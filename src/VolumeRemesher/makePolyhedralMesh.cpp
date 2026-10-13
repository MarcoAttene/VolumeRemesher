#include <stdarg.h>
#include <time.h>
#include <cassert>
#include <cfloat>
#include <cstdio>
#include "BSP.h"
#include "string.h"

namespace vol_rem {
#ifdef _MSC_VER
// clang-format off
#include <windows.h>
#include <psapi.h>
// clang-format on

// To ensure correct resolution of symbols, add Psapi.lib to TARGETLIBS
// and compile with -DPSAPI_VERSION=1

double getPeakMegabytesUsed()
{
    HANDLE hProcess =
        OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
    if (NULL == hProcess) return 0;

    PROCESS_MEMORY_COUNTERS pmc;
    double mem = 0;
    if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
        mem = pmc.PeakWorkingSetSize / 1048576.0;
    }

    CloseHandle(hProcess);
    return mem;
}

void PrintMemoryInfo()
{
    printf("Peak memory usage (Mb): %f\n", getPeakMegabytesUsed());
}
#else
void PrintMemoryInfo()
{
    printf("Peak memory usage unavailable. This feature is implemented on Windows/MSVC only.\n");
}
#endif

// Vertex type

struct vertex_t
{
    double coord[3]; // Coordinates
    uint64_t inc_tet; // Incident tetrahedron
    uint32_t original_index; // Index to support reordering
};

int vertex_compare(const void* void_v1, const void* void_v2)
{
    const vertex_t* v1 = (vertex_t*)void_v1;
    const vertex_t* v2 = (vertex_t*)void_v2;
    const double dx = v1->coord[0] - v2->coord[0];
    const double dy = v1->coord[1] - v2->coord[1];
    const double dz = v1->coord[2] - v2->coord[2];
    return (4 * ((dx > 0) - (dx < 0)) + 2 * ((dy > 0) - (dy < 0)) + ((dz > 0) - (dz < 0)));
}

inline bool coincident_points(const vertex_t* a, const vertex_t* b)
{
    return !vertex_compare(a->coord, b->coord);
}

// Report subdivision progress every this many cells created. Ordinary models finish well
// inside one step and so stay silent; only a subdivision that is actually running long
// says anything.
static const uint64_t SUBDIVISION_PROGRESS_STEP = (uint64_t)1 << 20;

// The cell-subdivision loop is otherwise completely silent, so a slow input is
// indistinguishable from a stuck one -- both look like a hang with growing memory.
//
// "retired" is the loop cursor: cells already known to be free of constraints and hence
// final. Its share of the total is the useful number to watch, because subdivision ends
// when it reaches 100%. A run that is converging shows it climbing; one that is not shows
// it flat.
static void report_subdivision_progress(uint64_t created, uint64_t retired, clock_t started)
{
    printf(
        "\tsubdividing: %llu cells, %llu retired (%.1f%%), %.0f s\n",
        (unsigned long long)created,
        (unsigned long long)retired,
        created ? 100.0 * (double)retired / (double)created : 0.0,
        (double)(clock() - started) / CLOCKS_PER_SEC);
    fflush(stdout);
}

void remove_duplicated_points(
    vertex_t** vertices_p,
    uint32_t* npts,
    vertex_t* vrts_copy,
    uint32_t* map,
    uint32_t* diff)
{
    // Sorting vertices by coordinates lexicographic order (x,y,z)
    qsort(vrts_copy, *npts, sizeof(vertex_t), vertex_compare);

    // Count and memory position of duplicated vertices.
    uint32_t vrts_counter = 0;
    for (uint32_t i = 1; i < (*npts); i++) {
        if (coincident_points(vrts_copy + i - 1, vrts_copy + i)) vrts_counter++;
        diff[i] = vrts_counter;
    }

    // Set original_index to follow vertices permutation.
    for (uint32_t i = 0; i < (*npts); i++) map[vrts_copy[i].original_index] = i;

    // Allocating memory to store uinque mesh vertices (vertices_p).
    *vertices_p = (vertex_t*)malloc(sizeof(vertex_t) * (*npts - vrts_counter));

    // Fill mesh vertices (vertices_p)
    memcpy(*vertices_p, vrts_copy, sizeof(vertex_t));
    for (uint32_t i = vrts_counter = 1; i < (*npts); i++)
        if (!coincident_points(vrts_copy + i - 1, vrts_copy + i))
            memcpy(*vertices_p + (vrts_counter++), vrts_copy + i, sizeof(vertex_t));

    // Update uinque vertices number.
    (*npts) = vrts_counter;
}


/// //////////////////////////////////////////////////////////////////////////////////////////

bool misAlignment(const double* p, const double* q, const double* r)
{
    explicitPoint3D ep(p[0], p[1], p[2]), eq(q[0], q[1], q[2]), er(r[0], r[1], r[2]);
    return explicitPoint3D::misaligned(ep, eq, er);
}

void read_nodes_and_constraints(
    const double* coords_A,
    uint32_t npts_A,
    const uint32_t* tri_idx_A,
    uint32_t ntri_A,
    vertex_t** vertices_p,
    uint32_t* npts,
    uint32_t** tri_vertices_p,
    uint32_t* ntri,
    uint32_t** tri_original_p,
    bool verbose)
{
    // Reading points coordinates.
    *npts = npts_A;
    *ntri = ntri_A;
    vertex_t* tmp = (vertex_t*)malloc(*npts * sizeof(vertex_t));
    uint32_t* diff = (uint32_t*)calloc(*npts, sizeof(uint32_t));
    uint32_t* map = (uint32_t*)malloc(*npts * sizeof(uint32_t));
    *tri_vertices_p = (uint32_t*)malloc(sizeof(uint32_t) * 3 * (*ntri));
    *tri_original_p = (uint32_t*)malloc(sizeof(uint32_t) * (*ntri));

    for (uint32_t i = 0; i < (*npts); i++) {
        tmp[i].coord[0] = coords_A[i * 3];
        tmp[i].coord[1] = coords_A[i * 3 + 1];
        tmp[i].coord[2] = coords_A[i * 3 + 2];
        tmp[i].original_index = i;
    }

    if (*npts > 1) remove_duplicated_points(vertices_p, npts, tmp, map, diff);
    if (verbose) printf("Using %u unique vertices\n", *npts);

    free(tmp);

    for (uint32_t i = 0, j = 0; i < (*ntri); j++) {
        const uint32_t i1 = tri_idx_A[j * 3];
        const uint32_t i2 = tri_idx_A[j * 3 + 1];
        const uint32_t i3 = tri_idx_A[j * 3 + 2];
        (*tri_vertices_p)[3 * i] = map[i1] - diff[map[i1]];
        (*tri_vertices_p)[3 * i + 1] = map[i2] - diff[map[i2]];
        (*tri_vertices_p)[3 * i + 2] = map[i3] - diff[map[i3]];

        const double* v1c = ((*vertices_p) + (*tri_vertices_p)[3 * i])->coord;
        const double* v2c = ((*vertices_p) + (*tri_vertices_p)[3 * i + 1])->coord;
        const double* v3c = ((*vertices_p) + (*tri_vertices_p)[3 * i + 2])->coord;

        if (!misAlignment(v1c, v2c, v3c))
            (*ntri)--;
        else {
            (*tri_original_p)[i] = j; // remember this constraint's input-file index
            i++;
        }
    }
    free(map);
    free(diff);

    if (verbose) printf("Using %u non-degenerate constraints\n", *ntri);

    //// DEBUG
    // FILE *file = fopen("prova.off", "w");
    // fprintf(file, "OFF\n%u %u 0\n", *npts, *ntri);
    // for (uint32_t i = 0; i < (*npts); i++) fprintf(file, "%f %f %f\n", (*vertices_p)[i].coord[0],
    // (*vertices_p)[i].coord[1], (*vertices_p)[i].coord[2]); for (uint32_t i = 0; i < (*ntri); i++)
    // fprintf(file, "3 %u %u %u\n", (*tri_vertices_p)[3 * i], (*tri_vertices_p)[3 * i + 1],
    // (*tri_vertices_p)[3 * i + 2]); fclose(file); exit(0);
}

void read_nodes_and_constraints_twoInput(
    double* coords_A,
    uint32_t npts_A,
    uint32_t* tri_idx_A,
    uint32_t ntri_A,
    double* coords_B,
    uint32_t npts_B,
    uint32_t* tri_idx_B,
    uint32_t ntri_B,
    vertex_t** vertices_p,
    uint32_t* npts,
    uint32_t** tri_vertices_p,
    uint32_t* ntri,
    CONSTR_GROUP_T** tri_group,
    bool verbose)
{
    // Global numbers of points and triangles.
    *npts = npts_A + npts_B;
    *ntri = ntri_A + ntri_B;

    // Reading points coordinates.
    vertex_t* tmp = (vertex_t*)malloc(*npts * sizeof(vertex_t));
    uint32_t* diff = (uint32_t*)calloc(*npts, sizeof(uint32_t));
    uint32_t* map = (uint32_t*)malloc(*npts * sizeof(uint32_t));

    for (uint32_t i = 0; i < npts_A; i++) {
        tmp[i].coord[0] = coords_A[i * 3];
        tmp[i].coord[1] = coords_A[i * 3 + 1];
        tmp[i].coord[2] = coords_A[i * 3 + 2];
        tmp[i].original_index = i;
    }
    for (uint32_t i = 0; i < npts_B; i++) {
        tmp[i + npts_A].coord[0] = coords_B[i * 3];
        tmp[i + npts_A].coord[1] = coords_B[i * 3 + 1];
        tmp[i + npts_A].coord[2] = coords_B[i * 3 + 2];
        tmp[i + npts_A].original_index = i + npts_A;
    }

    if (*npts > 1) remove_duplicated_points(vertices_p, npts, tmp, map, diff);
    if (verbose) printf("Using %u unique vertices\n", *npts);

    free(tmp);

    // Reading triangle vertices indices.
    *tri_vertices_p = (uint32_t*)malloc(*ntri * 3 * sizeof(uint32_t));
    *tri_group = (uint32_t*)malloc(*ntri * sizeof(uint32_t));

    uint32_t i1, i2, i3, base_ntri_A = ntri_A;
    for (uint32_t ti = 0, j = 0; ti < (*ntri); j++) {
        if (ti < ntri_A) {
            i1 = tri_idx_A[3 * j];
            i2 = tri_idx_A[3 * j + 1];
            i3 = tri_idx_A[3 * j + 2];

            (*tri_vertices_p)[3 * ti] = map[i1] - diff[map[i1]];
            (*tri_vertices_p)[3 * ti + 1] = map[i2] - diff[map[i2]];
            (*tri_vertices_p)[3 * ti + 2] = map[i3] - diff[map[i3]];
            (*tri_group)[ti] = CONSTR_A;

            const double* v1c = (*vertices_p + (*tri_vertices_p)[3 * ti])->coord;
            const double* v2c = (*vertices_p + (*tri_vertices_p)[3 * ti + 1])->coord;
            const double* v3c = (*vertices_p + (*tri_vertices_p)[3 * ti + 2])->coord;
            if (!misAlignment(v1c, v2c, v3c)) {
                (*ntri)--;
                ntri_A--;
            } else
                ti++;

        } else {
            if (ti == ntri_A) j = 0;
            i1 = tri_idx_B[3 * j];
            i2 = tri_idx_B[3 * j + 1];
            i3 = tri_idx_B[3 * j + 2];

            i1 += npts_A;
            i2 += npts_A;
            i3 += npts_A;
            (*tri_vertices_p)[3 * ti] = map[i1] - diff[map[i1]];
            (*tri_vertices_p)[3 * ti + 1] = map[i2] - diff[map[i2]];
            (*tri_vertices_p)[3 * ti + 2] = map[i3] - diff[map[i3]];
            (*tri_group)[ti] = CONSTR_B;

            const double* v1c = (*vertices_p + (*tri_vertices_p)[3 * ti])->coord;
            const double* v2c = (*vertices_p + (*tri_vertices_p)[3 * ti + 1])->coord;
            const double* v3c = (*vertices_p + (*tri_vertices_p)[3 * ti + 2])->coord;
            if (!misAlignment(v1c, v2c, v3c))
                (*ntri)--;
            else
                ti++;
        }
    }

    *tri_vertices_p = (uint32_t*)realloc(*tri_vertices_p, *ntri * 3 * sizeof(uint32_t));
    *tri_group = (uint32_t*)realloc(*tri_group, *ntri * sizeof(uint32_t));
    free(map);
    free(diff);

    if (verbose) printf("Using %u non-degenerate constraints\n", *ntri);
}


void logger(FILE* fp, const char* msg, ...)
{
    char fmt[2048], fms[4096];
    va_list ap;
    va_start(ap, msg);
    strcpy(fmt, msg);
    vsprintf(fms, fmt, ap);

    fprintf(fp, fms);

    va_end(ap);
}


// void saveVRML(const char* filename, double* coords_A, uint32_t npts_A, uint32_t* tri_idx_A,
// uint32_t ntri_A, bool wireframe)
//{
//     FILE* fp = fopen(filename, "w");
//     fprintf(fp,
//         "#VRML V1.0 ascii\n\n");
//     fprintf(fp,
//         "Separator {\n"
//         "ShapeHints{\n"
//         "vertexOrdering COUNTERCLOCKWISE creaseAngle 0.0\n"
//         "shapeType SOLID\n"
//         "}\n"
//         "Separator{\n"
//         "Coordinate3 {\n"
//         "point[\n");
//
//     for (uint32_t i = 0; i < npts_A; i++)
//         fprintf(fp, "%f %f %f,\n", coords_A[3 * i], coords_A[3 * i + 1], coords_A[3 * i + 2]);
//
//     fprintf(fp,
//         "]\n"
//         "}\n"
//         "Material{\n"
//         "diffuseColor 0.7 0.5 0.5\n"
//         "}\n"
//         "IndexedFaceSet{\n"
//         "coordIndex[\n");
//
//     for (uint32_t i = 0; i < ntri_A; i++)
//         fprintf(fp, "%u, %u, %u, -1,\n", tri_idx_A[3 * i], tri_idx_A[3 * i + 1], tri_idx_A[3 * i
//         + 2]);
//
//     fprintf(fp,
//         "]\n"
//         "}\n");
//
//     if (wireframe)
//     {
//
//         fprintf(fp,
//             "Material{\n"
//             "diffuseColor 0.0 0.0 0.0\n"
//             "}\n"
//             "DrawStyle{\n"
//             "lineWidth 2\n"
//             "}\n"
//             "IndexedLineSet{\n"
//             "coordIndex[\n");
//
//         for (uint32_t i = 0; i < ntri_A; i++)
//         {
//             uint32_t i1, i2;
//             i1 = tri_idx_A[3 * i];
//             i2 = tri_idx_A[3 * i + 1];
//             if (i1 < i2) fprintf(fp, "%u, %u, -1,\n", i1, i2);
//             i1 = tri_idx_A[3 * i + 1];
//             i2 = tri_idx_A[3 * i + 2];
//             if (i1 < i2) fprintf(fp, "%u, %u, -1,\n", i1, i2);
//             i1 = tri_idx_A[3 * i + 2];
//             i2 = tri_idx_A[3 * i];
//             if (i1 < i2) fprintf(fp, "%u, %u, -1,\n", i1, i2);
//         }
//
//         fprintf(fp, "]\n}\n");
//     }
//     fprintf(fp, "}\n}\n");
//
//     fclose(fp);
// }


/// <summary>
/// Main function - Create a polyhedral mesh out of the input
/// Input may be made of either one or two models to be combined into a boolean composition
/// </summary>
/// <param name="fileA_name">Name of first model</param>
/// <param name="coords_A">Serialized coordinates of first model vertices</param>
/// <param name="npts_A">Number of first model vertices</param>
/// <param name="tri_idx_A">Serialized indexes of first model triangles</param>
/// <param name="ntri_A">Number of first model triangles</param>
/// <param name="fileB_name">Name of second model</param>
/// <param name="coords_B">Serialized coordinates of second model vertices</param>
/// <param name="npts_B">Number of second model vertices</param>
/// <param name="tri_idx_B">Serialized indexes of second model triangles</param>
/// <param name="ntri_B">Number of second model triangles</param>
/// <param name="bool_opcode">Boolean operation (0 = no op, U = union, D = difference, I = intersection</param>
/// <param name="free_mem">Input models memory is released</param>
/// <param name="verbose">Print useful info during the process</param>
/// <param name="logging">Append info to mesh_generator.log</param>
/// <returns>The resulting BSPcomplex structure</returns>
BSPcomplex* makePolyhedralMesh(
    const char* fileA_name,
    double* coords_A,
    uint32_t npts_A,
    uint32_t* tri_idx_A,
    uint32_t ntri_A,
    const char* fileB_name,
    double* coords_B,
    uint32_t npts_B,
    uint32_t* tri_idx_B,
    uint32_t ntri_B,
    char bool_opcode,
    bool free_mem,
    bool verbose,
    bool logging,
    bool fill_holes,
    bool keep_all_cells,
    const extra_features_t* extra)
{
    // Inserted edges/points only survive if we keep every cell (the min-cut would
    // otherwise delete the cells they bound), so they imply keep_all_cells.
    const bool has_extra = extra && (extra->n_edges > 0 || extra->n_points > 0);
    if (has_extra) keep_all_cells = true;
    // saveVRML("input1.wrl", coords_A, npts_A, tri_idx_A, ntri_A, false);
    // saveVRML("input2.wrl", coords_A, npts_A, tri_idx_A, ntri_A, true);

    bool two_input = (bool_opcode != '0');

    if (verbose) {
        if (coords_B == NULL) {
            printf("\nResolve auto-intersections and/or repair.\n\n");
        } else {
            printf("\nBoolean operator: ");
            if (bool_opcode == 'U')
                printf(" union.\n\n");
            else if (bool_opcode == 'I')
                printf(" intersection.\n\n");
            else if (bool_opcode == 'D')
                printf(" difference.\n\n");
            else {
                ip_error("INVALID\n\n");
            }
        }
    }

    //--Initialization-----------------

    vertex_t* vertices;
    uint32_t num_vertices;
    constraints_t* constraints = new constraints_t;

    if (!two_input) {
        if (npts_A > 0) {
            read_nodes_and_constraints(
                coords_A,
                npts_A,
                tri_idx_A,
                ntri_A,
                &vertices,
                &num_vertices,
                &constraints->tri_vertices,
                &constraints->num_triangles,
                &constraints->tri_original_index,
                verbose);
        } else {
            // No triangle input: start with an empty surface. The edges/points inserted
            // below supply all vertices and constraints. (This is what makes the triangle
            // input optional -- only edges, only points, or any combination is valid.)
            vertices = NULL;
            num_vertices = 0;
            constraints->tri_vertices = NULL;
            constraints->num_triangles = 0;
            constraints->tri_original_index = NULL;
        }
    } else { // two input
        read_nodes_and_constraints_twoInput(
            coords_A,
            npts_A,
            tri_idx_A,
            ntri_A,
            coords_B,
            npts_B,
            tri_idx_B,
            ntri_B,
            &vertices,
            &num_vertices,
            &constraints->tri_vertices,
            &constraints->num_triangles,
            &constraints->constr_group,
            verbose);
    }

    TetMesh* mesh = new TetMesh;
    std::vector<double> coordinates(num_vertices * 3);
    for (uint32_t i = 0, j = 0; i < num_vertices; i++) {
        coordinates[j++] = vertices[i].coord[0];
        coordinates[j++] = vertices[i].coord[1];
        coordinates[j++] = vertices[i].coord[2];
    }
    mesh->init_vertices(coordinates.data(), num_vertices);

    free(vertices);

    if (!two_input) {
        // Cap open boundaries so that every input triangle bounds a closed volume
        // (see fill_holes_in_constraints). The cap triangles are excluded from the
        // surface tracking. Must run before the Delaunay permutation and virtual
        // constraints below, so the caps are remapped/handled just like real input.
        if (fill_holes) fill_holes_in_constraints(constraints, mesh, verbose);
        // Force extra edges/points into the output (appended after the caps). Same
        // timing requirement: before the Delaunay permutation and virtual constraints.
        if (has_extra) {
            insert_edges_and_points(constraints, mesh, *extra, verbose);
            // Edge/point forcing features can be the only input (an edges-only or points-only
            // run has no surface triangles or vertices of its own). init_vertices computed the
            // static predicate filters over just the surface vertices -- of which there may
            // have been none, leaving max_coord at its DBL_MAX sentinel -- so recompute them
            // now that the forcing vertices are present.
            mesh->init_static_filters();
        }
        constraints->constr_group = (uint32_t*)calloc(constraints->num_triangles, sizeof(uint32_t));
    }

    // Validate now that any edge/point forcing features have been added: an edges-only or
    // points-only input has no surface triangles or vertices, so these counts are only
    // meaningful here. (Before the Delaunay3D externalization this check sat after
    // insert_edges_and_points; the refactor moved it earlier, which rejected no-surface input.)
    if (mesh->numVertices() < 4) ip_error("Cannot mesh less than 4 vertices.");
    if (constraints->num_triangles < 1) ip_error("No non-degenerate constraints loaded.");

    if (free_mem) {
        free(coords_A);
        free(tri_idx_A);
        if (two_input) {
            free(coords_B);
            free(tri_idx_B);
        }
    }
    FILE* logfile;

    if (logging) {
        logfile = fopen("mesh_generator.log", "r");
        if (logfile == NULL) {
            logfile = fopen("mesh_generator.log", "a");
            logger(
                logfile,
                "Filename,Delaunay,Virtual,Mapping,Total Delaunay,Conversion,Subdivision,Face "
                "Coloring,In-Out,Total BSP,TOTAL,Mem\n");
        } else
            logfile = freopen("mesh_generator.log", "a", logfile);

        if (two_input)
            logger(logfile, "%s_%c_%s,", fileA_name, bool_opcode, fileB_name);
        else
            logger(logfile, "%s,", fileA_name);
        fflush(logfile);
    }

    clock_t time0 = clock();
    clock_t time1 = time0;

    //--Delaunay-Insertion-----------------

    // Setup for following vertex permutation in tetrahedrize
    //for (uint32_t i = 0; i < mesh->num_vertices; i++) mesh->vertices[i].original_index = i;

    // Create Delaunay tetrahedrization of the vertices
    mesh->tetrahedrize();
    // mesh->saveTET("delaunay.tet");

    //// Align constraint vertices after vertex permutation

    //if (mesh->vertices[2].original_index != 3) {
    //    for (uint32_t k = 0; k < 3 * constraints->num_triangles; k++)
    //        constraints->tri_vertices[k] =
    //            mesh->vertices[constraints->tri_vertices[k]].original_index;
    //} else {
    //    for (uint32_t k = 0; k < 3 * constraints->num_triangles; k++) {
    //        uint32_t l = mesh->vertices[3].original_index;
    //        uint32_t c = constraints->tri_vertices[k];
    //        if (c == 2)
    //            constraints->tri_vertices[k] = l;
    //        else if (c == 3)
    //            constraints->tri_vertices[k] = 2;
    //        else if (c == l)
    //            constraints->tri_vertices[k] = 3;
    //    }
    //}


    clock_t time2 = clock();
    if (verbose) printf("\tDelaunay insertion: %f s\n", (double)(time2 - time1) / CLOCKS_PER_SEC);

    //--Half-Edges-and-Virtual-Constraints----------------------

    half_edge_t* half_edges =
        (half_edge_t*)calloc(3 * constraints->num_triangles, sizeof(half_edge_t));
    fill_half_edges(constraints, half_edges);

    sort_half_edges(half_edges, 3 * constraints->num_triangles);

    uint32_t nvc = place_virtual_constraints(mesh, constraints, half_edges);
    if (verbose) printf("\t%u virtual constraints added\n", nvc);

    free(half_edges);

    clock_t time3 = clock();
    if (verbose) printf("\tHalf-edges: %f s\n", (double)(time3 - time2) / CLOCKS_PER_SEC);

    //--Map-Tetrahedra-Constraint-Intersections----------------

    // Initialize an array to map the tetrahedra improperly intersecated
    // by the constraints:
    //    the i-th element is a pointer (of tetrahedra index type) used to store
    //    the indices of the constraints that improper intersect the i-th
    //    tetrahedron,
    //    it points to an array (map[i] of lenght num_map[i]) whose elements
    //    are the indices of the constraints which improperly intersect the
    //    i-th tetrahedron.
    uint32_t* num_map = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));

    // Initalize other 4 maps to memory the tet-faces that are partially or
    // completely overlap with a constraint.
    uint32_t* num_map_f0 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f0 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));
    uint32_t* num_map_f1 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f1 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));
    uint32_t* num_map_f2 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f2 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));
    uint32_t* num_map_f3 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f3 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));

    insert_constraints(
        mesh,
        constraints,
        num_map,
        map,
        num_map_f0,
        map_f0,
        num_map_f1,
        map_f1,
        num_map_f2,
        map_f2,
        num_map_f3,
        map_f3);

    clock_t time4 = clock();
    if (verbose) printf("\tMap creation: %f s\n", (double)(time4 - time3) / CLOCKS_PER_SEC);

    double DEL_time = (double)(time4 - time0) / CLOCKS_PER_SEC;
    if (verbose) printf("TOTAL Delaunay+map: %f s\n\n", DEL_time);

    initFPU(); // From here on we need indirect predicates

    //-Init BSP with mesh and constraints---------------------------------------
    BSPcomplex* complex_p = new BSPcomplex(
        mesh,
        constraints,
        (const uint32_t**)map,
        num_map,
        (const uint32_t**)map_f0,
        num_map_f0,
        (const uint32_t**)map_f1,
        num_map_f1,
        (const uint32_t**)map_f2,
        num_map_f2,
        (const uint32_t**)map_f3,
        num_map_f3);
    BSPcomplex& complex = *complex_p;

    // Free the memory used by map(s) and num_map(s)
    for (uint64_t i = 0; i < mesh->numTets(); i++) {
        if (num_map[i]) free(map[i]);
        if (num_map_f0[i]) free(map_f0[i]);
        if (num_map_f1[i]) free(map_f1[i]);
        if (num_map_f2[i]) free(map_f2[i]);
        if (num_map_f3[i]) free(map_f3[i]);
    }
    free(num_map);
    free(map);
    free(num_map_f0);
    free(map_f0);
    free(num_map_f1);
    free(map_f1);
    free(num_map_f2);
    free(map_f2);
    free(num_map_f3);
    free(map_f3);

    delete mesh;
    delete constraints;

    clock_t time5 = clock();
    if (verbose) printf("\tDelaunay -> Complex %lf s\n", (double)(time5 - time4) / CLOCKS_PER_SEC);
    if (verbose) printf("\tInitial cells: %llu\n", complex.cells.size());

    //-Subdivision----------------------------------------------------------------

    // #define MEM_THRESHOLD   32768

    uint64_t cell_ind = 0;
    uint64_t next_progress = SUBDIVISION_PROGRESS_STEP;
    while (cell_ind < complex.cells.size()) {
#ifdef MEM_THRESHOLD
        if ((cell_ind & 0x0000000000000fff) == 0) {
            printf("\r%llu cells created                  ", cell_ind);
            fflush(stdout);
            if (getPeakMegabytesUsed() > MEM_THRESHOLD) ip_error("Maximum memory reached\n");
        }
#endif
        if (complex.cells[cell_ind].constraints.size() > 0)
            complex.splitCell(cell_ind);
        else
            cell_ind++;
        if (verbose && complex.cells.size() >= next_progress) {
            report_subdivision_progress(complex.cells.size(), cell_ind, time5);
            next_progress = complex.cells.size() + SUBDIVISION_PROGRESS_STEP;
        }
    }
    clock_t time6 = clock();
    if (verbose) printf("\tCell subdivision %f s\n", (double)(time6 - time5) / CLOCKS_PER_SEC);
    if (verbose) printf("\tFinal cells: %llu\n", complex.cells.size());


    // for (BSPcell& c : complex.cells) c.place = INTERNAL_A;
    // complex.saveMesh("splitch.tet", '0', true);
    // exit(0);


    //--Decide colour of GREY faces-----------------------------------------------
    for (uint64_t face_ind = 0; face_ind < complex.faces.size(); face_ind++) {
        BSPface& face = complex.faces[face_ind];
        if (face.colour == GREY)
            face.colour = complex.blackAB_or_white(face_ind, bool_opcode != '0');
    }

    clock_t time7 = clock();
    if (verbose) printf("\tFind black faces %f s\n", (double)(time7 - time6) / CLOCKS_PER_SEC);

    //-Classification:intrenal/external cells-------------------------------------
    if (keep_all_cells) {
        // Skip the min-cut in/out classification and keep the WHOLE tetrahedralized
        // domain. Every input triangle is then realized as a face shared by two kept
        // cells, so the BLACK-face surface tagging is exact: no cell is deleted, hence
        // no thin feature (e.g. an open-boundary flap) is shaved off by the
        // area-minimizing cut. Marking every cell INTERNAL_A makes makeTetrahedra keep
        // them all without any further plumbing.
        for (BSPcell& c : complex.cells) c.place = INTERNAL_A;
    } else
        complex.constraintsSurface_complexPartition(bool_opcode != '0');

    clock_t time8 = clock();
    if (verbose) printf("\tInt-ext class. %f s\n", (double)(time8 - time7) / CLOCKS_PER_SEC);


    clock_t time9 = clock();
    double BSP_time = (double)(time9 - time4) / CLOCKS_PER_SEC;
    if (verbose) printf("TOTAL BSP: %f s\n\n", BSP_time);
    if (verbose) printf("TOTAL time: %f s\n\n", BSP_time + DEL_time);
    if (verbose) PrintMemoryInfo();

    if (logging) {
        logger(
            logfile,
            "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
            (double)(time2 - time1) / CLOCKS_PER_SEC,
            (double)(time3 - time2) / CLOCKS_PER_SEC,
            (double)(time4 - time3) / CLOCKS_PER_SEC,
            DEL_time,
            (double)(time5 - time4) / CLOCKS_PER_SEC,
            (double)(time6 - time5) / CLOCKS_PER_SEC,
            (double)(time7 - time6) / CLOCKS_PER_SEC,
            (double)(time8 - time7) / CLOCKS_PER_SEC,
            BSP_time,
            BSP_time + DEL_time
#ifdef _MSC_VER
            ,
            getPeakMegabytesUsed()
#endif
        );
        fclose(logfile);
    }

    // uint32_t npts, ntri, * tri_idx;
    // double* coord;
    // complex.extractSkinTriMesh("input4.wrl", bool_opcode, &coord, &npts, &tri_idx, &ntri);
    // saveVRML("input4.wrl", coord, npts, tri_idx, ntri, false);
    // saveVRML("input5.wrl", coord, npts, tri_idx, ntri, true);
    // exit(0);

    return complex_p;
}


BSPcomplex* remakePolyhedralMesh(
    const double* coords_A,
    uint32_t npts_A,
    const uint32_t* tri_idx_A,
    uint32_t ntri_A,
    const double* coords_B,
    uint32_t npts_B,
    const uint32_t* tet_idx_B,
    uint32_t ntet_B,
    bool verbose,
    bool no_class,
    const extra_features_t* extra)
{
    if (verbose) printf("\nEmbed trimesh into tetmesh.\n\n");

    //--Initialization-----------------

    constraints_t* constraints = new constraints_t;
    vertex_t* vertices;
    uint32_t num_vertices;

    vertices = NULL;
    num_vertices = 0;
    if (npts_A > 0) {
        read_nodes_and_constraints(
            coords_A,
            npts_A,
            tri_idx_A,
            ntri_A,
            &vertices,
            &num_vertices,
            &constraints->tri_vertices,
            &constraints->num_triangles,
            &constraints->tri_original_index,
            verbose);
    } else {
        // No surface: the extra edges/points supply all constraints (optional surface).
        constraints->tri_vertices = NULL;
        constraints->num_triangles = 0;
        constraints->tri_original_index = NULL;
    }

    TetMesh* mesh = new TetMesh;

    std::vector<double> coordinates(npts_B * 3);
    for (size_t i = 0; i < coordinates.size(); i++) coordinates[i] = coords_B[i];
    mesh->init_vertices(coordinates.data(), npts_B);
    mesh->init_tets(tet_idx_B, ntet_B);

    uint32_t orig_numv = mesh->numVertices();
    if (verbose)
        printf("\nLoaded tetmesh made of %u vertices and %u tets.\n\n", orig_numv, ntet_B);

    if (npts_A > 0) {
        for (uint32_t i = 0; i < constraints->num_triangles * 3; i++)
            constraints->tri_vertices[i] += orig_numv;
  
        coordinates.resize(num_vertices * 3);
        for (uint32_t i = 0, j = 0; i < num_vertices; i++) {
            coordinates[j++] = vertices[i].coord[0];
            coordinates[j++] = vertices[i].coord[1];
            coordinates[j++] = vertices[i].coord[2];
        }

        mesh->vertices.reserve(mesh->numVertices() + num_vertices);
        for (uint32_t i = 0; i < num_vertices; i++)
            mesh->pushVertex(
                basicVec3d(coordinates[i * 3], coordinates[i * 3 + 1], coordinates[i * 3 + 2]));

        free(vertices);
    }


    // Force extra edges/points into the embedding: their vertices are appended here so
    // the insertExistingVerticesNonDelaunay call below inserts them along with the
    // surface, and their forcing triangles ride through the rest of the pipeline as
    // ordinary (non-surface) constraints. They must lie inside the tet mesh domain.
    const bool has_extra_embed = extra && (extra->n_edges > 0 || extra->n_points > 0);
    if (has_extra_embed)
        insert_edges_and_points(constraints, mesh, *extra, verbose, /*keep_corners_inside=*/true);

#ifndef NDEBUG
    // Debug guard: inserted edge/point vertices -- and their forcing-triangle apexes, which
    // are offset by the edge length (edges) or the average edge length (points) -- must lie
    // inside the tet mesh M. Otherwise the non-Delaunay insertion below cannot locate them
    // and crashes. This bounding-box test is a necessary condition (exact for a convex M)
    // and turns the common "apex outside M" mistake into a clear assertion.
    if (has_extra_embed) {
        double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
        for (uint32_t v = 0; v < orig_numv; v++)
            for (int k = 0; k < 3; k++) {
                const double c = mesh->vertices[v].coord[k];
                if (c < lo[k]) lo[k] = c;
                if (c > hi[k]) hi[k] = c;
            }
        for (uint32_t v = orig_numv; v < mesh->numVertices(); v++)
            for (int k = 0; k < 3; k++) {
                const double c = mesh->vertices[v].coord[k];
                if (c < lo[k] || c > hi[k]) {
                    fprintf(
                        stderr,
                        "remakePolyhedralMesh: an inserted edge/point vertex or forcing-triangle "
                        "apex lies outside the tet mesh's bounding box. Edges/points and their "
                        "apexes (offset by the edge length / average edge length) must be inside "
                        "the embedding tet mesh.\n");
                    assert(false && "inserted edge/point or apex is outside the tet mesh domain");
                }
            }
    }
#endif

    // Allocate constr_group after all constraints exist (surface + edge/point forcing).
    constraints->constr_group = (uint32_t*)calloc(constraints->num_triangles, sizeof(uint32_t));
    if (constraints->num_triangles < 1) ip_error("No constraints loaded (surface/edges/points).");

    // Insert constraint vertices into mesh
    //for (uint32_t i = 0; i < mesh->numVertices(); i++) mesh->vertices[i].original_index = i;
    uint64_t tet = 0;
    for (uint32_t vi = orig_numv; vi < mesh->numVertices(); vi++) {
        tet = mesh->searchTetrahedron(tet, vi);
        if (tet != UINT64_MAX) mesh->insertExistingVertexNonDelaunay(vi, tet);
    }
        

    //for (uint32_t i = 0; i < constraints->num_triangles * 3; i++) {
    //    constraints->tri_vertices[i] = mesh->vertices[constraints->tri_vertices[i]].original_index;
    //}

    // Here we may still have unused vertices (copies of coincident vertices that do not belong to
    // tets) Should clean... here below a sketch of the algorithm

    // uint32_t offset = 0;
    // std::vector<uint32_t> map;
    // map.resize(mesh->num_vertices);
    // for (uint32_t i = 0; i < mesh->num_vertices; i++) map[i] = i;

    // for (uint32_t i = 0; i < mesh->num_vertices; i++) {
    //     vertex_t& v = mesh->vertices[i];
    //     map[i] -= offset;
    //     if (v.original_index != i) {
    //         v.original_index = UINT32_MAX;
    //         offset++;
    //     }
    // }


    // 0) Init 'offset' to zero
    // 1) Create another vector<uint32_t> 'map' as large as mesh->vertices and init it with
    // increasing indexes 2) for each vertex v = vertices[i], decrease map[i] by offset
    //    Check whether 'v's original index corresponds to its position i.
    //    if not (the vertex must be eventually deleted):
    //      set original_index to UINT32_MAX
    //      increase 'offset'
    // 3) std::erase all vertices having original_index == UINT32_MAX
    // 4) replace each connstraints->tri_vertices[i] with map[constraints->tri_vertices[i]]

    clock_t time0 = clock();
    clock_t time2 = time0;

    //--Half-Edges-and-Virtual-Constraints----------------------

    half_edge_t* half_edges =
        (half_edge_t*)calloc(3 * constraints->num_triangles, sizeof(half_edge_t));
    fill_half_edges(constraints, half_edges);

    sort_half_edges(half_edges, 3 * constraints->num_triangles);

    uint32_t nvc = place_virtual_constraints(mesh, constraints, half_edges);
    if (verbose) printf("\t%u virtual constraints added\n", nvc);

    free(half_edges);

    clock_t time3 = clock();
    if (verbose) printf("\tHalf-edges: %f s\n", (double)(time3 - time2) / CLOCKS_PER_SEC);

    //--Map-Tetrahedra-Constraint-Intersections----------------

    // Initialize an array to map the tetrahedra improperly intersecated
    // by the constraints:
    //    the i-th element is a pointer (of tetrahedra index type) used to store
    //    the indices of the constraints that improper intersect the i-th
    //    tetrahedron,
    //    it points to an array (map[i] of lenght num_map[i]) whose elements
    //    are the indices of the constraints which improperly intersect the
    //    i-th tetrahedron.
    uint32_t* num_map = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));

    // Initalize other 4 maps to memory the tet-faces that are partially or
    // completely overlap with a constraint.
    uint32_t* num_map_f0 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f0 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));
    uint32_t* num_map_f1 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f1 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));
    uint32_t* num_map_f2 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f2 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));
    uint32_t* num_map_f3 = (uint32_t*)calloc(mesh->numTets(), sizeof(uint32_t));
    uint32_t** map_f3 = (uint32_t**)calloc(mesh->numTets(), sizeof(uint32_t*));

    insert_constraints(
        mesh,
        constraints,
        num_map,
        map,
        num_map_f0,
        map_f0,
        num_map_f1,
        map_f1,
        num_map_f2,
        map_f2,
        num_map_f3,
        map_f3);

    clock_t time4 = clock();
    if (verbose) printf("\tMap creation: %f s\n", (double)(time4 - time3) / CLOCKS_PER_SEC);

    double DEL_time = (double)(time4 - time0) / CLOCKS_PER_SEC;
    if (verbose) printf("TOTAL Delaunay+map: %f s\n\n", DEL_time);

    initFPU(); // From here on we need indirect predicates

    //-Init BSP with mesh and constraints---------------------------------------
    BSPcomplex* complex_p = new BSPcomplex(
        mesh,
        constraints,
        (const uint32_t**)map,
        num_map,
        (const uint32_t**)map_f0,
        num_map_f0,
        (const uint32_t**)map_f1,
        num_map_f1,
        (const uint32_t**)map_f2,
        num_map_f2,
        (const uint32_t**)map_f3,
        num_map_f3);
    BSPcomplex& complex = *complex_p;

    // Free the memory used by map(s) and num_map(s)
    for (uint64_t i = 0; i < mesh->numTets(); i++) {
        if (num_map[i]) free(map[i]);
        if (num_map_f0[i]) free(map_f0[i]);
        if (num_map_f1[i]) free(map_f1[i]);
        if (num_map_f2[i]) free(map_f2[i]);
        if (num_map_f3[i]) free(map_f3[i]);
    }
    free(num_map);
    free(map);
    free(num_map_f0);
    free(map_f0);
    free(num_map_f1);
    free(map_f1);
    free(num_map_f2);
    free(map_f2);
    free(num_map_f3);
    free(map_f3);

    delete mesh;
    delete constraints;

    clock_t time5 = clock();
    if (verbose) printf("\tDelaunay -> Complex %lf s\n", (double)(time5 - time4) / CLOCKS_PER_SEC);
    if (verbose) printf("\tInitial cells: %llu\n", complex.cells.size());

    //-Subdivision----------------------------------------------------------------

    // #define MEM_THRESHOLD   32768

    uint64_t cell_ind = 0;
    uint64_t next_progress = SUBDIVISION_PROGRESS_STEP;
    while (cell_ind < complex.cells.size()) {
#ifdef MEM_THRESHOLD
        if ((cell_ind & 0x0000000000000fff) == 0) {
            printf("\r%llu cells created                  ", cell_ind);
            fflush(stdout);
            if (getPeakMegabytesUsed() > MEM_THRESHOLD) ip_error("Maximum memory reached\n");
        }
#endif
        if (complex.cells[cell_ind].constraints.size() > 0)
            complex.splitCell(cell_ind);
        else
            cell_ind++;
        if (verbose && complex.cells.size() >= next_progress) {
            report_subdivision_progress(complex.cells.size(), cell_ind, time5);
            next_progress = complex.cells.size() + SUBDIVISION_PROGRESS_STEP;
        }
    }
    clock_t time6 = clock();
    if (verbose) printf("\tCell subdivision %f s\n", (double)(time6 - time5) / CLOCKS_PER_SEC);
    if (verbose) printf("\tFinal cells: %llu\n", complex.cells.size());


    // for (BSPcell& c : complex.cells) c.place = INTERNAL_A;
    // complex.saveMesh("splitch.tet", '0', true);
    // exit(0);


    //--Decide colour of GREY faces-----------------------------------------------
    for (uint64_t face_ind = 0; face_ind < complex.faces.size(); face_ind++) {
        BSPface& face = complex.faces[face_ind];
        if (face.colour == GREY) face.colour = complex.blackAB_or_white(face_ind, false);
    }

    clock_t time7 = clock();
    if (verbose) printf("\tFind black faces %f s\n", (double)(time7 - time6) / CLOCKS_PER_SEC);

    //-Classification:intrenal/external cells-------------------------------------
    if (!no_class) {
        complex.constraintsSurface_complexPartition(false);

        clock_t time8 = clock();
        if (verbose) printf("\tInt-ext class. %f s\n", (double)(time8 - time7) / CLOCKS_PER_SEC);
    }

    clock_t time9 = clock();
    double BSP_time = (double)(time9 - time4) / CLOCKS_PER_SEC;
    if (verbose) printf("TOTAL BSP: %f s\n\n", BSP_time);
    if (verbose) printf("TOTAL time: %f s\n\n", BSP_time + DEL_time);
    if (verbose) PrintMemoryInfo();

    return complex_p;
}
} // namespace vol_rem