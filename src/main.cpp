#include "BSP.h"

void read_OFF_file(const char* filename,
    double** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose) {

    FILE* file = fopen(filename, "r");
    if (file == NULL)
        ip_error("read_OFF_file: FATAL ERROR "
            "cannot open input file.\n");

    // Check OFF mark (1st line).
    char file_ext_read[3];
    char file_ext_target[] = { 'O','F','F' };
    if (fscanf(file, "%3c", file_ext_read) == 0)
        ip_error("read_OFF_file: FATAL ERROR "
            "cannot read 1st line of input file\n");

    for (uint32_t i = 0; i < 3; i++)
        if (file_ext_read[i] != file_ext_target[i])
            ip_error("read_OFF_file: FATAL ERROR "
                "1st line of input file is different from OFF\n");

    // Reading number of points and triangles.
    if (fscanf(file, " %d %d %*d ", npts, ntri) == 0)
        ip_error("read_OFF_file: FATAL ERROR 2st line of "
            "input file do not contanins point and triangles numbers.\n");

    if (verbose) printf("file %s contains %d vertices and %d constraints (triangles)\n",
        filename, *npts, *ntri);

    // Reading points coordinates.
    *vertices_p = (double*)malloc(sizeof(double) * 3 * (*npts));
    *tri_vertices_p = (uint32_t*)malloc(sizeof(uint32_t) * 3 * (*ntri));

    for (uint32_t i = 0; i < (*npts); i++) {
        if (fscanf(file, " %lf %lf %lf ",
            (*vertices_p) + (i * 3), (*vertices_p) + (i * 3 + 1), (*vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
    }

    uint32_t nv;
    for (uint32_t i = 0; i < (*ntri); i++) {
        if (fscanf(file, " %u %u %u %u ", &nv,
            (*tri_vertices_p) + (i * 3), (*tri_vertices_p) + (i * 3 + 1), (*tri_vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
        if (nv != 3) ip_error("Non-triangular faces not supported\n");
    }
    fclose(file);
}

void read_TET_file(const char* filename, double** vertices_p, uint32_t* npts,
    uint32_t** tet_vertices_p, uint32_t* ntet, bool verbose) {
    FILE* file = fopen(filename, "r");

    uint32_t nv, nt;
    if (!fscanf(file, "%d vertices\n", &nv)) ip_error("read_TET_file: FATAL ERROR\ncannot read 1st line of input file\n");
    if (!fscanf(file, "%d tets\n", &nt)) ip_error("read_TET_file: FATAL ERROR\ncannot read 2nd line of input file\n");

    if (verbose) printf("file %s contains %u vertices and %u tetrahedra\n",
        filename, nv, nt);

    double* v_crd = (double*)malloc(nv * sizeof(double) * 3);

    for (uint32_t i = 0; i < nv; i++) {
        double* crd = v_crd + i * 3;
        if (fscanf(file, "%lf %lf %lf\n", crd, crd + 1, crd + 2) != 3) ip_error("read_TET_file: FATAL ERROR\ncannot read vertex coordinates\n");
    }

    uint32_t* t_idx = (uint32_t*)malloc(nt * sizeof(uint32_t) * 4);

    for (uint32_t i = 0; i < nt; i++) {
        uint32_t* nodes = t_idx + i * 4;
        if (fscanf(file, "4 %u %u %u %u\n", nodes, nodes + 1, nodes + 2, nodes + 3) != 4) ip_error("read_TET_file: FATAL ERROR\ncannot read tet indexes\n");
    }

    fclose(file);

    *vertices_p = v_crd;
    *npts = nv;
    *tet_vertices_p = t_idx;
    *ntet = nt;
}

void read_MEDIT_file(const char* filename, double** vertices_p, uint32_t* npts,
    uint32_t** tet_vertices_p, uint32_t* ntet, bool verbose) {
    FILE* file = fopen(filename, "r");

    char line[1024];
    do {
        if (!fscanf(file, "%s\n", line)) ip_error("read_MEDIT_file: Could not find Vertices keyword\n");
    } while (strcmp(line, "Vertices"));

    int dummy, nv = 0;
    if (!fscanf(file, "%d\n", &nv)) ip_error("read_MEDIT_file: Could not read num vertices\n");
    printf("Reading %d vertices\n", nv);

    double* v_crd = (double*)malloc(nv * sizeof(double) * 3);

    for (int i = 0; i < nv; i++) {
        double* crd = v_crd + i * 3;
        if (fscanf(file, "%lf %lf %lf %d\n", crd, crd + 1, crd + 2, &dummy) != 4) ip_error("read_MEDIT_file: Could not read vertex coords\n");
    }

    int nt;
    if (!fscanf(file, "%s\n", line) || strcmp(line, "Tetrahedra")) ip_error("read_MEDIT_file: Could not find Tetrahedra keyword\n");
    if (!fscanf(file, "%d\n", &nt)) ip_error("read_MEDIT_file: Could not read num tet\n");
    printf("Reading %d tetrahedra\n", nt);

    uint32_t* t_idx = (uint32_t*)malloc(nt * sizeof(uint32_t) * 4);

    for (uint32_t i = 0; i < (uint32_t)nt; i++) {
        uint32_t* nodes = t_idx + i * 4;
        if (fscanf(file, "%d %d %d %d %d\n", nodes, nodes + 1, nodes + 2, nodes + 3, &dummy) != 5)
            ip_error("read_MEDIT_file: Could not read tet indexes\n");
    }
    for (uint32_t i = 0; i < (uint32_t)nt * 4; i++) t_idx[i]--; // Decrease all indexes by one

    fclose(file);

    if (verbose) printf("file %s contains %u vertices and %u tetrahedra\n",
        filename, nv, nt);

    *vertices_p = v_crd;
    *npts = nv;
    *tet_vertices_p = t_idx;
    *ntet = nt;
}

#include "embed.h"

/// <summary>
/// Main function
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(int argc, char** argv)
{
interval_number a(0);

#ifdef NDEBUG
    if (argc < 2) {
        printf("\nUsage: mesh_generator [-v | -l | -s | -b | -t] inputfile_A.off [bool_opcode inputfile_B.off]\n\n"
            "Defines the volume enclosed by the input OFF file(s) and saves a volume mesh to 'volume.msh'\n\n"
            "Command line arguments:\n"
            "-v = verbose mode\n"
            "-l = logging mode (appends a line to mesh_generator.log)\n"
            "-s = save the mesh bounding surface to 'skin.off'\n"
            "-b = save the subdivided constraints to 'black_faces.off'\n"
            "-t = triangulate/tetrahedrize output\n"
            "bool_opcode: {U, I, D}\n"
            "  U -> union (AuB),\n"
            "  I -> intersection (A^B),\n"
            "  D -> difference (A\\B)\n\n"
            "Example:\n"
            "mesh_generator ant.off U pig.off\n");
        return 0;
    }

    bool triangulate = false;
    bool verbose = false;
    bool logging = false;
    bool surfmesh = false;
    bool blackfaces = false;
    char* fileA_name = NULL;
    char* fileB_name = NULL;
    char bool_opcode = '0';

    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (argv[i][1] == 'v') verbose = true;
            else if (argv[i][1] == 't') triangulate = true;
            else if (argv[i][1] == 'l') logging = true;
            else if (argv[i][1] == 'b') blackfaces = true;
            else if (argv[i][1] == 's') surfmesh = true;
            else ip_error("Unknown option\n");
        }
        else if (fileA_name == NULL) fileA_name = argv[i];
        else if (bool_opcode == '0') bool_opcode = argv[i][0];
        else if (fileB_name == NULL) fileB_name = argv[i];
        else ip_error("Too many args passed\n");
    }

    bool two_input = (bool_opcode != '0');

    if (verbose)
    {
        if (fileB_name == NULL) {
            printf("\nResolve auto-intersections and/or repair.\n\n");
            printf("Loading %s\n", fileA_name);
        }
        else {
            printf("\nBoolean operator: ");
            if (bool_opcode == 'U') printf(" union.\n\n");
            else if (bool_opcode == 'I') printf(" intersection.\n\n");
            else if (bool_opcode == 'D') printf(" difference.\n\n");
            else { printf("INVALID\n\n"); return 0; }
            printf("Loading %s and %s.\n\n", fileA_name, fileB_name);
        }
    }
#else
    bool triangulate = true;
    bool verbose = true;
    bool logging = false;
    bool surfmesh = false;
    bool blackfaces = false;
    char* fileA_name = "D:\\SYNC_DATA\\Sviluppo_Software\\My_Software\\GIT_REPOS\\VolumeRemesher\\models\\Octocat-v1.off";
    char* fileB_name = "D:\\SYNC_DATA\\Sviluppo_Software\\My_Software\\GIT_REPOS\\VolumeRemesher\\models\\Octocat.bg.tet";
    char bool_opcode = 'U';
    bool two_input = (bool_opcode != '0');
#endif

    double* coords_A, * coords_B = NULL;
    uint32_t ncoords_A = 0, ncoords_B = 0;
    uint32_t* tri_idx_A = NULL, * tri_idx_B = NULL;
    uint32_t ntriidx_A = 0, ntriidx_B = 0;

    const bool file_B_is_tet = fileB_name != NULL && !strcmp(fileB_name + strlen(fileB_name) - 4, ".tet");
    const bool file_B_is_msh = fileB_name != NULL && !strcmp(fileB_name + strlen(fileB_name) - 4, ".msh");

    bool embedsurf = (file_B_is_tet || file_B_is_msh);
    if (embedsurf) two_input = false;

    read_OFF_file(fileA_name, &coords_A, &ncoords_A, &tri_idx_A, &ntriidx_A, verbose);
    if (two_input) read_OFF_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, verbose);

    BSPcomplex* complex;
    if (embedsurf) {
        if (file_B_is_tet) read_TET_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, verbose);
        else read_MEDIT_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, verbose);

        //const std::vector<double> tri_vrt_coords(coords_A, coords_A + ncoords_A * 3);
        //const std::vector<uint32_t> triangle_indexes(tri_idx_A, tri_idx_A + ntriidx_A * 3);
        //const std::vector<double> tet_vrt_coords(coords_B, coords_B + ncoords_B * 3);
        //const std::vector<uint32_t> tet_indexes(tri_idx_B, tri_idx_B + ntriidx_B * 4);
        //std::vector<bigrational> vertices;
        //std::vector<uint32_t> facets;
        //std::vector<uint32_t> cells;
        //std::vector<uint32_t> facets_on_input;

        //embed_tri_in_poly_mesh(
        //    tri_vrt_coords,
        //    triangle_indexes,
        //    tet_vrt_coords,
        //    tet_indexes,
        //    vertices,
        //    facets,
        //    cells,
        //    facets_on_input,
        //    verbose
        //);

        complex = remakePolyhedralMesh(
            coords_A, ncoords_A, tri_idx_A, ntriidx_A,
            coords_B, ncoords_B, tri_idx_B, ntriidx_B, verbose
        );
        free(coords_A);
        free(tri_idx_A);
    }
    else complex = makePolyhedralMesh(
        fileA_name, coords_A, ncoords_A, tri_idx_A, ntriidx_A,
        fileB_name, coords_B, ncoords_B, tri_idx_B, ntriidx_B,
        bool_opcode, true, verbose, logging
        );

    printf("Writing output file...\n");
    if (blackfaces) complex->saveBlackFaces("black_faces.off", triangulate);
    else if (surfmesh) complex->saveSkin("skin.off", bool_opcode, triangulate);
    else complex->saveMesh((triangulate)?("volume.tet"):("volume.msh"), bool_opcode, triangulate);
    printf("Done.\n");

    return 0;
}
