#ifndef _CONFORMING_MESH_
#define _CONFORMING_MESH_

#include "delaunay3d_wrapper.h"

#define CONSTR_GROUP_T uint32_t
#define CONSTR_A 0
#define CONSTR_B 1

namespace vol_rem {
// half-edges struct
struct half_edge_t
{
    uint32_t endpts[2]; // endpoints of a constraint(triangle) edge,
    // be sure that endpts[0] < endpts[1].
    uint32_t tri_ind; // id_number of the constraint: tri_vertices[i]/3, where
    // tri_vertices[i]=endpts[0] or tri_vertices[i]=endpts[1].
};

// constraints are triangles
class constraints_t
{
public:
    uint32_t* tri_vertices;
    uint32_t num_triangles;
    uint32_t* constr_group;
    uint32_t num_virtual_triangles; // virtual constraints are indexed
    // from num_triangles - num_virtual_triangles
    // to num_triangles-1.
    uint32_t* tri_original_index; // per (real) constraint: its index in the input
    // file's triangle list, before degenerate triangles were dropped. NULL if not
    // tracked. Length = number of real (non-virtual) constraints (input + hole caps).
    uint32_t num_input_triangles; // number of genuine input surface constraints. Extra
    // constraints are appended after these in a fixed order: hole caps, then edge-forcing
    // triangles (2 per input edge), then point-forcing triangles (3 per input point), then
    // the virtual (bounding-box) constraints. UINT32_MAX means "unset" (no caps/edges/points
    // added); callers then treat every real constraint as input surface. It can legitimately
    // be 0 (edges/points but no surface triangles), which is why the sentinel is not 0.
    uint32_t num_edge_triangles; // 2 per inserted edge (see insert_edges_and_points)
    uint32_t num_point_triangles; // 3 per inserted point

    constraints_t()
        : tri_vertices(NULL)
        , num_triangles(0)
        , constr_group(NULL)
        , num_virtual_triangles(0)
        , tri_original_index(NULL)
        , num_input_triangles(UINT32_MAX)
        , num_edge_triangles(0)
        , num_point_triangles(0) {};
    ~constraints_t()
    {
        if (tri_vertices) free(tri_vertices);
        if (constr_group) free(constr_group);
        if (tri_original_index) free(tri_original_index);
    }
};


void fill_half_edges(const constraints_t* constraints, half_edge_t* half_edges);
void sort_half_edges(half_edge_t* half_edges, uint32_t num_half_edges);

// Detect the open boundaries (holes) of the current input constraints and cap
// each one with an ear-clipped triangulation, so that every input triangle bounds
// a closed volume. The cap ("fake") triangles are appended right after the real
// input triangles (num_triangles grows); constraints->num_input_triangles is set
// to the original count so the surface-tracking can exclude the caps. Operates on
// the deduplicated mesh vertices; must be called before the Delaunay permutation
// and before place_virtual_constraints. Returns the number of cap triangles added.
uint32_t fill_holes_in_constraints(constraints_t* constraints, const TetMesh* mesh, bool verbose);

// Extra 1D/0D features to force into the tetrahedralization output.
//   edges:  edge_verts holds n_edge_verts vertices (x,y,z...); edge_idx holds n_edges
//           endpoint index pairs into edge_verts.
//   points: point_verts holds n_points vertices (x,y,z...).
struct extra_features_t
{
    const double* edge_verts = nullptr;
    uint32_t n_edge_verts = 0;
    const uint32_t* edge_idx = nullptr;
    uint32_t n_edges = 0;
    const double* point_verts = nullptr;
    uint32_t n_points = 0;
};

// Force the given edges and points to appear in the output tetrahedralization by
// adding, per edge, the two largest-area triangles among (A, B, A + |AB|*e_{x,y,z})
// (they share edge AB, pinning it), and per point, the three "corner" triangles
// (P, P + L*e_i, P + L*e_j) with L = average input-triangle edge length (they share
// vertex P, pinning it). New vertices are appended to the mesh and the forcing
// triangles to the constraints (num_edge_triangles / num_point_triangles record the
// counts). Must run before the Delaunay permutation and place_virtual_constraints.
//
// keep_corners_inside: the forcing-triangle corners must stay inside the bounding box of the
// vertices already in the mesh -- set when embedding into a given tet mesh, whose vertices are
// all the room there is. A corner that +L along an axis would leave that box is placed at -L
// instead; the triangles still meet only along AB (or only at P), so they pin it the same way.
void insert_edges_and_points(
    constraints_t* constraints,
    TetMesh* mesh,
    const extra_features_t& extra,
    bool verbose,
    bool keep_corners_inside = false);
uint32_t
place_virtual_constraints(TetMesh* mesh, constraints_t* constraints, half_edge_t* half_edges);
void insert_constraints(
    TetMesh*,
    constraints_t*,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**);
} // namespace vol_rem

#endif
