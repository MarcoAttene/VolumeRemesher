#include "BSP.h"
#include <vector>

/// <summary>
/// embed_tri_in_poly_mesh
/// Takes a set of triangles T and a tetrahedral mesh M.
/// Cuts tetrahedra in M so as to form a polyhedral mesh M' such that:
/// 1) each triangle in T is represented by the union of facets in M';
/// 2) each polyhedral cell in M' is (possibly weakly) convex.
/// Since cut points may be not representable using floating point
/// coordinates, M' is returned with rational coordinates.
/// </summary>
///
/// INPUT (T and M)
/// <param name="tri_vrt_coords">Vertex coordinates of T
/// (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param> <param
/// name="triangle_indexes">Triangle indexes of T (t1_v1, t1_v2, t1_v3, t2_v1,
/// t2_v2, t2_v3, ..., tn_v1, tn_v2, tn_v3)</param> <param
/// name="tet_vrt_coords">Vertex coordinates of M
/// (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param> <param name="tet_indexes">Tet
/// indexes of M (t1_v1, t1_v2, t1_v3, t1_v4, t2_v1, t2_v2, t2_v3, t2_v4,
/// ...)</param>
///
/// OUTPUT (M')
/// <param name="vertices">Vertex coordinates of M'
/// (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param> <param name="facets">Polygonal
/// facets in M'. A facet with n vertices is a sequence (n, p_v1, p_v2, ...,
/// p_vn). The first number in the sequence is the number of vertices in the
/// polygon, whereas the subsequent n numbers are its vertex indexes.</param>
/// <param name="cells">Polyhedral cells in M'. A polyhedron with n facets is a
/// sequence (n, p_f1, p_f2, ..., p_fn). The first number in the sequence is the
/// number of facets bounding the polyhedron, whereas the subsequent n numbers
/// are its facet indexes.</param> <param name="facets_on_input">Indexes of
/// facets that overlap with T.</param>
///
/// <param name="verbose">Set to TRUE to enable verbosity.</param>

void embed_tri_in_poly_mesh(const std::vector<double> &tri_vrt_coords,
                            const std::vector<uint32_t> &triangle_indexes,
                            const std::vector<double> &tet_vrt_coords,
                            const std::vector<uint32_t> &tet_indexes,
                            std::vector<bigrational> &vertices,
                            std::vector<uint32_t> &facets,
                            std::vector<uint32_t> &cells,
                            std::vector<uint32_t> &facets_on_input,
                            bool verbose);

//
// TO DO: CONSIDER INPUT POINTS AND SEGMENTS
