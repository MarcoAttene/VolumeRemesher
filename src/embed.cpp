#include "embed.h"
#include <BSP.h>
#include <vector>

/// <summary>
/// embed_tri_in_poly_mesh
/// Takes a set of triangles T and a tetrahedral mesh M.
/// Cuts tetrahedra in M so as to form a polyhedral mesh M' such that
/// each triangle in T is represented by the union of facets in M'.
/// Since cutting points may be not representable using floating point 
/// coordinates, M' is returned with rational coordinates.
/// </summary>
/// <param name="tri_vrt_coords">Vertex coordinates of T (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param>
/// <param name="triangle_indexes">Triangle indexes of T (t1_v1, t1_v2, t1_v3, t2_v1, t2_v2, t2_v3, ..., tn_v1, tn_v2, tn_v3)</param>
/// <param name="tet_vrt_coords">Vertex coordinates of M (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param>
/// <param name="tet_indexes">Tet indexes of M (t1_v1, t1_v2, t1_v3, t1_v4, t2_v1, t2_v2, t2_v3, t2_v4, ...)</param>
/// <param name="out_vrt_coords">Vertex coordinates of M' (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param>
/// <param name="out_poly_indexes">Polyhedra indexes of M'. A polyhedron with n vertices is a sequence (n, p_v1, p_v2, ..., p_vn).
/// The first number in the sequence is the number of vertices in the polyhedron, whereas the subsequent n numbers are its
/// vertex indexes.</param>
/// <param name="verbose">Set to TRUE to enable verbosity.</param>

void embed_tri_in_poly_mesh(
	const std::vector<double>& tri_vrt_coords,
	const std::vector<uint32_t>& triangle_indexes,
	const std::vector<double>& tet_vrt_coords,
	const std::vector<uint32_t>& tet_indexes,
	std::vector<bigrational>& out_vrt_coords,
	std::vector<uint32_t>& out_poly_vindexes,
	std::vector<uint32_t>& out_cell_findexes,
	std::vector<uint32_t>& facets_on_input,
	bool verbose
	)
{
	// Make a conformal polyhedralization
	BSPcomplex* complex = remakePolyhedralMesh(
		tri_vrt_coords.data(), (uint32_t)tri_vrt_coords.size(), triangle_indexes.data(), (uint32_t)triangle_indexes.size(),
		tet_vrt_coords.data(), (uint32_t)tet_vrt_coords.size(), tet_indexes.data(), (uint32_t)tet_indexes.size(), verbose,
		true);

	// Get exact vertex coordinates
	out_vrt_coords.reserve(complex->vertices.size() * 3);
	for (uint64_t v_id = 0; v_id < complex->vertices.size(); v_id++) {
		bigrational x, y, z;
		if (complex->vertices[v_id]->getExactXYZCoordinates(x, y, z)) {
			out_vrt_coords.push_back(x);
			out_vrt_coords.push_back(y);
			out_vrt_coords.push_back(z);
		}
		else ip_error("embed_tri_in_poly_mesh: could not compute exact coordinates. Should not happen!\n");
	}

	// Get facets
	for (size_t f_id = 0; f_id < complex->faces.size(); f_id++) {
		BSPface& face = complex->faces[f_id];
		std::vector<uint32_t> face_vrts;
		complex->list_faceVertices(face, face_vrts);
		out_poly_vindexes.push_back((uint32_t)face_vrts.size());
		for (uint32_t cvi : face_vrts) out_poly_vindexes.push_back(cvi);
		if (face.colour == BLACK_A) facets_on_input.push_back((uint32_t)f_id);
	}

	// Get polyhedra
	for (uint64_t c_id = 0; c_id < complex->cells.size(); c_id++) {
		BSPcell& cell = complex->cells[c_id];
		out_cell_findexes.push_back((uint32_t)cell.faces.size());
		for (uint64_t cfi : cell.faces) out_cell_findexes.push_back((uint32_t)cfi);
	}
}
