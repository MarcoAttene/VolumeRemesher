#include "embed.h"
#include "BSP.h"

#include <cassert>
#include <vector>

namespace vol_rem {
void embed_tri_in_poly_mesh(
    const std::vector<double> &tri_vrt_coords,
    const std::vector<uint32_t> &triangle_indexes,
    const std::vector<double> &tet_vrt_coords,
    const std::vector<uint32_t> &tet_indexes,
    std::vector<bigrational> &out_vrt_coords,
    std::vector<uint32_t> &out_poly_vindexes,
    std::vector<uint32_t> &out_cell_findexes,
    std::vector<std::array<uint32_t, 4>> &out_tets,
    std::vector<uint32_t> &final_tets_parent,
    std::vector<uint32_t> &facets_on_input,
    std::vector<bool> &cells_with_faces_on_input,
    std::vector<std::vector<uint32_t>> &final_tets_parent_faces, bool verbose) {
  // Make a conformal polyhedralization
  BSPcomplex *complex = remakePolyhedralMesh(
      tri_vrt_coords.data(), (uint32_t)tri_vrt_coords.size() / 3,
      triangle_indexes.data(), (uint32_t)triangle_indexes.size() / 3,
      tet_vrt_coords.data(), (uint32_t)tet_vrt_coords.size() / 3,
      tet_indexes.data(), (uint32_t)tet_indexes.size() / 4, verbose, true);

  for (size_t f_id = 0; f_id < complex->faces.size(); f_id++)
    complex->triangulateFace(f_id);

  complex->makeTetrahedra(verbose);

  if (verbose)
    printf("Producing vertices...\n");
  // Get exact vertex coordinates
  out_vrt_coords.resize(complex->vertices.size() * 3);
  for (uint64_t v_id = 0; v_id < complex->vertices.size(); v_id++) {
    if (!complex->vertices[v_id]->getExactXYZCoordinates(
            out_vrt_coords[v_id * 3], out_vrt_coords[v_id * 3 + 1],
            out_vrt_coords[v_id * 3 + 2]))
      ip_error("embed_tri_in_poly_mesh: could not compute exact coordinates. "
               "Should not happen!\n");
  }

  if (verbose)
    printf("Producing facets...\n");
  // Get facets
  cells_with_faces_on_input.resize(complex->cells.size(), false);
  for (size_t f_id = 0; f_id < complex->faces.size(); f_id++) {
    BSPface &face = complex->faces[f_id];
    std::vector<uint32_t> face_vrts(face.edges.size(), 0);
    complex->list_faceVertices(face, face_vrts);
    out_poly_vindexes.push_back((uint32_t)face_vrts.size());
    for (uint32_t cvi : face_vrts)
      out_poly_vindexes.push_back(cvi);
    if (face.colour == BLACK_A) {
      facets_on_input.push_back((uint32_t)f_id);
      if (face.conn_cells[0] >= 0)
        cells_with_faces_on_input[face.conn_cells[0]] = true;
      if (face.conn_cells[1] >= 0)
        cells_with_faces_on_input[face.conn_cells[1]] = true;
    }
  }

  if (verbose)
    printf("Producing cells...\n");
  // Get polyhedra
  for (uint64_t c_id = 0; c_id < complex->cells.size(); c_id++) {
    BSPcell &cell = complex->cells[c_id];
    out_cell_findexes.push_back((uint32_t)cell.faces.size());
    for (uint64_t cfi : cell.faces)
      out_cell_findexes.push_back((uint32_t)cfi);
  }
  if (verbose)
    printf("Producing tets...\n");
  // Get polyhedra
  final_tets_parent = complex->final_tets_parent;
  final_tets_parent_faces = complex->final_tets_parent_faces;
  assert(complex->final_tets.size() % 4 == 0);
  out_tets.resize(complex->final_tets.size() / 4);
  for (uint64_t t_id = 0; t_id < complex->final_tets.size(); t_id += 4) {
    for (uint64_t i = 0; i < 4; ++i)
      out_tets[t_id / 4][i] = complex->final_tets[t_id + i];
  }
  if (verbose)
    printf("Done\n");
}
} // namespace vol_rem
