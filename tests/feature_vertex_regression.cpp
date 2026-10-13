// Regression: forcing a point or an edge into an embedding whose background mesh already has a
// vertex at a surface vertex's position, and forcing one whose forcing-triangle corners would
// leave the background mesh.
//
// wildmeshing-toolkit's tetwild builds the background mesh from the input vertices plus box
// points, so every surface vertex has a coincident background vertex, and the box is padded by
// diag/15 -- less than a point's forcing displacement on small inputs. The background here is
// built the same way at a small scale: a grid whose lines pass through the cube's corners, with
// a padding of 0.2 against a point displacement of 0.1 x the mean surface edge (0.228).

#include <catch2/catch_test_macros.hpp>

#include "VolumeRemesher/embed.h"

#include <array>
#include <cstdint>
#include <vector>

namespace {

struct Background
{
    std::vector<double> coords;
    std::vector<uint32_t> tets;
};

// A grid with lines at -0.2, 0, 2, 2.2 on each axis, each cell split into the 6 Kuhn tets around
// its main diagonal, every tet oriented so that (v1-v0)x(v2-v0).(v3-v0) > 0.
Background grid_background()
{
    const double xs[4] = {-0.2, 0.0, 2.0, 2.2};
    Background b;
    const auto vid = [](int i, int j, int k) { return uint32_t(i + 4 * j + 16 * k); };
    for (int k = 0; k < 4; k++)
        for (int j = 0; j < 4; j++)
            for (int i = 0; i < 4; i++) {
                b.coords.push_back(xs[i]);
                b.coords.push_back(xs[j]);
                b.coords.push_back(xs[k]);
            }
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int k = 0; k < 3; k++)
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
                for (const auto& p : perms) {
                    int c[3] = {i, j, k};
                    std::array<uint32_t, 4> t;
                    t[0] = vid(c[0], c[1], c[2]);
                    for (int s = 0; s < 3; s++) {
                        c[p[s]]++;
                        t[s + 1] = vid(c[0], c[1], c[2]);
                    }
                    const double* a = &b.coords[3 * t[0]];
                    double e[3][3];
                    for (int s = 0; s < 3; s++)
                        for (int d = 0; d < 3; d++) e[s][d] = b.coords[3 * t[s + 1] + d] - a[d];
                    const double vol = e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1]) -
                                       e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0]) +
                                       e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]);
                    if (vol < 0) std::swap(t[2], t[3]);
                    b.tets.insert(b.tets.end(), t.begin(), t.end());
                }
    return b;
}

// The cube [0,2]^3, outward triangles.
const std::vector<double> cube_coords = {0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0,
                                         0, 0, 2, 2, 0, 2, 2, 2, 2, 0, 2, 2};
const std::vector<uint32_t> cube_tris = {0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4,
                                         2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3};

struct Result
{
    std::vector<vol_rem::bigrational> coords;
    std::vector<std::array<uint32_t, 4>> tets;
    std::vector<std::vector<std::array<uint32_t, 3>>> edge_provenance;
    std::vector<std::array<uint32_t, 2>> point_provenance;
};

Result embed(
    const std::vector<double>& edge_coords,
    const std::vector<uint32_t>& edge_indexes,
    const std::vector<double>& point_coords)
{
    const Background bg = grid_background();
    Result r;
    std::vector<uint32_t> poly_vindexes, cell_findexes, tets_parent, facets_on_input;
    std::vector<bool> cells_with_faces_on_input;
    std::vector<std::vector<uint32_t>> tets_parent_faces;
    std::vector<std::vector<std::array<uint32_t, 4>>> tri_provenance;
    std::vector<uint32_t> tri_group;
    vol_rem::embed_tri_in_poly_mesh(
        cube_coords,
        cube_tris,
        bg.coords,
        bg.tets,
        r.coords,
        poly_vindexes,
        cell_findexes,
        r.tets,
        tets_parent,
        facets_on_input,
        cells_with_faces_on_input,
        tets_parent_faces,
        edge_coords,
        edge_indexes,
        point_coords,
        tri_provenance,
        tri_group,
        r.edge_provenance,
        r.point_provenance,
        /*verbose=*/false);
    return r;
}

// The output vertex v sits exactly at p.
bool vertex_at(const Result& r, uint32_t v, const double* p)
{
    for (int d = 0; d < 3; d++)
        if (!(r.coords[3 * v + d] == vol_rem::bigrational(p[d]))) return false;
    return true;
}

} // namespace

TEST_CASE("a forced point on a surface vertex that the background also has", "[feature_vertex]")
{
    const double p[3] = {0, 0, 0};
    const Result r = embed({}, {}, std::vector<double>(p, p + 3));
    REQUIRE(!r.tets.empty());
    REQUIRE(r.point_provenance.size() == 1);
    CHECK(vertex_at(r, r.point_provenance[0][1], p));
}

TEST_CASE("a forced edge starting at a surface vertex that the background also has", "[feature_vertex]")
{
    const std::vector<double> e = {0, 0, 0, 0.5, 0.5, 0.5};
    const Result r = embed(e, {0, 1}, {});
    REQUIRE(!r.tets.empty());
    REQUIRE(r.edge_provenance.size() == 1);
    CHECK(!r.edge_provenance[0].empty());
}
