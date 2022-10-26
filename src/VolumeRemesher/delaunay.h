#ifndef _DELAUNAY_
#define _DELAUNAY_

#include <VolumeRemesher/implicit_point.h>
#include <cstring>

int orient2d(double p1x, double p1y, double p2x, double p2y, double p3x,
             double p3y);
int orient3d(double px, double py, double pz, double qx, double qy, double qz,
             double rx, double ry, double rz, double sx, double sy, double sz);
int inSphere(const double pax, const double pay, const double paz,
             const double pbx, const double pby, const double pbz,
             const double pcx, const double pcy, const double pcz,
             const double pdx, const double pdy, const double pdz,
             const double pex, const double pey, const double pez);

inline double orient2d(const double *p1, const double *p2, const double *p3) {
  return -orient2d(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
}

inline double orient3d(const double *p, const double *q, const double *r,
                       const double *s) {
  return -orient3d(p[0], p[1], p[2], q[0], q[1], q[2], r[0], r[1], r[2], s[0],
                   s[1], s[2]);
}

inline double insphere(const double *pa, const double *pb, const double *pc,
                       const double *pd, const double *pe) {
  return -inSphere(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0], pc[1],
                   pc[2], pd[0], pd[1], pd[2], pe[0], pe[1], pe[2]);
}

#pragma intrinsic(fabs)

#define INFINITE_VERTEX UINT32_MAX

// Vertex type

struct vertex_t {
  double coord[3];         // Coordinates
  uint64_t inc_tet;        // Incident tetrahedron
  uint32_t original_index; // Index to support reordering
};

// Tetrahedral mesh data structure

class TetMesh {
public:
  // Vertices
  uint32_t num_vertices; // Total number of vertices
  vertex_t *vertices;    // Vertex array

  // Tetrahedra
  uint32_t tet_num_vertices; // Number of vertices that belong to tetrahedra
  uint64_t tet_num;          // Number of tetrahedra
  uint64_t tet_size;         // Current capacity of the tetrahedron array
  uint32_t *tet_node;        // Tetrahedron array

  // Additional information
  uint64_t *tet_neigh;       // Adjacent tetrahedron array
  double *tet_subdet;        // Precomputed sub-determinants for insphere
  double o3d_static_filter;  // Static filter for orient3d
  double isp_static_filter;  // Static filter for insphere
  uint32_t *mark_tetrahedra; // General purpose tetrahedron marks

  // Constructor and destructor
  TetMesh();
  ~TetMesh();

  // Init the mesh with a tet connecting four non coplanar points in vertices
  void init();

  // Init the mesh from two arrays defining vertices and tets
  void init(const double *vertex_coords, size_t _num_vertices,
            const uint32_t *tet_indices, size_t _num_tets);

  // Create a Delaunay tetrahedrization by incremental insertion
  void tetrahedrize();

  // Save the mesh to a .tet file
  void saveTET(const char *filename);

  // Return an array containing incident tetrahedra at a given vertex v.
  // Store the array length in numtets
  uint64_t *incident_tetrahedra(const uint32_t v, uint64_t *numtets);

  // Return the i'th tet in neighbors 'n'
  inline uint64_t getIthNeighbor(const uint64_t *n, const uint64_t i) const {
    return n[i] & 0xFFFFFFFFFFFFFFFC;
  }

  // Return the i'th tet adjacent to 't'
  inline uint64_t getNeighbor(const uint64_t t, const uint64_t i) const {
    return getIthNeighbor(tet_neigh + t, i);
  }

  // Return an array containing incident tetrahedra at a given edge
  // edge_ends[2]. The first element of the array is first_tet_ind, which is
  // assumed to be incident at the edge. Store the array length in numtets
  uint64_t *ETrelation(const uint32_t *edge_ends, const uint64_t first_tet_ind,
                       uint64_t *numtets) const;

protected:
  struct DelTmp {
    uint32_t node[4];
    uint64_t bnd;
  } * Del_tmp;

  uint64_t Del_num_tmp;
  uint64_t Del_size_tmp;
  uint64_t *Del_deleted;
  uint64_t Del_num_deleted;
  uint64_t Del_size_deleted;
  uint32_t *Del_buffer;

  void allocTmpStruct(uint32_t num_v);
  void releaseTmpStruct();
  void bnd_push(uint32_t vta, uint32_t node1, uint32_t node2, uint32_t node3,
                uint64_t bnd);

  uint64_t searchTetrahedron(uint64_t tet, const uint32_t v_id);
  void deleteInSphereTets(uint64_t tet, const uint32_t v_id);
  void tetrahedrizeHole(uint64_t *tet);
  void removeDelTets();
  void compute_subDet(const uint64_t tet);

  double vertexInTetSphere(uint64_t tet, uint32_t v_id);

  // Pre-allocate memory to store tetrahedra
  void reserve(uint64_t numtet);

  // Inserts a vertex which is already in the vertices array
  // vi is the vertex index in 'vertices', ct is the previously inserted
  // tetrahedron
  void insertExistingVertex(const uint32_t vi, uint64_t &ct);

  // Same as above, but this simply splits tetrahedra without making them
  // Delaunay
  void insertExistingVertexNonDelaunay(const uint32_t vi, uint64_t &ct);
  void deleteInVertexTets(uint64_t tet, const uint32_t v_id);
  bool vertexNotCoplanar(uint32_t f1, uint32_t f2, uint32_t f3,
                         uint32_t v_id) const;

public:
  bool loadTET(const char *filename);
  bool loadMEDIT(const char *filename);
  void insertExistingVerticesNonDelaunay(const uint32_t first_v);

  void checkMesh() const {
    size_t i;
    // Check tet nodes
    for (i = 0; i < tet_num; i++) {
      const uint32_t *tn = tet_node + i * 4;
      if (tn[0] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[1] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[2] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[3] != INFINITE_VERTEX && tet_node[i * 4 + 3] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[3] != INFINITE_VERTEX &&
          orient3d(vertices[tn[0]].coord, vertices[tn[1]].coord,
                   vertices[tn[2]].coord, vertices[tn[3]].coord) <= 0)
        ip_error("Inverted/degn tet\n");
      if (tn[0] == tn[1] || tn[0] == tn[2] || tn[1] == tn[2])
        ip_error("Wrong tet node indexes!\n");
    }
    // Check neighbors
    for (i = 0; i < tet_num * 4; i++)
      if (tet_neigh[tet_neigh[i]] != i)
        ip_error("Wrong neighbor!\n");
    // Check vt*
    for (i = 0; i < num_vertices; i++) {
      const uint32_t *tn = tet_node + vertices[i].inc_tet * 4;
      if (tn[0] != i && tn[1] != i && tn[2] != i && tn[3] != i)
        ip_error("Wrong vt*!\n");
    }
  }
};

#endif // _DELAUNAY_
