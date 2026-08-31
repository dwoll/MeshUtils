// ----------------------------------------------------------------------- //
// Code adapted from packages
// https://github.com/stla/Boov/
// https://github.com/stla/PolygonSoup/
// https://github.com/stla/cgalMeshes/
// developed and copyright by
// Stéphane Laurent <laurent_step@outlook.fr>
// adapted by
// Daniel Wollschlaeger
// License: GPL-3
// ----------------------------------------------------------------------- //

#ifndef _CGALMESHHEADER_
#include "MeshUtils.h"
#endif

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
bool is_triangle_soup(const std::vector<std::vector<std::size_t>>& polygons) {
    for (const auto& poly : polygons) {
        if (poly.size() != 3) {
            return false;
        }
    }
    return true;
}

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
// currently unused
template <typename MeshInT, PointInT, typename MeshOutT, typename PointOutT>
MeshOutT switch_epeck_epick(const MeshInT &mesh_in) {
  const std::size_t nVertices = mesh_in.number_of_vertices();
  const std::size_t nEdges    = mesh_in.number_of_edges();
  const std::size_t nFaces    = mesh_in.number_of_faces();
  MeshOutT mesh_out;
  mesh_out.reserve(nVertices, nEdges, nFaces);
  for(MeshInT::Vertex_index vd : mesh_in.vertices()) {
    const PointInT vertex = mesh_in.point(vd);
    const double x = CGAL::to_double<KernelInT::FT>(vertex.x());
    const double y = CGAL::to_double<KernelInT::FT>(vertex.y());
    const double z = CGAL::to_double<KernelInT::FT>(vertex.z());
    PointOutT pt(x, y, z);
    mesh_out.add_vertex(pt);
  }
  for(MeshInT::Face_index fd : mesh_in.faces()) {
    std::vector<MeshOutT::Vertex_index> face;
    for(MeshInT::Vertex_index vd :
          vertices_around_face(mesh_in.halfedge(fd), mesh_in)) {
      face.push_back(vd);
    }
    mesh_out.add_face(face);
  }
  return mesh;
}
