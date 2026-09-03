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

#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List smoothShape_cpp(
  const Rcpp::List rmesh,
  const Rcpp::IntegerVector indices,
  const unsigned int nIter,
  const double time,
  const bool normals) {
    Mesh3 mesh = make_surf_mesh<K, Mesh3, Point3>(
        rmesh,
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false,       // remove_intersections
        1,           // remove_method
        false,       // fill_holes
        false,       // fair hole
        0);          // max_num_holes
    std::set<Mesh3::Vertex_index> constrained_vertices;
    for(Mesh3::Vertex_index v : vertices(mesh)) {
      if(is_border(v, mesh)) {
          constrained_vertices.insert(v);
      }
    }
    CGAL::Boolean_property_map<std::set<Mesh3::Vertex_index>> vcmap(constrained_vertices);
    const size_t nIdx = indices.size();
    if(nIdx == 0) {
        PMP::smooth_shape(mesh, time,
                          CGAL::parameters::number_of_iterations(nIter)
                          .vertex_is_constrained_map(vcmap));
    } else {
        std::list<fc_dscrptr> selectedFaces;
        const size_t nFaces = mesh.number_of_faces();
        for(std::size_t i = 0; i < nIdx; i++) {
          const size_t idx = indices(i);
          if(idx >= nFaces) {
            Rcpp::stop("Face index too large.");
          }
          // TODO use `CGAL::SM_Face_index(idx)`
          selectedFaces.push_back(*(mesh.faces().begin() + idx));
        }
        PMP::smooth_shape<Mesh3>(selectedFaces, mesh, time,
                                 PMP::parameters::number_of_iterations(nIter)
                                 .vertex_is_constrained_map(vcmap));
    }
    // subdivision requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}
