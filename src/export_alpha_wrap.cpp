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

#include <CGAL/alpha_wrap_3.h>

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List alphaWrapPoints_cpp(
    const Rcpp::NumericMatrix pts,
    const double alpha_rel,
    const double offset_rel,
    const bool normals) {
  std::vector<Point3> points = matrix_to_points3<Point3>(pts);
  CGAL::Bbox_3 bbox = CGAL::bbox_3(std::cbegin(points), std::cend(points));
  const double diag_len = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                    CGAL::square(bbox.ymax() - bbox.ymin()) +
                                    CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha  = diag_len / alpha_rel;
  const double offset = diag_len / offset_rel;
  Mesh3 mesh_wrap;
  CGAL::alpha_wrap_3(points, alpha, offset, mesh_wrap);
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh_wrap, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List alphaWrapMesh_cpp(
    const Rcpp::List rmesh,
    const double alpha_rel,
    const double offset_rel,
    const bool normals) {
  Mesh3 mesh = make_surf_mesh<K, Mesh3, Point3>(
      rmesh,
      true,           // triangulate - must be triangle
      false,          // repair_soup
      false,          // remove_intersections
      1,              // remove_method
      false,          // fill_holes
      false,          // fair hole
      0);             // max_num_holes

  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  const double diag_len = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                    CGAL::square(bbox.ymax() - bbox.ymin()) +
                                    CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha  = diag_len / alpha_rel;
  const double offset = diag_len / offset_rel;
  Mesh3 mesh_wrap;
  CGAL::alpha_wrap_3(mesh, alpha, offset, mesh_wrap);
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh_wrap, false, normals);
}
