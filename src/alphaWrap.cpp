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
// [[Rcpp::export]]
Rcpp::List alphaWrap_cpp(
    const Rcpp::NumericMatrix pts, const double alphaRel, const double offsetRel) {
  std::vector<Point3> points = matrix_to_points3<Point3>(pts);
  CGAL::Bbox_3 bbox = CGAL::bbox_3(std::cbegin(points), std::cend(points));
  const double diag_len = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                    CGAL::square(bbox.ymax() - bbox.ymin()) +
                                    CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha  = diag_len / alphaRel;
  const double offset = diag_len / offsetRel;
  Mesh3 meshWrap;
  CGAL::alpha_wrap_3(points, alpha, offset, meshWrap);
  return getRmesh(meshWrap, false);
}
