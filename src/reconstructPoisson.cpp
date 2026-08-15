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

// [[Rcpp::export]]
Rcpp::List reconstructPoisson_cpp(const Rcpp::NumericMatrix pts,
                                  const Rcpp::NumericMatrix normals,
                                  const double spacing,
                                  const double smAngle,
                                  const double smRadius,
                                  const double smDistance,
                                  const bool clean) {
  const size_t nPts = pts.ncol();
  std::vector<P3wn> points(nPts);
  for(size_t i = 0; i < nPts; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    const Rcpp::NumericVector nrml_i = normals(Rcpp::_, i);
    points[i] =
      std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)),
                     Vector3(nrml_i(0), nrml_i(1), nrml_i(2)));
  }

  double spacingUse;
  if(spacing == -1.0) {
    spacingUse = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
      points, 6, /* knn = 1 ring */
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>()));
  } else {
      spacingUse = spacing;
  }

  Polyhedron poly;
  const bool success = CGAL::poisson_surface_reconstruction_delaunay(
    points.begin(), points.end(),
    CGAL::First_of_pair_property_map<P3wn>(),
    CGAL::Second_of_pair_property_map<P3wn>(),
    poly,
    spacingUse, smAngle, smRadius, smDistance);

  if(!success) {
    Rcpp::stop("Poisson surface reconstruction has failed.");
  }

  size_t id = 1;
  for(Polyhedron::Vertex_iterator vit = poly.vertices_begin();
      vit != poly.vertices_end(); ++vit) {
    vit->id() = id;
    id++;
  }

  bool isTriangle = PMP::triangulate_faces(poly);
  if(!isTriangle) {
    Rcpp::stop("Could not triangluate polyhedron.");
  }

  Mesh3 mesh;
  CGAL::copy_face_graph(poly, mesh);
  Rcpp::List rmeshOut = getRmesh(mesh);
  return rmeshOut;
}
