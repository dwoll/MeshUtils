// ----------------------------------------------------------------------- //
// Code adapted from packages
// https://github.com/stla/Boov/
// https://github.com/stla/PolygonSoup/
// https://github.com/stla/cgalMeshes/
// developed and copyright by
// Stéphane Laurent <laurent_step@outlook.fr>
// License: GPL-3
// ----------------------------------------------------------------------- //

#ifndef _CGALMESHHEADER_
#include "MeshUtils.h"
#endif

// [[Rcpp::export]]
Rcpp::List reconstructPoisson_cpp(Rcpp::NumericMatrix pts,
                                  Rcpp::NumericMatrix normals,
                                  double spacing,
                                  double sm_angle,
                                  double sm_radius,
                                  double sm_distance) {
  const size_t npoints = pts.ncol();
  std::vector<P3wn> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    const Rcpp::NumericVector nrml_i = normals(Rcpp::_, i);
    points[i] =
      std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)),
                     Vector3(nrml_i(0), nrml_i(1), nrml_i(2)));
  }

  if(spacing == -1.0) {
    spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
      points, 6,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>()));
  }

  Polyhedron mesh;
  const bool psr = CGAL::poisson_surface_reconstruction_delaunay(
    points.begin(), points.end(), CGAL::First_of_pair_property_map<P3wn>(),
    CGAL::Second_of_pair_property_map<P3wn>(), mesh, spacing, sm_angle,
    sm_radius, sm_distance);

  if(!psr) {
    throw Rcpp::exception("Poisson surface reconstruction has failed.");
  }

  int id = 1;
  for(Polyhedron::Vertex_iterator vit = mesh.vertices_begin();
      vit != mesh.vertices_end(); ++vit) {
    vit->id() = id;
    id++;
  }

  const size_t nfacets = mesh.size_of_facets();
  const size_t nvertices = mesh.size_of_vertices();

  Rcpp::IntegerMatrix facets(3, nfacets);
  {
    size_t i = 0;
    for(Polyhedron::Facet_iterator fit = mesh.facets_begin();
        fit != mesh.facets_end(); fit++) {
      Rcpp::IntegerVector facet_i(3);
      facet_i(0) = fit->halfedge()->vertex()->id();
      facet_i(1) = fit->halfedge()->next()->vertex()->id();
      facet_i(2) = fit->halfedge()->opposite()->vertex()->id();
      facets(Rcpp::_, i) = facet_i;
      i++;
    }
  }

  Rcpp::NumericMatrix vertices(3, nvertices);
  {
    size_t i = 0;
    for(Polyhedron::Vertex_iterator vit = mesh.vertices_begin();
        vit != mesh.vertices_end(); vit++) {
      Rcpp::NumericVector vertex_i(3);
      vertex_i(0) = vit->point().x();
      vertex_i(1) = vit->point().y();
      vertex_i(2) = vit->point().z();
      vertices(Rcpp::_, i) = vertex_i;
      i++;
    }
  }

  // std::ofstream("out.off") << std::setprecision(17) << mesh;

  return Rcpp::List::create(Rcpp::Named("vertices") = vertices,
                            Rcpp::Named("faces")    = facets, // TODO why "facets"?
                            Rcpp::Named("spacing")  = spacing);
}
