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

#include <CGAL/Polygon_mesh_processing/distance.h>

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getHausdorffApprox_cpp(
    const Rcpp::List rmesh1,
    const Rcpp::List rmesh2,
    const bool symmetric,
    const unsigned int n) {
  Mesh3 mesh1 = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh1,
      true,        // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
  Mesh3 mesh2 = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh2,
      true,        // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
  if(CGAL::is_empty(mesh1)) {
    Rcpp::warning("Mesh 1 is empty.");
    return Rcpp::NumericVector::get_na();
  }
  if(CGAL::is_empty(mesh2)) {
    Rcpp::warning("Mesh 2 is empty.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_triangle_mesh(mesh1)) {
    Rcpp::warning("Mesh 1 is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_triangle_mesh(mesh2)) {
    Rcpp::warning("Mesh 2 is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  double d;
  // PMP::parameters::use_random_uniform_sampling(true)     // true
  // PMP::parameters::do_sample_vertices(true)              // true
  // PMP::parameters::do_sample_edges(true)                 // true
  // PMP::parameters::do_sample_faces(true)                 // true
  //
  // PMP::parameters::number_of_points_on_faces(n)          // unsigned int *
  // PMP::parameters::number_of_points_on_edges(n)          // unsigned int
  //
  // PMP::parameters::number_of_points_per_distance_unit(n) // double
  // PMP::parameters::number_of_points_per_edge(n)          // unsigned int
  // PMP::parameters::number_of_points_per_area_unit(n)     // double
  // PMP::parameters::number_of_points_per_face(n)          // unsigned int
  if(symmetric) {
    if(n > 0) {
        d = CGAL::to_double<K::FT>(PMP::approximate_symmetric_Hausdorff_distance<PIA_TAG>(
          mesh1, mesh2, PMP::parameters::number_of_points_on_faces(n)));
    } else {
        d = CGAL::to_double<K::FT>(PMP::approximate_symmetric_Hausdorff_distance<PIA_TAG>(
          mesh1, mesh2));
    }
  } else {
    if(n > 0) {
        d = CGAL::to_double<K::FT>(PMP::approximate_Hausdorff_distance<PIA_TAG>(
          mesh1, mesh2, PMP::parameters::number_of_points_on_faces(n)));
    } else {
        d = CGAL::to_double<K::FT>(PMP::approximate_Hausdorff_distance<PIA_TAG>(
          mesh1, mesh2));
    }
  }
  return d;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getHausdorffEst_cpp(
    const Rcpp::List rmesh1,
    const Rcpp::List rmesh2,
    const bool symmetric,
    const double error_bound) {
    Mesh3 mesh1 = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh1,
        true,        // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    Mesh3 mesh2 = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh2,
        true,        // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    if(CGAL::is_empty(mesh1)) {
      Rcpp::warning("Mesh 1 is empty.");
      return Rcpp::NumericVector::get_na();
    }
    if(CGAL::is_empty(mesh2)) {
      Rcpp::warning("Mesh 2 is empty.");
      return Rcpp::NumericVector::get_na();
    }
    if(!CGAL::is_triangle_mesh(mesh1)) {
      Rcpp::warning("Mesh 1 is not triangle.");
      return Rcpp::NumericVector::get_na();
    }
    if(!CGAL::is_triangle_mesh(mesh2)) {
      Rcpp::warning("Mesh 2 is not triangle.");
      return Rcpp::NumericVector::get_na();
    }
    double d;
    if(symmetric) {
        d = CGAL::to_double<K::FT>(PMP::bounded_error_symmetric_Hausdorff_distance<PIA_TAG>(
            mesh1, mesh2, error_bound));
    } else {
        d = CGAL::to_double<K::FT>(PMP::bounded_error_Hausdorff_distance<PIA_TAG>(
            mesh1, mesh2, error_bound));
    }
    return d;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List getMetro_cpp(
    const Rcpp::List rmesh1,
    const Rcpp::List rmesh2,
    const bool symmetric,
    const double p,
    const unsigned int n) {
  Mesh3 mesh1 = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh1,
      true,        // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
  Mesh3 mesh2 = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh2,
      true,        // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
  if(CGAL::is_empty(mesh1)) {
    Rcpp::warning("Mesh 1 is empty.");
    return Rcpp::NumericVector::get_na();
  }
  if(CGAL::is_empty(mesh2)) {
    Rcpp::warning("Mesh 2 is empty.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_triangle_mesh(mesh1)) {
    Rcpp::warning("Mesh 1 is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_triangle_mesh(mesh2)) {
    Rcpp::warning("Mesh 2 is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  std::tuple<double, double, double> metro = getMetro<K, Mesh3, Point3>(
    mesh1,
    mesh2,
    symmetric,
    p,
    n);
  double HDq  = get<0>(metro);
  double ASSD = get<1>(metro);
  double RMSE = get<2>(metro);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("HDq")  = HDq,
                                      Rcpp::Named("ASSD") = ASSD,
                                      Rcpp::Named("RMSE") = RMSE);
  return out;
}
