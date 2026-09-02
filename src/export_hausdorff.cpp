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
double getHausdorffApprox_cpp(
    const Rcpp::List rmesh1,
    const Rcpp::List rmesh2,
    const bool symmetric) {
  Mesh3 mesh1 = makeSurfMesh<K, Mesh3, Point3>(
      rmesh1,
      true,         // triangulate
      false,        // repair_soup
      false,        // remove_intersections
      1,            // remove_method
      false,        // fill_holes
      false,        // fair hole
      0);           // max_num_holes
  Mesh3 mesh2 = makeSurfMesh<K, Mesh3, Point3>(
      rmesh2,
      true,         // triangulate
      false,        // repair_soup
      false,        // remove_intersections
      1,            // remove_method
      false,        // fill_holes
      false,        // fair hole
      0);           // max_num_holes
  if(CGAL::is_empty(mesh1)) {
    Message("Mesh 1 is empty.");
    return Rcpp::NumericVector::get_na();
  }
  if(CGAL::is_empty(mesh2)) {
    Message("Mesh 2 is empty.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_triangle_mesh(mesh1)) {
    Message("Mesh 1 is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_triangle_mesh(mesh2)) {
    Message("Mesh 2 is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  double d;
  if(symmetric) {
    d = CGAL::to_double<K::FT>(PMP::approximate_symmetric_Hausdorff_distance<PIA_TAG>(mesh1, mesh2));
  } else {
    d = CGAL::to_double<K::FT>(PMP::approximate_Hausdorff_distance<PIA_TAG>(mesh1, mesh2));
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
    Mesh3 mesh1 = makeSurfMesh<K, Mesh3, Point3>(
        rmesh1,
        true,         // triangulate
        false,        // repair_soup
        false,        // remove_intersections
        1,            // remove_method
        false,        // fill_holes
        false,        // fair hole
        0);           // max_num_holes
    Mesh3 mesh2 = makeSurfMesh<K, Mesh3, Point3>(
        rmesh2,
        true,         // triangulate
        false,        // repair_soup
        false,        // remove_intersections
        1,            // remove_method
        false,        // fill_holes
        false,        // fair hole
        0);           // max_num_holes
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
