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

#include <CGAL/Subdivision_method_3/subdivision_methods_3.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// use EPEC kernel here
// [[Rcpp::export]]
Rcpp::List subdivideCatmullClark_cpp(
  const Rcpp::List rmesh, const unsigned int nIter, const bool triangulate) {
    EMesh3 mesh = makeSurfMesh<EK, EMesh3, EPoint3>(
        rmesh,
        triangulate, // triangulate
        false,       // repair_soup
        false,       // remove_intersections
        1,           // remove_method
        false,       // fill_holes
        false,       // fair hole
        0);          // max_num_holes
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties<EMesh3, EPoint3>(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::CatmullClark_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter));
    mesh.collect_garbage();
    return getRmesh<EK, EMesh3, EPoint3, EVector3>(mesh, triangulate);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// use EPEC kernel here
// [[Rcpp::export]]
Rcpp::List subdivideDooSabin_cpp(
  const Rcpp::List rmesh, const unsigned int nIter, const bool triangulate) {
    EMesh3 mesh = makeSurfMesh<EK, EMesh3, EPoint3>(
        rmesh,
        triangulate, // triangulate
        false,       // repair_soup
        false,       // remove_intersections
        1,           // remove_method
        false,       // fill_holes
        false,       // fair hole
        0);          // max_num_holes
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties<EMesh3, EPoint3>(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::DooSabin_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter));
    mesh.collect_garbage();
    return getRmesh<EK, EMesh3, EPoint3, EVector3>(mesh, triangulate);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// use EPEC kernel here
// [[Rcpp::export]]
Rcpp::List subdivideSqrt3_cpp(
  const Rcpp::List rmesh, const unsigned int nIter, const bool triangulate) {
    EMesh3 mesh = makeSurfMesh<EK, EMesh3, EPoint3>(
        rmesh,
        triangulate, // triangulate
        false,       // repair_soup
        false,       // remove_intersections
        1,           // remove_method
        false,       // fill_holes
        false,       // fair hole
        0);          // max_num_holes
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties<EMesh3, EPoint3>(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::Sqrt3_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter));
    mesh.collect_garbage();
    return getRmesh<EK, EMesh3, EPoint3, EVector3>(mesh, triangulate);
}
