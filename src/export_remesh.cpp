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

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshIsotropicUniform_cpp(
    const Rcpp::List rmesh,
    const double targetEdgeLen,
    const unsigned int nIter,
    const unsigned int nRelaxSteps,
    const bool normals) {
    Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
        rmesh,
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false,       // remove_intersections
        1,           // remove_method
        false,       // fill_holes
        false,       // fair hole
        0);          // max_num_holes
    std::vector<hlfdg_dscrptr> borderHalfEdges;
    // PMP::border_halfedges(faces(mesh), mesh, std::back_inserter(borderHalfEdges));
    CGAL::border_halfedges(faces(mesh), mesh, std::back_inserter(borderHalfEdges));  // requires CGAL 6.2
    std::vector<dg_dscrptr> border;
    std::size_t nheBorder = borderHalfEdges.size();
    border.reserve(nheBorder);
    // for(std::size_t i = 0; i < nheBorder; i++) {
    for(std::size_t i : borderHalfEdges) {
      border.emplace_back(mesh.edge(borderHalfEdges[i]));
    }
    PMP::split_long_edges(border, targetEdgeLen, mesh);
    PMP::Uniform_sizing_field<Mesh3> sizing_field(targetEdgeLen, mesh);
    PMP::isotropic_remeshing(
      faces(mesh),
      sizing_field,
      mesh,
      PMP::parameters::number_of_iterations(nIter)
                      .number_of_relaxation_steps(nRelaxSteps)
                      .protect_constraints(true));
    mesh.collect_garbage();
    // remeshing requires triangle mesh -> output is triangle
    return getRmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshIsotropicAdapt_cpp(
    const Rcpp::List rmesh,
    const double tol,
    const double edgeMin,
    const double edgeMax,
    const unsigned int nIter,
    const unsigned int nRelaxSteps,
    const bool normals) {
    Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
        rmesh,
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false,       // remove_intersections
        1,           // remove_method
        false,       // fill_holes
        false,       // fair hole
        0);          // max_num_holes
    const std::pair edge_min_max{ edgeMin, edgeMax };
    PMP::Adaptive_sizing_field<Mesh3> sizing_field(
        tol,
        edge_min_max,
        faces(mesh),
        mesh);
    PMP::isotropic_remeshing(
      faces(mesh),
      sizing_field,
      mesh,
      PMP::parameters::number_of_iterations(nIter)
                      .number_of_relaxation_steps(nRelaxSteps)
                      .protect_constraints(true));
    mesh.collect_garbage();
    // remeshing requires triangle mesh -> output is triangle
    return getRmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}
