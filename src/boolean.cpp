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
template <typename MeshT>
void checkMesh1(const MeshT &mesh, std::size_t i) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "Mesh n" + std::to_string(i) + " self-intersects.";
    Rcpp::stop(msg);
  }
}

// ----------------------------------------------------------------------- //
template <typename MeshT>
void checkMesh2(const MeshT &mesh, const std::string& what) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "The " + what + " self-intersects.";
    Rcpp::stop(msg);
  }
}

// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
MeshT boolIntersection(const Rcpp::List &rmeshes,
                       const Rcpp::LogicalVector triangulate,
                       const bool repairSoup) {
  const std::size_t nMeshes = rmeshes.size();
  std::vector<MeshT> meshes(nMeshes);
  Rcpp::List rmesh_0 = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh1");
  MeshT mesh_0 = makeSurfMesh<KernelT, MeshT, PointT>(
      rmesh_0,
      triangulate[0], // triangulate
      repairSoup,     // repair_soup
      false,          // remove_intersections
      1,              // remove_method
      false,          // fill_holes
      false,          // fair hole
      0);             // max_num_holes

  meshes[0] = mesh_0;
  for(std::size_t i = 1; i < nMeshes; i++) {
    if(i == 1) {
      checkMesh1<MeshT>(meshes[0], 1);
    } else {
      checkMesh2<MeshT>(meshes[i - 1], "intersection");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh" + meshnum);
    MeshT mesh_i = makeSurfMesh<KernelT, MeshT, PointT>(
        rmesh_i,
        triangulate[i], // triangulate
        repairSoup,     // repair_soup
        false,          // remove_intersections
        1,              // remove_method
        false,          // fill_holes
        false,          // fair hole
        0);             // max_num_holes
    checkMesh1<MeshT>(mesh_i, i + 1);
    Message("Before corefine_and_compute_intersection().");
    const bool ok = PMP::corefine_and_compute_intersection(
      meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    } else {
      Message("Intersection computation succeeded.");
    }
  }
  Message("Before return.");
  return meshes[nMeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List boolIntersectionEK_cpp(const Rcpp::List rmeshes,
                                  const Rcpp::LogicalVector triangulate,
                                  const bool repairSoup,
                                  const bool normals) {
  EMesh3 mesh = boolIntersection<EK, EMesh3, EPoint3>(
      rmeshes,
      triangulate,
      repairSoup);
  // PMP::corefine_and_compute_intersection() requires triangle mesh
  // -> output is triangle
  return RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(mesh, normals, 3);
}

// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
MeshT boolDifferenceSM(MeshT &mesh1,
                       MeshT &mesh2,
                       const bool triangulate1,
                       const bool triangulate2,
                       const bool repairSoup) {
  checkMesh1<MeshT>(mesh1, 1);
  checkMesh1<MeshT>(mesh2, 2);
  MeshT mesh_d;
  bool ok = PMP::corefine_and_compute_difference(mesh1, mesh2, mesh_d);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return mesh_d;
}

// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
MeshT boolDifference(const Rcpp::List &rmesh1,
                     const Rcpp::List &rmesh2,
                     const bool triangulate1,
                     const bool triangulate2,
                     const bool repairSoup) {
  Message("Processing mesh1");
  MeshT smesh1 = makeSurfMesh<KernelT, MeshT, PointT>(
      rmesh1,
      triangulate1,   // triangulate
      repairSoup,     // repair_soup
      false,          // remove_intersections
      1,              // remove_method
      false,          // fill_holes
      false,          // fair hole
      0);             // max_num_holes
  checkMesh1<MeshT>(smesh1, 1);
  Message("Processing mesh2");
  MeshT smesh2 = makeSurfMesh<KernelT, MeshT, PointT>(
      rmesh2,
      triangulate2,   // triangulate
      repairSoup,     // repair_soup
      false,          // remove_intersections
      1,              // remove_method
      false,          // fill_holes
      false,          // fair hole
      0);             // max_num_holes
  checkMesh1<MeshT>(smesh2, 2);
  MeshT mesh_d;
  bool ok = PMP::corefine_and_compute_difference(smesh1, smesh2, mesh_d);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return mesh_d;
}

// [[Rcpp::export]]
Rcpp::List boolDifferenceEK_cpp(const Rcpp::List rmesh1,
                                const Rcpp::List rmesh2,
                                const bool triangulate1,
                                const bool triangulate2,
                                const bool repairSoup,
                                const bool normals) {
  EMesh3 mesh = boolDifference<EK, EMesh3, EPoint3>(
    rmesh1, rmesh2, triangulate1, triangulate2, repairSoup);
  // PMP::corefine_and_compute_difference() requires triangle mesh
  // -> output is triangle
  return RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(mesh, normals, 3);
}

// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
MeshT boolUnion(const Rcpp::List &rmeshes,
                const Rcpp::LogicalVector triangulate,
                const bool repairSoup) {
  const std::size_t nMeshes = rmeshes.size();
  std::vector<MeshT> meshes(nMeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh1");
  MeshT mesh_0 = makeSurfMesh<KernelT, MeshT, PointT>(
      rmesh,
      triangulate[0], // triangulate
      repairSoup,     // repair_soup
      false,          // remove_intersections
      1,              // remove_method
      false,          // fill_holes
      false,          // fair hole
      0);             // max_num_holes
  meshes[0] = mesh_0;
  for(std::size_t i = 1; i < nMeshes; i++) {
    if(i == 1) {
      checkMesh1<MeshT>(meshes[0], 1);
    } else {
      checkMesh2<MeshT>(meshes[i - 1], "union");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh" + meshnum);
    MeshT mesh_i = makeSurfMesh<KernelT, MeshT, PointT>(
        rmesh_i,
        triangulate[i], // triangulate
        repairSoup,     // repair_soup
        false,          // remove_intersections
        1,              // remove_method
        false,          // fill_holes
        false,          // fair hole
        0);             // max_num_holes
    checkMesh1<MeshT>(mesh_i, i + 1);
    const bool ok =
        PMP::corefine_and_compute_union(meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Union computation has failed.");
    }
  }

  return meshes[nMeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List boolUnionEK_cpp(const Rcpp::List rmeshes,
                           const Rcpp::LogicalVector triangulate,
                           const bool repairSoup,
                           const bool normals) {
  EMesh3 mesh = boolUnion<EK, EMesh3, EPoint3>(
      rmeshes,
      triangulate,
      repairSoup);
  // PMP::corefine_and_compute_union() requires triangle mesh
  // -> output is triangle
  return RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(mesh, normals, 3);
}
