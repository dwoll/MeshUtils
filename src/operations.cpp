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
void checkMesh(MeshT mesh, size_t i) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "Mesh n\u00b0" + std::to_string(i) + " self-intersects.";
    Rcpp::stop(msg);
  }
}

// ----------------------------------------------------------------------- //
template <typename MeshT>
void checkMesh2(MeshT mesh, const std::string& what) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "The " + what + " self-intersects.";
    Rcpp::stop(msg);
  }
}

// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
MeshT boolIntersect(const Rcpp::List rmeshes,
                    const bool clean,
                    const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh1");
  MeshT mesh_0 = makeSurfMesh<MeshT, PointT>(rmesh, clean, triangulate[0]);
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    if(i == 1) {
      checkMesh<MeshT>(meshes[0], 1);
    } else {
      checkMesh2<MeshT>(meshes[i - 1], "intersection");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh" + meshnum);
    MeshT mesh_i = makeSurfMesh<MeshT, PointT>(rmesh_i, clean, triangulate[i]);
    checkMesh<MeshT>(mesh_i, i + 1);
    const bool ok = PMP::corefine_and_compute_intersection(
      meshes[i - 1], mesh_i, meshes[i]
    );
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List intersectionEK_cpp(const Rcpp::List rmeshes,
                              const bool clean,
                              const bool normals,
                              const Rcpp::LogicalVector triangulate) {
  EMesh3 mesh = boolIntersect<EK, EMesh3, EPoint3>(rmeshes, clean, triangulate);
  return RSurfTEKMesh(mesh, normals);
}

// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
MeshT boolDiff(const Rcpp::List rmesh1,
               const Rcpp::List rmesh2,
               const bool clean,
               const bool triangulate1,
               const bool triangulate2) {
  Message("Processing mesh1");
  MeshT smesh1 = makeSurfMesh<MeshT, PointT>(rmesh1, clean, triangulate1);
  checkMesh<MeshT>(smesh1, 1);
  Message("Processing mesh2");
  MeshT smesh2 = makeSurfMesh<MeshT, PointT>(rmesh2, clean, triangulate2);
  checkMesh<MeshT>(smesh2, 2);
  MeshT outmesh;
  bool ok = PMP::corefine_and_compute_difference(smesh1, smesh2, outmesh);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return outmesh;
}

// [[Rcpp::export]]
Rcpp::List differenceEK_cpp(const Rcpp::List rmesh1,
                            const Rcpp::List rmesh2,
                            const bool clean,
                            const bool normals,
                            const bool triangulate1,
                            const bool triangulate2) {
  EMesh3 mesh = boolDiff<EK, EMesh3, EPoint3>(
    rmesh1, rmesh2, clean, triangulate1, triangulate2);
  return RSurfTEKMesh(mesh, normals);
}

// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
MeshT boolUnion(const Rcpp::List rmeshes,
                const bool clean,
                const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("Processing mesh1");
  MeshT mesh_0 = makeSurfMesh<MeshT, PointT>(rmesh, clean, triangulate[0]);
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    if(i == 1) {
      checkMesh<MeshT>(meshes[0], 1);
    } else {
      checkMesh2<MeshT>(meshes[i - 1], "union");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("Processing mesh" + meshnum);
    MeshT mesh_i = makeSurfMesh<MeshT, PointT>(rmesh_i, clean, triangulate[i]);
    checkMesh<MeshT>(mesh_i, i + 1);
    const bool ok =
        PMP::corefine_and_compute_union(meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Union computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List unionEK_cpp(const Rcpp::List rmeshes,
                       const bool clean,
                       const bool normals,
                       const Rcpp::LogicalVector triangulate) {
  EMesh3 mesh = boolUnion<EK, EMesh3, EPoint3>(rmeshes, clean, triangulate);
  return RSurfTEKMesh(mesh, normals);
}
