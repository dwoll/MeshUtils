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

#include <CGAL/optimal_bounding_box.h>

// ----------------------------------------------------------------------- //
// initial mesh generation - EPIC kernel - TODO make parameter
// [[Rcpp::export]]
Rcpp::List makeMesh_cpp(const Rcpp::List rmesh,
                        const bool triangulate,
                        const bool repairSoup,
                        const bool removeIntersections,
                        const int removeMethod,
                        const bool fillHoles,
                        const bool fairHole,
                        const unsigned int maxNumHoles,
                        const bool normals) {
  Message("Processing mesh...");
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      triangulate,         // triangulate
      repairSoup,          // repair_soup
      removeIntersections, // remove_intersections
      removeMethod,        // remove_method
      fillHoles,           // fill_holes
      fairHole,            // fair hole
      maxNumHoles);        // max_num_holes
  return RSurfMesh1<K, Mesh3, Point3, Vector3>(mesh, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool isValid_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      false,       // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  return mesh.is_valid(false);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool hasGarbage_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      false,       // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  return mesh.has_garbage();
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesBoundVolume_cpp(
  const Rcpp::List rmesh, const bool triangulate = false) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      triangulate, // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::warning("The mesh is not triangle.");
    return Rcpp::LogicalVector::get_na();
  }
  return PMP::does_bound_a_volume(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesSelfIntersect_cpp(
  const Rcpp::List rmesh, const bool triangulate = false) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      triangulate, // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::warning("The mesh is not triangle.");
    return Rcpp::LogicalVector::get_na();
  }
  return PMP::does_self_intersect(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool isClosed_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      false,       // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  return CGAL::is_closed(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List orientToBoundVolume_cpp(
  const Rcpp::List rmesh, const bool triangulate = false) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
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
  PMP::orient_to_bound_a_volume(mesh);
  return getRmesh<K, Mesh3, Point3, Vector3>(mesh, false);
}

// ----------------------------------------------------------------------- //
// use EPEC kernel for autorefine_triangle_soup()
// [[Rcpp::export]]
Rcpp::List removeSelfIntersections_cpp(
  const Rcpp::List rmesh,
  const bool triangulate,
  const int method) {
  EMesh3 mesh = makeSurfMesh<EK, EMesh3, EPoint3>(
      rmesh,
      triangulate, // triangulate
      true,        // repair_soup
      true,        // remove_intersections
      method,      // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  return getRmesh<EK, EMesh3, EPoint3, EVector3>(mesh, false);
}

// ----------------------------------------------------------------------- //
// use EPEC kernel for fillBoundaryHoles()
// [[Rcpp::export]]
Rcpp::List fillBoundaryHoles_cpp(
  const Rcpp::List rmesh,
  const bool fairHole,
  const unsigned int maxNumHoles) {
  EMesh3 mesh = makeSurfMesh<EK, EMesh3, EPoint3>(
      rmesh,
      true,         // triangulate
      true,         // repair_soup
      false,        // remove_intersections
      1,            // remove_method
      true,         // fill_holes
      fairHole,     // fair hole
      maxNumHoles); // max_num_holes
  return getRmesh<EK, EMesh3, EPoint3, EVector3>(mesh, false);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getArea_cpp(const Rcpp::List rmesh, const bool triangulate = false) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      triangulate, // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::warning("The mesh is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  if(PMP::does_self_intersect(mesh)) {
    Rcpp::warning("The mesh self-intersects.");
    return Rcpp::NumericVector::get_na();
  }
  const K::FT a = PMP::area(mesh);
  return CGAL::to_double<K::FT>(a);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getVolume_cpp(
  const Rcpp::List rmesh, const bool triangulate = false) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      triangulate, // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  if(!CGAL::is_triangle_mesh(mesh)) {
    Message("The mesh is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_closed(mesh)) {
    Message("The mesh is not closed.");
    return Rcpp::NumericVector::get_na();
  }
  if(PMP::does_self_intersect(mesh)) {
    Message("The mesh self-intersects.");
    return Rcpp::NumericVector::get_na();
  }
  const K::FT vol = PMP::volume(mesh);
  return CGAL::to_double<K::FT>(vol);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::NumericVector getCentroid_cpp(
  const Rcpp::List rmesh, const bool triangulate = false) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      triangulate, // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  Rcpp::NumericVector ctr(3);
  if(!CGAL::is_triangle_mesh(mesh)) {
      Message("The mesh is not triangle.");
      ctr(0) = Rcpp::NumericVector::get_na();
      ctr(1) = Rcpp::NumericVector::get_na();
      ctr(2) = Rcpp::NumericVector::get_na();
  } else {
      const Point3 centroid = PMP::centroid(mesh);
      ctr(0) = CGAL::to_double<K::FT>(centroid.x());
      ctr(1) = CGAL::to_double<K::FT>(centroid.y());
      ctr(2) = CGAL::to_double<K::FT>(centroid.z());
  }
  return ctr;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List optimalBoundingBox_cpp(const Rcpp::List rmeshIn, const bool normals) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmeshIn,
      false,       // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  std::array<Point3, 8> obb_pts;
  CGAL::oriented_bounding_box(mesh, obb_pts,
                              CGAL::parameters::use_convex_hull(true));
  // make mesh out of oriented bounding box
  Mesh3 obb_mesh;
  CGAL::make_hexahedron(
    obb_pts[0], obb_pts[1], obb_pts[2], obb_pts[3],
    obb_pts[4], obb_pts[5], obb_pts[6], obb_pts[7],
    obb_mesh
  );
  Rcpp::List rmesh_obb = RSurfMesh2<K, Mesh3, Point3, Vector3>(obb_mesh, normals, 4);
  Rcpp::NumericMatrix hex_verts(3, 8);
  for(int i = 0; i < 8; i++) {
    Point3 pt = obb_pts[i];
    Rcpp::NumericVector v =
      Rcpp::NumericVector::create(pt.x(), pt.y(), pt.z());
    hex_verts(Rcpp::_, i) = v;
  }
  return Rcpp::List::create(
    Rcpp::Named("mesh") = rmesh_obb,
    Rcpp::Named("hxVertices") = hex_verts
  );
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List boundingBox_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      false,       // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  Rcpp::NumericVector lcorner = { bbox.xmin(), bbox.ymin(), bbox.zmin() };
  Rcpp::NumericVector ucorner = { bbox.xmax(), bbox.ymax(), bbox.zmax() };
  return Rcpp::List::create(
    Rcpp::Named("lcorner") = lcorner,
    Rcpp::Named("ucorner") = ucorner
  );
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::NumericVector getDistance_cpp(
    const Rcpp::List rmesh, const Rcpp::NumericMatrix points, const bool triangulate = false) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmesh,
      triangulate, // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  const std::size_t nPts = points.ncol();
  Rcpp::NumericVector distances(nPts);
  if(!CGAL::is_triangle_mesh(mesh)) {
      Message("The mesh is not triangle.");
      for(std::size_t i = 0; i < nPts; i++) {
          distances(i) = Rcpp::NumericVector::get_na();
      }
  } else {
      for(std::size_t i = 0; i < nPts; i++) {
          Rcpp::NumericVector point_i = points(Rcpp::_, i);
          std::vector<Point3> pt = { Point3(point_i(0), point_i(1), point_i(2)) };
          distances(i) = PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(pt, mesh);
      }
  }
  return distances;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List addVnormals_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
    rmesh,
    false,       // triangulate
    false,       // repair_soup
    false,       // remove_intersections
    1,           // remove_method
    false,       // fill_holes
    false,       // fair hole
    0);          // max_num_holes
  Rcpp::List mesh_out;
  if(CGAL::is_triangle_mesh(mesh)) {
    mesh_out = RSurfMesh2<K, Mesh3, Point3, Vector3>(mesh, true, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    mesh_out = RSurfMesh2<K, Mesh3, Point3, Vector3>(mesh, true, 4);
  } else {
    mesh_out = RSurfMesh1<K, Mesh3, Point3, Vector3>(mesh, true);
  }
  return mesh_out;
}
