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
// initial mesh generation -> use EPEC kernel to enable
// filling holes and removing self-intersections
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
  EMesh3 mesh = makeSurfMesh<EK, EMesh3, EPoint3>(
      rmesh,
      triangulate,         // triangulate
      repairSoup,          // repair_soup
      removeIntersections, // remove_intersections
      removeMethod,        // remove_method
      fillHoles,           // fill_holes
      fairHole,            // fair hole
      maxNumHoles);        // max_num_holes
  return RSurfMesh1<EK, EMesh3, EPoint3, EVector3>(mesh, normals);
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
  return getRmesh(mesh, false);
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
  return getRmesh(mesh, false);
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
  return getRmesh(mesh, false);
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
  Rcpp::NumericVector out(3);
  if(!CGAL::is_triangle_mesh(mesh)) {
      Message("The mesh is not triangle.");
      out(0) = Rcpp::NumericVector::get_na();
      out(1) = Rcpp::NumericVector::get_na();
      out(2) = Rcpp::NumericVector::get_na();
  } else {
      const Point3 centroid = PMP::centroid(mesh);
      out(0) = CGAL::to_double<K::FT>(centroid.x());
      out(1) = CGAL::to_double<K::FT>(centroid.y());
      out(2) = CGAL::to_double<K::FT>(centroid.z());
  }
  return out;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List optimalBoundingBox_cpp(const Rcpp::List rmeshIn) {
  Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
      rmeshIn,
      false,       // triangulate
      false,       // repair_soup
      false,       // remove_intersections
      1,           // remove_method
      false,       // fill_holes
      false,       // fair hole
      0);          // max_num_holes
  std::array<Point3, 8> obbPoints;
  CGAL::oriented_bounding_box(mesh, obbPoints,
                              CGAL::parameters::use_convex_hull(true));
  // Make a mesh out of the oriented bounding box
  Mesh3 obbMesh;
  CGAL::make_hexahedron(
    obbPoints[0], obbPoints[1], obbPoints[2], obbPoints[3],
    obbPoints[4], obbPoints[5], obbPoints[6], obbPoints[7],
    obbMesh
  );
  Rcpp::List rmeshOut = RSurfMesh2<K, Mesh3, Point3, Vector3>(obbMesh, false, 4);
  Rcpp::NumericMatrix hxVertices(3, 8);
  for(int i = 0; i < 8; i++) {
    Point3 pt = obbPoints[i];
    Rcpp::NumericVector v =
      Rcpp::NumericVector::create(pt.x(), pt.y(), pt.z());
    hxVertices(Rcpp::_, i) = v;
  }
  return Rcpp::List::create(
    Rcpp::Named("mesh") = rmeshOut,
    Rcpp::Named("hxVertices") = hxVertices
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
    for(std::size_t i = 0; i < nPts; i++){
      distances(i) = Rcpp::NumericVector::get_na();
    }
  } else {
      for(std::size_t i = 0; i < nPts; i++){
        Rcpp::NumericVector point_i = points(Rcpp::_, i);
        std::vector<Point3> pt = { Point3(point_i(0), point_i(1), point_i(2)) };
        distances(i) = PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(pt, mesh);
      }
  }
  return distances;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getHausdorffApprox_cpp(
    const Rcpp::List rmesh1,
    const Rcpp::List rmesh2,
    const bool symmetric,
    const bool triangulate1 = false,
    const bool triangulate2 = false) {
  Mesh3 mesh1 = makeSurfMesh<K, Mesh3, Point3>(
      rmesh1,
      triangulate1, // triangulate
      false,        // repair_soup
      false,        // remove_intersections
      1,            // remove_method
      false,        // fill_holes
      false,        // fair hole
      0);           // max_num_holes
  Mesh3 mesh2 = makeSurfMesh<K, Mesh3, Point3>(
      rmesh2,
      triangulate2, // triangulate
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
    const double errorBound,
    const bool triangulate1 = false,
    const bool triangulate2 = false) {
    Mesh3 mesh1 = makeSurfMesh<K, Mesh3, Point3>(
        rmesh1,
        triangulate1, // triangulate
        false,        // repair_soup
        false,        // remove_intersections
        1,            // remove_method
        false,        // fill_holes
        false,        // fair hole
        0);           // max_num_holes
    Mesh3 mesh2 = makeSurfMesh<K, Mesh3, Point3>(
        rmesh2,
        triangulate2, // triangulate
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
            mesh1, mesh2, errorBound));
    } else {
        d = CGAL::to_double<K::FT>(PMP::bounded_error_Hausdorff_distance<PIA_TAG>(
            mesh1, mesh2, errorBound));
    }
    return d;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshIsotropicUniform_cpp(
    const Rcpp::List rmesh,
    const double targetEdgeLen,
    const unsigned int nIter,
    const unsigned int nRelaxSteps) {
    Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
        rmesh,
        true,        // triangulate
        false,       // repair_soup
        false,       // remove_intersections
        1,           // remove_method
        false,       // fill_holes
        false,       // fair hole
        0);          // max_num_holes
    std::vector<hedgdescr> borderHalfEdges;
    PMP::border_halfedges(faces(mesh), mesh, std::back_inserter(borderHalfEdges));
    std::vector<edgdescr> border;
    std::size_t nheBorder = borderHalfEdges.size();
    border.reserve(nheBorder);
    for(std::size_t i = 0; i < nheBorder; i++) {
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
    return getRmesh(mesh, false);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshIsotropicAdapt_cpp(
    const Rcpp::List rmesh,
    const double tol,
    const double edgeMin,
    const double edgeMax,
    const unsigned int nIter,
    const unsigned int nRelaxSteps) {
    Mesh3 mesh = makeSurfMesh<K, Mesh3, Point3>(
        rmesh,
        true,        // triangulate
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
    return getRmesh(mesh, false);
}

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
    removeProperties(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::CatmullClark_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter));
    mesh.collect_garbage();
    return getRmesh(mesh, triangulate);
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
    removeProperties(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::DooSabin_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter));
    mesh.collect_garbage();
    return getRmesh(mesh, triangulate);
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
    removeProperties(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::Sqrt3_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter));
    return getRmesh(mesh, triangulate);
}
