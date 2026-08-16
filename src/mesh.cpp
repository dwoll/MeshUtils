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
Rcpp::List SurfMesh_cpp(const Rcpp::List rmeshIn,
                        const bool clean,
                        const bool triangulate,
                        const bool normals) {
  Message("Processing mesh...");
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmeshIn, clean, triangulate);
  Message("... done.\n");
  Rcpp::List rmeshOut = RSurfMesh1<K, Mesh3, Point3, Vector3>(mesh, normals);
  return rmeshOut;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesBoundVolume_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  return PMP::does_bound_a_volume(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesSelfIntersect_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  return PMP::does_self_intersect(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool isClosed_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  return CGAL::is_closed(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List orientToBoundVolume_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  PMP::orient_to_bound_a_volume(mesh);
  return getRmesh(mesh);
}

// ----------------------------------------------------------------------- //
// use EPEC kernel for autorefine_triangle_soup()
// [[Rcpp::export]]
Rcpp::List removeSelfIntersections_cpp(const Rcpp::List rmeshIn, const unsigned int method) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmeshIn, false, false);
  Rcpp::List rmeshOut;

  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.\n");
  }
  if(PMP::does_self_intersect(mesh)) {
      // use a polygon soup as container as the output will most likely be non-manifold
      std::vector<EPoint3> points;
      std::vector<std::vector<std::size_t>> polygons;
      PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
      bool success;
      if(method == 1) {
          success = PMP::autorefine_triangle_soup(points, polygons);
      } else if(method == 2) {
          const auto& snap = CGAL::parameters::apply_iterative_snap_rounding(true);
          success = PMP::autorefine_triangle_soup(points, polygons, snap);
      }
      if(success) {
          Message("Autorefine successful.\n");
      } else {
          Message("Autorefine not successful.\n");
      }

      // PMP::does_polygon_soup_self_intersect(points, polygons));
      CGAL::Conforming_constrained_Delaunay_triangulation_3<EK> ccdt;
      ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
      bool clean = true;
      bool triangulate = false;
      EMesh3 meshNSI = soup_to_mesh<EMesh3, EPoint3>(points, polygons, clean, triangulate);
      if(PMP::does_self_intersect(meshNSI)) {
        Message("Self intersections could not be removed.\n");
      } else {
        Message("Self intersections removed.\n");
      }

      rmeshOut = getRmesh(meshNSI);
  } else {
      Message("Mesh does not self-intersect. Nothing done.\n");
      rmeshOut = getRmesh(mesh);
  }

  return rmeshOut;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getArea_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Message("The mesh is not triangle.");
    return Rcpp::NumericVector::get_na();
  }
  if(PMP::does_self_intersect(mesh)) {
    Message("The mesh self-intersects.");
    return Rcpp::NumericVector::get_na();
  }
  const K::FT a = PMP::area(mesh);
  return CGAL::to_double<K::FT>(a);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getVolume_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  if(!CGAL::is_closed(mesh)) {
    Message("The mesh is not closed.");
    return Rcpp::NumericVector::get_na();
  }
  if(!CGAL::is_triangle_mesh(mesh)) {
    Message("The mesh is not triangle.");
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
Rcpp::NumericVector getCentroid_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  Rcpp::NumericVector out(3);
  if(!CGAL::is_triangle_mesh(mesh)) {
      Message("The mesh is not triangle.");
      out(0) = Rcpp::NumericVector::get_na();
      out(1) = Rcpp::NumericVector::get_na();
      out(2) = Rcpp::NumericVector::get_na();
  } else {
      const Point3 centroid = PMP::centroid(mesh);
      out(0) = centroid.x();
      out(1) = centroid.y();
      out(2) = centroid.z();
  }
  return out;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List optimalBoundingBox_cpp(const Rcpp::List rmeshIn) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmeshIn, false, false);
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
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
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
    const Rcpp::List rmesh, const Rcpp::NumericMatrix points) {
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, false, false);
  const size_t nPts = points.ncol();
  Rcpp::NumericVector distances(nPts);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Message("The mesh is not triangle.");
    for(size_t i = 0; i < nPts; i++){
      distances(i) = Rcpp::NumericVector::get_na();
    }
  } else {
      for(size_t i = 0; i < nPts; i++){
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
    const Rcpp::List rmesh1, const Rcpp::List rmesh2, bool symmetric) {
  Mesh3 mesh1 = makeSurfMesh<Mesh3, Point3>(rmesh1, false, false);
  Mesh3 mesh2 = makeSurfMesh<Mesh3, Point3>(rmesh2, false, false);
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
    d = PMP::approximate_symmetric_Hausdorff_distance<PIA_TAG>(mesh1, mesh2);
  } else {
    d = PMP::approximate_Hausdorff_distance<PIA_TAG>(mesh1, mesh2);
  }
  return d;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getHausdorffEst_cpp(
    const Rcpp::List rmesh1, const Rcpp::List rmesh2, bool symmetric, double errorBound) {
    Mesh3 mesh1 = makeSurfMesh<Mesh3, Point3>(rmesh1, false, false);
    Mesh3 mesh2 = makeSurfMesh<Mesh3, Point3>(rmesh2, false, false);
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
        d = PMP::bounded_error_symmetric_Hausdorff_distance<PIA_TAG>(
            mesh1, mesh2, errorBound);
  } else {
        d = PMP::bounded_error_Hausdorff_distance<PIA_TAG>(
            mesh1, mesh2, errorBound);
  }
  return d;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshIsotropic_cpp(
    const Rcpp::List rmeshIn,
    const double targetEdgeLen,
    const unsigned nIter,
    const unsigned nRelaxSteps) {
    Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmeshIn, false, false);
    std::vector<hedgdescr> borderHalfEdges;
    PMP::border_halfedges(mesh.faces(), mesh, std::back_inserter(borderHalfEdges));
    std::vector<edgdescr> border;
    unsigned int nhBorder = borderHalfEdges.size();
    border.reserve(nhBorder);
    for(unsigned int i = 0; i < nhBorder; i++) {
      border.emplace_back(mesh.edge(borderHalfEdges[i]));
    }
    PMP::split_long_edges(border, targetEdgeLen, mesh);
    PMP::isotropic_remeshing(
      mesh.faces(),
      targetEdgeLen,
      mesh,
      PMP::parameters::number_of_iterations(nIter)
                      .number_of_relaxation_steps(nRelaxSteps)
                      .protect_constraints(true));
    mesh.collect_garbage();
    Rcpp::List rmeshOut = getRmesh(mesh);
    return rmeshOut;
}

/*
// TODO need to fix removeProperties() first
// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// use EPEC kernel here
Rcpp::List subdivideCatmullClark_cpp(const Rcpp::List rmeshIn, unsigned int nIter) {
    EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmeshIn, false, false);
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::CatmullClark_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter)
    );
    mesh.collect_garbage();
    rmeshOut = getRmesh(mesh);
    return rmeshOut;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// use EPEC kernel here
Rcpp::List subdivideDooSabin_cpp(const Rcpp::List rmeshIn, unsigned int nIter) {
    EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmeshIn, false, false);
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::DooSabin_subdivision(
      mesh, CGAL::parameters::number_of_iterations(nIter));
    mesh.collect_garbage();
    rmeshOut = getRmesh(mesh);
    return rmeshOut;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// use EPEC kernel here
Rcpp::List subdivideSqrt3_cpp(const Rcpp::List rmeshIn, unsigned int nIter) {
    EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmeshIn, false, false);
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties(mesh, {"v:normal"});
    CGAL::Subdivision_method_3::Sqrt3_subdivision(
      mesh, CGAL::parameters::number_of_iterations(iterations));
      rmeshOut = getRmesh(mesh);
      return rmeshOut;
}
*/
