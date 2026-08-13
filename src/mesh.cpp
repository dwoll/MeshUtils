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
Rcpp::List SurfEMesh_cpp(const Rcpp::List rmesh,
                         const bool isTriangle,
                         const bool triangulate,
                         const bool clean,
                         const bool normals) {
  Message("\u2014 Processing mesh...");
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, clean, false);
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::DataFrame Edges0;
  Rcpp::NumericMatrix Normals0;
  if(really_triangulate) {
    Edges0 = getEdges<EK, EMesh3, EPoint3>(mesh);
    if(normals) {
      Normals0 = getEKNormals(mesh);
    }
    Message("Triangulation.");
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
    if(CGAL::is_closed(mesh)) {
      if(!PMP::is_outward_oriented(mesh)) {
        PMP::reverse_face_orientations(mesh);
      }
      const bool bv = PMP::does_bound_a_volume(mesh);
      std::string msg2;
      if(bv) {
        msg2 = "The mesh bounds a volume.";
      } else {
        msg2 = "The mesh does not bound a volume - reorienting.";
        PMP::orient_to_bound_a_volume(mesh);
      }
      Message(msg2);
    }
  }
  Message("... done.\n");
  Rcpp::List routmesh = RSurfEKMesh(mesh, normals);
  if(really_triangulate) {
    routmesh["edges0"] = Edges0;
    if(normals) {
      routmesh["normals0"] = Normals0;
    }
  }
  return routmesh;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesBoundVolume_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  return PMP::does_bound_a_volume(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesSelfIntersect_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  return PMP::does_self_intersect(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool isClosed_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
  return CGAL::is_closed(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List orientToBoundVolume_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
  if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  PMP::orient_to_bound_a_volume(mesh);
  // TODO necessary to update normals?
  return getRmesh(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List removeSelfIntersections_cpp(const Rcpp::List rmesh_in, const unsigned int method) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh_in, false, false);
  Rcpp::List rmesh_out;

  if(!CGAL::is_triangle_mesh(mesh)) {
    Message("The mesh is not triangle. Nothing done.\n");
    rmesh_out = getRmesh(mesh);
  }
  if(PMP::does_self_intersect(mesh)) {
      CGAL::Conforming_constrained_Delaunay_triangulation_3<EK> ccdt;
      if(method == 1) {
          // use a polygon soup as container as the output will most likely be non-manifold
          std::vector<EPoint3> points;
          std::vector<std::vector<std::size_t>> polygons;
          PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
          bool success = PMP::autorefine_triangle_soup(points, polygons);
          // PMP::does_polygon_soup_self_intersect(points, polygons));
          if(success) {
              Message("Autorefine successful.\n");
          } else {
              Message("Autorefine not successful.\n");
          }

          ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
          bool clean = true;
          bool triangulate = false;
          EMesh3 meshrsi = soup_to_mesh<EMesh3, EPoint3>(points, polygons, clean, triangulate);
          if(PMP::does_self_intersect(meshrsi)) {
            Message("Self intersections could not be removed.\n");
          } else {
            Message("Self intersections removed.\n");
          }

          rmesh_out = getRmesh(meshrsi);
      } else if(method == 2) {
          // use a polygon soup as container as the output will most likely be non-manifold
          std::vector<EPoint3> points;
          std::vector<std::vector<std::size_t>> polygons;
          const auto& snap = CGAL::parameters::apply_iterative_snap_rounding(true);
          PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
          bool success = PMP::autorefine_triangle_soup(points, polygons, snap);
          if(success) {
              Message("Autorefine (with iterative snap) successful.\n");
          } else {
              Message("Autorefine (with iterative snap) not successful.\n");
          }

          ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
          bool clean = true;
          bool triangulate = false;
          EMesh3 meshrsi = soup_to_mesh<EMesh3, EPoint3>(points, polygons, clean, triangulate);
          if(PMP::does_self_intersect(meshrsi)) {
            Message("Self intersections could not be removed (with iterative snap).\n");
          } else {
            Message("Self intersections removed (with iterative snap).\n");
          }

          rmesh_out = getRmesh(meshrsi);
      }
  } else {
    Message("Mesh does not self-intersect. Nothing done.\n");
    rmesh_out = getRmesh(mesh);
  }

  return rmesh_out;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getVolume_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
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
  const EK::FT vol = PMP::volume(mesh);
  return CGAL::to_double<EK::FT>(vol);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::NumericVector getCentroid_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
  Rcpp::NumericVector out(3);
  if(!CGAL::is_triangle_mesh(mesh)) {
      Message("The mesh is not triangle.");
      out(0) = Rcpp::NumericVector::get_na();
      out(1) = Rcpp::NumericVector::get_na();
      out(2) = Rcpp::NumericVector::get_na();
  } else {
      Mesh3 epickCopy;
      CGAL::copy_face_graph(mesh, epickCopy);
      const Point3 centroid = PMP::centroid(epickCopy);
      out(0) = centroid.x();
      out(1) = centroid.y();
      out(2) = centroid.z();
      // const EPoint3 centroid = PMP::centroid(mesh);
      // out(0) = CGAL::to_double<EK::FT>(centroid.x());
      // out(1) = CGAL::to_double<EK::FT>(centroid.y());
      // out(2) = CGAL::to_double<EK::FT>(centroid.z());
  }
  return out;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List optimalBoundingBox_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
  Mesh3 kmesh;
  CGAL::copy_face_graph(mesh, kmesh);
  std::array<Point3, 8> obb_points;
  CGAL::oriented_bounding_box(kmesh, obb_points,
                              CGAL::parameters::use_convex_hull(true));
  // Make a mesh out of the oriented bounding box
  Mesh3 obbMesh;
  CGAL::make_hexahedron(
    obb_points[0], obb_points[1], obb_points[2], obb_points[3],
    obb_points[4], obb_points[5], obb_points[6], obb_points[7],
    obbMesh
  );
  EMesh3 obbEMesh;
  CGAL::copy_face_graph(obbMesh, obbEMesh);
  Rcpp::List rmesh_out = RSurfEKMesh2(obbEMesh, false, 4);
  Rcpp::NumericMatrix hxVertices(3, 8);
  for(int i = 0; i < 8; i++) {
    Point3 pt = obb_points[i];
    Rcpp::NumericVector v =
      Rcpp::NumericVector::create(pt.x(), pt.y(), pt.z());
    hxVertices(Rcpp::_, i) = v;
  }
  return Rcpp::List::create(
    Rcpp::Named("mesh") = rmesh_out,
    Rcpp::Named("hxVertices") = hxVertices
  );
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List boundingBox_cpp(const Rcpp::List rmesh) {
  EMesh3 mesh = makeSurfMesh<EMesh3, EPoint3>(rmesh, false, false);
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  Rcpp::NumericVector lcorner = { bbox.xmin(), bbox.ymin(), bbox.zmin() };
  Rcpp::NumericVector ucorner = { bbox.xmax(), bbox.ymax(), bbox.zmax() };
  return Rcpp::List::create(
    Rcpp::Named("lcorner") = lcorner,
    Rcpp::Named("ucorner") = ucorner
  );
}
