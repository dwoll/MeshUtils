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
// ----------------------------------------------------------------------- //
void Message(std::string msg) {
  SEXP rmsg = Rcpp::wrap(msg);
  Rcpp::message(rmsg);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix &M) {
  const size_t nPts = M.ncol();
  std::vector<PointT> points;
  points.reserve(nPts);
  for(std::size_t i = 0; i < nPts; i++) {
    const Rcpp::NumericVector pt = M(Rcpp::_, i);
    points.emplace_back(PointT(pt(0), pt(1), pt(2)));
  }
  return points;
}

template std::vector<Point3>  matrix_to_points3<Point3>(const Rcpp::NumericMatrix&);
template std::vector<EPoint3> matrix_to_points3<EPoint3>(const Rcpp::NumericMatrix&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// create triangle faces
std::vector<std::vector<std::size_t>> matrix_to_Tfaces(
  const Rcpp::IntegerMatrix &face_mat) {
  const std::size_t nFaces = face_mat.ncol();
  std::vector<std::vector<std::size_t>> faces;
  faces.reserve(nFaces);
  for(std::size_t i = 0; i < nFaces; i++) {
    const Rcpp::IntegerVector face_rcpp = face_mat(Rcpp::_, i);
    // need static_cast here because initializing with {} instead of ()
    std::vector<std::size_t> face = { static_cast<std::size_t>(face_rcpp(0)),
                                      static_cast<std::size_t>(face_rcpp(1)),
                                      static_cast<std::size_t>(face_rcpp(2)) };
    faces.emplace_back(face);
  }
  return faces;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// create faces - may not be triangle
std::vector<std::vector<std::size_t>> list_to_faces1(const Rcpp::List &L) {
  const std::size_t nFaces = L.size();
  std::vector<std::vector<std::size_t>> faces;
  faces.reserve(nFaces);
  for(std::size_t i = 0; i < nFaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    std::vector<std::size_t> face(face_rcpp.begin(), face_rcpp.end());
    faces.emplace_back(face);
  }
  return faces;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
std::pair<std::vector<std::vector<std::size_t>>, bool> list_to_faces2(
    const Rcpp::List &L) {
  const std::size_t nFaces = L.size();
  std::vector<std::vector<std::size_t>> faces;
  faces.reserve(nFaces);
  bool triangle = true;
  for(std::size_t i = 0; i < nFaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    std::vector<std::size_t> face(face_rcpp.begin(), face_rcpp.end());
    // std::transform(
    //     face.begin(), face.end(), face.begin(),
    // 	   std::bind(std::minus<int>(), std::placeholders::_1, 1));
    faces.emplace_back(face);
    triangle = triangle && (face.size() == 3);
  }
  return std::make_pair(faces, triangle);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
bool is_triangle_soup(const std::vector<std::vector<std::size_t>>& polygons) {
    for (const auto& poly : polygons) {
        if (poly.size() != 3) {
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::polygon_soup_to_polygon_mesh() with a lot of checks
// hole filling, removing self-intersections
// points and faces are changed -> no const, no reference
template <typename KernelT, typename MeshT, typename PointT>
MeshT soup_to_mesh(std::vector<PointT> points,
                   std::vector<std::vector<std::size_t>> faces,
                   const bool triangulate,
                   const bool repair_soup,
                   const bool remove_intersections,
                   const int remove_method,
                   const bool fill_holes,
                   const bool fair_hole,
                   const unsigned int max_num_holes) {
    // determine EPEC vs. EPIC kernel
    /*
    bool is_epeck;
    if constexpr (std::is_same_v<MeshT, EMesh3>) {
        is_epeck = true;
        Message("EP_E_CK");
    } else {
        is_epeck = false;
        Message("EP_I_CK");
    }
    */
    if(repair_soup) {
        PMP::repair_polygon_soup(points, faces);
    }
    const bool is_oriented = PMP::orient_polygon_soup(points, faces);
    if(!is_oriented) {
        Message("Polygon orientation failed. Remove holes / self-intersections if present.");
    }
    // triangulate if necessary
    bool is_triangle = is_triangle_soup(faces);
    if(triangulate && !is_triangle) {
        is_triangle = PMP::triangulate_polygons(points, faces);
    }
    if(is_triangle) {
        Message("Mesh is triangle");
    } else {
        Message("Mesh is not triangle.\n  Cannot ensure it bounds a volume.\n  Cannot try to remove self-intersections.");
    }
    // check for self-intersections
    const bool soup_has_self_int = PMP::does_polygon_soup_self_intersect(points, faces);
    if(soup_has_self_int) {
        Message("Polygon soup has self-intersections.");
    }
    // remove self-intersections if necessary and possible
    MeshT mesh;
    if(soup_has_self_int && remove_intersections && is_triangle) {
        removeSelfIntSoup<KernelT, PointT>(points, faces, remove_method);
        PMP::orient_polygon_soup(points, faces);
        PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
    } else {
        PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
    }

    // filling boundary holes if necessary and possible
    std::string msg_closed;
    std::string msg_notclosed;
    if(!CGAL::is_closed(mesh) && fill_holes && (max_num_holes > 0)) {
        // mesh is passed by reference and modified in fillBoundaryHoles()
        // TODO also pass parameters related to hole size
        MeshT mesh_tmp = fillBoundaryHoles<MeshT, PointT>(mesh, fair_hole, -1, -1, max_num_holes);
        mesh = std::move(mesh_tmp);
        msg_closed    = "now ";
        msg_notclosed = "still ";
    }
    if(CGAL::is_closed(mesh)) {
        msg_closed    = "Mesh is " + msg_closed    + "closed.";
        Message(msg_closed);
    } else {
        msg_notclosed = "Mesh is " + msg_notclosed + "not closed.";
        Message(msg_notclosed);
    }
    if(PMP::does_self_intersect(mesh)) {
        Message("Mesh has self-intersections.");
    } else {
        Message("Mesh does not have self-intersections.");
    }
    if(mesh.is_valid(true)) {
        Message("Mesh is valid.\n");
    } else {
        Message("Mesh is not valid.\n");
    }
    return mesh;
}

template Mesh3 soup_to_mesh<K, Mesh3, Point3>(
    std::vector<Point3>,
    std::vector<std::vector<std::size_t>>,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

template EMesh3 soup_to_mesh<EK, EMesh3, EPoint3>(
    std::vector<EPoint3>,
    std::vector<std::vector<std::size_t>>,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// general conversion from R list to Surface_mesh_3
// turn R data structures for vertices and faces into C++ vectors
// then call soup_to_mesh()
// faces may not be triangle -> list
template <typename KernelT, typename MeshT, typename PointT>
MeshT makeSurfMesh(
  const Rcpp::List &rmesh,
  const bool triangulate,
  const bool repair_soup,
  const bool remove_intersections,
  const int remove_method,
  const bool fill_holes,
  const bool fair_hole,
  const unsigned int max_num_holes) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<std::size_t>> faces = list_to_faces1(rfaces);
  return soup_to_mesh<KernelT, MeshT, PointT>(
      points,
      faces,
      triangulate,
      repair_soup,
      remove_intersections,
      remove_method,
      fill_holes,
      fair_hole,
      max_num_holes);
}

template Mesh3 makeSurfMesh<K, Mesh3, Point3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

template EMesh3 makeSurfMesh<EK, EMesh3, EPoint3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// like makeSurfMesh() but for triangles -> rfaces is matrix
// currently unused
template <typename KernelT, typename MeshT, typename PointT>
MeshT makeSurfTMesh(
    const Rcpp::List &rmesh,
    const bool repair_soup,
    const bool remove_intersections,
    const int remove_method,
    const bool fill_holes,
    const bool fair_hole,
    const unsigned int max_num_holes) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::IntegerMatrix rfaces =
      Rcpp::as<Rcpp::IntegerMatrix>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<std::size_t>> faces = matrix_to_Tfaces(rfaces);
  return soup_to_mesh<KernelT, MeshT, PointT>(
      points,
      faces,
      false,               // triangulate
      repair_soup,
      remove_intersections,
      remove_method,
      fill_holes,
      fair_hole,
      max_num_holes);
}

template Mesh3 makeSurfTMesh<K, Mesh3, Point3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

template EMesh3 makeSurfTMesh<EK, EMesh3, EPoint3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);
