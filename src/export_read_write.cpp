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

#include <CGAL/IO/io.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <locale>     // tolower()
#include <filesystem> // path()

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
std::string toLower(std::string s) {
  for(char& c : s) {
    c = std::tolower(c);
  }
  return s;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List readFileSoup_cpp(const std::string filename, const bool verbose) {
    std::vector<Point3> points;
    std::vector<std::vector<std::size_t>> faces;
    const bool ok = CGAL::IO::read_polygon_soup(
        filename, points, faces, CGAL::parameters::verbose(verbose));
    if(!ok) {
      Rcpp::stop("Reading failure.");
    }
    const std::size_t nPts = points.size();
    Rcpp::NumericMatrix vertex_mat(3, nPts);
    for(std::size_t i = 0; i < nPts; i++) {
      const Point3 point_i = points[i];
      Rcpp::NumericVector col_i =
          Rcpp::NumericVector::create(point_i.x(), point_i.y(), point_i.z());
      vertex_mat(Rcpp::_, i) = col_i;
    }
    const std::size_t nFaces = faces.size();
    Rcpp::List face_list(nFaces);
    for(std::size_t i = 0; i < nFaces; i++) {
      const std::vector<std::size_t> face_i = faces[i];
      Rcpp::IntegerVector col_i(face_i.begin(), face_i.end());
      face_list(i) = col_i + 1;
    }
    Rcpp::List out;
    out["vertices"] = Rcpp::transpose(vertex_mat);
    out["faces"]    = face_list;
    return out;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List readFileMesh_cpp(
  const std::string filename, const bool normals, const bool verbose) {
  Mesh3 mesh;
  const bool ok = CGAL::IO::read_polygon_mesh(
      filename, mesh, CGAL::parameters::verbose(verbose));
  if(!ok) {
    Rcpp::stop("Reading failure.");
  }
  if(verbose) {
    run_mesh_checks<MeshT>(mesh);
  }
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
void writeFile_cpp(const std::string filename,
                   const bool binary,
                   const int precision,
                   const Rcpp::NumericMatrix vertices,
                   const Rcpp::List faceList) {
  const std::vector<Point3> points = matrix_to_points3<Point3>(vertices);
  const std::pair<std::vector<std::vector<std::size_t>>, bool> faces =
      list_to_faces2(faceList);
  if(filename.length() < 5) {
      Rcpp::stop("`filename` needs at least 5 characters, including dot file-extension.");
  }
  const std::string ext = toLower(filename.substr(filename.length() - 4, 4));
  bool ok = false;
  if(ext == ".ply") {
    ok = CGAL::IO::write_PLY(
      filename, points, faces.first,
      CGAL::parameters::use_binary_mode(binary).stream_precision(precision));
  } else if(ext == ".stl") {
    if(!faces.second) {
      Rcpp::stop("STL files only accept triangular faces.");
    }
    ok = CGAL::IO::write_STL(
      filename, points, faces.first,
      CGAL::parameters::use_binary_mode(binary).stream_precision(precision));
  } else if(ext == ".obj") {
    ok = CGAL::IO::write_OBJ(
      filename, points, faces.first,
      CGAL::parameters::stream_precision(precision));
  } else if(ext == ".off") {
    ok = CGAL::IO::write_OFF(
      filename, points, faces.first,
      CGAL::parameters::stream_precision(precision));
  } else {
    Rcpp::stop("Unknown file extension.");
  }
  if(!ok) {
    Rcpp::stop("Failed to write file.");
  }
}
