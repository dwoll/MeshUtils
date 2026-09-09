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

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/io.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
void rmessage(std::string msg) {
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
std::vector<std::vector<std::size_t>> matrix_to_tfaces(
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
                   const unsigned int max_num_holes,
                   const bool verbose) {
    if(repair_soup) {
        PMP::repair_polygon_soup(points, faces);
    }
    const bool is_oriented = PMP::orient_polygon_soup(points, faces);
    if(!is_oriented) {
        Rcpp::warning("Polygon soup orientation failed.\n  Remove holes / self-intersections if present.");
    }
    // triangulate if necessary
    bool is_triangle = is_triangle_soup(faces);
    if(triangulate && !is_triangle) {
        is_triangle = PMP::triangulate_polygons(points, faces);
    }
    // check for self-intersections
    // problem: may not have self-intersections in polygon soup,
    // but may have self-intersections after turning into mesh
    // const bool soup_has_self_int = PMP::does_polygon_soup_self_intersect(points, faces);
    // remove self-intersections if necessary and possible
    // if(soup_has_self_int && remove_intersections && is_triangle) {
    //     remove_selfint_soup<KernelT, PointT>(points, faces, remove_method);
    // }
    // create mesh from polygon soup
    MeshT mesh;
    PMP::orient_polygon_soup(points, faces);
    PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
    // filling boundary holes if necessary, requested, and possible
    if(!CGAL::is_closed(mesh) && fill_holes && (max_num_holes > 0)) {
        // mesh is passed by reference and modified in fill_boundary_holes()
        // TODO also pass parameters related to hole size
        MeshT mesh_tmp = fill_boundary_holes<MeshT, PointT>(mesh, fair_hole, -1, -1, max_num_holes, verbose);
        mesh = std::move(mesh_tmp);
    }
    // remove self-intersections if necessary, requested, and possible
    if(PMP::does_self_intersect(mesh) && remove_intersections && is_triangle) {
        remove_selfint_mesh<KernelT, MeshT, PointT>(mesh, remove_method, verbose);
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
    const unsigned int,
    const bool);

template EMesh3 soup_to_mesh<EK, EMesh3, EPoint3>(
    std::vector<EPoint3>,
    std::vector<std::vector<std::size_t>>,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int,
    const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::polygon_soup_to_polygon_mesh() with fewer checks
// points and faces are changed -> no const, no reference
template <typename MeshT, typename PointT>
MeshT soup_to_mesh_valid(std::vector<PointT> points,
                         std::vector<std::vector<std::size_t>> faces,
                         const bool triangulate,
                         const bool repair_soup) {  // repair_soup_currently ignored (performance)
  // if(repair_soup) {
  //     PMP::repair_polygon_soup(points, faces);
  // }
  // triangulate if necessary
  const bool is_triangle = is_triangle_soup(faces);
  if(triangulate && !is_triangle) {
      PMP::triangulate_polygons(points, faces);
  }
  const bool is_oriented = PMP::orient_polygon_soup(points, faces);
  if(!is_oriented) {
      Rcpp::warning("Polygon soup orientation failed.\n  Remove holes / self-intersections if present.");
  }
  MeshT mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
  return mesh;
}

template Mesh3 soup_to_mesh_valid<Mesh3, Point3>(
  std::vector<Point3>, std::vector<std::vector<std::size_t>>, const bool, const bool);

template EMesh3 soup_to_mesh_valid<EMesh3, EPoint3>(
  std::vector<EPoint3>, std::vector<std::vector<std::size_t>>, const bool, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT vf_to_mesh(const Rcpp::NumericMatrix &vertices,
                 const Rcpp::List &faces,
                 const bool triangulate) {
  MeshT mesh;
  using face_descriptor = typename boost::graph_traits<MeshT>::face_descriptor;

  const std::size_t nVerts = vertices.ncol();
  for(std::size_t j = 0; j < nVerts; j++) {
    Rcpp::NumericVector vertex = vertices(Rcpp::_, j);
    PointT pt(vertex(0), vertex(1), vertex(2));
    mesh.add_vertex(pt);
  }
  const std::size_t nFaces = faces.size();
  for(std::size_t i = 0; i < nFaces; i++) {
    Rcpp::IntegerVector intface = Rcpp::as<Rcpp::IntegerVector>(faces(i));
    const std::size_t face_size = intface.size();
    std::vector<typename MeshT::Vertex_index> face;
    face.reserve(face_size);
    for(std::size_t k = 0; k < face_size; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface(k)));
    }
    face_descriptor fd = mesh.add_face(face);
    if(fd == mesh.null_face()) {
      Rcpp::stop("Cannot add face " + std::to_string(i+1) + ".");
    }
  }
  // triangulate if necessary and requested
  if(!CGAL::is_triangle_mesh(mesh) && triangulate) {
      PMP::triangulate_faces(mesh);
  }
  return mesh;
}

template Mesh3  vf_to_mesh<Mesh3,  Point3>(const Rcpp::NumericMatrix&,  const Rcpp::List&, const bool);
template EMesh3 vf_to_mesh<EMesh3, EPoint3>(const Rcpp::NumericMatrix&, const Rcpp::List&, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT>
void run_mesh_checks(MeshT mesh) {
    const bool is_triangle  = CGAL::is_triangle_mesh(mesh);
    const bool has_self_int = PMP::does_self_intersect(mesh);
    const bool is_closed    = CGAL::is_closed(mesh);
    if(is_closed) {
        rmessage("Mesh is closed.");
    } else {
        rmessage("Mesh is not closed.");
    }
    if(is_triangle) {
        rmessage("Mesh is triangle.");
    } else {
        rmessage("Mesh is not triangle. Cannot ensure it bounds a volume.");
    }
    if(is_triangle && is_closed) {
        if(!PMP::is_outward_oriented(mesh)) {
            PMP::reverse_face_orientations(mesh);
        }
        if(!has_self_int) {
            if(PMP::does_bound_a_volume(mesh)) {
                rmessage("Mesh bounds a volume.");
            } else {
                PMP::orient_to_bound_a_volume(mesh);
                if(!PMP::does_bound_a_volume(mesh)) {
                    rmessage("Mesh does not bound a volume (after trying).");
                }
            }
        }
    }
    if(has_self_int) {
        rmessage("Mesh has self-intersections. Mesh does not bound a volume.");
    } else {
        rmessage("Mesh does not have self-intersections.");
    }
    if(mesh.is_valid()) {
        rmessage("Mesh is valid.\n");
    } else {
        rmessage("Mesh is not valid.\n");
    }
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// general conversion from R list to Surface_mesh_3
// turn R data structures for vertices and faces into C++ vectors
// then call soup_to_mesh()
// faces may not be triangle -> list
template <typename KernelT, typename MeshT, typename PointT>
MeshT make_surf_mesh(
  const Rcpp::List &rmesh,
  const bool triangulate,
  const bool repair_soup,
  const bool remove_intersections,
  const int remove_method,
  const bool fill_holes,
  const bool fair_hole,
  const unsigned int max_num_holes,
  const bool verbose) {
  const Rcpp::NumericMatrix             rvertices = Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::List                      rfaces    = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<PointT>                   points    = matrix_to_points3<PointT>(rvertices);
  std::vector<std::vector<std::size_t>> faces     = list_to_faces1(rfaces);
  MeshT mesh = soup_to_mesh<KernelT, MeshT, PointT>(
      points,
      faces,
      triangulate,
      repair_soup,
      remove_intersections,
      remove_method,
      fill_holes,
      fair_hole,
      max_num_holes,
      verbose);

  if(verbose) {
    run_mesh_checks<MeshT>(mesh);
  }
  return mesh;
}

template Mesh3 make_surf_mesh<K, Mesh3, Point3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int,
    const bool);

template EMesh3 make_surf_mesh<EK, EMesh3, EPoint3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int,
    const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// read mesh file (polygon soup) and convert to Surface_mesh_3
// faces may not be triangle -> list
template <typename KernelT, typename MeshT, typename PointT>
MeshT make_surf_mesh_ff(
  const Rcpp::String filename,
  const bool triangulate,
  const bool repair_soup,
  const bool remove_intersections,
  const int remove_method,
  const bool fill_holes,
  const bool fair_hole,
  const unsigned int max_num_holes,
  const bool verbose) {
  std::vector<PointT> points;
  std::vector<std::vector<std::size_t>> polygons;
  const bool ok = CGAL::IO::read_polygon_soup(
      filename, points, polygons, CGAL::parameters::verbose(verbose));
  if(!ok) {
    Rcpp::stop("Reading failure.");
  }
  MeshT mesh = soup_to_mesh<KernelT, MeshT, PointT>(
      points,
      polygons,
      triangulate,
      repair_soup,
      remove_intersections,
      remove_method,
      fill_holes,
      fair_hole,
      max_num_holes,
      verbose);
  if(verbose) {
    run_mesh_checks<MeshT>(mesh);
  }
  return mesh;
}

template Mesh3 make_surf_mesh_ff<K, Mesh3, Point3>(
    const Rcpp::String,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int,
    const bool);

template EMesh3 make_surf_mesh_ff<EK, EMesh3, EPoint3>(
    const Rcpp::String,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int,
    const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// make_surf_tmesh()
// like make_surf_mesh() but for triangles -> rfaces is matrix
// const Rcpp::IntegerMatrix             rfaces = Rcpp::as<Rcpp::IntegerMatrix>(rmesh["faces"]);
// std::vector<std::vector<std::size_t>> faces  = matrix_to_tfaces(rfaces);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT make_surf_mesh_valid(const Rcpp::List &rmesh,
                           const bool soup,
                           const bool triangulate,
                           const bool repair_soup,
                           const bool verbose) {
    const Rcpp::NumericMatrix rvertices = Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
    const Rcpp::List          rfaces    = Rcpp::as<Rcpp::List>(rmesh["faces"]);
    MeshT mesh;
    if(soup) {
        MeshT mesh_tmp = soup_to_mesh_valid<MeshT, PointT>(
            matrix_to_points3<PointT>(rvertices),
            list_to_faces1(rfaces),
            triangulate,
            repair_soup);
        mesh = std::move(mesh_tmp);
    } else {
        MeshT mesh_tmp = vf_to_mesh<MeshT, PointT>(rvertices, rfaces, triangulate);
        mesh = std::move(mesh_tmp);
    }

    if(verbose) {
      run_mesh_checks<MeshT>(mesh);
    }
    return mesh;
}

template Mesh3 make_surf_mesh_valid<Mesh3,  Point3>(
    const Rcpp::List&, const bool, const bool, const bool, const bool);

template EMesh3 make_surf_mesh_valid<EMesh3, EPoint3>(
    const Rcpp::List&, const bool, const bool, const bool, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT make_surf_mesh_valid_ff(const Rcpp::String filename,
                              const bool soup,
                              const bool triangulate,
                              const bool repair_soup,
                              const bool verbose) {
    MeshT mesh;
    if(soup) {
      std::vector<PointT> points;
      std::vector<std::vector<std::size_t>> polygons;
      const bool ok = CGAL::IO::read_polygon_soup(
          filename, points, polygons, CGAL::parameters::verbose(verbose));
      if(!ok) {
        Rcpp::stop("Reading failure.");
      }
      MeshT mesh_tmp = soup_to_mesh_valid<MeshT, PointT>(
          points, polygons, triangulate, repair_soup);
      mesh = std::move(mesh_tmp);
    } else {
      MeshT mesh_tmp;
      const bool ok = PMP::IO::read_polygon_mesh(filename, mesh_tmp);
      if(!ok) {
        Rcpp::stop("Reading failure.");
      }
      mesh = std::move(mesh_tmp);
    }
    if(verbose) {
      run_mesh_checks<MeshT>(mesh);
    }
    return mesh;
}

template Mesh3 make_surf_mesh_valid_ff<Mesh3,  Point3>(
    const Rcpp::String, const bool, const bool, const bool, const bool);

template EMesh3 make_surf_mesh_valid_ff<EMesh3, EPoint3>(
    const Rcpp::String, const bool, const bool, const bool, const bool);
