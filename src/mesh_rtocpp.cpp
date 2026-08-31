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
// PMP::polygon_soup_to_polygon_mesh() with fewer checks
// points and faces are changed -> no const, no reference
// currently unused (only in unused makeMesh())
template <typename MeshT, typename PointT>
MeshT csoup_to_mesh(std::vector<PointT> points,
                    std::vector<std::vector<std::size_t>> faces,
                    const bool repairSoup) {
  if(repairSoup) {
    PMP::repair_polygon_soup(points, faces);
  }
  const bool success = PMP::orient_polygon_soup(points, faces);
  if(!success) {
    Rcpp::warning("Polygon orientation failed.");
  }
  MeshT mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
  if(!mesh.is_valid(false)) {
    Rcpp::warning("Mesh is not valid.");
  }
  return mesh;
}

template Mesh3 csoup_to_mesh<Mesh3, Point3>(
  std::vector<Point3>, std::vector<std::vector<std::size_t>>, const bool);

template EMesh3 csoup_to_mesh<EMesh3, EPoint3>(
  std::vector<EPoint3>, std::vector<std::vector<std::size_t>>, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT vf_to_mesh(const Rcpp::NumericMatrix &vertices,
                 const Rcpp::List &faces) {
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
    const std::size_t sf = intface.size();
    std::vector<typename MeshT::Vertex_index> face;
    face.reserve(sf);
    for(std::size_t k = 0; k < sf; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface(k)));
    }
    face_descriptor fd = mesh.add_face(face);
    if(fd == mesh.null_face()) {
      Rcpp::stop("Cannot add face " + std::to_string(i+1) + ".");
    }
  }
  return mesh;
}

template Mesh3  vf_to_mesh<Mesh3,  Point3>(const Rcpp::NumericMatrix&,  const Rcpp::List&);
template EMesh3 vf_to_mesh<EMesh3, EPoint3>(const Rcpp::NumericMatrix&, const Rcpp::List&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
Rcpp::NumericVector defaultNormal() {
  Rcpp::NumericVector def =
    {
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na()
    };
  return def;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// currently unused
template <typename MeshT, typename PointT>
EMesh3 makeMesh(const Rcpp::NumericMatrix &vertices,
                const Rcpp::List &faces,
                const bool soup,
                const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_) {
  using v_descriptor = typename boost::graph_traits<MeshT>::vertex_descriptor;
  using norm_map_r   = typename MeshT::template Property_map<v_descriptor, Rcpp::NumericVector>;
  if(soup) {
    return csoup_to_mesh<MeshT, PointT>(
        matrix_to_points3<PointT>(vertices),
        list_to_faces1(faces),
        true);
  }

  MeshT mesh = vf_to_mesh<MeshT, PointT>(vertices, faces);
  if(normals_.isNotNull()) {
    Rcpp::NumericMatrix normals_mat(normals_);
    const unsigned int nNormals = static_cast<unsigned int>(normals_mat.ncol());
    if(mesh.number_of_vertices() != nNormals) {
      Rcpp::stop(
        "The number of normals does not match the number of vertices.");
    }
    Rcpp::NumericVector def = defaultNormal();
    norm_map_r normalsmap =
      mesh.add_property_map<v_descriptor, Rcpp::NumericVector>(
        "v:normal", def).first;
    for(std::size_t j = 0; j < nNormals; j++) {
      Rcpp::NumericVector normal = normals_mat(Rcpp::_, j);
      normalsmap[CGAL::SM_Vertex_index(j)] = normal;
    }
  }
  return mesh;
}

template Mesh3  makeMesh<Mesh3,  Point3>(
    const Rcpp::NumericMatrix&,
    const Rcpp::List&, const bool, const Rcpp::Nullable<Rcpp::NumericMatrix> &);

template EMesh3 makeMesh<EMesh3, EPoint3>(
    const Rcpp::NumericMatrix&,
    const Rcpp::List&, const bool, const Rcpp::Nullable<Rcpp::NumericMatrix> &);

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

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// EPEC kernel only
// TODO template
bool is_small_hole(hlfdg_descriptor h,
                   const EMesh3 &mesh,
                   const double max_hole_diam,
                   const int max_num_hole_edges) {
  int num_hole_edges = 0;
  CGAL::Bbox_3 hole_bbox;
  for (hlfdg_descriptor hc : CGAL::halfedges_around_face(h, mesh)) {
    const EPoint3& p = mesh.point(target(hc, mesh));

    hole_bbox += p.bbox();
    ++num_hole_edges;

    // exit early, to avoid unnecessary traversal of large holes
    if (num_hole_edges > max_num_hole_edges) return false;
    if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
    if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
    if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
  }

  return true;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh is changed in function -> not const
// CAVE: pass by reference, modifies mesh
template <typename MeshT, typename PointT>
MeshT fillBoundaryHoles(
    MeshT &mesh,
    const bool fair_hole,
    const double max_hole_diam,
    const int max_num_hole_edges,
    const unsigned int max_num_holes) {
  using face_descriptor     = typename boost::graph_traits<MeshT>::face_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<MeshT>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<MeshT>::halfedge_descriptor;
  Message("Attempting to fill hole(s).");
  if(max_num_holes == 0) {
    Message("'max_num_holes' is 0. Nothing done.");
    return mesh;
  }
  PMP::remove_almost_degenerate_faces(mesh);
  std::vector<halfedge_descriptor> border_cycles;
  unsigned int nb_holes_ok   = 0;
  unsigned int nb_holes_fail = 0;
  CGAL::extract_boundary_cycles(mesh, std::back_inserter(border_cycles)); // requires CGAL 6.2
  // PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
  size_t n_border = border_cycles.size();
  if(n_border == 0) {
    Message("There's no border in this mesh. Nothing done.");
    return mesh;
  }

  // collect one halfedge per boundary cycle
  for(halfedge_descriptor h : border_cycles) {
    if((nb_holes_ok + nb_holes_fail) >= max_num_holes) {
        break;
    }
    /*
    if((max_hole_diam > 0)      &&
       (max_num_hole_edges > 0) &&
       !is_small_hole(h, mesh, max_hole_diam, max_num_hole_edges)) {
        continue;
    }
    */

    std::vector<face_descriptor>   patch_facets;
    std::vector<vertex_descriptor> patch_vertices;
    bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(
        mesh, h,
        CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
                       .vertex_output_iterator(std::back_inserter(patch_vertices))));
    if(success) {
        nb_holes_ok++;
    } else {
        nb_holes_fail++
    }
  }

  std::string msg1;
  msg1 = "Filled " + std::to_string(nb_holes_ok) + " boundary hole(s).";
  Message(msg1);
  if(nb_holes_fail > 0) {
      std::string msg2;
      msg2 = "Failed to fill " + std::to_string(nb_holes_fail) + " boundary hole(s).";
      Message(msg2);
  }

  std::vector<PointT> points;
  std::vector<std::vector<std::size_t>> polygons;
  PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
  // PMP::repair_polygon_soup(points, polygons);
  PMP::merge_duplicate_polygons_in_polygon_soup(
      points, polygons,
      CGAL::parameters::erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD)
  );
  MeshT mesh_out;
  PMP::orient_polygon_soup(points, polygons);
  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
  if(!CGAL::is_closed(mesh_out)) {
      Message("mesh_out not closed immediately after polygon_soup_to_polygon_mesh()");
  }
  return mesh_out;
}

template Mesh3  fillBoundaryHoles<Mesh3,  Point3>(Mesh3&,   const bool, const double, const int, const unsigned int);
template EMesh3 fillBoundaryHoles<EMesh3, EPoint3>(EMesh3&, const bool, const double, const int, const unsigned int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// CAVE: points, polygons are passed by reference, and are modified in place
// uses a polygon soup as container as the output will most likely be non-manifold
template <typename KernelT, typename PointT>
bool removeSelfIntSoup(std::vector<PointT> &points,
                       std::vector<std::vector<std::size_t>> &polygons,
                       const int method) {
    Message("Attempting to remove self-intersections.");
    bool success;
    if(method == 1) {
        success = PMP::autorefine_triangle_soup(points, polygons,
            CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag())
                .erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD)  // requires CGAL version 6.2
        );
    } else if(method == 2) {
        // TODO .snap_grid_size(grid_size).number_of_iterations(15));
        success = PMP::autorefine_triangle_soup(points, polygons,
            CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag())
                .apply_iterative_snap_rounding(true)
                .erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD)  // requires CGAL version 6.2
        );
    } else {
        Message("Wrong method. Needs to be 1 or 2. Nothing done.");
        return false;
    }
    if(success) {
        Message("Autorefine successful.");
    } else {
        Message("Autorefine not successful.");
    }

    if(PMP::does_polygon_soup_self_intersect(points, polygons)) {
        Message("Polygon soup still self-intersects after autorefine.");
    }
    CGAL::Conforming_constrained_Delaunay_triangulation_3<KernelT> ccdt;
    ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
    return success;
}

template bool removeSelfIntSoup<K, Point3>(
    std::vector<Point3>&, std::vector<std::vector<std::size_t>>&, const int);

template bool removeSelfIntSoup<EK, EPoint3>(
    std::vector<EPoint3>&, std::vector<std::vector<std::size_t>>&, const int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh passed by const reference, then uses a polygon soup
// as container as the output will most likely be non-manifold
template <typename KernelT, typename MeshT, typename PointT>
MeshT removeSelfIntMesh(const MeshT &mesh, const int method) {
  std::vector<PointT> points;
  std::vector<std::vector<std::size_t>> polygons;
  PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
  const bool success = removeSelfIntSoup<KernelT, PointT>(points, polygons, method);
  MeshT mesh_out;
  if(PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
  } else {
    Message("Polygon soup not a polygon mesh after removing intersections. Nothing done.");
    return mesh;
  }
  if(PMP::does_self_intersect(mesh_out)) {
    Message("Mesh self-intersections could not be removed.");
  }
  return mesh_out;
}

template Mesh3  removeSelfIntMesh<K,  Mesh3,  Point3>(const Mesh3&,   const int);
template EMesh3 removeSelfIntMesh<EK, EMesh3, EPoint3>(const EMesh3&, const int);
