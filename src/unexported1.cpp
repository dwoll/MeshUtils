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
  for(size_t i = 0; i < nPts; i++) {
    const Rcpp::NumericVector pt = M(Rcpp::_, i);
    points.emplace_back(PointT(pt(0), pt(1), pt(2)));
  }
  return points;
}

template std::vector<Point3>  matrix_to_points3<Point3>(const Rcpp::NumericMatrix&);
template std::vector<EPoint3> matrix_to_points3<EPoint3>(const Rcpp::NumericMatrix&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename PointT>
Rcpp::NumericMatrix points3_to_matrix(const std::vector<PointT> &points) {
  const size_t nPts = points.size();
  Rcpp::NumericMatrix M(3, nPts);
  for(size_t i = 0; i != nPts; i++) {
    Rcpp::NumericVector col_i(3);
    const PointT point = points[i];
    col_i(0) = CGAL::to_double(point.x());
    col_i(1) = CGAL::to_double(point.y());
    col_i(2) = CGAL::to_double(point.z());
    M(Rcpp::_, i) = col_i;
  }
  return M;
}

template Rcpp::NumericMatrix points3_to_matrix<Point3>(const std::vector<Point3>&);
template Rcpp::NumericMatrix points3_to_matrix<EPoint3>(const std::vector<EPoint3>&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
std::vector<std::vector<size_t>> matrix_to_Tfaces(
  const Rcpp::IntegerMatrix &faceMat) {
  const size_t nFaces = faceMat.ncol();
  std::vector<std::vector<size_t>> faces;
  faces.reserve(nFaces);
  for(size_t i = 0; i < nFaces; i++) {
    const Rcpp::IntegerVector face_rcpp = faceMat(Rcpp::_, i);
    // need static cast here because initializing with {} instead of ()
    std::vector<size_t> face = { static_cast<size_t>(face_rcpp(0)),
                                 static_cast<size_t>(face_rcpp(1)),
                                 static_cast<size_t>(face_rcpp(2)) };
    faces.emplace_back(face);
  }
  return faces;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
std::vector<std::vector<size_t>> list_to_faces(const Rcpp::List &L) {
  const size_t nFaces = L.size();
  std::vector<std::vector<size_t>> faces;
  faces.reserve(nFaces);
  for(size_t i = 0; i < nFaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    std::vector<size_t> face(face_rcpp.begin(), face_rcpp.end());
    faces.emplace_back(face);
  }
  return faces;
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
                   const unsigned int remove_method,
                   const bool fill_holes,
                   const bool fair_hole,
                   const unsigned int max_num_holes) {
  if(repair_soup) {
    PMP::repair_polygon_soup(points, faces);
  }
  bool success = PMP::orient_polygon_soup(points, faces);
  if(success) {
    Message("Successful polygon orientation.");
  } else {
    Rcpp::warning("Polygon orientation failed.",
                  "Try to check for, and then remove self-intersections.");
  }
  bool hasSelfIntersections = PMP::does_polygon_soup_self_intersect(points, faces));
  if(hasSelfIntersections) {
      Message("Polygon soup has self-intersections.");
  }
  if(hasSelfIntersections && remove_intersections) {
    // check if polygons from the polygon soup are all triangles
    bool is_triangle_soup = std::all_of(
      faces.begin(),
      faces.end(),
      [](const auto& p) { return p.size() == 3; });
    if(is_triangle_soup) {
        Message("Attempting to remove self-intersections.");
        removeSelfIntersections<KernelT, MeshT, PointT>(points, faces, remove_method)
    } else {
        Rcpp::warning("Not a triangle mesh, could not try to remove self-intersections.");
    }
  }

  MeshT mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
  const bool isValid = mesh.is_valid(false);
  if(!isValid) {
    Rcpp::warning("The mesh is not valid.");
  }
  bool isTriangle = CGAL::is_triangle_mesh(mesh);
  if(triangulate && !isTriangle) {
    isTriangle = PMP::triangulate_faces(mesh);
  }
  if(isTriangle) {
    Message("The mesh is triangle.");
  } else {
    Message("The mesh is not triangle; no way to ensure it "
            "bounds a volume, and whether it is outward oriented.");
  }
  if(!CGAL::is_closed(mesh) && fill_holes && (max_num_holes > 0)) {
      Message("Attempting to fill hole(s).");
      // boundary hole filling with EPEC kernel
      // need to work around mesh being potentially based on EPIC kernel
      EMesh3 mesh_epeck;
      CGAL::copy_face_graph(mesh, mesh_epeck);
      // mesh_epeck is passed by reference and modified in fillBoundaryHoles()
      fillBoundaryHoles(mesh_epeck, fair_hole, -1, -1, max_num_holes);
      MeshT mesh_new;
      CGAL::copy_face_graph(mesh_epeck, mesh_new);
      mesh = mesh_new;
  }
  if(CGAL::is_closed(mesh)) {
      Message("The mesh is closed.");
  } else {
      Message("The mesh is not closed.");
  }
  if(isTriangle) {
    if(!PMP::is_outward_oriented(mesh)) {
      PMP::reverse_face_orientations(mesh);
    }
    const bool bv = PMP::does_bound_a_volume(mesh);
    if(bv) {
      Message("The mesh bounds a volume.");
    } else {
      Message("The mesh does not bound a volume - reorienting.");
      PMP::orient_to_bound_a_volume(mesh);
    }
  }
  return mesh;
}

template Mesh3 soup_to_mesh<Mesh3, Point3>(
    std::vector<Point3>,
    std::vector<std::vector<std::size_t>>,
    const bool,
    const bool,
    const bool,
    const unsigned int,
    const bool,
    const bool,
    const unsigned int);

template EMesh3 soup_to_mesh<EMesh3, EPoint3>(
    std::vector<EPoint3>,
    std::vector<std::vector<std::size_t>>,
    const bool,
    const bool,
    const bool,
    const unsigned int,
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
                    std::vector<std::vector<size_t>> faces,
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
  const bool valid = mesh.is_valid(false);
  if(!valid) {
    Rcpp::warning("The mesh is not valid.");
  }
  return mesh;
}

template Mesh3 csoup_to_mesh<Mesh3, Point3>(
  std::vector<Point3>, std::vector<std::vector<size_t>>, const bool);

template EMesh3 csoup_to_mesh<EMesh3, EPoint3>(
  std::vector<EPoint3>, std::vector<std::vector<size_t>>, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT vf_to_mesh(const Rcpp::NumericMatrix &vertices,
                 const Rcpp::List &faces) {
  MeshT mesh;
  const size_t nVerts = vertices.ncol();
  for(size_t j = 0; j < nVerts; j++) {
    Rcpp::NumericVector vertex = vertices(Rcpp::_, j);
    PointT pt(vertex(0), vertex(1), vertex(2));
    mesh.add_vertex(pt);
  }
  const size_t nFaces = faces.size();
  for(size_t i = 0; i < nFaces; i++) {
    Rcpp::IntegerVector intface = Rcpp::as<Rcpp::IntegerVector>(faces(i));
    const size_t sf = intface.size();
    std::vector<typename MeshT::Vertex_index> face;
    face.reserve(sf);
    for(size_t k = 0; k < sf; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface(k)));
    }
    typename boost::graph_traits<MeshT>::face_descriptor fd = mesh.add_face(face);
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
// TODO template, decide whether vertices and faces can be reference
// challenge: normals_map, vertex_descriptor vs. nrmlsmap, vxdescr
// currently unused
EMesh3 makeMesh(const Rcpp::NumericMatrix vertices,
                const Rcpp::List faces,
                bool soup,
                const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_) {
  if(soup) {
    return csoup_to_mesh<EMesh3, EPoint3>(matrix_to_points3<EPoint3>(vertices),
                                          list_to_faces(faces),
                                          true);
  }

  EMesh3 mesh = vf_to_mesh<EMesh3, EPoint3>(vertices, faces);
  if(normals_.isNotNull()) {
    Rcpp::NumericMatrix normals(normals_);
    if(mesh.number_of_vertices() != normals.ncol()) {
      Rcpp::stop(
        "The number of normals does not match the number of vertices.");
    }
    Rcpp::NumericVector def = defaultNormal();
    normals_map normalsmap =
      mesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
        "v:normal", def).first;
    for(int j = 0; j < normals.ncol(); j++) {
      Rcpp::NumericVector normal = normals(Rcpp::_, j);
      normalsmap[CGAL::SM_Vertex_index(j)] = normal;
    }
  }
  return mesh;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT makeSurfMesh(
  const Rcpp::List &rmesh,
  const bool triangulate,
  const bool repair_soup,
  const bool remove_intersections,
  const unsigned int remove_method,
  const bool fill_holes,
  const bool fair_hole,
  const unsigned int max_num_holes) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<size_t>> faces = list_to_faces(rfaces);
  return soup_to_mesh<MeshT, PointT>(
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

template Mesh3  makeSurfMesh<Mesh3,  Point3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const unsigned int,
    const bool,
    const bool,
    const unsigned int);

template EMesh3 makeSurfMesh<EMesh3, EPoint3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const unsigned int,
    const bool,
    const bool,
    const unsigned int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT makeSurfTMesh(
    const Rcpp::List &rmesh,
    const bool triangulate,
    const bool repair_soup,
    const bool remove_intersections,
    const unsigned int remove_method,
    const bool fill_holes,
    const bool fair_hole,
    const unsigned int max_num_holes) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::IntegerMatrix rfaces =
      Rcpp::as<Rcpp::IntegerMatrix>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<size_t>> faces = matrix_to_Tfaces(rfaces);
  return soup_to_mesh<MeshT, PointT>(
      points,
      triangulate,
      repair_soup,
      remove_intersections,
      remove_method,
      fill_holes,
      fair_hole,
      max_num_holes);
}

template Mesh3 makeSurfTMesh<Mesh3, Point3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const unsigned int,
    const bool,
    const bool,
    const unsigned int);

template EMesh3 makeSurfTMesh<EMesh3, EPoint3>(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const unsigned int,
    const bool,
    const bool,
    const unsigned int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// compatibility wrapper for CGAL property_map(std::string) API changes:
// older returned std::pair<Property_map, bool>
// newer returns  std::optional<Property_map>
// property_map_pair returns a std::pair<Property_map, bool> in both cases
template <typename KeyT, typename T, typename MeshT>
std::pair<typename MeshT::template Property_map<KeyT,T>, bool>
property_map_pair(MeshT mesh, const std::string name) {
  using RetType = decltype(mesh.template property_map<KeyT,T>(name));
  using Pmap = typename MeshT::template Property_map<KeyT,T>;
  if constexpr (std::is_same_v<RetType, std::pair<Pmap, bool>>) {
    return mesh.template property_map<KeyT,T>(name);
  } else {
    auto opt = mesh.template property_map<KeyT,T>(name);
    if(opt) {
        return std::make_pair(*opt, true);
    }
    return std::make_pair(Pmap(), false);
  }
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::triangulate_faces() modifies -> mesh not const
Rcpp::List getRmesh(Mesh3 &mesh, const bool triangulate) {
  std::pair<nrmlsmap, bool> normalsmap_ =
      property_map_pair<vxdescr, Rcpp::NumericVector, Mesh3>(mesh, "v:normal");
  const bool has_normals = normalsmap_.second;
  bool isTriangle = CGAL::is_triangle_mesh(mesh);
  if(triangulate && !isTriangle) {
    isTriangle = PMP::triangulate_faces(mesh);
  }
  Rcpp::List rmesh;
  if(isTriangle) {
    rmesh = RSurfMesh2<K, Mesh3, Point3, Vector3>(mesh, false, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    rmesh = RSurfMesh2<K, Mesh3, Point3, Vector3>(mesh, false, 4);
  } else {
    rmesh = RSurfMesh1<K, Mesh3, Point3, Vector3>(mesh, false);
  }
  if(has_normals) {
    nrmlsmap normalsmap = normalsmap_.first;
    Rcpp::NumericMatrix normals_mat(3, mesh.number_of_vertices());
    for(size_t i = 0; i < mesh.number_of_vertices(); i++) {
      normals_mat(Rcpp::_, i) = normalsmap[CGAL::SM_Vertex_index(i)];
    }
    rmesh["normals"] = normals_mat;
  }

  return rmesh;
}

// PMP::triangulate_faces() modifies -> mesh not const
Rcpp::List getRmesh(EMesh3 &mesh, const bool triangulate) {
  std::pair<normals_map, bool> normalsmap_ =
      property_map_pair<vertex_descriptor, Rcpp::NumericVector, EMesh3>(mesh, "v:normal");
  const bool has_normals = normalsmap_.second;
  bool isTriangle = CGAL::is_triangle_mesh(mesh);
  if(triangulate && !isTriangle) {
    isTriangle = PMP::triangulate_faces(mesh);
  }
  Rcpp::List rmesh;
  if(isTriangle) {
    rmesh = RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(mesh, false, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    rmesh = RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(mesh, false, 4);
  } else {
    rmesh = RSurfMesh1<EK, EMesh3, EPoint3, EVector3>(mesh, false);
  }
  if(has_normals) {
    normals_map normalsmap = normalsmap_.first;
    Rcpp::NumericMatrix normals_mat(3, mesh.number_of_vertices());
    for(size_t i = 0; i < mesh.number_of_vertices(); i++) {
      normals_mat(Rcpp::_, i) = normalsmap[CGAL::SM_Vertex_index(i)];
    }
    rmesh["normals"] = normals_mat;
  }

  return rmesh;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// EPEC kernel only
// max_num_hole_edges is int (and not unsigned int) because -1 is used
// as a flag in a calling function
bool is_small_hole(halfedge_descriptor h,
                   const EMesh3 &mesh,
                   const double max_hole_diam,
                   const int max_num_hole_edges) {
  int num_hole_edges = 0;
  CGAL::Bbox_3 hole_bbox;
  for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh))
  {
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

// mesh is changed in function -> not const
// CAVE: pass by reference, modifies mesh
void fillBoundaryHoles(
    EMesh3 &mesh,
    const bool fairHole,
    const double max_hole_diam,
    const int max_num_hole_edges,
    const unsigned int max_num_holes) {
  if(max_num_holes == 0) {
    Message("'max_num_holes' is 0 - nothing done.");
    return;
  }
  unsigned int nb_holes = 0;
  std::vector<halfedge_descriptor> border_cycles;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
  const size_t nBorders = border_cycles.size();
  if(nBorders == 0) {
    Rcpp::stop("There's no border in this mesh.");
  }

  for(halfedge_descriptor h : border_cycles) {
      if((max_hole_diam > 0)      &&
         (max_num_hole_edges > 0) &&
         !is_small_hole(h, mesh, max_hole_diam, max_num_hole_edges)) {
        continue;
      }

      std::vector<face_descriptor>   patch_faces;
      std::vector<vertex_descriptor> patch_vertices;
      if(fairHole) {
        const bool success = std::get<0>(
          PMP::triangulate_refine_and_fair_hole(
            mesh, h,
            CGAL::parameters::face_output_iterator(
              std::back_inserter(patch_faces)
            ).vertex_output_iterator(std::back_inserter(patch_vertices))
          )
        );
        if(!success) {
          Message("Fairing failed.");
        }
      } else {
        PMP::triangulate_and_refine_hole(
          mesh, h,
          CGAL::parameters::face_output_iterator(
            std::back_inserter(patch_faces)
          ).vertex_output_iterator(std::back_inserter(patch_vertices))
        );
      }
      ++nb_holes;
      if(nb_holes >= max_num_holes) {
        break;
      }
  }

  std::string msg;
  msg = "Filled " + std::to_string(nb_holes) + " boundary holes.";
  Message(msg);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// TODO template
// use property_map_pair() as defined above
// CAVE: pass by reference, modifies mesh
void removeProperties(EMesh3 &mesh, const std::vector<std::string> props) {
  for(size_t i = 0; i < props.size(); i++) {
    std::string prop = props[i];
    if(prop == "v:color") {
      std::pair<vcolors_map, bool> pmap_ =
         property_map_pair<vertex_descriptor, std::string, EMesh3>(mesh, "v:color");
        // mesh.property_map<vertex_descriptor, std::string>("v:color");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "f:color") {
      std::pair<fcolors_map, bool> pmap_ =
        property_map_pair<face_descriptor, std::string, EMesh3>(mesh, "f:color");
        // mesh.property_map<face_descriptor, std::string>("f:color");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:normal") {
      std::pair<normals_map, bool> pmap_ =
        property_map_pair<vertex_descriptor, Rcpp::NumericVector, EMesh3>(mesh, "v:normal");
        // mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:scalar") {
      std::pair<vscalars_map, bool> pmap_ =
        property_map_pair<vertex_descriptor, double, EMesh3>(mesh, "v:scalar");
        // mesh.property_map<vertex_descriptor, double>("v:scalar");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "f:scalar") {
      std::pair<fscalars_map, bool> pmap_ =
        property_map_pair<face_descriptor, double, EMesh3>(mesh, "f:scalar");
        // mesh.property_map<face_descriptor, double>("f:scalar");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    }
  }
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// CAVE: points, polygons are passed by reference, and are modified in place
// uses a polygon soup as container as the output will most likely be non-manifold
template <typename KernelT, typename MeshT, typename PointT>
void removeSelfIntersections(std::vector<PointT> &points,
                             std::vector<std::vector<std::size_t>> &polygons,
                             const unsigned int method) {
  bool success;
  if(method == 1) {
      success = PMP::autorefine_triangle_soup(points, polygons);
  } else if(method == 2) {
      const auto& snap = CGAL::parameters::apply_iterative_snap_rounding(true);
      success = PMP::autorefine_triangle_soup(points, polygons, snap);
  } else {
      Rcpp::stop("Wrong method");
  }
  if(success) {
      Message("Autorefine successful.\n");
  } else {
      Message("Autorefine not successful.\n");
  }

  CGAL::Conforming_constrained_Delaunay_triangulation_3<KernelT> ccdt;
  ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
  if(PMP::does_polygon_soup_self_intersect(points, polygons)) {
    Message("Self intersections could not be removed.\n");
  } else {
    Message("Self intersections removed.\n");
  }
}

template void removeSelfIntersections<K,  Mesh3,  Point3>(std::vector<PointT>&,  std::vector<std::vector<std::size_t>>&, const unsigned int);
template void removeSelfIntersections<EK, EMesh3, EPoint3>(std::vector<PointT>&, std::vector<std::vector<std::size_t>>&, const unsigned int);
