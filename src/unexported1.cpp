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
template <typename MeshT, typename PointT>
MeshT soup_to_mesh(std::vector<PointT> points,
                   std::vector<std::vector<size_t>> faces,
                   const bool clean,
                   const bool triangulate) {
  if(clean) {
    PMP::repair_polygon_soup(points, faces);
  }
  bool success = PMP::orient_polygon_soup(points, faces);
  if(success) {
    Message("Successful polygon orientation.");
  } else {
    Message("Polygon orientation failed.");
    // check if polygons from the polygon soup are all triangles
    bool is_triangle_soup = std::all_of(
        faces.begin(),
        faces.end(),
        [](const auto& p) { return p.size() == 3; }
    );

    if(clean && is_triangle_soup) {
        bool success = PMP::autorefine_triangle_soup(points, faces);
        if(success) {
            Message("Intersections removed from triangle soup");
        } else {
            Message("Intersections could not be removed from triangle soup");
        }
        PMP::repair_polygon_soup(points, faces);
    }
  }
  MeshT mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
  const bool valid = mesh.is_valid(false);
  if(!valid) {
    Message("The mesh is not valid.");
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
  std::string msg1;
  if(!CGAL::is_closed(mesh)) {
      if(clean) {
          // boundary hole filling with EPEC kernel
          // need to work around mesh being potentially based on EPIC kernel
          EMesh3 mesh_epeck;
          CGAL::copy_face_graph(mesh, mesh_epeck);
          EMesh3 mesh_filled = fillBoundaryHoles(mesh_epeck, true, -1, -1);
          if(!CGAL::is_closed(mesh_filled)) {
              Rcpp::stop("The mesh is still not closed after trying to fill holes.");
          } else {
              MeshT mesh_copy;
              CGAL::copy_face_graph(mesh_filled, mesh_copy);
              mesh = mesh_copy;
              msg1 = " after filling holes";
          }
      } else {
          Rcpp::stop("The mesh is not closed.");
      }
  }
  Message("The mesh is closed" + msg1 + ".");
  if(isTriangle) {
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
  return mesh;
}

template Mesh3 soup_to_mesh<Mesh3, Point3>(
    std::vector<Point3>, std::vector<std::vector<size_t>>, const bool, const bool);

template EMesh3 soup_to_mesh<EMesh3, EPoint3>(
    std::vector<EPoint3>, std::vector<std::vector<size_t>>, const bool, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::polygon_soup_to_polygon_mesh() with fewer checks
// points and faces are changed -> no const, no reference
// currently unused (only in unused makeMesh())
template <typename MeshT, typename PointT>
MeshT csoup_to_mesh(std::vector<PointT> points,
                    std::vector<std::vector<size_t>> faces,
                    const bool clean) {
  if(clean) {
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
  const Rcpp::List &rmesh, const bool clean, const bool triangulate) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<size_t>> faces = list_to_faces(rfaces);
  return soup_to_mesh<MeshT, PointT>(points, faces, clean, triangulate);
}

template Mesh3  makeSurfMesh<Mesh3,  Point3>(const Rcpp::List&,  const bool, const bool);
template EMesh3 makeSurfMesh<EMesh3, EPoint3>(const Rcpp::List&, const bool, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT makeSurfTMesh(
  const Rcpp::List &rmesh, const bool clean, const bool triangulate) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::IntegerMatrix rfaces =
      Rcpp::as<Rcpp::IntegerMatrix>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<size_t>> faces = matrix_to_Tfaces(rfaces);
  return soup_to_mesh<MeshT, PointT>(points, faces, clean, triangulate);
}

template Mesh3  makeSurfTMesh<Mesh3,  Point3>(const Rcpp::List&,  const bool, const bool);
template EMesh3 makeSurfTMesh<EMesh3, EPoint3>(const Rcpp::List&, const bool, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// Compatibility wrapper for CGAL property_map API changes:
// Older CGAL returned std::pair<Property_map, bool>; newer returns std::optional<Property_map>.
// property_map_pair returns a std::pair<Property_map, bool> in both cases.
template <typename KeyT, typename T, typename MeshT>
std::pair<typename MeshT::template Property_map<KeyT,T>, bool>
property_map_pair(MeshT mesh, const std::string name) {
  using RetType = decltype(mesh.template property_map<KeyT,T>(name));
  using Pmap = typename MeshT::template Property_map<KeyT,T>;
  if constexpr (std::is_same_v<RetType, std::pair<Pmap, bool>>) {
    return mesh.template property_map<KeyT,T>(name);
  } else {
    auto opt = mesh.template property_map<KeyT,T>(name);
    if(opt) return std::make_pair(*opt, true);
    return std::make_pair(Pmap(), false);
  }
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
Rcpp::List getRmesh(const Mesh3 &mesh) {
  std::pair<nrmlsmap, bool> normalsmap_ =
      property_map_pair<vxdescr, Rcpp::NumericVector>(mesh, "v:normal");
  const bool there_is_normals = normalsmap_.second;
  Rcpp::List rmesh;
  if(CGAL::is_triangle_mesh(mesh)) {
    rmesh = RSurfMesh2<K, Mesh3, Point3, Vector3>(mesh, false, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    rmesh = RSurfMesh2<K, Mesh3, Point3, Vector3>(mesh, false, 4);
  } else {
    rmesh = RSurfMesh1<K, Mesh3, Point3, Vector3>(mesh, false);
  }
  if(there_is_normals) {
    nrmlsmap normalsmap = normalsmap_.first;
    Rcpp::NumericMatrix Normals(3, mesh.number_of_vertices());
    for(size_t i = 0; i < mesh.number_of_vertices(); i++) {
      Normals(Rcpp::_, i) = normalsmap[CGAL::SM_Vertex_index(i)];
    }
    rmesh["normals"] = Normals;
  }

  return rmesh;
}

Rcpp::List getRmesh(const EMesh3 &mesh) {
  std::pair<normals_map, bool> normalsmap_ =
      property_map_pair<vertex_descriptor, Rcpp::NumericVector>(mesh, "v:normal");
  const bool there_is_normals = normalsmap_.second;
  Rcpp::List rmesh;
  if(CGAL::is_triangle_mesh(mesh)) {
    rmesh = RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(mesh, false, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    rmesh = RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(mesh, false, 4);
  } else {
    rmesh = RSurfMesh1<EK, EMesh3, EPoint3, EVector3>(mesh, false);
  }
  if(there_is_normals) {
    normals_map normalsmap = normalsmap_.first;
    Rcpp::NumericMatrix Normals(3, mesh.number_of_vertices());
    for(size_t i = 0; i < mesh.number_of_vertices(); i++) {
      Normals(Rcpp::_, i) = normalsmap[CGAL::SM_Vertex_index(i)];
    }
    rmesh["normals"] = Normals;
  }

  return rmesh;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// EPEC kernel only
bool is_small_hole(halfedge_descriptor h, const EMesh3 &mesh,
                   double max_hole_diam, int max_num_hole_edges) {
  int num_hole_edges = 0;
  CGAL::Bbox_3 hole_bbox;
  for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh))
  {
    const EPoint3& p = mesh.point(target(hc, mesh));

    hole_bbox += p.bbox();
    ++num_hole_edges;

    // Exit early, to avoid unnecessary traversal of large holes
    if (num_hole_edges > max_num_hole_edges) return false;
    if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
    if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
    if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
  }

  return true;
}

// mesh is changed in function -> not const, not reference
EMesh3 fillBoundaryHoles(
    EMesh3 mesh, bool fairhole, double max_hole_diam, int max_num_hole_edges) {
  unsigned int nb_holes = 0;
  std::vector<halfedge_descriptor> border_cycles;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
  const int nBorders = border_cycles.size();
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
      if(fairhole) {
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
  }

  std::string msg;
  msg = "Filled " + std::to_string(nb_holes) + " boundary holes.";
  Message(msg);

  return mesh;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// TODO use property_map_pair() defined above
// and write analogous function remove_property_map_pair()
/*
void removeProperties(EMesh3& mesh, std::vector<std::string> props) {
  for(unsigned int i = 0; i < props.size(); i++) {
    std::string prop = props[i];
    if(prop == "f:color") {
      std::pair<Fcolors_map, bool> pmap_ =
        mesh.property_map<face_descriptor, std::string>("f:color");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:color") {
      std::pair<Vcolors_map, bool> pmap_ =
        mesh.property_map<vertex_descriptor, std::string>("v:color");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:normal") {
      std::pair<normals_map, bool> pmap_ =
        mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:scalar") {
      std::pair<Vscalars_map, bool> pmap_ =
        mesh.property_map<vertex_descriptor, double>("v:scalar");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "f:scalar") {
      std::pair<Fscalars_map, bool> pmap_ =
        mesh.property_map<face_descriptor, double>("f:scalar");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    }
  }
}
*/
