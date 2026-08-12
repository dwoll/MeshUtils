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

void Message(std::string msg) {
  SEXP rmsg = Rcpp::wrap(msg);
  Rcpp::message(rmsg);
}

template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix M) {
  const size_t npoints = M.ncol();
  std::vector<PointT> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt = M(Rcpp::_, i);
    points.emplace_back(PointT(pt(0), pt(1), pt(2)));
  }
  return points;
}

template std::vector<Point3> matrix_to_points3<Point3>(
    const Rcpp::NumericMatrix);
template std::vector<EPoint3> matrix_to_points3<EPoint3>(
    const Rcpp::NumericMatrix);

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
template <typename PointT>
Rcpp::NumericMatrix points3_to_matrix(std::vector<PointT> points) {
  const size_t npoints = points.size();
  Rcpp::NumericMatrix M(3, npoints);
  for(size_t i = 0; i != npoints; i++) {
    Rcpp::NumericVector col_i(3);
    const PointT point = points[i];
    col_i(0) = CGAL::to_double(point.x());
    col_i(1) = CGAL::to_double(point.y());
    col_i(2) = CGAL::to_double(point.z());
    M(Rcpp::_, i) = col_i;
  }
  return M;
}

template Rcpp::NumericMatrix points3_to_matrix<EPoint3>(std::vector<EPoint3>);

std::vector<std::vector<size_t>> matrix_to_Tfaces(
    const Rcpp::IntegerMatrix Faces) {
  const size_t nfaces = Faces.ncol();
  std::vector<std::vector<size_t>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    const Rcpp::IntegerVector face_rcpp = Faces(Rcpp::_, i);
    std::vector<size_t> face = { face_rcpp(0), face_rcpp(1), face_rcpp(2) };
    faces.emplace_back(face);
  }
  return faces;
}

std::vector<std::vector<size_t>> list_to_faces(const Rcpp::List L) {
  const size_t nfaces = L.size();
  std::vector<std::vector<size_t>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    std::vector<size_t> face(face_rcpp.begin(), face_rcpp.end());
    faces.emplace_back(face);
  }
  return faces;
}

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
        Message("Trying to remove intersections from triangle soup.");
        bool success = PMP::autorefine_triangle_soup(points, faces);
        if(success) {
            Message("Intersections removed");
        } else {
            Message("Intersections could not be removed");
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
    Message("Triangulation.");
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Message("Triangulation has failed.");
    }
    isTriangle = true;
  }
  if(isTriangle) {
    Message("The mesh is triangle.");
  } else {
    Message(
        "The mesh is not triangle; no way to ensure it bounds a volume "
        "and whether it is outward oriented.");
  }
  if(!CGAL::is_closed(mesh)) {
      if(clean) {
          MeshT mesh_filled = fillBoundaryHoles(mesh, true, -1, -1);
          if(!CGAL::is_closed(mesh_filled)) {
              Rcpp::stop("The mesh is still not closed after trying to fill holes.");
          } else {
              mesh = mesh_filled;
              Message("The mesh is now closed after filling holes.");
          }
      } else {
          Rcpp::stop("The mesh is not closed.");
      }
  }
  Message("The mesh is closed.");
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

template EMesh3 csoup_to_mesh<EMesh3, EPoint3>(
  std::vector<EPoint3>, std::vector<std::vector<size_t>>, const bool
);

EMesh3 vf_to_emesh(const Rcpp::NumericMatrix vertices,
                   const Rcpp::List faces) {
  EMesh3 mesh;
  const int nv = vertices.ncol();
  for(int j = 0; j < nv ; j++) {
    Rcpp::NumericVector vertex = vertices(Rcpp::_, j);
    EPoint3 pt = EPoint3(vertex(0), vertex(1), vertex(2));
    mesh.add_vertex(pt);
  }
  const int nf = faces.size();
  for(int i = 0; i < nf; i++) {
    Rcpp::IntegerVector intface = Rcpp::as<Rcpp::IntegerVector>(faces(i));
    const int sf = intface.size();
    std::vector<EMesh3::Vertex_index> face;
    face.reserve(sf);
    for(int k = 0; k < sf; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface(k)));
    }
    face_descriptor fd = mesh.add_face(face);
    if(fd == EMesh3::null_face()) {
      Rcpp::stop("Cannot add face " + std::to_string(i+1) + ".");
    }
  }
  return mesh;
}

Rcpp::NumericVector defaultNormal() {
  Rcpp::NumericVector def =
    {
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na()
    };
  return def;
}

EMesh3 makeMesh(const Rcpp::NumericMatrix vertices,
                const Rcpp::List faces,
                bool soup,
                const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_) {
  if(soup) {
    return csoup_to_mesh<EMesh3, EPoint3>(
      matrix_to_points3<EPoint3>(vertices),
      list_to_faces(faces),
      true
    );
  }

  EMesh3 mesh = vf_to_emesh(vertices, faces);
  if(normals_.isNotNull()) {
    Rcpp::NumericMatrix normals(normals_);
    if(mesh.number_of_vertices() != normals.ncol()) {
      Rcpp::stop(
        "The number of normals does not match the number of vertices."
      );
    }
    Rcpp::NumericVector def = defaultNormal();
    Normals_map normalsmap =
      mesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
        "v:normal", def
      ).first;
    for(int j = 0; j < normals.ncol(); j++) {
      Rcpp::NumericVector normal = normals(Rcpp::_, j);
      normalsmap[CGAL::SM_Vertex_index(j)] = normal;
    }
  }
  return mesh;
}

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate
) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<size_t>> faces = list_to_faces(rfaces);
  return soup_to_mesh<MeshT, PointT>(points, faces, clean, triangulate);
}

template EMesh3 makeSurfMesh<EMesh3, EPoint3>(
  const Rcpp::List, const bool, const bool
);

template <typename MeshT, typename PointT>
MeshT makeSurfTMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate
) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::IntegerMatrix rfaces =
      Rcpp::as<Rcpp::IntegerMatrix>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<size_t>> faces = matrix_to_Tfaces(rfaces);
  return soup_to_mesh<MeshT, PointT>(points, faces, clean, triangulate);
}

template EMesh3 makeSurfTMesh<EMesh3, EPoint3>(
  const Rcpp::List, const bool, const bool
);

// Compatibility wrapper for CGAL property_map API changes:
// Older CGAL returned std::pair<Property_map, bool>; newer returns std::optional<Property_map>.
// property_map_pair returns a std::pair<Property_map, bool> in both cases.
template <typename Key, typename T, typename Mesh>
std::pair<typename Mesh::template Property_map<Key,T>, bool>
property_map_pair(Mesh mesh, const std::string name) {
  using RetType = decltype(mesh.template property_map<Key,T>(name));
  using Pmap = typename Mesh::template Property_map<Key,T>;
  if constexpr (std::is_same_v<RetType, std::pair<Pmap, bool>>) {
    return mesh.template property_map<Key,T>(name);
  } else {
    auto opt = mesh.template property_map<Key,T>(name);
    if(opt) return std::make_pair(*opt, true);
    return std::make_pair(Pmap(), false);
  }
}

// ----------------------------------------------------------------------- //
Rcpp::List getRmesh(EMesh3 mesh) {
  std::pair<Normals_map, bool> normalsmap_ =
      property_map_pair<vertex_descriptor, Rcpp::NumericVector>(mesh, "v:normal");
  const bool there_is_normals = normalsmap_.second;
  Rcpp::List rmesh;
  if(CGAL::is_triangle_mesh(mesh)) {
    rmesh = RSurfEKMesh2(mesh, false, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    rmesh = RSurfEKMesh2(mesh, false, 4);
  } else {
    rmesh = RSurfEKMesh(mesh, false);
  }
  if(there_is_normals) {
    Normals_map normalsmap = normalsmap_.first;
    Rcpp::NumericMatrix Normals(3, mesh.number_of_vertices());
    for(int i = 0; i < mesh.number_of_vertices(); i++) {
      Normals(Rcpp::_, i) = normalsmap[CGAL::SM_Vertex_index(i)];
    }
    rmesh["normals"] = Normals;
  }

  return rmesh;
}

// ----------------------------------------------------------------------- //
bool is_small_hole(halfedge_descriptor h, EMesh3 & mesh,
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

// TODO makes this work for Mesh3, template
EMesh3 fillBoundaryHoles(EMesh3 mesh, bool fairhole, double max_hole_diam, int max_num_hole_edges) {
  unsigned int nb_holes = 0;
  std::vector<halfedge_descriptor> border_cycles;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
  const int nborders = border_cycles.size();
  if(nborders == 0) {
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
EMesh3 removeSelfIntersections(const EMesh3 mesh, const unsigned int method) {
    EMesh3 mesh_out;
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
        if(method == 1) {
            CGAL::Conforming_constrained_Delaunay_triangulation_3<EK> ccdt;

            // use a polygon soup as container as the output will most likely be non-manifold
            std::vector<EPoint3> points;
            std::vector<std::vector<std::size_t>> polygons;
            PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
            bool success = PMP::autorefine_triangle_soup(points, polygons);
            // PMP::does_polygon_soup_self_intersect(points, polygons))
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
            mesh_out = meshrsi;
        } else if(method == 2) {
            CGAL::Conforming_constrained_Delaunay_triangulation_3<EK> ccdt;

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
              Message("Self intersections could not be removed (with itertive snap).\n");
            } else {
              Message("Self intersections removed (with iterative snap).\n");
            }
            mesh_out = meshrsi;
        }
    } else {
      Message("Mesh does not self-intersect\n");
      mesh_out = mesh;
    }

    return mesh_out;
}
