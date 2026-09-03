/*

typedef boost::graph_traits<Mesh3>::vertex_descriptor                vrtx_dscrptr;
typedef Mesh3::Property_map<vrtx_dscrptr, Rcpp::NumericVector>       nrmls_map_r;

typedef boost::graph_traits<EMesh3>::vertex_descriptor               vrtx_descriptor;
typedef EMesh3::Property_map<vrtx_descriptor, Rcpp::NumericVector>   normals_map_r;

typedef boost::graph_traits<EMesh3>::edge_descriptor                 dg_descriptor;
typedef boost::graph_traits<EMesh3>::halfedge_descriptor             hlfdg_descriptor;

// EPoint3 with normal EVector3
typedef std::pair<EPoint3, EVector3>                                 EP3EV3;
typedef boost::graph_traits<EMesh3>::face_descriptor                 fc_descriptor;

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT csoup_to_mesh(
    std::vector<PointT>, std::vector<std::vector<std::size_t>>, const bool);

template <typename KernelT, typename PointT>
Rcpp::NumericMatrix points3_to_matrix(const std::vector<PointT>&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename VectorT>
std::optional<Rcpp::NumericMatrix> getVNormals(const MeshT &mesh) {
    using vertex_descriptor = typename boost::graph_traits<MeshT>::vertex_descriptor;
    using vnormals_map      = typename MeshT::template Property_map<vertex_descriptor, VectorT>;

    std::optional<Rcpp::NumericMatrix> normals_mat;
    std::optional<vnormals_map> vnormals =
        mesh.template property_map<vertex_descriptor, VectorT>("v:normal");
    if(vnormals.has_value()) {
        Rcpp::NumericMatrix nm(3, mesh.number_of_vertices());
        std::size_t i = 0;
        for(vertex_descriptor vd : vertices(mesh)) {
          Rcpp::NumericVector col_i(3);
          const VectorT normal = vnormals.value()[vd];
          col_i(0) = CGAL::to_double<typename KernelT::FT>(normal.x());
          col_i(1) = CGAL::to_double<typename KernelT::FT>(normal.y());
          col_i(2) = CGAL::to_double<typename KernelT::FT>(normal.z());
          nm(Rcpp::_, i) = col_i;
          i++;
        }
        normals_mat = std::move(nm);
    }
    return normals_mat;
}

template std::optional<Rcpp::NumericMatrix> getVNormals<K,  Mesh3,  Vector3>(const  Mesh3&);
template std::optional<Rcpp::NumericMatrix> getVNormals<EK, EMesh3, EVector3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// currently not used
template <typename MeshT, typename PointT>
bool is_small_hole(typename boost::graph_traits<MeshT>::halfedge_descriptor h,
                   const MeshT &mesh,
                   const double max_hole_diam,
                   const int max_num_hole_edges) {
  using halfedge_descriptor = typename boost::graph_traits<MeshT>::halfedge_descriptor;
  int num_hole_edges = 0;
  CGAL::Bbox_3 hole_bbox;
  for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh)) {
    const PointT& p = mesh.point(target(hc, mesh));

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

template bool is_small_hole<Mesh3,  Point3>(typename  boost::graph_traits<Mesh3>::halfedge_descriptor,  const Mesh3&,  const double, const int);
template bool is_small_hole<EMesh3, EPoint3>(typename boost::graph_traits<EMesh3>::halfedge_descriptor, const EMesh3&, const double, const int);

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
// compatibility wrapper for CGAL property_map(std::string) API changes:
// older returned std::pair<Property_map, bool>
// newer returns  std::optional<Property_map>
// property_map_pair returns a std::pair<Property_map, bool> in both cases
template <typename KeyT, typename T, typename MeshT>
std::pair<typename MeshT::template Property_map<KeyT,T>, bool>
property_map_pair(MeshT &mesh, const std::string name) {
  using RetType = decltype(mesh.template property_map<KeyT,T>(name));
  using Pmap    = typename MeshT::template Property_map<KeyT,T>;
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
template <typename KernelT, typename PointT>
Rcpp::NumericMatrix points3_to_matrix(const std::vector<PointT> &points) {
  const std::size_t nPts = points.size();
  Rcpp::NumericMatrix M(3, nPts);
  for(std::size_t i = 0; i < nPts; i++) {
    Rcpp::NumericVector col_i(3);
    const PointT point = points[i];
    col_i(0) = CGAL::to_double<typename KernelT::FT>(point.x());
    col_i(1) = CGAL::to_double<typename KernelT::FT>(point.y());
    col_i(2) = CGAL::to_double<typename KernelT::FT>(point.z());
    M(Rcpp::_, i) = col_i;
  }
  return M;
}

template Rcpp::NumericMatrix points3_to_matrix<K,  Point3>(const  std::vector<Point3>&);
template Rcpp::NumericMatrix points3_to_matrix<EK, EPoint3>(const std::vector<EPoint3>&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::polygon_soup_to_polygon_mesh() with fewer checks
// points and faces are changed -> no const, no reference
// used in makeMesh())
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
template <typename MeshT, typename PointT>
MeshT makeMesh(const Rcpp::NumericMatrix &vertices,
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

  // TODO
  // normals_mat should be C++ vector
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
      mesh.template add_property_map<v_descriptor, Rcpp::NumericVector>(
        "v:normal", def).first;
    for(std::size_t j = 0; j < nNormals; j++) {
      Rcpp::NumericVector normal = normals_mat(Rcpp::_, j);
      normalsmap[CGAL::SM_Vertex_index(j)] = normal;
    }
  }
  return mesh;
}

template Mesh3 makeMesh<Mesh3,  Point3>(
    const Rcpp::NumericMatrix&,
    const Rcpp::List&, const bool, const Rcpp::Nullable<Rcpp::NumericMatrix> &);

template EMesh3 makeMesh<EMesh3, EPoint3>(
    const Rcpp::NumericMatrix&,
    const Rcpp::List&, const bool, const Rcpp::Nullable<Rcpp::NumericMatrix> &);
*/
