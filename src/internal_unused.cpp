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
*/
// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
/*
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
*/
// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
/*
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
*/
// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
/*
Rcpp::NumericVector defaultNormal() {
  Rcpp::NumericVector def =
    {
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na()
    };
  return def;
}
*/
// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// compatibility wrapper for CGAL property_map(std::string) API changes:
// older returned std::pair<Property_map, bool>
// newer returns  std::optional<Property_map>
// property_map_pair returns a std::pair<Property_map, bool> in both cases
/*
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
*/
// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
/*
(const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_)
using norm_map_r   = typename MeshT::template Property_map<v_descriptor, Rcpp::NumericVector>;
using vertex_descriptor = typename boost::graph_traits<MeshT>::vertex_descriptor;
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
*/
