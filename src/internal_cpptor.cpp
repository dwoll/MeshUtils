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

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
Rcpp::NumericMatrix getVertices(const MeshT &mesh) {
  const std::size_t nVerts = mesh.number_of_vertices();
  Rcpp::NumericMatrix Vertices(3, nVerts);
  {
    std::size_t i = 0;
    for(typename MeshT::Vertex_index vd : mesh.vertices()) {
      Rcpp::NumericVector col_i(3);
      const PointT vertex = mesh.point(vd);
      col_i(0) = CGAL::to_double<typename KernelT::FT>(vertex.x());
      col_i(1) = CGAL::to_double<typename KernelT::FT>(vertex.y());
      col_i(2) = CGAL::to_double<typename KernelT::FT>(vertex.z());
      Vertices(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Vertices;
}

template Rcpp::NumericMatrix getVertices<K,  Mesh3,  Point3>(const  Mesh3&);
template Rcpp::NumericMatrix getVertices<EK, EMesh3, EPoint3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges(const MeshT &mesh) {
  const std::size_t nEdges = mesh.number_of_edges();
  Rcpp::IntegerVector I1(nEdges);
  Rcpp::IntegerVector I2(nEdges);
  Rcpp::NumericVector Length(nEdges);
  Rcpp::NumericVector Angle(nEdges);
  Rcpp::LogicalVector Exterior(nEdges);
  Rcpp::LogicalVector Coplanar(nEdges);
  {
    std::size_t i = 0;
    for(typename MeshT::Edge_index ed : mesh.edges()) {
      typename MeshT::Vertex_index s = source(ed, mesh);
      typename MeshT::Vertex_index t = target(ed, mesh);
      I1(i) = static_cast<int>(s) + 1;
      I2(i) = static_cast<int>(t) + 1;
      std::vector<PointT> points(4);
      points[0] = mesh.point(s);
      points[1] = mesh.point(t);
      typename MeshT::Halfedge_index h0 = mesh.halfedge(ed, 0);
      points[2] = mesh.point(mesh.target(mesh.next(h0)));
      typename MeshT::Halfedge_index h1 = mesh.halfedge(ed, 1);
      points[3] = mesh.point(mesh.target(mesh.next(h1)));
      typename KernelT::FT angle = CGAL::abs(CGAL::approximate_dihedral_angle(
          points[0], points[1], points[2], points[3]));
      Angle(i) = CGAL::to_double<typename KernelT::FT>(angle);
      Exterior(i) = angle < 179.0 || angle > 181.0;
      Coplanar(i) = CGAL::coplanar(points[0], points[1], points[2], points[3]);
      typename KernelT::FT el = PMP::edge_length(h0, mesh);
      Length(i) = CGAL::to_double<typename KernelT::FT>(el);
      i++;
    }
  }
  Rcpp::DataFrame Edges = Rcpp::DataFrame::create(
    Rcpp::Named("i1")       = I1,
    Rcpp::Named("i2")       = I2,
    Rcpp::Named("length")   = Length,
    Rcpp::Named("angle")    = Angle,
    Rcpp::Named("exterior") = Exterior,
    Rcpp::Named("coplanar") = Coplanar
  );
  return Edges;
}

template Rcpp::DataFrame getEdges<K,  Mesh3,  Point3>(const  Mesh3&);
template Rcpp::DataFrame getEdges<EK, EMesh3, EPoint3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for possibly heterogeneous faces -> list
template <typename MeshT>
Rcpp::List getFaces1(const MeshT &mesh) {
  const std::size_t nFaces = mesh.number_of_faces();
  Rcpp::List face_list(nFaces);
  {
    std::size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      face_list(i) = col_i;
      i++;
    }
  }
  return face_list;
}

template Rcpp::List getFaces1<Mesh3>(const  Mesh3&);
template Rcpp::List getFaces1<EMesh3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for homogeneous faces -> matrix
template <typename MeshT>
Rcpp::IntegerMatrix getFaces2(const MeshT &mesh, const std::size_t nSides) {
  const std::size_t nFaces = mesh.number_of_faces();
  Rcpp::IntegerMatrix face_mat(nSides, nFaces);
  {
    std::size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      face_mat(Rcpp::_, i++) = col_i;
    }
  }
  return face_mat;
}

template Rcpp::IntegerMatrix getFaces2<Mesh3>(const Mesh3&,   const std::size_t);
template Rcpp::IntegerMatrix getFaces2<EMesh3>(const EMesh3&, const std::size_t);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh is changed (normals map added)
// -> no const, no ref (as const in calling function)
template <typename KernelT, typename MeshT, typename VectorT>
Rcpp::NumericMatrix computeVNormals(MeshT mesh) {
    using vertex_descriptor = typename boost::graph_traits<MeshT>::vertex_descriptor;
    Rcpp::NumericMatrix normals_mat(3, mesh.number_of_vertices());
    remove_properties<MeshT, VectorT>(mesh, {"v:normal"});
    auto vnormals = mesh.template add_property_map<vertex_descriptor, VectorT>(
                            "v:normal", CGAL::NULL_VECTOR).first;
    // PMP::compute_normals(mesh, vnormals, fnormals);
    PMP::compute_vertex_normals(mesh, vnormals);
    std::size_t i = 0;
    for(vertex_descriptor vd : vertices(mesh)) {
        Rcpp::NumericVector col_i(3);
        const VectorT normal = vnormals[vd];
        col_i(0) = CGAL::to_double<typename KernelT::FT>(normal.x());
        col_i(1) = CGAL::to_double<typename KernelT::FT>(normal.y());
        col_i(2) = CGAL::to_double<typename KernelT::FT>(normal.z());
        normals_mat(Rcpp::_, i) = col_i;
        i++;
    }
    return normals_mat;
}

template Rcpp::NumericMatrix computeVNormals<K,  Mesh3,  Vector3>(Mesh3);
template Rcpp::NumericMatrix computeVNormals<EK, EMesh3, EVector3>(EMesh3);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for possibly heterogeneous faces -> list
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh1(const MeshT &mesh, const bool normals) {
  Rcpp::DataFrame     Edges    = getEdges<KernelT, MeshT, PointT>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices<KernelT, MeshT, PointT>(mesh);
  Rcpp::List          Faces    = getFaces1<MeshT>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges")    = Edges,
                                      Rcpp::Named("faces")    = Faces);
  if(normals) {
    Rcpp::NumericMatrix vnormals_mat = computeVNormals<KernelT, MeshT, VectorT>(mesh);
    out["normals"] = vnormals_mat;
  }
  return out;
}

template Rcpp::List make_rmesh1<K,  Mesh3,  Point3,  Vector3>(const Mesh3&,   const bool);
template Rcpp::List make_rmesh1<EK, EMesh3, EPoint3, EVector3>(const EMesh3&, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for homogeneous faces -> matrix
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh2(const MeshT &mesh, const bool normals, const std::size_t nSides) {
  Rcpp::DataFrame     Edges    = getEdges<KernelT, MeshT, PointT>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices<KernelT, MeshT, PointT>(mesh);
  Rcpp::IntegerMatrix Faces    = getFaces2<MeshT>(mesh, nSides);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix vnormals_mat = computeVNormals<KernelT, MeshT, VectorT>(mesh);
    out["normals"] = vnormals_mat;
  }
  return out;
}

template Rcpp::List make_rmesh2<K,  Mesh3,  Point3,  Vector3>(const Mesh3&,   const bool, const std::size_t);
template Rcpp::List make_rmesh2<EK, EMesh3, EPoint3, EVector3>(const EMesh3&, const bool, const std::size_t);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::triangulate_faces() modifies -> mesh not const
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List get_rmesh(MeshT &mesh, const bool triangulate, const bool normals) {
  bool is_triangle = CGAL::is_triangle_mesh(mesh);
  if(triangulate && !is_triangle) {
    is_triangle = PMP::triangulate_faces(mesh);
  }
  Rcpp::List rmesh;
  if(is_triangle) {
    rmesh = make_rmesh2<KernelT, MeshT, PointT, VectorT>(mesh, normals, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    rmesh = make_rmesh2<KernelT, MeshT, PointT, VectorT>(mesh, normals, 4);
  } else {
    rmesh = make_rmesh1<KernelT, MeshT, PointT, VectorT>(mesh, normals);
  }

  return rmesh;
}

template Rcpp::List get_rmesh<K,  Mesh3,  Point3,  Vector3>(Mesh3&,   const bool, const bool);
template Rcpp::List get_rmesh<EK, EMesh3, EPoint3, EVector3>(EMesh3&, const bool, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// CAVE: pass by reference, modifies mesh
template <typename MeshT, typename VectorT>
void remove_properties(MeshT &mesh, const std::vector<std::string> &props) {
  using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
  using face_descriptor    = typename boost::graph_traits<MeshT>::face_descriptor;
  using vertex_colors_map  = typename MeshT::template Property_map<vertex_descriptor, std::string>;
  using face_colors_map    = typename MeshT::template Property_map<face_descriptor,   std::string>;
  using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;
  using vertex_scalars_map = typename MeshT::template Property_map<vertex_descriptor, double>;
  using face_scalars_map   = typename MeshT::template Property_map<face_descriptor,   double>;

  for(std::size_t i = 0; i < props.size(); i++) {
    std::string prop = props[i];
    if(prop == "v:color") {
      std::optional<vertex_colors_map> pmap_ =
        mesh.template property_map<vertex_descriptor, std::string>("v:color");
      if(pmap_.has_value()) {
        mesh.template remove_property_map(pmap_.value());
      }
    } else if(prop == "f:color") {
      std::optional<face_colors_map> pmap_ =
        mesh.template property_map<face_descriptor, std::string>("f:color");
      if(pmap_.has_value()) {
        mesh.template remove_property_map(pmap_.value());
      }
    } else if(prop == "v:normal") {
      std::optional<vertex_normals_map> pmap_ =
        mesh.template property_map<vertex_descriptor, VectorT>("v:normal");

      if(pmap_.has_value()) {
        mesh.template remove_property_map(pmap_.value());
      }
    } else if(prop == "v:scalar") {
      std::optional<vertex_scalars_map> pmap_ =
        mesh.template property_map<vertex_descriptor, double>("v:scalar");
      if(pmap_.has_value()) {
        mesh.template remove_property_map(pmap_.value());
      }
    } else if(prop == "f:scalar") {
      std::optional<face_scalars_map> pmap_ =
        mesh.template property_map<face_descriptor, double>("f:scalar");
      if(pmap_.has_value()) {
        mesh.template remove_property_map(pmap_.value());
      }
    }
  }
}

template void remove_properties<Mesh3,  Vector3>(Mesh3&,   const std::vector<std::string>&);
template void remove_properties<EMesh3, EVector3>(EMesh3&, const std::vector<std::string>&);
