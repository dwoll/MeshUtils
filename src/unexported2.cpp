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
template <typename KernelT, typename MeshT, typename PointT>
Rcpp::NumericMatrix getVertices(const MeshT &mesh) {
  const size_t nVerts = mesh.number_of_vertices();
  Rcpp::NumericMatrix Vertices(3, nVerts);
  {
    size_t i = 0;
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

template Rcpp::NumericMatrix getVertices<K,  Mesh3,  Point3>(const Mesh3&);
template Rcpp::NumericMatrix getVertices<EK, EMesh3, EPoint3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges(const MeshT &mesh) {
  const size_t nEdges = mesh.number_of_edges();
  Rcpp::IntegerVector I1(nEdges);
  Rcpp::IntegerVector I2(nEdges);
  Rcpp::NumericVector Length(nEdges);
  Rcpp::NumericVector Angle(nEdges);
  Rcpp::LogicalVector Exterior(nEdges);
  Rcpp::LogicalVector Coplanar(nEdges);
  {
    size_t i = 0;
    for(typename MeshT::Edge_index ed : mesh.edges()) {
      typename MeshT::Vertex_index s = source(ed, mesh);
      typename MeshT::Vertex_index t = target(ed, mesh);
      I1(i) = (int)s + 1;
      I2(i) = (int)t + 1;
      std::vector<PointT> points(4);
      points[0] = mesh.point(s);
      points[1] = mesh.point(t);
      typename MeshT::Halfedge_index h0 = mesh.halfedge(ed, 0);
      points[2] = mesh.point(mesh.target(mesh.next(h0)));
      typename MeshT::Halfedge_index h1 = mesh.halfedge(ed, 1);
      points[3] = mesh.point(mesh.target(mesh.next(h1)));
      typename KernelT::FT angle = CGAL::abs(CGAL::approximate_dihedral_angle(
          points[0], points[1], points[2], points[3]));
      Angle(i) = CGAL::to_double(angle);
      Exterior(i) = angle < 179.0 || angle > 181.0;
      Coplanar(i) = CGAL::coplanar(points[0], points[1], points[2], points[3]);
      typename KernelT::FT el = PMP::edge_length(h0, mesh);
      Length(i) = CGAL::to_double(el);
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

template Rcpp::DataFrame getEdges<K,  Mesh3,  Point3>(const Mesh3&);
template Rcpp::DataFrame getEdges<EK, EMesh3, EPoint3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT>
Rcpp::List getFaces1(const MeshT &mesh) {
  const size_t nFaces = mesh.number_of_faces();
  Rcpp::List Faces(nFaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(i) = col_i;
      i++;
    }
  }
  return Faces;
}

template Rcpp::List getFaces1<Mesh3>(const Mesh3&);
template Rcpp::List getFaces1<EMesh3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT>
Rcpp::IntegerMatrix getFaces2(const MeshT &mesh, const int nSides) {
  const size_t nFaces = mesh.number_of_faces();
  Rcpp::IntegerMatrix Faces(nSides, nFaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      // bool TEST = CGAL::is_triangle(mesh.halfedge(fd), mesh);
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(Rcpp::_, i++) = col_i;
    }
  }
  return Faces;
}

template Rcpp::IntegerMatrix getFaces2<Mesh3>(const Mesh3&,   const int);
template Rcpp::IntegerMatrix getFaces2<EMesh3>(const EMesh3&, const int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename MeshT>
Rcpp::IntegerMatrix getTFaces(const MeshT &mesh) {
  const size_t nFaces = mesh.number_of_faces();
  Rcpp::IntegerMatrix Faces(3, nFaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Faces;
}

template Rcpp::IntegerMatrix getTFaces<Mesh3>(const Mesh3&);
template Rcpp::IntegerMatrix getTFaces<EMesh3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh used as this -> no const, no reference
Rcpp::NumericMatrix getNormals(Mesh3 mesh) {
  const size_t nVerts = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nVerts);
  auto vnormals = mesh.add_property_map<Mesh3::Vertex_index, Vector3>(
                          "v:normals", CGAL::NULL_VECTOR)
                      .first;
  auto fnormals = mesh.add_property_map<Mesh3::Face_index, Vector3>(
                          "f:normals", CGAL::NULL_VECTOR)
                      .first;
  PMP::compute_normals(mesh, vnormals, fnormals);
  {
    size_t i = 0;
    for(Mesh3::Vertex_index vd : vertices(mesh)) {
      Rcpp::NumericVector col_i(3);
      const Vector3 normal = vnormals[vd];
      col_i(0) = CGAL::to_double<K::FT>(normal.x());
      col_i(1) = CGAL::to_double<K::FT>(normal.y());
      col_i(2) = CGAL::to_double<K::FT>(normal.z());
      Normals(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Normals;
}

Rcpp::NumericMatrix getNormals(EMesh3 mesh) {
  const size_t nVerts = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nVerts);
  auto vnormals = mesh.add_property_map<EMesh3::Vertex_index, EVector3>(
                          "v:normals", CGAL::NULL_VECTOR)
                      .first;
  auto fnormals = mesh.add_property_map<EMesh3::Face_index, EVector3>(
                          "f:normals", CGAL::NULL_VECTOR)
                      .first;
  PMP::compute_normals(mesh, vnormals, fnormals);
  {
    size_t i = 0;
    for(EMesh3::Vertex_index vd : vertices(mesh)) {
      Rcpp::NumericVector col_i(3);
      const EVector3 normal = vnormals[vd];
      col_i(0) = CGAL::to_double<EK::FT>(normal.x());
      col_i(1) = CGAL::to_double<EK::FT>(normal.y());
      col_i(2) = CGAL::to_double<EK::FT>(normal.z());
      Normals(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Normals;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List RSurfMesh1(const MeshT &mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges<KernelT, MeshT, PointT>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices<KernelT, MeshT, PointT>(mesh);
  Rcpp::List Faces = getFaces1<MeshT>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

template Rcpp::List RSurfMesh1<K,  Mesh3,  Point3,  Vector3>(const Mesh3&,   const bool);
template Rcpp::List RSurfMesh1<EK, EMesh3, EPoint3, EVector3>(const EMesh3&, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List RSurfMesh2(const MeshT &mesh, const bool normals, const int nSides) {
  Rcpp::DataFrame Edges = getEdges<KernelT, MeshT, PointT>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices<KernelT, MeshT, PointT>(mesh);
  Rcpp::IntegerMatrix Faces = getFaces2<MeshT>(mesh, nSides);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

template Rcpp::List RSurfMesh2<K,  Mesh3,  Point3,  Vector3>(const Mesh3&,   const bool, const int);
template Rcpp::List RSurfMesh2<EK, EMesh3, EPoint3, EVector3>(const EMesh3&, const bool, const int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List RSurfTMesh(const MeshT &mesh, const bool normals) {
  Rcpp::DataFrame Edges = getEdges<KernelT, MeshT, PointT>(mesh);
  Rcpp::NumericMatrix Vertices = getVertices<KernelT, MeshT, PointT>(mesh);
  Rcpp::IntegerMatrix Faces = getTFaces<MeshT>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("edges") = Edges,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}

template Rcpp::List RSurfTMesh<K,  Mesh3,  Point3,  Vector3>(const Mesh3&,   const bool);
template Rcpp::List RSurfTMesh<EK, EMesh3, EPoint3, EVector3>(const EMesh3&, const bool);
