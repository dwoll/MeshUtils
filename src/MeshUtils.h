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
#define _CGALMESHHEADER_

#include <Rcpp.h>

#define CGAL_EIGEN3_ENABLED 1
#define PIA_TAG CGAL::Parallel_if_available_tag

// ----------------------------------------------------------------------- //
#include "MeshUtils_types.h"

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

// TODO CatmullClark, DooSabin, Sqrt3 subdivision
#include <CGAL/Subdivision_method_3/subdivision_methods_3.h>

// -------------------------------------------------------------------------- //
typedef CGAL::Advancing_front_surface_reconstruction<>     AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3                AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2 AFS_Tds2;

typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>      Polyhedron;
typedef CGAL::Parallel_if_available_tag                              Concurrency_tag;

typedef CGAL::Scale_space_surface_reconstruction_3<K>                SSS_reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother<K> SSS_smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher<K>    SSS_mesher;
typedef SSS_reconstruction::Facet_const_iterator                     SSS_facet_iterator;
typedef SSS_reconstruction::Point_const_iterator                     SSS_point_iterator;

// Point3 with normal vector
typedef std::pair<Point3, Vector3>                          P3wn;
typedef boost::graph_traits<Mesh3>::vertex_descriptor       vxdescr;
typedef boost::graph_traits<Mesh3>::face_descriptor         fdescr;
typedef boost::graph_traits<Mesh3>::edge_descriptor         edgdescr;
typedef boost::graph_traits<Mesh3>::halfedge_descriptor     hedgdescr;
typedef Mesh3::Property_map<vxdescr, Rcpp::NumericVector>   nrmlsmap;

// EPoint3 with normal vector
typedef std::pair<EPoint3, EVector3>                        EP3wn;
typedef boost::graph_traits<EMesh3>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<EMesh3>::face_descriptor        face_descriptor;
typedef boost::graph_traits<EMesh3>::edge_descriptor        edge_descriptor;
typedef boost::graph_traits<EMesh3>::halfedge_descriptor    halfedge_descriptor;
typedef EMesh3::Property_map<vertex_descriptor, Rcpp::NumericVector> normals_map;
typedef EMesh3::Property_map<vertex_descriptor, std::string> vcolors_map;
typedef EMesh3::Property_map<face_descriptor, std::string>   fcolors_map;
typedef EMesh3::Property_map<vertex_descriptor, double>      vscalars_map;
typedef EMesh3::Property_map<face_descriptor, double>        fscalars_map;

// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;

// -------------------------------------------------------------------------- //
template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix&);

template <typename KernelT, typename PointT>
Rcpp::NumericMatrix points3_to_matrix(const std::vector<PointT>&);

template <typename MeshT, typename PointT>
MeshT soup_to_mesh(
    std::vector<PointT>,                    // points
    std::vector<std::vector<std::size_t>>,  // faces
    const bool,                             // triangulate
    const bool,                             // repair_soup
    const bool,                             // remove_intersections
    const int,                              // remove_method
    const bool,                             // fill_holes
    const bool,                             // fair_hole
    const unsigned int);                    // max_num_holes

template <typename MeshT, typename PointT>
MeshT csoup_to_mesh(std::vector<PointT>, std::vector<std::vector<std::size_t>>, const bool);

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

// no bool triangulate as triangle mesh assumed
template <typename MeshT, typename PointT>
MeshT makeSurfTMesh(
    const Rcpp::List&,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges(const MeshT&);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List RSurfMesh1(const MeshT&, const bool);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List RSurfMesh2(const MeshT&, const bool, const std::size_t);

// -------------------------------------------------------------------------- //
// TODO template
Rcpp::List getRmesh(Mesh3&,  const bool);
Rcpp::List getRmesh(EMesh3&, const bool);

// template <typename KernelT, typename MeshT, typename VectorT>
// Rcpp::NumericMatrix getNormals(MeshT);
Rcpp::NumericMatrix getNormals(Mesh3);
Rcpp::NumericMatrix getNormals(EMesh3);

// -------------------------------------------------------------------------- //
// no template
void Message(std::string);

std::vector<std::vector<std::size_t>> list_to_faces(const Rcpp::List&);

// use only EPEC kernel
void fillBoundaryHoles(EMesh3&, const bool, const double, const int, const unsigned int);

// use only EPEC kernel
void removeSelfIntSoup(
    std::vector<EPoint3>&, std::vector<std::vector<std::size_t>>&, const int);

// use only EPEC kernel
EMesh3 removeSelfIntMesh(const EMesh3&, const int);

// currently only used with EPEC kernel
void removeProperties(EMesh3&, const std::vector<std::string>);

// -------------------------------------------------------------------------- //
#endif
