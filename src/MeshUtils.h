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
// #include <CGAL/squared_distance_2.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

// -------------------------------------------------------------------------- //
typedef CGAL::Advancing_front_surface_reconstruction<>     AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3                AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2 AFS_Tds2;

typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>      Polyhedron;
typedef CGAL::Parallel_if_available_tag                              Concurrency_tag;
typedef CGAL::Face_filtered_graph<EMesh3>                            Filtered_graph;

typedef CGAL::Scale_space_surface_reconstruction_3<K>                SSS_reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother<K> SSS_smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher<K>    SSS_mesher;
typedef SSS_reconstruction::Facet_const_iterator                     SSS_facet_iterator;
typedef SSS_reconstruction::Point_const_iterator                     SSS_point_iterator;

typedef std::pair<Point3, Vector3> P3wn;
typedef boost::graph_traits<Mesh3>::vertex_descriptor       vxdescr;
typedef boost::graph_traits<EMesh3>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<Mesh3>::face_descriptor         fdescr;
typedef boost::graph_traits<EMesh3>::face_descriptor        face_descriptor;
typedef boost::graph_traits<EMesh3>::halfedge_descriptor    halfedge_descriptor;
typedef EMesh3::Property_map<vxdescr, Rcpp::NumericVector>  Normals_map;
typedef EMesh3::Property_map<face_descriptor, std::size_t>  Face_index_map;

// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;

// -------------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT soup_to_mesh(std::vector<PointT>, std::vector<std::vector<size_t>>, const bool, const bool triangulate);

template <typename MeshT, typename PointT>
MeshT csoup_to_mesh(std::vector<PointT>, std::vector<std::vector<size_t>>, const bool);

std::vector<std::vector<size_t>> list_to_faces(const Rcpp::List);

template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix);

template <typename PointT>
Rcpp::NumericMatrix points3_to_matrix(std::vector<PointT>);

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(const Rcpp::List, const bool, const bool);

template <typename MeshT, typename PointT>
MeshT  makeSurfTMesh(const Rcpp::List, const bool, const bool);

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges(MeshT);

Rcpp::NumericMatrix getEKNormals(EMesh3);

Rcpp::List getRmesh(EMesh3);
Rcpp::List RSurfEKMesh(EMesh3, const bool);
Rcpp::List RSurfEKMesh2(EMesh3, const bool, const int);
Rcpp::List RSurfTEKMesh(EMesh3, const bool);

void Message(std::string);

EMesh3 fillBoundaryHoles(EMesh3, bool, double, int);
EMesh3 removeSelfIntersections(const EMesh3, const unsigned int);

// -------------------------------------------------------------------------- //
#endif
