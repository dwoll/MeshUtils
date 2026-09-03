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

// ----------------------------------------------------------------------- //
#include "MeshUtils_types.h"

#include <CGAL/Vector_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;

// -------------------------------------------------------------------------- //
#define CGAL_EIGEN3_ENABLED 1
#define PIA_TAG CGAL::Parallel_if_available_tag

typedef std::pair<Point3, Vector3>                      P3V3;  // Point3 with normal Vector3
typedef boost::graph_traits<Mesh3>::face_descriptor     fc_dscrptr;
typedef boost::graph_traits<Mesh3>::edge_descriptor     dg_dscrptr;
typedef boost::graph_traits<Mesh3>::halfedge_descriptor hlfdg_dscrptr;

// -------------------------------------------------------------------------- //
template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix&);

template <typename KernelT, typename MeshT, typename PointT>
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

template <typename KernelT, typename MeshT, typename PointT>
MeshT make_surf_mesh(
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
MeshT make_surf_tmesh(
    const Rcpp::List&,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh1(const MeshT&, const bool);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh2(const MeshT&, const bool, const std::size_t);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List get_rmesh(MeshT&, const bool, const bool);

template <typename KernelT, typename MeshT, typename PointT>
MeshT remove_selfint_mesh(const MeshT&, const int);

template <typename MeshT, typename PointT>
MeshT fill_boundary_holes(
    MeshT&, const bool, const double, const int, const unsigned int);

template <typename MeshT, typename VectorT>
void remove_properties(MeshT&, const std::vector<std::string>&);

// -------------------------------------------------------------------------- //
// no template
void rmessage(std::string);
bool is_triangle_soup(const std::vector<std::vector<std::size_t>>&);

          std::vector<std::vector<std::size_t>>        list_to_faces1(const Rcpp::List&);
std::pair<std::vector<std::vector<std::size_t>>, bool> list_to_faces2(const Rcpp::List&);

// -------------------------------------------------------------------------- //
#endif
