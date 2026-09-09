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

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh is changed in function -> not const
// CAVE: pass by reference, modifies mesh
template <typename MeshT, typename PointT>
MeshT fill_boundary_holes(
    MeshT &mesh,
    const bool fair_hole,
    const double max_hole_diam,
    const int max_num_hole_edges,
    const unsigned int max_num_holes,
    const bool verbose) {
  using face_descriptor     = typename boost::graph_traits<MeshT>::face_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<MeshT>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<MeshT>::halfedge_descriptor;
  if(verbose) { rmessage("Trying to fill boundary holes."); }
  if(max_num_holes == 0) {
    if(verbose) { rmessage("'max_num_holes' is 0. Nothing done."); }
    return mesh;
  }
  // PMP::remove_almost_degenerate_faces(mesh);
  std::vector<halfedge_descriptor> border_cycles;
  unsigned int nb_holes_ok   = 0;
  unsigned int nb_holes_fail = 0;
  // requires CGAL 6.2 (was PMP::extract_...)
  CGAL::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
  size_t n_border = border_cycles.size();
  if(n_border == 0) {
    if(verbose) { rmessage("There's no border in this mesh. Nothing done."); }
    return mesh;
  }

  // collect one halfedge per boundary cycle
  for(halfedge_descriptor h : border_cycles) {
    if((nb_holes_ok + nb_holes_fail) >= max_num_holes) {
        break;
    }
    // if((max_hole_diam > 0)      &&
    //    (max_num_hole_edges > 0) &&
    //    !is_small_hole<MeshT, PointT>(h, mesh, max_hole_diam, max_num_hole_edges)) {
    //     continue;
    // }
    std::vector<face_descriptor>   patch_facets;
    std::vector<vertex_descriptor> patch_vertices;
    bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(
        mesh, h,
        CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
                       .vertex_output_iterator(std::back_inserter(patch_vertices))));
    if(success) {
        nb_holes_ok++;
    } else {
        nb_holes_fail++;
    }
  }

  std::string msg1;
  msg1 = "Filled " + std::to_string(nb_holes_ok) + " boundary hole(s).";
  if(verbose) { rmessage(msg1); }
  if(nb_holes_fail > 0) {
      std::string msg2;
      msg2 = "Failed to fill " + std::to_string(nb_holes_fail) + " boundary hole(s).";
      if(verbose) { rmessage(msg2); }
  }

  std::vector<PointT> points;
  std::vector<std::vector<std::size_t>> polygons;
  PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
  // could do autorefine_triangle_soup() here, but apparently not necessary
  PMP::merge_duplicate_polygons_in_polygon_soup(
      points, polygons,
      // without the following option, self-intersections are created
      // requires CGAL 6.2
      CGAL::parameters::erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD));
  PMP::remove_isolated_points_in_polygon_soup(points, polygons);
  MeshT mesh_out;
  const bool orient_ok = PMP::orient_polygon_soup(points, polygons);
  if(!orient_ok) {
      Rcpp::warning("Polygon orientation after filling holes failed.");
  }
  if(!PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
      Rcpp::warning("Polygon soup not a mesh - non-manifold vertex? - need PMP::autorefine_triangle_soup()?");
  }
  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
  // PMP::merge_duplicated_vertices_in_boundary_cycles(mesh);
  // PMP::remove_isolated_vertices(mesh);
  // std::size_t nNmv = PMP::duplicate_non_manifold_vertices(mesh_out);
  // std::string msg_nmv;
  // rmessage("Number of non-manifold vertices: " + std::to_string(nNmv) + ".");
  // if(!CGAL::is_closed(mesh_out)) {
  //     rmessage("Failed to fill boundary hole(s).");
  // }
  // if(PMP::does_self_intersect(mesh_out)) {
  //     rmessage("Mesh self-intersects after filling boundary hole(s).");
  // }
  return mesh_out;
}

template Mesh3  fill_boundary_holes<Mesh3,  Point3>(Mesh3&,   const bool, const double, const int, const unsigned int, const bool);
template EMesh3 fill_boundary_holes<EMesh3, EPoint3>(EMesh3&, const bool, const double, const int, const unsigned int, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// CAVE: points, polygons are passed by reference, and are modified in place
// uses a polygon soup as container as the output will most likely be non-manifold
template <typename KernelT, typename PointT>
bool remove_selfint_soup(std::vector<PointT> &points,
                       std::vector<std::vector<std::size_t>> &polygons,
                       const int method,
                       const bool verbose) {
    bool success;
    std::string msg;
    std::string msg_method;
    if(verbose) { rmessage("Trying to remove self-intersections."); }
    if(method == 1) {
        msg_method = "";
        success = PMP::autorefine_triangle_soup(points, polygons,
            CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag())
                // requires CGAL version 6.2
                .erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD));
    } else if(method == 2) {
        msg_method=" with snap rounding";
        // TODO .snap_grid_size(grid_size).number_of_iterations(15));
        success = PMP::autorefine_triangle_soup(points, polygons,
            CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag())
                .apply_iterative_snap_rounding(true)
                // requires CGAL version 6.2
                .erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD));
    } else {
        msg_method = "";
        Rcpp::warning("Wrong method. Needs to be 1 (`auto`) or 2 (`auto_snap`). Nothing done.");
        return false;
    }
    if(!success) {
        msg = "Autorefine " + msg_method + " not successful.";
        Rcpp::warning(msg);
    }

    // autorefine_triangle_soup() can remove edges, put isolated vertices may remain
    // result may be non-manifold
    // PMP::repair_polygon_soup(points, polygons);
    PMP::merge_duplicate_points_in_polygon_soup(points, polygons,
        // without the following option, self-intersections are created
        // requires CGAL 6.2
        CGAL::parameters::erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD)
    );
    PMP::merge_duplicate_polygons_in_polygon_soup(
        points, polygons,
        // without the following option, self-intersections are created
        // requires CGAL 6.2
        CGAL::parameters::erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD)
    );
    PMP::remove_isolated_points_in_polygon_soup(points, polygons);
    if(PMP::does_polygon_soup_self_intersect(points, polygons)) {
        Rcpp::warning("Polygon soup still self-intersects after autorefine.");
    }
    // TODO what to do with this?
    // CGAL::Conforming_constrained_Delaunay_triangulation_3<KernelT> ccdt;
    // ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
    return success;
}

template bool remove_selfint_soup<K, Point3>(
    std::vector<Point3>&, std::vector<std::vector<std::size_t>>&, const int, const bool);

template bool remove_selfint_soup<EK, EPoint3>(
    std::vector<EPoint3>&, std::vector<std::vector<std::size_t>>&, const int, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh passed by const reference, then uses a polygon soup
// as container as the output will most likely be non-manifold
template <typename KernelT, typename MeshT, typename PointT>
MeshT remove_selfint_mesh(const MeshT &mesh, const int method, const bool verbose) {
  std::vector<PointT> points;
  std::vector<std::vector<std::size_t>> polygons;
  PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
  const bool success = remove_selfint_soup<KernelT, PointT>(points, polygons, method, verbose);
  MeshT mesh_out;
  PMP::orient_polygon_soup(points, polygons);
  if(PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
  } else {
    Rcpp::warning("Polygon soup not a polygon mesh after removing intersections. Nothing done.");
    // PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
    // PMP::merge_duplicated_vertices_in_boundary_cycles(mesh_out);
    // PMP::duplicate_non_manifold_vertices(mesh_out);
    // PMP::stitch_borders(mesh_out);
    return mesh;
  }
  if(!mesh_out.is_valid()) {
      Rcpp::warning("Mesh not valid after removing intersections.");
  }
  if(PMP::does_self_intersect(mesh_out)) {
    Rcpp::warning("Mesh self-intersections could not be removed.");
  }
  return mesh_out;
}

template Mesh3  remove_selfint_mesh<K,  Mesh3,  Point3>(const Mesh3&,   const int, const bool);
template EMesh3 remove_selfint_mesh<EK, EMesh3, EPoint3>(const EMesh3&, const int, const bool);
