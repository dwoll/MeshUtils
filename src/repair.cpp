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

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh is changed in function -> not const
// CAVE: pass by reference, modifies mesh
template <typename MeshT, typename PointT>
MeshT fillBoundaryHoles(
    MeshT &mesh,
    const bool fair_hole,
    const double max_hole_diam,
    const int max_num_hole_edges,
    const unsigned int max_num_holes) {
  using face_descriptor     = typename boost::graph_traits<MeshT>::face_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<MeshT>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<MeshT>::halfedge_descriptor;
  Message("Attempting to fill hole(s).");
  if(max_num_holes == 0) {
    Message("'max_num_holes' is 0. Nothing done.");
    return mesh;
  }
  // PMP::remove_almost_degenerate_faces(mesh);
  std::vector<halfedge_descriptor> border_cycles;
  unsigned int nb_holes_ok   = 0;
  unsigned int nb_holes_fail = 0;
  CGAL::extract_boundary_cycles(mesh, std::back_inserter(border_cycles)); // requires CGAL 6.2
  // PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
  size_t n_border = border_cycles.size();
  if(n_border == 0) {
    Message("There's no border in this mesh. Nothing done.");
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
  Message(msg1);
  if(nb_holes_fail > 0) {
      std::string msg2;
      msg2 = "Failed to fill " + std::to_string(nb_holes_fail) + " boundary hole(s).";
      Message(msg2);
  }

  std::vector<PointT> points;
  std::vector<std::vector<std::size_t>> polygons;
  PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
  // PMP::repair_polygon_soup(points, polygons);
  PMP::merge_duplicate_polygons_in_polygon_soup(
      points, polygons,
      CGAL::parameters::erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD) // requires CGAL 6.2
  );
  MeshT mesh_out;
  PMP::orient_polygon_soup(points, polygons);
  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
  if(!CGAL::is_closed(mesh_out)) {
      Message("mesh_out not closed immediately after polygon_soup_to_polygon_mesh()");
  }
  return mesh_out;
}

template Mesh3  fillBoundaryHoles<Mesh3,  Point3>(Mesh3&,   const bool, const double, const int, const unsigned int);
template EMesh3 fillBoundaryHoles<EMesh3, EPoint3>(EMesh3&, const bool, const double, const int, const unsigned int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// CAVE: points, polygons are passed by reference, and are modified in place
// uses a polygon soup as container as the output will most likely be non-manifold
template <typename KernelT, typename PointT>
bool removeSelfIntSoup(std::vector<PointT> &points,
                       std::vector<std::vector<std::size_t>> &polygons,
                       const int method) {
    Message("Attempting to remove self-intersections.");
    bool success;
    if(method == 1) {
        success = PMP::autorefine_triangle_soup(points, polygons,
            CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag())
                .erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD)  // requires CGAL version 6.2
        );
    } else if(method == 2) {
        // TODO .snap_grid_size(grid_size).number_of_iterations(15));
        success = PMP::autorefine_triangle_soup(points, polygons,
            CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag())
                .apply_iterative_snap_rounding(true)
                .erase_policy(PMP::Duplicate_polygon_erase_policy::KEEP_ONE_IF_ODD)  // requires CGAL version 6.2
        );
    } else {
        Message("Wrong method. Needs to be 1 or 2. Nothing done.");
        return false;
    }
    if(success) {
        Message("Autorefine successful.");
    } else {
        Message("Autorefine not successful.");
    }

    if(PMP::does_polygon_soup_self_intersect(points, polygons)) {
        Message("Polygon soup still self-intersects after autorefine.");
    }
    // TODO what to do with this?
    // CGAL::Conforming_constrained_Delaunay_triangulation_3<KernelT> ccdt;
    // ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
    return success;
}

template bool removeSelfIntSoup<K, Point3>(
    std::vector<Point3>&, std::vector<std::vector<std::size_t>>&, const int);

template bool removeSelfIntSoup<EK, EPoint3>(
    std::vector<EPoint3>&, std::vector<std::vector<std::size_t>>&, const int);
