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
// EPEC kernel only
// currently not used
template <typename MeshT, typename PointT>
bool is_small_hole(//hlfdg_descriptor h,
                   //const EMesh3 &mesh,
                   typename boost::graph_traits<MeshT>::halfedge_descriptor h,
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
    /*
    if((max_hole_diam > 0)      &&
       (max_num_hole_edges > 0) &&
       !is_small_hole<MeshT, PointT>(h, mesh, max_hole_diam, max_num_hole_edges)) {
        continue;
    }
    */

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
    CGAL::Conforming_constrained_Delaunay_triangulation_3<KernelT> ccdt;
    ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);
    return success;
}

template bool removeSelfIntSoup<K, Point3>(
    std::vector<Point3>&, std::vector<std::vector<std::size_t>>&, const int);

template bool removeSelfIntSoup<EK, EPoint3>(
    std::vector<EPoint3>&, std::vector<std::vector<std::size_t>>&, const int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// mesh passed by const reference, then uses a polygon soup
// as container as the output will most likely be non-manifold
template <typename KernelT, typename MeshT, typename PointT>
MeshT removeSelfIntMesh(const MeshT &mesh, const int method) {
  std::vector<PointT> points;
  std::vector<std::vector<std::size_t>> polygons;
  PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
  const bool success = removeSelfIntSoup<KernelT, PointT>(points, polygons, method);
  MeshT mesh_out;
  if(PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
  } else {
    Message("Polygon soup not a polygon mesh after removing intersections. Nothing done.");
    return mesh;
  }
  if(PMP::does_self_intersect(mesh_out)) {
    Message("Mesh self-intersections could not be removed.");
  }
  return mesh_out;
}

template Mesh3  removeSelfIntMesh<K,  Mesh3,  Point3>(const Mesh3&,   const int);
template EMesh3 removeSelfIntMesh<EK, EMesh3, EPoint3>(const EMesh3&, const int);
