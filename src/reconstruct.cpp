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

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/jet_smooth_point_set.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>      Polyhedron;

typedef CGAL::Advancing_front_surface_reconstruction<>               AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3                          AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2           AFS_tds2;

typedef CGAL::Scale_space_surface_reconstruction_3<K>                SSS_reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother<K> SSS_smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher<K>    SSS_mesher;
typedef SSS_reconstruction::Facet_const_iterator                     SSS_facet_iterator;
typedef SSS_reconstruction::Point_const_iterator                     SSS_point_iterator;

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List reconstructAFS_cpp(const Rcpp::NumericMatrix pts,
                              const int nNeighbors,
                              const bool repairSoup) {
  std::vector<Point3> points = matrix_to_points3<Point3>(pts);
  if(nNeighbors >= 2) {
    CGAL::jet_smooth_point_set<CGAL::Sequential_tag>(points, nNeighbors);
  }

  AFS_triangulation3 dt(points.begin(), points.end());
  AFS_reconstruction reconstruction(dt);
  reconstruction.run();
  const AFS_tds2& tds = reconstruction.triangulation_data_structure_2();
  std::vector<Point3> vertices;
  vertices.reserve(pts.ncol());
  std::size_t counter = 0;
  for(AFS_tds2::Face_iterator fit = tds.faces_begin();
      fit != tds.faces_end(); ++fit) {
    if(reconstruction.has_on_surface(fit)) {
      counter++;
      AFS_triangulation3::Facet f = fit->facet();
      AFS_triangulation3::Cell_handle ch = f.first;
      int ci = f.second;
      Point3 vs[3];
      for(int i = 0, j = 0; i < 4; i++) {
        if(ci != i) {
          vs[j] = ch->vertex(i)->point();
          j++;
        }
      }
      for(int k = 0; k < 3; k++) {
        const Point3 p = vs[k];
        // const EPoint3 v = EPoint3(p.x(), p.y(), p.z());
        vertices.push_back(p);
      }
    }
  }
  std::vector<std::vector<std::size_t>> triangles;
  triangles.reserve(counter);
  for(std::size_t i = 0; i < counter; i++) {
    const std::size_t k = 3*i;
    const std::vector<std::size_t> triangle = { (k), (k+1), (k+2) };
    triangles.emplace_back(triangle);
  }

  // repair_polygon_soup() is called in soup_to_mesh()
  // PMP::merge_duplicate_points_in_polygon_soup(vertices, triangles);
  Mesh3 mesh = soup_to_mesh<K, Mesh3, Point3>(
      vertices,
      triangles,
      false,      // triangulate
      repairSoup, // repair_soup
      false,      // remove_intersections
      1,          // remove_method
      false,      // fill_holes
      false,      // fair hole
      0);         // max_num_holes
  return getRmesh<K, Mesh3, Point3, Vector3>(mesh, false);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List reconstructPoisson_cpp(const Rcpp::NumericMatrix pts,
                                  const Rcpp::NumericMatrix normals,
                                  const double spacing,
                                  const double smAngle,
                                  const double smRadius,
                                  const double smDistance) {
  const std::size_t nPts = pts.ncol();
  std::vector<P3V3> points(nPts);   // normals
  for(std::size_t i = 0; i < nPts; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    const Rcpp::NumericVector nrml_i = normals(Rcpp::_, i);
    points[i] =
      std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)),
                     Vector3(nrml_i(0), nrml_i(1), nrml_i(2)));
  }

  double spacingUse;
  if(spacing == -1.0) {
    spacingUse = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
      points, 6, /* knn = 1 ring */
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3V3>()));
  } else {
      spacingUse = spacing;
  }

  Polyhedron poly;
  const bool success = CGAL::poisson_surface_reconstruction_delaunay(
    points.begin(), points.end(),
    CGAL::First_of_pair_property_map<P3V3>(),
    CGAL::Second_of_pair_property_map<P3V3>(),
    poly,
    spacingUse, smAngle, smRadius, smDistance);

  if(!success) {
    Rcpp::stop("Poisson surface reconstruction has failed.");
  }

  std::size_t id = 1;
  for(Polyhedron::Vertex_iterator vit = poly.vertices_begin();
      vit != poly.vertices_end(); ++vit) {
    vit->id() = id;
    id++;
  }

  bool isTriangle = PMP::triangulate_faces(poly);
  if(!isTriangle) {
    Rcpp::stop("Could not triangluate polyhedron.");
  }

  Mesh3 mesh;
  CGAL::copy_face_graph(poly, mesh);
  return getRmesh<K, Mesh3, Point3, Vector3>(mesh, false);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List reconstructSSS_cpp(
  const Rcpp::NumericMatrix pts,
  const int scaleIterations,
  const int nNeighbors,
  const int nSamples,
  const bool separateShells,
  const bool forceManifold,
  const double borderAngle,
  const bool repairSoup) {
  std::vector<Point3> points = matrix_to_points3<Point3>(pts);
  SSS_reconstruction SSSR(points.begin(), points.end());
  SSS_smoother smoother(nNeighbors, nSamples);
  SSSR.increase_scale(scaleIterations, smoother);
  SSS_mesher mesher(
    smoother.squared_radius(), separateShells, forceManifold, borderAngle
  );
  SSSR.reconstruct_surface(mesher);
  SSS_reconstruction::Point_range smoothed(SSSR.points());
  SSS_reconstruction::Facet_range polygons(SSSR.facets());

  Mesh3 mesh;
  if(repairSoup) {
    PMP::repair_polygon_soup(smoothed, polygons);
  }
  PMP::orient_polygon_soup(smoothed, polygons);
  PMP::polygon_soup_to_polygon_mesh(smoothed, polygons, mesh);
  if(CGAL::is_triangle_mesh(mesh)) {
      if(!PMP::is_outward_oriented(mesh)) {
        PMP::reverse_face_orientations(mesh);
      }

      if(!PMP::does_bound_a_volume(mesh)) {
          PMP::orient_to_bound_a_volume(mesh);
      }
  }
  return getRmesh<K, Mesh3, Point3, Vector3>(mesh, false);
}
