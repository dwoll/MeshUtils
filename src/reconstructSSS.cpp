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
// [[Rcpp::export]]
Rcpp::List reconstructSSS_cpp(
  const Rcpp::NumericMatrix pts,
  const unsigned int scaleIterations,
  const unsigned int nNeighbors,
  const unsigned int nSamples,
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
  return getRmesh(mesh, false);
}
