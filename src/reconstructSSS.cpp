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

// [[Rcpp::export]]
Rcpp::List reconstructSSS_cpp( // Rcpp::XPtr<EMesh3>
  const Rcpp::NumericMatrix pts,
  const size_t scaleIterations,
  const unsigned nneighs,
  const unsigned nsamples,
  const bool separateShells,
  const bool forceManifold,
  const double borderAngle
) {

  std::vector<Point3> points = matrix_to_points3<Point3>(pts);

  SSS_reconstruction SSSR(points.begin(), points.end());
  SSS_smoother smoother(nneighs, nsamples);
  SSSR.increase_scale(scaleIterations, smoother);
  SSS_mesher mesher(
    smoother.squared_radius(), separateShells, forceManifold, borderAngle
  );
  SSSR.reconstruct_surface(mesher);

  Mesh3 mesh;
  for(
    SSS_point_iterator it = SSSR.points_begin(); it != SSSR.points_end(); ++it
  ) {
    mesh.add_vertex(*it);
  }
  for(
    SSS_facet_iterator it = SSSR.facets_begin(); it != SSSR.facets_end(); ++it
  ) {
    std::array<size_t, 3> face = *it;
    mesh.add_face(vxdescr(face[0]), vxdescr(face[1]), vxdescr(face[2]));
  }

  EMesh3 emesh;
  CGAL::copy_face_graph(mesh, emesh);

  // return Rcpp::XPtr<EMesh3>(new EMesh3(emesh), false);
  Rcpp::List out = getRmesh(emesh);
  return out;
}
