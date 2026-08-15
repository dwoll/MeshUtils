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
Rcpp::List reconstructAFS_cpp(const Rcpp::NumericMatrix pts,
                              const unsigned nNeighs,
                              const bool clean) {
  std::vector<Point3> points = matrix_to_points3<Point3>(pts);
  if(nNeighs >= 2) {
    CGAL::jet_smooth_point_set<CGAL::Sequential_tag>(points, nNeighs);
  }

  AFS_triangulation3 dt(points.begin(), points.end());
  AFS_reconstruction reconstruction(dt);
  reconstruction.run();
  const AFS_Tds2& tds = reconstruction.triangulation_data_structure_2();
  std::vector<Point3> vertices;
  vertices.reserve(pts.ncol());
  size_t counter = 0;
  for(AFS_Tds2::Face_iterator fit = tds.faces_begin();
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
  std::vector<std::vector<size_t>> triangles;
  triangles.reserve(counter);
  for(size_t i = 0; i < counter; i++) {
    const size_t k = 3*i;
    const std::vector<size_t> triangle = {k, k+1, k+2};
    triangles.emplace_back(triangle);
  }

  PMP::merge_duplicate_points_in_polygon_soup(vertices, triangles);
  Mesh3 mesh = soup_to_mesh<Mesh3, Point3>(vertices, triangles, clean, false);
  Rcpp::List rmeshOut = getRmesh(mesh);
  return rmeshOut;
}
