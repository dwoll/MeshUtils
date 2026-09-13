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

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

#include <algorithm>
#include <cmath>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
void rmessage(std::string msg) {
  SEXP rmsg = Rcpp::wrap(msg);
  Rcpp::message(rmsg);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
bool is_triangle_soup(const std::vector<std::vector<std::size_t>>& polygons) {
    for (const auto& poly : polygons) {
        if (poly.size() != 3) {
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// distances from a random sample of points on `mesh` to `mesh_target`,
// one point at a time so that max_distance_to_triangle_mesh() (which only
// reports the distance of the single furthest point) yields a per-point value
template <typename KernelT, typename MeshT, typename PointT>
std::vector<double> sampled_distances_to_mesh(
  const MeshT& mesh_source, const MeshT& mesh_target, const unsigned int n) {
  typedef CGAL::AABB_face_graph_triangle_primitive<MeshT> Primitive;
  typedef CGAL::AABB_traits_3<KernelT, Primitive> Tree_Traits;
  typedef CGAL::AABB_tree<Tree_Traits> Tree;
  std::vector<PointT> pts;
  if(n > 0) {
    PMP::sample_triangle_mesh(
      mesh_source, std::back_inserter(pts),
      PMP::parameters::number_of_points_on_faces(n));
  } else {
    PMP::sample_triangle_mesh(
      mesh_source, std::back_inserter(pts));
  }
  Tree tree(faces(mesh_target).first, faces(mesh_target).second, mesh_target);
  std::vector<double> dists;
  dists.reserve(pts.size());
  for(const PointT& p : pts) {
    std::vector<PointT> pt = { p };
    dists.push_back(std::sqrt(tree.squared_distance(p)));
  }
  return dists;
}

template std::vector<double> sampled_distances_to_mesh<K, Mesh3, Point3>(
    const Mesh3&, const Mesh3&, const unsigned int);

template std::vector<double> sampled_distances_to_mesh<EK, EMesh3, EPoint3>(
    const EMesh3&, const EMesh3&, const unsigned int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// use std::nth_element() to avoid full sorting
std::optional<double> get_quantile(std::vector<double> &data, double p) {
    if(data.empty() || !std::isfinite(p) || (p <= 0.0) || (p >= 1.0)) {
        return std::nullopt;
    }
    std::vector<double>::iterator it_b = data.begin();
    std::vector<double>::iterator it_e = data.end();
    double idx = p * (data.size() - 1); // index, may be fractional
    const std::size_t pos_lower = static_cast<std::size_t>(std::floor(idx));
    const std::size_t pos_upper = static_cast<std::size_t>(std::ceil( idx));
    std::vector<double>::iterator it_lower = it_b; // to be the quantile (lower)
    std::advance(it_lower, pos_lower);
    std::nth_element(it_b, it_lower, it_e);
    double q_lower = *it_lower;
    if(pos_lower == pos_upper) {
        return q_lower;
    } else {
        std::vector<double>::iterator it_upper = it_b; // to be the quantile (upper)
        std::advance(it_upper, pos_upper);
        std::nth_element(it_b, it_upper, it_e);
        double q_upper = *it_upper;
        // linear interpolation between lower and upper value
        double weight = idx - pos_lower;
        return (1.0 - weight)*q_lower + weight*q_upper;
    }
}
