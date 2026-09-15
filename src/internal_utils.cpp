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
#include <numeric>

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

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::parameters::use_random_uniform_sampling(true)     // true
// PMP::parameters::use_grid_sampling(false)              // false
// PMP::parameters::use_monte_carlo_sampling(false)       // false
//
// PMP::parameters::do_sample_vertices(true)              // true
// PMP::parameters::do_sample_edges(true)                 // true
// PMP::parameters::do_sample_faces(true)                 // true
//
// PMP::parameters::grid_spacing(n)                       // double, for grid sampling
//
// PMP::parameters::number_of_points_on_faces(n)          // *unsigned int, for random sampling
// PMP::parameters::number_of_points_on_edges(n)          //  unsigned int, for random sampling
//
// PMP::parameters::number_of_points_per_distance_unit(n) // double, for random sampling and Monte Carlo sampling
// PMP::parameters::number_of_points_per_area_unit(n)     // double, for random sampling and Monte Carlo sampling
// PMP::parameters::number_of_points_per_edge(n)          // unsigned int, for Monte-Carlo sampling
// PMP::parameters::number_of_points_per_face(n)          // unsigned int, for Monte-Carlo sampling
// distances from a random sample of points on `mesh_source` to `mesh_target`
template <typename KernelT, typename MeshT, typename PointT>
std::vector<double> sampled_distances_to_mesh(
  const MeshT& mesh_source,
  const MeshT& mesh_target,
  const sample_opts &opt) {
  typedef CGAL::AABB_face_graph_triangle_primitive<MeshT> Primitive;
  typedef CGAL::AABB_traits_3<KernelT, Primitive> Tree_Traits;
  typedef CGAL::AABB_tree<Tree_Traits> Tree;
  std::vector<PointT> pts;
  if(opt.ptsOnFaces > 0) {
    PMP::sample_triangle_mesh(
      mesh_source, std::back_inserter(pts),
      PMP::parameters::number_of_points_on_faces(opt.ptsOnFaces));
  } else {
    PMP::sample_triangle_mesh(
      mesh_source, std::back_inserter(pts));
  }
  Tree tree(faces(mesh_target).first, faces(mesh_target).second, mesh_target);
  std::vector<double> dsts;
  dsts.reserve(pts.size());
  for(const PointT& pt : pts) {
    double dsq = CGAL::to_double<typename KernelT::FT>(tree.squared_distance(pt));
    dsts.push_back(std::sqrt(dsq));
  }
  return dsts;
}

template std::vector<double> sampled_distances_to_mesh<K, Mesh3, Point3>(
    const Mesh3&,
    const Mesh3&,
    const sample_opts&);

template std::vector<double> sampled_distances_to_mesh<EK, EMesh3, EPoint3>(
    const EMesh3&,
    const EMesh3&,
    const sample_opts&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// Hausdorff distance quantile
// average symmetric surface distance
// root mean squared error
// https://github.com/zarquon42b/Rvcg/blob/master/src/metroSampling.h
template <typename KernelT, typename MeshT, typename PointT>
std::tuple<double, double, double> get_metro(
  const MeshT& mesh1,
  const MeshT& mesh2,
  const bool symmetric,
  const double p,
  const sample_opts &opt) {
  // distances of sampled points from mesh1 to mesh2
  std::vector<double> dsts12 = sampled_distances_to_mesh<KernelT, MeshT, PointT>(
    mesh1, mesh2, opt);
  // distances of sampled points from mesh2 to mesh1
  std::vector<double> dsts21 = sampled_distances_to_mesh<KernelT, MeshT, PointT>(
    mesh2, mesh1, opt);
  // number of samples may differ between meshes
  const std::size_t nDst = dsts12.size() + dsts21.size();
  // Hausdorff distance quantile
  const std::optional<double> dst12_q = get_quantile(dsts12, p);
  const std::optional<double> dst21_q = get_quantile(dsts21, p);
  double HDq;
  if(!dst12_q.has_value() || !dst21_q.has_value()) {
    HDq = std::nan("0");
  } else {
    if(symmetric) {
        HDq = std::max(dst12_q.value(), dst21_q.value());
    } else {
        HDq = dst12_q.value();
    }
  }

  // average symmetric surface distance
  const double sum_dsts12 = std::reduce(dsts12.begin(), dsts12.end()); // could be auto sum12
  const double sum_dsts21 = std::reduce(dsts21.begin(), dsts21.end()); // could be auto sum21
  const double assd = (sum_dsts12 + sum_dsts21) / static_cast<double>(nDst);
  // root mean squared error
  const double ssq_dsts12 = std::inner_product(dsts12.begin(), dsts12.end(), dsts12.begin(), 0.0); // sum of squared distances
  const double ssq_dsts21 = std::inner_product(dsts21.begin(), dsts21.end(), dsts21.begin(), 0.0); // sum of squared distances
  const double msq_dst    = (ssq_dsts12 + ssq_dsts21) / static_cast<double>(nDst);    // mean squared distance
  const double rmse       = std::sqrt(msq_dst);
  return std::tuple<double, double, double>(HDq, assd, rmse);
}

template std::tuple<double, double, double> get_metro<K, Mesh3, Point3>(
  const Mesh3&,
  const Mesh3&,
  const bool,
  const double,
  const sample_opts&);

template std::tuple<double, double, double> get_metro<EK, EMesh3, EPoint3>(
  const EMesh3&,
  const EMesh3&,
  const bool,
  const double,
  const sample_opts&);
