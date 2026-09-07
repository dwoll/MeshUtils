# TODO

  * use `make_surf_mesh_valid()` instead of `make_surf_mesh()` wherever reasonable
  * `readMeshFile()` and `makeMesh()` in one step - via `filename` option
  * `fill_boundary_hole()`
      * pass more parameters (small holes)
  * `remeshIsotropic()`
      * does not work well for pentagrammic prism - why? (`Rvcg::vcgIsotropicRemeshing()` works)

# Wishlist

  * messages
      * verbose/silent option to remove output during processing except for warnings()
      * `rmessage()` vs. `Rcpp::warning()`
  * triangulated surface mesh simplification
      * `edge_collapse()`
  * remeshing
      * https://www.cgal.org/2025/05/22/Surface_remeshing/
      * https://doc.cgal.org/latest/PMP_Remeshing/index.html#Chapter_PMPRemeshing
      * `approximated_centroidal_Voronoi_diagram_remeshing()`
      * `surface_Delaunay_remeshing()`
      * `PMP::remesh_planar_patches()`
      * `PMP::remesh_almost_planar_patches()`
  * smoothing
      * `PMP::angle_and_area_smoothing()`
      * `PMP::tangential_relaxation()`
  * bounding meshes
      * approximate bounding ellipsoid
      * bounding spheres

# CRAN

  * `cran-comments.md` file
