# TODO

  * `mesh.collect_garbage();` required in more places? (remeshing), `mesh.has_garbage()`
  * normals handling - get, calculate, and assign normals (`cgalMesh::get/compute/assignNormals()`)
  * `fillBoundaryHole()` as separate function
  * is it possible to make CGAL `Mesh_3` without polygon soup step?
  * `remeshIsotropic()` does not work well for pentagrammic prism -> because of `soup_to_mesh()`?
  * `readFile_cpp()` STL filename instead of `infile`?
  * configure script: gmp, mpfr necessary?

  # Wishlist

  * silent option to remove output during processing
  * fewer checks than in `soup_to_mesh()` when not required (via `SurfMesh_cpp` -> `makeSurfMesh()`)
      * use `vf_to_emesh()`, `csoup_to_mesh()`?
      * conversion from `mesh3d`
  * Remeshing
      * https://www.cgal.org/2025/05/22/Surface_remeshing/
      * https://doc.cgal.org/6.1/Polygon_mesh_processing/index.html#mesh3rem
      * Isotropic remeshing: `Uniform_sizing_field` vs. `Adaptive_sizing_field`
  * Smoothing with `PMP::angle_and_area_smoothing()`, `PMP::smooth_shape()`, `PMP::tangential_relaxation()`
  * Subdivision with `Sqrt3Subdivision()`, `CatmullClark()`, `DooSabin_subdivision()`
  * Triangulated Surface Mesh Simplification
  * Bounding meshes
      * Approximate Bounding Ellipsoid
      * Bounding Spheres
