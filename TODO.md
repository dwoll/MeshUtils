# TODO

  * `fillBoundaryHole()`
      * template
      * pass more parameters
  * `remeshIsotropic()` does not work well for pentagrammic prism
   * normals handling
      * get
      * calculate
      * assign (`cgalMesh::get/compute/assignNormals()`)
  * `readFile_cpp()` STL filename instead of `infile`?
  * configure script: gmp, mpfr necessary?

  # Wishlist

  * add CGAL documentation links to help pages for details
  * silent option to remove output during processing
  * fewer checks than in `soup_to_mesh()` when not required
      * use `vf_to_mesh()`, `csoup_to_mesh()`?
      * all conversions from existing `CGALmesh`
  * remeshing
      * https://www.cgal.org/2025/05/22/Surface_remeshing/
      * https://doc.cgal.org/6.1/Polygon_mesh_processing/index.html#mesh3rem
  * smoothing
      * `cgalMeshes`
      * `PMP::angle_and_area_smoothing()`, `PMP::smooth_shape()`, `PMP::tangential_relaxation()`
  * triangulated surface mesh simplification
  * bounding meshes
      * approximate bounding ellipsoid
      * bounding spheres
