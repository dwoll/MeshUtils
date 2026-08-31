# TODO

  * normals
      * R to C++: Normals are lost
      * C++ to R: Vertex normals via `RSurfMesh1/2()` -> `getVxNormals()`
  * usage of `getRmesh()` vs. `RSurfMesh()`
      * Call `RsurfMesh()`: `makeMesh()`, `boolInt/Union/Diff()`, `optimalBoundingBox_cpp()`, `getRmesh()`
          * Computes normals via `getVxNormals()`
      * Everything else: call `getRmesh()`
          * Takes normals that are already present in CGAL mesh
          * Triangulates if required
          * Calls `RSurfMesh1/2()`
  * `getVxNormals()`
      * remove old maps first?
      * also assign on C++ side? (`cgalMesh::get/compute/assignNormals()`)
  * `smoothShape()`
      * check
  * `fillBoundaryHole()`
      * pass more parameters (small holes)
  * `remeshIsotropic()`
      * does not work well for pentagrammic prism - why? (`Rvcg::vcgIsotropicRemeshing()` works)
  * `readFile_cpp()` STL filename instead of `infile`?

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
      * `PMP::angle_and_area_smoothing()`
      * `PMP::tangential_relaxation()`
  * triangulated surface mesh simplification
  * bounding meshes
      * approximate bounding ellipsoid
      * bounding spheres
