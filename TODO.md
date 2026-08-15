# TODO

  * normals handling
  * is it possible to make CGAL `Mesh_3` without polygon soup step?
  * `remeshIsotropic()` does not work well -> because of `soup_to_mesh()`?
  * `mesh.collect_garbage();` required in more places? (remeshing)
  * `storageMode("double")` required in more places?
  * `readFile_cpp()` STL filename instead of infile?

  # Wishlist
  
  * silent option to remove output during processing
  * fewer checks than in `soup_to_mesh()` when not required (via `SurfMesh_cpp` -> `makeSurfMesh()`)
      * use `vf_to_emesh()`, `csoup_to_mesh()`?
      * conversion from `mesh3d`
  * `Sqrt3Subdivision()`, `CatmullClark()`, `DooSabin_subdivision()`
  * C++: pass by reference where possible
  * configure script: gmp, mpfr necessary?
