# TODO

  * drop `data.table` requirement -> replace `uniqueN()`
  * silent option to remove output during processing
  * only deal with triangle meshes -> simplify code that is more general
  * `approximate_Hausdorff_distance()` from CGAL?
      * `cgalMeshes::HausdorffApproximate()`
  * `readFile_cpp()` STL filename instead of infile?
  * write documentation (roxygen)
  * `reconstructPoisson_cpp()` could make neighbors for spacing a parameter, currently fixed at 6
      * `reconstructPoisson()` does not work as expected
  * ckecks in `SurfEMesh_cpp` (via `makeMesh()`) and `soup_to_mesh()` via (`makeSurfMesh()`): redundant
      * use `vf_to_mesh()`?
      * conversion from `mesh3d`
  * `cgalMeshes`: `isotropicRemeshing()`, `Sqrt3Subdivision()`, `CatmullClark()`, `DooSabin_subdivision()`
  * `storageMode("double")` required in more places?
  * C++: pass by reference where possible
  * configure script: gmp, mpfr necessary?
  * `fillBoundaryHoles()` - template to make work for EPIC, EPEC
  * either template more with respect to EPIC vs. EPEC construction kernels, or only use EPEC, or use EPIC except for union / intersection methods
