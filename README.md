# MeshUtils: Utility functions for 3D meshes based on CGAL

`MeshUtils` ís an R package that supports basic processing of 3D meshes using as backend the 'C++' library [`CGAL`](https://www.cgal.org/) via R package [`RcppCGAL`](https://cran.r-project.org/package=RcppCGAL). Features:

  * Read in mesh files with common formats (STL, PLY, OBJ, OFF)
  * Conversion to / from class `mesh3d` from package [`rgl`](https://cran.r-project.org/package=rgl), also compatible with package [`Rvcg`](https://cran.r-project.org/package=Rvcg)
  * Mesh repair
      * Filling holes
      * Removing self intersections
  * Boolean mesh operations
      * Union
      * Dfference
      * Intersection
  * Isotropic remeshing
  * Smoothing
  * Subdivision
      * Catmull-Clark
      * Doo-Sabin
      * Sqrt3
      * Loop
  * Surface reconstruction
      * AFS
      * SSS
      * Poisson
      * Alpha wrapping
  * Bounding box
      * Axis-parallel bounding box
      * Optimal (oriented) bounding box
  * Hausdorff distance between two meshes
  * Distance from points to mesh
  * Centroid
  * Area
  * Volume
  * Vertex normals

For an application, see package [MeshAgreement](https://github.com/dwoll/MeshAgreement/) that calculates various distance and similarity metrics for two given 3D meshes.

## CAVE

Note that this package requires CGAL headers version 6.2. As of August 2026, `RcppCGAL` provides CGAL version 6.1, but it is possible to update as explained in the [`RcppCGAL::set_cgal()` documentation](https://ericdunipace.r-universe.dev/RcppCGAL/doc/manual.html#set_cgal).

See package [`Rmpfr`](https://cran.r-project.org/package=Rmpfr) for a note on how to install system requirements [MPFR](https://www.mpfr.org/) and [GMP](https://gmplib.org/).

## Implementation

This package includes code adapted from packages [`Boov`](https://github.com/stla/Boov/), [`PolygonSoup`](https://github.com/stla/PolygonSoup/), and [`cgalMeshes`](https://github.com/stla/cgalMeshes/) developed and copyright by [Stéphane Laurent](https://laustep.github.io/stlahblog/).  Currently, only a subset of the functionality of these packages is provided in `MeshUtils`.

A fork / adaptation of packages [`Boov`](https://github.com/stla/Boov/), [`PolygonSoup`](https://github.com/stla/PolygonSoup/), and [`cgalMeshes`](https://github.com/stla/cgalMeshes/) was carried out as upstream changes to CGAL introduced incompatibilities, and the packages were archived from [CRAN](https://cran.r-project.org/).

The design was chosen such that mesh data resides in R space. This means that for each mesh operation, data is first transferred to the C++ side (using `Rcpp`), converted to a CGAL surface mesh, subjected to CGAL functions, and then transferred back to R. The package does not maintain a pointer to a C++ data structure to keep the mesh data there - unlike packages such as [`terra`](https://cran.r-project.org/package=terra) or [`cgalMeshes`](https://github.com/stla/cgalMeshes/). This approach carries a performance penalty, but from an R perspective, it is more straightforward. In particular, there are no serialization issues (saving meshes). Furthermore, memory management is easier.

## License

This package is provided under the GPL-3 license, but it uses the C++ library [`CGAL`](https://www.cgal.org/). To use CGAL for commercial purposes, you must obtain a license from the [GeometryFactory](https://geometryfactory.com).
