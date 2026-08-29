# MeshUtils: Utility functions for 3D meshes based on CGAL

`MeshUtils` ís an R package that supports basic processing of 3D meshes using as backend the 'C++' library [`CGAL`](https://www.cgal.org/) via R package [`RcppCGAL`](https://cran.r-project.org/package=RcppCGAL). Features:

  * Read in mesh files with common formats (STL, PLY, OBJ, OFF)
  * Conversion to / from class `mesh3d` from package [`rgl`](https://cran.r-project.org/package=rgl), also compatible with package [`Rvcg`](https://cran.r-project.org/package=Rvcg)
  * Mesh repair
      * Filling holes
      * Removing self intersections
  * Boolean mesh operations: union, difference, intersection
  * Isotropic remeshing
  * Surface reconstruction
      * AFS
      * SSS
      * Poisson
      * Alpha wrapping
  * Bounding box
      * Axis-parallel bounding box
      * Optimal (oriented) bounding box
  * Hausdorff distance between two meshes
  * Centroid
  * Area
  * Volume
  * Distance from points to mesh

For an application, see package [MeshAgreement](https://github.com/dwoll/MeshAgreement) that calculates various distance and similarity metrics for two given 3D meshes.

## CAVE

Note that this package requires CGAL headers version 6.2. As of August 2026, `RcppCGAL` provides CGAL version 6.1, but it is possible to update as explained in the [`RcppCGAL::set_cgal()` documentation](https://cloud.r-project.org/web/packages/RcppCGAL/refman/RcppCGAL.html#set_cgal).

## Development

This package includes code adapted from packages [`Boov`](https://github.com/stla/Boov/), [`PolygonSoup`](https://github.com/stla/PolygonSoup/), and [`cgalMeshes`](https://github.com/stla/cgalMeshes/) developed and copyright by [Stéphane Laurent](https://laustep.github.io/stlahblog/).  Currently, only a subset of the functionality of these packages is provided in `MeshUtils`.

A fork / adaptation of packages [`Boov`](https://github.com/stla/Boov/), [`PolygonSoup`](https://github.com/stla/PolygonSoup/), and [`cgalMeshes`](https://github.com/stla/cgalMeshes/) was carried out as upstream changes to CGAL introduced incompatibilities, and the packages were archived from [CRAN](https://cran.r-project.org/).

## License

This package is provided under the GPL-3 license, but it uses the C++ library [`CGAL`](https://www.cgal.org/). To use CGAL for commercial purposes, you must obtain a license from the [GeometryFactory](https://geometryfactory.com).
