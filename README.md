# MeshUtils: Utility functions for meshes based on CGAL

`MeshUtils` supports basic processing of 3D meshes as well as Boolean operations (union, difference and intersection) using as backend the 'C++' library [`CGAL`](https://www.cgal.org/) via R package [`RcppCGAL`](https://cran.r-project.org/package=RcppCGAL). Conversion to / from class `mesh3d` from package [`rgl`](https://cran.r-project.org/package=rgl) is possible.

## Development

This package includes code adapted from packages [`Boov`](https://github.com/stla/Boov/), [`PolygonSoup`](https://github.com/stla/PolygonSoup/), and [`cgalMeshes`](https://github.com/stla/cgalMeshes/) developed and copyright by [Stéphane Laurent](https://laustep.github.io/stlahblog/).  Currently, only a subset of the functionality of these packages is provided in `MeshUtils`.

A fork / adaptation of packages [`Boov`](https://github.com/stla/Boov/), [`PolygonSoup`](https://github.com/stla/PolygonSoup/), and [`cgalMeshes`](https://github.com/stla/cgalMeshes/) was carried out as upstream changes to CGAL introduced incompatibilities, and the packages were archived from [CRAN](https://cloud.r-project.org/web/packages/index.html).

## License

This package is provided under the GPL-3 license, but it uses the C++ library [`CGAL`](https://www.cgal.org/). To use CGAL for commercial purposes, you must obtain a license from the [GeometryFactory](https://geometryfactory.com).
