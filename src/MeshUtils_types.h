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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// #include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

// -------------------------------------------------------------------------- //
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel   EK;

typedef K::Point_3  Point3;
typedef EK::Point_3 EPoint3;

typedef K::Vector_3  Vector3;
typedef EK::Vector_3 EVector3;

typedef CGAL::Surface_mesh<Point3>  Mesh3;
typedef CGAL::Surface_mesh<EPoint3> EMesh3;
