## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/Boov/
## https://github.com/stla/PolygonSoup/
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title Isotropic remeshing
#' @description Isotropic remeshing.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param method \code{character}: Either \code{"uniform"} for uniform sizing field
#'    or \code{"adaptive"} for adapative sizing field.
#' @param targetEdgeLen Positive number for \code{method="uniform"}:
#'   The target edge length of the remeshed mesh.
#' @param tol Positive number for \code{method="adaptive"}:
#'   Error tolerance. See details
#' @param edgeMin Positive number for \code{method="adaptive"}:
#'   Minimum edge length. See details.
#' @param edgeMax Positive number for \code{method="adaptive"}:
#'   Maximum edge length. See details.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param nRelaxSteps Positive \code{integer}: Number of relaxation steps.
#' @param out \code{character}: Indicate output mesh format.
#' @return A \code{CGALmesh} object or a \code{\link[rgl]{mesh3d}} object from package \strong{rgl},
#'   depending on option \code{out}.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/classCGAL_1_1Polygon__mesh__processing_1_1Adaptive__sizing__field.html},
#'   \url{https://doc.cgal.org/latest/PMP_Remeshing/group__PMP__local__remeshing__grp.html#ga412f696ec3009074bf957f1bba638248}.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh         <- makeMesh(mesh=dataTruncIcosahedron, triangulate=TRUE)
#' mesh_rgl     <- toRGL(mesh)
#' mesh_rem_rgl <- remeshIsotropic(mesh, targetEdgeLen=1, out="rgl")
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_rem_rgl)
#'
#' @export
remeshIsotropic <- function(
        x,
        method = c("uniform", "adaptive"),
        tol = 0.001,
        edgeMin = 0.001,
        edgeMax,
        targetEdgeLen,
        nIter = 1,
        nRelaxSteps = 1,
        out = c("CGALmesh", "rgl")) {
    method <- match.arg(method)
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
    }
    out <- match.arg(out)
    stopifnot(isStrictPositiveInteger(nIter))
    stopifnot(isStrictPositiveInteger(nRelaxSteps))
    meshCPP <- fromR(x)
    if(method == "uniform") {
      stopifnot(isPositiveNumber(targetEdgeLen))
      meshRe <- remeshIsotropicUniform_cpp(
          meshCPP,
          targetEdgeLen,
          as.integer(nIter),
          as.integer(nRelaxSteps))
    } else if(method == "adaptive") {
      stopifnot(isPositiveNumber(tol))
      stopifnot(isPositiveNumber(edgeMin))
      stopifnot(isPositiveNumber(edgeMax))

      meshRe <- remeshIsotropicAdapt_cpp(
          meshCPP,
          tol,
          edgeMin,
          edgeMax,
          as.integer(nIter),
          as.integer(nRelaxSteps))
    }

    meshReWN <- vcgUpdateNormals(
        tmesh3d(meshRe[["vertices"]],
                meshRe[["faces"]],
                normals=NULL))
    meshOut <- if(out == "rgl") {
        meshReWN
    } else {
        makeMesh(mesh=meshReWN)
    }
    meshOut
}
