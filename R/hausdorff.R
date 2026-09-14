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

#' @title Hausdorff distance between two meshes
#' @description Hausdorff distance between two meshes. Either the
#'   approximate distance, or the distance estimate with a given error bound.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param mesh2 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param symmetric Boolean. Whether to consider the symmetric Hausdorff
#'   distance.
#' @param n \code{integer}. Number of points sampled for the approximate Hausdorff distance.
#'   If missing and \code{errorBound} is missing as well, the number of vertices is used.
#' @param errorBound A positive number. Upper bound on the error of the
#'   estimate. If missing, the approximate distance is returned.
#' @returns A number. For the apprixmate distance, the algorithm uses
#'   simulation and thus the result can vary.
#' @details See \url{https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#PMPDistance} for details.
#'   The approximate distance uses random uniform vertex sampling, edge sampling, face sampling. See
#'   \url{https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__distance__grp.html}
#'   for details.
#' @seealso See \code{\link[MeshUtils]{getHausdorffQuantile}} for the quantile Hausdorff distance.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' ## approximate symmetric Hausdorff distance
#' getHausdorff(dataHeart1, dataHeart2, symmetric=TRUE)
#'
#' ## estimate with error bound
#' getHausdorff(dataHeart1, dataHeart2, symmetric=TRUE,
#'                      errorBound=0.001)
#'
#' @export
getHausdorff <- function(mesh1, mesh2, symmetric = TRUE, n, errorBound) {
  stopifnot(inherits(mesh1, "CGALmesh"))
  stopifnot(inherits(mesh2, "CGALmesh"))
  stopifnot(isBoolean(symmetric))
  meshCPP1 <- fromR(mesh1)
  meshCPP2 <- fromR(mesh2)
  if(!missing(errorBound)) {
    stopifnot(isPositiveNumber(errorBound))
    getHausdorffEst_cpp(meshCPP1, meshCPP2, symmetric, errorBound)
  } else {
    if(missing(n)) {
      n <- 0L
    } else {
      stopifnot(isStrictPositiveInteger(n))
    }
    getHausdorffApprox_cpp(meshCPP1, meshCPP2, symmetric, n)
  }
}

#' @title Quantile Hausdorff distance between two meshes
#' @description Approximate the Hausdorff distance between two meshes using a
#'   quantile (by default the 95th percentile, i.e. "HD95") of the
#'   distances from a random sample of points on one mesh to the other mesh,
#'   instead of their maximum. This is less sensitive to small, isolated
#'   outlying regions than the full Hausdorff distance.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param symmetric Boolean. Whether to pool the sampled distances from
#'   \code{mesh1} to \code{mesh2} with those from \code{mesh2} to \code{mesh1}
#'   before taking the quantile.
#' @param p A number in \eqn{[0, 1]}. The quantile probability, defaults
#'   to \code{0.95}.
#' @param n \code{integer}. Number of points sampled on each mesh.
#'   If missing, the number of vertices is used.
#' @returns A number: the requested quantile of the sampled point-to-mesh
#'   distances. The algorithm uses random sampling, so the result can vary.
#' @details See \url{https://metrics-reloaded.dkfz.de/metric-library/xhd} for details.
#'   Note: This function under the hood simply calls
#'   \code{\link[MeshUtils]{getMetro}} and extracts the quantile Hausdorff distance.
#'   See there for details on the sampling process.
#' @seealso See \code{\link[MeshUtils]{getHausdorff}} for the (approximate) Hausdorff distance.
#'   See \code{\link[MeshUtils]{getMetro}} for other surface metrics.
#' @author Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' ## symmetric 95\% Hausdorff distance ("HD95")
#' getHausdorffQuantile(dataHeart1, dataHeart2, n=1000L, p=0.95)
#'
#' @export
getHausdorffQuantile <- function(mesh1, mesh2, symmetric = TRUE, p = 0.95, n) {
  metroL <- getMetro(mesh1, mesh2, symmetric, p, n)
  metroL[["HDq"]]
}

#' @title Several distance metrics between two meshes
#' @description Quantile Hausdorff distance (HD), average symmetric surface distance (ASSD),
#'   and root mean squared error (RMSE) for the surface distance between two meshes.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param mesh2 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param symmetric Boolean. Whether to pool the sampled distances from
#'   \code{mesh1} to \code{mesh2} with those from \code{mesh2} to \code{mesh1}
#'   before taking the quantile.
#' @param p A number in \eqn{[0, 1]}. The quantile probability, defaults
#'   to \code{0.95}.
#' @param n \code{integer}. Number of points sampled on each mesh.
#'   If missing, the number of vertices is used.
#' @returns A number: the requested quantile of the sampled point-to-mesh
#'   distances. The algorithm uses random sampling, so the result can vary.
#' @details See \code{\link[MeshUtils]{getHausdorffQuantile}} for details on the quantile
#'   Hausdorff distance, \url{https://metrics-reloaded.dkfz.de/metric-library/assd} for ASSD.
#'   Inspired by \url{http://vcglib.net/metro.html}.
#'   Uses random uniform vertex sampling, edge sampling, face sampling. See
#'   \url{https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__distance__grp.html}
#'   for details.
#' @seealso See \code{\link[MeshUtils]{getHausdorffQuantile}} for the quantile Hausdorff distance.
#' @author Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' getMetro(dataHeart1, dataHeart2, n=1000L, p=0.95)
#'
#' @export
getMetro <- function(mesh1, mesh2, symmetric = TRUE, p = 0.95, n) {
  stopifnot(inherits(mesh1, "CGALmesh"))
  stopifnot(inherits(mesh2, "CGALmesh"))
  stopifnot(isBoolean(symmetric))
  if(missing(n)) {
    n <- 0L
  } else {
    stopifnot(isStrictPositiveInteger(n))
  }
  stopifnot(is.numeric(p), length(p) == 1L, p > 0, p < 1)
  meshCPP1 <- fromR(mesh1)
  meshCPP2 <- fromR(mesh2)
  getMetro_cpp(meshCPP1, meshCPP2, symmetric, p, as.integer(n))
}
