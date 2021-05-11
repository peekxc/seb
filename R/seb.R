#' @title Smallest Enclosing Ball
#' @description Computes the Smallest Enclosing Ball (SEB) of a point set. 
#' @param x matrix of points. 
#' @param validate whether to validate the resulting ball. Defaults to TRUE. 
#' @details This function computes the smallest enclosing ball (the \emph{seb}) of a point set \code{x} of 
#' arbitrary dimension. The smallest bounding sphere of a non-empty point set is guaranteed both 
#' to exist and to be unique. If you're unfamiliar with this problem, see [1]. For this implementation, the 
#' returned seb is exact, not an approximation. The code has been shown to work with dimensions upwards of 10,000 [2]. \cr 
#' \cr
#' The method for computing the \emph{seb} is an iterative procedure which continuously shrinks 
#' an enclosing ball towards the circumcenter of a set of intermediate boundary points. The circumcenter 
#' is continuously updated by sequences of orthogonal projections onto an affine subspace spanned by the 
#' intermediate boundary points, similar to the simplex algorithm from linear programming. Coefficients 
#' associated with the subspace projection inform the algorithm of which points to ignore as the algorithm
#' 'walks' along certain convex hulls of point subsets. Once the circumsphere of these intermediate boundary sets 
#' encapsulates the whole of the original point set, the algorithm terminates and the seb(x) is returned. See
#' [2] for more details. \cr
#' \cr
#' The orthogonal projection step which allows updating the circumcenter of 
#' the intermediate boundary points involves performing rank-1 updates to an 
#' underlying QR decomposition. 
#' 
#' \itemize{
#'    \item{\strong{min_lambda}}{ the minimum lambda value of the convex combination of boundary points with respect to 
#' the returned center point. If positive, this ensures the center point is in the convex hull of the 
#' boundary points used to compute the seb.}
#'    \item{\strong{max_overlength}}{ the ratio of the maximum distance of any point in the original set 
#' over the radius of the reported seb, or 0 if all points reported within a distance of 'radius' 
#' to the seb. Expected to be 0 or close to 0 in cases where points are in general position, but may be
#' positive if many points lie on the boundary of the seb.}
#'    \item{\strong{min_underlength}}{ the absolute-valued ratio between the distance of a point in the boundary set 
#' to the center point of the reported seb and the reported radius. Expected to be 0 or close to 0.} 
#'    \item{\strong{qr_inconsistency}}{ indicates the degree to which the QR factorization was 
#' inconconsistent. It is computed by measuring the maximum discrepancy found during the 
#' orthogonal projection phase of the algorithm.}
#' }
#' This implementation is a direct port of the C++ code given in [3].
#' @author The R interface was written by Matt Piekenbrock. The original C++ code was written by Martin Kutz, Kaspar Fischer, and Bernd Gärtner.
#' @return A list with components \emph{center} and \emph{radius} describing the seb of \code{x}. If \code{validate} was set to TRUE (the default), 
#' information about the iterative procedure is recorded in \emph{info} attribute of the list. 
#' @references 
#' \enumerate{
#'   \item https://en.wikipedia.org/wiki/Bounding_sphere
#'   \item Fischer, Kaspar, Bernd Gärtner, and Martin Kutz. "Fast smallest-enclosing-ball computation in high dimensions." European Symposium on Algorithms. Springer, Berlin, Heidelberg, 2003. 
#'   \item https://github.com/hbf/miniball
#' }
#' 
#' @examples
#' ## Random point set 
#' x <- replicate(2, rnorm(40))
#' 
#' ## Compute enclosing ball 
#' sebx <- seb::seb(x)
#' 
#' ## Visualize the ball
#' range <- apply(x, 2, function(dim){ range(dim) + diff(range(dim))*0.30*c(-1,1) })
#' plot(x, xlim = range[,1], ylim = range[,1], asp = 1)
#' points(sebx$center[1], sebx$center[2], pch = 20, col = "red")
#' 
#' \dontrun{ 
#' ## This uses plotrix to draw the ball
#' plotrix::draw.circle(x = sebx$center[1], y = sebx$center[2], radius = sebx$radius)
#' }
#' 
#' @useDynLib seb, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
seb <- function(x, validate = TRUE){
	stopifnot(is.matrix(x), is.logical(validate))
	results <- seb_rcpp(x, validate)
	if (validate){
		attr(results, "info") <- results$info
		results["info"] <- NULL
	}
	return(results)
}