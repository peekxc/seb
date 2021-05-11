#include <Rcpp.h>
using namespace Rcpp;

#include "seb.h"
typedef SEB_NAMESPACE::Point< double > Point;
typedef std::vector< Point > PointVector;
typedef SEB_NAMESPACE::Smallest_enclosing_ball< double > Miniball;

// [[Rcpp::export]]
List seb_rcpp(NumericMatrix x, bool validate) {
	
	const size_t n = x.nrow();
	const size_t d = x.ncol(); 
	
  // Construct n random points in dimension d
  PointVector S;
  for (size_t i = 0; i < n; ++i) {
  	NumericVector r = x.row(i);
    S.push_back(Point(d,r.begin()));
  }

  // Compute the miniball by inserting each value
  Miniball mb(d, S);
 
  // Output
  NumericVector center = NumericVector(mb.center_begin(), mb.center_end()); 
  List output = List::create(
  	_["radius"] = static_cast<double>(mb.radius()),
  	_["center"] = center
  );
  
  if (validate){
  	std::array< double, 4 > res_info = mb.verify();
  	NumericVector E = NumericVector(res_info.begin(), res_info.end());
  	output["info"] = E;
  }
  
  return(output);
}

/*** R
x <- replicate(2, rnorm(40))
sebx <- seb::seb(x)
range <- apply(x, 2, function(dim){ range(dim) + diff(range(dim))*cc*c(-1,1) })
cc <- 0.30
plot(x, xlim = range[,1], ylim = range[,1], asp = 1)
points(sebx$center[1], sebx$center[2], pch = 20, col = "red")
plotrix::draw.circle(x = sebx$center[1], y = sebx$center[2], radius = sebx$radius)
*/
