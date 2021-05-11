# seb - Smallest Enclosing Ball 
`seb` is an R package which provides an interface for efficiently computing the [smallest enclosing ball](https://en.wikipedia.org/wiki/Smallest-circle_problem) (*seb*) of a point set. The underlying code is ported from C++ code in the [miniball](https://github.com/hbf/miniball/) codebase which has been shown to be very efficient in low dimensions and practically efficient in dimensions upwards of 10,000. The technique for computing the *seb* is described in the following paper: 

>  Fischer, Kaspar, Bernd Gärtner, and Martin Kutz. "Fast smallest-enclosing-ball computation in high dimensions." *European Symposium on Algorithms*. Springer, Berlin, Heidelberg, 2003.

The returned *seb* is exact, not an approximation---note the smallest bounding sphere of a non-empty point set is guaranteed both to exist and to be unique.

The method for computing the *seb* used here is an iterative procedure which continuously shrinks an enclosing ball towards the circumcenter of a set of intermediate boundary points. The circumcenter is continuously updated by sequences of orthogonal projections onto an affine subspace spanned by the intermediate boundary points, similar to the simplex algorithm from linear programming. Coefficients associated with the subspace projection inform the algorithm of which points to ignore as the algorithm 'walks' along certain convex hulls of point subsets. Once the circumsphere of these intermediate boundary sets encapsulates the whole of the original point set, the algorithm terminates and the *seb* is returned. 

Some quick benchmarks using [microbenchmark](https://cran.r-project.org/web/packages/microbenchmark/) (100 runs each benchmark): 

```R
library(microbenchmark) 
setup_expr <- quote({ x <- replicate(n = d, expr = rnorm(1500)) })
avg_ms <- sapply(2:15, function(d) {
  bench <- microbenchmark({ seb::seb(x) }, times = 100L, setup = eval(setup_expr))
  mean(bench$time) * 1e-6 ## milliseconds
})
cat(sprintf("dim=%d, n=1500 took %.2g milliseconds\n", 2:15, avg_ms))
# dim=2, n=1500 took 0.45 milliseconds
# dim=3, n=1500 took 0.49 milliseconds
# dim=4, n=1500 took 0.49 milliseconds
# dim=5, n=1500 took 0.6 milliseconds
# dim=6, n=1500 took 0.99 milliseconds
# dim=7, n=1500 took 0.68 milliseconds
# dim=8, n=1500 took 0.7 milliseconds
# dim=9, n=1500 took 0.83 milliseconds
# dim=10, n=1500 took 0.87 milliseconds
# dim=11, n=1500 took 0.97 milliseconds
# dim=12, n=1500 took 0.99 milliseconds
# dim=13, n=1500 took 1.1 milliseconds
# dim=14, n=1500 took 1.2 milliseconds
# dim=15, n=1500 took 1.2 milliseconds
```

