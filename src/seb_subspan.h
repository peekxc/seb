// Synopsis: Class representing the affine hull of a point set.
//
// Authors: Martin Kutz <kutz@math.fu-berlin.de>,
//          Kaspar Fischer <kf@iaeth.ch>

#ifndef SEB_SUBSPAN_H
#define SEB_SUBSPAN_H

#include <vector>
#include <cmath>
#include <numeric>
#include <ostream>

namespace SEB_NAMESPACE {

  template<typename Float>
  inline Float sqr(const Float x)
  {
    return x * x;
  }

  template<typename Float, class Pt, class PointAccessor>
  class Subspan
  // An instance of this class represents the affine hull of a
  // non-empty set M of affinely independent points.  The set M is not
  // represented explicity; when an instance of this class is
  // constructed, you pass a list S of points to it (which the
  // instance will never change and which is assumed to stay fixed for
  // the lifetime of this instance): The set M is then a subset of S,
  // and its members are identified by their (zero-based) indices in
  // S.  The following routines are provided to query and change the
  // set M:
  //
  // - int size() returns the size of the instance's set M, a number
  //   between 0 and dim+1.
  //   Complexity: O(1).
  //
  // - bool is_member(int global_index) returns true iff S[global_index]
  //   is a member of M.
  //   Complexity: O(1)
  //
  // - global_index(int local_index) returns the index (into S) of the
  //   local_index-th point in M.  The points in M are internally
  //   ordered (in an arbitrary way) and this order only changes when
  //   add() or remove() (see below) is called.
  //   Complexity: O(1)
  //
  // - void add_point(int global_index) adds the global_index-th point
  //   of S to the instance's set M.
  //   Precondition: !is_member(global_index)
  //   Complexity: O(dim^2).
  //
  // - void remove_point(int local_index) removes the local_index-th
  //   point in M.
  //   Precondition: 0<=local_index<=size() and size()>1
  //   Complexity: O(dim^2)
  //
  // - int any_member() returns the global index (into S) of an
  //   arbitrary element of M.
  //   Precondition: size()>0
  //   Postcondition: is_member(any_member())
  //
  // The following routines are provided to query the affine hull of M:
  //
  // - void shortest_vector_to_span(p,w): Computes the vector w
  //   directed from point p to v, where v is the point in aff(M) that
  //   lies nearest to p.  Returned is the squared length of w.
  //   Precondition: size()>0
  //   Complexity: O(dim^2)
  //
  // - void find_affine_coefficients(c,coeffs):
  //   Preconditions: c lies in the affine hull aff(M) and size() > 0.
  //   Calculates the size()-many coefficients in the representation
  //   of c as an affine combination of the points M.  The i-th computed
  //   coefficient coeffs[i] corresponds to the i-th point in M, or,
  //   in other words, to the point in S with index global_index(i).
  //   Complexity: O(dim^2)
  {
  public: // construction and deletion:

    Subspan(unsigned int dim, const PointAccessor& S, int i);
    // Constructs an instance representing the affine hull aff(M) of M={p},
    // where p is the point S[i] from S.
    //
    // Notice that S must not changed as long as this instance of
    // Subspan<Float> is in use.

    ~Subspan();

  public: // modification:

    void add_point(int global_index);
    void remove_point(unsigned int local_index);

  public: // access:

    unsigned int size() const
    {
      return r+1;
    }

    bool is_member(unsigned int i) const
    {
      return membership[i];
    }

    unsigned int global_index(unsigned int i) const
    {
      return members[i];
    }

    unsigned int any_member() const {
      return members[r];
    }

    template<typename RandomAccessIterator1,
    typename RandomAccessIterator2>
    Float shortest_vector_to_span(RandomAccessIterator1 p,
                                  RandomAccessIterator2 w);

    template<typename RandomAccessIterator1,
    typename RandomAccessIterator2>
    void find_affine_coefficients(RandomAccessIterator1 c,
                                  RandomAccessIterator2 coeffs);

  public: // debugging routines:

    Float representation_error();
    // Computes the coefficient representations of all points in the
    // (internally used) system (Q) and returns the maximal deviation
    // from the theoretical values.
    // Warning: This routine has running time O(dim^3).

  private: // private helper routines:

    void append_column();
    // Appends the new column u (which is a member of this instance) to
    // the right of "A = QR", updating Q and R.  It assumes r to still
    // be the old value, i.e., the index of the column used now for
    // insertion; r is not altered by this routine and should be changed
    // by the caller afterwards.
    // Precondition: r<dim

    void hessenberg_clear(unsigned int start);
    // Given R in lower Hessenberg form with subdiagonal entries 0 to
    // pos-1 already all zero, clears the remaining subdiagonal entries
    // via Givens rotations.

    void special_rank_1_update();
    // Update current QR-decomposition "A = QR" to
    // A + u [1,...,1] = Q' R'.

  private: // member fields:
    const PointAccessor &S;            // a const-reference to the set S
    std::vector<bool> membership;      // S[i] in M iff membership[i]
    const unsigned int dim;            // ambient dimension (not to be
    // confused with the rank r,
    // see below)

    // Entry i of members contains the index into S of the i-th point
    // in M.  The point members[r] is called the "origin."
    std::vector<unsigned int> members;

  private: // member fields for maintaining the QR-decomposition:
    Float **Q, **R;                    // (dim x dim)-matrices Q
    // (orthogonal) and R (upper
    // triangular); notice that
    // e.g.  Q[j][i] is the element
    // in row i and column j
    Float *u,*w;                       // needed for rank-1 update
    unsigned int r;                    // the rank of R (i.e. #points - 1)
  };

} // namespace SEB_NAMESPACE


// The point members[r] is called the origin; we use the following macro
// to increase the readibilty of the code:
#define SEB_AFFINE_ORIGIN S[members[r]]

namespace SEB_NAMESPACE {

  template<typename Float>
  inline void givens(Float& c, Float& s, const Float a, const Float b)
  //  Determine the Givens coefficients (c,s) satisfying
  //
  //     c * a + s * b = +/- (a^2 + b^2)
  //     c * b - s * a = 0
  //
  //  We don't care about the signs here for efficiency,
  //  so make sure not to rely on them anywhere.
  //
  //  Source: taken from "Matrix Computations" (2nd edition) by Gene
  //  H. B. Golub & Charles F. B. Van Loan (Johns Hopkins University
  //  Press, 1989), p. 216.
  {
    using std::abs;
    using std::sqrt;

    if (b == 0) {
      c = 1;
      s = 0;
    } else if (abs(b) > abs(a)) {
      const Float t = a / b;
      s = 1 / sqrt (1 + sqr(t));
      c = s * t;
    } else {
      const Float t = b / a;
      c = 1 / sqrt (1 + sqr(t));
      s = c * t;
    }
  }

  template<typename Float, class Pt, class PointAccessor>
  Subspan<Float, Pt, PointAccessor>::Subspan(unsigned int dim, const PointAccessor& S, int index)
  : S(S), membership(S.size()), dim(dim), members(dim+1)
  {
    // allocate storage for Q, R, u, and w:
    Q = new Float *[dim];
    R = new Float *[dim];
    for (unsigned int i=0; i<dim; ++i) {
      Q[i] = new Float[dim];
      R[i] = new Float[dim];
    }
    u = new Float[dim];
    w = new Float[dim];

    // initialize Q to the identity matrix:
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
        Q[i][j] = (i==j)? 1 : 0;

    members[r = 0] = index;
    membership[index] = true;
  }

  template<typename Float, class Pt, class PointAccessor>
  Subspan<Float, Pt, PointAccessor>::~Subspan()
  {
    for (unsigned int i=0; i<dim; ++i) {
      delete[] Q[i];
      delete[] R[i];
    }
    delete[] Q;
    delete[] R;
    delete[] u;
    delete[] w;
  }

  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::add_point(int index) {

    // compute S[i] - origin into u:
    for (unsigned int i=0; i<dim; ++i)
      u[i] = S[index][i] - SEB_AFFINE_ORIGIN[i];

    // appends new column u to R and updates QR-decomposition,
    // routine work with old r:
    append_column();

    // move origin index and insert new index:
    membership[index] = true;
    members[r+1] = members[r];
    members[r]   = index;
    ++r;
  }

  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::remove_point(const unsigned int local_index) {

    membership[global_index(local_index)] = false;

    if (local_index == r) {
      // origin must be deleted

      // We choose the right-most member of Q, i.e., column r-1 of R,
      // as the new origin.  So all relative vectors (i.e., the
      // columns of "A = QR") have to be updated by u:= old origin -
      // S[global_index(r-1)]:
      for (unsigned int i=0; i<dim; ++i)
        u[i] = SEB_AFFINE_ORIGIN[i] - S[global_index(r-1)][i];

      --r;
      special_rank_1_update();

    } else {
      // general case: delete column from R

      //  shift higher columns of R one step to the left
      Float *dummy = R[local_index];
      for (unsigned int j = local_index+1; j < r; ++j) {
        R[j-1] = R[j];
        members[j-1] = members[j];
      }
      members[r-1] = members[r];  // shift down origin
      R[--r] = dummy;             // relink trash column

      // zero out subdiagonal entries in R
      hessenberg_clear(local_index);
    }
  }

  template<typename Float, class Pt, class PointAccessor>
  template<typename RandomAccessIterator1,
  typename RandomAccessIterator2>
  Float Subspan<Float, Pt, PointAccessor>::
  shortest_vector_to_span(RandomAccessIterator1 p,
                          RandomAccessIterator2 w)
  {
    using std::inner_product;

    // compute vector from p to origin, i.e., w = origin - p:
    for (unsigned int i=0; i<dim; ++i)
      w[i] = SEB_AFFINE_ORIGIN[i] - p[i];

    // remove projections of w onto the affine hull:
    for (unsigned int j = 0; j < r; ++j) {
      const Float scale = inner_product(w,w+dim,Q[j],Float(0));
      for (unsigned int i = 0; i < dim; ++i)
        w[i] -= scale * Q[j][i];
    }

    return inner_product(w,w+dim,w,Float(0));
  }

  template<typename Float, class Pt, class PointAccessor>
  Float Subspan<Float, Pt, PointAccessor>::representation_error()
  {
    using std::abs;

    std::vector<Float> lambdas(size());
    Float max = 0;
    Float error;

    // cycle through all points in hull
    for (unsigned int j = 0; j < size(); ++j) {
      // compute the affine representation:
      find_affine_coefficients(S[global_index(j)],lambdas.begin());

      // compare coefficient of point #j to 1.0
      error = abs(lambdas[j] - 1.0);
      if (error > max) max = error;

      // compare the other coefficients against 0.0
      for (unsigned int i = 0; i < j; ++i) {
        error = abs(lambdas[i] - 0.0);
        if (error > max) max = error;
      }
      for (unsigned int i = j+1; i < size(); ++i) {
        error = abs(lambdas[i] - 0.0);
        if (error > max) max = error;
      }
    }

    return max;
  }

  template<typename Float, class Pt, class PointAccessor>
  template<typename RandomAccessIterator1,
  typename RandomAccessIterator2>
  void Subspan<Float, Pt, PointAccessor>::
  find_affine_coefficients(RandomAccessIterator1 p,
                           RandomAccessIterator2 lambdas)
  {
    // compute relative position of p, i.e., u = p - origin:
    for (unsigned int i=0; i<dim; ++i)
      u[i] = p[i] - SEB_AFFINE_ORIGIN[i];

    // calculate Q^T u into w:
    for (unsigned int i = 0; i < dim; ++i) {
      w[i] = 0;
      for (unsigned int k = 0; k < dim; ++k)
        w[i] += Q[i][k] * u[k];
    }

    // We compute the coefficients by backsubstitution.  Notice that
    //
    //     c = \sum_{i\in M} \lambda_i (S[i] - origin)
    //       = \sum_{i\in M} \lambda_i S[i] + (1-s) origin
    //
    // where s = \sum_{i\in M} \lambda_i.-- We compute the coefficient
    // (1-s) of the origin in the variable origin_lambda:
    Float origin_lambda = 1;
    for (int j = r-1; j>=0; --j) {
      for (unsigned int k=j+1; k<r; ++k)
        w[j] -= *(lambdas+k) * R[k][j];
      origin_lambda -= *(lambdas+j) = w[j] / R[j][j];
    }
    // The r-th coefficient corresponds to the origin (cf. remove_point()):
    *(lambdas+r) = origin_lambda;
  }

  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::append_column()
  // Appends the new column u (which is a member of this instance) to
  // the right of "A = QR", updating Q and R.  It assumes r to still
  // be the old value, i.e., the index of the column used now for
  // insertion; r is not altered by this routine and should be changed
  // by the caller afterwards.
  // Precondition: r<dim
  {
    //  compute new column R[r] = Q^T * u
    for (unsigned int i = 0; i < dim; ++i) {
      R[r][i] = 0;
      for (unsigned int k = 0; k < dim; ++k)
        R[r][i] += Q[i][k] * u[k];
    }

    //  zero all entries R[r][dim-1] down to R[r][r+1]
    for (unsigned int j = dim-1; j > r; --j) {
      //  j is the index of the entry to be cleared
      //  with the help of entry j-1

      //  compute Givens coefficients c,s
      Float c, s;
      givens (c,s,R[r][j-1],R[r][j]);

      //  rotate one R-entry (the other one is an implicit zero)
      R[r][j-1] = c * R[r][j-1] + s * R[r][j];

      //  rotate two Q-columns
      for (unsigned int i = 0; i < dim; ++i) {
        const Float a = Q[j-1][i];
        const Float b = Q[j][i];
        Q[j-1][i] =  c * a + s * b;
        Q[j][i]   =  c * b - s * a;
      }
    }
  }

  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::hessenberg_clear (unsigned int pos)
  // Given R in lower Hessenberg form with subdiagonal entries 0 to
  // pos-1 already all zero, clears the remaining subdiagonal entries
  // via Givens rotations.
  {
    //  clear new subdiagonal entries
    for (; pos < r; ++pos) {
      //  pos is the column index of the entry to be cleared

      //  compute Givens coefficients c,s
      Float c, s;
      givens (c,s,R[pos][pos],R[pos][pos+1]);

      //  rotate partial R-rows (of the first pair, only one entry is
      //  needed, the other one is an implicit zero)
      R[pos][pos] = c * R[pos][pos] + s * R[pos][pos+1];
      //  (then begin at posumn pos+1)
      for (unsigned int j = pos+1; j < r; ++j) {
        const Float a = R[j][pos];
        const Float b = R[j][pos+1];
        R[j][pos]   =  c * a + s * b;
        R[j][pos+1] =  c * b - s * a;
      }

      //  rotate Q-columns
      for (unsigned int i = 0; i < dim; ++i) {
        const Float a = Q[pos][i];
        const Float b = Q[pos+1][i];
        Q[pos][i]   =  c * a + s * b;
        Q[pos+1][i] =  c * b - s * a;
      }
    }
  }

  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::special_rank_1_update ()
  // Update current QR-decomposition "A = QR" to
  // A + u * [1,...,1] = Q' R'.
  {
    //  compute w = Q^T * u
    for (unsigned int i = 0; i < dim; ++i) {
      w[i] = 0;
      for (unsigned int k = 0; k < dim; ++k)
        w[i] += Q[i][k] * u[k];
    }

    //  rotate w down to a multiple of the first unit vector;
    //  the operations have to be recorded in R and Q
    for (unsigned int k = dim-1; k > 0; --k) {
      //  k is the index of the entry to be cleared
      //  with the help of entry k-1

      //  compute Givens coefficients c,s
      Float c, s;
      givens (c,s,w[k-1],w[k]);

      //  rotate w-entry
      w[k-1] = c * w[k-1] + s * w[k];

      //  rotate two R-rows;
      //  the first column has to be treated separately
      //  in order to account for the implicit zero in R[k-1][k]
      R[k-1][k]    = -s * R[k-1][k-1];
      R[k-1][k-1] *=  c;
      for (unsigned int j = k; j < r; ++j) {
        const Float a = R[j][k-1];
        const Float b = R[j][k];
        R[j][k-1] =  c * a + s * b;
        R[j][k]   =  c * b - s * a;
      }

      //  rotate two Q-columns
      for (unsigned int i = 0; i < dim; ++i) {
        const Float a = Q[k-1][i];
        const Float b = Q[k][i];
        Q[k-1][i] =  c * a + s * b;
        Q[k][i]   =  c * b - s * a;
      }
    }

    //  add w * (1,...,1)^T to new R
    //  which means simply to add u[0] to each column
    //  since the other entries of u have just been eliminated
    for (unsigned int j = 0; j < r; ++j)
      R[j][0] += w[0];

    //  clear subdiagonal entries
    hessenberg_clear(0);
  }

} // namespace SEB_NAMESPACE

#endif // SEB_SUBSPAN_H