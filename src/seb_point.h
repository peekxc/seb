// Synopsis: Simple point class
//
// Authors: Martin Kutz <kutz@math.fu-berlin.de>,
//          Kaspar Fischer <kf@iaeth.ch>

#ifndef SEB_POINT_H
#define SEB_POINT_H

#include <vector>

namespace SEB_NAMESPACE {

  template<typename Float> // A simple class representing a d-dimensional point.
  class Point {
  public: // types:
    typedef typename std::vector<Float>::const_iterator Const_iterator;
    typedef typename std::vector<Float>::iterator Iterator;

  public: // construction and destruction:

  	// Constructs a d-dimensional point with undefined coordinates.
    Point(int d): c(d){ }

  	// Constructs a d-dimensional point with Cartesian center
    // coordinates [first,first+d).
    template<typename InputIterator>
    Point(int d, InputIterator first) : c(first,first+d) { }

  public: // access:

  	// Returns a const-reference to the i-th coordinate.
    const Float& operator[](unsigned int i) const {
      return c[i];
    }
  	
  	// Returns a reference to the i-th coordinate.
    Float& operator[](unsigned int i){
      return c[i];
    }

    // Returns a const-iterator to the first of the d Cartesian coordinates.
    Const_iterator begin() const{
      return c.begin();
    }

    // Returns the past-the-end iterator corresponding to begin().
    Const_iterator end() const {
      return c.end();
    }

  private: // member fields:
    std::vector<Float> c;       // Cartesian center coordinates
  };
} // namespace SEB_NAMESPACE

#endif // SEB_POINT_H

