/* scream_Vector.hpp
 * 
 * Header file for vectors used in scream routines.  Cartesion vectors in 3 D. 
 * 
 * Copyright (c) 2003.  Victor Wai Tak Kam.  
 */

#ifndef SCREAM_VECTOR_H_
#define SCREAM_VECTOR_H_

#include <iostream>
#include "scream_atom.hpp"

class ScreamVector {

public:
  ScreamVector();
  ScreamVector(double, double, double);
  ScreamVector(const SCREAM_ATOM*);
  ScreamVector(const ScreamVector&);
  ~ScreamVector();
  double & operator[](int i);
  const double & operator[](int i) const;
  ScreamVector & operator=(const ScreamVector &);

  bool operator==(const ScreamVector &) const; // If two vectors are identical.
  bool operator<(const ScreamVector& ) const; // For ordering purposes-- v1 < v2 if v1[0] < v2[0], v1[1] < v2[1], v1[2] < v2[2].

  //ostream & operator<<(ostream &, const ScreamVector &);
  
  ScreamVector operator+(const ScreamVector &) const;
  ScreamVector operator-(const ScreamVector &) const;
  ScreamVector operator*(double) const;
  ScreamVector operator/(double) const;
 
  ScreamVector& operator+=(const ScreamVector &);

  double magnitude() const;
  ScreamVector normalizedVector() const;     ///< returns a normalized vector
  void normalize();                ///< normalizes this vector
  
  ScreamVector cross(const ScreamVector&) const;
  double dot(const ScreamVector& ) const;
  
  double angleBtwn(const ScreamVector&) const;  // returns angle between the two vectors.  

  void printMe();

private:
  double* coords;  // really, coords[3]

};

#endif 
