/* scurvy_matrix.h
 * 
 * Header file for matrix manipulation code used in backbone matching in scream.  3D matrix operations.
 *
 * Coypright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 */

/* scream_matrix.hpp Header file for matrix manipulation functions */


#ifndef SCREAM_MATRIX_HPP
#define SCREAM_MATRIX_HPP

#include <iostream>
#include "scream_vector.hpp"

/** ScreamMatrix is the matrix class used by SCREAM.  
 * Most common arithmetic operations that involved 3x3 matrices involved are included.  This includes multiplying a matrix with a vector--SreamMatrix recognizes ScreamVector as a legal procedure.  
 */

class ScreamMatrix {

public: 
  ScreamMatrix();
  ScreamMatrix(const ScreamVector&, const ScreamVector&, const ScreamVector&);  ///< initialize by column vectors
  ScreamMatrix(const ScreamMatrix&);
  ~ScreamMatrix();
  ScreamMatrix & operator=(const ScreamMatrix &);
  
  ScreamVector& operator[](int);  ///< returns row vector
  const ScreamVector & operator[](int) const;
  ScreamVector& getCol(int);
  const ScreamVector& getCol(int) const;
  
  ScreamMatrix operator+(const ScreamMatrix &) const;	///< Addition operator.
  ScreamMatrix operator-(const ScreamMatrix &) const;	///< Substraction operator.
  ScreamMatrix operator*(double) const; ///< Scalar multiplication operator.
  ScreamMatrix operator/(double) const; ///< Scalar division operator.

  ScreamVector operator*(const ScreamVector &) const;	///< Left multiply a vector.  
  ScreamMatrix operator*(const ScreamMatrix &) const;	///< Left multiply a matrix.

  double det() const;		///< Returns value of determinant.
  ScreamMatrix inverse() const;	///< Returns inverse of matrix.
  ScreamMatrix transpose() const; ///< Returns transpose of matrix.

  /** Returns the matrix that aligns the second vector with the first one.
   * alignTwoVectors takes in Vectors fixed and tobealigned and generate the transformation matrix that takes the coords of vector tobealigned to vector fixed.  The two vectors are assumed to have the same origin.
   */
  ScreamMatrix alignTwoVectors(const ScreamVector& fixedVec, const ScreamVector& toBeAlignVec);
  
  /** alignZ takes in m vecZ.  returns the transformation matrix that takes the Cartesian coordinate system (meaning x = [1 0 0], y = [0 1 0], z = [0 0 1]) to a coordinate system with a z axis aligned with vecZ.
   */
  ScreamMatrix alignWithZ(const ScreamVector &) const;
  /** rotZ takes in angle.  returns the rotational matrix that rotates around the z = [1 0 0] axis by the specified angle. 
   * Right-handed rule rotation.  I.e. pointing thumb in Z = [0 0 1] direction, positive rotation is defined in the direction of curling fingers.
   */
  ScreamMatrix rotAboutZ(double angle) const;

  /** rotAboutV takes in a ScreamVector and an angle and returns the rotation matrix that rotates by angle about ScreamVector.
   * Directionality is taken into account.  The usual right-handed rule applies; i.e. positive rotation is defined as right-handed rotation with the vector pointing in the direction of the thumb.  I.e. clockwise when looking in from the start of the vector ray.
   */
  ScreamMatrix rotAboutV(const ScreamVector &, double) const;

  void printMe();


private:
  // need to update both member variables everytime an operation is done.
  ScreamVector* rowVectors;
  ScreamVector* colVectors;
} ;

#endif

