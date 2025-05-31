/* scurvy_matrix.h
 * Code for matrix routines used in MPSim-SCREAM.
 * Coypright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 */

/* This file contains some routines from fsm_matrix.h.
 */


/* scream_matrix.cpp Source file for matrix manipulation functions */

//#include <stdlib.h>
//#include <stdio.h>
//#include <stdarg.h>
//#include <math.h>
#include "stdlib.h"
#include "stdio.h"
#include "stdarg.h"
#include "math.h"

#include <iostream>
#include "scream_vector.hpp"
#include "scream_matrix.hpp"


ScreamMatrix::ScreamMatrix() {
  rowVectors = new ScreamVector[3];
  colVectors = new ScreamVector[3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j <3; ++j) {
      rowVectors[i][j] = 0.0;
      colVectors[i][j] = 0.0;
    }
  }
}

ScreamMatrix::ScreamMatrix(const ScreamVector& V1, const ScreamVector& V2, const ScreamVector& V3) { 
  
  rowVectors = new ScreamVector[3];
  colVectors = new ScreamVector[3];

  for (int i = 0; i<3; ++i) {
    colVectors[0][i] = V1[i];
    rowVectors[i][0] = V1[i];
  }
  for (int i = 0; i<3; ++i) {
    colVectors[1][i] = V2[i];
    rowVectors[i][1] = V2[i];
  }
  for (int i = 0; i<3; ++i) {
    colVectors[2][i] = V3[i];
    rowVectors[i][2] = V3[i];
  }
}

ScreamMatrix::ScreamMatrix(const ScreamMatrix& inMat) {
  rowVectors = new ScreamVector[3];
  colVectors = new ScreamVector[3];

  for (int i = 0; i<3; ++i) {
    for (int j = 0; j<3; ++j) {
      rowVectors[i][j] = inMat.rowVectors[i][j];
      colVectors[j][i] = rowVectors[i][j];
    }
  }
}

ScreamMatrix::~ScreamMatrix() {
  delete [] rowVectors;
  delete [] colVectors;
}


ScreamVector& ScreamMatrix::operator[](int i) {
  // check index before continuing
  if (i<0 || i>2) {
    cerr << "Error in ScreamMatrix limits: " << i << " is a bad index\n";
  }
  return rowVectors[i];
}

ScreamMatrix& ScreamMatrix::operator=(const ScreamMatrix & inMat) {
  if (this == &inMat) 
    return *this;
  for (int i =0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      rowVectors[i][j] = colVectors[i][j] = inMat.rowVectors[i][j];
    }
  }
  return *this;
}
      

const ScreamVector& ScreamMatrix::operator[](int i) const {
  // check index before continuing
  if (i<0 || i>2) {
    cerr << "Error in ScreamMatrix limits: " << i << " is a bad index\n";
  }
  return rowVectors[i];
}

ScreamVector& ScreamMatrix::getCol(int i) {
  // check index before continuing
  if (i<0 || i>2) {
    cerr << "Error in ScreamMatrix limits: "<< i << " is a bad index\n";
  }
  return colVectors[i];
}

const ScreamVector& ScreamMatrix::getCol(int i) const {
  // check index before continuing
  if (i<0 || i>2) {
    cerr << "Error in ScreamMatrix limits: "<< i << " is a bad index\n";
  }
  return colVectors[i];
}


ScreamMatrix ScreamMatrix::operator+(const ScreamMatrix & M ) const {
  ScreamMatrix newMat;
  for (int i = 0; i<3; ++i) {
    newMat[i] = rowVectors[i] + M[i];
    newMat.getCol(i) = colVectors[i] + M.getCol(i);
  }
  return newMat;

  /*  ScreamMatrix* newMat = new ScreamMatrix;
  for (int i = 0; i<3; ++i) {
    (*newMat)[i] = rowVectors[i] + M[i];
    (*newMat).getCol(i) = colVectors[i] + M.getCol(i);
  }	   
  return *newMat;
  */
}
  

ScreamMatrix ScreamMatrix::operator-(const ScreamMatrix & M ) const {
  ScreamMatrix newMat;
  for (int i = 0; i<3; ++i) {
    newMat[i] = rowVectors[i] - M.rowVectors[i];
    newMat.getCol(i) = colVectors[i] - M.getCol(i);
  }
  return newMat;

  /*  ScreamMatrix* newMat = new ScreamMatrix;
  for (int i = 0; i<3; ++i) {
    (*newMat)[i] = rowVectors[i] - M.rowVectors[i];
    (*newMat).getCol(i) = colVectors[i] - M.getCol(i);
  }	   
  return *newMat;
  */
}


ScreamMatrix ScreamMatrix::operator*(double factor) const {
/********  ScreamMatrix newMat;
  for (int i=0; i<3; ++i) {
    newMat[i] = rowVectors[i] * factor;
    
  }
*/
    ScreamMatrix* newMat = new ScreamMatrix;
  for (int i=0; i<3; ++i) {
    (*newMat)[i] = rowVectors[i] * factor;
    (*newMat).getCol(i) = colVectors[i] * factor;
  }
  return *newMat;
}

ScreamMatrix ScreamMatrix::operator/(double factor) const {
  if (factor == 0.0) {
    cerr << "Divisor is zero.\n";
  }
  return this->operator*(1/factor);
}
    
ScreamVector ScreamMatrix::operator*(const ScreamVector & V) const {
  return ScreamVector(rowVectors[0].dot(V),
		      rowVectors[1].dot(V),
		      rowVectors[2].dot(V));
}

ScreamMatrix ScreamMatrix::operator*(const ScreamMatrix & M ) const {
  ScreamMatrix newMat;
  for (int i =0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      newMat[i][j] = rowVectors[i].dot(M.getCol(j));   
      newMat.getCol(j)[i] = rowVectors[i].dot(M.getCol(j));
    }
  }

  return newMat;
  /*  ScreamMatrix* newMat = new ScreamMatrix;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      (*newMat)[i][j] = rowVectors[i].dot(M.getCol(j));   
      (*newMat).getCol(j)[i] = rowVectors[i].dot(M.getCol(j));
    }
  }
  return *newMat;
  */
  
}

double ScreamMatrix::det() const {

  return colVectors[0].dot(colVectors[1].cross(colVectors[2]));

}

ScreamMatrix ScreamMatrix::inverse() const {

  return ScreamMatrix(rowVectors[1].cross(rowVectors[2]) / det(),
		      rowVectors[2].cross(rowVectors[0]) / det(),
		      rowVectors[0].cross(rowVectors[1]) / det());
}


ScreamMatrix ScreamMatrix::transpose() const {
  return ScreamMatrix(rowVectors[0], rowVectors[1], rowVectors[2]);
}


ScreamMatrix ScreamMatrix::alignTwoVectors(const ScreamVector& fixedVec, const ScreamVector& toBeAlignVec) {
  
  double angle = fixedVec.angleBtwn(toBeAlignVec);

  ScreamVector newZAxis = fixedVec.cross(toBeAlignVec);
  /*
  ScreamMatrix transformToNewZ = alignWithZ(newZAxis);
  //  transformToNewZ.printMe();
  ScreamMatrix rotationAroundNewZ = rotAboutZ(angle);
  //  rotationAroundNewZ.printMe();
  */
  
  ScreamMatrix transform = rotAboutV(newZAxis, -angle);

  /*
  ScreamMatrix overallAlignmentMatrix = transformToNewZ * rotationAroundNewZ * transformToNewZ.inverse();
  
  return ScreamMatrix(overallAlignmentMatrix);
  */
  return transform;

}


ScreamMatrix ScreamMatrix::alignWithZ(const ScreamVector& vectorZ) const {
  
  ScreamVector outX;
  ScreamVector outY;
  ScreamVector outZ = vectorZ.normalizedVector();
  
  double angleBetween = outZ.angleBtwn(ScreamVector(0,0,1));

  if (angleBetween == 0) {
    outX = ScreamVector(1.0,0.0,0.0);
    outY = ScreamVector(0.0,1.0,0.0);
    outZ = ScreamVector(0.0,0.0,1.0);
  } else if ( fabs(fabs(angleBetween) - 180) < 0.00002) {
    outX = ScreamVector(-1.0,0.0,0.0);
    outY = ScreamVector(0.0,1.0,0.0);
    outZ = ScreamVector(0.0,0.0,-1.0);
  } else {
    outX = outZ.cross(ScreamVector(0.0,0.0,1.0));
    outY = outZ.cross(outX);
    outZ = outZ;
  }

  ScreamMatrix R(outX, outY, outZ);
  return R.inverse();

}


ScreamMatrix ScreamMatrix::rotAboutZ(double angle) const {

  angle = angle * 3.1415926535 / 180;

  return ScreamMatrix( ScreamVector(cos(angle), sin(angle), 0),
		       ScreamVector(-sin(angle), cos(angle), 0),
		       ScreamVector(0, 0, 1) );
}

ScreamMatrix ScreamMatrix::rotAboutV(const ScreamVector & V, double angle) const {

  // Algorithm: first generate rotation vector about Z.

  ScreamMatrix rotZ = rotAboutZ(angle); // angle is in degrees.  rotAboutZ converts it to radians.
  // Then do the T * rotZ * T_inv thing. 

  ScreamMatrix T = alignWithZ(V);
  
  return ScreamMatrix(T.inverse() * rotZ * T);

}


void ScreamMatrix::printMe() {
  rowVectors[0].printMe(); cout << "\n";
  rowVectors[1].printMe(); cout << "\n";
  rowVectors[2].printMe(); cout << "\n";
}
