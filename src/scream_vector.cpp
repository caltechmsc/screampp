/* scream_Vector.cpp
 *
 * Source code for ScreamVector class.
 *
 * Copyright (c) 2003.  Victor Wai Tak Kam.
 */


#include <math.h>
#include <iostream>
#include "scream_atom.hpp"
#include "scream_vector.hpp"

ScreamVector::ScreamVector() {
  coords = new double[3];
  coords[0] = coords[1] = coords[2] = 0.0;
}

ScreamVector::ScreamVector(double px, double py, double pz) {
  coords = new double[3];
  coords[0] = px;
  coords[1] = py;
  coords[2] = pz;
}
  
ScreamVector::ScreamVector(const SCREAM_ATOM* inAtom) {

  coords = new double[3];
  coords[0] = inAtom->x[0];
  coords[1] = inAtom->x[1];
  coords[2] = inAtom->x[2];
}

ScreamVector::ScreamVector(const ScreamVector& inVec) {
  
  coords = new double[3];
  for (int i = 0; i < 3; ++i) {
    coords[i] = inVec.coords[i];
  }

}

ScreamVector::~ScreamVector() {
 delete [] coords;
}

double & ScreamVector::operator[](int i) {
  // check index before continuing
  if (i<0 || i>2) {
    cerr << "Error in ScreamVector limits: " << i << " is a bad index\n";
  }

  return coords[i];
}

const double & ScreamVector::operator[](int i) const {
  // check index before continuing
  if (i<0 || i>2) {
    cerr << "Error in ScreamVector limits: " << i << " is a bad index\n";
  }

  return coords[i];
}



ScreamVector & ScreamVector::operator=(const ScreamVector& inVec) {
  if (this == &inVec)
    return *this;
  for (int i =0; i < 3; ++i) {
    coords[i] = inVec.coords[i];
  }
  return *this;
}

bool ScreamVector::operator==(const ScreamVector& v) const {
//   if ( (v[0] - 0.00001) <= this->coords[0] and this->coords[0] <= (v[0] + 0.00001) and
//        (v[1] - 0.00001) <= this->coords[1] and this->coords[1] <= (v[1] + 0.00001) and
//        (v[2] - 0.00001) <= this->coords[2] and this->coords[2] <= (v[2] + 0.00001) ) {
  if ( v[0] == this->coords[0] and
       v[1] == this->coords[1] and 
       v[2] == this->coords[2] ) {
    return true;
  } else {
    return false;
  }
}

bool ScreamVector::operator<(const ScreamVector& v) const {
  if (this->coords[0] < v[0]) {
    return true;
  } else if (this->coords[0] == v[0] and this->coords[1] < v[1]) {
    return true;
  } else if (this->coords[0] == v[0] and this->coords[1] == v[1] and this->coords[2] < v[2]) {
    return true;
  } else {
    return false;
  }
}

//ostream & operator<<(ostream & os, ScreamVector & inVec) {
//  os << "[";
//  for (int i = 0; i<3; ++i) {
//    os << inVec.coords[i] << " ";
//  }
//  os << "]";
//} 

ScreamVector ScreamVector::operator+(const ScreamVector& inVec) const {
  return ScreamVector(coords[0]+inVec.coords[0], coords[1]+inVec.coords[1], coords[2]+inVec.coords[2]);
}

ScreamVector ScreamVector::operator-(const ScreamVector& inVec) const {
  return ScreamVector(coords[0]-inVec.coords[0], coords[1]-inVec.coords[1], coords[2]-inVec.coords[2]);
}

ScreamVector ScreamVector::operator*(double factor) const {
  return ScreamVector(coords[0]*factor, coords[1]*factor, coords[2]*factor);
}

ScreamVector ScreamVector::operator/(double factor) const {
  return ScreamVector(coords[0]/factor, coords[1]/factor, coords[2]/factor);
}

ScreamVector& ScreamVector::operator+=(const ScreamVector& inVec) {
  this->coords[0] += inVec.coords[0];
  this->coords[1] += inVec.coords[1];
  this->coords[2] += inVec.coords[2];
  
  return *this;

}

double ScreamVector::magnitude() const {
  double sum = 0;
  for (int i = 0; i<3; ++i) {
    sum+= coords[i]*coords[i];
  }
  return sqrt(sum);
}

ScreamVector ScreamVector::normalizedVector() const {
  double magnitude1 = magnitude();
  return ScreamVector(coords[0] / magnitude1, coords[1] / magnitude1, coords[2] / magnitude1);
}

void ScreamVector::normalize() {
  double magnitude1 = magnitude();
  for (int i = 0; i<3; ++i) {
    coords[i] = coords[i] / magnitude1;
  }
}

ScreamVector ScreamVector::cross(const ScreamVector& inVec) const {
  return ScreamVector(coords[1]*inVec[2] - coords[2]*inVec[1],
		      -coords[0]*inVec[2] + coords[2]*inVec[0],
		      coords[0]*inVec[1] - coords[1]*inVec[0] );
}

double ScreamVector::dot(const ScreamVector& inVec) const {
  return coords[0]*inVec[0] + coords[1]*inVec[1] + coords[2]*inVec[2];
}

double ScreamVector::angleBtwn(const ScreamVector& V) const {
  // no directionality involved.  return value always between 0 and 180.
  double cos_value = this->dot(V) / (this->magnitude() * V.magnitude() );
  return double(acos(cos_value))*180/3.1415926535;

}



void ScreamVector::printMe() {
  cout << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
}

