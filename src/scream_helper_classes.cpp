#include "scream_helper_classes.hpp"

RotamerAxis::RotamerAxis(SCREAM_ATOM* a1, SCREAM_ATOM *a2) {

  this->a1 = a1;
  this->a2 = a2;

  double a1_x = a1->x[0];
  double a2_x = a2->x[0];
  double a1_y = a1->x[1];
  double a2_y = a2->x[1];
  double a1_z = a1->x[2];
  double a2_z = a2->x[2];

  this->axisDirection = new ScreamVector(a2_x - a1_x, a2_y - a1_y, a2_z - a1_z);

}

RotamerAxis::~RotamerAxis() {
  delete axisDirection;
}


SCREAM_ATOM* const RotamerAxis::getFirstAtom() {

  return this->a1;
  
}

SCREAM_ATOM* const RotamerAxis::getSecondAtom() {

  return this->a2;

}
