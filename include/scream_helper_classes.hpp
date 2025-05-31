//! scream_helper_classes

#ifndef SCREAMHELPER_CLASSES
#define SCREAMHELPER_CLASSES

#include <vector>
using namespace std;

#include "scream_atom.hpp"
#include "scream_vector.hpp"


class RotamerAxis {

public:
  RotamerAxis() {};
  RotamerAxis(SCREAM_ATOM*, SCREAM_ATOM*);
  ~RotamerAxis();
  SCREAM_ATOM* const getFirstAtom(); 
  SCREAM_ATOM* const getSecondAtom();
  
  
private:
  SCREAM_ATOM *a1, *a2;
  ScreamVector* axisDirection;


};



#endif /* SCREAMHELPER_CLASSES */
