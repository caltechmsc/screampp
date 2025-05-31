#ifndef SCREAM_ROTATION
#define SCREAM_ROTATION

#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_matrix.hpp"
#include "scream_vector.hpp"
#include <vector>


using namespace std;



class Scream_Rotation {

public:
  Scream_Rotation();
  Scream_Rotation(SCREAM_ATOM*, SCREAM_ATOM*, ScreamAtomV&);
  ~Scream_Rotation();

  SCREAM_ATOM* get_start_atom() { return start_atom;};
  SCREAM_ATOM* get_end_atom() {return end_atom;};

  ScreamAtomV get_atoms_to_be_rotated() { return this->atoms_to_be_rotated;};

private:
  SCREAM_ATOM* start_atom;
  SCREAM_ATOM* end_atom;

  ScreamAtomV atoms_to_be_rotated;
  ScreamAtomV atoms_to_be_rotated_tmp;

  /* Routines. */
  ScreamAtomV _return_atoms_to_be_rotated(SCREAM_ATOM*, SCREAM_ATOM*, ScreamAtomV&);
  ScreamAtomV _return_atoms_to_be_rotated_DFS(SCREAM_ATOM*, ScreamAtomV&, ScreamAtomV&);

};




class Scream_Rotation_List {
public:

  vector<Scream_Rotation> sc_rot_v;

private:

};

#endif /* SCREAM_ROTATION */

