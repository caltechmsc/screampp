/**\file scream_E_functionals_coulomb
 * 
 */

#ifndef SCREAM_E_FUNCTIONALS_COULOMB
#define SCREAM_E_FUNCTIONALS_COULOMB

#include <math.h>
#include <vector>
#include <string>

#include "scream_atom.hpp"

using namespace std;


/* Notes: 
   Need to take care of:
   1) Bad initial geometries for Morse potential, since E_Morse goes to -inf when R approaches 0.
*/

class SCREAM_Coulomb_OBJ {

public:
  SCREAM_Coulomb_OBJ();
  SCREAM_Coulomb_OBJ(string, double);
  ~SCREAM_Coulomb_OBJ();

  double calc_Coulomb(SCREAM_ATOM*, SCREAM_ATOM*, double = 1);
  double calc_Coulomb_normal(SCREAM_ATOM*, SCREAM_ATOM*, double = 1);
  double calc_Coulomb_distance_dielectric(SCREAM_ATOM*, SCREAM_ATOM*, double = 1);
  
  void read_param_line(string); 
  void read_param_file(string); ///< reads an mpsim style parameter file.

  double getEpsilon() {return this->epsilon;};
  void set_dielectric(double epsilon) {this->epsilon = epsilon;};
  
  void set_normal_mode() { this->mode = 1;};
  void set_distance_dependent_mode() { this->mode = 2;};


private:
  int mode; // {1: Normal  2: Distance dependent dielectric.}
  double epsilon; // applies to both normal coulomb or distance dependent dielectric
  double R_on, R_off;
};





#endif /* SCREAM_E_FUNCTIONALS_COULOMB */
