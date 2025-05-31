/**\file scream_E_functionals_hb
 *
 */


#ifndef SCREAM_E_FUNCTIONALS_HB
#define SCREAM_E_FUNCTIONALS_HB

#include <math.h>
#include <vector>
#include <string>

#include "scream_atom.hpp"
#include "AtomResInfo.hpp"

using namespace std;

class SCREAM_HB_fields {
public:
  SCREAM_HB_fields();
  SCREAM_HB_fields(double, double);
  ~SCREAM_HB_fields();

  double DE;
  double RE;
};

class HB_delta_fields {
public:
  HB_delta_fields();
  HB_delta_fields(double, double);
  ~HB_delta_fields();

  HB_delta_fields& operator=(const HB_delta_fields&);

  double mu, sigma;
};


class SCREAM_HB_OBJ {
public:

  SCREAM_HB_OBJ();
  ~SCREAM_HB_OBJ();

  double calc_HB_Dre(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*); ///< give: Acceptor, Hydrogen, Donor.  LJ 12-10.
  double calc_HB_CHARMM(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*);

  double calc_Scream_HB(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*); ///< give: Acceptor, Donor.  Hydrogen not necessary.  LJ 12-10, with buffer zone specified in HB fields.
  double calc_full_delta_HB(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*, double); 
  double calc_flat_delta_HB(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*, double);
  double calc_residue_delta_HB(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*);
  double calc_scaled_inner_wall_HB(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*, double);

  void read_param_line(string);
  void read_param_file(string); // now obsolete as all EE classes are organized by scream_EE
  void read_Scream_delta_file(string);

  SCREAM_HB_fields* get_HB_fields(SCREAM_ATOM*, SCREAM_ATOM*);
  double get_RE(SCREAM_ATOM*, SCREAM_ATOM*);
  double get_DE(SCREAM_ATOM*, SCREAM_ATOM*);

  double _calc_angle(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*); ///< A, H, D.

  /* Mapping of Atom Labels to Integers for fast access (integer cmp instead of string cmp)  */
  map<string, int> hb_atom_type_mapping; ///< For a heavy atom type that shows up in the HB list of a par file, like O_3, N_R, S_3, an associated integer is created.  
  map< pair<int, int>, SCREAM_HB_fields*> hb_dict; ///< pair<int, int>: hb atom type from mapping.
  //map< pair<string, string>, SCREAM_HB_fields*> hb_dict;

  vector<string> returnAllHBTypes() const;

  double R_on, R_off;
  double theta_on, theta_off;
  double _optimize_10_12(double, double) const;  ///< Returns the value which minimizes the 10-12 functional form for the domain between the two specified doubles.

private:
  map<string, map<AtomResInfo, HB_delta_fields*> * > hb_delta_library_dict;

};

class SCREAM_HB_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_HB_BASE_FUNCTIONAL_OBJ(SCREAM_HB_OBJ*);
  virtual ~SCREAM_HB_BASE_FUNCTIONAL_OBJ();

  virtual double operator()(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*) = 0;
  
  SCREAM_HB_OBJ* scream_hb_obj;


private:


};

class SCREAM_calc_full_delta_HB : public SCREAM_HB_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_HB(SCREAM_HB_OBJ*, double);
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*);
  
  double n_sigma;


};

class SCREAM_calc_flat_delta_HB : public SCREAM_HB_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_HB(SCREAM_HB_OBJ*, double);
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*);

  double delta;


};


#endif /* SCREAM_E_FUNCTIONALS_HB */
