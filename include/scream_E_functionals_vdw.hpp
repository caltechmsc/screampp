/**\file scream_E_functionals_vdw
 * 
 */

#ifndef SCREAM_E_FUNCTIONALS_VDW
#define SCREAM_E_FUNCTIONALS_VDW

#include <math.h>
#include <vector>
#include <string>

#include "scream_atom.hpp"
#include "AtomResInfo.hpp"

using namespace std;


/* Notes: 
   Need to take care of:
   1) Bad initial geometries for Morse potential, since E_Morse goes to -inf when R approaches 0.
*/
class VDW_fields {
public:
  VDW_fields();
  VDW_fields(double, double, double);
  ~VDW_fields();

  VDW_fields& operator=(const VDW_fields&);

  double RNB;
  double DENB;
  double SCALE;

};

class VDW_delta_fields {
public:
  VDW_delta_fields();
  VDW_delta_fields(double, double);
  ~VDW_delta_fields();

  VDW_delta_fields& operator=(const VDW_delta_fields&);

  double mu, sigma;

};


class SCREAM_VDW_OBJ {

public:
  SCREAM_VDW_OBJ();
  ~SCREAM_VDW_OBJ();

  // The following routines will initialize this object.

  void read_param_line(string);///< reads one line of parameter file
  void read_param_file(string); ///< reads an mpsim style parameter file.
  void read_Scream_delta_file(string); ///< reads a SCREAM format neighborhood delta file 
  

  double calc_VDW_6_8(SCREAM_ATOM*, SCREAM_ATOM*);
  double calc_VDW_6_9(SCREAM_ATOM*, SCREAM_ATOM*);
  double calc_VDW_6_10(SCREAM_ATOM*, SCREAM_ATOM*);
  double calc_VDW_6_12(SCREAM_ATOM*, SCREAM_ATOM*);
  double calc_VDW_X6(SCREAM_ATOM*, SCREAM_ATOM*);
  double calc_VDW_Morse(SCREAM_ATOM*, SCREAM_ATOM*);
  
  double calc_Scream_VDW_6_12(SCREAM_ATOM*, SCREAM_ATOM*); // experimental... will become obsolete.

  double calc_full_delta_VDW_6_12(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using the delta region.
  double calc_flat_delta_VDW_6_12(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using delta where only one flat delta value is used.

  double calc_VDW_6_12_scaled_inner_wall(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy with a scaled inner wall.

  double calc_full_delta_asym_VDW_6_12(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using the delta region.
  double calc_flat_delta_asym_VDW_6_12(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using delta where only one flat delta value is used.
  double calc_full_delta_asym_VDW_6_11(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using the delta region.
  double calc_flat_delta_asym_VDW_6_11(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using delta where only one flat delta value is used.
  double calc_full_delta_asym_VDW_6_10(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using the delta region.
  double calc_flat_delta_asym_VDW_6_10(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using delta where only one flat delta value is used.
  double calc_full_delta_asym_VDW_6_9(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using the delta region.
  double calc_flat_delta_asym_VDW_6_9(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using delta where only one flat delta value is used.
  double calc_full_delta_asym_VDW_6_8(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using the delta region.
  double calc_flat_delta_asym_VDW_6_8(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using delta where only one flat delta value is used.
  double calc_full_delta_asym_VDW_6_7(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using the delta region.
  double calc_flat_delta_asym_VDW_6_7(SCREAM_ATOM*, SCREAM_ATOM*, double); // calculation of energy using delta where only one flat delta value is used.

  double calc_full_delta_asym_X6(SCREAM_ATOM*, SCREAM_ATOM*, double);
  double calc_flat_delta_asym_X6(SCREAM_ATOM*, SCREAM_ATOM*, double);

  double get_RNB(string ff_type);
  double get_DENB(string ff_type);
  double get_SCALE(string ff_type);
  VDW_fields* get_VDW_fields(string ff_type);
  VDW_delta_fields* get_VDW_delta_fields(string, AtomResInfo) const;

  double _geom_mean(double i, double j) { return sqrt(i*j); };
  double _arith_mean(double i, double j) { return (i+j)/2; };
  double _harm_mean(double i, double j) { return 1/(1/i + 1/j);};

private:
  //  map<string, double> RNB; ///< RNB values for string FF_type.
  //  map<string, double> DENB; ///< DENB values for string FF_type.
  map<string, VDW_fields*> vdw_dict; 
  map<string, map<AtomResInfo, VDW_delta_fields*>* > vdw_delta_library_dict;

  //double _locally_optimize_vdw_functional(double, double); // returns value that optimizes (locally) a LJ type VDW function.  LJ type: -pow(r,-6) + pow(r,12), with just one minimum over [0,inf).
  double _symmetric_optimize_LJ_vdw_functional(double, double); // returns value that optimizes (locally) a LJ type VDW function.  LJ type: -pow(r,-6) + pow(r,12), with just one minimum over [0,inf).
  double _asymmetric_optimize_LJ_vdw_functional(double, double); // returns 

  void _checkAndInitFullAtomVdwAndDeltaFields(SCREAM_ATOM*, SCREAM_ATOM*, double);
  void _checkAndInitFlatAtomVdwAndDeltaFields(SCREAM_ATOM*, SCREAM_ATOM*, double);
};

class Morse_VDW : public SCREAM_VDW_OBJ {
public:
  double calc_VDW(SCREAM_ATOM*, SCREAM_ATOM*);
private:
  double scaling_factor;
};


class SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ(SCREAM_VDW_OBJ*);
  virtual ~SCREAM_VDW_BASE_FUNCTIONAL_OBJ();
  virtual double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const = 0;

  SCREAM_VDW_OBJ* scream_vdw_obj;

private:

};

class SCREAM_calc_full_delta_VDW_6_12 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_VDW_6_12(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_VDW_6_12();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

class SCREAM_calc_flat_delta_VDW_6_12 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_VDW_6_12(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_VDW_6_12();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;
};
 
class SCREAM_calc_flat_delta_asym_VDW_6_12 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_asym_VDW_6_12(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_asym_VDW_6_12();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;

};

class SCREAM_calc_full_delta_asym_VDW_6_12 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_asym_VDW_6_12(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_asym_VDW_6_12();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};


/* 6-11 */

class SCREAM_calc_full_delta_VDW_6_11 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_VDW_6_11(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_VDW_6_11();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

class SCREAM_calc_flat_delta_VDW_6_11 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_VDW_6_11(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_VDW_6_11();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;
};
 
class SCREAM_calc_flat_delta_asym_VDW_6_11 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_asym_VDW_6_11(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_asym_VDW_6_11();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;

};

class SCREAM_calc_full_delta_asym_VDW_6_11 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_asym_VDW_6_11(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_asym_VDW_6_11();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};



/* 6-10 */
class SCREAM_calc_full_delta_VDW_6_10 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_VDW_6_10(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_VDW_6_10();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

class SCREAM_calc_flat_delta_VDW_6_10 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_VDW_6_10(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_VDW_6_10();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;
};
 
class SCREAM_calc_flat_delta_asym_VDW_6_10 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_asym_VDW_6_10(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_asym_VDW_6_10();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;

};

class SCREAM_calc_full_delta_asym_VDW_6_10 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_asym_VDW_6_10(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_asym_VDW_6_10();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};


/* 6-9 */
class SCREAM_calc_full_delta_VDW_6_9 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_VDW_6_9(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_VDW_6_9();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

class SCREAM_calc_flat_delta_VDW_6_9 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_VDW_6_9(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_VDW_6_9();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;
};
 
class SCREAM_calc_flat_delta_asym_VDW_6_9 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_asym_VDW_6_9(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_asym_VDW_6_9();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;

};

class SCREAM_calc_full_delta_asym_VDW_6_9 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_asym_VDW_6_9(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_asym_VDW_6_9();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};


/* 6-8 */
class SCREAM_calc_full_delta_VDW_6_8 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_VDW_6_8(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_VDW_6_8();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

class SCREAM_calc_flat_delta_VDW_6_8 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_VDW_6_8(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_VDW_6_8();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;
};
 
class SCREAM_calc_flat_delta_asym_VDW_6_8 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_asym_VDW_6_8(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_asym_VDW_6_8();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;

};

class SCREAM_calc_full_delta_asym_VDW_6_8 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_asym_VDW_6_8(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_asym_VDW_6_8();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

/* 6-7 */
class SCREAM_calc_full_delta_VDW_6_7 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_VDW_6_7(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_VDW_6_7();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

class SCREAM_calc_flat_delta_VDW_6_7 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_VDW_6_7(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_VDW_6_7();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;
};
 
class SCREAM_calc_flat_delta_asym_VDW_6_7 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_asym_VDW_6_7(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_asym_VDW_6_7();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;

};

class SCREAM_calc_full_delta_asym_VDW_6_7 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_asym_VDW_6_7(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_asym_VDW_6_7();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};
/* end 6-7 */

class SCREAM_calc_flat_delta_asym_X6 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_flat_delta_asym_X6(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_flat_delta_asym_X6();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double delta;

};

class SCREAM_calc_full_delta_asym_X6 : public SCREAM_VDW_BASE_FUNCTIONAL_OBJ {
public:
  SCREAM_calc_full_delta_asym_X6(SCREAM_VDW_OBJ*, double);
  virtual ~SCREAM_calc_full_delta_asym_X6();
  double operator()(SCREAM_ATOM*, SCREAM_ATOM*) const;

  double n_sigma;

};

#endif /* SCREAM_E_FUNCTIONALS_VDW */
