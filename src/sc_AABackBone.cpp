/* sc_AABackBone.cpp
 *
 * Source file for classes relevant to AABackBone in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#include "sc_BackBone.hpp"
#include "sc_AABackBone.hpp"
#include "scream_atom.hpp"
#include "Rotlib.hpp"

#include "scream_vector.hpp"
#include "scream_tools.hpp"
using namespace scream_tools;

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace std;


AABackBone::AABackBone() : BackBone() {
  
  init_atom_label_maps();
  PROLINE_flag = false;
}

AABackBone::AABackBone(const SCREAM_ATOM* atom_list) {
  // not implemented
  init_atom_label_maps();
  PROLINE_flag = false;
}

AABackBone::AABackBone(const vector<SCREAM_ATOM*>& atom_list) : BackBone(atom_list) {
  init_atom_label_maps();

  /* Setting PROLINE_flag.
   */
  vector<SCREAM_ATOM*>::const_iterator itr = atom_list.begin();
  if ((*itr)->resName == string("PRO")) {
    PROLINE_flag = true;
  } else {
    PROLINE_flag = false;
  }

}

AABackBone::AABackBone(const AABackBone& bb) {

  //  this->operator=(bb);
  this->PROLINE_flag = bb.PROLINE_flag;
  this->bb_atom_mm = bb.get_bb_atom_mm();
}


AABackBone& AABackBone::operator=(const AABackBone& bb) {

  this->PROLINE_flag = bb.PROLINE_flag;
  if (this == &bb) {
    return *this;
  } else {    // copying member variables; shallow copy.
    this->bb_atom_mm = bb.get_bb_atom_mm(); 
  }

  return *this;
}



AABackBone::~AABackBone() {

}


void AABackBone::assign_atom_fftype() {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); itr++) {
    SCREAM_ATOM* this_atom = itr->second;
    this_atom->atomType = atom_label_fftype_map[scream_tools::strip_whitespace(this_atom->atomLabel)];
  }

}


ScreamVector AABackBone::calc_C_i_minus_one() const {

  // This routine should only be called when a AARotamer object needs a PHI-PSI angle without information from its preceeding residue (i.e. i minu one).

  //  if (this->PROLINE_flag) {                       // if backbone is part of proline, 
    
  //  }

  double N_C_BOND_DIST = 1.32;
  
  SCREAM_ATOM* HN_a;
  HN_a = this->HN();
  if (HN_a == NULL) {
    throw AtomNotFoundException("HN not found in AABackBone::calc_C_i_minus_one() call");
    return 0;
  }

  ScreamVector HN = ScreamVector(HN_a->x[0],
				 HN_a->x[1],
				 HN_a->x[2]);
  ScreamVector N  = ScreamVector(this->N()->x[0],
				 this->N()->x[1],
				 this->N()->x[2]);
  ScreamVector CA = ScreamVector(this->CA()->x[0],
				 this->CA()->x[1],
				 this->CA()->x[2]);

  /*
  ScreamVector HN = ScreamVector(this->bb_atom_mm.find(string("HN"))->second->x[0],
				 this->bb_atom_mm.find(string("HN"))->second->x[1],
				 this->bb_atom_mm.find(string("HN"))->second->x[2]);
  ScreamVector N  = ScreamVector(this->bb_atom_mm.find(string("N"))->second->x[0],
				 this->bb_atom_mm.find(string("N"))->second->x[1],
				 this->bb_atom_mm.find(string("N"))->second->x[2]);
  ScreamVector CA = ScreamVector(this->bb_atom_mm.find(string("CA"))->second->x[0],
				 this->bb_atom_mm.find(string("CA"))->second->x[1],
				 this->bb_atom_mm.find(string("CA"))->second->x[2]);
  */
  ScreamVector N_CA = (CA - N).normalizedVector();
  ScreamVector N_HN = (HN - N).normalizedVector();

  ScreamVector C_i_minus1 = N - (N_CA + N_HN).normalizedVector() * N_C_BOND_DIST;

  return C_i_minus1;

}


ScreamVector AABackBone::calc_N_i_plus_one() const {

  double N_C_BOND_DIST = 1.32;

  SCREAM_ATOM* CA_atom = this->CA();
  SCREAM_ATOM* C_atom  = this->C();
  SCREAM_ATOM* O_atom  = this->O();

  if (CA_atom == NULL or C_atom == NULL or O_atom == NULL) {
    cerr << "Missing atoms in backbone" << endl;
    throw AtomNotFoundException(" Missing either CA, C or O (likely O");
    return 0;
  }

  ScreamVector CA(this->CA());
  ScreamVector C(this->C());
  ScreamVector O(this->O());
  
  ScreamVector C_CA = (CA - C).normalizedVector();
  ScreamVector O_CA = (O - C).normalizedVector();

  ScreamVector N_i_plus1 = C - (C_CA + O_CA).normalizedVector() * N_C_BOND_DIST;

  return N_i_plus1;

}


double AABackBone::PHI() const {

  double phi_angle = 0;

  if (PROLINE_flag) {                                  // if is proline backbone, return zero.  can't calculate C(i-1) position.
    return phi_angle;
  } else {
    SCREAM_ATOM* N_atom = this->N();
    SCREAM_ATOM* CA_atom = this->CA();
    SCREAM_ATOM* C_atom = this->C();

    // throwing exceptions better, eliminate if statements.  learn exceptions.

    if (N_atom == NULL or CA_atom == NULL or C_atom == NULL) {
      cerr << "Missing atom for PHI angle calculation" << endl;
      return 0;
    }

    SCREAM_ATOM* C_i_minus_one_atom = new SCREAM_ATOM();
    ScreamVector v;
    try {
      v = this->calc_C_i_minus_one();
    }
    catch (AtomNotFoundException& notFoundE) {
      cerr << " AtomNotFoundException: " << notFoundE.what << endl;
      return -999.999; // return a phi_angle of -999.999
    }
    C_i_minus_one_atom->x[0] = v[0];
    C_i_minus_one_atom->x[1] = v[1];
    C_i_minus_one_atom->x[2] = v[2];
    
    phi_angle = calc_dihedral(C_i_minus_one_atom, N_atom, CA_atom, C_atom); // returns in degree
    delete C_i_minus_one_atom;

  }
  return phi_angle;
}

double AABackBone::PHI(SCREAM_ATOM* C_i_minus_one) const {
  /*
  cout << "C_i_minus_one coordinates: " << C_i_minus_one->x[0] << " " << C_i_minus_one->x[1] << " " << C_i_minus_one->x[2] << endl;
  cout << "N coords: " << this->N()->x[0] << " " << this->N()->x[1] << " " << this->N()->x[2] << endl;
  cout << "CA coords: " << this->CA()->x[0] << " " << this->CA()->x[1] << " " << this->CA()->x[2] << endl;
  cout << "C coords: " << this->C()->x[0] << " " << this->C()->x[1] << " " << this->C()->x[2] << endl;
  */
  SCREAM_ATOM* N_atom = this->N();
  SCREAM_ATOM* CA_atom = this->CA();
  SCREAM_ATOM* C_atom = this->C();
    

  if (N_atom == NULL or CA_atom == NULL or C_atom == NULL) {
    cerr << "Missing atom for PHI angle calculation" << endl;
    return 0;
  }
  
  return calc_dihedral(C_i_minus_one, this->N(), this->CA(), this->C()); // returns in degree

}

//ScreamMatrix AABackBone::set_PHI(double new_PHI, SCREAM_ATOM* C_i_minus_1 = NULL) { // syntax no logner allowed by gcc 3.0 or above
ScreamMatrix AABackBone::set_PHI(double new_PHI, SCREAM_ATOM* C_i_minus_1) {
  new_PHI = new_PHI * 3.1415926535 / 180;
  double initial_PHI;

  if (C_i_minus_1 == NULL) {
    initial_PHI = this->PHI() * 3.1415926535 / 180;
  } else {
    initial_PHI = this->PHI(C_i_minus_1) * 3.1415926535 / 180;
  }
  while ((new_PHI > 3.1415926535) or (new_PHI <= -3.1415926535)) {
    if (new_PHI > 3.1415926535) {
      new_PHI -= 3.1415926535 * 2;
    } else {
      new_PHI += 3.1415926535 * 2;
    }
  }

  // add 180 to both values so both values are > 0.  otherwise just too confusing for me.

  initial_PHI += 3.1415926535;
  new_PHI += 3.1415926535;

  // now take difference.  this values gives you the angle that you have to rotate right-handed the N-CA bond by.
  // if negative, just add 360.

  double angle_diff = new_PHI - initial_PHI;
  if (angle_diff < 0) {
    angle_diff += 2 * 3.1415926535;
  }

  // Calculate matrix that rotates by angle angle_diff about the N-CA bond.
  
  ScreamMatrix M;  // dummy ScreamMatrix
  ScreamVector N_CA = ScreamVector(this->CA()) - ScreamVector(this->N());

  ScreamMatrix R = M.rotAboutV(N_CA, angle_diff * 180 / 3.1415926535);

  // Now do the transformations for atoms.  In AABAckBone, setting PHI angle would change the positions of C and O.

  ScreamVector origin_offset = ScreamVector(this->N());

  ScreamVector HCA_atom = ScreamVector(this->HCA()) - origin_offset;
  ScreamVector C_atom = ScreamVector(this->C()) - origin_offset;
  ScreamVector O_atom = ScreamVector(this->O()) - origin_offset;

  HCA_atom = R * HCA_atom + origin_offset;
  C_atom = R * C_atom + origin_offset;
  O_atom = R * O_atom + origin_offset;
  
  for (int i = 0; i<=2; i++) {
    this->HCA()->x[i] = HCA_atom[i];
    this->C()->x[i] = C_atom[i];
    this->O()->x[i] = O_atom[i];
  }

  return R;

}

double AABackBone::PSI() const {

  SCREAM_ATOM* N_atom = this->N();
  SCREAM_ATOM* CA_atom = this->CA();
  SCREAM_ATOM* C_atom = this->C();
  
  if (N_atom == NULL or CA_atom == NULL or C_atom == NULL) {
    cerr << "Missing atom for PSI angle calculation" << endl;
    return 0;
  }

  
  try {
    return calc_dihedral(this->N(), this->CA(), this->C(), this->calc_N_i_plus_one()); // returns in degree
  }
  catch (AtomNotFoundException& notFoundE) {
    cerr << " AtomNotFoundException caught: " << notFoundE.what << endl;
    return -999.9999; // PSI if some atom not found.
  }

}

double AABackBone::PSI(SCREAM_ATOM* N_i_plus_one) const {

  SCREAM_ATOM* N_atom = this->N();
  SCREAM_ATOM* CA_atom = this->CA();
  SCREAM_ATOM* C_atom = this->C();
  
  if (N_atom == NULL or CA_atom == NULL or C_atom == NULL) {
    cerr << "Missing atom for PSI angle calculation" << endl;
    return 0;
  }

  return calc_dihedral(this->N(), this->CA(), this->C(), N_i_plus_one); // returns in degree

}

ScreamMatrix AABackBone::set_PSI(double new_PSI, SCREAM_ATOM* N_i_plus_1)  {

  new_PSI = new_PSI * 3.1415926535 / 180;
  double initial_PSI;

  if (N_i_plus_1 == NULL) {
    initial_PSI = this->PSI() * 3.1415926535 / 180;
  } else {
    initial_PSI = this->PSI(N_i_plus_1) * 3.1415926535 / 180;
  }

  while ((new_PSI > 3.1415926535) or (new_PSI <= -3.1415926535)) {
    if (new_PSI > 3.1415926535) {
      new_PSI -= 3.1415926535 * 2;
    } else {
      new_PSI += 3.1415926535 * 2;
    }
  }

  // add 180 to both values so both values are > 0.  otherwise just too confusing for me.

  initial_PSI += 3.1415926535;
  new_PSI += 3.1415926535;

  // now take difference.  this values gives you the angle that you have to rotate right-handed the N-CA bond by.
  // if negative, just add 360.

  double angle_diff = new_PSI - initial_PSI;
  if (angle_diff < 0) {
    angle_diff += 2 * 3.1415926535;
  }
  
  // Calculate matrix that rotates by angle angle_diff about the N-CA bond.
  
  ScreamMatrix M;  // dummy ScreamMatrix
  ScreamVector CA_C = ScreamVector(this->C()) - ScreamVector(this->CA());
  ScreamMatrix R = M.rotAboutV(CA_C, angle_diff * 180 / 3.1415926535);

  // Now do the transformations for atoms.  In AABAckBone, setting PSI angle would change the positions of only O.  Of course, atoms belonging to the next connected backbone mononer would be changed.  That is taken care of in other classes (AAChain), as one mononer should only know about itself and little else.
  ScreamVector origin_offset = ScreamVector(this->CA());

  ScreamVector O_atom = ScreamVector(this->O());

  O_atom = R * (O_atom - origin_offset) + origin_offset;
  
  for (int i = 0; i<=2; i++) {
    this->O()->x[i] = O_atom[i];
  }

  return R;
  
}

/* Specify atoms gets. */

SCREAM_ATOM* AABackBone::N() const {
  
  SCREAM_ATOM* atom = this->get(" N  ");
  if ( atom != NULL) {
    //    cout << "N() works! " << endl;
    return atom;
  } else {
    cerr << "Warning: missing N in res " << this->bb_atom_mm.begin()->second->resNum << endl;
    //  throw range_error("missing N");
    
    return NULL;
  }

}

SCREAM_ATOM* AABackBone::HN() const {
  
  SCREAM_ATOM* atom = this->get("HN  ");
  if ( atom != NULL) {
    //    cout << "N() works! " << endl;
    return atom;
  } else {
    cerr << "Warning: missing HN in res " << this->bb_atom_mm.begin()->second->resNum << endl;
    //    throw range_error("missing HN");
    return NULL;
  }

}

SCREAM_ATOM* AABackBone::CA() const {
  
  SCREAM_ATOM* atom = this->get("CA  ");
  if ( atom != NULL) {
    //    cout << "N() works! " << endl;
    return atom;
  } else {
    cerr << "Warning: missing CA in res " << this->bb_atom_mm.begin()->second->resNum << endl;
    //    throw range_error("missing CA");
    return NULL;
  }

}

SCREAM_ATOM* AABackBone::HCA() const {
  
  SCREAM_ATOM* atom = this->get("HCA ");
  if ( atom != NULL) {
    //    cout << "N() works! " << endl;
    return atom;
  } else {
    cerr << "Warning: missing HCA in res " << this->bb_atom_mm.begin()->second->resNum << endl;
    //    throw range_error("missing HCA");
    return NULL;
  }

}

SCREAM_ATOM* AABackBone::C() const {
  
  SCREAM_ATOM* atom = this->get(" C  ");
  if ( atom != NULL) {
    //    cout << "N() works! " << endl;
    return atom;
  } else {
    cerr << "Warning: missing C in res " << this->bb_atom_mm.begin()->second->resNum << endl;
    //    throw range_error("missing C");
    //    cout << "after throwing range_err" << endl; flush(cout);
    return NULL;
  }

}

SCREAM_ATOM* AABackBone::O() const {
  
  SCREAM_ATOM* atom = this->get(" O  ");
  if ( atom != NULL) {
    //    cout << "N() works! " << endl;
    return atom;
  } else {
    cerr << "Warning: missing O in res " << this->bb_atom_mm.begin()->second->resNum << endl;
    //    throw range_error("missing O");
    return NULL;
  }

}


/******* private functions **********/

map<string, string> AABackBone::atom_label_fftype_map;
map<string, double> AABackBone::atom_label_CHARM22_map;


void AABackBone::init_atom_label_maps() {

  if (atom_label_fftype_map.empty()) {
    atom_label_fftype_map.insert(make_pair(string( "N"  ),  string("N_R"  )));   // dreidii322-quanta.cnv is wrong here.  need more sophisticated conversion tool, because if residue is PRO, N is N_3.
    atom_label_fftype_map.insert(make_pair(string("HN"  ),  string("H___A")));
    atom_label_fftype_map.insert(make_pair(string( "CA" ),  string("C_3"  )));
    atom_label_fftype_map.insert(make_pair(string("HCA" ),  string("H_"   )));
    atom_label_fftype_map.insert(make_pair(string( "C"  ),  string("C_R"  )));   // same problem in dreidii322-quanta.cnv.  if residue following current is PRO, C atom type is C_2.
    atom_label_fftype_map.insert(make_pair(string( "O"  ),  string("O_2"  )));
  }

  if (atom_label_CHARM22_map.empty()) {
    atom_label_CHARM22_map.insert(make_pair(string( "N"  ),  -0.47000));
    atom_label_CHARM22_map.insert(make_pair(string("HN"  ),   0.31000));
    atom_label_CHARM22_map.insert(make_pair(string( "CA" ),   0.07000));
    atom_label_CHARM22_map.insert(make_pair(string("HCA" ),   0.09000));
    atom_label_CHARM22_map.insert(make_pair(string( "C"  ),   0.51000));
    atom_label_CHARM22_map.insert(make_pair(string( "O"  ),  -0.51000));
  }


}
