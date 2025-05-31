/* AARotamer.cpp
 *
 * Source file for classes relevant to rotamer in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */
#include "defs.hpp"

#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"
#include "scream_tools.hpp"

//#include "ProteinTools.hpp"
#include "Rotamer.hpp"
#include "AARotamer.hpp"
//#include "Rotlib.hpp"
#include "sc_BackBone.hpp"
#include "sc_AABackBone.hpp"
#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"

#include <typeinfo>
#include <exception>

#include <cmath>

AARotamer::AARotamer() : Rotamer() {
  // bb, sc initialized to NULL in Rotamer(
  this->resName = string("XXX");
}

AARotamer::AARotamer(const ScreamAtomV& scream_atom_v) {
  /* Note that this constructor is NOT INHERITED from Rotamer(const ScreamAtomV& atom_list)!!! */
  Debug debugInfo("AARotamer::AARotamer(const ScreamAtomV& scream_atom_v)");

  this->_setDefaults();

  ScreamAtomVConstItr itr;
  ScreamAtomV sc_atom_v;
  ScreamAtomV bb_atom_v;

  // copy contents of atom list passed in over to this->atom_list.
  assert(atom_list.size() == 0);
  this->atom_list.insert(atom_list.begin(), scream_atom_v.begin(), scream_atom_v.end());

  // populate sc_atom_v and bb_atom_v
  for (itr = scream_atom_v.begin(); itr != scream_atom_v.end(); ++itr) {
    string this_atom_label = scream_tools::strip_whitespace((*itr)->atomLabel);
    if (scream_tools::is_BB_atom(this_atom_label) ) {
      bb_atom_v.push_back(*itr);
    } else if (scream_tools::is_SC_atom(this_atom_label) ) {
      sc_atom_v.push_back(*itr);
    }
  }

  // deal with Glycine exception.
  if (scream_tools::strip_whitespace(bb_atom_v[0]->resName) == "GLY" ) {
    sc_atom_v.clear();
    sc_atom_v = this->_determine_and_fix_GLY_sidechain_HCA_issue(bb_atom_v);
  }
  
  
  // make new bb and sc.
  debugInfo.out(" Number of atoms in bb_atom_v: " + string(itoa(bb_atom_v.size())) );
  bb = new AABackBone(bb_atom_v);
  debugInfo.out("AABackBone new'ed");

  debugInfo.out(" Number of atoms in sc_atom_v: " + string(itoa(sc_atom_v.size())) );
  sc = new AASideChain(sc_atom_v);
  debugInfo.out("AASideChain new'ed");  


  this->resName = sc->get_SC_name();
  debugInfo.out("AARotamer constructor successful");
}


AARotamer::AARotamer(const vector<string>& string_v) {
  is_Original = false;
  same_backbone = false;
  self_E = 99999;
  this->allocatedScreamAtoms = true;

  // the constructor called by Rotlib constructor when reading in rotamer libraries
  // Format: Rotlib format

  this->initRotamerAtomList(string_v);
  this->PHI = 0;
  this->PSI = 0; // PHI, PSI not used.
  //this->PHI = ((AABackBone*)this->bb)->PHI();  // a) static casting should be fine.  
  //this->PSI = ((AABackBone*)this->bb)->PSI();  // b) PHI, PSI angles calculated using HN and O position respectively since original library don't contain C(i-1) or N(i+1) coordinates..

}

AARotamer::AARotamer(AABackBone* const bb_ps , AASideChain* const sc_ps) {
  is_Original = false;
  same_backbone = false;
  // ps stands for "passed in"
  // Note: an empty GlySC structure could be passed in.

  this->rotamer_n = 0;
  this->resName = string("GEN");

  this->bb = new AABackBone();
  //  this->sc = AASideChain::factory(sc_ps);    // comment: factory returns a pointer to an empty derived AASideChain object.
  this->sc = new AASideChain();          // comment: old code, before 20AA classes introduced.

  (*(AABackBone*)(this->bb)) = (*bb_ps);
  (*(AASideChain*)(this->sc)) = (*sc_ps);    // operator=().
  //obsolete comment.   // if AASideChain is one of 20 derived SC's, operator= overridden in derived SC classes.  here, residue type of sc_ps is guaranteed to be same as this->sc because of the this->sc = AASideChain::factory(sc_ps) line above.

  this->PHI = 0;
  this->PSI = 0;

  //this->PHI = ((AABackBone*)this->bb)->PHI();  // static casting should be fine
  //this->PSI = ((AABackBone*)this->bb)->PSI();

  this->resName = sc->get_SC_name();

}

AARotamer::~AARotamer() {
  // all taken care of in ~Rotamer() 

}


AARotamer& AARotamer::deepcopy(const AARotamer& rot) {
  // Deepcopy: recursively copies SCREAM_ATOM's.  new Backbone and Sidechain objects created.  new SCREAM_ATOMs are allocated.

  if (this == &rot) {
    return *this;
  }

  this->is_Original = rot.is_Original;
  this->same_backbone = rot.same_backbone;

  this->rotamer_n = rot.get_rotamer_n();
  this->resName = rot.resName;
  this->mult_H_n = rot.get_mult_H_n();
  this->rotlib_E = rot.get_rotlib_E();
  this->sc_valence_E = rot.get_sc_valence_E();
  this->sc_coulomb_E = rot.get_sc_coulomb_E();
  this->sc_vdw_E = rot.get_sc_vdw_E();
  this->sc_total_nb_E = rot.get_sc_total_nb_E();
  this->sc_solvation_E = rot.get_sc_solvation_E();
  this->sc_total_E = rot.get_sc_total_E();

  // pre-calc terms

  string preCalc_line = rot.get_preCalc_Energy_Line();
  this->populate_preCalc_Terms(preCalc_line);

  // Now do the real work.  First free up memory.
  if (this->allocatedScreamAtoms) {
    for (ScreamAtomVItr itr = this->atom_list.begin(); itr != atom_list.end(); ++itr)
          delete (*itr);
  }

  // Then clear the original atomlist, if it is occupied.
  this->atom_list.clear();

  // Then copy over entire atom_list from rotamer passed in (rot).
  this->allocatedScreamAtoms = true;
  map<SCREAM_ATOM*, SCREAM_ATOM*> old_to_new_atom_map;

  for (ScreamAtomVConstItr itr = rot.atom_list.begin(); itr != rot.atom_list.end(); ++itr) {
    SCREAM_ATOM* new_atom = new SCREAM_ATOM();
    new_atom->copy(*(*itr)); // shallow copies connectivities
    this->atom_list.push_back(new_atom);
    new_atom->connectivity_m.clear();
    old_to_new_atom_map.insert(make_pair(*itr, new_atom));
  }

  this->allocatedScreamAtoms = true;
  // fix connectivities.  for all atoms in this rotamer (i.e. not just the sidechain.
  for (map<SCREAM_ATOM*, SCREAM_ATOM*>::iterator o2nu = old_to_new_atom_map.begin();
       o2nu != old_to_new_atom_map.end(); ++o2nu) {
    SCREAM_ATOM* old_base_atom = o2nu->first;
    SCREAM_ATOM* new_base_atom = o2nu->second;
    new_base_atom->connectivity_m.clear();

    
    for (map<SCREAM_ATOM*, int>::const_iterator conn_itr = old_base_atom->connectivity_m.begin();
	 conn_itr != old_base_atom->connectivity_m.end(); ++conn_itr) {
      
      SCREAM_ATOM* old_connected_atom = conn_itr->first;
      if (old_connected_atom == NULL) {
	continue;
      }
      
      map<SCREAM_ATOM*, SCREAM_ATOM*>::iterator find_it = old_to_new_atom_map.find(old_connected_atom);
      if (find_it != old_to_new_atom_map.end()) {
	
	// do following only if atoms included in this map.
	SCREAM_ATOM* new_connected_atom = find_it->second;
	int old_order = conn_itr->second;
	new_base_atom->connectivity_m.insert(make_pair(new_connected_atom, old_order));
      } 
      else {
	continue; 
      }
    }

  }

  // Now take care of bb and sc.
  if (this->bb != NULL) delete bb;
  if (this->sc != NULL) delete sc;

  // Populate bb_atom_v and sc_atom_v.
  ScreamAtomV bb_atom_v, sc_atom_v;

  for (ScreamAtomVConstItr itr = this->atom_list.begin(); itr != this->atom_list.end(); ++itr) {
    string this_atom_label = scream_tools::strip_whitespace((*itr)->atomLabel);
    if (scream_tools::is_BB_atom(this_atom_label) ) {
      bb_atom_v.push_back(*itr);
    } else if (scream_tools::is_SC_atom(this_atom_label) ) {
      sc_atom_v.push_back(*itr);
    }
  }
  // deal with GLY exception.
  if (scream_tools::strip_whitespace(bb_atom_v[0]->resName) == "GLY" ) {
    sc_atom_v.clear();
    sc_atom_v = this->_determine_and_fix_GLY_sidechain_HCA_issue(bb_atom_v);
  }

  // make new bb
  this->bb = new AABackBone(bb_atom_v);
  // make new sc
  this->sc = new AASideChain(sc_atom_v);

  this->resName = sc->get_SC_name();

  this->PHI = 0;
  this->PSI = 0;

  //this->PHI = ((AABackBone*)this->bb)->PHI();  // static casting should be fine
  //this->PSI = ((AABackBone*)this->bb)->PSI();

  return *this;  

}


void AARotamer::initRotamerAtomList(const vector<string>& string_v) {

  // Initializes ScreamAtomV atom_list by passing in a vector of string.  Also initiliazes resName.
  
  string tmp, res;
  ScreamAtomV sc_atom_list;
  ScreamAtomV bb_atom_list;
 
  
 
  vector<string>::const_iterator pos;


  for (pos = string_v.begin(); pos != string_v.end(); ++pos) {
    
    string line = *pos;
    stringstream ss(line);

    if (line.substr(0,3) == "REM") {

      if (line.substr(6,7) == "rotamer") {
	ss >> tmp;  // REM text
	ss >> this->resName;
	ss >> tmp;
	ss >> tmp;  // rotamer number
	this->rotamer_n = atoi(tmp.c_str());
      }
      else if (line.substr(4, 6) == "energy") {
	// 3-13-06 new code: energy input format:
	// REM energy 1.750 0.050 0.253 0.044 0.000 0.000 1.402 0.000 0.000 
	// Legend:    TotE  Bonds Angle Torsion Inv Coul  VDW   HB     Solv
	// 

	vector<std::string> f;
	split(line, f); // standard " " as separator.
	if (f.size() == 4) {  
	  // means old style: eg REM energy 9.717384 kcal/mol
	  this->rotlib_E = atof(f[2].c_str());
	}
	else { // assume New All energy terms included style
	  try {
	    this->preCalc_TotE = atof(f[2].c_str());
	    this->preCalc_BondsE = atof(f[3].c_str());
	    this->preCalc_AnglesE = atof(f[4].c_str());
	    this->preCalc_TorsionsE = atof(f[5].c_str());
	    this->preCalc_InversionsE = atof(f[6].c_str());
	    this->preCalc_CoulombE = atof(f[7].c_str());
	    this->preCalc_vdwE = atof(f[8].c_str());
	    this->preCalc_HBondE = atof(f[9].c_str());
	    this->preCalc_SolvE = atof(f[10].c_str());
	  }
	  catch (exception& e) {
	    cerr << "Exception caught: probably out of bounds.  " << e.what() << endl;
	    // keep going.
	  }
	}

      }
      else if (line.substr(4,3) == "lib") {
	stringV f;
	split(line, " ", f);
	this->library_name = f[2];
      }
      
    } 
    else if (line.substr(0, 4) == "ATOM") {
      SCREAM_ATOM* atom = new SCREAM_ATOM(line);
      atom->library_name = this->library_name;
      this->atom_list.push_back(atom);
      
      if (scream_tools::is_BB_atom(scream_tools::atom_label(line))) {
	bb_atom_list.push_back(atom);
      } else {  // to be redundant add (scream_tools::is_SC_atom(screamtools::atom_label(line)))
	sc_atom_list.push_back(atom);
      } 
    } else if (line.substr(0,3) == "PHI") {
      ss >> tmp;
      ss >> tmp;
      
      sscanf(tmp.c_str(), "%lf", &(this->PHI));
    } else if (line.substr(0,3) == "PSI") {
      ss >> tmp;
      ss >> tmp;
      sscanf(tmp.c_str(), "%lf", &(this->PSI));
    }

  }
  
  // deal with GLY exception.
  if (scream_tools::strip_whitespace(bb_atom_list[0]->resName) == "GLY" ) {
    sc_atom_list.clear();
    sc_atom_list = this->_determine_and_fix_GLY_sidechain_HCA_issue(bb_atom_list);
  }
  
  this->bb = new AABackBone(bb_atom_list); 
  this->sc = new AASideChain(sc_atom_list);
  
  this->resName = sc->get_SC_name();
  
}

double AARotamer::calc_PHI() {

  this->PHI = ((AABackBone*)this->bb)->PHI();
  return this->PHI;

}

double AARotamer::calc_PSI() {
  
  this->PSI = ((AABackBone*)this->bb)->PSI();
  return this->PSI;

}

double AARotamer::chi1() {

  return private_chi(1);

}

double AARotamer::chi2() {

  return private_chi(2);

}


double AARotamer::chi3() {

  return private_chi(3);

}


double AARotamer::chi4() {

  return private_chi(4);

}


double AARotamer::chi5() {

  return private_chi(5);

}

void AARotamer::match_bb(const Rotamer* const rot_ps) {

  // for now: assumes Rotamer* is of AARotamer* type.  

  if (rot_ps->is_Original) {  // hack! big time
    cout << "Original rotamer placement." << endl;
    this->sc->copy_atom_positions(rot_ps->get_sc()); 
  } else if (rot_ps->same_backbone) {  // again, hack!
    this->sc->copy_atom_positions(rot_ps->get_sc());  
  } else {
    
    //okay if Rotamer* is passed in.

    // First, match backbone.  Multiple ways to do this.  The method used here aligns the direction in which the N-C bond points.
    ScreamVector this_CAN = ScreamVector( ((AABackBone*)this->bb)->N() ) - ScreamVector( ((AABackBone*)this->bb)->CA() );

    SCREAM_ATOM* N = ((AABackBone*)this->bb)->N();
    SCREAM_ATOM* CA = ((AABackBone*)this->bb)->CA();

    ScreamVector N_vector = ScreamVector(N);
    
    ScreamVector in_CAN = ScreamVector(( (AABackBone*)rot_ps->get_bb())->N()) - ScreamVector(((AABackBone*)rot_ps->get_bb())->CA());
    
    ScreamVector this_CAC = ScreamVector(((AABackBone*)this->bb)->C()) - ScreamVector(((AABackBone*)this->bb)->CA());
    ScreamVector in_CAC = ScreamVector(((AABackBone*)rot_ps->get_bb())->C()) - ScreamVector(((AABackBone*)rot_ps->get_bb())->CA());
    
    ScreamMatrix M; // dummy 
    ScreamMatrix R = M.alignTwoVectors(this_CAN, in_CAN); // aligns two N-C bonds.
    
    ScreamVector this_CAC_Red = this_CAC;
    ScreamVector in_CAC_Red = R * in_CAC;
    
    ScreamVector in_CAN_Red = R * in_CAN;
    
    double dihedral_angle_to_rotate = (this_CAN.cross(this_CAC_Red)).angleBtwn(this_CAN.cross(in_CAC_Red));
    
    if (this_CAN.cross(this_CAC).dot(in_CAC_Red) > 0) {
      dihedral_angle_to_rotate = dihedral_angle_to_rotate;
    } else {
      dihedral_angle_to_rotate = -dihedral_angle_to_rotate;
    }
    
    ScreamMatrix R_about_CAN = M.rotAboutV(ScreamVector(0,0,0) - this_CAN, dihedral_angle_to_rotate);
    
    AASideChain self_sc = *((AASideChain*)this->get_sc());
    self_sc.copy_atom_positions_and_library_name( rot_ps->get_sc() );

    self_sc.translate(ScreamVector(0,0,0) - ScreamVector(((AABackBone*)rot_ps->get_bb())->CA()));
    self_sc.transform(R);
    self_sc.transform(R_about_CAN);
    self_sc.translate(ScreamVector(((AABackBone*)this->bb)->CA()));

  }
  
}

void AARotamer::match_CB(const Rotamer* rot_ps, const SCREAM_ATOM* const CB, double lambda) {
  // remember to put is_Original_flag!!!!
  if (rot_ps->is_Original) {  // hack! big time
    cout << "Original rotamer placement." << endl;
    this->sc->copy_atom_positions(rot_ps->get_sc()); 
  } else if (rot_ps->same_backbone) {  // again, hack!
    this->sc->copy_atom_positions(rot_ps->get_sc());  
  } else {
    
    AASideChain* rot_sc = (AASideChain*)rot_ps->get_sc();
    AABackBone* rot_bb = (AABackBone*)rot_ps->get_bb();
    AASideChain* self_sc = this->get_sc();
    AABackBone* self_bb = this->get_bb();
    // 1. Match CA-CB axis.
    SCREAM_ATOM* CA = self_bb->CA();
    
    ScreamVector CA_v(CA); 
    ScreamVector CB_v(CB); // this CB atom is the new one that's passed in.  It is the one that's been created either through an algorithm or the original CB atom.
    ScreamVector in_CA_v(rot_bb->CA());
    ScreamVector in_CB_v(rot_sc->get_atom_by_greek_name("B"));

    ScreamVector this_CBCA_v = CA_v - CB_v;
    ScreamVector in_CBCA_v = in_CA_v - in_CB_v;

    ScreamMatrix M; // dummy
    ScreamMatrix R = M.alignTwoVectors(this_CBCA_v, in_CBCA_v); // aligns rotamer CB-CA axis to backbone CB-CA axis.
    
    // 2. Match the vector as specified by the lamdba parameter.
    ScreamVector in_C_v(rot_bb->C());
    ScreamVector in_N_v(rot_bb->N());

    ScreamVector in_C_v_Red = in_C_v - in_CB_v;
    ScreamVector in_N_v_Red = in_N_v - in_CB_v;
    ScreamVector in_CA_v_Red = in_CA_v - in_CB_v;

    in_C_v_Red = R * in_C_v_Red;
    in_N_v_Red = R * in_N_v_Red;
    in_CA_v_Red = R * in_CA_v_Red;

    ScreamVector in_CAC_v_Red = in_C_v_Red - in_CA_v_Red;
    ScreamVector in_CAN_v_Red = in_N_v_Red - in_CA_v_Red;
    ScreamVector in_lambda_v_Red = in_CAC_v_Red * (1-lambda) + in_CAN_v_Red * lambda; 
    
    
    ScreamVector this_C_v(self_bb->C());
    ScreamVector this_N_v(self_bb->N());
    
    ScreamVector this_CAC_v = this_C_v - CA_v;
    ScreamVector this_CAN_v = this_N_v - CA_v;

    ScreamVector this_lambda_v = this_CAC_v * (1-lambda) + this_CAN_v * lambda;
    // need to translate to CB as origin coordinate system.  Not that it matters for the subsequent angle calculation.
    
    double dihedral_angle_to_rotate = (this_CBCA_v.cross(this_lambda_v)).angleBtwn(this_CBCA_v.cross(in_lambda_v_Red));
    //cout << "Dihedral_Angle_To_Rotate is: " <<  dihedral_angle_to_rotate  << endl;
    
    //    cout << "Angle between the two CA-C vector is: " << this_lambda_v.angleBtwn(in_lambda_v_Red) << endl;

    if (this_CBCA_v.cross(this_lambda_v).dot(in_lambda_v_Red) > 0) 
      dihedral_angle_to_rotate = -dihedral_angle_to_rotate;
    else
      dihedral_angle_to_rotate = dihedral_angle_to_rotate;

    ScreamMatrix R_about_CBCA = M.rotAboutV(this_CBCA_v, dihedral_angle_to_rotate);
    
    // 3. Do the proper transformations.
    self_sc->copy_atom_positions_and_library_name( rot_sc );

    // Debug stuff:
//     for (int i=0; i<3; i++) {
//       self_bb->N()->x[i] = rot_bb->N()->x[i];
//       self_bb->C()->x[i] = rot_bb->C()->x[i];
//       self_bb->CA()->x[i] = rot_bb->CA()->x[i];

//     }

    self_sc->translate( ScreamVector(0,0,0) - in_CB_v ); // CB position is now pivot.
    self_sc->transform(R);
    self_sc->transform(R_about_CBCA);
    self_sc->translate( CB );  // new CB position is what's wanted on rotamer sidechain.


    // Debug stuff:
//     self_bb->translate( ScreamVector(0,0,0) - in_CB_v);
//     self_bb->transform(R);
//     self_bb->transform(R_about_CBCA);
//     self_bb->translate( CB );


  }
}

void AARotamer::create_CB(vector<double>& CreateCBParameters, SCREAM_ATOM* new_CB) {
  SCREAM_ATOM* C = ((AABackBone*)this->bb)->C();
  SCREAM_ATOM* N = ((AABackBone*)this->bb)->N();
  SCREAM_ATOM* CA = ((AABackBone*)this->bb)->CA();

  ScreamVector CAC_v = ScreamVector(C) - ScreamVector(CA); CAC_v.normalize();
  ScreamVector CAN_v = ScreamVector(N) - ScreamVector(CA); CAN_v.normalize();

  double offBisectorAngle = CreateCBParameters[0];
  double offPlaneAngle = CreateCBParameters[1];
  double bondLength = CreateCBParameters[2];
  double rotamerMatchLambda = CreateCBParameters[3];

  /* Step 1: make C-CA-N plane component of CA-CB vector. After some math, this vector (after normalizing) is:
     CACB_CCAN_projection_flip_v = lambda_A * A_v + lambda_B * B_v
     where: A_v is CAC_v, B_v is CAN_v.  lambda_A and lambda_b are values such that 1) the CACB_CCAN_projection_flip_v has specified angle wrt to CAC_v and CAN_v and 2) CACB_CCAN_projection_flip_v has unit length.
     turns out lambda_A = sin(beta-alpha) / sin(beta)   and lambda_B = sin(alpha)/sin(beta).  where alpha is the angle between CACB_CCAN_projection_flip_v and CAC_v and beta is angle between CAC_v and CAN_v.
     Two ways to derive this: One is by trigonometry.  The other is by vector arithmetic, which is easier although less elegant. */
  
  double beta = acos( CAC_v.dot(CAN_v) ); // * 180 / 3.1415926535;

  double alpha = beta/2 - (offBisectorAngle * 3.1415926535 / 180);

  ScreamVector CACB_CCAN_projection_flip_v = (CAC_v * sin(beta-alpha) + CAN_v * sin(alpha) ) * 1/sin(beta);
  ScreamVector CACB_CCAN_projection_v = ScreamVector(0,0,0) - CACB_CCAN_projection_flip_v;

  /* Step 2: make normal to C-CA-N plane component of this CA-CB vector. */
  ScreamVector CCAN_plane_normal_v = CAN_v.cross(CAC_v);
  ScreamVector CACB_v = CCAN_plane_normal_v * tan(offPlaneAngle * 3.1415926535 /180) + CACB_CCAN_projection_v;
  CACB_v.normalize();

  /* Step 3: bond length adjustment. */
  CACB_v = CACB_v * bondLength;
  
  for (int i=0; i<3; i++)
    new_CB->x[i] = CACB_v[i];

  /* Step 4: add CA coordinate to new_CB coordinates to get the final coordinates. */
  for (int i=0; i<3; i++)
    new_CB->x[i] += CA->x[i];
  
}


void AARotamer::assign_atom_fftype() {

  ((AABackBone*)this->bb)->assign_atom_fftype();
  ((AASideChain*)this->sc)->assign_atom_fftype();

  if (this->resName == string("PRO")) {
    multimap<string, SCREAM_ATOM*> bb_atom_mm = this->bb->get_bb_atom_mm();
    multimap<string, SCREAM_ATOM*>::const_iterator itr;
    
    for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
      SCREAM_ATOM* this_atom = itr->second;
      if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
	this_atom->atomType = string("C_2");         // this is actually wrong, but following mistakes in dreidii322-quanta.cnv.  as mentioned in comments in AABackBone::assign_atom_fftype(), need more sophisticated atom ff type assignment.  
      }
    }
  }
}

void AARotamer::assign_charges(string scheme, int aa_state) {
  
  //((AASideChain*)this->sc)->assign_charges(scheme);
  multimap<string, SCREAM_ATOM*> bb_atom_mm = this->bb->get_bb_atom_mm();
  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  char * SCREAM_NEW_CHG = getenv("SCREAM_NEW_CHG");
  if (SCREAM_NEW_CHG == NULL) {printf("error: enviromental variable SCREAM_NEW_CHG is not set \n"); exit(1);}

  // AMBER or AMBER_N charges: Different backbone charges for each residue
  if (string(SCREAM_NEW_CHG) == "amber" || string(SCREAM_NEW_CHG) == "amber_n") {
      // GLY
      if (this->resName == "GLY") {
          // GLY:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // GLY:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // GLY:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // GLY:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // GLY:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.13595;
                  }
              }

              // GLY:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // GLY:CA:amber/amber_n
                  this_atom->q[0] = -0.02520;
              }

              // GLY:HCA (backbone)
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // GLY:HCA:amber/amber_n
                  this_atom->q[0] = 0.06980;
              }

              // GLY:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // GLY:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // GLY:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // GLY:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
          
          // GLY:sidechain
          multimap<string, SCREAM_ATOM*> sc_atom_mm = this->sc->get_sc_atom_mm();
          for (itr = sc_atom_mm.begin() ; itr != sc_atom_mm.end(); ++itr) {
              // GLY:HCA (sidechain)
              SCREAM_ATOM* this_atom = itr->second;
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // GLY:HCA:amber/amber_n
                      this_atom->q[0] = 0.06980;
              }
          }
      }

      // PRO
      else if (this->resName == string("PRO")) {
          // PRO:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // PRO:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // PRO:N:amber/amber_n
                  this_atom->q[0] = -0.25480;
              }

              // PRO:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // PRO:HN:amber/amber_n/qeq/qeq_n/charmm/charmm_n all == 0
                  this_atom->q[0] =  0.00000;
              }

              // PRO:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // PRO:CA:amber/amber_n
                  this_atom->q[0] = -0.02660;
              }

              // PRO:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // PRO:HCA:amber/amber_n
                  this_atom->q[0] = 0.06410;
              }

              // PRO:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // PRO:C:amber/amber_n
                  this_atom->q[0] = 0.58960;
              }

              // PRO:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // PRO:O:amber/amber_n
                  this_atom->q[0] = -0.57480;
              }
          }
      }


      // ALA ARN ASN APP CYS CYX GLN GLP HIS HSE ILE LEU LYN MET PHE SER THR TRP TYR VAL: amber == amber_n
      // ALA
      else if (this->resName == string("ALA")) {
          // ALA:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // ALA:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // ALA:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // ALA:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // ALA:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // ALA:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // ALA:CA:amber/amber_n
                  this_atom->q[0] = 0.03370;
              }

              // ALA:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // ALA:HCA:amber/amber_n
                  this_atom->q[0] = 0.08230;
              }

              // ALA:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // ALA:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // ALA:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // ALA:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // ARN
      else if (this->resName == string("ARN")) {
          // ARN:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // ARN:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // ARN:N:amber/amber_n
                  this_atom->q[0] = -0.61425;
              }

              // ARN:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // ARN:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.37526;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.37526 / 2.00000;
                  }
              }

              // ARN:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // ARN:CA:amber/amber_n
                  this_atom->q[0] = -0.04834;
              }

              // ARN:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // ARN:HCA:amber/amber_n
                  this_atom->q[0] = 0.09546;
              }

              // ARN:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // ARN:C:amber/amber_n
                  this_atom->q[0] = 0.68923;
              }

              // ARN:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // ARN:O:amber/amber_n
                  this_atom->q[0] = -0.60332;
              }
          }
      }

      // ASN
      else if (this->resName == string("ASN")) {
          // ASN:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // ASN:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // ASN:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // ASN:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // ASN:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // ASN:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // ASN:CA:amber/amber_n
                  this_atom->q[0] = 0.01430;
              }

              // ASN:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // ASN:HCA:amber/amber_n
                  this_atom->q[0] = 0.10480;
              }

              // ASN:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // ASN:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // ASN:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // ASN:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // APP
      else if (this->resName == string("APP")) {
          // APP:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // APP:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // APP:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // APP:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // APP:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // APP:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // APP:CA:amber/amber_n
                  this_atom->q[0] = -0.03410;
              }

              // APP:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // APP:HCA:amber/amber_n
                  this_atom->q[0] = 0.08640;
              }

              // APP:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // APP:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // APP:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // APP:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // CYS
      else if (this->resName == string("CYS")) {
          // CYS:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // CYS:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // CYS:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // CYS:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // CYS:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // CYS:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // CYS:CA:amber/amber_n
                  this_atom->q[0] = 0.02130;
              }

              // CYS:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // CYS:HCA:amber/amber_n
                  this_atom->q[0] = 0.11240;
              }

              // CYS:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // CYS:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // CYS:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // CYS:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // CYX
      else if (this->resName == string("CYX")) {
          // CYX:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // CYX:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // CYX:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // CYX:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // CYX:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // CYX:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // CYX:CA:amber/amber_n
                  this_atom->q[0] = -0.04290;
              }

              // CYX:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // CYX:HCA:amber/amber_n
                  this_atom->q[0] = 0.07660;
              }

              // CYX:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // CYX:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // CYX:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // CYX:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // GLN
      else if (this->resName == string("GLN")) {
          // GLN:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // GLN:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // GLN:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // GLN:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // GLN:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // GLN:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // GLN:CA:amber/amber_n
                  this_atom->q[0] = -0.00310;
              }

              // GLN:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // GLN:HCA:amber/amber_n
                  this_atom->q[0] = 0.08500;
              }

              // GLN:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // GLN:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // GLN:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // GLN:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // GLP
      else if (this->resName == string("GLP")) {
          // GLP:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // GLP:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // GLP:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // GLP:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // GLP:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // GLP:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // GLP:CA:amber/amber_n
                  this_atom->q[0] = -0.01450;
              }

              // GLP:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // GLP:HCA:amber/amber_n
                  this_atom->q[0] = 0.07790;
              }
 
              // GLP:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // GLP:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // GLP:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // GLP:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
         }
      }

      // HIS
      else if (this->resName == string("HIS")) {
          // HIS:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // HIS:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // HIS:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // HIS:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // HIS:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // HIS:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // HIS:CA:amber/amber_n
                  this_atom->q[0] = 0.01880;
              }

              // HIS:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // HIS:HCA:amber/amber_n
                  this_atom->q[0] = 0.08810;
              }

              // HIS:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // HIS:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // HIS:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // HIS:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // HSE
      else if (this->resName == string("HSE")) {
          // HSE:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // HSE:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // HSE:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // HSE:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // HSE:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // HSE:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // HSE:CA:amber/amber_n
                  this_atom->q[0] = -0.05810;
              }

              // HSE:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // HSE:HCA:amber/amber_n
                  this_atom->q[0] = 0.13600;
              }

              // HSE:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // HSE:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // HSE:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // HSE:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // ILE
      else if (this->resName == string("ILE")) {
          // ILE:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // ILE:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // ILE:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // ILE:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // ILE:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // ILE:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // ILE:CA:amber/amber_n
                  this_atom->q[0] = -0.05970;
              }

              // ILE:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // ILE:HCA:amber/amber_n
                  this_atom->q[0] = 0.08690;
              }

              // ILE:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // ILE:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // ILE:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // ILE:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // LEU
      else if (this->resName == string("LEU")) {
          // LEU:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // LEU:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // LEU:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // LEU:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // LEU:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // LEU:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // LEU:CA:amber/amber_n
                  this_atom->q[0] = -0.05180;
              }

              // LEU:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // LEU:HCA:amber/amber_n
                  this_atom->q[0] = 0.09220;
              }

              // LEU:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // LEU:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // LEU:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // LEU:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // LYN
      else if (this->resName == string("LYN")) {
          // LYN:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // LYN:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // LYN:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // LYN:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // LYN:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // LYN:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // LYN:CA:amber/amber_n
                  this_atom->q[0] = -0.07206;
              }

              // LYN:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // LYN:HCA:amber/amber_n
                  this_atom->q[0] = 0.09940;
              }

              // LYN:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // LYN:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // LYN:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // LYN:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // MET
      else if (this->resName == string("MET")) {
          // MET:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // MET:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // MET:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // MET:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // MET:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // MET:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // MET:CA:amber/amber_n
                  this_atom->q[0] = -0.02370;
              }

              // MET:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // MET:HCA:amber/amber_n
                  this_atom->q[0] = 0.08800;
              }

              // MET:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // MET:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // MET:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // MET:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // PHE
      else if (this->resName == string("PHE")) {
          // PHE:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // PHE:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // PHE:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // PHE:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // PHE:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // PHE:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // PHE:CA:amber/amber_n
                  this_atom->q[0] = -0.00240;
              }

              // PHE:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // PHE:HCA:amber/amber_n
                  this_atom->q[0] = 0.09780;
              }

              // PHE:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // PHE:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // PHE:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // PHE:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // SER
      else if (this->resName == string("SER")) {
          // SER:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // SER:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // SER:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // SER:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // SER:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // SER:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // SER:CA:amber/amber_n
                  this_atom->q[0] = -0.02490;
              }

              // SER:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // SER:HCA:amber/amber_n
                  this_atom->q[0] = 0.08430;
              }

              // SER:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // SER:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // SER:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // SER:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // THR
      else if (this->resName == string("THR")) {
          // THR:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // THR:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // THR:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // THR:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // THR:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // THR:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // THR:CA:amber/amber_n
                  this_atom->q[0] = -0.03890;
              }

              // THR:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // THR:HCA:amber/amber_n
                  this_atom->q[0] = 0.10070;
              }

              // THR:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // THR:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // THR:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // THR:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // TRP
      else if (this->resName == string("TRP")) {
          // TRP:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // TRP:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // TRP:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // TRP:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // TRP:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // TRP:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // TRP:CA:amber/amber_n
                  this_atom->q[0] = -0.02750;
              }

              // TRP:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // TRP:HCA:amber/amber_n
                  this_atom->q[0] = 0.11230;
              }

              // TRP:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // TRP:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // TRP:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // TRP:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // TYR
      else if (this->resName == string("TYR")) {
          // TYR:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // TYR:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // TYR:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // TYR:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // TYR:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // TYR:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // TYR:CA:amber/amber_n
                  this_atom->q[0] = -0.00140;
              }

              // TYR:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // TYR:HCA:amber/amber_n
                  this_atom->q[0] = 0.08760;
              }

              // TYR:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // TYR:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // TYR:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // TYR:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // VAL
      else if (this->resName == string("VAL")) {
          // VAL:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // VAL:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // VAL:N:amber/amber_n
                  this_atom->q[0] = -0.41570;
              }

              // VAL:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // VAL:HN:amber/amber_n
                  if (aa_state == 0) {
                      this_atom->q[0] =  0.27190;
                  } else if (aa_state == 1) {
                      this_atom->q[0] =  0.27190 / 2.00000;
                  }
              }

              // VAL:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // VAL:CA:amber/amber_n
                  this_atom->q[0] = -0.08750;
              }

              // VAL:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // VAL:HCA:amber/amber_n
                  this_atom->q[0] = 0.09690;
              }

              // VAL:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // VAL:C:amber/amber_n
                  this_atom->q[0] = 0.59730;
              }

              // VAL:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // VAL:O:amber/amber_n
                  this_atom->q[0] = -0.56790;
              }
          }
      }

      // ARG ASP GLU HSP LYS: amber != amber_n
      // ARG
      else if (this->resName == string("ARG")) {
          // ARG:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // ARG:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // ARG:N:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.34790;
                  }
                  // ARG:N:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.38957;
                  }
              }

              // ARG:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // ARG:HN:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.27470;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.27470 / 2.00000;
                      }
                  }
                  // ARG:HN:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.23303;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.23303 / 2.00000;
                      }
                  }
              }

              // ARG:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // ARG:CA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.26370;
                  }
                  // ARG:CA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.30537;
                  }
              }

              // ARG:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // ARG:HCA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.15600;
                  }
                  // ARG:HCA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.11433;
                  }
              }

              // ARG:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // ARG:C:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.73410;
                  }
                  // ARG:C:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.69243;
                  }
              }

              // ARG:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // ARG:O:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.58940;
                  }
                  // ARG:O:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.63107;
                  }
              }
          }
      }

      // ASP
      else if (this->resName == string("ASP")) {
          // ASP:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // ASP:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // ASP:N:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.51630;
                  }
                  // ASP:N:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.43297;
                  }
              }

              // ASP:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // ASP:HN:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.29360;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.29360 / 2.00000;
                      }
                  }
                  // ASP:HN:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.37693;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.37693 / 2.00000;
                      }
                  }
              }

              // ASP:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // ASP:CA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.03810;
                  }
                  // ASP:CA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.12143;
                  }
              }

              // ASP:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // ASP:HCA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.08800;
                  }
                  // ASP:HCA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.17133;
                  }
              }

              // ASP:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // ASP:C:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.53660;
                  }
                  // ASP:C:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.61993;
                  }
              }

              // ASP:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // ASP:O:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.58190;
                  }
                  // ASP:O:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.49857;
                  }
              }
          }
      }

      // GLU
      else if (this->resName == string("GLU")) {
          // GLU:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // GLU:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // GLU:N:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.51630;
                  }
                  // GLU:N:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.44963;
                  }
              }

              // GLU:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // GLU:HN:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.29360;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.29360 / 2.00000;
                      }
                  }
                  // GLU:HN:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.36027;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.36027 / 2.00000;
                      }
                  }
              }

              // GLU:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // GLU:CA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.03970;
                  }
                  // GLU:CA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.10637;
                  }
              }

              // GLU:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // GLU:HCA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.11050;
                  }
                  // GLU:HCA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.17717;
                  }
              }

              // GLU:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // GLU:C:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.53660;
                  }
                  // GLU:C:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.60327;
                  }
              }

              // GLU:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // GLU:O:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.58190;
                  }
                  // GLU:O:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.51523;
                  }
              }
          }
      }

      // HSP
      else if (this->resName == string("HSP")) {
          // HSP:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // HSP:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // HSP:N:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.34790;
                  }
                  // HSP:N:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.40346;
                  }
              }

              // HSP:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // HSP:HN:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.27470;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.27470 / 2.00000;
                      }
                  }
                  // HSP:HN:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.21914;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.21914 / 2.00000;
                      }
                  }
              }

              // HSP:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // HSP:CA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.13540;
                  }
                  // HSP:CA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.19096;
                  }
              }

              // HSP:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // HSP:HCA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.12120;
                  }
                  // HSP:HCA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.06564;
                  }
              }

              // HSP:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // HSP:C:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.73410;
                  }
                  // HSP:C:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.67854;
                  }
              }

              // HSP:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // HSP:O:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.58940;
                  }
                  // HSP:O:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.64496;
                  }
              }
          }
      }

      // LYS
      else if (this->resName == string("LYS")) {
          // LYS:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // LYS:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // LYS:N:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.34790;
                  }
                  // LYS:N:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.39335;
                  }
              }

              // LYS:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // LYS:HN:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.27470;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.27470 / 2.00000;
                      }
                  }
                  // LYS:HN:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.22925;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.22925 / 2.00000;
                      }
                  }
              }

              // LYS:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // LYS:CA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.24000;
                  }
                  // LYS:CA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.28545;
                  }
              }

              // LYS:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // LYS:HCA:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.14260;
                  }
                  // LYS:HCA:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.09715;
                  }
              }

              // LYS:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // LYS:C:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = 0.73410;
                  }
                  // LYS:C:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = 0.68865;
                  }
              }

              // LYS:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // LYS:O:amber
                  if (string(SCREAM_NEW_CHG) == "amber") {
                      this_atom->q[0] = -0.58940;
                  }
                  // LYS:O:amber_n
                  else if (string(SCREAM_NEW_CHG) == "amber_n") {
                      this_atom->q[0] = -0.63485;
                  }
              }
          }
      }
  }

  // CHARMM or CHARMM_N or QEQ or QEQ_N charges: Same backbone charges
  else {
  
      // GLY
      if (this->resName == "GLY") {
          // GLY:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // GLY:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // GLY:N:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.36910;
                  }
                  // GLY:N:charmm/charmm_n
                  else {
                      this_atom->q[0] = -0.47000;
                  }
              }

              // GLY:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // GLY:HN:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.25090;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.12545;
                      }
                  }
                  // GLY:HN:charmm/charmm_n
                  else {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.31000;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.15500;
                      }
                  }
              }

              // GLY:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // GLY:CA:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.03165;
                  }
                  // GLY:CA:charmm/charmm_n
                  else {
                      this_atom->q[0] = -0.02000;
                  }
              }

              // GLY:HCA (backbone)
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // GLY:HCA:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.11213;
                  }
                  // GLY:HCA:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.09000;
                  }
              }
              
              // GLY:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // GLY:C:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.27455;
                  }
                  // GLY:C:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.51000;
                  }
              }

              // GLY:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // GLY:O:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.34896;
                  }
                  // GLY:O:charmm/charmm_n
                  else {
                      this_atom->q[0] = -0.51000;
                  }
              }
          }

          // GLY:sidechain
          multimap<string, SCREAM_ATOM*> sc_atom_mm = this->sc->get_sc_atom_mm();
          for (itr = sc_atom_mm.begin() ; itr != sc_atom_mm.end(); ++itr) {
              // GLY:HCA (sidechain)
              SCREAM_ATOM* this_atom = itr->second;
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // GLY:HCA:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.11213;
                  }
                  // GLY:HCA:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.09000;
                  }
              }
          }
      }

      // PRO
      else if (this->resName == string("PRO")) {
          // PRO:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // PRO:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // PRO:N:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.29671;
                  }
                  // PRO:N:charmm/charmm_n
                  else {
                      this_atom->q[0] = -0.29000;
                  }
              }

              // PRO:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // PRO:HN:qeq/qeq_n/charmm/charmm_n all == 0
                  this_atom->q[0] =  0.00000;
              }

              // PRO:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // PRO:CA:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.07225;
                  }
                  // PRO:CA:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.02000;
                  }
              }

              // PRO:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // PRO:HCA:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.10840;
                  }
                  // PRO:HCA:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.09000;
                  }
              }

              // PRO:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // PRO:C:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.28178;
                  }
                  // PRO:C:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.51000;
                  }
              }

              // PRO:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // PRO:O:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.34112;
                  }
                  // PRO:O:charmm/charmm_n
                  else {
                      this_atom->q[0] = -0.51000;
                  }
              }
          }
      }

      // ALL OTHERS
      else if (this->resName != string("PRO")) {
          // OTHERS:backbone
          for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
              SCREAM_ATOM* this_atom = itr->second;
              // RES:N
              if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
                  // RES:N:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.36910;
                  }
                  // RES:N:charmm/charmm_n
                  else {
                      this_atom->q[0] = -0.47000;
                  }
              }

              // RES:HN
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
                  // RES:HN:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.25090;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.12545;
                      }
                  }
                  // RES:HN:charmm/charmm_n
                  else {
                      if (aa_state == 0) {
                          this_atom->q[0] =  0.31000;
                      } else if (aa_state == 1) {
                          this_atom->q[0] =  0.15500;
                      }
                  }
              }

              // RES:CA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
                  // RES:CA:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.03165;
                  }
                  // RES:CA:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.07000;
                  }
              }

              // RES:HCA
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
                  // RES:HCA:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.22426;
                  }
                  // RES:HCA:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.09000;
                  }
              }

              // RES:C
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
                  // RES:C:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = 0.27455;
                  }
                  // RES:C:charmm/charmm_n
                  else {
                      this_atom->q[0] = 0.51000;
                  }
              }

              // RES:O
              else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
                  // RES:O:qeq/qeq_n
                  if (string(SCREAM_NEW_CHG) == "qeq" || string(SCREAM_NEW_CHG) == "qeq_n") {
                      this_atom->q[0] = -0.34896;
                  }
                  // RES:O:charmm/charmm_n
                  else {
                      this_atom->q[0] = -0.51000;
                  }
              }
          }
      }
  }
}     
  
void AARotamer::assign_lone_pair() {

  ((AASideChain*)this->sc)->assign_lone_pair();  

  multimap<string, SCREAM_ATOM*> bb_atom_mm = this->bb->get_bb_atom_mm();
  multimap<string, SCREAM_ATOM*>::const_iterator itr;

  if (this->resName != string("PRO")) {

    for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
      SCREAM_ATOM* this_atom = itr->second;
      if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 3;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HN")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 1;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 4;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 1;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 3;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
	this_atom->lone_pair = 2;
	this_atom->atoms_connected = 1;
      } 

    }

  } else if (this->resName == string("PRO")) {

    for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
      SCREAM_ATOM* this_atom = itr->second;
      if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("N")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 3;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("CA")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 4;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("HCA")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 1;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("C")) {
	this_atom->lone_pair = 0;
	this_atom->atoms_connected = 3;
      } else if (scream_tools::strip_whitespace(this_atom->atomLabel) == string("O")) {
	this_atom->lone_pair = 2;
	this_atom->atoms_connected = 1;
      } 
    }
  }
}


ScreamVector AARotamer::calc_C_i_minus_one() {

  double N_C_BOND_DIST = 1.32;
  
  return ((AABackBone*)this->bb)->calc_C_i_minus_one(); // casting safe here--make sure bb is of AABackBone* type.
  

}


void AARotamer::center_CA() {

  SCREAM_ATOM* CA_atom = this->bb->get("CA");
  ScreamVector CA_pstn = ScreamVector(CA_atom->x[0], CA_atom->x[1], CA_atom->x[2]);
  this->bb->translate(ScreamVector(0,0,0)-CA_pstn);
  this->sc->translate(ScreamVector(0,0,0)-CA_pstn);

}

void AARotamer::print_Me() const {


  cout << "REM " << this->resName << " rotamer " << this->rotamer_n << endl;
  cout << "REM " << "energy " << this->rotlib_E << endl;
  this->bb->print_Me();
  this->sc->print_Me();
  cout << "PHI          " << setprecision(15) << this->PHI * 180 / 3.1415926535 << endl;
  cout << "PSI          " << this->PSI * 180 / 3.1415926535 << endl;

}


void AARotamer::print_ordered_by_n() const {
  
  cout << "REM " << this->resName << " rotamer " << this->rotamer_n << endl;
  cout << "REM energy " << setprecision(7) << this->rotlib_E << " kcal/mol" << endl;
  this->bb->print_ordered_by_n();
  this->sc->print_ordered_by_n();
  //  cout << "PHI          " << this->PHI * 180 / 3.1415926535 << endl;
  //  cout << "PSI          " << this->PSI * 180 / 3.1415926535 << endl;
  cout <<  "PHI          " << setprecision(15) << this->PHI << endl;
  cout <<  "PSI          " << setprecision(15) << this->PSI << endl;
}

void AARotamer::append_to_filehandle(ostream* ofstream_p) const {

  *ofstream_p << "REM " << this->resName << " rotamer " << this->rotamer_n << endl;
  *ofstream_p << "REM energy " << setprecision(7) << this->rotlib_E << " kcal/mol" << endl;
  this->bb->append_to_filehandle(ofstream_p);
  this->sc->append_to_filehandle(ofstream_p);
  //  cout << "PHI          " << this->PHI * 180 / 3.1415926535 << endl;
  //  cout << "PSI          " << this->PSI * 180 / 3.1415926535 << endl;
  *ofstream_p <<  "PHI          " << setprecision(15) << this->PHI << endl;
  *ofstream_p <<  "PSI          " << setprecision(15) << this->PSI << endl;

}

void AARotamer::pdb_append_to_filehandle(ostream* ofstream_p) const {

  
  *ofstream_p << "REM " << this->resName << " rotamer " << this->rotamer_n << endl;
  this->bb->pdb_append_to_filehandle(ofstream_p);
  this->sc->pdb_append_to_filehandle(ofstream_p);
  //  cout << "PHI          " << this->PHI * 180 / 3.1415926535 << endl;
  //  cout << "PSI          " << this->PSI * 180 / 3.1415926535 << endl;
  *ofstream_p <<  "PHI          " << setprecision(15) << this->PHI << endl;
  *ofstream_p <<  "PSI          " << setprecision(15) << this->PSI << endl;


}

void AARotamer::append_to_ostream_connect_info(ostream* ofstream_p) const {

  cout << "connectivity info not printed for rotamer" << endl;

}

/* Private Functions */

double AARotamer::private_chi(int n) {

  AASideChain* sc_ptr_t = (AASideChain*)sc;
  string sc_name = sc_ptr_t->get_SC_name();

  SCREAM_ATOM* N = ((AABackBone*)bb)->N();
  SCREAM_ATOM* CA = ((AABackBone*)bb)->CA();
  SCREAM_ATOM* CB = sc_ptr_t->get_atom_by_greek_name("B");

  switch (n) {
  case 1:
    if (sc_name == "GLY" or 
	sc_name == "ALA") {
      return 1000;  // error 
    } else {
      SCREAM_ATOM* XG = sc_ptr_t->get_atom_by_greek_name("G");
      if (XG == NULL) {
	cout << "XG not found" << endl;
	return 1000;
      }
      return scream_tools::calc_dihedral(N, CA, CB, XG);
    }
    break;
  case 2:
    if (sc_name == "ARG" or
	sc_name == "ASN" or 
	sc_name == "ASP" or 
	sc_name == "GLN" or
	sc_name == "GLU" or
	sc_name == "HIS" or    
	sc_name == "HSE" or    
	sc_name == "HSP" or    
	sc_name == "ILE" or 	
	sc_name == "LEU" or
	sc_name == "LYS" or
	sc_name == "MET" or
	sc_name == "PHE" or
	sc_name == "PRO" or
	sc_name == "TRP" or
	sc_name == "TYR") {
      SCREAM_ATOM* XG = sc_ptr_t->get_atom_by_greek_name("G");
      SCREAM_ATOM* XD = sc_ptr_t->get_atom_by_greek_name("D");

      return scream_tools::calc_dihedral(CA, CB, XG, XD);
    } else {
      return 1000;//error
    }
    break;
  case 3:
    if (sc_name == "ARG" or
	sc_name == "GLN" or
	sc_name == "GLU" or
	sc_name == "HSP" or    
	sc_name == "LYS" or
	sc_name == "MET") {
      SCREAM_ATOM* XG = sc_ptr_t->get_atom_by_greek_name("G");
      SCREAM_ATOM* XD = sc_ptr_t->get_atom_by_greek_name("D");
      SCREAM_ATOM* XE = sc_ptr_t->get_atom_by_greek_name("E");

      return scream_tools::calc_dihedral(CB, XG, XD, XE);
    } else {
      return 1000; // error
    }
    break;
  case 4:
    if (sc_name == "ARG" or 
	sc_name == "LYS") {
      SCREAM_ATOM* XG = sc_ptr_t->get_atom_by_greek_name("G");
      SCREAM_ATOM* XD = sc_ptr_t->get_atom_by_greek_name("D");
      SCREAM_ATOM* XE = sc_ptr_t->get_atom_by_greek_name("E");
      SCREAM_ATOM* XZ = sc_ptr_t->get_atom_by_greek_name("Z");

      return scream_tools::calc_dihedral(XG, XD, XE, XZ);
    } else {
      return 1000; // error
    }
    break;
  case 5:
    if (sc_name == "ARG") {
      SCREAM_ATOM* XD = sc_ptr_t->get_atom_by_greek_name("D");
      SCREAM_ATOM* XE = sc_ptr_t->get_atom_by_greek_name("E");
      SCREAM_ATOM* XZ = sc_ptr_t->get_atom_by_greek_name("Z");
      SCREAM_ATOM* XH = sc_ptr_t->get_atom_by_greek_name("H");

      return scream_tools::calc_dihedral(XD, XE, XZ, XH);
    } else {
      return 1000; // error
    }
    break;
  default:
    return 1000;
  }

}


ScreamAtomV AARotamer::_determine_and_fix_GLY_sidechain_HCA_issue(ScreamAtomV& bb_atom_v) {
  /* the L-Amino Acid HCA is removed from bb_atom_v */

  ScreamAtomV sc_atom_v;  // Will only contain the L-AminoAcid HCA.
  
  SCREAM_ATOM *N, *CA, *C, *HCA_1, *HCA_2, *L_HCA;
  int HCA_1_flag = 0;
  // get those atoms.
  for (ScreamAtomVConstItr itr = bb_atom_v.begin(); itr != bb_atom_v.end(); ++itr) {
    string atom_label = scream_tools::strip_whitespace((*itr)->atomLabel);
    if (atom_label == "N") { N = *itr; } 
    else if (atom_label == "CA" ) { CA = *itr; }
    else if (atom_label == "C") { C = *itr;  }
    else if (scream_tools::is_HCA_atom(atom_label) and HCA_1_flag == 0) {  HCA_1 = *itr;       HCA_1_flag = 1;    }
    else if (scream_tools::is_HCA_atom(atom_label) and HCA_1_flag == 1) {      HCA_2 = *itr; }
  }

  // then determine left-right-ness.  use the cross product -- dot product test.

  ScreamVector CA_N = ScreamVector(N) - ScreamVector(CA);
  ScreamVector CA_C = ScreamVector(C) - ScreamVector(CA);
  ScreamVector CA_HCA_1 = ScreamVector(HCA_1) - ScreamVector(CA);  

  if (  (CA_N.cross(CA_C)).dot(CA_HCA_1) > 0 ) {  L_HCA = HCA_1; } 
  else {    L_HCA = HCA_2;   }

  sc_atom_v.push_back(L_HCA);

  ScreamAtomVItr L_HCA_itr = find(bb_atom_v.begin(), bb_atom_v.end(), L_HCA);
  bb_atom_v.erase(L_HCA_itr);

  return sc_atom_v;

}
