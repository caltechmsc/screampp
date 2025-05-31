/* sc_AminoAcid.cpp
 *
 * Header file for classes relevant to AminoAcids in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#include "defs.hpp"

#include "sc_AminoAcid.hpp"
#include "sc_BackBone.hpp"
#include "sc_AABackBone.hpp"
#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"
#include "Rotamer.hpp"
#include "AARotamer.hpp"

#include "scream_tools.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"


#include <iostream>

#include <typeinfo>


AminoAcid::AminoAcid() {

  this->rot = new AARotamer();

}

AminoAcid::AminoAcid(int state) {

}


AminoAcid::AminoAcid(const ScreamAtomV& atom_list) {
  Debug debugInfo("AminoAcid::AminoAcid(const ScreamAtomV& atom_list)");

  ScreamAtomV core_atoms_v;  // includes bb_atoms and sc_atoms.

  string this_atom_label; 

  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    core_atoms_v.push_back(*itr);

    if (debugInfo.active())  {
      this_atom_label = scream_tools::strip_whitespace((*itr)->atomLabel);
      if ( scream_tools::is_N_term_hydrogen(this_atom_label) ) {
      } else if (scream_tools::is_C_term_atom(this_atom_label)) {
	debugInfo.out("This chain has a C-Term, as evidenced by the existence of a \"OXT\" atom, or HC atom.");
      }
    }

  }
  
  this->rot = new AARotamer(core_atoms_v);  // need to take care of C-term and N-term extra atoms.

}


AminoAcid::AminoAcid(const AminoAcid&) {

}

AminoAcid::AminoAcid(const AABackBone&, const AASideChain&) {

}

AminoAcid::~AminoAcid() {
  delete rot;
}

AminoAcid & AminoAcid::operator=(const AminoAcid & in_aa) {

  if (this == &in_aa) {
    return *this;
  } else {
    delete this->rot;
    this->rot = new AARotamer(*(in_aa.rot));
    return *this;
  }
  
}


AminoAcid& AminoAcid::deepcopy(const AminoAcid & in_aa) {
  // not implemented yet
  if (this == &in_aa) {
    return *this;
  } 
  delete this->rot;
  ((AARotamer*)this->rot)->deepcopy(*(AARotamer*)(in_aa.rot));

  return *this;

}

SCREAM_ATOM* AminoAcid::operator[](int index) const {

  SCREAM_ATOM* a = new SCREAM_ATOM();
  return a;

}


void AminoAcid::SC_replacement(const AARotamer* const aa_rot, string placementMethod, vector<double>& CreateCBParameters) {
  if (placementMethod == "Default")
    this->rot->match_bb(aa_rot);
  else 
    if (placementMethod == "CreateCB") {
      SCREAM_ATOM* CB = new SCREAM_ATOM;
      this->rot->create_CB(CreateCBParameters, CB);
      this->rot->match_CB(aa_rot, CB, CreateCBParameters[3]); // CreateCBParameters[3] is rotamerMatchLambda

      delete CB;
    }
    else if (placementMethod == "UseExistingCB") {
      SCREAM_ATOM* CB = new SCREAM_ATOM;
      SCREAM_ATOM* this_CB = this->get_CB();
      for (int i=0; i<3; i++) CB->x[i] = this_CB->x[i];
      this->rot->match_CB(aa_rot, CB, CreateCBParameters[3]);

      delete CB;
    }

}

void AminoAcid::assign_atom_fftype() {

  //this->rot->get_bb()->assign_atom_fftype();
  //  this->rot->get_sc()->assign_atom_fftype();  
  this->rot->assign_atom_fftype();

}

void AminoAcid::assign_charges(string scheme, int AA_STATE) {
  
  this->rot->assign_charges(scheme, AA_STATE);

}

void AminoAcid::assign_lone_pair() {

  this->rot->assign_lone_pair();

}

void AminoAcid::translate(ScreamVector V) {

  rot->get_bb()->translate(V);
  rot->get_sc()->translate(V);

}

void AminoAcid::transform(ScreamMatrix M) {

  rot->get_bb()->transform(M);
  rot->get_sc()->transform(M);

}

ScreamAtomV AminoAcid::getAtomList() const {

  ScreamAtomV returnList;
  multimap<string, SCREAM_ATOM*> bb_atoms = this->rot->get_bb()->get_bb_atom_mm();
  multimap<string, SCREAM_ATOM*> sc_atoms = this->rot->get_sc()->get_sc_atom_mm();

  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = bb_atoms.begin(); itr != bb_atoms.end(); ++itr) {
    returnList.push_back(itr->second);
  }
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = sc_atoms.begin(); itr != sc_atoms.end(); ++itr) {
    returnList.push_back(itr->second);
  }

  return returnList;
}

void AminoAcid::fix_toggle(bool value) {
  this->rot->fix_toggle(value);
}


void AminoAcid::fix_sc_toggle(bool value) {
  this->rot->fix_sc_toggle(value);
}

void AminoAcid::fix_bb_toggle(bool value) {
  this->rot->fix_bb_toggle(value);
}

/*
void AminoAcid::SC_replacement(const AARotamer&) {

}
*/
double AminoAcid::PHI() const {
  return ((AABackBone*)this->rot->get_bb())->PHI();
}

double AminoAcid::PHI(const AminoAcid* const aa_prev) const {

  SCREAM_ATOM* C_i_minus_one = ((AABackBone*)aa_prev->rot->get_bb())->C();
  
  return ((AABackBone*)this->rot->get_bb())->PHI(C_i_minus_one);

}

ScreamMatrix AminoAcid::set_PHI(double angle, AminoAcid* prev_aa) {

  
  ScreamMatrix R;
  if (prev_aa == NULL ) {
    R = ((AABackBone*)this->rot->get_bb())->set_PHI(angle);
  } else {
    R = ((AABackBone*)this->rot->get_bb())->set_PHI(angle, ( (AABackBone*)(prev_aa->get_rot()->get_bb()) )->C()  );
  }
  ScreamVector origin_offset = ScreamVector( ( (AABackBone*) (this->rot->get_bb()) )->N());

  this->rot->get_sc()->translate(ScreamVector(0,0,0) - origin_offset);
  this->rot->get_sc()->transform(R);
  this->rot->get_sc()->translate(origin_offset);

  return R;

}

double AminoAcid::PSI() const {

  return ((AABackBone*)this->rot->get_bb())->PSI();

}

double AminoAcid::PSI(const AminoAcid* const aa_next) const {

  SCREAM_ATOM* N_i_plus_one = ((AABackBone*)aa_next->rot->get_bb())->N();
  return ((AABackBone*)this->rot->get_bb())->PSI(N_i_plus_one);

}

ScreamMatrix AminoAcid::set_PSI(double angle, AminoAcid* next_aa) {

  ScreamMatrix R;
  if (next_aa == NULL) {
    R = ((AABackBone*)this->rot->get_bb())->set_PSI(angle);
  } else {
    R = ((AABackBone*)this->rot->get_bb())->set_PSI(angle, ((AABackBone*)(next_aa->get_rot()->get_bb()))->N() );
  }

  ScreamVector origin_offset = ScreamVector( ( (AABackBone*)(this->rot->get_bb()) )->CA());

  return R;

}

double AminoAcid::OMEGA() {

  return 0;

}

ScreamMatrix AminoAcid::set_OMEGA(double angle) {
  return ScreamMatrix();
}

double AminoAcid::chi1() const {

  return this->rot->chi1();

}

double AminoAcid::chi2() const {

  return this->rot->chi2();

}


double AminoAcid::chi3() const {

  return this->rot->chi3();

}


double AminoAcid::chi4() const {

  return this->rot->chi4();

}

double AminoAcid::chi5() const {

  return this->rot->chi5();

}

SCREAM_ATOM* AminoAcid::get_N() const {

  BackBone* bb = this->get_rot()->get_bb();
  SCREAM_ATOM* N = bb->get("N");
  return N;

}

SCREAM_ATOM* AminoAcid::get_CA() const {

  multimap<string, SCREAM_ATOM*> sc_mm = this->rot->get_bb()->get_bb_atom_mm();
  if (sc_mm.find("CA") != sc_mm.end()) {
    return sc_mm.find("CA")->second;
  } else if (sc_mm.find(" CA ") != sc_mm.end()) {
    return sc_mm.find(" CA ")->second;
  } else {
    return NULL;
  }
  
}

SCREAM_ATOM* AminoAcid::get_CB() const {

  multimap<string, SCREAM_ATOM*> sc_mm = this->rot->get_sc()->get_sc_atom_mm();
  if (sc_mm.find("CB") != sc_mm.end()) {
    return sc_mm.find("CB")->second;
  } else if (sc_mm.find(" CB ") != sc_mm.end()) {
    return sc_mm.find(" CB ")->second;
  } else {
    return NULL;
  }
  

}


ScreamAtomV AminoAcid::get_sc_atoms() const {
  //cout << "getting sc_atoms in AminoAcid" << endl;
  return this->rot->get_sc()->get_atoms();

}

ScreamAtomV AminoAcid::get_bb_atoms() const {

  return this->rot->get_bb()->get_atoms();

}

int AminoAcid::number_of_atoms() const {

  //  int nterm_n = nterm_mm.size();
  //  int cterm_n = cterm_mm.size();
  int body_n = this->rot->number_of_atoms();

  //  return (nterm_n + cterm_n + body_n);
  return body_n;

}

double AminoAcid::total_charge() const {

  double total_charge = 0.0;
  //  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = nterm_mm.begin(); itr != nterm_mm.end(); ++itr) {
  //    total_charge += itr->second->q[0];
  //  }
  //  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = cterm_mm.begin(); itr != cterm_mm.end(); ++itr) {
  //    total_charge += itr->second->q[0];
  //  }
  total_charge += this->rot->total_charge();
  
  return total_charge;

}

void AminoAcid::print_Me() const {

  this->rot->get_bb()->print_Me();
  this->rot->get_sc()->print_Me();
  
}

void AminoAcid::print_ordered_by_n() const {

  this->rot->get_bb()->print_ordered_by_n();
  this->rot->get_sc()->print_ordered_by_n();

}
 
void AminoAcid::append_to_filehandle(ostream* ofstream_p) const {

  this->rot->get_bb()->append_to_filehandle(ofstream_p);
  this->rot->get_sc()->append_to_filehandle(ofstream_p);
  

}

void AminoAcid::pdb_append_to_filehandle(ostream* ofstream_p) const {
  this->rot->get_bb()->pdb_append_to_filehandle(ofstream_p);
  this->rot->get_sc()->pdb_append_to_filehandle(ofstream_p);
 
}


void AminoAcid::append_to_ostream_connect_info(ostream* ofstream_p) const {
  this->rot->get_bb()->append_to_ostream_connect_info(ofstream_p);
  this->rot->get_sc()->append_to_ostream_connect_info(ofstream_p);
 
}
