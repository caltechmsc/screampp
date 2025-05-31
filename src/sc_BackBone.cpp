/* sc_BackBone.cpp
 *
 * Header file for classes relevant to BackBone in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#include "sc_BackBone.hpp"
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
//using std::string;
using namespace std;

BackBone::BackBone() {
}

BackBone::BackBone(const ScreamAtomV& atom_list_a) {

  vector<SCREAM_ATOM*>::const_iterator itr = atom_list_a.begin();
  for (; itr != atom_list_a.end(); ++itr) {
    bb_atom_mm.insert(make_pair(scream_tools::strip_whitespace(string( (*itr)->atomLabel )), *itr) );
  }

}


void BackBone::print_Me() const {
  
  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  //  cout << this->bb_atom_mm.size() << endl;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    itr->second->dump();
  }
}

void BackBone::fix_toggle(bool value) {

  multimap<string, SCREAM_ATOM*>::iterator itr;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    itr->second->fix_atom(value);
  }

}


void BackBone::print_ordered_by_n() {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->dump();
  }

}

void BackBone::append_to_filehandle(ostream* ofstream_p) const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->append_to_filehandle(ofstream_p);
  }

}

void BackBone::pdb_append_to_filehandle(ostream* ofstream_p) const {


  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->pdb_append_to_filehandle(ofstream_p);
  }

}


void BackBone::append_to_ostream_connect_info(ostream* ofstream_p) const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->append_to_ostream_connect_info(ofstream_p);
  }

}

void BackBone::translate(const ScreamVector& V) {
  
  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = this->bb_atom_mm.begin(); itr != this->bb_atom_mm.end(); ++itr) {
    for (int i = 0; i <= 2; ++i) {
      itr->second->x[i] += V[i];
    }
  }

}


void BackBone::transform(const ScreamMatrix& M) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = this->bb_atom_mm.begin(); itr != this->bb_atom_mm.end(); ++itr) {
    ScreamVector this_atom(itr->second->x[0], itr->second->x[1], itr->second->x[2]);
    ScreamVector transformed_atom((M * this_atom));
    for (int i = 0; i<=2; ++i) {
      itr->second->x[i] = transformed_atom[i];
    }

  }
}


SCREAM_ATOM* BackBone::get(const string atom_l) const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  /*
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    itr->second->dump();
  }
  */
  string atomLabel = scream_tools::strip_whitespace(atom_l);
  itr = bb_atom_mm.find(atomLabel);
  
  if (itr != bb_atom_mm.end()) {
    return itr->second;
  } else {
    cerr << "Can't find atom " << atomLabel << " in BackBone::get(const string atom_l)" << endl;
    return NULL;
  }

}


double BackBone::total_charge() const {

  double chg = 0.0;
  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    chg += itr->second->q[0];
  }
  return chg;
}

vector<SCREAM_ATOM*> BackBone::get_atoms() const {

  vector<SCREAM_ATOM*> atoms;
  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = bb_atom_mm.begin(); itr != bb_atom_mm.end(); ++itr) {
    atoms.push_back(itr->second);
  }
  return atoms;

}

void BackBone::copy_atom_positions(const BackBone* in_bb) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  multimap<string, SCREAM_ATOM*>::iterator in_bb_itr;
  multimap<string, SCREAM_ATOM*> in_bb_copy = in_bb->get_bb_atom_mm();

  for (itr = this->bb_atom_mm.begin(); itr != this->bb_atom_mm.end(); ++itr) {
    in_bb_itr = in_bb_copy.find(itr->first);
    if (in_bb_itr == in_bb_copy.end()) {
      cerr << " Can't find atom " << itr->first << " in BackBone::copy_atom_positions(const BackBone* in_bb)! " << endl;
      cerr << " Info: in_bb->print_Me() " << endl;
      in_bb->print_Me();
    }
    for (int i = 0; i<=2; ++i) {
      itr->second->x[i] = in_bb_itr->second->x[i];
    }
    in_bb_copy.erase(in_bb_itr);
  }

}
