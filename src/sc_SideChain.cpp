/* sc_SideChain.cpp 
 *
 * Source file for classes relevant to SideChains in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#include <cstdlib>

#include "defs.hpp"
#include "scream_atom.hpp"
#include "Rotlib.hpp"
#include <math.h>
#include "scream_tools.hpp"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

#include "sc_SideChain.hpp"

SideChain::SideChain() {
  SC_name = "GEN";
}

SideChain::SideChain(const ScreamAtomV& atom_list_a) {
  
  SC_name = "GEN";
  if (atom_list_a.size() == 0) {
    // do nothing 
  } else {

    ScreamAtomVConstItr itr = atom_list_a.begin();
    for (itr = atom_list_a.begin(); itr != atom_list_a.end(); ++itr) {
      sc_atom_mm.insert(make_pair(scream_tools::strip_whitespace(string( (*itr)->atomLabel )), *itr) );
    }
    SC_name = this->sc_atom_mm.begin()->second->resName;
  }
}

SideChain::~SideChain() {
  
}

void SideChain::fix_toggle(bool value) {

  multimap<string, SCREAM_ATOM*>::iterator itr;
  for (itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    itr->second->fix_atom(value);
  }

}

void SideChain::print_Me() const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    itr->second->dump();
  }

}


void SideChain::print_ordered_by_n() const {
  
  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->dump();
  }
  
}

void SideChain::append_to_filehandle(ostream* ofstream_p) const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->append_to_filehandle(ofstream_p);
  }
  
}

void SideChain::pdb_append_to_filehandle(ostream* ofstream_p) const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->pdb_append_to_filehandle(ofstream_p);
  }
 
}


void SideChain::append_to_ostream_connect_info(ostream* ofstream_p) const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  map<int, SCREAM_ATOM*> order_m;
  for (itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    order_m.insert(make_pair(itr->second->n, itr->second));
  }
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = order_m.begin(); itr_m != order_m.end(); ++itr_m) {
    itr_m->second->append_to_ostream_connect_info(ofstream_p);
  }

}


  
//Void SideChain::translate(ScreamVector& V){
void SideChain::translate(const ScreamVector& V) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = this->sc_atom_mm.begin(); itr != this->sc_atom_mm.end(); ++itr) {
    for (int i = 0; i <= 2; ++i) {
      itr->second->x[i] += V[i];
    }
  }

}

void SideChain::transform(const ScreamMatrix& M) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = this->sc_atom_mm.begin(); itr != this->sc_atom_mm.end(); ++itr) {
    ScreamVector this_atom(itr->second->x[0], itr->second->x[1], itr->second->x[2]);
    ScreamVector transformed_atom((M * this_atom));
    for (int i = 0; i<=2; ++i) {
      itr->second->x[i] = transformed_atom[i];
    }
  }

}

double SideChain::rms(const SideChain* in_sc) const {

  //multimap<string, SCREAM_ATOM*> in_sc_atom_mm = in_sc->get_sc_atom_mm();
  return 0;

}

double SideChain::rms_heavyatom(const SideChain* in_sc) const {

  multimap<string, SCREAM_ATOM*> in_sc_atom_mm = in_sc->get_sc_atom_mm();

  double total_rms = 0;
  int n_of_atoms = 0;
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = in_sc_atom_mm.begin(); itr != in_sc_atom_mm.end(); ++itr) {
    if (scream_tools::is_heavy_atom(itr->first) != true) {
      continue;
    } else {
      SCREAM_ATOM* atom_counterpart = this->sc_atom_mm.find(scream_tools::strip_whitespace(itr->first))->second;
      total_rms += atom_counterpart->distance_squared(itr->second);
      ++n_of_atoms;
    }
  }
  total_rms = sqrt(total_rms);
  return (total_rms/n_of_atoms);


}

void SideChain::copy_atom_positions(const SideChain* new_sc) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  multimap<string, SCREAM_ATOM*>::iterator new_sc_itr;
  multimap<string, SCREAM_ATOM*> new_sc_copy = new_sc->get_sc_atom_mm();
  /*
  // sc_atom_mm 
  for (itr = this->sc_atom_mm.begin(); itr != this->sc_atom_mm.end(); ++itr) {
    cout << string(":::") << itr->first << string(":::") << endl;
  }
  
  // new_sc_copy
  for (itr = new_sc_copy.begin(); itr != new_sc_copy.end(); ++itr) {
    cout << string(":::") << itr->first << string(":::") << endl;
  }
  */
  for (itr = this->sc_atom_mm.begin(); itr != this->sc_atom_mm.end(); ++itr) {
    new_sc_itr = new_sc_copy.find(itr->first);    // finding is well defined because all atomLabel's have been stripped of whitespaces.
    if (new_sc_itr == new_sc_copy.end()) { // didn't find
      break;
    }
    // Copy coordinates.
    itr->second->copyJustCoords(*(new_sc_itr->second));

    // for (int i = 0; i<=2; ++i) {
    //       itr->second->x[i] = new_sc_itr->second->x[i];
    //     }
    new_sc_copy.erase(new_sc_itr);
  }

}


void SideChain::copy_atom_positions_and_library_name(const SideChain* in_sc) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  multimap<string, SCREAM_ATOM*>::iterator in_sc_itr;
  multimap<string, SCREAM_ATOM*> in_sc_copy = in_sc->get_sc_atom_mm();

  for (itr = this->sc_atom_mm.begin(); itr != this->sc_atom_mm.end(); ++itr) {
    string atom_label = itr->first;
    SCREAM_ATOM* self_atom = itr->second;

    in_sc_itr = in_sc_copy.find(atom_label);    // finding is well defined because all atomLabel's have been stripped of whitespaces.

    if (in_sc_itr == in_sc_copy.end()) {
      cout << "SIDECHAIN ATOM NOT FOUND! " << endl;
      cout << " the atom label is: " << atom_label << endl;
      cout << " quitting. " << endl;
      exit(2);
    }

    self_atom->library_name = in_sc_itr->second->library_name; // library name copied.
    for (int i = 0; i<=2; ++i) {
      self_atom->x[i] = in_sc_itr->second->x[i];
    }

    in_sc_copy.erase(in_sc_itr);

  }

}

void SideChain::copy_library_name(const SideChain* in_sc) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  multimap<string, SCREAM_ATOM*>::iterator in_sc_itr;
  multimap<string, SCREAM_ATOM*> in_sc_copy = in_sc->get_sc_atom_mm();

  for (itr = this->sc_atom_mm.begin(); itr != this->sc_atom_mm.end(); ++itr) {
    string atom_label = itr->first;
    SCREAM_ATOM* self_atom = itr->second;

    in_sc_itr = in_sc_copy.find(atom_label);    // finding is well defined because all atomLabel's have been stripped of whitespaces.

    if (in_sc_itr == in_sc_copy.end()) {
      cout << "SIDECHAIN ATOM NOT FOUND! " << endl;
      cout << " the atom label is: " << atom_label << endl;
      cout << " quitting. " << endl;
      exit(2);
    }

    self_atom->library_name = in_sc_itr->second->library_name; // library name copied.
    in_sc_copy.erase(in_sc_itr);

  }
  

}

double SideChain::worst_clash_distance(const vector<SCREAM_ATOM*>& atom_v) const {

  vector<SCREAM_ATOM*>::const_iterator itr_v;
  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  double worst_dist = 999999;
  SCREAM_ATOM *atom, *sc_atom;
  SCREAM_ATOM *worst_clash_atom;
  //  cout << "Size of atom_v passed in: " << atom_v.size() << endl;

  //  cout << "in worst_clash_distance: *************dumping ALL******************" << endl;

  for (itr_v = atom_v.begin(); itr_v != atom_v.end(); ++itr_v) {
    atom = *itr_v;
    //if (atom->resName == "TYR") {
    //atom->dump();
      //    }
    for (itr_mm = sc_atom_mm.begin(); itr_mm != sc_atom_mm.end(); ++itr_mm) {
      sc_atom = itr_mm->second;
      if ( ! (scream_tools::is_heavy_atom(sc_atom->atomLabel) )) {
	continue;
      }
      
      double tmp_dist = sc_atom->distance_squared(atom);
      //cout << tmp_dist << endl;
      if (tmp_dist < worst_dist) {
	//atom->dump();
	//sc_atom->dump();
	worst_clash_atom = atom;
	worst_dist = tmp_dist;
	  
      }
    }
  }
  cout << "Worst clash atom: bgf line:" << endl;
  worst_clash_atom->dump();
  return sqrt(worst_dist);

}

void SideChain::remove_atom(string atomLabel) {
  this->sc_atom_mm.erase(scream_tools::strip_whitespace(atomLabel));
  
}

void SideChain::add_atom(string atomLabel, SCREAM_ATOM* atom) {
  this->sc_atom_mm.insert(make_pair(scream_tools::strip_whitespace(atomLabel), atom) )  ;
}


SCREAM_ATOM* SideChain::get(const string atom_l) const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  string atomLabel = scream_tools::strip_whitespace(atom_l);

  itr = sc_atom_mm.find(atomLabel);
  //cout << "in SideChain::get " << atom_l << endl;
  for (map<string, SCREAM_ATOM*>::const_iterator itr = this->sc_atom_mm.begin(); 
       itr != this->sc_atom_mm.end(); ++itr) {
    //cout << ":::" << itr->first << ":::" << endl;
  }

  if (itr != sc_atom_mm.end()) {
    return itr->second;
  } else {
    itr = sc_atom_mm.find(scream_tools::strip_whitespace(atom_l));   // if input includes whitespace
    if (itr != sc_atom_mm.end()) {
      return itr->second;
    } else if (itr == sc_atom_mm.end()) {
      cerr << "Can't find " << atom_l << " in sc_SideChain::get()" << endl;
      //      exit(8);
      return NULL;
    }
  }


  

}

double SideChain::total_charge() const {

  double total_charge = 0.0;
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    total_charge += itr->second->q[0];
    
  }
  return total_charge;
}

vector<SCREAM_ATOM*> SideChain::get_atoms() const {

  vector<SCREAM_ATOM*> atoms;
  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = sc_atom_mm.begin(); itr != sc_atom_mm.end(); ++itr) {
    atoms.push_back(itr->second);
  }
  return atoms;
}


