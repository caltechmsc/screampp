#include "sc_Hetatm.hpp"
#include "scream_atom.hpp"
#include "scream_tools.hpp"
#include "sc_ProteinComponent.hpp"

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

Hetatm::Hetatm() {

  hetatm_atoms_on_free_store = false;

}

Hetatm::Hetatm(const vector<SCREAM_ATOM*>& atom_v) {

  hetatm_atoms_on_free_store = false;
  vector<SCREAM_ATOM*>::const_iterator itr;
  
  for (itr = atom_v.begin(); itr != atom_v.end(); ++itr) {
    hetatm_mm.insert(make_pair(scream_tools::strip_whitespace( (*itr)->atomLabel ), *itr));
  }

}

Hetatm::Hetatm(const Hetatm& hetatm) {
}

Hetatm::~Hetatm() {

  if (hetatm_atoms_on_free_store) {
  
    multimap<string, SCREAM_ATOM*>::const_iterator itr;
    for (itr = this->hetatm_mm.begin(); itr != hetatm_mm.end(); ++itr) {
      delete itr->second;
    } 
  }

}


vector<SCREAM_ATOM*> Hetatm::getAtomList() const {

  vector<SCREAM_ATOM*> returnList;
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = hetatm_mm.begin(); itr != hetatm_mm.end(); ++itr) {
    returnList.push_back(itr->second);
  }
  return returnList;

}

void Hetatm::print_Me() const {

  this->append_to_filehandle(&cout);

}

void Hetatm::append_to_filehandle(ostream* ofstream_p) const {

  map<int, SCREAM_ATOM*> ordered_m;
  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  for (itr_mm = this->hetatm_mm.begin(); itr_mm != this->hetatm_mm.end(); ++itr_mm) {
    ordered_m.insert(make_pair(itr_mm->second->n, itr_mm->second));
  }

  map<int, SCREAM_ATOM*>::const_iterator itr_m;
    for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_filehandle(ofstream_p);
  }
}

void Hetatm::append_to_ostream_connect_info(ostream* ofstream_p) const {

  map<int, SCREAM_ATOM*> ordered_m;
  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  for (itr_mm = this->hetatm_mm.begin(); itr_mm != this->hetatm_mm.end(); ++itr_mm) {
    ordered_m.insert(make_pair(itr_mm->second->n, itr_mm->second));
  }
  
  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_ostream_connect_info(ofstream_p);
  }
  
}

ProteinComponent* Hetatm::copy() const {

  // not yet implemented

  return new Hetatm(*this);

}
