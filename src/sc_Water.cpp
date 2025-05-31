#include "sc_Water.hpp"
#include "sc_ProteinComponent.hpp"
#include "scream_atom.hpp"

#include <vector>
#include <map>
#include <algorithm>

Water::Water() {

  water_atoms_on_free_store = false;

}

Water::Water(const vector<SCREAM_ATOM*>& atom_v) {
  
  water_atoms_on_free_store = false;
  
  vector<SCREAM_ATOM*>::const_iterator itr;
  for (itr = atom_v.begin(); itr != atom_v.end(); ++itr) {
    string atomLabel = (*itr)->atomLabel;
    SCREAM_ATOM* atom = (*itr);

    water_mm.insert(make_pair(atomLabel, atom));


  }

}

Water::~Water() {

  if (water_atoms_on_free_store) {
    for (multimap<string, SCREAM_ATOM*>::iterator itr = water_mm.begin(); itr != water_mm.end(); ++itr) {
      delete itr->second;
    }
  }
}

vector<SCREAM_ATOM*> Water::getAtomList() const {

  vector<SCREAM_ATOM*> returnList;
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = water_mm.begin(); itr != water_mm.end(); ++itr) {
    returnList.push_back(itr->second);
  }
  return returnList;

}

void Water::print_Me() const {

  this->append_to_filehandle(&cout);

}

void Water::append_to_filehandle(ostream* ofstream_p) const {
  
  map<int, SCREAM_ATOM*> ordered_m;
  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  for (itr_mm = this->water_mm.begin(); itr_mm != this->water_mm.end(); ++itr_mm) {
    ordered_m.insert(make_pair(itr_mm->second->n, itr_mm->second));
  }

  map<int, SCREAM_ATOM*>::const_iterator itr_m;
    for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_filehandle(ofstream_p);
  }

}


void Water::append_to_ostream_connect_info(ostream* ofstream_p) const {
  
  map<int, SCREAM_ATOM*> ordered_m;
  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  for (itr_mm = this->water_mm.begin(); itr_mm != this->water_mm.end(); ++itr_mm) {
    ordered_m.insert(make_pair(itr_mm->second->n, itr_mm->second));
  }

  map<int, SCREAM_ATOM*>::const_iterator itr_m;
    for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_ostream_connect_info(ofstream_p);
  }

}




ProteinComponent* Water::copy() const {

  // not yet implemented

  return new Water(*this);

}
