/** \file sc_Ligand.cpp
 * Ligand class functions defined in this file.
 */

#include "sc_Ligand.hpp"
#include "scream_atom.hpp"
#include "scream_tools.hpp"
#include "sc_ProteinComponent.hpp"
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

Ligand::Ligand() {

  ligand_atoms_on_free_store = false;

}

Ligand::Ligand(const vector<SCREAM_ATOM*>& atom_v)  {

  ligand_atoms_on_free_store = false;
  vector<SCREAM_ATOM*>::const_iterator itr;
  
  for (itr = atom_v.begin(); itr != atom_v.end(); itr++) {
    lig_mm.insert(make_pair(scream_tools::strip_whitespace( (*itr)->atomLabel ), *itr) );
  }

}

Ligand::Ligand(const Ligand& lig)   {

}

Ligand::~Ligand() {

  if (ligand_atoms_on_free_store) {

    multimap<string, SCREAM_ATOM*>::const_iterator itr;
    for (itr = this->lig_mm.begin(); itr != lig_mm.end(); itr++) {
      delete itr->second;
    }
  }

}

vector<SCREAM_ATOM*> Ligand::getAtomList() const {

  vector<SCREAM_ATOM*> returnList;
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = lig_mm.begin(); itr != lig_mm.end(); itr++) {
    returnList.push_back(itr->second);
  }
  return returnList;

}

void Ligand::print_Me() const {

  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  for (itr_mm = lig_mm.begin(); itr_mm != lig_mm.end(); itr_mm++) {
    itr_mm->second->dump();
  }

}

void Ligand::append_to_filehandle(ostream* ofstream_p) const {

  map<int, SCREAM_ATOM*> ordered_m;
  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  for (itr_mm = this->lig_mm.begin(); itr_mm != this->lig_mm.end(); itr_mm++) {
    ordered_m.insert(make_pair(itr_mm->second->n, itr_mm->second));
  }

  map<int, SCREAM_ATOM*>::const_iterator itr_m;
    for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); itr_m++) {
    itr_m->second->append_to_filehandle(ofstream_p);
  }

  
}

void Ligand::append_to_ostream_connect_info(ostream* ofstream_p) const {

  map<int, SCREAM_ATOM*> ordered_m;
  multimap<string, SCREAM_ATOM*>::const_iterator itr_mm;
  for (itr_mm = this->lig_mm.begin(); itr_mm != this->lig_mm.end(); itr_mm++) {
    ordered_m.insert(make_pair(itr_mm->second->n, itr_mm->second));
  }

  map<int, SCREAM_ATOM*>::const_iterator itr_m;
    for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); itr_m++) {
    itr_m->second->append_to_ostream_connect_info(ofstream_p);
  }

  
}


ProteinComponent* Ligand::copy() const {

  // not yet implemented 

  return new Ligand(*this);

}
