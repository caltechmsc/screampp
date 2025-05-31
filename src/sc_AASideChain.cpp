/* sc_AASideChain.cpp 
 *
 * Source file for classes relevant to AASideChains in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#include "scream_atom.hpp"
#include "Rotlib.hpp"
#include "defs.hpp"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"
#include "scream_tools.hpp"

#include <typeinfo>

AASideChain::AASideChain() {
  SideChain();
  /*
  SC_name = "GEN";
  */
}

AASideChain::AASideChain(const ScreamAtomV& atom_list) : SideChain(atom_list) {

  if (atom_list.size() == 1 and scream_tools::is_HCA_atom(atom_list[0]->getAtomLabel()) ) {
    SC_name = string("GLY");
  } else {
    SC_name = this->sc_atom_mm.begin()->second->resName;
  }

}

AASideChain::AASideChain(const AASideChain& sc) {

  // copy constructor: call operator=: no!  can't do this; memory not initialized.
  // copy constructor: deepcopy, or shallow copy?  shallow copy okay.

  this->sc_atom_mm = sc.get_sc_atom_mm();

}


AASideChain& AASideChain::operator=(const AASideChain& sc) {

  return dummy_assignment_operator(sc);

}


AASideChain::~AASideChain() {

}

/* private functions */

AASideChain& AASideChain::dummy_assignment_operator(const AASideChain& sc) {
  /* Check if it's itself. */
  if (this == &sc) {
    return *this;
  } 
  /* Cleans up allocated storage.  Most likely obsolete. */
  /* SC_name assignment */
  this->SC_name = sc.get_SC_name();

  this->sc_atom_mm = sc.get_sc_atom_mm();

  return *this;
}

/*
double AASideChain::private_chi1(int chi) {

  // no error checking involved!  make sure these atom labels exist.

  if (chi == 1) {

    // working on this!
    
    return calc_dihedral(find_atom_by_greek_name());

  }

  return calc_dihedral();

}
*/

SCREAM_ATOM* AASideChain::get_atom_by_greek_name(string alphabet) {

  multimap<string, SCREAM_ATOM*>::const_iterator mapItr;
  string key;

  if (alphabet == "B") {
    return sc_atom_mm.find("CB")->second;
  } 
  else if (alphabet == "G") {
    for (mapItr = sc_atom_mm.begin(); mapItr != sc_atom_mm.end(); mapItr++) {
      key = scream_tools::strip_whitespace(mapItr->first);
      if (key == "OG" or key == "OG1" or
	  key == "SG" or 
	  key == "CG" or key == "CG1") {
       	return mapItr->second;
      } 
    }
    if (mapItr == sc_atom_mm.end() ) {  // change this to exception handling later.
      return NULL;
    }

  } 
  else if (alphabet == "D") {
    for (mapItr = sc_atom_mm.begin(); mapItr != sc_atom_mm.end(); mapItr++) {
      key = scream_tools::strip_whitespace(mapItr->first);
      if (key == "OD" or key == "OD1" or
	  key == "SD" or
	  key == "ND1" or 
	  key == "CD" or key == "CD1") {
	return mapItr->second;
      } 
    }
    if (mapItr == sc_atom_mm.end() ) {  // change this to exception handling later.
      return NULL;
    }
  } 
  else if (alphabet == "E") {
    for (mapItr = sc_atom_mm.begin(); mapItr != sc_atom_mm.end(); mapItr++) {
      key = scream_tools::strip_whitespace(mapItr->first);
      if (key == "OE" or key == "OE1" or
	  key == "NE" or 
	  key == "CE" or key == "CE1") {
	return mapItr->second;
      }
    }
    if (mapItr == sc_atom_mm.end() ) {  // change this to exception handling later.
      return NULL;
    }
  } 
  else if (alphabet == "Z") {
    for (mapItr = sc_atom_mm.begin(); mapItr != sc_atom_mm.end(); mapItr++) {
      key = scream_tools::strip_whitespace(mapItr->first);
      if (key == "NZ" or
	  key == "CZ") {
	return mapItr->second;
      } 
    }
    if (mapItr == sc_atom_mm.end() ) {  // change this to exception handling later.
      return NULL;
    }
  }
  else if (alphabet == "H") {
    for (mapItr = sc_atom_mm.begin(); mapItr != sc_atom_mm.end(); mapItr++) {
      key = scream_tools::strip_whitespace(mapItr->first);
      if (key == "NH1") {
	return mapItr->second;
      } 
    }
    if (mapItr == sc_atom_mm.end() ) {  // change this to exception handling later.
      return NULL;
    }

  }
  else {
    return NULL;    		// change this to exception handling later.
  }

}
