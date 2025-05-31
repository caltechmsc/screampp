/* Rotamer.cpp
 *
 * Header file for classes relevant to rotamer in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#include "defs.hpp"

#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"
//#include "ProteinTools.hpp"
#include "Rotamer.hpp"
//b#include "Rotlib.hpp"
#include "sc_BackBone.hpp"
#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"


#include "scream_tools.hpp"

Rotamer::Rotamer() : allocatedScreamAtoms(false), sc(NULL), bb(NULL) {
  // other fields that aren't conceptually crucial
  this->_setDefaults();
}

Rotamer::Rotamer(const ScreamAtomV& atom_list_a, const RotConnInfo* rotConnInfo, bool atomListFromProtein) : allocatedScreamAtoms(false), sc(NULL), bb(NULL) {
  /* atomListFromProtein is set to false by default */

  this->_setDefaults();

  assert(rotConnInfo != NULL);
  this->atom_list.clear();

  /* Create a new list */

  for (ScreamAtomVConstItr itr = atom_list_a.begin(); itr != atom_list_a.end(); ++itr) {
    SCREAM_ATOM* atom = new SCREAM_ATOM();
    atom->copy(*(*itr));
    // connectivities not dealt with here!
    this->atom_list.push_back(atom);
  }

  //  scream_tools::fix_connectivities(this->atom_list, atom_list_a);
  
  if (atomListFromProtein == true) {
    // Then there is a need to do the Rotamer -->Protein atom_n mapping.
    for (ScreamAtomVItr itr = this->atom_list.begin(); itr != atom_list.end(); ++itr) {
      SCREAM_ATOM* atom = *itr;
      map<int, int> inverse_map;
      for (map<int, int>::const_iterator map_int = rotConnInfo->atom_n_map.begin();
	   map_int != rotConnInfo->atom_n_map.end(); ++map_int) {
	inverse_map[map_int->second] = map_int->first;
      }
      atom->n = inverse_map[atom->n]; // atom->n is original Protein atom N; inverse_map gives the rotamer map.
    }
  } else {
    // Then no need to do the Rotamer --> Protein atom_n mapping.     
  }
  
  this->_initSideChain(this->atom_list, rotConnInfo->side_chain_atoms);
  this->_initBackBone(this->atom_list, rotConnInfo->anchor_pts);
  

}

Rotamer::Rotamer(const stringV& rotamerCoords, const RotConnInfo* rotConnInfo) : allocatedScreamAtoms(true), sc(NULL), bb(NULL) {
  is_Original = false;
  same_backbone = false;
  self_E = 99999;

  assert(rotConnInfo != NULL);
  //  cout << "before _initAtomList" << endl;
  this->_initAtomList(rotamerCoords);
  //  cout << "before _initSideChain" << endl;
  this->_initSideChain(this->atom_list, rotConnInfo->side_chain_atoms);
  //  cout << "before _initBackBone" << endl;
  this->_initBackBone(this->atom_list, rotConnInfo->anchor_pts);
  // working on this!


}

Rotamer::~Rotamer() {

  //  cout << "Deleting Rotamer! " << endl;

  if (sc != NULL) {
    delete sc;  			// need to delete sc and bb before ScreamAtomV, i think, because otherwise the SideChain's and the BackBone's would have pointers pointing to nothing before they themselves are deleted.  probably doesn't matter, though.
  }
  if (bb != NULL) {
    delete bb;
  }
  // if Rotamer allocated Scream Atoms, need to free them.
  if (this->allocatedScreamAtoms) {
    ScreamAtomVItr itr;
    for (itr = this->atom_list.begin(); itr != atom_list.end(); ++itr) {
      delete (*itr);
    }
  }

  //  cout << "Done deleting rotamer!" << endl;
}


Rotamer& Rotamer::deepcopy(const Rotamer& rot) {
  /* Deecopies SCREAM_ATOM's. */
  if (this == &rot) {
    return *this;
  }
  
  this->is_Original = rot.is_Original;
  this->same_backbone = rot.same_backbone;

  this->rotamer_n = rot.get_rotamer_n();
  //  this->resName = rot.resName;
  this->mult_H_n = rot.get_mult_H_n();
  this->rotlib_E = rot.get_rotlib_E();
  this->sc_valence_E = rot.get_sc_valence_E();
  this->sc_coulomb_E = rot.get_sc_coulomb_E();
  this->sc_vdw_E = rot.get_sc_vdw_E();
  this->sc_total_nb_E = rot.get_sc_total_nb_E();
  this->sc_solvation_E = rot.get_sc_solvation_E();
  this->sc_total_E = rot.get_sc_total_E();


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
    new_atom->copy(*(*itr));
    this->atom_list.push_back(new_atom);
    old_to_new_atom_map.insert(make_pair(*itr, new_atom));
  }
  this->allocatedScreamAtoms = true;

  // fix connectivities.  for all atoms in this rotamer.  fixing "boundary conditions" (such as N(i) to C(i-1) ) unnecessary.  sidechain placement only affects connectivities that concern atoms immediaitely next to the sidechain portion of a residue.
  for (map<SCREAM_ATOM*, SCREAM_ATOM*>::iterator o2nu = old_to_new_atom_map.begin();
       o2nu != old_to_new_atom_map.end(); ++o2nu) {
    SCREAM_ATOM* old_base_atom = o2nu->first;
    SCREAM_ATOM* new_base_atom = o2nu->second;
    new_base_atom->connectivity_m.clear();



    for (map<SCREAM_ATOM*, int>::const_iterator conn_itr = old_base_atom->connectivity_m.begin();
	 conn_itr != old_base_atom->connectivity_m.end(); ++conn_itr) {
      SCREAM_ATOM* old_connected_atom = conn_itr->first;
      map<SCREAM_ATOM*, SCREAM_ATOM*>::iterator find_it = old_to_new_atom_map.find(old_connected_atom);
      if (find_it != old_to_new_atom_map.end()) {
	
	SCREAM_ATOM* old_connected_atom = find_it->first;

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

  // Figure out which atoms are on the backbone and which atoms are in the sidechain/variable parts. 
  vector<int> bb_atom_n, sc_atom_n;

  assert(rot.get_bb() != NULL);
  assert(rot.get_sc() != NULL);

  ScreamAtomV rot_bb_atoms = rot.get_bb()->get_atoms();
  ScreamAtomV rot_sc_atoms = rot.get_sc()->get_atoms();

  //  cout << "size of rot_sc_atoms " << rot_sc_atoms.size() << endl;

  for (ScreamAtomVConstItr itr = rot_bb_atoms.begin();
       itr != rot_bb_atoms.end(); ++itr) {
    bb_atom_n.push_back( (*itr)->n);
  }
  for (ScreamAtomVConstItr itr  = rot_sc_atoms.begin();
       itr != rot_sc_atoms.end(); ++itr) {
    sc_atom_n.push_back( (*itr)->n);
  }

  this->_initSideChain(atom_list, sc_atom_n);
  this->_initBackBone(atom_list, bb_atom_n);

  this->allocatedScreamAtoms = true; // repeated here from above!  for emphasis.

  return *this;

}

bool Rotamer::read_cnn_lines(const stringV atom_lines) {

}

void Rotamer::printEnergies() const {
  cout << "Rotamer: " << this->rotamer_n << endl;
  cout << "Total: " << this->empty_lattice_E_abs << endl;
  cout << "PreCalc: " << this->preCalc_TotE << endl;
  cout << "VDW: " << this->sc_vdw_E << endl;
  cout << "HB: " << this->sc_hb_E << endl;
  cout << "Coulomb: " << this->sc_coulomb_E << endl;

}


void Rotamer::print_Me() const {
  
  //  cout << "REM " << this->resName << " rotamer " << this->rotamer_n << endl;
  cout << "REM " << "ARB" << " rotamer " << this->rotamer_n << endl;
  cout << "REM " << "energy " << this->rotlib_E << endl;
  this->bb->print_Me();
  this->sc->print_Me();
  cout << "PHI          " << "N/A" << endl;
  cout << "PSI          " << "N/A" << endl;

}

void Rotamer::print_ordered_by_n() const {

  //  cout << "REM " << this->resName << " rotamer " << this->rotamer_n << endl;
  cout << "REM " << "ARB" << " rotamer " << this->rotamer_n << endl;
  cout << "REM energy " << setprecision(7) << this->rotlib_E << " kcal/mol" << endl;
  this->bb->print_ordered_by_n();
  this->sc->print_ordered_by_n();
  //  cout << "PHI          " << this->PHI * 180 / 3.1415926535 << endl;
  //  cout << "PSI          " << this->PSI * 180 / 3.1415926535 << endl;
  // cout <<  "PHI          " << setprecision(15) << this->PHI << endl;
  //  cout <<  "PSI          " << setprecision(15) << this->PSI << endl;
}



SCREAM_ATOM* Rotamer::getAtom(int wantedAtomN) {
  /* Returns SCREAM_ATOM with number n in Rotamer.*/
  SCREAM_ATOM* atomWanted;
  
  /* Find Atom */
  for (ScreamAtomVConstItr itr = this->atom_list.begin(); itr != this->atom_list.end(); ++itr) {
    int loopingAtomN = (*itr)->n;
    if (wantedAtomN == loopingAtomN) {
      atomWanted = (*itr);
      break;
    }
  }


  return atomWanted;

}

ScreamAtomV Rotamer::getTheseAtoms(vector<int> atomNumbers) {
    /* Returns a vector of SCREAM_ATOM*'s from the list provided. */
  ScreamAtomV wantedAtoms;

  for (vector<int>::const_iterator itr = atomNumbers.begin(); itr != atomNumbers.end(); ++itr) {
    SCREAM_ATOM* a = this->getAtom(*itr);
    if (a == NULL) {
      cout << "Warning: Atom not found! in Rotamer::getTheseAtoms." << endl;
    }
    wantedAtoms.push_back(this->getAtom(*itr));
  }
  return wantedAtoms;
}


void Rotamer::fix_toggle(bool value) {

  this->sc->fix_toggle(value);
  this->bb->fix_toggle(value);

}

void Rotamer::fix_sc_toggle(bool value) {

  this->sc->fix_toggle(value);

}

void Rotamer::fix_bb_toggle(bool value) {

  this->bb->fix_toggle(value);

}

void Rotamer::match_bb(const Rotamer*) {

}

void Rotamer::match_CB(const Rotamer* rot, const SCREAM_ATOM* const new_CB, double angle) {

}

int Rotamer::number_of_atoms() const {

  int sc_n = this->sc->number_of_atoms();
  int bb_n = this->bb->number_of_atoms();

  return (sc_n + bb_n);
    

}
double Rotamer::total_charge() const {

  double sc_chg = this->sc->total_charge();
  double bb_chg = this->bb->total_charge();

  return (sc_chg + bb_chg);

}

int Rotamer::sameResidueTypeAs(Rotamer* otherRot) {

  SCREAM_ATOM* thisAtom = this->atom_list[0];
  SCREAM_ATOM* otherAtom = otherRot->getAtom(1);

  if (thisAtom->getResName() == otherAtom->getResName()) {
    return 1;
  } else {
    return 0;
  }
  
}

std::string Rotamer::get_preCalc_Energy_Line() const {

  ostringstream s;
  s << preCalc_TotE << " ";
  s << preCalc_BondsE << " ";
  s << preCalc_AnglesE << " ";
  s << preCalc_TorsionsE << " ";
  s << preCalc_InversionsE << " ";
  s << preCalc_CoulombE << " ";
  s << preCalc_vdwE << " ";
  s << preCalc_HBondE << " ";
  s << preCalc_SolvE;

  return s.str();

}


void Rotamer::populate_preCalc_Terms(std::string l) {
  // line l like:
  //  1.750 0.050 0.253 0.044 0.000 0.000 1.402 0.000 0.000 
  //  TotE  Bonds Angle Torsion Inv Coul  VDW   HB     Solv

  // no bound check for now

  vector<std::string> f;
  split(l, f); // split by " "

  try {

    this->preCalc_TotE = atof(f[0].c_str());
    this->preCalc_BondsE = atof(f[1].c_str());
    this->preCalc_AnglesE = atof(f[2].c_str());
    this->preCalc_TorsionsE = atof(f[3].c_str());
    this->preCalc_InversionsE = atof(f[4].c_str());
    this->preCalc_CoulombE = atof(f[5].c_str());
    this->preCalc_vdwE = atof(f[6].c_str());
    this->preCalc_HBondE = atof(f[7].c_str());
    this->preCalc_SolvE = atof(f[8].c_str());

  }
  
  catch (exception& e) {
    cerr << "Exception caught in Rotamer::populate_preCalc_Terms.  Probably array out of bounds for input energy line string" << endl;
  }


}

/* Protected Function Definitions */


void Rotamer::_setDefaults() {

  this->is_Original = false;
  this->same_backbone = false;
  this->empty_lattice_energy_rank = -1;
  this->self_E = 99999;
  this->rotamer_n = 0;
  this->empty_lattice_E_abs = 99999;
  this->sc_valence_E = 99999;
  this->sc_coulomb_E = 99999;
  this->sc_hb_E = 99999;
  this->sc_vdw_E = 99999;
  this->sc_total_nb_E = 99999;
  this->sc_solvation_E = 99999;
  this->sc_total_E = 99999;

  this->preCalc_TotE = 99999;
  this->preCalc_BondsE = 99999;
  this->preCalc_AnglesE = 99999;
  this->preCalc_TorsionsE = 99999;
  this->preCalc_InversionsE = 99999;
  this->preCalc_CoulombE = 99999;
  this->preCalc_vdwE = 99999;
  this->preCalc_HBondE = 99999;
  this->preCalc_SolvE = 99999;
  


  this->library_name = "";

}

void Rotamer::_initAtomList(const stringV& rotLines) {
  /* Initializes atom_list and nothing else.  Compare AARotamer::initRotamerAtomList */
  stringVConstItr str_itr = rotLines.begin();

  //cout << "in initatomlist" << endl;

  for (; str_itr != rotLines.end(); ++str_itr) {
    string tmp;
    string line = *str_itr;
    if (line.length() <= 1) { // make format errors like an empty line or wrong REMARK lines don't creep in
      continue;
    }
    stringstream ss(line);

    if (line.substr(0,3) == "REM") {
      if (line.substr(6,7) == "rotamer") {
	ss >> tmp;  // REM text
	ss >> this->conformerName;
	ss >> tmp;
	ss >> tmp;  // rotamer number
	this->rotamer_n = atoi(tmp.c_str());
      }
      else if (line.substr(4,6) == "energy") {
	ss >> tmp;
	ss >> tmp;
	ss >> tmp;  // energy value
	if (tmp == string("N/A") or tmp == string("NA") ) {
	  this->rotlib_E = 0;
	} else {
	  std::sscanf(tmp.c_str(), "%lf", &(this->rotlib_E));
	}
      } 
      else if (line.substr(4,3) == "LIB") {
	ss >> tmp;
	ss >> tmp;
	ss >> tmp; // rotamer library name
	this->library_name = tmp;
      }
    } else if (line.substr(0,4) == "ATOM") {
      //cout << "before newing SCREAM_ATOM " << endl;
      SCREAM_ATOM* atom = new SCREAM_ATOM(line);
      atom->library_name = this->library_name;
      this->atom_list.push_back(atom);
      
    } 
  }
  // cout << "exiting _initAtomList" << endl;
}


void Rotamer::_initSideChain(ScreamAtomV& atom_list_a, const vector<int>& sc_int_list_a) {

  assert(atom_list_a.size() != 0);
  assert(sc_int_list_a.size() != 0);
  //cout << "atom_list_a. size  is " << atom_list_a.size() << endl;
  //cout << "sc_int_list_a size is " << sc_int_list_a.size() << endl;
  
  ScreamAtomV sc_atom_list;
  sc_atom_list.clear();

  for (ScreamAtomVItr itr = atom_list_a.begin(); itr != atom_list_a.end(); ++itr) {
    const int atom_n = (*itr)->n;
    /* determine whether this atom is a side chain atom */
    vector<int>::const_iterator pos = find(sc_int_list_a.begin(),
				     sc_int_list_a.end(),
				     atom_n);
    if (pos != sc_int_list_a.end() ){
      sc_atom_list.push_back(*itr); // if is side chain atom, add it to list.
    } else {
      // do nothing; continue;
    }

  } // end looping atom_list_a

  this->sc = new SideChain(sc_atom_list);
  //this->sc->print_Me();
}

void Rotamer::_initBackBone(ScreamAtomV& atom_list_a, const vector<int>& anchor_pts_a) {
  assert(atom_list_a.size() != 0 and anchor_pts_a.size() != 0);
  ScreamAtomV bb_atom_list;
  bb_atom_list.clear();

  for (ScreamAtomVItr itr = atom_list_a.begin(); itr != atom_list_a.end(); ++itr) {
    const int atom_n = (*itr)->n;
    /* determine whether this atom is a bb atom/anchor pt atom */
    vector<int>::const_iterator pos = find(anchor_pts_a.begin(),
					   anchor_pts_a.end(),
					   atom_n);
    if (pos != anchor_pts_a.end()) {
      bb_atom_list.push_back(*itr);
    } else {
      // no nothing; continue
    }
  } // end looping atom_list_a

  this->bb = new BackBone(bb_atom_list);
  //this->bb->print_Me();
}

