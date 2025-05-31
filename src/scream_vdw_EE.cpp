#include "defs.hpp"
#include "MutInfo.hpp"
#include "scream_vdw_EE.hpp"
#include <algorithm>

#include <time.h>
#include <stdio.h>
#include <cassert>

#include "ClashCollection.hpp"
#include "RotamerNeighborList.hpp"

using namespace std;

VDW_EE::VDW_EE() {
  this->clashCollection = NULL;
  this->rotamerNeighborList = NULL;
}

VDW_EE::VDW_EE(Protein* ptn, vector<MutInfo> mutInfo_V, SCREAM_VDW_OBJ* vdw_obj_) : vdw_obj(vdw_obj_) {

  this->ptn = ptn;

  this->_initVariableAndFixedAtomPairList(ptn, mutInfo_V);
  this->_initVariableAndVariableAtomPairList(ptn, mutInfo_V);

  this->clashCollection = NULL; // Instantiation.Initialization.
  this->ON_THE_FLY = 1;
}

VDW_EE::~VDW_EE() {

  

}

void VDW_EE::init_after_addedMutInfoRotConnInfo(Protein* ptn, SCREAM_VDW_OBJ* vdw_obj_) {
  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->vdw_obj = vdw_obj_;

  this->_initVariableAndFixedAtomPairListArb(ptn, this->mutInfo_rotConnInfo_map);
  this->_initVariableAndVariableAtomPairListArb(ptn, this->mutInfo_rotConnInfo_map);

  this->ptn = ptn;
  this->ON_THE_FLY = 0;

}

void VDW_EE::init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein* ptn, SCREAM_VDW_OBJ* vdw_obj_) {
  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->vdw_obj = vdw_obj_;

  //  this->_initFixedMoveableAtomsOnProtein(ptn, this->mutInfo_rotConnInfo_map);
  this->ptn = ptn;
  this->ON_THE_FLY = 1;

}

void VDW_EE::init_after_addedMutInfoRotConnInfo_neighbor_list(Protein* ptn, SCREAM_VDW_OBJ* vdw_obj_, RotamerNeighborList* rNL) {
  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->vdw_obj = vdw_obj_;
  
  this->ptn = ptn;
  this->ON_THE_FLY = 2;
  this->rotamerNeighborList = rNL;

}

void VDW_EE::addMutInfoRotConnInfo(MutInfo mutInfo, RotConnInfo* rotConnInfo = NULL) {
  /* adds a mutInfo, rotconninfo pair to mutInfo_rotConnInfo_map */
  this->mutInfo_rotConnInfo_map[mutInfo] = rotConnInfo;
}


void VDW_EE::addClashCollection(ClashCollection* cc) {

  this->clashCollection = cc;

}

void VDW_EE::cleanClashCollection() {
  //delete this->clashCollection; // should delete be here?  no... ClashCollection usually instantiated in upper level.
  this->clashCollection = NULL;
  
}


double VDW_EE::calc_empty_lattice_E(const MutInfo& mutInfo) {
  
  double total_E = 0;
  vector<AtomPair> atomPairList(this->variable_and_fixed[mutInfo]);
  
  //clock_t t1 = clock();
  if (! this->ON_THE_FLY) {
    for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
      total_E += this->vdw_obj->calc_VDW_6_12( (*ap_i).a1, (*ap_i).a2 );
    }
  }
  else {

    // first get atoms from this sidechain.
    ScreamAtomV sc_atoms = ptn->get_sc_atoms(mutInfo.chn, mutInfo.pstn);

    // then get all atoms on the protein.
    ScreamAtomV atom_list = ptn->getAtomList();

    // main double loop
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
	int ptn_flag = (*ptn_atom)->flags & 0x2; // Second least significant bit determines whether the atom is moveable or fixed with respect to this energy expression.
	if (ptn_flag == 0) { // i.e. part of sidechain/variable atoms
	  if (*sc_atom == *ptn_atom) continue;  // if same atom continue
	  if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {  // if two atoms on same sidechain
	    if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
		   or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) )
	      total_E += this->vdw_obj->calc_VDW_6_12( *ptn_atom, *sc_atom);   // then only calc energies if 1-2 and 1-3 exclusions are satisfied.
	  }
	  else 
	    continue; // if the two atoms are not on the same residue, don't calculate energies.
	}
	
	else {
	  total_E += this->vdw_obj->calc_VDW_6_12( *ptn_atom, *sc_atom); // calculate energies in al other instances.

	}
      
      }
    }
  }
  //  clock_t t2 = clock();
  //  cout << "that took: " << (double)(t2 - t1) / (double) CLOCKS_PER_SEC << endl;

  return total_E;

}

double VDW_EE::calc_empty_lattice_E_delta(const MutInfo& mutInfo, string mode, double r) {

  double total_E = 0;
  int  i = 0;
  int  double_count = 0;

  vector<AtomPair> atomPairList(this->variable_and_fixed[mutInfo]);
  //cout << " DEBUG: this->ON_THE_FLY value: " << this->ON_THE_FLY << endl;
  //  cout << "VDW EE calc size: " << atomPairList.size() << endl;
  //clock_t t1 = clock();  
  if (! this->ON_THE_FLY) {
    if (mode == "FULL") {
      for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
	total_E += this->vdw_obj->calc_full_delta_VDW_6_12( (*ap_i).a1, (*ap_i).a2, r );
      }
    }
    else if (mode == "FLAT") {
      for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
	total_E += this->vdw_obj->calc_flat_delta_VDW_6_12( (*ap_i).a1, (*ap_i).a2, r );
      }
    }
    else if (mode == "SCALED") {
      for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
	total_E += this->vdw_obj->calc_VDW_6_12_scaled_inner_wall( (*ap_i).a1, (*ap_i).a2, r );
      }
    }
  }
  else { // on the fly
    if ((this->ON_THE_FLY == 1)  or (this->ON_THE_FLY == 2) ) {
      total_E += this->_calc_empty_lattice_E_on_the_fly_loop(mutInfo, mode, r);
    }
  } // end on the fly calc loop

  return total_E;
}


// double VDW_EE::calc_empty_lattice_E_delta_asym(const MutInfo& mutInfo, string mode, double r) {
//   double total_E = 0;
//   vector<AtomPair> atomPairList(this->variable_and_fixed[mutInfo]);
//   //  cout << "VDW EE calc size: " << atomPairList.size() << endl;
//   //clock_t t1 = clock();  
//   if (! this->ON_THE_FLY) {
//     if (mode == "FULL") {
//       for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
// 	total_E += this->vdw_obj->calc_full_delta_asym_VDW_6_12( (*ap_i).a1, (*ap_i).a2, r );
//       }
//     }
//     else if (mode == "FLAT") {
//       for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
// 	total_E += this->vdw_obj->calc_flat_delta_asym_VDW_6_12( (*ap_i).a1, (*ap_i).a2, r );
//       }
//     }
//     else if (mode == "SCALED") {
//       for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
// 	total_E += this->vdw_obj->calc_VDW_6_12_scaled_inner_wall( (*ap_i).a1, (*ap_i).a2, r );
//       }
//     }
//     else if (mode == "RESIDUE") {
//       for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
// 	total_E += this->vdw_obj->calc_residue_delta_asym_VDW_6_12( (*ap_i).a1, (*ap_i).a2 );
//       }
//     }
//   }
//   else if ( this->ON_THE_FLY) {

//   }
//   //  clock_t t2 = clock();
//   //  cout << "that took: " << (double)(t2 - t1) / (double) CLOCKS_PER_SEC << endl;
//   return total_E;


// }

double VDW_EE::calc_residue_interaction_E(const MutInfo& mI) {
  // calculates interaction between MutInfo residue and rest of the variable atoms.
  double total_E = 0;

  if (! this->ON_THE_FLY) {

    for (map< MutInfoPair, vector<AtomPair> >::iterator itr = this->variable_and_variable.begin();
	 itr != this->variable_and_variable.end(); ++itr) {
      if ( mI == (itr->first).mutInfo1 or mI == (itr->first).mutInfo2 ) {
	for (vector<AtomPair>::iterator ap_i = (itr->second).begin(); ap_i != (itr->second).end(); ++ap_i) {
	  total_E += this->vdw_obj->calc_VDW_6_12( (*ap_i).a1, (*ap_i).a2 );
	}
      }
    }
  }

  else if (this->ON_THE_FLY) {

  }

  return total_E;
}


double VDW_EE::calc_residue_interaction_E(const MutInfo& mI1, const MutInfo& mI2, string method, double scale) {

  double interaction_E = 0;
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor;
  if (method == "FULL") {
    vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (method == "FLAT") {
    vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (method == "SCALE") {
    // not implemented yet.
    cout << "SCALE functor not implemented yet." << endl;
  } else {
    cout << " functor " << method  << " not implemented yet." << endl;
  }

  ScreamAtomV mI1_atoms, mI2_atoms; mI1_atoms.clear(); mI2_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI1_rCI_itr = this->mutInfo_rotConnInfo_map.find(mI1);
  map<MutInfo, RotConnInfo*>::const_iterator mI2_rCI_itr = this->mutInfo_rotConnInfo_map.find(mI2);

  if (mI1_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI1_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    mI1_atoms = this->ptn->get_sc_atoms(mI1);  
  }
  else {
    mI1_atoms = this->ptn->get_variable_atoms(mI1_rCI_itr->second);
  }
  
  if (mI2_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI2_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    mI2_atoms = this->ptn->get_sc_atoms(mI2);  
  }
  else {
    mI2_atoms = this->ptn->get_variable_atoms(mI2_rCI_itr->second);
  }
  
  // Actual calculations of energies.

  double total_E = 0;

  for (ScreamAtomVConstItr mI1_atom_itr = mI1_atoms.begin(); mI1_atom_itr != mI1_atoms.end(); ++mI1_atom_itr) {
    for (ScreamAtomVConstItr mI2_atom_itr = mI2_atoms.begin(); mI2_atom_itr != mI2_atoms.end(); ++mI2_atom_itr) {

      SCREAM_ATOM* mI1_atom = *mI1_atom_itr;
      SCREAM_ATOM* mI2_atom = *mI2_atom_itr;

      if (! (scream_tools::should_exclude_on_1_2(mI1_atom, mI2_atom)
	       or scream_tools::should_exclude_on_1_3(mI1_atom, mI2_atom) ) ) {
	double E = (*vdw_functor)(mI1_atom, mI2_atom);
	total_E += E;
	
      }
    }
  }

  delete vdw_functor;
  return total_E;

}

double VDW_EE::calc_all_interaction_E() {
  // calculates interaction energies between all atom pairs
  double total_E = 0;

  if (! this->ON_THE_FLY) {
    
    for (map< MutInfoPair, vector<AtomPair> >::iterator itr = this->variable_and_variable.begin(); 
	 itr != this->variable_and_variable.end(); ++itr) {
      for (vector<AtomPair>::iterator ap_i = (itr->second).begin(); ap_i != (itr->second).end(); ++ap_i) {
	total_E += this->vdw_obj->calc_VDW_6_12( (*ap_i).a1, (*ap_i).a2 );
	
      }
    }

  }

  else if (this->ON_THE_FLY) {
    total_E += this->_calc_all_interaction_E_on_the_fly_loop("FLAT", 0);
  }

  return total_E;
}

double VDW_EE::calc_all_interaction_E_delta(std::string mode, double r) {


  double total_E = 0;

  if (! this->ON_THE_FLY) {
    
    for (map< MutInfoPair, vector<AtomPair> >::iterator itr = this->variable_and_variable.begin(); 
	 itr != this->variable_and_variable.end(); ++itr) {
      for (vector<AtomPair>::iterator ap_i = (itr->second).begin(); ap_i != (itr->second).end(); ++ap_i) {
	// could optimize inner loop
	if (mode == "FULL") {
	  total_E += this->vdw_obj->calc_full_delta_VDW_6_12( (*ap_i).a1, (*ap_i).a2 , r);
	}
	if (mode == "FLAT") {
	  total_E += this->vdw_obj->calc_flat_delta_VDW_6_12((*ap_i).a1, (*ap_i).a2, r);
	}
	if (mode == "SCALED") {
	  total_E += this->vdw_obj->calc_VDW_6_12_scaled_inner_wall((*ap_i).a1, (*ap_i).a2, r);
	}
	
      }
    }
  }
  
  else if (this->ON_THE_FLY) {
    total_E += this->_calc_all_interaction_E_on_the_fly_loop(mode, r);
  }

  return total_E;

}

double VDW_EE::calc_EL_rot_selfBB(const MutInfo& mutInfo, std::string method, double scale) {
  Debug debugInfo("VDW_EE::calc_EL_rot_selfBB");
  
  /* method: FULL, FLAT, SCALE, FULL_NONONPOLARH, FLAT_NONONPOLARH */
  int i = 0;
  int double_count = 0;

  string base_method;
  stringV f;
  split(method, "_", f);
  base_method = f[0];

  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  if (f.size() == 2) {
    if (f[1] == "ASYM") {      base_method = base_method + "_ASYM";    }
    else if (f[1] == "NONONPOLARH" ) {      nonPolar_H_calc = 0;    }
  }

  if (f.size() == 3) {
    base_method += "_ASYM";
    nonPolar_H_calc = 0;
  }

  //debugInfo.out("Method used: " + method);
  //debugInfo.out("Base Method: " + base_method);
  //debugInfo.out("Nonpolar H calcalations? 1 = yes, 0 = no: " + string(itoa(nonPolar_H_calc)) );

  double total_E = 0;
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor, *no_scale_functor;
  if (base_method == "FULL") {
    vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "FLAT") {
    vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "FULL_ASYM") {
    vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "FLAT_ASYM") {
    vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
  } else {
    cout << " functor " << base_method << " not implemeneted yet." << endl;
    exit(2);
  }

  no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);

  // first get atoms from this sidechain.
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);

  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }

  // then get all relevant atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  atom_list = this->ptn->getAtomList();

  double dist_sq_cutoff = 10.0*10.0;
  int count_tot_E = 0;
  // Remake the two atom lists for nonPolar_H_calc.
  if (nonPolar_H_calc == 0) {
    ScreamAtomV noNonPolarH_atom_list; noNonPolarH_atom_list.clear();
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      if ( "H_" == (*ptn_atom)->stripped_atomType ) {
	continue;
      }
      else {
	noNonPolarH_atom_list.push_back(*ptn_atom);
      }
    }
    ScreamAtomV noNonPolarH_sc_atom_list; noNonPolarH_sc_atom_list.clear();
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      if ( "H_" == (*sc_atom)->stripped_atomType ) {
	continue;      }
      else {
	noNonPolarH_sc_atom_list.push_back(*sc_atom);      
      }
    }
    atom_list.clear();
    sc_atoms.clear();
    atom_list.insert(atom_list.end(), noNonPolarH_atom_list.begin(), noNonPolarH_atom_list.end());
    sc_atoms.insert(sc_atoms.end(), noNonPolarH_sc_atom_list.begin(), noNonPolarH_sc_atom_list.end());
  }
  //debugInfo.out(" Size of atom_list: " + string(itoa(atom_list.size()) ) );
  //debugInfo.out(" Size of sc_atoms: " + string(itoa(sc_atoms.size())) );


  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  /* For this particular one: calc_EL_rot_selfBB */
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }
  if (mutInfoChn_is_Z) {
    return 0; // If mutInfoChn is Z, selfBB energy undefined, since ligand 
  }
  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 1: moveable atoms.  point for moveable and immoveable atoms if for optimization.
      if (ptn_flag == 0) { // ptn_flag == 0: moveable atoms, i.e. those atoms on sidechains to be replaced.
	continue; // now that there's no need for internal SC-SC VDW E calculation, skip.
      }
      // Case 2: fixed atoms.
      else { // if ptn_flag ==1 
	// same residue atoms, including backbone atoms.
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  string scAtomLabel = (*sc_atom)->stripped_atomLabel;
	  string ptnAtomLabel = (*ptn_atom)->stripped_atomLabel;
	  // Note: the following lines need to change so that these atom names won't be hardcoded.
	  if ( ptnAtomLabel != "O" and
	       ptnAtomLabel != "HN" and
	       ptnAtomLabel != "OXT") {
	    // NEW SCHEME: SC-SC VDW not longer included in calculations, provided by rotamer library. 5-03-06.  HN, O, OXT's interactions with sidechain NOT provided by library because they are backbone dependent, therefore needs to be calulcated by this program.  Coulomb and HBond's still calculated because those terms were excluded in the calculations.
	    continue;
	  }
	  // following conditional: unnecessary, because sc atoms always away by at least 3 bonds from O, HN, or OXT.
	  if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
		 or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	    //debugInfo.out(" ptn_atom->stripped_atomLabel is: " + ptnAtomLabel);
	    //debugInfo.out(" sc_atom ->stripped_atomLabel is: " + scAtomLabel);

	    //double E = (*vdw_functor)(*ptn_atom, *sc_atom); // no need for /2.  not gonna bother with spline.
	    double E = (*no_scale_functor)(*ptn_atom, *sc_atom); // if same backbone, don't use flat-bottom potential function.  Or maybe i should, don't know.
	    total_E += E;

	  } // end should not exclude block
	} // end same residue block
      } // end ptn_flag == 1 block (fixed atoms)
      
      
    } // end inner loop
  } // end outer loop
  
  return total_E;

}

double VDW_EE::calc_EL_rot_otherBB(const MutInfo& mutInfo, std::string method, double scale) {
  Debug debugInfo("VDW_EE::calc_EL_rot_otherBB");
  
  /* method: FULL, FLAT, SCALE, FULL_NONONPOLARH, FLAT_NONONPOLARH */
  int i = 0;
  int double_count = 0;

  string base_method;
  stringV f;
  split(method, "_", f);
  base_method = f[0];

  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  if (f.size() == 2) {
    if (f[1] == "ASYM") {
      base_method = base_method + "_ASYM";
    }
    else if (f[1] == "NONONPOLARH" ) {
      nonPolar_H_calc = 0;
    }
  }

  if (f.size() == 3) {
    base_method += "_ASYM";
    nonPolar_H_calc = 0;
  }

  //debugInfo.out("Method used: " + method);
  //debugInfo.out("Base Method: " + base_method);
  //debugInfo.out("Nonpolar H calcalations? 1 = yes, 0 = no: " + string(itoa(nonPolar_H_calc)) );

  double total_E = 0;
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor, *no_scale_functor;
  if (base_method == "FULL") {
    vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "FLAT") {
    vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "SCALE") {
    // not implemented yet.
  } else if (base_method == "FULL_ASYM") {
    vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_full_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "FLAT_ASYM") {
    vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "SCALE_ASYM") {
    // no implemented yet
    cout << " This functor, SCALE_ASYM, not implemented yet." << endl;
    exit(2);
  } else {
    cout << " functor " << base_method << " not implemeneted yet." << endl;
    exit(2);
  }

  no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);
  //  no_scale_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, 0);

  // first get atoms from this sidechain.
  
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);


  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }

  // then get all relevant atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  //  if (this->ON_THE_FLY == 1) {
  atom_list = this->ptn->getAtomList();

  double dist_sq_cutoff = 10.0*10.0;
  
  int count_tot_E = 0;
  // Remake the two atom lists for nonPolar_H_calc.
  if (nonPolar_H_calc == 0) {
    ScreamAtomV noNonPolarH_atom_list; noNonPolarH_atom_list.clear();
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      //      if ( (*ptn_atom)->atomType.c_str()[0] != 'H') {
      //	noNonPolarH_atom_list.push_back(*ptn_atom);	
      //	continue;
      //      } // fast filter; string comparison too slow
      if ( "H_" == (*ptn_atom)->stripped_atomType ) {
	continue;
      }
      else {
	noNonPolarH_atom_list.push_back(*ptn_atom);
      }
    }
    ScreamAtomV noNonPolarH_sc_atom_list; noNonPolarH_sc_atom_list.clear();
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      //if ( (*sc_atom)->atomType.c_str()[0] != 'H') {
      //noNonPolarH_sc_atom_list.push_back(*sc_atom);
      //continue;
      //}
      if ( "H_" == (*sc_atom)->stripped_atomType ) {
	continue;      }
      else {
	noNonPolarH_sc_atom_list.push_back(*sc_atom);      
      }
    }
    atom_list.clear();
    sc_atoms.clear();
    atom_list.insert(atom_list.end(), noNonPolarH_atom_list.begin(), noNonPolarH_atom_list.end());
    sc_atoms.insert(sc_atoms.end(), noNonPolarH_sc_atom_list.begin(), noNonPolarH_sc_atom_list.end());
  }
  //debugInfo.out(" Size of atom_list: " + string(itoa(atom_list.size()) ) );
  //debugInfo.out(" Size of sc_atoms: " + string(itoa(sc_atoms.size())) );


  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }
  if (mutInfoChn_is_Z) {
    return 0;
  }

  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {

      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      // Case 1: moveable atoms.  point for moveable and immoveable atoms if for optimization.
      if (ptn_flag == 0) { // ptn_flag == 0: moveable atoms, i.e. those atoms on sidechains to be replaced.
	continue; // now that there's no need for internal SC-SC VDW E calculation, skip.
      }
      // Case 2: fixed atoms.
      else { // if ptn_flag ==1 
	// same residue atoms, including backbone atoms--don't do these
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  continue;
	} // end same residue block
	// if is a side chain atom, don't do the calculation.
	if ( scream_tools::is_SC_atom((*ptn_atom)->getAtomLabel()) ) {
	  continue;
	}
	else { // do calcultion only if it on backbone.
	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);
	  if (dist_sq > dist_sq_cutoff) {
	    continue;
	  } 
	  else {
	    double E = (*vdw_functor)(*ptn_atom, *sc_atom);  // no division by 2
	    total_E += E;
	  }
	}
      } // end ptn_flag == 1 block (fixed atoms)
    } // end inner loop
  } // end outer loop
  
  return total_E;

}

double VDW_EE::calc_EL_rot_fixedSC(const MutInfo& mutInfo, std::string method, double scale) {
  Debug debugInfo("VDW_EE::calc_EL_rot_fixedSC");
  
  /* method: FULL, FLAT, SCALE, FULL_NONONPOLARH, FLAT_NONONPOLARH */
  int i = 0;
  int double_count = 0;

  // Jan 16: if ASYM, use X6, the more accurate one.  If SYMM, the old one, X6 not written, hence LJ 6-12.
  // Mar 04: back to LJ.  Why did I even bother?  X6 much slower, need to take care of inflection point else E --> -inf.

  string base_method;
  stringV f;
  split(method, "_", f);
  base_method = f[0];

  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  if (f.size() == 2) {
    if (f[1] == "ASYM") {
      base_method = base_method + "_ASYM";
    }
    else if (f[1] == "NONONPOLARH" ) {
      nonPolar_H_calc = 0;
    }
  }

  if (f.size() == 3) {
    base_method += "_ASYM";
    nonPolar_H_calc = 0;
  }

  //debugInfo.out("Method used: " + method);
  //debugInfo.out("Base Method: " + base_method);
  //debugInfo.out("Nonpolar H calcalations? 1 = yes, 0 = no: " + string(itoa(nonPolar_H_calc)) );

  double total_E = 0;
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor, *no_scale_functor;
  if (base_method == "FULL") {
    vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "FLAT") {
    vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "SCALE") {
    // not implemented yet.
  } else if (base_method == "FULL_ASYM") {
    vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_full_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "FLAT_ASYM") {
    vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "SCALE_ASYM") {
    // no implemented yet
    cout << " This functor, SCALE_ASYM, not implemented yet." << endl;
    exit(2);
  } else {
    cout << " functor " << base_method << " not implemeneted yet." << endl;
    exit(2);
  }

  no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);
  //  no_scale_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, 0);

  // first get atoms from this sidechain.
  
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);


  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }

  // then get all relevant atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  //  if (this->ON_THE_FLY == 1) {
  atom_list = this->ptn->getAtomList();
  double dist_sq_cutoff = 10.0*10.0;
  int count_tot_E = 0;
  // Remake the two atom lists for nonPolar_H_calc.
  if (nonPolar_H_calc == 0) {
    ScreamAtomV noNonPolarH_atom_list; noNonPolarH_atom_list.clear();
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      if ( "H_" == (*ptn_atom)->stripped_atomType ) {
	continue;
      }
      else {
	noNonPolarH_atom_list.push_back(*ptn_atom);
      }
    }
    ScreamAtomV noNonPolarH_sc_atom_list; noNonPolarH_sc_atom_list.clear();
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      if ( "H_" == (*sc_atom)->stripped_atomType ) {
	continue;      }
      else {
	noNonPolarH_sc_atom_list.push_back(*sc_atom);      
      }
    }
    atom_list.clear();
    sc_atoms.clear();
    atom_list.insert(atom_list.end(), noNonPolarH_atom_list.begin(), noNonPolarH_atom_list.end());
    sc_atoms.insert(sc_atoms.end(), noNonPolarH_sc_atom_list.begin(), noNonPolarH_sc_atom_list.end());
  }
  //debugInfo.out(" Size of atom_list: " + string(itoa(atom_list.size()) ) );
  //debugInfo.out(" Size of sc_atoms: " + string(itoa(sc_atoms.size())) );


  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }
  if (mutInfoChn_is_Z) {
    return 0;
  }

  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {

      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      if (mutInfoChn_is_Z) {
	return 0;
      }
      // Case 1: moveable atoms.  point for moveable and immoveable atoms if for optimization.
      else if (ptn_flag == 0) { // ptn_flag == 0: moveable atoms, i.e. those atoms on sidechains to be replaced.
	continue; // now that there's no need for internal SC-SC VDW E calculation, skip.
      }
      // Case 2: fixed atoms.
      else { // if ptn_flag ==1 
	// same residue atoms, including backbone atoms.
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  continue;
	} // end same residue block
	if ( scream_tools::is_BB_atom((*ptn_atom)->getAtomLabel()) ) {
	  continue;
	}
	else { // iteration with other firxed SC's.
	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);
	  if (dist_sq > dist_sq_cutoff) {
	    continue;
	  } 
	  else {
	    double E = (*vdw_functor)(*ptn_atom, *sc_atom);  // no division by 2
	    total_E += E;
	  }
	}
      } // end ptn_flag == 1 block (fixed atoms)
    } // end inner loop
  } // end outer loop
  
  return total_E;

}

double VDW_EE::calc_EL_rot_fixedHET(const MutInfo& mutInfo, std::string method, double scale) {
  Debug debugInfo("VDW_EE::_calc_empty_lattice_E_on_the_fly_loop");
  
  /* method: FULL, FLAT, SCALE, FULL_NONONPOLARH, FLAT_NONONPOLARH */
  int i = 0;
  int double_count = 0;

  // Jan 16: if ASYM, use X6, the more accurate one.  If SYMM, the old one, X6 not written, hence LJ 6-12.
  // Mar 04: back to LJ.  Why did I even bother?  X6 much slower, need to take care of inflection point else E --> -inf.

  string base_method;
  stringV f;
  split(method, "_", f);
  base_method = f[0];

  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  if (f.size() == 2) {
    if (f[1] == "ASYM") {
      base_method = base_method + "_ASYM";
    }
    else if (f[1] == "NONONPOLARH" ) {
      nonPolar_H_calc = 0;
    }
  }

  if (f.size() == 3) {
    base_method += "_ASYM";
    nonPolar_H_calc = 0;
  }

  //debugInfo.out("Method used: " + method);
  //debugInfo.out("Base Method: " + base_method);
  //debugInfo.out("Nonpolar H calcalations? 1 = yes, 0 = no: " + string(itoa(nonPolar_H_calc)) );

  double total_E = 0;
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor, *no_scale_functor;
  if (base_method == "FULL") {
    vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "FLAT") {
    vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "SCALE") {
    // not implemented yet.
  } else if (base_method == "FULL_ASYM") {
    vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_full_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "FLAT_ASYM") {
    vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "SCALE_ASYM") {
    // no implemented yet
    cout << " This functor, SCALE_ASYM, not implemented yet." << endl;
    exit(2);
  } else {
    cout << " functor " << base_method << " not implemeneted yet." << endl;
    exit(2);
  }

  no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);
  //  no_scale_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, 0);

  // first get atoms from this sidechain.
  
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);


  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }

  

  // then get all relevant atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  //  if (this->ON_THE_FLY == 1) {
  atom_list = this->ptn->getAtomList();

    //}
//   else if (this->ON_THE_FLY == 2) {
//     atom_list = this->rotamerNeighborList->returnEmptyLatticeNeighborList(mutInfo);
//   }
  
  double dist_sq_cutoff = 10.0*10.0;

  
  int count_tot_E = 0;
  // Remake the two atom lists for nonPolar_H_calc.
  if (nonPolar_H_calc == 0) {
    ScreamAtomV noNonPolarH_atom_list; noNonPolarH_atom_list.clear();
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      //      if ( (*ptn_atom)->atomType.c_str()[0] != 'H') {
      //	noNonPolarH_atom_list.push_back(*ptn_atom);	
      //	continue;
      //      } // fast filter; string comparison too slow
      if ( "H_" == (*ptn_atom)->stripped_atomType ) {
	continue;
      }
      else {
	noNonPolarH_atom_list.push_back(*ptn_atom);
      }
    }
    ScreamAtomV noNonPolarH_sc_atom_list; noNonPolarH_sc_atom_list.clear();
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      //if ( (*sc_atom)->atomType.c_str()[0] != 'H') {
      //noNonPolarH_sc_atom_list.push_back(*sc_atom);
      //continue;
      //}
      if ( "H_" == (*sc_atom)->stripped_atomType ) {
	continue;      }
      else {
	noNonPolarH_sc_atom_list.push_back(*sc_atom);      
      }
    }
    atom_list.clear();
    sc_atoms.clear();
    atom_list.insert(atom_list.end(), noNonPolarH_atom_list.begin(), noNonPolarH_atom_list.end());
    sc_atoms.insert(sc_atoms.end(), noNonPolarH_sc_atom_list.begin(), noNonPolarH_sc_atom_list.end());
  }
  //debugInfo.out(" Size of atom_list: " + string(itoa(atom_list.size()) ) );
  //debugInfo.out(" Size of sc_atoms: " + string(itoa(sc_atoms.size())) );


  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }

  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {

      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      //if (mutInfo.chn == "Z") {
      if (mutInfoChn_is_Z) {
	if (*sc_atom == *ptn_atom) continue;
	if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
	       or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	  //if (this->ON_THE_FLY == 1) {
	  //total_E += (*vdw_functor)(*ptn_atom, *sc_atom)/2;
	  //count_tot_E++;
	  //} else if (this->ON_THE_FLY == 2) {
	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);
	  if (dist_sq > dist_sq_cutoff) {
	    continue;
	  } else {
	    double E = (*vdw_functor)(*ptn_atom, *sc_atom)/2;
	    total_E += E;

	    
	  }
	}
      }
      // Case 1: moveable atoms.  point for moveable and immoveable atoms if for optimization.
      else if (ptn_flag == 0) { // ptn_flag == 0: moveable atoms, i.e. those atoms on sidechains to be replaced.
	continue; // now that there's no need for internal SC-SC VDW E calculation, skip.
      }
      // Case 2: fixed atoms.
      else { // if ptn_flag ==1 
	// same residue atoms, including backbone atoms.
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  string scAtomLabel = (*sc_atom)->stripped_atomLabel;
	  string ptnAtomLabel = (*ptn_atom)->stripped_atomLabel;
	  // Note: the following lines need to change so that these atom names won't be hardcoded.
	  if ( ptnAtomLabel != "O" and
	       ptnAtomLabel != "HN" and
	       ptnAtomLabel != "OXT") {
	    // NEW SCHEME: SC-SC VDW not longer included in calculations, provided by rotamer library. 5-03-06.  HN, O, OXT's interactions with sidechain NOT provided by library because they are backbone dependent, therefore needs to be calulcated by this program.  Coulomb and HBond's still calculated because those terms were excluded in the calculations.
	    continue;
	  }
	  // following conditional: unnecessary, because sc atoms always away by at least 3 bonds from O, HN, or OXT.
	  if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
		 or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	    // CA <=> CG interactions already excluded by the previous condition.
	    //if (this->ON_THE_FLY == 1) {
	    //total_E += (*vdw_functor)(*ptn_atom, *sc_atom)/2;
	    //} else if (this->ON_THE_FLY == 2) {

	    
	    //double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);  // no longer checking distance; if on same residue, always calculate interaction
	    //if (dist_sq > dist_sq_cutoff) {
	    //continue;
	    //} else {
	    //debugInfo.out(" ptn_atom->stripped_atomLabel is: " + ptnAtomLabel);
	    //debugInfo.out(" sc_atom ->stripped_atomLabel is: " + scAtomLabel);

	    double E = (*vdw_functor)(*ptn_atom, *sc_atom); // no need for /2.  not gonna bother with spline.
	    //double E = (*no_scale_functor)(*ptn_atom, *sc_atom); // if same backbone, don't use flat-bottom potential function.  Or maybe i should, don't know.
	    total_E += E;

// 	    	      if (E > 5.0) {
// 	    		cout << total_E << " is > 5.0kcal/mol, clash " << endl;
// 	    		(*sc_atom)->dump();
// 	    		(*ptn_atom)->dump();
// 	    	      }

	      //	      }
	      
	      //}
	  } // end should not exclude block
	} // end same residue block
	else { // i.e. if not on the same residue
	  //if (this->ON_THE_FLY == 1) {
	  //total_E += (*vdw_functor)(*ptn_atom, *sc_atom)/2;

	  //} 
	  //else if (this->ON_THE_FLY == 2) {
	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);

	  if (dist_sq > dist_sq_cutoff) {
	    continue;
	  } 
	  else {
	    double E = (*vdw_functor)(*ptn_atom, *sc_atom);  // no division by 2
	    total_E += E;
	    
// 	    	      if (E > 5.0) {
// 	    		cout << total_E << " is > 5.0kcal/mol, clash " << endl;
// 	    		(*sc_atom)->dump();
// 	    		(*ptn_atom)->dump();
// 	    	      }
	    
	  }
	  //}
	}
	
      } // end ptn_flag == 1 block (fixed atoms)
      
      
    } // end inner loop
  } // end outer loop
  
    //  cout << "Operations: " << i << endl;
    //  cout << "Double counts: " << double_count << endl;
  //cout << "     DEBUG::: survived main double loop" << endl;
  //delete vdw_functor;
  return total_E;

}

double VDW_EE::calc_EL_rot_moveableHET(const MutInfo& mutInfo, std::string method, double scale) {
  Debug debugInfo("VDW_EE::_calc_empty_lattice_E_on_the_fly_loop");
  
  /* method: FULL, FLAT, SCALE, FULL_NONONPOLARH, FLAT_NONONPOLARH */
  int i = 0;
  int double_count = 0;

  // Jan 16: if ASYM, use X6, the more accurate one.  If SYMM, the old one, X6 not written, hence LJ 6-12.
  // Mar 04: back to LJ.  Why did I even bother?  X6 much slower, need to take care of inflection point else E --> -inf.

  string base_method;
  stringV f;
  split(method, "_", f);
  base_method = f[0];

  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  if (f.size() == 2) {
    if (f[1] == "ASYM") {
      base_method = base_method + "_ASYM";
    }
    else if (f[1] == "NONONPOLARH" ) {
      nonPolar_H_calc = 0;
    }
  }

  if (f.size() == 3) {
    base_method += "_ASYM";
    nonPolar_H_calc = 0;
  }

  //debugInfo.out("Method used: " + method);
  //debugInfo.out("Base Method: " + base_method);
  //debugInfo.out("Nonpolar H calcalations? 1 = yes, 0 = no: " + string(itoa(nonPolar_H_calc)) );

  double total_E = 0;
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor, *no_scale_functor;
  if (base_method == "FULL") {
    vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "FLAT") {
    vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
  } else if (base_method == "SCALE") {
    // not implemented yet.
  } else if (base_method == "FULL_ASYM") {
    vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_full_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "FLAT_ASYM") {
    vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
    //vdw_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "SCALE_ASYM") {
    // no implemented yet
    cout << " This functor, SCALE_ASYM, not implemented yet." << endl;
    exit(2);
  } else {
    cout << " functor " << base_method << " not implemeneted yet." << endl;
    exit(2);
  }

  no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);
  //  no_scale_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, 0);

  // first get atoms from this sidechain.
  
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);


  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }

  

  // then get all relevant atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  //  if (this->ON_THE_FLY == 1) {
  atom_list = this->ptn->getAtomList();

    //}
//   else if (this->ON_THE_FLY == 2) {
//     atom_list = this->rotamerNeighborList->returnEmptyLatticeNeighborList(mutInfo);
//   }
  
  double dist_sq_cutoff = 10.0*10.0;

  
  int count_tot_E = 0;
  // Remake the two atom lists for nonPolar_H_calc.
  if (nonPolar_H_calc == 0) {
    ScreamAtomV noNonPolarH_atom_list; noNonPolarH_atom_list.clear();
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      //      if ( (*ptn_atom)->atomType.c_str()[0] != 'H') {
      //	noNonPolarH_atom_list.push_back(*ptn_atom);	
      //	continue;
      //      } // fast filter; string comparison too slow
      if ( "H_" == (*ptn_atom)->stripped_atomType ) {
	continue;
      }
      else {
	noNonPolarH_atom_list.push_back(*ptn_atom);
      }
    }
    ScreamAtomV noNonPolarH_sc_atom_list; noNonPolarH_sc_atom_list.clear();
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      //if ( (*sc_atom)->atomType.c_str()[0] != 'H') {
      //noNonPolarH_sc_atom_list.push_back(*sc_atom);
      //continue;
      //}
      if ( "H_" == (*sc_atom)->stripped_atomType ) {
	continue;      }
      else {
	noNonPolarH_sc_atom_list.push_back(*sc_atom);      
      }
    }
    atom_list.clear();
    sc_atoms.clear();
    atom_list.insert(atom_list.end(), noNonPolarH_atom_list.begin(), noNonPolarH_atom_list.end());
    sc_atoms.insert(sc_atoms.end(), noNonPolarH_sc_atom_list.begin(), noNonPolarH_sc_atom_list.end());
  }
  //debugInfo.out(" Size of atom_list: " + string(itoa(atom_list.size()) ) );
  //debugInfo.out(" Size of sc_atoms: " + string(itoa(sc_atoms.size())) );


  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }

  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {

      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      //if (mutInfo.chn == "Z") {
      if (mutInfoChn_is_Z) {
	if (*sc_atom == *ptn_atom) continue;
	if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
	       or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	  //if (this->ON_THE_FLY == 1) {
	  //total_E += (*vdw_functor)(*ptn_atom, *sc_atom)/2;
	  //count_tot_E++;
	  //} else if (this->ON_THE_FLY == 2) {
	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);
	  if (dist_sq > dist_sq_cutoff) {
	    continue;
	  } else {
	    double E = (*vdw_functor)(*ptn_atom, *sc_atom)/2;
	    total_E += E;

	    
	  }
	}
      }
      // Case 1: moveable atoms.  point for moveable and immoveable atoms if for optimization.
      else if (ptn_flag == 0) { // ptn_flag == 0: moveable atoms, i.e. those atoms on sidechains to be replaced.
	continue; // now that there's no need for internal SC-SC VDW E calculation, skip.
      }
      // Case 2: fixed atoms.
      else { // if ptn_flag ==1 
	// same residue atoms, including backbone atoms.
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  string scAtomLabel = (*sc_atom)->stripped_atomLabel;
	  string ptnAtomLabel = (*ptn_atom)->stripped_atomLabel;
	  // Note: the following lines need to change so that these atom names won't be hardcoded.
	  if ( ptnAtomLabel != "O" and
	       ptnAtomLabel != "HN" and
	       ptnAtomLabel != "OXT") {
	    // NEW SCHEME: SC-SC VDW not longer included in calculations, provided by rotamer library. 5-03-06.  HN, O, OXT's interactions with sidechain NOT provided by library because they are backbone dependent, therefore needs to be calulcated by this program.  Coulomb and HBond's still calculated because those terms were excluded in the calculations.
	    continue;
	  }
	  // following conditional: unnecessary, because sc atoms always away by at least 3 bonds from O, HN, or OXT.
	  if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
		 or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	    // CA <=> CG interactions already excluded by the previous condition.
	    //if (this->ON_THE_FLY == 1) {
	    //total_E += (*vdw_functor)(*ptn_atom, *sc_atom)/2;
	    //} else if (this->ON_THE_FLY == 2) {

	    
	    //double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);  // no longer checking distance; if on same residue, always calculate interaction
	    //if (dist_sq > dist_sq_cutoff) {
	    //continue;
	    //} else {
	    //debugInfo.out(" ptn_atom->stripped_atomLabel is: " + ptnAtomLabel);
	    //debugInfo.out(" sc_atom ->stripped_atomLabel is: " + scAtomLabel);

	    double E = (*vdw_functor)(*ptn_atom, *sc_atom); // no need for /2.  not gonna bother with spline.
	    //double E = (*no_scale_functor)(*ptn_atom, *sc_atom); // if same backbone, don't use flat-bottom potential function.  Or maybe i should, don't know.
	    total_E += E;

// 	    	      if (E > 5.0) {
// 	    		cout << total_E << " is > 5.0kcal/mol, clash " << endl;
// 	    		(*sc_atom)->dump();
// 	    		(*ptn_atom)->dump();
// 	    	      }

	      //	      }
	      
	      //}
	  } // end should not exclude block
	} // end same residue block
	else { // i.e. if not on the same residue
	  //if (this->ON_THE_FLY == 1) {
	  //total_E += (*vdw_functor)(*ptn_atom, *sc_atom)/2;

	  //} 
	  //else if (this->ON_THE_FLY == 2) {
	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);

	  if (dist_sq > dist_sq_cutoff) {
	    continue;
	  } 
	  else {
	    double E = (*vdw_functor)(*ptn_atom, *sc_atom);  // no division by 2
	    total_E += E;
	    
// 	    	      if (E > 5.0) {
// 	    		cout << total_E << " is > 5.0kcal/mol, clash " << endl;
// 	    		(*sc_atom)->dump();
// 	    		(*ptn_atom)->dump();
// 	    	      }
	    
	  }
	  //}
	}
	
      } // end ptn_flag == 1 block (fixed atoms)
      
      
    } // end inner loop
  } // end outer loop
  
    //  cout << "Operations: " << i << endl;
    //  cout << "Double counts: " << double_count << endl;
  //cout << "     DEBUG::: survived main double loop" << endl;
  //delete vdw_functor;
  return total_E;

}

void VDW_EE::setup_variableAtomsOnEachSidechain() {
  // this->_variable_atoms_on_each_sidechain
  this->_variable_atoms_on_each_sidechain.clear();
  this->_mutInfo_n_map.clear();

  /* Then make mutinfo maps. */
  int mutInfo_n = 1; // index count

  map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = this->mutInfo_rotConnInfo_map.begin();
  for (; mIrotC_itr != this->mutInfo_rotConnInfo_map.end(); ++mIrotC_itr) {
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;

    ScreamAtomV tmp_sc_atoms;
    
    if (rotConnInfo == NULL) {
      tmp_sc_atoms = ptn->get_sc_atoms(mutInfo);

    } else {
      // need to get variable atoms 
      tmp_sc_atoms = ptn->get_variable_atoms(rotConnInfo);
    }

    ScreamAtomV relevantAtoms;
    for (ScreamAtomVItr itr = tmp_sc_atoms.begin(); itr != tmp_sc_atoms.end(); ++itr) {
      if ( (*itr)->flags & 0x1 )
	continue;
      else
	relevantAtoms.push_back(*itr);  // want to populate everything that's variable sidechains; even though that have their energy expression currently set to fixed.  

    }
    
    this->_variable_atoms_on_each_sidechain[mutInfo_n] = relevantAtoms;
    this->_mutInfo_n_map[mutInfo_n] = &(mIrotC_itr->first); // will this work?
    ++mutInfo_n;
  }

  



}



void VDW_EE::_initVariableAndFixedAtomPairList(Protein* ptn, const vector<MutInfo> mutInfoV) {
  /* This routine initializes vector<AtomPair> contacts_between_variable_and_fixed list. */
  
  ScreamAtomV fixed_atoms, all_variable_atoms, variable_atoms_for_one_MutInfo;
  fixed_atoms.clear(), all_variable_atoms.clear(), variable_atoms_for_one_MutInfo.clear();

  
  vector<MutInfo>::const_iterator mutInfoItr = mutInfoV.begin();
  
  
  for (; mutInfoItr != mutInfoV.end(); ++mutInfoItr) {
    /* initializing variable atoms */
    string chn = (*mutInfoItr).chn;
    int pstn = (*mutInfoItr).pstn;
    
    //ScreamAtomV tmp_sc_atoms = ptn->get_sc_atoms(chn, pstn);
    ScreamAtomV tmp_sc_atoms = ptn->get_sc_atoms(*mutInfoItr);
    all_variable_atoms.insert(all_variable_atoms.end(), tmp_sc_atoms.begin(), tmp_sc_atoms.end());
  }

  /* initialize fixed atoms.  all non-variable atoms are fixed atoms. */
  ScreamAtomV ptn_atom_list = ptn->getAtomList();
  fixed_atoms.insert(fixed_atoms.end(), ptn_atom_list.begin(), ptn_atom_list.end());
  for (ScreamAtomVItr itr = all_variable_atoms.begin(); itr != all_variable_atoms.end(); ++itr) {
    ScreamAtomVItr variable_atom_i = find(fixed_atoms.begin(), fixed_atoms.end(), *itr);
    fixed_atoms.erase(variable_atom_i);
  }

  // initialize map<MutInfo, vector<atom-pair> > dictionary.  This allows each rotamer to store calculated energies.
  this->variable_and_fixed.clear();
  for (vector<MutInfo>::const_iterator mutInfoItr = mutInfoV.begin(); mutInfoItr != mutInfoV.end(); ++mutInfoItr) {
    string chn = (*mutInfoItr).chn;
    int pstn = (*mutInfoItr).pstn;
    //ScreamAtomV this_sc_atoms = ptn->get_sc_atoms(chn, pstn);
    ScreamAtomV this_sc_atoms = ptn->get_sc_atoms(*mutInfoItr);
    
    vector<AtomPair> this_sc_atom_pair_list; // sidechain pair list with fixed atoms on the protein.
    
    /* make atom pair list */
    /* first, initialize interaction between sidechain and fixed atoms. */
    for (ScreamAtomVItr Var_i = this_sc_atoms.begin(); Var_i != this_sc_atoms.end(); ++Var_i) {
      for (ScreamAtomVItr Fixed_i = fixed_atoms.begin(); Fixed_i != fixed_atoms.end(); ++Fixed_i) {
	// 1-2 and 1-3 exclusion
	if (! (scream_tools::should_exclude_on_1_2(*Var_i, *Fixed_i)
	       or scream_tools::should_exclude_on_1_3(*Var_i, *Fixed_i) ) ) 
	  // no HB exclusion; term substracted later on.
	  this_sc_atom_pair_list.push_back(AtomPair(*Var_i, *Fixed_i));
	
      }
    }
    /* then, initialize interaction within the sidechain */
    for (ScreamAtomVItr Var_i1 = this_sc_atoms.begin(); Var_i1 != this_sc_atoms.end(); ++Var_i1) {
      for (ScreamAtomVItr Var_i2 = Var_i1; Var_i2 != this_sc_atoms.end(); ++Var_i2) {
	if (! (scream_tools::should_exclude_on_1_2(*Var_i1, *Var_i2)
	       or scream_tools::should_exclude_on_1_3(*Var_i1, *Var_i2) ) )
	  this_sc_atom_pair_list.push_back(AtomPair(*Var_i1, *Var_i2));
      }
    }

    this->variable_and_fixed[*mutInfoItr] = this_sc_atom_pair_list;
  }

  for (map<MutInfo, vector<AtomPair> >::iterator itr = this->variable_and_fixed.begin(); 
       itr != this->variable_and_fixed.end(); ++itr) {
    cout << (*itr).first.AA << (*itr).first.pstn << "_" << (*itr).first.chn;
    cout << " has size " << (*itr).second.size() << endl;
  }

}


void VDW_EE::_initVariableAndFixedAtomPairListArb(Protein* ptn, const map<MutInfo, RotConnInfo*> mIrotC_map) {
  /* This routine initializes variable_and_fixed atom pair lists */
  ScreamAtomV fixed_atoms, all_variable_atoms, variable_atoms_for_one_MutInfo;
  fixed_atoms.clear(), all_variable_atoms.clear(), variable_atoms_for_one_MutInfo.clear();

  map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = mIrotC_map.begin();
  
  for (; mIrotC_itr != mIrotC_map.end(); ++mIrotC_itr) {
    /* initializing variable atoms */
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;
    string chn = mutInfo.chn;
    int pstn = mutInfo.pstn;

    if (rotConnInfo == NULL) {
      int pstn = mutInfo.pstn;
      //ScreamAtomV tmp_sc_atoms = ptn->get_sc_atoms(chn, pstn);
      ScreamAtomV tmp_sc_atoms = ptn->get_sc_atoms(mutInfo);
      all_variable_atoms.insert(all_variable_atoms.end(), tmp_sc_atoms.begin(), tmp_sc_atoms.end());
    } else {
      // need to get variable atoms 
      ScreamAtomV tmp_sc_atoms = ptn->get_variable_atoms(rotConnInfo);
      all_variable_atoms.insert(all_variable_atoms.end(), tmp_sc_atoms.begin(), tmp_sc_atoms.end());
    }
  }
  /* initialize fixed atoms.  all non-variable atoms are fixed atoms. */
  ScreamAtomV ptn_atom_list = ptn->getAtomList();
  fixed_atoms.insert(fixed_atoms.end(), ptn_atom_list.begin(), ptn_atom_list.end());
  for (ScreamAtomVItr itr = all_variable_atoms.begin(); itr != all_variable_atoms.end(); ++itr) {
    ScreamAtomVItr variable_atom_i = find(fixed_atoms.begin(), fixed_atoms.end(), *itr);
    fixed_atoms.erase(variable_atom_i);
  }
  /* initialize map<MutInfo, vector<atom-pair> > dictionary.  This allows each rotamer to store calculated energies.*/
  this->variable_and_fixed.clear();
  for (map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = mIrotC_map.begin(); mIrotC_itr != mIrotC_map.end(); ++mIrotC_itr) {
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;
    string chn = mutInfo.chn;
    int pstn = mutInfo.pstn;

    ScreamAtomV this_sc_atoms;
    vector<AtomPair> this_sc_atom_pair_list; // variable atom pair list with fixed atoms on the protein.

    if (rotConnInfo == NULL ) {
      //      this_sc_atoms = ptn->get_sc_atoms(chn, pstn);
      this_sc_atoms = ptn->get_sc_atoms(mutInfo);

    } else {
      this_sc_atoms = ptn->get_variable_atoms(rotConnInfo);

    }
    /* make atom pair list */
    /* first, initialize interaction between sidechain and fixed atoms. */
    
    for (ScreamAtomVItr Var_i = this_sc_atoms.begin(); Var_i != this_sc_atoms.end(); ++Var_i) {
      for (ScreamAtomVItr Fixed_i = fixed_atoms.begin(); Fixed_i != fixed_atoms.end(); ++Fixed_i) {
	// 1-2 and 1-3 exclusion
	if (! (scream_tools::should_exclude_on_1_2(*Var_i, *Fixed_i)
	       or scream_tools::should_exclude_on_1_3(*Var_i, *Fixed_i) ) )
	  // no HB exclusion; term subtracted later on
	  this_sc_atom_pair_list.push_back(AtomPair(*Var_i, *Fixed_i));
      }
    }
 
    /* then, initialize interaction within the sidechain */
    for (ScreamAtomVItr Var_i1 = this_sc_atoms.begin(); Var_i1 != this_sc_atoms.end(); ++Var_i1) {
      for (ScreamAtomVItr Var_i2 = Var_i1; Var_i2 != this_sc_atoms.end(); ++Var_i2) {
	if (! (scream_tools::should_exclude_on_1_2(*Var_i1, *Var_i2)
	       or scream_tools::should_exclude_on_1_3(*Var_i1, *Var_i2) ) )
	  this_sc_atom_pair_list.push_back(AtomPair(*Var_i1, *Var_i2));
      }
    }
    this->variable_and_fixed[mutInfo] = this_sc_atom_pair_list;
  }
  
  for (map<MutInfo, vector<AtomPair> >::iterator itr = this->variable_and_fixed.begin(); 
       itr != this->variable_and_fixed.end(); ++itr) {
    cout << (*itr).first.AA << (*itr).first.pstn << "_" << (*itr).first.chn;
    cout << " has size " << (*itr).second.size() << endl;
  }
}


void VDW_EE::_initVariableAndVariableAtomPairList(Protein* ptn, const vector<MutInfo> mutInfoV) {
  /* initializes interaction between variable atoms and variable atoms only; for scream */

  map<MutInfo, ScreamAtomV> variable_atoms_on_each_sidechain; ///< each MutInfo has its own variable atoms.

  vector<MutInfo>::const_iterator mutInfoItr = mutInfoV.begin();

  for (; mutInfoItr != mutInfoV.end(); ++mutInfoItr) {
    /* initializing sc atom list */
    string chn = (*mutInfoItr).chn;
    int pstn = (*mutInfoItr).pstn;
    
    //ScreamAtomV tmp_sc_atoms = ptn->get_sc_atoms(chn, pstn);
    ScreamAtomV tmp_sc_atoms = ptn->get_sc_atoms(*mutInfoItr);
    this->each_sc_atom_list[*mutInfoItr] = tmp_sc_atoms;
  }

  // now make all possible pairs from this map<MutInfo, ScreamAtomV> list

  for (map<MutInfo, ScreamAtomV>::iterator itr1 = this->each_sc_atom_list.begin(); itr1 != this->each_sc_atom_list.end(); ++itr1) {
    MutInfo mutInfo1 = (*itr1).first;
    ScreamAtomV atomList1 = (*itr1).second;
    map<MutInfo, ScreamAtomV>::iterator itr2 = itr1;
    ++itr2;
    for (; itr2 != this->each_sc_atom_list.end(); ++itr2) {

      MutInfo mutInfo2 = (*itr2).first;

      ScreamAtomV atomList2 = (*itr2).second;
      // now make pairs: both pair<MutInfo, MutInfo> and list of AtomPairs.
      MutInfoPair mutInfo_pair(mutInfo1, mutInfo2);
      vector<AtomPair> atomPair_list;
      for (ScreamAtomVConstItr a1 = atomList1.begin(); a1 != atomList1.end(); ++a1) {
	for (ScreamAtomVConstItr a2 = atomList2.begin(); a2 != atomList2.end(); ++a2) {
	  atomPair_list.push_back(AtomPair(*a1, *a2));
	}
      }
      this->variable_and_variable[mutInfo_pair] = atomPair_list;


    }
  }

}


void VDW_EE::_initVariableAndVariableAtomPairListArb(Protein* ptn, const map<MutInfo, RotConnInfo*> mIrotC_map) {
  /* initializes interaction between variable atoms and variable atoms only; for scream */
  map<MutInfo, ScreamAtomV> variable_atoms_on_each_sidechain; ///< each MutInfo has its own variable atoms.

  map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = mIrotC_map.begin();
  for (; mIrotC_itr != mIrotC_map.end(); ++mIrotC_itr) {
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;
    string chn = mutInfo.chn;
    int pstn = mutInfo.pstn;
    ScreamAtomV tmp_sc_atoms;

    if (rotConnInfo == NULL) {
      //tmp_sc_atoms = ptn->get_sc_atoms(chn, pstn);
      tmp_sc_atoms = ptn->get_sc_atoms(mutInfo);
    } else {
      // need to get variable atoms 
      tmp_sc_atoms = ptn->get_variable_atoms(rotConnInfo);

    }
    this->each_sc_atom_list[mutInfo] = tmp_sc_atoms;
  }
  
  // now make all possible pairs from this map<MutInfo, ScreamAtomV> list

  for (map<MutInfo, ScreamAtomV>::iterator itr1 = this->each_sc_atom_list.begin(); itr1 != this->each_sc_atom_list.end(); ++itr1) {
    MutInfo mutInfo1 = (*itr1).first;
    ScreamAtomV atomList1 = (*itr1).second;
    map<MutInfo, ScreamAtomV>::iterator itr2 = itr1;
    ++itr2;
    for (; itr2 != this->each_sc_atom_list.end(); ++itr2) {
      
      MutInfo mutInfo2 = (*itr2).first;
      
      ScreamAtomV atomList2 = (*itr2).second;
      // now make pairs: both pair<MutInfo, MutInfo> and list of AtomPairs.
      MutInfoPair mutInfo_pair(mutInfo1, mutInfo2);
      vector<AtomPair> atomPair_list;
      for (ScreamAtomVConstItr a1 = atomList1.begin(); a1 != atomList1.end(); ++a1) {
	for (ScreamAtomVConstItr a2 = atomList2.begin(); a2 != atomList2.end(); ++a2) {
	  atomPair_list.push_back(AtomPair(*a1, *a2));
	}
      }
      this->variable_and_variable[mutInfo_pair] = atomPair_list;
    }
  }
}

double VDW_EE::_calc_empty_lattice_E_on_the_fly_loop(const MutInfo& mutInfo, string method, double scale) {
  Debug debugInfo("VDW_EE::_calc_empty_lattice_E_on_the_fly_loop");
  /* method: FULL, FLAT, SCALE, FULL_NONONPOLARH, FLAT_NONONPOLARH */
  // Jan 16: if ASYM, use X6, the more accurate one.  If SYMM, the old one, X6 not written, hence LJ 6-12.
  // Mar 04: back to LJ.  Why did I even bother?  X6 much slower, need to take care of inflection point else E --> -inf.

  string base_method;
  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  int CBCalc = 1; // default = 1: calculates ground spectrum with presence of CB on variable standard sidechains.
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor, *no_scale_functor;

  this->_figureOutvdwPotentialMethods(method, base_method, nonPolar_H_calc, CBCalc, scale, &vdw_functor);
  //no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);
  this->_returnNoScaleFunctor(method, &no_scale_functor);

  //debugInfo.out("Method used: " + method);
  //debugInfo.out("Base Method: " + base_method);
  //debugInfo.out("Nonpolar H calcalations? 1 = yes, 0 = no: " + string(itoa(nonPolar_H_calc)) );

 /* Get Sidechain atoms from this sidechain, the mutInfo passed in as first argument.  Reminder: CBCalc case. */
  ScreamAtomV atom_list; 
  ScreamAtomV sc_atoms;
  ScreamAtomV self_bb_atoms; 

  this->_setupAtomListStuff(mutInfo, atom_list, sc_atoms, self_bb_atoms);

  this->_remakeForNonPolarHCalc(nonPolar_H_calc, atom_list, sc_atoms); // need to include self_bb_atoms as well, do this later...

  if (sc_atoms.size() == 1) { // a Glycine
    delete vdw_functor;
    vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);
  }

  //  debugInfo.out(" Size of atom_list: " + string(itoa(atom_list.size()) ) );
  //  debugInfo.out(" Size of sc_atoms: " + string(itoa(sc_atoms.size())) );

  /* Start energy calculations.  First calculate energies between sidechain and self bb. */
  if (sc_atoms.size() == 0) return 0;
  double total_E = 0;
  double dist_sq_cutoff = 8.0*8.0;
  int count_tot_E = 0;

  for (ScreamAtomVConstItr bb_atom = self_bb_atoms.begin(); bb_atom != self_bb_atoms.end(); ++bb_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      total_E += (*no_scale_functor)(*bb_atom, *sc_atom);// No need for further conditionals since: only 1-4 VDWs for x-HN, x-O will be calculated
    }
  }
  //  cout << "Timing: Calculating Sidechain-Self-backbone time: " << (float(t4)-float(t3)) / (float)CLOCKS_PER_SEC << endl;

  // NOTE ON ABOVE: NEW SCHEME: SC-SC VDW not longer included in calculations, provided by rotamer library. 5-03-06.  HN, O, OXT's interactions with sidechain NOT provided by library because they are backbone dependent, therefore needs to be calulcated by this program.  Coulomb and HBond's still calculated because those terms were excluded in the calculations.
  
  /*  Main double loop below. */

  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }

  ScreamAtomVConstItr atom_list_end = atom_list.end();
  ScreamAtomVConstItr sc_atoms_end = sc_atoms.end();

  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {   // Note: atom_list only contains atoms that are relevent, 1-4 exclusion still needs to be checked however.
    //    int ptn_flag = (*ptn_atom)->flags & 0x2;  // No longer necessary: moveable atoms already been removed in earlier section.
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      //      cout << (*sc_atom)->return_bgf_line() <<
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      if (mutInfoChn_is_Z) {
	if (*sc_atom == *ptn_atom) continue;
	if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
	       or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {

	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);
	  if (dist_sq > dist_sq_cutoff)
	    continue;
	  else {
	    double E = (*vdw_functor)(*ptn_atom, *sc_atom)/2;
	    total_E += E;
	  }
	}
      }

      else { // if ptn_flag ==1 
	// Only Fixed atoms that don't belong to the same backbone are included here.
	//debugInfo.out(" ptn_atom->stripped_atomLabel is: " + (*sc_atom)->stripped_atomLabel); // WTF? compiler isn't compiling these lines out?
	//	debugInfo.out(" sc_atom ->stripped_atomLabel is: " + (*ptn_atom)->stripped_atomLabel);
	
	double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);
	if (dist_sq > dist_sq_cutoff)
	  continue;
	else {
	  double E = (*vdw_functor)(*ptn_atom, *sc_atom);  // no division by 2
	  total_E += E;
	}
      } // end ptn_flag == 1 block (fixed atoms)
      
    } // end inner loop
  } // end outer loop


  delete vdw_functor; delete no_scale_functor;
  return total_E;
}


double VDW_EE::_calc_all_interaction_E_on_the_fly_loop(string mode, double r) {
  /* this routine also takes care of arb rotamers */
  double total_E = 0;
 
  /* First, set up various stuff.*/
  string base_method;
  double scale = r;
  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  int CBCalc = 1; // default = 1: calculates ground spectrum with presence of CB on variable standard sidechains.
  SCREAM_VDW_BASE_FUNCTIONAL_OBJ* vdw_functor, *no_scale_functor;

  this->_figureOutvdwPotentialMethods(mode, base_method, nonPolar_H_calc, CBCalc, scale, &vdw_functor);
  this->_returnNoScaleFunctor(mode, &no_scale_functor);
  //no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);

  if (this->_variable_atoms_on_each_sidechain.empty())
    this->setup_variableAtomsOnEachSidechain();

  /* Then, setup moveable atoms and fixed atoms.*/
  map<int, ScreamAtomV> fixed_mutInfo_atoms, moveable_mutInfo_atoms;
  map<int, ScreamAtomV>::iterator v_end = this->_variable_atoms_on_each_sidechain.end();
  ScreamAtomV fixed_atoms, moveable_atoms;
  for (map<int, ScreamAtomV>::iterator itr = this->_variable_atoms_on_each_sidechain.begin(); itr != v_end; ++itr) {
    int mutInfo_n = itr->first;
    ScreamAtomV * atoms = &(itr->second);
    
    for (ScreamAtomVConstItr atom = atoms->begin(); atom != atoms->end(); ++atom) {
      if ( ( (*atom)->flags & 0x4 ) == 0) // i.e. invisible
	continue; 
      if (  (*atom)->flags & 0x2)  {// i.e. fixed 
	fixed_atoms.push_back(*atom);
      }
      else {
	moveable_atoms.push_back(*atom);
      }
    }
    if (fixed_atoms.size() != 0)
      fixed_mutInfo_atoms[mutInfo_n] = fixed_atoms;
    if (moveable_atoms.size() != 0)
      moveable_mutInfo_atoms[mutInfo_n] = moveable_atoms;
    fixed_atoms.clear(); 
    moveable_atoms.clear();
  }

//   cout << "Total Moveable atoms: ";
//   int n=0;
//   for (map<int, ScreamAtomV>::iterator jjj = moveable_mutInfo_atoms.begin(); jjj != moveable_mutInfo_atoms.end(); ++jjj)
//     n+=(jjj->second).size();
//   cout << n << endl;
//   cout << "Total Fixed atoms: ";
//   n=0;
//   for (map<int, ScreamAtomV>::iterator jjj = fixed_mutInfo_atoms.begin(); jjj != fixed_mutInfo_atoms.end(); ++jjj) {
//     n+=(jjj->second).size();
//   }
//   cout << n << endl;
  
  /* Then, first calculate moveable-moveable atom interaction E.*/
  double dist_cutoff_sq = 8*8; // previously: 10*10.  Now: 8*8.  8A cutoff should be good.

  for (map<int, ScreamAtomV>::iterator itr1 = moveable_mutInfo_atoms.begin(); itr1 != moveable_mutInfo_atoms.end(); ++itr1) {
    map<int, ScreamAtomV>::iterator itr2 = itr1;     ++itr2;
    for (; itr2 != moveable_mutInfo_atoms.end(); ++itr2) {
      double interactionE=0;
      ScreamAtomV* list1 = &(itr1->second);      
      ScreamAtomV* list2 = &(itr2->second);

      for (ScreamAtomVConstItr a1 = list1->begin(); a1 != list1->end(); ++a1)
	for (ScreamAtomVConstItr a2 = list2->begin(); a2 != list2->end(); ++a2) {
	  double dist_sq = (*a1)->distance_squared(*a2);
	  if (dist_sq < dist_cutoff_sq) {
	    double E = (*vdw_functor)(*a1, *a2);
	    interactionE += E;
	  }
	}
      if (this->clashCollection) {
	this->clashCollection->addClashPair( *(_mutInfo_n_map[itr1->first]), *(_mutInfo_n_map[itr2->first]), interactionE);
	if (interactionE > 250)
	  cout << "Clashing pair (E > 250) in on_the_fly_loop, VDW calculations: " <<  *(_mutInfo_n_map[itr1->first]) << " " << *(_mutInfo_n_map[itr2->first]) << " " <<  interactionE << endl;
      }
      total_E += interactionE;
    }
  }
  
  /* Then, main double loop, over moveable-fixed atoms. */
  for (map<int, ScreamAtomV>::iterator itr1 = moveable_mutInfo_atoms.begin(); itr1 != moveable_mutInfo_atoms.end(); ++itr1)
    for (map<int, ScreamAtomV>::iterator itr2 = fixed_mutInfo_atoms.begin(); itr2 != fixed_mutInfo_atoms.end(); ++itr2) {
      double interactionE = 0;
      ScreamAtomV* list1 = &(itr1->second);      
      ScreamAtomV* list2 = &(itr2->second);

      for (ScreamAtomVConstItr a1 = list1->begin(); a1 != list1->end(); ++a1)
	for (ScreamAtomVConstItr a2 = list2->begin(); a2 != list2->end(); ++a2) {
	  double dist_sq = (*a1)->distance_squared(*a2);
	  if (dist_sq < dist_cutoff_sq) {
	    double E = (*vdw_functor)(*a1, *a2);
	    interactionE += E;
	  }
	}
      if (this->clashCollection) 
	this->clashCollection->addClashPair( *(_mutInfo_n_map[itr1->first]), *(_mutInfo_n_map[itr2->first]), interactionE);
      total_E += interactionE;
    }

  return total_E;
}


void VDW_EE::_figureOutvdwPotentialMethods(const string method, string& base_method, int& nonPolar_H_calc, int& CBCalc, double scale, SCREAM_VDW_BASE_FUNCTIONAL_OBJ* *vdw_functor) {
  /* This function figures out potential methods.*/

  /* Input format for method: */
  stringV f;
  split(method, "_", f);
  base_method = f[0];
  string LJOption = "12-6";
  /* Fields that needed to be searched: "ASYM", "NOCB", "NONONPOLARH". */
  for (int i=0; i<f.size(); ++i) {
    string s = f[i];
    if (s == "ASYM")
      base_method = base_method + "_ASYM";
    else if (s == "NONONPOLARH")
      nonPolar_H_calc = 0;
    else if (s == "NOCB")
      CBCalc = 0; // NOTE: CBCalc setup no longer explicitly handled in this routine.  the initFixedMoveable routine in scream_EE sets up the fixed/moveable flag value in atom.
    else // LJ option: either 12-6, 11-6, etc.
      LJOption = s;
  }


  /* Now set up functors for different energy functions */
  if (base_method == "FULL") {
    *vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
    // Not implemented yet: below
//     if (LJOption == "12-6")
//       *vdw_functor = new SCREAM_calc_full_delta_VDW_6_12(this->vdw_obj, scale);
//     else if (LJOption == "11-6") 
//       *vdw_functor = new SCREAM_calc_full_delta_VDW_6_11(this->vdw_obj, scale);
//     else if (LJOption == "10-6") 
//       *vdw_functor = new SCREAM_calc_full_delta_VDW_6_10(this->vdw_obj, scale);
//     else if (LJOption == "9-6") 
//       *vdw_functor = new SCREAM_calc_full_delta_VDW_6_9(this->vdw_obj, scale);
//     else if (LJOption == "8-6") 
//       *vdw_functor = new SCREAM_calc_full_delta_VDW_6_8(this->vdw_obj, scale);
//     else if (LJOption == "7-6") 
//       *vdw_functor = new SCREAM_calc_full_delta_VDW_6_7(this->vdw_obj, scale);

  } else if (base_method == "FLAT") {
    *vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
    // Not implemented yet: below
//     if (LJOption == "12-6")
//       *vdw_functor = new SCREAM_calc_flat_delta_VDW_6_12(this->vdw_obj, scale);
//     else if (LJOption == "11-6") 
//       *vdw_functor = new SCREAM_calc_flat_delta_VDW_6_11(this->vdw_obj, scale);
//     else if (LJOption == "10-6") 
//       *vdw_functor = new SCREAM_calc_flat_delta_VDW_6_10(this->vdw_obj, scale);
//     else if (LJOption == "9-6") 
//       *vdw_functor = new SCREAM_calc_flat_delta_VDW_6_9(this->vdw_obj, scale);
//     else if (LJOption == "8-6") 
//       *vdw_functor = new SCREAM_calc_flat_delta_VDW_6_8(this->vdw_obj, scale);
//     else if (LJOption == "7-6") 
//       *vdw_functor = new SCREAM_calc_flat_delta_VDW_6_7(this->vdw_obj, scale);

  } else if (base_method == "SCALE") {
    // not implemented yet.
  } else if (base_method == "FULL_ASYM") {
    //*vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
    if (LJOption == "12-6")
      *vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
    else if (LJOption == "11-6") 
      *vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_11(this->vdw_obj, scale);
    else if (LJOption == "10-6") 
      *vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_10(this->vdw_obj, scale);
    else if (LJOption == "9-6") 
      *vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_9(this->vdw_obj, scale);
    else if (LJOption == "8-6") 
      *vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_8(this->vdw_obj, scale);
    else if (LJOption == "7-6") 
      *vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_7(this->vdw_obj, scale);
    else {
      cout << "LJ parameter not found.  Using default 12-6." << endl;
      *vdw_functor = new SCREAM_calc_full_delta_asym_VDW_6_12(this->vdw_obj, scale);
    }
      
    //vdw_functor = new SCREAM_calc_full_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "FLAT_ASYM") {
    //    *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
    if (LJOption == "12-6")
      *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
    else if (LJOption == "11-6") 
      *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_11(this->vdw_obj, scale);
    else if (LJOption == "10-6") 
      *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_10(this->vdw_obj, scale);
    else if (LJOption == "9-6") 
      *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_9(this->vdw_obj, scale);
    else if (LJOption == "8-6") 
      *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_8(this->vdw_obj, scale);
    else if (LJOption == "7-6") 
      *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_7(this->vdw_obj, scale);
    else {
      cout << "LJ parameter not found.  Using default 12-6." << endl;
      *vdw_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, scale);
    }
    //vdw_functor = new SCREAM_calc_flat_delta_asym_X6(this->vdw_obj, scale);
  } else if (base_method == "SCALE_ASYM") {
    // no implemented yet
    cout << " This functor, SCALE_ASYM, not implemented yet." << endl;
    exit(2);
  } else {
    cout << " functor " << base_method << " not implemeneted yet." << endl;
    exit(2);
  }

  
}


void VDW_EE::_returnNoScaleFunctor(const string& method, SCREAM_VDW_BASE_FUNCTIONAL_OBJ* *no_scale_functor) {
  /* Input format for method: */
  stringV f;
  split(method, "_", f);
  string LJOption = "12-6";

  /* Fields that needed to be searched: "ASYM", "NOCB", "NONONPOLARH". */
  for (int i=0; i<f.size(); ++i) {
    string s = f[i];
    if (s == "ASYM")
      continue;
    else if (s == "NONONPOLARH")
      continue;
    else if (s == "NOCB")
      continue; // NOTE: CBCalc setup no longer explicitly handled in this routine.  the initFixedMoveable routine in scream_EE sets up the fixed/moveable flag value in atom.
    else // LJ option: either 12-6, 11-6, etc.
      LJOption = s;
  }

  if (LJOption == "12-6")
    *no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_12(this->vdw_obj, 0);
  else if (LJOption == "11-6") 
    *no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_11(this->vdw_obj, 0);
  else if (LJOption == "10-6") 
    *no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_10(this->vdw_obj, 0);
  else if (LJOption == "9-6") 
    *no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_9(this->vdw_obj, 0);
  else if (LJOption == "8-6") 
    *no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_8(this->vdw_obj, 0);
  else if (LJOption == "7-6") 
    *no_scale_functor = new SCREAM_calc_flat_delta_asym_VDW_6_7(this->vdw_obj, 0);

}

void VDW_EE::_setupAtomListStuff(const MutInfo& mutInfo, ScreamAtomV& atom_list, ScreamAtomV& sc_atoms, ScreamAtomV& self_bb_atoms) {
  ScreamAtomV tmp_sc_atoms; tmp_sc_atoms.clear();
  atom_list.clear(); sc_atoms.clear();self_bb_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);
  
  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) )                      // if ->second == NULL, natural AA.
    tmp_sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  else tmp_sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  
  for (ScreamAtomVItr itr = tmp_sc_atoms.begin(); itr != tmp_sc_atoms.end(); ++itr)  // initialize sc_atom list, after checking if the sidechain is considered to be empty.
    if ( ((*itr)->flags & 0x2) == 0) sc_atoms.push_back(*itr);
  
  /* Prepare all relevant atoms on the protein.  Remember CBCalc case. */
  atom_list = this->ptn->getAtomList();
  ScreamAtomV relevent_aL; relevent_aL.clear();

  int pstn = mutInfo.getPstn();
  string chain = mutInfo.getChn();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    int ptn_flag = ( ( (*itr)->flags & 0x2 ) and ( (*itr)->flags & 0x8) ); // 0x2: determined.  0x8
    if ( ptn_flag ) { // 0 is moveable, don't include.  anything else is fixed, include. 
      if ( (*itr)->resNum == pstn and (*itr)->chain == chain) { // since no residue chain should even be chain Z, no arbitrary rotamers "backbones" would be included here.
	if ( (*itr)->stripped_atomLabel == "HN" or (*itr)->stripped_atomLabel == "O" or (*itr)->stripped_atomLabel == "OXT") // only want these atoms; other atoms, like "N", "C", "CA" are not wanted in this list.
	  self_bb_atoms.push_back(*itr);
      }
      else
	relevent_aL.push_back(*itr);
    }
  } // CB case handled by scream_EE setup.

  atom_list.clear();
  atom_list.insert(atom_list.end(), relevent_aL.begin(), relevent_aL.end());

}


void VDW_EE::_remakeForNonPolarHCalc(int nonPolar_H_calc, ScreamAtomV& atom_list, ScreamAtomV& sc_atoms) {

  /* Remake the two atom lists for nonPolar_H_calc. */
  if (nonPolar_H_calc == 0) {
    ScreamAtomV noNonPolarH_atom_list; noNonPolarH_atom_list.clear();
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      if ( "H_" == (*ptn_atom)->stripped_atomType )
	continue;
      else
	noNonPolarH_atom_list.push_back(*ptn_atom);
    }
    ScreamAtomV noNonPolarH_sc_atom_list; noNonPolarH_sc_atom_list.clear();
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      if ( "H_" == (*sc_atom)->stripped_atomType )
	continue;
      else 
	noNonPolarH_sc_atom_list.push_back(*sc_atom);      
    }
    atom_list.clear(); sc_atoms.clear();
    atom_list.insert(atom_list.end(), noNonPolarH_atom_list.begin(), noNonPolarH_atom_list.end());
    sc_atoms.insert(sc_atoms.end(), noNonPolarH_sc_atom_list.begin(), noNonPolarH_sc_atom_list.end());
  }



}
