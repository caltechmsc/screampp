#include "defs.hpp"
#include "MutInfo.hpp"
#include "scream_hb_EE.hpp"
#include <cassert>
#include <algorithm>

#include "RotamerNeighborList.hpp"

HB_EE::HB_EE() {
  this->rotamerNeighborList = NULL;
}

HB_EE::HB_EE(Protein* ptn, vector<MutInfo> mutInfo_V, SCREAM_HB_OBJ* hb_obj_) : hb_obj(hb_obj_) {

  this->ptn = ptn;

  this->_initVariableAndFixedAtomTripleList(ptn, mutInfo_V);
  this->_initVariableAndVariableAtomTripleList(ptn, mutInfo_V);

}

HB_EE::~HB_EE() {
  
}

void HB_EE::init_after_addedMutInfoRotConnInfo(Protein* ptn, SCREAM_HB_OBJ* hb_obj_) {
  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->hb_obj = hb_obj_;

  this->_initVariableAndFixedAtomTripleListArb(ptn, this->mutInfo_rotConnInfo_map);
  this->_initVariableAndVariableAtomTripleListArb(ptn, this->mutInfo_rotConnInfo_map);

  this->ptn = ptn;
  this->ON_THE_FLY = 0;

}

void HB_EE::init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein* ptn, SCREAM_HB_OBJ* hb_obj_) {

  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->hb_obj = hb_obj_;

  //  this->_initFixedMoveableAtomsOnProtein(ptn, this->mutInfo_rotConnInfo_map);
  this->ptn = ptn;
  this->ON_THE_FLY = 1;


}

void HB_EE::init_after_addedMutInfoRotConnInfo_neighbor_list(Protein* ptn, SCREAM_HB_OBJ* hb_obj_, RotamerNeighborList* rNL) {
  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->hb_obj = hb_obj_;
  
  this->ptn = ptn;
  this->ON_THE_FLY = 2;
  this->rotamerNeighborList = rNL;
  

}

void HB_EE::addMutInfoRotConnInfo(MutInfo mutInfo, RotConnInfo* rotConnInfo = NULL) {
  /* adds a mutInfo, rotconninfo pair to mutInfo_rotConnInfo_map */
  this->mutInfo_rotConnInfo_map[mutInfo] = rotConnInfo;
}


double HB_EE::calc_empty_lattice_E(const MutInfo mutInfo) {
  double total_E = 0;
  vector< vector<SCREAM_ATOM*> > atomTripleList(this->variable_and_fixed[mutInfo]);

  //  clock_t t1 = clock();

  if (! this->ON_THE_FLY) {

    for (vector< vector<SCREAM_ATOM*> >::iterator ap_i = atomTripleList.begin(); ap_i != atomTripleList.end(); ++ap_i) {
      total_E += this->hb_obj->calc_HB_Dre( (*ap_i)[0], (*ap_i)[1], (*ap_i)[2]);
    }
  }

  else {

    //total_E += this->_calc_empty_lattice_E_on_the_fly_loop(mutInfo, "FLAT", 0);

  }
  //  clock_t t2 = clock();
  //  cout << "that took: " << (double)(t2 - t1) / (double) CLOCKS_PER_SEC << endl;

  return total_E;

}

double HB_EE::calc_empty_lattice_E_delta(const MutInfo& mutInfo, string mode, double r) {
  double total_E = 0;
  vector< vector<SCREAM_ATOM*> > atomTripleList(this->variable_and_fixed[mutInfo]);

  if ( ! this->ON_THE_FLY) {

    if (mode == "FULL") {
      for (vector< vector<SCREAM_ATOM*> >::iterator ap_i = atomTripleList.begin(); ap_i != atomTripleList.end(); ++ap_i) {
	total_E += this->hb_obj->calc_full_delta_HB( (*ap_i)[0], (*ap_i)[1], (*ap_i)[2], r );
      }
    }
    if (mode == "FLAT") {
      for (vector< vector<SCREAM_ATOM*> >::iterator ap_i = atomTripleList.begin(); ap_i != atomTripleList.end(); ++ap_i) {
	total_E += this->hb_obj->calc_flat_delta_HB( (*ap_i)[0], (*ap_i)[1], (*ap_i)[2], r );
      }
    }
    if (mode == "SCALED") {
      for (vector< vector<SCREAM_ATOM*> >::iterator ap_i = atomTripleList.begin(); ap_i != atomTripleList.end(); ++ap_i) {
	total_E += this->hb_obj->calc_scaled_inner_wall_HB( (*ap_i)[0], (*ap_i)[1], (*ap_i)[2], r );
      }
    }
    if (mode == "RESIDUE") {
      for (vector< vector<SCREAM_ATOM*> >::iterator ap_i = atomTripleList.begin(); ap_i != atomTripleList.end(); ++ap_i) {
	total_E += this->hb_obj->calc_residue_delta_HB( (*ap_i)[0], (*ap_i)[1], (*ap_i)[2] );
      }
    }
  } 
  else {
    total_E = this->_calc_empty_lattice_E_on_the_fly_loop(mutInfo, mode, r);

  }
  return total_E;
}


double HB_EE::calc_residue_interaction_E(const MutInfo mI) {
// calculates interaction between MutInfo residue and rest of the variable atoms.
  double total_E = 0;
  for (map< MutInfoPair, vector< ScreamAtomV > >::iterator itr = this->variable_and_variable.begin();
       itr != this->variable_and_variable.end(); ++itr) {
    if ( mI == (itr->first).mutInfo1 or mI == (itr->first).mutInfo2 ) {
      for (vector< ScreamAtomV >::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); ++AHD_i) {
	 SCREAM_ATOM* A = (*AHD_i)[0];
	 SCREAM_ATOM* H = (*AHD_i)[1];
	 SCREAM_ATOM* D = (*AHD_i)[2];
	 //total_E += this->hb_obj->calc_HB_Dre(A, H, D);
	 total_E += this->hb_obj->calc_Scream_HB(A, H, D );
      }
    }
  }
    
  return total_E;

}

double HB_EE::calc_residue_interaction_E(const MutInfo mI1, const MutInfo mI2, string mode, double s) {
  
  double total_E = 0;
  SCREAM_HB_BASE_FUNCTIONAL_OBJ* hb_functor;
  if (mode == "FULL" or mode == "FULL_ASYM") {
    hb_functor = new SCREAM_calc_full_delta_HB(this->hb_obj, s);
    //cout << "FULL functor initiated" << endl;
  } else if (mode == "FLAT" or mode == "FLAT_ASYM") {
    hb_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, s);
    //cout << "FLAT functor initiated" << endl;
  } else if (mode == "SCALE") {
    // not implemented yet.
    cout << "SCALE functor not implemented" << endl;
  } else {
    cout << "no functors initiated" << endl;
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

  
  double R_on = this->hb_obj->R_on;
  double R_off = this->hb_obj->R_off;
  double theta_on = this->hb_obj->theta_on;
  double theta_off = this->hb_obj->theta_off;

  for (ScreamAtomVConstItr mI1_atom_itr = mI1_atoms.begin(); mI1_atom_itr != mI1_atoms.end(); ++mI1_atom_itr) {
    for (ScreamAtomVConstItr mI2_atom_itr = mI2_atoms.begin(); mI2_atom_itr != mI2_atoms.end(); ++mI2_atom_itr) {

      SCREAM_ATOM *A, *H, *D;
      
      if ( (*mI1_atom_itr)->atomType != "H___A" and (!(this->potential_HB_acceptor(*mI1_atom_itr))) ) continue;
      if ( (*mI1_atom_itr)->atomType == "H___A" ) {
	A = NULL;
	H = *mI1_atom_itr;
	D = H->connectivity_m.begin()->first; // H___A should be connected to only 1 O... will be exceptions, but ignore this for now.
      }
      else {
	A = *mI1_atom_itr;
	H = NULL;
	D = NULL;
      }
      if (*mI1_atom_itr == *mI2_atom_itr) continue; // self.
      // Define A, H and D.
      if (A == NULL) {
	if (!(this->potential_HB_acceptor(*mI2_atom_itr)) ) {
	  continue;
	} else {
	  A = *mI2_atom_itr; 
	}
      }
      else if (A != NULL) {
	if ( (*mI2_atom_itr)->atomType != "H___A") {
	  continue;
	} else {
	  H = *mI2_atom_itr;
	  D = H->connectivity_m.begin()->first;
	}
      }
      
      int ptn_flag = (*mI2_atom_itr)->flags & 0x2;
      // Case 0: arb lib cases, not optimized.
      if (! (scream_tools::should_exclude_on_1_2(A, D)
	     or scream_tools::should_exclude_on_1_3(A, D) ) ) {
	//total_E += this->hb_obj->calc_VDW_6_12( *mI2_atom_itr, *mI1_atom_itr);   // then only calc energies if 1-2 and 1-3 exclusions are satisfied.
	total_E += (*hb_functor)(A, H, D);
      }
	
    }
  }
  
  return total_E;

}


double HB_EE::calc_all_interaction_E() {
 // calculates interaction energies between all atom pairs
  double total_E = 0;

  if ( ! ON_THE_FLY) {

    for (map< MutInfoPair, vector<ScreamAtomV> >::iterator itr = this->variable_and_variable.begin(); 
	 itr != this->variable_and_variable.end(); ++itr) {
      for (vector<ScreamAtomV>::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); ++AHD_i) {
	SCREAM_ATOM* A = (*AHD_i)[0];
	SCREAM_ATOM* H = (*AHD_i)[1];
	SCREAM_ATOM* D = (*AHD_i)[2];
	
	//total_E += this->hb_obj->calc_HB_Dre(A, H, D);
	total_E += this->hb_obj->calc_Scream_HB(A, H, D);
      }
    }

  } else {
    total_E += this->_calc_all_interaction_E_on_the_fly_loop("FLAT", 0);

  }
  return total_E;

}



double HB_EE::calc_all_interaction_E_delta(string mode, double r) {
 // calculates interaction energies between all atom pairs
  double total_E = 0;
  
  if (! ON_THE_FLY) {
    
    if (mode == "FULL") {
      for (map< MutInfoPair, vector<ScreamAtomV> >::iterator itr = this->variable_and_variable.begin(); 
	   itr != this->variable_and_variable.end(); ++itr) {
	for (vector<ScreamAtomV>::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); ++AHD_i) {
	  SCREAM_ATOM* A = (*AHD_i)[0];
	  SCREAM_ATOM* H = (*AHD_i)[1];
	  SCREAM_ATOM* D = (*AHD_i)[2];
	  
	  total_E += this->hb_obj->calc_full_delta_HB(A, H, D, r);
	}
      }
    }
    else if (mode == "FLAT") {
      for (map< MutInfoPair, vector<ScreamAtomV> >::iterator itr = this->variable_and_variable.begin(); 
	   itr != this->variable_and_variable.end(); ++itr) {
	for (vector<ScreamAtomV>::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); ++AHD_i) {
	  SCREAM_ATOM* A = (*AHD_i)[0];
	  SCREAM_ATOM* H = (*AHD_i)[1];
	  SCREAM_ATOM* D = (*AHD_i)[2];
	  
	  total_E += this->hb_obj->calc_flat_delta_HB(A, H, D,r);
	}
      }
    }
    
    else if (mode == "RESIDUE") {
      for (map< MutInfoPair, vector<ScreamAtomV> >::iterator itr = this->variable_and_variable.begin(); 
	   itr != this->variable_and_variable.end(); ++itr) {
	for (vector<ScreamAtomV>::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); ++AHD_i) {
	  SCREAM_ATOM* A = (*AHD_i)[0];
	  SCREAM_ATOM* H = (*AHD_i)[1];
	  SCREAM_ATOM* D = (*AHD_i)[2];
	  
	  total_E += this->hb_obj->calc_residue_delta_HB(A, H, D);
	}
      }
      
    }
    else if (mode == "SCALED") {
      for (map< MutInfoPair, vector<ScreamAtomV> >::iterator itr = this->variable_and_variable.begin(); 
	   itr != this->variable_and_variable.end(); ++itr) {
	for (vector<ScreamAtomV>::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); ++AHD_i) {
	  SCREAM_ATOM* A = (*AHD_i)[0];
	  SCREAM_ATOM* H = (*AHD_i)[1];
	  SCREAM_ATOM* D = (*AHD_i)[2];
	  
	  total_E += this->hb_obj->calc_scaled_inner_wall_HB(A, H, D,r);
	}
      }
    }
  } 
  else { //  ON_THE_FLY
    total_E += this->_calc_all_interaction_E_on_the_fly_loop(mode, r);
  }
  return total_E;

}




void HB_EE::_initVariableAndFixedAtomTripleList(Protein* ptn, const vector<MutInfo> mutInfoV) {

  // first find list of H___A atoms and the oxygen attached to it.

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

  ScreamAtomV fixed_atom_H_Acceptors;
  for (ScreamAtomVItr itr = fixed_atoms.begin(); itr != fixed_atoms.end(); ++itr) {
    if (this->potential_HB_acceptor(*itr )) 
      fixed_atom_H_Acceptors.push_back(*itr);
  }
  
  
  ScreamAtomV fixed_atom_H___As;
  for (ScreamAtomVItr itr = fixed_atoms.begin(); itr != fixed_atoms.end(); ++itr) {
    if ( (*itr)->stripped_atomType == "H___A") {
      fixed_atom_H___As.push_back(*itr);
    }

  }
  
  // initialize map<MutInfo, vector< vector<SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*> > dictionary.  This allows each rotamer to store calculated energies.
  this->variable_and_fixed.clear();
  for (vector<MutInfo>::const_iterator mutInfoItr = mutInfoV.begin(); mutInfoItr != mutInfoV.end(); ++mutInfoItr) {
    string chn = (*mutInfoItr).chn;
    int pstn = (*mutInfoItr).pstn;
    //ScreamAtomV this_sc_atoms = ptn->get_sc_atoms(chn, pstn);
    ScreamAtomV this_sc_atoms = ptn->get_sc_atoms(*mutInfoItr);
    vector< vector<SCREAM_ATOM*> > this_sc_HB_atom_list;

    /* make atom pair list */
    for (ScreamAtomVItr Var_i = this_sc_atoms.begin(); Var_i != this_sc_atoms.end(); ++Var_i) {
      if  ( ((*Var_i)->atomType != "H___A") and (!(this->potential_HB_acceptor(*Var_i))) ) continue;
      /* initialize H___A on variable , acceptor on fixed */
      if ( (*Var_i)->atomType == "H___A" ) {
	SCREAM_ATOM* H = *Var_i;
	SCREAM_ATOM* D = H->connectivity_m.begin()->first; //  H___A should be connected to only 1 O... will be exceptions, but ignore this for now.
	for (ScreamAtomVItr Fixed_HB_A = fixed_atom_H_Acceptors.begin(); Fixed_HB_A != fixed_atom_H_Acceptors.end(); ++Fixed_HB_A) {
	  //if (D->distance(*Fixed_HB_A) > this->hb_obj->R_off) continue; // don't do this here
	  vector<SCREAM_ATOM*> HB_atoms;
	  HB_atoms.push_back(*Fixed_HB_A);
	  HB_atoms.push_back(H);
	  HB_atoms.push_back(D);
	  this_sc_HB_atom_list.push_back(HB_atoms);	  
	}
      }
      /* initialize accpetor on variable, H___A on donor */
      if ( this->potential_HB_acceptor(*Var_i) ) {
	SCREAM_ATOM* A = *Var_i;
	
	for (ScreamAtomVItr Fixed_HB_H = fixed_atom_H___As.begin(); Fixed_HB_H != fixed_atom_H___As.end(); ++Fixed_HB_H) {
	  vector<SCREAM_ATOM*> HB_atoms;
	  SCREAM_ATOM* H = *Fixed_HB_H;
	  SCREAM_ATOM* D = H->connectivity_m.begin()->first;
	  HB_atoms.push_back(A);
	  HB_atoms.push_back(H);
	  HB_atoms.push_back(D);
	  this_sc_HB_atom_list.push_back(HB_atoms);
	}
      }
      
      
    }
    cout << "Number of H-bonds: " << this_sc_HB_atom_list.size() << endl;
    this->variable_and_fixed[*mutInfoItr] = this_sc_HB_atom_list;
  }
  
  for (map<MutInfo, vector< vector<SCREAM_ATOM*> > >::iterator itr = this->variable_and_fixed.begin(); 
       itr != this->variable_and_fixed.end(); ++itr) {
    cout << (*itr).first.AA << (*itr).first.pstn << "_" << (*itr).first.chn;
    cout << " has size " << (*itr).second.size() << endl;
  }
  
  
}

void HB_EE::_initVariableAndFixedAtomTripleListArb(Protein* ptn, const map<MutInfo, RotConnInfo*> mIrotC_map) {

  // first find list of H___A atoms and the oxygen attached to it.

  ScreamAtomV fixed_atoms, all_variable_atoms, variable_atoms_for_one_MutInfo;
  fixed_atoms.clear(), all_variable_atoms.clear(), variable_atoms_for_one_MutInfo.clear();

  map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = mIrotC_map.begin();
  for (; mIrotC_itr != mIrotC_map.end(); ++mIrotC_itr) {
    /* initializing variable atoms */
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
    
    all_variable_atoms.insert(all_variable_atoms.end(), tmp_sc_atoms.begin(), tmp_sc_atoms.end());
  }
  
  /* initialize fixed atoms.  all non-variable atoms are fixed atoms. */
  ScreamAtomV ptn_atom_list = ptn->getAtomList();
  fixed_atoms.insert(fixed_atoms.end(), ptn_atom_list.begin(), ptn_atom_list.end());
  for (ScreamAtomVItr itr = all_variable_atoms.begin(); itr != all_variable_atoms.end(); ++itr) {
    ScreamAtomVItr variable_atom_i = find(fixed_atoms.begin(), fixed_atoms.end(), *itr);
    fixed_atoms.erase(variable_atom_i);
  }

  ScreamAtomV fixed_atom_H_Acceptors;
  for (ScreamAtomVItr itr = fixed_atoms.begin(); itr != fixed_atoms.end(); ++itr) {
    if (this->potential_HB_acceptor(*itr )) 
      fixed_atom_H_Acceptors.push_back(*itr);
  }
  
  ScreamAtomV fixed_atom_H___As;
  for (ScreamAtomVItr itr = fixed_atoms.begin(); itr != fixed_atoms.end(); ++itr) {
    if ( (*itr)->stripped_atomType == "H___A") {
      fixed_atom_H___As.push_back(*itr);
    }

  }

  // initialize map<MutInfo, vector< vector<SCREAM_ATOM*, SCREAM_ATOM*, SCREAM_ATOM*> > dictionary.  This allows each rotamer to store calculated energies.
  this->variable_and_fixed.clear();
  for (map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = mIrotC_map.begin(); mIrotC_itr != mIrotC_map.end(); ++mIrotC_itr) {
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;
    string chn = mutInfo.chn;
    int pstn = mutInfo.pstn;

    ScreamAtomV this_sc_atoms;
    vector< vector<SCREAM_ATOM*> > this_sc_HB_atom_list;

    if (rotConnInfo == NULL ) {
      //  this_sc_atoms = ptn->get_sc_atoms(chn, pstn);
      this_sc_atoms = ptn->get_sc_atoms(mutInfo);
    } 
    else { this_sc_atoms = ptn->get_variable_atoms(rotConnInfo); }

    /* make atom pair list */
    for (ScreamAtomVItr Var_i = this_sc_atoms.begin(); Var_i != this_sc_atoms.end(); ++Var_i) {
      if  ( ((*Var_i)->atomType != "H___A") and (!(this->potential_HB_acceptor(*Var_i))) ) continue;
      /* initialize H___A on variable , acceptor on fixed */
      if ( (*Var_i)->atomType == "H___A" ) {

	SCREAM_ATOM* H = *Var_i;
	SCREAM_ATOM* D = H->connectivity_m.begin()->first; //  H___A should be connected to only 1 O... will be exceptions, but ignore this for now.
	for (ScreamAtomVItr Fixed_HB_A = fixed_atom_H_Acceptors.begin(); Fixed_HB_A != fixed_atom_H_Acceptors.end(); ++Fixed_HB_A) {
	  //if (D->distance(*Fixed_HB_A) > this->hb_obj->R_off) continue; // don't do this here
	  vector<SCREAM_ATOM*> HB_atoms;
	  HB_atoms.push_back(*Fixed_HB_A);
	  HB_atoms.push_back(H);
	  HB_atoms.push_back(D);
	  this_sc_HB_atom_list.push_back(HB_atoms);	  
	}
      }
      /* initialize accpetor on variable, H___A on donor */
      if ( this->potential_HB_acceptor(*Var_i) ) {
	SCREAM_ATOM* A = *Var_i;

	for (ScreamAtomVItr Fixed_HB_H = fixed_atom_H___As.begin(); Fixed_HB_H != fixed_atom_H___As.end(); ++Fixed_HB_H) {
	  vector<SCREAM_ATOM*> HB_atoms;
	  SCREAM_ATOM* H = *Fixed_HB_H;
	  SCREAM_ATOM* D = H->connectivity_m.begin()->first;
	  HB_atoms.push_back(A);
	  HB_atoms.push_back(H);
	  HB_atoms.push_back(D);
	  this_sc_HB_atom_list.push_back(HB_atoms);
	}
      }
    }

    
    cout << "Number of H-bonds: " << this_sc_HB_atom_list.size() << endl;
    this->variable_and_fixed[mutInfo] = this_sc_HB_atom_list;
  }
  
  for (map<MutInfo, vector< vector<SCREAM_ATOM*> > >::iterator itr = this->variable_and_fixed.begin(); 
       itr != this->variable_and_fixed.end(); ++itr) {
    cout << (*itr).first.AA << (*itr).first.pstn << "_" << (*itr).first.chn;
    cout << " has size " << (*itr).second.size() << endl;
  }
  
  
}


void HB_EE::_initVariableAndVariableAtomTripleList(Protein* ptn, const vector<MutInfo> mutInfoV) {
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

  // instantiate list for <H___A, Donor> SCREAM_ATOM* pair on siechains.
  // instantiate list for potential H acceptor atoms on sidechains.  both done in same loop.
  map<MutInfo, vector<AtomPair> > HD_list;
  map<MutInfo, ScreamAtomV> A_list;
  vector<AtomPair> tmp_HD_list;
  ScreamAtomV tmp_A_list;
  for (map<MutInfo, ScreamAtomV>::iterator itr = this->each_sc_atom_list.begin(); itr != this->each_sc_atom_list.end(); ++itr) {
    MutInfo mutInfo = itr->first;
    tmp_HD_list.clear();
    tmp_A_list.clear();
    for (ScreamAtomV::iterator a = (*itr).second.begin(); a != (*itr).second.end(); ++a) {
      if ( (*a)->stripped_atomType == "H___A") {
	SCREAM_ATOM* H = *a;
	SCREAM_ATOM* D = H->connectivity_m.begin()->first; //  H___A should be connected to only 1 O... will be exceptions, but ignore this for now.
	AtomPair aP(H, D);
	tmp_HD_list.push_back(aP);
      }
      if ( this->potential_HB_acceptor(*a) ) {
	tmp_A_list.push_back(*a);
      }
    }

    HD_list[mutInfo] = tmp_HD_list;
    A_list[mutInfo] = tmp_A_list;

  }
  
  
  for (map<MutInfo, vector<AtomPair> >::iterator HD_itr = HD_list.begin(); HD_itr != HD_list.end(); ++HD_itr) {
    MutInfo mutInfo1 = (*HD_itr).first;
    vector<AtomPair> HD_v = (*HD_itr).second;
    for (map<MutInfo, ScreamAtomV>::iterator A_itr = A_list.begin(); A_itr != A_list.end(); ++A_itr) {
      MutInfo mutInfo2 = (*A_itr).first;
      ScreamAtomV A_v = (*A_itr).second;
      if (mutInfo1 == mutInfo2) {
	continue;
      }
      vector<ScreamAtomV> tmp_AHD_list;
      for (vector<AtomPair>::iterator HD_ap = HD_v.begin(); HD_ap != HD_v.end(); ++HD_ap) {
	for (ScreamAtomV::iterator A_a = A_v.begin(); A_a != A_v.end(); ++A_a) {
	  ScreamAtomV tmp_triple;
	  tmp_triple.push_back(*A_a);
	  tmp_triple.push_back((*HD_ap).a1);
	  tmp_triple.push_back((*HD_ap).a2);

	  tmp_AHD_list.push_back(tmp_triple);

	}
      }
      MutInfoPair mP(mutInfo1, mutInfo2);
      
      if (variable_and_variable.end() == variable_and_variable.find(mP) ) {
	this->variable_and_variable[mP] = tmp_AHD_list;
      } else {
	vector<ScreamAtomV> orig_AHD_list = this->variable_and_variable[mP];
	orig_AHD_list.insert(orig_AHD_list.end(), tmp_AHD_list.begin(), tmp_AHD_list.end());
	this->variable_and_variable[mP] = orig_AHD_list;
      }
    }

  }
  for (map<MutInfoPair, vector< vector<SCREAM_ATOM*> > >::iterator itr = this->variable_and_variable.begin(); itr != this->variable_and_variable.end(); ++itr) {
    cout << (*itr).first.mutInfo1 << " " << (*itr).first.mutInfo2 << " size: " << (*itr).second.size() << endl;;
  }
}



void HB_EE::_initVariableAndVariableAtomTripleListArb(Protein* ptn, const map<MutInfo, RotConnInfo*> mIrotC_map) {
    /* initializes interaction between variable atoms and variable atoms only; for scream */

  map<MutInfo, ScreamAtomV> variable_atoms_on_each_sidechain; ///< each MutInfo has its own variable atoms.
  map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = mIrotC_map.begin();

  for (; mIrotC_itr != mIrotC_map.end(); ++mIrotC_itr) {
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;
    string chn = mutInfo.chn;
    int pstn = mutInfo.pstn;
    ScreamAtomV tmp_sc_atoms;
    
    //if (rotConnInfo == NULL) { tmp_sc_atoms = ptn->get_sc_atoms(chn, pstn);} 
    if (rotConnInfo == NULL) { tmp_sc_atoms = ptn->get_sc_atoms(mutInfo);} 
    else { tmp_sc_atoms = ptn->get_variable_atoms(rotConnInfo); }

    this->each_sc_atom_list[mutInfo] = tmp_sc_atoms;
  }

  // instantiate list for <H___A, Donor> SCREAM_ATOM* pair on siechains.
  // instantiate list for potential H acceptor atoms on sidechains.  both done in same loop.
  map<MutInfo, vector<AtomPair> > HD_list;
  map<MutInfo, ScreamAtomV> A_list;
  vector<AtomPair> tmp_HD_list;
  ScreamAtomV tmp_A_list;
  for (map<MutInfo, ScreamAtomV>::iterator itr = this->each_sc_atom_list.begin(); itr != this->each_sc_atom_list.end(); ++itr) {
    MutInfo mutInfo = itr->first;
    tmp_HD_list.clear();
    tmp_A_list.clear();
    for (ScreamAtomV::iterator a = (*itr).second.begin(); a != (*itr).second.end(); ++a) {
      if ( (*a)->stripped_atomType == "H___A") {
	SCREAM_ATOM* H = *a;
	SCREAM_ATOM* D = H->connectivity_m.begin()->first; //  H___A should be connected to only 1 O... will be exceptions, but ignore this for now.
	AtomPair aP(H, D);
	tmp_HD_list.push_back(aP);
      }
      if ( this->potential_HB_acceptor(*a) ) {
	tmp_A_list.push_back(*a);
      }
    }

    HD_list[mutInfo] = tmp_HD_list;
    A_list[mutInfo] = tmp_A_list;

  }
  
  
  for (map<MutInfo, vector<AtomPair> >::iterator HD_itr = HD_list.begin(); HD_itr != HD_list.end(); ++HD_itr) {
    MutInfo mutInfo1 = (*HD_itr).first;
    vector<AtomPair> HD_v = (*HD_itr).second;
    for (map<MutInfo, ScreamAtomV>::iterator A_itr = A_list.begin(); A_itr != A_list.end(); ++A_itr) {
      MutInfo mutInfo2 = (*A_itr).first;
      ScreamAtomV A_v = (*A_itr).second;
      if (mutInfo1 == mutInfo2) {
	continue;
      }
      vector<ScreamAtomV> tmp_AHD_list;
      for (vector<AtomPair>::iterator HD_ap = HD_v.begin(); HD_ap != HD_v.end(); ++HD_ap) {
	for (ScreamAtomV::iterator A_a = A_v.begin(); A_a != A_v.end(); ++A_a) {
	  ScreamAtomV tmp_triple;
	  tmp_triple.push_back(*A_a);
	  tmp_triple.push_back((*HD_ap).a1);
	  tmp_triple.push_back((*HD_ap).a2);

	  tmp_AHD_list.push_back(tmp_triple);

	}
      }
      MutInfoPair mP(mutInfo1, mutInfo2);
      
      if (variable_and_variable.end() == variable_and_variable.find(mP) ) {
	this->variable_and_variable[mP] = tmp_AHD_list;
      } else {
	vector<ScreamAtomV> orig_AHD_list = this->variable_and_variable[mP];
	orig_AHD_list.insert(orig_AHD_list.end(), tmp_AHD_list.begin(), tmp_AHD_list.end());
	this->variable_and_variable[mP] = orig_AHD_list;
      }
    }

  }

  for (map<MutInfoPair, vector< vector<SCREAM_ATOM*> > >::iterator itr = this->variable_and_variable.begin(); itr != this->variable_and_variable.end(); ++itr) {
    cout << (*itr).first.mutInfo1 << " " << (*itr).first.mutInfo2 << " size: " << (*itr).second.size() << endl;;
  }
}


double HB_EE::_calc_empty_lattice_E_on_the_fly_loop(const MutInfo& mutInfo, string mode, double s) {
  Debug debugInfo("HB_EE::_calc_empty_lattice_E_on_the_fly_loop(const MutInfo& mutInfo, string mode, double s)");

  double total_E = 0;
  string base_method;  

  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  int CBCalc = 1; // default = 1: calculates ground spectrum with presence of CB on variable standard sidechains.  Does not affect H-Bond calculations.

  SCREAM_HB_BASE_FUNCTIONAL_OBJ* hb_functor, *zero_delta_functor;

  this->_figureOutHBPotentialMethods(mode, base_method,  nonPolar_H_calc, CBCalc, s, &hb_functor);
  zero_delta_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, 0);

  ScreamAtomV tmp_sc_atoms; tmp_sc_atoms.clear();
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);

  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) )                       // if ->second == NULL, natural AA.
    tmp_sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  else 
    tmp_sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);

  for (ScreamAtomVItr atom = tmp_sc_atoms.begin(); atom != tmp_sc_atoms.end(); ++atom) { // initialize sc_atom list, after checking if the sidechain is considered to be empty.
    if ( ((*atom)->flags & 0x2) == 0) { // moveable atom--to include if H___A or other positive value of atom->hb_da
      if ( (*atom)->hb_da == -600) {
	cout << "This should not come up! in on_the_fly_loop in hb_EE. 1. " << endl;
	map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*atom)->stripped_atomType );
	if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) {// non polar heavy atom 
	  (*atom)->hb_da = -1;
	}
	else 
	  (*atom)->hb_da = mStrInt_Itr->second;
      }
      if ( (*atom)->hb_da >=0 )
	sc_atoms.push_back(*atom);
    }
    
  }
  if (sc_atoms.size() == 0)
    return 0;
  //  cout << "Relevent sidechain atoms in empty lattice H-Bond calculations: " << sc_atoms.size() << endl;

  // then get all atoms on the protein. no need to take special care of NOCB case.
  ScreamAtomV * ptn_atom_list_ptr; 
  ScreamAtomV atom_list; // atom_list that are relevent to HBonding.
  //if (this->ON_THE_FLY ==1 ) {
  ptn_atom_list_ptr = &(this->ptn->getAtomList() ); ScreamAtomVConstItr ptn_atom_list_ptr_end = ptn_atom_list_ptr->end();
  ScreamAtomV self_bb_list;

  // Now organize atom_list.
  for (ScreamAtomVConstItr ptn_atom = ptn_atom_list_ptr->begin(); ptn_atom != ptn_atom_list_ptr_end; ++ptn_atom) {
    if ( ((*ptn_atom)->flags & 0x2) and ((*ptn_atom)->flags & 0x8) ) { // 0x2: fixed atoms 0x8: ones to be considered as visible in EL calculations (for flexibility).
      if ( (*ptn_atom)->hb_da == -600) {
	cout << "This should not come up! in on_the_fly_loop in hb_EE. 2. " << endl;
	map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*ptn_atom)->stripped_atomType );
	if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) // non polar heavy atom.
	  (*ptn_atom)->hb_da = -1;
	else {
	  (*ptn_atom)->hb_da = mStrInt_Itr->second;
	  
	}
      }
      if ((*ptn_atom)->resNum == mutInfo.pstn and (*ptn_atom)->chain == mutInfo.chn) { // self bb atoms
	if ( (*ptn_atom)->hb_da >= 0 )
	  self_bb_list.push_back(*ptn_atom);
      }
      else // all other atoms
	if ( (*ptn_atom)->hb_da >= 0 )
	  atom_list.push_back(*ptn_atom);
    }
  }

  double R_on = this->hb_obj->R_on;  double R_off = this->hb_obj->R_off;
  double theta_on = this->hb_obj->theta_on;  double theta_off = this->hb_obj->theta_off;

  SCREAM_ATOM *A, *H, *D;

  ScreamAtomVConstItr sc_atoms_end = sc_atoms.end();
  ScreamAtomVConstItr atom_list_end = atom_list.end();


  // First, do SC self-bb HBonds.  No need to check 1-3, 1-4 exclusion.
  for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms_end; ++sc_atom) {
    for (ScreamAtomVConstItr self_bb_atom = self_bb_list.begin(); self_bb_atom != self_bb_list.end(); ++self_bb_atom) {
      if ( (*sc_atom)->hb_da == 0) { // if sc atom is H___A
	if ( (*self_bb_atom)->hb_da == 0)
	  continue;
	else {
	  A = *self_bb_atom; 
	  H = *sc_atom;
	  D = H->connectivity_m.begin()->first;
	}
      }
      else { // if sc atom is heavy atom acceptor
	if ( (*self_bb_atom)->hb_da != 0) 
	  continue;
	else {
	  A = *sc_atom;
	  H = *self_bb_atom;
	  D = H->connectivity_m.begin()->first;
	}
	
      }
      total_E += (*zero_delta_functor)(A, H, D);
    }
  }

  /* Main double loop.*/

  for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms_end; ++sc_atom) {
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list_end; ++ptn_atom) {
      //if (*sc_atom == *ptn_atom) continue; // self. impossible.
      if ( (*sc_atom)->hb_da == 0 ) {
	if ( (*ptn_atom)->hb_da == 0) // if ptn atom is H___A
	  continue; // since both sc_atom and ptn atom are H___A
	else {
	  A = *ptn_atom;
	  H = *sc_atom;
	  D = H->connectivity_m.begin()->first; // H___A should be connected to only 1 atom... 
	}
      }
      else {
	if ( (*ptn_atom)->hb_da != 0) // if ptn atom is not H___A
	  continue; // since acceptor needs to accept H___A
	else {
	  A = *sc_atom;
	  H = *ptn_atom;
	  D = H->connectivity_m.begin()->first;
	}
      }

      // Case 0: arb lib cases, not optimized.
      if (mutInfo.chn == "Z") {
	if (! (scream_tools::should_exclude_on_1_2(A, D)
	       or scream_tools::should_exclude_on_1_3(A, D) ) ) {
	  total_E += (*hb_functor)(A, H, D);
	  //total_E += (*zero_delta_functor)(A,H,D);
	}
      }
      else
	  total_E += (*hb_functor)(A, H, D);
    }
  }

  //  clock_t t2 = clock();
    

  delete zero_delta_functor; delete hb_functor;
  return total_E;
}

double HB_EE::_calc_all_interaction_E_on_the_fly_loop(string mode, double r) {
  Debug debugInfo("HB_EE::_calc_all_interaction_E_on_the_fly_loop(string mode, double r)");

  double total_E = 0;
  string base_method;
  int nonPolar_H_calc = 1; // default = 1: calculates nonPolar H VDW's.
  int CBCalc = 1; // default = 1: calculates ground spectrum with presence of CB on variable standard sidechains.  doesn't matter for hb calculations

  SCREAM_HB_BASE_FUNCTIONAL_OBJ* hb_functor, *zero_delta_functor;

  this->_figureOutHBPotentialMethods(mode, base_method,  nonPolar_H_calc, CBCalc, r, &hb_functor);
  zero_delta_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, 0);

  if (this->_variable_atoms_on_each_sidechain.empty())
    this->setup_variableAtomsOnEachSidechain();

  /* Then, setup moveable atoms and fixed atoms.  Take care to only include residues and atoms that are HBond donors and acceptors.*/
  map<int, ScreamAtomV> fixed_mutInfo_atoms, moveable_mutInfo_atoms;
  map<int, ScreamAtomV>::iterator v_end = this->_variable_atoms_on_each_sidechain.end();
  ScreamAtomV fixed_polar_atoms, moveable_polar_atoms;

  for (map<int, ScreamAtomV>::iterator itr = this->_variable_atoms_on_each_sidechain.begin(); itr != v_end; ++itr) {
    int mutInfo_n = itr->first;
    ScreamAtomV * atoms = &(itr->second);

    for (ScreamAtomVConstItr atom = atoms->begin(); atom != atoms->end(); ++atom) {
      if ( ( (*atom)->flags & 0x4 ) == 0 ) // i.e. invisible
      	continue;
      if ( (*atom)->flags & 0x2)  { // i.e. fixed 
	if ( (*atom)->hb_da < -500) {
	  cout << "ERROR!  hb_da should not be -600 here! 1. " << endl;
	  (*atom)->dump();
	}
	
	if ( (*atom)->hb_da >= 0)  // i.e. donor or acceptor or H___A (H___A has hb_da value of 0)
	  fixed_polar_atoms.push_back(*atom);
      }
      else {
	if ( (*atom)->hb_da < -500) {
	  cout << "ERROR!  hb_da should not be -600 here! 2. " << endl;
	  (*atom)->dump();
	}
	if ( (*atom)->hb_da >= 0)
	  moveable_polar_atoms.push_back(*atom);
      }
    }
    if (fixed_polar_atoms.size() != 0)
      fixed_mutInfo_atoms[mutInfo_n] = fixed_polar_atoms;
    if (moveable_polar_atoms.size() != 0)
      moveable_mutInfo_atoms[mutInfo_n] = moveable_polar_atoms;
    fixed_polar_atoms.clear();
    moveable_polar_atoms.clear();

  }

  /* Then, first calculate moveable-moveable atom HB.*/
  SCREAM_ATOM *A(NULL), *H(NULL), *D(NULL);
  double dist_cutoff_sq = 6*6;
  map<int, ScreamAtomV>::iterator mm_end = moveable_mutInfo_atoms.end();
  for (map<int, ScreamAtomV>::iterator itr1 = moveable_mutInfo_atoms.begin(); itr1 != mm_end; ++itr1) {
    map<int, ScreamAtomV>::iterator  itr2 = itr1; ++itr2;
    for (; itr2 != mm_end; ++itr2) {
      ScreamAtomV* list1 = &(itr1->second); 
      ScreamAtomV* list2 = &(itr2->second);

      for (ScreamAtomVConstItr a1 = list1->begin(); a1 != list1->end(); ++a1)
	for (ScreamAtomVConstItr a2 = list2->begin(); a2 != list2->end(); ++a2) {
	  if ( (*a1)->hb_da == 0) { // if a1 is H___A
	    if ( (*a2)->hb_da == 0) continue;
	    else {
	      A = *a2;
	      H = *a1;
	      D = H->connectivity_m.begin()->first;
	    }
	  }
	  else {
	    if ( (*a2)->hb_da != 0) continue;
	    else {
	      A = *a1;
	      H = *a2;
	      D = H->connectivity_m.begin()->first;
	    }
	  }
	  // haven't put in distance cutoff yet.
	  total_E += (*hb_functor)(A,H,D);
	}
    }
  }
  
  /* Then, main loop; moveable-fixed atom HB. */
  map<int, ScreamAtomV>::iterator fm_end = fixed_mutInfo_atoms.end();
  for (map<int, ScreamAtomV>::iterator itr1 = moveable_mutInfo_atoms.begin(); itr1 != mm_end; ++itr1)
    for (map<int, ScreamAtomV>::iterator itr2 = fixed_mutInfo_atoms.begin(); itr2 != fm_end; ++itr2) {
      ScreamAtomV* list1 = &(itr1->second);      
      ScreamAtomV* list2 = &(itr2->second);
         
      for (ScreamAtomVConstItr a1 = list1->begin(); a1 != list1->end(); ++a1)
	for (ScreamAtomVConstItr a2 = list2->begin(); a2 != list2->end(); ++a2) {
	  if ( (*a1)->hb_da == 0) { // if a1 is H___A 

	    if ( (*a2)->hb_da == 0) continue;
	    else {
	      A = *a2;
	      H = *a1;
	      D = H->connectivity_m.begin()->first;
	    }
	  }
	  else {
	    if ( (*a2)->hb_da != 0) continue;
	    else {
	      A = *a1;
	      H = *a2;
	      D = H->connectivity_m.begin()->first;
	    }
	  }
	  // haven't put in distance cutoff yet.
	  total_E += (*hb_functor)(A,H,D);
	}
    }

  delete hb_functor; delete zero_delta_functor;
  return total_E;

}


void HB_EE::_figureOutHBPotentialMethods(const string method, string& base_method, int& nonPolar_H_calc, int& CBCalc, double s, SCREAM_HB_BASE_FUNCTIONAL_OBJ* *hb_functor) {
  stringV f;
  split(method, "_", f);
  base_method = f[0];

  /* Fields that needed to be searched: "ASYM", "NOCB", "NONONPOLARH". */
  for (int i=0; i<f.size(); ++i) {
    string s = f[i];
    if (s == "ASYM")
      base_method = base_method + "_ASYM";
    if (s == "NONONPOLARH")
      nonPolar_H_calc = 0;
    if (s == "NOCB")
      CBCalc = 0;
  }

  /* Initialize functors. */


  if (base_method == "FULL" or base_method == "FULL_ASYM") {
    *hb_functor = new SCREAM_calc_full_delta_HB(this->hb_obj, s);
  } else if (base_method == "FLAT" or base_method == "FLAT_ASYM") {
    *hb_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, s);
    //cout << "FLAT functor initiated" << endl;
  } else if (base_method == "SCALE") {
    // not implemented yet.
    cout << "SCALE functor not implemented" << endl;
  } else {
    cout << "no functors initiated" << endl;
  }


}


int HB_EE::potential_HB_acceptor(SCREAM_ATOM* atom) {
  string ff_label = atom->stripped_atomType;
  // ADAM: Added F_, Cl, Br, I_ as HB acceptor
  bool potential_acceptor = ((ff_label == "O_3") or 
			     (ff_label == "O_2") or
			     (ff_label == "O_R") or
			     (ff_label == "N_2") or
			     (ff_label == "N_3") or
			     (ff_label == "N_R") or 
			     (ff_label == "S_3") or
                             (ff_label == "F_")  or
                             (ff_label == "Cl")  or
                             (ff_label == "Br")  or
                             (ff_label == "I_"));

  if (potential_acceptor) { return 1; }
  else return 0;
  
  
}

double HB_EE::calc_EL_rot_selfBB(const MutInfo& mutInfo, std::string mode, double s) {
  Debug debugInfo("HB_EE::calc_EL_rot_selfBB(const MutInfo& mutInfo, string mode, double s)");

  double total_E = 0;
  SCREAM_HB_BASE_FUNCTIONAL_OBJ* hb_functor, *zero_delta_functor;
  if (mode == "FULL" or mode == "FULL_ASYM" or mode == "FULL_ASYM_NONONPOLAR") {
    hb_functor = new SCREAM_calc_full_delta_HB(this->hb_obj, s);
  } else if (mode == "FLAT" or mode == "FLAT_ASYM" or mode == "FLAT_ASYM_NONONPOLAR") {
    hb_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, s);
  } else if (mode == "SCALE") {
    cout << "SCALE functor not implemented" << endl;
  } else {
    cout << "no functors initiated because they're not implemented" << endl;
  }
  zero_delta_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, 0);
  // first get atoms from this sidechain.
  //ScreamAtomV sc_atoms = this->ptn->get_sc_atoms(mutInfo.chn, mutInfo.pstn);
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);


  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }


  // then get all atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();


  atom_list = this->ptn->getAtomList();

  // note: for ON_THE_FLY case 2, not going to enforce symmetry since the R_ON R_OFF implicitly enforces "sphere-like-ness".
  // also, there are so few HB bonding partners it is not going to be the bottleneck.  so saving myself coding time.
  
  // main double loop
 
  //debugInfo.out( "Protein atom list size: " + string(itoa(atom_list.size() )) );
  //debugInfo.out( "SC atom list size : " + string(itoa(sc_atoms.size() )) );

  double R_on = this->hb_obj->R_on;
  double R_off = this->hb_obj->R_off;
  double theta_on = this->hb_obj->theta_on;
  double theta_off = this->hb_obj->theta_off;

  SCREAM_ATOM *A, *H, *D;

  for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
    if ( (*sc_atom)->hb_da == -600) { // then assign hb_da.
      map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*sc_atom)->stripped_atomType );
      if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) { // non polar heavy atom.
	(*sc_atom)->hb_da = -1;
      }
      else {
	(*sc_atom)->hb_da = mStrInt_Itr->second;
      }
    }

    //if ( (*sc_atom)->atomType != "H___A" and (!(this->potential_HB_acceptor(*sc_atom))) ) continue;
    if ( (*sc_atom)->hb_da == -1 ) continue; // if current atom has nothing to do with HBond.  0: H___A.  other values, >1: base atom or acceptor.  -1: value if atom has nothing to do with such.

    // protein atom loop here!
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      if (*sc_atom == *ptn_atom) continue; // self.
      if ( (*ptn_atom)->hb_da == -600) {
	map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*ptn_atom)->stripped_atomType );
	if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) { // non polar heavy atom.
	  (*ptn_atom)->hb_da = -1;
	}
      }
      
      if ( (*sc_atom)->hb_da == 0 ) {
	A = NULL;
	H = *sc_atom;
	D = H->connectivity_m.begin()->first; // H___A should be connected to only 1 atom... 
      }
      else {
	A = *sc_atom;
	H = NULL;
	D = NULL;
      }

      // Define A, H and D.
      if (A == NULL) {
	//if (!(this->potential_HB_acceptor(*ptn_atom)) ) {
	if ( (*ptn_atom)->hb_da <= 0 ) { // if ptn atom has nothing to do with HBonding or it's an H___A
	  continue;
	} else {
	  A = *ptn_atom; 
	}
      }
      else if (A != NULL) {
	if ( (*ptn_atom)->hb_da != 0 ) { // if ptn atom is not H___A
	  continue;
	} else {
	  H = *ptn_atom;
	  D = H->connectivity_m.begin()->first;
	}
      }
      
      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 0: arb lib cases, not optimized.
      if (mutInfo.chn == "Z") {
	return 0; // self bb only
	if (! (scream_tools::should_exclude_on_1_2(A, D)
	       or scream_tools::should_exclude_on_1_3(A, D) ) ) {
	  total_E += (*zero_delta_functor)(A,H,D);
	}
      }

      // Case 1: looping over sidechain atoms in protein atom list.  
      else if (ptn_flag == 0) {
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {  // if two atoms on same sidechain
	  // This loop is actually never executed; no sidechain can make H-bond with self sidechain.
	  continue;
	}
	else {
	  // do nothing: empty lattice calculations.
	  continue;
	}
      } 
      // Case 2: fixed atoms (backbone atoms)
      else if (ptn_flag == 1) { // i.e. looping over fixed atoms.  additional loop because checking is slow.
	// This part: in testing.  Should sidechain make self backbone?  Running QM to test.
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
		 or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	    total_E += (*zero_delta_functor)(A,H,D); // uses zero-delta-functor
	  }
	}// end if backbone atom is same residue as siechain
      } // end ptn_flag == 1 block
      
    } // inner loop
  } // outer loop
  return total_E;


}

double HB_EE::calc_EL_rot_otherBB(const MutInfo& mutInfo, std::string mode, double s) {

  Debug debugInfo("HB_EE::calc_EL_rot_otherBB(const MutInfo& mutInfo, string mode, double s)");

  double total_E = 0;
  SCREAM_HB_BASE_FUNCTIONAL_OBJ* hb_functor, *zero_delta_functor;
  if (mode == "FULL" or mode == "FULL_ASYM" or mode == "FULL_ASYM_NONONPOLAR") {
    hb_functor = new SCREAM_calc_full_delta_HB(this->hb_obj, s);
  } else if (mode == "FLAT" or mode == "FLAT_ASYM" or mode == "FLAT_ASYM_NONONPOLAR") {
    hb_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, s);
  } else if (mode == "SCALE") {
    // not implemented yet.
    cout << "SCALE functor not implemented" << endl;
  } else {
    cout << "no functors initiated" << endl;
  }
  zero_delta_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, 0);
  // first get atoms from this sidechain.
  //ScreamAtomV sc_atoms = this->ptn->get_sc_atoms(mutInfo.chn, mutInfo.pstn);

  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);


  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }


  // then get all atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  atom_list = this->ptn->getAtomList();
  
  // note: for ON_THE_FLY case 2, not going to enforce symmetry since the R_ON R_OFF implicitly enforces "sphere-like-ness".
  // also, there are so few HB bonding partners it is not going to be the bottleneck.  so saving myself coding time.
  
  // main double loop
 
  //debugInfo.out( "Protein atom list size: " + string(itoa(atom_list.size() )) );
  //debugInfo.out( "SC atom list size : " + string(itoa(sc_atoms.size() )) );

  double R_on = this->hb_obj->R_on;
  double R_off = this->hb_obj->R_off;
  double theta_on = this->hb_obj->theta_on;
  double theta_off = this->hb_obj->theta_off;

  SCREAM_ATOM *A, *H, *D;

  for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
    if ( (*sc_atom)->hb_da == -600) { // then assign hb_da.
      map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*sc_atom)->stripped_atomType );
      if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) { // non polar heavy atom.
	(*sc_atom)->hb_da = -1;
      }
      else {
	(*sc_atom)->hb_da = mStrInt_Itr->second;
      }
    }

    //if ( (*sc_atom)->atomType != "H___A" and (!(this->potential_HB_acceptor(*sc_atom))) ) continue;
    if ( (*sc_atom)->hb_da == -1 ) continue; // if current atom has nothing to do with HBond.  0: H___A.  other values, >1: base atom or acceptor.  -1: value if atom has nothing to do with such.

    // protein atom loop here!
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      if (*sc_atom == *ptn_atom) continue; // self.
      if ( (*ptn_atom)->hb_da == -600) {
	map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*ptn_atom)->stripped_atomType );
	if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) { // non polar heavy atom.
	  (*ptn_atom)->hb_da = -1;
	}
      }
      

      if ( (*sc_atom)->hb_da == 0 ) {
	A = NULL;
	H = *sc_atom;
	D = H->connectivity_m.begin()->first; // H___A should be connected to only 1 atom... 
      }
      else {
	A = *sc_atom;
	H = NULL;
	D = NULL;
      }

      // Define A, H and D.
      if (A == NULL) {
	//if (!(this->potential_HB_acceptor(*ptn_atom)) ) {
	if ( (*ptn_atom)->hb_da <= 0 ) { // if ptn atom has nothing to do with HBonding or it's an H___A
	  continue;
	} else {
	  A = *ptn_atom; 
	}
      }
      else if (A != NULL) {
	//if ( (*ptn_atom)->atomType != "H___A") {
	if ( (*ptn_atom)->hb_da != 0 ) { // if ptn atom is not H___A
	  continue;
	} else {
	  H = *ptn_atom;
	  D = H->connectivity_m.begin()->first;
	}
      }
      
      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 0: arb lib cases, not optimized.
      if (mutInfo.chn == "Z") {
	return 0; // if "Z", return 0 for now.
	if (! (scream_tools::should_exclude_on_1_2(A, D)
	       or scream_tools::should_exclude_on_1_3(A, D) ) ) {
	  //total_E += this->hb_obj->calc_VDW_6_12( *ptn_atom, *sc_atom);   // then only calc energies if 1-2 and 1-3 exclusions are satisfied.
	  //total_E += (*hb_functor)(A, H, D);
	  total_E += (*zero_delta_functor)(A,H,D);
	}
	
      }

      // Case 1: looping over sidechain atoms in protein atom list.  
      else if (ptn_flag == 0) {
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {  // if two atoms on same sidechain
	  // This loop is actually never executed; no sidechain can make H-bond with self sidechain.  sometimes self-backbone, but no self sidechain.
	}
	else {
	  // do nothing: empty lattice calculations.
	}
      } 
      // Case 2: fixed atoms (backbone atoms)
      // OtherBB.
      else if (ptn_flag == 1) { // i.e. looping over fixed atoms.  additional loop because checking is slow.
	// This part: in testing.  Should sidechain make self backbone?  Running QM to test.
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  continue;
	}// end if backbone atom is same residue as siechain
	
	else {
	  if ( scream_tools::is_BB_atom((*ptn_atom)->stripped_atomLabel) ) {
	    total_E += (*hb_functor)(A, H, D);
	  }
	}
      } // end ptn_flag == 1 block
    } // inner loop
  } // outer loop
  return total_E;

}

double HB_EE::calc_EL_rot_fixedSC(const MutInfo& mutInfo, std::string mode, double s) {
  Debug debugInfo("HB_EE::calc_EL_rot_fixedSC(const MutInfo& mutInfo, string mode, double s)");

  double total_E = 0;
  SCREAM_HB_BASE_FUNCTIONAL_OBJ* hb_functor, *zero_delta_functor;
  if (mode == "FULL" or mode == "FULL_ASYM" or mode == "FULL_ASYM_NONONPOLAR") {
    hb_functor = new SCREAM_calc_full_delta_HB(this->hb_obj, s);
  } else if (mode == "FLAT" or mode == "FLAT_ASYM" or mode == "FLAT_ASYM_NONONPOLAR") {
    hb_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, s);
  } else if (mode == "SCALE") {
    // not implemented yet.
    cout << "SCALE functor not implemented" << endl;
  } else {
    cout << "no functors initiated" << endl;
  }
  zero_delta_functor = new SCREAM_calc_flat_delta_HB(this->hb_obj, 0);
  // first get atoms from this sidechain.
  //ScreamAtomV sc_atoms = this->ptn->get_sc_atoms(mutInfo.chn, mutInfo.pstn);

  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);


  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }


  // then get all atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();

  atom_list = this->ptn->getAtomList();
  
  // note: for ON_THE_FLY case 2, not going to enforce symmetry since the R_ON R_OFF implicitly enforces "sphere-like-ness".
  // also, there are so few HB bonding partners it is not going to be the bottleneck.  so saving myself coding time.
  
  // main double loop
 
  //debugInfo.out( "Protein atom list size: " + string(itoa(atom_list.size() )) );
  //debugInfo.out( "SC atom list size : " + string(itoa(sc_atoms.size() )) );

  double R_on = this->hb_obj->R_on;
  double R_off = this->hb_obj->R_off;
  double theta_on = this->hb_obj->theta_on;
  double theta_off = this->hb_obj->theta_off;

  SCREAM_ATOM *A, *H, *D;

  for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
    if ( (*sc_atom)->hb_da == -600) { // then assign hb_da.
      map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*sc_atom)->stripped_atomType );
      if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) { // non polar heavy atom.
	(*sc_atom)->hb_da = -1;
      }
      else {
	(*sc_atom)->hb_da = mStrInt_Itr->second;
      }
    }

    if ( (*sc_atom)->hb_da == -1 ) continue; // if current atom has nothing to do with HBond.  0: H___A.  other values, >1: base atom or acceptor.  -1: value if atom has nothing to do with such.

    // protein atom loop here!
    for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
      if (*sc_atom == *ptn_atom) continue; // self.
      if ( (*ptn_atom)->hb_da == -600) {
	map<string, int>::iterator mStrInt_Itr = this->hb_obj->hb_atom_type_mapping.find( (*ptn_atom)->stripped_atomType );
	if (mStrInt_Itr == this->hb_obj->hb_atom_type_mapping.end() ) { // non polar heavy atom.
	  (*ptn_atom)->hb_da = -1;
	}
      }
      
      if ( (*sc_atom)->hb_da == 0 ) {
	A = NULL;
	H = *sc_atom;
	D = H->connectivity_m.begin()->first; // H___A should be connected to only 1 atom... 
      }
      else {
	A = *sc_atom;
	H = NULL;
	D = NULL;
      }

      // Define A, H and D.
      if (A == NULL) {
	if ( (*ptn_atom)->hb_da <= 0 ) { // if ptn atom has nothing to do with HBonding or it's an H___A
	  continue;
	} else {
	  A = *ptn_atom; 
	}
      }
      else if (A != NULL) {
	if ( (*ptn_atom)->hb_da != 0 ) { // if ptn atom is not H___A
	  continue;
	} else {
	  H = *ptn_atom;
	  D = H->connectivity_m.begin()->first;
	}
      }
      
      int ptn_flag = (*ptn_atom)->flags & 0x2;
      // Case 0: arb lib cases, not optimized.
      if (mutInfo.chn == "Z") {
	return 0; // don't do this case for now.
      }

      // Case 1: looping over sidechain atoms in protein atom list.  
      else if (ptn_flag == 0) {
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {  // if two atoms on same sidechain
	  // This loop is actually never executed; no sidechain can make H-bond with self sidechain.  sometimes self-backbone, but no self sidechain.

	}
	else {
	  // do nothing: empty lattice calculations.
	}
      } 
      // Case 2: fixed atoms (backbone atoms)
      else if (ptn_flag == 1) { // i.e. looping over fixed atoms.  additional loop because checking is slow.
	// This part: in testing.  Should sidechain make self backbone?  Running QM to test.
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  continue;
	}// end if backbone atom is same residue as siechain
	
	else {
	  if ( scream_tools::is_SC_atom( (*ptn_atom)->stripped_atomLabel) ) {
	    total_E += (*hb_functor)(A, H, D);
	  }
	}
	
      } // end ptn_flag == 1 block
      
    } // inner loop
  } // outer loop
  return total_E;

}

void HB_EE::setup_variableAtomsOnEachSidechain() {
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
      // below: now handled by scream_EE initFixedMoveable routine in scream_EE.
      //if ( ! ( (*itr)->stripped_atomLabel == "CB" or (*itr)->stripped_atomLabel == "HCB" ) )
      //relevantAtoms.push_back(*itr);
    }
    //tmp_sc_atoms.clear();
    //tmp_sc_atoms.insert(tmp_sc_atoms.end(), relevantAtoms.begin(), relevantAtoms.end());
    
    this->_variable_atoms_on_each_sidechain[mutInfo_n] = relevantAtoms;
    this->_mutInfo_n_map[mutInfo_n] = &(mIrotC_itr->first); // will this work?
    ++mutInfo_n;
  }

}


double HB_EE::calc_EL_rot_fixedHET(const MutInfo& mutInfo, std::string mode, double s) {

}

double HB_EE::calc_EL_rot_moveableHET(const MutInfo& mutInfo, std::string mode, double s) {

}
