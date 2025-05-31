#include "defs.hpp"
#include "MutInfo.hpp"
#include "scream_coulomb_EE.hpp"
#include <algorithm>

#include <time.h>
#include <stdio.h>
#include <cassert>

#include "RotamerNeighborList.hpp"

Coulomb_EE::Coulomb_EE() {
  this->rotamerNeighborList = NULL;
}

Coulomb_EE::Coulomb_EE(Protein* ptn, vector<MutInfo> mutInfoV, SCREAM_Coulomb_OBJ* coulomb_obj_) : coulomb_obj(coulomb_obj_) {
  this->_initVariableAndFixedAtomPairList(ptn, mutInfoV);
  this->_initVariableAndVariableAtomPairList(ptn, mutInfoV);

  this->ptn = ptn;

}

Coulomb_EE::~Coulomb_EE() {

}

void Coulomb_EE::init_after_addedMutInfoRotConnInfo(Protein* ptn, SCREAM_Coulomb_OBJ* coulomb_obj_) {
  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->coulomb_obj = coulomb_obj_;

  this->_initVariableAndFixedAtomPairListArb(ptn, this->mutInfo_rotConnInfo_map);
  this->_initVariableAndVariableAtomPairListArb(ptn, this->mutInfo_rotConnInfo_map);

  this->ptn = ptn;
  this->ON_THE_FLY = 0;

}

void Coulomb_EE::init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein* ptn, SCREAM_Coulomb_OBJ* coulomb_obj_) {

  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->coulomb_obj = coulomb_obj_;

  //  this->_initFixedMoveableAtomsOnProtein(ptn, this->mutInfo_rotConnInfo_map);

  this->ptn = ptn;
  this->ON_THE_FLY = 1;

}

void Coulomb_EE::init_after_addedMutInfoRotConnInfo_neighbor_list(Protein* ptn, SCREAM_Coulomb_OBJ* coulomb_obj_, RotamerNeighborList* rNL) {

  assert(this->mutInfo_rotConnInfo_map.size() != 0);
  this->coulomb_obj = coulomb_obj_;

  this->ptn = ptn;
  this->ON_THE_FLY = 2;
  this->rotamerNeighborList = rNL;

}

void Coulomb_EE::addMutInfoRotConnInfo(MutInfo mutInfo, RotConnInfo* rotConnInfo = NULL) {
  /* adds a mutInfo, rotconninfo pair to mutInfo_rotConnInfo_map */
  this->mutInfo_rotConnInfo_map[mutInfo] = rotConnInfo;
}


double Coulomb_EE::calc_empty_lattice_E(const MutInfo mutInfo) {
  double total_E = 0;
  //map<const MutInfo, vector<AtomPair> >::iterator itr = this->variable_and_fixed[mutInfo];
  vector<AtomPair> atomPairList(this->variable_and_fixed[mutInfo]);
  //  cout << "coulomb EE calc size: " << atomPairList.size() << endl;
  //  clock_t t1 = clock();
  //  cout << atomPairList[0].a1->q[0] << endl;
  //  cout << atomPairList[0].a2->q[0] << endl;
  

  if (! this->ON_THE_FLY) {
    //    cout << "Total number of coulomb pairs: " << atomPairList.size() << endl;
    for (vector<AtomPair>::iterator ap_i = atomPairList.begin(); ap_i != atomPairList.end(); ++ap_i) {
      total_E += this->coulomb_obj->calc_Coulomb( (*ap_i).a1, (*ap_i).a2 );
    }
    
  }
  else {
    total_E += this->_calc_empty_lattice_E_on_the_fly_loop(mutInfo);

  }
  
  //  clock_t t2 = clock();
  //  cout << "that took: " << (double)(t2 - t1) / (double) CLOCKS_PER_SEC << endl;
  return total_E;

}

double Coulomb_EE::calc_residue_interaction_E(const MutInfo mI) {
   // calculates interaction between MutInfo residue and rest of the variable atoms.
  double total_E = 0;
  for (map< MutInfoPair, vector<AtomPair> >::iterator itr = this->variable_and_variable.begin();
       itr != this->variable_and_variable.end(); ++itr) {
    if ( mI == (itr->first).mutInfo1 or mI == (itr->first).mutInfo2 ) {
      for (vector<AtomPair>::iterator ap_i = (itr->second).begin(); ap_i != (itr->second).end(); ++ap_i) {
	total_E += this->coulomb_obj->calc_Coulomb( (*ap_i).a1, (*ap_i).a2 );
      }
    }
  }
    
  return total_E;

}

double Coulomb_EE::calc_residue_interaction_E(const MutInfo mI1, const MutInfo mI2) {

  double total_E = 0;
  //  int i = 0;
  //  int double_count  = 0;
  // first get atoms from this sidechain.
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

  
  for (ScreamAtomVConstItr mI1_atom_itr = mI1_atoms.begin(); mI1_atom_itr != mI1_atoms.end(); ++mI1_atom_itr) {
    for (ScreamAtomVConstItr mI2_atom_itr = mI2_atoms.begin(); mI2_atom_itr != mI2_atoms.end(); ++mI2_atom_itr) {

      SCREAM_ATOM* mI1_atom = *mI1_atom_itr;
      SCREAM_ATOM* mI2_atom = *mI2_atom_itr;

      if (! (scream_tools::should_exclude_on_1_2(mI1_atom, mI2_atom)
	       or scream_tools::should_exclude_on_1_3(mI1_atom, mI2_atom) ) ) {
	double E = this->coulomb_obj->calc_Coulomb(mI1_atom, mI2_atom);
	total_E += E;
	
      }
    }
  }

  return total_E;

  

}

double Coulomb_EE::calc_all_interaction_E() {
  // calculates interaction energies between all atom pairs
  double total_E = 0;

  if (! ON_THE_FLY) {

    for (map< MutInfoPair, vector<AtomPair> >::iterator itr = this->variable_and_variable.begin(); 
	 itr != this->variable_and_variable.end(); ++itr) {
      for (vector<AtomPair>::iterator ap_i = (itr->second).begin(); ap_i != (itr->second).end(); ++ap_i) {
	total_E += this->coulomb_obj->calc_Coulomb( (*ap_i).a1, (*ap_i).a2 );
      }
    }
  }
  else {
    total_E += this->_calc_all_interaction_E_on_the_fly_loop();

  }
  return total_E;
}

double Coulomb_EE::calc_all_interaction_E_delta() {
  double total_E = 0;
  total_E = this->calc_all_interaction_E();
  return total_E;
}

double Coulomb_EE::calc_EL_rot_selfBB(const MutInfo& mutInfo) {
  Debug debugInfo("Coulomb_EE::calc_EL_rot_selfBB(const MutInfo& mutInfo)");

  double total_E = 0;

  //  int i = 0;
  //  int double_count  = 0;
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
  
  
  // then get all atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  atom_list = this->ptn->getAtomList();

  double cutoff_on_sq = 8.5*8.5;
  double cutoff_off_sq = 10.5*10.5;
  double cutoff_diff = cutoff_off_sq - cutoff_on_sq;

  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }
  if (mutInfoChn_is_Z) {
    return 0;
  }

  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      int ptn_flag = ( (*ptn_atom)->flags & 0x2 );
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      // Case 1: moveable atoms.
      if (ptn_flag == 0) {
	continue;
      }
      // Case 2: fixed atoms.  (i.e. non-sidechain atoms)
      else { // if ptn_flag ==1 
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  // Comment Jul-29-06: put it back in; when when mutation isn't an issue, don't need to worry.
	  // The Following 1-2 and 1-3 exclusion necessary; sidechain can interact with self backbone.  if FF accurate no need to exclude sidechain backbone self interaction; i.e. GLN and ASN self-H-bond cases shouldn't come up (here, self strong dipoles).
	  string ptnAtomLabel = scream_tools::strip_whitespace((*ptn_atom)->getAtomLabel());
	  // Note: make the following not hardcoded at some point.
	  if ( ptnAtomLabel != "O" and 
	       ptnAtomLabel != "HN" and
	       ptnAtomLabel != "OXT") {
	    continue;
	  }
	  if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
		 or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	    double coul_E = this->coulomb_obj->calc_Coulomb(*ptn_atom, *sc_atom);
	    if (coul_E > 100) {
	      debugInfo.out("coul_E > 100" );
	      debugInfo.out( stringify(coul_E) );
	      debugInfo.out( (*ptn_atom)->return_bgf_line() );
	      debugInfo.out( (*sc_atom)->return_bgf_line() );
	    }
	    
	    total_E += coul_E;
	  }
	}// end if backbone atom is same residue as siechain
	
      } // end inner loop
    } // end outer loop
  }
  return total_E;
  
}
   
double Coulomb_EE::calc_EL_rot_otherBB(const MutInfo& mutInfo) {
  Debug debugInfo("Coulomb_EE::calc_EL_rot_otherBB(const MutInfo& mutInfo)");
     
  double total_E = 0;

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

  double cutoff_on_sq = 8.5*8.5;
  double cutoff_off_sq = 10.5*10.5;
  double cutoff_diff = cutoff_off_sq - cutoff_on_sq;

  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }
  if (mutInfoChn_is_Z) {
    return 0;
  }

   for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      int ptn_flag = ( (*ptn_atom)->flags & 0x2 );
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      if (mutInfoChn_is_Z) {
	return 0;
      }
      // Case 1: moveable atoms.
      else if (ptn_flag == 0) {
	continue;
      }
      // Case 2: fixed atoms.  (i.e. non-sidechain atoms)
      else { // if ptn_flag ==1 
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  continue;
	}// end if backbone atom is same residue as siechain
	// other BB interactions
	else { // no checking necessary; guaranteed to be at least 3 valence bonds away.  also functional groups _almost_ not being cutoff; even if so protected by spline and the fact they're far away.
	  if ( scream_tools::is_SC_atom((*ptn_atom)->getAtomLabel())) {
	    continue;
	  }
	  double distance_sq = (*sc_atom)->distance_squared(*ptn_atom);
	  if (distance_sq > cutoff_off_sq) {
	    continue;
	  }
	  else {
	    // electrostatic cutoffs.
	    double coul_E = this->coulomb_obj->calc_Coulomb(*ptn_atom, *sc_atom); // no doulbe counting here, no need for /2.
	    if (distance_sq > cutoff_on_sq) {
	      // use simple spline: linear interpolation (on squared distances) between R_on and R_off on Coul_E(R_ij).
	      double spline_value = this->_linearSpline(cutoff_on_sq, cutoff_diff, distance_sq);
	      coul_E *= spline_value;
	    }

	    if (coul_E > 100) {
	      debugInfo.out("coul_E > 100" );
	      debugInfo.out( stringify(coul_E) );
	      debugInfo.out( (*ptn_atom)->return_bgf_line() );
	      debugInfo.out( (*sc_atom)->return_bgf_line() );
	    }
	    
	    total_E += coul_E;
	  }
	}
      }
    } // end inner loop
  } // end outer loop

  return total_E;

}

double Coulomb_EE::calc_EL_rot_fixedSC(const MutInfo& mutInfo) {
  Debug debugInfo("Coulomb_EE::calc_EL_rot_fixedSC(const MutInfo& mutInfo)");

  double total_E = 0;

  //  int i = 0;
  //  int double_count  = 0;
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

  // then get all atoms on the protein.
  ScreamAtomV atom_list; atom_list.clear();
  atom_list = this->ptn->getAtomList();

  double cutoff_on_sq = 8.5*8.5;
  double cutoff_off_sq = 10.5*10.5;
  double cutoff_diff = cutoff_off_sq - cutoff_on_sq;

  //  main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }
  if (mutInfoChn_is_Z) {
    return 0;
  }

  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list.end(); ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms.end(); ++sc_atom) {
      int ptn_flag = ( (*ptn_atom)->flags & 0x2 );
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      // Case 1: moveable atoms.
      if (ptn_flag == 0) {
	continue;
      }
      // Case 2: fixed atoms.  (i.e. non-sidechain atoms)
      else { // if ptn_flag ==1 
	if ( (*ptn_atom)->chain == mutInfo.chn and (*ptn_atom)->resNum == mutInfo.pstn) {
	  continue;
	}
	if ( scream_tools::is_BB_atom((*ptn_atom)->getAtomLabel()) ) {
	  continue;
	}
	double distance_sq = (*sc_atom)->distance_squared(*ptn_atom);
	if (distance_sq > cutoff_off_sq) {
	  continue;
	}
	else {
	  // electrostatic cutoffs.
	  double coul_E = this->coulomb_obj->calc_Coulomb(*ptn_atom, *sc_atom); // no doulbe counting here, no need for /2.
	  if (distance_sq > cutoff_on_sq) {
	    // use simple spline: linear interpolation (on squared distances) between R_on and R_off on Coul_E(R_ij).
	    double spline_value = 1- (distance_sq - cutoff_on_sq)/cutoff_diff; // = 1-0 when dist = cutoff_on; = 1-1 when dist = cutoff_off
	    coul_E *= spline_value;
	  }
	  if (coul_E > 100) {
	    debugInfo.out("coul_E > 100" );
	    debugInfo.out( stringify(coul_E) );
	    debugInfo.out( (*ptn_atom)->return_bgf_line() );
	    debugInfo.out( (*sc_atom)->return_bgf_line() );
	  }
	  
	  total_E += coul_E;
	}
	
	
      }
      
    } // end inner loop
  } // end outer loop
  
  return total_E;

}

double Coulomb_EE::calc_EL_rot_fixedHET(const MutInfo& mutInfo) {
  return 0;
}

double Coulomb_EE::calc_EL_rot_moveableHET(const MutInfo& mutInfo) {
  return 0;
}


void Coulomb_EE::setup_variableAtomsOnEachSidechain() {
  // Sets up this->_variable_atoms_on_each_sidechain.

  this->_variable_atoms_on_each_sidechain.clear();

  int mutInfo_n = 1; // index count
  // First set up variable_atoms_on_each_sidechain map.
  map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = this->mutInfo_rotConnInfo_map.begin();
  for (; mIrotC_itr != this->mutInfo_rotConnInfo_map.end(); ++mIrotC_itr) {
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;
    ScreamAtomV tmp_sc_atoms;
    ScreamAtomV relevantAtoms;
    if (rotConnInfo == NULL)
      relevantAtoms = ptn->get_sc_atoms(mutInfo);
    else
      // need to get variable atoms 
      relevantAtoms = ptn->get_variable_atoms(rotConnInfo);

    for (ScreamAtomVItr itr = relevantAtoms.begin(); itr != relevantAtoms.end(); ++itr) {
      if ( (*itr)->flags & 0x1 == 1) continue; // CB atom setup handled by scream_EE _initFixedMoveableAtomsOnProtein routine.
      else 
	tmp_sc_atoms.push_back(*itr);
    //tmp_sc_atoms.push_back(*itr); // put all sc atoms on the list
    }
    
    this->_variable_atoms_on_each_sidechain[mutInfo_n] = tmp_sc_atoms;
    ++mutInfo_n;
  }

}

void Coulomb_EE::_initVariableAndFixedAtomPairList(Protein* ptn, const vector<MutInfo> mutInfoV) {
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

  //for_each(all_variable_atoms.begin(), all_variable_atoms.end(), dump);

  /* initialize fixed atoms.  all non-variable atoms are fixed atoms. */
  ScreamAtomV ptn_atom_list = ptn->getAtomList();
  fixed_atoms.insert(fixed_atoms.end(), ptn_atom_list.begin(), ptn_atom_list.end());
  for (ScreamAtomVItr itr = all_variable_atoms.begin(); itr != all_variable_atoms.end(); ++itr) {
    ScreamAtomVItr variable_atom_i = find(fixed_atoms.begin(), fixed_atoms.end(), *itr);
    fixed_atoms.erase(variable_atom_i);
  }
  //for_each(fixed_atoms.begin(), fixed_atoms.end(), dump);


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


void Coulomb_EE::_initVariableAndFixedAtomPairListArb(Protein* ptn, const map<MutInfo, RotConnInfo*> mIrotC_map) {
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
      //this_sc_atoms = ptn->get_sc_atoms(chn, pstn);
      this_sc_atoms = ptn->get_sc_atoms(mutInfo);

    } else {
      this_sc_atoms = ptn->get_variable_atoms(rotConnInfo);

    }
    /* first, initialize interaction between sidechain and fixed atoms. */
    for (ScreamAtomVItr Var_i = this_sc_atoms.begin(); Var_i != this_sc_atoms.end(); ++Var_i) {
      for (ScreamAtomVItr Fixed_i = fixed_atoms.begin(); Fixed_i != fixed_atoms.end(); ++Fixed_i) {
	// 1-2 and 1-3 exclusion
	if (! (scream_tools::should_exclude_on_1_2(*Var_i, *Fixed_i)
	       or scream_tools::should_exclude_on_1_3(*Var_i, *Fixed_i) ) )
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




void Coulomb_EE::_initVariableAndVariableAtomPairList(Protein* ptn, const vector<MutInfo> mutInfoV) {
  
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

  for (map<MutInfoPair, vector<AtomPair> >::iterator itr = this->variable_and_variable.begin(); itr != this->variable_and_variable.end(); ++itr) {
    cout << (*itr).first.mutInfo1 << " " << (*itr).first.mutInfo2 << " size: " << (*itr).second.size() ;
  }

}



void Coulomb_EE::_initVariableAndVariableAtomPairListArb(Protein* ptn, const map<MutInfo, RotConnInfo*> mIrotC_map) {
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
      int pstn = mutInfo.pstn;
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


double Coulomb_EE::_calc_empty_lattice_E_on_the_fly_loop(const MutInfo& mutInfo) {
  Debug debugInfo("Coulomb_EE::_calc_empty_lattice_E_on_the_fly_loop(const MutInfo& mutInfo)");

  double total_E = 0;

  ScreamAtomV tmp_sc_atoms; tmp_sc_atoms.clear();
  ScreamAtomV sc_atoms; sc_atoms.clear();
  map<MutInfo, RotConnInfo*>::const_iterator mI_rCI_itr = this->mutInfo_rotConnInfo_map.find(mutInfo);

  if (mI_rCI_itr == this->mutInfo_rotConnInfo_map.end() or   // if AminoAcid name does not match... happens in cases of mutation.
      (mI_rCI_itr->second == NULL) ) {                       // if ->second == NULL, natural AA.
    tmp_sc_atoms = this->ptn->get_sc_atoms(mutInfo);  
  }
  else {
    tmp_sc_atoms = this->ptn->get_variable_atoms(mI_rCI_itr->second);
  }

  for (ScreamAtomVItr itr = tmp_sc_atoms.begin(); itr != tmp_sc_atoms.end(); ++itr)  // initialize sc_atom list, after checking if the sidechain is considered to be empty.
    if ( ((*itr)->flags & 0x2) == 0 ) // 0x2: if atom is moveable (CB is fixed sometimes) 
      sc_atoms.push_back(*itr);
  
  if (sc_atoms.size() == 0)
    return 0;

  /* Then get atoms on protein.  Organize. */
  ScreamAtomV atom_list; atom_list.clear();
  ScreamAtomV * ptn_atom_list_ptr = &(this->ptn->getAtomList()); 
  ScreamAtomVConstItr ptn_atom_list_ptr_end = ptn_atom_list_ptr->end();
  ScreamAtomV self_bb_list;

  for (ScreamAtomVConstItr ptn_atom = ptn_atom_list_ptr->begin(); ptn_atom != ptn_atom_list_ptr_end; ++ptn_atom) {
    if ( ((*ptn_atom)->flags & 0x2) and ((*ptn_atom)->flags & 0x8)  )  { // 0x2: fixed/moveable; if moveable: is some other sc.  0x8: EL invisibility.  for more flexible sc-fixed calculations.
      if ( (*ptn_atom)->resNum == mutInfo.pstn and (*ptn_atom)->chain == mutInfo.chn) // self bb atoms
	self_bb_list.push_back(*ptn_atom);
      else
	atom_list.push_back(*ptn_atom);

    }
  }

  /* Calculate self-bb SC electrostatics. */  
  ScreamAtomVConstItr atom_list_end = atom_list.end();
  ScreamAtomVConstItr sc_atoms_end = sc_atoms.end();

  for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms_end; ++sc_atom)
    for (ScreamAtomVConstItr self_bb_atom = self_bb_list.begin(); self_bb_atom != self_bb_list.end(); ++self_bb_atom)
      if (! (scream_tools::should_exclude_on_1_2(*self_bb_atom, *sc_atom) 
	     or scream_tools::should_exclude_on_1_3(*self_bb_atom, *sc_atom) ) ) {
	string selfBBAtomLabel = scream_tools::strip_whitespace((*self_bb_atom)->getAtomLabel());
	// O, HN and OXT only self-bb atoms whose electrostatic interaction with sidechain not included in rotamer library.
	if ( selfBBAtomLabel != "O" and 
	     selfBBAtomLabel != "HN" and
	     selfBBAtomLabel != "OXT") {
	  continue;
	}
	total_E += this->coulomb_obj->calc_Coulomb(*self_bb_atom, *sc_atom);
      }

  /* Now main double loop. */
  double cutoff_on_sq = 8.5*8.5;
  double cutoff_off_sq = 10.5*10.5;
  double cutoff_diff = cutoff_off_sq - cutoff_on_sq;

  /* Main double loop below.  first pull out expensive string comparison stuff from double loop.  then main double loop.  All atom_list atoms are fixed atoms (i.e. backbone protein atoms or other fixed sidechains). */
  bool mutInfoChn_is_Z = false;
  if (mutInfo.chn == "Z") {   mutInfoChn_is_Z = true;  }
  
  for (ScreamAtomVConstItr ptn_atom = atom_list.begin(); ptn_atom != atom_list_end; ++ptn_atom) {
    for (ScreamAtomVConstItr sc_atom = sc_atoms.begin(); sc_atom != sc_atoms_end; ++sc_atom) {
      // Case 0: arblib atoms.  if arblib atom, currently no optimization.
      if (mutInfoChn_is_Z) {
	// Include electrostatics for arbitrary rotamer library guys.
	if (*sc_atom == *ptn_atom) continue; // should never happen.
	if (! (scream_tools::should_exclude_on_1_2(*ptn_atom, *sc_atom)
	       or scream_tools::should_exclude_on_1_3(*ptn_atom, *sc_atom) ) ) {
	  double dist_sq = (*sc_atom)->distance_squared(*ptn_atom);
	  if (dist_sq > cutoff_off_sq) {
	    continue;
	  } else {
	    double coul_E = this->coulomb_obj->calc_Coulomb(*ptn_atom, *sc_atom); // division by 2 not needed; variables never included in atom_list.
	    if (dist_sq > cutoff_on_sq) {
	      // use simple spline: linear interpolation (on squared distances) between R_on and R_off on Coul_E(R_ij).
	      double spline_value = 1- (dist_sq - cutoff_on_sq)/cutoff_diff; // = 1-0 when dist = cutoff_on; = 1-1 when dist = cutoff_off
	      coul_E *= spline_value;
	    }
	    total_E += coul_E;
	  } 
	}
      }

      else { 
	double distance_sq = (*sc_atom)->distance_squared(*ptn_atom);
	if (distance_sq > cutoff_off_sq) {
	  continue;
	}
	else {
	  // electrostatic cutoffs.
	  double coul_E = this->coulomb_obj->calc_Coulomb(*ptn_atom, *sc_atom); // no doulbe counting here, no need for /2.
	  if (distance_sq > cutoff_on_sq) { 
	    // use simple spline: linear interpolation (on squared distances) between R_on and R_off on Coul_E(R_ij).
	    double spline_value = 1- (distance_sq - cutoff_on_sq)/cutoff_diff; // = 1-0 when dist = cutoff_on; = 1-1 when dist = cutoff_off
	    coul_E *= spline_value;
	  }
	  
	  total_E += coul_E;
	}
	
      }
      
    } // end inner loop
  } // end outer loop
  
  return total_E;
}

double Coulomb_EE::_calc_all_interaction_E_on_the_fly_loop() {
  /* initializes interaction between variable atoms and variable atoms only; for scream */
  double total_E = 0;

  if (this->_variable_atoms_on_each_sidechain.empty())
    this->setup_variableAtomsOnEachSidechain();

  /* Then setup fixed atoms and non-fixed atoms. */
  ScreamAtomV fixed_atoms, moveable_atoms;
  map<int, ScreamAtomV>::iterator v_end = this->_variable_atoms_on_each_sidechain.end();
  for (map<int, ScreamAtomV>::iterator itr = this->_variable_atoms_on_each_sidechain.begin(); itr != v_end; ++itr) {
    ScreamAtomV * atoms = &(itr->second);
   
    for (ScreamAtomVConstItr atom = atoms->begin(); atom != atoms->end(); ++atom) {
      if ( ( (*atom)->flags & 0x4 ) == 0) // i.e. invisible
	continue; 
      if ( (*atom)->flags & 0x2)  // i.e. fixed
	fixed_atoms.push_back(*atom);
      else {
	(*atom)->a = itr->first; // stores temporarily variable--mutInfo_n assignment
	moveable_atoms.push_back(*atom);
      }
    }
  }
  
  /* Setup spline stuff. */
  double cutoff_on_sq = 8.5*8.5;
  double cutoff_off_sq = 10.5*10.5;
  double cutoff_diff = cutoff_off_sq - cutoff_on_sq;


  /* First do moveable-moveable loop.*/
  ScreamAtomVConstItr fixed_atoms_end=fixed_atoms.end();
  ScreamAtomVConstItr moveable_atoms_end=moveable_atoms.end();
  for (ScreamAtomVConstItr a1=moveable_atoms.begin(); a1 != moveable_atoms_end; ++a1) {
    ScreamAtomVConstItr a2=a1;
    ++a2;
    for (; a2 != moveable_atoms_end; ++a2) {
      if ( (*a1)->a != (*a2)->a ) { // calculate interaction only if a1 and a2 has differenct mutInfo_n assignment
	double coulombE = 0;
	double dist_sq = (*a1)->distance_squared(*a2);
	if (dist_sq < cutoff_off_sq) {
	  coulombE = this->coulomb_obj->calc_Coulomb(*a1, *a2);
	  if (dist_sq > cutoff_on_sq) 
	    coulombE *= this->_linearSpline(cutoff_on_sq, cutoff_diff, dist_sq);
	  total_E += coulombE;
	}
      }
    }
  }
    
  /* Then do main loop. */
  for (ScreamAtomVConstItr mv_atom = moveable_atoms.begin(); mv_atom != moveable_atoms_end; ++mv_atom)
    for (ScreamAtomVConstItr fx_atom = fixed_atoms.begin(); fx_atom != fixed_atoms_end; ++fx_atom) {
      double coulombE = 0;
      double dist_sq = (*mv_atom)->distance_squared(*fx_atom);
      if (dist_sq < cutoff_off_sq)
	coulombE = this->coulomb_obj->calc_Coulomb(*mv_atom, *fx_atom);
      if (dist_sq > cutoff_on_sq) 
	coulombE *= this->_linearSpline(cutoff_on_sq, cutoff_diff, dist_sq);
      total_E += coulombE;
    }
  
  return total_E;
}

double Coulomb_EE::_linearSpline(double cutoff_on_sq, double cutoff_diff_sq, double dist_sq) {
  // returns a ratio.
  return ( 1- (dist_sq - cutoff_on_sq)/cutoff_diff_sq ); // = 1-0 when dist = cutoff_on; = 1-1 when dist = cutoff_off

}
