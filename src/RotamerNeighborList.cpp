#include "defs.hpp"
#include "RotamerNeighborList.hpp"
#include "scream_vector.hpp"
//vcvicek
#include <algorithm>


RotamerNeighborList::RotamerNeighborList() {
  this->ptn = NULL;
}

RotamerNeighborList::RotamerNeighborList(Protein* ptn, map<MutInfo, RotConnInfo*> mI_rCI_map, double cutoff) {
  this->ptn = ptn;
  //  this->mutInfo_rotConnInfo = mI_rCI_map;
  this->cutoff = cutoff;
  this->mutInfo_rotConnInfo.clear();

  vector<MutInfo> HJ_list; HJ_list.clear();
  for (map<MutInfo, RotConnInfo*>::const_iterator itr = mI_rCI_map.begin(); itr != mI_rCI_map.end(); itr++) {
    MutInfo mI = itr->first;
    // His and Hse exception.  Need to make two neighbor lists.
    if (mI.getAA() == "H" or mI.getAA() == "J") {
      string other_AA = "J"; // Usually, H is the default one, J is the "other" one.
      if (mI.getAA() == "J") other_AA = "H";
      string HJ_mI_string = other_AA + string(itoa(mI.getPstn())) + "_" + mI.getChn();
      MutInfo HJ_mI(HJ_mI_string);
      
      vector<MutInfo>::iterator HJ_itr = find(HJ_list.begin(), HJ_list.end(), mI);
      
      if (HJ_itr == HJ_list.end()) {
	HJ_list.push_back(HJ_mI);
      }
      else {
	HJ_list.erase(HJ_itr);
      }
      
    }
    this->mutInfo_rotConnInfo[mI] = itr->second;
  }
  
  for (vector<MutInfo>::const_iterator itr = HJ_list.begin(); itr != HJ_list.end(); itr++) {
    this->mutInfo_rotConnInfo[*itr] = NULL;
  }
  
//   for (map<MutInfo, RotConnInfo*>::const_iterator itr = this->mutInfo_rotConnInfo.begin(); itr != this->mutInfo_rotConnInfo.end(); itr++) {
//     //    cout << "Printing map<MutInfo, RotConnInfo*> in RotamerNEighborList" << endl;
//     //    cout << itr->first << endl;
//   }
  //cout << "before initRotamerNeighborList" << endl;
  this->initRotamerNeighborList();

}

void RotamerNeighborList::setProtein(Protein* ptn) {
  this->ptn = ptn;
}

Protein* RotamerNeighborList::getProtein() {
  return this->ptn;
}

void RotamerNeighborList::addMutInfoRotConnInfo(MutInfo mI, RotConnInfo* rCI) {

  this->mutInfo_rotConnInfo.insert(make_pair(mI, rCI));

}

void RotamerNeighborList::initRotamerNeighborList() {
  /* Sets up a rotamer neighbor list for each MutInfo, RotConnInfo pair. If MutInfo NtrlAA, ignore associated RotConninfo */
  /* State variables */
  this->cutoff; // cutoff distance
  this->ptn; // protein.
  this->mutInfo_rotConnInfo; // residue positions and residues defined
  /* Final result */
  this->emptyLatticeNeighborLists.clear();

  /* First, prepare new relevant atom list, exclude all atoms in <MutInfo, RotConnInfo> pair */

  ScreamAtomV emptyLatticeAtomList = this->_prepareEmptyLatticeAtomList(this->ptn, this->mutInfo_rotConnInfo);
  //  for_each (emptyLatticeAtomList.begin(), emptyLatticeAtomList.end(), dump);


  /* Then iterate through all residues and make neighboring atom lists */
  for (map<MutInfo, RotConnInfo*>::const_iterator itr = this->mutInfo_rotConnInfo.begin(); 
       itr != this->mutInfo_rotConnInfo.end(); itr++) {
    MutInfo mI = itr->first;
    RotConnInfo* rCI = itr->second;
    //cout << mI << endl;
    //    cout << "about to _initOneRotamerNeighborList" << endl;
    //    cout << " before _initOneRotamerNeighborList" << endl;
    this->emptyLatticeNeighborLists[mI] = this->_initOneRotamerNeighborList(this->cutoff, emptyLatticeAtomList, mI, rCI);
    // cout << "done " << mI << endl;
  }
  
  //cout << "Exiting initRotamerNeighborList" << endl;

}

ScreamAtomV RotamerNeighborList::_prepareEmptyLatticeAtomList(Protein* ptn, map<MutInfo, RotConnInfo*>& mutInfo_rotConnInfo) {
  cout << " Preparing EmptyLatticeAtomList! " << endl;
  ScreamAtomV emptyLatticeAtomList = ptn->getAtomList(); // emptyLatticeAtomList is: atoms that form the "background" in emptylattice calculations.
  //  cout << "emptyLatticeAtomList size: " << emptyLatticeAtomList.size() << endl;

  /* First, determine which atoms are to be excluded from emptyLatticeAtomList. */
  ScreamAtomV atomsToBeRemoved;
  for (map<MutInfo, RotConnInfo*>::const_iterator itr = mutInfo_rotConnInfo.begin();
       itr != mutInfo_rotConnInfo.end(); itr++) {
    MutInfo mI = itr->first;
    RotConnInfo* rCI = itr->second;

    if (rCI) { // canonical rCI != NULL case; before clustering happened.
      ScreamAtomV mI_sc_atoms = ptn->get_variable_atoms(rCI);
      atomsToBeRemoved.insert(atomsToBeRemoved.end(), mI_sc_atoms.begin(), mI_sc_atoms.end());
    }

    else {
      ScreamAtomV mI_sc_atoms = ptn->get_sc_atoms(mI);
      atomsToBeRemoved.insert(atomsToBeRemoved.end(), mI_sc_atoms.begin(), mI_sc_atoms.end());
    }
    
    // else { // If NtrlAA's.
//       string chn = mI.getChn();
//       int pstn = mI.getPstn();
//       ScreamAtomV mI_sc_atoms = ptn->get_sc_atoms(chn, pstn);
//       atomsToBeRemoved.insert(atomsToBeRemoved.end(), mI_sc_atoms.begin(), mI_sc_atoms.end());
//     }

  }

  /* Then, remove those atoms from ptn atomlist */
  
  for (ScreamAtomVItr itr = atomsToBeRemoved.begin(); itr != atomsToBeRemoved.end(); itr++) {
    ScreamAtomVItr to_be_removed_atom_i = find(emptyLatticeAtomList.begin(), emptyLatticeAtomList.end(), *itr);
    if (to_be_removed_atom_i != emptyLatticeAtomList.end() )
      emptyLatticeAtomList.erase(to_be_removed_atom_i);
  }

  return emptyLatticeAtomList;

}

ScreamAtomV RotamerNeighborList::_initOneRotamerNeighborList(double cutoff, ScreamAtomV& emptyLatticeAtomList, MutInfo mI, RotConnInfo* rCI) {

  /** Algorithm: Core idea: Decide a point to be the center of the sphere.  Build atom list around it.
   *   Which point though?  Possibilities: center of current residue, or center of average reside rotamer (need external rotamer library information).
   */

  /* Pseudocode */
  //cout << "Printing MI in _initOneRotamerNeighborList" << endl;
  //cout << mI << endl;
  //cout << "before _determineCenter " << endl;
  ScreamVector CENTER = this->_determineCenter(emptyLatticeAtomList, mI, rCI);
  //cout << "after _determineCenter " << endl;
  //ScreamVector CENTER =  ScreamVector();
  SCREAM_ATOM* CENTER_ATOM = new SCREAM_ATOM(); // helper atom
  CENTER_ATOM->x[0] = CENTER[0];
  CENTER_ATOM->x[1] = CENTER[1];
  CENTER_ATOM->x[2] = CENTER[2];

  double cutoff_sq = cutoff * cutoff;
  ScreamAtomV neighborList; neighborList.clear();

  //  cout << " now entering building neighborList loop" << endl;

  for (ScreamAtomVConstItr itr = emptyLatticeAtomList.begin(); itr != emptyLatticeAtomList.end(); itr++) 
    {
      SCREAM_ATOM* atom = *itr;
      if (CENTER_ATOM->distance_squared(atom) > cutoff_sq) {
	continue;
      } 
      else {
	// also need to check for 1-3 and 1-4 exclusion.  but seems like this is better left for scream_vdw_EE and the likes.
	neighborList.push_back(atom);
      }
    }

  // Need to include the sidechain atoms back.  This is because of the way energies are calculated in scream_vdw_EE etc; may need to think about the optimal way to code this down the road.  No solid conclusion at this point.

  //ScreamAtomV sc_atoms = this->ptn->get_sc_atoms(mI.chn, mI.pstn);
  ScreamAtomV sc_atoms;
  if (rCI == NULL) {
    sc_atoms = this->ptn->get_sc_atoms(mI);
  }
  else {
    sc_atoms = this->ptn->get_variable_atoms(rCI);
  }

  neighborList.insert(neighborList.end(), sc_atoms.begin(), sc_atoms.end() );
  return neighborList;

}

ScreamAtomV& RotamerNeighborList::returnEmptyLatticeNeighborList(MutInfo mI) {

  return this->emptyLatticeNeighborLists[mI];

}

ScreamVector RotamerNeighborList::_determineCenter(ScreamAtomV& eLAtomList, MutInfo mI, RotConnInfo* rCI) {

  return ScreamVector(0,0,0);

  // need to rewrite this; need to incorporate tree-structure MutInfo's.
  if (rCI) {
    // Working on this! Need Protein information to be able to easily extract atoms from RotConnInfo.
    return ScreamVector(0,0,0);
  }
  else {

    string chn; 
    int pstn;

    if (mI.isClusterMutInfo()) {
      // for now, just pick one of the MutInfo's from list.  Need to improve on this.
      vector<MutInfo*> mI_list = mI.getAllMutInfos();
      chn = mI_list[0]->getChn();
      pstn = mI_list[0]->getPstn();
    }
    else {
      chn = mI.chn;
      pstn = mI.pstn;
    }
    
    
    for (ScreamAtomVConstItr itr = eLAtomList.begin(); itr != eLAtomList.end(); itr++) {
      string atom_chn = (*itr)->chain;
      int atom_pstn = (*itr)->resNum;
      if (chn == atom_chn && pstn == atom_pstn ) {
	if ( (*itr)->stripped_atomLabel == "CA" ) {
	  ScreamVector CA( (*itr)->x[0], (*itr)->x[1], (*itr)->x[2] );
	  return CA;
	}
      }
    }
  }

}
