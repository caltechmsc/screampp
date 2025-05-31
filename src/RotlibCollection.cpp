#include "RotlibCollection.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <time.h>

RotlibCollection::RotlibCollection() {

  // Set default values.
  HIGHEST_ALLOWED_ROTAMER_E = 350;
  //vcvicek ELEnergyTablePstn = NULL;
  this->_mutInfo_mapping.clear();
  this->_mutInfo_inverse_mapping.clear();
}

RotlibCollection::~RotlibCollection() {
  
  // doesn't do anything because Rotlib* 's should take care of themselves.
}

void RotlibCollection::addRotlib(const string mutInfo, Rotlib* rotlib) {
  ////this->RotamerLibraryMap.insert(make_pair(mutInfo, rotlib));
  int crntMutInfoMapping_size = this->_mutInfo_mapping.size();
  if ( this->_mutInfo_mapping.find(mutInfo) == this->_mutInfo_mapping.end() ) {
    this->_mutInfo_mapping[mutInfo] = crntMutInfoMapping_size + 1; // start from 1, not zero.  offset one.  arbitrary, seems like a good idea at the time.
    this->_mutInfo_inverse_mapping[crntMutInfoMapping_size + 1] = mutInfo;
    this->RotamerLibraryMap.insert(make_pair( crntMutInfoMapping_size +1, rotlib) );
  }
  else {
    cout << " Input mutInfo string cooresponds to a mutInfo that has already been added!  Rotamer library not added." << endl;
    exit(2);
  }

}

void RotlibCollection::addClashCollection(ClashCollection* cc) {

  this->clashCollection = cc;

}

void RotlibCollection::cleanClashCollection() {

  this->clashCollection = NULL;

}

void RotlibCollection::initEmptyLatticeDataStructures() {
  /* 8-07-06: This routine: after not having been used for ages, probably wrong */

  /* Comment: there smarter ways to do this.  See dynamic memory routine for one example of such.  This is now obsolete. */

  cout << "initializing Data Structures for Empty Lattice applications..." << endl;

  /* Sort by empty lattice E */
  this->_sortAllRotlibByEmptyLatticeEnergy();
  cout << "Done sorting All Rotamer Libraries by their Empty Lattice Energies!" << endl;

  /* Then PRE-ORDER EVERYTHING: theoretically very inefficient.  But in practice probably okay because you probably only have about 5 rotamers max at each site.  Ok when you have 200 residues we're gonna run into trouble.  I'll write the dynamic programming implementation some time down the road. */
  
  /* First determine cutoff. */
  double CUTOFF = this->HIGHEST_ALLOWED_ROTAMER_E;
  if (CUTOFF <= 0) {
    cout << "WARNING: CUTOFF (i.e. HIGHEST_ALLOWED_ROTAMER_E) is less than or equal to zero!" << endl;
  }
  
  /* Now enumerate all combinations allowed */
  map<int, Rotlib*>::iterator RotlibItr;
  set< ExcitationEnumeration_n > tempSetOfExcitationLevels;

  for (RotlibItr = RotamerLibraryMap.begin(); RotlibItr != RotamerLibraryMap.end(); RotlibItr++) {
    /* preparation */
    int mutInfo_n = RotlibItr->first;
    Rotlib* currentRotlib = RotlibItr->second;
    assert(currentRotlib != NULL);

    int RotamerInEmptyLatticeEnvironmentRank = -1;
    bool enoughRotamerSamplingFlag = false;  
    currentRotlib->reset_pstn();  // not necessary; already done previously.  just to be safe.

    Rotamer* aaRot = currentRotlib->get_next_rot_with_empty_lattice_E_below(CUTOFF);
    
    //set< ExcitationEnumeration_n > tempSetOfExcitationLevels = this->EnumerationList; // with n residue pstns on each ExcitationEnumeration_n.

    
    if (aaRot != NULL) {	// if aaRot == NULL, don't clear EnumerationList!  otherwise it will stay empty forever because the loop underneath it never gets executed.
      this->EnumerationList.clear();
    }
    else { continue;}
    
    /* Looping through all rotamers in a Rotlib */
    
    while (aaRot != NULL and
	   enoughRotamerSamplingFlag == false) {
      
      RotamerInEmptyLatticeEnvironmentRank = aaRot->get_empty_lattice_energy_rank();
      /* if first time through outer loop; base case */
      if (tempSetOfExcitationLevels.empty()) {
	////ExcitationEnumeration tempEE;
	ExcitationEnumeration_n tempEE;
	tempEE[mutInfo_n] = RotamerInEmptyLatticeEnvironmentRank;	// now 1 residue on each ExcitationEnumeration
	this->EnumerationList.insert(tempEE); 
      }
      /* if not first time through outer loop */
      else {
	//for (vector<ExcitationEnumeration>::const_iterator ExcitationLevelItr = tempVOfExcitationLevels.begin(); 
	for (set<ExcitationEnumeration_n >::const_iterator ExcitationLevelItr = tempSetOfExcitationLevels.begin(); 
	     ExcitationLevelItr != tempSetOfExcitationLevels.end() ;
	     ExcitationLevelItr++ ) {
	  ExcitationEnumeration_n tempEE = *ExcitationLevelItr;
	  tempEE[mutInfo_n] = RotamerInEmptyLatticeEnvironmentRank; // now n+1 residue one each ExcitationEnumeration
	  this->EnumerationList.insert(tempEE);
	}
      }
      
      aaRot = currentRotlib->get_next_rot_with_empty_lattice_E_below(CUTOFF);
      // Remark: this i will forget.
      // if CLASHes, empty_lattice_E set at 999999999, therefore will never be returned a certain cutoff value.
    } // end of one Rotlib loop

  } // end of outer iterating through all rotlib loop
  cout << "Number of total enumerations: " << this->EnumerationList.size() << endl;
  /* Initializes ELEnumerationToEnergy */
  this->_calcEmptyLatticeLinearSumEnergy();
  cout << "Done calculating linear sum of empty lattice energies!" << endl;
  
  /* Initializes ELEnergyToEnumeration */
  this->_copyELEnumerationToEnergyToEnergyMultimap();  
  //cout << "Done _copyELEnumerationToEnergyToEnergyMultimap!" << endl;
    
  /* Initializes Current Position Pointer */
  this->ELEnergyTablePstn = ELEnergyToEnumeration.begin();
  
  cout << "Done initializing Data Structures for Empty Lattice Calculation!" << endl;
    

}


void RotlibCollection::initDynamicMemoryDataStructures() {
  /* Remark: this is the dynamic memory datat structure initialization setup.  The members used in this scheme is very much alike the schemes used by initEmptyLatticeDataStructures(), but be careful not to call both initializations.  Design-wise, the init's warrant two subclasses of a more general class.  But there's no time for such elegance, unfortunately.  */

  cout << "Initializing Data Structure for Dynamic Memory Data Structures. " << endl;

  /* Sort by empty lattice E */
  this->_sortAllRotlibByEmptyLatticeEnergy();
  cout << "Done sorting All Rotamer Libraries by their Empty Lattice Energies!" << endl;

  /* Calculate max number of rotamer configurations allowed after clash test. */
  this->maxRotamerConfigurations = this->_calcMaxNoClashRotamerConfigurations();
  cout << "Maximum number of excited rotamer configurations: " << this->maxRotamerConfigurations << endl;

  /* cutoff.  is this obsolete? */
  double CUTOFF = this->HIGHEST_ALLOWED_ROTAMER_E;
  if (CUTOFF <= 0) {
    cout << "WARNING: CUTOFF (i.e. HIGHEST_ALLOWED_ROTAMER_E) is less than or equal to zero!" << endl;
  }

  /* Initialize only the ground state. */
  ////map<string, Rotlib*>::iterator RotlibItr;
  map<int, Rotlib*>::iterator RotlibItr;
  ExcitationEnumeration_n tempEE;
  for (RotlibItr = RotamerLibraryMap.begin(); RotlibItr != RotamerLibraryMap.end(); RotlibItr++) {
    /* preparation */
    ////string mutInfo = string(RotlibItr->first);
    int mutInfo_n = RotlibItr->first;
    Rotlib* currentRotlib = RotlibItr->second;
    assert(currentRotlib != NULL);

    currentRotlib->reset_pstn();  // not necessary; already done previously.

    Rotamer* aaRot = currentRotlib->get_next_rot_with_empty_lattice_E_below(CUTOFF);

    /* Take just the ground state rotamer from each Rotlib, store in EnumerationList. */
    tempEE[mutInfo_n] = 0; // ground state is guaranteed to be 0th rank.

  }
  this->EnumerationList.insert(tempEE);
  
  /* Now that EnumerationList has been initialized, populate other relevent info. */
  this->_calcEmptyLatticeLinearSumEnergy();
  this->_copyELEnumerationToEnergyToEnergyMultimap(); 
  this->_initMethod = "DynamicMemory";

  /* Pointer initialization.  But since pointer is no longer valid in an ever-expanding multi-map structure, simply use a rank. */
  
  this->currentEmptyLatticeExcitationN = 0;  // 0 is ground state, 1 is 1st excitation, so on.
  
  cout << "Done initializing Data structures for Dynamic Memory Energy Excitation Scheme Calculation." << endl;

}

void RotlibCollection::initAllocationUnderEnergyThreshold(double thresholdE) {
  /* 8-07-06 */
  /* Function: initializes data structure relevent to using the "allocation under energy threshold" rotamer allocation method.

     Description: The DynamicMemoryDataStructures makes no assumption about underlying energies for the singles, and as a result, requires a lot more space then is necessary.  This function uses the information for singles energies and builds all configurations under a threshold excited energies.  Since information about the singles spectrum and a threshold is provided, this scheme requires a lot less memory space.  

     Usage: when run out of excited configurations, one needs to increase the threshold and expand the list of configurations.

     Implementation details: Uses the following:
     
       multimap< double, ExcitationEnumeration_n > ELEnergyToEnumeration;
  */

  cout << "Initializing Data Structure for Allocation Under Energy Threshold Structures. " << endl;

  /* First sort. */
  
  this->EnumerationList.clear();
  this->ELEnergyToEnumeration.clear();
  this->ELEnumerationToEnergy.clear();
  this->_sortAllRotlibByEmptyLatticeEnergy();

  /* Init maximum rotamer configurations. */
  this->maxRotamerConfigurations = this->_calcMaxNoClashRotamerConfigurations();
  cout << "From AllocationUnderThreshold engine: maximum number of excited rotamer configurations: " << this->maxRotamerConfigurations << endl;

  /* Then define variables for building ELEnergyToEnumeration. */

  /* Pseudocode.  Builds lexiconigraphically, but addition of alphabets at all positions must be less than thresholdE.
  /* Variables: sizeOfRotlib, tmpBelowThresholdConfigs, thresholdE

    initialize ground state (0,0,...0)
    for i (1..sizeOfRotlib): // i = position
      maxN_i = maximum index of Rotlib i such that the singles energy is under thresholdE 
      for j (1..maxN_i): // j = excitation (note: starts from 1, not 0)
        energy_ij = energy of j'th excited rotamer of i'th Rotlib 
	availableE = thresholdE - energyN_j
        build(tmpBelowThresholdConfigs, availableE, i,j)
    
    routine build(tmpBelowThresholdConfigs, maximumE, i,j):
      if i = 1: // first Rotlib
        allowedSubConfigs = NULL // all zeros, the ground state, i.e. (0,0,0...0)
      else:
        allowedSubConfigs = all combinatorial configs in (i-1, i-2, ... 1) positions with energy < maximumE

      for subConf in allowedSubConfigs:
        newConfig = attach (i,j) to allowedSubConfigs
	tmpBelowThresholdConfigs.add(newConfig)  // take care of the 0's as well.

  */

  multimap< double, ExcitationEnumeration_n > allConfigs; allConfigs.clear();

  this->_buildAllUnderThresholdConfigs(allConfigs, thresholdE);
  //  cout << "Done _buildAllUnderThresholdConfigs! " << endl;
  // Now Copy allConfigs back to where it belongs by a swap.
  this->ELEnergyToEnumeration.swap(allConfigs);
  this->currentEmptyLatticeExcitationN = 0;
  this->ELEnergyTablePstn = this->ELEnergyToEnumeration.begin();
  this->currentThresholdEnergy = thresholdE; 

  // Now populate other data structures in this class.
  this->_copyEnergyMultimapToELEnumerationToEnergy(); 
  this->_copyEnergyMultimapToELEnumerationList();
  this->_initMethod = "AllocationUnderThreshold";

  cout << " Data structure initialized for allocation under threshold energy for a value of " << thresholdE << " kcal/mol.  A total of " << this->ELEnergyToEnumeration.size() << " configurations initialized." << endl;
  

}

ExcitedRotamers RotlibCollection::getNextRotamersByELEnergy() {
  if (this->_initMethod == "DynamicMemory") {
    return this->getNextDynamicMemoryRotamers_And_Expand();
  }
  else if (this->_initMethod == "AllocationUnderThreshold") {
    return this->getNextUnderEnergyThresholdRotamers();
  }


}

void RotlibCollection::_buildAllUnderThresholdConfigs(multimap< double, ExcitationEnumeration_n > & allConfigs, 
						      double thresholdE) {
  /* Builds up all configurations under thresholdE, and stores it in the passed in by reference variable allConfigs.
   *
   */

  //  Debug debugInfo("RotlibCollection::_buildAllUnderThresholdConfigs");
  // Variable init.
  int numberOfRotlibs = this->_mutInfo_mapping.size();
  multimap< double, ExcitationEnumeration_n > newConfigs; newConfigs.clear();
  
  // Ground State init.
  ExcitationEnumeration_n groundState;
  for (unsigned short i = 1; i <= numberOfRotlibs; i++) { // remember: offset 1.
    groundState.insert(make_pair(i,0));
  }

  typedef multimap< double, ExcitationEnumeration_n >::value_type eeType ;
  allConfigs.insert(eeType(0,groundState)); // value 0: ground state energy is zero.  

  // Main loops.
  double availableE = thresholdE - 0; // ground state: 0.

  for (int i = 1; i <= numberOfRotlibs; i++) {
    Rotlib* rotlib_i = this->RotamerLibraryMap.find(i)->second;
    int maxN_i = rotlib_i->n_rotamers_below_empty_lattice_energy(thresholdE);
    for (int j = 1; j <  maxN_i; j++) { // start from excitation 1, nt zero.  strictly less than: maxN_i counts ground state as well.
      double energy_ij = 0;
      Rotamer* rot = rotlib_i->get_empty_lattice_E_rot_after_sorted_by_empty_lattice_E(j);
      if (rot == NULL) {
	continue;
      }
      energy_ij = rot->get_empty_lattice_E();

      availableE = thresholdE - energy_ij;

      this->_buildOneSetUnderThresholdConfigs(allConfigs, newConfigs, availableE, i,j, energy_ij);

    }
    // Done with previousRoundConfigs.
    allConfigs.insert(newConfigs.begin(), newConfigs.end());
    newConfigs.clear();
    
  }
}


void RotlibCollection::resetEmptyLatticeCrntPstn() {
  assert(ELEnergyToEnumeration.size() != 0);
  this->ELEnergyTablePstn = ELEnergyToEnumeration.begin();
}

void RotlibCollection::resetTotalEnergyCrntPstn() {
  assert(TEEnergyToEnumeration.size() != 0);
  this->TEEnergyTablePstn = TEEnergyToEnumeration.begin();
}


ExcitedRotamers RotlibCollection::getNextEmptyLatticeExcitationRotamers() {
  /* This routine (3/13/04) extracts from member variable set<ExcitationEnumeration> EnumerationList and returns a ExcitedRotamers object (a map<std::string, AARotamer*> ). */
  assert(ELEnergyToEnumeration.size() != 0);
  assert(ELEnergyToEnumeration.size() == ELEnumerationToEnergy.size());

  ExcitedRotamers toBeReturned;
  if (this->ELEnergyTablePstn == ELEnergyToEnumeration.end()) {
    cout << "Reached end of ELEnergyToEnumeration! " << endl;
    assert (toBeReturned.size() == 0);
    return toBeReturned;


  } else {
    toBeReturned = this->getELExcitedRotamerFromEnumeration(this->ELEnergyTablePstn->second);
    ELEnergyTablePstn++;
  } 

  return toBeReturned;

}

ExcitedRotamers RotlibCollection::getNextTotalEnergyExcitationRotamers() {

  assert(TEEnergyToEnumeration.size() != 0);
  //  assert(ELEnergyToEnumeration.size() == ELEnumerationToEnergy.size());
  ExcitedRotamers toBeReturned;

  if (this->TEEnergyTablePstn == TEEnergyToEnumeration.end()) {
    cout << "Reached end of TEEnergyToEnumeration! " << endl;
    assert (toBeReturned.size() == 0);
    return toBeReturned;


  } else {
    toBeReturned = this->getELExcitedRotamerFromEnumeration(TEEnergyTablePstn->second);
    TEEnergyTablePstn++;
  } 

  return toBeReturned;

}
  

ExcitedRotamers RotlibCollection::getNthEmptyLatticeExcitationRotamers() {

}

////ExcitedRotamers RotlibCollection::getELExcitedRotamerFromEnumeration(ExcitationEnumeration& Enumeration)  {
ExcitedRotamers RotlibCollection::getELExcitedRotamerFromEnumeration(ExcitationEnumeration_n& Enumeration)  {
  /* Basically, transforms one representation into another. */
  /* Objects querried: RotamerLibraryMap. Looks at the Rotamer.get_empty_lattice_energy_rank() */
  
  bool nullRotamerCheck = false;
  //ExcitedRotamers_n toReturnExcitedRotamers_n;
  ExcitedRotamers toReturnExcitedRotamers;

  //ExcitationEnumeration::const_iterator EEItr = Enumeration.begin();
  ExcitationEnumeration_n::const_iterator EEItr = Enumeration.begin();
  for (; EEItr != Enumeration.end(); EEItr++) {
    //const string mutInfo = EEItr->first;
    unsigned short mutInfo_n = EEItr->first;
    const unsigned short Rank = EEItr->second;
    //    NtrlAARotlib* rotlib = RotamerLibraryMap[mutInfo];
    Rotlib* rotlib = RotamerLibraryMap[mutInfo_n];
    //AARotamer* aaRot = rotlib->get_empty_lattice_E_rot(Rank);
    Rotamer* aaRot = rotlib->get_empty_lattice_E_rot(Rank);

    //cout << "in getELExcitedRotamerFromEnumeration: " << Rank << " " << mutInfo << endl;
    

    if (aaRot == NULL) {
      nullRotamerCheck = true;
      cerr << "nullRotamerCheck failed, in ExcitedRotamers RotlibCollection::getELExcitedRotamerFromEnumeration(ExcitationEnumeration Enumeration" << endl;
      cerr << "This means that a non-existent rank for a rotamer is specified.  Usually this means there are fewer rotamers in the rotamer library than the number specified." << endl;
    }
    toReturnExcitedRotamers[ this->_mutInfo_inverse_mapping[mutInfo_n] ] = aaRot;
    //toReturnExcitedRotamers_n[mutInfo_n] = aaRot;
  }

  return toReturnExcitedRotamers;

}

ExcitedRotamers_n RotlibCollection::getELExcitedRotamer_nFromEnumeration_n(ExcitationEnumeration_n& Enumeration)  {
  /* Basically, transforms one representation into another. */
  /* Objects querried: RotamerLibraryMap. Looks at the Rotamer.get_empty_lattice_energy_rank() */
  
  bool nullRotamerCheck = false;
  ExcitedRotamers_n toReturnExcitedRotamers_n;

  //ExcitationEnumeration::const_iterator EEItr = Enumeration.begin();
  ExcitationEnumeration_n::const_iterator EEItr = Enumeration.begin();
  for (; EEItr != Enumeration.end(); EEItr++) {
    //const string mutInfo = EEItr->first;
    unsigned short mutInfo_n = EEItr->first;
    const unsigned short Rank = EEItr->second;
    //    NtrlAARotlib* rotlib = RotamerLibraryMap[mutInfo];
    Rotlib* rotlib = RotamerLibraryMap[mutInfo_n];
    //AARotamer* aaRot = rotlib->get_empty_lattice_E_rot(Rank);
    Rotamer* aaRot = rotlib->get_empty_lattice_E_rot(Rank);

    //cout << "in getELExcitedRotamerFromEnumeration: " << Rank << " " << mutInfo << endl;
    

    if (aaRot == NULL) {
      nullRotamerCheck = true;
      cerr << "nullRotamerCheck failed, in ExcitedRotamers RotlibCollection::getELExcitedRotamerFromEnumeration(ExcitationEnumeration Enumeration" << endl;
      cerr << "This means that a non-existent rank for a rotamer is specified.  Usually this means there are fewer rotamers in the rotamer library than the number specified." << endl;
    }
    //toReturnExcitedRotamers[ this->_mutInfo_inverse_mapping[mutInfo_n] ] = aaRot;
    toReturnExcitedRotamers_n[mutInfo_n] = aaRot;
  }

  return toReturnExcitedRotamers_n;

}



ExcitationEnumeration RotlibCollection::getELEnumerationFromExcitedRotamer(ExcitedRotamers& ER) {

  /* Transforms ExcitedRotamer representation to ExcitationEnumeration representation. */
  /* Objects querried: looks at rotamer.empty_lattice_energy_rank */
  

  ExcitationEnumeration tempEE;
  ExcitedRotamers::iterator ER_itr = ER.begin();
  for (; ER_itr != ER.end(); ER_itr++) {

    Rotamer* thisRot = ER_itr->second;
    assert(thisRot != NULL);
    
    string mutInfo = ER_itr->first;
    unsigned short rotELRank = thisRot->get_empty_lattice_energy_rank();
    
    tempEE[mutInfo] = rotELRank;
    
  }

  return tempEE;

}

ExcitationEnumeration_n RotlibCollection::getELEnumeration_nFromExcitedRotamer_n(ExcitedRotamers_n& ER) {

  ExcitationEnumeration_n tempEE;
  ExcitedRotamers_n::iterator ER_itr = ER.begin();
  for (; ER_itr != ER.end(); ER_itr++) {

    Rotamer* thisRot = ER_itr->second;
    assert(thisRot != NULL);

    unsigned short mutInfo_n = ER_itr->first;
    unsigned short rotELRank = thisRot->get_empty_lattice_energy_rank();
    
    tempEE[mutInfo_n] = rotELRank;

  }

  return tempEE;

}

ExcitationEnumeration_n RotlibCollection::_ExcitationEnumerationToExcitationEnumeration_n(ExcitationEnumeration& EE) {

  ExcitationEnumeration_n tempEE_n;
  ExcitationEnumeration::iterator EE_itr = EE.begin();
  for (; EE_itr != EE.end(); EE_itr++) {
    
    unsigned short excitation = EE_itr->second;
    string mutInfo = EE_itr->first;
    unsigned short mutInfo_n = this->_mutInfo_mapping[mutInfo];

    tempEE_n[mutInfo_n] = excitation;

  }
  
  return tempEE_n;
}

ExcitedRotamers RotlibCollection::_ExcitedRotamers_nToExcitedRotamers(ExcitedRotamers_n& ER_n) {

  ExcitedRotamers ER;
  ExcitedRotamers_n::iterator ER_n_itr = ER_n.begin();
  for (; ER_n_itr != ER_n.end(); ER_n_itr++) {

    Rotamer* thisRot = ER_n_itr->second;
    assert(thisRot != NULL);

    unsigned short mutInfo_n = ER_n_itr->first;
    string mutInfo = this->_mutInfo_inverse_mapping[mutInfo_n];

    ER[mutInfo] = thisRot;

  }

  return ER;

}

ClashCollection* RotlibCollection::getClashCollection() {

  return this->clashCollection;

}

ExcitedRotamers RotlibCollection::getNextDynamicMemoryRotamers_And_Expand() {
  /* This routine (8/19/04, start date) */
  /* Modification, add latest_round_EnumerationList and latest_round_ELEnumerationToEnergy time-saving helper structures. (8/25/05) */
  
  //  float t1 = (float) clock();
  ExcitedRotamers_n nextRotamers_n = this->_getNextRotamerDynamicMemory();
  //  cout << "after this->_getNextRotamerDynamicMemory()" << endl;

  //  float t2 = (float) clock();
  ExcitationEnumeration_n nextEnumeration = this->getELEnumeration_nFromExcitedRotamer_n(nextRotamers_n);
  //  cout << "after this->getELEnumeration_nFromExcitedRotamer_n(nextRotamers_n)"<< endl;

  //  float t3 = (float) clock();

  if (nextEnumeration.size() == 0) {
    // end of all possible rotamer configurations that passed clash test 
    cout << "Reached end of all possible rotamer configurations that pass individual rotamer clash test! " << endl;
    nextRotamers_n.clear();
    // do nothing
  } else {

    //    t4 = (float) clock();
    this->_expandDynamicMemory(nextEnumeration);
    //    t5 = (float) clock();
    // Increment.
    this->currentEmptyLatticeExcitationN++;

  }

//   cout << "Timing in getNextDynamicMemoryRotamers_And_Expand() " << endl;
//   cout << "Timing: _getNextRotamerDynamicMemory: " << (t2-t1)/(float) CLOCKS_PER_SEC << endl;
//   cout << "Timing: getELEnumerationFromExcitedRotamer: " << (t3-t2)/(float) CLOCKS_PER_SEC << endl;
//   cout << "Timing: in between... " << (t4-t3)/(float) CLOCKS_PER_SEC << endl;
//   cout << "Timing: _expandDynamicMemory: " << (t5-t4)/(float) CLOCKS_PER_SEC << endl;
//   cout << "TIming end." << endl;
  ExcitedRotamers nextRotamers = this->_ExcitedRotamers_nToExcitedRotamers(nextRotamers_n);
  
  return nextRotamers;
}


ExcitedRotamers RotlibCollection::getNextDynamicClashEliminatedRotamers_And_Expand() {

  //  float t1 = (float) clock() ;
  
  ExcitedRotamers eR = this->getNextDynamicMemoryRotamers_And_Expand();

  //  float t2 = (float) clock(); 

  //  cout << " Time for this->getNextDynamicMemoryRotamers_And_Expand(): " << (t2-t1)/(float)CLOCKS_PER_SEC << endl;

  if (eR.size() == 0 ) {
    return eR; ///< Means there are no more dnamic memory rotamers to be returned--reached end of all possible rotamer configurations.
  }

  ExcitationEnumeration eE = this->getELEnumerationFromExcitedRotamer(eR);
  int clash_eliminated_n = 0;
  assert(this->clashCollection);
  int clash_status = this->clashCollection->checkClash(eE);

  int clash_count = 0;
  int max_clash_count = 1000;
  while (clash_status == 1) {
    // if clash, iterate until reaches one without clashes.
    cout << "There is clash!" << endl;
    clash_eliminated_n++;
    this->clashCollection->increment_total_clashing_rotamers_eliminated();
    eR = this->getNextDynamicMemoryRotamers_And_Expand();
    eE = this->getELEnumerationFromExcitedRotamer(eR);
    clash_status = this->clashCollection->checkClash(eE);

    clash_count++;
    if (clash_count > max_clash_count) {
      cout << "Max clash count, 1000, reached!  Returning empty/NULL rotamer configuration. " << endl;
      eR.clear();
      eE.clear();
    }
  } 

  cout << " Eliminated " << clash_eliminated_n << " set of rotamer configurations. " << endl;
  this->clashCollection->storeCurrentRotamerConfiguration(eE);
  return eR;

}

void RotlibCollection::increaseConfigurationsUnderEnergyThreshold(double extraE) {

}

ExcitedRotamers RotlibCollection::getNextUnderEnergyThresholdRotamers() {

  //  Debug debugInfo("RotlibCollection::getNextUnderEnergyThresholdRotamers");

  ExcitedRotamers ER; ER.clear();
  //cout << " in getNextUnderEnergyThresholdRotamers " << endl;
  if (this->ELEnergyTablePstn == this->ELEnergyToEnumeration.end()) {
    if (this->ELEnergyToEnumeration.size() >= this->maxRotamerConfigurations) {
      cout << " All rotamer configurations visited! Returning empty ExcitedRotamer." << endl;
    } 
    else {
      // Expand--can't think of easy way, just rebuild the whole damn thing.
      cout << "Expanding allocated enumerations!" << endl;
      
      int sizeOfSystem = this->_mutInfo_mapping.size();
      int crntPstn = this->ELEnergyToEnumeration.size();
      int beforeExpandSize = crntPstn;
      double crntEnergy = this->ELEnergyToEnumeration.rbegin()->first;
      double multiplier = pow(3, (float) (1/ (float)sizeOfSystem) );
      double newEnergy = this->currentThresholdEnergy * multiplier; // approximately three times as much stuff allocated, i.e. can last approx twice as long as previously.
      int afterExpandSize = beforeExpandSize;
      
      while (afterExpandSize == beforeExpandSize) {

	this->ELEnergyToEnumeration.clear();
	this->ELEnumerationToEnergy.clear();
	// Rebuild--not writing special expand function.
	this->_buildAllUnderThresholdConfigs(this->ELEnergyToEnumeration, newEnergy);
	
	afterExpandSize = this->ELEnergyToEnumeration.size();
	newEnergy *= multiplier;


      }
	// Find lowest energy above previous threshold.
	//this->ELEnergyTablePstn = this->ELEnergyToEnumeration.lower_bound(crntEnergy);
      int nthPstn = this->ELEnergyToEnumeration.size() - crntPstn;
      std::multimap<double, ExcitationEnumeration_n >::iterator nthItr = this->ELEnergyToEnumeration.begin();
      for (int c = 1; c <= crntPstn; c++) {
	nthItr++;
      }

      ExcitationEnumeration_n ER_n = nthItr->second;
      ER = this->getELExcitedRotamerFromEnumeration(ER_n);

      // Misc state vars to update.
      this->currentEmptyLatticeExcitationN = crntPstn+1;
      this->currentThresholdEnergy = newEnergy; 
      nthItr++;
      this->ELEnergyTablePstn = nthItr;
      this->_copyEnergyMultimapToELEnumerationToEnergy();
    }
  }
  else {
    //cout << "size of this->ELEnergyToEnumeration " << this->ELEnergyToEnumeration.size() << endl;
    //    debugInfo.out(" Size of this->ELEnergyToEnumeration " + string(itoa(this->ELEnergyToEnumeration.size() ) ) );
    ExcitationEnumeration_n ER_n = this->ELEnergyTablePstn->second;
    this->ELEnergyTablePstn++;

    ER = this->getELExcitedRotamerFromEnumeration(ER_n);

  }
  return ER;
}


map<std::string, NtrlAARotlib*> RotlibCollection::getMutInfoRotlibMap() {

}

map<std::string, NtrlAARotlib*> RotlibCollection::getMutInfoRotlibDict() {

  //  return RotamerLibraryMap;

}

void RotlibCollection::setExcitationEnergy(ExcitationEnumeration EE, double energy) {
  //  assert(TEEnumerationToEnergy.size() == TEEnergyToEnumeration.size() );

  
  ExcitationEnumeration_n EE_n = this->_ExcitationEnumerationToExcitationEnumeration_n(EE);
  this->TEEnumerationToEnergy[EE_n] = energy;
  
  // also, populate TEEnergyToEnumeration
  this->TEEnergyToEnumeration.insert(make_pair(energy, EE_n));

  //cout << " Energies have been evaluated for " << TEEnumerationToEnergy.size() << " rotamer configurations." << endl;

//  cout << "TEEnumerationToEnergy.size(): " << TEEnumerationToEnergy.size() << endl;
//  cout << "TEEnergyToEnumeration.size(): " << TEEnergyToEnumeration.size() << endl;

  assert(TEEnumerationToEnergy.size() >= TEEnergyToEnumeration.size() );
}

double RotlibCollection::getExcitationEnergy(ExcitationEnumeration EE) {

  ExcitationEnumeration_n EE_n = this->_ExcitationEnumerationToExcitationEnumeration_n(EE);

  map< ExcitationEnumeration_n, double>::const_iterator ee_itr = this->TEEnumerationToEnergy.find(EE_n);
  if (ee_itr != this->TEEnumerationToEnergy.end()) {
    return ee_itr->second;
  }
  else {
    return 99999;
  }

}


void RotlibCollection::printExcitationEnergyTable() const {
  /* prints TEEnergyToEnumeration */
  multimap<double, ExcitationEnumeration_n>::const_iterator itr = TEEnergyToEnumeration.begin();

  /* prints header */
  cout << "Energy     ";
  for (ExcitationEnumeration_n::const_iterator citr = (itr->second).begin(); 
       citr != (itr->second).end();
       citr++) {
    //    cout << citr->first << "    "; // mutInfo
    //cout << this->_mutInfo_inverse_mapping[citr->first] << "    "; // mutInfo
    cout << this->_mutInfo_inverse_mapping.find(citr->first)->second << "    ";
  }
  cout << endl;


  /* prints content */
  for (; itr != TEEnergyToEnumeration.end(); itr++) {
    cout << itr->first << "   "; // itr->first is energy
    for (ExcitationEnumeration_n::const_iterator citr = (itr->second).begin();
	 citr != (itr->second).end();
	 citr++) {
      cout << citr->second << "     "; // Rank (enumeration)
    }
    cout << endl;
  }

}

void RotlibCollection::printExcitationEnergyTable(const string) const {

}

void RotlibCollection::printEmptyLatticeTable() const {
  /* Take EnumerationList and prints it! */
  cout << "Printing Empty Lattice Table! " << endl;
  cout << "Empty Lattice Table size: " << EnumerationList.size() << endl;

  // print EnumerationList

  for (set<ExcitationEnumeration_n>::const_iterator ExcitationEnumerationItr = this->EnumerationList.begin();
       ExcitationEnumerationItr != this->EnumerationList.end();
       ExcitationEnumerationItr++) {

    for (ExcitationEnumeration_n::const_iterator MutInfoRankItr = (*ExcitationEnumerationItr).begin();
	 MutInfoRankItr != (*ExcitationEnumerationItr).end();
	 MutInfoRankItr++) {

      //cout << "MutInfo is " << string(MutInfoRankItr->first) << " , Rank is " << MutInfoRankItr->second << ", ";
      cout << "MutInfo is " << this->_mutInfo_inverse_mapping.find(MutInfoRankItr->first)->second << " , Rank is " << MutInfoRankItr->second << ", ";
    }
    cout << endl;
  }
}

void RotlibCollection::printEmptyLatticeLinearEnergyTable() const {
  /* takes ELEnumerationToEnergy and prints it ! */

  map< ExcitationEnumeration_n, double >::const_iterator LevelItr = ELEnumerationToEnergy.begin(); // LevelItr is of pair(ExcitationEnumeration, double) type.
  ExcitationEnumeration_n::const_iterator ExiItr = ( LevelItr->first).begin();

  // First print header line
  for (; ExiItr != (*LevelItr).first.end(); ExiItr++) {
    cout << this->_mutInfo_inverse_mapping.find( ExiItr->first )->second << "  ";
  }
  cout << endl;

  // Now print meat
  for (; LevelItr != ELEnumerationToEnergy.end(); LevelItr++) {
    for (ExiItr = (LevelItr->first).begin(); ExiItr != (LevelItr->first).end(); ExiItr++) {
      cout << ExiItr->second << "    ";	// rank
    }
    cout << LevelItr->second << endl; // energy
    
  }

}

int RotlibCollection::cmpMaxRotamerConfigurations(int upperbound) {
  if (upperbound > this->maxRotamerConfigurations) {
    return 1;
  } else {
    return 0;
  }
  

}

void RotlibCollection::_sortAllRotlibByEmptyLatticeEnergy() {
  /* Sort all rotamer libraries by empty lattice E */
  ///  map<string, NtrlAARotlib*>::iterator RotlibItr;
  //map<string, Rotlib*>::iterator RotlibItr;
  map<int, Rotlib*>::iterator RotlibItr;
  assert(RotamerLibraryMap.size() != 0);
  for (RotlibItr = RotamerLibraryMap.begin(); RotlibItr != RotamerLibraryMap.end(); RotlibItr++) {
    assert(RotlibItr->second != NULL);
    RotlibItr->second->sort_by_empty_lattice_E();
    RotlibItr->second->reset_pstn();
  }

}

void RotlibCollection::_calcEmptyLatticeLinearSumEnergy() {
  /* calculates linear sum of indiivudal residue empty lattice energy and put them in ELEnumerationToEnergy */
  /* not particularly efficient, but works */

  assert(EnumerationList.size() != 0);
  ELEnumerationToEnergy.clear();
  int sizeOfSystem = RotamerLibraryMap.size();

  map<int, Rotlib*>::iterator RotlibItr;
  assert(RotamerLibraryMap.size() != 0);
  for (RotlibItr = RotamerLibraryMap.begin(); RotlibItr != RotamerLibraryMap.end(); RotlibItr++) {
    assert(RotlibItr->second != NULL);
    RotlibItr->second->reset_pstn();
  }  

  for (set< ExcitationEnumeration_n >::iterator LevelItr = EnumerationList.begin();
       LevelItr != EnumerationList.end();
       LevelItr++) {
    assert( (*LevelItr).size() == sizeOfSystem);

    
    double linearEnergyOfIndividualResidues = 0.0;
    /* Inner loop: loops over residues */
    ////for (ExcitationEnumeration::iterator ExiItr = (*LevelItr).begin();
    for (ExcitationEnumeration_n::const_iterator ExiItr = (*LevelItr).begin();
	 ExiItr != (*LevelItr).end();
	 ExiItr++) {

      unsigned short mutInfo_n = ExiItr->first;
      unsigned short Rank = ExiItr->second; // Rank starts from 0.
      //      NtrlAARotlib* thisRotlib = getNtrlAARotlib(mutInfo);
      Rotlib* thisRotlib = this->getNtrlAARotlib(mutInfo_n);
      
      //      AARotamer* aaRot = thisRotlib->get_empty_lattice_E_rot(Rank); // rotamer n starts from 1.  so add 1.
      Rotamer* aaRot = thisRotlib->get_empty_lattice_E_rot(Rank); // rotamer n starts from 1.  so add 1.
      // Below warning: should check to make sure aaRot isn't NULL!  or will seg fault.
      assert(aaRot != NULL);
      linearEnergyOfIndividualResidues += aaRot->get_empty_lattice_E();
      
    }
    ELEnumerationToEnergy.insert(make_pair(*LevelItr, linearEnergyOfIndividualResidues));
    
  }

  assert(ELEnumerationToEnergy.size() == EnumerationList.size());

}

void RotlibCollection::_calcLatestRoundEmptyLatticeLinearSumEnergy() {
  //  Debug debugInfo("RotlibCollection::_calcLatestRoundEmptyLatticeLinearSumEnergy()");
  /* calculates linear sum of indiivudal residue empty lattice energy and put them in ELEnumerationToEnergy */
  /* not particularly efficient, but works */
#ifdef DEBUG 
  assert(EnumerationList.size() != 0);
#endif

  this->latest_round_ELEnumerationToEnergy.clear();

  int sizeOfSystem = this->RotamerLibraryMap.size();

  //  cout << "ELEnumerationToEnergy size: " << ELEnumerationToEnergy.size() << endl;
  //  cout << "EnumerationList size: " << EnumerationList.size() << endl;

  map<int, Rotlib*>::iterator RotlibItr;
  assert(RotamerLibraryMap.size() != 0);
  for (RotlibItr = this->RotamerLibraryMap.begin(); RotlibItr != this->RotamerLibraryMap.end(); RotlibItr++) {
    assert(RotlibItr->second != NULL);
    RotlibItr->second->reset_pstn();
  }  

  for (list< ExcitationEnumeration_n >::iterator LevelItr = this->latest_round_EnumerationList.begin();
       LevelItr != this->latest_round_EnumerationList.end();
       LevelItr++) {
#ifdef DEBUG 
    assert( (*LevelItr).size() == sizeOfSystem);
#endif
    
    double linearEnergyOfIndividualResidues = 0.0;
    /* Inner loop: loops over residues */
    for (ExcitationEnumeration_n::const_iterator ExiItr = (*LevelItr).begin();
	 ExiItr != (*LevelItr).end();
	 ExiItr++) {
      unsigned short mutInfo_n = ExiItr->first;
      unsigned short Rank = ExiItr->second; // Rank starts from 0.

      Rotlib* thisRotlib = this->getNtrlAARotlib(mutInfo_n);
      
      Rotamer* aaRot = thisRotlib->get_empty_lattice_E_rot_after_sorted_by_empty_lattice_E(Rank); // miniscully faster.  so imperceptible i can't even time it.
      // Below warning: should check to make sure aaRot isn't NULL!  or will seg fault.
#ifdef DEBUG 
      assert(aaRot != NULL);
#endif
      linearEnergyOfIndividualResidues += aaRot->get_empty_lattice_E();
      
    }

    bool insertionSuccess;

    std::map<ExcitationEnumeration_n, double>::value_type tmpEEn_Energy_pair(*LevelItr, linearEnergyOfIndividualResidues);
    insertionSuccess = (ELEnumerationToEnergy.insert(tmpEEn_Energy_pair)).second;

    if (insertionSuccess) {
      this->EnumerationList.insert(*LevelItr);
      //this->latest_round_ELEnumerationToEnergy.insert(make_pair(*LevelItr, linearEnergyOfIndividualResidues));    
      this->latest_round_ELEnumerationToEnergy.insert(tmpEEn_Energy_pair);    
    }


  }
  //  cout << " Total number of sums done: " << i << endl;
  
  //debugInfo.out(" this->latest_round_ELEnumerationToEnergy.size() is: " + string(itoa(this->latest_round_ELEnumerationToEnergy.size()) ) );
  //debugInfo.out(" ELEnumerationToEnergy.size() after all's added is: " + string(itoa(ELEnumerationToEnergy.size()) ) );
  //debugInfo.out(" EnumerationList.size() is: " + string(itoa(EnumerationList.size()) ) );
  //#ifdef DEBUG 

  //  cout << " After : ELEnumerationToEnergy size: " << ELEnumerationToEnergy.size() << endl;
  //  cout << " After : EnumerationList size: " << EnumerationList.size() << endl;

  assert(ELEnumerationToEnergy.size() == EnumerationList.size());
  //#endif
}



 
void RotlibCollection::_copyELEnumerationToEnergyToEnergyMultimap() {
  /* Reverse key and value for ELEnumerationToEnergy and stores them in ELEnergyToEnumeration.  */

  assert(this->ELEnumerationToEnergy.size() != 0);
  this->ELEnergyToEnumeration.clear();

  //cout << " _copyELEnumerationToEnergyToEnergyMultimap()" << endl;
  
  map<ExcitationEnumeration_n, double>::const_iterator EnumItr = ELEnumerationToEnergy.begin();
  //  map<double, ExcitationEnumeration_n>::iterator EnergyItr;

  //cout << "entering for loop in _copyELEnumerationToEnergyToEnergyMultimap! " << endl;

  //  int i = 0;
  for (; EnumItr != ELEnumerationToEnergy.end(); EnumItr++) {
    //    i++;
    //this->ELEnergyToEnumeration.insert(make_pair(EnumItr->second, EnumItr->first)); // inverses (key, value).
    this->ELEnergyToEnumeration.insert(std::multimap< double, ExcitationEnumeration_n >::value_type(EnumItr->second, EnumItr->first) );
  }

  //  cout << "Total number of ELEnergyToEnumeration insertions done: " << i << endl;
  
  assert(ELEnumerationToEnergy.size() == ELEnergyToEnumeration.size() );

}



void RotlibCollection::_copyEnergyMultimapToELEnumerationToEnergy() {

  this->ELEnumerationToEnergy.clear();
  assert(this->ELEnergyToEnumeration.size() != 0);
  
  map<double, ExcitationEnumeration_n>::const_iterator enerItr = this->ELEnergyToEnumeration.begin();
  for (; enerItr != ELEnergyToEnumeration.end(); enerItr++)
    this->ELEnumerationToEnergy.insert(std::multimap< ExcitationEnumeration_n, double >::value_type(enerItr->second, enerItr->first) );
  
  assert(ELEnumerationToEnergy.size() == ELEnergyToEnumeration.size() );

}

void RotlibCollection::_copyEnergyMultimapToELEnumerationList() {
  
  this->EnumerationList.clear();
  assert(this->ELEnergyToEnumeration.size() != 0);

  map<double, ExcitationEnumeration_n>::const_iterator enerItr = this->ELEnergyToEnumeration.begin();
  for (; enerItr != ELEnergyToEnumeration.end(); enerItr++)
    this->EnumerationList.insert(enerItr->second);
  
  assert(ELEnumerationToEnergy.size() == ELEnergyToEnumeration.size() );


}


void RotlibCollection::_copyLatestRoundELEnumerationToEnergyToEnergyMultimap() {
  /* Reverse key and value for ELEnumerationToEnergy and stores them in ELEnergyToEnumeration.  */

  assert(this->ELEnumerationToEnergy.size() != 0);

  //  assert(this->latest_round_ELEnumerationToEnergy.size() != 0);
  //this->ELEnergyToEnumeration.clear();
  

  // cout << " _copyELEnumerationToEnergyToEnergyMultimap()" << endl;
  
  //map<ExcitationEnumeration, double>::const_iterator EnumItr = ELEnumerationToEnergy.begin();
  map<ExcitationEnumeration_n, double>::const_iterator EnumItr = this->latest_round_ELEnumerationToEnergy.begin();
  //  map<double, ExcitationEnumeration_n>::iterator EnergyItr;

  // cout << "entering for loop in _copyELEnumerationToEnergyToEnergyMultimap! " << endl;

  //  int i = 0;

  //  for (; EnumItr != ELEnumerationToEnergy.end(); EnumItr++) {
  for (; EnumItr != this->latest_round_ELEnumerationToEnergy.end(); EnumItr++) {
    //    i++;
    //ELEnergyToEnumeration.insert(make_pair(EnumItr->second, EnumItr->first)); // inverses (key, value). it's a multimap.
    this->ELEnergyToEnumeration.insert(std::multimap< double, ExcitationEnumeration_n >::value_type(EnumItr->second, EnumItr->first) );
  }

  //  cout << "Total number of ELEnergyToEnumeration insertions done: " << i << endl;
  //#ifdef DEBUG 
  assert(ELEnumerationToEnergy.size() == ELEnergyToEnumeration.size() );
  //#endif
}



long double RotlibCollection::_calcMaxNoClashRotamerConfigurations() {
  /* Calculates total number of passed clash-test passed rotamer configurations */
  //int totalN = 1;
  long double totalN = 1;
  ////  map<const std::string, Rotlib*>::const_iterator lib_itr = this->RotamerLibraryMap.begin();
  map<int, Rotlib*>::const_iterator lib_itr = this->RotamerLibraryMap.begin();
  for (; lib_itr != this->RotamerLibraryMap.end(); lib_itr++) {
    int libRotCount = lib_itr->second->n_rotamers_below_empty_lattice_energy(this->HIGHEST_ALLOWED_ROTAMER_E);
    totalN *= libRotCount;
          
  }
  return totalN;

}

////ExcitedRotamers RotlibCollection::_getNextRotamerDynamicMemory() {
ExcitedRotamers_n RotlibCollection::_getNextRotamerDynamicMemory() {

  /* (08/19/04) similar to getNextEmptyLatticeExcitationRotamers(), written with a different concept/purpose in mind.  Returns an empty ExcitationEnumeration map if this isn't possible (like when all steric-clash rotamers have been exhausted.) */
  
  assert(ELEnergyToEnumeration.size() != 0);
  assert(ELEnergyToEnumeration.size() == ELEnumerationToEnergy.size() );

  ////ExcitedRotamers toBeReturned;
  ExcitedRotamers_n toBeReturned;
  toBeReturned.clear();

  // first check if all clash-test-passed rotamer have already been exhausted. 
  if (this->currentEmptyLatticeExcitationN + 1 > this->maxRotamerConfigurations) {  // currentEmptyLatticeExcitationN starts from 0.
    return toBeReturned;
  }
  
  if (this->currentEmptyLatticeExcitationN >= ELEnergyToEnumeration.size()) {
    cout << "Reached end of ELEnergyToEnumeration! in _getNextRotamerDynamicMemory()." << endl;
    assert (toBeReturned.size() == 0);
    return toBeReturned;
  } 
  else {
    ////    multimap<double, ExcitationEnumeration >::iterator config_itr = this->ELEnergyToEnumeration.begin();
    multimap<double, ExcitationEnumeration_n >::iterator config_itr = this->ELEnergyToEnumeration.begin();
    if (this->currentEmptyLatticeExcitationN == 0) {
      
      //toBeReturned = this->getELExcitedRotamerFromEnumeration(config_itr->second);
      toBeReturned = this->getELExcitedRotamer_nFromEnumeration_n(config_itr->second);

    } 
    else {
      for (int count = 1; count <= this->currentEmptyLatticeExcitationN; count++) {  // currentEmptyLatticeExcitationN starts from 0.
	config_itr++;
      }
    }
    toBeReturned = this->getELExcitedRotamer_nFromEnumeration_n(config_itr->second);
  }


  return toBeReturned;

  // increment done in parent function: getNextDynamicMemoryRotamers_And_Expand().

}

void RotlibCollection::_expandDynamicMemory(ExcitationEnumeration_n& thisEE) {
  //  Debug debugInfo("RotlibCollection::_expandDynamicMemory(ExcitationEnumeration_n& thisEE)");

  assert(thisEE.size() != 0);

  // 1. makeExpansionVector, (if 26 residues, get 26 additional ExcitationEnumeration's.).  then, append this ExpansionVector to this->EnumerationList, which is a vector of ExcitationEnumerations.

  // make size(thisEE) copies of thisEE.  

  list< ExcitationEnumeration_n > expansionList;
  expansionList.clear();

  for (ExcitationEnumeration_n::iterator itr = thisEE.begin(); 
       itr != thisEE.end(); itr++) {

    ExcitationEnumeration_n a_copy_of_thisEE = thisEE;

    unsigned short mutInfo_n = itr->first;
    unsigned short currentN = a_copy_of_thisEE[mutInfo_n];
    // below: if currentN is greater than the total N of rotamer in mutInfo library, skip.  
    if (currentN + 1 + 1 > RotamerLibraryMap[mutInfo_n]->size()) {  // size(): starts from 1.  currentN: starts from 0.  hence need +1 +1.
      continue;
    }
    
    a_copy_of_thisEE[mutInfo_n]++;

    expansionList.push_back(a_copy_of_thisEE);
  }

  //  float t2 = (float) clock();

  //assert(thisEE.size() == expansionSet.size()); this is no longer true, because ones that exceed boundary will be truncated.  hence assertion no longer true.

  // Insertion: get rid of repeats, though not necessary because this vector will be processed into a map (not multimap) structure:
  //   i.e. this->ELEnumerationToEnergy is what's populated by EnumerationList through the function this->_calcEmptyLatticeLinearSumEnergy.  

  // Below: repeats not gotten rid of.
  //this->EnumerationList.insert(EnumerationList.end(), expansionV.begin(), expansionV.end()); // Insert expanded ones into set.
  //vector< ExcitationEnumeration > tempEnumerationList;
  
  this->latest_round_EnumerationList.clear();
  this->latest_round_ELEnumerationToEnergy.clear();

  //  debugInfo.out(" expansionSet size: " + string( itoa(expansionList.size()) ) );
  
  this->latest_round_EnumerationList.swap(expansionList);

  // 2. Convert this->EnumerationList to the two operatable multi-maps by the functions:
  
  //this->_calcEmptyLatticeLinearSumEnergy();
  this->_calcLatestRoundEmptyLatticeLinearSumEnergy();  // ALERT!  This subroutine now also takes care of proper element insertion of EnumerationList, ELEnumerationToEnergy, and latest_round_ELEnumerationToEnergy.

  //this->_copyELEnumerationToEnergyToEnergyMultimap(); 
  this->_copyLatestRoundELEnumerationToEnergyToEnergyMultimap(); // this->latest_round_ELEnumerationToEnergy updated from previous subroutine.


  // 3. that's it!  increment done in parent function: getNextDynamicMemoryRotamers_And_Expand()
  
}

void RotlibCollection::_buildOneSetUnderThresholdConfigs(multimap< double, ExcitationEnumeration_n > & previousConfigs, 
							 multimap< double, ExcitationEnumeration_n > & newConfigs, 
							 double availableE, int rotlib_n, int excitation_n, double crntRotE) {

  /* Subroutine for building just one set of configurations under:
     Rotlib index pstn, with the excitation number. */

  assert( availableE >= 0);

  typedef multimap< double, ExcitationEnumeration_n >::value_type eeType ;
  multimap< double, ExcitationEnumeration_n >::iterator upperBoundItr = previousConfigs.upper_bound(availableE); // first value past the last entry matching the argument
  
  for ( multimap< double, ExcitationEnumeration_n >::iterator itr = previousConfigs.begin(); itr != upperBoundItr; itr++) {

      double prevE = itr->first;
      double thisE = crntRotE + prevE;
    
      for (ExcitationEnumeration_n::iterator I = itr->second.begin(); I != itr->second.end(); I++) {
	//cout << I->first << endl;
	//cout << I->second << endl;
      }

      ExcitationEnumeration_n thisEE(itr->second);
      thisEE[rotlib_n] = excitation_n;
      
      newConfigs.insert(eeType(thisE, thisEE));
      //cout << "size of newConfigs: " << newConfigs.size() << endl;
    }
    //  }

  //    cout << " Exiting _buildOneSetUnderThresholdConfigs" << endl;

}


Rotlib* RotlibCollection::getNtrlAARotlib(const std::string mutInfo) {
  unsigned short mutInfo_n = this->_mutInfo_mapping.find(mutInfo)->second;
  return RotamerLibraryMap.find(mutInfo_n)->second;
}

Rotlib* RotlibCollection::getNtrlAARotlib(int mutInfo_n) {
  return RotamerLibraryMap.find(mutInfo_n)->second;

}

double RotlibCollection::getEmptyLatticeLinearSumEnergy( ExcitationEnumeration EE ) {
  ExcitationEnumeration_n EE_n = this->_ExcitationEnumerationToExcitationEnumeration_n(EE);
  return ELEnumerationToEnergy.find(EE_n)->second;
}


bool RotlibCollection::_handCheckEEIdentity(const ExcitationEnumeration_n& EE1, const ExcitationEnumeration_n& EE2) {

  // This helper function explicitly checks whether two ExcitationEnumeration objects are identical in value, and that includes how the values are internally ORDERS.  thankfully EE's keys are the mutInfo, so if mutInfo are wrongly ordered we're in deep sh*t anyway.
  // Returns true if the two are indeed identical, otherwise, if different, returns false.

  if (EE1.size() != EE2.size()) {
    return false;
  }

  //vector<string> EE1_keys, EE2_keys;
  vector<int> EE1_keys, EE2_keys;
  vector<int> EE1_values, EE2_values;

  ExcitationEnumeration_n::const_iterator itr1 = EE1.begin();
  ExcitationEnumeration_n::const_iterator itr2 = EE2.begin();

  for (; itr1 != EE1.end() and itr2 != EE2.end(); itr1++, itr2++) {
    if ((itr1->first != itr2->first) or (itr1->second != itr2->second)) {
      return false;
    }
  }

  return true;

}
