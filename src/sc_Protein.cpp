#include "defs.hpp"
#include "sc_bgf_handler.hpp"
//#include "GenericRotamer_defs.hpp"

#include "scream_vector.hpp"
#include "scream_matrix.hpp"

#include "sc_Protein.hpp"
#include "sc_ProteinComponent.hpp"
#include "scream_atom.hpp"
#include "sc_AAChain.hpp"
#include "sc_Ligand.hpp"
#include "sc_Water.hpp"
#include "sc_Hetatm.hpp"
#include "scream_tools.hpp"

#include "sc_SideChain.hpp"

#include "Rotamer.hpp"
#include "AARotamer.hpp"
#include "Rotlib.hpp"

#include "scream_rtf.hpp"

#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <iostream> 
#include <string>
#include <cassert>
//#include <typeinfo>
#include <algorithm>
#include <stdexcept>
#include <stdlib.h>

using namespace std;

//vcvicek
//Protein::Protein()
//: atom_v(),
//  new_mapping(),
//  rtf(NULL) 
//{
//  //ScreamAtomV temp;
//  //this->atom_v = temp;
//}

Protein::Protein(ScreamAtomV* atom_list) : atom_v(*atom_list), rtf(NULL){
  //cout << "Address in Protein(ScreamAtomV* atom_list) : &atom_v is: " << &atom_v << endl;
  this->InitDataStructures(*atom_list);
}

Protein::~Protein() {
  Debug debugInfo("in Protein::~Protein()");
  
  if (!this->rtf) {
    delete rtf;
  }

  vector<ProteinComponent*>::iterator itr;  // Components are always new'ed structures.  Hence, always deleted.
  for (itr = this->pc_v.begin(); itr != this->pc_v.end(); ++itr) {
    debugInfo.out("deleting ProteinComponent", "stdout");
    delete *itr;
  }
  
  //ScreamAtomVItr itr_v;     // If this list has non-zero size, then the structure is initialized by reading a bgf file.  Hence, SCREAM_ATOM's need to be deleted.

  /* deletion of atoms is no longer necessary as atom_v is merely storing a bunch of pointers.  */
  /* those pointers were NOT allocated in this Protein class, but allocated in a higher level class.  destruction of those pointers should, therefore, be take care of in that highlevel class. */
  /* specifically, in bgf_handler. */
  /* 7-23-06 */
  /* addHydrogens() function changes this.  However, atoms still need not be deleted at this level. */
  /* And that's because the this->atom_v shares the same address as whatever higher level construct that passes in the list of atoms. That higher level construct will take care of all atom deletion of that atom_v.  Whats created here, is deleted there; but whatever created here is immediately added to atom_v, and therefore safe. */

  debugInfo.out("Protein Destroyed.", "stdout");

}

AAChain* Protein::get_AAChain(string chain_desig) const {
  /* Need to assert that chain_design actually exists!!!! */
  Debug debugInfo("in Protein::get_AAChain(string chain_desig)");

  vector<ProteinComponent*>::const_iterator itr;
  bool CHAIN_FOUND = false;

  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    AAChain* temp = dynamic_cast<AAChain*>(*itr);     // temp = NULL if cannot be casted.
    if (temp != NULL) {

      if ( temp->get_chain_desig() == chain_desig) {
	CHAIN_FOUND = true;
	return dynamic_cast<AAChain*>(*itr);
      }
    } else {
      continue;
    }
  }
  if (CHAIN_FOUND == false) {
    string info = "CHAIN " + chain_desig +  " not found in protein.";
    debugInfo.out(info, "stdout");
    return NULL;
  }

}

Ligand* Protein::get_Ligand() const {
  /* Returns the ONLY ligand in system.  If there is no Ligand, returns a NULL */
  Debug dI("Protein::get_Ligand()");

  vector<ProteinComponent*>::const_iterator itr;
  bool LIGAND_FOUND = false;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    Ligand* temp = dynamic_cast<Ligand*>(*itr);     // temp = NULL if cannot be casted.
    if (temp != NULL) {
      LIGAND_FOUND = true;
      return dynamic_cast<Ligand*>(*itr);
    }
  }

  if (LIGAND_FOUND == false) {
    dI.out("No Ligands found in protein.", "stdout");
    return NULL;
  }
  

}

ProteinComponent* Protein::operator[](int n) const {
  Debug dI("Protein::operator[](int n)");

  try {
    return this->pc_v[n];
  } catch (out_of_range& e) {    // out of range error.  (out_of_range)
    dI.out(e.what(), "stdout");
    return NULL;
  }

}

ProteinComponent* Protein::get_Component_with_ChainName(string chn) const {
  Debug dI("Protein::get_Component_with_ChainName(string chn)");

  ProteinComponent* pc = NULL;
  vector<ProteinComponent*>::const_iterator itr = this->pc_v.begin();
  for (; itr != this->pc_v.end(); ++itr) {

    ScreamAtomV atom_list = (*itr)->getAtomList();
    if (atom_list[0]->getChain() == chn) return *itr;

  }

  if (pc == NULL) {
    string s = "No protein component with chain name " + chn + " found!  Quitting.";
    dI.out(s, "stdout");
    exit(2);
  }

}

Protein& Protein::operator=(const Protein& ptn) {

  if (this == &ptn) {
    return *this;
  } else {
    this->pc_v = ptn.pc_v;
    this->ptr_to_atom_v = ptn.ptr_to_atom_v;
    this->atom_v = *(this->ptr_to_atom_v);
  }
  return *this;
}


string Protein::get_res_type(string chn, int pstn) const {

  ScreamAtomV sc_atoms = this->get_sc_atoms(chn, pstn);
  string res_type = sc_atoms[0]->resName;
  return res_type;

}

ScreamAtomV Protein::get_sc_atoms(string chn, int pstn) const {

  return this->get_AAChain(chn)->get_sc_atoms(pstn);

}

ScreamAtomV Protein::get_sc_atoms(MutInfo mI) const {
  ScreamAtomV sc_atoms; sc_atoms.clear();

  if (mI.isClusterMutInfo()) {
    vector<MutInfo*> mI_list = mI.getAllMutInfos();
    
    for (vector<MutInfo*>::const_iterator mI_itr = mI_list.begin(); mI_itr != mI_list.end(); ++mI_itr) {

      if ((*mI_itr)->getRotConnInfo() == NULL) {
	string chn = (*mI_itr)->getChn();
	int pstn = (*mI_itr)->getPstn();
	ScreamAtomV mI_sc_atoms = this->get_sc_atoms(chn, pstn);
	sc_atoms.insert(sc_atoms.end(), mI_sc_atoms.begin(), mI_sc_atoms.end());
      }
      else {
	ScreamAtomV mI_sc_atoms = this->get_variable_atoms( (*mI_itr)->getRotConnInfo() );
	sc_atoms.insert(sc_atoms.end(), mI_sc_atoms.begin(), mI_sc_atoms.end());
      }
    }
  }
  else {
    if (mI.getRotConnInfo() == NULL) {
      sc_atoms = this->get_AAChain(mI.getChn())->get_sc_atoms(mI.getPstn());
    }
    else {
      sc_atoms = this->get_variable_atoms( mI.getRotConnInfo() );
    }
  }
  return sc_atoms;
}

ScreamAtomV Protein::get_variable_atoms(RotConnInfo* rotConnInfo) const { 

  map<int, int> atomMapping = rotConnInfo->atom_n_map;
  vector<int> sideChainAtoms =  rotConnInfo->side_chain_atoms;
  vector<int> variableAtomsNumInPtn = this->_getAtomNumbersInProtein(atomMapping, sideChainAtoms);
  
  ScreamAtomV variableAtomsInPtn = this->getTheseAtoms(variableAtomsNumInPtn);
  return variableAtomsInPtn;

}

ScreamAtomV Protein::get_visible_in_EL_mutInfo_atoms(MutInfo& mI, RotConnInfo* rotConnInfo) const {
  ScreamAtomV aL; aL.clear();

  if (mI.getChn() == "Z") {
    // what should i do here?  not do anything for now.
    //if (rotConnInfo != NULL)
    //aL = this->get_variable_atoms(rotConnInfo);
  }
  else {
    string chn = mI.getChn();
    ScreamAtomV chn_atom_v = this->get_AAChain(chn)->getAtomList();

    int pstn = mI.getPstn();
    ScreamAtomV mI_visible_atoms;

    for (ScreamAtomVConstItr itr=chn_atom_v.begin(); itr!=chn_atom_v.end(); ++itr)
      if ( ((*itr)->resNum == pstn) and ((*itr)->flags & 0x8) )  // i.e. correct pstn AND it's an already fixed atom.
	mI_visible_atoms.push_back(*itr);
    
  }
  return aL;

}

void Protein::addHydrogens() {
  Debug debugInfo("Protein::addHydrogens()");

  // 7-23-06: adding hydrogens by recognizing the C_33 type and what nots read in from forcefield.  Run this after connectivity has been added.
  // The variable in X_YZ (e.g. C_33).  
  // Y matters the most.  Than X, depends on C, S, O or N.  O, N, S are H-bond acceptors and donors.
  // This routine is simple; it merely makes sure the correct topologies.  H-Bonds are not necessarily optimized in this routine.

  // Atom creation is done in this function.  No deletion done, destruction of atoms done in high level construct; need to make sure whatever atom that is created is immediately added to this->atom_v else won't get deleted.
  ScreamAtomV newHydrogens_global; newHydrogens_global.clear();
  ScreamAtomV cysHydrogens_global; cysHydrogens_global.clear();

  for (vector<ProteinComponent*>::iterator itr = this->pc_v.begin(); itr != this->pc_v.end(); ++itr) {
    AAChain* chain = dynamic_cast<AAChain*>(*itr);      // temp = NULL if cannot be casted.
    if (chain != NULL) {
      debugInfo.out("Adding Hydrogens to a peptide chain..."); // could also be modifying forcefield type

      ScreamAtomV aaAtoms;

      vector<int> res_numbers = chain->getResidueNumbers();
      for (vector<int>::const_iterator itr = res_numbers.begin();
	   itr != res_numbers.end(); ++itr) {
	AminoAcid* aa = (*chain)[*itr];
	string resName = aa->residue_type();
	debugInfo.out(resName);

	ScreamAtomV aaAtoms = aa->getAtomList();

	for (ScreamAtomVItr itr2 = aaAtoms.begin(); itr2 != aaAtoms.end(); ++itr2) {
	  debugInfo.out( (*itr2)->return_bgf_line() );
	  /* HIS exceptions. */	  
	  if (resName == "HIS" or resName == "HSE" or resName == "HDD") {
	    if ( scream_tools::strip_whitespace((*itr2)->getAtomLabel()) == "ND1" ) {
	      debugInfo.out( (*itr2)->getAtomLabel() + " hydrogen excluded by default. ");
	      continue;
	    }
	  }
	  
	  else if (resName == "HSE" or resName == "HDD") {
	    if ( scream_tools::strip_whitespace( (*itr2)->getAtomLabel()) == "NE2") {
	      debugInfo.out( (*itr2)->getAtomLabel() + " hydrogen excluded by default. ");
	      continue;
	    }
	  }
	    
	  /* Calculate new hydrogen coordinates. */
	  vector<ScreamVector> newHydrogenCoords = scream_tools::generateHydrogenCoords( *itr2 );
	  debugInfo.out(" Hydrogens added: " + string(itoa(newHydrogenCoords.size()) ) );
	  
	  /* Create H SCREAM_ATOM's from hydrogen coordinates. */
	  ScreamAtomV newHydrogens = scream_tools::createHydrogens(newHydrogenCoords, *itr2); // includes Atom label, FF type (H_ or H___A), no charges.
	  debugInfo.out(" Hydrogens created: " + string(itoa(newHydrogens.size()) ) );
	  
	  /* Insert them into new hydrogen list for later insertion. CYS exception. */
	  if (resName == "CYS") {
	    cysHydrogens_global.insert(cysHydrogens_global.end(), newHydrogens.begin(), newHydrogens.end());
	    debugInfo.out(" CYS hydrogens inserted.");
	  }
	  else {
	    newHydrogens_global.insert(newHydrogens_global.end(), newHydrogens.begin(), newHydrogens.end());
	    debugInfo.out(" Hydrogens inserted. " );
	  }
	  
	} 
	
	
      }// end for (vector<int>::const_iterator itr = res_numbers.begin();	   itr != res_numbers.end(); ++itr)
    } else { // if chain == NULL, i.e. not AAchain.
      // do nothing for now.  If one day implement HETATM rtf stuff, just replace this block, do a dynamic cast<Hetatm*> or something.
    }
  }
  int beforeHAdd_atom_v_size = atom_v.size();
  int i = 1;
  /* Decide if CYS atoms should be added. For now, yes. */
  newHydrogens_global.insert(newHydrogens_global.end(), cysHydrogens_global.begin(), cysHydrogens_global.end());
  
  for (ScreamAtomVItr ii = newHydrogens_global.begin(); ii != newHydrogens_global.end(); ++ii, ++i) {
    (*ii)->n = beforeHAdd_atom_v_size + i;
  }
  this->atom_v.insert(atom_v.end(), newHydrogens_global.begin(), newHydrogens_global.end() );
  
  
  
  // need to sort atom_v after hydrogens have been added.
  this->_fix_entire_atom_list_ordering();
  
}


void Protein::addConnectivity() {
  Debug debugInfo("Protein::addConnectivity()");

  // vcvicek
  char * SCREAM_NEW_RTF = getenv("SCREAM_NEW_RTF");
  if (SCREAM_NEW_RTF == NULL) {printf("error: enviromental variable SCREAM_NEW_RTF is not set \n"); exit(1);}
  
  string rtf_file = string(SCREAM_NEW_RTF);
  //string rtf_file = "/ul/victor/SCREAM/SCREAM_v2.0/SCREAM/src/SCREAM.rtf";
  if (!this->rtf) {
    this->rtf = new SCREAM_RTF(rtf_file);
  }
  for (vector<ProteinComponent*>::iterator itr = this->pc_v.begin(); itr != this->pc_v.end(); ++itr) {
    AAChain* chain = dynamic_cast<AAChain*>(*itr);      // temp = NULL if cannot be casted.
    if (chain != NULL) {
      debugInfo.out("Adding bonds to a peptide chain...");

      SCREAM_ATOM* _prev_C = NULL;
      ScreamAtomV aaAtoms;
      map<string, SCREAM_ATOM*> aaAtom_map;

      vector<int> res_numbers = chain->getResidueNumbers();
      for (vector<int>::const_iterator itr = res_numbers.begin();
	   itr != res_numbers.end(); ++itr) {
	AminoAcid* aa = (*chain)[*itr];
	string resName = aa->residue_type();
	debugInfo.out(resName);
	// map of atom label to atoms

	ScreamAtomV aaAtoms = aa->getAtomList();
	aaAtom_map.clear();
	for (ScreamAtomVItr itr2 = aaAtoms.begin(); itr2 != aaAtoms.end(); ++itr2) {
	  //debugInfo.out( (*itr)->return_bgf_line() );
	  aaAtom_map[scream_tools::strip_whitespace((*itr2)->getAtomLabel())] = *itr2;
	}
	// prepping RTF stuff
	AminoAcid_RTF* aaRTF = this->rtf->get_AminoAcid_RTF(resName);
	multimap<string, string> aaBonds = aaRTF->return_bonds_table();
	debugInfo.out(string(itoa(aaBonds.size())));
	// adding bonds
	for (multimap<string, string>::iterator itr2 = aaBonds.begin(); itr2 != aaBonds.end(); ++itr2) {
	  string a1_label = itr2->first;
	  string a2_label = itr2->second;
	  debugInfo.out(a1_label);
	  debugInfo.out(a2_label);
	  SCREAM_ATOM* a1 = aaAtom_map[a1_label];
	  SCREAM_ATOM* a2 = aaAtom_map[a2_label]; // behavior of maps: if can't find key, default object is added to map and returned, which is NULL for SCREAM_ATOM*.
	  if (a1 and a2) { // both found; or not NULL.
	    a1->make_bond(a2);
	  }
	}

	// REMEMBER TO ADD +N--C and N-- -C bonds (the peptide bonds)  no boundary condition need be considered.
	if (_prev_C == NULL) { // first time through
	  _prev_C = aaAtom_map["C"];
	  if (_prev_C == NULL) {
	    cerr << " Mainchain C atom not found on a amino acid number " << *itr << "on chain " << chain->get_chain_desig() << endl;
	  }
	}
	else {
	  SCREAM_ATOM* _crnt_N = aaAtom_map["N"];
	  if (_crnt_N == NULL) {
	    cerr << " Mainchain N atom not found on a amino acid number " << *itr << "on chain " << chain->get_chain_desig() << endl;
	  }
	  _crnt_N->make_bond(_prev_C);
	  _prev_C = aaAtom_map["C"];
	}

      } // end for (vector<int>::const_iterator itr = res_numbers.begin();	   itr != res_numbers.end(); ++itr)
    } else { // if chain == NULL, i.e. not AAchain.
      // do nothing for now.  If one day implement HETATM rtf stuff, just replace this block, do a dynamic cast<Hetatm*> or something.
    }
  }

}

void Protein::assignFFType() {
  // Note: Should call this after assign connectivity.
  Debug debugInfo("Protein::assignFFType()");


  // 7-21-06: hardcoded path
  //vcvicek
  char * SCREAM_NEW_RTF = getenv("SCREAM_NEW_RTF");
  if (SCREAM_NEW_RTF == NULL) {printf("error: enviromental variable SCREAM_NEW_RTF is not set \n"); exit(1);}
  string rtf_file = string(SCREAM_NEW_RTF);
  //string rtf_file = "/ul/victor/SCREAM/SCREAM_v2.0/SCREAM/src/SCREAM.rtf";
  if (!this->rtf) {
    this->rtf = new SCREAM_RTF(rtf_file);
  }
  for (vector<ProteinComponent*>::iterator itr = this->pc_v.begin(); itr != this->pc_v.end(); ++itr) {
    AAChain* chain = dynamic_cast<AAChain*>(*itr);      // temp = NULL if cannot be casted.
    if (chain != NULL) {
      debugInfo.out("Adding FF type to a peptide chain..."); // could also be modifying forcefield type

      ScreamAtomV aaAtoms;

      vector<int> res_numbers = chain->getResidueNumbers();
      for (vector<int>::const_iterator itr = res_numbers.begin();
	   itr != res_numbers.end(); ++itr) {
	AminoAcid* aa = (*chain)[*itr];
	string resName = aa->residue_type();
	debugInfo.out(resName);

	// prepping RTF stuff
	AminoAcid_RTF* aaRTF = this->rtf->get_AminoAcid_RTF(resName);
	ScreamAtomV aaAtoms = aa->getAtomList();

	// adding ff type for each atom in amino acid
	for (ScreamAtomVItr itr2 = aaAtoms.begin(); itr2 != aaAtoms.end(); ++itr2) {
	  debugInfo.out( (*itr2)->return_bgf_line() );
	  string atom_label = scream_tools::strip_whitespace((*itr2)->getAtomLabel());
	  (*itr2)->setAtomType( aaRTF->get_ff_type(atom_label) );
	    
	}
	
	
      } // end for (vector<int>::const_iterator itr = res_numbers.begin();	   itr != res_numbers.end(); ++itr)
    } else { // if chain == NULL, i.e. not AAchain.
      // do nothing for now.  If one day implement HETATM rtf stuff, just replace this block, do a dynamic cast<Hetatm*> or something.
    }
  }
  
}

void Protein::assignFFType(SCREAM_RTF* rtf) {
  // nothing for now
}

vector<MutInfo> Protein::residuesAroundAtomN(vector<int> n_list, double distance, string inclusive_mode) const {
  // convert everything to atoms and the pass on to residuesAroundAtom
  ScreamAtomV atom_list = this->getTheseAtoms(n_list);
  return this->residuesAroundAtom(atom_list, distance, inclusive_mode);
  
}

vector<MutInfo> Protein::residuesAroundResidue(vector<MutInfo> mI_list, double distance, string inclusive_mode) const {
  // convert everything to atoms and then pass on to residuesAroundAtom.
  ScreamAtomV atom_list; atom_list.clear(); 
  for (vector<MutInfo>::const_iterator itr = mI_list.begin(); itr != mI_list.end(); ++itr) {
    ScreamAtomV mI_atoms = this->get_sc_atoms(*itr);
    atom_list.insert(atom_list.end(), mI_atoms.begin(), mI_atoms.end());
  }
  return this->residuesAroundAtom(atom_list, distance, inclusive_mode);
}

vector<MutInfo> Protein::residuesAroundChain(vector<string> chn_list, double distance, string inclusive_mode) const {
  // convert everything to atoms and then pass on to residuesAroundAtom.
  ScreamAtomV atom_list; atom_list.clear(); 
  for (vector<string>::const_iterator itr = chn_list.begin(); itr != chn_list.end(); ++itr ) {
    //ScreamAtomV chn_atoms = this->get_AAChain(*itr)->getAtomList();
    ScreamAtomV chn_atoms = this->get_Component_with_ChainName(*itr)->getAtomList();
    atom_list.insert(atom_list.end(), chn_atoms.begin(), chn_atoms.end());
  }
  return this->residuesAroundAtom(atom_list, distance, inclusive_mode);
}

vector<MutInfo> Protein::residuesAroundAtom(ScreamAtomV& ligand_atoms, double distance, string inclusive_mode) const {
  // named "ligand atoms" because most of the time want to get residues around a ligand.

  std::set<MutInfo> mI_set;
  double dist_sq = distance * distance;

  ScreamAtomV all_atoms_list = this->getAtomList();

  
  for (ScreamAtomVConstItr itr = all_atoms_list.begin(); itr != all_atoms_list.end(); ++itr) {
    if (inclusive_mode == "SideChainOnly") { // means only want sidechains.
      if ( ! scream_tools::is_SC_atom( (*itr)->atomLabel ) ) continue;
    }
    else if (inclusive_mode == "BackBoneOnly") {
      if ( ! scream_tools::is_BB_atom( (*itr)->atomLabel ) ) continue;
    }
    // else WholeResidue.
    for (ScreamAtomVConstItr itr2 = ligand_atoms.begin(); itr2 != ligand_atoms.end(); ++itr2) {
      if ( (*itr)->distance_squared(*itr2) > dist_sq ) continue;
     
      // excluding CYX.  9-10-05, temp measure.
      if ( (*itr)->getResName() == "CYX" ) {
	continue;
      }

      string chn = (*itr)->getChain();
      int pstn = (*itr)->getResNum();
      string AA = (*itr)->oneLetterResName;
      string mI_string = AA + itoa(pstn) + "_" + chn;
      MutInfo mI(mI_string);
      mI_set.insert(mI);

    }
  }

  vector<MutInfo> mI_list;
  for (std::set<MutInfo>::const_iterator itr = mI_set.begin(); itr != mI_set.end(); ++itr) {
    mI_list.push_back(*itr);
  }

  return mI_list;
}

ScreamAtomV& Protein::getAtomList() {
  
  return (this->atom_v);

}

const ScreamAtomV& Protein::getAtomList() const {

  return (this->atom_v);

}

ScreamAtomV Protein::getNewAtomList() {

  return this->atom_v;

}


SCREAM_ATOM* Protein::getAtom(int wantedAtomN) const {
  /* Returns SCREAM_ATOM with number n in Protein.*/
  SCREAM_ATOM* atomWanted;
  
  atomWanted = this->atom_v[wantedAtomN-1];
  
  /* Do the following if got atom not the one we wanted */
  if (atomWanted->n != wantedAtomN) {
    for (ScreamAtomVConstItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      int loopingAtomN = (*itr)->n;
      if (wantedAtomN == loopingAtomN) {
	atomWanted = (*itr);
	break;
      }
    }
  }

  return atomWanted;

}


SCREAM_ATOM* Protein::getAtom(MutInfo mI, string atom_label) const {
  /* Returns SCREAM_ATOM with MutInfo mI and atom_label in Protein.  If none found, return NULL..*/
  SCREAM_ATOM* atomWanted = NULL;
  atom_label = scream_tools::strip_whitespace(atom_label);
  string chn = mI.getChn();
  int pstn = mI.getPstn();
  for (int i=0; i<this->atom_v.size(); ++i) {
    atomWanted = this->atom_v[i];
    if ((atomWanted->getChain() == chn) and (atomWanted->getResNum() == pstn))
      if (atomWanted->stripped_atomLabel == atom_label)
	break;
  }
  return atomWanted;
}

ScreamAtomV Protein::getTheseAtoms(vector<int>& atomNumbers) const {
  /* Returns a vector of SCREAM_ATOM*'s from the list provided. */
  Debug dI("Protein::getTheseAtoms(vector<int>& atomNumbers)");

  ScreamAtomV wantedAtoms;
  for (vector<int>::const_iterator itr = atomNumbers.begin(); itr != atomNumbers.end(); ++itr) {
    SCREAM_ATOM* a = this->getAtom(*itr);
    if (a == NULL) {
      dI.out("Warning: Atom not found! in Protein::getTheseAtoms.", "stdout");
    }
    wantedAtoms.push_back(this->getAtom(*itr));
  }
  return wantedAtoms;
}

int Protein::totalComponents() const {

  return pc_v.size();

}

int Protein::mutationDone() const {
  // Done by checking size of new_mapping.
  if (this->new_mapping.size() == 0) {
    return 0;
  } else {
    return 1;
  }

}

void Protein::setMutInfoStrainEnergy(MutInfo mI, double strainE) {
  this->mutInfoStrainEnergies[mI] = strainE;
}

double Protein::getMutInfoStrainEnergy(MutInfo mI) {
  return this->mutInfoStrainEnergies[mI];
}

void Protein::printAtomFlagStates() {
  for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
    string l = (*itr)->return_bgf_line();
    cout << l.substr(0, 79) << " " << (*itr)->flags << endl;
  }
}

int Protein::ntrlRotamerPlacement(string chn, int pstn, AARotamer* rot) {
  /* rotamer placement, modifies underlying protein.*/
  /* rot: 100% of the time, this will be a AARotamer.*/
  /* AARotamer* doesn't need to be same as Aminoacid as current one--mutation facility */
  /* Return 0 if no mutations, returns 1 if mutations have occurred, may need to reset EE. */
  Debug dI("Protein::ntrlRotamerPlacement(string chn, int pstn, AARotamer* rot)");

#ifdef DEBUG
  assert(rot != NULL);
  // following test wrong; dynamic_cast can't cast a AARot* into an AARot*--AARotamer isn't virtual.
//   AARotamer* debugAArot = dynamic_cast<AARotamer*>(rot);
//   //  if (rot == NULL) {
//   if (debugAArot) {
//     cerr << "error in Protein::ntrlRotamerPlacement()" << endl;
//     cerr << "Rotamer passed in cannot be casted to an AARotamer" << endl;
//     exit(2);
//     return 0;
//   }
#endif
  
  AARotamer* AArot = rot; // code should be such that AARotamer* rot is guaranteed to be a AARotamer*.

  AAChain* scream_chn = this->get_AAChain(chn);

#ifdef DEBUG
  assert(scream_chn != NULL);
#endif

  AminoAcid* AA = scream_chn->get_aa_m()[pstn];
  if (AA == NULL) {
    cerr << "Pstn " << pstn << " does not exist on chain " << chn << " !" << endl;
    cerr << "Program exiting, from Protein::ntrlRotamerPlacement." << endl;
    exit(2);
  }
  
  string crntResName = this->get_res_type(chn, pstn);
  string newResName = AArot->get_resName();

  if (crntResName == newResName) { // no mutations, simple placement
    // no need to manage memory, but do need to reset new_mapping.
    this->new_mapping.clear();

  } else { // with mutations
    cout << crntResName << " ==> " << newResName << " Mutation! " << endl;

#ifdef DEBUG
    dI.out(" in sc_Protein PRINTING ORIGINAL_ROTAMER (deepcopied) sidechain atom connectivity info", "stdout");
    ScreamAtomV sc_atoms = AArot->get_sc_atoms();
    for (ScreamAtomVConstItr itr = sc_atoms.begin(); itr != sc_atoms.end(); ++itr) {
      (*itr)->dump();
      cout << "  Number of connected atoms: " << (*itr)->connectivity_m.size() << endl;
      cout << "  THEY ARE: " << endl;
      for (map<SCREAM_ATOM*, int>::const_iterator conn_itr = (*itr)->connectivity_m.begin();
	   conn_itr != (*itr)->connectivity_m.end(); ++conn_itr) {
	cout << "    address of connected atom: " << conn_itr->first << endl;
	(conn_itr->first)->dump(); 
      }
    }
#endif

    /* Part 1.  Memory management at the atoms level, allocation and fixing connectivities for atoms. */
    ScreamAtomV old_sc_atoms; // atoms to be replaced and deleted
    map<SCREAM_ATOM*, SCREAM_ATOM*> map_rot_to_new_sc_atoms; // map from rotamer atoms to would-be created new atoms, necessary to instantiate connectivities.

    ScreamAtomV new_sc_atoms = this->_mutationHelpers_alloc_new_sc_atoms(chn, pstn, AArot, old_sc_atoms, map_rot_to_new_sc_atoms);

    this->_mutationHelpers_connectivities_fix(chn, pstn, AArot, new_sc_atoms, map_rot_to_new_sc_atoms); // also fixes memory for the GLY and PRO special cases.

    this->_mutationHelpers_PRO_connectivities_and_numbering_fix(chn, pstn, AArot, new_sc_atoms, map_rot_to_new_sc_atoms);

    this->_mutationHelpers_GLY_connectivities_and_numbering_fix(chn, pstn, AArot, new_sc_atoms, map_rot_to_new_sc_atoms);

    /* Part 2.  Insertion of new atoms, deleting old sc atoms, renumbering, renaming resname on backbone.*/

    this->_mutationHelpers_init_new_mapping(this->atom_v);

    this->_mutationHelpers_insert_new_sc_atoms(chn, pstn, AArot, new_sc_atoms, old_sc_atoms);

    this->_mutationHelpers_delete_old_sc_atoms(old_sc_atoms);

    this->_mutationHelpers_renaming_mut_atoms(chn, pstn, newResName, new_sc_atoms);

    //this->_fix_entire_atom_list_ordering();
    this->_fix_residue_in_atom_list_ordering(chn, pstn);

    /* Part 3.  Memory management for high level constructs, such as AminoAcid. Delete old Aminoacid and make new AminoAcid.*/
    // these data structures need to be either reinstantiated or removed: AminoAcid, AARotamer, SideChain, Backbone.  easiest way to do this is probably just to delete AminoAcid, and re-instantiate AARotamer, SideChain and Backbone.  Backbone: from old backbone.  SideChain: from new_sc_atoms.  AARotamer: from Backbone and Sidechain.  AminoAcid: end terminals.  Keep track of end terminal guys.  To do this: in AAChain.
    ScreamAtomV new_atom_list;
    for (ScreamAtomVConstItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      if ((*itr)->chain == chn and (*itr)->resNum == pstn ) {
	new_atom_list.push_back(*itr);
      }
    }
    AminoAcid* new_AminoAcid = new AminoAcid(new_atom_list); 


    this->get_AAChain(chn)->replace_AminoAcid(pstn, new_AminoAcid); // this function deletes the old AminoAcid, thus achieving "conservation of memory"

    this->_fix_charges();  // need to fix charges as well

#ifdef DEBUG
    map<int, AminoAcid*> aa_m = this->get_AAChain(chn)->get_aa_m();
    for (map<int, AminoAcid*>::const_iterator itr = aa_m.begin();
	 itr != aa_m.end(); ++itr) {
      itr->second->print_Me();
    }
#endif

    /* Part 4.  Memory management for structures that created atom_lists in the first place, sc_bgf_handler, ModSimSCREAMInterface etc. */
    /* Further note: only need to do this when class Protein messes with atom_lists. */
    //*ptr_to_atom_v = this->atom_v;  

    dI.out(" Done mutation checkpoint.", "stdout");
  }
  scream_chn->SC_replacement(pstn, AArot, this->placementMethod, this->CreateCBParameters); // protein has been modified.

  if (crntResName == newResName)
    return 0; 
  else
    return 1; // returns 1 if there is a mutation
}

AARotamer* Protein::getAARotamer(std::string chn, int pstn) {
  
  AAChain* scream_chn = this->get_AAChain(chn);
  AARotamer* rotamerRequested = scream_chn->operator[](pstn)->get_rot();
  return rotamerRequested;
}


int Protein::conformerPlacement(Rotamer* conformer, RotConnInfo* rotConnInfo) {
  /* confomer placement. */
  assert(rotConnInfo != NULL and conformer != NULL);
  Rotamer tempConformer;
  tempConformer.deepcopy(*conformer);
  
  map<int, int> atomMapping = rotConnInfo->atom_n_map;
  vector<int> anchorAtoms = rotConnInfo->anchor_pts;
  vector<int> sideChainAtoms = rotConnInfo->side_chain_atoms;
  vector<int> exactMatchAtoms = rotConnInfo->atoms_of_exact_match;

  // cout << " Size of rotConnInfo.size() " << rotConnInfo->side_chain_atoms.size() << endl;

  /* 1) decide which are the anchor atoms and which ones are the sidechain/variable atoms */
  vector<int> anchorAtomsNumInProtein = this->_getAtomNumbersInProtein(atomMapping, anchorAtoms);
  vector<int> variableAtomsNumInProtein = this->_getAtomNumbersInProtein(atomMapping, sideChainAtoms);

  ScreamAtomV anchorAtomsInProtein = this->getTheseAtoms(anchorAtomsNumInProtein);
  //cout << "Size of anchorAtomsInProtein is: " << anchorAtomsInProtein.size() << " (Should be 3) " << endl;
  ScreamAtomV variableAtomsInProtein = this->getTheseAtoms(variableAtomsNumInProtein);
  //cout << "Size of variableAtomsInProtein is: " << variableAtomsInProtein.size() << " (Should be 6 for ASP) " << endl;

  /* 2) do match-anchor-atoms-matrix */
  anchorAtomsInProtein;
  ScreamAtomV anchorAtomsInConformer = tempConformer.getTheseAtoms(anchorAtoms);
  ScreamAtomV variableAtomsInConformer = tempConformer.getTheseAtoms(sideChainAtoms);

  //for_each (anchorAtomsInProtein.begin(), anchorAtomsInProtein.end(), dump);
  //for_each (anchorAtomsInConformer.begin(), anchorAtomsInConformer.end(), dump);

  // Match conformer atom to protein atom //
  bool exactMatch = 0;
  if (exactMatchAtoms.size() == 3) {
    exactMatch == 1;
  }
  // temp measure: for DEBUGGING PURPOSES
  exactMatch = 1;
  pair<ScreamMatrix, ScreamVector> T = scream_tools::getTransformationPairFromAtoms(anchorAtomsInProtein, anchorAtomsInConformer, exactMatch);

  //T.first.printMe();//first is Matrix
  //T.second.printMe();//second is the vector/translation
  
  /* 3) Transform variable sidechain atoms.  Note: should do translate and THEN transform. */
  //  cout << "Print VARIABLE ATOMS in CONFORMER" << endl;
  //  for_each (variableAtomsInConformer.begin(), variableAtomsInConformer.end(), dump);
  tempConformer.get_sc()->transform(T.first);
  //  cout << "Print VARIABLE ATOMS in CONFORMER after transform T.first (Identity Matrix)" << endl;
  //  for_each (variableAtomsInConformer.begin(), variableAtomsInConformer.end(), dump);
  tempConformer.get_sc()->translate(T.second);
  //cout << "Print VARIABLE ATOMS in CONFORMER after transform T.second (Translation Vector)" << endl;
  //  for_each (variableAtomsInConformer.begin(), variableAtomsInConformer.end(), dump);
  
  /* 4) Copy them over onto the protein. */
  assert(variableAtomsInConformer.size() == variableAtomsInProtein.size());
  ScreamAtomVConstItr conformer_itr = variableAtomsInConformer.begin();
  ScreamAtomVItr proteinAtom_itr = variableAtomsInProtein.begin();
  //cout << "This is COnformer Number: " << tempConformer.get_rotamer_n() << endl;
  for (; conformer_itr != variableAtomsInConformer.end(); ++conformer_itr, ++proteinAtom_itr) {
    //cout << "Conformer atom number: " << (*conformer_itr)->n << endl;
    (*proteinAtom_itr)->copyJustCoords( *(*conformer_itr) );   

  }
  // Done!!!
  //cout << "Should match: N position " << endl;
  //cout << "Print VARIABLE ATOMS in CONFORMER" << endl;
  //for_each (variableAtomsInConformer.begin(), variableAtomsInConformer.end(), dump);
  //for (ScreamAtomV::const_iterator itr = variableAtomsInConformer.begin(); itr != variableAtomsInConformer.end(); ++itr) {
  //(*itr)->dump();
  //  }
  //  cout << "Printe in VARIABLE ATOMS IN PRFOTEIN " << endl;
  //  for_each (variableAtomsInProtein.begin(),variableAtomsInProtein.end(), dump);
  
  return true;
}

Rotamer* Protein::conformerExtraction(RotConnInfo* rotConnInfo) {
  assert(rotConnInfo != NULL);

  /* creating an atom list with all the Rotamer atom's in it */

  map<int, int> atomMapping = rotConnInfo->atom_n_map;
  vector<int> anchorAtoms = rotConnInfo->anchor_pts;
  vector<int> sideChainAtoms = rotConnInfo->side_chain_atoms;

  vector<int> anchorAtomsNumInPtn = this->_getAtomNumbersInProtein(atomMapping, anchorAtoms);
  vector<int> variableAtomsNumInPtn = this->_getAtomNumbersInProtein(atomMapping, sideChainAtoms);
  ScreamAtomV anchorAtomsInPtn = this->getTheseAtoms(anchorAtomsNumInPtn);
  ScreamAtomV variableAtomsInPtn = this->getTheseAtoms(variableAtomsNumInPtn);

  ScreamAtomV allConformerAtomsInPtn(anchorAtomsInPtn);
  allConformerAtomsInPtn.insert(allConformerAtomsInPtn.end(), variableAtomsInPtn.begin(), variableAtomsInPtn.end());

  //cout << "dumping allConformerAtomsInPtn! " << endl;
  //  for_each(allConformerAtomsInPtn.begin(), allConformerAtomsInPtn.end(), dump);

  /* Actually create a conformer.  This conformer goes out of scope when this function exits.  */
  //cout << " !!!!!!!!Size of rotConnInfo->side_chain_atoms.size() " << rotConnInfo->side_chain_atoms.size() << endl;
  Rotamer tempConformer(allConformerAtomsInPtn, rotConnInfo, true);

  /* Now deep copy it.  This rotamer is returned; must be deleted outside!  Note: deepcopy sets allocatedScreamAtoms flag to true. */
  Rotamer* newConformer = new Rotamer();
  //cout << "after new Rotamer()" << endl;
  newConformer->deepcopy(tempConformer);
  //cout << "after deepcopy() " << endl;
  return newConformer;

  // Note: only reason to go through first creating a tempConformer object and then deepcopy it is because the Rotamer(ScreamAtomV, Rotconninfo*) constructor does not set allocatedScreamAtoms to true.  Though I can certainly set it to true here, it feels slightly less error prone this way because the less one needs to track memory the better.  The reason isn't that strong anyway because it is after all twice as first if I just do a Rotamer* newConformer = new Rotamer(atomlist, rotConnInfo) and then newConformer->allocatedScreamAtoms = true.

}

int Protein::rotamerClusterPlacement(Rotamer* cluster, MutInfo* mI) {
  /* Rotamer*, MutInfo* reasons: so that the base case will work */

  /* Note: assumes RotamerCluster and ClusterMutInfo have the exact same tree structure. */

  int mutationFlag = 0;

  vector<MutInfo*> mI_vector = mI->getAllMutInfos();
  vector<Rotamer*> rot_vector = cluster->getAllRotamers();

  assert(mI_vector.size() == rot_vector.size());
  
  for (int i = 0; i != mI_vector.size(); ++i) {
    string chn = mI_vector[i]->getChn();
    int pstn = mI_vector[i]->getPstn();

    if (mI_vector[i]->getRotConnInfo() == NULL) {
      mutationFlag += this->ntrlRotamerPlacement(chn, pstn, (AARotamer*)rot_vector[i] );
    } else {
      this->conformerPlacement(rot_vector[i], mI_vector[i]->getRotConnInfo());
    }
  }

  if (mutationFlag) 
    return 1;
  else 
    return 0;

}

void Protein::setRotamerClusterEmptyLatticeEnergy(Rotamer* cluster, MutInfo* mI, double E) {
  /* Rotamer*, MutInfo* reasons: so that the base case will work */

  /* Note: assumes RotamerCluster and ClusterMutInfo have the exact same tree structure. */

  vector<MutInfo*> mI_vector = mI->getAllMutInfos();
  vector<Rotamer*> rot_vector = cluster->getAllRotamers();

  assert(mI_vector.size() == rot_vector.size());
  
  for (int i = 0; i != mI_vector.size(); ++i) {
    string chn = mI_vector[i]->getChn();
    int pstn = mI_vector[i]->getPstn();

    if (mI_vector[i]->getRotConnInfo() == NULL) {
      //      this->ntrlRotamerPlacement(chn, pstn, (AARotamer*)rot_vector[i] );
      this->setEmptyLatticeEnergy(chn, pstn, rot_vector[i]->get_empty_lattice_E_abs());
      this->setPreCalcEnergy(chn, pstn, rot_vector[i]->get_preCalc_TotE() );
    } else {
      // will add back this line later on; when necessary
      //      this->conformerPlacement(rot_vector[i], mI_vector[i]->getRotConnInfo());
    }
  }

}

double Protein::getRotamerClusterEmptyLatticeEnergy(Rotamer* cluster, MutInfo* mI) {
  /* Rotamer*, MutInfo* reasons: so that the base case will work */

  /* Note: assumes RotamerCluster and ClusterMutInfo have the exact same tree structure. */

  vector<MutInfo*> mI_vector = mI->getAllMutInfos();
  vector<Rotamer*> rot_vector = cluster->getAllRotamers();

  assert(mI_vector.size() == rot_vector.size());
  double E=0;
  for (int i = 0; i != mI_vector.size(); ++i) {
    string chn = mI_vector[i]->getChn();
    int pstn = mI_vector[i]->getPstn();

    if (mI_vector[i]->getRotConnInfo() == NULL) {
      //this->ntrlRotamerPlacement(chn, pstn, (AARotamer*)rot_vector[i] );
      E+=this->getEmptyLatticeEnergy(chn, pstn);
    } else {
      //this->conformerPlacement(rot_vector[i], mI_vector[i]->getRotConnInfo());
      // do nothing currently; will modify this later on.
    }
  }
  return E;
}


int Protein::mutation(string chn, int pstn, string AA_3_letter_name) {
  // Calls ntrlRotamerPlacement.  Memory management done in ntrlRotamerPlacement.

  cout << "Performing mutation on chain " << chn << " position " << pstn << " to " << AA_3_letter_name << endl;

  string one_letter_AA = scream_tools::one_letter_AA(AA_3_letter_name);
  //vcvicek
  char * SCREAM_NEW_OLDLIB = getenv("SCREAM_NEW_OLDLIB");
  if (SCREAM_NEW_OLDLIB == NULL) {printf("error: enviromental variable SCREAM_NEW_OLDLIB is not set \n"); exit(1);}
  
  string default_lib_path = string(SCREAM_NEW_OLDLIB);
  //string library_file = default_lib_path + "G" + "/" + "G" + "_50.lib";
  string library_file = default_lib_path + one_letter_AA + "/" + one_letter_AA + "_50.lib";

  NtrlAARotlib library(library_file);
  library.reset_pstn();
  AARotamer *rot = library.get_next_rot();

  this->ntrlRotamerPlacement(chn, pstn, rot);

  return 1;


}

void Protein::setPreCalcEnergy(string chn, int pstn, double E) {
  this->get_AAChain(chn)->operator[](pstn)->get_rot()->set_preCalc_TotE(E);
  
}

double Protein::getPreCalcEnergy(string chn, int pstn) const {
  return this->get_AAChain(chn)->operator[](pstn)->get_rot()->get_preCalc_TotE();
}

void Protein::setEmptyLatticeEnergy(string chn, int pstn, double E) {
  this->get_AAChain(chn)->operator[](pstn)->get_rot()->set_empty_lattice_E(E);
}

double Protein::getEmptyLatticeEnergy(string chn, int pstn) const {
  return this->get_AAChain(chn)->operator[](pstn)->get_rot()->get_empty_lattice_E();
}


void Protein::setSideChainLibraryName(string chn, int pstn, string res) {
  ScreamAtomV atom_v = this->get_sc_atoms(chn, pstn);
  for (ScreamAtomVItr itr = atom_v.begin(); itr != atom_v.end(); ++itr) {
    (*itr)->library_name = res;
  }
}

void Protein::setProteinLibraryName(string res) {
  vector<ProteinComponent*>::iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    AAChain* tmp_chn = dynamic_cast<AAChain*>(*itr);      // temp = NULL if cannot be casted.
    if (tmp_chn != NULL) {
      string chn = tmp_chn->get_chain_desig();
      vector<int> residue_numbers = tmp_chn->getResidueNumbers();
      for (vector<int>::iterator i = residue_numbers.begin(); i != residue_numbers.end(); ++i) {
	this->setSideChainLibraryName(chn, *i, res);
      }
    }
  }

}

void Protein::fix_toggle(bool value) {

  vector<ProteinComponent*>::iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    (*itr)->fix_toggle(value);
  }

}


void Protein::fix_sc_toggle(bool value, int res_n, string chn) {

  AAChain* chain = this->get_AAChain(chn);
  chain->fix_toggle_sc_pstn(value, res_n);

}


void Protein::resetFlags() {

  for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr)
    (*itr)->resetFlag();
}


double Protein::sc_clash(string chain, int pstn) const {

  ScreamAtomV allAtomsInSystem;

  /* currently: entire list of atom_v is passed in */
  //  if (atom_v.size() != 0) {
  //    return get_AAChain(chain)->operator[](pstn)->get_rot()->get_sc()->worst_clash_distance(this->atom_v);
  //  } else {
  /* below: get all atoms.  This algorithm is obsolete because now protein NECESSARILY contains a copy/pointer to all atoms in system. */
  vector<ProteinComponent*>::const_iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    ScreamAtomV atomsInComponent = (*itr)->getAtomList();
    allAtomsInSystem.insert(allAtomsInSystem.end(), atomsInComponent.begin(), atomsInComponent.end());
  }


  /* First gather some info about this sidechain.  */
  SCREAM_ATOM* sample_atom;
  chain =  scream_tools::strip_whitespace(chain);
  for ( ScreamAtomVItr itr = allAtomsInSystem.begin(); itr != allAtomsInSystem.end(); ++itr) {
    if (  scream_tools::strip_whitespace((*itr)->chain) == chain and (*itr)->resNum == pstn ) {
      sample_atom = *itr;
      break;
    }
  }
  bool res_is_PRO = false;
  if (scream_tools::strip_whitespace(sample_atom->resName) == "PRO") {
    res_is_PRO = true;
  }


  /* need to remove sidechain currently being tested for clashes atoms */
  ScreamAtomV allButCurrentSC;
  
  //cout << "Chain:::" << chain << ":::" << endl;
  //cout << "Pstn:::" << pstn << ":::" << endl;

  for ( ScreamAtomVItr itr = allAtomsInSystem.begin(); itr != allAtomsInSystem.end(); ++itr) {
    // remove sidechain and CA that connects to CB of current residue sidechain.
    if ( scream_tools::strip_whitespace((*itr)->chain) == chain and (*itr)->resNum == pstn) {  // if on chain, pstn.
      if ( scream_tools::is_SC_atom((*itr)->atomLabel) or scream_tools::strip_whitespace((*itr)->atomLabel) == "CA")  {
	continue; 	// bypass
      }
      // getting rid of 1-3 interactions: N(i)-CB(i) and C(i)-CB(i)
      if ( scream_tools::strip_whitespace((*itr)->atomLabel) == "N" 
	   or scream_tools::strip_whitespace((*itr)->atomLabel) == "N" 
	   or scream_tools::strip_whitespace((*itr)->atomLabel) == "C") {
	continue; // bypass
      }
      // special case: Proline.  N conected to CD.
      if ( scream_tools::strip_whitespace((*itr)->resName) == "PRO" and 
	   (( scream_tools::strip_whitespace((*itr)->atomLabel) == "N") or 
	    ( scream_tools::strip_whitespace((*itr)->atomLabel) == "NT")) ) {   // Remark: also redundant.  1-3 interactions already took care of N.  But nature different, included here for clarity.
	continue; // bypass
      }

    }

    // more proline special case.  1-3 interaction between C of i-1 and CD of PRO.  Must find better way to place Proline!
    if (res_is_PRO) {
      if (  scream_tools::strip_whitespace((*itr)->chain) == chain and (*itr)->resNum == pstn -1 ) { 
	if ( scream_tools::strip_whitespace((*itr)->atomLabel) == "C") {
	  continue; // bypass.
	}
      }
      
    } 

    //else {
    // heavy atoms only
    if ( scream_tools::strip_whitespace((*itr)->atomLabel.substr(0,1)) != "H") {
      //(*itr)->dump();
      allButCurrentSC.push_back(*itr);
    }

    //}
  }
  //  cout << "Size of allButCurrentSC is " << allButCurrentSC.size() << endl;
  return get_AAChain(chain)->operator[](pstn)->get_rot()->get_sc()->worst_clash_distance(allButCurrentSC);
}

double Protein::sc_clash(string chain, int pstn, ScreamAtomV& atomsPassedIn) const {
  /* Returns the worst heavy atom distance from side chain specified by chain and pstn to atoms in atomsPassedIn */
  
  return get_AAChain(chain)->operator[](pstn)->get_rot()->get_sc()->worst_clash_distance(atomsPassedIn);
  
}

double Protein::conformer_clash(RotConnInfo* rotConnInfo) const {
  /* returns clash from of the variable region/side_chain atoms of this conformer to the rest of the protein. */

  /* First figure out which atoms are part of the conformer sidechain.*/
  vector<int> sideChain_heavyatoms_n;
  ScreamAtomV sideChain_heavyatoms;
  for (vector<int>::const_iterator itr = rotConnInfo->side_chain_atoms.begin();
       itr != rotConnInfo->side_chain_atoms.end(); ++itr) {
    int ptn_n = rotConnInfo->atom_n_map[*itr];
    SCREAM_ATOM* atom = this->getAtom(ptn_n);
    if ( scream_tools::strip_whitespace(atom->atomLabel.substr(0,1)) != "H") {
      sideChain_heavyatoms_n.push_back(ptn_n);      
      sideChain_heavyatoms.push_back(atom);// hope this works.  these pointers things always messed me up.
    }
  }
  //cout << "Now dumping sideChain_heavyatoms" << endl;
  //  for_each(sideChain_heavyatoms.begin(), sideChain_heavyatoms.end(), dump);
  //  cout << "end......" << endl;

  /* Secondly figure out which atoms belong to the non-conformer protein system. */
  vector<int> connection_point_atoms_ptn_n;
  for (vector<int>::const_iterator itr = rotConnInfo->connection_point_atoms.begin(); 
       itr != rotConnInfo->connection_point_atoms.end(); ++itr) {
    connection_point_atoms_ptn_n.push_back(rotConnInfo->atom_n_map[*itr]);
  }

  ScreamAtomV allButConformerSC_HeavyAtoms;
  allButConformerSC_HeavyAtoms.clear();
  for ( ScreamAtomVConstItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
    int n = (*itr)->n;
    /* insert non-conformer SC atoms into allButConformerSC ScreamAtomList. */
    vector<int>::const_iterator where_sc = find(sideChain_heavyatoms_n.begin(), sideChain_heavyatoms_n.end(), n);
    vector<int>::const_iterator where_connection_point = find(connection_point_atoms_ptn_n.begin(), 
							connection_point_atoms_ptn_n.end(), n);
    if (where_sc == sideChain_heavyatoms_n.end() and where_connection_point == connection_point_atoms_ptn_n.end() )  { // n is not a Conformer sc atom AND not a connection point (like CA)  ===> include in allButConformerSC
      
      if ( scream_tools::strip_whitespace((*itr)->atomLabel.substr(0,1)) != "H") {
	allButConformerSC_HeavyAtoms.push_back(*itr);
      }
    }
  }
  
  /* Finally calculate worst clash distances. */
  double worst_clash_dist = 999999;
  SCREAM_ATOM ** worst_clash_atom_on_Protein_p_tmp = &(*(sideChain_heavyatoms.begin()));
  SCREAM_ATOM* worst_clash_atom_on_Protein = NULL;
  SCREAM_ATOM* worst_clash_atom_on_conformer = NULL;
  //cout << "Dumping allButConformerSC_HeavyAtoms" << endl;
  //for_each(allButConformerSC_HeavyAtoms.begin(), allButConformerSC_HeavyAtoms.end(), dump);

  for (ScreamAtomVItr itr = sideChain_heavyatoms.begin(); itr != sideChain_heavyatoms.end(); ++itr) {
    double tmp_dist = (*itr)->worst_clash_dist(allButConformerSC_HeavyAtoms, worst_clash_atom_on_Protein_p_tmp);
    if (tmp_dist < worst_clash_dist) {
      //      cout <<"Updateing"<< endl;
      worst_clash_atom_on_conformer = *itr;
      worst_clash_atom_on_Protein = *worst_clash_atom_on_Protein_p_tmp;
      worst_clash_dist = tmp_dist;
    }
  }


  assert(worst_clash_atom_on_conformer != NULL);
  assert(worst_clash_atom_on_Protein != NULL);
  assert(*worst_clash_atom_on_Protein_p_tmp != NULL);
  if (worst_clash_atom_on_conformer != NULL and worst_clash_atom_on_Protein != NULL) {
    cout << "Worst clashing atom pair.  First, on confomer: " << endl;
    worst_clash_atom_on_conformer->dump();
    cout << "Now, the clashing atom on protein: " << endl;
    worst_clash_atom_on_Protein->dump();
  }
  return worst_clash_dist;

}

double Protein::conformer_clash(RotConnInfo* rotConnInfo, ScreamAtomV&) const {
  
}


//double Protein::sc_CRMS(string, int, Rotamer*) const {
  
//  return 0;

//}

double Protein::sc_CRMS(string chn, int pstn, Protein* other_ptn) const {
  
  ScreamAtomV self_atom_list = this->get_sc_atoms(chn, pstn);
  ScreamAtomV other_atom_list = other_ptn->get_sc_atoms(chn, pstn);
  
  ScreamAtomV self_heavy_atom_list = scream_tools::return_heavy_atoms(self_atom_list);
  ScreamAtomV other_heavy_atom_list = scream_tools::return_heavy_atoms(other_atom_list);

  // Turn them into maps.
  map<string, SCREAM_ATOM*> self_heavy_atom_map, other_heavy_atom_map;
  for (ScreamAtomVConstItr itr = self_heavy_atom_list.begin(); itr != self_heavy_atom_list.end(); ++itr) {
    string atom_label = scream_tools::strip_whitespace( (*itr)->atomLabel );
    self_heavy_atom_map.insert(make_pair(atom_label, *itr));
  }
  for (ScreamAtomVConstItr itr = other_heavy_atom_list.begin(); itr != other_heavy_atom_list.end(); ++itr) {
    string atom_label = scream_tools::strip_whitespace( (*itr)->atomLabel );
    other_heavy_atom_map.insert(make_pair(atom_label, *itr));
  }

  // Calculate unflipped CRMS.

  double CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
  double flipped_CRMS = CRMS;

  // Now deal with flipped residue cases: Asp, Glu, Phe, Tyr, Arg.  Flip:

  string res_name = scream_tools::strip_whitespace(self_atom_list[0]->resName);
  if (res_name == "PHE" or res_name == "TYR") {
    // CD1, CD2; CE1, CE2
    SCREAM_ATOM *CD1, *CD2, *CE1, *CE2;
    CD1 = other_heavy_atom_map["CD1"];
    CD2 = other_heavy_atom_map["CD2"];
    CE1 = other_heavy_atom_map["CE1"];
    CE2 = other_heavy_atom_map["CE2"];
    // flipping
    other_heavy_atom_map["CD1"] = CD2;
    other_heavy_atom_map["CD2"] = CD1;
    other_heavy_atom_map["CE1"] = CE2;
    other_heavy_atom_map["CE2"] = CE1;

    flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
  }
  else if (res_name == "ASP") {
    // OD1, OD2
    SCREAM_ATOM *OD1, *OD2;
    OD1 = other_heavy_atom_map["OD1"];
    OD2 = other_heavy_atom_map["OD2"];
    // flipping
    other_heavy_atom_map["OD1"] = OD2;
    other_heavy_atom_map["OD2"] = OD1;

    flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
  }
  else if (res_name == "GLU") {
    // OE1, OE2
    SCREAM_ATOM *OE1, *OE2;
    OE1 = other_heavy_atom_map["OE1"];
    OE2 = other_heavy_atom_map["OE2"];
    // flipping
    other_heavy_atom_map["OE1"] = OE2;
    other_heavy_atom_map["OE2"] = OE1;
    flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
  }
  else if (res_name == "ARG") {
    // NH1, NH2
    SCREAM_ATOM *NH1, *NH2;
    NH1 = other_heavy_atom_map["NH1"];
    NH2 = other_heavy_atom_map["NH2"];
    // flipping
    other_heavy_atom_map["NH1"] = NH2;
    other_heavy_atom_map["NH2"] = NH1;
    flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
  }

  CRMS = (CRMS < flipped_CRMS) ? CRMS : flipped_CRMS;
  return CRMS;

}

double Protein::conformer_CRMS(Protein*, RotConnInfo*) const {
  return 0;
}


pair<double, string> Protein::max_atom_dist_on_SC(string chn, int pstn, Protein* other_ptn) const {

  ScreamAtomV self_atom_list = this->get_sc_atoms(chn, pstn);
  ScreamAtomV other_atom_list = other_ptn->get_sc_atoms(chn, pstn);
  
  ScreamAtomV self_heavy_atom_list = scream_tools::return_heavy_atoms(self_atom_list);
  ScreamAtomV other_heavy_atom_list = scream_tools::return_heavy_atoms(other_atom_list);

  
  // Turn them into maps.
  map<string, SCREAM_ATOM*> self_heavy_atom_map, other_heavy_atom_map;
  for (ScreamAtomVConstItr itr = self_heavy_atom_list.begin(); itr != self_heavy_atom_list.end(); ++itr) {
    string atom_label = scream_tools::strip_whitespace( (*itr)->atomLabel );
    self_heavy_atom_map.insert(make_pair(atom_label, *itr));
  }
  for (ScreamAtomVConstItr itr = other_heavy_atom_list.begin(); itr != other_heavy_atom_list.end(); ++itr) {
    string atom_label = scream_tools::strip_whitespace( (*itr)->atomLabel );
    other_heavy_atom_map.insert(make_pair(atom_label, *itr));
  }

  // Calculate unflipped MAX dist.
  pair<double, string> Max_Dist_unflipped = scream_tools::max_equivalent_atom_dist(self_heavy_atom_map, other_heavy_atom_map);
  pair<double, string> Max_Dist_flipped = Max_Dist_unflipped;
  
  // Now deal with Flippants.

  string res_name = scream_tools::strip_whitespace(self_atom_list[0]->resName);
  if (res_name == "PHE" or res_name == "TYR") {
    // CD1, CD2; CE1, CE2
    SCREAM_ATOM *CD1, *CD2, *CE1, *CE2;
    CD1 = other_heavy_atom_map["CD1"];
    CD2 = other_heavy_atom_map["CD2"];
    CE1 = other_heavy_atom_map["CE1"];
    CE2 = other_heavy_atom_map["CE2"]; 
    // flipping
    other_heavy_atom_map["CD1"] = CD2;
    other_heavy_atom_map["CD2"] = CD1;
    other_heavy_atom_map["CE1"] = CE2;
    other_heavy_atom_map["CE2"] = CE1;

    Max_Dist_flipped = scream_tools::max_equivalent_atom_dist(self_heavy_atom_map, other_heavy_atom_map);
  }
  else if (res_name == "ASP") {
    // OD1, OD2
    SCREAM_ATOM *OD1, *OD2;
    OD1 = other_heavy_atom_map["OD1"];
    OD2 = other_heavy_atom_map["OD2"];
    // flipping
    other_heavy_atom_map["OD1"] = OD2;
    other_heavy_atom_map["OD2"] = OD1;

    Max_Dist_flipped = scream_tools::max_equivalent_atom_dist(self_heavy_atom_map, other_heavy_atom_map);
  }
  else if (res_name == "GLU") {
    // OE1, OE2
    SCREAM_ATOM *OE1, *OE2;
    OE1 = other_heavy_atom_map["OE1"];
    OE2 = other_heavy_atom_map["OE2"];
    // flipping
    other_heavy_atom_map["OE1"] = OE2;
    other_heavy_atom_map["OE2"] = OE1;
    Max_Dist_flipped = scream_tools::max_equivalent_atom_dist(self_heavy_atom_map, other_heavy_atom_map);
  }
  else if (res_name == "ARG") {
    // NH1, NH2
    SCREAM_ATOM *NH1, *NH2;
    NH1 = other_heavy_atom_map["NH1"];
    NH2 = other_heavy_atom_map["NH2"];
    // flipping
    other_heavy_atom_map["NH1"] = NH2;
    other_heavy_atom_map["NH2"] = NH1;
    Max_Dist_flipped = scream_tools::max_equivalent_atom_dist(self_heavy_atom_map, other_heavy_atom_map);
  }

  // So there are two Max_Dist_flipped values--one of them is an artifact.  Pick the smaller one.
  pair<double, string> final_Max_Dist = (Max_Dist_unflipped.first < Max_Dist_flipped.first) ? Max_Dist_unflipped : Max_Dist_flipped;
  return final_Max_Dist;

}

double Protein::sc_atom_CRMS(int atom_n, Protein* other_ptn) const {

  SCREAM_ATOM* a1 = this->getAtom(atom_n);
  SCREAM_ATOM* a2 = other_ptn->getAtom(atom_n);
  if (a1 == NULL or a2 == NULL) {
    cerr << "atom number " << atom_n << " not found!" << endl;
    return -999;
  }
  if (a1->resName != a2->resName or a1->chain != a2->chain or a1->atomLabel != a2->atomLabel) {
    cerr << "The two atoms do not match!" << endl;
    return -999;
  }

  string chn = a1->chain;
  int pstn = a1->resNum;

  double unflipped_distance = a1->distance(a2);
  double flipped_distance = 9999;

  string res_name = scream_tools::strip_whitespace(a1->resName);
  string a1_atom_label = scream_tools::strip_whitespace(a1->atomLabel);
  //ScreamAtomV other_heavy_atom_list;
  map<string, SCREAM_ATOM*> self_heavy_atom_map, other_heavy_atom_map;


  if (res_name == "PHE" or res_name == "TYR" or res_name == "ASP" or res_name == "GLU" or res_name == "ARG") {
    ScreamAtomV other_atom_list = other_ptn->get_sc_atoms(chn, pstn);
    ScreamAtomV other_heavy_atom_list = scream_tools::return_heavy_atoms(other_atom_list);

    ScreamAtomV self_atom_list = this->get_sc_atoms(chn,pstn);
    ScreamAtomV self_heavy_atom_list = scream_tools::return_heavy_atoms(self_atom_list);

    for (ScreamAtomVConstItr itr = self_heavy_atom_list.begin(); itr != self_heavy_atom_list.end(); ++itr) {
      string atom_label = scream_tools::strip_whitespace( (*itr)->atomLabel );
      self_heavy_atom_map.insert(make_pair(atom_label, *itr));
    }

    for (ScreamAtomVConstItr itr = other_heavy_atom_list.begin(); itr != other_heavy_atom_list.end(); ++itr) {
      string atom_label = scream_tools::strip_whitespace( (*itr)->atomLabel );
      other_heavy_atom_map.insert(make_pair(atom_label, *itr));
    }

  } 

  // done with the non-flipping residues.  now to the flipping residues.
  SCREAM_ATOM* flipped_atom;
  if (res_name == "PHE" or res_name == "TYR") {
    
    if (a1_atom_label == "CD1" or a1_atom_label == "CD2" or a1_atom_label == "CE1" or a1_atom_label == "CE2") {
      double CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
      SCREAM_ATOM *CD1, *CD2, *CE1, *CE2;
      CD1 = other_heavy_atom_map["CD1"];
      CD2 = other_heavy_atom_map["CD2"];
      CE1 = other_heavy_atom_map["CE1"];
      CE2 = other_heavy_atom_map["CE2"]; 
      // flipping
      other_heavy_atom_map["CD1"] = CD2;
      other_heavy_atom_map["CD2"] = CD1;
      other_heavy_atom_map["CE1"] = CE2;
      other_heavy_atom_map["CE2"] = CE1;
    
      double flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);

      if (flipped_CRMS < CRMS) {
	if (a1_atom_label == "CD1") flipped_distance = a1->distance(CD2);
	if (a1_atom_label == "CD2") flipped_distance = a1->distance(CD1);
	if (a1_atom_label == "CE1") flipped_distance = a1->distance(CE2);
	if (a1_atom_label == "CE2") flipped_distance = a1->distance(CE1);
      }

    }
  }
  else if (res_name == "ASP") {
    if (a1_atom_label == "OD1" or a1_atom_label == "OD2") {
      double CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
      SCREAM_ATOM *OD1, *OD2;
      OD1 = other_heavy_atom_map["OD1"];
      OD2 = other_heavy_atom_map["OD2"];
      // flipping
      other_heavy_atom_map["OD1"] = OD2;
      other_heavy_atom_map["OD2"] = OD1;
      
      double flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
      
      if (flipped_CRMS < CRMS) {
	if (a1_atom_label == "OD1") flipped_distance = a1->distance(OD2);
	if (a1_atom_label == "OD2") flipped_distance = a1->distance(OD1);
      }
    }
  }
  else if (res_name == "GLU") {
    if (a1_atom_label == "OE1" or a1_atom_label == "OE2") {
      double CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
      SCREAM_ATOM *OE1, *OE2;
      OE1 = other_heavy_atom_map["OE1"];
      OE2 = other_heavy_atom_map["OE2"];
      // flipping
      other_heavy_atom_map["OE1"] = OE2;
      other_heavy_atom_map["OE2"] = OE1;

      double flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);

      if (flipped_CRMS < CRMS) {
	if (a1_atom_label == "OE1") flipped_distance = a1->distance(OE2);
	if (a1_atom_label == "OE2") flipped_distance = a1->distance(OE1);
      }
    }
  }
  else if (res_name == "ARG") {
    if (a1_atom_label == "NH1" or a1_atom_label == "NH2") {
      double CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
      SCREAM_ATOM *NH1, *NH2;
      NH1 = other_heavy_atom_map["NH1"];
      NH2 = other_heavy_atom_map["NH2"];
      // flipping
      other_heavy_atom_map["NH1"] = NH2;
      other_heavy_atom_map["NH2"] = NH1;
            
      double flipped_CRMS = scream_tools::distance_by_atom_label_map(self_heavy_atom_map, other_heavy_atom_map);
      if (flipped_CRMS < CRMS) {
	if (a1_atom_label == "NH1") flipped_distance = a1->distance(NH2);
	if (a1_atom_label == "NH2") flipped_distance = a1->distance(NH1);
      }
    }
  }

  return (flipped_distance < unflipped_distance) ? flipped_distance : unflipped_distance;
}


void Protein::print_Me() const {
  // obsolete
  this->append_to_filehandle(&cout);

  /*
  vector<ProteinComponent*>::const_iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    (*itr)->print_Me();
  }
  */
}

void Protein::print_ordered_by_n() const {
  //obsolete

  this->append_to_filehandle(&cout);

  /*
  vector<ProteinComponent*>::const_iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    (*itr)->print_ordered_by_n();
  }
  */
}

void Protein::append_to_filehandle(ostream* ostream_p) const {

  fflush(stdout);

  vector<ProteinComponent*>::const_iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    (*itr)->append_to_filehandle(ostream_p);
  }

}


void Protein::print_bgf_file(ostream* ostream_p) const {
  /* Now creates a bgf_handler object and prints it using functions defined there. */
  //  bgf_handler TEMP
  //  TEMP.printfile(this->atom_v, ostream_p);
  //bgf_handler::printfile(this->atom_v, ostream_p);

}

void Protein::InitDataStructures(ScreamAtomV& atom_list) {
  /* This routine does not make use of atom connectivity records.  Atom connectivitity should have been initiated by this point. */

  new_mapping.clear();

  /* First Copy atom_list over */
  // Philosophy: most memory management (initialization, destruction) is done at a high level (except for deleting atoms when mutating stuff), hence the entire atom_list needs to be pointed to, not copied.  11/10/04.

  this->ptr_to_atom_v = &atom_list;
  this->atom_v = *ptr_to_atom_v; // by reference; i.e. it's an alias for *ptr_to_atom_v.  Why here?  Should be pointer, this way don't have to search and replace.  Shouldn't be so lazy.
  //cout << "Address: in InitDataStructures after this->atom_v = *ptr_to_atom_v: " << &atom_v << endl;
  // old code (prior to 11/10/04)
  //this->atom_v.insert(atom_v.begin(), atom_list.begin(), atom_list.end());
  this->placementMethod = "Default";
  this->CreateCBParameters = vector<double>(4,0);
  
  /* Initilziing Variables */
  ScreamAtomVConstItr itr;

  cout << "Total atoms in this system: " << atom_list.size()  << endl;

  string this_resName, last_resName;
  string this_chain, last_chain;
  int this_resNum, last_resNum;
  string this_keyw, last_keyw;

  this_resName = "beginning";
  this_chain = "beginning";
  this_resNum = -1;

  ScreamAtomV atoms_in_component;
  itr = atom_list.end();

  while (true) {
    if (itr == atom_list.end() ) {                       // initial condition
      itr = atom_list.begin();
      
      atoms_in_component.push_back(*itr);
      last_resName = this_resName;    last_chain = this_chain;    last_resNum = this_resNum;    
      this_resName = (*itr)->resName;    this_chain = (*itr)->chain;    this_resNum = (*itr)->resNum;
      continue;

    } else {                                          
      //      (*itr)->dump();
      ++itr;			// Iterating.
      
      if (itr == atom_list.end()) { // End of entire atom_list passed in.  Instantiate structures.
	/* Initialize last ProteinComponent Structure */
	SCREAM_ATOM* atom = *(atoms_in_component.begin());
	
	ProteinComponent* new_component;
	if (atom->resName == "RES") {
	  new_component = new Ligand(atoms_in_component);
	  cout << "New SCREAM Ligand Initiated... " << endl;
	} else if (atom->resName == "HOH") {
	  new_component = new Water(atoms_in_component);
	  cout << "Water added... " << endl;
	} else if (atom->keyw != "HETATM") {
	  cout << "Initializing chain..." << endl;
	  new_component = new AAChain(atoms_in_component);
	  cout << "New SCREAM Chain Initialized... " << endl;
	} else {
	  new_component = new Hetatm(atoms_in_component);
	  cout << "New SCREAM HETATM list Initialized... " << endl;
	}

	pc_v.push_back(new_component);
	atoms_in_component.clear();

	//	cout << "right before break" << endl;
	break;
      } // end end-condition block

      last_resName = this_resName;    last_chain = this_chain;    last_resNum = this_resNum;

      this_resName = (*itr)->resName;    
      this_chain = (*itr)->chain;    
      this_resNum = (*itr)->resNum;

    } // Closing "else" of any atom but-the-first-atom-loop.

    /* below: main loop.  a new Protein Component is initiated when one of the conditionals is visited. */

    if (last_chain != this_chain) {                   // the "canonical" case. different chain name guarantees a new protein componentn is encoutnered.
      SCREAM_ATOM* atom = *(atoms_in_component.begin());
      ProteinComponent* new_component;
      
      cout << "adding new component..." << endl;  flush(cout);
  
      if (atom->resName == "RES") {
	new_component = new Ligand(atoms_in_component);
	cout << "New SCREAM Ligand Initiated... " << endl;
      } else if (atom->resName == "HOH" or atom->resName == "TIP") {
	new_component = new Water(atoms_in_component);
	cout << "Water added... " << endl;

      } else if (atom->keyw != "HETATM") {
	new_component = new AAChain(atoms_in_component);
	cout << "Added new chain... " << endl; flush(cout);

      } else {
	new_component = new Hetatm(atoms_in_component);
	cout << "New SCREAM HETATM list Initialized... " << endl;
      }
      pc_v.push_back(new_component);
      //      cout << "atoms_in_component.size() is " << atoms_in_component.size() << endl;
      atoms_in_component.clear();


    } else {			// If chain designation same, need other criteria to determine the type of component encountered.
      SCREAM_ATOM* atom = *(atoms_in_component.begin());
      ProteinComponent* new_component;

      if (last_resName != this_resName and this_chain == " ") {            // ligand.
	if (last_resName == "RES") {
	  new_component = new Ligand(atoms_in_component);
	  cout << "New SCREAM Ligand Initiated... " << endl;

	  pc_v.push_back(new_component);
	  atoms_in_component.clear();
	} 
      } 
      else if (last_resNum != this_resNum and this_chain == " ") {         // water.
	if (last_resName == "HOH") {
	  new_component = new Water(atoms_in_component);
	  cout << "Water added... " << endl;

	  pc_v.push_back(new_component);
	  atoms_in_component.clear();
	} else if (scream_tools::strip_whitespace(atom->keyw) == "ATOM" and 
		   (*(atoms_in_component.begin()))->keyw == "HETATM") {       // HETATM's.
	  new_component = new Hetatm(atoms_in_component);
	  cout << "New SCREAM HETATM list Initialized... " << endl;
	  pc_v.push_back(new_component);
	  atoms_in_component.clear();
	}
      }
    }

    atoms_in_component.push_back(*itr);

    if (itr == atom_list.end()) {  // should never reach here.
      break;
    }
  } // end atom_v list while loop
  cout << "New Protein created with ";
  cout << pc_v.size();
  cout << " components." << endl;
  
}

void Protein::make_bond_from_connect_info(const vector<string>& connect_info_strv) {
  
  map<int, SCREAM_ATOM*> atom_n_SCREAM_ATOM_map;
  ScreamAtomVConstItr itr_scatom;
  for (itr_scatom = this->atom_v.begin(); itr_scatom != atom_v.end(); ++itr_scatom) {
    atom_n_SCREAM_ATOM_map.insert(make_pair((*itr_scatom)->n, (*itr_scatom)));
  }

  vector<string>::const_iterator itr;
  for (itr = connect_info_strv.begin(); itr != connect_info_strv.end(); ++itr) {
    vector<string> fields = scream_tools::split(*itr);
    if ( (*(fields.begin())) == string("ORDER")) {
      continue;
    }

    bool first_atom_FLAG = true;
    SCREAM_ATOM* cur_atom = NULL;
    
    for (vector<string>::iterator itr_s = fields.begin(); itr_s != fields.end(); ++itr_s) {
      if ( *itr_s == string("CONECT") ) {
	continue;
      } 

      if (first_atom_FLAG) {
	int cur_atom_n = atoi((*itr_s).c_str());
	cur_atom = atom_n_SCREAM_ATOM_map.find(cur_atom_n)->second;
	first_atom_FLAG = false;
	continue;
      }
      int atom_n = atoi((*itr_s).c_str());
      //      cout << "in sc_Protein::make_bond_from_connect_info atom_n is:::" << atom_n << endl;
      SCREAM_ATOM* atom_to_be_connected = atom_n_SCREAM_ATOM_map.find(atom_n)->second;
      //      cout << atom_to_be_connected->n << endl;
      cur_atom->make_bond(atom_to_be_connected);

    }
  }

}

void Protein::pdb_append_to_filehandle(ostream* ostream_p) const {

  vector<ProteinComponent*>::const_iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    (*itr)->pdb_append_to_filehandle(ostream_p);
  }
}

/* Helper functions */

class cmp_BGF_ptn_atom_ordering {
 public:
  bool operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
    bool order = false;
    // Comparision order: 0 means a1 < a2, 1 means a2 > a1.
    // 1. chain name
    // " " (space) is < "A", but want A > " " (space).
    if (a1->chain == " " and a2->chain != " ") { order = false; return order; }  
    else if (a2->chain == " " and a1->chain != " ") { order = true; return order;}
    else if (a1->chain < a2->chain) { order = true; return order; }
    else if (a1->chain > a2->chain) { order = false; return order; }
    // below: cases where a1->chain == a2->chain
    // 2. res_pstn
    else if ( a1->resNum != a2->resNum) { order = (a1->resNum < a2->resNum) ? true : false ;}
    // if the two atoms have same resName
    else if ( scream_tools::strip_whitespace(a1->resName) == scream_tools::strip_whitespace(a2->resName) ) {
	// 3.1 if res_name defined (defined means one of the 20 natural AMINO acids)
      if (scream_tools::is_natural_AA(a1->resName)) {
	return scream_tools::AA_atom_order(a1, a2);
      }
      // 3.2 maintain current order if res_name not defined (meaning not being of of the 20 natural AMINO acids)
      else // i.e. not a natural 
	{
	  //a1->dump(); a2->dump();
	  return (a1->n < a2->n);  // maintain original ordering
	}
      
    }
    else {  // if the above rules don't decide anything, use what's already in place.
      return (a1->n < a2->n);
    }
    return order;
  }
};

void Protein::_fix_entire_atom_list_ordering() {
  /* This function reorders atoms in the atom list so to conform to the BGF file protein atom ordering format. */
  //cout << "Before sorting: " << endl;
  //for_each (this->atom_v.begin(), this->atom_v.end(), dump);

  sort(this->atom_v.begin(), this->atom_v.end(), cmp_BGF_ptn_atom_ordering());  // very slow
  //cout << "After sorting: " << endl;
  //for_each(this->atom_v.begin(), this->atom_v.end(), dump);
  int c = 1;
  for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
    (*itr)->n = c;
    ++c;
  }
}

void Protein::_fix_residue_in_atom_list_ordering(string chn, int pstn) {
  // renumberring--then, all atoms after sidechain.
  
  //clock_t start = clock();

  SCREAM_ATOM* first_residue_ATOM, *last_residue_ATOM;

  int first_res_n = 99999999;
  int last_res_n = -1;

  ScreamAtomVItr first_residue_ATOM_itr, last_residue_ATOM_itr;

  for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
    if ( (*itr)->chain == chn and (*itr)->resNum == pstn) {
      int n = (*itr)->n;
      if (n < first_res_n) {
	first_res_n = n;
	first_residue_ATOM_itr = itr;
      }
      if (n> last_res_n) {
	last_res_n = n;
	last_residue_ATOM_itr = itr;
      }
    }
  }
  ++last_residue_ATOM_itr;

  // needed for new_mapping; 9-13-05
  int i = 0;
  map<SCREAM_ATOM*, int> old_residue_ordering_map; old_residue_ordering_map.clear();
  for (ScreamAtomVItr itr = first_residue_ATOM_itr; itr != last_residue_ATOM_itr; ++itr, ++i) {
    old_residue_ordering_map[*itr] = i; // original mapping
    //    (*itr)->dump();
    //    cout << "original mapping" << i << endl;
  }

  sort(first_residue_ATOM_itr, last_residue_ATOM_itr, cmp_BGF_ptn_atom_ordering());
  //sort(this->atom_v.begin(), this->atom_v.end(), cmp_BGF_ptn_atom_ordering());

  i = 0;
  map<SCREAM_ATOM*, int> new_residue_ordering_map; new_residue_ordering_map.clear();
  for (ScreamAtomVItr itr = first_residue_ATOM_itr; itr != last_residue_ATOM_itr; ++itr, ++i) {
    new_residue_ordering_map[*itr] = i; // new mapping
    //(*itr)->dump();
    //    cout << "new mapping" << i << endl;
  }

  map<int, int> after_ordering_mapping; after_ordering_mapping.clear();
  for (map<SCREAM_ATOM*, int>::const_iterator itr = old_residue_ordering_map.begin();
       itr != old_residue_ordering_map.end(); ++itr) {
    SCREAM_ATOM* atom = itr->first;
    int atom_n = itr->second;

    after_ordering_mapping[atom_n] = new_residue_ordering_map[atom];
  }

  //  clock_t   end = clock();
  //cout << "  sorting took (total) " << (double) (end - start) / CLOCKS_PER_SEC << " seconds." << endl;

  /*int c = 0;
  for (ScreamAtomVItr itr = first_residue_ATOM_itr; itr != last_residue_ATOM_itr; ++itr) {
    (*itr)->n = first_res_n + c;
    ++c;
    }*/

  int c = 1;
  for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
    (*itr)->n = c;
    ++c;
  }

  int after_ordering_first_n = (*first_residue_ATOM_itr)->n;

  for (map<int, int>::iterator itr = after_ordering_mapping.begin(); itr != after_ordering_mapping.end(); ++itr) {
    int atom_n = after_ordering_first_n + itr->first;

    //    cout << "before: atom_n: " << atom_n << "  new_mapping[atom_n]" << this->new_mapping[atom_n] << endl;
    //this->atom_v[atom_n-1]->dump();

    if ( this->new_mapping[atom_n] == -1) {
      continue;
    }
    else {
      // Note to self: need to keep clear head what's being mapped to what... this portion has nothing to do with the pre-mutated atom_n.  ONLY dealing with atom_n after mutation, and thus the indices for this->new_mapping (pre-mutation atom_n) are irrelevent.  only the values contained (post-mutation atom_n) 
      int new_atom_n = after_ordering_first_n + itr->second;
      replace(this->new_mapping.begin(), this->new_mapping.end(), atom_n, new_atom_n);

    }
    //cout << " after: atom_n: " << atom_n << "  new_mapping[atom_n]" << this->new_mapping[atom_n] << endl;
  }
}


void Protein::_fix_charges() {
  /* THis function (somewhat of a temp function) fixes charges. */

  vector<ProteinComponent*>::const_iterator itr;
  for (itr = pc_v.begin(); itr != pc_v.end(); ++itr) {
    AAChain* temp = dynamic_cast<AAChain*>(*itr);      // temp = NULL if cannot be casted.
    // Assumes that an AAChain contains residue numbers that are in sequence.

    if (temp != NULL) {
      // fix chaarges
      vector<int> res_numbers = temp->getResidueNumbers();

      for (int i = 0; i < res_numbers.size() ; ++i) {
      //for (vector<int>::const_iterator itr = res_numbers.begin(); itr != res_numbers.end(); ++itr) {
	if (i==0)
	  (*temp)[res_numbers[i]]->assign_charges("CHARMM22", AminoAcid::NTERM);
	else if (i== (res_numbers.size() - 1) ) 
	  (*temp)[res_numbers[i]]->assign_charges("CHARMM22", AminoAcid::CTERM); // assigns new charges for the backbone only.
	else
	  (*temp)[res_numbers[i]]->assign_charges("CHARMM22", AminoAcid::MAINCHAIN);
      }
      
    }
    
    else {
      // do nothing
      continue;
    }

  }

}


/* Private Readibility improvement functions */

vector<int> Protein::_getAtomNumbersInProtein(const map<int, int>& atomMapping, const vector<int>& atomNumbersInConformer) const {
  /* returns list of atom numbers that cooresp*/
  
  vector<int> atomsInProtein;
  atomsInProtein.clear();

  for (vector<int>::const_iterator itr = atomNumbersInConformer.begin(); 
       itr != atomNumbersInConformer.end(); ++itr) {
    const int proteinAtomNum = atomMapping.find(*itr)->second;
    atomsInProtein.push_back(proteinAtomNum);
  }

  return atomsInProtein;

}

ScreamAtomV Protein::_mutationHelpers_alloc_new_sc_atoms(string chn, int pstn, AARotamer* AArot, ScreamAtomV& old_sc_atoms, map<SCREAM_ATOM*, SCREAM_ATOM*>& map_rot_to_new_sc_atoms) {
  //  cout << "in _mutationHelpers_alloc_new_sc_atoms" << endl;

  old_sc_atoms = this->get_sc_atoms(chn, pstn);
  ScreamAtomV rot_sc_atoms = AArot->get_sc_atoms();

  // first allocate memory and create a map between rotamer atoms and new atoms to facilitate connectivity issues.
  ScreamAtomV new_sc_atoms;
  
  for (ScreamAtomVConstItr itr = rot_sc_atoms.begin(); itr != rot_sc_atoms.end(); ++itr) {
    SCREAM_ATOM* atom = new SCREAM_ATOM();
    new_sc_atoms.push_back(atom);
    atom->copy(*(*itr)); // Careful--needs to modify connecitivity.--do this after this loop
    map_rot_to_new_sc_atoms[*itr] = atom;
  }

  // take care of CB position, if CB position exists in current backbone.
  SCREAM_ATOM* ptn_bb_CB = this->getAtom(MutInfo("A"+stringify(pstn) + "_" + chn), "CB");
  if (ptn_bb_CB != NULL) {
    for (ScreamAtomVItr itr = new_sc_atoms.begin(); itr != new_sc_atoms.end(); ++itr)
      if ( (*itr)->stripped_atomLabel == "CB" )
	for (int i=0; i<3; ++i)
	  (*itr)->x[i] = ptn_bb_CB->x[i];
  }
  else { // == NULL, i.e. a GLY.
    SCREAM_ATOM* new_CB;
    this->get_AAChain(chn)->get_aa_m()[pstn]->get_rot()->create_CB(this->CreateCBParameters, new_CB);
    for (ScreamAtomVItr itr = new_sc_atoms.begin(); itr != new_sc_atoms.end(); ++itr)
      if ( (*itr)->stripped_atomLabel == "CB" )
	for (int i=0; i<3; ++i)
	  (*itr)->x[i] = new_CB->x[i];
  }
  return new_sc_atoms;
}

void Protein::_mutationHelpers_connectivities_fix(string chn, int pstn, AARotamer* AArot, ScreamAtomV& new_sc_atoms, map<SCREAM_ATOM*, SCREAM_ATOM*>& map_rot_to_new_sc_atoms) {
  //  cout << "in _mutationHelpers_connectivities_fix" << endl;

  string old_res_name = this->get_AAChain(chn)->get_aa_m()[pstn]->get_CA()->resName;
  string new_res_name = AArot->get_resName();

  // fix internal connecitivity--intensive debugging needed
  for (ScreamAtomVConstItr new_atom_itr = new_sc_atoms.begin(); new_atom_itr != new_sc_atoms.end(); ++new_atom_itr) {

    map<SCREAM_ATOM*, int> new_connectivity_m;
    for (map<SCREAM_ATOM*, int>::const_iterator conn_itr = (*new_atom_itr)->connectivity_m.begin(); 
	 conn_itr != (*new_atom_itr)->connectivity_m.end(); ++conn_itr) {
      SCREAM_ATOM* rot_atom = conn_itr->first; // if base atom == CB and connected atom == CA, rot_atom == NULL
      SCREAM_ATOM* new_atom = map_rot_to_new_sc_atoms[rot_atom];
      if (new_atom == NULL) { // sometimes, CA atom if NULL.
	continue;
      }
      else 
	new_connectivity_m[new_atom] = conn_itr->second;
    }
    (*new_atom_itr)->connectivity_m.clear();
    (*new_atom_itr)->connectivity_m.insert(new_connectivity_m.begin(), new_connectivity_m.end());

  }
  
  // fix boundary connecitivity--CA, CB atom connecitivites atom on old atoms and rotamer atoms.
  // get CA atom on both protein backbone and rotamer.
  SCREAM_ATOM* ptn_CA = this->get_AAChain(chn)->get_aa_m()[pstn]->get_CA();
  SCREAM_ATOM* rot_CA = AArot->get_bb()->get("CA");
  SCREAM_ATOM* ptn_CB = this->get_AAChain(chn)->get_aa_m()[pstn]->get_CB();
  SCREAM_ATOM* rot_CB;
  if (new_res_name != "GLY" ) { rot_CB = AArot->get_sc()->get("CB"); } else {rot_CB = AArot->get_sc()->get("HCA");} 
  if (ptn_CB == NULL) { // i.e. if old_res is GLY
    ptn_CB = this->get_AAChain(chn)->get_aa_m()[pstn]->get_rot()->get_sc()->get("HCA");
  }

//   cout << "!!!!! ptn_CB connectivity_m size: " << ptn_CB->connectivity_m.size() << endl;
//   for (map<SCREAM_ATOM*, int>::iterator ii = ptn_CB->connectivity_m.begin();
//        ii != ptn_CB->connectivity_m.end(); ++ii) {
//     ii->first->dump(); 
//   }

  SCREAM_ATOM* new_CB = map_rot_to_new_sc_atoms[rot_CB];

//   //  if (scream_tools::strip_whitespace(new_CB->atomLabel) == "CB" and new_CB->resNum == 178) {
//     cout << " in sc_Protein ntrl placement: " << endl;
//     new_CB->dump();
//     cout << "  Size of connectivity_m of new_cb: " << new_CB->connectivity_m.size() << endl;
//     for (map<SCREAM_ATOM*, int>::iterator ii = new_CB->connectivity_m.begin(); 
// 	 ii != new_CB->connectivity_m.end(); ++ii) {
//       cout << ii->first << endl;
//       ii->first->dump(); 
//     }
//     //  }
  
  int order = ptn_CA->connectivity_m[ptn_CB];
  ptn_CA->connectivity_m.erase(ptn_CB);
  ptn_CA->connectivity_m[new_CB] = order;  // from protein CA to new sidechain CB

//   if (scream_tools::strip_whitespace(new_CB->atomLabel) == "CB" and new_CB->resNum == 178) {
//     cout << " before .erase: " << endl;
//     new_CB->dump();
//     cout << "  Size of connectivity_m of new_cb: " << new_CB->connectivity_m.size() << endl;
//   }


  
  new_CB->connectivity_m.erase(rot_CA);    // should actually be unnecessary.
  new_CB->connectivity_m[ptn_CA] = order;  // from new sidechain CB to protein CA.

//   if (scream_tools::strip_whitespace(new_CB->atomLabel) == "CB" and new_CB->resNum == 178) {
//     cout << " after .erase: " << endl;
//     new_CB->dump();
//     cout << "  Size of connectivity_m of new_cb: " << new_CB->connectivity_m.size() << endl;
//   }



}

void Protein::_mutationHelpers_PRO_connectivities_and_numbering_fix(string chn, int pstn, AARotamer* AArot, ScreamAtomV& new_sc_atoms, map<SCREAM_ATOM*, SCREAM_ATOM*>& map_rot_to_new_sc_atoms) {
  //  cout << "in _mutationHelpers_PRO_connectivities_and_numbering_fix" << endl;
  
  /* This routine fixes Proline connectivities and what-nots. */

  string old_res_name = this->get_AAChain(chn)->get_aa_m()[pstn]->get_CA()->resName;
  string new_res_name = AArot->get_resName();

  /* First, PRO--> non-PRO */
    
  if (old_res_name == "PRO") {  
    // first delete N--CD bond.
    SCREAM_ATOM* ptn_N = this->get_AAChain(chn)->get_aa_m()[pstn]->get_N();
    SCREAM_ATOM* old_ptn_CD;
    
    for (map<SCREAM_ATOM*, int>::iterator itr = ptn_N->connectivity_m.begin(); itr != ptn_N->connectivity_m.end(); ++itr) {
      SCREAM_ATOM* a = itr->first;
      if (scream_tools::strip_whitespace(a->atomLabel) == "CD") {
	old_ptn_CD = a;
	break;
      }
    }
    ptn_N->connectivity_m.erase(old_ptn_CD);

    // then add HN atom.

    // First, initialize the CA_i, N_i, and C_i_minus_1 atoms.

    SCREAM_ATOM* new_HN_atom = new SCREAM_ATOM();
    SCREAM_ATOM *CA_i, *N_i, *C_i_minus_1;

    CA_i = this->get_AAChain(chn)->operator[](pstn)->get_CA();
    N_i = this->get_AAChain(chn)->operator[](pstn)->get_N();
    AminoAcid* AA = this->get_AAChain(chn)->operator[](pstn-1);
    if (AA == NULL) {
      cout << "Mutating PRO to a non-Pro not supported for N-term PRO!  Exiting. " << endl;
      exit(2);
    }
    C_i_minus_1 = AA->get_rot()->get_bb()->get("C");


    SCREAM_ATOM* template_HN_atom =  this->get_AAChain(chn)->get_aa_m().begin()->second->get_rot()->get_bb()->get("HN"); // assumes that the peptides chain actually has a HN atom at N-term.
    
    new_HN_atom->copy(*template_HN_atom);

    // Then make the HN atom on the backbone.

    scream_tools::calc_new_HN_atom_coords(CA_i, N_i, C_i_minus_1, new_HN_atom);
    new_HN_atom->n = 99999;
    new_HN_atom->resNum = pstn;
    new_HN_atom->resName = new_res_name;
    
    // make HN-N bond.
    new_HN_atom->connectivity_m.clear();
    new_HN_atom->connectivity_m[N_i] = 1;  // from HN to N
    N_i->connectivity_m[new_HN_atom] = 1;  // from N to HN

    // insert into atom_list and renumber everything behind it.
    ScreamAtomVItr one_past_N_itr;
    int atom_N_num;
    for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      if ( (*itr)->chain == chn and (*itr)->resNum == pstn and scream_tools::strip_whitespace((*itr)->atomLabel) == "N") {
	atom_N_num = (*itr)->n;
	one_past_N_itr = itr;
	one_past_N_itr++;
	break;
      }
    }
    new_HN_atom->n = atom_N_num;  // temp value
    this->atom_v.insert(one_past_N_itr, new_HN_atom);
    for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      if ( (*itr)->n > atom_N_num) {
	(*itr)->n++;  // pushing everything back by 1.
      }
    }
    // now there's room for new_HN_atom->n.
    new_HN_atom->n++;
  }

  /* Then, the non-PRO to PRO case */

  if (new_res_name == "PRO") { 
    // first remove HN atom and bond.
    SCREAM_ATOM* ptn_N = this->get_AAChain(chn)->get_aa_m()[pstn]->get_N();
    SCREAM_ATOM* HN_to_be_removed;

    int HN_removed = 0;
    int removed_HN_atom_n = 0;
    for ( map<SCREAM_ATOM*, int>::iterator itr = ptn_N->connectivity_m.begin();
	  itr != ptn_N->connectivity_m.end(); ++itr) {
      if (scream_tools::strip_whitespace(itr->first->atomLabel) == "HN" and HN_removed == 0) {
	HN_removed++;
	HN_to_be_removed = itr->first;
	removed_HN_atom_n = HN_to_be_removed->n;
	break;
      }
    }
    ptn_N->connectivity_m.erase(HN_to_be_removed);

    ScreamAtomVItr to_be_removed_HN_itr = find(this->atom_v.begin(), this->atom_v.end(), HN_to_be_removed);
    this->atom_v.erase(to_be_removed_HN_itr);

    delete HN_to_be_removed;

    // renumber everthing behind the removed HN atom.
    for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      if ( (*itr)->n > removed_HN_atom_n) {
	(*itr)->n--;  // pushing everything back by 1.
      }
    }


    // add N--CD bond.

    SCREAM_ATOM* rot_CD = AArot->get_sc()->get("CD");
    SCREAM_ATOM* new_CD = map_rot_to_new_sc_atoms[rot_CD];
    
    ptn_N->make_bond(new_CD);



  }


}

void Protein::_mutationHelpers_GLY_connectivities_and_numbering_fix(string chn, int pstn, AARotamer* AArot, ScreamAtomV& new_sc_atoms, map<SCREAM_ATOM*, SCREAM_ATOM*>& map_rot_to_new_sc_atoms) {

  //  cout << "in _mutationHelpers_GLY_connectivities_and_numbering_fix" << endl;

  /* This routine fixes Glycine connectivities and what-nots. */

  string old_res_name = this->get_AAChain(chn)->get_aa_m()[pstn]->get_CA()->resName;
  string new_res_name = AArot->get_resName();

  /* First, GLY--> non-PRO */
  // will code this later

  if (old_res_name == "GLY") {

    SCREAM_ATOM* backbone_O = this->get_AAChain(chn)->get_aa_m()[pstn]->get_rot()->get_bb()->get("O");
    SCREAM_ATOM* L_amino_HCA = this->get_AAChain(chn)->get_aa_m()[pstn]->get_rot()->get_sc()->get("HCA");  // only HCA is here, in the "sidechain" of GLYCINE.

    ScreamAtomVItr backbone_O_itr, L_amino_HCA_itr;
    // find stuff.
    
    for (ScreamAtomVItr itr  = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      if ( *itr == L_amino_HCA ) { L_amino_HCA_itr = itr; };
    }

    this->atom_v.erase(L_amino_HCA_itr);

    // find again, since vector content has changed.
    for (ScreamAtomVItr itr  = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      if ( *itr == backbone_O ) { backbone_O_itr = itr; };
    }

    ++backbone_O_itr;  // move one past, could be .end()
    this->atom_v.insert(backbone_O_itr, L_amino_HCA);
    
    // Finished moving the L-aminoacid HCA atom behind C and O in atom_list.  This is only temporary; since GLY is about to be mutated to something else.  (only for placeholder convenice).
    // Now, renumber.  Must as well renumber from 1 to whatever the size of atom_v is.
    int c = 1;
    for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      (*itr)->n = c;
      ++c;
    }
  }

  /* non-GLY --> GLY */
  // code this first, since i need this code is more urgent.

  if (new_res_name == "GLY") { 
    // Most of the work done in AminoAcid contructor already.

  }

}


void Protein::_mutationHelpers_init_new_mapping(ScreamAtomV& atom_list) {

  this->new_mapping.clear();
  this->new_mapping.push_back(-1); // zeroth index set to -1.  placeholder.
  
  int i = 1;
  ScreamAtomVConstItr itr = atom_list.begin(); 

  // GLY, PRO should be taken care of already.
  for (; itr != atom_list.end(); ++itr, ++i) {
    this->new_mapping.push_back(i);
  }

}

void Protein::_mutationHelpers_insert_new_sc_atoms(string chn, int pstn, AARotamer* AArot, ScreamAtomV& new_sc_atoms, ScreamAtomV& old_sc_atoms) {
  Debug _debugInfo("Protein::_mutationHelpers_insert_new_sc_atoms");
  _debugInfo.out(" Starting _mutationHelpers_insert_new_sc_atoms");

  //  ScreamAtomV old_sc_atoms = this->get_sc_atoms(chn, pstn);
  ScreamAtomV rot_sc_atoms = AArot->get_sc_atoms();
  
  // get start and end sc pstns for renumbering.  insertion assumes that system has the corect atom ordering.
  int start_sc_pstn = 99999999;
  int end_sc_pstn = -1;
  for (ScreamAtomVConstItr itr = old_sc_atoms.begin(); itr != old_sc_atoms.end(); ++itr) {
    int n = (*itr)->n;
    start_sc_pstn = (n < start_sc_pstn) ? n : start_sc_pstn;
    end_sc_pstn = (n > end_sc_pstn) ? n : end_sc_pstn;
  }
  int new_end_sc_pstn = end_sc_pstn + (rot_sc_atoms.size() - old_sc_atoms.size() );

  // take care of new_mapping.  9-12-05.
  for (int i = 1; i <= this->new_mapping.size(); ++i) {
    if ( (start_sc_pstn <= i) and (i <= end_sc_pstn) ) {
      _debugInfo.out("Gonner atoms: " + itoa(i));
      this->new_mapping[i] = 0; // these atoms are gonner's.
    }
    if (i > end_sc_pstn) {
      this->new_mapping[i] = new_end_sc_pstn + ( i - end_sc_pstn );
    }
  }
  
  // renumbering--first, just atoms on new sidechain.
  
  SCREAM_ATOM* last_sc_ATOM;
  int n_count = -1;
  for (ScreamAtomVConstItr itr = old_sc_atoms.begin(); itr != old_sc_atoms.end(); ++itr) {
    last_sc_ATOM = ( (*itr)->n > n_count) ? *itr : last_sc_ATOM; // gets the last_sc_ATOM
    n_count = last_sc_ATOM->n;
  }
  int c = 0;
  for (ScreamAtomVConstItr itr = new_sc_atoms.begin(); itr != new_sc_atoms.end(); ++itr, ++c) {
    (*itr)->n = start_sc_pstn + c;
  }
  
    
  // renumberring--then, all atoms after sidechain.
  ScreamAtomVItr last_sc_ATOM_itr =  find(this->atom_v.begin(), this->atom_v.end(), last_sc_ATOM);
  
  ++last_sc_ATOM_itr; // could become .end()
  for (; last_sc_ATOM_itr != this->atom_v.end(); ++last_sc_ATOM_itr, ++c) {
    (*last_sc_ATOM_itr)->n = start_sc_pstn + c;
  }
  
  // then remove and insert atoms
  // inserting
  last_sc_ATOM_itr =  find(this->atom_v.begin(), this->atom_v.end(), last_sc_ATOM);
  ++last_sc_ATOM_itr; // one past
  this->atom_v.insert(last_sc_ATOM_itr, new_sc_atoms.begin(), new_sc_atoms.end());

  // erasing original residue
  SCREAM_ATOM* first_sc_ATOM = *(old_sc_atoms.begin());
  ScreamAtomVItr first_sc_ATOM_itr = find(this->atom_v.begin(), this->atom_v.end(), first_sc_ATOM);
  last_sc_ATOM_itr =  find(this->atom_v.begin(), this->atom_v.end(), last_sc_ATOM);
  ++last_sc_ATOM_itr;
  this->atom_v.erase(first_sc_ATOM_itr, last_sc_ATOM_itr);

  // new: March 4th.  fixed atom flags needs to be taken care of.  before old sc atoms get deleted.
  //  SCREAM_ATOM* sample_old_sc_atom = *(old_sc_atoms.begin());
  //  int old_sc_atom_flag = sample_old_sc_atom->flags;
  //  for (ScreamAtomVItr itr = new_sc_atoms.begin(); itr != new_sc_atoms.end(); ++itr) {
  //    (*itr)->flags = old_sc_atom_flag;
  //  }

  // Modification: 7-17-07.  fixed atom flags more for first atom could be a CB, which could be fixed or moveable depending on previous assignment.  Thus, need to survey atom flag info a little bit deeper.
  //  cout << "7-17-07 start" << endl;
  int CB_flag_value = -1;
  int other_flag_value = -1;
  for (ScreamAtomVItr sample_old_sc_atom = old_sc_atoms.begin(); sample_old_sc_atom != old_sc_atoms.end(); ++sample_old_sc_atom) {

    if ( (*sample_old_sc_atom)->stripped_atomLabel == "CB") // or (*sample_old_sc_atom)->stripped_atomLabel == "HCB")
      other_flag_value = (*sample_old_sc_atom)->flags;
    else
      other_flag_value = (*sample_old_sc_atom)->flags;
    if ( CB_flag_value != -1 and other_flag_value != -1)
      break;
  }
  //cout << "7-17-07 b start" << endl;

  for (ScreamAtomVItr itr = new_sc_atoms.begin(); itr != new_sc_atoms.end(); ++itr) {
    //    cout << ":::" << (*itr)->stripped_atomLabel << ":::" << endl;
    //    (*itr)->dump();
    if ( (*itr)->stripped_atomLabel == "CB") // or (*itr)->stripped_atomLabel == "HCB")
      (*itr)->flags = CB_flag_value;
    else
      (*itr)->flags = other_flag_value;
  }
  //  cout << "7-17-07 c end" << endl;
}


void Protein::_mutationHelpers_delete_old_sc_atoms(ScreamAtomV& old_sc_atoms) {
  //  cout << "in _mutationHelpers_delete_old_sc_atoms" << endl;
  for (ScreamAtomVItr itr = old_sc_atoms.begin(); itr != old_sc_atoms.end(); ++itr) {
    delete (*itr);
  }
}

void Protein::_mutationHelpers_renaming_mut_atoms(string chn, int pstn, string newResName, ScreamAtomV& new_sc_atoms) {
  //  cout << "in _mutationHelpers_renaming_mut_atoms" << endl;
    // first for the sidechain atoms
    for (ScreamAtomVItr itr = new_sc_atoms.begin(); itr != new_sc_atoms.end(); ++itr) {
      (*itr)->chain = chn;
      (*itr)->resNum = pstn;
    }
    // then for the backbone atoms
    for (ScreamAtomVItr itr = this->atom_v.begin(); itr != this->atom_v.end(); ++itr) {
      if ( (*itr)->chain == chn and (*itr)->resNum == pstn) {
	(*itr)->resName = newResName;
      }
    }

}
