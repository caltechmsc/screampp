/* sc_Protein.hpp
 *
 * Header file for classes relevant to proteins in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#ifndef SC_PROTEIN_HPP
#define SC_PROTEIN_HPP

#include "defs.hpp"
#include "MutInfo.hpp"

#include "sc_ProteinComponent.hpp"
#include "scream_atom.hpp"
#include "sc_AAChain.hpp"
#include "sc_Ligand.hpp"
#include "sc_Hetatm.hpp"
#include "sc_Water.hpp"
#include <vector>

#include "scream_rtf.hpp"

#include "AARotamer.hpp"

using namespace std;

/** Class Protein.  More accurately, a System.  Chose the name protein because SCREAM deals primarily with protein systems.
 *  
 */

class Protein {

public:

  /* Contructors and Destructors. */
  
  //vcvicek commented this one line: Protein();
  //  Protein(const Protein&);
  //Protein(ScreamAtomV&); // not a const ScreamAtomV because mutations will mess with the underlying atom_list structure.
  Protein(ScreamAtomV*); // now a pointer to a list of scream atoms.  will this work?  3/8/05.
  ~Protein();

  /* Get SCREAM_ATOM, scream component functions. */

  AAChain* get_AAChain(string) const;   ///< Returns AAChain with provided designation.
  Ligand* get_Ligand() const; 	///< Returns the ONLY ligand in this structure.  Temporary measure!  For some proteins may have more than just one ligand.
  ProteinComponent* operator[](int) const; ///< Returns pointer to ProteinComponent object.
  ProteinComponent* get_Component_with_ChainName(string) const; ///< Returns pointer to Protein Component object.
  //vcvicek Protein& Protein::operator=(const Protein&); ///< overrides default assignment operator.
  Protein& operator=(const Protein&); ///< overrides default assignment operator.

  string get_res_type(string, int) const; ///< Returns the residue type of chain string postn int.  3 letter abbreviation.
  ScreamAtomV get_sc_atoms(string, int) const;  ///< Returns specified SC atoms in a vector. 
  ScreamAtomV get_sc_atoms(MutInfo) const; ///< Returns sidechain atoms that correspond for a MutInfo object.  Recall MutInfo objects could stored multiple number of MutInfo's.
  ScreamAtomV get_variable_atoms(RotConnInfo*) const; ///< Returns variables atoms in a RotConnInfo.
  ScreamAtomV get_visible_in_EL_mutInfo_atoms(MutInfo&, RotConnInfo*) const; ///< Returns atoms associated to a particular mutInfo that are visible in empty lattice calculations.
  ScreamAtomV& getAtomList();              ///< Returns atom_list.  This returns the actual vector<SCREAM_ATOM*> from this class.
  const ScreamAtomV& getAtomList() const;  ///< Returns atom_list. 
  ScreamAtomV getNewAtomList(); ///< Returns atom_list, but construct a new list of atoms. 
  SCREAM_ATOM* getAtom(int) const;	///< Returns pointer to atom numbered n.  Returns NULL if not found.
  SCREAM_ATOM* getAtom(MutInfo, string) const;  ///< Returns pointer to atom with MutInfo and atomlabel specified by string.  If multiple occurance of atom label, return "first" one (in atom list ordering).
  ScreamAtomV getTheseAtoms(vector<int>&) const; ///< Returns SCREAM_ATOM*'s with atom numbers corresponding to ones provided in list.

  /* Hydrogen adding functions, Connectivity adding functions. */

  void addHydrogens(); ///< 7-21-06: Adds hydrogens to AAChain (amino acid) only.  
  void addConnectivity(); ///< 7-21-06: Adds connectivities to AAChain (peptide) only.
  void assignFFType(); ///< adds FF Type using internal RTF object.
  void assignFFType(SCREAM_RTF*); ///< adds FF type using passed in SCREAM_RTF*.
  /* Binding site determination functions. */

  vector<MutInfo> residuesAroundAtomN(vector<int>, double, string) const; // double: distance.  string: SidechainOnly or BackBoneAsWell.
  vector<MutInfo> residuesAroundResidue(vector<MutInfo>, double, string) const;
  vector<MutInfo> residuesAroundChain(vector<string>, double, string) const;

  vector<MutInfo> residuesAroundAtom(ScreamAtomV&, double, string) const; // Calculates which residues are within distance to atoms specified on protein.

  /* Object State functions. */

  int totalComponents() const;	///< Returns total number of Components in this Protein object.
  int mutationDone() const; ///< Returns whether mutations were made to this protein from the last ntrlRotamerPlacement step.  

  void setMutInfoStrainEnergy(MutInfo, double); ///< Sets the energy assocaited with this mutInfo.
  double getMutInfoStrainEnergy(MutInfo); ///< Returns the strain energy associated with that MutInfo.

  void printAtomFlagStates(); ///< Prints atom moveable states.
  /* Insert Create CB parameters stuff here. */
  std::string getPlacementMethod() { return this->placementMethod;};
  void setPlacementMethod(std::string method) { this->placementMethod = method;};

  double getOffBisectorAngle() { return this->CreateCBParameters[0];};
  double getOffPlaneAngle() {return this->CreateCBParameters[1];};
  double getBondLength() { return this->CreateCBParameters[2];};
  double getRotamerMatchVectorLamdba() { return this->CreateCBParameters[3];};
  
  void setOffBisectorAngle(double d) { this->CreateCBParameters[0] = d;};
  void setOffPlaneAngle(double d) { this->CreateCBParameters[1] = d;};
  void setBondLength(double d) { this->CreateCBParameters[2] = d;};
  void setRotamerMatchVectorLamdba(double d) { this->CreateCBParameters[3] = d;};

  /* Sidechain maniuplation functions.  Meat! */

  int ntrlRotamerPlacement(string, int, AARotamer*); ///< Rotamer placement function.  Need to be careful: EE for all_interaction_calculation needs to be re-setup if mutation is done (returns 1 if there's mutation).  Otherwise, use the ntrlRotamerPlacement function in scream_EE, a short cut.  
  AARotamer* getAARotamer(std::string, int); ///< Returns an AARotamer on chain string and pstn int.

  int conformerPlacement(Rotamer* conformer, RotConnInfo* ); ///< replaces conformation give a Rotamer structure.  This is the more general version of SC_replacement, which is only used on a AARotamer*.  atomMapping gives the mapping between the atoms on rotamer and the atoms to be placed on the protein.  anchorAtoms gives the atom number (on rotamer) which one wishes to be overlayed on top.  sidechainAtoms stands for the atoms one actually wants to be placed; the "variable" part if you will.  For a AminoAcid this will the sidechain, for a random ligand this will simply be the entire ligand, since ligands in the library will have different conformers.  All these info are contained in RotConnInfo.
  Rotamer* conformerExtraction(RotConnInfo*); ///< returns a NEWed Rotamer* object. Remember to remove it from memory later on!

  int rotamerClusterPlacement(Rotamer* cluster, MutInfo* mI); ///< Placement of rotamers by passing in a RotamerCluster* and a tree-linked ClusterMutInfo* object.  Argument for based objects: more general.
  void setRotamerClusterEmptyLatticeEnergy(Rotamer* cluster, MutInfo* mI, double E);
  double getRotamerClusterEmptyLatticeEnergy(Rotamer* cluster, MutInfo* mI);

  int mutation(string, int, string); ///< (string Chain, int pstn, string Amino Acid, 3 letter name).  Does mutation on specfied Chain-Pstn to specified Amino Acid.  Returns whether or not successful.

  void setPreCalcEnergy(string, int, double); ///< Sets the PreCalc energy for a certain residue position.
  double getPreCalcEnergy(string, int) const; ///< Returns the PreCalc energy for that residues position.

  void setEmptyLatticeEnergy(string , int, double); ///< Sets the Empty Lattice energy for a certain residue position.
  double getEmptyLatticeEnergy(string , int) const; ///< Returns the Empty Lattice energy for a certain residue position.
  
  /* Functions that deal with setting resolution of atoms. */
  void setSideChainLibraryName(string , int, string); ///< Sets all atoms on sidechain to a certain library resolution.
  void setProteinLibraryName(string); ///< Sets all sidechain atoms on protein to a certain library resolution.


  /* Functions that act on the entire protein system.
     Mostly, these are functions that need the entire atom_list info.
   */

  void resetFlags(); ///< Resets fixed moveable flags.

  vector<int> getNewMapping() const { return this->new_mapping;}; ///< returns new mapping.

  double sc_clash(string, int) const;  ///< Returns distance of worst clash between sidechain and rest of protein.  1st arg: chain, 2nd arg: pstn.
  double sc_clash(string, int, ScreamAtomV&) const; ///< Returns distance of worst clash between sidechain and a list of atoms passed in.  1st arg: chain, 2nd arg: pstn, 3rd arg: list of atoms.
  
  double conformer_clash(RotConnInfo*) const;  ///< Checks clash between conformer and rest of protein system.
  double conformer_clash(RotConnInfo*, ScreamAtomV&) const; ///< Checks clash between conformer and list of atoms passed in.

  /* CRMS calculation functions. 
   */
  //  double sc_CRMS(string, int, Rotamer*) const;  ///< Calculates the CRMS distance between one position of the current protein to a specified Rotamer.
  double sc_CRMS(string, int, Protein*) const;  ///< Calculates the CRMS distance between one position of one Protein/system and current Protein system.

  double conformer_CRMS(Protein*, RotConnInfo*) const;

  pair<double, string> max_atom_dist_on_SC(string, int, Protein*) const; ///< Calculates the MAX distance between two equivalent atoms on two sidechains.
  double sc_atom_CRMS(int, Protein*) const; ///< Calculates the distance between the two equivalent atoms specified by int on the two Proteins.  Does a "flip" for residues like Phe, Asp, Glu, Tyr, Arg since they have equivalent atoms.


  /* Convenience helper functions. */
  void fix_entire_atom_list_ordering() { this->_fix_entire_atom_list_ordering();}; ///< Fixes the ordering of atoms according to standard bgf protein atom ordering format.

  /* MPSim functions.  Not sure if obsolete. */

  void fix_toggle(bool);                   ///< if TRUE, fixes all, if FALSE, all moveable.
  void fix_sc_toggle(bool, int, string);        ///<b makes a SC moveable.  string = CHAIN, int = pstn.

  /* Print functions. A pointer to an output stream is passed in in all of the following functions. */

  void append_to_filehandle(ostream*) const; ///< Prints to file; only prints atom bgf info, no connectivity, no header info, etc.
  void pdb_append_to_filehandle(ostream*) const; ///< Prints pdb to filehandle.

  void print_bgf_file(ostream*) const;  ///< The Print BGF File routine.

  /* Functions that are now obsolete.  Written for easy debugging purposes.*/

  void print_Me() const;           ///< Prints.  Obsolete.
  void print_ordered_by_n() const; ///< Prints.  Obsolete.

private:

  /* Private member varialbes.  Protein not meant to be inherited.
   */
  vector<ProteinComponent*> pc_v;                    ///< vector of ProteinComponents.
  map<MutInfo, double> mutInfoStrainEnergies; ///< Stores strain energies (basicaly PreCalc numbers) in Protein structure.

  ScreamAtomV *ptr_to_atom_v;	///< 11/10/04.  Memory management philosophy: all memories (except mutations) dealt at a high level.  Creation of atom_list not done in Protein class; either in sc_bgf_handler (bgf reader) or from ModSimSCREAMInterface (atom_list from elsewhere).  If class Protein never messes with the underlying SCREAM_ATOM*'s, all is fine.  But this is not the case because of mutations, which is only easily implementable in a class like Protein.  In other words, class Protein does mess with the underlying SCREAM_ATOM*'s.  This means there needs to be some way to tell higher level structures about these changes, otherwise when those higher level structures deallocates SCREAM_ATOM*'s, it would try to deallocate the SCREAM_ATOM*'s that have been mutated away (deleted by class Protein), leading to a seg fault.  The only solution I can think of at this point is to make a pointer ptr_to_atom_v to hold the external atom_list that is passed in to the constructor of the Protein class.  Memory management issues are complex, and I did not have the foresight to see this coming.  I have been using ScreamAtomV atom_v all along, and unwilling to change all instances of atom_v to *ptr_to_atom_v in my code, I am just going to make atom_v = (*ptr_to_atom_v).
  // Further note: this means I'll need to update this *ptr_to_atom_v value periodcially.  Alternatively, change all instances of atom_v to *ptr_to_atom_v.
  ScreamAtomV &atom_v;// = *ptr_to_atom_v;                       ///< vector of SCREAM_ATOM's.  see above.
  SCREAM_RTF* rtf; ///< the residue topology file structure.

  string placementMethod; ///< Method for placing ntrl rotamer sidechains.  Options: Default, CreateCB, and UseExistingCB.
  vector<double> CreateCBParameters; ///< Parameters for creating CB.
  int _CBCalc_flag; ///< if = 1, consider CB as part of protein backbone.  so that after placement updates of flag (protein fixed) value would be correct.

  vector<int> new_mapping; ///< If ntrlRotamerPlacement induces a mutation, this contains the new mapping from old to new atoms.  This is used for fixing RotConnInfo stuff (ArbLib).  Index == atom_n.  Note: when index == 0 entry is dummy, i.e. new_mapping[0] == 0, always.  The list starts at 1.  9-12-05.

  /* Private member functions. */

  void InitDataStructures(ScreamAtomV&);   ///< Initializes data structures from a vector of SCREAM_ATOM*'s.  Not a const because mutations will mess with underlying SCREAM_ATOM* structures.

  /* obsolete because of sc_bgf_handler */
  void make_bond_from_connect_info(const vector<string>&);  ///< makes bonds from bgf file style connectivity info.  

  /* Helper functions */

  void _fix_entire_atom_list_ordering(); ///< Fixes the ordering of atoms according to standard bgf protein atom ordering format.
  void _fix_charges(); 		///< Fixes charges on atoms using CHARMM2 charges.  temp scheme.
  void _fix_residue_in_atom_list_ordering(string, int); ///< Fixes the ordering atoms from only a portion of the atom list.


  /* Readibility improvement functions */
  vector<int> _getAtomNumbersInProtein(const map<int, int>&, const vector<int>&) const; ///< Returns a list of atoms with which 
  ScreamAtomV _mutationHelpers_alloc_new_sc_atoms(string, int, AARotamer*, ScreamAtomV&, map<SCREAM_ATOM*, SCREAM_ATOM*>&);
  void _mutationHelpers_connectivities_fix(string, int, AARotamer*, ScreamAtomV&, map<SCREAM_ATOM*, SCREAM_ATOM*>&);
  void _mutationHelpers_PRO_connectivities_and_numbering_fix(string, int, AARotamer*, ScreamAtomV&, map<SCREAM_ATOM*, SCREAM_ATOM*>&);
  void _mutationHelpers_GLY_connectivities_and_numbering_fix(string, int, AARotamer*, ScreamAtomV&, map<SCREAM_ATOM*, SCREAM_ATOM*>&);
  void _mutationHelpers_init_new_mapping(ScreamAtomV&);
  void _mutationHelpers_insert_new_sc_atoms(string, int, AARotamer*, ScreamAtomV&, ScreamAtomV&);
  void _mutationHelpers_delete_old_sc_atoms(ScreamAtomV&);
  void _mutationHelpers_renaming_mut_atoms(string, int, string, ScreamAtomV&);

  

};

#endif /* SC_PROTEIN_HPP */



