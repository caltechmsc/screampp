/** v2.0 updates:
 * removed: calc_empty_lattice_E().  (the no argument version)
 */

#ifndef SCREAM_VDW_EE
#define SCREAM_VDW_EE

#include "scream_E_functionals_vdw.hpp"
#include "sc_Protein.hpp"
#include "MutInfo.hpp"
#include "RotConnInfo.hpp"

#include "ClashCollection.hpp"
#include "RotamerNeighborList.hpp"

//class HB_VDW_Exclusion_EE;

class VDW_EE { // optiized
  /* excluded VDW interactions: 
   * 1. intra-sidechain VDW contacts
   * 2. backbone-backbone contacts
   * 3. 
   */
  friend class VDW_HB_Exclusion_EE; 

public:
  VDW_EE();
  VDW_EE(Protein*, vector<MutInfo>, SCREAM_VDW_OBJ*);
  ~VDW_EE();

  /* initialization routines for easy SWIGging. */
  void init_after_addedMutInfoRotConnInfo(Protein*, SCREAM_VDW_OBJ*);
  void init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein*, SCREAM_VDW_OBJ*);
  void init_after_addedMutInfoRotConnInfo_neighbor_list(Protein*, SCREAM_VDW_OBJ*, RotamerNeighborList*);
  void addMutInfoRotConnInfo(MutInfo, RotConnInfo*);
  void addClashCollection(ClashCollection* cc); 
  void cleanClashCollection();

  SCREAM_VDW_OBJ* vdw_obj;

  double calc_empty_lattice_E(const MutInfo&);
  double calc_empty_lattice_E_delta(const MutInfo&, string, double); // string: mode (FULL, FLAT, RESIDUE), double: for FULL, number of sigmas; for FLAT, value of delta, in angstroms.
  //  double calc_empty_lattice_E_delta_asym(const MutInfo&, string , double); // deprecated

  double calc_residue_interaction_E(const MutInfo&);
  double calc_residue_interaction_E(const MutInfo&, const MutInfo&, string = "FLAT", double = 0); // calculates inter-residue energy (only between sidechains)
  double calc_all_interaction_E(); // calculates the overall interaction energy between all variable atoms (i.e. sidechain interaction energies)

  double calc_all_interaction_E_delta(std::string, double);

  //double calc_residue_E(const MutInfo&); // calculates residue energy (only sidechain)
  map<MutInfo, RotConnInfo*> returnMutInfoRotConnInfoMap() { return this->mutInfo_rotConnInfo_map; };  //for SWIG ease


  /* Functions for pecific components of energy. */
  double calc_EL_rot_selfBB(const MutInfo&, std::string, double);
  double calc_EL_rot_otherBB(const MutInfo&, std::string, double);
  double calc_EL_rot_fixedSC(const MutInfo&, std::string, double);
  double calc_EL_rot_fixedHET(const MutInfo&, std::string, double);
  double calc_EL_rot_moveableHET(const MutInfo&, std::string, double);


  /* Setup 50% on the fly loop.  by 50% i mean it's still on the fly, but not like a permanent neighborlist. */
  void setup_variableAtomsOnEachSidechain(); ///< Sets up this->_variable_atoms_on_each_sidechain.  Call this everytime the energy expression has been updated, as in the fixed and moveable atoms have been modified.


  
private:
  int ON_THE_FLY; ///< if != 0, sets up lists.  if == 0, does not set up lists.
  Protein* ptn; ///< Needed to do energy calculations on the fly.
  ClashCollection* clashCollection; ///< Information about whether rotamers clash.
  RotamerNeighborList* rotamerNeighborList; ///< information about neighborlist for each rotamer.

  ScreamAtomV fixed_atoms;	///< includes the protein backbone, for sure.  ligand in some cases.  basically, atom coordinates that won't be modified by SCREAM.
  ScreamAtomV variable_atoms; 	///< includes: sidechains under question, perhaps all other sidechains.  ligand, maybe.  i.e. atoms whose coordinates will be modifed by SCREAM.
  
  vector<AASideChain*> sidechains_2b_SCREAMed; ///< sidechains that will be screamed; specified by user, in control file

  /* Private variables that save time for the on_the_fly_loops. */
  map<int, const MutInfo*> _mutInfo_n_map;  // lookup table for mutInfo_n to MutInfo.  integer operations faster than MutInfo operatrions.
  map<int, ScreamAtomV> _variable_atoms_on_each_sidechain; // // Stores ScreamAtomV for each MutInfo, using MutInfo index, i.e. mutInfo_n.  Used in calc_all_interactions stuff.


  /* Atom Pair info */

  // Need to save energy values in the Rotamer structure--remember this
  map<MutInfo, RotConnInfo*> mutInfo_rotConnInfo_map; ///< Stores info about which residues to be SCREAMed.  Necessary for initialization.
  map<const MutInfo, vector<AtomPair> > variable_and_fixed; ///< to get single energy excitation spectrum  
  void _initVariableAndFixedAtomPairList(Protein*, const vector<MutInfo>); ///< init the above structure
  void _initVariableAndFixedAtomPairListArb(Protein*, const map<MutInfo, RotConnInfo*>); ///< init the above structure, except use mutInfo_rotConnInfo_map

  // this is done on the fly--nothing is saved?  should i save some energy values?
  //vector<AtomPair> variable_and_variable; ///< add this energy to the single energy; save computational effort--
  void _initVariableAndVariableAtomPairList(Protein*, const vector<MutInfo>); ///< init the above structure
  void _initVariableAndVariableAtomPairListArb(Protein*, const map<MutInfo, RotConnInfo*>); ///< init the above structure, except use mutInfo_rotConnInfo_map

  double _calc_empty_lattice_E_on_the_fly_loop(const MutInfo&, string, double); ///< a loop that does on the fly energy calculations, takes in a pointer to an energy evaluation.
  double _calc_all_interaction_E_on_the_fly_loop(string, double); ///< a loop that does on the fly energy calculations for interactions between sidechains being replaced.

  /* Other help functions. */
  void _figureOutvdwPotentialMethods(const string, string&, int&, int&, double, SCREAM_VDW_BASE_FUNCTIONAL_OBJ**); ///< helper function.
  void _returnNoScaleFunctor(const string&, SCREAM_VDW_BASE_FUNCTIONAL_OBJ**); // Return an exact LJ function.  string is method string, contains a substring with 12-6, 11-6, 10-6, 9-6, 8-6, 7-6.
  void _setupAtomListStuff(const MutInfo&, ScreamAtomV&, ScreamAtomV&, ScreamAtomV&);
  void _remakeForNonPolarHCalc(int, ScreamAtomV&, ScreamAtomV&); ///< helper function 
  /* initialize fixed/moveable atom tags in scream_atom*/
  //void _initFixedMoveableAtomsOnProtein(Protein*, const map<MutInfo, RotConnInfo*>);

  map<MutInfo, ScreamAtomV> each_sc_atom_list; ///< stores rotamer atoms by a residue basis.  points to atoms in ptn*. 
  map< MutInfoPair, vector<AtomPair> > variable_and_variable;

};


#endif /* SCREAM_VDW_EE */
