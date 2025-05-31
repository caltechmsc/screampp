/** v2.0: updates:
 *  removed calc_empty_lattice_E().  (no argument version)
 */


#ifndef SCREAM_HB_EE
#define SCREAM_HB_EE

#include "defs.hpp"
#include "MutInfo.hpp"
#include "scream_E_functionals_hb.hpp"
#include "sc_Protein.hpp"

#include "RotamerNeighborList.hpp"

class HB_EE {
  /* need to take care of short distances--okay for LJ H-bond for now */
  
  friend class VDW_HB_Exclusion_EE;

public:


  HB_EE();
  HB_EE(Protein*, vector<MutInfo>, SCREAM_HB_OBJ*);
  ~HB_EE();

  SCREAM_HB_OBJ* hb_obj;

   /* initialization routines for easy SWIGging. */
  void init_after_addedMutInfoRotConnInfo(Protein*, SCREAM_HB_OBJ*);
  void init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein*, SCREAM_HB_OBJ*);
  void init_after_addedMutInfoRotConnInfo_neighbor_list(Protein*, SCREAM_HB_OBJ*, RotamerNeighborList*);
  void addMutInfoRotConnInfo(MutInfo, RotConnInfo*);

  double calc_empty_lattice_E(const MutInfo);

  double calc_empty_lattice_E_delta(const MutInfo&, string , double = 0);


  double calc_residue_interaction_E(const MutInfo);
  double calc_residue_interaction_E(const MutInfo, const MutInfo, string = "FLAT", double = 0); // calculates inter-residue energy (only between sidechains)
  double calc_all_interaction_E(); // calculates the overall interaction energy between all variable atoms (i.e. sidechain interaction energies)
  double calc_all_interaction_E_delta(string , double = 0);
  int potential_HB_acceptor(SCREAM_ATOM*);

  /* Functions for pecific components of energy. */
  double calc_EL_rot_selfBB(const MutInfo&, std::string, double);
  double calc_EL_rot_otherBB(const MutInfo&, std::string, double);
  double calc_EL_rot_fixedSC(const MutInfo&, std::string, double);
  double calc_EL_rot_fixedHET(const MutInfo&, std::string, double);
  double calc_EL_rot_moveableHET(const MutInfo&, std::string, double);

  /* Setup 50% on the fly loop.  by 50% i mean it's still on the fly, but not like a permanent neighborlist. */
  void setup_variableAtomsOnEachSidechain(); ///< Sets up this->_variable_atoms_on_each_sidechain.  Call this everytime the energy expression has been updated, as in the fixed and moveable atoms have been modified.


private:

  int ON_THE_FLY;
  Protein *ptn;
  RotamerNeighborList* rotamerNeighborList;

  ScreamAtomV fixed_atoms;
  ScreamAtomV variable_atoms;

  vector<AASideChain*> sidechains_2b_SCREAMed;

  /* Private variable that saves time for the on_the_fly_loops. */
  map<int, const MutInfo*> _mutInfo_n_map;
  map<int, ScreamAtomV> _variable_atoms_on_each_sidechain;

  /* MutInfo storage */
  map<MutInfo, RotConnInfo*> mutInfo_rotConnInfo_map; ///< Stores info about which residues to be SCREAMed.  Necessary for initialization.

  /* Atom triples; for H-bonding */
  map<const MutInfo, vector< vector<SCREAM_ATOM*> > > variable_and_fixed; // interaction between variable and fixed atoms. inner vector<SCREAM_ATOM*> has size 3: contains A*, H*, D*.
  void _initVariableAndFixedAtomTripleList(Protein*, const vector<MutInfo>); // init the above structure
  void _initVariableAndFixedAtomTripleListArb(Protein*, const map<MutInfo, RotConnInfo*>); ///< init the above structure, except use mutInfo_rotConnInfo_map
  //vector< vector<SCREAM_ATOM*> > variable_and_variable; // interactions between variable atoms, i.e. sidechains.  sometimes ligands.

  void _initVariableAndVariableAtomTripleList(Protein*, const vector<MutInfo>);
  void _initVariableAndVariableAtomTripleListArb(Protein*, const map<MutInfo, RotConnInfo*>); ///< init the above structure, except use mutInfo_rotConnInfo_map

  double _calc_empty_lattice_E_on_the_fly_loop(const MutInfo&, string, double); ///< loops over on the fly energy calculations.
  double _calc_all_interaction_E_on_the_fly_loop(string, double); ///< loops over overall energy on the fly calculations.
  
  /* Other help functions. */
  void _figureOutHBPotentialMethods(const string, string&, int&, int&, double, SCREAM_HB_BASE_FUNCTIONAL_OBJ**); ///< helper function.

  map<MutInfo, ScreamAtomV> each_sc_atom_list; ///< stores rotamer atoms by a residue basis.  points to atoms in ptn*. 
  map< MutInfoPair, vector< vector<SCREAM_ATOM*> > > variable_and_variable;
  

};







#endif /* SCREAM_HB_EE */
