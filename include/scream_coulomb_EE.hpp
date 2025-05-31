/** v2.0 updates
 *  removed calc_empty_lattice_E() (no argument version)
 */

#ifndef SCREAM_COULOMB_EE
#define SCREAM_COULOMB_EE

#include "defs.hpp"
#include "MutInfo.hpp"
#include "scream_E_functionals_coulomb.hpp"
#include "sc_Protein.hpp"

#include "RotamerNeighborList.hpp"

class Coulomb_EE {
  /* need to take care of short distances--okay for LJ H-bond for now */

public:
  Coulomb_EE();
  Coulomb_EE(Protein*, vector<MutInfo>, SCREAM_Coulomb_OBJ*);
  ~Coulomb_EE();

  SCREAM_Coulomb_OBJ* coulomb_obj;

  /* initialization routines for easy SWIGging. */
  void init_after_addedMutInfoRotConnInfo(Protein*, SCREAM_Coulomb_OBJ*);
  void init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein*, SCREAM_Coulomb_OBJ*);
  void init_after_addedMutInfoRotConnInfo_neighbor_list(Protein*, SCREAM_Coulomb_OBJ*, RotamerNeighborList*);
  void addMutInfoRotConnInfo(MutInfo, RotConnInfo*);
  
  double calc_empty_lattice_E(const MutInfo);  

  double calc_residue_interaction_E(const MutInfo);
  double calc_residue_interaction_E(const MutInfo, const MutInfo); // calculates inter-residue energy (only between sidechains)
  double calc_all_interaction_E(); // calculates the overall interaction energy between all variable atoms (i.e. sidechain interaction energies)
  double calc_all_interaction_E_delta(); // same thing as above for coulombic terms, for now.

  /* Functions for specific components of energy. */
  double calc_EL_rot_selfBB(const MutInfo&);
  double calc_EL_rot_otherBB(const MutInfo&);
  double calc_EL_rot_fixedSC(const MutInfo&);
  double calc_EL_rot_fixedHET(const MutInfo&);
  double calc_EL_rot_moveableHET(const MutInfo&);

  /* Setup 50% on the fly loop.  by 50% i mean it's still on the fly, but not like a permanent neighborlist. */
  void setup_variableAtomsOnEachSidechain(); ///< Sets up this->_variable_atoms_on_each_sidechain.  Call this everytime the energy expression has been updated, as in the fixed and moveable atoms have been modified.

private:

  int ON_THE_FLY; 
  Protein *ptn;
  RotamerNeighborList* rotamerNeighborList;

  ScreamAtomV fixed_atoms;
  ScreamAtomV variable_atoms;

  /* Private variables that save time for the on_the_fly_loops. */
  map<int, ScreamAtomV> _variable_atoms_on_each_sidechain; // Used in _calc_all_interaction_E_on_the_fly_loop(), a 50% neighborlist build up.  this is setup by calling 

  vector<AASideChain*> sidechains_2b_SCREAMed;

  /* MutInfo storage*/
  map<MutInfo, RotConnInfo*> mutInfo_rotConnInfo_map; ///< Stores info about which residues to be SCREAMed.  Necessary for initialization.

  map<const MutInfo, vector<AtomPair> > variable_and_fixed; // interaction between variable and fixed atoms. inner vector<SCREAM_ATOM*> has size 3: contains A*, H*, D*.
  void _initVariableAndFixedAtomPairList(Protein*, const vector<MutInfo>); // init the above structure
  void _initVariableAndFixedAtomPairListArb(Protein*, const map<MutInfo, RotConnInfo*>); ///< init the above structure, except use mutInfo_rotConnInfo_map

  //  vector<AtomPair> variable_and_variable; // interactions between variable atoms, i.e. sidechains.  sometimes ligands.
  void _initVariableAndVariableAtomPairList(Protein*, const vector<MutInfo>);
  void _initVariableAndVariableAtomPairListArb(Protein*, const map<MutInfo, RotConnInfo*>); ///< init the above structure, except use mutInfo_rotConnInfo_map

  double _calc_empty_lattice_E_on_the_fly_loop(const MutInfo&);
  double _calc_all_interaction_E_on_the_fly_loop();

  double _linearSpline(double, double, double); // 1st double: spline_start^2. 2nd double: (spline_end - spline_start)^2, 3rd double, dist^2.

  map<MutInfo, ScreamAtomV> each_sc_atom_list; ///< stores rotamer atoms by a residue basis.  points to atoms in ptn*. 
  map< MutInfoPair, vector<AtomPair> > variable_and_variable;

};

#endif /* SCREAM_COULOMB_EE */
