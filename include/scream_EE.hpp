/** SCREAM_EE class.
 *  Stands for SCREAM Energy Engine.
 */


#ifndef SCREAM_EE 
#define SCREAM_EE

#include "defs.hpp"
#include "MutInfo.hpp"

#include "sc_Protein.hpp"
#include "RotConnInfo.hpp"

#include "scream_coulomb_EE.hpp"
#include "scream_hb_EE.hpp"
#include "scream_vdw_EE.hpp"
#include "scream_vdw_hb_exclusion_EE.hpp"

#include "ClashCollection.hpp"
#include "RotamerNeighborList.hpp"

#include "scream_ctl_reader.hpp"

class Scream_EE {
public: 
  Scream_EE();
  //  Scream_EE(Protein*, vector<MutInfo>, string);
  //Scream_EE(Protein*, vector<MutInfo>, string, ClashCollection* = NULL);
  ~Scream_EE();
  
  // Note on ClashCollection* object: this class needs to comuunicate with RotlibCollection.  Putting a member that stores ClashCollection is A solution, though not a particularly elegant one.  

  /* Initialization routines; for ease of SWIGing */

  void init(Protein*, vector<std::string>, string, string); // Initializes.  vector<string>: list of MutInfo string.  1st string: FF_file.  2nd string: delta_file.
  void setCalcNonPolarHydrogen_flag(int flag) { this->_calcNonPolarHydrogen_flag=flag; } ; ///< Sets this->_calcNonPolarHydrogen_flag.  Must sure set this before the init functions!  else irrelevant.  1 = yes, calc non polar H.  0 = no, don't cal non polar H, also turns electrostatics off.
  int getCalcNonPolarHydrogen_flag() {return this->_calcNonPolarHydrogen_flag; }; ///< Returns this->_calcNonPolarHydrogen_flag.  Must sure set this before the init functions!  else irrelevant.  1 = yes, calc non polar H.  0 = no, don't cal non polar H, also turns electrostatics off.  Though user has control in Python high level control whether or not to include electrostatics calculation.

  void addMutInfoRotConnInfo(MutInfo, RotConnInfo* = NULL); // adds a MutInfo-RotConnInfo* pair to object.
  //  void init_after_addedMutInfoRotConnInfo(Protein*, std::string, std::string); // string1: FF_file.  string2: SCREAM_delta_file.
  //  void init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein*, std::string, std::string); 
  //  void init_after_addedMutInfoRotConnInfo_neighbor_list(Protein*, std::string, std::string); ///< after adding MutInfo RotConnInfo pairs, call this if want to do neighbor list.  
  void init_after_addedMutInfoRotConnInfo(Protein*, ScreamParameters*); 
  void init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein*, ScreamParameters*);
  void init_after_addedMutInfoRotConnInfo_neighbor_list(Protein*, ScreamParameters*);

  /* Flag manipulation stuff.*/

  void fix_mutInfo(MutInfo&, RotConnInfo* =NULL, int=0); ///< Make atoms on this mutinfo fixed.  last argument: whether or not to re-setup variable_atoms_on_each_sidechain for all_interactions_E calculations.
  void moveable_mutInfo(MutInfo&, RotConnInfo* = NULL, int=0); ///< Make atoms on this mutinfo moveable.  last argument: whether or not to re-setup variable_atoms_on_each_sidechain for all_interactions_E calculations.
  void fix_all(); ///< Fix all atoms in system.

  void visible_mutInfo(MutInfo&, RotConnInfo* = NULL, int=0);  ///< Make atoms on this mutinfo visible, sc-sc interaction calculations.
  void invisible_mutInfo(MutInfo&, RotConnInfo* = NULL, int=0);

  void visible_EL_mutInfo(MutInfo&, RotConnInfo* = NULL, int=0); ///< make atoms on this mutInfo visible, EL calculations.
  void invisible_EL_mutInfo(MutInfo&, RotConnInfo* = NULL, int=0); 
  
  void visible_all(); ///< Makes all atoms visible, in terms of sc-sc interations.  Also resets all bb sidechain atom visibility, hence, not very useful.
  void invisible_all(); ///< Makes all atom invisible, n terms of sc-sc interations.

  void visible_EL_all(); ///< Make all atoms EL visible, also resets EL visibility of sidechain atoms.  Hence, not very useful.
  void invisible_EL_all(); ///< Makes all atoms EL invisible.

  void make_chain_invisible(string); ///< Makes all atoms on specified chain invisible (in terms of sc-sc energies).  No undo functions; use resetFlags() to undo.
  void make_chain_EL_invisible(string); ///< Makes all atoms on specified chain EL invisible.  No undo functions; use resetFlags to undo.

  void resetFlags(int=0); ///< Reset atom moveable fixed flags to original.  argument: whether or not to re-setup variable_atoms_on_each_sidechain for all_interactions_E calculations.

  void setup_variableAtomsOnEachSidechain(); ///< Sets up pre calcalulated list in scream_vdw_EE, scream_coulomb_EE for faster all_interaction calculation.  Note: this sets up the variable sidechain according to the very original fixed/moveable flag.

  /* init function for delta fields, do this after init_after_addedMutInfoRotConnInfo_neighbor_list etc. */
  void initScreamAtomDeltaValue(string, string, double, string); ///< initializes mu value for atom.  first string: library name.  second string: method, i.e. "FLAT" or "FULL". last double: either the mu value itself in "FLAT" case, or the alpha value in "FULL" case, i.e. how much of sigma (stdev) to reduce.  Last string: file for EachAtomDeltaFile.
  void initScreamAtomVdwHbFields() { this->_initScreamAtomVdwHbFields(this->_calcNonPolarHydrogen_flag);}; // wrapper for private function that apparently's not as private as i though.

  void addClashCollection(ClashCollection*); // Stores a ClashCollection object in this class.  Needs to be pointer.
  void cleanClashCollection(); // Cleans out informations about ClashCollection.

  /* Misc. routines. */

  double getDistanceDielectricPrefactor();
  void setDistanceDielectricPrefactor(double);
  void setNormalDielectric(double);

  int ntrlRotamerPlacement(string, int, AARotamer*); ///< Wrapper for ptn->ntrlRotamerPlacement that also takes care of EE expression for calc_all_interaction_E.  Use this wrapper when calc_all_interaction_E function will be immediately employed, but not necessary when doing singles (empty lattice) calculations.
  //int conformerPlacement(std::string, int) {}; ///< Currently not implemented, since no need to modify EE expression for calc_all_itneraction_E.

  /* Scream Energy Calculation Objects. */

  SCREAM_Coulomb_OBJ coulomb_obj;
  SCREAM_VDW_OBJ vdw_obj;
  SCREAM_HB_OBJ hb_obj;

  Coulomb_EE* coulomb_EE;
  VDW_EE* vdw_EE;
  HB_EE* hb_EE;
  VDW_HB_Exclusion_EE* vdw_hb_exclusion_EE;

  /* Empty lattice Energy Functions, meant for easy SWIG access. */
  
  double calc_empty_lattice_E(const MutInfo&);
  double calc_empty_lattice_E_full_delta(const MutInfo&, double); // double is how much sigma to add
  double calc_empty_lattice_E_flat_delta(const MutInfo&, double); // double is value of delta
  double calc_empty_lattice_E_scaled_inner_wall(const MutInfo&, double); // double is value of the scaling factor.  no greater than 1 please.

  //double calc_empty_lattice_E_full_delta_asym(const MutInfo&, double); // only asymmetric for VDW, not HB.  deprecated
  //double calc_empty_lattice_E_flat_delta_asym(const MutInfo&, double); // deprecated


  /* Boring wrapper functions not meant for swig. */
  double calc_empty_lattice_coulomb_E_delta(const MutInfo&);
  double calc_empty_lattice_vdw_E_delta(const MutInfo&, std::string, double);  // string: mode.  double: n-sigma/flat delta value/scaling factor.
  double calc_empty_lattice_hb_E_delta(const MutInfo&, std::string, double); // similar to above.
  double calc_empty_lattice_vdw_hb_exclusion_E_delta(const MutInfo&, std::string, double);

  /* Functions for: getting specific energy terms, all wrappers. 09-22-06. */
  double calc_EL_vdw_rot_selfBB(const MutInfo&, std::string, double);
  double calc_EL_vdw_rot_otherBB(const MutInfo&, std::string, double);
  double calc_EL_vdw_rot_fixedSC(const MutInfo&, std::string, double);
  double calc_EL_vdw_rot_fixedHET(const MutInfo&, std::string, double);
  double calc_EL_vdw_rot_moveableHET(const MutInfo&, std::string, double);

  double calc_EL_coulomb_rot_selfBB(const MutInfo&);
  double calc_EL_coulomb_rot_otherBB(const MutInfo&);
  double calc_EL_coulomb_rot_fixedSC(const MutInfo&);
  double calc_EL_coulomb_rot_fixedHET(const MutInfo&);
  double calc_EL_coulomb_rot_moveableHET(const MutInfo&);
  
  double calc_EL_hb_rot_selfBB(const MutInfo&, std::string, double);
  double calc_EL_hb_rot_otherBB(const MutInfo&, std::string, double);
  double calc_EL_hb_rot_fixedSC(const MutInfo&, std::string, double);
  double calc_EL_hb_rot_fixedHET(const MutInfo&, std::string, double);
  double calc_EL_hb_rot_moveableHET(const MutInfo&, std::string, double);



  /* Easy access of energy evaluation for swig's sake. */

  double calc_all_interaction_E();
  double calc_all_interaction_E_full_delta(double);
  double calc_all_interaction_E_flat_delta(double);
  
  /* Boring wrapper functions not meant for swig. */
  
  double calc_all_interaction_coulomb_E_delta();
  double calc_all_interaction_vdw_E_delta(std::string, double);
  double calc_all_interaction_hb_E_delta(std::string, double);
  double calc_all_interaction_vdw_hb_exclusion_E_delta(std::string, double);

  /* Calculation of interaction energy */

  double calc_residue_interaction_E(const MutInfo);
  double calc_residue_interaction_E(const MutInfo, const MutInfo);

  double calc_residue_interaction_vdw_E(const MutInfo, const MutInfo, string, double);
  double calc_residue_interaction_hb_E(const MutInfo, const MutInfo, string, double);
  double calc_residue_interaction_coulumb_E(const MutInfo, const MutInfo);

  /* Calculation of energies for various cases. */
  // WORKING ON THIS!  WHAT DO I NEED?
  // ABSOLUTE NEEDS:
  // 1. selected sidechains--selected sidechains E
  // 2. selected sidechains--selected backbone E
  
  // NICE TO HAVE:
  // 1. list of H bonds, and contributing energy terms.
  // 2. energy of entire structure

  /* Calculation of energies on the fly 
   * Separated out from the rest for more efficient calculation purposes.
   */
  
  void setProtein(Protein* ptn) { this->ptn = ptn;};
  Protein* getProtein() { return this->ptn;};

  Protein *ptn;
  vector<MutInfo> mutInfoV;
 
private:
  int _CBCalc_flag; ///< if = 1, consider CB as part of protein backbone.
  int _calcNonPolarHydrogen_flag; ///< if = 1, calculates non polar hydrogen VDW.  if = 0, doesn't calculate non polar hydrogen VDW, plus no electrostatics at all.
  map<MutInfo, RotConnInfo*> mutInfo_rotConnInfo_map; 
  //map<MutInfo, ScreamAtomV > residue_neighbor_list_map;
  RotamerNeighborList* rotamerNeighborList; ///< stores information about neighborlist for each rotamer.

  ClashCollection* clashCollection;

  void _read_FF_param_file(string);
  void _read_SCREAM_delta_param_file(string);
  void _read_EachAtomDeltaFile(std::string file, map<int, double>& deltaMap);

  /* init function for setting neighbor list */
  //void _initTolerantNeighborList(Protein*, map<MutInfo, ScreamAtomV >);

  /* init function for adding vdw_r, vdw_s and vdw_d to SCREAM_ATOM */
  void _initScreamAtomVdwHbFields(int); ///< initialized VDW fields in SCREAM_ATOM* construct so far during run don't have to look them up by string search/map.  int: this->_calcNonPolarHydrogen_flag value.


  /* init function for on the fly E calculations */
  void _initFixedMoveableAtomsOnProtein(Protein*, map<MutInfo, RotConnInfo*>);

  /* helper functions */
  double _calc_empty_lattice_E_delta(const MutInfo&, string, double = 0); // called by the front end calc_empty_lattice_E_flat_delta series of functions.  the double = 0 value is either: # of standard d away from mean or the flat delta value to use.  string: "FULL", "FLAT" or "RESIDUE".
  //double _calc_empty_lattice_E_delta_asym(const MutInfo&, string, double = 0); // deprecated
  double _calc_all_interaction_E_delta(string, double = 0); // similar as above.  string: "FULL", "FLAT" or "RESIDUE".

};

#endif /* SCREAM_EE */
