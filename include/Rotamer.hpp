//! Rotamer.hpp
/** 
 * Rotamer.hpp
 *
 * Header file for classes relevant to rotamer in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.  
 *
 */

#ifndef ROTAMER_HPP
#define ROTAMER_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"
#include "sc_BackBone.hpp"
#include "sc_SideChain.hpp"
#include "RotConnInfo.hpp"


/**
 * Rotamer Class.  Base class for other types of rotamer classes.
 */

class Rotamer {

public:
  
  // Temporary measure: introduce friends.  Moving get_bb() and get_sc() to protected part of the class so that Don't need to worry about SWIG.
  friend class Protein;
  friend class AminoAcid;
  friend class AAChain;

  Rotamer();
  Rotamer(const ScreamAtomV&, const RotConnInfo*, bool = true); ///< Initializes a Rotamer object with a list of atoms passed in.  No need to delete these atoms afterwards if called through this method.  3rd arg: if true, atom list is from Protein, meaning needs to use atom_n_map info in rotConnInfo.
  Rotamer(const stringV& rotamerCoords, const RotConnInfo* = NULL); ///< Initializes a rotamer object with strings of rotamer coordinate lines.  Bgf Format.
  //  Rotamer(const string, const string); ///< Constructor with 1st arg being rotamer structure file, 2nd arg being rotamer conenctivity file.
  virtual ~Rotamer();

  Rotamer& deepcopy(const Rotamer&); 	///< Deepcopies a rotamer object.

  //  void initConnectivity(const stringV& connectivityLines); ///< Initializes the connectivity between atoms in a rotamer.

  bool read_cnn_lines(const stringV); ///< reads a connectivity file.
 
  virtual void print_Me() const;  ///< Prints atom in bgf style, unordered.
  virtual void print_ordered_by_n() const;  ///< Prints atom in bgf style, ordered by atom number.
  virtual void append_to_filehandle(ostream*) const {};
  virtual void pdb_append_to_filehandle(ostream*) const {}; 
  virtual void append_to_ostream_connect_info(ostream*) const {};
  
  void printEnergies() const; ///< Prints energies, with Total, PreCalc, VDW, HB and Coulomb.
  
  // double kai_angle(int) {};             // returns n'th kai angle value
  //  virtual void match_onto_bb() = 0;    // should put this in AminoAcid?
  // virtual void read_library(string lib_name) = 0;  // 
  //  virtual Rotamer& operator=(const Rotamer&) = 0; ///< Copy assignment.
  //virtual Rotamer& copy(const Rotamer&) = 0; ///< Copy expression.

  double self_E;		       ///< The self energy of the sidechain.  Backbone emptied environment.  It's here (public) because I'm lazy and don't want to bother with all the get/set functions.
  int is_Original;                    ///< hack!!! 1: true, 0, false
  
  int get_is_Original_flag() const { return this->is_Original;};
  void set_is_Original_flag(int FLAG) {this->is_Original = FLAG; };
  
  
  bool same_backbone;                  ///< hack!!!  no matrix matching necessary.

  string library_name; ///< library name from this rotamer came from

  virtual BackBone* get_bb() const { return this->bb;};
  virtual SideChain* get_sc() const { return this->sc;};

  virtual ScreamAtomV get_sc_atoms() const { return this->sc->get_atoms(); };
  virtual ScreamAtomV get_bb_atoms() const { return this->bb->get_atoms(); };

  virtual SCREAM_ATOM* getAtom(int); 	///< Returns a single SCREAM_ATOM* by its number.
  virtual ScreamAtomV getTheseAtoms(vector<int>); // returns a vector of SCREAM_ATOM*'s with corresponding atom numbers to specified.

  virtual vector<Rotamer*> getAllRotamers() { vector<Rotamer*> rot_list; rot_list.clear(); rot_list.push_back(this); return rot_list;}; ///< Mainly for purposes of RotamerCluster getAllRotamers() function.

  void fix_toggle(bool);  ///< Toggles fix(true)/moveable(false) flag.  
  void fix_sc_toggle(bool);   ///< Toggles fix/moveable flag for sc.
  void fix_bb_toggle(bool);   ///< Toggles fix/moveable flag for bb.

  virtual int number_of_atoms() const;
  virtual double total_charge() const;
  int get_rotamer_n() const {return rotamer_n;};
  int get_mult_H_n() const {return mult_H_n;};
  string get_library_name() {return this->library_name;}; 

  /* Get/set functions for rotamer empty lattice energies and current base-statement-relative energies. */
  double get_empty_lattice_E() const {return empty_lattice_E;};

  void set_empty_lattice_E(double E) {empty_lattice_E = E;};

  double get_empty_lattice_E_abs() const {return empty_lattice_E_abs;};

  void set_empty_lattice_E_abs(double E ) {empty_lattice_E_abs = E;};

  int get_empty_lattice_energy_rank() { return empty_lattice_energy_rank;};
  void set_empty_lattice_energy_rank(int rank) { empty_lattice_energy_rank = rank;};

  void setFailedDistanceCutoff() { this->failedCutoff = 1;};
  void setPassedDistanceCutoff() { this->failedCutoff = 0;};
  int failedDistanceCutoff() { return failedCutoff; } ;

  int sameResidueTypeAs(Rotamer*); ///< Compares this rotamer type with input rotamer.  checks by looking at the first atom in each rotamer and compare the RES type.
  
  /* Updated Rotamer get/set functions 3-13-06. */
  std::string get_preCalc_Energy_Line() const;
  void populate_preCalc_Terms(std::string);

  double get_preCalc_TotE() const {return preCalc_TotE;};
  double get_preCalc_BondsE() const {return preCalc_BondsE;};
  double get_preCalc_AnglesE() const {return preCalc_AnglesE;};
  double get_preCalc_TorsionsE() const {return preCalc_TorsionsE;};
  double get_preCalc_InversionsE() const {return preCalc_InversionsE;};
  double get_preCalc_CoulombE() const {return preCalc_CoulombE;};
  double get_preCalc_vdwE() const{return preCalc_vdwE;} ;
  double get_preCalc_HBondE() const {return preCalc_HBondE;};
  double get_preCalc_SolvE() const {return preCalc_SolvE;} ;
 
  void set_preCalc_TotE(double E) {preCalc_TotE = E;};
  void set_preCalc_BondsE(double E) {preCalc_BondsE = E;};
  void set_preCalc_AnglesE(double E) {preCalc_AnglesE = E;};
  void set_preCalc_TorsionsE(double E) {preCalc_TorsionsE = E;};
  void set_preCalc_InversionsE(double E) {preCalc_InversionsE = E;};
  void set_preCalc_CoulombE(double E) {preCalc_CoulombE = E;};
  void set_preCalc_vdwE(double E) {preCalc_vdwE = E;};
  void set_preCalc_HBondE(double E) {preCalc_HBondE = E;};
  void set_preCalc_SolvE(double E) {preCalc_SolvE = E;};
 

  /* Rotamer energy get/set functions */

  double get_rotlib_E() const {return rotlib_E;};
  double get_sc_valence_E() const {return sc_valence_E;};
  double get_sc_coulomb_E() const {return sc_coulomb_E;};
  double get_sc_vdw_E() const {return sc_vdw_E;};
  double get_sc_hb_E() const {return sc_hb_E;};
  double get_sc_total_nb_E() const {return sc_total_nb_E;};
  double get_sc_solvation_E() const {return sc_solvation_E;};
  double get_sc_total_E() const {return sc_total_E;};

  void set_rotamer_n(int n) {rotamer_n = n;};
  void set_rotlib_E(double E) {rotlib_E = E;};
  void set_sc_valence_E(double E) {sc_valence_E = E;};
  void set_sc_coulomb_E(double E) {sc_coulomb_E = E;};
  void set_sc_vdw_E(double E) {sc_vdw_E = E;};
  void set_sc_hb_E(double E) {sc_hb_E = E;};
  void set_sc_total_nb_E(double E) {sc_total_nb_E = E;};
  void set_sc_solvation_E(double E)  {sc_solvation_E = E;};
  void set_sc_total_E(double E)  {sc_total_E =E;};

  /** Reassigns the sidechain position of current rotamer to that of passed in rotamer. 
   * match_bb takes in a second rotamer of the same type (e.g. an AARotamer takes in an AARotamer); otherwise bad things could happen.  
   */

  virtual void match_bb(const Rotamer*);
  virtual void match_CB(const Rotamer*, const SCREAM_ATOM* const, double );

  virtual void assign_atom_fftype() {};
  virtual void assign_charges(string scheme) {};
  virtual void assign_lone_pair() {};

  bool declaredInRotlibScope() {return this->declaredInRotlibScope_m;};	///< returns declaredInRotlibScope member.
  void setDeclaredInRotlibScope(bool b) {this->declaredInRotlibScope_m = b;}; ///< sets declaredInRotlibScope_m to a value.

protected:

  bool allocatedScreamAtoms;	///< false if rotamer is generated from an atomList (memory precontructed) or from an original file (memory allocated within Rotamer).
  bool declaredInRotlibScope_m; ///< true if rotamer is instantiated within a Rotlib class.  this takes care of add_rotamer operation in Rotlib.

  ScreamAtomV atom_list;	///< Rotamer atoms.

  BackBone* bb;            	///< Pointer to BackBone structure.
  SideChain* sc;		///< Pointer to SideChain structure.

  vector<int> RotamerAxisAtoms;	///< The two atoms on the rotamer axis.  Labeled by atom number.  Read in from a .cnn file.
  vector<int> AnchorAtoms;	///< Reference atoms to which one "fixes" the rotamer.  Like (N, C, CA) or (N, C, CB) in AminoAcids.

  string conformerName;		///< Name of this Confomer.
  
  int rotamer_n;                       ///< Rotamer number from Rotlib file (heavy atom rotamer)
  int mult_H_n;                        ///< Certain AminoAcids (Ser, Thr, Tyr) allow multiple H positions.  This parameter denotes which H rotamer.

  /* Rotamer empty lattice energies and zeroth base state reference energies, plus current base reference state energies */
  double empty_lattice_E;	///< Energy of this rotamer in its empty lattice environment.  Relative to reference state, i.e. always >= 0.
  
  double empty_lattice_E_abs;	///< Energy of this rotamer in its empty lattice environment.  Absolute value, i.e. not relative to its reference state.

  int empty_lattice_energy_rank; ///< the rank of this rotamer relative to other rotamers in the empty lattice environment.  0 is best, 1 is next best... so on.
  
  int failedCutoff;		///< whether this rotamer passed or failed a distance cutoff test with the rest of the protein.

  /* Rotamer energy variables */
  // 3-13-06 New set: easier on my mind to use the following new names.
  // Note: 3-13-06: precalc energy only includes the valence and VDW.  no Coulomb or HBond.  
  double preCalc_TotE;
  double preCalc_BondsE;
  double preCalc_AnglesE;
  double preCalc_TorsionsE;
  double preCalc_InversionsE;
  double preCalc_CoulombE;
  double preCalc_vdwE;
  double preCalc_HBondE;
  double preCalc_SolvE;       



  double rotlib_E;                     ///< Energy value from rotamer library.  
                                       ///< This energy value is the vacuum energy of this rotamer in an Ala-X-Ala tripeptide with ends capped by an N-Methyl and methylaldehyde.  While it is acceptable to use this value to initially screen out bad rotamers, consideration of the rotamer's environment is integral to sidechain placement and hence such an initial approximation is of dubious value.
                                       ///< The above is now obsolete.  This is actually the energy of just the sidechain in vacuum w/o electrostatic effects.
  double sc_valence_E;                 ///< Valence terms of rotamer sidechain energy.
                                       ///< This is sum of sidechain energy_per_atom_b.
  double sc_coulomb_E;                 ///< Coulombic energy with the rest of the protein. 
  double sc_hb_E;                      ///< HB energy.
  double sc_vdw_E;		       ///< VDW energy of sidechain with rest of protein.
  double sc_total_nb_E;		       ///< sc_coulomb_E + sc_vdw_E + sc_hb_E

  double sc_solvation_E;	       ///< Solvation for just the sidechain.
                                       ///< Both polar solvation and cavity solvation included.  If doing just FSM in MPSim, cavity term is zero so is automatically ignored.
  double sc_total_E;                   ///< Total Energy of rotamer in its (protein) environment.  Includes valence terms, hydrogen bonds, VDW, etc.   
  
protected:
  /* Helper functions */
  void _setDefaults(); 		///< Initializes some member variables to their default values.  Used mostly from Constructors.
  void _initAtomList(const stringV&); ///< initializes an atom list from a string of vectors.
  void _initSideChain(ScreamAtomV&, const vector<int>&); ///< list of atoms numbers are provided; these atom numbers are taken as sidechain atoms (i.e. variable atoms).
  void _initBackBone(ScreamAtomV&, const vector<int>&); ///< list of ints (atom numbers) provided; these are the backbone atom numbers, i.e. the anchor points, aka non-variable region.


private:


};
 
#endif /* ROTAMER_HPP */
