/* sc_AminoAcid.hpp
 *
 * Header file for classes relevant to AminoAcids in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

/* Update Notes:
 * 2.0: 
 * removed xxx_on_free_store.  AminoAcid will never itself create atoms; this class will only be constructed through passing of atoms in.
 * removed aa_state.  Too convoluated, plus, now reading in PDB files, has to be made more general, i.e. more robust; can't detect CTERM NTERM merely by counting HT hydrogens or OXT labels.  This means the "assign_charges" series of functions need to updated.
 */


#ifndef SC_AMINO_ACID_HPP
#define SC_AMINO_ACID_HPP

#include "defs.hpp"

#include "sc_BackBone.hpp"
#include "sc_AABackBone.hpp"
#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"
#include "Rotamer.hpp"
#include "AARotamer.hpp"

#include <map>
#include <vector>
using namespace std;
/** AminoAcid.  
 *
 */
class AminoAcid {

public:

  /* Constructors and destructors. */

  AminoAcid();
  AminoAcid(const ScreamAtomV&);
  AminoAcid(int aa_state);
  AminoAcid(const AminoAcid&);
  AminoAcid(const AABackBone&, const AASideChain&);
  ~AminoAcid();

  /* A set of constants relevant to this class. */

  const static int MAINCHAIN = 0;
  const static int NTERM = 1;
  const static int CTERM= 2;
  const static int SINGLE = 3;


  /** The assignment operator does not always perform a deep copy--it copies everything up to the pointers to SCREAM_ATOMs.
   * Copies everything up to what the SCREAM_ATOMS pointers point to.  SCREAM_ATOM data are not deep copied if those SCREAM_ATOM's are not created on the heap.
   */

  AminoAcid & operator=(const AminoAcid &);

  /** The deepcopy routine performs a deep copy.  It copies everything including the SCREAM_ATOM structures.
   *  As the program runs, there will be two kind of SCREAM_ATOM pointers, memory-wise.  One kind of SCREAM_ATOMs points to a globally or externally stored SCREAM_ATOM list.  This happens, for instance, when a list of MPSim Atoms is passed in and this entire list is converted into a list of SCREAM_ATOM's, which would be global.  Subsequent data structures (AminoAcid, AASideChain etc.) would contain SCREAM_ATOM pointers that point to this external/global list of SCREAM_ATOMs.  The other kind of SCREAM_ATOM has its memory dynamically allocated.  For instance, a Rotlib object contains many Rotamers, each of which allocates memory for SCREAM_ATOMs as the rotamer library is being read in.  operator=() should be used only when pointer copying (shallow copy) is needed.  deepcopy() guarantees that the object is copied up to the data which the SCREAM_ATOM pointers point to.
   */

  AminoAcid & deepcopy(const AminoAcid &);

  /** operator[] returns pointer to SCREAM_ATOM.
   * The order of atoms is same as order of atoms in bgf files.  I.e. the backbone atoms first, then the sidechain atoms.  So operator[0] returns a pointer the SCREAM_ATOM labeled N.
   */
  SCREAM_ATOM* operator[](int) const;  // not implemented yet

  /** SC_replacement takes in a pointer to a rotamer.  
   * calls mathc_bb() in AARotamer
   */
  void SC_replacement(const AARotamer* const, string, vector<double>&);

  /** Replaces the sidechain of this AminoAcid with the rotamer passed in.
   * Calls function in AARotamer.
   */

  void translate(ScreamVector);
  void transform(ScreamMatrix);

  /** Returns an array of pointers to atoms belonging to this Component (chain, ligand, etc).
   *  For convenience's sake.  Good coding would not need this.
   */
  
  ScreamAtomV getAtomList() const;


  /* Toggle flags.  Written for MPSim-SCREAM integration. */

  void fix_toggle(bool); ///< Toggles fix(true)/moveable(false) flag.  
  void fix_sc_toggle(bool);   ///< Toggles fix/moveable flag for sc.
  void fix_bb_toggle(bool);   ///< Toggles fix/moveable flag for bb.


  string residue_type() const { return rot->get_resName();}; ///< Returns the type of this residue.

  //??? why is there a duplicate SC_replacement function?
  //  void SC_replacement(const AARotamer&); ///< Replaces sidechain.  See functions (match_bb) in AARotamer for additional info.

  void assign_atom_fftype(); ///< Assigns atom FF from atom label.  
                             ///< Anomolous atoms on Cterm and Nterm are not assigned by calling this function.  Those must be manually entered.
  void assign_charges(string, int);  ///< Assigns charges to atoms.  Argument can be "CHARM22" etc.  Second argument is AA_STATE value.
  void assign_lone_pair();           ///< Assigns lone pair and connectivity info to atoms.

  /* Angle calculation functions. */

  double PHI() const;			///< Calculates PHI angle.  
                                ///< C(i-1)-N-CA-C(i) dihedral, zero when C(i-1) and C(i) cis.  Right-handed rotation around N-CA bond.  Range: -180 to 180.
  double PHI(const AminoAcid* const) const; ///< Calculates PHI angle given next AminoAcid on chain..
  ScreamMatrix set_PHI(double, AminoAcid* = NULL);	///< Sets PHI angle.
                                ///< Rotates all atoms off the C end beginning C(i-1).  Returns Matrix that does the rotation.
  
  double PSI() const;			///< Calculates PSI angle.
                                ///< N(i)-CA-C-N(i+1) dihedral, zero when N(i) and N(i+1) cis.  Right-handed rotation around CA-Cbond.  Rnage: -180 to 180.
  double PSI(const AminoAcid* const) const; ///< Calculates PSI angle given previous AminoAcid on chain.
  ScreamMatrix set_PSI(double, AminoAcid* = NULL);	///< Sets PSI angle.
                                ///< Rotates all atoms off the C end beginning N(i+1).

  double OMEGA();		///< Calculates OMEGA angle.
                                ///< The OMEGA angle is the angle of right-handed rotation about C-N bond, the value being zero if CA-C bond of the preceding residue is cis to N-CA bond.  Most amino acid residues in a protein are involved in two peptide bonds.  The peptide bond formed by the residues i and i+1 is assigned to the residue i + 1.  The same applies to the OMEGA angle.  For that reason no omega angle is assigned to the first residue  (from GARLIC docs).
  
  ScreamMatrix set_OMEGA(double angle);     ///< Sets OMEGA angle.
                                ///< Rotates all atoms off the C end beginning CA(i+1). 
            
  double chi1() const;		///< returns chi1 angle.  returns the value 1000 if chiN() angle is not defined for this particular sidechain.
  double chi2() const;
  double chi3() const;
  double chi4() const;
  double chi5() const;
                    

  /* Misc. functions.  Some get functions. */
  SCREAM_ATOM* get_N() const;   ///< yeah, i should just write a general get function.
  SCREAM_ATOM* get_CA() const; 
  SCREAM_ATOM* get_CB() const;	///< returns pointer to CB atom.  returns NULL if CB not found (Gly).

  int number_of_atoms() const;                    ///< Returns the number of atoms contained in this AminoAcid/Residue.
  double total_charge() const;                    ///< Returns total charge on this AminoAcid/Residue.
  //int get_aa_state() const {return aa_state;};		///< Returns aa_state.
  //  multimap<string, SCREAM_ATOM*> get_nterm_mm() {return nterm_mm;}; ///< Returns nterm_mm.
  //  multimap<string, SCREAM_ATOM*> get_cterm_mm() {return cterm_mm;};///< Returns cterm_mm.

  AARotamer* get_rot() const {return rot;} ;

  ScreamAtomV get_sc_atoms() const;
  ScreamAtomV get_bb_atoms() const;

  /* Output to filehandle functions. */

  void append_to_filehandle(ostream*) const;      
  void pdb_append_to_filehandle(ostream*) const;
  void append_to_ostream_connect_info(ostream*) const;

  /* Obsolete Functions. */

  void print_Me() const;
  void print_ordered_by_n() const;


protected:

  AARotamer* rot;
 
  //int aa_state;

  //  multimap<string, SCREAM_ATOM*> nterm_mm; ///< Store atom info for N terminal atoms.
  //  multimap<string, SCREAM_ATOM*> cterm_mm; ///< Store atom info for C terminal atoms.


};

/*
class NTerm_AA : public AminoAcid {

public:

};


class CTerm_AA : public AminoAcid {

public:

};

*/
#endif  /* SC_AMINOACID_HPP */
