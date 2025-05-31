//!AARotamer.hpp
/**
 * AARotamer.hpp
 *
 */

#ifndef AAROTAMER_HPP
#define AAROTAMER_HPP

#include "defs.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"
#include "sc_BackBone.hpp"
#include "sc_AABackBone.hpp"
#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"
//#include "sc_AASideChain_20AA.hpp"
#include "Rotamer.hpp"

/** AARotamer class.  Derived from Rotamer.
 *
 */

class AARotamer : public Rotamer { 

public:
  // Temporary friend declaration for easy (and dirty SWIGGING).
  friend class Protein;
  friend class AminoAcid;
  friend class AAChain;

  
  ///< Constructors: like SideChain and BackBone classes, those that take in strings constructor new SCREAM_ATOM's on free store.  Those that pass in SCREAM_ATOM lists or AABackBone and AASideChain list simply directs pointers to the already existing structures.

  AARotamer();
  AARotamer(const ScreamAtomV&); ///< Construct an AARotamer by passing a list of SCREAM_ATOMs.  Note: if the passed in ScreamAtom goes out of scope, atoms defined in this AARotamer goes out of scope as well--meaning only pointers to ScreamAtoms in ScreamAtomV are stored and copied, i.e. no deepcopying SCREAM_ATOMs is done.  Caution Caution Caution!
  //  AARotamer(const string, const string); ///< 1st arg: Rotamer file.  2nd arg: Rotamer connectivity file.
  AARotamer(const stringV&); ///< This is not obsolete.  A AARotamer needs to capability to read strings when it's reading a Rotamer library.
  AARotamer(AABackBone* const, AASideChain* const);  ///< This is seldom used.
  //  AARotamer(const AARotamer&);	///< This again, is seldomly used.

  ~AARotamer(); 	///<  Not necessary to be virtual.  No derived classes; else too convoluated.

  /** Deepcopy.
   *  Recursively deepcopies all SCREAM_ATOM structures.  I.e. doesn't just copy references/pointers to SCREAM_ATOM's.
   */

  AARotamer& deepcopy(const AARotamer&);

  /* Misc. get/set Functions. */

  string get_resName() const {return this->resName;};
  void set_resName(std::string s) { this->resName = s;};
  void initRotamerAtomList(const vector<string>&);

  /* Misc. angle calculating routines. */

  double calc_PHI();  ///< Calculates PHI angle from positions of HN, N, CA, C.  Note: pstn of HN.  Sets PHI to calculated value and returns this value. 
  double calc_PSI();  ///< Calculates PSI angle from positions of N, CA, C, O.   Note: pstn of O.  Sets PSI to calculated value and returns this value.

  double get_PHI() {return PHI;};
  double get_PSI() {return PSI;};

  double chi1();        ///< Returns chi1 angle, if defined.  Else, returns 1000.  (error)
  double chi2();	///< Returns chi2 angle, if defined.   Else, returns 1000.  (error)
  double chi3();	///< Returns chi3 angle, if defined. Else, returns 1000.  (error)
  double chi4();	///< Returns chi4 angle, if defined. Else, returns 1000.  (error)
  double chi5();	///< Returns chi5 angle, if defined. Else, returns 1000.  (error)

  /** Overrides match_bb from base class.
   * Make sure a AARotamer is passed in, not just any Rotamer.
   */
  
  void match_bb(const Rotamer* const);

  /** Other rotamer placement methods. */
  void match_CB(const Rotamer*, const SCREAM_ATOM* const, double); ///< place rotamer by matching first CB atom.  Last double: rotamerMatchVector (0.5 means matching the C-CA-N bisector of the rotamer and the backbone).
  void create_CB(vector<double>&, SCREAM_ATOM*); ///< the 4 doubles: offBisectorAngle, offPlaneAngle, bondLength, rotamerMatchVectorLambda.  SCREAM_ATOM*: a SCREAM_ATOM* that gets passed in to store the x, y, z coordinates of the created CB.

  /** Assigns atom FF type from atom labels (Dreiding).
   * 
   */

  void assign_atom_fftype();

  /** Assigns charges on atoms from atom labels according to method.  String argument taken; can be "CHARM22".
   *  Int tells the class whether it's a C-term or N-term aminoacid.  aa_state from sc_AminoAcid.
   */

  void assign_charges(string, int = 0);

  /** Assigns lone pair and atom connection info on atoms.
   *
   */

  void assign_lone_pair();

  /** Calculates the position of C(i-1) given N(i), HN(i) and CA(i).
   * These atoms positions of included in bb.  Calls function in AABackBone.
   */

  ScreamVector calc_C_i_minus_one();


  void center_CA();           ///< Translates coordinates of rotamer such that CA = (0,0,0).

  /* Print functions. */

  void append_to_filehandle(ostream*) const;
  void pdb_append_to_filehandle(ostream*) const;
  void append_to_ostream_connect_info(ostream*) const;

  /* Obsolete Functions. */

  void print_Me() const;
  void print_ordered_by_n() const;


  //protected:

  AABackBone* get_bb() const { return (AABackBone*)(this->bb);};   ///< overriding with derived return type allowed in C++.  not allowed on MS Visual C++ though.
  AASideChain* get_sc() const { return (AASideChain*)(this->sc);}; ///< overriding with derived return type allowed in C++.  



  double PHI;                          ///< PHI angle from rotamer library.
  double PSI;                          ///< PSI angle from rotamer library.
  string resName;              ///< three letter abbrev. of residue name

  double private_chi(int); ///< Private function that calcs chi angles.  Int argument denotes which chi angle to calc.  
                           ///< Error checking handled here.


  ScreamAtomV _determine_and_fix_GLY_sidechain_HCA_issue(ScreamAtomV&); ///< a helper function that takes in a list of SCREAM_ATOM*'s of a Glycine residue (i.e. N, HN, CA, HCA x 2, C, O) and returns the HCA atom that occupies the L-amino acid position, wrapped in a ScreamAtomV.  Also removes the L-amino acid HCA from bb_atom_v.

};


#endif /* AAROTAMER_HPP */
