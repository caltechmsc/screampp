/*
 * scream_atom.hpp  first draft 1-19-03
 * made into a class 6-12-03.
 *
 * Defines the structure of the scream atom.
 * 
 * Copyright (c) 2003, Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#ifndef SCREAM_ATOM_HPP
#define SCREAM_ATOM_HPP

#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
using namespace std;

/** SCREAM atom definitions.  See inline comments.
 *
 */

/* Explanation on atom flags:
   bit 0x1: Original assignment of fixed and moveable.  1 == fixed (i.e. backbone or sidechains not being SCREAM'ed), 0 == moveable (sidechains or ligands designiated to be SCREAM'ed).

   All of the below bits only affects residue-residue interactions.

   bit 0x2: Energy interaction calculation (sc-sc): like energy expression setup in MD codes, but only applies to residue-residue interactions.  1 == fixed, 0 == moveable.  
   bit 0x4: Visibility during energy itneraction calculation (sc-sc).  0 == invisible, 1 == visible.  If an atom is invisible, it wouldn't be included residue-residue interaction calculations. 
   bit 0x8: Visibility during sidechain-fixed (i.e. empty lattice) calculations.  0 == invisible, 1 == visible.  If an atom is invisible, it wouldn't be included in sidechain-fixed interaction calculations.
 */

class SCREAM_ATOM {

public:


  SCREAM_ATOM();
  SCREAM_ATOM(const string);
  SCREAM_ATOM(SCREAM_ATOM*);
  ~SCREAM_ATOM();

  void pdb_init(string ); ///< initializes using a PDB line.

  /* get set functions */
  void set_x(double x) {this->x[0] = x;};
  void set_y(double y) {this->x[1] = y;};
  void set_z(double z) {this->x[2] = z;};

  string getAtomLabel() { return this->atomLabel; };
  void setAtomLabel(string);
  string getAtomType() { return this->atomType; };
  void setAtomType(string);
  double getX() { return this->x[0];};
  void setX(double x) { this->x[0]=x;};
  double getY() { return this->x[1];};
  void setY(double y) { this->x[1]=y;};
  double getZ() { return this->x[2];};
  void setZ(double z) { this->x[2]=z;};
  double getCharge() { return this->q[0];};
  void setCharge(double q) { this->q[0] = q;};
  string getResName() { return this->resName;};
  void setResName(string s) { this->resName = s;};
  string getChain() { return this->chain; };
  void setChain(string chn) { this->chain = chn; };
  int getResNum() { return this->resNum; };
  void setResNum(int resNum) { this->resNum = resNum; };

  /* Global and static counts to keep track of atoms. */

  int GLOBAL_ATOM_N;                  ///< global atom n.
  //  static int GLOBAL_ATOM_COUNT;       ///< total atom instances encountered

  /* attributes unique to SCREAM atom: */

  string keyw;       ///< indicates whether this is a HETATM or ATOM  .

  string atomLabel;         ///< atom label fof this atom.  must look up table!
  string stripped_atomLabel; ///< atom label stripped for this atom.
  int isSC_Flag;            ///< if ==0: not a SC atom
                            ///< if ==1: a SC atom (CB is a SC atom)
  int isAAResAtom;          ///< if ==1: a amino acid residue atom
  
  string atomType;          ///< atom type from forcefield
  string stripped_atomType; ///< whitespace stripped atom type.

  double occupancy;		///< occupancy field of pdb file
  double BFactor;		///< BFactor field of pdb file.

  string resName;           ///< 3 letter name of amino acid this atom belongs to
  string oneLetterResName;  ///< the 1 letter res name of the amino acid.
  string chain;             ///< chain that this atom belongs to 
  int resNum;            ///< the residue number this atom belongs to


  int atoms_connected;        ///< substr(68,1) of bgf line. 
  int lone_pair;            ///< substr(70,1) of bgf line.

  /* attributes with names that correspond directly to MPSim ATOM struct: */

  double x[3];               ///< coordinates of this atom
  double q[2];               ///< charge and VDW attraction "charge"
  int n;                     ///< global atom number
  int type;                  ///< type from forcefield
  int flags;                 ///< fixed atom flag.  See header of this file for explanations.
  double m;                  ///< mass; useful for calc. moments
  double vchg2;              ///< van der waals repulsion "charge"

  /* attributes that have to do with which rotamer library the atom came from */
  string library_name;       ///< where this atom came from.  if native, == "".  else, values like "SCWRL" or "10".

  /* attributes with names that correspond directly to MPSim FF_ATOM struct: */
  
  double vdw_r;                   ///< vdw radius
  double vdw_d;                   ///< vdw well depth
  double vdw_s;                   ///< vdw scale factor (not used in MPSim)
  double vachg;                   ///< van der waals attraction
  double vrchg;                   ///< van der waals repulsion
  int hb_da;                      ///< hydrogen bond donor acceptor integer association.  See HB_EE for more info.  0 if H___A, positive number if acceptor or donor base, -1 if nothing.
  int a; ///< for problems that need to store information, "a" is the variable for any information to store.
  
  /* attributes that are involved with full delta method for.  stored here for fast access instead of using external lookup table. */
  double delta; // delta is uncertainty of position, adjusted.
  //  double mu; // mu is avareage uncertainty of position; kind of.  or the "uncertainty" that you want to account for.  preinitialized.
  //double sigma; // by how much standard deviation of uncertainty to adjust

  /** Connectivity information stored in this atom.
   *  Corresponds to BOND struct in SCREAM.
   */
  
  map<SCREAM_ATOM*, int> connectivity_m;   ///< connecitivity info.  First element of the pair is pointer to connected atom, second element is bond type.  This is actually pretty fast.


  /** functions involved with this atom
   */
  
  void initFlag() { int least_sig=this->flags&1; least_sig += least_sig<<1; (this->flags|=0x3)&=least_sig; this->flags += 0x4; this->flags += 0x8;} ///< Initializes atom calculation flag flag value to original fixed/moveable value.  0x4 bit is set to 1: visible.  0x8 bit is set to 1: EL visible.
  void resetFlag() { initFlag(); }; ///< Resets flag value to original flag value, based on least significant bit.  Same as above, but different concept.

  void make_atom_moveable() { this->flags &= ~0x2;}; ///< Make atom moveable in residue-residue energy expression.
  void make_atom_fixed() { this->flags |= 0x2;}; ///< Make atom fixed in energy expression.

  void make_atom_invisible() { this->flags &= ~0x4;}; ///< Make atom visible in residue-residue energy expression.
  void make_atom_visible() { this->flags |= 0x4;}; ///< Make atom invisible in residue-residue energy expression.

  void make_atom_EL_invisible() { this->flags &= ~0x8;};
  void make_atom_EL_visible() { this->flags |= 0x8;};

  bool is_part_of_EE() { return this->flags & 0x2;}; ///< Returns true if this atom is part of energy expression.

  void fix_atom(bool value) {if (value) { this->flags = 1; } else { this->flags = 0;} };  ///< Used when setting atoms for the very first time.

  double distance(SCREAM_ATOM*); ///< Calculates Euclidean distance between two SCREAM_ATOM atoms.
  double distance_squared(SCREAM_ATOM*);  ///< Calculates distance squared between two SCREAM_ATOM atoms.  Just to save the time needed to take square roots.

  double worst_clash_dist(vector<SCREAM_ATOM*>&, SCREAM_ATOM** = NULL); ///< returns worst clash dist between this atom and a list of atoms.  The other SCREAM_ATOM* that gets passed in stores info for that worst clash atom.

  void feed_me(const string bgf_line); ///< populates this atom from bgf_line.
  void feed_me_pdb(const string pdb_line); ///< pdb equivalent of feed_me(string bgf_line).
  bool make_bond(SCREAM_ATOM*, int = 0);  ///< makes bond with SCREAM_ATOM with type int.  returns true if successful.  Bond type has default value zero--for now.  Meant to be bond order, like in BIOGRAF.
  bool delete_bond(SCREAM_ATOM*);

  /* Trying to avoid operator overloading to facilitate SWIGging. */
  SCREAM_ATOM& copy(const SCREAM_ATOM&); ///< Copies coordiantes, charges, etc from one atom over to another.  Bond connectivity is not correctedly treated!  So use with caution.
  //  SCREAM_ATOM& deepcopy(const SCREAM_ATOM&); ///< Copies coordinates, charges, etc from one atom over to another.
  SCREAM_ATOM& copyJustCoords(const SCREAM_ATOM&); ///< Copies just coordinates.


  /** Print routines.
   *
   */
  void dump() const;                                  ///< prints bgf line to cout.
  void pdb_dump() const;                              ///< prints pdb line to cout,
  string return_bgf_line() const;
  string return_pdb_line() const;
  void append_to_filehandle(ostream*) const;   ///< appends bgf output to filehandle for atom info.
  void pdb_append_to_filehandle(ostream*) const;   ///< appends pdb output to filehandle for atom info.
  void pdb_append_to_ostream_connect_info(ostream*) const; 
  void append_to_ostream_connect_info(ostream*) const;


private:
  static std::map<std::string, std::string> AA321_map;
  void _init_AA321_map();
};


#endif /* SCREAM_ATOM.HPP */

