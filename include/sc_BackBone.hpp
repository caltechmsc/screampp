/* sc_BackBone.hpp
 *
 * Header file for classes relevant to BackBone in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#ifndef SC_BACKBONE_HPP
#define SC_BACKBONE_HPP

#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"

#include <vector>
#include <map>
#include <string>
using namespace std;

/** BackBone class: Base class for AABackBone etc.
 * 
 */

class BackBone {

public:
  BackBone();
  BackBone(const ScreamAtomV&); 	///< Creates a BackBone object from a list of atoms.
  virtual ~BackBone() {};

  /** Toggles atom fix flag.  True, fixed.  False, moveable.
   *
   */

  void fix_toggle(bool);

  /** Appends output to filehandle.
   *
   */

  void append_to_filehandle(ostream* ofstream_p) const;

  /** Appends pdb output to filehandle.
   *
   */

  void pdb_append_to_filehandle(ostream* ofstream_p) const;

  /** Appends connectivity info of backbone to filehandle.
   *
   */

  void append_to_ostream_connect_info(ostream* ofstream_p) const;

  /** Translates the positions of backbone atoms by values provided by argument ScreamVector.
   * Positions of atoms in backbone are, of course, changed.  Since the atoms are not physically stored as member variables of this class but are instead referred to by pointers, modifications are made to variables OUTSIDE of this class.  Therefore care must be taken when performing this operation: only BackBones that are disposable should call this function. 
   */
  void translate(const ScreamVector&);

  /** Transforms the positions of all backbone atoms by the transformation matrix.
   * Positions of atoms in backbone are, of course, changed.  Since the atoms are not physically stored as member variables of this class but are instead referred to by pointers, modifications are made to variables OUTSIDE of this class.  Therefore care must be taken when performing this operation: only BackBones that are disposable should call this function. 
   */

  void transform(const ScreamMatrix&);

  /* Misc. functions. */


  /** Returns SCREAM_ATOM* of given atom label.
   * Shorthand for backbone->bb_atom_mm.find(string a)->second, which since bb_atom_mm is private is illegal anyway.
   */
  SCREAM_ATOM* get(const string) const;

  int number_of_atoms() const {return bb_atom_mm.size();};    ///< returns the number of atoms in this backbone structure.
  double total_charge() const;                                ///< return total charge on this backbone.

  multimap<string, SCREAM_ATOM*> get_bb_atom_mm() const {return bb_atom_mm;};
  vector<SCREAM_ATOM*> get_atoms() const;

  void copy_atom_positions(const BackBone*);   ///< returns a vector of pointers to atoms in this backbone.

  /* Obsolete Functions. */

  void print_Me() const;   ///<Prints bgf lines as ordered alphabetaically by atom labels.
  void print_ordered_by_n();  ///< Prints bgf lines as ordered by global atom number.  Different from print_Me(), which prints unsorted; or rather, sorted alphabetically by atom labels.
   
protected:

  multimap<string, SCREAM_ATOM*> bb_atom_mm; ///< key = atom_name, value = SCREAM_ATOM*

};


#endif /* SC_BACKBONE_HPP */
