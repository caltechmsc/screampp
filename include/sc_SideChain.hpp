/* sc_SideChain.hpp
 *
 * Header file for classes relevant to the SideChain class in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */


#ifndef SC_SIDECHAIN_HPP
#define SC_SIDECHAIN_HPP

#include <vector>
#include <map>
#include <string>
using namespace std;

#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"

/** Class SideChain is the base class of all other SideChain classes.
 * Class SideChain is the base class of all other SideChain classes.  Derived classes include: AASideChain.
 */


class SideChain {

public:
  SideChain();
  SideChain(const ScreamAtomV&);
  virtual ~SideChain();

  /** Toggles atom fixed flag.  If true, fixed; false, moveable.
   *  
   */

  void fix_toggle(bool);

  /** Returns name of SC.  this variable is stored in base classes.
   *
   */

  //  virtual string get_SC_name() const {return string(" ");};


  /** Prints bgf lines as ordered alphabetaically by atom labels.
   *
   */

  void print_Me() const;

  /** Prints bgf lines as ordered by global atom number.
   * Different from print_Me(), which prints unsorted; or rather, sorted alphabetically by atom labels.
   */
  void print_ordered_by_n() const;

  /** Appends output to filehandle.
   *
   */

  void append_to_filehandle(ostream* ofstream_p) const;

  /** Appends pdb output to filehandle.
   * 
   */

  void pdb_append_to_filehandle(ostream* ofstream_p) const;

  /** Appends connectivity info.
   *
   */

  void append_to_ostream_connect_info(ostream* ofstream_p) const;

  /** Translates the positions of sidechain atoms by values provided by argument ScreamVector.
   * Positions of atoms in sidechains are, of course, changed.  Since the atoms are not physically stored as member variables of this class but are instead referred to by pointers, modifications are made to variables OUTSIDE of this class.  Therefore care must be taken when performing this operation: only SideChains that are disposable should call this function. 
   */
  void translate(const ScreamVector&);

  /** Transforms positions of sidechain atoms by transformation matrix.
   * 
   */
  void transform(const ScreamMatrix&);

  /** Calculates rms with respect to another SideChain.  
   *  Care must be taken to ensure the two structures are "comparable".
   */

  double rms(const SideChain*) const;

  /** Calculates rms only for heavy atoms.
   *
   */

  double rms_heavyatom(const SideChain*) const;


  /** Copies coordinates of atoms from input SideChain to current SideChain.
   *
   */

  void copy_atom_positions(const SideChain*);
  
  /** Same as above, except the library name included in the atoms are copied over as well.
   *
   */
  
  void copy_atom_positions_and_library_name(const SideChain*);

  /** Copies just library name.
   *
   */

  void copy_library_name(const SideChain*);

  /** Takes in a vector of SCREAM_ATOM's and returns the distance with the worst clash.
   *
   */

  double worst_clash_distance(const vector<SCREAM_ATOM*>&) const;
  
  /* Misc. functions. */
  void remove_atom(string);
  void add_atom(string, SCREAM_ATOM*);

  SCREAM_ATOM* get(const string) const;

  int number_of_atoms() const {return sc_atom_mm.size();};  ///< returns number of atoms in this sidechain structure.
  double total_charge() const;                              ///< returns total electric charge on sidechain.
  
  multimap<string, SCREAM_ATOM*> get_sc_atom_mm() const {return sc_atom_mm;};
  vector<SCREAM_ATOM*> get_atoms() const;     ///< returns a vector of pointers to atoms in this sidechain.

  string get_SC_name() const {return SC_name;};

protected:

  string SC_name;
  multimap<string, SCREAM_ATOM*> sc_atom_mm; ///< Uses multimap instead of vector.  
				              ///< Want to do this to enable fast access of an atom with a particular label.  Note: string is the atom label of the prescribed SCREAM_ATOM*, with one exception: no whitespaces in atom label.  I.e. to access a SCREAM_ATOM with atom label " CB ", just provide "CB" as key, not " CB ".
  //  vector<string> orig_atom_l_seq; ///< Original atom lable sequence.  
				  ///<  Multimap loses ordering, so must store this info somewhere.  Is this necessary?

private:

};

#endif /* SC_SIDECHAIN_HPP */
