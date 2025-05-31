
/* sc_ProteinComponent.hpp
 *
 * Header file for classes relevant to ProteinComponent in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

/*****************************************
 ** P Protein component: *****************
 ** simply what a protein is made up *****
 ** made up of.  a chain of residues, ****
 ** waters, ligands, aminoacids. *********
 ****************************************/

#ifndef SC_PROTEINCOMPONENT_HPP
#define SC_PROTEINCOMPONENT_HPP

#include <fstream>
#include <vector>

#include "scream_atom.hpp"

using namespace std;

/** ProteinComponent is abstract base class for components of a protein, such as chain, ligand, water, etc.
 *
 */

class ProteinComponent {

public:

  //  virtual int number_of_atoms() const = 0;
  //  virtual ProteinComponent* copy() const =0;
  
  virtual ~ProteinComponent() { };

  /** Returns an array of pointers to atoms belonging to this Component (chain, ligand, etc).
   *
   */
  
  virtual vector<SCREAM_ATOM*> getAtomList() const = 0;

  ProteinComponent& operator=(const ProteinComponent&);
  ProteinComponent& merge(ProteinComponent&);

  
  virtual void fix_toggle(bool) {}; ///< toggles fixed atom flag (all atoms).  true == fixed, false == moveable.

  virtual void print_Me() const {};
  virtual void print_ordered_by_n() const {};
  virtual void append_to_filehandle(ostream* ofstream_p) const {};
  virtual void pdb_append_to_filehandle(ostream* ofstream_p) const {};
  virtual void append_to_ostream_connect_info(ostream* ofstream_p) const {};

  virtual string whatAmI() const = 0;  ///< I.e. which kind of proteinComponent am i?  I.e. AAChain, Ligand, Water, HETATMs, etc.

};


#endif /* SC_PROTEINCOMPONENT_HPP */
