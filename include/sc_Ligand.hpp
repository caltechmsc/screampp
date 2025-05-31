/** \file sc_Ligand.hpp
 * Ligand class declared here.  Derived from ProteinComponentAbstract.
 */

#ifndef SC_LIGAND_HPP
#define SC_LIGAND_HPP

#include "scream_atom.hpp"
#include "sc_ProteinComponent.hpp"

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

/** Ligand.  
 * Ligand class.  Non-water molecules are mostly ligands.  
 */

class Ligand : public ProteinComponent {

public:

  Ligand();
  Ligand(const Ligand&);	///< Copy constructor.
                                ///< Copy constructors of classes derived from ProteinComponentAbstract are special in that a deep copy is done, meaning that the object being constructed here have pointers that point to newly created copies of SCREAM_ATOM structs identical to those originally pointed to by reference object that is passed in.  
                                ///<  Note that this differs from the copy constructors of ProteinComponent, the surrogate class. 
  Ligand(const vector<SCREAM_ATOM*>&);
  ~Ligand();			///< Destructor.  Kills all objects pointed to by variables of this class.

  vector<SCREAM_ATOM*> getAtomList() const;
  
  void print_Me() const;
  void append_to_filehandle(ostream* ofstream_p) const;
  void append_to_ostream_connect_info(ostream* ofstream_p) const;

  ProteinComponent* copy() const; ///< Provides definition for the pure virtual copy() function in the base ProteinComponentAbstract class.
  string whatAmI() const { return string("Ligand"); };

private:

  bool ligand_atoms_on_free_store;
  multimap<string, SCREAM_ATOM*> lig_mm;

};

#endif /* SC_LIGAND_HPP */
