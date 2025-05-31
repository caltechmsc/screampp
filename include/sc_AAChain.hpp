#ifndef SC_AACHAIN_HPP
#define SC_AACHAIN_HPP

#include "defs.hpp"
#include "scream_atom.hpp"
#include "sc_ProteinComponent.hpp"
#include "sc_AminoAcid.hpp"

#include "scream_vector.hpp"
#include "scream_matrix.hpp"
#include "scream_tools.hpp"

#include <map>
#include <vector>
#include <iostream>
using namespace std;


/** AAChain.
 * AAChain class.  I.e. polypeptides.
 */

class AAChain : public ProteinComponent {

public:

  AAChain();
  AAChain(const ScreamAtomV&);
  AAChain(const AAChain&);
  ~AAChain();

  string whatAmI() const { return string("AAChain");}; 

  /* Get SCREAM_ATOM and get lower data structures functions. */

  ScreamAtomV getAtomList() const;        ///< Returns an array of pointers to atoms belonging to this Component (chain, ligand, etc).
  ScreamAtomV get_sc_atoms(int) const;   ///< Returns specified SC atoms in a vector. 
  AminoAcid* operator[](int) const;  ///< Returns pointer to AminoAcid at residue number (int).
                                     ///< If index > size of AAChain, returns NULL.

  /* Object State functions (i.e. gives you information about the state of AAChain; # of AA's, # of ATOM's, etc. */

  int length() const;                ///< Returns the length of this Chain.
  int number_of_atoms() const;       ///< Returns the number of atoms contained in this AAChain.
  double total_charge() const;       ///< Returns the total charge on this AAChain.
  vector<int> getResidueNumbers() const; ///< Returns all residues numbers that are defined on this chain.
  /* Sidechain replacement */

  void SC_replacement(int, const AARotamer* const, string, vector<double>& );  ///< replaces sidechain at position int with AARotamer sidechain.

  /* PHI, PSI, OMEGA functions. */

  double PHI(int) const; ///< Returns PHI value for residue number provided as argument.  Residue number is taken from bgf file.
  //double PHI_offset(int) const; ///< Returns PHI value for residue number specified, but residue number has been offset so that the initial residue is numbered zero.
  double PSI(int) const; ///< Returns PSI value for residue number provided as argument.  Residue number is taken from bgf file.
  //double PSI_offset(int) const; ///< Returns PSI value for residue number specified, but residue number has been offset so that the initial residue is numbered zero.
  //  double OMEGA(int) const; ///< Returns OMEGA value for residue number provided as argument.  Residue number is taken from bgf file.
  //double OMEGA_offset(int) const; ///< Returns OMEGA value for residue number specified, but residue number has been offset so that the initial residue is numbered zero.

  ScreamMatrix set_PHI(int, double); ///< Sets PHI value on residue int.  residue num is consistent with input; i.e. from bgf file or atom info.
                                     ///< Changes positions of all atoms downstream of the residue that has its PHI value set.
  ScreamMatrix set_PSI(int, double); ///< Sets PSI value on residue int.  residue num if consistent with input; i.e. from bgf file or atom info.
                                     ///< Changes positions of all atoms downstream of the residue that has its PSI value set.
  //  ScreamMatrix set_OMEGA(int, double); ///< Sets OMEGA value.  Changes positions of all atoms downstream of the residue that has its OMEGA value set.

  /* Misc. functions. */

  /* Functions related to mutations. */
  void replace_AminoAcid(int, AminoAcid*); ///< Memory management; deletes an AminoAcid and replaces it with another.

  /* MPSim functions.  Not sure if obsolete. */

  void fix_toggle(bool);     ///< makes all atoms on entire chain fixed (true) or moveable (false).
  void fix_toggle_sc_pstn(bool, int); ///< make all atoms on selected rotamer fixed (true) or moveable (false).

  /* Print functions.  Many are obsolete by now.  They are here for debugging purposes. */

  void print_Me() const;
  void print_ordered_by_n() const;
  void append_to_filehandle(ostream* ofstream_p) const;                 ///< prints to file specified.
  void pdb_append_to_filehandle(ostream* ofstream_p) const;
  void append_to_ostream_connect_info(ostream* ofstream_p) const;

  void print_chi_angle_spread(ostream* ofstream_p) const;

  map<int, AminoAcid*> get_aa_m() { return aa_m;};

  string get_chain_desig() { return chain_desig; };

  

private: 

  /* Private Member variables.  AAChain is not meant to be inherited. */

  /* map<int, AminoAcid*>: Bad idea.  Use a AminoAcid number map instead--like in ModSim get_atom. */
  map<int, AminoAcid*> aa_m;               ///< Meat.

  multimap<string, SCREAM_ATOM*> nterm_mm; ///< Store atom info for N terminal atoms.
  multimap<string, SCREAM_ATOM*> cterm_mm; ///< Store atom info for C terminal atoms.

  string chain_desig;                     ///< Chain designation, one letter.  From atom_list <== bgf file.

  /* free store variables.  possibly obsolete because this structure never actually new's SCREAM_ATOM's. */

  bool nterm_mm_on_free_store;            ///< True if nterm_mm SCREAM_ATOM*'s are created using operator new.
  bool cterm_mm_on_free_store;            ///< True if cterm_mm SCREAM_ATOM*'s are created using operator new.

  /* Private member functions. */

  void InitFromAtomList(const ScreamAtomV&);

  void mutation_replacement(int, const AARotamer* const); ///< private function called by SC_replacement to do mutations.

};

#endif /* SC_AACHAIN_HPP */
