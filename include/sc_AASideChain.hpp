/* sc_AASideChain.hpp
 *
 * Header file for classes relevant to the AASideChain class in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */


#ifndef SC_AASIDECHAIN_HPP
#define SC_AASIDECHAIN_HPP

#include <vector>
#include <map>
#include <string>
using namespace std;

#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"
#include "sc_SideChain.hpp"

/** AASideChain: base class for 20 different types of SideChains.
 * This design of this class is quite involved; probably unnecessarily complicated..  operator=() here is virtual, and in order to generate the right derived objected runtime (the AlaSC's etc) extensive use of the factory() function is applied.  Factories return a pointer to a newly constructed object, which can then be populated (if object is empty) by using operator=().  For careful details, see comments under operator=() and the static factory() functions.
 */


class AASideChain : public SideChain {

public: 
  AASideChain();

  /** Constructors. 
   * Constructors that take a list of SCREAM_ATOM's as argument do not allocate memory on free store.  sc_atoms_mm maps ATOM label to a pointer to the corresponding atom in the list.
   *  Constructors that take a list of strings as part of argument allocates memory from heap.  This means a list of SCREAM_ATOM structures is created, and sc_atoms_mm maps ATOM label to the corresponding atom in that newly created list.  When the object goes out of scope, the destructor destroys those SCREAM_ATOM's.
   */

  AASideChain(const ScreamAtomV&);
  AASideChain(const AASideChain&);
  virtual ~AASideChain();

  virtual AASideChain& operator=(const AASideChain&);
  //virtual AASideChain& deepcopy(const AASideChain&); // aasidechain has no business creating their own SCREAM_ATOM's.


  /* The following set of functions use atom label information.  A function will be written at some point to automatically determine atom labels from connectivity.  If structure read in from SCREAM_ATOM, no need to call these functions. */

  virtual void assign_charges(string) { cout << "I'm in non-overridden assign_charges" << endl;}; ///< Assigns charges according to a particular scheme, such as CHARM22.  
  virtual void assign_atom_fftype() {};   ///< Assigns atom type given atom labels.
  virtual void assign_lone_pair() {};     ///< Assigns lone pair values from atom labels.
  virtual void place_hb_hydrogen() {};    ///< Places hydrogens on hydrogen bonding hydrogens (donors).  E.g.: Ser, Thr, Cys, Tyr.

  /** find_atom_by_greek_name(string): Returns atom on the sidechain whose position is denoted by the greek alphabet argument.  Returns NULL if not found.
   *  This function is written for easy access of sidechain atoms when calculating chi angles.  The atom which defines the chi angles is specified as follows:
   *    - when there are two atoms at the same greek alphabet position (say OE1 and NE2), S has highest precedence, O next, N next, C last.
   *    - X1 always taken in favor of X2.  This actually implicitly takes care of the preceding rule because of nonmenclature rules.
   */

  SCREAM_ATOM* get_atom_by_greek_name(string);  

  


protected:


  AASideChain& dummy_assignment_operator(const AASideChain&); ///< Dummy assignment operator that operator= calls.  
                                                              ///< Reason: operator= overridden in derived class, yet want to use behavior of operator= from this class.  No way this can be done, hence introduce dummy_assignment_operator.  This is called by operator= in derived classes (e.g. AlaSC).
  
//  double private_chi(int, SCREAM_ATOM*, SCREAM_ATOM*);	///< Private function that calcs chi angles.  Int argument denotes which chi angle to calc.  The two SCREAM_ATOM*'s needed for backbone N and CA.

};

#endif /* SC_AASIDECHAIN_HPP */
