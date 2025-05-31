/* sc_AABackBone.hpp
 *
 * Header file for classes relevant to AABackBone in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#ifndef SC_AABACKBONE_HPP
#define SC_AABACKBONE_HPP

#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"
#include "sc_BackBone.hpp"

#include <vector>
#include <map>
#include <string>
using namespace std;


/** AABackBone.  Derived from BackBone.
 * Contains many basic operations on the BackBone of an Amino Acid chain.
 */


class AABackBone : public BackBone {
  
public: 

  AABackBone();
  AABackBone(const SCREAM_ATOM*);                ///< initializes by taking first encounters of bb atoms and populates multimap::bb_atoms.  No new SCREAM_ATOM 's created, only pointers pointing to those SCREAM_ATOM's are stored.
  AABackBone(const vector<SCREAM_ATOM*>&);   ///< same as above, except reads in a vector of SCREAM_ATOM*
  AABackBone(const AABackBone&);                    ///< copy constructor
  ~AABackBone();

  AABackBone& operator=(const AABackBone&); ///< copy assignment
  //  AABackBone& deepcopy(const AABackBone&); ///< deepcopy: copies SCREAM_ATOM's.  No longer allowed; bb has no business creating their own atoms.

  void assign_atom_fftype(); ///< assigns FF type from atom labels.

  /** Calculates the position of C(i-1) given N(i), HN(i) and CA(i).
   * These atoms positions of included in bb.  Algorithm: C(i-1) is placed on the bisector of the HN-N-CA angle 1.32A away from the N atom.
   */
  ScreamVector calc_C_i_minus_one() const; 
  ScreamVector calc_N_i_plus_one() const; ///< Calculates position of N(i+1) from positions of CA, C and O.  Look under calc_C_i_minus_one for more info.

  

  double PHI() const;			///< returns PHI angle in degrees, using C(i-1) position from calc_C_i_minus_one().
  double PHI(SCREAM_ATOM* C_i_minus_one) const; ///< returns PHI angle providing C(i-1) atom position.
                                          ///< An argument is needed because each backbone mononer does not know positions of atoms of other backbone monomers.  This function is therefore called primarily by the CHAIN class.

  ScreamMatrix set_PHI(double, SCREAM_ATOM* C_i_minus_1 = NULL);		///< Sets AABackBone to specified PHI angle and returns the transformation matrix that sets the PHI angle.  See AminoAcid class for more info.
  double PSI() const;			///< returns PSI angle in degrees, using N(i+1) position from calc_N_i_plus_one().
  double PSI(SCREAM_ATOM* N_i_plus_one) const; ///< returns PSI angle in degrees providing N(i+1) atom position.
                                         ///< An argument is needed because each backbone mononer does not know positions of atoms of other backbone monomers.  This function is therefore called primarily by the CAHIN class.
  ScreamMatrix set_PSI(double, SCREAM_ATOM* N_i_plus_1 = NULL);		///< Sets AABackBone to specified PSI angle and returns the transformation matrix the sets the PSI angle.
  
  /* Easy access to specific Atoms.  Interface. */


  SCREAM_ATOM* N() const;
  SCREAM_ATOM* HN() const;
  SCREAM_ATOM* CA() const; ///< returns CA SCREAM_ATOM*.
  SCREAM_ATOM* HCA() const; ///< returns HCA SCREAM_ATOM*.
  SCREAM_ATOM* C() const; ///< returns C SCREAM_ATOM*.
  SCREAM_ATOM* O() const; ///< returns O SCREAM_ATOM*.

private:

  static map<string, string> atom_label_fftype_map;
  static map<string, double> atom_label_CHARM22_map;

  void init_atom_label_maps();

  bool PROLINE_flag;                                ///< True if backbone is part of a proline residue.


};

#endif /* SC_AABACKBONE_HPP */
