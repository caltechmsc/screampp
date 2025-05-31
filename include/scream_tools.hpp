/**\file scream_tools
 * package includes the following helper functions:
 *   - calc_dihedral(arg1, arg2, arg3, arg4)
 *
 */

#ifndef SCREAM_TOOLS_HPP
#define SCREAM_TOOLS_HPP

#include <math.h>
#include <vector>
#include <map>
#include <string>

using namespace std;

#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"

namespace scream_tools {

  string atom_label(const string); ///< Returns atom type from a bgf line.
  int res_number(const string);	///< Returns residue number from a bgf line.
  string res_type(const string); ///< Returns the type of residue from a bgf line.
  string chain_desig(const string); ///< Returns the chain designation from a bgf line.
  string strip_whitespace(string); ///< Strips leading and trailing white spaces from string (please, one contiguous word).
  vector<string> split(string);    ///< Perl-style split string.
  vector<string> char_double_star_to_str_of_vector(char*[]);  ///< Used in scream_interface.h .cpp.  Takes in an array of pointers to char arrays and returns a vector of strings.
  vector<string> mut_list_string_to_str_of_vector(char*);    ///< Used in scream_interface.  Takes in a string of mutation/rerplacement candidates.


  bool is_bgf_header_line(string); ///< Returns true if line is a header line in a bgf file.
  bool is_bgf_atom_line(string);   ///< Returns true if line is an atom (or hetatm) line in bgf file.
  bool is_format_connect_line(string); ///< Returns true if line is a FORMAT CONECT or FORMAT ORDER line.
  bool is_connectivity_info(string);   ///< Returns true if line is a CONECT or ORDER line in bgf file.
  

  bool is_AA(const string);        ///< Returns true if res_type(bgf_line) is an AA atom.  Argument: 3 letter abbrev.
                                   ///< Inner working: if res_type is one of 20 AA's, return true.  else return false.

  bool is_metal_atom(const string);     ///< Returns true if atom is metal atom.
  bool is_SC_atom(const string);	///< Returns true if the atom is a (aminoacid) sidechain atom.  Takes in atom_label.
                                ///< This function draws on the fact that all aminoacid sidechain atoms are branched off the C alpha carbon.  Therefore all sidechain atoms are labeled something beta or something gamma or beyond.
  bool is_BB_atom(const string); ///< Returns true if atom is a backbone atom.  Takes in atom label.
  bool is_heavy_atom(const string); ///< Returns true if atom is a heavy atom.

  bool is_N_term_hydrogen(const string); ///< Returns true if atom is a N term hydrogen.
  bool is_C_term_atom(const string); ///< Returns true of atom unequivocaly indicates that this amino acid is a C-term.
  
  //  bool is_nonpolar_hydrogen(const string); ///< Returns true if query atom label 

  bool is_HCA_atom(const string); ///< REturns true if query string is atom label for a HCA atom.

  bool is_natural_AA(const string); ///< Returns true is string provided (3 letter AA) is a natural AminoAcid.
  
  bool AA_atom_order(SCREAM_ATOM*, SCREAM_ATOM*); ///< Returns true if SCREAM_ATOM* 1 comes before SCREAM_ATOM* 2.

  ScreamAtomV return_heavy_atoms(const ScreamAtomV& ); ///< Returns a list of heavy atoms from atoms passed in.


  /** deep copies a multimap with its values being pointers to SCREAM_ATOMs.  
   * This means the SCREAM_ATOMs that the pointers point to are reconstructed.  This is useful in many occasions: for instance, copying BackBone instances that were initialized by vector<string>.
   */


  string get_AA_from_mutationInfo(const string); ///< Returns the AA type to mutate to from mutationInfo string (eg C218_X)
  int get_pstn_from_mutationInfo(const string); ///< Returns the pstn to do mutation from mutationInfo string (eg C218_X)
  string get_chn_from_mutationInfo(const string); ///< Returns the chain on which to do mutation from mutationinfo string ( eg C218_X)

  multimap<string, SCREAM_ATOM*> deep_copy_str_SCREAM_ATOM_mm(const multimap<string, SCREAM_ATOM*>&, SCREAM_ATOM** = NULL);
  void deep_copy_ScreamAtomV(const ScreamAtomV&, ScreamAtomV&); // Does not take care of memory deletion!  User beware: when using this function, make sure you delete SCREAM_ATOM*'s at end.  Also, make sure you pass an empty new_atom_list in; otherwise, issues error.

  string three_letter_AA(const string);    ///< Returns the three letter AA designation.
  string one_letter_AA(const string); ///< Return the one letter AA designation.

  /** Calculates dihedral angle in degrees, given 4 ScreamVectors.  
   * The dihedral angle is defined as the angle of right-handed rotation around the ray starting from the 2nd ScreamVector argument to the 3rd ScreamVector argument, the value being zero if 1st ScreamVector argument and 4th ScreamVector argument are on the same side of the ray.  Range: -180 to 180 degrees.
   */
  double calc_dihedral(const ScreamVector&, const ScreamVector&, const ScreamVector&, const ScreamVector&);

  /** Calculates dihedral angle, given 4 SCREAM_ATOM* 's.
   *  The dihedral angle is the angle of right-handed rotation around the bond (or line) connecting the 2nd and 3rd atom, the value being zero if the bond (or line) between the 1st and 2nd atom is cis to the bond (or line) between the 3rd and 4th atom.  Range: -180 to 180 degrees.
   */
  double calc_dihedral(const SCREAM_ATOM* const, const SCREAM_ATOM* const, const SCREAM_ATOM* const, const SCREAM_ATOM* const);

  /** Translates position of an atom.  Pointer is passed in, position is updated.  Pointer is const, i.e. address doesn't change.
   * This function is overloaded in many ways.
   */

  void translation(SCREAM_ATOM* const, const ScreamVector&);

  /** Translates position of all atoms in a vector of atoms.
   * This function is overloaded in many ways.
   */

  void translation(const vector<SCREAM_ATOM*>&, const ScreamVector&);

  /** Translates position of all atoms in a multimap of atoms.
   * This function is overloaded in many ways.
   */

  void translation(const multimap<string, SCREAM_ATOM*>&, const ScreamVector&);
  
  /** Returns a ScreamMatrix object that transforms toBeMoved atom positions to stationary atom positions to the best of its ability.
   *
   */

  pair<ScreamMatrix, ScreamVector> getTransformationPairFromAtoms(const ScreamAtomV& stationary, const ScreamAtomV& toBeMoved, bool exactMatch = 0);
  
  /** Calculates the CRMS distance between one listt of atoms and a second list of atoms.
   */
  double distance(const ScreamAtomV&, const ScreamAtomV&); 

  /** Calculates the CRMS distance between one map of atoms and a second map of atoms.  Assuming: no 
   */
  double distance_by_atom_label_map(const map< string, SCREAM_ATOM* >&, const map< string, SCREAM_ATOM* >&);

  /** Calculates the maximum distance between two equivalent atoms on two lists of atoms.
   */
  pair<double, string>  max_equivalent_atom_dist(const map< string, SCREAM_ATOM* >&, const map< string, SCREAM_ATOM* >&);

  /** Adds a HN atom when doing a mutation from PRO to NON-PRO. 
   *  HN atom passed in; only coordinates are calculated.  Connectivities, ff_types, atom_n, etc need to be taken care of elsewhere.
   */
  void calc_new_HN_atom_coords(const SCREAM_ATOM*, const SCREAM_ATOM*, const SCREAM_ATOM*, SCREAM_ATOM*); 
  
  /** Makes connectivity_m map of a for a vector of SCREAM_ATOM objects from a connectivity map.
   *  Used primarily for making connecitivies for Rotamer atoms using information from RotConnInfo.
   */
  void make_connectivity(ScreamAtomV&, const map< int, vector<int> > &);

  /** Checks 1-2 and 1-3 connectivity. 
   */

  int should_exclude_on_1_2(SCREAM_ATOM*, SCREAM_ATOM*); // returns non-zero (1) if the two atoms are bonded (1-2)
  int should_exclude_on_1_3(SCREAM_ATOM*, SCREAM_ATOM*); // returns non-zero (1) if the two atoms are 1-3 bonded

  int should_exclude_on_1_2_or_1_3(SCREAM_ATOM*, SCREAM_ATOM*); // returns non-zero (1) if either bonded (1-2) or 1-3.

  /** Hydrogen generation routines. 
   */

  vector<ScreamVector> generateHydrogenCoords(SCREAM_ATOM*); ///< generates coord for new hydrogen.
  ScreamAtomV createHydrogens(vector<ScreamVector>&, SCREAM_ATOM*); ///< Creates (new's) hydrogens and make connectivity.  SCREAM_ATOM*: the base atom.

  //vector<ScreamVector> generateNHHydrogenCoords(ScreamVector&, ScreamVector&, ScreamVector&); ///< generates the mainchain NH hydrogen if missing.
  vector<ScreamVector> generate3SP3HydrogenCoords(ScreamVector&, ScreamVector&, ScreamVector&); ///< generates 3 sp3 positions.  first vector: base atom.  second: connected atom.  third: a 1-3 atom, to get the right staggered conformation, with the staggered angle always being 60 deg, sp2 or sp3 on the connected atom.  DREIDING: sp3-sp3 min 60deg, sp3-sp2 min also 60 deg.

  vector<ScreamVector> generate2SP3HydrogenCoords(ScreamVector&, ScreamVector&, ScreamVector&); ///< generates 2 sp3 position hydrogens.  First ScreamVector: base atom coord.  Second: a connected atom coords.  Third: a second connected atom coord.  generates 2 sp3 hydrogens with 2 of the 4 possible positions already in place.
  vector<ScreamVector> generate1SP3HydrogenCoords(ScreamVector&, ScreamVector&, ScreamVector&, ScreamVector&); ///< generates 1 sp3 position hydrogen.  first: base.  2nd through 4th: connected atoms.  Only for Proline.
  vector<ScreamVector> generate2SP2HydrogenCoords(ScreamVector&, ScreamVector&, ScreamVector&); ///< generates 2 sp2 position hydrogens.  First ScreamVector: base atom coord.  Second: a connected atom coords.  third: a 1-3 atoms.  120 deg assumed.
  vector<ScreamVector> generate1SP2HydrogenCoords(ScreamVector&, ScreamVector&, ScreamVector&); ///< generates 1 sp2 position hydrogen.  first: base.  second: first connected atom.  third: second connected atom.  120 degree assumed when possible.

};

class AtomNotFoundException {
public:
  const string what;
  
  AtomNotFoundException(const std::string& i_what): what(i_what){};

};
  
#endif /* SCREAM_TOOLS_HPP */
