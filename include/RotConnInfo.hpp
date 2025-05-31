#ifndef ROTCONNINFO_HPP
#define ROTCONNINFO_HPP

#include <vector>
#include <map>
#include <string>
using namespace std;


/** RotConnInfo.  Contains information about the connectivity of a rotamer conformer.
 *
 */

class RotConnInfo {
  
  /* Stuff relevent to connectivity info and anchor info.
   */
public:
  string targetRotamerLibFile;	///< Name of target Rotamer Library File as specified in rotamer_connectivity file.
  vector<int> anchor_pts;	///< Almost equivalent to sidechain atoms.  These points, under most circumstances, should match PERFECTLY with protein structure points.
  vector<int> atoms_of_exact_match; ///< These points will be exactly matched.
  map<int, int> atom_n_map; 	///< Key: atom number on rotamer provided.  Value: atom number on actual protein.
  map<int, string> atom_n_label_map; ///< Key: atom number on rotamer provided.  Value: atom label for that atom.
  vector<int> side_chain_atoms;	///< An array of atoms that are "sidechain" atoms.  I.e. these are the atoms that form the "sidechain" part of the "rotamer".  I.e. these are the atoms whose positions will be replaced.
  map< int, vector<int> > atom_connectivity_info; // basically a table that contains all connectivities for an atom.  Ints are in "rotamer n space/representation".
  vector<int> connection_point_atoms; ///< Atoms through which the conformer is connected to the rest of the system.  Like CA would be a connection point for CB.
  
  //string connectivityFileName;		///< How this RotConnInfo was created.

  void modifyMappingInProteinAtoms(vector<int>& ); ///< modifies map<int, int> after mutations have been done in Protein.

  void clear() {
    this->anchor_pts.clear();
    this->atoms_of_exact_match.clear();
    this->atom_n_map.clear();
    this->atom_n_label_map.clear();
    this->side_chain_atoms.clear();
    this->atom_connectivity_info.clear();

  }; ///< Deleted all information stored in this structure.

};

#endif /* ROTCONNINFO_HPP */
