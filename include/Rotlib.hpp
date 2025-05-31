/* Rotlib.hpp
 *
 * Header file for classes relevant to rotamer rotamer library in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#ifndef ROTLIB_HPP
#define ROTLIB_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "defs.hpp"
#include "RotConnInfo.hpp"
#include "scream_atom.hpp"
#include "Rotamer.hpp"
#include "AARotamer.hpp"
#include "RotamerCluster.hpp"
#include "sc_BackBone.hpp"
#include "sc_AABackBone.hpp"
#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"

/**
 *  Rotlib.  Non-virtual base class for all rotamer library classes.
 *  Remark (04-06-04): The concept of a rotamer has been generalized to mean a "Confomer".  The Backbone in a Rotamer is really the "AnchorPoints", whereas the SideChain of a Rotamer is the "Variable" part.  I.e., if you provide a dock structure, the protein itself is fixed with the ligand provided a number of different conformations/translocations in the binding site.  We can those these conformations/translocations a "Confomer", or as in the class "Rotamer", which is exactly implemented in this sense (it didn't start off this way but evolved to mean something like this).
 */

class Rotlib {   // originally a namespace

public:

  friend class HIS_NtrlAARotlib;
  friend class Multiple_NtrlAARotlib;
  /* Constructors and destructors */
  
  Rotlib();
  Rotlib(string);		///< Only argument: rotamer connectivity/anchor info file.  
  Rotlib(string, string);	///< First argument stands for the rotamer coordinate file, second argument stands for the rotamer connectivity file, which contains info about connectivity of all rotamers in the rotamer coordinate file and the anchoring information.
  virtual ~Rotlib();

  /* Readers. */

  void readConnectivityFile(string); ///< Reads in a connectivity file and adds appropriate connectivities. This automatically reads the rotamer library file, because rotamer library file is specified in the connectivity file.
  virtual void readRotamerLibrary(string); 	///< Reads in a coordinate library file and populates conformers.

  /* Accessors. */

  string get_library_name() { return this->library_name; } ///< Return the class (i.e. label) for this Rotamer library.

  virtual RotConnInfo* getRotConnInfo()  { return &rotConnInfo;}; ///< Returns rotamer connectivity info.

  virtual Rotamer* get_next_rot();     ///< returns pointer to next Rotamer.  returns NULL if reaches end of Rotlib.
  virtual Rotamer* get_current_rot();  ///< returns pointer to current Rotamer.  returns NULL if reaches end of Rotlib.
  virtual Rotamer* get_next_rot_with_empty_lattice_E_below(double); ///< returns next rotamer with an empty lattice E of below E provided as input.  E here is relative to base state, not absolute E.
  virtual Rotamer* get_empty_lattice_E_rot(int); ///< returns pointer to Rotamer with Empty lattice emtpy ordering of int provided.
  Rotamer* get_empty_lattice_E_rot_after_sorted_by_empty_lattice_E(int); ///< like above, but faster, assumes underlying vector of rotamers already sorted by empty lattice E.

  /* Simple Operations. */

  virtual void reset_pstn();           ///< Current rotamer position is reset to the beginning of rotamer library.

  int size();  ///< Returns the number of rotamers currently in this Rotlib.
  int n_rotamers_below_empty_lattice_energy(double E); ///< Returns the number of rotamers in this Rotlib with an empty lattice energy of below this value.  BY default, lowest state has energy zero.  

  /* External Rotamer adding functions */

  virtual void add_rotamer(Rotamer*); ///< Adds a rotamer to this rotamer library.  virtual: HIS_NtrlAARotlib overloads this.
  virtual Rotamer* new_rotamer(); ///< Returns a new Rotamer object that is bound to rot_v.
  virtual RotamerCluster* new_rotamer_cluster(); ///< Returns a new RotamerCluster object that is bound to rot_v.  needed for SWIG only, else OOJ takes care of thing in new_rotamer().
  /* Test print functions. */

  void print_Me();
  virtual void print_to_file(ofstream*) {};       ///< Outputs rotamer library to ofstream*.

  /** Sorting functions.
   *  Reorders vector<Rotamer*>.
   */

  void sort_by_rotlib_E(); ///< sorts by rotlib_E. Reorders vector<Rotamer*>.
  void sort_by_self_E(); ///< sorts by rotamer self energy in protein environment (backbone emptied).

  void sort_by_empty_lattice_E();   ///< As suggested by function name--sorts by increasing empty lattice energy of a rotamer.  Automatically initializes values of empty_lattice_E from empty_lattice_abs_E.

  double get_best_preCalc_E(); ///< returns the value for the best precalc energy.  used to set energy for existing rotamer from protein.


  /* Hierachical rotamer clustering i.e. family/parent/child approach to rotamer placement related functions. */
  
  

protected:

  string library_name;         ///< name of this library.
  
  int orig_rotamer_n;           ///< Original number of rotamers read from library.
  int current_rotamer_n;        ///< Current number of rotamers.  
                                ///< Rotamer number may increase or decrease because of diversity or bump check/energy cutoff  criteria.  

  vector<vector<string> > line_vv;   ///< Lines read in from rotamer library file.  A rotamer library contains many rotamers.  Each rotamer has its lines stored in vector<string>.  Hence, a vector of vector<string>.  Is this variable necessary?
  vector<Rotamer*> rot_v;       ///< Vector of pointers to rotamers.  .
  vector<Rotamer*>::const_iterator rot_itr;  ///< Iterator/pointer to current rotamer.

  /* Stuff relevent to connectivity info and anchor info.
   */

  RotConnInfo rotConnInfo;	// RotConnInfo structure.

  /**
  * Reader.  Private function that populates line_vv (vector<vector<string> > of read in lines.
  */
  
  void populate_line_vv(string in_Rotlib_File);  
  
  /**
   * Private function and helper functions that stores connectivity info as well as info like anchor info.  
   * These info are stored in rotConninfo.
   */
  void store_connectivity_info(string);
  void _populate_anchor_info(stringV); // pure readibilty function in store_connectivity_info.
  void _populate_atom_mapping_info(stringV); // same as above.
  void _populate_connectivity_info(stringV); // same as above.
    
};

/**
 *  AARotlib: AminoAcid rotamers.  Natural and Unnatural aminoacid rotamers are inherited from this class.
 */

class AARotlib : public Rotlib {
public:
  AARotlib();
  virtual ~AARotlib();

  AARotamer* get_next_rot();     ///< returns pointer to next Rotamer.  Returns NULL if reaches end of Rotlib.  Overrides base function.  Notice different return type.
  AARotamer* get_current_rot();  ///< returns pointer to current Rotamer.  Returns NULL if reaches end of Rotlib.  Overrides base function.  Notice different return type.
  AARotamer* get_rot(int);       ///< returns pointer to Rotamer with rotamer number specified.  REturns NULL if not found.
  AARotamer* reset_rot_pstn();   ///< resets internal rotamer pointer.  Returns pointer to first rotamer(in internal representation).
  AARotamer* set_rot_pstn(int);  ///< sets internal rotamer pointer to this rotamer number.  Returns NULL if not found.  Returns pointer to that rotamer if found.  If not found pointer would be reset.
  AARotamer* get_next_rot_with_empty_lattice_E_below(double); ///< returns next rotamer with an empty lattice E of below E provided as input.  E here is relative to base state, not absolute E.
  AARotamer* get_empty_lattice_E_rot(int); ///< returns pointer to Rotamer with Empty lattice emtpy ordering of int provided.

  void center_CA();                 ///< Translates coordinates of rotamers such that CA = (0,0,0).
  void calc_all_PHI();              ///< Calculates PHI angle for all rotamers in library.
  void calc_all_PSI();              ///< Calculates PSI angle for all rotamers in library.

  string resName;               ///< The type of AA this Rotlib contains.  3 letter abbreviation.

};

class NtrlAARotlib : public AARotlib {

public:


/**
 * Default constructor does nothing.
 */

  NtrlAARotlib();


  NtrlAARotlib(string);   ///< where string is rotamer library filename.  Rotamer connectivity file points to default location.
  NtrlAARotlib(string, string);   ///< like Rotlib, where the second specified rotamer connectivity file.  

  void setup_library();   ///< Initiates rotamer library.  Assigns charges, calculates PHI PSI angles, etc.

  //  void load_library(string);   ///< loads new library.

  /*  Below: loads rotamer library with specifications.  
   *  resName: three letter abbreviation.
   *  resolution: as usual, from 
   */
  //  NtrlAARotlib(string resName, int resolution, bool hydrogen);  

/**
 * Destructor for NtrlAARotlib.  Removes all objects pointed to from vector<Rotamer*> rot_v.
 *
 */

  virtual ~NtrlAARotlib();

  void assign_atom_fftype(); ///< assigns atom FFType for all rotamers in this rotamer library according to their atom labels.
  void assign_charges(string); ///< Assigns charges onto atoms according to charge assignment scheme.  Currently supported charge schemes: CHARM22.
  void assign_lone_pair();     ///< Assigns lone pair number and number of bonds that belong to each atom.
  
  void append_to_filehandle(ofstream*);       ///< Overrides fcn in Rotlib. Outputs rotamer library to ofstream*.

};

class Multiple_NtrlAARotlib : public NtrlAARotlib {

public:
  Multiple_NtrlAARotlib();
  Multiple_NtrlAARotlib(string, int); ///< string == path.  int == resolution.
  Multiple_NtrlAARotlib(string, int, vector<std::string>); ///< string == path, int == resolution, vector<string> == list of 1 letter amino acid libraries to include.
  ~Multiple_NtrlAARotlib();

  RotConnInfo* getRotConnInfo(); ///< Overloads base class getRotConnInfo.

  void add_library(std::string); ///< Appends rotamers from rotamer library to this rotlib.  No memory issues--all rotamer created and deleted by this rotlib.
  void add_library(NtrlAARotlib*);  ///< Appends rotamers from rotamer library to this rotlib.  Memory issues?

private:
  
  map<string, NtrlAARotlib*> rotlib_m;
  

};

class HIS_NtrlAARotlib : public NtrlAARotlib {

public:
  HIS_NtrlAARotlib();
  HIS_NtrlAARotlib(string , string, string, string); // the 2 singly protonated HIS libraries.  First string: H lib, second string, J lib.  Though order doesn't matter.
  ~HIS_NtrlAARotlib();

  void add_rotamer(Rotamer*); ///< Adds a rotamer, takes care of singly protonation states.

private:
  // use short cut: instantiate H_lib and J_lib and make them members.
  NtrlAARotlib H_AARotlib, J_AARotlib;


  void _make_other_singly_protonated_HIS(AARotamer*, AARotamer*); 
};

/* sample rotlib line:
          1         2         3         4         5
0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM       7  CB   SER      2    1.02194  -0.64038  -0.99924 C_3    4 0  0.05000 
*/


#endif /* ROTLIB_HPP */
