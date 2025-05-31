#ifndef BGF_READER_HPP
#define BGF_READER_HPP

#include "scream_atom.hpp"

#include "defs.hpp"

#include <iostream>
#include <string>
using namespace std;

/** bgf_handler:
 *  handles I/O aspest of bgf file.  
 *  New: now also handles PDB files!
 *  REMARK: printing anything on the protein level should be done by creating a bgf_handler object..
 */

class bgf_handler {
public:

  bgf_handler();
  bgf_handler(string);	///< Creates a bgf_handler object.  String is name of bgf file.
  bgf_handler(const bgf_handler&); ///< Copy constructor.  , in the sense that all SCREAM_ATOM*'s are contructued a second time.
  ~bgf_handler();		///< Destructor.  Destroys ScreamAtomV if that list is initialized.


  bool readfile(const string);                ///< Reads file.

  bool readfile(const string, ScreamAtomV&);  ///< Reads file and populates an atom list passed in by reference.

  bool readPDB(const string);	///< Reads PDB file.
  
  bool readPDB(const string, ScreamAtomV&); ///< Reads PDB file and populates an atom list passed in by reference.

  // Print commands.

  bool printToFile(std::string);                     ///< Prints a bgf file from atom_list stored in this class.
  bool printToFile(std::string, std::string);        ///< Second: additional remark.
  bool printfile(const ScreamAtomV&, ostream*, string = string(""), int = 332); ///< Prints a bgf file, taking a list of SCREAM ATOMS.
  
  bool printPDB(std::string);  ///< Prints structure in pdb format from atom_list 
  bool printToPDB(const ScreamAtomV&, ostream*, string = string("") ); ///< Prints a file in pdb format, taking a list of SCREAM ATOMS.

  // Print sequence commands.

  bool printSequenceToFile(std::string); ///< Prints a file that contains all the sequence of this bgf file amino acid.  Format: just a string of one-letter AA's.
  std::string returnSequence(); ///< Returns a string containing the sequence info.
  static bool printSequence(const ScreamAtomV&, ostream*);

  void pass_atomlist(ScreamAtomV*); ///< Populates a ScreamAtomV atomlist passed in by reference using the ScreamAtomV stored in this class.
  ScreamAtomV* getAtomList() { return &(this->atom_list); }; ///< Returns a pointer to the underlying atom list; so that modification is passible outside of this class.

  /* Variables */
  ScreamAtomV atom_list;

  /*! BGF header info */
  stringV header_lines;
  
  /*! BGF atom lines */
  stringV atom_lines;

  /*! BGF format conect lines */
  stringV conect_format_lines;

  /*! BGF connectivity lines */
  stringV connectivity_record_lines;
  

private: 

  void make_bonds(stringV&);
  void make_pdb_bonds(stringV&);

};

#endif /* BGF_READER_HPP */
