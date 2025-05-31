#include "scream_atom.hpp"
#include "scream_tools.hpp"

#include "sc_bgf_handler.hpp"
#include "defs.hpp"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include <algorithm>
using namespace std;

bgf_handler::bgf_handler() {

}

bgf_handler::bgf_handler(string filename) {

  this->readfile(filename);
  //  cout << "Address: in bgf_handler constructor: " << &atom_list << endl;
}

bgf_handler::bgf_handler(const bgf_handler& BH) {

  /* This copy constructor DEEP copies all attributes; but the most important attributes are coping SCREAM_ATOM's. */
  // First copy the string fields, easy.  Then again, these guys might no longer be relevant.

  this->header_lines = BH.header_lines;
  this->atom_lines = BH.atom_lines;
  this->conect_format_lines = BH.conect_format_lines;
  this->connectivity_record_lines = BH.connectivity_record_lines;

  // Then, deepcopy the ScreamAtomV field.
  // 1. First copy.  
  this->atom_list.clear();
  map<SCREAM_ATOM*, SCREAM_ATOM*> oldATOM_to_newATOM_map;

  for (ScreamAtomVConstItr itr = BH.atom_list.begin(); itr != BH.atom_list.end(); ++itr) {
    SCREAM_ATOM* oldAtom = *itr;
    SCREAM_ATOM* newAtom = new SCREAM_ATOM();

    newAtom->copy(*oldAtom);
    oldATOM_to_newATOM_map[oldAtom] = newAtom;

    this->atom_list.push_back(newAtom);
  
  }

  // 2. Then take care of connectivities.
  for (ScreamAtomVItr itr = this->atom_list.begin(); itr != this->atom_list.end(); ++itr) {
    SCREAM_ATOM* baseAtom = *itr;
    map<SCREAM_ATOM*, int> new_connectivity_m; new_connectivity_m.clear();
    for (map<SCREAM_ATOM*, int>::iterator mapItr = baseAtom->connectivity_m.begin(); mapItr != baseAtom->connectivity_m.end(); ++mapItr) {
      SCREAM_ATOM* oldConnectedAtom = mapItr->first;
      SCREAM_ATOM* newConnectedAtom = oldATOM_to_newATOM_map[oldConnectedAtom]; // guaranteed to find, though should be less lazy and put in an assert statement of sort.
      int order = mapItr->second;

      new_connectivity_m[newConnectedAtom] = order;
    }

    baseAtom->connectivity_m.clear();
    baseAtom->connectivity_m = new_connectivity_m;

  }

  // Done!
  cout << "new BGF_HANDLER constructed." << endl;

}


bgf_handler::~bgf_handler() {
  Debug info("bgf_handler::~bgf_handler");
  
  info.out("in bgf handler destructor");

  /* Iterate through all SCREAM_ATOM*'s in atom_list and delete them. */
  ScreamAtomVItr itr;
  if (atom_list.size() == 0) {
    info.out("atom_list was not allocated memory, no need to free up newed memory", "stdout");
  } else {
    for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
#ifdef DEBUG
      string s = (*itr)->return_bgf_line();
      info.out(s);
#endif
      delete *itr;
    }
    atom_list.clear();
  }

  info.out("Bgf_handler Destruction successful!");

}

bool bgf_handler::readfile(const string filename) {
  /* This is the main readfile in bgf_handler class. */

  if (this->atom_list.size() != 0) {
    cout << "Object already initialized!" << endl;
    cout << "To open new file, create a new sc_bgf_handler object." << endl;
    return 0;
  }

  ifstream FILE(filename.c_str());
  /* returns 0 if unable to open file */
  if ( !(FILE.good() ) ) {
    cerr << "Unable to open " << filename << endl;
    exit(2);
    return 0;
  }
  string line;
  char line_ch[256];
  while (!FILE.eof()) {

    FILE.getline(line_ch, sizeof(line_ch));
    line = string(line_ch);
    /* Conditions; populate xxxlines. */

    if (line.substr(0,3) == "END") {
        break;
    }
    if (scream_tools::is_bgf_header_line(line) ) {
      this->header_lines.push_back(line);
    } 
    else if (scream_tools::is_bgf_atom_line(line) ) {
      this->atom_lines.push_back(line);
    }
    else if (scream_tools::is_format_connect_line(line) ) {
      this->conect_format_lines.push_back(line);
    }
    else if (scream_tools::is_connectivity_info(line )) {
      this->connectivity_record_lines.push_back(line);
      }
  }
  
  FILE.close();

  /* Now populate atomlist */

  for (stringVConstItr itr = atom_lines.begin(); itr != atom_lines.end(); ++itr) {
    SCREAM_ATOM* atom = new SCREAM_ATOM(*itr);
    this->atom_list.push_back(atom);
  }

  cout << "Number of atoms in system: " << atom_list.size() << endl;

  /* Now make bonds. */
  this->make_bonds(connectivity_record_lines);

}


bool bgf_handler::readPDB(const string filename) {
  /* This is the readfile for PDB in bgf_handler class. */

  if (this->atom_list.size() != 0) {
    cout << "Object already initialized!" << endl;
    cout << "To open new file, create a new sc_bgf_handler object." << endl;
    return 0;
  }

  ifstream FILE(filename.c_str());
  /* returns 0 if unable to open file */
  if ( !(FILE.good() ) ) {
    cerr << "Unable to open " << filename << endl;
    exit(2);
    return 0;
  }
  string line;

  while (!FILE.eof()) {
    char line_ch[256];
    FILE.getline(line_ch, sizeof(line_ch));
    line = string(line_ch);
    /* Conditions; populate xxxlines. */

    if (line.substr(0,3) == "END") {
        break;
    }
    if (scream_tools::is_bgf_header_line(line) ) { 
        this->header_lines.push_back(line);
    } 
    else if (scream_tools::is_bgf_atom_line(line) ) { // same as bgf atom keyword
      this->atom_lines.push_back(line);
    }
    else if (scream_tools::is_format_connect_line(line) ) {
      this->conect_format_lines.push_back(line);
    }
    else if (scream_tools::is_connectivity_info(line )) { // same as bgf conn line
      this->connectivity_record_lines.push_back(line);
      }
  }
  
  /* Now populate atomlist */

  for (stringVConstItr itr = atom_lines.begin(); itr != atom_lines.end(); ++itr) {
    SCREAM_ATOM* atom = new SCREAM_ATOM();
    atom->pdb_init(*itr);
    this->atom_list.push_back(atom);
  }
  
  cout << "Number of atoms in system: " << atom_list.size() << endl;

  /* Now make bonds. */
  this->make_pdb_bonds(connectivity_record_lines); 

}


bool bgf_handler::readPDB(const string filename, ScreamAtomV& atomList) {
  // do nothing for now
}

bool bgf_handler::readfile(const string filename, ScreamAtomV& atom_list) {
  /* Calls readfile(const string filename).  Used when a ScreamAtomV is passed in to be populated. */

  this->readfile(filename);
  this->pass_atomlist(&atom_list);  // this is buggy.  &atom_list is not an lvalue, not assignable.  will debug this later--this is a run time error, will not caught by compiler.

}


bool bgf_handler::printToFile(std::string filename) {

  ostream* OUTPUT = new ofstream;
  if (filename == "stdout") {
    OUTPUT = &cout;
  } else {
    ((ofstream*)OUTPUT)->open(filename.c_str());
  }
  
  this->printfile(this->atom_list, OUTPUT, "", 332);
  if (filename != "stdout") {

    ((ofstream*)OUTPUT)->close();
  }

  
  delete OUTPUT;

  return 1;

}

bool bgf_handler::printToFile(std::string filename, std::string additionalRemark) {

  ostream* OUTPUT = new ofstream;
  if (filename == "stdout") {
    OUTPUT = &cout;
  } else {
    ((ofstream*)OUTPUT)->open(filename.c_str());
  }
  
  this->printfile(this->atom_list, OUTPUT, additionalRemark, 332);
  
  if (filename != "stdout") {

    ((ofstream*)OUTPUT)->close();
  }

  delete OUTPUT;
  return 1;

}




bool bgf_handler::printfile(const ScreamAtomV& atom_list, ostream* ostream_p, string additionalREMARK, int BIOGRF) {
  /* Prints an atom_list */
  //  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
  //    (*itr)->dump();
  //  }

  /** HEADER lines **/
  
  /* DESCRP BGF VERSION LINE, first find it. */
  string DESCRP_field = "";
  for (stringVItr itr = this->header_lines.begin(); itr != this->header_lines.end(); ++itr) {
    if ( (*itr).substr(0,6) == "DESCRP" ) {
      vector<string> fields = scream_tools::split(*itr);
      if (fields.size() >= 2) {
	for (int i = 1; i < fields.size(); ++i) {
	  DESCRP_field = DESCRP_field + fields[i];
	  DESCRP_field += " ";
	}
      }
      break;
    }
  }

  if (DESCRP_field == "") {
    DESCRP_field = "N/A";
  }

  *ostream_p << string("BIOGRF ") << BIOGRF << endl;
  /* REMARK LINES--THESE ARE PASSED IN */ 
  *ostream_p << string("DESCRP ") << DESCRP_field << endl;
  *ostream_p << string("REMARK") << additionalREMARK << endl;

  /* FORCEFIELD and ATOM FORMAT LINE */
  *ostream_p << string("FORCEFIELD DREIDING") << endl;
  *ostream_p << string("FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)") << endl;


  /* ATOM and HETATM lines */
  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    //cout << *itr << endl;
    //cout << (*itr)->resName << endl;
    //cout << (*itr)->n << endl;
    (*itr)->append_to_filehandle(ostream_p);
    //cout << "before dump" << endl;
    //(*itr)->dump();
  }
  /* FORMAT CONECT lines */
  *ostream_p << string("FORMAT CONECT (a6,14i6)") << endl;
  *ostream_p << string("FORMAT ORDER (a6,i6,13f6.3)") << endl;

  /* CONECT i.e. connectivity lines */
  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    (*itr)->append_to_ostream_connect_info(ostream_p);
  }

  /* END */
  *ostream_p << "END" << endl;

  return 1;

}

bool bgf_handler::printPDB(std::string filename) {

  ostream* OUTPUT = new ofstream;
  if (filename == "stdout") {
    OUTPUT = &cout;
  } else {
    ((ofstream*)OUTPUT)->open(filename.c_str());
  }
  
  this->printToPDB(this->atom_list, OUTPUT, "");
  
  if (filename != "stdout") {
    ((ofstream*)OUTPUT)->close();
  }

  delete OUTPUT;
  return 1;

}

bool bgf_handler::printToPDB(const ScreamAtomV& atom_list, ostream* ostream_p, string HEADER) {
  /** HEADER lines **/
  
  /* DESCRP PDB HEADER LINE */
  *ostream_p << string("HEADER    ") << HEADER << endl;
  /* REMARK LINES--THESE ARE PASSED IN */ 
  *ostream_p << string("TITLE") << endl;

  /* ATOM and HETATM lines */
  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    (*itr)->pdb_append_to_filehandle(ostream_p);
  }

  /* CONECT i.e. connectivity lines */
  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    (*itr)->pdb_append_to_ostream_connect_info(ostream_p);
  }

  /* END */
  *ostream_p << "END" << endl;

  return 1;

}

bool bgf_handler::printSequenceToFile(std::string filename) {
  /* Format: simply a one-liner with all the amino acids.  Amino acids in chains A will appear before those in chain B.  Those that are unnatural, or RES, will be outputted as X.  Water's will not be outputted. */
  /* Calls printSequenceToFile */

  ostream* OUTPUT = new ofstream;
  if (filename == "stdout") {
    OUTPUT = &cout;
  } else {
    ((ofstream*)OUTPUT)->open(filename.c_str());
  }
  
  this->printSequence(this->atom_list, OUTPUT);
  
  if (filename != "stdout") {
    ((ofstream*)OUTPUT)->close();
  }

  delete OUTPUT;
  return 1;

}

std::string bgf_handler::returnSequence() {
  /* Format: simply a one-liner with all the amino acids.  Amino acids will appear in the order they are encountered from the bgf file that was originally read in, regardless of chain name.  Those that are unnatural, or RES, will be outputted as X.  Water's will not be outputted. */
  string returnString = "";
  /* ATOM and HETATM lines */
  int resNum = 0;
  for (ScreamAtomVConstItr itr = this->atom_list.begin(); itr != this->atom_list.end(); ++itr) {
    int crntResNum = (*itr)->getResNum();
    if (crntResNum != resNum) {
      string crntResName = (*itr)->getResName();
      string oneLetter = scream_tools::one_letter_AA(crntResName);
      returnString += oneLetter;
      resNum = crntResNum;
    }
  }
  return returnString;

}

bool bgf_handler::printSequence(const ScreamAtomV& atom_list, ostream* ostream_p) {
  /* Format: simply a one-liner with all the amino acids.  Amino acids will appear in the order they are encountered from the bgf file that was originally read in, regardless of chain name.  Those that are unnatural, or RES, will be outputted as X.  Water's will not be outputted. */

  /* ATOM and HETATM lines */
  int resNum = 0;
  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    int crntResNum = (*itr)->getResNum();
    if (crntResNum != resNum) {
      string crntResName = (*itr)->getResName();
      string oneLetter = scream_tools::one_letter_AA(crntResName);
      *ostream_p << oneLetter;
      resNum = crntResNum;
    }
  }
  *ostream_p << endl;
  /* CONECT i.e. connectivity lines */
  // ignore.

  /* END */
  return 1;

}

void bgf_handler::pass_atomlist(ScreamAtomV* atom_list) {
  /* Reassigns the ScreamAtomV list with contents in this->atom_list. */

  //  atom_list.clear();
  //scream_tools::deep_copy_ScreamAtomV(this->atom_list, atom_list);
  atom_list = &(this->atom_list);   // this is just plain wrong.  Have to pass pointer to pointer.



  //  cout << "Size of this->atom_list: " << this->atom_list.size() << endl;
  //  cout << "Size of passed in atom_listt: " << (*atom_list).size() << endl;
  //  atom_list.assign(this->atom_list.begin(), this->atom_list.end());

}

/* Private member functions */

void bgf_handler::make_bonds(stringV& connectlines) {
  /* Makes a map/dictionary for easy access from atom number to atom structure.  Just like get_atom() function in ModulaSim. */
  /* Note: this is migrated from Protein class */

  map<int, SCREAM_ATOM*> atom_map;  // For easy access to SCREAM_ATOM* and its atom number.
  for (ScreamAtomVConstItr itr = this->atom_list.begin(); itr != atom_list.end(); ++itr) {
    atom_map.insert(make_pair( (*itr)->n, *itr ) );
  }

  /* Looping over CONNECTIVITY LINES*/
  for (stringVConstItr itr = connectivity_record_lines.begin(); itr != connectivity_record_lines.end(); ++itr) {
    
    bool base_atom_FLAG = true;   // If CONECT 34 51 52, base atom is 34, while 51 and 52 are connectivities to be added to base atom.
    SCREAM_ATOM* base_atom = NULL;
    SCREAM_ATOM* atom = NULL;     // Placeholder.

    stringV fields = scream_tools::split(*itr);
    map<int, SCREAM_ATOM*>::iterator atom_itr;    
    if ( *(fields.begin() ) == string("ORDER") )  {   // Skip ORDER lines for now.  
      continue;
    } else {  // else make connectivity

      for (stringVConstItr itr_f = fields.begin(); itr_f != fields.end(); ++itr_f) {
	if ( *itr_f == string("CONECT") ) 
	  continue;                  // if first entry is CONECT, skip.
	
	int atom_n = atoi( (*itr_f).c_str());
	atom_itr = atom_map.find(atom_n);
	if (atom_itr == atom_map.end() ) {
	  cout << "Can't find atom with atom number " << atom_n << "!  Please check your connectivities before proceeding.  Exiting.";
	  fflush(stdout);
	  exit(3);
	}
	if ( base_atom_FLAG == true) {  // FLAG to find base atom in a CONECT line
	  base_atom = atom_itr->second;
	  base_atom_FLAG = false;
	  continue;                      // keep looping.

	} else { 

	  atom = atom_itr->second;
	  base_atom->make_bond(atom);

	}

      }  /* end loop one line of CONECTIVITY record */

    } /* end else ORDER */
    
  }  /* end of loop CONECTIVITY LINES */
  
}

void bgf_handler::make_pdb_bonds(stringV& connectlines) {
  /* Makes a map/dictionary for easy access from atom number to atom structure.  Just like get_atom() function in ModulaSim. */
  /* Note: this is migrated from Protein class */

  Debug info("bgf_handler::make_pdb_bonds(stringV& connectlines)");
  map<int, SCREAM_ATOM*> atom_map;  // For easy access to SCREAM_ATOM* and its atom number.
  for (ScreamAtomVConstItr itr = this->atom_list.begin(); itr != atom_list.end(); ++itr) {
    atom_map.insert(make_pair( (*itr)->n, *itr ) );
  }

  /* Looping over CONNECTIVITY LINES*/
  for (stringVConstItr itr = connectivity_record_lines.begin(); itr != connectivity_record_lines.end(); ++itr) {
    
    bool base_atom_FLAG = true;   // If CONECT 34 51 52, base atom is 34, 51, 52 are connectivities to be added to base atom.
    SCREAM_ATOM* base_atom = NULL;
    SCREAM_ATOM* atom = NULL;     // Placeholder.

    //stringV fields = scream_tools::split(*itr); inpdb files, need to split by every 5 atom n.
    // like :Conect 5051 5052 5053 5054 5073
    stringV fields;
    string l = *itr;
    cout << l << endl;
    int i = 0;
    fields.push_back("CONECT");
    while (6+5*i+5 < l.length()) {
      string s = l.substr(6+5*i,5);
      info.out(s);
      if (s == "     ") { // if a line is like: CONECT 5051 5052 5053 5054 5073                                         1ADS2736
	break;
      }
      fields.push_back(s);
      i++;
    }
    
    if ( *(fields.begin() ) == string("ORDER") )  {   // Skip ORDER lines for now.  
      continue;
    } else {  // else make connectivity

      for (stringVConstItr itr_f = fields.begin(); itr_f != fields.end(); ++itr_f) {

	if ( *itr_f == string("CONECT") ) {
	  continue;                  // if first entry is CONECT, skip.
	} 

	// DEBUG LINE	cout << (*itr_f) << endl;

	int atom_n = atoi( (*itr_f).c_str());


	if ( base_atom_FLAG == true) {  // FLAG to find base atom in a CONECT line

	  base_atom = atom_map.find(atom_n)->second;
	  base_atom_FLAG = false;
	  continue;                      // keep looping.

	} else { 

	  atom = atom_map.find(atom_n)->second;
	  base_atom->make_bond(atom);

	}

      }  /* end loop one line of CONECTIVITY record */

    } /* end else ORDER */
    
  }  /* end of loop CONECTIVITY LINES */
  
}

