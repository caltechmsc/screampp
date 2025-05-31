/* Rotlib.cpp
 *
 * Source for classes relevant to rotamer library in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>

using std::string;
using namespace std;

#include "defs.hpp"
#include "scream_atom.hpp"
#include "Rotamer.hpp" 
#include "AARotamer.hpp"
#include "RotamerCluster.hpp"
#include "sc_BackBone.hpp"
#include "sc_AABackBone.hpp"
#include "sc_SideChain.hpp"
#include "sc_AASideChain.hpp"


#include "scream_tools.hpp"
#include "Rotlib.hpp"

#include "scream_vector.hpp"

Rotlib::Rotlib() {

  this->library_name = "";
  this->orig_rotamer_n = 0;
  this->rot_itr = this->rot_v.begin();  
  
  //  this->rot_itr = NULL;  not allowed by gcc 3.0 or above
}

Rotlib::Rotlib(std::string connectivity_file) {

  this->library_name = "";
  this->readConnectivityFile(connectivity_file);
  cout << "Done with reading ConnectivityFile.  Initializing rot iterators." << endl;
  this->rot_itr = this->rot_v.begin();  
}

Rotlib::Rotlib(string coordinate_file, string connectivity_file) {

  this->library_name = "";
  this->populate_line_vv(coordinate_file);
  this->store_connectivity_info(connectivity_file);

  int n = 1;
  /* Now create Rotamers* (really Conformers*). */
  vector<vector<string> >::const_iterator itr_vv = this->line_vv.begin();
  for (; itr_vv != this->line_vv.end(); ++itr_vv) {
    //cout << "Rotamer count: " << n << endl;
    //cout << "before newing Rotamer" << endl;
    Rotamer* new_r = new Rotamer(*itr_vv, &(this->rotConnInfo));
    //new_r->library_name = this->library_name; // vk, 01/31/05
    //cout << "after newing Rotamer" << endl;
    this->rot_v.push_back(new_r);
    new_r->setDeclaredInRotlibScope(true);
    ++n;
  }
  /* Initialize rot_itr. */
  rot_itr = rot_v.begin();
}

Rotlib::~Rotlib() {
  //cout << "Deleting Rotlib! " << endl;
  vector<Rotamer*>::const_iterator itr;
  //cout << "rot_v size: " << rot_v.size()  << endl;
  if (rot_v.size() != 0) 
    {
      for (itr = rot_v.begin(); itr != rot_v.end(); ++itr) {
	if ( (*itr)->declaredInRotlibScope() ) {
	  delete *itr;
	}
      }
    }

  //cout << "Done deleting Rotlib!" << endl;
} 


void Rotlib::readConnectivityFile(string connectivity_file) {
  /* If this is done, delete whatever is currerntly stored in this->rot_v. */
  /* Also, note that this routine also reads in the coordinate file. */

  cout << "Deleting Rotamers currently populated in this Rotamer library object..." << endl;
  vector<Rotamer*>::const_iterator itr;
  //cout << "rot_v size: " << rot_v.size()  << endl;
  if (rot_v.size() != 0) 
    {
      for (itr = rot_v.begin(); itr != rot_v.end(); ++itr) {
	if ( (*itr)->declaredInRotlibScope() ) {
	  cout << "in cleaning rotamers, declaredInRotlibScope true" << endl;
	  delete *itr;
	}
      }
    }
  this->rot_v.clear();
  this->line_vv.clear();
  this->rotConnInfo.clear();
  
  //  cout << "Breakpoint 1" << endl;

  this->store_connectivity_info(connectivity_file);

  //  cout << "breakpoint 1.5" << endl;

  /* find string that coordinate_file string in connectivity_file */
  string coordinate_file = this->rotConnInfo.targetRotamerLibFile;
  this->populate_line_vv(coordinate_file);

  //  cout << "Breakpoint 2 " << endl;

  /* Now create Rotamers* (really Conformers*). */
  int n = 1;
  vector<vector<string> >::const_iterator itr_vv = this->line_vv.begin();
  for (; itr_vv != this->line_vv.end(); ++itr_vv) {
    cout << "Rotamer count: " << n << endl;
    Rotamer* new_r = new Rotamer(*itr_vv, &(this->rotConnInfo));

    this->rot_v.push_back(new_r);
    new_r->setDeclaredInRotlibScope(true);
    ++n;
  }

  //  cout << "breakpoint 3" << endl;

  /* Initialize rot_itr. */
  rot_itr = rot_v.begin();

  cout << "readConnectivityFile done." << endl;

}




void Rotlib::readRotamerLibrary(string file) {
  cout << "I'm a Rotlib.  I don't directly read a rotamer library.  I read a connectivity file.  To read a rotamer library, you want NtrlAARotlib.  They assume something about the connecitivities of the rotamers, I don't.  I am more general." << endl;

}

Rotamer* Rotlib::get_next_rot() {
  Rotamer* to_be_returned = NULL;
  if (rot_itr == this->rot_v.end() ) {
    //return NULL;  // do nothing; to_be_returned is already NULL>
  } else {
    to_be_returned = *rot_itr;
    ++rot_itr;
  }
  return to_be_returned;

}

Rotamer* Rotlib::get_current_rot() {
  if (rot_itr != this->rot_v.end() ) {
    return *rot_itr;
  } else {
    return NULL;
  }
  
}

void Rotlib::reset_pstn() {

  rot_itr = rot_v.begin();

}

Rotamer* Rotlib::get_next_rot_with_empty_lattice_E_below(double CUTOFF_E) {

  assert(CUTOFF_E >= 0);
  Rotamer* to_be_returned;

  if (rot_itr == rot_v.begin() ) { // head of vector
    to_be_returned = (*rot_itr);
    ++rot_itr;
    return to_be_returned;
    
  } 
  else {			// any other position

    if (rot_itr != rot_v.end() ) {
      to_be_returned =  (*rot_itr) ;
      ++rot_itr;
      // if less than CUTOFF_E, return
      if ( to_be_returned->get_empty_lattice_E() <= CUTOFF_E) {
	return to_be_returned;
      } 
      // else, return NULL
      else {
	return NULL;
      }
    } 
    else {			// if at end of vector
      return NULL;
    }
    //    ++rot_itr;			// should never reach here
    
  }


}


Rotamer* Rotlib::get_empty_lattice_E_rot(int EL_rank) {

  vector<Rotamer*>::const_iterator temp_RotItr = this->rot_v.begin();

  Rotamer* to_be_returned;
  for (; temp_RotItr != rot_v.end(); ++temp_RotItr) {
    int rank = (*temp_RotItr)->get_empty_lattice_energy_rank();

    if ( (*temp_RotItr)->get_empty_lattice_energy_rank() == EL_rank) {
      to_be_returned = (*temp_RotItr);
      break;      
    }
  }
    
  if (rot_itr == rot_v.end()) {
    to_be_returned = NULL;
  }

  return to_be_returned;
}


Rotamer* Rotlib::get_empty_lattice_E_rot_after_sorted_by_empty_lattice_E(int EL_rank) {
  /* Assumes that rot_v is sorted by empty lattice energies.  Returns NULL if EL_rank is not found. */
  Rotamer* to_be_returned;
  if (EL_rank >= this->rot_v.size() ) { // >= because: if vector size 10, index is from 0 to 9.
    to_be_returned = NULL;    
  }
  else {
    to_be_returned = this->rot_v[EL_rank]; // EL_rank starts from 0, already sorted.
  }

  return to_be_returned;
}




void Rotlib::print_Me() {

  vector<Rotamer*>::const_iterator itr;

  for (itr = this->rot_v.begin(); itr != this->rot_v.end(); ++itr) {
    
    //(*itr)->print_ordered_by_n();
    (*itr)->print_Me();

  }

}

void Rotlib::add_rotamer(Rotamer* rot) {

  this->rot_v.push_back(rot);
  if (rot->library_name == "") {
    rot->library_name = this->library_name;
  }

}

Rotamer* Rotlib::new_rotamer() {

  Rotamer* new_rot = new Rotamer();
  new_rot->setDeclaredInRotlibScope(true);
  this->rot_v.push_back(new_rot);
  return new_rot;

}

RotamerCluster* Rotlib::new_rotamer_cluster() {
  RotamerCluster* new_rot = new RotamerCluster();
  new_rot->setDeclaredInRotlibScope(true);
  this->rot_v.push_back(new_rot);
  return new_rot;

}

int Rotlib::size() {

  return this->rot_v.size();

}


class cmp_empty_lattice_E_abs {
public:
  bool operator()(Rotamer* r1, Rotamer* r2) {
    return ( r1->get_empty_lattice_E_abs() < r2->get_empty_lattice_E_abs());
  }
};


int Rotlib::n_rotamers_below_empty_lattice_energy(double E) {

  assert(E > 0);
  int c = 0;
  vector<Rotamer*> tmp_rot_v = this->rot_v;
  sort(tmp_rot_v.begin(), tmp_rot_v.end(), cmp_empty_lattice_E_abs());

  for (vector<Rotamer*>::const_iterator itr = tmp_rot_v.begin(); 
       itr != tmp_rot_v.end(); ++itr) {

    if ( (*itr)->get_empty_lattice_E() > E) {
      break;
    } 
    else {
      ++c;
    }

  }
  return c;

}

class cmp_rotlib_E {
public:
  bool operator()(Rotamer* r1, Rotamer* r2) {
    return ( r1->get_rotlib_E() < r2->get_rotlib_E() );
  }
};

void Rotlib::sort_by_rotlib_E() {
  
  //  sort(rot_v.begin(), rot_v.end(), Rotlib::cmp_rotlib_E);
  sort(rot_v.begin(), rot_v.end(), cmp_rotlib_E());
  //  sort(rot_v.begin(), rot_v.end(), &(Rotlib::cmp_rotlib_E));
  int i = 1;
  for (vector<Rotamer*>::iterator itr = rot_v.begin(); itr != rot_v.end(); ++itr, ++i) {
    (*itr)->set_rotamer_n(i);
  }

}

class cmp_self_E {
public:
  bool operator()(Rotamer* r1, Rotamer* r2) {
    return ( r1->get_sc_total_E() < r2->get_sc_total_E() );
  }
};


void Rotlib::sort_by_self_E() {

  sort(rot_v.begin(), rot_v.end(), cmp_self_E());
  int i = 1;
  for (vector<Rotamer*>::iterator itr = rot_v.begin(); itr != rot_v.end(); ++itr, ++i) {
    (*itr)->set_rotamer_n(i);
  }

}

void Rotlib::sort_by_empty_lattice_E() {
  
  sort(rot_v.begin(), rot_v.end(), cmp_empty_lattice_E_abs());
  int i = 0;
  double lowestEnergy = (*(rot_v.begin()))->get_empty_lattice_E_abs(); // Rotamers have already been sorted, so head of list has lowest E.
  for (vector<Rotamer*>::iterator itr = rot_v.begin(); itr != rot_v.end(); ++itr, ++i) {

    double thisRotamerAbsRotamerE = (*itr)->get_empty_lattice_E_abs();
    (*itr)->set_empty_lattice_E(thisRotamerAbsRotamerE - lowestEnergy);
    (*itr)->set_empty_lattice_energy_rank(i);      
    //    (*itr)->set_rotamer_n(i);
  }
}

double Rotlib::get_best_preCalc_E() {
  double best_preCalc = 99999;
  for (vector<Rotamer*>::const_iterator i = this->rot_v.begin(); i != this->rot_v.end(); ++i) {
    double preCalc = (*i)->get_preCalc_TotE();
    if (preCalc < best_preCalc ) {
      best_preCalc = preCalc;
    }
  }
  return best_preCalc;
}


void Rotlib::populate_line_vv(string in_Rotlib_File) {

  vector<string> line_tmp;          // Temporary line storage container.
  
  std::ifstream IN_ROTLIB(in_Rotlib_File.c_str());

  //  cout << in_Rotlib_File << endl;

  if (!IN_ROTLIB.is_open()) {
    std::cerr << "Unable to open 1 " << in_Rotlib_File << endl;
    exit(8);
  }

  if (IN_ROTLIB.bad()) {
    std::cerr << "Unable to open 2 " << in_Rotlib_File << endl;
    exit (8);
  }
  
  // read in rotamer library and count number of rotamers
  bool REM_FLAG = false;
  bool first_FLAG = true;
  this->orig_rotamer_n = 0;
  //  int c = 1;

  while (true) {
    char line_c[25600];
    IN_ROTLIB.getline(line_c, sizeof(line_c));
    string line = string(line_c);
    //cout << "in Rotlib.cpp populate_line_vv" << line << endl;
    if (line.substr(0,1) == "#") {
      continue; // # comments are headed by #.
    }
    
    if (line.substr(0,3) == "LIB") {
      vector<string> fields = scream_tools::split(line);
      this->library_name = fields[1];
      continue;
    }
    
    if (line.substr(0, 3) == "REM" && REM_FLAG != true) {
      this->orig_rotamer_n++;
      REM_FLAG = true;

      // Initial conditions.
      if (first_FLAG == true) {
	first_FLAG = false;
      } 
      else if (first_FLAG == false) {
	vector<string> one_rot_tmp_v(line_tmp);
	this->line_vv.push_back(one_rot_tmp_v);
	line_tmp.clear();
      }

    } else if (line.substr(0,4) == "ATOM") {
      REM_FLAG = false;
    }

    if (line != "\n") {
      line_tmp.push_back(line);
    }
    // End conditions.

    if (IN_ROTLIB.eof()) {
      vector<string> one_rot_tmp_v(line_tmp);
      this->line_vv.push_back(one_rot_tmp_v);
      //      vector<string>* one_rot_tmp_v = new vector<string>(line_tmp);
      //      this->line_vv.push_back(*one_rot_tmp_v);
      break;
    }
  }
  IN_ROTLIB.close();

}

void Rotlib::store_connectivity_info(string connectivity_file) {

  //  cout << "STORE_CONNECTIVITY_INFO is : " << connectivity_file << endl;

  stringV anchor_info_lines;	// aptly named: a vector of strings that contain the anchor info lines.
  stringV atom_mapping_lines; 	// as above, but for atom mapping.
  stringV connectivity_info_lines; // as above, but for connectivity info.

  std::ifstream IN_FILE(connectivity_file.c_str());
  
  if ( (!IN_FILE.is_open())
       or IN_FILE.bad() ) {
    std::cerr << "Unable to open 3 " << connectivity_file << endl;
    std::cerr << "Failure in Rotlib::store_connectivity_info.  Quitting." << endl;
    exit(8);
  }

  /* Parsing File and Splitting lines_t into 3 parts: Anchor info, ATOM_MAPPING info and CONECT info. */
  int read_mode_flag = 0;	// 0: init.  1: reading ANCHOR_INFO.  2: reading ATOM_MAPPING.  3: reading FORMAT CONECT.
  int i = 0;

  while (!IN_FILE.eof()) {
    char line_c[25600];
    IN_FILE.getline(line_c, sizeof(line_c));
    if ( (line_c[0] == '#') or
	 (line_c[0] == '/') or
	 (line_c[0] == '%') or
	 (line_c[0] == ' ') or
	 (line_c[0] == '\n') or
	 (line_c[0] == '\t') or
	 (line_c[0] == '\0') ) {
      continue;
    }
    string line_str(line_c);
    stringV fields;
    split(line_str, string(" "), fields);

    if (fields[0].substr(0,3) == string("END")) {
      read_mode_flag = 0;
      continue;
    } else if (fields[0] == string("ANCHOR_INFO")) {
      read_mode_flag = 1;
      continue;
    } else if (fields[0] == string("ATOM_MAPPING")) {
      read_mode_flag = 2;
      continue;
    } else if (fields[0] == string("FORMAT")) {	// skips all lines that start with FORMAT; i.e. FORMAT CONECT and FORMAT ORDER
      read_mode_flag = 3;
      continue;
    }
    
    if (read_mode_flag == 0) {
      continue;
    } else if (read_mode_flag == 1) {
      anchor_info_lines.push_back(line_str);
    } else if (read_mode_flag == 2) {
      atom_mapping_lines.push_back(line_str);
    } else if (read_mode_flag == 3) {
      connectivity_info_lines.push_back(line_str);
    }
  }

//   cout << "ANCHOR_INFO lines: " << endl;
//   for_each (anchor_info_lines.begin(), anchor_info_lines.end(), print);
//   cout << "ATOM_MAPPING lines: " << endl;
//   for_each (atom_mapping_lines.begin(), atom_mapping_lines.end(), print) ;
//   cout << "CONNECTIVITY_INFO lines: " << endl;
//   for_each (connectivity_info_lines.begin(), connectivity_info_lines.end(), print);

  this->_populate_anchor_info(anchor_info_lines);
  //  cout << "passed this" << endl;
  this->_populate_atom_mapping_info(atom_mapping_lines);
  //  cout << "passed atom mapping" << endl;
  this->_populate_connectivity_info(connectivity_info_lines);
  //  cout << "passed populate connectivity_info_lines" << endl;

  //  cout << "connectivity_info.size: " << this->rotConnInfo.atom_connectivity_info.size() << endl;
  //  cout << "atom_n_map.szie: " << this->rotConnInfo.atom_n_map.size() << endl;

  if (this->rotConnInfo.atom_connectivity_info.size() > this->rotConnInfo.atom_n_map.size()) {
    cout << "Atom mapping size if greater than atom connectivity size.  Quitting." << endl;
    exit(8);
  }


}

void Rotlib::_populate_anchor_info(stringV anchor_info_lines) {
  /* populates TargetRotamerLibFile, anchor_pts, atom_n_map, atom_n_label_map, side_chain_atoms */
  for (stringVConstItr itr = anchor_info_lines.begin(); itr != anchor_info_lines.end(); ++itr) {
    stringV fields;
    split(*itr, fields);
    string keyw = fields[0];
    /*cout << keyw << endl;
    for (stringVConstItr itr = fields.begin(); itr != fields.end(); ++itr) {
      cout << *itr << ":::";
    }
    cout << endl;
    */

    if (keyw == "TargetRotamerLibFile") {
      this->rotConnInfo.targetRotamerLibFile = fields[1];
    } 
    else if (keyw == "AnchorAtoms") {
      this->rotConnInfo.anchor_pts.clear();
      for (int i = 1; i < fields.size(); ++i) {
	int anchorAtomN = atoi(fields[i].c_str());
	this->rotConnInfo.anchor_pts.push_back(anchorAtomN);
      }
    }
    else if (keyw == "SideChainAtoms") {
      this->rotConnInfo.side_chain_atoms.clear();
      for (int i = 1; i < fields.size(); ++i) {
	int sideChainAtomN = atoi(fields[i].c_str());
	this->rotConnInfo.side_chain_atoms.push_back(sideChainAtomN);
      }
    }
    else if (keyw == "RotamerAxis") {
      // do nothing
    }
    else if (keyw == "MatchChirality") {
      // do nothing
    }
    else if (keyw == "AtomsOfExactMatch") {
      this->rotConnInfo.atoms_of_exact_match.clear();
      for (int i = 1; i < fields.size(); ++i) {
	int exactMatchAtomN = atoi(fields[i].c_str());
	this->rotConnInfo.atoms_of_exact_match.push_back(exactMatchAtomN);
      }
    }
    else if (keyw == "ConnectionPointAtoms") {
      this->rotConnInfo.connection_point_atoms.clear();
      for (int i =1; i < fields.size(); ++i) {
	int connectionAtomN = atoi(fields[i].c_str());
	this->rotConnInfo.connection_point_atoms.push_back(connectionAtomN);
      }
    }
  }  

  /* Error conditions. */
  
  if (this->rotConnInfo.anchor_pts.size() == 0) {
    cerr << "No anchor points specified!  Program quitting." << endl;
    exit(8);
  }
  if (this->rotConnInfo.anchor_pts.size() != 3) {
    cerr << "Number of anchor points specified not equal to 3!  Program quitting." << endl;
    exit(8);
  }
  if (this->rotConnInfo.side_chain_atoms.size() == 0) {
    cerr << "No side_chain_atoms, i.e. atoms that correspond to different conformations of the structure, are specified.  Program quitting." << endl;
    exit(8);
  }
  if (this->rotConnInfo.atoms_of_exact_match.size() == 0) {
    cerr << "No atoms of exact match specified!  Program quitting" << endl;
    exit(8);
  }
  /* Insert other erro conditions here */


}

void Rotlib::_populate_atom_mapping_info(stringV atom_mapping_lines) {
  /* populates the following structures: atom_n_map, atom_n_label_map */
  this->rotConnInfo.atom_n_map.clear();
  this->rotConnInfo.atom_n_label_map.clear();

  for (stringVConstItr itr = atom_mapping_lines.begin(); itr != atom_mapping_lines.end(); ++itr) {
    stringV fields;
    split(*itr, fields);
    string atom_label = fields[0];
    int atom_n_in_rotamer = atoi(fields[1].c_str());
    int atom_n_in_protein = atoi(fields[2].c_str());

    this->rotConnInfo.atom_n_map.insert(make_pair(atom_n_in_rotamer, atom_n_in_protein));
    this->rotConnInfo.atom_n_label_map.insert(make_pair(atom_n_in_rotamer, atom_label));
    
  }
  assert(this->rotConnInfo.atom_n_map.size() == this->rotConnInfo.atom_n_label_map.size());

}

void Rotlib::_populate_connectivity_info(stringV connectivity_info_lines) {
  /* populates the map<int, vector<int> > atom_connectivity_info structure. */
  for (stringVConstItr itr = connectivity_info_lines.begin(); itr != connectivity_info_lines.end(); ++itr) {

    stringV fields;
    split(*itr, fields);

    if (fields[0] == "ORDER") { continue; } 
    else if (fields[0] == "CONECT")
      {
	int base_atom = atoi(fields[1].c_str());
	//cout << "base_atom " << base_atom << endl;
	vector<int> connected_atoms;
	connected_atoms.clear();
	if (fields.size() == 2) {
	  // do nothing; don't populate connected_atoms
	} else {
	  for (int i = 2; i < fields.size(); ++i) {
	    int atom_n = atoi(fields[i].c_str());
	    connected_atoms.push_back(atom_n);
	  }
	}
	this->rotConnInfo.atom_connectivity_info.insert(make_pair(base_atom, connected_atoms));
      }
    
  }

}


/* Hierachical rotamer clustering scheme related functions. */



/*
bool Rotlib::cmp_rotlib_E(Rotamer* r1, Rotamer* r2) {

  return ( r1->get_rotlib_E() < r2->get_rotlib_E() );

}
*/

/* *********** AARotlib **********************/

AARotlib::AARotlib() {

  Rotlib();
  

}

AARotlib::~AARotlib() {
  //cout << "~AARotlib successful" << endl;
}

AARotamer* AARotlib::get_next_rot() {

  AARotamer* to_be_returned;
  if (rot_itr == rot_v.begin() ) {
    to_be_returned = ( (AARotamer*)(*rot_itr) );
    ++rot_itr;
    return to_be_returned;
    
  } else {
    if (rot_itr != rot_v.end() ) {
      to_be_returned = ( (AARotamer*)(*rot_itr) );
      ++rot_itr;
      return to_be_returned;
    } else {
      return NULL;
    }
    ++rot_itr;
    
  }
}

AARotamer* AARotlib::get_current_rot() {

  if (rot_itr != this->rot_v.end() ) {
    return ( (AARotamer*)(*rot_itr) );
  } else {
    return NULL;
  }


}

AARotamer* AARotlib::get_rot(int rotamer_n) {

  int current_pstn = (*rot_itr)->get_rotamer_n();
  AARotamer* to_be_returned;
  for (this->rot_itr = this->rot_v.begin(); rot_itr != rot_v.end(); ++rot_itr) {
    if ( (*rot_itr)->get_rotamer_n() == rotamer_n) {
      //AARotamer* to_be_returned = (AARotamer*)(*rot_itr);   // now rot_itr is at right position.
      to_be_returned = (AARotamer*)(*rot_itr);
      break;      
    }
  }
    
  if (rot_itr == rot_v.end()) {
    to_be_returned = NULL;
  }

  this->set_rot_pstn(current_pstn);
  return to_be_returned;

}

AARotamer* AARotlib::reset_rot_pstn() {

  this->rot_itr = rot_v.begin();
  return (AARotamer*)(*rot_itr);

}

AARotamer* AARotlib::set_rot_pstn(int residue_n) {

  for (this->rot_itr = this->rot_v.begin(); rot_itr != rot_v.end(); ++rot_itr) {
    if ((*rot_itr)->get_rotamer_n() == residue_n) {
      break;                        // now at correct position.
    }
  }

  if (rot_itr == rot_v.end()) {
    this->reset_rot_pstn();
    return NULL;
  } else {
    return (AARotamer*)(*rot_itr);
  }

}

AARotamer* AARotlib::get_next_rot_with_empty_lattice_E_below(double CUTOFF_E) {
  
  assert(CUTOFF_E >= 0);
  AARotamer* to_be_returned;

  if (rot_itr == rot_v.begin() ) { // head of vector
    to_be_returned = ( (AARotamer*)(*rot_itr) );
    ++rot_itr;
    return to_be_returned;
    
  } 
  else {			// any other position

    if (rot_itr != rot_v.end() ) {
      to_be_returned = ( (AARotamer*)(*rot_itr) );
      ++rot_itr;
      // if less than CUTOFF_E, return
      if ( to_be_returned->get_empty_lattice_E() <= CUTOFF_E) {
	return to_be_returned;
      } 
      // else, return NULL
      else {
	return NULL;
      }
    } 
    else {			// if at end of vector
      return NULL;
    }
    //    ++rot_itr;			// should never reach here
    
  }

}

AARotamer* AARotlib::get_empty_lattice_E_rot(int EL_rank) {
  /* returns NULL if EL_rank is not found. */
  //  cout << "in get_empty_lattice_E_rot: EL_rank is " << EL_rank << endl;
  //  assert(rot_itr != rot_v.end() );
  //int current_pstn = (*rot_itr)->get_rotamer_n();
  //  cout << "in get_empty_lattice_E_rot: current_pstn is " << current_pstn << endl;

  vector<Rotamer*>::const_iterator temp_RotItr = this->rot_v.begin();

  AARotamer* to_be_returned;
  //  for (this->rot_itr = this->rot_v.begin(); rot_itr != rot_v.end(); ++rot_itr) {
  for (; temp_RotItr != rot_v.end(); ++temp_RotItr) {
    int rank = (*temp_RotItr)->get_empty_lattice_energy_rank();
    //    cout << "Rank for temp_RotItr (when i'm trying to find the right rotamer is: " << rank << " and the one i'm looking for is " << EL_rank << endl;
    //    cout << "  and, it has an energy of: " << (*temp_RotItr)->get_empty_lattice_E() << endl;

    //    if ( (*rot_itr)->get_empty_lattice_energy_rank() == rotamer_n) {
    if ( (*temp_RotItr)->get_empty_lattice_energy_rank() == EL_rank) {
      //AARotamer* to_be_returned = (AARotamer*)(*rot_itr);   // now rot_itr is at right position.
      to_be_returned = (AARotamer*)(*temp_RotItr);
      break;      
    }
  }
    
  if (rot_itr == rot_v.end()) {
    to_be_returned = NULL;
  }

  //  cout << "in get_empty_lattice_E_rot: this rotamer has energy " << to_be_returned->get_empty_lattice_E() << endl;
  //this->set_rot_pstn(current_pstn);
  return to_be_returned;
}

void AARotlib::center_CA() {

  vector<Rotamer*>::iterator itr;
  for (itr = this->rot_v.begin(); itr != rot_v.end(); ++itr) {

    ((AARotamer*)(*itr))->center_CA();

  }

}

void AARotlib::calc_all_PHI() {

  vector<Rotamer*>::iterator itr;
  for (itr = this->rot_v.begin(); itr != rot_v.end(); ++itr) {

    ((AARotamer*)(*itr))->calc_PHI();
    
  }

}

void AARotlib::calc_all_PSI() {

  vector<Rotamer*>::iterator itr;
  for (itr = this->rot_v.begin(); itr != rot_v.end(); ++itr) {

    ((AARotamer*)(*itr))->calc_PSI();
    
  }

}

/* ********** NaturalAARotlib ****************/

NtrlAARotlib::NtrlAARotlib() {



}

/**
 * Constructor with a string passed in as only argument initializes the Rotlib object of the rotamer library with that string name.
 */

NtrlAARotlib::NtrlAARotlib(string in_Rotlib_File) {
  
  AARotlib();
  this->library_name = "";

  /* populates the vector<vector<string> > line_vv member variable. */
  cout << "loading rotamer library " << in_Rotlib_File << "..." << endl;
  populate_line_vv(in_Rotlib_File);
  /* Populates rot_v. */
  vector<vector<string> >::const_iterator itr;
  for (itr = line_vv.begin(); itr != line_vv.end(); ++itr) {
    Rotamer* new_r = new AARotamer(*itr);
    this->rot_v.push_back(new_r);
    /* now set Rotamer declaredInRotlibScope to true */
    new_r->setDeclaredInRotlibScope(true);
  }

  /* Now determine which type of AA this Rotlib contains. */
  this->resName = static_cast<AARotamer*>(rot_v[0])->get_resName();

  /* Initialize rot_itr. */
  rot_itr = rot_v.begin();

  /* Make internal connectivities. */
  string AA = scream_tools::one_letter_AA(this->resName);
  //string default_rotConnInfo_path = "/project/Biogroup/Software/SCREAM/lib/NtrlAARotConn/";
  //vcvicek
  char * SCREAM_NEW_CNN = getenv("SCREAM_NEW_CNN");
  if (SCREAM_NEW_CNN == NULL) {printf("error: enviromental variable SCREAM_NEW_CNN is not set \n"); exit(1);}
  
  string default_rotConnInfo_path = string(SCREAM_NEW_CNN);
  string connectivity_file = default_rotConnInfo_path + AA + ".cnn";
  this->store_connectivity_info(connectivity_file);

  for (vector<Rotamer*>::iterator itr = this->rot_v.begin(); itr != this->rot_v.end(); ++itr) {
    ScreamAtomV sc_atoms = (*itr)->get_sc_atoms();
    scream_tools::make_connectivity(sc_atoms, this->rotConnInfo.atom_connectivity_info);
  }

}


NtrlAARotlib::NtrlAARotlib(string in_Rotlib_File, string cnn_file) {
  Debug debugInfo("NtrlAARotlib::NtrlAARotlib(string in_Rotlib_File, string cnn_file)");


  AARotlib();

  this->library_name = "";

  debugInfo.out("Connectivity file is: " + cnn_file);

  /* populates the vector<vector<string> > line_vv member variable. */
  debugInfo.out("loading rotamer library " + in_Rotlib_File + "...");
  populate_line_vv(in_Rotlib_File);
  /* Populates rot_v. */
  vector<vector<string> >::const_iterator itr;
  for (itr = line_vv.begin(); itr != line_vv.end(); ++itr) {
    debugInfo.out("before new Rotamer");
    Rotamer* new_r = new AARotamer(*itr);
    //new_r->library_name = this->library_name;
    debugInfo.out("after new Rotamer" );
    this->rot_v.push_back(new_r);
    /* now set Rotamer declaredInRotlibScope to true */
    new_r->setDeclaredInRotlibScope(true);
  }

  /* Now determine which type of AA this Rotlib contains. */
  this->resName = static_cast<AARotamer*>(rot_v[0])->get_resName();

  /* Initialize rot_itr. */
  rot_itr = rot_v.begin();

  /* Make internal connectivities. */
  string AA = scream_tools::one_letter_AA(this->resName);
  //string default_rotConnInfo_path = "/project/Biogroup/Software/SCREAM/lib/NtrlAARotConn/";
  //string connectivity_file = default_rotConnInfo_path + AA + ".cnn";
  //this->store_connectivity_info(connectivity_file);
  this->store_connectivity_info(cnn_file);


  for (vector<Rotamer*>::iterator itr = this->rot_v.begin(); itr != this->rot_v.end(); ++itr) {
    ScreamAtomV sc_atoms = (*itr)->get_sc_atoms();
    scream_tools::make_connectivity(sc_atoms, this->rotConnInfo.atom_connectivity_info);
  }

}

void NtrlAARotlib::setup_library() {

  this->center_CA();
  //  this->calc_all_PSI();
  //  this->calc_all_PHI();
  this->assign_atom_fftype();
  this->assign_charges(string("CHARM22"));
  this->assign_lone_pair();

}


NtrlAARotlib::~NtrlAARotlib() {
  // move this to Rotlib destructor
  /*vector<Rotamer*>::const_iterator itr;
  for (itr = rot_v.begin(); itr != rot_v.end(); ++itr) {
    delete *itr;
  }
  */
  //  cout << "~NtrlAARotlib successful" << endl;
}


void NtrlAARotlib::assign_atom_fftype() {

  vector<Rotamer*>::const_iterator itr;
  for (itr = rot_v.begin(); itr != rot_v.end(); ++itr) {
    (*itr)->assign_atom_fftype();
    
  }

}

void NtrlAARotlib::assign_charges(string scheme) {

  vector<Rotamer*>::const_iterator itr;
  for (itr = rot_v.begin(); itr != rot_v.end(); ++itr) {
    ((AARotamer*)(*itr))->assign_charges(scheme);
    
  }

}

void NtrlAARotlib::assign_lone_pair() {

  vector<Rotamer*>::const_iterator itr;
  for (itr = rot_v.begin(); itr != rot_v.end(); ++itr) {
    (*itr)->assign_lone_pair();
    
  }


}

void NtrlAARotlib::append_to_filehandle(ofstream* rotlib_s) {
  
  
  vector<Rotamer*>::const_iterator itr;
  for (itr = rot_v.begin(); itr != rot_v.end(); ++itr) {
    (*itr)->append_to_filehandle(rotlib_s);
    
  }

}


/* Multiple_NtrlAARotlib
 * Designed for doing mutant designs--all rotamers in one single library.
 */

Multiple_NtrlAARotlib::Multiple_NtrlAARotlib() {

}

Multiple_NtrlAARotlib::Multiple_NtrlAARotlib(string lib_path, int resolution) {

  /* Instantiate map<string, Rotlib*> rotlib_m */

  this->rotlib_m.clear();
  string AA_string("ACDERGHIJKLMNPQRSTVWY");

  string resolution_str = itoa(resolution);
  for (string::const_iterator itr = AA_string.begin(); itr != AA_string.end(); ++itr) {
    char AA_char[2];
    AA_char[0] = *itr;
    AA_char[1] = '\0';
    string AA(AA_char);
    string library_filename = lib_path + AA + "/" + AA + "_" + resolution_str + ".lib";
    NtrlAARotlib* rotlib = new NtrlAARotlib(library_filename);
    this->rotlib_m.insert(make_pair(AA, rotlib));
  }

  /* Then populate vector<Rotamer*> */
  this->library_name = "";
  this->rot_v.clear();
  AARotamer* rot = NULL;

  for (map<string, NtrlAARotlib*>::iterator itr = this->rotlib_m.begin(); itr != this->rotlib_m.end(); ++itr) {
    string AA = itr->first;
    NtrlAARotlib* rotlib = itr->second;
    rotlib->reset_pstn();
    rot = rotlib->get_next_rot();
    while (rot != NULL) {
      this->rot_v.push_back(rot);
      rot = rotlib->get_next_rot();
    }

  }
  
  /* Renumber rotamers. */
  int i =1;
  vector<Rotamer*>::iterator itr = this->rot_v.begin();
  for (; itr != this->rot_v.end(); ++i, ++itr) {
    (*itr)->set_rotamer_n(i);
  }

  /* Initialize rot_itr. */
  rot_itr = this->rot_v.begin();

}

Multiple_NtrlAARotlib::Multiple_NtrlAARotlib(string lib_path, int resolution, vector<std::string> amino_acids) {
  /* Instantiate map<string, Rotlib*> rotlib_m */

  this->rotlib_m.clear();
  //  string AA_string("ACDERGHIJKLMNPQRSTVWY");
  string AA_string("");
  for (vector<std::string>::const_iterator itr = amino_acids.begin(); itr != amino_acids.end(); ++itr) {
    AA_string.insert(0, *itr);
  }

  string resolution_str = itoa(resolution);
  for (string::const_iterator itr = AA_string.begin(); itr != AA_string.end(); ++itr) {
    char AA_char[2];
    AA_char[0] = *itr;
    AA_char[1] = '\0';
    string AA(AA_char);
    string library_filename = lib_path + AA + "/" + AA + "_" + resolution_str + ".lib";
    NtrlAARotlib* rotlib = new NtrlAARotlib(library_filename);
    this->rotlib_m.insert(make_pair(AA, rotlib));
  }

  /* Then populate vector<Rotamer*> */
  this->library_name = "";
  this->rot_v.clear();
  AARotamer* rot = NULL;

  for (map<string, NtrlAARotlib*>::iterator itr = this->rotlib_m.begin(); itr != this->rotlib_m.end(); ++itr) {
    string AA = itr->first;
    NtrlAARotlib* rotlib = itr->second;
    rotlib->reset_pstn();
    rot = rotlib->get_next_rot();
    while (rot != NULL) {
      //cout << "rot->library_name: " << rot->library_name << endl;
      this->rot_v.push_back(rot);
      rot = rotlib->get_next_rot();
    }

  }
  
  /* Renumber rotamers. */
  int i =1;
  vector<Rotamer*>::iterator itr = this->rot_v.begin();
  for (; itr != this->rot_v.end(); ++i, ++itr) {
    (*itr)->set_rotamer_n(i);
  }

  /* Initialize rot_itr. */
  rot_itr = this->rot_v.begin();

}


Multiple_NtrlAARotlib::~Multiple_NtrlAARotlib() {
  for (map<string, NtrlAARotlib*>::iterator itr = this->rotlib_m.begin(); itr != this->rotlib_m.end(); ++itr ) {
    itr->second->rot_v.clear(); // vector<Rotamer*> deleted in destructor ~Rotlib().  Don't want to double delete Rotamer*.
    delete itr->second; // Now that all vector<Rotamer*> members in NtrlAARotlib's are cleared, deleting wouldn't be deleting the same Rotamers a second time.
  }
  
}

RotConnInfo* Multiple_NtrlAARotlib::getRotConnInfo() {
  
  /* Returns the RotConnInfo* of the current Rotamer pointed to by rot_itr */

  /* Need to get the rotamer type of current rotamer */

  string AA_3_letter = static_cast<AARotamer*>(*(this->rot_itr))->get_resName();
  string AA_1_letter = scream_tools::one_letter_AA(AA_3_letter);
  
  return this->rotlib_m[AA_1_letter]->getRotConnInfo();

}

void Multiple_NtrlAARotlib::add_library(std::string library) {

}

void Multiple_NtrlAARotlib::add_library(NtrlAARotlib* rotlib) {

}


/* HIS_NtrlAARotlib 
 * A Specialized rotlib for doing HIS rotamers--singly protonated.
 */

HIS_NtrlAARotlib::HIS_NtrlAARotlib() {
  
}

HIS_NtrlAARotlib::HIS_NtrlAARotlib(string H_lib, string J_lib, string H_cnn_file, string J_cnn_file) : H_AARotlib(H_lib, H_cnn_file), J_AARotlib(J_lib, J_cnn_file) {
  /* Now populate vector<Rotamer*> */
  // Memory issue: set declaredInRotlibScope to false for Rotamer* 's in vector<Rotamer*>.  
  // Remark: not necessary; if i clear H_AARotlib and J_AARotlib's rot_v then Rotamer*'s in  vector<Rotamer*> deleted only once.

  this->library_name = "";

  this->rot_v.clear();
  this->H_AARotlib.reset_pstn();
  this->J_AARotlib.reset_pstn();
  
  AARotamer* rot = H_AARotlib.get_next_rot();
  while (rot != NULL) {
    this->rot_v.push_back(rot);
    rot = H_AARotlib.get_next_rot();
  }
  rot = J_AARotlib.get_next_rot();
  while (rot != NULL) {
    this->rot_v.push_back(rot);
    rot = J_AARotlib.get_next_rot();
  }

  /* Renumber rotamers. */
  int i = 1;
  vector<Rotamer*>::iterator itr = this->rot_v.begin();
  for ( ; itr != this->rot_v.end(); ++i, ++itr) {
    (*itr)->set_rotamer_n(i);
  }

  /* Initialize rot_itr. */
  rot_itr = this->rot_v.begin();

}

HIS_NtrlAARotlib::~HIS_NtrlAARotlib() {
  // need to clear these rot_v's, else double deletes Rotamer*.  Once in this destructor, once in H_AARotlib and J_AARotlib's destructors.
  // this creates a ton of orphan pointers if not careful... but these rot_v Rotamer*'s should have been copied into this->rot_v in the constructor, so should be okay.
  this->H_AARotlib.rot_v.clear();  
  this->J_AARotlib.rot_v.clear();
}

void HIS_NtrlAARotlib::add_rotamer(Rotamer* his_rot) {

  AARotamer* other_his_rot = new AARotamer();
  AARotamer* his_AArot = static_cast<AARotamer*>(his_rot);
  this->_make_other_singly_protonated_HIS(his_AArot, other_his_rot);
  
  if (other_his_rot == NULL) {
    cout << "Adding Histidine rotamer failure! " << endl;
  }

  other_his_rot->setDeclaredInRotlibScope(true);
  his_rot->setDeclaredInRotlibScope(false);

  if (his_AArot->get_resName() == "HIS") other_his_rot->set_resName("HSE");
  if (his_AArot->get_resName() == "HSE") other_his_rot->set_resName("HIS");
  
  this->rot_v.push_back(his_rot);
  this->rot_v.push_back(other_his_rot);

  if (his_rot->library_name == "") his_rot->library_name = this->library_name;
  if (other_his_rot->library_name == "") other_his_rot->library_name = this->library_name;

  //  his_rot->set_rotamer_n(0);
  other_his_rot->set_rotamer_n(this->rot_v.size() -1 );

  // also needs to set "is original rotamer flag to other_his_rot.

}

void HIS_NtrlAARotlib::_make_other_singly_protonated_HIS(AARotamer* his_rot, AARotamer* other_his_rot) {
  /* 1. deepcopy his_rot to other_his_rot (so don't have to worry about memory issues) */
  other_his_rot->deepcopy(*his_rot);
  
  /* 2. check HSE and HSD state, and get relevant atoms. */
  string old_HIS_state; // HSD and HSE.
  SCREAM_ATOM* old_H;
  ScreamAtomV sc_atoms = other_his_rot->get_sc_atoms();
  for (ScreamAtomVConstItr itr = sc_atoms.begin(); itr != sc_atoms.end(); ++itr) {
    if (scream_tools::strip_whitespace( (*itr)->atomLabel ) == "HND1") {
      old_HIS_state = "HSD";
      old_H = *itr;
      break;
    }
    // since these are singly protonated HIS, only one or other can exist in rotamer.
    if (scream_tools::strip_whitespace( (*itr)->atomLabel ) == "HNE2") {
      old_HIS_state = "HSE";
      old_H = *itr;
      break;
    }
  }

  /* 3. add other hydrogen position. 
   *   new_H
   *     |
   *     A
   *    / \ 
   *   B1  B2
   *
   */
   SCREAM_ATOM *A, *B1, *B2; // as illustrated above.
   SCREAM_ATOM* new_H = new SCREAM_ATOM();
   new_H->copy(*old_H);
   new_H->connectivity_m.clear();

   for (ScreamAtomVConstItr itr = sc_atoms.begin(); itr != sc_atoms.end(); ++itr) {
     string stripped_atomLabel = scream_tools::strip_whitespace( (*itr)->atomLabel );
     if (old_HIS_state == "HSD") { // then new atoms are: around HSE: A : NE2, B1: CE1, B2: CD2
       if (stripped_atomLabel == "NE2") A = *itr;
       if (stripped_atomLabel == "CE1") B1 = *itr;
       if (stripped_atomLabel == "CD2") B2 = *itr;
     } else { // old_HIS_state = "HSE".  the new atoms are: around HSD: A: ND1, B1: CB, B2: CE1
       if (stripped_atomLabel == "ND1") A = *itr;
       if (stripped_atomLabel == "CB") B1 = *itr;
       if (stripped_atomLabel == "CE1") B2 = *itr;
     }
   }

   ScreamVector new_H_position = ScreamVector(A) 
     - ( (ScreamVector(B1) - ScreamVector(A)).normalizedVector() + 
	 (ScreamVector(B2) - ScreamVector(A)).normalizedVector() ) * 0.97; // 0.97: distance from ND1 to new H.

   new_H->x[0] = new_H_position[0];
   new_H->x[1] = new_H_position[1];
   new_H->x[2] = new_H_position[2];
   if (old_HIS_state == "HSD") {
     new_H->resName = "HSE";
     new_H->atomLabel = "HNE2";
   }
   if (old_HIS_state == "HSE") {
     new_H->resName = "HSD";
     new_H->atomLabel = "HND1";
   }


   /* 4. Taking care of connectivies and sc_atom stuff. */
   SCREAM_ATOM* base_atom = old_H->connectivity_m.begin()->first;
   base_atom->connectivity_m.erase(old_H);

   new_H->make_bond(A);
   other_his_rot->get_sc()->remove_atom(old_H->atomLabel);
   other_his_rot->get_sc()->add_atom(new_H->atomLabel, new_H);
   
   /* 5. Renaming res names.  Also Renumber atom numbers. */

   //   other_his_rot->print_Me();
   ScreamAtomV post_bb_atoms = other_his_rot->get_bb_atoms();
   ScreamAtomV post_sc_atoms = other_his_rot->get_sc_atoms();
   

   for (ScreamAtomVItr itr = post_bb_atoms.begin(); itr != post_bb_atoms.end(); ++itr) {
     (*itr)->dump();
     if (old_HIS_state == "HSD") 
       (*itr)->resName = "HSE";
     else 
       (*itr)->resName = "HSD";
   }

   for (ScreamAtomVItr itr = post_sc_atoms.begin(); itr != post_sc_atoms.end(); ++itr) {
     if (old_HIS_state == "HSD") 
       (*itr)->resName = "HSE";
     else 
       (*itr)->resName = "HSD";
   }

   // then renumber atom numbers.

   /* 6.. Taking care of charges.  Working on this! */
   
   
}
