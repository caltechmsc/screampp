#include "scream_atom.hpp"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
//vcvicek
#include <stdlib.h>


std::map<std::string, std::string> SCREAM_ATOM::AA321_map;


SCREAM_ATOM::SCREAM_ATOM() {
  /* Default Constructor; should initialize everything to zero.  The following two items are set to zero for now. */
  this->lone_pair = 0;
  this->flags = 0;
  this->atoms_connected = 0;

  this->_init_AA321_map();

  this->vdw_r = -600;
  this->vdw_d = -600;
  this->vdw_s = -600;

  this->hb_da = -600;
  this->delta = -600;

}

SCREAM_ATOM::SCREAM_ATOM(const string bgf_line) {
  /*
           1        2         3         4         5         6         7         8         9  
0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM       1  N    VAL A   86    2.68796  16.15051 -19.20364 N_R    3 0 -0.47000 0  -1  11.99086
  */

  this->_init_AA321_map();

  this->feed_me(bgf_line);

  this->library_name = "";

  this->flags = 0;

  this->vdw_r = -600;
  this->vdw_d = -600;
  this->vdw_s = -600;

  this->hb_da = -600;
  this->delta = -600;
  // SCREAM delta values initialized elsewhere: in scream_EE object.
}

SCREAM_ATOM::SCREAM_ATOM(SCREAM_ATOM* atom) {

  

}

SCREAM_ATOM::~SCREAM_ATOM() {
  //  cout << "~SCREAM_ATOM()" << endl;
  // do nothing for now
}

void SCREAM_ATOM::pdb_init(string pdb_line) {

  this->feed_me_pdb(pdb_line);


  this->library_name = "";

  this->flags = 0;


}

SCREAM_ATOM& SCREAM_ATOM::copy(const SCREAM_ATOM& atom) {

  this->keyw = atom.keyw;
  this->atomLabel = atom.atomLabel;
  this->stripped_atomLabel = atom.stripped_atomLabel;
  this->isSC_Flag = atom.isSC_Flag;
  this->isAAResAtom = atom.isAAResAtom; 
  this->atomType = atom.atomType;
  this->stripped_atomType = atom.stripped_atomType;
  this->resName = atom.resName;
  this->oneLetterResName = atom.oneLetterResName;
  this->chain = atom.chain;
  this->resNum = atom.resNum;
  this->atoms_connected = atom.atoms_connected;
  this->lone_pair = atom.lone_pair;
  
  this->x[0] = atom.x[0];
  this->x[1] = atom.x[1];
  this->x[2] = atom.x[2];
  this->q[0] = atom.q[0];
  this->q[1] = atom.q[1];
  this->n = atom.n;
  this->type = atom.type;
  this->flags = atom.flags;
  this->m = atom.m;
  this->vchg2 = atom.vchg2;
  
  this->library_name = atom.library_name;

  this->vdw_r = atom.vdw_r;
  this->vdw_d = atom.vdw_d;
  this->vdw_s = atom.vdw_s;
  this->vachg = atom.vachg;
  this->vrchg = atom.vrchg;
  this->hb_da = atom.hb_da;

  this->delta = atom.delta;

  /* Careful with the following */
  this->connectivity_m.clear();
  for (map<SCREAM_ATOM*, int>::const_iterator itr = atom.connectivity_m.begin(); itr != atom.connectivity_m.end(); ++itr) {
    this->connectivity_m.insert(*itr);
  }

  return *this;

}

void SCREAM_ATOM::setAtomLabel(string atomLabel) {
  this->atomLabel = atomLabel;

  char al[256];
  strcpy(al, this->atomLabel.c_str());

  this->stripped_atomLabel = string(strtok(al, " "));
}

void SCREAM_ATOM::setAtomType(string atomType) {
  this->atomType = atomType;

  char al[256];
  strcpy(al, this->atomType.c_str());

  this->stripped_atomType = string(strtok(al, " "));
}


SCREAM_ATOM& SCREAM_ATOM::copyJustCoords(const SCREAM_ATOM& atom) {
  //cout << atom.x[0] << " " << atom.x[1] << " " << atom.x[2] << endl;

  this->x[0] = atom.x[0];
  this->x[1] = atom.x[1];
  this->x[2] = atom.x[2];

}

void SCREAM_ATOM::feed_me(const string bgf_line) {

  this->keyw = bgf_line.substr(0, 6);
  this->n = atoi(bgf_line.substr(7, 5).c_str());
  this->atomLabel = bgf_line.substr(13, 4);
  char al[256];
  strcpy(al, this->atomLabel.c_str());
  this->stripped_atomLabel = string(strtok(al, " "));
  
  this->resName = bgf_line.substr(19, 3);
  this->oneLetterResName = this->AA321_map[resName];
  this->chain = bgf_line.substr(23, 1);  
  this->resNum = atoi(bgf_line.substr(25, 4).c_str());

  std::sscanf(bgf_line.substr(30, 10).c_str(), "%lf", &this->x[0]);
  std::sscanf(bgf_line.substr(40, 10).c_str(), "%lf", &this->x[1]);
  std::sscanf(bgf_line.substr(50, 10).c_str(), "%lf", &this->x[2]);

  this->atomType = bgf_line.substr(61, 5);
  char s[256];
  strcpy(s, this->atomType.c_str());
  this->stripped_atomType = string(strtok(s, " "));
  this->atoms_connected = atoi(bgf_line.substr(68, 1).c_str());
  this->lone_pair = atoi(bgf_line.substr(70, 1).c_str());

  std::sscanf(bgf_line.substr(72, 8).c_str(), "%lf", &this->q[0]);

  

}

void SCREAM_ATOM::feed_me_pdb(const string pdb_line) {
  //ATOM      1  N   ALA     1      14.419  46.428  83.486  1.00 47.38
  this->keyw = pdb_line.substr(0, 6);
  this->n = atoi(pdb_line.substr(6, 5).c_str());
  this->atomLabel = pdb_line.substr(12, 4);
  char al[256];
  strcpy(al, this->atomLabel.c_str());
  this->stripped_atomLabel = string(strtok(al, " "));
  
  this->resName = pdb_line.substr(17, 3);
  this->oneLetterResName = this->AA321_map[resName];
  this->chain = pdb_line.substr(21, 1);  
  this->resNum = atoi(pdb_line.substr(22, 4).c_str());

  std::sscanf(pdb_line.substr(30, 8).c_str(), "%lf", &this->x[0]);
  std::sscanf(pdb_line.substr(38, 10).c_str(), "%lf", &this->x[1]);
  std::sscanf(pdb_line.substr(46, 10).c_str(), "%lf", &this->x[2]);

  if (pdb_line.length() >= 60) {
    std::sscanf(pdb_line.substr(55, 5).c_str(), "%lf", &this->occupancy);
  } 
  if (pdb_line.length() >= 66) {
    std::sscanf(pdb_line.substr(60, 6).c_str(), "%lf", &this->BFactor);
  }

}

/*
void feed_bond_info(const string bgf_connect_line) {
  // can't use scream_tools::split because don't want to get into dependency problems
  stringstream line_ss(bgf_connect_line);
  string entry;
  vector<string> after_split;
  while (!line_ss.eof()) {
    line_ss >> entry;
    after_split.push_back(entry);
  }  
  
}
*/

bool SCREAM_ATOM::make_bond(SCREAM_ATOM* in_atom, int type) {
  //    cout << "in SCREAM_ATOM::make_bond" << endl;
  //    cout << n << endl;
  //    cout << in_atom->n << endl;

  if (this == in_atom) {                                             // don't make self-bond
    cout << "making bond with self not allowed.  Printing passed in atom and self atom (in that order):" << endl;

    this->dump();
    in_atom->dump();
    return false;
  }

  if (this->connectivity_m.find(in_atom) == connectivity_m.end() or this->connectivity_m.empty() ) {  // check if connection is not already there
    this->connectivity_m.insert(make_pair(in_atom, type));
    in_atom->connectivity_m.insert(make_pair(this, type));
    //    cout << "here" << endl;
    return true;
  } else {
    return false;                                                    // returns false if connectivity already established.
  }
}

bool SCREAM_ATOM::delete_bond(SCREAM_ATOM* a) {
  if (a == this) {
    cout << "Deleting bond with self not allowed." << endl;
    return false;
  }
  else {
    map<SCREAM_ATOM*, int>::iterator itr = this->connectivity_m.find(a);
    if (itr == this->connectivity_m.end())
      return false;
    else {
      this->connectivity_m.erase(a);
      a->connectivity_m.erase(this);
      return true;
    }

  }
}


double SCREAM_ATOM::distance(SCREAM_ATOM* atom) {

  return sqrt(this->distance_squared(atom));

}

double SCREAM_ATOM::distance_squared(SCREAM_ATOM* atom) {

  //  cout << " in distance_squared: " << endl;
  //  atom->dump();
  //  cout << "Atom x[0] : " << atom->x[0];
  //  cout << "Atom x[1] : " << atom->x[1];
  //  cout << "Atom x[2] : " << atom->x[2];
  //    
  //  cout << "this x[0]: " << this->x[0];
  //  cout << "this x[1] : " << this->x[1];
  //  cout << "this x[2] : " << this->x[2];

  //return ( pow((this->x[0] - atom->x[0]),2) + pow((this->x[1] - atom->x[1]),2) + pow((this->x[2] - atom->x[2]),2) );
  double x_ = this->x[0] - atom->x[0];
  double y_ = this->x[1] - atom->x[1];
  double z_ = this->x[2] - atom->x[2];
  return x_ * x_ + y_ * y_ + z_ * z_;


}

double SCREAM_ATOM::worst_clash_dist(vector<SCREAM_ATOM*>& atom_list, SCREAM_ATOM** worst_clash_atom_p) {  
  //#double SCREAM_ATOM::worst_clash_dist(vector<SCREAM_ATOM*>& atom_list) {
  vector<SCREAM_ATOM*>::iterator itr = atom_list.begin();
  double worst_dist = 999999;
  SCREAM_ATOM* atom;
  *worst_clash_atom_p = NULL;  // this is the SCREAM_ATOM that gets passed in--keeps track of the atom that has the worst clash distance.  Potential memory leak?  This is a pointer to a pointer.
  

  for (; itr != atom_list.end(); ++itr) {
    atom = *itr;
    double tmp_dist = this->distance_squared(atom);
    if (tmp_dist < worst_dist) {
      worst_dist = tmp_dist;
      *worst_clash_atom_p = *itr;
    }
  }

  //  (*worst_clash_atom_p)->dump();
  return sqrt(worst_dist);
}

void SCREAM_ATOM::dump() const {
    
  //    printf("%6s %5d %4s  %3s %1s %4d %10.5lf%10.5lf%10.5lf %-5s  %1d %1d %8.5lf\n", keyw.c_str(), n, atomLabel.c_str(), resName.c_str(), chain.c_str(), resNum, x[0], x[1], x[2], atomType.c_str(), atoms_connected, lone_pair, q[0]); 

  setiosflags(std::ios::fixed);
  setprecision(5);

  cout << setiosflags(std::ios::fixed) << setw(6) << keyw << ' ' << setw(5) << n << ' ' << setw(4) << atomLabel << ' ' << ' ' << setw(3) << resName << ' ' << setw(1) << chain << ' ' << setw(4) << resNum << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10)  << x[0] << setw(10) << x[1] << setw(10) << setprecision(5) << x[2] << ' ' << setw(5) << atomType << ' ' << ' ' << setw(1) << atoms_connected << ' ' << setw(1) << lone_pair << setprecision(5) << setw(9) << q[0] << endl;

  flush(cout);
    //std::cout << setw(6) << keyw << ' ' << setw(5) << n << ' ' << setw(5) << atomLabel << ' ' << setw(3) << resName << ' ' << setw(10) << setprecision(5) << x[0] << ' ' << setw(10) << setprecision(5) << x[1] << ' ' << setw(10) << setprecision(5) << x[2] << ' ' << setw(5) << atomType << endl;

}

void SCREAM_ATOM::pdb_dump() const {

  //  cout << setw(6) << string("HETATM") << setw(5) << n << ' ' << setw(4) << atomLabel << ' ' << setw(3) << resName << ' ' << setw(1) << chain << setw(4) << resNum << setw(4) << string("    ") << setprecision(3) << setiosflags(std::ios::fixed) << setw(8) << x[0] << setw(8) << x[1] << setw(8) << x[2] << endl;

  cout << setw(6) << keyw << setw(5) << n << ' ' << setw(4) << atomLabel << ' ' << setw(3) << resName << ' ' << setw(1) << chain << setw(4) << resNum << setw(4) << string("    ") << setprecision(3) << setiosflags(std::ios::fixed) << setw(8) << x[0] << setw(8) << x[1] << setw(8) << x[2] << endl;

}

string SCREAM_ATOM::return_bgf_line() const{
  char buf[500];  // 500 > 80, otherwise arbitrary.
  sprintf(buf, "%-6s %5d %-4s  %-3s %1s %4d %10.5lf%10.5lf%10.5lf %-5s  %1d %1d %8.5lf\n", keyw.c_str(), n, atomLabel.c_str(), resName.c_str(), chain.c_str(), resNum, x[0], x[1], x[2], atomType.c_str(), atoms_connected, lone_pair, q[0]); 
  
  string buf_cpp(buf);
  return buf_cpp;
	  
}

string SCREAM_ATOM::return_pdb_line() const {
  //TATM12345  N   VAL A  86       2.688  16.151 -19.204
  char buf[500];  // 500 > 80, otherwise arbitrary.
  sprintf(buf, "%-6s%5d %-4s %-3s %1s%4d    %8.3lf%8.3lf%8.3lf\n", keyw.c_str(), n, atomLabel.c_str(), resName.c_str(), chain.c_str(), resNum, x[0], x[1], x[2]);
  string buf_cpp(buf);
  return buf_cpp;

}


void SCREAM_ATOM::append_to_filehandle(ostream* ostream_p) const {
  char buf[500];  // 500 > 80, otherwise arbitrary.
  sprintf(buf, "%-6s %5d %-4s  %-3s %1s %4d %10.5lf%10.5lf%10.5lf %-5s  %1d %1d %8.5lf\n", keyw.c_str(), n, atomLabel.c_str(), resName.c_str(), chain.c_str(), resNum, x[0], x[1], x[2], atomType.c_str(), atoms_connected, lone_pair, q[0]); 
  
  string buf_cpp(buf);
  *ostream_p << buf_cpp;

  // Wouldn't work on borg machines.

  // how set precision so that the trailing zero's are printed out as well?

  //setiosflags(std::ios::fixed);
  //setprecision(5);
  
  //if (keyw == "HETATM") {

  //  *ostream_p << setw(6) << ("HETATM") << ' ' << setw(5) << n << ' ' << setw(4) << atomLabel << ' ' << ' ' << setw(3) << resName << ' ' << setw(1) << chain << ' ' << setw(4) << resNum << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10)  << x[0] << setw(10) << x[1] << setw(10) << setprecision(5) << x[2] << ' ' << setw(5) << atomType << ' ' << ' ' << setw(1) << atoms_connected << ' ' << setw(1) << lone_pair << setprecision(5) << setw(9) << q[0] << endl;

  //  } else {
  //    *ostream_p << setw(6) << ("ATOM  ") << ' ' << setw(5) << n << ' ' << setw(4) << atomLabel << ' ' << ' ' << setw(3) << resName << ' ' << setw(1) << chain << ' ' << setw(4) << resNum << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10)  << x[0] << setw(10) << x[1] << setw(10) << setprecision(5) << x[2] << ' ' << setw(5) << atomType << ' ' << ' ' << setw(1) << atoms_connected << ' ' << setw(1) << lone_pair << setprecision(5) << setw(9) << q[0] << endl;
  //  }

  /* ATOM/HETATM keyw left aligned */

  //  ostream_p->setf(std::ios_base::left,
  //	  std::ios_base::adjustfield);
  //  *ostream_p << setiosflags(std::ios::fixed) << setw(6) << keyw;  

  /* atom number right aligned */
  //  ostream_p->setf(std::ios_base::right,
  //		  std::ios_base::adjustfield);
//  *ostream_p << ' ' << setw(5) << n << ' ' ;
  
  /* atomLabel, tricky.  right aligned for now. 2 spaces. */
  //ostream_p->setf(std::ios_base::left,
//	  std::ios_base::adjustfield);
//if (atomLabel.substr(0,1) == string("H")) {
//  *ostream_p << setw(4) << atomLabel << ' ' << ' ';
//} else {
//  *ostream_p << ' ' << setw(3) << atomLabel << ' ' << ' ';
//}

  /* everything else, right aligned as default */
  //ostream_p->setf(std::ios_base::right,
//	  std::ios_base::adjustfield);
//*ostream_p << setw(3) << resName << ' ' << setw(1) << chain << ' ' << setw(4) << resNum << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10)  << x[0] << setw(10) << x[1] << setw(10) << setprecision(5) << x[2] << ' ';

  /* atom type left aligned */
  //ostream_p->setf(std::ios_base::left,
//	  std::ios_base::adjustfield);
//*ostream_p << setw(5) << atomType << ' ' << ' ' ;

  /* back to right aligned */
  //ostream_p->setf(std::ios_base::right,
//	  std::ios_base::adjustfield);
//*ostream_p << setw(1) << atoms_connected << ' ' << setw(1) << lone_pair << setprecision(5) << setw(9) << q[0] << endl;

}

void SCREAM_ATOM::pdb_append_to_filehandle(ostream* ostream_p) const {
  /*                      12345
012345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM      1  N   MET A   1      62.090 101.800  60.911  1.00 30.47           N  
  */
  char buf[500];  // 500 > 80, otherwise arbitrary.
  sprintf(buf, "%-6s%5d %-4s %-3s %1s%4d    %8.3lf%8.3lf%8.3lf\n", keyw.c_str(), n, atomLabel.c_str(), resName.c_str(), chain.c_str(), resNum, x[0], x[1], x[2]);
  string buf_cpp(buf);
  *ostream_p << buf_cpp;

  // Reason for not using streams: choose lowest common denominator because borg doesn't have the right ++C STL.

  //  *ostream_p << setw(6) << string("HETATM") << setw(5) << n << ' ' << setw(4) << atomLabel << ' ' << setw(3) << resName << ' ' << setw(1) << chain << setw(4) << resNum << setw(4) << string("    ") << setprecision(3) << setiosflags(std::ios::fixed) << setw(8) << x[0] << setw(8) << x[1] << setw(8) << x[2] << endl;

  //*ostream_p << setiosflags(std::ios::fixed) << setw(6) << keyw << setw(5) << n << ' ' << setw(4) << atomLabel << ' ' << setw(3) << resName << ' ' << setw(1) << chain << setw(4) << resNum << setw(4) << string("    ") << setprecision(3) << setiosflags(std::ios::fixed) << setw(8) << x[0] << setw(8) << x[1] << setw(8) << x[2] << endl;

}

void SCREAM_ATOM::append_to_ostream_connect_info(ostream* ostream_p) const {

  *ostream_p << string("CONECT") << setw(6) << n;
  for (map<SCREAM_ATOM*, int>::const_iterator itr = connectivity_m.begin(); itr != connectivity_m.end(); ++itr) {
    *ostream_p << setw(6) << itr->first->n ;
  }
  *ostream_p << endl;

}

void SCREAM_ATOM::pdb_append_to_ostream_connect_info(ostream* ostream_p) const {
  char CONECT[500];
  strcpy(CONECT, "CONECT");
  *ostream_p << CONECT;
  
  char this_n[500];
  sprintf(this_n, "%5d", this->n);
  *ostream_p << string(this_n);

  for (map<SCREAM_ATOM*, int>::const_iterator itr = connectivity_m.begin(); itr != connectivity_m.end(); ++itr) {
    char buf[500];
    sprintf(buf, "%5d", itr->first->n);
    *ostream_p << string(buf);
  }
  *ostream_p << endl;


  /*
  *ostream_p << string("CONECT") << setw(6) << n;
  for (map<SCREAM_ATOM*, int>::const_iterator itr = connectivity_m.begin(); itr != connectivity_m.end(); ++itr) {
    *ostream_p << setw(5) << itr->first->n;
  }
  *ostream_p << endl;
  */
}

void SCREAM_ATOM::_init_AA321_map() {
  if (AA321_map.empty()) {
    AA321_map["ALA"] = "A";
    //AA321_map["HSP"] = "B"; // doesn't have B library anyway. // Adam: now it does
    AA321_map["HSP"] = "B";
    AA321_map["CYS"] = "C";
    AA321_map["CYX"] = "C";
    AA321_map["ASP"] = "D";
    AA321_map["GLU"] = "E";
    AA321_map["PHE"] = "F";
    AA321_map["GLY"] = "G";
    AA321_map["HIS"] = "H";
    AA321_map["ILE"] = "I";
    AA321_map["HSE"] = "H";
    AA321_map["HSD"] = "J"; 
    AA321_map["LYS"] = "K";
    AA321_map["LEU"] = "L";
    AA321_map["MET"] = "M";
    AA321_map["ASN"] = "N";
    AA321_map["PRO"] = "P";
    AA321_map["GLN"] = "Q";
    AA321_map["ARG"] = "R";
    AA321_map["SER"] = "S";
    AA321_map["THR"] = "T";
    AA321_map["VAL"] = "V";
    AA321_map["TRP"] = "W";
    AA321_map["TYR"] = "Y";
  }
}
