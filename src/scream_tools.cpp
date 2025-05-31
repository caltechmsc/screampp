#include <cstdlib>
#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_tools.hpp"
#include <math.h>
#include "scream_vector.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <string.h>

using namespace std;

string scream_tools::atom_label(const string bgf_line) {

  stringstream l_ss(bgf_line.substr(13, 5));
  string atom_label;
  //atom_label << l_ss;
  l_ss >> atom_label;
  return atom_label; 
  
}

int scream_tools::res_number(const string bgf_line)  {

  return atoi(bgf_line.substr(25, 4).c_str());

}

string scream_tools::res_type(const string bgf_line) {

  return bgf_line.substr(19,3);

}

string scream_tools::chain_desig(const string bgf_line) {
  string chain_desig = bgf_line.substr(23,1);
  if ( chain_desig != string(" ")) {
    return chain_desig;
  } else {
    return chain_desig;
  }

}

string scream_tools::strip_whitespace(string s) {

  size_t start = s.find_first_not_of(" ");  
  size_t end = s.find_last_not_of(" ");

  return s.substr(start, end-start+1);

}

vector<string> scream_tools::split(string line) {

  vector<string> after_split;
  const char* c_line = line.c_str();
  char cc_line[40000];
  strcpy(cc_line, c_line);
  
  //  cout << "here" << endl;

  char* w;
  w = strtok(cc_line, " \t\n");
  
  if (w != NULL)  {
    after_split.push_back(string(w));
  } else {
    return after_split;
  }

  //  printf("%s  after w!=NULL test\n", w);

  while ( (w = strtok(NULL, " \t\n")) != NULL) {
    after_split.push_back(string(w));
  }

  return after_split;

  // below: stringstream just won't work in mpsim-scream integration.  reason: unknown.  possibly compiler option related.
  // reason: caused by -malign-double option.  possibly conflicts with STL.
  /*

  cout << ":::" << line << ":::" << endl;
  //  line = line + "\n";
  stringstream line_ss(line);
  if (line_ss.good()) 
    cout << "ss is good " << endl;
  if (line_ss.eof()) 
    cout << "ss is eof  " << endl;
  if (line_ss.fail()) 
    cout << "ss is fail " << endl;
  if (line_ss.bad()) 
    cout << "ss. if bad " << endl;
  //  cout << " cout << line_ss << endl is ::: " << line_ss << endl;

  cout << "line_ss is :!!!: " << line_ss << endl;

  char test[256];
  while (line_ss >> test) {
    cout << test << endl;
  }

  string entry;
  //  char entry[50];
  vector<string> after_split;
  while (!line_ss.eof()) {
    cout << "does it even get here? " << endl;
    line_ss >> entry;
    cout << entry << endl;
    cout << "entry is " << entry << endl;
    after_split.push_back(entry);
    }
  cout << ":::" << line << ":::" << endl;
  */
  /*  while (line_ss >> entry) {
    cout << "do i get here? " << endl;
    cout << " entry is :::" << entry;
    after_split.push_back(string(entry));
    }*/
  /*
  for (vector<string>::const_iterator itr = after_split.begin(); itr != after_split.end(); ++itr) {
    cout << "printing fields in split:" << *itr;
  }
  */

}

vector<string> scream_tools::char_double_star_to_str_of_vector(char* array_of_ptrs_to_char_array[]) {

  vector<string> vector_of_strings;

  char* ptr_to_char_array = array_of_ptrs_to_char_array[0];
  int i = 0;
  //while (ptr_to_char_array != NULL and i < 5) {
  while (ptr_to_char_array != NULL) {

    vector_of_strings.push_back(string(ptr_to_char_array));
    ptr_to_char_array = array_of_ptrs_to_char_array[++i];
  }

  return vector_of_strings;

}

vector<string> scream_tools::mut_list_string_to_str_of_vector(char* mut_c_string) {

  vector<string> vector_of_strings;
  string mut_cpp_string(mut_c_string);

  cout << "here" << mut_cpp_string << endl; flush(cout);

  unsigned int i = 0;
  while (mut_cpp_string.substr(i,1) == " ") {
    ++i;
  }
  unsigned int j = i;

  if (i == 0) ++i;

  for (; i < mut_cpp_string.size(); ++i) {
    //    cout << "value of i is " << i << endl;
    //    cout << "value of j is " << j << endl;
    if (mut_cpp_string.substr(i,1) == " " and mut_cpp_string.substr(i-1, 1) != " ") {
      vector_of_strings.push_back(mut_cpp_string.substr(j, i-j));
      j = i+1;
      //      cout << "1st" << endl;
    }
    if (mut_cpp_string.substr(i,1) != " " and mut_cpp_string.substr(i-1, 1) == " ") {
      j = i;
      //      cout << "2nd" << endl;
    }
    if (i == mut_cpp_string.size() -1) {
      vector_of_strings.push_back(mut_cpp_string.substr(j, i-j+1));
      //      cout << "3rd" << endl;
    }
  }
  cout << "testing vector_of_strings" << endl;
  for (vector<string>::const_iterator itr = vector_of_strings.begin(); itr != vector_of_strings.end(); ++itr) {
    cout << *itr << endl;
  }

  return vector_of_strings;

}

bool scream_tools::is_bgf_header_line(string line) {

  if (line.substr(0,6) == "BIOGRF" ||
      line.substr(0,6) == "DESCRP" ||
      line.substr(0,6) == "REMARK" ||
      line.substr(0,10) == "FORCEFIELD" ||
      line.substr(0,11) == "FORMAT ATOM")  {
    return true;
  } else {
    return false;
  }

}

bool scream_tools::is_bgf_atom_line(string line) {

  if (line.substr(0,6) == "ATOM  " ||
      line.substr(0,6) == "HETATM" ) {
    return true;
  } else {
    return false;
  }

}

bool scream_tools::is_format_connect_line(string line) {

  if (line.substr(0, 13) == "FORMAT CONECT" || 
      line.substr(0, 12) == "FORMAT ORDER" ) {
    return true;
  } else {
    return false;
  }

}

bool scream_tools::is_connectivity_info(string line) {

  if (line.substr(0, 6) == "CONECT" ||
      line.substr(0,6)  == "ORDER " ) {
    return true;
  } else {
    return false;
  }

}

bool scream_tools::is_AA(const string s) {

  //string res_type_ = scream_tools::res_type(s);
  string res_type_= s;
  if ((res_type_ == "ALA") ||
      (res_type_ == "ARG") ||
      (res_type_ == "ASN") ||
      (res_type_ == "ASP") ||
      (res_type_ == "CYS") ||
      (res_type_ == "GLN") ||
      (res_type_ == "GLU") ||
      (res_type_ == "GLY") ||
      (res_type_ == "HIS") || (res_type_ == "HSE") ||
      (res_type_ == "HSP") ||
      (res_type_ == "ILE") ||
      (res_type_ == "LEU") ||
      (res_type_ == "LYS") ||
      (res_type_ == "MET") ||
      (res_type_ == "PHE") ||
      (res_type_ == "PRO") ||
      (res_type_ == "SER") ||
      (res_type_ == "THR") ||
      (res_type_ == "TRY") ||
      (res_type_ == "TYR") ||
      (res_type_ == "VAL")) {

    return true;

  } else {
    return false;
  }
}

bool scream_tools::is_metal_atom(const string s) {
  
  string atom_label = s;
  atom_label = scream_tools::strip_whitespace(atom_label);

  if (atom_label == "Mg" or atom_label == "MG" or
      atom_label == "Zn" or atom_label == "ZN" or
      atom_label == "Na" or atom_label == "NA"
      ) {
    return true;
  }
  else return false;

}

bool scream_tools::is_SC_atom(const string s) {
  
  //  string atom_label = scream_tools::atom_label(s);
  string atom_label = s;
  atom_label = scream_tools::strip_whitespace(atom_label);

  if (scream_tools::is_BB_atom(atom_label)) {
    return false;
  }
  
  if (scream_tools::is_metal_atom(atom_label)) {
    return false;
  }
  // what's below?  forgot.   whatever it is, it works.
  if (atom_label.substr(0,1) == string("H")) {
    if ( (atom_label.substr(2,1) != string("A")) and (atom_label.substr(2,1) != string(" ")) ) {
      return true;
    }
  } else if (atom_label.substr(0,1) != string("H")) {
    if ( (atom_label.substr(1,1) != string("A")) and (atom_label.substr(1,1) != string(" ")) ) {
      return true;
    }
  } else {
    return false;	
  }

}

bool scream_tools::is_BB_atom(const string s) {

  //string atom_label = scream_tools::atom_label(s);
  static map<std::string, bool> bbAtoms;
  if (bbAtoms.empty()) {
    bbAtoms["HN"] = true;
    bbAtoms["1H"] = true;
    bbAtoms["2H"] = true;
    bbAtoms["3H"] = true;
    bbAtoms["H1"] = true;
    bbAtoms["H2"] = true;
    bbAtoms["H3"] = true;
    bbAtoms["HT1"] = true;
    bbAtoms["HT2"] = true;
    bbAtoms["HT3"] = true;
    bbAtoms["H"] = true;
    bbAtoms["H NT"] = true;
    bbAtoms["N"] = true;
    bbAtoms["NT"] = true;
    bbAtoms["CA"] = true;
    bbAtoms["HCA"] = true;
    bbAtoms["HA"] = true;
    bbAtoms["1HA"] = true; // GLY
    bbAtoms["2HA"] = true; // GLY
    bbAtoms["C"] = true;
    bbAtoms["O"] = true;
    bbAtoms["O'"] = true;
    bbAtoms["O''"] = true;
    bbAtoms["OXT"] = true;
    bbAtoms["OT1"] = true;
    bbAtoms["OT2"] = true;
    bbAtoms["HC"] = true;
    bbAtoms["HOXT"] = true;
  }

  string atom_label = s;
  atom_label = scream_tools::strip_whitespace(atom_label);

  if (bbAtoms.find(atom_label) != bbAtoms.end()) {
    return true;
  }

  else {
    return false;
  }
}

bool scream_tools::is_heavy_atom(const string atomlabel) {

  if ( (scream_tools::strip_whitespace(atomlabel)).substr(0,1) != "H") {
    return true;
  } else {
    return false;
  }

}

bool scream_tools::is_N_term_hydrogen(const string atomLabel) {
  string label = strip_whitespace(atomLabel);
  // Below: Ideally, want to have a lookup table stored in a file than hardcoded like this.
  if (label.substr(0,1) != "H" and label.substr(1,1) != "H") {
    return false;
  }
  if ( label == "HN" or 
       label == "H NT" or
       label == "1H" or label == "2H" or label == "3H" or
       label.substr(0,2) == "HT" ) {
    return true;
  } else {
    return false;
  }
}

bool scream_tools::is_C_term_atom(const string atomLabel) {

  string label = strip_whitespace(atomLabel);
  if ( label == "OX" or
       label == "OXT" or
       label == "CT" or
       label == "HC") {
    return true;
  }
  else {
    return false;
  }

}

bool scream_tools::is_HCA_atom(const string atomLabel) {
  static map<string, bool> HCA_map;
  if (HCA_map.empty()) {
    HCA_map["HA"] = true;
    HCA_map["HCA"] = true;
    HCA_map["HA1"] = true;
    HCA_map["HA2"] = true;
    HCA_map["1HA"] = true;
    HCA_map["2HA"] = true;
  }
  
  if (HCA_map.find(atomLabel) == HCA_map.end()) {
    return false;
  } else {
    return true;
  }

}

bool scream_tools::is_natural_AA(const string resName) {
  
  string AA_name = scream_tools::strip_whitespace(resName);
  if (AA_name == "ALA" or 
      AA_name == "CYS" or 
      AA_name == "ASP" or 
      AA_name == "GLU" or 
      AA_name == "PHE" or 
      AA_name == "GLY" or 
      AA_name == "HIS" or 
      AA_name == "HSE" or 
      AA_name == "HSP" or 
      AA_name == "HSD" or
      AA_name == "ILE" or 
      AA_name == "LYS" or
      AA_name == "LEU" or 
      AA_name == "MET" or
      AA_name == "ASN" or 
      AA_name == "PRO" or
      AA_name == "GLN" or
      AA_name == "ARG" or
      AA_name == "SER" or
      AA_name == "THR" or
      AA_name == "VAL" or
      AA_name == "TRP" or
      AA_name == "TYR" or 
      AA_name == "APP" or
      AA_name == "GLP" or 
      AA_name == "ARN" or 
      AA_name == "LYN") {
    return true;
  } else {
    return false;
  }

}

bool scream_tools::AA_atom_order(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
  // N < HN < CA < HCA < C < O < CB < HCB < ... < OXT
  // returns true if a1 < a2.

  static map<string, double> value_map;
  // Backbone.
  if (value_map.empty()) {
    value_map[""] = 0;
    value_map["N"] = 1;
    value_map["NT"] = 1;
    value_map["HN"] = 2;
    value_map["HN1"] = 2.1;
    value_map["HN2"] = 2.2;
    value_map["HN3"] = 2.3;
    value_map["CA"] = 3;
    value_map["HCA"] = 4;
    value_map["HA1"] = 4.1;
    value_map["HA2"] = 4.2;
    value_map["C"] = 500;
    value_map["O"] = 600;
    
    // SideChain.
    value_map["CB"] = 11;
    value_map["HCB"] = 12;
    value_map["HB1"] = 12.1;
    value_map["HB2"] = 12.2;
    value_map["HB3"] = 12.3;
    
    value_map["CG"] = 21;
    value_map["HCG"] = 21.5;
    value_map["OG1"] = 22;
    value_map["HOG1"] = 22.5;
    value_map["CG1"] = 23;
    value_map["HCG1"] = 24;
    value_map["CG2"] = 25;
    value_map["HCG2"] = 26; // plus all those HG11, etc.  put these in later.
    value_map["SG"] = 27;
    value_map["HSG"] = 28;
    value_map["OG"] = 29;
    value_map["HOG"] = 30;
  
    value_map["CD"] = 31;
    value_map["HCD"] = 32;
    value_map["ND1"] = 33;
    value_map["HND1"] = 34;
    value_map["CD1"] = 35;
    value_map["HCD1"] = 36;
    value_map["CD2"] = 37;
    value_map["HCD2"] = 38;
    value_map["OD1"] = 39;
    value_map["HOD1"] = 39.1;
    value_map["OD2"] = 40;
    value_map["HOD2"] = 40.1;
    value_map["ND2"] = 41;
    value_map["HND2"] = 42;
    value_map["SD"] = 43;
    value_map["HSD"] = 43.1;
    
    value_map["CE"] = 51;
    value_map["HCE"] = 52;
    value_map["NE"] = 53;
    value_map["HNE"] = 54;
    value_map["NE1"] = 55;
    value_map["HNE1"] = 56;
    value_map["CE1"] = 57;
    value_map["HCE1"] = 58;
    value_map["OE1"] = 59;
    value_map["HOE1"] = 59.1;
    value_map["OE2"] = 60;
    value_map["HOE2"] = 60.1;
    value_map["NE2"] = 61;
    value_map["HNE2"] = 62;;
    value_map["CE2"] = 63;
    value_map["HCE2"] = 64;
    value_map["CE3"] = 65;
    value_map["HCE3"] = 66;
    
    value_map["CZ"] = 81;
    value_map["HCZ"] = 82;
    value_map["NZ"] = 83;
    value_map["HNZ"] = 84;
    value_map["CZ2"] = 85;
    value_map["HCZ2"] = 86;
    value_map["CZ3"] = 87;
    value_map["HCZ3"] = 88;
    
    value_map["CH2"] = 91;  // TRP
    value_map["HCH2"] = 92;
    value_map["OH"] = 93; // TYR
    value_map["HOH"] = 94;
    value_map["NH1"] = 95; // ARG
    value_map["HNH1"] = 96;
    value_map["HH11"] = 96.1; // ARN
    value_map["HH12"] = 96.2;
    value_map["NH2"] = 97;
    value_map["HNH2"] = 98;
    value_map["HH21"] = 98.1;
    value_map["HH22"] = 98.2;
    
    
    value_map["OX"] = 997;  
    value_map["HC"] = 998;
    value_map["OXT"] = 999;
    value_map["HOXT"] = 999.1;
    
    value_map["Na"] = 3001;
    value_map["Cl"] = 3002;
  }

  double a1_value=value_map[a1->stripped_atomLabel];
  double a2_value=value_map[a2->stripped_atomLabel];

  if (a1_value < a2_value)
    return true;
  else if ( a1_value > a2_value)
    return false;
  else
    return (a1->n < a2->n) ? true : false;
}

ScreamAtomV scream_tools::return_heavy_atoms(const ScreamAtomV& atom_list) {
  // Heavy atoms are atoms that do NOT start with an "H".  
  // Of course you could go wrong when you have HG, i.e. mercury.  But that doesn't happen too often.
  
  ScreamAtomV heavy_atoms;
  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    string atom_label = scream_tools::strip_whitespace((*itr)->atomLabel);
    if ( atom_label.substr(0, 1) != "H") {
      heavy_atoms.push_back(*itr);
    }
  }
  return heavy_atoms;

}

string scream_tools::get_AA_from_mutationInfo(const string mutInfo) {

  return mutInfo.substr(0,1);

}

int scream_tools::get_pstn_from_mutationInfo(const string mutInfo) {

  int num_end_pstn;
  
  for (num_end_pstn = 1; mutInfo.substr(num_end_pstn, 1) != "_"; ++num_end_pstn) {
  }

  string num_string = mutInfo.substr(1, num_end_pstn -1);
  return atoi(num_string.c_str());

}

string scream_tools::get_chn_from_mutationInfo(const string mutInfo) {

  return mutInfo.substr(mutInfo.size()-1, 1);

}

string scream_tools::three_letter_AA(const string AA_str) {

  map<string, string> AA123_map;
  AA123_map.insert(make_pair(string("A"), string("ALA")));
  AA123_map.insert(make_pair(string("B"), string("HSP")));
  AA123_map.insert(make_pair(string("C"), string("CYS")));
  AA123_map.insert(make_pair(string("D"), string("ASP")));
  AA123_map.insert(make_pair(string("E"), string("GLU")));
  AA123_map.insert(make_pair(string("F"), string("PHE")));
  AA123_map.insert(make_pair(string("G"), string("GLY")));
  AA123_map.insert(make_pair(string("H"), string("HIS")));
  AA123_map.insert(make_pair(string("I"), string("ILE")));
  AA123_map.insert(make_pair(string("J"), string("HSE")));
  AA123_map.insert(make_pair(string("K"), string("LYS")));
  AA123_map.insert(make_pair(string("L"), string("LEU")));
  AA123_map.insert(make_pair(string("M"), string("MET")));
  AA123_map.insert(make_pair(string("N"), string("ASN")));
  AA123_map.insert(make_pair(string("P"), string("PRO")));
  AA123_map.insert(make_pair(string("Q"), string("GLN")));
  AA123_map.insert(make_pair(string("R"), string("ARG")));
  AA123_map.insert(make_pair(string("S"), string("SER")));
  AA123_map.insert(make_pair(string("T"), string("THR")));
  AA123_map.insert(make_pair(string("V"), string("VAL")));
  AA123_map.insert(make_pair(string("W"), string("TRP")));
  AA123_map.insert(make_pair(string("Y"), string("TYR")));
  
  if (AA123_map.find(AA_str) != AA123_map.end()) {
    return AA123_map.find(AA_str)->second;
  }
  
  else {
    return AA_str;
  }

}

string scream_tools::one_letter_AA(const string AA_str) {

  map<string, string> AA123_map;
  AA123_map.insert(make_pair(string("A"), string("ALA")));
  AA123_map.insert(make_pair(string("B"), string("HSP")));
  AA123_map.insert(make_pair(string("C"), string("CYS")));
  AA123_map.insert(make_pair(string("D"), string("ASP")));
  AA123_map.insert(make_pair(string("E"), string("GLU")));
  AA123_map.insert(make_pair(string("F"), string("PHE")));
  AA123_map.insert(make_pair(string("G"), string("GLY")));
  AA123_map.insert(make_pair(string("H"), string("HIS")));
  AA123_map.insert(make_pair(string("I"), string("ILE")));
  AA123_map.insert(make_pair(string("J"), string("HSE")));
  AA123_map.insert(make_pair(string("K"), string("LYS")));
  AA123_map.insert(make_pair(string("L"), string("LEU")));
  AA123_map.insert(make_pair(string("M"), string("MET")));
  AA123_map.insert(make_pair(string("N"), string("ASN")));
  AA123_map.insert(make_pair(string("P"), string("PRO")));
  AA123_map.insert(make_pair(string("Q"), string("GLN")));
  AA123_map.insert(make_pair(string("R"), string("ARG")));
  AA123_map.insert(make_pair(string("S"), string("SER")));
  AA123_map.insert(make_pair(string("T"), string("THR")));
  AA123_map.insert(make_pair(string("V"), string("VAL")));
  AA123_map.insert(make_pair(string("W"), string("TRP")));
  AA123_map.insert(make_pair(string("Y"), string("TYR")));

  map<string, string> AA321_map;
  for (map<string, string>::const_iterator itr = AA123_map.begin(); itr != AA123_map.end(); ++itr) {
    AA321_map.insert(make_pair(itr->second, itr->first));
  }
  
  if (AA321_map.find(AA_str) != AA321_map.end()) {
    return AA321_map.find(AA_str)->second;
  } else {
    //return AA_str;
    return "X";
  }


}

multimap<string, SCREAM_ATOM*> scream_tools::deep_copy_str_SCREAM_ATOM_mm(const multimap<string, SCREAM_ATOM*>& in_mm, SCREAM_ATOM** head) {

  multimap<string, SCREAM_ATOM*> out_mm;
  multimap<string, SCREAM_ATOM*>::const_iterator ptr;

  int size_of_in_mm = in_mm.size();

  cout << "before new SCREAM_ATOM[]" << endl;
  //  cout << *head << "is head before assignment in deep copy" << endl;
  SCREAM_ATOM* head_to_new_SCREAM_array = new SCREAM_ATOM[size_of_in_mm];

  if (head != NULL) {
    *head = head_to_new_SCREAM_array;
  }

  int c = 0;

  for (ptr = in_mm.begin(); ptr != in_mm.end(); ++ptr) {
    //    head_to_new_SCREAM_array[c] = *(ptr->second);                           // using default assignment operator
    head_to_new_SCREAM_array[c].copy( *(ptr->second) ); // copy is default assignment operator; bonds not treated.
    out_mm.insert(make_pair(ptr->first, &(head_to_new_SCREAM_array[c])));
    ++c;
  }
  return out_mm;
}

void scream_tools::deep_copy_ScreamAtomV(const ScreamAtomV & in_atom_list, ScreamAtomV& new_atom_list) {

  if (new_atom_list.size() != 0) {
    cout << "new_atom_list.size() is not 0!  Exiting from scream_tools::deep_copy_ScreamAtomV. " << endl;
    exit(2);
  }
  
  map<SCREAM_ATOM*, SCREAM_ATOM*> old_to_new_atom_map;
  for (ScreamAtomVConstItr itr = in_atom_list.begin(); itr != in_atom_list.end(); ++itr) {
    SCREAM_ATOM* atom = new SCREAM_ATOM();
    atom->copy(*(*itr));
    new_atom_list.push_back(atom);
    old_to_new_atom_map[*itr] = atom;
  }
  // now update connectivities.
  
  for (ScreamAtomVConstItr new_atom_itr = new_atom_list.begin(); new_atom_itr != new_atom_list.end(); ++new_atom_itr) {
    
    map<SCREAM_ATOM*, int> new_connectivity_m;
    for (map<SCREAM_ATOM*, int>::const_iterator conn_itr = (*new_atom_itr)->connectivity_m.begin();
	 conn_itr != (*new_atom_itr)->connectivity_m.end(); ++conn_itr) {
      
      SCREAM_ATOM* old_atom = conn_itr->first;
      SCREAM_ATOM* new_atom = old_to_new_atom_map[old_atom];
      new_connectivity_m[new_atom] = conn_itr->second;
    }
    (*new_atom_itr)->connectivity_m.clear();
    (*new_atom_itr)->connectivity_m.insert(new_connectivity_m.begin(), new_connectivity_m.end());
      /*    for (map<SCREAM_ATOM*, int>::const_iterator conn_itr = new_connectivity_m.begin(); 
	 conn_itr != new_connectivity_m.end(); ++conn_itr) {
      
    } 
      */
  }

}


double scream_tools::calc_dihedral(const ScreamVector& V1, const ScreamVector& V2, const ScreamVector& V3, const ScreamVector& V4) {

  ScreamVector cross1 = (V1 - V2).cross(V2-V3);
  ScreamVector cross2 = (V4 - V3).cross(V2-V3);

  double angle = cross1.angleBtwn(cross2); // just angle between 2 vectors; need to append directionality.  returns in degrees.
  
  if (cross1.cross(V2-V3).dot(cross2) >= 0) {
    return angle;
  } else {
    return -angle;
  }
 
}
/*
int main() {

  double angle = scream_tools::calc_dihedral(ScreamVector(-1, 0, 0), ScreamVector(0, 0, 0), ScreamVector(0, 1, 0), ScreamVector(0, 1, 1));
  cout << "angle is " << angle * 180 / 3.1415926535 << endl;
  
}
*/
double scream_tools::calc_dihedral(const SCREAM_ATOM* const a1, const SCREAM_ATOM* const a2, const SCREAM_ATOM* const a3, const SCREAM_ATOM* const a4) {

  if (a1 == NULL or a2 == NULL or a3 == NULL or a4 == NULL) {
    cerr << "Cannot calculate dihedral angle because some atoms not defined" << endl;
    return 0;
  }

  return scream_tools::calc_dihedral(ScreamVector(a1->x[0], a1->x[1], a1->x[2]),
				     ScreamVector(a2->x[0], a2->x[1], a2->x[2]),
				     ScreamVector(a3->x[0], a3->x[1], a3->x[2]),
				     ScreamVector(a4->x[0], a4->x[1], a4->x[2]));

}


void scream_tools::translation(SCREAM_ATOM* const atom, const ScreamVector& V) {

  for (int i = 0; i <= 2; ++i) {
    atom->x[i] += V[i];
  }

}

void scream_tools::translation(const vector<SCREAM_ATOM*>& atom_v, const ScreamVector& V) {

  vector<SCREAM_ATOM*>::const_iterator itr;
  for (itr = atom_v.begin(); itr != atom_v.end(); ++itr) {
    for (int i = 0; i<= 2; ++i) {
      (*itr)->x[i] += V[i];
    }
  }

}

void scream_tools::translation(const multimap<string, SCREAM_ATOM*>& atom_mm, const ScreamVector& v) {

  multimap<string, SCREAM_ATOM*>::const_iterator itr;

  for (itr = atom_mm.begin(); itr != atom_mm.end(); ++itr) {
    for (int i = 0; i <= 2; ++i) {
      itr->second->x[i] += v[i];
    }
  }

}

pair<ScreamMatrix, ScreamVector> scream_tools::getTransformationPairFromAtoms(const ScreamAtomV& stationary, const ScreamAtomV& toBeMoved, bool exactMatch) {
  /* This functions returns a ScreamMatrix object that does the linear transformation that best matches the toBeMoved atoms to the stationary atoms. */
  /* If exactMatch is true, only a translation matrix will be returned. */

  ScreamMatrix M(ScreamVector(1,0,0),
		 ScreamVector(0,1,0),
		 ScreamVector(0,0,1));
  ScreamVector V;

  if (exactMatch)
    {
      /* M already is the Identity matrix */

      SCREAM_ATOM* s = *(stationary.begin() );
      SCREAM_ATOM* m = *(toBeMoved.begin() );
      
      ScreamVector s_pstn(s->x[0], s->x[1], s->x[2]);
      ScreamVector m_pstn(m->x[0], m->x[1], m->x[2]);

      V = s_pstn - m_pstn;
    }

  else 
    {
      // not implemented yet; will do eventually.
    }

  return make_pair(M, V);

}


double scream_tools::distance(const ScreamAtomV& list1, const ScreamAtomV& list2) {
  /* Calculates distance (Euclidean, CRMS) between two lists of atoms. */
  /* Remark: if list1.size() != list2.size() issues error, proceeds to compare just the first n atom CRMS's.*/

  int list1_size = list1.size();
  int list2_size = list2.size();
  int min_size = (list1_size < list2_size) ? list1_size : list2_size;

  if (list1_size == 0 or list2_size == 0) {
    cerr << "List size is zero! Exiting." << endl;
    exit(2);
    return -999;
  }

  if (list1_size != list2_size) {
    cerr << "Warning: Comparing lists of different sizes, in scream_tools::distance.  Proceed." << endl;
  }

  SCREAM_ATOM* atom1, *atom2;
  double CSquared = 0;

  for (int i = 0; i < min_size; ++i) {  // min_size: number of elements
    atom1 = list1[i];
    atom2 = list2[i];
    CSquared += atom1->distance_squared(atom2);
  }
  double CMS = CSquared / min_size;
  double CRMS = sqrt(CMS);

  return CRMS;
}


double scream_tools::distance_by_atom_label_map(const map< string, SCREAM_ATOM* >& map1, const map< string, SCREAM_ATOM* >& map2) {
  int map1_size = map1.size();
  int map2_size = map2.size();

  int min_size = (map1_size < map2_size) ? map1_size : map2_size;

  if (map1_size == 0 or map2_size == 0) {
    cerr << "Map size is zero! Exiting." << endl;
    exit(2);
    return -999;
  }
  

  if (map1_size != map2_size) {
    cerr << "Warning: Compariing maps of idfferent sizes, in scream_tools::distance_by_atom_label_map. Exiting." << endl;
    exit(2);
    return -999;
  }

  SCREAM_ATOM* atom1, *atom2;
  double CSquared = 0;
  for (map<string, SCREAM_ATOM*>::const_iterator itr = map1.begin(); itr != map1.end(); ++itr) {
    string atom_label = strip_whitespace( itr->first);
    atom1 = itr->second;
    atom2 = (map2.find(atom_label))->second;

    CSquared += atom1->distance_squared(atom2);

  }

  double CMS = CSquared / min_size;
  double CRMS = sqrt(CMS);

  return CRMS;

}


pair<double, string> scream_tools::max_equivalent_atom_dist(const map< string, SCREAM_ATOM* >& map1 , const map< string, SCREAM_ATOM* >& map2) {
    int map1_size = map1.size();
  int map2_size = map2.size();

  int min_size = (map1_size < map2_size) ? map1_size : map2_size;

  if (map1_size == 0 or map2_size == 0) {
    cerr << "Map size is zero! Exiting." << endl;
    exit(2);
    return make_pair(-999, " ");
  }
  

  if (map1_size != map2_size) {
    cerr << "Warning: Compariing maps of idfferent sizes, in scream_tools::distance_by_atom_label_map. Exiting." << endl;
    exit(2);
    return make_pair(-999, " ");
  }

  SCREAM_ATOM* atom1, *atom2;
  string Max_Dist_Atom_Label;
  double Max_Dist = -1;
  for (map<string, SCREAM_ATOM*>::const_iterator itr = map1.begin(); itr != map1.end(); ++itr) {
    string atom_label = strip_whitespace( itr->first);
    atom1 =  itr->second;
    atom2 = (map2.find(atom_label))->second;
    double dist = atom1->distance_squared(atom2);
    if (dist > Max_Dist) {
      Max_Dist = dist;
      Max_Dist_Atom_Label = atom_label;
    }

  }

  Max_Dist = sqrt(Max_Dist);
  return make_pair(Max_Dist, Max_Dist_Atom_Label);


}


void scream_tools::calc_new_HN_atom_coords(const SCREAM_ATOM* CA_i, const SCREAM_ATOM* N_i, const SCREAM_ATOM* C_i_minus_1, SCREAM_ATOM* HN_atom) {

  /* This functions allocates memory for a new SCREAM_ATOM*.  This needs to be deleted in a destructor elsewhere to maintain conservation of memory. */
  
  ScreamVector CA_i_v(CA_i);
  ScreamVector N_i_v(N_i);
  ScreamVector C_i_minus_1_v(C_i_minus_1);
  
  // now make coords for HN; assumes the angle CA-N-C is close to 120.  The positioning of HN is approximate.
  ScreamVector N_CA = (CA_i_v - N_i_v).normalizedVector();
  ScreamVector N_C  = (C_i_minus_1_v - N_i_v).normalizedVector();

  ScreamVector N_HN = (ScreamVector(0,0,0) - N_CA - N_C).normalizedVector();
  N_HN = N_HN * 0.97; // distance of HN-N bond: take 0.97A as bond length.
  ScreamVector HN = N_i_v + N_HN;  // position of HN atom.

  // Now make new HN atom.  bgf types assumed.  in the future, will make typing automatic, reading in a dreidii or CHARM typing file.
  
  HN_atom->x[0] = HN[0];
  HN_atom->x[1] = HN[1];
  HN_atom->x[2] = HN[2];
  
}


void scream_tools::make_connectivity(ScreamAtomV& atom_list, const map< int, vector<int> > & connect_map) {
  /* Assumes that atom_list has corresponding atom numbering information to connect_map. Also, erases previous connecitivity records.*/
  
  //  cout << "Entering make_connectivity" << endl;

  // First make a map<int, SCREAM_ATOM*> to make access easier. */
  map<int, SCREAM_ATOM*> n_to_atom_map;
  for (ScreamAtomVConstItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    int n = (*itr)->n;
    SCREAM_ATOM* a = *itr;
    n_to_atom_map[n] = a;
    //a->dump();
  }
  
  // Then populate atom->connectivity_m.
  for (map< int, vector<int> >::const_iterator itr = connect_map.begin(); itr != connect_map.end(); ++itr) {
    int base_n = itr->first;
    SCREAM_ATOM* base_atom = n_to_atom_map[base_n];

    //base_atom->dump();

    //    base_atom->connectivity_m.clear();
    for (vector<int>::const_iterator itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {
      int connect_n = *itr2;
      SCREAM_ATOM* connected_atom = n_to_atom_map[connect_n];
      if (connected_atom != NULL) {   // Possible that the connecitivity map contains atoms that aren't included in the atom_list passed in.
	base_atom->connectivity_m[connected_atom] = 1; // order assumed to be 1.
	//cout << "base_atom: " << base_atom->n << "  connected atom: " << connected_atom->n << endl;
      }

    }
  }
  
}


int scream_tools::should_exclude_on_1_2(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
  // returns 1 (true) if the two atoms are connected, i.e. should exclude based on 1-2.
  int should_exclude = 0;
  for (map<SCREAM_ATOM*, int>::const_iterator itr = a1->connectivity_m.begin();
       itr != a1->connectivity_m.end(); ++itr) {
    SCREAM_ATOM* connected_atom = (*itr).first;
    if (a2 == connected_atom) {
      should_exclude = 1;
      break;
    }
  }
  return should_exclude;
}

int scream_tools::should_exclude_on_1_3(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
  // should exclude based on 1-3.  one level of search tree deeper than should_exclude_on_1_2.
  int should_exclude = 0;
  for (map<SCREAM_ATOM*, int>::const_iterator itr = a1->connectivity_m.begin();
       itr != a1->connectivity_m.end(); ++itr) {
    SCREAM_ATOM* connected_atom = (*itr).first;
    if (scream_tools::should_exclude_on_1_2(connected_atom, a2) ) {   // call exclude on 1_2.
      should_exclude = 1;
      break;
    }
  }
  return should_exclude;
}

int scream_tools::should_exclude_on_1_2_or_1_3(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
  int should_exclude =0;
  for (map<SCREAM_ATOM*, int>::const_iterator itr = a1->connectivity_m.begin();
       itr != a1->connectivity_m.end(); ++itr) {
    SCREAM_ATOM* connected_atom = (*itr).first;
    // WORKING ON THIS!
  }
  return should_exclude;

}

vector<ScreamVector> scream_tools::generateHydrogenCoords(SCREAM_ATOM* base) {
  Debug debugInfo("scream_tools::generateHydrogenCoords(SCREAM_ATOM* base)");

  vector<ScreamVector> generatedCoords; generatedCoords.clear();
  string ff_type = base->getAtomType();
  // C_3n, if > 4 connected atoms, quit
  int connectedAtoms = base->connectivity_m.size();
  if (connectedAtoms >= 4) {
    return generatedCoords;
  }
  // C_22, C_R1.  if >= 3 connected atoms, quit
  if ( (ff_type[2] == 'R' or ff_type[2] == '2') and connectedAtoms >= 3) {
    return generatedCoords;
  }

  map<SCREAM_ATOM*, int>::iterator baseConnectedItr = base->connectivity_m.begin();
  
  SCREAM_ATOM* connectedAtom1 = baseConnectedItr->first;
  SCREAM_ATOM* connectedAtom2(NULL);
  SCREAM_ATOM* connectedAtom3(NULL);
  SCREAM_ATOM* connected13Atom(NULL); // 1-3 atom: used so that generated sp3 hydrogens are staggers, not eclipsed.
  
  ScreamVector baseV(base);
  ScreamVector connectedAtom1V(connectedAtom1);

  // Real work--calculating hydrogen positions.
  if ( ff_type[2] == '3') {
    if (connectedAtoms == 1) {
      // if terminal heavy atom: two cases: sp3 - sp3, and sp3 - sp2.  either way, make sure staggered.
      // find an appropriate connected13Atom--just pick the first atom, heavy or hydrogen.
      map<SCREAM_ATOM*, int>::iterator connectedAtom1Itr = connectedAtom1->connectivity_m.begin();
      while (connectedAtom1Itr->first == base) {
	connectedAtom1Itr++;
      }
      connected13Atom = connectedAtom1Itr->first;

      ScreamVector connected13AtomV(connected13Atom);
      debugInfo.out(" Calling generate3SP3HydrogenCoords");
      generatedCoords = scream_tools::generate3SP3HydrogenCoords(baseV,connectedAtom1V, connected13AtomV);
      if (ff_type[0] == 'O' or ff_type[0] == 'S') {
	generatedCoords.pop_back();
	generatedCoords.pop_back();
      }
    }
    else if (connectedAtoms == 2) {
      // then generate 2 hydrogens.  middle of an aliphatic chain.
      ++baseConnectedItr;
      connectedAtom2 = baseConnectedItr->first;
      ScreamVector connectedAtom2V(connectedAtom2);
      debugInfo.out(" Calling generate2SP3HydrogenCoords");
      generatedCoords = scream_tools::generate2SP3HydrogenCoords(baseV, connectedAtom1V, connectedAtom2V);
      if (ff_type[0] == 'O' or ff_type[0] == 'S') {
	generatedCoords.pop_back();
	generatedCoords.pop_back();
      }
    }
    else if (connectedAtoms == 3) {
      // PRO.  maybe something else too.
      connectedAtom2 = baseConnectedItr->first;
      ++baseConnectedItr;
      connectedAtom3 = baseConnectedItr->first;
      ScreamVector connectedAtom2V(connectedAtom2);
      ScreamVector connectedAtom3V(connectedAtom3);
      debugInfo.out(" Calling generate1SP3HydrogenCoords");
      generatedCoords = scream_tools::generate1SP3HydrogenCoords(baseV, connectedAtom1V, connectedAtom2V, connectedAtom3V);
    }
  }
  else if ( ff_type[2] == 'R' or ff_type[2] == '2') {
    if (ff_type[0] == 'O' or ff_type[0] == 'S') {
      return generatedCoords;
    }

    // two cases: sp2 - sp3, and sp2 - sp2.
    // don't think sp2-sp3 exists for proteins.
    if (connectedAtoms == 1) {
      // like amine in ASN, GLN.  NH2 - C == 0 -- sp2-sp2.  is sp2-sp3, taking 0 degree (hydrogens added in plane) is the DREIDING minimum anyway.
      map<SCREAM_ATOM*, int>::iterator connectedAtom1Itr = connectedAtom1->connectivity_m.begin();
      while (connectedAtom1Itr->first == base) {
	++connectedAtom1Itr;
      }
      connected13Atom = connectedAtom1Itr->first;

      ScreamVector connected13AtomV(connected13Atom);
      debugInfo.out(" Calling generate2SP2HydrogenCoords");
      generatedCoords = scream_tools::generate2SP2HydrogenCoords(baseV, connectedAtom1V, connected13AtomV);

    } 
    else if (connectedAtoms == 2) {
      // like backbone NH.
      ++baseConnectedItr;
      connectedAtom2 = baseConnectedItr->first;
      ScreamVector connectedAtom2V(connectedAtom2);
      debugInfo.out(" Calling generate1SP2HydrogenCoords");
      generatedCoords = scream_tools::generate1SP2HydrogenCoords(baseV, connectedAtom1V, connectedAtom2V);

    }
  }
  return generatedCoords;

}

ScreamAtomV scream_tools::createHydrogens(vector<ScreamVector>& H_coords, SCREAM_ATOM* base) {

  /* Fields to add in H:
   * 0. KEYWORD ATOM.
   * 1. H_ or H___A.  DREIDING ff type.
   * 2. Atom name.  HCB, or HB1/HB2/HB3?  
   * 3. Residue name, from base atom.
   * 4. Chain name, from base atom.
   * 5. Residue position, from base atom.
   * 6. Lone pair = 0.
   * 7. X,Y,Z, from H_coords.
   * 8. connectivity: to base atom.
   * Fields not to worry about:
   * 1. atom number.  take care of this in higher level routine.
   * 2. charges: high level routine.
   */

  string baseAtomLabel = scream_tools::strip_whitespace(base->getAtomLabel());
  // below: new hydrogen name.
  string H_dreidiiType = "H___A";
  if (baseAtomLabel[0] == 'C') {
    H_dreidiiType = "H_";
  }
  string H_atomLabel = "H" + baseAtomLabel;
  
  
  ScreamAtomV newHs;

  for (vector<ScreamVector>::iterator itr = H_coords.begin(); itr != H_coords.end(); ++itr) {
    SCREAM_ATOM* newH = new SCREAM_ATOM();
    newH->keyw = "ATOM";
    newH->setAtomType(H_dreidiiType);
    newH->setAtomLabel(H_atomLabel);
    newH->setResName(base->resName);
    newH->setChain(base->chain);
    newH->resNum = base->resNum;
    // lone pair in default constructor
    newH->setX( (*itr)[0] );
    newH->setY( (*itr)[1] );
    newH->setZ( (*itr)[2] );

    newH->make_bond(base);

    newHs.push_back(newH);

  }

  return newHs;

}

vector<ScreamVector> scream_tools::generate3SP3HydrogenCoords(ScreamVector& base, ScreamVector& connected, ScreamVector& samePlane) {
  // Uses HC or NH or OH or SH distance of 1.00Angstrom.
  /* Picture:
   CONN  (b)   BASE         (b): base_connected_axis
         --->               (i): in_plane_vector
  (o)|  ^    \ (h)          (h): in_plane_hydrogen
     | /      \             (o): in_plane_orthogonal_vector
     V/ (i)    V

  SAMEPLANE

  */

  ScreamVector base_connected_axis = (base - connected).normalizedVector();
  ScreamVector in_plane_vector = (samePlane - connected).normalizedVector();
  ScreamVector in_plane_orthogonal_vector = (in_plane_vector - base_connected_axis * base_connected_axis.dot(in_plane_vector) ).normalizedVector();

  // 0.35355: arc tan 19.47, which is 109.47 - 90 degrees.
  ScreamVector in_plane_hydrogen = base + (base_connected_axis * 0.35355 +  in_plane_orthogonal_vector).normalizedVector();
  
  // now figure out staggered hydrogen positions.  
  ScreamMatrix dummy;
  ScreamMatrix R = dummy.rotAboutV(base_connected_axis, 60);
  
  ScreamVector H1 = R * (in_plane_hydrogen -base) + base ;
  ScreamVector H2 = R * R * (H1 -base) + base;
  ScreamVector H3 = R * R * (H2 -base) + base;

  vector<ScreamVector> new_H_coords;
  new_H_coords.push_back(H1);
  new_H_coords.push_back(H2);
  new_H_coords.push_back(H3);

  return new_H_coords;

}

vector<ScreamVector> scream_tools::generate2SP3HydrogenCoords(ScreamVector& base, ScreamVector& connected1, ScreamVector& connected2) {
  /* Picture:
     Connected1                         mdp: midBasePoint
        \                                 X; base
	 \
    (mbp) \
	---X
	  /
         /
        /
     Connected2

     
  */
  
  ScreamVector normalizedConnected1V = base + (connected1 - base).normalizedVector();
  ScreamVector normalizedConnected2V = base + (connected2 - base).normalizedVector();

  ScreamVector midBasePoint = normalizedConnected1V + (normalizedConnected2V - normalizedConnected1V)/2;
  ScreamVector midBaseAxis = (base - midBasePoint).normalizedVector();
  
  ScreamVector connected2Vector = (connected2 - midBasePoint).normalizedVector();
  ScreamVector orthogonalVector = (connected2Vector - midBaseAxis * midBaseAxis.dot(connected2Vector) ).normalizedVector();
  
  // 0.5 and 0.7071: appropriate ratios for generating the tetrahedral.
  ScreamVector in_plane_hydrogen = base + (midBaseAxis * 0.5  + orthogonalVector * 0.7071).normalizedVector();
  // now figure out staggered hydrogen positions.
  ScreamMatrix dummy;
  ScreamMatrix R = dummy.rotAboutV(midBaseAxis, 90); // 90 degrees: located along points of a cube.
  
  ScreamVector H1 = R * (in_plane_hydrogen - base) + base;
  ScreamVector H2 = R * R * (H1 - base) + base;

  
  vector<ScreamVector> new_H_coords;
  new_H_coords.push_back(H1);
  new_H_coords.push_back(H2);

  return new_H_coords;


}

vector<ScreamVector> scream_tools::generate1SP3HydrogenCoords(ScreamVector& base, ScreamVector& connected1, ScreamVector& connected2, ScreamVector& connected3) {
  // No pictures needed, I assume?

  ScreamVector normalizedConnected1V = (connected1 - base).normalizedVector();
  ScreamVector normalizedConnected2V = (connected2 - base).normalizedVector();
  ScreamVector normalizedConnected3V = (connected3 - base).normalizedVector();
  
  ScreamVector oppositeDirection = (normalizedConnected1V + normalizedConnected2V + normalizedConnected3V).normalizedVector();

  ScreamVector H1 = base - oppositeDirection;

  H1.printMe();

  vector<ScreamVector> new_H_coords;
  new_H_coords.push_back(H1);

  return new_H_coords;

}

vector<ScreamVector> scream_tools::generate2SP2HydrogenCoords(ScreamVector& base, ScreamVector& connected, ScreamVector& samePlane) {
  // First part of this is same code as generate3SP3HydrogenCoords.
  // For sp2-sp3 and sp2-sp2 situations the same--put hydrogens in plane with one of the connected atoms.  for sp2-sp2, correct.  for sp2-sp3, this leads to one of 3 possible structures, each of which should be consistent with the DREIDING minima.  Ideally, should pick a heavy atom as the samePlane atom.
//   cout << " ??? generate2SP2HydrogenCoords " << endl;
//   base.printMe();
//   connected.printMe();
//   samePlane.printMe();
//   cout << " ??? over generate2SP2HydrogenCoords" << endl;

  ScreamVector base_connected_axis = (base - connected).normalizedVector();
  ScreamVector in_plane_vector = (samePlane - connected).normalizedVector();
  ScreamVector in_plane_orthogonal_vector = (in_plane_vector - base_connected_axis * base_connected_axis.dot(in_plane_vector) ).normalizedVector();

  // 0.57735: tan 30, which is 120 - 90 degrees.
  ScreamVector in_plane_hydrogen1 = base + (base_connected_axis * 0.57735 +  in_plane_orthogonal_vector).normalizedVector();
  ScreamVector in_plane_hydrogen2 = base + (base_connected_axis * 0.57735 -  in_plane_orthogonal_vector).normalizedVector();
  

  vector<ScreamVector> new_H_coords;
  new_H_coords.push_back(in_plane_hydrogen1);
  new_H_coords.push_back(in_plane_hydrogen2);


  return new_H_coords;

}

vector<ScreamVector> scream_tools::generate1SP2HydrogenCoords(ScreamVector& base, ScreamVector& connected1, ScreamVector& connected2) {
  // This is straightforward.  Take the two connected atom vectors and just add them to the base atom in the opposite direction, normalize.

  ScreamVector normalizedConnected1V = (connected1 - base).normalizedVector();
  ScreamVector normalizedConnected2V = (connected2 - base).normalizedVector();

  ScreamVector oppositeDirection = (normalizedConnected1V + normalizedConnected2V).normalizedVector();

  ScreamVector H1 = base - oppositeDirection;

  vector<ScreamVector> new_H_coords;
  new_H_coords.push_back(H1);

  return new_H_coords;

}




//AtomNotFoundException::AtomNotFoundException(string
