#include "defs.hpp"
#include "scream_E_functionals_vdw.hpp"
#include "scream_tools.hpp"
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

VDW_fields::VDW_fields() {

}

VDW_fields::VDW_fields(double RNB, double DENB, double SCALE) {

  this->RNB = RNB;
  this->DENB = DENB;
  this->SCALE = SCALE;

}

VDW_fields::~VDW_fields() {

}

VDW_fields& VDW_fields::operator=(const VDW_fields& in) {

  if (this == &in) return *this;

  this->RNB = in.RNB;
  this->DENB = in.DENB;
  this->SCALE = in.SCALE;
  return *this;

}

VDW_delta_fields::VDW_delta_fields() {

}

VDW_delta_fields::VDW_delta_fields(double mu, double sigma) {

  this->mu = mu;
  this->sigma = sigma;

}

VDW_delta_fields::~VDW_delta_fields() {

}




SCREAM_VDW_OBJ::SCREAM_VDW_OBJ() {

}


SCREAM_VDW_OBJ::~SCREAM_VDW_OBJ() {
  for (map<string, VDW_fields*>::iterator itr = this->vdw_dict.begin();
       itr != this->vdw_dict.end(); ++itr ) {
    delete (*itr).second;

  }
  for (map<string, map<AtomResInfo, VDW_delta_fields*>* >::iterator itr = this->vdw_delta_library_dict.begin();
       itr != this->vdw_delta_library_dict.begin(); ++itr) {
    for (map<AtomResInfo, VDW_delta_fields*>::iterator ii = itr->second->begin(); 
	 ii != itr->second->end(); ++ii) {

      delete ii->second;
    }
    delete itr->second;
  }

}

void SCREAM_VDW_OBJ::read_param_line(string line) {
  vector<string> fields;
  split(line, string(" "), fields);
  
  //      cout << line << endl;
  if (fields[0] == "VDW") { return; };
  double RNB = atof(fields[2].c_str());
  double DENB = atof(fields[3].c_str());
  double SCALE = atof(fields[4].c_str());
  this->vdw_dict[fields[0]] = new VDW_fields(RNB, DENB, SCALE); // fields[0]: ff_type
  
  if (fields[0] == "AUTOTYPE") return;

}

void SCREAM_VDW_OBJ::read_param_file(string ff_file) {

  ifstream FF_FILE;
  FF_FILE.open(ff_file.c_str());

  if (!FF_FILE.good()) {
    cerr << "Unable to open forcefield file: " << ff_file << endl;
    exit(8);
  }

  string line;
  string read_mode = "";


  while (!FF_FILE.eof()) {
    
    char line_ch[256];
    FF_FILE.getline(line_ch, sizeof(line_ch));
    line = string(line_ch);

    if (line.substr(0,1) == "*" or line.substr(0,1) == "#" or line.substr(0,1) == "!") continue;

    /* Parameter file keywords:
     * PARAMETER FORMET 
     * FORCEFIELD
     * DEFAULTS
     * RNB GEOMN
     * SCAL NB14
     * LCOULMB  
     * R*EPS    
     * DIELCTRIC
     * LHBOND   
     * USRENERGY
     * FFLABEL
     * ADDED H
     * LONE PAIRS
     * DEL RE ATOMS
     * GASTEIGER
     * VDW
     * AUTOTYPE
     * NONBOND-OFF
     * BONDSTRTCH 
     * ANGLE-(L-C-R)
     * TORSION
     * ...
     */
    vector<string> fields;
    split(line, string(" "), fields);
    if (fields[0] == "VDW") read_mode = "VDW";
    if (read_mode == "VDW") {
      //      cout << line << endl;
      if (fields[0] == "VDW") { continue; };
      double RNB = atof(fields[2].c_str());
      double DENB = atof(fields[3].c_str());
      double SCALE = atof(fields[4].c_str());
      this->vdw_dict[fields[0]] = new VDW_fields(RNB, DENB, SCALE); // fields[0]: ff_type
      
    }
    if (fields[0] == "AUTOTYPE") break;
  }
  FF_FILE.close();
}


void SCREAM_VDW_OBJ::read_Scream_delta_file(std::string file) {
  /* Format of Scream Neighborhood delta file: Library Type, Residue Name, AtomLabel, Delta value  */
  Debug debugInfo("SCREAM_VDW_OBJ::read_Scream_delta_file(std::string file) ");
  debugInfo.out("SCREAM DELTA FILENAME IS : " + file );
  ifstream FILE;
  FILE.open(file.c_str());

  if (!FILE.good()) {
    cerr << "Unable to open SCREAM delta file: " << file << endl;
    exit(8);
  }


  char line_ch[256];
  string library_type, resName, atomLabel;
  double atom_mu, atom_sigma;

  this->vdw_delta_library_dict.clear();  
  while (!FILE.eof()) {

    FILE.getline(line_ch, sizeof(line_ch));
    string line(line_ch);

    if (line == "" or int(line[0]) == 0) continue; 
    if (line.substr(0,1) == "*" or line.substr(0,1) == "#" or line.substr(0,1) == "!" or int(line[0]) == 0) continue;

    stringstream ss(line);

    ss >> library_type >> resName >> atomLabel >> atom_mu >> atom_sigma;

    //    vector<string> fields;
    //    split(line, fields);

//     string library_type = fields[0];
//     string resName = fields[1];
//     string atomLabel = fields[2];
//     double atom_mu = atof(fields[3].c_str());
//     double atom_sigma = atof(fields[4].c_str());

    // init one element of the map
    if (this->vdw_delta_library_dict.count(library_type) == 0) {
      map<AtomResInfo, VDW_delta_fields*> * new_aI_vdw_field_map = new map<AtomResInfo, VDW_delta_fields*>;
      new_aI_vdw_field_map->clear();
      this->vdw_delta_library_dict[library_type] = new_aI_vdw_field_map;
    }

    map<AtomResInfo, VDW_delta_fields*> * aI_vdw_map = this->vdw_delta_library_dict[library_type];
    AtomResInfo aRI(resName, atomLabel);

    if (aI_vdw_map->find(aRI) != aI_vdw_map->end()) { // value initialized earlier already; overwrite value and put in new values.
      delete (*aI_vdw_map)[aRI];
      VDW_delta_fields* vDF = new VDW_delta_fields(atom_mu, atom_sigma);
      aI_vdw_map->insert(make_pair(aRI, vDF) );
    }
    else {
      VDW_delta_fields* vDF = new VDW_delta_fields(atom_mu, atom_sigma);
      aI_vdw_map->insert(make_pair(aRI, vDF) );
    }



    //cout << this->vdw_delta_library_dict[library_type][AtomResInfo(resName, atomLabel)]->delta << endl;

  // init the default value for library_type NA, or empty string.  7-27-07: why do i want to do the below?  if NA, user should fix problem.
  //  assert( this->vdw_delta_library_dict.size() != 0);
  //  map<string, map<AtomResInfo, VDW_delta_fields*> >::iterator sample = this->vdw_delta_library_dict.begin();
  
  //  map<AtomResInfo, VDW_delta_fields*> new_aI_vdw_field_map;
  //  this->vdw_delta_library_dict[""] = new_aI_vdw_field_map;

  //  for (map<AtomResInfo, VDW_delta_fields*>::iterator itr = sample->second.begin(); 
  //       itr != sample->second.end(); ++itr ) {
  //    AtomResInfo arI = (*itr).first;
  //    this->vdw_delta_library_dict[""][arI] = new VDW_delta_fields(0,0); // init to 0
    // }
  }
  
  FILE.close();

//   for (map<string, map<AtomResInfo, VDW_delta_fields*> >::iterator iii = this->vdw_delta_library_dict.begin();
//        iii != this->vdw_delta_library_dict.end(); ++iii) {
//     cout << ":::" << iii->first << ":::" << endl;
//     for (map<AtomResInfo, VDW_delta_fields*>::iterator itr = iii->second.begin(); 
// 	 itr != iii->second.end(); ++itr) {
//       cout << itr->first << endl;
//       cout << itr->second->delta << endl;
//     }
//   }
  
 

}


double SCREAM_VDW_OBJ::get_RNB(string ff_type) {
  //  ff_type = scream_tools::strip_whitespace(ff_type);  // too slow
  map<string, VDW_fields*>::iterator itr = this->vdw_dict.find(ff_type);
  if (itr == this->vdw_dict.end()) {
    cerr << "ff_type " << ff_type << " not found in parameter file!  Exiting" << endl;
    exit(0);
  }

  return this->vdw_dict[ff_type]->RNB;
}

double SCREAM_VDW_OBJ::get_DENB(string ff_type) {
  // ff_type = scream_tools::strip_whitespace(ff_type); // too slow
  map<string, VDW_fields*>::iterator itr = this->vdw_dict.find(ff_type);
  if (itr == this->vdw_dict.end()) {
    cerr << "ff_type " << ff_type << " not found in parameter file!  Exiting" << endl;
    exit(0);
  }

  return this->vdw_dict[ff_type]->DENB;
}

double SCREAM_VDW_OBJ::get_SCALE(string ff_type) {
  // ff_type = scream_tools::strip_whitespace(ff_type); // too slow
  map<string, VDW_fields*>::iterator itr = this->vdw_dict.find(ff_type);
  if (itr == this->vdw_dict.end()) {
    cerr << "ff_type " << ff_type << " not found in parameter file!  Exiting" << endl;
    exit(0);
  }

  return this->vdw_dict[ff_type]->SCALE;
}

VDW_fields* SCREAM_VDW_OBJ::get_VDW_fields(string ff_type) {
  map<string, VDW_fields*>::iterator itr = this->vdw_dict.find(ff_type);
  if (itr == this->vdw_dict.end()) {
    cerr << "ff_type " << ff_type << " not found in parameter file!  Exiting" << endl;
    exit(0);
  }

  return (*itr).second;
}

VDW_delta_fields* SCREAM_VDW_OBJ::get_VDW_delta_fields(string library_name, AtomResInfo aRI) const {
  // Returns NULL if can't find a match.
  
  map<string, map<AtomResInfo, VDW_delta_fields*>* >::const_iterator itr = this->vdw_delta_library_dict.find(library_name);
  if (itr == this->vdw_delta_library_dict.end()) {
    cerr << "library_name " << library_name << " not found in SCREAM delta file!  Exiting" << endl;
    exit(0);
  }

//   cout << "size of VDW_delta_fields " << itr->second.size() << endl;
//   for (map<AtomResInfo, VDW_delta_fields*>::const_iterator ii = itr->second.begin(); ii != itr->second.end(); ++ii) {
//     cout << ii->first.resName << " " << ii->first.atomLabel << endl;
//   }

  
  map<AtomResInfo, VDW_delta_fields*>::const_iterator itr2 = itr->second->find(aRI);
  if (itr2 == itr->second->end()) {
    cerr << "SCREAM delta value doesn't exist for this sidechain: " << aRI.resName << " " << aRI.atomLabel << endl;
    return NULL;
  }
  
  return itr2->second;


}



double SCREAM_VDW_OBJ::calc_VDW_6_8(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
  /* old, don't call this function */
  string a1_type = a1->stripped_atomType;
  string a2_type = a2->stripped_atomType;

  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }

  //  double RNB = this->_geom_mean(this->get_RNB(a1_type), this->get_RNB(a2_type));
  //  double DENB = this->_geom_mean(this->get_DENB(a1_type), this->get_DENB(a2_type));
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double E =  DENB * ( 3 * pow(rho,-8) - 4 * pow(rho,-6) );
  return E;

}

double SCREAM_VDW_OBJ::calc_VDW_6_9(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
 /* 6-9 Functional form: 
   * DENB * { rho^(-9) - 2 * rho^(-6) }
   */

  string a1_type = a1->stripped_atomType;
  string a2_type = a2->stripped_atomType;

  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }

  //  double RNB = this->_geom_mean(this->get_RNB(a1_type), this->get_RNB(a2_type));
  //  double DENB = this->_geom_mean(this->get_DENB(a1_type), this->get_DENB(a2_type));
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double E =  DENB * ( 2 * pow(rho,-9) - 3 * pow(rho,-6) );
  return E;


}


double SCREAM_VDW_OBJ::calc_VDW_6_10(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
 /* 6-9 Functional form: 
   * DENB * { 1.5 * rho^(-10) - 2.5 * rho^(-6) }
   */

  string a1_type = a1->stripped_atomType;
  string a2_type = a2->stripped_atomType;
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }

  //  double RNB = this->_geom_mean(this->get_RNB(a1_type), this->get_RNB(a2_type));
  //  double DENB = this->_geom_mean(this->get_DENB(a1_type), this->get_DENB(a2_type));

  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);


  double R = a1->distance(a2);
  double rho = R/RNB;

  double E =  DENB * ( 1.5 * pow(rho,-10) - 2.5 * pow(rho,-6) );
  return E;


}



double SCREAM_VDW_OBJ::calc_VDW_6_12(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
 /* 6-12 Functional form: 
   * DENB * { rho^(-12) - 2 * rho^(-6) }
   */

  string a1_type = a1->stripped_atomType;
  string a2_type = a2->stripped_atomType;

  //double RNB = this->_geom_mean(this->get_RNB(a1_type), this->get_RNB(a2_type));
  //double DENB = this->_geom_mean(this->get_DENB(a1_type), this->get_DENB(a2_type));
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }

  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);


  double R = a1->distance(a2);
  double rho = R/RNB;

  double rho_m3 = pow(rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m12 = rho_m6 * rho_m6;

  //double E =  DENB * ( pow(rho,-12) - 2 * pow(rho,-6) );
  double E = DENB * (rho_m12 - 2 * rho_m6);
  return E;

}

double SCREAM_VDW_OBJ::calc_VDW_X6(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {

  /* X6 functional  form
   * DENB * { (6/(SCALE-6) ) * exp(SCALE*(1-rho)) - (SCALE/(SCALE-6)*rho^(-6) )}
   * where rho is R/RNB
   */

  // REMARK:  need to take care of close distances: when RHO -> 0, rho^(-6) term -> -inf and dominates.
  
  // default: geom 
  //string a1_type = scream_tools::strip_whitespace(a1->atomType);
  //  string a2_type = scream_tools::strip_whitespace(a2->atomType);

  string a1_type = a1->stripped_atomType;
  string a2_type = a2->stripped_atomType;

  // the following is slow

  //double RNB = this->_geom_mean(this->get_RNB(a1_type), this->get_RNB(a2_type));
  //double DENB = this->_geom_mean(this->get_DENB(a1_type), this->get_DENB(a2_type));
  //double SCALE = (this->get_SCALE(a1_type) + this->get_SCALE(a2_type))/2;
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
    a1->vdw_s = vdw_f->SCALE;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
    a2->vdw_s = vdw_f->SCALE;
  }

  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);
  double SCALE = this->_arith_mean(a1->vdw_s, a2->vdw_s);

  double R = a1->distance(a2);
    //  cout << "R, distance between the two atoms: " << R << endl;

  double rho = R/RNB;
  //  cout << "Standardized rho = R/RNB: " << endl;

  double E = DENB * ( (6/(SCALE-6) ) * exp(SCALE*(1-rho)) - (SCALE/(SCALE-6)* pow(rho,-6) ) )  ;
//   if (E < -10) {
//     cout << "E < -10 ! " << E << endl;
//     a1->dump();
//     a2->dump();
//   }

  // store atom energies

  //a1->energy_per_atom_nb += E / 2;
  //a2->energy_per_atom_nb += E / 2;

  return E;
  

}

double SCREAM_VDW_OBJ::calc_VDW_Morse(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {

}


double SCREAM_VDW_OBJ::calc_Scream_VDW_6_12(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {
  /* Scream functionals: added neighborhood_search_delta.
   * 6-12 Functional form: 
   * DENB * { rho^(-12) - 2 * rho^(-6) }
   */

  string a1_type = a1->stripped_atomType;
  string a2_type = a2->stripped_atomType;

  //double RNB = this->_geom_mean(this->get_RNB(a1_type), this->get_RNB(a2_type));
  //double DENB = this->_geom_mean(this->get_DENB(a1_type), this->get_DENB(a2_type));
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }

  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);


  double R = a1->distance(a2);
  double rho = R/RNB;
  double neighborhood_search_delta = 0.1; // 0.4/4 = 0.1.  0.4: angstrom; 4, typical RNB between two soft atoms (like C).
  // But Max(MaxMin) has a value closer to 0.69A.  I'll use 0.6/4= 0.15 instead.
  neighborhood_search_delta = 0.15;

  double opt_rho = this->_symmetric_optimize_LJ_vdw_functional(rho - neighborhood_search_delta, rho + neighborhood_search_delta );

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m12 = rho_m6 * rho_m6;

  //double E =  DENB * ( pow(opt_rho,-12) - 2 * pow(opt_rho,-6) );
  double E = DENB * ( rho_m12 - 2 * rho_m6 );
  return E;

}

double SCREAM_VDW_OBJ::calc_full_delta_VDW_6_12(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {

  /* calculates delta adjusted energy */

  /* 9-21-06 update: get mu value directly from scream atom; those fields should be instantiated before being used here or anywhere else. */

  /* First, initiliaize if not already initialized */
  this->_checkAndInitFullAtomVdwAndDeltaFields(a1, a2, t);

  /* Then do real calculations */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double a1_delta = a1->delta;
  double a2_delta = a2->delta;
  
  //double total_delta = a1_delta + a2_delta;
  double total_delta = sqrt(a1_delta * a1_delta + a2_delta * a2_delta);

  //if ((total_delta > a1_delta) and (total_delta > a2_delta)) total_delta /= 1.4142;  // divided by sqrt(2): from statistics, assuming 2 underlying iid distributions.
  
  total_delta /=  RNB;
  
  //cout << "total_delta is " << total_delta << endl;
  double opt_rho = this->_symmetric_optimize_LJ_vdw_functional(rho - total_delta, rho + total_delta);
  
  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m12 = rho_m6 * rho_m6;

  //double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  double E = DENB * ( rho_m12 - 2 * rho_m6 );
  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_VDW_6_12(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* calculates delta adjusted energy */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, delta);


  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10) {
    total_delta = delta;
  }
  total_delta /=  RNB;
  
  double opt_rho = this->_symmetric_optimize_LJ_vdw_functional(rho - total_delta, rho + total_delta);

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m12 = rho_m6 * rho_m6;

  //double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  double E = DENB * ( rho_m12 - 2 * rho_m6 );
  return E;

}

double SCREAM_VDW_OBJ::calc_VDW_6_12_scaled_inner_wall(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double scale) {
  /* calculates energy using a scaled inner wall for the 6_12 potential */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, scale);

  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }

  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;
  double rho_n6 = pow(rho, -6); 
  //double E = DENB * ( pow(rho, -12) - 2 * pow(rho, -6));
  double E = DENB * ( rho_n6 * rho_n6 - 2 * rho_n6);
  
  /* now adjust for the scaled potential */

  if (rho < 1) {
    E = (E + DENB)* scale - DENB;
  }

  return E;
}

double SCREAM_VDW_OBJ::calc_full_delta_asym_VDW_6_12(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  /* 9-21: first pre-initialize if not already initialized */
  this->_checkAndInitFullAtomVdwAndDeltaFields(a1, a2, t);

  /* Then do real calclations */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;
  
  double total_delta = (a1->delta > a2->delta) ? a1->delta : a2->delta; 

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  
  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m12 = rho_m6 * rho_m6;

  double E = DENB * ( rho_m12 - 2 * rho_m6 );
  //double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_asym_VDW_6_12(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* First check fields see if populated. */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, delta);

  /* Then do real work. */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10)
    total_delta = delta;

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  //double opt_rho = rho;

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m12 = rho_m6 * rho_m6;

  double E = DENB * ( rho_m12 - 2 * rho_m6 );
  //  double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  return E;

}

/* 6-11 */
double SCREAM_VDW_OBJ::calc_full_delta_asym_VDW_6_11(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  /* 9-21: first pre-initialize if not already initialized */
  this->_checkAndInitFullAtomVdwAndDeltaFields(a1, a2, t);

  /* Then do real calclations */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;
  
  double total_delta = (a1->delta > a2->delta) ? a1->delta : a2->delta; 

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  
  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m11 = rho_m6 * rho_m6 * opt_rho;

  double E = DENB * ( 1.2 * rho_m11 - 2.2 * rho_m6 );
  //double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_asym_VDW_6_11(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* First check fields see if populated. */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, delta);

  /* Then do real work. */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10)
    total_delta = delta;

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  //double opt_rho = rho;

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m11 = rho_m6 * rho_m6 * opt_rho;

  double E = DENB * ( 1.2 * rho_m11 - 2.2 * rho_m6 );

  return E;
}
 
/* 6-10 */
double SCREAM_VDW_OBJ::calc_full_delta_asym_VDW_6_10(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  /* 9-21: first pre-initialize if not already initialized */
  this->_checkAndInitFullAtomVdwAndDeltaFields(a1, a2, t);

  /* Then do real calclations */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;
  
  double total_delta = (a1->delta > a2->delta) ? a1->delta : a2->delta; 

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  
  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m10 = rho_m6 * rho_m6 * opt_rho * opt_rho;

  double E = DENB * ( 1.5 * rho_m10 - 2.5 * rho_m6 );
  //double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_asym_VDW_6_10(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* First check fields see if populated. */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, delta);

  /* Then do real work. */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10)
    total_delta = delta;

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  //double opt_rho = rho;

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m10 = rho_m6 * rho_m6 * opt_rho * opt_rho;

  double E = DENB * ( 1.5 * rho_m10 - 2.5 * rho_m6 );


  //  double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  return E;

}
/* 6-9 */
double SCREAM_VDW_OBJ::calc_full_delta_asym_VDW_6_9(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  /* 9-21: first pre-initialize if not already initialized */
  this->_checkAndInitFullAtomVdwAndDeltaFields(a1, a2, t);

  /* Then do real calclations */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;
  
  double total_delta = (a1->delta > a2->delta) ? a1->delta : a2->delta; 

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  
  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m9 = rho_m6 * rho_m3;

  double E = DENB * ( 2 * rho_m9 - 3 * rho_m6 );

  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_asym_VDW_6_9(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* First check fields see if populated. */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, delta);

  /* Then do real work. */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10)
    total_delta = delta;

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  //double opt_rho = rho;

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m9 = rho_m6 * rho_m3;

  double E = DENB * ( 2 * rho_m9 - 3 * rho_m6 );

  return E;

}
 
/* 6-8 */
double SCREAM_VDW_OBJ::calc_full_delta_asym_VDW_6_8(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  /* 9-21: first pre-initialize if not already initialized */
  this->_checkAndInitFullAtomVdwAndDeltaFields(a1, a2, t);

  /* Then do real calclations */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;
  
  double total_delta = (a1->delta > a2->delta) ? a1->delta : a2->delta; 

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  
  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m8 = rho_m6 * rho_m3 * opt_rho;

  double E = DENB * ( 3 * rho_m8 - 4 * rho_m6 );
  //double E = DENB * ( pow(opt_rho, -12) - 2 * pow(opt_rho, -6));
  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_asym_VDW_6_8(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* First check fields see if populated. */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, delta);

  /* Then do real work. */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10)
    total_delta = delta;

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  //double opt_rho = rho;

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m8 = rho_m6 * rho_m3 * opt_rho;

  double E = DENB * ( 3 * rho_m8 - 4 * rho_m6 );

  return E;

}

/* 6-7 */
double SCREAM_VDW_OBJ::calc_full_delta_asym_VDW_6_7(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  /* 9-21: first pre-initialize if not already initialized */
  this->_checkAndInitFullAtomVdwAndDeltaFields(a1, a2, t);

  /* Then do real calclations */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;
  
  double total_delta = (a1->delta > a2->delta) ? a1->delta : a2->delta; 

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  
  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m7 = rho_m6 / opt_rho;

  double E = DENB * ( 6 * rho_m7 - 7 * rho_m6 );

  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_asym_VDW_6_7(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* First check fields see if populated. */
  this->_checkAndInitFlatAtomVdwAndDeltaFields(a1, a2, delta);

  /* Then do real work. */
  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10)
    total_delta = delta;

  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);

  double rho_m3 = pow(opt_rho,-3);
  double rho_m6 = rho_m3 * rho_m3;
  double rho_m7 = rho_m6 / opt_rho;

  double E = DENB * ( 6 * rho_m7 - 7 * rho_m6 );

  return E;

}
 
/* end 6-7 */



double SCREAM_VDW_OBJ::calc_full_delta_asym_X6(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  
  /* calculates delta adjusted energy */
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
    a1->vdw_s = vdw_f->SCALE;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
    a2->vdw_s = vdw_f->SCALE;
  }
  if (a1->delta < -500) {
    VDW_delta_fields* v_f = (*(this->vdw_delta_library_dict[a1->library_name]))[AtomResInfo(a1->oneLetterResName, a1->stripped_atomLabel)];
    a1->delta = v_f->mu + t * v_f->sigma;
  }
  if (a2->delta < -500) {
    VDW_delta_fields* v_f = (*(this->vdw_delta_library_dict[a2->library_name]))[AtomResInfo(a2->oneLetterResName, a2->stripped_atomLabel)];
    a2->delta = v_f->mu + t * v_f->sigma;
  }

  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);
  double SCALE = this->_geom_mean(a1->vdw_s, a2->vdw_s);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double a1_delta = a1->delta;
  double a2_delta = a2->delta;
  
  double total_delta = a1_delta + a2_delta;
  if ((total_delta > a1_delta) and (total_delta > a2_delta)) total_delta /= 1.4142;  // divided by sqrt(2): from statistics, assuming 2 underlying iid distributions.
  total_delta /=  RNB;
  
  //cout << "total_delta is " << total_delta << endl;
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho);
  double E = DENB * ( (6/(SCALE-6) ) * exp(SCALE*(1-opt_rho)) - (SCALE/(SCALE-6)* pow(opt_rho,-6) ) );
  return E;

}

double SCREAM_VDW_OBJ::calc_flat_delta_asym_X6(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  /* calculates delta adjusted energy */
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
    a1->vdw_s = vdw_f->SCALE;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
    a2->vdw_s = vdw_f->SCALE;
  }

  double RNB = this->_geom_mean(a1->vdw_r, a2->vdw_r);
  double DENB = this->_geom_mean(a1->vdw_d, a2->vdw_d);
  double SCALE = this->_arith_mean(a1->vdw_s, a2->vdw_s);

  double R = a1->distance(a2);
  double rho = R/RNB;

  double total_delta = 0;
  if (R < 10) {
    total_delta = delta;
  }
  total_delta /=  RNB;
  
  double opt_rho = this->_asymmetric_optimize_LJ_vdw_functional(rho - total_delta, rho); // the asymmetric optimize function same for both LJ and X6.
  double E = DENB * ( (6/(SCALE-6) ) * exp(SCALE*(1-opt_rho)) - (SCALE/(SCALE-6)* pow(opt_rho,-6) ) )  ;
  return E;

}


double SCREAM_VDW_OBJ::_asymmetric_optimize_LJ_vdw_functional(double inner_wall, double rho) {
  /* temporarily: just testing "only inner wall delta" method. */
  /* return rho if rho > 1, i.e. lies on outer wall; if rho < 1, than as long as rho belongs in [1-delta, 1], return rho, else if rho < 1-delta, return rho + delta */
  double delta = rho - inner_wall;

  if (rho >= 1) {
    return rho;
  }
  else { // i.e. if rho < 1
    return (rho + delta > 1) ? 1 : rho + delta;
  }
}


void SCREAM_VDW_OBJ::_checkAndInitFullAtomVdwAndDeltaFields(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double t) {
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initialized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }
  if (a1->delta < -500) {
    VDW_delta_fields* v_f = (*(this->vdw_delta_library_dict[a1->library_name]))[AtomResInfo(a1->oneLetterResName, a1->stripped_atomLabel)];
    if (v_f == NULL) {
      cout << "Can't find v_f delta value for : " << a1->library_name << " " << a1->oneLetterResName << " " << a1->stripped_atomLabel << endl;
      a1->delta = 0;
    }
    else 
      a1->delta = v_f->mu + t * v_f->sigma;
  }
  if (a2->delta < -500) {
    VDW_delta_fields* v_f = (*(this->vdw_delta_library_dict[a2->library_name]))[AtomResInfo(a2->oneLetterResName, a2->stripped_atomLabel)];
    if (v_f == NULL) {
      cout << "Can't find v_f delta value for : " << a2->library_name << " " << a2->oneLetterResName << " " << a2->stripped_atomLabel << endl;
      a2->delta = 0;
    }
    else 
      a2->delta = v_f->mu + t * v_f->sigma;
  }

}

void SCREAM_VDW_OBJ::_checkAndInitFlatAtomVdwAndDeltaFields(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double delta) {
  if (a1->vdw_r < -500 ) { // in SCREAM_ATOM constructor, -600 is initislized value.
    VDW_fields* vdw_f = this->get_VDW_fields(a1->stripped_atomType);
    a1->vdw_r = vdw_f->RNB;
    a1->vdw_d = vdw_f->DENB;
  }
  if (a2->vdw_r < -500) {
    VDW_fields* vdw_f = this->get_VDW_fields(a2->stripped_atomType);
    a2->vdw_r = vdw_f->RNB;
    a2->vdw_d = vdw_f->DENB;
  }


}

double SCREAM_VDW_OBJ::_symmetric_optimize_LJ_vdw_functional(double inner_wall, double outer_wall) {
  /* 
   * Simply: if 1 belongs in [inner_wall, outer_wall], return 1.  else, return outer_wall if 1 > outer_wall, return inner_wall if 1 < inner_wall.
   */
  if (inner_wall <= 1 and 1 <= outer_wall) {
    return 1;
  }

  else {
    return (outer_wall < 1) ? outer_wall : inner_wall;
  }
}

SCREAM_VDW_BASE_FUNCTIONAL_OBJ::SCREAM_VDW_BASE_FUNCTIONAL_OBJ(SCREAM_VDW_OBJ* vdw_obj) {

  this->scream_vdw_obj = vdw_obj;

}

SCREAM_VDW_BASE_FUNCTIONAL_OBJ::~SCREAM_VDW_BASE_FUNCTIONAL_OBJ() {

}

SCREAM_calc_full_delta_VDW_6_12::SCREAM_calc_full_delta_VDW_6_12(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_VDW_6_12::~SCREAM_calc_full_delta_VDW_6_12() {
  
}

double SCREAM_calc_full_delta_VDW_6_12::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_VDW_6_12(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_VDW_6_12::SCREAM_calc_flat_delta_VDW_6_12(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_VDW_6_12::~SCREAM_calc_flat_delta_VDW_6_12() {

}

double SCREAM_calc_flat_delta_VDW_6_12::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_VDW_6_12(a1, a2, delta);
}


// The Asymmetric ones.
/* 6-12 */

SCREAM_calc_full_delta_asym_VDW_6_12::SCREAM_calc_full_delta_asym_VDW_6_12(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_asym_VDW_6_12::~SCREAM_calc_full_delta_asym_VDW_6_12() {
  
}

double SCREAM_calc_full_delta_asym_VDW_6_12::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_asym_VDW_6_12(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_asym_VDW_6_12::SCREAM_calc_flat_delta_asym_VDW_6_12(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_asym_VDW_6_12::~SCREAM_calc_flat_delta_asym_VDW_6_12() {

}

double SCREAM_calc_flat_delta_asym_VDW_6_12::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_asym_VDW_6_12(a1, a2, delta);
}

/* 6-11 */
SCREAM_calc_full_delta_asym_VDW_6_11::SCREAM_calc_full_delta_asym_VDW_6_11(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_asym_VDW_6_11::~SCREAM_calc_full_delta_asym_VDW_6_11() {
  
}

double SCREAM_calc_full_delta_asym_VDW_6_11::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_asym_VDW_6_11(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_asym_VDW_6_11::SCREAM_calc_flat_delta_asym_VDW_6_11(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_asym_VDW_6_11::~SCREAM_calc_flat_delta_asym_VDW_6_11() {

}

double SCREAM_calc_flat_delta_asym_VDW_6_11::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_asym_VDW_6_11(a1, a2, delta);
}



/* 6-10 */
SCREAM_calc_full_delta_asym_VDW_6_10::SCREAM_calc_full_delta_asym_VDW_6_10(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_asym_VDW_6_10::~SCREAM_calc_full_delta_asym_VDW_6_10() {
  
}

double SCREAM_calc_full_delta_asym_VDW_6_10::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_asym_VDW_6_10(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_asym_VDW_6_10::SCREAM_calc_flat_delta_asym_VDW_6_10(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_asym_VDW_6_10::~SCREAM_calc_flat_delta_asym_VDW_6_10() {

}

double SCREAM_calc_flat_delta_asym_VDW_6_10::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_asym_VDW_6_10(a1, a2, delta);
}

/* 6-9 */
SCREAM_calc_full_delta_asym_VDW_6_9::SCREAM_calc_full_delta_asym_VDW_6_9(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_asym_VDW_6_9::~SCREAM_calc_full_delta_asym_VDW_6_9() {
  
}

double SCREAM_calc_full_delta_asym_VDW_6_9::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_asym_VDW_6_9(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_asym_VDW_6_9::SCREAM_calc_flat_delta_asym_VDW_6_9(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_asym_VDW_6_9::~SCREAM_calc_flat_delta_asym_VDW_6_9() {

}

double SCREAM_calc_flat_delta_asym_VDW_6_9::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_asym_VDW_6_9(a1, a2, delta);
}

/* 6-8 */
SCREAM_calc_full_delta_asym_VDW_6_8::SCREAM_calc_full_delta_asym_VDW_6_8(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_asym_VDW_6_8::~SCREAM_calc_full_delta_asym_VDW_6_8() {
  
}

double SCREAM_calc_full_delta_asym_VDW_6_8::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_asym_VDW_6_8(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_asym_VDW_6_8::SCREAM_calc_flat_delta_asym_VDW_6_8(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_asym_VDW_6_8::~SCREAM_calc_flat_delta_asym_VDW_6_8() {

}

double SCREAM_calc_flat_delta_asym_VDW_6_8::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_asym_VDW_6_8(a1, a2, delta);
}

/* 6-7 */
SCREAM_calc_full_delta_asym_VDW_6_7::SCREAM_calc_full_delta_asym_VDW_6_7(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_asym_VDW_6_7::~SCREAM_calc_full_delta_asym_VDW_6_7() {
  
}

double SCREAM_calc_full_delta_asym_VDW_6_7::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_asym_VDW_6_7(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_asym_VDW_6_7::SCREAM_calc_flat_delta_asym_VDW_6_7(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_asym_VDW_6_7::~SCREAM_calc_flat_delta_asym_VDW_6_7() {

}

double SCREAM_calc_flat_delta_asym_VDW_6_7::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_asym_VDW_6_7(a1, a2, delta);
}
/* end 6-7 */




SCREAM_calc_full_delta_asym_X6::SCREAM_calc_full_delta_asym_X6(SCREAM_VDW_OBJ* vdw_obj, double n_sigma) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->n_sigma = n_sigma;
}

SCREAM_calc_full_delta_asym_X6::~SCREAM_calc_full_delta_asym_X6() {
  
}

double SCREAM_calc_full_delta_asym_X6::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_full_delta_asym_X6(a1, a2, n_sigma);
}

SCREAM_calc_flat_delta_asym_X6::SCREAM_calc_flat_delta_asym_X6(SCREAM_VDW_OBJ* vdw_obj, double delta) : SCREAM_VDW_BASE_FUNCTIONAL_OBJ(vdw_obj) {
  this->delta = delta;
}

SCREAM_calc_flat_delta_asym_X6::~SCREAM_calc_flat_delta_asym_X6() {

}

double SCREAM_calc_flat_delta_asym_X6::operator()(SCREAM_ATOM* a1, SCREAM_ATOM* a2) const {
  return this->scream_vdw_obj->calc_flat_delta_asym_X6(a1, a2, delta);
}


