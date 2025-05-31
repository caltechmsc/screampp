//vcvicek
#include <cstdlib>
#include <stdlib.h>

#include "defs.hpp"
#include "scream_E_functionals_coulomb.hpp"

SCREAM_Coulomb_OBJ::SCREAM_Coulomb_OBJ() {
  this->epsilon = 1;
}

SCREAM_Coulomb_OBJ::SCREAM_Coulomb_OBJ(string mode_str, double dielectric_) : epsilon(dielectric_){

  if (mode_str == "DistanceDependent")
    this->mode = 2;
  else
    this->mode = 1;

}

SCREAM_Coulomb_OBJ::~SCREAM_Coulomb_OBJ() {

}

double SCREAM_Coulomb_OBJ::calc_Coulomb(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double dielectric) {
  // Note: 7-17-07: dielectric passed in actually isn't used--even though 
  if (this->mode == 1)
    return this->calc_Coulomb_normal(a1, a2, this->epsilon);
  else 
    return this->calc_Coulomb_distance_dielectric(a1, a2, this->epsilon); 
}

double SCREAM_Coulomb_OBJ::calc_Coulomb_normal(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double dielectric) {

  /* Formula: conversion * q_i * q_j / dielectric R_ij * Spline(R_ij, R_on, R_off)
   */

  double q_i = a1->q[0];
  double q_j = a2->q[0];
  double R_ij = a1->distance(a2);

  double epsilon = this->epsilon;
  double conversion = 332.0637;

  double E = conversion * q_i * q_j / (epsilon * R_ij);

  return E;

}


double SCREAM_Coulomb_OBJ::calc_Coulomb_distance_dielectric(SCREAM_ATOM* a1, SCREAM_ATOM* a2, double e) {
  double q_i = a1->q[0];
  double q_j = a2->q[0];
  double R_ij = a1->distance(a2);

  if (R_ij <= 0.0001) { // if distance too small
    R_ij = 0.0001;
  }
  
  //double epsilon = this->dielectric;

  double epsilon = this->epsilon;
  double conversion = 332.0637;

  double E = conversion * q_i * q_j / (epsilon * R_ij);

  if (R_ij > 1.0)
    E /= R_ij; // distance dielectric

//   if (E < -1 or E > 1) {
//     cout << "Strong Electrostatic interactions between: " << E << endl;
//     a1->dump();
//     a2->dump();

//   }
  return E;

}

void SCREAM_Coulomb_OBJ::read_param_line(string line) {
   vector<string> fields;
   split(line, string(" "), fields);
   this->epsilon = atof(fields[1].c_str()) ;
   
}


void SCREAM_Coulomb_OBJ::read_param_file(string ff_file) {

  ifstream  FF_FILE;
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
    if (fields[0] == "DIELCTRIC") read_mode = "DIELCTRIC";
    if (read_mode == "DIELCTRIC") {
      cout << line << endl;
      this->epsilon = atof(fields[1].c_str()) ;
      this->epsilon = this->epsilon; // set these two equal for now. 12-9-05.  just so don't get an undefined value for epsilon.
      //cout << "this->epsilon: " << this->epsilon << endl;
      break;
    }
  }

}
