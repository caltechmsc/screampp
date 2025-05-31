//vcvicek
#include <cstdlib>

#include "defs.hpp"
#include <math.h>
#include "scream_vector.hpp"
#include "scream_tools.hpp"
#include "scream_E_functionals_hb.hpp"
#include <string>
#include <cassert>

using namespace std;

SCREAM_HB_fields::SCREAM_HB_fields() {

}

SCREAM_HB_fields::SCREAM_HB_fields(double DE_, double RE_) : DE(DE_), RE(RE_) {

}

SCREAM_HB_fields::~SCREAM_HB_fields() {
}

HB_delta_fields::HB_delta_fields() {

}

HB_delta_fields::HB_delta_fields(double mu, double sigma) {

  this->mu = mu;
  this->sigma = sigma;

}

HB_delta_fields::~HB_delta_fields() {

}



////// Start SCREAM_HB_OBJ stuff 


SCREAM_HB_OBJ::SCREAM_HB_OBJ() : R_on(0), R_off(5.0), theta_on(180), theta_off(90) {

}

SCREAM_HB_OBJ::~SCREAM_HB_OBJ() {
  // for each hb_fields in map, delete
  for (map< pair<int, int>, SCREAM_HB_fields*>::iterator itr = this->hb_dict.begin();
       itr != this->hb_dict.end(); ++itr) {
    delete itr->second;
  }

  for (map< string, map<AtomResInfo, HB_delta_fields*> * >::iterator itr = this->hb_delta_library_dict.begin(); 
       itr != this->hb_delta_library_dict.end(); ++itr) {
    for (map<AtomResInfo, HB_delta_fields*>::iterator itr2 = itr->second->begin(); itr2 != itr->second->end(); ++itr2) {
      delete itr2->second;
    }
    delete itr->second;
  }


}


double SCREAM_HB_OBJ::calc_HB_Dre(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D) {
  double angle_AHD = this->_calc_angle(A, H, D);
  //  cout << "angle_AHD: " << angle_AHD << endl;

  if (angle_AHD < theta_off) return 0;

  double R_AD = A->distance(D);
  //  cout << "distance: " << R_AD << endl;

  if (R_AD > this->R_off) return 0;

  double cos_AHD = cos(angle_AHD*3.1415926535897/180);

  //double RE = this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->RE;
  //double DE = fabs(this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->DE );

  
  double RE = this->hb_dict[make_pair(A->hb_da, D->hb_da)]->RE;
  double DE = fabs(this->hb_dict[make_pair(A->hb_da, D->hb_da)]->DE );


  // Formula: E_hb * cos_AHD^4 * Spline(R_AD^2, R_on^2, R_off^2) * Spline(cos_AHD^2, cos_on^2, cos_off^2)
  // E_hb = DE { 5* (RE/R)^12 - 6 * (RE/R)^10 }
  // 
  
  double R_o = RE / R_AD;
  
  double E_hb = DE * ( 5* pow(R_o, 12) - 6* pow(R_o, 10) );

  double E = E_hb * pow(cos_AHD , 4);  // put in spline for R later.
  //double E = E_hb * pow(cos_AHD, 2); 


//   A->dump();
//   D->dump();
//   cout << "HB Energy: " << E << endl;
//   cout << "Distance: " << A->distance(D) << endl;
//   cout << endl;

  return E;

}



double SCREAM_HB_OBJ::calc_HB_CHARMM(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D) {

}

double SCREAM_HB_OBJ::calc_Scream_HB(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D) {
  /* Philosophy:
   * 1. acceptor-hydrogen-donor angle ignored; assumed to take on best possible value.  !!! This hasn't been implemented yet.
   * 2. a "neighborhood_search_delta" is introduced for the acceptor and donor atoms on long sidechains.  This value is the radius around the distance between the Acceptor and Donor atom.which the optimal value for HB is calculated.  I.e.: 
   *    for r in { r_ij - r_delta < r_ij < r_ij + r_delta  } minimize HB_energy(r)   
   *            where r_ij is the distance between the acceptor and donor atoms.
   */
  double angle_AHD = this->_calc_angle(A, H, D);
  if (angle_AHD < theta_off) return 0;

  double R_AD = A->distance(D);
  if (R_AD > this->R_off) return 0;

  double cos_AHD = cos(angle_AHD*3.1415926535897/180);

  map< pair<int, int>, SCREAM_HB_fields*>::iterator mHbDict_Itr = this->hb_dict.find(make_pair(A->hb_da, D->hb_da));
  SCREAM_HB_fields* hb_f = mHbDict_Itr->second;

  double RE = hb_f->RE;
  double DE = hb_f->DE;

  //double RE = this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->RE;
  //double DE = fabs(this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->DE );

  // Formula: Min(E_hb) for R_AD in { R_AD - delta < R_AD < R_AD + delta }
  // where E_hb = DE { 5* (RE/R)^12 - 6 * (RE/R)^10 }
  // 
  
  double rho = R_AD / RE;
  double neighborhood_search_delta = 0.15; // typical value of RE: 2.75.  value wanted for neighborhood_search_delta: 0.4A.  0.4/2.75 = 14.5%.


  double opt_rho = this->_optimize_10_12(rho - neighborhood_search_delta, rho + neighborhood_search_delta );

  double E_hb = DE * ( 5* pow(opt_rho, -12) - 6* pow(opt_rho, -10) );
  double E = E_hb * pow(cos_AHD , 4); // put cosine term back in.
  return E;
}

double SCREAM_HB_OBJ::calc_full_delta_HB(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D, double t) {
  /* Philosophy:
   * 1. acceptor-hydrogen-donor angle ignored; assumed to take on best possible value.  !!! This hasn't been implemented yet.
   * 2. a "neighborhood_search_delta" is introduced for the acceptor and donor atoms on long sidechains.  This value is the radius around the distance between the Acceptor and Donor atom.which the optimal value for HB is calculated.  I.e.: 
   *    for r in { r_ij - r_delta < r_ij < r_ij + r_delta  } minimize HB_energy(r)   
   *            where r_ij is the distance between the acceptor and donor atoms.
   */

  double R_AD = A->distance(D);
  if (R_AD > this->R_off) return 0;

  double angle_AHD = this->_calc_angle(A, H, D);
  if (angle_AHD < theta_off) return 0;


  double cos_AHD = cos(angle_AHD*3.1415926535897/180);

  // Formula: Min(E_hb) for R_AD in { R_AD - delta < R_AD < R_AD + delta }
  // where E_hb = DE { 5* (RE/R)^12 - 6 * (RE/R)^10 }


  // Below: initialization if not already initialized (mutation cases) // 7-30-07: no longer needed, hb_da tests.
  map<string, int>::iterator mStrInt_Itr;
//   if (A->hb_da == -600) { // -600 means it's not populated
//     mStrInt_Itr = this->hb_atom_type_mapping.find(A->stripped_atomType);
//     if (mStrInt_Itr != this->hb_atom_type_mapping.end()) {
//       A->hb_da = mStrInt_Itr->second; // if exists, assign that int.
//     } else {
//       A->hb_da = -1; // otherwise, the atom doesn't play a role in HBonding.
//     }
//   }
//   if (D->hb_da == -600) { // -600 means it's not populated
//     mStrInt_Itr = this->hb_atom_type_mapping.find(D->stripped_atomType);
//     if (mStrInt_Itr != this->hb_atom_type_mapping.end()) {
//       D->hb_da = mStrInt_Itr->second; // if exists, assign that int.
//     } else {
//       D->hb_da = -1; // otherwise, the atom doesn't play a role in HBonding.
//     }
//   }
  
  if ( A->delta < -500) {
    HB_delta_fields *hb_f = (*(this->hb_delta_library_dict[A->library_name]))[AtomResInfo(A->oneLetterResName, A->stripped_atomLabel)];
    if (hb_f == NULL) {
      cout << "Can't find delta value for: " << A->resName << " " << A->atomLabel << " of library: " << A->library_name << endl;
      A->delta = 0;
      
    }
    else 
      A->delta = hb_f->mu + t * hb_f->sigma;
  }

  if ( D->delta < -500) {
    HB_delta_fields* hb_f = (*(this->hb_delta_library_dict[D->library_name]))[AtomResInfo(D->oneLetterResName, D->stripped_atomLabel)];
    if (hb_f == NULL) {
      cout << "Can't find delta value for: " << D->resName << " " << D->atomLabel << " of library: " << D->library_name << endl;
      D->delta = 0;
    }
    else
      D->delta = hb_f->mu + t * hb_f->sigma;
  }

  /* Then initialize Acceptor/Donor params if not already initialized */
  //map< pair<string, string>, SCREAM_HB_fields*>::iterator hb_dict_itr = this->hb_dict.find(make_pair(A->stripped_atomType, D->stripped_atomType) );
  map< pair<int, int>, SCREAM_HB_fields*>::iterator hb_dict_itr = this->hb_dict.find(make_pair(A->hb_da, D->hb_da) );
  if (hb_dict_itr == this->hb_dict.end()) {
    cout << " Hydrogen Bond Acceptor/Donor pair not found: " << endl;
    A->dump();
    D->dump();
    cout << "Exiting. " << endl;
    exit(2);

  }

  double RE = hb_dict_itr->second->RE;
  double DE = fabs(hb_dict_itr->second->DE);
  
  double rho = R_AD / RE;
  //double total_delta = (A->delta > D->delta) ? A->delta : D->delta;
  double total_delta = sqrt(A->delta * A->delta + D->delta * D->delta);

  //if ((total_delta > A_delta) and (total_delta > D_delta)) total_delta /= 1.4142;  // divided by sqrt(2): from statistics, assuming 2 underlying iid distributions.
  total_delta /= RE;

  double opt_rho = this->_optimize_10_12(rho - total_delta, rho + total_delta );

  double E_hb = DE * ( 5* pow(opt_rho, -12) - 6* pow(opt_rho, -10) );
  double E = E_hb * pow(cos_AHD , 4); // put cosine term back in.
  return E;
}


double SCREAM_HB_OBJ::calc_flat_delta_HB(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D, double delta) {
  /* Philosophy:
   * 1. acceptor-hydrogen-donor angle ignored; assumed to take on best possible value.  !!! This hasn't been implemented yet.
   * 2. a "neighborhood_search_delta" is introduced for the acceptor and donor atoms on long sidechains.  This value is the radius around the distance between the Acceptor and Donor atom.which the optimal value for HB is calculated.  I.e.: 
   *    for r in { r_ij - r_delta < r_ij < r_ij + r_delta  } minimize HB_energy(r)   
   *            where r_ij is the distance between the acceptor and donor atoms.
   */

  double R_AD = A->distance(D);
  if (R_AD > this->R_off) return 0;

  double angle_AHD = this->_calc_angle(A, H, D);
  if (angle_AHD < theta_off) return 0;

  double cos_AHD = cos(angle_AHD*3.1415926535897/180);

  //map< pair<string, string>, SCREAM_HB_fields*>::iterator hb_dict_itr = this->hb_dict.find(make_pair(A->stripped_atomType, D->stripped_atomType) );
  map<string, int>::iterator mStrInt_Itr;
  if (A->hb_da == -600) { // -600 means it's not populated
    mStrInt_Itr = this->hb_atom_type_mapping.find(A->stripped_atomType);
    if (mStrInt_Itr != this->hb_atom_type_mapping.end()) {
      A->hb_da = mStrInt_Itr->second; // if exists, assign that int.
    } else {
      A->hb_da = -1; // otherwise, the atom doesn't play a role in HBonding.
    }
  }
  if (D->hb_da == -600) { // -600 means it's not populated
    mStrInt_Itr = this->hb_atom_type_mapping.find(D->stripped_atomType);
    if (mStrInt_Itr != this->hb_atom_type_mapping.end()) {
      D->hb_da = mStrInt_Itr->second; // if exists, assign that int.
    } else {
      D->hb_da = -1; // otherwise, the atom doesn't play a role in HBonding.
    }
  }


  map< pair<int, int>, SCREAM_HB_fields*>::iterator hb_dict_itr = this->hb_dict.find(make_pair(A->hb_da, D->hb_da) );
  if (hb_dict_itr == this->hb_dict.end()) {
    cout << " Hydrogen Bond Acceptor/Donor pair not found: " << endl;
    A->dump();
    D->dump();
    cout << "Exiting. " << endl;
    exit(2);

  }

  double RE = hb_dict_itr->second->RE;
  double DE = fabs(hb_dict_itr->second->DE);


  //  double RE = this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->RE;
  //  double DE = fabs(this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->DE );

  // Formula: Min(E_hb) for R_AD in { R_AD - delta < R_AD < R_AD + delta }
  // where E_hb = DE { 5* (RE/R)^12 - 6 * (RE/R)^10 }
  // 
  
  double rho = R_AD / RE;
  double neighborhood_search_delta = delta / RE; // typical value of RE: 2.75.  value wanted for neighborhood_search_delta: 0.4A.  0.4/2.75 = 14.5%.

  double opt_rho = this->_optimize_10_12(rho - neighborhood_search_delta, rho + neighborhood_search_delta );

  

  double E_hb = DE * ( 5* pow(opt_rho, -12) - 6* pow(opt_rho, -10) );
  double E = E_hb * pow(cos_AHD , 4); // put cosine term back in.

  return E;
}


double SCREAM_HB_OBJ::calc_residue_delta_HB(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D) {
  return 0;
}

double SCREAM_HB_OBJ::calc_scaled_inner_wall_HB(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D, double scale) {
  
  double angle_AHD = this->_calc_angle(A, H, D);
  if (angle_AHD < theta_off) return 0;

  double R_AD = A->distance(D);
  if (R_AD > this->R_off) return 0;

  double cos_AHD = cos(angle_AHD*3.1415926535897/180);

  //map< pair<string, string>, SCREAM_HB_fields*>::iterator hb_dict_itr = this->hb_dict.find(make_pair(A->stripped_atomType, D->stripped_atomType) );
  map< pair<int, int>, SCREAM_HB_fields*>::iterator hb_dict_itr = this->hb_dict.find(make_pair(A->hb_da, D->hb_da) );
  if (hb_dict_itr == this->hb_dict.end()) {
    cout << " Hydrogen Bond Acceptor/Donor pair not found: " << endl;
    A->dump();
    D->dump();
    cout << "Exiting. " << endl;
    exit(2);

  }

  double RE = hb_dict_itr->second->RE;
  double DE = hb_dict_itr->second->DE;

  //  double RE = this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->RE;
  //  double DE = fabs(this->hb_dict[make_pair(A->stripped_atomType, D->stripped_atomType)]->DE );

  // Formula: Min(E_hb) for R_AD in { R_AD - delta < R_AD < R_AD + delta }
  // where E_hb = DE { 5* (RE/R)^12 - 6 * (RE/R)^10 }
  // 
  
  double rho = R_AD / RE;
  double E_hb = DE * ( 5* pow(rho, -12) - 6* pow(rho, -10) );
  double E = E_hb * pow(cos_AHD , 4); // put cosine term back in.

  // now scale.

  if (rho < 1) {
    E = (E + DE)* scale - DE;
  }

  return E;


}

double SCREAM_HB_OBJ::_calc_angle(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D) {
  // gives value between 0 and 180
  ScreamVector A_v(A);
  ScreamVector H_v(H);
  ScreamVector D_v(D);

  double angle = (A_v - H_v).angleBtwn(D_v - H_v);

  return fabs(angle);

}

void SCREAM_HB_OBJ::read_param_line(string line) {
  if (this->hb_atom_type_mapping.empty()) {
    hb_atom_type_mapping["H___A"] = 0;
  }

  vector<string> fields;
  split(line, string(" "), fields);
  string A_ff_type = fields[0];
  vector<string> f; f.clear();

  if (fields[0] == "MPSIM_HB") return;

  split(fields[1], "-", f);
  string D_ff_type = f[1];

  int A_ff_int, D_ff_int;
  /* Initiate map<string, int> hb_atom_type_mapping. */
  map<string, int>::iterator mStrInt_itr;
  mStrInt_itr = this->hb_atom_type_mapping.find(A_ff_type);
  if ( mStrInt_itr == this->hb_atom_type_mapping.end() ) {
    int crntSize = this->hb_atom_type_mapping.size();
    this->hb_atom_type_mapping[A_ff_type] = crntSize; // 0 is H___A.
    A_ff_int = crntSize;
  } else {
    A_ff_int = mStrInt_itr->second;
  }
  mStrInt_itr = this->hb_atom_type_mapping.find(D_ff_type);
  if ( mStrInt_itr == this->hb_atom_type_mapping.end() ) {
    int crntSize = this->hb_atom_type_mapping.size();
    this->hb_atom_type_mapping[D_ff_type] = crntSize;
    D_ff_int = crntSize;
  } else {
    D_ff_int = mStrInt_itr->second;
  }
  
  double DE = atof(fields[3].c_str());
  double RE = atof(fields[4].c_str());

  pair<int, int> HB_AD(A_ff_int, D_ff_int);
      
      
  this->hb_dict[HB_AD] = new SCREAM_HB_fields(DE, RE); // fields[0]: ff_type
  //cout << line << endl;
}

void SCREAM_HB_OBJ::read_param_file(string ff_file){
  Debug debugInfo("SCREAM_HB_OBJ::read_param_file(string ff_file)");
  // NOTE: works only on dreidii322-mpsim.par
  this->hb_atom_type_mapping.clear();
  this->hb_dict.clear();
  this->hb_atom_type_mapping["H___A"] = 0; // default value for H___A.

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

    if (line[0] == '*' or line[0] == '#' or line[0] == '!') continue;
    vector<string> fields;
    split(line, string(" "), fields);
    if (fields[0] == "MPSIM_HB") read_mode = "MPSIM_HB";
    if (fields[0] == "USER") break;
    if (read_mode == "MPSIM_HB") {

      if (fields[0] == "MPSIM_HB") { continue; };

      string A_ff_type = fields[0];
      vector<string> f; f.clear();
      split(fields[1], "-", f);
      string D_ff_type = f[1];

      int A_ff_int, D_ff_int;
      /* Initiate map<string, int> hb_atom_type_mapping. */
      map<string, int>::iterator mStrInt_itr;
      mStrInt_itr = this->hb_atom_type_mapping.find(A_ff_type);
      if ( mStrInt_itr == this->hb_atom_type_mapping.end() ) {
	int crntSize = this->hb_atom_type_mapping.size();
	this->hb_atom_type_mapping[A_ff_type] = crntSize; // 0 is H___A.
	A_ff_int = crntSize;
      } else {
	A_ff_int = mStrInt_itr->second;
      }
      mStrInt_itr = this->hb_atom_type_mapping.find(D_ff_type);
      if ( mStrInt_itr == this->hb_atom_type_mapping.end() ) {
	int crntSize = this->hb_atom_type_mapping.size();
	this->hb_atom_type_mapping[D_ff_type] = crntSize;
	D_ff_int = crntSize;
      } else {
	D_ff_int = mStrInt_itr->second;
      }
      
      /* Then, initiate hb_dict by those integers */

      //pair<string, string> HB_AD(A_ff_type, D_ff_type);
      pair<int, int> HB_AD(A_ff_int, D_ff_int);
      
      double DE = atof(fields[3].c_str());
      double RE = atof(fields[4].c_str());

      this->hb_dict[HB_AD] = new SCREAM_HB_fields(DE, RE); // fields[0]: ff_type
      
    }

  }
  FF_FILE.close();
}


void SCREAM_HB_OBJ::read_Scream_delta_file(std::string file) {
  Debug debugInfo("SCREAM_HB_OBJ::read_Scream_delta_file(std::string file)");
  debugInfo.out("SCREAM DELTA FILENAME IS: " + file );
  ifstream FILE;
  FILE.open(file.c_str());

  if (!FILE.good()) {
    cerr << "Unable to open Scream Delta file: " << file << endl;
    exit(8);
  }

  char line_ch[256];
  string library_type, resName, atomLabel;
  double atom_mu, atom_sigma;

  this->hb_delta_library_dict.clear();
  while (!FILE.eof()) {

    FILE.getline(line_ch, sizeof(line_ch));
    string line(line_ch);

    if (line == "" or int(line[0]) == 0) continue; 
    if (line.substr(0,1) == "*" or line.substr(0,1) == "#" or line.substr(0,1) == "!" or int(line[0]) == 0) continue;

    stringstream ss(line);

    ss >> library_type >> resName >> atomLabel >> atom_mu >> atom_sigma;


    // init one element of the map
    if (this->hb_delta_library_dict.count(library_type) == 0) {
      map<AtomResInfo, HB_delta_fields*> * new_aI_hb_field_map = new map<AtomResInfo, HB_delta_fields*>;
      new_aI_hb_field_map->clear();
      this->hb_delta_library_dict[library_type] = new_aI_hb_field_map;
    }

    map<AtomResInfo, HB_delta_fields*> * aI_hb_map = this->hb_delta_library_dict[library_type];
    AtomResInfo aRI(resName, atomLabel);

    if (aI_hb_map->find(aRI) != aI_hb_map->end() ) {
      delete (*aI_hb_map)[aRI];
      HB_delta_fields* hDF = new HB_delta_fields(atom_mu, atom_sigma);
      aI_hb_map->insert(make_pair(aRI, hDF));
    }
    else {
      HB_delta_fields* hDF = new HB_delta_fields(atom_mu, atom_sigma);
      aI_hb_map->insert(make_pair(aRI, hDF));

    }
	

  }

    // Init the default value for library_type NA, or empty string.  7-27-07: why do i want to do the below?  if NA, user should fix problem.
  //  assert( this->hb_delta_library_dict.size() != 0);
  //  map<string, map<AtomResInfo, HB_delta_fields*> >::iterator sample = this->hb_delta_library_dict.begin();
  
  //  map<AtomResInfo, HB_delta_fields*> new_aI_hb_field_map;
  //  this->hb_delta_library_dict[""] = new_aI_hb_field_map;

  //  for (map<AtomResInfo, HB_delta_fields*>::iterator itr = sample->second.begin(); 
  //       itr != sample->second.end(); ++itr ) {
  //    AtomResInfo arI = (*itr).first;
  //    this->hb_delta_library_dict[""][arI] = new HB_delta_fields(0,0); // init to 0
  //}

  FILE.close();

}

SCREAM_HB_fields* SCREAM_HB_OBJ::get_HB_fields(SCREAM_ATOM* A, SCREAM_ATOM* D) {

  
  string A_ff_type = A->stripped_atomType;
  string D_ff_type = D->stripped_atomType;

  int A_ff_int = this->hb_atom_type_mapping.find(A_ff_type)->second;
  int D_ff_int = this->hb_atom_type_mapping.find(D_ff_type)->second;

  // checking: comment this out for speed.  or if have separate DEBUG mode or something
  //map< pair<string, string>, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(make_pair(A_ff_type, D_ff_type));
  map< pair<int, int>, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(make_pair(A_ff_int, D_ff_int));
  if (itr == this->hb_dict.end()) {
    cerr << "Can't find this HB Donor Pair type!  " << A_ff_type << " " << D_ff_type << endl;
    exit(2);
  }
  
  // if no checking just return SCREAM_HB_fields*
  return (*itr).second;

}


double SCREAM_HB_OBJ::get_RE(SCREAM_ATOM* A, SCREAM_ATOM* D) {

  string A_ff_type = A->stripped_atomType;
  string D_ff_type = D->stripped_atomType;

  int A_ff_int = this->hb_atom_type_mapping.find(A_ff_type)->second;
  int D_ff_int = this->hb_atom_type_mapping.find(D_ff_type)->second;

  // checking: comment this out for speed.  or if have separate DEBUG mode or something
  //  map<HB_Acceptor_Donor_Pair_FF_type, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(HB_Acceptor_Donor_Pair_FF_type(A_ff_type, D_ff_type));
  //  map< pair<string, string>, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(make_pair(A_ff_type, D_ff_type));
  map< pair<int, int>, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(make_pair(A_ff_int, D_ff_int));

  if (itr == this->hb_dict.end()) {
    cerr << "Can't find this HB Donor Pair type!  " << A_ff_type << " " << D_ff_type << endl;
    exit(2);
  }
  
  // if no checking just return SCREAM_HB_fields->RE
  return (*itr).second->RE;
  

}

double SCREAM_HB_OBJ::get_DE(SCREAM_ATOM* A, SCREAM_ATOM* D) {

  string A_ff_type = A->stripped_atomType;
  string D_ff_type = D->stripped_atomType;

  int A_ff_int = this->hb_atom_type_mapping.find(A_ff_type)->second;
  int D_ff_int = this->hb_atom_type_mapping.find(D_ff_type)->second;

  // checking: comment this out for speed.  or if have separate DEBUG mode or something
  //  map<HB_Acceptor_Donor_Pair_FF_type, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(HB_Acceptor_Donor_Pair_FF_type(A_ff_type, D_ff_type));
  //  map< pair<string, string>, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(make_pair(A_ff_type, D_ff_type));
  map< pair<int, int>, SCREAM_HB_fields*>::iterator itr = this->hb_dict.find(make_pair(A_ff_int, D_ff_int));

  if (itr == this->hb_dict.end()) {
    cerr << "Can't find this HB Donor Pair type!  " << A_ff_type << " " << D_ff_type << endl;
    exit(2);
  }
  
  // if no checking just return SCREAM_HB_fields->RE
  return (*itr).second->DE;
}

vector<string> SCREAM_HB_OBJ::returnAllHBTypes() const {
  vector<string> allLabels;
  for (map<string, int>::const_iterator itr = this->hb_atom_type_mapping.begin();
       itr != this->hb_atom_type_mapping.end(); ++itr) 
    allLabels.push_back( itr->first );

  return allLabels;


}

double SCREAM_HB_OBJ::_optimize_10_12(double inner_wall, double outer_wall) const {
  /* 
   * Simply: if 1 belongs in [inner_wall, outer_wall], return 1.  else, return outer_wall if 1 > outer_wall, return inner_wall if 1 < inner_wall.
   */
  double opt;
  if (inner_wall <= 1 and 1 <= outer_wall) {
    opt = 1;
    return 1;
  }

  else {
    if (outer_wall < 1) {
      return outer_wall;
    }
    else {
      return inner_wall;
    }
  }

}


SCREAM_HB_BASE_FUNCTIONAL_OBJ::SCREAM_HB_BASE_FUNCTIONAL_OBJ(SCREAM_HB_OBJ* hb_obj) {

  this->scream_hb_obj = hb_obj;

}

SCREAM_HB_BASE_FUNCTIONAL_OBJ::~SCREAM_HB_BASE_FUNCTIONAL_OBJ() {

}

SCREAM_calc_full_delta_HB::SCREAM_calc_full_delta_HB(SCREAM_HB_OBJ* hb_obj, double n) : SCREAM_HB_BASE_FUNCTIONAL_OBJ(hb_obj) {
  this->n_sigma = n;
}

double SCREAM_calc_full_delta_HB::operator()(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D) {
  
  return this->scream_hb_obj->calc_full_delta_HB(A, H, D, this->n_sigma);

}

SCREAM_calc_flat_delta_HB::SCREAM_calc_flat_delta_HB(SCREAM_HB_OBJ* hb_obj, double delta) : SCREAM_HB_BASE_FUNCTIONAL_OBJ(hb_obj) {
  this->delta = delta;
}

double SCREAM_calc_flat_delta_HB::operator()(SCREAM_ATOM* A, SCREAM_ATOM* H, SCREAM_ATOM* D) {

  double E = this->scream_hb_obj->calc_flat_delta_HB(A, H, D, this->delta);
  return E;
}
