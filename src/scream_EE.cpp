#include "defs.hpp"
#include "RotConnInfo.hpp"
#include "scream_EE.hpp"

#include <sstream>
#include "ClashCollection.hpp"

using namespace std;

Scream_EE::Scream_EE() : ptn(NULL), coulomb_EE(NULL), hb_EE(NULL), vdw_EE(NULL), rotamerNeighborList(NULL), _calcNonPolarHydrogen_flag(1) {
  //this->coulomb_obj = new SCREAM_Coulomb_OBJ();
  //this->vdw_obj = new SCREAM_VDW_OBJ();
  //this->hb_obj = new SCREAM_HB_OBJ();
}


Scream_EE::~Scream_EE() {
  Debug debugInfo("~Scream_EE()");
  debugInfo.out(" destroying scream_EE ... ");

  if (this->coulomb_EE) delete this->coulomb_EE;
  if (this->vdw_EE) delete this->vdw_EE;
  if (this->hb_EE) delete this->hb_EE;
  //  if (this->vdw_hb_exclusion_EE) delete this->hb_vdw_exclusion;

  if (this->rotamerNeighborList) delete this->rotamerNeighborList;

  // Do not delete Protein* ptn.

  debugInfo.out(" Scream_EE destroyed.");

}




void Scream_EE::init(Protein* ptn, vector<std::string> stringV_, string FF_file, string delta_file) {

  this->ptn = ptn;
  for (vector<string>::const_iterator itr = stringV_.begin(); itr != stringV_.end(); ++itr) {
    MutInfo mutInfo(*itr);
    this->mutInfoV.push_back(mutInfo);
  }
  this->_read_FF_param_file(FF_file); // populates relevent fields in coulomb_obj, hb_obj, vdw_obj
  this->_read_SCREAM_delta_param_file(delta_file);

  this->coulomb_EE = new Coulomb_EE(ptn, this->mutInfoV, &(this->coulomb_obj));
  this->vdw_EE = new VDW_EE(ptn, this->mutInfoV, &(this->vdw_obj));
  this->hb_EE = new HB_EE(ptn, this->mutInfoV, &(this->hb_obj));

  // need to initialize vdw_r, vdw_d, vdw_s fields in SCREAM_ATOM.
  //this->_initFixedMoveableAtomsOnProtein(this->ptn, 
  this->_initScreamAtomVdwHbFields(this->_calcNonPolarHydrogen_flag);

}

void Scream_EE::addMutInfoRotConnInfo(MutInfo mutInfo, RotConnInfo* rotConnInfo) {
  this->mutInfo_rotConnInfo_map[mutInfo] = rotConnInfo;
}

void Scream_EE::init_after_addedMutInfoRotConnInfo(Protein* ptn, ScreamParameters* scream_param) {
  Debug debugInfo("Scream_EE::init_after_addedMutInfoRotConnInfo(Protein* ptn, string FF_file, string SCREAM_delta_file)");
  string FF_file = scream_param->getOneEnergyFFParFile();
  string SCREAM_delta_file = scream_param->getDeltaParFile();
  
  string CB_calc_state = scream_param->getCBGroundSpectrumCalc(); // YES/NO.
  string NonPolarHCalc = scream_param->getNonPolarHCalc();

  if (NonPolarHCalc == "YES")
    this->_calcNonPolarHydrogen_flag = 1;
  else
    this->_calcNonPolarHydrogen_flag = 0;
  
  this->ptn = ptn;
  this->_read_FF_param_file(FF_file); // populates relevent fields in coulomb_obj, hb_obj, vdw_obj

  debugInfo.out("Reading SCREAM delta param file... ");
  this->_read_SCREAM_delta_param_file(SCREAM_delta_file); // populates relevent fields in hb_obj, vdw_obj

  this->coulomb_EE = new Coulomb_EE();
  this->vdw_EE = new VDW_EE();
  this->hb_EE = new HB_EE();
  this->vdw_hb_exclusion_EE = new VDW_HB_Exclusion_EE(vdw_EE, hb_EE);

  for (map<MutInfo, RotConnInfo*>::const_iterator itr = this->mutInfo_rotConnInfo_map.begin();
       itr != this->mutInfo_rotConnInfo_map.end(); ++itr ) {
    MutInfo mutInfo = itr->first;
    RotConnInfo* rotConnInfo = itr->second;
    this->coulomb_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
    this->vdw_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
    this->hb_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
  }
  
  this->coulomb_EE->init_after_addedMutInfoRotConnInfo(ptn, &(this->coulomb_obj) );
  debugInfo.out("done coulomb_EE init! ");
  this->vdw_EE->init_after_addedMutInfoRotConnInfo(ptn, &(this->vdw_obj) );
  debugInfo.out("done vdw_EE init! ");
  this->hb_EE->init_after_addedMutInfoRotConnInfo(ptn, &(this->hb_obj) );
  debugInfo.out("done hb_EE init! ");
  
  // need to initialize vdw_r, vdw_d, vdw_s fields in SCREAM_ATOM.
  this->_initScreamAtomVdwHbFields(this->_calcNonPolarHydrogen_flag);


}

void Scream_EE::init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein* ptn, ScreamParameters* scream_param) {
  Debug debugInfo("Scream_EE::init_after_addedMutInfoRotConnInfo_on_the_fly_E(Protein* ptn, std::string  FF_file, std::string SCREAM_delta_file");

  string FF_file = scream_param->getOneEnergyFFParFile();
  string SCREAM_delta_file = scream_param->getDeltaParFile();
  
  string CB_calc_state = scream_param->getCBGroundSpectrumCalc(); // YES/NO.
  string NonPolarHCalc = scream_param->getNonPolarHCalc();

  if (NonPolarHCalc == "YES")
    this->_calcNonPolarHydrogen_flag = 1;
  else
    this->_calcNonPolarHydrogen_flag = 0;


  if (CB_calc_state == "YES")
    this->_CBCalc_flag = 1;
  else
    this->_CBCalc_flag = 0;

  this->ptn = ptn;

  this->_read_FF_param_file(FF_file);
  debugInfo.out("Setting up SCREAM for on-the-fly energy calculations.  No atom list setups. " );
  this->_read_SCREAM_delta_param_file(SCREAM_delta_file); // populates relevent fields in hb_obj, vdw_obj

  this->coulomb_EE = new Coulomb_EE();
  this->vdw_EE = new VDW_EE();
  this->hb_EE = new HB_EE();
  //  this->vdw_hb_exclusion_EE = new VDW_HB_Exclusion_EE(vdw_EE, hb_EE);

  for (map<MutInfo, RotConnInfo*>::const_iterator itr = this->mutInfo_rotConnInfo_map.begin();
       itr != this->mutInfo_rotConnInfo_map.end(); ++itr ) {
    MutInfo mutInfo = itr->first;
    RotConnInfo* rotConnInfo = itr->second;
    this->coulomb_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
    this->vdw_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
    this->hb_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
  }
  
  this->coulomb_EE->init_after_addedMutInfoRotConnInfo_on_the_fly_E(ptn, &(this->coulomb_obj) );
  debugInfo.out("done coulomb_EE! ");
  this->vdw_EE->init_after_addedMutInfoRotConnInfo_on_the_fly_E(ptn, &(this->vdw_obj) );
  debugInfo.out("done vdw_EE! ");
  this->hb_EE->init_after_addedMutInfoRotConnInfo_on_the_fly_E(ptn, &(this->hb_obj) );
  debugInfo.out("done hb_EE! ");
  
  if (scream_param->getCoulombMode() == "Normal")
    this->setNormalDielectric(scream_param->getDielectric());
  else if (scream_param->getCoulombMode() == "DistanceDependent")
    this->setDistanceDielectricPrefactor(scream_param->getDielectric());

  this->_initFixedMoveableAtomsOnProtein(ptn, this->mutInfo_rotConnInfo_map );
  debugInfo.out("done _initFixedMoveableAtomsOnProtein !");
  // Remark: need to add routine for hb_vdw_exclusion on the fly E.

  // need to initialize vdw_r, vdw_d, vdw_s fields in SCREAM_ATOM.
  this->_initScreamAtomVdwHbFields(this->_calcNonPolarHydrogen_flag);
  debugInfo.out("done _initScreamAtomVdwHbFields!");

}

void Scream_EE::init_after_addedMutInfoRotConnInfo_neighbor_list(Protein* ptn, ScreamParameters* scream_param) {

  Debug debugInfo("Scream_EE::init_after_addedMutInfoRotConnInfo_neighbor_list(Protein* ptn, std::string FF_file, std::string SCREAM_delta_file");

  this->ptn = ptn;

  string FF_file = scream_param->getOneEnergyFFParFile();
  string SCREAM_delta_file = scream_param->getDeltaParFile();
  
  string CB_calc_state = scream_param->getCBGroundSpectrumCalc(); // YES/NO.
  string NonPolarHCalc = scream_param->getNonPolarHCalc();

  if (NonPolarHCalc == "YES")
    this->_calcNonPolarHydrogen_flag = 1;
  else
    this->_calcNonPolarHydrogen_flag = 0;


  this->_read_FF_param_file(FF_file);
  debugInfo.out("Setting up rotamer neighbor lists for fast Empty Lattice Energy calculations." );
  this->_read_SCREAM_delta_param_file(SCREAM_delta_file); // populates relevant fields in hb_obj, vdw_obj.

  this->coulomb_EE = new Coulomb_EE();
  this->vdw_EE = new VDW_EE();
  this->hb_EE = new HB_EE();
  this->vdw_hb_exclusion_EE = new VDW_HB_Exclusion_EE(vdw_EE, hb_EE);

  this->rotamerNeighborList = new RotamerNeighborList(ptn, this->mutInfo_rotConnInfo_map, 10000);
  debugInfo.out("RotamerNeighborList initialized!");

  for (map<MutInfo, RotConnInfo*>::const_iterator itr = this->mutInfo_rotConnInfo_map.begin();
       itr != this->mutInfo_rotConnInfo_map.end(); ++itr ) {
    MutInfo mutInfo = itr->first;
    RotConnInfo* rotConnInfo = itr->second;
    this->coulomb_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
    this->vdw_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
    this->hb_EE->addMutInfoRotConnInfo(mutInfo, rotConnInfo);
  }
  
  this->coulomb_EE->init_after_addedMutInfoRotConnInfo_neighbor_list(ptn, &(this->coulomb_obj), this->rotamerNeighborList );
  this->vdw_EE->init_after_addedMutInfoRotConnInfo_neighbor_list(ptn, &(this->vdw_obj), this->rotamerNeighborList );
  this->hb_EE->init_after_addedMutInfoRotConnInfo_neighbor_list(ptn, &(this->hb_obj), this->rotamerNeighborList );
  
  if (scream_param->getCoulombMode() == "Normal")
    this->setNormalDielectric(scream_param->getDielectric());
  else if (scream_param->getCoulombMode() == "DistanceDependent")
    this->setDistanceDielectricPrefactor(scream_param->getDielectric());


  this->_initFixedMoveableAtomsOnProtein(ptn,this->mutInfo_rotConnInfo_map); // still needed; scream_vdw_EE energy calculations need info for  fixed nd moveable atoms.

  // need to initialize vdw_r, vdw_d, vdw_s fields in SCREAM_ATOM.
  this->_initScreamAtomVdwHbFields(this->_calcNonPolarHydrogen_flag);

}

void Scream_EE::fix_mutInfo(MutInfo& mI, RotConnInfo* rCI, int setup) {

  ScreamAtomV variable_atoms; variable_atoms.clear();
  if (rCI == NULL)
    variable_atoms = this->ptn->get_sc_atoms(mI);
  else 
    variable_atoms = this->ptn->get_variable_atoms(rCI);


  for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr)
    (*itr)->make_atom_fixed();

  
  if (setup)
    this->setup_variableAtomsOnEachSidechain();

}

void Scream_EE::moveable_mutInfo(MutInfo& mI, RotConnInfo* rCI, int setup) {  
  ScreamAtomV variable_atoms; variable_atoms.clear();
  if (rCI == NULL) {
    variable_atoms = this->ptn->get_sc_atoms(mI);
    for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr) {
      if (this->_CBCalc_flag == 1) {
	if ( (*itr)->stripped_atomLabel == "CB") // or (*itr)->stripped_atomLabel == "HCB")
	  continue;
	else
	  (*itr)->make_atom_moveable();
      }
      else
	(*itr)->make_atom_moveable();
    }
  }
  else {
    variable_atoms = this->ptn->get_variable_atoms(rCI);
    for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr)
      (*itr)->make_atom_moveable();
  }

  if (setup) 
    this->setup_variableAtomsOnEachSidechain();

}
 
 
void Scream_EE::fix_all() {
  ScreamAtomV atom_list = this->ptn->getAtomList();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    (*itr)->make_atom_fixed();
  }

}
 
void Scream_EE::visible_mutInfo(MutInfo& mI, RotConnInfo* rCI, int setup) {
  ScreamAtomV variable_atoms; variable_atoms.clear();
  if (rCI == NULL) {
    variable_atoms = this->ptn->get_sc_atoms(mI);
    for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr) 
      if (this->_CBCalc_flag == 1) {
	if ( (*itr)->stripped_atomLabel == "CB") // or (*itr)->stripped_atomLabel == "HCB")
	  continue;
	(*itr)->make_atom_visible();
	
      }
  }
  else {
    variable_atoms = this->ptn->get_variable_atoms(rCI);
    for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr)
      (*itr)->make_atom_visible();
  }

  if (setup) 
    this->setup_variableAtomsOnEachSidechain();

}


void Scream_EE::invisible_mutInfo(MutInfo& mI, RotConnInfo* rCI, int setup) {
  ScreamAtomV variable_atoms; variable_atoms.clear();
  if (rCI == NULL) {
    variable_atoms = this->ptn->get_sc_atoms(mI);
    for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr) 
      if (this->_CBCalc_flag == 1) {
	if ( (*itr)->stripped_atomLabel == "CB" ) //or (*itr)->stripped_atomLabel == "HCB")
	  continue;
	(*itr)->make_atom_invisible();
	
      }
  }
  else {
    variable_atoms = this->ptn->get_variable_atoms(rCI);
    for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr)
      (*itr)->make_atom_invisible();
  }

  if (setup) 
    this->setup_variableAtomsOnEachSidechain();

}

// void Scream_EE::visible_EL_mutInfo(MutInfo& mI, RotConnInfo* rCI, int setup) {
//   ScreamAtomV variable_atoms; variable_atoms.clear();
//   if (rCI == NULL) {
//     variable_atoms = this->ptn->get_visible_in_EL_mutInfo_atoms(mI);
//     for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr) 
//       if (this->_CBCalc_flag == 1) {
// 	if ( (*itr)->stripped_atomLabel == "CB") // or (*itr)->stripped_atomLabel == "HCB")
// 	  continue;
// 	(*itr)->make_atom_visible();
	
//       }
//   }
//   else {
//     variable_atoms = this->ptn->get_variable_atoms(rCI);
//     for (ScreamAtomVItr itr = variable_atoms.begin(); itr != variable_atoms.end(); ++itr)
//       (*itr)->make_atom_visible();
//   }

//   if (setup) 
//     this->setup_variableAtomsOnEachSidechain();

// }


void Scream_EE::invisible_EL_mutInfo(MutInfo& mI, RotConnInfo* rCI, int setup) {
  ScreamAtomV previously_visible_atoms = this->ptn->get_visible_in_EL_mutInfo_atoms(mI, rCI);
  for (ScreamAtomVItr itr = previously_visible_atoms.begin(); itr != previously_visible_atoms.end(); ++itr) 
    (*itr)->make_atom_EL_invisible();
  
  if (setup) 
    this->setup_variableAtomsOnEachSidechain();

}

void Scream_EE::visible_EL_mutInfo(MutInfo& mI, RotConnInfo* rCI, int setup) {
  // do nothing for now
}


void Scream_EE::visible_all() {
  // 
  ScreamAtomV atom_list = this->ptn->getAtomList();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr)
    (*itr)->make_atom_visible();
}

void Scream_EE::invisible_all() {
  ScreamAtomV atom_list = this->ptn->getAtomList();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr)
    (*itr)->make_atom_invisible();
}

void Scream_EE::visible_EL_all() {
  ScreamAtomV atom_list = this->ptn->getAtomList();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) 
    (*itr)->make_atom_EL_visible();
}

void Scream_EE::invisible_EL_all() {
  ScreamAtomV atom_list = this->ptn->getAtomList();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr) 
    (*itr)->make_atom_EL_invisible();
  
}

void Scream_EE::make_chain_invisible(string chn) {
  ScreamAtomV atom_list = this->ptn->get_AAChain(chn)->getAtomList();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr)
    (*itr)->make_atom_invisible();

}

void Scream_EE::make_chain_EL_invisible(string chn) {
  ScreamAtomV atom_list = this->ptn->get_AAChain(chn)->getAtomList();
  for (ScreamAtomVItr itr = atom_list.begin(); itr != atom_list.end(); ++itr)
    (*itr)->make_atom_EL_invisible();

}

void Scream_EE::resetFlags(int setup) {
  // Resets atom flag (fixed moveable) to original state (the least significant bit).
  this->ptn->resetFlags();

  if (setup)
    this->setup_variableAtomsOnEachSidechain();
}


void Scream_EE::setup_variableAtomsOnEachSidechain() {
  this->coulomb_EE->setup_variableAtomsOnEachSidechain();
  this->vdw_EE->setup_variableAtomsOnEachSidechain();
  this->hb_EE->setup_variableAtomsOnEachSidechain();
}

void Scream_EE::addClashCollection(ClashCollection* cc) {

  this->vdw_EE->addClashCollection(cc);
  //  this->coulomb_EE->addClashCollection(cc);
  //  this->scream_hb_EE->addClashCollection(cc);

}

void Scream_EE::cleanClashCollection() {
  this->vdw_EE->cleanClashCollection();
}

double Scream_EE::getDistanceDielectricPrefactor() {

  return this->coulomb_obj.getEpsilon();

}

void Scream_EE::setNormalDielectric(double e) {
  this->coulomb_obj.set_normal_mode();
  this->coulomb_obj.set_dielectric(e);
}

int Scream_EE::ntrlRotamerPlacement(string chn, int pstn, AARotamer* rotamer) {
  int mutationFlag = this->ptn->ntrlRotamerPlacement(chn, pstn, rotamer);
  if (mutationFlag) {
    this->setup_variableAtomsOnEachSidechain();
    // Below: can optimize such that only that one residue mutated is changed.
    this->_initScreamAtomVdwHbFields(this->_calcNonPolarHydrogen_flag);
    // No need to call _initFixedMoveableAtomsOnProtein.
  }

}

void Scream_EE::setDistanceDielectricPrefactor(double e) {
  this->coulomb_obj.set_distance_dependent_mode();
  this->coulomb_obj.set_dielectric(e);
}

double Scream_EE::calc_empty_lattice_E(const MutInfo& mutInfo) {
  double total_E = 0;
  //cout << "coulomb_E" << endl;
  double coulomb_E = this->coulomb_EE->calc_empty_lattice_E(mutInfo);
  //cout << "calc vdw" << endl;
  double vdw_E = this->vdw_EE->calc_empty_lattice_E(mutInfo);
  //cout << "calc hb" << endl;
  double hb_E = this->hb_EE->calc_empty_lattice_E(mutInfo);
  //cout << "after all" << endl;
  //double vdw_hb_exclusion_E = this->vdw_hb_exclusion_EE->calc_empty_lattice_E(mutInfo);


  total_E += coulomb_E;
  total_E += vdw_E;
  total_E += hb_E;
  //total_E += vdw_hb_exclusion_E;

  cout << "Self Energy: " << endl;
  cout << "  VDW:     " << vdw_E;
  cout << "  Coulomb: " << coulomb_E;
  cout << "  HB:      " << hb_E << endl;
  //  cout << "  VDW_HB_Ex: " << vdw_hb_exclusion_E << endl;
  cout << "  Total Energy: " << total_E << endl;
  return total_E;

}

double Scream_EE::calc_empty_lattice_E_full_delta(const MutInfo& mutInfo, double t) {

  double E = this->_calc_empty_lattice_E_delta(mutInfo, "FULL", t);
  return E;
}


double Scream_EE::calc_empty_lattice_E_flat_delta(const MutInfo& mutInfo, double d) {

  double E = this->_calc_empty_lattice_E_delta(mutInfo, "FLAT", d);
  return E;

}

double Scream_EE::calc_empty_lattice_E_scaled_inner_wall(const MutInfo& mutInfo, double s) {

  double E = this->_calc_empty_lattice_E_delta(mutInfo, "SCALED", s);
  return E;

}


double Scream_EE::calc_empty_lattice_coulomb_E_delta(const MutInfo& mutInfo) {

  double E = this->coulomb_EE->calc_empty_lattice_E(mutInfo);
  return E;

}

double Scream_EE::calc_empty_lattice_vdw_E_delta(const MutInfo& mutInfo, std::string mode, double r) {
  //  cout << "DEBUG::: in calc_empty_lattice_vdw_E_delta" << endl;
  double E = this->vdw_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);
  return E;
  
}

double Scream_EE::calc_empty_lattice_hb_E_delta(const MutInfo& mutInfo, std::string mode, double r) {
  double E = this->hb_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);
  return E;

}

double Scream_EE::calc_empty_lattice_vdw_hb_exclusion_E_delta(const MutInfo& mutInfo, std::string mode, double r) {

  double E = this->vdw_hb_exclusion_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);
  return E;

}

double Scream_EE::calc_EL_vdw_rot_selfBB(const MutInfo& mI, std::string method, double sigma) {
  return this->vdw_EE->calc_EL_rot_selfBB(mI, method, sigma);
}
 
double Scream_EE::calc_EL_vdw_rot_otherBB(const MutInfo& mI, std::string method, double sigma) {
  return this->vdw_EE->calc_EL_rot_otherBB(mI, method, sigma);
}
 
double Scream_EE::calc_EL_vdw_rot_fixedSC(const MutInfo& mI, std::string method, double sigma) {
  return this->vdw_EE->calc_EL_rot_fixedSC(mI, method, sigma);
}
 
double Scream_EE::calc_EL_vdw_rot_fixedHET(const MutInfo& mI, std::string method, double sigma) {
  return this->vdw_EE->calc_EL_rot_fixedHET(mI, method, sigma);
}

double Scream_EE::calc_EL_vdw_rot_moveableHET(const MutInfo& mI, std::string method, double sigma) {
  return this->vdw_EE->calc_EL_rot_moveableHET(mI, method, sigma);
}

double Scream_EE::calc_EL_coulomb_rot_selfBB(const MutInfo& mI) {
  return this->coulomb_EE->calc_EL_rot_selfBB(mI);
}
 
double Scream_EE::calc_EL_coulomb_rot_otherBB(const MutInfo& mI) {
  return this->coulomb_EE->calc_EL_rot_otherBB(mI);
}

double Scream_EE::calc_EL_coulomb_rot_fixedSC(const MutInfo& mI) {
  return this->coulomb_EE->calc_EL_rot_fixedSC(mI);
}

double Scream_EE::calc_EL_coulomb_rot_fixedHET(const MutInfo& mI) {
  return this->coulomb_EE->calc_EL_rot_fixedHET(mI);
}
 
double Scream_EE::calc_EL_coulomb_rot_moveableHET(const MutInfo& mI) {
  return this->coulomb_EE->calc_EL_rot_moveableHET(mI);
}
  
double Scream_EE::calc_EL_hb_rot_selfBB(const MutInfo& mI, std::string method, double sigma) {
  return this->hb_EE->calc_EL_rot_selfBB(mI, method, sigma);
}
 
double Scream_EE::calc_EL_hb_rot_otherBB(const MutInfo& mI, std::string method, double sigma) {
  return this->hb_EE->calc_EL_rot_otherBB(mI, method, sigma);
}
 
double Scream_EE::calc_EL_hb_rot_fixedSC(const MutInfo& mI, std::string method, double sigma) {
  return this->hb_EE->calc_EL_rot_fixedSC(mI, method, sigma);
}
 
double Scream_EE::calc_EL_hb_rot_fixedHET(const MutInfo& mI, std::string method, double sigma) {
  return this->hb_EE->calc_EL_rot_fixedHET(mI, method, sigma);
}

double Scream_EE::calc_EL_hb_rot_moveableHET(const MutInfo& mI, std::string method, double sigma) {
  return this->hb_EE->calc_EL_rot_moveableHET(mI, method, sigma);
}



double Scream_EE::calc_all_interaction_E() {
  // Called only when DeltaForInteraction == NO, which is currently the default (12-6-05).  Ooops.  In the process of changing this.  Python script should call all VDW, HB and Coulomb subroutines directly for greater control over displaying the energy values.
  double total_E = 0;
  double coulumb_E =  this->coulomb_EE->calc_all_interaction_E();
  double vdw_E = this->vdw_EE->calc_all_interaction_E();
  double hb_E = this->hb_EE->calc_all_interaction_E();
  //double vdw_hb_exclusion_E = this->vdw_hb_exclusion_EE->calc_all_interaction_E();


  total_E += coulumb_E;
  total_E += vdw_E;
  total_E += hb_E;
  //  total_E += vdw_hb_exclusion_E;

  cout << "Interaction Energy: " << endl;
  cout << "VDW:     " << vdw_E << endl;
  cout << "Coulomb: " << coulumb_E << endl;
  cout << "HB:      " << hb_E << endl;
  //cout << "VDW_HB_Ex: " << vdw_hb_exclusion_E << endl;

  return total_E;

}

double Scream_EE::calc_all_interaction_E_full_delta(double t) {

  double E = this->_calc_all_interaction_E_delta("FULL", t);
  return E;

}

double Scream_EE::calc_all_interaction_E_flat_delta(double d) {

  double E = this->_calc_all_interaction_E_delta("FLAT", d);
  return E;

}

double Scream_EE::calc_all_interaction_coulomb_E_delta() {

  double E = this->coulomb_EE->calc_all_interaction_E_delta();
  return E;

}

double Scream_EE::calc_all_interaction_vdw_E_delta(string mode, double r) {

  double E = this->vdw_EE->calc_all_interaction_E_delta(mode, r);
  return E;
}

double Scream_EE::calc_all_interaction_hb_E_delta(string mode, double r) {

  double E = this->hb_EE->calc_all_interaction_E_delta(mode, r);
  return E;
}

double Scream_EE::calc_all_interaction_vdw_hb_exclusion_E_delta(string mode, double r) {
  double E = this->hb_EE->calc_all_interaction_E_delta(mode, r);
  return E;
}


double Scream_EE::calc_residue_interaction_E(const MutInfo mI) {

  double total_E = 0;

  total_E += this->coulomb_EE->calc_residue_interaction_E(mI);
  total_E += this->vdw_EE->calc_residue_interaction_E(mI);
  total_E += this->hb_EE->calc_residue_interaction_E(mI);

  return total_E;

}


double Scream_EE::calc_residue_interaction_E(const MutInfo m1, const MutInfo m2) {
  double total_E = 0;
  total_E += this->coulomb_EE->calc_residue_interaction_E(m1, m2);
  total_E += this->vdw_EE->calc_residue_interaction_E(m1, m2);
  total_E += this->hb_EE->calc_residue_interaction_E(m1, m2);

  return total_E;

}

double Scream_EE::calc_residue_interaction_vdw_E(const MutInfo m1, const MutInfo m2, string method, double scale) {
  return this->vdw_EE->calc_residue_interaction_E(m1, m2, method, scale);
}

double Scream_EE::calc_residue_interaction_hb_E(const MutInfo m1, const MutInfo m2, string method, double scale) {
  return this->hb_EE->calc_residue_interaction_E(m1, m2, method, scale);
}

double Scream_EE::calc_residue_interaction_coulumb_E(const MutInfo m1, const MutInfo m2) {
  return this->coulomb_EE->calc_residue_interaction_E(m1, m2);
}

void Scream_EE::_read_FF_param_file(string ff_file) {
  
  ifstream FF_FILE;
  FF_FILE.open(ff_file.c_str());
  
  if (!FF_FILE.good()) {
    cerr << "Unable to open forcefield file: " << ff_file << endl;
    exit(8);
  }
  
  string line;
  string read_mode = "";

  bool done_VDW = 0;
  bool done_Coulomb = 0;
  bool done_HB = 0;
  
  while (!FF_FILE.eof()) {
    
    char line_ch[256];
    FF_FILE.getline(line_ch, sizeof(line_ch));
    line = string(line_ch);

    if (line == "") continue;
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
    if (fields[0] == "MPSIM_HB") read_mode = "MPSIM_HB";
    if (fields[0] == "DIELCTRIC") read_mode = "DIELCTRIC";
    if (fields[0] == "AUTOTYPE") read_mode = "NULL"; // ends VDW section
    if (fields[0] == "USER") read_mode = "NULL"; // ends MPSIM_HB section

    if (read_mode == "VDW") {
      this->vdw_obj.read_param_line(line);
      done_VDW = 1;
    }
    if (read_mode == "MPSIM_HB") {
      this->hb_obj.read_param_line(line);
      done_HB = 1;
    }
    if (read_mode == "DIELCTRIC") {
      this->coulomb_obj.read_param_line(line);
      done_Coulomb = 1;
      read_mode = "NULL"; // only one line of dielectric info
    }

  }

  if (!done_HB) cout << "HB parameters not parsed!" << endl;
  if (!done_VDW) cout << "VDW parameters not parsed!" << endl;
  if (!done_Coulomb) cout << "Coulomb parameters not parsed! " << endl;

}

void Scream_EE::_read_SCREAM_delta_param_file(string SCREAM_delta_file) {

  this->vdw_obj.read_Scream_delta_file(SCREAM_delta_file);
  this->hb_obj.read_Scream_delta_file(SCREAM_delta_file); // could cut down on time here

}

void Scream_EE::_initScreamAtomVdwHbFields(int flag) {
  /* Remark: flag normally == this->_calcNonPolarHydrogen_flag */
  Debug debugInfo("Scream_EE::_initScreamAtomVdwHbFields(int flag)");
  
  VDW_fields* vdw_f;
  SCREAM_HB_fields* hb_f;
  ScreamAtomV aL = this->ptn->getAtomList();// getAtomList returns reference.

  /* flag == 0: instead of using original atom FF type, figure out FF type on heavy atoms, i.e. C_3 becomes C_33 or C_31. */
  debugInfo.out(" flag value: " + string(itoa(flag)) );
  if (flag == 0) {
    //debugInfo.out("flag == 0");

    for (ScreamAtomVItr itr = aL.begin(); itr != aL.end(); ++itr ) {
      string atom_type = (*itr)->stripped_atomType;
      debugInfo.out(" Atom type of current atom: " + atom_type);
      if ( atom_type[0] == 'C' and atom_type[1] == '_') {
	SCREAM_ATOM* heavyNonPolarAtom = *itr;
	int connectedH_c = 0;
	for (map<SCREAM_ATOM*, int>::iterator connItr = heavyNonPolarAtom->connectivity_m.begin();
	     connItr != heavyNonPolarAtom->connectivity_m.end(); ++connItr) {
	  SCREAM_ATOM* connAtom = (*connItr).first;
	  if (connAtom->stripped_atomType[0] == 'H') { ++connectedH_c; }
	}

	if (connectedH_c) { // not == 0
	  string base_atom_type = atom_type.substr(0,3); // get C_3 or C_R out of C_31 or C_R2
	  string new_atom_type = base_atom_type + string(itoa(connectedH_c));
	  heavyNonPolarAtom->stripped_atomType = new_atom_type; // heavyNonPolarAtom->stripped_atomType attributes.
	}
	// Done with heavyNonPolarAtom FF typing.
      }
    }
    
  }
  /* Then, set the vdw_r etc fields.
  /* flag == 1: calculate everything, use original atom FF type. */
  //  else if (flag == 1) {
  //debugInfo.out("flag == 1");

  for (ScreamAtomVItr itr = aL.begin(); itr != aL.end(); ++itr ) {
    string atom_type = (*itr)->stripped_atomType;
    //    debugInfo.out(" Atom type of current atom: " + atom_type);

    vdw_f = this->vdw_EE->vdw_obj->get_VDW_fields(atom_type);
    
    (*itr)->vdw_r = vdw_f->RNB;
    //    debugInfo.out( stringify((*itr)->vdw_r) );
    (*itr)->vdw_d = vdw_f->DENB;
    //    debugInfo.out( stringify((*itr)->vdw_d) );
    (*itr)->vdw_s = vdw_f->SCALE;
    //    debugInfo.out( stringify((*itr)->vdw_s) );

    // then deal with HB
    map< string, int>::iterator mStrInt_itr;
    mStrInt_itr = this->hb_obj.hb_atom_type_mapping.find(atom_type);
    if (mStrInt_itr != this->hb_obj.hb_atom_type_mapping.end() ) {
      (*itr)->hb_da = mStrInt_itr->second;
      //      debugInfo.out(" atom involved in HBonding. " + atom_type + " " + string(itoa( (*itr)->hb_da) ) );
    }
    else {
      //      debugInfo.out(" atom not involved in HBonding. " + atom_type );
      (*itr)->hb_da = -1; // H___A --> 0, included in hb_atom_type_mapping.  for those not sure, not involved in HBonding, assign value -1.
    }
  }
  //} 

}

void Scream_EE::initScreamAtomDeltaValue(string library_name, string method, double alpha, string eachAtomDeltaFile) {
  /* Remark: run this after other init functions */
  
  /* alpha does double duty: if method is "FLAT" alpha --> FLAT value for all ie mu value, if method if "FULL", alpha --> sigma value, where the mu value would be populated from existing data structures. */

  Debug debugInfo("Scream_EE::_initScreamAtomDeltaFields()");
  
  ScreamAtomV aL = this->ptn->getAtomList();
  for (ScreamAtomVItr itr = aL.begin(); itr != aL.end(); ++itr ) {
    if ( (*itr)->keyw == "HETATM" ) {
      (*itr)->delta = 0;
      continue;
    }
    if ( ((*itr)->flags & 0x1) == 1 ) { // fixed , including backbone atoms.
      (*itr)->delta = 0;
      continue;
    }
    
    if (method.substr(0,4) == "FLAT") {
      (*itr)->delta = alpha;
      continue;
    }

    // otherwise, proceed FULL method initialization.
    //string res_name = scream_tools::one_letter_AA((*itr)->resName);
    string label = (*itr)->stripped_atomLabel;
    string one_letter_name = (*itr)->oneLetterResName;
    VDW_delta_fields* v_f = this->vdw_EE->vdw_obj->get_VDW_delta_fields(library_name, AtomResInfo(one_letter_name, label));

    if (v_f == NULL) { // i.e. no match found 
      cerr << "No SCREAM delta value found for atom n: " << endl;
      cerr << (*itr)->return_bgf_line() << endl;
      cerr << "Delta value set to 0.0." << endl;
      (*itr)->delta = 0; // if not found, assume delta value of 0.
    }
    // below: initialize value of mu.  alpha: is the amount of sigma to reduce/increase
    else 
      (*itr)->delta = v_f->mu + alpha * v_f->sigma;
  }
  
  /* Now read and populate atom fields using values from eachAtomDeltaFile. */
  /* Remark: This is placed after all the standard values have already been populated--i.e. specifications from this file overwrites values from earlier. */
  std::map<int, double> eachAtomDeltaMap;
  this->_read_EachAtomDeltaFile(eachAtomDeltaFile, eachAtomDeltaMap);
  for (int i=0; i<aL.size(); ++i) {
    map<int, double>::iterator _find = eachAtomDeltaMap.find(aL[i]->n);
    if (_find == eachAtomDeltaMap.end() )
	continue;
    else {
      aL[i]->delta = _find->second;
      //cout << "Delta value updated for to " << aL[i]->delta << " for atom " << aL[i]->n << endl;
    }
    
  }
   
  cout << "Delta value updated for " << eachAtomDeltaMap.size() << " atoms." << endl;

}

void Scream_EE::_initFixedMoveableAtomsOnProtein(Protein* ptn, map<MutInfo, RotConnInfo*>) {
  // Non Polar hydrogen calculation cases not handled here yet.
  ScreamAtomV fixed_atoms, all_variable_atoms, variable_atoms_for_one_MutInfo;
  fixed_atoms.clear(), all_variable_atoms.clear(), variable_atoms_for_one_MutInfo.clear();
  
  map<MutInfo, RotConnInfo*>::const_iterator mIrotC_itr = mutInfo_rotConnInfo_map.begin();
  
  for (; mIrotC_itr != mutInfo_rotConnInfo_map.end(); ++mIrotC_itr) {
    /* initializing variable atoms */
    MutInfo mutInfo = mIrotC_itr->first;
    RotConnInfo* rotConnInfo = mIrotC_itr->second;
    string chn = mutInfo.chn;
    int pstn = mutInfo.pstn;

    if (rotConnInfo == NULL) {
      ScreamAtomV tmp_sc_atoms = ptn->get_sc_atoms(mutInfo);
      ScreamAtomV relevantAtoms;
      if (this->_CBCalc_flag == 1) {
	for (ScreamAtomVItr itr = tmp_sc_atoms.begin(); itr != tmp_sc_atoms.end(); ++itr)
	  if ( (*itr)->stripped_atomLabel == "CB") //  or (*itr)->stripped_atomLabel == "HCB")
	    continue;
	  else
	    relevantAtoms.push_back(*itr);
      }
      else
	relevantAtoms.insert(relevantAtoms.end(), tmp_sc_atoms.begin(), tmp_sc_atoms.end());
	
      all_variable_atoms.insert(all_variable_atoms.end(), relevantAtoms.begin(), relevantAtoms.end());
    } else {
      // need to get variable atoms 
      ScreamAtomV tmp_sc_atoms = ptn->get_variable_atoms(rotConnInfo);
      all_variable_atoms.insert(all_variable_atoms.end(), tmp_sc_atoms.begin(), tmp_sc_atoms.end());
    }
  }
  /* initialize fixed atoms.  all non-variable atoms are fixed atoms. */
  ScreamAtomV ptn_atom_list = ptn->getAtomList();
  fixed_atoms.insert(fixed_atoms.end(), ptn_atom_list.begin(), ptn_atom_list.end());
  for (ScreamAtomVItr itr = all_variable_atoms.begin(); itr != all_variable_atoms.end(); ++itr) {
    ScreamAtomVItr variable_atom_i = find(fixed_atoms.begin(), fixed_atoms.end(), *itr);
    fixed_atoms.erase(variable_atom_i);
  }
  /* set fixed flags in SCREAM_ATOM*.  if == 0, moveable.  if ==1, fixed. */
  //  cout << "all_variable_atoms.size(): " << all_variable_atoms.size() << endl;
  for (ScreamAtomVItr itr = all_variable_atoms.begin(); itr != all_variable_atoms.end(); ++itr) {
    (*itr)->flags = 0;
    (*itr)->initFlag();
  }
  for (ScreamAtomVItr itr = fixed_atoms.begin(); itr != fixed_atoms.end(); ++itr ) {
    (*itr)->flags = 1;
    (*itr)->initFlag();
  }
  
}

double Scream_EE::_calc_empty_lattice_E_delta(const MutInfo& mutInfo, string mode, double r) {

  double total_E = 0;
  //  cout << "before coubomb_E" << endl;
  double coulomb_E = this->coulomb_EE->calc_empty_lattice_E(mutInfo);
  //  cout << "before vdw_E" << endl;
  double vdw_E = this->vdw_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);
  //  cout << "before hb_E" << endl;
  double hb_E = this->hb_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);
  // cout << "before vdw_hb_exclusion_E" << endl;
  //double vdw_hb_exclusion_E = this->vdw_hb_exclusion_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);

  total_E += coulomb_E;
  total_E += vdw_E;
  total_E += hb_E;
  //  total_E += vdw_hb_exclusion_E;

  return total_E;
  
}

// double Scream_EE::_calc_empty_lattice_E_delta_asym(const MutInfo& mutInfo, string mode, double r) {
  
//   double total_E = 0;
//   //  cout << "before coubomb_E" << endl;
//   double coulomb_E = this->coulomb_EE->calc_empty_lattice_E(mutInfo);
//   //  cout << "before vdw_E" << endl;
//   double vdw_E = this->vdw_EE->calc_empty_lattice_E_delta_asym(mutInfo, mode, r);
//   //  cout << "before hb_E" << endl;
//   double hb_E = this->hb_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);
//   // cout << "before vdw_hb_exclusion_E" << endl;
//   //double vdw_hb_exclusion_E = this->vdw_hb_exclusion_EE->calc_empty_lattice_E_delta(mutInfo, mode, r);

//   total_E += coulomb_E;
//   total_E += vdw_E;
//   total_E += hb_E;
//   //  total_E += vdw_hb_exclusion_E;

//   return total_E;
 

// } 


double Scream_EE::_calc_all_interaction_E_delta(string mode, double r) {
  
  double total_E = 0;

  double vdw_E = this->vdw_EE->calc_all_interaction_E_delta(mode, r);
  double hb_E = this->hb_EE->calc_all_interaction_E_delta(mode, r);
  double coulomb_E = this->coulomb_EE->calc_all_interaction_E_delta();
  //double vdw_hb_exclusion = this->vdw_hb_exclusion_EE->calc_all_interaction_E_delta(mode, r);
  total_E += vdw_E;
  total_E += hb_E;
  total_E += coulomb_E;

  return total_E;

}

void Scream_EE::_read_EachAtomDeltaFile(std::string file, map<int, double>& deltaMap) {
  /* Simply populates deltaMap map<int, int> structure, where first int is atom n, double is delta value specified by input file.*/

  deltaMap.clear();

  ifstream FILE;
  FILE.open(file.c_str());

  if (!FILE.good()) {
    cerr << "Warning: Unable to open SCREAM each atom delta file: " << file << endl;
    cerr << "If not specified this is probably not a problem.  Proceeding. " << endl;
    return;
  }

  string line;
  int n;
  double delta;

  while (!FILE.eof()) {
    char line_ch[256];
    FILE.getline(line_ch, sizeof(line_ch));
    line = string(line_ch);
    stringstream ss("");
    if (line[0] == '#' or int(line[0]) == 0)
      continue;
    ss << line;
    ss >> n >> delta;
    
    deltaMap[n] = delta;
  }

}
