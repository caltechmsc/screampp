#include "scream_vdw_hb_exclusion_EE.hpp"
#include "scream_vdw_EE.hpp"
#include "scream_hb_EE.hpp"

#include "defs.hpp"
#include "MutInfo.hpp"

#include <vector>

VDW_HB_Exclusion_EE::VDW_HB_Exclusion_EE() {

}

VDW_HB_Exclusion_EE::VDW_HB_Exclusion_EE(VDW_EE* vdw_EE, HB_EE* hb_EE) {
  this->vdw_EE = vdw_EE;
  this->hb_EE = hb_EE;
}

VDW_HB_Exclusion_EE::~VDW_HB_Exclusion_EE() {

}

double VDW_HB_Exclusion_EE::calc_empty_lattice_E(const MutInfo& mI) {
  double total_E = 0;
  vector< vector<SCREAM_ATOM*> > atomTripleList(this->hb_EE->variable_and_fixed[mI]);
  vector<AtomPair> already_subtracted_atomPairs;
  
  double theta_off = this->hb_EE->hb_obj->theta_off;
  double R_off = this->hb_EE->hb_obj->R_off;

  vector< vector<SCREAM_ATOM*> >::const_iterator AHD_itr = atomTripleList.begin();
  for (; AHD_itr != atomTripleList.end(); AHD_itr++) {
    SCREAM_ATOM* A = (*AHD_itr)[0];
    SCREAM_ATOM* H = (*AHD_itr)[1];
    SCREAM_ATOM* D = (*AHD_itr)[2];
    double angle_AHD = this->hb_EE->hb_obj->_calc_angle(A, H, D);
    if (angle_AHD < theta_off) continue;
    
    double R_AD = A->distance(D);
    if (R_AD > R_off) continue;

    AtomPair tmpAD(A, D);
    bool already_subtracted = this->find_atompair_in_list(tmpAD, already_subtracted_atomPairs);

    // If reach here, means need to calc VDW-HB-exclusion.
    if ( ! already_subtracted) {
      already_subtracted_atomPairs.push_back(tmpAD);
      total_E += this->vdw_EE->vdw_obj->calc_VDW_6_12(A,D);

    }
  }

  cout << "size of already subtracted atomPairs in VDW_HB exclusion code: " << already_subtracted_atomPairs.size() << endl;

  return (-total_E);
  
  
}

double VDW_HB_Exclusion_EE::calc_empty_lattice_E_delta(const MutInfo& mI, string mode, double r) {
  double total_E = 0;
  vector< vector<SCREAM_ATOM*> > atomTripleList(this->hb_EE->variable_and_fixed[mI]);
  vector<AtomPair> already_subtracted_atomPairs;
  
  double theta_off = this->hb_EE->hb_obj->theta_off;
  double R_off = this->hb_EE->hb_obj->R_off;

  vector< vector<SCREAM_ATOM*> >::const_iterator AHD_itr = atomTripleList.begin();
  for (; AHD_itr != atomTripleList.end(); AHD_itr++) {
    SCREAM_ATOM* A = (*AHD_itr)[0];
    SCREAM_ATOM* H = (*AHD_itr)[1];
    SCREAM_ATOM* D = (*AHD_itr)[2];
    double angle_AHD = this->hb_EE->hb_obj->_calc_angle(A, H, D);
    if (angle_AHD < theta_off) continue;
    
    double R_AD = A->distance(D);
    if (R_AD > R_off) continue;

    AtomPair tmpAD(A, D);
    bool already_subtracted = this->find_atompair_in_list(tmpAD, already_subtracted_atomPairs);

    // If reach here, means need to calc VDW-HB-exclusion.
    if ( ! already_subtracted) {
      already_subtracted_atomPairs.push_back(tmpAD);
      if (mode == "FULL") {
	total_E += this->vdw_EE->vdw_obj->calc_full_delta_VDW_6_12(A, D, r);
      } else if (mode == "FLAT") {
	total_E += this->vdw_EE->vdw_obj->calc_flat_delta_VDW_6_12(A,D,r);
      } else if (mode == "SCALED") {
	total_E += this->vdw_EE->vdw_obj->calc_VDW_6_12_scaled_inner_wall(A,D,r);
      }

    }
  }

  return (-total_E);

}


double VDW_HB_Exclusion_EE::calc_residue_interaction_E(const MutInfo& mI) {
  // will implement this later on.


}

double VDW_HB_Exclusion_EE::calc_residue_interaction_E(const MutInfo& mI1, const MutInfo& mI2) {
  // will implement this later on.
}

double VDW_HB_Exclusion_EE::calc_all_interaction_E() {
  double total_E = 0;
  vector<AtomPair> already_subtracted_atomPairs;
  
  double theta_off = this->hb_EE->hb_obj->theta_off;
  double R_off = this->hb_EE->hb_obj->R_off;
  
  for (map< MutInfoPair, vector<ScreamAtomV> >::iterator itr = this->hb_EE->variable_and_variable.begin(); 
       itr != this->hb_EE->variable_and_variable.end(); itr++) {
    for (vector<ScreamAtomV>::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); AHD_i++) {
      SCREAM_ATOM* A = (*AHD_i)[0];
      SCREAM_ATOM* H = (*AHD_i)[1];
      SCREAM_ATOM* D = (*AHD_i)[2];
      
      double angle_AHD = this->hb_EE->hb_obj->_calc_angle(A, H, D);
      if (angle_AHD < theta_off) continue;
      
      double R_AD = A->distance(D);
      if (R_AD > R_off) continue;
      
      AtomPair tmpAD(A, D);
      bool already_subtracted = this->find_atompair_in_list(tmpAD, already_subtracted_atomPairs);
    
      if ( ! already_subtracted) {
	already_subtracted_atomPairs.push_back(tmpAD);
// 	if (mode == "FULL") {
// 	  total_E += this->vdw_EE->vdw_obj->calc_full_delta_VDW_6_12(A, D, r);
// 	} else if (mode == "FLAT") {
// 	  total_E += this->vdw_EE->vdw_obj->calc_flat_delta_VDW_6_12(A,D,r);
// 	} else if (mode == "RESIDUE") {
// 	  total_E += this->vdw_EE->vdw_obj->calc_residue_delta_VDW_6_12(A,D);
// 	} else if (mode == "SCALED") {
// 	  total_E += this->vdw_EE->vdw_obj->calc_VDW_6_12_scaled_inner_wall(A,D,r);
// 	}
	total_E += this->vdw_EE->vdw_obj->calc_VDW_6_12(A, D);
      }
    }
    

  }
  return (-total_E);

}

double VDW_HB_Exclusion_EE::calc_all_interaction_E_delta(string mode, double r) {
  double total_E = 0;
  vector<AtomPair> already_subtracted_atomPairs;
  
  double theta_off = this->hb_EE->hb_obj->theta_off;
  double R_off = this->hb_EE->hb_obj->R_off;
  
  for (map< MutInfoPair, vector<ScreamAtomV> >::iterator itr = this->hb_EE->variable_and_variable.begin(); 
       itr != this->hb_EE->variable_and_variable.end(); itr++) {
    for (vector<ScreamAtomV>::iterator AHD_i = (itr->second).begin(); AHD_i != (itr->second).end(); AHD_i++) {
      SCREAM_ATOM* A = (*AHD_i)[0];
      SCREAM_ATOM* H = (*AHD_i)[1];
      SCREAM_ATOM* D = (*AHD_i)[2];
      
      double angle_AHD = this->hb_EE->hb_obj->_calc_angle(A, H, D);
      if (angle_AHD < theta_off) continue;
      
      double R_AD = A->distance(D);
      if (R_AD > R_off) continue;
      
      AtomPair tmpAD(A, D);
      bool already_subtracted = this->find_atompair_in_list(tmpAD, already_subtracted_atomPairs);
    
      if ( ! already_subtracted) {
	already_subtracted_atomPairs.push_back(tmpAD);
	if (mode == "FULL") {
	  total_E += this->vdw_EE->vdw_obj->calc_full_delta_VDW_6_12(A, D, r);
	} else if (mode == "FLAT") {
	  total_E += this->vdw_EE->vdw_obj->calc_flat_delta_VDW_6_12(A,D,r);
	} else if (mode == "SCALED") {
	  total_E += this->vdw_EE->vdw_obj->calc_VDW_6_12_scaled_inner_wall(A,D,r);
	}
      }
    }
    
    
  }
  return (-total_E);
  
}


bool VDW_HB_Exclusion_EE::find_atompair_in_list(AtomPair& aP, vector<AtomPair>& aP_v) {
  // returns true if atomPair is in list, else return false.

  for (vector<AtomPair>::const_iterator aP_i = aP_v.begin(); aP_i != aP_v.end(); aP_i++) {
    SCREAM_ATOM* aP_a1 = aP.a1;
    SCREAM_ATOM* aP_a2 = aP.a2;

    SCREAM_ATOM* aP_i_a1 = (*aP_i).a1;
    SCREAM_ATOM* aP_i_a2 = (*aP_i).a2;
    // pointer address comparison is good enough for this.
    if (((aP_a1 == aP_i_a1) and (aP_a2 == aP_i_a2)) or 
	((aP_a1 == aP_i_a2) and (aP_a2 == aP_i_a1)) ) {
      return true;
    }
  }
  return false;
}
