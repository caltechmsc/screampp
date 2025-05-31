#ifndef SCREAM_VDW_HB_EXCLUSION_EE
#define SCREAM_VDW_HB_EXCLUSION_EE

#include "sc_Protein.hpp"
#include "MutInfo.hpp"
#include "RotConnInfo.hpp"

#include "scream_vdw_EE.hpp"
#include "scream_hb_EE.hpp"

#include <vector>

class VDW_HB_Exclusion_EE {
public:
  VDW_HB_Exclusion_EE();
  VDW_HB_Exclusion_EE(VDW_EE*, HB_EE*); // do this after init_after_addedMutInfoRotConnInfo in both vdw ee and hb ee.
  ~VDW_HB_Exclusion_EE();

  double calc_empty_lattice_E(const MutInfo&);  // does the same thing as scream_vdw_EE, scream_hb_EE etc EXCEPT it returns the negative of how much VDW energy is calculated in VDW_EE which shouldn't have been due to the HBonding VDW exclusion.
  double calc_empty_lattice_E_delta(const MutInfo&, string , double); // again, returns the negative of over-contribution of VDW due to HBonding VDW exclusion.

  double calc_residue_interaction_E(const MutInfo&); 
  double calc_residue_interaction_E(const MutInfo&, const MutInfo&);
  double calc_all_interaction_E(); 
  double calc_all_interaction_E_delta(string , double = 0);

private:

  int ON_THE_FLY;

  VDW_EE* vdw_EE; // stores a pointer to VDW_EE.
  HB_EE* hb_EE; // stores a pointer to HB_EE.

  bool find_atompair_in_list(AtomPair&, vector<AtomPair>&);

};

#endif /* SCREAM_VDW_HB_EXCLUSION_EE */
