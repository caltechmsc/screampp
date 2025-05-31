#include <assert.h>
#include "ClashCollection.hpp"
#include "MutInfo.hpp"

#include <map>
#include <vector>

ClashCollection::ClashCollection() {
  this->threshold_E = 300; // in absolute terms.
  this->total_clashing_rotamers_eliminated = 0;
  this->currentRotamerConfiguration.clear();
}

ClashCollection::ClashCollection(double threshold_E) {

  this->threshold_E = threshold_E;
  this->total_clashing_rotamers_eliminated = 0;
  this->currentRotamerConfiguration.clear();

}

ClashCollection::~ClashCollection() {
  this->currentRotamerConfiguration.clear();
}

void ClashCollection::setThresholdE(double E) {

  this->threshold_E = E;

}

double ClashCollection::getThresholdE() {
  return this->threshold_E;
}

int ClashCollection::getNumberOfClashes() {
  return this->clashList.size();
}


void ClashCollection::addClashPair(MutInfoPair mP, double energy) {

  if (energy > this->threshold_E) {
    mP.setClashE(energy);
    this->clashList.push_back(mP);
  }

}

void ClashCollection::addClashPair(const MutInfo& m1, const MutInfo& m2, double E) {

  if (E < this->threshold_E) {
    return;
  }

  string m1_string = m1.getString();
  string m2_string = m2.getString();

  int m1_index = this->currentRotamerConfiguration[m1_string];
  int m2_index = this->currentRotamerConfiguration[m2_string];

  MutInfo indexed_m1(m1);
  indexed_m1.setIndex(m1_index);
  MutInfo indexed_m2(m2);
  indexed_m2.setIndex(m2_index);

  /* Below: checking that this clash does not already exsit would be nice; but should be unnecessary.  */
  /* Why? Proof: if current clash pair in clashMap, current configuration would have been discarded (clashing pair).  Therefore, cannot happen. */

  
  MutInfoPair mP(indexed_m1, indexed_m2);
  
  mP.setClashE(E);
  //cout << " Original Clash List Size: " << this->clashList.size() << endl;
  this->clashList.push_back(mP);
  //  cout << " Now Clash List Size: " << this->clashList.size() << endl;
 

}

int ClashCollection::checkClash(ExcitationEnumeration& ee) {


  //  for (multimap<double, MutInfoPair >::const_iterator d_mP_itr = this->clashMap.begin();
  //       d_mP_itr != clashMap.end(); d_mP_itr++) {

  for (vector<MutInfoPair>::const_iterator itr = this->clashList.begin();
       itr != this->clashList.end(); itr++) {

    MutInfoPair tmp_mI = *itr;

    string m1_string = tmp_mI.getMutInfo1().getString();
    string m2_string = tmp_mI.getMutInfo2().getString();
    
    int m1_index = tmp_mI.getMutInfo1().getIndex();
    int m2_index = tmp_mI.getMutInfo2().getIndex();
      
    assert (m1_index != 0);
    assert (m2_index != 0);

    if (ee[m1_string] == m1_index and ee[m2_string] == m2_index)
    {
      return 1;
    }
  }
  
  return 0;

}

void ClashCollection::storeCurrentRotamerConfiguration(ExcitationEnumeration& ee) {

  this->currentRotamerConfiguration = ee;

}

void ClashCollection::increment_total_clashing_rotamers_eliminated() {

  this->total_clashing_rotamers_eliminated++;

}

void ClashCollection::set_total_clashing_rotamers_eliminated(int i) {

  this->total_clashing_rotamers_eliminated = i;

}

int ClashCollection::get_total_clashing_rotamers_eliminated() {

  return this->total_clashing_rotamers_eliminated;

}

class cmp_mutInfo_E_descend {
public:
  bool operator()(const MutInfoPair& mip1, const MutInfoPair& mip2) {
    return (mip1.getClashE() > mip2.getClashE());
    }
};

vector< MutInfoPair > ClashCollection::getDiscreteClashPairList() {
  Debug debugInfo("ClashCollection::getDiscreteClashPairList()");

  map<string, int> mutInfo_appearance; mutInfo_appearance.clear();
  vector< MutInfoPair > clashing_MutInfos; clashing_MutInfos.clear();
  debugInfo.out("Before sorting this->clashList ");
  sort(this->clashList.begin(), this->clashList.end());
  
  debugInfo.out("After sorting this->clashList " );

  for (vector<MutInfoPair >::const_iterator itr = this->clashList.begin(); itr != this->clashList.end(); itr++) {
    MutInfo mI1 = itr->getMutInfo1();
    MutInfo mI2 = itr->getMutInfo2();

    string mI1_string = mI1.getString();
    string mI2_string = mI2.getString();

    //cout << mI1_string << endl;
    //cout << mI2_string << endl;

    bool flag1 = true;
    bool flag2 = true;

    if (mutInfo_appearance.find(mI1_string) == mutInfo_appearance.end() ) {
      mutInfo_appearance[mI1_string] = 0;
      flag1 = true;
    } else {
      flag1 = false;
    }

    if (mutInfo_appearance.find(mI2_string) == mutInfo_appearance.end() ) {
      mutInfo_appearance[mI2_string] = 0;
      flag2 = true;
    } else {
      flag2 = false;
    }
    
    if (flag1 and flag2) {
      clashing_MutInfos.push_back(MutInfoPair(mI1, mI2, itr->getClashE()) );
    }

  }

  // Then sort clashing_MutInfos by energies stored in MutInfoPair objects.

  sort(clashing_MutInfos.begin(), clashing_MutInfos.end(), cmp_mutInfo_E_descend());
  return clashing_MutInfos;
}
