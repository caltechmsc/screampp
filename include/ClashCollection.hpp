#ifndef CLASHCOLLECTION_HPP
#define CLASHCOLLECTION_HPP

#include "defs.hpp"
#include "MutInfo.hpp"

#include <algorithm>
#include <map>
#include <vector>

using namespace std;

typedef std::map<std::string, unsigned short> ExcitationEnumeration;

class ClashCollection {
public:
  ClashCollection();
  ClashCollection(double); // double: threshold clashing energy.  if energy > double, clash. else, no clash.
  ~ClashCollection();

  void setThresholdE(double); ///< defines the clash threshold energy.
  void addClashPair(MutInfoPair, double); // adds a clashing pair with MutInfo.
  void addClashPair(const MutInfo&, const MutInfo&, double); // adds a clashing pair of rotamers.  Indices are not cinluded; uses information from this->currentRotamerConfiguration, which is information passed in from RotlibCollection.

  int checkClash(ExcitationEnumeration&); // returns 1 if the enumeration contains a pair of clashing rotamers.
  double getThresholdE(); ///< returns the clash threshold energy.
  int getNumberOfClashes(); ///< returns the number of clashing rotamer pairs.


  void storeCurrentRotamerConfiguration(ExcitationEnumeration&); ///< Stores the current rotamer configuration which energy energy will be evaluated/tested.
  void increment_total_clashing_rotamers_eliminated(); ///< this->total_clashing_rotamers_eliminated++.
  void set_total_clashing_rotamers_eliminated(int); ///< Set the number of total clashing rotamers eliminated.
  int get_total_clashing_rotamers_eliminated(); ///< Returns the number of total clashing rotamers eliminated.
  
  //multimap<double, MutInfoPair> getClashMap() { return this->clashMap; } ; ///< Return clashList.
  vector<MutInfoPair> getClashList() { return this->clashList;};

  vector< MutInfoPair > getDiscreteClashPairList(); ///< Returns a list of Clashing Pair that do not share a common node.

private:
  //multimap<double, MutInfoPair> clashMap; ///< When all the clashes are stored.  Not doing this: for some reason, keep crashing.  Don't know why.
  vector<MutInfoPair> clashList;  ///< Stores all the clashing pairs.


  double threshold_E; // above which rotamer-rotamer interaction energy would be considered clashing.
  double bestConfiguration_E; // energy value for the set of rotamers (i.e. configuration) that has the best energy.  absolute energy is used here. 
  int total_clashing_rotamers_eliminated; ///< Stores the number of total clashing rotamers eliminated.

  ExcitationEnumeration currentRotamerConfiguration; ///< state variable.  Updated by: storeCurrentRotamerConfiguration, which is called ONLY by RotlibCollection.  Energy evaluation functions, such as those called by scream_vdw_EE, will input to this class MutInfo and energies.  Those classes lack rotamer indices, so they need to be provided from elsewhere.  RotlibCollection is really the only place where these indices can be accessed.

};

#endif /* CLASHCOLLECTION_HPP */
