#ifndef ROTAMERNEIGHBORLIST_HPP
#define ROTAMERNEIGHBORLIST_HPP

#include "defs.hpp"
#include "MutInfo.hpp"
#include "RotConnInfo.hpp"
#include "sc_Protein.hpp"

#include <vector>
#include <map>

using namespace std;

class RotamerNeighborList {
 public:
  RotamerNeighborList();
  RotamerNeighborList(Protein*, map<MutInfo, RotConnInfo*>, double); ///< Contructor: need to know which residues are being screamed.  ScreamAtomV: list of all atoms to be considered.  double, double, double: ine biograf 8.0, 8.5, 9.0 style.

  /* necessary info for neighborlist setup */

  void setCutoff(double cutoff) { this->cutoff = cutoff; };
  double getCutoff() {return this->cutoff;}; 
  void setProtein(Protein*);
  Protein* getProtein();
  void addMutInfoRotConnInfo(MutInfo, RotConnInfo* = NULL);

  /* actual work */

  void initRotamerNeighborList(); ///< Use in conjunctionwith addMutInfoRotConninfo; after all MutInfo added, call this method.

  ScreamAtomV& returnEmptyLatticeNeighborList(MutInfo);

 private:
  double cutoff; ///< Neighbor list cutoff distnce.
  //ScreamAtomV relevantAtomList; ///< Atoms that are relevent to this structure.
  Protein* ptn; ///< Protein underlying.
  map<MutInfo, RotConnInfo*> mutInfo_rotConnInfo; ///< contains information about MutInfo and corresponding RotConnInfo.
  map<MutInfo, ScreamAtomV> emptyLatticeNeighborLists;   ///< the neighbor list for just bb atoms.
  //  map<MutInfo, ScreamAtomV> allAtomsNeighborList; ///< neighborlist for all atoms

  /* Helper functions */
  ScreamAtomV _prepareEmptyLatticeAtomList(Protein*, map<MutInfo, RotConnInfo*>&); ///< Takes in an atom list, removes sidechain portion of residues specified in all <MutInfo, RotConnInfo*>.
  ScreamAtomV _initOneRotamerNeighborList(double, ScreamAtomV&, MutInfo, RotConnInfo*); ///< Returns list of ScreamAtomV that will is the neighbor list for MutInfo.
  ScreamVector _determineCenter(ScreamAtomV&, MutInfo, RotConnInfo*);

};

#endif /* ROTAMERNEIGHBORLIST_HPP */
