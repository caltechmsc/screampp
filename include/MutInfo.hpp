#ifndef MUTINFO_HPP
#define MUTINFO_HPP

#include <string>
#include <iostream>
#include "RotConnInfo.hpp"

using namespace std;

/* SCREAM mutation info */

class MutInfo {
public:
  MutInfo();
  MutInfo(string);
  //MutInfo(int); ///< Constructor from a long to a mutInfo string.  e.g., C218_X <-> 3021824, T219_F <-> 20021906.  I.e.: int in: 2-4-2 format.  INT_MAX	2147483647, in limits.h.
  MutInfo(const MutInfo&);
  MutInfo(const MutInfo&, const MutInfo&);
  virtual ~MutInfo();

  void init(string);
  //  void initInt(int);
  void addMutInfo(const string);
  //  void addMutInfoInt(const int);
  void addMutInfo(const MutInfo&);


  virtual MutInfo& operator=(const MutInfo&);
  virtual bool operator==(const MutInfo&) const;
  virtual bool operator<(const MutInfo&) const;
  friend ostream& operator<<(ostream&, const MutInfo&);

  string getChn() const { return chn;};
  int getPstn() const { return pstn; };
  string getAA() const { return AA; };

  string chn;
  int pstn;
  string AA;

  string str;
  int mIInt;


  virtual void print_Me() const;
  virtual string getString() const; ///< Returns str.
  virtual vector<MutInfo*> getAllMutInfos();

  void setIndex(int i) { this->index = i; };
  int getIndex() const { return this->index; };

  void setRotConnInfo(RotConnInfo* rCI) { this->rCI = rCI;};
  RotConnInfo* getRotConnInfo() const { return this->rCI;};

  map<std::string, RotConnInfo*> getMutInfoStringWithRotConnInfo() const; ///< Returns all MutInfo's with RotConnInfos.  Used for copy constructor and operator=().
  void searchAndAddRotConnInfo(MutInfo, RotConnInfo*); ///< This routine searches the relevant MutInfo (most likely an ArbLib MutInfo) and adds relevant RotConnInfo info to it for easier placement.

  int isClusterMutInfo() const {return this->childMutInfo.size(); } ; 

private:
  vector< MutInfo* > childMutInfo;
  int index; ///< Functions as rotamer # most of the time, when applicable.

  RotConnInfo* rCI; ///< For dealing with arbitrary rotamers.

  void _init_base_MutInfo(string); ///< Initializes a C218_A string.
  int _determineMaxLevel(string) const; ///< Returns the max depth from a string of A||||B||C|B|||D etc.
  bool _cmp_base_MutInfo(const MutInfo&) const; ///< helper function for operator<.
  string _getString_for_Cluster() const; ///< get String for Cluster.
  void _copy_rotConnInfo(const MutInfo&); ///< helper function to copy rotConnInfo from MutInfo.

};

class MutInfoPair {
public:
  MutInfoPair();
  MutInfoPair(const MutInfoPair&);
  MutInfoPair(MutInfo, MutInfo);
  MutInfoPair(MutInfo, MutInfo, double);
  ~MutInfoPair();

  void init(MutInfo, MutInfo);
  MutInfoPair& operator=(const MutInfoPair&);
  bool operator==(const MutInfoPair&) const;
  bool operator<(const MutInfoPair&) const;

  string getString() { string str(""); str += this->mutInfo1.getString(); str += " "; str += this->mutInfo2.getString(); return str; } ;

  MutInfo mutInfo1, mutInfo2;

  MutInfo getMutInfo1() const { return this->mutInfo1;};
  MutInfo getMutInfo2() const { return this->mutInfo2;};

  void setClashE(double E) { this->clashE = E;};
  double getClashE() const { return this->clashE;};

private:
  double clashE;

};

#endif /* MUTINFO_HPP */
