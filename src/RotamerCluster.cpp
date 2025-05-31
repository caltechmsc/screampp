#include "RotamerCluster.hpp"

RotamerCluster::RotamerCluster() :  rCI(NULL) {
  // this->allocatedScreamAtoms = false; 
  this->childRotamerList.clear();
}

RotamerCluster::RotamerCluster(Rotamer* rot1, Rotamer* rot2) : rCI(NULL) {
  this->childRotamerList.clear();
  this->childRotamerList.push_back(rot1);
  this->childRotamerList.push_back(rot2);

}

RotamerCluster::~RotamerCluster() {
  
}

void RotamerCluster::addRotamerCluster(Rotamer* newRot) {
  this->childRotamerList.push_back(newRot);
}

vector<Rotamer*> RotamerCluster::getAllRotamers() {
  
  vector<Rotamer*> rot_list; rot_list.clear();

  /* base case */
  if (this->childRotamerList.size() == 1) {
    rot_list.push_back(this->childRotamerList[0]);
    return rot_list;
  }
  /* traversing */
  else {
    for (vector<Rotamer*>::iterator itr = this->childRotamerList.begin();
	 itr != this->childRotamerList.end(); itr++) {
      vector <Rotamer*> childList = (*itr)->getAllRotamers();
      rot_list.insert(rot_list.end(), childList.begin(), childList.end());
    }
  }
  return rot_list;

}

void RotamerCluster::print_Me() const {
  vector<Rotamer*>::const_iterator itr = this->childRotamerList.begin();
  for (; itr != this->childRotamerList.end(); itr++) {
    (*itr)->print_Me(); ///< Recursively called.
  }
}

ScreamAtomV RotamerCluster::get_sc_atoms() {

  ScreamAtomV parentSCAtoms; parentSCAtoms.clear();
  vector<Rotamer*>::const_iterator itr = this->childRotamerList.begin();
  for (; itr != this->childRotamerList.end(); itr++) {
    ScreamAtomV childSCAtoms = (*itr)->get_sc_atoms();
    parentSCAtoms.insert(parentSCAtoms.end(), childSCAtoms.begin(), childSCAtoms.end());
  }
  return parentSCAtoms;

}

ScreamAtomV RotamerCluster::get_bb_atoms() {
  ScreamAtomV parentBBAtoms; parentBBAtoms.clear();
  vector<Rotamer*>::const_iterator itr = this->childRotamerList.begin();
  for (; itr != this->childRotamerList.end(); itr++) {
    ScreamAtomV childBBAtoms = (*itr)->get_bb_atoms();
    parentBBAtoms.insert(parentBBAtoms.end(), childBBAtoms.begin(), childBBAtoms.end() );
  }
  return parentBBAtoms;
}

