#include <cstdlib>
#include "AtomResInfo.hpp"

std::map<std::string, double> AtomResInfo::atomLabel_value_map;

AtomResInfo::AtomResInfo() : resName(""), atomLabel("") {
  if (atomLabel_value_map.empty()) {
    atomLabel_value_map[""] = 0;
    atomLabel_value_map["N"] = 1;
    atomLabel_value_map["NT"] = 1;
    atomLabel_value_map["HN"] = 2;
    atomLabel_value_map["HN1"] = 2.1;
    atomLabel_value_map["HN2"] = 2.2;
    atomLabel_value_map["HN3"] = 2.3;
    atomLabel_value_map["CA"] = 3;
    atomLabel_value_map["HCA"] = 4;
    atomLabel_value_map["HA1"] = 4.1;
    atomLabel_value_map["HA2"] = 4.2;
    atomLabel_value_map["C"] = 500;
    atomLabel_value_map["O"] = 600;
    
    // SideChain.
    atomLabel_value_map["CB"] = 11;
    atomLabel_value_map["HCB"] = 12;
    atomLabel_value_map["HB1"] = 12.1;
    atomLabel_value_map["HB2"] = 12.2;
    atomLabel_value_map["HB3"] = 12.3;
    
    atomLabel_value_map["CG"] = 21;
    atomLabel_value_map["HCG"] = 21.5;
    atomLabel_value_map["OG1"] = 22;
    atomLabel_value_map["HOG1"] = 22.5;
    atomLabel_value_map["CG1"] = 23;
    atomLabel_value_map["HCG1"] = 24;
    atomLabel_value_map["CG2"] = 25;
    atomLabel_value_map["HCG2"] = 26; // plus all those HG11, etc.  put these in later.
    atomLabel_value_map["SG"] = 27;
    atomLabel_value_map["HSG"] = 28;
    atomLabel_value_map["OG"] = 29;
    atomLabel_value_map["HOG"] = 30;
  
    atomLabel_value_map["CD"] = 31;
    atomLabel_value_map["HCD"] = 32;
    atomLabel_value_map["ND1"] = 33;
    atomLabel_value_map["HND1"] = 34;
    atomLabel_value_map["CD1"] = 35;
    atomLabel_value_map["HCD1"] = 36;
    atomLabel_value_map["CD2"] = 37;
    atomLabel_value_map["HCD2"] = 38;
    atomLabel_value_map["OD1"] = 39;
    atomLabel_value_map["HOD1"] = 39.1;
    atomLabel_value_map["OD2"] = 40;
    atomLabel_value_map["HOD2"] = 40.1;
    atomLabel_value_map["ND2"] = 41;
    atomLabel_value_map["HND2"] = 42;
    atomLabel_value_map["SD"] = 43;
    atomLabel_value_map["HSD"] = 43.1;

    atomLabel_value_map["CE"] = 51;
    atomLabel_value_map["HCE"] = 52;
    atomLabel_value_map["NE"] = 53;
    atomLabel_value_map["HNE"] = 54;
    atomLabel_value_map["NE1"] = 55;
    atomLabel_value_map["HNE1"] = 56;
    atomLabel_value_map["CE1"] = 57;
    atomLabel_value_map["HCE1"] = 58;
    atomLabel_value_map["OE1"] = 59;
    atomLabel_value_map["HOE1"] = 59.1;
    atomLabel_value_map["OE2"] = 60;
    atomLabel_value_map["HOE2"] = 60.1;
    atomLabel_value_map["NE2"] = 61;
    atomLabel_value_map["HNE2"] = 62;
    atomLabel_value_map["CE2"] = 63;
    atomLabel_value_map["HCE2"] = 64;
    atomLabel_value_map["CE3"] = 65;
    atomLabel_value_map["HCE3"] = 66;
    
    atomLabel_value_map["CZ"] = 81;
    atomLabel_value_map["HCZ"] = 82;
    atomLabel_value_map["NZ"] = 83;
    atomLabel_value_map["HNZ"] = 84;
    atomLabel_value_map["CZ2"] = 85;
    atomLabel_value_map["HCZ2"] = 86;
    atomLabel_value_map["CZ3"] = 87;
    atomLabel_value_map["HCZ3"] = 88;
    
    atomLabel_value_map["CH2"] = 91;  // TRP
    atomLabel_value_map["HCH2"] = 92;
    atomLabel_value_map["OH"] = 93; // TYR
    atomLabel_value_map["HOH"] = 94;
    atomLabel_value_map["NH1"] = 95; // ARG
    atomLabel_value_map["HNH1"] = 96;
    atomLabel_value_map["HH11"] = 96.1; // ARN
    atomLabel_value_map["HH12"] = 96.2;
    atomLabel_value_map["NH2"] = 97;
    atomLabel_value_map["HNH2"] = 98;
    atomLabel_value_map["HH21"] = 98.1;
    atomLabel_value_map["HH22"] = 98.2;
    
    
    atomLabel_value_map["OX"] = 997;  
    atomLabel_value_map["HC"] = 998;
    atomLabel_value_map["OXT"] = 999;
    atomLabel_value_map["HOXT"] = 999.1;
    
    atomLabel_value_map["Na"] = 3001;
    atomLabel_value_map["Cl"] = 3002;
  }
}

AtomResInfo::AtomResInfo(std::string resName, std::string atomLabel) {

  this->resName = resName;
  this->atomLabel = atomLabel;

  if (atomLabel_value_map.empty()) {
    atomLabel_value_map[""] = 0;
    atomLabel_value_map["N"] = 1;
    atomLabel_value_map["NT"] = 1;
    atomLabel_value_map["HN"] = 2;
    atomLabel_value_map["HN1"] = 2.1;
    atomLabel_value_map["HN2"] = 2.2;
    atomLabel_value_map["HN3"] = 2.3;
    atomLabel_value_map["CA"] = 3;
    atomLabel_value_map["HCA"] = 4;
    atomLabel_value_map["HA1"] = 4.1;
    atomLabel_value_map["HA2"] = 4.2;
    atomLabel_value_map["C"] = 500;
    atomLabel_value_map["O"] = 600;
    
    // SideChain.
    atomLabel_value_map["CB"] = 11;
    atomLabel_value_map["HCB"] = 12;
    atomLabel_value_map["HB1"] = 12.1;
    atomLabel_value_map["HB2"] = 12.2;
    atomLabel_value_map["HB3"] = 12.3;
    
    atomLabel_value_map["CG"] = 21;
    atomLabel_value_map["HCG"] = 21.5;
    atomLabel_value_map["OG1"] = 22;
    atomLabel_value_map["HOG1"] = 22.5;
    atomLabel_value_map["CG1"] = 23;
    atomLabel_value_map["HCG1"] = 24;
    atomLabel_value_map["CG2"] = 25;
    atomLabel_value_map["HCG2"] = 26; // plus all those HG11, etc.  put these in later.
    atomLabel_value_map["SG"] = 27;
    atomLabel_value_map["HSG"] = 28;
    atomLabel_value_map["OG"] = 29;
    atomLabel_value_map["HOG"] = 30;
  
    atomLabel_value_map["CD"] = 31;
    atomLabel_value_map["HCD"] = 32;
    atomLabel_value_map["ND1"] = 33;
    atomLabel_value_map["HND1"] = 34;
    atomLabel_value_map["CD1"] = 35;
    atomLabel_value_map["HCD1"] = 36;
    atomLabel_value_map["CD2"] = 37;
    atomLabel_value_map["HCD2"] = 38;
    atomLabel_value_map["OD1"] = 39;
    atomLabel_value_map["HOD1"] = 39.1;
    atomLabel_value_map["OD2"] = 40;
    atomLabel_value_map["HOD2"] = 40.1;
    atomLabel_value_map["ND2"] = 41;
    atomLabel_value_map["HND2"] = 42;
    atomLabel_value_map["SD"] = 43;
    atomLabel_value_map["HSD"] = 43.1;
    
    atomLabel_value_map["CE"] = 51;
    atomLabel_value_map["HCE"] = 52;
    atomLabel_value_map["NE"] = 53;
    atomLabel_value_map["HNE"] = 54;
    atomLabel_value_map["NE1"] = 55;
    atomLabel_value_map["HNE1"] = 56;
    atomLabel_value_map["CE1"] = 57;
    atomLabel_value_map["HCE1"] = 58;
    atomLabel_value_map["OE1"] = 59;
    atomLabel_value_map["HOE1"] = 59.1;
    atomLabel_value_map["OE2"] = 60;
    atomLabel_value_map["HOE2"] = 60.1;
    atomLabel_value_map["NE2"] = 61;
    atomLabel_value_map["HNE2"] = 62;
    atomLabel_value_map["CE2"] = 63;
    atomLabel_value_map["HCE2"] = 64;
    atomLabel_value_map["CE3"] = 65;
    atomLabel_value_map["HCE3"] = 66;
    
    atomLabel_value_map["CZ"] = 81;
    atomLabel_value_map["HCZ"] = 82;
    atomLabel_value_map["NZ"] = 83;
    atomLabel_value_map["HNZ"] = 84;
    atomLabel_value_map["CZ2"] = 85;
    atomLabel_value_map["HCZ2"] = 86;
    atomLabel_value_map["CZ3"] = 87;
    atomLabel_value_map["HCZ3"] = 88;
    
    atomLabel_value_map["CH2"] = 91;  // TRP
    atomLabel_value_map["HCH2"] = 92;
    atomLabel_value_map["OH"] = 93; // TYR
    atomLabel_value_map["HOH"] = 94;
    atomLabel_value_map["NH1"] = 95; // ARG
    atomLabel_value_map["HNH1"] = 96;
    atomLabel_value_map["HH11"] = 96.1; // ARN
    atomLabel_value_map["HH12"] = 96.2;
    atomLabel_value_map["NH2"] = 97;
    atomLabel_value_map["HNH2"] = 98;
    atomLabel_value_map["HH21"] = 98.1;
    atomLabel_value_map["HH22"] = 98.2;
    
    
    atomLabel_value_map["OX"] = 997;  
    atomLabel_value_map["HC"] = 998;
    atomLabel_value_map["OXT"] = 999;
    atomLabel_value_map["HOXT"] = 999.1;
    
    atomLabel_value_map["Na"] = 3001;
    atomLabel_value_map["Cl"] = 3002;
  }
}

AtomResInfo::~AtomResInfo() {

}

AtomResInfo& AtomResInfo::operator=(const AtomResInfo& atomResInfo) {
  if (this == &atomResInfo) 
    return *this;
  this->resName = atomResInfo.resName;
  this->atomLabel = atomResInfo.atomLabel;

  return *this;

}

bool AtomResInfo::operator==(const AtomResInfo& atomResInfo) const {
  if (this->resName == atomResInfo.resName and
      this->atomLabel == atomResInfo.atomLabel) {
    return true;
  }
  else return false;

}


bool AtomResInfo::operator<(const AtomResInfo& atomResInfo) const {
  // ALA comes before CYS, 

  if (this->resName < atomResInfo.resName) {
    return true;
  } else if (this->resName == atomResInfo.resName) {
    std::map<std::string, double>::iterator this_itr = this->atomLabel_value_map.find(this->atomLabel);
    std::map<std::string, double>::iterator other_itr = this->atomLabel_value_map.find(atomResInfo.atomLabel);

    if (this_itr == this->atomLabel_value_map.end()) {
      std::cerr << "Atom label value unassigned!  Exiting.  (AtomResInfo operator<() value unassigned for " << this->atomLabel << std::endl;
      exit(2);
    }
    if (other_itr == this->atomLabel_value_map.end()) {
      std::cerr << "Atom label value unassigned!  Exiting.  (AtomResInfo operator<() value unassigned for " << atomResInfo.atomLabel << std::endl;
    }
    
    double self_value = this_itr->second;
    double other_value = other_itr->second;
    //double self_value = this->atomLabel_value_map[this->atomLabel];
    //double other_value = this->atomLabel_value_map[atomResInfo.atomLabel];
    return ( (self_value < other_value) ? true : false );
  } else if (this->resName > atomResInfo.resName) {
    return false;
  }



}


std::ostream& operator<<(std::ostream& os, const AtomResInfo& aI) {

  os << aI.resName << " " << aI.atomLabel << std::endl;
  return os;

}


