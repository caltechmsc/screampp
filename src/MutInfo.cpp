#include <cstdlib>
#include "defs.hpp"
#include "MutInfo.hpp"


MutInfo::MutInfo() : index(0), rCI(NULL) {

}

MutInfo::MutInfo(string mutInfo) : index(0), rCI(NULL) {
 
  this->init(mutInfo);

}

// MutInfo::MutInfo(int mI_int) : index(0), rCI(NULL) {
//   this->initInt(mI_int);
// }

MutInfo::MutInfo(const MutInfo& mI) {
  string mI_string = mI.getString();
  this->init(mI_string);
  this->index = mI.getIndex();
  this->_copy_rotConnInfo(mI);
}


MutInfo::MutInfo(const MutInfo& mI1, const MutInfo& mI2) : index(0), rCI(NULL) {

  string mI1_string = mI1.getString();
  string mI2_string = mI2.getString();
  
  int max_level_1 = this->_determineMaxLevel(mI1_string);
  int max_level_2 = this->_determineMaxLevel(mI2_string);

  int max_level = (max_level_1 > max_level_2) ? max_level_1 : max_level_2;
  max_level += 1;
  string sep = string(max_level, '|');


  string new_string = mI1_string + sep + mI2_string;

  this->init(new_string);

}

MutInfo::~MutInfo() {
  if (this->isClusterMutInfo() ) {
    for (vector<MutInfo*>::iterator itr = this->childMutInfo.begin();
	 itr != this->childMutInfo.end(); itr++) {
      delete (*itr);
    }
  }
}


void MutInfo::init(string s) {

  this->childMutInfo.clear();

  /* for a string of: A|B||C|||D|E|F||G,
     the object is to setup a tree that represents this structure.
     Recursive algorithm is the easiest here.
   */

  int max_level = this->_determineMaxLevel(s);
  
  /* Base case. */

  if (max_level == 0) {
    this->_init_base_MutInfo(s);
    return;
  }

  /* Take care of AA, pstn, settings if not base case */
  this->chn = "~"; // why "~": greatest value in ASCII code.
  this->pstn = 9999999;
  this->AA = "~";
  
  /* Make string ||||...|   <- max_level times. */
  string sep = string(max_level, '|');
  vector<string> split_s; split_s.clear();

  int last_mark = 0;
  for (int i = 0; i != s.size(); i++) {
    if (s.substr(i,max_level) == sep) {
      split_s.push_back(s.substr(last_mark,i-last_mark));
      last_mark = i+max_level;
    }
    if (i == s.size()-1) {
      split_s.push_back(s.substr(last_mark, i-last_mark+1));
    }
  }

  
  /* Recurse. */
  for (vector<string>::iterator itr = split_s.begin();
       itr != split_s.end(); itr++) {
    MutInfo* newMI = new MutInfo(*itr);
    this->childMutInfo.push_back(newMI);
  }

  return;

}

// void MutInfo::initInt(int mIInt) {

// }

void MutInfo::addMutInfo(const string mI_s) {

  /* Take care of AA, pstn, settings if not base case */
  this->chn = "~"; // why "~": greatest value in ASCII code.
  this->pstn = 9999999;
  this->AA = "~";

  MutInfo* newMI = new MutInfo(mI_s);
  this->childMutInfo.push_back(newMI);

}

void MutInfo::addMutInfo(const MutInfo& mI) {

  /* Take care of AA, pstn, settings if not base case */
  this->chn = "~"; // why "~": greatest value in ASCII code.
  this->pstn = 9999999;
  this->AA = "~";

  MutInfo* newMI = new MutInfo(mI);
  this->childMutInfo.push_back(newMI);

}

MutInfo& MutInfo::operator=(const MutInfo& mutInfo) {
  if (this == &mutInfo) {
    return *this;
  }
  else {
    
    for (vector<MutInfo*>::iterator itr = this->childMutInfo.begin();
	 itr != this->childMutInfo.end(); itr++) {
      delete *itr;
    }

    string mI_string = mutInfo.getString();
    this->init(mI_string);
    this->index = mutInfo.index;

    // Now add rCI info.
    this->_copy_rotConnInfo(mutInfo);


  } 

  return *this;
}

bool MutInfo::operator==(const MutInfo& mI) const {
  // If base MutInfo
  if (! this->isClusterMutInfo() ){

    if (this->chn == mI.chn and
	this->pstn == mI.pstn and
	this->AA == mI.AA) {
      return true;
    }
    else return false;
  }
  
  // If ClusterMutInfo
  else {
    
    string mI_string = mI.getString();
    string my_string = this->getString();
    
    if (mI_string == my_string) return true;
    else return false;

  }
}

bool MutInfo::operator<(const MutInfo& mI) const {
  if ( (! this->isClusterMutInfo()) and (! mI.isClusterMutInfo() ) ) {
    return this->_cmp_base_MutInfo(mI);
  }
  else {
      string mI_string = mI.getString();
      string my_string = this->getString();
      
      /* Below: isn't what i want for the order of my strings, but will do for now.  will recode later */
      
      int mI_level = this->_determineMaxLevel(mI_string);
      int my_level = this->_determineMaxLevel(my_string);
      
      if (my_level < mI_level) {
	return true;
      } 
      else if (mI_level == my_level) { 
	if (my_string < mI_string) {
	  return true;
	} 
	else 
	  return false;
      }
      else {
	return false;
      }
      
  }
}

void MutInfo::print_Me() const {
  cout << *this;
}

string MutInfo::getString() const {
  if ( ! this->isClusterMutInfo() ) {
    return this->str;
  }
  else {
    return this->_getString_for_Cluster();
  }

}

vector<MutInfo*> MutInfo::getAllMutInfos() {
  vector<MutInfo*> mI_list; mI_list.clear();

  /* base case */
  if (! this->isClusterMutInfo()) {
     mI_list.push_back(this); 
     return mI_list ;
  }

 /* traversing */
  else {
    for (vector<MutInfo*>::iterator itr = this->childMutInfo.begin(); 
	 itr != this->childMutInfo.end(); itr++) {
      vector< MutInfo* > childList = (*itr)->getAllMutInfos();
      mI_list.insert(mI_list.end(), childList.begin(), childList.end());
    }
  }

  return mI_list;

}


map<std::string, RotConnInfo*> MutInfo::getMutInfoStringWithRotConnInfo() const {

  /* Used for copy constructor and operator=() */

  map<std::string, RotConnInfo*> str_rCI_map; str_rCI_map.clear();

  /* base case */
  if (! this->isClusterMutInfo()) {
    if (this->getRotConnInfo() == NULL) {
      return str_rCI_map;
    }
    else {
      string self_string = this->getString();
      str_rCI_map[self_string] = this->getRotConnInfo();
    }
  }

  /* traversing */
  else {
    for (vector<MutInfo*>::const_iterator itr = this->childMutInfo.begin(); 
	 itr != this->childMutInfo.end(); itr++) {
      map<string, RotConnInfo*> childMap = (*itr)->getMutInfoStringWithRotConnInfo();
      str_rCI_map.insert(childMap.begin(), childMap.end());

    }
  }

  return str_rCI_map;

}


void MutInfo::searchAndAddRotConnInfo(MutInfo mI, RotConnInfo* rCI) {

  vector<MutInfo*> mI_list = this->getAllMutInfos();
  for (vector<MutInfo*>::iterator itr = mI_list.begin(); itr != mI_list.end(); itr++) {
    if ( (*(*itr)) == mI ) {
      (*itr)->setRotConnInfo(rCI);
    }
  }

}


void MutInfo::_init_base_MutInfo(string s) {
  
  this->str = s;

  this->AA = s.substr(0,1);
  int underscore_i = s.find("_");
  this->pstn = atoi(s.substr(1,underscore_i-1).c_str());

  this->chn = s.substr(underscore_i+1,1);

  this->childMutInfo.clear();

}


int MutInfo::_determineMaxLevel(string s) const {
  int max_level = 0;
  int current_level = 0;
  /* determine max_level */
  for (int i = 1; i != s.size()-1; i++) {  // i =1: starts from second character.
    if (s.substr(i,1) == "|" and s.substr(i-1,1) != "|" ) {
      current_level = 1;
    } else if (s.substr(i,1) == "|" and s.substr(i-1,1) == "|") {
      current_level++;
    } else {
      current_level = 0;
    }

    if (max_level < current_level) {
      max_level = current_level;
    }

  }

  return max_level;

}


bool MutInfo::_cmp_base_MutInfo(const MutInfo& mutInfo) const {
  // e.g.: C123_A comes before D123_A, after S122_A, before A1_B.  chn name is first checked.
  if (this->chn < mutInfo.chn) {  // chn
    return true;
  } else if (this->chn == mutInfo.chn) 
    { // chn
      if (this->pstn < mutInfo.pstn) { // pstn
	return true;
      } else if (this->pstn == mutInfo.pstn) 
	{  // pstns
	  if (this->AA < mutInfo.AA) {  // AA
	    return true;
	  } else { // AA
	    return false;
	  }
	} 
      else {   // pstn
	return false;
      }
    } 
  else 
    return false; // chn
}


string MutInfo::_getString_for_Cluster() const {
  
  /* returns in A1_A|B1_A||C1_A format */

  vector<string> same_level_strings;
  int max_level = 0;

  if (this->childMutInfo.size() == 1) {
    /* base case */
    max_level = 0;
    same_level_strings.push_back(this->childMutInfo[0]->getString());
  }
  else {
    for (vector<MutInfo*>::const_iterator itr = this->childMutInfo.begin();
	 itr != this->childMutInfo.end(); itr++) {
      
      string str = (*itr)->getString();
      same_level_strings.push_back(str);
      
      int new_string_level = this->_determineMaxLevel(str);
      
      max_level = (new_string_level > max_level) ? new_string_level : max_level;
      
    }
    max_level++; // need to increment max_level if not base case.
  }

  /* concatenate the strings */
  string sep(max_level, '|');
  string str;
  for (vector<string>::const_iterator itr = same_level_strings.begin(); itr != same_level_strings.end(); itr++) {
       str += sep;
       str += *itr;
  }
  str = str.substr(max_level,str.size());

  return str;

}

void MutInfo::_copy_rotConnInfo(const MutInfo& mutInfo) {
  // copy over RotConnInfo* stuff.

  map<string, RotConnInfo*> str_rCI_map = mutInfo.getMutInfoStringWithRotConnInfo();
  for (map<string, RotConnInfo*>::const_iterator itr = str_rCI_map.begin();
       itr != str_rCI_map.end(); itr++) {
    MutInfo mI(itr->first);
    this->searchAndAddRotConnInfo(mI, itr->second );
  
  }
  this->rCI = mutInfo.rCI; // should be NULL if MutInfo is a ClusterMutInfo.  else, already taken care of.
}

MutInfoPair::MutInfoPair() {

}

MutInfoPair::MutInfoPair(const MutInfoPair& mP) {
  //  cout << "MutInfoPair MutInfoPair& constructor" << endl;
  //this->init(mP.getMutInfo1(), mP.getMutInfo2());
  this->mutInfo1 = mP.getMutInfo1(); 
  this->mutInfo2 = mP.getMutInfo2();
  this->clashE = mP.getClashE();
  //cout << "after MutInfoPair & constructor" << endl;

}

MutInfoPair::MutInfoPair(MutInfo m1, MutInfo m2) {
  this->init(m1, m2);

}

MutInfoPair::MutInfoPair(MutInfo m1, MutInfo m2, double E) : clashE(E) {
  this->init(m1, m2);

}


MutInfoPair::~MutInfoPair() {

}

void MutInfoPair::init(MutInfo m1, MutInfo m2) {
  if (m1 == m2) {
    cerr << "the two MutInfos are identical! exiting." << endl;
    exit(2);
  }
  if (m1 < m2) {
    this->mutInfo1 = m1; 
    this->mutInfo2 = m2;
  } else {

    this->mutInfo2 = m1;
    this->mutInfo1 = m2;

  }

}

MutInfoPair& MutInfoPair::operator=(const MutInfoPair& mP) {

  if (this == &mP) 
    return *this;
  this->mutInfo1 = mP.getMutInfo1();
  this->mutInfo2 = mP.getMutInfo2();
  this->clashE = mP.getClashE();
  return *this;

}

bool MutInfoPair::operator==(const MutInfoPair& mP) const {
  if ((this->mutInfo1 == mP.mutInfo1 and this->mutInfo2 == mP.mutInfo2) or
      (this->mutInfo1 == mP.mutInfo2 and this->mutInfo2 == mP.mutInfo1) ) {  // this is unnecessary as each MutInfoPair structure is ordered.
    return true;
  } else {
    return false;
  }
      

}

bool MutInfoPair::operator<(const MutInfoPair& mP) const {

  if (this->getClashE() < mP.getClashE()) {
    return true;
  } else {
    return false;
  }

//   if (this->mutInfo1 < mP.mutInfo1) return true; // no need to cross check since MutInfoPairs are ordered by mutInfo1 and mutInfo2
//   if (this->mutInfo2 < mP.mutInfo2) return true;
//   return false;

}

ostream& operator<<(ostream& os, const MutInfo& m) {
  string m_s = m.getString();
  os << m_s;
  return os;
}

/*
ClusterMutInfo::ClusterMutInfo() {
  // Note: making ~ and 99999 to ensure when it comes to ordering ClusterMutInfo will alwys come behind a simple MutInfo
  this->chn = "~";
  this->pstn = 99999999;
  this->AA = "~";
}

ClusterMutInfo::ClusterMutInfo(string s) {

  this->init(s);

}

ClusterMutInfo::ClusterMutInfo(ClusterMutInfo& cMI) {
  
  string cMI_string = cMI.getString();
  this->init(cMI_string);

}

ClusterMutInfo::ClusterMutInfo(vector<string> v_s) {

  this->childClusterMutInfo.clear();
  for (vector<string>::iterator itr = v_s.begin(); itr != v_s.end(); itr++) {
    ClusterMutInfo* new_mI = new ClusterMutInfo(*itr);
    this->childClusterMutInfo.push_back(new_mI);
  }
  
}

ClusterMutInfo::~ClusterMutInfo() {

  for (vector<MutInfo*>::iterator itr = this->childClusterMutInfo.begin();
       itr != this->childClusterMutInfo.end(); itr++) {
    delete (*itr);
  }

}

ClusterMutInfo& ClusterMutInfo::operator=(const MutInfo& mI) {

  string mI_string = mI.getString();
  this->init(mI_string);

}

bool ClusterMutInfo::operator==(const MutInfo& mI) const {

  string mI_string = mI.getString();
  string my_string = this->getString();
  
  if (mI_string == my_string) return true;
  else return false;

}

bool ClusterMutInfo::operator<(const MutInfo& mI) const {

  string mI_string = mI.getString();
  string my_string = this->getString();
 
  // Below: isn't what i want for the order of my strings, but will do for now.  will recode later

  int mI_level = this->_determineMaxLevel(mI_string);
  int my_level = this->_determineMaxLevel(my_string);

  if (my_level < mI_level) {
    return true;
  } 
  else if (mI_level == my_level) { 
    if (my_string < mI_string) {
      return true;
    } 
    else 
      return false;
  }
  else {
    return false;
  }

}

void ClusterMutInfo::init(string s) {
  this->chn = "~";
  this->pstn = 9999999;
  this->AA = "~";

  // for a string of: A|B||C|||D|E|F||G,
   //  the object is to setup a tree that represents this structure.
//     Recursive algorithm is the easiest here.
   
  
  int max_level = this->_determineMaxLevel(s);

  // Base case. 

  if (max_level == 0) {
    MutInfo* mutInfo = new MutInfo(s);
    this->childClusterMutInfo.push_back(mutInfo);
    return;
  }


  // Make string ||||...|   <- max_level times. 
  string sep = string(max_level, '|');
  vector<string> split_s; split_s.clear();

  int last_mark = 0;
  for (int i = 0; i != s.size(); i++) {
    if (s.substr(i,max_level) == sep) {
      split_s.push_back(s.substr(last_mark,i-last_mark));
      last_mark = i+max_level;
    }
    if (i == s.size()-1) {
      split_s.push_back(s.substr(last_mark, i-last_mark+1));
    }
  }

  // Recurse.
  for (vector<string>::iterator itr = split_s.begin();
       itr != split_s.end(); itr++) {
    MutInfo* newMI = new ClusterMutInfo(*itr);
    this->childClusterMutInfo.push_back(newMI);
  }

  return;

}

void ClusterMutInfo::print_Me() const {

  // If input is: A|B||C|D
  // prints like:
  // [[A][B]][[C][D]]
  // {{A}{B}}{{C}{D}}
  //

  // Base case
//   if (this->childClusterMutInfo.size() == 1) {
//     cout << "[" << endl;
//     cout << *(this->childClusterMutInfo[0]) << endl;
//     cout << "" << endl;
//   }

  // Recurse 
  //cout << "size of list" << this->childClusterMutInfo.size() << endl;
  cout << "[";
  for (vector<MutInfo*>::const_iterator itr = this->childClusterMutInfo.begin();
       itr != this->childClusterMutInfo.end(); itr++) {
    (*itr)->print_Me();
  }
  cout << "]";


}

string ClusterMutInfo::getString() const {

  // returns in A1_A|B1_A||C1_A format

  vector<string> same_level_strings;
  int max_level = 0;

  if (this->childClusterMutInfo.size() == 1) {
  // base case 
    max_level = 0;
    same_level_strings.push_back(this->childClusterMutInfo[0]->getString());
  }
  else {
    for (vector<MutInfo*>::const_iterator itr = this->childClusterMutInfo.begin();
	 itr != this->childClusterMutInfo.end(); itr++) {
      
      string str = (*itr)->getString();
      same_level_strings.push_back(str);
      
      int new_string_level = this->_determineMaxLevel(str);
      
      max_level = (new_string_level > max_level) ? new_string_level : max_level;
      
    }
    max_level++; // need to increment max_level if not base case.
  }

  // concatenate the strings
  string sep(max_level, '|');
  string str;
  for (vector<string>::const_iterator itr = same_level_strings.begin(); itr != same_level_strings.end(); itr++) {
       str += sep;
       str += *itr;
  }
  str = str.substr(max_level,str.size());

  return str;
  

}


vector< MutInfo* > ClusterMutInfo::getAllMutInfos() {

  vector< MutInfo* > mI_list; mI_list.clear();

  // base case 
  if (this->childClusterMutInfo.size() == 1 ) {
    mI_list.push_back(this->childClusterMutInfo[0]);
    return mI_list;
  }
  // traversing 
  else {
    for (vector<MutInfo*>::iterator itr = this->childClusterMutInfo.begin(); 
	 itr != this->childClusterMutInfo.end(); itr++) {
      vector< MutInfo* > childList = (*itr)->getAllMutInfos();
      mI_list.insert(mI_list.end(), childList.begin(), childList.end());
    }
  }
  return mI_list;

}

*/
