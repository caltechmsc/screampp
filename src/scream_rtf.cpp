//vcvicek
#include <cstdlib>

#include "scream_rtf.hpp"
#include <iostream>
#include <fstream>

SCREAM_RTF::SCREAM_RTF() {

}

SCREAM_RTF::SCREAM_RTF(string rtf_file) {
  this->_init(rtf_file);
}

SCREAM_RTF::~SCREAM_RTF() {
  Debug debugInfo("SCREAM_RTF::~SCREAM_RTF()");
  for (map<string,AminoAcid_RTF*>::iterator itr = this->_rtfTable.begin();
       itr != this->_rtfTable.end(); itr++ ) {
    debugInfo.out(" Deleting a AminoAcid_RTF.");
    delete (*itr).second;
  }

}

AminoAcid_RTF* SCREAM_RTF::get_AminoAcid_RTF(string resName) {
  return this->_rtfTable[resName];
}

string SCREAM_RTF::get_ff_type(string resName, string atomLabel) {
  return this->_rtfTable[resName]->get_ff_type(atomLabel);
}

void SCREAM_RTF::_init(string rtf_file) {
  // FORMAT OF a SCREAM rtf file: exactly identical to a NAMD/CHARMM rtf file.
  Debug debugInfo("SCREAM_RTF::_init(string rtf_file)");

  ifstream FILE(rtf_file.c_str());
  if (!(FILE.good())) {
    cerr << "Unable to open " << rtf_file << endl;
    exit(2); return;
  }
  
  stringV resiLines;
  string crntRes = "";
  bool flag = false;

  string l;

  while(getline(FILE,l)) {

    /* Skip over comments. */
    if (l[0] == '!' or l[0] == '#' or l[0] == '*' or l[0] == ' ' or l == "") {
      continue;
    }
    debugInfo.out(l);
    /* read real data. */
    if (l.substr(0,4) == "RESI" ) {
      if (flag == false) {  // first time through
	flag = true;
	stringV fields;
	split(l, " ", fields);
	crntRes = fields[1];
      }
      else {
	AminoAcid_RTF* aaRTF = new AminoAcid_RTF(resiLines);
	debugInfo.out("newed aaRTF--" + crntRes);
	this->_rtfTable[crntRes] = aaRTF;

	resiLines.clear();
	stringV fields;
	split(l, " ", fields);
	crntRes = fields[1];
      }
    }
     
    if (flag == true) {
      resiLines.push_back(l);
    }
    
  }
  /* End condition */
  AminoAcid_RTF* aaRTF = new AminoAcid_RTF(resiLines);
  debugInfo.out("newed aaRTF--" + crntRes);
  this->_rtfTable[crntRes] = aaRTF;
   
  debugInfo.out("Done!");
}

AminoAcid_RTF::AminoAcid_RTF() {

}

AminoAcid_RTF::AminoAcid_RTF(stringV& lines) {
  this->_init(lines);
}

AminoAcid_RTF::~AminoAcid_RTF() {

  // nothing is needed to do.  yay!

}

string AminoAcid_RTF::get_ff_type(string atom_label) {
  return this->ff_type[atom_label];
}

void AminoAcid_RTF::_init(stringV& lines) {
  Debug debugInfo("AminoAcid_RTF::_init");
  for (stringVItr itr = lines.begin(); itr != lines.end(); itr++) {
    /* first line: RESI ALA    something */
    debugInfo.out(*itr);
    stringV fields;
    split(*itr, " ", fields);
    if (fields[0] == "RESI") {
      this->resName = fields[1];
    }
    /* ALIAS label */
    
    /* ignore "GROUP" label */
    
    /* ATOM label read */
    if (fields[0] == "ATOM") {
      this->ff_type[fields[1]] = fields[2];
      // charges: maybe?
    } else
    
    /* BOND or DOUBLE label read. */
    if (fields[0] == "BOND" or fields[0] == "DOUBLE") {
      // read two at a time.  should assert fields.size() >= 1.
      int i = 0;
      int n = fields.size();
      while (1+ 2*(i+1) <= n) {
	//this->bonds[fields[1+2*i]] = fields[1+(2*i+1)];
	this->bonds.insert(make_pair(fields[1+2*i], fields[1+2*i+1]));
	i++;
      }

    } 
    /* Ignore DONOR, ACCEPTOR, IC for now.*/
    // else if (fields[0] == "GROUP" or fields[0] == "IMPR" or fields[0] == "DONOR" or fields[0] == "ACCEPTOR" or fields[0] == "IC") {
//       continue;
//     }
    
    if (fields[0] == "PRES") {
      break;
    }
  }
  
}
