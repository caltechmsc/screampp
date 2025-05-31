#include "defs.hpp"
#include "MutInfo.hpp"
#include "scream_ctl_reader.hpp"
#include "scream_tools.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

ScreamParameters::ScreamParameters(string ctl_name) {

  this->_init_default_params();
  this->read_scream_par_file(ctl_name);

}

std::string ScreamParameters::getKeepOriginalRotamer() {
  return this->KeepOriginalRotamer;
}

std::string ScreamParameters::getUseScreamEnergyFunction() {
  return this->UseScreamEnergyFunction;
}

vector<std::string> ScreamParameters::getDesignPositions() {
  
  vector<std::string> DesignPositions;
  for (stringVConstItr itr = this->DesignPositionAndClass.begin(); itr != this->DesignPositionAndClass.end(); ++itr) {
    vector<std::string> ss;
    split(*itr, "_", ss);
    DesignPositions.push_back(ss[0]);
  }

  return DesignPositions;

}

std::string ScreamParameters::getDesignClassFromPosition(std::string Chain_and_Position) {
  Chain_and_Position = scream_tools::strip_whitespace(Chain_and_Position);
  for (stringVConstItr itr = this->DesignPositionAndClass.begin(); itr != this->DesignPositionAndClass.end(); ++itr) {
    stringV ss;
    split(*itr, "_", ss);
    if (ss[0] == Chain_and_Position) {
      return ss[1];
    }
  }
}

vector<std::string> ScreamParameters::getDesignClassAAs(std::string Class) {

  Class = scream_tools::strip_whitespace(Class);
  vector<std::string> AAs;
  for (stringVConstItr itr = this->DesignAAClassDefns.begin(); itr != this->DesignAAClassDefns.end(); ++itr) {
    stringV ss;
    split(*itr, "_", ss);
    if (ss[0] == Class) {
      for (int i = 0; i != ss[1].length(); ++i) {
	AAs.push_back(string(ss[1].substr(i,1)));
      }
    }
  }
  return AAs;
}

void ScreamParameters::read_scream_par_file(string par_file) {

  ifstream PAR_FILE;
  PAR_FILE.open(par_file.c_str());
  if (PAR_FILE.bad()) {
    cerr << "in read_scream_par_file: Unable to open " << par_file << endl;
    exit(8);
    // return false;
  }

  string line_holder;

  if (PAR_FILE.fail()) {
    cout << "Can't open " << par_file << endl;
    cout << "will fail" << endl;
    exit(0);
  }

  while (!PAR_FILE.eof()) {

    //    PAR_FILE.getline(line_char, sizeof(line_char));
    getline(PAR_FILE, line_holder);

    string line = line_holder;

    vector<string> fields;
    fields = scream_tools::split(line);

    if (fields.empty()) 
      continue;

    if (fields[0] == "InputFileName") {
      if (fields.size() == 2) 
	this->InputFileName = fields[1];

    } else if (fields[0] == "MutateResidueInfo") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->MutateResidueInfo.push_back(*itr);
	  ++itr;
	}
      }	  

    } else if (fields[0] == "AdditionalLibraryInfo") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->AdditionalLibraryInfo.push_back(*itr);
	  ++itr;
	}
      }
    } else if (fields[0] == "Library") {
      if (fields.size() == 2) 
	this->Library = fields[1];
    
    } else if (fields[0] == "PlacementMethod") {
      if (fields.size() == 2)
	this->PlacementMethod = fields[1];

    } else if (fields[0] == "CreateCBParameters") {
      if (fields.size() == 5) {
	this->CreateCBParameters.clear();
	for (int i=1; i<5; i++) {
	  double value;
	  from_string(value, fields[i], std::dec);
	  this->CreateCBParameters.push_back(value);
	}
      }
      else
	cerr << "CreateCBParameters has wrong number of parameters!  Using default values." << endl;

    } else if (fields[0] == "KeepOriginalRotamer") {
      if (fields.size() == 2)
	this->KeepOriginalRotamer = fields[1];

    } else if (fields[0] == "UseScreamEnergyFunction") {
      if (fields.size() == 2)
	this->UseScreamEnergyFunction = fields[1];

    } else if (fields[0] == "UseDeltaMethod") {
      if (fields.size() == 2)
	this->UseDeltaMethod = fields[1];

    } else if (fields[0] == "UseRotamerNeighborList") {
      if (fields.size() == 2) 
	this->UseRotamerNeighborList = fields[1];

    } else if (fields[0] == "UseAsymmetricDelta") {
      if (fields.size() == 2)
	this->UseAsymmetricDelta = fields[1];

    } else if (fields[0] == "UseDeltaForInterResiE") {
      if (fields.size() == 2)
	this->UseDeltaForInterResiE = fields[1];

    } else if (fields[0] == "FlatDeltaValue") {
      if (fields.size() == 2)
	this->FlatDeltaValue = atof(fields[1].c_str());

    } else if (fields[0] == "DeltaStandardDevs") {
      if (fields.size() == 2)
	this->DeltaStandardDevs = atof(fields[1].c_str());

    } else if (fields[0] == "InnerWallScalingFactor") {
      if (fields.size() == 2)
	this->InnerWallScalingFactor = atof(fields[1].c_str());

    } else if (fields[0] == "NonPolarHCalc") {
      if (fields.size() == 2)
	this->NonPolarHCalc = fields[1];

    } else if (fields[0] == "ScoringFunction") {
      if (fields.size() == 2) 
	this->ScoringFunction = fields[1];

    } else if (fields[0] == "MultiplePlacementMethod") {
      if (fields.size() == 2)
	this->MultiplePlacementMethod = fields[1];

    } else if (fields[0] == "CBGroundSpectrumCalc") {
      if (fields.size() == 2)
	this->CBGroundSpectrumCalc = fields[1];

    } else if (fields[0] == "OneEnergyFFParFile") {
      if (fields.size() == 2) 
	this->OneEnergyFFParFile = fields[1];
    
    } else if (fields[0] == "DeltaParFile") {
      if (fields.size() == 2)
	this->DeltaParFile = fields[1];

    } else if (fields[0] == "EachAtomDeltaFile") {
      if (fields.size() == 2)
	this->EachAtomDeltaFile = fields[1];
      
    } else if (fields[0] == "PolarOptimizationExclusions") {
      if (fields.size() == 2)
	this->PolarOptimizationExclusions = fields[1];

    } else if (fields[0] == "LJOption") {
      if (fields.size() == 2)
	this->LJOption = fields[1];

    } else if (fields[0] == "CoulombMode") {
      if (fields.size() == 2)
	this->CoulombMode = fields[1];

    } else if (fields[0] == "Dielectric") {
      if (fields.size() == 2)
	this->Dielectric = atof(fields[1].c_str());

    } else if (fields[0] == "Selections") {
      if (fields.size() == 2)
	this->Selections = atoi(fields[1].c_str());

    } else if (fields[0] == "MaxSearchNumber") {
      if (fields.size() == 2)
	this->MaxSearchNumber = atoi(fields[1].c_str());

    } else if (fields[0] == "AbsStericClashCutoffEL") {
      if (fields.size() == 2)
	this->AbsStericClashCutoffEL = atof(fields[1].c_str());

    } else if (fields[0] == "StericClashCutoffEnergy") {
      if (fields.size() == 2)
	this->StericClashCutoffEnergy = atof(fields[1].c_str());

    } else if (fields[0] == "StericClashCutoffDist") {
      if (fields.size() == 2) 
	this->StericClashCutoffDist = atof(fields[1].c_str());
    
    } else if (fields[0] == "MaxFinalStepRunTime") {
      if (fields.size() == 2)
	this->MaxFinalStepRunTime = atoi(fields[1].c_str());

    } else if (fields[0] == "LibPath") {
      if (fields.size() == 2) {
	this->LibPath = fields[1];
      }

    } else if (fields[0] == "Verbosity") {
      if (fields.size() == 2)
	this->Verbosity = atoi(fields[1].c_str());

    } else if (fields[0] == "DesignPositionAndClass") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->DesignPositionAndClass.push_back(*itr);
	  ++itr;
	}
      }	  

    } else if (fields[0] == "DesignAAClassDefns") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->DesignAAClassDefns.push_back(*itr);
	  ++itr;
	}
      }	  

    } else if (fields[0] == "JustOutputSequence") {
      if (fields.size() == 2)
	this->JustOutputSequence = fields[1];

    } else if (fields[0] == "BindingSiteMode") {
      if (fields.size() == 2) 
	this->BindingSiteMode = fields[1];

    } else if (fields[0] == "FixedResidues") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->FixedResidues.push_back(*itr);
	  ++itr;
	}
      }	  

    } else if (fields[0] == "AroundAtom") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->AroundAtom.push_back(atoi((*itr).c_str()));
	  ++itr;
	}
      }

    } else if (fields[0] == "AroundResidue") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->AroundResidue.push_back(MutInfo(*itr));
	  ++itr;
	}
      }	 

    } else if (fields[0] == "AroundChain") {
      vector<string>::const_iterator itr = fields.begin();
      if (fields.size() >= 2) {
	++itr;
	while (itr != fields.end() ) {
	  this->AroundChain.push_back(*itr);
	  ++itr;
	}
      }	 

    } else if (fields[0] == "AroundDistance") {
      if (fields.size() == 2) 
	this->AroundDistance = atof(fields[1].c_str());

    } else if (fields[0] == "AroundDistanceDefn") {
      if (fields.size() == 2)
	this->AroundDistanceDefn = fields[1];

    } else {
      // do nothing
    }
  }

  this->print_to_output(&cout);
  
  PAR_FILE.close();

}

void ScreamParameters::print_to_output(ostream* ostream_p) const {

  *ostream_p << "InputFileName           " << this->InputFileName << endl;

  *ostream_p << "MutateResidueInfo       ";
  for (vector<string>::const_iterator itr = this->MutateResidueInfo.begin(); 
       itr != this->MutateResidueInfo.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  *ostream_p << "AdditionalLibraryInfo   ";
  for (vector<string>::const_iterator itr = this->AdditionalLibraryInfo.begin(); 
       itr != this->AdditionalLibraryInfo.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  string CBCreateParamString = "";
  for (int i=0; i < this->CreateCBParameters.size(); i++) {
    CBCreateParamString += stringify(this->CreateCBParameters[i]);
    CBCreateParamString += " ";
  }


  *ostream_p << "Library                 " << this->Library << endl;
  *ostream_p << "PlacementMethod         " << this->PlacementMethod << endl;
  *ostream_p << "CreateCBParameters      " << CBCreateParamString << endl; // CBCreateParamString variable defined a few lines above.
  *ostream_p << "UseScreamEnergyFunction " << this->UseScreamEnergyFunction << endl;
  *ostream_p << "UseDeltaMethod          " << this->UseDeltaMethod << endl;
  *ostream_p << "UseRotamerNeighborList  " << this->UseRotamerNeighborList << endl;
  *ostream_p << "UseAsymmetricDelta      " << this->UseAsymmetricDelta << endl;
  *ostream_p << "UseDeltaForInterResiE   " << this->UseDeltaForInterResiE << endl;
  *ostream_p << "FlatDeltaValue          " << this->FlatDeltaValue << endl;
  *ostream_p << "DeltaStandardDevs       " << this->DeltaStandardDevs << endl;
  *ostream_p << "InnerWallScalingFactor  " << this->InnerWallScalingFactor << endl;
  *ostream_p << "NonPolarHCalc           " << this->NonPolarHCalc << endl;
  *ostream_p << "ScoringFunction         " << this->ScoringFunction << endl;
  *ostream_p << "OneEnergyFFParFile      " << this->OneEnergyFFParFile << endl;
  *ostream_p << "DeltaParFile            " << this->DeltaParFile << endl;
  *ostream_p << "EachAtomDeltaFile       " << this->EachAtomDeltaFile << endl;
  *ostream_p << "PolarOptimizationExclusions " << this->PolarOptimizationExclusions << endl;
  *ostream_p << "LJOption                " << this->LJOption << endl;
  *ostream_p << "CoulombMode             " << this->CoulombMode << endl;
  *ostream_p << "Dielectric              " << this->Dielectric << endl;
  *ostream_p << "MultiplePlacementMethod " << this->MultiplePlacementMethod << endl;
  *ostream_p << "CBGroundSpectrumCalc    " << this->CBGroundSpectrumCalc << endl;
  *ostream_p << "Selections              " << this->Selections << endl;
  *ostream_p << "MaxSearchNumber         " << this->MaxSearchNumber << endl;
  *ostream_p << "AbsStericClashCutoffEL  " << this->AbsStericClashCutoffEL << endl;
  *ostream_p << "StericClashCutoffEnergy " << this->StericClashCutoffEnergy << endl;
  *ostream_p << "StericClashCutoffDist   " << this->StericClashCutoffDist  << endl;
  *ostream_p << "MaxFinalStepRunTime     " << this->MaxFinalStepRunTime << endl;
  *ostream_p << "LibPath                 " << this->LibPath << endl;
  *ostream_p << "Verbosity               " << this->Verbosity << endl;
  
  *ostream_p << "DesignPositionAndClass  " ;
  for (vector<string>::const_iterator itr = this->DesignPositionAndClass.begin(); 
       itr != this->DesignPositionAndClass.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  *ostream_p << "DesignAAClassDefns      ";
  for (vector<string>::const_iterator itr = this->DesignAAClassDefns.begin(); 
       itr != this->DesignAAClassDefns.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  *ostream_p << "JustOutputSequence      " << this->JustOutputSequence << endl;

  *ostream_p << "BindingSiteMode         " << this->BindingSiteMode << endl;

  *ostream_p << "FixedResidues           " ;
  for (vector<string>::const_iterator itr = this->FixedResidues.begin(); 
       itr != this->FixedResidues.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  *ostream_p << "AroundAtom              " ;
  for (vector<int>::const_iterator itr = this->AroundAtom.begin(); 
       itr != this->AroundAtom.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  *ostream_p << "AroundResidue           " ;
  for (vector<MutInfo>::const_iterator itr = this->AroundResidue.begin(); 
       itr != this->AroundResidue.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  *ostream_p << "AroundChain             " ;
  for (vector<string>::const_iterator itr = this->AroundChain.begin(); 
       itr != this->AroundChain.end();
       ++itr) {
    *ostream_p << (*itr) << " ";
  }
  *ostream_p << endl;

  *ostream_p << "AroundDistance          " << this->AroundDistance << endl;
  *ostream_p << "AroundDistanceDefn      " << this->AroundDistanceDefn << endl;

  *ostream_p << endl;



}

string ScreamParameters::minimizationMethod() const {

  if (ScoringFunction == "Default") {
    return "V";
  } else {
    return ScoringFunction.substr(0,1);
  }

}

int ScreamParameters::minimizationSteps() const {

  if (ScoringFunction == "Default") {
    return 5;
  } else {
    int i = ScoringFunction.rfind("OE");
    return atoi(ScoringFunction.substr(1, i-1).c_str());
  }

}


string ScreamParameters::oneEMethod() const {

  if (ScoringFunction == "Default") {
    return "V";
  }
  else {
    //    string::size_type i = ScoringFunction.rfind("OE");
    int i = ScoringFunction.rfind("OE");
    return ScoringFunction.substr(i+2, 1);
  }

}

string ScreamParameters::residueToScreamName(string mutInfo) const {
  
  return mutInfo.substr(0,1);

}

int ScreamParameters::residueToScreamPstn(string mutInfo) const {

  int i = mutInfo.rfind("_");
  return atoi(mutInfo.substr(1, i-1).c_str());
  
}

string ScreamParameters::residueToScreamChn(string mutInfo) const {

  int i = mutInfo.rfind("_");
  return mutInfo.substr(i+1, 1);

}

string ScreamParameters::determineLibDirPath() const {
  // Reads in environment variables set in SCREAM_OLD_LIB_PATH, SCREAM_LIB_PATH, SCWRL_LIB_PATH

  //vcvicek
  //char * SCREAM_OLD_LIB_PATH_c = getenv("SCREAM_OLD_LIB_PATH");
  //char * SCREAM_LIB_PATH_c = getenv("SCREAM_LIB_PATH");
  //char * SCWRL_LIB_PATH_c = getenv("SCWRL_LIB_PATH");
  char * SCREAM_NEW = getenv("SCREAM_NEW");
  if (SCREAM_NEW == NULL) {printf("error: enviromental variable SCREAM_NEW is not set \n"); exit(1);}

  char * SCREAM_NEW_LIB = getenv("SCREAM_NEW_LIB");
  if (SCREAM_NEW_LIB == NULL) {printf("error: environmental variable SCREAM_NEW_LIB is not set \n"); exit(1);}
  string SCREAM_LIB_PATH=string(SCREAM_NEW_LIB);
  
  //string SCREAM_OLD_LIB_PATH=string(SCREAM_NEW)+"/scream/lib/libClusteredOld/";
  //string SCWRL_LIB_PATH=string(SCREAM_NEW)+"/scream/lib/libScwrl";
 
  //if (this->Library == "SCWRL") {
  //  return SCWRL_LIB_PATH;
  //} 
  //if ((this->Library).substr(0,1) == "H") {
  //  // If Library is of name H01, H02, etc. up to H50.
  //  return SCREAM_OLD_LIB_PATH;
  //}

  if ((this->Library).substr(0,1) == "V") {
    // This set of libraries has the "right energies" (internal energies, like torsion, included)
    if (this->LibPath != "")  // default value for LibPath is "". 
      return this->LibPath;
    else
      return SCREAM_LIB_PATH;
  }

  if (this->Library == "USER") {
    return this->LibPath;
  }

  // Else, user specifies the library path.
  return this->LibPath;

}

string ScreamParameters::determineLibDirFileNameSuffix() const {
  // This is determined by the field this->Library.  Follows original perl SCP format: I.e. H10 is the 1.0A library, etc.
  // But now, only the the number is read in and used, until the string is "SCWRL".

  string suffix = "";

  if (this->Library == "SCWRL") {
    suffix = string("_10.lib");
  }
  else if (this->Library == "USER") {
    suffix = string(".lib");
  }
  else if (this->Library == "") {
    suffix = string("_10.lib");
  }
  else {   // Only SCWRL is doesn't have any numbers in it.
    string resolution = this->Library.substr(1,2);
    suffix = "_" + resolution + ".lib";
  }

  return suffix;

}

string ScreamParameters::determineCnnDirPath() const {
  // Reads in environment variables set in SCREAM_OLD_LIB_PATH, SCREAM_LIB_PATH, SCWRL_LIB_PATH

  //vcvicek
  //char * SCREAM_OLD_CNN_PATH_c = getenv("SCREAM_OLD_CNN_PATH");
  //char * SCREAM_CNN_PATH_c = getenv("SCREAM_CNN_PATH");
  char * SCREAM_NEW_CNN = getenv("SCREAM_NEW_CNN");
  if (SCREAM_NEW_CNN == NULL) {printf("error: enviromental variable SCREAM_NEW_CNN is not set \n"); exit(1);}

  //string SCREAM_OLD_CNN_PATH = string(SCREAM_NEW)+"/lib/NtrlAARotConnOld/";
  string SCREAM_CNN_PATH = string(SCREAM_NEW_CNN);
 
  //if (this->Library == "SCWRL") {
  //  return SCREAM_CNN_PATH;
  //} 
  //if ((this->Library).substr(0,1) == "H") {
  //  // If Library is of name H01, H02, etc. up to H50.
  //  return SCREAM_OLD_CNN_PATH;
  //}

  if ((this->Library).substr(0,1) == "V") {
    // This set of libraries has the "right energies" (internal energies, like torsion, included)
    if (this->LibPath != "")
      return this->LibPath;
    else
      return SCREAM_CNN_PATH;
  }

}

int ScreamParameters::getLibResolution() const {

  if (this->Library == "SCWRL" ) {
    return 0;
  }
  else {
    int resolution = 0;
    resolution = atoi(this->Library.substr(1,2).c_str());
    return resolution;
  }

}

string ScreamParameters::returnEnergyMethod() const {
  string method = this->UseDeltaMethod;
  if ( this->UseAsymmetricDelta == "YES" )
    method += "_ASYM";
  if ( this->CBGroundSpectrumCalc == "NO" )
    method += "_NOCB";

  return method;

}

double ScreamParameters::returnEnergyMethodTValue() const {
  if ( this->UseDeltaMethod == "FLAT" )
    return this->FlatDeltaValue;
  else if ( this->UseDeltaMethod == "FULL" )
    return this->DeltaStandardDevs;
  else if ( this->UseDeltaMethod == "SCALE" )
    return this->InnerWallScalingFactor;

}



void ScreamParameters::_init_default_params() {

  this->InputFileName = string("");
  //this->AdditionalLibraryInfo = string("");
  this->Library = string("V10");

  this->PlacementMethod = "Default";
  this->CreateCBParameters.clear();
  CreateCBParameters.push_back(1.81); // offBisectorAngle
  CreateCBParameters.push_back(51.1); // offPlaceAngle
  CreateCBParameters.push_back(1.55); // bond length
  CreateCBParameters.push_back(0.5); // rotamerMatchVector

  this->KeepOriginalRotamer = string("YES");
  this->UseScreamEnergyFunction = string("YES");
  this->UseDeltaMethod = string("FLAT");
  this->UseRotamerNeighborList = string ("NO");
  this->UseAsymmetricDelta = string ("NO");
  this->UseDeltaForInterResiE = string("YES");
  this->FlatDeltaValue = 0.3;
  this->DeltaStandardDevs = -0.8;
  this->InnerWallScalingFactor = 0.3;
  this->NonPolarHCalc = string("YES");
  this->MultiplePlacementMethod = "ExcitationWithClustering";
  this->CBGroundSpectrumCalc = "YES";

  //vcvicek
  char * SCREAM_NEW = getenv("SCREAM_NEW");
  if (SCREAM_NEW == NULL) {printf("erorr: enviromental variable SCREAM_NEW is not set\n"); exit(1);}

  this->OneEnergyFFParFile = string(SCREAM_NEW)+"/lib/SCREAM_delta_par_files/dreidii322-scream-EE.par";
  this->DeltaParFile = string(SCREAM_NEW)+"/lib/SCREAM_delta_par_files/SCREAM_delta.par";
  this->EachAtomDeltaFile = string(SCREAM_NEW)+"/lib/SCREAM_delta_par_files/SCREAM_EachAtomDeltaFileStub.par";
  this->PolarOptimizationExclusions = string(SCREAM_NEW)+"/lib/SCREAM_delta_par_files/SCREAM_PolarOptimizationExclusionsStub.par";

  this->LJOption = "12-6";
  this->CoulombMode = "Normal";
  this->Dielectric = 2.5;
  
  this->Selections = 1;
  this->MaxSearchNumber = 99999999;
  
  this->AbsStericClashCutoffEL = 5000;
  this->StericClashCutoffEnergy = 300;
  this->StericClashCutoffDist = 2.0;

  this->MaxFinalStepRunTime = 36000;
  //this->LibPath = string("/project/Biogroup/Software/SCREAM/lib/v3.0/Clustered/lib/");
  //this->LibPath = string("project/Biogroup/Software/SCREAM/lib/new_clustering/sorted_lib/");
  this->LibPath = string(""); // empty string
  this->Verbosity = 0;
  this->DesignPositionAndClass.clear();
  this->DesignAAClassDefns.clear();
  
  this->JustOutputSequence = "NO";

  this->BindingSiteMode = "";
  this->FixedResidues.clear();
  this->AroundAtom.clear();
  this->AroundResidue.clear();
  this->AroundChain.clear();
  this->AroundDistance = 4.0;
  this->AroundDistanceDefn = "SideChainOnly";

  // Below: Obsolete stuff or stuff that must be provided in ctl file.
  this->ScoringFunction = "";

}

