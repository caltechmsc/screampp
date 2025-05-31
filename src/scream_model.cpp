#include "scream_model.hpp"
#include "defs.hpp"
#include "MutInfo.hpp"
#include <string>

using namespace std;

//vcvicek
//ScreamModel::ScreamModel() {
//
//}

ScreamModel::ScreamModel(string scream_ctl_file) : scream_parameters(scream_ctl_file), HANDLER(scream_parameters.InputFileName), ptn(&(HANDLER.atom_list)) {

  //  this->_initScreamEE();  // now obsolete.  PYTHON script does this now... with all those init functions from Scream_EE.
  
  ptn.setPlacementMethod(scream_parameters.getPlacementMethod());
  vector<double> CBParameters = scream_parameters.getCreateCBParameters();
  ptn.setOffBisectorAngle(CBParameters[0]);
  ptn.setOffPlaneAngle(CBParameters[1]);
  ptn.setBondLength(CBParameters[2]);
  ptn.setRotamerMatchVectorLamdba(CBParameters[3]);
  cout << "THIS LINE INDICATES YOU USED THE GRIFFITH SCREAM VERSION" << endl;
}

ScreamModel::~ScreamModel() {
  cout << "Destroying ScreamModel! " << endl;

  for (vector<Scream_EE* >::iterator itr = this->scream_EE_list.begin();
       itr != this->scream_EE_list.end(); itr++) {
    delete (*itr);
  }

  for (vector<Rotlib*>::iterator itr = this->rotlib_list.begin();
       itr != this->rotlib_list.end(); itr++) {
    delete (*itr);
  }

  cout << "THIS LINE INDICATES YOU USED THE GRIFFITH SCREAM VERSION" << endl;
  cout << "ScreamModel destroyed." << endl;

}

Scream_EE* ScreamModel::new_ScreamEE() {

  Scream_EE* new_scream_EE = new Scream_EE();
  this->scream_EE_list.push_back(new_scream_EE);

  return new_scream_EE;

}

Rotlib* ScreamModel::new_Rotlib() {

  Rotlib* new_rotlib = new Rotlib();
  this->rotlib_list.push_back(new_rotlib);
  
  return new_rotlib;
  
}

void ScreamModel::_initScreamEE() {

  // This routine does not load the rotamer libraries.  Rotamer libraries are loaded in an upper level construct (usually in python, though this extra level of indirection seems more like bad design than anything else.)

  this->scream_EE = Scream_EE();
  stringV mutInfoList = this->scream_parameters.getMutateResidueInfoList();
  stringV additionalLib_list = this->scream_parameters.getAdditionalLibraryInfo();


  if (mutInfoList.size() == 1 and mutInfoList[0] == "DESIGN") {
    // setup protein design routines
    stringV designPositions = scream_parameters.getDesignPositionAndClass();
    for (stringVConstItr designInfo = designPositions.begin();
	 designInfo != designPositions.end(); designInfo++) {
      string mutInfoName = this->_convertDesignPositionToMutInfoName(*designInfo);
      cout << mutInfoName << endl;
      MutInfo mI(mutInfoName);
      this->scream_EE.addMutInfoRotConnInfo(mI);
    }

  } 

  else {

    if (mutInfoList.size() == 1 and mutInfoList[0] == "BINDING_SITE") {
      vector<MutInfo> mI_list;
      double dist = scream_parameters.getAroundDistance();
      string inclusiveMode = scream_parameters.getAroundDistanceDefn();
      if (scream_parameters.getBindingSiteMode() == "AroundAtom") {
	mI_list = this->ptn.residuesAroundAtomN(scream_parameters.getAroundAtom(), dist, inclusiveMode);
      } else if (scream_parameters.getBindingSiteMode() == "AroundResidue") {
	mI_list = this->ptn.residuesAroundResidue(scream_parameters.getAroundResidue(), dist, inclusiveMode);
      } else if (scream_parameters.getBindingSiteMode() == "AroundChain") {
	mI_list = this->ptn.residuesAroundChain(scream_parameters.getAroundChain(), dist, inclusiveMode);
      }
      // unfinished; unnecessary since scream_model no longer a black box for Python front.
    }

    // just residue placement.
    for (stringVConstItr mutInfo = mutInfoList.begin(); mutInfo != mutInfoList.end(); mutInfo++) {
      MutInfo mI(*mutInfo);
      this->scream_EE.addMutInfoRotConnInfo(mI);
      
    }
  }
  
  /* Remark: Following: energy calculations done on the fly.  If no mutation, slightly more efficient if do energy calculations by first setting up a list (the function scream_EE.init_after_addedMutInfoRotConnInfo() ) */
  if (this->scream_parameters.getUseRotamerNeighborList() == "YES") {
    this->scream_EE.init_after_addedMutInfoRotConnInfo_neighbor_list( &(this->ptn), &scream_parameters);
  } else {
    this->scream_EE.init_after_addedMutInfoRotConnInfo_on_the_fly_E( &(this->ptn), &scream_parameters);
  }

  /* Then, need to initialize DELTA values for FULL implementation. */
//   if (this->scream_parameters.UseDeltaMethod == "FULL") {
//     string lib_name = this->scream_parameters.Library;
//     string method = "FULL";
//     double alpha = this->scream_parameters.DeltaStandardDevs;
    
//     if (lib_name[0] == 'V') {
//       lib_name = lib_name.substr(1,2);
//     }
    
//     this->scream_EE.initScreamAtomDeltaValue(lib_name, method, alpha);
//   }

}

string ScreamModel::_convertDesignPositionToMutInfoName(string DesignName) {

  vector<string> splitDesignName;
  split(DesignName, "_", splitDesignName);

  string chain = splitDesignName[0].substr(0,1);
  string position = splitDesignName[0].substr(1);
  string mutInfo = "A" + position + "_" + chain;

  return mutInfo;
  

}
