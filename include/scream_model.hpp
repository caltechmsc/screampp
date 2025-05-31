/**\file scream_model
 *
 */

#ifndef SCREAM_MODEL
#define SCREAM_MODEL

#include "sc_Protein.hpp"
#include "scream_EE.hpp"
#include "Rotlib.hpp"
#include "scream_ctl_reader.hpp"
#include "sc_bgf_handler.hpp"

#include <string>

using namespace std;

class ScreamModel {

public:
  //vcvicek commented ScreamModel();
  ScreamModel(std::string); ///< Specify scream_ctl file.
  ~ScreamModel();

  
  ScreamParameters scream_parameters;
  bgf_handler HANDLER;
  Protein ptn;
  Scream_EE scream_EE;

  //! Add new constructs to scream_model world.

  Scream_EE* new_ScreamEE();
  Rotlib* new_Rotlib();

  vector<Scream_EE*> scream_EE_list;
  vector<Rotlib*> rotlib_list;

private:

  void _initScreamEE(); ///< not implemented to completion yet.  probably gonna leave this to python, for flexibility.
  string _convertDesignPositionToMutInfoName(string); ///< takes in a DesignName (like A32, which means chain A position 32), returns a legal mutInfo name (like A32_A, alanine, position 32, chain A.  Alanine is set as the default.)



};

#endif /* SCREAM_MODEL */
