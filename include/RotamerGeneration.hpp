//! RotamerGeneration.hpp
/**
 * RotamerGeneration.hpp
 *
 * Header file for classes relevant to RotamerGeneration in SCREAM module.
 *
 * Copyright (c) 2004 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */
 
#ifndef ROTAMERGENERATION_HPP
#define ROTAMERGENERATION_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "defs.hpp"
#include "scream_atom.hpp"
#include "scream_vector.hpp"
#include "scream_matrix.hpp"
#include "sc_Protein.hpp"

#include "scream_helper_classes.hpp"

/**
 * RotamerGeneration: generates rotamer in the environment of a Protein (Protein: synonym for a system).
 */

class RotamerGeneration  {


public:
  RotamerGeneration();
  RotamerGeneration(Protein*, RotamerAxis*);


private:
  Protein* ptn;			///< Needs a pointer to Protein/system to know something about its environment.
  RotamerAxis* axis; ///< The axis that the rotamer is allowed to spin around on.  The two pointers to SCREAM_ATOMs here specifies the direction as well, and are always atoms in the protein system.
  

};


#endif /* ROTAMERGENERATION_HPP */

