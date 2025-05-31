/* sc_ProteinComponent.cpp
 *
 * Header file for classes relevant to proteins in SCREAM module.  
 * 
 * Copyright (c) 2003 Victor Wai Tak Kam.  All Rights Reserved.
 *
 */


#include "sc_ProteinComponent.hpp"

/************************\
ProteinComponent Defitions
\************************/
/*0
ProteinComponent* ProteinComponent::copy() const {

  return this;

}
*/
ProteinComponent& ProteinComponent::operator=(const ProteinComponent& pc) {

  return *this;

}

ProteinComponent& ProteinComponent::merge(ProteinComponent& pc) {

  return *this;

}

