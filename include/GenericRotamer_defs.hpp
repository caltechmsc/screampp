#ifndef GENERICROTAMER_DEFS_HPP
#define GENERICROTAMER_DEFS_HPP

#include "AARotamer.hpp"

namespace GenericRotamer {

  string GenericRotamerPath = "../lib/GenericRotamers/";
  const AARotamer* A = new AARotamer(GenericRotamerPath + "A.lib", GenericRotamerPath + "A.cnn");
  const AARotamer* C = new AARotamer(GenericRotamerPath + "C.lib", GenericRotamerPath + "C.cnn");
  
};


#endif /* GENERICROTAMER_DEFS_HPP */
