#ifndef SCREAM_RTF_HPP
#define SCREAM_RTF_HPP

#include "defs.hpp"
#include <map>

class AminoAcid_RTF {
public:
  AminoAcid_RTF();
  AminoAcid_RTF(stringV&);
  ~AminoAcid_RTF();

  string get_ff_type(string); ///< string: atom label of this residue.
  multimap<string, string> return_bonds_table() const {return this->bonds;}; ///< returns the bonds table.

private:
  
  string resName;
  void _init(stringV&);
  multimap<string, string> bonds; ///< first string: first atom label. second string: second atom label.  in labels: all unique.
  map<string, string> ff_type; ///< first string: first atom label.  second string: forcefield type.  

};

class SCREAM_RTF {
public:
  SCREAM_RTF();
  SCREAM_RTF(string); // reader
  ~SCREAM_RTF();

  AminoAcid_RTF* get_AminoAcid_RTF(string); ///< string: resName.
  string get_ff_type(string, string); ///< string1: resName. string2: atomLabel.

private:
  void _init(string);
  map<string, AminoAcid_RTF*> _rtfTable; 	///< Table.

};


#endif /* ifndef SCREAM_RTF_HPP */
