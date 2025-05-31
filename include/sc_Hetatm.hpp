#ifndef SC_HETATM_HPP
#define SC_HETATM_HPP


#include <map>
#include <vector>

using namespace std;

#include "scream_atom.hpp"
#include "sc_ProteinComponent.hpp"

class Hetatm : public ProteinComponent {

public:

  Hetatm();
  Hetatm(const Hetatm&);

  Hetatm(const vector<SCREAM_ATOM*>&);
  ~Hetatm();

  vector<SCREAM_ATOM*> getAtomList() const;

  void print_Me() const;
  void append_to_filehandle(ostream* ofstream_p) const;
  void append_to_ostream_connect_info(ostream* ofstream_p) const;

  ProteinComponent* copy() const;

  string whatAmI() const { return string("Hetatm"); };
  
private:

  bool hetatm_atoms_on_free_store;
  multimap<string, SCREAM_ATOM*> hetatm_mm;

};


#endif /* SC_HETATM_HPP */
