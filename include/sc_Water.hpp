#ifndef SC_WATER_HPP
#define SC_WATER_HPP

#include "sc_ProteinComponent.hpp"
#include "scream_atom.hpp"

#include <vector>
#include <map>

using namespace std;

class Water : public ProteinComponent {

public:

  Water();
  Water(const vector<SCREAM_ATOM*>&);
  ~Water();
  
  vector<SCREAM_ATOM*> getAtomList() const;

  void print_Me() const;
  void append_to_filehandle(ostream* ofstream_p) const;
  void append_to_ostream_connect_info(ostream* ofstream_p) const;

  ProteinComponent* copy() const;
  string whatAmI() const { return string("Water"); };

private:
  
  bool water_atoms_on_free_store;
  multimap<string, SCREAM_ATOM*> water_mm;


};




#endif /* SC_WATER_HPP */
