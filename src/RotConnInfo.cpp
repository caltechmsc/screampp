#include "RotConnInfo.hpp" 

#include <iostream>
//vcvicek
#include <cstdlib>


using namespace std;

void RotConnInfo::modifyMappingInProteinAtoms(vector<int>& int_list) {

    ///< If protein atoms have their mappings modified, need to update RotConnInfo mappings.  This function is called to do that; vector<int> stores info about what has happened to protein atoms; with index being atom_n.  when index == 0, vector[0] == 0, i.e. a placeholder.

    // Should do more checks; skip checkings for now.

    for (map<int, int>::iterator itr = this->atom_n_map.begin(); itr != this->atom_n_map.end(); itr++) {
      // checks
      int rot_rep = itr->first;
      int ptn_rep = itr->second;

      if (ptn_rep >= int_list.size() ) {
	cout << "Serious error: mutation in protein has wiped out certain atoms .cnn files specification.  Please fix your .cnn files." << endl;
	cerr << "Serious error: mutation in protein has wiped out certain atoms .cnn files specification.  Please fix your .cnn files." << endl;
	exit(2);
      }

      int new_ptn_rep = int_list[ptn_rep];
      if (new_ptn_rep == 0 or new_ptn_rep == -1) {
	cout << "Serious error: mutation in protein has wiped out certain atoms .cnn files specification.  Please fix your .cnn files." << endl;
	cerr << "Serious error: mutation in protein has wiped out certain atoms .cnn files specification.  Please fix your .cnn files." << endl;
	exit(2);
      }
      this->atom_n_map[rot_rep] = new_ptn_rep;

    }


}
