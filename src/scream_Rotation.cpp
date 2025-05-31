#include "scream_Rotation.hpp"
#include <algorithm>

Scream_Rotation::Scream_Rotation() {

}

Scream_Rotation::Scream_Rotation(SCREAM_ATOM* start_atom, SCREAM_ATOM* end_atom, ScreamAtomV& atom_list) {

  this->start_atom = start_atom;
  this->end_atom = end_atom;
  
  cout << " 4.1" << endl;

  this->atoms_to_be_rotated = this->_return_atoms_to_be_rotated(this->start_atom, this->end_atom, atom_list);

}

Scream_Rotation::~Scream_Rotation() {

}

ScreamAtomV Scream_Rotation::_return_atoms_to_be_rotated(SCREAM_ATOM* start_atom, SCREAM_ATOM* end_atom, ScreamAtomV& atom_list) {

  // Currently, since atoms are connected by pointers to connected atoms, atom_list is not strictly needed to be passed in.

  // Error case.
  // if either connectivity_m.size == 0.

  // Initialization.
  ScreamAtomV to_be_rotated_atom_list; to_be_rotated_atom_list.clear();
  ScreamAtomV visited_atoms; visited_atoms.clear();
  //  visited_atoms.push_back(start_atom);  // Should be commented because if start_atom gets visited ever, there is a cycle, and wouldn't want to miss atoms connected to start_atom. 
  visited_atoms.push_back(end_atom); // Should this be commented?  Doesn't matter too much.  Unlike above, all connected atoms of this atom are guaranteed to have been visited even in the event of having a cycle/ring.  If commented out, means this atom will be revisited when there's a ring.  If not commented out, this atom wouldn't be revisited.  But either way, this atom will be removed from the final to_be_rotated_atom_list so it doesn't matter.  


  // Depth first search.

  for (map<SCREAM_ATOM*, int>::const_iterator itr = end_atom->connectivity_m.begin();
       itr != end_atom->connectivity_m.end(); 
       itr++) {

    SCREAM_ATOM* new_start_atom = end_atom;
    SCREAM_ATOM* connected_atom = itr->first;

    cout << "4.2" << endl;

    if (connected_atom != start_atom) { // Don't want to traverse back to original start atom; only for the bond defining rotation axis.
      ScreamAtomV append_atoms = this->_return_atoms_to_be_rotated_DFS(connected_atom, atom_list, visited_atoms);
      to_be_rotated_atom_list.insert(to_be_rotated_atom_list.begin(), append_atoms.begin(), append_atoms.end());
    }

  }

  cout << "4.8" << endl;

  // Below: may not be necessary.
  ScreamAtomVItr startAtomItr = std::find(to_be_rotated_atom_list.begin(), to_be_rotated_atom_list.end(), start_atom);
  if (startAtomItr != to_be_rotated_atom_list.end() ) {
    to_be_rotated_atom_list.erase(startAtomItr);
  }
  ScreamAtomVItr endAtomItr = std::find(to_be_rotated_atom_list.begin(), to_be_rotated_atom_list.end(), end_atom);
  if (endAtomItr != to_be_rotated_atom_list.end() ) {
    to_be_rotated_atom_list.erase(endAtomItr);
  }

  return to_be_rotated_atom_list;

}


ScreamAtomV Scream_Rotation::_return_atoms_to_be_rotated_DFS(SCREAM_ATOM* v, ScreamAtomV& atom_list, ScreamAtomV& visited_atoms) {

  visited_atoms.push_back(v);  
  ScreamAtomV to_be_rotated_atom_list; to_be_rotated_atom_list.clear();

  cout << "4.4" << endl;

  for (map<SCREAM_ATOM*, int>::const_iterator itr = v->connectivity_m.begin();
       itr != v->connectivity_m.end(); 
       itr++) {

    SCREAM_ATOM* connected_atom = itr->first;

    cout << " 4.5 " << endl;

    // Note: if there is a ring that goes back to original start atom, what ends up happening is that the entire molecule will get rotated.
    if ( find(visited_atoms.begin(), visited_atoms.end(), connected_atom ) == visited_atoms.end() ) { // Not visited
      ScreamAtomV append_atoms = this->_return_atoms_to_be_rotated_DFS(connected_atom, atom_list, visited_atoms);
      to_be_rotated_atom_list.insert(to_be_rotated_atom_list.end(), append_atoms.begin(), append_atoms.end());

    }
    else {
      continue;
    }
  }

  
  cout << " 4.6 " << endl;

  to_be_rotated_atom_list.push_back(v);
  return to_be_rotated_atom_list;

}
