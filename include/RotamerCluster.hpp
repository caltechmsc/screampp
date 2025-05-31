/* RotamerCluster.hpp class: 7-30-05 */

#ifndef ROTAMER_CLUSTER_HPP
#define ROTAMER_CLUSTER_HPP

#include "defs.hpp"

#include "Rotamer.hpp"
#include "RotConnInfo.hpp"

#include <vector>

using namespace std;

/* Justification for RotamerCluster being a subclass of Rotamer because:
1. Makes coding easier.  I don't have to rewrite everything for RotlibCollection.  This is the primary reason: while not conceptually clean, it's practical and that counts a lot.
2. RotamerCluster IS a "generalized/specialized" rotamer structure.  Therefore, many functions provided by the Rotamer class comes free.
3. Otherwise, RotamerCluster should really be a separate class.  Well, it's impressive how long this code has survived without restructuring anyway.  

Cons:
1. Many functions, like get_sc() etc becomes quite useless.  In fact, violates principle of OOP.
2. Unclean.
3. Hopefully, memory allocation wouldn't be too much of an issue.
*/

/*
  Principle: Rotamer Cluster does not actually INITIALIZE ANY atoms.  All accesses to those actual atoms are done through references of Rotamer*.
 */

class RotamerCluster : public Rotamer {
 public:
  RotamerCluster();
  RotamerCluster(Rotamer*, Rotamer*); ///< Convenient constructor for two rotamers/rotamerclusters.
  ~RotamerCluster();
  
  /* RotamerCluster specific functions */
  void addRotamerCluster(Rotamer*); ///< Function that adds a rotamer cluster into the rotamer child tree.
  vector<Rotamer*> getAllRotamers(); ///< Returns list of Rotamers, depth-first search sequence.


  /* Below: overridden functions */
  void print_Me() const; ///< Prints atom in bgf style,unordered.
  virtual void append_to_filehandle(ostream*) const {};
  virtual void pdb_append_to_filehandle(ostream*) const {}; 
  virtual void append_to_ostream_connect_info(ostream*) const {};

  ScreamAtomV get_sc_atoms();
  ScreamAtomV get_bb_atoms();

  //int number_of_atoms() const; 
  //double total_charge() const;
  

 private:
  RotConnInfo* rCI;
  //Rotamer* Left, *Right; ///< Binary tree representation of Rotamer Cluster. Recursively search for its contents in its left and right tree.
  vector<Rotamer*> childRotamerList; ///< Or, a multi tree representation.  Cleaner, easier to implement.
  string Idx; ///< a string to identify self.

};

#endif /* ROTAMER_CLUSTER_HPP */
