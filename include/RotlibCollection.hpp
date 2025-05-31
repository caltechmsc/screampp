/* RotlibCollection.hpp
 * 
 * Header file for class that contains a collection of Rotamer Libraries.
 *
 */

#ifndef ROTLIBCOLLECTION_HPP
#define ROTLIBCOLLECTION_HPP

#include <string>
#include <map>
#include <set>
//using namespace std;

#include "Rotamer.hpp"
#include "AARotamer.hpp"
#include "Rotlib.hpp"

#include "ClashCollection.hpp"

typedef std::map<std::string, unsigned short> ExcitationEnumeration; ///< Stands for a Particular Excitation enumeration, like 13112.  char* stands for the name of the residue, like C218_X, whereas int stands for the level of the excitation.
//typedef std::map<std::string, AARotamer*> ExcitedRotamers; ///< A dictionary of Excitated Rotamers.
typedef std::map<std::string, Rotamer*> ExcitedRotamers;

/* Two new map structs, string --> n, a mapping for string--so that comparison would be much faster. */
typedef std::map<unsigned short, unsigned short> ExcitationEnumeration_n; ///< same as ExcitationEnumeration, except key is the mapped MutInfo.
typedef std::map<unsigned short, Rotamer*> ExcitedRotamers_n; ///< same as ExcitedRotamers, except key is the mapped MutInfo.

class RotlibCollection {

public:
  /* Constructors, Destructors and related functions. */
  RotlibCollection();		///< Constructs an empty RotlibCollection object.
  ~RotlibCollection();		///< RotlibCollection destructor.
  //  void addRotlib(const std::string, NtrlAARotlib*); ///< Adds a Rotamer Library to map<char*, NtrlAARotlib*> mutInfo_Rotlib_Map variable.
  void addRotlib(const std::string, Rotlib*);

  void addClashCollection(ClashCollection*); ///< Should be in constructor for cleanliness purposes; but rotlibcollection does not necessarily need this if not doing this particular Excitation energy-eliminate clash algorithm.

  /* These functions can all be accessed by swigged Python classes, with return types as dictionaries. */
  //  void initEmptyLatticeExcitationRotamers(); ///< Initializes a getNextEmptyLatticeExcitationRotamers() calls.


  /* initialization functions. use only one of the following in any run! */

  void initEmptyLatticeDataStructures(); ///< Initializes EmptyLattice related Date structures.

  void initDynamicMemoryDataStructures(); ///< Initializes data structure relevant to Dynamic Memory scheme.  These schemes use overlapping data structures.  
  void initAllocationUnderEnergyThreshold(double); ///< Initializes data structure relevant to using the memory allocation under a certain energy threshold for the sum of singles energies.  double: is the threshold.


  /* interface functions. */
  ExcitedRotamers getNextRotamersByELEnergy(); ///< One stop get rotamer routine.



  /* Other functions. */

  void resetEmptyLatticeCrntPstn(); ///< Resets the position of pointer to current Empty Lattice excitation.
  void resetTotalEnergyCrntPstn(); ///< Resets the position of pointer to current Total Energy position.
  ExcitedRotamers getNextEmptyLatticeExcitationRotamers(); ///< Returns the next EL excitation.
  ExcitedRotamers getNextTotalEnergyExcitationRotamers(); ///< Returns the next Excitation.

  ExcitedRotamers getNthEmptyLatticeExcitationRotamers(); ///< Returns the nth excitation.
  ExcitedRotamers getELExcitedRotamerFromEnumeration(ExcitationEnumeration_n&); ///< Returns a ExcitedRotamers structure from a ExcitationEnumeration structure.
  ExcitedRotamers_n getELExcitedRotamer_nFromEnumeration_n(ExcitationEnumeration_n&); ///< As name of fcn suggests.
  ExcitationEnumeration getELEnumerationFromExcitedRotamer(ExcitedRotamers&); // ER2EE: Acronym for ExcitationRotamer to ExcitationEnumeration.
  ExcitationEnumeration_n getELEnumeration_nFromExcitedRotamer_n(ExcitedRotamers_n&); // counterpart for _n for above.

  /* Converting to and from _n rep to string rep */
  ExcitationEnumeration_n _ExcitationEnumerationToExcitationEnumeration_n(ExcitationEnumeration&); // conversion from excitation enumeration to _n rep.
  ExcitedRotamers _ExcitedRotamers_nToExcitedRotamers(ExcitedRotamers_n&); // conversion from _n to strings.
  

  ClashCollection* getClashCollection(); ///< return the ClashCollection pointer.
  void cleanClashCollection(); ///< Clean pointer to ClashCollection.

  map<std::string, NtrlAARotlib*> getMutInfoRotlibMap(); ///< Gets the mutInfo_Rotlib_Map.
  map<std::string, NtrlAARotlib*> getMutInfoRotlibDict(); ///< Gets the mutInfo_Rotlib_Map. Called "dict" to be consistent in swigged functions.

  /* Functions involved in Dynamic memory scheme */

  ExcitedRotamers getNextDynamicMemoryRotamers_And_Expand(); ///< Returns the next EL excitation and expands collection.
  // A couple of private functions defined below in private section.
  ExcitedRotamers getNextDynamicClashEliminatedRotamers_And_Expand(); ///< Returns the next EL excitation, not including ones that fail the clash test, and expands collection.

  /* Functions involved in Allocate memory lower than to Energy threshold scheme. */
  void increaseConfigurationsUnderEnergyThreshold(double); ///< double: the amount of energy threshold to expand under.
  ExcitedRotamers getNextUnderEnergyThresholdRotamers(); ///< Gets the next rotamer while using the threshold rotamer scheme.
  

  /* Energy Query Functions. */
  void setExcitationEnergy(ExcitationEnumeration, double); ///< Sets ExcitationEnumeration's Energy to energy.
  double getExcitationEnergy(ExcitationEnumeration  EE) ; ///< Returns the energy of a particular ExcitationEnumeration.
  void printExcitationEnergyTable() const; ///< Prints a table of excitation levels of residues and the energies associated with them.  Output: stdout.
  void printEmptyLatticeLinearEnergyTable() const; ///< Prints a table of Empty Lattice Linear Sum energies (energies used to order the excitations).
  void printExcitationEnergyTable(std::string) const;
  void printEmptyLatticeTable() const; ///< Prints a table of Empty Lattice orderings, including energies.

  /* Misc get/set functions */
  std::string getInitMethod() const { return this->_initMethod;}; ///< Returns the initialization method, dynamicMemory or allocationUnderThresholdEnergy.
  int sizeOfSystem() { return RotamerLibraryMap.size(); }; ///< Returns the number of pstns that undergo sidechain replacement.
  double getHighestAllowedRotamerE() {return HIGHEST_ALLOWED_ROTAMER_E;}; ///< Get Highest Allowed Rotamer Energy to be considered.  Currently, this is used in initEmptyLatticeExcitationRotamers because a pre-ordering is being done.  When Dynamic programming for this process is implemented, HIGHEST_ALLOWED_ROTAMER_E would not be quite as useful, though will still improve speed.
  void setHighestAllowedRotamerE(double E) {HIGHEST_ALLOWED_ROTAMER_E = E;}; ///< Set Highest_Allowed_Rotamer_E.

  /* Public integer value */
  long double maxRotamerConfigurations; // maximum number of rotamer configurations that passed the clash test.
  int cmpMaxRotamerConfigurations(int); // Returns true if maxRotamerConfigurations < int.

private:
  /* Initialization method: dynamicMemory on the fly expansion, or allocationUnderThresholdEnergy. */
  string _initMethod; // Possible values: "DynamicMemory" and "AllocationUnderThreshold".
  

  /* where rotlibs are stored. */
  map<const std::string, unsigned short> _mutInfo_mapping; ///< Lookup table for a MutInfo string to an internal int for fast add/find/cmp etc.  
  map<unsigned short, std::string> _mutInfo_inverse_mapping; ///< Inverse mapping for above lookup table.
  //  map<const std::string, Rotlib*> RotamerLibraryMap;
  map<int, Rotlib*> RotamerLibraryMap;

  /* where empty lattice ordering is stored */
  set < ExcitationEnumeration_n > EnumerationList;  ///< EmptyLattices states.  Not ordered according to linear sum energy.
  list < ExcitationEnumeration_n > latest_round_EnumerationList;    ///< Track-keeping device to save time on calculations.  Represents the latest round, or newest round, of configurations that represent the very frontier of our configuration space.

  map < ExcitationEnumeration_n, double> ELEnumerationToEnergy; //< Linear combination of the energies of the spectrum energies on residues.
  map < ExcitationEnumeration_n, double> latest_round_ELEnumerationToEnergy; ///< Time saving structure.

  multimap< double, ExcitationEnumeration_n > ELEnergyToEnumeration;

  /* where calculated energies (i.e. energies including interactions) are stored */
  map< ExcitationEnumeration_n, double > TEEnumerationToEnergy; ///< The Energy corresponding to an Excitation State.
  multimap< double, ExcitationEnumeration_n > TEEnergyToEnumeration; ///< The Excitation state corredsponding to an energy.

  /* State Variables */
  multimap< double, ExcitationEnumeration_n >::iterator ELEnergyTablePstn; // Keeps track of where I'm currently at.
  int ELCurrentPnst_int; // Keeps track of where I'm currently at by counting.  Starts from 0.
  multimap< double, ExcitationEnumeration_n >::iterator TEEnergyTablePstn; // Keeps track of where I'm currently at.  TE stands for Total Energy.
  int TECurrentPstn_int; // Keeps track of where I'm currently at by counting.  Starts from 0.


  /* ClashCollection for eliminating configuration sets that have a pair of clashing rotamers.*/
  ClashCollection* clashCollection; 

  /* Helper functions */
  void _sortAllRotlibByEmptyLatticeEnergy(); ///< Sorts rotamers in all Rotlib's by Empty Lattice Energy.
  void _calcEmptyLatticeLinearSumEnergy(); ///< iterates through EnumerationList and calculates the linear sum of empty lattice energies and stores them in ELEnergyToEnumeration.
  void _calcLatestRoundEmptyLatticeLinearSumEnergy(); ///< like above, but a time saving devise, doesn't need to step through all items in EnumerationList.
  void _copyELEnumerationToEnergyToEnergyMultimap(); ///< Reverses key and value for ELEnumerationToEnergy and stores them in ELEnergyToEnumeration.  Now, key is Energy and value is ExcitedRotamers.
  void _copyEnergyMultimapToELEnumerationToEnergy(); ///< Reverses key and value for ELEnergyToEnumeration and stroes them in ELEnumerationToEnergy.
  void _copyLatestRoundELEnumerationToEnergyToEnergyMultimap(); ///< Like above, but only the latest round.
  void _copyEnergyMultimapToELEnumerationList(); ///< Copies 
  long double _calcMaxNoClashRotamerConfigurations(); ///< Calculates the total number of rotamer configurations that passed the clash test.

  bool _handCheckEEIdentity(const ExcitationEnumeration_n&, const ExcitationEnumeration_n&); ///< Checks whether two EE's are identical by explicitly comparison the contents of each multimap.

  /* Dynamic memory involved help functions. */

  ExcitedRotamers_n _getNextRotamerDynamicMemory(); //  these functions should do what they should do.
  void _expandDynamicMemory(ExcitationEnumeration_n&); // 

  //NtrlAARotlib* getNtrlAARotlib(const std::string mutInfo) { return RotamerLibraryMap.find(mutInfo)->second;}; ///< returns pointer to NtrlAARotlib by specifying mutInfo.
  Rotlib* getNtrlAARotlib(const std::string mutInfo); ///< returns pointer to Rotlib by specifying mutInfo.
  Rotlib* getNtrlAARotlib(int); ///< returns pointer to Rotlib by specifying the index for mutInfo (from lookup table).

  /* Helper functions for underThreshold implementation. */
  void _buildAllUnderThresholdConfigs(multimap< double, ExcitationEnumeration_n > &, 
				      double); ///< double: available E.  first int: index of Rotlib.  second int: index of excitation of that Rotlib meant for expansion.  Second double: energy of the passed in rotamer specified by the two ints.  The two multimaps: first one, stores the configurations from the previous round/lexiconigraphic entries.  Second one: all entries, including new ones that are being created.

  void _buildOneSetUnderThresholdConfigs(multimap< double, ExcitationEnumeration_n > &, 
					 multimap< double, ExcitationEnumeration_n > &, 
					 double, int, int, double); ///< double: available E.  first int: index of Rotlib.  second int: index of excitation of that Rotlib meant for expansion.  Second double: energy of the passed in rotamer specified by the two ints.  The two multimaps: first one, stores the configurations from the previous round/lexiconigraphic entries.  Second one: all entries, including new ones that are being created.

  /* Other State Variables */

  double HIGHEST_ALLOWED_ROTAMER_E;		///< Cutoff in energy.  Don't use any rotamers with a reference state energy > cutoff.
  double currentThresholdEnergy; ///< Cutoff in total energy for Allocation under energy threshold scheme.

  int currentEmptyLatticeExcitationN; ///< current Empty Lattice Excitation number, i.e. specifies the position.
  double getEmptyLatticeLinearSumEnergy( ExcitationEnumeration EE ) ;    	///< a conveniece function that takes in a Excitation Enumeration and returns the linear sum of Empty lattice energy.

  /* Zeroth Excitation Reference Variables */
  
  double bestEnergy;		///< bestEnergy So Far.
  int nSinceLastDescent;	///< number of excitations since last decrease in energy.
  int currentExcitationN;	///< current excitation number.

  /* Updated Reference State Variables, for iteration instead of storing all configurations in memory and sort them that way */
  /* Remark: of course, this saves disk space, but obviously slower */
  int crntExcitationN_itr_scheme; ///< Keeps tracks of current excitation number for use in iteration scheme.
  double crntExcitationN_itr_scheme_E; ///< Current energy.

  // multimap< double, ExcitationEnumeration > Energy2ExcitationEnumeration_itr_scheme; ///< The Excitation state corredsponding to energy.
//   // remark: the two members above should coincide exactly with currentExcitationN and TEEnergyToEnumeration.  There only there to illustrate the parallelism.

//   multimap< double, ExcitationEnumeration > temp_iteration_configurations; ///< Stores all configurations up to crnexcitationN_itr_sheme.
  multimap< double, ExcitationEnumeration_n > Energy2ExcitationEnumeration_itr_scheme; ///< The Excitation state corredsponding to energy.
//   // remark: the two members above should coincide exactly with currentExcitationN and TEEnergyToEnumeration.  There only there to illustrate the parallelism.

  multimap< double, ExcitationEnumeration_n > temp_iteration_configurations; ///< Stores all configurations up to crnexcitationN_itr_sheme.

  

};

#endif /* ROTLIBCOLLECTION_HPP */
