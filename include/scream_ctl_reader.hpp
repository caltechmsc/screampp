/*! \file scream_ctl_reader.hpp
  Header file containing routines for reading a SCREAM ctl file.
*/

#ifndef SCREAM_CTL_READER_HPP
#define SCREAM_CTL_READER_HPP

#include "MutInfo.hpp"
#include <string>
#include <vector>
using namespace std;

class ScreamParameters {
public:
  ScreamParameters() {};
  ScreamParameters(string);	///< Constructor.
  ~ScreamParameters() {};

  string InputFileName;		// no longer necessary
  vector<string> MutateResidueInfo; // specify the residues
  vector<string> AdditionalLibraryInfo; // Additional library info.  Stores connectivity/anchor info filenames.
  std::string Library;		// Name of the library.
  std::string PlacementMethod; // Method for placing a standard amino acid rotamer onto the backbone.
  std::vector<double> CreateCBParameters; // Parameters for how to create a CB atom on the backbone.  Order: offBisectorAngle, offPlaceAngle, bondLength, rotamerMatchVector.
  std::string KeepOriginalRotamer;	// Self Explanatory.  would make this a bool, but python doesn't have bools so make into string, easier to handle.
  std::string UseScreamEnergyFunction;  // Uses energy functions implemented in SCREAM, not ModSim.
  std::string UseDeltaMethod; // Specified which delta method to use: flat delta value, residue wide delta value, or full empirical delta values.
  std::string UseRotamerNeighborList; // Whether or not to use rotamer neighbor list.
  std::string UseAsymmetricDelta; // YES/NO.  whether or not to use Asymmetric delta.
  std::string UseDeltaForInterResiE; // YES/NO.  whether or not to use delta scheme on inter residue energy calculations
  double FlatDeltaValue;     // Delta value for all atoms on all rotamers to be placed.
  double DeltaStandardDevs;  // For Full Empricial delta values, number of standard devs away from mean error to use.
  double InnerWallScalingFactor; // Inner wall scaling factor for the "scaled vdw function".

  string NonPolarHCalc; // Whether non-polar H's are used for energy calculations.  YES: included.  NO: not included, also turns eletrostatics off.

  string ScoringFunction;	// Obsolete?  Maybe not yet.
  string MultiplePlacementMethod; // either Brute or EnergyExcitation
  string CBGroundSpectrumCalc;  // Whether the CB atom on canonical sidechains are included for the calculation of ground state spectrum.

  string OneEnergyFFParFile;  // Forcefield Parameter Par file for calculating One Energy.
  string DeltaParFile;        ///< Delta parameter file for using delta atom parameters.
  string EachAtomDeltaFile;   ///< Delta parameter for each atom in the bgf file.  Format: 2 columns, first column corresponds to atom number in bgf file.  Second column correspond to surface area.
  string PolarOptimizationExclusions;  ///< File that includes a list of all residues excluded in the polar optimization step.  Mainly, these guys are the surface exposed residues.  

  string LJOption; ///< 12-6, 11-6, 10-6, 9-6, 8-6 or 7-6.
  string CoulombMode; ///< Either normal coulomb or distance dependent coulomb.
  double Dielectric; ///< Dielectric value; if distance dependent, still the dielectric.

  
  int Selections;		///< Number of structures to return.
  int MaxSearchNumber; ///< Maximum number of rot conf. searches to be done.
  double AbsStericClashCutoffEL; ///< Absolute steric clash cutoff; if absolute energy > AbsStericClashCutoffEL in Empty Lattice Evaluation, Exit.
  double StericClashCutoffEnergy; ///< Relative steric clash cutoff.  Relative to ground state rotamer.
  double StericClashCutoffDist; ///< Temporarily disabled. (8-28-05)

  double MaxFinalStepRunTime;  ///< Maximum runtime on final combinatorial/optimization step.
  string LibPath;

  int Verbosity; ///< How much output is printed.

  //! Parameters relevant to DESIGN.
  vector<std::string> DesignPositionAndClass; ///< Defines positions at which you perform protein design on; uses multiple libraries for one site and of course includes mutations.

  vector<std::string> DesignAAClassDefns; ///< Defines the Classes defined in DesignPositionAndClass.  For more information, see Scream_par_file.doc.

  std::string JustOutputSequence; ///< If YES, outputs only the sequence in a file called Sequences.out.
  int StructuresPerSequence; ///< Defines the number of structures to output per unique sequence.


  //! Parameters relevant to BINDING_SITE_MODE.
  std::string BindingSiteMode; ///< AroundAtom, AroundResidue, AroundChain
  vector<std::string> FixedResidues; ///< Ntrl AA residues that would be fixed.  If self not included, self will be SCREAM'ed.
  vector<int> AroundAtom; ///< Around which atoms?
  vector<MutInfo> AroundResidue; ///< Around which residues?
  vector<std::string> AroundChain; ///< Around which chains?
  double AroundDistance; ///< Residues within this distance would be included in Binding Site SCREAM.
  std::string AroundDistanceDefn; ///< Distance to sidechains, or distance to backbone?  Possible values: SideChainOnly, BackBoneOnly and WholeResidue. 

  /* get, set functions for general use.  Especially with SWIG in mind. */

  vector<std::string> getMutateResidueInfoList() {return MutateResidueInfo; };
  vector<std::string> getAdditionalLibraryInfo() {return AdditionalLibraryInfo; };
  std::string getKeepOriginalRotamer();
  std::string getUseScreamEnergyFunction();

  std::string getPlacementMethod() { return PlacementMethod; };
  vector<double> getCreateCBParameters() { return CreateCBParameters; } ;

  std::string getUseDeltaMethod() { return UseDeltaMethod; };
  std::string getUseRotamerNeighborList() { return UseRotamerNeighborList;};
  std::string getUseAsymmetricDelta() { return UseAsymmetricDelta; };
  std::string getUseDeltaForInterResiE() { return UseDeltaForInterResiE; };
  double getFlatDeltaValue() {return FlatDeltaValue; };
  double getDeltaStandardDevs() {return DeltaStandardDevs; };
  double getInnerWallScalingFactor() { return InnerWallScalingFactor;};
  std::string getNonPolarHCalc() { return NonPolarHCalc;};
  std::string getOneEnergyFFParFile() {return OneEnergyFFParFile; };
  std::string getDeltaParFile() {return DeltaParFile; } ;
  std::string getEachAtomDeltaFile() { return EachAtomDeltaFile; };
  std::string getPolarOptimizationExclusions() { return PolarOptimizationExclusions; };
  int getSelections() {return Selections; };
  int getMaxSearchNumber() { return MaxSearchNumber; } ;
  double getAbsStericClashCutoffEL() { return AbsStericClashCutoffEL; } ; 
  double getStericClashCutoffEnergy() {return StericClashCutoffEnergy; };
  double getStericClashCutoffDist() { return StericClashCutoffDist;};
  
  double getMaxFinalStepRunTime() { return MaxFinalStepRunTime; };

  vector<std::string> getDesignPositionAndClass() { return DesignPositionAndClass; };
  vector<std::string> getDesignAAClassDefns() { return DesignAAClassDefns; } ;
  std::string getJustOutputSequence() { return JustOutputSequence; };

  std::string getLJOption() { return LJOption; };
  std::string getCoulombMode() { return CoulombMode; };
  double getDielectric() { return Dielectric; };

  std::string getBindingSiteMode() { return BindingSiteMode; } ;
  vector<std::string> getFixedResidues() { return FixedResidues; };
  vector<int> getAroundAtom() { return AroundAtom; };
  vector<MutInfo> getAroundResidue() { return AroundResidue; };
  vector<std::string> getAroundChain() { return AroundChain; };
  double getAroundDistance() { return AroundDistance; };
  std::string getAroundDistanceDefn() { return AroundDistanceDefn; } ;

  vector<std::string> getDesignPositions(); // returns a list of design positions.

  std::string getDesignClassFromPosition(std::string); ///< string: DesignPosition string.
  vector<std::string> getDesignClassAAs(std::string); ///< string: class name.

  string multiplePlacementMethod() { return MultiplePlacementMethod; };
  string getCBGroundSpectrumCalc() { return CBGroundSpectrumCalc; } ;
  
  void read_scream_par_file(string);
  void print_to_output(ostream*) const;

  string minimizationMethod() const;
  int minimizationSteps() const;
  string oneEMethod() const;

  string residueToScreamName(string) const;
  int residueToScreamPstn(string) const;
  string residueToScreamChn(string) const;

  string determineLibDirPath() const;  ///< returns library path to load rotamer library.  2 parameter consulted: Library and LibPath.
  string determineLibDirFileNameSuffix() const; ///< returns the suffix to the file name of a library, like _10.lib, etc.
  string determineCnnDirPath() const;  ///< returns consistent cnn path for rotamer library to be loaded.
  int getLibResolution() const; ///< Returns the resolution of library specified in .par file.

  string returnEnergyMethod() const; ///< Returns a string that contains all information what method is being used.  Final string output looks like FULL_ASYM_NOCB.
  double returnEnergyMethodTValue() const; ///< Returns a double that contains the value for the method being used.

  int getVerbosity() const {return this->Verbosity;}; ///< Returns the verbosity level.

  void _init_default_params();


};


#endif /* SCREAM_CTL_READER_HPP */
