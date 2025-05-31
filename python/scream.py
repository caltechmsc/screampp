#!/usr/bin/env python

import sys, os, random, math, re
from packages import timing
#import timing
from packages.py_scream_ee import *

def usage():
    print '''
Does SCREAM.
Usage: %s <SCREAM ctl file>
E.g.: %s PheRS_F258_F260.ctl
''' % (sys.argv[0], sys.argv[0])
    sys.exit(0)

def main():
    if len(sys.argv) < 2:
        usage()

    # Control file
    ctl_file = sys.argv[1]

    # Required SCREAM variables
    SCREAM_NEW = os.getenv("SCREAM_NEW")
    if SCREAM_NEW == None:
        print ' Error!  enviromental variable SCREAM_NEW not set '
        sys.exit(2)
    print 'SCREAM_NEW     = ' + SCREAM_NEW

    SCREAM_NEW_LIB = os.getenv("SCREAM_NEW_LIB")
    if SCREAM_NEW_LIB == None:
        print ' Error!  enviromental variable SCREAM_NEW_LIB not set '
        sys.exit(2)
    print 'SCREAM_NEW_LIB = ' + SCREAM_NEW_LIB

    SCREAM_NEW_CNN = os.getenv("SCREAM_NEW_CNN")
    if SCREAM_NEW_CNN == None:
        print ' Error!  enviromental variable SCREAM_NEW_CNN not set '
        sys.exit(2)
    print 'SCREAM_NEW_CNN = ' + SCREAM_NEW_CNN

    SCREAM_NEW_RTF = os.getenv("SCREAM_NEW_RTF")
    if SCREAM_NEW_RTF == None:
        print ' Error!  enviromental variable SCREAM_NEW_RTF not set '
        sys.exit(2)
    print 'SCREAM_NEW_RTF = ' + SCREAM_NEW_RTF

    SCREAM_NEW_CHG = os.getenv("SCREAM_NEW_CHG")
    if SCREAM_NEW_CHG == None:
        print ' Error!  enviromental variable SCREAM_NEW_CHG not set '
        sys.exit(2)
    print 'SCREAM_NEW_CHG = ' + SCREAM_NEW_CHG

    # SCREAM_NEW_SCALE_COU = os.getenv("SCREAM_NEW_SCALE_COU")
    # if SCREAM_NEW_SCALE_COU == None:
    #     print ' Error!  environmental variable SCREAM_NEW_SCALE_COU not set '
    #     sys.exit(2)
    # print 'SCREAM_NEW_SCALE_COU = ' + SCREAM_NEW_SCALE_COU

    # SCREAM_NEW_SCALE_HB = os.getenv("SCREAM_NEW_SCALE_HB")
    # if SCREAM_NEW_SCALE_HB == None:
    #     print ' Error!  environmental variable SCREAM_NEW_SCALE_HB not set '
    #     sys.exit(2)
    # print 'SCREAM_NEW_SCALE_HB = ' + SCREAM_NEW_SCALE_HB

    # SCREAM_NEW_SCALE_VDW = os.getenv("SCREAM_NEW_SCALE_VDW")
    # if SCREAM_NEW_SCALE_VDW == None:
    #     print ' Error!  environmental variable SCREAM_NEW_SCALE_VDW not set '
    #     sys.exit(2)
    # print 'SCREAM_NEW_SCALE_VDW = ' + SCREAM_NEW_SCALE_VDW

    print

    # AssignProteinCharges
    f = open(ctl_file, 'r')
    ctl_text = f.read()
    m = re.search('InputFileName\s+(\S+)',ctl_text)
    apcbgf = m.group(1)
    assign_protein_charges = '/project/Biogroup/scripts/perl/AssignProteinCharges.pl -b ' + apcbgf + ' --' + SCREAM_NEW_CHG + ' -o'
    os.system(assign_protein_charges)

    # Initialization.
    timing.start()
    SCREAM_MODEL = ScreamModel(ctl_file)
    (SCREAM_PARAMS, BGF_HANDLER, ptn, scream_EE) = (SCREAM_MODEL.scream_parameters, SCREAM_MODEL.HANDLER, SCREAM_MODEL.ptn, SCREAM_MODEL.scream_EE)
    TIMING = open('timing.txt', 'w')
    
    # Rotamer Library init, and relevant info.
    (Libraries_Dict, mutInfo_rotConnInfo_Dict, ntrl_orig_rots) = initRotlibs(SCREAM_PARAMS, ptn)  # need ntrl_orig_rots so doesn't get deleted

    ntrlMutInfo_list = SCREAM_PARAMS.getMutateResidueInfoList()
    additionalLib_list = SCREAM_PARAMS.getAdditionalLibraryInfo()

    # scream_EE init.
    init_Scream_EE(scream_EE, SCREAM_PARAMS, ptn, Libraries_Dict, mutInfo_rotConnInfo_Dict)
    timing.finish()
    print >>TIMING, 'SCREAM setup took ' + str(timing.micro() / 1000000.00) + ' seconds.'
    TIMING.flush()
    
    # Stage 1. Empty Lattice Energy calculations.
    timing.start()
    Calc_EL_Energies(ptn, scream_EE, SCREAM_PARAMS, Libraries_Dict, mutInfo_rotConnInfo_Dict, BGF_HANDLER) # temp added BGF_HANDLER
    timing.finish()
    #Print_Library_EL_Energies(Libraries_Dict)
    print >>TIMING, 'SCREAM singles took ' + str(timing.micro() / 1000000.00) + ' seconds.'
    TIMING.flush()

    # print best structure after this stage.

    EL_RC = RotlibCollectionPy()
    initialize_rotlibCollection(EL_RC, Libraries_Dict)

    #printBestELStructure(ptn, BGF_HANDLER, EL_RC, mutInfo_rotConnInfo_Dict, 'BestEL.bgf')

    # Stage 2. Clash Resolution.  I.e. Clustering clashing rotamers.  Iteration until no longer detect clustering rotamers; i.e. final Libraries_Dict
    print 'Start Clash Resolution (or DeterministicClashEliminationWithSnapshotField)'
    
    tmp_Libraries_Dict = Libraries_Dict.copy()
    tmp_mutInfo_rotConnInfo_Dict = mutInfo_rotConnInfo_Dict.copy()

    clash_Libraries_Dict = {}

    if SCREAM_PARAMS.multiplePlacementMethod() == 'ExcitationWithClustering':
        print 'Doing ExcitationWithClustering!'
        (clash_scream_EE, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict) = ground_state_calc_and_cluster(SCREAM_MODEL, scream_EE, tmp_Libraries_Dict, tmp_mutInfo_rotConnInfo_Dict, TIMING) # retains scream_EE
        #output_file = 'Residue-E-NoClash.txt'
        #BGF_HANDLER.printToFile('NoClash.bgf')
        #Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file)

        #sys.exit(2)
        #rotlibCollection = excitation_after_clash_resolution(SCREAM_MODEL, scream_EE, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict, TIMING)
        rotlibCollection = excitation_after_clash_resolution(SCREAM_MODEL, clash_scream_EE, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict, TIMING)
        printBestStructures(SCREAM_MODEL, BGF_HANDLER, scream_EE, rotlibCollection, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict, Libraries_Dict, mutInfo_rotConnInfo_Dict)
        print 'Done ExcitationWithClustering!'
        sys.exit(2)
        
    elif SCREAM_PARAMS.multiplePlacementMethod() == 'ClusteringThenDoubletExcitation':
        print 'Doing ClusteringThenDoubleExcitation!'
        (clash_scream_EE, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict) = ground_state_calc_and_cluster(SCREAM_MODEL, scream_EE, tmp_Libraries_Dict, tmp_mutInfo_rotConnInfo_Dict, TIMING) # retains scream_EE
        #BGF_HANDLER.printToFile('NoClash.bgf')
        #output_file = 'Residue-E-NoClash.txt'
        #Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file)
        #sys.exit(2)
        doublet_excitation(SCREAM_MODEL, scream_EE, Libraries_Dict, mutInfo_rotConnInfo_Dict, TIMING)

        #BGF_HANDLER.printToFile('Doublets-Optimize.bgf')
        #output_file = 'Residue-E-Doublets.txt'
        #Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file)
    elif SCREAM_PARAMS.multiplePlacementMethod() == 'ClusteringThenDoubletsThenSinglets':
        print 'Doing ClusteringThenDoubletsThenSinglets!'
        (clash_scream_EE, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict) = ground_state_calc_and_cluster(SCREAM_MODEL, scream_EE, tmp_Libraries_Dict, tmp_mutInfo_rotConnInfo_Dict, TIMING) # retains scream_EE
        #BGF_HANDLER.printToFile('NoClash.bgf')
        #output_file = 'Residue-E-NoClash.txt'
        #Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file)
        #sys.exit(2)
        doublet_excitation(SCREAM_MODEL, scream_EE, Libraries_Dict, mutInfo_rotConnInfo_Dict, TIMING)

        #BGF_HANDLER.printToFile('Doublets-Optimize.bgf')
        #output_file = 'Residue-E-Doublets.txt'
        #Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file)

        T_list = [0.001]
        deterministicClashEliminationWithSnapshotField(SCREAM_MODEL, scream_EE, tmp_Libraries_Dict, tmp_mutInfo_rotConnInfo_Dict, T_list, TIMING, BGF_HANDLER)
        #BGF_HANDLER.printToFile('Singlets-Optimize.bgf')
        #output_file = 'Residue-E-Singlets.txt'
        #Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file)
        

        #printBestStructures(SCREAM_MODEL, BGF_HANDLER, scream_EE, rotlibCollection, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict, Libraries_Dict, mutInfo_rotConnInfo_Dict)
        
    elif SCREAM_PARAMS.multiplePlacementMethod() == 'DeterministicClashEliminationWithProgressiveField':
        pass
    
    elif SCREAM_PARAMS.multiplePlacementMethod() == 'DeterministicClashEliminationWithSnapshotField':
        T_list = [600, 300, 100, 0.1]
        #deterministicClashEliminationWithSnapshotField(SCREAM_MODEL, scream_EE, tmp_Libraries_Dict, tmp_mutInfo_rotConnInfo_Dict, T_list, TIMING, BGF_HANDLER)
        #BGF_HANDLER.printToFile('Stochastic.bgf')
        print 'Done Deterministicclasheliminationwithsnapshotfield!'
        #sys.exit(2)
        
    elif SCREAM_PARAMS.multiplePlacementMethod() == 'StochasticClashEliminationWithSnapshotField':
        pass
    
    elif SCREAM_PARAMS.multiplePlacementMethod() == 'ClusteringAndThenAnneal':
        (clash_scream_EE, clash_Libraries_Dict, clash_mutInfo_rotConnInfo_Dict) = ground_state_calc_and_cluster(SCREAM_MODEL, scream_EE, tmp_Libraries_Dict, tmp_mutInfo_rotConnInfo_Dict, TIMING) # retains scream_EE
        #BGF_HANDLER.printToFile('NoClash.bgf')
        print 'Clustering done, now anneal.'
        T_list = [600, 300, 100, 0.1]
        deterministicClashEliminationWithSnapshotField(SCREAM_MODEL, scream_EE, tmp_Libraries_Dict, tmp_mutInfo_rotConnInfo_Dict, T_list, TIMING, BGF_HANDLER)
        #BGF_HANDLER.printToFile('Stochastic.bgf')
        print 'Annealing done!'

    # Brand new and final RotlibCollection, after all clashes have been eliminiated.
    timing.start()
    rotlibCollection = RotlibCollectionPy()
    rotlibCollection.setHighestAllowedRotamerE(SCREAM_PARAMS.StericClashCutoffEnergy)

    #Print_Library_EL_Energies(clash_Libraries_Dict,0)
    clash_Libraries_Dict = tmp_Libraries_Dict
    initialize_rotlibCollection(rotlibCollection, clash_Libraries_Dict)

    # Need to remember to add clashcollection!  Else crashes.  Bug.
    clashCollection = ClashCollection(15)
    rotlibCollection.addClashCollection(clashCollection)
    scream_EE.addClashCollection(clashCollection)

    print 'ClashCollection Added.'
    
    timing.finish()
    print >>TIMING, 'SCREAM setup for Opportunity stage took ' + str( timing.micro() / 1000000.00) + ' seconds.'
    TIMING.flush()

    # Print best structure after Clash Resolution.
    #printBestELStructure(ptn, BGF_HANDLER, rotlibCollection, mutInfo_rotConnInfo_Dict, 'NoClashEL.bgf')

    # Stage 3.  Hydrogen bond network optimization.
    # vcvicek
    # ResidueReachFile = '/project/Biogroup/Software/SCREAM/lib/SCREAM_delta_par_files/Residue_Reach.par'
    SCREAM_NEW = os.getenv("SCREAM_NEW")
    if SCREAM_NEW == None:
        print ' Error!  enviromental variable SCREAM_NEW not set '
        sys.exit(2)
    ResidueReachFile = SCREAM_NEW+'/lib/delta/Residue_Reach.par'
    
    ResidueExclusionFile = SCREAM_PARAMS.getPolarOptimizationExclusions()
    print 'Starting HBNetworkOptimization!'
    #HBNetworkOptimization(SCREAM_MODEL, Libraries_Dict, mutInfo_rotConnInfo_Dict, ResidueReachFile, ResidueExclusionFile, TIMING)
    print 'Done HBNetworkOptimization.'
    TIMING.close()

    # Last: Clean up: print energies.
    output_file = 'Residue-E.txt'
    BGF_HANDLER.printToFile('best_1.bgf')
    Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file)
    
    print 'SCREAM all done!  Exiting.'


    

def HBNetworkOptimization(SCREAM_MODEL, Libraries_Dict, mutInfo_rotConnInfo_Dict, ResidueReachFile, ResidueExclusionFile, TIMING):
    # Variable setup.
    BGF_HANDLER = SCREAM_MODEL.HANDLER
    SCREAM_PARAMS = SCREAM_MODEL.scream_parameters
    ptn = SCREAM_MODEL.ptn
    scream_EE = SCREAM_MODEL.scream_EE
    scream_EE.resetFlags(1)
    scream_EE.initScreamAtomVdwHbFields()
    
    HB_OUTPUT = open('HBOptimization.txt', 'w')

    # 1) Find polar pairs that are close to each other, and optimize them.  Setup list of close feasible polar contacts.  Use "Residue_Reach.par" table to determine this.  Go through list in increasing distance.
    (polar_res_pair_list, polar_Libraries_Dict) = polar_res_pair_priority(ptn, Libraries_Dict, ResidueReachFile, ResidueExclusionFile)

    print >>HB_OUTPUT, '-----------------------------------------'
    print >>HB_OUTPUT, 'Possibly coupled polar residues (Distance of CB atom, includes repeats): '
    print >>HB_OUTPUT, '%7s %7s %7s' % ('Dist', 'Res1', 'Res2')
    for i in polar_res_pair_list:
        print >>HB_OUTPUT, '%7.3f %7s %7s' % (i[0], i[1], i[2])
    print >>HB_OUTPUT, '-----------------------------------------'
    # 2) Optimize all pairs of H-Bonding residues.
    timing.start()
    print >>HB_OUTPUT, 'Hydrogen bond network optimization: Stage 1.  Polar Residue pair optimization.  These residues pairs are optimized in the presence of other sidechains.'
    solved_flag = [] # contains those residues that have already been touched
    for i in polar_res_pair_list:
        print >>HB_OUTPUT, '-----------------------------------------'
        (mI1_str, mI2_str) = (i[1], i[2])
        if (mI1_str in solved_flag) or (mI2_str in solved_flag):
            continue

        mI1mI2_Str = mI1_str + ' ' + mI2_str
        
        (mI1, mI2) = (MutInfo(mI1_str), MutInfo(mI2_str))
        print >>HB_OUTPUT, 'Optimizing polar residue pair ' , mI1_str, mI2_str
        HB_OUTPUT.flush()
        # Reset scream_EE fixed-moveable info.
        scream_EE.resetFlags(1)
        scream_EE.initScreamAtomVdwHbFields()

        #(all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)
        #print >>HB_OUTPUT, 'Energies before pair optimization: %7.3f %7.3f %7.3f %7.3f' % (all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E)
        
        scream_EE.fix_all()
        scream_EE.moveable_mutInfo(mI1)
        scream_EE.moveable_mutInfo(mI2, None, 1) # Last argument == 1 means also do setup_variableAtomsOnEachSidechain().
        
        crnt_polar_pair_dict = {}
        crnt_polar_pair_dict[mI1_str] = polar_Libraries_Dict[mI1_str]
        crnt_polar_pair_dict[mI2_str] = polar_Libraries_Dict[mI2_str]
        
        crnt_pair_rotlibCollection = RotlibCollectionPy()
        crnt_pair_rotlibCollection.setHighestAllowedRotamerE(SCREAM_PARAMS.StericClashCutoffEnergy)
        initialize_rotlibCollection(crnt_pair_rotlibCollection, crnt_polar_pair_dict)
        
        pair_rotamer = crnt_pair_rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()
        (c, max_c, lowestE) = (0, SCREAM_MODEL.scream_parameters.getMaxSearchNumber(), 999999999)
        max_c = 5000
        isLessThanMaxCount = crnt_pair_rotlibCollection.cmpMaxRotamerConfigurations(max_c)

        while len(pair_rotamer) != 0:
            c += 1
            if c > max_c:
                break
            print >>HB_OUTPUT, ''
            print >>HB_OUTPUT, 'Rotamer pair ', c, ' (Residue breakdown: ',
            for mI_str in pair_rotamer.keys():
                print >>HB_OUTPUT, mI_str, pair_rotamer[mI_str].get_rotamer_n(),
            print >>HB_OUTPUT, ')'
            placeMultipleConformers(ptn, pair_rotamer, mutInfo_rotConnInfo_Dict, scream_EE) # Adjusts setup_variableAtomsOnEachSidechain if there's a mutation.
            (all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)
            crntEnergy = all_interactions_E
            EL_E = 0
            for mI in pair_rotamer.keys():
                rot = pair_rotamer[mI]
                EL_E += rot.get_empty_lattice_E()
                #print mI, pair_rotamer[mI].get_rotamer_n(), ' EL energy:', rot.get_empty_lattice_E()

            crntEnergy += EL_E
            print >>HB_OUTPUT, 'Total:         ' + str(crntEnergy)
            print >>HB_OUTPUT, 'EmptyLattice:  ' + str(EL_E)
            print >>HB_OUTPUT, 'All Else VDW:  ' + str(all_vdw_E)
            print >>HB_OUTPUT, 'All Else HB:   ' + str(all_hb_E)
            print >>HB_OUTPUT, 'All Else Coul: ' + str(all_coulomb_E)
            
            crnt_pair_rotlibCollection.setEnergyForExcitedRotamers(pair_rotamer, crntEnergy) # set
            pair_rotamer = crnt_pair_rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()

        crnt_pair_rotlibCollection.resetTotalEnergyCrntPstn()
        (best_E, best_pair_rotamers) = placeBestStructureSoFar(ptn, crnt_pair_rotlibCollection, mutInfo_rotConnInfo_Dict)
        print >>HB_OUTPUT,  'Best pair has local energy: ', best_E, ' ' , mI1_str, mI2_str, ' (self + doublet + interaction with rest of protein)'
        HB_OUTPUT.flush()
        solved_flag.append(mI1_str)
        solved_flag.append(mI2_str)

        scream_EE.resetFlags(1)
        scream_EE.initScreamAtomVdwHbFields()
        
        (total_E, vdw_E, hb_E, coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)
        EL_E = 0
        for mI_tmp in Libraries_Dict.keys():
            mutInfo = MutInfo(mI_tmp)
            EL_E += ptn.getEmptyLatticeEnergy(mutInfo.getChn(), mutInfo.getPstn())
        total_E += EL_E

        print >>HB_OUTPUT, 'Overall Protein Energy after optimizing %s pair: (TotalE EL_E vdwE hbE coulombE):  %7.3f %7.3f %7.3f %7.3f %7.3f' % (mI1mI2_Str, total_E, EL_E, vdw_E, hb_E, coulomb_E)
        print >>HB_OUTPUT, '-----------------------------------------'
        HB_OUTPUT.flush()
    
    timing.finish()
    print >>TIMING, 'SCREAM Hydrogen bond Stage 1 (pairs) network optimization took ' + str( timing.micro() / 1000000.00) + ' seconds.'
    TIMING.flush()
    
    BGF_HANDLER.printToFile('PolarOptimized-Stage1.bgf')
    # 3) Second Stage of Polar Optimization: singles polar.  basically a iterations step.
    timing.start()
    print >>HB_OUTPUT, 'Hydrogen bond network optimization: Stage 2.  Optimizing polar residues, singles!'



    for mI in polar_Libraries_Dict.keys():
        (c, total_E) = (0,0)

        if mI in solved_flag:
            continue
        else:
            scream_EE.fix_all()
            scream_EE.moveable_mutInfo(MutInfo(mI), None, 1)
            solved_flag.append(mI)
            print >>HB_OUTPUT, '-----------------------------------------'
            print >>HB_OUTPUT, 'Optimizing single residue: ', mI
            c+=1
            crntLib = polar_Libraries_Dict[mI]
            crntLib.reset_pstn()
            crntRotamer = crntLib.get_next_rot()
            (best_E, best_rot) = (999999999,None)

            while crntRotamer != None:
                print >>HB_OUTPUT, 'Rotamer %d' % c
                c+=1
                if crntRotamer.get_empty_lattice_E() > 500:
                    print >>HB_OUTPUT, ' EL E > 500, not sampling.'
                    crntRotamer = crntLib.get_next_rot()
                    continue
                placeSingleConformer(ptn, mI, crntRotamer, mutInfo_rotConnInfo_Dict, scream_EE)
                (interaction_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)

                EL_E = crntRotamer.get_empty_lattice_E_abs()
                total_E = interaction_E + EL_E
                
                print >>HB_OUTPUT, 'Total:         ' + str(total_E)
                if SCREAM_PARAMS.getVerbosity() == 1:
                    print >>HB_OUTPUT, 'EmptyLattice:  ' + str(EL_E)
                    print >>HB_OUTPUT, 'All Else VDW:  ' + str(all_vdw_E)
                    print >>HB_OUTPUT, 'All Else HB:   ' + str(all_hb_E)
                    print >>HB_OUTPUT, 'All Else Coul: ' + str(all_coulomb_E)

                if total_E < best_E:
                    (best_E, best_rot) = (total_E, crntRotamer)
                crntRotamer = crntLib.get_next_rot()

        placeSingleConformer(ptn, mI, best_rot, mutInfo_rotConnInfo_Dict, scream_EE)
        scream_EE.resetFlags(1)
        scream_EE.initScreamAtomVdwHbFields()

        (total_E, vdw_E, hb_E, coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)

        EL_E = 0
        for mI_tmp in Libraries_Dict.keys():
            mutInfo = MutInfo(mI_tmp)
            EL_E += ptn.getEmptyLatticeEnergy(mutInfo.getChn(), mutInfo.getPstn())
        total_E += EL_E
        
        print >>HB_OUTPUT, 'Overall Protein Energy after optimizing %s (TotalE EL_E vdwE hbE coulombE):  %7.3f %7.3f %7.3f %7.3f %7.3f' % (mI, total_E, EL_E, vdw_E, hb_E, coulomb_E)
        print >>HB_OUTPUT, '-----------------------------------------'

    timing.finish()
    print >>TIMING, 'SCREAM Hydrogen bond Stage 2 (singles) network optimization took ' + str ( timing.micro() / 1000000.00) + ' seconds.'
    TIMING.flush()
    

    BGF_HANDLER.printToFile('PolarOptimized-Stage2.bgf')
    BGF_HANDLER.printToFile('best_1.bgf')
    #printBestELStructure(ptn, BGF_HANDLER, rotlibCollection, mutInfo_rotConnInfo_Dict, 'AfterPolarOptimization.bgf')


    HB_OUTPUT.close()

def Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, output_file):
    ptn = SCREAM_MODEL.ptn
    SCREAM_PARAMS = SCREAM_MODEL.scream_parameters

    OUTPUT = open(output_file, 'w')
    #print >>OUTPUT, '%7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s' % ('Res', 'Intern', 'V_sBB', 'V_oBB', 'V_fSC', 'C_sBB', 'C_oBB', 'C_fSC', 'H_sBB', 'H_oBB', 'H_fSC', 'V_Tot_fm', 'C_Tot_fm', 'H_Tot_fm', 'Tot_fm', 'V_ScSc', 'C_ScSc', 'H_ScSc', 'Tot_ScSc', 'Total' )
    print >>OUTPUT, '%7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s' % ('Res', 'Intern', 'Tot_EL', 'V_EL', 'C_EL', 'H_EL', 'Tot_ScSc', 'V_ScSc', 'C_ScSc', 'H_ScSc', 'Total' )

    clashCollection = ClashCollection(15)
    scream_EE.addClashCollection(clashCollection)
    
    scream_EE.resetFlags(1)
    scream_EE.initScreamAtomVdwHbFields()
    print 'before'
    (total_E, vdw_E, hb_E, coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)

    #print 'after Calc_This_Interaction_Energy in Print_Final_Residue_Energy_Breakdown'
    
    ntrlAA_list = SCREAM_PARAMS.getMutateResidueInfoList()
    if len(ntrlAA_list) != 0:
        if ntrlAA_list[0] == 'BINDING_SITE':
            ntrlAA_list = return_ntrlAA_list_from_BINDING_SITE_mode(SCREAM_PARAMS, ptn)
        # nothing to do for now for DESIGN mode.
    (Intern,  Tot_EL,    V_EL,    C_EL,    H_EL) = (0,0,0,0,0)
    mI_E_Dict = {}
    for mI in ntrlAA_list:
        E_list = MutInfoEnergies(ptn, scream_EE, SCREAM_PARAMS, mI)
        mI_E_Dict[mI] = E_list
        Intern += E_list[0]
        Tot_EL += E_list[1]
        V_EL += E_list[2]
        C_EL += E_list[3]
        H_EL += E_list[4]

    Total_E = Intern + Tot_EL + total_E
    print >>OUTPUT, '%7s' % 'All', '%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f' % (Intern, Tot_EL, V_EL, C_EL, H_EL, total_E, vdw_E, coulomb_E, hb_E, Total_E)    

    for mI in ntrlAA_list:
        print >>OUTPUT, '%7s' % mI, '%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f' % mI_E_Dict[mI]


    OUTPUT.close()
    pass


def make_new_rotlib(SCREAM_MODEL, Libraries_Dict, mutInfo_rotConnInfo_Dict = None):
    # Makes a new rotamer library, usually combining 2 rotamer libraries.
    print '<<<<<<<<<<<<<<<<<<<<<<<<<<<New rotlib>>>>>>>>>>>>>>>>>>>>>>>>>'
    
    # First, Setups.
    new_Rotlib = SCREAM_MODEL.new_Rotlib()
    ptn = SCREAM_MODEL.ptn
    SCREAM_PARAMS = SCREAM_MODEL.scream_parameters
    rotlibCollection = RotlibCollectionPy()
    scream_EE = Scream_EE() # not adding to SCREAM_MODEL because it's only being used in this routine.
    
    for mutInfoString in Libraries_Dict.keys():
        mI = MutInfo(mutInfoString)
        for arblib_mI in mutInfo_rotConnInfo_Dict:
            mI.searchAndAddRotConnInfo(MutInfo(arblib_mI), mutInfo_rotConnInfo_Dict[arblib_mI])
        scream_EE.addMutInfoRotConnInfo(mI)
    #scream_EE.init_after_addedMutInfoRotConnInfo_neighbor_list(ptn, SCREAM_PARAMS.getOneEnergyFFParFile(), SCREAM_PARAMS.getDeltaParFile() ) # why neighbor_list? 5-14-07.  ahh... neighbor_lists actually are never used, even though they are set up.
    if SCREAM_PARAMS.getCoulombMode() == 'Normal':
        scream_EE.setNormalDielectric(SCREAM_PARAMS.getDielectric())
    elif SCREAM_PARAMS.getCoulombMode() == 'DistanceDependent':
        scream_EE.setDistanceDielectricPrefactor(SCREAM_PARAMS.getDielectric())
        
    scream_EE.init_after_addedMutInfoRotConnInfo_on_the_fly_E(ptn, SCREAM_PARAMS) # is this necessary?
    _initScreamAtomDeltaValuePy(scream_EE, SCREAM_PARAMS) # added 9-23-06

    # First do excitations.
    initialize_rotlibCollection(rotlibCollection, Libraries_Dict)
    theseExcitedConformers = rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()
    (total_count, count_since_last_best) = (0,0)
    max_count = SCREAM_PARAMS.getMaxSearchNumber()
    lowestE = 999999999

    #maxRotamerConfigurations = rotlibCollection.maxRotamerConfigurations
    max_count = 500
    isLessThanMaxCount = rotlibCollection.cmpMaxRotamerConfigurations(max_count)
    #print "Is the maximum number of rotamer configurations less than 1000?  If so, just evaluate all rotamers in positions in this cluster (1==YES, 0==NO): " + str(isLessThanMaxCount)
    #sys.exit(2)

    new_mutInfo = MutInfo() # 11-29-05

    while len(theseExcitedConformers) != 0:
        #print total_count 
        total_count += 1
        if total_count > max_count:
            break
        count_since_last_best = count_since_last_best +1
        placeMultipleConformers(ptn, theseExcitedConformers, mutInfo_rotConnInfo_Dict, scream_EE)
        (all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)
        crntEnergy = all_interactions_E
    
        for mutInfo in theseExcitedConformers.keys():
            crntConformer = theseExcitedConformers[mutInfo]
            crntEnergy += crntConformer.get_empty_lattice_E()

        # Now add new rotamer to new rotlib.
        new_cluster_rotamer = new_Rotlib.new_rotamer_cluster()
        new_mutInfo = MutInfo()
        for mutInfo in theseExcitedConformers.keys():
            crntConformer = theseExcitedConformers[mutInfo]
            new_cluster_rotamer.addRotamerCluster(crntConformer)
            new_mutInfo.addMutInfo(mutInfo)
            
        
        new_cluster_rotamer.set_empty_lattice_E_abs(crntEnergy)
        new_cluster_rotamer.set_rotamer_n(total_count)

        # Keep going.
        rotlibCollection.setEnergyForExcitedRotamers(theseExcitedConformers, crntEnergy)
                
        if crntEnergy < lowestE:
            lowestE = crntEnergy
            count_since_last_best = 0


        if rotlibCollection._shouldKeepGoing(count_since_last_best, 0, 1) or isLessThanMaxCount == 1:
            theseExcitedConformers = rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()
        else:
            break

    # cleaning up new_Rotlib.
    new_Rotlib.sort_by_empty_lattice_E()
    new_Rotlib.reset_pstn()

    new_Library_Dict = {}
    new_Library_Dict[new_mutInfo.getString()] = new_Rotlib

    t = new_Library_Dict

    return new_Library_Dict



def clash_resolution(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict):
    # Might as well implement clash resolution for arbitrary rotamer libraries as well.
    # Current
    pass



def deterministicClashEliminationWithSnapshotField(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, T_list, TIMING, BGF_HANDLER):

    #T_list = [0.1, 500, 400, 300, 200, 100, 0.1]
    #T_list = [600, 300, 100, 0.01]
    #T_list = [0.001]
    OUTPUT = open('Anneal-Energies.txt', 'w')
    OUTPUT.close()
    for c in range(1,len(T_list)+1):
        timing.start()
        print 'Starting round ', c, ' of deterministicClashEliminationWithSnapshotField procedure.'
        _oneRoundDeterministicClashEliminationWithSnapshotField(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, T_list[c-1])
        #print 'Energies after this round, round %d, is %7.3f: ' % (c, energy)
        print 'End round ', c, ' of deterministicClashEliminationWithSnapshotField procedure.'
        timing.finish()
        print >>TIMING, 'SCREAM deterministicClashEliminationWithSnapshotField round', c, 'took ' + str(timing.micro() / 1000000.00) + ' seconds.'
        TIMING.flush()
        BGF_HANDLER.printToFile('Field%d.bgf' % c)

def _oneRoundDeterministicClashEliminationWithSnapshotField(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, T):
    # 1. Protein will be continuously updated.
    # 2. The background (snapshot) of the sidechains are visible throughout.  Because of this, no new energy expression needs to be set up.  Yay!
    # 3. Don't eliminate 2 clashing sidechains at a time, just 1 clashing sidechain at a time.  Because of this, no need to use clustered rotamer libraries.  Yay!  That was a serius pain in the ass...
    # 4. I.E. the algorithm is essentially iterative.  Biased Monte Carlo, i.e. stochastic.  "Biased" comes from clash energy calculations.
    
    # Step 1.  Find clash list.  Do a 1 point energy calculating and return a list that contains all the clashes.
    clashCollection = ClashCollection(15)

    scream_EE.addClashCollection(clashCollection)

    scream_EE.resetFlags(1)
    scream_EE.initScreamAtomVdwHbFields()



    print 'Calculating all interaction energies before this round begins:'
    (all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_MODEL.scream_parameters)
    
    list = clashCollection.getDiscreteClashPairList()

    if len(list) != 0:
        print 'Clashes identified.  Listing them:'
        for mutInfoPair in list:
            mI1 = mutInfoPair.getMutInfo1()
            mI2 = mutInfoPair.getMutInfo2()
            print mI1.getString(), mI2.getString()
    else:
        print 'No Clashes identified in this round.'


    # Step 2.  Determine surface residues.
    surfaceResidueList = get_polar_optimization_exclusions(SCREAM_MODEL)
    if len(surfaceResidueList) != 0:
        print 'Polar Optimization/Surface Residues definitions given.  They are:'
        for mI in surfaceResidueList:
            print mI
    else:
        print 'No Polar Optimization/Surface Residues definitions given.  Proceeding.'

    OUTPUT = open('Anneal-Energies.txt', 'a')
    # Step 3. Remove clashes that involve the polar residues.
    Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict
    solved_mi_list = []
    for mutInforPair in list:
        (mI1, mI2) = (mutInforPair.getMutInfo1().getString(), mutInforPair.getMutInfo2().getString())
        if mI1 in surfaceResidueList:
            SCList = [mI1]
            print 'Solving Surface Clashing Sidechain: ', mI1
            _placeSideChainWithFullField_and_print(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, T, OUTPUT)
            solved_mi_list.append(mI1)
        if mI2 in surfaceResidueList:
            SCList = [mI2] 
            print 'Solving Surface Clashing Sidechain: ', mI2
            _placeSideChainWithFullField_and_print(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, T, OUTPUT)
            solved_mi_list.append(mI2)

    # Step 4.  Remove clashes that involve everything else.
    for mutInfoPair in list:
        (mI1, mI2) = (mutInforPair.getMutInfo1().getString(), mutInforPair.getMutInfo2().getString())
        if not (mI1 in solved_mi_list):
            SCList = [mI1]
            print 'Solving core sidechain: ', mI1
            _placeSideChainWithFullField_and_print(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, T, OUTPUT)
            solved_mi_list.append(mI1)
        if not (mI2 in solved_mi_list):
            print 'Solving core sidechain: ', mI2
            SCList = [mI2] 
            _placeSideChainWithFullField_and_print(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, T, OUTPUT)
            solved_mi_list.append(mI2)

    # Step 5.  SCREAM everything else.
    rL_list = Primieval_Libraries_Dict.keys()
    #rL_list = randomizeList(rL_list)
    #for mI in Primieval_Libraries_Dict.keys():
    rL_E_pair = []#, rL_new =[]
    for mI in rL_list:
        if mI in solved_mi_list:
            continue
        elif mI in surfaceResidueList:
            continue
        else:
            (PreCalc_E, EL_all, EL_rot_vdw_E, EL_rot_coulomb_E,  EL_rot_hb_E, All_E, All_rot_vdw_E, All_rot_coulomb_E, All_rot_hb_E, Total_E) = MutInfoEnergies(SCREAM_MODEL.ptn, scream_EE, SCREAM_MODEL.scream_parameters, mI)
            rL_E_pair.append([Total_E, mI])

    rL_E_pair.sort(lambda x, y: cmp(x[0], y[0]) )

    #     for mI in rL_list:
    #         if mI in solved_mi_list:
    #             continue
    #         elif mI in surfaceResidueList:
    #             continue
    #         else:
    #             print 'Solving rest of protein.'
    #             SCList = [mI]
    #             _placeSideChainWithFullField_and_print(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, T, OUTPUT)
    #             solved_mi_list.append(mI)
    print 'Solving rest of protein.'
    for rL_E in rL_E_pair:
        SCList = [rL_E[1]]
        _placeSideChainWithFullField_and_print(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, T, OUTPUT)

    # Step 6. cleanup.
    scream_EE.cleanClashCollection()

    OUTPUT.close()


def _placeSideChainWithFullField_and_print(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, temperature, OUTPUT):
    (total_E, EL_E, vdw_E, hb_E, coulomb_E) = _placeSideChainWithFullField(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, temperature)
    mI = SCList[0]
    print >>OUTPUT, 'After annealing %s (totE, EL, vdw, hb, coulomb): %7.3f %7.3f %7.3f %7.3f %7.3f' % (mI, total_E, EL_E, vdw_E, hb_E, coulomb_E)

def randomizeList(l):
    c = len(l)
    for i in range(0,c):
        rand = random.randint(0,c-i-1) # inclusive
        (l[(c-1)-i], l[rand]) = (l[rand], l[(c-1)-i])
    return l


    

def _placeSideChainWithFullField(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, SCList, temperature):
    # Calculates the sidechain energy interacting with other originally variable sidechains (interaction energy with backbone already calculated).
    # Note: places only 1 Sidechain at a time.
    ptn = SCREAM_MODEL.ptn
    SCREAM_PARAMS = SCREAM_MODEL.scream_parameters

    (method, t) = _returnEnergyMethod(SCREAM_PARAMS)
    scream_EE.resetFlags(1)
    scream_EE.initScreamAtomVdwHbFields()
    rot_vdw_E = scream_EE.calc_all_interaction_vdw_E_delta(method, t)
    rot_hb_E = scream_EE.calc_all_interaction_hb_E_delta(method, t)
    rot_coulomb_E = scream_EE.calc_all_interaction_coulomb_E_delta()

    #rot_vdw_hb_exclusion_E = scream_EE.calc_empty_lattice_vdw_hb_exclusion_E_delta(mI, method, t)

    EL_E = 0
    for mI_tmp in Primieval_Libraries_Dict.keys():
        mutInfo = MutInfo(mI_tmp)
        EL_E += ptn.getEmptyLatticeEnergy(mutInfo.getChn(), mutInfo.getPstn())

    rot_E = rot_vdw_E + rot_hb_E + rot_coulomb_E + EL_E
    #print 'Checking energies before trying out difference sidechains, should match previous lowest E values: VDW: ', rot_vdw_E, ' HB: ', rot_hb_E, ' rot_coulomb_E: ', rot_coulomb_E, ' EL: ', EL_E, 'Total (including precalc) ', rot_E
    EL_E = 0

    #print ''
    #print 'Placing sidechain (in order of Empty Lattice energy, only below a certain cutoff threshold): ', 
    #for mI in SCList:
    #print mI, ' '
    for mI in SCList:
        energy_map = {}
        library = Primieval_Libraries_Dict[mI]
        library.sort_by_empty_lattice_E()
        library.reset_pstn()
        currentRotamer = library.get_next_rot()
        (mutAA, mutPstn, mutChn) = unpackMutInfo(mI)
        
        scream_EE.fix_all()
        scream_EE.moveable_mutInfo(MutInfo(mI), None, 1)
        
        while currentRotamer != None:
            
            if currentRotamer.get_empty_lattice_E_abs() > SCREAM_MODEL.scream_parameters.getAbsStericClashCutoffEL():
                currentRotamer = library.get_next_rot()
                continue
            if mutAA != 'Z':
                #ptn.ntrlRotamerPlacement(mutChn, mutPstn, currentRotamer)
                scream_EE.ntrlRotamerPlacement(mutChn, mutPstn, currentRotamer)
                if ptn.mutationDone():
                    int_map = ptn.getNewMapping()
                    for i in mutInfo_rotConnInfo_Dict.keys():
                        mutInfo_rotConnInfo_Dict[i].modifyMappingInProteinAtoms(int_map)
                if currentRotamer.get_is_Original_flag() == 1:
                    pass
                else:
                    rot_preCalc_E = currentRotamer.get_preCalc_TotE() 

            else:
                rotConnInfo_Z = mutInfo_rotConnInfo_Dict[mutInfo]
                ptn.conformerPlacement(currentRotamer, rotConnInfo_Z)

            # Now decide which SCREAM energy function to use
            rot_EL_E = currentRotamer.get_empty_lattice_E_abs()
            rot_vdw_E = scream_EE.calc_all_interaction_vdw_E_delta(method, t)
            rot_hb_E = scream_EE.calc_all_interaction_hb_E_delta(method, t)
            rot_coulomb_E = scream_EE.calc_all_interaction_coulomb_E_delta()
            rot_E = rot_vdw_E + rot_hb_E + rot_coulomb_E + rot_EL_E


            if SCREAM_PARAMS.getVerbosity() == 1:
                print ''
                print 'Rotamer ', currentRotamer.get_rotamer_n()
                print 'Total E:         ' + str(rot_E)
                print 'PreCalc:         ' + str(rot_preCalc_E)
                print 'Empty Lattice E: ' + str(rot_EL_E) + ' (includes PreCalc)'
                print 'VDW:             ' + str(rot_vdw_E)
                print 'HB:              ' + str(rot_hb_E)
                print 'Coulomb:         ' + str(rot_coulomb_E)

            energy_map[currentRotamer] = rot_E
            currentRotamer = library.get_next_rot()

        # Now pick a sidechain to place on protein
        (k, T, R)  = (1.38*pow(10,-23)/4.18/1000, temperature, 6.023*pow(10,23))
        adjusted_beta = 1/(k*T*R)
        lowest_E = ''
        lowest_rot_n = ''
        partition_function = 0
        rotamerList = [] # list of rotamers
        prob_point = [] # probability entry point; for instance, line from 0 - 1 divided into 0.4, 0.457, 0.667, 0.83, 0.92, 1.0.  Then first rotamer has prob of 0.4 of getting picked.
        for rotamer in energy_map.keys():
            if lowest_E == '':
                lowest_E = energy_map[rotamer]
                lowest_rot_n = rotamer.get_rotamer_n()
            if energy_map[rotamer] < lowest_E:
                lowest_E = energy_map[rotamer]
                lowest_rot_n = rotamer.get_rotamer_n()
        for rotamer in energy_map.keys():
            E = energy_map[rotamer]
            boltzmann_ratio = math.exp(-(E-lowest_E) * adjusted_beta)
            partition_function += boltzmann_ratio
            if boltzmann_ratio > 0.99:
                pass
                #print 'Dominant rotamer: ', rotamer.get_rotamer_n()
            rotamerList.append(rotamer)
            prob_point.append(boltzmann_ratio)
        #print 'Lowest Energy value is: ', lowest_E
        #print 'Lowest Energy Rotamer number is: ', lowest_rot_n
        #print 'Boltzmann ratio: ',
        #for ratio in prob_point:
        #print ratio
        #print 'Partition function value: ', partition_function
        prob_point = map(lambda x: x/partition_function, prob_point)
        for i in range(1,len(prob_point)):
            prob_point[i] += prob_point[i-1]
        # print prob_point
        #print 'Prob point printing!'
        #for i in range(0,len(prob_point)):
        #print i, ' ', prob_point[i], ' ', energy_map[rotamerList[i]] - lowest_E
        #print 'Rotamer n s:'
        #for r in rotamerList:
        #print r.get_rotamer_n()
        X = random.uniform(0,1)
        rotamerToPlace = ''
        for i in range(0,len(prob_point)):
            if X < prob_point[i]:
                rotamerToPlace = rotamerList[i]
                break
        #print 'Rotamer number picked: ' , rotamerToPlace.get_rotamer_n()
        # Now place that rotamer
        if mutAA != 'Z':
            #ptn.ntrlRotamerPlacement(mutChn, mutPstn, rotamerToPlace)
            scream_EE.ntrlRotamerPlacement(mutChn, mutPstn, rotamerToPlace)
            if ptn.mutationDone():
                int_map = ptn.getNewMapping()
                for i in mutInfo_rotConnInfo_Dict.keys():
                    mutInfo_rotConnInfo_Dict[i].modifyMappingInProteinAtoms(int_map)
            if rotamerToPlace.get_is_Original_flag() == 1:
                pass
            else:
                #rot_preCalc_E = rotamerToPlace.get_preCalc_TotE()
                rot_EL_E = rotamerToPlace.get_empty_lattice_E_abs()

        else:
            rotConnInfo_Z = mutInfo_rotConnInfo_Dict[mutInfo]
            ptn.conformerPlacement(rotamerToPlace, rotConnInfo_Z)

        
        ptn.setEmptyLatticeEnergy(MutInfo(mI).getChn(), MutInfo(mI).getPstn(), rotamerToPlace.get_empty_lattice_E_abs() )
        ptn.setPreCalcEnergy(MutInfo(mI).getChn(), MutInfo(mI).getPstn(), rotamerToPlace.get_preCalc_TotE() )

        # Now calculate overall energy.
        scream_EE.resetFlags(1)
        vdw_E = scream_EE.calc_all_interaction_vdw_E_delta(method, t)
        hb_E = scream_EE.calc_all_interaction_hb_E_delta(method, t)
        coulomb_E = scream_EE.calc_all_interaction_coulomb_E_delta()
        EL_E = 0
        for mI_tmp in Primieval_Libraries_Dict.keys():
            mutInfo = MutInfo(mI_tmp)
            EL_E += ptn.getEmptyLatticeEnergy(mutInfo.getChn(), mutInfo.getPstn())
        
        total_E = EL_E + vdw_E + hb_E + coulomb_E
        #print "Overall energy after %s is tried (totE, EL, vdw, hb, coulomb): %7.3f %7.3f %7.3f %7.3f %7.3f" % (mI, total_E, EL_E, vdw_E, hb_E, coulomb_E)
        
        return (total_E, EL_E, vdw_E, hb_E, coulomb_E)
    
            
def _returnEnergyMethod(SCREAM_PARAMS):
    # Format: e.g.: FLAT_ASYM_NOCB_96  (96 == LJ 9-6.  if 126, == LJ 12-6.)

    method = SCREAM_PARAMS.getUseDeltaMethod()
    t = 0.0
    if SCREAM_PARAMS.getUseAsymmetricDelta() == 'YES':
        method += '_ASYM'
    if method[0:4] == 'FLAT':
        t = SCREAM_PARAMS.getFlatDeltaValue()
    elif method[0:4] == 'FULL':
        t = SCREAM_PARAMS.getDeltaStandardDevs()
    elif method[0:6] == 'SCALED':
        t = SCREAM_PARAMS.getInnerWallScalingFactor()
    if SCREAM_PARAMS.getCBGroundSpectrumCalc() == 'NO':
        method += '_NOCB'

    LJOption = SCREAM_PARAMS.getLJOption()
    method = method + '_' + LJOption

    return (method, t)


def excitation_after_clash_resolution(SCREAM_MODEL, scream_EE, clash_resolved_Libraries, clash_resolved_mutInfo_rotConnInfo_Dict, TIMING):
    # Excitation step after clashes has been resolved (or kind of resolved).

    ptn = SCREAM_MODEL.ptn
    SCREAM_PARAMS = SCREAM_MODEL.scream_parameters

    scream_EE.resetFlags(1)
    scream_EE.initScreamAtomVdwHbFields()
    #scream_EE.

    
    rotlibCollection = RotlibCollectionPy()
    rotlibCollection.setHighestAllowedRotamerE(SCREAM_PARAMS.getStericClashCutoffEnergy())
    # initialize rotlibcollection
    initialize_rotlibCollection(rotlibCollection, clash_resolved_Libraries)
    # add clashcollection, or crashes, kind of a bug
    clashCollection = ClashCollection(15)
    rotlibCollection.addClashCollection(clashCollection)
    scream_EE.addClashCollection(clashCollection)
    
    theseExcitedConformers = rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()
    total_count = 0
    count_since_last_best = 0
    max_count = SCREAM_PARAMS.getMaxSearchNumber()
    lowestE = 99999999

    timing.start()
    Combinatorial_Time = 0.0
    allowed_time = SCREAM_PARAMS.getMaxFinalStepRunTime()


    while len(theseExcitedConformers) != 0:
        total_count += 1
        if total_count > max_count:
            break
        count_since_last_best += 1
        placeMultipleConformers(ptn, theseExcitedConformers, clash_resolved_mutInfo_rotConnInfo_Dict, scream_EE)
        (E, vdw_E, hb_E, coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)
        total_E = E
        for mutInfo in theseExcitedConformers.keys():
            crntConformer = theseExcitedConformers[mutInfo]
            total_E += crntConformer.get_empty_lattice_E_abs()

        rotlibCollection.setEnergyForExcitedRotamers(theseExcitedConformers, total_E)

        if total_E < lowestE:
            lowestE = total_E
            count_since_last_best=0

        timing.finish()
        Combinatorial_Time = Combinatorial_Time + (timing.micro() / 1000000.00)
        timing.start()

        if rotlibCollection._shouldKeepGoing(count_since_last_best, Combinatorial_Time, allowed_time):
            theseExcitedConformers = rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()
        else:
            break

    timing.finish()
    Excitation_Time = timing.micro() / 1000000.00


    return rotlibCollection

    #printEnergyEvolution(rotlibCollection)

    #print "Empty Lattice Calculation time: " + str(EL_calculation_time)
    #print "Combinatorial Calculation time: " + str(Combinatorial_Time)
    #print "Total number of rotamer configuration sets evaluated: " + str(total_count)

def SCREAM_resolve_clashes(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, TIMING):
    # WORKING ON THIS!  May nee to rewrite this whole routine; for step #1 below, might need to move into main body when fully implemented.
    
    # SCREAM algorithm for ground state and cluster:
    # 1. For residues that only admit 1 rotamer, keep it as fixed for the rest of the calculation.  By "admit" I mean energies that are relative energies above a specific threshold, say 20 kcal/mol.
    # 2. Do same thing as ground_state_calc_and_cluster.
    # 3. If ground state cluster does not converge (i.e. still a clash after that stage), do double excitation on the clashing residues.
    timing.start()

    print 'Doing SCREAM_resolve_clashes!  An update of old SCREAM clash resolution.'

    ptn = SCREAM_MODEL.ptn
    SCREAM_PARAMS.scream_parameters

    new_scream_EE = scream_EE
    # 1. Decide if there are any rotamers that will be fixed for the rest of this calculation.
    thres_E = 20 # Threshold energy for not accepting rotamer.
    (new_Library_Dict, new_mutInfo_rotConnInfo_Dict) = Interesting_Ground_States(Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, thres_E)
    
    rotlibCollection = RotlibCollectionPy()
    rotlibCollection.setHighestAllowedRotamerE(SCREAM_PARAMETERS.getStericClashCutoffEnergy())

    initialize_rotlibCollection(rotlibCollection, new_Library_Dict())

    # ClashCollection Init.
    clashCollection = ClashCollection(15)
    new_scream_EE.addClashCollection(clashCollection)
    rotlibCollection.addClashCollection(clashCollection)

    # Reinitialize scream_EE: using info from new_Library_Dict (i.e. the on
    new_scream_EE.fix_all()
    for mI in new_mutInfo_rotConnInfo_Dict.keys():
        new_scream_EE.moveable_mutInfo(mI, new_mutInfo_rotConnInfo_Dict[mI])
    new_scream_EE.setup_variableAtomsOnEachSidechain() 
    new_scream_EE.initScreamAtomVdwHbFields()
    new_scream_EE._initScreamAtomDeltaValuePy()

    # Step 2: do clustering.  Very similar to ground_state_calc_and_cluster, but with a twist: i.e. if not converging, will do double excitation instead.
    while (1):
        theseExcitedConformers = rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()
        placeClusters(ptn, theseExcitedConformers, new_mutInfo_rotConnInfo_Dict, new_scream_EE)
        (all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(new_scream_EE, SCREAM_PARAMS)

        print 'Number of Clashes: ' + str(clashCollection.getNumberOfClashes())
        if clashCollection.getNumberOfClashes() == 0:
            print 'All "All ground state clashes eliminated!'
            break
        else:
            print 'Eliminating ground state clashes.'
            list = clashCollection.getDiscreteClashPairList() # returns top clashing pair
            # For each clash component, do local multiple excitation.
            for mutInfoPair in list:
                mI1 = mutInfoPair.getMutInfo1()
                mI2 = mutInfoPair.getMutInfo2()

                local_Library_Dict = {}
                local_Library_Dict[mI1.getString()] = new_Libraries_Dict[mI1.getString()]
                local_Library_Dict[mI2.getString()] = new_Libraries_Dict[mI2.getString()]

                # build new Rotamers for this set of RotlibCollection and scream_EE
                clustered_Library_Dict = make_new_rotlib(SCREAM_MODEL, local_Library_Dict, mutInfo_rotConnInfo_Dict)
                
                # update Libraries_Dict, mutInfo_rotConnInfo_Dict
                del new_Libraries_Dict[mI1.getString()]
                del new_Libraries_Dict[mI2.getString()]


def Interesting_Ground_States(Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, E):
    # Returns a new Library Dict that contains only those libraries that have multiple rotamers acceptable below a threshold E.
    new_Library_Dict = {}
    new_mutInfo_rotConnInfo_Dict = {}
    for mI in Primieval_Libraries_Dict.keys():
        Lib = Primieval_Libraries_Dict[mI]
        Lib.sort_by_empty_lattice_E()
        (c, rot) = (0,Lib.get_next_rot() )
        while rot != None:
            if rot.get_empty_lattice_E() > E:
                break
            (c, rot) = (c+1, Lib.get_next_rot() )
            if c >= 2:
                new_Library_Dict[mI] = Lib
                new_mutInfo_rotConnInfo_Dict[mI] = mutInfo_rotConnInfo_Dict(mI)
                break

    return new_Library_Dict



def doublet_excitation(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, TIMING):
    # Doublet excitation: hopefully, resolves some clashes that are unresolvable in the clash resolution stage.
    # Algorithm:
    # 1. Go through list of clashes; from high to low, do doublet excitation for highest energy pair.
    # 2. Repeast step 1, until energy no longer improves.
    timing.start()
    ptn = SCREAM_MODEL.ptn
    SCREAM_PARAMS = SCREAM_MODEL.scream_parameters

    # Reinitialize scream_EE for safety, in case earlier operations didn't initialize scream_EE right.

    #scream_EE.setup_variableAtomsOnEachSidechain() # NEW ADDITION
    #scream_EE.initScreamAtomVdwHbFields() # Needed for HSE/HIS "mutations".

    print '+++++++++++++++++++++++++++++++++++++++++++++++'
    print 'Doing SCREAM doublet excitation.'
    visited_clashing_mutInfoPairs = []
    while (1):
        # ClashCollection Init.
        clashCollection = ClashCollection(15)
        scream_EE.addClashCollection(clashCollection)
        #rotlibCollection.addClashCollection(clashCollection) # is this necessary?
        scream_EE.resetFlags(1)
        scream_EE.initScreamAtomVdwHbFields()
        (all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS)
        if clashCollection.getNumberOfClashes() == 0:
            print 'No further clashes needed to be resolved by doublet excitation.'
            break
        list = clashCollection.getDiscreteClashPairList() # Returns a list of clashing pairs, in descending clashing E.
        print 'Number of clashes: ', len(list), '  Printing them:'
        for l in list:
            print l.getMutInfo1().getString(), l.getMutInfo2().getString(), str(l.getClashE() )
        mutInfoPair = None
        for i in range(0, len(list)-1):
            l = list[i]
            if l in visited_clashing_mutInfoPairs:
                continue
            else:
                mutInfoPair = l
                break

        if not mutInfoPair:
            print 'All venue explored, doublets can\'t improve energies any further.'
            break

        (mI1, mI2) = (mutInfoPair.getMutInfo1(), mutInfoPair.getMutInfo2() )
        visited_clashing_mutInfoPairs.append(MutInfoPair(mI1, mI2, mutInfoPair.getClashE()))
        print 'Exploring doublet possibility for ', mI1.getString(), mI2.getString()
        
        # For each clash component, do local multiple excitation.

        doublet_Library_Dict = {}
        doublet_Library_Dict[mI1.getString()] = Primieval_Libraries_Dict[mI1.getString()]
        doublet_Library_Dict[mI2.getString()] = Primieval_Libraries_Dict[mI2.getString()]
        
        # build new Rotamers for this set of RotlibCollection and scream_EE
        rotlibCollection = RotlibCollectionPy()
        initialize_rotlibCollection(rotlibCollection, doublet_Library_Dict)
        
        doubletConformers = rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()

        scream_EE.resetFlags()
        scream_EE.fix_all()

        (mI1_rCI, mI2_rCI) = (None, None)
        try:
            (mI1_rCI, mI2_rCI) = (mutInfo_rotConnInfo_Dict[mI1.getString()], mutInfo_rotConnInfo_Dict[mI2.getString()])
        except KeyError:
            (mI1_rCI, mI2_rCI) = (None, None)
        scream_EE.moveable_mutInfo(mI1, mI1_rCI, 0)
        scream_EE.moveable_mutInfo(mI2, mI2_rCI, 1)

        # Later on: put in calculation of original energies.
        #         crnt_EL_E = ptn.getEmptyLatticeEnergy(mI1.getChn(), mI1.
        #         vdw_E = scream_EE.calc_all_interaction_vdw_E_delta(method, t)
        #         hb_E = scream_EE.calc_all_interaction_hb_E_delta(method, t)
        #         coulomb_E = scream_EE.calc_all_interaction_coulomb_E_delta()
        #         E = vdw_E + hb_E + coulomb_E + EL_E
        (crnt_E, best_doublet, c) = ('', '', 0)
        
        while len(doubletConformers) != 0 and c <= 500:
            placeMultipleConformers(ptn, doubletConformers, mutInfo_rotConnInfo_Dict, scream_EE)
            # Now calulate energies
            EL_E = 0
            for key in doubletConformers:
                EL_E += doubletConformers[key].get_empty_lattice_E_abs()
            (method, t) = _returnEnergyMethod(SCREAM_PARAMS)
            vdw_E = scream_EE.calc_all_interaction_vdw_E_delta(method, t)
            hb_E = scream_EE.calc_all_interaction_hb_E_delta(method, t)
            coulomb_E = scream_EE.calc_all_interaction_coulomb_E_delta()
            E = vdw_E + hb_E + coulomb_E + EL_E
            if E < crnt_E or crnt_E == '':
                best_doublet = doubletConformers
                crnt_E = E
            # Increment.
            (c, doubletConformers) = (c+1, rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy() )
        print 'New energy: ' , crnt_E
        placeMultipleConformers(ptn, best_doublet, mutInfo_rotConnInfo_Dict, scream_EE)
    
            
    timing.finish()
    print >>TIMING, 'SCREAM doublet excitation (after clustering) took ' + str(timing.micro() / 1000000.00) + ' seconds.'


def ground_state_calc_and_cluster(SCREAM_MODEL, scream_EE, Primieval_Libraries_Dict, mutInfo_rotConnInfo_Dict, TIMING):
    # For now, mutInfo_rotConnInfo_Dict says nothing new.  all Libraries_Dict for ground_state_calc and cluster MUST be NtrlAARotlib.
    # I.e., all rotConnInfo == NULL.
    # Here, apparently, neighbor_list is used to do clustering.

    # RotlibCollection Init.
    timing.start()
    
    ptn = SCREAM_MODEL.ptn
    SCREAM_PARAMS = SCREAM_MODEL.scream_parameters

    new_scream_EE = scream_EE
    new_Libraries_Dict = Primieval_Libraries_Dict.copy()
    new_mutInfo_rotConnInfo_Dict = mutInfo_rotConnInfo_Dict.copy()

    rotlibCollection = RotlibCollectionPy()
    rotlibCollection.setHighestAllowedRotamerE(SCREAM_PARAMS.StericClashCutoffEnergy)

    initialize_rotlibCollection(rotlibCollection, new_Libraries_Dict)

    # ClashCollection Init.
    clashCollection = ClashCollection(15)
    new_scream_EE.addClashCollection(clashCollection)
    rotlibCollection.addClashCollection(clashCollection)

    # Reinitialize scream_EE for safety, in case earlier operations didn't initialize scream_EE right.
    new_scream_EE.setup_variableAtomsOnEachSidechain() # NEW ADDITION
    new_scream_EE.initScreamAtomVdwHbFields() # NEW ADDITION

    print '------------------------------------'
    (last_clash_n, crnt_clash_n, rounds) = (0,0,0)
    # Keep going until the ClashCollection returns no clashes
    while (1):
        theseExcitedConformers = rotlibCollection.getNextDynamicMemoryRotamers_And_ExpandPy()
        placeClusters(ptn, theseExcitedConformers, mutInfo_rotConnInfo_Dict, new_scream_EE)
        (all_interactions_E, all_vdw_E, all_hb_E, all_coulomb_E) = Calc_This_Interaction_Energy(new_scream_EE, SCREAM_PARAMS)

        #for mI_tmp in Primieval_Libraries_Dict.keys():
        #mutInfo = MutInfo(mI_tmp)
        #print mutInfo.getString(), ' '
        #print ptn.getEmptyLatticeEnergy(mutInfo.getChn(), mutInfo.getPstn())
        
        # get clash collection info
        (last_clash_n, crnt_clash_n, rounds) = (crnt_clash_n, clashCollection.getNumberOfClashes(), rounds+1)
        print 'Number of Clashes identified: ' + str(crnt_clash_n)
        if crnt_clash_n == 0:
            print 'All ground state clashes eliminated!'
            break
        #elif crnt_clash_n == last_clash_n and rounds >= 20:
            #print 'Reached limit of clashing algorithm, quiting, perserving current best structure.'
        else:
            print "Eliminating ground state clashes."
            # cluster rotlibs that correspond to rotamer, calc new cluster energies.
            list = clashCollection.getDiscreteClashPairList() # returns top clashing pair
            # For each clash component, do local multiple excitation.
            for mutInfoPair in list:
                mI1 = mutInfoPair.getMutInfo1()
                mI2 = mutInfoPair.getMutInfo2()

                local_Library_Dict = {}
                local_Library_Dict[mI1.getString()] = new_Libraries_Dict[mI1.getString()]
                local_Library_Dict[mI2.getString()] = new_Libraries_Dict[mI2.getString()]

                # build new Rotamers for this set of RotlibCollection and scream_EE
                clustered_Library_Dict = make_new_rotlib(SCREAM_MODEL, local_Library_Dict, mutInfo_rotConnInfo_Dict)
                
                # update Libraries_Dict, mutInfo_rotConnInfo_Dict
                del new_Libraries_Dict[mI1.getString()]
                del new_Libraries_Dict[mI2.getString()]

                new_Libraries_Dict.update(clustered_Library_Dict)
                
            # new round of RotlibCollection and scream_EE
            print "End round of clash elimination.  Examining if another round needed." 
            rotlibCollection = RotlibCollectionPy()
            new_scream_EE = Scream_EE()
            initialize_rotlibCollection(rotlibCollection, new_Libraries_Dict)

            initialize_scream_EE(ptn, SCREAM_PARAMS, new_scream_EE, new_Libraries_Dict, mutInfo_rotConnInfo_Dict)
            
            # ClashCollection Init.
            clashCollection = ClashCollection(15)
            new_scream_EE.addClashCollection(clashCollection)
            rotlibCollection.addClashCollection(clashCollection)
            scream_EE.addClashCollection(clashCollection)

    print 'Clash resolution finished!'

    timing.finish()
    print >>TIMING, 'SCREAM clash resolution (clustering) took ' + str(timing.micro() / 1000000.00) + ' seconds.'
    TIMING.flush()
    
    return (new_scream_EE, new_Libraries_Dict, new_mutInfo_rotConnInfo_Dict)

def setup_polar_libraries_dict(Libraries_Dict):
    # Returns a dictionary of polar residues, in dictionary format: {mutInfo : RotLib}
    polar_Libraries_Dict = {}
    for mI in Libraries_Dict.keys():
        mI_string = mI
#ADAM        if mI_string[0] in ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y']:
        if mI_string[0] in ['B', 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y']:
            polar_Libraries_Dict[mI] = Libraries_Dict[mI]
            print mI, ' in setup_polar_libraries_dict'
    return polar_Libraries_Dict

def get_polar_optimization_exclusions(SCREAM_MODEL):
    polarOptExclusion = SCREAM_MODEL.scream_parameters.getPolarOptimizationExclusions()
    polarExcludeList = []
    POLAR_EXCLUDE = open(polarOptExclusion, 'r')
    for l in POLAR_EXCLUDE.xreadlines():
        if l[0] == '#':
            continue
        f = l.split()
        polarExcludeList.append(f[0])

    return polarExcludeList

def init_core_polar_EE(SCREAM_MODEL, tmpPtn, Libraries_Dict, mutInfo_rotConnInfo_Dict):
    polarOptExclusion = SCREAM_MODEL.scream_parameters.getPolarOptimizationExclusions()
    polarExcludeList = []
    POLAR_EXCLUDE = open(polarOptExclusion, 'r')
    for l in POLAR_EXCLUDE.xreadlines():
        if l[0] == '#':
            continue
        f = l.split()
        polarExcludeList.append(f[0])
    
    (polar_EE, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict) = ('', {}, {})
    
    # 1. Figure out which residues are polar and which aren't.  Currently: assuming all ArbLib residues are non polar.
    for mI in Libraries_Dict.keys():
        mI_string = mI
        if mI_string in polarExcludeList:
            continue
        #if mI_string[0] in ['C', 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y']:
#ADAM        if mI_string[0] in ['D', 'E', 'H', 'J', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y']:
        if mI_string[0] in ['B', 'D', 'E', 'H', 'J', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y']:
            polar_Libraries_Dict[mI] = Libraries_Dict[mI]
            print mI, ' is added to polar EE.'
            # if mutInfo_rotConnInfo_Dict[mI] exists, add to polar_mutInfo_rotConnInfo_Dict (only exists for ArbRotlibs)
            #polar_mutInfo_rotConnInfo_Dict[mI] = mutInfo_rotConnInfo_Dict[mI]

    if len(polar_Libraries_Dict) == 0:
        print ' No Polar residues!'
    else:
        polar_EE = Scream_EE()
        #initialize_scream_EE(tmpPtn, SCREAM_MODEL.scream_parameters, polar_EE, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict)
        init_Scream_EE(polar_EE, SCREAM_MODEL.scream_parameters, tmpPtn, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict)

    return (polar_EE, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict)



def init_polar_EE(SCREAM_MODEL, tmpPtn, Libraries_Dict, mutInfo_rotConnInfo_Dict):
    # Sets up polar_EE, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict.
    (polar_EE, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict) = ('', {}, {})
    
    # 1. Figure out which residues are polar and which aren't.  Currently: assuming all ArbLib residues are non polar.
    for mI in Libraries_Dict.keys():
        mI_string = mI
        #if mI_string[0] in ['C', 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y']:
#ADAM        if mI_string[0] in ['D', 'E', 'H', 'J', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y']:
        if mI_string[0] in ['B', 'D', 'E', 'H', 'J', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y']:
            polar_Libraries_Dict[mI] = Libraries_Dict[mI]
            print mI, ' in init_polar_EE'
            #polar_mutInfo_rotConnInfo_Dict[mI] = mutInfo_rotConnInfo_Dict[mI]

    if len(polar_Libraries_Dict) == 0:
        print ' No Polar residues!'
    else:
        polar_EE = Scream_EE()
        initialize_scream_EE(tmpPtn, SCREAM_MODEL.scream_parameters, polar_EE, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict)

    return (polar_EE, polar_Libraries_Dict, polar_mutInfo_rotConnInfo_Dict)
    
def polar_res_pair_priority(ptn, Libraries_Dict, ResidueReachFile, ResidueExclusionFile):
    # Returns a {} with {distance: [residue 1, residue 2] } structure.
    # Version 0.1: Just order from small to big: priority_dist = total_dist - residue1_dist - residue2_dist.  if priority_dist > 0, skip over.

#ADAM    AA123_map = {'S':'SER', 'T':'THR', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'H':'HIS', 'J': 'HIS', 'K':'LYS', 'N':'ASN', 'Q':'GLN', 'R':'ARG', 'Y':'TYR'}
    AA123_map = {'S':'SER', 'T':'THR', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'H':'HIS', 'J':'HIS', 'B':'HSP', 'K':'LYS', 'N':'ASN', 'Q':'GLN', 'R':'ARG', 'Y':'TYR'}
    non_polar_list = ['A', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W'] # Tyr, His, Cys's are considered polar.  Trp is not.
    # 1. Setup
    reachDict = read_residue_reach_file(ResidueReachFile) #  { S : [2.8, 'CA', 'CB'] }
    excludeList = read_residue_polar_exclusion_file(ResidueExclusionFile)
    polar_map = {} # {mI: [2.8, 'CA', 'CB', CA_Atom, CB_Atom]}
    mI_array = [] # pstn --> mI
    polar_Libraries_Dict = {}
    
    for mI in Libraries_Dict.keys():
        if mI in excludeList or mI[0:1] in non_polar_list:
            continue
        polar_Libraries_Dict[mI] = Libraries_Dict[mI]
        AA3letter = ''
        try:
            AA3letter = AA123_map[mI[0:1]]
        except:
            print ' Error!  No AA123 found', sys.exc_info()[0]
            print mI, ' is the culprit.'
            sys.exit(2)

            
        f = reachDict[AA3letter]
        (atom1, atom2) = (ptn.getAtom(MutInfo(mI), f[1]), ptn.getAtom(MutInfo(mI), f[2]) )
        polar_map[mI] = [ f[0], f[1], f[2], atom1, atom2]
        mI_array.append(mI)
        
    # 2. Build N x N matrix.  Just use a double python array for now--easy.  Also make increasing order list of residues--but only those with value < 0.0 (i.e. sufficiently close to each other)
    close_dist = [] # [[dist, mI1, mI2]]
    dist_matrix = init_matrix(len(mI_array), len(mI_array)) # distance transformed by the the 2.8 in [2.8, 'CA', 'CB'] for each residue (the reach)
    for i in range(0,len(mI_array)):
        mI_i = mI_array[i]
        (mI_i_atom2, mI_i_reach) = (polar_map[mI_i][4], polar_map[mI_i][0])
        for j in range(i+1,len(mI_array)):
            if i == j:
                dist_matrix[i][j] = 999 # i.e. some big number--big number means no need to consider
                continue
            mI_j = mI_array[j]
            (mI_j_atom2, mI_j_reach) = (polar_map[mI_j][4], polar_map[mI_j][0])
            abs_dist = mI_i_atom2.distance(mI_j_atom2)
            dist_minus_reach = abs_dist - mI_i_reach - mI_j_reach
            dist_matrix[i][j] = dist_minus_reach
            dist_matrix[j][i] = dist_minus_reach
            if dist_minus_reach < 0:
                close_dist.append([abs_dist, mI_i, mI_j])

    # dist_matrix data currently not used
    #print_matrix(dist_matrix)
    close_dist.sort(lambda x, y: cmp(x[0], y[0]) ) # sort by abs dist. x, y are lists that look like [dist, mI1, mI2]

    # 3. Filter the close_dist list.  Don't want + charge with + charge residues, dont' want - charge with - charge residues.  Later, don't want surface residues.
    plus_charge_list = ['B', 'K', 'R']
    minus_charge_list = ['E', 'D']
    new_close_dist = []
    for i in range(0,len(close_dist)):
        if (close_dist[i][1][0] in plus_charge_list) and (close_dist[i][2][0] in plus_charge_list):
            continue
        if (close_dist[i][1][0] in minus_charge_list) and (close_dist[i][2][0] in minus_charge_list):
            continue
        new_close_dist.append(close_dist[i])
        
    return (new_close_dist, polar_Libraries_Dict)

    
def print_matrix(mat):
    for row in mat:
        for col in row:
            print col, ' '
        print ''


def init_matrix(m,n):
    # An m by n matrix.  m: rows.  n: columns.
    mat = []
    row = []
    for i in range(0,m):
        for j in range(0,n):
            row.append(0)
        mat.append(row)
        row = []
    return mat

def read_residue_reach_file(filename):
    # Returns: {RES_name_1_letter : [dist, first atom, second atom] }
    # E.g.: { S : [2.8, 'CA', 'CB'] }   so CA-CB is the vector we're gonna be using (in an updated version).

    dict = {}
    INPUT = ''
    try:
        INPUT = open(filename, 'r')
    except e :
        print filename, 'does not exist!'
        return {}
        
    for l in INPUT.readlines():
        if l[0] == '#':
            continue
        f = l.split()
        (res, reach, atom1, atom2) = (f[0].strip(), float(f[1]), f[2].strip(), f[3].strip())
        dict[res] = [reach, atom1, atom2]
    
    INPUT.close()
    return dict

def read_residue_polar_exclusion_file(filename):
    # Input file format:
    # in MutInfo's to begin with!

    list = []
    INPUT = ''
    try:
        INPUT = open(filename, 'r')
    except:
        print filename , 'does not exist!'
        return []
    for l in INPUT.readlines():
        if l[0] == '#':
            continue
        f = l.split()
        mI = f[0].strip() # There should be only 1 field per line anyway.
        list.append(mI)
    
    INPUT.close()
    return list


def initialize_rotlibCollection(rotlibCollection, Libraries_Dict):
    
    for mutInfo in Libraries_Dict.keys():
        library = Libraries_Dict[mutInfo]
        rotlibCollection.addRotlib(mutInfo, library)
    numberOfRotlib = len(Libraries_Dict)

    #if numberOfRotlib <= 10:
    rotlibCollection.initDynamicMemoryDataStructures()
    #else:
        #   rotlibCollection.initAllocationUnderEnergyThreshold(0.1) # 0.1 kcal.
    #print 'Done Initializing RotlibCollection!'
    
def initialize_scream_EE(ptn, SCREAM_PARAMS, scream_EE, Libraries_Dict, mutInfo_rotConnInfo_Dict):
    for mutInfoString in Libraries_Dict.keys():
        mI = MutInfo(mutInfoString)
        for arblib_mI in mutInfo_rotConnInfo_Dict:
            mI.searchAndAddRotConnInfo(MutInfo(arblib_mI), mutInfo_rotConnInfo_Dict[arblib_mI])
        scream_EE.addMutInfoRotConnInfo(mI)

    scream_EE.init_after_addedMutInfoRotConnInfo_on_the_fly_E(ptn, SCREAM_PARAMS)
    _initScreamAtomDeltaValuePy(scream_EE, SCREAM_PARAMS)
    
def init_Scream_EE(scream_EE, SCREAM_PARAMS, ptn, Libraries_Dict, mutInfo_rotConnInfo_Dict):
    '''
    Inits scream_EE using a protein structure and a SCREAM_PARAMs structure.
    '''
    ntrlMutInfo_list = SCREAM_PARAMS.getMutateResidueInfoList()
    additionalLib_list = SCREAM_PARAMS.getAdditionalLibraryInfo()

    if ntrlMutInfo_list.size() == 1 and ntrlMutInfo_list[0] == 'DESIGN':
        ntrlMutInfo_list = []
    elif ntrlMutInfo_list.size() == 1 and ntrlMutInfo_list[0] == 'BINDING_SITE':
        ntrlMutInfo_list = []

    for mutInfo in Libraries_Dict:
        ntrlMutInfo_list.append(mutInfo)

    for mutInfo in ntrlMutInfo_list:
        mI = MutInfo(mutInfo)
        scream_EE.addMutInfoRotConnInfo(mI)
    for mutInfo in mutInfo_rotConnInfo_Dict.keys():
        mI = MutInfo(mutInfo)
        rotConnInfo_z = mutInfo_rotConnInfo_Dict[mutInfo]
        mI.setRotConnInfo(rotConnInfo_z) # 9-08-05: MutInfo() now contains info about RotConnInfo.
        scream_EE.addMutInfoRotConnInfo(mI, rotConnInfo_z)

    ff_file = SCREAM_PARAMS.getOneEnergyFFParFile()
    scream_delta_file = SCREAM_PARAMS.getDeltaParFile()

    #ff_file = 'dreidii322-mpsim-dielectric-2.5.par'
    #scream_delta_file = '/project/Biogroup/Software/SCREAM/lib/SCREAM_delta_par_files/SCREAM_delta.par'
    #scream_EE.init_after_addedMutInfoRotConnInfo(ptn, ff_file, scream_delta_file)
    if SCREAM_PARAMS.getUseRotamerNeighborList() == 'YES':
        scream_EE.init_after_addedMutInfoRotConnInfo_neighbor_list(ptn, SCREAM_PARAMS)
    else:
        scream_EE.init_after_addedMutInfoRotConnInfo_on_the_fly_E(ptn, SCREAM_PARAMS)

    # Then, init FULL pre-initialization.
    _initScreamAtomDeltaValuePy(scream_EE, SCREAM_PARAMS)


def _initScreamAtomDeltaValuePy(scream_EE, SCREAM_PARAMS):
    if SCREAM_PARAMS.getUseDeltaMethod().strip()[0:4] == 'FULL':
        lib_name = SCREAM_PARAMS.getLibResolution()
        if lib_name == 0:
            lib_name = 'SCWRL'
        else:
            if lib_name <= 9:
                lib_name = '0' + str(lib_name)
            else:
                lib_name = str(lib_name)
        method = 'FULL'
        alpha = SCREAM_PARAMS.getDeltaStandardDevs()
        eachAtomDeltaFile = SCREAM_PARAMS.getEachAtomDeltaFile()
        scream_EE.initScreamAtomDeltaValue(lib_name, method, alpha, eachAtomDeltaFile)

def placeBestStructureSoFar(ptn, rotlibCollection, Conformer_RotConnInfo_Dict):
    rotlibCollection.resetTotalEnergyCrntPstn()
    bestRotamers = rotlibCollection.getNextTotalEnergyExcitationRotamersPy()
    placeMultipleConformers_no_EE_reset(ptn, bestRotamers, Conformer_RotConnInfo_Dict)
    groundStateE = rotlibCollection.getEnergyForExcitedRotamers(bestRotamers)
    return (groundStateE, bestRotamers)
    

def printBestELStructure(ptn, BGF_HANDLER, rotlibCollection, Conformer_RotConnInfo_Dict, filename):

    rotlibCollection.resetEmptyLatticeCrntPstn()

    bestELRotamers = rotlibCollection.getNextEmptyLatticeExcitationRotamersPy()
    placeMultipleConformers_no_EE_reset(ptn, bestELRotamers, Conformer_RotConnInfo_Dict)
    BGF_HANDLER.printToFile(filename)
    

def printBestStructures(SCREAM_MODEL, BGF_HANDLER, scream_EE, rotlibCollection, Libraries_Dict, Conformer_RotConnInfo_Dict, Orig_Libraries_Dict, Orig_mutInfo_rotConnInfo_Dict):
      """Prints best structures give a rotlibCollection, with the total number specfied by the Selections parameter in SCREAM param file."""
      print "Printing Best Structures!"
      (ptn, SCREAM_PARAMS) = (SCREAM_MODEL.ptn, SCREAM_MODEL.scream_parameters)
      
      count = 1
      filename_prefix = 'best_'
      rotlibCollection.resetTotalEnergyCrntPstn()
      theseExcitedRotamers = rotlibCollection.getNextTotalEnergyExcitationRotamersPy()

      seqList = []
      For_Residue_E_scream_EE = Scream_EE()
      init_Scream_EE(For_Residue_E_scream_EE, SCREAM_PARAMS, ptn, Orig_Libraries_Dict, Orig_mutInfo_rotConnInfo_Dict)
      
      # First print the one with the best energy.
      filename = filename_prefix + `count` + '.bgf'
      
      placeMultipleConformers_no_EE_reset(ptn, theseExcitedRotamers, Conformer_RotConnInfo_Dict)
      groundStateHandler = bgf_handler(BGF_HANDLER)
      groundStatePtn = Protein(groundStateHandler.getAtomList())
      groundStateE = 0
      
      if SCREAM_PARAMS.getJustOutputSequence() != 'YES':
          # Get energies, which positions have their rotamers changed, and CRMS among the changes.
          # 1. energies.
          groundStateE = rotlibCollection.getEnergyForExcitedRotamers(theseExcitedRotamers)
          # 2. Which ones changed?
          #groundStatePtn.fix_toggle(0)
          groundStateHandler.printToFile(filename, ' Energy: ' + str(groundStateE) + ' kcal/mol Relative to ground state: 0 kcal/mol')
          residue_e_file = 'Residue-E-best_1.txt'
          #Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, scream_EE, residue_e_file)
          Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, For_Residue_E_scream_EE, residue_e_file)
      else:
          seqString = str(count) + " " + BGF_HANDLER.returnSequence()
          seqList.append(seqString)
      theseExcitedRotamers = rotlibCollection.getNextTotalEnergyExcitationRotamersPy()

      # Then proceed with the others
      while len(theseExcitedRotamers) != 0:
          count = count+1
          filename = filename_prefix + `count` + '.bgf'
          residue_e_file = 'Residue-E-best_' + `count` + '.txt'
          if count > SCREAM_PARAMS.Selections:
              break
          #placeMultipleConformers(ptn, theseExcitedRotamers, Conformer_RotConnInfo_Dict)
          placeMultipleConformers_no_EE_reset(ptn, theseExcitedRotamers, Conformer_RotConnInfo_Dict)

          excitedRotamerInfo = rotlibCollection.emptyLatticeRankInfo(theseExcitedRotamers)
          if SCREAM_PARAMS.getJustOutputSequence() != 'YES':
              # Get energies, which positions have their rotamers changed, and CRMS among the changes.
              # 1. energies.
              E = rotlibCollection.getEnergyForExcitedRotamers(theseExcitedRotamers)
              E_string = '%10.4f' % E
              # 2. Which ones changed and CRMS.
              whichOnesChangedAndCRMSString = returnWhichSideChainsChangedAndCRMSString(groundStatePtn, ptn, rotlibCollection, theseExcitedRotamers, Conformer_RotConnInfo_Dict)
              differenceFromGroundStateE = E - groundStateE
              differenceFromGroundStateE_string = '%10.4f' % differenceFromGroundStateE
              #ptn.fix_toggle(0) # screws up results!
              BGF_HANDLER.printToFile(filename, ' Energy: ' + E_string + ' kcal/mol Relative to ground state: ' + differenceFromGroundStateE_string + ' kcal/mol CRMS From Ground State is different for the following SideChains: ' + whichOnesChangedAndCRMSString)
              Print_Final_Residue_Energy_Breakdown(SCREAM_MODEL, For_Residue_E_scream_EE, residue_e_file)

          else:
              seqString = str(count) + " " + BGF_HANDLER.returnSequence()
              seqList.append(seqString)
          theseExcitedRotamers = rotlibCollection.getNextTotalEnergyExcitationRotamersPy()

      # Printing seqList to a file.
      if SCREAM_PARAMS.getJustOutputSequence() == 'YES':
          SeqOutput = open('Sequences.txt', 'w')
          for seq in seqList:
              SeqOutput.write(seq)
              SeqOutput.write('\n')
          SeqOutput.close()
          
      print 'Done printing Best Structures.'

def returnWhichSideChainsChangedAndCRMSString(groundStatePtn, excitedPtn, rotlibCollection, excitedRotamers, Conformer_RotConnInfo_Dict):
    completeBreakdownDict = rotlibCollection.emptyLatticeRankInfo_completeBreakdown(excitedRotamers, Conformer_RotConnInfo_Dict)
    SC_PlacementInfo_String = ''
    
    for k in completeBreakdownDict.keys():
        CRMS = 0
        if k[0] == 'Z':
            # arb lib CRMS calculation
            rotConnInfo = Conformer_RotConnInfo_Dict[k]
            CRMS = groundStatePtn.conformer_CRMS(excitedPtn, rotConnInfo)
        else:
            (AA_type, pstn, chain) = unpackMutInfo(k)
            CRMS = groundStatePtn.sc_CRMS(chain, pstn, excitedPtn)
        
        if CRMS <= 0.001:
            continue
        else:
            CRMS_string = '%5.3f' % CRMS
            SC_PlacementInfo_String = SC_PlacementInfo_String + k + ' ' + CRMS_string + ' '
    return SC_PlacementInfo_String

    


def placeClusters(ptn, clustersDict, clusterRotConnInfoDict = None, scream_EE = None):
    '''
    Places rotamer clusters.
    '''
    # Need to add arblib info to MutInfo, else can't place arb libs.
    mutationFlag = 0
    for mutInfo in clustersDict.keys():
        cluster = clustersDict[mutInfo]
        
        mI = MutInfo(mutInfo)

        for arblib_mI in clusterRotConnInfoDict:
            mI.searchAndAddRotConnInfo(MutInfo(arblib_mI), clusterRotConnInfoDict[arblib_mI])

        mutationFlag += ptn.rotamerClusterPlacement(cluster, mI)
        ptn.setRotamerClusterEmptyLatticeEnergy(cluster, mI, 0) # last field is actually irrelevent

    if scream_EE != None:
        if mutationFlag != 0:
            scream_EE.setup_variableAtomsOnEachSidechain()
            scream_EE.initScreamAtomVdwHbFields()
        

def placeSingleConformer(ptn, mutInfo, conformer, mutInfo_rotConnInfo_Dict, scream_EE):
    '''
    Places a single conformer on protein.  In the case of a mutation, also sets up energy expression if scream_EE is passed in.
    '''

    (mutAA, mutPstn, mutChn) = unpackMutInfo(mutInfo)
    mutation_flag = 0
    
    if mutAA != 'Z':
        mutation_flag = scream_EE.ntrlRotamerPlacement(mutChn, mutPstn, conformer)
        ptn.setEmptyLatticeEnergy(mutChn, mutPstn, conformer.get_empty_lattice_E_abs())
        ptn.setPreCalcEnergy(mutChn, mutPstn, conformer.get_preCalc_TotE() )
        
        if ptn.mutationDone():
            int_map = ptn.getNewMapping()
            for i in mutInfo_rotConnInfo_Dict.keys():
                mutInfo_rotConnInfo_Dict[i].modifyMappingInProteinAtoms(int_map)
    else:
        rotConnInfo_Z = mutInfo_rotConnInfo_Dict[mutInfo]
        ptn.conformerPlacement(conformer, rotConnInfo_Z) 
        
    if scream_EE != None:
        if mutation_flag:
            scream_EE.setup_variableAtomsOnEachSidechain()
            scream_EE.initScreamAtomVdwHbFields()

def placeMultipleConformers_no_EE_reset(ptn, conformersDict, conformerRotConnInfoDict = None):
    '''
    Wrapper for placeMultipleConformers routine that only changes structure but not SC-SC energy expression.
    '''
    placeMultipleConformers(ptn, conformersDict, conformerRotConnInfoDict, None)

def placeMultipleConformers(ptn, conformersDict, conformerRotConnInfoDict = None, scream_EE = None):
    '''
    Places multiple conformers on protein.  In the case of a mutation, also sets up energy expression if scream_EE is passed in.
    Note to self: More general than placeMultipleSideChains.  Should replace all instances of placeMultipleSideChains.
    '''

    mutation_flag = 0
    
    for mutInfo in conformersDict.keys():
        conformer = conformersDict[mutInfo]
        (mutAA, mutPstn, mutChn) = unpackMutInfo(mutInfo)
        
        if mutAA == 'Cluster' and mutChn == 'Cluster':
            clusterDict = {}
            clusterDict[mutInfo] = conformer
            placeClusters(ptn, clusterDict, conformerRotConnInfoDict, scream_EE)
            
        elif mutAA == 'Z' and mutChn == 'Z': # should not enter this conditional if conformersDict contains only rotamers.
            rotConnInfo = conformerRotConnInfoDict[mutInfo]
            ptn.conformerPlacement(conformer, rotConnInfo)
            
        else:
            mutation_flag += ptn.ntrlRotamerPlacement(mutChn, mutPstn, castRotamerToAARotamer(conformer))
            ptn.setPreCalcEnergy(mutChn, mutPstn, conformer.get_preCalc_TotE() )
            ptn.setEmptyLatticeEnergy(mutChn, mutPstn, conformer.get_empty_lattice_E_abs())
            if ptn.mutationDone():
                int_map = ptn.getNewMapping()
                for i in conformerRotConnInfoDict.keys():
                    conformerRotConnInfoDict[i].modifyMappingInProteinAtoms(int_map)

    if scream_EE != None:
        if mutation_flag != 0: # i.e. there has been at least 1 mutation
            scream_EE.setup_variableAtomsOnEachSidechain()
            scream_EE.initScreamAtomVdwHbFields()



def Calc_EL_Energies(ptn, scream_EE, SCREAM_PARAMS, Libraries_Dict, mutInfo_rotConnInfo_Dict, BGF_HANDLER): # temp added BGF_HANDLER
    # This subroutine invokes the desired energy calculation routines.
    # First print which method is being used.
    print '##################################'
    print 'Calculating Empty Lattice Energies'
    print 'Method used: ', SCREAM_PARAMS.getUseDeltaMethod()
    print '##################################'
    if scream_EE == '':
        print 'no scream_EE defined!  Returning.'
        return

    for mutInfo in Libraries_Dict.keys():
        (mutAA, mutPstn, mutChn) = unpackMutInfo(mutInfo)
        mI = MutInfo()
        mI.init(mutInfo)

        library = Libraries_Dict[mutInfo]
        library.reset_pstn()
        currentRotamer = library.get_next_rot()
        count = 1



        print ''
        print mutInfo
        while currentRotamer != None:

            print '' 
            print 'Rotamer', count
            count += 1

            rot_E = 0
            rot_vdw_E  = 0
            rot_hb_E = 0
            rot_coulomb_E = 0
            rot_preCalc_E = 0
            
            if mutAA != 'Z':
                #ptn.ntrlRotamerPlacement(mutChn, mutPstn, currentRotamer)
                scream_EE.ntrlRotamerPlacement(mutChn, mutPstn, currentRotamer)
                if ptn.mutationDone():
                    int_map = ptn.getNewMapping()
                    for i in mutInfo_rotConnInfo_Dict.keys():
                        mutInfo_rotConnInfo_Dict[i].modifyMappingInProteinAtoms(int_map)
                if currentRotamer.get_is_Original_flag() == 1:
                    pass
                else:
                    rot_preCalc_E = currentRotamer.get_preCalc_TotE() 
                    #print 'rot_preCalc_E is : ' + str(rot_preCalc_E)

            else:
                rotConnInfo_Z = mutInfo_rotConnInfo_Dict[mutInfo]
                ptn.conformerPlacement(currentRotamer, rotConnInfo_Z)
                
             # Now decide which SCREAM energy function to use
            if SCREAM_PARAMS.getUseScreamEnergyFunction() == 'YES':
                (method, t) = _returnEnergyMethod(SCREAM_PARAMS)
                
                rot_vdw_E = scream_EE.calc_empty_lattice_vdw_E_delta(mI, method, t)
                rot_hb_E = scream_EE.calc_empty_lattice_hb_E_delta(mI, method, t)
                rot_coulomb_E = scream_EE.calc_empty_lattice_coulomb_E_delta(mI)

                #rot_vdw_hb_exclusion_E = scream_EE.calc_empty_lattice_vdw_hb_exclusion_E_delta(mI, method, t)
                rot_E = rot_vdw_E + rot_hb_E + rot_coulomb_E + rot_preCalc_E
                #+ rot_vdw_hb_exclusion_E
                
                print 'Total:   ' + str(rot_E)
                if SCREAM_PARAMS.getVerbosity() == 1:
                    print 'PreCalc: ' + str(rot_preCalc_E)
                    print 'VDW:     ' + str(rot_vdw_E)
                    print 'HB:      ' + str(rot_hb_E)
                    print 'Coulomb: ' + str(rot_coulomb_E)

                #rot_E = scream_EE.calc_empty_lattice_E(mI)  # obsolete.
                 
            else:
                print 'SCREAM standalone version can only use SCREAM energy functions.'
                sys.exit(2)
            currentRotamer.set_empty_lattice_E_abs(rot_E)
            currentRotamer.set_sc_coulomb_E(rot_coulomb_E)
            #currentRotamer.set_sc_vdw_E(rot_vdw_E + rot_vdw_hb_exclusion_E)
            currentRotamer.set_sc_vdw_E(rot_vdw_E)
            currentRotamer.set_sc_hb_E(rot_hb_E)
            currentRotamer = library.get_next_rot()

    # Then check to see if Absolute Energy > AbsStericClashCutoffEL.
    overall_steric_clash_flag = 1 # 1 means everything passes.  Needs all 1's to pass.
    AbsStericClashCutoffEL = SCREAM_PARAMS.getAbsStericClashCutoffEL()

    for mutInfo in Libraries_Dict.keys():
        (mutAA, mutPstn, mutChn) = unpackMutInfo(mutInfo)
        mI = MutInfo()
        mI.init(mutInfo)
        library = Libraries_Dict[mutInfo]
        library.reset_pstn()
        rot = library.get_next_rot()

        steric_clash_flag = 0 # 0 means none < AbsStericClashCutoffEL.  1 false is enough to fail.
        
        while rot != None:
            if rot.get_empty_lattice_E_abs() < AbsStericClashCutoffEL:
                steric_clash_flag = 1 # if there exist 1 candidate, return true.
                break
            rot = library.get_next_rot()

        if steric_clash_flag == 0:
            print ' No rotamers pass AbsStericClashCutoffEL test for library ' + mutInfo
            
        overall_steric_clash_flag = overall_steric_clash_flag and steric_clash_flag

    if overall_steric_clash_flag == 0:
        print ' Quitting! There are rotamer libraries for which no rotamers pass the AbsStericClashCutoffEL test. '
        sys.exit()

def Calc_This_Interaction_Energy(scream_EE, SCREAM_PARAMS):
    #deltaForInteractionE = SCREAM_PARAMS.getUseDeltaForInterResiE() # YES/NO
    (method, t) = _returnEnergyMethod(SCREAM_PARAMS)

    
    #if deltaForInteractionE == 'YES':
    vdw_E = scream_EE.calc_all_interaction_vdw_E_delta(method, t)
    hb_E = scream_EE.calc_all_interaction_hb_E_delta(method, t)
    coulomb_E = scream_EE.calc_all_interaction_coulomb_E_delta()
    
    #print "VDW Energy: " + str(vdw_E)
    #print "HB Energy: " + str(hb_E)
    #print "Coulomb Energy: " + str(coulomb_E)
    
    total_E = vdw_E + hb_E + coulomb_E
    #print "Total Energy: " + str(total_E)

    return (total_E, vdw_E, hb_E, coulomb_E)



def Print_Library_EL_Energies(Libraries_Dict, Flag=0):
    filename = ''
    i = 1
    c = 1
    MAPFILE = open('MapFile.tbl', 'w')
    for mutInfo in Libraries_Dict.keys():
        filename = ''
        if Flag == 0:
            filename = mutInfo + '.lib.EL'
        if Flag == 1:
            filename = mutInfo + '.HB.lib.EL'
        library = Libraries_Dict[mutInfo]
        library.sort_by_empty_lattice_E()
        library.reset_pstn()
        crntRotamer = library.get_next_rot()

        mI = MutInfo(mutInfo)

        print >>MAPFILE, '%d %s' % (c, mutInfo)

        OUTPUT = ''
        if len(filename) <= 255:
            OUTPUT = open(filename, 'w')
        else:
            if Flag == 0:
                OUTPUT = open('tooLongName%d.lib.EL' % c, 'w')
            if Flag == 1:
                OUTPUT = open('tooLongName%d.HB.lib.E' % c, 'w')
        print >>OUTPUT, '%5s %3s %10s %10s %10s %10s %10s' % ('Rank', 'RES', 'Tot Energy', 'VDW', 'Coulomb', 'HB', 'PreCalc')
        while crntRotamer != None:
            resName = ''
            if mI.isClusterMutInfo():
                resName = ' | '
            elif mutInfo[0] == 'Z':
                resName = mutInfo[0:3]
            else:
                resName = crntRotamer.get_resName()            

            EL_energy = crntRotamer.get_empty_lattice_E()
            EL_coulomb_E = crntRotamer.get_sc_coulomb_E()
            EL_vdw_E = crntRotamer.get_sc_vdw_E()
            EL_hb_E = crntRotamer.get_sc_hb_E()
            EL_preCalc_E = crntRotamer.get_preCalc_TotE() 
            
            print >>OUTPUT, '%5d %3s %10.5f %10.5f %10.5f %10.5f %10.5f' % (i, resName, EL_energy, EL_vdw_E, EL_coulomb_E, EL_hb_E, EL_preCalc_E )
            i = i+1
            crntRotamer = library.get_next_rot()
        OUTPUT.close()
        i = 1
        c += 1
    MAPFILE.close()

def EmptyLatticeExcitationWithArbLib_scream_EE():
    pass

def printEnergyEvolution(rotlibCollection):
    """Prints the estimated and actual energy evolution.  Sorted by estimated energy."""
    Evolution_File = open('Evolution.dat', 'w')
    #rotlibCollection.resetTotalEnergyCrntPstn()
    rotlibCollection.resetEmptyLatticeCrntPstn()
    theseExcitedRotamers = rotlibCollection.getNextEmptyLatticeExcitationRotamersPy()

    info = rotlibCollection.emptyLatticeRankInfo(theseExcitedRotamers)
    mutInfo_ordering = []
    
    header_line = 'Index    Estimated     Actual'
    cc = 1
    for i in info.split():
        if cc%2 == 1:
            header_line = header_line + ' ' + i
            mutInfo_ordering.append(i)
        cc = cc+1
    
    print >>Evolution_File, header_line

    c = 1
    while len(theseExcitedRotamers) != 0:
        rankInfo = rotlibCollection.emptyLatticeRankInfo(theseExcitedRotamers)
        estimated_E = rotlibCollection.getEstimatedEnergyForExcitedRotamers(theseExcitedRotamers)
        actual_E = rotlibCollection.getEnergyForExcitedRotamers(theseExcitedRotamers)
        
        if estimated_E >= 99999999:
            breakmor
        ranking_info_dict = rotlibCollection.emptyLatticeRankInfo_in_dict_form(theseExcitedRotamers)
        ranks_line = ' '
        for mutInfo in mutInfo_ordering:
            n = ranking_info_dict[mutInfo]
            ranks_line = ranks_line + n + ' '

        print >>Evolution_File, '%5d %10.5f %10.5f' % ( c ,estimated_E, actual_E) + ranks_line
        theseExcitedRotamers = rotlibCollection.getNextEmptyLatticeExcitationRotamersPy()
        c += 1
        
    Evolution_File.close()
    

def initRotlibs(SCREAM_PARAMS, ptn):
    '''
    Initializes the following objects:
    self.Libraries_Dict                    # Dict { mutInfo_str : Rotlib* }
    self.Original_Conformers_Dict          # Dict { mutInfo_str : Rotamer* }
    self.Conformer_RotConnInfo_Dict        # Dict { mutInfo_str : rotConnInfo* } for Conformers.
    '''
    (Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict) = return_Ntrl_Library_Dicts(SCREAM_PARAMS, ptn)

    if SCREAM_PARAMS.getKeepOriginalRotamer() == 'YES':
        appendOriginalRotamersToRotlibs(Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict)

    # Next initialize the Arbitrary rotamer libraries.
    additional_list = SCREAM_PARAMS.getAdditionalLibraryInfo()
    (Add_Library_Dict, Add_Orig_Conformer_Dict, mutInfo_rotConnInfo_Dict) = ( {}, {}, {} )

    c = 1
    if len(additional_list) != 0:
        for lib in additional_list:
            # first load library
            library_file = lib
            library = Rotlib(library_file)
            print 'Done loading ' + library_file + '.'
            mutInfo = 'Z' + str(c) + '_Z'
            Add_Library_Dict[mutInfo] = library
            mutInfo_rotConnInfo_Dict[mutInfo] = library.getRotConnInfo()
            # then add original conformer.
            originalConformer = ptn.conformerExtraction(library.getRotConnInfo())
            Add_Orig_Conformer_Dict[mutInfo] = originalConformer
            c += 1

    if SCREAM_PARAMS.getKeepOriginalRotamer() == 'YES':
        appendOriginalConformersToRotlibs(Add_Library_Dict, Add_Orig_Conformer_Dict)
        
    # Put them together.
    Ntrl_Library_Dict.update(Add_Library_Dict)
    Libraries_Dict = Ntrl_Library_Dict.copy()     # Dict { mutInfo : Rotlib* }
    Ntrl_Orig_Rotamer_Dict.update(Add_Orig_Conformer_Dict)
    
    print  'Total Number of Libraries loaded: ' , len(Libraries_Dict)
    return (Libraries_Dict, mutInfo_rotConnInfo_Dict, Ntrl_Orig_Rotamer_Dict)
    
def return_Ntrl_Library_Dicts(SCREAM_PARAMS, ptn):
    # Returns two dictionaries: Ntrl_Library_Dict and Ntrl_Orig_Rotamer_Dict, according to entries in SCREAM_PARAMS.
    ntrlAA_list = SCREAM_PARAMS.getMutateResidueInfoList()
    (Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict) = ( {} , {} )
    
    # Populating Ntrl library and original rotamers.
    if len(ntrlAA_list) != 0:
        if ntrlAA_list[0] == 'DESIGN':
            # setup multiple AA rotamer libraries from getDesignPositionAndClass(), getDesignAAClassDefns(), getDesignClassFromPosition(string), getDesignClassAAs(string).
            # First get list of Design positions, then get their corresponding mutation AA lists, set them up, done.
            library_path_first_half = SCREAM_PARAMS.determineLibDirPath()
            
            DesignPositions = SCREAM_PARAMS.getDesignPositions()
            for dP in DesignPositions:
                DesignClass = SCREAM_PARAMS.getDesignClassFromPosition(dP)
                DesignAAs = SCREAM_PARAMS.getDesignClassAAs(DesignClass)
                Resolution = SCREAM_PARAMS.getLibResolution()
                multipleAARotlib = Multiple_NtrlAARotlib(library_path_first_half, Resolution, DesignAAs)
                mutInfo = convertDesignPositionToMutInfoName(dP)
                Ntrl_Library_Dict[mutInfo] = multipleAARotlib

                (mutAA, mutPstn, mutChn) = unpackMutInfo(mutInfo) 
                originalRotamer = AARotamer()
                originalRotamer.deepcopy(ptn.getAARotamer(mutChn, mutPstn))

                originalRotamer.is_Original = 1
                Ntrl_Orig_Rotamer_Dict[mutInfo] = originalRotamer
                return (Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict)
                
        elif ntrlAA_list[0] == 'BINDING_SITE':
            ntrlAA_list = return_ntrlAA_list_from_BINDING_SITE_mode(SCREAM_PARAMS, ptn)

        (Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict) = Load_Rotamer_Libraries(ntrlAA_list, ptn, SCREAM_PARAMS)
        
    else:
        print 'Warning: No Natural Amino Acids are specified to be SCREAM\'ed!  This is possible when only a ligand rotamer library is to be SCREAM\ed.'
            
    return (Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict)


def return_ntrlAA_list_from_BINDING_SITE_mode(SCREAM_PARAMS, ptn):
    ntrlAA_list = []
    aroundWhatEntity = SCREAM_PARAMS.getBindingSiteMode()
    aroundDistance = SCREAM_PARAMS.getAroundDistance()
    aroundDistanceDefn = SCREAM_PARAMS.getAroundDistanceDefn()
    
    mI_list = []

    if aroundWhatEntity == 'AroundAtom':
        list = SCREAM_PARAMS.getAroundAtom()
        if len(list) == 0:
            print 'AroundAtom not defined!  Please specify AroundAtom parameters.'
            sys.exit()
        mI_list = ptn.residuesAroundAtomN(list, aroundDistance, aroundDistanceDefn)

    if aroundWhatEntity == 'AroundResidue':
        list = SCREAM_PARAMS.getAroundResidue() # vector<MutInfo>
        if len(list) == 0:
            print 'AroundResidue not defined!  Please specify AroundResidue parameters.'
            sys.exit()
        mI_list = ptn.residuesAroundResidue(list, aroundDistance, aroundDistanceDefn) # vector<MutInfo>
        # Now, modify mI_list so that if a mutation is specified in AroundResidue, mI_list reflects this.
        getChnAndPstn = lambda mutInfo : mutInfo.getChn() + str(mutInfo.getPstn())
        mI_list_new = []
        print "##################"
        print mI_list
        map(lambda mI : mI.print_Me(), mI_list)
        print "\n"
        flag = 0 # if already added, has flag 1.  else flag 0.
        for mI in mI_list:
            for aroundResidueMI in list:
                #aroundResidueMI = MutInfo(aroundResidueStr)
                if getChnAndPstn(aroundResidueMI) == getChnAndPstn(mI):
                    print aroundResidueMI.getString()
                    mI_list_new.append(MutInfo(aroundResidueMI.getString()))
                    flag = 1
                    break
            if flag == 0:
                mI_list_new.append(MutInfo(mI.getString()))
            if flag == 1:
                flag = 0
        map(lambda mI : mI.print_Me(), mI_list_new)
        print "#########"
        mI_list = [MutInfo(mI.getString()) for mI in mI_list_new]

    if aroundWhatEntity == 'AroundChain':
        list = SCREAM_PARAMS.getAroundChain()
        if len(list) == 0:
            print 'AroundChain ot defined! Please specify AroundChain parameters.'
        mI_list = ptn.residuesAroundChain(list, aroundDistance, aroundDistanceDefn)

    fixedResidue_list = SCREAM_PARAMS.getFixedResidues()

    for mI in mI_list:
        mI_string = mI.getString()
        mutChn = mI_string[-1]
        if mI_string in fixedResidue_list:
            continue
        if mutChn in fixedResidue_list:
            continue
        if mI_string[0] == 'A' or mI_string[0] == 'G':
            continue
        ntrlAA_list.append(mI_string)

    print 'List of Residues that would be replaced: ',
    for mIStr in ntrlAA_list:
        print mIStr,
    print ''

    return ntrlAA_list
                

def Load_Rotamer_Libraries(ntrlAA_list, ptn, SCREAM_PARAMS):
    (Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict) = ( {}, {} )
    for mutInfo in ntrlAA_list:
        (mutAA, mutPstn, mutChn) = unpackMutInfo(mutInfo)
        #if mutAA == 'G' or mutAA == 'A':
        #continue
        print mutInfo

        library_path_first_half = SCREAM_PARAMS.determineLibDirPath()
        library_file_suffix = SCREAM_PARAMS.determineLibDirFileNameSuffix()
        cnn_file_path = SCREAM_PARAMS.determineCnnDirPath()
        library_file = ''
        # first populate the libraries
        if mutAA == 'J' or mutAA == 'H':
            H_library_file = library_path_first_half + '/' + 'H' + '/' + 'H' + library_file_suffix
            J_library_file = library_path_first_half + '/' + 'J' + '/' + 'J' + library_file_suffix
            H_cnn_file = cnn_file_path + '/H.cnn'
            J_cnn_file = cnn_file_path + '/J.cnn'
            library = HIS_NtrlAARotlib(H_library_file, J_library_file, H_cnn_file, J_cnn_file)
        else:
            library_file = library_path_first_half + '/' + mutAA + '/' + mutAA + library_file_suffix
            cnn_file = cnn_file_path + '/' + mutAA + '.cnn'
            library = NtrlAARotlib(library_file, cnn_file)

        Ntrl_Library_Dict[mutInfo] = library

        # then the original rotamers.
        originalRotamer = AARotamer()
        originalRotamer.deepcopy(ptn.getAARotamer(mutChn, mutPstn))

        originalRotamer.is_Original = 1
        Ntrl_Orig_Rotamer_Dict[mutInfo] = originalRotamer

    return (Ntrl_Library_Dict, Ntrl_Orig_Rotamer_Dict)



def unpackMutInfo(mutInfo):
    """_unpackMutInfo: unpacks a string like C123_X, i.e. returns a (C, 123, X) tuple."""

    if mutInfo.find('|') != -1:
        return ('Cluster', 0, 'Cluster')
    else:
        mutAA = mutInfo[0]
        (mutPstn, mutChn) = mutInfo[1:].split('_')
        mutPstn = int(mutPstn)
        return (mutAA, mutPstn, mutChn)

def convertDesignPositionToMutInfoName(DesignName):
    '''convertDesignPositionToMutInfoName: takes in a DesignName (like A32, which means chain A position 32), returns a legal mutInfo name (like A32_A, alanine, position 32, chain A.  Alanine is set as the default.)'''
    chain = DesignName[0]
    position = DesignName[1:]
    mutInfo = 'A' + position + '_' + chain
    return mutInfo

    
# Two functions that help setup Rotamer libraries.

def appendOriginalRotamersToRotlibs(Mutation_Library_Dict, Original_Rotamer_Dict):
    assert(len(Mutation_Library_Dict) == len(Original_Rotamer_Dict))
    for mutInfo in Original_Rotamer_Dict.keys():
        crntRotlib = Mutation_Library_Dict[mutInfo]
        #crntRotlib.sort_by_rotlib_E()
        crntRotlib.reset_pstn()
        lowestRotlibERot = crntRotlib.get_current_rot()
        #lowest_E = lowestRotlibERot.get_rotlib_E()
        lowest_preCalc_E  = crntRotlib.get_best_preCalc_E()
        
        # need to set an energy for original rotamer... right now just takes the lowset energy from rotlib_E.  need to be improved.
        crntRotamer = Original_Rotamer_Dict[mutInfo]
        crntRotamer.setDeclaredInRotlibScope(0)
        #crntRotamer.set_rotlib_E(lowest_E)
        crntRotamer.set_preCalc_TotE(lowest_preCalc_E)

        # if there is mutation, don't add rotamer.
        if crntRotamer.sameResidueTypeAs(lowestRotlibERot):
            crntRotlib.add_rotamer(crntRotamer) # if not, there is mutation , don't add rotamer.  There might be problems with this in design cases where you may possibly want the original rotamer.  That would have a more complicated decision tree.
    
    
def appendOriginalConformersToRotlibs(Conformer_Library_Dict, Original_Conformer_Dict):
    # Actually identical code to appendOriginalRotamersToRotlibs
    appendOriginalRotamersToRotlibs(Conformer_Library_Dict, Original_Conformer_Dict)


def cavityAnalysis(ptn, scream_EE, probeMutInfo, mutInfo_rotConnInfo_Dict, SCREAM_PARAM):
    # Does simple cavity analysis.
    mI_list = []
    if len(probeMutInfo) == 1:
        mI_list = ptn.residuesAroundChain(probeMutInfo, 3.5, 'SideChainOnly')
    else:
        mI_list = ptn.residuesAroundResidue(probeMutInfo, 3.5, 'SideChainOnly')

    for mI in mI_list:
        scream_EE
    pass

    

# RotlibCollectionPy class, derived from RotlibCollection, added convenience functions.

class RotlibCollectionPy(RotlibCollection):
    "RotlibCollection clas.  Convenience class derived from C++ RotlibCollection to make data transfer easier and cleaner."
    def __init__(self):
        RotlibCollection.__init__(self)
    
    def _shouldKeepGoing(self, count_since_last_best, elapsed_time, allowed_time):
        """Primitive decision whether or not to keep going."""
        #print 'Timing data (elapsed/allowed): ' + str(elapsed_time) + '/' + str(allowed_time)
        if elapsed_time > allowed_time:
            return 0
        
        sizeOfSystem = self.sizeOfSystem()
        testBlockSize = 1.8**(sizeOfSystem+1)

        if testBlockSize < 100:
            testBlockSize = 100
        if testBlockSize > 5000:
            testBlockSize = 5000
        if count_since_last_best > testBlockSize:
            return 0
        else:
            return 1

    def emptyLatticeRankInfo(self, excitedRotamers):
        """Returns a string in the format C3_A 2 D4_A 0 E5_A 5, i.e. name of residue, Empty lattice rank.  ExcitedRotamers is in a format compatible to usage in Python, i.e. PyStr and <C_instance_AARotamer_p>."""
        outputString = ''
        for mutInfo in excitedRotamers.keys():
            rotamer = excitedRotamers[mutInfo]
            Rank = rotamer.get_empty_lattice_energy_rank()
            outputString = outputString + mutInfo + " " + `Rank` + " "
        return outputString

    def emptyLatticeRankInfo_in_dict_form(self, excitedRotamers):
        '''returns a dict with: {'C3_A', 2; 'D4_A', 14; ...} etc.
        '''
        rankInfo_line = self.emptyLatticeRankInfo(excitedRotamers)
        mutInfo_rank_dict = {}
        mutInfo_list = rankInfo_line.split()
        c = 0
        while c < len(mutInfo_list):
            mutInfo_rank_dict[mutInfo_list[c]] = mutInfo_list[c+1]
            c += 2
        return mutInfo_rank_dict

    def emptyLatticeRankInfo_completeBreakdown(self, excitedRotamers, Conformer_RotConnInfo_Dict):
        '''Returns a A1_A B2_B C3_C string instead of a A1_A|B2_B||C3_C string, for instance.'''
        mutInfo_rank_dict = self.emptyLatticeRankInfo_in_dict_form(excitedRotamers)
        completeBreakdownDict = {}
        for k in mutInfo_rank_dict.keys():
            if '|' in k:
                splitList = k.split('|')
                for s in splitList:
                    if s == '':
                        continue
                    else:
                        if s in Conformer_RotConnInfo_Dict.keys():
                            completeBreakdownDict[s] = Conformer_RotConnInfo_Dict[s]
                        else:
                            completeBreakdownDict[s] = None
            else:
                #print 'size of Conformer_RotConnInfo_Dict: ' + str(len(Conformer_RotConnInfo_Dict))
                if k in Conformer_RotConnInfo_Dict.keys():
                    completeBreakdownDict[k] = Conformer_RotConnInfo_Dict[k]
                else:
                    completeBreakdownDict[k] = None
                
        return completeBreakdownDict
    
    def getNextEmptyLatticeExcitationRotamersPy(self ):
        """Returns a dictionary of Python useable {string, AARotamer*} pairs instead of the C++ returned {_p_string, _p_p_AARotamer} format."""
        excitedRotamersRightFormat = dict()
        excitedRotamersWrongFormat = self.getNextEmptyLatticeExcitationRotamers()
        if excitedRotamersWrongFormat.size() == 0:
            return {}
        
        for wrongFormatString in excitedRotamersWrongFormat.keys():
            #rightString = derefString(wrongFormatString)
            rightString = wrongFormatString
            wrongFormatRotamer = excitedRotamersWrongFormat[wrongFormatString]
            #rightRotamer = derefRotamer(wrongFormatRotamer)
            rightRotamer = wrongFormatRotamer
            excitedRotamersRightFormat[rightString] = rightRotamer
        return excitedRotamersRightFormat

    def getNextRotamersByELEnergy_Py(self):
        '''Wrapper for getNextRotamersByELEnergy in RotlibCollection.'''
        excitedRotamersRightFormat = dict()
        excitedRotamersRightFormat = self.getNextRotamersByELEnergy()
        if excitedRotamersRightFormat.size() == 0:
            return {}
        return excitedRotamersRightFormat

    def getNextDynamicMemoryRotamers_And_ExpandPy(self):
        '''Returns a distionary of Python useable {string, Rotamer*} pairs instead of the C++ returned types {_p_string, _p_pAARotamer} or _p_pRotamer format.'''
        excitedRotamersRightFormat = dict()

        #timing.start()
        
        #excitedRotamersWrongFormat = self.getNextDynamicMemoryRotamers_And_Expand()
        excitedRotamersRightFormat = self.getNextDynamicMemoryRotamers_And_Expand()

        #timing.finish()

        #print ' getNextDynamicMemoryRotamers_And_Expand timing: ' + str(timing.micro() / 1000000.00)
        
        #if excitedRotamersWrongFormat.size() == 0:
        if excitedRotamersRightFormat.size() == 0:
            return {}
        # Don't go through following loop if size == 0.

        #timing.start()
        
        #for wrongFormatString in excitedRotamersWrongFormat.keys():
            #rightString = derefString(wrongFormatString)
            #   rightString = wrongFormatString
            #  wrongFormatRotamer = excitedRotamersWrongFormat[wrongFormatString]
            #rightRotamer = derefRotamer(wrongFormatRotamer)
            # rightRotamer = wrongFormatRotamer
            #excitedRotamersRightFormat[rightString] = rightRotamer

        #timing.finish()

        #print ' getNextDynamicMemoryRotamers_And_Expand python string conversion part timing: ' + str(timing.micro() / 1000000.00)
            
        return excitedRotamersRightFormat



    def getNextTotalEnergyExcitationRotamersPy(self):
        """Returns a dictionary of Python useable {string, AARotamer*} pairs instead of the C++ returned {_p_string, _p_p_AARotamer} format from Total Energy dict."""
        excitedRotamersRightFormat = {}
        excitedRotamersWrongFormat = self.getNextTotalEnergyExcitationRotamers()
        if excitedRotamersWrongFormat.size() == 0:
            return {}
        for wrongFormatString in excitedRotamersWrongFormat.keys():
            #rightString = derefString(wrongFormatString)
            rightString = wrongFormatString
            wrongFormatRotamer = excitedRotamersWrongFormat[wrongFormatString]
            #rightRotamer = derefRotamer(wrongFormatRotamer)
            rightRotamer = wrongFormatRotamer
            excitedRotamersRightFormat[rightString] = rightRotamer
        return excitedRotamersRightFormat

    
    def setEnergyForExcitedRotamers(self, excitedRotamers, energy):
        """Takes a ExcitedRotamer structure ( dict{mutInfo, AARotamer*} ) and stores its energy in RotlibCollection."""
        #EE = {}
        EE = ExcitationEnumeration({})
        for mutInfo in excitedRotamers.keys():
            rotamer = excitedRotamers[mutInfo]
            EmptyLatticeRank = rotamer.get_empty_lattice_energy_rank()
            EE[mutInfo] = EmptyLatticeRank
        self.setExcitationEnergy(EE, energy)
        
    def getEnergyForExcitedRotamers(self, excitedRotamers):
        EE = ExcitationEnumeration({})
        for mutInfo in excitedRotamers.keys():
            rotamer = excitedRotamers[mutInfo]
            EmptyLatticeRank = rotamer.get_empty_lattice_energy_rank()
            EE[mutInfo] = EmptyLatticeRank

        energy = self.getExcitationEnergy(EE)
        return energy

    def getEstimatedEnergyForExcitedRotamers(self, excitedRotamers):
        total_estimated_E = 0
        for mutInfo in excitedRotamers.keys():
            rotamer = excitedRotamers[mutInfo]
            this_estimated_E = rotamer.get_empty_lattice_E()
            total_estimated_E += this_estimated_E
        return total_estimated_E


########## Routines that calculate energies on relevent to a single sidechain. ###########

def MutInfoEnergies(ptn, scream_EE, SCREAM_PARAMS, mutInfo):
    mI = MutInfo(mutInfo)
    (method, t) = _returnEnergyMethod(SCREAM_PARAMS)

    chain = mI.getChn()
    pstn = mI.getPstn()
    AA = mI.getAA()

    scream_EE.resetFlags(1)
    scream_EE.initScreamAtomVdwHbFields()
    _initScreamAtomDeltaValuePy(scream_EE, SCREAM_PARAMS) # not 100% sure if necessary

    # Want: PreCalc E, EL_E, Interaction_E.
    # PreCalc Energy.
    PreCalc_E = ptn.getPreCalcEnergy(chain, pstn)
    
    # First, Empty Lattice Energies.
    EL_rot_vdw_E = scream_EE.calc_empty_lattice_vdw_E_delta(mI, method, t)
    EL_rot_hb_E = scream_EE.calc_empty_lattice_hb_E_delta(mI, method, t)
    EL_rot_coulomb_E = scream_EE.calc_empty_lattice_coulomb_E_delta(mI)
    EL_all = EL_rot_vdw_E + EL_rot_hb_E + EL_rot_coulomb_E
    # Second, Interaction Energies.
    scream_EE.fix_all()
    scream_EE.moveable_mutInfo(mI, None, 1)

    All_rot_vdw_E = scream_EE.calc_all_interaction_vdw_E_delta(method, t)
    All_rot_hb_E = scream_EE.calc_all_interaction_hb_E_delta(method, t)
    All_rot_coulomb_E = scream_EE.calc_all_interaction_coulomb_E_delta()
    All_E = All_rot_vdw_E + All_rot_hb_E + All_rot_coulomb_E

    # Total.
    Total_E = PreCalc_E + EL_all + All_E

    return (PreCalc_E, EL_all, EL_rot_vdw_E, EL_rot_coulomb_E,  EL_rot_hb_E, All_E, All_rot_vdw_E, All_rot_coulomb_E, All_rot_hb_E, Total_E)

if __name__ == '__main__':
    main()
