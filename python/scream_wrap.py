#!/usr/bin/env python

import sys, os

def usage():
    print '''
This is a wrapper for SCREAM.
Usage: %s <bgffile> <rotlib_specs> <ff> <residue list>
E.g.:  %s PheRS.bgf 05 F258_A F260_A
Function: Runs SCREAM using the specified residues.
Rotamer library suggestions:
   05 = 0.5 Angstrom (slower, more fine)
   10 = 1.0 Angstrom (faster, more coarse)
FF: suggestion: /project/Biogroup/FF/dreiding-0.3.par
Selections: 1 (use scream_multi.py for multiple selections)
Delta/value: Non-adjustable: FULL 0.0
''' % (sys.argv[0], sys.argv[0])
    sys.exit()

def main():

    if len(sys.argv) <= 2:
        usage()

    bgf_file          = sys.argv[1]
    rotlib_specs      = sys.argv[2]
    ff_file           = sys.argv[3]
    mutInfo_list      = sys.argv[4:]
    delta_method      = 'FULL'
    final_delta_value = 0.0

    SCREAM_NEW = os.getenv("SCREAM_NEW")
    if SCREAM_NEW == None:
        print ' Error!  enviromental variable SCREAM_NEW not set '
        sys.exit(2)

#    ff_file           = SCREAM_NEW+'/lib/ff/dreiding-0.3.par'
    scream_delta_file = SCREAM_NEW+'/lib/delta/SCREAM_delta_Total_Min.par'
    each_atom_delta   = SCREAM_NEW+'/lib/delta/SCREAM_EachAtomDeltaFileStub.par'
    polar_opt_excl    = SCREAM_NEW+'/lib/delta/SCREAM_PolarOptimizationExclusionsStub.par'
    multi_placement   = 'ClusteringThenDoubletsThenSinglets'

    scream_par_filename = 'scream.par'
    SCREAM_PAR_FILE     = open(scream_par_filename, 'w')
    mutInfo_string      = ''
    extra_lib_string    = ''
    for mutInfo in mutInfo_list:
        import re
        CNN_REGEX = re.compile("cnn")
        if CNN_REGEX.search(mutInfo):\
               extra_lib_string = mutInfo
        else:
            mutInfo_string += mutInfo
            mutInfo_string += ' '
    rotlib_specs_string = ''
    if rotlib_specs == 'SCWRL':
        rotlib_specs_string = rotlib_specs
    else:
        rotlib_specs_string = 'V' + rotlib_specs

    print >>SCREAM_PAR_FILE, '%25s %s' % ('InputFileName',               bgf_file)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('MutateResidueInfo',           mutInfo_string)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('AdditionalLibraryInfo',       extra_lib_string)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('Library',                     rotlib_specs_string)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('PlacementMethod',             'CreateCB')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('CreateCBParameters',          '1.81 51.1 1.55 0.5')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('UseScreamEnergyFunction',     'YES')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('UseDeltaMethod',              delta_method)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('UseDeltaForInterResiE',       'YES')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('UseAsymmetricDelta',          'YES')
    print >>SCREAM_PAR_FILE, '%25s %f' % ('FlatDeltaValue',              final_delta_value)
    print >>SCREAM_PAR_FILE, '%25s %f' % ('DeltaStandardDevs',           final_delta_value)
    print >>SCREAM_PAR_FILE, '%25s %f' % ('InnerWallScalingFactor',      final_delta_value)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('MultiplePlacementMethod',     multi_placement)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('CBGroundSpectrumCalc',        'YES')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('KeepOriginalRotamer',         'NO')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('OneEnergyFFParFile',          ff_file)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('DeltaParFile',                scream_delta_file)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('EachAtomDeltaFile',           each_atom_delta)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('PolarOptimizationExclusions', polar_opt_excl)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('LJOption',                    '12-6')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('CoulombMode',                 'Normal')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('Dielectric',                  '2.5')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('Selections',                  1)
    print >>SCREAM_PAR_FILE, '%25s %s' % ('AbsStericClashCutoffEL',      '999999999')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('StericClashCutoffEnergy',     '250')
    print >>SCREAM_PAR_FILE, '%25s %s' % ('StericClashCutoffDist',       '0.5')
    SCREAM_PAR_FILE.close()

    # Required SCREAM variables
    SCREAM_NEW = os.getenv("SCREAM_NEW")
    if SCREAM_NEW == None:
        print ' Error!  enviromental variable SCREAM_NEW not set '
        sys.exit(2)

    SCREAM_NEW_LIB = os.getenv("SCREAM_NEW_LIB")
    if SCREAM_NEW_LIB == None:
        print ' Error!  enviromental variable SCREAM_NEW_LIB not set '
        sys.exit(2)

    SCREAM_NEW_CNN = os.getenv("SCREAM_NEW_CNN")
    if SCREAM_NEW_CNN == None:
        print ' Error!  enviromental variable SCREAM_NEW_CNN not set '
        sys.exit(2)

    SCREAM_NEW_RTF = os.getenv("SCREAM_NEW_RTF")
    if SCREAM_NEW_RTF == None:
        print ' Error!  enviromental variable SCREAM_NEW_RTF not set '
        sys.exit(2)

    SCREAM_NEW_CHG = os.getenv("SCREAM_NEW_CHG")
    if SCREAM_NEW_CHG == None:
        print ' Error!  enviromental variable SCREAM_NEW_CHG not set '
        sys.exit(2)

    # SCREAM_NEW_SCALE_COU = os.getenv("SCREAM_NEW_SCALE_COU")
    # if SCREAM_NEW_SCALE_COU == None:
    #     print ' Error!  environmental variable SCREAM_NEW_SCALE_COU not set '
    #     sys.exit(2)

    # SCREAM_NEW_SCALE_HB = os.getenv("SCREAM_NEW_SCALE_HB")
    # if SCREAM_NEW_SCALE_HB == None:
    #     print ' Error!  environmental variable SCREAM_NEW_SCALE_HB not set '
    #     sys.exit(2)

    # SCREAM_NEW_SCALE_VDW = os.getenv("SCREAM_NEW_SCALE_VDW")
    # if SCREAM_NEW_SCALE_VDW == None:
    #     print ' Error!  environmental variable SCREAM_NEW_SCALE_VDW not set '
    #     sys.exit(2)

    # Run SCREAM
    run_line_str = SCREAM_NEW + '/python/scream.py' + ' ' + scream_par_filename + ' > ' + 'scream.out'
    os.system(run_line_str)
    sys.exit()

def unpackMutInfo(mutInfo):
    """_unpackMutInfo: unpacks a string like C123_X, i.e. returns a (C, 123, X) tuple."""
    mutAA = mutInfo[0]
    (mutPstn, mutChn) = mutInfo[1:].split('_')
    mutPstn = int(mutPstn)
    return (mutAA, mutPstn, mutChn)


if __name__ == '__main__':
    main()
