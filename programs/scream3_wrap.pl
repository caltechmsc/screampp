#!/usr/bin/perl

BEGIN {
    use FindBin qw($Bin);
    use lib $FindBin::Bin;
}

use Cwd;
use File::Basename;
use File::Copy;
use File::Path;
use File::Spec;
use Getopt::Long qw(:config no_ignore_case pass_through);
use List::Util qw(min max);
use POSIX qw(ceil floor);
use Sys::Hostname;
use Time::Local;

# Variables
my $ff           = "/project/Biogroup/FF/dreiding-0.4-x6.par";
my $rotlib       = 0.5;
my $chargetype   = "charmm";

if (@ARGV == 0) { help(); }

#ABCDEFGHIJKLMNOPQRSTUVWXYZ
#abcdefghijklmnopqrstuvwxyz

#ABCDEFGHIJKLMNOPQRSTUVWXYZ
#abcdefg ijklmnopqrstuvwxyz
GetOptions ('h|help'          => \$help,
	    'b|bgf=s'         => \$bgf,
	    'l|rotlib=f'      => \$rotlib,
	    'r|residues=s'    => \$residues,
	    'c|charge=s'      => \$chargetype,
	    'f|ff=s'          => \$ff,
	    'debug'           => \$debug);

if ($help) { help(); }

# Scream Environment
my $screamdir  = "/project/Biogroup/Software/scream3/";
if ($debug) {
    $screamdir = "${Bin}/../"; }
my $scream       = "${screamdir}/python/scream.py";
my $screamsingle = "${screamdir}/python/scream_wrap.py";
my $screammulti  = "${screamdir}/python/scream_multi.py";

############################################################
### Main Routine                                         ###
############################################################

# BGF
if (!$bgf) {
    die "scream3 :: Must provide BGF to scream!\n";
} elsif (! -e $bgf) {
    die "scream3 :: Could not find BGF to scream: $bgf\n";
}

# Rotlib
if (($rotlib < 0.1) || ($rotlib > 5.0)) {
    die "scream3 :: Rotamer library must be: 0.1 <= rotlib <= 5.0\n";
}
$rotlib = sprintf "%.1f", $rotlib;
if ($rotlib < 1) {
    $rotlib = sprintf "%d", $rotlib * 10;
    $rotlib = "0" . $rotlib;
} else {
    $rotlib = sprintf "%d", $rotlib * 10;
}

# Residues
if (!$residues) {
    die "scream3 :: No residues specified for SCREAMing!\n";
}
@reslist = split(/\s+/, $residues);
if (@reslist < 1) {
    die "scream3 :: No residues specified for SCREAMing!\n";
}

# Charge
if ($chargetype !~ /^(charmm|charmm_n|amber|amber_n|qeq|qeq_n)$/) {
    die "scream3 :: Charge type must be: charmm, charmm_n, amber, amber_n\n".
	"        :: qeq, or qeq_n!\n";
}

# FF
if (! -e $ff) {
    die "scream3 :: Could not find forcefield file: $ff\n";
}

# Scream setup
loadscream();

# Print
printf
    "scream3_wrap.pl " . localtime() . "\n".
    " :: BGF File    :: $bgf\n".
    " :: Rotlib      :: %.1f\n".
    " :: FF File     :: $ff\n".
    " :: Version     :: $screamdir\n".
    " :: Charge type :: $chargetype\n\n",
    $rotlib / 10;
$nres = @reslist;
print "Residues: $nres";
@sortreslist = sort sortres @reslist;
for ($i = 0; $i < @sortreslist; $i++) {
    if ($i % 8 == 0) {
	print "\n  ";
    }
    printf "%8s", $sortreslist[$i];
}
print "\n\n";

# Run scream
$command = "$screamsingle $bgf $rotlib $ff " . join(" ", @reslist);
`$command >& scream.out`;

# Success
if (-e "best_1.bgf") {
    `/project/Biogroup/scripts/perl/BGFFormat.pl -b best_1.bgf`;
    print "Scream completed successfully!\n\n";
}

# Fail
else {
    print "ERROR :: Scream failed!\n\n";
}

# Cleanup
unlink "timing.txt";
unlink "Anneal-Energies.txt";
unlink "Field1.bgf";
unlink "Residue-E.txt";

exit;

sub sortres {
    $a =~ /(\d+)\_(\S)/;
    $a_num = $1;
    $a_chn = $2;

    $b =~ /(\d+)\_(\S)/;
    $b_num = $1;
    $b_chn = $2;

    if (($a_chn <=> $b_chn) != 0) {
	return $a_chn <=> $b_chn;
    }

    return $a_num <=> $b_num;
}

sub loadscream {
    my $arch = `uname -m`; chomp $arch;
    if ($arch =~ /\_64$/) {
	$ENV{"python"} = "/exec/python/pythonEPD-7.0-2-rh5-x86_64/bin/python";
	$ENV{"PATH"} =
	    "${screamdir}/python:/exec/python/pythonEPD-7.0-2-rh5-x86_64/bin:" . $ENV{"PATH"};
    } else {
	$ENV{"python"} = "/exec/python/pythonEPD-7.0-2-rh3-x86/bin/python";
	$ENV{"PATH"} =
	    "${screamdir}/python:/exec/python/pythonEPD-7.0-2-rh3-x86/bin:" . $ENV{"PATH"};
    }

    $ENV{"SCREAM_NEW"}          = "${screamdir}";
    $ENV{"SCREAM_NEW_CHG"}      = "${chargetype}";
    $ENV{"SCREAM_NEW_LIB"}      = "${screamdir}/lib/${chargetype}/";
    $ENV{"SCREAM_NEW_CNN"}      = "${screamdir}/lib/cnn/";
    $ENV{"SCREAM_NEW_RTF"}      = "${screamdir}/lib/rft/";

    if ($arch =~ /\_64$/) {
	$ENV{"PYTHONPATH"} =
	    "${screamdir}/build/lib.linux-x86_64-2.7".
	    ":${screamdir}/python/packages".
	    ":" . $ENV{"PYTHONPATH"};
    } else {
	$ENV{"PYTHONPATH"} =
	    "${screamdir}/build/lib.linux-i686-2.7".
	    ":${screamdir}/python/packages".
	    ":" . $ENV{"PYTHONPATH"};
    }
}

############################################################
### Help                                                 ###
############################################################

sub help {

    my $help = "
Program:
 :: scream3_wrap.pl

Author:
 :: Adam R. Griffith (griffith\@wag.caltech.edu)

Usage:
 :: scream3_wrap.pl -b {bgf} -l {rotlib} -r {residues} -c {charge}

Input:
 :: -b | --bgf         :: Filename
 :: BGF file to scream.

 :: -l | --rotlib      :: Rotamer Library
 :: Rotamer library diversity: eg. 0.5, 1.0, 2.0, etc.
 :: Libraries range from 0.1 to 5.0
 :: Default = 0.5

 :: -r | --residues    :: List of SCREAM residues
 :: Residues to SCREAM.
 :: Example: -r 'D103_C F288_6'

 :: -f | --ff          :: Forcefield File
 :: Alternate forcefield file for energies.
 :: User specified FF must match format of default
 :: Default: /project/Biogroup/FF/dreiding-0.4-x6.par

 :: -c | --charge      :: Charge keyword
 :: charmm, charmm_n, amber, amber_n, qeq, or qeq_n
 :: Default = charmm

 :: -h | --help        :: No Input
 :: Displays this help message.

Description:
 :: Runs SCREAM using the provided input.
 :: SCREAM version: $screamdir

";

    die "$help";
}
