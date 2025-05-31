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
my $alatypes     = "FILMYVW";
my $rotlib       = 0.5;
my $chargetype   = "charmm";
my $forceres     = "";
my $skipres      = "";

%resmap = ();
$resmap{ALA} = "A";
$resmap{CYS} = "C";
$resmap{CYX} = "C";
$resmap{ASP} = "D";
$resmap{APP} = "D";
$resmap{GLU} = "E";
$resmap{GLP} = "E";
$resmap{PHE} = "F";
$resmap{GLY} = "G";
$resmap{HIS} = "H";
$resmap{HSE} = "H";
$resmap{HSP} = "B";
$resmap{ILE} = "I";
$resmap{LYS} = "K";
$resmap{LYN} = "K";
$resmap{LEU} = "L";
$resmap{MET} = "M";
$resmap{ASN} = "N";
$resmap{PRO} = "P";
$resmap{GLN} = "Q";
$resmap{ARG} = "R";
$resmap{ARN} = "R";
$resmap{SER} = "S";
$resmap{THR} = "T";
$resmap{VAL} = "V";
$resmap{TRP} = "W";
$resmap{TYR} = "Y";

$screamablelist =
    "CYS|ASP|APP|GLU|GLP|PHE|HIS|HSE|HSP|".
    "ILE|LYS|LYN|LEU|MET|ASN|GLN|ARG|ARN|SER|".
    "THR|VAL|TRP|TYR";

if (@ARGV == 0) { help(); }

#ABCDEFGHIJKLMNOPQRSTUVWXYZ
#abcdefghijklmnopqrstuvwxyz

#ABCDEFGHIJKLMNOPQRSTUVWXYZ
#abcdefg ijklmnopqrstuvwxyz
GetOptions ('h|help'          => \$help,
	    'b|bgf=s'         => \$bgf,
	    'a|alatypes=s'    => \$alatypes,
	    'f|force=s'       => \$forceres,
	    's|skip=s'        => \$skipres,
	    'c|charge=s'      => \$chargetype,
	    'ff=s'            => \$ff,
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

# Residues
@skiplist  = split(/\s+/, $skipres);
@forcelist = split(/\s+/, $forceres);

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

# Read BGF to get residues
@reslist = ();
@alalist = ();
open BGF, "$bgf";
while (<BGF>) {
    $line = $_;
    if ($line =~ /^ATOM\s+\d+\s+N\s+(\S+)\s+(\S+)\s+(\S+)/) {
	$res3 = $1;
	$chn  = $2;
	$num  = $3;
	$res1 = res3to1($res3);

	my $skip = 0;
	foreach $r (@skiplist) {
	    if ("${res1}${num}_${chn}" eq $r) { $skip = 1; }
	}

	my $force = 0;
	foreach $r (@forcelist) {
	    if ("${res1}${num}_${chn}" eq $r) { $force = 1; }
	}

	# Alanizeable & Non-skip
	if (($res1 =~ /[$alatypes]/) && (!$skip)) {
	    push @reslist, "${res1}${num}_${chn}";
	    push @alalist, "A${num}_${chn}";
	}

	# Force
	elsif ($force) {
	    push @reslist, "${res1}${num}_${chn}";
	    push @alalist, "A${num}_${chn}";
	}
    }
}
close BGF;

# Print
printf
    "scream3_ala.pl " . localtime() . "\n".
    " :: BGF File    :: $bgf\n".
    " :: FF File     :: $ff\n".
    " :: Ala Types   :: $alatypes\n".
    " :: Version     :: $screamdir\n".
    " :: Charge type :: $chargetype\n\n",
    $rotlib / 10;

$nskip = @skiplist;
print "Residues to Skip:";
if (@skiplist == 0) {
    print " none\n\n";
} else {
    print " $nskip";
    @sortreslist = sort sortres @skiplist;
    for ($i = 0; $i < @sortreslist; $i++) {
	if ($i % 8 == 0) {
	    print "\n  ";
	}
	printf "%8s", $sortreslist[$i];
    }
    print "\n\n";
}

$nforce = @forcelist;
print "Residues to Force:";
if (@forcelist == 0) {
    print " none\n\n";
} else {
    print " $nforce";
    @sortreslist = sort sortres @forcelist;
    for ($i = 0; $i < @sortreslist; $i++) {
	if ($i % 8 == 0) {
	    print "\n  ";
	}
	printf "%8s", $sortreslist[$i];
    }
    print "\n\n";
}

$nres = @reslist;
print "Residues to Alanize: $nres";
@sortreslist = sort sortres @reslist;
for ($i = 0; $i < @sortreslist; $i++) {
    if ($i % 8 == 0) {
	print "\n  ";
    }
    printf "%8s", $sortreslist[$i];
}
print "\n\n";

# Run scream
$command = "$screamsingle $bgf 10 $ff " . join(" ", @alalist);
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
unlink glob "Residue-E*.txt";

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
### res3to1
############################################################
sub res3to1 {
    my $res = $_[0];
    if (defined $resmap{$res}) {
	return $resmap{$res};
    } else {
	return $res;
    }
}

############################################################
### Help                                                 ###
############################################################

sub help {

    my $help = "
Program:
 :: scream3_ala.pl

Author:
 :: Adam R. Griffith (griffith\@wag.caltech.edu)

Usage:
 :: scream3_ala.pl -b {bgf} -a {ala types} -f {force residues}
 ::                -s {skip residues} -c {chargetype}

Input:
 :: -b | --bgf         :: Filename
 :: BGF file to scream.

 :: -a | --alatypes    :: String
 :: Types of residues to alanize.  Default: FILMYVW

 :: -f | --force       :: SCREAM residue list
 :: Force residues to be alanized.  Useful if you want to alanize
 :: a specific polar residue.
 :: Example: -f 'R103_4 R228_5'

 :: -s | --skip        :: SCREAM residue list
 :: Prevent residues from being alanized.  Useful if you want to
 :: keep a particular non-polar residue.
 :: Example: -s 'W220_6 F288_6 F289_6'

 :: --ff               :: Forcefield File
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
 :: Alanizes nonpolar residues.
 :: SCREAM version: $screamdir

";

    die "$help";
}
