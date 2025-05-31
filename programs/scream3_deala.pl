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
my $radius       = 5;
my $chargetype   = "charmm";

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
    "ALA|CYS|ASP|APP|GLU|GLP|PHE|GLY|HIS|HSE|HSP|".
    "ILE|LYS|LYN|LEU|MET|ASN|PRO|GLN|ARG|ARN|SER|".
    "THR|VAL|TRP|TYR";

if (@ARGV == 0) { help(); }

#ABCDEFGHIJKLMNOPQRSTUVWXYZ
#abcdefghijklmnopqrstuvwxyz

#ABCDEFGHIJKLMNOPQRSTUVWXYZ
#abcdefg ijklmnopqrstuvwxyz
GetOptions ('h|help'          => \$help,
	    'a|alabgf=s'      => \$alabgf,
	    'r|refbgf=s'      => \$refbgf,
	    'l|rotlib=f'      => \$rotlib,
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
if (!$alabgf) {
    die "scream3 :: Must provide BGF to scream!\n";
} elsif (! -e $alabgf) {
    die "scream3 :: Could not find BGF to scream: $alabgf\n";
}
if (!$refbgf) {
    die "scream3 :: Must provide reference BGF!\n";
} elsif (! -e $refbgf) {
    die "scream3 :: Could not find reference BGF: $refbgf\n";
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
    "scream3_deala.pl " . localtime() . "\n".
    " :: Alanized BGF  :: $alabgf\n".
    " :: Reference BGF :: $refbgf\n".
    " :: Rotlib        :: %.1f\n".
    " :: FF File       :: $ff\n".
    " :: Version       :: $screamdir\n".
    " :: Charge type   :: $chargetype\n\n",
    $rotlib / 10;

# Alanized BGF
my %ala = ();
open BGF, "$alabgf";
while (<BGF>) {
    my $line = $_;
    if ($line =~ /^ATOM\s+\d+\s+N\s+(\S+)\s+(\S)\s+(\d+)\s+/) {
	my $res3 = $1;
	my $chn  = $2;
	my $num  = $3;
	my $res1 = res3to1($res3);
	$ala{"${num}_${chn}"} = "$res1";
    }
}
close BGF;

# Reference BGF
my %dat = ();
my @reslist = ();
open BGF, "$refbgf";
while (<BGF>) {
    my $line = $_;
    if ($line =~ /^ATOM\s+\d+\s+N\s+(\S+)\s+(\S)\s+(\d+)\s+/) {
	my $res3 = $1;
	my $chn  = $2;
	my $num  = $3;
	my $res1 = res3to1($res3);

	# If different residue
	if ($ala{"${num}_${chn}"} ne $res1) {
	    push @reslist, "${res1}${num}_${chn}";
	    if ($ala{"${num}_${chn}"} ne "A") {
		$dat{"${res1}${num}_${chn}"} = "[".$ala{"${num}_${chn}"} . " to ${res1}]";
	    } else {
		$dat{"${res1}${num}_${chn}"} = "";
	    }
	}
    }
}
close BGF;

$nres = @reslist;
print "Residues to Dealanize: $nres";
@sortreslist = sort sortres @reslist;
for ($i = 0; $i < @sortreslist; $i++) {
    if ($i % 5 == 0) {
	print "\n  ";
    }
    printf " %8s %8s", $sortreslist[$i], $dat{$sortreslist[$i]};
}
print "\n\n";

# Run scream
$command = "$screamsingle $alabgf $rotlib $ff " . join(" ", @reslist);
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

sub is_screamable {
    my $res = $_[0];

    # Screamable amino acids
    if ($res =~ /^($screamablelist)$/) {
	return 1;
    }

    # Not screamable
    return 0;
}

sub res3to1 {
    my $res = $_[0];
    if (defined $resmap{$res}) {
	return $resmap{$res};
    } else {
	return $res;
    }
}

sub distance {
    my $a = $_[0];
    my $b = $_[1];
    my $c = $_[2];
    my $x = $_[3];
    my $y = $_[4];
    my $z = $_[5];

    return sqrt( ($a - $x)**2 + ($b - $y)**2 + ($c - $z)** 2 );
}

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

############################################################
### Help                                                 ###
############################################################

sub help {

    my $help = "
Program:
 :: scream3_deala.pl

Author:
 :: Adam R. Griffith (griffith\@wag.caltech.edu)

Usage:
 :: scream3_deala.pl -a {alabgf} -r {refbgf} -l {rotlib} -c {charge}

Input:
 :: -a | --alabgf      :: Filename
 :: Alanized BGF

 :: -r | --refbgf      :: Binding site radius
 :: Reference BGF to determine residues to dealanize

 :: -l | --rotlib      :: Rotamer Library
 :: Rotamer library diversity: eg. 0.5, 1.0, 2.0, etc.
 :: Libraries range from 0.1 to 5.0
 :: Default = 0.5

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
 :: Dealanizes the provided BGF based on the residues found
 :: in the reference BGF.
 ::
 :: Note: will also work in cases where valine was used instead
 :: of alanine.
 ::
 :: Warning: Will undo any mutations that are not present in
 :: both the alanized and reference BGFs.
 ::
 :: SCREAM version: $screamdir

";

    die "$help";
}
