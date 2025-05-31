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
	    'b|bgf=s'         => \$bgf,
	    'l|rotlib=f'      => \$rotlib,
	    'r|radius=s'      => \$radius,
	    'c|charge=s'      => \$chargetype,
	    'f|ff=s'          => \$ff,
	    'debug'           => \$debug);

if ($help) { help(); }

# Scream Environment
my $screamdir  = "/project/Biogroup/Software/scream3/";
if ($debug) {
    $screamdir = "${Bin}/../";
}
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

# Radius
if ($radius < 3) {
    die "scream3 :: Radius must be >= 3\n";
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
    "scream3_bindsite.pl " . localtime() . "\n".
    " :: BGF File    :: $bgf\n".
    " :: Rotlib      :: %.1f\n".
    " :: FF File     :: $ff\n".
    " :: Site Radius :: $radius\n".
    " :: Version     :: $screamdir\n".
    " :: Charge type :: $chargetype\n\n",
    $rotlib / 10;

@reslist = getbindsite($bgf,$radius);
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

sub getbindsite {
    my $bgf    = $_[0];
    my $radius = $_[1];

    # Read structure
    my %atoms    = ();
    my %residues = ();
    my @ligatoms = ();
    open BGF, "$bgf";
    my @lines = <BGF>;
    close BGF;
    foreach my $line (@lines) {

	# Atom line
	if ($line =~ /^(ATOM|HETATM)/) {
	    my @split = split(/\s+/, $line);
	    $atoms{$split[1]}{x} = $split[6];
	    $atoms{$split[1]}{y} = $split[7];
	    $atoms{$split[1]}{z} = $split[8];

	    # Ligand
	    if ($split[4] eq "X") {
		push @ligatoms, $split[1];
	    }

	    # Not Ligand
	    else {
		push @{$residues{"$split[3]_$split[5]_$split[4]"}}, $split[1];
	    }
	}

	# Format conect
	elsif ($line =~ /^FORMAT CONECT/) {
	    last;
	}
    }

    # Bindsite
    my @bindsite = ();

    # Search by residue
    foreach my $res (keys %residues) {
	$res =~ /(\S+)\_(\d+)\_(\S)/;
	my $res3 = $1; 	my $num  = $2; 	my $chn  = $3;
	my $res1 = res3to1($res3);

	# If not screamable, skip
	if (! is_screamable($res3)) {
	    next;
	}
	
	# Atoms in this residue
	my @resatoms = reverse @{$residues{$res}};

	# Cycle through residue atoms
	foreach my $r_atom (@resatoms) {
	    my $ratom_last = 0;

	    # Cycle through protein atoms
	    foreach my $l_atom (@ligatoms) {

		# Distance
		my $d = sqrt( ($atoms{$r_atom}{x} - $atoms{$l_atom}{x}) ** 2 +
			      ($atoms{$r_atom}{y} - $atoms{$l_atom}{y}) ** 2 +
			      ($atoms{$r_atom}{z} - $atoms{$l_atom}{z}) ** 2 );

		# Check distance
		if ($d < $radius) {
		    push @bindsite, "${res1}${num}_${chn}";

		    # Only need one r_atom <=> $l_atom distance < $radius
		    # Exit @resatoms and @ligatoms loop
		    $ratom_last = 1;
		    last;
		}
	    }

	    # Exit @resatoms loop if one r_atom <=> $l_atom distance < $radius
	    if ($ratom_last) {
		last;
	    }
	}
    }

    # Bindsite
    return @bindsite;
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

    $ENV{"LD_LIBRARY_PATH"} = ":/ul/griffith/ld_lib/libstdc++.so.6".$ENV{"LD_LIBRARY_PATH"};

    die $ENV{"PYTHONPATH"} . "\n";
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
 :: scream3_bindsite.pl

Author:
 :: Adam R. Griffith (griffith\@wag.caltech.edu)

Usage:
 :: scream3_bindsite.pl -b {bgf} -l {rotlib} -r {radius} -c {charge}

Input:
 :: -b | --bgf         :: Filename
 :: BGF file to scream.

 :: -l | --rotlib      :: Rotamer Library
 :: Rotamer library diversity: eg. 0.5, 1.0, 2.0, etc.
 :: Libraries range from 0.1 to 5.0
 :: Default = 0.5

 :: -r | --radius      :: Binding site radius
 :: Default = 5

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
