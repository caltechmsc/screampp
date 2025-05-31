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
use Getopt::Long qw(:config no_ignore_case);
use List::Util qw(min max);
use POSIX qw(ceil floor);
use Sys::Hostname;
use Time::Local;

@libs = glob "*.lib";
mkdir "bak";

foreach $lib (@libs) {
    print "$lib\n";
    copy($lib, "bak");

    $txt = "";
    open LIB, $lib;
    while (<LIB>) {
	$line = $_;

	# Should be 2 bonds, 1 lp
	if ($line =~ /\sND1\s/) {
	    substr($line,68,3) = "2 1";
	}

	# Should be 3 bonds, 0 lp
	elsif ($line =~ /\sNE2\s/) {
	    substr($line,68,3) = "3 0";
	}

	$txt .= $line;
    }
    close LIB;
    open LIB, ">$lib";
    print LIB $txt;
    close LIB;
}


#                                                                    68
#0         10        20        30        40        50        60        70
#01234567890123456789012345678901234567890123456789012345678901234567890
#ATOM       9  ND1  HSE A    2    2.22981  -0.63633  -2.18325 N_R    2 1 -0.70000
#ATOM      14  NE2  HSE A    2    2.53966   0.95537  -3.67575 N_R    3 0 -0.36000
