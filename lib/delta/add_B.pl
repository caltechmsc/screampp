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

$file = "SCREAM_delta_Total_Min.par";
copy($file, "${file}.bak");

open FILE, "$file";
@lines = <FILE>;
close FILE;

$txt = "";
for ($i = 0; $i < @lines; $i++) {
    $line = $lines[$i];
    if ($line =~ /^\d+ J/) {
	($tmp = $line) =~ s/J/B/;
	$lines[$i] = $line . $tmp;
    }
    $txt .= $lines[$i];
}




open FILE, ">$file";
print FILE $txt;
close FILE;







