#!/usr/bin/perl -w

# This script will calculate, from an output snapshot file, the 3d mass density
# The script is kind of slow in Perl, but I don't plan on using it much...

use Getopt::Std;
use strict;

# declare variables and subroutines
my ($usage, %Options, $line, @vals, @r, @m);
my ($n, $mtot, $i, $j, $p, $r, $imin);
my ($imax, $rmin, $rmax, $m, $rho, $radbin);

# base 10 logarithm
sub log10 {
    my $inval = shift(@_);
    return(log($inval)/log(10));
}

# default options
$p = 10;
$radbin = 0;

# the usage
$usage = "Calculates 3D density from snapshot file on stdin and writes it stdout.

USAGE:
  $0 [options...]

OPTIONS:
  -p <p> : bin by star number and set 1/2 of averaging kernel to this value [$p]
  -b <b> : bin by radius and use this many bins [$radbin]
  -h     : display this help text\n";

# get options
if (!getopts('p:b:h', \%Options)) {
    die("$usage");
}

# no required arguments
if ($#ARGV != -1) {
    die("$usage");
}

# print help if necessary
if ($Options{h}) {
    die("$usage");
}

# these options are mutually exclusive
if ($Options{p} && $Options{b}) {
    die("You must specify only one of -p or -b.\n");
}

# specified options
if ($Options{p}) {
    $p = $Options{p};
    $radbin = 0;
}

if ($Options{b}) {
    $radbin = $Options{b};
    $p = 0;
}

# read in snapfile
while ($line = <STDIN>) {
    
    if ($line !~ /^\#/) {
	@vals = split(/[\s]+/, $line);
	
	push(@r, $vals[2]);
	push(@m, $vals[1]);
    }
}

# this is the total number of stars
$n = $#r + 1;

# normalize to N-body units (radii should already be in N-body units)
$mtot = 0.0;
for ($i=0; $i<=$n-1; $i++) {
    $mtot += $m[$i];
}
for ($i=0; $i<=$n-1; $i++) {
    $m[$i] /= $mtot;
}

# calculate 3d density by averaging
if ($p) {
    for ($i=0; $i<=$n-1; $i++) {
	$r = $r[$i];
	
	$imin = $i - $p;
	if ($imin < 0) {
	    $imin = 0;
	}
	
	$imax = $i + $p;
	if ($imax > $n-1) {
	    $imax = $n-1;
	}
	
	$rmin = $r[$imin];
	$rmax = $r[$imax];
	
	# use only half the mass of the stars on the edges
	$m = $m[$imin]/2.0 + $m[$imax]/2.0;
	for ($j=$imin+1; $j<=$imax-1; $j++) {
	    $m += $m[$j];
	}
	
	$rho = $m/(4.0/3.0 * 3.14159265358979312 * ($rmax*$rmax*$rmax - $rmin*$rmin*$rmin));
	
	print("$r $rho\n");
    }
} elsif ($radbin) {
    $n = 0;
    for ($i=0; $i<$radbin; $i++) {
	$rmin = 10**($i/$radbin * (log10($r[$#r])-log10($r[0])) + log10($r[0]));
	$rmax = 10**(($i+1)/$radbin * (log10($r[$#r])-log10($r[0])) + log10($r[0]));
	
	$m = 0;
	while ($n <= $#r && $r[$n] <= $rmax) {
	    $m += $m[$n];
	    $n++;
	}
	
	$r = 10**(0.5 * (log10($rmin) + log10($rmax)));
	$rho = $m/(4.0/3.0 * 3.14159265358979312 * ($rmax*$rmax*$rmax - $rmin*$rmin*$rmin));
	
	print("$r $rho\n");
    }
} else {
    die("Unknown binning scheme.\n");
}
