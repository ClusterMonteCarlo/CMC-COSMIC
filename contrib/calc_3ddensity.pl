#!/usr/bin/perl -w

# This script will calculate, from an output snapshot file, the 3d mass density
# The script is kind of slow in Perl, but I don't plan on using it much...

use strict;

# declare variables and subroutines
my ($usage, $line, @vals, @r, @m);
my ($n, $mtot, $i, $j, $p, $r, $imin);
my ($imax, $rmin, $rmax, $m, $rho);

# the usage
$usage = 
"Usage: $0
Calculates the 3d mass density profile and writes it to stdout.
You'll have to cat (or zcat) the snapshot file to stdin.
";

# wrong number of arguments?
if ($#ARGV != -1) {
    die("$usage");
}

# read in snapfile
while ($line = <STDIN>) {
    
    @vals = split(/[\s]+/, $line);
    
    push(@r, $vals[2]);
    push(@m, $vals[3]);
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
$p = 20;
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
