#!/usr/bin/perl -w

# This script will calculate, from an output snapshot file, the distribution
# function dN/dEdJ.
# The script is kind of slow in Perl, but I don't plan on using it much...

use strict;

# declare variables and subroutines
my ($usage, $line, @vals, @lines, @sortedlines);
my ($n, $i, $j, $k, $l, $p, $N, $Emax, $Jmax, $Emin, $Jmin, $dE, $dJ);
my ($Eminus, $Eplus, $E, $Jminus, $Jplus, $J, @Jslice, @sortedJslice, $distfunc);

# the usage
$usage = 
"Usage: $0
Calculates the distribution function dN/dEdJ and writes it to stdout.
You'll have to cat (or zcat) the snapshot file to stdin.
";

# wrong number of arguments?
if ($#ARGV != -1) {
    die("$usage");
}

# read in snapfile
while ($line = <STDIN>) {
    
    @vals = split(/[\s]+/, $line);
    
    push(@lines, "$vals[6]:$vals[7]");
}

# this is the total number of stars
$n = $#lines + 1;

# sort lines by E
@sortedlines = sort { (split ':', $a, 2)[0] <=> (split ':', $b, 2)[0] } @lines;

# find Emax, Jmax, Emin, and Jmin
$Emin = (split ':', $sortedlines[0], 2)[0];
$Emax = (split ':', $sortedlines[$n-1], 2)[0];
$Jmin = 100.0;
$Jmax = 0.0;
for ($i=0; $i<$n; $i++) {
    if ((split ':', $sortedlines[$i], 2)[1] >= $Jmax) {
	$Jmax = (split ':', $sortedlines[$i], 2)[1];
    }
    if ((split ':', $sortedlines[$i], 2)[1] <= $Jmin) {
	$Jmin = (split ':', $sortedlines[$i], 2)[1];
    }
}

# create histogram
$p = 50;
$dE = ($Emax-$Emin)/$p;
$dJ = ($Jmax-$Jmin)/$p;
$k = 0;
for ($i=0; $i<$p; $i++) {
    $Eminus = $Emin + $i * ($Emax-$Emin)/$p;
    $Eplus = $Emin + ($i+1) * ($Emax-$Emin)/$p;
    $E = ($Eminus + $Eplus)/2.0;

    # create slice
    @Jslice = ();
    while ($k <= $n-1 && (split ':', $sortedlines[$k], 2)[0] <= $Eplus) {
	push(@Jslice, (split ':', $sortedlines[$k], 2)[1]);
	$k++;
    }

    # sort slice
    @sortedJslice = sort {$a <=> $b} @Jslice;

    # go through slice
    $l = 0;
    for ($j=0; $j<$p; $j++) {
	$Jminus = $Jmin + $j * ($Jmax-$Jmin)/$p;
	$Jplus = $Jmin + ($j+1) * ($Jmax-$Jmin)/$p;
	$J = ($Jminus + $Jplus)/2.0;

	$N = 0;
	while ($l <= $#sortedJslice && $sortedJslice[$l] <= $Jplus) {
	    $N++;
	    $l++;
	}
	
	$distfunc = $N/($dE*$dJ);

	print("$E $J $distfunc\n");
    }
    
    # newline required for Gnuplot
    print("\n");
}
