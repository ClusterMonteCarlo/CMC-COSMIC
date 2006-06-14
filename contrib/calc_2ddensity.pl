#!/usr/bin/perl -w

# This script will calculate, from an output snapshot file, the 2d (projected) mass density
# The script is kind of slow in Perl, but I don't plan on using it much...

use strict;

# declare variables and subroutines
my ($usage, $rh, $line, @vals, @r, @m, @L);
my ($n, $i, $j, $p, $r, $imin);
my ($imax, $rmin, $rmax, $m, $L, @rho, @rhoL);
my ($R, $Rmin, $Rmax, $Rrh, $Sigma, $SigmaL, $z, $zminus);

# base 10 logarithm
sub log10 {
    my $inval = shift(@_);
    return(log($inval)/log(10.0));
}

# the usage
$usage = 
"Usage: $0 <rh>
Calculates the 2d (projected) mass and luminosity density profiles and writes them to stdout.
You'll have to cat (or zcat) the snapshot file to stdin.
";

# wrong number of arguments?
if ($#ARGV != 0) {
    die("$usage");
}

$rh = $ARGV[0];

# read in snapfile
while ($line = <STDIN>) {
    if ($line !~ /^\#/) {
	@vals = split(/[\s]+/, $line);
	
	# mass should already be in MSUN
	push(@m, $vals[1]);
	# this is the Luminosity in LSUN
	push(@L, ($vals[1])**3);
	# this is the radius in units of r_h
	push(@r, $vals[2]/$rh);
    }
}

# this is the total number of stars
$n = $#r + 1;

# calculate 3d density by averaging, and save in array
$p = 10;
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
    $L = $L[$imin]/2.0 + $L[$imax]/2.0;
    for ($j=$imin+1; $j<=$imax-1; $j++) {
	$m += $m[$j];
	$L += $L[$j];
    }
    
    $rho[$i] = $m/(4.0/3.0 * 3.14159265358979312 * ($rmax*$rmax*$rmax - $rmin*$rmin*$rmin));
    $rhoL[$i] = $L/(4.0/3.0 * 3.14159265358979312 * ($rmax*$rmax*$rmax - $rmin*$rmin*$rmin));
}

# calculate 2d density
$p = 100;
$Rmin = 1.0e-4;
$Rmax = $r[$n-1];
for ($i=0; $i<=$p; $i++) {
    $R = $Rmin * 10**($i/$p*log10($Rmax/$Rmin));
    $Sigma = 0.0;
    $SigmaL = 0.0;
    
    # calculate numerical integral Sigma = \int rho(r) dz
    $j = $n-1;
    while ($j >= 0 && $r[$j] >= $R) {
	$z = sqrt($r[$j]*$r[$j] - $R*$R);
	if ($j == 0 || $r[$j-1] < $R) {
	    $zminus = 0.0;
	} else {
	    $zminus = sqrt($r[$j-1]*$r[$j-1] - $R*$R);
	}

	$Sigma += $rho[$j] * ($z - $zminus);
	$SigmaL += $rhoL[$j] * ($z - $zminus);
	$j--;
    }
    
    $Sigma *= 2.0;
    $SigmaL *= 2.0;
    
    print("$R $Sigma $SigmaL\n");
}
