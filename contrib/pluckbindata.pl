#!/usr/bin/perl -w

use strict;

# declare variables and subroutines
my ($usage, $zerofile, $threefile, $binfile, $zeroline, $threeline, $binline, @threevals, @binvals);
my ($Ec0, $trh0, $tthree, $tbin, $Mbin, $Mbin0, $Ec, $T, $W, $M, $Ecesc, $rc, $Ebesc);
my ($Ebin, $rhsingle, $rhbin, $DEbb, $DEbs, $Ebb, $Ebs, $fb_core, $fb);
my ($Econs, $vir);

# a subroutine for extracting keyword-value pairs from text
sub grabval {
    my $key = shift(@_);
    my $line = shift(@_);
    $line =~ /$key[\s]*=[\s]*([0-9.\-+eE]+)/;
    return($1);
}

# the usage
$usage = "Usage: $0 <out_prefix>\nDistills the data required for a Gao-type plot.\n";

# wrong number of arguments?
if ($#ARGV != 0) {
    die("$usage");
}

# the relevant files
$zerofile = $ARGV[0]."_0";
$threefile = $ARGV[0]."_3";
$binfile = $ARGV[0]."_binary";

# calculate initial quantities
open(FP, "$zerofile");
$zeroline = <FP>;

# determine initial total energy of cluster
while ($zeroline !~ /^[\s]*Etotal/) {
    $zeroline = <FP>;
}
$Ec0 = grabval("Etotal", $zeroline);

# determine initial trh
while ($zeroline !~ /^[\s]*trh/) {
    $zeroline = <FP>;
}
$trh0 = grabval("trh", $zeroline);

close(FP);

# open files
open(FPthree, "$threefile");
$threeline = <FPthree>;
open(FPbin, "$binfile");
$binline = <FPbin>;

@threevals = split(/[\s]+/, $threeline);
@binvals = split(/[\s]+/, $binline);

$tthree = $threevals[0];
$tbin = $binvals[0];
$Mbin = $binvals[2];

if (abs($tbin-$tthree) >= 1.0e-5*(($tbin+$tthree)/2 + 1)) {
    die("Error: tthree!=tbin in first line: tthree=$tthree tbin=$tbin.\n");
}

$Mbin0 = $Mbin;

# print header
print("#t3/trh0 tbin/trh0 vir Econs Ec T W rc rh,s rh,b M Mbin/Mbin0 Eb/Ec0 DEbb/Ec0 DEbs/Ec0 Ebb/Ec0 Ebs/Ec0 f_b,core f_b\n");

while (($threeline = <FPthree>)&&($binline = <FPbin>)) {
    @threevals = split(/[\s]+/, $threeline);
    @binvals = split(/[\s]+/, $binline);
    
    $tthree = $threevals[0];
    $Ec = $threevals[2];
    $T = $threevals[3];
    $W = $threevals[4];
    $M = $threevals[5];
    $Ecesc = $threevals[7];
    $rc = $threevals[11];
    $Ebesc = $threevals[15];
    
    $tbin = $binvals[0];
    $Mbin = $binvals[2];
    $Ebin = $binvals[3];
    $rhsingle = $binvals[4];
    $rhbin = $binvals[5];
    $DEbb = $binvals[8];
    $DEbs = $binvals[9];
    $Ebb = $binvals[10];
    $Ebs = $binvals[11];
    $fb_core = $binvals[14];
    $fb = $binvals[15];

    # this is the conserved energy in the absence of collisions
    $Econs = $Ec + $Ecesc - $Ebin - $Ebesc;
    $vir = -2 * $T / $W;

    if (abs($tbin-$tthree) >= 1.0e-5*(($tbin+$tthree)/2 + 1)) {
	warn("Warning: tthree!=tbin: tthree=$tthree tbin=$tbin.  Skipping line...\n");
    } else {
	printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
	       $tthree/$trh0, #1
	       $tbin/$trh0, #2
	       $vir, #3
	       $Econs, #4
	       $Ec, #5
	       $T, #6
	       $W, #7
	       $rc, #8
	       $rhsingle, #9
	       $rhbin, #10
	       $M, #11
	       $Mbin/$Mbin0, #12
	       $Ebin/$Ec0, #13
	       $DEbb/$Ec0, #14
	       $DEbs/$Ec0, #15
	       $Ebb/$Ec0, #16
	       $Ebs/$Ec0, #17
	       $fb_core, #18
	       $fb); #19
	       
    }
}

close(FPthree);
close(FPbin);
