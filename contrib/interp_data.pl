#!/usr/bin/perl -w

use Getopt::Std;
use strict;

# "declare" variables
my ($usage, %Options, $col, $dx, $xstart, $xfinish, $outfile, $datafile);
my ($firstline, $line, @vals, $xin, @valsprev, $xinprev, $xout, $i);

# base 10 logarithm
sub log10 {
    my $inval = shift(@_);
    return(log($inval)/log(10));
}

# default options
$col = 1;
$dx = 1;
$xstart = 0;
$xfinish = 10;
$outfile = "-";

# define usage
$usage = "Interpolates data to a fixed grid in the x variable.

USAGE:
  $0 [options...] <datafile>

OPTIONS:
  -x <col>     : x column in data file [$col]
  -d <dx>      : output interval [$dx]
  -s <xstart>  : start value of x [$xstart]
  -f <xfinish> : finish value of x [$xfinish]
  -o <outfile> : output file [$outfile]
  -h           : display this help text\n";

# get options
if (!getopts('x:d:s:f:o:h', \%Options)) {
    die("$usage");
}

# need to specify output file
if ($#ARGV != 0) {
    die("$usage");
}

# print help if necessary
if ($Options{h}) {
    die("$usage");
}

$datafile = $ARGV[0];

# specified options
if ($Options{x}) {
    $col = $Options{x};
}

if ($Options{d}) {
    $dx = $Options{d};
}

if ($Options{s}) {
    $xstart = $Options{s};
}

if ($Options{f}) {
    $xfinish = $Options{f};
}

if ($Options{o}) {
    $outfile = $Options{o};
}

# open files
if ($outfile =~ "/^-$/") {
    *OP = *STDOUT;
} else {
    open(OP, ">$outfile") or die("Can't open outfile \"$outfile\" for writing.\n");
}

open(FP, "$datafile") or die("Can't open datafile \"$datafile\" for reading.\n");

# read through data file
$firstline = 1;
$xout = $xstart;
while ($line = <FP>) {
    if ($line !~ /^[\s]*\#/) { # skip commented lines
	# special stuff for first line, to force reading of the second line (if it exists)
	if ($firstline) {
	    @vals = split(/[\s]+/, $line);
	    $xin = $vals[$col-1];
	    @valsprev = @vals;
	    $xinprev = $xin;
	    $firstline = 0;
	} else {
	    @valsprev = @vals;
	    $xinprev = $xin;
	    @vals = split(/[\s]+/, $line);
	    $xin = $vals[$col-1];
	}

	while ($xout >= $xinprev && $xout < $xin && $xout <= $xfinish) {
	    if ($xin == $xinprev) {
		die("Error: xin==xprev! xin=$xin xinprev=$xinprev\n");
	    }
	    for ($i=0; $i<=$#vals; $i++) {
		if ($i != $col-1) {
		    printf(OP "%.5g ", ($xout-$xinprev)*($vals[$i]-$valsprev[$i])/($xin-$xinprev)+$valsprev[$i]);
		} else {
		    printf(OP "%.5g ", $xout);
		}
	    }
	    printf(OP "\n");
	    $xout += $dx;
	}
    } else { # print commented lines
	print(OP $line);
    }
}

# close files
if ($outfile !~ "/^-$/") {
    close(OP);
}

close(FP);
