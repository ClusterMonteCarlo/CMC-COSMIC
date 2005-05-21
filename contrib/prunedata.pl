#!/usr/bin/perl -w

use Getopt::Std;
use strict;

# "declare" variables
my ($usage, %Options, $col, $dx, $dlogx, $n, $outfile, $datafile);
my ($firstline, $line, @vals, $x, $ctr, $xnext, $lastline, $lastlineprinted);

# base 10 logarithm
sub log10 {
    my $inval = shift(@_);
    return(log($inval)/log(10));
}

# default options
$col = 1;
$dx = 0;
$dlogx = 0;
$n = 1;
$outfile = "-";

# define usage
$usage = "Prunes data in a data file.

USAGE:
  $0 [options...] <datafile>

OPTIONS:
  -x <col>     : x column in data file [$col]
  -d <dx>      : approximate output interval [$dx]
  -l <dlogx>   : approximate logarithmic (base 10) ouput interval [$dlogx]
  -n <n>       : print every n lines [$n]
  -o <outfile> : output file [$outfile]
  -h           : display this help text\n";

# get options
if (!getopts('x:d:l:n:o:h', \%Options)) {
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

# these options are mutually exclusive
if (($Options{d} && ($Options{l} || $Options{n})) || 
    ($Options{l} && $Options{n})) {
    die("You must specify only one of -d, -l, or -n.\n");
}

$datafile = $ARGV[0];

# specified options
if ($Options{x}) {
    $col = $Options{x};
}

if ($Options{d}) {
    $dx = $Options{d};
}

if ($Options{l}) {
    $dlogx = $Options{l};
}

if ($Options{n}) {
    $n = $Options{n};
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
while ($line = <FP>) {
    if ($line !~ /^[\s]*\#/) { # skip commented lines
	@vals = split(/[\s]+/, $line);
	$x = $vals[$col-1];
	
	# make sure first line is printed
	if ($firstline) {
	    $ctr = $n-1;
	    $xnext = $x;
	    $firstline = 0;
	}
	
	# and last line, too
	$lastlineprinted = 0;
	
	# three different modes
	if ($Options{n}) {
	    # number mode
	    $ctr++;
	    if ($ctr == $n) {
		$ctr = 0;
		print(OP $line);
		$lastlineprinted = 1;
	    }
	} elsif ($Options{d}) {
	    # dx mode
	    if ($x >= $xnext) {
		$xnext = $x + $dx;
		print(OP $line);
		$lastlineprinted = 1;
	    }
	} else {
	    # dlogx mode
	    if ($x >= $xnext) {
		# can only work with positive values
		if ($x >= 0.0) {
		    $xnext = 10**(log10($x) + $dlogx);
		} else {
		    $xnext = $x;
		}
		print(OP $line);
		$lastlineprinted = 1;
	    }
	}
	
	# save last line
	$lastline = $line;
    } else { # print commented lines
	print(OP $line);
    }
}

# print last line if necessary
if (!$lastlineprinted) {
    print(OP $lastline);
}

# close files
if ($outfile !~ "/^-$/") {
    close(OP);
}

close(FP);
