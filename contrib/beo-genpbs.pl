#!/usr/bin/perl -w

use strict;
use Getopt::Std;

# declare variables and subroutines
my (@chars, $randstring, $i, $invokedas, $date);
my ($srcdir, $jobname, $cmdline, $files, $nodes, $cput, $ncpus, $queue, $usage);
my ($outfile, @file);

# generate random string, used for job name
@chars=(0..9, 'A'..'Z', 'a'..'z');
$randstring = '';
for ($i=0; $i<8; $i++) {
    $randstring = $randstring.$chars[int(rand($#chars+1))];
}

# default options
$srcdir = `pwd`;
$srcdir =~ s/\n$//;
$jobname = $randstring;
$cmdline = "echo \"this is a test\"";
$files = "";
$nodes = 1;
$cput = "150:00:00";
$ncpus = 1;
$queue = "mediumjob";

# define usage
$usage = "Generates a PBS file using a few simple rules.

USAGE:
  $0 [options...] <outfile>

OPTIONS:
  -d <srcdir>  : source directory [$srcdir]
  -j <jobname> : job name (random)
  -c <cmdline> : command line [$cmdline]
  -f <files>   : comma-separated list of required files [$files]
  -m <nodes>   : number of nodes required [$nodes]
  -t <cput>    : cput [$cput]
  -n <ncpus>   : ncpus [$ncpus]
  -q <queue>   : queue [$queue]
  -h           : display this help text

EXAMPLE:
  $0 -f cmc,test.input,test.fits -c \"./cmc test.input test.out\" -t 300:00:00 -q longjob test.pbs\n";

# store command line before argument parsing, quoting anything special
# note that the quoting doesn't always work properly
$invokedas = "$0";
for ($i=0; $i<=$#ARGV; $i++) {
    if ($ARGV[$i] =~ /[^a-zA-Z0-9_\-\.\/,]/) {
	$invokedas = $invokedas." "."\"".$ARGV[$i]."\"";
    } else {
	$invokedas = $invokedas." ".$ARGV[$i];
    }
}

# store date
$date = `date`;
$date =~ s/\n$//;

# get options
my %Options;
if (!getopts('d:j:c:f:m:t:n:q:h', \%Options)) {
    die("$usage");
}

# need to specify output file
if ($#ARGV != 0) {
    die("$usage");
}

# print help if necessary
if ($Options{h}) {die("$usage");}

$outfile = $ARGV[0];

# specified options
if ($Options{d}) {$srcdir = $Options{d};}
if ($Options{j}) {$jobname = $Options{j};}
if ($Options{c}) {$cmdline = $Options{c};}
if ($Options{f}) {$files = $Options{f};}
if ($Options{m}) {$nodes = $Options{m};}
if ($Options{t}) {$cput = $Options{t};}
if ($Options{n}) {$ncpus = $Options{n};}
if ($Options{q}) {$queue = $Options{q};}

# open outfile
open(OP, ">$outfile") or die("Can't open outfile \"$outfile\" for writing.\n");

print(OP "#!/bin/sh\n");
print(OP "#\n");
print(OP "# created with $0 on $date\n");
print(OP "# commandline: $invokedas\n");
print(OP "#\n");
print(OP "\n");
print(OP "# PBS directives can't use shell variables, so have to be done by hand.\n");
print(OP "#PBS -N $jobname\n");
print(OP "#PBS -r n\n");
print(OP "#PBS -e $srcdir/$jobname.stderr\n");
print(OP "#PBS -o $srcdir/$jobname.stdout\n");
print(OP "#PBS -l nodes=$nodes,cput=$cput,ncpus=$ncpus\n");
print(OP "#PBS -q $queue\n");
print(OP "\n");
print(OP "# Define variables to make PBS script more flexible.\n");
print(OP "SRCDIR=$srcdir\n");
print(OP "LOCDIR=/usr/local/cluster/users/`whoami`\n");
print(OP "JOBNAME=$jobname\n");
print(OP "\n");
print(OP "# Create working directory based on job name.\n");
print(OP "cd \$LOCDIR\n");
print(OP "mkdir \$JOBNAME\n");
print(OP "cd \$JOBNAME\n");
print(OP "\n");
if ($files !~ /^$/) {
    print(OP "# Copy needed files to working directory.\n");
    @file = split(/,/, $files);
    for ($i=0; $i<=$#file; $i++) {
	$file[$i] =~ s/^([^\/])/\$SRCDIR\/$1/;
	print(OP "cp $file[$i] .\n");
    }
    print(OP "\n");
}
print(OP "# Run code.\n");
print(OP "$cmdline\n");
print(OP "\n");
print(OP "# Copy all files back.\n");
print(OP "rcp * master:\$SRCDIR && rm -f *\n");
print(OP "\n");
print(OP "# And clean up.\n");
print(OP "cd ..\n");
print(OP "rmdir \$JOBNAME\n");
print(OP "\n");
print(OP "# This is just in case the previous command fails for some reason.\n");
print(OP "exit 0\n");

close(OP);

system("chmod 755 $outfile");
