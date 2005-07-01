#!/usr/bin/perl -w

use strict;

# declare variables and subroutines
my ($usage, $logfile, $tmpfile, $epsfile, $smfile, $megayear, $line, @chars, $randstring);
my (%event, $ncoll, $t, $m, @idi, @mi, $i, $mmax, $idmax, $type);
my ($xmin, $xmax, $ymin, $ymax);

# a subroutine for extracting keyword-value pairs from text
sub grabval {
    my $key = shift(@_);
    my $line = shift(@_);
    $line =~ /$key[\s]*=[\s]*([0-9.\-+eE]+)/;
    return($1);
}

# base 10 logarithm
sub log10 {
    my $inval = shift(@_);
    return(log($inval)/log(10));
}

# draws a point of a given type
sub draw_point {
    my ($x, $y, $size, $type);
    $x = shift(@_);
    $y = shift(@_);
    $size = 1.0 * log10(shift(@_));
    $type = shift(@_);
    
    # keep track of min and max for bounds
    if ($x <= $xmin) {
	$xmin = $x;
    }
    if ($x >= $xmax) {
	$xmax = $x;
    }
    if ($y <= $ymin) {
	$ymin = $y;
    }
    if ($y >= $ymax) {
	$ymax = $y;
    }
    
    # clamp size if too small
    if ($size <= 0.6) {
	$size = 0.6;
    }

    printf(TMPFP "expand %g\n", $size);
    printf(TMPFP "ptype 30 3\n");
    if ($type eq "s") {
	printf(TMPFP "ctype yellow\n");
    } elsif ($type eq "ss") {
	printf(TMPFP "ctype red\n");
    } elsif ($type eq "bs") {
	printf(TMPFP "ctype green\n");
    } elsif ($type eq "bb") {
	printf(TMPFP "ctype blue\n");
    }
    printf(TMPFP "points (%g) (%g)\n", $x, $y);
    printf(TMPFP "ptype 30 0\n");
    printf(TMPFP "ctype black\n");
    printf(TMPFP "points (%g) (%g)\n", $x, $y);
    printf(TMPFP "expand 1.0\n");
}

# draws points and lines for a given event
sub draw_event {
    my ($id, $t, $m, $type, $ncoll, @idi, @mi, $i);
    $id = shift(@_);
    
    if (!exists($event{$id})) {
	die("Event $id does not exist!\n");
    }

    $t = $event{$id}{'t'};
    $m = $event{$id}{'m'};
    $type = $event{$id}{'type'};
    $ncoll = $event{$id}{'ncoll'};
    $idi[1] = $event{$id}{'id1'};
    $idi[2] = $event{$id}{'id2'};
    $idi[3] = $event{$id}{'id3'};
    $idi[4] = $event{$id}{'id4'};
    $mi[1] = $event{$id}{'m1'};
    $mi[2] = $event{$id}{'m2'};
    $mi[3] = $event{$id}{'m3'};
    $mi[4] = $event{$id}{'m4'};
    
    # draw this event's lines, and points of stars without their own events
    for ($i=1; $i<=$ncoll; $i++) {
	printf(TMPFP "relocate %g %g\n", $t, log10($m));
	if (exists($event{$idi[$i]})) {
	    printf(TMPFP "draw %g %g\n", $event{$idi[$i]}{'t'}, log10($mi[$i]));
	} else {
	    printf(TMPFP "draw %g %g\n", $t, log10($mi[$i]));
	    draw_point($t, log10($mi[$i]), $mi[$i], "s");
	}
    }
    
    # here comes the recursion
    for ($i=1; $i<=$ncoll; $i++) {
	if (exists($event{$idi[$i]})) {
	    draw_event($idi[$i]);
	}
    }

    # draw this point last so it is on top
    draw_point($t, log10($m), $m, $type);
}

# the usage
$usage = "Usage: $0 <collision log file> <eps file> <Myr in code time units>\nExtracts merger tree from the collision log file.\n";

# wrong number of arguments?
if ($#ARGV != 2) {
    die("$usage");
}

# the relevant files
$logfile = $ARGV[0];
$epsfile = $ARGV[1];
$smfile = $epsfile;
$smfile =~ s/\.[a-zA-Z]*$/\.sm/;
$megayear = $ARGV[2];

# generate random string, used as part of job directory
@chars=(0..9, 'A'..'Z', 'a'..'z');
for ($i=0; $i<8; $i++) {
    $randstring = $randstring.$chars[int(rand($#chars+1))];
}
$tmpfile = $randstring.".sm";

# read in file and compile events
open(FP, "$logfile");
while ($line = <FP>) {
    if ($line =~ /single-single/) {
	$type = "ss";
    } elsif ($line =~ /binary-single/) {
	$type = "bs";
    } elsif ($line =~ /binary-binary/) {
	$type = "bb";
    }

    if ($line !~ /^\#/) {
	if ($line =~ /m4=/) {
	    $ncoll = 4;
	} elsif ($line =~ /m3=/) {
	    $ncoll = 3;
	} elsif ($line =~ /m2=/) {
	    $ncoll = 2;
	} else {
	    die("Can't determine number of stars merging!\n");
	}
	
	$t = grabval("t", $line) / $megayear;
	$idi[0] = grabval("idm", $line);
	$idi[1] = grabval("id1", $line);
	$idi[2] = grabval("id2", $line);
	$idi[3] = grabval("id3", $line);
	$idi[4] = grabval("id4", $line);
	$mi[0] = grabval("mm", $line);
	$mi[1] = grabval("m1", $line);
	$mi[2] = grabval("m2", $line);
	$mi[3] = grabval("m3", $line);
	$mi[4] = grabval("m4", $line);
	
	$event{$idi[0]}{'t'} = $t;
	$event{$idi[0]}{'m'} = $mi[0];
	$event{$idi[0]}{'type'} = $type;
	$event{$idi[0]}{'ncoll'} = $ncoll;
	$event{$idi[0]}{'id1'} = $idi[1];
	$event{$idi[0]}{'id2'} = $idi[2];
	$event{$idi[0]}{'id3'} = $idi[3];
	$event{$idi[0]}{'id4'} = $idi[4];
	$event{$idi[0]}{'m1'} = $mi[1];
	$event{$idi[0]}{'m2'} = $mi[2];
	$event{$idi[0]}{'m3'} = $mi[3];
	$event{$idi[0]}{'m4'} = $mi[4];
    }
}
close(FP);

# find likely runaway (largest mass merger)
$mmax = -1.0e300;
for my $key ( keys %event ) {
    if ($event{$key}{'m'} >= $mmax) {
	$idmax = $key;
	$mmax = $event{$key}{'m'};
    }
}

# sm file
open(TMPFP, ">$tmpfile");

# this does everything
$xmin = 1.0e300;
$xmax = -1.0e300;
$ymin = 1.0e300;
$ymax = -1.0e300;
draw_event($idmax);

printf(TMPFP "hardcopy\n");
close(TMPFP);

# print sm file
open(SMFP, ">$smfile");

printf(SMFP "device postfile \"%s\"\n", $epsfile);
printf(SMFP "erase\n");
printf(SMFP "expand 1.0\n");
printf(SMFP "lweight 1\n");
printf(SMFP "ticksize 0 0 -1 10\n");
printf(SMFP "limits %g %g %g %g\n", 
       $xmin-($xmax-$xmin)*0.05, $xmax+($xmax-$xmin)*0.05, 
       $ymin-($ymax-$ymin)*0.05, $ymax+($ymax-$ymin)*0.05);
printf(SMFP "xlabel t [Myr]\n");
printf(SMFP "ylabel M [M_{sun}]\n");
printf(SMFP "box\n");
printf(SMFP "ltype 0\n");

close(SMFP);

system("cat $tmpfile >> $smfile");
system("rm $tmpfile");
system("sm < $smfile");

