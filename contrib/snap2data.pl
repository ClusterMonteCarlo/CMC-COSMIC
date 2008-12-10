#!/usr/bin/perl -w

use Math::Trig;

my ($PI, $SIGMA, $LSUN, $RSUN);
my (@m, @r, @rproj, @vr, @vt, @startype, @L, @rad);
my (@binflag, @binm0, @binm1, @binstartype0, @binstartype1, @binstarlum0, @binstarlum1, @binstarrad0, @binstarrad1);

# some physical constants
$PI = 3.14159265358;
$SIGMA = 5.67051e-5;
$LSUN = 3.826e33;
$RSUN = 69599000000.0;

# base 10 logarithm
sub log10 {
    my $inval = shift(@_);
    return(log($inval)/log(10.0));
}

# the usage
$usage = 
"Usage: $0 <out.conv.sh filename> <command name> [command parameters]
  Commands:
    extractwds <L_min/L_sun> <L_max/L_sun>
      extract WDs in luminosity range [L_min,L_max]
    extracthrdiag
      extract HR diagram
    extract3dvrms <p>
      extract 3D velocity dispersion as a function of radius in cluster, 
      averaging over p stars for each data point
    extract2dvrms <p>
      extract 3D velocity dispersion as a function of projected radius in
      cluster, averaging over p stars for each data point
    extractwd2dvrms <p> <L_min/L_sun> <L_max/L_sun>
      extract 3D velocity dispersion as a function of projected radius in
      cluster for single WDs, averaging over p WDs for each data point
Note that the out.conv.sh filename must be supplied for physical parameters, and
that the snapshot file must be fed as STDIN to this script.
";

# wrong number of arguments?
if ($#ARGV+1 < 2) {
    die("$usage");
}

$convfile = $ARGV[0];
$commandname = $ARGV[1];

# parse conv file some physical parameters
$lengthunitparsec = 0.0;
$lengthunitcgs = 0.0;
$nbtimeunitcgs = 0.0;
open(FP, "$convfile") or die("Can't open \"$convfile\" for reading.\n");
while ($line = <FP>) {
    if ($line !~ /^[\s]*\#/) { # skip commented lines
	if ($line =~ /^[\s]*([\w]+)=([0-9.\-+eE]+)/) {
	    $token = $1;
	    $value = $2;
	    if ($token =~ /^lengthunitparsec$/) {
		$lengthunitparsec = $value;
	    } elsif ($token =~ /^lengthunitcgs$/) {
		$lengthunitcgs = $value;
	    } elsif ($token =~ /^nbtimeunitcgs$/) {
		$nbtimeunitcgs = $value;
	    }
	}
    }
}
close(FP);
if ($lengthunitparsec == 0.0 || $lengthunitcgs == 0.0 || $nbtimeunitcgs == 0.0) {
    die("Couldn't set all physical parameter values from \"$convfile\".\n");
}

# read in snapfile
while ($line = <STDIN>) {
    if ($line !~ /^\#/) {
	@vals = split(/[\s]+/, $line);
	
	# mass in M_sun
	push(@m, $vals[1]);
	# radius in parsecs
	push(@r, $vals[2] * $lengthunitparsec);
	# projected radius in parsecs
	#$phi = 2.0 * $PI * rand();
	$theta = acos(2.0*rand()-1.0);
	push(@rproj, $vals[2] * $lengthunitparsec * sin($theta));
	# velocities in km/s
	push(@vr, $vals[3] * $lengthunitcgs / $nbtimeunitcgs / 1.0e5);
	push(@vt, $vals[4] * $lengthunitcgs / $nbtimeunitcgs / 1.0e5);
	# stellar type
	push(@startype, $vals[14]);
	# luminosity in L_sun
	push(@L, $vals[15]);
	# radius in R_sun
	push(@rad, $vals[16]);

	# binary flag (0=single, 1=binary)
	push(@binflag, $vals[7]);
	# mass in M_sun
	push(@binm0, $vals[8]);
	push(@binm1, $vals[9]);
	# stellar type
	push(@binstartype0, $vals[17]);
	push(@binstartype1, $vals[18]);
	# luminosity in L_sun
	push(@binstarlum0, $vals[19]);
	push(@binstarlum1, $vals[20]);
	# radius in R_sun
	push(@binstarrad0, $vals[21]);
	push(@binstarrad1, $vals[22]);
    }
}

# this is the total number of stars
$n = $#r + 1;

# extract white dwarves
sub extractwds {
    # wrong number of arguments?
    if ($#ARGV+1 != 4) {
	die("$usage");
    }
    
    $Lmin = $ARGV[2];
    $Lmax = $ARGV[3];

    $i=0;
    while ($i < $n) {
	if ($binflag[$i] == 0 && ($startype[$i] == 10 || $startype[$i] == 11 || $startype[$i] == 12)) {
	    if ($L[$i] >= $Lmin && $L[$i] <= $Lmax) {
		printf("%.8g %.8g %.8g %.8g\n", $r[$i], $rproj[$i], $m[$i], $L[$i]);
	    }
	}
	$i++;
    }
}
# extract HR diagram
sub extracthrdiag {
    # wrong number of arguments?
    if ($#ARGV+1 != 2) {
	die("$usage");
    }
    
    $i = 0;
    while ($i < $n) {
	if ($binflag[$i] == 0) {
	    $temp = ($L[$i]*$LSUN/(4.0*$PI*($rad[$i]*$RSUN)**2*$SIGMA))**0.25;
	    printf("%g %g\n", log10($temp), log10($L[$i]));
	} else {
	    $temp0 = ($binstarlum0[$i]*$LSUN/(4.0*$PI*($binstarrad0[$i]*$RSUN)**2*$SIGMA))**0.25;
	    $temp1 = ($binstarlum1[$i]*$LSUN/(4.0*$PI*($binstarrad1[$i]*$RSUN)**2*$SIGMA))**0.25;
	    $lum0 = $binstarlum0[$i];
	    $lum1 = $binstarlum1[$i];
	    $temp = ($lum0 * $temp0 + $lum1 * $temp1) / ($lum0 + $lum1);
	    printf("%g %g\n", log10($temp), log10($lum0+$lum1));
	}
	$i++;
    }
}

# extract 3D velocity dispersion
sub extract3dvrms {
    # wrong number of arguments?
    if ($#ARGV+1 != 3) {
	die("$usage");
    }
    
    $p = $ARGV[2];
    $i = 0;
    $vrms = 0.0;
    $rave = 0.0;
    $pctr = 0;
    while ($i < $n) {
	if ($binflag[$i] == 0 && $startype[$i] < 10) {
	    $pctr++;
	    $vrms += $vr[$i] * $vr[$i] + $vt[$i] * $vt[$i];
	    $rave += $r[$i];
	} elsif ($binflag[$i] == 1 && ($binstartype0[$i] < 10 || $binstartype1[$i] < 10)) {
	    $pctr++;
	    $vrms += $vr[$i] * $vr[$i] + $vt[$i] * $vt[$i];
	    $rave += $r[$i];
	} 
	
	if ($pctr == $p) {
	    $vrms /= $p;
	    $vrms = sqrt($vrms);
	    $rave /= $p;
	    printf("%g %g\n", $rave, $vrms);
	    $pctr = 0;
	    $vrms = 0.0;
	    $rave = 0.0;
	}
	
	$i++;
    }
}

# extract projected velocity dispersion
sub extract2dvrms {
    # wrong number of arguments?
    if ($#ARGV+1 != 3) {
	die("$usage");
    }
    
    # sort on rproj
    @idx = sort {$rproj[$a] <=> $rproj[$b]} 0 .. $#rproj;

    $p = $ARGV[2];
    $i = 0;
    $vrms = 0.0;
    $rave = 0.0;
    $pctr = 0;
    while ($i < $n) {
	if ($binflag[$idx[$i]] == 0 && $startype[$idx[$i]] < 10) {
	    $pctr++;
	    $vrms += $vr[$idx[$i]] * $vr[$idx[$i]] + $vt[$idx[$i]] * $vt[$idx[$i]];
	    $rave += $rproj[$idx[$i]];
	} elsif ($binflag[$idx[$i]] == 1 && ($binstartype0[$idx[$i]] < 10 || $binstartype1[$idx[$i]] < 10)) {
	    $pctr++;
	    $vrms += $vr[$idx[$i]] * $vr[$idx[$i]] + $vt[$idx[$i]] * $vt[$idx[$i]];
	    $rave += $rproj[$idx[$i]];
	} 
	
	if ($pctr == $p) {
	    $vrms /= $p;
	    $vrms = sqrt($vrms);
	    $rave /= $p;
	    printf("%g %g\n", $rave, $vrms);
	    $pctr = 0;
	    $vrms = 0.0;
	    $rave = 0.0;
	}
	
	$i++;
    }
}

# extract projected velocity dispersion of single WDs
sub extractwd2dvrms {
    # wrong number of arguments?
    if ($#ARGV+1 != 5) {
	die("$usage");
    }
    
    # sort on rproj
    @idx = sort {$rproj[$a] <=> $rproj[$b]} 0 .. $#rproj;

    $p = $ARGV[2];
    $Lmin = $ARGV[3];
    $Lmax = $ARGV[4];
    $i = 0;
    $vrms = 0.0;
    $rave = 0.0;
    $pctr = 0;
    while ($i < $n) {
	if ($binflag[$idx[$i]] == 0 && 
	    ($startype[$idx[$i]] == 10 || $startype[$idx[$i]] == 11 || $startype[$idx[$i]] == 12) &&
	    ($L[$idx[$i]] >= $Lmin && $L[$idx[$i]] <= $Lmax)) {
	    $pctr++;
	    $vrms += $vr[$idx[$i]] * $vr[$idx[$i]] + $vt[$idx[$i]] * $vt[$idx[$i]];
	    $rave += $rproj[$idx[$i]];
	}
	
	if ($pctr == $p) {
	    $vrms /= $p;
	    $vrms = sqrt($vrms);
	    $rave /= $p;
	    printf("%g %g\n", $rave, $vrms);
	    $pctr = 0;
	    $vrms = 0.0;
	    $rave = 0.0;
	}
	
	$i++;
    }
}

# the main attraction
if ($commandname =~ /^extractwds$/) {
    extractwds();
} elsif ($commandname =~ /^extracthrdiag$/) {
    extracthrdiag();
} elsif ($commandname =~ /^extract3dvrms$/) {
    extract3dvrms();
} elsif ($commandname =~ /^extract2dvrms$/) {
    extract2dvrms();
} elsif ($commandname =~ /^extractwd2dvrms$/) {
    extractwd2dvrms();
} else {
    die("$usage");
}


