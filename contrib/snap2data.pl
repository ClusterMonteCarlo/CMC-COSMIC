#!/usr/bin/perl -w

use Math::Trig;
use List::Util qw[min max];

my ($PI, $SIGMA, $LSUN, $RSUN);
my (@m, @r, @rproj, @vr, @vt, @startype, @L, @rad, @id, @id0, @id1);
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

# minimum of two quantities
sub mymin {
    my $x = shift(@_);
    my $y = shift(@_);
    if ($x <= $y) {
	return($x);
    } else {
	return($y);
    }
}

# maximum of two quantities
sub mymax {
    my $x = shift(@_);
    my $y = shift(@_);
    if ($x >= $y) {
	return($x);
    } else {
	return($y);
    }
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
    extractobsbinfrac <q_crit> <r_min> <r_max>
      extract observable binary fraction, estimated as number of MS-MS
      binaries with q>q_crit relative to number of stars appearing on or near MS, 
      for cluster overall, within r_min<r<r_max (with r in parsecs), and within
      the 10% Lagrange radius
    extractmss <M_min/M_sun> <M_max/M_sun>
      extract MS stars in mass range [M_min,M_max]
Note that the out.conv.sh filename must be supplied for physical parameters, and
that the snapshot file must be fed as STDIN to this script.
";

# wrong number of arguments?
if ($#ARGV+1 < 2) {
    die("$usage");
}

$convfile = $ARGV[0];
$commandname = $ARGV[1];

# parse conv file for some physical parameters
$lengthunitparsec = 0.0;
$lengthunitcgs = 0.0;
$massunitcgs = 0.0;
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
	    } elsif ($token =~ /^massunitcgs$/) {
		$massunitcgs = $value;
	    } elsif ($token =~ /^nbtimeunitcgs$/) {
		$nbtimeunitcgs = $value;
	    }
	}
    }
}
close(FP);

if ($lengthunitparsec == 0.0 || $lengthunitcgs == 0.0 || $massunitcgs == 0.0) {
    die("Couldn't set all physical parameter values from \"$convfile\".\n");
}
if ($nbtimeunitcgs == 0.0) {
    $nbtimeunitcgs = $lengthunitcgs**1.5 / sqrt(6.67259e-08 * $massunitcgs);
}

# read in snapfile
while ($line = <STDIN>) {
    if ($line !~ /^\#/) {
	@vals = split(/[\s]+/, $line);
	
	#ids
	push(@id, $vals[0]);
	push(@id0,$vals[10]);
	push(@id1, $vals[11]);	
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
    my @m_more;
    
    # wrong number of arguments?
    if ($#ARGV+1 != 2) {
	die("$usage");
    }
    
    $i = 0;
    while ($i < $n) {
	if ($binflag[$i] == 0) {
		if ($rad[$i]>0.0){
	    		$temp = ($L[$i]*$LSUN/(4.0*$PI*($rad[$i]*$RSUN)**2*$SIGMA))**0.25;
	    		printf("%g %g\n",log10($temp), log10($L[$i]));
		} 
	} else {
		if ($binstarrad0[$i]>0.0 && $binstarrad1[$i]>0.0){
	    		$temp0 = ($binstarlum0[$i]*$LSUN/(4.0*$PI*($binstarrad0[$i]*$RSUN)**2*$SIGMA))**0.25;
	    		$temp1 = ($binstarlum1[$i]*$LSUN/(4.0*$PI*($binstarrad1[$i]*$RSUN)**2*$SIGMA))**0.25;
	    		$lum0 = $binstarlum0[$i];
	    		$lum1 = $binstarlum1[$i];
	    		$temp = ($lum0 * $temp0 + $lum1 * $temp1) / ($lum0 + $lum1);
	    		$m_more[$i] = max($binm0[$i],$binm1[$i]);
	    		printf("%g %g\n",log10($temp), log10($lum0+$lum1));
		}
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

# extract observed binary fraction
sub extractobsbinfrac {
    # wrong number of arguments?
    if ($#ARGV+1 != 5) {
	die("$usage");
    }
    
    $qcrit = $ARGV[2];
    $rmin = $ARGV[3];
    $rmax = $ARGV[4];

    # first calculate total cluster mass
    $mtotal = 0.0;
    $i = 0;
    while ($i < $n) {
	$mtotal += $m[$i];
	$i++;
    }

    # calculate binary fractions for whole cluster
    $i = 0;
    $mcount = 0.0;
    $nbin = 0;
    $ntot = 0;
    $nbinobs = 0;
    $ntotobs = 0;
    while ($i < $n) {
	if ($binflag[$i] == 1) {
	    $nbin++;
	    # a binary with a large enough mass ratio gets counted as a binary
	    if (($binstartype0[$i] == 0 || $binstartype0[$i] == 1) && ($binstartype1[$i] == 0 || $binstartype1[$i] == 1)) {
		$ntotobs++;
		$q = mymin($binm0[$i], $binm1[$i])/mymax($binm0[$i], $binm1[$i]);
		if ($q >= $qcrit) {
		    $nbinobs++;
		}
	    }
	    # a MS-compact object binary gets counted as a single MS object if it's not too far away from the MS
	    if ( (($binstartype0[$i] == 0 || $binstartype0[$i] == 1) && $binstartype1[$i] >= 10) ||
		 (($binstartype1[$i] == 0 || $binstartype1[$i] == 1) && $binstartype0[$i] >= 10) ) {
		# check for special zero radius case
		if ($binstarrad0[$i] != 0.0) {
		    $temp0 = ($binstarlum0[$i]*$LSUN/(4.0*$PI*($binstarrad0[$i]*$RSUN)**2*$SIGMA))**0.25;
		} else {
		    $temp0 = 0.0;
		}
		$lum0 = $binstarlum0[$i];
		if ($binstarrad1[$i] != 0.0) {
		    $temp1 = ($binstarlum1[$i]*$LSUN/(4.0*$PI*($binstarrad1[$i]*$RSUN)**2*$SIGMA))**0.25;
		} else {
		    $temp1 = 0.0;
		}
		$lum1 = $binstarlum1[$i];
		$temp = ($lum0 * $temp0 + $lum1 * $temp1) / ($lum0 + $lum1);
		$lum = $lum0 + $lum1;
		if ($binstartype0[$i] == 0 || $binstartype0[$i] == 1) {
		    $reftemp = $temp0;
		    $reflum = $lum0;
		} else {
		    $reftemp = $temp1;
		    $reflum = $lum1;
		}
		if (2.0*abs($temp-$reftemp)/($temp+$reftemp) <= 0.1 &&
		    2.0*abs($lum-$reflum)/($lum+$reflum) <= 0.1) {
		    $ntotobs++;
		}
	    }
	} else {
	    if ($startype[$i] == 0 || $startype[$i] == 1) {
		$ntotobs++;
	    }
	}
	$ntot++;
	$i++;
    }
    printf("total cluster: f_b=%g+/-%g f_b,obs=%g+/-%g\n", $nbin/$ntot, sqrt($nbin)/$ntot, $nbinobs/$ntotobs, sqrt($nbinobs)/$ntotobs);

    $i = 0;
    $mcount = 0.0;
    $nbin = 0;
    $ntot = 0;
    $nbinobs = 0;
    $ntotobs = 0;
    while ($mcount < 0.1 * $mtotal) {
	if ($binflag[$i] == 1) {
	    $nbin++;
	    # a binary with a large enough mass ratio gets counted as a binary
	    if (($binstartype0[$i] == 0 || $binstartype0[$i] == 1) && ($binstartype1[$i] == 0 || $binstartype1[$i] == 1)) {
		$ntotobs++;
		$q = mymin($binm0[$i], $binm1[$i])/mymax($binm0[$i], $binm1[$i]);
		if ($q >= $qcrit) {
		    $nbinobs++;
		}
	    }
	    # a MS-compact object binary gets counted as a single MS object if it's not too far away from the MS
	    if ( (($binstartype0[$i] == 0 || $binstartype0[$i] == 1) && $binstartype1[$i] >= 10) ||
		 (($binstartype1[$i] == 0 || $binstartype1[$i] == 1) && $binstartype0[$i] >= 10) ) {
		# check for special zero radius case
		if ($binstarrad0[$i] != 0.0) {
		    $temp0 = ($binstarlum0[$i]*$LSUN/(4.0*$PI*($binstarrad0[$i]*$RSUN)**2*$SIGMA))**0.25;
		} else {
		    $temp0 = 0.0;
		}
		$lum0 = $binstarlum0[$i];
		if ($binstarrad1[$i] != 0.0) {
		    $temp1 = ($binstarlum1[$i]*$LSUN/(4.0*$PI*($binstarrad1[$i]*$RSUN)**2*$SIGMA))**0.25;
		} else {
		    $temp1 = 0.0;
		}
		$lum1 = $binstarlum1[$i];
		$temp = ($lum0 * $temp0 + $lum1 * $temp1) / ($lum0 + $lum1);
		$lum = $lum0 + $lum1;
		if ($binstartype0[$i] == 0 || $binstartype0[$i] == 1) {
		    $reftemp = $temp0;
		    $reflum = $lum0;
		} else {
		    $reftemp = $temp1;
		    $reflum = $lum1;
		}
		if (2.0*abs($temp-$reftemp)/($temp+$reftemp) <= 0.1 &&
		    2.0*abs($lum-$reflum)/($lum+$reflum) <= 0.1) {
		    $ntotobs++;
		}
	    }
	} else {
	    if ($startype[$i] == 0 || $startype[$i] == 1) {
		$ntotobs++;
	    }
	}
	$ntot++;
	$mcount += $m[$i];
	$i++;
    }
    printf("10%% Lag. rad.: f_b=%g+/-%g f_b,obs=%g+/-%g\n", $nbin/$ntot, sqrt($nbin)/$ntot, $nbinobs/$ntotobs, sqrt($nbinobs)/$ntotobs);

    $i = 0;
    $mcount = 0.0;
    $nbin = 0;
    $ntot = 0;
    $nbinobs = 0;
    $ntotobs = 0;
    while ($r[$i] <= $rmin) {
	$i++;
    }
    while ($i < $n && $r[$i] < $rmax) {
	if ($binflag[$i] == 1) {
	    $nbin++;
	    # a binary with a large enough mass ratio gets counted as a binary
	    if (($binstartype0[$i] == 0 || $binstartype0[$i] == 1) && ($binstartype1[$i] == 0 || $binstartype1[$i] == 1)) {
		$ntotobs++;
		$q = mymin($binm0[$i], $binm1[$i])/mymax($binm0[$i], $binm1[$i]);
		if ($q >= $qcrit) {
		    $nbinobs++;
		}
	    }
	    # a MS-compact object binary gets counted as a single MS object if it's not too far away from the MS
	    if ( (($binstartype0[$i] == 0 || $binstartype0[$i] == 1) && $binstartype1[$i] >= 10) ||
		 (($binstartype1[$i] == 0 || $binstartype1[$i] == 1) && $binstartype0[$i] >= 10) ) {
		if ($binstarrad0[$i] != 0.0) {
		    $temp0 = ($binstarlum0[$i]*$LSUN/(4.0*$PI*($binstarrad0[$i]*$RSUN)**2*$SIGMA))**0.25;
		} else {
		    $temp0 = 0.0;
		}
		$lum0 = $binstarlum0[$i];
		if ($binstarrad1[$i] != 0.0) {
		    $temp1 = ($binstarlum1[$i]*$LSUN/(4.0*$PI*($binstarrad1[$i]*$RSUN)**2*$SIGMA))**0.25;
		} else {
		    $temp1 = 0.0;
		}
		$lum1 = $binstarlum1[$i];
		$temp = ($lum0 * $temp0 + $lum1 * $temp1) / ($lum0 + $lum1);
		$lum = $lum0 + $lum1;
		if ($binstartype0[$i] == 0 || $binstartype0[$i] == 1) {
		    $reftemp = $temp0;
		    $reflum = $lum0;
		} else {
		    $reftemp = $temp1;
		    $reflum = $lum1;
		}
		if (2.0*abs($temp-$reftemp)/($temp+$reftemp) <= 0.1 &&
		    2.0*abs($lum-$reflum)/($lum+$reflum) <= 0.1) {
		    $ntotobs++;
		}
	    }
	} else {
	    if ($startype[$i] == 0 || $startype[$i] == 1) {
		$ntotobs++;
	    }
	}
	$ntot++;
	$i++;
    }
    printf("r_min<r<r_max: f_b=%g+/-%g f_b,obs=%g+/-%g\n", $nbin/$ntot, sqrt($nbin)/$ntot, $nbinobs/$ntotobs, sqrt($nbinobs)/$ntotobs);
}

# extract main sequence stars
sub extractmss {
    # wrong number of arguments?
    if ($#ARGV+1 != 4) {
	die("$usage");
    }
    
    $Mmin = $ARGV[2];
    $Mmax = $ARGV[3];

    $i=0;
    while ($i < $n) {
	if ($binflag[$i] == 0 && ($startype[$i] == 0 || $startype[$i] == 1)) {
	    if ($m[$i] >= $Mmin && $m[$i] <= $Mmax) {
		printf("%.8g %.8g %.8g %.8g\n", $r[$i], $rproj[$i], $m[$i], $L[$i]);
	    }
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
} elsif ($commandname =~ /^extractobsbinfrac$/) {
    extractobsbinfrac();
} elsif ($commandname =~ /^extractmss$/) {
    extractmss();
} else {
    die("$usage");
}
