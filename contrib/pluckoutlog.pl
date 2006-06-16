#!/usr/bin/perl -w

use strict;

# declare variables and subroutines
my ($usage, $PI, $trh0set, $firstline, $line);
my ($tcount, $t, $dt);
my ($e, $max_r, $nbound, $rtidal);
my ($m, $pe, $ke, $vratio);
my ($tidalmassloss);
my ($rc, $rhoc, $v_core, $trc, $conc, $Ncore);
my ($trh, $trh0, $rh, $rh_single, $rh_binary);
my ($nb, $mb, $eb);

sub grabval {
    my $key = shift(@_);
    my $line = shift(@_);
    $line =~ /$key[\s]*=[\s]*([0-9.\-+eE]+)/;
    return($1);
}

# define the usage
$usage = "Usage: $0\nPlucks data from the out.log file on stdin and sends it to stdout.\n";

# give usage for wrong number of arguments
if ($#ARGV != -1) {
    die("$usage");
}

# constants
$PI = 3.141592653589793116;

# have we already set these initial values?
$trh0set = 0;
$firstline = 0;

# and go...
while ($line = <STDIN>) {
    if ($line =~ /^[\s]*$/) {
	# do nothing for empty lines
    } elsif ($line =~ /^\*+$/) {
	# do nothing for these lines
    } elsif ($line =~ /^[\s]*tcount/) {
	$tcount = grabval("tcount", $line);
	$t = grabval("TotalTime", $line);
	$dt = grabval("Dt", $line);
    } elsif ($line =~ /^[\s]*sub\.count/) {
	# do nothing for these lines
    } elsif ($line =~ /^[\s]*Etotal/) {
	$e = grabval("Etotal", $line);
	$max_r = grabval("max_r", $line);
	$nbound = grabval("N_bound", $line);
	$rtidal = grabval("Rtidal", $line);
    } elsif ($line =~ /^[\s]*Mtotal/) {
	$m = grabval("Mtotal", $line);
	$pe = grabval("Etotal\.P", $line);
	$ke = grabval("Etotal\.K", $line);
	$vratio = grabval("VRatio", $line);
    } elsif ($line =~ /^[\s]*TidalMassLoss/) {
	$tidalmassloss = grabval("TidalMassLoss", $line);
    } elsif ($line =~ /^[\s]*core_radius/) {
	$rc = grabval("core_radius", $line);
	$rhoc = grabval("rho_core", $line);
	$v_core = grabval("v_core", $line);
	$trc = grabval("Trc", $line);
	$conc = grabval("conc_param", $line);
	$Ncore = grabval("N_core", $line);
    } elsif ($line =~ /^[\s]*trh/) {
	$trh = grabval("trh", $line);
	if ($trh0set == 0) {
	    $trh0 = $trh;
	    $trh0set = 1;
	}
	$rh = grabval("rh", $line);
	$rh_single = grabval("rh_single", $line);
	$rh_binary = grabval("rh_binary", $line);
    } elsif ($line =~ /^[\s]*N_b/) {
	$nb = grabval("N_b", $line);
	$mb = grabval("M_b", $line);
	$eb = grabval("E_b", $line);
	
	if ($rc != 0) {
	    # this is the last line output per timestep, so print now
	    if ($firstline == 0) { # print header
		print("#1:t/trh0 2:M 3:nbound 4:E 5:PE 6:KE 7:VRatio 8:rhoc 9:v_core 10:trc 11:trh 12:conc 13:ncore 14:rc 15:rh 16:rh_single 17:rh_binary 18:nb 19:mb 20:eb\n");
		$firstline = 1;
	    }
	    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
		   $t/$trh0, #1
		   
		   $m, #2
		   $nbound, #3
		   
		   $e, #4
		   $pe, #5
		   $ke, #6
		   $vratio, #7
		   
		   $rhoc, #8
		   $v_core, #9
		   $trc, #10
		   $trh, #11
		   $conc, #12
		   $Ncore, #13
		   
		   $rc, #14
		   $rh, #15
		   $rh_single, #16
		   $rh_binary, #17
		   
		   $nb, #18
		   $mb, #19
		   $eb); #20
	}
    } elsif ($line =~ /^[\s]*perturb_stars/) {
	# do nothing for these lines
    } elsif ($line =~ /^[\s]*sniff_stars/) {
	# do nothing for these lines
    } elsif ($line =~ /^[\s]*get_positions/) {
	# do nothing for these lines
    } elsif ($line =~ /^[\s]*dynamics_apply/) {
	# do nothing for these lines
    } else {
	warn("warning: unmatched line: $line");
    }
}
