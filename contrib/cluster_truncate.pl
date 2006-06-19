#!/usr/bin/perl
#
# This script truncates output files (made with cmc) after a certain input time 

# Usage: 

if($#ARGV != 1){
die("Usage: $0 <file prefix> <stop time>\n e.g., $0 w11_n1e5_fb0.1.out 0.43\n");
}

$file_prefix = $ARGV[0];
$stoptime = $ARGV[1];

# Store all suffixes in an array  

@file_suffixes = (".avemass_lagrad.dat",".centmass.dat", ".dyn.dat", ".ke_rad_lagrad.dat", ".ke_tan_lagrad.dat", ".lagrad0-0.1-1.dat", ".lagrad1-1-10.dat", ".lagrad2-10-100.dat", ".lagrad3-100-1000.dat", ".lagrad.dat", ".nostar_lagrad.dat", ".relaxation.dat", ".rho_lagrad.dat", ".v2_rad_lagrad.dat", ".v2_tan_lagrad.dat");  

# Shift loop - process each array element in turn

for($i=1;$i<=15;$i++){
$file_suffix = shift(@file_suffixes); 
&truncate;
}

# Upon end...

print "Complete. \n";


# Generic truncation subroutine

sub truncate
{

print ("Truncating $file_prefix"."$file_suffix . . .\n");
  $kept  = 0;
  $total = 0;

# Check successful open

if (!open(FP, $file_prefix.$file_suffix)){
    print "File $file_prefix"."$file_suffix not found. . . skipping . . .\n";
    end;
}  else {

# If so, proceed

    rename("$file_prefix$file_suffix", "$file_prefix$file_suffix.bak");
    open(OUT, ">$file_prefix$file_suffix");

# Read in line by line

    while($lines = <FP>){
	$total++;    
	$_ = $lines;

# Check if first entry is a number (with or without  exponent) If so, set to time
   
    if (/^([0-9]*\.*[0-9]+([eE][-+][0-9]+)*)/ == 1){
        $time = $1;
}

# Check if time is within inputted time constraint. If so, print to file; otherwise do nothing
  
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

}
