#!/usr/bin/perl
if($#ARGV != 1){
    die("Usage: $0 <file prefix> <stop time>\nfor example, $0 cluster_out 0.78\n");
}

$prefix = $ARGV[0];
$stoptime = $ARGV[1];

# crop the _0 file
print("Cropping $prefix"."_0 . . .\n");
open(FP, "$prefix"."_0");
@lines = <FP>;
close(FP);
rename("$prefix"."_0", "$prefix"."_0.uncropped");

open(FP, ">$prefix"."_0");
$_ = shift(@lines);
while(!/^>>>Tcount +=/){
    printf(FP);
    $_ = shift(@lines);
}

/^>>>.+ +Time += +([0-9.]+)/;
$time = $1;
while($time <= $stoptime){
    printf(FP);
    $_ = shift(@lines);
    if(/^>>>.+ +Time += +([0-9.]+)/){
	$time = $1;
    }
}

close(FP);
print(" . . . done\n");

# crop the _1 file
print("Cropping $prefix"."_1 . . .\n");
open(FP, "$prefix"."_1");
@lines = <FP>;
close(FP);
rename("$prefix"."_1", "$prefix"."_1.uncropped");

open(FP, ">$prefix"."_1");
$_ = shift(@lines);
/^[0-9]+ +([0-9.]+)/;
$time = $1;
while($time <= $stoptime){
    printf(FP);
    $_ = shift(@lines);
    /^[0-9]+ +([0-9.]+)/;
    $time = $1;
}

close(FP);
print(" . . . done\n");

#crop the _3 file
print("Cropping $prefix"."_3 . . .\n");
open(FP, "$prefix"."_3");
@lines = <FP>;
close(FP);
rename("$prefix"."_3", "$prefix"."_3.uncropped");

open(FP, ">$prefix"."_3");
$_ = shift(@lines);
/^([0-9]+)/;
$time = $1;
while($time <= $stoptime){
    printf(FP);
    $_ = shift(@lines);
    /^([0-9.]+)/;
    $time = $1;
}

close(FP);
print(" . . . done\n");

#crop the _4 file
print("Cropping $prefix"."_4 . . .\n");
open(FP, "$prefix"."_4");
@lines = <FP>;
close(FP);
rename("$prefix"."_4", "$prefix"."_4.uncropped");

open(FP, ">$prefix"."_4");
$_ = shift(@lines);
/^[0-9]+ +([0-9.]+)/;
$time = $1;
while($time <= $stoptime){
    printf(FP);
    $_ = shift(@lines);
    /^[0-9]+ +([0-9.]+)/;
    $time = $1;
}

close(FP);
print(" . . . done\n");

#crop the _5 file
print("Cropping $prefix"."_5 . . .\n");
open(FP, "$prefix"."_5");
@lines = <FP>;
close(FP);
rename("$prefix"."_5", "$prefix"."_5.uncropped");

open(FP, ">$prefix"."_5");
$_ = shift(@lines);
/^([0-9]+)/;
$time = $1;
while($time <= $stoptime){
    printf(FP);
    $_ = shift(@lines);
    /^([0-9.]+)/;
    $time = $1;
}

close(FP);
print(" . . . done\n");

#crop the _esc file
print("Cropping $prefix"."_esc . . .\n");
open(FP, "$prefix"."_esc");
@lines = <FP>;
close(FP);
rename("$prefix"."_esc", "$prefix"."_esc.uncropped");

open(FP, ">$prefix"."_esc");
$_ = shift(@lines);
/^[0-9]+ +([0-9.]+)/;
$time = $1;
while($time <= $stoptime){
    printf(FP);
    $_ = shift(@lines);
    /^[0-9]+ +([0-9.]+)/;
    $time = $1;
}

close(FP);
print(" . . . done\n");

#crop the _stellar file
print("Cropping $prefix"."_stellar . . .\n");
open(FP, "$prefix"."_stellar");
@lines = <FP>;
close(FP);
rename("$prefix"."_stellar", "$prefix"."_stellar.uncropped");

open(FP, ">$prefix"."_stellar");
$_ = shift(@lines);
/^[0-9]+ +([0-9.]+)/;
$time = $1;
while($time <= $stoptime){
    printf(FP);
    $_ = shift(@lines);
    /^[0-9]+ +([0-9.]+)/;
    $time = $1;
}

close(FP);
print(" . . . done\n");

