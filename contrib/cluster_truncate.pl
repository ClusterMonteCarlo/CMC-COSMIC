#!/usr/bin/perl

if($#ARGV != 1){
die("Usage: $0 <file prefix> <stop time>\n e.g., $0 w11_n1e5_fb0.1.out 0.43\n");
}

$file_prefix = $ARGV[0];
$stoptime = $ARGV[1];

# truncate the .avemass_lagrad.dat file
print ("Truncating $file_prefix".".avemass_lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".avemass_lagrad.dat")){
    print "File $file_prefix".".avemass_lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".avemass_lagrad.dat", "$file_prefix".".avemass_lagrad.dat.bak");
    open(OUT, ">$file_prefix.avemass_lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
        $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .centmass.dat file
print ("Truncating $file_prefix".".centmass.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".centmass.dat")){
    print "File $file_prefix".".centmass.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".centmass.dat", "$file_prefix".".centmass.dat.bak");
    open(OUT, ">$file_prefix.centmass.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
	    $time = $1;
	}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .dyn.dat file
print ("Truncating $file_prefix".".dyn.dat . . .\n");
   $kept  = 0;
   $total = 0;
if (!open(FP, $file_prefix.".dyn.dat")){
    print "File $file_prefix".".dyn.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".dyn.dat", "$file_prefix".".dyn.dat.bak");
    open(OUT, ">$file_prefix.dyn.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+)/ == 1){
	    $time = $1;
	}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .ke_rad_lagrad.dat file
print ("Truncating $file_prefix".".ke_rad_lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".ke_rad_lagrad.dat")){
    print "File $file_prefix".".ke_rad_lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".ke_rad_lagrad.dat", "$file_prefix".".ke_rad_lagrad.dat.bak");
    open(OUT, ">$file_prefix.ke_rad_lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
        $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .ke_tan_lagrad.dat file
print ("Truncating $file_prefix".".ke_tan_lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".ke_tan_lagrad.dat")){
    print "File $file_prefix".".ke_tan_lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".ke_tan_lagrad.dat", "$file_prefix".".ke_tan_lagrad.dat.bak");
    open(OUT, ">$file_prefix.ke_tan_lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
        $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .lagrad0-0.1-1.dat file
print ("Truncating $file_prefix".".lagrad0-0.1-1.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".lagrad0-0.1-1.dat")){
    print "File $file_prefix".".lagrad0-0.1-1.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".lagrad0-0.1-1.dat", "$file_prefix".".lagrad0-0.1-1.dat.bak");
    open(OUT, ">$file_prefix.lagrad0-0.1-1.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
	    $time = $1;
	}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .lagrad1-1-10.dat file
print ("Truncating $file_prefix".".lagrad1-1-10.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".lagrad1-1-10.dat")){
    print "File $file_prefix".".lagrad1-1-10.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".lagrad1-1-10.dat", "$file_prefix".".lagrad1-1-10.dat.bak");
    open(OUT, ">$file_prefix.lagrad1-1-10.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
    $time = $1;
}
if($time <= $stoptime){
            print OUT "$lines";
    $kept++; 
}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .lagrad2-10-100.dat file
print ("Truncating $file_prefix".".lagrad2-10-100.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".lagrad2-10-100.dat")){
    print "File $file_prefix".".lagrad2-10-100.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".lagrad2-10-100.dat", "$file_prefix".".lagrad2-10-100.dat.bak");
    open(OUT, ">$file_prefix.lagrad2-10-100.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
    $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .lagrad3-100-1000.dat file
print ("Truncating $file_prefix".".lagrad3-100-1000.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".lagrad3-100-1000.dat")){
    print "File $file_prefix".".lagrad3-100-1000.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".lagrad3-100-1000.dat", "$file_prefix".".lagrad3-100-1000.dat.bak");
    open(OUT, ">$file_prefix.lagrad3-100-1000.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
	    $time = $1;
	}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .lagrad.dat file
print ("Truncating $file_prefix".".lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".lagrad.dat")){
    print "File $file_prefix".".lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".lagrad.dat", "$file_prefix".".lagrad.dat.bak");
    open(OUT, ">$file_prefix.lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
    $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .nostar_lagrad.dat file
print ("Truncating $file_prefix".".nostar_lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".nostar_lagrad.dat")){
    print "File $file_prefix".".nostar_lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".nostar_lagrad.dat", "$file_prefix".".nostar_lagrad.dat.bak");
    open(OUT, ">$file_prefix.nostar_lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
	    $time = $1;
	}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .relaxation.dat file
print ("Truncating $file_prefix".".relaxation.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".relaxation.dat")){
    print "File $file_prefix".".relaxation.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".relaxation.dat", "$file_prefix".".relaxation.dat.bak");
    open(OUT, ">$file_prefix.relaxation.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+)/ == 1){
    $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .rho_lagrad.dat file
print ("Truncating $file_prefix".".rho_lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".rho_lagrad.dat")){
    print "File $file_prefix".".rho_lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".rho_lagrad.dat", "$file_prefix".".rho_lagrad.dat.bak");
    open(OUT, ">$file_prefix.rho_lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
    $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .v2_rad_lagrad.dat file
print ("Truncating $file_prefix".".v2_rad_lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".v2_rad_lagrad.dat")){
    print "File $file_prefix".".v2_rad_lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".v2_rad_lagrad.dat", "$file_prefix".".v2_rad_lagrad.dat.bak");
    open(OUT, ">$file_prefix.v2_rad_lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
	    $time = $1;
	}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

# truncate the .v2_tan_lagrad.dat file
print ("Truncating $file_prefix".".v2_tan_lagrad.dat . . .\n");
  $kept  = 0;
  $total = 0;
if (!open(FP, $file_prefix.".v2_tan_lagrad.dat")){
    print "File $file_prefix".".v2_tan_lagrad.dat not found. . . skipping . . .\n";
    end;
}  else {
    rename("$file_prefix".".v2_tan_lagrad.dat", "$file_prefix".".v2_tan_lagrad.dat.bak");
    open(OUT, ">$file_prefix.v2_tan_lagrad.dat");
    while($lines = <FP>){
	$total++;    
	$_ = $lines;
	if (/^([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)/ == 1){
    $time = $1;
}
	if($time <= $stoptime){
            print OUT "$lines";
	    $kept++; 
	}
    }
    print "Done. (Kept $kept / $total lines) \n"; 
}

print "Complete.\n";
