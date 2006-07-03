#!/usr/bin/ruby

require 'zlib'
include Math

if ARGV.length != 3
print "Usage: ./3d_density.rb <file_prefix> <snapshot #> <half mass radius>\n e.g. ./3d_density.rb w6_n1e5_fb0.03.out 0076 0.804557 \n"
exit
end

prefix = ARGV[0]
snap_num = ARGV[1]
r_h = ARGV[2].to_f

mass = []
lum = []
rad = []

# Read in file and assign relative information to arrays
begin
input = open(prefix + "." + "snap" + snap_num + ".dat.gz")
rescue Errno::ENOENT
  puts "File not found, aborting..."
exit
end

gz = Zlib::GzipReader.new(input)

gz.each_line { |line|
	if line =~ /\A[0-9]+/
	    	info_array = [ ]
	   	info_array = line.split(" ")
   		mass.push(info_array[1].to_f)
		lum.push(info_array[1].to_f**3)
		rad.push(info_array[2].to_f/r_h)
	else 	
	end
}

n_stars = mass.length - 1
bin_width = 10
i = 0

write = File.new(prefix + "." + "snap" + snap_num + ".3d", "w")

n_stars.times {
	
	r = rad[i] 
	i_m = i - bin_width
	
	if i_m < 0
		i_m = 0
	end

	i_M = i + bin_width
	
	if i_M > n_stars
		i_M = n_stars
	end

	rad_m = rad[i_m]
	rad_M = rad[i_M]

	tot_mass = mass[i_m]/2.0 + mass[i_M]/2.0
	tot_lum  = lum[i_m]/2.0 + lum[i_M]/2.0
	
	((i_m+1)..(i_M-1)).each { |x|
		tot_mass += mass[x]
		tot_lum += lum[x]
	}

	rho = tot_mass/(4.0/3.0 * 3.14159265358979312 * (rad_M**3 - rad_m**3))
	
	printf(write, "%f %f\n", r, rho) 
	
	i = i + 1
}

puts "Done."
