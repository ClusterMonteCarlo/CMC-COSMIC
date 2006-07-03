#!/usr/bin/ruby

require 'zlib'
require '~/bin/D_compute'
include Math

if ARGV.length != 3
print "Usage: ./2d_density.rb <file_prefix> <snapshot #> <half mass radius>\n e.g. ./2d_density.rb w6_n1e5_fb0.03.out 0076 0.804557 \n"
exit
end
  
prefix = ARGV[0]
snap_num = ARGV[1]
r_h = ARGV[2].to_f


# Read in file and assign relative information to arrays                                                                                              
begin
input = open(prefix + "." + "snap" + snap_num + ".dat.gz")

# Handle errors                                                                                                                                       
rescue Errno::ENOENT
puts "File not found, aborting..."
exit
end

class D_compute
       attr_accessor :mass, :lum, :rad, :rho, :rhoL, :rRad_m, :rRad_M, :rRad_ratio
       def initialize(mass = [], lum = [], rad = [], rho = [], rhoL = [], rRad_m = 0, rRad_M = 0, rRad_ratio = 0)
       	   @mass = mass
	   @lum = lum
	   @rad = rad
	   @rho = rho
	   @rhoL = rhoL
	   @rRad_m = rRad_m
	   @rRad_M = rRad_M
	   @rRad_ratio = rRad_ratio
       end
end
compute = D_compute.new

gz = Zlib::GzipReader.new(input)

gz.each_line { |line|
        if line =~ /\A[0-9]+/
               info_array = [ ]
               info_array = line.split(" ")
               compute.mass.push(info_array[1].to_f)
          if info_array[7].to_i == 0
               compute.lum.push(info_array[1].to_f**3)
          elsif info_array[7].to_i == 1
               compute.lum.push((info_array[8].to_f**3)+(info_array[9].to_f**3))
          else
                puts "Binaries mislabeled."
          end
               compute.rad.push(info_array[2].to_f/r_h)
        else
        end
}

n_stars = compute.mass.length - 1
bin_width = 10
i = 0

filename = (prefix + "." + "snap" + snap_num + ".2d")
write = File.new(filename, "w")

n_stars.times {

        r = compute.rad[i]
        i_m = i - bin_width

        if i_m < 0
                i_m = 0
        end

        i_M = i + bin_width

        if i_M > n_stars
                i_M = n_stars
        end

        rad_m = compute.rad[i_m]
        rad_M = compute.rad[i_M]

        tot_mass = compute.mass[i_m]/2.0 + compute.mass[i_M]/2.0
        tot_lum  = compute.lum[i_m]/2.0 + compute.lum[i_M]/2.0

  ((i_m+1)..(i_M-1)).each { |x|
                tot_mass += compute.mass[x]
                tot_lum += compute.lum[x]
  }

        compute.rho[i] = tot_mass/(4.0/3.0 * 3.14159265358979312 * (rad_M**3 - rad_m**3))
        compute.rhoL[i] = tot_lum/(4.0/3.0 * 3.14159265358979312 * (rad_M**3 - rad_m**3))
        i = i + 1
}

p = 1000
compute.rRad_m = 1.0e-4
compute.rRad_M = compute.rad[n_stars]

compute.rRad_ratio = compute.rRad_M/compute.rRad_m
compute.do(n_stars, p)
File.rename("temp.dat", filename) 
puts "Done."
