#!/usr/bin/ruby

# This script can be used to recover the original, non-truncated files stored at *.bak, located from a inputted file prefix  

if ARGV.length != 1
print "Usage: ./recover_truncated <file_prefix>
      e.g. ./recover_truncated w11_1e5_fb0.1\n"
exit
end

prefix = ARGV[0]


suffixes = [".avemass_lagrad.dat", ".centmass.dat", ".dyn.dat", ".ke_rad_lagrad.dat", ".ke_tan_lagrad.dat", ".lagrad0-0.1-1.dat", ".lagrad1-1-10.dat"\
,".lagrad2-10-100.dat", ".lagrad3-100-1000.dat", ".lagrad.dat", ".nostar_lagrad.dat", ".relaxation.dat", ".rho_lagrad.dat", ".v2_rad_lagrad.dat",   \
".v2_tan_lagrad.dat"]

15.times {
begin
filename = prefix + suffixes.shift
File.rename(filename + '.bak', filename)
rescue SystemCallError
       print "Cannot find " + filename + ".bak...\n"
end
}
print "Done.\n"
 