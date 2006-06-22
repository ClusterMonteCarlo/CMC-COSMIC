#!/usr/bin/ruby
#
# This script truncates various files associated with the cmc output with a user inputted file prefix and stop time

if ARGV.length != 2
print "Usage: ./truncate.rb <file_prefix> <stop time>\n i.e. ./truncate.rb w11_n1e5_fb0.1 0.75\n"
exit
end

# set variables

prefix = ARGV[0]
stoptime = ARGV[1].to_f

suffixes = [".avemass_lagrad.dat", ".centmass.dat", ".dyn.dat", ".ke_rad_lagrad.dat", ".ke_tan_lagrad.dat", ".lagrad0-0.1-1.dat", ".lagrad1-1-10.dat",".lagrad2-10-100.dat", ".lagrad3-100-1000.dat", ".lagrad.dat", ".nostar_lagrad.dat", ".relaxation.dat", ".rho_lagrad.dat", ".v2_rad_lagrad.dat",    ".v2_tan_lagrad.dat"]

# truncation loop

15.times {
filename = prefix + suffixes.shift
puts 'Truncating ' + filename + '...'
kept = 0
begin
	File.rename(filename, filename + ".bak")
	write = File.new(filename, "w")
	File.open(filename + ".bak").each { |line|
		if line =~ /^([0-9.]+[eE]*[-+]*[0-9]*)/
			if $&.to_f <= stoptime
				printf(write, line)
				kept = kept + 1
			       	
			else
				break
				end
			else	
				printf(write, line)
			end

}
print "Done. (Kept " + kept.to_s + " lines)\n"
rescue SystemCallError
       print "Cannot find file '" + filename + "', skipping...\n"
end
}
print "Complete.\n"