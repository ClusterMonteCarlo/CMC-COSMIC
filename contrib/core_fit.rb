#!/usr/bin/ruby

include Math

# This program takes a .2d extracted density/luminosity file and fits the core (defined where the luminosity drops by half,
# corresponds to a magnitude drop of 0.7525275) to a logarithmic linear least squares line.

if ARGV.length != 2
puts "Usage: ./core_fit.rb <file prefix> <snapshot #>"
puts "e.g. ./core_fit.rb sp_w6_n1e5_fb0.03.out 0954"
exit
end

filename = ARGV[0]
snapnum = ARGV[1]

x = 0.0
y = 0.0
xx = 0.0
xy = 0.0

checkflag = 0
corelimit = 0.0
xstore = []
ystore = []
n = 0
filename = filename + '.snap' + snapnum

write = File.new(filename + '.core', "w")
info = []

File.open(filename + '.2d', "r").each { |line|

  info = line.split(" ")
  if checkflag == 0
    puts "central mag = " + info[2] + " mag/arcsec^2"
    corelimit = info[2].to_f + 0.7525275
    puts "core limit = " + corelimit.to_s + " mag/arcsec^2"
    checkflag = 1
  end
  if info[2].to_f <= corelimit
    xstore.push(info[0].to_f)
    ystore.push(info[2].to_f)
    x +=  log10(info[0].to_f)
    y +=  (info[2].to_f)
    xy += (log10(info[0].to_f) * (info[2].to_f))
    xx += (log10(info[0].to_f) * log10(info[0].to_f))
    n += 1
  else
    puts "core radius = " + info[0] + " r_h"
    break
    end
}

  # Calculate least squares linear fit
  x /= n
  y /= n
  xx /= n
  xy /= n

  a = ((xy - x*y) / (xx - x*x))
  b = (xx*y - x*xy) / (xx - x*x)

  puts "\nLEAST SQUARES INFO (Y = AX + B)"
  puts "a = " + a.to_s + " b = " + b.to_s

xstore.each_index { |xval|
  yval = ystore.shift
  fitval = (a)*log10(xstore[xval]) + b
  printf(write, "%f %f %f\n", xstore[xval], yval, fitval)
}
