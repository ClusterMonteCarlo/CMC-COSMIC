#!/usr/bin/python2.2

# a simple code to average a given data file
# usage: ./smoother.py [-s<sigma>] [<infile> [<outfile>]]
# hopping not implemented yet

import sys, string, getopt
from math import exp

# Defaults
defa_sigma = 20
defa_hop_length = 1
defa_hop_offset = 0

# === begin parsing and setting options ===
args = '-ssig -hhop -ooff'
optlist, args = getopt.getopt(sys.argv[1:], 's:h:o:', '')
#print optlist
#print args, len(args)
#print '******************'

def dete_value(optlist, optiden, optname, default, lowlim) :
	for opt in optlist :
		if opt[0] == optiden :
			value = int(opt[1])
			if value < lowlim :
#				print optname, 'has to be at least', lowlim
#				print 'using default value', default
				value = default
			break
	else :
#		print 'using default value', default, 'for', optname
		value = default
	return value

sigma = dete_value(optlist, '-s', 'sigma', defa_sigma, 1)
hop_length = dete_value(optlist, '-h', 'hop length', defa_hop_length, 1)
hop_offset = dete_value(optlist, '-o', 'hop offset', defa_hop_offset, 0)

#print 'sigma is', sigma
#print 'hop length is', hop_length
#print 'hop offset is', hop_offset

# Determining the input and output file, if any
if len(args) == 0 :
	instream = sys.stdin
	# print 'reading from stdin'
	outstream = sys.stdout
	# print 'writing to stdout'
elif len(args) == 1: 
	if args[0] == '--' :
		instream = sys.stdin
		# print 'reading from stdin'
	else :
		instream = file(args[0])
	outstream = sys.stdout
	# print 'writing to stdout'
elif len(args) == 2:
	if args[0] == '--' :
		instream = sys.stdin
		# print 'reading from stdin'
	else :
		instream = file(args[0])
	if args[1] == '--' :
		outstream = sys.stdout
		# print 'writing to stdout'
	else :
		outstream = file(args[1], 'w')
else :
	print 'illegal number of arguments: ', len(args)
	sys.exit(1)
# === end parsing and setting options ===

rows = [line.split() for line in instream]

# create the list 'goodrows' while 
# 1. discarding empty lines and lines starting with '#'
# 2. checking if the other lines have same number of fields
goodrows = []
no_field = 0
for i in range(len(rows)):
	if (rows[i] != []) and (rows[i][0][0] != '#') :
		goodrows.append(rows[i])
		if no_field == 0:
			no_field = len(rows[i])
		else :
			if no_field != len(rows[i]):
				print 'number of fields in line', i,
				print 'does not match previous lines'
				sys.exit(2)

if no_field == 0 :
	print 'no lines found!'
	sys.exit(3)
#elif no_field > 1 :
#	print 'handling multiple columns in data files is not implemented yet'
#	sys.exit(4)

norows = len(goodrows)
nocols = len(goodrows[0])

# gooddata is transpose of goodrows
gooddata = [[r[col] for r in goodrows] for col in range(len(goodrows[0]))]

def calc_aver(data, first, point, last) :
	sum = 0
	wsum = 0
	for i in range(first, last) :
		sum += float(data[i]) * exp(-(i-point)**2/(sigma**2))
		wsum += exp(-(i-point)**2/sigma**2) 
	return sum/wsum

# main loop of the program, calculate and print averages
for i in range(hop_offset, norows, hop_length) :
	if i<2*sigma :		# we are too close to beginning
		first = 0;
	else:
		first = i-2*sigma

	if i+2*sigma+1>norows :	# we are too close to the end
		last = norows
	else:
		last = i+2*sigma+1

	if no_field==1 :	
		aver = calc_aver(gooddata[0], first, i, last)
		print >>outstream, aver
		#print >>outstream, gooddata[0][i]
	else :
		print >>outstream, gooddata[0][i],
		for j in range(1,nocols) :
			aver = calc_aver(gooddata[j], first, i, last)
			print >>outstream, aver,
		print >>outstream
