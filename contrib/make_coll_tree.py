#!/usr/bin/env python2.4

import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--infile", dest="filename",
                  help="name of the input file", metavar="FILE")
parser.add_option("-o", "--outfile", dest="psfilename",
                  help="name of the output file", metavar="FILE",
			default="coll_tree")

(options, args) = parser.parse_args()

# The data file has two kinds of lines
# 1. comments, which start with "#", to be ignored
# 2. data, with the following format
# t=<float> <interaction type> idm=<long>(mm=<float>) id1=<long>(m1=<float>):id2=<long>(m2=<float>):...
# -> steps for parsing (to avoid using regular expressions <wink>):
# 1. make substitution to create space seperated strings,
# 2. split using space,
# 3. assign from result of splitting

# reading our files
def cnp(line) : # calculate number of parents for a given line
    if line.rfind("id4")>-1 : nops = 4
    elif line.rfind("id3")>-1 : nops = 3
    elif line.rfind("id2")>-1 : nops = 2
    else : raise UserWarning('the following line does not have legal '+
				     'parent IDs:\n %s' % line)
    return "%d " % nops

def readcollfile(filename) :
    print 'input file:', filename
    return [line.replace("t=","")
		    .replace("single-single idm=",cnp(line))
		    .replace("binary-binary idm=",cnp(line))
		    .replace("binary-single idm=",cnp(line))
		    .replace("(mm="," ")
		    .replace(") id1="," ").replace("(m1="," ")
		    .replace("):id2="," ").replace("(m2="," ")
		    .replace("):id3="," ").replace("(m3="," ")
		    .replace("):id4="," ").replace("(m4="," ")
		    .replace(")\n","") 
		for line in file(filename) if line[0]!="#"]

if options.filename==None and args == []: # nothing is given, 
						      # use default name
    rows = readcollfile('cmc_debug.collision.log')
elif options.filename != None : # if option is given, use it
    rows = readcollfile(options.filename)
else : # option is not given but there are arguments
       # which we use as input filenames
    rows = []
    for x in args : rows += readcollfile(x)

# filling our dictionary/associative array
# { ID : { 'time'     : time of collision/formation,
#          'mass'     : mass of the collision product
#          'parents'  : dictionary for containing information about parents
#                       'nopar'  : number of parents
#                       'IDs'    : list containing IDs of parents
#                       'masses' : list containing masses of parents
#          'fin_prod' : flag to see if this is a final product or went into
#                       another collision
#         }
# }
coll_prod = {}
for i in range(len(rows)):
	row_elems = rows[i].split()
	ID = long(row_elems[2])
	nopars = int(row_elems[1])
	parents = {'nopar': nopars,
		     'IDs':    [long(row_elems[4+i*2]) for i in range(nopars)],
		     'masses': [float(row_elems[5+i*2]) for i in range(nopars)]
		    }
	coll_prod[ID] = {'time': float(row_elems[0]), 
			     'mass': float(row_elems[3]), 
			     'parents' : parents,  
			     'fin_prod' : 1
			     }
	for parid in parents['IDs']:
		if coll_prod.has_key(parid):
			coll_prod[parid]['fin_prod'] = 0

# sort the final products with respect to their masses
fin_prods = [(k, coll_prod[k]['mass']) for k,v in coll_prod.iteritems() 
				if coll_prod[k]['fin_prod']==1]
fin_prods.sort(lambda x, y: -cmp(x[1],y[1]))
mxm = fin_prods[0][1]

# but determine the maximum time (for axis) by simple a max()
mxtime =  max([coll_prod[x[0]]['time'] for x in fin_prods])

# drawing the tree
from pyx import *
from pyx.graph import graphxy, data, axis
from pyx.graph.axis import painter, tick
from pyx.deco import earrow

from math import pow

p = painter.regular(basepathattrs=[earrow.normal])
g = graphxy(width=15, height=8,
		x2 = None, y2 = None,
		x = axis.linear(title='time (code units)', 
				    min = mxtime*(-0.025), 
				    max = mxtime*1.05,
				    painter=p),
		y = axis.logarithmic(title=r'mass $(M_\odot)$', 
					   min = 0.1, max = mxm*2.0,
					   painter=p, density=1.0)
	     )
def draw_star(time, mass, con_color) :
    x, y= g.pos(time, mass)
    rad = pow(mass/mxm, 1.0/3)*0.1
    p = path.circle(x, y, rad)
    g.stroke(p, [deco.filled([color.cmyk.Yellow]),con_color])

def connect_star(time1, mass1, time2, mass2, con_color) :
    x1, y1 = g.pos(time1, mass1)
    x2, y2 = g.pos(time2, mass2)
    p = path.line(x1, y1, x2, y2)
    g.stroke(p, [con_color])

def draw_tree(id, con_color):
    star = coll_prod[id]
    for parno in range(star['parents']['nopar']) : 
	  parid = star['parents']['IDs'][parno]
	  parmass = star['parents']['masses'][parno]
	  if coll_prod.has_key(parid) :
		connect_star(coll_prod[parid]['time'], parmass,
				 star['time'], star['mass'], con_color)
		draw_tree(parid,  con_color)
	  else :
		connect_star(star['time']*0.9, parmass,
				 star['time'], star['mass'], con_color)
		draw_star(star['time']*0.9, parmass, con_color)
# alternatives for the tails of collision products:			
#			connect_star(star['time'], parent['mass'],
#					 star['time'], star['mass'])
#			draw_star(star['time'], parent['mass'])
#			connect_star(0.0, parent['mass'],
#					 star['time'], star['mass'])
#			draw_star(0.0, parent['mass'])
    draw_star(star['time'], star['mass'], con_color)

for i in range (5, len(fin_prods)) :
#for i in range (5, 10) :
    draw_tree(fin_prods[i][0], color.cmyk.Cyan)
#draw_tree(fin_prods[4][0], color.cmyk.Periwinkle)
#draw_tree(fin_prods[3][0], color.rgb.blue)
#draw_tree(fin_prods[2][0], color.cmyk.Brown)
#draw_tree(fin_prods[1][0], color.cmyk.Tan)
#draw_tree(fin_prods[0][0], color.rgb.red)
for i in range (4,-1,-1) :
    draw_tree(fin_prods[i][0], color.palette.ReverseHue.getcolor(float(i)/9))

print "%ld collision products" % len(fin_prods)
for i in range(5,-1,-1) :
    print fin_prods[i]

g.writeEPSfile(options.psfilename)
