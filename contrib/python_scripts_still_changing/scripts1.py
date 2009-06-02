#!/opt/local/bin//python
from numpy import *
import gzip

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
		    .replace("idm=",cnp(line))
		    .replace("idm=",cnp(line))
		    .replace("idm=",cnp(line))
		    .replace("(mm="," ")
		    .replace(") id1="," ").replace("(m1="," ")
		    .replace("):id2="," ").replace("(m2="," ")
		    .replace("):id3="," ").replace("(m3="," ")
		    .replace("):id4="," ").replace("(m4="," ")
		    .replace(") (r="," ")
		    .replace(")\n","") 
		for line in file(filename) if line[0]!="#"]

#rows=readcollfile('coll.log')


def collision(collfile): 
	rows = readcollfile(collfile)
	#print rows
	coll_prod = {}
	for i in range(len(rows)):
		each_row = rows[i].split()
		ID = long(each_row[3])
		nopars = int(each_row[2])
		parents = {'nopar': nopars,
			'IDs': [long(each_row[5+i*2]) for i in range(nopars)],
			'masses': [float(each_row[6+i*2]) for i in range(nopars)]
			}
		coll_prod[ID] = {'time': float(each_row[0]),
				'interaction': each_row[1],
				'coll_id': long(each_row[3]),
				'mass': float(each_row[4]),
				'parents': parents,
				'position': float(each_row[4+nopars*2+1]),
				'fin_prod': 1
				}
		for parid in parents['IDs']:
			parid
			if coll_prod.has_key(parid):
				coll_prod[parid]['fin_prod'] = 0
	return coll_prod



def id_tracker_coll(coll,id,count,coll_history):
	"""takes: 
	coll: collision dictionary made by the collision module, 
	id: star id in question, 
	count: counter for the level of branching of the collision tree
	coll_history: another initialized dictionary where the tree is written out
	the reason for the need of collision_history is to be able to call this as a module externally from another module"""
	try:
		#print id
		if coll.has_key(id):
			coll_history[count]=coll[id]
			#coll_history
		else:
			#print 'star %ld did not have a collision' % (id,)
			raise StopIteration()
	
		#print coll_history[count]['parents']['nopar']
		for i in range(coll_history[count]['parents']['nopar']):
			parid = coll_history[count]['parents']['IDs'][i]
			print parid, count
			if coll.has_key(parid):
				#print 'parid found'
				count += 1
				id_tracker_coll(coll,parid,count,coll_history)
			else:
				#print 'not found'
				raise StopIteration()
	except StopIteration:
		pass
	return coll_history



def call_collision_tree(coll,id):
	"""This one basically calls the tree maker above
	collfile: the collisionfile is fed to collision routine to populate coll dictionary
	id: id in question
	coll_history: initialize the dictionary where the tree will be written out"""
	#coll=collision(collfile)
	collision_history={}
	count=0 # this is the first level in the tree so always zero, recursion is in id_tracker_coll
	coll_tree = id_tracker_coll(coll,id,count,collision_history)
	return coll_tree





