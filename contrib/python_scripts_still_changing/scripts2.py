#!/opt/local/bin//python
from numpy import *
import gzip
#scripts needed for extracting BSE merger and disruption history

def check_case(line) : # check different merger cases
	"""checks and sorts different merger cases"""
    	if line.rfind("disrupt1")>-1 : case = 1
    	elif line.rfind("disrupt2")>-1 : case = 2
    	elif line.rfind("disruptboth")>-1 : case = 3
    	else : raise UserWarning("this line does not have legal merger flag: %s\n" % line)
    	return "%d " % case

def readmergefile(filename) :
	"""reads the full merger file one line at a time and stores it in a convenient manner"""
	print 'input file:', filename
	return [line.replace("t=","")
		.replace("disrupt1 ",check_case(line))
		.replace("disrupt2 ",check_case(line))
		.replace("disruptboth ",check_case(line)+' -100 -100 ')
		.replace("idr=","").replace("(mr="," ")
		.replace(") id1="," ").replace("(m1="," ")
		.replace("):id2="," ").replace("(m2="," ")
		.replace(" id1="," ").replace(" id2="," ")
		.replace(") (r="," ").replace(")","")
	for line in file(filename) if line[0]!="#"]

def bse_int(filename):
	"""reads the full merger file and creates the merger dictionary"""
	rows = readmergefile(filename)
	#print rows
	se_dict = {}
	#print se_dict
	for i in range(len(rows)):
		each_row = rows[i].split()
		#print each_row
		#print (i,len(rows),each_row[1])
		if int(each_row[1])<3:
			ID = long(each_row[2])
			case = int(each_row[1])
			parents = {'nopar': 2,
				'IDs': [long(each_row[4+i*2]) for i in range(2)],
				'masses': [float(each_row[5+i*2]) for i in range(2)]
				}
			se_dict[ID] = {'time': float(each_row[0]),
					'interaction': each_row[1],
					'id': long(each_row[2]),
					'mass': float(each_row[3]),
					'parents': parents,
					'position': float(each_row[3+2*2+1]),
					'fin_prod': 1,
					'description':'mergers with one star gone'
					}
					
		elif int(each_row[1])==3:
			ID1 = long(each_row[4])
			ID2 = long(each_row[6])
			parents = {}
			se_dict[ID1] = {'time': float(each_row[0]),
					'interaction': each_row[1],
					'id': long(each_row[4]),
					'mass': float(each_row[5]),
					'parents': parents,
					'position': float(each_row[3+2*2+1]),
					'fin_prod': 1,
					'description': 'both stars intact'
					}
			se_dict[ID2] = {'time': float(each_row[0]),
					'interaction': each_row[1],
					'id': long(each_row[6]),
					'mass': float(each_row[7]),
					'parents': parents,
					'position': float(each_row[3+2*2+1]),
					'fin_prod': 1,
					'description': 'both stars intact'
					}
		else:
			raise UserWarning("totally confused")
	
	return se_dict	
	
		
			
	

