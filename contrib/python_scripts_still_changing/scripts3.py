#!/opt/local/bin//python
from numpy import *
import gzip

def find_positions(f):
	"""Takes binint file's file pointer to find the byte positions for interactions"""
	#f=open(s1,'r')
	seek_line=()
	try:
    		while f:
        		line=f.readline()
        		#print line
        		if line=='********************************************************************************\n':
            			pos=f.tell()
            			#print pos
            			seek_line+=(pos,)
        		elif line=='':
				posend=f.tell()
				seek_line+=(posend,)
            			raise StopIteration()
	except StopIteration:
    		pass
		#print 'EOF'
	#f.close()
	return seek_line


def read_segment(f,position):
	"""reads each segment of the binint file that are contained within the **************s"""
	goodrun=1
	try:
		na=-100
		f.seek(position)
		line=f.readline()
		parsed=line.replace("type=","").replace("BB","BB 4").replace("BS","BS 3").replace("t="," ")
		values=parsed.split()
		#print values
		type={'type': values[0],
			'nopars': int(values[1]),
			'time': float(values[2])
			}
		line=f.readline()
		parsed=line.replace("params: b=","").replace("v=","")
		values=parsed.split()
		params={'b': float(values[0]),
			'v': float(values[1])
			}
		input={}	
		for i in range(2):
			line=f.readline()
			parsed=line.replace("input: type=","")\
				.replace("binary","binary 2").replace("single","single 1")\
				.replace("m0=","").replace("m1=","").replace("m=","")\
				.replace("R0=","").replace("R1=","").replace("R=","")\
				.replace("Eint1=","").replace("Eint2=","").replace("Eint=","")\
				.replace("id0=","").replace("id1=","").replace("id=","")\
				.replace("a=","").replace("e=","")
			values=parsed.split()
			#print values
			if values[0]=='binary':
				input[i]={'type': values[0],
					'no': int(values[1]),
					'm': [float(values[2]),float(values[3])],
					'R': [float(values[4]),float(values[5])],
					'Eint': [float(values[6]),float(values[7])],
					'ids': [long(values[8]),long(values[9])],
					'a': float(values[10]),
					'e': float(values[11])
					}
			
			elif values[0]=='single':
				input[i]={'type': values[0],
					'no': int(values[1]),
					'm': [float(values[2])],
					'R': [float(values[3])],
					'Eint': [float(values[4])],
					'ids': [long(values[5])],
					'a': na,
					'e': na
					}	
		
			else:
				raise UserWarning('the following line does not have legal '+
					     'input:\n %s' % line)
			
		
		line=f.readline()
		parsed=line.replace("status: DE/E=","").replace("DE=","")\
				.replace("DL/L=","").replace("DL=","").replace("tcpu=","")
		values=parsed.split()
		#print values
		status={'DE/E': float(values[0]),
			'DE': float(values[1]),
			'DL/L': float(values[2]),
			'DL': float(values[3]),
			'tcpu': float(values[4])
			}
		
		line=f.readline()
		if line.rfind('stopped')>-1:
			raise StopIteration()

		#print line+' look here %ld ' % position
		pos1=line.rfind('nstar=')+6
		pos2=line.rfind('nobj=')+5
		pos3=line.rfind('(')
		pos4=line.rfind(')')
		#print line[pos1]+'and here'
		outcome={'nstar': int(line[pos1]),
			'nobj': int(line[pos2]),
			'outcome': line[pos2+4:pos3-1], 
			'objects': line[pos3+1:pos4],
			}
		output={}	
		for i in range(outcome['nobj']):
			line=f.readline()
			parsed=line.replace("output: type=","").replace("binary","binary 2").replace("single","single 1")\
				.replace("m0=","").replace("m1=","").replace("m=","")\
				.replace("R0=","").replace("R1=","").replace("R=","")\
				.replace("Eint1=","").replace("Eint2=","").replace("Eint=","")\
				.replace("id0=","").replace("id1=","").replace("id=","")\
				.replace("a=","").replace("e=","")\
				.replace("triple","triple 3")\
				.replace("min0=","").replace("min1=","").replace("min=","")\
				.replace("mout0=","").replace("mout1=","").replace("mout=","")\
				.replace("Rin0=","").replace("Rin1=","").replace("Rin=","")\
				.replace("Rout0=","").replace("Rout1=","").replace("Rout=","")\
				.replace("Eintin0=","").replace("Eintin1=","").replace("Eintin=","")\
				.replace("Eintout0=","").replace("Eintout1=","").replace("Eintout=","")\
				.replace("idin1=","").replace("idin2=","").replace("idin=","")\
				.replace("idout0=","").replace("idout1=","").replace("idout=","")\
				.replace("ain=","").replace("aout=","").replace("ein=","").replace("eout=","")
			values=parsed.split()
			#print values
			
			no=int(values[1])

			tmp={}
			merger=[]	
			for j in range(no):
				if values[2 + 3*no + j].rfind(':')>-1:
					tmp[j]=values[2 + 3*no + j].split(':')
					merger+=[1]
				else:
					merger+=[0]
			merge={'merge': merger,
				'ids': tmp
				}
			a=[]
			e=[]
			if no>1:
				for j in range(no-1):
						a+=[float(values[-2* (no-1) +j])]
						e+=[float(values[-(no-1) +j])]
					
			output[i]={'type': values[0],
				'no': no,
				'm': [float(values[2+ j]) for j in range(no)],
				'R': [float(values[2+ no + j]) for j in range(no)],
				'Eint': [float(values[2+ 2*no +j]) for j in range(no)],
				'ids': [values[2+ 3*no +j] for j in range(no)],
				'a': a,
				'e': e,
				'merge': merge 
				}
	except StopIteration:
		goodrun=0
		outcome={}
		output={}
		pass
	
	binint_each={'success': goodrun,
		'type': type,
		'params': params,
		'input': input,
		'status': status,
		'outcome': outcome,
		'output': output
		}
	
	return binint_each


def read_binint(filename):
	"""calls the above scripts and stores all binint information in a dictionary called binint"""
	f=open(filename,'r')
	print 'inputfile %s' %(filename)
	positions=find_positions(f)
	binint={}
	for i in range(len(positions)-1):
		#print positions[i]
		binint[i]=read_segment(f,positions[i])
	#print binint.keys(),len(binint.keys()),len(positions)
	
	return binint
	

def binary_interaction(binint,id):
	"""asks for id of a star and then stores the binint info/history for that id
	first need to call read_binint to populate the binint dictionary
	then using the binint dictionary make id binint tree"""
	#binint=read_binint(filename)
	count=0
	binint_id={}
	for i in binint.keys():
		for j in binint[i]['output'].keys():
			for k in range(len(binint[i]['output'][j]['ids'])):
				if binint[i]['output'][j]['merge']['merge'][k]==0:
					if int(binint[i]['output'][j]['ids'][k])==id:
						binint_id[count]={'merge': 0,
								'id': id,
								'mergeids': 0,
								'interaction': binint[i]
								}
						count+=1
				elif binint[i]['output'][j]['merge']['merge'][k]==1:
					for l in range(len(binint[i]['output'][j]['merge']['ids'][k])):
						if int(binint[i]['output'][j]['merge']['ids'][k][l])==id:
							binint_id[count]={'merge': 1,
									'id': id,
									'mergeids': binint[i]['output'][j]['merge']['ids'][k],
									'interaction': binint[i]
									}
							count+=1
				else:
					raise UserWarning('the following binint entry does not have a legal merge no '+
					     'input:\n'+binint[i])

	return binint_id	


	






def read_binintfile(id):
	f=open('binint.log','r')
	seek=find_positions('binint.log')
	interaction={}
	for i in range(len(seek)):
		f.seek(seek[i])
		value={}
		try:
			k=0
			while f.tell() < seek[i+1]:
				line=f.readline()
				a=line.split()
				value[k]=a
				k += 1
			int_id={}
			out_id={}
			out_obj_type=()
			#single_array = ()
			for k in value.keys():
				if value[k][0]=='output:':
					b1=value[k][1]
					b2=b1.replace("type=","")
					single_array=()
					if b2=='single':
						c1 = value[k][5].replace("id=","")
						if int(c1)==id:
							print 'success'
							int_id[id]={'typein': value[0][0].replace("type=",""),
									'time': value[0][1].replace("type",""),
									'objin1': [value[2][j].replace("type=","").replace("m=","").replace("R=","").replace("Eint=","").replace("id=","").replace("m0=","").replace("m1=","").replace("id0=","").replace("id1=","").replace("R1=","").replace("R2=","").replace("Eint1=","").replace("Eint2=","").replace("a=","").replace("e=","") for j in range(1,len(value[2]))],
									'objin2': [value[3][j].replace("type=","").replace("m=","").replace("R=","").replace("Eint=","").replace("id=","").replace("m0=","").replace("m1=","").replace("id0=","").replace("id1=","").replace("R1=","").replace("R2=","").replace("Eint1=","").replace("Eint2=","").replace("a=","").replace("e=","").replace("R0=","") for j in range(1,len(value[3]))],
									'objoutcome': [value[5][j].replace("type=","").replace("m=","").replace("R=","").replace("Eint=","").replace("id=","").replace("m0=","").replace("m1=","").replace(":","").replace("nstar=","").replace("nobj=","").replace("(=","").replace(")=","") for j in range(1,len(value[5]))]
								}

							#for l in range(int(int_id[id]['objoutcome'][1])):
							#	print l
							#	out_id[id'+'l]={l: 'dada'}
								#out_id[id]={'objout': [value[6+l][j].replace("type=","").replace("m=","").replace("R=","").replace("Eint=","").replace("id=","").replace("m0=","").replace("m1=","").replace("id0=","").replace("id1=","").replace("R1=","").replace("R2=","").replace("Eint1=","").replace("Eint2=","").replace("a=","").replace("e=","") for j in range(1,len(value[3]))]
										#}
							
							
							#for l in range(int(int_id[id]['objoutcome'][1])):
							#	int_id[id]['objout'][l]=[value[6+l][j].replace("type=","").replace("m=","").replace("R=","").replace("Eint=","").replace("id=","").replace("m0=","").replace("m1=","").replace("id0=","").replace("id1=","").replace("R1=","").replace("R2=","").replace("Eint1=","").replace("Eint2=","").replace("a=","").replace("e=","") for j in range(1,len(value[3]))]							
							print int_id
							print out_id
							
							



							#for j in range(2,len(value)):
							#print value[k][j]
							#c2=value[k][j].replace("m=","").replace("R=","").replace("Eint=","").replace("id=","")
							#single_array += (c2,)
						


						#print len(single_array)
						#if int(single_array[3])==14023:
							#print 'here'
			
					
						
				#elif value[k][0]=='outcome:':
				#	b1=value[k][1]
				#	b2=b1.replace("nstar=","")
				#	nstar=int(b2)
				#	b1=value[k][2]
				#	b2=b1.replace("nobj=","").replace(":","")
				#	nobj=int(b2)
				#	b1=value[k][3+nstar]
				#	b2=b1.replace("(","").replace(")","")
				#	outtype=b2

			#print (nstar,nobj,outtype,single_array)
		except IndexError:
			pass




