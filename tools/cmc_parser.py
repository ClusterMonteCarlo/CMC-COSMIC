import re

class collision:
	"""A class to store information from the collisions file
	
	the __init__ function takes the string of a single line of the collision
	file and converts it to useful numeric values
	
	Note that the instances can also be accessed as a dictionary by using
	collision_intance.__dict__
	
	"""
		
	def parse_one_coll(self,args):
		_,out_id,_,out_m,_ = re.split('=|\(|\)',args[2])
		_,id1,_,m1,_,_,id2,_,m2,_ = re.split('=|\(|\)|:',args[3])
		self.out_id = int(out_id)
		self.out_mass = float(out_m)
		self.in_ids = (int(id1),int(id2))
		self.in_masses =(float(m1),float(m2))
		km = args[5].split('=')[1]
		k1 = args[6].split('=')[1]
		k2 = args[7].split('=')[1].rstrip()
		self.out_type = int(km)
		self.in_types = (int(k1),int(k2))
		
	def parse_two_coll(self,args):
		_,out_id,_,out_m,_ = re.split('=|\(|\)',args[2])
		_,id1,_,m1,_,_,id2,_,m2,_,_,id3,_,m3,_ = re.split('=|\(|\)|:',args[3])
		self.out_id = int(out_id)
		self.out_mass = float(out_m)
		self.in_ids = (int(id1),int(id2),int(id3))
		self.in_masses = (float(m1),float(m2),float(m3))
		km = args[5].split('=')[1]
		k1 = args[6].split('=')[1]
		k2 = args[7].split('=')[1]
		k3 = args[8].split('=')[1].rstrip()
		self.out_type = int(km)
		self.in_types = (int(k1),int(k2),int(k3))

	def parse_three_coll(self,args):
		_,out_id,_,out_m,_ = re.split('=|\(|\)',args[2])
		_,id1,_,m1,_,_,id2,_,m2,_,_,id3,_,m3,_,_,id4,_,m4,_ = re.split('=|\(|\)|:',args[3])
		self.out_id = int(out_id)
		self.out_mass = float(out_m)
		self.in_ids = (int(id1),int(id2),int(id3),int(id4))
		self.in_masses = (float(m1),float(m2),float(m3),float(m4))
		km = args[5].split('=')[1]
		k1 = args[6].split('=')[1]
		k2 = args[7].split('=')[1]
		k3 = args[8].split('=')[1]
		k4 = args[9].split('=')[1].rstrip()
		self.out_type = int(km)
		self.in_types = (int(k1),int(k2),int(k3),int(k4))

	def __init__(self, string):
		args = string.split()
		
		self.time = float(args[0].split('=')[1])
		self.radius = float(args[4].split('=')[1].rstrip(')'))
		self.type = args[1]
		
		if self.type == 'single-single':
			self.parse_one_coll(args)
		elif self.type == 'binary-single':
			if len(args) == 8:
				self.parse_one_coll(args)
			else:
				self.parse_two_coll(args)
		elif self.type == 'binary-binary':
			if len(args) == 8:
				self.parse_one_coll(args)
			elif len(args) == 9:
				self.parse_two_coll(args)
			else:
				self.parse_three_coll(args)


class se_merger:
	"""
	A class to store information from the se_merge file
	
	the __init__ function takes the string of a single line of the merger
	file and converts it to useful numeric values
	
	Note that the instances can also be accessed as a dictionary by using
	collision_intance.__dict__
	
	"""
		
	def parse_disrupt_1(self,args): #same as the collision file
		_,out_id,_,out_m,_ = re.split('=|\(|\)',args[2])
		_,id1,_,m1,_,_,id2,_,m2,_ = re.split('=|\(|\)|:',args[3])
		self.out_id = int(out_id)
		self.out_mass = float(out_m)
		self.in_ids = (int(id1),int(id2))
		self.in_masses =(float(m1),float(m2))
		km = args[5].split('=')[1]
		k1 = args[6].split('=')[1]
		k2 = args[7].split('=')[1].rstrip()
		self.out_type = int(km)
		self.in_types = (int(k1),int(k2))
				
	def parse_disrupt_both(self,args):
		_,id1,_,m1,_ = re.split('=|\(|\)',args[2])
		_,id2,_,m2,_ = re.split('=|\(|\)',args[3])

		self.out_id = -1
		self.out_mass = -1
		self.in_ids = (int(id1),int(id2))
		self.in_masses = (float(m1),float(m2))
		k1 = args[5].split('=')[1]
		k2 = args[6].split('=')[1].rstrip()
		self.out_type = -1
		self.in_types = (int(k1),int(k2))

	def __init__(self, string):		
		args = string.split()
		
		self.time = float(args[0].split('=')[1])
		self.radius = float(args[4].split('=')[1].rstrip(')'))
		self.type = args[1]
		
		if self.type == 'disrupt1' or self.type == 'disrupt2':
			self.parse_disrupt_1(args)
		else:
			self.parse_disrupt_both(args)

class star:
	"""
	A container for keeping track of stars
	"""	
	def __init__(self,m,r,id,ktype=-1):
		self.m_MSUN = m
		self.r_RSUN = r
		self.id = id
		self.ktype = ktype
		
class binary:
	"""
	A container for keeping track of binaries
	"""	
	def __init__(self,m1,r1,id1,m2,r2,id2,a,e,k1=-1,k2=-1):
		self.star1 = star(m1,r1,id1,k1)
		self.star2 = star(m2,r2,id2,k2)
		self.a_AU = a
		self.e = e
		
	@property
	def m1_MSUN(self):
		return self.star1.m_MSUN
	@property
	def r1_RSUN(self):
		return self.star1.r_RSUN
	@property
	def id1(self):
		return self.star1.id
	@property
	def k1(self):
		return self.star1.ktype
	@property
	def m2_MSUN(self):
		return self.star2.m_MSUN
	@property
	def r2_RSUN(self):
		return self.star2.r_RSUN
	@property
	def id2(self):
		return self.star2.id
	@property
	def k2(self):
		return self.star2.ktype
	
	def same_binary(self,b1):
		if(b1.id1 == self.id1 and b1.id2 == self.id2 or
		   b1.id1 == self.id2 and b1.id2 == self.id1):
			return True
		else:
			return False
	
class triple:
	"""
	A container for keeping track of triples
	"""	
	def __init__(self,m1,r1,id1,m2,r2,id2,m3,r3,id3,a_in,e_in,a_out,e_out,k1=-1,k2=-1,k3=-1):
		self.star = star(m3,r3,id3,k3)
		self.binary = binary(m1,r1,id1,m2,r2,id2,a_in,e_in,k1,k2)
		self.a_out_AU = a_out
		self.e_out = e_out
		
	@property
	def m1_MSUN(self):
		return self.binary.m1_MSUN
	@property
	def r1_RSUN(self):
		return self.binary.r1_RSUN
	@property
	def id1(self):
		return self.binary.id1
	@property
	def k1(self):
		return self.binary.k1
	@property
	def m2_MSUN(self):
		return self.binary.m2_MSUN
	@property
	def r2_RSUN(self):
		return self.binary.r2_RSUN
	@property
	def id2(self):
		return self.binary.id2
	@property
	def k2(self):
		return self.binary.k2
	@property
	def a_in_AU(self):
		return self.binary.a_AU
	@property
	def e_in(self):
		return self.binary.e
	@property
	def m3_MSUN(self):
		return self.star.m_MSUN
	@property
	def r3(self):
		return self.star.r_RSUN
	@property
	def id3(self):
		return self.star.id
	@property
	def k3(self):
		return self.star.ktype

	
	
class binint:
	"""
	A class to store information from the binary interaction file
	
	the __init__ function takes the FILENAME, and pulls in lines till its a
	the next interaction, then converts it to useful numeric values
	
	Note that the instances can also be accessed as a dictionary by using
	collision_intance.__dict__
	
	"""
		
	def parse_single(self,args):
		m = float(args[2].split('=')[1])
		R = float(args[3].split('=')[1])
		id_s = args[5].split('=')[1].rstrip('\n')
		if ':' in id_s:
			self.collision = 0
		else:
			id_s = int(id_s)
		if len(args) > 6:
			ktype = int(args[6].split('=')[1].rstrip('\n'))
		else:
			ktype = -1
		return star(m,R,id_s,ktype)
	
	def parse_binary(self,args):
		m1 = float(args[2].split('=')[1])
		R1 = float(args[4].split('=')[1])
		id1 = args[8].split('=')[1]
		if ':' in id1:
			self.collision = 0
		else:
			id1 = int(id1)

		m2 = float(args[3].split('=')[1])
		R2 = float(args[5].split('=')[1])
		id2 = args[9].split('=')[1]
		if ':' in id2:
			self.collision = 0
		else:
			id2 = int(id2)

		if len(args) > 12:
			k1 = int(args[12].split('=')[1])
			k2 = int(args[13].split('=')[1].rstrip('\n'))
		else:
			k1 = -1
			k2 = -1
			
		a = float(args[10].split('=')[1])
		e = float(args[11].split('=')[1].rstrip('\n'))
		
		return binary(m1,R1,id1,m2,R2,id2,a,e,k1,k2)
	
	def parse_triple(self,args):
		m1 = float(args[2].split('=')[1])
		R1 = float(args[5].split('=')[1])
		id1 = args[11].split('=')[1]
		if ':' in id1:
			self.collision = 0
		else:
			id1 = int(id1)

		m2 = float(args[3].split('=')[1])
		R2 = float(args[6].split('=')[1])
		id2 = args[12].split('=')[1]
		if ':' in id2:
			self.collision = 0
		else:
			id2 = int(id2)

		m3 = float(args[4].split('=')[1])
		R3 = float(args[7].split('=')[1])
		id3 = args[13].split('=')[1]
		if ':' in id3:
			self.collision = 0
		else:
			id3 = int(id3)

		if len(args) > 18:
			k1 = int(args[18].split('=')[1])
			k2 = int(args[19].split('=')[1])
			k3 = int(args[20].split('=')[1].rstrip('\n'))
		else:
			k1 = -1
			k2 = -1
			k3 = -1

		a_in = float(args[14].split('=')[1])
		e_in = float(args[16].split('=')[1])

		a_out = float(args[15].split('=')[1])
		e_out = float(args[17].split('=')[1].rstrip('\n'))

		return triple(m1,R1,id1,m2,R2,id2,m3,R3,id3,a_in,e_in,a_out,e_out,k1,k2,k3)

	def parse_status(self,args):
		self.vesc_KMS = float(args[7].split('=')[1])
		self.de_gw = float(args[5].split('=')[1].rstrip())
		
	def parse_params(self,args):
		self.b = float(args[1].split('=')[1])
		self.v = float(args[2].split('=')[1].rstrip())
		
	def parse_input(self,args):
		if args[1] == 'type=single':
			self.in_singles.append(self.parse_single(args))
		else:
			self.in_binaries.append(self.parse_binary(args))
			
	def parse_output(self,args):
		if args[1] == 'type=single':
			self.out_singles.append(self.parse_single(args))
		elif args[1] == 'type=binary':
			self.out_binaries.append(self.parse_binary(args))
		else:
			self.out_triple.append(self.parse_triple(args))

			
	def __init__(self, binint_file): 
		self.in_singles  = []
		self.in_binaries = []

		self.out_singles  = []
		self.out_binaries = []
		
		self.out_triple = []
		
		self.collision = -1
		
		line = binint_file.readline()
		if '******' in line:
			line = binint_file.readline()


		args = line.split()
		
		if 'type=BS' == args[0]:
			self.type = 'binary-single'
		else:
			self.type = 'binary-binary'
			
		self.time = float(args[1].split('=')[1].rstrip())
		
		while '******' not in line and line != '':
			args = line.split()
			
			if args[0] == 'params:':
				self.parse_params(args)
			elif args[0] == 'input:':
				self.parse_input(args)
			elif args[0] == 'status:':
				self.parse_status(args)
			elif args[0] == 'output:':
				self.parse_output(args)		 
				
			line = binint_file.readline()
			
			
	def outputs(self):
		num_singles = len(self.out_singles)
		num_binaries = len(self.out_binaries)
		num_triples = len(self.out_triple)
		
		return [num_singles,num_binaries,num_triples]
	
	
	def number_exchanges(self):
		num_exchange = 0
		for bo in self.out_binaries:
			same = False
			for bi in self.in_binaries:
				if bo.same_binary(bi):
					same = True
			if not same:
				num_exchange += 1
		return num_exchange
	
	def in_ids(self):
		ids = []
		for star in self.in_singles:
			ids.append(star.id)
		for binary in self.in_binaries:
			ids.append(binary.id1)
			ids.append(binary.id2)
		return ids
	
	def out_ids(self):
		ids = []
		for star in self.out_singles:
			ids.append(star.id)
		for binary in self.out_binaries:
			ids.append(binary.id1)
			ids.append(binary.id2)
		for triple in self.out_triple:
			ids.append(triple.id1)
			ids.append(triple.id2)
			ids.append(triple.id3)
		return ids

class conversion_file:
	"""
	A class to store information from the unit conversion script 

	the __init__ function takes the FILENAME, and saves all the units as class objects
	"""
	
	def parse_conv_file(self,filename):
		with open(filename) as f:
			lines = f.readlines()
		for line in lines:
			line_array = line.split('=')
			if len(line_array) == 2:
				if line_array[0] == 'massunitcgs':
					self.mass_cgs = float(line_array[1])
				if line_array[0] == 'massunitmsun':
					self.mass_msun = float(line_array[1])
				if line_array[0] == 'mstarunitcgs':
					self.mass_cgs_mstar = float(line_array[1])
				if line_array[0] == 'mstarunitmsun':
					self.mass_msun_mstar = float(line_array[1])
				if line_array[0] == 'lengthunitcgs':
					self.length_cgs = float(line_array[1])
				if line_array[0] == 'lengthunitparsec':
					self.length_parsec = float(line_array[1])
				if line_array[0] == 'timeunitcgs':
					self.time_cgs = float(line_array[1])
				if line_array[0] == 'timeunitsmyr':
					self.time_myr = float(line_array[1])
				if line_array[0] == 'nbtimeunitcgs':
					self.time_cgs_nb = float(line_array[1])
				if line_array[0] == 'nbtimeunitsmyr':
					self.time_myr_nb = float(line_array[1])
					break


	def __init__(self,filename):
	   self.parse_conv_file(filename)

	@property
	def time_nb_relax(self):
		return self.time_myr / self.time_myr_nb

def load_ineraction_files(file_prefix):
	"""
	Actually does the work of loading the binint, collision, and se_merger files into
	the pythonic classes we've developed. Takes the file prefix as an input and returns

	(units,[list of binints],[list of collisions],[list of mergers]) 
    """
	binintfile=file_prefix+'.binint.log'
	collfile=file_prefix+'.collision.log'
	mergefile=file_prefix+'.semergedisrupt.log'
	unitsfile=file_prefix+'.conv.sh'

	with open(binintfile, 'r') as f:
		num_int = sum([1 for line in f if '****' in line])
		f.seek(0) #need to count the encounters then reset the file handle
		binints = [binint(f) for __ in range(num_int)]

	with open(mergefile,'r') as f:
		next(f) #skip the header
		mergers = [se_merger(line) for line in f]
		
	with open(collfile,'r') as f:
		next(f) #skip the header
		collisions = [collision(line) for line in f]

	units = conversion_file(unitsfile) 

	return (units,binints,mergers,collisions)
