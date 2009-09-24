#!/opt/local/bin//python

from numpy import *
import gzip
import math
import constants
import scripts1
import scripts2
import scripts3
import history_cmc
import subprocess

def read_units(filestring):
	"""reads from the supplied conv file and stores the physical units"""
	s = filestring+'.conv.sh'
	f=open(s,'r')
	count=0
	while count<10:
		a=f.readline()
		b=a.split('=')
		if b[0]=='massunitcgs':
			massunitcgs = float(b[1])
			count+=1
		elif b[0]=='massunitmsun':
			massunitmsun = float(b[1])
			count+=1
		elif b[0]=='mstarunitcgs':
			mstarunitcgs = float(b[1])
			count+=1
		elif b[0]=='mstarunitmsun':
			mstarunitmsun = float(b[1])
			count+=1
		elif b[0]=='lengthunitcgs':
			lengthunitcgs = float(b[1])
			count+=1
		elif b[0]=='lengthunitparsec':
			lengthunitparsec = float(b[1])
			count+=1
		elif b[0]=='timeunitcgs':
			timeunitcgs = float(b[1])
			count+=1
		elif b[0]=='timeunitsmyr':
			timeunitsmyr = float(b[1])
			count+=1
		elif b[0]=='nbtimeunitcgs':
			nbtimeunitcgs = float(b[1])
			count+=1
		elif b[0]=='nbtimeunitsmyr':
			nbtimeunitsmyr = float(b[1])
			count+=1
	f.close()
	units = []
	unittype = [('m_cgs', float), ('m_msun', float), ('mstar_cgs', float), ('mstar_msun', float), ('l_cgs', float), ('l_pc', float), ('t_cgs', float),('t_myr', float), ('nbt_cgs', float), ('nbt_myr', float)]
	units.append((massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr))
	units = array(units, dtype=unittype)
	#return (massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr)
	return units

def find_correlations(string,snapno, binary):
	"""finds anything that may be interesting to see correlations:
		string: file string
		snapno: snapno as a string
		binary: whether there were primordial binaries or not.  0: not, 1: yes"""

	snap = string+'.snap'+snapno+'.dat.gz'
	dyn = string+'.dyn.dat'
	from glob import glob
	out = string+'.out.dat'
	bin = string+'.bin.dat'

	units = read_units(string)
	convert_v = units[0]['l_cgs']/units[0]['nbt_cgs']/1e5
	convert_rho = units[0]['m_msun']/(units[0]['l_pc'])**3.
	t_myr = find_t_myr(string, snapno)
	m_to = find_MS_turnoff(t_myr*1000000.)

	f_snap = gzip.open(snap,'rb')
	line_snap = f_snap.readline()
	a_snap = line_snap.split()
	b_snap = a_snap[1]
	c_snap = b_snap.split('=')
	t_snap = float(c_snap[1])
	print t_snap
	
	dyndata = loadtxt(dyn)
	for i in range(len(dyndata)):
		if (dyndata[i,0])==t_snap:
			m = dyndata[i,4]
			rc = dyndata[i,7]
			rh = dyndata[i,20]
			rho0 = dyndata[i,21]
			#lo behold: time unit in central v calculation is in nbody units
			v0_rms = dyndata[i,23]
			print 'v=',v0_rms,convert_v	
			m_msun,rc_pc,rh_pc,rho0_msun_pc3,v0_rms_km_s = m*units[0]['m_msun'],rc*units[0]['l_pc'],rh*units[0]['l_pc'],rho0*convert_rho,v0_rms*convert_v
	
	try:
		bindata = loadtxt(bin)
		f_bc_array=[]
		f_b_array=[]
		rhb_pc_array=[]
		rhs_pc_array=[]
		print 'opened bindata'
		for i in range(len(bindata)):
			if (bindata[i,0])-t_snap < 0.00001 and (bindata[i,0])-t_snap > -0.00001:
				print 'found time'
				f_bc_array.append(bindata[i,10])
				f_b_array.append(bindata[i,11])
				rhb_pc_array.append(bindata[i,5]*units[0]['l_pc'])
				rhs_pc_array.append(bindata[i,4]*units[0]['l_pc'])
			f_bc=mean(f_bc_array)
			f_b=mean(f_b_array)
			rhb_pc=mean(rhb_pc_array)
			rhs_pc=mean(rhs_pc_array)
	except IOError:
		f_bc = 0.0
		f_b = 0.0
		rhb_pc = 0.0
		rhs_pc = rh_pc	
	
	BSS_count_core = find_BSS_core(string, snapno, rc)
	BSS_count = find_BSS(string,snapno)
	total_BSS = BSS_count[0]+BSS_count[1]
	print 'look here', rc, rc_pc
	
	total_BSS_core = BSS_count_core[0]+BSS_count_core[1]
	history = history_cmc.history_maker(BSS_count[3], BSS_count[4], string, binary)
	branching_ratio = history_cmc.branching_ratio(history, binary)
	#print branching_ratio
	BSS_binint = total_BSS*branching_ratio['binint']
	BSS_pure_binint = total_BSS*branching_ratio['pure_binint']
	BSS_coll = total_BSS*branching_ratio['coll']
	BSS_pure_coll = total_BSS*branching_ratio['pure_coll']
	BSS_merger = total_BSS*branching_ratio['merger']
	BSS_pure_merger = total_BSS*branching_ratio['pure_merger']
	

	print "got BSS"
	
	giants_count = find_giants(string,snapno)
	giants_core_count = find_giants_core(string, snapno, rc)
	print "got giants"
	total_giants = (giants_count[0]+giants_count[1])
	total_giants_core = (giants_core_count[0] + giants_core_count[1])

	try:
		outfile=open(out,'r')
	except IOError:
		print 'could not find out.dat'
		out=glob('*.o*')[0]
		outfile = open(out,'r')
	line='initial'
	tcoll_array=[]
	while len(line)>0:
		line=outfile.readline()
		if line.rfind('Tcoll')>-1:
			tcoll=line.split()
			tcoll
			if float(tcoll[2])>(0.001*t_myr - 1.) and float(tcoll[2])<=(0.001*t_myr) and tcoll[6]!='nan':
				tcoll_array.append(float(tcoll[6]))
	mean_gamma=(mean(tcoll_array))**(-1.)

	dtype = [('m',float), ('rc', float), ('rh',float), ('rho',float), ('gamma',float), ('v0',float), ('t',float), ('mto', float), 
			('fbc',float), ('fb',float), ('rhb',float), ('rhs',float), ('bss_sing',int), ('bss_bin',int), ('bssc_sing',int),
			('bssc_bin',int), ('bss_binint',float), ('bss_pure_binint',float), ('bss_coll',float), ('bss_pure_coll',float), 
			('bss_merger',float), ('bss_pure_merger',float), ('rgb_sing',int), ('rgb_bin',int), ('rgbc_sing',int), 
			('rgbc_bin',int)]
	dummy=[]
	dummy.append( (m_msun,rc_pc,rh_pc,rho0_msun_pc3,mean_gamma,v0_rms_km_s,t_myr,m_to,f_bc,f_b,rhb_pc,rhs_pc,BSS_count[0],BSS_count[1],BSS_count_core[0], BSS_count_core[1], BSS_binint, BSS_pure_binint, BSS_coll, BSS_pure_coll, BSS_merger, BSS_pure_merger, giants_count[0], giants_count[1], giants_core_count[0], giants_core_count[1]) )
	dummy=array(dummy,dtype=dtype)
	

	#for i in range(len(dummy)):
	#	print "%.3f " %(dummy[i],)
#
#	print "\n"

	return dummy
	
	

def find_t_myr(filestring, snapno):
	"""goes in the given snapshot and finds the physical time corresponding to that snap"""
	snapfile = filestring+'.snap'+snapno+'.dat.gz'

	f=gzip.open(snapfile,'rb')
	line=f.readline()
	a=line.split()
	b=a[1]
	c=b.split('=')
	t=float(c[1])
	d=read_units(filestring)
	t_myr=t*d['t_myr'][0]
	f.close()
	return (t_myr)

def find_MS_turnoff(t):
	"""given the time in year it finds the MS turn-off mass in Solar masses.  Very simple now.  Need to make the MS lifetime formula better. """
	lm = (9.921 - log10(t))/3.6648
	m = 10**lm
	return(m)	

def find_MS_TO(z):
	eta = log10(z/0.02)
	m_hook = 1.0185 + 0.16015*eta + 0.0892*eta**2.
	m_HeF = 1.995 + 0.25*eta + 0.087*eta**2.
	m_FGB = 13.048*(z/0.02)**0.06/(1+0.0012*(0.02/z)**1.27)

	x = max([0.95,min([0.95-0.03*(eta+0.30103)]),0.99])

	return (eta,m_hook,m_HeF,m_FGB,x)
	


def find_giants(file_string,snapno):
	"""takes the file string for the runs.  Returns the number of giants in single and the 
	number of giants in binaries."""
	snapfile=file_string+'.snap'+snapno+'.dat.gz'
	convfile=file_string+'.conv.sh'
	wfile1=file_string+'.snap'+snapno+'.giants.dat'

	t_yr = (find_t_myr(file_string, snapno))*10**6

	f=gzip.open(snapfile,'rb')
	f1 = open(wfile1, 'w')
	giants_sing_count = 0
	giants_bin_count = 0

	#print header
	f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	#f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	for line in f:	
		a=line.split()
		if a[0]!='#' and a[0]!='#1:id':
			if a[7]=='0' and float(a[14])==3.:
				for i in range(len(a)):
					f1.write("%s " % (a[i]))
				giants_sing_count += 1
				f1.write("\n")
				
			
			if a[7]!='0' and a[8]!='0' and a[9]!='0': 
				if float(a[17])==3.:
					for i in range(len(a)):
						f1.write("%s " % (a[i]))
					giants_bin_count += 1
					f1.write("\n")
				elif float(a[18])==3.:
					for i in range(len(a)):
						f1.write("%s " % (a[i]))
					giants_bin_count += 1
					f1.write("\n")
	f1.close()
	return (giants_sing_count, giants_bin_count,t_yr/10**6)


def find_giants_core(file_string,snapno, rc):
	"""takes the snapshot file, the conv.sh file, the file name to write the BSS in singles and the file name to write the BSS in binaries at that snapshot.  Returns the number of BSS in single within the core and the number of BSS in binaries within the core."""
	snapfile=file_string+'.snap'+snapno+'.dat.gz'
	convfile=file_string+'.conv.sh'
	
	t_yr = (find_t_myr(file_string, snapno))*10**6
	t_gyr = '%.1f' %(t_yr/10.**9.)

	#wfile1=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.single_BSS.dat'
	#wfile2=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.binary_BSS.dat'

	m = find_MS_turnoff(float(t_yr))
	m_cut = (1.+0.1)*m
	print "looking at a snap at t = %f Gyr, MS turnoff is %f MSun" % (t_yr/10**9, m)

	f=gzip.open(snapfile,'rb')
	f.seek(0)
	#f1 = open(wfile1,'w')
	#f2 = open(wfile2,'w')
	giants_sing_count = 0
	giants_bin_count = 0

	#print header
	#f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	#f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	giants_ids = []
	giants_positions = []
	for line in f:	
		a=line.split()
		if a[0]!='#' and a[0]!='#1:id':
			if float(a[2]) <= float(rc):
				if a[7]=='0' and float(a[14])==3.:
					#for i in range(len(a)):
						#f1.write("%s " % (a[i]))
					giants_sing_count += 1
					#f1.write("\n")
					
				
				if a[7]!='0' and a[8]!='0' and a[9]!='0': 
					if float(a[17])==3.:
						#for i in range(len(a)):
							#f2.write("%s " % (a[i]))
						giants_bin_count += 1
						#f2.write("\n")
					if float(a[18])==3.:
						#for i in range(len(a)):
							#f2.write("%s " % (a[i]))
						giants_bin_count += 1
						#f2.write("\n")
	f.close()
	#f1.close()
	#f2.close()
	return (giants_sing_count, giants_bin_count,t_yr/10**6)


	

def find_BSS(file_string,snapno):
	"""takes the snapshot file, the conv.sh file, the file name to write the BSS in singles and the file name to write the BSS in binaries at that snapshot.  Returns the number of BSS in single and the number of BSS in binaries."""
	from glob import glob
	binint=glob('*.binint.log')
	binary=1
	if len(binint)==0:
		binary=0
	print 'binary=',binary
	snapfile=file_string+'.snap'+snapno+'.dat.gz'
	units = read_units(file_string)
	
	t_yr = (find_t_myr(file_string, snapno))*10**6
	t_gyr = '%.1f' %(t_yr/10.**9.)

	wfile1=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.single_BSS.dat'
	wfile2=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.binary_BSS.dat'

	m = find_MS_turnoff(float(t_yr))
	m_cut = (1.+0.05)*m
	print "looking at a snap at t = %f Gyr, MS turnoff is %f MSun" % (t_yr/10**9, m)

	f=gzip.open(snapfile,'rb')
	f.seek(0)
	f1 = open(wfile1,'w')
	f2 = open(wfile2,'w')
	bss_sing_count = 0
	bss_bin_count = 0

	print 'header'
	f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN] 24.binint? 25.pure_binint? 26.coll? 27.pure_coll? 28.se? 29.pure_se? 30.Leff(LSUN) 31.Teff(K)\n")
	f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN] 24.binint? 25.pure_binint? 26.coll? 27.pure_coll? 28.se? 29.pure_se? 30.Leff(LSUN) 31.Teff(K)\n")
	BSS_ids = []
	BSS_positions = []
	for line in f:
		binint, pure_binint, coll, pure_coll, se, pure_se = 0, 0, 0, 0, 0, 0
		coll_time = []
		a=line.split()
		if a[0]!='#' and a[0]!='#1:id':
			if float(a[1])>m_cut and a[7]=='0' and float(a[14])<2.:
				for i in range(len(a)):
					f1.write("%s " % (a[i]))
				bss_sing_count += 1
				dummy_id = int(a[0])
				BSS_ids.append(dummy_id)
				BSS_positions.append(float(a[2]))
				#f1.write("\n")
				#now get info about the BSS
				hist=history_cmc.history_maker([dummy_id], [float(a[2])], file_string, binary)
				if binary==1:
					if len(hist[dummy_id]['binint']['binint'].keys())>0:
						binint = 1
						if len(hist[dummy_id]['coll'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
							pure_binint = 1
				if len(hist[dummy_id]['coll'].keys())>0:
					coll = 1
					for i in hist[dummy_id]['coll'].keys():
						print 'key', i
						coll_time.append(hist[dummy_id]['coll'][i]['coll_params']['time'])
						coll_time.sort()
						last_coll=coll_time[-1]

					if len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
						pure_coll = 1
				if len(hist[dummy_id]['se'].keys())>0:
					se = 1
					if len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['coll'].keys())==0:
						pure_se = 1

				f1.write("%d %d %d %d %d %d " %(binint, pure_binint, coll, pure_coll, se, pure_se))
				Teff = find_T(float(a[15]), float(a[16]))
				Leff = float(a[15])
				f1.write("%f %f\n" %(Leff, Teff))


				
			
			if a[7]!='0' and a[8]!='0' and a[9]!='0': 
				if (float(a[8]))>m_cut and float(a[17])<2.:
					for i in range(len(a)):
						f2.write("%s " % (a[i]))
					bss_bin_count += 1
					dummy_id=int(a[10])
					BSS_ids.append(dummy_id)
					BSS_positions.append(float(a[2]))

					#now get info about the BSS
					hist=history_cmc.history_maker([dummy_id], [float(a[2])], file_string, binary)
					if binary==1 and len(hist[dummy_id]['binint']['binint'].keys())>0:
						binint = 1
						if len(hist[dummy_id]['coll'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
							pure_binint = 1
					if len(hist[dummy_id]['coll'].keys())>0:
						coll = 1
						for i in hist[dummy_id]['coll'].keys():
							coll_time.append(hist[dummy_id]['coll'][i]['coll_params']['time'])
							coll_time.sort()
							last_coll=coll_time[-1]
						if len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
							pure_coll = 1
					if len(hist[dummy_id]['se'].keys())>0:
						se = 1
						if len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['coll'].keys())==0:
							pure_se = 1

					f2.write("%d %d %d %d %d %d " %(binint, pure_binint, coll, pure_coll, se, pure_se))
					L0,R0,L1,R1=float(a[19]), float(a[21]), float(a[20]), float(a[22])
					T0=find_T(L0, R0)
					T1=find_T(L1, R1)
					Leff=L0+L1
					Teff=(L0*T0+L1*T1)/Leff
					f2.write("%f %f\n" %(Leff, Teff))


					#f2.write("\n")
				if (float(a[9]))>m_cut and float(a[18])<2.:
					for i in range(len(a)):
						f2.write("%s " % (a[i]))
					bss_bin_count += 1
					dummy_id = int(a[11])
					BSS_ids.append(dummy_id)
					BSS_positions.append(float(a[2]))
					
					#now get info about the BSS
					hist=history_cmc.history_maker([dummy_id], [float(a[2])], file_string, binary)
					if binary==1 and len(hist[dummy_id]['binint']['binint'].keys())>0:
						binint = 1
						if len(hist[dummy_id]['coll'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
							pure_binint = 1
					if len(hist[dummy_id]['coll'].keys())>0:
						coll = 1
						for i in hist[dummy_id]['coll'].keys():
							coll_time.append(hist[dummy_id]['coll'][i]['coll_params']['time'])
							coll_time.sort()
							last_coll=coll_time[-1]
						if len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['se'].keys())==0:
							pure_coll = 1
					if len(hist[dummy_id]['se'].keys())>0:
						se = 1
						if len(hist[dummy_id]['binint'].keys())==0 and len(hist[dummy_id]['coll'].keys())==0:
							pure_se = 1

					f2.write("%d %d %d %d %d %d " %(binint, pure_binint, coll, pure_coll, se, pure_se))

					L0,R0,L1,R1=float(a[19]), float(a[21]), float(a[20]), float(a[22])
					T0=find_T(L0, R0)
					T1=find_T(L1, R1)
					Leff=L0+L1
					Teff=(L0*T0+L1*T1)/Leff
					f2.write("%f %f\n" %(Leff, Teff))

					

					#f2.write("\n")
	f.close()
	f1.close()
	f2.close()
	return (bss_sing_count, bss_bin_count,t_yr/10**6, BSS_ids, BSS_positions)



def find_BSS_core(file_string,snapno, rc):
	"""takes the snapshot file, the conv.sh file, the file name to write the BSS in singles and the file name to write the BSS in binaries at that snapshot.  Returns the number of BSS in single within the core and the number of BSS in binaries within the core."""
	snapfile=file_string+'.snap'+snapno+'.dat.gz'
	convfile=file_string+'.conv.sh'
	
	t_yr = (find_t_myr(file_string, snapno))*10**6
	t_gyr = '%.1f' %(t_yr/10.**9.)

	#wfile1=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.single_BSS.dat'
	#wfile2=file_string+'.snap'+snapno+'.'+str(t_gyr)+'.core.binary_BSS.dat'

	m = find_MS_turnoff(float(t_yr))
	m_cut = (1.+0.05)*m
	print "looking at a snap at t = %f Gyr, MS turnoff is %f MSun" % (t_yr/10**9, m)

	f=gzip.open(snapfile,'rb')
	f.seek(0)
	#f1 = open(wfile1,'w')
	#f2 = open(wfile2,'w')
	bss_sing_count = 0
	bss_bin_count = 0

	#print header
	#f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	#f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	BSS_ids = []
	BSS_positions = []
	for line in f:	
		a=line.split()
		if a[0]!='#' and a[0]!='#1:id':
			if float(a[2]) <= float(rc):
				if float(a[1])>m_cut and a[7]=='0' and float(a[14])<2.:
					#for i in range(len(a)):
						#f1.write("%s " % (a[i]))
					bss_sing_count += 1
					BSS_ids.append(int(a[0]))
					BSS_positions.append(float(a[2]))
					#f1.write("\n")
					
				
				if a[7]!='0' and a[8]!='0' and a[9]!='0': 
					if (float(a[8]))>m_cut and float(a[17])<2.:
						#for i in range(len(a)):
							#f2.write("%s " % (a[i]))
						bss_bin_count += 1
						BSS_ids.append(int(a[10]))
						BSS_positions.append(float(a[2]))
						#f2.write("\n")
					if (float(a[9]))>m_cut and float(a[18])<2.:
						#for i in range(len(a)):
							#f2.write("%s " % (a[i]))
						bss_bin_count += 1
						BSS_ids.append(int(a[11]))
						BSS_positions.append(float(a[2]))
						#f2.write("\n")
	f.close()
	#f1.close()
	#f2.close()
	return (bss_sing_count, bss_bin_count,t_yr/10**6, BSS_ids, BSS_positions)



def BS_r_cum(singlefile, binaryfile, giantfile, rc, writefile, nbin):
	sing_data=loadtxt(singlefile)
	bin_data=loadtxt(binaryfile)
	giant_data=loadtxt(giantfile)
	f=open(writefile,'w')
	
	#initialize
	sing_roverrc=zeros(len(sing_data))
	bin_roverrc=zeros(len(bin_data))
	both_roverrc=zeros(len(bin_data)+len(sing_data))
	g_roverrc=zeros(len(giant_data))
	
	#getting the scaled radial distribution length scaled to the core radius 
	for i in range(len(sing_data)):
		sing_roverrc[i]=sing_data[i,2]/rc
	for i in range(len(bin_data)):
		bin_roverrc[i]=bin_data[i,2]/rc
	for i in range(len(giant_data)):
		g_roverrc[i]=giant_data[i,2]/rc
	for i in range(len(sing_data)):
		both_roverrc[i]=sing_roverrc[i]
	for i in range(len(bin_data)):
		both_roverrc[i+len(sing_data)]=bin_roverrc[i]


	nbinmem = len(g_roverrc)/nbin
	sing_cum, bin_cum, all_cum = 0., 0., 0.
	for i in range(1,nbin):
		temp_sing, temp_bin, temp_all = 0., 0., 0.
		rmin, rmax = g_roverrc[(i-1)*nbinmem], g_roverrc[i*nbinmem]
		print rmin, rmax, (rmin+rmax)/2.
		for j in range(len(sing_roverrc)):
			if sing_roverrc[j] < rmax and sing_roverrc[j]>=rmin:
				temp_sing += 1
		for j in range(len(bin_roverrc)):
			if bin_roverrc[j] < rmax and bin_roverrc[j]>=rmin:
				temp_bin += 1
		for j in range(len(both_roverrc)):
			if both_roverrc[j] < rmax and both_roverrc[j]>=rmin:
				temp_all += 1
		r = (rmin+rmax)/2.
		temp_sing = temp_sing/float(nbinmem)
		temp_bin = temp_bin/float(nbinmem)
		temp_all = temp_all/float(nbinmem)
		sing_cum += temp_sing
		bin_cum += temp_bin
		all_cum += temp_all

		f.write("%f %f %f %f %f %f %f\n" %(r, temp_sing, temp_bin, temp_all, sing_cum, bin_cum, all_cum) )
		#print r, temp_sing, temp_bin, temp_all, sing_cum, bin_cum, all_cum
	f.close()

	
				
def BS_r(singlefile, binaryfile, giantfile, rc, writefile):
	"""gives the normalized radial distribution of the BSS population
	to get this there is some preprocessing needs to be done
	1. run hrdiag on your favorite snapshot
	2. from hrdiag file filter out the singles with the BSS cut as you like in one file
	3. from hrdiag file filter out the binaries with the BSS cut as you like in another file
	4. from hrdiag file filter out the Giants 
	5. feed these 3 files and the value of rc at that time from *.dyn.dat"""
	sing_data=loadtxt(singlefile)
	bin_data=loadtxt(binaryfile)
	giant_data=loadtxt(giantfile)
	f=open(writefile,'w')
	
	#initialize
	sing_roverrc=zeros(len(sing_data))
	bin_roverrc=zeros(len(bin_data))
	both_roverrc=zeros(len(bin_data)+len(sing_data))
	g_roverrc=zeros(len(giant_data))
	
	#getting the scaled radial distribution length scaled to the core radius 
	for i in range(len(sing_data)):
		sing_roverrc[i]=sing_data[i,8]/rc
	for i in range(len(bin_data)):
		bin_roverrc[i]=bin_data[i,8]/rc
	for i in range(len(giant_data)):
		g_roverrc[i]=giant_data[i,8]/rc
	for i in range(len(sing_data)):
		both_roverrc[i]=sing_roverrc[i]
	for i in range(len(bin_data)):
		both_roverrc[i+len(sing_data)]=bin_roverrc[i]

	#print sing_roverrc, bin_roverrc, both_roverrc
	#print len(sing_roverrc), len(bin_roverrc), len(both_roverrc)

	
	
	#get histogram for the scaled radii
	amin, amax = min(g_roverrc[:]), max(g_roverrc[:])
	hd_both=histogram(both_roverrc,bins=50, range=(amin, amax))
	hd_sing=histogram(sing_roverrc,bins=50, range=(amin, amax))
	hd_bin=histogram(bin_roverrc,bins=50, range=(amin, amax))
	hd_g=histogram(g_roverrc,bins=50, range=(amin, amax))

	print hd_both, hd_sing, hd_bin
	print len(hd_both[0]), len(hd_sing[0]), len(hd_bin[0])

	#normalize the numbers with the giant population
	normalized_both=zeros((len(hd_both[0]),3))
	normalized_sing=zeros((len(hd_sing[0]),3))
	normalized_bin=zeros((len(hd_sing[0]),3))
	#print normalized_both, normalized_sing, normalized_bin
	#print len(normalized_both), len(normalized_sing), len(normalized_bin)
	for i in range(1,len(hd_both[1])):
		count=0
		for j in range(len(g_roverrc)):
			if g_roverrc[j]<hd_both[1][i] and g_roverrc[j]>hd_both[1][i-1]:
				count += 1
		print count
		normalized_both[i-1]=[ float(hd_both[0][i-1])/count , (hd_both[0][i-1]**(-0.5) + count**-0.5)*float(hd_both[0][i-1])/count , hd_both[1][i-1] ]
		f.write ("%f %f %f\n" % (normalized_both[i-1,2],normalized_both[i-1,0],normalized_both[i-1,1]))
	
	for i in range(1,len(hd_sing[1])):
		count=0
		for j in range(len(g_roverrc)):
			if g_roverrc[j]<hd_sing[1][i] and g_roverrc[j]>hd_sing[1][i-1]:
				count += 1
		print count
		normalized_sing[i-1]=[ float(hd_sing[0][i-1])/count , (hd_sing[0][i-1]**(-0.5) + count**-0.5)*float(hd_sing[0][i-1])/count , hd_sing[1][i-1] ]
	
	for i in range(1,len(hd_bin[1])):
		count=0
		for j in range(len(g_roverrc)):
			if g_roverrc[j]<hd_bin[1][i] and g_roverrc[j]>hd_bin[1][i-1]:
				count += 1
		print count
		normalized_bin[i-1]=[ float(hd_bin[0][i-1])/count , (hd_bin[0][i-1]**(-0.5) + count**-0.5)*float(hd_bin[0][i-1])/count, hd_bin[1][i-1] ]
	
	#print normalized_both, normalized_sing, normalized_bin
	#print len(normalized_both), len(normalized_sing), len(normalized_bin)
	

	#now plot
	import gracePlot
	gpl=gracePlot.gracePlot()
	gpl.hold()
	gpl.plot(normalized_both[:,2],normalized_both[:,0],normalized_both[:,1])

	f.close()

	#gpl.plot(normalized_sing[:,2],normalized_sing[:,0],normalized_sing[:,1])
	#gpl.plot(normalized_bin[:,2],normalized_bin[:,0],normalized_bin[:,1])	


#
#	return (normalized_both, normalized_sing, normalized_bin)
		

def find_problematic(id,filename):
	"""takes the id and then returns if that id was somewhere in the given file"""	
	f=open(filename,'r')
	f.seek(0)
	try:
		while f:
			line=f.readline()
			if line.rfind(id)>-1:
				#print "%s\n" % (line)
				problem=line
			elif line=='':
				raise StopIteration()
	except StopIteration:
		f.seek(0)
	f.close()
	return (problem)

def c_w0(c):
	"""taking the concentration parameter gives the w0 value within range c (0.5,3.5)"""
	w0 = -2.47882 + 9.53235*c - 0.843304*c**2 - 2.34742*c**3 + 1.22119*c**4 - 0.167361*c**5
	return w0 

def w0_c(w0):
	"""taking the w0 gives the concentration parameter within range w0 (2,15)"""
	c=0.247266 + 0.2079*w0 - 0.057241*w0**2 + 0.0145234*w0**3 - 0.00120862*w0**4 -0.0000330051*w0**5
	return c

def metallicity(m,s):
	"""give the metalicity and tell how to convert
	options are:
		ztofe/h
		fe/htoz"""

	if s=='ztofe/h':
		m1 = log10(m/0.02)
	elif s=='fe/htoz':
		m1 = 0.02*10**m
	return m1

def find_rc(file_string,snapno):
	"""finds r_c at this snapfile"""
	snapfile=file_string+'.snap'+snapno+'.dat.gz'
	dynfile=file_string+'.dyn.dat'

	snap=gzip.open(snapfile,'rb')
	line=snap.readline()
	data = line.split()
	data1=data[1].split('=')
	time=data1[1]
	
	dyndata=loadtxt(dynfile)
	for i in range(len(dyndata)):
		if dyndata[i,0]==float(time):
			rc = dyndata[i,7]

	return rc
	


def filter_single(filter,T):
	data=loadtxt(filter)
	sum=0.0
	normalize=0.0
	const=2*constants.h*constants.c**2
	bin_min = data[0,0] - (data[1,0] - data[0,0])
	for i in range (len(data)):
		sum+= const / (data[i,0]*constants.Angstrom)**5 / (exp(constants.c * constants.h/(data[i,0]*constants.Angstrom)/constants.k / T) -1) * (data[i,0]-bin_min)*constants.Angstrom * data[i,1] 
		normalize+= (data[i,0]-bin_min)*constants.Angstrom * data[i,1]
		bin_min = data[i,0]
	normalized = (sum/normalize)/4/pi/85000./constants.PC
	mag = 1.2 -2.5*log10(normalized)
	return (const,normalized,normalize,sum,mag)

def find_T(L,R):
	"""L in LSUN and R in RSUN: returns T in K"""
	T= (L*constants.Lsun / (4*pi*(R*constants.Rsun)**2) / constants.sigma)**0.25
	return (T)

def find_g(M,R):
	"""M in MSun and R in RSun: returns g in CGS"""
	g = constants.G*M*constants.Msun/(R*constants.Rsun)**2.
	return(g)


def filter_singles(filterdata,L,R):
	#filterdata = loadtxt(filter)
	sum = 0.0
	sun_sum = 0.0 
	normalize = 0.0 
	const = 2*constants.h*constants.c**2
	T = (L*constants.Lsun / (4*pi*(R*constants.Rsun)**2) / constants.sigma)**0.25
	
	bin_min = filterdata[0,0] - (filterdata[1,0] - filterdata[0,0])
	
	for i in range(len(filterdata)):
		a1 = const/(filterdata[i,0]*constants.Angstrom)**5
		a2 = ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/T) -1 )**-1.
		sun_a2 = ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/constants.Tsun) -1 )**-1.
		a3 = (filterdata[i,0] - bin_min)*constants.Angstrom

		sum += a1*a2*a3*filterdata[i,1]
		sun_sum += a1*sun_a2*a3*filterdata[i,1]
		normalize += a3*filterdata[i,1]	
		
		#sum += const / (filterdata[i,0]*constants.Angstrom)**5. / (exp(constants.c * constants.h / constants.k / filterdata[i,0]*constants.Angstrom / T) - 1) * (filterdata[i,0] - bin_min)*constants.Angstrom * filterdata[i,1]
	
		#normalize += (filterdata[i,0] - bin_min)*constants.Angstrom * filterdata[i,1]		
		
		bin_min = filterdata[i,0]
		#print (sum, normalize, bin_min)
	#print (sum,sun_sum, T)
	#normalized = (sum/normalize) * (4*pi*R**2) * pi 
	#sun_normalized = (sun_sum/normalize) * (4*pi*R**2) * pi
	#fraction = normalized/L/constants.Lsun
	#mag = -2.5*log10(normalized/sun_normalized)
	lum = sum * 4*pi*(R*constants.Rsun)**2 * pi
	absolute_mag = -2.5*log10(lum) + 92.0218
	#mag = absolute_mag + 5*log10(dist_mod) - 5
	return (absolute_mag)


def filter_binaries(filterdata,L1,R1,L2,R2):
	sum = 0.0
	const = 2*constants.h*constants.c**2
	T1 = (L1*constants.Lsun / (4*pi*R1**2*constants.Rsun**2) / constants.sigma)**0.25
	T2 = (L2*constants.Lsun / (4*pi*R2**2*constants.Rsun**2) / constants.sigma)**0.25
	bin_min = filterdata[0,0] - (filterdata[1,0] - filterdata[0,0])
	for i in range(len(filterdata)):
		a1 = const / (filterdata[i,0]*constants.Angstrom)**5
		a2_1 = 4*pi*(R1*constants.Rsun)**2 * pi * ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/T1) -1 )**-1.
		a2_2 = 4*pi*(R2*constants.Rsun)**2 * pi * ( exp(constants.c * constants.h / filterdata[i,0]/constants.Angstrom/constants.k/T2) -1 )**-1.
		a3 = (filterdata[i,0] - bin_min)*constants.Angstrom

		sum += a1*a2_1*a3*filterdata[i,1]
		sum += a1*a2_2*a3*filterdata[i,1]
		
		bin_min = filterdata[i,0]
	
	lum = sum 
	absolute_mag = -2.5*log10(lum) + 92.0218
	#mag = absolute_mag + 5*log10(dist_mod) - 5
	return (absolute_mag)



def hrdiag(snap,filter1,filter2,writefile):
	print 'started'
	snap=loadtxt(snap)
	print 'snap done'
	f1=loadtxt(filter1)
	print 'filter1 done'
	f2=loadtxt(filter2)
	print 'filter2 done'
	na=-100
	fwrite=open(writefile,'w')
	print 'open file done'
	fwrite.write("#1.binflag 2.startype0 3.startype1 4.id 5.id0 6.id1 7.m0[MSUN] 8.m1[MSUN] 9.r 10.rad0[RSUN] 11.rad1[RSUN] 12.L0[erg/s] 13.L1[erg/s] 14.B-I 15.B\n")
	print 'header done'

	for i in range(len(snap)):
		if snap[i,7]==0.:
			magb=filter_singles(f1,snap[i,15],snap[i,16])
			magi=filter_singles(f2,snap[i,15],snap[i,16])
			fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,14],na,snap[i,0],na,na,snap[i,1],na,snap[i,2],snap[i,16],na,snap[i,15],na,magb-magi,magb ))
			print (i,snap[i,7])
		elif snap[i,7]==1.:
			magb=filter_binaries(f1,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
			magi=filter_binaries(f2,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
			fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,17],snap[i,18],snap[i,0],snap[i,10],snap[i,11],snap[i,8],snap[i,9],snap[i,2],snap[i,21],snap[i,22],snap[i,19],snap[i,20],magb-magi,magb ))
			print (i,snap[i,7])

	fwrite.close()


def hrdiag_L_T(filestring, snapno):
	"""Give it a snap file and it gives you L and T_eff for all stars
	for binaries the temperature is a luminosity averaged temperature"""
	print 'started' 
	snapfile=filestring+'.snap'+snapno+'.dat.gz'
	snap=loadtxt(snapfile)
	print 'snap done'
	writefile=filestring+'.snap'+snapno+'.hrdiag.dat'
	fwrite=open(writefile,'w')
	na=-100
	fwrite.write("#1.binflag 2.startype0 3.startype1 4.id 5.id0 6.id1 7.m0[MSUN] 8.m1[MSUN] 9.r 10.rad0[RSUN] 11.rad1[RSUN] 12.L0[LSUN] 13.L1[LSUN] 14.T0(K) 15.T1(K) 16.log10(T_eff/K) 17.log10(Leff/LSUN)\n")  #T_eff and T0 are the same for singles. For binaries Teff is L averaged T
	print 'header done'

	for i in range(len(snap)):
		if snap[i,7]==0:
			T0=find_T(snap[i,15],snap[i,16])
			fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %d %g %g\n" % ( snap[i,7],snap[i,14],na,snap[i,0],na,na,snap[i,1],na,snap[i,2],snap[i,16],na,snap[i,15],na,T0,na,log10(T0),log10(snap[i,15])))
			print (i,snap[i,7])
		if snap[i,7]==1:
			T0=find_T(snap[i,19],snap[i,21])
			T1=find_T(snap[i,20],snap[i,22])
			Teff=(snap[i,19]*T0 + snap[i,20]*T1)/(snap[i,19]+snap[i,20])
			fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %d %g %g\n" % ( snap[i,7],snap[i,17],snap[i,18],snap[i,0],snap[i,10],snap[i,11],snap[i,8],snap[i,9],snap[i,2],snap[i,21],snap[i,22],snap[i,19],snap[i,20],T0,T1,log10(Teff),log10(snap[i,19]+snap[i,20])))
			print (i,snap[i,7])
	fwrite.close()


def prune_hrdiag_L_T(filestring, snapno, smoothlength):
	hrdfile=filestring+'.snap'+snapno+'.hrdiag.dat'
	writefile=filestring+'.snap'+snapno+'.hrdiag_pruned.dat'
	f=open(writefile, 'w')

	data=loadtxt(hrdfile)
	for i in range(len(data)):
		if data[i,0]==0:  #single
			if data[i,1]<2 and i%smoothlength==0:  #prune some single MSSs
				for j in range(len(data[i,:])):
					f.write("%f " %(data[i,j],))
				f.write("\n")
			elif data[i,1]>=2:  #all other stars don't prune
				for j in range(len(data[i,:])):
					f.write("%f " %(data[i,j],))
				f.write("\n")
		elif data[i,0]==1:   #binaries
			if data[i,1]<2 and data[i,2]<2 and i%smoothlength==0:   #prune some MS binaries
				for j in range(len(data[i,:])):
					f.write("%f " %(data[i,j],))
				f.write("\n")
			else:						#print everything else
				for j in range(len(data[i,:])):
					f.write("%f " %(data[i,j],))
				f.write("\n")
	f.close()	


def color_comp(snap,writefile,filter1,filter2,filter3,filter4):
	print 'started'
	snap=loadtxt(snap)
	print 'snap done'
	f1=loadtxt(filter1)
	print 'filter1 done'
	f2=loadtxt(filter2)
	print 'filter2 done'
	f3=loadtxt(filter3)
	print 'filter3 done'
	f4=loadtxt(filter4)
	print 'filter4 done'
	na=-100
	fwrite=open(writefile,'w')
	print 'open file done'
	fwrite.write("#1.binflag 2.startype0 3.startype1 4.id 5.id0 6.id1 7.m0[MSUN] 8.m1[MSUN] 9.r 10.rad0[RSUN] 11.rad1[RSUN] 12.L0[erg/s] 13.L1[erg/s] 14.g1(cm/s^2) 15.g2(cm/s^2) 16.T1 (K) 17.T2(K) 18.f435w 19.f555w 20.f606w 21.f814w\n")
	print 'header done'

	for i in range(len(snap)):
		if snap[i,7]==0.:
			g1=constants.G*snap[i,1]*constants.Msun/(snap[i,16]*constants.Rsun)**2.
			T1=find_T(snap[i,15],snap[i,16])
			magb=filter_singles(f1,snap[i,15],snap[i,16])
			magv=filter_singles(f2,snap[i,15],snap[i,16])
			magr=filter_singles(f3,snap[i,15],snap[i,16])
			magi=filter_singles(f4,snap[i,15],snap[i,16])
			fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,14],na,snap[i,0],na,na,snap[i,1],na,snap[i,2],snap[i,16],na,snap[i,15],na,g1,na,T1,na,magb,magv,magr,magi ))
			print (i,snap[i,7],len(snap))
		elif snap[i,7]==1.:
			g1=constants.G*snap[i,8]*constants.Msun/(snap[i,21]*constants.Rsun)**2.
			g2=constants.G*snap[i,9]*constants.Msun/(snap[i,22]*constants.Rsun)**2.
			T1=find_T(snap[i,19],snap[i,21])
			T2=find_T(snap[i,20],snap[i,22])
			magb=filter_binaries(f1,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
			magv=filter_binaries(f2,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
			magr=filter_binaries(f3,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
			magi=filter_binaries(f4,snap[i,19],snap[i,21],snap[i,20],snap[i,22])
			fwrite.write("%d %d %d %ld %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n" % ( snap[i,7],snap[i,17],snap[i,18],snap[i,0],snap[i,10],snap[i,11],snap[i,8],snap[i,9],snap[i,2],snap[i,21],snap[i,22],snap[i,19],snap[i,20],g1,g2,T1,T2,magb,magv,magr,magi ))
			print (i,snap[i,7],len(snap))

	fwrite.close()

def find_rtidal(mc,vg,rg):
	"""takes 
	the cluster mass (mc) in solar mass, 
	the galactic circular velocity (vg) in km/s, 
	and 
	the galactocentric distance (rg) in pc
	returns: the tidal radius of the cluster in pc"""
	r_tidal = (constants.G * mc*constants.Msun / 2 / (vg*constants.km)**2.)**(1./3.) * (rg*constants.PC)**(2./3.) / constants.PC
	return r_tidal


def find_t_vs_rgb(filestring):
	from glob import glob
	files = glob('*.dat.gz')
	dynfile = filestring+'.dyn.dat'
	dyndata = loadtxt(dynfile)
	convfile = filestring+'.conv.sh'
	units = read_units(convfile)
	wfile = filestring+'.t_rgb_bss.dat'
	wfile = open(wfile, 'w')
	wfile.write("#1.t(Myr) 2.SBSS 3.BBSS 4.tot_BSS 5.SRGB 6.BRGB 7.tot_RGB 8.SBSS_c 9.BBSS_c 10.tot_BSS_c\n")
	for i in range(len(files)):
		filestring = files[i].split('.snap')[0]
		snapno = files[i].split('.snap')[1].split('.')[0]
		snapfile = gzip.open(str(files[i]),'rb')
		t_snap = float(snapfile.readline().split('=')[1].split('[')[0])
		m_to = find_MS_turnoff(t_snap*units[7]*1e6)
		m_cut = 1.1*m_to
		snapdata = loadtxt(snapfile)

		for j in range(len(dyndata)):
			if t_snap == dyndata[j,0]:
				rc_snap = dyndata[j,7]
				rh_snap = dyndata[j,20]

		#get bss
		bss_sing_count = 0
		bss_bin_count = 0
		bss_sing_count_core = 0
		bss_bin_count_core = 0
		giants_sing_count = 0
		giants_bin_count = 0
		for j in range(len(snapdata)):
			if snapdata[j,0]!='#' and snapdata[j,0]!='#1:id':
				if snapdata[j,1]>m_cut and snapdata[j,7]==0 and snapdata[j,14]<2.:
					bss_sing_count += 1
					if snapdata[j,2]<=rc_snap:
						bss_sing_count_core += 1
				if snapdata[j,7]==0 and snapdata[j,14]==3:
					giants_sing_count += 1
			
				if snapdata[j,7]!=0. and snapdata[j,8]!=0. and snapdata[j,9]!=0.: 
					if snapdata[j,8]>m_cut and snapdata[j,17]<2.:
						bss_bin_count += 1
						if snapdata[j,2]<=rc_snap:
							bss_bin_count_core += 1
					if snapdata[j,9]>m_cut and snapdata[j,18]<2.:
						bss_bin_count += 1
						if snapdata[j,2]<=rc_snap:
							bss_bin_count_core += 1
					if snapdata[j,17]==3 or snapdata[j,18]==3:
						giants_bin_count += 1


		print t_snap*units[7]/10.**3., 'Gyr'
		tot_bss_snap = bss_sing_count + bss_bin_count 
		tot_bss_core_snap = bss_sing_count_core + bss_bin_count_core 
		tot_giants_snap = giants_sing_count + giants_bin_count 

		dummy = ( t_snap*units[7], bss_sing_count, bss_bin_count, tot_bss_snap, giants_sing_count, giants_bin_count, tot_giants_snap, bss_sing_count_core, bss_bin_count_core, tot_bss_core_snap )
		for j in range(len(dummy)):
			wfile.write("%f " %(dummy[j]))
		wfile.write("\n")

		snapfile.close()

	wfile.close()



def two_point_correlation(a1, a2):
	"""given two arrays of numbers of equal dimension gives the 2-point correlation coefficient.  """

	a1_bar = mean(a1)
	a2_bar = mean(a2)
	a1_sd = std(a1)
	a2_sd = std(a2)

	sum=0
	try:
		if len(a1)!=len(a2):
			raise StopIteration()
		else:
			for i in range(len(a1)):
				sum += (a1[i]-a1_bar)*(a2[i]-a2_bar)/(a1_sd*a2_sd)
	except StopIteration:
		print "bad arrays: arrays have different dimensions"
		pass
		
	cor = float(sum)/float(len(a1))

	return cor

def tidal_mass_loss(filestring):
	"""creates file with tidal mass loss and stellar evolution mass loss
	takes the filestring and creates a file"""

	dynfile=filestring+'.dyn.dat'
	outfile=filestring+'.mass_loss.dat'
	convfile=filestring+'.conv.sh'

	data=loadtxt(dynfile)
	units=read_units(convfile)
	f=open(outfile, 'w')

	sum_Mse=0.
	f.write("#1.t(Myr) 2.M/M(0) 3.M_SE/M(0) 4.dM_SE/M(0)\n")
	for i in range(len(data)):
		sum_Mse += data[i,25]/units[1]
		f.write("%f %f %f %f\n" %(data[i,0]*units[7], data[i,4], data[0,4]-sum_Mse, data[i,25]/units[1]))
	
	f.close()

def mass_vs_tdiss(filestring):

	dynfile=filestring+'.dyn.dat'
	convfile=filestring+'.conv.sh'

	units=read_units(convfile)
	data=loadtxt(dynfile, usecols=(0,4))
	tdiss=[-100, -100, -100, -100, -100, -100, -100, -100]
	count=0
	for j in (0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.05):
		success=0
		for i in range(len(data)):
			if data[i,1]<j and success==0:
				t1=data[i-1,0]
				t2=data[i,0]
				m1=data[i-1,1]
				m2=data[i,1]
				t=t1 - (j-m1)*(t2-t1)/(m2-m1)
				tdiss[count]=t*units[7]
				success=1
				count += 1
	print "%f %f %f %f %f %f %f %f %f" %(units[1], tdiss[0], tdiss[1], tdiss[2], tdiss[3], tdiss[4], tdiss[5], tdiss[6], tdiss[7])


def bss_rdist_cumulative(hrdiagfile, mcut, nbin):
	"""take the file created by hrdiag_L_T and create a cumulative 
	radial distribution normalized with the RGBs"""
	data = loadtxt(hrdiagfile)
	bss_sing_r = []
	bss_bin_r = []
	bss_all_r = []
	rgb_r = []
	for i in range(len(data)):
		if data[i,0]==0:  #singles
			if data[i,6]>mcut and data[i,1]<2:  #BSSs
				bss_all_r.append(data[i,8])
				bss_sing_r.append(data[i,8])
		if data[i,0]==1: #binaries
			if data[i,6]>mcut and data[i,1]<2:  #BSSs
				bss_all_r.append(data[i,8])
				bss_bin_r.append(data[i,8])
			if data[i,7]>mcut and data[i,2]<2:  #BSSs
				bss_all_r.append(data[i,8])
				bss_bin_r.append(data[i,8])
		if data[i,1]==3 or data[i,2]==3:
			rgb_r.append(data[i,8])

	maxr_g, minr_g = max(rgb_r), min(rgb_r)
	maxr_bss = max(bss_all_r)
	if maxr_g >= maxr_bss:
		maxr = maxr_g
	else:
		maxr = maxr_bss
	minr = minr_g

	print len(bss_all_r), len(bss_sing_r), len(bss_bin_r), len(rgb_r)
	
	cum_bss_all = cumulate_array(bss_all_r, minr, maxr, nbin)
	cum_bss_sing = cumulate_array(bss_sing_r, minr, maxr, nbin)
	cum_bss_bin = cumulate_array(bss_bin_r, minr, maxr, nbin)
	cum_rgb = cumulate_array(rgb_r, minr, maxr, nbin)

	binsize = (maxr-minr)/nbin
	for i in range(len(cum_rgb)):
		print cum_rgb[i][0], cum_bss_all[i][1], cum_bss_sing[i][1], cum_bss_bin[i][1], cum_rgb[i][1]
		

def cumulate_array(a, minvalue, maxvalue, nbin):
	"""Cumulates from array a from minvalue with steps of 
	binsize=(maxvalue-minvalue)/nbin and prints out the values"""
	binsize = (maxvalue-minvalue)/nbin
	c_array = []
	for i in range(nbin+1):
		c = 0
		a_temp = minvalue+(i+1)*binsize
		for j in range(len(a)):
			if a[j] < a_temp:
				c += 1
		c_array.append([a_temp, c])
	return c_array
		

def nbss_mcut(hrdiagfile, mcut):
	data = loadtxt(hrdiagfile)
	for j in range(100):
		n, ns, nbin = 0, 0, 0
		for i in range(len(data)):
        		if data[i,0]==0 and data[i,1]<2 and data[i,6]>mcut:
            			n += 1
				ns += 1
        		elif data[i,1]==1:
            			if data[i,1]<2 and data[i,6]>mcut:
                			n += 1
					nbin += 1
            			elif data[i,2]<2 and data[i,7]>mcut:
                			n += 1
					nbin += 1
    		print mcut, n, ns, nbin
    		mcut += 0.01


def make_3D(r):
	"""takes some quantity, like the radial position of a star and creates the 3D position 
	by randomly creating the angles \theta and \phi.  Assumes spherical symmetry.  """
	import random
	costheta = random.uniform(-1, 1)
	sintheta = (1-costheta**2.)**0.5
	phi = random.uniform(0, 4*pi)
	rz = r*costheta
	rx = r*sintheta*cos(phi)
	ry = r*sintheta*sin(phi)

	return (rx, ry, rz)

def surface_density_profile_L(filestring, snapno, seedy, proj, nbins, bintype, writefile):
	"""takes the filestring for the run, snapno of interest, seed the random generator, 
	projections like (0, 1, 2) gives 3D projection, (0, 1) gives 2d suppressing z
	no of bins wanted and the binning type, 0:equal member in each bin 1:bins are equidistant 
	in log10(r2D)"""
	#first learn the physical units
	units = read_units(filestring)
	lpc = units[0]['l_pc']

	writefile=open(writefile, 'w')

	#read the snapfile
	snapfile = filestring+'.snap'+snapno+'.dat.gz'
	colnos = (2, 7, 15, 19, 20, 14, 17, 18)
	data = loadtxt(snapfile, usecols=colnos)

	dtype = [('r', float), ('x', float), ('y', float), ('z', float), ('r2D', float), ('logr2D', float), ('binflag', int), ('L', float), ('type0', int), ('bintype0', int), ('bintype1', int)]
	a = []
	import random
	random.seed(seedy)
	for i in range(len(data)):
		rs = make_3D(data[i,0])
		r2d = 0
		for j in proj:
			r2d += rs[j]**2. 
		r2d = r2d**0.5

		if data[i,1]==0:		#singles
			a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,2], data[i,5], data[i,6], data[i,7]))
		elif data[i,1]==1:		#binaries
			a.append((data[i,0]*lpc, rs[0]*lpc, rs[1]*lpc, rs[2]*lpc, r2d*lpc, log10(r2d*lpc), data[i,1], data[i,3]+data[i,4], data[i,5], data[i,6], data[i,7]))

	a = array(a, dtype=dtype)
	a = sort(a, order='logr2D')
	print a[:]['logr2D']
	#for equal member binning
	if bintype==0:
		binsize = len(a)/nbins  
		for i in range(binsize, len(a), binsize):
			if a[i]['type0']!=3 and a[i]['bintyppe0']!=3 and a[i]['bintype1']!=3:
				lsum = sum(a[i-binsize:i]['L'])
			#print lsum, i, binsize, i-binsize
			#for j in range(i-binsize, i):
			#	print a[j]['L']
			area = 4.*pi*(a[i]['r2D']**2. - a[i-binsize]['r2D']**2.)
			lsum, lsum_err = lsum/area, lsum/float(binsize**0.5)/area
			n2d, n2d_err = binsize/area, float(binsize**0.5)/area
			r2d_av = (a[i]['r2D']+a[i-binsize]['r2D'])/2.
			writefile.write("%f %f %f %f %f\n" %(r2d_av, lsum, lsum_err, n2d, n2d_err))
			#print r2d_av, lsum, lsum_err, n2d, n2d_err


	#for binning equidistant in logr2D
	elif bintype==1:
		lbinsize = abs((a[-1]['logr2D']-a[0]['logr2D'])/nbins)
		n2d_prev=0	
		for i in range(1, nbins):
			lsum, lsum_err, rsum, n2d, n2d_err = 0., 0., 0., 0., 0.
			lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err = 0., 0., 0., 0.
			lr_upperbound, lr_lowerbound = a[0]['logr2D']+i*lbinsize, a[0]['logr2D']+(i-1)*lbinsize
			lr2D_av = (lr_upperbound+lr_lowerbound)/2.
			area = pi*( (10.**lr_upperbound)**2. - (10**lr_lowerbound)**2. )
			try:
				for j in range(n2d_prev, len(a)):
					if a[j]['logr2D']<lr_upperbound and a[j]['logr2D']>=lr_lowerbound:
						lsum += a[j]['L'] 
						n2d += 1
						if a[j]['L']<=20.:
							lsum_cut += a[j]['L']
							n2d_cut += 1
					else:
						raise StopIteration()
			except StopIteration:
				print n2d, n2d_prev, 'got the values'
			n2d_prev += n2d
			#print lr2D_av, n2d
			if n2d>1:	#just to avoid the points with 100% poisson error
				lsum, lsum_err = lsum/area, lsum/float(n2d)**(0.5)/area
				lsum_cut, lsum_cut_err = lsum_cut/area, lsum_cut/float(n2d_cut)**(0.5)/area
				n2d, n2d_err = n2d/area, float(n2d)**0.5/area
				n2d_cut, n2d_cut_err = n2d_cut/area, float(n2d_cut)**0.5/area
				writefile.write("%f %f %f %f %f %f %f %f %f\n" %(10**lr2D_av, lsum, lsum_err, n2d, n2d_err, lsum_cut, lsum_cut_err, n2d_cut, n2d_cut_err))
				#print 10**lr2D_av, lsum, lsum_err, n2d, n2d_err


	writefile.close()





def king_fit(r, n, rc_guess, rt_guess):
	"""fits a single mass King profile given the r array, the n array and 
	the initial guesses for rc and rt"""
	rc = []
	rt = []
	for i in range(100):
		rc.append(rc_guess*(1.-50*0.01)+0.01*rc_guess*i)
		rt.append(rt_guess*(1.-50*0.01)+0.01*rt_guess*i)
	a = []
	dtype = [('rc', float), ('rt', float), ('norm', float), ('sqsum', float)]
	for i in range(len(rc)):
		for j in range(len(rt)):
			norm = n[len(n)/2]/( ( 1/ (1+ (r[len(r)/2]/rc[i])**2.)**0.5 ) - ( 1/ (1+ (rt[j]/rc[i])**2.)**0.5 ) )**2.
			sqsum = 0.
			for k in range(len(n)):
				n_k = norm*( 1/( 1+ (r[k]/rc[i])**2. )**0.5 - 1/( 1+ (rt[j]/rc[i])**2. )**0.5 )**2.
				sqsum += (n[k] - n_k)**2.
			a.append( (rc[i], rt[j], norm, sqsum) ) 
	a = array(a, dtype=dtype)
	a = sort(a, order='sqsum')

	return a

def smooth_data(filename, smoothlength, writefile):
	"""smoothens data with the smoothing length smoothlength and writes in writefile"""
	data=loadtxt(filename)
	writefile=open(writefile, 'w')
	for i in range(smoothlength, len(data)):
    		if i%smoothlength==0:
        		for j in range(len(data[i,:])):
            			temp=mean(data[i-smoothlength:i,j])
            			writefile.write("%f " %(temp))
        		writefile.write("\n")
	writefile.close()


def get_roche_r(q, a):
	R = a*0.49*q**(2./3.) / (0.6*q**(2./3.) + log(1+q**(1./3.)))

	return R

def get_period(a,m1,m2): #all in CGS
	p = (4.*pi**2.*(a**3.)/constants.G/(m1+m2))**0.5 /3600./24.	#p in days

	return p

def plot_P_e(bsfilename):
	data = loadtxt(bsfilename)
	p = []
	for i in range(len(data)):
		p.append(get_period(data[i,12]*constants.AU, data[i,8]*constants.Msun, data[i,9]*constants.Msun))
	
	
	#separate
	msms=[]
	msg=[]
	for i in range(len(data)):
		if data[i,17]<2 and data[i,18]<2:
			msms.append((p[i],data[i,13]))
		else:
			msg.append((p[i],data[i,13]))




	dtype = [('p', float), ('e', float)] 
        msms = array(msms, dtype=dtype)
	msg = array(msg, dtype=dtype)

	print msms[:]['p'], msms[:]['e']
	#now plot
	import gracePlot
	gpl=gracePlot.gracePlot()
	gpl.hold()
	gpl.plot(msms[:]['p'], msms[:]['e'], symbols=1)
	gpl.plot(msg[:]['p'], msg[:]['e'], symbols=1)

