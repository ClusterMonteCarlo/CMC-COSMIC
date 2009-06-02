#!/opt/local/bin//python

from numpy import *
import gzip
import math
import constants

def read_units(s):
	"""reads from the supplied conv file and stores the physical units"""
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
	return (massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr)


def find_correlations(string,snapno):
	snap = string+'.snap'+snapno+'.dat.gz'
	dyn = string+'.dyn.dat'
	conv = string+'.conv.sh'
	out = string+'.out.dat'
	bin = string+'.bin.dat'

	units = read_units(conv)
	t_myr = find_t_myr(snap,conv)
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
			m_msun = dyndata[i,4]*units[1]
			rc_pc = dyndata[i,7]*units[5]
			rh_pc = dyndata[i,20]*units[5]
			rho_msun_pc3 = dyndata[i,21]*units[1]/units[5]**3.
			#lo behold: time unit in central v calculation is in nbody units
			v0_rms_km_s = dyndata[i,23]*units[4]/units[8]/100000.
			print (m_msun,rc_pc,rh_pc,rho_msun_pc3,v0_rms_km_s,t_myr,m_to)
	
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
				f_bc_array.append(bindata[i,11])
				f_b_array.append(bindata[i,12])
				rhb_pc_array.append(bindata[i,6]*units[5])
				rhs_pc_array.append(bindata[i,5]*units[5])
			f_bc=mean(f_bc_array)
			f_b=mean(f_b_array)
			rhb_pc=mean(rhb_pc_array)
			rhs_pc=mean(rhs_pc_array)
	except IOError:
		f_bc = 0.0
		f_b = 0.0
		rhb_pc = 0.0
		rhs_pc = rh_pc	

	BSS_count = find_BSS(string,snapno)
	print "got BSS"
	giants_count = find_giants(string,snapno)
	print "got giants"
	total_giants = (giants_count[0]+giants_count[1])

	outfile=open(out,'r')
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
				
	print "%f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %f %d %d %d %f %f %f\n" %(m_msun,rc_pc,rh_pc,rc_pc/rh_pc,rho_msun_pc3,mean_gamma,v0_rms_km_s,t_myr,m_to,f_bc,f_b,rhb_pc,rhs_pc,BSS_count[0],BSS_count[1],BSS_count[0]+BSS_count[1],rhb_pc/rhs_pc,giants_count[0], giants_count[1], giants_count[0]+giants_count[1], float(BSS_count[0])/float(total_giants), float(BSS_count[1])/float(total_giants), float(BSS_count[0]+BSS_count[1])/float(total_giants))

	return (m_msun,rc_pc,rh_pc,rc_pc/rh_pc,rho_msun_pc3,mean_gamma,v0_rms_km_s,t_myr,m_to,f_bc,f_b,rhb_pc,rhs_pc,BSS_count[0],BSS_count[1],BSS_count[0]+BSS_count[1], rhb_pc/rhs_pc,giants_count[0], giants_count[1], giants_count[0]+giants_count[1], float(BSS_count[0])/float(total_giants), float(BSS_count[1])/float(total_giants), float(BSS_count[0]+BSS_count[1])/float(total_giants)) 

	
	
	

def find_t_myr(s1,s2):
	"""goes in the given snapshot and finds the physical time corresponding to that snap"""
	f=gzip.open(s1,'rb')
	line=f.readline()
	a=line.split()
	b=a[1]
	c=b.split('=')
	t=float(c[1])
	d=read_units(s2)
	t_myr=t*d[7]
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
	wfile1=file_string+'.snap'+snapno+'single_giants.dat'
	wfile2=file_string+'.snap'+snapno+'binary_giants.dat'

	t_yr = (find_t_myr(snapfile,convfile))*10**6

	f=gzip.open(snapfile,'rb')
	f1 = open(wfile1,'w')
	f2 = open(wfile2,'w')
	giants_sing_count = 0
	giants_bin_count = 0

	#print header
	f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
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
						f2.write("%s " % (a[i]))
					giants_bin_count += 1
					f2.write("\n")
				elif float(a[18])==3.:
					for i in range(len(a)):
						f2.write("%s " % (a[i]))
					giants_bin_count += 1
					f2.write("\n")
	f1.close()
	f2.close()
	return (giants_sing_count, giants_bin_count,t_yr/10**6)

	

def find_BSS(file_string,snapno):
	"""takes the snapshot file, the conv.sh file, the file name to write the BSS in singles and the file name to write the BSS in binaries at that snapshot.  Returns the number of BSS in single and the number of BSS in binaries."""
	snapfile=file_string+'.snap'+snapno+'.dat.gz'
	convfile=file_string+'.conv.sh'
	wfile1=file_string+'.snap'+snapno+'single_BSS.dat'
	wfile2=file_string+'.snap'+snapno+'binary_BSS.dat'

	t_yr = (find_t_myr(snapfile,convfile))*10**6
	m = find_MS_turnoff(float(t_yr))
	m_cut = (1.+0.1)*m
	print "looking at a snap at t = %f Gyr, MS turnoff is %f MSun" % (t_yr/10**9, m)

	f=gzip.open(snapfile,'rb')
	f1 = open(wfile1,'w')
	f2 = open(wfile2,'w')
	bss_sing_count = 0
	bss_bin_count = 0

	#print header
	f1.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	f2.write("#1:id #2:m[MSUN] #3:r #4:vr #5:vt #6:E #7:J #8:binflag #9:m0[MSUN] #10:m1[MSUN] #11:id0 #12:id1 #13:a[AU] #14:e #15:startype #16:luminosity[LSUN] #17:radius[RSUN]  #18:bin_startype0 #19:bin_startype1 #20:bin_star_lum0[LSUN] #21:bin_star_lum1[LSUN] #22:bin_star_radius0[RSUN] #23:bin_star_radius1[RSUN]\n")
	for line in f:	
		a=line.split()
		if a[0]!='#' and a[0]!='#1:id':
			if float(a[1])>m_cut and a[7]=='0' and float(a[14])<2.:
				for i in range(len(a)):
					f1.write("%s " % (a[i]))
				bss_sing_count += 1
				f1.write("\n")
				
			
			if a[7]!='0' and a[8]!='0' and a[9]!='0': 
				if (float(a[8]))>m_cut and float(a[17])<2.:
					for i in range(len(a)):
						f2.write("%s " % (a[i]))
					bss_bin_count += 1
					f2.write("\n")
				if (float(a[9]))>m_cut and float(a[18])<2.:
					for i in range(len(a)):
						f2.write("%s " % (a[i]))
					bss_bin_count += 1
					f2.write("\n")
	f1.close()
	f2.close()
	return (bss_sing_count, bss_bin_count,t_yr/10**6)




				
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
	hd_both=histogram(both_roverrc,5)
	hd_sing=histogram(sing_roverrc,5)
	hd_bin=histogram(bin_roverrc,5)

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
	"""L in LSUN and R in RSUN"""
	T= (L*constants.Lsun / (4*pi*(R*constants.Rsun)**2) / constants.sigma)**0.25
	return (T)


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


def hrdiag_L_T(snap,writefile):
	"""Give it a snap file and it gives you L and T_eff for all stars
	for binaries the temperature is a luminosity averaged temperature"""
	print 'started' 
	snap=loadtxt(snap)
	print 'snap done'
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


	
	




	
	
	

					
				
	
		
