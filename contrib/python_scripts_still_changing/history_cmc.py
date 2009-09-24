#!/opt/local/bin//python
from numpy import *
import gzip
import scripts
import scripts1
import scripts2
import scripts3

def history_maker(ids,positions,file_string,binary):
	"""creates the full history of interesting stars
	ids: array containing the ids of interesting stars
	binintfile: binint file
	collfile: collision log file
	mergefile: merger file
	filestring: the filestring from cmc
	binary: whether there was any primordial binaries or not -> 0: no binary, 1: yes binary"""
	binintfile=file_string+'.binint.log'
	collfile=file_string+'.collision.log'
	mergefile=file_string+'.semergedisrupt.log'
	history_interactions={}
	#make the collision dictionary by reading and sorting the collision log file
	coll=scripts1.collision(collfile)
	#make the pure stellar evolution dictionary
	se=scripts2.bse_int(mergefile)

	if binary:
		#make the binint dictionary by reading from the binint log file
		binint=scripts3.read_binint(binintfile) 


	#now go through the ids of interest and populate the history file with every kind of interactions
	for i in range(len(ids)):
		if binary:
			#read the binint history of the id in this dictionary
			binint_id=scripts3.binary_interaction(binint,ids[i])
			bininteract={'binint': binint_id,
				'binint_coll': {}
				}
			#cross_correlate with collisions: if binint resulted in a merger this gives if any of the stars had any collisions before this
			if len(binint_id.keys())>0:
				for j in binint_id.keys():
					print j
					if binint_id[j]['merge']==1:
						for k in range(len(binint_id[j]['mergeids'])):
							cross_id=int(binint_id[j]['mergeids'][k])
							print cross_id
							binint_coll=scripts1.call_collision_tree(coll,cross_id)
							bininteract['binint_coll'][k]=binint_coll
		else:
			bininteract = {}
			
		collision=scripts1.call_collision_tree(coll,ids[i])

		if se.has_key(ids[i]):
			se_id=se[ids[i]]
		else:
			se_id={}

		history_interactions[ids[i]]={'binint': bininteract,
				'coll': collision,
				'se': se_id,
				'position': positions[i]
				}
	
	return history_interactions



def branching_ratio(history, binary):
	all=[]
	binint=[]
	pure_binint=[]
	coll=[]
	pure_coll=[]
	merger=[]
	pure_merger=[]
	binint_coll=[]
	binint_merger=[]
	
	for i in history.keys():
		all+=[history[i]['position']]
		if binary and len(history[i]['binint']['binint'].keys())>0:
		#if len(history[i]['binint']['binint'].keys())>0:
			binint+=[history[i]['position']]
		
		if binary==1 and len(history[i]['binint']['binint'].keys())>0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())==0:
			pure_binint+=[history[i]['position']]

		if len(history[i]['coll'].keys())>0:
			coll+=[history[i]['position']]

		if binary==1 and len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())>0 and len(history[i]['se'].keys())==0:
			pure_coll+=[history[i]['position']]
			
		if binary==1 and len(history[i]['se'].keys())>0:
			merger+=[history[i]['position']]

		if binary==1 and len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())>0:
			pure_merger+=[history[i]['position']]

		if binary==1 and len(history[i]['coll'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_coll+=[history[i]['position']]

		if binary==1 and len(history[i]['se'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_merger+=[history[i]['position']]

	if binary==0:
		pure_coll = coll

	branching_ratio={'binint': float(len(binint))/float(len(history.keys())),
			'pure_binint': float(len(pure_binint))/float(len(history.keys())),
			'coll': float(len(coll))/float(len(history.keys())),
			'pure_coll': float(len(pure_coll))/float(len(history.keys())),
			'merger': float(len(merger))/float(len(history.keys())),
			'pure_merger': float(len(pure_merger))/float(len(history.keys()))
			}

	return branching_ratio



def branching_ratio_plot(history,r_reference):
	all=[]
	binint=[]
	pure_binint=[]
	coll=[]
	pure_coll=[]
	merger=[]
	pure_merger=[]
	binint_coll=[]
	binint_merger=[]
	
	for i in history.keys():
		all+=[history[i]['position']]
		if len(history[i]['binint']['binint'].keys())>0:
			binint+=[history[i]['position']]
		
		if len(history[i]['binint']['binint'].keys())>0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())==0:
			pure_binint+=[history[i]['position']]

		if len(history[i]['coll'].keys())>0:
			coll+=[history[i]['position']]

		if len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())>0 and len(history[i]['se'].keys())==0:
			pure_coll+=[history[i]['position']]
			
		if len(history[i]['se'].keys())>0:
			merger+=[history[i]['position']]

		if len(history[i]['binint']['binint'].keys())==0 and len(history[i]['coll'].keys())==0 and len(history[i]['se'].keys())>0:
			pure_merger+=[history[i]['position']]

		if len(history[i]['coll'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_coll+=[history[i]['position']]

		if len(history[i]['se'].keys())>0 and len(history[i]['binint']['binint'].keys())>0:
			binint_merger+=[history[i]['position']]

	branching_ratio={'binint': float(len(binint))/float(len(history.keys())),
			'pure_binint': float(len(pure_binint))/float(len(history.keys())),
			'coll': float(len(coll))/float(len(history.keys())),
			'pure_coll': float(len(pure_coll))/float(len(history.keys())),
			'merger': float(len(merger))/float(len(history.keys())),
			'pure_merger': float(len(pure_merger))/float(len(history.keys()))
			}
	branching={'all': all,
		'binint': binint,
		'pure_binint': pure_binint,
		'coll': coll,
		'pure_coll': pure_coll,
		'merger': merger,
		'pure_merger': pure_merger,
		'branching_ratio': branching_ratio,
		'r_dist': r_reference
		}


	#now make plot	
	nobins=20
	min=0
	max=5

	r_dist_hist, lower_edges = histogram(branching['r_dist'],bins=nobins,range=(min,max),normed=False)
	binwidth = lower_edges[-1]/len(lower_edges)
	lower_edges=lower_edges[:(len(lower_edges)-1)]
	r_dist_cumhist = r_dist_hist.cumsum()

	all_hist, lower_edges = histogram(branching['all'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_all_hist = zeros(nobins)
	for i in range(nobins):
		norm_all_hist[i] = float(all_hist[i]) / float(r_dist_hist[i])

	all_cumhist = all_hist.cumsum()
	norm_all_cum = zeros(nobins)
	for i in range(nobins):
		norm_all_cum[i] = float(all_cumhist[i]) / float(r_dist_cumhist[i])
	
	binint_hist, lower_edges = histogram(branching['binint'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_binint_hist = zeros(nobins)
	for i in range(nobins):
		norm_binint_hist[i] = float(binint_hist[i]) / float(r_dist_hist[i])

	binint_cumhist = binint_hist.cumsum()
	norm_binint_cum = zeros(nobins)
	for i in range(nobins):
		norm_binint_cum[i] = float(binint_cumhist[i]) / float(r_dist_cumhist[i])
	
	pure_binint_hist, lower_edges = histogram(branching['pure_binint'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]
	
	norm_pure_binint_hist=zeros(nobins)
	for i in range(nobins):
		norm_pure_binint_hist[i] = float(pure_binint_hist[i]) / float(r_dist_hist[i])

	pure_binint_cumhist = pure_binint_hist.cumsum()
	norm_pure_binint_cum=zeros(nobins)
	for i in range(nobins):
		norm_pure_binint_cum[i] = float(pure_binint_cumhist[i]) / float(r_dist_cumhist[i])

	coll_hist, lower_edges = histogram(branching['coll'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_coll_hist = zeros(nobins)
	for i in range(nobins):
		norm_coll_hist[i] = float(coll_hist[i]) / float(r_dist_hist[i])

	coll_cumhist = coll_hist.cumsum()
	norm_coll_cum = zeros(nobins)
	for i in range(nobins):
		norm_coll_cum[i] = float(coll_cumhist[i]) / float(r_dist_cumhist[i])
	
	pure_coll_hist, lower_edges = histogram(branching['pure_coll'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_pure_coll_hist = zeros(nobins)
	for i in range(nobins):
		norm_pure_coll_hist[i] = float(pure_coll_hist[i]) / float(r_dist_hist[i])

	pure_coll_cumhist = pure_coll_hist.cumsum()
	norm_pure_coll_cum = zeros(nobins)
	for i in range(nobins):
		norm_pure_coll_cum[i] = float(pure_coll_cumhist[i]) / float(r_dist_cumhist[i])
	
	merger_hist, lower_edges = histogram(branching['merger'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_merger_hist = zeros(nobins)
	for i in range(nobins):
		norm_merger_hist[i] = float(merger_hist[i]) / float(r_dist_hist[i])

	merger_cumhist = merger_hist.cumsum()
	norm_merger_cum = zeros(nobins)
	for i in range(nobins):
		norm_merger_cum[i] = float(merger_cumhist[i]) / float(r_dist_cumhist[i])
	
	pure_merger_hist, lower_edges = histogram(branching['pure_merger'],bins=nobins,range=(min,max),normed=False)
	lower_edges=lower_edges[:(len(lower_edges)-1)]

	norm_pure_merger_hist = zeros(nobins)
	for i in range(nobins):
		norm_pure_merger_hist[i] = float(pure_merger_hist[i]) / float(r_dist_hist[i])

	pure_merger_cumhist = pure_merger_hist.cumsum()
	norm_pure_merger_cum = zeros(nobins)
	for i in range(nobins):
		norm_pure_merger_cum[i] = float(pure_merger_cumhist[i]) / float(r_dist_cumhist[i])


	import gracePlot
	gpl=gracePlot.gracePlot()
	gpl.hold()

	gpl.plot(lower_edges+0.5*binwidth, norm_all_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_binint_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_pure_binint_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_coll_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_pure_coll_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_merger_hist, symbols=1, styles=1)
	gpl.plot(lower_edges+0.5*binwidth, norm_pure_merger_hist, symbols=1, styles=1)
	
	return branching['branching_ratio']


def find_giants(snapfile):
	snapdata=loadtxt(snapfile)
	r_giants=[]
	for i in range(len(snapdata)):
		#single RGBs
		if snapdata[i,7]==0 and snapdata[i,14]==4:
			r_giants.append(snapdata[i,2])
			print len(r_giants)
		elif snapdata[i,7]==1:
			if snapdata[i,17]==4 or snapdata[i,18]==4:
				r_giants.append(snapdata[i,2])
				print len(r_giants)
	
	return r_giants

def find_r_dist(snapfile):
	snapdata=loadtxt(snapfile)
	r=[]
	for i in range(len(snapdata)):
		r.append(snapdata[i,2])

	return r

			

		
		
