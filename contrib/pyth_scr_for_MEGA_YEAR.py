#!/usr/bin/python2.2
# a python script to calculate MEGA_YEAR for cmc, 
# for collision runs, which are Plummer models, and whose rh(0) are given in pc
#             trh(0) in FP units   <-- 0.0931 for Plummer model
# MEGA_YEAR = ------------------
#               trh(0) in Myr

from math import *

G_cgs = 6.672e-8
N = 5e5
rh_pc = 100.0
rh_cgs = rh_pc*3.086e18
mave_Msun = 0.6893375799
m_Msun = N * mave_Msun
m_cgs = m_Msun * 1.989e33
gamma = 0.04

trh_cgs = 0.138*N/log(gamma*N)*sqrt(rh_cgs**3/(G_cgs*m_cgs))
trh_myr = trh_cgs / (365*24*3600*1e6)

MY = 0.0931/trh_myr

print 'trh(0) =', trh_myr, 'Myr'
print 'MEGA_YEAR =', MY
