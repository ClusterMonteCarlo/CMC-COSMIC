noinst_LTLIBRARIES = libbse.la
# Since this is just a convenience library, we set the -static flag to avoid 
# double compilation.
libbse_la_FFLAGS = -static
libbse_la_SOURCES= comenv.f corerd.f deltat.f dgcore.f evolv1.f evolv2.f gntage.f hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f rl.f star.f zcnsts.f zdata.h zfuncs.f

if USE_TAUS
  libbse_la_SOURCES+=int64.f tausworth.f taus113-ran3.f
else
  libbse_la_SOURCES+=ran3.f
endif
