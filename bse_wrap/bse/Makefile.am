noinst_LTLIBRARIES = libbse.la
# Since this is just a convenience library, we set the -static flag to avoid 
# double compilation.
libbse_la_FFLAGS = -static
libbse_la_SOURCES= ./COSMIC/cosmic/src/comenv.f ./COSMIC/cosmic/src/corerd.f ./COSMIC/cosmic/src/deltat.f ./COSMIC/cosmic/src/dgcore.f ./COSMIC/cosmic/src/evolv1.f ./COSMIC/cosmic/src/evolv2.f ./COSMIC/cosmic/src/gntage.f ./COSMIC/cosmic/src/hrdiag.f ./COSMIC/cosmic/src/instar.f ./COSMIC/cosmic/src/kick.f ./COSMIC/cosmic/src/mix.f ./COSMIC/cosmic/src/mlwind.f ./COSMIC/cosmic/src/mrenv.f ./COSMIC/cosmic/src/rl.f ./COSMIC/cosmic/src/star.f ./COSMIC/cosmic/src/zcnsts.f ./COSMIC/cosmic/src/zdata.h ./COSMIC/cosmic/src/zfuncs.f ./COSMIC/cosmic/src/concatkstars.f ./COSMIC/cosmic/src/bpp_array.f

if USE_TAUS
  libbse_la_SOURCES+=./COSMIC/cosmic/src/int64.f ./COSMIC/cosmic/src/tausworth.f ./COSMIC/cosmic/src/taus113-ran3.f
else
  libbse_la_SOURCES+=./COSMIC/cosmic/src/ran3.f
endif
