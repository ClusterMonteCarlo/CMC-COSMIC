ACLOCAL_AMFLAGS = -I m4 --install
AUTOMAKE_OPTIONS = foreign
SUBDIRS = libs bse_wrap fewbody-0.24 

bin_PROGRAMS = cmc
cmc_SOURCES = cmc_bhlosscone.c cmc_binbin.c cmc_binsingle.c cmc.c cmc_core.c \
	      cmc_core.h cmc_dynamics.c cmc_dynamics_helper.c \
	      cmc_evolution_thr.c cmc_fits.c cmc_funcs.c cmc.h cmc_init.c \
	      cmc_io.c cmc_nr.c cmc_orbit.c \
	      cmc_remove_star.c cmc_search_grid.c cmc_sort.c cmc_sscollision.c\
	      cmc_stellar_evolution.c cmc_trace.c cmc_utils.c cmc_vars.h \
	      common/fitslib.h common/gensearch.h common/gensort.h \
              common/taus113-v2.h cmc_bse_utils.c cmc_bse_utils.h
cmc_CPPFLAGS = -DCMCVERSION=\"$(CMCVERSION)\" -DCMCDATE=\"$(CMCDATE)\"
cmc_LDADD = libs/libsupport.la fewbody-0.24/libfewbody.la bse_wrap/libbsewrap.la
cmc_LDFLAGS = $(LIBTOOL_LDFLAGS)

