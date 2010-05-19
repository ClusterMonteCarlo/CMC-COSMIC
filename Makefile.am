ACLOCAL_AMFLAGS = -I m4 
AUTOMAKE_OPTIONS = foreign
SUBDIRS = libs bse_wrap fewbody-0.24 

bin_PROGRAMS = cmc cmc_mkplummer cmc_mkking cmc_setimf cmc_setunits cmc_setstellar cmc_addbinaries

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

cmc_mkplummer_SOURCES = utils/mkplummer/mkplummer.c
cmc_mkplummer_LDADD = libs/libsupport.la
cmc_mkplummer_LDFLAGS = $(LIBTOOL_LDFLAGS)

cmc_mkking_SOURCES = utils/mkking/mkking.c
cmc_mkking_LDADD = libs/libsupport.la
cmc_mkking_LDFLAGS = $(LIBTOOL_LDFLAGS)

cmc_setimf_SOURCES = utils/setimf/setimf.c utils/setimf/genimf.c utils/setimf/genimf.h
cmc_setimf_LDADD = libs/libsupport.la
cmc_setimf_LDFLAGS = $(LIBTOOL_LDFLAGS)

cmc_setunits_SOURCES = utils/setunits/setunits.c
cmc_setunits_LDADD = libs/libsupport.la
cmc_setunits_LDFLAGS = $(LIBTOOL_LDFLAGS)

cmc_setstellar_SOURCES = utils/setstellar/setstellar.c
cmc_setstellar_LDADD = libs/libsupport.la bse_wrap/libbsewrap.la
cmc_setstellar_LDFLAGS = $(LIBTOOL_LDFLAGS)

cmc_addbinaries_SOURCES = utils/addbinaries/addbinaries.c
cmc_addbinaries_LDADD = libs/libsupport.la bse_wrap/libbsewrap.la
cmc_addbinaries_LDFLAGS = $(LIBTOOL_LDFLAGS)

