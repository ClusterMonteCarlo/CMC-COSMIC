ACLOCAL_AMFLAGS = -I m4 
AUTOMAKE_OPTIONS = foreign
SUBDIRS = libs bse_wrap fewbody-0.24 

.cu.cpp:
	${NVCC} --cuda -arch ${CUDAARCH} ${CPPFLAGS} -o $@ $<

bin_PROGRAMS = cmc

if BUILD_UTILS
bin_PROGRAMS+= cmc_mkplummer cmc_mkking cmc_setimf cmc_setunits cmc_setstellar cmc_addbinaries
endif

cmc_SOURCES= cmc.c
cmc_CPPFLAGS = -DCMCVERSION=\"$(CMCVERSION)\" -DCMCDATE=\"$(CMCDATE)\"
cmc_LDADD = libcmc.la libs/libsupport.la fewbody-0.24/libfewbody.la bse_wrap/libbsewrap.la
cmc_LDFLAGS = $(LIBTOOL_LDFLAGS)

if USE_CUDA
cmc_SOURCES += cuda/cmc_cuda.cpp cuda/cmc_cuda.h
cmc_CPPFLAGS += $(CUDA_CONFIG_INCLUDES) -DUSE_CUDA
cmc_LDFLAGS += $(CUDA_CONFIG_LIBRARIES)
endif

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

noinst_LTLIBRARIES = libcmc.la
# Since this is just a convenience library, we set the -static flag to avoid
# double compilation.
libcmc_la_CFLAGS = -static
libcmc_la_CPPFLAGS = -DCMCVERSION=\"$(CMCVERSION)\" -DCMCDATE=\"$(CMCDATE)\"
libcmc_la_SOURCES = cmc_bhlosscone.c cmc_binbin.c cmc_binsingle.c cmc_core.c \
	      cmc_core.h cmc_dynamics.c cmc_dynamics_helper.c \
	      cmc_evolution_thr.c cmc_fits.c cmc_funcs.c cmc.h cmc_init.c \
	      cmc_io.c cmc_nr.c cmc_orbit.c \
	      cmc_remove_star.c cmc_search_grid.c cmc_sort.c cmc_sscollision.c\
	      cmc_stellar_evolution.c cmc_trace.c cmc_utils.c cmc_vars.h \
	      common/fitslib.h common/gensearch.h common/gensort.h \
              common/taus113-v2.h cmc_bse_utils.c cmc_bse_utils.h

