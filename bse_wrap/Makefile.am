SUBDIRS= bse
noinst_LTLIBRARIES = libbsewrap.la
libbsewrap_la_SOURCES = bse_wrap.c bse_wrap.h
libbsewrap_la_LIBADD = bse/libbse.la $(FLIBS)

bin_PROGRAMS = c_sse c_bse popsynth
AM_LDFLAGS = $(LIBTOOL_LDFLAGS)
c_sse_SOURCES = c_sse.c
c_sse_LDADD = libbsewrap.la
c_bse_SOURCES = c_bse.c
c_bse_LDADD = libbsewrap.la
popsynth_SOURCES = popsynth.c
popsynth_LDADD = libbsewrap.la

