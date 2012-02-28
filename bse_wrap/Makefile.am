SUBDIRS= bse
noinst_LTLIBRARIES = libbsewrap.la
# Since this is just a convenience library, we set the -static flag to avoid 
# double compilation.
libbsewrap_la_CFLAGS=-static
# We cannot set AM_CFLAGS as this would interfere with the compilation of c_sse
# etc. as AM_CFLAGS are also used for linking. The meaning of -static in link
# mode is different and overrides the -all-static LDFLAG, breaking static
# linkning of these programs.
libbsewrap_la_SOURCES = bse_wrap.c bse_wrap.h
# FLIBS is for the standard libraries of the Fortran compiler.
libbsewrap_la_LIBADD = bse/libbse.la $(FLIBS)

bin_PROGRAMS = c_sse c_bse popsynth
AM_LDFLAGS = $(LIBTOOL_LDFLAGS)
c_sse_SOURCES = c_sse.c
c_sse_LDADD = libbsewrap.la
c_bse_SOURCES = c_bse.c
c_bse_LDADD = libbsewrap.la
popsynth_SOURCES = popsynth.c
popsynth_LDADD = libbsewrap.la

