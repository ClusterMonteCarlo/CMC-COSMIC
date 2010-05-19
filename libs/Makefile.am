noinst_LTLIBRARIES = libsupport.la
# Since this is just a convenience library, we set the -static flag to avoid 
# double compilation.
libsupport_la_CFLAGS = -static
libsupport_la_SOURCES = fitslib.c taus113-v2.c

bin_PROGRAMS= fitslibtester
fitslibtester_SOURCES = fitslibtester.c
fitslibtester_LDADD = libsupport.la

