# installation prefix; executables are installed in $(PREFIX)/bin
PREFIX = $(HOME)

# test for ccache
CCACHE = $(shell which ccache 2>/dev/null)

ifneq ($(CCACHE),)
CC = ccache gcc
else
CC = gcc
endif

# test for architecture
UNAME = $(shell uname)

ifeq ($(UNAME),Linux)
CFLAGS = -Wall -O3
LIBFLAGS = -lgsl -lgslcblas -lm
else
ifeq ($(UNAME),Darwin)

# Test if gsl and getopt are installed via darwin ports or fink
FINK_DIR = /sw
DARWIN_DIR = /opt/local

ifeq ($(wildcard $(FINK_DIR)/include/gsl_sys.h), )
ifeq ($(wildcard $(DARWIN_DIR)/include/gsl/gsl_sys.h), )
$(error "Cannot find GSL headers!")
else
GSL_INCLUDE = $(DARWIN_DIR)/include
GSL_LIB = $(DARWIN_DIR)/lib
endif
else
GSL_INCLUDE = $(FINK_DIR)/include
GSL_LIB = $(FINK_DIR)/lib
endif
ifeq ($(wildcard $(FINK_DIR)/include/gnugetopt/getopt.h), )
ifeq ($(wildcard $(DARWIN_DIR)/include/gnugetopt/getopt.h), )
$(error "Cannot find GNU getopt headers!")
else
GETOPT_INCLUDE = $(DARWIN_DIR)/include/gnugetopt
GETOPT_LIB = $(DARWIN_DIR)/lib
endif
else
GETOPT_INCLUDE = $(FINK_DIR)/include/gnugetopt
GETOPT_LIB = $(FINK_DIR)/lib
endif

CFLAGS = -Wall -O3 -I$(GETOPT_INCLUDE) -I$(GSL_INCLUDE) -L$(GETOPT_LIB) -L$(GSL_LIB)
LIBFLAGS = -lgsl -lgslcblas -lgnugetopt -lm

else
CFLAGS = -Wall -O3
LIBFLAGS = -lgsl -lgslcblas -lm
endif
endif

# the core fewbody objects
FEWBODY_OBJS = fewbody.o fewbody_classify.o fewbody_coll.o fewbody_hier.o \
	fewbody_int.o fewbody_io.o fewbody_isolate.o fewbody_ks.o \
	fewbody_nonks.o fewbody_scat.o fewbody_utils.o

all: cluster triplebin binbin binsingle sigma_binsingle triple

cluster: cluster.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

triplebin: triplebin.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

binbin: binbin.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

binsingle: binsingle.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

sigma_binsingle: sigma_binsingle.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

triple: triple.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

cluster.o: cluster.c cluster.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

triplebin.o: triplebin.c triplebin.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

binbin.o: binbin.c binbin.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

binsingle.o: binsingle.c binsingle.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

sigma_binsingle.o: sigma_binsingle.c sigma_binsingle.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

triple.o: triple.c triple.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.c fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

install: cluster triplebin binbin binsingle sigma_binsingle triple
	mkdir -p $(PREFIX)/bin/
	install -m 0755 cluster triplebin binbin binsingle sigma_binsingle \
	triple $(PREFIX)/bin/

uninstall:
	rm -f $(PREFIX)/bin/cluster $(PREFIX)/bin/triplebin $(PREFIX)/bin/binbin \
	$(PREFIX)/bin/binsingle $(PREFIX)/bin/sigma_binsingle $(PREFIX)/bin/triple

clean:
	rm -f $(FEWBODY_OBJS) cluster.o triplebin.o binbin.o binsingle.o \
	sigma_binsingle.o cluster triplebin binbin binsingle sigma_binsingle \
	triple

mrproper: clean
	rm -f *~ *.bak *.dat ChangeLog
	rm -f */*~ */*.bak */*.dat
