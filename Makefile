##############################################################################
### installation prefix; executables are installed in $(PREFIX)/bin
##############################################################################
PREFIX = $(HOME)

##############################################################################
### the name of the code can be changed (with some work)
##############################################################################
NAME = cmc
PRETTYNAME = ClusterMonteCarlo

##############################################################################
### C compiler
##############################################################################
CC = gcc
# Intel compiler
#CC = icc
# Portland Group compiler
#CC = pgcc

##############################################################################
### test for ccache
##############################################################################
CCACHE = $(shell which ccache 2>/dev/null)

ifneq ($(CCACHE),)
CC := ccache $(CC)
endif

##############################################################################
### test for architecture
##############################################################################
UNAME = $(shell uname)

ifeq ($(UNAME),Linux)
CFLAGS = -Wall -O3 -DPRETTYNAME=\"$(PRETTYNAME)\"
LIBFLAGS = -lpthread -lz -lgsl -lgslcblas -lcfitsio -lm
else
ifeq ($(UNAME),Darwin)
CFLAGS = -Wall -O3 -fast -I/sw/include -I/sw/include/gnugetopt -L/sw/lib -DPRETTYNAME=\"$(PRETTYNAME)\"
LIBFLAGS = -lz -lgsl -lgslcblas -lcfitsio -lm
else
CFLAGS = -Wall -O3 -DPRETTYNAME=\"$(PRETTYNAME)\"
LIBFLAGS = -lpthread -lz -lgsl -lgslcblas -lcfitsio -lm
endif
endif

##############################################################################
### extra C flags
##############################################################################
CHRISCFLAGS = -Wno-uninitialized -Wno-unused
# possibilities: 
# -DUSE_THREADS	(using threads, results will change due different RNGs)
# -DUSE_THREADS_SORT	(using threads, for sorting only, results should
# 					remain the same)
# -DUSE_FIND	(using Hoare's FIND for parallel sorting)
# -DUSE_MODIFIND	(using Zabrodsky's MODIFIND)
# 	(the default for the above two is to use Floyd and Rivest's SELECT)
# some gcc flags
#CFLAGS := $(CFLAGS) -malign-double -fomit-frame-pointer -fmove-all-movables
# Ato's debugging
#CFLAGS = -Wall -pedantic -ggdb
# Intel compiler
#CFLAGS = -O2 -mp -unroll -pc80 -tpp7 -xiMKW
# Portland Group compiler
#CFLAGS = -B -O4 -Minform,warn

##############################################################################
### special hosts
##############################################################################
ifeq ($(HOSTNAME),chinook.astro.northwestern.edu)
CC = gcc3.3.2
#CFLAGS := $(CFLAGS) -mcpu=pentium4 -mmmx -msse -msse2
else
ifeq ($(HOSTNAME),master.cluster)
CFLAGS := $(CFLAGS) -mcpu=athlon-mp -mmmx -msse -m3dnow
LIBFLAGS := $(LIBFLAGS) -static
else
endif
endif

##############################################################################
### the dependencies
##############################################################################
EXE = $(NAME)
OBJS = $(NAME)_binbin.o $(NAME)_binsingle.o \
	$(NAME)_dynamics.o $(NAME)_evolution_thr.o $(NAME)_funcs.o \
	$(NAME)_init.o $(NAME)_io.o $(NAME).o $(NAME)_nr.o \
	$(NAME)_utils.o taus113-v2.o $(NAME)_fits.o singl.o belgy.o \
	$(NAME)_fits_sshot.o $(NAME)_sort.o
FEWBODYDIR = fewbody
FEWBODYOBJS = $(FEWBODYDIR)/fewbody.o $(FEWBODYDIR)/fewbody_classify.o \
	$(FEWBODYDIR)/fewbody_coll.o $(FEWBODYDIR)/fewbody_hier.o \
	$(FEWBODYDIR)/fewbody_int.o $(FEWBODYDIR)/fewbody_io.o \
	$(FEWBODYDIR)/fewbody_isolate.o $(FEWBODYDIR)/fewbody_ks.o \
        $(FEWBODYDIR)/fewbody_nonks.o $(FEWBODYDIR)/fewbody_scat.o \
	$(FEWBODYDIR)/fewbody_utils.o
CONTRIBS = contrib/calc_2ddensity.pl contrib/calc_3ddensity.pl \
	contrib/calc_distfunc.pl contrib/pluck0.pl contrib/pluckbindata.pl \
	contrib/prunedata.pl contrib/beo-genpbs.pl
EXTRAS = FITS

all: $(EXE) $(EXTRAS)

$(EXE): $(OBJS) $(FEWBODYOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

singl.o: singl.c Makefile
	$(CC) $(CFLAGS) $(CHRISCFLAGS) -c $< -o $@

$(FEWBODYDIR)/%.o: $(FEWBODYDIR)/%.c $(FEWBODYDIR)/fewbody.h Makefile
	$(CC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

%.o: %.c $(NAME).h $(NAME)_vars.h Makefile
	$(CC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

.PHONY: FITS install clean fewbodyclean mrproper

FITS:
	cd ato-fits && $(MAKE) -s

install: $(EXE) $(CONTRIBS)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/

clean:
	rm -f $(OBJS) $(FEWBODYOBJS) $(EXE)
	cd ato-fits && $(MAKE) -s clean

$(NAME)clean:
	rm -f $(OBJS) $(EXE)

fewbodyclean:
	rm -f $(FEWBODYOBJS) $(EXE)

mrproper: clean
	rm -f   *~   .smhist   *.dat   *out_*   *.stdout   *.stderr
	rm -f */*~ */.smhist */*.dat */*out_* */*.stdout */*.stderr
