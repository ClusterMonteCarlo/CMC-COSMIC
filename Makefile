##############################################################################
### installation prefix; executables are installed in $(PREFIX)/bin
##############################################################################
PREFIX = $(HOME)

##############################################################################
### the name of the code can be changed (with some work)
##############################################################################
NAME = cmc

##############################################################################
### C compiler
##############################################################################
CC = gcc

##############################################################################
### test for condor
##############################################################################
CONDOR = $(shell type condor_compile 2>/dev/null)

ifneq ($(CONDOR),)
CONDORCC := condor_compile $(CC)
endif

##############################################################################
### test for ccache
##############################################################################
CCACHE = $(shell type ccache 2>/dev/null)

ifneq ($(CCACHE),)
CC := ccache $(CC)
endif

##############################################################################
### set versioning
##############################################################################
VERSION = $(shell grep revision .svn/entries | cut -d \" -f 2)
DATE = $(shell date | sed -e 's|[[:space:]]|_|g')

##############################################################################
### test for architecture to set CFLAGS and LIBFLAGS
##############################################################################
UNAME = $(shell uname)

ifeq ($(UNAME),Linux)
CFLAGS = -Wall -O3 -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\""
LIBFLAGS = -lpthread -lz -lgsl -lgslcblas -lcfitsio -lm
else
ifeq ($(UNAME),Darwin)
CC = gcc
#CFLAGS = -Wall -O3 -fast -I/sw/include -I/sw/include/gnugetopt -L/sw/lib -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\""
CFLAGS = -Wall -O3 -I/opt/local/include -L/opt/local/lib -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\""
LIBFLAGS = -lz -lgsl -lgslcblas -lcfitsio -lm
else
ifeq ($(UNAME),AIX)
CFLAGS = -Wall -O3 -I/u/ac/fregeau/local/include -L/u/ac/fregeau/local/lib -I/usr/local/include -L/usr/local/lib -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\""
LIBFLAGS = -lz -lgsl -lgslcblas -lcfitsio -liberty -lm
else
CFLAGS = -Wall -O3 -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\""
LIBFLAGS = -lpthread -lz -lgsl -lgslcblas -lcfitsio -lm
endif
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
ifeq ($(HOSTNAME),master.cluster)
CFLAGS := $(CFLAGS) -march=athlon-mp -I/opt/gsl/include -L/opt/gsl/lib
LIBFLAGS := $(LIBFLAGS) -static
endif

ifeq ($(HOSTNAME),fugu.phys.northwestern.edu)
#CC = pathcc
#CFLAGS := -Wall -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\"" -Ofast -OPT:fast_math=on -LNO:fu=9:full_unroll_size=7000 -static-data -I/usr/include/cfitsio
#CFLAGS := $(CFLAGS) -march=opteron -I/usr/include/cfitsio
#CC = gcc
CFLAGS := $(CFLAGS) -march=k8 -I/usr/include/cfitsio
LIBFLAGS := $(LIBFLAGS) -static
endif

DOMNAME = $(shell hostname | cut -d . -f 2-)
ifeq ($(DOMNAME),ncsa.uiuc.edu)
CC = icc
CFLAGS := -wd864,1188 -I $(HOME)/libs_et_al/include -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\""
# redefine libflags, leave out -lm to link with intel math library
# turn of diagn. 864: extern inline function ... was referenced but not defined
#           and 1188: floating-point value cannot be represented exactly
LIBFLAGS = -lpthread -lz -lgsl -lgslcblas -lcfitsio
LIBFLAGS := $(LIBFLAGS) -L $(HOME)/libs_et_al/lib -static
CHRISCFLAGS = 
endif

##############################################################################
### the dependencies
##############################################################################
FEWBODYDIR = fewbody

# standard executable
EXE = $(NAME)
OBJS = $(NAME)_binbin.o $(NAME)_binsingle.o $(NAME)_dynamics.o \
	$(NAME)_dynamics_helper.o $(NAME)_evolution_thr.o $(NAME)_funcs.o \
	$(NAME)_init.o $(NAME)_io.o $(NAME).o $(NAME)_nr.o \
	$(NAME)_utils.o taus113-v2.o $(NAME)_fits.o startrack/singl.o \
	$(NAME)_stellar_evolution.o $(NAME)_fits_sshot.o \
	$(NAME)_sort.o $(NAME)_sscollision.o $(NAME)_bhlosscone.o \
	$(NAME)_search_grid.o
FEWBODYOBJS = $(FEWBODYDIR)/fewbody.o $(FEWBODYDIR)/fewbody_classify.o \
	$(FEWBODYDIR)/fewbody_coll.o $(FEWBODYDIR)/fewbody_hier.o \
	$(FEWBODYDIR)/fewbody_int.o $(FEWBODYDIR)/fewbody_io.o \
	$(FEWBODYDIR)/fewbody_isolate.o $(FEWBODYDIR)/fewbody_ks.o \
        $(FEWBODYDIR)/fewbody_nonks.o $(FEWBODYDIR)/fewbody_scat.o \
	$(FEWBODYDIR)/fewbody_utils.o

# condor executable
CONDOREXE := $(NAME).condor
COBJS = $(NAME)_binbin.co $(NAME)_binsingle.co $(NAME)_dynamics.co \
	$(NAME)_dynamics_helper.co $(NAME)_evolution_thr.co $(NAME)_funcs.co \
	$(NAME)_init.co $(NAME)_io.co $(NAME).co $(NAME)_nr.co \
	$(NAME)_utils.co taus113-v2.co $(NAME)_fits.co startrack/singl.co \
	$(NAME)_stellar_evolution.co $(NAME)_fits_sshot.co \
	$(NAME)_sort.co $(NAME)_sscollision.co $(NAME)_bhlosscone.co \
	$(NAME)_search_grid.co
FEWBODYCOBJS = $(FEWBODYDIR)/fewbody.co $(FEWBODYDIR)/fewbody_classify.co \
	$(FEWBODYDIR)/fewbody_coll.co $(FEWBODYDIR)/fewbody_hier.co \
	$(FEWBODYDIR)/fewbody_int.co $(FEWBODYDIR)/fewbody_io.co \
	$(FEWBODYDIR)/fewbody_isolate.co $(FEWBODYDIR)/fewbody_ks.co \
        $(FEWBODYDIR)/fewbody_nonks.co $(FEWBODYDIR)/fewbody_scat.co \
	$(FEWBODYDIR)/fewbody_utils.co

# extra stuff
CONTRIBS = contrib/calc_2ddensity.pl contrib/calc_3ddensity.pl \
	contrib/calc_distfunc.pl contrib/pluckoutlog.pl contrib/pluckbindata.pl \
	contrib/prunedata.pl contrib/beo-genpbs.pl \
	contrib/quick_cluster_plot.sh contrib/quick_rel_plot.sh \
	contrib/extract_merger_tree.pl contrib/quick_binary_plot.sh \
	contrib/extract_bins.sh contrib/cluster_truncate.pl \
	contrib/quick_energy_plot.sh contrib/quick_binary_plot_noe.sh
EXTRAS = FITS

# everything available
ifneq ($(CONDOR),)
ALLEXES = $(EXE) $(CONDOREXE)
else
ALLEXES = $(EXE)
endif

all: $(ALLEXES) $(EXTRAS)

# the standard executable
$(EXE): $(OBJS) $(FEWBODYOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

startrack/singl.o: startrack/singl.c Makefile
	$(CC) $(CFLAGS) $(CHRISCFLAGS) -c $< -o $@

$(FEWBODYDIR)/%.o: $(FEWBODYDIR)/%.c $(FEWBODYDIR)/fewbody.h Makefile
	$(CC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

%.o: %.c $(NAME).h $(NAME)_vars.h Makefile
	$(CC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

# the condor executable
$(CONDOREXE): $(COBJS) $(FEWBODYCOBJS)
	$(CONDORCC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

startrack/singl.co: startrack/singl.c Makefile
	$(CONDORCC) $(CFLAGS) $(CHRISCFLAGS) -c $< -o $@

$(FEWBODYDIR)/%.co: $(FEWBODYDIR)/%.c $(FEWBODYDIR)/fewbody.h Makefile
	$(CONDORCC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

%.co: %.c $(NAME).h $(NAME)_vars.h Makefile
	$(CONDORCC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

# fake dependency
.PHONY: FITS install clean fewbodyclean mrproper

# Ato's FITS stuff
FITS:
	cd ato-fits && $(MAKE)

install: $(ALLEXES) $(CONTRIBS)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/
	cd ato-fits && $(MAKE) install

clean:
	rm -f $(OBJS) $(FEWBODYOBJS) $(EXE) $(COBJS) $(FEWBODYCOBJS) $(CONDOREXE)
	cd ato-fits && $(MAKE) clean

mrproper: clean
	rm -f   *~   .smhist   *.dat   *.dat.gz    *out_*   *.stdout   *.stderr   *.log
	rm -f */*~ */.smhist */*.dat */*.dat.gz  */*out_* */*.stdout */*.stderr */*.log
