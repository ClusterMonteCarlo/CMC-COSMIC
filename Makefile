include common/common.mk

# standard executable
EXE = cmc
OBJS = cmc_binbin.o cmc_binsingle.o cmc_dynamics.o \
	cmc_dynamics_helper.o cmc_evolution_thr.o cmc_funcs.o \
	cmc_init.o cmc_io.o cmc.o cmc_nr.o \
	cmc_utils.o cmc_fits.o \
	cmc_stellar_evolution.o cmc_fits_sshot.o \
	cmc_sort.o cmc_sscollision.o cmc_bhlosscone.o \
	cmc_search_grid.o cmc_trace.o \
	startrack/singl.o \
	libs/fitslib.o libs/taus113-v2.o
FEWBODYOBJS = $(FEWBODYDIR)/fewbody.o $(FEWBODYDIR)/fewbody_classify.o \
	$(FEWBODYDIR)/fewbody_coll.o $(FEWBODYDIR)/fewbody_hier.o \
	$(FEWBODYDIR)/fewbody_int.o $(FEWBODYDIR)/fewbody_io.o \
	$(FEWBODYDIR)/fewbody_isolate.o $(FEWBODYDIR)/fewbody_ks.o \
        $(FEWBODYDIR)/fewbody_nonks.o $(FEWBODYDIR)/fewbody_scat.o \
	$(FEWBODYDIR)/fewbody_utils.o

# condor executable
CONDOREXE := cmc.condor
COBJS = cmc_binbin.co cmc_binsingle.co cmc_dynamics.co \
	cmc_dynamics_helper.co cmc_evolution_thr.co cmc_funcs.co \
	cmc_init.co cmc_io.co cmc.co cmc_nr.co \
	cmc_utils.co cmc_fits.co startrack/singl.co \
	cmc_stellar_evolution.co cmc_fits_sshot.co \
	cmc_sort.co cmc_sscollision.co cmc_bhlosscone.co \
	cmc_search_grid.co cmc_trace.co \
	libs/fitslib.co libs/taus113-v2.co
FEWBODYCOBJS = $(FEWBODYDIR)/fewbody.co $(FEWBODYDIR)/fewbody_classify.co \
	$(FEWBODYDIR)/fewbody_coll.co $(FEWBODYDIR)/fewbody_hier.co \
	$(FEWBODYDIR)/fewbody_int.co $(FEWBODYDIR)/fewbody_io.co \
	$(FEWBODYDIR)/fewbody_isolate.co $(FEWBODYDIR)/fewbody_ks.co \
        $(FEWBODYDIR)/fewbody_nonks.co $(FEWBODYDIR)/fewbody_scat.co \
	$(FEWBODYDIR)/fewbody_utils.co

# everything available
ifneq ($(CONDOR),)
ALLEXES = $(EXE) $(CONDOREXE)
else
ALLEXES = $(EXE)
endif

# default super-target
all: $(ALLEXES) UTILS CONTRIBS

# peripheral stuff
libs/fitslib.o: LIBS

libs/taus113-v2.o: LIBS

LIBS:
	cd libs && $(MAKE)

UTILS: LIBS
	cd utils && $(MAKE)

CONTRIBS: LIBS
	cd contrib && $(MAKE)

# the standard executable
$(EXE): $(OBJS) $(FEWBODYOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

startrack/singl.o: startrack/singl.c Makefile
	$(CC) $(CFLAGS) $(CHRISCFLAGS) -c $< -o $@

$(FEWBODYDIR)/%.o: $(FEWBODYDIR)/%.c $(FEWBODYDIR)/fewbody.h Makefile
	$(CC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

%.o: %.c cmc.h cmc_vars.h Makefile
	$(CC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

# the condor executable
$(CONDOREXE): $(COBJS) $(FEWBODYCOBJS)
	$(CONDORCC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

startrack/singl.co: startrack/singl.c Makefile
	$(CONDORCC) $(CFLAGS) $(CHRISCFLAGS) -c $< -o $@

$(FEWBODYDIR)/%.co: $(FEWBODYDIR)/%.c $(FEWBODYDIR)/fewbody.h Makefile
	$(CONDORCC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

%.co: %.c cmc.h cmc_vars.h Makefile
	$(CONDORCC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

# fake dependency to tell make that these targets do not produce files
.PHONY: install clean mrproper

install: $(ALLEXES)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/
	cd utils && $(MAKE) install
	cd contrib && $(MAKE) install

clean:
	rm -f $(OBJS) $(FEWBODYOBJS) $(EXE) $(COBJS) $(FEWBODYCOBJS) $(CONDOREXE)
	cd utils && $(MAKE) clean
	cd libs && $(MAKE) clean

mrproper: clean
	rm -f   *~   .smhist   *.dat   *.dat.gz    *out_*   *.stdout   *.stderr   *.log
	rm -f */*~ */.smhist */*.dat */*.dat.gz  */*out_* */*.stdout */*.stderr */*.log
