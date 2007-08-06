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

# default super-target
all: $(EXE) UTILS CONTRIBS

# peripheral stuff
libs/fitslib.o: 
	cd libs && $(MAKE)

libs/taus113-v2.o:
	cd libs && $(MAKE)

UTILS: libs/fitslib.o libs/taus113-v2.o startrack/singl.o
	cd utils && $(MAKE)

CONTRIBS:
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

# fake dependency to tell make that these targets do not produce files
.PHONY: install clean mrproper

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/
	cd utils && $(MAKE) install
	cd contrib && $(MAKE) install

clean:
	rm -f $(OBJS) $(FEWBODYOBJS) $(EXE) $(COBJS) $(FEWBODYCOBJS)
	cd utils && $(MAKE) clean
	cd libs && $(MAKE) clean

mrproper: clean
	rm -f   *~   .smhist   *.dat   *.dat.gz    *out_*   *.stdout   *.stderr   *.log
	rm -f */*~ */.smhist */*.dat */*.dat.gz  */*out_* */*.stdout */*.stderr */*.log
