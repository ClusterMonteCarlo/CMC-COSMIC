include common/common.mk

#----------- CUDA suppport ----------------#
#determines if cuda is compiled and if emulation mode
use_cuda=0
emu=0

ifeq ($(use_cuda), 1)
CFLAGS   += -DUSE_CUDA
CUDAOBJS = $(CUDADIR)/cmc_cuda.cu_o
CUDAINC  = -I./$(CUDADIR)/common/inc
CFLAGS	 += -L/usr/local/lib -L./$(CUDADIR)/common/lib -L./$(CUDADIR)/lib
CUDALIB  = -fPIC -lcuda -lcudart
ifeq ($(emu), 1)
CFLAGS 	+= -D__DEVICE_EMULATION__
endif 
endif
#------------------------------------------#

# standard executable
EXE = cmc
OBJS = cmc_binbin.o cmc_binsingle.o cmc_dynamics.o \
	cmc_dynamics_helper.o cmc_evolution_thr.o cmc_funcs.o \
	cmc_init.o cmc_io.o cmc.o cmc_nr.o \
	cmc_utils.o cmc_fits.o cmc_stellar_evolution.o \
	cmc_sort.o cmc_sscollision.o cmc_remove_star.o cmc_bhlosscone.o \
	cmc_search_grid.o cmc_trace.o cmc_orbit.o cmc_core.o \
	libs/fitslib.o libs/taus113-v2.o cmc_bse_utils.o
FEWBODYOBJS = $(FEWBODYDIR)/fewbody.o $(FEWBODYDIR)/fewbody_classify.o \
	$(FEWBODYDIR)/fewbody_coll.o $(FEWBODYDIR)/fewbody_hier.o \
	$(FEWBODYDIR)/fewbody_int.o $(FEWBODYDIR)/fewbody_io.o \
	$(FEWBODYDIR)/fewbody_isolate.o $(FEWBODYDIR)/fewbody_ks.o \
        $(FEWBODYDIR)/fewbody_nonks.o $(FEWBODYDIR)/fewbody_scat.o \
	$(FEWBODYDIR)/fewbody_utils.o
BSEOBJS = $(BSEWRAPDIR)/bse_wrap.o $(BSEDIR)/comenv.o $(BSEDIR)/corerd.o $(BSEDIR)/deltat.o \
        $(BSEDIR)/dgcore.o $(BSEDIR)/evolv1.o $(BSEDIR)/evolv2.o \
        $(BSEDIR)/gntage.o $(BSEDIR)/hrdiag.o $(BSEDIR)/instar.o \
        $(BSEDIR)/kick.o $(BSEDIR)/mix.o $(BSEDIR)/mlwind.o \
        $(BSEDIR)/mrenv.o $(BSEDIR)/ran3.o $(BSEDIR)/rl.o $(BSEDIR)/star.o \
        $(BSEDIR)/zcnsts.o $(BSEDIR)/zfuncs.o

# default super-target
all: $(EXE) UTILS CONTRIBS $(CUDAOBJS)

# peripheral stuff
$(BSEWRAPDIR)/bse_wrap.o: $(BSEWRAPDIR)/bse_wrap.c $(BSEWRAPDIR)/bse_wrap.h
	cd $(BSEWRAPDIR) && make && rm -f c_sse c_bse popsynth

libs/fitslib.o: 
	cd libs && $(MAKE)

libs/taus113-v2.o:
	cd libs && $(MAKE)

UTILS: libs/fitslib.o libs/taus113-v2.o
	cd utils && $(MAKE)

CONTRIBS:
	cd contrib && $(MAKE)

$(CUDAOBJS):
	cd $(CUDADIR) && $(MAKE)

# the standard executable
$(EXE): $(CUDAOBJS) $(OBJS) $(FEWBODYOBJS) $(BSEOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS) $(CUDALIB)

$(FEWBODYDIR)/%.o: $(FEWBODYDIR)/%.c $(FEWBODYDIR)/fewbody.h Makefile
	$(CC) $(CFLAGS) -I$(FEWBODYDIR) -c $< -o $@

%.o: %.c cmc.h cmc_vars.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

# fake dependency to tell make that these targets do not produce files
.PHONY: install clean mrproper

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/
	cd utils && $(MAKE) install
	cd contrib && $(MAKE) install

clean:
	rm -f $(OBJS) $(FEWBODYOBJS) $(BSEOBJS) $(EXE)
	cd utils && $(MAKE) clean
	cd libs && $(MAKE) clean
	cd bse_wrap && $(MAKE) clean
	cd $(CUDADIR) && $(MAKE) clean

mrproper: clean
	find . -name \*~ | xargs rm -f
	find . -name .smhist | xargs rm -f
	find . -name \*.dat | xargs rm -f
	find . -name \*.dat.gz | xargs rm -f
	find . -name \*.log | xargs rm -f
	find . -name \*.stdout | xargs rm -f
	find . -name \*.stderr | xargs rm -f
	find . -name \*.cmc.parsed | xargs rm -f
	find . -name \*.conv.sh | xargs rm -f
	find . -name fort.99 | xargs rm -f
	find . -name \*.fits | xargs rm -f
	find . -name \*.fit | xargs rm -f

