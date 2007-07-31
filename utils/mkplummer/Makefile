PREFIX = $(HOME)

CC = gcc

##############################################################################
### test for architecture to set CFLAGS and LIBFLAGS
##############################################################################
UNAME = $(shell uname)

ifeq ($(UNAME),Linux)
CFLAGS = -Wall -O3 -I../gensort -I../taus113-v2
LIBFLAGS = -lcfitsio -lgsl -lgslcblas -lm
else
ifeq ($(UNAME),Darwin)
CC = gcc
CFLAGS = -Wall -O3 -I/opt/local/include -L/opt/local/lib -I../gensort -I../taus113-v2
LIBFLAGS = -lcfitsio -lgsl -lgslcblas -lm
else
CFLAGS = -Wall -O3 -I../gensort -I../taus113-v2
LIBFLAGS = -lcfitsio -lgsl -lgslcblas -lm
endif
endif

##############################################################################
### test for ccache
##############################################################################
CCACHE = $(shell type ccache 2>/dev/null)

ifneq ($(CCACHE),)
CC := ccache $(CC)
endif

##############################################################################
### special hosts
##############################################################################
ifeq ($(HOSTNAME),master.cluster)
CFLAGS := $(CFLAGS) -mcpu=athlon-mp -mmmx -msse -m3dnow -I/opt/gsl/include -L/opt/gsl/lib
LIBFLAGS := $(LIBFLAGS) -static
endif

ifeq ($(HOSTNAME),fugu.phys.northwestern.edu)
#CC = pathcc
#CFLAGS := -Wall -DCMCVERSION="\"$(VERSION)\"" -DCMCDATE="\"$(DATE)\"" -Ofast -OPT:fast_math=on -LNO:fu=9:full_unroll_size=7000 -static-data -I/usr/include/cfitsio
#CFLAGS := $(CFLAGS) -march=opteron -I/usr/include/cfitsio
#CC = gcc
#CFLAGS := $(CFLAGS) -march=k8 -I/usr/include/cfitsio
CFLAGS := $(CFLAGS) -m32 -march=k8 -I/share/apps/gsl/include -L/share/apps/gsl/lib -I/share/apps/cfitsio/include -L/share/apps/cfitsio/lib
LIBFLAGS := $(LIBFLAGS) -static
endif

DOMNAME = $(shell hostname | cut -d . -f 2-)
ifeq ($(DOMNAME),ncsa.uiuc.edu)
CC = icc
CFLAGS := $(CFLAGS) -I $(HOME)/libs_et_al/include
LIBFLAGS := $(LIBFLAGS) -L $(HOME)/libs_et_al/lib -static
endif

EXE = plumcik
OBJS = $(EXE).o ../taus113-v2/taus113-v2.o

#pattern rule to compile object files from C files
%.o : %.c mkplummer
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) $(LIBFLAGS)

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/cmc_$^

clean:
	rm -f $(OBJS) $(EXE)
