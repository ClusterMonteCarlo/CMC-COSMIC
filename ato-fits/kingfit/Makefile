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

DOMNAME = $(shell hostname | cut -d . -f 2-)
ifeq ($(DOMNAME),ncsa.uiuc.edu)
CC = icc
CFLAGS := $(CFLAGS) -I $(HOME)/libs_et_al/include
LIBFLAGS := $(LIBFLAGS) -L $(HOME)/libs_et_al/lib -static
endif


EXE = kingfit
OBJS = $(EXE).o ../taus113-v2/taus113-v2.o

#pattern rule to compile object files from C files
%.o : %.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) $(LIBFLAGS)

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/

clean:
	rm -f $(OBJS) $(EXE)
