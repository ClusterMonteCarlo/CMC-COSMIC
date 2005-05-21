PREFIX = $(HOME)

CC = gcc
CFLAGS = -Wall -O3 -I../gensort -I../taus113-v2
LIBFLAGS = -lcfitsio -lgsl -lgslcblas -lm

##############################################################################
### test for ccache
##############################################################################
CCACHE = $(shell which ccache 2>/dev/null)

ifneq ($(CCACHE),)
CC := ccache $(CC)
endif

DOMNAME = $(shell dnsdomainname)
ifeq ($(DOMNAME),ncsa.uiuc.edu)
CC = icc
CFLAGS := $(CFLAGS) -I $(HOME)/libs_et_al/include
LIBFLAGS := $(LIBFLAGS) -L $(HOME)/libs_et_al/lib -static
endif

EXE = massspec
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
