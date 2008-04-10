PREFIX = $(HOME)

CC = gcc

##############################################################################
### test for ccache
##############################################################################
CCACHE = $(shell which ccache 2>/dev/null)

ifneq ($(CCACHE),)
CC := ccache $(CC)
endif

COPTS = -Wall -O3 -I../gensort -I../taus113-v2
LIBS = -lcfitsio -lgsl -lgslcblas -lm

EXE = nbody2fits
OBJS = $(EXE).o ../taus113-v2/taus113-v2.o

#pattern rule to compile object files from C files
%.o : %.c Makefile
	$(CC) $(COPTS) -c $< -o $@

$(EXE): $(OBJS) Makefile
	$(CC) $(COPTS) $(OBJS) -o $(EXE) $(LIBS)

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/

clean:
	rm -f $(OBJS) $(EXE)
