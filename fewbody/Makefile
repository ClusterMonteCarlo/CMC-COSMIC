# installation prefix; executables are installed in $(PREFIX)/bin
PREFIX = $(HOME)

# normal compiling
CC = gcc
CFLAGS = -Wall -O3
LIBFLAGS = -lgsl -lgslcblas -lm

# profiling
#CC = gcc
#CFLAGS = -Wall -Winline -pedantic -O3 -pg
#LIBFLAGS = -lgsl -lgslcblas -lm

# hardcore debugging, a la GSL
#CC = gcc
#CFLAGS = -ansi -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes \
#	 -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual \
#	 -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums \
#	 -Dinline= -g -O4
#LIBFLAGS = -lgsl -lgslcblas -lm -lefence

FEWBODY_OBJS = fewbody.o fewbody_classify.o fewbody_coll.o fewbody_hier.o \
	fewbody_int.o fewbody_io.o fewbody_isolate.o fewbody_ks.o \
	fewbody_nonks.o fewbody_scat.o fewbody_utils.o

all: cluster triplebin binbin binsingle

cluster: cluster.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

triplebin: triplebin.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

binbin: binbin.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

binsingle: binsingle.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

cluster.o: cluster.c cluster.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

triplebin.o: triplebin.c triplebin.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

binbin.o: binbin.c binbin.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

binsingle.o: binsingle.c binsingle.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.c fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

install: cluster triplebin binbin binsingle
	mkdir -p $(PREFIX)/bin/
	install -m 0755 cluster triplebin binbin binsingle $(PREFIX)/bin/

uninstall:
	rm -f $(PREFIX)/bin/cluster $(PREFIX)/bin/triplebin $(PREFIX)/bin/binbin \
	$(PREFIX)/bin/binsingle

clean:
	rm -f $(FEWBODY_OBJS) cluster.o triplebin.o binbin.o binsingle.o \
	cluster triplebin binbin binsingle

mrproper: clean
	rm -f *~ *.bak *.dat
	rm -f */*~ */*.bak */*.dat
