include ../common/common.mk

all: taus113-v2.o fitslib.o fitslibtester taus113tester

%.o : %.c Makefile ../common/taus113-v2.h ../common/fitslib.h
	$(CC) $(CFLAGS) -c $< -o $@

fitslibtester: fitslibtester.o fitslib.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

taus113tester: taus113tester.o taus113-v2.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

clean:
	rm -f taus113-v2.o fitslib.o fitslibtester.o fitslibtester taus113tester.o taus113tester
