include ../common/common.mk

all: taus113-v2.o fitslib.o fitslibtester

%.o : %.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

fitslibtester: fitslibtester.o fitslib.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

clean:
	rm -f taus113-v2.o fitslib.o fitslibtester.o fitslibtester
