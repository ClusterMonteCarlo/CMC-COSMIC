include ../common/common.mk

all: fitslib.o fitslibtester

%.o : %.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

fitslibtester: fitslibtester.o fitslib.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

clean:
	rm -f fitslib.o fitslibtester.o fitslibtester
