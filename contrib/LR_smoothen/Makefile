SHELL = /bin/sh
# decide on compiler and flags depending on host
ifeq ($(HOSTNAME),chinook.astro.northwestern.edu)
CC = /usr/local/bin/gcc3.3.2
F77 = /usr/local/bin/g773.3.2
CFLAGS = -mcpu=pentium4 -mmmx -msse -msse2
else
ifeq ($(HOSTNAME),master.cluster)
CC = /usr/local/bin/gcc3
F77 = /usr/local/bin/g773
CFLAGS = -mcpu=athlon-mp -mmmx -msse -m3dnow
else
CC = gcc
F77 = g77
CFLAGS =
endif
endif
ALL_CFLAGS = -Wall -O6 $(CFLAGS)

OBJS = LRproc.o cubgcv.o
LIBS = -lm
EXE = LRproc

# pattern rule to compile object files from C files
# might not work with make programs other than GNU make
%.o : %.c Makefile
	$(CC) $(ALL_CFLAGS) -c $< -o $@
%.o : %.f Makefile
	$(F77) $(ALL_CFLAGS) -c $< -o $@

all: $(EXE)

$(EXE): $(OBJS) Makefile
	$(F77) $(ALL_CFLAGS) $(OBJS) -o $(EXE) $(LIBS)

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXE)
