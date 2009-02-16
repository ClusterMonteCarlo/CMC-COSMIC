include ../common/common.mk

#CC = gcc
#CFLAGS = -Wall -O3 -L/sw/lib
#LIBFLAGS = -lm -lg2c

BSEOBJS = bse/comenv.o bse/corerd.o bse/deltat.o \
	bse/dgcore.o bse/evolv1.o bse/evolv2.o \
	bse/gntage.o bse/hrdiag.o bse/instar.o \
	bse/kick.o bse/mix.o bse/mlwind.o \
	bse/mrenv.o bse/ran3.o bse/rl.o bse/star.o \
	bse/zcnsts.o bse/zfuncs.o

all: c_sse c_bse popsynth

c_sse: c_sse.o bse_wrap.o $(BSEOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

c_sse.o: c_sse.c bse_wrap.h
	$(CC) $(CFLAGS) -c $< -o $@

c_bse: c_bse.o bse_wrap.o $(BSEOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

c_bse.o: c_bse.c bse_wrap.h
	$(CC) $(CFLAGS) -c $< -o $@

popsynth: popsynth.o bse_wrap.o $(BSEOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

popsynth.o: popsynth.c popsynth.h bse_wrap.h
	$(CC) $(CFLAGS) -c $< -o $@

bse_wrap.o: bse_wrap.c bse_wrap.h
	$(CC) $(CFLAGS) -c $< -o $@

bse/%.o: bse/%.f bse/const_bse.h bse/zdata.h
	cd bse && $(MAKE) `echo $@ | cut -d / -f 2`

bseclean:
	rm -f $(BSEOBJS)

localclean:
	rm -f bse_wrap.o c_sse.o c_sse c_bse.o c_bse popsynth.o popsynth

clean: bseclean localclean
