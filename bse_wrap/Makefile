CC = gcc
CFLAGS = -Wall -O3 -L/sw/lib
LIBFLAGS = -lm -lg2c

BSEDIR = bse
BSEOBJS = $(BSEDIR)/comenv.o $(BSEDIR)/corerd.o $(BSEDIR)/deltat.o \
	$(BSEDIR)/dgcore.o $(BSEDIR)/evolv1.o $(BSEDIR)/evolv2.o \
	$(BSEDIR)/gntage.o $(BSEDIR)/hrdiag.o $(BSEDIR)/instar.o \
	$(BSEDIR)/kick.o $(BSEDIR)/mix.o $(BSEDIR)/mlwind.o \
	$(BSEDIR)/mrenv.o $(BSEDIR)/ran3.o $(BSEDIR)/rl.o $(BSEDIR)/star.o \
	$(BSEDIR)/zcnsts.o $(BSEDIR)/zfuncs.o

all: c_sse c_bse

c_sse: c_sse.o bse_wrap.o $(BSEOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

c_sse.o: c_sse.c bse_wrap.h
	$(CC) $(CFLAGS) -c $< -o $@

c_bse: c_bse.o bse_wrap.o $(BSEOBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

c_bse.o: c_bse.c bse_wrap.h
	$(CC) $(CFLAGS) -c $< -o $@

bse_wrap.o: bse_wrap.c bse_wrap.h
	$(CC) $(CFLAGS) -c $< -o $@

$(BSEDIR)/%.o: $(BSEDIR)/%.f $(BSEDIR)/const_bse.h $(BSEDIR)/zdata.h
	cd $(BSEDIR) && $(MAKE) `echo $@ | cut -d / -f 2`

bseclean:
	rm -f $(BSEOBJS)

localclean:
	rm -f bse_wrap.o c_sse.o c_sse c_bse.o c_bse

clean: bseclean localclean
