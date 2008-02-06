include ../../common/common.mk

EXE = addbinaries
BSEOBJS = ../../$(BSEWRAPDIR)/bse_wrap.o ../../$(BSEDIR)/comenv.o ../../$(BSEDIR)/corerd.o ../../$(BSEDIR)/deltat.o \
        ../../$(BSEDIR)/dgcore.o ../../$(BSEDIR)/evolv1.o ../../$(BSEDIR)/evolv2.o \
        ../../$(BSEDIR)/gntage.o ../../$(BSEDIR)/hrdiag.o ../../$(BSEDIR)/instar.o \
        ../../$(BSEDIR)/kick.o ../../$(BSEDIR)/mix.o ../../$(BSEDIR)/mlwind.o \
        ../../$(BSEDIR)/mrenv.o ../../$(BSEDIR)/ran3.o ../../$(BSEDIR)/rl.o ../../$(BSEDIR)/star.o \
        ../../$(BSEDIR)/zcnsts.o ../../$(BSEDIR)/zfuncs.o
OBJS = $(EXE).o ../../libs/taus113-v2.o ../../libs/fitslib.o $(BSEOBJS)

#pattern rule to compile object files from C files
%.o : %.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) $(LIBFLAGS)

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/cmc_$^

clean:
	rm -f $(EXE).o $(EXE)
