include ../../common/common.mk

# unset FLIBS as here it is not needed
FLIBS=

EXE = setimf
OBJS = $(EXE).o genimf.o ../../libs/taus113-v2.o ../../libs/fitslib.o

#pattern rule to compile object files from C files
%.o : %.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) $(LIBFLAGS)

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/cmc_$^

clean:
	rm -f $(EXE).o $(EXE) genimf.o
