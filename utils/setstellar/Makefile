include ../../common/common.mk

EXE = setstellar
OBJS = $(EXE).o ../../libs/taus113-v2.o ../../libs/fitslib.o ../../startrack/singl.o

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
