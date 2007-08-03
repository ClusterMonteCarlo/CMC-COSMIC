include ../../common/common.mk

EXE = setimf
OBJS = $(EXE).o ../taus113-v2/taus113-v2.o

#pattern rule to compile object files from C files
%.o : %.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) $(LIBFLAGS)

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/cmc_$^

clean:
	rm -f $(OBJS) $(EXE)
