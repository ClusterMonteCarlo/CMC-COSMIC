PREFIX = $(HOME)

CC = gcc
COPTS = -Wall -pedantic -O3 -I../gensort -I../taus113-v2
LIBS = -lcfitsio -lgsl -lgslcblas -lm

EXE = kingfit
OBJS = $(EXE).o ../taus113-v2/taus113-v2.o

#pattern rule to compile object files from C files
%.o : %.c Makefile
	$(CC) $(COPTS) -c $< -o $@

$(EXE): $(OBJS) Makefile
	$(CC) $(COPTS) $(OBJS) -o $(EXE) $(LIBS)

install: $(EXE)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/

clean:
	rm -f $(OBJS) $(EXE)
