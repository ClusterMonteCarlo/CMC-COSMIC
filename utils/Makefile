##############################################################################
### A simple Makefile to execute make in subdirectories
##############################################################################

all: mkking/mkking setimf/setimf mkplummer/mkplummer

.PHONY: clean mkking setimf mkplummer
mkking/mkking: 
	cd mkking && $(MAKE) mkking

setimf/setimf: 
	cd setimf && $(MAKE) setimf

mkplummer/mkplummer:
	cd mkplummer && $(MAKE) mkplummer

install: mkking/mkking setimf/setimf mkplummer/mkplummer
	cd mkking && $(MAKE) install
	cd setimf && $(MAKE) install
	cd mkplummer && $(MAKE) install

clean:
	cd mkking && $(MAKE) clean
	cd setimf && $(MAKE) clean
	cd mkplummer && $(MAKE) clean
