include ../common/common.mk

all: mkking/mkking mkplummer/mkplummer setimf/setimf setunits/setunits

.PHONY: clean mkking mkplummer setunits

mkking/mkking: 
	cd mkking && $(MAKE) mkking

mkplummer/mkplummer:
	cd mkplummer && $(MAKE) mkplummer

setimf/setimf: 
	cd setimf && $(MAKE) setimf

setunits/setunits: 
	cd setunits && $(MAKE) setunits

install: mkking/mkking mkplummer/mkplummer setimf/setimf setunits/setunits
	cd mkking && $(MAKE) install
	cd mkplummer && $(MAKE) install
	cd setimf && $(MAKE) install
	cd setunits && $(MAKE) install

clean:
	cd mkking && $(MAKE) clean
	cd mkplummer && $(MAKE) clean
	cd setimf && $(MAKE) clean
	cd setunits && $(MAKE) clean
