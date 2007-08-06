include ../common/common.mk

all: mkking/mkking mkplummer/mkplummer setimf/setimf setstellar/setstellar setunits/setunits

.PHONY: clean mkking mkplummer setimf setstellar setunits

mkking/mkking: 
	cd mkking && $(MAKE) mkking

mkplummer/mkplummer:
	cd mkplummer && $(MAKE) mkplummer

setimf/setimf: 
	cd setimf && $(MAKE) setimf

setstellar/setstellar: 
	cd setstellar && $(MAKE) setstellar

setunits/setunits: 
	cd setunits && $(MAKE) setunits

install: mkking/mkking mkplummer/mkplummer setimf/setimf setstellar/setstellar setunits/setunits
	cd mkking && $(MAKE) install
	cd mkplummer && $(MAKE) install
	cd setimf && $(MAKE) install
	cd setstellar && $(MAKE) install
	cd setunits && $(MAKE) install

clean:
	cd mkking && $(MAKE) clean
	cd mkplummer && $(MAKE) clean
	cd setimf && $(MAKE) clean
	cd setstellar && $(MAKE) clean
	cd setunits && $(MAKE) clean
