include ../common/common.mk

all: addbinaries/addbinaries mkking/mkking mkplummer/mkplummer \
	mkelson/mkelson setimf/setimf setstellar/setstellar setunits/setunits

.PHONY: clean addbinaries mkking mkplummer mkelson setimf setstellar setunits

addbinaries/addbinaries: 
	cd addbinaries && $(MAKE) addbinaries

mkking/mkking: 
	cd mkking && $(MAKE) mkking

mkplummer/mkplummer:
	cd mkplummer && $(MAKE) mkplummer

mkelson/mkelson:
	cd mkelson && $(MAKE) mkelson

setimf/setimf: 
	cd setimf && $(MAKE) setimf

setstellar/setstellar: 
	cd setstellar && $(MAKE) setstellar

setunits/setunits: 
	cd setunits && $(MAKE) setunits

install: addbinaries/addbinaries mkking/mkking mkplummer/mkplummer \
	 mkelson/mkelson setimf/setimf setstellar/setstellar setunits/setunits
	cd addbinaries && $(MAKE) install
	cd mkking && $(MAKE) install
	cd mkplummer && $(MAKE) install
	cd mkelson && $(MAKE) install
	cd setimf && $(MAKE) install
	cd setstellar && $(MAKE) install
	cd setunits && $(MAKE) install

clean:
	cd addbinaries && $(MAKE) clean
	cd mkking && $(MAKE) clean
	cd mkplummer && $(MAKE) clean
	cd mkelson && $(MAKE) install
	cd setimf && $(MAKE) clean
	cd setstellar && $(MAKE) clean
	cd setunits && $(MAKE) clean
