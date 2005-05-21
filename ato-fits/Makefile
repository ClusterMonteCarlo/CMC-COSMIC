##############################################################################
### A simple Makefile to execute make in subdirectories
##############################################################################

all: kingfit/kingfit massspec/massspec plumcik/plumcik

.PHONY: clean kingfit massspec plumcik
kingfit/kingfit: 
	cd kingfit && $(MAKE) kingfit

massspec/massspec: 
	cd massspec && $(MAKE) massspec

plumcik/plumcik:
	cd plumcik && $(MAKE) plumcik

install: kingfit/kingfit massspec/massspec plumcik/plumcik
	cd kingfit && $(MAKE) install
	cd massspec && $(MAKE) install
	cd plumcik && $(MAKE) install

clean:
	cd kingfit && $(MAKE) clean
	cd massspec && $(MAKE) clean
	cd plumcik && $(MAKE) clean
