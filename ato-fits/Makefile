##############################################################################
### A simple Makefile to execute make in subdirectories
##############################################################################

all: kingfit massspec plumcik

.PHONY : clean kingfit massspec plumcik
kingfit: 
	@echo "*** making kingfit ***"
	cd kingfit && $(MAKE) kingfit

massspec: 
	@echo "*** making massspec ***"
	cd massspec && $(MAKE) massspec

plumcik:
	@echo "*** making plumcik ***"
	cd plumcik && $(MAKE) plumcik

clean:
	cd kingfit && $(MAKE) clean
	cd massspec && $(MAKE) clean
	cd plumcik && $(MAKE) clean
