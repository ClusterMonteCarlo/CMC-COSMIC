cmc_mkking -w 5 -N 100000 -o /tmp/test.fits
cmc_setimf -i /tmp/test.fits -o /tmp/test2.fits -I 0 -m 0.05 -M 5 && rm /tmp/test.fits
cmc_setunits -i /tmp/test2.fits -o /tmp/test3.fits -R 3.5 && rm /tmp/test2.fits
cmc_setstellar -i /tmp/test3.fits -o /tmp/test4.fits && rm /tmp/test3.fits
cmc_addbinaries -i /tmp/test4.fits -o binaries.fits -N 10000 -l 0 && rm /tmp/test4.fits
