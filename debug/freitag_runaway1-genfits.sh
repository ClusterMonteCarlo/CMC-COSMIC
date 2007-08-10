cmc_mkking -N 1250000 -w 3 -o /tmp/test.fits
cmc_setimf -i /tmp/test.fits -o /tmp/test2.fits -I 1 -p -2.35 -m 0.2 -M 120 && rm /tmp/test.fits
cmc_setunits -i /tmp/test2.fits -o /tmp/test3.fits -R 0.1 -T 1000000 && rm /tmp/test2.fits
cmc_setstellar -i /tmp/test3.fits -o freitag_runaway1.fits && rm /tmp/test3.fits
