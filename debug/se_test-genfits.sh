cmc_mkplummer -N 100000 -o /tmp/test.fits
cmc_setimf -i /tmp/test.fits -o /tmp/test2.fits -I 0 -m 0.7 -M 50 -s 1 && rm /tmp/test.fits
cmc_setunits -i /tmp/test2.fits -o /tmp/test3.fits -R 1.5 && rm /tmp/test2.fits
cmc_setstellar -i /tmp/test3.fits -o se_test.fits && rm /tmp/test3.fits
