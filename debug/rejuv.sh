cmc_mkking -w 4 -N 1000000 -o /tmp/test1.fits -s 659139000
cmc_setimf -i /tmp/test1.fits -o /tmp/test2.fits -I 0 -m 0.1 -M 5.0 && rm /tmp/test1.fits
cmc_setunits -i /tmp/test2.fits -Z 0.02 -o /tmp/test3.fits -R 1.5 && rm /tmp/test2.fits
cmc_setstellar -i /tmp/test3.fits -o /tmp/test4.fits && rm /tmp/test3.fits
cmc_addbinaries -i /tmp/test4.fits -o ../k4k_0.1_5_b10_n1e6.fits -N 100000 -l 0 && rm /tmp/test4.fits
