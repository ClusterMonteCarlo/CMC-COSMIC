cmc_mkking -w 4 -N 100000 -o /tmp/test1.fits -s 100
cmc_setimf -i /tmp/test1.fits -o /tmp/test2.fits -I 0 -m 0.1 -M 4.0 && rm /tmp/test1.fits
cmc_setunits -i /tmp/test2.fits -o /tmp/test3.fits -R 1.5 && rm /tmp/test2.fits
cmc_setstellar -i /tmp/test3.fits -o /tmp/test4.fits && rm /tmp/test3.fits
cmc_addbinaries -i /tmp/test4.fits -o ../k4k_0.1_4_b10_n1e5.fits -N 10000 -l 0 && rm /tmp/test4.fits
