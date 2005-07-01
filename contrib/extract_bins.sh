zcat $1 | awk '{if ($8 == 1) { print $9 " " $10 " " $11 " " $12 " " $13 " " $14}}' > `echo $1 | sed -e 's|dat.gz|bin.dat|'`
