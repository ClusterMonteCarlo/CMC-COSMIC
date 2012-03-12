#echo $1
for ((  i = 0 ;  i < $1;  i++  ))
do
#	echo "test_out_par$i.dat"
	cat "test_out_par$i.dat" >> test_out_par.dat
	rm "test_out_par$i.dat"
done
