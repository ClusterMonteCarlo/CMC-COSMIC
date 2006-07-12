# Usage check_snaptime.sh <file_prefix> <snap #>
# e.g. check_snaptime.sh w6_n1e5_fb0.03.out 0133
export time=`zcat $1.snap$2.dat.gz | head -n 1 | awk '{print $2}' | cut -d = -f 2`
echo snapshot time = $time
