#!/bin/sh

echo "Creating $1.etc.dat..."
echo "# t(Gyr) N_NS(core)/N_NS N_NS_in_bin(core)/N_NS N_bin(core)/N_bin" > $1.etc.dat
cat $1.out_etc | awk '{print $2 * 6.03088e+10 / 1e9 " " $7/$5 " " $10/$5 " " $6/$4}' >> $1.etc.dat

echo "Creating $1.escallbin.dat..."
echo "# t(Gyr) log_10(r(vir)) m_1(Msun) m_2(Msun) a(AU) log_10(P_orb(day)) e^2" > $1.escallbin.dat
cat $1.out_esc | awk '{if ($15==5) print $2 * 6.03088e+10 / 1e9 " " log($5)/log(10) " " $16 " " $17 " " $18 " " log(sqrt($18^3/($16+$17))*365.25)/log(10) " " $19^2}' >> $1.escallbin.dat

echo "Creating $1.escnsbin.dat..."
echo "# t(Gyr) log_10(r(vir)) m_1(Msun) m_2(Msun) a(AU) log_10(P_orb(day)) e^2" > $1.escnsbin.dat
cat $1.out_esc | awk '{if (($15==5)&&(($16==1.4)||($17==1.4))) print $2 * 6.03088e+10 / 1e9 " " log($5)/log(10) " " $16 " " $17 " " $18 " " log(sqrt($18^3/($16+$17))*365.25)/log(10) " " $19^2}' >> $1.escnsbin.dat

echo "Creating $1.escbinns.dat..."
echo "# t(Gyr) log_10(r(vir)) m_1(Msun) m_2(Msun) a(AU) log_10(P_orb(day)) e^2" > $1.escbinns.dat
cat $1.out_esc | awk '{if (($15==5)&&(($16==1.4)&&($17==1.4))) print $2 * 6.03088e+10 / 1e9 " " log($5)/log(10) " " $16 " " $17 " " $18 " " log(sqrt($18^3/($16+$17))*365.25)/log(10) " " $19^2}' >> $1.escbinns.dat

for i in $1.out_etc????.gz; do echo "Creating $i.allbin.dat..."; \
    echo "# [line from etc file] log_10(P_orb(day)) e^2" > $i.allbin.dat; \
    zcat $i | awk '{print $line " " log(sqrt($8^3/($6+$7))*365.25)/log(10) " " $9^2}' >> \
    $i.allbin.dat; done

for i in $1.out_etc????.gz; do echo "Creating $i.nsbin.dat..."; \
    echo "# [line from etc file] log_10(P_orb(day)) e^2" > $i.nsbin.dat; \
    zcat $i | awk '{if(($6==1.4)||($7==1.4)) print $line " " log(sqrt($8^3/($6+$7))*365.25)/log(10) " " $9^2}' >> \
    $i.nsbin.dat; done

for i in $1.out_etc????.gz; do echo "Creating $i.binns.dat..."; 
    echo "# [line from etc file] log_10(P_orb(day)) e^2" > $i.binns.dat; \
    zcat $i | awk '{if(($6==1.4)&&($7==1.4)) print $line " " log(sqrt($8^3/($6+$7))*365.25)/log(10) " " $9^2}' >> \
    $i.binns.dat; done

for i in $1.out_snap????.gz; do echo "Creating $i.bg.logr.dat..."; 
    echo "# log_10(r(vir))" > $i.bg.logr.dat; \
    zcat $i | awk '{if (($11==1)&&($8==0.7)) print log($3)/log(10)}' > \
    $i.bg.logr.dat; done

for i in $1.out_snap????.gz; do echo "Creating $i.ns.logr.dat..."; 
    echo "# log_10(r(vir))" > $i.ns.logr.dat; \
    zcat $i | awk '{if (($11==1)&&($8==1.4)) print log($3)/log(10)}' > \
    $i.ns.logr.dat; done

for i in $1.out_snap????.gz; do echo "Creating $i.bin.logr.dat..."; 
    echo "# log_10(r(vir))" > $i.bin.logr.dat; \
    zcat $i | awk '{if ($11==5) print log($3)/log(10)}' > \
    $i.bin.logr.dat; done
