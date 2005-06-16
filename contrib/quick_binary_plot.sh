xmgrace -pexec "arrange(3,1,0.15,0,0)" \
    -graph 0 -block $1.bin.dat -settype xy -bxy 1:3 \
    -pexec "s0.y = s0.y / s0.y[0]" \
    -graph 0 -block $1.dyn.dat -settype xy -bxy 1:5 \
    -pexec "s1.y = s1.y / s1.y[0]" \
    -pexec "yaxis label \"M\sb\N/M\sb\N(0), M/M(0)\"" \
    -pexec "xaxis ticklabel off" \
    -graph 1 -block $1.bin.dat -settype xy -bxy 1:13 \
    -pexec "s0.y = s0.y / 0.25" \
    -graph 1 -block $1.bin.dat -settype xy -bxy 1:14 \
    -pexec "s1.y = s1.y / 0.25" \
    -pexec "yaxis label \"E\sbb\N, E\sbs\N [|E\sc\N(0)|]\"" \
    -pexec "xaxis ticklabel off" \
    -graph 2 -log y -block $1.dyn.dat -settype xy -bxy 1:8 \
    -graph 2 -log y -block $1.bin.dat -settype xy -bxy 1:5 \
    -graph 2 -log y -block $1.bin.dat -settype xy -bxy 1:6 \
    -pexec "yaxis label \"r\sc\N, r\sh,b\N, r\sh,s\N [r\sNB\N]\"" \
    -pexec "xaxis label \"t [FP units]\""
