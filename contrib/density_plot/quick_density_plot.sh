xmgrace -pexec "arrange(3,1,0.15,0,0.2)" -graph 0 -log y -log x -block $1.snap$2.3d -settype xy -bxy 1:2 -pexec "world xmin 0.0001" -pexec "xaxis ticklabel off" -pexec "yaxis label \"r\sh\N\S-3\"" -graph 1 -log y -log x -block $1.snap$2.2d -settype xy -bxy 1:2 -pexec "yaxis label \"r\sh\N\S-2\"" -pexec "xaxis ticklabel off" -graph 2 -log y -log x -block $1.snap$2.2d -settype xy -bxy 1:3 -pexec "yaxis label \"L(L\ssun\N)/r\sh\N\S2\"" -pexec "xaxis label \"r [r\sh\N]\"" &

echo "Done."
