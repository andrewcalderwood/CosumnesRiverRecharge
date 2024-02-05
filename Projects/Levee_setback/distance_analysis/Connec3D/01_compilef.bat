gfortran -o Connec3DLarge Connec3D.FOR -fforce-addr 
rem ifx Connec3DLarge_ifort Connec3D.FOR
PAUSE
rem optimization -o, -Qx cpu optimizatoin
rem alternate run type to stop when overflow occurs -ffpe-trap=overflow
rem -O -Wall -fcheck=all -g -fbacktrace rem helps check for errors in code
rem -g -ffpe-summary=zero,invalid,overflow,underflow 