Rsun = 6.95510e10

se ter po eps enhanced color
se tics font "Times-Roman"
se size 0.7,0.7
se output "vr.eps"
se log
se xrange [1.e8:]
se yrange [:1.e10]
se format x "10^{%L}"
se format y "10^{%L}"
se xlabel "radius  r [cm]" font "Times-Roman, 16"
se ylabel "wind radial velocity  v_r [cm/s]" font "Times-Roman, 16"
plot \
"../inshot.dat" usi 6:7 w l title "",\
"../outshot.dat" usi 6:7 w l title ""

unse yrange
se output "T.eps"
se ylabel "temperature  T [K]" font "Times-Roman, 16"
plot \
"../inshot.dat" usi 6:8 w l title "",\
"../outshot.dat" usi 6:8 w l title ""

unse yrange
se output "rho.eps"
se ylabel "density  {/Symbol r} [g/cm^3]" font "Times-Roman, 16"
plot \
"../inshot.dat" usi 6:9 w l title "",\
"../outshot.dat" usi 6:9 w l title ""

unse yrange
se output "vphi.eps"
se ylabel "wind tangential velocity  v_{/Symbol j} [cm/s]" font "Times-Roman, 16"
plot \
"../inshot.dat" usi 6:($10) w l title "",\
"../outshot.dat" usi 6:($10) w l title ""

unse yrange
se output "Bphi.eps"
se ylabel "tangential magnetic field  B_{/Symbol j} [G]" font "Times-Roman, 16"
unse yrange 
plot \
"../inshot.dat" usi 6:(-$12) w l title "",\
"../outshot.dat" usi 6:(-$12) w l title ""

unse yrange
unse format y
se output "Lr.eps"
se ylabel "luminosity  L_{rad} [erg/s]" font "Times-Roman, 16"
plot \
"../inshot.dat" usi 6:13 w l title "",\
"../outshot.dat" usi 6:13 w l title ""

unse yrange
se output "kappa.eps"
se ylabel "opacity  {/Symbol k} [g/cm^2]" font "Times-Roman, 16"
plot \
"../inshot.dat" usi 6:14 w l title "",\
"../outshot.dat" usi 6:14 w l title ""