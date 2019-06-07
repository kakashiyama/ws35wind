Rsun = 6.95510e10

se ter po eps enhanced color
se tics font "Times-Roman"
se size 0.7,0.7
se output "Lr.eps"
se log
se xrange [1.e8:10*Rsun]
se format x "10^{%L}"
se format y "10^{%L}"
se xlabel "radius  r [cm]" font "Times-Roman, 16"
se ylabel "luminosity  L_{rad} [erg/s]" font "Times-Roman, 16"
plot \
"../inshot.dat" usi 6:13 w l title "",\
"../outshot.dat" usi 6:13 w l title ""