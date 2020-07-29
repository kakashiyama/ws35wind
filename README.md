# ws35wind 
Finding a rotating magnetic wind solution from a magnetic white dwarf

## input
- input.dat : input parameters <br>
- tabXX.dat : opacity table 

## code
shooting.c : finding a wind solution using a shooting method <br>
ref) Kashiyama, Fujisawa, and Shigeyama 2019

## how to run the code
`gcc shooting.c` 

## output
- output.dat : eigenvalues of the wind solution  
- inshot.dat : wind profile inside the Alfven point
- outshot.dat : wind profile outside the Alfven point

## how to plot the results
`cd plt` <br>
`gnuplot plt/hydro.plt`

