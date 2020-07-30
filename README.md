# ws35wind 
Finding a rotating magnetic wind solution from a magnetic white dwarf <br>
ref) [Kashiyama, Fujisawa, and Shigeyama 2019](https://iopscience.iop.org/article/10.3847/1538-4357/ab4e97)

*TODO: opacity table needs to be refined for high temperature regions*

## input
- **input.dat** : input parameters 
- **tabXX.dat** : opacity table 

## code
- **shooting.c** : finding a wind solution using a shooting method <br>

## how to run the code
`gcc shooting.c` 

## output
When you run the code, you will see something like <br>
> In-shot No. 1 with rA = 7.05280819262871e+08 [cm] --> end with << flag 3 >> <br>
> In-shot No. 2 with rA = 2.71255376703998e+09 [cm] --> end with << flag 2 >> <br>
> ...

It is searching a wind solution inside the Alfven radius. 

If you get a << flag 0 >>, or the "In-shot No." exceeds the maximum value (50 by default), it proceeds to the outshot sequence.

Then, you will see something like <br>
> Out-shot No. 1 with dudxA = 4.45000000000000e-01 --> end with << flag 1 >> <br>
> Out-shot No. 2 with dudxA = 4.34275532462097e-01 --> end with << flag 2 >> <br>
> ...

It is searching a wind solution outside the Alfven radius. 

If you get a << flag 0 >>, or the "Out-shot No." exceeds the maximum value (50 by default), it finishes the calculation.

The results of the final trials of the inshot and outshot sequences are outputted in 

- **output.dat** : eigenvalues of the wind solution  
- **inshot.dat** : wind profile inside the Alfven point
- **outshot.dat** : wind profile outside the Alfven point

**When you get << flag 0 >> for both the inshot and outshot, it is a wind solution**.

## how to plot the results
`cd plt` <br>
`gnuplot plt/hydro.plt`

