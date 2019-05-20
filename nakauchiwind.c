#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double v_beta(double v_inf, double R0, double r, double beta);
void set_r(double R0, double Rmax, int rbin, double r[]);
void set_Mr(double Mc, double Mtot, int rbin, double Mr[]);


/* inner boundary condition */
const double rhoc = 1.0e9;
const double Tc = 1.0e8;


int main(){
    FILE *op;
    
    double v_inf = 1.e9;
    double R0 = 1.e8;
    double beta = 0.7;
    double r_tmp = 1.e10;
    //printf("%12.3e \n",v_beta(v_inf,R0,r_tmp,beta));


    double Rmax = 1.e13;
    int rbin = 100;
    double r[rbin];
    set_r(R0,Rmax,rbin,r);
    
    double drc = 1.e5;
    double Mc = 4.*M_PI/3.*pow(drc,3.);
    double Mtot = 2.e33;
    double Mr[rbin];
    set_Mr(Mc,Mtot,rbin,Mr);

    
    op = fopen("test.dat","w+");
    int i;
    for (i=0; i<rbin; i++) {
        //fprintf(op,"%12.3e %12.3e \n",r[i],v_beta(v_inf,R0,r[i],beta));
        fprintf(op,"%12.3e %12.3e \n",r[i],Mr[i]);
    }
    return 0;
    fclose(op);
}

/* let me first reproduce Nakauchi et al.18 */

double v_beta(double v_inf, double R0, double r, double beta){
    return v_inf*pow(1.-R0/r,beta);
}

void set_r(double R0, double Rmax, int rbin, double r[]){
    double del_ln_r = log(Rmax/R0)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        r[i] = R0*exp(del_ln_r*(double)i);
    }
}

void set_Mr(double Mc, double Mtot, int rbin, double Mr[]){
    double del_ln_r = log(Mtot/Mc)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        Mr[i] = Mc*exp(del_ln_r*(double)i);
    }
}





