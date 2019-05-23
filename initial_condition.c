#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* physical constant */
const double mol = 6.02e23;
const double joul_to_erg = 1.e7;
const double Rgas = 8.314462618*joul_to_erg/mol;
const double arad = 7.5657e-15;
const double G = 6.67408e-8;
const double C = 2.99792458e10;
const double Msun = 2.e33;
const double Rsun = 6.9551e10;
const double Yr = 365.*24.*60.*60.;

/* model parameters */
const double mu = 2.0;
const double Mwd = 1.0*Msun;
const double Rwd = 1.0e8;
const double Bwd = 1.0e9;
const double Omega = 0.01;
const double Mdot = 1.0e6*Msun/Yr;

void set_initial_guess_at_rA(double rA_tmp, double vA_tmp, double dudxA_tmp, double *BrA, double *rhoA, double *vphiA);

int main()
{
    double rA_tmp=0.1*Rsun;
    double fac_W=1.5;
    double dudxA_tmp=1.5;
    
    //double TA_tmp=1.0e5;

    double vA_tmp=rA_tmp*Omega/fac_W;
    
    double BrA,rhoA,vphiA;
    
    set_initial_guess_at_rA(rA_tmp,vA_tmp,dudxA_tmp,&BrA,&rhoA,&vphiA);
    printf("%12.3e %12.3e %12.3e %12.3e \n",vA_tmp,BrA,rhoA,vphiA);
    return 0;
}

/* rAとdudxAとLambda、Tを振って... */

void set_initial_guess_at_rA(double rA_tmp, double vA_tmp, double dudxA_tmp, double *BrA, double *rhoA, double *vphiA)
{
    double BrA_tmp = pow(rA_tmp/Rwd,-2.0)*Bwd;
    double rhoA_tmp = Mdot/4./M_PI/vA_tmp/pow(rA_tmp,2.0);
    double vphiA_tmp = rA_tmp*Omega*dudxA_tmp/(2.+dudxA_tmp);
    
    *BrA = BrA_tmp;
    *rhoA = rhoA_tmp;
    *vphiA = vphiA_tmp;
    
}
