#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double v_beta(double v_inf, double R0, double r, double beta);
void set_r(double R0, double Rmax, int rbin, double r[]);
void set_Mr(double Mc, double Mtot, int rbin, double Mr[]);
double Pgas(double mu, double rho, double T);
double Prad(double T);
void dr_dP_and_dLr(double dMr, double Mr, double r, double rho, double epsnu, double *dr, double *dP, double *dLr);
void dT(double r, double T, double rho, double dr, double Lrad, double kappa, double *dT);
void drho(double mu, double T, double rho, double dP, double dT, double *drho);


/* physical constant */
const double mol = 6.02e23;
const double joul_to_erg = 1.e7;
const double Rgas = 8.314462618*joul_to_erg/mol;
const double arad = 7.5657e-15;
const double G = 6.67408e-8;
const double C = 2.99792458e10;
const double Msun = 2.e33;


/* inner boundary condition */
const double rhoc = 1.0e9;
const double Tc = 1.0e8;


/* simplification */
const double mu = 2.;
const double kappa_fix = 0.4;
const double epsnu_fix = 0.;


int main()
{
    FILE *op;
    
    double v_inf = 1.e9;
    double R0 = 1.e8;
    double beta = 0.7;
    double r_tmp = 1.e10;
    //printf("%12.3e \n",v_beta(v_inf,R0,r_tmp,beta));


    double Rmax = 1.e12;
    int rbin = 100;
    double r[rbin];
    set_r(R0,Rmax,rbin,r);
    
    double drc = 1.e5;
    double Mc = 4.*M_PI/3.*pow(drc,3.);
    double Mtot = 10.*Msun;
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

/* static core */
void set_r(double R0, double Rmax, int rbin, double r[])
{
    double del_ln_r = log(Rmax/R0)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        r[i] = R0*exp(del_ln_r*(double)i);
    }
}

void set_Mr(double Mc, double Mtot, int rbin, double Mr[])
{
    double del_ln_r = log(Mtot/Mc)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        Mr[i] = Mc*exp(del_ln_r*(double)i);
    }
}


double Pgas(double mu, double rho, double T)
{
    return Rgas/mu*rho*T;
}

double Prad(double T)
{
    return arad*pow(T,4.)/3.;
}

void dr_dP_and_dLr(double dMr, double Mr, double r, double rho, double epsnu, double *dr, double *dP, double *dLr)
{
    *dr = 1./4./M_PI/r/r/rho*dMr;
    *dP = -G*Mr/4./M_PI/pow(r,4.);
    *dLr = epsnu*dMr;
}

void dT(double r, double T, double rho, double dr, double Lrad, double kappa, double *dT)
{
    *dT = -Lrad/(16.*M_PI*arad*C*r*r*pow(T,3.)/3./kappa/rho);
}

void drho(double mu, double T, double rho, double dP, double dT, double *drho)
{
    *drho = mu/Rgas/T*(dP-dT*(Rgas/mu*rho+4./3.*arad*pow(T,3.)));
}




/* wind region */
double v_beta(double v_inf, double R0, double r, double beta)
{
    return v_inf*pow(1.-R0/r,beta);
}
