#include <stdio.h>
#include <math.h>

const double arad = 7.5657e-15;
const double C = 2.99792458e10;

double calc_Df_simple(double R);
double calc_Df(double R);
double find_Rfld(double r, double T, double Lr);
double find_dTdr(double rho, double kappa, double T, double Rfld);


int main(){
    double r = 6.0e9;
    double Lr = 2.0e38;
    double rho = 1.0e-10;
    double kappa = 0.2;
    double T = 3.0e5;
    
    double Rfld = find_Rfld(r,T,Lr);
    double dTdr = find_dTdr(rho,kappa,T,Rfld);

    printf("%12.3e %12.3e \n",Rfld,dTdr);
}

double calc_Df_simple(double R)
{
    return (2. + R)/(6. + 3.*R + R*R);
}

double calc_Df(double R)
{
    return (1./tanh(R) - 1./R)/R;
}

double find_Rfld(double r, double T, double Lr)
{
    /* analytic solution of Rfld for given r, T, and Lr */
    double s = Lr/(4.*M_PI*pow(r,2.)*arad*pow(T,4.)*C);
    
    return (3.*s-2+sqrt(4.+12.*s-15.*s*s))/2./(1.-s);
}

double find_dTdr(double rho, double kappa, double T, double Rfld)
{
    return -Rfld*kappa*rho*T/4.;
}
