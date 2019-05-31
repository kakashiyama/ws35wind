#include <stdio.h>
#include <math.h>

const double arad = 7.5657e-15;
const double C = 2.99792458e10;

double calc_Rfld(double T, double dTdr, double kappa, double rho);
double calc_Df_simple(double R);
double calc_Df(double R);
double calc_dTdr(double r, double T, double kappa, double rho, double Lr, double Rfld);
double find_Rfld(double r, double T, double Lr);

int main(){
    double r = 6.0e9;
    double Lr = 2.0e38;
    double rho = 1.0e-10;
    double kappa = 0.2;
    double T = 3.0e5;
    double dTdr = -1.0e-4;
    
    double Rfld = calc_Rfld(T,dTdr,kappa,rho);
    double Df_simple = calc_Df_simple(Rfld);
    double Df = calc_Df(Rfld);
    double dTdr_tmp = calc_dTdr(r,T,kappa,rho,Lr,1./3.);
    double tmp = find_Rfld(r,T,Lr);

    printf("%12.3e %12.3e %12.3e %12.3e %12.3e \n",Rfld,Df_simple,Df,dTdr_tmp,tmp);
}

double calc_Rfld(double T, double dTdr, double kappa, double rho)
{
    return -dTdr/T/kappa/rho;
}

double calc_Df_simple(double R)
{
    return (2. + R)/(6. + 3.*R + R*R);
}

double calc_Df(double R)
{
    return (1./tanh(R) - 1./R)/R;
}

double calc_dTdr(double r, double T, double kappa, double rho, double Lr, double Rfld)
{
    return -kappa*rho/16./M_PI/arad/C/pow(T,3.)/calc_Df(Rfld)*Lr/pow(r,2.);
}

double find_Rfld(double r, double T, double Lr)
{
    double flux_density = Lr/(4.*M_PI*pow(r,2.)*arad*pow(T,4.)*C);
    double flux_density_min = 0.;
    double flux_density_max = 1.;
    double flux_density_tmp = 0.5*(flux_density_min + flux_density_max);
    
    /* under construction */
    
    return flux_density;
}
