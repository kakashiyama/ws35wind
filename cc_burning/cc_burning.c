#include <stdio.h>
#include <math.h>

double calc_epsCC(double T, double rho);

int main()
{
  double rho = 1.0e5;
  double T = 1.0e9;
  double epsCC = calc_epsCC(T,rho);
  printf("%12.3e \n",epsCC);

  return 0;
}

double calc_epsCC(double T, double rho)
{
  double f_CC = 1.;
  double X12 = 0.5;
  double T9 = T/1.e9;
  double T9alpha = T9/(1.+0.067*T9);

  return 5.49e43*f_CC*rho*pow(X12,2.)*pow(T9,-3./2.)*pow(T9alpha,5./6.)*exp(-84.165/pow(T9alpha,1./3.))/(exp(-0.01*pow(T9alpha,4.)) + 5.56e-3*exp(1.685*pow(T9alpha,2./3.)));
}
