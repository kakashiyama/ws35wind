#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "const.h"

const int n_input_para = 8;

struct _input load_inputfile();
struct _fixed calc_fixed_para(struct _input in, double rA, double dudxA);

struct _input{
    double mu_mol;
    double Mwd;
    double Rwd;
    double Bwd;
    double Omega;
    double Mdot;
    double TA;
    double LrA;
};

struct _fixed{
    double Fm;
    double FB;
    double vA;
    double rhoA;
    double BrA;
    double vphiA;
    double BphiA;
    double Lang;
    double kA;
    double hA;
    double etotA;
};

int main()
{
    /* load input file */
    struct _input in = load_inputfile();
    printf("%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",
           in.mu_mol,in.Mwd,in.Rwd,in.Bwd,in.Omega,in.Mdot,in.TA,in.LrA);
    
    /* trial parameters for each shooting */
    double rA = 0.1*Rsun;
    double dudxA = 1.;
    
    /* calculate parameters fixed during the shooting */
    struct _fixed fix = calc_fixed_para(in,rA,dudxA);
    printf("%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",
           fix.Fm,fix.FB,fix.vA,fix.rhoA,fix.BrA,fix.vphiA,fix.BphiA,fix.Lang,fix.kA,fix.hA,fix.etotA);
    
    return 0;
}

struct _input load_inputfile()
{
    double dmy[n_input_para];
    FILE *ip;
    ip = fopen("input.dat","r");
    int i=0;
    while (fscanf(ip,"%*s %lf %*[\n]",&dmy[i])!=EOF){
        i++;
    }
    fclose(ip);
    
    struct _input in;
    in.mu_mol = dmy[0];
    in.Mwd = dmy[1]*Msun;
    in.Rwd = dmy[2];
    in.Bwd = dmy[3];
    in.Omega = dmy[4];
    in.Mdot = dmy[5]*Msun/Yr;
    in.TA = dmy[6];
    in.LrA = dmy[7];
    
    return in;
}

struct _fixed calc_fixed_para(struct _input in, double rA, double dudxA)
{
    struct _fixed fix;
    
    fix.Fm = in.Mdot/4./M_PI;
    fix.FB = in.Bwd*pow(in.Rwd,2.);
    fix.vA = pow(in.Bwd,2.)*pow(in.Rwd,4.)/in.Mdot/pow(rA,2.);
    fix.rhoA = fix.Fm/fix.vA/rA/rA;
    fix.BrA = fix.FB/rA/rA;
    fix.vphiA = rA*in.Omega*dudxA/(2.+ dudxA);
    fix.BphiA = -fix.BrA*rA*in.Omega/fix.vA*(2./(2.+ dudxA));
    fix.Lang = rA*rA*in.Omega;
    fix.kA = 0.5*(pow(fix.vA,2.) + pow(fix.vphiA,2.));
    fix.hA = 5./2.*kB*in.TA/in.mu_mol/Mu;
    fix.etotA = in.LrA/4./M_PI/fix.Fm + fix.kA + fix.hA - G*in.Mwd/rA - rA*in.Omega*fix.vphiA + fix.Lang*in.Omega;

    return fix;
}

