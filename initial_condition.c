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
const double kB = 1.38064852e-16;
const double Mu = 1.6726219e-24;
const double Msun = 2.e33;
const double Rsun = 6.9551e10;
const double Yr = 365.*24.*60.*60.;

/* model parameters */
const double mu_mol = 2.0;
const double Mwd = 1.0*Msun;
const double Rwd = 3.0e8;
const double Bwd = 1.0e9;
const double Omega = 0.01;
const double Mdot = 1.0e-6*Msun/Yr;

/* see table #46 */
const int rbin = 1000;
const int index_T=58;
const int index_R=20;

double michel_wind_velocity(double r, double br, double omega, double mdot);
void set_r_from_rA_to_rWD(double rA, double r[]);
void set_para_at_rA(double rA_tmp, double vA_tmp, double dudxA_tmp, double TA_tmp, double LrA_tmp, double *BrA, double *rhoA, double *vphiA, double *BphiA, double *Fm, double *FB, double *Lang, double *etot);
void solve_constraint_eqs(double r, double rA, double vA, double vr, double T, double Fm, double FB, double Lang, double etot, double *Br, double *Bphi, double *vphi, double *Lr, double *rho);

void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T6, double log10rho, double kappa_tab[index_T][index_R]);
void calc_dTdr(double Tinput, double *Toutput, double dr, double r, double rho, double Lr);
void calc_dVrdr_1ststep(double vrinput, double *vroutput, double dr, double dudxA, double vA, double rA);



int main()
{
    double VM = michel_wind_velocity(Rwd,Bwd,Omega,Mdot);
    printf("VM = %12.3e cm/s \n",VM);
    
    /* trial parameters */
    double rA_tmp = 0.1*Rsun;
    double vA_tmp = 0.1*VM;
    double dudxA_tmp = 1.2;
    double TA_tmp = 1.0e5;
    double LrA_tmp = 1.e38;
    
    /* to-be-caluclated parameters */
    double BrA,rhoA,vphiA,BphiA,Fm,FB,Lang,etot;
    set_para_at_rA(rA_tmp,vA_tmp,dudxA_tmp,TA_tmp,LrA_tmp,&BrA,&rhoA,&vphiA,&BphiA,&Fm,&FB,&Lang,&etot);
    printf("%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",BrA,rhoA,vphiA,BphiA,Fm,FB,Lang,etot);
    
    double radius[rbin];
    set_r_from_rA_to_rWD(rA_tmp,radius);
    
    
    /* calculate the 1st dTdr & dVrdr */
    double Toutput,vroutput;
    double dr = radius[1]-radius[0];
    calc_dTdr(TA_tmp,&Toutput,dr,rA_tmp,rhoA,LrA_tmp);
    calc_dVrdr_1ststep(vA_tmp,&vroutput,dr,dudxA_tmp,vA_tmp,rA_tmp);
    
    printf("%12.3e %12.3e %12.3e %12.3e \n",TA_tmp,Toutput,vA_tmp,vroutput);

    
    
    return 0;
}

double michel_wind_velocity(double r, double br, double omega, double mdot)
{
    return pow(pow(r,4.)*pow(br,2.)*pow(omega,2.)/mdot,1./3.);
}

void set_r_from_rA_to_rWD(double rA, double r[])
{
    double del_ln_r = log(Rwd/rA)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        r[i] = rA*exp(del_ln_r*(double)i);
    }
}

void set_para_at_rA(double rA_tmp, double vA_tmp, double dudxA_tmp, double TA_tmp, double LrA_tmp, double *BrA, double *rhoA, double *vphiA, double *BphiA, double *Fm, double *FB, double *Lang, double *etot)
{
    double BrA_tmp = pow(rA_tmp/Rwd,-2.0)*Bwd;
    double rhoA_tmp = Mdot/4./M_PI/vA_tmp/pow(rA_tmp,2.0);
    double vphiA_tmp = rA_tmp*Omega*dudxA_tmp/(2.+dudxA_tmp);
    double BphiA_tmp = -BrA_tmp*rA_tmp*Omega/vA_tmp*(2./(2.+dudxA_tmp));
    
    double hA_tmp = 5./2.*kB*TA_tmp/mu_mol/Mu + 4.*arad*pow(TA_tmp,4.)/3./rhoA_tmp;
    double kA_tmp = 0.5*(vA_tmp*vA_tmp+vphiA_tmp*vphiA_tmp);
    double Fm_tmp = rhoA_tmp*vA_tmp*rA_tmp*rA_tmp;
    double Lang_tmp = rA_tmp*rA_tmp*Omega;
    double etot_tmp = LrA_tmp/4./M_PI/Fm_tmp + kA_tmp + hA_tmp - G*Mwd/rA_tmp - rA_tmp*Omega*vphiA_tmp + Lang_tmp*Omega;
   
    *BrA = BrA_tmp;
    *rhoA = rhoA_tmp;
    *vphiA = vphiA_tmp;
    *BphiA = BphiA_tmp;
    
    *FB = Rwd*Rwd*Bwd;
    *Fm = Fm_tmp;
    *Lang = Lang_tmp;
    *etot = etot_tmp;
}

void calc_dTdr(double Tinput, double *Toutput, double dr, double r, double rho, double Lr)
{
    double kappa_tab[index_T][index_R],kappa;
    load_kappa_table(kappa_tab);
    kappa = kappa_fit(log10(Tinput),log10(rho),kappa_tab);
    
    *Toutput = Tinput - 3.*kappa*rho*Lr/16./M_PI/arad/C/pow(Tinput,3.)/pow(r,2.)*dr;
}

void calc_dVrdr_1ststep(double vrinput, double *vroutput, double dr, double dudxA, double vA, double rA)
{
    *vroutput = vrinput + dudxA*vA/rA*dr;
}

void solve_constraint_eqs(double r, double rA, double vA, double vr, double T, double Fm, double FB, double Lang, double etot, double *Br, double *Bphi, double *vphi, double *Lr, double *rho)
{
    double Br_tmp = FB/r/r;
    double rho_tmp = Fm/vr/r/r;
    double u = vr/vA;
    double x = r/rA;
    double vphi_tmp = rA*Omega*x*(1.-u)/(1.-x*x*u);
    double Bphi_tmp = -Br_tmp*rA*Omega/vA*(1.+x*x)/(1.-x*x*u);
    
    double h_tmp = 5./2.*kB*T/mu_mol/Mu + 4.*arad*pow(T,4.)/3./rho_tmp;
    double k_tmp = 0.5*(vr*vr+vphi_tmp*vphi_tmp);
    double Lr_tmp = 4.*M_PI*Fm*(etot-k_tmp-h_tmp+G*Mwd/r+r*Omega*vphi_tmp-Lang*Omega);
    
    *Br = Br_tmp;
    *Bphi = Bphi_tmp;
    *vphi = vphi_tmp;
    *Lr = Lr_tmp;
    *rho = rho_tmp;
}


void load_kappa_table(double kappa_tab[index_T][index_R])
{
    int i,j;
    
    FILE *ip;
    ip = fopen("table46.dat","r");
    for (i=0; i<index_T; i++) {
        for (j=0; j<index_R; j++) {
            fscanf(ip,"%le ",&kappa_tab[i][j]);
        }
    }
    fclose(ip);
}


double kappa_fit(double log10T, double log10rho, double kappa_tab[index_T][index_R])
{
    /* bilinear extrapolation of the OPAL opacity table */
    int i,j,index_T_fit=0,index_R_fit=0;
    double log10R = log10rho-3.*(log10T-6.);
    
    i=1;
    while(kappa_tab[i][0]<log10T && i<index_T-1)
        i++;
    index_T_fit = i;
    if (index_T_fit == 1)
        index_T_fit = 2;
    
    j=1;
    while(kappa_tab[0][j]<log10R && j<index_R-1)
        j++;
    index_R_fit = j;
    if (index_R_fit == 1)
        index_R_fit = 2;
    /* tableの外側にはバカ外挿 */
    
    double kappa_Tdirect_min = (kappa_tab[index_T_fit-1][index_R_fit]-kappa_tab[index_T_fit-1][index_R_fit-1])/(kappa_tab[0][index_R_fit]-kappa_tab[0][index_R_fit-1])*(log10R-kappa_tab[0][index_R_fit-1])+kappa_tab[index_T_fit-1][index_R_fit-1];
    double kappa_Tdirect_max = (kappa_tab[index_T_fit][index_R_fit]-kappa_tab[index_T_fit][index_R_fit-1])/(kappa_tab[0][index_R_fit]-kappa_tab[0][index_R_fit-1])*(log10R-kappa_tab[0][index_R_fit-1])+kappa_tab[index_T_fit][index_R_fit-1];
    double log10kappa = (kappa_Tdirect_max-kappa_Tdirect_min)/(kappa_tab[index_T_fit][0]-kappa_tab[index_T_fit-1][0])*(log10T-kappa_tab[index_T_fit-1][0])+kappa_Tdirect_min;
    
    return pow(10.,log10kappa);
}
