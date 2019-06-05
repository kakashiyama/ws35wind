# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "rkf45.h"
# include "rkf45.c"

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
const double mu_mol = .5;
const double Mwd = 1.3*Msun;
const double Rwd = 1.e8;
const double Bwd = 3.e8;
const double Omega = .1;
const double Mdot = 3.e-6*Msun/Yr;
const double Fm = Mdot/4./M_PI;
const double FB = Rwd*Rwd*Bwd;

/* optcaity table #46 */
const int index_T=71;
const int index_R=20;

/* trial parameters at rA */
const double dudxA = 3.;
const double rA = 8.2e9;
const double TA = 3.e5;
const double LrA = 2.e38;

/* calculate other parameters at rA */
const double vA = Bwd*Bwd*Rwd*Rwd*Rwd*Rwd/Mdot/rA/rA;
const double rhoA = Fm/vA/rA/rA;
const double BrA = FB/rA/rA;
const double vphiA = rA*Omega*dudxA/(2.+ dudxA);
const double BphiA = -BrA*rA*Omega/vA*(2./(2.+ dudxA));
const double Lang = rA*rA*Omega;
const double kA = 0.5*(vA*vA + vphiA*vphiA);
const double hA = 5./2.*kB*TA/mu_mol/Mu + 4.*arad*TA*TA*TA*TA/3./rhoA;
const double etot = LrA/4./M_PI/Fm + kA + hA - G*Mwd/rA - rA*Omega*vphiA + Lang*Omega;


int main(void);
void test(void);
void r8_f2(double x, double y[], double yp[]);
void solve_constraint_eqs(double r, double vr, double T, double *rho, double *vphi, double *Br, double *Bphi, double *Lr, double *kappa);
void calc_derivatives(double r, double vr, double T, double rho, double vphi, double Br, double Bphi, double Lr, double kappa, double *dvrdr, double *dTdr);
void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T, double log10rho, double kappa_tab[index_T][index_R]);
double solve_Rfld(double r, double T, double Lr);
double solve_dTdr(double rho, double kappa, double T, double Rfld);


int main()
{
  timestamp();

  test();
/*
  Terminate.
*/
  printf("\n");
  printf("RKF45\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  timestamp();

  return 0;
}


void test( )
{
# define NEQN 2

    double abserr;
    int flag;
    int i_step;
    int n_step;
    double relerr;
    double x;
    double x_out;
    double x_start;
    double x_stop;
    double y[NEQN];
    double yp[NEQN];
    
    double r,vr,T;
    double rho,vphi,Br,Bphi,Lr,kappa;
    double dvrdr,dTdr;

    abserr = 1.;//sqrt(r8_epsilon());
    relerr = 1.;//sqrt(r8_epsilon());
    
    flag = 1;
    
    x_start = 1.;
    x_stop = 1.5;
    
    n_step = 1000;
    
    x = 1.;
    x_out = 1.;
    
    y[0] = 1.;
    y[1] = 1.;
    r8_f2(x,y,yp);
    r = x*rA;
    vr = y[0]*vA;
    T = y[1]*TA;
    solve_constraint_eqs(r,vr,T,&rho,&vphi,&Br,&Bphi,&Lr,&kappa);
    calc_derivatives(r,vr,T,rho,vphi,Br,Bphi,Lr,kappa,&dvrdr,&dTdr);
    
    FILE *op;
    op = fopen("test.dat","w");
    fprintf(op,"%4d %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e \n",
            flag,x,y[0],y[1],yp[0],yp[1],r,vr,T,rho,vphi,Br,Bphi,Lr,kappa,dvrdr,dTdr);
    
    for ( i_step = 1; i_step <= n_step; i_step++ )
    {
        x = ( ( double ) ( n_step - i_step + 1 ) * x_start
            + ( double ) (          i_step - 1 ) * x_stop )
            / ( double ) ( n_step              );
        
        x_out = ( ( double ) ( n_step - i_step ) * x_start
                + ( double ) (          i_step ) * x_stop )
                / ( double ) ( n_step          );

        flag = r8_rkf45(r8_f2,NEQN,y,yp,&x,x_out,&relerr,abserr,flag);

        r = x*rA;
        vr = y[0]*vA;
        T = y[1]*TA;
        solve_constraint_eqs(r,vr,T,&rho,&vphi,&Br,&Bphi,&Lr,&kappa);
        calc_derivatives(r,vr,T,rho,vphi,Br,Bphi,Lr,kappa,&dvrdr,&dTdr);
        
        fprintf(op,"%4d %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e \n",
                flag,x,y[0],y[1],yp[0],yp[1],r,vr,T,rho,vphi,Br,Bphi,Lr,kappa,dvrdr,dTdr);
    }
    
    fclose(op);
    
    return;
# undef NEQN
}


void r8_f2(double x, double y[], double yp[])
{
    double u = y[0];
    double t = y[1];
    
    double r = x*rA;
    double vr = u*vA;
    double T = t*TA;
    
    double rho,vphi,Br,Bphi,Lr,kappa;
    solve_constraint_eqs(r,vr,T,&rho,&vphi,&Br,&Bphi,&Lr,&kappa);

    double dvrdr,dTdr;
    calc_derivatives(r,vr,T,rho,vphi,Br,Bphi,Lr,kappa,&dvrdr,&dTdr);
    
    if (x == 1.){
        yp[0] = dudxA;
    } else {
        yp[0] = dvrdr*(rA/vA);
    }
    yp[1] = dTdr*(rA/TA);

    return;
}


void solve_constraint_eqs(double r, double vr, double T, double *rho, double *vphi, double *Br, double *Bphi, double *Lr, double *kappa)
{
    double x = r/rA;
    double u = vr/vA;
    double t = T/TA;
    
    double rho_tmp = rhoA/u/x/x;
    double Br_tmp = BrA/x/x;
    double vphi_tmp,Bphi_tmp;
    if (x == 1.){
        vphi_tmp = vphiA;
        Bphi_tmp = BphiA;
    } else {
        vphi_tmp = rA*Omega*x*(1.-u)/(1.-x*x*u);
        Bphi_tmp = -BrA*rA*Omega/vA*x*(1.-x*x)/(1.-x*x*u);
    }
    
    double h = 5./2.*kB*T/mu_mol/Mu + 4.*arad*pow(T,4.)/3./rho_tmp;
    double k = 0.5*(vr*vr + vphi_tmp*vphi_tmp);
    double Lr_tmp = 4.*M_PI*Fm*(etot - k - h + G*Mwd/r + r*Omega*vphi_tmp - Lang*Omega);
    
    /* load kappa table */
    double kappa_tab[index_T][index_R];
    load_kappa_table(kappa_tab);
    double kappa_tmp = kappa_fit(log10(T),log10(rho_tmp),kappa_tab);

    *rho = rho_tmp;
    *vphi = vphi_tmp;
    *Br = Br_tmp;
    *Bphi = Bphi_tmp;
    *Lr = Lr_tmp;
    *kappa = kappa_tmp;
    
    return ;
}


void calc_derivatives(double r, double vr, double T, double rho, double vphi, double Br, double Bphi, double Lr, double kappa, double *dvrdr, double *dTdr)
{
    double Ar = Br/sqrt(4.*M_PI*rho);
    double Aphi = Bphi/sqrt(4.*M_PI*rho);
    
    double denominator_of_dvrdr = (vr*vr - kB*T/mu_mol/Mu - Aphi*Aphi*vr*vr/(vr*vr-Ar*Ar))*r/vr;
    
    double dPdr_term = 3.*kappa*Lr/16./M_PI/arad/C/pow(T,3.)/r*(rho*kB/mu_mol/Mu + 4./3.*arad*pow(T,3.)) + 2.*kB*T/mu_mol/Mu;
    double gravity_term = -G*Mwd/r;
    double centrifugal_force_term = vphi*vphi;
    double magnetic_term = 2.*vr*vphi*Ar*Aphi/(vr*vr-Ar*Ar);
    double numerator_of_dvrdr = dPdr_term + gravity_term + centrifugal_force_term + magnetic_term;
    
    double Rfld = solve_Rfld(r,T,Lr);
    
    if (numerator_of_dvrdr*denominator_of_dvrdr < 0.)
        exit(1);
    
    *dvrdr = numerator_of_dvrdr/denominator_of_dvrdr;
    *dTdr = solve_dTdr(rho,kappa,T,Rfld);
    
    return ;
}


void load_kappa_table(double kappa_tab[index_T][index_R])
{
    int i,j;
    
    FILE *ip;
    ip = fopen("../tab46.dat","r");
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
    while(kappa_tab[i][0]<log10T && i<index_T)
        i++;
    index_T_fit = i;
    
    if (index_T_fit == 1){
        //printf("Warning (T too small for opal kappa): log10T = %le --> %le , log10R = %le \n",log10T,kappa_tab[1][0],log10R);
        index_T_fit = 2;
        log10T = kappa_tab[1][0];
    } else if (index_T_fit == index_T){
        //printf("Warning (T too large for opal kappa): log10T = %le --> %le , log10R = %le \n",log10T,kappa_tab[index_T-1][0],log10R);
        log10T = kappa_tab[index_T-1][0];
    }
    
    j=1;
    while(kappa_tab[0][j]<log10R && j<index_R)
        j++;
    index_R_fit = j;
    
    if (index_R_fit == 1){
        //printf("Warning (R too small for opal kappa): log10R = %le --> %le , log10T = %le, log10rho = %le \n",log10R,kappa_tab[0][1],log10T,log10rho);
        index_R_fit = 2;
        log10R = kappa_tab[0][1];
    } else if (index_R_fit == index_R){
        //printf("Warning (R too large got opsl kappa): log10R = %le --> %le , log10T = %le \n",log10R,kappa_tab[0][index_R-1],log10T);
        log10R = kappa_tab[0][index_R-1];
    }
    
    double kappa_Tdirect_min = (kappa_tab[index_T_fit-1][index_R_fit]-kappa_tab[index_T_fit-1][index_R_fit-1])/(kappa_tab[0][index_R_fit]-kappa_tab[0][index_R_fit-1])*(log10R-kappa_tab[0][index_R_fit-1])+kappa_tab[index_T_fit-1][index_R_fit-1];
    double kappa_Tdirect_max = (kappa_tab[index_T_fit][index_R_fit]-kappa_tab[index_T_fit][index_R_fit-1])/(kappa_tab[0][index_R_fit]-kappa_tab[0][index_R_fit-1])*(log10R-kappa_tab[0][index_R_fit-1])+kappa_tab[index_T_fit][index_R_fit-1];
    double log10kappa = (kappa_Tdirect_max-kappa_Tdirect_min)/(kappa_tab[index_T_fit][0]-kappa_tab[index_T_fit-1][0])*(log10T-kappa_tab[index_T_fit-1][0])+kappa_Tdirect_min;
    
    //return 0.2;
    return pow(10.,log10kappa);
}


double solve_Rfld(double r, double T, double Lr)
{
    /* analytic solution of Rfld for given r, T, and Lr */
    double s = Lr/(4.*M_PI*pow(r,2.)*arad*pow(T,4.)*C);
    
    return (3.*s-2+sqrt(4.+12.*s-15.*s*s))/2./(1.-s);
}


double solve_dTdr(double rho, double kappa, double T, double Rfld)
{
    return -Rfld*kappa*rho*T/4.;
}
