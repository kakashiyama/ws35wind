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
const double mu_mol = 0.5;
const double Mwd = 1.*Msun;
const double Rwd = 3.0e8;
const double Bwd = 3.0e8;
const double Omega = 0.1;
const double Mdot = 3.0e-6*Msun/Yr;

/* number of radial bin */
const int rbin = 10000;

/* optcaity table #46 */
const int index_T=58;
const int index_R=20;

void set_r_from_rA_to_rWD(double rA, double r[]);
void set_para_at_rA(double rA_tmp, double vA_tmp, double dudxA_tmp, double TA_tmp, double LrA_tmp,
                    double *BrA, double *rhoA, double *vphiA, double *BphiA, double *Fm, double *FB, double *Lang, double *etot);
void calc_dTdr(double Tinput, double *Toutput, double dr, double r, double rho, double Lr, double kappa);
void calc_dVrdr_1ststep(double vrinput, double *vroutput, double dr, double dudxA, double vA, double rA);
void calc_dVrdr_2ststep_and_more(double vrinput, double r, double dr, double rho, double Br, double Bphi, double vphi, double T, double Lr, double kappa,
                                 double *vroutput, double *denominator_of_dvrdr, double *numerator_of_dvrdr);
void solve_constraint_eqs(double r, double rA, double vA, double vr, double T, double Fm, double FB, double Lang, double etot,
                          double *Br, double *Bphi, double *vphi, double *Lr, double *rho, double *a);

double michel_wind_velocity(double r, double br, double omega, double mdot);
void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T, double log10rho, double kappa_tab[index_T][index_R]);


int main()
{
    /* load kappa table */
    double kappa_tab[index_T][index_R],kappa[rbin];
    load_kappa_table(kappa_tab);
    
    /* fixed parameters */
    double vA = .15*michel_wind_velocity(Rwd,Bwd,Omega,Mdot); /* normalized by the Michel velocity */
    double dudxA = 1.0;
    double TA = 3.e5;
    double LrA = 2.e38;
    
    /* to-be-calculated parameters */
    double BrA,rhoA,vphiA,BphiA,Fm,FB,Lang,etot;
    double r[rbin],dr;
    double T[rbin],vr[rbin],Br[rbin],Bphi[rbin],vphi[rbin],Lr[rbin],rho[rbin],a,denominator_of_dvrdr[rbin],numerator_of_dvrdr[rbin];

    
    /* trial parameters at rA */
    double rAmax = 2.*Rsun;
    double rAmin = 0.01*Rsun;
    double rA = .5*(rAmax+rAmin);
    int max_trial = 50,j=0;
    
    while (j<max_trial){
        printf("%d: rA = %lf \n",j,rA/Rsun);
        
        /* set radial coordinate */
        set_r_from_rA_to_rWD(rA,r);
        
        /* calculate other parameters at rA */
        set_para_at_rA(rA,vA,dudxA,TA,LrA,&Br[0],&rho[0],&vphi[0],&Bphi[0],&Fm,&FB,&Lang,&etot);
        vr[0] = vA;
        T[0] = TA;
        Lr[0] = LrA;
        kappa[0] = kappa_fit(log10(T[0]),log10(rho[0]),kappa_tab);
        
        /* calculate the 1st step */
        dr = r[1]-r[0];
        calc_dTdr(T[0],&T[1],dr,r[0],rho[0],Lr[0],kappa[0]);
        calc_dVrdr_1ststep(vr[0],&vr[1],dr,dudxA,vA,rA);
        solve_constraint_eqs(r[1],rA,vA,vr[1],T[1],Fm,FB,Lang,etot,&Br[1],&Bphi[1],&vphi[1],&Lr[1],&rho[1],&a);
        kappa[1] = kappa_fit(log10(T[1]),log10(rho[1]),kappa_tab);
        
        /* calculate the 2nd step and more */
        int i=1;
        while (i<rbin-1 && vr[i] > 0.){
            dr = r[i+1]-r[i];
            calc_dTdr(T[i],&T[i+1],dr,r[i],rho[i],Lr[i],kappa[i]);
            calc_dVrdr_2ststep_and_more(vr[i],r[i],dr,rho[i],Br[i],Bphi[i],vphi[i],T[i],Lr[i],kappa[i],&vr[i+1],&denominator_of_dvrdr[i+1],&numerator_of_dvrdr[i+1]);
            solve_constraint_eqs(r[i+1],rA,vA,vr[i+1],T[i+1],Fm,FB,Lang,etot,&Br[i+1],&Bphi[i+1],&vphi[i+1],&Lr[i+1],&rho[i+1],&a);
            kappa[i+1] = kappa_fit(log10(T[i+1]),log10(rho[i+1]),kappa_tab);
            i++;
        }
        
        
        if (i >= rbin-1) {
            if (numerator_of_dvrdr[i]/denominator_of_dvrdr[i] > 0.){
                break;
            } else {
                rAmin = rA;
                rA = 0.5*(rA+rAmax);
            }
        } else {
            rAmax = rA;
            rA = 0.5*(rA+rAmin);
        }
        
        j++;
    
    }
    
    FILE *op;
    op = fopen("test.dat","w");
    for (j=0; j<rbin; j++) {
        fprintf(op,"%lf %lf \n",r[j],vr[j]);
    }
    fclose(op);
    
    if (j<max_trial){
        printf("find a solution \n");
    } else {
        printf("could not find a solution\n");
    }
    
    return 0;
}


void set_r_from_rA_to_rWD(double rA, double r[])
{
    double del_ln_r = log(Rwd/rA)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        r[i] = rA*exp(del_ln_r*(double)i);
    }
}

void set_para_at_rA(double rA, double vA, double dudxA, double TA, double LrA,
                    double *BrA, double *rhoA, double *vphiA, double *BphiA, double *Fm, double *FB, double *Lang, double *etot)
{
    double BrA_tmp = pow(rA/Rwd,-2.0)*Bwd;
    double rhoA_tmp = Mdot/4./M_PI/vA/pow(rA,2.0);
    double vphiA_tmp = rA*Omega*dudxA/(2.+ dudxA);
    double BphiA_tmp = -BrA_tmp*rA*Omega/vA*(2./(2.+ dudxA));
    
    double hA = 5./2.*kB*TA/mu_mol/Mu + 4.*arad*pow(TA,4.)/3./rhoA_tmp;
    double kA = 0.5*(vA*vA + vphiA_tmp*vphiA_tmp);
    double Fm_tmp = rhoA_tmp*vA*rA*rA;
    double Lang_tmp = rA*rA*Omega;
    double etot_tmp = LrA/4./M_PI/Fm_tmp + kA + hA - G*Mwd/rA - rA*Omega*vphiA_tmp + Lang_tmp*Omega;
   
    *BrA = BrA_tmp;
    *rhoA = rhoA_tmp;
    *vphiA = vphiA_tmp;
    *BphiA = BphiA_tmp;
    
    *FB = Rwd*Rwd*Bwd;
    *Fm = Fm_tmp;
    *Lang = Lang_tmp;
    *etot = etot_tmp;
}


void calc_dTdr(double Tinput, double *Toutput, double dr, double r, double rho, double Lr, double kappa)
{
    *Toutput = Tinput - 3.*kappa*rho*Lr/16./M_PI/arad/C/pow(Tinput,3.)/pow(r,2.)*dr;
}


void calc_dVrdr_1ststep(double vrinput, double *vroutput, double dr, double dudxA, double vA, double rA)
{
    *vroutput = vrinput + dr*dudxA*(vA/rA);
}


void calc_dVrdr_2ststep_and_more(double vrinput, double r, double dr, double rho, double Br, double Bphi, double vphi, double T, double Lr, double kappa,
                                 double *vroutput, double *denominator_of_dvrdr, double *numerator_of_dvrdr)
{
    double Ar = Br/sqrt(4.*M_PI*rho);
    double Aphi = Bphi/sqrt(4.*M_PI*rho);
    
    double denominator_of_dvrdr_tmp = (vrinput*vrinput -kB*T/mu_mol/Mu - Aphi*Aphi*vrinput*vrinput/(vrinput*vrinput-Ar*Ar))*r/vrinput;
    double dPdr_term = 3.*kappa*Lr/16./M_PI/arad/C/pow(T,3.)/r*(rho*kB/mu_mol/Mu + 4./3.*arad*pow(T,3.)) + 2.*kB*T/mu_mol/Mu;
    double gravity_term = -G*Mwd/r;
    double centrifugal_force_term = vphi*vphi;
    double magnetic_term = 2.*vrinput*vphi*Ar*Aphi/(vrinput*vrinput-Ar*Ar);
    double numerator_of_dvrdr_tmp = dPdr_term + gravity_term + centrifugal_force_term + magnetic_term;
    
    *vroutput = vrinput + dr*numerator_of_dvrdr_tmp/denominator_of_dvrdr_tmp;
    *denominator_of_dvrdr = denominator_of_dvrdr_tmp;
    *numerator_of_dvrdr = numerator_of_dvrdr_tmp;
}


void solve_constraint_eqs(double r, double rA, double vA, double vr, double T, double Fm, double FB, double Lang, double etot,
                          double *Br, double *Bphi, double *vphi, double *Lr, double *rho, double *a)
{
    double Br_tmp = FB/r/r;
    double rho_tmp = Fm/vr/r/r;
    double u = vr/vA;
    double x = r/rA;
    double vphi_tmp = rA*Omega*x*(1.-u)/(1.-x*x*u);
    double Bphi_tmp = -Br_tmp*rA*Omega/vA*x*(1.-x*x)/(1.-x*x*u);
    double a_tmp = sqrt(kB*T/mu_mol/Mu);

    
    double h_tmp = 5./2.*kB*T/mu_mol/Mu + 4.*arad*pow(T,4.)/3./rho_tmp;
    double k_tmp = 0.5*(vr*vr+vphi_tmp*vphi_tmp);
    double Lr_tmp = 0.;//4.*M_PI*Fm*(etot-k_tmp-h_tmp+G*Mwd/r+r*Omega*vphi_tmp-Lang*Omega);
    
    *Br = Br_tmp;
    *Bphi = Bphi_tmp;
    *vphi = vphi_tmp;
    *Lr = Lr_tmp;
    *rho = rho_tmp;
    *a = a_tmp;
}


double michel_wind_velocity(double r, double br, double omega, double mdot)
{
    return pow(pow(r,4.)*pow(br,2.)*pow(omega,2.)/mdot,1./3.);
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
