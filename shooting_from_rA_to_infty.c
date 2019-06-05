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
const double Mwd = 1.3*Msun;
const double Rwd = 1.0e8;
const double Bwd = 3.0e8;
const double Omega = 0.1;
const double Mdot = 3.0e-6*Msun/Yr;

/* number of radial bin */
const int rbin = 1000000;

/* shooting trial number */
const int max_1sttrial = 50;
const int max_2ndtrial = 50;

/* optcaity table #46 */
const int index_T=71;
const int index_R=20;

void trial_1st_stage(double rA, double rmax, double TA, double LrA, double ln_dudxAmax, double ln_dudxAmin, double *dudxA_1st, double *ddudxA_1st);

void set_r_from_rA_to_infty(double rA, double rmax, double r[]);
void set_para_at_rA(double rA, double dudxA, double TA, double LrA,
                    double *vA, double *BrA, double *rhoA, double *vphiA, double *BphiA, double *Fm, double *FB, double *Lang, double *etot);
void calc_dTdr(double Tinput, double dr, double dTdr,
               double *Toutput);
void calc_dVrdr_1ststep(double vrinput, double dr, double dudxA, double vA, double rA,
                        double *vroutput);
void calc_dVrdr_2ndstep_and_more(double vrinput, double r, double dr, double rho, double Br, double Bphi, double vphi, double T, double Lr, double kappa,
                                 double *vroutput, double *denominator_of_dvrdr, double *numerator_of_dvrdr);
void solve_constraint_eqs(double r, double vr, double T,
                          double rA, double vA, double Fm, double FB, double Lang, double etot,
                          double *Br, double *Bphi, double *vphi, double *Lr, double *rho);

void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T, double log10rho, double kappa_tab[index_T][index_R]);
double solve_Rfld(double r, double T, double Lr);
double solve_dTdr(double rho, double kappa, double T, double Rfld);

double michel_wind_velocity(double r, double br, double omega, double mdot);


int main()
{
    /* input parameters */
<<<<<<< Updated upstream
    double rA = 8.2e9;
    double rmax = rA*100.0;
=======
    double rA = 8185856364.370256;
    double rmax = rA*1.05;
>>>>>>> Stashed changes
    double TA = 3.e5;
    double LrA = 2.e38;
    
    /* trial 1st stage */
    double ln_dudxAmax = log(10.);
    double ln_dudxAmin = log(1.e-3);
    double dudxA,ddudxA;
    trial_1st_stage(rA,rmax,TA,LrA,ln_dudxAmax,ln_dudxAmin,&dudxA,&ddudxA);
    
    return 0;
}


void trial_1st_stage(double rA, double rmax, double TA, double LrA, double ln_dudxAmax, double ln_dudxAmin, double *dudxA_1st, double *ddudxA_1st)
{
    /* initial guess */
    double dudxA = exp(.5*(ln_dudxAmax+ln_dudxAmin));
    
    /* to-be-calculated qunatities */
    double vA,BrA,rhoA,vphiA,BphiA; /* at rA */
    double Fm,FB,Lang,etot; /* conserved quantities */
    double dr,ddudxA,Rfld,dTdr;
    double *r,*vr,*T,*Br,*Bphi,*vphi,*Lr,*rho,*kappa,*denominator_of_dvrdr,*numerator_of_dvrdr;
    r = (double *)malloc(rbin * sizeof(double));
    vr = (double *)malloc(rbin * sizeof(double));
    T = (double *)malloc(rbin * sizeof(double));
    Br = (double *)malloc(rbin * sizeof(double));
    Bphi = (double *)malloc(rbin * sizeof(double));
    vphi = (double *)malloc(rbin * sizeof(double));
    Lr = (double *)malloc(rbin * sizeof(double));
    rho = (double *)malloc(rbin * sizeof(double));
    kappa = (double *)malloc(rbin * sizeof(double));
    denominator_of_dvrdr = (double *)malloc(rbin * sizeof(double));
    numerator_of_dvrdr = (double *)malloc(rbin * sizeof(double));
    
    /* load kappa table */
    double kappa_tab[index_T][index_R];
    load_kappa_table(kappa_tab);
    
    int j=0;
    while (j<max_1sttrial){
        printf("1st trial No. %d: dudxA = %lf \n",j,dudxA);
        
        /* set radial coordinate */
        set_r_from_rA_to_infty(rA,rmax,r);
        
        /* calculate other parameters at rA */
        set_para_at_rA(rA,dudxA,TA,LrA,&vA,&Br[0],&rho[0],&vphi[0],&Bphi[0],&Fm,&FB,&Lang,&etot);
        vr[0] = vA;
        T[0] = TA;
        Lr[0] = LrA;
        kappa[0] = kappa_fit(log10(T[0]),log10(rho[0]),kappa_tab);
        
        /* calculate the 1st step */
        dr = r[1]-r[0];
        Rfld = solve_Rfld(r[0],T[0],Lr[0]);
        dTdr = solve_dTdr(rho[0],kappa[0],T[0],Rfld);
        calc_dTdr(T[0],dr,dTdr,&T[1]);
        calc_dVrdr_1ststep(vr[0],dr,dudxA,vA,rA,&vr[1]);
        denominator_of_dvrdr[1] = -dudxA*vA*vA;
        numerator_of_dvrdr[1] = -vA*rA;
        solve_constraint_eqs(r[1],vr[1],T[1],rA,vA,Fm,FB,Lang,etot,&Br[1],&Bphi[1],&vphi[1],&Lr[1],&rho[1]);
        kappa[1] = kappa_fit(log10(T[1]),log10(rho[1]),kappa_tab);
        
        /* calculate the 2nd step and more */
        int i=1;
        while (i<rbin-1 && denominator_of_dvrdr[i]*numerator_of_dvrdr[i] > 0.){
            dr = r[i+1]-r[i];
            Rfld = solve_Rfld(r[i],T[i],Lr[i]);
            dTdr = solve_dTdr(rho[i],kappa[i],T[i],Rfld);
            calc_dTdr(T[i],dr,dTdr,&T[i+1]);
            calc_dVrdr_2ndstep_and_more(vr[i],r[i],dr,rho[i],Br[i],Bphi[i],vphi[i],T[i],Lr[i],kappa[i],&vr[i+1],&denominator_of_dvrdr[i+1],&numerator_of_dvrdr[i+1]);
            solve_constraint_eqs(r[i+1],vr[i+1],T[i+1],rA,vA,Fm,FB,Lang,etot,&Br[i+1],&Bphi[i+1],&vphi[i+1],&Lr[i+1],&rho[i+1]);
            kappa[i+1] = kappa_fit(log10(T[i+1]),log10(rho[i+1]),kappa_tab);
            i++;
        }
        
        /* set next trial rA */
        if (i >= rbin-1) {
            break;
        } else if (denominator_of_dvrdr[i] < 0.){
            ln_dudxAmin = log(dudxA);
            dudxA = exp(.5*(log(dudxA)+ln_dudxAmax));
            ddudxA = dudxA - exp(ln_dudxAmin);
            printf("increase dudxA \n");
        } else if (numerator_of_dvrdr[i] < 0.) {
            ln_dudxAmax = log(dudxA);
            dudxA = exp(.5*(log(dudxA)+ln_dudxAmin));
            ddudxA = exp(ln_dudxAmin)-dudxA;
            printf("decrease dudxA \n");
        }
        
        j++;
        
    }
    
    if (j<max_1sttrial){
        printf("find a rough solution \n\n");
    } else {
        printf("could not find a rough solution \n\n");
    }
    
    FILE *op;
    op = fopen("rA2infty1.dat","w");
    for (j=0; j<rbin; j++) {
        //if (j % 100 == 0){
            fprintf(op,"%12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e \n",
                    r[j],vr[j],T[j],Br[j],Bphi[j],vphi[j],Lr[j],rho[j],kappa[j],denominator_of_dvrdr[j],numerator_of_dvrdr[j]);
        //}
    }
    fclose(op);
    
    *dudxA_1st = dudxA;
    *ddudxA_1st = ddudxA/(double)max_2ndtrial;
    
    free(r);
    free(vr);
    free(T);
    free(Br);
    free(Bphi);
    free(vphi);
    free(Lr);
    free(rho);
    free(denominator_of_dvrdr);
    free(numerator_of_dvrdr);
}


void set_para_at_rA(double rA, double dudxA, double TA, double LrA,
                    double *vA, double *BrA, double *rhoA, double *vphiA, double *BphiA, double *Fm, double *FB, double *Lang, double *etot)
{
    double vA_tmp = pow(Bwd,2.)*pow(Rwd,4.)/Mdot/pow(rA,2.0);
    double BrA_tmp = pow(rA/Rwd,-2.0)*Bwd;
    double rhoA_tmp = Mdot/4./M_PI/vA_tmp/pow(rA,2.0);
    double vphiA_tmp = rA*Omega*dudxA/(2.+ dudxA);
    double BphiA_tmp = -BrA_tmp*rA*Omega/vA_tmp*(2./(2.+ dudxA));
    
    double hA = 5./2.*kB*TA/mu_mol/Mu + 4.*arad*pow(TA,4.)/3./rhoA_tmp;
    double kA = 0.5*(vA_tmp*vA_tmp + vphiA_tmp*vphiA_tmp);
    double Fm_tmp = rhoA_tmp*vA_tmp*rA*rA;
    double Lang_tmp = rA*rA*Omega;
    double etot_tmp = LrA/4./M_PI/Fm_tmp + kA + hA - G*Mwd/rA - rA*Omega*vphiA_tmp + Lang_tmp*Omega;
   
    *vA = vA_tmp;
    *BrA = BrA_tmp;
    *rhoA = rhoA_tmp;
    *vphiA = vphiA_tmp;
    *BphiA = BphiA_tmp;
    
    *FB = Rwd*Rwd*Bwd;
    *Fm = Fm_tmp;
    *Lang = Lang_tmp;
    *etot = etot_tmp;
}


void set_r_from_rA_to_infty(double rA, double rmax, double r[])
{
    double del_ln_r = log(rmax/rA)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        r[i] = rA*exp(del_ln_r*(double)i);
    }
}


void calc_dTdr(double Tinput, double dr, double dTdr,
               double *Toutput)
{
    *Toutput = Tinput + dr*dTdr;
}


void calc_dVrdr_1ststep(double vrinput, double dr, double dudxA, double vA, double rA,
                        double *vroutput)
{
    *vroutput = vrinput + dr*dudxA*(vA/rA);
}


void calc_dVrdr_2ndstep_and_more(double vrinput, double r, double dr, double rho, double Br, double Bphi, double vphi, double T, double Lr, double kappa,
                                 double *vroutput, double *denominator_of_dvrdr, double *numerator_of_dvrdr)
{
    double Ar = Br/sqrt(4.*M_PI*rho);
    double Aphi = Bphi/sqrt(4.*M_PI*rho);
    
    double denominator_of_dvrdr_tmp = (vrinput*vrinput - kB*T/mu_mol/Mu - Aphi*Aphi*vrinput*vrinput/(vrinput*vrinput-Ar*Ar))*r/vrinput;
    double dPdr_term = 3.*kappa*Lr/16./M_PI/arad/C/pow(T,3.)/r*(rho*kB/mu_mol/Mu + 4./3.*arad*pow(T,3.)) + 2.*kB*T/mu_mol/Mu;
    double gravity_term = -G*Mwd/r;
    double centrifugal_force_term = vphi*vphi;
    double magnetic_term = 2.*vrinput*vphi*Ar*Aphi/(vrinput*vrinput-Ar*Ar);
    double numerator_of_dvrdr_tmp = dPdr_term + gravity_term + centrifugal_force_term + magnetic_term;
    
    *vroutput = vrinput + dr*numerator_of_dvrdr_tmp/denominator_of_dvrdr_tmp;
    *denominator_of_dvrdr = denominator_of_dvrdr_tmp;
    *numerator_of_dvrdr = numerator_of_dvrdr_tmp;
}


void solve_constraint_eqs(double r, double vr, double T,
                          double rA, double vA, double Fm, double FB, double Lang, double etot,
                          double *Br, double *Bphi, double *vphi, double *Lr, double *rho)
{
    double Br_tmp = FB/r/r;
    double rho_tmp = Fm/vr/r/r;
    double u = vr/vA;
    double x = r/rA;
    double vphi_tmp = rA*Omega*x*(1.-u)/(1.-x*x*u);
    double Bphi_tmp = -Br_tmp*rA*Omega/vA*x*(1.-x*x)/(1.-x*x*u);
    double a_tmp = sqrt(kB*T/mu_mol/Mu);

    double h_tmp = 5./2.*kB*T/mu_mol/Mu + 4.*arad*pow(T,4.)/3./rho_tmp;
    double k_tmp = 0.5*(vr*vr + vphi_tmp*vphi_tmp);
    double Lr_tmp = 4.*M_PI*Fm*(etot - k_tmp - h_tmp + G*Mwd/r + r*Omega*vphi_tmp - Lang*Omega);
    
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
    ip = fopen("tab46.dat","r");
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


double michel_wind_velocity(double r, double br, double omega, double mdot)
{
    return pow(pow(r,4.)*pow(br,2.)*pow(omega,2.)/mdot,1./3.);
}

