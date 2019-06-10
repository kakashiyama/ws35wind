#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "const.h"
#include "shooting.h"

/* number of radial bin */

int main()
{
    /* load input file */
    struct _input in = load_inputfile();
    struct _fixed fix;

    /* trial parameters for each shot */
    int flag;
    double rstop = in.Rwd;
    int rbin = 10000;
    int nshot = 1;
    int nshot_max = 50;
    //int nshot_max = 1;
    //double dudxA = 2.080621739e+00;
    double dudxA = 1.;
    
    double logrAmax = log(1.*Rsun);
    double logrAmin = log(2.*sqrt(in.LrA/(4.*M_PI*arad*pow(in.TA,4.)*C)));
    double rA = exp(0.5*(logrAmax+logrAmin));
    //double rA = 5.440052523e+09;
    
    flag = 99;
    while (flag != 0 && nshot < nshot_max){
        fix = calc_fixed_para(in,rA,dudxA);
        flag = inshot(in,fix,rA,dudxA,rstop,rbin);

        printf("In-shot No. %d with rA = %12.9e [cm] --> end with << flag %d >> \n",nshot,rA,flag);

        if (flag == 1) {
            logrAmax = log(rA);
            rA = exp(0.5*(logrAmin+log(rA)));
        } else if (flag == 2){
            logrAmax = log(rA);
            rA = exp(0.5*(logrAmin+log(rA)));
        } else if (flag == 3){
            logrAmin = log(rA);
            rA = exp(0.5*(log(rA)+logrAmax));
        } else {
            break;
        }
        nshot++;
    }
    
    rstop = 100.*rA;
    flag = outshot(in,fix,rA,dudxA,rstop,rbin);
    printf("Out-shot with rA = %12.9e [cm] --> end with << flag %d >> \n",rA,flag);
    
    return 0;
}

struct _input load_inputfile()
{
    double dmy[8];
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

int inshot(struct _input in, struct _fixed fix, double rA, double dudxA, double rstop, int rbin)
{
    int flag = 0;
    
    int i;
    double x,dx;
    double x_start;
    double x_stop;
    double y[3];
    double yp[3];
    
    double r[rbin],vr,T,etot;
    double rho,vphi,Br,Bphi,Lr,kappa;
    double dvrdr,dTdr,detotdr,nume,deno;
    set_r_from_rA_to_infty(rA,rstop,rbin,r);
    
    y[0] = 1.;
    y[1] = 1.;
    y[2] = 1.;
    
    FILE *op;
    op = fopen("inshot.dat","w");
    for (i=0; i<rbin-1; i++) {
        x = r[i]/rA;
        dx = (r[i+1]-r[i])/rA;
        
        /* only for output */
        vr = y[0]*fix.vA;
        T = y[1]*in.TA;
        etot = y[2]*fix.etotA;
        solve_constraint_eqs(in,fix,rA,r[i],vr,T,etot,&rho,&vphi,&Br,&Bphi,&Lr,&kappa);
        calc_derivatives(in,fix,rA,dudxA,r[i],vr,T,rho,vphi,Br,Bphi,Lr,kappa,&dvrdr,&dTdr,&detotdr,&nume,&deno);
        fprintf(op,"%12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e \n",
                x,y[0],y[1],yp[0],yp[1],r[i],vr,T,rho,vphi,Br,Bphi,Lr,kappa,etot,dvrdr,dTdr,detotdr);
        
        rk(in,fix,rA,dudxA,x,dx,y,yp);
        
        
        if ( y[0] < 0.){
            flag = 1;
            break;
        } else if ( isnan(y[0])) {
            flag = 2;
            break;
        }
        
    }
    fclose(op);

    if (i>=rbin-1 && nume*deno < 0.)
        flag = 3;
    
    return flag;
}


int outshot(struct _input in, struct _fixed fix, double rA, double dudxA, double rstop, int rbin)
{
    int flag = 0;
    
    int i;
    double x,dx;
    double x_start;
    double x_stop;
    double y[3];
    double yp[3];
    
    double r[rbin],vr,T,etot;
    double rho,vphi,Br,Bphi,Lr,kappa;
    double dvrdr,dTdr,detotdr,nume,deno;
    set_r_from_rA_to_infty(rA,rstop,rbin,r);
    
    y[0] = 1.;
    y[1] = 1.;
    y[2] = 1.;
    
    FILE *op;
    op = fopen("outshot.dat","w");
    for (i=0; i<rbin-1; i++) {
        x = r[i]/rA;
        dx = (r[i+1]-r[i])/rA;
        
        /* only for output */
        vr = y[0]*fix.vA;
        T = y[1]*in.TA;
        etot = y[2]*fix.etotA;
        solve_constraint_eqs(in,fix,rA,r[i],vr,T,etot,&rho,&vphi,&Br,&Bphi,&Lr,&kappa);
        calc_derivatives(in,fix,rA,dudxA,r[i],vr,T,rho,vphi,Br,Bphi,Lr,kappa,&dvrdr,&dTdr,&detotdr,&nume,&deno);
        fprintf(op,"%12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e \n",
                x,y[0],y[1],yp[0],yp[1],r[i],vr,T,rho,vphi,Br,Bphi,Lr,kappa,etot,dvrdr,dTdr,detotdr);
        
        rk(in,fix,rA,dudxA,x,dx,y,yp);
        
        if (nume*deno < 0.){
            if (nume < 0.)
                flag = 1;
            else
                flag = 2;
            break;
        }
        
    }
    fclose(op);
    
    return flag;
}


void set_r_from_rA_to_infty(double rA, double rstop, int rbin, double r[])
{
    double del_ln_r = log(rstop/rA)/(double)(rbin-1);
    int i;
    for (i=0; i<rbin; i++) {
        r[i] = rA*exp(del_ln_r*(double)i);
    }
}


void radial_step(struct _input in, struct _fixed fix, double rA, double dudxA, double x, double y[], double yp[])
{
    double u = y[0];
    double t = y[1];
    double e = y[2];
    
    double r = x*rA;
    double vr = u*fix.vA;
    double T = t*in.TA;
    double etot = e*fix.etotA;
    
    double rho,vphi,Br,Bphi,Lr,kappa;
    solve_constraint_eqs(in,fix,rA,r,vr,T,etot,&rho,&vphi,&Br,&Bphi,&Lr,&kappa);
    
    double dvrdr,dTdr,detotdr,nume,deno;
    calc_derivatives(in,fix,rA,dudxA,r,vr,T,rho,vphi,Br,Bphi,Lr,kappa,&dvrdr,&dTdr,&detotdr,&nume,&deno);
    
    if (x == 1.){
        yp[0] = dudxA;
    } else {
        yp[0] = dvrdr*(rA/fix.vA);
    }
    yp[1] = dTdr*(rA/in.TA);
    yp[2] = detotdr*(rA/fix.etotA);
    
    return;
}


void rk(struct _input in, struct _fixed fix, double rA, double dudxA, double x, double dx, double y[], double yp[])
{
    double y_tmp[3];
    double yp_tmp[3];
    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    
    y_tmp[0] = y[0];
    y_tmp[1] = y[1];
    y_tmp[2] = y[2];
    radial_step(in,fix,rA,dudxA,x,y,yp_tmp);
    k1[0] = dx*yp_tmp[0];
    k1[1] = dx*yp_tmp[1];
    k1[2] = dx*yp_tmp[2];
    
    yp[0] = yp_tmp[0];
    yp[1] = yp_tmp[1];
    yp[2] = yp_tmp[2];
    
    y_tmp[0] = y[0] + .5*k1[0];
    y_tmp[1] = y[1] + .5*k1[1];
    y_tmp[2] = y[2] + .5*k1[2];
    radial_step(in,fix,rA,dudxA,x+.5*dx,y_tmp,yp_tmp);
    k2[0] = dx*yp_tmp[0];
    k2[1] = dx*yp_tmp[1];
    k2[2] = dx*yp_tmp[2];
    
    y_tmp[0] = y[0] + .5*k2[0];
    y_tmp[1] = y[1] + .5*k2[1];
    y_tmp[2] = y[2] + .5*k2[2];
    radial_step(in,fix,rA,dudxA,x+.5*dx,y_tmp,yp_tmp);
    k3[0] = dx*yp_tmp[0];
    k3[1] = dx*yp_tmp[1];
    k3[2] = dx*yp_tmp[2];
    
    y_tmp[0] = y[0] + k3[0];
    y_tmp[1] = y[1] + k3[1];
    y_tmp[2] = y[2] + k3[2];
    radial_step(in,fix,rA,dudxA,x+dx,y_tmp,yp_tmp);
    k4[0] = dx*yp_tmp[0];
    k4[1] = dx*yp_tmp[1];
    k4[2] = dx*yp_tmp[2];
    
    y[0] = y[0] + (k1[0]+2.*k2[0]+2.*k3[0]+k4[0])/6.;
    y[1] = y[1] + (k1[1]+2.*k2[1]+2.*k3[1]+k4[1])/6.;
    y[2] = y[2] + (k1[2]+2.*k2[2]+2.*k3[2]+k4[2])/6.;
    
    return;
}


void solve_constraint_eqs(struct _input in, struct _fixed fix, double rA,
                          double r, double vr, double T, double etot,
                          double *rho, double *vphi, double *Br, double *Bphi, double *Lr, double *kappa)
{
    double x = r/rA;
    double u = vr/fix.vA;
    double t = T/in.TA;
    
    double rho_tmp = fix.rhoA/u/x/x;
    double Br_tmp = fix.BrA/x/x;
    double vphi_tmp,Bphi_tmp;
    if (x == 1.){
        vphi_tmp = fix.vphiA;
        Bphi_tmp = fix.BphiA;
    } else {
        vphi_tmp = rA*in.Omega*x*(1.-u)/(1.-x*x*u);
        Bphi_tmp = -Br_tmp*rA*in.Omega/fix.vA*x*(1.-x*x)/(1.-x*x*u);
    }
    
    double h = 5./2.*kB*T/in.mu_mol/Mu;
    double k = 0.5*(vr*vr + vphi_tmp*vphi_tmp);
    double Lr_tmp = 4.*M_PI*fix.Fm*(etot - k - h + G*in.Mwd/r + r*in.Omega*vphi_tmp - fix.Lang*in.Omega);
    
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

void calc_derivatives(struct _input in, struct _fixed fix, double rA, double dudxA,
                      double r, double vr, double T, double rho, double vphi, double Br, double Bphi, double Lr, double kappa,
                      double *dvrdr, double *dTdr, double *detotdr, double *nume, double *deno)
{
    double Rfld = solve_Rfld(r,T,Lr);
    double lambda = calc_lambda(Rfld);
    double dTdr_tmp = solve_dTdr(rho,kappa,T,Rfld);
    
    double Ar = Br/sqrt(4.*M_PI*rho);
    double Aphi = Bphi/sqrt(4.*M_PI*rho);
    double denominator_of_dvrdr = (vr*vr - kB*T/in.mu_mol/Mu - Aphi*Aphi*vr*vr/(vr*vr-Ar*Ar))*r/vr;
    
    //double dPdr_term = dTdr_tmp*(rho*kB/in.mu_mol/Mu + 4.*lambda*arad*pow(T,3.)) + 2.*kB*T/in.mu_mol/Mu;
    //double dPdr_term = rho*kappa*Lr/(4.*M_PI*r*r*C) + dTdr_tmp*rho*kB/in.mu_mol/Mu + 2.*kB*T/in.mu_mol/Mu;
    double dPdr_term = kappa*Lr/(4.*M_PI*r*C) + dTdr_tmp*r*kB/in.mu_mol/Mu + 2.*kB*T/in.mu_mol/Mu;
    double gravity_term = -G*in.Mwd/r;
    double centrifugal_force_term = vphi*vphi;
    double magnetic_term = 2.*vr*vphi*Ar*Aphi/(vr*vr-Ar*Ar);
    double numerator_of_dvrdr = dPdr_term + gravity_term + centrifugal_force_term + magnetic_term;
    
    
    if (vr != Ar)
        *dvrdr = numerator_of_dvrdr/denominator_of_dvrdr;
    else
        *dvrdr = dudxA*fix.vA/rA;
    *dTdr = dTdr_tmp;
    *detotdr = kappa*Lr/(4.*M_PI*r*r*C);//4.*lambda*arad*pow(T,3.)/rho*dTdr_tmp;
    *nume = numerator_of_dvrdr;
    *deno = denominator_of_dvrdr;
    
    return ;
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
    
    //return .2;
    return pow(10.,log10kappa);

}


double solve_Rfld(double r, double T, double Lr)
{
    /* analytic solution of Rfld for given r, T, and Lr */
    double s = Lr/(4.*M_PI*pow(r,2.)*arad*pow(T,4.)*C);
    //printf("%12.3e ",s);
    return (3.*s-2+sqrt(4.+12.*s-15.*s*s))/2./(1.-s);
}


double calc_lambda(double Rfld)
{
    //printf("%12.3e ",Rfld);
    return (2.+Rfld)/(6.+3.*Rfld+Rfld*Rfld);
}


double solve_dTdr(double rho, double kappa, double T, double Rfld)
{
    return -Rfld*kappa*rho*T/4.;
}
