#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void load_data();


int main()
{
    load_data();
    return 0;
}

void load_data()
{
    int i,rbin = 100000;
    
    double *r,*vr,*T,*Br,*Bphi,*vphi,*Lr,*rho,*kappa;
    r = (double *)malloc(2*rbin * sizeof(double));
    vr = (double *)malloc(2*rbin * sizeof(double));
    T = (double *)malloc(2*rbin * sizeof(double));
    Br = (double *)malloc(2*rbin * sizeof(double));
    Bphi = (double *)malloc(2*rbin * sizeof(double));
    vphi = (double *)malloc(2*rbin * sizeof(double));
    Lr = (double *)malloc(2*rbin * sizeof(double));
    rho = (double *)malloc(2*rbin * sizeof(double));
    kappa = (double *)malloc(2*rbin * sizeof(double));
    
    FILE *ip;
    ip = fopen("../inshot.dat","r");
    for (i = rbin-1; i==0; i--){
        fscanf(ip,"%*lf %*lf %*lf %*lf %*lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %*lf %*lf %*lf %*lf \n",
               &r[i],&vr[i],&T[i],&rho[i],&vphi[i],&Br[i],&Bphi[i],&Lr[i],&kappa[i]);
    }
    fclose(ip);
//    ip = fopen("../outshot.dat","r");
//    for (i = rbin; i<2*rbin-2; i++){
//        fscanf(ip,"%*lf %*lf %*lf %*lf %*lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %*lf %*lf %*lf %*lf \n",
//               &r[i],&vr[i],&T[i],&rho[i],&vphi[i],&Br[i],&Bphi[i],&Lr[i],&kappa[i]);
//    }
//    fclose(ip);
    
    printf("v_inf = %lf [cm/s] \n",vr[i]);
    
//    double tau=0.;
//    while (tau < 1.){
//        tau += rho[i]*kappa[i]*(r[i]-r[i-1]);
//        printf("%12.7e \n",tau);
//        i --;
//    }
//    printf("r_ph = %lf [cm] \n",r[i]);
    
    
    free(r);
    free(vr);
    free(T);
    free(Br);
    free(Bphi);
    free(vphi);
    free(Lr);
    free(rho);
    
    return ;
}
