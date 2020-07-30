#include <stdio.h>
#include <math.h>

const int index_T=71;
const int index_R=20;
const int index_table_tot = 71*20;

void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T6, double log10rho, double kappa_tab[index_T][index_R]);


int main()
{
  FILE* op;
  double kappa_tab[index_T][index_R],kappa;
  load_kappa_table(kappa_tab);
  
  double log10Tmin = 6.0,log10Tmax = 9.5;
  double log10rhomin = -7.,log10rhomax = 6.0;
  int i,j,bin=99;
  double log10T[bin+1],log10rho[bin+1];

  op = fopen("kappa_fit.dat","w");
  for (i=0; i<bin; i++){
    log10T[i] = log10Tmin + (log10Tmax-log10Tmin)*i/(double)(bin-1);
    log10rho[i] = log10rhomin + (log10rhomax-log10rhomin)*i/(double)(bin-1);
    printf("%12.3e %12.3e \n",log10T[i], log10rho[i]);
  } 

  for (i=0;i<bin;i++){
    fprintf(op,"%le ",log10rho[i]);
    for (j=0;j<bin;j++){
      fprintf(op,"%le ",kappa_fit(log10T[j],log10rho[i],kappa_tab));
    }
    fprintf(op,"\n");
  }
  fclose(op);

  return 0;
}


void load_kappa_table(double kappa_tab[index_T][index_R])
{
    int i,j;
    
    FILE *ip;
    ip = fopen("tab56.dat","r");
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
        index_T_fit = 2;
        log10T = kappa_tab[1][0];
    } else if (index_T_fit == index_T){
        /* extrapolating the table for higher temperatures */;
        index_T_fit --;
    }

    j=1;
    while(kappa_tab[0][j]<log10R && j<index_R)
        j++;
    index_R_fit = j;

    if (index_R_fit == 1){
        index_R_fit = 2;
        log10R = kappa_tab[0][1];
    } else if (index_R_fit == index_R){
        log10R = kappa_tab[0][index_R-1];
    }

    double kappa_Tdirect_min = (kappa_tab[index_T_fit-1][index_R_fit]-kappa_tab[index_T_fit-1][index_R_fit-1])/(kappa_tab[0][index_R_fit]-kappa_tab[0][index_R_fit-1])*(log10R-kappa_tab[0][index_R_fit-1])+kappa_tab[index_T_fit-1][index_R_fit-1];
    double kappa_Tdirect_max = (kappa_tab[index_T_fit][index_R_fit]-kappa_tab[index_T_fit][index_R_fit-1])/(kappa_tab[0][index_R_fit]-kappa_tab[0][index_R_fit-1])*(log10R-kappa_tab[0][index_R_fit-1])+kappa_tab[index_T_fit][index_R_fit-1];
    double log10kappa = (kappa_Tdirect_max-kappa_Tdirect_min)/(kappa_tab[index_T_fit][0]-kappa_tab[index_T_fit-1][0])*(log10T-kappa_tab[index_T_fit-1][0])+kappa_Tdirect_min;
    
    return pow(10.,log10kappa);

}
