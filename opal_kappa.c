#include <stdio.h>
#include <math.h>

/* see table #46 */
const int index_T=58;
const int index_R=20;

void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T6, double log10rho, double kappa_tab[index_T][index_R]);


int main()
{
    double kappa_tab[index_T][index_R],kappa;
    load_kappa_table(kappa_tab);
    double log10T=5.31,log10rho=-5.;
    kappa = kappa_fit(log10T,log10rho,kappa_tab);

    printf("%le \n",kappa);

    return 0;
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

