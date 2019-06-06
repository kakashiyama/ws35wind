#include <stdio.h>

const int n_para = 6;

void read_inputfile(struct _para_in);

struct _para_in{
    double mu_mol;
    double Mwd;
    double Rwd;
    double Bwd;
    double Omega;
    double Mdot;
};

int main()
{
    struct _para_in in_para;
    double input[n_para];
    read_inputfile(in_para);
    return 0;
}

void read_inputfile(struct _para_in)
{
    double input[n_para];
    int i=0;
    FILE *ip;
    ip = fopen("input.dat","r");
    while (fscanf(ip,"%*s %lf %*[\n]",&input[i])!=EOF){
        i++;
    }
    fclose(ip);
    
    struct _para_in in_para;
    in_para.mu_mol = input[0];
    in_para.Mwd = input[1];
    in_para.Rwd = input[2];
    in_para.Bwd = input[3];
    in_para.Omega = input[4];
    in_para.Mdot = input[5];
    
    printf("%12.3e \n",in_para.mu_mol);
    
    return ;
}
