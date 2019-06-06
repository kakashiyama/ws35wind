#ifndef TEST
#define TEST

void shooting();
void radial_step(double x, double y[], double yp[]);
void rk(double x, double dx, double y[], double yp[]);
void set_r_from_rA_to_infty(double rA, double r[]);
void solve_constraint_eqs(double r, double vr, double T, double etot, double *rho, double *vphi, double *Br, double *Bphi, double *Lr, double *kappa);
void calc_derivatives(double r, double vr, double T, double rho, double vphi, double Br, double Bphi, double Lr, double kappa, double *dvrdr, double *dTdr, double *detotdr);
void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T, double log10rho, double kappa_tab[index_T][index_R]);
double solve_Rfld(double r, double T, double Lr);
double calc_lambda(double Rfld);
double solve_dTdr(double rho, double kappa, double T, double Rfld);

#endif
