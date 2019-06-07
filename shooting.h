struct _input load_inputfile();
struct _fixed calc_fixed_para(struct _input in, double rA, double dudxA);
int inshot(struct _input in, struct _fixed fix, double rA, double dudxA, double rstop, int rbin);
int outshot(struct _input in, struct _fixed fix, double rA, double dudxA, double rstop, int rbin);
void set_r_from_rA_to_infty(double rA, double rstop, int rbin, double r[]);
void radial_step(struct _input in, struct _fixed fix, double rA, double dudxA, double x, double y[], double yp[]);
void rk(struct _input in, struct _fixed fix, double rA, double dudxA, double x, double dx, double y[], double yp[]);
void solve_constraint_eqs(struct _input in, struct _fixed fix, double rA,
                          double r, double vr, double T, double etot,
                          double *rho, double *vphi, double *Br, double *Bphi, double *Lr, double *kappa);
void calc_derivatives(struct _input in,
                      double r, double vr, double T, double rho, double vphi, double Br, double Bphi, double Lr, double kappa,
                      double *dvrdr, double *dTdr, double *detotdr, double *nume, double *deno);
void load_kappa_table(double kappa_tab[index_T][index_R]);
double kappa_fit(double log10T, double log10rho, double kappa_tab[index_T][index_R]);
double solve_Rfld(double r, double T, double Lr);
double calc_lambda(double Rfld);
double solve_dTdr(double rho, double kappa, double T, double Rfld);

struct _input{
  double mu_mol;
  double Mwd;
  double Rwd;
  double Bwd;
  double Omega;
  double Mdot;
  double TA;
  double LrA;
};

struct _fixed{
  double Fm;
  double FB;
  double vA;
  double rhoA;
  double BrA;
  double vphiA;
  double BphiA;
  double Lang;
  double kA;
  double hA;
  double etotA;
};
