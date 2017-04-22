#ifndef _COMMON_
#define _COMMON_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include "params.h"
#include "cosmo_mad.h"

typedef struct {
  gsl_interp_accel *intacc;
  gsl_spline *spline;
  double x0,xf;
  double y0,yf;
} SplPar;

typedef struct {
  double om,ol,ob;
  double w0,wa,h0;
  double ns,s8;
  char **fname_window;
  char fname_bias[256];
  char fname_sbias[256];
  char fname_abias[256];
  char fname_pk[256];
  char prefix_out[256];
  int lmax;
  Csm_params *cpar;
  double chi_horizon;
  double z_kappa;
  double z_isw;
  double chi_kappa;
  double chi_isw;
  double chimax_nc[2];
  double chimin_nc[2];
  double prefac_lensing;
  double dchi;
  double r_smooth;
  int l_limber_min;
  int do_nc;
  int do_shear;
  int do_cmblens;
  int do_isw;
  int has_bg;
  int has_dens;
  int has_rsd;
  int has_lensing;
  int has_intrinsic_alignment;
  SplPar *aofchi;
  SplPar *zofchi;
  SplPar *hofchi;
  SplPar *gfofchi;
  SplPar *fgofchi;
  SplPar **wind_0;
  SplPar **wind_M;
  SplPar **wind_L;
  SplPar *bias;
  SplPar *sbias;
  SplPar *abias;
  double *cl_dd;
  double *cl_d1l2;
  double *cl_d2l1;
  double *cl_dc;
  double *cl_di;
  double *cl_ll;
  double *cl_lc;
  double *cl_li;
  double *cl_cc;
  double *cl_ci;
  double *cl_ii;
  int do_w_theta;
  int do_w_theta_logbin;
  double th_min;
  double th_max;
  int n_th;
  int n_th_logint;
  double *wt_dd;
  double *wt_d1l2;
  double *wt_d2l1;
  double *wt_dc;
  double *wt_di;
  double *wt_ll_pp;
  double *wt_ll_mm;
  double *wt_lc;
  double *wt_li;
  double *wt_cc;
  double *wt_ci;
  double *wt_ii;
} RunParams;

//Defined in common.c
int my_linecount(FILE *f);
void report_error(int level,char *fmt,...);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
FILE *my_fopen(const char *path,const char *mode);
SplPar *spline_init(int n,double *x,double *y,double y0,double yf);
double spline_eval(double x,SplPar *spl);
void spline_free(SplPar *spl);
RunParams *param_new(void);
void param_free(RunParams *par);

//Defined in cosmo.c
RunParams *init_params(char *fname_ini);

//Defined in transfers.c
double transfer_wrap(int l,double k,RunParams *par,char *trtype,int ibin);

//Defined in spectra.c
void compute_spectra(RunParams *par);
void compute_w_theta(RunParams *par);

//Defined in io.c
int read_parameter_file(char *fname,RunParams *par);
void write_output(RunParams *par);

#endif //_COMMON_
