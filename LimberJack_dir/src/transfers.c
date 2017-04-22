#include "common.h"

//TODO: check non-limber rsd
//TODO: check spin  prefactors (sqrt?)
//TODO: check lensing signs

#define N_CHI 1000
static double transfer_dens_nolim(int l,double k,RunParams *par,int ibin)
{
  int i;
  double ret=0;
  double dchi=(par->chimax_nc[ibin]-par->chimin_nc[ibin])/N_CHI;
  for(i=0;i<N_CHI;i++) {
    double chi=par->chimin_nc[ibin]+(i+0.5)*dchi;
    double gf=spline_eval(chi,par->gfofchi);
    double h=spline_eval(chi,par->hofchi);
    double z=spline_eval(chi,par->zofchi);
    double pz=spline_eval(z,par->wind_0[ibin]);
    double bz=spline_eval(z,par->bias);
    double jl=csm_j_bessel(l,k*chi);

    ret+=h*pz*bz*gf*jl;
  }

  return ret*dchi;
}

static double transfer_rsd_nolim(int l,double k,RunParams *par,int ibin)
{
  int i;
  double ret=0;
  double dchi=(par->chimax_nc[ibin]-par->chimin_nc[ibin])/N_CHI;
  for(i=0;i<N_CHI;i++) {
    double chi=par->chimin_nc[ibin]+(i+0.5)*dchi;
    double x=k*chi;
    double gf=spline_eval(chi,par->gfofchi);
    double fg=spline_eval(chi,par->fgofchi);
    double h=spline_eval(chi,par->hofchi);
    double z=spline_eval(chi,par->zofchi);
    double pz=spline_eval(z,par->wind_0[ibin]);
    double jl=csm_j_bessel(l,x);
    double jlp1=csm_j_bessel(l+1,x);
    double jlpp=((x*x-l*(l-1))*jl-2*x*jlp1)/(x*x);

    ret+=h*pz*fg*gf*jlpp;
  }

  return ret*dchi;
}

static double transfer_magnification_nolim(int l,double k,RunParams *par,int ibin)
{
  int i;
  double ret=0;
  double dchi=par->chimax_nc[ibin]/N_CHI;
  for(i=0;i<N_CHI;i++) {
    double chi=(i+0.5)*dchi;
    double gf=spline_eval(chi,par->gfofchi);
    double a=spline_eval(chi,par->aofchi);
    double w=spline_eval(chi,par->wind_M[ibin]);
    double jl=csm_j_bessel(l,k*chi);
    
    ret+=gf*w*jl/(a*chi);
  }

  ret*=-2*par->prefac_lensing*l*(l+1)/(k*k);
  return ret*dchi;
}

static double transfer_lensing_nolim(int l,double k,RunParams *par,int ibin)
{
  int i;
  double ret=0;
  double dchi=par->chimax_nc[ibin]/N_CHI;
  for(i=0;i<N_CHI;i++) {
    double chi=(i+0.5)*dchi;
    double gf=spline_eval(chi,par->gfofchi);
    double a=spline_eval(chi,par->aofchi);
    double w=spline_eval(chi,par->wind_L[ibin]);
    double jl=csm_j_bessel(l,k*chi);
    
    ret+=gf*w*jl/(a*chi);
  }

  ret*=par->prefac_lensing*sqrt((l+2.)*(l+1.)*l*(l-1.))/(k*k);
  return ret*dchi;
}

static double transfer_IA_NLA_nolim(int l,double k,RunParams *par,int ibin)
{
  int i;
  double ret=0;
  double dchi=(par->chimax_nc[ibin]-par->chimin_nc[ibin])/N_CHI;
  for(i=0;i<N_CHI;i++) {
    double chi=par->chimin_nc[ibin]+(i+0.5)*dchi;
    double x=k*chi;
    double gf=spline_eval(chi,par->gfofchi);
    double h=spline_eval(chi,par->hofchi);
    double z=spline_eval(chi,par->zofchi);
    double pz=spline_eval(z,par->wind_0[ibin]);
    double az=spline_eval(z,par->abias);
    double jl=csm_j_bessel(l,x);
    
    ret+=h*pz*az*gf*jl/(x*x);
  }

  ret*=sqrt((l+2.)*(l+1.)*l*(l-1.));
  return ret*dchi;
}

static double transfer_cmblens_nolim(int l,double k,RunParams *par)
{
  int i;
  double ret=0;
  double dchi=par->chi_kappa/N_CHI;
  for(i=0;i<N_CHI;i++) {
    double chi=(i+0.5)*dchi;
    double w=1-chi/par->chi_kappa;
    double a=spline_eval(chi,par->aofchi);
    double gf=spline_eval(chi,par->gfofchi);
    double jl=csm_j_bessel(l,k*chi);

    ret+=gf*w*jl/(chi*a);
  }
  
  ret*=(par->prefac_lensing)*l*(l+1)/(k*k);
  return ret*dchi;
}

static double transfer_isw_nolim(int l,double k,RunParams *par)
{
  int i;
  double ret=0;
  double dchi=par->chi_isw/N_CHI;
  for(i=0;i<N_CHI;i++) {
    double chi=(i+0.5)*dchi;
    double gf=spline_eval(chi,par->gfofchi);
    double fg=spline_eval(chi,par->fgofchi);
    double h=spline_eval(chi,par->hofchi);
    double jl=csm_j_bessel(l,k*chi);

    ret+=h*gf*(1-fg)*jl;
  }
  
  ret*=2*par->prefac_lensing/(k*k);
  return ret*dchi;
}

static double transfer_dens(int l,double k,RunParams *par,int ibin)
{
  double chi=(l+0.5)/k;
  double gf=spline_eval(chi,par->gfofchi);
  double h=spline_eval(chi,par->hofchi);
  double z=spline_eval(chi,par->zofchi);
  double pz=spline_eval(z,par->wind_0[ibin]);
  double bz=spline_eval(z,par->bias);

  return pz*bz*gf*h;
}

static double transfer_rsd(int l,double k,RunParams *par,int ibin)
{
  double chi0=(l+0.5)/k;
  double chi1=(l+1.5)/k;
  double z0=spline_eval(chi0,par->zofchi);
  double z1=spline_eval(chi1,par->zofchi);
  double gf0=spline_eval(chi0,par->gfofchi);
  double gf1=spline_eval(chi1,par->gfofchi);
  double fg0=spline_eval(chi0,par->fgofchi);
  double fg1=spline_eval(chi1,par->fgofchi);
  double h0=spline_eval(chi0,par->hofchi);
  double h1=spline_eval(chi1,par->hofchi);
  double pz0=spline_eval(z0,par->wind_0[ibin]);
  double pz1=spline_eval(z1,par->wind_0[ibin]);
  double term0=pz0*fg0*gf0*h0*(1+8.*l)/((2*l+1.)*(2*l+1.));
  double term1=pz1*fg1*gf1*h1*sqrt((l+0.5)/(l+1.5))*4./(2*l+3);

  return term0-term1;
}

static double transfer_magnification(int l,double k,RunParams *par,int ibin)
{
  double chi=(l+0.5)/k;
  double gf=spline_eval(chi,par->gfofchi);
  double a=spline_eval(chi,par->aofchi);
  double w=spline_eval(chi,par->wind_M[ibin]);

  if(w<=0)
    return 0;
  else
    return -2*par->prefac_lensing*l*(l+1)*gf*w/(a*chi*k*k);
}

static double transfer_lensing(int l,double k,RunParams *par,int ibin)
{
  double chi=(l+0.5)/k;
  double gf=spline_eval(chi,par->gfofchi);
  double a=spline_eval(chi,par->aofchi);
  double w=spline_eval(chi,par->wind_L[ibin]);

  if(w<=0)
    return 0;
  else
    return par->prefac_lensing*sqrt((l+2.)*(l+1.)*l*(l-1.))*gf*w/(a*chi*k*k);
  //    return par->prefac_lensing*l*(l+1)*gf*w/(a*chi*k*k);
}

static double transfer_IA_NLA(int l,double k,RunParams *par,int ibin)
{
  double chi=(l+0.5)/k;
  double gf=spline_eval(chi,par->gfofchi);
  double h=spline_eval(chi,par->hofchi);
  double z=spline_eval(chi,par->zofchi);
  double pz=spline_eval(z,par->wind_0[ibin]);
  double az=spline_eval(z,par->abias);

  return pz*az*gf*h*sqrt((l+2.)*(l+1.)*l*(l-1.))/((l+0.5)*(l+0.5));
}

static double transfer_cmblens(int l,double k,RunParams *par)
{
  double chi=(l+0.5)/k;

  if(chi>=par->chi_kappa)
    return 0;
  else {
    double gf=spline_eval(chi,par->gfofchi);
    double a=spline_eval(chi,par->aofchi);
    double w=1-chi/par->chi_kappa;

    return par->prefac_lensing*l*(l+1)*gf*w/(a*chi*k*k);
  }
}

static double transfer_isw(int l,double k,RunParams *par)
{
  double chi=(l+0.5)/k;

  if(chi>=par->chi_isw)
    return 0;
  else {
    double gf=spline_eval(chi,par->gfofchi);
    double h=spline_eval(chi,par->hofchi);
    double fg=spline_eval(chi,par->fgofchi);

    return 2*par->prefac_lensing*h*gf*(1-fg)/(k*k);
  }
}

double transfer_wrap(int l,double k,RunParams *par,char *trtype,int ibin)
{
  double tr=-1;
  if(!strcmp(trtype,"nc")) {
    if(par->do_nc!=1)
      report_error(1,"Asked to calculate NC transfer function, but can't!\n");
    else {
      tr=0;
      if(par->has_dens) {
	if(l<=par->l_limber_min)
	  tr+=transfer_dens_nolim(l,k,par,ibin);
	else
	  tr+=transfer_dens(l,k,par,ibin);
      }
      if(par->has_rsd) {
	if(l<=par->l_limber_min)
	  tr+=transfer_rsd_nolim(l,k,par,ibin);
	else
	  tr+=transfer_rsd(l,k,par,ibin);
      }
      if(par->has_lensing) {
	if(l<=par->l_limber_min)
	  tr+=transfer_magnification_nolim(l,k,par,ibin);
	else
	  tr+=transfer_magnification(l,k,par,ibin);
      }
    }
  }
  else if(!strcmp(trtype,"shear")) {
    if(par->do_shear!=1)
      report_error(1,"Asked to calculate shear transfer function, but can't!\n");
    else {
      if(l<=par->l_limber_min)
	tr=transfer_lensing_nolim(l,k,par,ibin);
      else
	tr=transfer_lensing(l,k,par,ibin);
      if(par->has_intrinsic_alignment) {
	if(l<=par->l_limber_min)
	  tr+=transfer_IA_NLA_nolim(l,k,par,ibin);
	else
	  tr+=transfer_IA_NLA(l,k,par,ibin);
      }
    }
  }
  else if(!strcmp(trtype,"cmblens")) {
    if(par->do_cmblens!=1)
      report_error(1,"Asked to calculate shear transfer function, but can't!\n");
    else {
      if(l<=par->l_limber_min)
	tr=transfer_cmblens_nolim(l,k,par);
      else
	tr=transfer_cmblens(l,k,par);
    }
  }
  else if(!strcmp(trtype,"isw")) {
    if(par->do_isw!=1)
      report_error(1,"Asked to calculate ISW transfer function, but can't!\n");
    else {
      if(l<=par->l_limber_min)
	tr=transfer_isw_nolim(l,k,par);
      else
	tr=transfer_isw(l,k,par);
    }
  }
  else {
    report_error(1,"Unknown transfer type %s\n",trtype);
    tr=-1.;
  }

  if(l>par->l_limber_min)
    tr*=sqrt(M_PI/(2*l+1.))/k;

  return tr;
}
