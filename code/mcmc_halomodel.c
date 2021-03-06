/*  Copyright (C) <2016-2019>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */
 
/** @file   mcmc_halomodel.c
 *  @date   2013-2019
 *  @author Marcel van Daalen
 *  @author Stefan Hilbert
 *
 *  @brief  a halo model for clustering 
 *          in the MCMC parameter estimation
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>

#include "allvars.h"
#include "proto.h"
#include "mcmc_vars.h"
#include "mcmc_halomodel.h"
#include "mcmc_mpfit.h"

#define WORKSIZE 100000

const double xmin    = 0.0        ;
#ifdef MCRIT
const double xmax    = 6.         ;
#else
const double xmax    = 5.         ;
#endif

const double pow_ns  = 1.00       ;
const double kmax    = 1e4        ;

#ifdef PROJLIMITS
const double pimin   =  0.0       ;
const double pimax   = 40.0       ;
#endif

const double pa_low  = -2.0       ;
const double pa_high =  0.47712125;//log10(3)
const double pb_low  = -1.0       ;
const double pb_high =  1.0       ;
const double pc_low  = -1.0       ;
const double pc_high =  1.5       ;

void halomodel(const int snap, const float masslimit_low, const float masslimit_high, const int NR, double* r_arr, double* proj_arr) 
{
  const int NK=60;
  const int extralin=16;
  const double kbin=6./79.;
  const double kstart=-2;
  
  int    i;
  double k,p1h,p2h,m10,m11;
  double r,corrtmp;
  double *rCorrTable,*CorrTable;
  double Rhalo_cc,Rhalo_cs,Rhalo_ss,pcorr_cc,pcorr_cs,pcorr_ss;
  double *P2HTable_cc;
  double *P2HTable_cs;
  double *P2HTable_ss;
  double *kPcorrTable;
  double *PcorrTable;
  
  gsl_set_error_handler_off();
  setbuf(stdout,NULL);
  
  cutoff_low=malloc(6*sizeof(double));
  cutoff_high=malloc(6*sizeof(double));
  
  init_numgal(masslimit_low, masslimit_high, snap);
  ngal_mean=ngal_mean_calc(0);
  PowerTable=malloc((NK+extralin)*sizeof(double));
  kPowerTable=malloc((NK+extralin)*sizeof(double));
  m10=Mcensat(1,1,0);
  P2HTable_cc=malloc((NK+extralin)*sizeof(double));
  P2HTable_cs=malloc((NK+extralin)*sizeof(double));
  P2HTable_ss=malloc((NK+extralin)*sizeof(double));
  Rhalo_cc=2*calc_mean_rhalo_simple(0,0)/calc_mean_rhalo_simple(0,1);
  Rhalo_ss=2*calc_mean_rhalo_simple(1,0)/calc_mean_rhalo_simple(1,1);
  Rhalo_cs=0.5*(Rhalo_cc+Rhalo_ss);
  for (i=0; i<NK+extralin; i++)
  {
    if (i>2*extralin) k=pow(10.,(i-extralin)*kbin+kstart);
    else k=pow(10.,i*0.5*kbin+kstart);
    if (i>0 && kPowerTable[i-1]>2) k=pow(10.,kPowerTable[i-1]+2.044/6.);
    if (k>=0.01) p1h=2*Mcensat(k,0,1)+Mcensat(k,0,2);
    else p1h=0.0;
    m11=Mcensat(k,1,1);
    p2h=PowerSpec(k)*(m10*m10+2*m10*m11+m11*m11);
    P2HTable_cc[i]=log10(PowerSpec(k)*m10*m10);
    P2HTable_cs[i]=log10(PowerSpec(k)*2*m10*m11);
    P2HTable_ss[i]=log10(PowerSpec(k)*m11*m11);
    kPowerTable[i]=log10(k);
    PowerTable[i]=log10((p1h+p2h)/(gsl_spline_eval(ellipSpline,kPowerTable[i],ellipAcc)+1)/(gsl_spline_eval(alignSpline,kPowerTable[i],alignAcc)+1));
  } //for
  Twopow_ccAcc=gsl_interp_accel_alloc();
  Twopow_ccSpline=gsl_spline_alloc(gsl_interp_cspline,(NK+extralin));
  gsl_spline_init(Twopow_ccSpline,kPowerTable,P2HTable_cc,(NK+extralin));
  Twopow_csAcc=gsl_interp_accel_alloc();
  Twopow_csSpline=gsl_spline_alloc(gsl_interp_cspline,(NK+extralin));
  gsl_spline_init(Twopow_csSpline,kPowerTable,P2HTable_cs,(NK+extralin));
  Twopow_ssAcc=gsl_interp_accel_alloc();
  Twopow_ssSpline=gsl_spline_alloc(gsl_interp_cspline,(NK+extralin));
  gsl_spline_init(Twopow_ssSpline,kPowerTable,P2HTable_ss,(NK+extralin));
  PcorrTable=malloc((NK+extralin)*sizeof(double));
  kPcorrTable=malloc((NK+extralin)*sizeof(double));
  for (i=0; i<NK+extralin; ++i)
  {
    if (i>2*extralin) k=pow(10.,(i-extralin)*kbin+kstart);
    else k=pow(10.,i*0.5*kbin+kstart);
    if (i>0 && kPowerTable[i-1]>2) k=pow(10.,kPowerTable[i-1]+2.044/6.);
    kPcorrTable[i]=log10(k);
    if (k<4*M_PI/Rhalo_cc) 
    {
      pcorr_cc=min(4*M_PI/3.*pow(Rhalo_cc,3)*(pconv_W(k,Rhalo_cc,0)+TopHatWindow(k*Rhalo_cc)),TwoPowerSpec(k,0));
      if (Rhalo_cs>0)
        pcorr_cs=min((mugal_qawo(k,4*M_PI/3.*Delta*rho_mean*pow(Rhalo_cs,3))-TopHatWindow(k*Rhalo_cs/Delta_invth))*4*M_PI/3.*pow(Rhalo_cs,3)*(pconv_W(k,Rhalo_cs,1)+TopHatWindow(k*Rhalo_cs)),TwoPowerSpec(k,1));
      else
        pcorr_cs=0;
      if (Rhalo_ss>0)
        pcorr_ss=min((pow(mugal_qawo(k,4*M_PI/3.*Delta*rho_mean*pow(Rhalo_ss,3)),2)-pow(TopHatWindow(k*Rhalo_ss/Delta_invth),2))*4*M_PI/3.*pow(Rhalo_ss,3)*(pconv_W(k,Rhalo_ss,2)+TopHatWindow(k*Rhalo_ss)),TwoPowerSpec(k,2));
      else 
        pcorr_ss=0;
    } //if
    else
    {
      pcorr_cc=TwoPowerSpec(k,0);
      if (Rhalo_cs>0) pcorr_cs=TwoPowerSpec(k,1);
      else pcorr_cs=0;
      if (Rhalo_ss>0) pcorr_ss=TwoPowerSpec(k,2);
      else pcorr_ss=0;
    } //else
    PcorrTable[i]=pcorr_cc+pcorr_cs+pcorr_ss;
  } //for
  TwopcorrAcc=gsl_interp_accel_alloc();
  TwopcorrSpline=gsl_spline_alloc(gsl_interp_cspline,(NK+extralin));
  gsl_spline_init(TwopcorrSpline,kPcorrTable,PcorrTable,(NK+extralin));
  gsl_spline_free(pcSpline);
  gsl_interp_accel_free(pcAcc);
  gsl_spline_free(pbSpline);
  gsl_interp_accel_free(pbAcc);
  gsl_spline_free(paSpline);
  gsl_interp_accel_free(paAcc);
  for (i=5; i>=0; --i)
  {
    gsl_spline_free(NgalSpline[i]);
    gsl_interp_accel_free(NgalAcc[i]);
  } //for
  free(NgalSpline);
  free(NgalAcc);
  NewpowAcc=gsl_interp_accel_alloc();
  NewpowSpline=gsl_spline_alloc(gsl_interp_cspline,(NK+extralin));
  gsl_spline_init(NewpowSpline,kPowerTable,PowerTable,(NK+extralin));
#if defined(OUTPUTCORR) || defined(OUTPUTPOW)
  FILE *file_;
  char file_name_[512];
  const float mingalmass = 8.77 + (ThisTask%6  ) * 0.5;
  const float maxgalmass = 8.77 + (ThisTask%6+1) * 0.5;
#endif
#ifdef OUTPUTPOW
  sprintf(file_name_,"pow_%.2f-%.2f_%d.dat",mingalmass,maxgalmass,snap);
  file_=fopen(file_name_,"w");
  for (i=0; i<1000; ++i)
  {
    fprintf(file_,"%g %g\n",pow(10.,0.001*i*5-2),pow(10.,gsl_spline_eval(NewpowSpline,0.001*i*5-2,NewpowAcc))-gsl_spline_eval(TwopcorrSpline,0.001*i*5-2,TwopcorrAcc));
  } //for
  fclose(file_);
#endif
  CorrTable =malloc(NR*10*sizeof(double));
  rCorrTable=malloc(NR*10*sizeof(double));
  for (i=0; i<NR*10; ++i)
  {
    rCorrTable[i]=pow(10.,(i+0.5)*(log10(610.)+3.1)/(float)(NR*10)-3.1);
    CorrTable[i]=corr_qawo(rCorrTable[i],0.01,0.99)+corr_qawo(rCorrTable[i],1.,99.)+corr_qawo(rCorrTable[i],100.,900.)+corr_qawo(rCorrTable[i],1e3,5e3)+corr_qawo(rCorrTable[i],6e3,5e3);
  } //for
  CorrAcc=gsl_interp_accel_alloc();
  CorrSpline=gsl_spline_alloc(gsl_interp_cspline,NR*10);
  gsl_spline_init(CorrSpline,rCorrTable,CorrTable,NR*10);
#ifdef OUTPUTCORR
  sprintf(file_name_,"corr_%.2f-%.2f_%d.dat",mingalmass,maxgalmass,snap);
  file_=fopen(file_name_,"w");
  for (i=0; i<1000; ++i)
  {
    fprintf(file_,"%g %g\n",pow(10.,0.001*i*5.785-3),gsl_spline_eval(CorrSpline,pow(10.,0.001*i*5.785-3),CorrAcc));
  } //for
  fclose(file_);
#endif
  for (i=0; i<NR; ++i)
  {
    r=pow(10.,(i+0.5)*(log10(80.)+2.1)/(float)(NR-1)-2.3);
    corrtmp=proj_corr(r);
    r_arr[i]=r * inv_Hubble_h;
    proj_arr[i]=corrtmp * inv_Hubble_h;
  } //for
#ifdef OUTPUTPROJ
  sprintf(file_name_,"proj_%.2f-%.2f_%d.dat",mingalmass,maxgalmass,snap);
  file_=fopen(file_name_,"w");
  for (i=0; i<NR; ++i)
  {
    fprintf(file_,"%g   %g\n",r_arr[i],proj_arr[i]);
  } //for
  fclose(file_);
#endif
  gsl_spline_free(CorrSpline);
  gsl_interp_accel_free(CorrAcc);
  gsl_spline_free(NewpowSpline);
  gsl_interp_accel_free(NewpowAcc);
  gsl_spline_free(TwopcorrSpline);
  gsl_interp_accel_free(TwopcorrAcc);
  gsl_spline_free(Twopow_ssSpline);
  gsl_interp_accel_free(Twopow_ssAcc);
  gsl_spline_free(Twopow_csSpline);
  gsl_interp_accel_free(Twopow_csAcc);
  gsl_spline_free(Twopow_ccSpline);
  gsl_interp_accel_free(Twopow_ccAcc);
  free(kPcorrTable);
  free(PcorrTable);
  free(P2HTable_ss);
  free(P2HTable_cs);
  free(P2HTable_cc);
  free(kPowerTable);
  free(PowerTable);
  free(cutoff_high);
  free(cutoff_low);
} //halomodel


double TwoPowerSpec(const double k, const int censat)
{
  if (censat==0)
  {
    if (k>0.01) 
      return pow(10.,gsl_spline_eval(Twopow_ccSpline,log10(k),Twopow_ccAcc));
    else if (k>1e-4) 
      return PowerSpec(k)/PowerSpec(0.01)*pow(10.,gsl_spline_eval(Twopow_ccSpline,-2.,Twopow_ccAcc));
    else 
      return pow(k,pow_ns)/pow(1e-4,pow_ns)*PowerSpec(1e-4)/PowerSpec(0.01)*pow(10.,gsl_spline_eval(Twopow_ccSpline,-2.,Twopow_ccAcc));
  } //if
  else if (censat==1)
  {
    if (k>0.01) 
      return pow(10.,gsl_spline_eval(Twopow_csSpline,log10(k),Twopow_csAcc));
    else if (k>1e-4)
      return PowerSpec(k)/PowerSpec(0.01)*pow(10.,gsl_spline_eval(Twopow_csSpline,-2.,Twopow_csAcc));
    else 
      return pow(k,pow_ns)/pow(1e-4,pow_ns)*PowerSpec(1e-4)/PowerSpec(0.01)*pow(10.,gsl_spline_eval(Twopow_csSpline,-2.,Twopow_csAcc));
  } //else if
  else /* if (censat == 2) */
  {
    if (k>0.01)
      return pow(10.,gsl_spline_eval(Twopow_ssSpline,log10(k),Twopow_ssAcc));
    else if (k>1e-4)
      return PowerSpec(k)/PowerSpec(0.01)*pow(10.,gsl_spline_eval(Twopow_ssSpline,-2.,Twopow_ssAcc));
    else
      return pow(k,pow_ns)/pow(1e-4,pow_ns)*PowerSpec(1e-4)/PowerSpec(0.01)*pow(10.,gsl_spline_eval(Twopow_ssSpline,-2.,Twopow_ssAcc));
  } //else
} //TwoPowerSpec


double pconv_W_P_func(const double theta_, void *p)
{
  const struct conv_W_P_params *params=(struct conv_W_P_params *)p;
  const double k=(params->k);
  const double q=(params->q);
  const int censat=(params->censat);
  return TwoPowerSpec(sqrt(k*k-2*q*k*cos(theta_)+q*q),censat)*sin(theta_);
} //pconv_W_P_func


double pconv_W_func(const double lq, void *p)
{
  const struct conv_W_params *params=(struct conv_W_params *)p;
  const double k=(params->k);
  const double R=(params->R);
  const int censat=(params->censat);
  const double q=exp(lq);
  const double arg=(k*k+q*q-kmax*kmax)/(2.*k*q);
  
  double result=0,abserr,thetamax;
  
  if (arg>-1 && arg<1) thetamax=acos(arg);
  else if (arg<=-1) thetamax=M_PI;
  else thetamax=0.;
  
  gsl_function F;
  // int status;
  struct conv_W_P_params params2={ k,q,censat };
  F.function=&pconv_W_P_func;
  F.params=&params2;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  gsl_integration_qag(&F,0.,thetamax,0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return TopHatWindow(q*R)*result*q*q*q;
} //pconv_W_func


double pconv_W(const double k, const double R, const int censat)
{
  double result=0,abserr;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  struct conv_W_params params={ k,R,censat };
  F.function=&pconv_W_func;
  F.params=&params;
  gsl_integration_qag(&F,log(pow(10.,-8.)),log(k+kmax),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result/(4*M_PI*M_PI);
} //pconv_W


double calc_mean_rhalo_simple_func(const double lm, void *p)
{
  const struct rhalo_params *params=(struct rhalo_params *)p;
  const int censat=(params->censat);
  const int norm=(params->norm);
  const double m=exp(lm);
  
  double rterm=1.;
  if (norm==0) 
  {
    rterm=Delta_invth*Radius(m);
    if (censat>0)
      rterm-=min(calc_mean_rsat(log10(m),1)/calc_mean_rsat(log10(m),0),1.)*Delta_invth*Radius(m);
  } //if
  return nbargal(m)*NgalF(m,2+censat)*m*rterm;
} //calc_mean_rhalo_simple_func


double calc_mean_rhalo_simple(const int censat, const int norm) 
{
  double result=0,abserr;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  struct rhalo_params params={ censat,norm };
  F.function=&calc_mean_rhalo_simple_func;
  F.params=&params;
  gsl_integration_qag(&F,log(pow(10.,cutoff_low[2+censat])),log(pow(10.,cutoff_high[2+censat])),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //calc_mean_rhalo_simple


double calc_mean_rsat_func(const double lr, void *p) 
{
  const struct rsat_params *params=(struct rsat_params *)p;
  const double lm=(params->lm);
  const double pa=pow(10.,pa_eval(lm));
  const double pb=pow(10.,pa_eval(lm));
  const double pc=pow(10.,pa_eval(lm));
  const int extrar=(params->extrar);
  const double r=pow(10.,lr);
  
  // double nm= pow(r/pb,pa)/r * exp(-pow(r/pb, pc));
  // if (extrar)
  //   nm*=r;
  // return r*nm;
  
  const double nmr = pow(r/pb,pa) * exp(-pow(r/pb, pc));
  if (extrar)
    return nmr*r;
  else
    return nmr;
} //calc_mean_rsat_func


double calc_mean_rsat(const double lm, const int extrar)
{
  struct rsat_params params={ lm,extrar };
  double result=0,abserr;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&calc_mean_rsat_func;
  F.params=&params;
  gsl_integration_qag(&F,-4.,log10(xmax),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS51,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //calc_mean_rsat


double NewPowerSpec(const double k)
{
  //return pow(10.,gsl_spline_eval(NewpowSpline,log10(k),NewpowAcc));
  return pow(10.,gsl_spline_eval(NewpowSpline,log10(k),NewpowAcc))-gsl_spline_eval(TwopcorrSpline,log10(k),TwopcorrAcc);
} //NewPowerSpec


double corr_qawo_func(const double k,void *params)
{
  return k*NewPowerSpec(k)*exp(-k*k/1e7); //exponential needed for convergence
  (void)params; /* suppress unused-parameter warning */
} //corr_qawo_func


double corr_qawo(const double r, const double a, const double L)
{
  double result=0,abserr;
  gsl_function F;
  //int status;
  F.function=&corr_qawo_func;
  F.params=0;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  gsl_integration_qawo_table *t=gsl_integration_qawo_table_alloc(r,L,GSL_INTEG_SINE,40);
  /* status = */ gsl_integration_qawo(&F,a,0,1.0e-3,WORKSIZE,w,t,&result,&abserr);
  gsl_integration_qawo_table_free(t);
  gsl_integration_workspace_free(w);
  return result/(2.*M_PI*M_PI*r);
} //corr_qawo


double proj_corr_func(const double r, void *params) 
{
  const double sigma= *(double *) params;
  return 2*r*gsl_spline_eval(CorrSpline,r,CorrAcc)/sqrt(pow(r,2)-pow(sigma,2));
} //proj_corr_func


double proj_corr(double sigma)
{
  double result=0,abserr;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&proj_corr_func;
  F.params=&sigma;
#ifndef PROJLIMITS
  gsl_integration_qag(&F,sigma,600.,0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS51,w,&result,&abserr);
#else
  gsl_integration_qag(&F,sqrt(pimin*pimin+sigma*sigma),sqrt(pimax*pimax+sigma*sigma),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS51,w,&result,&abserr);
#endif
  gsl_integration_workspace_free(w);
  return result;
} //proj_corr


double Radius(const double m)
{
  return pow(m/rho_mean*3./4./M_PI,1./3.);
} //Radius


double Sigma2(const double m) 
{
  return pow(10.,gsl_spline_eval(SigmaSpline,log10(m),SigmaAcc));
} //Sigma2


double PowerSpec(const double k) 
{
  const double lk=log10(k);
  const double lp=gsl_spline_eval(PowSpline,lk,PowAcc);
  return Norm*pow(10.,lp-3*lk)*2.*M_PI*M_PI;
} //PowerSpec


double TopHatWindow(const double kr)
{
  const double kr2=kr*kr;
  const double kr3=kr2*kr;
  if (kr<1e-8) return 1.;
  return 3.*(sin(kr)/kr3-cos(kr)/kr2);
} //TopHatWindow


double sigma2_func(const double k,void *params)
{
  const double r_tophat=*(double *) params;
  const double W=TopHatWindow(r_tophat*k);
  return 1./(2.*M_PI*M_PI)*k*k*W*W*PowerSpec(k);
} //sigma2_func


double TopHatSigma2(double R)
{
  double result=0,abserr;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&sigma2_func;
  F.params=&R;
  gsl_integration_qag(&F,2*M_PI/500.,500.,0,1.0e-5,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //TopHatSigma2


double Mass(const double R)
{
  return 4.0*M_PI*R*R*R/3.0*Omega*rho_c;
} //Mass


double nbargal(const double m)
{
  const double lm=log10(m);
  if (lm<cutoff_fof_low || lm>cutoff_fof_high) 
    return 0.;
  else
  {
    const double value=gsl_spline_eval(FofSpline,lm,FofAcc);
    if (value>0) 
      return value;
    else 
      return 0.;
  } //else
} //nbargal


double b(const double m, const int i)
{
  switch(i) 
  {
    // case 0: 
    // { return 1.0; }
    case 1:
    {
      const double delta_c2 = delta_c*delta_c;
      const double nu2=delta_c2/Sigma2(m);
      //Tinker et al. (2010)
      const double y=log10(Delta);
      const double AT=1.0+0.24*y*exp(-pow(4./y,4));
      const double aT=0.44*y-0.88;
      const double BT=0.183;
      const double bT=1.5;
      const double CT=0.019+0.107*y+0.19*exp(-pow(4./y,4));
      const double cT=2.4;
      const double b1=1-AT*pow(nu2,0.5*aT)/(pow(nu2,0.5*aT)+pow(delta_c2,0.5*aT))+BT*pow(nu2,0.5*bT)+CT*pow(nu2,0.5*cT);
      return b1;
    }
    default:
    { return 1.0; }
  } //switch
} //b


double mugal_qawo_func(const double r,void *p)
{
  const struct mugal_qawo_params *params=(struct mugal_qawo_params *)p;
  const double pa=(params->pa);
  const double pb=(params->pb);
  const double pc=(params->pc);
  const double nm=pow(r/pb,pa)/pow(r,3)*exp(-pow(r/pb,pc));
  return r*nm;
} //mugal_qawo_func


double mugal_qawo(const double k,double m) 
{
  const double rvir=Radius(m)*Delta_invth;
  const double rvir3=rvir*rvir*rvir;
  const double pa=pow(10.,pa_eval(log10(m)));
  const double pb=pow(10.,pb_eval(log10(m)));
  const double pc=pow(10.,pc_eval(log10(m)));
  const double norm=pc/(rvir3*4*M_PI*exp(gsl_sf_lngamma(pa/pc)+log(gsl_sf_gamma_inc_P(pa/pc,pow(xmax/pb,pc)))));
  struct mugal_qawo_params params={ pa,pb,pc };
  double result=0,abserr;
  gsl_function F;
  // int status;
  F.function=&mugal_qawo_func;
  F.params=&params;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  gsl_integration_qawo_table
*t=gsl_integration_qawo_table_alloc(k*rvir,xmax,GSL_INTEG_SINE,25);
  /*status=*/ gsl_integration_qawo(&F,0.0,0,1.0e-3,WORKSIZE,w,t,&result,&abserr);
  gsl_integration_qawo_table_free(t);
  gsl_integration_workspace_free(w);
  return norm*4*M_PI*rvir*rvir*result/k;
} //mugal_qawo


double NgalF(const double m, const int j)
{ //0: ngal, 1: ngal*(ngal-1), 2: ncen, 3: nsat, 4:ncen*nsat, 5: nsat*(nsat-1)
  const double lm=log10(m);
  if (lm<cutoff_low[j] || lm>cutoff_high[j]) 
    return 0.;
  else
  {
    const double value=gsl_spline_eval(NgalSpline[j],lm,NgalAcc[j]);
    if (value>0) 
      return value;
    else
      return 0.;
  } //else
} //NgalF


double pa_eval(const double m) 
{
  const double value=gsl_spline_eval(paSpline,min(max(m,parscutoff_low),parscutoff_high),paAcc);
  if (value>=pa_low && value<=pa_high)
    return value;
  else if (value<pa_low) 
    return pa_low;
  else
    return pa_high;
} //pa_eval


double pb_eval(const double m)
{
  const double value=gsl_spline_eval(pbSpline,min(max(m,parscutoff_low),parscutoff_high),pbAcc);
  if (value>=pb_low && value<=pb_high)
    return value;
  else if (value<pb_low) 
    return pb_low;
  else 
    return pb_high;
} //pb_eval


double pc_eval(const double m)
{
  const double value=gsl_spline_eval(pcSpline,min(max(m,parscutoff_low),parscutoff_high),pcAcc);
  if (value>=pc_low && value<=pc_high)
    return value;
  else if (value<pc_low) 
    return pc_low;
  else 
    return pc_high;
} //pc_eval


double Mcensat_func(const double lm,void *p)
{
  const struct M_params *params=(struct M_params *) p;
  const double k=(params->k);
  const int i=(params->i);
  const int j=(params->j);
  const double m=exp(lm);
  if (i==0) 
  { //counterterm due to Valageas & Nishimichi (2011)
    if (j==0)
      return 0.;
#ifdef MCRIT
    else if (j==1)
      return nbargal(m)*NgalF(m,4)/pow(ngal_mean,2)*mugal_qawo(k,m)*m;
    else 
      return nbargal(m)*NgalF(m,5)/pow(ngal_mean,2)*pow(mugal_qawo(k,m),2)*m;
#else
    else if (j==1)
      return nbargal(m)*NgalF(m,4)/pow(ngal_mean,2)*(mugal_qawo(k,m)-TopHatWindow(k*Radius(m)))*m;
    else 
      return nbargal(m)*NgalF(m,5)/pow(ngal_mean,2)*(pow(mugal_qawo(k,m),2)-pow(TopHatWindow(k*Radius(m)),2))*m;
#endif
  } //if
  else
  {
    if (j==0)
      return nbargal(m)*b(m,1)*NgalF(m,2)/ngal_mean*m;
    else 
      return nbargal(m)*b(m,1)*NgalF(m,3)/ngal_mean*mugal_qawo(k,m)*m;
  } //else
} //Mcensat_func


double Mcensat(const double k, const int i, const int j)
{ //k=wavenumber, i=which haloterm [0/1], j=central/satellite [0/1/2]
  double result=0,abserr;
  struct M_params params={ k,i,j };
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&Mcensat_func;
  F.params=&params;
  gsl_integration_qag(&F,log(pow(10.,cutoff_low[3-i+j])),log(pow(10.,cutoff_high[3-i+j])),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //Mcensat


double ngal_mean_func(const double lm,void *p)
{
  const struct N_params *params=(struct N_params *) p;
  const int j=(params->j);
  const double m=exp(lm);
  return nbargal(m)*NgalF(m,j)*m; //extra m for logarithmic integration
} //ngal_mean_func


double ngal_mean_calc(const int j)
{
  double result=0,abserr;
  struct N_params params={ j };
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&ngal_mean_func;
  F.params=&params;
  gsl_integration_qag(&F,log(pow(10.,cutoff_low[j])),log(pow(10.,cutoff_high[j])),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //ngal_mean_calc


void init_power()
{
  FILE *file_;
  char file_name_[1024];
  double k,p;
  int NPowerTable=0;
  sprintf(file_name_,"%s/powrealized_rebin_corrected.dat",MCMCHaloModelDir);
  if (!(file_=fopen(file_name_,"r"))) 
  {
    char error_message_[2048];
    sprintf(error_message_,"Can't read input spectrum in file '%s'.\n",file_name_);
    terminate(error_message_);
  } //if
  do 
  {
    if (fscanf(file_," %lg %lg ",&k,&p)==2)
      NPowerTable++;
    else break;
  } //do
  while(1);
  fclose(file_);
  PowerTable=malloc(NPowerTable*sizeof(double));
  kPowerTable=malloc(NPowerTable*sizeof(double));
  file_=fopen(file_name_,"r");
  NPowerTable=0;
  do
  { //k and Delta
    if (fscanf(file_," %lg %lg ",&k,&p)==2) 
    {
      kPowerTable[NPowerTable]=k-log10(ScalePos);
      PowerTable[NPowerTable]=p;
      NPowerTable++;
    } //if
    else break;
  } //do
  while(1);
  fclose(file_);
  PowAcc=gsl_interp_accel_alloc();
  PowSpline=gsl_spline_alloc(gsl_interp_cspline,NPowerTable);
  gsl_spline_init(PowSpline,kPowerTable,PowerTable,NPowerTable);
  free(kPowerTable);
  free(PowerTable);
  NPowerTable=0;
  sprintf(file_name_,"%s/ellip_corr.dat",MCMCHaloModelDir);
  if (!(file_=fopen(file_name_,"r"))) 
  {
    char error_message_[2048];
    sprintf(error_message_,"Can't read correction spectrum in file '%s'.\n",file_name_);
    terminate(error_message_);
  } //if
  do 
  {
    if (fscanf(file_," %lg %lg ",&k,&p)==2)
      NPowerTable++;
    else break;
  } //do
  while(1);
  fclose(file_);
  PowerTable=malloc(NPowerTable*sizeof(double));
  kPowerTable=malloc(NPowerTable*sizeof(double));
  file_=fopen(file_name_,"r");
  NPowerTable=0;
  do 
  { //k and Delta
    if (fscanf(file_," %lg %lg ",&k,&p)==2)
    {
      kPowerTable[NPowerTable]=k-log10(ScalePos);
      PowerTable[NPowerTable]=p;
      NPowerTable++;
    } //if
    else break;
  } //do
  while(1);
  fclose(file_);
  ellipAcc=gsl_interp_accel_alloc();
  ellipSpline=gsl_spline_alloc(gsl_interp_cspline,NPowerTable);
  gsl_spline_init(ellipSpline,kPowerTable,PowerTable,NPowerTable);
  free(kPowerTable);
  free(PowerTable);
  NPowerTable=0;
  sprintf(file_name_,"%s/align_corr.dat",MCMCHaloModelDir);
  if(!(file_=fopen(file_name_,"r"))) 
  {
    printf("Can't read correction spectrum in file '%s'.\n",file_name_);
    exit(0);
  } //if
  do
  {
    if (fscanf(file_," %lg %lg ",&k,&p)==2) 
      NPowerTable++;
    else break;
  } //do
  while(1);
  fclose(file_);
  PowerTable=malloc(NPowerTable*sizeof(double));
  kPowerTable=malloc(NPowerTable*sizeof(double));
  file_=fopen(file_name_,"r");
  NPowerTable=0;
  do 
  { //k and Delta
    if (fscanf(file_," %lg %lg ",&k,&p)==2) 
    {
      kPowerTable[NPowerTable]=k-log10(ScalePos);
      PowerTable[NPowerTable]=p;
      NPowerTable++;
    } //if
    else break;
  } //do
  while(1);
  fclose(file_);
  alignAcc=gsl_interp_accel_alloc();
  alignSpline=gsl_spline_alloc(gsl_interp_cspline,NPowerTable);
  gsl_spline_init(alignSpline,kPowerTable,PowerTable,NPowerTable);
  free(kPowerTable);
  free(PowerTable);
} //init_power


void init_sigma() 
{
  int i,NSigmaTable=500;
  double MSigmaTable[500],SigmaTable[500];
  double R;
  for (i=0; i<NSigmaTable; i++) 
  {
    MSigmaTable[i]=i*15./(float)NSigmaTable+5.;
    R=Radius(pow(10.,MSigmaTable[i]));
    SigmaTable[i]=log10(TopHatSigma2(R));
  } //for
  SigmaAcc=gsl_interp_accel_alloc();
  SigmaSpline=gsl_spline_alloc(gsl_interp_cspline,NSigmaTable);
  gsl_spline_init(SigmaSpline,MSigmaTable,SigmaTable,NSigmaTable);
} //init_sigma


void init_numgal(const float masslimit_low, const float masslimit_high, const int snap)
{
  const int NMassTable2=6;
  const double massoffset=log10(2*(Gravity/1e10)/(Delta*Omega*1e4));
  
  int i,j,jj;
  int mbin,mbin2,k,found,ncen,nsat;
  int nsat_tot=0;
  
  int *NgalTotal,*NfofTotal,*NgalInRange,*usedbymass;
  int status;
  double *MassTable;
  double *masstmp,*patmp,*pbtmp,*pctmp;
  double **NgalTable;
  double *borders,*radii,*p,*perror,*pberr,*paerr,*pcerr;
  double boxsize=BoxSize;
  double r,rvir,relx,rely,relz;
  double *MassTable2,*pa_m,*pb_m,*pc_m;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *abc;
  gsl_multimin_function_fdf my_func;
  size_t iter;
  for (i=0; i<6; ++i) 
  {
    cutoff_low[i]=0;
    cutoff_high[i]=0;
  } //for
  NgalAcc=malloc(6*sizeof(gsl_interp_accel*));
  NgalSpline=malloc(6*sizeof(gsl_spline*));
  usedbymass=malloc(massbins*sizeof(int));
  for (i=0; i<massbins; ++i)
  { usedbymass[i]=0; }
  for (i=0; i<6; ++i)
  {
    NgalAcc[i]=gsl_interp_accel_alloc();
    NgalSpline[i]=gsl_spline_alloc(gsl_interp_cspline,massbins);
  } //for
  borders=malloc((NMassTable2+1)*sizeof(double));
  MassTable=malloc(massbins*sizeof(double));
  NgalTable=malloc(6*sizeof(double*));
  for (i=0; i<6; ++i)
  { NgalTable[i]=malloc(massbins*sizeof(double)); }
  for (j=0; j<massbins; ++j)
  {
    MassTable[j]=0;
    for (i=0; i<6; ++i) NgalTable[i][j]=0;
  } //for
  for (mbin=0; mbin<massbins; ++mbin)
  {
    for (j=0; j<UsedFofsInSample[snap]; ++j)
    {
      ncen=0;
      nsat=0;
      if (MCMC_FOF2[j].M_Mean200 > minfofmass+mbin*(maxfofmass-minfofmass)/(double)massbins && MCMC_FOF2[j].M_Mean200 < minfofmass+(mbin+1)*(maxfofmass-minfofmass)/(double)massbins)
      {
        i=MCMC_FOF2[j].IndexOfCentralGal;
        if (i>=0) 
        {
          if (MCMC_GAL[snap][HashTable[i]].StellarMass>masslimit_low && MCMC_GAL[snap][HashTable[i]].StellarMass<=masslimit_high)
            ncen=1;
          for (k=1; k<MCMC_GAL[snap][HashTable[i]].ngal; ++k)
          {
            if (MCMC_GAL[snap][HashTable[i+k]].StellarMass>masslimit_low && MCMC_GAL[snap][HashTable[i+k]].StellarMass<=masslimit_high) 
              nsat++;
          } //for
        } //if
        MassTable[mbin]+=MCMC_FOF2[j].M_Mean200;
        usedbymass[mbin]++;
      } //if
      NgalTable[0][mbin]+=ncen+nsat;
      NgalTable[1][mbin]+=(ncen+nsat)*max((ncen+nsat)-1,0);
      NgalTable[2][mbin]+=ncen;
      NgalTable[3][mbin]+=nsat;
      NgalTable[4][mbin]+=ncen*nsat;
      NgalTable[5][mbin]+=nsat*max(nsat-1,0);
      nsat_tot+=nsat;
    } //for
  } //for
  for (j=0; j<massbins; ++j)
  {
    if (usedbymass[j]>0) 
      MassTable[j]/=usedbymass[j];
    else 
      MassTable[j]=minfofmass+(j+0.5)*(maxfofmass-minfofmass)/(double)massbins;
  } //for
  nsat=NgalTable[3][0];
  if (NgalTable[0][0]>0) 
  {
    char error_message_[2048];
    sprintf(error_message_,"First Ngal mass bin not equal to zero, divergence will ensue.\n");
    terminate(error_message_);
  } //if
  borders[0]=minfofmass;
  borders[NMassTable2]=maxfofmass;
  for (i=1; i<6; ++i)
  { borders[i]=10+i; }
  for (j=0; j<massbins; ++j)
  {
    if (usedbymass[j]>0) 
    {
      for (i=0; i<6; ++i)
      { NgalTable[i][j]/=usedbymass[j]; }
    } //if
  } //for
  for (i=0; i<6; ++i) 
  {
    for (mbin=0; mbin<massbins-1 && cutoff_low[i]==0; ++mbin)
    {
      if (NgalTable[i][mbin]==0 && NgalTable[i][mbin+1]>0)
      { cutoff_low[i]=MassTable[mbin]; } 
    } //for
    for (mbin=massbins-1; mbin>0 && cutoff_high[i]==0; --mbin)
    {
      if (NgalTable[i][mbin-1]>0 && NgalTable[i][mbin]==0) 
      { cutoff_high[i]=MassTable[mbin]; }
    } //for
  } //for
  for (i=0; i<6; ++i) 
  { gsl_spline_init(NgalSpline[i],MassTable,NgalTable[i],massbins); }
  for (i=5; i>=0; --i)
  { free(NgalTable[i]); }
  free(NgalTable);
  free(MassTable);
  MassTable2=malloc(NMassTable2*sizeof(double));
  NgalTotal=malloc(NMassTable2*sizeof(int));
  NfofTotal=malloc(NMassTable2*sizeof(int));
  NgalInRange=malloc(NMassTable2*sizeof(int));
  for (i=0; i<NMassTable2; ++i)
  {
    MassTable2[i]=0;
    NgalTotal[i]=0;
    NfofTotal[i]=0;
    NgalInRange[i]=0;
  } //for
  for (j=0; j<UsedFofsInSample[snap]; ++j)
  {
    i=MCMC_FOF2[j].IndexOfCentralGal;
    if (i>=0) 
    {
      mbin2=-1;
      for (jj=1; jj<=NMassTable2; ++jj)
      {
        if (MCMC_GAL[snap][HashTable[i]].M_Mean200<borders[jj])
        {
          mbin2=jj-1;
          break;
        } //if
      } //for
      if (mbin2>=0 && mbin2<NMassTable2 && MCMC_GAL[snap][HashTable[i]].ngal>1)
      {
        found=0;
        rvir=pow(10.,(MCMC_GAL[snap][HashTable[i]].M_Mean200+massoffset)/3.);
        for (jj=1; jj<MCMC_GAL[snap][HashTable[i]].ngal; ++jj)
        {
          if (MCMC_GAL[snap][HashTable[i+jj]].StellarMass>masslimit_low && MCMC_GAL[snap][HashTable[i+jj]].StellarMass<=masslimit_high) {
            found++;
            relx=min(fabs(MCMC_GAL[snap][HashTable[i+jj]].x-MCMC_GAL[snap][HashTable[i]].x),fabs(boxsize-fabs(MCMC_GAL[snap][HashTable[i+jj]].x-MCMC_GAL[snap][HashTable[i]].x)));
            rely=min(fabs(MCMC_GAL[snap][HashTable[i+jj]].y-MCMC_GAL[snap][HashTable[i]].y),fabs(boxsize-fabs(MCMC_GAL[snap][HashTable[i+jj]].y-MCMC_GAL[snap][HashTable[i]].y)));
            relz=min(fabs(MCMC_GAL[snap][HashTable[i+jj]].z-MCMC_GAL[snap][HashTable[i]].z),fabs(boxsize-fabs(MCMC_GAL[snap][HashTable[i+jj]].z-MCMC_GAL[snap][HashTable[i]].z)));
            r=sqrt(relx*relx+rely*rely+relz*relz);
            if (r/rvir>xmin && r/rvir<=xmax) 
            { NgalInRange[mbin2]++; }
          } //if
        } //for
        if (found>0) 
        {
          NgalTotal[mbin2]+=found;
          NfofTotal[mbin2]++;
          MassTable2[mbin2]+=found*MCMC_GAL[snap][HashTable[i]].M_Mean200;
        } //if
      } //if
    } //if
  } //for
  for (mbin2=0; mbin2<NMassTable2; ++mbin2) 
  {
    if (NgalTotal[mbin2]>0) 
    { MassTable2[mbin2]/=NgalTotal[mbin2]; }
    else 
    { MassTable2[mbin2]=0.5*(borders[mbin2]+borders[mbin2+1]); }
  } //for
  pa_m=malloc(NMassTable2*sizeof(double));
  pb_m=malloc(NMassTable2*sizeof(double));
  pc_m=malloc(NMassTable2*sizeof(double));
  paerr=malloc(NMassTable2*sizeof(double));
  pberr=malloc(NMassTable2*sizeof(double));
  pcerr=malloc(NMassTable2*sizeof(double));
  p=malloc(3*sizeof(double));
  perror=malloc(3*sizeof(double));
  my_func.n=3;
  my_func.f=&my_f;
  my_func.df=&my_df;
  my_func.fdf=&my_fdf;
  abc=gsl_vector_alloc(3);
  T=gsl_multimin_fdfminimizer_conjugate_pr;
  s=gsl_multimin_fdfminimizer_alloc(T,3);
  for (mbin2=0; mbin2<NMassTable2; ++mbin2)
  {
    if (NgalInRange[mbin2]>2)
    {
      numrad=NgalInRange[mbin2];
      radii=malloc(numrad*sizeof(double));
      NgalInRange[mbin2]=0;
      for (j=0; j<UsedFofsInSample[snap]; ++j) 
      {
        i=MCMC_FOF2[j].IndexOfCentralGal;
        if (i>=0) 
        {
          if (MCMC_GAL[snap][HashTable[i]].M_Mean200 >= borders[mbin2]   && 
              MCMC_GAL[snap][HashTable[i]].M_Mean200 <  borders[mbin2+1] &&
              MCMC_GAL[snap][HashTable[i]].ngal      >  1)
          {
            rvir=pow(10.,(MCMC_GAL[snap][HashTable[i]].M_Mean200+massoffset)/3.);
            for (jj=1; jj<MCMC_GAL[snap][HashTable[i]].ngal; ++jj) 
            {
              if (MCMC_GAL[snap][HashTable[i+jj]].StellarMass>masslimit_low && MCMC_GAL[snap][HashTable[i+jj]].StellarMass<=masslimit_high)
              {
                relx=min(fabs(MCMC_GAL[snap][HashTable[i+jj]].x-MCMC_GAL[snap][HashTable[i]].x),fabs(boxsize-fabs(MCMC_GAL[snap][HashTable[i+jj]].x-MCMC_GAL[snap][HashTable[i]].x)));
                rely=min(fabs(MCMC_GAL[snap][HashTable[i+jj]].y-MCMC_GAL[snap][HashTable[i]].y),fabs(boxsize-fabs(MCMC_GAL[snap][HashTable[i+jj]].y-MCMC_GAL[snap][HashTable[i]].y)));
                relz=min(fabs(MCMC_GAL[snap][HashTable[i+jj]].z-MCMC_GAL[snap][HashTable[i]].z),fabs(boxsize-fabs(MCMC_GAL[snap][HashTable[i+jj]].z-MCMC_GAL[snap][HashTable[i]].z)));
                r=sqrt(relx*relx+rely*rely+relz*relz);
                if (r/rvir>xmin && r/rvir<=xmax) 
                {
                  radii[NgalInRange[mbin2]]=r/rvir;
                  NgalInRange[mbin2]++;
                } //if
              } //if
            } //for
          } //if
        } //if
      } //for
      p[0]=0.0;
      p[1]=-0.4+0.12*(ThisTask%6);
      p[2]=0.1*(ThisTask%6);
      gsl_vector_set(abc,0,p[0]);
      gsl_vector_set(abc,1,p[1]);
      gsl_vector_set(abc,2,p[2]);
      my_func.params=(void *)radii;
      gsl_multimin_fdfminimizer_set(s,&my_func,abc,1e-4,1e-2); //stepsize, gradient tolerance
      iter=0;
      do
      {
        iter++;
        status=gsl_multimin_fdfminimizer_iterate(s);
        if (status) break;
        status=gsl_multimin_test_gradient(s->gradient,1e-2);
        gsl_vector_set(s->x,0,max(min(gsl_vector_get(s->x,0),pa_high),pa_low));
        gsl_vector_set(s->x,1,max(min(gsl_vector_get(s->x,1),pb_high),pb_low));
        gsl_vector_set(s->x,2,max(min(gsl_vector_get(s->x,2),pc_high),pc_low));
      } //do
      while (status==GSL_CONTINUE && iter<500);
      p[0]=max(min(gsl_vector_get(s->x,0),pa_high),pa_low);
      p[1]=max(min(gsl_vector_get(s->x,1),pb_high),pb_low);
      p[2]=max(min(gsl_vector_get(s->x,2),pc_high),pc_low);
      perror[0]=0.1;
      perror[1]=0.1;
      perror[2]=0.1;
      paramerror(radii,p,perror);
      pa_m[mbin2]=p[0];
      pb_m[mbin2]=p[1];
      pc_m[mbin2]=p[2];
      paerr[mbin2]=perror[0];
      pberr[mbin2]=perror[1];
      pcerr[mbin2]=perror[2];
      if (paerr[mbin2]<0) paerr[mbin2]=0.1;
      if (pberr[mbin2]<0) pberr[mbin2]=0.1;
      if (pcerr[mbin2]<0) pcerr[mbin2]=0.1;
      free(radii);
    } //if
  } //for
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(abc);
  free(perror);
  free(p);
  masstmp=malloc((NMassTable2+2)*sizeof(double));
  patmp=malloc((NMassTable2+2)*sizeof(double));
  pbtmp=malloc((NMassTable2+2)*sizeof(double));
  pctmp=malloc((NMassTable2+2)*sizeof(double));
  masstmp[0]=minfofmass;
  masstmp[NMassTable2+1]=maxfofmass;
  patmp[0]=pa_m[0];
  pbtmp[0]=pb_m[0];
  pctmp[0]=pc_m[0];
  patmp[NMassTable2+1]=pa_m[NMassTable2-1];
  pbtmp[NMassTable2+1]=pb_m[NMassTable2-1];
  pctmp[NMassTable2+1]=pc_m[NMassTable2-1];
  for (i=1; i<=NMassTable2; ++i) 
  {
    masstmp[i]=MassTable2[i-1];
    patmp[i]=pa_m[i-1];
    pbtmp[i]=pb_m[i-1];
    pctmp[i]=pc_m[i-1];
  } //for
  paAcc=gsl_interp_accel_alloc();
  if (NMassTable2>2)
    paSpline=gsl_spline_alloc(gsl_interp_akima,NMassTable2+2);
  else 
    paSpline=gsl_spline_alloc(gsl_interp_cspline,NMassTable2+2);
  pbAcc=gsl_interp_accel_alloc();
  if (NMassTable2>2) 
    pbSpline=gsl_spline_alloc(gsl_interp_akima,NMassTable2+2);
  else
    pbSpline=gsl_spline_alloc(gsl_interp_cspline,NMassTable2+2);
  pcAcc=gsl_interp_accel_alloc();
  if (NMassTable2>2) 
    pcSpline=gsl_spline_alloc(gsl_interp_akima,NMassTable2+2);
  else 
    pcSpline=gsl_spline_alloc(gsl_interp_cspline,NMassTable2+2);
  gsl_spline_init(paSpline,masstmp,patmp,NMassTable2+2);
  gsl_spline_init(pbSpline,masstmp,pbtmp,NMassTable2+2);
  gsl_spline_init(pcSpline,masstmp,pctmp,NMassTable2+2);
  parscutoff_low=minfofmass;
  parscutoff_high=maxfofmass;
  free(pctmp);
  free(pbtmp);
  free(patmp);
  free(masstmp);
  free(pcerr);
  free(pberr);
  free(paerr);
  free(pc_m);
  free(pb_m);
  free(pa_m);
  free(NgalInRange);
  free(NfofTotal);
  free(NgalTotal);
  free(MassTable2);
  free(borders);
} //init_numgal


void initialize_halomodel(void)
 {
  int i,NnuTable=100;
  double res;
  double MnuTable[100],nuTable[100];
  FILE *file_;
  char file_name_[1024];
  double m,h;
  double FofmassTable[massbins],FofnumTable[massbins];
  
  rho_c=RhoCrit*1e10;
  delta_c=(3./5.*pow(0.5*3*M_PI,2./3.))*(1-0.0123*log10(1+(1./Omega-1)));
  Delta_invth=1./pow(Delta,1./3.);
  rho_mean=Omega*rho_c;
  Norm=1.0;
  init_power();
  res=TopHatSigma2(8.);
  Norm=Sigma8*Sigma8/res;
  res=TopHatSigma2(8.);
  init_sigma();
  for (i=0; i<NnuTable; i++) 
  {
    MnuTable[i]=log(pow(10,i*3./(float)NnuTable+11.));
    nuTable[i]=delta_c/sqrt(Sigma2(exp(MnuTable[i])));
  } //for
  nuAcc=gsl_interp_accel_alloc();
  nuSpline=gsl_spline_alloc(gsl_interp_cspline,NnuTable);
  gsl_spline_init(nuSpline,nuTable,MnuTable,NnuTable);
  Mstar=exp(gsl_spline_eval(nuSpline,1.0,nuAcc));
  gsl_spline_free(nuSpline);
  gsl_interp_accel_free(nuAcc);
  Rstar=Radius(Mstar);
  //sprintf(file_name_,"%s/fofnum.dat",MCMCHaloModelDir);
#ifdef MCRIT
  sprintf(file_name_,"%s/fofnum_m200_z0.02.dat",MCMCHaloModelDir);
#else
  sprintf(file_name_,"%s/fofnum_m200mean.dat",MCMCHaloModelDir);
#endif
  if (!(file_=fopen(file_name_,"r"))) 
  {
    char error_message_[2048];
    sprintf(error_message_,"Can't read input FoF mass function in file '%s'.\n",file_name_);
    terminate(error_message_);
  } //if
  i=0;
  do
  {
    if (fscanf(file_," %lg %lg ",&m,&h)==2)
      i++;
    else break;
    FofmassTable[i-1]=m+log10(ScaleMass);
    FofnumTable[i-1]=h/pow(ScalePos,3);
  } //do
  while(1);
  fclose(file_);
  if (FofnumTable[0]>0) 
  {
    char error_message_[2048];
    sprintf(error_message_,"First FoF mass bin not equal to zero, divergence will ensue.\n");
    terminate(error_message_);
  } //if
  cutoff_fof_low=0.;
  cutoff_fof_high=0.;
  for (i=0; i<massbins && cutoff_fof_low==0; ++i)
  {
    if (FofnumTable[i]==0 && FofnumTable[i+1]>0)
    { cutoff_fof_low=FofmassTable[i]; }
  } //for
  for (i=massbins-1; i>0 && cutoff_fof_high==0; --i)
  {
    if (FofnumTable[i-1]>0 && FofnumTable[i]==0) 
    { cutoff_fof_high=FofmassTable[i]; }
  } //for
  for (i=0; i<massbins; ++i) 
  {  FofnumTable[i] /= pow(10.,FofmassTable[i])*log(10.); }
  FofAcc=gsl_interp_accel_alloc();
  FofSpline=gsl_spline_alloc(gsl_interp_cspline,massbins);
  gsl_spline_init(FofSpline,FofmassTable,FofnumTable,massbins);
} //initialize_halomodel


double my_f(const gsl_vector *v,void *params)
 {
  int i;
  double a,b,c,tot;
  double *p=(double *)params;
  if (gsl_vector_get(v,0)<pa_low) a=pow(10.,pa_low);
  else if (gsl_vector_get(v,0)>pa_high) a=pow(10.,pa_high);
  else a=pow(10.,gsl_vector_get(v,0));
  if (gsl_vector_get(v,1)<pb_low) b=pow(10.,pb_low);
  else if (gsl_vector_get(v,1)>pb_high) b=pow(10.,pb_high);
  else b=pow(10.,gsl_vector_get(v,1));
  if (gsl_vector_get(v,2)<pc_low) c=pow(10.,pc_low);
  else if (gsl_vector_get(v,2)>pc_high) c=pow(10.,pc_high);
  else c=pow(10.,gsl_vector_get(v,2));
  tot=numrad*(log(c)-a*log(b)-gsl_sf_lngamma(a/c));
  for (i=0; i<numrad; ++i) tot+=(a-3)*log(p[i])-pow(p[i]/b,c);
  return -tot;
} //my_f


/* The gradient of f, df = (df/dx, df/dy). */
void my_df(const gsl_vector *v,void *params,gsl_vector *df)
 {
  int i;
  double a,b,c,tot;
  double *p=(double *)params;
  if (gsl_vector_get(v,0)<pa_low) a=pow(10.,pa_low);
  else if (gsl_vector_get(v,0)>pa_high) a=pow(10.,pa_high);
  else a=pow(10.,gsl_vector_get(v,0));
  if (gsl_vector_get(v,1)<pb_low) b=pow(10.,pb_low);
  else if (gsl_vector_get(v,1)>pb_high) b=pow(10.,pb_high);
  else b=pow(10.,gsl_vector_get(v,1));
  if (gsl_vector_get(v,2)<pc_low) c=pow(10.,pc_low);
  else if (gsl_vector_get(v,2)>pc_high) c=pow(10.,pc_high);
  else c=pow(10.,gsl_vector_get(v,2));
  tot=-numrad*(a/c*gsl_sf_psi(a/c));
  for (i=0; i<numrad; ++i) tot+=a*log(p[i]/b);
  gsl_vector_set(df,0,-tot);
  tot=-numrad*a;
  for (i=0; i<numrad; ++i) tot+=c*pow(p[i]/b,c);
  gsl_vector_set(df,1,-tot);
  tot=numrad*(1+a/c*gsl_sf_psi(a/c));
  for (i=0; i<numrad; ++i) tot-=c*log(p[i]/b)*pow(p[i]/b,c);
  gsl_vector_set(df,2,-tot);
} //my_df


/* Compute both f and df together. */
void my_fdf(const gsl_vector *x,void *params,double *f, gsl_vector *df)
 {
  *f = my_f(x,params); 
  my_df(x, params, df);
} //my_fdf


void paramerror(double *x,double *p,double *perror)
 {
  int i;
  
  const double a=pow(10.,p[0]);
  const double b=pow(10.,p[1]);
  const double c=pow(10.,p[2]);
  const double polygamma0=gsl_sf_psi(a/c);
  
  double dblderiv=0;
  for (i=0; i<numrad; ++i) dblderiv+=pow(a/c*polygamma0-a*log(x[i]/b),2);
  perror[0]=pow(1./(sqrt(2*M_PI)*dblderiv),1./3.);
  dblderiv=0;
  for (i=0; i<numrad; ++i) dblderiv+=pow(a-c*pow(x[i]/b,c),2);
  perror[1]=pow(1./(sqrt(2*M_PI)*dblderiv),1./3.);
  dblderiv=0;
  for (i=0; i<numrad; ++i) dblderiv+=pow(1-c*log(x[i]/b)*pow(x[i]/b,c)+a/c*polygamma0,2);
  perror[2]=pow(1./(sqrt(2*M_PI)*dblderiv),1./3.);
} //paramerror


int poissonfit(int m,int n,double *p,double *dy,double **dvec,void *vars)
 {
  (void) m; /* suppress unused-parameter warning */
  (void) dvec; /* suppress unused-parameter warning */
  
  const double a=pow(10.,p[0]);
  const double b=pow(10.,p[1]);
  const double c=pow(10.,p[2]);
  const double *x=(double *) vars; //should be r/rvir for every satellite with xmin<r/rvir<xmax
  const double digam=gsl_sf_psi(a/c)/c;
  
  int i;
  for (i=0; i<n; ++i)
  {  dy[i]=0; }
  for (i=0; i<numrad; ++i)
  {
    dy[0]+=log(x[i]/b)-digam;
    dy[1]+=c/b*pow(x[i]/b,c)-a/b;
    dy[2]+=1./c+a/c*digam-log(x[i]/b)*pow(x[i]/b,c);
  } //for
  return 0;
} //poissonfit
