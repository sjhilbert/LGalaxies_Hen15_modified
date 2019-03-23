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

/** @file   mcmc_halomodel.h
 *  @date   2013-2019
 *  @author Marcel van Daalen
 *  @author Stefan Hilbert
 *
 *  @brief  prototypes for functions for halo model for mcmc.
 */

void halomodel(const int snap, const float masslimit_low, const float masslimit_high, const int number_of_r, double* r_arr,double* proj_arr);
double TwoPowerSpec(const double k_, const int censat);
double pconv_W_P_func(const double theta_,void *p_);
double pconv_W_func(const double lq,void *p);
double pconv_W(const double k, const double R, const int censat);
double calc_mean_rhalo_simple_func(const double lm, void *p);
double calc_mean_rhalo_simple(const int censat, const int norm);
double calc_mean_rsat_func(const double lr,void *p);
double calc_mean_rsat(const double lm, const int extrar);
double NewPowerSpec(const double k);
double corr_qawo_func(const double k, void *params);
double corr_qawo(const double r, const double a, const double L);
double proj_corr_func(const double r,void *p);
double proj_corr(double sigma);
double Radius(const double m);
double Sigma2(const double m);
double PowerSpec(const double k);
double sigma2_func(const double k,void *params);
double TopHatSigma2(const double R);
double TopHatWindow(const double kr);
double nbargal(const double m);
double b(const double m,const int i);
double mugal_qawo_func(const double r,void *p);
double mugal_qawo(const double k, const double m);
double NgalF(const double m, const int j);
double pa_eval(const double m);
double pb_eval(const double m);
double pc_eval(const double m);
double Mcensat_func(const double lm,void *p);
double Mcensat(const double k, const int i, const int j);
double ngal_mean_func(const double lm,void *p);
double ngal_mean_calc(const int j);
void init_power(void);
void init_sigma(void);
void init_numgal(const float masslimit_low, const float masslimit_high, const int snap);
void initialize_halomodel(void);
double my_f(const gsl_vector *v,void *params);
void my_df(const gsl_vector *v,void *params,gsl_vector *df);
void my_fdf(const gsl_vector *x,void *params,double *f,gsl_vector *df);
void paramerror(double *x,double *p,double *perror);
int poissonfit(int m,int n,double *p,double *dy,double **dvec,void *vars);
