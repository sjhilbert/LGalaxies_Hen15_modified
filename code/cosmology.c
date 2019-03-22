/*  Copyright (C) <2016+>  <L-Galaxies>
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

/** @file cosmology.c
 *  @brief   code for cosmological distances etc.
 * 
 *  @author  Stefan Hilbert
 *  @author  Rachel Asquith
 *  @author  auhtor of time_to_present from age.c
 * 
 *  @date    2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "proto.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
// #include <gsl/gsl_interp.h>
// #include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_hyperg.h>

#include "allvars.h"

/* workspace for gsl integration:*/
#ifndef COSMOLOGY_GSL_INTEGRATION_WORKSPACE_SIZE
#define COSMOLOGY_GSL_INTEGRATION_WORKSPACE_SIZE 1000
#endif /* not defined COSMOLOGY_GSL_INTEGRATION_WORKSPACE_SIZE */

static gsl_integration_workspace *cosmology_gsl_integration_workspace_;


#ifdef ASSUME_FLAT_LCDM
/** @brief checks for flatness, i.e. 1 = Omega + OmegaLambda */
void assert_flat_LCDM(void)
{
  if(fabs(1 - Omega - OmegaLambda) > 0.01)
  {
    printf("code compiled with ASSUME_FLAT_LCDM, but input density paramters Omega = %f and OmegaLambda = %f don't add up to 1.", Omega, OmegaLambda);
    terminate("code compiled with ASSUME_FLAT_LCDM, but input density paramters Omega and OmegaLambda don't add up to 1.");
  }
}
#endif /* defined ASSUME_FLAT_LCDM */


/** @brief inits arrays for computing redshifts from distances
 *         by interpolation
 */
void init_redshift_for_comoving_distance(void)
{
  int i_;
  set_interpolation_tables(i_, 0, N_REDSHIFTS_FOR_INTERPOLATION, 0., DELTA_REDSHIFT_FOR_INTERPOLATION, comoving_los_distance_for_redshift, redshift_table_for_interpolation, distance_table_for_interpolation); 
}


/** @brief inits workspace for gsl integration */
void init_cosmology_gsl_integration(void)
{
  cosmology_gsl_integration_workspace_ = gsl_integration_workspace_alloc(COSMOLOGY_GSL_INTEGRATION_WORKSPACE_SIZE);
}


/** @brief inits stuff for cosmological calculations */
void init_cosmology(void)
{
  init_redshift_for_comoving_distance();
  init_cosmology_gsl_integration();
}


/** @brief aux. function for computing cosmic distances */
static inline double 
integrand_comoving_los_distance_for_redshift_(double inv_a_, void *param_)
{
  (void)param_;  /* avoid unused-parameter warning */
  return 1 / sqrt(Omega * inv_a_ * inv_a_ * inv_a_ + (1 - Omega - OmegaLambda) * inv_a_ * inv_a_ + OmegaLambda);
}


/** @brief compute comoving los distance for given cosmological redshift
 *         using numerical integrationg */
static inline double 
comoving_los_distance_for_redshift_ni(const double redshift_)
{
  gsl_function F_;
  double result_, abserr_;
 
  F_.function = &integrand_comoving_los_distance_for_redshift_;

  gsl_integration_qag(&F_, 1.0, 1.0 + redshift_, 1e-5, 1.0e-6, COSMOLOGY_GSL_INTEGRATION_WORKSPACE_SIZE, GSL_INTEG_GAUSS21, cosmology_gsl_integration_workspace_, &result_, &abserr_);

  return D_HUBBLE * result_;
}


/** @brief   compute comoving los distance for given cosmological redshift 
 *           using hypergeometric functions
 *
 *  @param [in] redshift_ cosmological redshift for which distance to compute
 *  @return l.o.s. comoving distance in simulation units
 *         (i.e. units of D_HUBBLE, usually Mpc/h or kpc/h)
 *
 *  @warning assumes flat LCDM with Omega + OmegaLambda = 1
 */
static inline double
comoving_los_distance_for_redshift_flat_LCDM(const double redshift_)
{
  const double inv_a_    =  1. + redshift_;
  const double inv_omlf_ =  1. / (OmegaLambda + inv_a_ * inv_a_ * inv_a_ * Omega);
  const double result_   = 
  (0.99 < OmegaLambda * inv_omlf_) ?
    D_HUBBLE * redshift_ :
    D_HUBBLE * (  2.                            * gsl_sf_hyperg_2F1(1./2., 1., 7./6., OmegaLambda            )
                - 2. * inv_a_ * sqrt(inv_omlf_) * gsl_sf_hyperg_2F1(1./2., 1., 7./6., OmegaLambda * inv_omlf_));
  return result_;
}


/** @brief compute comoving los distance for given cosmological redshift
 *
 *  @param [in] redshift_ cosmological redshift for which distance to compute
 *  @return l.o.s. comoving distance in simulation units
 *         (i.e. units of D_HUBBLE, usually Mpc/h or kpc/h)
 */
double
comoving_los_distance_for_redshift(const double redshift_)
{
#ifdef ASSUME_FLAT_LCDM
  return comoving_los_distance_for_redshift_flat_LCDM(redshift_);
#else /* not defined ASSUME_FLAT_LCDM */
 return comoving_los_distance_for_redshift_ni(redshift_);
#endif /* not defined ASSUME_FLAT_LCDM */
}


/** @brief   compute luminosity distance for given 
 *           cosmological redshift  */
double
luminosity_distance_for_redshift(const double redshift_)
{
#ifdef ASSUME_FLAT_LCDM
  return (1. + redshift_) * comoving_los_distance_for_redshift_flat_LCDM(redshift_);
#else /* not defined ASSUME_FLAT_LCDM */
  const double OmegaK_ = 1. - Omega - OmegaLambda;
  if(OmegaK_ < -0.01)
  {
    const double R_K_ = D_HUBBLE / sqrt(-OmegaK_);
    return (1. + redshift_) * R_K_ * sinh(comoving_los_distance_for_redshift_ni(redshift_) / R_K_);
  }
  else if(OmegaK_ > 0.01)
  {
    const double R_K_ = D_HUBBLE / sqrt(OmegaK_);
    return (1. + redshift_) * R_K_ * sin(comoving_los_distance_for_redshift_ni(redshift_) / R_K_);
  }
  else
  { return (1. + redshift_) * comoving_los_distance_for_redshift_flat_LCDM(redshift_); }
#endif /* not defined ASSUME_FLAT_LCDM */
}


/** @brief compute cosmological redshift for given 
 *         comoving los distance
 * 
 *  @details currently uses linear interpolation and 
 *  zeroth-order extrapolation
 * 
 *  @warning requires init_redshift_for_comoving_distance before
 */
double redshift_for_comoving_los_distance(const double d_)
{
  double res_;
  int i_;
  linear_interpolate(i_, 0, N_REDSHIFTS_FOR_INTERPOLATION, d_, distance_table_for_interpolation, redshift_table_for_interpolation, res_, <, linear);    
  return res_;
}


/** @brief compute doppler redshift for given 
 *         los peculiar velocity
 */
double redshift_for_radial_velocity(const double v_)
{
  const double beta_ = v_ * (1. / SpeedOfLight); /* velocity in units of speed of light */
  if(beta_ < 0.03)
  { return beta_ + 0.5 * beta_ * beta_ + 0.5 * beta_ * beta_ * beta_; }
  else
  { return sqrt((1. + beta_) / (1. - beta_)) - 1.; }
}

/** @brief combines two redshifts into single redshift
 *
 * assumes as composition law:
 * (1 + z_res) = (1 + z_1) (1 + z_2)
 */
double combine_redshifts(const double z_1_, const double z_2_)
{ return z_1_ + z_2_ + z_1_ * z_2_; }


/** @brief returns combined redshift from cosmological distance and peculiar velocity
 */
double redshift_for_comoving_los_distance_and_radial_velocity(const double d_, const double v_)
{
  const double cosmological_redshift_ = redshift_for_comoving_los_distance(d_);
  const double velocity_redshift_     = redshift_for_radial_velocity(v_);
  const double combined_redshift_     = combine_redshifts(cosmological_redshift_, velocity_redshift_);
  return combined_redshift_;
}


/** @brief aux. function for computing cosmic ages */
static inline double 
integrand_time_to_present(double a_, void *param_)
{
  (void)param_;  /* avoid unused-parameter warning */
  return 1 / sqrt(Omega / a_ + (1 - Omega - OmegaLambda) + OmegaLambda * a_ * a_);
}


/** @brief For a given redshift, returns time to present for a given redshift
 *  using numerical integration
 *
 * @eturn time to present in code units/Hubble_h
 */
static inline double
time_to_present_ni(const double redshift_)
{
  gsl_function F_;
  double result_, abserr_;
 
  F_.function = &integrand_time_to_present;

  /* I think that the 1.0/Hubble here should be any small number, such
     as 1e-6 - it is the aboslute error accuracy required.  PAT */
  gsl_integration_qag(&F_, 1.0 / (redshift_ + 1), 1.0, 1.0 / Hubble,
                      1.0e-8, COSMOLOGY_GSL_INTEGRATION_WORKSPACE_SIZE, GSL_INTEG_GAUSS21, cosmology_gsl_integration_workspace_, &result_, &abserr_);

  return 1 / Hubble * result_;
}


/** @brief For a given redshift, returns time to present for a given redshift
 *  assuming flat LCDM using analytic formula
 *
 * @return time to present (lookback time) in code units/Hubble_h
 *
 * \f$ t(z) = \frac{2}{3 H_0 sqrt{\Omega_{\Lambda}}} \left[\mathrm{asinh}\sqrt{\frac{\Omega_{\Lambda}}{\Omega_{m}}} - \mathrm{asinh}\sqrt{\frac{\Omega_{\Lambda}}{\Omega_{m}(1+z)^3}}\right] \f$
 *
 * @bug (corrected on 2018-02.17) forgot sqrt around first OmegaLambda 
 */
static inline double
time_to_present_flat_LCDM(const double redshift_)
{
  const double a_ =  1. / (1. + redshift_);
  return (1. / Hubble) * (2./ 3.) / sqrt(OmegaLambda) * (asinh(sqrt(OmegaLambda / Omega)) - asinh(sqrt(OmegaLambda / Omega * a_ * a_ * a_)));
}


/** @brief For a given redshift, returns time to present for a given redshift
 *  
 * returns the time from that redshift until redshift zero:
 * \f$H_0t_0=\int_0^z\frac{dz}{(1+z)\sqrt{(1+z)^2(1+z\Omega_m)-z(2+x)\Omega_{\Lambda}}}\f$
 *
 * @return time to present (lookback time) in code units/Hubble_h
 */
double time_to_present(const double redshift_)
{
#ifdef ASSUME_FLAT_LCDM
  return time_to_present_flat_LCDM(redshift_);
#else /* not defined ASSUME_FLAT_LCDM */
  return time_to_present_ni(redshift_);
#endif /* not defined ASSUME_FLAT_LCDM */       
}
