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

/* Created on: April 10, 2012
 *         by: Bruno Henriques (bmh20)
 *
 * modified on 2018 Jan 15+ by Stefan Hilbert (hilbert)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

#include "model_dust_extinction_inline.h"

/** @file post_process_spec_mags.c
 *  @brief post_process_spec_mags.c can be used to compute mags or spectra from
 *        star formation histories. It also applies dust corrections.
 *
 *  When STAR_FORMATION_HISTORY option is ON in the Makefile the code
 *  stores the star formation history in the different components of
 *  the galaxy in log bins. If POST_PROCESS_MAGS is ON these star
 *  formation histories are used in this routine to compute magnitudes
 *  just before output.
 * */

/* standard #ifdef pattern for dealing with GALAXY_OUTPUT magnitudes: */
#ifdef OUTPUT_REST_MAGS
#ifdef ICL
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */                                                               
#ifdef OUTPUT_OBS_MAGS                                                                              
#ifdef ICL
#endif /* defined ICL */
#ifdef OUTPUT_FB_OBS_MAGS                                                                          
#ifdef ICL
#endif /* defined ICL */
#ifdef ICL
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */


/** @brief computes lumsinosities for galaxies */
static inline void 
//compute_post_process_luminosities (const double mass_, const double age_, const double metal_fraction_, const int snapshot_number_, const int filter_number_,
compute_post_process_luminosities (const double mass_,
                                   const int met_index_, const double f_met_1_, const double f_met_2_,
                                   const int age_index_, const double f_age_1_, const double f_age_2_,
                                   const int filter_number_, 
#ifdef OUTPUT_OBS_MAGS
                                   const int snapshot_number_, 
#endif /* defined OUTPUT_OBS_MAGS */
#ifdef OUTPUT_REST_MAGS
                                   double *Lum_            , double *YLum_            , 
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
                                   double *ObsLum_         , double *ObsYLum_         , 
#ifdef OUTPUT_FB_OBS_MAGS
                                   double *backward_ObsLum_        , double *backward_ObsYLum_        , 
                                   double *forward_ObsLum_, double *forward_ObsYLum_,
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
                                   const bool have_to_add_to_y_lum_)
{
  int    redshift_index_;
  double lum_to_add_;
  
#ifdef OUTPUT_REST_MAGS
  redshift_index_ = 0;
  lum_to_add_     = mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][redshift_index_][filter_number_]  +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_    ][redshift_index_][filter_number_]) +
                             f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][redshift_index_][filter_number_]  +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][redshift_index_][filter_number_])  );
  *Lum_ += lum_to_add_;
  if(have_to_add_to_y_lum_)
  { *YLum_ += lum_to_add_; }

#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
  redshift_index_ = LastDarkMatterSnapShot - snapshot_number_;
  lum_to_add_     = mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][redshift_index_][filter_number_] +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_    ][redshift_index_][filter_number_]) +
                             f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][redshift_index_][filter_number_] +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][redshift_index_][filter_number_]));
  *ObsLum_ += lum_to_add_;
  if(have_to_add_to_y_lum_)
  { *ObsYLum_ += lum_to_add_; }

#ifdef OUTPUT_FB_OBS_MAGS
  redshift_index_ = LastDarkMatterSnapShot - (snapshot_number_ > 0 ? snapshot_number_ - 1 : 0);
  lum_to_add_     = mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][redshift_index_][filter_number_] +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_    ][redshift_index_][filter_number_]) +
                             f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][redshift_index_][filter_number_] +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][redshift_index_][filter_number_]));
  *backward_ObsLum_ += lum_to_add_;
  if(have_to_add_to_y_lum_)
  { *backward_ObsYLum_ += lum_to_add_; }

  redshift_index_ = LastDarkMatterSnapShot - (snapshot_number_ < LastDarkMatterSnapShot ? snapshot_number_ + 1 : LastDarkMatterSnapShot);
  if(redshift_index_ < 0) redshift_index_ = 0;
  lum_to_add_     = mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][redshift_index_][filter_number_] +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_    ][redshift_index_][filter_number_]) +
                             f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][redshift_index_][filter_number_] +
                                         f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][redshift_index_][filter_number_]));
  *forward_ObsLum_ += lum_to_add_;
  if(have_to_add_to_y_lum_)
  { *forward_ObsYLum_ += lum_to_add_; }

#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
}


/** @brief computes dust corrections for lums. */
static inline void
make_dust_correction_for_disk_luminosities(const int filter_number_, const int snapshot_number_, const double Z_g_, const double cold_gas_, const double gas_disk_radius_, const double cos_inclination_ 
#ifdef OUTPUT_REST_MAGS
                                         , double *LumDisk_            , const double YLumDisk_
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
                                         , double *ObsLumDisk_         , const double ObsYLumDisk_
#ifdef OUTPUT_FB_OBS_MAGS
                                         , double *backward_ObsLumDisk_, const double backward_ObsYLumDisk_
                                         , double *forward_ObsLumDisk_ , const double forward_ObsYLumDisk_
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
                                          )
{
  double tau_, alam_;
//  double taubc_,
  double a_bc_;
  double mu_;
  
  const double FilterLambda_ratio_ = FilterLambda[filter_number_] * (1. / 0.55); /* filter wavelength / V band wavelength (reference) */

  /* 0.94 = 2.83/3. - 3 to get scale lenght and 2.83 = 1.68^2 */
  /* n_h_ = ColdGas_ / (M_PI * pow(GasDiskRadius_ * 0.56, 2) * 1.4); */
  /* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (2.1 10^21 atoms/cm^2) */
  /* note: 3252.37 = 10^(3.5122) */
  /* add  redshift dependence of dust-to-ColdGas_ ratio */
//  const double n_h_ = ColdGas_ / (M_PI * pow(GasDiskRadius_ * 0.94, 2) * 1.4 * 3252.37 * (1 + ZZ[snapshot_number_]));
  const double n_h_ = cold_gas_ / ((M_PI * 0.94 * 0.94 * 1.4 * 3252.37) * gas_disk_radius_ * gas_disk_radius_ * (1 + ZZ[snapshot_number_]));

 /* minimum inclination ~80 degrees, i.e. cos_inclination_ == 0.2 */
  const double n_h_sec_ = (cos_inclination_ < 0.2) ? (5. * n_h_) : (n_h_ / cos_inclination_);

  /* mu_ for YS extinction, given by a Gaussian with centre 0.3 (MUCENTER)
   * and width 0.2 (MUWIDTH), truncated at 0.1 and 1.  */
  do { mu_ = gsl_ran_gaussian(random_generator, MUWIDTH) + MUCENTER; }
  while (mu_ < 0.1 || mu_ > 1.0);
  
  // for testing:
  mu_ = MUCENTER;

  // extinction on Vband used as reference for the BC extinction
  const double tauvbc_ = get_extinction(NMAG, Z_g_, 0) * n_h_ * (1. / mu_ - 1.);

#ifdef OUTPUT_REST_MAGS
  //extinction due to ISM and light from young stars absorbed by birth clouds
//   taubc_ = tauvbc_ * pow(FilterLambda_ratio_, -0.7);
  a_bc_  = 1. - exp(-tauvbc_ * pow(FilterLambda_ratio_, -0.7));
  tau_   = get_extinction(filter_number_, Z_g_, 0) * n_h_sec_;
  alam_  = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;

  *LumDisk_ -= YLumDisk_ * a_bc_;
  *LumDisk_ *= alam_;
  
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
  //extinction due to ISM and light from young stars absorbed by birth clouds
//   taubc_ = tauvbc_ * pow((1. + ZZ[snapshot_number_]) * FilterLambda_ratio_, -0.7);
  a_bc_  = 1. - exp(-tauvbc_ * pow((1. + ZZ[snapshot_number_]) * FilterLambda_ratio_, -0.7));
  tau_   = get_extinction(filter_number_, Z_g_, ZZ[snapshot_number_]) * n_h_sec_;
  alam_  = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;

  *ObsLumDisk_ -= ObsYLumDisk_ * a_bc_;
  *ObsLumDisk_ *= alam_;

#ifdef OUTPUT_FB_OBS_MAGS   // compute same thing at z + 1
  const int earlier_snapshot_number_ = (snapshot_number_ > 0) ? (snapshot_number_ - 1) : 0;
//   taubc_ = tauvbc_ * pow((1. + ZZ[earlier_snapshot_number_]) * FilterLambda_ratio_, -0.7);
  a_bc_  = 1. - exp(-tauvbc_ * pow((1. + ZZ[earlier_snapshot_number_]) * FilterLambda_ratio_, -0.7));
  tau_   = get_extinction(filter_number_, Z_g_, ZZ[earlier_snapshot_number_]) * n_h_sec_;
  alam_  = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;

  *backward_ObsLumDisk_ -= backward_ObsYLumDisk_ * a_bc_;
  *backward_ObsLumDisk_ *= alam_;

  const int later_snapshot_number_ = (snapshot_number_ < LastDarkMatterSnapShot) ? (snapshot_number_ + 1) : LastDarkMatterSnapShot;
//   taubc_ = tauvbc_ * pow((1. + ZZ[later_snapshot_number_]) * FilterLambda_ratio_, -0.7);
  a_bc_  = 1. - exp(-tauvbc_ * pow((1. + ZZ[later_snapshot_number_]) * FilterLambda_ratio_, -0.7));
  tau_   = get_extinction(filter_number_, Z_g_, ZZ[later_snapshot_number_]) * n_h_sec_;
  alam_  = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;

  *forward_ObsLumDisk_ -= forward_ObsYLumDisk_ * a_bc_;
  *forward_ObsLumDisk_ *= alam_;
  
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
}


/** @brief computes lums./mags. for galaxies */
void post_process_spec_mags(struct GALAXY_OUTPUT *galaxy_)
{
  /** if one wants to have finner bins for the star formation then the STEPS
    * of the calculation, N_FINE_AGE_BINS should be set to > 1 */
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
  int fine_age_step_;
#endif /* defined N_FINE_AGE_BINS > 1 */

  /* Time below which the luminosities are corrected for extinction due to
   * molecular birth clouds.  */
  const double birthcloud_age_ = 10.0 / UnitTime_in_Megayears * Hubble_h;

  /* used for dust corrections: */
  const double Z_g_        = 50. * metals_total(galaxy_->MetalsColdGas) / galaxy_->ColdGas;
  const double t_snapshot_ = NumToTime(galaxy_->SnapNum);

  /* used for interpolating luminosities: */
#if !((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
  /* static int snapshot_number_of_previous_galaxy_ = -1; */
#endif /* not defined N_FINE_AGE_BINS > 1 */  
  /* static */ int       age_index_[SFH_NBIN]; /* static */ double       f_age_1_[SFH_NBIN]; /* static */ double       f_age_2_[SFH_NBIN];
         int  disk_met_index_[SFH_NBIN];        double  disk_f_met_1_[SFH_NBIN];        double  disk_f_met_2_[SFH_NBIN];
         int bulge_met_index_[SFH_NBIN];        double bulge_f_met_1_[SFH_NBIN];        double bulge_f_met_2_[SFH_NBIN];
#ifdef ICL
         int   icm_met_index_[SFH_NBIN];        double   icm_f_met_1_[SFH_NBIN];        double   icm_f_met_2_[SFH_NBIN];
#endif /* defined ICL */

  int sfh_bin_number_, filter_number_;
  /* precompute parameters for interpolating luminosities: */
#if !((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
  /* if(snapshot_number_of_previous_galaxy_ != galaxy_->SnapNum)  */
  {
    /* snapshot_number_of_previous_galaxy_ = galaxy_->SnapNum; */
    for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
    {
      const double age_       = SFH_t[galaxy_->SnapNum][0][sfh_bin_number_] + 0.5 * SFH_dt[galaxy_->SnapNum][0][sfh_bin_number_] - t_snapshot_;
      
      if(age_ <= 0)
      {
        age_index_[sfh_bin_number_] = 0; f_age_1_[sfh_bin_number_] = 0.; f_age_2_[sfh_bin_number_] = 0.;
      } 
      else
      {
        const double log10_age_ = log10(age_);
        find_age_luminosity_interpolation_parameters(log10_age_ , age_index_[sfh_bin_number_], f_age_1_[sfh_bin_number_], f_age_2_[sfh_bin_number_]);
      }
    }
  }
#endif /* not defined N_FINE_AGE_BINS > 1 */
  
  /** @bug (fixed by Stefan Hilbert) before, for MetallicityOption == 0 (only solar mets.) just the met. indices be set, but not the f_met_ factors,
  *       now also  f_met_ factors are set.
  */
  if(MetallicityOption == 0)
  {
    for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
    {
      disk_met_index_ [sfh_bin_number_] = 4;  disk_f_met_1_ [sfh_bin_number_] = 1.;  disk_f_met_2_ [sfh_bin_number_] = 0.;
      bulge_met_index_[sfh_bin_number_] = 4; bulge_f_met_1_ [sfh_bin_number_] = 1.; bulge_f_met_2_ [sfh_bin_number_] = 0.;
#ifdef ICL
      icm_met_index_  [sfh_bin_number_] = 4;   icm_f_met_1_ [sfh_bin_number_] = 1.;   icm_f_met_2_ [sfh_bin_number_] = 0.;
#endif /* defined ICL */
    }
  }
  else
  {  
    for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
    {
      if(galaxy_->sfh_DiskMass [sfh_bin_number_] > 0.)
      {
        const double disk_metal_fraction_  = metals_total(galaxy_->sfh_MetalsDiskMass [sfh_bin_number_]) / galaxy_->sfh_DiskMass [sfh_bin_number_];
        if(disk_metal_fraction_ <= 0.)
        {
          disk_met_index_ [sfh_bin_number_] = 0.; disk_f_met_1_ [sfh_bin_number_] = 1.;  disk_f_met_2_ [sfh_bin_number_] = 0.;
        }
        else
        {
          const double log10_disk_metal_fraction_  = log10(disk_metal_fraction_);
          find_metallicity_luminosity_interpolation_parameters(log10_disk_metal_fraction_ , disk_met_index_ [sfh_bin_number_], disk_f_met_1_ [sfh_bin_number_], disk_f_met_2_ [sfh_bin_number_]);
        }
      }
      if(galaxy_->sfh_BulgeMass[sfh_bin_number_] > 0.)
      {
        const double bulge_metal_fraction_ = metals_total(galaxy_->sfh_MetalsBulgeMass[sfh_bin_number_]) / galaxy_->sfh_BulgeMass[sfh_bin_number_];
        if(bulge_metal_fraction_ <= 0.)
        {
          bulge_met_index_ [sfh_bin_number_] = 0.; bulge_f_met_1_ [sfh_bin_number_] = 1.;  bulge_f_met_2_ [sfh_bin_number_] = 0.;
        }
        else
        {
          const double log10_bulge_metal_fraction_ = log10(bulge_metal_fraction_);
          find_metallicity_luminosity_interpolation_parameters(log10_bulge_metal_fraction_, bulge_met_index_[sfh_bin_number_], bulge_f_met_1_[sfh_bin_number_], bulge_f_met_2_[sfh_bin_number_]);
        }
      }
#ifdef ICL
      if(galaxy_->sfh_ICM      [sfh_bin_number_] > 0.)
      {
        const double icm_metal_fraction_   = metals_total(galaxy_->sfh_MetalsICM      [sfh_bin_number_]) / galaxy_->sfh_ICM      [sfh_bin_number_];
        if(icm_metal_fraction_ <= 0.)
        {
          icm_met_index_ [sfh_bin_number_] = 0.; icm_f_met_1_ [sfh_bin_number_] = 1.;  icm_f_met_2_ [sfh_bin_number_] = 0.;
        }
        else
        {
          const double log10_icm_metal_fraction_   = log10(icm_metal_fraction_);
          find_metallicity_luminosity_interpolation_parameters(log10_icm_metal_fraction_  , icm_met_index_  [sfh_bin_number_], icm_f_met_1_  [sfh_bin_number_], icm_f_met_2_  [sfh_bin_number_]);
        }
      }
#endif /* defined ICL */
    }
  }

  /** r-band weighted ages
   * 
   *  @bug (fixed by Stefan Hilbert)
   *       if not defined OUTPUT_REST_MAGS,
   *       but 0 < (galaxy_->DiskMass + galaxy_->BulgeMass),
   *       galaxy_->rbandWeightAge would compute to 0./0. (i.e. NaN).
   *       nowg alaxy_->rbandWeightAg is set to 0 and stays 0 
   *       if not fefined OUTPUT_REST_MAG
   * 
   * @bug  filter_number_==17 is not necessarily r-band,
   *       since filters are assigned from info in parameter file
   *
   * @todo parametrize filter band number for computing weighted stellar age
   */
#ifdef OUTPUT_REST_MAGS
  const int r_band_filter_number_ = 17;
#endif /* defined OUTPUT_REST_MAGS */
  galaxy_->rbandWeightAge = 0.;

  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  {
#ifdef OUTPUT_REST_MAGS
    double LumDisk_  = 0.; double YLumDisk_  = 0.;
    double LumBulge_ = 0.; double YLumBulge_ = 0.;
#ifdef ICL
    double LumICL_   = 0.; double YLumICL_   = 0.;
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS                                                                              
    double ObsLumDisk_  = 0.; double ObsYLumDisk_  = 0.;
    double ObsLumBulge_ = 0.; double ObsYLumBulge_ = 0.;
#ifdef ICL
    double ObsLumICL_   = 0.; double ObsYLumICL_   = 0.;
#endif /* defined ICL */
#ifdef OUTPUT_FB_OBS_MAGS                                                                          
    double backward_ObsLumDisk_  = 0.; double backward_ObsYLumDisk_  = 0.;
    double backward_ObsLumBulge_ = 0.; double backward_ObsYLumBulge_ = 0.;
#ifdef ICL
    double backward_ObsLumICL_   = 0.; double backward_ObsYLumICL_   = 0.;
#endif /* defined ICL */
    double forward_ObsLumDisk_  = 0.; double forward_ObsYLumDisk_  = 0.;
    double forward_ObsLumBulge_ = 0.; double forward_ObsYLumBulge_ = 0.;
#ifdef ICL
    double forward_ObsLumICL_   = 0.; double forward_ObsYLumICL_   = 0.;
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */

#ifdef OUTPUT_REST_MAGS
    double previous_r_band_luminosity_ = 0.;
#endif /* defined OUTPUT_REST_MAGS */                                                               

    //loop on SFH bins
    for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
    {
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
      for(fine_age_step_ = 0; fine_age_step_< N_FINE_AGE_BINS; fine_age_step_++)
      {
        const double age_ = SFH_t[galaxy_->SnapNum][0][sfh_bin_number_] + (0.5 + fine_age_step_) * SFH_dt[galaxy_->SnapNum][0][sfh_bin_number_] / N_FINE_AGE_BINS - t_snapshot_;
#else  /* not defined N_FINE_AGE_BINS > 1 */
        const double age_ = SFH_t[galaxy_->SnapNum][0][sfh_bin_number_] + 0.5 * SFH_dt[galaxy_->SnapNum][0][sfh_bin_number_] - t_snapshot_;
#endif /* not defined N_FINE_AGE_BINS > 1 */
      
        if(age_ <= 0.) continue;
        
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
        const double log10_age_ = log10(age_);
        find_age_luminosity_interpolation_parameters(log10_age_ , age_index_[sfh_bin_number_], f_age_1_[sfh_bin_number_], f_age_2_[sfh_bin_number_]);
#endif /* defined N_FINE_AGE_BINS > 1 */

        const bool   is_affected_by_birthclould_ = (age_ <= birthcloud_age_);

        /* The stellar populations tables have magnitudes for all the mass
         * formed in stars including what will be shortly lost by SNII   */
#ifdef DETAILED_METALS_AND_MASS_RETURN
        const double disk_mass_  = galaxy_->sfh_DiskMass [sfh_bin_number_] * inv_Hubble_h * 0.1;
        const double bulge_mass_ = galaxy_->sfh_BulgeMass[sfh_bin_number_] * inv_Hubble_h * 0.1;
#ifdef ICL
        const double icm_mass_   = galaxy_->sfh_ICM      [sfh_bin_number_] * inv_Hubble_h * 0.1;
#endif /* defined ICL */
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
        const double disk_mass_  = galaxy_->sfh_DiskMass [sfh_bin_number_] * inv_Hubble_h * 0.1 / (1 - RecycleFraction);
        const double bulge_mass_ = galaxy_->sfh_BulgeMass[sfh_bin_number_] * inv_Hubble_h * 0.1 / (1 - RecycleFraction);
#ifdef ICL              
        const double icm_mass_   = galaxy_->sfh_ICM      [sfh_bin_number_] * inv_Hubble_h * 0.1 / (1 - RecycleFraction);
#endif /* defined ICL */
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */

        /* disk luminosities */
        if(galaxy_->sfh_DiskMass[sfh_bin_number_] > 0.)
          compute_post_process_luminosities (disk_mass_,
                                            disk_met_index_[sfh_bin_number_], disk_f_met_1_[sfh_bin_number_], disk_f_met_2_[sfh_bin_number_], 
                                            age_index_     [sfh_bin_number_], f_age_1_     [sfh_bin_number_], f_age_2_     [sfh_bin_number_], 
                                            filter_number_,
#ifdef OUTPUT_OBS_MAGS
                                            galaxy_->SnapNum, 
#endif /* defined OUTPUT_OBS_MAGS */
#ifdef OUTPUT_REST_MAGS
                                            &LumDisk_            , &YLumDisk_            , 
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
                                            &ObsLumDisk_         , &ObsYLumDisk_         , 
#ifdef OUTPUT_FB_OBS_MAGS
                                            &backward_ObsLumDisk_        , &backward_ObsYLumDisk_        , 
                                            &forward_ObsLumDisk_, &forward_ObsYLumDisk_,
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
                                            is_affected_by_birthclould_);

        /* bulge luminosities */
        if(galaxy_->sfh_BulgeMass[sfh_bin_number_] > 0.)
          compute_post_process_luminosities (bulge_mass_,  
                                            bulge_met_index_[sfh_bin_number_], bulge_f_met_1_[sfh_bin_number_], bulge_f_met_2_[sfh_bin_number_], 
                                            age_index_     [sfh_bin_number_], f_age_1_     [sfh_bin_number_], f_age_2_     [sfh_bin_number_], 
                                            filter_number_,
#ifdef OUTPUT_OBS_MAGS
                                            galaxy_->SnapNum, 
#endif /* defined OUTPUT_OBS_MAGS */
#ifdef OUTPUT_REST_MAGS
                                            &LumBulge_            , &YLumBulge_            , 
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
                                            &ObsLumBulge_         , &ObsYLumBulge_         , 
#ifdef OUTPUT_FB_OBS_MAGS
                                            &backward_ObsLumBulge_        , &backward_ObsYLumBulge_        , 
                                            &forward_ObsLumBulge_, &forward_ObsYLumBulge_,
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
                                            is_affected_by_birthclould_);

#ifdef ICL
        /** @note ICL does not seem to get correction  by dust or birth clouds */
        /* ICL luminosities */
        if(galaxy_->sfh_ICM[sfh_bin_number_] > 0.)
          compute_post_process_luminosities (icm_mass_,  
                                            icm_met_index_ [sfh_bin_number_], icm_f_met_1_ [sfh_bin_number_], icm_f_met_2_ [sfh_bin_number_], 
                                            age_index_     [sfh_bin_number_], f_age_1_     [sfh_bin_number_], f_age_2_     [sfh_bin_number_], 
                                            filter_number_,
#ifdef OUTPUT_OBS_MAGS
                                            galaxy_->SnapNum, 
#endif /* defined OUTPUT_OBS_MAGS */
#ifdef OUTPUT_REST_MAGS
                                            &LumICL_            , &YLumICL_            , 
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
                                            &ObsLumICL_         , &ObsYLumICL_         , 
#ifdef OUTPUT_FB_OBS_MAGS
                                            &backward_ObsLumICL_        , &backward_ObsYLumICL_        , 
                                            &forward_ObsLumICL_, &forward_ObsYLumICL_,
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
                                            false);
#endif /* defined ICL */

#ifndef LIGHT_OUTPUT
#ifdef OUTPUT_REST_MAGS
        if((galaxy_->DiskMass+galaxy_->BulgeMass) > 0. && filter_number_ == r_band_filter_number_)
        {
          galaxy_->rbandWeightAge += age_ * (LumDisk_ + LumBulge_ - previous_r_band_luminosity_);
          previous_r_band_luminosity_= LumDisk_ + LumBulge_;
        }
#endif /* defined OUTPUT_REST_MAGS */
#endif /* not defined LIGHT_OUTPUT */


#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
      } // end of loop on fine bins
#endif /* defined N_FINE_AGE_BINS > 1 */
    }//end of loop on sfh bins

#ifndef LIGHT_OUTPUT
#ifdef OUTPUT_REST_MAGS
    if((galaxy_->DiskMass+galaxy_->BulgeMass) > 0. && filter_number_ == r_band_filter_number_)
    {
      galaxy_->rbandWeightAge /= (LumDisk_ + LumBulge_);
      galaxy_->rbandWeightAge = galaxy_->rbandWeightAge * 1.e-3 * UnitTime_in_Megayears * inv_Hubble_h; //conversion in age_ from code units/h -> Gyr
    }
#endif /* defined OUTPUT_REST_MAGS */

#ifndef LIGHT_OUTPUT
#ifdef OUTPUT_REST_MAGS
    galaxy_->Mag                 [filter_number_] = lum_to_mag(LumDisk_ + LumBulge_);
    galaxy_->MagBulge            [filter_number_] = lum_to_mag(LumBulge_);
#ifdef ICL
    galaxy_->MagICL              [filter_number_] = lum_to_mag(LumICL_);
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
    galaxy_->ObsMag              [filter_number_] = lum_to_mag(ObsLumDisk_ + ObsLumBulge_);
    galaxy_->ObsMagBulge         [filter_number_] = lum_to_mag(ObsLumBulge_);
#ifdef ICL
    galaxy_->ObsMagICL           [filter_number_] = lum_to_mag(ObsLumICL_);
#endif /* defined ICL */

#ifdef OUTPUT_FB_OBS_MAGS
    galaxy_->backward_ObsMag     [filter_number_] = lum_to_mag(backward_ObsLumDisk_ + backward_ObsLumBulge_);
    galaxy_->backward_ObsMagBulge[filter_number_] = lum_to_mag(backward_ObsLumBulge_);
#ifdef ICL
    galaxy_->backward_ObsMagICL  [filter_number_] = lum_to_mag(backward_ObsLumICL_);
#endif /* defined ICL */
    galaxy_->forward_ObsMag      [filter_number_] = lum_to_mag(forward_ObsLumDisk_ + forward_ObsLumBulge_);
    galaxy_->forward_ObsMagBulge [filter_number_] = lum_to_mag(forward_ObsLumBulge_);
#ifdef ICL                       
    galaxy_->forward_ObsMagICL   [filter_number_] = lum_to_mag(forward_ObsLumICL_);
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
#endif /* not defined LIGHT_OUTPUT */

    /* inclination is needed for disk dust correction,
     * but not sure, inclination shouldn't already be computed elsewhere */
    galaxy_->CosInclination = fabs(galaxy_->StellarSpin[2]) /
        sqrt(galaxy_->StellarSpin[0] * galaxy_->StellarSpin[0] +
             galaxy_->StellarSpin[1] * galaxy_->StellarSpin[1] +
             galaxy_->StellarSpin[2] * galaxy_->StellarSpin[2]);

    /* dust correction for disk (remove light from dust and young stars absorbed by birth clouds) */
    if(galaxy_->ColdGas > 0.0)
      make_dust_correction_for_disk_luminosities(filter_number_, galaxy_->SnapNum, Z_g_, galaxy_->ColdGas, galaxy_->GasDiskRadius, galaxy_->CosInclination
#ifdef OUTPUT_REST_MAGS
                                               , &LumDisk_            , YLumDisk_
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
                                               , &ObsLumDisk_         , ObsYLumDisk_
#ifdef OUTPUT_FB_OBS_MAGS
                                               , &backward_ObsLumDisk_        , backward_ObsYLumDisk_
                                               , &forward_ObsLumDisk_, forward_ObsYLumDisk_
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
                                                );

    /* dust correction for bulges (remove light from young stars absorbed by birth clouds) */
#ifdef OUTPUT_REST_MAGS
    LumBulge_             -= YLumBulge_             * (1. - ExpTauBCBulge);
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
    ObsLumBulge_          -= ObsYLumBulge_          * (1. - ExpTauBCBulge);
#ifdef OUTPUT_FB_OBS_MAGS
    backward_ObsLumBulge_         -= backward_ObsYLumBulge_         * (1. - ExpTauBCBulge);
    forward_ObsLumBulge_ -= forward_ObsYLumBulge_ * (1. - ExpTauBCBulge);
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */

#ifdef OUTPUT_REST_MAGS
    galaxy_->MagDust            [filter_number_] = lum_to_mag(LumDisk_+LumBulge_);
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
    galaxy_->ObsMagDust         [filter_number_] = lum_to_mag(ObsLumDisk_+ObsLumBulge_);
#ifdef OUTPUT_FB_OBS_MAGS
    galaxy_->backward_ObsMagDust        [filter_number_] = lum_to_mag(backward_ObsLumDisk_+backward_ObsLumBulge_);
    galaxy_->forward_ObsMagDust[filter_number_] = lum_to_mag(forward_ObsLumDisk_+forward_ObsLumBulge_);
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
  }//end of loop on bands
}
