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

/* Created on 2018 Jan 15+ by Stefan Hilbert (hilbert)
 * orignially based on version by Bruno Henriques (bmh20)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

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


// standard #ifdef pattern for dealing with GALAXY_OUTPUT magnitudes:
#ifdef OUTPUT_REST_MAGS
#ifdef ICL
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */                                                               
#ifdef OUTPUT_OBS_MAGS                                                                              
#ifdef ICL
#endif /* defined ICL */
#ifdef OUTPUT_MOMAF_INPUTS                                                                          
#ifdef ICL
#endif /* defined ICL */
#ifdef KITZBICHLER                                                                                  
#ifdef ICL
#endif /* defined ICL */
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined OUTPUT_OBS_MAGS */


/** @brief Converts luminosities into magnitudes
 *
 * Converts luminosities into magnitudes:
 * \f$ M=-2.5\mathrm{log}_{10}(L) \f$ */
static inline double 
lum_to_lum_or_mag(const double lum_)
{
#ifdef FULL_SPECTRA  
  return lum_;
#else  /* not defined FULL_SPECTRA */
  if(lum_ > 0)
    return -2.5 * log10(lum_);
  else
    return 99.0;
#endif /* not defined FULL_SPECTRA */
}


/** @brief computes lumsinosities for galaxies */
static inline double
interpolated_luminosity(const double mass_, const int filter_number_, 
                        const int met_index_, const double f_met_1_, const double f_met_2_,
                        const int age_index_, const double f_age_1_, const double f_age_2_, 
                        const int z_index_)
{
  return mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][z_index_][filter_number_]  +
                              f_age_2_ * LumTables[age_index_ + 1][met_index_    ][z_index_][filter_number_]) +
                  f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][z_index_][filter_number_]  +
                              f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][z_index_][filter_number_])   );
}


/** @brief computes lumsinosities for galaxies and adds them to fields */
static inline void
add_interpolated_luminosities(const double mass_,
                              const int met_index_, const double f_met_1_, const double f_met_2_,
                              const int age_index_, const double f_age_1_, const double f_age_2_, 
                              const int z_index_,
                              double (*lum_) [NMAG], double (*y_lum_) [NMAG],
                              const bool have_to_add_to_y_lum_
                             )
{
  int filter_number_;
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  { 
    const double lum_to_add_ = interpolated_luminosity(mass_, filter_number_, 
                                                      met_index_, f_met_1_, f_met_2_,
                                                      age_index_, f_age_1_, f_age_2_,
                                                      z_index_);
    (*lum_) [filter_number_] += lum_to_add_;
    if(have_to_add_to_y_lum_)
    { (*y_lum_) [filter_number_] += lum_to_add_; }
  }
}


/** @brief computes dust corrections for lums. */
static inline void
make_dust_correction_for_post_processing(const int snapshot_number_, const double Z_g_, const double cold_gas_, const double gas_disk_radius_, const double cos_inclination_
#ifdef OUTPUT_REST_MAGS
                                       , double (*LumDisk_            )[NMAG], const double YLumDisk_            [NMAG]
#endif /* defined OUTPUT_REST_MAGS */
#ifdef COMPUTE_OBS_MAGS
                                       , double (*ObsLumDisk_         )[NMAG], const double ObsYLumDisk_         [NMAG]
#ifdef OUTPUT_MOMAF_INPUTS
                                       , double (*dObsLumDisk_        )[NMAG], const double dObsYLumDisk_        [NMAG]
#ifdef KITZBICHLER
                                       , double (*dObsLumDisk_forward_)[NMAG], const double dObsYLumDisk_forward_[NMAG]
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined COMPUTE_OBS_MAGS */
                                        )
{
  int filter_number_;
  double tau_, alam_, a_bc_, mu_;
  
  const double inv_VBand_WaveLength_ = 1./ 0.55;

  /* 0.94 = 2.83/3. - 3 to get scale lenght and 2.83 = 1.68^2 */
  /* n_h_ = cold_gas_ / (M_PI * pow(gas_disk_radius_ * 0.56, 2) * 1.4); */
  /* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (2.1 10^21 atoms/cm^2) */
  /* note: 3252.37 = 10^(3.5122) */
  /* add redshift dependence of dust-to-cold_gas ratio */
  //  const double n_h_ = ColdGas_ / (M_PI * pow(GasDiskRadius_ * 0.94, 2) * 1.4 * 3252.37 * (1 + ZZ[snapshot_number_]));
  const double n_h_ = cold_gas_ / ((M_PI * 0.94 * 0.94 * 1.4 * 3252.37) * gas_disk_radius_ * gas_disk_radius_ * (1 + ZZ[snapshot_number_]));

 // minimum inclination ~80 degrees, i.e. cos(inclination) == 0.2
  const double n_h_sec_ = (cos_inclination_ < 0.2) ? (5. * n_h_) : (n_h_ / cos_inclination_);

  /* mu_ for YS extinction, given by a Gaussian with centre 0.3 (MUCENTER)
   * and width 0.2 (MUWIDTH), truncated at 0.1 and 1.  */
  do { mu_ = gsl_ran_gaussian(random_generator, MUWIDTH) + MUCENTER; }
  while (mu_ < 0.1 || mu_ > 1.0);
  
//  // for testing:
   mu_ = MUCENTER;

  // extinction on Vband used as reference for the BC extinction
  const double tauvbc_ = get_extinction(NMAG, Z_g_, 0) * n_h_ * (1. / mu_ - 1.);
  
  //extinction due to ISM and light from young stars absorbed by birth clouds
#ifdef OUTPUT_REST_MAGS
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  {
    tau_  = get_extinction(filter_number_, Z_g_, 0) * n_h_sec_;
    alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
    a_bc_ = 1. - exp(-tauvbc_ * pow(FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

    (*LumDisk_)[filter_number_] -= YLumDisk_[filter_number_] * a_bc_;
    (*LumDisk_)[filter_number_] *= alam_;
  }
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  {
    tau_  = get_extinction(filter_number_, Z_g_, ZZ[snapshot_number_]) * n_h_sec_;
    alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
    a_bc_ = 1. - exp(-tauvbc_ * pow((1. + ZZ[snapshot_number_]) * FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

    (*ObsLumDisk_)[filter_number_] -= ObsYLumDisk_[filter_number_] * a_bc_;
    (*ObsLumDisk_)[filter_number_] *= alam_;
  }
#ifdef OUTPUT_MOMAF_INPUTS   // compute same thing at z + 1
  const int earlier_snapshot_number_ = (snapshot_number_ > 0) ? (snapshot_number_ - 1) : 0;
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  {
    tau_  = get_extinction(filter_number_, Z_g_, ZZ[earlier_snapshot_number_]) * n_h_sec_;
    alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
    a_bc_ = 1. - exp(-tauvbc_ * pow((1. + ZZ[earlier_snapshot_number_]) * FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

    (*dObsLumDisk_)[filter_number_] -= dObsYLumDisk_[filter_number_] * a_bc_;
    (*dObsLumDisk_)[filter_number_] *= alam_;
  }
#ifdef KITZBICHLER
  const int later_snapshot_number_ = (snapshot_number_ < LastDarkMatterSnapShot) ? (snapshot_number_ + 1) : LastDarkMatterSnapShot;
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  {
    tau_  = get_extinction(filter_number_, Z_g_, ZZ[later_snapshot_number_]) * n_h_sec_;
    alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
    a_bc_ = 1. - exp(-tauvbc_ * pow((1. + ZZ[later_snapshot_number_]) * FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

    (*dObsLumDisk_forward_)[filter_number_] -= dObsYLumDisk_forward_[filter_number_] * a_bc_;
    (*dObsLumDisk_forward_)[filter_number_] *= alam_;
  }
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined OUTPUT_OBS_MAGS */
}


/** @brief computes lums./mags. for galaxies */
void
post_process_spec_mags(struct GALAXY_OUTPUT *galaxy_)
{
   /* Time below which the luminosities are corrected for extinction due to
   * molecular birth clouds.  */
  const double birthcloud_age_ = 10.0 / UnitTime_in_Megayears * Hubble_h;

  /* metallicity parameter for dust corrections: */
  const double Z_g_ = metals_total(galaxy_->MetalsColdGas)/galaxy_->ColdGas/0.02;
  
  /* for computing stellar ages: */
  const double t_snapshot_ = NumToTime(galaxy_->SnapNum);
  
  /** if one wants to have finner bins for the star formation then the STEPS
    * of the calculation, N_FINE_AGE_BINS should be set to > 1 */
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
  int fine_age_step_;
#endif /* defined N_FINE_AGE_BINS > 1 */

#ifdef OUTPUT_REST_MAGS
  /** @bug  filter_number_ = 17 is not necessarily r-band, since filters are assigned from info in parameter file
   *
   *  @todo  parametrize filter band number for computing weighted stellar age_
   */
  const int r_band_filter_number_ = 17;
  double previous_r_band_luminosity_ = 0.;
#endif /* defined OUTPUT_REST_MAGS */                                                               
 
  galaxy_->rbandWeightAge = 0.;
 
  /* for accumulating luminosities, 
   * incl. that of young stellar populations affected by molecular birth cloud.
   * note: consider going to float arrays and possibly using galaxy_'s fields */
#ifdef OUTPUT_REST_MAGS
  double LumDisk_               [NMAG]; set_array_to(LumDisk_              , 0.);
  double YLumDisk_              [NMAG]; set_array_to(YLumDisk_             , 0.);
  double LumBulge_              [NMAG]; set_array_to(LumBulge_             , 0.);
  double YLumBulge_             [NMAG]; set_array_to(YLumBulge_            , 0.);
#ifdef ICL
  double LumICL_                [NMAG]; set_array_to(LumICL_               , 0.);
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */                                                               
#ifdef OUTPUT_OBS_MAGS                                                                              
  double ObsLumDisk_            [NMAG]; set_array_to(ObsLumDisk_           , 0.);
  double ObsYLumDisk_           [NMAG]; set_array_to(ObsYLumDisk_          , 0.);
  double ObsLumBulge_           [NMAG]; set_array_to(ObsLumBulge_          , 0.);
  double ObsYLumBulge_          [NMAG]; set_array_to(ObsYLumBulge_         , 0.);
#ifdef ICL
  double ObsLumICL_             [NMAG]; set_array_to(ObsLumICL_            , 0.);
#endif /* defined ICL */
#ifdef OUTPUT_MOMAF_INPUTS                                                                          
  double dObsLumDisk_           [NMAG]; set_array_to(dObsLumDisk_          , 0.);
  double dObsYLumDisk_          [NMAG]; set_array_to(dObsYLumDisk_         , 0.);
  double dObsLumBulge_          [NMAG]; set_array_to(dObsLumBulge_         , 0.);
  double dObsYLumBulge_         [NMAG]; set_array_to(dObsYLumBulge_        , 0.);
#ifdef ICL
  double dObsLumICL_            [NMAG]; set_array_to(dObsLumICL_           , 0.);
#endif /* defined ICL */
#ifdef KITZBICHLER                                                                                  
  double dObsLumDisk_forward_   [NMAG]; set_array_to(dObsLumDisk_forward_  , 0.);
  double dObsYLumDisk_forward_  [NMAG]; set_array_to(dObsYLumDisk_forward_ , 0.);
  double dObsLumBulge_forward_  [NMAG]; set_array_to(dObsLumBulge_forward_ , 0.);
  double dObsYLumBulge_forward_ [NMAG]; set_array_to(dObsYLumBulge_forward_, 0.);
#ifdef ICL
  double dObsLumICL_forward_    [NMAG]; set_array_to(dObsLumICL_forward_   , 0.);
#endif /* defined ICL */
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined OUTPUT_OBS_MAGS */

  /* ICL does not seem to get corrected for m.b.c. yet,
   * but we use a dummy field to meet function parameter signature */
#ifdef ICL
  double YLumICL_dummy_ [NMAG];
#endif /* defined ICL */

  /* indices and fractions for interpolating luminosities: */
  int age_index_      ; double f_age_1_      ; double f_age_2_      ;
  int disk_met_index_ ; double disk_f_met_1_ ; double disk_f_met_2_ ;
  int bulge_met_index_; double bulge_f_met_1_; double bulge_f_met_2_;
#ifdef ICL  
  int icm_met_index_  ; double icm_f_met_1_  ; double icm_f_met_2_  ;
 #endif /* defined ICL */

#ifdef OUTPUT_REST_MAGS
  const int redshift_index_              = 0;
#endif /* defined OUTPUT_REST_MAGS */                                                               
#ifdef OUTPUT_OBS_MAGS  
  const int redshift_index_obs_          = LastDarkMatterSnapShot -  galaxy_->SnapNum;
#ifdef OUTPUT_MOMAF_INPUTS                                                                          
  const int redshift_index_dobs_         = LastDarkMatterSnapShot - (galaxy_->SnapNum > 0                      ? galaxy_->SnapNum - 1 : 0                     );
#ifdef KITZBICHLER                                                                                  
  const int redshift_index_dobs_forward_ = LastDarkMatterSnapShot - (galaxy_->SnapNum < LastDarkMatterSnapShot ? galaxy_->SnapNum + 1 : LastDarkMatterSnapShot);
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined OUTPUT_OBS_MAGS */

  int sfh_bin_number_;
  for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
  {
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
    for(fine_age_step_ = 0; fine_age_step_< N_FINE_AGE_BINS; fine_age_step_++)
    {
      const double age_ = SFH_t[galaxy_->SnapNum][0][sfh_bin_number_] + (0.5 + fine_age_step_) * SFH_dt[galaxy_->SnapNum][0][sfh_bin_number_] / N_FINE_AGE_BINS - t_snapshot_;

#else  /* not defined N_FINE_AGE_BINS > 1 */
      const double age_ = SFH_t[galaxy_->SnapNum][0][sfh_bin_number_] + 0.5 * SFH_dt[galaxy_->SnapNum][0][sfh_bin_number_] - t_snapshot_;
#endif /* not defined N_FINE_AGE_BINS > 1 */
    
      if(age_ <= 0.)
        continue;
      
      const double log10_age_                  = log10(age_);
      const bool   is_affected_by_birthclould_ = (age_ <= birthcloud_age_);
  
      find_age_luminosity_interpolation_parameters(log10_age_, age_index_, f_age_1_, f_age_2_);

      /* The stellar populations tables have magnitudes for all the mass
      * formed in stars including what will be shortly lost by SNII   */
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
#ifdef DETAILED_METALS_AND_MASS_RETURN
      const double disk_mass_  = galaxy_->sfh_DiskMass [sfh_bin_number_] * 0.1 / N_FINE_AGE_BINS / Hubble_h;
      const double bulge_mass_ = galaxy_->sfh_BulgeMass[sfh_bin_number_] * 0.1 / N_FINE_AGE_BINS / Hubble_h;
#ifdef ICL
      const double icm_mass_   = galaxy_->sfh_ICM      [sfh_bin_number_] * 0.1 / N_FINE_AGE_BINS / Hubble_h;
#endif /* defined ICL */
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */        
      const double disk_mass_  = galaxy_->sfh_DiskMass [sfh_bin_number_] * 0.1 / N_FINE_AGE_BINS / (Hubble_h * (1 - RecycleFraction) );
      const double bulge_mass_ = galaxy_->sfh_BulgeMass[sfh_bin_number_] * 0.1 / N_FINE_AGE_BINS / (Hubble_h * (1 - RecycleFraction) );
#ifdef ICL
      const double icm_mass_   = galaxy_->sfh_ICM      [sfh_bin_number_] * 0.1 / N_FINE_AGE_BINS / (Hubble_h * (1 - RecycleFraction) );
#endif /* defined ICL */
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */
#else  /* not defined N_FINE_AGE_BINS > 1 */
#ifdef DETAILED_METALS_AND_MASS_RETURN
      const double disk_mass_  = galaxy_->sfh_DiskMass [sfh_bin_number_] * 0.1 / Hubble_h;
      const double bulge_mass_ = galaxy_->sfh_BulgeMass[sfh_bin_number_] * 0.1 / Hubble_h;
#ifdef ICL
      const double icm_mass_   = galaxy_->sfh_ICM      [sfh_bin_number_] * 0.1 / Hubble_h;
#endif /* defined ICL */
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */        
      const double disk_mass_  = galaxy_->sfh_DiskMass [sfh_bin_number_] * 0.1 / (Hubble_h * (1 - RecycleFraction) );
      const double bulge_mass_ = galaxy_->sfh_BulgeMass[sfh_bin_number_] * 0.1 / (Hubble_h * (1 - RecycleFraction) );
#ifdef ICL
      const double icm_mass_   = galaxy_->sfh_ICM      [sfh_bin_number_] * 0.1 / (Hubble_h * (1 - RecycleFraction) );
#endif /* defined ICL */
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */
#endif /* not defined N_FINE_AGE_BINS > 1 */

      if(disk_mass_ > 0.)
      {
        if(MetallicityOption == 0)
        { disk_met_index_  = 4;  disk_f_met_1_ = 1.;  disk_f_met_2_ = 0.; }
        else if(metals_total(galaxy_->sfh_MetalsDiskMass [sfh_bin_number_]) <= 0.)
        { disk_met_index_  = 0;  disk_f_met_1_ = 1.;  disk_f_met_2_ = 0.; }
        else  
        {
          const double log10_disk_metal_fraction_  = log10(metals_total(galaxy_->sfh_MetalsDiskMass [sfh_bin_number_]) / galaxy_->sfh_DiskMass [sfh_bin_number_]);
          find_metallicity_luminosity_interpolation_parameters(log10_disk_metal_fraction_ , disk_met_index_ , disk_f_met_1_ , disk_f_met_2_);
        }
      
#ifdef OUTPUT_REST_MAGS
        add_interpolated_luminosities(disk_mass_, disk_met_index_, disk_f_met_1_, disk_f_met_2_,
                                                      age_index_,      f_age_1_,      f_age_2_, 
                                                  redshift_index_, &LumDisk_, &YLumDisk_, is_affected_by_birthclould_);
#endif /* defined OUTPUT_REST_MAGS */     
#ifdef OUTPUT_OBS_MAGS           
        add_interpolated_luminosities(disk_mass_, disk_met_index_, disk_f_met_1_, disk_f_met_2_,
                                                      age_index_,      f_age_1_,      f_age_2_, 
                                                  redshift_index_obs_, &ObsLumDisk_, &ObsYLumDisk_, is_affected_by_birthclould_);
#ifdef OUTPUT_MOMAF_INPUTS   
        add_interpolated_luminosities(disk_mass_, disk_met_index_, disk_f_met_1_, disk_f_met_2_,
                                                      age_index_,      f_age_1_,      f_age_2_, 
                                                  redshift_index_dobs_, &dObsLumDisk_, &dObsYLumDisk_, is_affected_by_birthclould_);
#ifdef KITZBICHLER       
        add_interpolated_luminosities(disk_mass_, disk_met_index_, disk_f_met_1_, disk_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                  redshift_index_dobs_forward_, &dObsLumDisk_forward_, &dObsYLumDisk_forward_, is_affected_by_birthclould_);
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined OUTPUT_OBS_MAGS */
      }
    
      if(bulge_mass_ > 0.)
      {
        if(MetallicityOption == 0)
        { bulge_met_index_ = 4; bulge_f_met_1_ = 1.; bulge_f_met_2_ = 0.; }
        else if(metals_total(galaxy_->sfh_MetalsBulgeMass [sfh_bin_number_]) <= 0.)
        { bulge_met_index_ = 0; bulge_f_met_1_ = 1.; bulge_f_met_2_ = 0.; }
        else
        {
          const double log10_bulge_metal_fraction_ = log10(metals_total(galaxy_->sfh_MetalsBulgeMass[sfh_bin_number_]) / galaxy_->sfh_BulgeMass[sfh_bin_number_]);
          find_metallicity_luminosity_interpolation_parameters(log10_bulge_metal_fraction_, bulge_met_index_, bulge_f_met_1_, bulge_f_met_2_);
        }
      
#ifdef OUTPUT_REST_MAGS
        add_interpolated_luminosities(bulge_mass_, bulge_met_index_, bulge_f_met_1_, bulge_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_, &LumBulge_, &YLumBulge_, is_affected_by_birthclould_);
#endif /* defined OUTPUT_REST_MAGS */     
#ifdef OUTPUT_OBS_MAGS           
        add_interpolated_luminosities(bulge_mass_, bulge_met_index_, bulge_f_met_1_, bulge_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_obs_, &ObsLumBulge_, &ObsYLumBulge_, is_affected_by_birthclould_);
#ifdef OUTPUT_MOMAF_INPUTS   
        add_interpolated_luminosities(bulge_mass_, bulge_met_index_, bulge_f_met_1_, bulge_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_dobs_, &dObsLumBulge_, &dObsYLumBulge_, is_affected_by_birthclould_);
#ifdef KITZBICHLER       
        add_interpolated_luminosities(bulge_mass_, bulge_met_index_, bulge_f_met_1_, bulge_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_dobs_forward_, &dObsLumBulge_forward_, &dObsYLumBulge_forward_, is_affected_by_birthclould_);
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined OUTPUT_OBS_MAGS */
      }

#ifdef ICL
      if(icm_mass_ > 0.)
      {
        if(MetallicityOption == 0)
        { icm_met_index_   = 4;   icm_f_met_1_ = 1.;   icm_f_met_2_ = 0.; }
        else if(metals_total(galaxy_->sfh_MetalsICMMass [sfh_bin_number_]) <= 0.)
        { icm_met_index_   = 0;   icm_f_met_1_ = 1.;   icm_f_met_2_ = 0.; }
        else
        {
          const double log10_icm_metal_fraction_   = log10(metals_total(galaxy_->sfh_MetalsICM      [sfh_bin_number_]) / galaxy_->sfh_ICM      [sfh_bin_number_]);
          find_metallicity_luminosity_interpolation_parameters(log10_icm_metal_fraction_  , icm_met_index_  , icm_f_met_1_  , icm_f_met_2_  );
        }
        
#ifdef OUTPUT_REST_MAGS
        add_interpolated_luminosities(icm_mass_, icm_met_index_, icm_f_met_1_, icm_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_, &LumICL_, &YLumICL_dummy_, false);
#endif /* defined OUTPUT_REST_MAGS */     
#ifdef OUTPUT_OBS_MAGS           
        add_interpolated_luminosities(icm_mass_, icm_met_index_, icm_f_met_1_, icm_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_obs_, &ObsLumICL_, &YLumICL_dummy_, false);
#ifdef OUTPUT_MOMAF_INPUTS   
        add_interpolated_luminosities(icm_mass_, icm_met_index_, icm_f_met_1_, icm_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_dobs_, &dObsLumICL_, &YLumICL_dummy_, false);
#ifdef KITZBICHLER       
        add_interpolated_luminosities(icm_mass_, icm_met_index_, icm_f_met_1_, icm_f_met_2_,
                                                        age_index_,      f_age_1_,      f_age_2_, 
                                                    redshift_index_dobs_forward_, &dObsLumICL_forward_, &YLumICL_dummy_, false);
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined OUTPUT_OBS_MAGS */
      }
#endif /* defined ICL */ 
   
#ifdef OUTPUT_REST_MAGS
      galaxy_->rbandWeightAge += (age_ * (LumDisk_[r_band_filter_number_] + LumBulge_[r_band_filter_number_] - previous_r_band_luminosity_));
      previous_r_band_luminosity_ = LumDisk_[r_band_filter_number_] + LumBulge_[r_band_filter_number_];
#endif /* defined OUTPUT_REST_MAGS */

#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
    }
#endif /* defined N_FINE_AGE_BINS > 1 */
  }
  
#ifdef OUTPUT_REST_MAGS
  if(LumDisk_[r_band_filter_number_] + LumBulge_[r_band_filter_number_] > 0.)
  { 
    galaxy_->rbandWeightAge /= LumDisk_[r_band_filter_number_] + LumBulge_[r_band_filter_number_];
    galaxy_->rbandWeightAge *= UnitTime_in_Megayears / Hubble_h / 1000.; //conversion in age_ from code units/h -> Gyr
  }
#endif /* defined OUTPUT_REST_MAGS */    
  
  int i_;
#ifdef OUTPUT_REST_MAGS
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->Mag     [i_] = lum_to_lum_or_mag(LumBulge_[i_] + LumDisk_[i_]); }
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->MagBulge[i_] = lum_to_lum_or_mag(LumBulge_[i_]               ); }
#ifdef ICL
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->MagICL  [i_] = lum_to_lum_or_mag(LumICL_  [i_]               ); }
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */
#ifdef COMPUTE_OBS_MAGS
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->ObsMag     [i_] = lum_to_lum_or_mag(ObsLumBulge_[i_] + ObsLumDisk_[i_]); }
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->ObsMagBulge[i_] = lum_to_lum_or_mag(ObsLumBulge_[i_]                  ); }
#ifdef ICL
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->ObsMagICL  [i_] = lum_to_lum_or_mag(ObsLumICL_  [i_]                  ); }
#endif /* defined ICL */
#ifdef OUTPUT_MOMAF_INPUTS
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMag     [i_] = lum_to_lum_or_mag(dObsLumBulge_[i_] + dObsLumDisk_[i_]); }
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMagBulge[i_] = lum_to_lum_or_mag(dObsLumBulge_[i_]                   ); }
#ifdef ICL
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMagICL  [i_] = lum_to_lum_or_mag(dObsLumICL_  [i_]                   ); }
#endif /* defined ICL */
#ifdef KITZBICHLER
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMag_forward     [i_] = lum_to_lum_or_mag(dObsLumBulge_forward_[i_] + dObsLumDisk_forward_[i_]); }
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMagBulge_forward[i_] = lum_to_lum_or_mag(dObsLumBulge_forward_[i_]                           ); }
#ifdef ICL
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMagICL_forward  [i_] = lum_to_lum_or_mag(dObsLumICL_forward_  [i_]                           ); }
#endif /* defined ICL */
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined COMPUTE_OBS_MAGS */

  /* inclination is needed for disk dust correction,
   * but not sure, inclination shouldn't already be computed elsewhere */
  galaxy_->CosInclination = fabs(galaxy_->StellarSpin[2]) / sqrt(galaxy_->StellarSpin[0]*galaxy_->StellarSpin[0]+
                                                                 galaxy_->StellarSpin[1]*galaxy_->StellarSpin[1]+
                                                                 galaxy_->StellarSpin[2]*galaxy_->StellarSpin[2]);

  /* Dust correction for disks (Inter-stellar Medium  + birth clouds) */
  if(galaxy_->ColdGas > 0.0)
  {
    make_dust_correction_for_post_processing(galaxy_->SnapNum, Z_g_, galaxy_->ColdGas, galaxy_->GasDiskRadius, galaxy_->CosInclination
#ifdef OUTPUT_REST_MAGS
                                           , &LumDisk_            , YLumDisk_           
#endif /* defined OUTPUT_REST_MAGS */
#ifdef COMPUTE_OBS_MAGS
                                           , &ObsLumDisk_         , ObsYLumDisk_        
#ifdef OUTPUT_MOMAF_INPUTS
                                           , &dObsLumDisk_        , dObsYLumDisk_       
#ifdef KITZBICHLER
                                           , &dObsLumDisk_forward_, dObsYLumDisk_forward_
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined COMPUTE_OBS_MAGS */
                                            );
  }

   /* dust correction for bulges (only remove light from young stars absorbed by birth clouds) */
#ifdef OUTPUT_REST_MAGS
  for(i_= 0; i_ < NMAG; i_++) { LumBulge_            [i_] -= YLumBulge_            [i_] * (1. - ExpTauBCBulge); }
#endif /* defined OUTPUT_REST_MAGS */
#ifdef COMPUTE_OBS_MAGS
  for(i_= 0; i_ < NMAG; i_++) { ObsLumBulge_         [i_] -= ObsYLumBulge_         [i_] * (1. - ExpTauBCBulge); }
#ifdef OUTPUT_MOMAF_INPUTS
  for(i_= 0; i_ < NMAG; i_++) { dObsLumBulge_        [i_] -= dObsYLumBulge_        [i_] * (1. - ExpTauBCBulge); }
#ifdef KITZBICHLER
  for(i_= 0; i_ < NMAG; i_++) { dObsLumBulge_forward_[i_] -= dObsYLumBulge_forward_[i_] * (1. - ExpTauBCBulge); }
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined COMPUTE_OBS_MAGS */
                                         
#ifdef OUTPUT_REST_MAGS
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->MagDust             [i_] = lum_to_lum_or_mag(LumBulge_            [i_] + LumDisk_            [i_]); }
#endif /* defined OUTPUT_REST_MAGS */                        
#ifdef COMPUTE_OBS_MAGS                                      
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->ObsMagDust          [i_] = lum_to_lum_or_mag(ObsLumBulge_         [i_] + ObsLumDisk_         [i_]); }
#ifdef OUTPUT_MOMAF_INPUTS                                   
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMagDust         [i_] = lum_to_lum_or_mag(dObsLumBulge_        [i_] + dObsLumDisk_        [i_]); }
#ifdef KITZBICHLER                                           
  for(i_= 0; i_ < NMAG; i_++) { galaxy_->dObsMagDust_forward [i_] = lum_to_lum_or_mag(dObsLumBulge_forward_[i_] + dObsLumDisk_forward_[i_]); }
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined COMPUTE_OBS_MAGS */
}
