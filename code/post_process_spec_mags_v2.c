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

#include "allvars.h"
#include "proto.h"

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


/** @brief adds metals if detailed metals on */
 #ifdef DETAILED_METALS_AND_MASS_RETURN
static inline float 
get_metals_total(struct metals metals_)
{ return(metals_.type1a+metals_.type2+metals_.agb); }
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
static inline float 
get_metals_total(const float metals_)
{ return(metals_); }
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */


/** @brief Used by add_to_luminosities() to interpolates into
 *         the age in the SSP tables.*/
static inline void 
find_age_luminosity_interpolation_parameters(double log10_age_, int *age_index_, double *f_age_1_, double *f_age_2_)
{
  if(log10_age_ > SSP_logAgeTab[SSP_NAGES - 1])        /* beyond table, take latest entry */
  {
    *age_index_ = SSP_NAGES - 2;
    *f_age_1_ = 0;
    *f_age_2_ = 1;
  }
  else if(log10_age_ < SSP_logAgeTab[1])        /* age younger than 1st enty, take 1st entry */
  {
    *age_index_ = 0;
    *f_age_1_ = 0;
    *f_age_2_ = 1;
  }
  else
  {
    int idx_ = get_jump_index(log10_age_);
    while(SSP_logAgeTab[idx_ + 1] < log10_age_) idx_++;
    *age_index_ = idx_;
    const double frac_ = (log10_age_ - SSP_logAgeTab[idx_]) / (SSP_logAgeTab[idx_ + 1] - SSP_logAgeTab[idx_]);
    *f_age_1_ = 1 - frac_;
    *f_age_2_ =     frac_;
  }
}


/** @brief Used by add_to_luminosities() to interpolates into
 *         the metallicity in the SSP tables.*/
static inline void 
find_metallicity_luminosity_interpolation_parameters(const double log10_metallicity_, int *met_index_, double *f_met_1_, double *f_met_2_)
{
  if(log10_metallicity_ > SSP_logMetalTab[SSP_NMETALLICITES - 1])        /* beyond table, take latest entry */
  {
    *met_index_ = SSP_NMETALLICITES - 2;
    *f_met_1_ = 0;
    *f_met_2_ = 1;
  }
  else if(log10_metallicity_ < SSP_logMetalTab[0])        /* mettallicity smaller 1st enty, take 1st entry */
  {
    *met_index_ = 0;
    *f_met_1_ = 1;
    *f_met_2_ = 0;
  }
  else
  {
    int idx_ = 0;
    while(SSP_logMetalTab[idx_ + 1] < log10_metallicity_) idx_++;
    *met_index_ = idx_;
    const double frac_ = (log10_metallicity_ - SSP_logMetalTab[idx_]) / (SSP_logMetalTab[idx_ + 1] - SSP_logMetalTab[idx_]);
    *f_met_1_ = 1 - frac_;
    *f_met_2_ =     frac_;
  }
}


/** @brief computes lumsinosities for galaxies */
static inline void 
//compute_post_process_luminosities (const double mass_, const double age_, const double metal_fraction_, const int snapshot_number_, const int filter_number_,
compute_post_process_luminosities (const double mass_, const double age_, 
                                   const int age_index_, const double f_age_1_, const double f_age_2_,
                                   const int met_index_, const double f_met_1_, const double f_met_2_,
                                   const int snapshot_number_, const int filter_number_,
                                   double *Lum, double *ObsLum, double *dObsLum, double *dObsLum_forward,
                                   double *YLum, double *ObsYLum, double *dObsYLum, double *dObsYLum_forward)
{
  int  redshift_index_;

  /* Time below which the luminosities are corrected for extinction due to
   * molecular birth clouds.  */
  const double birthcloud_age_ = 10.0 / UnitTime_in_Megayears * Hubble_h;

#ifdef OUTPUT_REST_MAGS
  redshift_index_   = 0;
  const double LumToAdd  = mass_ * (f_met_1_ * (f_age_1_ * LumTables[filter_number_][met_index_    ][redshift_index_][age_index_    ]  +
                                                f_age_2_ * LumTables[filter_number_][met_index_    ][redshift_index_][age_index_ + 1]) +
                                    f_met_2_ * (f_age_1_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_    ]  +
                                                f_age_2_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_ + 1])  );

  *Lum += LumToAdd;
  if(age_<=birthcloud_age_)
    *YLum += LumToAdd;

#else /* not defined OUTPUT_REST_MAGS */
  (void)Lum; /* suppress unused-parameter warning */
  (void)YLum; /* suppress unused-parameter warning */
#endif /* not defined OUTPUT_REST_MAGS */

#ifdef COMPUTE_OBS_MAGS
  redshift_index_   = (LastDarkMatterSnapShot+1) - 1 - snapshot_number_;
  const double ObsLumToAdd = mass_ * (f_met_1_ * (f_age_1_ * LumTables[filter_number_][met_index_][redshift_index_][age_index_] +
                                 f_age_2_ * LumTables[filter_number_][met_index_][redshift_index_][age_index_ + 1]) +
                        f_met_2_ * (f_age_1_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_] +
                                 f_age_2_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_ + 1]));

  *ObsLum += ObsLumToAdd;
  if(age_<=birthcloud_age_)
    *ObsYLum += ObsLumToAdd;

#ifdef OUTPUT_MOMAF_INPUTS
  double dObsLumToAdd;

  redshift_index_   = (LastDarkMatterSnapShot+1) - 1 - (snapshot_number_-1);
  dObsLumToAdd = mass_ * (f_met_1_ * (f_age_1_ * LumTables[filter_number_][met_index_][redshift_index_][age_index_] +
                                  f_age_2_ * LumTables[filter_number_][met_index_][redshift_index_][age_index_ + 1]) +
                         f_met_2_ * (f_age_1_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_] +
                                  f_age_2_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_ + 1]));

  *dObsLum += dObsLumToAdd;
  if(age_<=birthcloud_age_)
    *dObsYLum += dObsLumToAdd;

#ifdef KITZBICHLER
  redshift_index_   = (LastDarkMatterSnapShot+1) - 1 - (snapshot_number_+1);
  if(redshift_index_ < 0) redshift_index_ = 0;
  dObsLumToAdd = mass_ * (f_met_1_ * (f_age_1_ * LumTables[filter_number_][met_index_][redshift_index_][age_index_] +
                                  f_age_2_ * LumTables[filter_number_][met_index_][redshift_index_][age_index_ + 1]) +
                         f_met_2_ * (f_age_1_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_] +
                                  f_age_2_ * LumTables[filter_number_][met_index_ + 1][redshift_index_][age_index_ + 1]));

  *dObsLum_forward += dObsLumToAdd;
  if(age_<=birthcloud_age_)
    *dObsYLum_forward += dObsLumToAdd;

#else  /* not defined KITZBICHLER */
  (void)dObsLum_forward; /* suppress unused-parameter warning */
  (void)dObsYLum_forward; /* suppress unused-parameter warning */
#endif /* not defined KITZBICHLER */

#else  /* not defined OUTPUT_MOMAF_INPUTS */
  (void)dObsLum; /* suppress unused-parameter warning */
  (void)dObsYLum; /* suppress unused-parameter warning */
#endif /* not defined OUTPUT_MOMAF_INPUTS */

#else  /* not defined COMPUTE_OBS_MAGS */
  (void)ObsLum; /* suppress unused-parameter warning */
  (void)ObsYLum; /* suppress unused-parameter warning */
#endif /* not defined COMPUTE_OBS_MAGS */
}


/** @brief computes dust corrections for lums. */
static inline void
make_dust_correction_for_post_processing(int filter_number_, int snapshot_number_, double Z_g_, double ColdGas, double GasDiskRadius, double CosInclination, 
                                         double *LumDisk            , double YLumDisk,
                                         double *ObsLumDisk         , double ObsYLumDisk,
                                         double *dObsLumDisk        , double dObsYLumDisk,
                                         double *dObsLumDisk_forward, double dObsYLumDisk_forward)
{
  double tau, alam;
  double taubc, mu;
  
  const double VBand_WaveLength = 0.55;

  /* 0.94 = 2.83/3. - 3 to get scale lenght and 2.83 = 1.68^2 */
  /* nh = ColdGas / (M_PI * pow(GasDiskRadius * 0.56, 2) * 1.4); */
  /* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (2.1 10^21 atoms/cm^2) */
  /* note: 3252.37 = 10^(3.5122) */
  /* add  redshift dependence of dust-to-ColdGas ratio */
 const double nh = ColdGas / (M_PI * pow(GasDiskRadius * 0.94, 2) * 1.4 * 3252.37 * (1 + ZZ[snapshot_number_]));

 // minimum inclination ~80 degrees
  const double sec = (CosInclination < 0.2) ? (1. / 0.2) : (1. / CosInclination);

  /* mu for YS extinction, given by a Gaussian with centre 0.3 (MUCENTER)
   * and width 0.2 (MUWIDTH), truncated at 0.1 and 1.  */
  do { mu = gasdev(&mu_seed) * MUWIDTH + MUCENTER; }
  while (mu < 0.1 || mu > 1.0);
  
  // for testing:
  mu = MUCENTER;

  // extinction on Vband used as reference for the BC extinction
  const double tauvbc = get_extinction(NMAG, Z_g_, 0) * nh * (1. / mu - 1.);

#ifdef OUTPUT_REST_MAGS

  //extinction due to ISM and light from young stars absorbed by birth clouds
  taubc = tauvbc * pow(FilterLambda[filter_number_] / VBand_WaveLength, -0.7);
  tau   = get_extinction(filter_number_, Z_g_, 0) * nh * sec;
  alam  = (tau > 0.0) ? (1.0 - exp(-tau)) / tau : 1.;

  *LumDisk -= (YLumDisk) * (1. - exp(-taubc));
  *LumDisk *= alam;
  
#else  /* not defined OUTPUT_REST_MAGS */
  (void)LumDisk; /* suppress unused-parameter warning */
  (void)YLumDisk; /* suppress unused-parameter warning */
#endif /* not defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS

  //extinction due to ISM and light from young stars absorbed by birth clouds
  taubc = tauvbc * pow((FilterLambda[filter_number_] * (1. + ZZ[snapshot_number_])) / VBand_WaveLength, -0.7);
  tau   = get_extinction(filter_number_, Z_g_, ZZ[snapshot_number_]) * nh * sec;
  alam  = (tau > 0.0) ? (1.0 - exp(-tau)) / tau : 1.;

  *ObsLumDisk -= (ObsYLumDisk) * (1. - exp(-taubc));
  *ObsLumDisk *= alam;

#ifdef OUTPUT_MOMAF_INPUTS   // compute same thing at z + 1
  const int earlier_snapshot_number_ = (snapshot_number_ > 0) ? (snapshot_number_ - 1) : 0;
  taubc = tauvbc * pow((FilterLambda[filter_number_] * (1. + ZZ[earlier_snapshot_number_])) / VBand_WaveLength, -0.7);
  tau   = get_extinction(filter_number_, Z_g_, ZZ[earlier_snapshot_number_]) * nh * sec;
  alam  = (tau > 0.0) ? (1.0 - exp(-tau)) / tau : 1.;

  *dObsLumDisk -= (dObsYLumDisk) * (1. - exp(-taubc));
  *dObsLumDisk *= alam;

#ifdef KITZBICHLER
  const int later_snapshot_number_ = (snapshot_number_ < LastDarkMatterSnapShot) ? (snapshot_number_ + 1) : LastDarkMatterSnapShot;
  taubc = tauvbc * pow((FilterLambda[filter_number_] * (1. + ZZ[later_snapshot_number_])) / VBand_WaveLength, -0.7);
  tau   = get_extinction(filter_number_, Z_g_, ZZ[later_snapshot_number_]) * nh * sec;
  alam  = (tau > 0.0) ? (1.0 - exp(-tau)) / tau : 1.;

  *dObsLumDisk_forward -= (dObsYLumDisk_forward) * (1. - exp(-taubc));
  *dObsLumDisk_forward *= alam;
  
#else  /* not defined KITZBICHLER */
  (void)dObsLumDisk_forward; /* suppress unused-parameter warning */
  (void)dObsYLumDisk_forward; /* suppress unused-parameter warning */
#endif /* not defined KITZBICHLER */

#else  /* not defined OUTPUT_MOMAF_INPUTS */
  (void)dObsLumDisk; /* suppress unused-parameter warning */
  (void)dObsYLumDisk; /* suppress unused-parameter warning */
#endif /* not defined OUTPUT_MOMAF_INPUTS */

#else  /* not defined OUTPUT_OBS_MAGS */
  (void)ObsLumDisk; /* suppress unused-parameter warning */
  (void)ObsYLumDisk; /* suppress unused-parameter warning */
#endif /* not defined OUTPUT_OBS_MAGS */
}


/** @brief computes lums./mags. for galaxies */
void post_process_spec_mags(struct GALAXY_OUTPUT *galaxy_)
{
  int sfh_bin_number_, filter_number_;

  //used for dust corrections
  const double Z_g_ = metals_total(galaxy_->MetalsColdGas)/galaxy_->ColdGas/0.02;
  const double t_snapshot_ = NumToTime(galaxy_->SnapNum);

  int age_index_       [SFH_NBIN];
  double f_age_1_      [SFH_NBIN];
  double f_age_2_      [SFH_NBIN];
  int disk_met_index_  [SFH_NBIN];
  double disk_f_met_1_ [SFH_NBIN];
  double disk_f_met_2_ [SFH_NBIN];
  int bulge_met_index_ [SFH_NBIN];
  double bulge_f_met_1_[SFH_NBIN];
  double bulge_f_met_2_[SFH_NBIN];
#ifdef ICL
  int icm_met_index_   [SFH_NBIN];
  double icm_f_met_1_  [SFH_NBIN];
  double icm_f_met_2_  [SFH_NBIN];
#endif /* defined ICL */

  for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
  {
    const double age_                  = SFH_t[galaxy_->SnapNum][0][sfh_bin_number_]+SFH_dt[galaxy_->SnapNum][0][sfh_bin_number_]/2 - t_snapshot_;
    const double disk_metal_fraction_  = get_metals_total(galaxy_->sfh_MetalsDiskMass [sfh_bin_number_]) / galaxy_->sfh_DiskMass [sfh_bin_number_];
    const double bulge_metal_fraction_ = get_metals_total(galaxy_->sfh_MetalsBulgeMass[sfh_bin_number_]) / galaxy_->sfh_BulgeMass[sfh_bin_number_];
#ifdef ICL
    const double icm_metal_fraction_   = get_metals_total(galaxy_->sfh_MetalsICM      [sfh_bin_number_]) / galaxy_->sfh_ICM      [sfh_bin_number_];
#endif /* defined ICL */
    
    find_age_luminosity_interpolation_parameters        (log10(age_)                 , &(age_index_      [sfh_bin_number_]), &(f_age_1_      [sfh_bin_number_]), &(f_age_2_      [sfh_bin_number_]));
    find_metallicity_luminosity_interpolation_parameters(log10(disk_metal_fraction_ ), &(disk_met_index_ [sfh_bin_number_]), &(disk_f_met_1_ [sfh_bin_number_]), &(disk_f_met_2_ [sfh_bin_number_]));
    find_metallicity_luminosity_interpolation_parameters(log10(bulge_metal_fraction_), &(bulge_met_index_[sfh_bin_number_]), &(bulge_f_met_1_[sfh_bin_number_]), &(bulge_f_met_2_[sfh_bin_number_]));
#ifdef ICL
    find_metallicity_luminosity_interpolation_parameters(log10(icm_metal_fraction_  ), &(icm_met_index_  [sfh_bin_number_]), &(icm_f_met_1_  [sfh_bin_number_]), &(icm_f_met_2_  [sfh_bin_number_]));
#endif /* defined ICL */
  }

  // reset met index to use only solar metallicity
  if(MetallicityOption == 0)
  {
    for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
    {
      disk_met_index_ [sfh_bin_number_] = 4;
      bulge_met_index_[sfh_bin_number_] = 4;
#ifdef ICL
      icm_met_index_  [sfh_bin_number_] = 4;
#endif /* defined ICL */
    }
  }
    
  galaxy_->rbandWeightAge=0.0;

  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  {
    double LumDisk  = 0.; double YLumDisk  = 0.;
    double LumBulge = 0.; double YLumBulge = 0.;
    
    double ObsLumDisk  = 0.; double ObsYLumDisk  = 0.;
    double ObsLumBulge = 0.; double ObsYLumBulge = 0.;
    
    double dObsLumDisk  = 0.; double dObsYLumDisk  = 0.;
    double dObsLumBulge = 0.; double dObsYLumBulge = 0.;
    
    double dObsLumDisk_forward  = 0.; double dObsYLumDisk_forward  = 0.;
    double dObsLumBulge_forward = 0.; double dObsYLumBulge_forward = 0.;

#ifdef ICL
   double LumICL   = 0.; double YLumICL   = 0.;
   double ObsLumICL   = 0.; double ObsYLumICL   = 0.;
   double dObsLumICL   = 0.; double dObsYLumICL   = 0.;
   double dObsLumICL_forward   = 0.; double dObsYLumICL_forward   = 0.;
#endif /* defined ICL */

    double previous_lum_=0.;

    //loop on SFH bins
    for(sfh_bin_number_ = 0; sfh_bin_number_ <= galaxy_->sfh_ibin; sfh_bin_number_++)
    {
      const double age_ = SFH_t[galaxy_->SnapNum][0][sfh_bin_number_]+SFH_dt[galaxy_->SnapNum][0][sfh_bin_number_]/2 - t_snapshot_;

      /* The stellar populations tables have magnitudes for all the mass
       * formed in stars including what will be shortly lost by SNII   */
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

      //Disk MAGS
      if(galaxy_->sfh_DiskMass[sfh_bin_number_] > 0.0)
        compute_post_process_luminosities (disk_mass_, age_,
                                           age_index_     [sfh_bin_number_], f_age_1_     [sfh_bin_number_], f_age_2_     [sfh_bin_number_], 
                                           disk_met_index_[sfh_bin_number_], disk_f_met_1_[sfh_bin_number_], disk_f_met_2_[sfh_bin_number_], 
                                           galaxy_->SnapNum, filter_number_,
                                           &LumDisk, &ObsLumDisk, &dObsLumDisk, &dObsLumDisk_forward,
                                           &YLumDisk, &ObsYLumDisk, &dObsYLumDisk, &dObsYLumDisk_forward);

      //Bulge MAGS
      if(galaxy_->sfh_BulgeMass[sfh_bin_number_] > 0.0)
        compute_post_process_luminosities (bulge_mass_, age_, 
                                           age_index_     [sfh_bin_number_], f_age_1_     [sfh_bin_number_], f_age_2_     [sfh_bin_number_], 
                                           bulge_met_index_[sfh_bin_number_], bulge_f_met_1_[sfh_bin_number_], bulge_f_met_2_[sfh_bin_number_], 
                                           galaxy_->SnapNum, filter_number_,
                                           &LumBulge, &ObsLumBulge, &dObsLumBulge, &dObsLumBulge_forward,
                                           &YLumBulge, &ObsYLumBulge, &dObsYLumBulge, &dObsYLumBulge_forward);
#ifdef ICL
      //ICL MAGS
      if(galaxy_->sfh_ICM[sfh_bin_number_] > 0.0)
        compute_post_process_luminosities (icm_mass_, age_, 
                                           age_index_     [sfh_bin_number_], f_age_1_     [sfh_bin_number_], f_age_2_     [sfh_bin_number_], 
                                           icm_met_index_ [sfh_bin_number_], icm_f_met_1_ [sfh_bin_number_], icm_f_met_2_ [sfh_bin_number_], 
                                           galaxy_->SnapNum, filter_number_,
                                           &LumICL, &ObsLumICL, &dObsLumICL, &dObsLumICL_forward,
                                           &YLumICL, &ObsYLumICL, &dObsYLumICL, &dObsYLumICL_forward);
#endif /* defined ICL */

      //r-band weighted ages and mass_ weighted
      /** @bug  filter_number_==17 is not necessarily r-band, since filters are assigned from info in parameter file
      *
      *  @todo  parametrize filter band number for computing weighted stellar age_
      */
      if((galaxy_->DiskMass+galaxy_->BulgeMass)>0.0 && filter_number_==17)
      {
        galaxy_->rbandWeightAge += (age_ * (LumDisk+LumBulge-previous_lum_));
        previous_lum_=LumDisk+LumBulge;
      }
    }//end of loop on sfh bins

    if((galaxy_->DiskMass+galaxy_->BulgeMass)>0.0 && filter_number_==17)
    {
      //LumDisk & LumBulge are sdss r-band luminosities (filter_number_==17)
      galaxy_->rbandWeightAge /= (LumDisk+LumBulge);
      galaxy_->rbandWeightAge = galaxy_->rbandWeightAge / 1000. * UnitTime_in_Megayears / Hubble_h; //conversion in age_ from code units/h -> Gyr
    }

#ifdef OUTPUT_REST_MAGS
    galaxy_->Mag[filter_number_]=lum_to_lum_or_mag(LumDisk+LumBulge);
    galaxy_->MagBulge[filter_number_]=lum_to_lum_or_mag(LumBulge);
#ifdef ICL
    galaxy_->MagICL[filter_number_]=lum_to_lum_or_mag(LumICL);
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */
#ifdef COMPUTE_OBS_MAGS
    galaxy_->ObsMag[filter_number_]=lum_to_lum_or_mag(ObsLumDisk+ObsLumBulge);
    galaxy_->ObsMagBulge[filter_number_]=lum_to_lum_or_mag(ObsLumBulge);
#ifdef ICL
    galaxy_->ObsMagICL[filter_number_]=lum_to_lum_or_mag(ObsLumICL);
#endif /* defined ICL */

#ifdef OUTPUT_MOMAF_INPUTS
    galaxy_->dObsMag[filter_number_]=lum_to_lum_or_mag(dObsLumDisk+dObsLumBulge);
    galaxy_->dObsMagBulge[filter_number_]=lum_to_lum_or_mag(dObsLumBulge);
#ifdef ICL
    galaxy_->dObsMagICL[filter_number_]=lum_to_lum_or_mag(dObsLumICL);
#endif /* defined ICL */
#ifdef KITZBICHLER
    galaxy_->dObsMag_forward[filter_number_]=lum_to_lum_or_mag(dObsLumDisk_forward+dObsLumBulge_forward);
    galaxy_->dObsMagBulge_forward[filter_number_]=lum_to_lum_or_mag(dObsLumBulge_forward);
#ifdef ICL
    galaxy_->dObsMagICL_forward[filter_number_]=lum_to_lum_or_mag(dObsLumICL_forward);
#endif /* defined ICL */
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined COMPUTE_OBS_MAGS */

    galaxy_->CosInclination = fabs(galaxy_->StellarSpin[2]) /
        sqrt(galaxy_->StellarSpin[0]*galaxy_->StellarSpin[0]+
             galaxy_->StellarSpin[1]*galaxy_->StellarSpin[1]+
             galaxy_->StellarSpin[2]*galaxy_->StellarSpin[2]);
                                          
    if(galaxy_->ColdGas > 0.0)
      make_dust_correction_for_post_processing(filter_number_, galaxy_->SnapNum, Z_g_, galaxy_->ColdGas, galaxy_->GasDiskRadius, galaxy_->CosInclination,
                                               &LumDisk            , YLumDisk,
                                               &ObsLumDisk         , ObsYLumDisk,
                                               &dObsLumDisk        , dObsYLumDisk,
                                               &dObsLumDisk_forward, dObsYLumDisk_forward);

    //Dust correction for bulges (remove light from young stars absorbed by birth clouds)
    LumBulge             -= YLumBulge             * (1. - ExpTauBCBulge);
    ObsLumBulge          -= ObsYLumBulge          * (1. - ExpTauBCBulge);
    dObsLumBulge         -= dObsYLumBulge         * (1. - ExpTauBCBulge);
    dObsLumBulge_forward -= dObsYLumBulge_forward * (1. - ExpTauBCBulge);

#ifdef OUTPUT_REST_MAGS
    galaxy_->MagDust            [filter_number_] = lum_to_lum_or_mag(LumDisk+LumBulge);
#endif /* defined OUTPUT_REST_MAGS */
#ifdef COMPUTE_OBS_MAGS
    galaxy_->ObsMagDust         [filter_number_] = lum_to_lum_or_mag(ObsLumDisk+ObsLumBulge);
#ifdef OUTPUT_MOMAF_INPUTS
    galaxy_->dObsMagDust        [filter_number_] = lum_to_lum_or_mag(dObsLumDisk+dObsLumBulge);
#ifdef KITZBICHLER
    galaxy_->dObsMagDust_forward[filter_number_] = lum_to_lum_or_mag(dObsLumDisk_forward+dObsLumBulge_forward);
#endif /* defined KITZBICHLER */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined COMPUTE_OBS_MAGS */
      
  }//end of loop on bands
}
