/*  Copyright (C) <2016>  <L-Galaxies>
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

#include "model_dust_extinction_inline.h"


/** @file model_dust.c
 *  @brief model_dust.c is used to compute dust extinction as described
 *         in Delucia2007 + redshift dependence as Kitzbichler & White 2007.

 *  There are 2 extinction sources:
 *  Extinction from a diffuse inter-stellar medium (ISM) (Devriendt1999);
 *  Extinction from molecular clouds in young stars (YS) (Charlot2000);
 *  Both were introduced in Delucia2007.

 *  The optical depth of dust in each component
 *  \f$\tau^z_{\lambda}\f$(ISM)
 *  and \f$\tau_{\lambda}^{\rm{BC}}\f$(YS)
 *  is used to compute extinction assuming a slab geometry for the dust and
 *  a random inclination of the disk to the line of sight.

 *  Extinction curves for the ISM:
 *  \f$\left(\frac{A_{\lambda}}{A_{\rm{v}}}\right)_{Z_{\odot}}
    \left(\frac{Z_{\rm{gas}}}{Z_{\odot}}\right)^s\f$
 *  are computed in get_extinction
 *
 *  The optical depth for the ISM at a given \f$\lambda\f$ can be written as:

    \f$\tau_{\lambda}^{ISM}=\left(\frac{A_{\lambda}}{A_v}\right)_{Z_{\odot}}
    \left(\frac{Z_{\rm{gas}}}{Z_{\odot}}\right)^s\left(\frac{\langle N_H\rangle}
    {2.1 \times10^{21}{\rm{atoms}} \,{\rm{cm}}^{-2}}\right)\f$,

    where the mean column density of Hydrogen is:

    \f$\langle N_H\rangle=\frac{M_{\rm{cold}}}{1.4\,m_p\pi
    (a R_{\mathrm{D}})^2}{\rm{atoms}}\, {\rm{cm}}^{-2}.\f$


 *  The optical depth for YS (\f$\tau_{\lambda}^{\rm{BC}}\f$) is calibrated
 *  from the ISM optical depth in the V-band:
 *
 *  \f$\tau_{\lambda}^{BC}=\tau_{\rm{v}}^{\rm{ISM}}\left(\frac{1}{\mu}-1\right)
 *  \left(\frac{\lambda}{5500 \AA}\right)^{-0.7}\f$,
 */


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS

/** @brief main routine where the extinction is calculated */
void dust_model(const int galaxy_number_, const int output_number_)
{
#ifndef OUTPUT_REST_MAGS
#ifndef OUTPUT_OBS_MAGS
  return; /* early return if no mags to compute dust for */
#endif /* not defined OUTPUT_OBS_MAGS */
#endif /* not defined OUTPUT_REST_MAGS */ 
  
  double mu_, tau_, alam_, a_bc_;
  double Lum_disk_, dly_;
  int filter_number_;
  
  const double inv_VBand_WaveLength_ = 1. / 0.55;

  if(Gal[galaxy_number_].ColdGas > 0.0)
  {
    // /* 0.94 = 2.83/3. -> 3 to get scale lenght and 2.83 = 1.68^2 */
    // n_h_ = Gal[galaxy_number_].ColdGas / (M_PI * pow(Gal[galaxy_number_].GasDiskRadius * 0.94, 2) * 1.4);
    // /* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (2.1 10^21 atoms/cm^2) */
    // n_h_ = n_h_ / 3252.37;        // 3252.37 = 10^(3.5122) ... ha ha ! 
    // /*redshift dependence */
    // n_h_ = n_h_ * pow(1 + ZZ[Halo[halo_number_].SnapNum], -1.0);
    
    const int snapshot_number_ = ListOutputSnaps[output_number_];
    
    const double n_h_ = Gal[galaxy_number_].ColdGas / ((M_PI * 0.94 * 0.94 * 1.4 * 3252.37) * Gal[galaxy_number_].GasDiskRadius * Gal[galaxy_number_].GasDiskRadius * (1 + ZZ[snapshot_number_]));

    Gal[galaxy_number_].CosInclination = fabs(Gal[galaxy_number_].StellarSpin[2]) /
                                sqrt(Gal[galaxy_number_].StellarSpin[0]*Gal[galaxy_number_].StellarSpin[0]+
                                     Gal[galaxy_number_].StellarSpin[1]*Gal[galaxy_number_].StellarSpin[1]+
                                     Gal[galaxy_number_].StellarSpin[2]*Gal[galaxy_number_].StellarSpin[2]);

    /* minimum inclination ~80 degrees, i.e. cos(inclination) == 0.2 */
    const double n_h_sec_ = (Gal[galaxy_number_].CosInclination < 0.2) ? (5. * n_h_) : (n_h_ / Gal[galaxy_number_].CosInclination);

    const double Z_g_ = metals_total(Gal[galaxy_number_].MetalsColdGas) / Gal[galaxy_number_].ColdGas / 0.02;

    /* mu_ for YS extinction, given by a Gaussian with centre 0.3 (MUCENTER)
     * and width 0.2 (MUWIDTH), truncated at 0.1 and 1.  */
    do { mu_ = gsl_ran_gaussian(random_generator, MUWIDTH) + MUCENTER; }
    while (mu_ < 0.1 || mu_ > 1.0);
    
    //for testing:
    mu_ = MUCENTER;

    const double tauvbc_ = get_extinction(NMAG, Z_g_, 0) * n_h_ * (1. / mu_ - 1.);

#ifdef OUTPUT_REST_MAGS
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      tau_  = get_extinction(filter_number_, Z_g_, 0) * n_h_sec_;
      alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
      a_bc_ = 1. - exp(-tauvbc_ * pow(FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

      Lum_disk_ = Gal[galaxy_number_].Lum[output_number_][filter_number_] - Gal[galaxy_number_].LumBulge[output_number_][filter_number_];
      Gal[galaxy_number_].LumDust[output_number_][filter_number_] = Gal[galaxy_number_].LumBulge[output_number_][filter_number_] + Lum_disk_ * alam_;

      // now remove light from young stars absorbed by birth clouds

      dly_ = (Gal[galaxy_number_].LumY     [output_number_][filter_number_] - Gal[galaxy_number_].LumBulgeY[output_number_][filter_number_]) * alam_ * a_bc_  +
              Gal[galaxy_number_].LumBulgeY[output_number_][filter_number_] * (1. - ExpTauBCBulge);

      Gal[galaxy_number_].LumDust[output_number_][filter_number_] -= dly_;
      
    /*
      Gal[galaxy_number_].LumDust[output_number_][filter_number_] = 
      Gal[galaxy_number_].LumBulge[output_number_][filter_number_] + (Gal[galaxy_number_]. Lum[output_number_][filter_number_] - Gal[galaxy_number_]. LumBulge[output_number_][filter_number_]) * alam_
             -                                                       (Gal[galaxy_number_].LumY[output_number_][filter_number_] - Gal[galaxy_number_].LumBulgeY[output_number_][filter_number_]) * alam_ * a_bc_
             - Gal[galaxy_number_].LumBulgeY[output_number_][filter_number_] * (1. - ExpTauBCBulge)
             
     = Gal[galaxy_number_]. Lum     [output_number_][filter_number_] * alam_
     + Gal[galaxy_number_]. LumBulge[output_number_][filter_number_] * (1 - alam_)
     - Gal[galaxy_number_].LumY     [output_number_][filter_number_] * alam_ * a_bc_
     - Gal[galaxy_number_].LumBulgeY[output_number_][filter_number_] * (1. - ExpTauBCBulge - alam_ * a_bc_)
        */     
    }
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      tau_  = get_extinction(filter_number_, Z_g_, ZZ[snapshot_number_]) * n_h_sec_;
      alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
      a_bc_ = 1. - exp(-tauvbc_ * pow((1. + ZZ[snapshot_number_]) * FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

      Lum_disk_ = Gal[galaxy_number_].ObsLum[output_number_][filter_number_] - Gal[galaxy_number_].ObsLumBulge[output_number_][filter_number_];
      Gal[galaxy_number_].ObsLumDust[output_number_][filter_number_] = Gal[galaxy_number_].ObsLumBulge[output_number_][filter_number_] + Lum_disk_ * alam_;

      // now remove light from young stars absorbed by birth clouds

      dly_ = (Gal[galaxy_number_].ObsLumY[output_number_][filter_number_] - Gal[galaxy_number_].ObsLumBulgeY[output_number_][filter_number_]) * alam_ * a_bc_  +
                      Gal[galaxy_number_].ObsLumBulgeY[output_number_][filter_number_] * (1. - ExpTauBCBulge);

      Gal[galaxy_number_].ObsLumDust[output_number_][filter_number_] -= dly_;
    }

#ifdef OUTPUT_FB_OBS_MAGS   // compute same thing at z (snapshot_number_ +/- 1)
    const int earlier_snapshot_number_ = (snapshot_number_ > 0) ? (snapshot_number_ - 1) : 0;
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      tau_  = get_extinction(filter_number_, Z_g_, ZZ[earlier_snapshot_number_]) * n_h_sec_;
      alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
      a_bc_ = 1. - exp(-tauvbc_ * pow((1. + ZZ[earlier_snapshot_number_]) * FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

      Lum_disk_ = Gal[galaxy_number_].backward_ObsLum[output_number_][filter_number_] - Gal[galaxy_number_].backward_ObsLumBulge[output_number_][filter_number_];
      Gal[galaxy_number_].backward_ObsLumDust[output_number_][filter_number_] = Gal[galaxy_number_].backward_ObsLumBulge[output_number_][filter_number_] + Lum_disk_ * alam_;

      dly_ = (Gal[galaxy_number_].backward_ObsLumY[output_number_][filter_number_] - Gal[galaxy_number_].backward_ObsLumBulgeY[output_number_][filter_number_]) * alam_ * a_bc_ +
                      Gal[galaxy_number_].backward_ObsLumBulgeY[output_number_][filter_number_] * (1. - ExpTauBCBulge);

      Gal[galaxy_number_].backward_ObsLumDust[output_number_][filter_number_] -= dly_;
    }

    const int later_snapshot_number_ = (snapshot_number_ < LastDarkMatterSnapShot) ? (snapshot_number_ + 1) : LastDarkMatterSnapShot;
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      tau_  = get_extinction(filter_number_, Z_g_, ZZ[later_snapshot_number_]) * n_h_sec_;
      alam_ = (tau_ > 0.0) ? (1.0 - exp(-tau_)) / tau_ : 1.;
      a_bc_ = 1. - exp(-tauvbc_ * pow((1. + ZZ[later_snapshot_number_]) * FilterLambda[filter_number_] * inv_VBand_WaveLength_, -0.7));

      Lum_disk_ = Gal[galaxy_number_].forward_ObsLum[output_number_][filter_number_] - Gal[galaxy_number_].forward_ObsLumBulge[output_number_][filter_number_];
      Gal[galaxy_number_].forward_ObsLumDust[output_number_][filter_number_] = Gal[galaxy_number_].forward_ObsLumBulge[output_number_][filter_number_] + Lum_disk_ * alam_;

      dly_ = (Gal[galaxy_number_].forward_ObsLumY[output_number_][filter_number_] - Gal[galaxy_number_].forward_ObsLumBulgeY[output_number_][filter_number_]) * alam_ * a_bc_ +
                      Gal[galaxy_number_].forward_ObsLumBulgeY[output_number_][filter_number_] * (1. - ExpTauBCBulge);

      Gal[galaxy_number_].forward_ObsLumDust[output_number_][filter_number_] -= dly_;
    }
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
  }
}
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
