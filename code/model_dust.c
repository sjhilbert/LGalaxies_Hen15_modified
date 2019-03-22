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
void dust_model(const int p, const int snap, const int halonr)
{
  double nh, tau, alam, sec, Lum_disk, cosinc, Zg;
  double tauv, taubc, tauvbc, mu, dly;
  int k;
  
  const double VBand_WaveLength = 0.55;

  if(Gal[p].ColdGas > 0.0)
    {

      /* 0.94 = 2.83/3. -> 3 to get scale lenght and 2.83 = 1.68^2 */
      nh = Gal[p].ColdGas / (M_PI * pow(Gal[p].GasDiskRadius * 0.94, 2) * 1.4);
      /* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (2.1 10^21 atoms/cm^2) */
      nh = nh / 3252.37;        // 3252.37 = 10^(3.5122) ... ha ha ! 

      /*redshift dependence */
      nh = nh * pow(1 + ZZ[Halo[halonr].SnapNum], -1.0);


      Gal[p].CosInclination = fabs(Gal[p].StellarSpin[2]) /
                                  sqrt(Gal[p].StellarSpin[0]*Gal[p].StellarSpin[0]+
                                       Gal[p].StellarSpin[1]*Gal[p].StellarSpin[1]+
                                       Gal[p].StellarSpin[2]*Gal[p].StellarSpin[2]);
      cosinc = Gal[p].CosInclination;
      if(cosinc < 0.2)
            cosinc = 0.2;                // minimum inclination ~80 degrees
      sec = 1.0 / cosinc;

      /* mu for YS extinction, given by a Gaussian with centre 0.3 (MUCENTER)
       * and width 0.2 (MUWIDTH), truncated at 0.1 and 1.  */
      do { mu_ = gsl_ran_gaussian(random_generator, MUWIDTH) + MUCENTER; }
      while (mu < 0.1 || mu > 1.0);

      Zg = metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas/0.02;


#ifdef OUTPUT_REST_MAGS

      tauv = get_extinction(NMAG, Zg, 0) * nh;

      for(k = 0; k < NMAG; k++)
        {
          tau = get_extinction(k, Zg, 0) * nh;
          tau = tau * sec;

          if(tau > 0.0)
                alam = (1.0 - exp(-tau)) / tau;
          else
                alam = 1.;

          Lum_disk = Gal[p].Lum[k][snap] - Gal[p].LumBulge[k][snap];
          Gal[p].LumDust[k][snap] = Gal[p].LumBulge[k][snap] + Lum_disk * alam;

          // now remove light from young stars absorbed by birth clouds
          tauvbc = tauv * (1. / mu - 1.);
          taubc = tauvbc * pow(FilterLambda[k] / VBand_WaveLength, -0.7);

          dly = (Gal[p].YLum[k][snap] - Gal[p].YLumBulge[k][snap]) * alam * (1. - exp(-taubc)) +
                          Gal[p].YLumBulge[k][snap] * (1. - ExpTauBCBulge);

          Gal[p].LumDust[k][snap] -= dly;
        }
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS

      tauv = get_extinction(NMAG, Zg, 0) * nh;

      for(k = 0; k < NMAG; k++)
        {
          tau = get_extinction(k, Zg, ZZ[ListOutputSnaps[snap]]) * nh;
          tau = tau * sec;
          if(tau > 0.0)
                alam = (1.0 - exp(-tau)) / tau;
          else
                alam = 1.;

          Lum_disk = Gal[p].ObsLum[k][snap] - Gal[p].ObsLumBulge[k][snap];
          Gal[p].ObsLumDust[k][snap] = Gal[p].ObsLumBulge[k][snap] + Lum_disk * alam;

          // now remove light from young stars absorbed by birth clouds
          tauvbc = tauv * (1. / mu - 1.);
          taubc = tauvbc * pow((FilterLambda[k] * (1. + ZZ[ListOutputSnaps[snap]])) / VBand_WaveLength, -0.7);

          dly = (Gal[p].ObsYLum[k][snap] - Gal[p].ObsYLumBulge[k][snap]) * alam * (1. - exp(-taubc)) +
                          Gal[p].ObsYLumBulge[k][snap] * (1. - ExpTauBCBulge);

          Gal[p].ObsLumDust[k][snap] -= dly;


#ifdef OUTPUT_MOMAF_INPUTS   // compute same thing at z + 1

          if(snap < (LastDarkMatterSnapShot+1) - 1)
                tau = get_extinction(k, Zg, ZZ[ListOutputSnaps[snap] + 1]) * nh;
          else
                tau = get_extinction(k, Zg, ZZ[ListOutputSnaps[snap]]) * nh;
          tau = tau * sec;
          if(tau > 0.0)
                alam = (1.0 - exp(-tau)) / tau;
          else
                alam = 1.;

          Lum_disk = Gal[p].dObsLum[k][snap] - Gal[p].dObsLumBulge[k][snap];
          Gal[p].dObsLumDust[k][snap] = Gal[p].dObsLumBulge[k][snap] + Lum_disk * alam;

          // now remove light from young stars absorbed by birth clouds
          if(snap < (LastDarkMatterSnapShot+1) - 1)
                taubc = tauvbc * pow((FilterLambda[k] * (1. + ZZ[ListOutputSnaps[snap] + 1])) / VBand_WaveLength, -0.7);
          else
                taubc = tauvbc * pow((FilterLambda[k] * (1. + ZZ[ListOutputSnaps[snap]])) / VBand_WaveLength, -0.7);

          dly = (Gal[p].dObsYLum[k][snap] - Gal[p].dObsYLumBulge[k][snap]) * alam * (1. - exp(-taubc)) +
                          Gal[p].dObsYLumBulge[k][snap] * (1. - ExpTauBCBulge);

          Gal[p].dObsLumDust[k][snap] -= dly;

#endif /* defined OUTPUT_MOMAF_INPUTS */

        }//end for loop on mags (k)

#endif /* defined OUTPUT_OBS_MAGS */
    }
}
#endif /* not defined POST_PROCESS_MAGS */

#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
