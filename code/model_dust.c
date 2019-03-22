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

#include "allvars.h"
#include "proto.h"

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


/** @brief hand made rng (replace by gsl version?) */
static inline float 
ran1(long *idum)
{
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
  
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if(*idum <= 0 || !iy)
    {
      if(-(*idum) < 1)
                *idum = 1;
      else
                *idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--)
                {
                  k = (*idum) / IQ;
                  *idum = IA * (*idum - k * IQ) - IR * k;
                  if(*idum < 0)
                    *idum += IM;
                  if(j < NTAB)
                    iv[j] = *idum;
                }
      iy = iv[0];
    }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if(*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
  
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX  
}


/** @brief computes a gaussian random deviate to calculate a random
 *         inclination for extinction. */
float gasdev(long *idum)
{
  static int iset = 0;
  static float gset;
  float fac, rsq, v1, v2;

  if(iset == 0)
    {
      do
        {
          v1 = 2.0 * ran1(idum) - 1.0;
          v2 = 2.0 * ran1(idum) - 1.0;
          rsq = v1 * v1 + v2 * v2;
        }
      while(rsq >= 1.0 || rsq == 0.0);
      fac = sqrt(-2.0 * log(rsq) / rsq);
      gset = v1 * fac;
      iset = 1;
      return v2 * fac;
    }
  else
    {
      iset = 0;
      return gset;
    }
}


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
      do { mu = gasdev(&mu_seed) * MUWIDTH + MUCENTER; } while (mu < 0.1 || mu > 1.0);

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


#ifndef CARDELLI_DUST
  /** @brief Dust extinction from Guiderdoni & Rocca-Volmerange 1987
   *         Mathis,Mezger&Panagia (1983) dust law is used */
double get_extinction(const int mag, const double Zg, const double redshift)
{
#ifdef MATHIS_N
#if MATHIS_N != 42
#error "already defined MATHIS_N != 42"
#endif /* MATHIS_N != 42 */
#else  /* not defined MATHIS_N */
#define MATHIS_N 42
#endif /* not defined MATHIS_N */

  static const float MathisLambda[MATHIS_N] = {0.091, 0.10, 0.13, 0.143, 0.18, 0.20, 0.21, 0.216, 0.23, 0.25,
    0.346, 0.435, 0.55, 0.7, 0.9, 1.2, 1.8, 2.2, 2.4, 3.4,
    4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 20.0, 25.0, 30.0, 40.0,
    50.0, 60.0, 70.0, 80.0, 100.0, 150.0, 200.0, 300.0, 400.0, 600.0,
    800.0, 1000.0
  };
  static const float MathisAv[MATHIS_N] = { 5.720, 4.650, 2.960, 2.700, 2.490, 2.780, 3.000, 3.120, 2.860, 2.350,
    1.580, 1.320, 1.000, 0.750, 0.480, 0.280, 0.160, 0.122, 0.093, 0.038,
    0.024, 0.018, 0.014, 0.013, 0.072, 0.030, 0.065, 0.062, 0.032, 0.017,
    0.014, 0.012, 9.7e-3, 8.5e-3, 6.5e-3, 3.7e-3, 2.5e-3, 1.1e-3, 6.7e-4, 2.5e-4,
    1.4e-4, 7.3e-5
  };
  static const float MathisAlbedo[MATHIS_N] = { 0.42, 0.43, 0.45, 0.45, 0.53, 0.56, 0.56, 0.56, 0.63, 0.63,
    0.71, 0.67, 0.63, 0.56, 0.50, 0.37, 0.25, 0.22, 0.15, 0.058,
    0.046, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00
  };
  
//   FilterLambda[NMAG] = 0.55;      //to use by the dust model for birth clouds, the wavelength of the V-band
   /* this belongs in an init routine, 
   * thus was moved to setup_LumTables_precomputed()
   * and setup_Spec_LumTables_onthefly() */

  const float ObsFrameLambda = FilterLambda[mag] / (1. + redshift);

  /* If ObsFrameLambda less then the first MathisLambda, higher then the last,
   * no interpolation is done. The first or value is taken respectively. */
  double A_Av, Albedo;
  if(ObsFrameLambda <= MathisLambda[0])
  {
    A_Av   = MathisAv    [0];
    Albedo = MathisAlbedo[0];
  }
  else if(ObsFrameLambda >=  MathisLambda[MATHIS_N - 1])
  {
    A_Av   = MathisAv    [MATHIS_N - 1];
    Albedo = MathisAlbedo[MATHIS_N - 1];
  }
  else
  {
   /** @bug (corrected by Stefan Hilbert)
   *       looks like in case of interpolation, i ended up after the while loop 
   *       such that ObsFrameLambda < MathisLambda[i], but interpolation
   *       assumes MathisLambda[i] <= ObsFrameLambda <= MathisLambda[i + 1].
   *       this was the old version:
   *       int i = 0; while(ObsFrameLambda > MathisLambda[i] && i < MathisN) i++;
   */
    int i = 0;
    while((i < (MATHIS_N - 1)) && (MathisLambda[i + 1] < ObsFrameLambda)) { i++; }
    
    const double f_ = (ObsFrameLambda - MathisLambda[i]) / (MathisLambda[i + 1] - MathisLambda[i]);
    A_Av   = (1. - f_) * MathisAv    [i] + f_ * MathisAv    [i + 1];
    Albedo = (1. - f_) * MathisAlbedo[i] + f_ * MathisAlbedo[i + 1];
  }

  if(ObsFrameLambda < 0.2)
  { A_Av *= pow(Zg, 1.35) * pow(1. - Albedo, 0.5); }
  else
  { A_Av *= pow(Zg, 1.6 ) * pow(1. - Albedo, 0.5); }

  return A_Av;
}


#else /* defined CARDELLI_DUST */
/** @brief Dust extinction Cardelli et al. */
double get_extinction(const int mag, const double Zg, const double redshift)
{
  static const float a_opt[]={0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999};
  static const float b_opt[]={1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002};
  static const float a_IR=0.574;
  static const float b_IR=-0.527;
  static const float a_UV[]={1.752,-0.316,-0.104,-4.67,0.341};
  static const float b_UV[]={-3.090,1.825,1.206,-4.62,0.263};
  static const float R_V=3.1;

  int k;
  float x,y,a,b;

  x=1./FilterLambda[mag];

  a=0.0;
  b=0.0;

  if(x<=1.1)
    {
      a=a_IR*pow(x,1.61f);
      b=b_IR*pow(x,1.61f);
      //a=a_IR*pow(x,-1.14f);
      //b=b_IR*pow(x,-1.14f);
    }
  else if(x>=1.1 && x<=3.3)
    {
      y=x-1.82;
      x=y;
      a=1.;
      b=0.;
      for(k=0;k<7;k++)
        {
          a+=a_opt[k]*y;
            b+=b_opt[k]*y;
            y*=x;
        }
    }
  else if(x>=3.3 && x<= 8.0)
    {
      a=a_UV[0]+a_UV[1]*x+a_UV[2]/((x+a_UV[3])*(x+a_UV[3])+a_UV[4]);
      b=b_UV[0]+b_UV[1]*x+b_UV[2]/((x+b_UV[3])*(x+b_UV[3])+b_UV[4]);
    }
  else if(x>=8.0 && x<=10.0)
    {
      a=1.0;
      b=13.415;
    }
  else if(x>=10.0)
    {
      a=1.0;
      b=16.732;
    }

  A_Av=(a+b/R_V);

  if(ObsFrameLambda < 0.2)
  { A_Av = A_Av * pow(Zg, 1.35); }
  else
  { A_Av = A_Av * pow(Zg, 1.6); }
  
  return A_Av;
}
#endif /* defined CARDELLI_DUST */

#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

