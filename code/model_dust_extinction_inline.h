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
 *  You should have received a_ copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

#ifndef MODEL_DUST_EXTICTION_INLINE_H
#define MODEL_DUST_EXTICTION_INLINE_H

#include <math.h>

#include "allvars.h"
#include "proto.h"

/** @file model_dust_extiction_inline.h
 *  @brief used to compute dust extinction as described
 *         in Delucia2007 + redshift_ dependence as Kitzbichler & White 2007.

 *  There are 2 extinction sources:
 *  Extinction from a_ diffuse inter-stellar medium (ISM) (Devriendt1999);
 *  Extinction from molecular clouds in young stars (YS) (Charlot2000);
 *  Both were introduced in Delucia2007.

 *  The optical depth of dust in each component
 *  \f$\tau^z_{\lambda}\f$(ISM)
 *  and \f$\tau_{\lambda}^{\rm{BC}}\f$(YS)
 *  is used to compute extinction assuming a_ slab geometry for the dust and
 *  a_ random inclination of the disk to the line of sight.

 *  Extinction curves for the ISM:
 *  \f$\left(\frac{A_{\lambda}}{A_{\rm{v}}}\right)_{Z_{\odot}}
    \left(\frac{Z_{\rm{gas}}}{Z_{\odot}}\right)^s\f$
 *  are computed in get_extinction
 *
 *  The optical depth for the ISM at a_ given \f$\lambda\f$ can be written as:

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

#ifdef CARDELLI_DUST

/** @brief Dust extinction Cardelli et al. */
static inline double 
get_extinction(const int filter_number_, const double Z_g_, const double redshift_)
{
  static const float a_opt_[]={0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999};
  static const float b_opt_[]={1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002};
  static const float a_IR_=0.574;
  static const float b_IR_=-0.527;
  static const float a_UV_[]={1.752,-0.316,-0.104,-4.67,0.341};
  static const float b_UV_[]={-3.090,1.825,1.206,-4.62,0.263};
  static const float R_V_=3.1;

  int k_;
  float x_,y_,a_,b_;

  x_=1./FilterLambda[filter_number_];

  a_=0.0;
  b_=0.0;

  if(x_<=1.1)
    {
      a_=a_IR_*pow(x_,1.61f);
      b_=b_IR_*pow(x_,1.61f);
      //a_=a_IR_*pow(x_,-1.14f);
      //b_=b_IR_*pow(x_,-1.14f);
    }
  else if(x_>=1.1 && x_<=3.3)
    {
      y_=x_-1.82;
      x_=y_;
      a_=1.;
      b_=0.;
      for(k_=0;k_<7;k_++)
        {
          a_+=a_opt_[k_]*y_;
            b_+=b_opt_[k_]*y_;
            y_*=x_;
        }
    }
  else if(x_>=3.3 && x_<= 8.0)
    {
      a_=a_UV_[0]+a_UV_[1]*x_+a_UV_[2]/((x_+a_UV_[3])*(x_+a_UV_[3])+a_UV_[4]);
      b_=b_UV_[0]+b_UV_[1]*x_+b_UV_[2]/((x_+b_UV_[3])*(x_+b_UV_[3])+b_UV_[4]);
    }
  else if(x_>=8.0 && x_<=10.0)
    {
      a_=1.0;
      b_=13.415;
    }
  else if(x_>=10.0)
    {
      a_=1.0;
      b_=16.732;
    }

  A_Av_=(a_+b_/R_V_);

  if(ObsFrameLambda_ < 0.2)
  { A_Av_ = A_Av_ * pow(Z_g_, 1.35); }
  else
  { A_Av_ = A_Av_ * pow(Z_g_, 1.6 ) ; }
  
  return A_Av_;
}

#else /* not defined CARDELLI_DUST */

#ifdef MATHIS_N_
#if MATHIS_N_ != 42
#error "already defined MATHIS_N_ != 42"
#endif /* MATHIS_N_ != 42 */
#else  /* not defined MATHIS_N_ */
#define MATHIS_N_ 42
#endif /* not defined MATHIS_N_ */

static const float MathisLambda_[MATHIS_N_] = {0.091, 0.10, 0.13, 0.143, 0.18, 0.20, 0.21, 0.216, 0.23, 0.25,
    0.346, 0.435, 0.55, 0.7, 0.9, 1.2, 1.8, 2.2, 2.4, 3.4,
    4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 20.0, 25.0, 30.0, 40.0,
    50.0, 60.0, 70.0, 80.0, 100.0, 150.0, 200.0, 300.0, 400.0, 600.0,
    800.0, 1000.0};
static const float MathisAv_[MATHIS_N_] = { 5.720, 4.650, 2.960, 2.700, 2.490, 2.780, 3.000, 3.120, 2.860, 2.350,
    1.580, 1.320, 1.000, 0.750, 0.480, 0.280, 0.160, 0.122, 0.093, 0.038,
    0.024, 0.018, 0.014, 0.013, 0.072, 0.030, 0.065, 0.062, 0.032, 0.017,
    0.014, 0.012, 9.7e-3, 8.5e-3, 6.5e-3, 3.7e-3, 2.5e-3, 1.1e-3, 6.7e-4, 2.5e-4,
    1.4e-4, 7.3e-5};
static const float MathisAlbedo_[MATHIS_N_] = { 0.42, 0.43, 0.45, 0.45, 0.53, 0.56, 0.56, 0.56, 0.63, 0.63,
    0.71, 0.67, 0.63, 0.56, 0.50, 0.37, 0.25, 0.22, 0.15, 0.058,
    0.046, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00};


  /** @brief Dust extinction from Guiderdoni & Rocca-Volmerange 1987
   *         Mathis,Mezger&Panagia (1983) dust law is used */
static inline double 
get_extinction(const int filter_number_, const double Z_g_, const double redshift_)
{
//   FilterLambda[NMAG] = 0.55;      //to use by the dust model for birth clouds, the wavelength of the V-band
   /* this belongs in an init routine, 
   * thus was moved to setup_LumTables_precomputed()
   * and setup_Spec_LumTables_onthefly() */

  const float ObsFrameLambda_ = FilterLambda[filter_number_] / (1. + redshift_);

  /* If ObsFrameLambda_ less then the first MathisLambda_, higher then the last,
   * no interpolation is done. The first or value is taken respectively. */
  double A_Av_, Albedo_;
  if(ObsFrameLambda_ <= MathisLambda_[0])
  {
    A_Av_   = MathisAv_    [0];
    Albedo_ = MathisAlbedo_[0];
  }
  else if(ObsFrameLambda_ >=  MathisLambda_[MATHIS_N_ - 1])
  {
    A_Av_   = MathisAv_    [MATHIS_N_ - 1];
    Albedo_ = MathisAlbedo_[MATHIS_N_ - 1];
  }
  else
  {
   /** @bug (corrected by Stefan Hilbert)
   *       looks like in case of interpolation, i_ ended up after the while loop 
   *       such that ObsFrameLambda_ < MathisLambda_[i_], but interpolation
   *       assumes MathisLambda_[i_] <= ObsFrameLambda_ <= MathisLambda_[i_ + 1].
   *       this was the old version:
   *       int i_ = 0; while(ObsFrameLambda_ > MathisLambda_[i_] && i_ < MathisN) i_++;
   */
    int i_ = 0;
    while((i_ < (MATHIS_N_ - 1)) && (MathisLambda_[i_ + 1] < ObsFrameLambda_)) { i_++; }
    
    const double f_ = (ObsFrameLambda_ - MathisLambda_[i_]) / (MathisLambda_[i_ + 1] - MathisLambda_[i_]);
    A_Av_   = (1. - f_) * MathisAv_    [i_] + f_ * MathisAv_    [i_ + 1];
    Albedo_ = (1. - f_) * MathisAlbedo_[i_] + f_ * MathisAlbedo_[i_ + 1];
  }

  if(ObsFrameLambda_ < 0.2)
  { A_Av_ *= pow(Z_g_, 1.35) * sqrt(1. - Albedo_); }
  else
  { A_Av_ *= pow(Z_g_, 1.6 ) * sqrt(1. - Albedo_); }

  return A_Av_;
}

#undef MATHIS_N_

#endif /* not defined CARDELLI_DUST */

#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

#endif /* header guard */
