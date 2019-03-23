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
 *  You should have received a_ copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

 /** @file   scale_cosmology.c
  *  @date   2010-2019
  *  @author Qi Guo (intitial version)
  *  @author Bruno Henriques
  *  @author Stefan Hilbert
  * 
  *  @brief  functions used to scale to a different cosmology
  */

#include <math.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

#ifdef MCMC
#include "mcmc_vars.h"
#include "mcmc_proto.h"
#endif


/** When we have tables for the scaling parameter in any cosmology only the
 * cosmological parameters will be inputs. Then this function will read the
 * scaling parameters from the tables. */
void read_scaling_parameters(void)
{
  const double om_min_=0.1;
  const double om_max_=0.6;
  const double s8_min_=0.5;
  const double s8_max_=1.0;

  const double om_binsize_=(om_max_-om_min_)/50;
  const double s8_binsize_=(s8_max_-s8_min_)/50;
//0.369798725 0.589

  Omega=0.370;
  const int om_n_bins_=(int)((0.370-om_min_)/om_binsize_);
  const int s8_n_bins_=(int)((0.589-s8_min_)/s8_binsize_);  

  double dummy_growth_, dummy_snap63_;
  char file_name_[1000], buffer_[1000];
  FILE *file_;

  sprintf(FileWithZList, "/galformod/scratch/bmh20/Workspace/CosmologyTables/zlist_%04d_%04d.txt", om_n_bins_, s8_n_bins_);

  read_zlist_new();
  read_output_snaps();

  sprintf(file_name_, "/galformod/scratch/bmh20/Workspace/CosmologyTables/fit_%04d_%04d.txt", om_n_bins_, s8_n_bins_);
  if(!(file_ = fopen(file_name_, "r")))
  {
    char error_message_[2048];
    sprintf(error_message_, "file `%s' not found.\n", file_name_);
    terminate(error_message_);
  }

  fgets(buffer_, 300, file_);
  fgets(buffer_, 300, file_);
  fgets(buffer_, 300, file_);

  if(fscanf(file_, "%lf %lf %lf %lf", &ScaleMass, &dummy_growth_, &ScalePos, &dummy_snap63_)!=4)
  {
    char error_message_[2048];
    sprintf(error_message_, "Wrong format of values in %s.\n", file_name_);
    terminate(error_message_);
  }

  fclose(file_);

  ScaleMass = 1./ScaleMass;
  ScalePos  = 1./ScalePos;

  PartMass  = PartMass_OriginalCosm * ScaleMass;
  BoxSize   = BoxSize_OriginalCosm  * ScalePos;

  printf("Boxsize=%f\n",BoxSize);
}


/** @brief scales all halos to new cosmology */
void scale_cosmology(const int n_halos_)
{
  int halo_number_, j_;
  double Scale_V_, Cen_Vel_[3], dv_;

#ifdef ALLOW_UNSCALE_COSMOLOGY
  //Save unscaled properties
  for(halo_number_ = 0; halo_number_ < n_halos_ ; halo_number_++)
  {
    //will make sure haloes in the future are not scaled/un_scaled
    if(Halo[halo_number_].SnapNum<=LastSnapShotNr)
    {
      HaloAux[halo_number_].M_Crit200_Unscaled = Halo[halo_number_].M_Crit200;
      HaloAux[halo_number_].M_Mean200_Unscaled = Halo[halo_number_].M_Mean200;
      HaloAux[halo_number_].Vmax_Unscaled = Halo[halo_number_].Vmax;
      for (j_ = 0; j_ < 3 ; j_++)
      {
        HaloAux[halo_number_].Pos_Unscaled [j_] = Halo[halo_number_].Pos [j_];
        HaloAux[halo_number_].Vel_Unscaled [j_] = Halo[halo_number_].Vel [j_];
        HaloAux[halo_number_].Spin_Unscaled[j_] = Halo[halo_number_].Spin[j_];
      }
    }
  }
#endif /* defined ALLOW_UNSCALE_COSMOLOGY */

  for (halo_number_ = 0; halo_number_ < n_halos_ ; halo_number_++)
  {
    Scale_V_ = scale_v_cen(Halo[Halo[halo_number_].FirstHaloInFOFgroup].SnapNum);

    //will make sure haloes in the future are not scaled/un_scaled
    if(Halo[halo_number_].SnapNum<=LastSnapShotNr)
    {
      if(Halo[halo_number_].M_Crit200 > 1.e-8)
        Halo[halo_number_].M_Crit200 = Halo[halo_number_].M_Crit200 * ScaleMass * c_correction(Halo[halo_number_].M_Crit200,Halo[halo_number_].SnapNum);
      if(Halo[halo_number_].M_Mean200 > 1.e-8)
        Halo[halo_number_].M_Mean200 = Halo[halo_number_].M_Mean200 * ScaleMass * c_correction(Halo[halo_number_].M_Mean200,Halo[halo_number_].SnapNum);
      Halo[halo_number_].Vmax = Halo[halo_number_].Vmax * sqrt(ScaleMass/ScalePos) * sqrt(AA_OriginalCosm[Halo[halo_number_].SnapNum]/AA[Halo[halo_number_].SnapNum]);

      for (j_ = 0; j_ < 3 ; j_++)
      {
        Halo[halo_number_].Pos[j_] = Halo[halo_number_].Pos[j_] * ScalePos;
        Halo[halo_number_].Spin[j_] *= ScalePos * sqrt(ScaleMass/ScalePos) * sqrt(AA[Halo[halo_number_].SnapNum]/AA_OriginalCosm[Halo[halo_number_].SnapNum]);

        Cen_Vel_[j_] = Halo[Halo[halo_number_].FirstHaloInFOFgroup].Vel[j_] * Scale_V_ ;
        if(halo_number_ !=  Halo[halo_number_].FirstHaloInFOFgroup) // subhalos
        {
          dv_ = Halo[halo_number_].Vel[j_] - Halo[Halo[halo_number_].FirstHaloInFOFgroup].Vel[j_];
          dv_ *=sqrt(ScaleMass/ScalePos) * sqrt(AA_OriginalCosm[Halo[halo_number_].SnapNum]/AA[Halo[halo_number_].SnapNum]);
          Halo[halo_number_].Vel[j_] = Cen_Vel_[j_] + dv_;
        }
        else //central halos
          Halo[halo_number_].Vel[j_] = Halo[halo_number_].Vel[j_] * Scale_V_  ;
      }
    }
  }
}

#ifdef ALLOW_UNSCALE_COSMOLOGY
/** @brief scales all halos back to old cosmology */
void un_scale_cosmology(const int n_halos_)
{
  int halo_number_, j_;

  for(halo_number_ = 0; halo_number_ < n_halos_ ; halo_number_++)
  {
    //will make sure haloes in the future are not scaled/un_scaled
    if(Halo[halo_number_].SnapNum<=LastSnapShotNr)
    {
      Halo[halo_number_].M_Crit200 = HaloAux[halo_number_].M_Crit200_Unscaled;
      Halo[halo_number_].M_Mean200 = HaloAux[halo_number_].M_Mean200_Unscaled;
      Halo[halo_number_].Vmax = HaloAux[halo_number_].Vmax_Unscaled;

      for (j_ = 0; j_ < 3 ; j_++)
      {
        Halo[halo_number_].Pos [j_] = HaloAux[halo_number_].Pos_Unscaled [j_];
        Halo[halo_number_].Vel [j_] = HaloAux[halo_number_].Vel_Unscaled [j_];
        Halo[halo_number_].Spin[j_] = HaloAux[halo_number_].Spin_Unscaled[j_];
      }
    }
  }
}
#endif /* defined ALLOW_UNSCALE_COSMOLOGY */


static inline double 
func_c(const double c_)
{ return log(1 + c_) - c_ / (1 + c_); }

static inline double 
func_c_p(const double c_)
{ return (log(1 + c_) - c_ / (1 + c_)) / (c_ * c_ * c_); }

/** @brief finds new concentration parameter
 * 
 * finds c_ using bisection
 * 
 * since the original version of this function showed up
 * surprisingly high on profile, a more optimized bisection
 * version was implemented
 * 
 * @warning assumes that initial values bracket the result 
 */
static inline double
find_c(const double c_ori_, const double ratio_)
{
  const double constant_ = ratio_ * func_c_p(c_ori_);
  
  double x_1_ = 1.;
  double x_2_ = 50.;
  
  double f_1_ = func_c_p(x_1_) - constant_;
  double x_m_, f_m_;
  while(
    x_m_ = 0.5 * (x_1_ + x_2_),
    f_m_ = func_c_p(x_m_) - constant_,
    fabs(f_m_) > 1.e-8
  ) 
  {
    if(f_m_ * f_1_ > 0)
    {       
      x_1_ = x_m_;
      f_1_ = f_m_;
    }
    else
    { x_2_ = x_m_; }
  }
  return x_m_;
}


/** @brief computes c_ correction */
double c_correction(const float halo_mass_, const int snapshot_number_)
{
  double c_original_, c_new_, Omega_new_, Omega_original_, ratio_;

  c_original_ = 5 * pow(0.0001 * halo_mass_, -0.1);
  
  Omega_new_ = Omega * 1./pow3(AA[snapshot_number_]) / 
              (Omega * 1./pow3(AA[snapshot_number_]) + OmegaLambda);
  
  Omega_original_ = Omega_OriginalCosm * 1./pow3(AA_OriginalCosm[snapshot_number_]) /
                   (Omega_OriginalCosm * 1./pow3(AA_OriginalCosm[snapshot_number_]) + OmegaLambda_OriginalCosm);
                  
  ratio_ = Omega_original_/ Omega_new_;
  
  c_new_ = find_c(c_original_, ratio_);

  return func_c(c_new_) / func_c(c_original_);
}


double dgrowth_factor_dt(const double a_, const double omega_m_, const double omega_l_)
{
  // const double g0 = 2.5 * omega_m_ / (pow(omega_m_, 4./7.) - omega_l_ + (1.0 + 0.5 * omega_m_) * (1.0 + (1./70.) *omega_l_));
  const double inv_g0_ = (pow(omega_m_, 4./7.) - omega_l_ + (1.0 + 0.5 * omega_m_) * (1.0 + (1./70.) *omega_l_)) / (2.5 * omega_m_);
  
  // const double o_m_ = omega_m_ * 1./pow3(a_)/(omega_m_ * 1./pow3(a_) + omega_l_);
  // const double o_l_ = omega_l_ / (omega_m_ * 1./pow3(a_) + omega_l_);
  // const double do_m = -3 * omega_m_ * omega_l_ / (a_ * a_ * a_ * a_) / (omega_m_ / (a_ * a_ * a_) + omega_l_) / (omega_m_ / (a_ * a_ * a_) + omega_l_);
  // const double do_l = -do_m;

  // const double o_m_ =                         omega_m_ * 1. /     (omega_m_ + pow3(a_) * omega_l_);
  // const double o_l_ =               pow3(a_) * omega_l_ * 1. /     (omega_m_ + pow3(a_) * omega_l_);
  // const double do_m = -3 * a_ * a_ * omega_m_ * omega_l_ * 1. / pow2(omega_m_ + pow3(a_) * omega_l_);
  // const double do_l = -do_m;

  const double inv_omega_m_plus_a_a_a_omega_l_ = 1. / (omega_m_ + pow3(a_) * omega_l_);
  const double o_m_ =           omega_m_ * inv_omega_m_plus_a_a_a_omega_l_;
  const double o_l_ = pow3(a_) * omega_l_ * inv_omega_m_plus_a_a_a_omega_l_;
  // const double do_m = -3 * o_m_ * o_l_ / a_;
  // const double do_l = -do_m;

  // const double hubble_a = sqrt(omega_m_ / pow3(a_) + omega_l_);

  //da_dtau = sqrt(1 + o_m_ * (1 / a_ - 1) + o_l_ * (a_ * a_ - 1));   //tau = H0*t
  // const double extra_fac = - ( 4/7.* pow(o_m_, -3./7) * do_m - do_l
  //                 +(do_m /2. *(1 + (1./70.) * o_l_) - (1 + o_m_ / 2.) * do_l / 70)
  //                   /(1 + (1./70.) * o_l_)/(1 + o_l_ / 70))/(pow(o_m_, 4./7.) - o_l_ + (1.0 + 0.5*o_m_)*(1.0 + (1./70.) * o_l_))/(pow(o_m_, 4./7.) - o_l_ + (1.0 + 0.5*o_m_)*(1.0 + (1./70.) * o_l_));
  
  // const double extra_fac = -((4./7.)* pow(o_m_, -3./7.) * do_m - do_l  + (0.5 * do_m * (1.0 + (1./70.) * o_l_) - (1.0 + 0.5 * o_m_) * (1./70.) * do_l) / pow2(1.0 + (1./70.) * o_l_))
  //                              * 1. / pow2(pow(o_m_, 4./7.) - o_l_ + (1.0 + 0.5 * o_m_) * (1.0 + (1./70.) * o_l_));
  // const double g_  = 2.5 *  o_m_ * 1. /     (pow(o_m_, 4./7.) - o_l_ + (1.0 + 0.5 * o_m_) * (1.0 + (1./70.) * o_l_));
  // const double dg = 2.5 * do_m * 1. /     (pow(o_m_, 4./7.) - o_l_ + (1.0 + 0.5 * o_m_) * (1.0 + (1./70.) * o_l_)) + 2.5 * o_m_ * extra_fac;
  
  // const double pow_o_m_4_7_ = pow(o_m_, 4./7.);
  // const double den_         = 1. / (pow_o_m_4_7_ - o_l_ + (1.0 + 0.5 * o_m_) * (1.0 + (1./70.) * o_l_));
  // const double extra_fac = -((4./7.) * pow_o_m_4_7_ / o_m_ * do_m - do_l  + (0.5 * do_m * (1.0 + (1./70.) * o_l_) - (1.0 + 0.5 * o_m_) * (1./70.) * do_l) / pow2(1.0 + (1./70.) * o_l_)) * pow2(den_);
  // const double g_  = 2.5 *  o_m_ * den_ ;
  // const double dg = 2.5 * do_m * den_ + 2.5 * o_m_ * extra_fac;
  
  const double pow_o_m_4_7_ = pow(o_m_, 4./7.);
  const double den_         = 1. / (pow_o_m_4_7_ - o_l_ + (1.0 + 0.5 * o_m_) * (1.0 + (1./70.) * o_l_));
  const double g_           =  2.5 * o_m_ * den_ ;
  const double dg_a_        = -7.5 * o_m_ * den_ * o_l_ * ( 1. - den_ * ((4./7.) * pow_o_m_4_7_ + o_m_ +  ((0.5 + (1./70.)) + (0.5/70.) * (o_l_ + o_m_)) * o_m_ / pow2(1.0 + (1./70.) * o_l_)));
  
  //  const double dD_dt_ = hubble_a * a_ * (dg * a_ +  g_ ) / g0;
  const double dD_dt_ = sqrt(omega_m_ / pow3(a_) + omega_l_) * a_ * inv_g0_ * (dg_a_ + g_);
  return dD_dt_;
}


double scale_v_cen(const int snapshot_number_)
{
  return ScalePos * dgrowth_factor_dt(AA[snapshot_number_],Omega, OmegaLambda) /
                        dgrowth_factor_dt(AA_OriginalCosm[snapshot_number_],Omega_OriginalCosm,OmegaLambda_OriginalCosm) *
                         AA[snapshot_number_]/AA_OriginalCosm[snapshot_number_] * Hubble_h / Hubble_h_OriginalCosm;
}

