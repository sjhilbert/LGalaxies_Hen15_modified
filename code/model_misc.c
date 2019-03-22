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
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

/**@file model_misc.c
 * @brief model_misc.c contains a mix of recipes used to: calculate disk
 *        sizes, initiate a galaxy structure, update central galaxy,
 *        update type 1 and type2, transfer stars and gas between galaxies,
 *        etc. */


/** @brief Calculates the separation of galaxies galaxy_number_a_ and galaxy_number_b_, allowing for wrapping */
double separation_gal(const int galaxy_number_a_, const int galaxy_number_b_)
 {
  int i_;
  double sep_1_,sep_squared_;

  sep_squared_=0.;
  for (i_=0; i_<3; i_++)
    {
      sep_1_ =  Gal[galaxy_number_a_].Pos[i_] - Gal[galaxy_number_b_].Pos[i_];
      sep_1_ =  wrap(sep_1_, BoxSize);
      sep_squared_+= sep_1_*sep_1_;
    }
  return sqrt(sep_squared_);
}


/** @brief Calculates the separation of galaxies halo_number_a_ and halo_number_b_, allowing for wrapping */
double separation_halo(const int halo_number_a_, const  int halo_number_b_) 
{
  int i_;
  double sep_1_, sep_squared_;

  sep_squared_=0.;
  for (i_=0; i_<3; i_++)
    {
          sep_1_ = Halo[halo_number_a_].Pos[i_] - Halo[halo_number_b_].Pos[i_];
          sep_1_ = wrap(sep_1_,BoxSize);
          sep_squared_+=sep_1_*sep_1_;
    }
  return sqrt(sep_squared_);
}


/** @brief computes disk radius
 * 
 *        computes disk radius from halo spin
 * 
 * @note This routine is no longer called (after Guo2010) */
double get_disk_radius(const int halo_number_, const int galaxy_number_)
{
  if(DiskRadiusModel == 1)
  {
    /*  See Mo, Mao & White (1998) eq12, and using a Bullock style lambda.  Since this is the scale length
        we take the typical star forming region as 3 times this using the Milky Way as an approximate guide */
    return 1.5 * sqrt(Halo[halo_number_].Spin[0] * Halo[halo_number_].Spin[0] + 
                      Halo[halo_number_].Spin[1] * Halo[halo_number_].Spin[1] +
                      Halo[halo_number_].Spin[2] * Halo[halo_number_].Spin[2]) / Gal[galaxy_number_].Vvir;
  }
  else
    /*  simple prescription */
    return 0.1 * Gal[galaxy_number_].Rvir;
}


/** @brief Initiates the value of the disk radius.
 *
 *  First determination of radius in Guo2010 (same as in Delucia2007), after this,
 *  the disks are updated using set_gas_disk_radius and set_stellar_disk_radius.
 *  Two options are available:
 *
 *    If DiskRadiusModel = 2 then \f$R_{\rm{disk}}=\frac{R_{\rm{vir}}}{10}\f$
 *
 *    If DiskRadiusModel = 0 or 1 then the Mo, Mao & White (1998) formalism is
 *    used with a Bullock style \f$\lambda\f$:
 *
 *    \f$ R_d=\frac{1}{\sqrt{2}}\frac{j_d}{m_d}\lambda r_{200}\f$
 *
 *    and using the Milky Way as an approximate guide \f$R_{\rm{disk}}=3R_d\f$. 
 * 
 * @bug possible bug: should really test for GasSpin != 0 instead of HaloSpin != 0?
 */
double get_initial_disk_radius(const int halo_number_, const int galaxy_number_)
{
  if(DiskRadiusModel == 0 || DiskRadiusModel == 1)
  {
    const double spin_squared_ = Halo[halo_number_].Spin[0] * Halo[halo_number_].Spin[0] +
                                 Halo[halo_number_].Spin[1] * Halo[halo_number_].Spin[1] +
                                 Halo[halo_number_].Spin[2] * Halo[halo_number_].Spin[2];  

    if(Gal[galaxy_number_].GasSpin[0]==0 && Gal[galaxy_number_].GasSpin[1]==0 && Gal[galaxy_number_].GasSpin[2]==0)
      return 0.1 * Gal[galaxy_number_].Rvir;
    else
      return 1.5 * sqrt(spin_squared_) / Gal[galaxy_number_].Vvir;
  }
  else
    /*  simpler prescription */
    return 0.1 * Gal[galaxy_number_].Rvir;
}


/** @brief Calculates the virial mass: \f$M_{\rm{crit200}}\f$ for central halos
 *        with \f$M_{\rm{crit200}}\f$ or Len*PartMass for central halos without. */
double get_virial_mass(const int halo_number_)
{
  if(halo_number_ == Halo[halo_number_].FirstHaloInFOFgroup && Halo[halo_number_].M_Crit200)
    return Halo[halo_number_].M_Crit200;        /* take spherical overdensity mass estimate */
  else
    return Halo[halo_number_].Len * PartMass;
}


/** @brief Calculates the virial velocity from the virial mass and radius.
 *
 * Calculates virial velocity:
 *    \f$ V_{\rm{vir}}=\frac{GM_{\rm{vir}}}{R_{\rm{vir}}} \f$*/

double get_virial_velocity(const int halo_number_)
{
  return sqrt(Gravity * get_virial_mass(halo_number_) / get_virial_radius(halo_number_));
}


/** @brief computes Hubbble parameter for halo redshift  */
double get_hubble_parameter_for_halo(const int halo_number_)
{
  double zplus1;

  zplus1 = 1 + ZZ[Halo[halo_number_].SnapNum];

  /*get H for current z*/
  return Hubble * sqrt(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 + OmegaLambda);
}


/** @brief computes Hubbble parameter for halo redshift  */
double get_hubble_parameter_squared_for_halo(const int halo_number_)
{
  const double zplus1 = 1 + ZZ[Halo[halo_number_].SnapNum];

  return Hubble *  Hubble * (Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 + OmegaLambda);
}


/** @brief Calculates virial radius from a critical overdensity
 *
 * Calculates virial radius using:
 * \f$ M_{\rm{vir}}=\frac{4}{3}\pi R_{\rm{vir}}^3 \rho_c \Delta_c\f$.
 *
 * From which, assuming \f$ \Delta_c=200\f$, *
 * \f$ R_{\rm{vir}}=\left( \frac{3M_{\rm{vir}}}{4\pi 200 \rho_c}\right)^{1/3}\f$
 */
double get_virial_radius(const int halo_number_)
// { return pow(get_virial_mass(halo_number_) / get_hubble_parameter_squared_for_halo(halo_number_) * (0.01 * Gravity), 1. / 3.); }
{ return cbrt(get_virial_mass(halo_number_) / get_hubble_parameter_squared_for_halo(halo_number_) * (0.01 * Gravity)); }


/** @brief Updates the gas disk radius.
 *
 *  The gas disk is assumed to be thin, in centrifugal equilibrium and to have
 *  an exponential density profile:
 *
 *  \f$ \Sigma(R_{\rm{gas}})=
 *      \Sigma_{\rm{gas0}}e^{-\frac{R_{\rm{gas}}}{R_{\rm{gas,d}}}}, \f$
 *
 *  where \f$R_{\rm{gas,d}}\f$ is the scale length of the gas disk and
 *  \f$\Sigma_{\rm{gas0}}\f$ is the corresponding central surface density.
 *
 *  Assuming a flat circular velocity curve (galaxy with a negligible self-gravity)
 *  in an isothermal dark matter halo, the gas disk scale length is given by:
 *
 *  \f$ R_{\rm{gas,d}}=\frac{J_{\rm{gas}}/M_{\rm{gas}}}{2V_{\rm{max}}}\f$,
 *
 *  assuming conservation of the angular momentum of the cooling gas and that the
 *  maximum circular velocity of satellite galaxies does not change after infall
 *  (inner dark matter regions are compact and don't change). */
void set_gas_disk_radius(const int galaxy_number_)
{
  const double spin_squared_ = Gal[galaxy_number_].GasSpin[0] * Gal[galaxy_number_].GasSpin[0] +
                               Gal[galaxy_number_].GasSpin[1] * Gal[galaxy_number_].GasSpin[1] +
                               Gal[galaxy_number_].GasSpin[2] * Gal[galaxy_number_].GasSpin[2];
  if(spin_squared_ > 0.)
  {
    Gal[galaxy_number_].GasDiskRadius = (Gal[galaxy_number_].Type == 0) ? 1.5 * sqrt(spin_squared_) / Gal[galaxy_number_].Vmax  :
                                                                          1.5 * sqrt(spin_squared_) / Gal[galaxy_number_].InfallVmax;
  }
  else
  { Gal[galaxy_number_].GasDiskRadius = 0.1 * Gal[galaxy_number_].Rvir;  }
}


/** @brief Updates the stellar disk radius.
 *
 *  The stellar disk is assumed to be thin, in centrifugal equilibrium and to have
 *  an exponential density profile:
 *
 *  \f$ \Sigma(R_{\star})=
 *      \Sigma_{\star\rm{0}}e^{-\frac{R_{\star}}{R_{\rm{\star,d}}}}, \f$
 *
 *  where \f$R_{\rm{\star,d}}\f$ is the scale length of the gas disk and
 *  \f$\Sigma_{\star\rm{0}}\f$ is the corresponding central surface density.
 *
 *  Assuming a flat circular velocity curve (galaxy with a negligible self-gravity)
 *  in an isothermal dark matter halo, the gas disk scale length is given by:
 *
 *  \f$ R_{\rm{\star,d}}=\frac{J_{\star}/M_{\star}}{2V_{\rm{max}}}\f$,
 *
 *  assuming that the maximum circular velocity of satellite galaxies does not
 *  change after infall (inner dark matter regions are compact and don't change). */
void set_stellar_disk_radius(const int galaxy_number_)
{
  const double spin_squared_ = Gal[galaxy_number_].StellarSpin[0] * Gal[galaxy_number_].StellarSpin[0] +
                               Gal[galaxy_number_].StellarSpin[1] * Gal[galaxy_number_].StellarSpin[1] +
                               Gal[galaxy_number_].StellarSpin[2] * Gal[galaxy_number_].StellarSpin[2];
  if(spin_squared_ > 0.)
  {
    Gal[galaxy_number_].StellarDiskRadius = (Gal[galaxy_number_].Type == 0) ? 1.5 * sqrt(spin_squared_) / Gal[galaxy_number_].Vmax : 
                                                                              1.5 * sqrt(spin_squared_) / Gal[galaxy_number_].InfallVmax;
  }
  else
  { Gal[galaxy_number_].StellarDiskRadius = 0.1 * Gal[galaxy_number_].Rvir;  }
}


/** @brief Initializes the Galaxy Structure by setting all its
 *         elements to zero. */
void init_galaxy(const int galaxy_number_, const int halo_number_)
{
  int j_, output_number_;
  
#ifndef POST_PROCESS_MAGS
  int filter_number_;
#endif /* not defined POST_PROCESS_MAGS */

  /* make explicitly sure that the whole galaxy structure has defined 0 values */
  memset(&Gal[galaxy_number_], 0, sizeof(struct GALAXY));

  Gal[galaxy_number_].NextGalaxy = -1;
#ifdef GALAXYTREE
  Gal[galaxy_number_].FirstProgGal = -1;
#endif /* defined GALAXYTREE */

  if(halo_number_ != Halo[halo_number_].FirstHaloInFOFgroup)
  {
    terminate("Hah?\n");
  }

  Gal[galaxy_number_].Type = 0;

  Gal[galaxy_number_].HaloNr = halo_number_;
  Gal[galaxy_number_].MostBoundID = Halo[halo_number_].MostBoundID;
  Gal[galaxy_number_].SnapNum = Halo[halo_number_].SnapNum - 1;
#ifdef HALOPROPERTIES
  Gal[galaxy_number_].HaloM_Mean200 = Halo[halo_number_].M_Mean200;
  Gal[galaxy_number_].HaloM_Crit200 = Halo[halo_number_].M_Crit200;
  Gal[galaxy_number_].HaloM_TopHat = Halo[halo_number_].M_TopHat;
  Gal[galaxy_number_].HaloVelDisp = Halo[halo_number_].VelDisp;
  Gal[galaxy_number_].HaloVmax = Halo[halo_number_].Vmax;
#endif /* defined HALOPROPERTIES */

  for(j_ = 0; j_ < 3; j_++)
  {
    Gal[galaxy_number_].Pos[j_] = Halo[halo_number_].Pos[j_];
    Gal[galaxy_number_].Vel[j_] = Halo[halo_number_].Vel[j_];
    Gal[galaxy_number_].GasSpin[j_] = Halo[halo_number_].Spin[j_];
    Gal[galaxy_number_].StellarSpin[j_] = Halo[halo_number_].Spin[j_];
    Gal[galaxy_number_].HaloSpin[j_] = Halo[halo_number_].Spin[j_];
    Gal[galaxy_number_].MergCentralPos[j_] = Gal[galaxy_number_].Pos[j_];
    Gal[galaxy_number_].DistanceToCentralGal[j_]=0.0;
#ifdef HALOPROPERTIES
    Gal[galaxy_number_].HaloPos[j_] = Halo[halo_number_].Pos[j_];
    Gal[galaxy_number_].HaloVel[j_] = Halo[halo_number_].Vel[j_];
#endif /* defined HALOPROPERTIES */
  }

  Gal[galaxy_number_].Len = Halo[halo_number_].Len;
  Gal[galaxy_number_].Vmax = Halo[halo_number_].Vmax;
  Gal[galaxy_number_].InfallVmax = Halo[halo_number_].Vmax;
  Gal[galaxy_number_].InfallVmaxPeak = Gal[galaxy_number_].InfallVmax;
  Gal[galaxy_number_].Vvir = get_virial_velocity(halo_number_);
  Gal[galaxy_number_].Mvir = get_virial_mass(halo_number_);
  Gal[galaxy_number_].Rvir = get_virial_radius(halo_number_);
  Gal[galaxy_number_].MergeSat = 0.0;
  Gal[galaxy_number_].InfallSnap = Halo[halo_number_].SnapNum;

  Gal[galaxy_number_].ColdGas = 0.0;
  Gal[galaxy_number_].DiskMass = 0.0;
  Gal[galaxy_number_].BulgeMass = 0.0;
  Gal[galaxy_number_].HotGas = 0.0;
  Gal[galaxy_number_].EjectedMass = 0.0;
  Gal[galaxy_number_].ICM = 0.0;
#ifdef TRACK_BURST
  Gal[galaxy_number_].BurstMass = 0.0;
#endif /* defined TRACK_BURST */

  Gal[galaxy_number_].BlackHoleMass = 0.0;
  /*ram pressure*/
  Gal[galaxy_number_].HotRadius=Gal[galaxy_number_].Rvir;
#ifdef GALAXYTREE
  Gal[galaxy_number_].DisruptOn = 0;
#endif /* defined GALAXYTREE */

#ifdef DETAILED_METALS_AND_MASS_RETURN
  Gal[galaxy_number_].MetalsColdGas = metals_init();
  Gal[galaxy_number_].MetalsDiskMass = metals_init();
  Gal[galaxy_number_].MetalsBulgeMass = metals_init();
  Gal[galaxy_number_].MetalsHotGas = metals_init();
  Gal[galaxy_number_].MetalsEjectedMass = metals_init();
  Gal[galaxy_number_].MetalsICM = metals_init();
#ifdef METALS_SELF
  Gal[galaxy_number_].MetalsHotGasSelf = metals_init();
#endif /* defined METALS_SELF */
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
  Gal[galaxy_number_].MetalsColdGas = 0.0;
  Gal[galaxy_number_].MetalsDiskMass = 0.0;
  Gal[galaxy_number_].MetalsBulgeMass = 0.0;
  Gal[galaxy_number_].MetalsHotGas = 0.0;
  Gal[galaxy_number_].MetalsEjectedMass = 0.0;
  Gal[galaxy_number_].MetalsICM = 0.0;
#ifdef METALS_SELF
  Gal[galaxy_number_].MetalsHotGasSelf = 0.0;
#endif /* defined METALS_SELF */
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */

  //inclination defined as the angle between galaxy spin and the z-axis
  Gal[galaxy_number_].CosInclination = 0.0;

  Gal[galaxy_number_].PrimordialAccretionRate = 0.0;
  Gal[galaxy_number_].CoolingRate = 0.0;
  Gal[galaxy_number_].CoolingRate_beforeAGN = 0.0;
  Gal[galaxy_number_].CoolingRadius = 0.0;
  Gal[galaxy_number_].CoolingGas = 0.0;
  Gal[galaxy_number_].QuasarAccretionRate=0.0;
  Gal[galaxy_number_].RadioAccretionRate=0.0;
  Gal[galaxy_number_].AGNheatingFromCentral = 0.0;

  Gal[galaxy_number_].Sfr = 0.0;
  Gal[galaxy_number_].SfrBulge = 0.0;

  Gal[galaxy_number_].StarMerge=0.0;

  Gal[galaxy_number_].XrayLum = 0.0;
  Gal[galaxy_number_].GasDiskRadius = get_initial_disk_radius(halo_number_, galaxy_number_);
  Gal[galaxy_number_].StellarDiskRadius = Gal[galaxy_number_].GasDiskRadius;
  Gal[galaxy_number_].BulgeSize = 0.0;

  Gal[galaxy_number_].OriMergTime = 0.0;
  Gal[galaxy_number_].MergTime = 0.0;

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    Gal[galaxy_number_].MassWeightAge[output_number_] = 0.0;

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    Gal[galaxy_number_].rbandWeightAge[output_number_] = 0.0;

#ifndef POST_PROCESS_MAGS
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
#ifdef OUTPUT_REST_MAGS
      Gal[galaxy_number_].Lum[output_number_][filter_number_]         = 0.0;
      Gal[galaxy_number_].LumY[output_number_][filter_number_]        = 0.0;
      Gal[galaxy_number_].LumBulge[output_number_][filter_number_]    = 0.0;
      Gal[galaxy_number_].LumBulgeY[output_number_][filter_number_]   = 0.0;
      Gal[galaxy_number_].LumDust[output_number_][filter_number_]     = 0.0;
#ifdef ICL
      Gal[galaxy_number_].LumICL[output_number_][filter_number_]      = 0.0;
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
      Gal[galaxy_number_].ObsLum[output_number_][filter_number_]        = 0.0;
      Gal[galaxy_number_].ObsLumY[output_number_][filter_number_]       = 0.0;
      Gal[galaxy_number_].ObsLumBulge[output_number_][filter_number_]   = 0.0;
      Gal[galaxy_number_].ObsLumBulgeY[output_number_][filter_number_]  = 0.0;
      Gal[galaxy_number_].ObsLumDust[output_number_][filter_number_]    = 0.0;
#ifdef ICL
      Gal[galaxy_number_].ObsLumICL[output_number_][filter_number_]     = 0.0;
#endif /* defined ICL */
          
#ifdef OUTPUT_FB_OBS_MAGS
      Gal[galaxy_number_].backward_ObsLum[output_number_][filter_number_]       = 0.0;
      Gal[galaxy_number_].backward_ObsLumY[output_number_][filter_number_]      = 0.0;
      Gal[galaxy_number_].backward_ObsLumBulge[output_number_][filter_number_]  = 0.0;
      Gal[galaxy_number_].backward_ObsLumBulgeY[output_number_][filter_number_] = 0.0;
      Gal[galaxy_number_].backward_ObsLumDust[output_number_][filter_number_]   = 0.0;
#ifdef ICL
      Gal[galaxy_number_].backward_ObsLumICL[output_number_][filter_number_]    = 0.0;
#endif /* defined ICL */

      Gal[galaxy_number_].forward_ObsLum[output_number_][filter_number_]       = 0.0;
      Gal[galaxy_number_].forward_ObsLumY[output_number_][filter_number_]      = 0.0;
      Gal[galaxy_number_].forward_ObsLumBulge[output_number_][filter_number_]  = 0.0;
      Gal[galaxy_number_].forward_ObsLumBulgeY[output_number_][filter_number_] = 0.0;
      Gal[galaxy_number_].forward_ObsLumDust[output_number_][filter_number_]   = 0.0;
#ifdef ICL
      Gal[galaxy_number_].forward_ObsLumICL[output_number_][filter_number_]    = 0.0;
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
    }
  }
#endif /* not defined POST_PROCESS_MAGS */

#ifdef GALAXYTREE
  Gal[galaxy_number_].FirstProgGal = -1;
#endif /* defined GALAXYTREE */

#ifdef STAR_FORMATION_HISTORY
  sfh_initialise(galaxy_number_);
#endif /* defined STAR_FORMATION_HISTORY */

#ifdef INDIVIDUAL_ELEMENTS
  Gal[galaxy_number_].ColdGas_elements = elements_init();
  Gal[galaxy_number_].DiskMass_elements = elements_init();
  Gal[galaxy_number_].BulgeMass_elements = elements_init();
  Gal[galaxy_number_].HotGas_elements = elements_init();
  Gal[galaxy_number_].EjectedMass_elements = elements_init();
  Gal[galaxy_number_].ICM_elements = elements_init();
#endif /* defined INDIVIDUAL_ELEMENTS */
}


/** @brief Updates properties of central galaxies.
 *
 *   \f$M_{\rm{vir}}\f$, \f$R_{\rm{vir}}\f$ and \f$V_{\rm{vir}}\f$ are only
 *   updated for type 0's. Once galaxies become satellites these quantities
 *   stay unchanged, so will be the values at infall.
 *
 *   If type = 0 then the HotRadius is the Viral Radius, which will be used in
 *   the cooling recipe.
 *
 *   Other infall information will not be used for central galaxies so we do not
 *   care whether they carry the correct values. */
void update_centralgal(const int galaxy_number_, const int halo_number_)
{
  int j_;
  Gal[galaxy_number_].Type = 0;
 
  Gal[galaxy_number_].InfallVmax = Halo[halo_number_].Vmax;
  if(Gal[galaxy_number_].InfallVmaxPeak < Gal[galaxy_number_].InfallVmax)
          Gal[galaxy_number_].InfallVmaxPeak = Gal[galaxy_number_].InfallVmax;
  Gal[galaxy_number_].Rvir = get_virial_radius(halo_number_);
  Gal[galaxy_number_].Vvir = get_virial_velocity(halo_number_);
  Gal[galaxy_number_].Mvir = get_virial_mass(halo_number_);
  Gal[galaxy_number_].InfallSnap = Halo[halo_number_].SnapNum;
  Gal[galaxy_number_].InfallHotGas=Gal[galaxy_number_].HotGas;
  //Gal[galaxy_number_].InfallHotGasRadius=Gal[galaxy_number_].Rvir;
  
  /* if type =0 then hotradius =viral radius, this will be used in the cooling recipe; */
  Gal[galaxy_number_].HotRadius = Gal[galaxy_number_].Rvir;
  Gal[galaxy_number_].MergeOn= 0;
  for (j_=0;j_<3;j_++)
  { Gal[galaxy_number_].HaloSpin[j_] = Halo[halo_number_].Spin[j_]; }
}


/** @brief Updates properties of type 1 galaxies.
 *
 * A dynamical friction decay time scale is calculated
 * for type 1's (as is done for type 2 - introduced for millennium II where the
 * increased resolution means type 1 always retain some dark matter and orbit
 * around for a long time). This is only calculated when the baryonic mass of
 * the type 1 becomes larger than its dark matter mass. The code finds the type
 * 0 to which this galaxy should merge and then sets up the merging clock.
 * */
void update_type_1(const int galaxy_number_, const int halo_number_, const int progenitor_halo_number_)
{
  int current_halo_number_, descendant_halo_number_,first_desc_halo_number_;

  Gal[galaxy_number_].Type = 1;

  if(Gal[galaxy_number_].MergeOn == 0)
  {
    /*If baryonic mass > dark matter mass*/
    if (Gal[galaxy_number_].ColdGas+Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass > Halo[halo_number_].Len*PartMass)
    {
      current_halo_number_ = halo_number_;
      descendant_halo_number_ = Halo[halo_number_].Descendant;
      first_desc_halo_number_ = Halo[Halo[halo_number_].FirstHaloInFOFgroup].Descendant;

      /* In case this is the last snapnum (first_desc_halo_number_ == -1), it means that we tracked all
       * the way down to redshift =0 and mergeon should be trun on. Otherwise, it is the
       * case that the current_halo_number_ halo and the corresponding fof central subhalo are
       * "mysteriously" lost in the dark matter simulation at an intermediate redshift
       * and this galaxy would not be treated further anyway further. Thus the mergeon
       * value is irrelevant. Here mergeon is set to 1. */
      if(descendant_halo_number_ == -1)
      { Gal[galaxy_number_].MergeOn = 1; }

          /* checks when the galaxy "disappears" (when it merges with the type 0) in order to get
           * the type 0 ID into which this type 1 will be merged. */
      while(descendant_halo_number_ >= 0 && first_desc_halo_number_ >= 0)
      {
        if(first_desc_halo_number_ != Halo[first_desc_halo_number_].FirstHaloInFOFgroup)
          break;

        if(Halo[descendant_halo_number_].FirstHaloInFOFgroup != Halo[first_desc_halo_number_].FirstHaloInFOFgroup)
          break;

        if(descendant_halo_number_ != Halo[descendant_halo_number_].FirstHaloInFOFgroup && current_halo_number_ == Halo[descendant_halo_number_].FirstProgenitor)
          if(Gal[galaxy_number_].ColdGas + Gal[galaxy_number_].DiskMass + Gal[galaxy_number_].BulgeMass < Halo[descendant_halo_number_].Len * PartMass)
            break;


        if(descendant_halo_number_ == Halo[descendant_halo_number_].FirstHaloInFOFgroup && current_halo_number_ == Halo[descendant_halo_number_].FirstProgenitor)
          break;

        if(descendant_halo_number_ == Halo[descendant_halo_number_].FirstHaloInFOFgroup && current_halo_number_ != Halo[descendant_halo_number_].FirstProgenitor)
        {
          Gal[galaxy_number_].MergeOn = 1;
          break;
        }

        if(descendant_halo_number_ != Halo[descendant_halo_number_].FirstHaloInFOFgroup && current_halo_number_ != Halo[descendant_halo_number_].FirstProgenitor)
          break;
    
        current_halo_number_=descendant_halo_number_;
        first_desc_halo_number_ = Halo[first_desc_halo_number_].Descendant;
        descendant_halo_number_=Halo[descendant_halo_number_].Descendant;
        
        /* In case this is the last snapnum (first_desc_halo_number_ == -1), it means that we tracked all
          * the way down to redshift =0 and mergeon should be trun on. Otherwise, it is the
          * case that the current_halo_number_ halo and the corresponding fof central subhalo are
          * "mysteriously" lost in the dark matter simulation at an intermediate redshift
          * and this galaxy would not be treated further anyway further. Thus the mergeon
          * value is irrelevant. Here mergeon is set to 1. */
        if(first_desc_halo_number_ == -1)
        {
          if (descendant_halo_number_ == -1)
            Gal[galaxy_number_].MergeOn = 1;
          break;
        }
      }
   
      /*Sets up the dynamical friction decay merging clock as for type 2 galaxies. */
      if (descendant_halo_number_ < 0 || Gal[galaxy_number_].MergeOn == 1)
      {
        Gal[galaxy_number_].MergeOn = 1;
        //In case central galaxy has no progenitor
        if (Halo[Halo[halo_number_].FirstHaloInFOFgroup].FirstProgenitor == -1 )
          Gal[galaxy_number_].MergTime = estimate_merging_time(progenitor_halo_number_,Halo[halo_number_].FirstHaloInFOFgroup,galaxy_number_);
        else
          Gal[galaxy_number_].MergTime = estimate_merging_time(progenitor_halo_number_,Halo[Halo[halo_number_].FirstHaloInFOFgroup].FirstProgenitor,galaxy_number_);
        Gal[galaxy_number_].MergTime -= NumToTime(Halo[halo_number_].SnapNum) - NumToTime(Halo[progenitor_halo_number_].SnapNum);
              //to calculate the position of type 2
        Gal[galaxy_number_].OriMergTime=Gal[galaxy_number_].MergTime;
        Gal[galaxy_number_].OriMvir = get_virial_mass(progenitor_halo_number_);
        Gal[galaxy_number_].OriRvir = get_virial_radius(progenitor_halo_number_);
      }
    }
  }
  /*Mvir, Rvir and Vvir keep their value fixed after infall*/
}         

  
/** @brief Updates properties of type 2 galaxies.
 *
 *  Sets Hot Radius to 0, since all the hot gas has been stripped.
 *  Calls estimate_merging_time to get the merging time scale, calculated for
 *  the orbital decay due to dynamical friction, since this galaxy has lost its
 *  dark matter halo and its position cannot be tracked. */
void update_type_2(const int galaxy_number_, const int halo_number_, const int progenitor_halo_number_, int most_massive_halo_number_)
{
 if(Gal[galaxy_number_].Type != 2)
  {
    int j_;
    for(j_=0; j_<3; j_++)
    {
      Gal[galaxy_number_].Pos_notupdated[j_] = Gal[galaxy_number_].Pos[j_];
      Gal[galaxy_number_].Vel_notupdated[j_] = Gal[galaxy_number_].Vel[j_];
    }
  }

  Gal[galaxy_number_].Type = 2;

  Gal[galaxy_number_].HotRadius = 0.0;

  /* Estimate remaining merging timescale. */
  if (Gal[galaxy_number_].MergeOn == 0)
  {
    //if central galaxy has no progenitor
    if (most_massive_halo_number_ == -1)
      most_massive_halo_number_ = halo_number_;

    Gal[galaxy_number_].MergTime = estimate_merging_time(progenitor_halo_number_,most_massive_halo_number_,galaxy_number_);
    Gal[galaxy_number_].MergTime -= NumToTime(Halo[halo_number_].SnapNum) - NumToTime(Halo[progenitor_halo_number_].SnapNum);
    //to calculate the position of type 2
    Gal[galaxy_number_].OriMergTime=Gal[galaxy_number_].MergTime;
    Gal[galaxy_number_].OriMvir = get_virial_mass(progenitor_halo_number_);
    Gal[galaxy_number_].OriRvir = get_virial_radius(progenitor_halo_number_);
  }
}


/** @brief transfer part of stars between reservoirs and galaxies
 * 
 * Transfers a fraction_ of component cq of galaxy q_ onto component cp of galaxy p_.
 * cp and cq must each be one of: 
 * DiskComponent, BulgeComponent, ICMComponent, BurstComponent
 *
 * If -DTRACK_BURST is set then can also specify BurstComponent as an option.
 * This is different in that it is not a separate component, so it should only
 * be transferred if both cq and cp are BurstComponent
 * 
 * @bug (correct by Stefan Hilbert)
 *       it was tested that "BurstMass" could only be transferred to "BurstMass",
 *       but it was "Burst" to indicate BurstMass transfer and "BurstMass" was
 *       never mentioned elsewhere.
 *       So, "BurstMass" was replaced by "Burst" (and then by BurstComponent)
 *      
 * @bug (correct by Stefan Hilbert)
 *       BurstMass was given to receiving galaxy, but not taken from donor galaxy
 *       code used to be: Gal[q_].BurstMass -= 0.;
 *       now replaced by: Gal[q_].BurstMass -= Mass;
 */
void transfer_stars(const int to_galaxy_number_, const StellarComponentType to_component_, const int from_galaxy_number_, const StellarComponentType from_component_, const double fraction_)
{
  float mass_transferred_;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals metals_transferred_;
#ifdef INDIVIDUAL_ELEMENTS
  struct elements elements_transferred_;
  struct elements sfh_elements_transferred_[SFH_NBIN];
#endif
#else
  float metals_transferred_;
#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef STAR_FORMATION_HISTORY
  int i_;
  float sfh_mass_transferred_[SFH_NBIN];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_metals_transferred_[SFH_NBIN];
#else
  float sfh_metals_transferred_[SFH_NBIN];
#endif
#endif

  /* Sanity checks */
  if (fraction_ > 1.)
  {
    printf("\n*** transfer_stars: fraction_>1 ***\n");
    exit(1);
  }
#ifdef STAR_FORMATION_HISTORY
  if (Gal[to_galaxy_number_].sfh_ibin != Gal[from_galaxy_number_].sfh_ibin)
  {
    printf("\n*** transfer_stars: inconsistent itimes ***\n");
    for(i_=0;i_<SFH_NBIN;i_++)
    {
      printf("Bin[%d] time_1=%e dt_1=%e Nbins_1=%d time_2=%e dt_2=%e Nbins_2=%d\n",i_,
                        Gal[to_galaxy_number_].sfh_t[i_],Gal[to_galaxy_number_].sfh_dt[i_],Gal[to_galaxy_number_].sfh_Nbins[i_],
                        Gal[from_galaxy_number_].sfh_t[i_],Gal[from_galaxy_number_].sfh_dt[i_],Gal[from_galaxy_number_].sfh_Nbins[i_]);
    }
    exit(1);
  }
#endif
#ifdef TRACK_BURST
  if ((to_component_ == BurstComponent && from_component_ != BurstComponent) ||
      (from_component_ == BurstComponent && to_component_ != BurstComponent)   )
  {
    printf("\n*** transfer_stars: used incorrectly with BurstMass ***\n");
    exit(1);
  }
#endif

  //mass_transferred_ and metals to be transfered and taken from galaxy from_galaxy_number_
  if (from_component_ == DiskComponent)
  {
    mass_transferred_ = fraction_*Gal[from_galaxy_number_].DiskMass;
    Gal[from_galaxy_number_].DiskMass -= mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[from_galaxy_number_].sfh_ibin; i_++)
    {
      sfh_mass_transferred_[i_]             = fraction_*Gal[from_galaxy_number_].sfh_DiskMass[i_];
      Gal[from_galaxy_number_].sfh_DiskMass[i_] -= sfh_mass_transferred_[i_];
    }
#endif
    metals_transferred_ = metals_fraction(Gal[from_galaxy_number_].MetalsDiskMass, fraction_);
             metals_deduct_from(&Gal[from_galaxy_number_].MetalsDiskMass, metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_transferred_  = elements_add_fraction(elements_init(),Gal[from_galaxy_number_].DiskMass_elements,fraction_);
             elements_deduct_from(&Gal[from_galaxy_number_].DiskMass_elements, elements_transferred_);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[from_galaxy_number_].sfh_ibin; i_++)
    {
      sfh_metals_transferred_[i_] = metals_fraction(Gal[from_galaxy_number_].sfh_MetalsDiskMass[i_],fraction_);
                      metals_deduct_from(&Gal[from_galaxy_number_].sfh_MetalsDiskMass[i_],sfh_metals_transferred_[i_]);
#ifdef INDIVIDUAL_ELEMENTS
      sfh_elements_transferred_[i_] = elements_add_fraction(elements_init(),Gal[from_galaxy_number_].sfh_ElementsDiskMass[i_],fraction_);
                        elements_deduct_from(&Gal[from_galaxy_number_].sfh_ElementsDiskMass[i_], sfh_elements_transferred_[i_]);
#endif
    }
#endif  
  }
  else if (from_component_ == BulgeComponent)
  {
    mass_transferred_              = fraction_*Gal[from_galaxy_number_].BulgeMass;
    Gal[from_galaxy_number_].BulgeMass -= mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[from_galaxy_number_].sfh_ibin; i_++) 
    { 
      sfh_mass_transferred_[i_]              = fraction_*Gal[from_galaxy_number_].sfh_BulgeMass[i_];
      Gal[from_galaxy_number_].sfh_BulgeMass[i_] -= sfh_mass_transferred_[i_];
    }
#endif
    metals_transferred_ = metals_fraction(Gal[from_galaxy_number_].MetalsBulgeMass,fraction_);
             metals_deduct_from(&Gal[from_galaxy_number_].MetalsBulgeMass,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_transferred_ = elements_add_fraction(elements_init(),Gal[from_galaxy_number_].BulgeMass_elements,fraction_);
            elements_deduct_from(&Gal[from_galaxy_number_].BulgeMass_elements,elements_transferred_);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[from_galaxy_number_].sfh_ibin; i_++)
    {
      sfh_metals_transferred_[i_] = metals_fraction(Gal[from_galaxy_number_].sfh_MetalsBulgeMass[i_],fraction_);
                      metals_deduct_from(&Gal[from_galaxy_number_].sfh_MetalsBulgeMass[i_],sfh_metals_transferred_[i_]);
#ifdef INDIVIDUAL_ELEMENTS
      sfh_elements_transferred_[i_] = elements_add_fraction(elements_init(),Gal[from_galaxy_number_].sfh_ElementsBulgeMass[i_],fraction_);
                        elements_deduct_from(&Gal[from_galaxy_number_].sfh_ElementsBulgeMass[i_],sfh_elements_transferred_[i_]);
#endif
    }
#endif  
  }
  else if (from_component_ == ICMComponent)
  {
    mass_transferred_        = fraction_*Gal[from_galaxy_number_].ICM;
    Gal[from_galaxy_number_].ICM -= mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[from_galaxy_number_].sfh_ibin; i_++)
    {
      sfh_mass_transferred_[i_]=fraction_*Gal[from_galaxy_number_].sfh_ICM[i_];
      Gal[from_galaxy_number_].sfh_ICM[i_] -= sfh_mass_transferred_[i_];
    }
#endif
    metals_transferred_ = metals_fraction(Gal[from_galaxy_number_].MetalsICM, fraction_);
             metals_deduct_from(&Gal[from_galaxy_number_].MetalsICM,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_transferred_ = elements_add_fraction(elements_init(),Gal[from_galaxy_number_].ICM_elements,fraction_);
            elements_deduct_from(&Gal[from_galaxy_number_].ICM_elements,elements_transferred_);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[from_galaxy_number_].sfh_ibin; i_++)
    {
      sfh_metals_transferred_[i_] = metals_fraction(Gal[from_galaxy_number_].sfh_MetalsICM[i_],fraction_);
                      metals_deduct_from(&Gal[from_galaxy_number_].sfh_MetalsICM[i_],sfh_metals_transferred_[i_]);
#ifdef INDIVIDUAL_ELEMENTS
      sfh_elements_transferred_[i_] = elements_add_fraction(elements_init(),Gal[from_galaxy_number_].sfh_ElementsICM[i_],fraction_);
                        elements_deduct_from(&Gal[from_galaxy_number_].sfh_ElementsICM[i_],sfh_elements_transferred_[i_],-1.);
#endif
    }
#endif  
  }
#ifdef TRACK_BURST
  else if (from_component_ == BurstComponent)
  {
    mass_transferred_              = fraction_*Gal[from_galaxy_number_].BurstMass;
    Gal[from_galaxy_number_].BurstMass -= mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[from_galaxy_number_].sfh_ibin; i_++) 
    { 
      sfh_mass_transferred_[i_]=fraction_*Gal[from_galaxy_number_].sfh_BurstMass[i_];
      Gal[from_galaxy_number_].sfh_BurstMass[i_] -= sfh_mass_transferred_[i_];
    }
#endif
  }
#endif
  else
  {
    printf("Unknown component type %s for galaxy from_galaxy_number_ in call to transfer_stars\n", StellarComponentStr[from_component_]);
    exit(1);
  }

  //Add to galaxy to_galaxy_number_
  if (to_component_ == DiskComponent)
  {
    Gal[to_galaxy_number_].DiskMass += mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[to_galaxy_number_].sfh_ibin; i_++) Gal[to_galaxy_number_].sfh_DiskMass[i_] += sfh_mass_transferred_[i_];
#endif
    metals_add_to(&Gal[to_galaxy_number_].MetalsDiskMass,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[to_galaxy_number_].DiskMass_elements, elements_transferred_);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[to_galaxy_number_].sfh_ibin; i_++)
    {
      metals_add_to(&Gal[to_galaxy_number_].sfh_MetalsDiskMass[i_],sfh_metals_transferred_[i_]);
#ifdef INDIVIDUAL_ELEMENTS
      elements_add_to(&Gal[to_galaxy_number_].sfh_ElementsDiskMass[i_],sfh_elements_transferred_[i_]);
#endif
    }
#endif  
  }
  else if (to_component_ == BulgeComponent)
  {
    Gal[to_galaxy_number_].BulgeMass += mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[to_galaxy_number_].sfh_ibin; i_++) 
    { Gal[to_galaxy_number_].sfh_BulgeMass[i_] += sfh_mass_transferred_[i_]; }
#endif
    metals_add_to(&Gal[to_galaxy_number_].MetalsBulgeMass,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[to_galaxy_number_].BulgeMass_elements,elements_transferred_);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[to_galaxy_number_].sfh_ibin; i_++)
    {
      metals_add_to(&Gal[to_galaxy_number_].sfh_MetalsBulgeMass[i_],sfh_metals_transferred_[i_]);
#ifdef INDIVIDUAL_ELEMENTS
      elements_add_to(&Gal[to_galaxy_number_].sfh_ElementsBulgeMass[i_],sfh_elements_transferred_[i_]);
#endif
    }
#endif  
  }
  else if (to_component_ == ICMComponent)
  {
    Gal[to_galaxy_number_].ICM += mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[to_galaxy_number_].sfh_ibin; i_++) { Gal[to_galaxy_number_].sfh_ICM[i_] += sfh_mass_transferred_[i_]; }
#endif
    metals_add_to(&Gal[to_galaxy_number_].MetalsICM,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[to_galaxy_number_].ICM_elements,elements_transferred_);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[to_galaxy_number_].sfh_ibin; i_++)
    {
      metals_add_to(&Gal[to_galaxy_number_].sfh_MetalsICM[i_],sfh_metals_transferred_[i_]);
#ifdef INDIVIDUAL_ELEMENTS
      elements_add_to(&Gal[to_galaxy_number_].sfh_ElementsICM[i_],sfh_elements_transferred_[i_]);
#endif
    }
#endif  
#ifdef TRACK_BURST
  }
  else if (to_component_ == BurstComponent)
  {
    Gal[to_galaxy_number_].BurstMass += mass_transferred_;
#ifdef STAR_FORMATION_HISTORY
    for (i_=0; i_<=Gal[to_galaxy_number_].sfh_ibin; i_++) Gal[to_galaxy_number_].sfh_BurstMass[i_] += sfh_mass_transferred_[i_];
#endif
#endif
  }
  else
  {
    printf("Unknown component type %s for galaxy to_galaxy_number_ in call to transfer_stars\n", StellarComponentStr[to_component_]);
    exit(1);
  }
}


/** @brief transfer part of gas between reservoirs and galaxies
 *
 * Transfers a fraction_ of component from_component_ of galaxy from_galaxy_number_ onto component to_component_ of galaxy to_galaxy_number_.
 * to_component_ and from_component_ must each be one of:
 *   ColdGasComponent
 *   HotGasComponent
 *   EjectedGasComponent
 * 
 * @note used to get info about which function was calling,
 *       but this is (likely by far) the most called function,
 *       so it should be kept as tidy as possible
 */
void transfer_gas(const int to_galaxy_number_, const GasComponentType to_component_, const int from_galaxy_number_, const GasComponentType from_component_, const double fraction_)
{
  float mass_transferred_;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals metals_transferred_;
#ifdef INDIVIDUAL_ELEMENTS
  struct elements elements_transferred_;
#endif
#else
  float metals_transferred_;
#endif

  /* Sanity check */
  if (fraction_ > 1.)
  {
    char error_message_[1000];
    sprintf(error_message_, "transfer_gas: fraction_>1\nfraction = %.11f\nFrom '%s' to '%s\n", fraction_, GasComponentStr[from_component_], GasComponentStr[to_component_]);
    terminate(error_message_);
  }

  /* mass_transferred_ and metals_transferred_ to be transfered taken from galaxy from_galaxy_number_ */
  if (from_component_ == ColdGasComponent)
  {
    mass_transferred_            = fraction_*Gal[from_galaxy_number_].ColdGas;
    Gal[from_galaxy_number_].ColdGas -= mass_transferred_;
    metals_transferred_          = metals_fraction    (Gal[from_galaxy_number_].MetalsColdGas, fraction_);
                      metals_deduct_from(&Gal[from_galaxy_number_].MetalsColdGas, metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_transferred_           = elements_fraction    (Gal[from_galaxy_number_].ColdGas_elements,fraction_);
                      elements_deduct_from(&Gal[from_galaxy_number_].ColdGas_elements,elements_transferred_);
#endif
  }
  else if (from_component_ == HotGasComponent)
  {
    mass_transferred_           = fraction_*Gal[from_galaxy_number_].HotGas;
    Gal[from_galaxy_number_].HotGas -= mass_transferred_;
    metals_transferred_         = metals_fraction    (Gal[from_galaxy_number_].MetalsHotGas,fraction_);
                     metals_deduct_from(&Gal[from_galaxy_number_].MetalsHotGas,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_transferred_          = elements_fraction    (Gal[from_galaxy_number_].HotGas_elements,fraction_);
                     elements_deduct_from(&Gal[from_galaxy_number_].HotGas_elements,elements_transferred_);
#endif
#ifdef METALS_SELF
   metals_deduct_from(&Gal[from_galaxy_number_].MetalsHotGasSelf, metals_transferred_);
#endif
  }
  else if (from_component_ == EjectedGasComponent)
  {
    mass_transferred_                = fraction_*Gal[from_galaxy_number_].EjectedMass;
    Gal[from_galaxy_number_].EjectedMass -= mass_transferred_;
    metals_transferred_              = metals_fraction    (Gal[from_galaxy_number_].MetalsEjectedMass, fraction_);
                          metals_deduct_from(&Gal[from_galaxy_number_].MetalsEjectedMass, metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_transferred_               = elements_fraction  (Gal[from_galaxy_number_].EjectedMass_elements, fraction_);
                          metals_deduct_from(&Gal[from_galaxy_number_].EjectedMass_elements, elements_transferred_);
#endif
  }
  else
  {
    printf("Unknown component type %s in call to transfer_gas\n", GasComponentStr[from_component_]);
    exit(1);
  }
 
 //Add to galaxy to_galaxy_number_
  if (to_component_ == ColdGasComponent)
  {
    Gal[to_galaxy_number_].ColdGas += mass_transferred_;
    metals_add_to(&Gal[to_galaxy_number_].MetalsColdGas,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[to_galaxy_number_].ColdGas_elements,elements_transferred_);
#endif
  }
  else if (to_component_ == HotGasComponent)
  {
    Gal[to_galaxy_number_].HotGas += mass_transferred_;
    metals_add_to(&Gal[to_galaxy_number_].MetalsHotGas,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[to_galaxy_number_].HotGas_elements,elements_transferred_);
#endif
#ifdef METALS_SELF
    if (to_galaxy_number_==from_galaxy_number_) metals_add_to(&Gal[to_galaxy_number_].MetalsHotGasSelf,metals_transferred_);
#endif
  }
  else if (to_component_ == EjectedGasComponent)
  {
    Gal[to_galaxy_number_].EjectedMass += mass_transferred_;
    metals_add_to(&Gal[to_galaxy_number_].MetalsEjectedMass,metals_transferred_);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[to_galaxy_number_].EjectedMass_elements,elements_transferred_);
#endif
  }
  else
  {
    printf("Unknown component type %s in call to transfer_gas\n", GasComponentStr[to_component_]);
    exit(1);
  }
  return;
}


/** @brief sanity checks on the masses of different components
 *
 *
 * Some sanity checks on the masses of different components. 
 * If due to rounding error, then apply a correction;
 * otherwise print error message and exit
 *
 * @note ROB: Should probably make some of these for the elements
 */
void perform_mass_checks(char tag_string_[], const int galaxy_number_)
{
#ifdef STAR_FORMATION_HISTORY
  int i_;
  double sfh_sum;
#endif

  if(Gal[galaxy_number_].ColdGas < 1.e-8)
  {
    Gal[galaxy_number_].ColdGas = 0.;
    Gal[galaxy_number_].MetalsColdGas = metals_init();
  }

  //check if the gas mass is less than 0
  if(Gal[galaxy_number_].ColdGas < 0.0)
  {
    if (Gal[galaxy_number_].ColdGas > -1e-7)
      Gal[galaxy_number_].ColdGas = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, ColdGas < 0. ***\n",tag_string_);
      printf("                ColdGas[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].ColdGas);
      terminate("");
    }
  }

  //check if the mass in metals is less than 0
  if(metals_total(Gal[galaxy_number_].MetalsColdGas) < 0.0)
  {
    if (metals_total(Gal[galaxy_number_].MetalsColdGas) > -1e-7)
      Gal[galaxy_number_].MetalsColdGas = metals_init();
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsColdGas < 0. ***\n",tag_string_);
      printf("                MetalsColdGas[%d] = %g\n",galaxy_number_,metals_total(Gal[galaxy_number_].MetalsColdGas));
      terminate("");
    }
  }

  //check if the mass in metals is greater than the gas mass
  if(metals_total(Gal[galaxy_number_].MetalsColdGas) > Gal[galaxy_number_].ColdGas)
  {
    if (metals_total(Gal[galaxy_number_].MetalsColdGas) < 1e-7)
      Gal[galaxy_number_].MetalsColdGas = metals_fraction(Gal[galaxy_number_].MetalsColdGas,
                                           Gal[galaxy_number_].ColdGas/metals_total(Gal[galaxy_number_].MetalsColdGas));
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsColdGas > ColdGas ***\n",tag_string_);
      printf("          MetalsColdGas[%d] = %g\n",galaxy_number_,metals_total(Gal[galaxy_number_].MetalsColdGas));
      printf("                ColdGas[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].ColdGas);
      terminate("");
    }
  }

  if(Gal[galaxy_number_].HotGas < 0.0)
  {
    if (Gal[galaxy_number_].HotGas > -1e-7)
      Gal[galaxy_number_].HotGas = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, HotGas < 0. ***\n",tag_string_);
      printf("                HotGas[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].HotGas);
      terminate("");
    }
  }

  if(metals_total(Gal[galaxy_number_].MetalsHotGas) < 0.0)
  {
    if (metals_total(Gal[galaxy_number_].MetalsHotGas) > -1e-7)
      Gal[galaxy_number_].MetalsHotGas = metals_init();
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsHotGas < 0. ***\n",tag_string_);
      printf("                MetalsHotGas[%d] = %g\n",galaxy_number_,metals_total(Gal[galaxy_number_].MetalsHotGas));
      terminate("");
    }
  }

  if(metals_total(Gal[galaxy_number_].MetalsHotGas) > Gal[galaxy_number_].HotGas)
  {
    if (metals_total(Gal[galaxy_number_].MetalsHotGas) < 1e-7)
      Gal[galaxy_number_].MetalsHotGas = metals_fraction(Gal[galaxy_number_].MetalsHotGas,
                                           Gal[galaxy_number_].HotGas/metals_total(Gal[galaxy_number_].MetalsHotGas));
   else
   {
      printf("\n***  Mass check error, called from: %s, MetalsHotGas > HotGas ***\n",tag_string_);
      printf("          MetalsHotGas[%d] = %g\n",galaxy_number_,metals_total(Gal[galaxy_number_].MetalsHotGas));
      printf("                HotGas[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].HotGas);
      printf("          MetalsHotGas[%d] = %.11f\n",galaxy_number_,metals_total(Gal[galaxy_number_].MetalsHotGas));
      printf("                HotGas[%d] = %.11f\n",galaxy_number_,Gal[galaxy_number_].HotGas);
      printf("             BulgeMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].BulgeMass);
      printf("           EjectedMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].EjectedMass);
      printf("                  Snapnum = %i_\n",Gal[galaxy_number_].SnapNum);
      terminate("");
    }
  }

  if(Gal[galaxy_number_].EjectedMass < 0.0)
  {
    if (Gal[galaxy_number_].EjectedMass > -1e-7)
      Gal[galaxy_number_].EjectedMass = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, EjectedMass < 0. ***\n",tag_string_);
      printf("                EjectedMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].EjectedMass);
      terminate("");
    }
  }

  if(metals_total(Gal[galaxy_number_].MetalsEjectedMass) < 0.0)
  {
    if (metals_total(Gal[galaxy_number_].MetalsEjectedMass) > -1e-7)
      Gal[galaxy_number_].MetalsEjectedMass = metals_init();
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsEjectedMass < 0. ***\n",tag_string_);
      printf("                MetalsEjectedMass[%d] = %g\n",galaxy_number_,metals_total(Gal[galaxy_number_].MetalsEjectedMass));
      terminate("");
    }
  }

  if(metals_total(Gal[galaxy_number_].MetalsEjectedMass) > Gal[galaxy_number_].EjectedMass)
  {
    if (metals_total(Gal[galaxy_number_].MetalsEjectedMass) < 1e-7)
      Gal[galaxy_number_].MetalsEjectedMass = metals_fraction(Gal[galaxy_number_].MetalsEjectedMass,
                                           Gal[galaxy_number_].EjectedMass/metals_total(Gal[galaxy_number_].MetalsEjectedMass));
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsEjectedMass > EjectedMass ***\n",tag_string_);
      printf("          MetalsEjectedMass[%d] = %g\n",galaxy_number_,metals_total(Gal[galaxy_number_].MetalsEjectedMass));
      printf("                EjectedMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].EjectedMass);
      terminate("");
    }
  }

  if(Gal[galaxy_number_].DiskMass < 0.0)
  {
    if (Gal[galaxy_number_].DiskMass > -1e-7)
      Gal[galaxy_number_].DiskMass = 0.;
    else 
    {
      printf("\n*** Mass check error, called from: %s, DiskMass < 0. ***\n",tag_string_);
      printf("                DiskMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].DiskMass);
      terminate("");
    }
  }

  if(Gal[galaxy_number_].BulgeMass < 0.0)
  {
    if (Gal[galaxy_number_].BulgeMass > -1e-7)
      Gal[galaxy_number_].BulgeMass = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, BulgeMass < 0. ***\n",tag_string_);
      printf("                BulgeMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].BulgeMass);
      terminate("");
    }
  }

  if(Gal[galaxy_number_].ICM < 0.0)
  {
    if (Gal[galaxy_number_].ICM > -1e-7)
      Gal[galaxy_number_].ICM = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, ICM < 0. ***\n",tag_string_);
      printf("                ICM[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].ICM);
      terminate("");
    }
  }

#ifdef TRACK_BURST
  if(Gal[galaxy_number_].BurstMass < 0.0) 
  {
    if (Gal[galaxy_number_].BurstMass > -1e-7)
      Gal[galaxy_number_].BurstMass = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, BurstMass < 0. ***\n",tag_string_);
      printf("                BurstMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].BurstMass);
      terminate("");
    }
  }
#endif

  /* If DETAILED_METALS_AND_MASS_RETURN, sfh stores accumulation of 'stars', not 'stars-recycFrac'.
   * Therefore, it's sum doesn't equal DiskMass any more.*/
#ifndef DETAILED_METALS_AND_MASS_RETURN
#ifdef STAR_FORMATION_HISTORY
  sfh_sum=-Gal[galaxy_number_].DiskMass;
  for (i_=0; i_<=Gal[galaxy_number_].sfh_ibin; i_++) sfh_sum+=Gal[galaxy_number_].sfh_DiskMass[i_];
  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[galaxy_number_].DiskMass) ||
     (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[galaxy_number_].DiskMass))
  {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                DiskMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].DiskMass);
    printf("            sfh_DiskMass[%d] = %g\n",galaxy_number_,sfh_sum+Gal[galaxy_number_].DiskMass);
    char error_message_[1000];
    sprintf(error_message_, "\n*** Mass check error, called from: %s, Inconsistent sfh for DiskMass.*** \n",tag_string_);
    terminate(error_message_);
  }

  sfh_sum=-Gal[galaxy_number_].BulgeMass;
  for (i_=0; i_<=Gal[galaxy_number_].sfh_ibin; i_++) sfh_sum+=Gal[galaxy_number_].sfh_BulgeMass[i_];
  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[galaxy_number_].BulgeMass) ||
     (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[galaxy_number_].BulgeMass))
  {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                BulgeMass[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].BulgeMass);
    printf("            sfh_BulgeMass[%d] = %g\n",galaxy_number_,sfh_sum+Gal[galaxy_number_].BulgeMass);
    char error_message_[1000];
    sprintf(error_message_, "\n*** Mass check error, called from: %s, Inconsistent sfh for BulgeMass. ***\n",tag_string_);
    terminate(error_message_);
  }

  sfh_sum=-Gal[galaxy_number_].ICM;
  for (i_=0; i_<=Gal[galaxy_number_].sfh_ibin; i_++) sfh_sum+=Gal[galaxy_number_].sfh_ICM[i_];
  if(sfh_sum < -1e-4 || sfh_sum > 1e-4)
  {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                ICM[%d] = %g\n",galaxy_number_,Gal[galaxy_number_].ICM);
    printf("            sfh_ICM[%d] = %g\n",galaxy_number_,sfh_sum+Gal[galaxy_number_].ICM);
    for (i_=0; i_<=Gal[galaxy_number_].sfh_ibin; i_++)
      printf("%d %f\n",i_,Gal[galaxy_number_].sfh_ICM[i_]);
    char error_message_[1000];
    sprintf(error_message_, "\n*** Mass check error, called from: %s, Inconsistent sfh for ICM. ***\n",tag_string_);
    terminate(error_message_);
  }
#endif //STAR_FORMATION_HISTORY
#endif //DETAILED_ENRICHEMENT

  return;
}


/** @brief print galaxy properties */
void print_galaxy(char tag_string_[], const int galaxy_number_, const int halo_number_)
{
/*        printf("%s Hnr=%d firstinFOF=%d prog=%d nestprog=%d Descendant=%d gal=%d Type=%d\n",
                        tag_string_, Gal[galaxy_number_].HaloNr, Halo[halo_number_].FirstHaloInFOFgroup, Halo[halo_number_].FirstProgenitor,
                        Halo[halo_number_].NextProgenitor, Halo[halo_number_].Descendant, galaxy_number_, Gal[galaxy_number_].Type);
        printf("     Mvir=%0.3e Vvir=%0.3e Hot=%0.3e Cold=%0.3e Eject=%0.3e disk=%0.3e bulge=%0.3e  GasDiskRadius=%0.3e\n",
                        Gal[galaxy_number_].Mvir*1.e10, Gal[galaxy_number_].Vvir, Gal[galaxy_number_].HotGas*1.e10, Gal[galaxy_number_].ColdGas*1.e10, Gal[galaxy_number_].EjectedMass*1.e10,
                        Gal[galaxy_number_].DiskMass*1.e10, Gal[galaxy_number_].BulgeMass*1.e10, Gal[galaxy_number_].GasDiskRadius);

        printf("     HotMetals=%0.3e ColdMetals=%0.3e diskMetals=%0.3e bulgeMetals=%0.3e\n",
                        Gal[galaxy_number_].MetalsHotGas, Gal[galaxy_number_].MetalsColdGas,Gal[galaxy_number_].MetalsDiskMass, Gal[galaxy_number_].MetalsBulgeMass);


        printf("     x=%0.3f y=%0.3f z=%0.3f vx=%0.3f vy=%0.3f vz=%0.3f\n",
                        Gal[galaxy_number_].Pos[0],Gal[galaxy_number_].Pos[1],Gal[galaxy_number_].Pos[2],Gal[galaxy_number_].Vel[0],Gal[galaxy_number_].Vel[1],Gal[galaxy_number_].Vel[2]);*/

  printf("%s Hnr=%d gal=%d Gal[gal].Hnr=%d Gal[gal].snap=%d\n",tag_string_, halo_number_, galaxy_number_, Gal[galaxy_number_].HaloNr, Gal[galaxy_number_].SnapNum);
  printf(" Hot=%0.3e Cold=%0.3e Eject=%0.3e disk=%0.3e bulge=%0.3e  GasDiskRadius=%0.3e StellarDiskRadius=%0.3e BulgeSize=%0.3e\n",
         Gal[galaxy_number_].HotGas*1.e10, Gal[galaxy_number_].ColdGas*1.e10, Gal[galaxy_number_].EjectedMass*1.e10,
         Gal[galaxy_number_].DiskMass*1.e10, Gal[galaxy_number_].BulgeMass*1.e10, Gal[galaxy_number_].GasDiskRadius, Gal[galaxy_number_].StellarDiskRadius, Gal[galaxy_number_].BulgeSize);

  if(isnan(Gal[galaxy_number_].GasDiskRadius))
  { exit(1); }
}
