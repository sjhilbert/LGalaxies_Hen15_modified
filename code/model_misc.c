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


/** @brief Calculates the separation of galaxies p_ and q, allowing for wrapping */
double separation_gal(const int p_, const int q)
 {
  int i;
  double sep1,sep2;

  sep2=0.;
  for (i=0; i<3; i++)
    {
      sep1 =  Gal[p_].Pos[i] - Gal[q].Pos[i];
      sep1 =  wrap(sep1, BoxSize);
      sep2+= sep1*sep1;
    }
  return sqrt(sep2);
}


/** @brief Calculates the separation of galaxies p_ and q, allowing for wrapping */
double separation_halo(const int p_, const  int q) 
{
  int i;
  double sep1, sep2;

  sep2=0.;
  for (i=0; i<3; i++)
    {
          sep1 = Halo[p_].Pos[i] - Halo[q].Pos[i];
          sep1 = wrap(sep1,BoxSize);
          sep2+=sep1*sep1;
    }
  return sqrt(sep2);
}


/** @brief computes disk radius
 * 
 *        computes disk radius from halo spin
 * 
 * @note This routine is no longer called (after Guo2010) */
double get_disk_radius(const int halo_number_, const int p_)
{
  if(DiskRadiusModel == 1)
  {
    /*  See Mo, Mao & White (1998) eq12, and using a Bullock style lambda.  Since this is the scale length
        we take the typical star forming region as 3 times this using the Milky Way as an approximate guide */
    return 1.5 * sqrt(Halo[halo_number_].Spin[0] * Halo[halo_number_].Spin[0] + 
                      Halo[halo_number_].Spin[1] * Halo[halo_number_].Spin[1] +
                      Halo[halo_number_].Spin[2] * Halo[halo_number_].Spin[2]) / Gal[p_].Vvir;
  }
  else
    /*  simple prescription */
    return 0.1 * Gal[p_].Rvir;
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
double get_initial_disk_radius(const int halonr_, const int p_)
{
  if(DiskRadiusModel == 0 || DiskRadiusModel == 1)
  {
    const double spin_squared_ = Halo[halonr_].Spin[0] * Halo[halonr_].Spin[0] +
                                 Halo[halonr_].Spin[1] * Halo[halonr_].Spin[1] +
                                 Halo[halonr_].Spin[2] * Halo[halonr_].Spin[2];  

    if(Gal[p_].GasSpin[0]==0 && Gal[p_].GasSpin[1]==0 && Gal[p_].GasSpin[2]==0)
      return 0.1 * Gal[p_].Rvir;
    else
      return 1.5 * sqrt(spin_squared_) / Gal[p_].Vvir;
  }
  else
    /*  simpler prescription */
    return 0.1 * Gal[p_].Rvir;
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
void set_gas_disk_radius(const int p_)
{
  const double spin_squared_ = Gal[p_].GasSpin[0] * Gal[p_].GasSpin[0] +
                               Gal[p_].GasSpin[1] * Gal[p_].GasSpin[1] +
                               Gal[p_].GasSpin[2] * Gal[p_].GasSpin[2];
  if(spin_squared_ > 0.)
  {
    Gal[p_].GasDiskRadius = (Gal[p_].Type == 0) ? 1.5 * sqrt(spin_squared_) / Gal[p_].Vmax  :
                                                  1.5 * sqrt(spin_squared_) / Gal[p_].InfallVmax;
  }
  else
  { Gal[p_].GasDiskRadius = 0.1 * Gal[p_].Rvir;  }
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
void set_stellar_disk_radius(const int p_)
{
  const double spin_squared_ = Gal[p_].StellarSpin[0] * Gal[p_].StellarSpin[0] +
                               Gal[p_].StellarSpin[1] * Gal[p_].StellarSpin[1] +
                               Gal[p_].StellarSpin[2] * Gal[p_].StellarSpin[2];
  if(spin_squared_ > 0.)
  {
    Gal[p_].StellarDiskRadius = (Gal[p_].Type == 0) ? 1.5 * sqrt(spin_squared_) / Gal[p_].Vmax : 
                                                      1.5 * sqrt(spin_squared_) / Gal[p_].InfallVmax;
  }
  else
  { Gal[p_].StellarDiskRadius = 0.1 * Gal[p_].Rvir;  }
}


/** @brief Initializes the Galaxy Structure by setting all its
 *         elements to zero. */
void init_galaxy(const int p_, const int halo_number_)
{
  int j_, output_number_;
  
#ifndef POST_PROCESS_MAGS
  int filter_number_;
#endif /* not defined POST_PROCESS_MAGS */

  /* make explicitly sure that the whole galaxy structure has defined 0 values */
  memset(&Gal[p_], 0, sizeof(struct GALAXY));

  Gal[p_].NextGalaxy = -1;
#ifdef GALAXYTREE
  Gal[p_].FirstProgGal = -1;
#endif /* defined GALAXYTREE */

  if(halo_number_ != Halo[halo_number_].FirstHaloInFOFgroup)
  {
    terminate("Hah?\n");
  }

  Gal[p_].Type = 0;

  Gal[p_].HaloNr = halo_number_;
  Gal[p_].MostBoundID = Halo[halo_number_].MostBoundID;
  Gal[p_].SnapNum = Halo[halo_number_].SnapNum - 1;
#ifdef HALOPROPERTIES
  Gal[p_].HaloM_Mean200 = Halo[halo_number_].M_Mean200;
  Gal[p_].HaloM_Crit200 = Halo[halo_number_].M_Crit200;
  Gal[p_].HaloM_TopHat = Halo[halo_number_].M_TopHat;
  Gal[p_].HaloVelDisp = Halo[halo_number_].VelDisp;
  Gal[p_].HaloVmax = Halo[halo_number_].Vmax;
#endif /* defined HALOPROPERTIES */

  for(j_ = 0; j_ < 3; j_++)
  {
    Gal[p_].Pos[j_] = Halo[halo_number_].Pos[j_];
    Gal[p_].Vel[j_] = Halo[halo_number_].Vel[j_];
    Gal[p_].GasSpin[j_] = Halo[halo_number_].Spin[j_];
    Gal[p_].StellarSpin[j_] = Halo[halo_number_].Spin[j_];
    Gal[p_].HaloSpin[j_] = Halo[halo_number_].Spin[j_];
    Gal[p_].MergCentralPos[j_] = Gal[p_].Pos[j_];
    Gal[p_].DistanceToCentralGal[j_]=0.0;
#ifdef HALOPROPERTIES
    Gal[p_].HaloPos[j_] = Halo[halo_number_].Pos[j_];
    Gal[p_].HaloVel[j_] = Halo[halo_number_].Vel[j_];
#endif /* defined HALOPROPERTIES */
  }

  Gal[p_].Len = Halo[halo_number_].Len;
  Gal[p_].Vmax = Halo[halo_number_].Vmax;
  Gal[p_].InfallVmax = Halo[halo_number_].Vmax;
  Gal[p_].InfallVmaxPeak = Gal[p_].InfallVmax;
  Gal[p_].Vvir = get_virial_velocity(halo_number_);
  Gal[p_].Mvir = get_virial_mass(halo_number_);
  Gal[p_].Rvir = get_virial_radius(halo_number_);
  Gal[p_].MergeSat = 0.0;
  Gal[p_].InfallSnap = Halo[halo_number_].SnapNum;

  Gal[p_].ColdGas = 0.0;
  Gal[p_].DiskMass = 0.0;
  Gal[p_].BulgeMass = 0.0;
  Gal[p_].HotGas = 0.0;
  Gal[p_].EjectedMass = 0.0;
  Gal[p_].ICM = 0.0;
#ifdef TRACK_BURST
  Gal[p_].BurstMass = 0.0;
#endif /* defined TRACK_BURST */

  Gal[p_].BlackHoleMass = 0.0;
  /*ram pressure*/
  Gal[p_].HotRadius=Gal[p_].Rvir;
#ifdef GALAXYTREE
  Gal[p_].DisruptOn = 0;
#endif /* defined GALAXYTREE */

#ifdef DETAILED_METALS_AND_MASS_RETURN
  Gal[p_].MetalsColdGas = metals_init();
  Gal[p_].MetalsDiskMass = metals_init();
  Gal[p_].MetalsBulgeMass = metals_init();
  Gal[p_].MetalsHotGas = metals_init();
  Gal[p_].MetalsEjectedMass = metals_init();
  Gal[p_].MetalsICM = metals_init();
#ifdef METALS_SELF
  Gal[p_].MetalsHotGasSelf = metals_init();
#endif /* defined METALS_SELF */
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
  Gal[p_].MetalsColdGas = 0.0;
  Gal[p_].MetalsDiskMass = 0.0;
  Gal[p_].MetalsBulgeMass = 0.0;
  Gal[p_].MetalsHotGas = 0.0;
  Gal[p_].MetalsEjectedMass = 0.0;
  Gal[p_].MetalsICM = 0.0;
#ifdef METALS_SELF
  Gal[p_].MetalsHotGasSelf = 0.0;
#endif /* defined METALS_SELF */
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */

  //inclination defined as the angle between galaxy spin and the z-axis
  Gal[p_].CosInclination = 0.0;

  Gal[p_].PrimordialAccretionRate = 0.0;
  Gal[p_].CoolingRate = 0.0;
  Gal[p_].CoolingRate_beforeAGN = 0.0;
  Gal[p_].CoolingRadius = 0.0;
  Gal[p_].CoolingGas = 0.0;
  Gal[p_].QuasarAccretionRate=0.0;
  Gal[p_].RadioAccretionRate=0.0;
  Gal[p_].AGNheatingFromCentral = 0.0;

  Gal[p_].Sfr = 0.0;
  Gal[p_].SfrBulge = 0.0;

  Gal[p_].StarMerge=0.0;

  Gal[p_].XrayLum = 0.0;
  Gal[p_].GasDiskRadius = get_initial_disk_radius(halo_number_, p_);
  Gal[p_].StellarDiskRadius = Gal[p_].GasDiskRadius;
  Gal[p_].BulgeSize = 0.0;

  Gal[p_].OriMergTime = 0.0;
  Gal[p_].MergTime = 0.0;

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    Gal[p_].MassWeightAge[output_number_] = 0.0;

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    Gal[p_].rbandWeightAge[output_number_] = 0.0;

#ifndef POST_PROCESS_MAGS
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
#ifdef OUTPUT_REST_MAGS
      Gal[p_].Lum[output_number_][filter_number_]         = 0.0;
      Gal[p_].YLum[output_number_][filter_number_]        = 0.0;
      Gal[p_].LumBulge[output_number_][filter_number_]    = 0.0;
      Gal[p_].YLumBulge[output_number_][filter_number_]   = 0.0;
      Gal[p_].LumDust[output_number_][filter_number_]     = 0.0;
#ifdef ICL
      Gal[p_].ICLLum[output_number_][filter_number_]      = 0.0;
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
      Gal[p_].ObsLum[output_number_][filter_number_]        = 0.0;
      Gal[p_].ObsYLum[output_number_][filter_number_]       = 0.0;
      Gal[p_].ObsLumBulge[output_number_][filter_number_]   = 0.0;
      Gal[p_].ObsYLumBulge[output_number_][filter_number_]  = 0.0;
      Gal[p_].ObsLumDust[output_number_][filter_number_]    = 0.0;
#ifdef ICL
      Gal[p_].ObsICL[output_number_][filter_number_]        = 0.0;
#endif /* defined ICL */
          
#ifdef OUTPUT_FB_OBS_MAGS
      Gal[p_].backward_ObsLum[output_number_][filter_number_]       = 0.0;
      Gal[p_].backward_ObsYLum[output_number_][filter_number_]      = 0.0;
      Gal[p_].backward_ObsLumBulge[output_number_][filter_number_]  = 0.0;
      Gal[p_].backward_ObsYLumBulge[output_number_][filter_number_] = 0.0;
      Gal[p_].backward_ObsLumDust[output_number_][filter_number_]   = 0.0;
#ifdef ICL
      Gal[p_].backward_ObsICL[output_number_][filter_number_]       = 0.0;
#endif /* defined ICL */

      Gal[p_].forward_ObsLum[output_number_][filter_number_]       = 0.0;
      Gal[p_].forward_ObsYLum[output_number_][filter_number_]      = 0.0;
      Gal[p_].forward_ObsLumBulge[output_number_][filter_number_]  = 0.0;
      Gal[p_].forward_ObsYLumBulge[output_number_][filter_number_] = 0.0;
      Gal[p_].forward_ObsLumDust[output_number_][filter_number_]   = 0.0;
#ifdef ICL
      Gal[p_].forward_ObsICL[output_number_][filter_number_]       = 0.0;
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
    }
  }
#endif /* not defined POST_PROCESS_MAGS */

#ifdef GALAXYTREE
  Gal[p_].FirstProgGal = -1;
#endif /* defined GALAXYTREE */

#ifdef STAR_FORMATION_HISTORY
  sfh_initialise(p_);
#endif /* defined STAR_FORMATION_HISTORY */

#ifdef INDIVIDUAL_ELEMENTS
  Gal[p_].ColdGas_elements = elements_init();
  Gal[p_].DiskMass_elements = elements_init();
  Gal[p_].BulgeMass_elements = elements_init();
  Gal[p_].HotGas_elements = elements_init();
  Gal[p_].EjectedMass_elements = elements_init();
  Gal[p_].ICM_elements = elements_init();
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
void update_centralgal(const int ngal, const int halo_number_)
{
  int j_;
  Gal[ngal].Type = 0;
 
  Gal[ngal].InfallVmax = Halo[halo_number_].Vmax;
  if(Gal[ngal].InfallVmaxPeak < Gal[ngal].InfallVmax)
          Gal[ngal].InfallVmaxPeak = Gal[ngal].InfallVmax;
  Gal[ngal].Rvir = get_virial_radius(halo_number_);
  Gal[ngal].Vvir = get_virial_velocity(halo_number_);
  Gal[ngal].Mvir = get_virial_mass(halo_number_);
  Gal[ngal].InfallSnap = Halo[halo_number_].SnapNum;
  Gal[ngal].InfallHotGas=Gal[ngal].HotGas;
  //Gal[ngal].InfallHotGasRadius=Gal[ngal].Rvir;
  
  /* if type =0 then hotradius =viral radius, this will be used in the cooling recipe; */
  Gal[ngal].HotRadius = Gal[ngal].Rvir;
  Gal[ngal].MergeOn= 0;
  for (j_=0;j_<3;j_++)
  { Gal[ngal].HaloSpin[j_] = Halo[halo_number_].Spin[j_]; }
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
void update_type_1(const int ngal, const int halo_number_, const int prog)
{
  int current,descendant,firstdes;

  Gal[ngal].Type = 1;

  if(Gal[ngal].MergeOn == 0)
  {
    /*If baryonic mass > dark matter mass*/
    if (Gal[ngal].ColdGas+Gal[ngal].DiskMass+Gal[ngal].BulgeMass > Halo[halo_number_].Len*PartMass)
    {
      current = halo_number_;
      descendant = Halo[halo_number_].Descendant;
      firstdes = Halo[Halo[halo_number_].FirstHaloInFOFgroup].Descendant;

      /* In case this is the last snapnum (firstdes == -1), it means that we tracked all
       * the way down to redshift =0 and mergeon should be trun on. Otherwise, it is the
       * case that the current halo and the corresponding fof central subhalo are
       * "mysteriously" lost in the dark matter simulation at an intermediate redshift
       * and this galaxy would not be treated further anyway further. Thus the mergeon
       * value is irrelevant. Here mergeon is set to 1. */
      if(descendant == -1)
      { Gal[ngal].MergeOn = 1; }

          /* checks when the galaxy "disappears" (when it merges with the type 0) in order to get
           * the type 0 ID into which this type 1 will be merged. */
      while(descendant >= 0 && firstdes >= 0)
      {
        if(firstdes != Halo[firstdes].FirstHaloInFOFgroup)
          break;

        if(Halo[descendant].FirstHaloInFOFgroup != Halo[firstdes].FirstHaloInFOFgroup)
          break;

        if(descendant != Halo[descendant].FirstHaloInFOFgroup && current == Halo[descendant].FirstProgenitor)
          if(Gal[ngal].ColdGas + Gal[ngal].DiskMass + Gal[ngal].BulgeMass < Halo[descendant].Len * PartMass)
            break;


        if(descendant == Halo[descendant].FirstHaloInFOFgroup && current == Halo[descendant].FirstProgenitor)
          break;

        if(descendant == Halo[descendant].FirstHaloInFOFgroup && current != Halo[descendant].FirstProgenitor)
        {
          Gal[ngal].MergeOn = 1;
          break;
        }

        if(descendant != Halo[descendant].FirstHaloInFOFgroup && current != Halo[descendant].FirstProgenitor)
          break;
    
        current=descendant;
        firstdes = Halo[firstdes].Descendant;
        descendant=Halo[descendant].Descendant;
        
        /* In case this is the last snapnum (firstdes == -1), it means that we tracked all
          * the way down to redshift =0 and mergeon should be trun on. Otherwise, it is the
          * case that the current halo and the corresponding fof central subhalo are
          * "mysteriously" lost in the dark matter simulation at an intermediate redshift
          * and this galaxy would not be treated further anyway further. Thus the mergeon
          * value is irrelevant. Here mergeon is set to 1. */
        if(firstdes == -1)
        {
          if (descendant == -1)
            Gal[ngal].MergeOn = 1;
          break;
        }
      }
   
      /*Sets up the dynamical friction decay merging clock as for type 2 galaxies. */
      if (descendant < 0 || Gal[ngal].MergeOn == 1)
      {
        Gal[ngal].MergeOn = 1;
        //In case central galaxy has no progenitor
        if (Halo[Halo[halo_number_].FirstHaloInFOFgroup].FirstProgenitor == -1 )
          Gal[ngal].MergTime = estimate_merging_time(prog,Halo[halo_number_].FirstHaloInFOFgroup,ngal);
        else
          Gal[ngal].MergTime = estimate_merging_time(prog,Halo[Halo[halo_number_].FirstHaloInFOFgroup].FirstProgenitor,ngal);
        Gal[ngal].MergTime -= NumToTime(Halo[halo_number_].SnapNum) - NumToTime(Halo[prog].SnapNum);
              //to calculate the position of type 2
        Gal[ngal].OriMergTime=Gal[ngal].MergTime;
        Gal[ngal].OriMvir = get_virial_mass(prog);
        Gal[ngal].OriRvir = get_virial_radius(prog);
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
void update_type_2(const int ngal, const int halo_number_, const int prog, int mostmassive)
{
 if(Gal[ngal].Type != 2)
  {
    int j_;
    for(j_=0; j_<3; j_++)
    {
      Gal[ngal].Pos_notupdated[j_] = Gal[ngal].Pos[j_];
      Gal[ngal].Vel_notupdated[j_] = Gal[ngal].Vel[j_];
    }
  }

  Gal[ngal].Type = 2;

  Gal[ngal].HotRadius = 0.0;

  /* Estimate remaining merging timescale. */
  if (Gal[ngal].MergeOn == 0)
  {
    //if central galaxy has no progenitor
    if (mostmassive == -1)
      mostmassive = halo_number_;

    Gal[ngal].MergTime = estimate_merging_time(prog,mostmassive,ngal);
    Gal[ngal].MergTime -= NumToTime(Halo[halo_number_].SnapNum) - NumToTime(Halo[prog].SnapNum);
    //to calculate the position of type 2
    Gal[ngal].OriMergTime=Gal[ngal].MergTime;
    Gal[ngal].OriMvir = get_virial_mass(prog);
    Gal[ngal].OriRvir = get_virial_radius(prog);
  }
}


/** @brief transfer part of stars between reservoirs and galaxies
 * 
 * Transfers a fraction of component cq of galaxy q onto component cp of galaxy p_.
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
 *       code used to be: Gal[q].BurstMass -= 0.;
 *       now replaced by: Gal[q].BurstMass -= Mass;
 */
void transfer_stars(const int p_, const StellarComponentType cp, const int q, const StellarComponentType cq, const double fraction) 
{
  float Mass;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals Metals;
#ifdef INDIVIDUAL_ELEMENTS
  struct elements Yield;
  struct elements sfh_Elements[SFH_NBIN];
#endif
#else
  float Metals;
#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef STAR_FORMATION_HISTORY
  int i;
  float sfh_Mass[SFH_NBIN];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_Metals[SFH_NBIN];
#else
  float sfh_Metals[SFH_NBIN];
#endif
#endif

  /* Sanity checks */
  if (fraction > 1.)
  {
    printf("\n*** transfer_stars: fraction>1 ***\n");
    exit(1);
  }
#ifdef STAR_FORMATION_HISTORY
  if (Gal[p_].sfh_ibin != Gal[q].sfh_ibin)
  {
    printf("\n*** transfer_stars: inconsistent itimes ***\n");
    for(i=0;i<SFH_NBIN;i++)
    {
      printf("Bin[%d] time_1=%e dt_1=%e Nbins_1=%d time_2=%e dt_2=%e Nbins_2=%d\n",i,
                        Gal[p_].sfh_t[i],Gal[p_].sfh_dt[i],Gal[p_].sfh_Nbins[i],
                        Gal[q].sfh_t[i],Gal[q].sfh_dt[i],Gal[q].sfh_Nbins[i]);
    }
    exit(1);
  }
#endif
#ifdef TRACK_BURST
  if ((cp == BurstComponent && cq != BurstComponent) ||
      (cq == BurstComponent && cp != BurstComponent)   )
  {
    printf("\n*** transfer_stars: used incorrectly with BurstMass ***\n");
    exit(1);
  }
#endif

  //Mass and metals to be transfered and taken from galaxy q
  if (cq == DiskComponent)
  {
    Mass = fraction*Gal[q].DiskMass;
    Gal[q].DiskMass -= Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
    {
      sfh_Mass[i]             = fraction*Gal[q].sfh_DiskMass[i];
      Gal[q].sfh_DiskMass[i] -= sfh_Mass[i];
    }
#endif
    Metals = metals_fraction(Gal[q].MetalsDiskMass, fraction);
             metals_deduct_from(&Gal[q].MetalsDiskMass, Metals);
#ifdef INDIVIDUAL_ELEMENTS
    Yield  = elements_add_fraction(elements_init(),Gal[q].DiskMass_elements,fraction);
             elements_deduct_from(&Gal[q].DiskMass_elements, Yield);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
    {
      sfh_Metals[i] = metals_fraction(Gal[q].sfh_MetalsDiskMass[i],fraction);
                      metals_deduct_from(&Gal[q].sfh_MetalsDiskMass[i],sfh_Metals[i]);
#ifdef INDIVIDUAL_ELEMENTS
      sfh_Elements[i] = elements_add_fraction(elements_init(),Gal[q].sfh_ElementsDiskMass[i],fraction);
                        elements_deduct_from(&Gal[q].sfh_ElementsDiskMass[i], sfh_Elements[i]);
#endif
    }
#endif  
  }
  else if (cq == BulgeComponent)
  {
    Mass              = fraction*Gal[q].BulgeMass;
    Gal[q].BulgeMass -= Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) 
    { 
      sfh_Mass[i]              = fraction*Gal[q].sfh_BulgeMass[i];
      Gal[q].sfh_BulgeMass[i] -= sfh_Mass[i];
    }
#endif
    Metals = metals_fraction(Gal[q].MetalsBulgeMass,fraction);
             metals_deduct_from(&Gal[q].MetalsBulgeMass,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    Yield = elements_add_fraction(elements_init(),Gal[q].BulgeMass_elements,fraction);
            elements_deduct_from(&Gal[q].BulgeMass_elements,Yield);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
    {
      sfh_Metals[i] = metals_fraction(Gal[q].sfh_MetalsBulgeMass[i],fraction);
                      metals_deduct_from(&Gal[q].sfh_MetalsBulgeMass[i],sfh_Metals[i]);
#ifdef INDIVIDUAL_ELEMENTS
      sfh_Elements[i] = elements_add_fraction(elements_init(),Gal[q].sfh_ElementsBulgeMass[i],fraction);
                        elements_deduct_from(&Gal[q].sfh_ElementsBulgeMass[i],sfh_Elements[i]);
#endif
    }
#endif  
  }
  else if (cq == ICMComponent)
  {
    Mass        = fraction*Gal[q].ICM;
    Gal[q].ICM -= Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
    {
      sfh_Mass[i]=fraction*Gal[q].sfh_ICM[i];
      Gal[q].sfh_ICM[i] -= sfh_Mass[i];
    }
#endif
    Metals = metals_fraction(Gal[q].MetalsICM, fraction);
             metals_deduct_from(&Gal[q].MetalsICM,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    Yield = elements_add_fraction(elements_init(),Gal[q].ICM_elements,fraction);
            elements_deduct_from(&Gal[q].ICM_elements,Yield);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
    {
      sfh_Metals[i] = metals_fraction(Gal[q].sfh_MetalsICM[i],fraction);
                      metals_deduct_from(&Gal[q].sfh_MetalsICM[i],sfh_Metals[i]);
#ifdef INDIVIDUAL_ELEMENTS
      sfh_Elements[i] = elements_add_fraction(elements_init(),Gal[q].sfh_ElementsICM[i],fraction);
                        elements_deduct_from(&Gal[q].sfh_ElementsICM[i],sfh_Elements[i],-1.);
#endif
    }
#endif  
  }
#ifdef TRACK_BURST
  else if (cq == BurstComponent)
  {
    Mass              = fraction*Gal[q].BurstMass;
    Gal[q].BurstMass -= Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) 
    { 
      sfh_Mass[i]=fraction*Gal[q].sfh_BurstMass[i];
      Gal[q].sfh_BurstMass[i] -= sfh_Mass[i];
    }
#endif
  }
#endif
  else
  {
    printf("Unknown component type %s for galaxy q in call to transfer_stars\n", StellarComponentStr[cq]);
    exit(1);
  }

  //Add to galaxy p_
  if (cp == DiskComponent)
  {
    Gal[p_].DiskMass += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p_].sfh_ibin; i++) Gal[p_].sfh_DiskMass[i] += sfh_Mass[i];
#endif
    metals_add_to(&Gal[p_].MetalsDiskMass,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[p_].DiskMass_elements, Yield);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p_].sfh_ibin; i++)
    {
      metals_add_to(&Gal[p_].sfh_MetalsDiskMass[i],sfh_Metals[i]);
#ifdef INDIVIDUAL_ELEMENTS
      elements_add_to(&Gal[p_].sfh_ElementsDiskMass[i],sfh_Elements[i]);
#endif
    }
#endif  
  }
  else if (cp == BulgeComponent)
  {
    Gal[p_].BulgeMass += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p_].sfh_ibin; i++) 
    { Gal[p_].sfh_BulgeMass[i] += sfh_Mass[i]; }
#endif
    metals_add_to(&Gal[p_].MetalsBulgeMass,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[p_].BulgeMass_elements,Yield);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p_].sfh_ibin; i++)
    {
      metals_add_to(&Gal[p_].sfh_MetalsBulgeMass[i],sfh_Metals[i]);
#ifdef INDIVIDUAL_ELEMENTS
      elements_add_to(&Gal[p_].sfh_ElementsBulgeMass[i],sfh_Elements[i]);
#endif
    }
#endif  
  }
  else if (cp == ICMComponent)
  {
    Gal[p_].ICM += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p_].sfh_ibin; i++) { Gal[p_].sfh_ICM[i] += sfh_Mass[i]; }
#endif
    metals_add_to(&Gal[p_].MetalsICM,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[p_].ICM_elements,Yield);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p_].sfh_ibin; i++)
    {
      metals_add_to(&Gal[p_].sfh_MetalsICM[i],sfh_Metals[i]);
#ifdef INDIVIDUAL_ELEMENTS
      elements_add_to(&Gal[p_].sfh_ElementsICM[i],sfh_Elements[i]);
#endif
    }
#endif  
#ifdef TRACK_BURST
  }
  else if (cp == BurstComponent)
  {
    Gal[p_].BurstMass += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p_].sfh_ibin; i++) Gal[p_].sfh_BurstMass[i] += sfh_Mass[i];
#endif
#endif
  }
  else
  {
    printf("Unknown component type %s for galaxy p_ in call to transfer_stars\n", StellarComponentStr[cp]);
    exit(1);
  }
}


/** @brief transfer part of gas between reservoirs and galaxies
 *
 * Transfers a fraction of component cq of galaxy q onto component cp of galaxy p_.
 * cp and cq must each be one of:
 *   ColdGasComponent
 *   HotGasComponent
 *   EjectedGasComponent
 * 
 * @note used to get info about which function was calling,
 *       but this is (likely by far) the most called function,
 *       so it should be kept as tidy as possible
 */
void transfer_gas(const int p_, const GasComponentType cp, const int q, const GasComponentType  cq, const double fraction)
{
  float Mass;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals Metals;
#ifdef INDIVIDUAL_ELEMENTS
  struct elements Yield;
#endif
#else
  float Metals;
#endif

  /* Sanity check */
  if (fraction > 1.)
  {
    char sbuf[1000];
    sprintf(sbuf, "transfer_gas: fraction>1\nfraction = %.11f\nFrom '%s' to '%s\n", fraction, GasComponentStr[cq], GasComponentStr[cp]);
    terminate(sbuf);
  }

  /* Mass and Metals to be transfered taken from galaxy q */
  if (cq == ColdGasComponent)
  {
    Mass            = fraction*Gal[q].ColdGas;
    Gal[q].ColdGas -= Mass;
    Metals          = metals_fraction    (Gal[q].MetalsColdGas, fraction);
                      metals_deduct_from(&Gal[q].MetalsColdGas, Metals);
#ifdef INDIVIDUAL_ELEMENTS
    Yield           = elements_fraction    (Gal[q].ColdGas_elements,fraction);
                      elements_deduct_from(&Gal[q].ColdGas_elements,Yield);
#endif
  }
  else if (cq == HotGasComponent)
  {
    Mass           = fraction*Gal[q].HotGas;
    Gal[q].HotGas -= Mass;
    Metals         = metals_fraction    (Gal[q].MetalsHotGas,fraction);
                     metals_deduct_from(&Gal[q].MetalsHotGas,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    Yield          = elements_fraction    (Gal[q].HotGas_elements,fraction);
                     elements_deduct_from(&Gal[q].HotGas_elements,Yield);
#endif
#ifdef METALS_SELF
   metals_deduct_from(&Gal[q].MetalsHotGasSelf, Metals);
#endif
  }
  else if (cq == EjectedGasComponent)
  {
    Mass                = fraction*Gal[q].EjectedMass;
    Gal[q].EjectedMass -= Mass;
    Metals              = metals_fraction    (Gal[q].MetalsEjectedMass, fraction);
                          metals_deduct_from(&Gal[q].MetalsEjectedMass, Metals);
#ifdef INDIVIDUAL_ELEMENTS
    Yield               = elements_fraction  (Gal[q].EjectedMass_elements, fraction);
                          metals_deduct_from(&Gal[q].EjectedMass_elements, Yield);
#endif
  }
  else
  {
    printf("Unknown component type %s in call to transfer_gas\n", GasComponentStr[cq]);
    exit(1);
  }
 
 //Add to galaxy p_
  if (cp == ColdGasComponent)
  {
    Gal[p_].ColdGas += Mass;
    metals_add_to(&Gal[p_].MetalsColdGas,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[p_].ColdGas_elements,Yield);
#endif
  }
  else if (cp == HotGasComponent)
  {
    Gal[p_].HotGas += Mass;
    metals_add_to(&Gal[p_].MetalsHotGas,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[p_].HotGas_elements,Yield);
#endif
#ifdef METALS_SELF
    if (p_==q) metals_add_to(&Gal[p_].MetalsHotGasSelf,Metals);
#endif
  }
  else if (cp == EjectedGasComponent)
  {
    Gal[p_].EjectedMass += Mass;
    metals_add_to(&Gal[p_].MetalsEjectedMass,Metals);
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[p_].EjectedMass_elements,Yield);
#endif
  }
  else
  {
    printf("Unknown component type %s in call to transfer_gas\n", GasComponentStr[cp]);
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
void perform_mass_checks(char string[], const int igal)
{
  (void) string; /* suppress unused-parameter warning */
  (void) igal;   /* suppress unused-parameter warning */
  return;
  
#ifdef STAR_FORMATION_HISTORY
  int i;
  double sfh_sum;
#endif

  if(Gal[igal].ColdGas < 1.e-8)
  {
    Gal[igal].ColdGas = 0.;
    Gal[igal].MetalsColdGas = metals_init();
  }

  //check if the gas mass is less than 0
  if(Gal[igal].ColdGas < 0.0)
  {
    if (Gal[igal].ColdGas > -1e-7)
      Gal[igal].ColdGas = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, ColdGas < 0. ***\n",string);
      printf("                ColdGas[%d] = %g\n",igal,Gal[igal].ColdGas);
      terminate("");
    }
  }

  //check if the mass in metals is less than 0
  if(metals_total(Gal[igal].MetalsColdGas) < 0.0)
  {
    if (metals_total(Gal[igal].MetalsColdGas) > -1e-7)
      Gal[igal].MetalsColdGas = metals_init();
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsColdGas < 0. ***\n",string);
      printf("                MetalsColdGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsColdGas));
      terminate("");
    }
  }

  //check if the mass in metals is greater than the gas mass
  if(metals_total(Gal[igal].MetalsColdGas) > Gal[igal].ColdGas)
  {
    if (metals_total(Gal[igal].MetalsColdGas) < 1e-7)
      Gal[igal].MetalsColdGas = metals_fraction(Gal[igal].MetalsColdGas,
                                           Gal[igal].ColdGas/metals_total(Gal[igal].MetalsColdGas));
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsColdGas > ColdGas ***\n",string);
      printf("          MetalsColdGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsColdGas));
      printf("                ColdGas[%d] = %g\n",igal,Gal[igal].ColdGas);
      terminate("");
    }
  }

  if(Gal[igal].HotGas < 0.0)
  {
    if (Gal[igal].HotGas > -1e-7)
      Gal[igal].HotGas = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, HotGas < 0. ***\n",string);
      printf("                HotGas[%d] = %g\n",igal,Gal[igal].HotGas);
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsHotGas) < 0.0)
  {
    if (metals_total(Gal[igal].MetalsHotGas) > -1e-7)
      Gal[igal].MetalsHotGas = metals_init();
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsHotGas < 0. ***\n",string);
      printf("                MetalsHotGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsHotGas));
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsHotGas) > Gal[igal].HotGas)
  {
    if (metals_total(Gal[igal].MetalsHotGas) < 1e-7)
      Gal[igal].MetalsHotGas = metals_fraction(Gal[igal].MetalsHotGas,
                                           Gal[igal].HotGas/metals_total(Gal[igal].MetalsHotGas));
   else
   {
      printf("\n***  Mass check error, called from: %s, MetalsHotGas > HotGas ***\n",string);
      printf("          MetalsHotGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsHotGas));
      printf("                HotGas[%d] = %g\n",igal,Gal[igal].HotGas);
      printf("          MetalsHotGas[%d] = %.11f\n",igal,metals_total(Gal[igal].MetalsHotGas));
      printf("                HotGas[%d] = %.11f\n",igal,Gal[igal].HotGas);
      printf("             BulgeMass[%d] = %g\n",igal,Gal[igal].BulgeMass);
      printf("           EjectedMass[%d] = %g\n",igal,Gal[igal].EjectedMass);
      printf("                  Snapnum = %i\n",Gal[igal].SnapNum);
      terminate("");
    }
  }

  if(Gal[igal].EjectedMass < 0.0)
  {
    if (Gal[igal].EjectedMass > -1e-7)
      Gal[igal].EjectedMass = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, EjectedMass < 0. ***\n",string);
      printf("                EjectedMass[%d] = %g\n",igal,Gal[igal].EjectedMass);
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsEjectedMass) < 0.0)
  {
    if (metals_total(Gal[igal].MetalsEjectedMass) > -1e-7)
      Gal[igal].MetalsEjectedMass = metals_init();
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsEjectedMass < 0. ***\n",string);
      printf("                MetalsEjectedMass[%d] = %g\n",igal,metals_total(Gal[igal].MetalsEjectedMass));
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsEjectedMass) > Gal[igal].EjectedMass)
  {
    if (metals_total(Gal[igal].MetalsEjectedMass) < 1e-7)
      Gal[igal].MetalsEjectedMass = metals_fraction(Gal[igal].MetalsEjectedMass,
                                           Gal[igal].EjectedMass/metals_total(Gal[igal].MetalsEjectedMass));
    else
    {
      printf("\n*** Mass check error, called from: %s, MetalsEjectedMass > EjectedMass ***\n",string);
      printf("          MetalsEjectedMass[%d] = %g\n",igal,metals_total(Gal[igal].MetalsEjectedMass));
      printf("                EjectedMass[%d] = %g\n",igal,Gal[igal].EjectedMass);
      terminate("");
    }
  }

  if(Gal[igal].DiskMass < 0.0)
  {
    if (Gal[igal].DiskMass > -1e-7)
      Gal[igal].DiskMass = 0.;
    else 
    {
      printf("\n*** Mass check error, called from: %s, DiskMass < 0. ***\n",string);
      printf("                DiskMass[%d] = %g\n",igal,Gal[igal].DiskMass);
      terminate("");
    }
  }

  if(Gal[igal].BulgeMass < 0.0)
  {
    if (Gal[igal].BulgeMass > -1e-7)
      Gal[igal].BulgeMass = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, BulgeMass < 0. ***\n",string);
      printf("                BulgeMass[%d] = %g\n",igal,Gal[igal].BulgeMass);
      terminate("");
    }
  }

  if(Gal[igal].ICM < 0.0)
  {
    if (Gal[igal].ICM > -1e-7)
      Gal[igal].ICM = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, ICM < 0. ***\n",string);
      printf("                ICM[%d] = %g\n",igal,Gal[igal].ICM);
      terminate("");
    }
  }

#ifdef TRACK_BURST
  if(Gal[igal].BurstMass < 0.0) 
  {
    if (Gal[igal].BurstMass > -1e-7)
      Gal[igal].BurstMass = 0.;
    else
    {
      printf("\n*** Mass check error, called from: %s, BurstMass < 0. ***\n",string);
      printf("                BurstMass[%d] = %g\n",igal,Gal[igal].BurstMass);
      terminate("");
    }
  }
#endif

  /* If DETAILED_METALS_AND_MASS_RETURN, sfh stores accumulation of 'stars', not 'stars-recycFrac'.
   * Therefore, it's sum doesn't equal DiskMass any more.*/
#ifndef DETAILED_METALS_AND_MASS_RETURN
#ifdef STAR_FORMATION_HISTORY
  sfh_sum=-Gal[igal].DiskMass;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_DiskMass[i];
  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[igal].DiskMass) ||
     (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[igal].DiskMass))
  {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                DiskMass[%d] = %g\n",igal,Gal[igal].DiskMass);
    printf("            sfh_DiskMass[%d] = %g\n",igal,sfh_sum+Gal[igal].DiskMass);
    char sbuf[1000];
    sprintf(sbuf, "\n*** Mass check error, called from: %s, Inconsistent sfh for DiskMass.*** \n",string);
    terminate(sbuf);
  }

  sfh_sum=-Gal[igal].BulgeMass;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_BulgeMass[i];
  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[igal].BulgeMass) ||
     (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[igal].BulgeMass))
  {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                BulgeMass[%d] = %g\n",igal,Gal[igal].BulgeMass);
    printf("            sfh_BulgeMass[%d] = %g\n",igal,sfh_sum+Gal[igal].BulgeMass);
    char sbuf[1000];
    sprintf(sbuf, "\n*** Mass check error, called from: %s, Inconsistent sfh for BulgeMass. ***\n",string);
    terminate(sbuf);
  }

  sfh_sum=-Gal[igal].ICM;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_ICM[i];
  if(sfh_sum < -1e-4 || sfh_sum > 1e-4)
  {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                ICM[%d] = %g\n",igal,Gal[igal].ICM);
    printf("            sfh_ICM[%d] = %g\n",igal,sfh_sum+Gal[igal].ICM);
    for (i=0; i<=Gal[igal].sfh_ibin; i++)
      printf("%d %f\n",i,Gal[igal].sfh_ICM[i]);
    char sbuf[1000];
    sprintf(sbuf, "\n*** Mass check error, called from: %s, Inconsistent sfh for ICM. ***\n",string);
    terminate(sbuf);
  }
#endif //STAR_FORMATION_HISTORY
#endif //DETAILED_ENRICHEMENT

  return;
}


/** @brief print galaxy properties */
void print_galaxy(char string[], const int p_, const int halo_number_)
{
/*        printf("%s Hnr=%d firstinFOF=%d prog=%d nestprog=%d Descendant=%d gal=%d Type=%d\n",
                        string, Gal[p_].HaloNr, Halo[halo_number_].FirstHaloInFOFgroup, Halo[halo_number_].FirstProgenitor,
                        Halo[halo_number_].NextProgenitor, Halo[halo_number_].Descendant, p_, Gal[p_].Type);
        printf("     Mvir=%0.3e Vvir=%0.3e Hot=%0.3e Cold=%0.3e Eject=%0.3e disk=%0.3e bulge=%0.3e  GasDiskRadius=%0.3e\n",
                        Gal[p_].Mvir*1.e10, Gal[p_].Vvir, Gal[p_].HotGas*1.e10, Gal[p_].ColdGas*1.e10, Gal[p_].EjectedMass*1.e10,
                        Gal[p_].DiskMass*1.e10, Gal[p_].BulgeMass*1.e10, Gal[p_].GasDiskRadius);

        printf("     HotMetals=%0.3e ColdMetals=%0.3e diskMetals=%0.3e bulgeMetals=%0.3e\n",
                        Gal[p_].MetalsHotGas, Gal[p_].MetalsColdGas,Gal[p_].MetalsDiskMass, Gal[p_].MetalsBulgeMass);


        printf("     x=%0.3f y=%0.3f z=%0.3f vx=%0.3f vy=%0.3f vz=%0.3f\n",
                        Gal[p_].Pos[0],Gal[p_].Pos[1],Gal[p_].Pos[2],Gal[p_].Vel[0],Gal[p_].Vel[1],Gal[p_].Vel[2]);*/

  printf("%s Hnr=%d gal=%d Gal[gal].Hnr=%d Gal[gal].snap=%d\n",string, halo_number_, p_, Gal[p_].HaloNr, Gal[p_].SnapNum);
  printf(" Hot=%0.3e Cold=%0.3e Eject=%0.3e disk=%0.3e bulge=%0.3e  GasDiskRadius=%0.3e StellarDiskRadius=%0.3e BulgeSize=%0.3e\n",
         Gal[p_].HotGas*1.e10, Gal[p_].ColdGas*1.e10, Gal[p_].EjectedMass*1.e10,
         Gal[p_].DiskMass*1.e10, Gal[p_].BulgeMass*1.e10, Gal[p_].GasDiskRadius, Gal[p_].StellarDiskRadius, Gal[p_].BulgeSize);

  if(isnan(Gal[p_].GasDiskRadius))
  { exit(1); }
}
