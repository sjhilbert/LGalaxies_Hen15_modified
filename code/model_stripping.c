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
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

/** @file   model_stripping.c
 *  @date   2016-2019
 *  @author ?
 *  @author Stefan Hilbert
 *
 *  @brief  models for stripping of gas and stars for satellites
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @brief This is where the gas and ICM components of newly accreted satellites
 *         are treated.
 *
 *         There are basically 2 options for the way satellite components are
 *         added into centrals:
 *
 *         HotGasStrippingModel ==1
 *         If inside Rvir, hot and ejected gas from satellites of both types 1 and 2
 *         is instantaneously striped and added to type 0.
 *
 *         HotGasStrippingModel ==0
 *         Type 1's keep an ejected component.
 *         Type 1's are stripped of hot and ejected gas gradually and later in the code.
 *         A fraction of the hot and ejected gas in the type 2's is
 *         added to the type 1 and the rest to the type 0.
 *         If satellites are outside Rvir, type 1 keeps all its components and receives
 *         everything from type 2's.
 *
 *         In these routines, central_galaxy_number_ is the type 0 at the centre of the halo;
 *         Gal[galaxy_number_].CentralGal is the galaxy around which each satellite orbits;
 *         For simplicity of reference in the comments in the code below, the latter
 *         will be called the type 1, even though it may be the same galaxy as the type 0.
 **/
void deal_with_satellites(const int central_galaxy_number_, const int n_galaxies_)
{
  int galaxy_number_, merger_centre_;
  double distance_, gas_fraction_into_type_1_, stripped_fraction_;

  for(galaxy_number_ = 0; galaxy_number_ < n_galaxies_; galaxy_number_++)    /* Loop over all galaxies in the FoF-halo */
  {
    mass_checks("Top of deal_with_satellites a",galaxy_number_);
    mass_checks("Top of deal_with_satellites b",central_galaxy_number_);
    mass_checks("Top of deal_with_satellites c",Gal[galaxy_number_].CentralGal);

    /* distance_ is the separation of the type 1 from the type 0 */
    distance_=separation_gal(central_galaxy_number_,Gal[galaxy_number_].CentralGal)/(1+ZZ[Halo[Gal[central_galaxy_number_].HaloNr].SnapNum]);


    /* HotGasStrippingModel ==  0=> Guo2010 non instantaneous treatment of gas stripping in type 1's
     *
     * if the galaxy is a type 2 and still has hot and ejected gas it is removed at this point
     * (meaning that the halo was fully stripped in previous step) and split between type 0 and type 1
     *
     * If type 2 is orbiting a type 1
     * if the type 2 is inside Rvir of the type 0  (distance_ < Gal[central_galaxy_number_].Rvir) the gas is split between 0 and 1
     * if type 2 is outside Rvir of type 0, all the gas goes to type 1
     *
     * If the type 2 is orbiting a type 0 central_galaxy_number_ and Gal[galaxy_number_].CentralGal both refer to the type 0*/

    if(HotGasStrippingModel == 0)
    {
      /* All gas Stripped from Type 2 galaxies */
      if (Gal[galaxy_number_].Type ==2)
      {
        //if type 2 is inside Rvir of type 0 split between type 0 and type 1
        if (distance_ < Gal[central_galaxy_number_].Rvir)
          gas_fraction_into_type_1_=Gal[Gal[galaxy_number_].CentralGal].HotRadius / Gal[Gal[galaxy_number_].CentralGal].Rvir;
        //if type 2 is outside Rvir of type 0, all goes to type 1
        else
          gas_fraction_into_type_1_=1.;

        Gal[galaxy_number_].HotRadius = 0.0;
        if(Gal[galaxy_number_].HotGas > 0.0)
          transfer_gas(Gal[galaxy_number_].CentralGal,HotGasComponent,galaxy_number_,HotGasComponent,gas_fraction_into_type_1_);
        if(Gal[galaxy_number_].EjectedMass > 0.0)
          transfer_gas(Gal[galaxy_number_].CentralGal,EjectedGasComponent,galaxy_number_,EjectedGasComponent,gas_fraction_into_type_1_);

        mass_checks("deal_with_satellites galaxy_number_ #0",galaxy_number_);
        mass_checks("deal_with_satellites Gal[galaxy_number_].CentraGal #0",Gal[galaxy_number_].CentralGal);
#ifdef TRACK_BURST
        /* Transfer burst component first */
        transfer_stars(Gal[galaxy_number_].CentralGal,BurstComponent,galaxy_number_,BurstComponent,
                       GasFraction_intotype1*Gal[galaxy_number_].ICM/(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass+Gal[galaxy_number_].ICM));
#endif
        transfer_stars(Gal[galaxy_number_].CentralGal,ICMComponent,galaxy_number_,ICMComponent,gas_fraction_into_type_1_);
        mass_checks("deal_with_satellites galaxy_number_ #1",galaxy_number_);
        mass_checks("deal_with_satellites Gal[galaxy_number_].CentraGal #1",Gal[galaxy_number_].CentralGal);
#ifndef POST_PROCESS_MAGS
#ifdef ICL
        transfer_ICL(Gal[galaxy_number_].CentralGal,galaxy_number_,gas_fraction_into_type_1_);
#endif
#endif
              //All the gas not moved to the type 1 yet goes to the type 0
        if (gas_fraction_into_type_1_ < 1.)
        {
          if(Gal[galaxy_number_].HotGas > 0.0)
            transfer_gas(central_galaxy_number_,HotGasComponent,galaxy_number_,HotGasComponent,1.);
          if(Gal[galaxy_number_].EjectedMass > 0.0)
            transfer_gas(central_galaxy_number_,EjectedGasComponent,galaxy_number_,EjectedGasComponent,1.);

#ifdef TRACK_BURST
          /* Transfer burst component first */
          transfer_stars(central_galaxy_number_,BurstComponent,galaxy_number_,BurstComponent,
                         Gal[galaxy_number_].ICM/(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass+Gal[galaxy_number_].ICM));
#endif
          transfer_stars(central_galaxy_number_,ICMComponent,galaxy_number_,ICMComponent,1.);
          mass_checks("deal_with_satellites #2",galaxy_number_);
          mass_checks("deal_with_satellites #2",central_galaxy_number_);
#ifndef POST_PROCESS_MAGS
#ifdef ICL
          transfer_ICL(central_galaxy_number_,galaxy_number_,1.);
#endif
#endif
        }
      }
      //Type 1 galaxies (or type 2's for the modified stripping) - stripping if galaxy inside Rvir of merger_centre_
      else if ( Gal[galaxy_number_].Type == 1 && distance_ < Gal[central_galaxy_number_].Rvir && Gal[galaxy_number_].HotGas > 0.0 )
      {
        merger_centre_ = central_galaxy_number_;

        //hot_retain_sat also re-evaluates HotRadius
        stripped_fraction_=1.-(float)(hot_retain_sat(galaxy_number_,merger_centre_))/Gal[galaxy_number_].HotGas;
        if (stripped_fraction_ < 0.)
          {
            printf("***Error in hot_retain_sat - returns value larger than HotGas***\n");
            exit(1);
          }

        transfer_gas(merger_centre_,HotGasComponent,galaxy_number_,HotGasComponent,stripped_fraction_);
        transfer_gas(merger_centre_,EjectedGasComponent,galaxy_number_,EjectedGasComponent,stripped_fraction_);
        mass_checks("deal_with_satellites #3",galaxy_number_);
        mass_checks("deal_with_satellites #3",merger_centre_);
#ifdef TRACK_BURST
        /* Transfer burst component first */
        transfer_stars(merger_centre_,BurstComponent,galaxy_number_,BurstComponent,
                       stripped_fraction_*Gal[galaxy_number_].ICM/(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass+Gal[galaxy_number_].ICM));
#endif
        transfer_stars(merger_centre_,ICMComponent,galaxy_number_,ICMComponent,stripped_fraction_);
        mass_checks("deal_with_satellites #4",galaxy_number_);
        mass_checks("deal_with_satellites #4",merger_centre_);
#ifndef POST_PROCESS_MAGS
#ifdef ICL
        transfer_ICL(merger_centre_,galaxy_number_,stripped_fraction_);
#endif
#endif
      }

      mass_checks("deal_with_satellites #5",galaxy_number_);
      mass_checks("deal_with_satellites #5",central_galaxy_number_);
    }
      /* Instantaneous stripping of gas from satellites and no ejection of type 2 into type 1,
       * still there is the condition on Rvir that determines that if a galaxy is a newly
       * accreted type 2 outside Rvir of type 0, its gas will go into the type 1. If it's
       * a type 1 outside Rvir of type 0, it retains all its gas. -> DeLucia2007*/
    else if (HotGasStrippingModel == 1)
    {
    /* If galaxy is a satellite inside Rvir it will lose its hot and
     * ejected gas into the hot gas component of the central_galaxy_number_.
     * Only galaxies within Rvir contribute to the central halo.*/
      if ( distance_ < Gal[central_galaxy_number_].Rvir && galaxy_number_ != central_galaxy_number_)
      {
        transfer_gas(central_galaxy_number_,HotGasComponent,galaxy_number_,HotGasComponent,1.);
        transfer_gas(central_galaxy_number_,EjectedGasComponent,galaxy_number_,EjectedGasComponent,1.);
#ifdef TRACK_BURST
        /* Transfer burst component first */
        transfer_stars(central_galaxy_number_,BurstComponent,galaxy_number_,BurstComponent,
                       Gal[galaxy_number_].ICM/(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass+Gal[galaxy_number_].ICM));
#endif
        transfer_stars(central_galaxy_number_,ICMComponent,galaxy_number_,ICMComponent,1.);
#ifndef POST_PROCESS_MAGS
#ifdef ICL
        transfer_ICL(central_galaxy_number_, galaxy_number_, 1.);
#endif
#endif
        Gal[galaxy_number_].HotRadius =0.;
      }
      /* If its a type 1 outside Rvir it retains all its gas components, so do nothing
       * else if (Gal[galaxy_number_].Type ==1) {}
       * If galaxy is a type 2 outside Rvir of type 0, then all its gas components
       * will be added to the type 1. */
      else  if (Gal[galaxy_number_].Type == 2)
      {
        transfer_gas(Gal[galaxy_number_].CentralGal,HotGasComponent,galaxy_number_,HotGasComponent,1.);
        transfer_gas(Gal[galaxy_number_].CentralGal,EjectedGasComponent,galaxy_number_,EjectedGasComponent,1.);
#ifdef TRACK_BURST
        /* Transfer burst component first */
        transfer_stars(Gal[galaxy_number_].CentralGal,BurstComponent,galaxy_number_,BurstComponent,
                       Gal[galaxy_number_].ICM/(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass+Gal[galaxy_number_].ICM));
#endif
        transfer_stars(Gal[galaxy_number_].CentralGal,ICMComponent,galaxy_number_,ICMComponent,1.);
#ifndef POST_PROCESS_MAGS
#ifdef ICL
        transfer_ICL(Gal[galaxy_number_].CentralGal,galaxy_number_,1.);
#endif
#endif
        Gal[galaxy_number_].HotRadius =0.;
      }
    }//end of HotGasStrippingModel == 1

    mass_checks("Bottom of deal_with_satellites galaxy_number_",galaxy_number_);
    mass_checks("Bottom of deal_with_satellites central_galaxy_number_",central_galaxy_number_);

  } /* End of HotGasStrippingModel choice */

  return;
}


/** Gradual stripping of hot and ejected gas from type 1 satellites. 
 *  This is caused both by tidal and ram-pressure stripping.
 *  This function returns the actual mass of hot gas that the
 *  type 1 retains.
 *
 *  TIDAL STRIPPING
 *  Hot gas is tidally stripped at the same rate at which dark matter is
 *  stripped:
 *
 * \f$ \frac{M_{\rm{hot}}(R_{\rm{tidal}})}{M_{\rm{hot,infall}}}=
 *  \frac{M_{\rm{DM}}}{M_{\rm{DM,infall}}}\f$
 *
 *  Since the hot gas distribution is assumed to be \f$ \rho \propto r^{-2}\f$
 *  this means \f$ M_{\rm{hot}}(r) \propto r.\f$ Therefore, the tidal
 *  radius beyond gas is stripped is given by:
 *
 *  \f$ R_{\rm{tidal}}=
 *  \left(\frac{M_{\rm{DM}}}{M_{\rm{DM,infall}}}\right)R_{\rm{DM,infall}}\f$
 *
 *  RAM PRESSURE STRIPING
 *  Let \f$R_{r.p.}\f$ represent the distance from the centre of the satellite
 *  at which ram pressure striping equals its self-gravity. Then:
 *
 *  \f$ \rho_{\rm{sat}}(R_{\rm{r.p.}})V^2_{\rm{sat}}=
 *      \rho_{\rm{par}}(R_{\rm{orbit}})V^2_{\rm{orbit}}\f$
 *  Where the four terms represent respectively the density of the satellite
 *  at \f$R_{\rm{r.p.}}\f$, the virial velocity of the satellite at infall,
 *  the density of the parent halo at the radius of the satellite and the
 *  orbit velocity of the satellite (given by \f$V_{\rm{c}} of the parent halo\f$)
 *
 *  The stripping radius is given by
 *
 *  \f$R_{\rm{strip}}=min(R_{\rm{tidal}},R_{\rm{r.p.}})\f$
 *
 * */
double hot_retain_sat(const int galaxy_number_, const int central_galaxy_number_)
{
  double hot_remain_;
  double R_stripping_, R_tidal_, R_ram_pressure_, R_orbit_, total_mass_sat_, V_orbit_;

  if (Gal[central_galaxy_number_].Type != 0)
    exit(0);

  /*Calculate tidal stripping radius*/
   R_tidal_=Gal[galaxy_number_].Len*PartMass/Gal[galaxy_number_].Mvir*Gal[galaxy_number_].Rvir;

  /*Ram pressure stripping radius calculation*/

  /*First calculate the orbital radius of the satellite R_orbit*/
   R_orbit_=separation_gal(central_galaxy_number_,galaxy_number_)/(1+ZZ[Halo[Gal[central_galaxy_number_].HaloNr].SnapNum]);


  /*If the central galaxy has no hot gas, it exerts no ram pressure stripping on the
   * satellite. */
  if (Gal[central_galaxy_number_].HotGas<1.e-6 || Gal[central_galaxy_number_].Mvir<RamPressureStrip_CutOffMass)
          R_ram_pressure_=Gal[galaxy_number_].HotRadius;
  else
  {
    total_mass_sat_=Gal[galaxy_number_].Mvir;
    V_orbit_=sqrt((Gravity*Gal[central_galaxy_number_].Mvir)/Gal[central_galaxy_number_].Rvir);
    R_ram_pressure_= sqrt(Gal[galaxy_number_].HotGas/Gal[galaxy_number_].HotRadius) * sqrt(Gravity * total_mass_sat_/Gal[galaxy_number_].Rvir) *
        sqrt(Gal[central_galaxy_number_].Rvir/Gal[central_galaxy_number_].HotGas)*R_orbit_ * 1./V_orbit_;
  }

  /*Get the smaller of tidal and ram pressure stripping radii.*/
  R_stripping_=min(R_tidal_, R_ram_pressure_);

  /*if the stripping radius is larger then hot radius there is
   * no stripping*/
  if (R_stripping_>Gal[galaxy_number_].HotRadius || Gal[galaxy_number_].HotGas < 1.e-8)
    hot_remain_=Gal[galaxy_number_].HotGas;         
  // If stripping radius is smaller than the hot radius
  else 
  {
    //Assuming M_hot(r) proportional to r, the remaining hot gas is given by:
    hot_remain_=Gal[galaxy_number_].HotGas*R_stripping_/Gal[galaxy_number_].HotRadius;
    // hot radius is updated to the stripping radius
    Gal[galaxy_number_].HotRadius=R_stripping_;

    // Check that HotRadius has sensible values
    if (Gal[galaxy_number_].HotRadius < 1.e-8)
      Gal[galaxy_number_].HotRadius = Gal[galaxy_number_].Len*PartMass/Gal[galaxy_number_].Mvir*Gal[galaxy_number_].Rvir;
    if (Gal[galaxy_number_].HotRadius > Gal[galaxy_number_].Rvir)          
      Gal[galaxy_number_].HotRadius = Gal[galaxy_number_].Rvir;
  }

  if(hot_remain_>Gal[galaxy_number_].HotGas)
    hot_remain_=Gal[galaxy_number_].HotGas;

  return hot_remain_;
}
