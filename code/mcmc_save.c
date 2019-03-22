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

/*
 *  Created in: 2008
 *      Author: Bruno Henriques
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"

// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>

#ifdef MCMC
/** @brief Writes galaxies into a structure to be used by the MCMC */
void save_galaxy_for_mcmc(const int galaxy_index_)
{
  int output_number_, fof_number_;
  
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS
  struct GALAXY_OUTPUT galaxy_output;
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES 

#ifdef MR_PLUS_MRII
  const float log10_stellar_mass_lower_limit_ = (Switch_MR_MRII==1) ?  9.5 : 6.0;
  const float log10_stellar_mass_upper_limit_ = (Switch_MR_MRII==1) ? 14.0 : 9.5;
#else
  const float log10_stellar_mass_lower_limit_ = 7.27;
  const float log10_stellar_mass_upper_limit_ = 13.0;
#ifdef MRII
  const float log10_stellar_mass_lower_limit_ = 6.0;
  const float log10_stellar_mass_upper_limit_ = 11.27;
#endif
#endif
 
  //THE ERROR IS NOW INCLUDED IN mcmc_likelihood.c
  //StellarMass+=gsl_ran_ugaussian(MCMC_rng)*0.08*(1+MCMCConstraintsZZ[output_number_]);
#ifndef HALOMODEL
  const float log10_stellar_mass_ = log10(1E10 * (HaloGal[galaxy_index_].DiskMass + HaloGal[galaxy_index_].BulgeMass) * Hubble_h);
#else
    //no h factor in masses for OPT+=DHALOMODEL
  const float log10_stellar_mass_ = log10(1E10 * (HaloGal[galaxy_index_].DiskMass + HaloGal[galaxy_index_].BulgeMass) * inv_Hubble_h);
#endif

  if(log10_stellar_mass_lower_limit_ < log10_stellar_mass_ && log10_stellar_mass_ < log10_stellar_mass_upper_limit_)
  {
    const double    log10_Hubble_h_ = log10(Hubble_h);
    
    const int halo_number_ = HaloGal[galaxy_index_].HaloNr;
    
    // const long long first_halo_in_fof_group_number_ = HaloIDs[halo_number_].FirstHaloInFOFgroup;

    for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    {
      /*
      // linear search:
      for(fof_number_ = 0; fof_number_ < NFofsInSample[output_number_]; fof_number_++)
      {
        if(first_halo_in_fof_group_number_ == MCMC_FOF[output_number_][fof_number_].FoFID)
      */  
       
       /*
       //binary search:
      if((NFofsInSample[output_number_] > 0) && (MCMC_FOF[output_number_][0].FoFID <= first_halo_in_fof_group_number_) && (first_halo_in_fof_group_number_ <= MCMC_FOF[output_number_][NFofsInSample[output_number_] - 1].FoFID))
      {
        unsigned int fof_number_lower_bound_ = 0;
        unsigned int fof_number_upper_bound_ = NFofsInSample[output_number_] - 1;
        while(fof_number_ = (fof_number_lower_bound_ + fof_number_upper_bound_) / 2, fof_number_lower_bound_ < fof_number_upper_bound_)
        {
          if(MCMC_FOF[output_number_][fof_number_].FoFID < first_halo_in_fof_group_number_)
            fof_number_lower_bound_ = fof_number_ + 1;
          else
            fof_number_upper_bound_ = fof_number_;
        }
        if(first_halo_in_fof_group_number_ == MCMC_FOF[output_number_][fof_number_].FoFID)
          */
        if(HaloAux[halo_number_].halo_is_in_MCMC_sample_for_output[output_number_])
        {
          fof_number_ = HaloAux[halo_number_].fof_number_in_MCMC_sample_for_output[output_number_];
          
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].StellarMass   = log10_stellar_mass_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].ColdGas       = log10(1E10 * (HaloGal[galaxy_index_].ColdGas*Hubble_h));
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].BulgeMass     = log10(1E10 * HaloGal[galaxy_index_].BulgeMass*Hubble_h);
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].BlackHoleMass = log10(1E10 * HaloGal[galaxy_index_].BlackHoleMass); //black hole in units of h^-1
          //in units of Solar Masses yr^-1 h-2
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].Sfr           = log10(HaloGal[galaxy_index_].Sfr * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS) * Hubble_h * Hubble_h);

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS

          //in case of postprocess magnitudes they are only calculates here, inside prepare
          prepare_galaxy_for_output(output_number_, &HaloGal[galaxy_index_], &galaxy_output);

          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagU = galaxy_output.MagDust[0]-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagB = galaxy_output.MagDust[1]-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagV = galaxy_output.MagDust[2]-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagJ = galaxy_output.MagDust[3]-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagK = galaxy_output.MagDust[4]-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].Magu = galaxy_output.MagDust[5]-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].Magr = galaxy_output.MagDust[6]-5.*log10_Hubble_h_;
#else
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagU = lum_to_mag(HaloGal[galaxy_index_].LumDust[output_number_][0])-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagB = lum_to_mag(HaloGal[galaxy_index_].LumDust[output_number_][1])-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagV = lum_to_mag(HaloGal[galaxy_index_].LumDust[output_number_][2])-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagJ = lum_to_mag(HaloGal[galaxy_index_].LumDust[output_number_][3])-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].MagK = lum_to_mag(HaloGal[galaxy_index_].LumDust[output_number_][4])-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].Magu = lum_to_mag(HaloGal[galaxy_index_].LumDust[output_number_][5])-5.*log10_Hubble_h_;
          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].Magr = lum_to_mag(HaloGal[galaxy_index_].LumDust[output_number_][6])-5.*log10_Hubble_h_;
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
                   
#ifdef HALOMODEL
          if(output_number_ == 0)
          {
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].fofid     = fof_number_;
            //MCMC_GAL[output_number_][TotMCMCGals[output_number_]].M_Crit200 = log10(Halo[HaloGal[galaxy_index_].HaloNr].Len*PartMass*1.e10);
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].M_Crit200 = log10(Halo[HaloGal[galaxy_index_].HaloNr].M_Crit200*1.e10);
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].M_Mean200 = log10(Halo[HaloGal[galaxy_index_].HaloNr].M_Mean200*1.e10);
#ifdef MCRIT
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].M_Mean200 = log10(Halo[HaloGal[galaxy_index_].HaloNr].M_Crit200*1.e10);
#endif
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].x         = HaloGal[galaxy_index_].Pos[0];
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].y         = HaloGal[galaxy_index_].Pos[1];
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].z         = HaloGal[galaxy_index_].Pos[2];
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].Type      = HaloGal[galaxy_index_].Type;
            MCMC_GAL[output_number_][TotMCMCGals[output_number_]].ngal      = 0;
          }
#endif

          MCMC_GAL[output_number_][TotMCMCGals[output_number_]].Weight      = MCMC_FOF[output_number_][fof_number_].Weight;

#ifdef HALOMODEL
          //NOW GET PROPERTIES FOR FOF GROUPS - done for the particular fof_number_ that current galaxy resides in
          ++MCMC_FOF[output_number_][fof_number_].NGalsInFoF;

          if(HaloGal[galaxy_index_].Type==0)
          {
            MCMC_FOF[output_number_][fof_number_].IndexOfCentralGal = TotMCMCGals[output_number_];
          //MCMC_FOF[output_number_][fof_number_].M_Crit200         = log10(Halo[HaloGal[galaxy_index_].HaloNr].Len*PartMass*1.e10);
            MCMC_FOF[output_number_][fof_number_].M_Crit200         = log10(Halo[HaloGal[galaxy_index_].HaloNr].M_Crit200*1.e10);
            MCMC_FOF[output_number_][fof_number_].M_Mean200         = log10(Halo[HaloGal[galaxy_index_].HaloNr].M_Mean200*1.e10);
#ifdef MCRIT                               
            MCMC_FOF[output_number_][fof_number_].M_Mean200         = log10(Halo[HaloGal[galaxy_index_].HaloNr].M_Crit200*1.e10);
#endif
          }
#endif //HALOMODEL

          ++TotMCMCGals[output_number_];

          if(TotMCMCGals[output_number_] > MCMC_MAX_NGALS_PER_OUTPUT)
            terminate("Maximum number of galaxies in MCMC structure reached. Increase MCMCSmartFactor\n");
        }
     /* } */
    }
  }
}
#endif


