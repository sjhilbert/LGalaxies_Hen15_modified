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

/**@file save.c
 * @brief Copies the relevant properties in Galaxy structure into
 *        Galaxy_Output structure and saves them into the output
 *        files (SA_z**_**) - redshift/file_number_.
 *        */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

#ifdef OUTPUT_BUFFERING
#include "dynamic_array_inline.h"
#endif /* defined OUTPUT_BUFFERING */


#ifdef OUTPUT_BUFFERING

/** initial buffer capacity */
#ifndef OUTPUT_BUFFERING_INIT_CAPACITY
#if OUTPUT_BUFFERING == 1
#define OUTPUT_BUFFERING_INIT_CAPACITY 10 * BYTES_PER_MB
#else /* OUTPUT_BUFFERING == 2 */
#define OUTPUT_BUFFERING_INIT_CAPACITY 10 * BYTES_PER_MB
#endif /* OUTPUT_BUFFERING == 2 */
#endif /* not defined OUTPUT_BUFFERING_INIT_CAPACITY_IN_MB */


/** @brief the buffer used to buffer galaxy output to disk */
static dynamic_array_type galaxy_output_buffer[NOUT];


/** @brief init the buffer used to buffer galaxy output to disk */
void save_galaxy_init_output_buffer(void)
{ 
  int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; ++output_number_)
  { dynamic_array_init(&galaxy_output_buffer[output_number_], 0, OUTPUT_BUFFERING_INIT_CAPACITY); }
}


/** @brief writes output buffer content to file and clears buffer  */
void save_galaxy_flush_output_buffer(void)
{
  int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    if(galaxy_output_buffer[output_number_].size > 0)
    {
      myfwrite_large_data(galaxy_output_buffer[output_number_].data, 1, galaxy_output_buffer[output_number_].size, FdGalDumps[output_number_]);
      dynamic_array_clear(&galaxy_output_buffer[output_number_]);
    }
}


/** @brief shows output buffer stats */
void save_galaxy_show_output_buffer_statistics(void)
{
  int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; ++output_number_)
  { printf("galaxy output buffer[%d]: capacity = %lu (%f MB)\n", output_number_, galaxy_output_buffer[output_number_].capacity, galaxy_output_buffer[output_number_].capacity / (1024. * 1024.)); }
}

#endif /* defined OUTPUT_BUFFERING */


 /** @brief create galaxy output files for redshifts */
void create_galaxy_files(const int file_number_)
{
  // create output files - snapshot option
  int output_number_, tree_number_;
  char file_name_[1000];

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    for(tree_number_ = 0; tree_number_ < Ntrees; tree_number_++)
    { TreeNgals[output_number_][tree_number_] = 0; }

    sprintf(file_name_, "%s/%s_z%1.2f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[output_number_]], file_number_);
    if(!(FdGalDumps[output_number_] = fopen(file_name_, "wb+")))
    {
      char error_message_[1000];
      sprintf(error_message_, "can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }

    fseek(FdGalDumps[output_number_], (2 + Ntrees) * sizeof(int), SEEK_SET);        /* skip the space for the header */
    TotGalaxies[output_number_] = 0;
    
#ifdef OUTPUT_BUFFERING
    dynamic_array_clear(&galaxy_output_buffer[output_number_]);
#endif /* defined OUTPUT_BUFFERING */     
  }
}


 /** @brief write header and close galaxy output files */
void close_galaxy_files(void)
{
  int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    fseek(FdGalDumps[output_number_], 0, SEEK_SET);
    myfwrite(&Ntrees, sizeof(int), 1, FdGalDumps[output_number_]);        //Number of trees
    myfwrite(&TotGalaxies[output_number_], sizeof(int), 1, FdGalDumps[output_number_]);        // total number of galaxies
    myfwrite(TreeNgals[output_number_], sizeof(int), Ntrees, FdGalDumps[output_number_]);        // Number of galaxies in each tree
#ifdef OUTPUT_BUFFERING
#if OUTPUT_BUFFERING == 2
    if(galaxy_output_buffer[output_number_].size > 0)
    {
      myfwrite_large_data(galaxy_output_buffer[output_number_].data, galaxy_output_buffer[output_number_].size, 1, FdGalDumps[output_number_]);
      dynamic_array_clear(&galaxy_output_buffer[output_number_]);
    }
#endif /* OUTPUT_BUFFERING == 2 */
#endif /* defined OUTPUT_BUFFERING */    
    fclose(FdGalDumps[output_number_]);
  }
}


/** @brief Saves the Galaxy_Output structure for all the galaxies in
 *        the current tree into the current output file (one for each
 *        input dark matter file) for the chosen snapshots.
 *
 *        If UPDATETYPETWO=1 then the positions and velocities of type 2
 *        galaxies are updated from the most bound dark matter particle.
 *        After that the GALAXY_OUTPUT structure is created and written.
 *        input: int file number (current file where the output is
 *        being written), int tree number (tree being currently treated).
 */
void save_galaxy_append(const int tree_number_, const int galaxy_number_, const int output_number_)
{
  struct GALAXY_OUTPUT galaxy_output_;
  prepare_galaxy_for_output(output_number_, &HaloGal[galaxy_number_], &galaxy_output_);
  
#ifdef OUTPUT_BUFFERING
  dynamic_array_push_back(&galaxy_output_buffer[output_number_], &galaxy_output_, sizeof(struct GALAXY_OUTPUT));
#else  /* not defined OUTPUT_BUFFERING */
  myfwrite(&galaxy_output_, sizeof(struct GALAXY_OUTPUT), 1, FdGalDumps[output_number_]);
#endif /* not defined OUTPUT_BUFFERING */  

  ++TotGalaxies[output_number_];                 //this will be written later
  ++TreeNgals  [output_number_][tree_number_];   //this will be written later (Number of galaxies in each tree)
}


/** @brief Copies all the relevant properties from the Galaxy structure
 *         into the Galaxy output structure, some units are corrected.
 *
 * @bug (corrected by Stefan Hilbert) #ifndef GALAXYTREE, DisruptOn was not set,
 *       but left to garbage. Now it will be set to zero in that case.
 *      
 * @bug (corrected by Stefan Hilbert) sfh_numbins was not set. it only was
 *       set in save_galaxy_tree_append(), where it was set to sfh_ibin. The 
 *       correct value, however, should be (sfh_ibin + 1), and it should be
 *       set here to work for all kinds of output (snapshot, tree, lightcone).
**/
void prepare_galaxy_for_output(const int output_number_, const struct GALAXY *galaxy_, struct GALAXY_OUTPUT *output_galaxy_)
{
  int j_;

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#if defined OUTPUT_REST_MAGS || defined OUTPUT_OBS_MAGS || defined OUTPUT_FB_OBS_MAGS
  int filter_number_;
#endif /* defined OUTPUT_REST_MAGS || defined OUTPUT_OBS_MAGS || defined OUTPUT_FB_OBS_MAGS */
#ifndef LIGHT_OUTPUT
#ifdef OUTPUT_REST_MAGS
  const int r_band_filter_number_ = 17;
#endif /* defined OUTPUT_REST_MAGS */
#endif /* not defined LIGHT_OUTPUT */
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
 
#ifdef LIGHT_OUTPUT 
  output_galaxy_->Type                       = galaxy_->Type;
  output_galaxy_->SnapNum                    = galaxy_->SnapNum;
  output_galaxy_->CentralMvir                = get_virial_mass  (Halo[galaxy_->HaloNr].FirstHaloInFOFgroup);
  output_galaxy_->CentralRvir                = get_virial_radius(Halo[galaxy_->HaloNr].FirstHaloInFOFgroup);

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->Pos[j_]                    = galaxy_->Pos[j_];

  output_galaxy_->Mvir                       = galaxy_->Mvir;
  output_galaxy_->Rvir                       = galaxy_->Rvir;
  output_galaxy_->Vvir                       = galaxy_->Vvir;

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->DistanceToCentralGal[j_]   = wrap(Halo[Halo[galaxy_->HaloNr].FirstHaloInFOFgroup].Pos[j_] - galaxy_->Pos[j_], BoxSize);

  output_galaxy_->ColdGas                    = galaxy_->ColdGas;
  output_galaxy_->BulgeMass                  = galaxy_->BulgeMass;
  output_galaxy_->DiskMass                   = galaxy_->DiskMass;
  output_galaxy_->HotGas                     = galaxy_->HotGas;
  output_galaxy_->BlackHoleMass              = galaxy_->BlackHoleMass;

  output_galaxy_->GasDiskRadius              = galaxy_->GasDiskRadius;
  output_galaxy_->CosInclination             = galaxy_->CosInclination;

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS
  post_process_spec_mags(output_galaxy_);
#else /* not defined POST_PROCESS_MAGS */
#ifdef OUTPUT_OBS_MAGS
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->ObsMagDust[filter_number_] = lum_to_mag(galaxy_->ObsLumDust[output_number_][filter_number_]);
#endif /* defined OUTPUT_OBS_MAGS */
#ifdef OUTPUT_REST_MAGS
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->MagDust[filter_number_]    = lum_to_mag(galaxy_->LumDust[output_number_][filter_number_]);
#endif /* defined OUTPUT_REST_MAGS */ 
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
 
#else  /* not defined LIGHT_OUTPUT */

#ifdef GALAXYTREE
  // for calculating MMSubID:
  int tmpfirst_ = Halo[galaxy_->HaloNr].FirstHaloInFOFgroup;
  int lenmax_ = 0;
  int next_ = tmpfirst_;
  while(next_ != -1)
  {
    if(Halo[next_].Len > lenmax_)
    {
      lenmax_ = Halo[next_].Len;
      tmpfirst_ = next_;
    }
    next_ = Halo[next_].NextHaloInFOFgroup;
  }
#endif /* defined GALAXYTREE */

#ifdef GALAXYTREE
#ifdef MCMC
//   output_galaxy_->GalID                                = -1;
//   output_galaxy_->HaloID                               = -1;
#else /* not defined MCMC */
//   output_galaxy_->GalID                                = -1;
  output_galaxy_->HaloID                               = HaloIDs[galaxy_->HaloNr].HaloID;
#endif /* not defined MCMC */ 
#endif /* defined GALAXYTREE */

#ifdef MBPID
  output_galaxy_->MostBoundID                          =  galaxy_->MostBoundID;
#endif /* defined MBPID */

#ifdef GALAXYTREE
//   output_galaxy_->FirstProgGal                         = -1;
//   output_galaxy_->NextProgGal                          = -1;
//   output_galaxy_->LastProgGal                          = -1;
//   output_galaxy_->FOFCentralGal                        = -1;
//   output_galaxy_->FileTreeNr                           = -1;
//   output_galaxy_->DescendantGal                        = -1;
//   output_galaxy_->MainLeafId                           = -1;
//   output_galaxy_->TreeRootId                           = -1;
  output_galaxy_->SubID                                = calc_big_db_subid_index(galaxy_->SnapNum, Halo[galaxy_->HaloNr].FileNr, Halo[galaxy_->HaloNr].SubhaloIndex);
  output_galaxy_->MMSubID                              = calc_big_db_subid_index(galaxy_->SnapNum, Halo[tmpfirst_].FileNr, Halo[tmpfirst_].SubhaloIndex);
  output_galaxy_->PeanoKey                             = peano_hilbert_key((int) floor(galaxy_->Pos[0] * ScaleFactor), (int) floor(galaxy_->Pos[1] * ScaleFactor), (int) floor(galaxy_->Pos[2] * ScaleFactor), Hashbits);
#endif /* defined GALAXYTREE */

  output_galaxy_->Redshift                             = ZZ[galaxy_->SnapNum];
#ifdef LIGHTCONE_OUTPUT
  output_galaxy_->ObsRedshift                          = ZZ[galaxy_->SnapNum];
  output_galaxy_->CubeShiftIndex                       = 0;
#endif /* defined LIGHTCONE_OUTPUT */

  output_galaxy_->Type                                 = galaxy_->Type;
#ifndef GALAXYTREE
  output_galaxy_->HaloIndex                            = galaxy_->HaloNr;
#endif /* not defined GALAXYTREE */

#ifdef HALOPROPERTIES
  output_galaxy_->HaloM_Mean200                        = galaxy_->HaloM_Mean200;
  output_galaxy_->HaloM_Crit200                        = galaxy_->HaloM_Crit200;
  output_galaxy_->HaloM_TopHat                         = galaxy_->HaloM_TopHat;

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->HaloPos[j_]                          = galaxy_->HaloPos[j_];

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->HaloVel[j_]                          = galaxy_->HaloVel[j_];

  output_galaxy_->HaloVelDisp                          = galaxy_->HaloVelDisp;
  output_galaxy_->HaloVmax                             = galaxy_->HaloVmax;

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->HaloSpin[j_]                         = galaxy_->HaloSpin[j_];
#endif /* defined HALOPROPERTIES */

  output_galaxy_->SnapNum                              = galaxy_->SnapNum;
  output_galaxy_->LookBackTimeToSnap                   = NumToTime(galaxy_->SnapNum) * UnitTime_in_years * inv_Hubble_h;
  output_galaxy_->CentralMvir                          = get_virial_mass(Halo[galaxy_->HaloNr].FirstHaloInFOFgroup);
  output_galaxy_->CentralRvir                          = get_virial_radius(Halo[galaxy_->HaloNr].FirstHaloInFOFgroup);

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->DistanceToCentralGal[j_]             = wrap(Halo[Halo[galaxy_->HaloNr].FirstHaloInFOFgroup].Pos[j_] - galaxy_->Pos[j_],BoxSize);

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->Pos[j_]                              = galaxy_->Pos[j_];

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->Vel[j_]                              = galaxy_->Vel[j_];

  output_galaxy_->Len                                  = galaxy_->Len;
  output_galaxy_->Mvir                                 = galaxy_->Mvir;
  output_galaxy_->Rvir                                 = galaxy_->Rvir;
  output_galaxy_->Vvir                                 = galaxy_->Vvir;
  output_galaxy_->Vmax                                 = galaxy_->Vmax;

  for(j_ = 0; j_ < 3; ++j_)
  output_galaxy_->GasSpin[j_]                          = galaxy_->GasSpin[j_];

  for(j_ = 0; j_ < 3; ++j_)                              
  output_galaxy_->StellarSpin[j_]                      = galaxy_->StellarSpin[j_];

  output_galaxy_->InfallVmax                           = galaxy_->InfallVmax;
  output_galaxy_->InfallVmaxPeak                       = galaxy_->InfallVmaxPeak;
  output_galaxy_->InfallSnap                           = galaxy_->InfallSnap;
  output_galaxy_->InfallHotGas                         = galaxy_-> InfallHotGas;
  output_galaxy_->HotRadius                            = galaxy_->HotRadius;
  output_galaxy_->OriMergTime                          = (galaxy_->Type == 2 || (galaxy_->Type == 1 && galaxy_->MergeOn == 1)) ? galaxy_->OriMergTime * UnitTime_in_years * inv_Hubble_h : 0.;
  output_galaxy_->MergTime                             = (galaxy_->Type == 2 || (galaxy_->Type == 1 && galaxy_->MergeOn == 1)) ? galaxy_->MergTime    * UnitTime_in_years * inv_Hubble_h : 0.;
  output_galaxy_->ColdGas                              = galaxy_->ColdGas;
  output_galaxy_->StellarMass                          = galaxy_->BulgeMass + galaxy_->DiskMass;
  output_galaxy_->BulgeMass                            = galaxy_->BulgeMass;
  output_galaxy_->DiskMass                             = galaxy_->DiskMass;
  output_galaxy_->HotGas                               = galaxy_->HotGas;
  output_galaxy_->EjectedMass                          = CORRECTDBFLOAT(galaxy_->EjectedMass);
  output_galaxy_->BlackHoleMass                        = galaxy_->BlackHoleMass;
  output_galaxy_->ICM                                  = galaxy_->ICM;
#ifdef DETAILED_METALS_AND_MASS_RETURN                  
  output_galaxy_->MetalsColdGas                        = galaxy_->MetalsColdGas;
  output_galaxy_->MetalsBulgeMass                      = galaxy_->MetalsBulgeMass;
  output_galaxy_->MetalsDiskMass                       = galaxy_->MetalsDiskMass;
  output_galaxy_->MetalsHotGas                         = galaxy_->MetalsHotGas;
  output_galaxy_->MetalsEjectedMass                    = galaxy_->MetalsEjectedMass;
  output_galaxy_->MetalsICM                            = galaxy_->MetalsICM;
#ifdef METALS_SELF                                       
  output_galaxy_->MetalsHotGasSelf                     = galaxy_->MetalsHotGasSelf;
#endif /* defined METALS_SELF */
#else /* not defined DETAILED_METALS_AND_MASS_RETURN */
  output_galaxy_->MetalsColdGas                        = CORRECTDBFLOAT(galaxy_->MetalsColdGas);
  output_galaxy_->MetalsStellarMass                    = CORRECTDBFLOAT(galaxy_->MetalsDiskMass) + CORRECTDBFLOAT(galaxy_->MetalsBulgeMass);
  output_galaxy_->MetalsBulgeMass                      = CORRECTDBFLOAT(galaxy_->MetalsBulgeMass);
  output_galaxy_->MetalsDiskMass                       = CORRECTDBFLOAT(galaxy_->MetalsDiskMass);
  output_galaxy_->MetalsHotGas                         = CORRECTDBFLOAT(galaxy_->MetalsHotGas);
  output_galaxy_->MetalsEjectedMass                    = CORRECTDBFLOAT(galaxy_->MetalsEjectedMass);
  output_galaxy_->MetalsICM                            = CORRECTDBFLOAT(galaxy_->MetalsICM);
#ifdef METALS_SELF                                       
  output_galaxy_->MetalsHotGasSelf                     = CORRECTDBFLOAT(galaxy_->MetalsHotGasSelf);
#endif /* defined METALS_SELF */                         
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */
#ifdef TRACK_BURST 
  output_galaxy_->BurstMass                            = galaxy_->BurstMass;
#endif /* defined TRACK_BURST */                         
  output_galaxy_->PrimordialAccretionRate              = CORRECTDBFLOAT(galaxy_->PrimordialAccretionRate * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS));
  output_galaxy_->CoolingRadius                        =                galaxy_->CoolingRadius;
  output_galaxy_->CoolingRate                          = CORRECTDBFLOAT(galaxy_->CoolingRate             * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS));
  output_galaxy_->CoolingRate_beforeAGN                = CORRECTDBFLOAT(galaxy_->CoolingRate_beforeAGN   * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)); 
  output_galaxy_->QuasarAccretionRate                  =                galaxy_->QuasarAccretionRate     * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  output_galaxy_->RadioAccretionRate                   =                galaxy_->RadioAccretionRate      * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  output_galaxy_->Sfr                                  = CORRECTDBFLOAT(galaxy_->Sfr                     * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS));
  output_galaxy_->SfrBulge                             = CORRECTDBFLOAT(galaxy_->SfrBulge                * (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS));
  output_galaxy_->XrayLum                              = galaxy_->XrayLum;
  output_galaxy_->BulgeSize                            = galaxy_->BulgeSize;
  output_galaxy_->StellarDiskRadius                    = galaxy_->StellarDiskRadius;
  output_galaxy_->GasDiskRadius                        = galaxy_->GasDiskRadius;
  output_galaxy_->CosInclination                       = galaxy_->CosInclination;
  
#ifdef GALAXYTREE
  output_galaxy_->DisruptOn                            = galaxy_->DisruptOn;
#else  /* not defined GALAXYTREE */
  output_galaxy_->DisruptOn                            = 0;
#endif /* not defined GALAXYTREE */

  output_galaxy_->MergeOn                              = galaxy_->MergeOn;
  
  output_galaxy_->MassWeightAge                        = (galaxy_->DiskMass + galaxy_->BulgeMass > 0.) ? galaxy_->MassWeightAge[output_number_] / (galaxy_->DiskMass + galaxy_->BulgeMass) * UnitTime_in_Gigayears * inv_Hubble_h : 0.;
  
#ifdef STAR_FORMATION_HISTORY
  output_galaxy_->sfh_ibin                             = galaxy_->sfh_ibin;
  output_galaxy_->sfh_numbins                          = galaxy_->sfh_ibin + 1;

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_DiskMass[j_]                     = galaxy_->sfh_DiskMass[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_DiskMass[j_]                     = 0.;

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_BulgeMass[j_]                    = galaxy_->sfh_BulgeMass[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_BulgeMass[j_]                    = 0.;

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_ICM[j_]                          = galaxy_->sfh_ICM[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_ICM[j_]                          = 0.;

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_MetalsDiskMass[j_]               = galaxy_->sfh_MetalsDiskMass[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_MetalsDiskMass[j_]               = metals_init();

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_MetalsBulgeMass[j_]              = galaxy_->sfh_MetalsBulgeMass[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_MetalsBulgeMass[j_]              = metals_init();

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_MetalsICM[j_]                    = galaxy_->sfh_MetalsICM[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_MetalsICM[j_]                    = metals_init();

#ifdef TRACK_BURST
  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_BurstMass[j_]                    = galaxy_->sfh_BurstMass[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_BurstMass[j_]                    = 0.;
#endif /* defined TRACK_BURST */

#endif /* defined STAR_FORMATION_HISTORY */

#ifdef INDIVIDUAL_ELEMENTS
  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_ElementsDiskMass[j_]             = galaxy_->sfh_ElementsDiskMass[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_ElementsDiskMass[j_]             = elements_init();

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_ElementsBulgeMass[j_]            = galaxy_->sfh_ElementsBulgeMass[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_ElementsBulgeMass[j_]            = elements_init();

  for (j_ = 0; j_ <= galaxy_->sfh_ibin; j_++)
  output_galaxy_->sfh_ElementsICM[j_]                  = galaxy_->sfh_ElementsICM[j_];
  for (j_ = galaxy_->sfh_ibin + 1; j_ < SFH_NBIN; j_++)
  output_galaxy_->sfh_ElementsICM[j_]                  = elements_init();

  output_galaxy_->DiskMass_elements                    = galaxy_->DiskMass_elements;
  output_galaxy_->BulgeMass_elements                   = galaxy_->BulgeMass_elements;
  output_galaxy_->ColdGas_elements                     = galaxy_->ColdGas_elements;
  output_galaxy_->HotGas_elements                      = galaxy_->HotGas_elements;
  output_galaxy_->ICM_elements                         = galaxy_->ICM_elements;
  output_galaxy_->EjectedMass_elements                 = galaxy_->EjectedMass_elements;
#endif /* defined INDIVIDUAL_ELEMENTS */
  
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS
  post_process_spec_mags(output_galaxy_);
#else  /* not defined POST_PROCESS_MAGS */

#ifdef OUTPUT_REST_MAGS
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->MagDust             [filter_number_] = lum_to_mag(galaxy_->LumDust             [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->Mag                 [filter_number_] = lum_to_mag(galaxy_->Lum                 [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->MagBulge            [filter_number_] = lum_to_mag(galaxy_->LumBulge            [output_number_][filter_number_]);
#ifdef ICL
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->MagICL              [filter_number_] = lum_to_mag(galaxy_->LumICL              [output_number_][filter_number_]);
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS 
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->ObsMagDust          [filter_number_] = lum_to_mag(galaxy_->ObsLumDust          [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->ObsMag              [filter_number_] = lum_to_mag(galaxy_->ObsLum              [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->ObsMagBulge         [filter_number_] = lum_to_mag(galaxy_->ObsLumBulge         [output_number_][filter_number_]);
#ifdef ICL
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->ObsMagICL           [filter_number_] = lum_to_mag(galaxy_->ObsLumICL           [output_number_][filter_number_]);
#endif /* defined ICL */
#ifdef OUTPUT_FB_OBS_MAGS
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->backward_ObsMagDust [filter_number_] = lum_to_mag(galaxy_->backward_ObsLumDust [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->backward_ObsMag     [filter_number_] = lum_to_mag(galaxy_->backward_ObsLum     [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->backward_ObsMagBulge[filter_number_] = lum_to_mag(galaxy_->backward_ObsLumBulge[output_number_][filter_number_]);
#ifdef ICL 
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->backward_ObsMagICL  [filter_number_] = lum_to_mag(galaxy_->backward_ObsLumICL  [output_number_][filter_number_]);
#endif /* defined ICL */
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->forward_ObsMagDust  [filter_number_] = lum_to_mag(galaxy_->forward_ObsLumDust  [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->forward_ObsMag      [filter_number_] = lum_to_mag(galaxy_->forward_ObsLum      [output_number_][filter_number_]);
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->forward_ObsMagBulge [filter_number_] = lum_to_mag(galaxy_->forward_ObsLumBulge [output_number_][filter_number_]);
#ifdef ICL
  for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
  output_galaxy_->forward_ObsMagICL   [filter_number_] = lum_to_mag(galaxy_->forward_ObsLumICL   [output_number_][filter_number_]);
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */

#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
  output_galaxy_->rbandWeightAge                       = (galaxy_->DiskMass + galaxy_->BulgeMass > 0.) ? galaxy_->rbandWeightAge[output_number_] / (galaxy_->Lum[output_number_][r_band_filter_number_]) * UnitTime_in_Gigayears * inv_Hubble_h;
#else  /* not defined OUTPUT_REST_MAGS */ 
  output_galaxy_->rbandWeightAge                       = 0.; 
#endif /* not defined OUTPUT_REST_MAGS */ 
#endif /* not defined POST_PROCESS_MAGS */
#else  /* not defined COMPUTE_SPECPHOT_PROPERTIES */
  output_galaxy_->rbandWeightAge                       = 0.; 
#endif /* not defined COMPUTE_SPECPHOT_PROPERTIES */
 
#endif /* not defined LIGHT_OUTPUT */

#ifdef FIX_OUTPUT_UNITS
  fix_units_for_ouput(output_galaxy_);
#endif /* defined FIX_OUTPUT_UNITS */
}


#ifdef FIX_OUTPUT_UNITS
/** @brief Removes h from units of galaxy properties. If desired (makefile option
 * FIX_OUTPUT_UNITS is set), the output properties of the galaxies can be scaled
 * to physical units excluding any factors of h == Hubble/100 km/s/Mpc. */
void fix_units_for_ouput(struct GALAXY_OUTPUT *output_galaxy_)
{
#ifdef LIGHT_OUTPUT
  output_galaxy_->Pos[0]            *= inv_Hubble_h;
  output_galaxy_->Pos[1]            *= inv_Hubble_h;
  output_galaxy_->Pos[2]            *= inv_Hubble_h;
  output_galaxy_->Mvir              *= inv_Hubble_h;
  output_galaxy_->Rvir              *= inv_Hubble_h;
  output_galaxy_->ColdGas           *= inv_Hubble_h;
  output_galaxy_->DiskMass          *= inv_Hubble_h;
  output_galaxy_->BulgeMass         *= inv_Hubble_h;
  output_galaxy_->HotGas            *= inv_Hubble_h;
  output_galaxy_->BlackHoleMass     *= inv_Hubble_h;

#else /* not defined LIGHT_OUTPUT */
  int j_;

  output_galaxy_->Pos[0]            *= inv_Hubble_h;
  output_galaxy_->Pos[1]            *= inv_Hubble_h;
  output_galaxy_->Pos[2]            *= inv_Hubble_h;
  output_galaxy_->CentralMvir       *= inv_Hubble_h;
  output_galaxy_->Mvir              *= inv_Hubble_h;
  output_galaxy_->Rvir              *= inv_Hubble_h;
#ifdef HALOPROPERTIES
  output_galaxy_->HaloM_Mean200     *= inv_Hubble_h;
  output_galaxy_->HaloM_Crit200     *= inv_Hubble_h;
  output_galaxy_->HaloM_TopHat      *= inv_Hubble_h;
  output_galaxy_->HaloPos[0]        *= inv_Hubble_h;
  output_galaxy_->HaloPos[1]        *= inv_Hubble_h;
  output_galaxy_->HaloPos[2]        *= inv_Hubble_h;
  output_galaxy_->HaloSpin[0]       *= inv_Hubble_h;
  output_galaxy_->HaloSpin[1]       *= inv_Hubble_h;
  output_galaxy_->HaloSpin[2]       *= inv_Hubble_h;
#endif /* defined HALOPROPERTIES */
  output_galaxy_->GasSpin[0]        *= inv_Hubble_h;
  output_galaxy_->GasSpin[1]        *= inv_Hubble_h;
  output_galaxy_->GasSpin[2]        *= inv_Hubble_h;
  output_galaxy_->StellarSpin[0]    *= inv_Hubble_h;
  output_galaxy_->StellarSpin[1]    *= inv_Hubble_h;
  output_galaxy_->StellarSpin[2]    *= inv_Hubble_h;
  output_galaxy_->HotRadius         *= inv_Hubble_h;
  output_galaxy_->ColdGas           *= inv_Hubble_h;
  output_galaxy_->DiskMass          *= inv_Hubble_h;
  output_galaxy_->BulgeMass         *= inv_Hubble_h;
  output_galaxy_->HotGas            *= inv_Hubble_h;
  output_galaxy_->EjectedMass       *= inv_Hubble_h;
  output_galaxy_->BlackHoleMass     *= inv_Hubble_h;
  output_galaxy_->ICM               *= inv_Hubble_h;
  output_galaxy_->BulgeSize         *= inv_Hubble_h;
  output_galaxy_->StellarDiskRadius *= inv_Hubble_h;
  output_galaxy_->GasDiskRadius     *= inv_Hubble_h;
  output_galaxy_->CoolingRadius     *= sqrt(inv_Hubble_h);
#ifdef TRACK_BURST
  output_galaxy_->BurstMass         *= inv_Hubble_h;
#endif /* defined TRACK_BURST */

  metals_multiply_by(&(output_galaxy_->MetalsColdGas    ), inv_Hubble_h);
  metals_multiply_by(&(output_galaxy_->MetalsDiskMass   ), inv_Hubble_h);
  metals_multiply_by(&(output_galaxy_->MetalsBulgeMass  ), inv_Hubble_h);
  metals_multiply_by(&(output_galaxy_->MetalsHotGas     ), inv_Hubble_h);
  metals_multiply_by(&(output_galaxy_->MetalsEjectedMass), inv_Hubble_h);
  metals_multiply_by(&(output_galaxy_->MetalsICM        ), inv_Hubble_h);

#ifdef STAR_FORMATION_HISTORY
  for(j_ = 0; j_ <= output_galaxy_->sfh_ibin; j_++)
  {
    output_galaxy_->sfh_DiskMass [j_] *= inv_Hubble_h;
    output_galaxy_->sfh_BulgeMass[j_] *= inv_Hubble_h;
    output_galaxy_->sfh_ICM      [j_] *= inv_Hubble_h;
#ifdef TRACK_BURST
    output_galaxy_->sfh_BurstMass[j_] *= inv_Hubble_h;
#endif /* defined TRACK_BURST */
    metals_multiply_by(&(output_galaxy_->sfh_MetalsDiskMass [j_]), inv_Hubble_h);
    metals_multiply_by(&(output_galaxy_->sfh_MetalsBulgeMass[j_]), inv_Hubble_h);
    metals_multiply_by(&(output_galaxy_->sfh_MetalsICM      [j_]), inv_Hubble_h);
  }
#endif /* defined STAR_FORMATION_HISTORY */
#endif /* not defined LIGHT_OUTPUT */
}
#endif /* defined FIX_OUTPUT_UNITS */


/** @brief calculated galaxy ID offset from file and tree no. */
long long calc_big_db_offset(const int file_number_, const int tree_number_)
{
  long long offset_;
#ifdef MRII
  offset_ = (((file_number_ * (long long) 1000000) + tree_number_) * (long long) 1000000000);
#else /* defined  */
  offset_ = (((file_number_ * (long long) 1000000) + tree_number_) * (long long) 1000000);
#endif /* defined  */
  return offset_;
}


/** @brief calculated subhalo ID offset from file and tree no. */
long long calc_big_db_subid_index(const int snapnum, const int file_number_, const int subhalo_index_)
{
  long long offset_;
#ifdef MRII
  offset_ = snapnum * (long long) 10000000000000 + file_number_ * (long long) 1000000000 + subhalo_index_;
#else /* not defined MRII */
  offset_ = snapnum * (long long) 1000000000000 + file_number_ * (long long) 100000000 + subhalo_index_;
#endif /* not defined MRII */
  return offset_;
}


#ifdef STAR_FORMATION_HISTORY
/** @brief writes SFH bins to separate file */
void write_sfh_bins()
{
  int snapshot_number_, j_;
  int n_bins_ = 0;
  for(snapshot_number_ = 0; snapshot_number_ < MAXSNAPS; snapshot_number_++)
    for(j_=0;j_ < SFH_ibin[snapshot_number_][0];j_++)
      n_bins_++;

  struct SFH_Time *sfh_times_ = (struct SFH_Time *) mymalloc("sfh_times_", sizeof(struct SFH_Time) * n_bins_);
  int i_bin_ = 0;
  for(snapshot_number_ = 0; snapshot_number_ < MAXSNAPS; snapshot_number_++)
  {
    for(j_=0;j_ < SFH_ibin[snapshot_number_][0];j_++)
    {
      sfh_times_[i_bin_].snapnum = snapshot_number_;
      sfh_times_[i_bin_].bin = j_;
      sfh_times_[i_bin_].lookbacktime = (SFH_t[snapshot_number_][0][j_]+SFH_dt[snapshot_number_][0][j_]/2.-NumToTime(snapshot_number_))*UnitTime_in_years * inv_Hubble_h;
      sfh_times_[i_bin_].dt=SFH_dt[snapshot_number_][0][j_]*UnitTime_in_years * inv_Hubble_h;
      sfh_times_[i_bin_].nbins = SFH_Nbins[snapshot_number_][0][j_];
      i_bin_++;
    }
  }

  FILE* SFH_Bins_File_;
  char file_name_[1000];
  sprintf(file_name_, "%s/SFH_Bins", FinalOutputDir);
  if(!(SFH_Bins_File_ = fopen(file_name_, "w")))
  {
    char error_message_[1000];
    sprintf(error_message_, "can't open file `%s'\n", file_name_);
    terminate(error_message_);
  }
  else
  { printf("writing sfh bins to %s\n", file_name_); }

  // write # bins
  myfwrite(&n_bins_, sizeof(int), 1, SFH_Bins_File_);        // write 1
  myfwrite(sfh_times_, sizeof(struct SFH_Time), n_bins_, SFH_Bins_File_);        // size of an output structure (Galaxy_Output)
  fflush(SFH_Bins_File_);
  fclose(SFH_Bins_File_);
}
#endif /* defined STAR_FORMATION_HISTORY */

