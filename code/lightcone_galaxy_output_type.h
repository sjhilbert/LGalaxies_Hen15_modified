/*  Copyright (C) <2016+>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in_ the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

/** @file lightcone_galaxy_output_type.h
 *  @brief   type for outputting galaxies on lightcone to disk
 *
 *  @author  Stefan Hilbert
 * 
 *  @date    2018
 */

#ifndef LIGHTCONE_GALAXY_OUTPUT_TYPE_H
#define LIGHTCONE_GALAXY_OUTPUT_TYPE_H

#include "allvars.h"


/** @brief struct for lightcone galaxies on disk
 * 
 * will be used as lightcone_galaxy_output_type,
 * if #defined LIGHTCONE_CUSTOM_OUTPUT
 * (otherwise, GALAXY_OUTPUT will be used instead)
 */
typedef struct lightcone_galaxy_custom_output_type_
{
#ifdef LIGHT_OUTPUT
  int       Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
  int       SnapNum; // The snapshot number where this galaxy was identified.
  float     CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float     CentralRvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float     Pos[3]; // 1/h Mpc - Galaxy Positions
  float     Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
  float     Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
  float     Vvir; // km/s - Virial velocity of the subhalo the galaxy is/was the center of.
  float     DistanceToCentralGal[3];

  /* baryonic reservoirs */
  float     ColdGas; // 10^10/h Msun - Mass in cold gas.
  float     BulgeMass; // 10^10/h Msun - Mass in the bulge
  float     DiskMass;
  float     HotGas; // 10^10/h Msun - Mass in hot gas
  float     BlackHoleMass; // 10^10/h Msun - Mass in black hole

  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_REST_MAGS
  float     MagDust[NMAG]; // dust corrected, rest-frame absolute mags
#endif /* defined OUTPUT_REST_MAGS */ 
#ifdef OUTPUT_OBS_MAGS
  float     ObsMagDust[NMAG]; // dust corrected, rest-frame absolute mags
#endif /* defined OUTPUT_OBS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

#else /* not defined LIGHT_OUTPUT */
#ifdef GALAXYTREE
  long long GalID; /** ID of galaxy, unique within simulation and SAM run.*/
  long long HaloID; // Unique ID of MPA halo containing this galaxy
  long long FOFCentralGal;
  long long SubID;
  long long MMSubID; // fofId, the subhaloid of the subhalo at the center of the fof group
  int       PeanoKey; // Peano-Hilbert key, (bits=8), for position in 500/h Mpc box
#endif /* defined GALAXYTREE */

  float     Redshift; // cosmological redshift of the galaxy (if GALAXYTREE then reuse Redshift defined there)
  float     ObsRedshift; // observed (l.galaxy_.s. cosmological + peculial velocity) redshift of the galaxy
  int       CubeShiftIndex; // index identifying periodic copy of simulation box the galaxy is in

  int       Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.

#ifdef HALOPROPERTIES
  float     HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
  float     HaloPos[3];
  float     HaloVel[3];
  float     HaloVelDisp;
  float     HaloVmax;
  float     HaloSpin[3];
#endif /* defined HALOPROPERTIES */
  
  int       SnapNum; // The snapshot number where this galaxy was identified.
  float     CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float     DistanceToCentralGal[3];
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float     Pos[3]; // 1/h Mpc - Galaxy Positions
  float     Vel[3]; // km/s - Galaxy Velocities
  int       Len; //Number of particles in the associated subhalo  
  float     Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
  float     Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
  float     Vvir; // km/s - Virial velocity of the subhalo the galaxy is/was the center of.
  float     Vmax; // km/s - Maximum rotational velocity of the subhalo, or the last value for type 2's galaxies.
  float     GasSpin[3]; // Gas Spin
  float     StellarSpin[3]; // Stellar Spin
  float     InfallVmax; // km/s - Vmax at infall
  float     InfallVmaxPeak; // km/s - Max previous Vmax at infall
  int       InfallSnap; // Snapnum at infall

  /* baryonic reservoirs */
  float     ColdGas; // 10^10/h Msun - Mass in cold gas.
  float     StellarMass; // 10^10/h Msun - Disk+Bulge
  float     BulgeMass; // 10^10/h Msun - Mass in the bulge
  float     DiskMass;
  float     HotGas; // 10^10/h Msun - Mass in hot gas
  float     EjectedMass; // 10^10/h Msun - Mass in ejected gas
  float     BlackHoleMass; // 10^10/h Msun - Mass in black hole
  float     ICM;            // mass in intra-cluster stars, for type 0,1

  /* misc */
  float     BulgeSize;
  float     StellarDiskRadius;
  float     GasDiskRadius;
  float     CosInclination; // cos(angle) between galaxy spin and the z-axis

  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_REST_MAGS
  float     MagDust[NMAG]; // dust corrected, rest-frame absolute mags
  float     Mag[NMAG]; // rest-frame absolute mags
  float     MagBulge[NMAG]; // rest-frame absolute mags for the bulge
#ifdef ICL
  float     MagICL[NMAG];          // rest-frame absolute mags of ICL
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
  float     ObsMagDust[NMAG]; // dust-corrected, obs-frame absolute mags
  float     ObsMag[NMAG]; // obs-frame absolute mags
  float     ObsMagBulge[NMAG]; // obs-frame absolute mags for the bulge
#ifdef ICL
  float     ObsMagICL[NMAG];  // observer-frame absolute mags for intra-cluster light
#endif /* defined ICL */
#endif /* defined OUTPUT_OBS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

  float     MassWeightAge;
  float     rbandWeightAge;
#endif /* not defined LIGHT_OUTPUT */
} lightcone_galaxy_custom_output_type;


#ifdef LIGHTCONE_CUSTOM_OUTPUT
typedef lightcone_galaxy_custom_output_type lightcone_galaxy_output_type;
#else /* not defined LIGHTCONE_CUSTOM_OUTPUT */
/* just use GALAXY_OUTPUT: */
typedef struct GALAXY_OUTPUT lightcone_galaxy_output_type;
#endif /* not defined LIGHTCONE_CUSTOM_OUTPUT */


/** @brief copies GALAXY_OUTPUT to lightcone_galaxy_output_type
*
*  @note  should only be needed (i.e. called) if 
*         defined LIGHTCONE_CUSTOM_OUTPUT.
*/
static inline void 
galaxy_output_to_lightcone_galaxy_output_type(const struct GALAXY_OUTPUT *galaxy_, lightcone_galaxy_output_type *lightcone_galaxy_)
{
#ifdef LIGHTCONE_CUSTOM_OUTPUT
  
  int i_;
  
#ifdef LIGHT_OUTPUT                           
  lightcone_galaxy_->Type                        = galaxy_->Type                       ;
  lightcone_galaxy_->SnapNum                     = galaxy_->SnapNum                    ;
  lightcone_galaxy_->CentralMvir                 = galaxy_->CentralMvir                ;
  lightcone_galaxy_->CentralRvir                 = galaxy_->CentralRvir                ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->Pos[i_]                     = galaxy_->Pos[i_]                    ;
  lightcone_galaxy_->Mvir                        = galaxy_->Mvir                       ;
  lightcone_galaxy_->Rvir                        = galaxy_->Rvir                       ;
  lightcone_galaxy_->Vvir                        = galaxy_->Vvir                       ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->DistanceToCentralGal[i_]    = galaxy_->DistanceToCentralGal[i_]   ;

  /* baryonic reservoirs */
  lightcone_galaxy_->ColdGas                     = galaxy_->ColdGas                    ;
  lightcone_galaxy_->BulgeMass                   = galaxy_->BulgeMass                  ;
  lightcone_galaxy_->DiskMass                    = galaxy_->DiskMass                   ;
  lightcone_galaxy_->HotGas                      = galaxy_->HotGas                     ;
  lightcone_galaxy_->BlackHoleMass               = galaxy_->BlackHoleMass              ;

  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_REST_MAGS
  for(i_ = 0; i_ < NMAG; i_++) 
  lightcone_galaxy_->MagDust[i_]                 = galaxy_->MagDust[i_]                ;
#endif /* defined OUTPUT_REST_MAGS */
#ifdef OUTPUT_OBS_MAGS
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->ObsMagDust[i_]              = galaxy_->ObsMagDust[i_]             ;
#endif /* defined OUTPUT_OBS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
#else /* not defined LIGHT_OUTPUT */

#ifdef GALAXYTREE
  lightcone_galaxy_->GalID                       = galaxy_->GalID                      ;
  lightcone_galaxy_->HaloID                      = galaxy_->HaloID                     ;
  lightcone_galaxy_->FOFCentralGal               = galaxy_->FOFCentralGal              ;
  lightcone_galaxy_->SubID                       = galaxy_->SubID                      ;
  lightcone_galaxy_->MMSubID                     = galaxy_->MMSubID                    ;
  lightcone_galaxy_->PeanoKey                    = galaxy_->PeanoKey                   ;
#endif /* defined GALAXYTREE */

  lightcone_galaxy_->Redshift                    = galaxy_->Redshift                   ;
  lightcone_galaxy_->ObsRedshift                 = galaxy_->ObsRedshift                ;
  lightcone_galaxy_->CubeShiftIndex              = galaxy_->CubeShiftIndex             ;
  lightcone_galaxy_->Type                        = galaxy_->Type                       ;

#ifdef HALOPROPERTIES
  lightcone_galaxy_->HaloM_Mean200               = galaxy_->HaloM_Mean200              ;
  lightcone_galaxy_->HaloM_Crit200               = galaxy_->HaloM_Crit200              ;
  lightcone_galaxy_->HaloM_TopHat                = galaxy_->HaloM_TopHat               ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->HaloPos[i_]                 = galaxy_->HaloPos[i_]                ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->HaloVel[i_]                 = galaxy_->HaloVel[i_]                ;
  lightcone_galaxy_->HaloVelDisp                 = galaxy_->HaloVelDisp                ;
  lightcone_galaxy_->HaloVmax                    = galaxy_->HaloVmax                   ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->HaloSpin[i_]                = galaxy_->HaloSpin[i_]               ;
#endif /* defined HALOPROPERTIES */

  lightcone_galaxy_->SnapNum                     = galaxy_->SnapNum                    ;
  lightcone_galaxy_->CentralMvir                 = galaxy_->CentralMvir                ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->DistanceToCentralGal[i_]    = galaxy_->DistanceToCentralGal[i_]   ;

  /* properties of subhalo at the last time...*/
  for(i_ = 0; i_ < 3; i_++)   
  lightcone_galaxy_->Pos[i_]                     = galaxy_->Pos[i_]                    ;
  for(i_ = 0; i_ < 3; i_++)            
  lightcone_galaxy_->Vel[i_]                     = galaxy_->Vel[i_]                    ;
  lightcone_galaxy_->Len                         = galaxy_->Len                        ;
  lightcone_galaxy_->Mvir                        = galaxy_->Mvir                       ;
  lightcone_galaxy_->Rvir                        = galaxy_->Rvir                       ;
  lightcone_galaxy_->Vvir                        = galaxy_->Vvir                       ;
  lightcone_galaxy_->Vmax                        = galaxy_->Vmax                       ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->GasSpin[i_]                 = galaxy_->GasSpin[i_]                ;
  for(i_ = 0; i_ < 3; i_++)
  lightcone_galaxy_->StellarSpin[i_]             = galaxy_->StellarSpin[i_]            ;
  lightcone_galaxy_->InfallVmax                  = galaxy_->InfallVmax                 ;
  lightcone_galaxy_->InfallVmaxPeak              = galaxy_->InfallVmaxPeak             ;
  lightcone_galaxy_->InfallSnap                  = galaxy_->InfallSnap                 ;

  /* baryonic reservoirs */
  lightcone_galaxy_->ColdGas                     = galaxy_->ColdGas                    ;
  lightcone_galaxy_->StellarMass                 = galaxy_->StellarMass                ;
  lightcone_galaxy_->BulgeMass                   = galaxy_->BulgeMass                  ;
  lightcone_galaxy_->DiskMass                    = galaxy_->DiskMass                   ;
  lightcone_galaxy_->HotGas                      = galaxy_->HotGas                     ;
  lightcone_galaxy_->EjectedMass                 = galaxy_->EjectedMass                ;
  lightcone_galaxy_->BlackHoleMass               = galaxy_->BlackHoleMass              ;
  lightcone_galaxy_->ICM                         = galaxy_->ICM                        ;

  /* misc */
  lightcone_galaxy_->BulgeSize                   = galaxy_->BulgeSize                  ;
  lightcone_galaxy_->StellarDiskRadius           = galaxy_->StellarDiskRadius          ;
  lightcone_galaxy_->GasDiskRadius               = galaxy_->GasDiskRadius              ;
  lightcone_galaxy_->CosInclination              = galaxy_->CosInclination             ;

  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_REST_MAGS           
  for(i_ = 0; i_ < NMAG; i_++)    
  lightcone_galaxy_->MagDust[i_]                 = galaxy_->MagDust[i_]                ;
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->Mag[i_]                     = galaxy_->Mag[i_]                    ;
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->MagBulge[i_]                = galaxy_->MagBulge[i_]               ;
#ifdef ICL
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->MagICL[i_]                  = galaxy_->MagICL[i_]                 ;
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->ObsMagDust[i_]              = galaxy_->ObsMagDust[i_]             ;
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->ObsMag[i_]                  = galaxy_->ObsMag[i_]                 ;
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->ObsMagBulge[i_]             = galaxy_->ObsMagBulge[i_]            ;
#ifdef ICL
  for(i_ = 0; i_ < NMAG; i_++)
  lightcone_galaxy_->ObsMagICL[i_]               = galaxy_->ObsMagICL[i_]              ;
#endif /* defined ICL */
#endif /* defined OUTPUT_OBS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

  lightcone_galaxy_->MassWeightAge               = galaxy_->MassWeightAge              ;
  lightcone_galaxy_->rbandWeightAge              = galaxy_->rbandWeightAge             ;
#endif /* not defined LIGHT_OUTPUT */
  
#else /* not defined LIGHTCONE_CUSTOM_OUTPUT */
  *lightcone_galaxy_ = *galaxy_; 
#endif /* not defined LIGHTCONE_CUSTOM_OUTPUT */
}


#ifdef GALAXYTREE
/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) for a galaxy in memory */
static inline void 
prepare_galaxy_tree_info_for_lightcone_output(const int file_number_, const int tree_number_, const struct galaxy_tree_data *tree_gal_, lightcone_galaxy_output_type *galaxy_)
{
  const long long big_db_offset_ = calc_big_db_offset(file_number_, tree_number_);

  galaxy_->GalID = tree_gal_->GalID;
  galaxy_->FOFCentralGal = tree_gal_->FOFCentralGal;
  
#ifndef LIGHTCONE_CUSTOM_OUTPUT
  galaxy_->FirstProgGal = tree_gal_->FirstProgGal;
  galaxy_->NextProgGal = tree_gal_->NextProgGal;
  galaxy_->LastProgGal = tree_gal_->LastProgGal;
  galaxy_->MainLeafId = tree_gal_->MainLeaf;
  galaxy_->TreeRootId = tree_gal_->TreeRoot;
  galaxy_->DescendantGal = tree_gal_->DescendantGal;
  galaxy_->FileTreeNr = big_db_offset_;
#endif /* not defined LIGHTCONE_CUSTOM_OUTPUT */

#ifdef CONTINUOUS_TREES
  // Reset big_offset_ (so only FileTreeNr has original value)
  // Then new values should coincide with positions in the file
  const long long big_offset_ = TotGalCount;
#else  /* not defined CONTINUOUS_TREES */
  const long long big_offset_ = big_db_offset_;
#endif /* not defined CONTINUOUS_TREES */

  galaxy_->GalID += big_offset_;
  galaxy_->FOFCentralGal += big_offset_;

#ifndef LIGHTCONE_CUSTOM_OUTPUT
  if(galaxy_->FirstProgGal >= 0)
    galaxy_->FirstProgGal += big_offset_;

  if(galaxy_->LastProgGal >= 0)
    galaxy_->LastProgGal += big_offset_;
  else
    galaxy_->LastProgGal = galaxy_->GalID;

  if(galaxy_->MainLeafId >= 0)
    galaxy_->MainLeafId += big_offset_;
  else
    galaxy_->MainLeafId = galaxy_->GalID;

  if(galaxy_->TreeRootId >= 0)
    galaxy_->TreeRootId += big_offset_;
  else
  {
    terminate("galaxy_->TreeRootId < 0");
    galaxy_->TreeRootId = -1;
  }

  if(galaxy_->NextProgGal >= 0)
    galaxy_->NextProgGal += big_offset_;

  if(galaxy_->DescendantGal >= 0)
    galaxy_->DescendantGal += big_offset_;
#endif /* not defined LIGHTCONE_CUSTOM_OUTPUT */
}
#endif /* defined GALAXYTREE */


/** @brief Reading routine for galaxies in_ lightcone files  */
static inline size_t 
myfread_lightcone_galaxy(lightcone_galaxy_output_type *galaxy_, const size_t n_to_read_, FILE * stream_)
{ return myfread(galaxy_, sizeof(lightcone_galaxy_output_type), n_to_read_, stream_); }


/** @brief writing routine for galaxies in_ lightcone files  */
static inline size_t 
myfwrite_lightcone_galaxy(lightcone_galaxy_output_type *galaxy_, const size_t n_to_write_, FILE * stream_)
{ return myfwrite(galaxy_, sizeof(lightcone_galaxy_output_type), n_to_write_, stream_); }


/** @brief writing routine for galaxies in_ lightcone files  */
static inline size_t 
myfwrite_lightcone_galaxy_from_galaxy_output(struct GALAXY_OUTPUT *galaxy_, const size_t n_to_write_, FILE * stream_)
{ 
#ifdef LIGHTCONE_CUSTOM_OUTPUT
  lightcone_galaxy_output_type lightcone_galaxy_; 
  size_t i_;
  for(i_ = 0; i_ < n_to_write_; ++i_)
  {
    galaxy_output_to_lightcone_galaxy_output_type(galaxy_ + i_, &lightcone_galaxy_);
    myfwrite(&lightcone_galaxy_, sizeof(lightcone_galaxy_output_type), 1, stream_);
  }
  return n_to_write_;
#else  /* not defined LIGHTCONE_CUSTOM_OUTPUT */
  return myfwrite(&galaxy_, sizeof(struct GALAXY_OUTPUT), 1, stream_);
#endif /* not defined LIGHTCONE_CUSTOM_OUTPUT */
}

#endif /* header guard */
