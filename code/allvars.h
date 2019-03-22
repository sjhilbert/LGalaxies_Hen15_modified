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

/** @file allvars.h declares all globally defined structs and variables 
  *                 and some constants.
  *                
  * Declares all globally defined structs and variables and some MACRO constants.
  *
  * Common mathematical constants (M_PI etc.) are provided by "mathematical_constants.h"
  * (unless already provided by the compliler).
  *
  * Phyical constants and units are provided by "physical_constants_and_units.h".
  */
#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdbool.h>

#include "mathematical_constants.h"
#include "physical_constants_and_units.h"


#define MIN_ALLOC_NUMBER       1000
#define ALLOC_INCREASE_FACTOR  1.1
#define ALLOC_DECREASE_FACTOR  0.7

//WATCH OUT! In the case of MCMC running both MR and MRII the larger value is used to "allocate" all the arrays
//inside the code its LastDarkMatterSnapShot+1 that defines the extent of the loops
//(in MCMC MR_plus_MRII mode this are not always identical)
#if defined MRII
#define  MAXSNAPS  68     /* Number of snapshots in the dark matter simulation */
#elif defined PHOENIX
#define  MAXSNAPS  72
#elif defined CATERPILLAR
#define  MAXSNAPS  256
#else  /* not defined any of MRII, PHOENIX, CATERPILLAR */
#define  MAXSNAPS  64  //NORMAL MILLENNIUM
#endif /* not defined any of MRII, PHOENIX, CATERPILLAR */


#define  MAXGALFAC 2.3 /*1.5/2.3 - maximum fraction of satellite without a halo (for memory allocation)*/

#define  STEPS 20		/* Number of integration intervals between two snapshots */

#define  ALLOCPARAMETER 50.  /* new definition !!! THIS HAS TO BE 50 !!! DONT EVER EVER EVER CHANGE !!! */


#ifdef GALAXYTREE
#ifdef NOUT
#undef NOUT
#endif /* defined NOUT */
#define NOUT MAXSNAPS
#endif /* defined GALAXYTREE */


#ifdef STAR_FORMATION_HISTORY

#define SFH_NMERGE 3  //  SFH_NMERGE=Nmax+1 (Nmax used in Shamshiri2014)

#ifdef CATERPILLAR
#define SFH_NBIN 24 //  CATERPILLAR - 256 snapshots
#else  /* not defined CATERPILLAR */
#define SFH_NBIN 20
#endif /* not defined CATERPILLAR */

#endif /* defined STAR_FORMATION_HISTORY */


/** @brief r-band filter number for computing rbandWeightAge
 *
 * @warning filter_number_==17 or 6 is not necessarily r-band,
 *          since filters are assigned from info in parameter file.
 *
 * @todo turn macro into program parameter read from file
 */
#ifndef R_BAND_FILTER_NUMBER
#ifdef MCMC
#define R_BAND_FILTER_NUMBER 6
#else  /* not defined MCMC */
#define R_BAND_FILTER_NUMBER 17
#endif /* not defined MCMC */
#endif /* not defined R_BAND_FILTER_NUMBER */

/* currently three levels of output buffering are recognized (0 = no, 1 = per tree, and 2 = per tree file): */
#ifdef OUTPUT_BUFFERING
#if OUTPUT_BUFFERING == 0
/* interpret as no buffering, so make equivalent to not defined OUTPUT_BUFFERING: */
#undef OUTPUT_BUFFERING
#elif OUTPUT_BUFFERING > 2
/* interpret as highest available level of buffering (currently two): */
#undef OUTPUT_BUFFERING
#define OUTPUT_BUFFERING 2
#endif /* OUTPUT_BUFFERING > 2 */
#endif /* defined OUTPUT_BUFFERING */


typedef enum GasComponentType_ { HotGasComponent, ColdGasComponent, EjectedGasComponent } GasComponentType;
extern const char* GasComponentStr[];


typedef enum StellarComponentType_ { DiskComponent, BulgeComponent, ICMComponent, BurstComponent } StellarComponentType;
extern const char* StellarComponentStr[];


#ifdef DETAILED_METALS_AND_MASS_RETURN
struct metals
{
  float type1a;
  float type2;
  float agb;
};

#ifdef INDIVIDUAL_ELEMENTS
/** @brief individual element content */
struct elements
{
  float H;
  float He;
#ifndef MAINELEMENTS
  float Cb; //NOTE: Carbon (C) is stored as Cb here
  float N;
#endif /* not defined MAINELEMENTS */
  float O;
#ifndef MAINELEMENTS
  float Ne;
#endif /* not defined MAINELEMENTS */
  float Mg;
#ifndef MAINELEMENTS
  float Si;
  float S;
  float Ca;
#endif /* not defined MAINELEMENTS */
  float Fe;
};


//Number of chemical elements tracked:
#ifdef MAINELEMENTS
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe]
#else  /* not defined MAINELEMENTS */
#define NUM_ELEMENTS 11 //All: [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe]
#endif /* not defined MAINELEMENTS */

#endif /* defined INDIVIDUAL_ELEMENTS */
#endif /* defined DETAILED_METALS_AND_MASS_RETURN */


/** @brief Galaxy structure for output */
#ifdef LIGHT_OUTPUT

struct GALAXY_OUTPUT
{
  int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
  int   SnapNum; // The snapshot number where this galaxy was identified.
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float CentralRvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float Pos[3]; // 1/h Mpc - Galaxy Positions
  float Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
  float Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
  float Vvir; // km/s -	Virial velocity of the subhalo the galaxy is/was the center of.
  float DistanceToCentralGal[3];

  /* baryonic reservoirs */
  float ColdGas; // 10^10/h Msun - Mass in cold gas.
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass;
  float HotGas; // 10^10/h Msun - Mass in hot gas
  float BlackHoleMass; // 10^10/h Msun - Mass in black hole
  
  float GasDiskRadius;
  float CosInclination;

  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_OBS_MAGS
  float ObsMagDust[NMAG]; // dust corrected, rest-frame absolute mags
#endif /* defined OUTPUT_OBS_MAGS */
#ifdef OUTPUT_REST_MAGS
  float MagDust[NMAG]; // dust corrected, rest-frame absolute mags
#endif /* defined OUTPUT_REST_MAGS */ 
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
};

#else /* not defined LIGHT_OUTPUT */

#ifdef PACK_OUTPUT
#pragma pack(1)  //structure alignment for 1 Byte(s).
#endif /* defined PACK_OUTPUT */
struct GALAXY_OUTPUT
{
#ifdef GALAXYTREE
  long long GalID; /** ID of galaxy, unique within simulation and SAM run.*/
  long long HaloID; // Unique ID of MPA halo containing this galaxy
#endif /* defined GALAXYTREE */
#ifdef MBPID
  long long MostBoundID; // Most bound particle at centre of subhalo last associated with this galaxy.  Put here as want all 8-byte blocks together at top of output record.
#endif /* defined MBPID */
#ifdef GALAXYTREE
  long long FirstProgGal;	// Main progenitor of this galaxy. Also the first progenitor in a linked list representation of the merger tree.
  long long NextProgGal;	// Next progenitor of this galaxy in linked list representation of merger tree
  long long LastProgGal;	// Galaxies with id between this galaxyId and this lastProgenitorId form the merger tree rooted in this galaxy.
  long long FOFCentralGal;
  long long FileTreeNr;
  long long DescendantGal;	// Pointer to the descendant of this galaxy in its merger tree; -1 if there is no descendant
  long long MainLeafId;
  long long TreeRootId;
  long long SubID;
  long long MMSubID; // fofId, the subhaloid of the subhalo at the center of the fof group
  int   PeanoKey; // Peano-Hilbert key, (bits=8), for position in 500/h Mpc box
#endif /* defined GALAXYTREE */
  float Redshift; // redshift of the snapshot where this galaxy resides
#ifdef LIGHTCONE_OUTPUT
  float ObsRedshift; // observed (l.o.s. cosmological + peculial velocity) redshift of the galaxy
  int   CubeShiftIndex; // index identifying periodic copy of simulation box the galaxy is in
#endif /* defined LIGHTCONE_OUTPUT */
  int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
#ifndef GALAXYTREE
  int   HaloIndex;
 // long long SubID;
 // long long FirstHaloInFOFgroup;
#endif /* not defined GALAXYTREE */
#ifdef HALOPROPERTIES
  float HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
  float HaloPos[3];
  float HaloVel[3];
  float HaloVelDisp;
  float HaloVmax;
  float HaloSpin[3];
#endif /* defined HALOPROPERTIES */
  int   SnapNum; // The snapshot number where this galaxy was identified.
  float LookBackTimeToSnap; //The time from a given snapshot to z=0, in years
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float CentralRvir; // Rvir of background (FOF) halo containing this galaxy
  float DistanceToCentralGal[3];
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float Pos[3]; // 1/h Mpc - Galaxy Positions
  float Vel[3]; // km/s - Galaxy Velocities
  int   Len; //Number of particles in the associated subhalo  
  float Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
  float Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
  float Vvir; // km/s -	Virial velocity of the subhalo the galaxy is/was the center of.
  float Vmax; // km/s - Maximum rotational velocity of the subhalo, or the last value for type 2's galaxies.
  float GasSpin[3]; // Gas Spin
  float StellarSpin[3]; // Stellar Spin
  float InfallVmax; // km/s - Vmax at infall
  float InfallVmaxPeak; // km/s - Max previous Vmax at infall
  int InfallSnap; // Snapnum at infall
  float InfallHotGas;
  float HotRadius; //Mpc/h - Radius of the hot gas
  /*dynamical friction merger time*/
  float OriMergTime;
  float MergTime;
  /* baryonic reservoirs */
  float ColdGas; // 10^10/h Msun - Mass in cold gas.
  float StellarMass; // 10^10/h Msun - Disk+Bulge
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass;
  float HotGas; // 10^10/h Msun - Mass in hot gas
  float EjectedMass; // 10^10/h Msun - Mass in ejected gas
  float BlackHoleMass; // 10^10/h Msun - Mass in black hole
  /* ICL magnitude and mass*/
  float ICM;            // mass in intra-cluster stars, for type 0,1
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals MetalsColdGas; // 10^10/h Msun -	Mass in metals in cold gas.
  struct metals MetalsBulgeMass; // 10^10/h Msun -	Mass in metals in the bulge
  struct metals MetalsDiskMass; // 10^10/h Msun -       Mass in metals in the disk
  struct metals MetalsHotGas; // 10^10/h Msun -	Mass in metals in the hot gas
  struct metals MetalsEjectedMass; // 10^10/h Msun -	Mass in metals in the ejected gas
  struct metals MetalsICM;  // total mass in metals in intra-cluster stars, for type 0,1
#ifdef METALS_SELF
  struct metals MetalsHotGasSelf; // hot gas metals that come from self
#endif /* defined METALS_SELF */
#else /* not defined DETAILED_METALS_AND_MASS_RETURN */
  float MetalsColdGas; // 10^10/h Msun -	Mass in metals in cold gas.
  float MetalsStellarMass; // 10^10/h Msun -	Mass in metals in the bulge+disk
  float MetalsBulgeMass; // 10^10/h Msun -	Mass in metals in the bulge
  float MetalsDiskMass; // 10^10/h Msun -       Mass in metals in the disk
  float MetalsHotGas; // 10^10/h Msun -	Mass in metals in the hot gas
  float MetalsEjectedMass; // 10^10/h Msun -	Mass in metals in the ejected gas
  float MetalsICM;  // total mass in metals in intra-cluster stars, for type 0,1
#ifdef METALS_SELF
  float MetalsHotGasSelf; // hot gas metals that come from self
#endif /* defined METALS_SELF */
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */
#ifdef TRACK_BURST
  float BurstMass; // Mass formed in starbursts
#endif /* defined TRACK_BURST */
  /* misc */
  float PrimordialAccretionRate;
  float CoolingRadius;  // Q: store this ? (was stored in Delucia20006a)
  float CoolingRate;
  float CoolingRate_beforeAGN;
  float QuasarAccretionRate;
  float RadioAccretionRate;
  float Sfr;
  float SfrBulge;
  float XrayLum;
  float BulgeSize;
  float StellarDiskRadius;
  float GasDiskRadius;
  float CosInclination; // cos(angle) between galaxy spin and the z-axis
  int   DisruptOn; // 0: galaxy merged onto merger center; 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center
  int   MergeOn;   // 0: standard delucia-like merger behaviour for type 1 galaxy; 1: galaxy mass > halo mass, separate dynamical friction time calculated ....
  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_REST_MAGS
  float MagDust[NMAG]; // dust corrected, rest-frame absolute mags
  float Mag[NMAG]; // rest-frame absolute mags
  float MagBulge[NMAG]; // rest-frame absolute mags for the bulge
#ifdef ICL
  float MagICL[NMAG];          // rest-frame absolute mags of ICL
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
  float ObsMagDust[NMAG]; // dust-corrected, obs-frame absolute mags
  float ObsMag[NMAG]; // obs-frame absolute mags
  float ObsMagBulge[NMAG]; // obs-frame absolute mags for the bulge
#ifdef ICL
  float ObsMagICL[NMAG];  // observer-frame absolute mags for intra-cluster light
#endif /* defined ICL */
#ifdef OUTPUT_FB_OBS_MAGS
  // define luminosities as if the galaxy were one snapshot earlier, i.e. higher redshift, than its actual snapshot
  float backward_ObsMagDust [NMAG];
  float backward_ObsMag     [NMAG];
  float backward_ObsMagBulge[NMAG];
#ifdef ICL
  float backward_ObsMagICL  [NMAG];
#endif /* defined ICL */
  // define luminosities as if the galaxy were one snapshot later, i.e. lower redshift, than its actual snapshot
  float forward_ObsMagDust  [NMAG];
  float forward_ObsMag      [NMAG];
  float forward_ObsMagBulge [NMAG];
#ifdef ICL
  float forward_ObsMagICL   [NMAG];
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */

#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

  float MassWeightAge;
  float rbandWeightAge;

#ifdef STAR_FORMATION_HISTORY
  int sfh_ibin; //Index of highest bin currently in use
  int sfh_numbins; // number of non empty bins
  float sfh_DiskMass[SFH_NBIN];
  float sfh_BulgeMass[SFH_NBIN];
  float sfh_ICM[SFH_NBIN];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass[SFH_NBIN]; // Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
#else /* not defined DETAILED_METALS_AND_MASS_RETURN */
  float sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  float sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */
#ifdef TRACK_BURST
  float sfh_BurstMass[SFH_NBIN]; // Mass formed in starbursts
#endif /* defined TRACK_BURST */
#endif /* defined STAR_FORMATION_HISTORY */

#ifdef INDIVIDUAL_ELEMENTS
  struct elements sfh_ElementsDiskMass[SFH_NBIN];
  struct elements sfh_ElementsBulgeMass[SFH_NBIN];
  struct elements sfh_ElementsICM[SFH_NBIN];

  //float DiskMass_elements[ELEMENT_NUM];
  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
#endif /* defined INDIVIDUAL_ELEMENTS */
};

// next only of interest to DB output, which generally requires complete tree
#ifdef STAR_FORMATION_HISTORY
struct SFH_BIN 
{
	long long GalID; // ID of the galaxy
	short snapnum; // snapnum of the galaxy, repeated here for faster lookups of times etc
    short sfh_ibin; //Index of highest bin currently in use
//    float sfh_time; //time to present at the middle of bin in years.
//    float sfh_dt; //time width of bin in years.
  float sfh_DiskMass;
  float sfh_BulgeMass;
  float sfh_ICM;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_MetalsDiskMass; // Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass; // Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM; // Metals locked up in stars in ICM.
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
  float sfh_MetalsDiskMass; // Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass; //Metals locked up in stars in bulge.
  float sfh_MetalsICM; // Metals locked up in stars in ICM.
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */
#ifdef TRACK_BURST
  float sfh_BurstMass; // Mass formed in starbursts
#endif /* defined TRACK_BURST */

#ifdef INDIVIDUAL_ELEMENTS
  struct elements sfh_ElementsDiskMass;
  struct elements sfh_ElementsBulgeMass;
  struct elements sfh_ElementsICM;

  //float DiskMass_elements[ELEMENT_NUM];
  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
#endif  /* defined INDIVIDUAL_ELEMENTS */
};


struct SFH_Time
{
 	int snapnum; // snapnum
 	int bin; // index of current bin
 	double lookbacktime; // lookback time in years (???) to center of current bin
 						// proposal: in output write the start of the bin and its end, rather than center and dt
 	double dt; // width of the current bin in years (???)
 	int nbins; // # of highest resolution bins used to create current bin
};
#endif  /* defined STAR_FORMATION_HISTORY */
#pragma pack()   //structure alignment ends.
#endif /* not defined LIGHT_OUTPUT */


extern struct galaxy_tree_data
{
  int HaloGalIndex;
  int IndexStored;
  int SnapNum;
  int GalID;
  int FirstProgGal;
  int NextProgGal;
  int LastProgGal;
  int DescendantGal;
  int MainLeaf;
  int TreeRoot;
  int FOFCentralGal;
  int Done;
#ifdef LIGHTCONE_OUTPUT
  long long lightcone_galaxy_number_in_file_begin;
  long long lightcone_galaxy_number_in_file_end;
#endif /* defined LIGHTCONE_OUTPUT */
}
 *GalTree;
 
 
/** @brief Structure with all the data associated with galaxies (this is not the same as the output!) */
extern struct GALAXY
{
  int HeapIndex;
  int GalTreeIndex;
  int NextGalaxy;
#ifdef GALAXYTREE
  int FirstProgGal;
#endif /* defined GALAXYTREE */
  int Type;
  int HaloNr;
  long long MostBoundID;
  int SnapNum;
  int CentralGal;  //own ID for types 0 and 1, unless 1's have little dark matter and already merging to type 0. For 2's its the merger centre
  float CentralMvir;
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float Pos[3];
  float MergCentralPos[3];
  float Vel[3];
  float Pos_notupdated[3];
  float Vel_notupdated[3];
#ifdef HALOPROPERTIES
  float HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
  float HaloPos[3];
  float HaloVel[3];
  float HaloVelDisp;
  float HaloVmax;
#endif /* defined HALOPROPERTIES */
  float HaloSpin[3];
  float GasSpin[3];
  float StellarSpin[3];
  int   Len;   
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float InfallVmax;
  float InfallVmaxPeak; // km/s - Max previous Vmax at infall
  int InfallSnap;
  float InfallHotGas;
  float InfallHotGasRadius;
  float HotRadius;
  /* baryonic reservoirs */
  float ColdGas;
  float BulgeMass;
  float DiskMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals MetalsColdGas;
  struct metals MetalsBulgeMass;
  struct metals MetalsDiskMass;
  struct metals MetalsHotGas;
  struct metals MetalsEjectedMass;
#ifdef METALS_SELF
  struct metals MetalsHotGasSelf;
#endif /* defined METALS_SELF */ 
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */ 
  float MetalsColdGas;
  float MetalsBulgeMass;
  float MetalsDiskMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
#ifdef METALS_SELF
  float MetalsHotGasSelf;
#endif /* defined METALS_SELF */ 
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */ 
#ifdef TRACK_BURST
  float BurstMass;
#endif /* defined TRACK_BURST */ 

  /* misc */
  float PrimordialAccretionRate;
  float CoolingRate;
  float CoolingRate_beforeAGN;
  float CoolingRadius;
  float CoolingGas;
  float QuasarAccretionRate;
  float RadioAccretionRate;
  float AGNheatingFromCentral;
  float Sfr;
  float SfrBulge;
  float StarMerge;
  float XrayLum;
  float BulgeSize;
  float StellarDiskRadius;
  float GasDiskRadius;
#ifdef GALAXYTREE
  int   DisruptOn;
#endif /* defined GALAXYTREE */
  // float halfradius;
  //float periradius;
  float CosInclination; //angle between galaxy spin and the z-axis
  float OriMergTime;
  float MergTime;
  float OriMvir;
  float OriRvir;
  float MergeSat;
  float DistanceToCentralGal[3];
  int MergeOn;
  float ICM;
#ifdef DETAILED_METALS_AND_MASS_RETURN
   struct metals MetalsICM;
 #else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
   float MetalsICM;
 #endif /* not defined DETAILED_METALS_AND_MASS_RETURN */

  /* luminosities in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
  float Lum      [NOUT][NMAG];
  float LumY     [NOUT][NMAG];
  float LumBulge [NOUT][NMAG];
  float LumBulgeY[NOUT][NMAG];
  float LumDust  [NOUT][NMAG];
#ifdef ICL             
  float LumICL   [NOUT][NMAG];
#endif /* defined ICL */
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
  float ObsLum      [NOUT][NMAG];
  float ObsLumY     [NOUT][NMAG];
  float ObsLumBulge [NOUT][NMAG];
  float ObsLumBulgeY[NOUT][NMAG];
  float ObsLumDust  [NOUT][NMAG];
#ifdef ICL                
  float ObsLumICL   [NOUT][NMAG];
#endif /* defined ICL */

#ifdef OUTPUT_FB_OBS_MAGS
  float backward_ObsLum      [NOUT][NMAG];
  float backward_ObsLumY     [NOUT][NMAG];
  float backward_ObsLumBulge [NOUT][NMAG];
  float backward_ObsLumBulgeY[NOUT][NMAG];
  float backward_ObsLumDust  [NOUT][NMAG];
#ifdef ICL                 
  float backward_ObsLumICL  [NOUT][NMAG];
#endif /* defined ICL */

  float forward_ObsLum      [NOUT][NMAG];
  float forward_ObsLumY     [NOUT][NMAG];
  float forward_ObsLumBulge [NOUT][NMAG];
  float forward_ObsLumBulgeY[NOUT][NMAG];
  float forward_ObsLumDust  [NOUT][NMAG];
#ifdef ICL                         
  float forward_ObsLumICL   [NOUT][NMAG];
#endif /* defined ICL */
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */

#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

  float MassWeightAge[NOUT];
  float rbandWeightAge[NOUT];
 
#ifdef STAR_FORMATION_HISTORY
  int sfh_ibin; //Index of highest bin currently in use
  double sfh_age; //Time in years of last call to sph_update_bins
  int sfh_flag[SFH_NBIN];
  float sfh_dt[SFH_NBIN]; //Size of time interval in units of years
  float sfh_t[SFH_NBIN]; //Time at low-redshift edge of bin in same units
  int sfh_Nbins[SFH_NBIN]; //Number of bins on the time interval
  float sfh_DiskMass[SFH_NBIN]; //Stellar mass in disk, in bin in standard units
  float sfh_BulgeMass[SFH_NBIN]; //Stellar mass in bulge, in bin in standard units
  float sfh_ICM[SFH_NBIN]; //Stellar mass in ICM, in bin in standard units
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_MetalsDiskMass[SFH_NBIN]; //Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM[SFH_NBIN]; //Metals locked up in stars in ICM.
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
  float sfh_MetalsDiskMass[SFH_NBIN]; //Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  float sfh_MetalsICM[SFH_NBIN]; //Metals locked up in stars in ICM.
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */
#ifdef TRACK_BURST
  float sfh_BurstMass[SFH_NBIN]; //Stellar mass formed in bursts, in standard units.
#endif /* defined TRACK_BURST */
#endif /* defined STAR_FORMATION_HISTORY */

#ifdef INDIVIDUAL_ELEMENTS
  struct elements sfh_ElementsDiskMass[SFH_NBIN];
  struct elements sfh_ElementsBulgeMass[SFH_NBIN];
  struct elements sfh_ElementsICM[SFH_NBIN];

  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
#endif /* defined INDIVIDUAL_ELEMENTS */
} *Gal, *HaloGal;


// Documentation can be found in the database
extern struct halo_data
{
	/* merger tree pointers */
	int Descendant;
	int FirstProgenitor;
	int NextProgenitor;
	int FirstHaloInFOFgroup;
	int NextHaloInFOFgroup;

  /* properties of halo */
	int Len;
	float M_Mean200, M_Crit200, M_TopHat;
	float Pos[3];
	float Vel[3];
	float VelDisp;
	float Vmax;
	float Spin[3];
	long long MostBoundID;

  /* original position in subfind output */
	int SnapNum;
	int FileNr;
	int SubhaloIndex;
	float SubHalfMass;
}
  *Halo, *Halo_Data;


// Documentation can be found in the database
#ifndef MCMC
extern struct halo_ids_data
{
  long long HaloID;
  long long FileTreeNr;
  long long FirstProgenitor;
  long long LastProgenitor;
  long long NextProgenitor;
  long long Descendant;
  long long FirstHaloInFOFgroup;
  long long NextHaloInFOFgroup;
#ifdef MRII
  long long MainLeafID; 
#endif /* defined MRII */
  double    Redshift;
  int       PeanoKey;
  int       dummy;      /* need to use this padding for 64bit alignment */
} *HaloIDs, *HaloIDs_Data;
#else  /* defined MCMC */
extern struct  halo_ids_data
{
  long long FirstHaloInFOFgroup;
  int       MCMC_FOF_number[NOUT];
} *HaloIDs, *HaloIDs_Data;
#endif /* defined MCMC */


// Documentation can be found in the database
extern struct halo_aux_data  /* auxiliary halo data */
{
	int DoneFlag;
	int HaloFlag;
	int NGalaxies;
	int FirstGalaxy;
#ifdef ALLOW_UNSCALE_COSMOLOGY
	float M_Mean200_Unscaled;
	float M_Crit200_Unscaled;
	float Pos_Unscaled[3];
	float Vel_Unscaled[3];
	float Vmax_Unscaled;
	float Spin_Unscaled[3];
#endif /* defined ALLOW_UNSCALE_COSMOLOGY */
}
 *HaloAux;

extern time_t GlobalStartingTime;

extern int FirstFile;		/* first and last file for processing */
extern int LastFile;

extern int Ntrees;		/* number of trees in current file */
extern double AllocValue_MaxHaloGal;
extern double AllocValue_MaxGal;
extern double AllocValue_MaxGalTree;

extern int MaxGal;		/* Maximum number of galaxies allowed for Gal[] array */
extern int NHaloGal, MaxHaloGal;
extern int NGalTree, MaxGalTree;
extern int *HaloGalHeap;
extern int IndexStored;

extern int LastSnapShotNr;

extern int LastDarkMatterSnapShot;
#ifdef MR_PLUS_MRII //OPTION for MCMC
extern int LastDarkMatterSnapShot_MR;
extern int LastDarkMatterSnapShot_MRII;
#endif /* defined MR_PLUS_MRII */


extern char SpecPhotDir[512];
extern char PhotPrefix[50];
extern char SpecPhotIMF[50];
extern char McFile[512];
extern char FileWithFilterNames[512];
extern char CoolFunctionsDir[512];
extern char OutputDir[512];
/* in case a second parameter is given as argument to the code, this will be taken as a
 * temporary outputdir to allow fast I/O. OutputDir will be replaced by this directory
 * and in the end everything will be moved to the FinalOutputDir (original OutputDir
 * given in input.par )*/
extern char FinalOutputDir[512];
extern char FileNameGalaxies[512];
extern char SimulationDir[512];
extern char FileWithOutputRedshifts[512];

extern char FileWithZList[512];
//variable used to scale to a different cosmology
extern char FileWithZList_OriginalCosm[512];
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern char FileWithZList_MR[512];
extern char FileWithZList_OriginalCosm_MR[512];
extern char FileWithZList_MRII[512];
extern char FileWithZList_OriginalCosm_MRII[512];
#endif /* defined MR_PLUS_MRII */

extern double ScalePos;
extern double ScaleMass;

#ifdef SPECIFYFILENR
extern char   FileNrDir[512];
extern int    ListInputFileNr[111];
#endif /* defined SPECIFYFILENR */

extern int TotHalos;
extern int TotGalaxies[NOUT];
extern int *TreeNgals[NOUT];

extern int *FirstHaloInSnap;

extern int *TreeNHalos;
extern int *TreeFirstHalo;

extern void *TreeAuxData;


extern double MaxMemSize;

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;

#ifdef PARALLEL
extern int ThisTask, NTask;
#else /* not defined PARALLEL */
#define NTask 1
#define ThisTask 0
#endif /* not defined PARALLEL */

extern int TotGalCount;

/* Cosmological parameters */
extern double BaryonFrac;
extern double Sigma8;
extern double Omega;
extern double OmegaLambda;
extern double Hubble_h;
extern double inv_Hubble_h;
extern double Omega_OriginalCosm;
extern double OmegaLambda_OriginalCosm;
extern double Hubble_h_OriginalCosm;
//SIMULATION RELATED
extern double PartMass;
extern double BoxSize;
extern double PartMass_OriginalCosm;
extern double BoxSize_OriginalCosm;
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern double PartMass_MR;
extern double BoxSize_MR;
extern double PartMass_OriginalCosm_MR;
extern double BoxSize_OriginalCosm_MR;
extern double PartMass_MRII;
extern double BoxSize_MRII;
extern double PartMass_OriginalCosm_MRII;
extern double BoxSize_OriginalCosm_MRII;
#endif /* defined MR_PLUS_MRII */


/* flags */
extern int ReionizationModel;
extern int DiskRadiusModel;
extern int StarFormationModel;
extern int FeedbackReheatingModel;
extern int FeedbackEjectionModel;
extern int FateOfSatellitesGas;
extern int ReIncorporationModel;
extern int AGNRadioModeModel;
extern int DiskInstabilityModel;
extern int BHGrowthInDiskInstabilityModel;
extern int HotGasStrippingModel;
extern int DisruptionModel;
extern int StarBurstModel;
extern int BulgeFormationInMinorMergersOn;
extern int MetallicityOption;

/* parameters */
extern double Reionization_z0;
extern double Reionization_zr;
extern double Yield;
extern double RecycleFraction;
extern double ThreshMajorMerger;
extern double MergerTimeMultiplier;
extern double RamPressureStrip_CutOffMass;
extern double SfrEfficiency;
extern double SfrColdCrit;
extern double SfrBurstEfficiency;
extern double SfrBurstSlope;
extern double AgnEfficiency;
extern double BlackHoleGrowthRate;
extern double BlackHoleSeedMass;
extern double BlackHoleCutoffVelocity;
extern double FeedbackReheatingEpsilon;
extern double ReheatPreVelocity;
extern double ReheatSlope;
extern double FeedbackEjectionEfficiency;
extern double EjectPreVelocity;
extern double EjectSlope;
extern double ReIncorporationFactor;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;

// internal units and phys. consts. in internal units moved to "physical_constants_and_units.h"
// extern double UnitTime_in_s;
// extern double UnitPressure_in_cgs;
// extern double UnitDensity_in_cgs;
// extern double UnitCoolingRate_in_cgs;
// extern double UnitEnergy_in_cgs;
// extern double UnitTime_in_Megayears; //Using time as stored in the code, this gives Myr/h
// extern double UnitTime_in_years;
// 
// extern double Gravity;
// extern double SpeedOfLight;
// extern double Hubble;
// extern double RhoCrit;
  
extern double	a0, ar;

extern int ListOutputSnaps[NOUT];
extern float ListOutputRedshifts[NOUT];

extern int ListOutputNumberOfSnapshot[MAXSNAPS];

extern double ZZ[MAXSNAPS];
extern double AA[MAXSNAPS];
//variable used to scale to a different cosmology
extern double AA_OriginalCosm[MAXSNAPS];

extern double Age[MAXSNAPS];

extern int    Zlistlen;

extern gsl_rng *random_generator;


extern int    NumMergers;


/*  tabulated stuff */

#ifdef STAR_FORMATION_HISTORY
/* SFH_ is the reference structure for storing the star formation histories in
 * logarithmic bins. It is computed in init.c generating a binning structure for
 * each snapshot/time step. In the code galaxy structures are adjusted with respect
 * to this structure at each step. */
extern double SFH_t[MAXSNAPS][STEPS][SFH_NBIN]; //Time to present at the lower edge of the bin (code units)
extern double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN]; //Time width of the bin (code units)
extern int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN]; //Number of bins merged in each bin (only useful for the merging algorithm)
extern int SFH_ibin[MAXSNAPS][STEPS]; //Last active bin
#ifdef DETAILED_METALS_AND_MASS_RETURN
extern double tau_t[STEPS*MAXSNAPS]; //Time-to-z=0 of every timestep in the code. (Used for SNe rates in yield_integrals.c)
extern double tau_dt[STEPS*MAXSNAPS];//Width of every timestep in the code. (Used for SNe rates in yield_integrals.c)
#endif /* defined DETAILED_METALS_AND_MASS_RETURN */
#endif /* defined STAR_FORMATION_HISTORY */


#ifdef DETAILED_METALS_AND_MASS_RETURN

//Number of interpolated points within the mass ranges for the four types of yield table:
#define LIFETIME_MASS_NUM 150
#define LIFETIME_Z_NUM 6
#define AGB_MASS_NUM 59 //55 //ROB: 59, when going from 0.85 to 7 Msun
#define AGB_Z_NUM 3
#ifdef PORTINARI
#define SNII_MASS_NUM 85  //ROB: 85, from 6 <= M[Msun] <= 120. Change SNII_MIN_MASS and SNII_MAX_MASS for shorter ranges.
#define SNII_Z_NUM 5
#endif /* defined PORTINARI */
#ifdef CHIEFFI
#define SNII_MASS_NUM 81 //ROB: 56 if 7 <= M[Msun] <= 50. 81 if 7 <= M[Msun] <= 120. (NB: You can set SNII_MASS_NUM 81, and SNII_MAX_MASS 50. But DON"T put SNII_MASS_NUM > 81 ever!)
#define SNII_Z_NUM 6
#endif /* defined CHIEFFI */
#define SNIA_MASS_NUM 83 //48 //Number increased after extending range to cover M2 masses (07-02-12)

//Mass ranges for the different modes of ejection:
#define AGB_MIN_MASS 0.85
#define AGB_MAX_MASS 7.0 //6.0
#define SNIA_MIN_MASS 3.0
#define SNIA_MAX_MASS 16.0
#ifdef PORTINARI
#define SNII_MIN_MASS 7.0 //6.0
#define SNII_MAX_MASS 120.0
#endif /* defined PORTINARI */
#ifdef CHIEFFI
#define SNII_MIN_MASS 7.0
#define SNII_MAX_MASS 120.0 //50.0
#endif /* defined CHIEFFI */

int ELETOBIGCOUNTA;
int FRACCOUNTA;

//Arrays that yield tables are written to:
float lifetimeMasses[LIFETIME_MASS_NUM];
float lifetimeMetallicities[LIFETIME_Z_NUM];
float lifetimes[LIFETIME_Z_NUM][LIFETIME_MASS_NUM];
float AGBMasses[AGB_MASS_NUM]; //Initial star masses [Msun]
float AGBMetallicities[AGB_Z_NUM]; //Initial star metallicities [Msun]
float AGBEjectedMasses[AGB_Z_NUM][AGB_MASS_NUM]; //Total mass ejected [Msun]
float AGBTotalMetals[AGB_Z_NUM][AGB_MASS_NUM]; //Total metal YIELD ejected [Msun]
float AGBYields[AGB_Z_NUM][11][AGB_MASS_NUM]; //YIELD ejected, for each element [Msun]
float SNIIMasses[SNII_MASS_NUM];
float SNIIMetallicities[SNII_Z_NUM];
float SNIIEjectedMasses[SNII_Z_NUM][SNII_MASS_NUM];
float SNIITotalMetals[SNII_Z_NUM][SNII_MASS_NUM];
float SNIIYields[SNII_Z_NUM][11][SNII_MASS_NUM];
#ifndef DTD
float SNIaMasses[SNIA_MASS_NUM];
float SNIaEjectedMasses[SNIA_MASS_NUM];
float SNIaTotalMetals[SNIA_MASS_NUM];
float SNIaYields[42][SNIA_MASS_NUM];
#else
float SNIaYields[42];
#endif

//Integrated yields arrays:
float NormSNIIMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIIMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
float NormSNIIYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif
float NormAGBMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormAGBMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
float NormAGBYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif
float NormSNIaMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIaMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
float NormSNIaYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif

//Arrays used to plot SNe rates from SFH bins (yield_integrals.c):
float TheSFH[SFH_NBIN];
float SNIIRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float SNIaRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float AGBRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
//Arrays used to plot SNe rates from SFH-timesteps (calc_SNe_rates.c):
float TheSFH2[STEPS*MAXSNAPS];
float SNIIRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float SNIaRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float AGBRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];

//SNIa parameters:
#define A_FACTOR 0.04 //Fraction of mass from all objects between SNIa_MIN_MASS and SNIA_MAX_MASS that comes from SN-Ia.
//#define FRAC2HOT 0.9 //Fraction of material released by disk stars that goes straight into the HotGas. Res goes in ColdGas.
#ifdef DTD
//#define KALPHA 1.4765 //1.59203 //Now set in yield_integrals.c
//#define	F316 0.0384 //Integral of the IMF (by number) from 3.0 - 16.0 Msun //Now set in yield_integrals.c
#define SNIAEJECMASS 1.2300971 //Total mass (and total metals) ejected by a SNIa explosion in Msun //Value form original yield table (42 elements): 1.3740855. //Value when only considering 11 elements: 1.2300971
#ifdef BIMODALDTD
	#define DTD_NORM 0.903206 //For P98 Z=0.02 lifetimes (26Myrs - 21Gyrs)
#endif
#ifdef CUSTOMDTD
	#define DTD_NORM 0.524836 //For P98 Z=0.02 lifetimes (26Myrs - 21Gyrs)
#endif
#ifdef GAUSSIANDTD
	#define DTD_NORM = 1.0
	#define TAUCHARAC 1.0 //Characteristic delay time for SNe-Ia (i.e. peak of Gaussian distribution) in Gyrs //default: 2.0
	#define SIGMA_TD 0.2*TAUCHARAC //0.2 for narrow-DTD, 0.5 for wide_DTD
#endif
#ifdef POWERLAWDTD
	#define DTD_NORM 7.21863 //For P98 Z=0.02 lifetimes (26Myrs - 21Gyrs)
	#define DTD_SLOPE -1.12 //Slope of power law, according to Maoz et al. (2012)
#endif
#ifdef RUITERDTD
	#define DTD_NORM 1.09545 //For P98 Z=0.02 lifetimes (26Myrs - 21Gyrs)
	#define TAUCHARAC 0.5 //Peak of Gaussian (prompt) component [in Gyrs]
	#define SIGMA_TD 0.2*TAUCHARAC //Width of Gaussian (prompt) component
	#define DTD_SLOPE -2.0 //Slope of power law (delayed) component (see Ruiter et al. 2012)
#endif
#endif

#endif //DETAILED_METALS_AND_MASS_RETURN


#ifdef COMPUTE_SPECPHOT_PROPERTIES
// SSP PHOT_TABLES - magnitues of starburst population as a function of age

#ifdef M05
#define SSP_NAGES 220		// Age grid of the SSP tables
#define SSP_NMETALLICITES 4			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef BC03
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef CB07
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1238
#endif
#endif

//table containing the Metallicity grid of the SSP tables (converted to log10)
extern float SSP_logMetalTab[SSP_NMETALLICITES];
//table containing the Age grid of the SSP tables (originally in years, converted to log10(internal time units 1e12 Yrs/h))
extern float SSP_logAgeTab[SSP_NAGES];
//table containing redshift (different from the one in the code when scaling to future times)
extern float RedshiftTab[MAXSNAPS];
extern float LumTables[SSP_NAGES][SSP_NMETALLICITES][MAXSNAPS][NMAG];
extern float FilterLambda[NMAG+1];//wavelength of each filter + 1 for V-band

//for speeding up lookup in table:
#define SSP_NJUMPTAB 1000
extern int SSP_log_age_jump_table[SSP_NJUMPTAB];
extern double SSP_log_age_jump_factor;


#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define MAX_NLambdaFilter 1000
extern int NLambdaFilter[NMAG];
//VEGA
#define NLambdaVega 3303
#endif

//DUST EXTINCTION
#define ExpTauBCBulge 0.5	// constant extinction for young stars in bulges.
#define MUWIDTH  0.2
#define MUCENTER 0.3
extern long mu_seed;

#endif //COMPUTE_SPECPHOT_PROPERTIES

/*For H2 formation recipe - Not Supported*/
#define RHO_LEN 101
#define Z_LEN 13

extern size_t HighMark;

#ifdef UPDATETYPETWO
extern int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
extern int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
extern long long *IdList;
extern float *PosList, *VelList;
#endif

extern int Hashbits;
extern int NumWrittenInParallel;
extern double ScaleFactor;	// factor by which to multiply a position to get its ph index (after floring)


#ifdef USE_MEMORY_TO_MINIMIZE_IO
extern char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
extern size_t offset_auxdata, offset_treedata, offset_dbids;
extern size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
extern size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#endif

#define N_REION_Z 45
extern float Reion_z[N_REION_Z];
extern float Reion_log10_Mc[N_REION_Z];

extern FILE *tree_file;
extern FILE *treeaux_file;
extern FILE *treedbids_file;
extern FILE *FdGalTree;
extern FILE *FdGalTreeSFH;
extern FILE *FdGalDumps[NOUT];

/* used to compute redshift for given l.o.s. distance by interpolation: */
#ifndef N_REDSHIFTS_FOR_INTERPOLATION
#define N_REDSHIFTS_FOR_INTERPOLATION 150
#endif /* not defined N_REDSHIFTS_FOR_INTERPOLATION */
#ifndef DELTA_REDSHIFT_FOR_INTERPOLATION
#define DELTA_REDSHIFT_FOR_INTERPOLATION 0.1
#endif /* not defined DELTA_REDSHIFT_FOR_INTERPOLATION */
extern double distance_table_for_interpolation[N_REDSHIFTS_FOR_INTERPOLATION];
extern double redshift_table_for_interpolation[N_REDSHIFTS_FOR_INTERPOLATION];

#ifdef LIGHTCONE_OUTPUT
/* file pointers for lightcone output files: */
extern FILE *FdLightconeGalTree;
extern FILE *FdLightconeGalDumps[NOUT];

/* keeping track of #galaxies in lightcone output files: */
extern long long TotLightconeGalCount;
extern long long TotLightconeGalaxies[NOUT];

/* redshift and distance boundaries for lightcone slices: */
extern float lightcone_slice_lower_redshift [NOUT];
extern float lightcone_slice_upper_redshift [NOUT];
  
extern float lightcone_slice_lower_los_distance [NOUT];
extern float lightcone_slice_upper_los_distance [NOUT];

extern float lightcone_radius_for_snapshot       [MAXSNAPS];
extern bool  is_outside_lightcone_for_snapshot   [MAXSNAPS];
extern bool  check_outside_lightcone_for_snapshot[MAXSNAPS];

/* position of observer w.r.t. simulation box: */
extern float lightcone_observer_position[3];
extern float lightcone_observer_distance_from_origin;

/* simple lightcone geometry and stellar mass selection criteria: */
extern float lightcone_lower_redshift;
extern float lightcone_upper_redshift;
extern float lightcone_lower_ra;
extern float lightcone_upper_ra;
extern float lightcone_lower_dec;
extern float lightcone_upper_dec;
extern float lightcone_lower_stellar_mass;

extern long long lightcone_N_fof_groups_skipped_construction;
extern long long lightcone_N_galaxies_skipped_construction;
extern long long lightcone_N_galaxies_skipped_output_early;
extern long long lightcone_N_galaxies_for_output;

#endif /* defined LIGHTCONE_OUTPUT */

#endif /* header guard */
