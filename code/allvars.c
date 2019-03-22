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


#include "allvars.h"

const char* GasComponentStr    [] = {"HotGas", "ColdGas", "EjectedGas"};
const char* StellarComponentStr[] = {"Disk", "Bulge", "ICM", "Burst"};

struct GALAXY			/* Galaxy data */
 *Gal, *HaloGal;

struct halo_data *Halo, *Halo_Data;

struct halo_aux_data		/* auxiliary halo data */
 *HaloAux;

struct halo_ids_data *HaloIDs, *HaloIDs_Data;

int FirstFile;			/* first and last file for processing */
int LastFile;


double AllocValue_MaxHaloGal;
double AllocValue_MaxGal;
double AllocValue_MaxGalTree;

int Ntrees;			/* number of trees in current file */

int MaxGal;
int NHaloGal, MaxHaloGal;
int NGalTree, MaxGalTree;
int *HaloGalHeap;
int IndexStored;

char SpecPhotDir[512];
char PhotPrefix[50];
char SpecPhotIMF[50];
char McFile[512];
char FileWithFilterNames[512];
char CoolFunctionsDir[512];
char OutputDir[512];
char FinalOutputDir[512];
char FileNameGalaxies[512];
char SimulationDir[512];
char FileWithOutputRedshifts[512];
char FileWithZList[512];
//variables used to scale to a different cosmology
char FileWithZList_OriginalCosm[512];
#ifdef MR_PLUS_MRII  //OPTION for MCMC
char FileWithZList_MR[512];
char FileWithZList_OriginalCosm_MR[512];
char FileWithZList_MRII[512];
char FileWithZList_OriginalCosm_MRII[512];
#endif

double ScalePos, ScaleMass;

#ifdef SPECIFYFILENR
char FileNrDir[512];
int ListInputFileNr[111];
#endif


int TotHalos;
int TotGalaxies[NOUT];
int *TreeNgals[NOUT];

int LastSnapShotNr;
int LastDarkMatterSnapShot;
#ifdef MR_PLUS_MRII  //OPTION for MCMC
int LastDarkMatterSnapShot_MR;
int LastDarkMatterSnapShot_MRII;
#endif


int *FirstHaloInSnap;
int *TreeNHalos;
int *TreeFirstHalo;

double MaxMemSize;

size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;


#ifdef PARALLEL
int ThisTask, NTask;
#endif /* defined PARALLEL */

#ifdef GALAXYTREE
int GalCount;
int TotGalCount;
struct galaxy_tree_data *GalTree;
#endif

size_t HighMark;


/* cosmological parameters */
double BaryonFrac;
double Sigma8;
double Omega;
double OmegaLambda;
double Hubble_h;
double inv_Hubble_h;
double Omega_OriginalCosm;
double OmegaLambda_OriginalCosm;
double Hubble_h_OriginalCosm;
//SIMULATION RELATED
double PartMass;
double BoxSize;
double PartMass_OriginalCosm;
double BoxSize_OriginalCosm;

#ifdef MR_PLUS_MRII  //OPTION for MCMC
double PartMass_MR;
double BoxSize_MR;
double PartMass_OriginalCosm_MR;
double BoxSize_OriginalCosm_MR;
double PartMass_MRII;
double BoxSize_MRII;
double PartMass_OriginalCosm_MRII;
double BoxSize_OriginalCosm_MRII;
#endif


/* flags */
int ReionizationModel;
int DiskRadiusModel;
int StarFormationModel;
int FeedbackReheatingModel;
int FeedbackEjectionModel;
int FateOfSatellitesGas;
int ReIncorporationModel;
int AGNRadioModeModel;
int DiskInstabilityModel;
int BHGrowthInDiskInstabilityModel;
int HotGasStrippingModel;
int DisruptionModel;
int StarBurstModel;
int BulgeFormationInMinorMergersOn;
int MetallicityOption;

/* parameters */
double Reionization_z0;
double Reionization_zr;
double RamPressureStrip_CutOffMass;
double SfrEfficiency;
double SfrColdCrit;
double SfrBurstEfficiency;
double SfrBurstSlope;
double Yield;
double RecycleFraction;
double ThreshMajorMerger;
double MergerTimeMultiplier;
double AgnEfficiency;
double BlackHoleGrowthRate;
double BlackHoleSeedMass;
double BlackHoleCutoffVelocity;
double FeedbackReheatingEpsilon;
double ReheatPreVelocity;
double ReheatSlope;
double FeedbackEjectionEfficiency;
double EjectPreVelocity;
double EjectSlope;
double ReIncorporationFactor;
double EnergySNcode, EnergySN;
double EtaSNcode, EtaSN;

// internal units and phys. consts. in internal units moved to "physical_constants_and_units.h"
// double UnitTime_in_s;
// double UnitPressure_in_cgs;
// double UnitDensity_in_cgs;
// double UnitCoolingRate_in_cgs;
// double UnitEnergy_in_cgs;
// double UnitTime_in_Megayears; //Using time as stored in the code, this gives Myr/h
// double UnitTime_in_years;
// 
// double Gravity;
// double SpeedOfLight;
// double Hubble;
// double RhoCrit;

double  a0, ar;

int ListOutputSnaps[NOUT];
float ListOutputRedshifts[NOUT];

int ListOutputNumberOfSnapshot[MAXSNAPS];

double ZZ[MAXSNAPS];
double AA[MAXSNAPS];
//variable used to scale to a different cosmology
double AA_OriginalCosm[MAXSNAPS];

double Age[MAXSNAPS];

int Zlistlen;

gsl_rng *random_generator;

int NumMergers;

/*  tabulated stuff */

/* fixed-metallicity spectrophotometric model */
/* tables hold magnitues of starburst population as a function of age */

#ifdef STAR_FORMATION_HISTORY
double SFH_t[MAXSNAPS][STEPS][SFH_NBIN];
double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN];
int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN];
int SFH_ibin[MAXSNAPS][STEPS];
#ifdef DETAILED_METALS_AND_MASS_RETURN
double tau_t[STEPS*MAXSNAPS]; //Time-to-z=0 of every timestep in the code. (Used for SNe rates in yield_integrals.c)
double tau_dt[STEPS*MAXSNAPS];//Width of every timestep in the code. (Used for SNe rates in yield_integrals.c)
#endif
#endif //STAR_FORMATION_HISTORY

#ifdef COMPUTE_SPECPHOT_PROPERTIES
//SSP PHOT TABLES
float SSP_logMetalTab[SSP_NMETALLICITES];
float SSP_logAgeTab[SSP_NAGES];
float RedshiftTab[MAXSNAPS];
float LumTables[SSP_NAGES][SSP_NMETALLICITES][MAXSNAPS][NMAG];
float FilterLambda[NMAG+1];	//wavelength of each filter + 1 for V-band
#ifdef SPEC_PHOTABLES_ON_THE_FLY
int NLambdaFilter[NMAG];
#endif

//for speeding up lookup in table: 
int SSP_log_age_jump_table[SSP_NJUMPTAB];
double SSP_log_age_jump_factor;

// dust
long mu_seed;
#endif

void *TreeAuxData;

#ifdef UPDATETYPETWO
int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
long long *IdList;
float *PosList, *VelList;
#endif


int Hashbits;
int NumWrittenInParallel;
double ScaleFactor;


#ifdef USE_MEMORY_TO_MINIMIZE_IO
char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
size_t offset_auxdata, offset_treedata, offset_dbids;
size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#endif

/* reionization Okamoto et al. 2008*/
float Reion_z[N_REION_Z];
float Reion_log10_Mc[N_REION_Z];

FILE *tree_file;
FILE *treeaux_file;
FILE *treedbids_file;
FILE *FdGalTree;
FILE *FdGalTreeSFH;
FILE *FdGalDumps[NOUT];

/* used to compute redshift for given l.o.s. distance by interpolation */
double distance_table_for_interpolation[N_REDSHIFTS_FOR_INTERPOLATION];
double redshift_table_for_interpolation[N_REDSHIFTS_FOR_INTERPOLATION];

/* for lightcone output */
#ifdef LIGHTCONE_OUTPUT
/* file pointers for lightcone output files: */
FILE *FdLightconeGalTree;
FILE *FdLightconeGalDumps[NOUT];

/* keeping track of #galaxies in lightcone output files: */
long long TotLightconeGalCount;
long long TotLightconeGalaxies[NOUT];

/* redshift and distance boundaries for lightcone slices: */
float lightcone_slice_lower_redshift [NOUT];
float lightcone_slice_upper_redshift [NOUT];
  
float lightcone_slice_lower_los_distance [NOUT];
float lightcone_slice_upper_los_distance [NOUT];

float lightcone_radius_for_snapshot       [MAXSNAPS];
bool  is_outside_lightcone_for_snapshot   [MAXSNAPS];
bool  check_outside_lightcone_for_snapshot[MAXSNAPS];

/* position of observer w.r.t. simulation box: */
float lightcone_observer_position[3];
float lightcone_observer_distance_from_origin;

/* simple lightcone geometry and stellar mass selection criteria: */
float lightcone_lower_redshift;
float lightcone_upper_redshift;
float lightcone_lower_ra;
float lightcone_upper_ra;
float lightcone_lower_dec;
float lightcone_upper_dec;
float lightcone_lower_stellar_mass;

long long lightcone_N_fof_groups_skipped_construction;
long long lightcone_N_galaxies_skipped_construction;
long long lightcone_N_galaxies_skipped_output_early;
long long lightcone_N_galaxies_for_output;
long long lightcone_N_galaxies_remaining_for_output_past_construct_galaxies;


#endif /* defined LIGHTCONE_OUTPUT */
