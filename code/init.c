/*  Copyright (C) <2016+>  <L-Galaxies>
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
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#ifdef PARALLEL
#include <mpi.h>
#endif /* defined PARALLEL */

#include "allvars.h"
#include "proto.h"
#include "mcmc_vars.h"

/** @file init.c
 *  @brief Sets up some unit conversion variables; converts SN and AGN feedback
 *         variables into internal units; and reads in input tables, including
 *         the desired output snapshots, the photometric and dust tables, the
 *         cooling functions and reionization tables.
 *    
 *         <B>set_units()</B> - THIS IS FUNDAMENTAL TO UNDERSTAND THE UNITS IN
 *         THE CODE! sets ups some variables used to convert from physical to
 *         internal units (as UnitDensity_in_cgs). These are obtained from
 *         UNITLENGTH_IN_CM, UNITMASS_IN_G and UNITVELOCITY_IN_CM_PER_S.
 *    
 *         <B>read_output_snaps()</B> - reads in the list of output redshifts from
 *         file ./input/desired_output_redshifts.txt and converts them into snapsshots
 *         for the given cosmology.
 *    
 *         <B>read_zlist()</B> - reads in 1/(z+1) from FileWithZList defined
 *         in ./input/input.par and creates a table with output
 *         redshift ZZ[] and ages Age[].
 *    
 *         <B>read_file_nrs()</B> - Done if SPECIFYFILENR OFF - the dark matter files
 *         to read can be defined in a file, instead of being read sequentially.
 *         These are defined in FileNrDir, in input.par, and read into
 *         ListInputFileNr[].
 *    
 *         <B>read_reionization()</B> - Reads in Reion_z[] and Reion_Mc[] from
 *         ./input/Mc.txt. These are used if UPDATEREIONIZATION ON to get Okamoto(2008)
 *         fitting parameters (Mc), instead of Kravtsov(2004)  for the Gnedin (2000)
 *         reionization formalism.
 *    
 *    
 *         <B>read_dust_tables()</B> - Reads in the dust extinction for the same bands
 *         read from the spectrophotometric tables both for the inter-galactic medium
 *         (**_Set_Ext_table.dat) and young birth clouds (**_Set_YSExt_table.dat).
 *         Detailed description at recipe_dust.c
 *    
 *         <B>read_cooling_functions()</B> - calls the functions that read in the
 *         cooling functions defined in cool_func.c
 *    
 *         In init.c, but called from other files, are the function to interpolate
 *         properties into different tables. They are all the same but interpolate
 *         into different properties (have different inputs):
 *    
 *    
 *         <B>find_interpolated_lum()</B> - Used by add_to_luminosities() to
 *         interpolates into the age and metallicity in the SSP tables.
 *    
 *          <B>find_interpolate_reionization()</B> - Called from recipe_infall
 *         interpolates into the Reion_z[], Reion_Mc[] table, giving a value of Mc
 *         for a given redshift.
 *    
 *         SuperNovae and AGN feedback parameters are converted into internal units:
 *    
 *         \f$ AgnEfficiency = \frac{UnitMass_{\rm{g}}}{1.58e^{-26}UnitTime_{\rm{s}}}\f$
 *    
 *         \f$ EnergySNcode = \frac{EnergySN}{UnitEnergy_{\rm{cgs}}} h; \f$
           \f$ EtaSNcode = EtaSN \frac{UnitMass_{\rm{g}}}{M_\odot h}. \f$
      
 *         */


/** @brief controlling recipe for init.c, calls functions to read in tables and
 *         defines some SN and AGN feedback parameters.
 *       
 *         Calls read_output_snaps(), read_zlist(), read_recgas(),
 *         read_file_nrs(), read_sfrz(), read_reionization(), read_dust_tables() and
 *         read_cooling_functions().
 *       
 *         Converts EnergySN (->EnergySNcode), EtaSN (->EtaSNcode) and AgnEfficiency
 *         into internal units.
 *        */
void init(void)
{
  struct rlimit rlim_;

  getrlimit(RLIMIT_CORE, &rlim_);
  rlim_.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim_);

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

#ifdef GALAXYTREE
  ScaleFactor = pow(2, Hashbits) / BoxSize;
#endif

  //reads in the redshifts for the used Cosmology
  read_zlist();
  //reads in the redshifts for Original Cosmology
  read_zlist_original_cosm();

  //reads in the desired output snapshots
  read_output_snaps();

#ifdef SPECIFYFILENR
  /* read in the number of the files to be processed */
  read_file_nrs();
#endif

  if(ReionizationModel == 0)
  { read_reionization(); }

  //Values of a for the beginning and end of reionization
  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);

  read_cooling_functions();

//CREATE ARRAYS OF SFH TIME STRUCTURE:
#ifdef  STAR_FORMATION_HISTORY
  create_sfh_bins();
#endif

#ifndef MR_PLUS_MRII //If this option (to run with MCMC) is on, the tables are read later
#ifdef COMPUTE_SPECPHOT_PROPERTIES
//read in photometric tables
#ifdef PHOTTABLES_PRECOMPUTED
#ifdef MRII
  setup_LumTables_precomputed("MRII");
#else
  setup_LumTables_precomputed("MR");
#endif
#endif
#ifdef SPEC_PHOTABLES_ON_THE_FLY
  setup_Spec_LumTables_onthefly();
#endif
#endif //COMPUTE_SPECPHOT_PROPERTIES
#endif //MR_PLUS_MRII

//READ IN THE YIELD TABLES, AND FORM NORMALISED YIELD ARRAYS:
#ifdef DETAILED_METALS_AND_MASS_RETURN
  read_yield_tables();
  init_integrated_yields();
  integrate_yields();
#endif

#ifdef ASSUME_FLAT_LCDM
  assert_flat_LCDM();
#endif /* defined ASSUME_FLAT_LCDM */

  init_cosmology();
  
#ifdef LIGHTCONE_OUTPUT
  init_lightcone();
#endif /* defined LIGHTCONE_OUTPUT */
}


/** @brief Reads in 1/(z+1) from FileWithZList defined
 *         in ./input/input.par for the list of output snapshots.
 *         Creates a table with redshift ZZ[] and ages Age[].
 */
void read_zlist(void)
{
	int i_;
  FILE *file_;
  char file_name_[1000];
  char error_message_[1000];
  //double dumb2_, dumb_;

  sprintf(file_name_, "%s", FileWithZList);

  if(!(file_ = fopen(file_name_, "r")))
  {
    sprintf(error_message_,"can't read output list in file '%s'\n", file_name_);
    terminate(error_message_);
  }

  Zlistlen = 0;
  do
  {
    //if(fscanf(file_, "%d %lg %lf\n", &dumb_, &ZZ[Zlistlen], &dumb2_) == 3)
    if(fscanf(file_, "%lg\n", &AA[Zlistlen]) == 1)
      Zlistlen++;
    else
      break;
  }
  while(Zlistlen < LastDarkMatterSnapShot+1);

  fclose(file_);

  for(i_ = 0; i_ < Zlistlen; i_++)
  {
    //convert AA[] into A_inv[]
    /* AA_inv[i_] = 1 / AA[i_]; */
    //convert AA[] into redshift - ZZ[]
    ZZ[i_] = 1 / AA[i_] - 1;
    //table with time in internal units (Mpc/Km/s/h)
    Age[i_] = time_to_present(ZZ[i_]);
  }

#ifndef MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif

}


/** @brief reads redshifts from file */
void read_zlist_new(void)
{
  int i_, dummy_snap_;
  double dummy_a_;
  FILE *file_;
  char file_name_[1000];

  sprintf(file_name_, "%s", FileWithZList);
  if(!(file_ = fopen(file_name_, "r")))
  {
  	char error_message_[1000];
  	sprintf(error_message_, "can't open file `%s'\n", file_name_);
  	terminate(error_message_);
  }

  Zlistlen = 0;
  do
  {
    if(fscanf(file_, "%d %lg %lg", &dummy_snap_, &ZZ[Zlistlen], &dummy_a_) == 3)
      Zlistlen++;
    else
      break;
  }
  while(Zlistlen < LastDarkMatterSnapShot+1);
  fclose(file_);

  for(i_ = 0; i_ < Zlistlen; i_++)
  {
  	//convert redshift - ZZ[] into AA[]
                /* AA_inv[i_] = ZZ[i_] + 1; */
  	AA[i_] = 1/(ZZ[i_] + 1);
  	//printf("z[%d]=%f\n",i_,ZZ[i_]);
  	//table with time in internal units (Mpc/Km/s/h)
  	if(ZZ[i_]>=0.0)
  		Age[i_] = time_to_present(ZZ[i_]);
  	else
  		Age[i_] = 0.0;
  		//break;
  }

#ifndef MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif
}


/** @brief reads original cosmology redshifts from file */
void read_zlist_original_cosm(void)
{
  FILE *file_;
  char file_name_[1000];

  sprintf(file_name_, "%s", FileWithZList_OriginalCosm);
  if(!(file_ = fopen(file_name_, "r")))
  {
    printf("can't read output list in file '%s'\n", file_name_);
    terminate("in read_zlist_original_cosm");
  }

  Zlistlen = 0;
  do
  {
    if(fscanf(file_, " %lg ", &AA_OriginalCosm[Zlistlen]) == 1)
    {  Zlistlen++;	}
    else
      break;
  }
  while(Zlistlen < LastDarkMatterSnapShot+1);

  fclose(file_);
  
/*  int i_;
  for(i_ = 0; i_ < Zlistlen; i_++)
  { AA_OriginalCosm_inv[i_] = 1 / AA_OriginalCosm[i_]; } */

#ifndef MR_PLUS_MRII //option for MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif
}


/** @brief Reads in the list of output snapshots from
 *         file /input/desired_output_snaps.txt
 */
void read_output_snaps(void)
{
#ifndef GALAXYTREE
  int output_number_, snapshot_number_;
  char file_name_[1000];
  FILE *file_;

  LastSnapShotNr=0;

  sprintf(file_name_, "%s", FileWithOutputRedshifts);

  if(!(file_ = fopen(file_name_, "r")))
  {
    char error_message_[1000];
    sprintf(error_message_, "file `%s' not found.\n", file_name_);
    terminate(error_message_);
  }

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    if(fscanf(file_, " %f ", &ListOutputRedshifts[output_number_]) != 1)
    {
      char error_message_[1000];
      sprintf(error_message_, "I/O error in file '%s'\n", file_name_);
      terminate(error_message_);
    }

    //find the snapshot corresponding to the desired output redshift in ListOutputRedshifts[]
    for(snapshot_number_ = 0; snapshot_number_ <= LastDarkMatterSnapShot; snapshot_number_++)
      if(ListOutputRedshifts[output_number_]>=ZZ[snapshot_number_])
      {
        if((snapshot_number_ > 0) && ((ZZ[snapshot_number_-1] - ListOutputRedshifts[output_number_]) < (ListOutputRedshifts[output_number_]-ZZ[snapshot_number_]) || (ZZ[snapshot_number_] < 0.)))
          ListOutputSnaps[output_number_]=snapshot_number_-1;
        else
          ListOutputSnaps[output_number_]=snapshot_number_;
        break;
      }

#ifdef MCMC
        if (ThisTask == 0 && CurrentMCMCStep==1)
#else /* not defined MCMC */
        if (ThisTask == 0)
#endif /* not defined MCMC */
          printf("output %d: requested z=%0.2f, available snap[%d] z=%f & snap[%d] z=%f, use snap[%d]\n",
              output_number_, ListOutputRedshifts[output_number_], snapshot_number_-1, ZZ[snapshot_number_-1], snapshot_number_, ZZ[snapshot_number_], ListOutputSnaps[output_number_]);

    //define the LastSnapShotNr as the highest snapshot need to be computed
    if(LastSnapShotNr<ListOutputSnaps[output_number_])
    LastSnapShotNr=ListOutputSnaps[output_number_];
  }
  fclose(file_);
  
  for(snapshot_number_ = 0; snapshot_number_ < MAXSNAPS; snapshot_number_++)
  {
    for( output_number_ = NOUT; (output_number_ > 0) && (output_number_--) && (ListOutputSnaps[output_number_] < snapshot_number_); ) {}
    ListOutputNumberOfSnapshot[snapshot_number_] = output_number_;
  }
  
// #ifdef MCMC
//   if (ThisTask == 0 && CurrentMCMCStep==1)
// #else /* not defined MCMC */
//   if (ThisTask == 0)
// #endif /* not defined MCMC */
//   {
//     printf("\n----------------\n");
//     for(snapshot_number_ = 0; snapshot_number_ < MAXSNAPS; snapshot_number_++)
//       printf("ListOutputNumberOfSnapshot[%d] = %d\n", snapshot_number_, ListOutputNumberOfSnapshot[snapshot_number_]);
//     printf("----------------\n\n");
//   }
  

#else /* defined GALAXYTREE */
  int output_number_, snapshot_number_;
  
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    ListOutputSnaps[output_number_] = output_number_;
  
  for(snapshot_number_ = 0; snapshot_number_ < MAXSNAPS; snapshot_number_++)
    ListOutputNumberOfSnapshot[snapshot_number_] = (snapshot_number_ < NOUT) ? snapshot_number_ : NOUT - 1;
  
  LastSnapShotNr=LastDarkMatterSnapShot;
#endif /* defined GALAXYTREE */
}


#ifdef SPECIFYFILENR
/** @brief reads file numbers to process from file */
void read_file_nrs(void)
{
  int i_;
  char file_name_[1000];
  FILE *file_;
  sprintf(file_name_, "%s", FileNrDir);
  if(!(file_ = fopen(file_name_, "r")))
  {
    char error_message_[1000];
    sprintf(error_message_, "file `%s' not found.\n", file_name_);
    terminate(error_message_);
  }

  for(i_ = 0; i_ < 111; i_++) //only 111files in ../input/filenrdir.txt are read in
  {
    if(fscanf(file_, " %d ", &ListInputFileNr[i_]) != 1)
    {
      char error_message_[1000];
      sprintf(error_message_, "I/O error in file '%s'\n", file_name_);
      terminate(error_message_);
    }
  }
  fclose(file_);
}
#endif


/** @brief reads reionization parameters from file */
void read_reionization(void)
{
  FILE *file_;
  int p_;
  float dumb_;

  if(!(file_ = fopen(McFile, "r")))
  {
    char error_message_[1000];
    sprintf(error_message_, "file `%s' not found.\n", McFile);
    terminate(error_message_);
  }

  for(p_ = 0; p_ < N_REION_Z; p_++)
  {
    fscanf(file_, "%f", &dumb_);
    fscanf(file_, "%f", &Reion_z[p_]);
    fscanf(file_, "%f", &Reion_log10_Mc[p_]); 
    Reion_log10_Mc[p_] = log10(Reion_log10_Mc[p_]);  /* file should contain reion. Mc, not log10(Mc) */
  }
  fclose(file_);
}


/** @brief determines number of files for process */
int get_nr_files_to_process()
{
#ifdef OVERWRITE_OUTPUT
  return LastFile - FirstFile + 1;
#else /* not defined OVERWRITE_OUTPUT */
  int n_files_ = 0;
  if(ThisTask==0)
  {  
    int file_number_;
    for(file_number_ = FirstFile; file_number_ <= LastFile; file_number_++)
    {
#ifdef SPECIFYFILENR
      const int input_file_number_ = ListInputFileNr[file_number_];
#else /* not defined SPECIFYFILENR */
      const int input_file_number_ = file_number_;
#endif /* not defined SPECIFYFILENR */

      char file_name_[1000];
#ifdef GALAXYTREE
      sprintf(file_name_, "%s/%s_galtree_%d", FinalOutputDir, FileNameGalaxies, input_file_number_);
#else /* not defined GALAXYTREE */
      sprintf(file_name_, "%s/%s_z%1.2f_%d", FinalOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], input_file_number_);
#endif /* not defined GALAXYTREE */
      struct stat file_status_;
      if(stat(file_name_, &file_status_) != 0)  // seems not to exist  
      { ++n_files_; }
    }
  }
 #ifdef PARALLEL
  MPI_Bcast(&n_files_,1, MPI_INT, 0, MPI_COMM_WORLD);
#endif /* defined PARALLEL */
  return n_files_;
#endif /* not defined OVERWRITE_OUTPUT */
}


/** @brief distributes tree files over available tasks */
void assign_files_to_tasks(int *FileToProcess_, int *TaskToProcess_, int n_files_)
{
  int i_,j_, file_number_, input_file_number_;

  if(ThisTask==0)
  {
    i_=0;
    j_=0;
    for(file_number_ = FirstFile; file_number_ <= LastFile; file_number_++)
    {
#ifdef SPECIFYFILENR
      input_file_number_ = ListInputFileNr[file_number_];
#else
      input_file_number_ = file_number_;
#endif
#ifndef OVERWRITE_OUTPUT
      char file_name_[1000];
#ifdef GALAXYTREE
      sprintf(file_name_, "%s/%s_galtree_%d", FinalOutputDir, FileNameGalaxies, input_file_number_);
#else
      sprintf(file_name_, "%s/%s_z%1.2f_%d", FinalOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], input_file_number_);
#endif
      struct stat file_status_;
      if(stat(file_name_, &file_status_) != 0)        // doesn't exist
      {
#endif
        FileToProcess_[i_]=input_file_number_;
#ifdef PARALLEL
        TaskToProcess_[i_]=j_;
#else
        TaskToProcess_[i_]=0;
#endif
        ++i_;
        ++j_;
        if(j_==NTask)
        { j_=0; }
#ifndef OVERWRITE_OUTPUT
      }
#endif
    }
  }
#ifdef PARALLEL
  MPI_Bcast(FileToProcess, sizeof(int) * n_files_, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(TaskToProcess, sizeof(int) * n_files_, MPI_BYTE, 0, MPI_COMM_WORLD);
#else /* not defined PARALLEL */
  (void)n_files_; /* avoid unused-parameter warning */
#endif /* not defined PARALLEL */
}
