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
 *         ListInputFilrNr[].
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
  struct rlimit rlim;

  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

#ifdef GALAXYTREE
  ScaleFactor = pow(2, Hashbits) / BoxSize;
#endif

  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;
  EtaSNcode = EtaSN * (UNITMASS_IN_G / SOLAR_MASS) / Hubble_h;

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
	int i;
  FILE *fd;
  char fname[1000];
  char sbuf[1000];
  //double dumb2, dumb;

  sprintf(fname, "%s", FileWithZList);

  if(!(fd = fopen(fname, "r")))
  {
    sprintf(sbuf,"can't read output list in file '%s'\n", fname);
    terminate(sbuf);
  }

  Zlistlen = 0;
  do
  {
    //if(fscanf(fd, "%d %lg %lf\n", &dumb, &ZZ[Zlistlen], &dumb2) == 3)
    if(fscanf(fd, "%lg\n", &AA[Zlistlen]) == 1)
      Zlistlen++;
    else
      break;
  }
  while(Zlistlen < LastDarkMatterSnapShot+1);

  fclose(fd);

  for(i = 0; i < Zlistlen; i++)
  {
    //convert AA[] into A_inv[]
    /* AA_inv[i] = 1 / AA[i]; */
    //convert AA[] into redshift - ZZ[]
    ZZ[i] = 1 / AA[i] - 1;
    //table with time in internal units (Mpc/Km/s/h)
    Age[i] = time_to_present(ZZ[i]);
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
  int i, dummy_snap;
  double dummy_a;
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithZList);
  if(!(fd = fopen(fname, "r")))
  {
  	char sbuf[1000];
  	sprintf(sbuf, "can't open file `%s'\n", fname);
  	terminate(sbuf);
  }

  Zlistlen = 0;
  do
  {
    if(fscanf(fd, "%d %lg %lg", &dummy_snap, &ZZ[Zlistlen], &dummy_a) == 3)
      Zlistlen++;
    else
      break;
  }
  while(Zlistlen < LastDarkMatterSnapShot+1);
  fclose(fd);

  for(i = 0; i < Zlistlen; i++)
  {
  	//convert redshift - ZZ[] into AA[]
                /* AA_inv[i] = ZZ[i] + 1; */
  	AA[i] = 1/(ZZ[i] + 1);
  	//printf("z[%d]=%f\n",i,ZZ[i]);
  	//table with time in internal units (Mpc/Km/s/h)
  	if(ZZ[i]>=0.0)
  		Age[i] = time_to_present(ZZ[i]);
  	else
  		Age[i] = 0.0;
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
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithZList_OriginalCosm);
  if(!(fd = fopen(fname, "r")))
  {
    printf("can't read output list in file '%s'\n", fname);
    terminate("in read_zlist_original_cosm");
  }

  Zlistlen = 0;
  do
  {
    if(fscanf(fd, " %lg ", &AA_OriginalCosm[Zlistlen]) == 1)
    {  Zlistlen++;	}
    else
      break;
  }
  while(Zlistlen < LastDarkMatterSnapShot+1);

  fclose(fd);
  
/*  int i;
  for(i = 0; i < Zlistlen; i++)
  { AA_OriginalCosm_inv[i] = 1 / AA_OriginalCosm[i]; } */

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
  int i, j;
  char buf[1000];
  FILE *fd;

  LastSnapShotNr=0;

  sprintf(buf, "%s", FileWithOutputRedshifts);

  if(!(fd = fopen(buf, "r")))
  {
    char sbuf[1000];
    sprintf(sbuf, "file `%s' not found.\n", buf);
    terminate(sbuf);
  }

  for(i = 0; i < NOUT; i++)
  {
    if(fscanf(fd, " %f ", &ListOutputRedshifts[i]) != 1)
    {
      char sbuf[1000];
      sprintf(sbuf, "I/O error in file '%s'\n", buf);
      terminate(sbuf);
    }

    //find the snapshot corresponding to the desired output redshift in ListOutputRedshifts[]
    for(j = 0; j <= LastDarkMatterSnapShot; j++)
      if(ListOutputRedshifts[i]>=ZZ[j])
      {
        if((j > 0) && ((ZZ[j-1] - ListOutputRedshifts[i]) < (ListOutputRedshifts[i]-ZZ[j]) || (ZZ[j] < 0.)))
          ListOutputSnaps[i]=j-1;
        else
          ListOutputSnaps[i]=j;
        break;
      }

#ifdef MCMC
        if (ThisTask == 0 && CurrentMCMCStep==1)
#else /* not defined MCMC */
        if (ThisTask == 0)
#endif /* not defined MCMC */
          printf("output %d: requested z=%0.2f, available snap[%d] z=%f & snap[%d] z=%f, use snap[%d]\n",
              i, ListOutputRedshifts[i], j-1, ZZ[j-1], j, ZZ[j], ListOutputSnaps[i]);

    //define the LastSnapShotNr as the highest snapshot need to be computed
    if(LastSnapShotNr<ListOutputSnaps[i])
    LastSnapShotNr=ListOutputSnaps[i];
  }
  fclose(fd);
  
  for(j = 0; j < MAXSNAPS; j++)
  {
    for( i = NOUT; (i > 0) && (i--) && (ListOutputSnaps[i] < j); ) {}
    ListOutputNumberOfSnapshot[j] = i;
  }
  
// #ifdef MCMC
//   if (ThisTask == 0 && CurrentMCMCStep==1)
// #else /* not defined MCMC */
//   if (ThisTask == 0)
// #endif /* not defined MCMC */
//   {
//     printf("\n----------------\n");
//     for(j = 0; j < MAXSNAPS; j++)
//       printf("ListOutputNumberOfSnapshot[%d] = %d\n", j, ListOutputNumberOfSnapshot[j]);
//     printf("----------------\n\n");
//   }
  

#else /* defined GALAXYTREE */
  int i;
  for(i = 0; i < NOUT; i++)
    ListOutputSnaps[i] = i;
  for(i = 0; i < MAXSNAPS; i++)
    ListOutputNumberOfSnapshot[i] = (i < NOUT) ? i : NOUT - 1;
  LastSnapShotNr=LastDarkMatterSnapShot;
#endif /* defined GALAXYTREE */
}


#ifdef SPECIFYFILENR
/** @brief reads file numbers to process from file */
void read_file_nrs(void)
{
  int i;
  char buf[1000];
  FILE *fd;
  sprintf(buf, "%s", FileNrDir);
  if(!(fd = fopen(buf, "r")))
  {
    char sbuf[1000];
    sprintf(sbuf, "file `%s' not found.\n", buf);
    terminate(sbuf);
  }

  for(i = 0; i < 111; i++) //only 111files in ../input/filenrdir.txt are read in
  {
    if(fscanf(fd, " %d ", &ListInputFilrNr[i]) != 1)
    {
      char sbuf[1000];
      sprintf(sbuf, "I/O error in file '%s'\n", buf);
      terminate(sbuf);
    }
  }
  fclose(fd);
}
#endif


/** @brief reads reionization parameters from file */
void read_reionization(void)
{
  FILE *fd;
  int p;
  float dumb;

  if(!(fd = fopen(McFile, "r")))
  {
    char sbuf[1000];
    sprintf(sbuf, "file `%s' not found.\n", McFile);
    terminate(sbuf);
  }

  for(p = 0; p < N_REION_Z; p++)
  {
    fscanf(fd, "%f", &dumb);
    fscanf(fd, "%f", &Reion_z[p]);
    fscanf(fd, "%f", &Reion_log10_Mc[p]); 
    Reion_log10_Mc[p] = log10(Reion_log10_Mc[p]);  /* file should contain reion. Mc, not log10(Mc) */

  }
  fclose(fd);
}


/** @brief determines number of files for process */
int get_nr_files_to_process()
{
#ifdef OVERWRITE_OUTPUT
  return LastFile - FirstFile + 1;
#else /* not defined OVERWRITE_OUTPUT */
  int nfiles = 0;
  if(ThisTask==0)
  {  
    int filenr;
    for(filenr = FirstFile; filenr <= LastFile; filenr++)
    {
#ifdef SPECIFYFILENR
      const int file = ListInputFilrNr[filenr];
#else /* not defined SPECIFYFILENR */
      const int file = filenr;
#endif /* not defined SPECIFYFILENR */

      char buf[1000];
#ifdef GALAXYTREE
      sprintf(buf, "%s/%s_galtree_%d", FinalOutputDir, FileNameGalaxies, file);
#else /* not defined GALAXYTREE */
      sprintf(buf, "%s/%s_z%1.2f_%d", FinalOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], file);
#endif /* not defined GALAXYTREE */
      struct stat filestatus;
      if(stat(buf, &filestatus) != 0)  // seems not to exist  
      { ++nfiles; }
    }
  }
 #ifdef PARALLEL
  MPI_Bcast(&nfiles,1, MPI_INT, 0, MPI_COMM_WORLD);
#endif /* defined PARALLEL */
  return nfiles;
#endif /* not defined OVERWRITE_OUTPUT */
}


/** @brief distributes tree files over available tasks */
void assign_files_to_tasks(int *FileToProcess, int *TaskToProcess, int nfiles)
{
  int i,j, filenr, file;

  if(ThisTask==0)
  {
    i=0;
    j=0;
    for(filenr = FirstFile; filenr <= LastFile; filenr++)
    {
#ifdef SPECIFYFILENR
      file = ListInputFilrNr[filenr];
#else
      file=filenr;
#endif
#ifndef OVERWRITE_OUTPUT
      char buf[1000];
#ifdef GALAXYTREE
      sprintf(buf, "%s/%s_galtree_%d", FinalOutputDir, FileNameGalaxies, file);
#else
      sprintf(buf, "%s/%s_z%1.2f_%d", FinalOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], file);
#endif
      struct stat filestatus;
      if(stat(buf, &filestatus) != 0)        // doesn't exist
      {
#endif
        FileToProcess[i]=file;
#ifdef PARALLEL
        TaskToProcess[i]=j;
#else
        TaskToProcess[i]=0;
#endif
        ++i;
        ++j;
        if(j==NTask)
        { j=0; }
#ifndef OVERWRITE_OUTPUT
      }
#endif
    }
  }
#ifdef PARALLEL
  MPI_Bcast(FileToProcess, sizeof(int) * nfiles, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(TaskToProcess, sizeof(int) * nfiles, MPI_BYTE, 0, MPI_COMM_WORLD);
#else /* not defined PARALLEL */
  (void)nfiles; /* avoid unused-parameter warning */
#endif /* not defined PARALLEL */
}
