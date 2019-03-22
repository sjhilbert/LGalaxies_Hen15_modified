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
 *  Created in: 2009
 *      Author: Chiara Tonini & Bruno Henriques
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"


#ifdef COMPUTE_SPECPHOT_PROPERTIES
/**@brief Reads in the look up tables from Stellar Population Synthesis Models.
 *
 * Reads in the look up tables from a given Stellar Population Synthesis Model.
 * There are different files, each one corresponding to a different metallicity
 * and PhotBand. On each file the tables have the Magnitudes for single bursts of
 * 1\f$M_{\odot}\f$ for a range of ages and "inversely" k-corrected in case the
 * code needs to compute observed frame magnitudes. The number of filters is
 * givens by NMAGS and the names by input file Filter_Names.txt
 *
 * On each file, the three first numbers correspond to: the number of snapshots,
 * the number of filter bands and the number of age bins. Then the age grid is listed.
 * After that for each snapshot, for each age the value corresponding to a burst
 * inversely k-corrected to that redshift, with that age (and metallicity) is listed.
 *
 * The basic structure is number of columns  corresponding to number of mags,
 * repeated over the number of ages on the initial listed grid, multiplied by the
 * number of snapshots, with the snapshot number listed in between.
 *
 * agTableZz[Mag][mettallicity][Snapshot][Age]. */
#ifdef PHOTTABLES_PRECOMPUTED
void setup_LumTables_precomputed(char SimName[])
{
  FilterLambda[NMAG] = 0.55;        //to use by the dust model for birth clouds, the wavelength of the V-filter_number_
  
  FILE *fa, *fb;
  int MetalLoop, AgeLoop, filter_number_, snap;
  char buf[1000], FilterName[100], dummy[100], SSP[1000];
  char dumb_FilterFile[100];
  float dumb_filterlambda;
  int dumb_ssp_nsnaps, dumb_ssp_nage, dumb_ssp_nmetallicites, dumb_nmag;

#ifdef BC03
  sprintf(SSP, "BC03");
#endif
#ifdef M05
  sprintf(SSP, "M05");
#endif
#ifdef CB07
  sprintf(SSP, "CB07");
#endif

  /*Read list of metallicities available from SSP*/
  sprintf(buf, "%s/PhotTables/%s_%s_Metallicity_list.dat", SpecPhotDir, SSP, SpecPhotIMF);
  if(!(fa = fopen(buf, "r")))
  {
    char sbuf[1000];
    sprintf(sbuf, "file `%s' not found.\n", buf);
    terminate(sbuf);
  }

  fscanf(fa, "%d", &dumb_ssp_nmetallicites);
  if(dumb_ssp_nmetallicites != SSP_NMETALLICITES)
  {
    terminate("nmetallicites on file not equal to SSP_NMETALLICITES");
  }

  for(MetalLoop=0;MetalLoop<SSP_NMETALLICITES;MetalLoop++)
  {
    fscanf(fa, "%f", &SSP_logMetalTab[MetalLoop]);
    SSP_logMetalTab[MetalLoop]=log10(SSP_logMetalTab[MetalLoop]);
  }
  fclose(fa);

  /*Loop over the different files corresponding to different metallicities */
  for(MetalLoop = 0; MetalLoop < SSP_NMETALLICITES; MetalLoop++)
  {
    sprintf(buf, "%s", FileWithFilterNames);
    if((fa = fopen(buf, "r")) == NULL)
    {
      printf("\n**Can't open file \"%s\" **\n", buf);
      char sbuf[1000];
      sprintf(sbuf, "Can't open file %s\n", buf);
      terminate(sbuf);
    }

    fscanf(fa, "%d", &dumb_nmag);
    if(dumb_nmag != NMAG)
    {
      char sbuf[1000];
      sprintf(sbuf,"nmag = %d on file %s not equal to NMAG = %d",dumb_nmag, buf, NMAG);
      terminate(sbuf);
    }

    //There is a different file for each filter_number_
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      fscanf(fa,"%s %f %s" ,dumb_FilterFile, &dumb_filterlambda, FilterName);
      //READ TABLES
      sprintf(buf, "%s/PhotTables/%s_%s_Phot_Table_%s_Mag%s_m%0.4f.dat", SpecPhotDir, PhotPrefix, SpecPhotIMF, SimName, FilterName,
              pow(10,SSP_logMetalTab[MetalLoop]));
      if(!(fb = fopen(buf, "r")))
      {
        char sbuf[1000];
        sprintf(sbuf, "file `%s' not found.\n", buf);
        terminate(sbuf);
      }

      if(ThisTask == 0)
        printf("reading file %s \n", buf);

      fscanf(fb, "%s %f %d %d", dummy, &FilterLambda[filter_number_], &dumb_ssp_nsnaps, &dumb_ssp_nage);
      /* check that the numbers on top of the file correspond to (LastDarkMatterSnapShot+1) and SSP_NAGES */
      if(FilterLambda[filter_number_] != dumb_filterlambda)
      {
        terminate("filterlambdas dont match");
      }
      if(dumb_ssp_nsnaps != (LastDarkMatterSnapShot+1))
      {
        terminate("nsnaps not equal to Snaplistlen");
      }
      if(dumb_ssp_nage != SSP_NAGES)
      {
        terminate("n_age_bins not equal to SSP_NAGES");
      }

      /* read ages of SSPs (done for all filter_number_s and metallicities)
        * last call stays on global array  SSP_logAgeTab*/
      for(AgeLoop = 0; AgeLoop < SSP_NAGES; AgeLoop++)
      {
        fscanf(fb, " %e ", &SSP_logAgeTab[AgeLoop]);

        if(SSP_logAgeTab[AgeLoop] > 0.0)        // avoid taking a log of 0 ...
        {
          /* converts SSP_AgeTab from years to log10(internal time units) */
          SSP_logAgeTab[AgeLoop] = SSP_logAgeTab[AgeLoop] / 1.0e6 / UnitTime_in_Megayears * Hubble_h;
          SSP_logAgeTab[AgeLoop] = log10(SSP_logAgeTab[AgeLoop]);
        }
        else
          SSP_logAgeTab[AgeLoop] = 0.;
      }

      //read luminosities at each output redshift
      for(snap = 0; snap < (LastDarkMatterSnapShot+1); snap++)
      {
        fscanf(fb, " %f ", &RedshiftTab[snap]);
        //for each age
        for(AgeLoop = 0; AgeLoop < SSP_NAGES; AgeLoop++)
        {
          fscanf(fb, "%e", &LumTables[filter_number_][MetalLoop][snap][AgeLoop]);
          LumTables[filter_number_][MetalLoop][snap][AgeLoop] = pow(10., -LumTables[filter_number_][MetalLoop][snap][AgeLoop] / 2.5);
        }                //end loop on age
      }                //end loop on redshift (everything done for current filter_number_)

      fclose(fb);
    }//end loop on filter_number_s

    fclose(fa);
  }//end loop on metallicity

  init_SSP_log_age_jump_index();
}
#endif //PHOTTABLES_PRECOMPUTED


/** Reads in the Full SEDs from a given stellar population and filter curves and
 * computes PhotTables on the fly. Just needs the SEDs in SpecPhotDir/FullSEDs/,
 * the file with filter names and wave_lengths "FileWithFilterNames" and the
 * filter curves in SpecPhotDir/Filters/
 *
 * Developed by Chiara Tonini, adapted by Bruno Henriques
 * */
#ifdef SPEC_PHOTABLES_ON_THE_FLY
void setup_Spec_LumTables_onthefly(void)
{
  FilterLambda[NMAG] = 0.55;        // used by the dust model for birth clouds, the wavelength of the V-filter_number_
  
  double AbsMAG;
  //FILTERS
  double LambdaFilter[NMAG][MAX_NLambdaFilter], FluxFilter[NMAG][MAX_NLambdaFilter];
  double *FluxFilterOnGrid, FluxFilterInt;
  //InputSSP spectra
  double LambdaInputSSP[SSP_NAGES][SSP_NLambda], FluxInputSSP[SSP_NAGES][SSP_NLambda];
  double *FluxInputSSPConv, *FluxInputSSPOnGrid, FluxInputSSPInt;
  //VEGA
  double LambdaVega[NLambdaVega], FluxVega[NLambdaVega];
  double *FluxVegaConv, *FluxVegaOnGrid, FluxVegaInt;
  //AUX ARRAYS for FILTERS and InputSSP
  double LambdaFilter_SingleFilter[MAX_NLambdaFilter], FluxFilter_SingleFilter[MAX_NLambdaFilter], LambdaInputSSP_SingleAge[SSP_NLambda];
  double redshift;
  //Loops in the code
  int MetalLoop, AgeLoop, snap, filter_number_, i;
  FILE *File_PhotTables[NMAG];

#ifdef PARALLEL
  if(ThisTask == 0)
          printf("\n\nComputing PhotTables on the fly...\n\n");
#else
  printf("\n\nComputing PhotTables on the fly...\n\n");
#endif

  read_vega_spectra(LambdaVega, FluxVega);
#ifndef FULL_SPECTRA
  read_filters(LambdaFilter, FluxFilter);
#endif
  setup_RedshiftTab();
  read_MetalTab();

  //1st loop on the mettalicity files
  for (MetalLoop=0;MetalLoop<SSP_NMETALLICITES;MetalLoop++)
  {
#ifdef PARALLEL
    if(ThisTask == 0)
      printf("Doing Metallicity File %d of %d\n",MetalLoop+1, SSP_NMETALLICITES);
#else
    printf("Doing Metallicity File %d of %d\n",MetalLoop+1, SSP_NMETALLICITES);
#endif
    //READ FULL INPUT SPECTRA into units of erg.s^-1.Hz^-1
    read_InputSSP_spectra(LambdaInputSSP, FluxInputSSP, MetalLoop);

    //2nd Loop on redshift
    for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++)
    {
      redshift=RedshiftTab[(LastDarkMatterSnapShot+1)-snap-1];

      //3rd loop on Age
      for(AgeLoop=0;AgeLoop<SSP_NAGES;AgeLoop++)
      {
        //4th loop on Bands
        //IF FULL_SPECTRA defined a filter_number_ correspond to a wavelength on the SSP spectra
        for(filter_number_=0;filter_number_<NMAG;filter_number_++)
        {
#ifndef FULL_SPECTRA
          //ALLOCATE GRID - size of filter, binning of the spectra
          int Min_Wave_Grid=0, Max_Wave_Grid=0, Grid_Length=0;

          double *lgrid=create_grid(LambdaFilter[filter_number_][0], LambdaFilter[filter_number_][NLambdaFilter[filter_number_]-1], AgeLoop, redshift,
                                        LambdaInputSSP, &Min_Wave_Grid, &Max_Wave_Grid, &Grid_Length);

          if(Grid_Length>0)
          {
            for(i=0;i<Grid_Length;i++)
              lgrid[i]=(1+redshift)*LambdaInputSSP[AgeLoop][Min_Wave_Grid+i];

            //VEGA - interpolate spectrum on integral grid
            FluxVegaOnGrid = malloc(sizeof(double) * Grid_Length);
            interpolate(lgrid, Grid_Length, LambdaVega, NLambdaVega, FluxVega, FluxVegaOnGrid);


            /*SSP - multiply by (1+z) to go from rest to observed SSP flux*/
            FluxInputSSPOnGrid = malloc(sizeof(double) * Grid_Length);
            for (i=0;i<Grid_Length;i++)
              FluxInputSSPOnGrid[i]=(1.+redshift)*FluxInputSSP[AgeLoop][Min_Wave_Grid+i];

            //FILTERS - interpolate on integral grid
            for(i=0;i<NLambdaFilter[filter_number_];i++)
            {
              LambdaFilter_SingleFilter[i]=LambdaFilter[filter_number_][i];
              FluxFilter_SingleFilter[i]=FluxFilter[filter_number_][i];
            }
            FluxFilterOnGrid = malloc(sizeof(double) * Grid_Length);
            interpolate(lgrid, Grid_Length, LambdaFilter_SingleFilter, NLambdaFilter[filter_number_], FluxFilter_SingleFilter, FluxFilterOnGrid) ;

            /* spectrum and filters are now defined on same grid
              * CONVOLUTION: direct (configuration) space
              * simply multiply filter*spectrum it's a convolution in Fourier space */
            FluxInputSSPConv = malloc(sizeof(double) * Grid_Length);
            for(i=0;i<Grid_Length;i++)
              FluxInputSSPConv[i]=FluxInputSSPOnGrid[i]*FluxFilterOnGrid[i];

            FluxVegaConv = malloc(sizeof(double) * Grid_Length);
            for(i=0;i<Grid_Length;i++) FluxVegaConv[i]=FluxVegaOnGrid[i]*FluxFilterOnGrid[i];

            //INTEGRATE
            FluxFilterInt=integrate(FluxFilterOnGrid, Grid_Length);
            FluxVegaInt=integrate(FluxVegaConv, Grid_Length);
            FluxInputSSPInt=integrate(FluxInputSSPConv, Grid_Length);

            //Absolute Observed Frame Magnitudes
            if (FluxInputSSPInt == 0. || FluxFilterInt == 0.)
              AbsMAG=99.;
            else
            {
#ifdef AB
              AbsMAG=get_AbsAB_magnitude(FluxInputSSPInt, FluxFilterInt, redshift);
#endif
#ifdef VEGA
              AbsMAG=get_AbsAB_magnitude(FluxInputSSPInt, FluxFilterInt, redshift);
              AbsMAG=AbsMAG+2.5*(log10(FluxVegaInt)-log10(FluxFilterInt))+48.6;
              //MagABVega[filter_number_]=-2.5*(log10(FluxVegaInt)
              //                -log10(FluxFilterInt))-48.6;
#endif
            }
            free(FluxVegaOnGrid);
            free(FluxVegaConv);
            free(FluxInputSSPOnGrid);
            free(FluxInputSSPConv);
            free(FluxFilterOnGrid);
          }
          else //if Grid_Length=0 (filter outside the spectra, can happen for observed frame)
            AbsMAG=99.;

          LumTables[filter_number_][MetalLoop][snap][AgeLoop] = pow(10.,-AbsMAG/2.5);
          free(lgrid);
#else //ifdef FULL_SPECTRA
          //FULL_SPECTRA defined -> a filter_number_ corresponds to a wavelength on the SSP spectra
          FilterLambda[filter_number_]=(1+redshift)*LambdaInputSSP[AgeLoop][filter_number_];
          LumTables[filter_number_][MetalLoop][snap][AgeLoop] = (1.+redshift)*FluxInputSSP[AgeLoop][filter_number_];
#endif
        }        // end age loop
      }//end snap loop
    }//end Band loop

  }  //end loop on metallicities
  printf("\nPhotTables Computed.\n\n");
}
#endif //SPEC_PHOTABLES_ON_THE_FLY


void init_SSP_log_age_jump_index(void)
{
  double age;
  int i, idx;

  SSP_log_age_jump_factor = SSP_NJUMPTAB / (SSP_logAgeTab[SSP_NAGES - 1] - SSP_logAgeTab[1]);

  for(i = 0; i < SSP_NJUMPTAB; i++)
  {
    age = SSP_logAgeTab[1] + i / SSP_log_age_jump_factor;
    idx = 1; while(SSP_logAgeTab[idx + 1] < age) { ++idx; };
    SSP_log_age_jump_table[i] = idx;
  }
}

#endif // COMPUTE_SPECPHOT_PROPERTIES

