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
#include <sys/stat.h>
#include <stddef.h>

#include "allvars.h"
#include "proto.h"

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef SPEC_PHOTABLES_ON_THE_FLY

/** @brief Setup table with redshifts to calculate observed frame magnitudes.
 *
 * The list is different from the one on the code (ZZ[]) since no negative
 * times are allowed, which can happen for scaled cosmologies. */
void setup_RedshiftTab()
{
  int snap;

  for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++)
  {
    RedshiftTab[snap]=ZZ[snap];
    //for the photometry calculation do not allow negative times (it can happen for scaled cosmologies)
    if(ZZ[snap]<0.0)
      RedshiftTab[snap]=0.;
  }
}


void read_vega_spectra(double *LambdaVega, double *FluxVega)
{
  int i;
  char buf[1000];
  FILE *fa;

  sprintf(buf, "%s/FullSEDs/VEGA_A0V_KUR_BB.SED",SpecPhotDir);
  if((fa=fopen(buf,"r"))==NULL)
  {
    char sbuf[1000];
    sprintf(sbuf, "Can't open file %s\n", buf);
    terminate(sbuf);
  }

  for (i=0;i<NLambdaVega;i++)
  {
    fscanf(fa,"%lf %lf" , &LambdaVega[i], &FluxVega[i]);
    //convert to ergs.cm^-2.s^-1.AA-1
    FluxVega[i]*=2e-17;
    FluxVega[i]=FluxVega[i]*LambdaVega[i]*LambdaVega[i]/(SPEED_OF_LIGHT*1e8);
  }
  fclose(fa);
}


void read_filters(double LambdaFilter[NMAG][MAX_NLambdaFilter], double FluxFilter[NMAG][MAX_NLambdaFilter])
{
  int j, bandn, NFilters;
  FILE *fa, *fb;
  char buf[1000], buf2[1000], FilterFile[1000], FilterName[1000];

  sprintf(buf, "%s",FileWithFilterNames);
  if((fa=fopen(buf,"r"))==NULL)
  {
    char sbuf[1000];
    sprintf(sbuf, "Can't open file %s\n", buf);
    terminate(sbuf);
  }

  fscanf(fa,"%d" ,&NFilters);
  if (NFilters != NMAG) {printf("NFilters not equal to  NMAG, line %d of read_filters.c!!! ",__LINE__);exit(0);}

  for(bandn=0;bandn<NMAG;bandn++)
  {
    fscanf(fa,"%s %f %s" ,FilterFile, &FilterLambda[bandn],FilterName);
    sprintf(buf2, "%s/Filters/%s",SpecPhotDir,FilterFile);

    if((fb=fopen(buf2,"r"))==NULL)
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file %s\n", buf2);
      terminate(sbuf);
    }

    fscanf(fb,"%d" ,&NLambdaFilter[bandn]);
    if(NLambdaFilter[bandn]>MAX_NLambdaFilter)
    {
      char sbuf[1000];
      sprintf(sbuf, "NLambdaFilter[%d]>MAX_NLambdaFilter \n", bandn);
      terminate(sbuf);
    }

    for(j=0;j<NLambdaFilter[bandn];j++)
      fscanf(fb,"%lf %lf" ,&LambdaFilter[bandn][j], &FluxFilter[bandn][j]);

    fclose(fb);
  }

  fclose(fa);

}


void read_MetalTab()
{
  int i,dumb_ssp_nmetallicites;
  FILE *fa;
  char buf[1000];
#ifdef M05
  char *SSP = {"M05"};
#endif
#ifdef CB07
  char *SSP = {"CB07"};
#endif
#ifdef BC03
  char *SSP = {"BC03"};
#endif

  sprintf(buf, "%s/FullSEDs/%s_%s_Metallicity_list.dat", SpecPhotDir, SSP, SpecPhotIMF);
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

  for(i=0;i<SSP_NMETALLICITES;i++)
  {
    fscanf(fa, "%f", &SSP_logMetalTab[i]);
    SSP_logMetalTab[i]=log10(SSP_logMetalTab[i]);
  }

  fclose(fa);
}


/** @brief Reads in the SSP full spectra for the current metallicity
 *
 * LambdaInputSSP[NAGE][NLambdaInputSSP] &  FluxInputSSP[NAGE][NLambdaInputSSP];
 * original units of FluxInputSSP[i] are (erg.s^-1.AA^-1);
 * Flux*Lambda^2/Clight*1e8 converts it to (erg.s^-1.Hz^-1) */
void read_InputSSP_spectra(double LambdaInputSSP[SSP_NAGES][SSP_NLambda], double FluxInputSSP[SSP_NAGES][SSP_NLambda], int MetalLoop)
{
  double Dumb1, age;
  int i, ageloop;
  FILE *fa;
  char buf1[1000];
#ifdef M05
  char *SSP = {"M05"};
#endif
#ifdef CB07
  char *SSP = {"CB07"};
#endif
#ifdef BC03
  char *SSP = {"BC03"};
#endif

  sprintf(buf1, "%s/FullSEDs/%s_%s_FullSED_m%0.4f.dat",SpecPhotDir, SSP, SpecPhotIMF, pow(10,SSP_logMetalTab[MetalLoop]));

  if((fa=fopen(buf1,"r"))==NULL)
  {
    char sbuf[1000];
    sprintf(sbuf, "Can't open file %s\n", buf1);
    terminate(sbuf);
  }

  for(ageloop=0;ageloop<SSP_NAGES;ageloop++)
    for(i=0;i<SSP_NLambda;i++)
    {
      LambdaInputSSP[ageloop][i]=0.0;
      FluxInputSSP[ageloop][i]=0.0;
    }

  for(ageloop=0;ageloop<SSP_NAGES;ageloop++)
  {
    for (i=0;i<SSP_NLambda;i++)
    {
      fscanf(fa,"%lf %lf %lf %lf\n" , &age, &Dumb1,&LambdaInputSSP[ageloop][i], &FluxInputSSP[ageloop][i]);
      FluxInputSSP[ageloop][i]=1e11*FluxInputSSP[ageloop][i]*LambdaInputSSP[ageloop][i]*LambdaInputSSP[ageloop][i]/(SPEED_OF_LIGHT*1.e8);
    }
    if(MetalLoop==0) //only read age table once
      if(age>0.)
        SSP_logAgeTab[ageloop]=log10(age / 1.0e6 / UnitTime_in_Megayears * Hubble_h);
      else
        SSP_logAgeTab[ageloop]=0.;
  }
  
  fclose(fa);
}


double get_AbsAB_magnitude(const double FluxInputSSPInt, const double FluxFilterInt, const double redshift)
{
//it needs to be converted to cm since the units of
//the original InputSSP spectra are (erg.s^-1.AA^-1) which
//was converted to (erg.s^-1.Hz^-1) by doing Flux*Lambda^2/Clight*1e8
//we need 3631Jy or 3631*10^-23 erg.s^-1.Hz^-1 cm-2
//4*pi*10pc^2*3631jy*erg.s^-1.Hz^-1

#ifdef APP
  const double area      = get_area(redshift);
  const double zeropoint = +48.6-2.5*area;
#else
  const double distance_cm = 1 0.0*3.08568025e18;
  const double zeropoint   = -2.5*log10(4.0*M_PI*distance_cm*distance_cm*3631.0*1.0e-23);
#endif

  const double AbsAB = -2.5*(log10(FluxInputSSPInt) -log10(FluxFilterInt)) - zeropoint;

  return AbsAB;
}


//use for apparent magnitudes
double get_area (const double redshift)
{
  double dist, area;

  //when calculating apparent magnitudes, minimum dist set to 10pc
  if (redshift<0.00000001)
  { dist=1e-5; }
  else
  { dist=lum_distance(redshift); }

//  if(dist > 0.0)
    area=log10(4.*M_PI)+2.*log10(dist*3.08568025e24);  //in cm (dl in Mpc)
//  else area=0.0;

  return area;
}


//LUMINOSITY DISTANCE
double lum_distance(const double redshift)
{
  int i, k, Npoints=1000;
  double x[1000];
  double sum[2], I[3], f[4];
  double h, integral, dl;

  for(i=0;i<2;i++)sum[i]=0.0;
  for(i=0;i<3;i++)I[i]=0.0;
  for(i=0;i<4;i++)f[i]=0.0;

  h=redshift/(Npoints-1);
  for(i=0;i<Npoints;i++) x[i]=h*(i-1);


  for (i=0;i<Npoints/2;i++)
    {
      k=2*i-1;
      f[2]=1./sqrt((1.+x[k])*(1.+x[k])*(1.+Omega*x[k])-x[k]*OmegaLambda*(2.+x[k]));
      sum[0]=sum[0]+f[2];
    }
  I[1]=sum[0]*4./3.;

  for (i=0;i<Npoints/2-1;i++)
    {
      k=2*i;
      f[3]=1./sqrt((1.+x[k])*(1.+x[k])*(1.+Omega*x[k])-x[k]*OmegaLambda*(2.+x[k]));
      sum[1]=sum[1]+f[3];
    }
  I[2]=sum[1]*2./3.;

  f[1]=1./sqrt((1.+x[0])*(1.+x[0])*(1.+Omega*x[0])-x[0]*OmegaLambda*(2.+x[0]));
  f[2]=1./sqrt((1.+x[Npoints-1])*(1.+x[Npoints-1])*(1.+Omega*x[Npoints-1])
               -x[Npoints-1]*OmegaLambda*(2.+x[Npoints-1]));

  I[0]=(f[0]+f[1])/3.;

  integral=h*(I[0]+I[1]+I[2]);

  dl=integral/1000.0;    //!Mpc

  dl*=(1.+redshift)*(SPEED_OF_LIGHT/100.)/(Hubble_h*100.);   //in Mpc

  return dl;
}


/** @brief numerical integration for fluxes
 * 
 * numerical integration for fluxes with
 * Simpsons quadratures for the signal
 * 
 * erg/s/A --> erg/s/Hz --> erg/s
 */
double integrate_flux(double *flux, const int Grid_Length)
{
  double sum[3], I[3], f[4];
  double integral=0.0;
  int i,k;

  for(i=0;i<3;i++)sum[i]=0.0;
  for(i=0;i<3;i++)I[i]=0.0;
  for(i=0;i<4;i++)f[i]=0.0;

  for(i=0;i<Grid_Length/2-2;i++)
    {
      k=2*i+1;                    //odd indexes
      f[2]=flux[k];
      sum[1]=sum[1]+f[2];
    }
  I[1]=sum[1]*2./3.;

  for(i=0;i<Grid_Length/2-1;i++)
    {
      k=2*i;
      f[3]=flux[k] ;     //even indexes
      sum[2]=sum[2]+f[3];
    }
  I[2]=sum[2]*4./3.;

  f[0]=flux[0];
  f[1]=flux[Grid_Length-1];
  I[0]=(f[0]+f[1])/3.;

  integral=I[0]+I[1]+I[2];

//if(Grid_Length==0)
//  printf("Integral=%e\n",integral);
  return integral;
}


/** @brief find interpolation point
 * 
 * Finds interpolation point
 * the value j so that xx[j] < x < xx[jj+1]
 * 
 * @note  MATH MISC - PROBABLY SHOULD GO INTO SEPARATE FILE
 */
void locate(double *xx, const int n, const double x, int *j)
{
  unsigned long ju,jm,jl;
  int ascnd;

  jl=0;
  ju=n+1;
  ascnd=(xx[n] >= xx[1]);

  while (ju-jl > 1)
    {
      jm=(ju+jl) >> 1;
      if ((x >= xx[jm]) == ascnd)
        jl=jm;
      else
        ju=jm;
    }

  if (x == xx[1]) *j=1;
  else if(x == xx[n]) *j=n-1;
  else *j=jl;

}


/** @brief interpolates filters on integral grid */
void interpolate_flux(double *lgrid, const int Grid_Length, double *lambda, const int nlambda, double *flux, double *FluxOnGrid)
{
  int kk=0, nn=0, m=2, i;

  for(i=0;i<Grid_Length;i++)
  {
    if (lgrid[i] < lambda[0] || lgrid[i] > lambda[nlambda-1]) FluxOnGrid[i]=0;
    //outside filter range, transmission is 0
  else
    {
      //finds where wavelenght is in respect to the grid
      locate(lambda,nlambda-1,lgrid[i],&nn);
      kk=min(max(nn-(m-1)/2,1),nlambda+1-m);
      FluxOnGrid[i]=flux[kk];
    }
  }
}


/** @brief creates a grid of points based on the input SSP.
 *
 * The first point on the grid is the first point for which 
 * LambdaInputSSP>FilterWaveMin. The last point is the largest
 * wavelength for which LambdaInputSSP<FilterWaveMax 
 */
double* create_grid (const double WaveMin, const double WaveMax, const int AgeLoop, const double redshift, double LambdaInputSSP[SSP_NAGES][SSP_NLambda],
                                      int *Min_Wave_Grid, int *Max_Wave_Grid, int *Grid_Length)
{
  double x0, x1, h;
  int i, min, max;
  double *grid;

  min=min(WaveMin, WaveMax);
  max=max(WaveMin, WaveMax);

  *Grid_Length=0;

  //get minimum of grid
  for(i=0;i<SSP_NLambda;i++)
    if((1+redshift)*LambdaInputSSP[AgeLoop][i]>=min)
      {
        *Min_Wave_Grid=i;
        break;
      }

  for(i=0;i<SSP_NLambda;i++)
    if((1+redshift)*LambdaInputSSP[AgeLoop][i]>=min)
      {
        *Grid_Length+=1;

        //point at maximum range or out of it, set maximum
        if((1+redshift)*LambdaInputSSP[AgeLoop][i]>=max)
          {
            if((1+redshift)*LambdaInputSSP[AgeLoop][i]==max)
              *Max_Wave_Grid=i;
            else //if point out of range, set max to previous
              {
                *Max_Wave_Grid=i-1;
                *Grid_Length-=1;
              }
            break;
          }
      }

  grid = malloc(sizeof(double) * *Grid_Length);
  for(i=0;i<*Grid_Length;i++)
    grid[i]=(1+redshift)*LambdaInputSSP[AgeLoop][*Min_Wave_Grid+i];

  return grid;
}


/** Reads in the Full SEDs from a given stellar population and filter curves and
 * computes PhotTables on the fly. Just needs the SEDs in SpecPhotDir/FullSEDs/,
 * the file with filter names and wave_lengths "FileWithFilterNames" and the
 * filter curves in SpecPhotDir/Filters/
 *
 * Developed by Chiara Tonini, adapted by Bruno Henriques
 *
 * 
 *
 * */
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
            interpolate_flux(lgrid, Grid_Length, LambdaVega, NLambdaVega, FluxVega, FluxVegaOnGrid);


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
            interpolate_flux(lgrid, Grid_Length, LambdaFilter_SingleFilter, NLambdaFilter[filter_number_], FluxFilter_SingleFilter, FluxFilterOnGrid) ;

            /* spectrum and filters are now defined on same grid
              * CONVOLUTION: direct (configuration) space
              * simply multiply filter*spectrum it's a convolution in Fourier space */
            FluxInputSSPConv = malloc(sizeof(double) * Grid_Length);
            for(i=0;i<Grid_Length;i++)
              FluxInputSSPConv[i]=FluxInputSSPOnGrid[i]*FluxFilterOnGrid[i];

            FluxVegaConv = malloc(sizeof(double) * Grid_Length);
            for(i=0;i<Grid_Length;i++) FluxVegaConv[i]=FluxVegaOnGrid[i]*FluxFilterOnGrid[i];

            //INTEGRATE
            FluxFilterInt=integrate_flux(FluxFilterOnGrid, Grid_Length);
            FluxVegaInt=integrate_flux(FluxVegaConv, Grid_Length);
            FluxInputSSPInt=integrate_flux(FluxInputSSPConv, Grid_Length);

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

          LumTables[AgeLoop][MetalLoop][snap][filter_number_] = pow(10.,-AbsMAG/2.5);
          free(lgrid);
#else //ifdef FULL_SPECTRA
          //FULL_SPECTRA defined -> a filter_number_ corresponds to a wavelength on the SSP spectra
          FilterLambda[filter_number_]=(1+redshift)*LambdaInputSSP[AgeLoop][filter_number_];
          LumTables[AgeLoop][MetalLoop][snap][filter_number_] = (1.+redshift)*FluxInputSSP[AgeLoop][filter_number_];
#endif
        }        // end age loop
      }//end snap loop
    }//end Band loop

  }  //end loop on metallicities
  printf("\nPhotTables Computed.\n\n");
}
#endif /* defined SPEC_PHOTABLES_ON_THE_FLY */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
