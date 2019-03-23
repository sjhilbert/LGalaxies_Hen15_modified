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

/** @file   model_spectro_photometric_onthefly.c
 *  @date   2009-2019
 *  @author Chiara Tonini
 *  @author Bruno Henriques
 *  @author Stefan Hilbert 
 *
 *  @brief  spectral and photometric properties of galaxies on the fly
 **/

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
  int snapshot_number_;

  for(snapshot_number_=0;snapshot_number_<(LastDarkMatterSnapShot+1);snapshot_number_++)
  {
    RedshiftTab[snapshot_number_]=ZZ[snapshot_number_];
    //for the photometry calculation do not allow negative times (it can happen for scaled cosmologies)
    if(ZZ[snapshot_number_]<0.0)
      RedshiftTab[snapshot_number_]=0.;
  }
}


void read_vega_spectra(double *LambdaVega_, double *FluxVega_)
{
  int i_;
  char buf[1000];
  FILE *fa;

  sprintf(buf, "%s/FullSEDs/VEGA_A0V_KUR_BB.SED",SpecPhotDir);
  if((fa=fopen(buf,"r"))==NULL)
  {
    char sbuf[1000];
    sprintf(sbuf, "Can't open file %s\n", buf);
    terminate(sbuf);
  }

  for (i_=0;i_<NLambdaVega;i_++)
  {
    fscanf(fa,"%lf %lf" , &LambdaVega_[i_], &FluxVega_[i_]);
    //convert to ergs.cm^-2.s^-1.AA-1
    FluxVega_[i_]*=2e-17;
    FluxVega_[i_]=FluxVega_[i_]*LambdaVega_[i_]*LambdaVega_[i_]/(SPEED_OF_LIGHT*1e8);
  }
  fclose(fa);
}


void read_filters(double LambdaFilter_[NMAG][MAX_NLambdaFilter], double FluxFilter_[NMAG][MAX_NLambdaFilter])
{
  int j_, bandn, NFilters;
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

    for(j_=0;j_<NLambdaFilter[bandn];j_++)
      fscanf(fb,"%lf %lf" ,&LambdaFilter_[bandn][j_], &FluxFilter_[bandn][j_]);

    fclose(fb);
  }
  fclose(fa);
}


void read_MetalTab()
{
  int i_,dumb_ssp_nmetallicites;
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
    sprintf(sbuf, "file `%s' not found.\n_", buf);
    terminate(sbuf);
  }

  fscanf(fa, "%d", &dumb_ssp_nmetallicites);
  if(dumb_ssp_nmetallicites != SSP_NMETALLICITES)
  {
    terminate("nmetallicites on file not equal to SSP_NMETALLICITES");
  }

  for(i_=0;i_<SSP_NMETALLICITES;i_++)
  {
    fscanf(fa, "%f", &SSP_logMetalTab[i_]);
    SSP_logMetalTab[i_]=log10(SSP_logMetalTab[i_]);
  }

  fclose(fa);
}


/** @brief Reads in the SSP full spectra for the current metallicity
 *
 * LambdaInputSSP_[NAGE][NLambdaInputSSP] &  FluxInputSSP_[NAGE][NLambdaInputSSP];
 * original units of FluxInputSSP_[i_] are (erg.s^-1.AA^-1);
 * Flux*Lambda^2/Clight*1e8 converts it to (erg.s^-1.Hz^-1) */
void read_InputSSP_spectra(double LambdaInputSSP_[SSP_NAGES][SSP_NLambda], double FluxInputSSP_[SSP_NAGES][SSP_NLambda], const int met_index_)
{
  double Dumb1, age;
  int i_, ageloop;
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

  sprintf(buf1, "%s/FullSEDs/%s_%s_FullSED_m%0.4f.dat",SpecPhotDir, SSP, SpecPhotIMF, pow(10,SSP_logMetalTab[met_index_]));

  if((fa=fopen(buf1,"r"))==NULL)
  {
    char sbuf[1000];
    sprintf(sbuf, "Can't open file %s\n", buf1);
    terminate(sbuf);
  }

  for(ageloop=0;ageloop<SSP_NAGES;ageloop++)
    for(i_=0;i_<SSP_NLambda;i_++)
    {
      LambdaInputSSP_[ageloop][i_]=0.0;
      FluxInputSSP_[ageloop][i_]=0.0;
    }

  for(ageloop=0;ageloop<SSP_NAGES;ageloop++)
  {
    for (i_=0;i_<SSP_NLambda;i_++)
    {
      fscanf(fa,"%lf %lf %lf %lf\n" , &age, &Dumb1,&LambdaInputSSP_[ageloop][i_], &FluxInputSSP_[ageloop][i_]);
      FluxInputSSP_[ageloop][i_]=1e11*FluxInputSSP_[ageloop][i_]*LambdaInputSSP_[ageloop][i_]*LambdaInputSSP_[ageloop][i_]/(SPEED_OF_LIGHT*1.e8);
    }
    if(met_index_==0) //only read age table once
    {
      if(age>0.)
        SSP_logAgeTab[ageloop]=log10(age * 1.0e-6 / UnitTime_in_Megayears * Hubble_h);
      else
        SSP_logAgeTab[ageloop]=0.;
    }
  }
  
  fclose(fa);
}

/** @brief get reference AB luminosity 
 *
 * it needs to be converted to cm since the units of
 * the original InputSSP spectra are (erg.s^-1.AA^-1) which
 * was converted to (erg.s^-1.Hz^-1) by doing Flux*Lambda^2/Clight*1e8
 * we need 3631Jy or 3631*10^-23 (= 3.631*10^-20) erg.s^-1.Hz^-1 cm-2
 * 4*pi*10pc^2*3631jy*erg.s^-1.Hz^-1
 */
static inline double
get_reference_AB_luminosity(const double FluxInputSSPInt_, const double FluxFilterInt_)
{
//it needs to be converted to cm since the units of
//the original InputSSP spectra are (erg.s^-1.AA^-1) which
//was converted to (erg.s^-1.Hz^-1) by doing Flux*Lambda^2/Clight*1e8
//we need 3631Jy or 3631*10^-23 (= 3.631*10^-20) erg.s^-1.Hz^-1 cm-2
//4*pi*10pc^2*3631jy*erg.s^-1.Hz^-1
  return FluxInputSSPInt_ / (FluxFilterInt_ * (4. * M_PI * 3.631e-20 * LENGTH_10_PC_IN_CM * LENGTH_10_PC_IN_CM));
}


/** @brief numerical integration for fluxes
 * 
 * numerical integration for fluxes with
 * Simpsons quadratures for the signal
 * 
 * erg/s/A --> erg/s/Hz --> erg/s
 */
static double integrate_flux(double *flux_, const int grid_length_)
{
  double sum[3], I[3], f[4];
  double integral=0.0;
  int i_,k;

  for(i_=0;i_<3;i_++)sum[i_]=0.0;
  for(i_=0;i_<3;i_++)I[i_]=0.0;
  for(i_=0;i_<4;i_++)f[i_]=0.0;

  for(i_=0;i_<grid_length_/2-2;i_++)
    {
      k=2*i_+1;                    //odd indexes
      f[2]=flux_[k];
      sum[1]=sum[1]+f[2];
    }
  I[1]=sum[1]*2./3.;

  for(i_=0;i_<grid_length_/2-1;i_++)
    {
      k=2*i_;
      f[3]=flux_[k] ;     //even indexes
      sum[2]=sum[2]+f[3];
    }
  I[2]=sum[2]*4./3.;

  f[0]=flux_[0];
  f[1]=flux_[grid_length_-1];
  I[0]=(f[0]+f[1])/3.;

  integral=I[0]+I[1]+I[2];

//if(grid_length_==0)
//  printf("Integral=%e\n",integral);
  return integral;
}


/** @brief find interpolation point
 * 
 * Finds interpolation point
 * the value j_ so that xx_[j_] < x_ < xx_[j_+1]
 */
static inline void 
locate(double *xx_, const int n_, const double x_, int *j_)
{
  unsigned long ju_,jm_,jl_;
  int ascnd_;

  jl_=0;
  ju_=n_+1;
  ascnd_=(xx_[n_] >= xx_[1]);
  while (ju_-jl_ > 1)
    {
      jm_=(ju_+jl_) >> 1;
      if ((x_ >= xx_[jm_]) == ascnd_)
        jl_=jm_;
      else
        ju_=jm_;
    }
  if (x_ == xx_[1]) *j_=1;
  else if(x_ == xx_[n_]) *j_=n_-1;
  else *j_=jl_;
}


/** @brief interpolates filters on integral grid */
static void interpolate_flux(double *luminosity_grid_, const int grid_length_, double *lambda_, const int n_lambda_, double *flux_, double *flux_on_grid_)
{
  int kk_=0, nn_=0, m_=2, i_;

  for(i_=0;i_<grid_length_;i_++)
  {
    if (luminosity_grid_[i_] < lambda_[0] || luminosity_grid_[i_] > lambda_[n_lambda_-1]) flux_on_grid_[i_]=0;
    //outside filter range, transmission is 0
  else
    {
      //finds where wavelenght is in respect to the grid
      locate(lambda_,n_lambda_-1,luminosity_grid_[i_],&nn_);
      kk_=min(max(nn_-(m_-1)/2,1),n_lambda_+1-m_);
      flux_on_grid_[i_]=flux_[kk_];
    }
  }
}


/** @brief creates a grid of points based on the input SSP.
 *
 * The first point on the grid is the first point for which 
 * LambdaInputSSP_>FilterWaveMin. The last point is the largest
 * wavelength for which LambdaInputSSP_<FilterWaveMax 
 */
static double* create_grid (const double wave_min_, const double wave_max_, const int age_index_, const double redshift_, double LambdaInputSSP_[SSP_NAGES][SSP_NLambda],
                                      int *min_wave_grid_, int *max_wave_grid_, int *grid_length_)
{
  int i_, min_, max_;
  double *grid;

  min_= min(wave_min_, wave_max_);
  max_= max(wave_min_, wave_max_);

  *grid_length_=0;

  //get minimum of grid
  for(i_=0;i_<SSP_NLambda;i_++)
    if((1+redshift_)*LambdaInputSSP_[age_index_][i_]>=min_)
    {
      *min_wave_grid_=i_;
      break;
    }

  for(i_=0;i_<SSP_NLambda;i_++)
    if((1+redshift_)*LambdaInputSSP_[age_index_][i_]>=min_)
    {
      *grid_length_+=1;

      //point at maximum range or out of it, set maximum
      if((1+redshift_)*LambdaInputSSP_[age_index_][i_]>=max_)
      {
        if((1+redshift_)*LambdaInputSSP_[age_index_][i_]==max_)
          *max_wave_grid_=i_;
        else //if point out of range, set max to previous
        {
          *max_wave_grid_=i_-1;
          *grid_length_-=1;
        }
        break;
      }
    }

  grid = malloc(sizeof(double) * *grid_length_);
  for(i_=0;i_<*grid_length_;i_++)
    grid[i_]=(1+redshift_)*LambdaInputSSP_[age_index_][*min_wave_grid_+i_];

  return grid;
}


/** Reads in the Full SEDs from given stellar population and filter curves and
 * computes PhotTables on the fly. Just needs the SEDs in SpecPhotDir/FullSEDs/,
 * the file with filter names and wave_lengths "FileWithFilterNames" and the
 * filter curves in SpecPhotDir/Filters/
 *
 * Developed by Chiara Tonini, adapted by Bruno Henriques, Stefan Hilbert
 **/
void setup_Spec_LumTables_onthefly(void)
{
  FilterLambda[NMAG] = 0.55;        // used by the dust model for birth clouds, the wavelength of the V-filter_number_
  
  double luminosity_;
  //FILTERS
  double LambdaFilter_[NMAG][MAX_NLambdaFilter], FluxFilter_[NMAG][MAX_NLambdaFilter];
  double *FluxFilterOnGrid_, FluxFilterInt_;
  //InputSSP spectra
  double LambdaInputSSP_[SSP_NAGES][SSP_NLambda], FluxInputSSP_[SSP_NAGES][SSP_NLambda];
  double *FluxInputSSPConv_, *FluxInputSSPOnGrid_, FluxInputSSPInt_;
#ifdef VEGA
  //VEGA
  double LambdaVega_[NLambdaVega], FluxVega_[NLambdaVega];
  double *FluxVegaConv_, *FluxVegaOnGrid_;
  double FluxVegaInt_;
#endif /* defined VEGA */

  //AUX ARRAYS for FILTERS and InputSSP
  double LambdaFilter_SingleFilter_[MAX_NLambdaFilter], FluxFilter_SingleFilter[MAX_NLambdaFilter];
  /* double LambdaInputSSP_SingleAge[SSP_NLambda];*/
  double redshift_;
  //Loops in the code
  int met_index_, age_index_, snapshot_number_, filter_number_, i_;

  if(ThisTask == 0)
    printf("\n_\nComputing PhotTables on the fly...\n_\n");
#ifdef VEGA
  read_vega_spectra(LambdaVega_, FluxVega_);
#endif /* defined VEGA */
#ifndef FULL_SPECTRA
  read_filters(LambdaFilter_, FluxFilter_);
#endif /* not defined FULL_SPECTRA */
  setup_RedshiftTab();
  read_MetalTab();

  //1st loop on the mettalicity files
  for (met_index_=0; met_index_ < SSP_NMETALLICITES; met_index_++)
  {
    if(ThisTask == 0)
      printf("Doing Metallicity File %d of %d\n",met_index_+1, SSP_NMETALLICITES);

    //READ FULL INPUT SPECTRA into units of erg.s^-1.Hz^-1
    read_InputSSP_spectra(LambdaInputSSP_, FluxInputSSP_, met_index_);

    //2nd Loop on redshift_
    for(snapshot_number_ = 0; snapshot_number_ <= LastDarkMatterSnapShot; snapshot_number_++)
    {
      redshift_ = RedshiftTab[LastDarkMatterSnapShot - snapshot_number_];

      //3rd loop on Age
      for(age_index_ = 0; age_index_ < SSP_NAGES; age_index_++)
      {
        //4th loop on Bands
        //IF FULL_SPECTRA defined a filter_number_ correspond to a wavelength on the SSP spectra
        for(filter_number_=0; filter_number_ < NMAG; filter_number_++)
        {
#ifndef FULL_SPECTRA
          //ALLOCATE GRID - size of filter, binning of the spectra
          int min_wave_grid_=0, max_wave_grid_=0, grid_length_=0;

          double *luminosity_grid_=create_grid(LambdaFilter_[filter_number_][0], LambdaFilter_[filter_number_][NLambdaFilter[filter_number_]-1], age_index_, redshift_,
                                        LambdaInputSSP_, &min_wave_grid_, &max_wave_grid_, &grid_length_);

           //if grid_length_=0 (filter outside the spectra, can happen for observed frame)
          if(grid_length_ <= 0)
          { luminosity_ = 0.; }   
          else
          {
            for(i_=0;i_<grid_length_;i_++)
              luminosity_grid_[i_]=(1+redshift_)*LambdaInputSSP_[age_index_][min_wave_grid_+i_];

#ifdef VEGA
            //VEGA - interpolate spectrum on integral grid
            FluxVegaOnGrid_ = malloc(sizeof(double) * grid_length_);
            interpolate_flux(luminosity_grid_, grid_length_, LambdaVega_, NLambdaVega, FluxVega_, FluxVegaOnGrid_);
#endif /* defined VEGA */

            /*SSP - multiply by (1+z) to go from rest to observed SSP flux_*/
            FluxInputSSPOnGrid_ = malloc(sizeof(double) * grid_length_);
            for (i_=0;i_<grid_length_;i_++)
              FluxInputSSPOnGrid_[i_]=(1.+redshift_)*FluxInputSSP_[age_index_][min_wave_grid_+i_];

            //FILTERS - interpolate on integral grid
            for(i_=0;i_<NLambdaFilter[filter_number_];i_++)
            {
              LambdaFilter_SingleFilter_[i_]=LambdaFilter_[filter_number_][i_];
              FluxFilter_SingleFilter[i_]=FluxFilter_[filter_number_][i_];
            }
            FluxFilterOnGrid_ = malloc(sizeof(double) * grid_length_);
            interpolate_flux(luminosity_grid_, grid_length_, LambdaFilter_SingleFilter_, NLambdaFilter[filter_number_], FluxFilter_SingleFilter, FluxFilterOnGrid_) ;

            /* spectrum and filters are now defined on same grid
              * CONVOLUTION: direct (configuration) space
              * simply multiply filter*spectrum it's a convolution in Fourier space */
            FluxInputSSPConv_ = malloc(sizeof(double) * grid_length_);
            for(i_=0;i_<grid_length_;i_++)
              FluxInputSSPConv_[i_]=FluxInputSSPOnGrid_[i_]*FluxFilterOnGrid_[i_];
#ifdef VEGA
            FluxVegaConv_ = malloc(sizeof(double) * grid_length_);
            for(i_=0;i_<grid_length_;i_++) FluxVegaConv_[i_]=FluxVegaOnGrid_[i_]*FluxFilterOnGrid_[i_];
#endif /* defined VEGA */
            //INTEGRATE
            FluxFilterInt_ = integrate_flux(FluxFilterOnGrid_, grid_length_);
#ifdef VEGA
            FluxVegaInt_=integrate_flux(FluxVegaConv_, grid_length_);
#endif /* defined VEGA */
            FluxInputSSPInt_ = integrate_flux(FluxInputSSPConv_, grid_length_);

            //Absolute Observed Frame Magnitudes
            if (FluxInputSSPInt_ == 0. || FluxFilterInt_ == 0.)
              luminosity_ = 0.;
            else
            {
#ifdef AB
              // AbsMAG=get_AbsAB_magnitude(FluxInputSSPInt_, FluxFilterInt_, redshift_);
             //  AbsMAG = -2.5 * log10(get_reference_AB_luminosity(FluxInputSSPInt_, FluxFilterInt_));
              luminosity_ = get_reference_AB_luminosity(FluxInputSSPInt_, FluxFilterInt_);
#endif /* defined AB */
#ifdef VEGA
              // AbsMAG=get_AbsAB_magnitude(FluxInputSSPInt_, FluxFilterInt_, redshift_);
              // AbsMAG = -2.5 * log10(get_reference_AB_luminosity(FluxInputSSPInt_, FluxFilterInt_));
              luminosity_ = get_reference_AB_luminosity(FluxInputSSPInt_, FluxFilterInt_);
             
              // AbsMAG=AbsMAG+2.5*(log10(FluxVegaInt_)-log10(FluxFilterInt_))+48.6;
              luminosity_ *= (FluxFilterInt_ * 3.631e-20) / FluxVegaInt_);

              //MagABVega[filter_number_]=-2.5*(log10(FluxVegaInt_)
              //                -log10(FluxFilterInt_))-48.6;
#endif /* defined VEGA */
            }
#ifdef VEGA
            free(FluxVegaOnGrid_);
            free(FluxVegaConv_);
#endif /* defined VEGA */
            free(FluxInputSSPOnGrid_);
            free(FluxInputSSPConv_);
            free(FluxFilterOnGrid_);
          }

          // LumTables[age_index_][met_index_][snapshot_number_][filter_number_] = pow(10.,-AbsMAG/2.5);
          LumTables[age_index_][met_index_][snapshot_number_][filter_number_] = luminosity_;
          free(luminosity_grid_);
#else /* defined FULL_SPECTRA */
          //FULL_SPECTRA defined -> a filter_number_ corresponds to a wavelength on the SSP spectra
          FilterLambda[filter_number_]=(1+redshift_)*LambdaInputSSP_[age_index_][filter_number_];
          LumTables[age_index_][met_index_][snapshot_number_][filter_number_] = (1.+redshift_)*FluxInputSSP_[age_index_][filter_number_];
#endif /* defined FULL_SPECTRA */
        }        // end age loop
      }//end snapshot_number_ loop
    }//end Band loop

  }  //end loop on metallicities
  printf("\nPhotTables Computed.\n_\n");
}
#endif /* defined SPEC_PHOTABLES_ON_THE_FLY */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
