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
 *  You should have received a_ copy of the GNU General Public License
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"
#ifdef HALOMODEL
#include "mcmc_halomodel.h"
#endif


//////////////
//LIKELIHOOD//
//////////////


/** @file mcmc_likelihood.c
 *  @brief This function computes the likelihhod of SAM galaxies in respect
 *         to a_ number of observational constraints (Henriques2009)
 *
 *  This function computes the likelihood between the SAM and a_ set of observational
 *  data sets. The galaxy properties for the SAM are computed using a_ given set
 *  of input parameters and in the end of each run this function is called to get
 *  the likelihood of properties with the new set of parameters and compare it with
 *  the previous run.
 *
 * Three different tests are available: Chi^2, Maximum Likelihhod Method and Binomial
 * probability. These are used to compare the model with the different observational
 * properties read at read_observations and stored in struct MCMC_Obs.
 *
 * prob_[0:13] are computed and the likelihood is given by a_ desired
 * combination of these at the exit of the get_likelihho() function.
 *
 * 1st Test - Chisquare test for the SMF
 * MCMC_Obs[].Obs[0] & MCMC_Obs[].Error[0] - observations from baldry2008
 * MCMC_GAL[].StellarMass - Masses from SAM
 *
 * 8th Test - MLM test for the colours
 * MCMC_Obs[].Obs[1] & MCMC_Obs[].Error[1] - observations from baldry2004
 * binredfraction - fraction of red galaxies in each mass bin
 * fract_error_ - error in the SAM red fraction (assumed to have the same value in all bins
 * - 0.05 to account for the wiggles in the 2 color luminosity function, due to the fact
 * that only one file is being used )
 *
 * 2nd Test - Chisquare test for the K-Band LF
 * MCMC_Obs[].Obs[2] & MCMC_Obs[].Error[2] - observations from Cole2003 + Bell2003+ Jones2006
 * MCMC_GAL[].MagK - SAM MagK
 *
 *
 * 10th Test - Binomial test for the Black Hole-Bulge Mass relation
 * MCMC_Obs[].ObsUp[0] & MCMC_Obs[].ObsDown[0] - observations for the BHBM from Haring & Rix 2004
 * MCMC_GAL[].Bulge & MCMC_GAL[].BlackHoleMass - masses from SAM
 * bin_black_hole_up_ - numbers of galaxies in the two upper bins on the bulge-blackhole mass relation
 * bin_black_hole_down_ - numbers of galaxies in the two lower bins on the bulge-blackhole mass relation */

#ifdef MCMC
double get_likelihood()
{
  //Variables for all the LF/SM function like tests
  double *binned_sam_data_, *sam_data_;
  //variables for the bulge-blackhole mass test using binomial
  double bin_black_hole_up_  [2] = {0.0, 0.0};
  double bin_black_hole_down_[2] = {0.0, 0.0};
  double final_probability_, redshift_probability_, current_probability_;
  // double prob_SMF_z0;
  int observation_number_, galaxy_number_, output_number_;
  double color_U_minus_V_, color_V_minus_J_;
  //BEST FIT TO MODEL CUT
#ifndef HALOMODEL
  const double offset_color_cut_[NOUT]={0.00, 1.085, 1.1, 1.0, 1.15}; //not used at z=0
  const double slope_color_cut_[NOUT]={0.00, 0.5, 0.48, 0.38, 0.18};
#else
  const double offset_color_cut_[NOUT]={0.00}; //not used at z=0
  const double slope_color_cut_[NOUT]={0.00};
#endif
  FILE *file_;
  char file_name_[1000];

  /* Bin sam_data_ into binned_sam_data_ according to observational constraints.
   * The first argument of the bin functions and of the likelihood
   * function (Chi^2 ot MLM) indicates the observational data set to use.
   * ALL Luminosity Function to be compared in units of M-5logh and phi(Mpc-3h-3mag-1) */

  final_probability_ = 1.;
  for(output_number_=0;output_number_<NOUT;output_number_++)
  {
    printf("output_number = %d\n",output_number_);
    redshift_probability_=1.;

    if((sam_data_ = malloc(sizeof(double) * TotMCMCGals[output_number_])) == NULL)
      terminate("get_likelihood: malloc for TotMCMCGals failed");

#ifdef HALOMODEL
    correct_for_correlation(output_number_);
#endif

    for(observation_number_=0;observation_number_<MCMCNConstraints;observation_number_++)
    {
      if(MCMC_Obs[observation_number_].ObsTest_Switch_z[output_number_]==1)
      {
        if((binned_sam_data_ = malloc(sizeof(double) * Nbins[output_number_][observation_number_])) == NULL)
          terminate("get_likelihood: malloc for binned_sam_data_ failed.");
        
        /* default prob_. for unreckognized test: */
        current_probability_ = 1;
        // current_probability_ = 0;

/******************************
**     Chi_Sq TESTS         **
******************************/
        if(strcmp(MCMC_Obs[observation_number_].TestType,"chi_sq")==0)
        {
          //bin all the galaxies into binned_sam_data_ for properties that just require normal histograms
          for(galaxy_number_ = 0; galaxy_number_ < TotMCMCGals[output_number_]; galaxy_number_++)
          {
            sam_data_[galaxy_number_] = 0.; //initialize
            if(strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunction")==0)
            { sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].StellarMass; }
            else if(strcmp(MCMC_Obs[observation_number_].Name,"KBandLF")==0)
            { sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].MagK; }
            else if(strcmp(MCMC_Obs[observation_number_].Name,"BBandLF")==0)
            {
              if((double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.)<0.2) //Bj band at z=0
                sam_data_[galaxy_number_]=MCMC_GAL[output_number_][galaxy_number_].MagB-0.28*(MCMC_GAL[output_number_][galaxy_number_].MagB-MCMC_GAL[output_number_][galaxy_number_].MagV);
              else //BBand at all other z
                sam_data_[galaxy_number_]=MCMC_GAL[output_number_][galaxy_number_].MagB;
            }
            else if(strcmp(MCMC_Obs[observation_number_].Name,"uBandLF")==0)
            { sam_data_[galaxy_number_] = MCMC_GAL[output_number_][observation_number_].Magu; }
            else if(strcmp(MCMC_Obs[observation_number_].Name,"ColdGasMassFunction")==0)
            {  sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].ColdGas*0.54; }

            //Stellar Mass Function of Passive Galaxies
            else if(strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionPassive")==0)
            {
              if( MCMC_GAL[output_number_][galaxy_number_].Sfr* 1.e9 /pow(10.,MCMC_GAL[output_number_][galaxy_number_].StellarMass) < 0.01 )
                sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].StellarMass;
            }
            //Stellar Mass Function of Active Galaxies
            else if(strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionActive")==0)
            {
              if( MCMC_GAL[output_number_][galaxy_number_].Sfr* 1.e9 /pow(10.,MCMC_GAL[output_number_][galaxy_number_].StellarMass) > 0.3 )
                sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].StellarMass;
            }
            //Stellar Mass Function of Red Galaxies
            //original cut in bladry 2004 (2.06-0.244*tanh((MCMC_GAL[output_number_][galaxy_number_].Magr+20.07)/1.09))
            else if(strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionRed")==0)
            {
              if((double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.)<0.2) //z=0
              {
                if( (MCMC_GAL[output_number_][galaxy_number_].Magu-MCMC_GAL[output_number_][galaxy_number_].Magr) > (1.9-0.244*tanh((MCMC_GAL[output_number_][galaxy_number_].Magr+20.07)/1.09)))
                  sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].StellarMass;
              }
              else //z>0
              {
                color_U_minus_V_=(MCMC_GAL[output_number_][galaxy_number_].MagU-MCMC_GAL[output_number_][galaxy_number_].MagV);
                color_V_minus_J_=(MCMC_GAL[output_number_][galaxy_number_].MagV-MCMC_GAL[output_number_][galaxy_number_].MagJ);
                if( (color_V_minus_J_ < (1.3-offset_color_cut_[output_number_])/slope_color_cut_[output_number_] && color_U_minus_V_ > 1.3) ||
                    (color_V_minus_J_ > (1.3-offset_color_cut_[output_number_])/slope_color_cut_[output_number_] && color_U_minus_V_ > color_V_minus_J_*slope_color_cut_[output_number_]+offset_color_cut_[output_number_]) )
                  sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].StellarMass;
              }
            }

            //Stellar Mass Function of Blue Galaxies
            //original cut in bladry 2004 (2.06-0.244*tanh((MCMC_GAL[output_number_][galaxy_number_].Magr+20.07)/1.09))
            else if(strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionBlue")==0)
            {
              if((double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.)<0.2) //z=0
              {
                if( (MCMC_GAL[output_number_][galaxy_number_].Magu-MCMC_GAL[output_number_][galaxy_number_].Magr) < (1.9-0.244*tanh((MCMC_GAL[output_number_][galaxy_number_].Magr+20.07)/1.09)))
                  sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].StellarMass;
              }
              else //z>0
              {
                color_U_minus_V_=(MCMC_GAL[output_number_][galaxy_number_].MagU-MCMC_GAL[output_number_][galaxy_number_].MagV);
                color_V_minus_J_=(MCMC_GAL[output_number_][galaxy_number_].MagV-MCMC_GAL[output_number_][galaxy_number_].MagJ);

                if( (color_V_minus_J_ < (1.3-offset_color_cut_[output_number_])/slope_color_cut_[output_number_] && color_U_minus_V_ < 1.3) ||
                    (color_V_minus_J_ > (1.3-offset_color_cut_[output_number_])/slope_color_cut_[output_number_] && color_U_minus_V_ < color_V_minus_J_*slope_color_cut_[output_number_]+offset_color_cut_[output_number_]) )
                  sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].StellarMass;
              }
            }

            //SFRF
            else if(strcmp(MCMC_Obs[observation_number_].Name,"SFRF")==0)
            { sam_data_[galaxy_number_] = MCMC_GAL[output_number_][galaxy_number_].Sfr; }

          }//end loop on number of galaxies
          bin_function(output_number_, observation_number_, sam_data_, binned_sam_data_);

          //SFRD only has one bin
          if(strcmp(MCMC_Obs[observation_number_].Name,"SFRD")==0)
          {
            binned_sam_data_[0]=0;
            for(galaxy_number_ = 0; galaxy_number_ < TotMCMCGals[output_number_]; galaxy_number_++)
              binned_sam_data_[0]+=MCMC_GAL[output_number_][galaxy_number_].Sfr*MCMC_GAL[output_number_][galaxy_number_].Weight;
          }

#ifdef HALOMODEL
          //Correlation Function - requires more than just binning
          else if(strncmp(MCMC_Obs[observation_number_].Name,"Clustering_MassBins_8.77_9.27",28)==0)
            compute_correlation_func(output_number_, observation_number_,  8.77,  9.27, binned_sam_data_);
          else if(strncmp(MCMC_Obs[observation_number_].Name,"Clustering_MassBins_9.27_9.77",28)==0)
            compute_correlation_func(output_number_, observation_number_,  9.27,  9.77, binned_sam_data_);
          else if(strncmp(MCMC_Obs[observation_number_].Name,"Clustering_MassBins_9.77_10.27",28)==0)
            compute_correlation_func(output_number_, observation_number_,  9.77, 10.27, binned_sam_data_);
          else if(strncmp(MCMC_Obs[observation_number_].Name,"Clustering_MassBins_10.27_10.77",28)==0)
            compute_correlation_func(output_number_, observation_number_, 10.27, 10.77, binned_sam_data_);
          else if(strncmp(MCMC_Obs[observation_number_].Name,"Clustering_MassBins_10.77_11.27",28)==0)
            compute_correlation_func(output_number_, observation_number_, 10.77, 11.27, binned_sam_data_);
          else if(strncmp(MCMC_Obs[observation_number_].Name,"Clustering_MassBins_11.27_11.47",28)==0)
            compute_correlation_func(output_number_, observation_number_, 11.27, 11.77, binned_sam_data_);
#endif

          current_probability_=chi_square_probability(output_number_, observation_number_, binned_sam_data_);

#ifndef PARALLEL
          //print comparison with observations into file (over-write, only makes sense if not PARALLEL)
          sprintf(file_name_, "%s/mcmc_plus_obs%dz%1.2f.txt",OutputDir,observation_number_,(double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.));
          if(!(file_ = fopen(file_name_, "w_")))
          {
            char error_message_[1000];
            sprintf(error_message_, "can't_ open file `%s'\n", file_name_);
            terminate(error_message_);
          }
          int ii;
          for(ii = 0; ii < Nbins[output_number_][observation_number_]; ii++)
            fprintf(file_, "%g %g %g %g\n",
                    MCMC_Obs[observation_number_].Bin_low[output_number_][ii]+(MCMC_Obs[observation_number_].Bin_high[output_number_][ii]-MCMC_Obs[observation_number_].Bin_low[output_number_][ii])/2.,
                    MCMC_Obs[observation_number_].Obs[output_number_][ii], MCMC_Obs[observation_number_].Error[output_number_][ii], binned_sam_data_[ii]);
          fclose(file_);
#endif

        }//end chi_sq tests


/******************************
** Maximum Likelihood TESTS **
******************************/
        else if(strcmp(MCMC_Obs[observation_number_].TestType,"maxlike")==0)
        {
          if(strcmp(MCMC_Obs[observation_number_].Name,"RedFraction")==0)
          {
            bin_red_fraction(output_number_, observation_number_, binned_sam_data_);
            current_probability_=maximum_likelihood_probability(output_number_, observation_number_, binned_sam_data_);
          }

          if(strcmp(MCMC_Obs[observation_number_].Name,"ColdGasFractionvsStellarMass")==0)
          {
            bin_ColdGasFractionvsStellarMass(output_number_, observation_number_, binned_sam_data_);
            current_probability_=maximum_likelihood_probability(output_number_, observation_number_, binned_sam_data_);
          }

          if(strcmp(MCMC_Obs[observation_number_].Name,"PassiveFraction")==0)
          {
            bin_passive_fraction(output_number_, observation_number_, binned_sam_data_);
            current_probability_=maximum_likelihood_probability(output_number_, observation_number_, binned_sam_data_);
          }

          if(strcmp(MCMC_Obs[observation_number_].Name,"BulgeFraction")==0)
          {
            bin_bulge_fraction(output_number_, observation_number_, binned_sam_data_);
            current_probability_=maximum_likelihood_probability(output_number_, observation_number_, binned_sam_data_);
          }
        }

/******************************
**    Binomial TESTS        **
******************************/
        else if(strcmp(MCMC_Obs[observation_number_].TestType,"binomial")==0)
        {
          if(strcmp(MCMC_Obs[observation_number_].Name,"BlackHoleBulgeMass")==0)
          {
            bin_bhbm(output_number_, bin_black_hole_up_, bin_black_hole_down_);
            current_probability_=binomial_probability(output_number_, observation_number_, bin_black_hole_up_, bin_black_hole_down_);
          }
        }
//END BINOMIAL TESTS

        free(binned_sam_data_);

        redshift_probability_ *= current_probability_;
        printf("prob[output_number = %d, observation_number = %d] = %0.5e\n",output_number_, observation_number_, current_probability_);
        //write likelihood for each constraint at each MCMC step
        fprintf(FILE_MCMC_PredictionsPerStep[output_number_][observation_number_], " %0.5e", current_probability_);

      }//end if(MCMC_Obs[observation_number_].ObsTest_Switch_z[output_number_]==1)
    }//end loop on MCMCNConstraints

    printf("prob[output_number = %d] = %0.5e\n",output_number_, redshift_probability_);

    final_probability_ *= redshift_probability_;
    free(sam_data_);
#ifdef HALOMODEL
    free(MCMC_FOF2);
    free(HashTable);
#endif
  }//end loop on output_number_s


  printf("final prob_=%0.5e\nlog_like=%0.8g\n",final_probability_, -log10(final_probability_));

  //write total likelihood into separate files and into
  //the end of the line on the file with the comparison to each constraint
  fprintf(FILE_MCMC_LIKELIHOOD,"%0.8g\n",-log10(final_probability_));
  fflush(FILE_MCMC_LIKELIHOOD);
  for(observation_number_=0;observation_number_<MCMCNConstraints;observation_number_++)
    for(output_number_=0;output_number_<NOUT;output_number_++)
      if(MCMC_Obs[observation_number_].ObsTest_Switch_z[output_number_]==1)
      {
        fprintf(FILE_MCMC_PredictionsPerStep[output_number_][observation_number_], " %0.8g\n", -log10(final_probability_));
        fflush(FILE_MCMC_PredictionsPerStep[output_number_][observation_number_]);
      }

  return final_probability_;
}


/**@brief Bin Luminosity or Stellar Mass Function*/
void bin_function(const int output_number_, const int observation_number_, const double *sam_data_, double *binned_sam_data_)
{
  const int N_samples_ = 100;
  int bin_number_, sample_number_, galaxy_number_;
  // FILE *file_;
  
  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  { binned_sam_data_[bin_number_] = 0.; }

  //because of the convolution with a_ random error, when masses are involved
  //we would get a_ different output for the same model parameters. Therefore
  //the mass function is computed for N_samples_ with convolved with different
  //random errors
  if(strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunction"       )==0 ||
     strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionPassive")==0 ||
     strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionActive" )==0 ||
     strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionRed"    )==0 ||
     strcmp(MCMC_Obs[observation_number_].Name,"StellarMassFunctionBlue"   )==0)
  {
    for(galaxy_number_=0;galaxy_number_<TotMCMCGals[output_number_];galaxy_number_++)
    {
      for(sample_number_ = 0; sample_number_ < N_samples_; sample_number_++)
      {  
        const double obs_sam_data_ = sam_data_[galaxy_number_] + gsl_ran_ugaussian(MCMC_rng)*AddedErrOnMass*(1+MCMCConstraintsZZ[output_number_]);
        for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
        {  
          if(MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_] <= obs_sam_data_ && obs_sam_data_ < MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_])
          { binned_sam_data_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight; }
        }
      }
    }
    for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
    { binned_sam_data_[bin_number_] /= ((MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_] - MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]) * N_samples_); } 
  }
  else
  {
    for(galaxy_number_ = 0; galaxy_number_ < TotMCMCGals[output_number_]; galaxy_number_++)
    {
      for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
      {
        if(MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_] <= sam_data_[galaxy_number_] && sam_data_[galaxy_number_] < MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_])
        { binned_sam_data_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight; }
      }
    } 
    for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
    { binned_sam_data_[bin_number_] /= (MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_] - MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]); } 
  }

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    
    fprintf(FILE_MCMC_PredictionsPerStep[output_number_][observation_number_], " %0.5e", binned_sam_data_[bin_number_]);
//#ifndef PARALLEL
  /*   fprintf(file_, "%g %g %g %g\n", MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]+(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2.,
            MCMC_Obs[observation_number_].Obs[output_number_][bin_number_], MCMC_Obs[observation_number_].Error[output_number_][bin_number_], binned_sam_data_[bin_number_]);*/
//#endif
  }
}


/**@brief Bin fraction of red galaxies.
 *    The separation is done based on g-r_ color
 *    at z=0 and and U-V vs V-J at higher z
 *
 * @note changed statistical treatment of averaging over stellar mass errors:
 *       changed from: red_fraction = <red/(red + blue)> 
 *       changed   to: red_fraction = <red>/(<red> + <blue>)
 */
void bin_red_fraction(const int output_number_, const int observation_number_, double *binned_red_fraction_)
{
  const int N_samples_ = 100;
  
  int bin_number_, sample_number_, galaxy_number_;
  double color_, color_U_minus_V_, color_V_minus_J_;
  bool is_red_galaxy_;
  
  double red_[MCMCMaxObsBins];
  double all_[MCMCMaxObsBins];

  //OBSERVATIONAL CUT
  //double offset_color_cut_[output_number_]={0.00, 0.69, 0.59, 0.59, 0.59}; //not used at z=0
  //double slope_color_cut_[output_number_]={0.00, 0.88, 0.88, 0.88, 0.88};
  //BEST FIT TO MODEL CUT
#ifndef HALOMODEL
  const double offset_color_cut_[NOUT]={0.00, 1.085, 1.1, 1.0, 1.15}; //not used at z=0
  const double slope_color_cut_[NOUT]={0.00, 0.5, 0.48, 0.38, 0.18};
#else
  const double offset_color_cut_[NOUT]={0.00}; //not used at z=0
  const double slope_color_cut_[NOUT]={0.00};
#endif

  FILE *file_;
  char file_name_[1000];

#ifndef PARALLEL
  sprintf(file_name_, "%s/mcmc_plus_obs%dz%1.2f.txt",OutputDir,observation_number_,(double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.));
  if(!(file_ = fopen(file_name_, "w_")))
  {
    char error_message_[1000];
    sprintf(error_message_, "can't_ open file `%s'\n", file_name_);
    terminate(error_message_);
  }
#endif

  const bool is_low_z_ = (double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.) < 0.2;

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    red_[bin_number_] = 0.;
    all_[bin_number_] = 0.;
  }

  for(galaxy_number_ = 0; galaxy_number_ < TotMCMCGals[output_number_]; galaxy_number_++)
  {
    //color g-r_ cut for z=0  - baldry 2004
    //original cut in bladry 2004 (2.06-0.244*tanh((MCMC_GAL[output_number_][j_].Magr+20.07)/1.09))
//       if((double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.)<0.2)
    if(is_low_z_)
    {
      color_         = MCMC_GAL[output_number_][galaxy_number_].Magu-MCMC_GAL[output_number_][galaxy_number_].Magr;
      is_red_galaxy_ = (color_ > (1.9-0.244* tanh((MCMC_GAL[output_number_][galaxy_number_].Magr+20.07)*(1./1.09))));
    }
    //U-V vs V-J at higher redshift
    else
    {
      color_U_minus_V_ = (MCMC_GAL[output_number_][galaxy_number_].MagU-MCMC_GAL[output_number_][galaxy_number_].MagV);
      color_V_minus_J_ = (MCMC_GAL[output_number_][galaxy_number_].MagV-MCMC_GAL[output_number_][galaxy_number_].MagJ);
      is_red_galaxy_   = ( (color_V_minus_J_ < (1.3-offset_color_cut_[output_number_])/slope_color_cut_[output_number_] && color_U_minus_V_ > 1.3 ) ||
                           (color_V_minus_J_ > (1.3-offset_color_cut_[output_number_])/slope_color_cut_[output_number_] && color_U_minus_V_ > color_V_minus_J_*slope_color_cut_[output_number_]+offset_color_cut_[output_number_]) );
    }

    for(sample_number_ = 0; sample_number_ < N_samples_; sample_number_++)
    {
      const double obs_stellar_mass_ = MCMC_GAL[output_number_][galaxy_number_].StellarMass + gsl_ran_ugaussian(MCMC_rng)*0.08*(1+MCMCConstraintsZZ[output_number_]);
      for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
      {
        if(MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_] <= obs_stellar_mass_ && obs_stellar_mass_ < MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_])
        {
          all_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight; 
          if(is_red_galaxy_)
          { red_ [bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight; } 
        }
      }
    }              
  }
  
  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    binned_red_fraction_[bin_number_] = (all_[bin_number_] > 0.) ? red_[bin_number_] / all_[bin_number_] : 0.; 

    fprintf(FILE_MCMC_PredictionsPerStep[output_number_][observation_number_], " %0.5e", binned_red_fraction_[bin_number_]);

#ifndef PARALLEL
    fprintf(file_, "%g %g %g %g\n",  MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]+(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2., MCMC_Obs[observation_number_].Obs[output_number_][bin_number_],
                                                MCMC_Obs[observation_number_].Error[output_number_][bin_number_], binned_red_fraction_[bin_number_]);
#endif
  }
#ifndef PARALLEL
  fclose(file_);
#endif
}


/**@brief Bin fraction of passive galaxies.
 *    The separation is done based on sfr
 *
 * @note no stellar mass error?
 */
void bin_passive_fraction(const int output_number_, const int observation_number_, double *binned_passive_fraction_)
{
  int bin_number_, galaxy_number_;
  
  double passive_[MCMCMaxObsBins];
  double all_    [MCMCMaxObsBins];

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    passive_[bin_number_] = 0.;
    all_    [bin_number_] = 0.;
  }

  for(galaxy_number_ = 0 ; galaxy_number_ < TotMCMCGals[output_number_]; galaxy_number_++)
  {
    const bool is_passive_galaxy_ = (MCMC_GAL[output_number_][galaxy_number_].Sfr* 1.e9 /pow(10.,MCMC_GAL[output_number_][galaxy_number_].StellarMass) < 0.01);

    for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
    {
      if(MCMC_GAL[output_number_][galaxy_number_].StellarMass >= MCMC_Obs[observation_number_].Bin_low [output_number_][bin_number_] &&
         MCMC_GAL[output_number_][galaxy_number_].StellarMass <= MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_])
      {
        all_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight;
        if(is_passive_galaxy_)        
        { passive_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight; }
      }
    }
  }

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {  
    binned_passive_fraction_[bin_number_] = (all_[bin_number_] > 0.) ? passive_[bin_number_]  / all_[bin_number_] : 0.;
 
    //printf("mass=%f obs_ =%f samcolors=%f\n",(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2., MCMC_Obs[observation_number_].Obs[output_number_][bin_number_], binned_passive__fraction_[bin_number_]);
  }
}


/**@brief Bin fraction of bulge galaxies.*/
void bin_bulge_fraction(const int output_number_, const int observation_number_, double *binned_bulge_fraction_)
{
  int bin_number_, galaxy_number_;
  
  double bulge_[MCMCMaxObsBins];
  double all_  [MCMCMaxObsBins];

  FILE *file_;
  char file_name_[1000];

#ifndef PARALLEL
  sprintf(file_name_, "%s/mcmc_plus_obs%dz%1.2f.txt",OutputDir,observation_number_,(double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.));
  if(!(file_ = fopen(file_name_, "w_")))
  {
    char error_message_[1000];
    sprintf(error_message_, "can't_ open file `%s'\n", file_name_);
    terminate(error_message_);
  }
#endif

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    bulge_[bin_number_] = 0.;
    all_  [bin_number_] = 0.;
  }

  for(galaxy_number_ = 0; galaxy_number_ < TotMCMCGals[output_number_]; galaxy_number_++)
  {
    //log10(0.7)=-0.154902
    const bool is_bulge_galaxy_ = ((MCMC_GAL[output_number_][galaxy_number_].BulgeMass-MCMC_GAL[output_number_][galaxy_number_].StellarMass) > -0.154902);
    for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
    {
      if(MCMC_GAL[output_number_][galaxy_number_].StellarMass-2.*log10(Hubble_h) >= MCMC_Obs[observation_number_].Bin_low [output_number_][bin_number_] &&
         MCMC_GAL[output_number_][galaxy_number_].StellarMass-2.*log10(Hubble_h) <= MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_])
      {
        all_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight; 
        if(is_bulge_galaxy_)           
        { bulge_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight; }    
      }
    }
  }
  
  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    binned_bulge_fraction_[bin_number_] = (all_[bin_number_] > 0.) ? bulge_[bin_number_]  / all_[bin_number_] : 0.;
 
   //printf("mass=%f obs_ =%f samcolors=%f\n",(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2., MCMC_Obs[observation_number_].Obs[output_number_][bin_number_], binpassivefraction[bin_number_]);

#ifndef PARALLEL
    fprintf(file_, "%g %g %g %g\n",  MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]+(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2., MCMC_Obs[observation_number_].Obs[output_number_][bin_number_],
                                                MCMC_Obs[observation_number_].Error[output_number_][bin_number_], binned_bulge_fraction_[bin_number_]);
#endif
  }
#ifndef PARALLEL
  fclose(file_);
#endif
}


/**@brief Bin fraction of bulge galaxies.
*
* @note no stellar mass error?
*/
void bin_ColdGasFractionvsStellarMass(const int output_number_, const int observation_number_, double *binned_gas_fraction_)
{
  int bin_number_, galaxy_number_;
  
  double gas_[MCMCMaxObsBins];
  double all_[MCMCMaxObsBins];

  FILE *file_;
  char file_name_[1000];

#ifndef PARALLEL
  sprintf(file_name_, "%s/mcmc_plus_obs%dz%1.2f.txt",OutputDir,observation_number_,(double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.));
  if(!(file_ = fopen(file_name_, "w_")))
  {
    char error_message_[1000];
    sprintf(error_message_, "can't_ open file `%s'\n", file_name_);
    terminate(error_message_);
  }
#endif

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    gas_[bin_number_] = 0.;
    all_[bin_number_] = 0.;
  }
  
  for(galaxy_number_=0;galaxy_number_<TotMCMCGals[output_number_];galaxy_number_++)
  {
    for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
    {
      if(MCMC_GAL[output_number_][galaxy_number_].StellarMass >= MCMC_Obs[observation_number_].Bin_low [output_number_][bin_number_] &&
         MCMC_GAL[output_number_][galaxy_number_].StellarMass <= MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_])
      {
        gas_[bin_number_] += pow(10., MCMC_GAL[output_number_][galaxy_number_].ColdGas - MCMC_GAL[output_number_][galaxy_number_].StellarMass) * MCMC_GAL[output_number_][galaxy_number_].Weight;
        all_[bin_number_] += MCMC_GAL[output_number_][galaxy_number_].Weight;
      }
    }
  }
  
  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {  
    binned_gas_fraction_[bin_number_] =  (all_[bin_number_] > 0.) ? gas_[bin_number_] / all_[bin_number_] : 0.;
#ifndef PARALLEL
    fprintf(file_, "%g %g %g %g\n",  MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]+(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2., MCMC_Obs[observation_number_].Obs[output_number_][bin_number_],
                                                MCMC_Obs[observation_number_].Error[output_number_][bin_number_], binned_gas_fraction_[bin_number_]);
#endif
  }
#ifndef PARALLEL
  fclose(file_);
#endif
}


/**@brief Bin SAM colours according to observations of Baldry2004
*
* @note no stellar mass error?
*
*/
void bin_color_hist(const int output_number_, const int observation_number_, double *binned_color_hist_)
{
  int bin_number_, galaxy_number_;
  double color_;

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  { binned_color_hist_[bin_number_] = 0; }
  
  for(galaxy_number_=0;galaxy_number_<TotMCMCGals[output_number_];galaxy_number_++)
  {
    color_=MCMC_GAL[output_number_][galaxy_number_].Magg-MCMC_GAL[output_number_][galaxy_number_].Magr;

    for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
    {
      if(color_ >= MCMC_Obs[observation_number_].Bin_low [output_number_][bin_number_] && 
         color_ <= MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_] && 
         MCMC_GAL[output_number_][galaxy_number_].StellarMass > 9.0 && MCMC_GAL[output_number_][galaxy_number_].StellarMass < 9.5)
      
      { binned_color_hist_[bin_number_]=binned_color_hist_[bin_number_]+MCMC_GAL[output_number_][galaxy_number_].Weight; }
    }
  }
  
  // for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  // {
  //   printf("%f %f %f\n",(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2., MCMC_Obs[observation_number_].Obs[output_number_][bin_number_],binned_color_hist_[bin_number_]);
  // }
}


/**@brief Bin SAM BH and BM according to observations of Haring and Rix 2004*/
void bin_bhbm(const int output_number_, double *binned_black_hole_up_, double *binned_black_hole_down_)
{
  int galaxy_number_;

  //for MCconnell & Ma 2012
  for(galaxy_number_ = 0; galaxy_number_ < TotMCMCGals[output_number_]; galaxy_number_++)
  {
    if(MCMC_GAL[output_number_][galaxy_number_].BlackHoleMass>1.05*MCMC_GAL[output_number_][galaxy_number_].BulgeMass-2.91961)
    {
      if((MCMC_GAL[output_number_][galaxy_number_].BlackHoleMass>-0.952381*MCMC_GAL[output_number_][galaxy_number_].BulgeMass+17.2) &&
         (MCMC_GAL[output_number_][galaxy_number_].BlackHoleMass<-0.952381*MCMC_GAL[output_number_][galaxy_number_].BulgeMass+18.88))
      { binned_black_hole_up_[0]+=MCMC_GAL[output_number_][galaxy_number_].Weight; }         
      else if (MCMC_GAL[output_number_][galaxy_number_].BlackHoleMass>-0.952381*MCMC_GAL[output_number_][galaxy_number_].BulgeMass+18.88)
      { binned_black_hole_up_[1]+=MCMC_GAL[output_number_][galaxy_number_].Weight; }
    }
    else if((MCMC_GAL[output_number_][galaxy_number_].BlackHoleMass>-0.952381*MCMC_GAL[output_number_][galaxy_number_].BulgeMass+17.2) &&
            (MCMC_GAL[output_number_][galaxy_number_].BlackHoleMass<-0.952381*MCMC_GAL[output_number_][galaxy_number_].BulgeMass+18.88))
    { binned_black_hole_down_[0]+=MCMC_GAL[output_number_][galaxy_number_].Weight; }
    else if(MCMC_GAL[output_number_][galaxy_number_].BlackHoleMass>-0.952381*MCMC_GAL[output_number_][galaxy_number_].BulgeMass+18.88)
    { binned_black_hole_down_[1]+=MCMC_GAL[output_number_][galaxy_number_].Weight; }
  }
  //printf("binup0=%f bindown0=%f binup1=%f bindown1=%f\n",
  //          binned_black_hole_up_[0], binned_black_hole_down_[0], binned_black_hole_up_[1], binned_black_hole_down_[1]);
}


/** @brief compute chi^2 (uncorrelated multivariate normal) probability for sam data */
double chi_square_probability(const int output_number_, const int observation_number_, const double *sam_data_)
{
  int df_, bin_number_;
  const int knstrn_ = 0;
  double chi_squared_ = 0.0, temp_ = 0.0, prob_ = 0.0;
  double obs_, obs_error_;

  df_ = Nbins[output_number_][observation_number_] - knstrn_;

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    obs_ = MCMC_Obs[observation_number_].Obs[output_number_][bin_number_];
    obs_error_ = MCMC_Obs[observation_number_].Error[output_number_][bin_number_];

    if(obs_ < 0.0 || (obs_ == 0 && sam_data_[bin_number_] > 0))
          printf("Bad expected number in chsone\n");

    if(obs_ == 0 && sam_data_[bin_number_] == 0)
      --df_;
    else
    {
      if((obs_error_/obs_)< MCMC_Minimum_Obs_Error)
        obs_error_= MCMC_Minimum_Obs_Error*obs_;
    }

    temp_=sam_data_[bin_number_]-obs_;
    chi_squared_ += (temp_*temp_)/(obs_error_*obs_error_)/MCMC_Obs[observation_number_].ObsTest_Weight_z[output_number_];
    //if(observation_number_==0)
    // printf("OBS[%d] output_number_=%d bin=%f sam=%f obs_=%f error=%f chi_squared_=%f\n",
    //                 observation_number_, output_number_, MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]+(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2.,
    //                 sam_data_[bin_number_], obs_, obs_error_,(temp_*temp_)/(obs_error_*obs_error_));
  }

  prob_=exp(-chi_squared_/2.);
  return prob_;
}


/** @brief compute maximum likelihood probability for sam data */
double maximum_likelihood_probability(const int output_number_, const int observation_number_, const double *sam_fract_)
{
  double fract_error_ = 0.025, aux_prob_ = 99.0, prob_ = 99.0;
  int bin_number_;
  double obs_, obs_error_;

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    obs_ = MCMC_Obs[observation_number_].Obs[output_number_][bin_number_];
    obs_error_ = MCMC_Obs[observation_number_].Error[output_number_][bin_number_];
    aux_prob_=exp(-(pow2(sam_fract_[bin_number_]-obs_)/(2.*(fract_error_*fract_error_+obs_error_*obs_error_)))/MCMC_Obs[observation_number_].ObsTest_Weight_z[output_number_]);
    if(bin_number_==0) prob_=aux_prob_;
    else prob_=prob_*aux_prob_;
  }

  return prob_;
}


/** @brief compute binomial probability for sam data */
double binomial_probability(const int output_number_, const int observation_number_, const double *sam_up_, const double *sam_down_)
{
  int bin_number_;
  double prob_ = 99.0, aux_prob_ = 99.0;
  double obs_up_, obs_down_;

  for(bin_number_ = 0; bin_number_ < Nbins[output_number_][observation_number_]; bin_number_++)
  {
    obs_up_ = MCMC_Obs[observation_number_].ObsUp[output_number_][bin_number_];
    obs_down_ = MCMC_Obs[observation_number_].ObsDown[output_number_][bin_number_];

    if((1.0-betai(obs_up_,obs_down_+1.0,sam_up_[bin_number_]/(sam_up_[bin_number_]+sam_down_[bin_number_])))>betai(obs_up_,obs_down_+1.0,sam_up_[bin_number_]/(sam_up_[bin_number_]+sam_down_[bin_number_])))
      aux_prob_=2.0*betai(obs_up_,obs_down_+1.0,sam_up_[bin_number_]/(sam_up_[bin_number_]+sam_down_[bin_number_]))+1e-20;
    else
      aux_prob_=2.0*(1.0-betai(obs_up_,obs_down_+1.0,sam_up_[bin_number_]/(sam_up_[bin_number_]+sam_down_[bin_number_])))+1e-20;

    if(bin_number_==0) prob_=aux_prob_;
    else prob_=prob_*aux_prob_;
  }
  //printf("\nBinomial Probability=%f\n",prob_);
  return prob_;
}


#ifdef HALOMODEL
/** @brief compute correction for correlation */
void correct_for_correlation(const int output_number_)
{
  int fof_number_, galaxy_number_, other_galaxy_number_, *cumulative_number_of_galaxies_;

  if((HashTable = malloc(sizeof(int) * TotMCMCGals[output_number_])) == NULL)
    terminate("correct_for_correlation: malloc for HashTable failed");

  MCMC_FOF2 = malloc(sizeof(struct MCMC_FOF_struct) * NFofsInSample[output_number_]);
  cumulative_number_of_galaxies_ = malloc(sizeof(int) * NFofsInSample[output_number_]);

  cumulative_number_of_galaxies_[0] = 0;

  for(fof_number_ = 0; fof_number_ < NFofsInSample[output_number_]; fof_number_++)
  {
    if (fof_number_<NFofsInSample[output_number_]-1)
    { cumulative_number_of_galaxies_[fof_number_+1] = cumulative_number_of_galaxies_[fof_number_]+MCMC_FOF[output_number_][fof_number_].NGalsInFoF; }
    MCMC_FOF2[fof_number_].M_Crit200 = MCMC_FOF[output_number_][fof_number_].M_Crit200;
    MCMC_FOF2[fof_number_].M_Mean200 = MCMC_FOF[output_number_][fof_number_].M_Mean200;
    if(MCMC_FOF[output_number_][fof_number_].NGalsInFoF>0)
    {
      MCMC_FOF2[fof_number_].NGalsInFoF        = MCMC_FOF[output_number_][fof_number_].NGalsInFoF;
      MCMC_FOF2[fof_number_].IndexOfCentralGal = MCMC_FOF[output_number_][fof_number_].IndexOfCentralGal;
    }
    else
    {
      MCMC_FOF2[fof_number_].NGalsInFoF = 0;
      MCMC_FOF2[fof_number_].IndexOfCentralGal = -1;
    }
  }
  UsedFofsInSample[output_number_]=NFofsInSample[output_number_];

  //BUILD HASHTABLE
  for(galaxy_number_=0;galaxy_number_<TotMCMCGals[output_number_];galaxy_number_++)
    HashTable[galaxy_number_]=-1;

  for(fof_number_=0;fof_number_<NFofsInSample[output_number_]; fof_number_++)
  {
    for(galaxy_number_=0;galaxy_number_<TotMCMCGals[output_number_];galaxy_number_++)
    {
      if(MCMC_GAL[output_number_][galaxy_number_].fofid==fof_number_)
      {
        MCMC_GAL[output_number_][galaxy_number_].ngal=MCMC_FOF2[fof_number_].NGalsInFoF;
        //if type=0 it gets the first place in hashtable for this group (cumulative_number_of_galaxies_[fof_number_])
        if(MCMC_GAL[output_number_][galaxy_number_].Type==0)
        {
          HashTable[cumulative_number_of_galaxies_[fof_number_]]=galaxy_number_;
          MCMC_FOF2[fof_number_].IndexOfCentralGal = cumulative_number_of_galaxies_[fof_number_];
        }
        //if not gets the first available place in hashtable for this group (cumulative_number_of_galaxies_[fof_number_]+other_galaxy_number_)
        else
        {
          for(other_galaxy_number_ = 1; other_galaxy_number_ < MCMC_FOF2[fof_number_].NGalsInFoF; other_galaxy_number_++)
          {
            if(HashTable[cumulative_number_of_galaxies_[fof_number_]+other_galaxy_number_]==-1)
            {
              HashTable[cumulative_number_of_galaxies_[fof_number_]+other_galaxy_number_]=galaxy_number_;
              break;
            }
          }
        }
      }
    }
  }//loop on fof_number_ to get hashtable

  free(cumulative_number_of_galaxies_);
}


/** @brief compute the correlation function for a_ given galaxy mass range */
void compute_correlation_func(const int output_number_, const int observation_number_, const float min_galaxy_mass_, const float max_galaxy_mass_, double *binned_sam_data_)
{
  double *r_, *proj_;
  int bin_number_;
#ifndef PARALLEL
  char file_name_[1000];
  FILE *file_;
#endif
  gsl_spline *Proj_Spline_;
  gsl_interp_accel *Proj_SplineAcc_;

  NR=60;
  r_=malloc(NR*sizeof(double));
  proj_=malloc(NR*sizeof(double));

  halomodel(r_,proj_,min_galaxy_mass_,max_galaxy_mass_,output_number_);

  Proj_Spline_    = gsl_spline_alloc(gsl_interp_cspline,NR);
  Proj_SplineAcc_ = gsl_interp_accel_alloc();
  gsl_spline_init(Proj_Spline_,r_,proj_,NR);

  for(bin_number_=0;bin_number_<Nbins[output_number_][observation_number_]-1;bin_number_++)
    binned_sam_data_[bin_number_]=gsl_spline_eval(Proj_Spline_,MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]+(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2.,Proj_SplineAcc);

#ifndef PARALLEL
  //full3 - full PLANCK
  sprintf(file_name_, "%s/correlation_guo10_bug_fix_Mmean_z0.00_%0.2f_%0.2f.txt",OutputDir, min_galaxy_mass_,max_galaxy_mass_);
  if(!(file_ = fopen(file_name_, "w_")))
  {
    char error_message_[1000];
    sprintf(error_message_, "can't_ open file `%s'\n", file_name_);
    terminate(error_message_);
  }
  for(bin_number_=0;bin_number_<Nbins[output_number_][observation_number_]-1;bin_number_++)
    fprintf(file_, "%g %g %g\n", MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_]+(MCMC_Obs[observation_number_].Bin_high[output_number_][bin_number_]-MCMC_Obs[observation_number_].Bin_low[output_number_][bin_number_])/2.,
            binned_sam_data_[bin_number_],binned_sam_data_[bin_number_]*0.1);
  fclose(file_);
#endif
  free(r_);
  free(proj_);
  gsl_spline_free(Proj_Spline_);
  gsl_interp_accel_free(Proj_SplineAcc_);
}
#endif //HALOMODEL


#define FPMIN 1.0e-30
#define ASWITCH 100
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MAXIT 100

double gammp(const double a_, const double x_)
{
  if(x_ < 0.0 || a_ <= 0.0)
    printf("bad args in gammp\n");
  if(x_ == 0.0)
    return 0.0;
  else if((int) a_ >= ASWITCH)
    return gammpapprox(a_, x_, 1);
  else if(x_ < a_ + 1.0)
    return gser(a_, x_);
  else
    return 1.0 - gcf(a_, x_);
}


double gammq(const double a_, const double x_)
{
  if(x_ < 0.0 || a_ <= 0.0)
    printf("bad args in gammq\n");
  if(x_ == 0.0)
    return 1.0;
  else if((int) a_ >= ASWITCH)
    return gammpapprox(a_, x_, 0);
  else if(x_ < a_ + 1.0)
    return 1.0 - gser(a_, x_);
  else
    return gcf(a_, x_);
}


double gser(const double a_, const double x_)
{
  double sum_, del_, ap_;
  int n_;

  const double gammln_a_ = gammln(a_);
  
  ap_ = a_;
  del_ = sum_ = 1.0 / a_;

  for(n_ = 1; n_ <= ASWITCH; n_++)
  {
    ++ap_;
    del_ *= x_ / ap_;
    sum_ += del_;
    if(fabs(del_) < fabs(sum_)*EPS) return sum_*exp(-x_+a_*log(x_)-gammln_a_);
  }
  // return regardless(?):
  return sum_*exp(-x_+a_*log(x_)-gammln_a_);
}


double gcf(const double a_, const double x_)
{
  int i_;
  double an, b_, c_, d_, del_, h_;

  const double gammln_a_ = gammln(a_);
  
  b_ = x_ + 1.0 - a_;
  c_ = 1.0 / FPMIN;
  d_ = 1.0 / b_;
  h_ = d_;

  for(i_ = 1; i_ <= ASWITCH; i_++)
  {
    an = -i_ * (i_ - a_);
    b_ += 2.0;
    d_ = an * d_ + b_;
    if(fabs(d_) < FPMIN)
          d_ = FPMIN;
    c_ = b_ + an / c_;
    if(fabs(c_) < FPMIN)
          c_ = FPMIN;
    d_ = 1.0 / d_;
    del_ = d_ * c_;
    h_ *= del_;
    if(fabs(del_-1.0) <= EPS) break;
  }
  return exp(-x_ + a_ * log(x_) - gammln_a_) * h_;
}


double gammpapprox(const double a_, const double x_, const int psig_)
{
  const double y_[18] = { 0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
    0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
    0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483,
    0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
    0.87126389619061517, 0.95698180152629142
  };
  const double w_[18] = { 0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
    0.034213810770299537, 0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
    0.058860144245324798, 0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
    0.076598410645870640, 0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
    0.085346685739338721, 0.085983275670394821
  };
  int j_;
  double xu_, t_, sum_, ans_;
  
  const int ngau_ = 18;
  const double a_minus_1_ = a_ - 1.0;
  const double ln_a_minus_1_ = log(a_minus_1_);
  const double sqrt_a_minus_1_ = sqrt(a_minus_1_);
  const double gammln_a_ = gammln(a_);

  if(x_ > a_minus_1_)
    if((a_minus_1_ + 11.5 * sqrt_a_minus_1_) > (x_ + 6.0 * sqrt_a_minus_1_))
      xu_ = a_minus_1_ + 11.5 * sqrt_a_minus_1_;
    else
      xu_ = x_ + 6.0 * sqrt_a_minus_1_;
  else if((a_minus_1_ - 7.5 * sqrt_a_minus_1_) < (x_ - 5.0 * sqrt_a_minus_1_))
    if(0. > (a_minus_1_ - 7.5 * sqrt_a_minus_1_))
      xu_ = 0;
    else
      xu_ = (a_minus_1_ - 7.5 * sqrt_a_minus_1_);
  else if(0. > (x_ - 5.0 * sqrt_a_minus_1_))
    xu_ = 0;
  else
    if(0.>( x_ - 5.0*sqrt_a_minus_1_)) xu_ =0;
    else xu_=( x_ - 5.0*sqrt_a_minus_1_);

  sum_ = 0;
  for(j_ = 0; j_ < ngau_; j_++)
  {
    t_ = x_ + (xu_ - x_) * y_[j_];
    sum_ += w_[j_] * exp(-(t_ - a_minus_1_) + a_minus_1_ * (log(t_) - ln_a_minus_1_));
  }
  ans_ = sum_ * (xu_ - x_) * exp(a_minus_1_ * (ln_a_minus_1_ - 1.) - gammln_a_);

  if(psig_)
  {
    if(ans_ > 0.0)
    { return ans_ + 1.; }
    else
    { return -ans_; }
  }
  else
  {
    if(ans_ >= 0.0)
    { return ans_; }
    else
    { return ans_ + 1.; }
  }
}


double gammln(const double xx_)
{
  int j_;
  double x_, tmp_, y_ = 0, ser_;

  double cof[14] = { 57.1562356658629235, -59.5979603554754912,
    14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
    .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
    -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
    .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5
  };
  if(xx_ <= 0)
    printf("bad arg in gammln\n");
  y_ = x_ = xx_;
  tmp_ = x_ + 5.24218750000000000;
  tmp_ = (x_ + 0.5) * log(tmp_) - tmp_;
  ser_ = 0.999999999999997092;
  for(j_ = 0; j_ < 14; j_++)
    ser_ += cof[j_] / ++y_;
  return tmp_ + log(2.5066282746310005 * ser_ / x_);
}


double betai(const double a_, const double b_, const double x_)
{
  double bt_;
  if(x_ < 0.0 || x_ > 1.0)
    printf("Bad x_ in routine betai\n");
  if(x_ == 0.0 || x_ == 1.0)
    bt_ = 0.0;
  else
    bt_ = exp(gammln(a_ + b_) - gammln(a_) - gammln(b_) + a_ * log(x_) + b_ * log(1.0 - x_));
  if(x_ < (a_ + 1.0) / (a_ + b_ + 2.0))
    return bt_ * betacf(a_, b_, x_) / a_;
  else
    return 1.0 - bt_ * betacf(b_, a_, 1.0 - x_) / b_;
}


double betacf(const double a_, const double b_, const double x_)
{
  int m_, m2_;
  double aa_, c_, d_, del_, h_, qab_, qam_, qap_;

  qab_=a_+b_;
  qap_=a_+1.0;
  qam_=a_-1.0;
  c_=1.0;
  d_=1.0-qab_*x_/qap_;
  if (fabs(d_) < FPMIN) d_=FPMIN;
  d_=1.0/d_;
  h_=d_;
  for (m_=1;m_<=MAXIT;m_++)
  {
    m2_=2*m_;
    aa_=m_*(b_-m_)*x_/((qam_+m2_)*(a_+m2_));
    d_=1.0+aa_*d_;
    if (fabs(d_) < FPMIN) d_=FPMIN;
    c_=1.0+aa_/c_;
    if (fabs(c_) < FPMIN) c_=FPMIN;
    d_=1.0/d_;
    h_ *= d_*c_;
    aa_ = -(a_+m_)*(qab_+m_)*x_/((a_+m2_)*(qap_+m2_));
    d_=1.0+aa_*d_;
    if (fabs(d_) < FPMIN) d_=FPMIN;
    c_=1.0+aa_/c_;
    if (fabs(c_) < FPMIN) c_=FPMIN;
    d_=1.0/d_;
    del_=d_*c_;
    h_ *= del_;
    if (fabs(del_-1.0) < EPS) break;
  }
  if (m_ > MAXIT) printf("a_ or b_ too big, or MAXIT too small in betacf\n");
  return h_;
}
#endif //MCMC

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#undef FPMIN
#undef ASWITCH
#undef MAXIT
