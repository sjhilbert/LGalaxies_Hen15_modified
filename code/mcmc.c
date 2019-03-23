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

/** @file   mcmc.c
 *  @date   2008, 2018
 *  @author Bruno Henriques
 *  @author Stefan Hilbert
 *
 *  @brief  This is the main mcmc file, it reads in observational tests,
 *          starts a chain, proposes new sets of parameters, calls the
 *          SAM, gets the likelihood for galaxy properties obtained with
 *          the current parameters, compares the likelihood with the
 *          previous run and decides on whether or not to accept the
 *          proposed parameters.
 *
 * To run you need to chose input_mcmc_wmp7.par and select the
 * desired sample of halos. Also choose the correct desired_output_snap
 * according to the redshift at which you want to run the constraints.
 * Turn on MCMC on the Makefile.
 *
 * Senna is the controlling routine. First is calls read_observations,
 * to read in all the observational data sets into the struct MCMC_Obs.
 * Then, an initial SAM is run, to obtain a first likelihood value.
 * This will be the starting point of the chain, given by the initial
 * values of p[] (this should be some non 0 likelihood region to avoid
 * a long burn in). The likelihood for each set of parameters is computed
 * in mcmc_likelihood.c (get_likelihood()) at the end of each SAM run.
 *
 * After this first step, a chain is started with mcmc() being called
 * for each step to control all the processes. Namely, propose_new_parameters(),
 * SAM() and then at each step the likelihood for the given set of
 * parameters compared with previous and accepted or not. Two options
 * available according to the value of MCMCMode in input.par. If MCMCMode
 * =0, normal MCMC is done. If MCMCMode = 1, the new set of parameters
 * is only accepted if the likelihood is higher than for the previous.
 * This makes the chain going up in likelihood very quickly and its
 * useful to find a region of non 0 likelihood when the parameter is
 * still unknown.
 *
 * The name of the parameters sampled can be seen in the beginning of
 * function SAM() where the values of the internal variables of the code
 * containing the parameters are changed according to the new set of
 * parameters proposed. At each step p[] will contain the previously
 * accepted set of parameters and prop[] will contain the newly proposed
 * set of parameters (given by propose_new_parameters()). If desired,
 * less parameters then the default number, given by MCMCNpar can be sampled.
 * This can be done in propose_new_parameters() by never assigning a new
 * value to a given parameter. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"


#ifdef HALOMODEL
void initialize_halomodel(void);
#endif /* defined HALOMODEL */


#ifdef MCMC
/** @brief main routine for controlling MCMC algorithm */
void Senna(void)
{
  int parameter_number_, output_number_, N_MCMC_steps_accepted_ = 0;
  FILE *mcmc_file_;
  char buf_[1000];
  char dir_num_;
  time_t local_initial_time_, local_final_time_;
  double proposed_likelihood_, q_ratio_, a_ratio_;
  bool proposed_MCMC_step_is_accepted_;

  if(ThisTask==0)
  {
        printf("\n\n\n");
        printf("**********************************************************\n");
        printf("*                                                        *\n");
        printf("*                   Starting Senna                       *\n");
        printf("*                                                        *\n");
        printf("*             MCMC parameter estimation                  *\n");
        printf("*  Applied to a Semi-Analytic Model of Galaxy Formation  *\n");
        printf("*                                                        *\n");
        printf("**********************************************************\n\n");
  }

  MCMC_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
  MCMCseed = -((ThisTask + FirstChainNumber) * 100 + 25);
  gsl_rng_set(MCMC_rng, MCMCseed);	/* start-up seed */ 

#ifdef HALOMODEL //to compute correlation function for MCMC
  initialize_halomodel();
  printf("halo model initialized\n");
#endif

  //This file will have the values for the parameters
  //It will be the output from the mcmc sampling used in the analysis
  //it is also where the parameters are read from, using previous chains
  getcwd(buf_, sizeof(buf_));
  dir_num_=buf_[strlen(buf_)-1];
  sprintf(buf_, "%s/senna_g%c_%d.txt", OutputDir, dir_num_, ThisTask+FirstChainNumber);
  //open to read and write. if file doesn't exists open just to write and read from backup
  if((mcmc_file_ = fopen(buf_, "r+")) == NULL)
    if((mcmc_file_ = fopen(buf_, "w")) == NULL)
    {
      char error_message_[1000];
      sprintf(error_message_, "can't open file `%s'\n", buf_);
      terminate(error_message_);
    }

  //read (from previous output) and initialize mcmc parameters, also MCMC_Likelihood
  initialize_mcmc_par_and_lhood (mcmc_file_);

  //Read observational constraints (SMF, LFs, BHBM relation, etc)
  read_observations();

  //open files that will contain the comparison with observations at each step (used to compute the likelihoods)
  open_files_with_comparison_to_observations();

  for(CurrentMCMCStep = 0; CurrentMCMCStep < ChainLength; ++CurrentMCMCStep, ++GlobalMCMCStep)
  {
    time(&local_initial_time_);
    printf("\n\n\nMCMC %d STARTED on Task %d\n\n\n", CurrentMCMCStep, ThisTask + FirstChainNumber);

    //get a new set of parameters and return q_ratio_ - the prior
    q_ratio_ = propose_new_parameters();
        
    //runs the SAM with the new parameters and gives the likelyhood for them
    proposed_likelihood_ = SAM(MCMCSampleFile);

    printf("LIKELY1=%0.5e\nLIKELY2=%0.5e\n", MCMC_Likelihood, proposed_likelihood_);

    if(isnan(MCMC_Likelihood))      MCMC_Likelihood      = 0.;
    if(isnan(proposed_likelihood_)) proposed_likelihood_ = 0.;

    // accepts the proposed parameters with probability=acceptance rate:
    if(MCMCMode == 0)
    {
      /* By default q_ratio_ = 1, meaning we assume a flat prior. Therefore, the acceptance
       * probability is just given by the ratio of the likelihoods from two steps. */
      a_ratio_ = q_ratio_ * (proposed_likelihood_ / MCMC_Likelihood);
      proposed_MCMC_step_is_accepted_ = (MCMC_Likelihood <= 0.) || (1. <= a_ratio_) || (gsl_rng_uniform(MCMC_rng) < a_ratio_);
    }
    // only accepts new parameters if likelihood increases:
    else if(MCMCMode == 1)
    { proposed_MCMC_step_is_accepted_ = (proposed_likelihood_ >= MCMC_Likelihood); }
    // don't accept in unknown mode:
    else
    { proposed_MCMC_step_is_accepted_ = false; }

    //if new set of parameters is accepted change current set of parameter and lhood
    if(proposed_MCMC_step_is_accepted_)
    {
      N_MCMC_steps_accepted_++;
      MCMC_Likelihood = proposed_likelihood_;
      for(parameter_number_ = 0; parameter_number_ < MCMCNpar; parameter_number_++)
        for(output_number_ = 0; output_number_ < NOUT; output_number_++)
          MCMC_PAR[parameter_number_].Value[output_number_] = MCMC_PAR[parameter_number_].PropValue[output_number_];
    }

    //print the values to file and to screen (to screen only if proposed_MCMC_step_is_accepted_)
    print_parameters(proposed_MCMC_step_is_accepted_, mcmc_file_);

    printf("\nCurrent acceptance rate of this chain=%0.2f%%\n", (N_MCMC_steps_accepted_ / (float)(CurrentMCMCStep + 1)) * 100.);

    time(&local_final_time_);
    printf("Task %d chain %d (%d) took %lds\n", ThisTask+FirstChainNumber, CurrentMCMCStep, GlobalMCMCStep, local_final_time_ - local_initial_time_);
    printf("Global Time Elapsed %lds, %f hours\n", local_final_time_ - GlobalStartingTime, (local_final_time_ - GlobalStartingTime) / 3600.);

#ifdef PARALLEL
    if(ThisTask==0)
      if((local_final_time_ - GlobalStartingTime)/3600. > MachineTimeOut)
      {
        sprintf(buf_, "%s %s %s%d.bash", JobSubmitCommand, JobSubmitPipe, JobSubmitFile, FirstChainNumber);
        printf("resubmit command: %s\n",buf_);
        system(buf_);
        fflush(stdout);
        terminate("\n\nMachine TimeOut reached, restarting\n\n");
      }
#endif
  }// ChainLength (MCMC steps)

  fclose(mcmc_file_);

  /* For each set of parameters writes the semi-analytic
   * predictions used for the likelihood calculation*/
  //write_comparison_to_obs();
  close_files_with_comparison_to_observations();

  myfree(MCMC_Obs);

#ifdef HALOMODEL //to compute correlation function for MCMC
    gsl_spline_free(FofSpline);
    gsl_interp_accel_free(FofAcc);
    gsl_spline_free(SigmaSpline);
    gsl_interp_accel_free(SigmaAcc);
    gsl_spline_free(ellipSpline);
    gsl_interp_accel_free(ellipAcc);
    gsl_spline_free(PowSpline);
#endif

  gsl_rng_free(MCMC_rng);

  printf("\nFinal acceptance rate of this chain=%f%%\n", ((float) N_MCMC_steps_accepted_ / ChainLength) * 100);
  printf("\n\nMCMC OVER\n\n");
}


////////
//MCMC//
////////

/** @brief Function to print parameters whenever
 *       a step is accepted in the MCMC */
void print_parameters (const bool proposed_MCMC_step_is_accepted_, FILE *mcmc_file_)
{
  int parameter_number_, output_number_;
  const int chainweight_ = 1;

  fprintf(mcmc_file_,"%d %0.8g", chainweight_, -(log10(MCMC_Likelihood)));

  //print parameters into output file
  for(parameter_number_=0;parameter_number_<MCMCNpar;++parameter_number_)
  {
    if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
    {
      if(Time_Dependent_PhysPar==1)
        for(output_number_ = 0; output_number_ < NOUT; output_number_++)
          fprintf(mcmc_file_," %0.6f", log10(MCMC_PAR[parameter_number_].Value[output_number_]));
      else
        fprintf(mcmc_file_," %0.6f", log10(MCMC_PAR[parameter_number_].Value[0]));
    }
  }
  fprintf(mcmc_file_,"\n");
  fflush(mcmc_file_);
  fflush(stdout);

  //print to screen
  if(proposed_MCMC_step_is_accepted_)
  {
    printf("\n******************************************************\n");
    printf("Accepted!!!\n");
    for(parameter_number_=0;parameter_number_<MCMCNpar;++parameter_number_)
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
      {
        if(Time_Dependent_PhysPar==1)
          for(output_number_=0;output_number_<NOUT;output_number_++)
            printf("%0.2g ",MCMC_PAR[parameter_number_].PropValue[output_number_]);
        else
          printf("%0.2g ",MCMC_PAR[parameter_number_].PropValue[0]);
      }
    printf("\n");


    //print log of parameters to screen
    printf("%d %0.8g ",chainweight_, -(log10(MCMC_Likelihood)));
    for(parameter_number_=0;parameter_number_<MCMCNpar;++parameter_number_)
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
      {
        if(Time_Dependent_PhysPar==1)
          for(output_number_=0;output_number_<NOUT;output_number_++)
            printf("%0.6f ",log10(MCMC_PAR[parameter_number_].PropValue[output_number_]));
        else
          printf("%0.6f ",log10(MCMC_PAR[parameter_number_].PropValue[0]));
      }

    printf("\n******************************************************\n\n\n\n");
  } //end if(proposed_MCMC_step_is_accepted_==1) - print ot screen
}


/** @brief initialize MCMC_PAR.Value and MCMC_PAR.PropValue with the same values
 *
 * the number of parameters, limits and switches (whitch to sample) are
 * read from MCMCParaPriorsAndSwitchesFile while the actual values are read from
 * previous output or MCMCStartingParFile
 *  */
void initialize_mcmc_par_and_lhood (FILE *mcmc_file_)
{
  int parameter_number_, output_number_, dumb_weight_;
  double aux_p_;
  char buf_[1000];
  FILE *file_;

  sprintf(buf_, "%s", MCMCParPriorsAndSwitchesFile);
  if(!(file_ = fopen(buf_, "r")))
  {
    char error_message_[1000];
    sprintf(error_message_, "can't open file `%s'\n", buf_);
    terminate(error_message_);
  }

  fgets(buf_, 300, file_);
  fscanf(file_,"%d\n",&MCMCNpar);

  //initialize structure to contain parameter names, values, priors and other properties
  MCMC_PAR = mymalloc("MCMC_PAR", sizeof(struct MCMC_PAR) * MCMCNpar);

  //read names and switches
  fgets(buf_, 300, file_);
  for(parameter_number_ = 0; parameter_number_ < MCMCNpar; parameter_number_++)
  {
    fscanf(file_,"%s %lg %lg %lg %d\n",MCMC_PAR[parameter_number_].Name, &MCMC_PAR[parameter_number_].PropValue[0],
           &MCMC_PAR[parameter_number_].PriorMin, &MCMC_PAR[parameter_number_].PriorMax, &MCMC_PAR[parameter_number_].Sampling_Switch);
  }

  fclose(file_);  //done reading from MCMCParPriorsAndSwitchesFile

  //read actual values from previous outputs and check if are inside priors
  GlobalMCMCStep = 0;
  while(2 == fscanf(mcmc_file_,"%d %lg ", &dumb_weight_, &MCMC_Likelihood))
  {
    //read parameter values
    for(parameter_number_=0;parameter_number_<MCMCNpar;++parameter_number_)
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
        fscanf(mcmc_file_,"%lg",&MCMC_PAR[parameter_number_].Value[0]);
    //printf("par[%d]=%f\n",parameter_number_, MCMC_PAR[parameter_number_].Value[0]);
    ++GlobalMCMCStep;
  }

  //if there is something in that file,this is a re-start so:
  if(GlobalMCMCStep > 0)
  { MCMC_Initial_Par_Displacement=0.; }
  else
  //if there is nothing in that file read from backup file 00
  {
    sprintf(buf_, "%s", MCMCStartingParFile);
    if((file_ = fopen(buf_, "r")) == NULL)
    {
      char error_message_[1000];
      sprintf(error_message_, "can't open file `%s'\n", buf_);
      terminate(error_message_);
    }
    fscanf(file_,"%d %lg ", &dumb_weight_, &MCMC_Likelihood);
    for(parameter_number_=0;parameter_number_<MCMCNpar;++parameter_number_)
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
        fscanf(file_,"%lg",&MCMC_PAR[parameter_number_].Value[0]);
    fclose(file_);
  }

  //convert from log
  MCMC_Likelihood=pow(10,-MCMC_Likelihood);
  for(parameter_number_=0;parameter_number_<MCMCNpar;++parameter_number_)
    if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
      MCMC_PAR[parameter_number_].Value[0]=pow(10,MCMC_PAR[parameter_number_].Value[0]);

  printf("Initial Parameter Values:\n");
  for(parameter_number_=0;parameter_number_<MCMCNpar;parameter_number_++)
  {
    if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
    { printf("%0.6f ",MCMC_PAR[parameter_number_].Value[0]); }
    if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
      if(MCMC_PAR[parameter_number_].Value[0]<MCMC_PAR[parameter_number_].PriorMin ||
         MCMC_PAR[parameter_number_].Value[0]>MCMC_PAR[parameter_number_].PriorMax)
      {
        printf("value=%0.6f priormin=%0.4f priormax=%0.4f\n",MCMC_PAR[parameter_number_].Value[0],MCMC_PAR[parameter_number_].PriorMin,MCMC_PAR[parameter_number_].PriorMax);
        char error_message_[1000];
        sprintf(error_message_, "parameter '%s' outside prior range \n", MCMC_PAR[parameter_number_].Name);
        terminate(error_message_);
      }
  }
  printf("\n");

  //if MCMC_Initial_Par_Displacement>0. introduce displacement in parameter values
  if(Time_Dependent_PhysPar==0)
  {
    for(parameter_number_ = 0; parameter_number_ < MCMCNpar; parameter_number_++)
    {
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
      {
        aux_p_ = MCMC_PAR[parameter_number_].Value[0];
        do
        { MCMC_PAR[parameter_number_].Value[0] = aux_p_ * exp(MCMC_Initial_Par_Displacement * gsl_ran_ugaussian(MCMC_rng)); }
        while(MCMC_PAR[parameter_number_].Value[0] < MCMC_PAR[parameter_number_].PriorMin
           || MCMC_PAR[parameter_number_].Value[0] > MCMC_PAR[parameter_number_].PriorMax);
      }
      else
        MCMC_PAR[parameter_number_].Value[0] = MCMC_PAR[parameter_number_].Value[0];

      MCMC_PAR[parameter_number_].PropValue[0] = MCMC_PAR[parameter_number_].Value[0];

      for(output_number_ = 1; output_number_ < NOUT; output_number_++)
      {
        MCMC_PAR[parameter_number_].Value    [output_number_] = MCMC_PAR[parameter_number_].Value[0];
        MCMC_PAR[parameter_number_].PropValue[output_number_] = MCMC_PAR[parameter_number_].Value[0];
      }
    }
  }
  else if(Time_Dependent_PhysPar==1)
  {
    for(parameter_number_ = 0; parameter_number_ < MCMCNpar; parameter_number_++)
    {
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
      {   
        for(output_number_ = 0; output_number_ < NOUT; output_number_++)
        {
          aux_p_=MCMC_PAR[parameter_number_].Value[0];
          do
          { MCMC_PAR[parameter_number_].Value[output_number_] = aux_p_ * exp(MCMC_Initial_Par_Displacement * gsl_ran_ugaussian(MCMC_rng)); }
          while(MCMC_PAR[parameter_number_].Value[output_number_] < MCMC_PAR[parameter_number_].PriorMin
             || MCMC_PAR[parameter_number_].Value[output_number_] > MCMC_PAR[parameter_number_].PriorMax);

          MCMC_PAR[parameter_number_].PropValue[output_number_] = MCMC_PAR[parameter_number_].Value[output_number_];
        }
      }
      else
      {
        for(output_number_ = 0; output_number_ < NOUT; output_number_++)
        {
          MCMC_PAR[parameter_number_].Value    [output_number_] = MCMC_PAR[parameter_number_].Value[0];
          MCMC_PAR[parameter_number_].PropValue[output_number_] = MCMC_PAR[parameter_number_].Value[0];
        }
      }
    }

    //LOAD INTITIAL VALUES FOR ALL PARAMETERS FROM A DIFFERENT FILE
    sprintf(buf_, "./input/mcmc_allz_par.txt");
    if(!(file_ = fopen(buf_, "r")))
    {
      char error_message_[1000];
      sprintf(error_message_, "can't open file `%s'\n", buf_);
      terminate(error_message_);
    }

    for(parameter_number_ = 0; parameter_number_ < MCMCNpar; parameter_number_++)
      for(output_number_ = 0; output_number_ < NOUT; output_number_++)
        if(parameter_number_ < 1 || (parameter_number_ > 2 && parameter_number_ < 13))
        {
          fscanf(file_,"%lg",&MCMC_PAR[parameter_number_].Value[output_number_]);
          aux_p_=MCMC_PAR[parameter_number_].Value[output_number_];
          do
          { MCMC_PAR[parameter_number_].Value[output_number_] = aux_p_ * exp(MCMC_Initial_Par_Displacement * gsl_ran_ugaussian(MCMC_rng)); }
          while(MCMC_PAR[parameter_number_].Value[output_number_] < MCMC_PAR[parameter_number_].PriorMin || MCMC_PAR[parameter_number_].Value[output_number_] > MCMC_PAR[parameter_number_].PriorMax);
          MCMC_PAR[parameter_number_].PropValue[output_number_] = MCMC_PAR[parameter_number_].Value[output_number_];
        }
    fclose(file_);
  }
}


///////////
//PROPOSE//
///////////

/** @brief Function to propose new parameters given by
 *       a random normal distribution gassdev*/
double propose_new_parameters(void)
{
  double q_ratio_;
  int parameter_number_, output_number_;
 
  if(Time_Dependent_PhysPar == 0)  // only sample for output_number_ == 0, then copy for other output_number_:
  {
    for(parameter_number_ = 0; parameter_number_ < MCMCNpar; parameter_number_++)
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
      {
        do
        { MCMC_PAR[parameter_number_].PropValue[0] = MCMC_PAR[parameter_number_].Value[0] * exp(MCMC_LogStep_Size * gsl_ran_ugaussian(MCMC_rng)); }
        while(MCMC_PAR[parameter_number_].PropValue[0] < MCMC_PAR[parameter_number_].PriorMin ||
              MCMC_PAR[parameter_number_].PropValue[0] > MCMC_PAR[parameter_number_].PriorMax);

        for(output_number_ = 1; output_number_ < NOUT; output_number_++)
        { MCMC_PAR[parameter_number_].PropValue[output_number_] = MCMC_PAR[parameter_number_].PropValue[0]; }
      }
  }
  else // sample for all output_number_:
  {
    for(parameter_number_ = 0; parameter_number_ < MCMCNpar; parameter_number_++)
      if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
        for(output_number_ = 0; output_number_ < NOUT; output_number_++)
        {
          do
          { MCMC_PAR[parameter_number_].PropValue[output_number_] = MCMC_PAR[parameter_number_].Value[output_number_] * exp(MCMC_LogStep_Size * gsl_ran_ugaussian(MCMC_rng)); }
          while(MCMC_PAR[parameter_number_].PropValue[output_number_] < MCMC_PAR[parameter_number_].PriorMin ||
                MCMC_PAR[parameter_number_].PropValue[output_number_] > MCMC_PAR[parameter_number_].PriorMax);
        }
  }

  q_ratio_ = 1;
  //q_ratio_= prop[0]/p[0]*prop[1]/p[1]*prop[2]/p[2]*prop[3]/p[3]*prop[4]/p[4];
  return q_ratio_;
}


/** @brief read in mcmc parameters */
void read_mcmc_par (const int snapshot_number_)
{
  int output_number_, parameter_number_;
 
  // //for snapshot_number_ between output_number_[i] and output_number_[i+1], parameters have the values of output_number_[i+1]
  // for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  //   if(snapshot_number_ < ListOutputSnaps[NOUT-output_number_-1]+1)
  //     break;
  // output_number_=NOUT-output_number_-1;
  // 
  // if(output_number_  <  0)
  //   output_number_ = 0;
  
  output_number_ = ListOutputNumberOfSnapshot[snapshot_number_];

  for(parameter_number_=0;parameter_number_<MCMCNpar;parameter_number_++)
    if(MCMC_PAR[parameter_number_].Sampling_Switch==1)
    {
      if(strcmp(MCMC_PAR[parameter_number_].Name,"SfrEfficiency")==0)
        SfrEfficiency = MCMC_PAR[parameter_number_].PropValue[output_number_];
      if(strcmp(MCMC_PAR[parameter_number_].Name,"SfrColdCrit")==0)
        SfrColdCrit = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"SfrBurstEfficiency")==0)
        SfrBurstEfficiency = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"SfrBurstSlope")==0)
        SfrBurstSlope = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"AgnEfficiency")==0)
        AgnEfficiency = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"BlackHoleGrowthRate")==0)
        BlackHoleGrowthRate = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"BlackHoleCutoffVelocity")==0)
        BlackHoleCutoffVelocity = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"FeedbackReheatingEpsilon")==0)
        FeedbackReheatingEpsilon = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"ReheatPreVelocity")==0)
        ReheatPreVelocity = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"ReheatSlope")==0)
        ReheatSlope = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"FeedbackEjectionEfficiency")==0)
        FeedbackEjectionEfficiency = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"EjectPreVelocity")==0)
        EjectPreVelocity = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"EjectSlope")==0)
        EjectSlope = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"Yield")==0)
        Yield = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"ThreshMajorMerger")==0)
        ThreshMajorMerger = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"MergerTimeMultiplier")==0)
        MergerTimeMultiplier = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"RamPressureStrip_CutOffMass")==0)
        RamPressureStrip_CutOffMass = MCMC_PAR[parameter_number_].PropValue[output_number_];

      else if(strcmp(MCMC_PAR[parameter_number_].Name,"Reionization_z0")==0)
        Reionization_z0 = MCMC_PAR[parameter_number_].PropValue[output_number_];
      else if(strcmp(MCMC_PAR[parameter_number_].Name,"Reionization_zr")==0)
      {
        Reionization_zr = MCMC_PAR[parameter_number_].PropValue[output_number_];
        if(Reionization_zr>(Reionization_z0-0.5))
        {
          Reionization_zr=Reionization_z0-0.5;
          MCMC_PAR[parameter_number_].PropValue[output_number_]=Reionization_z0-0.5;
        }
      }
                      //printf("EjectSlope=%g\n",EjectSlope);
    }
}


/** @brief Read in the IDs and Weights of the FOFs groups that
 *       constitute the sample for which galaxies will be
 *       compared with observations
 *
 *       There is a single file with merger trees that needs to be created before
 *       running the code. The file contains the trees corresponding to the IDs
 *       read here. It contains both MR and MRII trees. The MR trees are stored first.
 *       NOTE THAT THERE IS MORE THAN 1 FOF GROUP PER TREE. Another file is
 *       created with the number of trees in MR. After the code (in main()) has
 *       looped over that number of trees the necessary variables are changed to MRII.*/
void read_sample_info (void)
{
  int DumbTreeNrColector_, DumbFileNrColector_;
  int fof_number_, output_number_;
  FILE *file_;
  char file_name_[1000];
  char error_message_[1000];

#ifdef MR_PLUS_MRII
  if(Switch_MR_MRII==1)
  {
    sprintf(file_name_, "%s/%ssample_allz_nh_Switch_MR_MRII_%d.dat", MCMCSampleDir, MCMCSampleFilePrefix, MCMCSampleFile);
    if(!(file_ = fopen(file_name_, "r")))
    {
      sprintf(error_message_, "can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }
    fscanf(file_, "%d \n", &NTrees_Switch_MR_MRII);
    fclose(file_);
  }
#endif

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    sprintf(file_name_, "%s/%ssample_allz_nh_%d%d.dat", MCMCSampleDir, MCMCSampleFilePrefix, MCMCSampleFile, ListOutputSnaps[output_number_]);
    if(!(file_ = fopen(file_name_, "r")))
    {
      sprintf(error_message_, "can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }

    fscanf(file_, "%d \n", &NFofsInSample[output_number_]);

    MCMC_FOF[output_number_] = malloc(sizeof(struct MCMC_FOF_struct) * NFofsInSample[output_number_]);
    
    for(fof_number_ = 0; fof_number_ < NFofsInSample[output_number_]; fof_number_++)
    {
      fscanf(file_, "%lld %d %d %lg\n", &MCMC_FOF[output_number_][fof_number_].FoFID, &DumbTreeNrColector_, &DumbFileNrColector_, &MCMC_FOF[output_number_][fof_number_].Weight);
      MCMC_FOF[output_number_][fof_number_].Weight /= BoxSize*BoxSize*BoxSize;
#ifdef HALOMODEL
      MCMC_FOF[output_number_][fof_number_].NGalsInFoF=0;
      MCMC_FOF[output_number_][fof_number_].IndexOfCentralGal=-1;
#endif
    }

    fclose(file_);

    qsort(MCMC_FOF[output_number_], NFofsInSample[output_number_], sizeof(struct MCMC_FOF_struct), MCMC_FOF_compare_FoFID);

    // // debugging: checks ordering by fof id (order would ease work later):
    // for(fof_number_ = 1; fof_number_ < NFofsInSample[output_number_]; fof_number_++)
    // {
    //   if(!(MCMC_FOF[output_number_][fof_number_ - 1].FoFID < MCMC_FOF[output_number_][fof_number_].FoFID))
    //   {
    //     printf("fof ids not ordered in file `%s': output_number_ = %d, fof_number_ = %d, !(%lld < %lld)\n",
    //       file_name_, output_number_, fof_number_, MCMC_FOF[output_number_][fof_number_ - 1].FoFID, MCMC_FOF[output_number_][fof_number_].FoFID);        
    //     // sprintf(error_message_, "fof ids not ordered in file `%s': output_number_ = %d, fof_number_ = %d, !(%lld < %lld)\n",
    //     //   file_name_, output_number_, fof_number_, MCMC_FOF[output_number_][fof_number_ - 1].FoFID, MCMC_FOF[output_number_][fof_number_].FoFID);
    //     // terminate(error_message_);
    //     break;
    //   }
    // }
  }
}


/** @brief Read in the arrays of observational data. They will be compared
 *       with the outputs from the SAM On the get_likelihood routine*/
void read_observations (void)
{
  int constraint_number_, bin_number_, output_number_, mcmc_constraints_output_number_, number_of_tests, n_mcmc_constraints_outputs_;
  //FILE *f[MCMCNConstraints];
  FILE *file_;
  float BinValueColector_;
  char buf_[1000], error_message_[1000], aux_test_name_[1000], aux_test_type_[1000];
  bool match_found_;

  //allocate structure to contain observational data
  MCMC_Obs = mymalloc("MCMC_Obs", sizeof(struct MCMC_OBSCONSTRAINTS) * MCMCNConstraints);

  /*Read file MCMCObsConstraints.txt that contains a header with
   * number_of_tests, number_of_chi_tests, number_of_binom_tests,
   * n_mcmc_constraints_outputs_ and the redshifts that will be used to constrain the MCMC
   * Then there is a list of all the observational constraints, with
   * the type of test to be used and a switch controlling which redshifts
   * will be used for each constraint */
  sprintf(buf_, "%s", MCMCObsConstraints);
  if(!(file_ = fopen(buf_, "r")))
  {
    sprintf(error_message_, "can't open file `%s'\n", buf_);
    terminate(error_message_);
  }

  fgets(buf_, 500, file_);
  fscanf(file_,"%d\n",&number_of_tests);

  if(number_of_tests!=MCMCNConstraints)
  {
    sprintf(error_message_, "check MCMCNConstraints, NChiTests & NBinomTests\n");
    terminate(error_message_);
  }

  fgets(buf_, 500, file_);
  fgets(buf_, 500, file_);
  fscanf(file_,"%d\n",&n_mcmc_constraints_outputs_);
  if(n_mcmc_constraints_outputs_ > NOUT)
  {
    sprintf(error_message_, "n_mcmc_constraints_outputs_ > NOUT\n");
    terminate(error_message_);
  }

   //Check if all the required snaps are in the current ListOutputSnaps
  for(mcmc_constraints_output_number_ = 0; mcmc_constraints_output_number_ < n_mcmc_constraints_outputs_; mcmc_constraints_output_number_++)
  {
    fscanf(file_,"%lg\n",&MCMCConstraintsZZ[mcmc_constraints_output_number_]);
    match_found_ = false;
    for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    {
      if(MCMCConstraintsZZ[mcmc_constraints_output_number_] >= (double)((int)((ZZ[ListOutputSnaps[output_number_]]*10)+0.5)/10.)-0.1 &&
         MCMCConstraintsZZ[mcmc_constraints_output_number_] <= (double)((int)((ZZ[ListOutputSnaps[output_number_]]*10)+0.5)/10.)+0.1)
      {
        match_found_ = true;
        break;
      }
    }
    if(!match_found_)
    {
      sprintf(error_message_, "redshift %0.2f required for MCMC not in outputlist \n", MCMCConstraintsZZ[mcmc_constraints_output_number_]);
      terminate(error_message_);
    }
  }

  fgets(buf_, 500, file_);

  //Scan test names, types and redshift switches
  for(constraint_number_ = 0; constraint_number_ < MCMCNConstraints; constraint_number_++)
  {
    fscanf(file_,"%s %s",MCMC_Obs[constraint_number_].Name, MCMC_Obs[constraint_number_].TestType);
    for(output_number_=0;output_number_<NOUT;output_number_++)
      fscanf(file_,"%d",&MCMC_Obs[constraint_number_].ObsTest_Switch_z[output_number_]);
  }
  fclose(file_);

  //now read weights for different constraints
  sprintf(buf_, "%s", MCMCWeightsObsConstraints);
  if(!(file_ = fopen(buf_, "r")))
  {
    sprintf(error_message_, "can't open file `%s'\n", buf_);
    terminate(error_message_);
  }

  fgets(buf_, 500, file_);
  fscanf(file_,"%d\n",&number_of_tests);

  if(number_of_tests!=MCMCNConstraints)
  {
    sprintf(error_message_, "check MCMCNConstraints, NChiTests & NBinomTests\n");
    terminate(error_message_);
  }
  fgets(buf_, 500, file_);
  for(constraint_number_=0;constraint_number_<MCMCNConstraints;constraint_number_++)
  {
    fscanf(file_,"%s %s",aux_test_name_, aux_test_type_);
    for(output_number_=0;output_number_<NOUT;output_number_++)
      fscanf(file_,"%lf",&MCMC_Obs[constraint_number_].ObsTest_Weight_z[output_number_]);
  }
  fclose(file_);

  //now read the observations
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    //the round of ZZ[] is to assure that independently of the cosmology used you still
    //round to z=1.0 or 2.0,etc...
    for(constraint_number_=0;constraint_number_<MCMCNConstraints;++constraint_number_)
    {
      if(MCMC_Obs[constraint_number_].ObsTest_Switch_z[output_number_]==1)
      {
        sprintf(buf_, "%s/%s_z%1.2f.txt",ObsConstraintsDir,MCMC_Obs[constraint_number_].Name,
                (double)((int)((MCMCConstraintsZZ[output_number_]*10)+0.5)/10.) );
        if((file_=fopen(buf_,"r"))==NULL)
        {
          sprintf(error_message_, "can't open file `%s'\n", buf_);
          terminate(error_message_);
        }

        /* This values will give the size of the observational arrays
         * They will be used in get_likelihood. */
        fscanf(file_, "%d", &Nbins[output_number_][constraint_number_]);
        if(Nbins[output_number_][constraint_number_]>MCMCMaxObsBins)
        {
          sprintf(error_message_, "NBins for Obs[%d] Snap[%d] > MCMCMaxObsBins", constraint_number_,output_number_);
          terminate(error_message_);
        }

        //Read observational data
        for(bin_number_ = 0; bin_number_ < Nbins[output_number_][constraint_number_]; bin_number_++)
        {
          //Chi_Sq and Maximum Likelihood TESTS
          if(strcmp(MCMC_Obs[constraint_number_].TestType,"chi_sq")==0 || strcmp(MCMC_Obs[constraint_number_].TestType,"maxlike")==0)
            fscanf(file_, "%lg %lg %lg %lg", &MCMC_Obs[constraint_number_].Bin_low[output_number_][bin_number_], &MCMC_Obs[constraint_number_].Bin_high[output_number_][bin_number_],
                   &MCMC_Obs[constraint_number_].Obs[output_number_][bin_number_], &MCMC_Obs[constraint_number_].Error[output_number_][bin_number_]);
          //Binomial TESTS
          else if(strcmp(MCMC_Obs[constraint_number_].TestType,"binomial")==0)
            fscanf(file_, "%f %lg %lg", &BinValueColector_, &MCMC_Obs[constraint_number_].ObsUp[output_number_][bin_number_], &MCMC_Obs[constraint_number_].ObsDown[output_number_][bin_number_]);
        }
        fclose(file_);

      }//end if(MCMC_Obs[constraint_number_].ObsTest_Switch_z[output_number_]==1)
    }//end loop on tests
  }//end loop on snaps
}


/** @brief prepares files to contain comparison to observations 
 * 
 * A different file is written for each observational constraint and for each redshift
 * Each file as the following structure: *
 * int ChainLength+1
 * int Nbins[snapshot_number_][constraint]
 * phi[0],phi[1],phi[Nbins[snapshot_number_][constraint]],likelihood(of this constraint),likelihood(of this step)
 *
 * Plus one file to contain the total likelihood at each step*/
void open_files_with_comparison_to_observations()
{
  int constraint_number_, snapshot_number_;
  char file_name_[1000], error_message_[1000];

  sprintf(file_name_, "%sMCMC_LIKELIHOOD_%d.txt", OutputDir, ThisTask+FirstChainNumber);
  if((FILE_MCMC_LIKELIHOOD = fopen(file_name_, "w")) == NULL)
  {
    sprintf(error_message_, "can't open file `%s'\n", file_name_);
    terminate(error_message_);
  }

  for(constraint_number_ = 0; constraint_number_ < MCMCNConstraints; constraint_number_++)
    for(snapshot_number_=0;snapshot_number_<NOUT;snapshot_number_++)
    {
      if(MCMC_Obs[constraint_number_].ObsTest_Switch_z[snapshot_number_]==1)
      {
        sprintf(file_name_, "%sPredictionsPerStep_%s_z%1.2f_%d.txt", OutputDir, MCMC_Obs[constraint_number_].Name,
                (double)((int)((MCMCConstraintsZZ[snapshot_number_]*10)+0.5)/10.), ThisTask+FirstChainNumber);
        if((FILE_MCMC_PredictionsPerStep[snapshot_number_][constraint_number_] = fopen(file_name_, "w")) == NULL)
        {
          sprintf(error_message_, "can't open file `%s'\n", file_name_);
          terminate(error_message_);
        }
        fprintf(FILE_MCMC_PredictionsPerStep[snapshot_number_][constraint_number_], " %d\n", ChainLength+1);
        fprintf(FILE_MCMC_PredictionsPerStep[snapshot_number_][constraint_number_], " %d\n", Nbins[snapshot_number_][constraint_number_]);
      }
    }
}


/** @brief close files with comparison to observations */
void close_files_with_comparison_to_observations()
{
  int constraint_number_, snapshot_number_;

  fclose(FILE_MCMC_LIKELIHOOD);

  for(constraint_number_ = 0; constraint_number_ < MCMCNConstraints; constraint_number_++)
    for(snapshot_number_=0;snapshot_number_<NOUT;snapshot_number_++)
      if(MCMC_Obs[constraint_number_].ObsTest_Switch_z[snapshot_number_]==1)
        fclose(FILE_MCMC_PredictionsPerStep[snapshot_number_][constraint_number_]);
}


#ifdef MR_PLUS_MRII
/** @brief change dark matter simulation */
void change_dark_matter_sim(const char SimName[])
{
  if (strcmp(SimName,"MR")==0)
  {
    Switch_MR_MRII=1;
    sprintf(FileWithZList, "%s", FileWithZList_MR);
    PartMass=PartMass_MR;
    BoxSize=BoxSize_MR;

    sprintf(FileWithZList_OriginalCosm, "%s", FileWithZList_OriginalCosm_MR);
    PartMass_OriginalCosm=PartMass_OriginalCosm_MR;
    BoxSize_OriginalCosm=BoxSize_OriginalCosm_MR;

    LastDarkMatterSnapShot=LastDarkMatterSnapShot_MR;
  }
  else if (strcmp(SimName,"MRII")==0)
  {
    Switch_MR_MRII=2;
    sprintf(FileWithZList, "%s", FileWithZList_MRII);
    PartMass=PartMass_MRII;
    BoxSize=BoxSize_MRII;

    sprintf(FileWithZList_OriginalCosm, "%s", FileWithZList_OriginalCosm_MRII);
    PartMass_OriginalCosm=PartMass_OriginalCosm_MRII;
    BoxSize_OriginalCosm=BoxSize_OriginalCosm_MRII;

    LastDarkMatterSnapShot=LastDarkMatterSnapShot_MRII;
  }

  //do part of init() again.
  //Can't do the all functions because of metallicity arrays in read_cooling_functions()
  read_zlist();
  read_zlist_original_cosm();
  read_output_snaps();

  //CREATE ARRAYS OF SFH TIME STRUCTURE:
#ifdef  STAR_FORMATION_HISTORY
  create_sfh_bins();
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
  //read in photometric tables
#ifdef PHOTTABLES_PRECOMPUTED
  setup_LumTables_precomputed(SimName);
#endif
#ifdef SPEC_PHOTABLES_ON_THE_FLY
  setup_Spec_LumTables_onthefly();
#endif
#endif

#ifdef DETAILED_METALS_AND_MASS_RETURN
  init_integrated_yields();
  integrate_yields();
#endif

  /*After all the millennium trees have been done read sample for MRII trees
   * from treenr=NTrees_MR to NTrees_MR+NTrees_MRII=Ntrees */
  if (strcmp(SimName,"MR")==0)
  {
    sprintf(MCMCSampleFilePrefix,"%s",MCMCSampleFilePrefix_MR);
    MCMCSampleFile=MCMCSampleFile_MR;
  }
  else if (strcmp(SimName,"MRII")==0)
  {
    sprintf(MCMCSampleFilePrefix,"%s",MCMCSampleFilePrefix_MRII);
    MCMCSampleFile=MCMCSampleFile_MRII;
  }

  if (strcmp(SimName,"MRII")==0)
  { free_MCMC_FOF(); }

  read_sample_info();
}
#endif


#ifdef HALOMODEL
/** @brief asign FOF halo masses */
void assign_FOF_masses(const int snapshot_number_, const int tree_number_)
{
  int fof_group_number_, output_number_, halo_number_;

  for(halo_number_ = 0; halo_number_ < TreeNHalos[tree_number_]; halo_number_++)
    if(HaloAux[halo_number_].DoneFlag == 0 && Halo[halo_number_].SnapNum == snapshot_number_)
    {
      for(output_number_=0;output_number_<NOUT;output_number_++)
      {
        if(snapshot_number_==ListOutputSnaps[output_number_])
        {
          for(fof_group_number_=0;fof_group_number_<NFofsInSample[output_number_]; fof_group_number_++)
            if(HaloIDs[halo_number_].FirstHaloInFOFgroup == MCMC_FOF[output_number_][fof_group_number_].FoFID)
            {
              MCMC_FOF[output_number_][fof_group_number_].M_Crit200 = log10(Halo[halo_number_].M_Crit200*1.e10);
              MCMC_FOF[output_number_][fof_group_number_].M_Mean200 = log10(Halo[halo_number_].M_Mean200*1.e10);
#ifdef MCRIT
              MCMC_FOF[output_number_][fof_group_number_].M_Mean200 = log10(Halo[halo_number_].M_Crit200*1.e10);
#endif
            }
        }
      }
    }
}
#endif


/**@brief marks halos that will be used in MCMC sampling 
 *         and (if halo model is used) copies FOF masses 
 *         (because then we need halo masses even for FOFs
 *          with no galaxies)
 *
 * - requires MCMC_FOF[].FoFID[], NFofsInSample[] to be up to date, 
 *   which is currently read from file in read_sample_info().
 * - requires MCMC_FOF[].FoFID[] ordered for binary search.
 * - requires Halo, HaloIDs, HaloAux, TreeNHalos to be up to date,
 *   which are currently allocated/read in load_tree(tree_number_).
 */
void link_halos_and_MCMC_FOF(const int tree_number_)
{
  int halo_number_, output_number_, fof_number_ /*, fof_number_lower_bound_, fof_number_upper_bound_*/;
  
  // printf("marking halos for tree %d.\n", tree_number_);
  
  for(halo_number_ = 0; halo_number_ < TreeNHalos[tree_number_]; halo_number_++)
  {
    const long long fof_id_of_halo_ = HaloIDs[halo_number_].FirstHaloInFOFgroup;
    
    HaloAux[halo_number_].halo_is_in_MCMC_sample_for_any_output = false;

    for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    {
      HaloAux[halo_number_].fof_number_in_MCMC_sample_for_output[output_number_] = -1;
      HaloAux[halo_number_].halo_is_in_MCMC_sample_for_output   [output_number_] = false;

      //binary search:
      if((Halo[halo_number_].SnapNum == ListOutputSnaps[output_number_]) && (NFofsInSample[output_number_] > 0) && (MCMC_FOF[output_number_][0].FoFID <= fof_id_of_halo_) && (fof_id_of_halo_ <= MCMC_FOF[output_number_][NFofsInSample[output_number_] - 1].FoFID))
      {
        unsigned int fof_number_lower_bound_ = 0;
        unsigned int fof_number_upper_bound_ = NFofsInSample[output_number_] - 1;
        while(fof_number_ = (fof_number_lower_bound_ + fof_number_upper_bound_) / 2, fof_number_lower_bound_ < fof_number_upper_bound_)
        {
          if(MCMC_FOF[output_number_][fof_number_].FoFID < fof_id_of_halo_)
            fof_number_lower_bound_ = fof_number_ + 1;
          else
            fof_number_upper_bound_ = fof_number_;
        }
        if(fof_id_of_halo_ == MCMC_FOF[output_number_][fof_number_].FoFID)
        {
          HaloAux[halo_number_].fof_number_in_MCMC_sample_for_output[output_number_] = fof_number_;
          HaloAux[halo_number_].halo_is_in_MCMC_sample_for_output   [output_number_] = true;
          HaloAux[halo_number_].halo_is_in_MCMC_sample_for_any_output                = true;
          
#ifdef HALOMODEL
          MCMC_FOF[output_number_][fof_number_].M_Crit200 = log10(Halo[halo_number_].M_Crit200*1.e10);
          MCMC_FOF[output_number_][fof_number_].M_Mean200 = log10(Halo[halo_number_].M_Mean200*1.e10);
#ifdef MCRIT
          MCMC_FOF[output_number_][fof_number_].M_Mean200 = log10(Halo[halo_number_].M_Crit200*1.e10);
#endif /* defined MCRIT */
#endif /* defined HALOMODEL */
        }
      }
    }
  }
  // printf("done marking halos for tree %d.\n", tree_number_);
}


/** @brief free memory for MCMC_FOF **/
void free_MCMC_FOF (void)
{
  unsigned int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  { free(MCMC_FOF[output_number_]); }
}


/** @brief compares  for sorting MCMC_FOF_*/
int MCMC_FOF_compare_FoFID(const void *MCMC_FOF_a_, const void *MCMC_FOF_b_)
{
  if     (((struct MCMC_FOF_struct*)MCMC_FOF_a_)->FoFID < ((struct MCMC_FOF_struct*)MCMC_FOF_b_)->FoFID)
    return -1;   
  else if(((struct MCMC_FOF_struct*)MCMC_FOF_a_)->FoFID > ((struct MCMC_FOF_struct*)MCMC_FOF_b_)->FoFID)
    return +1;
  else
    return 0;
}

#endif /* defined MCMC */
