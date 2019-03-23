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

/** @file   model_spectro_photometric.c
 *  @date   2009-2019
 *  @author Chiara Tonini
 *  @author Bruno Henriques
 *  @author Stefan Hilbert 
 *
 *  @brief  spectral and photometric properties of galaxies
 **/

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
void setup_LumTables_precomputed(const char sim_name_[])
{
  FilterLambda[NMAG] = 0.55;        //to use by the dust model for birth clouds, the wavelength of the V-filter_number_
  
  FILE *metallicity_list_file_, *filter_file_;
  int met_index_, age_index_, filter_number_, snapshot_number_;
  char file_name_[1000], filter_name_[100], dummy_[100], SSP_name_[100];
  char dumb_filter_file_[100];
  float dumb_filter_lambda_;
  int dumb_ssp_n_snapshots_, dumb_ssp_n_ages_, dumb_ssp_n_metallicites_, dumb_n_mags_;

#ifdef BC03
  sprintf(SSP_name_, "BC03");
#endif
#ifdef M05
  sprintf(SSP_name_, "M05");
#endif
#ifdef CB07
  sprintf(SSP_name_, "CB07");
#endif

  /*Read list of metallicities available from SSP_name_*/
  sprintf(file_name_, "%s/PhotTables/%s_%s_Metallicity_list.dat", SpecPhotDir, SSP_name_, SpecPhotIMF);
  if(!(metallicity_list_file_ = fopen(file_name_, "r")))
  {
    char error_message_[2048];
    sprintf(error_message_, "file `%s' not found.\n", file_name_);
    terminate(error_message_);
  }

  fscanf(metallicity_list_file_, "%d", &dumb_ssp_n_metallicites_);
  if(dumb_ssp_n_metallicites_ != SSP_NMETALLICITES)
  {
    terminate("nmetallicites on file not equal to SSP_NMETALLICITES");
  }

  for(met_index_=0;met_index_<SSP_NMETALLICITES;met_index_++)
  {
    fscanf(metallicity_list_file_, "%f", &SSP_logMetalTab[met_index_]);
    SSP_logMetalTab[met_index_]=log10(SSP_logMetalTab[met_index_]);
  }
  fclose(metallicity_list_file_);

  /*Loop over the different files corresponding to different metallicities */
  for(met_index_ = 0; met_index_ < SSP_NMETALLICITES; met_index_++)
  {
    sprintf(file_name_, "%s", FileWithFilterNames);
    if((metallicity_list_file_ = fopen(file_name_, "r")) == NULL)
    {
      printf("\n**Can't open file \"%s\" **\n", file_name_);
      char error_message_[2048];
      sprintf(error_message_, "Can't open file %s\n", file_name_);
      terminate(error_message_);
    }

    fscanf(metallicity_list_file_, "%d", &dumb_n_mags_);
    if(dumb_n_mags_ != NMAG)
    {
      char error_message_[2048];
      sprintf(error_message_,"nmag = %d on file %s not equal to NMAG = %d",dumb_n_mags_, file_name_, NMAG);
      terminate(error_message_);
    }

    //There is a different file for each filter_number_
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      fscanf(metallicity_list_file_,"%s %f %s" ,dumb_filter_file_, &dumb_filter_lambda_, filter_name_);
      //READ TABLES
      sprintf(file_name_, "%s/PhotTables/%s_%s_Phot_Table_%s_Mag%s_m%0.4f.dat", SpecPhotDir, PhotPrefix, SpecPhotIMF, sim_name_, filter_name_,
              pow(10,SSP_logMetalTab[met_index_]));
      if(!(filter_file_ = fopen(file_name_, "r")))
      {
        char error_message_[2048];
        sprintf(error_message_, "file `%s' not found.\n", file_name_);
        terminate(error_message_);
      }

      if(ThisTask == 0)
        printf("reading file %s \n", file_name_);

      fscanf(filter_file_, "%s %f %d %d", dummy_, &FilterLambda[filter_number_], &dumb_ssp_n_snapshots_, &dumb_ssp_n_ages_);
      /* check that the numbers on top of the file correspond to (LastDarkMatterSnapShot+1) and SSP_NAGES */
      if(FilterLambda[filter_number_] != dumb_filter_lambda_)
      {
        terminate("filterlambdas dont match");
      }
      if(dumb_ssp_n_snapshots_ != (LastDarkMatterSnapShot+1))
      {
        terminate("nsnaps not equal to Snaplistlen");
      }
      if(dumb_ssp_n_ages_ != SSP_NAGES)
      {
        terminate("n_age_bins not equal to SSP_NAGES");
      }

      /* read ages of SSPs (done for all filter_number_s and metallicities)
        * last call stays on global array  SSP_logAgeTab*/
      for(age_index_ = 0; age_index_ < SSP_NAGES; age_index_++)
      {
        fscanf(filter_file_, " %e ", &SSP_logAgeTab[age_index_]);

        if(SSP_logAgeTab[age_index_] > 0.0)        // avoid taking a log of 0 ...
        {
          /* converts SSP_AgeTab from years to log10(internal time units) */
          SSP_logAgeTab[age_index_] = SSP_logAgeTab[age_index_] / UnitTime_in_years * Hubble_h;
          SSP_logAgeTab[age_index_] = log10(SSP_logAgeTab[age_index_]);
        }
        else
          SSP_logAgeTab[age_index_] = 0.;
      }

      //read luminosities at each output redshift
      for(snapshot_number_ = 0; snapshot_number_ < (LastDarkMatterSnapShot+1); snapshot_number_++)
      {
        fscanf(filter_file_, " %f ", &RedshiftTab[snapshot_number_]);
        //for each age
        for(age_index_ = 0; age_index_ < SSP_NAGES; age_index_++)
        {
          fscanf(filter_file_, "%e", &LumTables[age_index_][met_index_][snapshot_number_][filter_number_]);
          LumTables[age_index_][met_index_][snapshot_number_][filter_number_] = pow(10., -0.4 * LumTables[age_index_][met_index_][snapshot_number_][filter_number_]);
        }                //end loop on age
      }                //end loop on redshift (everything done for current filter_number_)

      fclose(filter_file_);
    }//end loop on filter_number_s

    fclose(metallicity_list_file_);
  }//end loop on metallicity

  init_SSP_log_age_jump_index();
}
#endif /* defined PHOTTABLES_PRECOMPUTED */


/** @brief init table for shortcut lookup into LumTables */
void init_SSP_log_age_jump_index(void)
{
  double log_age_;
  int i_, idx_;

  SSP_log_age_jump_factor = SSP_NJUMPTAB / (SSP_logAgeTab[SSP_NAGES - 1] - SSP_logAgeTab[1]);

  for(i_ = 0; i_ < SSP_NJUMPTAB; i_++)
  {
    log_age_ = SSP_logAgeTab[1] + i_ / SSP_log_age_jump_factor;
    idx_ = 1; while(SSP_logAgeTab[idx_ + 1] < log_age_) { ++idx_; };
    SSP_log_age_jump_table[i_] = idx_;
  }
}


/** @brief Whenever star formation occurs, calculates the luminosity corresponding
  *        to the mass of stars formed, considering the metallicity and age of the
  *        material.
  *
  * The semi-analytic code uses look up tables produced by Evolutionary Population
  * Synthesis Models to convert the mass formed on every star formation episode
  * into a luminosity. Each of These tables corresponds to a simple stellar
  * population i_.e, a population with a single metallicity. For a given IMF,
  * metatillicty and age, the tables give the luminosity for a
  * \f$ 10^{11}M_\odot\f$ burst. The default model uses a Chabrier IMF and
  * stellar populations from Bruzual & Charlot 2003 with 6 different metallicites.
  *
  * The magnitudes are immediately calculated for each output bin, so that we know
  * the age of each population that contributed to a galaxy total population: the
  * age between creation and output. Apart from the different ages of the populations
  * at a given output bin, if the option OUTPUT_OBS_MAGS is turned on, then we also
  * need to know the K-corrections (going the opposite directions as in observations)
  * that will affect each population.
  *
  * For each metallicity there is a look up table which has the different magnitudes
  * for each age and then this is k-corrected to all the snapshots.
  *
  * If MetallicityOption = 0 -> only solar metallicity.
  * If MetallicityOption = 1 -> 6 metallicities.
  *
  * @bug (corrected by Stefan Hilbert) 
  *      MetallicityOption = 0 (-> only solar metallicity) used to only set tabindex,
  *      but not fractions, now corrected (by Stefan Hilbert)
  **/
#ifndef  POST_PROCESS_MAGS
void add_to_luminosities(const int galaxy_number_, double stellar_mass_, double time_, double dt_, const double metallicity_)
{
#ifndef OUTPUT_REST_MAGS
#ifndef OUTPUT_OBS_MAGS
  return; /* early return if no mags to compute */
#endif /* not defined OUTPUT_OBS_MAGS */
#endif /* not defined OUTPUT_REST_MAGS */

  int output_number_, filter_number_;
  double luminosity_to_add_;
   
  int age_index_;  double f_age_1_, f_age_2_; 
  int met_index_;  double f_met_1_, f_met_2_;
  int redshift_index_;

  //if one wants to have finner bins for the star formation then the STEPS
  //of the calculation, N_FINE_AGE_BINS should be set to > 1
  
#if !((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
  (void)dt_; /* suppress unused-parameter warning */
#endif /* not defined N_FINE_AGE_BINS > 1 */

  /* Time below which the luminosities are corrected for extinction due to
   * molecular birth clouds.  */
  const double birthcloud_age_ = 10.0 / UnitTime_in_Megayears * Hubble_h;

  /* mstars converted from 1.e10Msun/h to 1.e11 Msun */
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
  stellar_mass_ *= 0.1 / N_FINE_AGE_BINS * inv_Hubble_h;
#else  /* not defined FINE_AGE_BINS > 1 */
  stellar_mass_ *= 0.1 * inv_Hubble_h;
#endif /* not defined FINE_AGE_BINS > 1 */

  /* now we have to change the luminosities accordingly. */
  /* note: we already know at which place we have to look up the tables,
   * since we know the output times, the current time_ and the metallicity.
   * find_interpolated_lum() finds the 2 closest points in the SPS table
   * in terms of age and metallicity. Time gives the time_to_present for
   * the current step while NumToTime(ListOutputSnaps[output_number_]) gives
   * the time_ of the output snapshot_number_ - units Mpc/Km/s/h */
  
  if(MetallicityOption == 0) // reset met index to use only solar metallicity
  { met_index_ = 4; f_met_1_ = 1., f_met_2_ = 0.; } 
  else if(metallicity_ <= 0.)
  { met_index_ = 0; f_met_1_ = 1., f_met_2_ = 0.; } 
  else   
  {
    const double log10_metallicity_ = log10(metallicity_);
    find_metallicity_luminosity_interpolation_parameters(log10_metallicity_, met_index_, f_met_1_, f_met_2_);
  }
  
#if ((defined N_FINE_AGE_BINS) && (N_FINE_AGE_BINS > 1))
  const double lower_time_ = time_ - 0.5 * dt;
  dt_  /= N_FINE_AGE_BINS;
  int fine_age_step_;
  for(fine_age_step_ = 0; fine_age_step_< N_FINE_AGE_BINS; fine_age_step_++)
  {
    time_ = lower_time_ + (fine_age_step_ + 0.5) * dt_;
#endif /* defined N_FINE_AGE_BINS > 1 */

#ifdef GALAXYTREE
    const int output_number_beg_ = Gal[galaxy_number_].SnapNum;
    const int output_number_end_ = NOUT;
#else  /* not defined GALAXYTREE */
    const int output_number_beg_ = 0;
    const int output_number_end_ = ListOutputNumberOfSnapshot[Gal[galaxy_number_].SnapNum] + 1;
#endif /* not defined GALAXYTREE */

    for(output_number_ = output_number_beg_; output_number_ < output_number_end_; output_number_++)
    {
      const double age_                        = time_ - NumToTime(ListOutputSnaps[output_number_]);
      
      if(age_ <= 0) continue;
      
      const double log10_age_                  = log10(age_);
      const bool   is_affected_by_birthclould_ = (age_ <= birthcloud_age_);
      find_age_luminosity_interpolation_parameters(log10_age_, age_index_, f_age_1_, f_age_2_);
      
#ifdef OUTPUT_REST_MAGS
      /* For rest-frame, there is no K-correction on magnitudes,
       * hence the 0 in LumTables[filter_number_][met_index_][0][age_index_] */
      for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
      {
        //interpolation between the points found by find_interpolated_lum
        luminosity_to_add_ = stellar_mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][0][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_    ][0][filter_number_]) +
                                              f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][0][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][0][filter_number_]));
                                         
        Gal[galaxy_number_].Lum[output_number_][filter_number_] += luminosity_to_add_;

        /*luminosity used for extinction due to young birth clouds */
        if(is_affected_by_birthclould_)
          Gal[galaxy_number_].LumY[output_number_][filter_number_] += luminosity_to_add_;
        
        if(filter_number_ == R_BAND_FILTER_NUMBER)
          Gal[galaxy_number_].rbandWeightAge[output_number_] += age_ * luminosity_to_add_;
      }
#endif /* defined OUTPUT_REST_MAGS */

#ifdef OUTPUT_OBS_MAGS
      redshift_index_ = LastDarkMatterSnapShot - ListOutputSnaps[output_number_];

      /* Note the zindex in LumTables[][][][] meaning the magnitudes are now
        * "inversely k-corrected to get observed frame at output bins" */
      for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
      {
        luminosity_to_add_ = stellar_mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][redshift_index_][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_    ][redshift_index_][filter_number_]) +
                                              f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][redshift_index_][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][redshift_index_][filter_number_]));
                                                       
        Gal[galaxy_number_].ObsLum[output_number_][filter_number_] += luminosity_to_add_;
  
        if(is_affected_by_birthclould_)
          Gal[galaxy_number_].ObsLumY[output_number_][filter_number_] += luminosity_to_add_;
      }
      
#ifdef OUTPUT_FB_OBS_MAGS
      redshift_index_ = (ListOutputSnaps[output_number_] > 0 ? LastDarkMatterSnapShot - ListOutputSnaps[output_number_] - 1 : LastDarkMatterSnapShot);
      for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
      {
        luminosity_to_add_ = stellar_mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][redshift_index_][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_    ][redshift_index_][filter_number_]) +
                                              f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][redshift_index_][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][redshift_index_][filter_number_]));
                                                       
        Gal[galaxy_number_].backward_ObsLum[output_number_][filter_number_] += luminosity_to_add_;

        if(is_affected_by_birthclould_)
          Gal[galaxy_number_].backward_ObsLumY[output_number_][filter_number_] += luminosity_to_add_;
      }
      
      redshift_index_ = (ListOutputSnaps[output_number_] < LastDarkMatterSnapShot ? LastDarkMatterSnapShot - ListOutputSnaps[output_number_] + 1 : 0);
      for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
      {
        luminosity_to_add_ = stellar_mass_ * (f_met_1_ * (f_age_1_ * LumTables[age_index_    ][met_index_    ][redshift_index_][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_    ][redshift_index_][filter_number_]) +
                                              f_met_2_ * (f_age_1_ * LumTables[age_index_    ][met_index_ + 1][redshift_index_][filter_number_]  +
                                                          f_age_2_ * LumTables[age_index_ + 1][met_index_ + 1][redshift_index_][filter_number_]));
                                                       
        Gal[galaxy_number_].forward_ObsLum[output_number_][filter_number_] += luminosity_to_add_;

        if(is_affected_by_birthclould_)
          Gal[galaxy_number_].forward_ObsLumY[output_number_][filter_number_] += luminosity_to_add_;
      }
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
    }

#ifdef FINE_AGE_BINS  
  }//end loop on small age bins
#endif /* defined FINE_AGE_BINS */
}
#endif  //POST_PROCESS_MAGS

#endif // COMPUTE_SPECPHOT_PROPERTIES

