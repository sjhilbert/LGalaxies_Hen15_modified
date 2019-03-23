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

/** @file   main.c
 *  @date   1978-2019
 *  @author Gabriella De Lucia, Qi Guo, Bruno Henriques, Stefan Hilbert, Guinevere Kauffmann,
 *          Volker Springel, Peter Thomas, Marcel van Daalen, Simon White, Rob Yates, et al.
 *
 *  @brief just the main function of L-Galaxies
 *
 *  Created originally in 1978
 *  Re-written in 2001 by Volker Springel
 * 
 *  Major Contributors: Qi Guo, Gabriella De Lucia, Bruno Henriques, Guinevere Kauffmann,
 *                      Volker Springel, Peter Thomas, Marcel van Daalen, Rob Yates, Simon White
 *
 *  Contributors:       Stefan Hilbert, Rachel Asquith
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#ifdef PARALLEL
#include <mpi.h>
#endif /* defined PARALLEL */

#include "allvars.h"
#include "proto.h"
#include "check_compile_time_options.h"


#ifdef MCMC
#include "mcmc_vars.h"
#include "mcmc_proto.h"
#endif /* defined MCMC */


/** @file main.c
 * @brief Controlling function of L-Galaxies plus SAM Construct Galaxies,
 *        Join Galaxies of progenitors, Evolve Galaxies and check
 *        Makefile options.
 *
 * Reads the parameter file given as an argument at run time:
 * read_parameter_file(argv_[1]).
 *
 * Checks for consistency between some Makefile options: check_options().
 *
 * Runs init() to initialize some stuff.
 *
 * Ifdef MCMC - it calls Senna which will start the MCMC parameter estimation,
 * calling SAM for each step of the chain to run the semi-analytic model on
 * a set of merger trees and compute a likelihood for the model with those
 * parameters with respect to a set of observations - see mcmc.c
 *
 * Otherwise, SAM is called for each of the chosen dark matter tree files and
 * for each tree, on each treefile: reads tree, constructs the galaxies, saves
 * the galaxies, frees memory for galaxies and tree.
 *
 * @returns       0
 * */

/**@brief Main routine of L-Galaxies*/
int main(int argc_, char **argv_)
{
#ifndef MCMC
  char buffer_[1000];
  time_t start_time_, current_time_;
#endif /* not defined MCMC */

#ifdef PARALLEL
  MPI_Init(&argc_, &argv_);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#endif /* defined PARALLEL */

#ifdef MCMC
  time(&GlobalStartingTime);
#endif /* defined MCMC */

  if(ThisTask==0)
  {
    printf("\n\n\n");
    printf("**************************************************************************\n");
    printf("*                                                                        *\n");
    printf("*                   Copyright (C) <2016+>  <L-Galaxies>                  *\n");
    printf("*                                                                        *\n");
    printf("*  This program is free software: you can redistribute it and/or modify  *\n");
    printf("*  it under the terms of the GNU General Public License as published by  *\n");
    printf("*  the Free Software Foundation, either version 3 of the License, or     *\n");
    printf("*  (at your option) any later version.                                   *\n");
    printf("*                                                                        *\n");
    printf("*  This program is distributed in the hope that it will be useful,       *\n");
    printf("*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *\n");
    printf("*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *\n");
    printf("*  GNU General Public License for more details.                          *\n");
    printf("*                                                                        *\n");
    printf("*  You should have received a copy of the GNU General Public License     *\n");
    printf("*  along with this program.  If not, see <http://www.gnu.org/licenses/>  *\n");
    printf("*                                                                        *\n");
    printf("**************************************************************************\n\n\n");
  }

  if(2 > argc_ || argc_ > 3)
  {
    printf("\n  Wrong number of runtime arguments\n\n");
    printf("\n  usage: ./L-Galaxies <parameterfile>\n\n");
    endrun(0);
  }

 if (ThisTask == 0)
 { printf("%s\n",COMPILETIMESETTINGS); }

 /* check compatibility of some Makefile Options*/
  check_compile_time_options();

  /*Reads the parameter file, given as an argument at run time. */
  read_parameter_file(argv_[1]);
  check_program_parameters();
  compute_derived_program_parameters();

#ifdef MR_PLUS_MRII
  //Start with MR files and later change to MRII
  LastDarkMatterSnapShot=LastDarkMatterSnapShot_MR;
  sprintf(FileWithZList, "%s", FileWithZList_MR);
  sprintf(FileWithZList_OriginalCosm, "%s", FileWithZList_OriginalCosm_MR);
#endif /* defined MR_PLUS_MRII */

  mymalloc_init();

  sprintf(FinalOutputDir, "%s", OutputDir);
#ifndef MCMC
  if(argc_ == 3)
  { sprintf(OutputDir, "%s", argv_[2]); }
#else /* defined MCMC */
  FirstChainNumber=0;
  if(argc_ == 3)
  { FirstChainNumber=atoi(argv_[2]); }
#endif /* defined MCMC */

  //time(&start_time_);

#ifdef COMPUTE_SPECPHOT_PROPERTIES
  //for dust_model
  mu_seed = -150;
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

  init();

#ifdef STAR_FORMATION_HISTORY
  if(ThisTask == 0)
  { write_sfh_bins(); }
#endif /* defined STAR_FORMATION_HISTORY */

#ifdef MCMC

  /* In MCMC mode only one file is loaded into memory
   * and the sampling for all the steps is done on it */
  sprintf(SimulationDir, "%s/", SimulationDir);
  // time(&start_time_);
  load_tree_table(MCMCTreeSampleFile);
  Senna(); // run the model in MCMC MODE

#else  /* not defined MCMC */

  const int n_tree_files_ = get_nr_files_to_process();
  int *tree_file_number_for_file_ = mymalloc("tree_file_number_for_file_", sizeof(int) * n_tree_files_);
  int *task_for_file_ = mymalloc("task_for_file_", sizeof(int) * n_tree_files_);
  assign_files_to_tasks(tree_file_number_for_file_, task_for_file_, n_tree_files_);

  int file_index_;
  for(file_index_ = 0; file_index_ < n_tree_files_; file_index_++)
  {
    if(ThisTask == task_for_file_[file_index_])
    {
      const int tree_file_number_ = tree_file_number_for_file_[file_index_];

      time(&start_time_);

#ifdef PARALLEL
      do
      { time(&current_time_); }
      while(difftime(current_time_, start_time_) < 1.0 * ThisTask);
#endif /* defined PARALLEL */

      // printf("loading tree table %d.\n", tree_file_number_);

      load_tree_table(tree_file_number_);
      
     // printf("running SAM for tree %d.\n", tree_file_number_);
     
      SAM(tree_file_number_); // run the model in NORMAL MODE

      time(&current_time_);
      printf("\ndone tree file %d in %ldmin and %lds\n\n", tree_file_number_, (current_time_ - start_time_) / 60, (current_time_ - start_time_) % 60);

      free_tree_table();
      //if temporary directory given as argument
      if(argc_ == 3)
      {
#ifdef GALAXYTREE
        sprintf(buffer_, "mv %s/%s_galtree_%d %s", OutputDir,FileNameGalaxies, tree_file_number_, FinalOutputDir);
#else /* not defined GALAXYTREE */
        sprintf(buffer_, "mv %s/%s_z*_%d %s", OutputDir,FileNameGalaxies, tree_file_number_, FinalOutputDir);
#endif /* not defined GALAXYTREE */
        system(buffer_);
      }
    }
  }
  myfree(task_for_file_);
  myfree(tree_file_number_for_file_);
  
#endif /* not defined MCMC */

  gsl_rng_free(random_generator);

#ifdef PARALLEL
  MPI_Finalize();
#endif /* defined PARALLEL */
  return 0;
}
