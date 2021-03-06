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

/** @file   sam.c
 *  @date   2016-2019
 *  @author ?
 *  @author Bruno Henriques
 *  @author Stefan Hilbert
 *
 * @brief SAM Construct Galaxies, Join Galaxies of progenitors, Evolve Galaxies
 *
 * SAM is to be called for each of the chosen dark matter tree files and
 * for each tree, on each treefile: reads tree, constructs the galaxies, saves
 * the galaxies, frees memory for galaxies and tree.
 **/
 
/* library headers: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#ifdef PARALLEL
#include <mpi.h>
#endif /* defined PARALLEL */
 
 
/* other program headers: */
#include "allvars.h"
#include "proto.h"

#ifdef MCMC
#include "mcmc_vars.h"
#include "mcmc_proto.h"
#endif /* defined MCMC */


/** @brief SAM() loops on trees and calls construct_galaxies. 
  *
  *
  * @bug possible bug: if PRELOAD_TREES but not MCMC,
  *      scale_cosmology() is called every time
 */
double SAM(const int tree_file_number_)
{
#ifdef MCMC
  MCMC_GAL = mymalloc("MCMC_Gal", sizeof(struct MCMC_GALAXY[MCMC_MAX_NGALS_PER_OUTPUT]) * NOUT);
  {
    int output_number_;
    for(output_number_ = 0; output_number_<NOUT; output_number_++)
    { TotMCMCGals[output_number_] = 0; }
  }

#ifdef MR_PLUS_MRII
  change_dark_matter_sim("MR");
#else /* not defined MR_PLUS_MRII */
  if(CurrentMCMCStep == 0)
  { read_sample_info(); }
#ifdef HALOMODEL
  else
  {
    int output_number_, fof_number_;
    for(output_number_=0;output_number_<NOUT;output_number_++)
      for(fof_number_=0;fof_number_<NFofsInSample[output_number_];fof_number_++)
        MCMC_FOF[output_number_][fof_number_].NGalsInFoF=0;
  }
#endif /* defined HALOMODEL */
#endif /* not defined MR_PLUS_MRII */
#endif /* defined MCMC */

  //to be used when we have tables for the scaling in any cosmology
  //read_scaling_parameters();

#ifndef MCMC

#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE
  create_galaxy_tree_file(tree_file_number_);
#else /* not defined GALAXYTREE */
  create_galaxy_files(tree_file_number_);
#endif /* not defined GALAXYTREE */
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */
#ifdef LIGHTCONE_OUTPUT
  create_lightcone_galaxy_files(tree_file_number_);
#endif /* defined LIGHTCONE_OUTPUT */ 

#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE
  FILE *tree_n_gal_file_ = fopen("treengal.dat", "w");
#endif /* defined GALAXYTREE */
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */
 
#endif /* not defined MCMC */

//***************************************************************************************
//***************************************************************************************
  TotGalCount = 0;

  //for(tree_number_ = 0; tree_number_ < NTrees_Switch_MR_MRII; tree_number_++)
  int tree_number_;  
  for(tree_number_ = 0; tree_number_ < Ntrees; tree_number_++)
  {
  //printf("doing tree %d of %d\n", tree_number_, Ntrees);
#ifdef MR_PLUS_MRII
    if(tree_number_ == NTrees_Switch_MR_MRII)
    { change_dark_matter_sim("MRII"); }
#endif /* defined MR_PLUS_MRII */

    load_tree(tree_number_);
   
#ifdef MCMC
#ifdef PRELOAD_TREES
    if(CurrentMCMCStep == 0)
#endif /* defined PRELOAD_TREES */
#endif /* defined MCMC */
    { scale_cosmology(TreeNHalos[tree_number_]); }

    gsl_rng_set(random_generator, tree_file_number_ * 100000 + tree_number_);
    NumMergers = 0;
    NHaloGal = 0;
#ifdef GALAXYTREE
    NGalTree = 0;
    IndexStored = 0;
#endif/* defined GALAXYTREE */

#ifdef MCMC
    link_halos_and_MCMC_FOF(tree_number_);
#endif /* defined MCMC */

    /* we process the snapshots now in temporal order 
     * (as a means to reduce peak memory usage).
     * LastSnapShotNr is the highest output snapshot. */
    int snapshot_number_;
    for(snapshot_number_ = 0; snapshot_number_ <= LastSnapShotNr; snapshot_number_++)
    {
#ifdef MCMC
      /* read the appropriate parameter list for current snapshot_number
       * into the parameter variables to be used in construct_galaxies */
      read_mcmc_par(snapshot_number_);
#endif /* defined MCMC */

      int halo_number_;
      for(halo_number_ = 0; halo_number_ < TreeNHalos[tree_number_]; halo_number_++)
      {
#ifdef FIRST_HALO_FIRST
/* construct halos starting with "official" first halo in fof group,
 * which then also takes care of all other halos in that fof group,
 * assuming, all progenitors have been constructed already
 * (i.e. acknowledge loop over snapshots) */
        if(halo_number_ == Halo[halo_number_].FirstHaloInFOFgroup && Halo[halo_number_].SnapNum == snapshot_number_)
          construct_galaxies_in_fof(tree_number_, halo_number_);
#else  /* not defined FIRST_HALO_FIRST */
/* construct halos starting with whatever halos come first in input file,
 * which then also takes care of all other halos in that fof group,
 * not assuming, all progenitors have been constructed already
 * (though they should have because of loop over snapshots) */
         if(HaloAux[halo_number_].DoneFlag == 0 && Halo[halo_number_].SnapNum == snapshot_number_)
           construct_galaxies(tree_number_, halo_number_);
#endif /* not defined FIRST_HALO_FIRST */
      }
    }

    /* output remaining galaxies as needed */
    while(NHaloGal)
    {
#ifdef MCMC 
      if(HaloAux[HaloGal[HaloGalHeap[0]].HaloNr].halo_is_in_MCMC_sample_for_any_output)
#endif /* defined MCMC */  
      output_galaxy(tree_number_, 0);
      pop_galaxy_from_heap(0);
    }

#ifndef MCMC

#ifdef GALAXYTREE
    update_galaxy_tree_ids();
#ifndef LIGHTCONE_OUTPUT_ONLY
    save_galaxy_tree_finalize(tree_file_number_, tree_number_);
#endif /*not defined LIGHTCONE_OUTPUT_ONLY */
#ifdef LIGHTCONE_OUTPUT
    save_lightcone_galaxy_finalize(tree_file_number_, tree_number_);
#endif /* defined LIGHTCONE_OUTPUT */  
#endif /* defined GALAXYTREE */

#if defined OUTPUT_BUFFERING && OUTPUT_BUFFERING == 1
#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE
    save_galaxy_tree_flush_output_buffer();
#else  /* not defined GALAXYTREE */
    save_galaxy_flush_output_buffer();
#endif /* not defined GALAXYTREE */
#endif /*not defined LIGHTCONE_OUTPUT_ONLY */
#ifdef LIGHTCONE_OUTPUT
    save_lightcone_galaxy_flush_output_buffer();
#endif /* defined LIGHTCONE_OUTPUT */  
#endif /* defined OUTPUT_BUFFERING && OUTPUT_BUFFERING == 1 */

#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE /* defined GALAXYTREE */
    fprintf(tree_n_gal_file_, "%d\n", NGalTree);
#endif /* defined GALAXYTREE */
#endif /* defined LIGHTCONE_OUTPUT */ 

#endif /* not defined MCMC */

#ifdef GALAXYTREE /* defined GALAXYTREE */
    TotGalCount += NGalTree;
#else  /* not defined GALAXYTREE */
    {
      int output_number_;
      TotGalCount = 0;
      for(output_number_ = 0;  output_number_ < NOUT; ++output_number_)
      { TotGalCount += TotGalaxies[output_number_]; }
#ifdef LIGHTCONE_OUTPUT 
      TotLightconeGalCount = 0;
      for(output_number_ = 0;  output_number_ < NOUT; ++output_number_)
      { TotLightconeGalCount += TotLightconeGalaxies[output_number_]; }
#endif /* defined LIGHTCONE_OUTPUT */  
    }
#endif /* not defined GALAXYTREE */

#ifndef MCMC
#ifndef PARALLEL
#ifdef LIGHTCONE_OUTPUT
    if(0 == non_roundness(tree_number_)) printf("tree_number_=%d  TotGalCount=%d   TotLightconeGalCount=%lld\n", tree_number_, TotGalCount, TotLightconeGalCount);
#else /* not defined LIGHTCONE_OUTPUT */
    if(0 == non_roundness(tree_number_)) printf("tree_number_=%d  TotGalCount=%d\n", tree_number_, TotGalCount);
#endif /* not defined LIGHTCONE_OUTPUT */
#endif /* not defined PARALLEL */
#endif /* not defined MCMC */

    fflush(stdout);

    free_galaxies_and_tree();
  }//loop on trees

#ifdef MCMC
  const double likelihood_ = get_likelihood();

#ifdef MR_PLUS_MRII
  free_MCMC_FOF();
#else /* not defined MR_PLUS_MRII */
  if(CurrentMCMCStep==ChainLength)
  { free_MCMC_FOF(); }
#endif /* not defined MR_PLUS_MRII */

  myfree(MCMC_GAL);
  return likelihood_;

#else /* not defined MCMC */

#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE
  close_galaxy_tree_file();
#else /* not defined GALAXYTREE */
  close_galaxy_files();
#endif /* not defined GALAXYTREE */
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */
#ifdef LIGHTCONE_OUTPUT  
  close_lightcone_galaxy_files();
#endif /* defined LIGHTCONE_OUTPUT */

#if defined OUTPUT_BUFFERING
#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE
    save_galaxy_tree_show_output_buffer_statistics();
#else  /* not defined GALAXYTREE */
    save_galaxy_show_output_buffer_statistics();
#endif /* not defined GALAXYTREE */
#endif /*not defined LIGHTCONE_OUTPUT_ONLY */
#ifdef LIGHTCONE_OUTPUT
    save_lightcone_galaxy_show_output_buffer_statistics();
#endif /* defined LIGHTCONE_OUTPUT */  
#endif /* defined OUTPUT_BUFFERING */

#ifdef LIGHTCONE_OUTPUT  
  show_lightcone_statistics();
#endif /* defined LIGHTCONE_OUTPUT */

  return 0.;
#endif /* not defined MCMC */
}


/** @brief  construct_galaxies() recursively runs the semi-analytic model.
  *         For each halo it checks if its main progenitor has been done, then
  *         if all the halos in the FOF of its main progenitor have been
  *         done and then runs the SAM in the current halo. This means that
  *         for the first time its called it will walk up the tree into the
  *         haloes in the earliest snapshot.
  *
  *         When it finds a halo that needs to be done it calls
  *         join_galaxies_of_progenitors and evolve_galaxies. */
void construct_galaxies(const int tree_number_, const int halo_number_)
{
  int progenitor_halo_number_, same_fof_halo_number_, n_galaxies_in_fof_group_, fof_merger_center_;
  
#ifdef LIGHTCONE_OUTPUT_ONLY
#ifdef LIGHTCONE_MAY_SKIP_CONSTRUCT_GALAXY
  /* if halos is outside lightcone, don't construct, but return early */
  if(is_outside_lightcone_for_snapshot[Halo[halo_number_].SnapNum])
  {
    lightcone_N_galaxies_skipped_construction++;
    return;  
  }
  
  if(check_outside_lightcone_for_snapshot[Halo[halo_number_].SnapNum])
  {
    float d_0_ = Halo[halo_number_].Pos[0] - lightcone_observer_position[0]; d_0_ = wrap(d_0_, BoxSize);
    float d_1_ = Halo[halo_number_].Pos[1] - lightcone_observer_position[1]; d_1_ = wrap(d_1_, BoxSize);
    float d_2_ = Halo[halo_number_].Pos[2] - lightcone_observer_position[2]; d_2_ = wrap(d_2_, BoxSize);

    if(d_0_ * d_0_ + d_1_ * d_1_ + d_2_ * d_2_ > pow2(lightcone_radius_for_snapshot[Halo[halo_number_].SnapNum]))
    {
      lightcone_N_galaxies_skipped_construction++;
      return;
    }
  }
#endif /* defined LIGHTCONE_MAY_SKIP_CONSTRUCT_GALAXY */
#endif /* defined LIGHTCONE_OUTPUT_ONLY */

  HaloAux[halo_number_].DoneFlag = 1;

  /* construct all galaxies in all progenitors of this halo (if not already done): */
  for(progenitor_halo_number_ = Halo[halo_number_].FirstProgenitor; progenitor_halo_number_ >= 0; progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor)
  {
    if(HaloAux[progenitor_halo_number_].DoneFlag == 0) /* progenitor hasn't been done yet */
    { construct_galaxies(tree_number_, progenitor_halo_number_); }/* so construct progenitor's galaxies now */
  }

  /* Construct all galaxies in all progenitors of all the halos in the
   * current FOF group, unless laready done.
   * First call to construct_galaxies for a halo of a fof group will construct
   * all progenitors of all halos in the group. So, later calls to
   * construct_galaxies for other halos of the same fof group will find that
   * work already done. 
   */
  const int first_in_fof_halo_number_ = Halo[halo_number_].FirstHaloInFOFgroup;
  
  if(HaloAux[first_in_fof_halo_number_].HaloFlag == 0) /* halo hasn't been done */
  {
    HaloAux[first_in_fof_halo_number_].HaloFlag = 1;  /* mark as to do */
    for(same_fof_halo_number_ = first_in_fof_halo_number_ ; same_fof_halo_number_ >= 0; same_fof_halo_number_ = Halo[same_fof_halo_number_].NextHaloInFOFgroup)  /* fof still has a halo to process */
    {
      for(progenitor_halo_number_ = Halo[same_fof_halo_number_].FirstProgenitor; progenitor_halo_number_ >= 0;  progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor) //build its progenitors
      {
        if(HaloAux[progenitor_halo_number_].DoneFlag == 0)
        { construct_galaxies(tree_number_, progenitor_halo_number_); }
      }
    }

    /* At this point, the galaxies for all progenitors of this halo have been
     * properly constructed. Also, the galaxies of the progenitors of all other 
     * halos in the same FOF-group have been constructed as well. We can hence go
     * ahead and construct all galaxies for the subhalos in this FOF halo, and
     * evolve them in time. */
    HaloAux[first_in_fof_halo_number_].HaloFlag = 2;

    /*For all the halos in the current FOF join all the progenitor galaxies together
      * n_galaxies_in_fof_group_ will be the total number of galaxies in the current FOF*/
    n_galaxies_in_fof_group_ = 0;
    fof_merger_center_ = get_merger_center(first_in_fof_halo_number_);     //Find type 0 for type 1 to merge into
    for(same_fof_halo_number_ = first_in_fof_halo_number_ ; same_fof_halo_number_ >= 0; same_fof_halo_number_ = Halo[same_fof_halo_number_].NextHaloInFOFgroup) 
    { join_galaxies_of_progenitors(same_fof_halo_number_, &n_galaxies_in_fof_group_, &fof_merger_center_); }

    /*Evolve the Galaxies -> SAM! */
    evolve_galaxies(first_in_fof_halo_number_, n_galaxies_in_fof_group_, tree_number_, fof_merger_center_);

#ifdef MASS_CHECKS
    {
      int galaxy_number_;
      for (galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
      { mass_checks("Construct_galaxies #1", galaxy_number_); }
    }
#endif /* defined MASS_CHECKS */
  }
}


/** @brief  construct_galaxies_in_fof() runs the semi-analytic model.
  *         For each FOF group it assumes that all the halos of all its
  *         progenitors have been done and then runs the SAM on all halos
  *         of the current FOF (i.e. runs join_galaxies_of_progenitors
  *         and evolve_galaxies). */
void construct_galaxies_in_fof(const int tree_number_, const int first_in_fof_halo_number_)
{
  int same_fof_halo_number_, n_galaxies_in_fof_group_, fof_merger_center_;

#ifdef LIGHTCONE_OUTPUT_ONLY
#ifdef LIGHTCONE_MAY_SKIP_CONSTRUCT_GALAXY
  /* if halos is outside lightcone, don't construct, but return early */
  if(is_outside_lightcone_for_snapshot[Halo[first_in_fof_halo_number_].SnapNum])
  {
    lightcone_N_fof_groups_skipped_construction++;
    return;  
  }
  
  if(check_outside_lightcone_for_snapshot[Halo[first_in_fof_halo_number_].SnapNum])
  {
    bool all_outside_ = true;
    for(same_fof_halo_number_ = first_in_fof_halo_number_; same_fof_halo_number_ >= 0 && all_outside_ == true; same_fof_halo_number_ = Halo[same_fof_halo_number_].NextHaloInFOFgroup)
    {
      float d_0_ = Halo[same_fof_halo_number_].Pos[0] - lightcone_observer_position[0]; d_0_ = wrap(d_0_, BoxSize);
      float d_1_ = Halo[same_fof_halo_number_].Pos[1] - lightcone_observer_position[1]; d_1_ = wrap(d_1_, BoxSize);
      float d_2_ = Halo[same_fof_halo_number_].Pos[2] - lightcone_observer_position[2]; d_2_ = wrap(d_2_, BoxSize);

      all_outside_ &= (d_0_ * d_0_ + d_1_ * d_1_ + d_2_ * d_2_ > pow2(lightcone_radius_for_snapshot[Halo[same_fof_halo_number_].SnapNum]));
    }
    if(all_outside_)
    {
      lightcone_N_fof_groups_skipped_construction++;
      return;      
    }
  }
#endif /* defined LIGHTCONE_MAY_SKIP_CONSTRUCT_GALAXY */
#endif /* defined LIGHTCONE_OUTPUT_ONLY */

  for(same_fof_halo_number_ = first_in_fof_halo_number_; same_fof_halo_number_ >= 0; same_fof_halo_number_ = Halo[same_fof_halo_number_].NextHaloInFOFgroup)
  { HaloAux[same_fof_halo_number_].DoneFlag = 1; }

  /* At this point, the galaxies for all progenitors of this halo have been
   * properly constructed. Also, the galaxies of the progenitors of all other 
   * halos in the same FOF-group have been constructed as well. We can hence go
   * ahead and construct all galaxies for the subhalos in this FOF halo, and
   * evolve them in time. */

  HaloAux[first_in_fof_halo_number_].HaloFlag = 2;

  /*For all the halos in the current FOF join all the progenitor galaxies together
    * n_galaxies_in_fof_group_ will be the total number of galaxies in the current FOF*/
  n_galaxies_in_fof_group_ = 0;
  fof_merger_center_ = get_merger_center(first_in_fof_halo_number_);     //Find type 0 for type 1 to merge into
  for(same_fof_halo_number_ = first_in_fof_halo_number_; same_fof_halo_number_ >= 0; same_fof_halo_number_ = Halo[same_fof_halo_number_].NextHaloInFOFgroup)
  { join_galaxies_of_progenitors(same_fof_halo_number_, &n_galaxies_in_fof_group_, &fof_merger_center_); }

  /*Evolve the Galaxies -> SAM! */
  evolve_galaxies(first_in_fof_halo_number_, n_galaxies_in_fof_group_, tree_number_, fof_merger_center_);

#ifdef MASS_CHECKS
  {
    int galaxy_number_;
    for (galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
    { mass_checks("Construct_galaxies #1", galaxy_number_); }
  }
#endif /* defined MASS_CHECKS */
}


/** @brief join_galaxies_of_progenitors() updates the properties of the
 *         galaxy from the dark matter halo properties and deals with
 *         merging clocks. This routine is called by construct_galaxies
 *         for every halo in the FOF being constructed. When there is no
 *         galaxy in the Halo of FirstProgenitor, the first_occupied_halo_number_
 *         pointer is changed to a subhalo which have the maximum mass.
 *    
 *         For a central galaxy it just updates its properties. For
 *         satellites it needs to know its most massive (or only progenitor)
 *         to keep track of the merging clock. It also finds the central
 *         galaxies into which galaxies should merge. Type 1's
 *         can merge if their baryonic mass is bigger than the dark matter
 *         mass and type 2's can merge into them. Once the type 1's merge
 *         into a type 0 all its satellites will have the merging clock
 *         into the type 0 reset .
 * */
void join_galaxies_of_progenitors(const int halo_number_, int *n_galaxies_in_fof_group_, int *merger_center_)
{
  int i_, j_, progenitor_halo_number_, first_occupied_halo_number_, central_galaxy_number_, most_massive_halo_number_, most_massive_halo_length_, galaxy_number_;
  
  const int galaxy_number_beg_ = *n_galaxies_in_fof_group_;
        int galaxy_number_end_ = *n_galaxies_in_fof_group_;

  /* If there is no galaxy in the Halo of FirstProgenitor,
   * but there are other progenitors with galaxies of type 0 or 1,
   * the first_occupied_halo_number_ pointer is changed 
   * to a subhalo which has the maximum mass.
   * This should only happen in the case that
   * the leaf on the firstprogenitor branch occurs
   * as a subhalo, in that case no galaxy would be assigned to it. */
  progenitor_halo_number_ = Halo[halo_number_].FirstProgenitor;
  first_occupied_halo_number_ = Halo[halo_number_].FirstProgenitor;
  if(progenitor_halo_number_ >= 0 && HaloAux[progenitor_halo_number_].NGalaxies == 0) //If halo has progenitors, but progenitor has no galaxies
  {
    most_massive_halo_length_ = 0;
    for(;progenitor_halo_number_ >= 0; progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor)
    {
      for(i_ = 0, galaxy_number_ = HaloAux[progenitor_halo_number_].FirstGalaxy; i_ < HaloAux[progenitor_halo_number_].NGalaxies; ++i_, galaxy_number_ = HaloGal[galaxy_number_].NextGalaxy)
        if((HaloGal[galaxy_number_].Type == 0 || HaloGal[galaxy_number_].Type == 1) && (Halo[progenitor_halo_number_].Len > most_massive_halo_length_))
        {
          most_massive_halo_length_ = Halo[progenitor_halo_number_].Len;
          first_occupied_halo_number_ = progenitor_halo_number_;  //define the new first_occupied_halo_number_
          break;
        }
    }
  }

  /* loop through all the progenitors and get the halo mass and ID of the most massive */
  most_massive_halo_length_ = 0;
  most_massive_halo_number_ = Halo[halo_number_].FirstProgenitor;
  for(progenitor_halo_number_ = Halo[halo_number_].FirstProgenitor; progenitor_halo_number_ >= 0; progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor)
  {
    if(Halo[progenitor_halo_number_].Len > most_massive_halo_length_)
    {
      most_massive_halo_length_ = Halo[progenitor_halo_number_].Len;
      most_massive_halo_number_ = progenitor_halo_number_;
    }
  }

  for(progenitor_halo_number_ = Halo[halo_number_].FirstProgenitor; progenitor_halo_number_ >= 0;  progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor)
  {
    for(i_ = 0, galaxy_number_ = HaloAux[progenitor_halo_number_].FirstGalaxy; i_ < HaloAux[progenitor_halo_number_].NGalaxies; ++i_, galaxy_number_ = HaloGal[galaxy_number_].NextGalaxy, ++galaxy_number_end_)
    {
      if(galaxy_number_end_ >= MaxGal)
      {
        AllocValue_MaxGal *= ALLOC_INCREASE_FACTOR;
        MaxGal = AllocValue_MaxGal;
        if(MaxGal < galaxy_number_end_ + 1) { MaxGal = galaxy_number_end_ + 1; }
        Gal = myrealloc_movable(Gal, sizeof(struct GALAXY) * MaxGal);
      }
      if(*merger_center_ == galaxy_number_)
      {  *merger_center_ =  galaxy_number_end_; }

      /* Copy galaxy properties from progenitor,
        * except for those that need initialising */
      Gal[galaxy_number_end_]                         = HaloGal[galaxy_number_];
      Gal[galaxy_number_end_].HaloNr                  = halo_number_;
      Gal[galaxy_number_end_].CoolingRadius           = 0.0;
      Gal[galaxy_number_end_].CoolingGas              = 0.0;
      Gal[galaxy_number_end_].PrimordialAccretionRate = 0.0;
      Gal[galaxy_number_end_].CoolingRate             = 0.0;
      Gal[galaxy_number_end_].CoolingRate_beforeAGN   = 0.0;
      Gal[galaxy_number_end_].Sfr                     = 0.0;
      Gal[galaxy_number_end_].SfrBulge                = 0.0;
      Gal[galaxy_number_end_].QuasarAccretionRate     = 0.0;
      Gal[galaxy_number_end_].RadioAccretionRate      = 0.0;
#ifdef GALAXYTREE
      Gal[galaxy_number_end_].FirstProgGal            = HaloGal[galaxy_number_].GalTreeIndex;    /* CHECK */
#endif /* defined GALAXYTREE */
 
#ifdef MASS_CHECKS
      // To fail this check means that we copy in a failed galaxy
      mass_checks("Middle of join_galaxies_of_progenitors", galaxy_number_end_);
#endif /* defined MASS_CHECKS */

      /* Update Properties of this galaxy with physical properties of halo */
      /* this deals with the central galaxies of subhalos */
      if(Gal[galaxy_number_end_].Type == 0 || Gal[galaxy_number_end_].Type == 1)
      {
        if(progenitor_halo_number_ == first_occupied_halo_number_)
        {
#ifdef HALOPROPERTIES
          Gal[galaxy_number_end_].HaloM_Mean200 = Halo[halo_number_].M_Mean200;
          Gal[galaxy_number_end_].HaloM_Crit200 = Halo[halo_number_].M_Crit200;
          Gal[galaxy_number_end_].HaloM_TopHat  = Halo[halo_number_].M_TopHat;
          Gal[galaxy_number_end_].HaloVelDisp   = Halo[halo_number_].VelDisp;
          Gal[galaxy_number_end_].HaloVmax      = Halo[halo_number_].Vmax;
#endif /* defined HALOPROPERTIES */
          Gal[galaxy_number_end_].MostBoundID = Halo[halo_number_].MostBoundID;
          for(j_ = 0; j_ < 3; j_++)
          {
            Gal[galaxy_number_end_].Pos[j_]     = Halo[halo_number_].Pos[j_];
            Gal[galaxy_number_end_].Vel[j_]     = Halo[halo_number_].Vel[j_];
#ifdef HALOPROPERTIES
            Gal[galaxy_number_end_].HaloPos[j_] = Halo[halo_number_].Pos[j_];
            Gal[galaxy_number_end_].HaloVel[j_] = Halo[halo_number_].Vel[j_];
#endif /* defined HALOPROPERTIES */
          }
          Gal[galaxy_number_end_].Len           = Halo[halo_number_].Len;
          Gal[galaxy_number_end_].Vmax          = Halo[halo_number_].Vmax;

          // FOFCentralGal property in case that is different from FirstGalaxy
          if(halo_number_ == Halo[halo_number_].FirstHaloInFOFgroup)
          { update_centralgal(galaxy_number_end_, halo_number_); }
          else
          { update_type_1(galaxy_number_end_, halo_number_, progenitor_halo_number_); }

          if(DiskRadiusModel == 1 || DiskRadiusModel == 2)
          {
            Gal[galaxy_number_end_].GasDiskRadius     = get_disk_radius(halo_number_, galaxy_number_end_);
            Gal[galaxy_number_end_].StellarDiskRadius = Gal[galaxy_number_end_].GasDiskRadius;
          }
        }
        else //type 2 galaxies
        { update_type_2(galaxy_number_end_, halo_number_, progenitor_halo_number_, most_massive_halo_number_); }
      }

      /* Note: Galaxies that are already type=2 do not need a special treatment at this point */
      if(Gal[galaxy_number_end_].Type < 0 || Gal[galaxy_number_end_].Type > 2)
      { terminate("Unknown galaxy type\n"); }
    }
  }

  /* If there are no progenitors with galaxies, a new galaxy is created.
   * However, if it's a subhalo, no galaxy is placed, since it would stay
   * at zero luminosity. */
  if(galaxy_number_end_ == 0)
  {
    *merger_center_ = 0;
    if(Halo[halo_number_].FirstHaloInFOFgroup == halo_number_)
    {
      init_galaxy(galaxy_number_end_, halo_number_);
      galaxy_number_end_++;
    }
  }

  /* satelites (type 2's) will preferably merge onto this type 1 rather than the type 0 */
  central_galaxy_number_ = -1;
  for(i_ = galaxy_number_beg_; i_ < galaxy_number_end_; i_++)
    if(Gal[i_].Type == 0 || Gal[i_].Type == 1)
    {
      if(central_galaxy_number_ != -1)
      { terminate("Subhalo hosts more than one Type 0/1\n"); }

      central_galaxy_number_ = i_;
    }

  for(i_ = galaxy_number_beg_; i_ < galaxy_number_end_; i_++)
  {
    Gal[i_].CentralGal = central_galaxy_number_;
    if(central_galaxy_number_ != -1)
      for(j_ = 0; j_ < 3; j_++)
      { Gal[i_].MergCentralPos[j_] = Gal[central_galaxy_number_].Pos[j_]; }
  }

  /* Satellites whose type 1 has merged into type 0, will be reset to merge
   * into the type 0. */
  if(central_galaxy_number_ == -1 && galaxy_number_end_ != galaxy_number_beg_)
  {
    for(i_ = galaxy_number_beg_; i_ < galaxy_number_end_; i_++)
    {
      Gal[i_].CentralGal = *merger_center_;
      for(j_ = 0; j_ < 3; j_++)
      { Gal[i_].MergCentralPos[j_] = Gal[*merger_center_].Pos[j_]; }
    }
  }
    
#ifdef MASS_CHECKS
  for (i_ = galaxy_number_beg_; i_<galaxy_number_end_; i_++)
  { mass_checks("Bottom of join_galaxies_of_progenitors",i_); }
#endif /* defined MASS_CHECKS */

  report_memory_usage(&HighMark, "join_galaxies");

  *n_galaxies_in_fof_group_ = galaxy_number_end_;
}


/** @brief evolve_galaxies() deals with most of the SA recipes.
  *
  *        Here most of the physical recipes are called, including:
  *        infall_recipe() (gets the fraction of primordial gas that infalled),
  *        add_infall_to_hot() (adds the infalled gas to the hot phase),
  *        reincorporate_gas() (reincorporates gas ejected by SN),
  *        cooling_recipe() (gets the amount of gas that cooled - takes into
  *        account AGN feedback), cool_gas_onto_galaxy() (adds the gas that
  *        cooled into the cold phase, starformation_and_feedback() (normal
  *        SF and SN feedback), deal_with_galaxy_merger() (adds components,
  *        grows black hole and deals with SF burst), disruption() (total and
  *        instantaneous disruption of type 2 satellites) and dust() (if galaxy
  *        is in an output time, dust extinction is computed).
  * 
  *        All these calculations are done in time steps of 1/STEPS the time
  *        between each snapshot (STEPS=20).
  * 
  * @note  halo_number_ is here the FOF-background subhalo (i.e. main halo) 
  */
void evolve_galaxies(const int halo_number_, const int n_galaxies_in_fof_group_, const int tree_number_, const int galaxy_number_for_merger_center_of_type_1_galaxies_)
{
  int galaxy_number_, other_galaxy_number_, step_number_, merger_centralgal_, current_halo_number_, progenitor_halo_number_, fof_central_galaxy_number_, previous_galaxy_number_, i_;
  double current_time_;

#ifdef STAR_FORMATION_HISTORY
  double age_in_years_;
#endif /* defined STAR_FORMATION_HISTORY */

  if(n_galaxies_in_fof_group_ <= 0)
    return;

  //previous_time_         = NumToTime(Gal[0].SnapNum);
  const double previous_time_ = NumToTime(Halo[halo_number_].SnapNum-1);
  const double new_time_      = NumToTime(Halo[halo_number_].SnapNum  );
  const double deltaT_        = previous_time_ - new_time_; /* Time between snapshots */
  const double dt_            = (previous_time_ - new_time_) / STEPS; 

#ifdef MASS_CHECKS
  for (galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
  { mass_checks("Evolve_galaxies #0", galaxy_number_); }
#endif /* defined MASS_CHECKS */

  fof_central_galaxy_number_ = Gal[0].CentralGal;
  //print_galaxy("\n\ncheck1", fof_central_galaxy_number_, halo_number_);
  if(Gal[fof_central_galaxy_number_].Type != 0 || Gal[fof_central_galaxy_number_].HaloNr != halo_number_)
  { terminate("Something wrong here ..... \n"); }

  /* Update all galaxies to same star-formation history time-bins.
   * Needed in case some galaxy has skipped a snapshot. */
#ifdef STAR_FORMATION_HISTORY
  age_in_years_ = (Age[0] - previous_time_) * UnitTime_in_years * inv_Hubble_h; //ROB: age_in_years_ is in units of "real years"!
  step_number_  = 0;
  for(galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
  { sfh_update_bins(galaxy_number_, Halo[halo_number_].SnapNum - 1, step_number_, age_in_years_); }
#endif /* defined STAR_FORMATION_HISTORY */

  /* Handle the transfer of mass between satellites and central galaxies */
  deal_with_satellites(fof_central_galaxy_number_, n_galaxies_in_fof_group_);

  /* Delete inconsequential galaxies */
  for (galaxy_number_ =0; galaxy_number_< n_galaxies_in_fof_group_; galaxy_number_++)
  {
    if (Gal[galaxy_number_].Type == 2 && Gal[galaxy_number_].ColdGas+Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass < 1.e-8)
    { Gal[galaxy_number_].Type = 3; }
#ifdef MASS_CHECKS
    else
    { mass_checks("Evolve_galaxies #0.1",galaxy_number_); }
#endif /* defined MASS_CHECKS */
  }
   
  /* Calculate how much hot gas needs to be accreted to give the correct baryon fraction
   * in the main halo. This is the universal fraction, less any reduction due to reionization. */
   
  const double infalling_gas_per_step_ = infall_recipe(fof_central_galaxy_number_, n_galaxies_in_fof_group_, ZZ[Halo[halo_number_].SnapNum]) / STEPS;
  Gal[fof_central_galaxy_number_].PrimordialAccretionRate = infalling_gas_per_step_ / dt_;
  
  /* All the physics are computed in a number of intervals between snapshots
   * equal to STEPS */
  for (step_number_ = 0; step_number_ < STEPS; step_number_++)
  {
    /* time to present of the current step */
    current_time_ = previous_time_ - (step_number_ + 0.5) * dt_;

    /* Update all galaxies to the star-formation history time-bins of current step*/
#ifdef STAR_FORMATION_HISTORY
    age_in_years_=(Age[0]-current_time_)*UnitTime_in_years * inv_Hubble_h;
    for (galaxy_number_=0; galaxy_number_<n_galaxies_in_fof_group_; galaxy_number_++)
    { sfh_update_bins(galaxy_number_,Halo[halo_number_].SnapNum - 1, step_number_, age_in_years_); }

#endif /* defined STAR_FORMATION_HISTORY */

    /* Infall onto central galaxy only, if required to make up a baryon deficit */
#ifndef GUO10
#ifndef GUO13
    if (infalling_gas_per_step_ > 0.)
#endif /* not defined GUO13 */
#endif /* not defined GUO10 */
    { add_infall_to_hot(fof_central_galaxy_number_, infalling_gas_per_step_); }

#ifdef MASS_CHECKS
    mass_checks("Evolve_galaxies #0.5",fof_central_galaxy_number_);
#endif /* defined MASS_CHECKS */

    for(galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
    {
      /* don't treat galaxies that have already merged */
      if(Gal[galaxy_number_].Type == 3)
      { continue; }
      
#ifdef MASS_CHECKS
      mass_checks("Evolve_galaxies #1",galaxy_number_);
#endif /* defined MASS_CHECKS */

      if (Gal[galaxy_number_].Type == 0 || Gal[galaxy_number_].Type == 1)
      {
        reincorporate_gas(galaxy_number_, dt_);
        /* determine cooling gas given halo properties and add it to the cold phase*/
#ifdef MASS_CHECKS
        mass_checks("Evolve_galaxies #1.5",galaxy_number_);
#endif /* defined MASS_CHECKS */

        compute_cooling(galaxy_number_, dt_);
      }
    }

    //this must be separated as now satellite AGN can heat central galaxies
    //therefore the AGN from all satellites must be computed, in a loop inside this function,
    //before gas is cooled into central galaxies (only suppress cooling, the gas is not actually heated)
    if(AGNRadioModeModel != 5)
    { do_AGN_heating(dt_, n_galaxies_in_fof_group_); }

    for (galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
    {
      cool_gas_onto_galaxy(galaxy_number_, dt_);
      
#ifdef MASS_CHECKS
      mass_checks("Evolve_galaxies #2",galaxy_number_);
#endif /* defined MASS_CHECKS */

      starformation(galaxy_number_, fof_central_galaxy_number_, current_time_, dt_);
      
#ifdef MASS_CHECKS
      mass_checks("Evolve_galaxies #3",galaxy_number_);
#endif /* defined MASS_CHECKS */

      //print_galaxy("check3", fof_central_galaxy_number_, halo_number_);
    }

    /* Check for merger events */
    for(galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
    {
      if(Gal[galaxy_number_].Type == 2 || (Gal[galaxy_number_].Type == 1 && Gal[galaxy_number_].MergeOn == 1))             /* satellite galaxy */
      {
        Gal[galaxy_number_].MergTime -= dt_;
        if(Gal[galaxy_number_].MergTime < 0.0)
        {
          NumMergers++;

          if(Gal[galaxy_number_].Type == 1)
            for(other_galaxy_number_ = 0; other_galaxy_number_ < n_galaxies_in_fof_group_; other_galaxy_number_++)
              if(Gal[other_galaxy_number_].Type == 2 && Gal[galaxy_number_].CentralGal == galaxy_number_)
              { Gal[other_galaxy_number_].CentralGal = galaxy_number_for_merger_center_of_type_1_galaxies_; }

          if(Gal[galaxy_number_].Type == 2)
          { merger_centralgal_ = Gal[galaxy_number_].CentralGal; }
          else
          { merger_centralgal_ = galaxy_number_for_merger_center_of_type_1_galaxies_; }
        
#ifdef MASS_CHECKS
          mass_checks("Evolve_galaxies #4",galaxy_number_);
          mass_checks("Evolve_galaxies #4",merger_centralgal_);
          mass_checks("Evolve_galaxies #4",fof_central_galaxy_number_);
#endif /* defined MASS_CHECKS */
          
           /* note: call to deal_with_galaxy_merger uses deltaT_ as parameter, and internally divides by STEPS where needed */
          deal_with_galaxy_merger(galaxy_number_, merger_centralgal_, fof_central_galaxy_number_, current_time_, deltaT_);
          
#ifdef MASS_CHECKS
          mass_checks("Evolve_galaxies #5",galaxy_number_);
          mass_checks("Evolve_galaxies #5",merger_centralgal_);
          mass_checks("Evolve_galaxies #5",fof_central_galaxy_number_);
#endif /* defined MASS_CHECKS */
        }
      }
    }//loop on all galaxies to detect mergers

#ifdef DETAILED_METALS_AND_MASS_RETURN
    //DELAYED ENRICHMENT AND MASS RETURN + FEEDBACK: No fixed yield or recycling fraction anymore. FB synced with enrichment
    for (galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
    { update_yields_and_return_mass(galaxy_number_, fof_central_galaxy_number_, dt_, step_number_); }
#endif /* defined DETAILED_METALS_AND_MASS_RETURN */

  }/* end move forward in interval STEPS */

  for(galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
  {
    if(Gal[galaxy_number_].Type == 2)
    {
#ifndef UPDATETYPETWO
      int jj;
      float tmppos;
      for(jj = 0; jj < 3; jj++)
      {
        tmppos = wrap(Gal[galaxy_number_].DistanceToCentralGal[jj],BoxSize);
        tmppos *=  2.*sqrt(Gal[galaxy_number_].MergTime/Gal[galaxy_number_].OriMergTime);
        Gal[galaxy_number_].Pos[jj] = Gal[galaxy_number_].MergCentralPos[jj] + tmppos;

        if(Gal[galaxy_number_].Pos[jj] < 0)
                  Gal[galaxy_number_].Pos[jj] = BoxSize + Gal[galaxy_number_].Pos[jj];
        if(Gal[galaxy_number_].Pos[jj] > BoxSize)
                  Gal[galaxy_number_].Pos[jj] = Gal[galaxy_number_].Pos[jj] - BoxSize;
      }
#endif /* not defined UPDATETYPETWO */
      /* Disruption of type 2 galaxies. Type 1 galaxies are not disrupted since usually
        * bayonic component is more compact than dark matter.*/

      if(DisruptionModel==0)
      { disrupt(galaxy_number_); }
    }
  }
  
#ifdef MASS_CHECKS
  for (galaxy_number_ =0;galaxy_number_<n_galaxies_in_fof_group_;galaxy_number_++)
  { mass_checks("Evolve_galaxies #6",galaxy_number_); }
#endif /* defined MASS_CHECKS */

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef  POST_PROCESS_MAGS
  /* If this is an output snapshot apply the dust model to each galaxy */
  int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    if(Halo[halo_number_].SnapNum == ListOutputSnaps[output_number_])
    {
      for(galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
      { dust_model(galaxy_number_, output_number_); }
      break;
    }
  }
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

  /* now save the galaxies of all the progenitors (and free the associated storage) */
  for(progenitor_halo_number_ = Halo[halo_number_].FirstProgenitor; progenitor_halo_number_ >= 0;  progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor)
  {
    for(i_ = 0, galaxy_number_ = HaloAux[progenitor_halo_number_].FirstGalaxy; i_ < HaloAux[progenitor_halo_number_].NGalaxies; i_++)
    {
      int next_galaxy_number_ = HaloGal[galaxy_number_].NextGalaxy;
#ifdef MCMC
      if(HaloAux[progenitor_halo_number_].halo_is_in_MCMC_sample_for_any_output)
#endif /* defined MCMC */     
      output_galaxy(tree_number_, HaloGal[galaxy_number_].HeapIndex);
      pop_galaxy_from_heap(HaloGal[galaxy_number_].HeapIndex);
      galaxy_number_ = next_galaxy_number_;
    }
  }
  
#ifdef GALAXYTREE
  const int evolved_galaxies_begin_ = NGalTree;
#endif /* defined GALAXYTREE */
  for(galaxy_number_ = 0, previous_galaxy_number_ = -1, current_halo_number_ = -1, fof_central_galaxy_number_ = -1; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
  {
    if(Gal[galaxy_number_].HaloNr != current_halo_number_)
    {
      current_halo_number_ = Gal[galaxy_number_].HaloNr;
      HaloAux[current_halo_number_].FirstGalaxy = -1;
      HaloAux[current_halo_number_].NGalaxies = 0;
    }

#ifdef MASS_CHECKS
    mass_checks("Evolve_galaxies #7",galaxy_number_);
#endif /* defined MASS_CHECKS */

    if(Gal[galaxy_number_].Type != 3)
    {
      if(NHaloGal >= MaxHaloGal)
      {
        int oldmax = MaxHaloGal;
        AllocValue_MaxHaloGal *= ALLOC_INCREASE_FACTOR;
        MaxHaloGal = AllocValue_MaxHaloGal;
        if(MaxHaloGal<NHaloGal+1)
        { MaxHaloGal=NHaloGal+1; }
        HaloGal = myrealloc_movable(HaloGal, sizeof(struct GALAXY) * MaxHaloGal);
        HaloGalHeap = myrealloc_movable(HaloGalHeap, sizeof(int) * MaxHaloGal);
        for(i_ = oldmax; i_ < MaxHaloGal; i_++)
        { HaloGalHeap[i_] = i_; }
      }

      Gal[galaxy_number_].SnapNum = Halo[current_halo_number_].SnapNum;

#ifndef GUO10
#ifdef UPDATETYPETWO
      update_type_two_coordinate_and_velocity(tree_number_, galaxy_number_, Gal[0].CentralGal);
#endif /* defined UPDATETYPETWO */
#endif /* not defined GUO10 */

      /* when galaxies are outputed, the slot is filled with the
        * last galaxy in the heap. New galaxies always take the last spot */
      int next_galaxy_number_ = HaloGalHeap[NHaloGal];
      HaloGal[next_galaxy_number_] = Gal[galaxy_number_];
      HaloGal[next_galaxy_number_].HeapIndex = NHaloGal;

      if(HaloAux[current_halo_number_].FirstGalaxy < 0)
      { HaloAux[current_halo_number_].FirstGalaxy = next_galaxy_number_; }

      if(previous_galaxy_number_ >= 0)
      { HaloGal[previous_galaxy_number_].NextGalaxy = next_galaxy_number_; }
      previous_galaxy_number_ = next_galaxy_number_;

      HaloAux[current_halo_number_].NGalaxies++;
      NHaloGal++;

#ifdef GALAXYTREE
      if(NGalTree >= MaxGalTree)
      {
        AllocValue_MaxGalTree *= ALLOC_INCREASE_FACTOR;
        MaxGalTree = AllocValue_MaxGalTree;
        if(MaxGalTree<NGalTree+1) MaxGalTree=NGalTree+1;
        GalTree = myrealloc_movable(GalTree, sizeof(struct galaxy_tree_data) * MaxGalTree);
      }
      HaloGal[next_galaxy_number_].GalTreeIndex = NGalTree;

      memset(&GalTree[NGalTree], 0, sizeof(struct galaxy_tree_data));
      GalTree[NGalTree].HaloGalIndex = next_galaxy_number_;
      GalTree[NGalTree].SnapNum = Halo[current_halo_number_].SnapNum;
      GalTree[NGalTree].NextProgGal = -1;
      GalTree[NGalTree].DescendantGal = -1;

      GalTree[NGalTree].FirstProgGal = Gal[galaxy_number_].FirstProgGal;
      if(Gal[galaxy_number_].Type == 0)
      { fof_central_galaxy_number_ = NGalTree; }
      NGalTree++;
#endif /* defined GALAXYTREE */
    }
  }

#ifdef GALAXYTREE
  if(fof_central_galaxy_number_ < 0)
  { terminate("fof_central_galaxy_number_ < 0, i_.e. did not find central galaxy"); }

  for(galaxy_number_ = evolved_galaxies_begin_; galaxy_number_ < NGalTree; galaxy_number_++)
  { GalTree[galaxy_number_].FOFCentralGal = fof_central_galaxy_number_; }
#endif /* defined GALAXYTREE */

  report_memory_usage(&HighMark, "evolve_galaxies");
}


/** @brief output_galaxy() outputs galaxy to various files on disk
 * 
 *  outputs galaxy to various files on disk,
 *  alos does some book-keeping needed for later updating tree info
 *  on disk,
 *  after writing, pops galaxy from heap (temporary memory)
 */
void output_galaxy(const int tree_number_, const int heap_index_)
{
  (void) tree_number_; /* avoid unused-parameter warning for certain compiler options */
  
#ifndef MCMC
#ifndef GALAXYTREE
  int output_number_;
#endif /* not defined GALAXYTREE */
#endif /* not defined MCMC */
  
  const int galaxy_index_ = HaloGalHeap[heap_index_];

  if(heap_index_ >= NHaloGal)
  { terminate("heap_index_ >= NHaloGal"); }
 
  if(HaloGal[galaxy_index_].HeapIndex != heap_index_)                // consistency check
  { terminate("this should not happen"); }

#ifdef GUO10
#ifdef UPDATETYPETWO
  update_type_two_coordinate_and_velocity(tree_number_, galaxy_index_, HaloGal[0].CentralGal);
#endif /* defined GUO10 */
#endif /* defined UPDATETYPETWO */

#ifdef MCMC
  /* if MCMC galaxies are saved into a structure to then be compared
   * with observations. The output will come in form of a file with
   * parameters, not galaxy properties */

  save_galaxy_for_mcmc(galaxy_index_);

#else  /* not defined MCMC */

#ifndef LIGHTCONE_OUTPUT_ONLY  
#ifdef GALAXYTREE
  GalTree[HaloGal[galaxy_index_].GalTreeIndex].IndexStored = IndexStored++;
  save_galaxy_tree_append(galaxy_index_);
#else /* not defined GALAXYTREE */
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    if(ListOutputSnaps[output_number_] == HaloGal[galaxy_index_].SnapNum)
    { save_galaxy_append(tree_number_, galaxy_index_, output_number_); }
#endif /* not defined GALAXYTREE */
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */

#ifdef LIGHTCONE_OUTPUT
#ifdef GALAXYTREE
  GalTree[HaloGal[galaxy_index_].GalTreeIndex].lightcone_galaxy_number_in_file_begin = TotLightconeGalCount;       
  save_lightcone_galaxy_append(galaxy_index_, HaloGal[galaxy_index_].SnapNum);
  GalTree[HaloGal[galaxy_index_].GalTreeIndex].lightcone_galaxy_number_in_file_end = TotLightconeGalCount;
#else  /* not defined GALAXYTREE */
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    if(ListOutputSnaps[output_number_] == HaloGal[galaxy_index_].SnapNum)
    { save_lightcone_galaxy_append(galaxy_index_, output_number_); }
#endif /* not defined GALAXYTREE */
#endif /* defined LIGHTCONE_OUTPUT */

#endif /* not defined MCMC */
}


/** @brief pops galaxy from heap (temporary memory) */
void pop_galaxy_from_heap(const int heap_index_)
{
  /* fill the gap in the heap with the galaxy in the last occupied slot */
  const int galaxy_index_ = HaloGalHeap[heap_index_];
  const int last_ = NHaloGal - 1;
  const int last_gal_index_ = HaloGalHeap[last_];
  HaloGalHeap[last_] = galaxy_index_;
  HaloGalHeap[heap_index_] = last_gal_index_;

  /* make sure that the back-pointer of the last galaxy is updated */
  HaloGal[last_gal_index_].HeapIndex = heap_index_;

  NHaloGal--;
}
