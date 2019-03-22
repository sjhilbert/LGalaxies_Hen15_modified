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

/** @file sam.c
 * @brief SAM Construct Galaxies, Join Galaxies of progenitors, Evolve Galaxies
 *
 * SAM is to be called for each of the chosen dark matter tree files and
 * for each tree, on each treefile: reads tree, constructs the galaxies, saves
 * the galaxies, frees memory for galaxies and tree.
 * */

 /* own header for consistency: */
#include "sam.h"
 
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



/**@brief SAM() loops on trees and calls construct_galaxies.*/
double SAM(const int filenr)
{
  int treenr, halonr;

#ifdef MCMC
  int ii;
  MCMC_GAL = mymalloc("MCMC_Gal", sizeof(struct MCMC_GALAXY) * MCMCAllocFactor);
  for(ii=0;ii<NOUT;ii++)
  { TotMCMCGals[ii] = 0; }

#ifdef MR_PLUS_MRII
  change_dark_matter_sim("MR");
#else /* not defined MR_PLUS_MRII */
  if(CurrentMCMCStep==1)
  { read_sample_info(); }
#ifdef HALOMODEL
  else
  {
    int snap, ii;
    for(snap=0;snap<NOUT;snap++)
      for(ii=0;ii<NFofsInSample[snap];ii++)
        MCMC_FOF[ii].NGalsInFoF[snap]=0;
  }
#endif /* defined HALOMODEL */
#endif /* not defined MR_PLUS_MRII */
#endif /* defined MCMC */

  //to be used when we have tables for the scaling in any cosmology
  //read_scaling_parameters();

#ifndef MCMC

#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE
  create_galaxy_tree_file(filenr);
#else /* not defined GALAXYTREE */
  create_galaxy_files(filenr);
#endif /* not defined GALAXYTREE */
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */
  
#ifdef LIGHTCONE_OUTPUT
  create_lightcone_galaxy_files(filenr);
#endif /* defined LIGHTCONE_OUTPUT */ 
  
#endif /* not defined MCMC */
  
#ifndef LIGHTCONE_OUTPUT_ONLY
#ifdef GALAXYTREE /* defined GALAXYTREE */
  FILE *fdg = fopen("treengal.dat", "w");
#endif /* defined GALAXYTREE */
#endif /* defined LIGHTCONE_OUTPUT */ 

//***************************************************************************************
//***************************************************************************************

  //for(treenr = 0; treenr < NTrees_Switch_MR_MRII; treenr++)
  for(treenr = 0; treenr < Ntrees; treenr++)
  {
  //printf("doing tree %d of %d\n", treenr, Ntrees);
#ifdef MR_PLUS_MRII
    if(treenr == NTrees_Switch_MR_MRII)
    { change_dark_matter_sim("MRII"); }
#endif /* defined MR_PLUS_MRII */

    load_tree(treenr);
#ifdef MCMC
#ifdef PRELOAD_TREES
    if(CurrentMCMCStep==1)
#endif /* defined PRELOAD_TREES */
#endif /* defined MCMC */
    { scale_cosmology(TreeNHalos[treenr]); }

    gsl_rng_set(random_generator, filenr * 100000 + treenr);
    NumMergers = 0;
    NHaloGal = 0;
#ifdef GALAXYTREE
    NGalTree = 0;
    IndexStored = 0;
#endif/* defined GALAXYTREE */
    int snapnum;
    //LastSnapShotNr is the highest output snapshot
    /* we process the snapshots now in temporal order 
      * (as a means to reduce peak memory usage) */
    for(snapnum = 0; snapnum <= LastSnapShotNr; snapnum++)
    {
#ifdef MCMC
      /* read the appropriate parameter list for current snapnum
       * into the parameter variables to be used in construct_galaxies */
      read_mcmc_par(snapnum);
#ifdef HALOMODEL
      //because we need halo masses even for FOFs
      //with no galaxies it needs to be done here
      assign_FOF_masses(snapnum, treenr);
#endif /* defined HALOMODEL */
#endif /* defined MCMC */
      for(halonr = 0; halonr < TreeNHalos[treenr]; halonr++)
        if(HaloAux[halonr].DoneFlag == 0 && Halo[halonr].SnapNum == snapnum)
          construct_galaxies(treenr, halonr);
    }

    /* output remaining galaxies as needed */
    while(NHaloGal)
    {
#ifdef LIGHTCONE_OUTPUT
      lightcone_N_galaxies_remaining_for_output_past_construct_galaxies++;
#endif /* defined LIGHTCONE_OUTPUT */              
      output_galaxy(treenr, 0);
    }

#ifndef MCMC
#ifdef GALAXYTREE
    update_galaxy_tree_ids();
#ifndef LIGHTCONE_OUTPUT_ONLY
    save_galaxy_tree_finalize(filenr, treenr);
#endif /*not defined LIGHTCONE_OUTPUT_ONLY */
#ifdef LIGHTCONE_OUTPUT
    save_lightcone_galaxy_finalize(filenr, treenr);
#endif /* defined LIGHTCONE_OUTPUT */
   
#ifndef PARALLEL
#ifdef LIGHTCONE_OUTPUT
    if((treenr/100)*100==treenr) printf("treenr=%d  TotGalCount=%d   TotLightconeGalCount=%lld\n", treenr, TotGalCount, TotLightconeGalCount);
#else /* not defined LIGHTCONE_OUTPUT */
    if((treenr/100)*100==treenr) printf("treenr=%d  TotGalCount=%d\n", treenr, TotGalCount);
#endif /* not defined LIGHTCONE_OUTPUT */
#endif /* note defined PARALLEL*/
    fflush(stdout);
#ifndef LIGHTCONE_OUTPUT_ONLY
    fprintf(fdg, "%d\n", NGalTree);
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */
#endif /* defined GALAXYTREE */
#else /* defined MCMC */
#endif /* defined MCMC */
    free_galaxies_and_tree();
  }//loop on trees

#ifdef MCMC
  double lhood = get_likelihood();

#ifdef MR_PLUS_MRII
  free(MCMC_FOF);
#else /* not defined MR_PLUS_MRII */
  if(CurrentMCMCStep==ChainLength)
  { free(MCMC_FOF); }
#endif /* not defined MR_PLUS_MRII */

  myfree(MCMC_GAL);
  return lhood;

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
void construct_galaxies(const int treenr, const int halonr)
{
  static int halosdone = 0;
  int prog, fofhalo, ngal, cenngal;
  
#ifdef LIGHTCONE_OUTPUT_ONLY
#ifdef LIGHTCONE_MAY_SKIP_CONSTRUCT_GALAXY
  /* if halos is outside lightcone, don't construct, but return early */
  if(is_outside_lightcone_for_snapshot[Halo[halonr].SnapNum])
  {
    lightcone_N_galaxies_skipped_construction++;
    return;  
  }
  
  if(check_outside_lightcone_for_snapshot[Halo[halonr].SnapNum])
  {
    float d_0_ = Halo[halonr].Pos[0] - lightcone_observer_position[0]; d_0_ = wrap(d_0_, BoxSize);
    float d_1_ = Halo[halonr].Pos[1] - lightcone_observer_position[1]; d_1_ = wrap(d_1_, BoxSize);
    float d_2_ = Halo[halonr].Pos[2] - lightcone_observer_position[2]; d_2_ = wrap(d_2_, BoxSize);

    if(d_0_ * d_0_ + d_1_ * d_1_ + d_2_ * d_2_ > pow2(lightcone_radius_for_snapshot[Halo[halonr].SnapNum]))
    {
      lightcone_N_galaxies_skipped_construction++;
      return;
    }
  }
#endif /* defined LIGHTCONE_MAY_SKIP_CONSTRUCT_GALAXY */
#endif /* defined LIGHTCONE_OUTPUT_ONLY */

  HaloAux[halonr].DoneFlag = 1;
  halosdone++;

  prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0) //If halo has a progenitor
  {
    if(HaloAux[prog].DoneFlag == 0) //If progenitor hasn't been done yet
    { construct_galaxies(treenr, prog); }
    prog = Halo[prog].NextProgenitor;         //Jump to next halo in progenitors FOF
  }

  //Now check for the progenitors of all the halos in the current FOF group
  fofhalo = Halo[halonr].FirstHaloInFOFgroup;   //Starting at the first halo in current FOF
  if(HaloAux[fofhalo].HaloFlag == 0)            //If it hasn't been done
  {
    HaloAux[fofhalo].HaloFlag = 1;            //mark as to do
    while(fofhalo >= 0)       //go through all the halos in current FOF
    {
      prog = Halo[fofhalo].FirstProgenitor;
      while(prog >= 0)              //build its progenitors
      {
        if(HaloAux[prog].DoneFlag == 0)
        { construct_galaxies(treenr, prog); }
        prog = Halo[prog].NextProgenitor;
      }

      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;   //Jump to next halo in FOF
    }
  }

  /* At this point, the galaxies for all progenitors of this halo have been
   * properly constructed. Also, the galaxies of the progenitors of all other 
   * halos in the same FOF-group have been constructed as well. We can hence go
   * ahead and construct all galaxies for the subhalos in this FOF halo, and
   * evolve them in time. */

  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  if(HaloAux[fofhalo].HaloFlag == 1)            //If it is marked as an halo to do
  {
    ngal = 0;
    HaloAux[fofhalo].HaloFlag = 2;

    cenngal = get_merger_center(fofhalo);     //Find type 0 for type 1 to merge into

    /*For all the halos in the current FOF join all the progenitor galaxies together
      * ngals will be the total number of galaxies in the current FOF*/
    while(fofhalo >= 0)
    {
      ngal = join_galaxies_of_progenitors(fofhalo, ngal, &cenngal);
      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }

    /*Evolve the Galaxies -> SAM! */
    evolve_galaxies(Halo[halonr].FirstHaloInFOFgroup, ngal, treenr, cenngal);

    int p;
    for (p =0;p<ngal;p++)
    { mass_checks("Construct_galaxies #1",p); }
  }
}


/** @brief join_galaxies_of_progenitors() updates the properties of the
 *         galaxy from the dark matter halo properties and deals with
 *         merging clocks. This routine is called by construct_galaxies
 *         for every halo in the FOF being constructed. When there is no
 *         galaxy in the Halo of FirstProgenitor, the first_occupied
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
int join_galaxies_of_progenitors(const int halonr, const int ngalstart, int *cenngal)
{
  int ngal, i, j, prog, first_occupied, lenmax, centralgal, mostmassive;

  /* When there is no galaxy in the Halo of FirstProgenitor, the first_occupied
   * pointer is changed to a subhalo which have the maximum mass (This should
   * only happen in the case that the leaf on the firstprogenitor branch occurs
   * as a subhalo, in that case no galaxy would be assigned to it). */
  prog = Halo[halonr].FirstProgenitor;
  first_occupied = Halo[halonr].FirstProgenitor;
  lenmax = 0;
  if(prog >= 0)                                 //If halo has progenitors
  {
    if(HaloAux[prog].NGalaxies == 0)          //if progenitor has no galaxies
      while(prog >= 0)
      {
        int currentgal;
        for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++)
        {
          if(HaloGal[currentgal].Type == 0 || HaloGal[currentgal].Type == 1)
          {
            if(Halo[prog].Len > lenmax)
            {
              lenmax = Halo[prog].Len;
              first_occupied = prog;          //define the new first_occupied
            }
          }
          currentgal = HaloGal[currentgal].NextGalaxy;
        }
        prog = Halo[prog].NextProgenitor;
      }
  }
   
  lenmax = 0;
  prog = Halo[halonr].FirstProgenitor;
  mostmassive = Halo[halonr].FirstProgenitor;

  /* loop through all the progenitors and get the halo mass and ID
   * of the most massive*/
  while(prog >= 0)
  {
    if(Halo[prog].Len > lenmax)
    {
      lenmax = Halo[prog].Len;
      mostmassive = prog;
    }
    prog = Halo[prog].NextProgenitor;
  }

  ngal = ngalstart;
  prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0)
  {
    int currentgal;
    for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++)
    {
      if(ngal >= MaxGal)
      {
        AllocValue_MaxGal *= ALLOC_INCREASE_FACTOR;
        MaxGal = AllocValue_MaxGal;
        if(MaxGal<ngal+1) MaxGal=ngal+1;
        Gal = myrealloc_movable(Gal, sizeof(struct GALAXY) * MaxGal);
      }
      if(*cenngal==currentgal)
      { *cenngal=ngal; }

      /* Copy galaxy properties from progenitor,
        * except for those that need initialising */
      Gal[ngal] = HaloGal[currentgal];

      Gal[ngal].HaloNr = halonr;
      Gal[ngal].CoolingRadius = 0.0;
      Gal[ngal].CoolingGas = 0.0;

      Gal[ngal].PrimordialAccretionRate = 0.0;
      Gal[ngal].CoolingRate = 0.0;
      Gal[ngal].CoolingRate_beforeAGN = 0.0;
      Gal[ngal].Sfr = 0.0;
      Gal[ngal].SfrBulge = 0.0;
      Gal[ngal].QuasarAccretionRate=0.0;
      Gal[ngal].RadioAccretionRate=0.0;
#ifdef GALAXYTREE
      Gal[ngal].FirstProgGal = HaloGal[currentgal].GalTreeIndex;    /* CHECK */
#endif /* defined GALAXYTREE */
      // To fail this check means that we copy in a failed galaxy
      mass_checks("Middle of join_galaxies_of_progenitors", ngal);

      /* Update Properties of this galaxy with physical properties of halo */
      /* this deals with the central galaxies of subhalos */
      if(Gal[ngal].Type == 0 || Gal[ngal].Type == 1)
      {
        if(prog == first_occupied)
        {
#ifdef HALOPROPERTIES
          Gal[ngal].HaloM_Mean200 = Halo[halonr].M_Mean200;
          Gal[ngal].HaloM_Crit200 = Halo[halonr].M_Crit200;
          Gal[ngal].HaloM_TopHat = Halo[halonr].M_TopHat;
          Gal[ngal].HaloVelDisp = Halo[halonr].VelDisp;
          Gal[ngal].HaloVmax = Halo[halonr].Vmax;
#endif /* defined HALOPROPERTIES */
          Gal[ngal].MostBoundID = Halo[halonr].MostBoundID;
          for(j = 0; j < 3; j++)
          {
            Gal[ngal].Pos[j] = Halo[halonr].Pos[j];
            Gal[ngal].Vel[j] = Halo[halonr].Vel[j];
#ifdef HALOPROPERTIES
            Gal[ngal].HaloPos[j] = Halo[halonr].Pos[j];
            Gal[ngal].HaloVel[j] = Halo[halonr].Vel[j];
#endif /* defined HALOPROPERTIES */
          }

          Gal[ngal].Len = Halo[halonr].Len;

          // FOFCentralGal property in case that is different from FirstGalaxy
          if(halonr == Halo[halonr].FirstHaloInFOFgroup)
          { update_centralgal(ngal, halonr); }
          else
          { update_type_1(ngal, halonr, prog); }

          if(DiskRadiusModel == 1 || DiskRadiusModel == 2)
          {
            Gal[ngal].GasDiskRadius     = get_disk_radius(halonr, ngal);
            Gal[ngal].StellarDiskRadius = Gal[ngal].GasDiskRadius;
          }
          Gal[ngal].Vmax = Halo[halonr].Vmax;
        }
        else //type 2 galaxies
        { update_type_2(ngal, halonr, prog, mostmassive); }
      }

      /* Note: Galaxies that are already type=2 do not need a special treatment at this point */
      if(Gal[ngal].Type < 0 || Gal[ngal].Type > 2)
      { terminate("Unknown galaxy type\n"); }

      ngal++;
      currentgal = HaloGal[currentgal].NextGalaxy;
    }
    prog = Halo[prog].NextProgenitor;
  }

  /* If there are no progenitors with galaxies, a new galaxy is created.
   * However, if it's a subhalo, no galaxy is placed, since it would stay
   * at zero luminosity. */
  if(ngal == 0)
  {
    *cenngal=0;
    if(Halo[halonr].FirstHaloInFOFgroup == halonr)
    {
      init_galaxy(ngal, halonr);
      ngal++;
    }
  }

  /* satelites (type 2's) will preferably merge onto this type 1 rather than the type 0 */
  for(i = ngalstart, centralgal = -1; i < ngal; i++)
    if(Gal[i].Type == 0 || Gal[i].Type == 1)
    {
      if(centralgal != -1)
      { terminate("Subhalo hosts more than one Type 0/1\n"); }

      centralgal = i;
    }

  for(i = ngalstart; i < ngal; i++)
  {
    Gal[i].CentralGal = centralgal;
    if(centralgal != -1)
      for(j = 0; j < 3; j++)
      { Gal[i].MergCentralPos[j] = Gal[centralgal].Pos[j]; }
  }

  /* Satellites whose type 1 has merged into type 0, will be reset to merge
   * into the type 0. */
  if(centralgal == -1 && ngal != ngalstart)
  {
    for(i = ngalstart; i < ngal; i++)
    {
      Gal[i].CentralGal = *cenngal;
      for(j = 0; j < 3; j++)
      { Gal[i].MergCentralPos[j] = Gal[*cenngal].Pos[j]; }
    }
  }
    
  for (i = ngalstart; i<ngal; i++)
  { mass_checks("Bottom of join_galaxies_of_progenitors",i); }

  report_memory_usage(&HighMark, "join_galaxies");

  return ngal;
}


/** @brief evolve_galaxies() deals with most of the SA recipes. This is
  *        where most of the physical recipes are called, including:
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
  * Note: halonr is here the FOF-background subhalo (i.e. main halo) 
  */
void evolve_galaxies(const int halonr, const int ngal, const int treenr, const int cenngal)
{
  int p, q, nstep, centralgal, merger_centralgal, currenthalo, prevgal, /* start,*/ i;
  double infallingGas, deltaT, Zcurr;
  double time, previoustime, newtime;

#ifdef STAR_FORMATION_HISTORY
  double age_in_years;
#endif /* defined STAR_FORMATION_HISTORY */

//   // Eddington time in code units
//   // code units are UnitTime_in_s/Hubble_h
//   const double t_Edd=1.42e16*Hubble_h/UnitTime_in_s;

  //previoustime = NumToTime(Gal[0].SnapNum);
  previoustime = NumToTime(Halo[halonr].SnapNum-1);
  newtime = NumToTime(Halo[halonr].SnapNum);

  /* Time between snapshots */
  deltaT = previoustime - newtime;
  /* Redshift of current Snapnum */
  Zcurr = ZZ[Halo[halonr].SnapNum];

  centralgal = Gal[0].CentralGal;

  for (p =0;p<ngal;p++)
  { mass_checks("Evolve_galaxies #0",p); }

  //print_galaxy("\n\ncheck1", centralgal, halonr);

  if(Gal[centralgal].Type != 0 || Gal[centralgal].HaloNr != halonr)
  { terminate("Something wrong here ..... \n"); }

  /* Update all galaxies to same star-formation history time-bins.
   * Needed in case some galaxy has skipped a snapshot. */
#ifdef STAR_FORMATION_HISTORY
  age_in_years=(Age[0]-previoustime)*UnitTime_in_years/Hubble_h; //ROB: age_in_years is in units of "real years"!
  nstep=0;
  for (p=0; p<ngal; p++)
  { sfh_update_bins(p,Halo[halonr].SnapNum-1,nstep,age_in_years); }
#endif /* defined STAR_FORMATION_HISTORY */

  /* Handle the transfer of mass between satellites and central galaxies */
  deal_with_satellites(centralgal, ngal);

  /* Delete inconsequential galaxies */
  for (p =0;p<ngal;p++)
  {
    if (Gal[p].Type ==2 && Gal[p].ColdGas+Gal[p].DiskMass+Gal[p].BulgeMass <1.e-8)
    { Gal[p].Type = 3; }
    else
    { mass_checks("Evolve_galaxies #0.1",p); }
  }
   
  /* Calculate how much hot gas needs to be accreted to give the correct baryon fraction
   * in the main halo. This is the universal fraction, less any reduction due to reionization. */
  infallingGas = infall_recipe(centralgal, ngal, Zcurr);
  Gal[centralgal].PrimordialAccretionRate=infallingGas/deltaT;
  
  /* All the physics are computed in a number of intervals between snapshots
   * equal to STEPS */
  for (nstep = 0; nstep < STEPS; nstep++)
  {
    /* time to present of the current step */
    time = previoustime - (nstep + 0.5) * (deltaT / STEPS);

    /* Update all galaxies to the star-formation history time-bins of current step*/
#ifdef STAR_FORMATION_HISTORY
    age_in_years=(Age[0]-time)*UnitTime_in_years/Hubble_h;
    for (p=0; p<ngal; p++)
    { sfh_update_bins(p,Halo[halonr].SnapNum-1,nstep,age_in_years); }

#endif /* defined STAR_FORMATION_HISTORY */

    /* Infall onto central galaxy only, if required to make up a baryon deficit */
#ifndef GUO10
#ifndef GUO13
    if (infallingGas > 0.)
#endif /* not defined GUO13 */
#endif /* not defined GUO10 */
    { add_infall_to_hot(centralgal, infallingGas / STEPS); }

    mass_checks("Evolve_galaxies #0.5",centralgal);

    for (p = 0; p < ngal; p++)
    {
      /* don't treat galaxies that have already merged */
      if(Gal[p].Type == 3)
      { continue; }
      
      mass_checks("Evolve_galaxies #1",p);

      if (Gal[p].Type == 0 || Gal[p].Type == 1)
      {
        reincorporate_gas(p, deltaT / STEPS);
        /* determine cooling gas given halo properties and add it to the cold phase*/
        mass_checks("Evolve_galaxies #1.5",p);
        compute_cooling(p, deltaT / STEPS);
      }
    }

    //this must be separated as now satellite AGN can heat central galaxies
    //therefore the AGN from all satellites must be computed, in a loop inside this function,
    //before gas is cooled into central galaxies (only suppress cooling, the gas is not actually heated)
    if(AGNRadioModeModel != 5)
    { do_AGN_heating(deltaT / STEPS, ngal); }

    for (p = 0; p < ngal; p++)
    {
      cool_gas_onto_galaxy(p, deltaT / STEPS);
      mass_checks("Evolve_galaxies #2",p);
      starformation(p, centralgal, time, deltaT / STEPS, nstep);
      mass_checks("Evolve_galaxies #3",p);
      //print_galaxy("check3", centralgal, halonr);
    }

    /* Check for merger events */
    for(p = 0; p < ngal; p++)
    {
      if(Gal[p].Type == 2 || (Gal[p].Type == 1 && Gal[p].MergeOn == 1))             /* satellite galaxy */
      {
        Gal[p].MergTime -= deltaT / STEPS;
        if(Gal[p].MergTime < 0.0)
        {
          NumMergers++;

          if(Gal[p].Type == 1)
            for(q = 0; q < ngal; q++)
              if(Gal[q].Type == 2 && Gal[p].CentralGal == p)
              { Gal[q].CentralGal = cenngal; }

          if(Gal[p].Type == 2)
          { merger_centralgal = Gal[p].CentralGal; }
          else
          { merger_centralgal = cenngal; }
          mass_checks("Evolve_galaxies #4",p);
          mass_checks("Evolve_galaxies #4",merger_centralgal);
          mass_checks("Evolve_galaxies #4",centralgal);
          deal_with_galaxy_merger(p, merger_centralgal, centralgal, time, deltaT, nstep);
          mass_checks("Evolve_galaxies #5",p);
          mass_checks("Evolve_galaxies #5",merger_centralgal);
          mass_checks("Evolve_galaxies #5",centralgal);
        }
      }
    }//loop on all galaxies to detect mergers

#ifdef DETAILED_METALS_AND_MASS_RETURN
    //DELAYED ENRICHMENT AND MASS RETURN + FEEDBACK: No fixed yield or recycling fraction anymore. FB synced with enrichment
    for (p = 0; p < ngal; p++)
    { update_yields_and_return_mass(p, centralgal, deltaT/STEPS, nstep); }
#endif /* defined DETAILED_METALS_AND_MASS_RETURN */

  }/* end move forward in interval STEPS */

  for(p = 0; p < ngal; p++)
  {
    if(Gal[p].Type == 2)
    {
#ifndef UPDATETYPETWO
      int jj;
      float tmppos;
      for(jj = 0; jj < 3; jj++)
      {
        tmppos = wrap(Gal[p].DistanceToCentralGal[jj],BoxSize);
        tmppos *=  2.*sqrt(Gal[p].MergTime/Gal[p].OriMergTime);
        Gal[p].Pos[jj] = Gal[p].MergCentralPos[jj] + tmppos;

        if(Gal[p].Pos[jj] < 0)
                  Gal[p].Pos[jj] = BoxSize + Gal[p].Pos[jj];
        if(Gal[p].Pos[jj] > BoxSize)
                  Gal[p].Pos[jj] = Gal[p].Pos[jj] - BoxSize;
      }
#endif /* not defined UPDATETYPETWO */
      /* Disruption of type 2 galaxies. Type 1 galaxies are not disrupted since usually
        * bayonic component is more compact than dark matter.*/

      if(DisruptionModel==0)
      { disrupt(p); }
    }
  }
    
  for (p =0;p<ngal;p++)
  { mass_checks("Evolve_galaxies #6",p); }

  
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef  POST_PROCESS_MAGS
  int n;
  /* If this is an output snapshot apply the dust model to each galaxy */
  for(n = 0; n < NOUT; n++)
    {
      if(Halo[halonr].SnapNum == ListOutputSnaps[n])
      {
        for(p = 0; p < ngal; p++)
        { dust_model(p, n, halonr); }
        break;
      }
    }
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

  /* now save the galaxies of all the progenitors (and free the associated storage) */
  int prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0)
  {
    int currentgal;
    for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++)
    {
      int nextgal = HaloGal[currentgal].NextGalaxy;
      /* this will write this galaxy to an output file and free the storage associate with it */
      output_galaxy(treenr, HaloGal[currentgal].HeapIndex);
      currentgal = nextgal;
    }
    prog = Halo[prog].NextProgenitor;
  }
    
  /** @todo check why start was set here, and why it was set to start = NGalTree */
  for(p = 0, prevgal = -1, currenthalo = -1, centralgal = -1 /*, start = NGalTree */; p < ngal; p++)
  {
    if(Gal[p].HaloNr != currenthalo)
    {
      currenthalo = Gal[p].HaloNr;
      HaloAux[currenthalo].FirstGalaxy = -1;
      HaloAux[currenthalo].NGalaxies = 0;
    }

    mass_checks("Evolve_galaxies #7",p);

    if(Gal[p].Type != 3)
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
        for(i = oldmax; i < MaxHaloGal; i++)
        { HaloGalHeap[i] = i; }
      }

      Gal[p].SnapNum = Halo[currenthalo].SnapNum;

#ifndef GUO10
#ifdef UPDATETYPETWO
      update_type_two_coordinate_and_velocity(treenr, p, Gal[0].CentralGal);
#endif /* defined UPDATETYPETWO */
#endif /* not defined GUO10 */

      /* when galaxies are outputed, the slot is filled with the
        * last galaxy in the heap. New galaxies always take the last spot */
      int nextgal = HaloGalHeap[NHaloGal];
      HaloGal[nextgal] = Gal[p];
      HaloGal[nextgal].HeapIndex = NHaloGal;

      if(HaloAux[currenthalo].FirstGalaxy < 0)
      { HaloAux[currenthalo].FirstGalaxy = nextgal; }

      if(prevgal >= 0)
      { HaloGal[prevgal].NextGalaxy = nextgal; }
      prevgal = nextgal;

      HaloAux[currenthalo].NGalaxies++;
      NHaloGal++;


#ifdef GALAXYTREE
      if(NGalTree >= MaxGalTree)
      {
        AllocValue_MaxGalTree *= ALLOC_INCREASE_FACTOR;
        MaxGalTree = AllocValue_MaxGalTree;
        if(MaxGalTree<NGalTree+1) MaxGalTree=NGalTree+1;
        GalTree = myrealloc_movable(GalTree, sizeof(struct galaxy_tree_data) * MaxGalTree);
      }
      HaloGal[nextgal].GalTreeIndex = NGalTree;

      memset(&GalTree[NGalTree], 0, sizeof(struct galaxy_tree_data));
      GalTree[NGalTree].HaloGalIndex = nextgal;
      GalTree[NGalTree].SnapNum = Halo[currenthalo].SnapNum;
      GalTree[NGalTree].NextProgGal = -1;
      GalTree[NGalTree].DescendantGal = -1;

      GalTree[NGalTree].FirstProgGal = Gal[p].FirstProgGal;
      if(Gal[p].Type == 0)
      { centralgal = NGalTree; }
      NGalTree++;
#endif /* defined GALAXYTREE */
    }
  }

#ifdef GALAXYTREE
/** @bug  possible BUG: start was set to start = NGalTree above,
 *        start seemed to be never updated after that, and below, 
 *        start is used to init loop with p = start and termination p < NGalTree,
 *        so the loop seems to never run.
 *        now I put  start = NGalTree here to reproduce current (buggy?) behaviour
 */
    int start = NGalTree;
    for(p = start; p < NGalTree; p++)
    {
      if(centralgal < 0)
      { terminate("centralgal < 0"); }
      GalTree[p].FOFCentralGal = centralgal;
    }
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
void output_galaxy(const int treenr, const int heap_index)
{
  (void) treenr; /* avoid unused-parameter warning for certain compiler options */
  int gal_index = HaloGalHeap[heap_index];

  if(heap_index >= NHaloGal)
  { terminate("heap_index >= NHaloGal"); }
 
  if(HaloGal[gal_index].HeapIndex != heap_index)                // consistency check
  { terminate("this should not happen"); }

#ifdef GUO10
#ifdef UPDATETYPETWO
  update_type_two_coordinate_and_velocity(treenr, gal_index, HaloGal[0].CentralGal);
#endif
#endif

#ifdef GALAXYTREE
  GalTree[HaloGal[gal_index].GalTreeIndex].IndexStored = IndexStored++;
#ifndef LIGHTCONE_OUTPUT_ONLY  
  save_galaxy_tree_append(gal_index);
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */
#else /* not defined GALAXYTREE */
#ifdef MCMC
  /* if MCMC galaxies are saved into a structure to then be compared
   * with observations. The output will come in form of a file with
   * parameters, not galaxy properties */

  save_galaxy_for_mcmc(gal_index);

#else  /* not defined MCMC */
#ifndef LIGHTCONE_OUTPUT_ONLY  
  {
    int n;
    for(n = 0; n < NOUT; n++)
      if(ListOutputSnaps[n] == HaloGal[gal_index].SnapNum)
      {
        save_galaxy_append(treenr, gal_index, n);
      }
  }
#endif /* not defined LIGHTCONE_OUTPUT_ONLY */
#endif /* not defined MCMC */
#endif /* not defined GALAXYTREE */

#ifndef MCMC
#ifdef LIGHTCONE_OUTPUT
#ifdef GALAXYTREE
  GalTree[HaloGal[gal_index].GalTreeIndex].lightcone_galaxy_number_in_file_begin = TotLightconeGalCount;       
  save_lightcone_galaxy_append(gal_index, HaloGal[gal_index].SnapNum);
  GalTree[HaloGal[gal_index].GalTreeIndex].lightcone_galaxy_number_in_file_end = TotLightconeGalCount;
#else  /* not defined GALAXYTREE */
  {
    int n;
    for(n = 0; n < NOUT; n++)
      if(ListOutputSnaps[n] == HaloGal[gal_index].SnapNum)
      {
        save_lightcone_galaxy_append(gal_index, n);
      }
  }
#endif /* not defined GALAXYTREE */

#endif /* defined LIGHTCONE_OUTPUT */
#endif /* not defined MCMC */ 

  /* fill the gap in the heap with the galaxy in the last occupied slot */

  int last = NHaloGal - 1;
  int last_gal_index = HaloGalHeap[last];
  HaloGalHeap[last] = gal_index;
  HaloGalHeap[heap_index] = last_gal_index;

  /* make sure that the back-pointer of the last galaxy is updated */
  HaloGal[last_gal_index].HeapIndex = heap_index;

  NHaloGal--;
}
