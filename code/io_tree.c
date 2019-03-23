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

/** @file   io_tree.c
 *  @date   2016-2019 
 *  @author ?
 *  @author Stefan Hilbert
 *
 *  @brief  Reads in the data from the dark matter simulation merger
 *          trees, creates output files and after calculations frees
 *          the allocated memory.
 *
 * There are three different input files: trees_** - normal tree
 * files; tree_dbids - file containing halo IDs (used if GALAXYTREE
 * ON); tree_aux - file containing the particle data (used if
 * UPDATETYPETWO ON).
 **/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

#ifdef PARALLEL
#ifdef MCMC
#include <mpi.h>
#endif
#endif

/** @file io_tree.c
 *  @brief Reads in the data from the dark matter simulation merger
 *         trees, creates output files and after calculations frees
 *         the allocated memory.
 *
 * There are three different input files: trees_** - normal tree
 * files; tree_dbids - file containing halo IDs (used if GALAXYTREE
 * ON); tree_aux - file containing the particle data (used if
 * UPDATETYPETWO ON).
 */

 
/** @brief Reads tree data from tree files.
 *
 *  Modified versions of myfread/myfwrite/myfseek are called to either
 *  read individual trees from the files into structures or from the
 *  pointers with all the data for the current file into structures.
 *
 *  For each tree the code reads in the header in trees_**: Ntrees -
 *  number of trees in the current file (int); totNHalos - total number
 *  of halos summed over all the trees in current file (int);
 *  TreeNHalos - number of halos in each tree (Ntrees).
 *
 *  Then output files are opened SA_z**_** - for snapshot output;
 *  SA_galtree_** for GALAXYTREE option and SA_**.** for MOMAF.
 *
 *  If UPDATETYPETWO ON the header in tree_aux is also read: NtotHalos -
 *  total number of halos in the file (int); TotIds - total number of
 *  particle IDs (int); Ntrees - total number of trees (int); TotSnaps -
 *  total number of snapshots (int). Define some other quantities:
 *  CountIDs_snaptree - number of Ids for each tree at each snapshot
 *  (TotSnaps * Ntrees); OffsetIDs_snaptree (TotSnaps * Ntrees);
 *  CountIDs_halo - Number of Ids per halo (NtotHalos); OffsetIDs_halo
 *  (int).
 */
void load_tree_table(const int file_number_)
{
  int i_, n_, totNHalos, SnapShotInFileName;
  char file_name_[1024];

#ifdef UPDATETYPETWO
  load_all_auxdata(file_number_);
#endif /* defined UPDATETYPETWO */

#ifdef PARALLEL
#ifndef MCMC
  printf("\nTask %d reading tree file %d\n", ThisTask, file_number_);
#else  /* defined MCMC */
  //if def MCMC and PARALLEL only task 0 reads the representative treefile, then broadcasts
  if(ThisTask==0)
  {
    printf("\nTask %d reading trees_%d\n", ThisTask, file_number_);
#endif /* defined MCMC */
#endif /* defined PARALLEL */

    SnapShotInFileName=LastDarkMatterSnapShot;

#ifdef MCMC
#ifdef MR_PLUS_MRII
    SnapShotInFileName=LastDarkMatterSnapShot_MRII;
#endif /* defined MR_PLUS_MRII */
#endif /* defined MCMC */

#ifdef LOADIDS
#ifndef MRII
    sprintf(file_name_, "%s/treedata/tree_dbids_%03d.%d", SimulationDir, SnapShotInFileName, file_number_);
#else
    sprintf(file_name_, "%s/treedata/tree_sf1_dbids_%03d.%d", SimulationDir, SnapShotInFileName, file_number_);
#endif
    if(!(treedbids_file = fopen(file_name_, "r")))
    {
      char error_message_[2048];
      sprintf(error_message_, "can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }
#endif

#ifndef MRII
    sprintf(file_name_, "%s/treedata/trees_%03d.%d", SimulationDir, SnapShotInFileName, file_number_);
#else
    sprintf(file_name_, "%s/treedata/trees_sf1_%03d.%d", SimulationDir, SnapShotInFileName, file_number_);
#endif

    if(!(tree_file = fopen(file_name_, "r")))
    {
      char error_message_[2048];
      sprintf(error_message_, "can't open file place `%s'\n", file_name_);
      terminate(error_message_);
    }

    //read header on trees_** file
    myfread(&Ntrees, 1, sizeof(int), tree_file);
    myfread(&totNHalos, 1, sizeof(int), tree_file);

    TreeNHalos = mymalloc("TreeNHalos", sizeof(int) * Ntrees);
    TreeFirstHalo = mymalloc("TreeFirstHalo", sizeof(int) * Ntrees);
    TreeNgals[0] = mymalloc("TreeNgals[n_]", NOUT * sizeof(int) * Ntrees);
    for(n_ = 1; n_ < NOUT; n_++)
      TreeNgals[n_] = TreeNgals[n_ - 1] + Ntrees;

    myfread(TreeNHalos, Ntrees, sizeof(int), tree_file);

    if(Ntrees)
      TreeFirstHalo[0] = 0;
    /*Define a variable containing the number you have to jump to
     * get from one firshalo to the next. */
    for(i_ = 1; i_ < Ntrees; i_++)
      TreeFirstHalo[i_] = TreeFirstHalo[i_ - 1] + TreeNHalos[i_ - 1];

#ifdef PRELOAD_TREES
    Halo_Data = mymalloc("Halo_Data", sizeof(struct halo_data) * totNHalos);
    myfseek(tree_file, sizeof(int) * (2 + Ntrees), SEEK_SET);
    myfread(Halo_Data, totNHalos, sizeof(struct halo_data), tree_file);
#ifdef PARALLEL
    printf("\nTask %d done loading trees_%d\n", ThisTask, file_number_);
#endif /* defined PARALLEL */

#ifdef LOADIDS
    HaloIDs_Data = mymalloc("HaloIDs_Data", sizeof(struct halo_ids_data) * totNHalos);
    myfseek(treedbids_file, 0, SEEK_SET);
    myfread(HaloIDs_Data, totNHalos, sizeof(struct halo_ids_data), treedbids_file);
#ifdef PARALLEL
    printf("\nTask %d done loading tree_dbids_%d\n", ThisTask, file_number_);
#endif /* defined PARALLEL */
#endif
#endif

  //if MCMC is turned only Task 0 reads the file and then broadcasts
#ifdef PARALLEL
#ifdef MCMC
  } // end if ThisTask==0

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(&Ntrees,sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&totNHalos,sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask>0)
    TreeNHalos = mymalloc("TreeNHalos", sizeof(int) * Ntrees);

  MPI_Bcast(TreeNHalos,sizeof(int)*Ntrees, MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask>0)
  {
    TreeNgals[0] = mymalloc("TreeNgals[n_]", NOUT * sizeof(int) * Ntrees);
    for(n_ = 1; n_ < NOUT; n_++)
      TreeNgals[n_] = TreeNgals[n_ - 1] + Ntrees;

    TreeFirstHalo = mymalloc("TreeFirstHalo", sizeof(int) * Ntrees);
    if(Ntrees)
      TreeFirstHalo[0] = 0;
    /* Define a variable containing the number you have to jump to
     * get from one firshalo to the next. */
    for(i_ = 1; i_ < Ntrees; i_++)
      TreeFirstHalo[i_] = TreeFirstHalo[i_ - 1] + TreeNHalos[i_ - 1];

    Halo_Data = mymalloc("Halo_Data", sizeof(struct halo_data) * totNHalos);
    HaloIDs_Data = mymalloc("HaloIDs_Data", sizeof(struct halo_ids_data) * totNHalos);
  }

  MPI_Bcast(HaloIDs_Data,totNHalos* sizeof(struct halo_ids_data), MPI_BYTE, 0, MPI_COMM_WORLD);
  
  MPI_Bcast_large_data(Halo_Data, totNHalos, sizeof(struct halo_data), 0, MPI_COMM_WORLD);

//  size_t bytes=totNHalos* sizeof(struct halo_data);
//  int ii, Nmessages=10000;
//  int HaloChunks=1000000;
//
//  //MPI_BCast has a limit of 2Gb so everything needs to be passed in smaller chunks
//  for(ii=0;ii<Nmessages;ii++)
//  {
//    //if next chunk is outside of array size_, just pass whats left and then exit the loop
//    if((ii+1)*HaloChunks>totNHalos)
//    {
//      MPI_Bcast(&Halo_Data[ii*HaloChunks],bytes-ii*HaloChunks* sizeof(struct halo_data), MPI_BYTE, 0, MPI_COMM_WORLD);
//      break;
//    }
//    else
//      MPI_Bcast(&Halo_Data[ii*HaloChunks],HaloChunks* sizeof(struct halo_data), MPI_BYTE, 0, MPI_COMM_WORLD);
//  }

  if(ThisTask==0)
    printf("all tree data has now been broadcasted\n");
#endif
#endif
}


/**@brief Deallocates the arrays used to read in the headers of
 *        trees_** and tree_aux (the ones with more than one
 *        element); if PRELOAD_TREES ON, deallocates
 *        the pointers containing all the input and output data.*/
void free_tree_table(void)
{
#ifdef PRELOAD_TREES
#ifdef LOADIDS
  myfree(HaloIDs_Data);
#endif
  myfree(Halo_Data);
#endif

  myfree(TreeNgals[0]);

  //deallocates header from trees_**
  myfree(TreeFirstHalo);//derived from the header of trees_**
  myfree(TreeNHalos);

#ifdef UPDATETYPETWO
  myfree(TreeAuxData);
#endif

#ifdef LOADIDS
  fclose(treedbids_file);
#endif

  fclose(tree_file);
}


/**@brief Reads the actual trees into structures to be used in the
 *        code: Halo and HaloIDs; the galaxy structures (HaloGal
 *        and Gal) and the HaloAux are also allocated
 *
 *  HaloAux structure is allocated =
 *  (sizeof(struct halo_aux_data) * TreeNHalos[])
 *
 *  Considering the number of halos allocated for the current tree
 *  the size of the structures containing the galaxies
 *  with a halo and the galaxies with and without a halo to be
 *  allocated is estimated:
 *
 *  For galaxies with a halo - MaxGals = MAXGALFAC * TreeNHalos[] and
 *  HaloGal = (sizeof(struct GALAXY) * MaxGals).
 *
 *  For all galaxies - FoF_MaxGals = 10000*15 and
 *  Gal = (sizeof(struct GALAXY) * FoF_MaxGals)
 *
 *  If GALAXYTREE ON, HaloIDs structure is read from tree_dbids =
 *  sizeof(struct halo_ids_data) * TreeNHalos[] */
void load_tree(const int tree_number_)
{
  int i_;

#ifdef PRELOAD_TREES
  Halo = Halo_Data + TreeFirstHalo[tree_number_];
  /*for(i_=0;i_<TreeNHalos[tree_number_];i_++)
          printf("vel=%f\n",Halo[i_].Vel[1]);*/
#ifdef LOADIDS
  HaloIDs = HaloIDs_Data + TreeFirstHalo[tree_number_];
#endif
#else

  Halo = mymalloc("Halo", sizeof(struct halo_data) * TreeNHalos[tree_number_]);
  myfseek(tree_file, sizeof(int) * (2 + Ntrees) + sizeof(struct halo_data) * TreeFirstHalo[tree_number_], SEEK_SET);
  myfread(Halo, TreeNHalos[tree_number_], sizeof(struct halo_data), tree_file);
#ifdef LOADIDS
  HaloIDs = mymalloc("HaloIDs", sizeof(struct halo_ids_data) * TreeNHalos[tree_number_]);
  myfseek(treedbids_file, sizeof(struct halo_ids_data) * TreeFirstHalo[tree_number_], SEEK_SET);
  myfread(HaloIDs, TreeNHalos[tree_number_], sizeof(struct halo_ids_data), treedbids_file);
#endif

#endif

  //Allocate HaloAux and Galaxy structures.
  HaloAux = mymalloc("HaloAux", sizeof(struct halo_aux_data) * TreeNHalos[tree_number_]);

  for(i_ = 0; i_ < TreeNHalos[tree_number_]; i_++)
  {
    HaloAux[i_].DoneFlag = 0;
    HaloAux[i_].HaloFlag = 0;
    HaloAux[i_].NGalaxies = 0;
  }

  if(AllocValue_MaxHaloGal == 0)
    AllocValue_MaxHaloGal = 1 + TreeNHalos[tree_number_] / (0.25 * (LastDarkMatterSnapShot+1));

  if(AllocValue_MaxGal == 0)
    AllocValue_MaxGal = 2000;

  MaxHaloGal = AllocValue_MaxHaloGal;
  NHaloGal = 0;
  HaloGal = mymalloc_movable(&HaloGal, "HaloGal", sizeof(struct GALAXY) * MaxHaloGal);
  HaloGalHeap = mymalloc_movable(&HaloGalHeap, "HaloGalHeap", sizeof(int) * MaxHaloGal);
  for(i_ = 0; i_ < MaxHaloGal; i_++)
    HaloGalHeap[i_] = i_;

  MaxGal = AllocValue_MaxGal;
  Gal = mymalloc_movable(&Gal, "Gal", sizeof(struct GALAXY) * MaxGal);

#ifdef GALAXYTREE
  if(AllocValue_MaxGalTree == 0)
    AllocValue_MaxGalTree = 1.5 * TreeNHalos[tree_number_];

  MaxGalTree = AllocValue_MaxGalTree;
  GalTree = mymalloc_movable(&GalTree, "GalTree", sizeof(struct galaxy_tree_data) * MaxGalTree);
#endif
}


/**@brief Frees all the Halo and Galaxy structures in the code. */
void free_galaxies_and_tree(void)
{
#ifdef GALAXYTREE
  myfree(GalTree);
#endif
  myfree(Gal);
  myfree(HaloGalHeap);
  myfree(HaloGal);
  myfree(HaloAux);

#ifndef PRELOAD_TREES
#ifdef LOADIDS
  myfree(HaloIDs);
#endif
  myfree(Halo);
#endif
}
