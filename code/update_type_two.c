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

/** @file   update_type_two.c
 *  @date   2016-2019
 *  @author ???
 *  @author Stefan Hilbert
 *
 *  @brief  
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"

#ifdef PARALLEL
#ifdef MCMC
#include <mpi.h>
#endif
#endif

#ifdef UPDATETYPETWO

void update_type_two_coordinate_and_velocity(const int tree_number_, const int galaxy_number_, const int central_galaxy_number_)
{
  int j_;
  float tmp_pos_;
  double Scale_V_, dv_;
  const int p = galaxy_number_;
//printf("updating type 2 treenr =%d\n",tree_number_);
#ifdef GUO10
  if(HaloGal[galaxy_number_].Type == 2)  /* Update positions of type 2's */
#else
  if(Gal[galaxy_number_].Type == 2)  /* Update positions of type 2's */
#endif
  {
#ifdef GUO10
    int snapshot_number_ = HaloGal[galaxy_number_].SnapNum;
#else
    int snapshot_number_ = Gal[galaxy_number_].SnapNum;
#endif
    Nids = CountIDs_snaptree[snapshot_number_ * Ntrees + tree_number_];
    OffsetIDs = OffsetIDs_snaptree[snapshot_number_ * Ntrees + tree_number_];

    size_t header_offset = 4 * sizeof(int) + 2 * TotSnaps * sizeof(int) + 2 * TotSnaps * Ntrees * sizeof(int) + 2 * sizeof(int) * NtotHalos;

    IdList  = (long long *) (TreeAuxData + header_offset);
    PosList = (float     *) (TreeAuxData + header_offset + TotIds * sizeof(long long));
    VelList = (float     *) (TreeAuxData + header_offset + TotIds * sizeof(long long) + TotIds * 3 * sizeof(float));

    IdList  +=     OffsetIDs;
    PosList += 3 * OffsetIDs;
    VelList += 3 * OffsetIDs;
#ifdef GUO10
    get_coordinates(HaloGal[galaxy_number_].Pos, HaloGal[galaxy_number_].Vel, HaloGal[galaxy_number_].MostBoundID, tree_number_, HaloGal[galaxy_number_].HaloNr, HaloGal[galaxy_number_].SnapNum);
#else
    get_coordinates(Gal[galaxy_number_].Pos, Gal[galaxy_number_].Vel, Gal[galaxy_number_].MostBoundID, tree_number_, Gal[galaxy_number_].HaloNr, Gal[galaxy_number_].SnapNum);
#endif

//#ifdef SCALE_COSMOLOGY
#ifdef GUO10
    for(j_ = 0; j_ < 3; j_++)
      HaloGal[galaxy_number_].Pos[j_] *= ScalePos;
#else
    for(j_ = 0; j_ < 3; j_++)
      Gal[galaxy_number_].Pos[j_] *= ScalePos;
#endif
//#endif

#ifdef GUO10
    for(j_ = 0; j_ < 3; j_++)
    {
      tmp_pos_ = wrap(-HaloGal[p].MergCentralPos[j_] + HaloGal[p].Pos[j_],BoxSize);
      tmp_pos_ *=  sqrt(HaloGal[p].MergTime/HaloGal[p].OriMergTime);

      HaloGal[p].Pos[j_]=HaloGal[p].MergCentralPos[j_] + tmp_pos_;

      if(HaloGal[p].Pos[j_] < 0)
        HaloGal[p].Pos[j_] = BoxSize + HaloGal[p].Pos[j_];
      if(HaloGal[p].Pos[j_] > BoxSize)
        HaloGal[p].Pos[j_] = HaloGal[p].Pos[j_] - BoxSize;
    }
#else
    for(j_ = 0; j_ < 3; j_++)
    {
      tmp_pos_ = wrap(-Gal[p].MergCentralPos[j_] + Gal[p].Pos[j_],BoxSize);
#ifdef GUO13
      tmp_pos_ *=  sqrt(Gal[p].MergTime/Gal[p].OriMergTime);
#else
      tmp_pos_ *=  (Gal[p].MergTime/Gal[p].OriMergTime);
#endif
      Gal[p].Pos[j_]=Gal[p].MergCentralPos[j_] + tmp_pos_;

      if(Gal[p].Pos[j_] < 0)
        Gal[p].Pos[j_] = BoxSize + Gal[p].Pos[j_];
      if(Gal[p].Pos[j_] > BoxSize)
        Gal[p].Pos[j_] = Gal[p].Pos[j_] - BoxSize;
    }
#endif

#ifdef GUO10
    //#ifdef SCALE_COSMOLOGY
    //add by Qi. 06/04/2012 to account for the scale of velocity field
    Scale_V_ = scale_v_cen(Halo[HaloGal[central_galaxy_number_].HaloNr].SnapNum);

    for (j_ = 0; j_ < 3 ; j_++)
    {
      dv_ = HaloGal[p].Vel[j_] - HaloGal[central_galaxy_number_].Vel[j_]/Scale_V_;
      dv_ *=sqrt(ScaleMass/ScalePos) * sqrt(AA_OriginalCosm[Halo[HaloGal[central_galaxy_number_].HaloNr].SnapNum]/AA[Halo[HaloGal[central_galaxy_number_].HaloNr].SnapNum]);
      HaloGal[p].Vel[j_] = HaloGal[central_galaxy_number_].Vel[j_] + dv_;
    }
//#endif

#else
    Scale_V_ = scale_v_cen(Halo[Gal[central_galaxy_number_].HaloNr].SnapNum);

    for (j_ = 0; j_ < 3 ; j_++)
    {
      dv_ = Gal[p].Vel[j_] - Gal[central_galaxy_number_].Vel[j_]/Scale_V_;
      dv_ *=sqrt(ScaleMass/ScalePos) * sqrt(AA_OriginalCosm[Halo[Gal[central_galaxy_number_].HaloNr].SnapNum]/AA[Halo[Gal[central_galaxy_number_].HaloNr].SnapNum]);
      Gal[p].Vel[j_] = Gal[central_galaxy_number_].Vel[j_] + dv_;
    }
#endif
  }
}


void get_coordinates(float *pos_, float *vel_, const long long ID_, const int tree_number_, const int halo_number_, const int snapshot_number_)
{
  int m_, k_, start_, n_ids_;

  start_ = OffsetIDs_halo[TreeFirstHalo[tree_number_] + halo_number_] - OffsetIDs;
  n_ids_ = CountIDs_halo[TreeFirstHalo[tree_number_] + halo_number_];

  while(n_ids_ > 0)
  {
    m_ = n_ids_ / 2;
    if(IdList[start_ + m_] == ID_)
    {
      for(k_ = 0; k_ < 3; k_++)
      {
        pos_[k_] = PosList[3 * (start_ + m_) + k_];
        vel_[k_] = sqrt(AA[snapshot_number_]) * VelList[3 * (start_ + m_) + k_];        /* to convert to peculiar velocity */
      }

      if(pos_[0] == 0 && pos_[1] == 0 && pos_[2] == 0)
      {
        terminate("This treeaux-files does not (yet) contain the coordinates\n for the desired output time!\n");
      }

      return;
    }

    if(IdList[start_ + m_] < ID_)
    {
      n_ids_ -= m_;
      start_ += m_;
    }
    else
    {
      n_ids_ = m_;
    }
  }

  terminate("ID_ not found! - What's going on?");
  return;
}


/**@brief If USE_MEMORY_TO_MINIMIZE_IO ON - Routine to read in all
 *        the tree_aux data for one file into a pointer (ptr_auxdata) -
 *        later myfread() will pass the data onto each tree_number_ structure*/

void load_all_auxdata(const int file_number_)
{
  char file_name_[1024];
  FILE *file_;
  struct stat file_status_;
  int SnapShotInFileName_;
  size_t bytes_;

  //if def MCMC and PARALLEL only task 0 reads the representative treefile, then broadcasts
#ifdef PARALLEL
#ifdef MCMC
  if(ThisTask==0)
  {
          printf("Task 0 reading aux data\n");
#endif
#endif
  SnapShotInFileName_=LastDarkMatterSnapShot;

#ifdef MCMC
#ifdef MR_PLUS_MRII
  SnapShotInFileName_=LastDarkMatterSnapShot_MRII;
#endif
#endif

#ifndef MRII
  sprintf(file_name_, "%s/treedata/treeaux_%03d.%d", SimulationDir, SnapShotInFileName_, file_number_);
#else
  sprintf(file_name_, "%s/treedata/treeaux_sf1_%03d.%d", SimulationDir, SnapShotInFileName_, file_number_);
#endif


  if(stat(file_name_, &file_status_) != 0)                  /* seems not to exist */
    {
      char error_message_[2048];
      sprintf(error_message_, "Can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }

  if(!(file_ = fopen(file_name_, "r")))
    {
      char error_message_[2048];
      sprintf(error_message_, "Can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }

  bytes_ = file_status_.st_size;
  TreeAuxData = mymalloc("TreeAuxData", bytes_);

  myfread(TreeAuxData, 1, bytes_, file_);

  fclose(file_);

#ifdef PARALLEL
#ifdef MCMC
  } //end if ThisTask==0

  if(ThisTask==0)
          printf("aux data read  by task %d\n", ThisTask);

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(&bytes_, sizeof(size_t), MPI_BYTE, 0, MPI_COMM_WORLD);
  if(ThisTask>0)
          TreeAuxData = mymalloc("TreeAuxData", bytes_);

  if(ThisTask==0)
          printf("broadcasting aux data\n", ThisTask);

  //MPI_BCast has a limit of 2Gb so everything needs to be passed in smaller chunks
  int ii, Nmessages=10000;
  long long  MsgSizeInBytes=10000000; //chunks of 10MsgSizeInBytes
  for(ii=0;ii<Nmessages;ii++)
    {
      //if next chunk is outside of array size, just pass whats left and then exit the loop
      if((ii+1)*MsgSizeInBytes>bytes_)
        {
          MPI_Bcast(&TreeAuxData[ii*MsgSizeInBytes],bytes_-ii*MsgSizeInBytes, MPI_BYTE, 0, MPI_COMM_WORLD);
          break;
        }
      else
        MPI_Bcast(&TreeAuxData[ii*MsgSizeInBytes],MsgSizeInBytes, MPI_BYTE, 0, MPI_COMM_WORLD);
  }

  if(ThisTask==0)
    printf("done broadcasting aux data\n", ThisTask);

  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  int *header = TreeAuxData;
  NtotHalos = header[0];
  TotIds = header[1];
  Ntrees = header[2];
  TotSnaps = header[3];

  CountIDs_snaptree = header + 4 + 2 * TotSnaps;
  OffsetIDs_snaptree = CountIDs_snaptree +  TotSnaps * Ntrees;
  CountIDs_halo = OffsetIDs_snaptree + TotSnaps * Ntrees;
  OffsetIDs_halo = CountIDs_halo + NtotHalos;
}

#endif

