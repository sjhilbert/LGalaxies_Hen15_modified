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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

#ifdef GALAXYTREE

/** @brief create galaxy output files for galaxy trees */
void create_galaxy_tree_file(const int file_number_)
{
  char file_name_[1000];

  sprintf(file_name_, "%s/%s_galtree_%d", OutputDir, FileNameGalaxies, file_number_);
  if(!(FdGalTree = fopen(file_name_, "wb+")))
  {
    char error_message_[1000];
    sprintf(error_message_, "can't open file `%s'\n", file_name_);
    terminate(error_message_);
  }

  /* skip one block to make room for header */
  myfseek(FdGalTree, sizeof(struct GALAXY_OUTPUT), SEEK_SET);
  TotGalCount = 0;
}


void close_galaxy_tree_file(void)
{
  int one_ = 1;
  int size_of_struct_ = sizeof(struct GALAXY_OUTPUT);

  /* write header information  */
  myfseek(FdGalTree, 0, SEEK_SET);
  myfwrite(&one_, sizeof(int), 1, FdGalTree);        // write 1
  myfwrite(&size_of_struct_, sizeof(int), 1, FdGalTree);        // size of an output structure (Galaxy_Output)
  myfwrite(&TotGalCount, sizeof(int), 1, FdGalTree);        // the total number of galaxies

  fclose(FdGalTree);
}


/** @brief  Saves the Galaxy_Output structure for galaxies in
 *        the current tree into the current output file (one for each
 *        input dark matter file).
 *
 * for now store SFH bins boht in GALAXY_OUTPUT and SFH_OUTPUT to check things are ok.
 * Later choose only latter mode.
 */
void save_galaxy_tree_append(const int galaxy_number_)
{
  struct GALAXY_OUTPUT galaxy_output_;

  prepare_galaxy_for_output(HaloGal[galaxy_number_].SnapNum, &HaloGal[galaxy_number_], &galaxy_output_);

#ifdef STAR_FORMATION_HISTORY
  galaxy_output_.sfh_numbins = galaxy_output_.sfh_ibin;
#endif

  myfwrite(&galaxy_output_, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
}


/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) in galaxy tree_number_  */
void update_galaxy_tree_ids(void)
{
  int galaxy_number_, progenitor_galaxy_number_, snapshot_number_;

  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
  {
    GalTree[galaxy_number_].Done = 0;
    GalTree[galaxy_number_].LastProgGal = -1;
    GalTree[galaxy_number_].MainLeaf = -1;
    GalTree[galaxy_number_].TreeRoot = -1;
  }

  GalCount = 0;

  for(snapshot_number_ = LastDarkMatterSnapShot; snapshot_number_ >= 0; snapshot_number_--)
  {
    for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
    {
      if(GalTree[galaxy_number_].SnapNum == snapshot_number_)
        if(GalTree[galaxy_number_].Done == 0)
            walk_galaxy_tree(galaxy_number_);
    }
  }

  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
  {
    progenitor_galaxy_number_ = GalTree[galaxy_number_].FirstProgGal;
    while(progenitor_galaxy_number_ >= 0)
    {
      GalTree[progenitor_galaxy_number_].DescendantGal = galaxy_number_;
      progenitor_galaxy_number_ = GalTree[progenitor_galaxy_number_].NextProgGal;
    }
  }

  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
  {
    if(GalTree[galaxy_number_].FirstProgGal >= 0)
      GalTree[galaxy_number_].FirstProgGal = GalTree[GalTree[galaxy_number_].FirstProgGal].GalID;

    if(GalTree[galaxy_number_].LastProgGal >= 0)
      GalTree[galaxy_number_].LastProgGal = GalTree[GalTree[galaxy_number_].LastProgGal].GalID;

    if(GalTree[galaxy_number_].MainLeaf >= 0)
      GalTree[galaxy_number_].MainLeaf = GalTree[GalTree[galaxy_number_].MainLeaf].GalID;

    if(GalTree[galaxy_number_].TreeRoot >= 0)
      GalTree[galaxy_number_].TreeRoot = GalTree[GalTree[galaxy_number_].TreeRoot].GalID;

    if(GalTree[galaxy_number_].NextProgGal >= 0)
      GalTree[galaxy_number_].NextProgGal = GalTree[GalTree[galaxy_number_].NextProgGal].GalID;

    if(GalTree[galaxy_number_].DescendantGal >= 0)
      GalTree[galaxy_number_].DescendantGal = GalTree[GalTree[galaxy_number_].DescendantGal].GalID;

    if(GalTree[galaxy_number_].FOFCentralGal >= 0)
      GalTree[galaxy_number_].FOFCentralGal = GalTree[GalTree[galaxy_number_].FOFCentralGal].GalID;
  }
}


/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) in galaxy output files  */
void save_galaxy_tree_finalize(const int file_number_, const int tree_number_)
{
  
  // order GalTree by current order of storage in file (IndexStored)
  qsort(GalTree, NGalTree, sizeof(struct galaxy_tree_data), save_galaxy_tree_compare);

  /* Before, the header was a simple integer for number of galaxies. So, the
     code had to jump over an int (used to store the number of galaxies) and
     and over all the galaxies written so far */ 
  // for DB compatible output, pad the first line with the size of one struct.
  
#ifdef UPDATE_GALAXY_OUTPUT_IN_MEM

  struct GALAXY_OUTPUT *galaxy_output_ = (struct GALAXY_OUTPUT*) mymalloc("tree_file_gal", sizeof(struct GALAXY_OUTPUT) * NGalTree);
    
  myfseek(FdGalTree, (1 + TotGalCount) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
  myfread(&(galaxy_output_[0]), sizeof(struct GALAXY_OUTPUT), NGalTree, FdGalTree);

  int galaxy_number_;
  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
    prepare_galaxy_tree_info_for_output(file_number_, tree_number_, &GalTree[galaxy_number_], &galaxy_output_[galaxy_number_]);

  
#ifdef SORT_GALAXY_OUTPUT
  qsort(galaxy_output_, NGalTree, sizeof(struct GALAXY_OUTPUT), output_galaxy_compare);
#endif /* not defined SORT_GALAXY_OUTPUT */  

  myfseek(FdGalTree, (1 + TotGalCount) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
  myfwrite(&galaxy_output_[0], sizeof(struct GALAXY_OUTPUT), NGalTree, FdGalTree);
  myfree(galaxy_output_);

#else /* not defined UPDATE_GALAXY_OUTPUT_IN_MEM */
  
  struct GALAXY_OUTPUT galaxy_output_;
  
  int galaxy_number_;
  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
  { 
    myfseek(FdGalTree, (1 + TotGalCount + galaxy_number_) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
    myfread(&galaxy_output_, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);

    prepare_galaxy_tree_info_for_output(file_number_, tree_number_, &GalTree[galaxy_number_], &galaxy_output_);

    myfseek(FdGalTree, (1 + TotGalCount + galaxy_number_) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
    myfwrite(&galaxy_output_, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
  }
  
#ifdef SORT_GALAXY_OUTPUT
  save_galaxy_tree_reorder_on_disk();
#endif /* defined SORT_GALAXY_OUTPUT */
#endif /* not defined UPDATE_GALAXY_OUTPUT_IN_MEM */

  TotGalCount += NGalTree;
}


/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) for a galaxy in memory */
void prepare_galaxy_tree_info_for_output(const int file_number_, const int tree_number_, const struct galaxy_tree_data *galaxy_, struct GALAXY_OUTPUT *output_galaxy_)
{
  const long long big_db_offset_ = calc_big_db_offset(file_number_, tree_number_);

  output_galaxy_->GalID = galaxy_->GalID;
  output_galaxy_->FOFCentralGal = galaxy_->FOFCentralGal; 
  output_galaxy_->FirstProgGal = galaxy_->FirstProgGal;
  output_galaxy_->NextProgGal = galaxy_->NextProgGal;
  output_galaxy_->LastProgGal = galaxy_->LastProgGal;
  output_galaxy_->MainLeafId = galaxy_->MainLeaf;
  output_galaxy_->TreeRootId = galaxy_->TreeRoot;
  output_galaxy_->DescendantGal = galaxy_->DescendantGal;
  output_galaxy_->FileTreeNr = big_db_offset_;

#ifdef CONTINUOUS_TREES
  // Reset big_offset_ (so only FileTreeNr has original value)
  // Then new values should coincide with positions in the file
  const long long big_offset_ = TotGalCount;
#else  /* not defined CONTINUOUS_TREES */
  const long long big_offset_ = big_db_offset_;
#endif /* not defined CONTINUOUS_TREES */

  output_galaxy_->GalID += big_offset_;
  output_galaxy_->FOFCentralGal += big_offset_;

  if(output_galaxy_->FirstProgGal >= 0)
    output_galaxy_->FirstProgGal += big_offset_;

  if(output_galaxy_->LastProgGal >= 0)
    output_galaxy_->LastProgGal += big_offset_;
  else
    output_galaxy_->LastProgGal = output_galaxy_->GalID;

  if(output_galaxy_->MainLeafId >= 0)
    output_galaxy_->MainLeafId += big_offset_;
  else
    output_galaxy_->MainLeafId = output_galaxy_->GalID;

  if(output_galaxy_->TreeRootId >= 0)
    output_galaxy_->TreeRootId += big_offset_;
  else
  {
    terminate("output_galaxy_->TreeRootId < 0");
    output_galaxy_->TreeRootId = -1;
  }

  if(output_galaxy_->NextProgGal >= 0)
    output_galaxy_->NextProgGal += big_offset_;

  if(output_galaxy_->DescendantGal >= 0)
    output_galaxy_->DescendantGal += big_offset_;
}


/** @brief walks galaxy tree_number_ to compute foreign gal ID keys */
int walk_galaxy_tree(const int galaxy_number_)
{
  int last_galaxy_number_;

  last_galaxy_number_ = galaxy_number_;

  if(GalTree[galaxy_number_].Done == 0)
  {
    GalTree[galaxy_number_].Done = 1;
    GalTree[galaxy_number_].GalID = GalCount++;

    if(GalTree[galaxy_number_].TreeRoot == -1)
      GalTree[galaxy_number_].TreeRoot = galaxy_number_;

    if(GalTree[galaxy_number_].FirstProgGal >= 0)
    {
      GalTree[GalTree[galaxy_number_].FirstProgGal].TreeRoot = GalTree[galaxy_number_].TreeRoot;
      last_galaxy_number_ = walk_galaxy_tree(GalTree[galaxy_number_].FirstProgGal);
      GalTree[galaxy_number_].MainLeaf = GalTree[GalTree[galaxy_number_].FirstProgGal].MainLeaf;
    }
    else
      GalTree[galaxy_number_].MainLeaf = galaxy_number_;

    GalTree[galaxy_number_].LastProgGal = last_galaxy_number_;

    if(GalTree[galaxy_number_].NextProgGal >= 0)
    {
      if(GalTree[galaxy_number_].NextProgGal >= NGalTree)
      {
        printf("\n galaxy_number_=%d NGalTree=%d GalTree[galaxy_number_].NextProgGal=%d\n", 
                galaxy_number_, NGalTree, GalTree[galaxy_number_].NextProgGal);
        terminate("GalTree[galaxy_number_].NextProgGal >= NGalTree");
      }

      GalTree[GalTree[galaxy_number_].NextProgGal].TreeRoot = GalTree[galaxy_number_].TreeRoot;
      last_galaxy_number_ = walk_galaxy_tree(GalTree[galaxy_number_].NextProgGal);
    }
  }

  return last_galaxy_number_;
}


struct mp_tree_data
{
  int index;
  long long key;
};


/** @brief sorts galaxies in output files  */
void save_galaxy_tree_reorder_on_disk(void)
{
  int galaxy_number_, id_source_, id_save_, dest_;
  struct GALAXY_OUTPUT galaxy_save_, galaxy_source_;
  
  struct mp_tree_data *mp_ = (struct mp_tree_data *) mymalloc("mp_", sizeof(struct mp_tree_data) * NGalTree);
  int *id_ = (int *) mymalloc("id_", sizeof(int) * NGalTree);
  
  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
  {
    mp_[galaxy_number_].index =         galaxy_number_;
    mp_[galaxy_number_].key   = GalTree[galaxy_number_].GalID;
  }
  
  qsort(mp_, NGalTree, sizeof(struct mp_tree_data), save_galaxy_tree_mp_comp);
  
  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
    id_[mp_[galaxy_number_].index] = galaxy_number_;

  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
  {
    if(id_[galaxy_number_] != galaxy_number_)
    {
      myfseek(FdGalTree, (1 + TotGalCount + galaxy_number_) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
      myfread(&galaxy_source_, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
      id_source_ = id_[galaxy_number_];
      dest_ = id_[galaxy_number_];

      do
      {
        myfseek(FdGalTree, (1 + TotGalCount + dest_) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
        myfread(&galaxy_save_, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
        id_save_ = id_[dest_];

        myfseek(FdGalTree, (1 + TotGalCount + dest_) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
        myfwrite(&galaxy_source_, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
        id_[dest_] = id_source_;

        if(dest_ == galaxy_number_)
          break;

        galaxy_source_ = galaxy_save_;
        id_source_ = id_save_;

        dest_ = id_source_;
      }
      while(1);
    }
  }

  myfree(id_);
  myfree(mp_);

  myfseek(FdGalTree, (1 + TotGalCount + NGalTree) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
}


/** @brief compares lightcone galaxy entries for sorting
 * 
 *  compares GALAXY_OUTPUT entries for sorting output on disk
 */
int output_galaxy_compare(const void *output_galaxy_a_, const void *output_galaxy_b_)
{
  if(((struct GALAXY_OUTPUT*) output_galaxy_a_)->GalID < ((struct GALAXY_OUTPUT*) output_galaxy_b_)->GalID)
    return -1;

  else if(((struct GALAXY_OUTPUT*) output_galaxy_a_)->GalID > ((struct GALAXY_OUTPUT*) output_galaxy_b_)->GalID)
    return +1;

  else 
    return 0;
}


/** @brief comparison function for sorting tree_number_ data for updating tree_number_ info for galaxie on disk */
int save_galaxy_tree_compare(const void *a_, const void *b_)
{
  if(((struct galaxy_tree_data *) a_)->IndexStored < ((struct galaxy_tree_data *) b_)->IndexStored)
    return -1;

  else if(((struct galaxy_tree_data *) a_)->IndexStored > ((struct galaxy_tree_data *) b_)->IndexStored)
    return +1;

  else
    return 0;
}


/** @brief comparison function for sorting galaxy data on disk */
int save_galaxy_tree_mp_comp(const void *mp_tree_data_a_, const void *mp_tree_data_b_)
{
  if(((struct mp_tree_data *) mp_tree_data_a_)->key < ((struct mp_tree_data *) mp_tree_data_b_)->key)
    return -1;

  else if(((struct mp_tree_data *) mp_tree_data_a_)->key > ((struct mp_tree_data *) mp_tree_data_b_)->key)
    return +1;

  else
    return 0;
}


#endif /* defined GALAXYTREE */
