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

/** @file   save_galtree.c
 *  @date   2016-2019
 *  @author ?
 *  @author Stefan Hilbert
 *
 * @brief   output to galaxy tree files
 **/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

#ifdef OUTPUT_BUFFERING
#include "dynamic_array_inline.h"
#endif /* defined OUTPUT_BUFFERING */

#ifdef GALAXYTREE

#ifdef OUTPUT_BUFFERING

/** initial buffer capacity */
#ifndef OUTPUT_BUFFERING_INIT_CAPACITY
#if OUTPUT_BUFFERING == 1
#define OUTPUT_BUFFERING_INIT_CAPACITY 200 * BYTES_PER_MB
#else /* OUTPUT_BUFFERING == 2 */
#define OUTPUT_BUFFERING_INIT_CAPACITY 1000 * BYTES_PER_MB
#endif /* OUTPUT_BUFFERING == 2 */
#endif /* not defined OUTPUT_BUFFERING_INIT_CAPACITY_IN_MB */


/** @brief the buffer used to buffer galaxy output to disk */
static dynamic_array_type galaxy_output_buffer;


/** @brief init the buffer used to buffer galaxy output to disk */
void save_galaxy_tree_init_output_buffer(void)
{ dynamic_array_init(&galaxy_output_buffer, 0, OUTPUT_BUFFERING_INIT_CAPACITY); }


/** @brief writes output buffer content to file and clears buffer  */
void save_galaxy_tree_flush_output_buffer(void)
{
  myfwrite_large_data(galaxy_output_buffer.data, 1, galaxy_output_buffer.size, FdGalTree);
  dynamic_array_clear(&galaxy_output_buffer);
}


/** @brief shows output buffer stats */
void save_galaxy_tree_show_output_buffer_statistics(void)
{ printf("galaxy tree output buffer: capacity = %lu (%f MB)\n", galaxy_output_buffer.capacity, galaxy_output_buffer.capacity / (1024. * 1024.)); }

#endif /* defined OUTPUT_BUFFERING */


/** @brief create galaxy output files for galaxy trees */
void create_galaxy_tree_file(const int file_number_)
{
  char file_name_[1536];

  sprintf(file_name_, "%s/%s_galtree_%d", OutputDir, FileNameGalaxies, file_number_);
  if(!(FdGalTree = fopen(file_name_, "wb+")))
  {
    char error_message_[2048];
    sprintf(error_message_, "can't open file `%s'\n", file_name_);
    terminate(error_message_);
  }

  /* skip one block to make room for header */
  myfseek(FdGalTree, sizeof(struct GALAXY_OUTPUT), SEEK_SET);
  
#ifdef OUTPUT_BUFFERING
  dynamic_array_clear(&galaxy_output_buffer);
#endif /* defined OUTPUT_BUFFERING */
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
  
#ifdef OUTPUT_BUFFERING
#if OUTPUT_BUFFERING == 2
  /* skip one block to make room for header */
  myfseek(FdGalTree, sizeof(struct GALAXY_OUTPUT), SEEK_SET);
  save_galaxy_tree_flush_output_buffer();
#endif /* OUTPUT_BUFFERING == 2 */
#endif /* defined OUTPUT_BUFFERING */      

  fclose(FdGalTree);
}


/** @brief  Saves the Galaxy_Output structure for galaxies in
 *        the current tree into the current output file (one for each
 *        input dark matter file).
 *
 * for now store SFH bins boht in GALAXY_OUTPUT and SFH_OUTPUT to check things are ok.
 * Later choose only latter mode.
 *
 * @bug (corrected by Stefan Hilbert) sfh_numbins is now set to sfh_ibin + 1 (was sfh_ibin),
 *      this is now done in prepare_galaxy_for_output().
 */
void save_galaxy_tree_append(const int galaxy_number_)
{
  struct GALAXY_OUTPUT galaxy_output_;

  prepare_galaxy_for_output(HaloGal[galaxy_number_].SnapNum, &HaloGal[galaxy_number_], &galaxy_output_);

#ifdef OUTPUT_BUFFERING
  dynamic_array_push_back(&galaxy_output_buffer, &galaxy_output_, sizeof(struct GALAXY_OUTPUT));
#else  /* not defined OUTPUT_BUFFERING */
  myfwrite(&galaxy_output_, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
#endif /* not defined OUTPUT_BUFFERING */  
}


/** @brief all things to be done on output files before next tree arrives 
  *
  * updates galaxy ids (GalID, FirstProgGal, etc.) in galaxy output files
  */
void save_galaxy_tree_finalize(const int file_number_, const int tree_number_)
{
  if(NGalTree > 0)
  {
    // order GalTree by current order of storage in file (IndexStored)
    qsort(GalTree, NGalTree, sizeof(struct galaxy_tree_data), save_galaxy_tree_compare);
    
#ifdef OUTPUT_BUFFERING
    // //debugging:
    // if(galaxy_output_buffer.size < NGalTree * sizeof(struct GALAXY_OUTPUT))
    // {
    //   printf("error: in save_galaxy_tree_finalize(): galaxy_output_buffer.size = %lu < %lu = NGalTree * sizeof(struct GALAXY_OUTPUT).\n"
    //          "NGalTree = %d, sizeof(struct GALAXY_OUTPUT) = %lu\n",
    //          galaxy_output_buffer.size, NGalTree * sizeof(struct GALAXY_OUTPUT), NGalTree, sizeof(struct GALAXY_OUTPUT));
    //   terminate("error: in save_galaxy_tree_finalize(): galaxy_output_buffer.size < NGalTree * sizeof(struct GALAXY_OUTPUT).\n");
    // }

    struct GALAXY_OUTPUT *galaxy_output_ = (struct GALAXY_OUTPUT*) (galaxy_output_buffer.data + galaxy_output_buffer.size - NGalTree * sizeof(struct GALAXY_OUTPUT));

#else  /* not defined OUTPUT_BUFFERING */  
    struct GALAXY_OUTPUT *galaxy_output_ = (struct GALAXY_OUTPUT*) mymalloc("tree_file_gal", sizeof(struct GALAXY_OUTPUT) * NGalTree);
    myfseek(FdGalTree, (1 + TotGalCount) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
    myfread(&(galaxy_output_[0]), sizeof(struct GALAXY_OUTPUT), NGalTree, FdGalTree);
#endif /* not defined OUTPUT_BUFFERING */  
    
    int galaxy_number_;
    for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
      prepare_galaxy_tree_info_for_output(file_number_, tree_number_, &GalTree[galaxy_number_], &galaxy_output_[galaxy_number_]);
    
#ifdef SORT_GALAXY_OUTPUT
    qsort(galaxy_output_, NGalTree, sizeof(struct GALAXY_OUTPUT), output_galaxy_compare);
#endif /* not defined SORT_GALAXY_OUTPUT */  

#ifndef OUTPUT_BUFFERING
    myfseek(FdGalTree, (1 + TotGalCount) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
    myfwrite(&galaxy_output_[0], sizeof(struct GALAXY_OUTPUT), NGalTree, FdGalTree);
    myfree(galaxy_output_);
#endif /* not defined OUTPUT_BUFFERING */
  }
}


static int GalCount; /* used to assign galaxy ids via walking galaxy tree */

/** @brief walks galaxy tree to compute foreign gal ID keys */
int walk_galaxy_tree(const int galaxy_number_)
{
  int last_galaxy_number_ = galaxy_number_;

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


/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) in galaxy tree_number_  */
void update_galaxy_tree_ids(void)
{
  int galaxy_number_, progenitor_galaxy_number_, snapshot_number_;

  for(galaxy_number_ = 0; galaxy_number_ < NGalTree; galaxy_number_++)
  {
    GalTree[galaxy_number_].Done         = 0;
    GalTree[galaxy_number_].LastProgGal = -1;
    GalTree[galaxy_number_].MainLeaf    = -1;
    GalTree[galaxy_number_].TreeRoot    = -1;
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


/** @brief comparison function for sorting tree_number_ data for updating tree_number_ info for galaxy on disk */
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
