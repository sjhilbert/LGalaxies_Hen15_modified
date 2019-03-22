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
void create_galaxy_tree_file(const int filenr)
{
  char buf[1000];

  sprintf(buf, "%s/%s_galtree_%d", OutputDir, FileNameGalaxies, filenr);
  if(!(FdGalTree = fopen(buf, "wb+")))
  {
    char sbuf[1000];
    sprintf(sbuf, "can't open file `%s'\n", buf);
    terminate(sbuf);
  }

  /* skip one block to make room for header */
  myfseek(FdGalTree, sizeof(struct GALAXY_OUTPUT), SEEK_SET);
  TotGalCount = 0;
}


void close_galaxy_tree_file(void)
{
  int one = 1;
  int size_of_struct = sizeof(struct GALAXY_OUTPUT);

  /* write header information  */
  myfseek(FdGalTree, 0, SEEK_SET);
  myfwrite(&one, sizeof(int), 1, FdGalTree);        // write 1
  myfwrite(&size_of_struct, sizeof(int), 1, FdGalTree);        // size of an output structure (Galaxy_Output)
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
void save_galaxy_tree_append(int i)
{
  struct GALAXY_OUTPUT galaxy_output;

  prepare_galaxy_for_output(HaloGal[i].SnapNum, &HaloGal[i], &galaxy_output);

#ifdef STAR_FORMATION_HISTORY
  galaxy_output.sfh_numbins = galaxy_output.sfh_ibin;
#endif

  myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
}


/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) in galaxy tree  */
void update_galaxy_tree_ids(void)
{
  int i, p, num;

  for(i = 0; i < NGalTree; i++)
  {
    GalTree[i].Done = 0;
    GalTree[i].LastProgGal = -1;
    GalTree[i].MainLeaf = -1;
    GalTree[i].TreeRoot = -1;
  }

  GalCount = 0;

  for(num = LastDarkMatterSnapShot; num >= 0; num--)
  {
    for(i = 0; i < NGalTree; i++)
    {
      if(GalTree[i].SnapNum == num)
        if(GalTree[i].Done == 0)
            walk_galaxy_tree(i);
    }
  }

  for(i = 0; i < NGalTree; i++)
  {
    p = GalTree[i].FirstProgGal;
    while(p >= 0)
    {
      GalTree[p].DescendantGal = i;
      p = GalTree[p].NextProgGal;
    }
  }

  for(i = 0; i < NGalTree; i++)
  {
    if(GalTree[i].FirstProgGal >= 0)
      GalTree[i].FirstProgGal = GalTree[GalTree[i].FirstProgGal].GalID;

    if(GalTree[i].LastProgGal >= 0)
      GalTree[i].LastProgGal = GalTree[GalTree[i].LastProgGal].GalID;

    if(GalTree[i].MainLeaf >= 0)
      GalTree[i].MainLeaf = GalTree[GalTree[i].MainLeaf].GalID;

    if(GalTree[i].TreeRoot >= 0)
      GalTree[i].TreeRoot = GalTree[GalTree[i].TreeRoot].GalID;

    if(GalTree[i].NextProgGal >= 0)
      GalTree[i].NextProgGal = GalTree[GalTree[i].NextProgGal].GalID;

    if(GalTree[i].DescendantGal >= 0)
      GalTree[i].DescendantGal = GalTree[GalTree[i].DescendantGal].GalID;

    if(GalTree[i].FOFCentralGal >= 0)
      GalTree[i].FOFCentralGal = GalTree[GalTree[i].FOFCentralGal].GalID;
  }
}


/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) in galaxy output files  */
void save_galaxy_tree_finalize(const int filenr, const int tree)
{
  int i;
  struct GALAXY_OUTPUT galaxy_output;
  
  // order GalTree by current order of storage in file (IndexStored)
  qsort(GalTree, NGalTree, sizeof(struct galaxy_tree_data), save_galaxy_tree_compare);

  /* Before, the header was a simple integer for number of galaxies. So, the
     code had to jump over an int (used to store the number of galaxies) and
     and over all the galaxies written so far */ 
  // for DB compatible output, pad the first line with the size of one struct.
 
  for(i = 0; i < NGalTree; i++)
  { 
    myfseek(FdGalTree, (1 + TotGalCount + i) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
    myfread(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);

    prepare_galaxy_tree_info_for_output(filenr, tree, &GalTree[i], &galaxy_output);

    myfseek(FdGalTree, (1 + TotGalCount + i) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
    myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
  }

  // GL: propose to only reorder file on disk based on input (or Makefile) parameter
  // if DB is fast in ordering, eg using SSDs, we could leave it out here.
  // BTW it is BETTER not to reorder galaxies if SFHBins are not also reordered,
  // if at least both are to be used in light cone post processing.
  // for easier to read
#ifdef SORT_GALAXYTREE_OUTPUT
  save_galaxy_tree_reorder_on_disk();
#else  /* not defined SORT_GALAXYTREE_OUTPUT */
  myfseek(FdGalTree, (1 + TotGalCount + NGalTree) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
#endif /* not defined SORT_GALAXYTREE_OUTPUT */

  TotGalCount += NGalTree;
}


/** @brief updates galaxy ids (GalID, FirstProgGal, etc.) for a galaxy in memory */
void prepare_galaxy_tree_info_for_output(const int filenr, const int tree, const struct galaxy_tree_data *g, struct GALAXY_OUTPUT *o)
{
  long long big = calc_big_db_offset(filenr, tree);

  o->GalID = g->GalID;
  o->FOFCentralGal = g->FOFCentralGal; 
  o->FirstProgGal = g->FirstProgGal;
  o->NextProgGal = g->NextProgGal;
  o->LastProgGal = g->LastProgGal;
  o->MainLeafId = g->MainLeaf;
  o->TreeRootId = g->TreeRoot;
  o->DescendantGal = g->DescendantGal;
  o->FileTreeNr = big;

#ifdef CONTINUOUS_TREES
  // Reset big (so only FileTreeNr has original value)
  // Then new values should coincide with positions in the file
  big = TotGalCount;
#endif

  o->GalID += big;
  o->FOFCentralGal += big;

  if(o->FirstProgGal >= 0)
    o->FirstProgGal += big;

  if(o->LastProgGal >= 0)
    o->LastProgGal += big;
  else
    o->LastProgGal = o->GalID;

  if(o->MainLeafId >= 0)
    o->MainLeafId += big;
  else
    o->MainLeafId = o->GalID;

  if(o->TreeRootId >= 0)
    o->TreeRootId += big;
  else
  {
    terminate("o->TreeRootId < 0");
    o->TreeRootId = -1;
  }

  if(o->NextProgGal >= 0)
    o->NextProgGal += big;

  if(o->DescendantGal >= 0)
    o->DescendantGal += big;
}


/** @brief walks galaxy tree to compute foreign gal ID keys */
int walk_galaxy_tree(const int nr)
{
  int last;

  last = nr;

  if(GalTree[nr].Done == 0)
  {
    GalTree[nr].Done = 1;
    GalTree[nr].GalID = GalCount++;

    if(GalTree[nr].TreeRoot == -1)
      GalTree[nr].TreeRoot = nr;

    if(GalTree[nr].FirstProgGal >= 0)
    {
      GalTree[GalTree[nr].FirstProgGal].TreeRoot = GalTree[nr].TreeRoot;
      last = walk_galaxy_tree(GalTree[nr].FirstProgGal);
      GalTree[nr].MainLeaf = GalTree[GalTree[nr].FirstProgGal].MainLeaf;
    }
    else
      GalTree[nr].MainLeaf = nr;

    GalTree[nr].LastProgGal = last;

    if(GalTree[nr].NextProgGal >= 0)
    {
      if(GalTree[nr].NextProgGal >= NGalTree)
      {
        printf("\n nr=%d NGalTree=%d GalTree[nr].NextProgGal=%d\n", 
                nr, NGalTree, GalTree[nr].NextProgGal);
        terminate("GalTree[nr].NextProgGal >= NGalTree");
      }

      GalTree[GalTree[nr].NextProgGal].TreeRoot = GalTree[nr].TreeRoot;
      last = walk_galaxy_tree(GalTree[nr].NextProgGal);
    }
  }

  return last;
}


struct mp_tree_data
{
  int index;
  long long key;
};


/** @brief sorts galaxies in output files  */
void save_galaxy_tree_reorder_on_disk(void)
{
  int i, idsource, idsave, dest;
  struct GALAXY_OUTPUT galaxy_save, galaxy_source;
  
  struct mp_tree_data *mp = (struct mp_tree_data *) mymalloc("mp", sizeof(struct mp_tree_data) * NGalTree);
  int *id = (int *) mymalloc("id", sizeof(int) * NGalTree);
  
  for(i = 0; i < NGalTree; i++)
  {
    mp[i].index       = i;
    mp[i].key = GalTree[i].GalID;
  }

  
  qsort(mp, NGalTree, sizeof(struct mp_tree_data), save_galaxy_tree_mp_comp);
  
  for(i = 0; i < NGalTree; i++)
    id[mp[i].index] = i;

  for(i = 0; i < NGalTree; i++)
  {
    if(id[i] != i)
    {
      myfseek(FdGalTree, (1 + TotGalCount + i) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
      myfread(&galaxy_source, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
      idsource = id[i];
      dest = id[i];

      do
      {
        myfseek(FdGalTree, (1 + TotGalCount + dest) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
        myfread(&galaxy_save, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
        idsave = id[dest];

        myfseek(FdGalTree, (1 + TotGalCount + dest) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
        myfwrite(&galaxy_source, sizeof(struct GALAXY_OUTPUT), 1, FdGalTree);
        id[dest] = idsource;

        if(dest == i)
          break;

        galaxy_source = galaxy_save;
        idsource = idsave;

        dest = idsource;
      }
      while(1);
    }
  }

  myfree(id);
  myfree(mp);

  myfseek(FdGalTree, (1 + TotGalCount + NGalTree) * sizeof(struct GALAXY_OUTPUT), SEEK_SET);
}


/** @brief comparison function for sorting tree data for updating tree info for galaxie on disk */
int save_galaxy_tree_compare(const void *a, const void *b)
{
  if(((struct galaxy_tree_data *) a)->IndexStored < ((struct galaxy_tree_data *) b)->IndexStored)
    return -1;

  else if(((struct galaxy_tree_data *) a)->IndexStored > ((struct galaxy_tree_data *) b)->IndexStored)
    return +1;

  else
    return 0;
}


/** @brief comparison function for sorting galaxy data on disk */
int save_galaxy_tree_mp_comp(const void *a, const void *b)
{
  if(((struct mp_tree_data *) a)->key < ((struct mp_tree_data *) b)->key)
    return -1;

  else if(((struct mp_tree_data *) a)->key > ((struct mp_tree_data *) b)->key)
    return +1;

  else
    return 0;
}


#endif
