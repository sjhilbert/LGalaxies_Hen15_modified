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

#ifdef PARALLEL
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


#define MAXBLOCKS 5000
#define MAXCHARS  32

static size_t TotBytes;
static void *Base;

static unsigned long Nblocks;

static void **Table;
static size_t *BlockSize;
static char *MovableFlag;
static void ***BasePointers;

static char *VarName;
static char *FunctionName;
static char *FileName;
static int *LineNumber;


void mymalloc_init(void)
{
  size_t n_;

  BlockSize = (size_t *) malloc(MAXBLOCKS * sizeof(size_t));
  Table = (void **) malloc(MAXBLOCKS * sizeof(void *));
  MovableFlag = (char *) malloc(MAXBLOCKS * sizeof(char));
  BasePointers = (void ***) malloc(MAXBLOCKS * sizeof(void **));
  VarName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FunctionName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FileName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  LineNumber = (int *) malloc(MAXBLOCKS * sizeof(int));

  memset(VarName     , 0, MAXBLOCKS * MAXCHARS);
  memset(FunctionName, 0, MAXBLOCKS * MAXCHARS);
  memset(FileName    , 0, MAXBLOCKS * MAXCHARS);

  n_ = MaxMemSize * ((size_t) 1024 * 1024);

  if(!(Base = malloc(n_)))
  {
    printf("Failed to allocate memory for `Base' (%g Mbytes).\n", MaxMemSize);
    terminate("failure to allocate memory");
  }

  TotBytes = FreeBytes = n_;

  AllocatedBytes = 0;
  Nblocks = 0;
  HighMarkBytes = 0;
}


void report_detailed_memory_usage_of_largest_task(size_t * OldHighMarkBytes_, const char *label_, const char *func_, const char *file_, const int line_)
{
  if(AllocatedBytes > 1.1 * (*OldHighMarkBytes_))
  {
    *OldHighMarkBytes_ = AllocatedBytes;

//#ifndef MCMC
    printf("\nAt '%s', %s()/%s/%d: Allocation = %g Mbyte (on task=%d)\n\n",
         label_, func_, file_, line_, AllocatedBytes / (1024.0 * 1024.0), ThisTask);
  dump_memory_table();
//#endif
  
    fflush(stdout);
  }
}



void dump_memory_table(void)
{
  unsigned int block_number_;
  size_t totBlocksize_ = 0;

  printf("------------------------ Allocated Memory Blocks-------------------------------------------------------\n");
  printf("Task   Nr F                          Variable      MBytes   Cumulative         Function/File/Linenumber\n");
  printf("-------------------------------------------------------------------------------------------------------\n");
  for(block_number_ = 0; block_number_ < Nblocks; block_number_++)
  {
    totBlocksize_ += BlockSize[block_number_];

    printf("%4d %4d %d  %32s  %10.4f   %10.4f  %s()/%s/%d\n",
     ThisTask, block_number_, MovableFlag[block_number_], VarName + block_number_ * MAXCHARS, BlockSize[block_number_] / (1024.0 * 1024.0),
     totBlocksize_ / (1024.0 * 1024.0), FunctionName + block_number_ * MAXCHARS,
     FileName + block_number_ * MAXCHARS, LineNumber[block_number_]);
  }
  printf("-------------------------------------------------------------------------------------------------------\n");
}


void *mymalloc_fullinfo(const char *var_name_, size_t n_, const char *func_, const char *file_, const int line_)
{
  if((n_ % 8) > 0)
    n_ = (n_ / 8 + 1) * 8;

  if(n_ < 8)
    n_ = 8;

  if(Nblocks >= MAXBLOCKS)
  {
    char error_message_[1000];
    sprintf(error_message_, "Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line_ %d. MAXBLOCKS=%d\n",
      ThisTask, func_, file_, line_, MAXBLOCKS);
    terminate(error_message_);
  }

  if(n_ > FreeBytes)
  {
    dump_memory_table();
    char error_message_[1000];
    sprintf(error_message_,"\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line_ %d (FreeBytes=%g MB).\n",
           ThisTask, n_ / (1024.0 * 1024.0), var_name_, func_, file_, line_, FreeBytes / (1024.0 * 1024.0));
    terminate(error_message_);
  }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n_;

  strncpy(VarName + Nblocks * MAXCHARS, var_name_, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func_, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file_, MAXCHARS - 1);
  LineNumber[Nblocks] = line_;

  AllocatedBytes += n_;
  BlockSize[Nblocks] = n_;
  MovableFlag[Nblocks] = 0;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}


void *mymalloc_movable_fullinfo(void *ptr_, const char *var_name_, size_t n_, const char *func_, const char *file_, const int line_)
{
  if((n_ % 8) > 0)
    n_ = (n_ / 8 + 1) * 8;

  if(n_ < 8)
    n_ = 8;

  if(Nblocks >= MAXBLOCKS)
  {
    char error_message_[1000];
    sprintf(error_message_, "Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line_ %d. MAXBLOCKS=%d\n",
      ThisTask, func_, file_, line_, MAXBLOCKS);
    terminate(error_message_);
  }

  if(n_ > FreeBytes)
  {
    dump_memory_table();
    char error_message_[1000];
    sprintf(error_message_,"\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line_ %d (FreeBytes=%g MB).\n",
                ThisTask, n_ / (1024.0 * 1024.0), var_name_, func_, file_, line_, FreeBytes / (1024.0 * 1024.0));
    terminate(error_message_);
  }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n_;

  strncpy(VarName + Nblocks * MAXCHARS, var_name_, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func_, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file_, MAXCHARS - 1);
  LineNumber[Nblocks] = line_;

  AllocatedBytes += n_;
  BlockSize[Nblocks] = n_;
  MovableFlag[Nblocks] = 1;
  BasePointers[Nblocks] = ptr_;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}


void myfree_fullinfo(void *ptr_, const char *func_, const char *file_, const int line_)
{
  if(Nblocks == 0)
    terminate("Nblocks == 0");

  if(ptr_ != Table[Nblocks - 1])
  {
    dump_memory_table();
    char error_message_[1000];
    sprintf(error_message_, "Task=%d: Wrong call of myfree() at %s()/%s/line_ %d: not the last allocated block!\n",
      ThisTask, func_, file_, line_);
    terminate(error_message_);
  }

  Nblocks -= 1;
  AllocatedBytes -= BlockSize[Nblocks];
  FreeBytes += BlockSize[Nblocks];
}


void myfree_movable_fullinfo(void *ptr_, const char *func_, const char *file_, const int line_)
{
  unsigned int i, nr;

  if(Nblocks == 0)
    terminate("no allocated blocks that could be freed");

  /* first, let's find the block */
  bool found = false;
  for(nr = Nblocks; nr--;)
    if(ptr_ == Table[nr])
    {
      found = true;
      break;
    }
    
  if(!found)
  {
    dump_memory_table();
    char error_message_[1000];
    sprintf(error_message_,
         "Task=%d: Wrong call of myfree_movable() from %s()/%s/line_ %d - this block has not been allocated!\n",
         ThisTask, func_, file_, line_);
    terminate(error_message_);
  }

  if(nr < Nblocks - 1)                /* the block is not the last allocated block */
  {
    /* check that all subsequent blocks are actually movable */
    for(i = nr + 1; i < Nblocks; i++)
      if(MovableFlag[i] == 0)
      {
        dump_memory_table();
        char error_message_[1000];
        sprintf(error_message_,"Task=%d: Wrong call of myfree_movable() from %s()/%s/line_ %d - behind block=%d there are subsequent non-movable allocated blocks\n",
                ThisTask, func_, file_, line_, nr);
        terminate(error_message_);
      }
  }

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  size_t offset = -BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
  {
    Table[i] += offset;
    *BasePointers[i] = *BasePointers[i] + offset;
  }

  for(i = nr + 1; i < Nblocks; i++)
  {
    Table[i - 1] = Table[i];
    BasePointers[i - 1] = BasePointers[i];
    BlockSize[i - 1] = BlockSize[i];
    MovableFlag[i - 1] = MovableFlag[i];

    strncpy(VarName + (i - 1) * MAXCHARS, VarName + i * MAXCHARS, MAXCHARS - 1);
    strncpy(FunctionName + (i - 1) * MAXCHARS, FunctionName + i * MAXCHARS, MAXCHARS - 1);
    strncpy(FileName + (i - 1) * MAXCHARS, FileName + i * MAXCHARS, MAXCHARS - 1);
    LineNumber[i - 1] = LineNumber[i];
  }

  Nblocks -= 1;
}


void *myrealloc_fullinfo(void *ptr_, size_t n_, const char *func_, const char *file_, const int line_)
{
  if((n_ % 8) > 0)
    n_ = (n_ / 8 + 1) * 8;

  if(n_ < 8)
    n_ = 8;

  if(Nblocks == 0)
    terminate("no allocated blocks that could be reallocated");

  if(ptr_ != Table[Nblocks - 1])
  {
    char error_message_[1000];
    dump_memory_table();
    sprintf(error_message_, "Task=%d: Wrong call of myrealloc() at %s()/%s/line_ %d - not the last allocated block!\n",
            ThisTask, func_, file_, line_);
    terminate(error_message_);
  }

  AllocatedBytes -= BlockSize[Nblocks - 1];
  FreeBytes += BlockSize[Nblocks - 1];

  if(n_ > FreeBytes)
  {
    dump_memory_table();
    char error_message_[1000];
    sprintf(error_message_, "Task=%d: Not enough memory in myremalloc(n_=%g MB) at %s()/%s/line_ %d. previous=%g FreeBytes=%g MB\n",
            ThisTask, n_ / (1024.0 * 1024.0), func_, file_, line_, BlockSize[Nblocks - 1] / (1024.0 * 1024.0),
            FreeBytes / (1024.0 * 1024.0));
    terminate(error_message_);
  }
  Table[Nblocks - 1] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n_;

  AllocatedBytes += n_;
  BlockSize[Nblocks - 1] = n_;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}


void *myrealloc_movable_fullinfo(void *ptr_, size_t n_, const char *func_, const char *file_, const int line_)
{
  unsigned int i, nr;

  if((n_ % 8) > 0)
    n_ = (n_ / 8 + 1) * 8;

  if(n_ < 8)
    n_ = 8;

  if(Nblocks == 0)
    terminate("no allocated blocks that could be reallocated");

  /* first, let's find the block */
  bool found = false;
  for(nr = Nblocks; nr--; )
    if(ptr_ == Table[nr])
    {
      found = true;
      break;
    }
 
  if(!found)
  {
    dump_memory_table();
    char error_message_[1000];
    sprintf(error_message_, "Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line_ %d - this block has not been allocated!\n",
            ThisTask, func_, file_, line_);
    terminate(error_message_);
  }
 
  if(nr < Nblocks - 1)                /* the block is not the last allocated block */
  {
    /* check that all subsequent blocks are actually movable */
    for(i = nr + 1; i < Nblocks; i++)
      if(MovableFlag[i] == 0)
        {
        dump_memory_table();
        char error_message_[1000];
        sprintf(error_message_, "Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line_ %d - behind block=%d there are subsequent non-movable allocated blocks\n",
                ThisTask, func_, file_, line_, nr);
        terminate(error_message_);
      }
  }

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  if(n_ > FreeBytes)
  {
    dump_memory_table();
    char error_message_[1000];
    sprintf(error_message_, "Task=%d: at %s()/%s/line_ %d: Not enough memory in myremalloc_movable(n_=%g MB). previous=%g FreeBytes=%g MB\n",
            ThisTask, func_, file_, line_, n_ / (1024.0 * 1024.0), BlockSize[nr] / (1024.0 * 1024.0),
            FreeBytes / (1024.0 * 1024.0));
    terminate(error_message_);
  }

  size_t offset = n_ - BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
  {
    Table[i] += offset;

    *BasePointers[i] = *BasePointers[i] + offset;
  }

  FreeBytes -= n_;
  AllocatedBytes += n_;
  BlockSize[nr] = n_;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[nr];
}

void endrun(const int error_number_)
{
  if(error_number_)
  {
#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD, error_number_);
#endif
    exit(error_number_);
  }

#ifdef PARALLEL
  MPI_Finalize();
#endif
  exit(0);
}
