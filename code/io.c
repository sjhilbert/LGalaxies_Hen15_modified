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

/** @file   io.c
 *  @date   2016-2019 
 *  @author ?
 *  @author Stefan Hilbert
 *
 *  @brief basic file io helper functions
 */
 
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
 
/** @brief Reading routine, either from a file into a structure or
 *         from a pointer to a structure.
 **/
size_t myfread(void *ptr, const size_t size_, const size_t n_memb_, FILE * stream_)
{
  size_t n_read_ = 0;
  if(size_ * n_memb_ > 0)
  {
    if((n_read_ = fread(ptr, size_, n_memb_, stream_)) != n_memb_)
    {
      if(feof(stream_))
        printf("I/O error (fread) has occured: end of file\n");
      else
        printf("I/O error (fread) has occured: %s\n", strerror(errno));
      fflush(stdout);
      terminate("read error");
    }
  }
  return n_read_;
}


/** @brief Writing routine, either to a file from a structure or
 *         to a pointer from a structure.
 **/
size_t myfwrite(void *ptr, const size_t size_, const size_t n_memb_, FILE * stream_)
{
  size_t n_written_ = 0;

  if(size_ * n_memb_ > 0)
  {
    if((n_written_ = fwrite(ptr, size_, n_memb_, stream_)) != n_memb_)
    {
      printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
      fflush(stdout);
      terminate("write error");
    }
  }
  return n_written_;
}


#ifdef WRITE_LARGE_DATA_IN_CHUNKS
#ifndef MAXIMUM_WRITE_SIZE_IN_BYTES  
#define MAXIMUM_WRITE_SIZE_IN_BYTES 2147483647ull
/* = 2^31 - 1, i.e. the largest possible 32-bit signed int value */
#endif /* not defined MAXIMUM_WRITE_SIZE_IN_BYTES */
#endif /* defined WRITE_LARGE_DATA_IN_CHUNKS */


/** @brief Writing routine, either to a file from a structure or
 *         to a pointer from a structure.
 *
 * @todo  integer multiply overflow check
 **/
size_t myfwrite_large_data(void *ptr, const size_t size_, const size_t n_memb_, FILE * stream_)
{
  if(size_ == 0 || n_memb_ == 0)
  { return 0; }
#ifdef WRITE_LARGE_DATA_IN_CHUNKS
  else if(size_ * n_memb_ <= MAXIMUM_WRITE_SIZE_IN_BYTES)
#endif /* defined WRITE_LARGE_DATA_IN_CHUNKS */
  {
    size_t n_written_ = 0;
    if((n_written_ = fwrite(ptr, size_, n_memb_, stream_)) != n_memb_)
    {
      printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
      fflush(stdout);
      terminate("write error");
    }
    return n_written_;
  }
#ifdef WRITE_LARGE_DATA_IN_CHUNKS    
  else /* now we write in chunks */
  {
    const size_t tot_n_bytes_to_write_ = size_* n_memb_;
    size_t curr_, n_bytes_to_write_, tot_n_bytes_written_ = 0;
    for(curr_ = 0; curr_ < tot_n_bytes_to_write_; curr_ += MAXIMUM_WRITE_SIZE_IN_BYTES)
    {
      n_bytes_to_write_ = min(MAXIMUM_WRITE_SIZE_IN_BYTES, tot_n_bytes_to_write_ - curr_);
      tot_n_bytes_written_ += fwrite(ptr , 1, n_bytes_to_write_, stream_);
      ptr += n_bytes_to_write_;
      /* fflush(stream_); */
    }
    
    if(tot_n_bytes_written_ != tot_n_bytes_to_write_)
    {
      printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
      fflush(stdout);
      terminate("write error");
    }
    return n_memb_;
  }
#endif /* defined WRITE_LARGE_DATA_IN_CHUNKS */ 
}


/** @brief Writing routine (repeating one value) **/
size_t myffill(void *ptr_, const size_t size_, const size_t n_memb_, FILE * stream_)
{
  size_t n_written_ = 0;
  for(n_written_ = 0; n_written_ < n_memb_; n_written_++)
    if(1 != fwrite(ptr_, size_, 1, stream_))
    {
      printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
      fflush(stdout);
      terminate("write error");
    }
  return n_written_;
}


/** @brief moving stream i/o position **/
int myfseek(FILE * stream_, const long offset_, const int whence_)
{
  if(fseek(stream_, offset_, whence_))
  {
    printf("I/O error (fseek) has occured: %s\n", strerror(errno));
    fflush(stdout);
    terminate("fseek error");
  }
  return 0;
}
