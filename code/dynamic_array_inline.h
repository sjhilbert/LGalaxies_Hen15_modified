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

/** @file   dynamic_array_inline.h
 *  @date   2018-2019
 *  @author Stefan Hilbert (hilbert)
 *
 *  @brief  dynamically allocated array (inspired by C++ STL vector)
 */
 
#ifndef DYNAMIC_ARRAY_INLINE_H
#define DYNAMIC_ARRAY_INLINE_H

#include <stddef.h>

#ifndef DYNAMIC_ARRAY_BLOCK_SIZE
#define DYNAMIC_ARRAY_BLOCK_SIZE (((size_t)1024) * ((size_t)1024))
#endif /* not defined DYNAMIC_ARRAY_BLOCK_SIZE */


#ifndef DYNAMIC_ARRAY_REALLOC_FACTOR
#define DYNAMIC_ARRAY_REALLOC_FACTOR ((size_t)2)
#endif /* not defined DYNAMIC_ARRAY_REALLOC_FACTOR */


typedef struct dynamic_array_type_
{
  size_t capacity;
  size_t size;
  void* data;
} dynamic_array_type;


static inline void 
dynamic_array_free(dynamic_array_type *array_)
{
  array_->size     = 0;
  array_->capacity = 0;
  free(array_->data);
}


static inline void 
dynamic_array_init(dynamic_array_type *array_, const size_t size_, const size_t suggested_capacity_)
{
  array_->size     = size_;
  if(size_ > suggested_capacity_)
  {  
    array_->capacity = ((size_ / DYNAMIC_ARRAY_BLOCK_SIZE) + (size_ % DYNAMIC_ARRAY_BLOCK_SIZE > 0)) * DYNAMIC_ARRAY_BLOCK_SIZE;
    if(!(array_->data = malloc(array_->capacity)))
    { terminate("dynamic_array_init: malloc failed"); }  }
  else if(suggested_capacity_ > 0)
  {
    array_->capacity = ((suggested_capacity_ / DYNAMIC_ARRAY_BLOCK_SIZE) + (suggested_capacity_ % DYNAMIC_ARRAY_BLOCK_SIZE > 0)) * DYNAMIC_ARRAY_BLOCK_SIZE;
    if(!(array_->data = malloc(array_->capacity)))
    { terminate("dynamic_array_init: malloc failed"); }
  }
  else
  { 
    array_->capacity = 0;
    array_->data     = NULL;
  }
}


static inline void 
dynamic_array_reserve(dynamic_array_type *array_, const size_t suggested_capacity_)
{
  if(suggested_capacity_ < array_->size)
  { /* do nothing */ }
  else if(suggested_capacity_  > 0)
  {
    array_->capacity = ((suggested_capacity_ / DYNAMIC_ARRAY_BLOCK_SIZE) + (suggested_capacity_ % DYNAMIC_ARRAY_BLOCK_SIZE > 0)) * DYNAMIC_ARRAY_BLOCK_SIZE;
    if(!(array_->data = realloc(array_->data, array_->capacity)))
    { terminate("dynamic_array_init: malloc failed"); }
  }
  else
  { 
    array_->capacity = 0;
    array_->data     = NULL;
  }
}


static inline void 
dynamic_array_resize(dynamic_array_type *array_, const size_t size_)
{
  array_->size = size_;
  if(array_->size > array_->capacity)
  {
    array_->capacity = (((size_t)(array_->size * DYNAMIC_ARRAY_REALLOC_FACTOR)) / DYNAMIC_ARRAY_BLOCK_SIZE + 1) * DYNAMIC_ARRAY_BLOCK_SIZE;
    if(!(array_->data = realloc(array_->data, array_->capacity)))
      terminate("dynamic_array_resize: realloc failed");
  }
}


static inline void  
dynamic_array_clear(dynamic_array_type *array_)
{ array_->size  = 0; }


static inline void  
dynamic_array_push_back(dynamic_array_type *array_, void *in_, size_t in_size_)
{ 
  if(array_->size + in_size_ > array_->capacity)
  {
    array_->capacity = (((size_t)((array_->size + in_size_) * DYNAMIC_ARRAY_REALLOC_FACTOR)) / DYNAMIC_ARRAY_BLOCK_SIZE + 1) * DYNAMIC_ARRAY_BLOCK_SIZE;
    if(!(array_->data = realloc(array_->data, array_->capacity)))
      terminate("dynamic_array_push_back: realloc failed");
  }
  memcpy(array_->data + array_->size, in_, in_size_);
  array_->size += in_size_;
}

  
#endif /* header guard */
