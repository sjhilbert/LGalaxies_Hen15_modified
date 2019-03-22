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

/** @file typed_dynamic_array_inline.h
 * 
 * @author Stefan Hilbert (hilbert)
 *
 * 
 * 
 */
#ifndef TYPED_DYNAMIC_ARRAY_INLINE_H
#define TYPED_DYNAMIC_ARRAY_INLINE_H

#include <stddef.h>

#ifndef TYPED_DYNAMIC_ARRAY_BLOCK_SIZE
#define TYPED_DYNAMIC_ARRAY_BLOCK_SIZE (((size_t)1024) * ((size_t)1024))
#endif /* not defined TYPED_DYNAMIC_ARRAY_BLOCK_SIZE */


#ifndef TYPED_DYNAMIC_ARRAY_REALLOC_FACTOR
#define TYPED_DYNAMIC_ARRAY_REALLOC_FACTOR ((size_t)2)
#endif /* not defined TYPED_DYNAMIC_ARRAY_REALLOC_FACTOR */


#define typed_dynamic_array_define_type(data_type_name_, data_type_)
typedef struct typed_dynamic_array_of_##data_type_name_##type_
{
  size_t capacity;
  size_t size;
  data_type_* data;
} typed_dynamic_array_of_##data_type_name_##type


#define typed_dynamic_array_type(data_type_name_)
typed_dynamic_array_of_##data_type_name_##type


#define typed_dynamic_array_free(array_)
{
  array_->size     = 0;
  array_->capacity = 0;
  free(array_->data);
}


#define typed_dynamic_array_init(array_, size_, suggested_capacity_)
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


#define typed_dynamic_array_reserve(array_, suggested_capacity_)
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


#define typed_dynamic_array_resize(array_, size_)
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
