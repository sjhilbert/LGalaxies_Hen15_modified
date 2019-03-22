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

/** @file aux_macro_functions.h
 * 
 * @author Stefan Hilbert (hilbert)
 * 
 */
#ifndef AUX_FUNCTIONS_INLINE_H
#define AUX_FUNCTIONS_INLINE_H

#include <stddef.h>

/** @brief sets all entries in array to given value */
#define set_array_to(arr_, value_) \
do { size_t arr_##idx_; for(arr_##idx_ = 0; arr_##idx_ < (sizeof arr_ / sizeof *arr_); arr_##idx_++) { arr_[arr_##idx_] = value_; } } while (0)
  

/** @brief finds interpolation index for interpolation from table
 * 
 * finds interpolation index i_ for interpolation based on values
 * in array arr_:
 * - if value_ <= arr_[i_beg_    ], return i_ = i_beg_    ,
 * - if value_ >= arr_[i_end_ - 2], return i_ = i_end_ - 2,
 * - else, return i_ such that arr_[i_] <= value_ < arr_[i_ + 1].
 * 
 * requirements:
 * - i_beg_ >= 0,
 * - i_end_ >= 2,
 * - arr_[i_] <=  arr_[i_ + 1] , i.e. assumes ascending order.
 * 
 *  uses ascending linear scan.
 */
#define locate_interpolation_index_ascend_linear(i_, i_beg_, i_end_, value_, arr_) \
do { i_ = i_beg_; while((i_ < (i_end_ - 2)) && (arr_[i_ + 1] < value_)) ++i_; } while (0)
  

/** @brief finds interpolation index for interpolation from table
 * 
 * finds interpolation index i_ for interpolation based on values
 * in array arr_:
 * - if value_ <= arr_[i_beg_    ], return i_ = i_beg_    ,
 * - if value_ >= arr_[i_end_ - 2], return i_ = i_end_ - 2,
 * - else, return i_ such that arr_[i_] <= value_ < arr_[i_ + 1].
 * 
 * requirements:
 * - i_beg_ >= 0,
 * - i_end_ >= 2,
 * - arr_[i_] <=  arr_[i_ + 1] , i.e. assumes ascending order.
 * 
 *  uses bisection
 */
#define locate_interpolation_index_ascend_bisect(i_, i_beg_, i_end_, value_, arr_) \
do                                                                                 \
{                                                                                  \
  size_t i_1_ = i_beg_;                                                            \
  size_t i_2_ = i_end_ - 1;                                                        \
  while(i_ = (i_1_ + i_2_) / 2, i_1_ + 1 < i_2_)                                   \
  {                                                                                \
    if(arr_[i_] > value_)                                                          \
    { i_2_ = i_; }                                                                 \
    else                                                                           \
    { i_1_ = i_; }                                                                 \
  }                                                                                \
} while (0)


#define locate_interpolation_index_ascend(i_, i_beg_, i_end_, value_, arr_) \
locate_interpolation_index_ascend_linear(i_, i_beg_, i_end_, value_, arr_)





#endif /* header guard */
