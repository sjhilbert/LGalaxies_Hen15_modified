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

/** @file aux_functions_inline.h
 * 
 * @author Stefan Hilbert (hilbert)
 *
 *  so far mostly macro magic for some basic compile time polymorphism 
 *  and simple generic algorithms 
 * 
 */
#ifndef AUX_FUNCTIONS_INLINE_H
#define AUX_FUNCTIONS_INLINE_H

#include <stddef.h>

/** @brief sets all entries in mem range to given value
 *
 * @param [in]    ptr_   pointer to begin of mem range
 * @param [in]    value_ value the elements in mem will be set to
 * @param [in]    count_ number of elements to be set
 */
#define set_mem_to(ptr_, value_, count_) \
do { size_t ptr_##idx_; for(ptr_##idx_ = 0; ptr_##idx_ < count_; ++ptr_##idx_) { ptr_[ptr_##idx_] = value_; } } while(0)
  

/** @brief sets all entries in array to given value
 *
 * @param [inout] arr_   array to be set
 * @param [in]    value_ value the array elements will be set to
 *
 * @warning won't work with pointers or arrays decayed to pointers
 */
#define set_array_to(arr_, value_) set_mem_to(arr_, value_, (sizeof arr_ / sizeof *arr_))


/** @brief compute interpolation tables */ 
#define set_interpolation_tables(i_, i_beg_, i_end_, arg_beg_, arg_step_, f_, arg_arr_, res_arr_) \
do{                                                                                               \
  for(i_ = i_beg_; i_ < i_end_; ++i_)                                                             \
  {                                                                                               \
    arg_arr_[i_] =    arg_beg_ + arg_step_ * i_;                                                  \
    res_arr_[i_] = f_(arg_beg_ + arg_step_ * i_);                                                 \
  }                                                                                               \
}while(0)  
  
  
/** @brief finds interpolation index for interpolation from table
 *
 * finds interpolation index i_ for interpolation based on values
 * in array arg_arr_:
 * - if arg_value_ <= arg_arr_[i_beg_    ], return i_ = i_beg_    ,
 * - if arg_value_ >= arg_arr_[i_end_ - 2], return i_ = i_end_ - 2,
 * - else, return i_ such that arg_arr_[i_] <= arg_value_ < arg_arr_[i_ + 1].
 * 
 * @param [out]   i_         interpolation index returned 
 * @param [in]    i_beg_     lower bound for interpolation index 
 * @param [in]    i_end_     upper bound for interpolation index 
 * @param [in]    arg_value_ function argument value for which interpolation index is requested
 * @param [in]    arg_arr_   array containing function argument values for interpolation
 * @param [in]    order_rel_ ordering in array values, e.g. '<' or '<=' for ascending order
 *
 * @pre
 * - i_beg_ >= 0,
 * - i_end_ >= 2,
 * - arg_arr_[i_] <=  arg_arr_[i_ + 1] for ascending  order (order_rel_ = '<'),
 * - arg_arr_[i_] >=  arg_arr_[i_ + 1] for descending order (order_rel_ = '>').
 * 
 *  uses ascending linear scan.
 */
#define locate_interpolation_index_linear(i_, i_beg_, i_end_, arg_value_, arg_arr_, order_rel_) \
do { i_ = i_beg_; while((i_ < (i_end_ - 2)) && (arg_arr_[i_ + 1] order_rel_##= arg_value_)) ++i_; } while(0)
  

/** @brief finds interpolation index for interpolation from table
 *
 * finds interpolation index i_ for interpolation based on values
 * in array arg_arr_:
 * - if arg_value_ <= arg_arr_[i_beg_    ], return i_ = i_beg_    ,
 * - if arg_value_ >= arg_arr_[i_end_ - 2], return i_ = i_end_ - 2,
 * - else, return i_ such that arg_arr_[i_] <= arg_value_ < arg_arr_[i_ + 1].
 * 
 * @param [out]   i_         interpolation index returned 
 * @param [in]    i_beg_     lower bound for interpolation index 
 * @param [in]    i_end_     upper bound for interpolation index 
 * @param [in]    arg_value_ function argument value for which interpolation index is requested
 * @param [in]    arg_arr_   array containing function argument values for interpolation
 * @param [in]    order_rel_ ordering in array values, e.g. '<' or '<=' for ascending order
 *
 * @pre
 * - i_beg_ >= 0,
 * - i_end_ >= 2,
 * - arg_arr_[i_] <=  arg_arr_[i_ + 1] for ascending  order (order_rel_ = '<'),
 * - arg_arr_[i_] >=  arg_arr_[i_ + 1] for descending order (order_rel_ = '>').
 * 
 *  uses bisection
 */
#define locate_interpolation_index_bisect(i_, i_beg_, i_end_, arg_value_, arg_arr_, order_rel_) \
do{                                                                                             \
  size_t i_1_ = i_beg_;                                                                         \
  size_t i_2_ = i_end_ - 1;                                                                     \
  while(i_ = (i_1_ + i_2_) / 2, i_1_ + 1 < i_2_)                                                \
  {                                                                                             \
    if(arg_arr_[i_] order_rel_##= arg_value_)                                                   \
    { i_1_ = i_; }                                                                              \
    else                                                                                        \
    { i_2_ = i_; }                                                                              \
  }                                                                                             \
}while(0)

  
/** @brief finds interpolation index for interpolation from table */
#define locate_interpolation_index(i_, i_beg_, i_end_, value_, arr_, order_rel_, method_) \
locate_interpolation_index_##method_(i_, i_beg_, i_end_, value_, arr_, order_rel_)


/** @brief finds interpolation and fractions for interpolation from table  */
#define locate_interpolation_index_and_fraction(i_, i_beg_, i_end_, arg_value_, arg_arr_, f_i_, f_i_p_1_, order_rel_, method_)              \
do{                                                                                                                                         \
  if(arg_value_ order_rel_##= arg_arr_[i_beg_])                                                                                             \
  { i_ = i_beg_    ; f_i_ = 1, f_i_p_1_ = 0; }                                                                                              \
  else if(!(arg_value_ order_rel_ arg_arr_[i_end_ - 1]))                                                                                    \
  { i_ = i_end_ - 2; f_i_ = 0, f_i_p_1_ = 1; }                                                                                              \
  else                                                                                                                                      \
  {                                                                                                                                         \
    locate_interpolation_index(i_, i_beg_, i_end_, arg_value_, arg_arr_, order_rel_, method_);                                              \
    f_i_p_1_ = (arg_value_ - arg_arr_[i_]) / (arg_arr_[i_ + 1] - arg_arr_[i_]);                                                             \
    f_i_     = 1 - f_i_p_1_;                                                                                                                \
  }                                                                                                                                         \
}while(0)
  

/** @brief finds interpolation and fractions for interpolation from table  */
#define locate_interpolation_index_and_fraction_bf(i_, i_beg_, i_end_, arg_value_, arg_arr_, f_i_, f_i_p_1_, order_rel_, method_, i_beg_f_) \
do{                                                                                                                                         \
  if(arg_value_ order_rel_##= arg_arr_[i_beg_])                                                                                             \
  { i_ = i_beg_    ; f_i_ = 1, f_i_p_1_ = 0; }                                                                                              \
  else if(!(arg_value_ order_rel_ arg_arr_[i_end_ - 1]))                                                                                    \
  { i_ = i_end_ - 2; f_i_ = 0, f_i_p_1_ = 1; }                                                                                              \
  else                                                                                                                                      \
  {                                                                                                                                         \
    locate_interpolation_index(i_, i_beg_f_(arg_value_), i_end_, arg_value_, arg_arr_, order_rel_, method_);                                \
    f_i_p_1_ = (arg_value_ - arg_arr_[i_]) / (arg_arr_[i_ + 1] - arg_arr_[i_]);                                                             \
    f_i_     = 1 - f_i_p_1_;                                                                                                                \
  }                                                                                                                                         \
}while(0)                                                                                                                                   

  
/** @brief linear interpolation from tables  */                                                                                             
#define linear_interpolate(i_, i_beg_, i_end_, arg_value_, arg_arr_, res_arr_, res_value_, order_rel_, method_)                             \
do{                                                                                                                                         \
  if(arg_value_ order_rel_##= arg_arr_[i_beg_])                                                                                             \
  { res_value_ = res_arr_[i_beg_]; }                                                                                                        \
  else if(!(arg_value_ order_rel_ arg_arr_[i_end_ - 1]))                                                                                    \
  { res_value_ = res_arr_[i_end_ - 1]; }                                                                                                    \
  else                                                                                                                                      \
  {                                                                                                                                         \
    locate_interpolation_index(i_, i_beg_, i_end_, arg_value_, arg_arr_, order_rel_, method_);                                              \
    res_value_ = res_arr_[i_] + (res_arr_[i_ + 1] - res_arr_[i_]) * (arg_value_ - arg_arr_[i_]) / (arg_arr_[i_ + 1] - arg_arr_[i_]);        \
  }                                                                                                                                         \
}while(0) 
  

/** @brief linear interpolation from tables  */                                                                                             
#define linear_interpolate_bf(i_, i_beg_, i_end_, arg_value_, arg_arr_, res_arr_, res_value_, order_rel_, method_, i_beg_f_)                \
do{                                                                                                                                         \
  if(arg_value_ order_rel_##= arg_arr_[i_beg_])                                                                                             \
  { res_value_ = res_arr_[i_beg_]; }                                                                                                        \
  else if(!(arg_value_ order_rel_ arg_arr_[i_end_ - 1]))                                                                                    \
  { res_value_ = res_arr_[i_end_ - 1]; }                                                                                                    \
  else                                                                                                                                      \
  {                                                                                                                                         \
    locate_interpolation_index(i_, i_beg_f_(arg_value_), i_end_, arg_value_, arg_arr_, order_rel_, method_);                                \
    res_value_ = res_arr_[i_] + (res_arr_[i_ + 1] - res_arr_[i_]) * (arg_value_ - arg_arr_[i_]) / (arg_arr_[i_ + 1] - arg_arr_[i_]);        \
  }                                                                                                                                         \
}while(0)

  
#endif /* header guard */
