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

/** @file   peano.c
 *  @date   ????-2019
 *  @author Volker Springel
 *  @author Stefan Hilbert
 *
 *  @brief  Peano-Hilbert keys
 **/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


static const char quadrants[24][2][2][2] = {
  /* rot_x_=0, rot_y_=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rot_x_=1, rot_y_=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rot_x_=2, rot_y_=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rot_x_=3, rot_y_=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rot_x_=4, rot_y_=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rot_x_=5, rot_y_=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static const char rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static const char rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
  11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static const char rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static const char roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static const char sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };


int peano_hilbert_key(const int x_, const int y_, const int z_, const int bits_)
{
  int i_, bit_x_, bit_y_, bit_z_, mask_, quad_, rotation_;
  char sense_, rot_x_, rot_y_;
  int key_;

  mask_     = 1 << (bits_ - 1);
  key_      = 0;
  rotation_ = 0;
  sense_    = 1;

  for(i_ = 0; i_ < bits_; i_++, mask_ >>= 1)
  {
    bit_x_ = (x_ & mask_) ? 1 : 0;
    bit_y_ = (y_ & mask_) ? 1 : 0;
    bit_z_ = (z_ & mask_) ? 1 : 0;

    quad_ = quadrants[rotation_][bit_x_][bit_y_][bit_z_];

    key_ <<= 3;
    key_ += (sense_ == 1) ? (quad_) : (7 - quad_);

    rot_x_ = rotx_table[quad_];
    rot_y_ = roty_table[quad_];
    sense_ *= sense_table[quad_];

    while(rot_x_ > 0)
    {
      rotation_ = rotxmap_table[rotation_];
      rot_x_--;
    }

    while(rot_y_ > 0)
    {
      rotation_ = rotymap_table[rotation_];
      rot_y_--;
    }
  }

  return key_;
}
