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

#ifndef GASDEV_INLINE_H
#define GASDEV_INLINE_H

#include <math.h>

#include "allvars.h"
#include "proto.h"

/** @file gasdev_inline.h
 *  @brief hand-made rng for normal distributed random numbers 
 */

/** @brief hand made rng (replace by gsl version?) */
static inline float 
ran1(long *idum)
{
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
  
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if(*idum <= 0 || !iy)
    {
      if(-(*idum) < 1)
                *idum = 1;
      else
                *idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--)
                {
                  k = (*idum) / IQ;
                  *idum = IA * (*idum - k * IQ) - IR * k;
                  if(*idum < 0)
                    *idum += IM;
                  if(j < NTAB)
                    iv[j] = *idum;
                }
      iy = iv[0];
    }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if(*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
  
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX  
}


/** @brief computes a gaussian random deviate to calculate a random
 *         inclination for extinction. */
static inline float
gasdev(long *idum)
{
  static int iset = 0;
  static float gset;
  float fac, rsq, v1, v2;

  if(iset == 0)
    {
      do
        {
          v1 = 2.0 * ran1(idum) - 1.0;
          v2 = 2.0 * ran1(idum) - 1.0;
          rsq = v1 * v1 + v2 * v2;
        }
      while(rsq >= 1.0 || rsq == 0.0);
      fac = sqrt(-2.0 * log(rsq) / rsq);
      gset = v1 * fac;
      iset = 1;
      return v2 * fac;
    }
  else
    {
      iset = 0;
      return gset;
    }
}


#endif /* header guard */
