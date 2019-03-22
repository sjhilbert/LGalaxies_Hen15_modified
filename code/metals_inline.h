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

/*
 * metals.c
 *
 * Routines to handle the creation and addition of containers for metals.
 *
 * The routines are needed even if the metallicity is a single float,
 * because the lines for creation and addition of floats in the code
 * have been replaced by function calls (this to avoid a multitude of
 * #ifdef DETAILED_METALS_AND_MASS_RETURN statements throughout the code).
 *
 * In the description below metal is either:
 *   float
 *   struct metals
 * It would be relatively straight forward to add a metallicity array option also.
 *
 * metal metals_add(metal m1, metal m2, float fraction);
 *  Function returns m1+fraction*m2 for each metal component.
 *
 * metal metals_init();
 *  Function returns 0.0 for each metal component.
 *
 * void metals_print(char [s] ,metal m);
 *  Prints out the value of each metal component, preceded by string s.
 *
 * float metals_total(metal m);
 *  Function returns the total of all components of metal.
 */

#include <stdio.h>

#include "allvars.h"
#include "proto.h"

#ifndef METALS_INLINE_H
#define METALS_INLINE_H

#ifdef DETAILED_METALS_AND_MASS_RETURN

static inline
struct metals metals_init()
{
  struct metals m;
  m.type1a=0.;
  m.type2=0.;
  m.agb=0.;
  return m;
}

static inline
struct metals metals_add(const struct metals m1,
                         const struct metals m2)
{
  struct metals m;
  m.type1a=m1.type1a+m2.type1a;
  m.type2=m1.type2+m2.type2;
  m.agb=m1.agb+m2.agb;
  return m;
}


static inline
struct metals metals_fraction(const struct metals m2,
                              const float fraction)
{
  struct metals m;
  m.type1a=fraction*m2.type1a;
  m.type2=fraction*m2.type2;
  m.agb=fraction*m2.agb;
  return m;
}


static inline
struct metals metals_add_fraction(const struct metals m1,
                                 const struct metals m2,
                                 const float fraction)
{
  struct metals m;
  m.type1a=m1.type1a+fraction*m2.type1a;
  m.type2=m1.type2+fraction*m2.type2;
  m.agb=m1.agb+fraction*m2.agb;
  return m;
}


static inline
void metals_add_to(struct *metals m,
                   const struct metals m2)
{
  m->type1a += m2.type1a;
  m->type2  += m2.type2;
  m->agb    += m2.agb;
}


static inline
void metals_add_fraction_to(struct *metals m,
                            const struct metals m2,
                            const float fraction)
{
  m->type1a += fraction * m2.type1a;
  m->type2  += fraction * m2.type2;
  m->agb    += fraction * m2.agb;
}


static inline
void metals_deduct_from(struct *metals m,
                        const struct metals m2)
{
  m->type1a -= m2.type1a;
  m->type2  -= m2.type2;
  m->agb    -= m2.agb;
}


static inline
void metals_deduct_fraction_from(struct *metals m,
                                 const struct metals m2,
                                 const float fraction)
{
  m->type1a -= fraction * m2.type1a;
  m->type2  -= fraction * m2.type2;
  m->agb    -= fraction * m2.agb;
}


static inline
void metals_multiply_by(struct *metals m,
                       const float factor)
{
  m->type1a *= factor;
  m->type2  *= factor;
  m->agb    *= factor;
}


static inline
void metals_divide_by(struct *metals m,
                      const float factor)
{
  m->type1a /= factor;
  m->type2  /= factor;
  m->agb    /= factor;
}


static inline
float metals_total(struct metals m)
{ return(m.type1a+m.type2+m.agb); }


static inline
void metals_print(char s[], const struct metals m)
{
  printf("%s.type1a [Msun] = %.2f\n",s,m.type1a*1.0e10 * inv_Hubble_h);
  printf("%s.type2 [Msun]  = %.2f\n",s,m.type2*1.0e10 * inv_Hubble_h);
  printf("%s.agb  [Msun]   = %.2f\n",s,m.agb*1.0e10 * inv_Hubble_h);
}

#else /* not defined DETAILED_METALS_AND_MASS_RETURN */

// The following mimics the original code with a single metallicity

static inline
float metals_init()
{ return(0.); }


static inline
float metals_add(const float m1,
                 const float m2)
{ return m1 + m2; }


static inline
float metals_fraction(const float m2,
                      const float fraction)
{ return fraction * m2; }


static inline
float metals_add_fraction(const float m1,
                          const float m2,
                          const float fraction)
{ return m1 + fraction * m2; }


static inline
void metals_add_to(float *m,
                   const float m2)
{ *m += m2; }


static inline
void metals_add_fraction_to(float *m,
                            const float m2,
                            const float fraction)
{ *m += fraction * m2; }


static inline
void metals_deduct_from(float *m,
                       const float m2)
{ *m -= m2; }


static inline
void metals_deduct_fraction_from(float *m,
                                 const float m2,
                                 const float fraction)
{ *m -= fraction * m2; }


static inline
void metals_multiply_by(float *m,
                        const float factor)
{ *m *= factor; }


static inline
void metals_divide_by(float *m,
                      const float factor)
{ *m /= factor; }


static inline
float metals_total(const float m)
{ return m; }


static inline
void metals_print(char s[], const float m)
{ printf("%s=%f\n",s, m); }

#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */

#endif /* header guard */

