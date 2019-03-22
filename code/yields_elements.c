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
 * elements.c
 *
 *  Created on: 20.01.2012
 *      Author: robyates
 *
 *  Where individual chemical element history arrays are created and dealt with (like metals.c for metal history arrays)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


#ifdef INDIVIDUAL_ELEMENTS

/* Defined in allvars.h
struct elements
{
  float H;
  float He;
#ifndef MAINELEMENTS
  float Cb; //NOTE: Carbon (C) is stored as Cb here
  float N;
#endif
  float O;
#ifndef MAINELEMENTS
  float Ne;
#endif
  float Mg;
#ifndef MAINELEMENTS
  float Si;
  float S;
  float Ca;
#endif
  float Fe;
};
*/


struct elements elements_init()
{
  struct elements ele;
  ele.H=0.0;
  ele.He=0.0;
#ifndef MAINELEMENTS
  ele.Cb=0.0;
  ele.N=0.0;
#endif
  ele.O=0.0;
#ifndef MAINELEMENTS
  ele.Ne=0.0;
#endif
  ele.Mg=0.0;
#ifndef MAINELEMENTS
  ele.Si=0.0;
  ele.S=0.0;
  ele.Ca=0.0;
#endif
  ele.Fe=0.0;
  return(ele);
}


struct elements elements_add(const struct elements ele1, const struct elements ele2)
{
  struct elements ele;
  ele.H=ele1.H+ele2.H;
  ele.He=ele1.He+ele2.He;
#ifndef MAINELEMENTS
  ele.Cb=ele1.Cb+ele2.Cb;
  ele.N=ele1.N+ele2.N;
#endif
  ele.O=ele1.O+ele2.O;
#ifndef MAINELEMENTS
  ele.Ne=ele1.Ne+ele2.Ne;
#endif
  ele.Mg=ele1.Mg+ele2.Mg;
#ifndef MAINELEMENTS
  ele.Si=ele1.Si+ele2.Si;
  ele.S=ele1.S+ele2.S;
  ele.Ca=ele1.Ca+ele2.Ca;
#endif
  ele.Fe=ele1.Fe+ele2.Fe;

  return(ele);
}


struct elements elements_fraction(const struct elements ele2, const float fraction)
{
  struct elements ele;
  ele.H= fraction*ele2.H;
  ele.He= fraction*ele2.He;
#ifndef MAINELEMENTS
  ele.Cb= fraction*ele2.Cb;
  ele.N= fraction*ele2.N;
#endif
  ele.O= fraction*ele2.O;
#ifndef MAINELEMENTS
  ele.Ne= fraction*ele2.Ne;
#endif
  ele.Mg= fraction*ele2.Mg;
#ifndef MAINELEMENTS
  ele.Si= fraction*ele2.Si;
  ele.S= fraction*ele2.S;
  ele.Ca= fraction*ele2.Ca;
#endif
  ele.Fe= fraction*ele2.Fe;

  return(ele);
}


struct elements elements_add_fraction(const struct elements ele1, const struct elements ele2, const float fraction)
{
  struct elements ele;
  ele.H=ele1.H+fraction*ele2.H;
  ele.He=ele1.He+fraction*ele2.He;
#ifndef MAINELEMENTS
  ele.Cb=ele1.Cb+fraction*ele2.Cb;
  ele.N=ele1.N+fraction*ele2.N;
#endif
  ele.O=ele1.O+fraction*ele2.O;
#ifndef MAINELEMENTS
  ele.Ne=ele1.Ne+fraction*ele2.Ne;
#endif
  ele.Mg=ele1.Mg+fraction*ele2.Mg;
#ifndef MAINELEMENTS
  ele.Si=ele1.Si+fraction*ele2.Si;
  ele.S=ele1.S+fraction*ele2.S;
  ele.Ca=ele1.Ca+fraction*ele2.Ca;
#endif
  ele.Fe=ele1.Fe+fraction*ele2.Fe;

  return(ele);
}


void elements_add_to(struct elements *ele, const struct elements ele2)
{
  ele->H  += *ele2.H;
  ele->He += *ele2.He;
#ifndef MAINELEMENTS
  ele->Cb += *ele2.Cb;
  ele->N  += *ele2.N;
#endif
  ele->O  += *ele2.O;
#ifndef MAINELEMENTS
  ele->Ne += *ele2.Ne;
#endif
  ele->Mg += *ele2.Mg;
#ifndef MAINELEMENTS
  ele->Si += *ele2.Si;
  ele->S  += *ele2.S;
  ele->Ca += *ele2.Ca;
#endif
  ele->Fe += *ele2.Fe;
}


void elements_add_fraction_to(struct elements *ele, const struct elements ele2, const float fraction)
{
  ele->H  += fraction*ele2.H;
  ele->He += fraction*ele2.He;
#ifndef MAINELEMENTS
  ele->Cb += fraction*ele2.Cb;
  ele->N  += fraction*ele2.N;
#endif
  ele->O  += fraction*ele2.O;
#ifndef MAINELEMENTS
  ele->Ne += fraction*ele2.Ne;
#endif
  ele->Mg += fraction*ele2.Mg;
#ifndef MAINELEMENTS
  ele->Si += fraction*ele2.Si;
  ele->S  += fraction*ele2.S;
  ele->Ca += fraction*ele2.Ca;
#endif
  ele->Fe += fraction*ele2.Fe;
}


void elements_deduct_from(struct elements *ele, const struct elements ele2)
{
  ele->H  -= *ele2.H;
  ele->He -= *ele2.He;
#ifndef MAINELEMENTS
  ele->Cb -= *ele2.Cb;
  ele->N  -= *ele2.N;
#endif
  ele->O  -= *ele2.O;
#ifndef MAINELEMENTS
  ele->Ne -= *ele2.Ne;
#endif
  ele->Mg -= *ele2.Mg;
#ifndef MAINELEMENTS
  ele->Si -= *ele2.Si;
  ele->S  -= *ele2.S;
  ele->Ca -= *ele2.Ca;
#endif
  ele->Fe -= *ele2.Fe;
}


void elements_deduct_fraction_from(struct elements *ele, const struct elements ele2, const float fraction)
{
  ele->H  -= fraction*ele2.H;
  ele->He -= fraction*ele2.He;
#ifndef MAINELEMENTS
  ele->Cb -= fraction*ele2.Cb;
  ele->N  -= fraction*ele2.N;
#endif
  ele->O  -= fraction*ele2.O;
#ifndef MAINELEMENTS
  ele->Ne -= fraction*ele2.Ne;
#endif
  ele->Mg -= fraction*ele2.Mg;
#ifndef MAINELEMENTS
  ele->Si -= fraction*ele2.Si;
  ele->S  -= fraction*ele2.S;
  ele->Ca -= fraction*ele2.Ca;
#endif
  ele->Fe -= fraction*ele2.Fe;
}

void elements_print(char s[],struct elements ele)
{
  printf("%s.H [Msun]  = %.2f\n",s,ele.H);
  printf("%s.He [Msun] = %.2f\n",s,ele.He);
#ifndef MAINELEMENTS
  printf("%s.Cb [Msun] = %.2f\n",s,ele.Cb);
  printf("%s.N [Msun]  = %.2f\n",s,ele.N);
#endif
  printf("%s.O [Msun]  = %.2f\n",s,ele.O);
#ifndef MAINELEMENTS
  printf("%s.Ne [Msun] = %.2f\n",s,ele.Ne);
#endif
  printf("%s.Mg [Msun] = %.2f\n",s,ele.Mg);
#ifndef MAINELEMENTS
  printf("%s.Si [Msun] = %.2f\n",s,ele.Si);
  printf("%s.S [Msun]  = %.2f\n",s,ele.S);
  printf("%s.Ca [Msun] = %.2f\n",s,ele.Ca);
#endif
  printf("%s.Fe [Msun] = %.2f\n",s,ele.Fe);
  return;
}


double elements_total(struct elements ele)
{
#ifndef MAINELEMENTS
  return(ele.H+ele.He+ele.Cb+ele.N+ele.O+ele.Ne+ele.Mg+ele.Si+ele.S+ele.Ca+ele.Fe);
#else
  return(ele.H+ele.He+ele.O+ele.Mg+ele.Fe);
#endif
}


double metal_elements_total(struct elements ele)
{
#ifndef MAINELEMENTS
  return(ele.Cb+ele.N+ele.O+ele.Ne+ele.Mg+ele.Si+ele.S+ele.Ca+ele.Fe);
#else
  return(ele.O+ele.Mg+ele.Fe);
#endif
}

#endif //INDIVIDUAL_ELEMENTS
