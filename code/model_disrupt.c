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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file model_disrupt.c
 *  @brief model_disrupt.c checks if a type 2 satellite galaxy should
 *         or not be disrupted due to tidal forces
 *
 *  This routine takes into account the tidal effects that satellite
 *  galaxies experience while orbiting a central companion. Since the
 *  baryonic component is more compact and denser that the dark matter
 *  it assumes that only type 2 satellites are affected (those that have
 *  already lost their dark matter halo). It assumes that the disruption
 *  is complete and instantaneous with all the satellite material being
 *  transferred into the central galaxy.
 *
 *  The satellite is assumed to orbit a singular isothermal potential:
 *
 *  \f$\phi(R)=V^2_{\rm{vir}}\rm{ln} R\f$ (Eq. 28 Guo2010)
 *
 *  Assuming conservation of energy and angular momentum along the
 *  orbit, its pericentric distance (closest point to the centre along
 *  the orbit) can be estimated from:
 *
 *  \f$\left(\frac{R}{R_{\rm{peri}}}\right)^2=
 *  \frac{lnR/R_{\rm{peri}}+\frac{1}{2}(V/V_{\rm{vir}})^2}
 *  {\frac{1}{2}(V_{\rm{t}}/V_{\rm{vir}})^2}\f$
 *  (Eq. 29 Guo2010).
 *
 *  The main halo density at this point is compared with the baryonic
 *  mass (cold gas + stellar) density of the satellite within its half
 *  mass radius. If
 *
 *  \f$ \frac{M_{\rm{DM,halo}}(R_{\rm{peri}})}{R^3_{\rm{peri}}}\equiv
 *  \rho_{\rm{DM,halo}}>
 *  \rho_{\rm{sat}}\equiv\frac{M_{\rm{sat}}}{R^3_{\rm{sat,half}}}\f$
 *  (Eq. 30 Guo2010)
 *
 *  the galaxy is disrupted.
 *
 *  */


// static inline double get_isothermal_mass(const double Mvir, const double Rvir, const double dr)
// {
//   return Mvir/Rvir * dr;
// }


/** @brief Returns the mass of a disk within a given radius in units of the scale length
 *         Disk profile -> exponential */
static inline double
get_disk_mass_for_radius(const double x)
{ return 1.-(1.+x)*exp(-x); }


/** @brief Returns the mass of a bulge at a certain radius.
 *         Bulge profile -> de Vaucouleurs type r^{1/4} law */
static inline double 
get_bulge_mass_for_radius(const double x)
{ return x/(1.+x); }


/** @brief Calculates the distance of the satellite to the pericentre of the
  *        main dark matter halo. */
static inline double
get_peri_radius_for_galaxy(const int p, const int centralgal)
{
  int i;
  double a, b, v[3], r[3], x, x0;
  for(i = 0; i < 3; i++)
    {
      r[i] = wrap(Gal[p].Pos[i]-Gal[centralgal].Pos[i],BoxSize);
      r[i] /= (1 + ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
      v[i] = Gal[p].Vel[i] - Gal[centralgal].Vel[i];
    }

  b = 1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / pow2(Gal[centralgal].Vvir);
  a = 1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] -
                pow2(r[0] * v[0] + r[1] * v[1] + r[2] * v[2])/
                    (r[0] * r[0] + r[1] * r[1] + r[2] * r[2])) / pow2(Gal[centralgal].Vvir);

  x = sqrt(b / a);
  x0 = 1000;
  while(abs(x0 - x) >= 1.e-8)
    {
      x0 = x;
      x = sqrt((log(x0) + b) / a);
    }
  if(x == 0)
    {
      terminate("wrong in peri_radius \n");
    }

  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/x;
}


/** @brief Calculates the half mass radius of satellite galaxies */
static inline double 
get_sat_radius_for_galaxy(const int p)
{
  double r, rd, rb, Mdisk, rmax, M;
  double Mgas, Mbulge, rgd, dr, totmass;
  int ii;
  #define SAT_RADIUS_RMIN 5e-7
  #define SAT_RADIUS_N 100

  r=0.;
  rgd = Gal[p].GasDiskRadius/3.;
  rd=Gal[p].StellarDiskRadius/3.;
  rb=Gal[p].BulgeSize;
  Mgas = Gal[p].ColdGas;
  Mdisk=Gal[p].DiskMass;
  Mbulge = Gal[p].BulgeMass;
  totmass = Mgas+Mdisk+Mbulge;

  rmax=max(rb,1.68*max(rd,rgd));
  if (rmax < 2.*SAT_RADIUS_RMIN)
        return(rmax);
  dr=(rmax-SAT_RADIUS_RMIN)/(float)SAT_RADIUS_N;

    /* increases the search radius until it encompasses half the total mass taking
     * into account the stellar disk, stellar bulge and cold gas disk. */
  ii = 0;
  do {
      // Not sure that we need the 0.5 here - it's all a matter of definition
      r = (SAT_RADIUS_RMIN) + (ii+0.5)* dr;
      M = Mgas*get_disk_mass_for_radius(r/rgd)+Mdisk*get_disk_mass_for_radius(r/rd);

#ifndef GUO10
#ifndef GUO13
#ifndef HENRIQUES13
      if(Mbulge>0.)
#endif
#endif
#endif
        M +=Mbulge*get_bulge_mass_for_radius(r/rb);

      ii++;
      if(ii > 1000) terminate ("couldn't find half mass radius");
  }
  while(M < 0.5*totmass);

  return (r);
}


/**  @brief checks if a type 2 satellite galaxy should be disrupted
 *   due to tidal forces */
void disrupt(const int p)
{
  double rho_sat, rho_cen;
  double cen_mass, r_sat, radius;
  int centralgal;

  /* If the main halo density at the pericentre (closest point in the orbit
   * to the central galaxy)is larger than the satellite's density at the
   * half mass radius, the satellite is completely disrupted. Note that since
   * the satellite is a type 2 the only mass components remaining and
   * contributing to the density are the cold gas and stellar mass. */

  centralgal=Gal[p].CentralGal;
 
  mass_checks("Top of disrupt",centralgal);
  mass_checks("Top of disrupt",p);

  /* Radius calculated at the peri-point */
  radius = get_peri_radius_for_galaxy(p, centralgal);
  if (radius < 0) {
   terminate("must be wrong \n");
  }

  /* Calculate the density of the main central halo at radius (the peri-centre).
   * The tidal forces are caused by the dark matter of the main halo, hence Mvir
   * is used. Assume isothermal. */
  cen_mass=Gal[centralgal].Mvir*radius/Gal[centralgal].Rvir;
  rho_cen=cen_mass/pow3(radius);

  /* Calculate the density of the satellite's baryonic material */
  if (Gal[p].DiskMass+Gal[p].BulgeMass>0)
  {
    /* Calculate the rho according to the real geometry */
    r_sat = get_sat_radius_for_galaxy(p);
    rho_sat=(Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ColdGas)/pow3(r_sat);
  }
  else
    rho_sat=0.0;

  /* If density of the main halo is larger than that of the satellite baryonic
   * component, complete and instantaneous disruption is assumed. Galaxy becomes
   * a type 3 and all its material is transferred to the central galaxy. */
  if (rho_cen > rho_sat)
  {
    Gal[p].Type = 3;
#ifdef GALAXYTREE
    int q;
    q = Gal[Gal[p].CentralGal].FirstProgGal;
    if (q >= 0)
    {
      // add progenitors of Gal[p] to the list of progentitors of Gal[p].CentralGal
      while (GalTree[q].NextProgGal >= 0)
      	q = GalTree[q].NextProgGal;
	
      GalTree[q].NextProgGal = Gal[p].FirstProgGal;
      
      if(GalTree[q].NextProgGal >= NGalTree)
	{
	  printf("q=%d p=%d GalTree[q].NextProgGal=%d NGalTree=%d\n",
		 q, p, GalTree[q].NextProgGal, NGalTree);
	  terminate("problem");
	}
    }

    if(q < 0)
    	terminate("this shouldn't happen");
	
    q = GalTree[q].NextProgGal;

    if(q < 0)
    	terminate("inconsistency");

    if(HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q)
    	terminate("inconsistency");

    HaloGal[GalTree[q].HaloGalIndex].DisruptOn = 1;
#endif
    /* Put gas component to the central galaxy hot gas and stellar material into the ICM.
     * Note that the satellite should have no extended components. */

    transfer_gas(centralgal,HotGasComponent,p,ColdGasComponent,1.);
    transfer_gas(centralgal,HotGasComponent,p,HotGasComponent,1.);
#ifdef TRACK_BURST
    /* Transfer burst component first */
    transfer_stars(centralgal,BurstComponent,p,BurstComponent,
		   (Gal[p].DiskMass+Gal[p].BulgeMass)/(Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ICM));
#endif
    transfer_stars(centralgal,ICMComponent,p,DiskComponent,1.);
    transfer_stars(centralgal,ICMComponent,p,BulgeComponent,1.);
    /* Add satellite's luminosity into the luminosity of the ICL
     * component of the central galaxy. */

#ifndef POST_PROCESS_MAGS
#ifdef ICL
    int outputbin, j;
    for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(j = 0; j < NMAG; j++)
      {
#ifdef OUTPUT_REST_MAGS 
      	Gal[centralgal].ICLLum[j][outputbin] += Gal[p].Lum[j][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
      	Gal[centralgal].ObsICL[j][outputbin] += Gal[p].ObsLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
      	Gal[centralgal].dObsICL[j][outputbin] += Gal[p].dObsLum[j][outputbin];
#endif
#endif
      }  
    }
#endif //ICL
#endif //POST_PROCESS_MAGS

  } //if (rho_cen > rho_sat)
  mass_checks("Bottom of disrupt",centralgal);
  mass_checks("Bottom of disrupt",p);
}
