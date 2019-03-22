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
get_peri_radius_for_galaxy(const int galaxy_number_, const int central_galaxy_number_)
{
  int i;
  double a, b, v[3], r[3], x, x0;
  for(i = 0; i < 3; i++)
    {
      r[i] = wrap(Gal[galaxy_number_].Pos[i]-Gal[central_galaxy_number_].Pos[i],BoxSize);
      r[i] /= (1 + ZZ[Halo[Gal[central_galaxy_number_].HaloNr].SnapNum]);
      v[i] = Gal[galaxy_number_].Vel[i] - Gal[central_galaxy_number_].Vel[i];
    }

  b = 1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / pow2(Gal[central_galaxy_number_].Vvir);
  a = 1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] -
                pow2(r[0] * v[0] + r[1] * v[1] + r[2] * v[2])/
                    (r[0] * r[0] + r[1] * r[1] + r[2] * r[2])) / pow2(Gal[central_galaxy_number_].Vvir);

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
get_sat_radius_for_galaxy(const int galaxy_number_)
{
  double r, rd, rb, Mdisk, rmax, M;
  double Mgas, Mbulge, rgd, dr, totmass;
  int ii;
  #define SAT_RADIUS_RMIN 5e-7
  #define SAT_RADIUS_N 100

  r=0.;
  rgd = Gal[galaxy_number_].GasDiskRadius/3.;
  rd=Gal[galaxy_number_].StellarDiskRadius/3.;
  rb=Gal[galaxy_number_].BulgeSize;
  Mgas = Gal[galaxy_number_].ColdGas;
  Mdisk=Gal[galaxy_number_].DiskMass;
  Mbulge = Gal[galaxy_number_].BulgeMass;
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
void disrupt(const int galaxy_number_)
{
  double rho_sat, rho_cen;
  double cen_mass, r_sat, radius;
  int central_galaxy_number_;

  /* If the main halo density at the pericentre (closest point in the orbit
   * to the central galaxy)is larger than the satellite's density at the
   * half mass radius, the satellite is completely disrupted. Note that since
   * the satellite is a type 2 the only mass components remaining and
   * contributing to the density are the cold gas and stellar mass. */

  central_galaxy_number_=Gal[galaxy_number_].CentralGal;
 
  mass_checks("Top of disrupt",central_galaxy_number_);
  mass_checks("Top of disrupt",galaxy_number_);

  /* Radius calculated at the peri-point */
  radius = get_peri_radius_for_galaxy(galaxy_number_, central_galaxy_number_);
  if (radius < 0) {
   terminate("must be wrong \n");
  }

  /* Calculate the density of the main central halo at radius (the peri-centre).
   * The tidal forces are caused by the dark matter of the main halo, hence Mvir
   * is used. Assume isothermal. */
  cen_mass=Gal[central_galaxy_number_].Mvir*radius/Gal[central_galaxy_number_].Rvir;
  rho_cen=cen_mass/pow3(radius);

  /* Calculate the density of the satellite's baryonic material */
  if (Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass>0)
  {
    /* Calculate the rho according to the real geometry */
    r_sat = get_sat_radius_for_galaxy(galaxy_number_);
    rho_sat=(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass+Gal[galaxy_number_].ColdGas)/pow3(r_sat);
  }
  else
    rho_sat=0.0;

  /* If density of the main halo is larger than that of the satellite baryonic
   * component, complete and instantaneous disruption is assumed. Galaxy becomes
   * a type 3 and all its material is transferred to the central galaxy. */
  if (rho_cen > rho_sat)
  {
    Gal[galaxy_number_].Type = 3;
#ifdef GALAXYTREE
    int q;
    q = Gal[Gal[galaxy_number_].CentralGal].FirstProgGal;
    if (q >= 0)
    {
      // add progenitors of Gal[galaxy_number_] to the list of progentitors of Gal[galaxy_number_].CentralGal
      while (GalTree[q].NextProgGal >= 0)
      	q = GalTree[q].NextProgGal;
	
      GalTree[q].NextProgGal = Gal[galaxy_number_].FirstProgGal;
      
      if(GalTree[q].NextProgGal >= NGalTree)
	{
	  printf("q=%d galaxy_number_=%d GalTree[q].NextProgGal=%d NGalTree=%d\n",
		 q, galaxy_number_, GalTree[q].NextProgGal, NGalTree);
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

    transfer_gas(central_galaxy_number_,HotGasComponent,galaxy_number_,ColdGasComponent,1.);
    transfer_gas(central_galaxy_number_,HotGasComponent,galaxy_number_,HotGasComponent,1.);
#ifdef TRACK_BURST
    /* Transfer burst component first */
    transfer_stars(central_galaxy_number_,BurstComponent,galaxy_number_,BurstComponent,
		   (Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass)/(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass+Gal[galaxy_number_].ICM));
#endif
    transfer_stars(central_galaxy_number_,ICMComponent,galaxy_number_,DiskComponent,1.);
    transfer_stars(central_galaxy_number_,ICMComponent,galaxy_number_,BulgeComponent,1.);
    
    /* Add satellite's luminosity into the luminosity of the ICL
     * component of the central galaxy. */
#ifndef POST_PROCESS_MAGS
#ifdef ICL
    int output_number_, filter_number_;
    for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    {
      for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
      {
#ifdef OUTPUT_REST_MAGS 
      	Gal[central_galaxy_number_].ICLLum [output_number_][filter_number_] += Gal[galaxy_number_].Lum    [output_number_][filter_number_];
#endif                                                                                           
#ifdef OUTPUT_OBS_MAGS                                                                          
      	Gal[central_galaxy_number_].ObsICL [output_number_][filter_number_] += Gal[galaxy_number_].ObsLum [output_number_][filter_number_];
#ifdef OUTPUT_FB_OBS_MAGS                                                                       
      	Gal[central_galaxy_number_].backward_ObsICL[output_number_][filter_number_] += Gal[galaxy_number_].backward_ObsLum[output_number_][filter_number_];
      	Gal[central_galaxy_number_].forward_ObsICL[output_number_][filter_number_] += Gal[galaxy_number_].forward_ObsLum[output_number_][filter_number_];
#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
      }  
    }
#endif /* defined ICL */
#endif /* not defined POST_PROCESS_MAGS */

  } //if (rho_cen > rho_sat)
  mass_checks("Bottom of disrupt",central_galaxy_number_);
  mass_checks("Bottom of disrupt",galaxy_number_);
}
