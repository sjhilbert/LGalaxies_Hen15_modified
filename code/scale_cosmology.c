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
 *  Created in: 2010
 *      Author: Qi Guo
 */

#include <math.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

#ifdef MCMC
#include "mcmc_vars.h"
#include "mcmc_proto.h"
#endif


/* When we have tables for the scaling parameter in any cosmology only the
 * cosmological parameters will be inputs. Then this function will read the
 * scaling parameters from the tables. */

 void read_scaling_parameters()
 {
  double om_min=0.1, om_max=0.6, s8_min=0.5, s8_max=1.0;
  double om_binsize, s8_binsize;
  int om_Nbin, s8_Nbin;
  double dummy_growth, dummy_snap63;
  char buf[1000], buf1[1000];
  FILE *fd;

  om_binsize=(om_max-om_min)/50;
  s8_binsize=(s8_max-s8_min)/50;
//0.369798725 0.589

  Omega=0.370;
  om_Nbin=(int)((0.370-om_min)/om_binsize);
  s8_Nbin=(int)((0.589-s8_min)/s8_binsize);

  sprintf(FileWithZList, "/galformod/scratch/bmh20/Workspace/CosmologyTables/zlist_%04d_%04d.txt", om_Nbin, s8_Nbin);

  read_zlist_new();
  read_output_snaps();

  sprintf(buf, "/galformod/scratch/bmh20/Workspace/CosmologyTables/fit_%04d_%04d.txt", om_Nbin, s8_Nbin);
  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", buf);
      terminate(sbuf);
    }

  fgets(buf1, 300, fd);
  fgets(buf1, 300, fd);
  fgets(buf1, 300, fd);

  if(fscanf(fd, "%lf %lf %lf %lf", &ScaleMass, &dummy_growth, &ScalePos, &dummy_snap63)!=4)
    {
      char sbuf[1000];
      sprintf(sbuf, "Wrong format of values in %s.\n", buf);
      terminate(sbuf);
    }

  fclose(fd);

  ScaleMass=1./ScaleMass;
  ScalePos=1./ScalePos;

  PartMass = 	PartMass_OriginalCosm * ScaleMass;
  BoxSize  =  BoxSize_OriginalCosm * ScalePos;

  printf("Boxsize=%f\n",BoxSize);

 }


 /** @file scale_cosmology.c
  *  @brief functions used to scale to a different cosmology
  *
  *  Add by Qi Guo adapted by Bruno Henriques */


void scale_cosmology(int nhalos)
{
  int i, j;
  double Scale_V,CenVel[3],dv;


  //Save unscaled properties
  for(i = 0; i < nhalos ; i++)
    {
      //will make sure haloes in the future are not scaled/un_scaled
      if(Halo[i].SnapNum<=LastSnapShotNr)
	{
	  HaloAux[i].M_Crit200_Unscaled = Halo[i].M_Crit200;
	  HaloAux[i].M_Mean200_Unscaled = Halo[i].M_Mean200;
	  HaloAux[i].Vmax_Unscaled = Halo[i].Vmax;
	  for (j = 0; j < 3 ; j++)
	    {
	      HaloAux[i].Pos_Unscaled[j] = Halo[i].Pos[j];
	      HaloAux[i].Vel_Unscaled[j] = Halo[i].Vel[j];
	      HaloAux[i].Spin_Unscaled[j] = Halo[i].Spin[j];
	    }
	}
    }

  for (i = 0; i < nhalos ; i++)
    {
      Scale_V = scale_v_cen(Halo[Halo[i].FirstHaloInFOFgroup].SnapNum);

      //will make sure haloes in the future are not scaled/un_scaled
      if(Halo[i].SnapNum<=LastSnapShotNr)
	{
	  if(Halo[i].M_Crit200 > 1.e-8)
	    Halo[i].M_Crit200 = Halo[i].M_Crit200 * ScaleMass * c_correction(Halo[i].M_Crit200,Halo[i].SnapNum);
	  if(Halo[i].M_Mean200 > 1.e-8)
	    Halo[i].M_Mean200 = Halo[i].M_Mean200 * ScaleMass * c_correction(Halo[i].M_Mean200,Halo[i].SnapNum);
	  Halo[i].Vmax = Halo[i].Vmax * sqrt(ScaleMass/ScalePos) * sqrt(AA_OriginalCosm[Halo[i].SnapNum]/AA[Halo[i].SnapNum]);

	  for (j = 0; j < 3 ; j++)
	    {
	      Halo[i].Pos[j] = Halo[i].Pos[j] * ScalePos;
	      Halo[i].Spin[j] *= ScalePos * sqrt(ScaleMass/ScalePos) * sqrt(AA[Halo[i].SnapNum]/AA_OriginalCosm[Halo[i].SnapNum]);

	      CenVel[j] = Halo[Halo[i].FirstHaloInFOFgroup].Vel[j] * Scale_V ;
	      if(i !=  Halo[i].FirstHaloInFOFgroup) // subhalos
		{
		  dv = Halo[i].Vel[j] - Halo[Halo[i].FirstHaloInFOFgroup].Vel[j];
		  dv *=sqrt(ScaleMass/ScalePos) * sqrt(AA_OriginalCosm[Halo[i].SnapNum]/AA[Halo[i].SnapNum]);
		  Halo[i].Vel[j] = CenVel[j] + dv;
		}
	      else //central halos
		Halo[i].Vel[j] = Halo[i].Vel[j] * Scale_V  ;
	    }
	}
    }
}

void un_scale_cosmology(int nhalos)
{
  int i, j;

  for(i = 0; i < nhalos ; i++)
    {
      //will make sure haloes in the future are not scaled/un_scaled
      if(Halo[i].SnapNum<=LastSnapShotNr)
	{
	  Halo[i].M_Crit200 = HaloAux[i].M_Crit200_Unscaled;
	  Halo[i].M_Mean200 = HaloAux[i].M_Mean200_Unscaled;
	  Halo[i].Vmax = HaloAux[i].Vmax_Unscaled;

	  for (j = 0; j < 3 ; j++)
	    {
	      Halo[i].Pos[j] = HaloAux[i].Pos_Unscaled[j];
	      Halo[i].Vel[j] = HaloAux[i].Vel_Unscaled[j];
	      Halo[i].Spin[j] = HaloAux[i].Spin_Unscaled[j];
	    }
	}
    }
}



/** @brief computes c correction */
double c_correction(const float mass, const int snapnum)
{
  double c_original, c_new, Omega_new, Omega_original, ratio;

  c_original = 5 * pow(0.0001 * mass, -0.1);
  
  Omega_new = Omega * 1./pow3(AA[snapnum]) / 
             (Omega * 1./pow3(AA[snapnum]) + OmegaLambda);
  
  Omega_original = Omega_OriginalCosm * 1./pow3(AA_OriginalCosm[snapnum]) /
                  (Omega_OriginalCosm * 1./pow3(AA_OriginalCosm[snapnum]) + OmegaLambda_OriginalCosm);
                  
  ratio = Omega_original/ Omega_new;
  
  c_new = find_c(c_original, ratio);

  return func_c(c_new) / func_c(c_original);
}


/** @brief finds c
 * 
 * finds c using bisection
 * 
 * since the original version of this function showed up
 * surprisingly high on profile, a more optimized bisection
 * version was implemented
 * 
 * @warning assumes that initial values bracket the result 
 */
double find_c(const double c_ori_, const double ratio_)
{
  const double constant_ = ratio_ * func_c_p(c_ori_);
  
  double x_1_ = 1.;
  double x_2_ = 50.;
  
  double f_1_ = func_c_p(x_1_) - constant_;
  double x_m_, f_m_;
  while(
    x_m_ = 0.5 * (x_1_ + x_2_),
    f_m_ = func_c_p(x_m_) - constant_,
    fabs(f_m_) > 1.e-8
  ) 
  {
    if(f_m_ * f_1_ > 0)
    {       
      x_1_ = x_m_;
      f_1_ = f_m_;
    }
    else
    { x_2_ = x_m_; }
  }
  return x_m_;
}


double func_c(const double c)
{
  return log(1 + c) - c / (1 + c);
}


double func_c_p(const double c)
{
  return (log(1 + c) - c / (1 + c)) / (c * c * c);
}


double dgrowth_factor_dt(double a, double omega_m, double omega_l)
{
  double o_m, o_l, do_m, do_l,  hubble_a, g0, g, extra_fac, dg, dDdt;
  o_m = omega_m * 1./pow3(a)/(omega_m * 1./pow3(a) + omega_l);
  o_l = omega_l / (omega_m * 1./pow3(a) + omega_l);
  do_m = -3 * omega_m * omega_l / (a * a * a * a) / (omega_m / (a * a * a) + omega_l) / (omega_m / (a * a * a) + omega_l);
  do_l = -do_m;

  hubble_a = sqrt(omega_m/(a*a*a) + omega_l);

  g0 = 2.5 * omega_m / (pow(omega_m, 4.0/7.0) - omega_l + (1.0 + 0.5*omega_m)*(1.0 + omega_l/70.0));
  g = 2.5 * o_m / (pow(o_m, 4.0/7.0) - o_l + (1.0 + 0.5*o_m)*(1.0 + o_l/70.0));

  //da_dtau = sqrt(1 + o_m * (1 / a - 1) + o_l * (a * a - 1));   //tau = H0*t
  extra_fac = - ( 4/7.* pow(o_m, -3./7) * do_m - do_l
		  +(do_m /2. *(1 + o_l / 70.) - (1 + o_m / 2.) * do_l / 70)
		    /(1 + o_l/70.)/(1 + o_l / 70))/(pow(o_m, 4.0/7.0) - o_l + (1.0 + 0.5*o_m)*(1.0 + o_l/70.0))/(pow(o_m, 4.0/7.0) - o_l + (1.0 + 0.5*o_m)*(1.0 + o_l/70.0));

  dg = 2.5 * do_m / (pow(o_m, 4.0/7.0) - o_l + (1.0 + 0.5*o_m)*(1.0 + o_l/70.0)) + 2.5 * o_m * extra_fac;

  dDdt = a * hubble_a * (dg * a +  g )/ g0;

  return dDdt;

}


 double scale_v_cen(int snapnum)
{
  double Scale_V;

  Scale_V= ScalePos * dgrowth_factor_dt(AA[snapnum],Omega, OmegaLambda) /
        		dgrowth_factor_dt(AA_OriginalCosm[snapnum],Omega_OriginalCosm,OmegaLambda_OriginalCosm) *
  	       	AA[snapnum]/AA_OriginalCosm[snapnum] *  Hubble_h/Hubble_h_OriginalCosm;
  return Scale_V;

}

