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

/** @file  model_reincorporation.c
 *  @brief model_reincorporation.c calculates the fraction of ejected gas that
 *         gets reincorporated into the hot fraction per timestep.
 *
 *         23options are available to reincorporate the gas from the external
 *         reservoir:
 *           -\f$\dot{M}_{\rm{eject}}=
 *                -\gamma \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *         (Eq. 3 Delucia2004) (ReIncorporationModel == 2);
 *           -\f$\dot{M}_{\rm{eject}}=
 *                -\gamma \left(\frac{V_{\rm{vir}}}{\rm{220km/s}}\right)
 *                \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *         (Eq. 23 Guo2010) (ReIncorporationModel == 1)
 *
 **/
/** @brief reincorporates ejected gas back into the central galaxy hot halo */

void reincorporate_gas(int p, double dt)
{
  double reincorporated, fraction, reinc_time;

  reincorporated = 0.;

  mass_checks("reincorporate_gas #1",p);


  if(FeedbackEjectionModel == 0)
    {
      if(ReIncorporationModel == 0)
	{
	  reinc_time= (Hubble_h/Gal[p].Mvir)*(ReIncorporationFactor/UnitTime_in_years);
	  reincorporated = Gal[p].EjectedMass / reinc_time * dt;
	  /* Henriques2013 Mdot_eject=-gama_ej*M_ejected*M_vir Mvir should be in units of 1e12, but inside the
	   * code Mvir is already in units of 1.e10*/
	}
      else
	if(ReIncorporationModel == 1)
	  reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * Gal[p].Vvir/220. *dt ;
      /* Guo2010 -> Mdot_eject=-gama_ej * M_ejected/tdyn * Vvir/220 */
	else
	  if(ReIncorporationModel == 2)
	    reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * dt;
    }
  else if(FeedbackEjectionModel == 1)
    {
      reincorporated = ReIncorporationFactor * Gal[p].EjectedMass /
	  (Gal[p].Rvir * min(FeedbackEjectionEfficiency,1.)*sqrt(EtaSNcode * EnergySNcode)/(Gal[p].Vvir*Gal[p].Vvir))
	  * Gal[p].Vvir/220. * 1.e-6* dt ;
    }

  if (reincorporated > Gal[p].EjectedMass)
    reincorporated = Gal[p].EjectedMass;
	
  mass_checks("reincorporate_gas #1.5",p);

  /*Update ejected and hot gas contents*/
  if (Gal[p].EjectedMass > 0.)
    {
      fraction=((float)reincorporated)/Gal[p].EjectedMass;
      transfer_gas(p,HotGasComponent,p,EjectedGasComponent,fraction);
    }

  mass_checks("reincorporate_gas #2",p);

}

