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
 
/** @file   model_reincorporation.c
 *  @date   2016-2019
 *  @author ?
 *  @author Stefan Hilbert
 *
 *  @brief  model_reincorporation.c calculates the fraction of ejected gas that
 *          gets reincorporated into the hot fraction per timestep.
 *          
 *          options are available to reincorporate the gas from the external
 *          reservoir:
 *            -\f$\dot{M}_{\rm{eject}}=
 *                 -\gamma \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *          (Eq. 3 Delucia2004) (ReIncorporationModel == 2);
 *            -\f$\dot{M}_{\rm{eject}}=
 *                 -\gamma \left(\frac{V_{\rm{vir}}}{\rm{220km/s}}\right)
 *                 \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *          (Eq. 23 Guo2010) (ReIncorporationModel == 1)
 **/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @brief reincorporates ejected gas back into the central galaxy hot halo */
void reincorporate_gas(const int galaxy_number_, const double dt_)
{
  double reincorporated = 0.;

  mass_checks("reincorporate_gas #1",galaxy_number_);

  if(FeedbackEjectionModel == 0)
  {
    if(ReIncorporationModel == 0)
    {
      // const double reinc_time     = (Hubble_h * ReIncorporationFactor) / ( Gal[galaxy_number_].Mvir * UnitTime_in_years);
     // reincorporated = Gal[galaxy_number_].EjectedMass / reinc_time * dt_;
      reincorporated = Gal[galaxy_number_].EjectedMass * Gal[galaxy_number_].Mvir * UnitTime_in_years * dt_ / (Hubble_h * ReIncorporationFactor);
      /* Henriques2013 Mdot_eject=-gama_ej*M_ejected*M_vir Mvir should be in units of 1e12, but inside the
        * code Mvir is already in units of 1.e10*/
    }
    else if(ReIncorporationModel == 1)
    { reincorporated = ReIncorporationFactor * Gal[galaxy_number_].EjectedMass / (Gal[galaxy_number_].Rvir / Gal[galaxy_number_].Vvir) * Gal[galaxy_number_].Vvir/220. *dt_; }
    /* Guo2010 -> Mdot_eject=-gama_ej * M_ejected/tdyn * Vvir/220 */
    else if(ReIncorporationModel == 2)
    { reincorporated = ReIncorporationFactor * Gal[galaxy_number_].EjectedMass / (Gal[galaxy_number_].Rvir / Gal[galaxy_number_].Vvir) * dt_; }
  }
  else if(FeedbackEjectionModel == 1)
  {
    reincorporated = ReIncorporationFactor * Gal[galaxy_number_].EjectedMass /
        (Gal[galaxy_number_].Rvir * min(FeedbackEjectionEfficiency,1.)*sqrt(EtaSNcode * EnergySNcode)/(Gal[galaxy_number_].Vvir*Gal[galaxy_number_].Vvir))
        * Gal[galaxy_number_].Vvir/220. * 1.e-6* dt_ ;
  }

  if (reincorporated > Gal[galaxy_number_].EjectedMass)
  { reincorporated = Gal[galaxy_number_].EjectedMass; }
        
  mass_checks("reincorporate_gas #1.5",galaxy_number_);

  /*Update ejected and hot gas contents*/
  if (Gal[galaxy_number_].EjectedMass > 0.)
  {
    const double fraction=((float)reincorporated)/Gal[galaxy_number_].EjectedMass;
    transfer_gas(galaxy_number_,HotGasComponent,galaxy_number_,EjectedGasComponent,fraction);
  }

  mass_checks("reincorporate_gas #2",galaxy_number_);
}

