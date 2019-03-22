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


/**@file model_infall.c
 * @brief model_infall.c calculates the amount of gas that infalls
 *       into the galaxy hot gas component at each time step. This
 *       is derived from the baryonic fraction taking reionization
 *       into account.
 */

double infall_recipe(int centralgal, int ngal, double Zcurr)
{
  int i;
  double tot_mass, reionization_modifier, infallingMass;
  double dis;

  /*  need to add up all the baryonic mass asociated with the full halo to check
   *  what baryonic fraction is missing/in excess. That will give the mass of gas
   *  that need to be added/subtracted to the hot phase, gas that infalled.*/
  tot_mass = 0.0;

  for(i = 0; i < ngal; i++) 
  {    /* Loop over all galaxies in the FoF-halo */
    /* dis is the separation of the galaxy which i orbits from the type 0 */
    dis=separation_gal(centralgal,Gal[i].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

    /* If galaxy is orbiting a galaxy inside Rvir of the type 0 it will contribute
     * to the baryon sum */
    if ( dis < Gal[centralgal].Rvir ) {
      tot_mass += Gal[i].DiskMass + Gal[i].BulgeMass + Gal[i].ICM + Gal[i].BlackHoleMass;
      tot_mass += Gal[i].ColdGas + Gal[i].HotGas + Gal[i].EjectedMass;
    }
  }

  /* The logic here seems to be that the infalling mass is supposed to be
   * newly-accreted diffuse gas.  So we take the total mass and subtract all
   * the components that we already know about.  In principle, this could lead
   * to negative infall. */
  /* The reionization modifier is applied to the whole baryon mass, not just the 
   * diffuse component.  It is not obvious that this is the correct thing to do. */
  /* The baryonic fraction is conserved by adding/subtracting the infallingMass
   * calculated here to/from the hot gas of the central galaxy of the FOF
   * This is done in main.c where the infall recipe is called.
   * If ReionizationModel<2, the impact of reonization on the fraction of infalling
   * gas is computed, this is done using the Gnedin formalism with a choice
   * of fitting parameters to the formulas proposed by these authors.
   * There are two options. In both cases reionization has the effect of reducing 
   * the fraction of baryons that collapse into dark matter halos, reducing the
   * amount of infalling gas. */
  if(ReionizationModel == 2)
    reionization_modifier = 1.0;
  else
    reionization_modifier = get_reionization_modifier(Gal[centralgal].Mvir, Zcurr);

  infallingMass = reionization_modifier * BaryonFrac * Gal[centralgal].Mvir - tot_mass;

  return infallingMass;
}


/** @brief computes reionization modifier  */
double get_reionization_modifier(float Mvir, double Zcurr)
{
  double modifier=1.;
  if (ReionizationModel == 2) 
  {
    printf("Should not be called with this option\n");
    exit(0);
  }
  else if (ReionizationModel == 0)
  {
    /* reionization recipie described in Gnedin (2000), with the fitting
     * from Okamoto et al. 2008 -> Qi(2010)*/

    //  const double a = 1. / (1 + Zcurr);
    //  // a0 = 1.;
    //  
    //  // /* if not use Okamoto et al. 2008*/
    //  const double x       = - (1 - Omega) * a * a * a / (Omega + (1-Omega) * a * a * a);
    //  const double delta_c = (178 + 82 *x - 39 * x * x) / (1. + x);
    //
    //  // const double x0 = - (1 - Omega) * a0 * a0 * a0 / (Omega + (1-Omega) * a0 * a0 * a0);
    //  const double x0      = - (1 - Omega);
    //  const double delta_0 = (178 + 82 *x0 - 39 * x0 * x0) / (1. + x0);
    //
    //  const double tau = 0.73 * pow(Zcurr + 1,0.18)* exp(-pow(0.25 * Zcurr, 2.1));
    //  
    //  const double Mc = sqrt(tau * tau * tau * a * a * a * delta_0 / delta_c);

    /* if use Okamoto et al. 2008*/
    int tabindex; 
    double Mc;
    linear_interpolate(tabindex, 0, N_REION_Z, Zcurr, Reion_z, Reion_log10_Mc, Mc, >, linear); 
    Mc = pow(10, Mc-10);
 
   // const double alpha = 2.0;
   // modifier = pow(1 + (pow(2, alpha/3.) -1) * pow(Mc / Mvir, alpha), -3./alpha);
   // 0.5874010519681995 = (pow(2, alpha/3.) -1)  for alpha = 2.0
    modifier = pow(1 + 0.5874010519681995 * Mc * Mc / (Mvir * Mvir), -1.5);
  }
  else if (ReionizationModel ==1)
  {
    /** reionization recipie described in Gnedin (2000), using the fitting */
    /*  formulas given by Kravtsov et al (2004) Appendix B, used after Delucia 2004*/

    /*  here are two parameters that Kravtsov et al keep fixed. */
    /*  alpha gives the best fit to the Gnedin data */
    
    //Gnedin (2000)
    double a, alpha;
    double f_of_a, a_on_a0, a_on_ar, Mfiltering, Mjeans, Mchar, mass_to_use;
    double Tvir, Vchar, omegaZ, xZ, deltacritZ, HubbleZ;
    
    alpha = 6.0;
    Tvir = 1e4;

    /*  calculate the filtering mass */

    a = 1.0 / (1.0 + Zcurr);
    a_on_a0 = a / a0;
    a_on_ar = a / ar;

    if(a <= a0)
      f_of_a = 3.0 * a / ((2.0 * alpha) * (5.0 + 2.0 * alpha)) * pow(a_on_a0, alpha);
    else if((a > a0) && (a < ar))
      f_of_a =
        (3.0 / a) * a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
        a * a / 10.0 - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5));
    else
      f_of_a = (3.0 / a) * (a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
            (ar * ar / 10.0) * (5.0 - 4.0 * pow(a_on_ar, -0.5)) - (a0 * a0 / 10.0) * (5.0 -
                4.0 * pow(a_on_a0, -0.5)) + a * ar / 3.0 - (ar * ar / 3.0) * (3.0 - 2.0 * pow(a_on_ar, -0.5)));

    /*  this is in units of 10^10Msun/h, note mu=0.59 and mu^-1.5 = 2.21 */
    Mjeans = 25.0 * pow(Omega, -0.5) * 2.21;
    Mfiltering = Mjeans * pow(f_of_a, 1.5);


    /*  calculate the characteristic mass coresponding to a halo temperature of 10^4K */
    Vchar = sqrt(Tvir / 36.0);
    omegaZ = Omega * (pow3(1.0 + Zcurr) / (Omega * pow3(1.0 + Zcurr) + OmegaLambda));
    xZ = omegaZ - 1.0;
    deltacritZ = 18.0 * M_PI * M_PI + 82.0 * xZ - 39.0 * xZ * xZ;
    HubbleZ = Hubble * sqrt(Omega * pow3(1.0 + Zcurr) + OmegaLambda);

    Mchar = Vchar * Vchar * Vchar / (G * HubbleZ * sqrt(0.5 * deltacritZ));

    /*  we use the maximum of Mfiltering and Mchar */
    mass_to_use = max(Mfiltering, Mchar);
    modifier = 1.0 / pow3(1.0 + 0.26 * (mass_to_use / Mvir));
  }
  return modifier;
}


/** @brief adds infalling gas to the hot gas of the central galaxy.  */
void add_infall_to_hot(int centralgal, double infallingGas) {

  /*  Add the infalling gas to the central galaxy hot component */
  Gal[centralgal].HotGas += infallingGas;

#ifdef INDIVIDUAL_ELEMENTS
  //Gal[centralgal].HotGas_elements.H += 0.75*infallingGas*1.0e10;
  Gal[centralgal].HotGas_elements.H += 0.75*(infallingGas/Hubble_h)*1.0e10;
  Gal[centralgal].HotGas_elements.He += 0.25*(infallingGas/Hubble_h)*1.0e10;
#endif
}
