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
 *  You should have received a_ copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

/**@file    model_infall.c
 *  @date   2016-2019
 *  @author ?
 *  @author Stefan Hilbert
 *
 * @brief   calculates the amount of gas that infalls
 *          into the galaxy hot gas component at each time step. This
 *          is derived from the baryonic fraction taking reionization
 *          into account.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

double infall_recipe(const int central_galaxy_number_, const int n_galaxies_, double current_redshift_)
{
  int galaxy_number_;
  double total_mass_, reionization_modifier_, infalling_mass_;
  double distance_;

  /*  need to add up all the baryonic mass asociated with the full halo to check
   *  what baryonic fraction is missing/in excess. That will give the mass of gas
   *  that need to be added/subtracted to the hot phase, gas that infalled.*/
  total_mass_ = 0.0;

  for(galaxy_number_ = 0; galaxy_number_ < n_galaxies_; galaxy_number_++) 
  {    /* Loop over all galaxies in the FoF-halo */
    /* distance_ is the separation of the galaxy which galaxy_number_ orbits from the type 0 */
    distance_=separation_gal(central_galaxy_number_,Gal[galaxy_number_].CentralGal)/(1+ZZ[Halo[Gal[central_galaxy_number_].HaloNr].SnapNum]);

    /* If galaxy is orbiting a_ galaxy inside Rvir of the type 0 it will contribute
     * to the baryon sum */
    if ( distance_ < Gal[central_galaxy_number_].Rvir )
    {
      total_mass_ += Gal[galaxy_number_].DiskMass + Gal[galaxy_number_].BulgeMass + Gal[galaxy_number_].ICM + Gal[galaxy_number_].BlackHoleMass;
      total_mass_ += Gal[galaxy_number_].ColdGas + Gal[galaxy_number_].HotGas + Gal[galaxy_number_].EjectedMass;
    }
  }

  /* The logic here seems to be that the infalling mass is supposed to be
   * newly-accreted diffuse gas.  So we take the total mass and subtract all
   * the components that we already know about.  In principle, this could lead
   * to negative infall. */
  /* The reionization modifier_ is applied to the whole baryon mass, not just the 
   * diffuse component.  It is not obvious that this is the correct thing to do. */
  /* The baryonic fraction is conserved by adding/subtracting the infalling_mass_
   * calculated here to/from the hot gas of the central galaxy of the FOF
   * This is done in sam.c where the infall recipe is called.
   * If ReionizationModel<2, the impact of reonization on the fraction of infalling
   * gas is computed, this is done using the Gnedin formalism with a_ choice
   * of fitting parameters to the formulas proposed by these authors.
   * There are two options. In both cases reionization has the effect of reducing 
   * the fraction of baryons that collapse into dark matter halos, reducing the
   * amount of infalling gas. */
  if(ReionizationModel == 2)
    reionization_modifier_ = 1.0;
  else
    reionization_modifier_ = get_reionization_modifier(Gal[central_galaxy_number_].Mvir, current_redshift_);

  infalling_mass_ = reionization_modifier_ * BaryonFrac * Gal[central_galaxy_number_].Mvir - total_mass_;

  return infalling_mass_;
}


/** @brief computes reionization modifier */
double get_reionization_modifier(const float M_vir_, const double current_redshift_)
{
  double modifier_ = 1.;
  if (ReionizationModel == 2) 
  {
    printf("Should not be called with this option\n");
    exit(0);
  }
  else if (ReionizationModel == 0)
  {
    /* reionization recipie described in Gnedin (2000), with the fitting
     * from Okamoto et al. 2008 -> Qi(2010)*/

    //  const double a_ = 1. / (1 + current_redshift_);
    //  // a0 = 1.;
    //  
    //  // /* if not use Okamoto et al. 2008*/
    //  const double x       = - (1 - Omega) * a_ * a_ * a_ / (Omega + (1-Omega) * a_ * a_ * a_);
    //  const double delta_c = (178 + 82 *x - 39 * x * x) / (1. + x);
    //
    //  // const double x0 = - (1 - Omega) * a0 * a0 * a0 / (Omega + (1-Omega) * a0 * a0 * a0);
    //  const double x0      = - (1 - Omega);
    //  const double delta_0 = (178 + 82 *x0 - 39 * x0 * x0) / (1. + x0);
    //
    //  const double tau = 0.73 * pow(current_redshift_ + 1,0.18)* exp(-pow(0.25 * current_redshift_, 2.1));
    //  
    //  const double M_c_ = sqrt(tau * tau * tau * a_ * a_ * a_ * delta_0 / delta_c);

    /* if use Okamoto et al. 2008*/
    int tab_index_; 
    double M_c_;
    linear_interpolate(tab_index_, 0, N_REION_Z, current_redshift_, Reion_z, Reion_log10_Mc, M_c_, >, linear); 
    M_c_ = pow(10, M_c_-10);
 
   // const double alpha_ = 2.0;
   // modifier_ = pow(1 + (pow(2, alpha_/3.) -1) * pow(M_c_ / M_vir_, alpha_), -3./alpha_);
   // 0.5874010519681995 = (pow(2, alpha_/3.) -1)  for alpha_ = 2.0
    modifier_ = pow(1 + 0.5874010519681995 * M_c_ * M_c_ / (M_vir_ * M_vir_), -1.5);
  }
  else if (ReionizationModel ==1)
  {
    /** reionization recipie described in Gnedin (2000), using the fitting */
    /*  formulas given by Kravtsov et al (2004) Appendix B, used after Delucia 2004*/

    /*  here are two parameters that Kravtsov et al keep fixed. */
    /*  alpha_ gives the best fit to the Gnedin data */
    
    //Gnedin (2000)
    
    const double alpha_ = 6.0;
    // const double Tvir  = 1e4;

    /*  calculate the filtering mass */

    const double a_ = 1.0 / (1.0 + current_redshift_);
    const double a0_on_a_ = a0 / a_;
    const double ar_on_a_ = ar / a_;

    // double f_of_a_;
    // if(a_ <= a0)
    //   f_of_a_ = 3.0 * a_ / ((2.0 * alpha_) * (5.0 + 2.0 * alpha_)) * pow(a_on_a0, alpha_);
    // else if((a_ > a0) && (a_ < ar))
    //   f_of_a_ =
    //     (3.0 / a_) * a0 * a0 * (1.0 / (2.0 + alpha_) - 2.0 * sqrt(a0_on_a_) / (5.0 + 2.0 * alpha_)) +
    //     a_ * a_ / 10.0 - (a0 * a0 / 10.0) * (5.0 - 4.0 * sqrt(a0_on_a_));
    // else
    //   f_of_a_ = (3.0 / a_) * (a0 * a0 * (1.0 / (2.0 + alpha_) - 2.0 * sqrt(a0_on_a_) / (5.0 + 2.0 * alpha_)) +
    //         (ar * ar / 10.0) * (5.0 - 4.0 * sqrt(ar_on_a_)) - (a0 * a0 / 10.0) * (5.0 -
    //             4.0 * sqrt(a0_on_a_)) + a_ * ar / 3.0 - (ar * ar / 3.0) * (3.0 - 2.0 * sqrt(ar_on_a_)));

    double f_of_a_;
    if(a_ <= a0)
      f_of_a_ = (3.0 / ((2.0 * alpha_) * (5.0 + 2.0 * alpha_))) * a_ * pow(a0_on_a_, -alpha_);
    else if((a_ > a0) && (a_ < ar))
      f_of_a_ = a0 * a0_on_a_ * (3.0 / (2.0 + alpha_) - (6.0 / (5.0 + 2.0 * alpha_)) * sqrt(a0_on_a_)) + 0.1 * a_ * a_ - 0.1 * a0 * a0 * (5.0 - 4.0 * sqrt(a0_on_a_));
    else
      f_of_a_ = (a0 * a0_on_a_ * (3.0 / (2.0 + alpha_) - (6.0 / (5.0 + 2.0 * alpha_)) * sqrt(a0_on_a_)) +
            0.3 * ar * ar_on_a_ * (5.0 - 4.0 * sqrt(ar_on_a_)) - 0.3 * a0 * a0_on_a_ * (5.0 - 4.0 * sqrt(a0_on_a_)) + ar  -  ar * ar_on_a_ * (3.0 - 2.0 * sqrt(ar_on_a_)));

    /*  this is in units of 10^10Msun/h, note mu=0.59 and mu^-1.5 = 2.21 */
    // Mjeans = 25.0 * pow(Omega, -0.5) * 2.21;
    // M_filtering_ = Mjeans * pow(f_of_a_, 1.5);
    const double M_filtering_ =  (25.0 * 2.21) * sqrt(f_of_a_ * f_of_a_ * f_of_a_ / Omega);

    /*  calculate the characteristic mass coresponding to a_ halo temperature of 10^4K */
    // V_char_ = sqrt(Tvir / 36.0);
    const double V_char_ = 100. / 6.;
    // const double omegaZ = Omega * (pow3(1.0 + current_redshift_) / (Omega * pow3(1.0 + current_redshift_) + OmegaLambda));
    // const double x_of_z_ = omegaZ - 1.0;
    const double x_of_z_ = Omega / (Omega  + OmegaLambda * pow3(a_)) - 1.0;
    const double deltacrit_of_z_ = 18.0 * M_PI * M_PI + 82.0 * x_of_z_ - 39.0 * x_of_z_ * x_of_z_;
    // const double HubbleZ = Hubble * sqrt(Omega * pow3(1.0 + current_redshift_) + OmegaLambda);

    // const double M_char_ = V_char_ * V_char_ * V_char_ / (Gravity * HubbleZ * sqrt(0.5 * deltacrit_of_z_));
    const double M_char_ = (V_char_ * V_char_ * V_char_) / (Gravity * Hubble * sqrt((Omega * pow3(1.0 + current_redshift_) + OmegaLambda) * 0.5 * deltacrit_of_z_));

    /*  we use the maximum of M_filtering_ and M_char_ */
    // const double mass_to_use = max(M_filtering_, M_char_);
    // modifier_ = 1.0 / pow3(1.0 + 0.26 * (mass_to_use / M_vir_));
    modifier_ = 1.0 / pow3(1.0 + 0.26 * (max(M_filtering_, M_char_) / M_vir_));
  }
  return modifier_;
}


/** @brief adds infalling gas to the hot gas of the central galaxy.  */
void add_infall_to_hot(const int central_galaxy_number_, const double infalling_gas_) 
{
  /*  Add the infalling gas to the central galaxy hot component */
  Gal[central_galaxy_number_].HotGas += infalling_gas_;

#ifdef INDIVIDUAL_ELEMENTS
  //Gal[central_galaxy_number_].HotGas_elements.H += 0.75 * infalling_gas_ * 1.0e10;
  Gal[central_galaxy_number_].HotGas_elements.H  += 0.75 * 1.0e10 * (infalling_gas_ * inv_Hubble_h);
  Gal[central_galaxy_number_].HotGas_elements.He += 0.25 * 1.0e10 * (infalling_gas_ * inv_Hubble_h);
#endif
}
