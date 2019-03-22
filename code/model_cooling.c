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

/** @file model_cooling.c
 *  @brief model_cooling.c calculates the amount of mass that cools
 *  from the hot to the cold phase at each timestep (this fraction
 *  is then reduced by AGN heating) and updates the hot and cold gas
 *  fractions accordingly.
 *
 * This recipe calculates the amount of mass that cools from the hot
 * to the cold phase at each timestep. Two different infalling regimes
 * are assumed depending on the redshift and mass of the halos. This is
 * the recipe proposed by (White & Rees, 1978) and that has since been
 * used in most SA.
 *
 * Before only central galaxies could cool gas, so R_hot was just set
 * equal to Rvir. After Guo2010, satellites can have cooling, but since
 * they are loosing dark matter, Rvir is not a good approximation of R_Hot
 * so Gal.HotRadius was introduced.
 *
 * -> At early times and for low mass halos, the cooling radius can be
 * larger then the virial radius. In this case, despite the infalling
 * gas being shock heated to the virial temperature, it condenses within
 * a halo dynamical time and a quasi-static atmosphere cannot form.
 * Single line on the code, with no calculations needed, just that
 * \f$t_{\rm{dyn,h}}=R_{\rm{vir}}/V_{\rm{vir}}\f$
 * and therefore
 * \f$\dot{M}_{\rm{cool}}=0.5\frac{m_{\rm{hot}}}{t_{\rm{dyn}}}
 * =0.5m_{\rm{hot}}\frac{V_{\rm{vir}}}{R_{\rm{vir}}}\f$
 *
 *
 * -> For massive halos and at late times, the cooling radius lies within
 * the virial radius and the infalling gas does form a quasi-static hot
 * atmosphere that extends to the virial radius. This gas can cool at
 * later times and its accretion into central regions is modeled through
 * a cooling flow.
 *
 * \f$\dot{M}_{\rm{cool}}=
 *  m_{\rm{hot}}\frac{r_{\rm{cool}}}{R_{\rm{vir}}}\frac{1}{t_{\rm{cool}}}\f$
 *  eq5 (no 0.5 factor...)
 *
 *
 * In both cases the condensed gas settles into a central disk where
 * star formation occurs. The cooling rate is given by the cooling-AGN
 * heating.
 *
 * BH quiet accretion rate:
 *  \f$\dot{m}_{\rm BH,R}=k_{\rm AGN}\left(\frac{m_{\rm
    BH}}{10^8\,M_{\odot}}\right)\left(\frac{f_{\rm
    hot}}{0.1}\right)\left(\frac{V_{\rm
    vir}}{200\,\mathrm{km\,s}^{-1}}\right)^3\f$

   Luminosity from quiet accretion:
   \f$L_{\rm BH}=\eta\,\dot{m}_{\rm BH,R}\,c^2,\f$

   Corresponding reduction in cooling:
    \f$\dot{m}'_{\rm{cool}}=
    \dot{m}_{\rm{cool}}-\frac{L_{BH}}{\frac{1}{2}V^2_{\rm{vir}}}\f$
 *
 *
 * The metal dependent gas cooling rates are interpolated from tables
 * read from files.
*/


#define COOLING_FUNC_N_METALS_ 8
#define COOLING_FUNC_N_TEMPRS_ 91
#define COOLING_FUNC_MIN_LOG10_T_ 4.0
#define COOLING_FUNC_MAX_LOG10_T_ 8.5
// #define COOLING_FUNC_INV_DLOG10_T_ ((COOLING_FUNC_N_TEMPRS_ - 1) / (COOLING_FUNC_MAX_LOG10_T_ - COOLING_FUNC_MIN_LOG10_T_))
#define COOLING_FUNC_INV_DLOG10_T_ 20.


static const char *cooling_function_file_name_[COOLING_FUNC_N_METALS_] = {
  "stripped_mzero.cie",
  "stripped_m-30.cie",
  "stripped_m-20.cie",
  "stripped_m-15.cie",
  "stripped_m-10.cie",
  "stripped_m-05.cie",
  "stripped_m-00.cie",
  "stripped_m+05.cie"
};


/* metallicies with repect to solar converted to
   absolut metallicity by adding  log10(Z_sun) = -1.6989700043360188 for Zsun=0.02 */
static const double cooling_function_metallicity_[COOLING_FUNC_N_METALS_] = {
  -5.0 -1.6989700043360188,                         /* actually primordial -> -infinity */
  -3.0 -1.6989700043360188,
  -2.0 -1.6989700043360188,
  -1.5 -1.6989700043360188,
  -1.0 -1.6989700043360188,
  -0.5 -1.6989700043360188,
  +0.0 -1.6989700043360188,
  +0.5 -1.6989700043360188
};


static double cooling_function_table_[COOLING_FUNC_N_METALS_][COOLING_FUNC_N_TEMPRS_];


/**@brief reads the gas cooling functions from files*/
void read_cooling_functions(void)
{
  FILE *cooling_function_file_;
  char full_cooling_function_file_name_[200];
  int i_metal_, i_temp_;
  float sd_logT_, sd_ne_, sd_nh_, sd_nt_, sd_logLnet_, sd_logLnorm_, sd_logU_, sd_logTau_, sd_logP12_, sd_logRho24_, sd_ci_, sd_mubar_;

  for(i_metal_ = 0; i_metal_ < COOLING_FUNC_N_METALS_; i_metal_++)
  {
    sprintf(full_cooling_function_file_name_, "%s/%s",CoolFunctionsDir, cooling_function_file_name_[i_metal_]);

    if(!(cooling_function_file_ = fopen(full_cooling_function_file_name_, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", full_cooling_function_file_name_);
      terminate(sbuf);
    }
    for(i_temp_ = 0; i_temp_ < COOLING_FUNC_N_TEMPRS_; i_temp_++)
    {
      fscanf(cooling_function_file_, " %f %f %f %f %f %f %f %f %f %f %f %f ",
             &sd_logT_, &sd_ne_, &sd_nh_, &sd_nt_,
             &sd_logLnet_, &sd_logLnorm_, &sd_logU_,
             &sd_logTau_, &sd_logP12_, &sd_logRho24_, &sd_ci_, &sd_mubar_);

      cooling_function_table_[i_metal_][i_temp_] = sd_logLnorm_;
    }
    fclose(cooling_function_file_);
  }

  if(ThisTask == 0)
    printf("cooling functions read.\n\n");
}


/** @brief compute metal dependent cooling rate temperature part
 * 
 *  @param [in] i_metal_ index for metallicity
 *  @param [in] log10_Z_ log10(metal fraction) 
 */
static inline double 
get_metaldependent_cooling_rate_temperature_part(const int i_metal_, const double log10_T_)
{
  if(log10_T_ <= COOLING_FUNC_MIN_LOG10_T_)
  { return cooling_function_table_[i_metal_][0]; }
  else if(log10_T_ >= COOLING_FUNC_MAX_LOG10_T_)
  { return cooling_function_table_[i_metal_][COOLING_FUNC_N_TEMPRS_ - 1]; }
  else
  {
    const int    i_temp_ = COOLING_FUNC_INV_DLOG10_T_ * (log10_T_ - COOLING_FUNC_MIN_LOG10_T_);
    const double f_      = COOLING_FUNC_INV_DLOG10_T_ * (log10_T_ - COOLING_FUNC_MIN_LOG10_T_) - i_temp_;
    return (1. - f_) * cooling_function_table_[i_metal_][i_temp_] + f_ * cooling_function_table_[i_metal_][i_temp_ + 1];
  }
}


/** @brief compute metal dependent cooling rate
 * 
 *  @param [in] log10_T_ log10(temperature/Kelvin)
 *  @param [in] log10_Z_ log10(metal fraction) 
 */
double get_metaldependent_cooling_rate(const double log10_T_, const double log10_Z_)
{
  if(log10_Z_ <= cooling_function_metallicity_[0])
  { return pow(10., get_metaldependent_cooling_rate_temperature_part(0, log10_T_)); }
  else if(log10_Z_ >= cooling_function_metallicity_[COOLING_FUNC_N_METALS_ - 1])
  { return pow(10., get_metaldependent_cooling_rate_temperature_part(COOLING_FUNC_N_METALS_ - 1, log10_T_)); }
  else
  {
    int i_metal_ = 0;
    while(log10_Z_ > cooling_function_metallicity_[i_metal_ + 1]) ++i_metal_;

    const double rate1 = get_metaldependent_cooling_rate_temperature_part(i_metal_    , log10_T_);
    const double rate2 = get_metaldependent_cooling_rate_temperature_part(i_metal_ + 1, log10_T_);
 
    return pow(10., rate1 + (rate2 - rate1) * (log10_Z_ - cooling_function_metallicity_[i_metal_]) / (cooling_function_metallicity_[i_metal_ + 1] - cooling_function_metallicity_[i_metal_]));
  }
}


/** @brief print cooling rates */
void test_metaldependent_cooling_rate(void)
{
  read_cooling_functions();

  double z_;
  for(z_ = -8.0; z_ < 2.0; z_ += 0.5)
    printf("z=%g\t rate= %g\n", z_, log10(get_metaldependent_cooling_rate(1.0, z_)));
}


/** @brief main cooling recipe, where the cooling rates are calculated */
void compute_cooling(const int galaxy_number_, const double dt_)
{
  mass_checks("cooling_recipe #1",galaxy_number_);

  Gal[galaxy_number_].CoolingGas = 0.;
  
  if(dt_ <= 0)
    return;

  if(Gal[galaxy_number_].HotGas <= 1.0e-6) /* no hot gas */
    return;
 
  /* get temp -> Temperature of the Gas in Kelvin, obtained from
   * hydrostatic equilibrium KT= 0.5*mu_p*(Vc)^2 assuming Vvir~Vc */  
  const double Vvir        = Gal[galaxy_number_].Vvir;
  const double temp        = 35.9 * Vvir * Vvir;

  if(temp < 1.e4) /* below photoionizing background */
   return;
  
  const double Rvir        = Gal[galaxy_number_].Rvir;
  const double tot_hotMass = Gal[galaxy_number_].HotGas;
  const double tot_metals  = metals_total(Gal[galaxy_number_].MetalsHotGas);

// const double tcool = Rvir / Vvir; // tcool = t_dynamical = Rvir/Vvir
// const double HotRadius = (Gal[galaxy_number_].Type == 0) ? Rvir : Gal[galaxy_number_].HotRadius;
// const double logZ      = (tot_metals > 0) ? log10(tot_metals / tot_hotMass) : -10.0;
// 
//   //eq. 3 and 4 Guo2010
//   const double lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
//   const double x = PROTONMASS * BOLTZMANN * temp / (lambda * UnitDensity_in_cgs * UnitTime_in_s);  // in internal units
//   const double rho_rcool = x / (0.28086 * tcool);
//   /* an isothermal density profile for the hot gas is assumed here */
//   const double rho0 = tot_hotMass / (4 * M_PI * HotRadius);
//   const double rcool = sqrt(rho0 / rho_rcool);
//   
//   if (Gal[galaxy_number_].CoolingRadius < rcool)
//     Gal[galaxy_number_].CoolingRadius = rcool;
//       
//   double coolingGas;
//   //if Hotradius is used, when galaxies become type 1's there will be a discontinuity in the cooling
//   if(rcool > Rvir) // INFALL DOMINATED REGIME
//     //coolingGas = tot_hotMass; - Delucia 2007
//     /*comes in to keep the continuity (Delucia2004) */
//     // coolingGas = tot_hotMass / (HotRadius / Vvir) * dt_;
//     coolingGas = tot_hotMass / HotRadius * Vvir * dt_;
//   else // HOT PHASE REGIME
//     /*coolingGas = (tot_hotMass / Rvir) * (rcool / tcool) * dt_ */
//     // coolingGas = (tot_hotMass / HotRadius) * (rcool / tcool) * dt_ ;
//     coolingGas = tot_hotMass * rcool / (HotRadius * Rvir) * Vvir * dt_ ;
//  
//   if(coolingGas > tot_hotMass)
//     coolingGas = tot_hotMass;
//   else if(coolingGas < 0.0)
//     coolingGas = 0.0;      

  // const double tcool = Rvir / Vvir; // tcool = t_dynamical = Rvir/Vvir
  const double HotRadius = (Gal[galaxy_number_].Type == 0) ? Rvir : Gal[galaxy_number_].HotRadius;
  const double logZ      = (tot_metals > 0) ? log10(tot_metals / tot_hotMass) : -10.0;

  //eq. 3 and 4 Guo2010
  const double lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
  // const double x = PROTONMASS * BOLTZMANN * temp / (lambda * UnitDensity_in_cgs * UnitTime_in_s);  // in internal units
  // const double x = PROTONMASS * BOLTZMANN * 35.9 * Vvir * Vvir / (lambda * UnitDensity_in_cgs * UnitTime_in_s);  // in internal units
  // const double rho_rcool = x / (0.28086 * tcool);
  // const double rho_rcool = (PROTONMASS * BOLTZMANN * 35.9 / (0.28086 * UnitDensity_in_cgs * UnitTime_in_s)) * Vvir * Vvir * Vvir / (lambda * Rvir);

  /* an isothermal density profile for the hot gas is assumed here */
  // const double rho0 = tot_hotMass / (4 * M_PI * HotRadius);
  // const double rcool = sqrt(rho0 / rho_rcool);
  const double rcool = sqrt(((0.28086 * UnitDensity_in_cgs * UnitTime_in_s) / (4 * M_PI * PROTONMASS * BOLTZMANN * 35.9)) * tot_hotMass * lambda * Rvir / (HotRadius * Vvir * Vvir * Vvir));
  
  if (Gal[galaxy_number_].CoolingRadius < rcool)
    Gal[galaxy_number_].CoolingRadius = rcool;
      
  double coolingGas;
  //if Hotradius is used, when galaxies become type 1's there will be a discontinuity in the cooling
  if(rcool > Rvir) // INFALL DOMINATED REGIME
    //coolingGas = tot_hotMass; - Delucia 2007
    /*comes in to keep the continuity (Delucia2004) */
    // coolingGas = tot_hotMass / (HotRadius / Vvir) * dt_;
    coolingGas = tot_hotMass / HotRadius * Vvir * dt_;
  else // HOT PHASE REGIME
    /*coolingGas = (tot_hotMass / Rvir) * (rcool / tcool) * dt_ */
    // coolingGas = (tot_hotMass / HotRadius) * (rcool / tcool) * dt_ ;
    coolingGas = tot_hotMass * rcool / (HotRadius * Rvir) * Vvir * dt_ ;

  if(coolingGas > tot_hotMass)
    coolingGas = tot_hotMass;

  Gal[galaxy_number_].CoolingGas = coolingGas;

  mass_checks("cooling_recipe #1.5", galaxy_number_);
}   


/** @brief calculates the energy released by black holes due to passive accretion,
  *
  * do_AGN_heating calculates the amount of energy released by
  * black holes due to passive accretion, Which is then used to reduce
  * the cooling.
  *
  * There is one parameter, AgnEfficiency, which is the efficiency of AGN
  * passive accretion and consequently of cooling flow reheating, Eqs. 10,
  * 11 & 12 in Croton 2006. There is a AGNRadioModeModel =2 (empirical)
  * with options 3 and 4 representing Bondi-Hoyle and cold cloud accretion.
  * The three should be identical and the use of empirical avoids people
  * shouting about duty-cycles being inconsistent with Bondi-Hoyle & etc.
  *
  * AGN accretion and heating is assumed to go from the gas from the type 0
  * galaxy at the centre of the FOF group to the most massive black hole
  * inside Rvir. If the most massive black hole is in a type 1 this produces
  * all the heating affecting the type 0 galaxy and will add to the accretion
  * occuring on the type 1,
  * this is done to account for the fact that the centres of FOF groups can
  * switch between galaxies. As a result a very small galaxy without a black hole
  * might be assigned the centre of a cluster leading to huge cooling. It is therefore
  * not necessary to do the same correction for satellites of subhalos.
  * 
  * @bug (fixed?) dist was used in comparison, but never set,
  *      now dist is computed as distance between galaxy and FoF center
  */
void do_AGN_heating(double dt_, const int n_galaxies_in_fof_group_)
{
  double AGNrate, AGNheating, AGNaccreted, fraction, FreeFallRadius;
  double HotGas, HotRadius, Rvir, Vvir, Mvir;
  double LeftOverEnergy, CoolingGas;
  
  /* caution: this is set to 0 here to avoid -Wmaybe-uninitialized for gcc,
   * which can't figure out that it's only relevant for AGNRadioModeModel == 0,
   * in which case right below, it is indeed set
   * (and the the program terminates, if setting to a valid value fails) */
  int FoFCentralGal = 0; 
  if(AGNRadioModeModel == 0)
  {
    for (FoFCentralGal = 0; FoFCentralGal < n_galaxies_in_fof_group_ && Gal[FoFCentralGal].Type != 0; FoFCentralGal++)  {}
    if(FoFCentralGal >= n_galaxies_in_fof_group_) { terminate("FoFCentralGal not found."); }
  }

  int galaxy_number_;
  for (galaxy_number_ = 0; galaxy_number_ < n_galaxies_in_fof_group_; galaxy_number_++)
  {
    Gal[galaxy_number_].CoolingRate_beforeAGN += Gal[galaxy_number_].CoolingGas / (dt_*STEPS);

    AGNrate=0.;
    AGNheating = 0.;
    LeftOverEnergy = 0.;

    HotGas     = Gal[galaxy_number_].HotGas;
    HotRadius  = Gal[galaxy_number_].HotRadius;
    CoolingGas = Gal[galaxy_number_].CoolingGas;
    Mvir       = Gal[galaxy_number_].Mvir;
    Rvir       = Gal[galaxy_number_].Rvir;
    Vvir       = Gal[galaxy_number_].Vvir;

    if(HotGas > 0.0)
    {
      if(AGNRadioModeModel == 0)
      {
        // AGNrate = AgnEfficiency * UnitTime_in_s * Gal[galaxy_number_].BlackHoleMass / Hubble_h * (HotGas/Hubble_h) *(10. * SOLAR_MASS / (UNITMASS_IN_G * SEC_PER_YEAR));
        AGNrate = AgnEfficiency * Gal[galaxy_number_].BlackHoleMass * HotGas *  inv_Hubble_h * inv_Hubble_h * (10. * SolarMass * UnitTime_in_s / SEC_PER_YEAR);
      }
      else if(AGNRadioModeModel == 2)
      {
        //empirical (standard) accretion recipe - Eq. 10 in Croton 2006
//         AGNrate = AgnEfficiency / (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
//                       * (Gal[galaxy_number_].BlackHoleMass / 0.01) * pow3(Vvir / 200.0)
//                       * ((HotGas / HotRadius * Rvir / Mvir) / 0.1);
                      
        // AGNrate = AgnEfficiency * UnitTime_in_s * Gal[galaxy_number_].BlackHoleMass * pow3(Vvir) * HotGas * Rvir / (HotRadius * Mvir) * (0.000125 * SOLAR_MASS / (UNITMASS_IN_G * SEC_PER_YEAR));
        AGNrate = AgnEfficiency * Gal[galaxy_number_].BlackHoleMass * pow3(Vvir) * HotGas * Rvir / (HotRadius * Mvir) * (0.000125 * SolarMass * UnitTime_in_s / SEC_PER_YEAR);
      }
      else if(AGNRadioModeModel == 3 || AGNRadioModeModel == 4)
      {
        const double tot_metals = metals_total(Gal[galaxy_number_].MetalsHotGas);

        /* temp -> Temperature of the Gas in Kelvin, obtained from
         * hydrostatic equilibrium KT=0.5*mu_p*(Vc)^2 assuming Vvir~Vc */
        const double temp   = 35.9 * Vvir * Vvir;
        const double logZ   = (tot_metals > 0) ? log10(tot_metals / HotGas): -10.0;
        const double lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
        // const double x      = PROTONMASS * BOLTZMANN * temp / (lambda * UnitDensity_in_cgs * UnitTime_in_s);  // in internal units
        const double x      = (PROTONMASS * BOLTZMANN / (UnitDensity_in_cgs * UnitTime_in_s)) * temp / lambda;  // in internal units

        /* Bondi-Hoyle accretion recipe -- efficiency = 0.15
         * Eq. 29 in Croton 2006 */
        if(AGNRadioModeModel == 3)
//        { AGNrate = (2.5 * M_PI * Gravity) * (0.75 * 0.6 * x) * Gal[galaxy_number_].BlackHoleMass * 0.15; }
        { AGNrate = (2.5 * M_PI * 0.75 * 0.6 * 0.15 * Gravity) * x * Gal[galaxy_number_].BlackHoleMass; }
        else if(AGNRadioModeModel == 4)
        {
          /* Cold cloud accretion recipe -- trigger: Rff = 50 Rdisk,
           * and accretion rate = 0.01% cooling rate
           * Eq. 25 in Croton 2006 */
          // FreeFallRadius = HotGas / (6.0 * 0.6 * x * Rvir * Vvir) / HotRadius * Rvir;
          FreeFallRadius = HotGas * Rvir / (6.0 * 0.6 * x * Rvir * Vvir * HotRadius);
          if(Gal[galaxy_number_].BlackHoleMass > 0.0 && FreeFallRadius < Gal[galaxy_number_].GasDiskRadius * 50.0)
          { AGNrate = 0.0001 * CoolingGas / dt_; }
          else
          { AGNrate = 0.0; }
        }
      }

      /* Eddington rate */
      /* Note that this assumes an efficiency of 50%
       * - it ignores the e/(1-e) factor in L = e/(1-e) Mdot c^2 */
      // const double EDDrate = 1.3e48 * Gal[galaxy_number_].BlackHoleMass / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;
      const double EDDrate = (1.3e48 / 9e10 * UnitTime_in_s / UnitEnergy_in_cgs) * Gal[galaxy_number_].BlackHoleMass;

      /* accretion onto BH is always limited by the Eddington rate */
      if(AGNrate > EDDrate)
      { AGNrate = EDDrate; }

      /*  accreted mass onto black hole the value of dt_ puts an h factor into AGNaccreted as required for code units */
      AGNaccreted = AGNrate * dt_;

      /* cannot accrete more mass than is available! */
      if(AGNaccreted > HotGas)
      { AGNaccreted = HotGas; }

      /*  coefficient to heat the cooling gas back to the virial temperature of the halo */
      /*  1.34e5 = sqrt(2*eta*c^2), eta=0.1 (standard efficiency) and c in km/s
       *  Eqs. 11 & 12 in Croton 2006 */
      // const double AGNcoeff = (1.34e5 / Vvir) * (1.34e5 / Vvir);
      const double AGNcoeff = (1.34e5 * 1.34e5) / (Vvir * Vvir);

      /*  cooling mass that can be suppressed from AGN heating */
      AGNheating = AGNcoeff * AGNaccreted;

      if(AGNRadioModeModel == 0 && Gal[galaxy_number_].Type==1)
      {
        if(separation_gal(galaxy_number_, FoFCentralGal) < Gal[FoFCentralGal].Rvir)
        {
          if(AGNheating > (Gal[galaxy_number_].CoolingGas + Gal[FoFCentralGal].CoolingGas))
          {
            AGNheating  = (Gal[galaxy_number_].CoolingGas + Gal[FoFCentralGal].CoolingGas);
            AGNaccreted = (Gal[galaxy_number_].CoolingGas + Gal[FoFCentralGal].CoolingGas) / AGNcoeff;
          }
          if(AGNheating > Gal[galaxy_number_].CoolingGas)
          { LeftOverEnergy = AGNheating - Gal[galaxy_number_].CoolingGas; }
        }
      }
      else if(AGNheating > Gal[galaxy_number_].CoolingGas)
      { AGNaccreted = Gal[galaxy_number_].CoolingGas / AGNcoeff; }

      /*  accreted mass onto black hole */
      Gal[galaxy_number_].BlackHoleMass += AGNaccreted; //ROB: transfer_mass functions should be used here
      Gal[galaxy_number_].RadioAccretionRate += AGNaccreted / (dt_*STEPS);
      fraction=AGNaccreted/Gal[galaxy_number_].HotGas;
      Gal[galaxy_number_].HotGas -= AGNaccreted;
      metals_add_fraction_to(&Gal[galaxy_number_].MetalsHotGas, Gal[galaxy_number_].MetalsHotGas, -fraction);

#ifdef INDIVIDUAL_ELEMENTS
      elements_add_fraction_to(&Gal[galaxy_number_].HotGas_elements,Gal[galaxy_number_].HotGas_elements,-fraction);
#endif
#ifdef METALS_SELF
      metals_add_fraction_to(&Gal[galaxy_number_].MetalsHotGasSelf,Gal[galaxy_number_].MetalsHotGasSelf,-fraction);
#endif 
    }
    
    /* limit heating to cooling rate */
    if(Gal[galaxy_number_].CoolingGas < AGNheating)
    { Gal[galaxy_number_].CoolingGas = 0.0; }
    else
    {
      Gal[galaxy_number_].CoolingGas -= AGNheating; 
      Gal[galaxy_number_].CoolingRate += Gal[galaxy_number_].CoolingGas / (dt_*STEPS);
    }

    if(AGNRadioModeModel == 0 && LeftOverEnergy > 0.)
    {
      if(Gal[FoFCentralGal].CoolingGas < LeftOverEnergy)
      { Gal[FoFCentralGal].CoolingGas = 0.0; }
      else
      { 
        Gal[FoFCentralGal].CoolingGas  -= LeftOverEnergy;
        Gal[FoFCentralGal].CoolingRate -= LeftOverEnergy / (dt_*STEPS);
      }
    }
    mass_checks("cooling_recipe #2.",galaxy_number_);
  }
}


/** @brief updates the fractions of hot and cold gas due to cooling.
 * 
 * cool_gas_onto_galaxy updates the fractions of hot and cold gas
 * due to cooling. This is done for the mass, metals and, after Guo2010,
 * spin components.
 */
void cool_gas_onto_galaxy(int galaxy_number_, double dt_)
{
  int i;

  const double Mdisk = Gal[galaxy_number_].ColdGas;
  const double Mcool = min(Gal[galaxy_number_].CoolingGas, Gal[galaxy_number_].HotGas);

  /*  add the fraction 1/STEPS of the total cooling gas to the cold disk */
  if(Mcool > 0.0)
  {
    //determine the xray luminosity of any cooling gas in this snapshot (White & Frenk 1991 eq21)
    // Gal[galaxy_number_].XrayLum = log10(2.5 * (Mcool / dt_) * 6.31 * Gal[galaxy_number_].Vvir * Gal[galaxy_number_].Vvir) + 35.0;
    Gal[galaxy_number_].XrayLum = log10(15.775 * Mcool * Gal[galaxy_number_].Vvir * Gal[galaxy_number_].Vvir / dt_) + 35.0;

    // We already know that 0<mcool<=Gal[galaxy_number_].HotGas
    const double fraction = ((float)Mcool) / Gal[galaxy_number_].HotGas;
    transfer_gas(galaxy_number_, ColdGasComponent, galaxy_number_, HotGasComponent, fraction);

    if (DiskRadiusModel == 0)
    {
      if (Gal[galaxy_number_].ColdGas != 0.0)
      {
        for (i=0;i<3;i++) 
        { 
          Gal[galaxy_number_].GasSpin[i] = (Gal[galaxy_number_].GasSpin[i] * Mdisk + Gal[galaxy_number_].HaloSpin[i] * Mcool) / Gal[galaxy_number_].ColdGas;
        }
      }
      set_gas_disk_radius(galaxy_number_);
    }
  }
  else
  { Gal[galaxy_number_].XrayLum = 0.0; }
}

