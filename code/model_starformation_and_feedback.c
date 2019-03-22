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
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_starformation_and_feedback.c
 *  @brief recipe_starformation_and_feedback.c computes the amount of stars_
 *         formed from the cold gas, the amount of gas reheated from cold to hot
 *         and the amount of gas ejected from hot to external.
 *
 * The routine is divided in two parts, star formation and SN feedback, with a
 * number of different implementations controlled by input parameters.
 *
 *
 *  0 -\f$M_{\rm{crit}}=3.8\times 10^9
 *     \left(\frac{V_{\rm{max}}}{200\,\rm{km s}^{-1}}\right)
 *     \left(\frac{r_{\rm{disk}}}{10\,\rm{kpc}}\right)M_{\odot}\f$
 *     (Eq. 16 Guo2010) (StarFormationModel = 0), \n
 *        - same as 1 but using \f$V_{\rm{max}}\f$ or \f$V_{\rm{max,infall}}\f$
 *          instead of \f$V_{\rm{vir}}\f$ and allowing SF in satellites. *

 *
 * There are 2 options for the <B>SN Feedback Recipe</B>:
 *
 * 0 - \f$\epsilon_{\rm{disk}}=\epsilon
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_1}\biggr]\f$,
 *     \f$\epsilon_{\rm{halo}}=\eta
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_2}\biggr]\f$
 *     (Eqs. 19 & 21 Guo2010)(FeedbackEjectionModel = 2)
 *     same as FeedbackEjectionModel = 1 * Vmax dependence.
 *
 * Also, Guo2010 alowed for type 1 satellite to have gas cycles and receive
 * gas from their own satellites when these are outside Rvir of the type 0.
 * */


/** @brief Main recipe, calculates the fraction_ of cold gas turned into stars due
  *        to star formation; the fraction_ of mass instantaneously recycled and
  *        returned to the cold gas; the fraction_ of gas reheated from cold to hot,
  *        ejected from hot to external and returned from ejected to hot due to
  *        SN feedback.   
  *
  * @note changing v_max_ from float to double seemed to impact t_dyn_ 
  */
void starformation(const int galaxy_number_, const int central_galaxy_number_, const double time_, const double dt_)
{
        /** Variables: reff-Rdisk, t_dyn_=Rdisk/Vmax, strdot=Mstar_dot, stars_=strdot*dt_ */
  /* Note that Units of dynamical time_ are Mpc/Km/s - no conversion on dt_ needed */
  double stars_;
 
  // if(Gal[galaxy_number_].Type == 0)
  // {
  //   t_dyn_ = Gal[galaxy_number_].GasDiskRadius / Gal[galaxy_number_].Vmax;
  //   cold_crit_ = SfrColdCrit * Gal[galaxy_number_].Vmax/200. * Gal[galaxy_number_].GasDiskRadius*100.;
  // }
  // else
  // {
  //    t_dyn_ = Gal[galaxy_number_].GasDiskRadius / Gal[galaxy_number_].InfallVmax;
  //   cold_crit_ = SfrColdCrit * Gal[galaxy_number_].InfallVmax/200. * Gal[galaxy_number_].GasDiskRadius*100.;
  // }
  
  const double v_max_     = (Gal[galaxy_number_].Type == 0) ? Gal[galaxy_number_].Vmax : Gal[galaxy_number_].InfallVmax;
  const double inv_t_dyn_ = v_max_ / Gal[galaxy_number_].GasDiskRadius;
  const double cold_crit_ = 0.5 * SfrColdCrit * v_max_ * Gal[galaxy_number_].GasDiskRadius;
  
  //standard star formation law (Croton2006, Delucia2007, Guo2010)
  /* if(StarFormationModel == 0) */
  { stars_ = (Gal[galaxy_number_].ColdGas > cold_crit_) ? SfrEfficiency * (Gal[galaxy_number_].ColdGas - cold_crit_) * dt_ * inv_t_dyn_ : 0.0; }
  /*else if(StarFormationModel == 1)
  {
          ALTERNATIVE STAR FORMATION LAW 
  }*/

  if(stars_ > 0.)
  {
//otherwise cold gas and stars_ share material in update_stars_due_to_reheat
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
    if(stars_ > Gal[galaxy_number_].ColdGas)
    { stars_ = Gal[galaxy_number_].ColdGas; }
#endif

    mass_checks("recipe_starform #1",galaxy_number_);
    mass_checks("recipe_starform #1.1",central_galaxy_number_);

  /* update for star formation
   * updates Mcold, StellarMass, MetalsMcold and MetalsStellarMass
   * in Guo2010 case updates the stellar spin -> hardwired, not an option */

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
    /* Store the value of the metallicity of the cold phase when SF occurs */
    const double metallicity_ = metals_total(Gal[galaxy_number_].MetalsColdGas) / Gal[galaxy_number_].ColdGas;
#endif //NDEF POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

    update_stars_due_to_reheat(galaxy_number_, &stars_);

    mass_checks("recipe_starform #2",galaxy_number_);
    mass_checks("recipe_starform #2.1",central_galaxy_number_);

    /*  update the star formation rate */
    /*Sfr=stars_/(dt_*steps)=strdot*dt_/(dt_*steps)=strdot/steps -> average over the STEPS*/
    Gal[galaxy_number_].Sfr += stars_ / (dt_ * STEPS);

    // update_from_star_formation can only be called
    // after SD_feeedback recipe since stars_ need to be reset once the reheated mass is known
    // (star formation and feedback share the same fraction_ of cold gas)
    update_from_star_formation(galaxy_number_, stars_, false); // false indicates not a burst

    update_massweightage(galaxy_number_, stars_, time_);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  /* ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN feedback is only called
   * when stars_ die, inside DETAILED_METALS_AND_MASS_RETURN */
    SN_feedback(galaxy_number_, central_galaxy_number_, stars_, ColdGasComponent);
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  /*  Update the luminosities due to the stars_ formed */
    add_to_luminosities(galaxy_number_, stars_, time_, dt_, metallicity_);
#endif //NDEF POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
  }

  if(Gal[galaxy_number_].DiskMass > 0.0)
  { check_disk_instability(galaxy_number_); }

  if (DiskRadiusModel== 0)
  { set_stellar_disk_radius(galaxy_number_); }
}


 /** @brief update mass of stars_ formed due to reheating 
  * 
  * @bug Which is the true central galaxy? Gal[galaxy_number_].CentralGal or Gal[central_galaxy_number_]?
  *       Tests show, these are not always the same.
  */
void update_stars_due_to_reheat(const int galaxy_number_, double *stars_)
{
//   double MergeCentralVvir=0.;
//   double CentralVvir=0.;
//   doble ejected_mass_=0.;
  /* SN FEEDBACK RECIPES */

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */

  //REHEAT
//   CentralVvir = Gal[central_galaxy_number_].Vvir; // main halo Vvir
//   MergeCentralVvir = Gal[Gal[galaxy_number_].CentralGal].Vvir; //central subhalo Vvir

  // Feedback depends on the circular velocity of the host halo
  // Guo2010 - eq 18 & 19
  if(FeedbackReheatingModel == 0)
  {
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
    double reheated_mass_;
    if(Gal[Gal[galaxy_number_].CentralGal].Type == 0)
    { reheated_mass_ = FeedbackReheatingEpsilon * (*stars_) * (0.5 + pow(Gal[Gal[galaxy_number_].CentralGal].Vmax      / ReheatPreVelocity, -ReheatSlope)); }
    else
    { reheated_mass_ = FeedbackReheatingEpsilon * (*stars_) * (0.5 + pow(Gal[Gal[galaxy_number_].CentralGal].InfallVmax / ReheatPreVelocity, -ReheatSlope)); }
  
    if(reheated_mass_ * Gal[Gal[galaxy_number_].CentralGal].Vvir * Gal[Gal[galaxy_number_].CentralGal].Vvir > (*stars_) * EtaSNcode * EnergySNcode)
    { reheated_mass_ = (*stars_) * EtaSNcode * EnergySNcode / (Gal[Gal[galaxy_number_].CentralGal].Vvir * Gal[Gal[galaxy_number_].CentralGal].Vvir); }

    if((*stars_ + reheated_mass_) > Gal[galaxy_number_].ColdGas)
    { *stars_ *= Gal[galaxy_number_].ColdGas / (*stars_ + reheated_mass_); }
#endif
  }
  /*else if(FeedbackReheatingModel == 1)
  {
    reheated_mass_ =  ALTERNATIVE Reheating LAW ;
  }*/
}


/** @brief Updates the different components due to star formation: mass
  *        and metals in stars_ and cold gas and stellar spin. */
//void update_from_star_formation(int galaxy_number_, double time_, double stars_, double metallicity)
void update_from_star_formation(const int galaxy_number_, const double stars_, const bool flag_burst_)
{
  if(Gal[galaxy_number_].ColdGas <= 0. || stars_ <= 0.)
  {
    printf("update_from_star_formation: Gal[galaxy_number_].ColdGas <= 0. || stars_ <= 0.\n");
    exit(0);
  }

  /* If DETAILED_METALS_AND_MASS_RETURN, no longer an assumed instantaneous
   * recycled fraction_. Mass is returned over time_ via SNe and AGB winds.
   * Update the Stellar Spin when forming stars_ */
#ifndef DETAILED_METALS_AND_MASS_RETURN
  const double stars_to_add_ = (1 - RecycleFraction) * stars_;
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
  const double stars_to_add_ = stars_;
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */

  if (Gal[galaxy_number_].DiskMass+stars_to_add_ > 1.e-8)
  {
    const double inv_new_disk_mass_ = 1. / (Gal[galaxy_number_].DiskMass + stars_to_add_);
    int i;
    for (i = 0; i < 3; i++)
    { Gal[galaxy_number_].StellarSpin[i]= (Gal[galaxy_number_].StellarSpin[i] * Gal[galaxy_number_].DiskMass + stars_to_add_ * Gal[galaxy_number_].GasSpin[i]) * inv_new_disk_mass_; }
  }
    /*  Update Gas and Metals from star formation */
  mass_checks("update_from_star_formation #0",galaxy_number_);

  const double fraction_ = stars_to_add_ / Gal[galaxy_number_].ColdGas;

#ifdef STAR_FORMATION_HISTORY
  Gal[galaxy_number_].sfh_DiskMass[Gal[galaxy_number_].sfh_ibin] += stars_to_add_; //ROB: Now, all SF gas is put in SFH array ("recycled' mass will return to gas phase over time_)
  metals_add_fraction_to(&Gal[galaxy_number_].sfh_MetalsDiskMass[Gal[galaxy_number_].sfh_ibin], Gal[galaxy_number_].MetalsColdGas, fraction_);
#ifdef INDIVIDUAL_ELEMENTS
  elements_add_fraction_to(&Gal[galaxy_number_].sfh_ElementsDiskMass[Gal[galaxy_number_].sfh_ibin], Gal[galaxy_number_].ColdGas_elements, fraction_);
#endif /* defined INDIVIDUAL_ELEMENTS */
#ifdef TRACK_BURST
  if (flag_burst_)
  { Gal[galaxy_number_].sfh_BurstMass[Gal[galaxy_number_].sfh_ibin] += stars_to_add_; }
#else  /* not defined TRACK_BURST */
  (void) flag_burst_; /* avoid unused-parameter warning */
#endif /* not defined TRACK_BURST */
#else  /* not defined  STAR_FORMATION_HISTORY */
  (void) flag_burst_; /* avoid unused-parameter warning */
#endif /* not defined STAR_FORMATION_HISTORY */

  metals_add_fraction_to(&Gal[galaxy_number_].MetalsDiskMass, Gal[galaxy_number_].MetalsColdGas,  fraction_);
  metals_add_fraction_to(&Gal[galaxy_number_].MetalsColdGas , Gal[galaxy_number_].MetalsColdGas, -fraction_);

  //GLOBAL PROPERTIES
  Gal[galaxy_number_].DiskMass += stars_to_add_;
  Gal[galaxy_number_].ColdGas  -= stars_to_add_;
#ifdef INDIVIDUAL_ELEMENTS
  elements_add_fraction_to(&Gal[galaxy_number_].DiskMass_elements, Gal[galaxy_number_].ColdGas_elements,  fraction_);
  elements_add_fraction_to(&Gal[galaxy_number_].ColdGas_elements , Gal[galaxy_number_].ColdGas_elements, -fraction_);
#endif
#ifdef TRACK_BURST
  if (flag_burst_) Gal[galaxy_number_].BurstMass += stars_to_add_;
#endif

  mass_checks("update_from_star_formation #1",galaxy_number_);

  /* Formation of new metals - instantaneous recycling approximation - only SNII
   * Also recompute the metallicity of the cold phase.*/
#ifndef DETAILED_METALS_AND_MASS_RETURN
  /* stars_ used because the Yield is defined as a fraction_ of
   * all stars_ formed, not just long lived */
  Gal[galaxy_number_].MetalsColdGas += Yield * stars_;
#endif

  if (DiskRadiusModel == 0)
    set_stellar_disk_radius(galaxy_number_);
}


/* there are two modes for supernova feedback corresponding to when the mass returning
 * by dying stars_ is returned to the cold gas - reheat and ejection; and when the mass
 * is returned to the hot gas - onle ejection.*/
void SN_feedback(const int galaxy_number_, const int central_galaxy_number_, const double stars_, const GasComponentType feedback_location_)
{
//   double CentralVvir, MergeCentralVvir=0.;
  double EjectVmax, EjectVvir, SN_energy_, reheat_energy_;
  double reheated_mass_ = 0., ejected_mass_ = 0.;
  
  
  /* SN FEEDBACK MODEL */

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */
   
  if (feedback_location_ == HotGasComponent)
  { reheated_mass_ = 0.; }
  else if(feedback_location_ == ColdGasComponent)
  {
//   CentralVvir = Gal[central_galaxy_number_].Vvir; // main halo Vvir
//   MergeCentralVvir = Gal[Gal[galaxy_number_].CentralGal].Vvir; //central subhalo Vvir

    mass_checks("recipe_starform #0",galaxy_number_);
    mass_checks("recipe_starform #0.1",central_galaxy_number_);

    // Feedback depends on the circular velocity of the host halo
    // Guo2010 - eq 18 & 19
    /* if(FeedbackReheatingModel == 0) */
    {
      if (Gal[Gal[galaxy_number_].CentralGal].Type == 0)
      { reheated_mass_ = FeedbackReheatingEpsilon * stars_ * (0.5 + pow(Gal[Gal[galaxy_number_].CentralGal].Vmax       / ReheatPreVelocity, -ReheatSlope)); }
      else
      { reheated_mass_ = FeedbackReheatingEpsilon * stars_ * (0.5 + pow(Gal[Gal[galaxy_number_].CentralGal].InfallVmax / ReheatPreVelocity, -ReheatSlope)); }
    
      if (reheated_mass_ * Gal[Gal[galaxy_number_].CentralGal].Vvir * Gal[Gal[galaxy_number_].CentralGal].Vvir > stars_ * (EtaSNcode * EnergySNcode))
      { reheated_mass_ = stars_ * (EtaSNcode * EnergySNcode) / (Gal[Gal[galaxy_number_].CentralGal].Vvir * Gal[Gal[galaxy_number_].CentralGal].Vvir); }
    }
   /* else if(FeedbackReheatingModel == 1)
    {
      reheated_mass_ =  ALTERNATIVE Reheating LAW ;
    } */

    if(reheated_mass_ > Gal[galaxy_number_].ColdGas)
    { reheated_mass_ = Gal[galaxy_number_].ColdGas; }
  }// end if feedback_location_

  /* Determine ejection (for FeedbackEjectionModel 2 we have the dependence on Vmax)
   * Guo2010 - eq 22
   * Note that satellites can now retain gas and have their own gas cycle*/
  if (Gal[Gal[galaxy_number_].CentralGal].Type == 0)
  {
    EjectVmax=Gal[central_galaxy_number_].Vmax;
    EjectVvir=Gal[central_galaxy_number_].Vvir;// main halo Vvir
  }
  else
  {
    EjectVmax=Gal[Gal[galaxy_number_].CentralGal].InfallVmax;
    EjectVvir=Gal[Gal[galaxy_number_].CentralGal].Vvir; //central subhalo Vvir
  }

  if(FeedbackEjectionModel == 0)
  {
    // ejected_mass_ = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) * stars_ *  min(1./FeedbackEjectionEfficiency, 0.5 + pow(EjectVmax / EjectPreVelocity, -EjectSlope)) - reheated_mass_ * EjectVvir * EjectVvir) / (EjectVvir * EjectVvir);

    const double ef_1_ = 1. / FeedbackEjectionEfficiency;
    const double ef_2_ = 0.5 + pow(EjectVmax / EjectPreVelocity, -EjectSlope);
    ejected_mass_ = FeedbackEjectionEfficiency * EtaSNcode * EnergySNcode * stars_ * min(ef_1_, ef_2_) / (EjectVvir * EjectVvir) - reheated_mass_;
  }
  else if(FeedbackEjectionModel == 1)//the ejected material is assumed to have V_SN
  {
    SN_energy_     = 0.5 * stars_ * (EtaSNcode * EnergySNcode);
    reheat_energy_ = 0.5 * reheated_mass_ * EjectVvir * EjectVvir;

    ejected_mass_ = (SN_energy_ - reheat_energy_)/(0.5 * FeedbackEjectionEfficiency * EtaSNcode * EnergySNcode);

    //if VSN^2<Vvir^2 nothing is ejected
    if(FeedbackEjectionEfficiency * EtaSNcode * EnergySNcode < EjectVvir * EjectVvir)
    { ejected_mass_ = 0.0; }
  }

  // Finished calculating mass exchanges, so just check that none are negative
  if (reheated_mass_ < 0.0) reheated_mass_ = 0.0;
  if (ejected_mass_ < 0.0) ejected_mass_ = 0.0;

  /* Update For Feedback */
  /* update cold, hot, ejected gas fractions and respective metallicities
   * there are a number of changes introduced by Guo2010 concerning where
   * the gas ends up */

  //ejected_mass_ = 0.01*Gal[central_galaxy_number_].HotGas;
  if (reheated_mass_ + ejected_mass_ > 0.)
    update_from_feedback(galaxy_number_, central_galaxy_number_, reheated_mass_, ejected_mass_);
}


/** @brief Updates cold, hot and external gas components due to SN
 *         reheating and ejection. */
void update_from_feedback(const int galaxy_number_, const int central_galaxy_number_, const double reheated_mass_, double ejected_mass_)
{
  double separation_=0.;
  double mass_remaining_;
  double fraction_;
//  int merger_centre;

  //mass_checks("update_from_feedback #1",galaxy_number_);

  if(Gal[galaxy_number_].ColdGas > 0.)
  {
    //REHEAT
    // if galaxy is a type 1 or a type 2 orbiting a type 1 with hot gas being striped,
    //some of the reheated and ejected masses goes to the type 0 and some stays in the type 1

    if(Gal[galaxy_number_].Type ==0)
    {
      transfer_gas(galaxy_number_,HotGasComponent,galaxy_number_,ColdGasComponent,((float)reheated_mass_)/Gal[galaxy_number_].ColdGas);
    }
    else if(Gal[galaxy_number_].Type <3)
    {
//           if(Gal[galaxy_number_].Type ==1)
//             merger_centre=central_galaxy_number_;
//           else if(Gal[galaxy_number_].Type ==2)
//             merger_centre=Gal[galaxy_number_].CentralGal;

      separation_=separation_gal(central_galaxy_number_,Gal[galaxy_number_].CentralGal)/(1+ZZ[Halo[Gal[central_galaxy_number_].HaloNr].SnapNum]);

      //compute share of reheated mass
      if (separation_<Gal[central_galaxy_number_].Rvir && Gal[Gal[galaxy_number_].CentralGal].Type == 1)
      {
        //mass that remains on type1 (the rest goes to type 0) for reheat - mass_remaining_, for eject - ejected mass
        mass_remaining_=reheated_mass_*Gal[galaxy_number_].HotRadius/Gal[galaxy_number_].Rvir;
        ejected_mass_ = ejected_mass_*Gal[galaxy_number_].HotRadius/Gal[galaxy_number_].Rvir;

        if (mass_remaining_ > reheated_mass_)
        { mass_remaining_ = reheated_mass_; }
      }
      else
      { mass_remaining_=reheated_mass_; }

      //needed due to precision issues, since we first remove mass_remaining_ and then (reheated_mass_-mass_remaining_)
      //from the satellite into the type 0 and type 1 the fraction_ might not add up on the second call
      //since Gal[galaxy_number_].ColdGas is a float and reheated_mass_ & mass_remaining_ are doubles
      if((mass_remaining_ + reheated_mass_)>Gal[galaxy_number_].ColdGas)
        mass_remaining_=Gal[galaxy_number_].ColdGas-reheated_mass_;

      //transfer mass_remaining_
      transfer_gas(Gal[galaxy_number_].CentralGal,HotGasComponent,galaxy_number_,ColdGasComponent,mass_remaining_/Gal[galaxy_number_].ColdGas);

      //transfer reheated_mass_-mass_remaining_ from galaxy to the type 0
      if (reheated_mass_ > mass_remaining_)
        if(Gal[galaxy_number_].ColdGas > 0.) //if the reheat to itself, left cold gas below limit do not reheat to central
          transfer_gas(central_galaxy_number_,HotGasComponent,galaxy_number_,ColdGasComponent,(reheated_mass_-mass_remaining_)/Gal[galaxy_number_].ColdGas);
    }//types
  }//if(Gal[galaxy_number_].ColdGas > 0.)

  mass_checks("update_from_feedback #2",galaxy_number_);

  //DO EJECTION OF GAS
  if (Gal[Gal[galaxy_number_].CentralGal].HotGas > 0.)
  {
    if (ejected_mass_ > Gal[Gal[galaxy_number_].CentralGal].HotGas)
      ejected_mass_ = Gal[Gal[galaxy_number_].CentralGal].HotGas;  //either eject own gas or merger_centre gas for ttype 2's

    fraction_=((float)ejected_mass_) / Gal[Gal[galaxy_number_].CentralGal].HotGas;

    if (Gal[Gal[galaxy_number_].CentralGal].Type == 1)
    {
      /* If type 1, or type 2 orbiting type 1 near type 0 */
      if (FateOfSatellitesGas == 0)
      {   transfer_gas(Gal[galaxy_number_].CentralGal, EjectedGasComponent, Gal[galaxy_number_].CentralGal, HotGasComponent, fraction_); }
      else if (FateOfSatellitesGas == 1)
      {
        if (separation_ < Gal[central_galaxy_number_].Rvir)
        { transfer_gas(central_galaxy_number_       , HotGasComponent    , Gal[galaxy_number_].CentralGal, HotGasComponent, fraction_); }
        else
        { transfer_gas(Gal[galaxy_number_].CentralGal, EjectedGasComponent, Gal[galaxy_number_].CentralGal, HotGasComponent, fraction_); }
      }
    }
    else // If galaxy type 0 or type 2 merging into type 0
    { transfer_gas(central_galaxy_number_,EjectedGasComponent,Gal[galaxy_number_].CentralGal,HotGasComponent,fraction_); }
  }//(Gal[Gal[galaxy_number_].CentralGal].HotGas > 0.)
}


//Age in Mpc/Km/s/h - code units
void update_massweightage(const int galaxy_number_, const double stars_, const double time_)
{
  int output_number_;
  double age_;

  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    age_ = time_ - NumToTime(ListOutputSnaps[output_number_]);
#ifdef DETAILED_METALS_AND_MASS_RETURN
                 Gal[galaxy_number_].MassWeightAge[output_number_] += age_ * stars_;
#else
                 Gal[galaxy_number_].MassWeightAge[output_number_] += age_ * stars_ * (1. - RecycleFraction);
#endif
  }
}


/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria */
void check_disk_instability(const int galaxy_number_)
{
  double M_crit_, fraction_, stars_, disk_mass_;

/** @brief Calculates the stability of the stellar disk as discussed
 *         in Mo, Mao & White (1998). For unstable stars_, the required
 *         amount is transfered to the bulge to make the disk stable again.
 *         Mass, metals and luminosities updated. After Guo2010 the bulge
 *         size is followed and needs to be updated.
 *         Eq 34 & 35 in Guo2010 are used.
 */
//  if (DiskInstabilityModel == 1)
//  { return; }
//  else
  if(DiskInstabilityModel == 0)
  {
    /* check stellar disk -> eq 34 Guo2010*/
    M_crit_ = (Gal[galaxy_number_].Type == 0) ? Gal[galaxy_number_].Vmax * Gal[galaxy_number_].Vmax * Gal[galaxy_number_].StellarDiskRadius / Gravity : Gal[galaxy_number_].InfallVmax * Gal[galaxy_number_].InfallVmax * Gal[galaxy_number_].StellarDiskRadius / Gravity; 
    disk_mass_ = Gal[galaxy_number_].DiskMass;
    stars_ = disk_mass_ - M_crit_;

    /* add excess stars_ to the bulge */
    if(stars_ > 0.0)
    {
      fraction_ = stars_ / disk_mass_;
      
      /* to calculate the bulge size */
      update_bulge_from_disk(galaxy_number_, stars_);
      transfer_stars(galaxy_number_, BulgeComponent, galaxy_number_, DiskComponent, fraction_);

      if((BHGrowthInDiskInstabilityModel == 1) && (Gal[galaxy_number_].ColdGas > 0.))
      {
        Gal[galaxy_number_].BlackHoleMass += Gal[galaxy_number_].ColdGas * fraction_;
        Gal[galaxy_number_].ColdGas       -= Gal[galaxy_number_].ColdGas * fraction_;
      }

#ifndef POST_PROCESS_MAGS
      double Lumdisk;
      int output_number_, j;
#ifdef OUTPUT_REST_MAGS
      for(output_number_ = 0; output_number_ < NOUT; output_number_++)
      {
        for(j = 0; j < NMAG; j++)
        {
          Lumdisk = Gal[galaxy_number_].Lum[j][output_number_]-Gal[galaxy_number_].LumBulge[j][output_number_];
          Gal[galaxy_number_].LumBulge[j][output_number_] += fraction_ * Lumdisk;
          Lumdisk = Gal[galaxy_number_].YLum[j][output_number_]-Gal[galaxy_number_].YLumBulge[j][output_number_];
          Gal[galaxy_number_].YLumBulge[j][output_number_] += fraction_ * Lumdisk;
        }
      }
#endif
#ifdef COMPUTE_OBS_MAGS
      for(output_number_ = 0; output_number_ < NOUT; output_number_++)
      {
        for(j = 0; j < NMAG; j++)
        {
          Lumdisk = Gal[galaxy_number_].ObsLum[j][output_number_]-Gal[galaxy_number_].ObsLumBulge[j][output_number_];
          Gal[galaxy_number_].ObsLumBulge[j][output_number_] += fraction_ * Lumdisk;
          Lumdisk = Gal[galaxy_number_].ObsYLum[j][output_number_]-Gal[galaxy_number_].ObsYLumBulge[j][output_number_];
          Gal[galaxy_number_].ObsYLumBulge[j][output_number_] += fraction_ * Lumdisk;
#ifdef OUTPUT_MOMAF_INPUTS
          Lumdisk = Gal[galaxy_number_].dObsLum[j][output_number_]-Gal[galaxy_number_].dObsLumBulge[j][output_number_];
          Gal[galaxy_number_].dObsLumBulge[j][output_number_] += fraction_ * Lumdisk;
          Lumdisk = Gal[galaxy_number_].dObsYLum[j][output_number_]-Gal[galaxy_number_].dObsYLumBulge[j][output_number_];
          Gal[galaxy_number_].dObsYLumBulge[j][output_number_] += fraction_ * Lumdisk;
#endif
        }
      }
#endif
#endif

//       // debugging:
//       if ((Gal[galaxy_number_].BulgeMass > 1e-9 && Gal[galaxy_number_].BulgeSize == 0.0) ||
//           (Gal[galaxy_number_].BulgeMass == 0.0 && Gal[galaxy_number_].BulgeSize >1e-9))
//       { terminate( "bulgesize wrong in diskinstablility.c \n"); }
    }
  }
}


/** @brief Updates bulge from disk instability -> stars_ represents the mass
  *        transfered to the bulge, which occupies a size in the bulge equal
  *        to the occupied in the disk.
  *
  *  Introduced in Guo2010 to track the change in size of bulges
  *         after their growth due to disk instabilities. */
void update_bulge_from_disk(const int galaxy_number_, const double stars_)
{      
  const double original_bulge_size_  =  Gal[galaxy_number_].BulgeSize; //remove, not used

  /* alpha_inter=2.0/C=0.5 (alpha larger than in mergers since
   * the existing and newly formed bulges are concentric)*/
  const double f_int_ = 4.0;

  /* update the stellardisk spin due to the angular momentum transfer
   * from disk to bulge changing the specific angular momentum for disk stars_.
   * This should be done on the main routine, as this is update bulge.*/
  const double mass_fraction_ = stars_ /  Gal[galaxy_number_].DiskMass;
  
  if(mass_fraction_ >= 1) //everything transferred to the bulge
  { 
    Gal[galaxy_number_].StellarSpin[0] = 0;  
    Gal[galaxy_number_].StellarSpin[1] = 0;  
    Gal[galaxy_number_].StellarSpin[2] = 0;
  }
  else
  {
    const double speedup_ = 1. / (1. - mass_fraction_);
    Gal[galaxy_number_].StellarSpin[0] *= speedup_;
    Gal[galaxy_number_].StellarSpin[1] *= speedup_;
    Gal[galaxy_number_].StellarSpin[2] *= speedup_;
  }

  /* update disksize done, disk mass is automatically given by total-bulge*/

//GET BULGE SIZE - Eq. 35 in Guo2010
  /* if previous Bulge Mass = 0
     -> bulge size is given directly from newly formed bulge */
  const double bulgesize = bulge_from_disk(mass_fraction_) * Gal[galaxy_number_].StellarDiskRadius/3.;
  if(Gal[galaxy_number_].BulgeMass <1.e-9) 
  {
    /* size of newly formed bulge, which consists of the stellar mass
     * transfered from the disk. This is calculated using bulge_from_disk
     * which receives Delta_M/DiskMass and returns Rb/Rd. From eq 35 and
     * since DiskMass=2PISigma(Rd)^2 we see that Delta_M/DiskMass=1-(1+Rb/Rd)*exp(-Rb/Rd),
     * so function bulge_from_disk avoids calculating the slow "ln" function */
    Gal[galaxy_number_].BulgeSize=bulgesize;
  }      
  else
  {
    /* combine the old with newly formed bulge and calculate the
     * bulge size assuming energy conservation as for mergers but
     * using alpha=2. - eq 33 */
    Gal[galaxy_number_].BulgeSize=(Gal[galaxy_number_].BulgeMass+stars_)*(Gal[galaxy_number_].BulgeMass+stars_)/
      (Gal[galaxy_number_].BulgeMass*Gal[galaxy_number_].BulgeMass/Gal[galaxy_number_].BulgeSize+stars_*stars_/bulgesize+f_int_*Gal[galaxy_number_].BulgeMass*stars_/(Gal[galaxy_number_].BulgeSize+bulgesize));
  }

  if((Gal[galaxy_number_].BulgeMass + stars_ > 1.e-9 && Gal[galaxy_number_].BulgeSize == 0.0) || 
     (Gal[galaxy_number_].BulgeMass + stars_ == 0    && Gal[galaxy_number_].BulgeSize > 1.e-9))
  {
    printf("\nerror:\nGasDiskMass=%e GasDiskSize=%e \nStellarDiskMass=%e StellarDiskSize=%e \nBulgeMass=%e Bulgesize=%e\n",
                    Gal[galaxy_number_].ColdGas, Gal[galaxy_number_].GasDiskRadius, Gal[galaxy_number_].DiskMass, Gal[galaxy_number_].StellarDiskRadius,
                    Gal[galaxy_number_].BulgeMass, Gal[galaxy_number_].BulgeSize);
    printf("TransferSize=%e, OriBulgeSize=%e\n", bulgesize, original_bulge_size_);
    terminate("bulgesize or mass wrong in disk instablility");
  }
}


/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. */
double bulge_from_disk(const double frac)
{
  double x1,x2,x0,value;
/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. The bulge is assumed
 *         to form with the same size. avoid doing "ln" from eq 35*/
  x1=0.0;
  x2=1.;
  while ((func_size(x2,frac) * func_size(x1,frac))>0) 
  {
    x1=x2;
    x2=x2*2;
  }
  x0=x1+(x2-x1)/2.;
  value=func_size(x0,frac);
  if (value < 0) 
    value = -value;

  while(value>0.00001)
  {
    if(func_size(x0,frac)*func_size(x2,frac)>0)
      x2=x0;
    else
      x1=x0;
    x0=x1+(x2-x1)/2.;
    value=func_size(x0,frac);
    if (value < 0) 
      value = -value;
  }
    
  return x0;
}


double func_size(const double x, const double a)
{  return  exp(-x)*(1+x)-(1-a);}  





