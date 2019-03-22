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
 *  @brief recipe_starformation_and_feedback.c computes the amount of stars
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


/** @brief Main recipe, calculates the fraction of cold gas turned into stars due
  *        to star formation; the fraction of mass instantaneously recycled and
  *        returned to the cold gas; the fraction of gas reheated from cold to hot,
  *        ejected from hot to external and returned from ejected to hot due to
  *        SN feedback.   */
void starformation(const int p, const int centralgal, const double time, const double dt, const int nstep)
{
        /** Variables: reff-Rdisk, tdyn=Rdisk/Vmax, strdot=Mstar_dot, stars=strdot*dt */
  /* Note that Units of dynamical time are Mpc/Km/s - no conversion on dt needed */
  double tdyn, stars, cold_crit;
  
  if(Gal[p].Type == 0)
  {
    tdyn = Gal[p].GasDiskRadius / Gal[p].Vmax;
//     cold_crit = SfrColdCrit * Gal[p].Vmax/200. * Gal[p].GasDiskRadius*100.;
    cold_crit = 0.5 * SfrColdCrit * Gal[p].Vmax * Gal[p].GasDiskRadius;
  }
  else
  {
    tdyn = Gal[p].GasDiskRadius / Gal[p].InfallVmax;
//     cold_crit = SfrColdCrit * Gal[p].InfallVmax/200. * Gal[p].GasDiskRadius*100.;
    cold_crit = 0.5 * SfrColdCrit * Gal[p].InfallVmax * Gal[p].GasDiskRadius;
  }

  //standard star formation law (Croton2006, Delucia2007, Guo2010)
  /* if(StarFormationModel == 0) */
  { stars = (Gal[p].ColdGas > cold_crit) ? SfrEfficiency * dt * (Gal[p].ColdGas - cold_crit) / tdyn : 0.0; }
  /*else if(StarFormationModel == 1)
  {
          ALTERNATIVE STAR FORMATION LAW 
  }*/

  if (stars > 0.)
  {
//otherwise cold gas and stars share material in update_stars_due_to_reheat
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
    if(stars > Gal[p].ColdGas)
    { stars = Gal[p].ColdGas; }
#endif

    mass_checks("recipe_starform #1",p);
    mass_checks("recipe_starform #1.1",centralgal);

  /* update for star formation
   * updates Mcold, StellarMass, MetalsMcold and MetalsStellarMass
   * in Guo2010 case updates the stellar spin -> hardwired, not an option */

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
    /* Store the value of the metallicity of the cold phase when SF occurs */
    const double metallicitySF = (Gal[p].ColdGas > 0.) ? metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas : 0;
#endif //NDEF POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

    update_stars_due_to_reheat(p, centralgal, &stars);

    mass_checks("recipe_starform #2",p);
    mass_checks("recipe_starform #2.1",centralgal);

    /*  update the star formation rate */
    /*Sfr=stars/(dt*steps)=strdot*dt/(dt*steps)=strdot/steps -> average over the STEPS*/
    Gal[p].Sfr += stars / (dt * STEPS);

    // update_from_star_formation can only be called
    // after SD_feeedback recipe since stars need to be re_set once the reheated mass is known
    // (star formation and feedback share the same fraction of cold gas)
    update_from_star_formation(p, stars, false, nstep); // false indicates not a burst

    update_massweightage(p, stars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  /* ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN feedback is only called
   * when stars die, inside DETAILED_METALS_AND_MASS_RETURN */
    SN_feedback(p, centralgal, stars, ColdGasComponent);
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  /*  Update the luminosities due to the stars formed */
    add_to_luminosities(p, stars, time, dt, metallicitySF);
#endif //NDEF POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
  }

  if(Gal[p].DiskMass > 0.0)
  { check_disk_instability(p); }

  if (DiskRadiusModel== 0)
  { set_stellar_disk_radius(p); }
}


 /** @brief update mass of stars formed due to reheating 
  * 
  * @bug Which is the true central galaxy? Gal[p].CentralGal or Gal[centralgal]?
  *       Tests show, these are not always the same.
  */
void update_stars_due_to_reheat(const int p, const int centralgal, double *stars)
{
  (void) centralgal;
  
//   double MergeCentralVvir=0.;
//   double CentralVvir=0.;
//   doble ejected_mass=0.;
  /* SN FEEDBACK RECIPES */

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */

  //REHEAT
//   CentralVvir = Gal[centralgal].Vvir; // main halo Vvir
//   MergeCentralVvir = Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir

  // Feedback depends on the circular velocity of the host halo
  // Guo2010 - eq 18 & 19
  if(FeedbackReheatingModel == 0)
  {
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
    double reheated_mass;
    if (Gal[Gal[p].CentralGal].Type == 0)
    { reheated_mass = FeedbackReheatingEpsilon * (*stars) * (0.5 + pow(Gal[Gal[p].CentralGal].Vmax      / ReheatPreVelocity, -ReheatSlope)); }
    else
    { reheated_mass = FeedbackReheatingEpsilon * (*stars) * (0.5 + pow(Gal[Gal[p].CentralGal].InfallVmax / ReheatPreVelocity, -ReheatSlope)); }
    if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > (*stars) * (EtaSNcode * EnergySNcode))
    { reheated_mass = (*stars) * (EtaSNcode * EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir); }

    if((*stars + reheated_mass) > Gal[p].ColdGas)
    { *stars *= Gal[p].ColdGas / (*stars + reheated_mass); }
#endif
  }
  /*else if(FeedbackReheatingModel == 1)
  {
    reheated_mass =  ALTERNATIVE Reheating LAW ;
  }*/
}


/** @brief Updates the different components due to star formation: mass
  *        and metals in stars and cold gas and stellar spin. */
//void update_from_star_formation(int p, double time, double stars, double metallicity)
void update_from_star_formation(const int p, const double stars, const bool flag_burst, const int nstep)
{
  (void) nstep; /* avoid unused-parameter warning */

  if(Gal[p].ColdGas <= 0. || stars <= 0.)
  {
    printf("update_from_star_formation: Gal[p].ColdGas <= 0. || stars <= 0.\n");
    exit(0);
  }

  /* If DETAILED_METALS_AND_MASS_RETURN, no longer an assumed instantaneous
   * recycled fraction. Mass is returned over time via SNe and AGB winds.
   * Update the Stellar Spin when forming stars */
#ifndef DETAILED_METALS_AND_MASS_RETURN
  const double stars_to_add = (1 - RecycleFraction) * stars;
#else  /* not defined DETAILED_METALS_AND_MASS_RETURN */
  const double stars_to_add = stars;
#endif /* not defined DETAILED_METALS_AND_MASS_RETURN */

  if (Gal[p].DiskMass+stars_to_add > 1.e-8)
  {
    int i;
    for (i = 0; i < 3; i++)
    { Gal[p].StellarSpin[i]= (Gal[p].StellarSpin[i] * Gal[p].DiskMass + stars_to_add * Gal[p].GasSpin[i]) / (Gal[p].DiskMass + stars_to_add); }
  }
    /*  Update Gas and Metals from star formation */
  mass_checks("update_from_star_formation #0",p);

  const double fraction = stars_to_add / Gal[p].ColdGas;

#ifdef STAR_FORMATION_HISTORY
  Gal[p].sfh_DiskMass[Gal[p].sfh_ibin] += stars_to_add; //ROB: Now, all SF gas is put in SFH array ("recycled' mass will return to gas phase over time)
  metals_add_fraction_to(&Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin], Gal[p].MetalsColdGas, fraction);
#ifdef INDIVIDUAL_ELEMENTS
  elements_add_fraction_to(&Gal[p].sfh_ElementsDiskMass[Gal[p].sfh_ibin], Gal[p].ColdGas_elements, fraction);
#endif /* defined INDIVIDUAL_ELEMENTS */
#ifdef TRACK_BURST
  if (flag_burst)
  { Gal[p].sfh_BurstMass[Gal[p].sfh_ibin] += stars_to_add; }
#else  /* not defined TRACK_BURST */
  (void) flag_burst; /* avoid unused-parameter warning */
#endif /* not defined TRACK_BURST */
#else  /* not defined  STAR_FORMATION_HISTORY */
  (void) flag_burst; /* avoid unused-parameter warning */
#endif /* not defined STAR_FORMATION_HISTORY */

  metals_add_fraction_to(&Gal[p].MetalsDiskMass, Gal[p].MetalsColdGas,  fraction);
  metals_add_fraction_to(&Gal[p].MetalsColdGas , Gal[p].MetalsColdGas, -fraction);

  //GLOBAL PROPERTIES
  Gal[p].DiskMass += stars_to_add;
  Gal[p].ColdGas  -= stars_to_add;
#ifdef INDIVIDUAL_ELEMENTS
  elements_add_fraction_to(&Gal[p].DiskMass_elements, Gal[p].ColdGas_elements,  fraction);
  Gelements_add_fraction_to(&Gal[p].ColdGas_elements , Gal[p].ColdGas_elements, -fraction);
#endif
#ifdef TRACK_BURST
  if (flag_burst) Gal[p].BurstMass += stars_to_add;
#endif

  mass_checks("update_from_star_formation #1",p);

  /* Formation of new metals - instantaneous recycling approximation - only SNII
   * Also recompute the metallicity of the cold phase.*/
#ifndef DETAILED_METALS_AND_MASS_RETURN
  /* stars used because the Yield is defined as a fraction of
   * all stars formed, not just long lived */
  Gal[p].MetalsColdGas += Yield * stars;
#endif

  if (DiskRadiusModel == 0)
    set_stellar_disk_radius(p);
}


/* there are two modes for supernova feedback corresponding to when the mass returning
 * by dying stars is returned to the cold gas - reheat and ejection; and when the mass
 * is returned to the hot gas - onle ejection.*/
void SN_feedback(const int p, const int centralgal, const double stars, const GasComponentType feedback_location)
{
//   double CentralVvir, MergeCentralVvir=0.;
  double EjectVmax, EjectVvir, SN_Energy, Reheat_Energy;
  double reheated_mass=0., ejected_mass=0.;
  /* SN FEEDBACK MODEL */

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */
   
  if (feedback_location == HotGasComponent)
  { reheated_mass = 0.; }
  else if(feedback_location == ColdGasComponent)
  {
//   CentralVvir = Gal[centralgal].Vvir; // main halo Vvir
//   MergeCentralVvir = Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir

    mass_checks("recipe_starform #0",p);
    mass_checks("recipe_starform #0.1",centralgal);

    // Feedback depends on the circular velocity of the host halo
    // Guo2010 - eq 18 & 19
    /* if(FeedbackReheatingModel == 0) */
    {
      if (Gal[Gal[p].CentralGal].Type == 0)
      { reheated_mass = FeedbackReheatingEpsilon * stars * (0.5 + pow(Gal[Gal[p].CentralGal].Vmax       / ReheatPreVelocity, -ReheatSlope)); }
      else
      { reheated_mass = FeedbackReheatingEpsilon * stars * (0.5 + pow(Gal[Gal[p].CentralGal].InfallVmax / ReheatPreVelocity, -ReheatSlope)); }
      if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > stars * (EtaSNcode * EnergySNcode))
      { reheated_mass = stars * (EtaSNcode * EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir); }
    }
   /* else if(FeedbackReheatingModel == 1)
    {
      reheated_mass =  ALTERNATIVE Reheating LAW ;
    } */

    if(reheated_mass > Gal[p].ColdGas)
    { reheated_mass = Gal[p].ColdGas; }
  }// end if feedback_location

  /* Determine ejection (for FeedbackEjectionModel 2 we have the dependence on Vmax)
   * Guo2010 - eq 22
   * Note that satellites can now retain gas and have their own gas cycle*/
  if (Gal[Gal[p].CentralGal].Type == 0)
  {
    EjectVmax=Gal[centralgal].Vmax;
    EjectVvir=Gal[centralgal].Vvir;// main halo Vvir
  }
  else
  {
    EjectVmax=Gal[Gal[p].CentralGal].InfallVmax;
    EjectVvir=Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir
  }

  if(FeedbackEjectionModel == 0)
  {
    ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) * stars *
            min(1./FeedbackEjectionEfficiency, 0.5 + pow(EjectVmax / EjectPreVelocity, -EjectSlope)) -
            reheated_mass * EjectVvir * EjectVvir) / (EjectVvir * EjectVvir);
  }
  else if(FeedbackEjectionModel == 1)//the ejected material is assumed to have V_SN
  {
    SN_Energy = 0.5 * stars * (EtaSNcode * EnergySNcode);
    Reheat_Energy = 0.5 * reheated_mass * EjectVvir * EjectVvir;

    ejected_mass = (SN_Energy - Reheat_Energy)/(0.5 * FeedbackEjectionEfficiency*(EtaSNcode * EnergySNcode));

    //if VSN^2<Vvir^2 nothing is ejected
    if(FeedbackEjectionEfficiency*(EtaSNcode * EnergySNcode)<EjectVvir*EjectVvir)
    { ejected_mass = 0.0; }
  }

  // Finished calculating mass exchanges, so just check that none are negative
  if (reheated_mass < 0.0) reheated_mass = 0.0;
  if (ejected_mass < 0.0) ejected_mass = 0.0;

  /* Update For Feedback */
  /* update cold, hot, ejected gas fractions and respective metallicities
   * there are a number of changes introduced by Guo2010 concerning where
   * the gas ends up */

  //ejected_mass = 0.01*Gal[centralgal].HotGas;
  if (reheated_mass + ejected_mass > 0.)
    update_from_feedback(p, centralgal, reheated_mass, ejected_mass);
}


/** @brief Updates cold, hot and external gas components due to SN
 *         reheating and ejection. */
void update_from_feedback(const int p, const int centralgal, const double reheated_mass, double ejected_mass)
{
  double dis=0.;
  double massremain;
  double fraction;
//  int merger_centre;

  //mass_checks("update_from_feedback #1",p);

  if(Gal[p].ColdGas > 0.)
  {
    //REHEAT
    // if galaxy is a type 1 or a type 2 orbiting a type 1 with hot gas being striped,
    //some of the reheated and ejected masses goes to the type 0 and some stays in the type 1

    if(Gal[p].Type ==0)
    {
      transfer_gas(p,HotGasComponent,p,ColdGasComponent,((float)reheated_mass)/Gal[p].ColdGas);
    }
    else if(Gal[p].Type <3)
    {
//           if(Gal[p].Type ==1)
//             merger_centre=centralgal;
//           else if(Gal[p].Type ==2)
//             merger_centre=Gal[p].CentralGal;

      dis=separation_gal(centralgal,Gal[p].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

      //compute share of reheated mass
      if (dis<Gal[centralgal].Rvir && Gal[Gal[p].CentralGal].Type == 1)
      {
        //mass that remains on type1 (the rest goes to type 0) for reheat - massremain, for eject - ejected mass
        massremain=reheated_mass*Gal[p].HotRadius/Gal[p].Rvir;
        ejected_mass = ejected_mass*Gal[p].HotRadius/Gal[p].Rvir;

        if (massremain > reheated_mass)
        { massremain = reheated_mass; }
      }
      else
      { massremain=reheated_mass; }

      //needed due to precision issues, since we first remove massremain and then (reheated_mass-massremain)
      //from the satellite into the type 0 and type 1 the fraction might not add up on the second call
      //since Gal[p].ColdGas is a float and reheated_mass & massremain are doubles
      if((massremain + reheated_mass)>Gal[p].ColdGas)
        massremain=Gal[p].ColdGas-reheated_mass;

      //transfer massremain
      transfer_gas(Gal[p].CentralGal,HotGasComponent,p,ColdGasComponent,massremain/Gal[p].ColdGas);

      //transfer reheated_mass-massremain from galaxy to the type 0
      if (reheated_mass > massremain)
        if(Gal[p].ColdGas > 0.) //if the reheat to itself, left cold gas below limit do not reheat to central
          transfer_gas(centralgal,HotGasComponent,p,ColdGasComponent,(reheated_mass-massremain)/Gal[p].ColdGas);
    }//types
  }//if(Gal[p].ColdGas > 0.)

  mass_checks("update_from_feedback #2",p);

  //DO EJECTION OF GAS
  if (Gal[Gal[p].CentralGal].HotGas > 0.)
  {
    if (ejected_mass > Gal[Gal[p].CentralGal].HotGas)
      ejected_mass = Gal[Gal[p].CentralGal].HotGas;  //either eject own gas or merger_centre gas for ttype 2's

    fraction=((float)ejected_mass) / Gal[Gal[p].CentralGal].HotGas;

    if (Gal[Gal[p].CentralGal].Type == 1)
    {
      /* If type 1, or type 2 orbiting type 1 near type 0 */
      if (FateOfSatellitesGas == 0)
      {   transfer_gas(Gal[p].CentralGal, EjectedGasComponent, Gal[p].CentralGal, HotGasComponent, fraction); }
      else if (FateOfSatellitesGas == 1)
      {
        if (dis < Gal[centralgal].Rvir)
        { transfer_gas(centralgal       , HotGasComponent    , Gal[p].CentralGal, HotGasComponent, fraction); }
        else
        { transfer_gas(Gal[p].CentralGal, EjectedGasComponent, Gal[p].CentralGal, HotGasComponent, fraction); }
      }
    }
    else // If galaxy type 0 or type 2 merging into type 0
    { transfer_gas(centralgal,EjectedGasComponent,Gal[p].CentralGal,HotGasComponent,fraction); }
  }//(Gal[Gal[p].CentralGal].HotGas > 0.)
}


//Age in Mpc/Km/s/h - code units
void update_massweightage(const int p, const double stars, const double time)
{
  int outputbin;
  double age;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
  {
    age = time - NumToTime(ListOutputSnaps[outputbin]);
#ifdef DETAILED_METALS_AND_MASS_RETURN
                 Gal[p].MassWeightAge[outputbin] += age * stars;
#else
                 Gal[p].MassWeightAge[outputbin] += age * stars * (1. - RecycleFraction);
#endif
  }
}


/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria */
void check_disk_instability(const int p)
{
  double Mcrit, fraction, stars, diskmass;

/** @brief Calculates the stability of the stellar disk as discussed
 *         in Mo, Mao & White (1998). For unstable stars, the required
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
    Mcrit = (Gal[p].Type == 0) ? Gal[p].Vmax * Gal[p].Vmax * Gal[p].StellarDiskRadius / G : Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].StellarDiskRadius / G; 
    diskmass = Gal[p].DiskMass;
    stars = diskmass - Mcrit;

    /* add excess stars to the bulge */
    if(stars > 0.0)
    {
      fraction = stars / diskmass;
      
      /* to calculate the bulge size */
      update_bulge_from_disk(p, stars);
      transfer_stars(p, BulgeComponent, p, DiskComponent, fraction);

      if((BHGrowthInDiskInstabilityModel == 1) && (Gal[p].ColdGas > 0.))
      {
        Gal[p].BlackHoleMass += Gal[p].ColdGas * fraction;
        Gal[p].ColdGas       -= Gal[p].ColdGas * fraction;
      }

#ifndef POST_PROCESS_MAGS
      double Lumdisk;
      int outputbin, j;
#ifdef OUTPUT_REST_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
      {
        for(j = 0; j < NMAG; j++)
        {
          Lumdisk = Gal[p].Lum[j][outputbin]-Gal[p].LumBulge[j][outputbin];
          Gal[p].LumBulge[j][outputbin] += fraction * Lumdisk;
          Lumdisk = Gal[p].YLum[j][outputbin]-Gal[p].YLumBulge[j][outputbin];
          Gal[p].YLumBulge[j][outputbin] += fraction * Lumdisk;
        }
      }
#endif
#ifdef COMPUTE_OBS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
      {
        for(j = 0; j < NMAG; j++)
        {
          Lumdisk = Gal[p].ObsLum[j][outputbin]-Gal[p].ObsLumBulge[j][outputbin];
          Gal[p].ObsLumBulge[j][outputbin] += fraction * Lumdisk;
          Lumdisk = Gal[p].ObsYLum[j][outputbin]-Gal[p].ObsYLumBulge[j][outputbin];
          Gal[p].ObsYLumBulge[j][outputbin] += fraction * Lumdisk;
#ifdef OUTPUT_MOMAF_INPUTS
          Lumdisk = Gal[p].dObsLum[j][outputbin]-Gal[p].dObsLumBulge[j][outputbin];
          Gal[p].dObsLumBulge[j][outputbin] += fraction * Lumdisk;
          Lumdisk = Gal[p].dObsYLum[j][outputbin]-Gal[p].dObsYLumBulge[j][outputbin];
          Gal[p].dObsYLumBulge[j][outputbin] += fraction * Lumdisk;
#endif
        }
      }
#endif
#endif

//       // debugging:
//       if ((Gal[p].BulgeMass > 1e-9 && Gal[p].BulgeSize == 0.0) ||
//           (Gal[p].BulgeMass == 0.0 && Gal[p].BulgeSize >1e-9))
//       { terminate( "bulgesize wrong in diskinstablility.c \n"); }
    }
  }
}


/** @brief Updates bulge from disk instability -> stars represents the mass
  *        transfered to the bulge, which occupies a size in the bulge equal
  *        to the occupied in the disk.
  *
  *  Introduced in Guo2010 to track the change in size of bulges
  *         after their growth due to disk instabilities. */
void update_bulge_from_disk(const int p, const double stars)
{      
  const double orisize  =  Gal[p].BulgeSize; //remove, not used

  /* alpha_inter=2.0/C=0.5 (alpha larger than in mergers since
   * the existing and newly formed bulges are concentric)*/
  const double fint = 4.0;

  /* update the stellardisk spin due to the angular momentum transfer
   * from disk to bulge changing the specific angular momentum for disk stars.
   * This should be done on the main routine, as this is update bulge.*/
  const double massfrac = stars /  Gal[p].DiskMass;
  
  if(massfrac >= 1) //everything transferred to the bulge
  { 
    Gal[p].StellarSpin[0] = 0;  
    Gal[p].StellarSpin[1] = 0;  
    Gal[p].StellarSpin[2] = 0;
  }
  else
  {
    const double speedup = 1. / (1. - massfrac);
    Gal[p].StellarSpin[0] *= speedup;
    Gal[p].StellarSpin[1] *= speedup;
    Gal[p].StellarSpin[2] *= speedup;
  }

  /* update disksize done, disk mass is automatically given by total-bulge*/

//GET BULGE SIZE - Eq. 35 in Guo2010
  /* if previous Bulge Mass = 0
     -> bulge size is given directly from newly formed bulge */
  const double bulgesize = bulge_from_disk(massfrac) * Gal[p].StellarDiskRadius/3.;
  if(Gal[p].BulgeMass <1.e-9) 
  {
    /* size of newly formed bulge, which consists of the stellar mass
     * transfered from the disk. This is calculated using bulge_from_disk
     * which receives Delta_M/DiskMass and returns Rb/Rd. From eq 35 and
     * since DiskMass=2PISigma(Rd)^2 we see that Delta_M/DiskMass=1-(1+Rb/Rd)*exp(-Rb/Rd),
     * so function bulge_from_disk avoids calculating the slow "ln" function */
    Gal[p].BulgeSize=bulgesize;
  }      
  else
  {
    /* combine the old with newly formed bulge and calculate the
     * bulge size assuming energy conservation as for mergers but
     * using alpha=2. - eq 33 */
    Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars)/
      (Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize+stars*stars/bulgesize+fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));
  }

  if((Gal[p].BulgeMass + stars > 1.e-9 && Gal[p].BulgeSize == 0.0) || 
     (Gal[p].BulgeMass + stars == 0    && Gal[p].BulgeSize > 1.e-9))
  {
    printf("\nerror:\nGasDiskMass=%e GasDiskSize=%e \nStellarDiskMass=%e StellarDiskSize=%e \nBulgeMass=%e Bulgesize=%e\n",
                    Gal[p].ColdGas, Gal[p].GasDiskRadius, Gal[p].DiskMass, Gal[p].StellarDiskRadius,
                    Gal[p].BulgeMass, Gal[p].BulgeSize);
    printf("TransferSize=%e, OriBulgeSize=%e\n", bulgesize, orisize);
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





