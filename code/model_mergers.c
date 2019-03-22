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

/** @file model_mergers.c
 *  @brief Calculates the merging time_, the central galaxy (for type_ 1's),
 *         adds galaxies together, calculates SF from bursts and grows
 *         black holes.
 *
 *
 *       <B>set_merger_center</B> - calculates the central galaxy for type_ 1's,
 *       since type_ 1's can also merge. Therefore,
 *       they need a merger central galaxy and will also have a merger clock
 *       (needed for millennium two, since due to the high resolution, haloes
 *       are very difficult to disrupt and orbit around forever).
 *
 *
 *
 *       <B>estimate_merging_time</B> sets up a merger clock. Originally this
 *       was done only for type_ 2 galaxies. The positions of galaxies in the
 *       model are normally given by the position of their dark matter halo.
 *       However, when galaxies become satellites, their dark matter haloes
 *       are stripped by the central object to the point where there is none
 *       left. At this point, the dark matter of the satellite becomes part of
 *       the main halo, but the galaxy's position should continue to evolve
 *       due to the dynamical friction force caused by the dark matter around
 *       it. \n
 *       The code does keep track of the position of the most bounded particle
 *       when the satellite's halo was disrupted, but this is not used to track
 *       the galaxy position. Instead a clock is set, giving the time_ left until
 *       the satellite mergers with the central galaxy. Before, this was only
 *       done for type_ 2's (satellites that lost a halo). Guo2010 included the
 *       option to set the clock also for type_ 1's (satellites with
 *       a halo), since for the highest resolution millennium 2, their haloes
 *       can be very small and orbit forever around the central companion.\n
 *       This time_ is computed using the Chandrasekhar's formula for dynamical
 *       friction, as in Binney & Tremaine 1987:
 *
 *       \f$F_{\rm{df}}=
 *               -\frac{4\pi {\rm{G}}^2 m^2_{\rm{sat}} \ln(\Lambda) \rho B(x)}
 *                     {v^2_{\rm{rel}}}.\f$
 *
 *       Which gives (B&T Eq. 7.26):
 *
 *       \f$t_{df}\approx1.17\frac{V_{\rm{vir}}r^2_{\rm{sat}}}
 *                                {{\rm{G}}\,m_{\rm{sat}}\ln(\Lambda)},\f$
 *
 *       that is afterwards multiplied by 2 (after Delucia2007 to fit the
 *       data). When the merging time_ reaches zero, the satellite is assumed to
 *       merge with the central galaxy.
 *
 *
 *
 *       <B>deal_with_galaxy_merger</B> deals with the process, according to
 *       the mass fraction of the merger. Major if
 *       \f$M_{\rm{sat}}/M_{\rm{central}}>0.3\f$ and minor otherwise. It calls
 *       - add_galaxies_together - Add the cold and stellar phase of the merged
 *       galaxy to the central one. Also form a bulge at the central galaxy
 *       with the stars from the satellite in a minor merger if
 *       BulgeFormationInMinorMergersOn=1 (Major mergers are dealt later).
 *       - Then calls grow_black_hole - Grows black hole through accretion from
 *       cold gas during mergers (due to the instabilities triggered), as in
 *       Kauffmann & Haehnelt (2000). This is commonly referred as the quasar
 *       mode, main responsible for the black hole growth. After Croton2006 this
 *       mode is active even in minor mergers:
 *       \f$\Delta m_{\rm BH,Q}=M_{\rm{BH,min}}
 *          \frac{f_{\rm BH}(m_{\rm sat}/m_{\rm central})\,m_{\rm cold}}
 *           {1+(280\,\mathrm{km\,s}^{-1}/V_{\rm vir})^2}.\f$

 *       - Finally the burst of star formation due to the merger is treated.
 *           - If StarBurstModel = 0 (since Croton2006), the Somerville 2001
 *           model of bursts is used collisional_starburst_recipe(). The burst
 *           can happen for both major and minor mergers, with a fraction of
 *           the added cold gas from the satellite and central being consumed:
 *           \f$\dot{m}_{\star}^{\rm{burst}}
 *             = 0.56 \left(\frac{m_{\rm{sat}}}{m_{\rm{central}}}\right)^{0.7}
 *               m_{\rm{gas}}\f$.
 *           SN Feedback from starformation is computed and the sizes of bulge
 *           and disk followed.
 *
 *       - When a major merger occurs, the disk of both merging galaxies is
 *       completely destroyed to form a bulge. In either type_ of mergers, the
 *       bulge size is updated using Eq. 33 in Guo2010:
 *       \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
 *          C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2}+\alpha_{\rm{inter}}
 *          \frac{GM_1M_2}{R_1+R_2}\f$*/

/** @brief Calculates the central galaxies for type_ 1's. */

int get_merger_center(const int halo_number_)
{
  /** @brief Get id of central galaxy, since type 1's
   *         can have their merger clock started before they become type 2's
   *         if M_star>M_vir. Introduced for millennium 2, where they can have
   *         very small masses, still be followed and never merge. At this
   *         moment the centre is still not known, reason why this function is
   *         needed. Also, if the type 1 merges, all the type 2's associated with
   *         it will need to know the "new" central galaxy they are merging into.
   *
   * @note now assumes that halos with at least one galaxy also have at least one
   *       type-0 or type-1 galaxy
   *
   * @bug  possible bug: the branch returning 0 seems wrong to depend on i_ == 0
   */

  int same_fof_halo_number_, progenitor_halo_number_, first_occupied_progenitor_halo_number_, type_;
  double most_massive_halo_length_;

  //loop on all the haloes in current FOF group - to find a merger centre
  int current_galaxy_number_ = -1;
  int i_ = 0;
  int n_galaxies_ = 0;
  for(same_fof_halo_number_ = halo_number_; same_fof_halo_number_ >= 0; same_fof_halo_number_ = Halo[same_fof_halo_number_].NextHaloInFOFgroup)
  {
    most_massive_halo_length_ = 0;
    first_occupied_progenitor_halo_number_ = Halo[same_fof_halo_number_].FirstProgenitor;
    progenitor_halo_number_ = Halo[same_fof_halo_number_].FirstProgenitor;

    /* If the main progenitor of the current halo had no galaxies,
     * set first_occupied_progenitor_halo_number_ to the most massive progenitor. */
    if(progenitor_halo_number_ >= 0 && HaloAux[progenitor_halo_number_].NGalaxies == 0)
    {
      for(; progenitor_halo_number_ >= 0; progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor)
      {
        // for(i_ = 0, current_galaxy_number_ = HaloAux[progenitor_halo_number_].FirstGalaxy; i_ < HaloAux[progenitor_halo_number_].NGalaxies; ++i_, current_galaxy_number_ = HaloGal[current_galaxy_number_].NextGalaxy)
        // {
        //   type_ = HaloGal[current_galaxy_number_].Type;
        //   if((type_ == 0 || type_ == 1) && (Halo[progenitor_halo_number_].Len > most_massive_halo_length_))
        //   {
        //     most_massive_halo_length_ = Halo[progenitor_halo_number_].Len;
        //     first_occupied_progenitor_halo_number_ = progenitor_halo_number_;
        //   }
        // }
        if(HaloAux[progenitor_halo_number_].NGalaxies > 0 && Halo[progenitor_halo_number_].Len > most_massive_halo_length_)
        {
          most_massive_halo_length_ = Halo[progenitor_halo_number_].Len;
          first_occupied_progenitor_halo_number_ = progenitor_halo_number_;          
        }
      }
    }

    for(progenitor_halo_number_ = Halo[same_fof_halo_number_].FirstProgenitor; progenitor_halo_number_ >= 0; progenitor_halo_number_ = Halo[progenitor_halo_number_].NextProgenitor)//loop over all the progenitors
    {
      for(i_ = 0, current_galaxy_number_ = HaloAux[progenitor_halo_number_].FirstGalaxy; i_ < HaloAux[progenitor_halo_number_].NGalaxies; ++i_, current_galaxy_number_ = HaloGal[current_galaxy_number_].NextGalaxy)//loop over all the galaxies in a given progenitor
      {
        n_galaxies_++;
        
        type_ = HaloGal[current_galaxy_number_].Type;

        if(type_ == 0 || type_ == 1) // the galaxy is a type_ 0 or 1?
          if(progenitor_halo_number_ == first_occupied_progenitor_halo_number_) //is the main progenitor?
            if(same_fof_halo_number_ == Halo[same_fof_halo_number_].FirstHaloInFOFgroup) //is the main halo?
              return current_galaxy_number_;
      }
    }

  //if the halo has no galaxies, return 0
  if(i_ == 0)
    if(Halo[same_fof_halo_number_].FirstHaloInFOFgroup == same_fof_halo_number_)
      return 0;
  }
  // //if the halo has no galaxies, return 0
  // if(n_galaxies_ == 0)
  //   return 0;

  char error_message_[1000];
  sprintf(error_message_, "wrong in finding the central fof %d, gal %d, snapshot %d, n_galaxies_ %d\n", halo_number_, current_galaxy_number_, Halo[halo_number_].SnapNum, n_galaxies_);
  terminate(error_message_);
  return -1;
}


/** @brief Calculates the merging time_ whenever a galaxy becomes a satellite*/
double estimate_merging_time(const int halo_number_, const int mother_halo_number_, const int galaxy_number_)
{
  int central_halo_number_;
  double coulomb_, merging_time_, satellite_mass_, satellite_radius_, mother_halo_Rvir_;

  /** @brief Binney & Tremaine 1987 - 7.26 merging time_ for satellites due to
   *         dynamical friction. After Delucia2007 *2, shown to agree with
   *         Kolchin2008 simulations in Delucia2010. This is set when a galaxy
   *         becomes a type_ 2 or being a type_ 1 \f$M_{\rm{star}}>M_{\rm{vir}}\f$.
   *         In DeLucia2007 they could only merge into a type_ 0, now (after
   *         guo2010) they can merge into a type_ 1. */


  /*  recipe updated for more accurate merging time_ (see BT eq 7.26),
     now satellite radius at previous timestep is included */
  central_halo_number_ = Halo[Halo[halo_number_].Descendant].FirstProgenitor;
  if(Gal[galaxy_number_].Type == 1)
    central_halo_number_ = mother_halo_number_;
  if(central_halo_number_ == halo_number_)
    {
      terminate("can't be...!\n");
    }

  coulomb_ = log(Halo[mother_halo_number_].Len / ((double) Halo[halo_number_].Len) + 1);

  /*  should include stellar+cold gas in satellite_mass_! */
  satellite_mass_ = get_virial_mass(halo_number_)+(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass);

  satellite_radius_ = separation_halo(central_halo_number_,halo_number_)/(1 + ZZ[Halo[halo_number_].SnapNum]);

  int j;
  for (j = 0; j < 3; j++)
   Gal[galaxy_number_].DistanceToCentralGal[j] =  wrap(Halo[central_halo_number_].Pos[j] - Halo[halo_number_].Pos[j], BoxSize);


  mother_halo_Rvir_ = get_virial_radius(mother_halo_number_);
  if(satellite_radius_ > mother_halo_Rvir_)
    satellite_radius_ = mother_halo_Rvir_;

  if(satellite_mass_ > 0.0) 
  {
    merging_time_ = 1.17 * satellite_radius_ * satellite_radius_ * get_virial_velocity(mother_halo_number_) / (coulomb_ * Gravity * satellite_mass_); // Binney & Tremaine Eq.7.26

    /* change introduced by Delucia2007 to fit observations */
    merging_time_ = MergerTimeMultiplier * merging_time_;
    //merging_time_ = 2.*merging_time_;
  }
  else
    merging_time_ = -99999.9;

  return merging_time_;
}


/** @brief Deals with all the physics triggered by mergers */
void deal_with_galaxy_merger(const int galaxy_number_, const int merger_centralgal_, const int centralgal_, const double time_, const double deltaT_)
{

/** @brief Deals with the physics triggered by mergers, according to the mass
 *         fraction of the merger \f$(M_{\rm{sat}}/M_{\rm{central}}><0.3)\f$.
 *         Add the cold and stellar phase of the satellite galaxy to the central
 *         one, form a bulge at the central galaxy with the stars from the
 *         satellite in a minor merger if BulgeFormationInMinorMergersOn=1.
 *         Grows black hole through accretion from cold gas "quasar mode".
 *         If StarBurstModel = 0, the Somerville 2001 model
 *         of bursts is used, SN Feedback from starformation is computed and
 *         the sizes of bulge and disk followed. When a major merger occurs,
 *         the disk of both merging galaxies is completely destroyed to form
 *         a bulge. New stars form of to the bulge*/


#ifdef GALAXYTREE
  mass_checks("deal_with_galaxy_merger #0",galaxy_number_);
  mass_checks("deal_with_galaxy_merger #0",merger_centralgal_);
  mass_checks("deal_with_galaxy_merger #0",centralgal_);

  int q = Gal[merger_centralgal_].FirstProgGal;
  if(q >= 0)
  {
    while(GalTree[q].NextProgGal >= 0)
    { q = GalTree[q].NextProgGal; }

    GalTree[q].NextProgGal = Gal[galaxy_number_].FirstProgGal;

    if(GalTree[q].NextProgGal >= NGalTree)
    {
      printf("q=%d galaxy_number_=%d GalTree[q].NextProgGal=%d NGalTree=%d\n",
             q, galaxy_number_, GalTree[q].NextProgGal, NGalTree);
      terminate("problem walking GalTree");
    }
  }

  if(q < 0)
  { terminate("may not happen"); }

  q = GalTree[q].NextProgGal;

  if(HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q)
  { terminate("inconsistency between HaloGal and GalTree indices"); }
 
  HaloGal[GalTree[q].HaloGalIndex].MergeOn = 2;

  if(Gal[galaxy_number_].Type == 1)
  { HaloGal[GalTree[q].HaloGalIndex].MergeOn = 3; }
#endif

  /* flag galaxy as finished */
  Gal[galaxy_number_].Type = 3;

  /*  calculate mass ratio of merging galaxies */
  const double mi         = Gal[galaxy_number_                ].DiskMass + Gal[galaxy_number_                ].BulgeMass + Gal[galaxy_number_                ].ColdGas;
  const double ma         = Gal[merger_centralgal_].DiskMass + Gal[merger_centralgal_].BulgeMass + Gal[merger_centralgal_].ColdGas;
  const double mass_ratio = (mi > 0 || ma > 0) ? ((mi < ma) ? (mi / ma) : (ma / mi)) : 1.0;

  /* record the gas and stellar component  mass of merger central and satellite
   * galaxies before the merger */
  const double Mcstar=(Gal[merger_centralgal_].DiskMass+Gal[merger_centralgal_].BulgeMass);
  const double Mcbulge=Gal[merger_centralgal_].BulgeMass;
  const double Mcgas=Gal[merger_centralgal_].ColdGas;
  const double Mpstar=(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass);
  const double Mpbulge=Gal[galaxy_number_].BulgeMass;
  const double Mpgas=Gal[galaxy_number_].ColdGas;

  mass_checks("deal_with_galaxy_merger #1",galaxy_number_);
  mass_checks("deal_with_galaxy_merger #1",merger_centralgal_);
  mass_checks("deal_with_galaxy_merger #1",centralgal_);

  /* Add the cold and stellar phase of the merged galaxy to the central one.
     Also form a bulge if BulgeFormationInMinorMergersOn is set on (transfer stars
     from satellite disk to central bulge). In a major merger (dealt at the
     make_bulge_from_burst) the disk of the central (now made up of central and
     satellite will be moved to the bulge). Any new stars formed will go to the bulge */
  add_galaxies_together(merger_centralgal_, galaxy_number_);

  mass_checks("deal_with_galaxy_merger #2",galaxy_number_);
  mass_checks("deal_with_galaxy_merger #2",merger_centralgal_);
  mass_checks("deal_with_galaxy_merger #2",centralgal_);

  /* grow black hole through accretion from cold disk during mergers, as in
   * Kauffmann & Haehnelt (2000) + minor mergers - Quasar Mode */
  if(AGNRadioModeModel != 5)
    grow_black_hole(merger_centralgal_, mass_ratio, deltaT_);

  mass_checks("deal_with_galaxy_merger #3",galaxy_number_);
  mass_checks("deal_with_galaxy_merger #3",merger_centralgal_);
  mass_checks("deal_with_galaxy_merger #3",centralgal_);

  if (StarBurstModel == 0)
  {
  /* Starburst as in Somerville 2001, with feedback computed inside. */
        /* All star formation happens in the disk, but in a major merger this will then
         * be destroyed with everything moved to the bulge. */
    const double frac = collisional_starburst_recipe(mass_ratio, merger_centralgal_, centralgal_, time_, deltaT_);
    bulgesize_from_merger(mass_ratio,merger_centralgal_,galaxy_number_, Mcstar, Mcbulge, Mcgas, Mpstar, Mpbulge, Mpgas, frac);

    mass_checks("deal_with_galaxy_merger #3.5",galaxy_number_);
    mass_checks("deal_with_galaxy_merger #3.5",merger_centralgal_);
    mass_checks("deal_with_galaxy_merger #3.5",centralgal_);

    if(mass_ratio > ThreshMajorMerger)
      make_bulge_from_burst(merger_centralgal_);
  }

  mass_checks("deal_with_galaxy_merger #4",galaxy_number_);
  mass_checks("deal_with_galaxy_merger #4",merger_centralgal_);
  mass_checks("deal_with_galaxy_merger #4",centralgal_);

  /* If we are in the presence of a minor merger, check disk stability (the disk
   * is completely destroyed in major mergers) */
  if(mass_ratio < ThreshMajorMerger && (Gal[merger_centralgal_].DiskMass+Gal[merger_centralgal_].BulgeMass) > 0.0)
  { check_disk_instability(merger_centralgal_); }

  if ((Gal[merger_centralgal_].BulgeMass > 1.e-6 && Gal[merger_centralgal_].BulgeSize == 0.0) ||
      (Gal[merger_centralgal_].BulgeMass == 0.0 && Gal[merger_centralgal_].BulgeSize >1.e-6))
  {
    char error_message_[1000];
    sprintf(error_message_, "central: stellarmass %f, bulgemass %f, bulgesize %f, coldgas %f,gasdisk %f,stellardisk %f \n",
                    (Gal[merger_centralgal_].DiskMass+Gal[merger_centralgal_].BulgeMass),Gal[merger_centralgal_].BulgeMass,
                    Gal[merger_centralgal_].BulgeSize,Gal[merger_centralgal_].ColdGas,Gal[merger_centralgal_].GasDiskRadius,
                    Gal[merger_centralgal_].StellarDiskRadius);
    terminate(error_message_);
  }

  if (DiskRadiusModel == 0) 
  {
    set_gas_disk_radius(merger_centralgal_);
    set_stellar_disk_radius(merger_centralgal_);
  }

  mass_checks("deal_with_galaxy_merger #5",galaxy_number_);
  mass_checks("deal_with_galaxy_merger #5",merger_centralgal_);
  mass_checks("deal_with_galaxy_merger #5",centralgal_);
}


/** @brief Grows black holes, through accretion from cold gas during mergers,
 *          as in Kauffmann & Haehnelt (2000) - Quasar Mode.  */

void grow_black_hole(const int merger_centralgal_, const double mass_ratio, const double deltaT_)
{
  double BHaccrete, fraction;

  /** @brief Grows black hole through accretion from cold gas during mergers,
   *         as in Kauffmann & Haehnelt (2000).
   *                              in main.c */

  if(Gal[merger_centralgal_].ColdGas > 0.0)
  {
    BHaccrete = BlackHoleGrowthRate * mass_ratio
        / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[merger_centralgal_].Vvir))) * Gal[merger_centralgal_].ColdGas;
    /* redshift dependent accretion, not published */
    /* BHaccrete = BlackHoleGrowthRate * (1.0 + ZZ[Halo[halo_number_].SnapNum]) * mass_ratio */

    /* cannot accrete more gas than is available! */
    if(BHaccrete > Gal[merger_centralgal_].ColdGas)
      BHaccrete = Gal[merger_centralgal_].ColdGas;

    fraction=BHaccrete/Gal[merger_centralgal_].ColdGas;

    Gal[merger_centralgal_].BlackHoleMass += BHaccrete;
    Gal[merger_centralgal_].QuasarAccretionRate += BHaccrete / deltaT_;

    Gal[merger_centralgal_].ColdGas -= BHaccrete;
    metals_add_fraction_to(&Gal[merger_centralgal_].MetalsColdGas, Gal[merger_centralgal_].MetalsColdGas,-fraction);

#ifdef INDIVIDUAL_ELEMENTS
    elements_add_fraction_to(&Gal[merger_centralgal_].ColdGas_elements, Gal[merger_centralgal_].ColdGas_elements,-fraction);
#endif
  }
}


/** @brief Adds all the components of the satellite galaxy into its
 *         central companion. */
void add_galaxies_together(const int t, const int galaxy_number_)
{
  /** @brief All the components of the satellite galaxy are added to the
   *         correspondent component of the central galaxy. Cold gas spin
   *         is updated and a bulge is formed at the central galaxy, with
   *         the stars of the satellite if  BulgeFormationInMinorMergersOn=1.
   *         In case of a major merger, everything that was put in the disk of
   *         the central galaxy will be moved into the bulge
   */
  int outputbin, j;
  float tspin[3],tmass,pmass;

  /* t central, galaxy_number_ satellite */

  mass_checks("add_galaxies_together #0",galaxy_number_);
  mass_checks("add_galaxies_together #0.1",t);

  /* angular momentum transfer between gas*/
  tmass= Gal[t].ColdGas;
  pmass= Gal[galaxy_number_].ColdGas;

  Gal[t].MergeSat +=(Gal[galaxy_number_].DiskMass+Gal[galaxy_number_].BulgeMass);
  Gal[galaxy_number_].MergeSat=0.;

  transfer_gas(t,ColdGasComponent,galaxy_number_,ColdGasComponent,1.);
  //transfer_gas(t,EjectedGasComponent,galaxy_number_,ColdGasComponent,1.);
  transfer_gas(t,HotGasComponent,galaxy_number_,HotGasComponent,1.);
  transfer_gas(t,EjectedGasComponent,galaxy_number_,EjectedGasComponent,1.);
#ifdef TRACK_BURST
    /* The whole burst component gets transferred */
  transfer_stars(t,BurstComponent,galaxy_number_,BurstComponent,1.);
#endif
  if(BulgeFormationInMinorMergersOn)
    transfer_stars(t,BulgeComponent,galaxy_number_,DiskComponent,1.);
  else
    transfer_stars(t,DiskComponent,galaxy_number_,DiskComponent,1.);
  transfer_stars(t,BulgeComponent,galaxy_number_,BulgeComponent,1.);
  transfer_stars(t,ICMComponent,galaxy_number_,ICMComponent,1.);

  Gal[t].BlackHoleMass += Gal[galaxy_number_].BlackHoleMass;
  Gal[galaxy_number_].BlackHoleMass=0.;
  Gal[t].StarMerge += Gal[galaxy_number_].StarMerge;
  Gal[galaxy_number_].StarMerge=0.;

  mass_checks("add_galaxies_together #1",galaxy_number_);
  mass_checks("add_galaxies_together #1.1",t);

  /*update the gas spin*/
  for(j=0;j<3;j++)
    tspin[j]=Gal[t].GasSpin[j]*tmass+Gal[t].HaloSpin[j]*pmass;
  if (Gal[t].ColdGas != 0)
    for (j=0;j<3;j++)
     Gal[t].GasSpin[j]=tspin[j]/(Gal[t].ColdGas);

  Gal[t].Sfr += Gal[galaxy_number_].Sfr;

  if(BulgeFormationInMinorMergersOn)
    Gal[t].SfrBulge += Gal[galaxy_number_].Sfr;

  Gal[t].QuasarAccretionRate += Gal[galaxy_number_].QuasarAccretionRate;
  Gal[t].RadioAccretionRate += Gal[galaxy_number_].RadioAccretionRate;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
           Gal[t].MassWeightAge[outputbin] += Gal[galaxy_number_].MassWeightAge[outputbin];

#ifndef  POST_PROCESS_MAGS

/* Add the luminosities of the satellite and central galaxy */
#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
  {
    for(j = 0; j < NMAG; j++)
    {
      Gal[t].Lum[j][outputbin] += Gal[galaxy_number_].Lum[j][outputbin];
      Gal[t].YLum[j][outputbin] += Gal[galaxy_number_].YLum[j][outputbin];
#ifdef ICL
      Gal[t].ICLLum[j][outputbin] += Gal[galaxy_number_].ICLLum[j][outputbin];
#endif
    }
    if(BulgeFormationInMinorMergersOn)
    {
      for(j = 0; j < NMAG; j++)
      {
        Gal[t].LumBulge[j][outputbin] += Gal[galaxy_number_].Lum[j][outputbin];
        Gal[t].YLumBulge[j][outputbin] += Gal[galaxy_number_].YLum[j][outputbin];
      }
    }
    else
    {
      for(j = 0; j < NMAG; j++)
      {
        Gal[t].LumBulge[j][outputbin]  += Gal[galaxy_number_].LumBulge[j][outputbin];
        Gal[t].YLumBulge[j][outputbin] += Gal[galaxy_number_].YLumBulge[j][outputbin];
      }
    }
  }
#endif // OUTPUT_REST_MAGS

#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(j = 0; j < NMAG; j++) {
      Gal[t].ObsLum[j][outputbin]   += Gal[galaxy_number_].ObsLum[j][outputbin];
      Gal[t].ObsYLum[j][outputbin]  += Gal[galaxy_number_].ObsYLum[j][outputbin];
#ifdef ICL
      Gal[t].ObsICL[j][outputbin]  += Gal[galaxy_number_].ObsICL[j][outputbin];
#endif

#ifdef OUTPUT_MOMAF_INPUTS
      Gal[t].dObsLum[j][outputbin] += Gal[galaxy_number_].dObsLum[j][outputbin];
      Gal[t].dObsYLum[j][outputbin] += Gal[galaxy_number_].dObsYLum[j][outputbin];
#ifdef ICL
      Gal[t].dObsICL[j][outputbin]  += Gal[galaxy_number_].dObsICL[j][outputbin];
#endif
#endif
    }
    if(BulgeFormationInMinorMergersOn) {
      for(j = 0; j < NMAG; j++) {
        Gal[t].ObsLumBulge[j][outputbin]   += Gal[galaxy_number_].ObsLum[j][outputbin];
        Gal[t].ObsYLumBulge[j][outputbin]  += Gal[galaxy_number_].ObsYLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
        Gal[t].dObsLumBulge[j][outputbin]  += Gal[galaxy_number_].dObsLum[j][outputbin];
        Gal[t].dObsYLumBulge[j][outputbin] += Gal[galaxy_number_].dObsYLum[j][outputbin];
#endif
      }
    }
    else
    {
      for(j = 0; j < NMAG; j++) {
        Gal[t].ObsLumBulge[j][outputbin]   += Gal[galaxy_number_].ObsLumBulge[j][outputbin];
        Gal[t].ObsYLumBulge[j][outputbin]  += Gal[galaxy_number_].ObsYLumBulge[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
        Gal[t].dObsLumBulge[j][outputbin]  += Gal[galaxy_number_].dObsLumBulge[j][outputbin];
        Gal[t].dObsYLumBulge[j][outputbin] += Gal[galaxy_number_].dObsYLumBulge[j][outputbin];
#endif
      }
    }
  }
#endif //COMPUTE_OBS_MAGS
#endif  //POST_PROCESS_MAGS
}


/** @brief In a major merger, both disks are destroyed and all the mass transferred
 *         to the bulge. The galaxies have already been merged, so all we need to do here
 *         is transfer stars from disk to bulge. */

void make_bulge_from_burst(const int galaxy_number_)
{
  /* generate bulge */
  transfer_stars(galaxy_number_,BulgeComponent,galaxy_number_,DiskComponent,1.);

  /*  update the star formation rate */
  Gal[galaxy_number_].SfrBulge  = Gal[galaxy_number_].Sfr;

#ifndef  POST_PROCESS_MAGS
        int outputbin, j;
  
#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[galaxy_number_].LumBulge[j][outputbin]  = Gal[galaxy_number_].Lum[j][outputbin];
      Gal[galaxy_number_].YLumBulge[j][outputbin] = Gal[galaxy_number_].YLum[j][outputbin];
    }
  }
#endif /* defined OUTPUT_REST_MAGS */
#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[galaxy_number_].ObsLumBulge[j][outputbin]   = Gal[galaxy_number_].ObsLum[j][outputbin];
      Gal[galaxy_number_].ObsYLumBulge[j][outputbin]  = Gal[galaxy_number_].ObsYLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
      Gal[galaxy_number_].dObsLumBulge[j][outputbin]  = Gal[galaxy_number_].dObsLum[j][outputbin];
      Gal[galaxy_number_].dObsYLumBulge[j][outputbin] = Gal[galaxy_number_].dObsYLum[j][outputbin];
#endif /* defined OUTPUT_MOMAF_INPUTS */
    }
  }
#endif /* defined COMPUTE_OBS_MAGS */
#endif /* not defined POST_PROCESS_MAGS */
}

/** @brief Merger burst recipe from Somerville 2001 (used after Croton2006) */

double collisional_starburst_recipe(const double mass_ratio, const int merger_centralgal_, const int centralgal_,
                                    const double time_, const double deltaT_)
{
  /** @brief If StarBurstModel = 0 (since Croton2006), the Somerville 2001
   *         model of bursts is used. The burst can happen for both major
   *         and minor mergers, with a fraction of the added cold gas from
   *         the satellite and central being consumed. SN Feedback from
   *         starformation is computed and the sizes of bulge and disk
   *         followed (not done for the other burst mode).*/

  /* This is the major and minor merger starburst recipe of Somerville 2001.
   * The coefficients in eburst are taken from TJ Cox's PhD thesis and should
   * be more accurate then previous. */
  const double Ggas = Gal[merger_centralgal_].ColdGas;
  if(Ggas > 0.)
  {
    /* the bursting fraction given the mass ratio */
    /* m_dot = 0.56*(m_sat/m_central)^0.7*m_gas */
    // const double eburst = SfrBurstEfficiency * pow(mass_ratio, SfrBurstSlope);
    //eburst = 0.56 * pow(mass_ratio, 0.7);
    double mstars = SfrBurstEfficiency * pow(mass_ratio, SfrBurstSlope) * Gal[merger_centralgal_].ColdGas;

    //otherwise there is another check inside SN_feedback
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
    if(mstars > Gal[merger_centralgal_].ColdGas)
    { mstars = Gal[merger_centralgal_].ColdGas; }
#endif /* defined FEEDBACK_COUPLED_WITH_MASS_RETURN */

    /*  update the star formation rate */
    Gal[merger_centralgal_].Sfr += mstars / deltaT_;

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
    /* Store the value of the metallicity of the cold phase when SF occurs.
     * Used to update luminosities below */
    const double metallicitySF = metals_total(Gal[merger_centralgal_].MetalsColdGas) / Gal[merger_centralgal_].ColdGas;
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

    update_stars_due_to_reheat(merger_centralgal_, &mstars);

    // update_from_star_formation can only be called
    // after SD_feeedback recipe since stars need to be re_set once the reheated mass is known
    // (star formation and feedback share the same fraction of cold gas)
    update_from_star_formation(merger_centralgal_, mstars, true); // true indicates starburst

    mass_checks("collisional_starburst_recipe #2",merger_centralgal_);
        
    update_massweightage(merger_centralgal_, mstars, time_);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
    SN_feedback(merger_centralgal_, centralgal_, mstars, ColdGasComponent);
#endif /* not defined FEEDBACK_COUPLED_WITH_MASS_RETURN */

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
    /*  update the luminosities due to the stars formed */
    add_to_luminosities(merger_centralgal_, mstars, time_, deltaT_ / STEPS, metallicitySF);
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

    return mstars/Ggas;
  }
  else
    return 0.0;
}


/** @brief Calculates the bulge size after a merger. */
void  bulgesize_from_merger(const double mass_ratio, const int merger_centralgal_, const int galaxy_number_,
                            const double Mcstar, const double Mcbulge, const double Mcgas,
                            const double Mpstar, const double Mpbulge, const double Mpgas, double frac)
{
  /** @brief For any type_ of merger calculates the new bulge size using
   *         Eq. 33 in Guo2010:
   *
   *         \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
   *             C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2}+                \
   *             alpha_{\rm{inter}}\frac{GM_1M_2}{R_1+R_2}\f$.
   *
   *         This implementation assumed that the new bulge occupies the
   *         same space as the components that formed it. */

  double Mc,Rc;
  double Mp,Rp;
  double fint,c;

  fint=0.5;
  c=0.5;

  /* calculate radius for the object that will form the new bulge - Rc and Rp */
  /* Minor Merger */
  if(mass_ratio < ThreshMajorMerger)
  {
    /* In a minor merger only the stars of the satellite galaxy are moved
     * to the bulge of the central galaxy, therefore only stellar
     * components are used to compute radius and masses. */
    frac=0.0;
    /* in a minor merger only consider the bulge mass of the central galaxy */
    Mc=Mcbulge;
    Rc=Gal[merger_centralgal_].BulgeSize;
    /* and stellarmass of satellite*/
    Mp=Mpstar;
    if (Mp >0.0)
      Rp=(Gal[galaxy_number_].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[galaxy_number_].BulgeSize*Mpbulge)/Mpstar;
    else
      Rp=0.0;
  }
  /* Major Merger */
  else 
  {
    /* on a major merger both stellar and gas (after a burst) components
     * from the final bulge and need to be considered */
    /* Mc = bulge mass + burst of central*/
    Mc=Mcstar+frac*Mcgas;
    if (Mc > 0.0)
      Rc=(Gal[merger_centralgal_].StellarDiskRadius/3.*1.68*(Mcstar-Mcbulge)+Gal[merger_centralgal_].BulgeSize*Mcbulge+Gal[merger_centralgal_].GasDiskRadius*frac*Mcgas/3.*1.68)/(Mcgas*frac+Mcstar);
    else
      Rc=0.0;
    /* and satellite Mp */
    Mp=Mpstar+frac*Mpgas;
    if (Mp > 0.0)
      Rp=(Gal[galaxy_number_].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[galaxy_number_].BulgeSize*Mpbulge+Gal[galaxy_number_].GasDiskRadius*frac*Mpgas/3.*1.68)/(Mpgas*frac+Mpstar);
    else
      Rp=0.0;
  }

  if(Rp>0. && Rp<1.e-8)
          Rp=1.e-8;
  if(Rc>0. && Rc<1.e-8)
            Rc=1.e-8;
  /* If both original radius are bigger then 0 then this is Eq. 33 in Guo 2010
   * solved for R_new,bulge with all terms divided by G and C. */
  if(Rc >= 1.e-8 && Rp >= 1.e-8)
          Gal[merger_centralgal_].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mp*Mp/Rp+Mc*Mc/Rc+fint/c*Mp*Mc/(Rp+Rc));

  if(Rc >= 1.e-8 && Rp <= 1.e-8)
    Gal[merger_centralgal_].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mc*Mc/Rc+fint/c*Mp*Mc/(Rp+Rc));

  if(Rc <= 1.e-8 && Rp <=1.e-8)
    Gal[merger_centralgal_].BulgeSize=0.0;

  if(Rc <= 1.e-8 && Rp >= 1.e-8)
    Gal[merger_centralgal_].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mp*Mp/Rp+fint/c*Mp*Mc/(Rp+Rc));


  if ((Mp+Mc > 0.0 && Gal[merger_centralgal_].BulgeSize == 0.0 )||(Mp+Mc == 0.0 && Gal[merger_centralgal_].BulgeSize> 0.0))
    {
      char error_message_[1000];
      printf("halo_number_=%d, merger_centralgal_ %d\n\n", Gal[merger_centralgal_].HaloNr, merger_centralgal_);

      printf("New Bulge Mass from Central (MC)=%e\n New Bulge Mass from Satellite (Mp)=%e\n NewBulge size=%e\n\n",
             Mp,Mc,Gal[merger_centralgal_].BulgeSize);

      if(mass_ratio < ThreshMajorMerger)
        {
          printf("minor merger, new mass=original mass\n");
          printf("New Bulge Mass From Central (Mc)   = Mcbulge = %f\n",Mcbulge);
          printf("New Bulge Mass From Satellite (Mp) = Mpstar  = %f\n",Mpstar);
        }
      else
        {
          printf("New Bulge From Central (Mc)   = Mcstar+frac*Mcgas = %f+%f*%f\n", Mcstar, frac, Mcgas);
          printf("New Bulge From Satellite (Mp) = Mpstar+frac*Mpgas = %f+%f*%f\n", Mpstar, frac, Mpgas);
        }

      printf("BulgeSize from Central (Rc)=%e\nBulgeSize from Satellite (Rp)=%e\nmass ratio=%f\n\n",Rc,Rp, mass_ratio);

      printf("the following masses don't tell a lot because things have been merged already!!!\n");
      printf("    sat: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
             Gal[galaxy_number_].BulgeMass, Gal[galaxy_number_].BulgeSize, Gal[galaxy_number_].ColdGas, Gal[galaxy_number_].GasDiskRadius,
             Gal[galaxy_number_].DiskMass, Gal[galaxy_number_].StellarDiskRadius);
      printf(        "central: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
                Gal[merger_centralgal_].BulgeMass,  Gal[merger_centralgal_].BulgeSize, Gal[merger_centralgal_].ColdGas,
                Gal[merger_centralgal_].GasDiskRadius, Gal[merger_centralgal_].DiskMass, Gal[merger_centralgal_].StellarDiskRadius);

      sprintf(error_message_,"\n bulgesize wrong in merger");
      terminate(error_message_);
      exit(0);
  }

}
