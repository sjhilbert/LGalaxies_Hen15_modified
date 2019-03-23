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

/** @file   star_formation_history.c
 *  @date   2010-2019
 *  @author Peter Thomas (original code)
 *  @author Stefan Hilbert
 *
 *  @brief  Routines to track the star-formation history of galaxies.
 *          Keeps bins of the size of the current step number of the semi-analytic.
 *          As a galaxy ages, bins are merged in factors of two.
 *
 * This version assumes that all stars formed on a particular timestep
 * go into a single time history bin with star formation time_ equal to
 * the mid point of the bin. All times are code units (Mpc/Km/s/h) and
 * are the time_ till present.
 *
 * Number of time_ bins needed:
 * #define SFH_NMERGE 3
 * #define SFH_NBIN 19
 * At worst there are SFH_NMERGE-1 bins of each size.
 * So set SFH_NBIN to the smallest integer greater than
 * (SFH_NMERGE-1)log_2(MAXSNAPS*STEPS/(SFH_NMERGE-1)+1)
 * MAXSNAPS=61, STEPS=20, NSH_NMERGE=3 --> SFH_NBIN=19 (actually uses 19)
 * MAXSNAPS=61, STEPS=20, NSH_NMERGE=5 --> SFH_NBIN=34 (actually uses 33)
 *
 * Usage:
 *
 * On init.c:
 *   create_sfh_bins();
 * Generates the reference structure for storing the star formation histories in
 * logarithmic bins (for each snapshot/time_ step_number_ combination). In the code galaxy
 * structures are adjusted with respect to this structure at each step_number_.
 * double SFH_t[MAXSNAPS][STEPS][SFH_NBIN]; //Time to present (sfh_bin_number_.e. z=0 ?) at the low-z edge of the bin (code units)
 * double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN]; //Time width of the bin (code units)
 * int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN]; //Number of bins merged in each bin (only useful for the merging algorithm)
 * int SFH_ibin[MAXSNAPS][STEPS]; //Last active bin
 *
 * On initialising galaxy p:
 *           sfh_initialise(p);
 * Whenever new stars are created:
 *    sfh_update_bins(p, time_);
 *    Gal[p].sfh_DiskMass[Gal[p].sfh_ibin_]+=added_mass;
 * (could make the above a function call but it hardly seems worth
 * it as would need a separate call for each property to be updated).
 * To merge galaxy from_galaxy_number_ into galaxy p:
 *    sfh_merge(p,from_galaxy_number_);
 * For debugging purposes:
 *    sfh_print(p);
 * will print the time_ bin structure for galaxy p.
 *
 * The following variables are all defined as part of the galaxy structures:
 * int sfh_ibin_
 *    Index of highest bin are currently in use
 * double sfh_age
 *    Time in code units of last call to sph_update_bins.
 *    (Not strictly required, but useful in L-Galaxies to prevent having to
 *    pass time_ explicitly to the save subroutine.)
 * int sfh_dt[SFH_NBIN]
 *    Width of each time_ bin in code units
 * int sfh_t_[SFH_NBIN]
 *    Time at low-z edge of bin in code units
 * float sfh_t_ime[SFH_NBIN]
 *    time_ till output from the middle of the bin in years.
 *    Used only in save.c
 * float sfh_DiskMass[SFH_NBIN]
 *    Disk mass in bin in standard mass units.
 * float sfh_MetalsDiskMass[SFH_NBIN];
 *    Metals locked up in stars in the disk ditto.
 **/

#include "allvars.h"
#include "proto.h"

void sfh_initialise(const int galaxy_number_)
{
  /* Initialises the sfh-variables for a galaxy */
  int sfh_bin_number_;

  for (sfh_bin_number_=0;sfh_bin_number_<SFH_NBIN;sfh_bin_number_++){
    Gal[galaxy_number_].sfh_dt[sfh_bin_number_]=0.;
    Gal[galaxy_number_].sfh_t[sfh_bin_number_]=0.;
    Gal[galaxy_number_].sfh_flag[sfh_bin_number_]=0;
    Gal[galaxy_number_].sfh_Nbins[sfh_bin_number_]=0;
    Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_]=0.;
    Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_]=0.;
    Gal[galaxy_number_].sfh_ICM[sfh_bin_number_]=0.;
    Gal[galaxy_number_].sfh_MetalsDiskMass[sfh_bin_number_]=metals_init();
    Gal[galaxy_number_].sfh_MetalsBulgeMass[sfh_bin_number_]=metals_init();
    Gal[galaxy_number_].sfh_MetalsICM[sfh_bin_number_]=metals_init();
#ifdef INDIVIDUAL_ELEMENTS
    Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_]=elements_init();
    Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_]=elements_init();
    Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_]=elements_init();
#endif
#ifdef TRACK_BURST
    Gal[galaxy_number_].sfh_BurstMass[sfh_bin_number_]=0.;
#endif
  }

  /* Create first bin */
  Gal[galaxy_number_].sfh_ibin = 0;

  /* Age is used for comparing galaxies during mergers, 
   * so needs to have a value set in case a merger happens before stars 
   * form (which can happen). */
  Gal[galaxy_number_].sfh_age = 0.;
}

void sfh_merge(const int to_galaxy_number_, const int from_galaxy_number_)
{
  /* Merge galaxy from_galaxy_number_ into galaxy to_galaxy_number_ */
  int sfh_bin_number_;

  /* Perform minimal test that the two galaxies have the same time_ structure */
  if (Gal[from_galaxy_number_].sfh_ibin != Gal[to_galaxy_number_].sfh_ibin)
  {
    printf("sfh_merge: trying to merge galaxies with different sfh bins\n");
    sfh_print(to_galaxy_number_);
    sfh_print(from_galaxy_number_);
    exit(1);
  }

  /* The zero-ing of galaxy from_galaxy_number_ here is not strictly necessary as galaxy from_galaxy_number_ should
   * cease to exist after merging, but helps to make mass conservation explicit. */
  for(sfh_bin_number_ = 0 ; sfh_bin_number_ <= Gal[to_galaxy_number_].sfh_ibin; sfh_bin_number_++) 
  {
    Gal[to_galaxy_number_].sfh_DiskMass[sfh_bin_number_]+=Gal[from_galaxy_number_].sfh_DiskMass[sfh_bin_number_];
    Gal[to_galaxy_number_].sfh_BulgeMass[sfh_bin_number_]+=Gal[from_galaxy_number_].sfh_BulgeMass[sfh_bin_number_];
    Gal[to_galaxy_number_].sfh_ICM[sfh_bin_number_]+=Gal[from_galaxy_number_].sfh_ICM[sfh_bin_number_];
    Gal[from_galaxy_number_].sfh_DiskMass[sfh_bin_number_]=0.;
    Gal[from_galaxy_number_].sfh_BulgeMass[sfh_bin_number_]=0.;
    Gal[from_galaxy_number_].sfh_ICM[sfh_bin_number_]=0.;
    metals_add_to(&Gal[to_galaxy_number_].sfh_MetalsDiskMass [sfh_bin_number_], Gal[from_galaxy_number_].sfh_MetalsDiskMass [sfh_bin_number_]);
    metals_add_to(&Gal[to_galaxy_number_].sfh_MetalsBulgeMass[sfh_bin_number_], Gal[from_galaxy_number_].sfh_MetalsBulgeMass[sfh_bin_number_]);
    metals_add_to(&Gal[to_galaxy_number_].sfh_MetalsICM      [sfh_bin_number_], Gal[from_galaxy_number_].sfh_MetalsICM      [sfh_bin_number_]);
    Gal[from_galaxy_number_].sfh_MetalsDiskMass[sfh_bin_number_]  = metals_init();
    Gal[from_galaxy_number_].sfh_MetalsBulgeMass[sfh_bin_number_] = metals_init();
    Gal[from_galaxy_number_].sfh_MetalsICM[sfh_bin_number_]       = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
    elements_add_to(&Gal[to_galaxy_number_].sfh_ElementsDiskMass [sfh_bin_number_],Gal[from_galaxy_number_].sfh_ElementsDiskMass [sfh_bin_number_]);
    elements_add_to(&Gal[to_galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_],Gal[from_galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_]);
    elements_add_to(&Gal[to_galaxy_number_].sfh_ElementsICM      [sfh_bin_number_],Gal[from_galaxy_number_].sfh_ElementsICM      [sfh_bin_number_]);
    Gal[from_galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_]=elements_init();
    Gal[from_galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_]=elements_init();
    Gal[from_galaxy_number_].sfh_ElementsICM[sfh_bin_number_]=elements_init();
#endif
#ifdef TRACK_BURST
    Gal[to_galaxy_number_].sfh_BurstMass[sfh_bin_number_]+=Gal[from_galaxy_number_].sfh_BurstMass[sfh_bin_number_];
    Gal[from_galaxy_number_].sfh_BurstMass[sfh_bin_number_]=0.;
#endif
  }
  /* Again, not strictly necessary, but safe. */
  Gal[from_galaxy_number_].sfh_ibin = 0;
  Gal[from_galaxy_number_].sfh_age = 0.;
}


void sfh_print(const int galaxy_number_)
{
  /* Prints out populated sfh_structure.
   * Does sum of Disk + Bulge only. */
  int sfh_bin_number_;

  printf("For galaxy %d:\n",galaxy_number_);
  printf("sfh_ibin=%d\n",Gal[galaxy_number_].sfh_ibin);
  printf("sfh_age=%f\n",Gal[galaxy_number_].sfh_age);
  printf("  sfh_bin    dt   t      Stars      Metals\n");
  for(sfh_bin_number_ = 0; sfh_bin_number_ < SFH_NBIN; sfh_bin_number_++)
    if (Gal[galaxy_number_].sfh_dt[sfh_bin_number_]!=0)
    {
      printf("%5d %5e %5e %12f\n",sfh_bin_number_,Gal[galaxy_number_].sfh_dt[sfh_bin_number_],Gal[galaxy_number_].sfh_t[sfh_bin_number_],(Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_]+Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_]));
      metals_print("..",metals_add_fraction(Gal[galaxy_number_].sfh_MetalsDiskMass[sfh_bin_number_],Gal[galaxy_number_].sfh_MetalsBulgeMass[sfh_bin_number_],1.));
#ifdef INDIVIDUAL_ELEMENTS
      elements_print("..",elements_add_fraction(Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_],Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_],1.));
#endif
      printf(".......................\n");
    }
}


void create_sfh_bins()
{
  double previous_time_, new_time_, deltaT_, time_;
  int snapshot_number_, step_number_, sfh_ibin_, sfh_bin_number_, sfh_Nbins_[SFH_NBIN];
  int ibin_max_=0;
  double sfh_t_[SFH_NBIN];

  for(snapshot_number_ = 0; snapshot_number_ < MAXSNAPS; snapshot_number_++) 
  {
    for(step_number_ = 0; step_number_ < STEPS; step_number_++) 
    {
      for(sfh_bin_number_=0; sfh_bin_number_ < SFH_NBIN; sfh_bin_number_++)
      {
        SFH_t    [snapshot_number_][step_number_][sfh_bin_number_] = 0;
        SFH_dt   [snapshot_number_][step_number_][sfh_bin_number_] = 0;
        SFH_Nbins[snapshot_number_][step_number_][sfh_bin_number_] = 0;
      }
      SFH_ibin   [snapshot_number_][step_number_]                  = 0;
    }
  }

  for(sfh_bin_number_ = 0; sfh_bin_number_ < SFH_NBIN; sfh_bin_number_++)
  {
    sfh_Nbins_[sfh_bin_number_] = 0;
    sfh_t_    [sfh_bin_number_] = 0.;
  }
  sfh_ibin_=0;

        //for(snapshot_number_=0;snapshot_number_<(LastDarkMatterSnapShot+1)-1;snapshot_number_++) {
  for(snapshot_number_=0;snapshot_number_<(LastDarkMatterSnapShot+1);snapshot_number_++) 
  {
    previous_time_ = NumToTime(snapshot_number_    );
    new_time_      = NumToTime(snapshot_number_ + 1);
    deltaT_        = previous_time_ - new_time_;

    for(step_number_ = 0; step_number_ < STEPS; ++step_number_)
    {
      int ibin_;
      int flag_merged_bins_; // Boolean used to check whether have merged bins
      int dt_merge_; // Size of bins that we are checking for merging
      int n_merge_; // Number of bins of this size

      time_ = previous_time_ - (step_number_ + 1.0) * (deltaT_ / STEPS);
      ibin_ = sfh_ibin_;

      //printf("sna=%d step_number_=%d step_number_ time_=%f time_ low=%f\n",
      //                snapshot_number_,step_number_,(previous_time_ - (step_number_ + 0.5) * (deltaT_ / STEPS))*UnitTime_in_years * inv_Hubble_h/1.e9,
      //                (time_)*UnitTime_in_years * inv_Hubble_h/1.e9);
      //Add one extra bin
      if(snapshot_number_==0 && step_number_==0) 
      {
        sfh_t_    [0] = time_;
        sfh_Nbins_[0] = 1;
      }
      else
      {
        ++ibin_;
        if(ibin_ == SFH_NBIN)
        { terminate("sfh_update_bins: too many bins required\n"); }
        ibin_max_ = max(ibin_max_, ibin_);
        sfh_Nbins_[ibin_] = 1;
        sfh_t_    [ibin_] = time_;
      }

      /* Now merge bins where we have SFH_NMERGE bins of the same size.
       * Need to do this iteratively. */
      flag_merged_bins_ = 1;
      while(flag_merged_bins_)
      {
        flag_merged_bins_ = 0;
        dt_merge_         = sfh_Nbins_[0];
        sfh_bin_number_   = 0;
        // Will have checked all bins once dt_merge_ drops to zero
        while(!flag_merged_bins_ && dt_merge_> 0 )
        {
          // Count number of bins of this size
          n_merge_ = 0;
          // The sfh_bin_number_=sfh_bin_number_ below is to suppress a warning message
          for(sfh_bin_number_= sfh_bin_number_; sfh_Nbins_[sfh_bin_number_]==dt_merge_; ++sfh_bin_number_) {++n_merge_; };
          /* If fewer than SFH_NMERGE bins then do nothing
           * (SFH_NMERGE+1 bins if dt_merge_=1)
           * else exit loop and flag for merging */
          if (n_merge_<SFH_NMERGE || (n_merge_==SFH_NMERGE && dt_merge_==1))
          {
            /* In new version of the code, treat smallest bins just like any others */
            //if (n_merge_<SFH_NMERGE) {
            dt_merge_/=2;
            n_merge_=0;
          }
          else 
          {
            flag_merged_bins_=1;
            sfh_bin_number_=sfh_bin_number_ - n_merge_;
          }
        }
      
        /* At this point, if flag_merged_bins_ is set then
         * we have to merge SFH_NMERGE bins into SFH_NMERGE-1. */
        if(flag_merged_bins_)
        {
          /* Merge bins sfh_bin_number_ and sfh_bin_number_+1 */
          sfh_Nbins_[sfh_bin_number_] *= 2;
          sfh_t_[sfh_bin_number_]=sfh_t_[sfh_bin_number_ + 1];
          /* Relabel all the other bins */
          for(sfh_bin_number_ = sfh_bin_number_ + 1; sfh_bin_number_ < ibin_; ++sfh_bin_number_)
          {
            sfh_Nbins_[sfh_bin_number_]=sfh_Nbins_[sfh_bin_number_+1];
            sfh_t_[sfh_bin_number_]=sfh_t_[sfh_bin_number_+1];
          }
          sfh_Nbins_[sfh_bin_number_]=0;
          sfh_t_[sfh_bin_number_]=0.;
          ibin_=sfh_bin_number_-1;
        }
      } // End loop over bin merging

      sfh_ibin_=ibin_;

      //if(step_number_==19)
      //printf("snapshot_number_=%d\n",snapshot_number_+1);
      for(sfh_bin_number_=0;sfh_bin_number_<=sfh_ibin_;sfh_bin_number_++)
      {
        SFH_t    [snapshot_number_][step_number_][sfh_bin_number_] = sfh_t_    [sfh_bin_number_]; //Time to present at the low-z edge of the bin (code units)
        SFH_Nbins[snapshot_number_][step_number_][sfh_bin_number_] = sfh_Nbins_[sfh_bin_number_];//Number of bins merged in each bin (only useful for the merging algorithm)
        SFH_dt   [snapshot_number_][step_number_][sfh_bin_number_] = (sfh_bin_number_==0) ? NumToTime(0) - sfh_t_[sfh_bin_number_] : sfh_t_[sfh_bin_number_ - 1] - sfh_t_[sfh_bin_number_];//Time width of the bin (code units)
      }        
      SFH_ibin   [snapshot_number_][step_number_]                  = sfh_ibin_; //Last active bin
    }//end loop on steps
  }//end loop on snaps

  if(ThisTask==0)
    printf("Max number of SFH bins used = %d\n",ibin_max_+1);
}


void sfh_update_bins(const int galaxy_number_, const int snapshot_number_, const int step_number_, const double time_)
{
  /* Adds new bins as required.
   * Then merges bins whenever you have three or more of the same size.
   * Assumes that time_ counts from zero at the big bang. */
  int sfh_bin_number_, higher_sfh_bin_number_; // loop index

  Gal[galaxy_number_].sfh_age = time_;

  //t=time_/SFH_TIME_INTERVAL;
  //ibin=Gal[galaxy_number_].sfh_ibin;
  int sfh_ibin_=0; // find highest currently active bin in SFH_Nbins (sfh_bin_number_.e. bin in question in for loop below)
  while((sfh_ibin_ < SFH_NBIN - 1) && (SFH_Nbins[snapshot_number_][step_number_][sfh_ibin_ + 1] > 0)) { ++sfh_ibin_; }

  if (Gal[galaxy_number_].sfh_ibin == 0) //sfh_bin_number_.e. If highest active bin is bin 0...
  {
    for(sfh_bin_number_ = 0; sfh_bin_number_ <= sfh_ibin_; ++sfh_bin_number_) 
    {
      Gal[galaxy_number_].sfh_t    [sfh_bin_number_] = SFH_t    [snapshot_number_][step_number_][sfh_bin_number_];
      Gal[galaxy_number_].sfh_Nbins[sfh_bin_number_] = SFH_Nbins[snapshot_number_][step_number_][sfh_bin_number_];
    }
    Gal[galaxy_number_].sfh_ibin = sfh_ibin_;
  }
  else //sfh_bin_number_.e. If highest active bin is > bin 0...
  {
    sfh_bin_number_ = 0;
    while(sfh_bin_number_ <= sfh_ibin_ && sfh_bin_number_ <= Gal[galaxy_number_].sfh_ibin) //Up to 'bin in question'...until highest active bin is reached...
    {
      if(Gal[galaxy_number_].sfh_Nbins[sfh_bin_number_] == SFH_Nbins[snapshot_number_][step_number_][sfh_bin_number_])
      { ++sfh_bin_number_; }        //...and until bin has grown to required size.
      else
      {
        // Merge bins sfh_bin_number_ and sfh_bin_number_+1
        Gal[galaxy_number_].sfh_Nbins                             [sfh_bin_number_] += Gal[galaxy_number_].sfh_Nbins            [sfh_bin_number_ + 1];
        Gal[galaxy_number_].sfh_t                                 [sfh_bin_number_]  = Gal[galaxy_number_].sfh_t                [sfh_bin_number_ + 1];
        Gal[galaxy_number_].sfh_DiskMass                          [sfh_bin_number_] += Gal[galaxy_number_].sfh_DiskMass         [sfh_bin_number_ + 1];
        Gal[galaxy_number_].sfh_BulgeMass                         [sfh_bin_number_] += Gal[galaxy_number_].sfh_BulgeMass        [sfh_bin_number_ + 1];
        Gal[galaxy_number_].sfh_ICM                               [sfh_bin_number_] += Gal[galaxy_number_].sfh_ICM              [sfh_bin_number_ + 1];
        metals_add_to  (&Gal[galaxy_number_].sfh_MetalsDiskMass   [sfh_bin_number_],   Gal[galaxy_number_].sfh_MetalsDiskMass   [sfh_bin_number_ + 1]);
        metals_add_to  (&Gal[galaxy_number_].sfh_MetalsBulgeMass  [sfh_bin_number_],   Gal[galaxy_number_].sfh_MetalsBulgeMass  [sfh_bin_number_ + 1]);
        metals_add_to  (&Gal[galaxy_number_].sfh_MetalsICM        [sfh_bin_number_],   Gal[galaxy_number_].sfh_MetalsICM        [sfh_bin_number_ + 1]);
#ifdef INDIVIDUAL_ELEMENTS                                                             
        elements_add_to(&Gal[galaxy_number_].sfh_ElementsDiskMass [sfh_bin_number_],   Gal[galaxy_number_].sfh_ElementsDiskMass [sfh_bin_number_ + 1]);
        elements_add_to(&Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_],   Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_ + 1]);
        elements_add_to(&Gal[galaxy_number_].sfh_ElementsICM      [sfh_bin_number_],   Gal[galaxy_number_].sfh_ElementsICM      [sfh_bin_number_ + 1]);
#endif /* defined INDIVIDUAL_ELEMENTS */
#ifdef TRACK_BURST
        Gal[galaxy_number_].sfh_BurstMass                         [sfh_bin_number_] += Gal[galaxy_number_].sfh_BurstMass        [sfh_bin_number_ + 1];
#endif /* defined TRACK_BURST */
        // Relabel all the other bins
        for(higher_sfh_bin_number_ = sfh_bin_number_ + 1; higher_sfh_bin_number_<Gal[galaxy_number_].sfh_ibin; ++higher_sfh_bin_number_) 
        {
          Gal[galaxy_number_].sfh_Nbins            [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_Nbins            [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_t                [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_t                [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_DiskMass         [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_DiskMass         [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_BulgeMass        [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_BulgeMass        [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_ICM              [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_ICM              [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_MetalsDiskMass   [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_MetalsDiskMass   [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_MetalsBulgeMass  [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_MetalsBulgeMass  [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_MetalsICM        [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_MetalsICM        [higher_sfh_bin_number_ + 1];
#ifdef INDIVIDUAL_ELEMENTS
          Gal[galaxy_number_].sfh_ElementsDiskMass [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_ElementsDiskMass [higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_ElementsBulgeMass[higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_ElementsBulgeMass[higher_sfh_bin_number_ + 1];
          Gal[galaxy_number_].sfh_ElementsICM      [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_ElementsICM      [higher_sfh_bin_number_ + 1];
#endif /* defined INDIVIDUAL_ELEMENTS */
#ifdef TRACK_BURST
          Gal[galaxy_number_].sfh_BurstMass        [higher_sfh_bin_number_] = Gal[galaxy_number_].sfh_BurstMass[higher_sfh_bin_number_+1];
#endif /* defined TRACK_BURST */
        }

        //set last bin to zero
        Gal[galaxy_number_].sfh_flag             [higher_sfh_bin_number_] = 0;
        Gal[galaxy_number_].sfh_Nbins            [higher_sfh_bin_number_] = 0;
        Gal[galaxy_number_].sfh_t                [higher_sfh_bin_number_] = 0.;
        Gal[galaxy_number_].sfh_DiskMass         [higher_sfh_bin_number_] = 0.;
        Gal[galaxy_number_].sfh_BulgeMass        [higher_sfh_bin_number_] = 0.;
        Gal[galaxy_number_].sfh_ICM              [higher_sfh_bin_number_] = 0.;
        Gal[galaxy_number_].sfh_MetalsDiskMass   [higher_sfh_bin_number_] = metals_init();
        Gal[galaxy_number_].sfh_MetalsBulgeMass  [higher_sfh_bin_number_] = metals_init();
        Gal[galaxy_number_].sfh_MetalsICM        [higher_sfh_bin_number_] = metals_init();
#ifdef INDIVIDUAL_ELEMENTS                                                  
        Gal[galaxy_number_].sfh_ElementsDiskMass [higher_sfh_bin_number_] = elements_init();
        Gal[galaxy_number_].sfh_ElementsBulgeMass[higher_sfh_bin_number_] = elements_init();
        Gal[galaxy_number_].sfh_ElementsICM      [higher_sfh_bin_number_] = elements_init();
#endif /* defined INDIVIDUAL_ELEMENTS */
#ifdef TRACK_BURST
        Gal[galaxy_number_].sfh_BurstMass        [higher_sfh_bin_number_] = 0.;
#endif /* defined TRACK_BURST */
        Gal[galaxy_number_].sfh_ibin                                      = higher_sfh_bin_number_ - 1;
        
        /* If there are no more time_ bins in the galaxy to merge and
         * the last bin still doesn't have the required size
         * re-size it according to the reference structure SFH */
        if(Gal[galaxy_number_].sfh_Nbins[sfh_bin_number_ + 1] == 0) 
        {
          Gal[galaxy_number_].sfh_Nbins[sfh_bin_number_] = SFH_Nbins[snapshot_number_][step_number_][sfh_bin_number_];
          Gal[galaxy_number_].sfh_t    [sfh_bin_number_] = SFH_t    [snapshot_number_][step_number_][sfh_bin_number_];
          ++sfh_bin_number_;
        }
      }
    }

    //no more bins available in the galaxy, fill the rest times from SFH array
    for(; sfh_bin_number_ <= sfh_ibin_; ++sfh_bin_number_)
    {
      Gal[galaxy_number_].sfh_Nbins            [sfh_bin_number_] = SFH_Nbins[snapshot_number_][step_number_][sfh_bin_number_];
      Gal[galaxy_number_].sfh_t                [sfh_bin_number_] = SFH_t    [snapshot_number_][step_number_][sfh_bin_number_];
      Gal[galaxy_number_].sfh_DiskMass         [sfh_bin_number_] = 0.;
      Gal[galaxy_number_].sfh_BulgeMass        [sfh_bin_number_] = 0.;
      Gal[galaxy_number_].sfh_ICM              [sfh_bin_number_] = 0.;
      Gal[galaxy_number_].sfh_MetalsDiskMass   [sfh_bin_number_] = metals_init();
      Gal[galaxy_number_].sfh_MetalsBulgeMass  [sfh_bin_number_] = metals_init();
      Gal[galaxy_number_].sfh_MetalsICM        [sfh_bin_number_] = metals_init();
#ifdef INDIVIDUAL_ELEMENTS                                         
      Gal[galaxy_number_].sfh_ElementsDiskMass [sfh_bin_number_] = elements_init();
      Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_] = elements_init();
      Gal[galaxy_number_].sfh_ElementsICM      [sfh_bin_number_] = elements_init();
#endif /* defined INDIVIDUAL_ELEMENTS */
#ifdef TRACK_BURST
      Gal[galaxy_number_].sfh_BurstMass        [sfh_bin_number_] = 0.;
#endif /* defined TRACK_BURST */
      Gal[galaxy_number_].sfh_ibin                               = sfh_bin_number_;
    }
  }//end else

  Gal[galaxy_number_].sfh_dt[0] = NumToTime(0) - Gal[galaxy_number_].sfh_t[0];
  for(sfh_bin_number_ = 1; sfh_bin_number_ < SFH_NBIN; sfh_bin_number_++)
  { Gal[galaxy_number_].sfh_dt[sfh_bin_number_] = Gal[galaxy_number_].sfh_t[sfh_bin_number_ - 1] - Gal[galaxy_number_].sfh_t[sfh_bin_number_]; }
}
