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
 * recipe_yields.c
 *
 *  Created on: 18.11.2011
 *      Author: robyates
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


void update_yields_and_return_mass(int galaxy_number_, int central_galaxy_number_, double dt_, int step_number_)
{
	int Zi_;
	double timestep_width_; //Width of current timestep in CODE UNITS
	int time_bin_; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp_, NormSNIIMassEjecRate_actual_, NormSNIaMassEjecRate_actual_, NormAGBMassEjecRate_actual_, NormSNIIMetalEjecRate_actual_, NormSNIaMetalEjecRate_actual_, NormAGBMetalEjecRate_actual_;
#ifdef INDIVIDUAL_ELEMENTS
	double NormSNIIYieldRate_actual_[NUM_ELEMENTS], NormSNIaYieldRate_actual_[NUM_ELEMENTS], NormAGBYieldRate_actual_[NUM_ELEMENTS];
#endif
	double mass_diff_;
	double time_t_, sfh_time_;
	//double time_to_ts; //Time from high-z (upper) edge of SFH bin to middle of current timestep (used for massive SNII to hot) [in Myrs]
	//double tcut; //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]
	double ColdGasSurfaceDensity_, f_wind_, SNIIEjectaToHot_; //Required for metal-rich wind implementation
	double DiskSFR_, step_width_times_DiskSFR_, DiskSFR_physical_units_, step_width_times_DiskSFR_physical_units_, inverse_DiskMass_physical_units_;
	double BulgeSFR_, step_width_times_BulgeSFR_, BulgeSFR_physical_units_, step_width_times_BulgeSFR_physical_units_, inverse_BulgeMass_physical_units_;
	double ICMSFR_, step_width_times_ICMSFR_, ICMSFR_physical_units_, step_width_times_ICMSFR_physical_units_, inverse_ICM_physical_units_;
	double Disk_total_metallicity_, Bulge_total_metallicity_, ICM_total_metallicity_;
	double NormMassEjecRateSumAllTypes_;
	double TotalMassReturnedToColdDiskGas_, TotalMassReturnedToHotGas_;
	int output_number_; //Iterator used for loop over NOUT when updating MassWeightedAge
	double AgeCorrectionDisk_[NOUT];
	double AgeCorrectionBulge_[NOUT];

	TotalMassReturnedToColdDiskGas_=0.0;
	TotalMassReturnedToHotGas_=0.0;

	for(output_number_=0;output_number_<NOUT;output_number_++)
	{
		AgeCorrectionDisk_[output_number_] = 0.0;
		AgeCorrectionBulge_[output_number_] = 0.0;
	}

	timestep_width_ = dt_; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	time_bin_ = (STEPS*Gal[galaxy_number_].SnapNum)+step_number_;//Bin in Yield tables corresponding to current timestep
	time_t_ = NumToTime(Gal[galaxy_number_].SnapNum) - (step_number_ + 0.5) * dt_; //Time from middle of the current timestep to z=0 (used here for MassWeightAge corrections)
	//NB: NumToTime(Gal[galaxy_number_].SnapNum) is the time to z=0 from start of current snapshot
	//    step_number_ is the number of the current timestep (0-19)
	//    dt_ is the width of one timestep within current snapshot
#ifdef METALRICHWIND
	ColdGasSurfaceDensity_ = max(0.0, (Gal[galaxy_number_].ColdGas*(1.0e10 * inv_Hubble_h))/(4.0*3.14159265*Gal[galaxy_number_].GasDiskRadius*Gal[galaxy_number_].GasDiskRadius * inv_Hubble_h));
	f_wind_ = min(1.0, max(0.0, 1.0/(ColdGasSurfaceDensity_/5.0e12))); //Fraction of SN-II ejecta put directly into HotGas
	if (Gal[galaxy_number_].ColdGas != (float)Gal[galaxy_number_].ColdGas) {f_wind_ = 1.0;}
#endif
#ifndef METALRICHWIND
	f_wind_ = 0.0; //For all stellar ejecta (from disk) to ColdGas
#endif

    int sfh_bin_number_;
    for (sfh_bin_number_=0;sfh_bin_number_<=Gal[galaxy_number_].sfh_ibin;sfh_bin_number_++) //LOOP OVER SFH BINS
    {
    	sfh_time_=Gal[galaxy_number_].sfh_t[sfh_bin_number_]+(0.5*Gal[galaxy_number_].sfh_dt[sfh_bin_number_]);
    	//time_to_ts = ((sfh_time_+(0.5*Gal[galaxy_number_].sfh_dt[sfh_bin_number_])) - time_t_)*(UnitTime_in_years * inv_Hubble_h)/1.0e6; //Time from high-z (upper) edge of SFH bin to middle of current timestep [in Myrs]
    	//tcut = 2.0*((Gal[galaxy_number_].Rvir/Gal[galaxy_number_].Vvir)/0.0001); //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]



    //*****************************************
    //ENRICHMENT FROM DISK STARS INTO COLD GAS:
    //*****************************************
    if (Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_] > 0.0)
    {
     	//pre-calculations to speed up the code
    	DiskSFR_ = Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_]/Gal[galaxy_number_].sfh_dt[sfh_bin_number_];
    	step_width_times_DiskSFR_ = timestep_width_ * DiskSFR_;
    	DiskSFR_physical_units_ = DiskSFR_ * (1.0e10 * inv_Hubble_h);
    	step_width_times_DiskSFR_physical_units_ = timestep_width_ * DiskSFR_physical_units_;
    	inverse_DiskMass_physical_units_=Hubble_h/(Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_]*1.0e10);
    	Disk_total_metallicity_=metals_total(Gal[galaxy_number_].sfh_MetalsDiskMass[sfh_bin_number_])/Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_];


    	Zi_ = find_initial_metallicity(galaxy_number_, sfh_bin_number_, 1, 1);
    	//Interpolate the disk metallicity on the lifetimeMetallicities tables:
    	Zi_disp_ = (Disk_total_metallicity_ - lifetimeMetallicities[Zi_])/(lifetimeMetallicities[Zi_+1] - lifetimeMetallicities[Zi_]);
    	if (Zi_disp_ < 0.0) Zi_disp_ = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

    	NormSNIIMassEjecRate_actual_ = NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIaMassEjecRate_actual_ = NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormAGBMassEjecRate_actual_ = NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIIMetalEjecRate_actual_ = NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIaMetalEjecRate_actual_ = NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormAGBMetalEjecRate_actual_ = NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
    	reset_ejection_rates(sfh_bin_number_, Gal[galaxy_number_].sfh_ibin,
    			&NormSNIIMassEjecRate_actual_, &NormSNIIMetalEjecRate_actual_,
    			&NormSNIaMassEjecRate_actual_, &NormAGBMassEjecRate_actual_,
    			&NormSNIaMetalEjecRate_actual_, &NormAGBMetalEjecRate_actual_);
#endif //INSTANTANEOUS_RECYCLE

    	//pre-calculations to speed up the code
     	NormMassEjecRateSumAllTypes_ = NormSNIIMassEjecRate_actual_ + NormSNIaMassEjecRate_actual_ + NormAGBMassEjecRate_actual_;

#ifdef INDIVIDUAL_ELEMENTS
    	int element_number_;
	    for (element_number_=0;element_number_<NUM_ELEMENTS;element_number_++)
	    {
	    	NormSNIIYieldRate_actual_[element_number_] = NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    	NormSNIaYieldRate_actual_[element_number_] = NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    	NormAGBYieldRate_actual_[element_number_] = NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    }
#endif

#ifdef PORTINARI
	    SNIIEjectaToHot_ = max(0.0, f_wind_ * step_width_times_DiskSFR_ * (NormSNIIMetalEjecRate_actual_ + (Disk_total_metallicity_ * NormSNIIMassEjecRate_actual_)));
	    Gal[galaxy_number_].MetalsHotGas.type2 += SNIIEjectaToHot_;
	    Gal[galaxy_number_].MetalsColdGas.type2 += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_ * (NormSNIIMetalEjecRate_actual_ + (Disk_total_metallicity_ * NormSNIIMassEjecRate_actual_)));
#endif
#ifdef CHIEFFI
	    SNIIEjectaToHot_ = max(0.0, f_wind_ * step_width_times_DiskSFR_ * NormSNIIMetalEjecRate_actual_);
	    Gal[galaxy_number_].MetalsHotGas.type2 += SNIIEjectaToHot_;
	    Gal[galaxy_number_].MetalsColdGas.type2 += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_ * NormSNIIMetalEjecRate_actual_);
#endif

#ifndef SNIATOHOT
	    Gal[galaxy_number_].HotGas += SNIIEjectaToHot_;
    	Gal[galaxy_number_].ColdGas += max(0.0, (step_width_times_DiskSFR_ * NormMassEjecRateSumAllTypes_)-SNIIEjectaToHot_);
    	TotalMassReturnedToColdDiskGas_ += max(0.0, (step_width_times_DiskSFR_ * NormMassEjecRateSumAllTypes_)-SNIIEjectaToHot_); //Only use energy from SNe that eject into ColdGas to reheat
	    TotalMassReturnedToHotGas_ += SNIIEjectaToHot_;
	    //TotalMassReturnedToColdDiskGas_ += max(0.0, step_width_times_DiskSFR_ * NormMassEjecRateSumAllTypes_); //Use energy from ALL SNe (that eject into ColdGas and HotGas) to reheat
#else
	    Gal[galaxy_number_].HotGas += max(0.0, step_width_times_DiskSFR_ * NormSNIaMassEjecRate_actual_) + SNIIEjectaToHot_;
	    Gal[galaxy_number_].ColdGas += max(0.0, step_width_times_DiskSFR_ * (NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)-SNIIEjectaToHot_);
	    TotalMassReturnedToColdDiskGas_ += max(0.0, step_width_times_DiskSFR_ * (NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)-SNIIEjectaToHot_); //Only use energy from SNe that eject into ColdGas to reheat
	    //TotalMassReturnedToColdDiskGas_ += max(0.0, step_width_times_DiskSFR_ * (NormMassEjecRateSumAllTypes_)); //Use energy from ALL SNe (that eject into ColdGas and HotGas) to reheat
	    TotalMassReturnedToHotGas_ += max(0.0, step_width_times_DiskSFR_ * NormSNIaMassEjecRate_actual_) + SNIIEjectaToHot_;
#endif

#ifndef SNIATOHOT
	    Gal[galaxy_number_].MetalsColdGas.type1a += max(0.0, step_width_times_DiskSFR_ * NormSNIaMetalEjecRate_actual_);
#else
    	Gal[galaxy_number_].MetalsHotGas.type1a += max(0.0, step_width_times_DiskSFR_ * NormSNIaMetalEjecRate_actual_);
#endif
	    Gal[galaxy_number_].MetalsColdGas.agb += max(0.0, step_width_times_DiskSFR_ * (NormAGBMetalEjecRate_actual_ + (Disk_total_metallicity_ * NormAGBMassEjecRate_actual_)));


#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
#ifndef SNIATOHOT
    		Gal[galaxy_number_].HotGas_elements.H += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[0] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_)); //SN-II ejecta to HotGas in metal-rich wind (f_wind_)
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[0] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_)); //SN-II ejecta to ColdGas (1.0-f_wind_)
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_)); //SN-Ia and AGB ejecta to ColdGas
    		Gal[galaxy_number_].HotGas_elements.He += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[1] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[1] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#ifndef MAINELEMENTS
    		Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.N += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[5] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[5] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[6] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[6] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Si += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[7] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[7] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.S += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[8] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[8] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[9] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[9] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[10] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[10] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#else
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#endif //MAINELEMENTS
#endif //SNIATOHOT
#ifdef SNIATOHOT
    		Gal[galaxy_number_].HotGas_elements.H += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[0] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_)); //SN-II ejecta to HotGas in metal-rich wind (f_wind_)
    		Gal[galaxy_number_].HotGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[0]); //SN-Ia ejecta to HotGas
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[0] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_)); //SN-II ejecta to ColdGas (1.0-f_wind_)
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[0] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_)); //AGB ejecta to ColdGas
    		Gal[galaxy_number_].HotGas_elements.He += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[1] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[1]);
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[1] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[1] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#ifndef MAINELEMENTS
    		Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[2]);
    		Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.N += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[5] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[5]);
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[5] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[5] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[6] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[6]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[6] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[6] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Si += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[7] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[7]);
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[7] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[7] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.S += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[8] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[8]);
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[8] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[8] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[9] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[9]);
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[9] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[9] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[10] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[10]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[10] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[10] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#else //MAINELEMENTS
       		Gal[galaxy_number_].HotGas_elements.O += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
        	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[2]);
        	Gal[galaxy_number_].ColdGas_elements.O += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
        	Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
       		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
        	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[3]);
        	Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
        	Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
       		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
        	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[4]);
        	Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * (NormSNIIYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe * inverse_DiskMass_physical_units_) * NormSNIIMassEjecRate_actual_));
        	Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#endif //MAINELEMENTS
#endif //SNIATOHOT
#endif //PORTINARI
#ifdef CHIEFFI
#ifndef SNIATOHOT
    		Gal[galaxy_number_].HotGas_elements.H += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[0]); //SN-II ejecta to HotGas in metal-rich wind (f_wind_) //NB: No unsynth component required for SN-II ejecta when using the CL04 SN-II yields
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[0]); //SN-II ejecta to ColdGas (1.0-f_wind_)
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_)); //SN-Ia and AGB ejecta to ColdGas
    		Gal[galaxy_number_].HotGas_elements.He += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[1]);
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[1]);
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#ifndef MAINELEMENTS
    		Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
    		Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
    		Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.N += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[5]);
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[5]);
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[6]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[6]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Si += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[7]);
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[7]);
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.S += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[8]);
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[8]);
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[9]);
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[9]);
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[10]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[10]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#else //MAINELEMENTS
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * ((NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#endif //MAINELEMENTS
#endif //SNIATOHOT
#ifdef SNIATOHOT
    		Gal[galaxy_number_].HotGas_elements.H += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[0]); //SN-II ejecta to HotGas in metal-rich wind (f_wind_) //NB: No unsynth component required for SN-II ejecta when using the CL04 SN-II yields
    		Gal[galaxy_number_].HotGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[0]); //SN-Ia ejecta to HotGas
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[0]); //SN-II ejecta to ColdGas (1.0-f_wind_)
    		Gal[galaxy_number_].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[0] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_)); //AGB ejecta to ColdGas
    		Gal[galaxy_number_].HotGas_elements.He += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[1]);
    		Gal[galaxy_number_].HotGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[1]);
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[1]);
    		Gal[galaxy_number_].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[1] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#ifndef MAINELEMENTS
    		Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
        	Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[2]);
        	Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
        	Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.N += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].HotGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[5]);
    		Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[5]);
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[5]);
    		Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[5] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[6]);
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[6]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[6]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[6] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Si += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[7]);
    		Gal[galaxy_number_].HotGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[7]);
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[7]);
    		Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[7] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.S += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[8]);
    		Gal[galaxy_number_].HotGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[8]);
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[8]);
    		Gal[galaxy_number_].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[8] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[9]);
    		Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[9]);
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[9]);
    		Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[9] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[10]);
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[10]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[10]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[10] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#else //MAINELEMENTS
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
    		Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[2]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[2]);
    		Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[2] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[3]);
    		Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[3] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, f_wind_ * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * NormSNIaYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, (1.0-f_wind_) * step_width_times_DiskSFR_physical_units_ * NormSNIIYieldRate_actual_[4]);
    		Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units_ * (NormAGBYieldRate_actual_[4] + (Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormAGBMassEjecRate_actual_));
#endif //MAINELEMENTS
#endif //SNIATOHOT
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS

    	//UPDATE DISK MASS COMPONENTS:
    	/*ROB (13-02-13): All the mass/metals/elements in the stars that die in this timestep are lost from the stellar component.
    	//sfh_bin_number_.e. All the mass/metals/elements in the stars at birth are removed...
    	//...Some goes to the gas (+ newly synthesised component), the rest goes into the 'stellar remnants' which are not tracked and do not contribute to the stellar component's mass/metals/elements budget.*/
    	Gal[galaxy_number_].DiskMass -= max(0.0, step_width_times_DiskSFR_ * NormMassEjecRateSumAllTypes_);
    	Gal[galaxy_number_].MetalsDiskMass.type2 -= max(0.0, step_width_times_DiskSFR_ * (Disk_total_metallicity_ * NormSNIIMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsDiskMass.type1a -= max(0.0, step_width_times_DiskSFR_ * (Disk_total_metallicity_ * NormSNIaMassEjecRate_actual_));
	    Gal[galaxy_number_].MetalsDiskMass.agb -= max(0.0, step_width_times_DiskSFR_ * (Disk_total_metallicity_ * NormAGBMassEjecRate_actual_));

#ifdef INDIVIDUAL_ELEMENTS
	    Gal[galaxy_number_].DiskMass_elements.H -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].H*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
	    Gal[galaxy_number_].DiskMass_elements.He -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].He*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
	    Gal[galaxy_number_].DiskMass_elements.Cb -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Cb*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
	    Gal[galaxy_number_].DiskMass_elements.N -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].N*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
	    Gal[galaxy_number_].DiskMass_elements.O -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].O*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
	    Gal[galaxy_number_].DiskMass_elements.Ne -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ne*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
	    Gal[galaxy_number_].DiskMass_elements.Mg -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Mg*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
	    Gal[galaxy_number_].DiskMass_elements.Si -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Si*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
	    Gal[galaxy_number_].DiskMass_elements.S -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].S*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
	    Gal[galaxy_number_].DiskMass_elements.Ca -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Ca*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
	    Gal[galaxy_number_].DiskMass_elements.Fe -= max(0.0, step_width_times_DiskSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsDiskMass[sfh_bin_number_].Fe*inverse_DiskMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif //INDIVIDUAL_ELEMENTS

	    //Update ages:
    	for(output_number_=0;output_number_<NOUT;output_number_++)
    	{
    		AgeCorrectionDisk_[output_number_] += max(0.0, (sfh_time_-NumToTime(ListOutputSnaps[output_number_]))*(step_width_times_DiskSFR_ * NormMassEjecRateSumAllTypes_));
    		if (AgeCorrectionDisk_[output_number_] < 0.0) AgeCorrectionDisk_[output_number_] = 0.0;
    	}
    } //if (Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_] > 0.0)

    //*****************************************
    //ENRICHMENT FROM BULGE STARS INTO HOT GAS:
    //*****************************************
    if (Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_] > 0.0)
    {
    	//pre-calculations to speed up the code
    	BulgeSFR_ = Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_]/Gal[galaxy_number_].sfh_dt[sfh_bin_number_];
    	step_width_times_BulgeSFR_ = timestep_width_ * BulgeSFR_;
    	BulgeSFR_physical_units_ = BulgeSFR_ * (1.0e10 * inv_Hubble_h);
    	step_width_times_BulgeSFR_physical_units_ = timestep_width_ * BulgeSFR_physical_units_;
    	inverse_BulgeMass_physical_units_=Hubble_h/(Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_]*1.0e10);
    	Bulge_total_metallicity_=metals_total(Gal[galaxy_number_].sfh_MetalsBulgeMass[sfh_bin_number_])/Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_];

    	Zi_ = find_initial_metallicity(galaxy_number_, sfh_bin_number_, 1, 2);
    	//Interpolate the bulge luminosity on the lifetimeMetallicities tables:
    	Zi_disp_ = (Bulge_total_metallicity_ - lifetimeMetallicities[Zi_])/(lifetimeMetallicities[Zi_+1] - lifetimeMetallicities[Zi_]);
    	if (Zi_disp_ < 0.0) Zi_disp_ = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

    	NormSNIIMassEjecRate_actual_ = NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIaMassEjecRate_actual_ = NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormAGBMassEjecRate_actual_ = NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIIMetalEjecRate_actual_ = NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIaMetalEjecRate_actual_ = NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormAGBMetalEjecRate_actual_ = NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);

    	//pre-calculations to speed up the code
    	NormMassEjecRateSumAllTypes_ = NormSNIIMassEjecRate_actual_ + NormSNIaMassEjecRate_actual_ + NormAGBMassEjecRate_actual_;

#ifdef INDIVIDUAL_ELEMENTS
    	int element_number_;
	    for (element_number_=0;element_number_<NUM_ELEMENTS;element_number_++)
	    {
	    	NormSNIIYieldRate_actual_[element_number_] = NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    	NormSNIaYieldRate_actual_[element_number_] = NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    	NormAGBYieldRate_actual_[element_number_] = NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    }
#endif //INDIVIDUAL_ELEMENTS

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
    	reset_ejection_rates(sfh_bin_number_, Gal[galaxy_number_].sfh_ibin,
    			&NormSNIIMassEjecRate_actual_, &NormSNIIMetalEjecRate_actual_,
    			&NormSNIaMassEjecRate_actual_, &NormAGBMassEjecRate_actual_,
    			&NormSNIaMetalEjecRate_actual_, &NormAGBMetalEjecRate_actual_);
#endif //INSTANTANEOUS_RECYCLE

    	//UPDATE HOT GAS COMPONENTS:
#ifndef BULGE_TO_COLD
    	Gal[galaxy_number_].HotGas += max(0.0, step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_);
    	TotalMassReturnedToHotGas_ += max(0.0, step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_);
#ifdef PORTINARI
    	Gal[galaxy_number_].MetalsHotGas.type2 += max(0.0, step_width_times_BulgeSFR_ * (NormSNIIMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormSNIIMassEjecRate_actual_)));
#endif
#ifdef CHIEFFI
    	Gal[galaxy_number_].MetalsHotGas.type2 += max(0.0, step_width_times_BulgeSFR_ * NormSNIIMetalEjecRate_actual_);
#endif
    	//Gal[galaxy_number_].MetalsHotGas.type1a += step_width_times_BulgeSFR_ * (NormSNIaMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormSNIaMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsHotGas.type1a += max(0.0, step_width_times_BulgeSFR_ * NormSNIaMetalEjecRate_actual_);
    	Gal[galaxy_number_].MetalsHotGas.agb += max(0.0, step_width_times_BulgeSFR_ * (NormAGBMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormAGBMassEjecRate_actual_)));

    	/*if (galaxy_number_==0 && Gal[galaxy_number_].sfh_ICM[sfh_bin_number_] > 0.0) {printf("Bulge:\output_number_");}
    	if (galaxy_number_==0 && Gal[galaxy_number_].sfh_ICM[sfh_bin_number_] > 0.0) {printf("%.11f | %.11f %.11f %.11f\output_number_", NormMassEjecRateSumAllTypes_, NormSNIIMassEjecRate_actual_, NormSNIaMassEjecRate_actual_, NormAGBMassEjecRate_actual_);}
    	if (galaxy_number_==0 && Gal[galaxy_number_].sfh_ICM[sfh_bin_number_] > 0.0) {printf("%.11f | %.11f %.11f %.11f | %.11f\output_number_", max(0.0, step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_), max(0.0, step_width_times_BulgeSFR_ * (NormSNIIMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormSNIIMassEjecRate_actual_))), max(0.0, step_width_times_BulgeSFR_ * NormSNIaMetalEjecRate_actual_), max(0.0, step_width_times_BulgeSFR_ * (NormAGBMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormAGBMassEjecRate_actual_))), max(0.0, step_width_times_BulgeSFR_ * NormSNIaMetalEjecRate_actual_) + max(0.0, step_width_times_BulgeSFR_ * (NormAGBMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormAGBMassEjecRate_actual_))));}*/

#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
    	Gal[galaxy_number_].HotGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[0] + NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].H*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[1] + NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].He*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Cb*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].N*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[5] + NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ne*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[6] + NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[7] + NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Si*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[8] + NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].S*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[9] + NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ca*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[10] + NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#else
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[galaxy_number_].HotGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[0] + NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].H*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[galaxy_number_].HotGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[1] + NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].He*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Cb*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].N*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[5] + NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ne*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[6] + NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[7] + NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Si*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[8] + NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].S*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[9] + NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ca*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[10] + NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
#else
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
#endif //MAINELEMENTS
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS
#else //BULGE_TO_COLD
    	Gal[galaxy_number_].ColdGas += max(0.0, step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_);
    	TotalMassReturnedToColdDiskGas_ += max(0.0, step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_);
    	//TotalMassReturnedToHotGas_ += 0.0;
    	//printf("BulgeToCold = %f\output_number_\output_number_", step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_);
#ifdef PORTINARI
    	Gal[galaxy_number_].MetalsColdGas.type2 += max(0.0, step_width_times_BulgeSFR_ * (NormSNIIMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormSNIIMassEjecRate_actual_)));
#endif
#ifdef CHIEFFI
    	Gal[galaxy_number_].MetalsColdGas.type2 += max(0.0, step_width_times_BulgeSFR_ * NormSNIIMetalEjecRate_actual_);
#endif
    	//Gal[galaxy_number_].MetalsColdGas.type1a += step_width_times_BulgeSFR_ * (NormSNIaMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormSNIaMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsColdGas.type1a += max(0.0, step_width_times_BulgeSFR_ * NormSNIaMetalEjecRate_actual_);
    	Gal[galaxy_number_].MetalsColdGas.agb += max(0.0, step_width_times_BulgeSFR_ * (NormAGBMetalEjecRate_actual_ + (Bulge_total_metallicity_ * NormAGBMassEjecRate_actual_)));

#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
    	Gal[galaxy_number_].ColdGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[0] + NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].H*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[1] + NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].He*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Cb*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].N*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[5] + NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ne*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[6] + NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[7] + NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Si*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[8] + NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].S*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[9] + NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ca*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[10] + NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#else
    	Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[galaxy_number_].ColdGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[0] + NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].H*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[galaxy_number_].ColdGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[1] + NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].He*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].ColdGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Cb*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].N*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[5] + NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ne*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[6] + NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[7] + NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Si*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[8] + NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].S*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[9] + NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ca*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[10] + NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
#else
    	Gal[galaxy_number_].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*(NormAGBMassEjecRate_actual_)));
#endif //MAINELEMENTS
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS
#endif //BULGE_TO_COLD

    	//UPDATE BULGE MASS COMPONENTS:
    	Gal[galaxy_number_].BulgeMass -= max(0.0, step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_);
    	Gal[galaxy_number_].MetalsBulgeMass.type2 -= max(0.0, step_width_times_BulgeSFR_ * (Bulge_total_metallicity_ * NormSNIIMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsBulgeMass.type1a -= max(0.0, step_width_times_BulgeSFR_ * (Bulge_total_metallicity_ * NormSNIaMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsBulgeMass.agb -= max(0.0, step_width_times_BulgeSFR_ * (Bulge_total_metallicity_ * NormAGBMassEjecRate_actual_));

#ifdef INDIVIDUAL_ELEMENTS
    	Gal[galaxy_number_].BulgeMass_elements.H -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].H*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].BulgeMass_elements.He -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].He*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].BulgeMass_elements.Cb -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Cb*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].BulgeMass_elements.N -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].N*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
    	Gal[galaxy_number_].BulgeMass_elements.O -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].O*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].BulgeMass_elements.Ne -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ne*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
    	Gal[galaxy_number_].BulgeMass_elements.Mg -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Mg*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].BulgeMass_elements.Si -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Si*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].BulgeMass_elements.S -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].S*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].BulgeMass_elements.Ca -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Ca*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
    	Gal[galaxy_number_].BulgeMass_elements.Fe -= max(0.0, step_width_times_BulgeSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsBulgeMass[sfh_bin_number_].Fe*inverse_BulgeMass_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif //INDIVIDUAL_ELEMENTS

    	//Update ages:
        for(output_number_=0;output_number_<NOUT;output_number_++)
        {
        	AgeCorrectionBulge_[output_number_] += max(0.0, (sfh_time_-NumToTime(ListOutputSnaps[output_number_]))*(step_width_times_BulgeSFR_ * NormMassEjecRateSumAllTypes_));
        	if (AgeCorrectionBulge_[output_number_] < 0.0) AgeCorrectionBulge_[output_number_] = 0.0;
        }
    } //if (Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_] > 0.0) //BULGE


    //*****************************************
    //ENRICHMENT FROM ICL STARS INTO HOT GAS:
    //*****************************************

    if (Gal[galaxy_number_].sfh_ICM[sfh_bin_number_] > 0.0)
    {
    	//pre-calculations to speed up the code
    	ICMSFR_ = Gal[galaxy_number_].sfh_ICM[sfh_bin_number_]/Gal[galaxy_number_].sfh_dt[sfh_bin_number_];
    	step_width_times_ICMSFR_ = timestep_width_ * ICMSFR_;
    	ICMSFR_physical_units_ = ICMSFR_ * (1.0e10 * inv_Hubble_h);
    	step_width_times_ICMSFR_physical_units_ = timestep_width_ * ICMSFR_physical_units_;
    	inverse_ICM_physical_units_=Hubble_h/(Gal[galaxy_number_].sfh_ICM[sfh_bin_number_]*1.0e10);
    	ICM_total_metallicity_=metals_total(Gal[galaxy_number_].sfh_MetalsICM[sfh_bin_number_])/Gal[galaxy_number_].sfh_ICM[sfh_bin_number_];

    	Zi_ = find_initial_metallicity(galaxy_number_, sfh_bin_number_, 1, 3);
    	//Interpolate the ICM metallicity on the lifetimeMetallicities tables:
    	Zi_disp_ = (ICM_total_metallicity_ - lifetimeMetallicities[Zi_])/(lifetimeMetallicities[Zi_+1] - lifetimeMetallicities[Zi_]);
    	if (Zi_disp_ < 0.0) Zi_disp_ = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

    	NormSNIIMassEjecRate_actual_ = NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIIMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIaMassEjecRate_actual_ = NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIaMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormAGBMassEjecRate_actual_ = NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormAGBMassEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIIMetalEjecRate_actual_ = NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIIMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormSNIaMetalEjecRate_actual_ = NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormSNIaMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);
    	NormAGBMetalEjecRate_actual_ = NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_] + ((NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_+1] - NormAGBMetalEjecRate[time_bin_][sfh_bin_number_][Zi_])*Zi_disp_);

    	//pre-calculations to speed up the code
    	NormMassEjecRateSumAllTypes_ = NormSNIIMassEjecRate_actual_ + NormSNIaMassEjecRate_actual_ + NormAGBMassEjecRate_actual_;

#ifdef INDIVIDUAL_ELEMENTS
    	int element_number_;
	    for (element_number_=0;element_number_<NUM_ELEMENTS;element_number_++)
	    {
	    	NormSNIIYieldRate_actual_[element_number_] = NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormSNIIYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    	NormSNIaYieldRate_actual_[element_number_] = NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormSNIaYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    	NormAGBYieldRate_actual_[element_number_] = NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_] + ((NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_+1][element_number_] - NormAGBYieldRate[time_bin_][sfh_bin_number_][Zi_][element_number_])*Zi_disp_);
	    }
#endif //INDIVIDUAL_ELEMENTS

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
    	reset_ejection_rates(sfh_bin_number_, Gal[galaxy_number_].sfh_ibin,
    			&NormSNIIMassEjecRate_actual_, &NormSNIIMetalEjecRate_actual_,
    			&NormSNIaMassEjecRate_actual_, &NormAGBMassEjecRate_actual_,
    			&NormSNIaMetalEjecRate_actual_, &NormAGBMetalEjecRate_actual_);
#endif //INSTANTANEOUS_RECYCLE

    	//UPDATE HOT GAS COMPONENTS:
    	Gal[galaxy_number_].HotGas += max(0.0, step_width_times_ICMSFR_ * NormMassEjecRateSumAllTypes_);
    	TotalMassReturnedToHotGas_ += max(0.0, step_width_times_ICMSFR_ * NormMassEjecRateSumAllTypes_);
#ifdef PORTINARI
    	Gal[galaxy_number_].MetalsHotGas.type2 += max(0.0, step_width_times_ICMSFR_ * (NormSNIIMetalEjecRate_actual_ + (ICM_total_metallicity_ * NormSNIIMassEjecRate_actual_)));
#endif
#ifdef CHIEFFI
    	Gal[galaxy_number_].MetalsHotGas.type2 += max(0.0, step_width_times_ICMSFR_ * NormSNIIMetalEjecRate_actual_);
#endif
    	//Gal[galaxy_number_].MetalsHotGas.type1a += step_width_times_ICMSFR_ * (NormSNIaMetalEjecRate_actual_ + (ICM_total_metallicity_ * NormSNIaMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsHotGas.type1a += max(0.0, step_width_times_ICMSFR_ * NormSNIaMetalEjecRate_actual_);
    	Gal[galaxy_number_].MetalsHotGas.agb += max(0.0, step_width_times_ICMSFR_ * (NormAGBMetalEjecRate_actual_ + (ICM_total_metallicity_ * NormAGBMassEjecRate_actual_)));

    	/*if (galaxy_number_==0) {printf("ICL:\output_number_");}
    	if (galaxy_number_==0) {printf("%.11f | %.11f %.11f %.11f\output_number_", NormMassEjecRateSumAllTypes_, NormSNIIMassEjecRate_actual_, NormSNIaMassEjecRate_actual_, NormAGBMassEjecRate_actual_);}
    	if (galaxy_number_==0) {printf("%.11f | %.11f %.11f %.11f | %.11f\output_number_\output_number_", max(0.0, step_width_times_ICMSFR_ * NormMassEjecRateSumAllTypes_), max(0.0, step_width_times_ICMSFR_ * (NormSNIIMetalEjecRate_actual_ + (ICM_total_metallicity_ * NormSNIIMassEjecRate_actual_))), max(0.0, step_width_times_ICMSFR_ * NormSNIaMetalEjecRate_actual_), max(0.0, step_width_times_ICMSFR_ * (NormAGBMetalEjecRate_actual_ + (ICM_total_metallicity_ * NormAGBMassEjecRate_actual_))), max(0.0, step_width_times_ICMSFR_ * NormSNIaMetalEjecRate_actual_) + max(0.0, step_width_times_ICMSFR_ * (NormAGBMetalEjecRate_actual_ + (ICM_total_metallicity_ * NormAGBMassEjecRate_actual_))));}*/

#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
    	Gal[galaxy_number_].HotGas_elements.H += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[0] + NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].H*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.He += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[1] + NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].He*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Cb*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.N += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].N*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].O*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[5] + NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Ne*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[6] + NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Mg*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Si += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[7] + NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Si*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.S += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[8] + NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].S*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[9] + NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Ca*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[10] + NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Fe*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#else
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].O*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Mg*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Fe*inverse_ICM_physical_units_)*(NormSNIIMassEjecRate_actual_ + NormAGBMassEjecRate_actual_)));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[galaxy_number_].HotGas_elements.H += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[0] + NormSNIaYieldRate_actual_[0] + NormAGBYieldRate_actual_[0]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].H*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[galaxy_number_].HotGas_elements.He += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[1] + NormSNIaYieldRate_actual_[1] + NormAGBYieldRate_actual_[1]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].He*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].HotGas_elements.Cb += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Cb*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.N += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].N*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].O*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ne += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[5] + NormSNIaYieldRate_actual_[5] + NormAGBYieldRate_actual_[5]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Ne*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[6] + NormSNIaYieldRate_actual_[6] + NormAGBYieldRate_actual_[6]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Mg*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Si += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[7] + NormSNIaYieldRate_actual_[7] + NormAGBYieldRate_actual_[7]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Si*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.S += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[8] + NormSNIaYieldRate_actual_[8] + NormAGBYieldRate_actual_[8]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].S*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Ca += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[9] + NormSNIaYieldRate_actual_[9] + NormAGBYieldRate_actual_[9]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Ca*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[10] + NormSNIaYieldRate_actual_[10] + NormAGBYieldRate_actual_[10]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Fe*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
#else
    	Gal[galaxy_number_].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[2] + NormSNIaYieldRate_actual_[2] + NormAGBYieldRate_actual_[2]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].O*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[3] + NormSNIaYieldRate_actual_[3] + NormAGBYieldRate_actual_[3]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Mg*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
    	Gal[galaxy_number_].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units_ * ((NormSNIIYieldRate_actual_[4] + NormSNIaYieldRate_actual_[4] + NormAGBYieldRate_actual_[4]) + (Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Fe*inverse_ICM_physical_units_)*(NormAGBMassEjecRate_actual_)));
#endif //MAINELEMENTS
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS

    	//UPDATE ICL COMPONENTS:
    	Gal[galaxy_number_].ICM -= max(0.0, step_width_times_ICMSFR_ * NormMassEjecRateSumAllTypes_);
    	Gal[galaxy_number_].MetalsICM.type2 -= max(0.0, step_width_times_ICMSFR_ * (Bulge_total_metallicity_ * NormSNIIMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsICM.type1a -= max(0.0, step_width_times_ICMSFR_ * (Bulge_total_metallicity_ * NormSNIaMassEjecRate_actual_));
    	Gal[galaxy_number_].MetalsICM.agb -= max(0.0, step_width_times_ICMSFR_ * (Bulge_total_metallicity_ * NormAGBMassEjecRate_actual_));

#ifdef INDIVIDUAL_ELEMENTS
    	Gal[galaxy_number_].ICM_elements.H -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].H*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].ICM_elements.He -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].He*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].ICM_elements.Cb -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Cb*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].ICM_elements.N -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].N*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
    	Gal[galaxy_number_].ICM_elements.O -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].O*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].ICM_elements.Ne -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Ne*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
    	Gal[galaxy_number_].ICM_elements.Mg -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Mg*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
#ifndef MAINELEMENTS
    	Gal[galaxy_number_].ICM_elements.Si -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Si*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].ICM_elements.S -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].S*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
    	Gal[galaxy_number_].ICM_elements.Ca -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Ca*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif
    	Gal[galaxy_number_].ICM_elements.Fe -= max(0.0, step_width_times_ICMSFR_physical_units_ * ((Gal[galaxy_number_].sfh_ElementsICM[sfh_bin_number_].Fe*inverse_ICM_physical_units_)*NormMassEjecRateSumAllTypes_));
#endif //INDIVIDUAL_ELEMENTS

    	/*//Update ages:
        for(output_number_=0;output_number_<NOUT;output_number_++)
        {
        	AgeCorrectionICM[output_number_] += max(0.0, (sfh_time_-NumToTime(ListOutputSnaps[output_number_]))*(step_width_times_ICMSFR_ * NormMassEjecRateSumAllTypes_));
        }*/
    } //if (Gal[galaxy_number_].sfh_ICM[sfh_bin_number_] > 0.0) //ICM

    } //for (sfh_bin_number_=0;sfh_bin_number_<=Gal[galaxy_number_].sfh_ibin;sfh_bin_number_++) //MAIN LOOP OVER SFH BINS

    /*//CALL SN-FEEDBACK RECIPE: Sending total mass returned to ColdGas to calculate FB energy:
    SN_feedback(galaxy_number_, central_galaxy_number_, TotalMassReturnedToColdDiskGas_);*/

    //Update Mass-weighted ages:
    for(output_number_=0;output_number_<NOUT;output_number_++)
    {
    	Gal[galaxy_number_].MassWeightAge[output_number_] -= (AgeCorrectionDisk_[output_number_]+AgeCorrectionBulge_[output_number_]);
    }
    
#ifdef H2_AND_RINGS
    double TotalMassReturnedToColdDiskGasr[RNUM], TotalMassReturnedToHotGasr[RNUM];
    double Coldmetallicityr[RNUM], Hotmetallicity[RNUM];
    int ii;
    for(ii=0;ii<RNUM;ii++)
    {
    	TotalMassReturnedToColdDiskGasr[ii]= TotalMassReturnedToColdDiskGas_/((float)RNUM);
    	TotalMassReturnedToHotGasr[ii]=TotalMassReturnedToHotGasr/((float)RNUM);
    	Coldmetallicityr[ii]=metals_total(Gal[galaxy_number_].MetalsColdGas)/Gal[galaxy_number_].ColdGas/((float)RNUM);
    	Hotmetallicity[ii]=metals_total(Gal[galaxy_number_].MetalsHotGas)/Gal[galaxy_number_].HotGas/((float)RNUM);
    }
#endif

#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
    if(TotalMassReturnedToColdDiskGas_>0.)
#ifndef H2_AND_RINGS
    	SN_feedback(galaxy_number_, central_galaxy_number_, TotalMassReturnedToColdDiskGas_, ColdGasComponent);
#else
    SN_feedback(galaxy_number_, central_galaxy_number_, TotalMassReturnedToColdDiskGas_, TotalMassReturnedToColdDiskGasr, ColdGasComponent, Coldmetallicityr);
#endif
    if(TotalMassReturnedToHotGas_>0.)
#ifndef H2_AND_RINGS
    	SN_feedback(galaxy_number_, central_galaxy_number_, TotalMassReturnedToHotGas_, HotGasComponent);
#else
    SN_feedback(galaxy_number_, central_galaxy_number_, TotalMassReturnedToHotGas_, TotalMassReturnedToHotGasr, HotGasComponent, Hotmetallicity);
#endif
#endif
}


int find_initial_metallicity(int galaxy_number_, int sfh_bin_number_, int table_type_, StellarComponentType component_)
{
	if (component_ == DiskComponent) //Disk stars
	{
	int sfh_bin_number_, Zi_bin_;
	double initMetals, Z_disk;

	initMetals = metals_total(Gal[galaxy_number_].sfh_MetalsDiskMass[sfh_bin_number_]); //IN [10^10/h Msun]
	Zi_bin_ = -1;
	sfh_bin_number_ = 0;
	if (initMetals == 0.0 || Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_] == 0.0)
	{
		Z_disk = 0.0;
	}
	else Z_disk = initMetals/Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_];

	switch (table_type_)
	{
		case 1: //Lifetime metallicity table
			while (Zi_bin_ == -1)
			{
				if (lifetimeMetallicities[sfh_bin_number_] < Z_disk)
				{
					sfh_bin_number_++;
					if (sfh_bin_number_ == LIFETIME_Z_NUM) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin_ = sfh_bin_number_;
			}
			break;
		case 2: //SN-II metallicity table
			while (Zi_bin_ == -1)
			{
				if (SNIIMetallicities[sfh_bin_number_] < Z_disk)
				{
					sfh_bin_number_++;
					if (sfh_bin_number_ == 5) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin_ = sfh_bin_number_;
			}
			break;
		//case 3 //SNIa yields are NOT metallicity dependent
		case 4: //AGB metallicity table
			while (Zi_bin_ == -1)
			{
				if (AGBMetallicities[sfh_bin_number_] < Z_disk)
				{
					sfh_bin_number_++;
					if (sfh_bin_number_ == 3) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin_ = sfh_bin_number_;
			}
			break;
	}

	if (Zi_bin_ == 0 ) return Zi_bin_;
	else return Zi_bin_-1;
	}
	else if (component_ == BulgeComponent) //Bulge stars
	{
		int sfh_bin_number_, Zi_bin_;
		double initMetals, Z_bulge;

		initMetals = metals_total(Gal[galaxy_number_].sfh_MetalsBulgeMass[sfh_bin_number_]); //IN [10^10/h Msun]
		Zi_bin_ = -1;
		sfh_bin_number_ = 0;
		if (initMetals == 0.0 || Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_] == 0.0)
		{
			Z_bulge = 0.0;
		}
		else Z_bulge = initMetals/Gal[galaxy_number_].sfh_BulgeMass[sfh_bin_number_];

		switch (table_type_)
		{
			case 1: //Lifetime metallicity table
				while (Zi_bin_ == -1)
				{
					if (lifetimeMetallicities[sfh_bin_number_] < Z_bulge) //Gal[galaxy_number_].sfh_MetalsDiskMass[sfh_bin_number_].type2/Gal[galaxy_number_].sfh_DiskMass[sfh_bin_number_])
					{
						sfh_bin_number_++;
						if (sfh_bin_number_ == LIFETIME_Z_NUM) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin_ = sfh_bin_number_;
				}
				break;
			case 2: //SN-II metallicity table
				while (Zi_bin_ == -1)
				{
					if (SNIIMetallicities[sfh_bin_number_] < Z_bulge)
					{
						sfh_bin_number_++;
						if (sfh_bin_number_ == 5) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin_ = sfh_bin_number_;
				}
				break;
			//case 3 //SNIa yields are NOT metallicity dependent
			case 4: //AGB metallicity table
				while (Zi_bin_ == -1)
				{
					if (AGBMetallicities[sfh_bin_number_] < Z_bulge)
					{
						sfh_bin_number_++;
						if (sfh_bin_number_ == 3) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin_ = sfh_bin_number_;
				}
				break;
		}
		if (Zi_bin_ == 0 ) return Zi_bin_;
		else return Zi_bin_-1;
	}
	else if (component_ == ICMComponent) //ICL stars
		{
			int sfh_bin_number_, Zi_bin_;
			double initMetals, Z_ICM;

			initMetals = metals_total(Gal[galaxy_number_].sfh_MetalsICM[sfh_bin_number_]); //IN [10^10/h Msun]
			Zi_bin_ = -1;
			sfh_bin_number_ = 0;
			if (initMetals == 0.0 || Gal[galaxy_number_].sfh_ICM[sfh_bin_number_] == 0.0)
			{
				Z_ICM = 0.0;
			}
			else Z_ICM = initMetals/Gal[galaxy_number_].sfh_ICM[sfh_bin_number_];

			switch (table_type_)
			{
				case 1: //Lifetime metallicity table
					while (Zi_bin_ == -1)
					{
						if (lifetimeMetallicities[sfh_bin_number_] < Z_ICM)
						{
							sfh_bin_number_++;
							if (sfh_bin_number_ == LIFETIME_Z_NUM) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
						}
						else Zi_bin_ = sfh_bin_number_;
					}
					break;
				case 2: //SN-II metallicity table
					while (Zi_bin_ == -1)
					{
						if (SNIIMetallicities[sfh_bin_number_] < Z_ICM)
						{
							sfh_bin_number_++;
							if (sfh_bin_number_ == 5) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
						}
						else Zi_bin_ = sfh_bin_number_;
					}
					break;
				//case 3 //SNIa yields are NOT metallicity dependent
				case 4: //AGB metallicity table
					while (Zi_bin_ == -1)
					{
						if (AGBMetallicities[sfh_bin_number_] < Z_ICM)
						{
							sfh_bin_number_++;
							if (sfh_bin_number_ == 3) Zi_bin_ = sfh_bin_number_; //If galaxy's Z is higher than max Z from table, then just take max Z from table
						}
						else Zi_bin_ = sfh_bin_number_;
					}
					break;
			}
			if (Zi_bin_ == 0 ) return Zi_bin_;
			else return Zi_bin_-1;
		}
	else { printf("Wrong stellar component type for Z_init calculation: Use either DiskComponent (disk), BulgeComponent (bulge) or ICMComponent (ICL)"); exit(1);}
}


#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
 void reset_ejection_rates(int sfh_bin_number_, int sfh_ibin_,
		 double *NormSNIIMassEjecRate_actual_, double *NormSNIIMetalEjecRate_actual_,
		 double *NormSNIaMassEjecRate_actual_, double *NormAGBMassEjecRate_actual_,
		 double *NormSNIaMetalEjecRate_actual_, double *NormAGBMetalEjecRate_actual_)
 {
    	if(sfh_bin_number_ == sfh_ibin_)
    	{
    		*NormSNIIMassEjecRate_actual_ = 0.43;
    		*NormSNIIMetalEjecRate_actual_ = 0.03;
    	}
    	else
    	{
    		*NormSNIIMassEjecRate_actual_ = 0.0;
    		*NormSNIIMetalEjecRate_actual_ = 0.0;
    	}
    	*NormSNIaMassEjecRate_actual_ = 0.0;
    	*NormAGBMassEjecRate_actual_ =  0.0;
    	*NormSNIaMetalEjecRate_actual_ = 0.0;
    	*NormAGBMetalEjecRate_actual_ =  0.0;
 }
#endif //INSTANTANEOUS_RECYCLE

