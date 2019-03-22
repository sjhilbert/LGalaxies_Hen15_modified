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

#include "allvars.h"
#include "proto.h"

#ifdef MCMC
#include "mcmc_vars.h"
#include "mcmc_proto.h"
#endif

#define PARAMETER_TYPE_IS_INT 1
#define PARAMETER_TYPE_IS_FLOAT 2
#define PARAMETER_TYPE_IS_DOUBLE 3
#define PARAMETER_TYPE_IS_STRING 4


#define MAX_N_TAGS 300


/** @file read_parameters.c reads all the parameters in input.par into global variables
 *       that can be used by the code. */

void read_parameter_file(char *file_name_)
{
  FILE *file_;
  char buf_0_[400], buf_1_[400], buf_2_[400], buf_3_[400];
  int i_, j_, n_parameters_ = 0;
  int parameter_type_[MAX_N_TAGS];
  void *parameter_address_[MAX_N_TAGS];
  char parameter_tag_[MAX_N_TAGS][50];
  int warning_flag_ = 0;
  int error_flag_ = 0;

#ifdef PARALLEL
  if(ThisTask == 0)
    printf("\nreading parameter file:\n\n");
#else
  printf("\nreading parameter file:\n\n");
#endif

  strcpy(parameter_tag_[n_parameters_], "OutputDir");
  parameter_address_[n_parameters_] = OutputDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "FileNameGalaxies");
  parameter_address_[n_parameters_] = FileNameGalaxies;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "McFile");
  parameter_address_[n_parameters_] = McFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "SimulationDir");
  parameter_address_[n_parameters_] = SimulationDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "FileWithOutputRedshifts");
  parameter_address_[n_parameters_] = FileWithOutputRedshifts;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

#ifdef SPECIFYFILENR
  strcpy(parameter_tag_[n_parameters_], "FileNrDir");
  parameter_address_[n_parameters_] = FileNrDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;
#endif

  strcpy(parameter_tag_[n_parameters_], "SpecPhotDir");
  parameter_address_[n_parameters_] = SpecPhotDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "PhotPrefix");
  parameter_address_[n_parameters_] = PhotPrefix;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "SpecPhotIMF");
  parameter_address_[n_parameters_] = SpecPhotIMF;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "FileWithFilterNames");
  parameter_address_[n_parameters_] = FileWithFilterNames;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;


  strcpy(parameter_tag_[n_parameters_], "CoolFunctionsDir");
  parameter_address_[n_parameters_] = CoolFunctionsDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MaxMemSize");
  parameter_address_[n_parameters_] = &MaxMemSize;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Hashbits");
  parameter_address_[n_parameters_] = &Hashbits;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

//Variables used in the MCMC
#ifdef MCMC
  strcpy(parameter_tag_[n_parameters_], "MCMCStartingParFile");
  parameter_address_[n_parameters_] = MCMCStartingParFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MCMCParPriorsAndSwitchesFile");
  parameter_address_[n_parameters_] = MCMCParPriorsAndSwitchesFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MCMCObsConstraints");
  parameter_address_[n_parameters_] = MCMCObsConstraints;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MCMCWeightsObsConstraints");
  parameter_address_[n_parameters_] = MCMCWeightsObsConstraints;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "ObsConstraintsDir");
  parameter_address_[n_parameters_] = ObsConstraintsDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MCMCSampleDir");
  parameter_address_[n_parameters_] = MCMCSampleDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

#ifdef MR_PLUS_MRII
  strcpy(parameter_tag_[n_parameters_], "MCMCSampleFilePrefix_MR");
  parameter_address_[n_parameters_] = MCMCSampleFilePrefix_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MCMCSampleFilePrefix_MRII");
  parameter_address_[n_parameters_] = MCMCSampleFilePrefix_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MCMCSampleFile_MR");
  parameter_address_[n_parameters_] = &MCMCSampleFile_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "MCMCSampleFile_MRII");
  parameter_address_[n_parameters_] = &MCMCSampleFile_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;
#else
  strcpy(parameter_tag_[n_parameters_], "MCMCSampleFilePrefix");
  parameter_address_[n_parameters_] = MCMCSampleFilePrefix;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "MCMCSampleFile");
  parameter_address_[n_parameters_] = &MCMCSampleFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;
#endif

#ifdef HALOMODEL
  strcpy(parameter_tag_[n_parameters_], "MCMCHaloModelDir");
  parameter_address_[n_parameters_] = MCMCHaloModelDir;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;
#endif

  strcpy(parameter_tag_[n_parameters_], "MCMCTreeSampleFile");
  parameter_address_[n_parameters_] = &MCMCTreeSampleFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "ChainLength");
  parameter_address_[n_parameters_] = &ChainLength;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "Time_Dependent_PhysPar");
  parameter_address_[n_parameters_] = &Time_Dependent_PhysPar;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "MCMCMode");
  parameter_address_[n_parameters_] = &MCMCMode;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "MCMC_LogStep_Size");
  parameter_address_[n_parameters_] = &MCMC_LogStep_Size;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "MCMC_Initial_Par_Displacement");
  parameter_address_[n_parameters_] = &MCMC_Initial_Par_Displacement;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "MCMC_Minimum_Obs_Error");
  parameter_address_[n_parameters_] = &MCMC_Minimum_Obs_Error;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "AddedErrOnMass");
  parameter_address_[n_parameters_] = &AddedErrOnMass;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "MachineTimeOut");
  parameter_address_[n_parameters_] = &MachineTimeOut;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "JobSubmitCommand");
  parameter_address_[n_parameters_] = JobSubmitCommand;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "JobSubmitFile");
  parameter_address_[n_parameters_] = JobSubmitFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "JobSubmitPipe");
  parameter_address_[n_parameters_] = JobSubmitPipe;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;
#endif
 
  //Variables for the Scaling & Cosmological Parameters

  strcpy(parameter_tag_[n_parameters_], "ScalePos");
  parameter_address_[n_parameters_] = &ScalePos;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "ScaleMass");
  parameter_address_[n_parameters_] = &ScaleMass;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BaryonFrac");
  parameter_address_[n_parameters_] = &BaryonFrac;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Sigma8");
  parameter_address_[n_parameters_] = &Sigma8;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Omega");
  parameter_address_[n_parameters_] = &Omega;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "OmegaLambda");
  parameter_address_[n_parameters_] = &OmegaLambda;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Hubble_h");
  parameter_address_[n_parameters_] = &Hubble_h;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Omega_OriginalCosm");
  parameter_address_[n_parameters_] = &Omega_OriginalCosm;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "OmegaLambda_OriginalCosm");
  parameter_address_[n_parameters_] = &OmegaLambda_OriginalCosm;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Hubble_h_OriginalCosm");
  parameter_address_[n_parameters_] = &Hubble_h_OriginalCosm;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

#ifdef MR_PLUS_MRII  //OPTION for MCMC
  //MR
  strcpy(parameter_tag_[n_parameters_], "FileWithZList_MR");
  parameter_address_[n_parameters_] = FileWithZList_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "PartMass_MR");
  parameter_address_[n_parameters_] = &PartMass_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BoxSize_MR");
  parameter_address_[n_parameters_] = &BoxSize_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "FileWithZList_OriginalCosm_MR");
  parameter_address_[n_parameters_] = FileWithZList_OriginalCosm_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "PartMass_OriginalCosm_MR");
  parameter_address_[n_parameters_] = &PartMass_OriginalCosm_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BoxSize_OriginalCosm_MR");
  parameter_address_[n_parameters_] = &BoxSize_OriginalCosm_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  //MRII
  strcpy(parameter_tag_[n_parameters_], "FileWithZList_MRII");
  parameter_address_[n_parameters_] = FileWithZList_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "PartMass_MRII");
  parameter_address_[n_parameters_] = &PartMass_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BoxSize_MRII");
  parameter_address_[n_parameters_] = &BoxSize_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "FileWithZList_OriginalCosm_MRII");
  parameter_address_[n_parameters_] = FileWithZList_OriginalCosm_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "PartMass_OriginalCosm_MRII");
  parameter_address_[n_parameters_] = &PartMass_OriginalCosm_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BoxSize_OriginalCosm_MRII");
  parameter_address_[n_parameters_] = &BoxSize_OriginalCosm_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;
#else
  strcpy(parameter_tag_[n_parameters_], "FileWithZList");
  parameter_address_[n_parameters_] = FileWithZList;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag_[n_parameters_], "PartMass");
  parameter_address_[n_parameters_] = &PartMass;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BoxSize");
  parameter_address_[n_parameters_] = &BoxSize;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "PartMass_OriginalCosm");
  parameter_address_[n_parameters_] = &PartMass_OriginalCosm;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BoxSize_OriginalCosm");
  parameter_address_[n_parameters_] = &BoxSize_OriginalCosm;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "FileWithZList_OriginalCosm");
  parameter_address_[n_parameters_] = FileWithZList_OriginalCosm;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_STRING;
#endif


#ifdef MR_PLUS_MRII  //OPTION for MCMC
  strcpy(parameter_tag_[n_parameters_], "LastDarkMatterSnapShot_MR");
  parameter_address_[n_parameters_] = &LastDarkMatterSnapShot_MR;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "LastDarkMatterSnapShot_MRII");
  parameter_address_[n_parameters_] = &LastDarkMatterSnapShot_MRII;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;
#else
  strcpy(parameter_tag_[n_parameters_], "LastDarkMatterSnapShot");
  parameter_address_[n_parameters_] = &LastDarkMatterSnapShot;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;
#endif

#ifndef MCMC
  strcpy(parameter_tag_[n_parameters_], "FirstFile");
  parameter_address_[n_parameters_] = &FirstFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "LastFile");
  parameter_address_[n_parameters_] = &LastFile;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;
#endif

  //Physical Recipes
  strcpy(parameter_tag_[n_parameters_], "ReionizationModel");
  parameter_address_[n_parameters_] = &ReionizationModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "DiskRadiusModel");
  parameter_address_[n_parameters_] = &DiskRadiusModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "StarFormationModel");
  parameter_address_[n_parameters_] = &StarFormationModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "FeedbackReheatingModel");
  parameter_address_[n_parameters_] = &FeedbackReheatingModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "FeedbackEjectionModel");
  parameter_address_[n_parameters_] = &FeedbackEjectionModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "FateOfSatellitesGas");
  parameter_address_[n_parameters_] = &FateOfSatellitesGas;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "ReIncorporationModel");
  parameter_address_[n_parameters_] = &ReIncorporationModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "AGNRadioModeModel");
  parameter_address_[n_parameters_] = &AGNRadioModeModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "DiskInstabilityModel");
  parameter_address_[n_parameters_] = &DiskInstabilityModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "BHGrowthInDiskInstabilityModel");
  parameter_address_[n_parameters_] = &BHGrowthInDiskInstabilityModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "HotGasStrippingModel");
  parameter_address_[n_parameters_] = &HotGasStrippingModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "DisruptionModel");
  parameter_address_[n_parameters_] = &DisruptionModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "StarBurstModel");
  parameter_address_[n_parameters_] = &StarBurstModel;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "BulgeFormationInMinorMergersOn");
  parameter_address_[n_parameters_] = &BulgeFormationInMinorMergersOn;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag_[n_parameters_], "MetallicityOption");
  parameter_address_[n_parameters_] = &MetallicityOption;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_INT;

#ifdef LIGHTCONE_OUTPUT
  strcpy(parameter_tag_[n_parameters_], "lightcone_observer_position_0");
  parameter_address_[n_parameters_] = &lightcone_observer_position[0];
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_observer_position_1");
  parameter_address_[n_parameters_] = &lightcone_observer_position[1];
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_observer_position_2");
  parameter_address_[n_parameters_] = &lightcone_observer_position[2];
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_lower_redshift");
  parameter_address_[n_parameters_] = &lightcone_lower_redshift;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_upper_redshift");
  parameter_address_[n_parameters_] = &lightcone_upper_redshift;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_lower_ra");
  parameter_address_[n_parameters_] = &lightcone_lower_ra;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_upper_ra");
  parameter_address_[n_parameters_] = &lightcone_upper_ra;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_lower_dec");
  parameter_address_[n_parameters_] = &lightcone_lower_dec;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag_[n_parameters_], "lightcone_upper_dec");
  parameter_address_[n_parameters_] = &lightcone_upper_dec;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;

  strcpy(parameter_tag_[n_parameters_], "lightcone_lower_stellar_mass");
  parameter_address_[n_parameters_] = &lightcone_lower_stellar_mass;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_FLOAT;

#endif /* defined LIGHTCONE_OUTPUT */

  //Physical Parameters

  strcpy(parameter_tag_[n_parameters_], "Reionization_z0");
  parameter_address_[n_parameters_] = &Reionization_z0;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Reionization_zr");
  parameter_address_[n_parameters_] = &Reionization_zr;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "Yield");
  parameter_address_[n_parameters_] = &Yield;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "RecycleFraction");
  parameter_address_[n_parameters_] = &RecycleFraction;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "ThreshMajorMerger");
  parameter_address_[n_parameters_] = &ThreshMajorMerger;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "MergerTimeMultiplier");
  parameter_address_[n_parameters_] = &MergerTimeMultiplier;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "RamPressureStrip_CutOffMass");
  parameter_address_[n_parameters_] = &RamPressureStrip_CutOffMass;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "SfrEfficiency");
  parameter_address_[n_parameters_] = &SfrEfficiency;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "SfrColdCrit");
  parameter_address_[n_parameters_] = &SfrColdCrit;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "SfrBurstEfficiency");
  parameter_address_[n_parameters_] = &SfrBurstEfficiency;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "SfrBurstSlope");
  parameter_address_[n_parameters_] = &SfrBurstSlope;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "AgnEfficiency");
  parameter_address_[n_parameters_] = &AgnEfficiency;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BlackHoleGrowthRate");
  parameter_address_[n_parameters_] = &BlackHoleGrowthRate;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BlackHoleSeedMass");
  parameter_address_[n_parameters_] = &BlackHoleSeedMass;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "BlackHoleCutoffVelocity");
  parameter_address_[n_parameters_] = &BlackHoleCutoffVelocity;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "FeedbackReheatingEpsilon");
  parameter_address_[n_parameters_] = &FeedbackReheatingEpsilon;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "ReheatPreVelocity");
  parameter_address_[n_parameters_] = &ReheatPreVelocity;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "ReheatSlope");
  parameter_address_[n_parameters_] = &ReheatSlope;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "FeedbackEjectionEfficiency");
  parameter_address_[n_parameters_] = &FeedbackEjectionEfficiency;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "EjectPreVelocity");
  parameter_address_[n_parameters_] = &EjectPreVelocity;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "EjectSlope");
  parameter_address_[n_parameters_] = &EjectSlope;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "ReIncorporationFactor");
  parameter_address_[n_parameters_] = &ReIncorporationFactor;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "EnergySN");
  parameter_address_[n_parameters_] = &EnergySN;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag_[n_parameters_], "EtaSN");
  parameter_address_[n_parameters_] = &EtaSN;
  parameter_type_ [n_parameters_++] = PARAMETER_TYPE_IS_DOUBLE;

  if((file_ = fopen(file_name_, "r")))
  {
    while(!feof(file_))
    {
      *buf_0_ = 0;
      fgets(buf_0_, 200, file_);
      if(sscanf(buf_0_, "%s%s%s", buf_1_, buf_2_, buf_3_) < 2)
        continue;

      if(buf_1_[0] == '%')
        continue;

      for(i_ = 0, j_ = -1; i_ < n_parameters_; i_++)
        if(strcmp(buf_1_, parameter_tag_[i_]) == 0)
        {
          j_ = i_;
          parameter_tag_[i_][0] = 0;
          break;
        }

      if(j_ >= 0)
      {
#ifdef PARALLEL
        if(ThisTask == 0)
          printf("%35s\t%10s\n", buf_1_, buf_2_);
#else
        printf("%35s\t%10s\n", buf_1_, buf_2_);
#endif
        switch (parameter_type_[j_])
        {
          case PARAMETER_TYPE_IS_DOUBLE:
            *((double *) parameter_address_[j_]) = atof(buf_2_);
            break;
          case PARAMETER_TYPE_IS_FLOAT:
            *((float *) parameter_address_[j_]) = atof(buf_2_);
            break;
          case PARAMETER_TYPE_IS_STRING:
            strcpy(parameter_address_[j_], buf_2_);
            break;
          case PARAMETER_TYPE_IS_INT:
            *((int *) parameter_address_[j_]) = atoi(buf_2_);
            break;
          default:
            printf("Error: unrecognized parameter type for parameter %s\n", parameter_tag_[i_]);
            terminate("Error: unrecognized parameter type for parameter\n");
        }
      }
      else
      {
        printf("Warning in file %s:   Tag '%s' not recognized or multiply defined.\n", file_name_, buf_1_);
        warning_flag_ = 1;
      }
    }
    fclose(file_);

    i_ = strlen(OutputDir);
    if(i_ > 0)
      if(OutputDir[i_ - 1] != '/')
        strcat(OutputDir, "/");
  }
  else
  {
    printf("Parameter file %s not found.\n", file_name_);
    error_flag_ = 1;
  }


  for(i_ = 0; i_ < n_parameters_; i_++)
  {
    if(*parameter_tag_[i_])
    {
      printf("Error. I miss a value for parameter_tag_ '%s' in parameter file '%s'.\n", parameter_tag_[i_], file_name_);
      error_flag_ = 1;
    }
  }

  if(warning_flag_)
    printf("Warning: parameter file %s: encountered unrecognized or multiply defined parameter_tags\n", file_name_);

  if(error_flag_)
  {
    printf("Error: parameter file \"%s\" missing or missing values for mandatory parameters\n", file_name_);
    terminate("parameterfile incorrect");
  }
}




/**
 * @brief checks some program parameters for valid values
 *
 * if used properly, may reduce number of sanity checks later in program
 */
void check_program_parameters()
{
  if(!(
   (ReionizationModel == 0) ||
   (ReionizationModel == 1) ||
   (ReionizationModel == 2)
    ))
  {
    printf("invalid value for program parameter encountered:\n  ReionizationModel = %d\n  error: unknown reionization model\n", ReionizationModel);
    terminate("invalid value program parameter encounted (ReionizationModel).");
  }
  
  if(!(
   (DiskRadiusModel == 0) ||
   (DiskRadiusModel == 1) ||
   (DiskRadiusModel == 2)
    ))
  {
    printf("invalid value for program parameter encountered:\n  DiskRadiusModel = %d\n  error: unknown disk radius model\n", DiskRadiusModel);
    terminate("invalid value program parameter encounted (DiskRadiusModel).");
  }

  if(!(
   (StarFormationModel == 0)
    ))
  {
    printf("invalid value for program parameter encountered:\n  StarFormationModel = %d\n  error: unknown star formation model\n", StarFormationModel);
    terminate("invalid value program parameter encounted (StarFormationModel).");
  }

  if(!(
   (FeedbackReheatingModel == 0)
    ))
  {
    printf("invalid value for program parameter encountered:\n  FeedbackReheatingModel = %d\n  error: unknown feedback reheating model\n", FeedbackReheatingModel);
    terminate("invalid value program parameter encounted (FeedbackReheatingModel).");
  }
  
  if(!(
   (FeedbackEjectionModel == 0) ||
   (FeedbackEjectionModel == 1)
    ))
  {
    printf("invalid value for program parameter encountered:\n  FeedbackEjectionModel = %d\n  error: unknown feedback ejection model\n", FeedbackEjectionModel);
    terminate("invalid value program parameter encounted (FeedbackEjectionModel).");
  }
  
  if(!(
   (FateOfSatellitesGas == 0) ||
   (FateOfSatellitesGas == 1)
    ))
  {
    printf("invalid value for program parameter encountered:\n  FateOfSatellitesGas = %d\n  error: unknown model of ejection of satellite galaxies\n", FateOfSatellitesGas);
    terminate("invalid value program parameter encounted (FateOfSatellitesGas).");
  }
   
  if(!(
   (ReIncorporationModel == 0) ||
   (ReIncorporationModel == 1) ||
   (ReIncorporationModel == 2)
    ))
  {
    printf("invalid value for program parameter encountered:\n  ReIncorporationModel = %d\n  error: unknown model for reincorporation time scale\n", ReIncorporationModel);
    terminate("invalid value program parameter encounted (ReIncorporationModel).");
  }
   
  if(!(
   (AGNRadioModeModel == 0) ||
//   (AGNRadioModeModel == 1) ||
   (AGNRadioModeModel == 2) ||
   (AGNRadioModeModel == 3) ||
   (AGNRadioModeModel == 4) ||
   (AGNRadioModeModel == 5)
    ))
  {
    printf("invalid value for program parameter encountered:\n  AGNRadioModeModel = %d\n  error: unknown AGN radio model\n", AGNRadioModeModel);
    terminate("invalid value program parameter encounted (AGNRadioModeModel).");
  }
   
  if(!(
   (DiskInstabilityModel == 0) ||
   (DiskInstabilityModel == 1) 
    ))
  {
    printf("invalid value for program parameter encountered:\n  DiskInstabilityModel = %d\n  error: unknown disk instability model\n", DiskInstabilityModel);
    terminate("invalid value program parameter encounted (DiskInstabilityModel).");
  }
   
  if(!(
   (BHGrowthInDiskInstabilityModel == 0) ||
   (BHGrowthInDiskInstabilityModel == 1) 
    ))
  {
    printf("invalid value for program parameter encountered:\n  BHGrowthInDiskInstabilityModel = %d\n  error: unknown disk instability feeds BH model\n", BHGrowthInDiskInstabilityModel);
    terminate("invalid value program parameter encounted (BHGrowthInDiskInstabilityModel).");
  }
    
  if(!(
   (HotGasStrippingModel == 0) ||
   (HotGasStrippingModel == 1) 
    ))
  {
    printf("invalid value for program parameter encountered:\n  HotGasStrippingModel = %d\n  error: unknown model for stripping gas off satellite galaxies\n", HotGasStrippingModel);
    terminate("invalid value program parameter encounted (HotGasStrippingModel).");
  }
    
  if(!(
   (DisruptionModel == 0) ||
   (DisruptionModel == 1) 
    ))
  {
    printf("invalid value for program parameter encountered:\n  DisruptionModel = %d\n  error: unknown model for tidal disruption of satellite galaxies\n", DisruptionModel);
    terminate("invalid value program parameter encounted (DisruptionModel).");
  }
    
  if(!(
   (StarBurstModel == 0) ||
   (StarBurstModel == 1) 
    ))
  {
    printf("invalid value for program parameter encountered:\n  StarBurstModel = %d\n  error: unknown star burst model\n", StarBurstModel);
    terminate("invalid value program parameter encounted (StarBurstModel).");
  }
    
  if(!(
   (BulgeFormationInMinorMergersOn == 0) ||
   (BulgeFormationInMinorMergersOn == 1) 
    ))
  {
    printf("invalid value for program parameter encountered:\n  BulgeFormationInMinorMergersOn = %d\n  error: unkown model for bulge formation in mergers\n", BulgeFormationInMinorMergersOn);
    terminate("invalid value program parameter encounted (BulgeFormationInMinorMergersOn).");
  }
    
  if(!(
   (MetallicityOption == 0) ||
   (MetallicityOption == 1) 
    ))
  {
    printf("invalid value for program parameter encountered:\n  MetallicityOption = %d\n  error: unkown choice for metallicities in photometric tables from SPS models\n", MetallicityOption);
    terminate("invalid value program parameter encounted (MetallicityOption).");
  }
}


/** @brief compute derived program parameters from primary parameters */
void compute_derived_program_parameters(void)
{
  inv_Hubble_h = 1. / Hubble_h;
  
  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;
  EtaSNcode    = EtaSN * (UNITMASS_IN_G / SOLAR_MASS) * inv_Hubble_h;
}
