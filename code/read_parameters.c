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

#define PARAMETER_TYPE_IS_DOUBLE 1
#define PARAMETER_TYPE_IS_FLOAT 2
#define PARAMETER_TYPE_IS_STRING 3
#define PARAMETER_TYPE_IS_INT 4
#define MAX_N_TAGS 300

/** @file read_parameters.c reads all the parameters in input.par into global variables
 *       that can be used by the code. */

void read_parameter_file(char *file_name)
{
  FILE *fd;
  char buf[400], buf1[400], buf2[400], buf3[400];
  int i, j, nt = 0;
  int parameter_type[MAX_N_TAGS];
  void *parameter_address[MAX_N_TAGS];
  char parameter_tag[MAX_N_TAGS][50];
  int warningFlag = 0;
  int errorFlag = 0;

#ifdef PARALLEL
  if(ThisTask == 0)
    printf("\nreading parameter file:\n\n");
#else
  printf("\nreading parameter file:\n\n");
#endif

  strcpy(parameter_tag[nt], "OutputDir");
  parameter_address[nt] = OutputDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "FileNameGalaxies");
  parameter_address[nt] = FileNameGalaxies;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "McFile");
  parameter_address[nt] = McFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "SimulationDir");
  parameter_address[nt] = SimulationDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "FileWithOutputRedshifts");
  parameter_address[nt] = FileWithOutputRedshifts;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

#ifdef SPECIFYFILENR
  strcpy(parameter_tag[nt], "FileNrDir");
  parameter_address[nt] = FileNrDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;
#endif

  strcpy(parameter_tag[nt], "SpecPhotDir");
  parameter_address[nt] = SpecPhotDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "PhotPrefix");
  parameter_address[nt] = PhotPrefix;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "SpecPhotIMF");
  parameter_address[nt] = SpecPhotIMF;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "FileWithFilterNames");
  parameter_address[nt] = FileWithFilterNames;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;


  strcpy(parameter_tag[nt], "CoolFunctionsDir");
  parameter_address[nt] = CoolFunctionsDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MaxMemSize");
  parameter_address[nt] = &MaxMemSize;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Hashbits");
  parameter_address[nt] = &Hashbits;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

//Variables used in the MCMC
#ifdef MCMC
  strcpy(parameter_tag[nt], "MCMCStartingParFile");
  parameter_address[nt] = MCMCStartingParFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MCMCParPriorsAndSwitchesFile");
  parameter_address[nt] = MCMCParPriorsAndSwitchesFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MCMCObsConstraints");
  parameter_address[nt] = MCMCObsConstraints;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MCMCWeightsObsConstraints");
  parameter_address[nt] = MCMCWeightsObsConstraints;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "ObsConstraintsDir");
  parameter_address[nt] = ObsConstraintsDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MCMCSampleDir");
  parameter_address[nt] = MCMCSampleDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

#ifdef MR_PLUS_MRII
  strcpy(parameter_tag[nt], "MCMCSampleFilePrefix_MR");
  parameter_address[nt] = MCMCSampleFilePrefix_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MCMCSampleFilePrefix_MRII");
  parameter_address[nt] = MCMCSampleFilePrefix_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MCMCSampleFile_MR");
  parameter_address[nt] = &MCMCSampleFile_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "MCMCSampleFile_MRII");
  parameter_address[nt] = &MCMCSampleFile_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;
#else
  strcpy(parameter_tag[nt], "MCMCSampleFilePrefix");
  parameter_address[nt] = MCMCSampleFilePrefix;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "MCMCSampleFile");
  parameter_address[nt] = &MCMCSampleFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;
#endif

#ifdef HALOMODEL
  strcpy(parameter_tag[nt], "MCMCHaloModelDir");
  parameter_address[nt] = MCMCHaloModelDir;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;
#endif

  strcpy(parameter_tag[nt], "MCMCTreeSampleFile");
  parameter_address[nt] = &MCMCTreeSampleFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "ChainLength");
  parameter_address[nt] = &ChainLength;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "Time_Dependent_PhysPar");
  parameter_address[nt] = &Time_Dependent_PhysPar;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "MCMCMode");
  parameter_address[nt] = &MCMCMode;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "MCMC_LogStep_Size");
  parameter_address[nt] = &MCMC_LogStep_Size;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "MCMC_Initial_Par_Displacement");
  parameter_address[nt] = &MCMC_Initial_Par_Displacement;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "MCMC_Minimum_Obs_Error");
  parameter_address[nt] = &MCMC_Minimum_Obs_Error;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "AddedErrOnMass");
  parameter_address[nt] = &AddedErrOnMass;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "MachineTimeOut");
  parameter_address[nt] = &MachineTimeOut;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "JobSubmitCommand");
  parameter_address[nt] = JobSubmitCommand;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "JobSubmitFile");
  parameter_address[nt] = JobSubmitFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "JobSubmitPipe");
  parameter_address[nt] = JobSubmitPipe;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;
#endif
 
  //Variables for the Scaling & Cosmological Parameters

  strcpy(parameter_tag[nt], "ScalePos");
  parameter_address[nt] = &ScalePos;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "ScaleMass");
  parameter_address[nt] = &ScaleMass;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BaryonFrac");
  parameter_address[nt] = &BaryonFrac;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Sigma8");
  parameter_address[nt] = &Sigma8;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Omega");
  parameter_address[nt] = &Omega;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "OmegaLambda");
  parameter_address[nt] = &OmegaLambda;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Hubble_h");
  parameter_address[nt] = &Hubble_h;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Omega_OriginalCosm");
  parameter_address[nt] = &Omega_OriginalCosm;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "OmegaLambda_OriginalCosm");
  parameter_address[nt] = &OmegaLambda_OriginalCosm;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Hubble_h_OriginalCosm");
  parameter_address[nt] = &Hubble_h_OriginalCosm;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

#ifdef MR_PLUS_MRII  //OPTION for MCMC
  //MR
  strcpy(parameter_tag[nt], "FileWithZList_MR");
  parameter_address[nt] = FileWithZList_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "PartMass_MR");
  parameter_address[nt] = &PartMass_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BoxSize_MR");
  parameter_address[nt] = &BoxSize_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "FileWithZList_OriginalCosm_MR");
  parameter_address[nt] = FileWithZList_OriginalCosm_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "PartMass_OriginalCosm_MR");
  parameter_address[nt] = &PartMass_OriginalCosm_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BoxSize_OriginalCosm_MR");
  parameter_address[nt] = &BoxSize_OriginalCosm_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  //MRII
  strcpy(parameter_tag[nt], "FileWithZList_MRII");
  parameter_address[nt] = FileWithZList_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "PartMass_MRII");
  parameter_address[nt] = &PartMass_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BoxSize_MRII");
  parameter_address[nt] = &BoxSize_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "FileWithZList_OriginalCosm_MRII");
  parameter_address[nt] = FileWithZList_OriginalCosm_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "PartMass_OriginalCosm_MRII");
  parameter_address[nt] = &PartMass_OriginalCosm_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BoxSize_OriginalCosm_MRII");
  parameter_address[nt] = &BoxSize_OriginalCosm_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;
#else
  strcpy(parameter_tag[nt], "FileWithZList");
  parameter_address[nt] = FileWithZList;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;

  strcpy(parameter_tag[nt], "PartMass");
  parameter_address[nt] = &PartMass;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BoxSize");
  parameter_address[nt] = &BoxSize;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "PartMass_OriginalCosm");
  parameter_address[nt] = &PartMass_OriginalCosm;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BoxSize_OriginalCosm");
  parameter_address[nt] = &BoxSize_OriginalCosm;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "FileWithZList_OriginalCosm");
  parameter_address[nt] = FileWithZList_OriginalCosm;
  parameter_type [nt++] = PARAMETER_TYPE_IS_STRING;
#endif


#ifdef MR_PLUS_MRII  //OPTION for MCMC
  strcpy(parameter_tag[nt], "LastDarkMatterSnapShot_MR");
  parameter_address[nt] = &LastDarkMatterSnapShot_MR;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "LastDarkMatterSnapShot_MRII");
  parameter_address[nt] = &LastDarkMatterSnapShot_MRII;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;
#else
  strcpy(parameter_tag[nt], "LastDarkMatterSnapShot");
  parameter_address[nt] = &LastDarkMatterSnapShot;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;
#endif

#ifndef MCMC
  strcpy(parameter_tag[nt], "FirstFile");
  parameter_address[nt] = &FirstFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "LastFile");
  parameter_address[nt] = &LastFile;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;
#endif

  //Physical Recipes
  strcpy(parameter_tag[nt], "ReionizationModel");
  parameter_address[nt] = &ReionizationModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "DiskRadiusModel");
  parameter_address[nt] = &DiskRadiusModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "StarFormationModel");
  parameter_address[nt] = &StarFormationModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "FeedbackReheatingModel");
  parameter_address[nt] = &FeedbackReheatingModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "FeedbackEjectionModel");
  parameter_address[nt] = &FeedbackEjectionModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "FateOfSatellitesGas");
  parameter_address[nt] = &FateOfSatellitesGas;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "ReIncorporationModel");
  parameter_address[nt] = &ReIncorporationModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "AGNRadioModeModel");
  parameter_address[nt] = &AGNRadioModeModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "DiskInstabilityModel");
  parameter_address[nt] = &DiskInstabilityModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "BHGrowthInDiskInstabilityModel");
  parameter_address[nt] = &BHGrowthInDiskInstabilityModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "HotGasStrippingModel");
  parameter_address[nt] = &HotGasStrippingModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "DisruptionModel");
  parameter_address[nt] = &DisruptionModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "StarBurstModel");
  parameter_address[nt] = &StarBurstModel;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "BulgeFormationInMinorMergersOn");
  parameter_address[nt] = &BulgeFormationInMinorMergersOn;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

  strcpy(parameter_tag[nt], "MetallicityOption");
  parameter_address[nt] = &MetallicityOption;
  parameter_type [nt++] = PARAMETER_TYPE_IS_INT;

#ifdef LIGHTCONE_OUTPUT
  strcpy(parameter_tag[nt], "lightcone_observer_position_0");
  parameter_address[nt] = &lightcone_observer_position[0];
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_observer_position_1");
  parameter_address[nt] = &lightcone_observer_position[1];
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_observer_position_2");
  parameter_address[nt] = &lightcone_observer_position[2];
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_lower_redshift");
  parameter_address[nt] = &lightcone_lower_redshift;
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_upper_redshift");
  parameter_address[nt] = &lightcone_upper_redshift;
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_lower_ra");
  parameter_address[nt] = &lightcone_lower_ra;
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_upper_ra");
  parameter_address[nt] = &lightcone_upper_ra;
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_lower_dec");
  parameter_address[nt] = &lightcone_lower_dec;
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;
  
  strcpy(parameter_tag[nt], "lightcone_upper_dec");
  parameter_address[nt] = &lightcone_upper_dec;
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;

  strcpy(parameter_tag[nt], "lightcone_lower_stellar_mass");
  parameter_address[nt] = &lightcone_lower_stellar_mass;
  parameter_type [nt++] = PARAMETER_TYPE_IS_FLOAT;

#endif /* defined LIGHTCONE_OUTPUT */

  //Physical Parameters

  strcpy(parameter_tag[nt], "Reionization_z0");
  parameter_address[nt] = &Reionization_z0;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Reionization_zr");
  parameter_address[nt] = &Reionization_zr;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "Yield");
  parameter_address[nt] = &Yield;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "RecycleFraction");
  parameter_address[nt] = &RecycleFraction;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "ThreshMajorMerger");
  parameter_address[nt] = &ThreshMajorMerger;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "MergerTimeMultiplier");
  parameter_address[nt] = &MergerTimeMultiplier;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "RamPressureStrip_CutOffMass");
  parameter_address[nt] = &RamPressureStrip_CutOffMass;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "SfrEfficiency");
  parameter_address[nt] = &SfrEfficiency;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "SfrColdCrit");
  parameter_address[nt] = &SfrColdCrit;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "SfrBurstEfficiency");
  parameter_address[nt] = &SfrBurstEfficiency;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "SfrBurstSlope");
  parameter_address[nt] = &SfrBurstSlope;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "AgnEfficiency");
  parameter_address[nt] = &AgnEfficiency;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BlackHoleGrowthRate");
  parameter_address[nt] = &BlackHoleGrowthRate;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BlackHoleSeedMass");
  parameter_address[nt] = &BlackHoleSeedMass;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "BlackHoleCutoffVelocity");
  parameter_address[nt] = &BlackHoleCutoffVelocity;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "FeedbackReheatingEpsilon");
  parameter_address[nt] = &FeedbackReheatingEpsilon;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "ReheatPreVelocity");
  parameter_address[nt] = &ReheatPreVelocity;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "ReheatSlope");
  parameter_address[nt] = &ReheatSlope;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "FeedbackEjectionEfficiency");
  parameter_address[nt] = &FeedbackEjectionEfficiency;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "EjectPreVelocity");
  parameter_address[nt] = &EjectPreVelocity;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "EjectSlope");
  parameter_address[nt] = &EjectSlope;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "ReIncorporationFactor");
  parameter_address[nt] = &ReIncorporationFactor;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "EnergySN");
  parameter_address[nt] = &EnergySN;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  strcpy(parameter_tag[nt], "EtaSN");
  parameter_address[nt] = &EtaSN;
  parameter_type [nt++] = PARAMETER_TYPE_IS_DOUBLE;

  if((fd = fopen(file_name, "r")))
    {
      while(!feof(fd))
        {
          *buf = 0;
          fgets(buf, 200, fd);
          if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
            continue;

          if(buf1[0] == '%')
            continue;

          for(i = 0, j = -1; i < nt; i++)
            if(strcmp(buf1, parameter_tag[i]) == 0)
              {
                j = i;
                parameter_tag[i][0] = 0;
                break;
              }

          if(j >= 0)
            {
#ifdef PARALLEL
              if(ThisTask == 0)
                printf("%35s\t%10s\n", buf1, buf2);
#else
              printf("%35s\t%10s\n", buf1, buf2);
#endif
              switch (parameter_type[j])
                {
                case PARAMETER_TYPE_IS_DOUBLE:
                  *((double *) parameter_address[j]) = atof(buf2);
                  break;
                case PARAMETER_TYPE_IS_FLOAT:
                  *((float *) parameter_address[j]) = atof(buf2);
                  break;
                case PARAMETER_TYPE_IS_STRING:
                  strcpy(parameter_address[j], buf2);
                  break;
                case PARAMETER_TYPE_IS_INT:
                  *((int *) parameter_address[j]) = atoi(buf2);
                  break;
                default:
                  printf("Error: unrecognized parameter type for parameter %s\n", parameter_tag[i]);
                  terminate("Error: unrecognized parameter type for parameter\n");
                }
            }
          else
            {
              printf("Warning in file %s:   Tag '%s' not recognized or multiply defined.\n", file_name, buf1);
              warningFlag = 1;
            }
        }
      fclose(fd);

      i = strlen(OutputDir);
      if(i > 0)
        if(OutputDir[i - 1] != '/')
          strcat(OutputDir, "/");
    }
  else
    {
      printf("Parameter file %s not found.\n", file_name);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*parameter_tag[i])
        {
          printf("Error. I miss a value for parameter_tag '%s' in parameter file '%s'.\n", parameter_tag[i], file_name);
          errorFlag = 1;
        }
    }

  if(warningFlag)
    printf("Warning: parameter file %s: encountered unrecognized or multiply defined parameter_tags\n", file_name);

  if(errorFlag)
  {
    printf("Error: parameter file \"%s\" missing or missing values for mandatory parameters\n", file_name);
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
