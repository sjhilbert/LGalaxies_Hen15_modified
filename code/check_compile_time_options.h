/*  Copyright (C) <2016+>  <L-Galaxies>
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

 /** @file check_compile_time_options.h
 *
 * @brief check  Makefile options.
 *
 * compile time (and runtime) check of
 * compile time options
 * */

#ifndef CHECK_COMPILE_TIME_OPTIONS_H
#define CHECK_COMPILE_TIME_OPTIONS_H
  
#include <stdio.h>
#include <stdlib.h>

#ifdef PARALLEL
#include <mpi.h>
#endif /* defined PARALLEL */

#include "allvars.h"
#include "proto.h"


/**
 * @brief Check whether makefile options are compatible.
 *
 * employs primarily compile time checks by preprocessor.
 */
static inline void 
check_compile_time_options(void)
{
#ifdef OUTPUT_OBS_MAGS
#ifndef COMPUTE_OBS_MAGS
#error "Makefile option OUTPUT_OBS MAGS requires option COMPUTE_OBS_MAGS"
/*  terminate("\n\n> Error : Makefile option OUTPUT_OBS MAGS requires option COMPUTE_OBS_MAGS \n"); */
#endif /* not defined COMPUTE_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */

#ifdef OUTPUT_MOMAF_INPUTS
#ifndef COMPUTE_OBS_MAGS 
#error "Makefile option OUTPUT_MOMAF_INPUTS requires option COMPUTE_OBS_MAGS"
/*   terminate("\n\n> Error : Makefile option OUTPUT_MOMAF_INPUTS requires option COMPUTE_OBS_MAGS \n"); */
#endif /* not defined COMPUTE_OBS_MAGS */
#endif /* defined OUTPUT_MOMAF_INPUTS */

#ifdef KITZBICHLER
#ifndef OUTPUT_MOMAF_INPUTS
#error "Makefile option KITZBICHLER requires option OUTPUT_MOMAF_INPUTS"
/*   terminate("\n\n> Error : Makefile option KITZBICHLER requires option OUTPUT_MOMAF_INPUTS \n"); */
#endif /* not defined OUTPUT_MOMAF_INPUTS */
#endif /* defined KITZBICHLER */

#ifdef GALAXYTREE
#ifndef LOADIDS
#error "Error : Makefile option GALAXYTREE requires LOADIDS"
/* terminate("\n\n> Error : Makefile option GALAXYTREE requires LOADIDS \n"); */
#endif /* not defined LOADIDS */
#endif /* defined GALAXYTREE */

#ifdef POST_PROCESS_MAGS
#ifndef STAR_FORMATION_HISTORY
#error "Makefile option POST_PROCESS_MAGS  requires STAR_FORMATION_HISTORY"
/* terminate("\n\n> Error : Makefile option POST_PROCESS_MAGS  requires STAR_FORMATION_HISTORY \n"); */
#endif /* not defined STAR_FORMATION_HISTORY */
#endif /* defined POST_PROCESS_MAGS */

#ifdef MCMC
#ifndef LOADIDS
#error "Makefile option MCMC requires LOADIDS"
/* terminate("\n\n> Error : Makefile option MCMC requires LOADIDS \n"); */
#endif /* not defined LOADIDS */
#endif /* defined MCMC */

#ifdef HALOMODEL
#ifdef MR_PLUS_MRII
#error "Makefile option HALOMODEL doesn't work yet with MR_PLUS_MRII"
/*   terminate("\n\n> Error : Makefile option HALOMODEL doesn't work yet with MR_PLUS_MRII\n"); */
#endif /* defined MR_PLUS_MRII */
#ifdef MRII
#error "Error : Makefile option HALOMODEL doesn't work yet with MRII"
/*   terminate("\n\n> Error : Makefile option HALOMODEL doesn't work yet with MRII\n"); */
#endif /* defined MRII */
#endif /* defined HALOMODEL */

#ifdef PHOTTABLES_PRECOMPUTED
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#error "Makefile option PHOTTABLES_PRECOMPUTED cannot be run with SPEC_PHOTABLES_ON_THE_FLY"
/* terminate("\n\n> Error : Makefile option PHOTTABLES_PRECOMPUTED cannot run with SPEC_PHOTABLES_ON_THE_FLY\n"); */
#endif /* defined PHOTTABLES_PRECOMPUTED */
#endif /* defined SPEC_PHOTABLES_ON_THE_FLY */

#ifdef LIGHT_OUTPUT
#ifdef POST_PROCESS_MAGS
#error "Makefile option LIGHT_OUTPUT cannot be run with POST_PROCESS_MAGS"
/* terminate("\n\n> Error : Makefile option LIGHT_OUTPUT cannot run with POST_PROCESS_MAGS \n"); */
#endif /* defined POST_PROCESS_MAGS */
#ifdef OUTPUT_MOMAF_INPUTS
#error "Makefile option LIGHT_OUTPUT cannot be run with OUTPUT_MOMAF_INPUTS"
/* terminate("\n\n> Error : Makefile option LIGHT_OUTPUT cannot run with OUTPUT_MOMAF_INPUTS \n"); */
#endif /* defined OUTPUT_MOMAF_INPUTS */
#endif /* defined LIGHT_OUTPUT */

#ifdef LIGHTCONE_OUTPUT
#ifdef MCMC
#error "Makefile option LIGHTCONE_OUTPUT cannot be run with MCMC"
/* terminate("\n\n> Error : Makefile option LIGHTCONE_OUTPUT cannot run with MCMC\n"); */
#endif /* defined MCMC */
#endif /* defined LIGHTCONE_OUTPUT */

#ifdef LIGHTCONE_OUTPUT
#ifdef FIX_OUTPUT_UNITS
#error "Makefile option LIGHTCONE_OUTPUT cannot be run with FIX_OUTPUT_UNITS"
/* terminate("\n\n> Error : Makefile option LIGHTCONE_OUTPUT cannot run with FIX_OUTPUT_UNITS\n"); */
#endif /* defined FIX_OUTPUT_UNITS */
#endif /* defined LIGHTCONE_OUTPUT */

#ifndef LIGHTCONE_OUTPUT
#ifdef LIGHTCONE_OUTPUT_ONLY
#error "Makefile option LIGHTCONE_OUTPUT_ONLY requires LIGHTCONE_OUTPUT"
/* terminate("\n\n> Error : Makefile option LIGHTCONE_OUTPUT_ONLY requires LIGHTCONE_OUTPUT \n"); */
#endif /* defined LIGHTCONE_OUTPUT_ONLY */
#endif /* not defined LIGHTCONE_OUTPUT */
}

#endif /* header guard */
