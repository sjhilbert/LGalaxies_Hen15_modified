/*  Copyright (C) <2016+>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in_ the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

/** @file save_lightcone.c
 *  @brief   code for outputting galaxies on lightcone to disk
 *
 *  @details copies the relevant properties in_ Galaxy structure into
 *           Galaxy_Output structure and saves them into the output
 *           files (lightcone_*).
 *  
 *  @author  Stefan Hilbert
 *  @author  Rachel Asquith
 * 
 *  @date    2018
 */

#include <stdio.h>
#include <stddef.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

/*** nothing here unless defined LIGHTCONE_OUTPUT ***/
#ifdef LIGHTCONE_OUTPUT
/****************************************************/

#include "lightcone_galaxy_output_type.h"


/** @brief proper integer modulo */
static inline int
modulo(const int i_, const int base_)
{ return (i_ % base_ < 0) ?  i_ % base_ + base_ : i_ % base_; }


/** @brief scalar product for 3d vectors cartesian vectors */
static inline float 
scalar_product(const float v_a_[3], const float v_b_[3])
{ return v_a_[0] * v_b_[0] + v_a_[1] * v_b_[1] + v_a_[2] * v_b_[2]; }


/** @brief returns square of euclidian norm of cartesian vector */
static inline float
euclidian_norm_squared(const float v_[3])
{ return v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]; }


/** @brief returns euclidian norm of cartesian vector */
static inline float
euclidian_norm(const float v_[3])
{ return sqrt(v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]); }


/** @brief returns euclidian distance between two cartesian positions  */
static inline float
euclidian_distance_between(const float pos_a_[3], const float pos_b_[3])
{ return sqrt((pos_b_[0] - pos_a_[0]) * (pos_b_[0] - pos_a_[0]) + (pos_b_[1] - pos_a_[1]) * (pos_b_[1] - pos_a_[1]) + (pos_b_[2] - pos_a_[2]) * (pos_b_[2] - pos_a_[2])); }


/** @brief returns component parallel to a given direction */
static inline float
parallel_component(const float v_[3], const float dir_[3])
{ return scalar_product(v_, dir_) / euclidian_norm(dir_); }


/** @brief shifts cartesian vector by shift_ vecor
 *
 * @param[in]      shift_ 3d shift_ vector
 * @param[in,out]  pos_   3d vector to be shifted
 */
static inline void
apply_shift(const float shift_[3], float (*pos_)[3])
{
  (*pos_)[0] += shift_[0];
  (*pos_)[1] += shift_[1];
  (*pos_)[2] += shift_[2];
}


/** @brief rotates cartesian vector by rotation matrix
 *
 * @param[in]      rot_m_ 3x3 rotation matrix
 * @param[in,out]  v     3d vector to be rotated
 */
static inline void
apply_rotation(/* const */ float rot_m_[3][3], float (*v_)[3])
{
  const float tmp_0_ = (*v_)[0];
  const float tmp_1_ = (*v_)[1];
  const float tmp_2_ = (*v_)[2];
   
  (*v_)[0] = rot_m_[0][0] * tmp_0_ + rot_m_[0][1] * tmp_1_ + rot_m_[0][2] * tmp_2_;
  (*v_)[1] = rot_m_[1][0] * tmp_0_ + rot_m_[1][1] * tmp_1_ + rot_m_[1][2] * tmp_2_;
  (*v_)[2] = rot_m_[2][0] * tmp_0_ + rot_m_[2][1] * tmp_1_ + rot_m_[2][2] * tmp_2_;
}


/** @brief transforms position from cartesian to spherical (ra,dec,r) coords.
 *
 * @param[in,out]  v     3d vector to be transformed
 */
static inline void
apply_cartesian_to_ra_dec_r(float (*pos_)[3])
{
  const float r   = euclidian_norm(*pos_);
  const float ra  = atan2((*pos_)[1], (*pos_)[0]);
  const float dec = asin((*pos_)[2] / r);

  (*pos_)[0] = ra;
  (*pos_)[1] = dec;
  (*pos_)[2] = r;
}


/** @brief computes rotation matrix for rotating cartesian vectors to local orthonormal frame 
 *         defined by unit tangent vectors of spherical (ra, dec, r) coord. system 
 *         using cartesian pos
 * 
 * @param[in]  pos_0_, pos_1_, pos_2_  3d cartesian position around which to set up local orthonormal system
 * @param[out] rot_m_                  rotation matrix to be computed
 *
 * @todo check for correctness
 */
static inline void
get_cartesian_to_ra_dec_r_local_orthogonal_rotation_matrix_from_cartesian_position(const float pos_0_, const float pos_1_, const float pos_2_, float (*rot_m_)[3][3])
{
  const float q_       = sqrt(pos_0_ * pos_0_ + pos_1_ * pos_1_);
  const float r_       = sqrt(pos_0_ * pos_0_ + pos_1_ * pos_1_ + pos_2_ * pos_2_);

  (*rot_m_)[0][0] = - pos_1_ / q_;
  (*rot_m_)[0][1] =   pos_0_ / q_;
  (*rot_m_)[0][2] =   0;
  (*rot_m_)[1][0] = - pos_0_ / r_ * pos_2_ / q_;
  (*rot_m_)[1][1] = - pos_1_ / r_ * pos_2_ / q_;
  (*rot_m_)[1][2] =       q_ / r_;
  (*rot_m_)[2][0] =   pos_0_ / r_;
  (*rot_m_)[2][1] =   pos_1_ / r_;
  (*rot_m_)[2][2] =   pos_2_ / r_;   
}


/** @brief computes rotation matrix for rotating cartesian vectors to local orthonormal frame 
 *         defined by unit tangent vectors of spherical (ra, dec, r) coord. system 
 *         using ra and dec
 * 
 * @param[in]  ra_     right ascension [rad] of position
 * @param[in]  dec_    declination [rad] of position
 * @param[out] rot_m_  rotation matrix to be computed
 */
static inline void
get_cartesian_to_ra_dec_r_local_orthogonal_rotation_matrix_from_ra_dec(const float ra_, const float dec_, float (*rot_m_)[3][3])
{
  const float sin_ra_  = sin(ra_);  /* in old cartesian coords, this should be: Pos[1] / sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]) */
  const float cos_ra_  = cos(ra_);  /* in old cartesian coords, this should be: Pos[0] / sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]) */
  const float sin_dec_ = sin(dec_); /* in old cartesian coords, this should be: Pos[2] / sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1] + Pos[2] * Pos[2]) */ 
  const float cos_dec_ = cos(dec_); /* in old cartesian coords, this should be: sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]) / sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1] + Pos[2] * Pos[2]) */ 

  (*rot_m_)[0][0] = - sin_ra_;
  (*rot_m_)[0][1] =   cos_ra_;
  (*rot_m_)[0][2] = 0;
  (*rot_m_)[1][0] = - sin_dec_ * cos_ra_;
  (*rot_m_)[1][1] = - sin_dec_ * sin_ra_;
  (*rot_m_)[1][2] = cos_dec_;
  (*rot_m_)[2][0] = cos_dec_ * cos_ra_;
  (*rot_m_)[2][1] = cos_dec_ * sin_ra_;
  (*rot_m_)[2][2] = sin_dec_;   
}


/** @brief combined 3d shift index into single integer */
static inline int
composed_cube_shift_index(const int shift_index_[3])
{ 
  const int base_ = 1000;
  return base_ * base_ * modulo(shift_index_[0], base_) + base_ * modulo(shift_index_[1], base_) + modulo(shift_index_[2], base_);  
}


/** @brief sets up global variables needed for lightcone geometry
 */
void 
init_lightcone(void)
{
  lightcone_observer_distance_from_origin = euclidian_norm(lightcone_observer_position);
  
  /* compute slice boundaries for lightcone */
  /* note: the output order is reversed between w/ and w/o GALAXYTREE */
  int output_number_;
  int snapshot_number_;
  
#ifdef GALAXYTREE
  for(output_number_ = 0; output_number_ < NOUT - 1; output_number_++)
  { lightcone_slice_lower_redshift[output_number_] = 0.5 * (ZZ[ListOutputSnaps[output_number_]] + ZZ[ListOutputSnaps[output_number_ + 1]]); };
  lightcone_slice_lower_redshift[NOUT - 1] = min(0., ZZ[ListOutputSnaps[NOUT - 1]]);

  lightcone_slice_upper_redshift[0] = ZZ[ListOutputSnaps[0]];
  for(output_number_ = 1; output_number_ < NOUT; output_number_++)
  { lightcone_slice_upper_redshift[output_number_] = 0.5 * (ZZ[ListOutputSnaps[output_number_ - 1]] + ZZ[ListOutputSnaps[output_number_]]); };

#else /* not defined GALAXYTREE */
  lightcone_slice_lower_redshift[0] = min(0., ZZ[ListOutputSnaps[0]]);
  for(output_number_ = 1; output_number_ < NOUT; output_number_++)
  { lightcone_slice_lower_redshift[output_number_] = 0.5 * (ZZ[ListOutputSnaps[output_number_ - 1]] + ZZ[ListOutputSnaps[output_number_]]); };

  for(output_number_ = 0; output_number_ < NOUT - 1; output_number_++)
  { lightcone_slice_upper_redshift[output_number_] = 0.5 * (ZZ[ListOutputSnaps[output_number_]] + ZZ[ListOutputSnaps[output_number_ + 1]]); };
  lightcone_slice_upper_redshift[NOUT - 1] = ZZ[ListOutputSnaps[NOUT - 1]];
  
#endif /* not defined GALAXYTREE */
  
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    lightcone_slice_lower_los_distance[output_number_] = comoving_los_distance_for_redshift(lightcone_slice_lower_redshift[output_number_]);
    lightcone_slice_upper_los_distance[output_number_] = comoving_los_distance_for_redshift(lightcone_slice_upper_redshift[output_number_]);
  };
  
  printf("\n--- lightcone output slice parameters: --------------------------------------------------\n"
         "slice no.: \t snap no., \t     z_lo, \t   z_snap, \t     z_hi, \t   chi_lo, \t   chi_hi\n"
         "-----------------------------------------------------------------------------------------\n");
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  { printf("slice %d: \t%d, \t%f, \t%f, \t%f, \t%f, \t%f\n", output_number_, ListOutputSnaps[output_number_], lightcone_slice_lower_redshift[output_number_],  ZZ[ListOutputSnaps[output_number_]], lightcone_slice_upper_redshift[output_number_], lightcone_slice_lower_los_distance[output_number_], lightcone_slice_upper_los_distance[output_number_]); }
  printf("------------------------------------------------------------------------------------------\n");

#ifdef GALAXYTREE   
  for(snapshot_number_ = MAXSNAPS; snapshot_number_--; )
  { 
    lightcone_radius_for_snapshot       [snapshot_number_] = lightcone_slice_upper_los_distance[snapshot_number_];
    is_outside_lightcone_for_snapshot   [snapshot_number_] = (lightcone_radius_for_snapshot[snapshot_number_] <= 0.);
    check_outside_lightcone_for_snapshot[snapshot_number_] = 0.866 * BoxSize > lightcone_radius_for_snapshot[snapshot_number_];
  }

#else /* not defined GALAXYTREE */ 
  output_number_ = -1;
  for(snapshot_number_ = MAXSNAPS; snapshot_number_--; )
  {
    while(output_number_ < NOUT - 1 && ListOutputSnaps[output_number_ + 1] >= snapshot_number_)
    { output_number_++; }
    lightcone_radius_for_snapshot       [snapshot_number_] = (output_number_ < 0) ? 0    : lightcone_slice_upper_los_distance[output_number_] + comoving_los_distance_for_redshift(ZZ[snapshot_number_]) - comoving_los_distance_for_redshift(ZZ[ListOutputSnaps[output_number_]]);
    is_outside_lightcone_for_snapshot   [snapshot_number_] = (output_number_ < 0) ? true : (lightcone_radius_for_snapshot[snapshot_number_] <= 0.);
    check_outside_lightcone_for_snapshot[snapshot_number_] = 0.866 * BoxSize > lightcone_radius_for_snapshot[snapshot_number_];
  }
 #endif /* not defined GALAXYTREE */
  
  /*
  printf("\n--- lightcone output snapshot parameters: ------------------------------------------------\n"
         "snap no.: \t     l.c.r., \t   is_outide, \t     check_outside\n"
         "-----------------------------------------------------------------------------------------\n");
  for(snapshot_number_ = MAXSNAPS; snapshot_number_--; )
  { printf("snap %d: \t     %f, \t   %d, \t     %d\n", snapshot_number_, lightcone_radius_for_snapshot[snapshot_number_], is_outside_lightcone_for_snapshot[snapshot_number_], check_outside_lightcone_for_snapshot[snapshot_number_]); }
  printf("------------------------------------\n");
  */
  
  lightcone_N_fof_groups_skipped_construction                       = 0;
  lightcone_N_galaxies_skipped_construction                         = 0;
  lightcone_N_galaxies_skipped_output_early                         = 0;
  lightcone_N_galaxies_for_output                                   = 0;
  lightcone_N_galaxies_remaining_for_output_past_construct_galaxies = 0;
}


/** @brief outputs statistics about lightcone galaxies
 */
void show_lightcone_statistics(void)
{
  long long TotLightconeGalaxies_sum_ = 0;
#ifdef GALAXYTREE
  TotLightconeGalaxies_sum_ += TotLightconeGalCount;
#else /* not defined GALAXYTREE */
  {
    int output_number_;
    for(output_number_ = 0; output_number_ < NOUT; output_number_++)
    {  TotLightconeGalaxies_sum_+= TotLightconeGalaxies[output_number_];   }
  }
#endif /* not defined GALAXYTREE */

  printf("\n--- statistics for lightcone files from last tree file: ---");
  printf("\n  lightcone_N_fof_groups_skipped_construction = %lld", lightcone_N_fof_groups_skipped_construction                                            );
  printf("\n  lightcone_N_galaxies_skipped_construction = %lld", lightcone_N_galaxies_skipped_construction                                                );
  printf("\n  lightcone_N_galaxies_remaining_for_output_past_construct_galaxies = %lld", lightcone_N_galaxies_remaining_for_output_past_construct_galaxies);
  printf("\n  lightcone_N_galaxies_skipped_output_early = %lld", lightcone_N_galaxies_skipped_output_early                                                );
  printf("\n  lightcone_N_galaxies_for_output = %lld", lightcone_N_galaxies_for_output                                                                    );
  printf("\n  lightcone_N_galaxies_for_output_not_skipped_early = %lld", lightcone_N_galaxies_for_output - lightcone_N_galaxies_skipped_output_early      );
  printf("\n  lightcone galaxies written = %lld", TotLightconeGalaxies_sum_                                                                               );
  printf("\n  size of GALAXY_OUTPUT = %lu", sizeof(struct GALAXY_OUTPUT)                                                                                  );
  printf("\n  size of lightcone galaxy on disk = %lu", sizeof(lightcone_galaxy_output_type)                                                               );
#ifdef GALAXYTREE
  const unsigned long long expected_file_size_ = (unsigned long long)(sizeof(lightcone_galaxy_output_type)) * (unsigned long long)(1 + TotLightconeGalCount);
  printf("\n  expected size of lightcone galaxy file on disk = %llu", expected_file_size_                                                                 );
#endif /* defined GALAXYTREE */
  printf("\n");
}


/** @brief creates and opens files for outputting galaxies on lightcone
 */
void 
create_lightcone_galaxy_files(int file_number_)
{
  // create output files 
#ifdef GALAXYTREE
  char file_name_[1000];
  sprintf(file_name_, "%s/lightcone_%s_galtree_%d", OutputDir, FileNameGalaxies, file_number_);
  if(!(FdLightconeGalTree = fopen(file_name_, "wb+")))
    {
      char error_message_[1000];
      sprintf(error_message_, "can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }
  /* fill one block to make room for header */
  char zero = 0;
  myffill(&zero, 1, sizeof(lightcone_galaxy_output_type), FdLightconeGalTree);
  TotLightconeGalCount = 0;
  
#else /* not defined GALAXYTREE */
  char file_name_[1000];
  int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    sprintf(file_name_, "%s/lightcone_%s_z%1.2f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[output_number_]], file_number_);
    if(!(FdLightconeGalDumps[output_number_] = fopen(file_name_, "wb+")))
    {
      char error_message_[1000];
      sprintf(error_message_, "can't open file `%s'\n", file_name_);
      terminate(error_message_);
    }
    
    char zero = 0;
    myffill(&zero, 1, sizeof(long long int), FdLightconeGalDumps[output_number_]);  /* space for the header */
    TotLightconeGalaxies[output_number_] = 0;
  }
#endif /* not defined GALAXYTREE */
}


/** @brief writes header and closes files for outputting galaxies on lightcone
 */
void 
close_lightcone_galaxy_files(void)
{
#ifdef GALAXYTREE
  
  long long one_ = 1;
  long long size_of_struct_ = sizeof(lightcone_galaxy_output_type);

  myfseek(FdLightconeGalTree, 0, SEEK_SET);
  myfwrite(&one_, sizeof(long long), 1, FdLightconeGalTree);                   // write 1 (to determine endianess?)
  myfwrite(&size_of_struct_, sizeof(long long), 1, FdLightconeGalTree);        // size of an output structure (Galaxy_Output)
  myfwrite(&TotLightconeGalCount, sizeof(long long), 1, FdLightconeGalTree);   // the total number of galaxies
#ifdef SORT_GALAXY_OUTPUT
  //  sort_lightcone_galaxy_in_file(FdLightconeGalTree); /* sorting is already done during save_lightcone_galaxy_finalize */
#endif /* defined SORT_GALAXY_OUTPUT */
  fclose(FdLightconeGalTree);

#else /* not defined GALAXYTREE */

  int output_number_;
  for(output_number_ = 0; output_number_ < NOUT; output_number_++)
  {
    fseek(FdLightconeGalDumps[output_number_], 0, SEEK_SET);
    myfwrite(&TotLightconeGalaxies[output_number_], sizeof(long long), 1, FdLightconeGalDumps[output_number_]);  // total number of galaxies
#ifdef SORT_GALAXY_OUTPUT
    sort_lightcone_galaxy_in_file(FdLightconeGalDumps[output_number_]);
#endif /* defined SORT_GALAXY_OUTPUT */
    fclose(FdLightconeGalDumps[output_number_]);
  }
#endif /* not defined GALAXYTREE */
}


/** @brief seeks to position of an entry in lightcone galaxy file
 */
static inline int
myfseek_lightcone_galaxy(FILE* lightcone_galaxy_file_, const long long galaxy_in_file_number_)
{
#ifdef GALAXYTREE
  myfseek(lightcone_galaxy_file_, (long)(1 + galaxy_in_file_number_) * (long)(sizeof(lightcone_galaxy_output_type)), SEEK_SET);
#else  /* not defined GALAXYTREE */
  myfseek(lightcone_galaxy_file_, (long)(sizeof(long long int)) + (long)(galaxy_in_file_number_) * (long)(sizeof(lightcone_galaxy_output_type)), SEEK_SET);
#endif /* not defined GALAXYTREE */
  return 0;
}


/** @brief reads header of lightcone galaxy file to get number of galaxies in file
 */
static inline size_t
myfread_lightcone_galaxy_number_of_entries(FILE* lightcone_galaxy_file_, long long* n_lightcone_galaxies_in_file_)
{
  *n_lightcone_galaxies_in_file_ = 0;
  fseek(lightcone_galaxy_file_, 0, SEEK_SET);
  
#ifdef GALAXYTREE
  long long one_;
  long long size_of_struct_;

  myfread(&one_, sizeof(long long), 1, lightcone_galaxy_file_);                   // 1 (to determine endianess?)
  myfread(&size_of_struct_, sizeof(long long), 1, lightcone_galaxy_file_);        // size of an output structure (Galaxy_Output)
  myfread(n_lightcone_galaxies_in_file_, sizeof(long long), 1, lightcone_galaxy_file_); 
 
  if(one_ != 1)
  { 
    printf("\nerror: in get_number_of_lightcone_galaxies_in_file(FILE* ): error reading file header: one on disk = %llu != 1 (supposed value).\n", one_);
    terminate("error reading file header");
  }
  if(size_of_struct_ != sizeof(lightcone_galaxy_output_type))
  { 
    printf("\nerror: in get_number_of_lightcone_galaxies_in_file(FILE* ): error reading file header: size of struct on disk = %llu != %lu (supposed value).\n", size_of_struct_, sizeof(lightcone_galaxy_output_type));
    terminate("error reading file header");
  }
  
  return 3;
#else /* not defined GALAXYTREE */
  myfread(n_lightcone_galaxies_in_file_, sizeof(long long), 1, lightcone_galaxy_file_); 
  return 1;
#endif /* not defined GALAXYTREE */
}


/** @brief adjusts galaxy properties for lightcone
 * 
 * Pos will be transformed from comoving cartesian to 
 * spherical coordinates relative to lightcone_observer_position
 * and stored as
 * Pos[0] = right ascension [rad]
 * Pos[1] = declination [rad]
 * Pos[2] = comoving radial l.o.s. distance [Mpc/h]
 * 
 * Redshift will be adjusted to cosmological redshift for l.o.s. of galaxy
 * ObsRedshift will be redshift accounting for both cosmological and l.o.s. peculiar 
 * motion (transverse motion neglected)
 *
 * Vel, DistanceToCentralGal, GasSpin, StellarSpin will be rotated into local 
 * orthonormal coordinates defined by spherical coordinate unit tangent vectors 
 * at galaxy position
 * 
 * CosInclination will be adjusted to cos(inclination) of stellar spin w.r.t. l.o.s.
 *
 * @warning dust extinction was computed using old CosInclination (relative to z axis),
 *          but MagDust, etc. are not updated to reflect extiction 
 *          for actual inclination w.r.t. l.o.s.
 * 
 * if defined OUTPUT_FB_OBS_MAGS
 * ObsMag, etc. will be adjusted to absolute observer frame mags for galaxy
 * redshifted to ObsRedshift
 */
void
adjust_galaxy_for_lightcone(struct GALAXY_OUTPUT *galaxy_, const float shift_[3], const int shift_index_[3], const int output_number_)
{
#ifndef OUTPUT_FB_OBS_MAGS 
  (void)output_number_;  /* suppress unused-parameter warning */
#endif /* not defined OUTPUT_FB_OBS_MAGS */
 
  apply_shift(shift_, &(galaxy_->Pos));
  apply_cartesian_to_ra_dec_r(&(galaxy_->Pos));
  
  /** @todo express sin(ra) etc. in terms of cartesian Pos?  */
  float rot_m_[3][3];
  get_cartesian_to_ra_dec_r_local_orthogonal_rotation_matrix_from_ra_dec(galaxy_->Pos[0], galaxy_->Pos[1], &rot_m_);
  apply_rotation(rot_m_, &(galaxy_->DistanceToCentralGal));
  apply_rotation(rot_m_, &(galaxy_->Vel));
  apply_rotation(rot_m_, &(galaxy_->GasSpin));
  apply_rotation(rot_m_, &(galaxy_->StellarSpin));
   
  galaxy_->Redshift    = redshift_for_comoving_los_distance(galaxy_->Pos[2]);
  galaxy_->ObsRedshift = combine_redshifts(galaxy_->Redshift, redshift_for_radial_velocity(galaxy_->Vel[2]));
  
  galaxy_->CubeShiftIndex = composed_cube_shift_index(shift_index_);

  galaxy_->CosInclination = galaxy_->StellarSpin[2] / euclidian_norm(galaxy_->StellarSpin);

#ifndef LIGHT_OUTPUT
#ifdef HALOPROPERTIES
  apply_shift(shift_, &(galaxy_->HaloPos));
  apply_cartesian_to_ra_dec_r(&(galaxy_->HaloPos));
#endif  /* defined HALOPROPERTIES  */ 
#endif /* not defined LIGHT_OUTPUT  */ 

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_OBS_MAGS
#ifdef OUTPUT_FB_OBS_MAGS
  int filter_number_;

  const int   current_snapshot_number_   = ListOutputSnaps[output_number_];
  const float current_snapshot_redshift_ = ZZ[current_snapshot_number_];
  
  const int   earlier_snapshot_number_   = current_snapshot_number_ > 0 ? current_snapshot_number_ - 1 : 0;
  const float earlier_snapshot_redshift_ = ZZ[earlier_snapshot_number_];


  const int   later_snapshot_number_   = current_snapshot_number_ < LastDarkMatterSnapShot ? current_snapshot_number_ + 1 : LastDarkMatterSnapShot;
  const float later_snapshot_redshift_ = ZZ[later_snapshot_number_];

  if(galaxy_->ObsRedshift < current_snapshot_redshift_)
  {
    const float obs_z_rel_dev_ = (current_snapshot_redshift_ - galaxy_->ObsRedshift) / (current_snapshot_redshift_ - later_snapshot_redshift_);

    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      galaxy_->ObsMagDust [filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMagDust [filter_number_] + obs_z_rel_dev_ * galaxy_->forward_ObsMagDust [filter_number_]; 
      galaxy_->ObsMag     [filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMag     [filter_number_] + obs_z_rel_dev_ * galaxy_->forward_ObsMag     [filter_number_]; 
      galaxy_->ObsMagBulge[filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMagBulge[filter_number_] + obs_z_rel_dev_ * galaxy_->forward_ObsMagBulge[filter_number_]; 
#ifdef ICL
      galaxy_->ObsMagICL  [filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMagICL  [filter_number_] + obs_z_rel_dev_ * galaxy_->forward_ObsMagICL  [filter_number_]; 
#endif /* ICL */
    }
  }
  else /* galaxy_->ObsRedshift >= current_snapshot_redshift_ */
  {
    const float obs_z_rel_dev_ = (galaxy_->ObsRedshift - current_snapshot_redshift_) / (earlier_snapshot_redshift_ - current_snapshot_redshift_);
    for(filter_number_ = 0; filter_number_ < NMAG; filter_number_++)
    {
      galaxy_->ObsMagDust [filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMagDust [filter_number_] + obs_z_rel_dev_ * galaxy_->backward_ObsMagDust [filter_number_]; 
      galaxy_->ObsMag     [filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMag     [filter_number_] + obs_z_rel_dev_ * galaxy_->backward_ObsMag     [filter_number_]; 
      galaxy_->ObsMagBulge[filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMagBulge[filter_number_] + obs_z_rel_dev_ * galaxy_->backward_ObsMagBulge[filter_number_]; 
#ifdef ICL
      galaxy_->ObsMagICL  [filter_number_] = (1. - obs_z_rel_dev_) * galaxy_->ObsMagICL  [filter_number_] + obs_z_rel_dev_ * galaxy_->backward_ObsMagICL  [filter_number_]; 
#endif /* ICL */
    }
  }

#endif /* defined OUTPUT_FB_OBS_MAGS */
#endif /* defined OUTPUT_OBS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */
}


/** @brief outputs galaxy on lightcone
 * 
 * checks if galaxy in_ any of the periodic copies of the simulation box
 * is on the lightcone,
 * then transforms those periodic copies of the galaxy on the lightcone
 * and writes them to file
 */
void 
save_lightcone_galaxy_append(int galaxy_number_, int output_number_)
{
  lightcone_N_galaxies_for_output++;      
        
// #ifndef LIGHTCONE_OUTPUT_ONLY
  /* if closest copy of galaxy is outside lightcone, return early */
  /* if defined LIGHTCONE_OUTPUT_ONLY, check may have already been done earlier */
  if(is_outside_lightcone_for_snapshot[HaloGal[galaxy_number_].SnapNum])
  {
    lightcone_N_galaxies_skipped_output_early++;
    return;  
  }
  
  if(check_outside_lightcone_for_snapshot[HaloGal[galaxy_number_].SnapNum])
  {
    float d_0_ = HaloGal[galaxy_number_].Pos[0] - lightcone_observer_position[0]; d_0_ = wrap(d_0_, BoxSize);
    float d_1_ = HaloGal[galaxy_number_].Pos[1] - lightcone_observer_position[1]; d_1_ = wrap(d_1_, BoxSize);
    float d_2_ = HaloGal[galaxy_number_].Pos[2] - lightcone_observer_position[2]; d_2_ = wrap(d_2_, BoxSize);

    if(d_0_ * d_0_ + d_1_ * d_1_ + d_2_ * d_2_ > pow2(lightcone_radius_for_snapshot[HaloGal[galaxy_number_].SnapNum]))
    {
      lightcone_N_galaxies_skipped_output_early++;
      return;
    }
  }
// #endif /* not defined LIGHTCONE_OUTPUT_ONLY */     
  
  /* ignore galaxy not fitting stellar mass selection */
  if(HaloGal[galaxy_number_].BulgeMass + HaloGal[galaxy_number_].DiskMass < lightcone_lower_stellar_mass)
    return; 
  
  /** @note prepare_galaxy_for_output() is expensive, so seems worth avoiding multiple calls to this function
   *  by keeping a caopy after first required call for the same input galaxy */
  struct GALAXY_OUTPUT galaxy_output_, prepared_galaxy_output_;
  bool prepared_galaxy_output_is_valid_ = false;
  
  int box_shift_i_end_ = ceil((lightcone_slice_upper_los_distance[output_number_] + lightcone_observer_distance_from_origin) / BoxSize);
  
  int box_shift_i_[3];
  for(box_shift_i_[0] = -box_shift_i_end_; box_shift_i_[0] < box_shift_i_end_; box_shift_i_[0]++)
  for(box_shift_i_[1] = -box_shift_i_end_; box_shift_i_[1] < box_shift_i_end_; box_shift_i_[1]++)
  for(box_shift_i_[2] = -box_shift_i_end_; box_shift_i_[2] < box_shift_i_end_; box_shift_i_[2]++)
  {      
    const float shift_[3] = {box_shift_i_[0] * BoxSize - lightcone_observer_position[0],
                             box_shift_i_[1] * BoxSize - lightcone_observer_position[1],
                             box_shift_i_[2] * BoxSize - lightcone_observer_position[2]};

    /* even though testing for lightcone geometry as early as possible means
       possibly computing transformed positions and redshifts twice,
       this may save us from a lot of other computations done in 
       prepare_galaxy_for_output() and adjust_galaxy_for_lightcone().
      
       if this double computation becomes significant, 
       one may think about passing the already computed values to 
       adjust_galaxy_for_lightcone()
    */
    float lightcone_galaxy_position_[3];
    lightcone_galaxy_position_[0] = HaloGal[galaxy_number_].Pos[0] + shift_[0];
    lightcone_galaxy_position_[1] = HaloGal[galaxy_number_].Pos[1] + shift_[1];
    lightcone_galaxy_position_[2] = HaloGal[galaxy_number_].Pos[2] + shift_[2];
  
    const float los_distance_to_lightcone_observer_ = euclidian_norm(lightcone_galaxy_position_);
    
    if(los_distance_to_lightcone_observer_ < lightcone_slice_lower_los_distance[output_number_] ||
       los_distance_to_lightcone_observer_ > lightcone_slice_upper_los_distance[output_number_])
      continue;

    const float cosmological_redshift_ = redshift_for_comoving_los_distance(los_distance_to_lightcone_observer_);
    const float peculiar_redshift_     = redshift_for_radial_velocity(parallel_component(HaloGal[galaxy_number_].Vel, lightcone_galaxy_position_));
    const float observed_redshift_     = combine_redshifts(cosmological_redshift_, peculiar_redshift_);
    
    if(observed_redshift_ < lightcone_lower_redshift ||
       observed_redshift_ > lightcone_upper_redshift   ) 
      continue;
      
    apply_cartesian_to_ra_dec_r(&lightcone_galaxy_position_);
    
    if(lightcone_galaxy_position_[0] < lightcone_lower_ra  ||
       lightcone_galaxy_position_[0] > lightcone_upper_ra  ||
       lightcone_galaxy_position_[1] < lightcone_lower_dec ||
       lightcone_galaxy_position_[1] > lightcone_upper_dec   )
      continue;  
      
    if(!prepared_galaxy_output_is_valid_)
    { 
      prepare_galaxy_for_output(output_number_, &HaloGal[galaxy_number_], &prepared_galaxy_output_);
      prepared_galaxy_output_is_valid_ = true;
    }
    galaxy_output_ = prepared_galaxy_output_;
    
    /** @todo consider passing transformed position and redshifts to avoid recomputation */
    adjust_galaxy_for_lightcone(&galaxy_output_, shift_, box_shift_i_, output_number_);

#ifdef GALAXYTREE

#ifdef STAR_FORMATION_HISTORY
    galaxy_output_.sfh_numbins = galaxy_output_.sfh_ibin;
#endif /* defined STAR_FORMATION_HISTORY */

    myfwrite_lightcone_galaxy_from_galaxy_output(&galaxy_output_, 1, FdLightconeGalTree);
    TotLightconeGalCount++; //this will be written later
  
#else /* not defined GALAXYTREE */
    myfwrite_lightcone_galaxy_from_galaxy_output(&galaxy_output_, 1, FdLightconeGalDumps[output_number_]);
    TotLightconeGalaxies[output_number_]++;  //this will be written later
#endif /* not defined GALAXYTREE */
  }
}


/** @brief updates galaxies in lightcone galaxy output file to contain proper tree info
 * 
 * requires that all galaxies in_ tree have been constructed and that tree info in_ memory
 * has been updated
 * 
 * then uses tree info in memory to update lightcone galaxies on disk
 */
void save_lightcone_galaxy_finalize(int file_number_, int tree_number_)
{
#ifdef GALAXYTREE
  if(NGalTree > 0)
  {
    // order GalTree by current order of storage in_ file (lightcone_galaxy_number_in_file_begin)
    qsort(GalTree, NGalTree, sizeof(struct galaxy_tree_data), save_lightcone_galaxy_tree_compare);
  
#ifdef UPDATE_GALAXY_OUTPUT_IN_MEM
/* updating done in memory */

    const long long galaxy_in_file_number_begin_ = GalTree[0           ].lightcone_galaxy_number_in_file_begin;
    const long long galaxy_in_file_number_end_   = GalTree[NGalTree - 1].lightcone_galaxy_number_in_file_end;
    const long long N_galaxies_                  = galaxy_in_file_number_end_ - galaxy_in_file_number_begin_;

    lightcone_galaxy_output_type *galaxy_output_ = (lightcone_galaxy_output_type*) mymalloc("lc_file_gal", sizeof(lightcone_galaxy_output_type) * N_galaxies_);
    
    myfseek_lightcone_galaxy(FdLightconeGalTree, galaxy_in_file_number_begin_);
    myfread_lightcone_galaxy(galaxy_output_, N_galaxies_, FdLightconeGalTree);

    int galaxy_in_tree_number_;
    long long galaxy_in_file_number_;
    for(galaxy_in_tree_number_ = 0; galaxy_in_tree_number_ < NGalTree; galaxy_in_tree_number_++)
      for(galaxy_in_file_number_ = GalTree[galaxy_in_tree_number_].lightcone_galaxy_number_in_file_begin; galaxy_in_file_number_ < GalTree[galaxy_in_tree_number_].lightcone_galaxy_number_in_file_end; galaxy_in_file_number_++)
      { prepare_galaxy_tree_info_for_lightcone_output(file_number_, tree_number_, &GalTree[galaxy_in_tree_number_], &galaxy_output_[galaxy_in_file_number_ - galaxy_in_file_number_begin_]); }

#ifdef SORT_GALAXY_OUTPUT
    qsort(galaxy_output_, N_galaxies_, sizeof(lightcone_galaxy_output_type), lightcone_galaxy_compare);
#endif /* not defined SORT_GALAXY_OUTPUT */  
 
    myfseek_lightcone_galaxy(FdLightconeGalTree, galaxy_in_file_number_begin_);
    myfwrite_lightcone_galaxy(galaxy_output_, N_galaxies_, FdLightconeGalTree);

    myfree(galaxy_output_);

#else /* not defined UPDATE_GALAXY_OUTPUT_IN_MEM */
 /* updating done on disk with little memory overhead */

    lightcone_galaxy_output_type galaxy_output_;
    int galaxy_in_tree_number_;
    long long galaxy_in_file_number_;
    for(galaxy_in_tree_number_ = 0; galaxy_in_tree_number_ < NGalTree; galaxy_in_tree_number_++)
    { 
      for(galaxy_in_file_number_ = GalTree[galaxy_in_tree_number_].lightcone_galaxy_number_in_file_begin; galaxy_in_file_number_ < GalTree[galaxy_in_tree_number_].lightcone_galaxy_number_in_file_end; galaxy_in_file_number_++)
      {
        myfseek_lightcone_galaxy(FdLightconeGalTree, galaxy_in_file_number_);
        myfread_lightcone_galaxy(&galaxy_output_, 1, FdLightconeGalTree);
        prepare_galaxy_tree_info_for_lightcone_output(file_number_, tree_number_, &GalTree[galaxy_in_tree_number_], &galaxy_output_);
        myfseek_lightcone_galaxy(FdLightconeGalTree, galaxy_in_file_number_);
        myfwrite_lightcone_galaxy(&galaxy_output_, 1, FdLightconeGalTree);
      }
    }
  
#ifdef SORT_GALAXY_OUTPUT
#warning "Sorry, save_lightcone_galaxy_finalize() not yet implemented for SORT_GALAXY_OUTPUT without UPDATE_GALAXY_OUTPUT_IN_MEM."
/** @todo implement on-disk version for sorting galaxy output,
 *        taking inspiration from save_galaxy_tree_reorder_on_disk() */

#endif /* not defined SORT_GALAXY_OUTPUT */  
#endif /* not defined UPDATE_GALAXY_OUTPUT_IN_MEM */  
   /* after updating tree info in file, make sure that file position indicator is put back to end of file.
    * note: file position indicator should be at end of file already,
    * if GalTree[].lightcone_galaxy_number_in_file_begin/end work as expected an are propery sorted)
    */
    // myfseek_lightcone_galaxy(FdLightconeGalTree, TotLightconeGalCount);
  }
#else /* not defined GALAXYTREE */
 (void) file_number_;  /* avoid unused-parameter warning */
 (void) tree_number_; /* avoid unused-parameter warning */
#endif /* not defined GALAXYTREE */
}


/** @brief sorts lightcone galaxies in output files
  *
  * currently, this is a straightforward version reading all galaxies, sorting them in memory,
  * then writing back to disk.
  */
void sort_lightcone_galaxy_in_file(FILE * lightcone_galaxy_file_)
{
#ifdef UPDATE_GALAXY_OUTPUT_IN_MEM
  long long N_galaxies_ = 0;
  myfread_lightcone_galaxy_number_of_entries(lightcone_galaxy_file_, &N_galaxies_);

  lightcone_galaxy_output_type *galaxy_ = (lightcone_galaxy_output_type*) mymalloc("lc_file_gal", sizeof(lightcone_galaxy_output_type) * N_galaxies_);
  
  myfseek_lightcone_galaxy(lightcone_galaxy_file_, 0);
  myfread_lightcone_galaxy(galaxy_, N_galaxies_, lightcone_galaxy_file_);

  qsort(galaxy_, N_galaxies_, sizeof(lightcone_galaxy_output_type), lightcone_galaxy_compare);
  
  myfseek_lightcone_galaxy(lightcone_galaxy_file_, 0);
  myfwrite_lightcone_galaxy(galaxy_, N_galaxies_, lightcone_galaxy_file_);
  
#else /* not defined UPDATE_GALAXY_OUTPUT_IN_MEM */    
#warning "Sorry, sort_lightcone_galaxy_in_file() not yet implemented without UPDATE_GALAXY_OUTPUT_IN_MEM."
  (void)lightcone_galaxy_file_; /* avoid unused-parameter warning */
#endif /* not defined UPDATE_GALAXY_OUTPUT_IN_MEM */  
}


/** @brief compares galaxy_tree_data entries for sorting
 * 
 *  compares galaxy_tree_data entries for sorting for writing tree info data for lightcone output
 */
int save_lightcone_galaxy_tree_compare(const void *galaxy_tree_data_a_, const void *galaxy_tree_data_b_)
{
  if(((struct galaxy_tree_data *) galaxy_tree_data_a_)->lightcone_galaxy_number_in_file_begin < ((struct galaxy_tree_data *) galaxy_tree_data_b_)->lightcone_galaxy_number_in_file_begin)
    return -1;

  else if(((struct galaxy_tree_data *) galaxy_tree_data_a_)->lightcone_galaxy_number_in_file_begin > ((struct galaxy_tree_data *) galaxy_tree_data_b_)->lightcone_galaxy_number_in_file_begin)
    return +1;
  
  else if(((struct galaxy_tree_data *) galaxy_tree_data_a_)->lightcone_galaxy_number_in_file_end < ((struct galaxy_tree_data *) galaxy_tree_data_b_)->lightcone_galaxy_number_in_file_end)
    return -1;

  else if(((struct galaxy_tree_data *) galaxy_tree_data_a_)->lightcone_galaxy_number_in_file_end > ((struct galaxy_tree_data *) galaxy_tree_data_b_)->lightcone_galaxy_number_in_file_end)
    return +1;
  
  else
    return 0;
}


/** @brief compares lightcone galaxy entries for sorting
 * 
 *  compares lightcone_galaxy_output_type entries for sorting output on disk
 */
int lightcone_galaxy_compare(const void *lightcone_galaxy_a_, const void *lightcone_galaxy_b_)
{
#ifdef GALAXYTREE 
  /* if GalID available, use GalID for sorting */
  if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->GalID < ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->GalID)
    return -1;

  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->GalID > ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->GalID)
    return +1;
#endif /* defined GALAXYTREE */

  /* next, use CubeShiftIndex for sorting (if GalID is available, (GalID,CubeShiftIndex) provide total order) */
  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->CubeShiftIndex < ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->CubeShiftIndex)
    return -1;

  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->CubeShiftIndex > ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->CubeShiftIndex)
    return +1;
  
  /* next, use SnapNum for sorting (reverse order) */
  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->SnapNum > ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->SnapNum)
    return -1;

  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->SnapNum < ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->SnapNum)
    return +1;

  /* next, use obs. redshift */
  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->ObsRedshift < ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->ObsRedshift)
    return -1;

  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->ObsRedshift > ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->ObsRedshift)
    return +1;

  /* next, use StellarMass (reverse order)*/
  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->StellarMass > ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->StellarMass)
    return -1;

  else if(((lightcone_galaxy_output_type*) lightcone_galaxy_a_)->StellarMass < ((lightcone_galaxy_output_type*) lightcone_galaxy_b_)->StellarMass)
    return +1;

  else 
    return 0;
}

/****************************************************/
#endif /* defined LIGHTCONE_OUTPUT */
/*** nothing here unless defined LIGHTCONE_OUTPUT ***/
