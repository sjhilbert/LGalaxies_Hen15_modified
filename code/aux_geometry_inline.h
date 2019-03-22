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

/** @file aux_geometry_inline.h 
 * 
 * @author Stefan Hilbert (hilbert)
 * 
 */
#ifndef AUX_GEOMETRY_INLINE_H
#define AUX_GEOMETRY_INLINE_H


/** @brief scalar product for 3d vectors cartesian vectors */
static inline float 
scalar_product_3d(const float v_a_[3], const float v_b_[3])
{ return v_a_[0] * v_b_[0] + v_a_[1] * v_b_[1] + v_a_[2] * v_b_[2]; }


/** @brief returns square of euclidian norm of cartesian vector */
static inline float
euclidian_norm_squared_3d(const float v_[3])
{ return v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]; }


/** @brief returns euclidian norm of cartesian vector */
static inline float
euclidian_norm_3d(const float v_[3])
{ return sqrt(v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]); }


/** @brief returns euclidian distance between two cartesian positions  */
static inline float
euclidian_distance_between_3d(const float pos_a_[3], const float pos_b_[3])
{ return sqrt((pos_b_[0] - pos_a_[0]) * (pos_b_[0] - pos_a_[0]) + (pos_b_[1] - pos_a_[1]) * (pos_b_[1] - pos_a_[1]) + (pos_b_[2] - pos_a_[2]) * (pos_b_[2] - pos_a_[2])); }


/** @brief returns component parallel to a given direction */
static inline float
parallel_component_3d(const float v_[3], const float dir_[3])
{ return scalar_product_3d(v_, dir_) / euclidian_norm_3d(dir_); }


/** @brief shifts cartesian vector by shift_ vecor
 *
 * @param[in]      shift_ 3d shift_ vector
 * @param[in,out]  pos_   3d vector to be shifted
 */
static inline void
apply_shift_3d(const float shift_[3], float (*pos_)[3])
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
apply_rotation_3d(/* const */ float rot_m_[3][3], float (*v_)[3])
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
  const float r   = euclidian_norm_3d(*pos_);
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

  
#endif /* header guard */
