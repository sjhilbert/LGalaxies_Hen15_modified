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
/** @file sam.h
 * @brief SAM Construct Galaxies, Join Galaxies of progenitors, Evolve Galaxies
 * */


#ifndef SAM_H
#define SAM_H

/**@brief SAM() loops on trees and calls construct_galaxies.*/
double SAM(const int filenr);


/** @brief  construct_galaxies() recursively runs the semi-analytic model. */
void construct_galaxies(const int treenr, const int halonr);


/** @brief updates the properties of the galaxy from the dark matter halo
 *  properties and deals with merging clocks. */
int join_galaxies_of_progenitors(const int halonr, const int ngalstart, int *cenngal)


/** @brief evolve_galaxies() deals with most of the SA recipes.
  * 
  * @note halonr is here the FOF-background subhalo (i.e. main halo) */
void evolve_galaxies(const int halonr, const int ngal, const int treenr, const int cenngal);


/** @brief output_galaxy() outputs galaxy to various files on disk */
void output_galaxy(const int treenr, const int heap_index);

#endif /* header guard */
