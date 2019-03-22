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

#ifndef PROTO_H
#define PROTO_H

#include "allvars.h"

#include "aux_functions_inline.h"

#include "metals_inline.h"

#include <math.h>
#include <stdbool.h>

double SAM(const int filenr);

void construct_galaxies(const int tree, const int halonr);
int join_galaxies_of_progenitors(const int halonr, const int ngalstart, int *cenngal);
void evolve_galaxies(const int halonr, const int ngal, const int tree, const int cenngal);

size_t myfread (void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t myfwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t myffill (void *ptr, size_t size, size_t nmemb, FILE * stream);
int myfseek(FILE * stream, long offset, int whence);

void load_all_auxdata(const int filenr);

int peano_hilbert_key(const int x, const int y, const int z,const  int bits);
void update_type_two_coordinate_and_velocity(const int tree, const int i, const int centralgal);
void output_galaxy(const int treenr, const int heap_index);

void close_galaxy_tree_file(void);
void create_galaxy_tree_file(const int filenr);
void save_galaxy_tree_append(const int i);
void update_galaxy_tree_ids(void);
void save_galaxy_tree_finalize(const int filenr, const int tree);
void prepare_galaxy_tree_info_for_output(const int filenr, const int tree, const struct galaxy_tree_data *g, struct GALAXY_OUTPUT *o);
int walk_galaxy_tree(const int nr);
int save_galaxy_tree_compare(const void *a, const void *b);

void save_galaxy_append(const int tree, const int i, const int n);
void close_galaxy_files(void);
void create_galaxy_files(const int filenr);

long long calc_big_db_subid_index(const int snapnum, const int filenr, const int subhaloindex);

void dump_memory_table(void);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, const int line);
void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, const int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, const int line);
void myfree_fullinfo(void *p, const char *func, const char *file, const int line);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, const int line);
void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, const int line);
void report_detailed_memory_usage_of_largest_task(size_t * OldHighMarkBytes, const char *label, const char *func, const char *file, const int line);
void mymalloc_init(void);

void save_galaxy_tree_reorder_on_disk(void);
int save_galaxy_tree_mp_comp(const void *a, const void *b);

void get_coordinates(float *pos, float *vel, const long long ID, const int tree, const int halonr, const int snapnum);


//functions used to scale to a different cosmology
void read_scaling_parameters();
void scale_cosmology(const int nhalos);
void un_scale_cosmology(const int nhalos);
void read_zlist_original_cosm(void);
void read_zlist_new(void);
double c_correction(const float mass, const int snapnum);
double find_c(const double c_ori, const double ratio);
double func_c(const double c);
double func_c_p(const double c);
double dgrowth_factor_dt(double a, double omega_m, double omega_l);
double scale_v_cen(const int snapnum);

long long calc_big_db_offset(const int filenr, const int treenr);

void init(void);
void read_parameter_file(char *fname);
void check_program_parameters();
// void check_compile_time_options();
void set_units(void);
#ifdef SPECIFYFILENR
void read_file_nrs(void);
#endif /* defined SPECIFYFILENR */
int get_nr_files_to_process(void);
void assign_files_to_tasks(int *FileToProcess, int *TaskToProcess, const int nfiles);
void read_reionization(void);


void load_tree_table(const int filenr);
void load_tree(const int nr);
void prepare_galaxy_for_output(const int n, const struct GALAXY *g, struct GALAXY_OUTPUT *o);
void fix_units_for_ouput(struct GALAXY_OUTPUT *o);

void free_galaxies_and_tree(void);
void free_tree_table(void);
void endrun(const int ierr);

void finalize_galaxy_file(const int filenr);

void starformation(const int p, const int centralgal, const double time, const double dt, const int nstep);
void update_stars_due_to_reheat(const int p, const int centralgal, double *stars);
void update_from_star_formation(const int p, const double stars, const bool flag_burst, const int nstep);
void SN_feedback(const int p, const int centralgal, const double stars, const GasComponentType feedback_location);
void update_from_feedback(const int p, const int centralgal, const double reheated_mass, double ejected_mass);
void update_massweightage(const int p, const double stars, const double time);

void add_galaxies_together(const int t, const int p);
void init_galaxy(const int p, const int halonr);
double infall_recipe(const int centralgal, const int ngal, const double Zcurr);
void add_infall_to_hot(const int centralgal, const double infallingGas);
void compute_cooling(const int p, const double dt);
void do_AGN_heating(const double dt, const int ngal);
void cool_gas_onto_galaxy(const int p, const double dt);
void reincorporate_gas(const int p, const double dt);
void deal_with_galaxy_merger(const int p, const int merger_centralgal, const int centralgal, const double time, const double deltaT, const int nstep);
double get_reionization_modifier(const float Mvir, const double Zcurr);

//MATH MISC
void locate(double *xx, const int n, const double x, int *j);
double integrate(double *flux, const int Grid_Length);
void polint(double xa[], double ya[], const int n, const double x, double *y, double *dy);
void nrerror(char error_text[]);
double* new_vector(const long nl, const long nh);
void free_vector(double *v, const long nl, const long nh);

//SPECTRO/PHOTOMETRY PROPERTIES
#ifdef COMPUTE_SPECPHOT_PROPERTIES

#ifdef PHOTTABLES_PRECOMPUTED
void setup_LumTables_precomputed(char SimName[]);
#endif /* defined PHOTTABLES_PRECOMPUTED */

#ifdef SPEC_PHOTABLES_ON_THE_FLY

void setup_Spec_LumTables_onthefly(void);

//initialize
void setup_RedshiftTab();
void read_vega_spectra(double *LambdaVega, double *FluxVega);
void read_filters(double LambdaFilter[NMAG][MAX_NLambdaFilter], double FluxFilter[NMAG][MAX_NLambdaFilter]);
void read_MetalTab();
void read_InputSSP_spectra(double LambdaInputSSP[SSP_NAGES][SSP_NLambda], double FluxInputSSP[SSP_NAGES][SSP_NLambda], int MetalLoop);
void free_input_spectra(double LambdaInputSSP[SSP_NAGES][SSP_NLambda], double FluxInputSSP[SSP_NAGES][SSP_NLambda]);

//misc
double get_AbsAB_magnitude(const double FluxInputSSPInt, const double FluxFilterInt, const double redshift);
double get_area (const double redshift);
double lum_distance(const double redshift);

//numerical
double* create_grid (const double WaveMin, const double WaveMax, const int AgeLoop, const double redshift, double LambdaInputSSP[SSP_NAGES][SSP_NLambda],
		                      int *Min_Wave_Grid, int *Max_Wave_Grid, int *Grid_Length);
void interpolate(double *lgrid, const int Grid_Length, double *lambda, const int nlambda, double *flux, double *FluxOnGrid);
#endif /* defined SPEC_PHOTABLES_ON_THE_FLY */

#ifdef POST_PROCESS_MAGS
void post_process_spec_mags(struct GALAXY_OUTPUT *o);
#else /* not defined POST_PROCESS_MAGS */
void add_to_luminosities(const int p, double mstars, double time, double dt, const double metallicity);
void dust_model(const int p, const int snap, const int halonr);
#endif /* not defined POST_PROCESS_MAGS */

// dust model
void read_dust_tables(void);
double get_extinction(const int mag, const double Zg, const double redshift);

// luminosity table lookup speedup:
void init_SSP_log_age_jump_index(void);

static inline int 
get_SSP_log_age_jump_index(const double log10_age_)
{ return SSP_log_age_jump_table[(int) ((log10_age_ - SSP_logAgeTab[1]) * SSP_log_age_jump_factor)]; }

#define find_age_luminosity_interpolation_parameters(log10_age_, age_index_, f_age_1_, f_age_2_) \
locate_interpolation_index_and_fraction_bf(age_index_, 0, SSP_NAGES, log10_age_, SSP_logAgeTab, f_age_1_, f_age_2_, <, linear, get_SSP_log_age_jump_index)

#define find_metallicity_luminosity_interpolation_parameters(log10_metallicity_, met_index_, f_met_1_, f_met_2_) \
locate_interpolation_index_and_fraction(met_index_, 0, SSP_NMETALLICITES, log10_metallicity_, SSP_logMetalTab, f_met_1_, f_met_2_, <, linear)

#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

float gasdev(long *idum);

double estimate_merging_time(const int halonr, const int mostmassive, const int p);

double get_hubble_parameter_for_halo(const int halonr);
double get_hubble_parameter_squared_for_halo(const int halonr);
double get_virial_velocity(const int halonr);
double get_virial_radius(const int halonr);
double get_virial_mass(const int halonr);
double collisional_starburst_recipe(const double mass_ratio, const int merger_centralgal, const int centralgal, const double time, const double deltaT);

void make_bulge_from_burst(const int p);
void grow_black_hole(const int merger_centralgal, const double mass_ratio, const double deltaT);
void check_disk_instability(const int p);

double get_disk_radius(const int halonr, const int p);
void get_gas_disk_radius(const int p);
void get_stellar_disk_radius(const int p);

void read_output_snaps(void);
void read_zlist(void);


void read_cooling_functions(void);
double get_metaldependent_cooling_rate(const double logTemp, const double logZ);


void disrupt(const int p);

double hot_retain_sat(const int i, const int centralgal);

double get_initial_disk_radius(const int halonr, const int p);
void update_bulge_from_disk(const int p, const double stars);
double bulge_from_disk(const double frac);
double func_size(const double x, const double a);
void bulgesize_from_merger(const double mass_ratio, const int merger_centralgal, const int p, const double Mcstar, const double Mcbulge, const double Mcgas, const  double Mpstar, const double Mpbulge, const double Mpgas, double frac);

void update_type_2(const int ngal, const int halonr, const int prog, int mostmassive);
void update_centralgal(const  int ngal, const int halonr);
void update_hotgas(const int ngal, const int centralgal);
void update_type_1(const int ngal, const int halonr,const int prog);
double separation_gal(const int p, const int q);
double separation_halo(const int p, const int q);

/** @bug can't find implementation for transfer_ICL() */
void transfer_ICL(const int p, const int q, const double fraction); 

/* check funs */
void update_hot_frac(const int p, const double reincorporated, const float HotGas);
int get_merger_center(const int fofhalo);

void transfer_stars(const int p, const StellarComponentType cp, const int q, StellarComponentType cq, const double fraction);
void transfer_gas(const int p, const GasComponentType cp, const int q, const GasComponentType cq, const double fraction);
void deal_with_satellites(int centralgal, int ngal);
void perform_mass_checks(char string[], int igal) ;

#ifdef STAR_FORMATION_HISTORY
void sfh_initialise(const int p);
void sfh_merge(const int p, const int p1);
void sfh_print(const int p);
void sfh_update_bins(const int p, const int snap, const int step, const double time);
void create_sfh_bins();
void write_sfh_bins();
#endif /* defined STAR_FORMATION_HISTORY */ 

#ifdef DETAILED_METALS_AND_MASS_RETURN
//in read_yield_tables.c:
void read_yield_tables();
double Chabrier_IMF(const double M);

//in calc_SNe_rates.c:
#ifdef INDIVIDUAL_ELEMENTS
void SNe_rates();
#endif /* defined INDIVIDUAL_ELEMENTS */

//in yield_integrals.c:
void init_integrated_yields();
void integrate_yields();
int find_initial_metallicity_comp(int Zi, int sfh_bin, int table_type);
int find_initial_mass(double lifetime, int Zi_bin);
int max_Mi_lower(int Mi_lower, int channel_type);
int min_Mi_upper(int Mi_upper, int channel_type);
int find_SNII_mass_bin(double masslimit);
int find_agb_mass_bin(double masslimit);
#ifdef DTD
double DTDcalc (double timevalue);
#endif /* defined DTD */
#ifdef INDIVIDUAL_ELEMENTS
void find_actual_ejecta_limits(const int channel_type, const double Mi_lower_actual, const double Mi_upper_actual, const int Mi_lower, const int Mi_upper, const int Zi,
		double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual,
		double* Yields_lower_actual, double* Yields_upper_actual);
#else /* not defined INDIVIDUAL_ELEMENTS */
void find_actual_ejecta_limits(const int channel_type, const double Mi_lower_actual, const double Mi_upper_actual, const int Mi_lower, const int Mi_upper, const int Zi,
		double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual);
#endif /* not defined INDIVIDUAL_ELEMENTS */

void print_galaxy(char string[], const int p, const int halonr);

//in elements.c:
struct elements elements_init();
struct elements elements_add(const struct elements ele1, const struct elements ele2);
struct elements elements_fraction(const struct elements ele2, const float fraction);
struct elements elements_add_fraction(const struct elements ele1, const struct elements ele2, const float fraction);
void elements_add_to(struct elements *ele1, const struct elements ele2);
void elements_add_fraction_to(struct elements *ele1, const struct elements ele2, const float fraction);
void elements_deduct_from(struct elements *ele, const struct elements ele2);
void elements_deduct_fraction_from(struct elements *ele, const struct elements ele2, const float fraction);
void elements_print(char s[],const struct elements ele);
double elements_total(const struct elements ele);
double metal_elements_total(const struct elements ele);

//in recipe_yields.c:
void update_yields_and_return_mass(const int p, const int centralgal, const double dt, const int nstep);
int find_initial_metallicity(const int p, const int sfh_bin, const int table_type, const int component_type);
#ifdef INSTANTANEOUS_RECYCLE
void reset_ejection_rates(const int i, const int sfh_ibin,
		 double *NormSNIIMassEjecRate_actual, double *NormSNIIMetalEjecRate_actual,
		 double *NormSNIaMassEjecRate_actual, double *NormAGBMassEjecRate_actual,
		 double *NormSNIaMetalEjecRate_actual, double *NormAGBMetalEjecRate_actual);
#endif /* defined INSTANTANEOUS_RECYCLE */

#endif /* defined DETAILED_METALS_AND_MASS_RETURN */

/* aux. functions for cosmology */
#ifdef ASSUME_FLAT_LCDM
void assert_flat_LCDM(void);
#endif /* defined ASSUME_FLAT_LCDM */
void init_redshift_for_comoving_distance(void);
double comoving_los_distance_for_redshift(const double redshift_);
double redshift_for_comoving_los_distance(const double d_);
double redshift_for_radial_velocity(const double v_);
double combine_redshifts(const double z_1_, const double z_2_);
double redshift_for_comoving_los_distance_and_radial_velocity(const double d_, const double v_);
double time_to_present(const double redshift_);

/** @brief snapshot number to cosmic time */
static inline double 
NumToTime(int snapnum)
{ return Age[snapnum]; }

/** @brief Converts luminosities into magnitudes
 *
 * Converts luminosities into magnitudes:
 * \f$ M=-2.5\mathrm{log}_{10}(L) \f$ */
static inline double 
lum_to_mag(double lum_)
{ return (lum_ > 0) ? -2.5 * log10(lum_) : 99.0; }

#ifdef LIGHTCONE_OUTPUT
void init_lightcone(void);
void show_lightcone_statistics(void);
void create_lightcone_galaxy_files(const int file_number_);
void close_lightcone_galaxy_files(void);
void adjust_galaxy_for_lightcone(struct GALAXY_OUTPUT *galaxy_, const float shift_[3], const int shift_index_[3], const int output_number_);
void save_lightcone_galaxy_append(const int galaxy_number_, const int output_number_);
void save_lightcone_galaxy_finalize(const int file_number_, const int tree_number_);
void sort_lightcone_galaxy_in_file(FILE *lightcone_galaxy_file_);
int save_lightcone_galaxy_tree_compare(const void *galaxy_tree_data_a_, const void *galaxy_tree_data_b_);
int lightcone_galaxy_compare(const void *lightcone_galaxy_a_, const void *lightcone_galaxy_b_);
#endif /* defined LIGHTCONE_OUTPUT */


#endif /* header guard */
