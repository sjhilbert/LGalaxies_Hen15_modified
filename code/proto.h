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

/* cosmology.c */
#ifdef ASSUME_FLAT_LCDM
void assert_flat_LCDM(void);
#endif /* defined ASSUME_FLAT_LCDM */

void init_redshift_for_comoving_distance(void);
void init_cosmology_gsl_integration(void);
void init_cosmology(void);

double comoving_los_distance_for_redshift(const double redshift_);
double luminosity_distance_for_redshift(const double redshift_);
double redshift_for_comoving_los_distance(const double distance_);
double redshift_for_radial_velocity(const double velocity_);
double combine_redshifts(const double redshift_1_, const double redshift_2_);
double redshift_for_comoving_los_distance_and_radial_velocity(const double distance_, const double velocity_);
double time_to_present(const double redshift_);


/** @brief snapshot number to cosmic time */
#define NumToTime(snapshot_number_) Age[snapshot_number_]


/** @brief Converts luminosities into magnitudes 
 *         (unless #defined FULL_SPECTRA).
 *
 * Converts luminosities into magnitudes:
 * \f$ M=-2.5\mathrm{log}_{10}(L) \f$ */
static inline double 
lum_to_mag(const double lum_)
{
#ifdef FULL_SPECTRA  
  return lum_;
#else  /* not defined FULL_SPECTRA */
  return (lum_ > 2.511886431509581e-40) ? -2.5 * log10(lum_) : 99.0; 
#endif /* not defined FULL_SPECTRA */
}


/* init.c */
void init(void);

void read_zlist(void);
void read_zlist_new(void);
void read_zlist_original_cosm(void);
void read_output_snaps(void);
void read_reionization(void);

#ifdef SPECIFYFILENR
void read_file_nrs(void);
#endif /* defined SPECIFYFILENR */

int get_nr_files_to_process();
void assign_files_to_tasks(int *FileToProcess_, int *TaskToProcess_, int n_files_);

/* io_tree.c */
void load_tree_table(const int file_number_);
void free_tree_table(void);
void load_tree(const int tree_number_);
void free_galaxies_and_tree(void);

size_t myfread(void *ptr, const size_t size_, const size_t n_memb_, FILE * stream_);
size_t myfwrite(void *ptr, const size_t size_, const size_t n_memb_, FILE * stream_);
size_t myfwrite_large_data(void *ptr, const size_t size_, const size_t n_memb_, FILE * stream_);
size_t myffill(void *ptr_, const size_t size_, const size_t n_memb_, FILE * stream_);
int myfseek(FILE * stream_, const long offset_, const int whence_);

/* model.c */
double SAM(const int tree_file_number_);
void construct_galaxies(const int tree_number_, const int halo_number_);
void construct_galaxies_in_fof(const int tree_number_, const int first_in_fof_halo_number_);
void join_galaxies_of_progenitors(const int halo_number_, int *n_galaxies_in_fof_group_, int *merger_center_);
void evolve_galaxies(const int halo_number_, const int n_galaxies_in_fof_group_, const int tree_number_, const int merger_center_);
void output_galaxy(const int tree_number_, const int heap_index_);

/* model_cooling.c */
void read_cooling_functions(void);

double get_metaldependent_cooling_rate(const double log10_T_, const double log10_Z_);
void test_metaldependent_cooling_rate(void);

void compute_cooling(const int galaxy_number_, const double dt_);
void do_AGN_heating(double dt_, const int n_galaxies_in_fof_group_);
void cool_gas_onto_galaxy(int galaxy_number_, double dt_);

/* model_disrupt.c */
void disrupt(const int galaxy_number_);

/* model_dust.c */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
void dust_model(const int galaxy_number_, const int output_number_);
#endif /* not defined POST_PROCESS_MAGS */
#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

/* model_infall.c */
double infall_recipe(const int central_galaxy_number_, const int n_galaxies_, double current_redshift_);
double get_reionization_modifier(const float M_vir_, const double current_redshift_);
void add_infall_to_hot(const int central_galaxy_number_, const double infalling_gas_);

/* model_mergers.c */
int get_merger_center(const int halo_number_);

double estimate_merging_time(const int halo_number_, const int mother_halo_number_, const int galaxy_number_);
void deal_with_galaxy_merger(const int galaxy_number_, const int merger_centralgal_, const int centralgal_, const double time_, const double deltaT_);
void grow_black_hole(const int merger_centralgal_, const double mass_ratio_, const double deltaT_);
void add_galaxies_together(const int central_galaxy_number_, const int satellite_galaxy_number_);
void make_bulge_from_burst(const int galaxy_number_);
double collisional_starburst_recipe(const double mass_ratio_, const int merger_centralgal_, const int centralgal_,  const double time_, const double deltaT_);

void  bulgesize_from_merger(const double mass_ratio_, const int merger_centralgal_, const int galaxy_number_,
                            const double M_c_star_, const double M_c_bulge_, const double M_c_gas_,
                            const double M_p_star_, const double M_p_bulge_, const double M_p_gas_, double fraction_);

/* model_misc.c */
double separation_gal(const int galaxy_number_a_, const int galaxy_number_b_);
double separation_halo(const int halo_number_a_, const  int halo_number_b_);
double get_disk_radius(const int halo_number_, const int galaxy_number_);
double get_initial_disk_radius(const int halo_number_, const int galaxy_number_);
double get_virial_mass(const int halo_number_);
double get_virial_velocity(const int halo_number_);
double get_hubble_parameter_for_halo(const int halo_number_);
double get_hubble_parameter_squared_for_halo(const int halo_number_);
double get_virial_radius(const int halo_number_);

void set_gas_disk_radius(const int galaxy_number_);
void set_stellar_disk_radius(const int galaxy_number_);

void init_galaxy(const int galaxy_number_, const int halo_number_);

void update_centralgal(const int galaxy_number_, const int halo_number_);
void update_type_1(const int galaxy_number_, const int halo_number_, const int progenitor_halo_number_);
void update_type_2(const int galaxy_number_, const int halo_number_, const int progenitor_halo_number_, int most_massive_halo_number_);

void transfer_stars(const int to_galaxy_number_, const StellarComponentType to_component_, const int from_galaxy_number_, const StellarComponentType from_component_, const double fraction_);
void transfer_gas(const int to_galaxy_number_, const GasComponentType to_component_, const int from_galaxy_number_, const GasComponentType from_component_, const double fraction_);

/** @bug can't find implementation for transfer_ICL() */
void transfer_ICL(const int p, const int q, const double fraction); 

void perform_mass_checks(char tag_string_[], const int galaxy_number_);
void print_galaxy(char tag_string_[], const int galaxy_number_, const int halo_number_);

/* model_reincorporation.h */
void reincorporate_gas(const int galaxy_number_, const double dt_);

#ifdef COMPUTE_SPECPHOT_PROPERTIES
/* model_spectro_photometric.c */

#ifdef PHOTTABLES_PRECOMPUTED
void setup_LumTables_precomputed(const char sim_name_[]);
#endif /* defined PHOTTABLES_PRECOMPUTED */

void init_SSP_log_age_jump_index(void);

// static inline int get_SSP_log_age_jump_index(const double log10_age_)
// { return SSP_log_age_jump_table[(int) ((log10_age_ - SSP_logAgeTab[1]) * SSP_log_age_jump_factor)]; }

#define get_SSP_log_age_jump_index(log10_age_) \
SSP_log_age_jump_table[(int) ((log10_age_ - SSP_logAgeTab[1]) * SSP_log_age_jump_factor)]

#define find_age_luminosity_interpolation_parameters(log10_age_, age_index_, f_age_1_, f_age_2_) \
locate_interpolation_index_and_fraction_bf(age_index_, 0, SSP_NAGES, log10_age_, SSP_logAgeTab, f_age_1_, f_age_2_, <, linear, get_SSP_log_age_jump_index)

#define find_metallicity_luminosity_interpolation_parameters(log10_metallicity_, met_index_, f_met_1_, f_met_2_) \
locate_interpolation_index_and_fraction(met_index_, 0, SSP_NMETALLICITES, log10_metallicity_, SSP_logMetalTab, f_met_1_, f_met_2_, <, linear)

#ifndef  POST_PROCESS_MAGS
void add_to_luminosities(const int galaxy_number_, double stellar_mass_, double time_, double dt_, const double metallicity_);
#endif /* POST_PROCESS_MAGS */

/* model_spectro_photometric_on_the_fly.c */
#ifdef SPEC_PHOTABLES_ON_THE_FLY
void setup_RedshiftTab();
void read_vega_spectra(double *LambdaVega_, double *FluxVega_);
void read_filters(double LambdaFilter_[NMAG][MAX_NLambdaFilter], double FluxFilter_[NMAG][MAX_NLambdaFilter]);
void read_MetalTab();
void read_InputSSP_spectra(double LambdaInputSSP_[SSP_NAGES][SSP_NLambda], double FluxInputSSP_[SSP_NAGES][SSP_NLambda], const int met_index_);
void setup_Spec_LumTables_onthefly(void);
#endif /* defined SPEC_PHOTABLES_ON_THE_FLY */

#endif /* defined COMPUTE_SPECPHOT_PROPERTIES */

/* model_starformation_and_feedback.c */
void starformation(const int galaxy_number_, const int central_galaxy_number_, const double time_, const double dt_);
void update_stars_due_to_reheat(const int galaxy_number_, double *stars_);
void update_from_star_formation(const int galaxy_number_, const double stars_, const bool flag_burst_);
void SN_feedback(const int galaxy_number_, const int fof_central_galaxy_number_, const double stars_, const GasComponentType feedback_location_);
void update_from_feedback(const int galaxy_number_, const int central_galaxy_number_, const double reheated_mass_, double ejected_mass_);
void update_massweightage(const int galaxy_number_, const double stars_, const double time_);
void check_disk_instability(const int galaxy_number_);
void update_bulge_from_disk(const int galaxy_number_, const double stars_);

/* model_stripping.c */
double hot_retain_sat(const int galaxy_number_, const int central_galaxy_number_);
void deal_with_satellites(const int central_galaxy_number_, const int n_galaxies_);

/* model_yields.c */
void update_yields_and_return_mass(int galaxy_number_, int central_galaxy_number_, double dt_, int step_number_);
int find_initial_metallicity(int galaxy_number_, int sfh_bin_number_, int table_type_, StellarComponentType component_);

#ifdef INSTANTANEOUS_RECYCLE
void reset_ejection_rates(int sfh_bin_number_, int sfh_ibin_,
		                      double *NormSNIIMassEjecRate_actual_, double *NormSNIIMetalEjecRate_actual_,
		                      double *NormSNIaMassEjecRate_actual_, double *NormAGBMassEjecRate_actual_,
		                      double *NormSNIaMetalEjecRate_actual_, double *NormAGBMetalEjecRate_actual_);
#endif /* defined INSTANTANEOUS_RECYCLE */

/* mymalloc.c */
void mymalloc_init(void);

void report_detailed_memory_usage_of_largest_task(size_t * OldHighMarkBytes_, const char *label_, const char *func_, const char *file_, const int line_);
void dump_memory_table(void);

void *mymalloc_fullinfo(const char *var_name_, size_t n_, const char *func_, const char *file_, const int line_);
void *mymalloc_movable_fullinfo(void *ptr_, const char *var_name_, size_t n_, const char *func_, const char *file_, const int line_);
void myfree_fullinfo(void *ptr_, const char *func_, const char *file_, const int line_);
void myfree_movable_fullinfo(void *ptr_, const char *func_, const char *file_, const int line_);
void *myrealloc_fullinfo(void *ptr_, size_t n_, const char *func_, const char *file_, const int line_);
void *myrealloc_movable_fullinfo(void *ptr_, size_t n_, const char *func_, const char *file_, const int line_);

void endrun(const int error_number_);

/* peano.c */
int peano_hilbert_key(const int x, const int y, const int z, const int bits);

/* post_process_spec_mags.c */
void post_process_spec_mags(struct GALAXY_OUTPUT *galaxy_);

/* read_parameters.c */
void read_parameter_file(char *file_name_);
void check_program_parameters();
void compute_derived_program_parameters(void);

/* save.c */
void create_galaxy_files(const int file_number_);
void save_galaxy_append(const int tree_number_, const int galaxy_number_, const int output_number_);
void close_galaxy_files(void);

void prepare_galaxy_for_output(const int output_number_, const struct GALAXY *galaxy_, struct GALAXY_OUTPUT *output_galaxy_);

long long calc_big_db_offset(const int file_number_, const int tree_number_);
long long calc_big_db_subid_index(const int snapnum, const int file_number_, const int subhalo_index_);

#ifdef OUTPUT_BUFFERING
void save_galaxy_init_output_buffer(void);
void save_galaxy_flush_output_buffer(void);
void save_galaxy_show_output_buffer_statistics(void);
#endif /* defined OUTPUT_BUFFERING */

#ifdef FIX_OUTPUT_UNITS
void fix_units_for_ouput(struct GALAXY_OUTPUT *output_galaxy_);
#endif /* defined FIX_OUTPUT_UNITS */

#ifdef STAR_FORMATION_HISTORY
void write_sfh_bins();
#endif /* defined STAR_FORMATION_HISTORY */

#ifdef GALAXYTREE
/* save_galtree.c */
void create_galaxy_tree_file(const int file_number_);
void save_galaxy_tree_append(const int galaxy_number_);
void save_galaxy_tree_finalize(const int file_number_, const int tree_number_);
void close_galaxy_tree_file(void);
void save_galaxy_tree_reorder_on_disk(void);

void prepare_galaxy_tree_info_for_output(const int file_number_, const int tree_number_, const struct galaxy_tree_data *galaxy_, struct GALAXY_OUTPUT *output_galaxy_);

void update_galaxy_tree_ids(void);
int walk_galaxy_tree(const int galaxy_number_);

int output_galaxy_compare(const void *output_galaxy_a_, const void *output_galaxy_b_);
int save_galaxy_tree_compare(const void *a_, const void *b_);
int save_galaxy_tree_mp_comp(const void *mp_tree_data_a_, const void *mp_tree_data_b_);

#ifdef OUTPUT_BUFFERING
void save_galaxy_tree_init_output_buffer(void);
void save_galaxy_tree_flush_output_buffer(void);
void save_galaxy_tree_show_output_buffer_statistics(void);
#endif /* defined OUTPUT_BUFFERING */

#endif /* defined GALAXYTREE */

#ifdef LIGHTCONE_OUTPUT
/* save_lightcone.c */
void init_lightcone(void);
void show_lightcone_statistics(void);

void create_lightcone_galaxy_files(int file_number_);
void close_lightcone_galaxy_files(void);
void save_lightcone_galaxy_append(int galaxy_number_, int output_number_);
void save_lightcone_galaxy_finalize(int file_number_, int tree_number_);
void sort_lightcone_galaxy_in_file(FILE * lightcone_galaxy_file_);

void adjust_galaxy_for_lightcone(struct GALAXY_OUTPUT *galaxy_, const float shift_[3], const int shift_index_[3], const int output_number_);

int lightcone_galaxy_compare(const void *lightcone_galaxy_a_, const void *lightcone_galaxy_b_);
int save_lightcone_galaxy_tree_compare(const void *galaxy_tree_data_a_, const void *galaxy_tree_data_b_);

#ifdef OUTPUT_BUFFERING
void save_lightcone_galaxy_init_output_buffer(void);
void save_lightcone_galaxy_flush_output_buffer(void);
void save_lightcone_galaxy_show_output_buffer_statistics(void);
#endif /* defined OUTPUT_BUFFERING */

#endif /* defined LIGHTCONE_OUTPUT */

/* scale_cosmology.c */
void read_scaling_parameters(void);

void scale_cosmology(const int n_halos_);

#ifdef ALLOW_UNSCALE_COSMOLOGY
void un_scale_cosmology(const int n_halos_);
#endif /* defined ALLOW_UNSCALE_COSMOLOGY */

double c_correction(const float halo_mass_, const int snapshot_number_);
double dgrowth_factor_dt(const double a_, const double omega_m_, const double omega_l_);
double scale_v_cen(const int snapshot_number_);

#ifdef STAR_FORMATION_HISTORY
/* star_formation_history.c */
void create_sfh_bins();
void sfh_initialise(const int galaxy_number_);
void sfh_update_bins(const int galaxy_number_, const int snapshot_number_, const int step_number_, const double time_);
void sfh_merge(const int to_galaxy_number_, const int from_galaxy_number_);
void sfh_print(const int galaxy_number_);
#endif /* defined STAR_FORMATION_HISTORY */

#ifdef UPDATETYPETWO
/* update_type_two.c */
void update_type_two_coordinate_and_velocity(const int tree_number_, const int galaxy_number_, const int central_galaxy_number_);
void get_coordinates(float *pos_, float *vel_, const long long ID_, const int tree_number_, const int halo_number_, const int snapshot_number_);
void load_all_auxdata(const int file_number_);
#endif

#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef INDIVIDUAL_ELEMENTS
/* yields_calc_SNe_rates.c */
void SNe_rates();
int find_initial_metallicity_comp2(int Zi, int sfh_bin, int table_type);
int find_initial_mass2(double lifetime, int Zi_bin);

/* yields_elements.c */
struct elements elements_init();
struct elements elements_add(const struct elements ele1, const struct elements ele2);
struct elements elements_fraction(const struct elements ele2, const float fraction);
struct elements elements_add_fraction(const struct elements ele1, const struct elements ele2, const float fraction);

void elements_add_to(struct elements *ele, const struct elements ele2);
void elements_add_fraction_to(struct elements *ele, const struct elements ele2, const float fraction);
void elements_deduct_from(struct elements *ele, const struct elements ele2);
void elements_deduct_fraction_from(struct elements *ele, const struct elements ele2, const float fraction);

void elements_print(char s[],struct elements ele);
double elements_total(struct elements ele);
double metal_elements_total(struct elements ele);
#endif /* defined INDIVIDUAL_ELEMENTS */

/* yields_integrals.c */
void init_integrated_yields();
void integrate_yields();
int find_initial_metallicity_comp(int Zi, int sfh_bin, int table_type);
int find_initial_mass(const double lifetime, const int Zi_bin);
int max_Mi_lower(int Mi_lower, int channel_type);
int min_Mi_upper(const int Mi_upper, const int channel_type);
int find_SNII_mass_bin(const double masslimit);
int find_agb_mass_bin(const double masslimit);

#ifdef DTD
double DTDcalc (const double timevalue);
#endif /* defined DTD */

#ifdef INDIVIDUAL_ELEMENTS
void find_actual_ejecta_limits(const int channel_type, const double Mi_lower_actual, const double Mi_upper_actual, const int Mi_lower, const int Mi_upper, const int Zi,
                               double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual,
                               double* Yields_lower_actual, double* Yields_upper_actual);
#else  /* not defined INDIVIDUAL_ELEMENTS */
void find_actual_ejecta_limits(const int channel_type, const double Mi_lower_actual, const double Mi_upper_actual, const int Mi_lower, const int Mi_upper, const int Zi,
                               double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual);
#endif /* not defined INDIVIDUAL_ELEMENTS */

/* yields_read_tables.c */
void read_yield_tables(void);
double Chabrier_IMF(const double M);

#endif /* defined DETAILED_METALS_AND_MASS_RETURN */

#endif /* header guard */
