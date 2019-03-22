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

 
/** @file mcmc_proto.h
 * 
 * @author Bruno Henriques
 * @author Stefan Hilbert
 * @date   2008, 2018
 */
 
/* mcmc.c */
void Senna (void);

void print_parameters (const bool current_MCMC_step_is_accepted_, FILE *fmcmc);
void initialize_mcmc_par_and_lhood (FILE *fmcmc);
void read_mcmc_par (const int snapshot_number_);
void read_sample_info(void);
void read_observations(void);
void open_files_with_comparison_to_observations(void);
void close_files_with_comparison_to_observations(void);
double propose_new_parameters(void);
void mark_halos_in_MCMC_sample(const int tree_number_);
void free_MCMC_FOF(void);
int MCMC_FOF_compare_FoFID(const void *MCMC_FOF_a_, const void *MCMC_FOF_b_);
#ifdef MR_PLUS_MRII
void change_dark_matter_sim(const char SimName[]);
#endif

#ifdef HALOMODEL
void assign_FOF_masses(const int snapshot_number_, const int tree_number_);
#endif

/* mcmc_save.c */
// void save_galaxy_for_mcmc(const int gal_index_, const int output_is_expected_);
void save_galaxy_for_mcmc(const int gal_index_);

/* mcmc_likelihood.c */
double get_likelihood (void);

void bin_function(const int output_number_, const int observation_number_, const double *sam_data_, double *binned_sam_data_);
void bin_red_fraction(const int output_number_, const int observation_number_, double *binned_red_fraction_);
void bin_passive_fraction(const int output_number_, const int observation_number_, double *binned_passive_fraction_);
void bin_bulge_fraction(const int output_number_, const int observation_number_, double *binned_bulge_fraction_);
void bin_ColdGasFractionvsStellarMass(const int output_number_, const int observation_number_, double *binned_gas_fraction_);
void bin_color_hist(const int output_number_, const int observation_number_, double *binned_color_hist_);
void bin_bhbm(const int output_number_, double *binned_blackhole_up_, double *binned_blackhole_down_);

void correct_for_correlation(const int output_number_);

void compute_correlation_func(const int output_number_, const int observation_number_, const float mingalmass, const float maxgalmass, double *binsamdata);

double chi_square_probability(const int output_number_, const int observation_number_, const double *sam_data_);
double maximum_likelihood_probability(const int output_number, const int observation_number_, const double *sam_fract_);
double binomial_probability(const int output_number_, const int observation_number_, const double *sam_up_, const double *sam_down_);

double gammp (const double a, const double x);
double gammq (const double a, const double x);
double gser (const double a, const double x);
double gcf(const double a, const double x);
double gammpapprox(const double a, const double x, const int psig);
double gammln(const double xx);
double betai(const double a, const double b, const double x);
double betacf(const double a, const double  b, const double x);
