# Changes

The following non-exhaustive list describes the changes made to the 2015 public-release version of L-Galaxies to create the current version:

## 2019-03-26 (and before):

### general:

- The directory layout of the code has been adjusted:
  - Object files are placed into their own ./obj/ subdirectory.
  - Executables are placed into their own ./bin/ subdirectory.
  - Dependency files (see below) are placed into their own ./dep/ subdirectory.
  
- The MakeFile and related files have been updated:
  - The compile-time options are now all in one MakeFile_options file that should work for both standard and MCMC mode.
  - Source code dependencies are now deduced automatically by make (and gcc's C-preprocessor) and stored in ./dep/*.dep files created in that process (at least if GNU make and gcc are used).
  
### naming/style/coding conventions:

- Functions, e.g., for the galaxy model have been regrouped into several source files.

- Certain function names have been adjusted to better reflect their function.

- Function declarations (aka 'prototypes') in proto.h are now sorted by source file name and order in source file.

- Local variable names (local to a function) can now be usually recognized by ending with an underscore. Conversely, a variable name without an underscore at the end likely denotes a global variable (note: this change has not been applied to all parts of the code yet).

- Variable names have been made more descriptive, e.g. most of the i, j, k, etc. have been changed to, e.g., galaxy_number_, halo_number_, sfh_bin_number_, etc. 

- Local variables now are sometimes not declared at the beginning of a block (usually a function) any more, but where they can or must be set in the code. This improves locality and readability, but violates pre-C99 ANSI/ISO C standard.

- Local variables (including most function arguments) are declared constant where possible. 

- Boolean variables are now declared with type bool, and the 'true' and 'false' are used instead of 1 or 0 for boolean values. This improves readability, but violates pre-C99 ANSI/ISO C standard.

- Tabs have been replaced by spaces. The indention scheme, bracket placing, spacing around operators, etc. have been updated to improve readability.

- #elif and #endif now always come with info on the corresponding #if to improve readability for the many nested #ifdef's.

- Code parts causing warnings about unused variables, uninitialized use, set-but-not-used variables, etc. have been adjusted (e.g. by elimination unused variables) to not cause such warnings. Now turning on such warning for compilation should help catching bugs in the code.

### regressions, optimizations, and new features:

- The mutual compatibility of compile-time options is now checked at compile-time.

- The handling of runtime options (parameters read from the input parameter file) has been updated:
  - Functionality for optional parameters has been added.
  - Certain checks on the parameters are performed (e.g. range chacks) after reading them from file. For some parameters that currently only have a single valid value, the checks and comparisons in the galaxy physics parts of the code have been replaced by a check after reading the parameter file. If one wants to add another possible value to such a parameter, those checks and comparisons in the galaxy physics code have to be put back again.

- Treatment of physical constants and units has been updated:
  - All physical constants and most units are now defined as macro/static constants in a single header file physical_constants_and_units.h.
  - Functionality to use updated physical constants has been added.
  
- Random numbers are now gerenated using GSL.

- The cosmological distance and time calculations have been updated:
  - There is now a separate source file cosmology.c for handling cosmological functions.
  - Functionality for computing various distances and redshifts has been added/extended.
  - For flat LCDM, distances are now computed using hypergeometric functions from GSL.
  - For non-flat LCDM, distances are computed by numerical intergation with a new kernel, which speeds up integration for given accuracy.
  - There are several other improvements on convenience/speed of functions for cosmic times or magnitudes.
  
- Functionality to output galaxies on a lightcone has been added:
  - Galaxy lightcones (with simple geometries defined by ra, dec, redshift) can now be output in addition or instead of standard output.
  - As a side effect, the treatment of magnitudes for forward and backward interpolation has been changed. For example, compile-time options KITZBICHLER and MOMAF have been replaced by one option OUTPUT_FB_OBS_MAG.
  
- Functionality to sort the galaxy output has been added/improved.

- Functionality for internally buffering the galaxy output has been added. (This option is likely not needed for most current computer systems, which are good at buffering output on their own.)
  
- The algorithm for turning halo trees into galaxies have been updated:
  - Once upon a time, the code probably generated galaxies by recursively calling the function 'construct_galaxies' to generate the galaxy progenitors (wich called construct_galaxies to create their progenitors, and so on), and then assembling the galaxies. In the current code, galaxies are explicitly created in chronological order, but construct_galaxies still has several checks and branches left so it could be called in a recursive manner. There is now a compile-time option to employ a simplified version of construct_galaxies without some of these checks (which has the side effect of sometimes changing the order of galaxies in the output).
  
- Saving galaxies for MCMC has been updated:
  - Galaxies to be used for statistical analysis in MCMC mode are now identified based on their host halo in a more efficient way, improving speed on this part of the code.
  - As a side effect, halo numbers for halos to be used in MCMC analysis are sorted upon reading them from file.
  
- The treatment of observational errors on, e.g., stellar masses in MCMC has been updated to improve speed.

- The code has been further optimized for speed by various methods:
  - Certain functions have been inlined or declared inline and put near there use to help compilers that don't optimize across units.
  - String comparisons to switch between a finite number of known alternatives in speed-critical functions (e.g. transfer_gas and tranfer_stars) have been replaced by comparisons on enums.
  - The order of certain nested loops has been changed to both explicitly reduce the amount of repeated computations, and to help optimizing compilers to move computations out of inner loops.
  - The order of dimensions of certain multi-dimensional arrays has been adjusted to improve caching when accessed within nested loops.
  - Certain functions employing numerical root finding have been replaced by standardized root-finding algorithms or by calls to special functions.
  - Certain functions employing interpolation have been replaced by standardized interpolation algorithms.
  - Certain physics calculations (that often employed multiple intermediate computations) have been simplified into mathematical equivalent, but faster calculations. In most cases, this produces exactly the same output as before, but there are few instances where outputs between the old and new version of the code differ slightly due to finite numerical accuracy.
  
  
### bug fixes:

- bugfix:
  - file: various files
  - function: various functions
  - description: the sizes of various char[] (those used for, e.g., file names or error messages) have been adjusted to avoid potential write past end when used in sprintf.

- bugfix:
  - file: model_dust_extinction_inline.c
  - function: get_extinction
  - description: corrected off-by-one error on table index for interpolation of dust extiction.

- bugfix:
  - file: save_galtree.c
  - function: save_galaxy_tree_append
  - description: sfh_numbins was set to sfh_ibin. sfh_numbins is now set to sfh_ibin + 1, and this is now done in prepare_galaxy_for_output().

- bugfix:
  - file: save.c
  - function: prepare_galaxy_for_output.c
  - description: sfh_numbins was not set. it only was set in save_galaxy_tree_append(), where it was set to sfh_ibin. The correct value, however, should be (sfh_ibin + 1), and it should be set here to work for all kinds of output (snapshot, tree, lightcone).
  
- bugfix:
  - file: save.c
  - function: prepare_galaxy_for_output.c
  - description: #ifndef GALAXYTREE, DisruptOn was not set, but left to garbage. now it will be set to zero in that case.

- bugfix:
  - file: post_process_spec_mags.c
  - function: post_process_spec_mags
  - description: 
    before: filter number for r band was hard-coded as 17.
    now   : R_BAND_FILTER_NUMBER defines r band filter number.

- bugfix:
  - file: post_process_spec_mags.c
  - function: post_process_spec_mags
  - description:
    before: if not defined OUTPUT_REST_MAGS, but 0 <= (galaxy_->DiskMass + galaxy_->BulgeMass), galaxy_->rbandWeightAge would compute to 0./0. (i.e. NaN).
    now   : if not defined OUTPUT_REST_MAGS, galaxy_->rbandWeightAg is set to 0 and stays 0.

- bugfix:
  - file: model_spectro_photometric.c
  - function: add_to_luminosities
  - description: for MetallicityOption == 0 (means only solar metallicity) only the table index for interpolation was set, but not the fractions for interpolation. fractions are now computed correctly, too.

- bugfix:
  - file: model_misc.c
  - function: transfer_stars
  - description: it was tested that "BurstMass" could only be transferred to "BurstMass", but it was "Burst" to indicate BurstMass transfer and "BurstMass" was never mentioned elsewhere. now "BurstMass" is replaced by "Burst" (and then by BurstComponent).

- bugfix:
  - file: model_cooling.c
  - function: do_AGN_heating
  - description: a variable 'dist' was used in a comparison to decide on the details of AGN heating for satellites, but never set. now it is set to the distance between galaxy and FOF central galaxy
 
 
### possible bugs (and fixes):

- bugfix(?):
  - file: sam.c
  - function: evolve_galaxies
  - description: moved "int start = NGalTree;" back to before loop on finding central and satellite galaxies, thus fixing calculation of "GalTree[p].FOFCentralGal".

- bugfix (?):
  - file: model_misc.c
  - function: transfer_stars
  - description: BurstMass was given to receiving galaxy, but not taken from donor galaxy. code used to be: Gal[q_].BurstMass -= 0.; now replaced by: Gal[q_].BurstMass -= Mass;

- bugfix (?):
  - file: model_mergers.c
  - function: get_merger_center
  - description: the branch returning 0 used to depend on i_ == 0. the check 'if(i_ == 0)' has now been removed there.
 
- bug (?):
  - file: model_misc.c
  - function: get_initial_disk_radius
  - description: code checks for zero gas spin, but perhaps should check for zero halo spin instead
 
- bug (?):
  - file: proto.h
  - function: transfer_ICL
  - description: missing function definition (implementation)

 
