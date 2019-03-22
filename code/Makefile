#  Copyright (C) <2016>  <L-Galaxies>

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>


.SUFFIXES:
.SUFFIXES: .o .c .h .dep
# suffixes of the file names to be used in the makefile

VPATH = .
# paths to be looked into by make while processing the makefile

##################################################################
# setting path variables for sources,deps,etc.:                  #
##################################################################

SRCDIR = ./code
OBJDIR = ./obj
BINDIR = ./bin
DEPDIR = ./dep

##################################################################
# the executable's name:                                         #
##################################################################

EXEC   = L-Galaxies

EXEC  := $(BINDIR)/$(EXEC)

##################################################################
# the source and object files:                                   #
##################################################################

SRCS   = main.c sam.c io_tree.c init.c cool_func.c \
         save.c save_galtree.c save_lightcone.c \
         cosmology.c \
         mymalloc.c read_parameters.c \
         peano.c allvars.c update_type_two.c \
         model_infall.c \
         model_cooling.c \
         model_starformation_and_feedback.c \
         model_reincorporation.c \
         model_mergers.c \
         model_dust.c \
         model_misc.c \
         model_disrupt.c \
         model_stripping.c \
         scale_cosmology.c

OBJS = $(SRCS:%.c=$(OBJDIR)/%.o)

##################################################################
# files for all (not found by dep files below):                 #
##################################################################

DEPSFORALL = Makefile My_Makefile_options Makefile_compilers 

# ./code/allvars.h  ./code/proto.h  ./code/aux_functions_inline.h  ./code/metals_inline.h

##################################################################
# include files for compile time options:                        #
##################################################################

# Either include the default set of Makefile options, or define your own
include My_Makefile_options
#include My_Makefile_options_MCMC
#include My_Makefile_options_MCMC_HaloModel

##################################################################
# include files for compiler choice:                             #
##################################################################

# Choose your system type (needs to match an entry in Makefile_compilers)
SYSTYPE = "MyMachine"
include Makefile_compilers

##################################################################
# now for real:                                                  #
##################################################################

LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas 

CFLAGS =   -g $(OPTIONS) $(OPT) -DCOMPILETIMESETTINGS=\""$(OPT)"\" $(OPTIMIZE) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

# the rule for making the object files %.o from the source files %.cpp:
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPSFORALL)
	$(CC) -c -o $@ $(CFLAGS) $<

# files containing the dependencies of the object files:
DEPS = $(SRCS:%.c=$(DEPDIR)/%.dep)

$(DEPDIR)/%.dep: $(SRCDIR)/%.c
	@set -e; rm -f $@; \
         $(CC) -MM $(CFLAGS) $< > $@.$$$$; \
         sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' < $@.$$$$ > $@;\
         rm -f $@.$$$$

include $(DEPS)

##################################################################
# phony targets:
##################################################################
.PHONY: clean tidy metadata metadata_db

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	${CC_MD} ${OPT} ${CFLAGS} -E ./code/allvars.h -o ./code/allvars.i
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_TypeString.awk > ./AuxCode/awk/L-Galaxies_Types.txt
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_DDL.awk > ./AuxCode/awk/L-Galaxies_DDL.sql	
ifeq (NORMALIZEDDB,$(findstring NORMALIZEDDB,$(OPT)))
	awk -f ./AuxCode/awk/extractSFH_BIN.awk ./code/allvars.i |awk -f ./AuxCode/awk/SFH_BIN_2_DDL.awk >> ./AuxCode/awk/L-Galaxies_DDL.sql
else
	awk -f ./AuxCode/awk/extractSFH_Time.awk ./code/allvars.i |awk -f ./AuxCode/awk/SFH_Time_2_DDL.awk >> ./AuxCode/awk/L-Galaxies_DDL.sql
endif	
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_struct.awk >  ./AuxCode/awk/idl/LGalaxy.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_hists.awk > ./AuxCode/awk/idl/LGalaxy_plot.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_testfloats.awk > ./AuxCode/awk/idl/LGalaxy_testfloats.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_zerofloats.awk > ./AuxCode/awk/idl/LGalaxy_zerofloats.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_python_struct.awk > ./AuxCode/awk/LGalaxy.py
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_LGalaxy.awk > ./AuxCode/awk/L-Galaxies.h
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_FileFormat.awk > ./AuxCode/awk/L-Galaxies_FileFormat.csv
	awk -f ./AuxCode/awk/extractSFH_BIN.awk ./code/allvars.i |awk -f ./AuxCode/awk/MOMAF_INPUT_2_MoMaFGalaxy.awk >> ./AuxCode/awk/L-Galaxies.h

metadata_db:
	awk -f ./AuxCode/awk/extract_struct_metals.awk ./code/allvars.i > ./AuxCode/awk/structs.dat
	awk -f ./AuxCode/awk/extract_struct_elements.awk ./code/allvars.i >> ./AuxCode/awk/structs.dat
	awk -f ./AuxCode/awk/extract_struct_GALAXY_OUTPUT.awk ./code/allvars.i >> ./AuxCode/awk/structs.dat
