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
# SRCDIR = ./version_2018_02_06_code
OBJDIR = ./obj
BINDIR = ./bin
DEPDIR = ./dep

AWKDIR = ./AuxCode/awk

##################################################################
# the executable's name:                                         #
##################################################################

EXEC   = L-Galaxies

EXEC  := $(BINDIR)/$(EXEC)

##################################################################
# the source and object files:                                   #
##################################################################

SRCS   = main.c \
         allvars.c \
	 mymalloc.c \
	 io.c \
	 io_tree.c \
	 init.c \
	 read_parameters.c \
	 cosmology.c \
         peano.c \
	 update_type_two.c \
	 model.c \
         model_cooling.c \
         model_disrupt.c \
         model_dust.c \
         model_infall.c \
         model_mergers.c \
         model_misc.c \
         model_reincorporation.c \
         model_starformation_and_feedback.c \
         model_stripping.c \
         save.c \
         scale_cosmology.c

OBJS = $(SRCS:%.c=$(OBJDIR)/%.o)


##################################################################
# files for all (not found by dep files below):                 #
##################################################################

DEPSFORALL = Makefile Makefile_options Makefile_compilers 

# ./code/allvars.h  ./code/proto.h  ./code/aux_functions_inline.h  ./code/metals_inline.h

##################################################################
# include files for compile time options:                        #
##################################################################

# Either include the default set of Makefile options, or define your own
include Makefile_options

##################################################################
# include files for compiler choice:                             #
##################################################################

# Choose your system type in Makefile_compilers:
include Makefile_compilers

##################################################################
# removing excessive whitespace in options:                      #
##################################################################
# (comment if your make does not like this)
OPT := $(strip $(OPT))

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
	rm -f $(OBJDIR)/*.o $(DEPDIR)/*.dep

tidy:
	rm -f $(OBJS) $(DEPS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	${CC_MD} ${OPT} ${CFLAGS} -E $(SRCDIR)/allvars.h -o $(SRCDIR)/allvars.i
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/GALAXY_OUTPUT_2_TypeString.awk > $(AWKDIR)/L-Galaxies_Types.txt
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/GALAXY_OUTPUT_2_DDL.awk > $(AWKDIR)/L-Galaxies_DDL.sql	
	${CC_MD} ${OPT} ${CFLAGS} -E $(SRCDIR)/lightcone_galaxy_output_type.h -o $(SRCDIR)/lightcone_galaxy_output_type.i
ifeq (NORMALIZEDDB,$(findstring NORMALIZEDDB,$(OPT)))
	awk -f $(AWKDIR)/extractSFH_BIN.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/SFH_BIN_2_DDL.awk >> $(AWKDIR)/L-Galaxies_DDL.sql
else
	awk -f $(AWKDIR)/extractSFH_Time.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/SFH_Time_2_DDL.awk >> $(AWKDIR)/L-Galaxies_DDL.sql
endif	
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/idl/GALAXY_OUTPUT_2_IDL_struct.awk >  $(AWKDIR)/idl/LGalaxy.pro
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/idl/GALAXY_OUTPUT_2_IDL_hists.awk > $(AWKDIR)/idl/LGalaxy_plot.pro
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/idl/GALAXY_OUTPUT_2_IDL_testfloats.awk > $(AWKDIR)/idl/LGalaxy_testfloats.pro
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/idl/GALAXY_OUTPUT_2_IDL_zerofloats.awk > $(AWKDIR)/idl/LGalaxy_zerofloats.pro
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/GALAXY_OUTPUT_2_python_struct.awk > $(AWKDIR)/LGalaxy.py
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/GALAXY_OUTPUT_2_LGalaxy.awk > $(AWKDIR)/L-Galaxies.h
	awk -f $(AWKDIR)/extractGALAXY_OUTPUT.awk $(SRCDIR)/allvars.i |awk -f $(AWKDIR)/GALAXY_OUTPUT_2_FileFormat.awk > $(AWKDIR)/L-Galaxies_FileFormat.csv
	awk -f $(AWKDIR)/extract_lightcone_galaxy_output_type.awk $(SRCDIR)/lightcone_galaxy_output_type.i |awk -f $(AWKDIR)/lightcone_galaxy_output_type_2_lightcone_galaxy_output_type.awk > $(AWKDIR)/lightcone_galaxy_output_type.h
	awk -f $(AWKDIR)/extract_lightcone_galaxy_output_type.awk $(SRCDIR)/lightcone_galaxy_output_type.i |awk -f $(AWKDIR)/lightcone_galaxy_output_type_2_python_struct.awk > $(AWKDIR)/lightcone_galaxy_output_type.py

metadata_db:
	awk -f $(AWKDIR)/extract_struct_metals.awk $(SRCDIR)/allvars.i > $(AWKDIR)/structs.dat
	awk -f $(AWKDIR)/extract_struct_elements.awk $(SRCDIR)/allvars.i >> $(AWKDIR)/structs.dat
	awk -f $(AWKDIR)/extract_struct_GALAXY_OUTPUT.awk $(SRCDIR)/allvars.i >> $(AWKDIR)/structs.dat
