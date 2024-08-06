#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2012-2014 Clemson University.
# 
#  This file was originally written by Bradley S. Meyer and Michael J. Bojazi.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#/////////////////////////////////////////////////////////////////////////////*/

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate codes for the simple snII project.
#//!
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# Here are lines to be edited, if desired.
#///////////////////////////////////////////////////////////////////////////////

SVNURL = http://svn.code.sf.net/p/nucnet-tools/code/trunk

NUCNET_TARGET = ../nucnet-tools-code
VENDORDIR = $(NUCNET_TARGET)/vendor
OBJDIR = $(NUCNET_TARGET)/obj

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
MY_USER_DIR = $(NUCNET_TARGET)/my_user
BUILD_DIR = $(NUCNET_TARGET)/build

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Svn.
#===============================================================================

SVN_CHECKOUT := $(shell if [ ! -d $(NUCNET_TARGET) ]; then svn co $(SVNURL) $(NUCNET_TARGET); else svn update $(NUCNET_TARGET); fi )

#===============================================================================
# My user.
#===============================================================================

MY_USER := $(shell mkdir -p $(NUCNET_TARGET)/my_user)

MY_USER_COPY := $(shell if [ ! -f $(NUCNET_TARGET)/my_user/my_rates.cpp ]; then cp my_user_default/Makefile.inc $(NUCNET_TARGET)/my_user/Makefile.inc; cp my_user_default/my_rates.cpp $(NUCNET_TARGET)/my_user/my_rates.cpp; cp my_user_default/my_rates.h $(NUCNET_TARGET)/my_user/my_rates.h; fi)

#===============================================================================
# Includes.
#===============================================================================

include Makefile_data.inc

include $(BUILD_DIR)/Makefile
include $(BUILD_DIR)/Makefile.sparse
include $(MY_USER_DIR)/Makefile.inc
include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR):$(MY_USER_DIR)

#===============================================================================
# Objects.
#===============================================================================

UTIL_OBJ = $(OBJDIR)/simple_snII_utilities.o

$(UTIL_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

SNII_OBJ = $(WN_OBJ)        \
           $(MY_USER_OBJ)   \
           $(NNT_OBJ)       \
           $(SOLVE_OBJ)     \
           $(USER_OBJ)      \
           $(UTIL_OBJ)      \

ifeq ($(NNT_USE_SPARSE_SOLVER),yes)
	SNII_OBJ += $(SP_OBJ)
        CFLAGS += -DSPARSE_SOLVER
	SNII_DEP= sparse
	FLIBS= -L$(SPARSKITDIR) -lskit -lstdc++
	MC = $(FF)
else
	MC = $(CC)
	SNII_DEP=
endif

SNII_DEP += $(SNII_OBJ)

#===============================================================================
# Executables.
#===============================================================================

SNII_EXEC = simple_snII			\
        run_single_zone                 \
        run_multiple_zone_omp           \
        t9_rho_vs_mr                    \
	zone_t9_rho			\

$(SNII_EXEC): $(SNII_DEP)
	$(MC) $(SNII_OBJ) -o $(BINDIR)/$@ $@.cpp $(CLIBS) $(FLIBS)

.PHONY all_simple_snII: $(SNII_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_snII cleanall_snII

clean_snII:
	rm -f $(SNII_OBJ)

cleanall_snII: clean_snII
	rm -f $(BINDIR)/$(SNII_EXEC) $(BINDIR)/$(SNII_EXEC).exe
