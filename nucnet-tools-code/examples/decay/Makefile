#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2017 Clemson University.
# 
#  This file was originally written by Bradley S. Meyer.
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
#//! \brief A makefile to generate codes for the multi-zone project.
#//!
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# Here are lines to be edited, if desired.
#///////////////////////////////////////////////////////////////////////////////

NUCNET_TARGET = ../..
VENDORDIR = $(NUCNET_TARGET)/vendor
OBJDIR = $(NUCNET_TARGET)/obj

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build

FC = gfortran

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Includes.
#===============================================================================

include $(BUILD_DIR)/Makefile
include $(BUILD_DIR)/Makefile.sparse
include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

#===============================================================================
# Debugging, if desired.
#===============================================================================

ifdef DEBUG
  CFLAGS += -DDEBUG
  FFLAGS += -DDEBUG
endif

#===============================================================================
# Objects.
#===============================================================================

DECAY_OBJS = $(SP_OBJ)          \
             $(ILU_OBJ)         \
             $(WN_OBJ)          \
             $(MY_USER_OBJ)     \
             $(NNT_OBJ)         \
             $(SOLVE_OBJ)       \
             $(USER_OBJ)        \
             $(DECAY_OBJ)       \

DECAY_DEP = $(DECAY_OBJS)

CFLAGS += -DSPARSKIT2 -lgfortran
DECAY_DEP += sparse
FLIBS= -L$(SPARSKITDIR) -lskit

#===============================================================================
# Executables.
#===============================================================================

DECAY_EXEC = decay_zone

$(DECAY_EXEC): $(DECAY_DEP)
	$(CC) -c -o $(OBJDIR)/$@.o $@.cpp
	$(CC) $(DECAY_OBJS) -o $(BINDIR)/$@ $(OBJDIR)/$@.o $(CLIBS) $(FLIBS)

.PHONY all_decay: $(DECAY_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_decay cleanall_decay

clean_decay:
	rm -f $(DECAY_OBJ)

cleanall_decay: clean_decay
	rm -f $(BINDIR)/$(DECAY_EXEC) $(BINDIR)/$(DECAY_EXEC).exe
