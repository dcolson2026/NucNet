#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2013 Clemson University.
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
#//! \brief A makefile to generate codes for jina->webnucleo xml conversions.
#//!
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# Here are lines to be edited, if desired.
#///////////////////////////////////////////////////////////////////////////////

SVNURL = svn://svn.code.sf.net/p/nucnet-tools/code/trunk

NUCNET_TARGET = ../nucnet-tools-code
VENDORDIR = $(NUCNET_TARGET)/vendor
OBJDIR = $(NUCNET_TARGET)/obj

NNT_DIR = $(NUCNET_TARGET)/nnt
BUILD_DIR = $(NUCNET_TARGET)/build

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Svn.
#===============================================================================

SVN_CHECKOUT := $(shell if [ ! -d $(NUCNET_TARGET) ]; then svn co $(SVNURL) $(NUCNET_TARGET); else svn update $(NUCNET_TARGET); fi )

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Includes.
#===============================================================================

include Makefile_data.inc

include $(BUILD_DIR)/Makefile

#===============================================================================
# Paths and objects.
#===============================================================================

VPATH = $(BUILD_DIR):$(NNT_DIR)

JINA_XML_OBJS = $(WN_OBJ)           	\
                $(NNT_OBJ)

#===============================================================================
# Executables.
#===============================================================================

JINA_XML_EXEC = nuclide_converter       \
                reaction_converter      \

$(JINA_XML_EXEC): $(JINA_XML_OBJS)
	$(CC) $(JINA_XML_OBJS) -o $(BINDIR)/$@ $@.cpp $(CLIBS)

.PHONY all_jina_to_webnucleo: $(JINA_XML_EXEC)

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean_jina_xml

clean_jina_xml:
	rm -f $(OBJ)
