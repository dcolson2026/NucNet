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
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to download data.  It is included from the main
#//!        makefile.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef SIMPLE_SNII_DATA_DEF

SIMPLE_SNII_DATA_URL=http://nucnet-projects.sourceforge.net/data_pub

SIMPLE_SNII_DATA_DIR=$(SIMPLE_SNII_DATA_URL)/2013-10-23

#===============================================================================
# Get data.
#===============================================================================

ifndef DATA_DIR
DATA_DIR=$(NUCNET_TARGET)/data_pub
endif

SIMPLE_SNII_DATA=simple_snII__zone.xml
N15_STUDY_DATA=n15_study.tar
ZONE_NU=zone_nu.xml
ZONE_NU_INCLUDE=zone_nu_include.xml

simple_snII_data:
	mkdir -p $(DATA_DIR)
	wget $(SIMPLE_SNII_DATA_DIR)/$(SIMPLE_SNII_DATA)
	mv $(SIMPLE_SNII_DATA) $(DATA_DIR)

n15_study_data:
	mkdir -p $(DATA_DIR)
	wget $(SIMPLE_SNII_DATA_DIR)/$(N15_STUDY_DATA).gz
	mv ${N15_STUDY_DATA}.gz $(DATA_DIR)
	tar xz -C ${DATA_DIR} -f ${DATA_DIR}/${N15_STUDY_DATA}.gz
	rm ${DATA_DIR}/${N15_STUDY_DATA}.gz

simple_snII_neutrino_data:
	make neutrino_data
	rm ${DATA_DIR}/${ZONE_NU}
	wget $(SIMPLE_SNII_DATA_DIR)/${ZONE_NU_INCLUDE}
	mv ${ZONE_NU_INCLUDE} ${DATA_DIR}

#===============================================================================
# End include.
#===============================================================================

SIMPLE_SNII_DATA_DEF = yes

endif
