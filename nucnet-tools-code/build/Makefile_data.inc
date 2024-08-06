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

ifndef DATA_DEF

SF_DATA_URL=http://nucnet-tools.sourceforge.net/data_pub

#===============================================================================
# Get data.
#===============================================================================

ifndef DATA_DIR
DATA_DIR=$(NUCNET_TARGET)/data_pub
endif

WN_ZONE_DIR=2018-12-01
WN_NET_DIR=2018-12-01

WN_ZONE_DATA_URL=$(SF_DATA_URL)/$(WN_ZONE_DIR)/common
WN_ZONE_DATA=zone_data

ifdef WN_DATA_VERSION
  WN_NET_DATA_URL=$(SF_DATA_URL)/$(WN_NET_DIR)/$(WN_DATA_VERSION)
else
  WN_NET_DATA_URL=$(SF_DATA_URL)/$(WN_NET_DIR)/default
endif

WN_NET_PRE=my_net

ifdef WN_DATA_LEAVE_VERSION_NAME
  ifdef WN_DATA_VERSION
     WN_NET_DEST_PRE=${WN_NET_PRE}_${WN_DATA_VERSION}
  else
     WN_NET_DEST_PRE=${WN_NET_PRE}_default
  endif
else
  WN_NET_DEST_PRE=${WN_NET_PRE}
endif

ifdef WN_DATA_FULL
  WN_NET_DATA=${WN_NET_PRE}_full.xml
  WN_NET_DEST=${WN_NET_DEST_PRE}_full.xml
else
  WN_NET_DATA=${WN_NET_PRE}.xml
  WN_NET_DEST=${WN_NET_DEST_PRE}.xml
endif

data:
	mkdir -p $(DATA_DIR)
	$(WGET) $(WN_ZONE_DATA_URL)/$(WN_ZONE_DATA).tar.gz
	tar xz -C ${DATA_DIR} -f ${WN_ZONE_DATA}.tar.gz
	rm ${WN_ZONE_DATA}.tar.gz
	$(WGET) $(WN_NET_DATA_URL)/$(WN_NET_DATA).gz
	gunzip ${WN_NET_DATA}.gz
	mv ${WN_NET_DATA} ${DATA_DIR}/${WN_NET_DEST}
	rm -f ${WN_NET_DATA}*

#===============================================================================
# Get multi-zone data.
#===============================================================================

MULTI_ZONE_URL = $(SF_DATA_URL)/2018-06-14
MULTI_ZONE_DATA = multi_zone_data
MULTI_ZONE_DIR = $(NUCNET_TARGET)/data_pub

multi_zone_data:
	$(WGET) ${MULTI_ZONE_URL}/${MULTI_ZONE_DATA}.tar.gz
	mkdir -p ${MULTI_ZONE_DIR}
	tar xz -C ${MULTI_ZONE_DIR}/ -f ${MULTI_ZONE_DATA}.tar.gz
	rm ${MULTI_ZONE_DATA}.tar.gz

#===============================================================================
# Get neutrino data.
#===============================================================================

NEUTRINO_URL = $(SF_DATA_URL)/2013-06-20
NEUTRINO_DATA = neutrino
NEUTRINO_DIR = $(DATA_DIR)

neutrino_data:
	$(WGET) ${NEUTRINO_URL}/${NEUTRINO_DATA}.tar.gz
	mkdir -p ${NEUTRINO_DIR}
	tar xz -C ${NEUTRINO_DIR} -f ${NEUTRINO_DATA}.tar.gz
	rm ${NEUTRINO_DATA}.tar.gz
 
#===============================================================================
# Get weak data.
#===============================================================================

WEAK_URL = $(SF_DATA_URL)/2016-08-15
WEAK_DATA = 2d_weak_data
WEAK_DIR = $(DATA_DIR)

weak_data:
	$(WGET) ${WEAK_URL}/${WEAK_DATA}.tar.gz
	mkdir -p ${WEAK_DIR}
	tar xz -C ${WEAK_DIR} -f ${WEAK_DATA}.tar.gz
	rm ${WEAK_DATA}.tar.gz

#===============================================================================
# End include.
#===============================================================================

DATA_DEF = yes

endif
