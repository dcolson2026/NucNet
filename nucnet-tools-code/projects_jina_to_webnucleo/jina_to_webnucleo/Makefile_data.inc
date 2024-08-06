#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2014 Clemson University.
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

ifndef JINA_TO_WEBNUCLEO_DATA_DEF

JINA_TO_WEBNUCLEO_DATA_URL=http://nucnet-projects.sourceforge.net/data_pub

JINA_TO_WEBNUCLEO_DATA_DIR=$(JINA_TO_WEBNUCLEO_DATA_URL)/2014-12-22

#===============================================================================
# Get data.
#===============================================================================

ifndef DATA_DIR
DATA_DIR=$(NUCNET_TARGET)/data_pub
endif

JINA_TO_WEBNUCLEO_DATA=jina_to_webnucleo.tar

jina_to_webnucleo_data:
	mkdir -p $(DATA_DIR)
	wget $(JINA_TO_WEBNUCLEO_DATA_DIR)/$(JINA_TO_WEBNUCLEO_DATA).gz
	mv ${JINA_TO_WEBNUCLEO_DATA}.gz $(DATA_DIR)
	tar xz -C ${DATA_DIR} -f ${DATA_DIR}/${JINA_TO_WEBNUCLEO_DATA}.gz
	rm ${DATA_DIR}/${JINA_TO_WEBNUCLEO_DATA}.gz

#===============================================================================
# End include.
#===============================================================================

JINA_TO_WEBNUCLEO_DATA_DEF = yes

endif
