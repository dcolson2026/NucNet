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
#//! \brief An example makefile to generate user-defined routine objects.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef MY_USER_DEF

#===============================================================================
# User object files.  Add routines to MY_USER_OBJ.
#===============================================================================

MY_USER_OBJ = $(OBJDIR)/my_hydro.o   \
              $(OBJDIR)/my_screen.o  \

$(MY_USER_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

#===============================================================================
# Clean up.
#===============================================================================

.phony: clean_my_user cleanall_my_user

clean_my_user:
	rm -f $(MY_USER_OBJ)

cleanall_my_user:
	rm -fr $(MY_USER_DIR)

#===============================================================================
# End include.
#===============================================================================

MY_USER_DEF = yes

endif
