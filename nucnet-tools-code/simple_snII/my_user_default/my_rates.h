////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2012-2013 Clemson University.
//
// This file was originally written by Michael J. Bojazi.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file simple_snII_utilities.h
//! \brief A header file for user-defined rates for the simple_snII
//!   project.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef MY_RATES_H
#define MY_RATES_H

#include <iostream>

#include <Libnucnet.h>

namespace my_user
{

double
hfcz_n15_p_a_c12_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
fcz2_n14_a_g_f18_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
fcz2_n15_p_g_o16_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
dh03_o18_a_g_ne22_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
ra94_o18_n_g_o19_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
wk82_f18_p_g_ne19_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
wk82_f18_p_a_o15_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
la90_o17_p_g_f18_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

void
register_rate_functions( Libnucnet__Reac * );

} // namespace my_user

#endif  // MY_RATES_H
