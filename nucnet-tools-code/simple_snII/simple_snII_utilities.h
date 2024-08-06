////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2012-2013 Clemson University.
//
// This file was originally written by Bradley S. Meyer and Michael J. Bojazi.
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
//! \brief A header file to define useful utilities for the simple snII
//!   project.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_SNII_UTILITIES_H
#define SIMPLE_SNII_UTILITIES_H

#include <boost/tuple/tuple.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/network_utilities.h"
#include "user/thermo.h"

#define S_P0    "p0"
#define S_P     "p"
#define S_E0    "e0"
#define S_E     "e"
#define S_U     "u"
#define S_F1    "f1"
#define S_F2    "f2"
#define S_F3    "f3"
#define S_V     "v"

#define S_SHOCK_ENERGY   "shock energy"
#define S_SHOCK_SPEED    "shock speed"
#define S_PARTICLE_TYPE  "particle type"
#define S_SHOCK_TIME     "shock time"
#define S_SHOCK_RHO      "shock rho"
#define S_SHOCK_T9       "shock t9"
#define S_VELOCITY_FACTOR "velocity factor"
#define S_XML_OUTPUT_FILE  "output file"

#define S_LOW_T_NETWORK_T9 "low T network T9"
#define S_LOW_T_NETWORK_REAC_XPATH "low T network reaction xpath"

int
zone_compare_by_mass( const Libnucnet__Zone *, const Libnucnet__Zone * );

boost::tuple<std::vector<double>, std::vector<double>, std::vector<double> >
get_zone_vectors(
  nnt::Zone&,
  Libnucnet *
);

double T9_post( double, nnt::Zone& );

double rho_post( double, nnt::Zone& );

void
update_zone_properties(
  nnt::Zone&,
  Libnucnet *
);

std::pair<Libnucnet *, std::vector<nnt::Zone> >
get_evolution_network_and_zones( Libnucnet *, char ** );

int
copy_properties_to_zone(
  const char *,
  const char *,
  const char *,
  const char *,
  Libnucnet__Zone *
);

void
set_low_T_network( nnt::Zone& );

#endif // SIMPLE_SNII_UTILITIES_H
