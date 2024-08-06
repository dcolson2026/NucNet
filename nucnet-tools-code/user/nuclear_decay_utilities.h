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
//! \file
//! \brief A header file to define useful utilities for nuclear decay.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NUCLEAR_DECAY_UTILITIES_H
#define NUCLEAR_DECAY_UTILITIES_H

#include <map>
#include <omp.h>

#include <Libnucnet.h>
#include <WnSparseSolve.h>
#include <boost/lexical_cast.hpp>

#include <gsl/gsl_blas.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "network_limiter.h"

namespace user
{

//##############################################################################
// Defines.
//##############################################################################

#define S_DECAY_XPATH "[count(reactant) = 1]"
#define S_RADIOACTIVE_MASS_FRACTION   "radioactive mass"

//##############################################################################
// Prototypes.
//##############################################################################

std::pair< WnMatrix *, WnMatrix * >
get_decay_matrices( Libnucnet__Net *, double );
 
void
convert_to_mass_fractions(
  Libnucnet__Nuc *,
  gsl_vector *
);

void
convert_to_abundances(
  Libnucnet__Nuc *,
  gsl_vector *
);

void
decay_abundances(
  Libnucnet *,
  double,
  const char *,
  std::string
);

void
decay_abundances(
  Libnucnet *,
  double,
  const char *
);

void
decay_abundances(
  Libnucnet *,
  double
);

void
push_abundances_to_daughters(
  Libnucnet *,
  double
);

} // namespace user

#endif // NUCLEAR_DECAY_UTILITIES_H
