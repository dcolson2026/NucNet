////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Code for user-defined rate routines.
////////////////////////////////////////////////////////////////////////////////

#include "my_rates.h"

namespace my_user
{

//##############################################################################
// hfcz_n15_p_a_c12_rate_function().
//##############################################################################

double
hfcz_n15_p_a_c12_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    1.08e12 / pow( d_t9, 2./3. ) * exp( -15.251 / pow( d_t9, 1./3. ) - 
      pow( d_t9/0.522, 2. ) ) * 
    ( 1. + 0.027 * pow( d_t9, 1./3. ) + 2.62 * 
      pow( d_t9, 2./3. ) + 0.501 * d_t9 + 5.36 * pow( d_t9, 4./3. ) + 2.60 * 
      pow( d_t9, 5./3. ) ) + 
    1.19e8 / pow( d_t9, 3./2. ) * exp( -3.676/d_t9 ) + 5.41e8 / 
    pow( d_t9, 1./2. ) * exp(-8.926/d_t9) + 0.1 * 4.72e8 / pow( d_t9, 3./2. ) *     exp( -7.721/d_t9 ) + 2.20e9 / pow( d_t9, 3./2. ) * exp( -11.418/d_t9 );

}

//##############################################################################
// fcz2_n14_a_g_f18_rate_function().
//##############################################################################

double
fcz2_n14_a_g_f18_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    7.78e9 / pow( d_t9, 2./3. ) * exp( -36.031 / pow( d_t9, 1./3. ) - 
      pow( d_t9/0.881, 2. ) ) * 
    ( 1. + 0.012 * pow( d_t9, 1./3. ) + 1.45 * 
      pow( d_t9, 2./3. ) + 0.117 * d_t9 + 1.97 * pow( d_t9, 4./3. ) + 0.406 *
      pow( d_t9, 5./3. ) ) + 
    2.36e-10 / pow( d_t9, 3./2. ) * exp( -2.798/d_t9) + 2.03e0 / 
    pow( d_t9, 3./2. ) * exp( -5.054/d_t9 ) + 1.15e4 / pow( d_t9, 2./3. ) *
    exp( -12.310/d_t9 );

}

//##############################################################################
// fcz2_n15_p_g_o16_rate_function().
//##############################################################################

double
fcz2_n15_p_g_o16_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    9.78e8 / pow( d_t9, 2./3. ) * exp( -15.251 / pow( d_t9, 1./3. ) - 
      pow( d_t9/0.450, 2. ) ) * 
    ( 1. + 0.027 * pow( d_t9, 1./3. ) + 0.219 * 
      pow( d_t9, 2./3. ) + 0.042 * d_t9 + 6.83 * pow( d_t9, 4./3. ) + 3.32 *
      pow( d_t9, 5./3. ) ) + 
    1.11e4 / pow( d_t9, 3./2. ) * exp( -3.328/d_t9) + 1.49e4 / 
    pow( d_t9, 3./2. ) * exp( -4.665/d_t9 ) + 3.80e6 / pow( d_t9, 3./2. ) *
    exp( -11.048/d_t9 );

}

//##############################################################################
// dh03_o18_a_g_ne22_rate_function().
//##############################################################################

double
dh03_o18_a_g_ne22_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    1.95e-13 / pow( d_t9, 3./2. ) * exp( -2.069/d_t9 ) + 1.248e-2 / 
    pow( d_t9, 3./2. ) * exp( -4.462/d_t9 ) + 1.845e-2 / pow( d_t9, 3./2. ) * 
    exp( -5.373/d_t9 ) + 5.95 / pow( d_t9, 3./2. ) * exp( -6.285/d_t9 ) + 
    12.74 / pow( d_t9, 3./2. ) * exp( -7.119/d_t9 ) + 31.19 / 
    pow( d_t9, 3./2. ) * exp( -7.287/d_t9 ) + 3.22e5 / pow( d_t9, 1./2. ) * 
    exp( -21.801/d_t9 );

}

//##############################################################################
// ra94_o18_n_g_o19_rate_function().
//##############################################################################

double
ra94_o18_n_g_o19_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    2.12e1 + 5.87e3 * d_t9 + 4.08e4 / pow( d_t9, 3./2. ) * exp( -1.769/d_t9 ) +
    9.70e4 / pow( d_t9, 3./2. ) * exp( -4.305/d_t9 ) + 5.11e5 / 
    pow( d_t9, 3./2. ) * exp( -7.265/d_t9 );

}

//##############################################################################
// wk82_f18_p_g_ne19_rate_function().
//##############################################################################

double
wk82_f18_p_g_ne19_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    1.658e7 / pow( d_t9, 2./3. ) * exp( -18.06/pow( d_t9, 1./3. ) ) * 
    ( 4.604 + 0.106 * pow( d_t9, 1./3. ) + 0.053 * pow( d_t9, 2./3. ) + 0.009 *
      d_t9 - 0.036 *  pow( d_t9, 4./3. ) - 0.015 * pow( d_t9, 5./3. ) ) +
    4.55e-14 / pow( d_t9, 3./2. ) * exp( -0.302/d_t9 ) + 
    3.27e2 / pow( d_t9, 3./2. ) * exp( -3.84/d_t9 ) +
    1.32e4 / pow( d_t9, 3./2. ) * exp( -5.22/d_t9 ) +
    93. / pow( d_t9, 3./2. ) * exp( -4.29/d_t9 ); 

}

//##############################################################################
// wk82_f18_p_a_o15_rate_function().
//##############################################################################

double
wk82_f18_p_a_o15_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    1.66e-10 / pow( d_t9, 3./2. ) * exp( -0.302/d_t9 ) + 
    1.56e5 / pow( d_t9, 3./2. ) * exp( -3.84/d_t9 ) +
    1.36e6 / pow( d_t9, 3./2. ) * exp( -5.22/d_t9 ) +
    8.1e-5 / pow( d_t9, 3./2. ) * exp( -1.05/d_t9 ) +
    8.9e-4 / pow( d_t9, 3./2. ) * exp( -1.51/d_t9 ) +
    3.0e5 / pow( d_t9, 3./2. ) * exp( -4.29/d_t9 );

}

//##############################################################################
// la90_o17_p_g_f18_rate_function().
//##############################################################################

double
la90_o17_p_g_f18_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr << "Routine should not have extra data." << std::endl;
    exit( EXIT_FAILURE );
  }

  return 
    7.97e7 * pow( d_t9/( 1. + 2.69 * d_t9 ), 5./6. ) / pow( d_t9, 3./2. ) *
    exp( -16.712/ pow( d_t9/( 1. + 2.69 * d_t9 ), 1./3. ) ) + 1.51e8 / 
    pow( d_t9, 2./3. ) * exp( -16.712/ pow( d_t9, 1./3. ) ) * ( 1. + 0.025 * 
       pow( d_t9, 1./3. ) - 0.051 *  pow( d_t9, 2./3. ) - 8.82e-3 * d_t9 ) +
    1.56e5 / d_t9 * exp( -6.272/d_t9 ) + 0.2 * 3.16e-5 / pow( d_t9, 3./2. ) *
    exp( -0.767/d_t9 );

}

//##############################################################################
// register_rate_functions().
//##############################################################################

void
register_rate_functions( Libnucnet__Reac * p_reac )
{

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my n15(p,a)c12 rate",
    (Libnucnet__Reaction__userRateFunction) hfcz_n15_p_a_c12_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my n14(a,g)f18 rate",
    (Libnucnet__Reaction__userRateFunction) fcz2_n14_a_g_f18_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my n15(p,g)o16 rate",
    (Libnucnet__Reaction__userRateFunction) fcz2_n15_p_g_o16_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my o18(a,g)ne22 rate",
    (Libnucnet__Reaction__userRateFunction) dh03_o18_a_g_ne22_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my o18(n,g)o19 rate",
    (Libnucnet__Reaction__userRateFunction) ra94_o18_n_g_o19_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my f18(p,g)ne19 rate",
    (Libnucnet__Reaction__userRateFunction) wk82_f18_p_g_ne19_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my f18(p,a)o15 rate",
    (Libnucnet__Reaction__userRateFunction) wk82_f18_p_a_o15_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    "my o17(p,g)f18 rate",
    (Libnucnet__Reaction__userRateFunction) la90_o17_p_g_f18_rate_function
  );

}

} // namespace my_user
