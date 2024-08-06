//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer and Michael Bojazi.
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
//! \brief Code for useful utilities for nuclear_decay.
////////////////////////////////////////////////////////////////////////////////
  
#include "nuclear_decay_utilities.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */

namespace user
{

//##############################################################################
// get_decay_matrices().
//##############################################################################

std::pair< WnMatrix *, WnMatrix * >
get_decay_matrices(
  Libnucnet__Net * p_net,
  const char * s_reac_xpath,
  double d_cutoff_decay_rate
)
{

  Libnucnet__NetView * p_net_view;
  Libnucnet__Species * p_reactant, * p_product;
  double d_rate;
  std::pair< WnMatrix *, WnMatrix * > my_matrices;
  WnMatrix * p_work;
  gsl_vector * p_vector;
  
  my_matrices.first =
    WnMatrix__new(
      Libnucnet__Nuc__getNumberOfSpecies( Libnucnet__Net__getNuc( p_net ) ),
      Libnucnet__Nuc__getNumberOfSpecies( Libnucnet__Net__getNuc( p_net ) )
    );

  p_work =
    WnMatrix__new(
      Libnucnet__Nuc__getNumberOfSpecies( Libnucnet__Net__getNuc( p_net ) ),
      Libnucnet__Nuc__getNumberOfSpecies( Libnucnet__Net__getNuc( p_net ) )
    );

  p_net_view =
    Libnucnet__NetView__new(
      p_net,
      "",
      s_reac_xpath
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    d_rate =
      Libnucnet__Reaction__computeRate(
        reaction.getNucnetReaction(),
        1.e-10,
        NULL
      );

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list(
        reaction.getNucnetReaction()
      );

    if( reactant_list.size() == 1 && d_rate > d_cutoff_decay_rate )
    {

      p_reactant =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( p_net ),
          Libnucnet__Reaction__Element__getName(
            (reactant_list.begin())->getNucnetReactionElement()
          )
        );

      WnMatrix__assignElement(
        p_work,
        Libnucnet__Species__getIndex( p_reactant ) + 1,
        Libnucnet__Species__getIndex( p_reactant ) + 1,
        d_rate
      );

      nnt::reaction_element_list_t product_list =
        nnt::make_reaction_nuclide_product_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement product, product_list )
      {

        p_product =
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc( p_net ),
            Libnucnet__Reaction__Element__getName(
              product.getNucnetReactionElement()
            )
          );

        WnMatrix__assignElement(
          p_work,
          Libnucnet__Species__getIndex( p_product ) + 1,
          Libnucnet__Species__getIndex( p_reactant ) + 1,
          -d_rate *
             (double) Libnucnet__Species__getA( p_product ) /
             (double) Libnucnet__Species__getA( p_reactant )
        );

      }

    }

  }

  Libnucnet__NetView__free( p_net_view );

  p_vector = WnMatrix__getDiagonalElements( p_work );

  for( size_t i = 0; i < p_vector->size; i++ )
  {

    if( WnMatrix__value_is_zero( gsl_vector_get( p_vector, i ) ) )
      WnMatrix__assignElement( my_matrices.first, i + 1, i + 1, 1. );

  }

  gsl_vector_free( p_vector );

  WnMatrix__addValueToDiagonals( p_work, 1.e-300 );

  my_matrices.second = WnMatrix__getTransferMatrix( p_work );

  WnMatrix__free( p_work );

  WnMatrix__scaleMatrix( my_matrices.second, -1. );

  return my_matrices;

} 

//##############################################################################
// decay_abundances().
//##############################################################################

/**
 \brief Routine to decay the abundances in the input network over the input time
        using the decay rates included in the network.

 \param p_nucnet A libnucnet reaction network.
 \param d_decay_time The time interval (in s) over which to decay the
                     abundances.
 \param s_reac_xpath A reaction XPath expression selecting the decay reactions
                     to use (optional but required if debug
                     present--default = "[count(reactant) = 1]").
 \param debug A string ("yes" or "no") determining whether to print out
              solution diagnostics (optional--default = "no").
 \return On successful return, the abundances have been decayed.
*/

void
decay_abundances(
  Libnucnet * p_nucnet,
  double d_decay_time,
  const char * s_reac_xpath,
  std::string debug
)
{

  WnMatrix * p_matrix;
  std::vector<nnt::Zone> zone_vector;

  Libnucnet__setZoneCompareFunction(
    p_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
    zone_vector.push_back( zone );
  }

  Libnucnet__Zone__updateNetView(
    zone_vector[0].getNucnetZone(),
    EVOLUTION_NETWORK,
    NULL,
    NULL,
    Libnucnet__NetView__copy(
      zone_vector[0].getNetView(
        "",
        s_reac_xpath
      )
    )
  ); 

  Libnucnet__Zone__computeRates(
    zone_vector[0].getNucnetZone(),
    1.e-10,
    1.e-10
  );

  p_matrix =
    Libnucnet__Zone__computeJacobian( zone_vector[0].getNucnetZone() );

  // Return if matrix empty--no decay

  if( WnMatrix__getNumberOfElements( p_matrix ) == 0 )
  {
    return;
  }
  
#ifndef NO_OPENMP
  #pragma omp parallel for schedule( dynamic, 1 )
#endif
    for( size_t i = 0; i < zone_vector.size(); i++ )
    {

      gsl_vector * p_vector =
        Libnucnet__Zone__getAbundances( zone_vector[i].getNucnetZone() );

      gsl_vector * p_const = gsl_vector_calloc( p_vector->size );
          
      WnSparseSolve__Phi * p_exp_solver = WnSparseSolve__Phi__new();

      if( debug == "yes" ) WnSparseSolve__Phi__setDebug( p_exp_solver );

      WnSparseSolve__Phi__updateMaximumIterations( p_exp_solver, 100000 );

      gsl_vector * p_sol =
	WnSparseSolve__Phi__solve(
	  p_exp_solver,
	  p_matrix,
	  p_vector,
	  p_const,
	  d_decay_time
	);

      if( !p_sol )
      {
	std::cerr << "Solution not found!" << std::endl;
	exit( EXIT_FAILURE );
      }

      Libnucnet__Zone__updateAbundances(
	zone_vector[i].getNucnetZone(),
	p_sol
      );

      gsl_vector_free( p_sol );
      
      if( debug == "yes" )
      {
	std::cout <<
	  Libnucnet__Zone__getLabel(
	    zone_vector[i].getNucnetZone(),
	    1
	  ) << std::endl;
      }

      gsl_vector_free( p_const );
      gsl_vector_free( p_vector );

      WnSparseSolve__Phi__free( p_exp_solver );

      zero_out_small_abundances( zone_vector[i], 0. );

    }

  WnMatrix__free( p_matrix );

} 

void
decay_abundances(
  Libnucnet * p_nucnet,
  double d_decay_time,
  const char * s_reac_xpath
)
{

  decay_abundances( p_nucnet, d_decay_time, s_reac_xpath, "no" );

}

void
decay_abundances( Libnucnet * p_nucnet, double d_decay_time )
{

  decay_abundances( p_nucnet, d_decay_time, S_DECAY_XPATH, "no" );

}

//##############################################################################
// push_abundances_to_daughters().
//##############################################################################

/**
 \brief Routine to push abundances in zones in a network with rates greater
        than the cutoff rate to their daughters with decay rate less than
        the cutoff rate.

 \param p_nucnet A libnucnet reaction network.
 \param d_cutoff_decay_rate The cutoff rate.
 \param s_reac_xpath A reaction XPath expression selecting the decay reactions
                     to use (optional--default = "[count(reactant) = 1]").
 \return On successful return, the abundances have been pushed to their
         appropriate daughters.
*/
void
push_abundances_to_daughters(
  Libnucnet * p_nucnet,
  double d_cutoff_decay_rate,
  const char * s_reac_xpath
)
{

  std::pair< WnMatrix *, WnMatrix * > my_matrices;
  std::vector<nnt::Zone> zone_vector;
  double d_mass;

  my_matrices =
    get_decay_matrices(
      Libnucnet__getNet( p_nucnet ),
      s_reac_xpath,
      d_cutoff_decay_rate
    );

  Libnucnet__setZoneCompareFunction(
    p_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_nucnet ) )
  {
    zone_vector.push_back( zone );
  }

#ifndef NO_OPENMP
  #pragma omp parallel for schedule( dynamic, 1 )
#endif
    for( size_t i = 0; i < zone_vector.size(); i++ )
    {

      gsl_vector * p_a =
        Libnucnet__Zone__getMassFractions( zone_vector[i].getNucnetZone() );

      gsl_vector * p_b =
        gsl_vector_alloc(
          WnMatrix__get_gsl_vector_size( p_a )
        );

      gsl_vector_memcpy( p_b, p_a );

      gsl_vector * p_undecayed =
        gsl_vector_alloc(
          WnMatrix__get_gsl_vector_size( p_a )
        );

      gsl_vector_memcpy( p_undecayed, p_a );

      for( size_t j = 0; j < 1000; j++ )
      {

        gsl_vector * p_c =
          WnMatrix__computeMatrixTimesVector( my_matrices.second, p_b );
        gsl_vector_add( p_a, p_c );
        gsl_vector_memcpy( p_b, p_c );
        gsl_vector_free( p_c );

        if( WnMatrix__value_is_zero( gsl_blas_dnrm2( p_b ) ) ) break;

      }

      gsl_vector * p_c =
        WnMatrix__computeMatrixTimesVector( my_matrices.first, p_a );

      d_mass = 0.;
      for( size_t j = 0; j < WnMatrix__get_gsl_vector_size( p_c ); j++ )
      {
        if( gsl_vector_get( p_c, j ) > 1.e-20 )
          d_mass +=
            gsl_vector_get( p_c, j ) -
            gsl_vector_get( p_undecayed, j );
      }

      zone_vector[i].updateProperty(
        S_RADIOACTIVE_MASS_FRACTION,
        boost::lexical_cast<std::string>( d_mass )
      );

      Libnucnet__Zone__updateMassFractions(
        zone_vector[i].getNucnetZone(), p_c
      );

      gsl_vector_free( p_a );
      gsl_vector_free( p_b );
      gsl_vector_free( p_c );
      gsl_vector_free( p_undecayed );

    }

  WnMatrix__free( my_matrices.first );
  WnMatrix__free( my_matrices.second );

}

void
push_abundances_to_daughters(
  Libnucnet * p_nucnet,
  double d_cutoff_decay_rate
)
{
  push_abundances_to_daughters( p_nucnet, d_cutoff_decay_rate, S_DECAY_XPATH );
}

} // namespace user
