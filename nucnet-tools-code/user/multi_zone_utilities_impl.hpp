////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2013-2014 Clemson University.
//
// This file was originally written by Bradley S. Meyer.
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
//! \file multi_zone_utilities_impl.hpp
//! \brief A header file to define useful multi-zone utilities.
//!
////////////////////////////////////////////////////////////////////////////////

#include <boost/any.hpp>
#include <boost/function.hpp>
#include <boost/bind/bind.hpp>

#ifndef MULTI_ZONE_UTILITIES_IMPL_HPP
#define MULTI_ZONE_UTILITIES_IMPL_HPP

namespace user
{

//##############################################################################
// multi_zone_vertex_writer().
//##############################################################################

struct
multi_zone_vertex_writer
{
  multi_zone_vertex_writer( user::zone_link_graph_t& g_ ) :          g( g_ ) {};

  template <class Vertex>
  void operator()( std::ostream& out, Vertex v )
  {
    out <<
      boost::format( "[label=\"(%s, %s, %s)\" ]" ) %
        Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 1 ) %
        Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 2 ) %
        Libnucnet__Zone__getLabel( g[v].getNucnetZone(), 3 );
  }

  user::zone_link_graph_t& g;
};

//##############################################################################
// multi_zone_edge_writer().
//##############################################################################

struct
multi_zone_edge_writer
{
  multi_zone_edge_writer( user::zone_link_graph_t& g_ ) :           g( g_ ) {};

  template <class Edge>
  void operator()( std::ostream& out, Edge e )
  {
    out << boost::format( "[label=\"%.2e\" ]" ) % g[e].getWeight();
  }

  user::zone_link_graph_t& g;
};

//##############################################################################
// exp_multi_zone().
//##############################################################################

template<class Function1, class Function2, class Function3, class Function4>
bool
exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 phi_fn,
  Function2 param_fn,
  Function3 check_fn,
  Function4 const_vec_fn,
  double d_t
)
{

  WnSparseSolve__Phi * p_phi;
  gsl_vector * p_sol, * p_old, * p_prev, * p_constant;
  std::vector<WnMatrix *> mix_matrices;
  WnMatrix * p_jacobian_matrix;
  double d_dt_e, d_dt_cum;
  size_t i, i_species, i_offset, i_delta;
  int i_check_fn;

  //============================================================================
  // Create solver and set parameters.
  //============================================================================

  p_phi = WnSparseSolve__Phi__new();

  //============================================================================
  // Get parameters.
  //============================================================================

  param_fn( p_phi );

  //============================================================================
  // Get number of species and abundance vectors.
  //============================================================================

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  i_offset = i_species;

  if( is_multi_mass_calculation( zones ) )
    i_delta = 1;
  else
    i_delta = 0;

  i_offset += i_delta;

  p_old = create_full_vector( zones );

  try{
    p_constant = const_vec_fn();
  }
  catch( boost::bad_function_call& ex )
  {
    p_constant = gsl_vector_calloc( i_offset * zones.size() );
  }

  //============================================================================
  // Exponential solution.
  //============================================================================

  mix_matrices = phi_fn();

  p_jacobian_matrix = get_evolution_matrix( zones[0] );
  WnMatrix__scaleMatrix( p_jacobian_matrix, -1. );

#ifndef NO_OPENMP
  #pragma omp parallel for schedule( dynamic, 1 )
#endif
    for( i = 0; i < mix_matrices.size(); i++ )
    {
      WnMatrix__insertMatrix(
        mix_matrices[i],
        p_jacobian_matrix,
        (i * i_offset) + i_delta + 1,
        (i * i_offset ) + i_delta + 1
      );

    }

  WnMatrix__free( p_jacobian_matrix );

  //============================================================================
  // Exponential solution.
  //============================================================================

  d_dt_e = d_t;

  d_dt_cum = 0;

  p_prev = gsl_vector_alloc( p_old->size );

  gsl_vector_memcpy( p_prev, p_old );

  while( gsl_fcmp( d_dt_cum, d_t, D_DT_E_CMP ) )
  {

    p_sol =
      phi__solve__parallel(
        p_phi,
        mix_matrices,
        p_prev,
        p_constant,
        d_dt_e
      );

    if( p_sol )
    {

      update_from_full_vector(
        zones,
        p_sol,
        "abundances"
      );

      try
      {
        i_check_fn = check_fn();
      }
      catch( boost::bad_function_call& ex )
      {
        i_check_fn = 1;
      }

    }

    if( !p_sol || i_check_fn != 1 )
    {
      update_from_full_vector(
        zones,
        p_prev,
        "abundances"
      );
      d_dt_e /= D_DT_E_STEP;
    }
    else
    {
      d_dt_cum += d_dt_e;
      gsl_vector_memcpy( p_prev, p_sol );
      d_dt_e = d_t - d_dt_cum;
      std::cout << d_dt_cum << "  " << d_t << std::endl;
    }
      
    gsl_vector_free( p_sol );

  }

  for( size_t i = 0; i < mix_matrices.size(); i++ )
    WnMatrix__free( mix_matrices[i] );

  WnSparseSolve__Phi__free( p_phi );

  gsl_vector_free( p_constant );
  gsl_vector_free( p_old );
  gsl_vector_free( p_prev );

  return true;

}

template<class Function1, class Function2, class Function3>
int
exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 phi_fn,
  Function2 param_fn,
  Function3 check_fn,
  double d_t
)
{

  boost::function< gsl_vector*( ) > const_vec_fn;
  const_vec_fn = 0;

  return
    exp_multi_zone(
      zones,
      phi_fn,
      param_fn,
      check_fn,
      boost::bind( const_vec_fn ),
      d_t
    );

}

template<class Function, class Function2>
int
exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function phi_fn,
  Function2 param_fn,
  double d_t
)
{

  boost::function< int() > check_fn;
  boost::function< gsl_vector*( ) > const_vec_fn;
  check_fn = 0;
  const_vec_fn = 0;

  return
    exp_multi_zone(
      zones,
      phi_fn,
      param_fn,
      boost::bind( check_fn ),
      boost::bind( const_vec_fn ),
      d_t
    );

}

template<class Function>
int
exp_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function phi_fn,
  double d_t
)
{

  boost::function< void() > param_fn;
  boost::function< int() > check_fn;
  boost::function< gsl_vector*( ) > const_vec_fn;
  param_fn = 0;
  check_fn = 0;
  const_vec_fn = 0;

  return
    exp_multi_zone(
      zones,
      phi_fn,
      boost::bind( param_fn ),
      boost::bind( check_fn ),
      boost::bind( const_vec_fn ),
      d_t
    );

}

//##############################################################################
// evolve_multi_zone().
//##############################################################################

template<class Function1, class Function2, class Function3>
bool
evolve_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 mat_vec_fn,
  Function2 param_fn,
  Function3 check_fn,
  double d_dt
)
{

  WnMatrix * p_matrix;
  gsl_vector * p_vector, * p_sol, * p_old, * p_current, * p_work;
  gsl_vector_view view;
  std::vector<WnMatrix *> jacobians;
  std::vector<gsl_vector *> rhs_vectors;
  size_t i, i_species, i_iter, i_offset, i_delta;
  int i_check;
  double d_check;

  //============================================================================
  // Get number of species and abundance vectors.
  //============================================================================

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  i_offset = i_species;

  if( is_multi_mass_calculation( zones ) )
    i_delta = 1;
  else
    i_delta = 0;

  i_offset += i_delta;

  p_old = create_full_vector( zones );

  p_current = gsl_vector_alloc( WnMatrix__get_gsl_vector_size( p_old ) );

  gsl_vector_memcpy( p_current, p_old );

  //============================================================================
  // Newton-Raphson loop.
  //============================================================================

  for( i_iter = 0; i_iter < IT_MAX; i_iter++ )
  {

    boost::tie( p_matrix, p_vector ) = mat_vec_fn( p_current );

    boost::tie( jacobians, rhs_vectors ) =
      get_zone_jacobians_and_rhs_vectors( zones );

    i = 0;
    BOOST_FOREACH( WnMatrix * p_sub_matrix, jacobians )
    {
      
      WnMatrix__insertMatrix(
	p_matrix,
	p_sub_matrix,
	(i * i_offset) + i_delta + 1,
	(i * i_offset ) + i_delta + 1
      );

      WnMatrix__free( p_sub_matrix );

      i++;

    }

    WnMatrix__addValueToDiagonals( p_matrix, 1. / d_dt );

    i = 0;
    BOOST_FOREACH( gsl_vector * p_sub_vector, rhs_vectors )
    {
      
      view =
        gsl_vector_subvector( p_vector, ( i * i_offset ) + i_delta, i_species );
      gsl_vector_add( &view.vector, p_sub_vector );
      gsl_vector_free( p_sub_vector );
      i++;

    }

    p_work = gsl_vector_alloc( WnMatrix__get_gsl_vector_size( p_current ) );
    gsl_vector_memcpy( p_work, p_current );
    gsl_vector_sub( p_work, p_old );
    gsl_vector_scale( p_work, 1. / d_dt );
    gsl_vector_sub( p_vector, p_work );
    gsl_vector_free( p_work );

    p_sol =
      solve_sparse_matrix_with_ilu_preconditioner(
        p_matrix,
        p_vector,
        param_fn
      ); 

    WnMatrix__free( p_matrix );
    gsl_vector_free( p_vector );

    if( !p_sol ) 
    {
      update_from_full_vector(
        zones,
        p_old,
        "abundances"
      );
      gsl_vector_free( p_old );
      gsl_vector_free( p_current );
      return false;
    }

    gsl_vector_add( p_current, p_sol );

    update_from_full_vector(
      zones,
      p_current,
      "abundances"
    );

    d_check =
      sqrt( gsl_blas_dnrm2( p_sol ) ) / sqrt( gsl_blas_dnrm2( p_current ) );

    gsl_vector_free( p_sol );

    if( d_check < 1.e-4 ) break;

  }

  try
  {
    i_check = check_fn();
  }
  catch( boost::bad_function_call& ex )
  {
    i_check = 1;
  }

  if( i_check == 0 )
  {
    update_from_full_vector(
      zones,
      p_old,
      "abundances"
    );
    gsl_vector_free( p_old );
    gsl_vector_free( p_current );
    return false;
  }

  gsl_vector_free( p_old );
  gsl_vector_free( p_current );

  return true;

}

template<class Function1, class Function2>
bool
evolve_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 mat_vec_fn,
  Function2 param_fn,
  double d_dt
)
{

  boost::function< int() > f;
  f = 0;

  return
    evolve_multi_zone(
      zones,
      mat_vec_fn,
      param_fn,
      boost::bind( f ),
      d_dt
    );

}

//##############################################################################
// safe_evolve_multi_zone().
//##############################################################################

template<class Function1, class Function2, class Function3>
void
safe_evolve_multi_zone(
  std::vector<nnt::Zone>& zones,
  Function1 mat_vec_fn,
  Function2 param_fn,
  Function3 check_fn,
  double d_dt
)
{

  double d_dt_e, d_dt_cum = 0.;

  gsl_vector * p_old, * p_new;

  p_old = create_full_vector( zones );

  if( !evolve_multi_zone( zones, mat_vec_fn, param_fn, check_fn, d_dt ) )
  {

    d_dt_e = d_dt / D_DT_E_STEP;

    while( !evolve_multi_zone( zones, mat_vec_fn, param_fn, check_fn, d_dt_e ) )
    {
      d_dt_e /= D_DT_E_STEP;
    }

    d_dt_cum = d_dt_e;

    while( gsl_fcmp( d_dt_cum, d_dt, D_DT_E_CMP ) )
    { 

      std::cout << d_dt_cum << "  " << d_dt << std::endl;

      if( evolve_multi_zone( zones, mat_vec_fn, param_fn, d_dt_e ) )
      {
	d_dt_cum += d_dt_e;
	d_dt_e *= D_DT_E_REG;
	if( d_dt_cum + d_dt_e > d_dt ) d_dt_e = d_dt - d_dt_cum;
      }
      else
      {
	d_dt_e /= D_DT_E_REG;
      }

    }

  }

  p_new = create_full_vector( zones );

  gsl_vector_sub( p_new, p_old );

  update_from_full_vector(
    zones,
    p_new,
    "abundance changes"
  );

  gsl_vector_free( p_old );
  gsl_vector_free( p_new );

}


} // namespace user

#endif // MULTI_ZONE_UTILITIES_IMPL_HPP
