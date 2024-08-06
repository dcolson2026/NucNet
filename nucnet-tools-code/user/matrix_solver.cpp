//////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Code for solving the matrix for a zone.
////////////////////////////////////////////////////////////////////////////////

#include "matrix_solver.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */
namespace user
{

#ifdef SPARSKIT2

//##############################################################################
// default_sparse_solver_properties().
//##############################################################################

std::pair<size_t, double>
default_sparse_solver_properties(
  WnSparseSolve__Mat * p_solver,
  nnt::Zone& zone
)
{

  size_t i_delta;
  double d_droptol;

  if( zone.hasProperty( nnt::s_ITER_SOLVER ) )
  {
    WnSparseSolve__Mat__updateSolverMethod(
      p_solver,
      zone.getProperty<std::string>( nnt::s_ITER_SOLVER ).c_str()
    );
  }
  else
  {
    std::cerr << "No sparse solver set." << std::endl;
    exit( EXIT_FAILURE );
  }

  if( zone.hasProperty( nnt::s_ITER_SOLVER_MAX_ITERATIONS ) )
  {
    WnSparseSolve__Mat__updateMaximumIterations(
      p_solver,
      zone.getProperty<int>( nnt::s_ITER_SOLVER_MAX_ITERATIONS )
    );
  }

  if( zone.hasProperty( nnt::s_ITER_SOLVER_REL_TOL ) )
  {
    WnSparseSolve__Mat__updateRelativeTolerance(
      p_solver,
      zone.getProperty<double>( nnt::s_ITER_SOLVER_REL_TOL )
    );
  }
  else
  {
    WnSparseSolve__Mat__updateRelativeTolerance(
      p_solver,
      1.e-16
    );
  }

  if( zone.hasProperty( nnt::s_ITER_SOLVER_ABS_TOL ) )
  {
    WnSparseSolve__Mat__updateAbsoluteTolerance(
      p_solver,
      zone.getProperty<double>( nnt::s_ITER_SOLVER_ABS_TOL )
    );
  }
  else
  {
    WnSparseSolve__Mat__updateAbsoluteTolerance(
      p_solver,
      0.
    );
  }

  if( zone.hasProperty( nnt::s_ITER_SOLVER_CONVERGENCE_METHOD ) )
  {
    WnSparseSolve__Mat__updateConvergenceMethod(
      p_solver,
      zone.getProperty<std::string>(
        nnt::s_ITER_SOLVER_CONVERGENCE_METHOD
      ).c_str()
    );
  }

  if(
    zone.hasProperty( nnt::s_ITER_SOLVER_DEBUG ) &&
    zone.getProperty<std::string>( nnt::s_ITER_SOLVER_DEBUG ) == "yes"
  )
  {
    WnSparseSolve__Mat__setDebug( p_solver );
  }

  //============================================================================
  // Set the preconditioner data.
  //============================================================================

  if( zone.hasProperty( nnt::s_ILU_DELTA ) )
  {
    i_delta = zone.getProperty<size_t>( nnt::s_ILU_DELTA );
  }
  else
  {
    i_delta = 1;
  }

  if( zone.hasProperty( nnt::s_ILU_DROP_TOL ) )
  {
    d_droptol = zone.getProperty<double>( nnt::s_ILU_DROP_TOL );
  }
  else
  {
    d_droptol = 0;
  }

  return std::make_pair( i_delta, d_droptol );

}

#endif

//##############################################################################
// solve_matrix_for_zone().
//##############################################################################

gsl_vector *
solve_matrix_for_zone(
  nnt::Zone& zone,
  WnMatrix * p_matrix,
  gsl_vector * p_rhs
)
{

  gsl_vector * p_sol;
  WnMatrix__Arrow * p_arrow;

#ifndef SPARSKIT2
  if(
    zone.hasProperty( nnt::s_ITER_SOLVER ) ||
    zone.hasProperty( nnt::s_ITER_SOLVER_T9 )
  )
  {
    std::cerr << std::endl;
    std::cerr << "Code was not compiled against Sparskit2. Either recompile"
       << std::endl;
    std::cerr << "the code with Sparskit2 or remove the iterative solver zone"
       << std::endl;
    std::cerr << "properties." << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }
#endif

  //============================================================================
  // Call optional matrix function.
  //============================================================================

  if( zone.hasFunction( nnt::s_MATRIX_MODIFICATION_FUNCTION ) )
  {
    boost::any_cast<boost::function<void(WnMatrix*, gsl_vector* )> >(
      zone.getFunction( nnt::s_MATRIX_MODIFICATION_FUNCTION )
    )( p_matrix, p_rhs );
  }

#ifdef SPARSKIT2

  //============================================================================
  // Try iterative solver, if set.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_ITER_SOLVER ) &&
    zone.hasProperty( nnt::s_ITER_SOLVER_T9 ) &&
    zone.getProperty<double>( nnt::s_T9 )
    <
    zone.getProperty<double>( nnt::s_ITER_SOLVER_T9 )
  )
  {

    //--------------------------------------------------------------------------
    // Get solution.
    //--------------------------------------------------------------------------
    
    #pragma omp critical
    {

    if( zone.hasProperty( nnt::s_SOLVER_PARAMETER_FUNCTION ) )
    {
      p_sol =
        solve_sparse_matrix_with_ilu_preconditioner(
          p_matrix,
          p_rhs,
          boost::bind(
            boost::any_cast<
              boost::function<std::pair<size_t, double>( WnSparseSolve__Mat * )>
            >(
                zone.getFunction( nnt::s_SOLVER_PARAMETER_FUNCTION )
            ),
            _1
          )
        );
    }
    else
    {
      p_sol =
        solve_sparse_matrix_with_ilu_preconditioner(
          p_matrix,
          p_rhs,
          boost::bind(
            default_sparse_solver_properties, _1, boost::ref( zone )
          )
        );
    }

    }

    if( p_sol ) return p_sol;

  }

#endif

  //============================================================================
  // Normal solver.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_SOLVER ) &&
    zone.getProperty<std::string>( nnt::s_SOLVER ) == nnt::s_ARROW
  )
  {
    p_arrow =
      WnMatrix__getArrow(
        p_matrix,
        zone.getProperty<unsigned long>( nnt::s_ARROW_WIDTH )
      );
    p_sol = WnMatrix__Arrow__solve( p_arrow, p_rhs );
    WnMatrix__Arrow__free( p_arrow );
  }
  else
  {
    p_sol = WnMatrix__solve( p_matrix, p_rhs );
  }

  return p_sol;

}

#ifdef SPARSKIT2

//##############################################################################
// phi__solve__parallel().
//##############################################################################

gsl_vector *
phi__solve__parallel(
  WnSparseSolve__Phi *self,
  std::vector<WnMatrix *>& matrices,
  const gsl_vector *p_input_vector,
  const gsl_vector *p_const_vector,
  double d_tn
)
{

  size_t i, i_rows, n_fpar = 6;
  int i_indic, i_ierr, i_iter, i_workspace;
  int i_ipar[5] = {0,0,0,0,0};
  double d_tol;
  gsl_vector *p_x, *p_y, *p_u, *p_output_vector, *p_fpar;

  //============================================================================
  // Check input.
  //============================================================================

  if( matrices.size() == 0 || !p_input_vector || !p_const_vector )
    WNSPARSESOLVE__ERROR( "Invalid input" );

  //============================================================================
  // Find debugging information.
  //============================================================================

  if( self->iDebug == 1 )
      printf( "\nUsing phipro\n\n" );

  //============================================================================
  // Set debugging flag.
  //============================================================================

  i_ipar[3] = self->iDebug;

  //============================================================================
  // Get number of rows.
  //============================================================================

  i_rows = WnMatrix__getNumberOfRows( matrices[0] );

  //============================================================================
  // Set local variables from structure.
  //============================================================================

  i_workspace = self->iWorkSpace;
  d_tol = self->dTolerance;

  //============================================================================
  // Vector to handle multiplication results.
  //============================================================================

  std::vector<gsl_vector *> xx( matrices.size() );

  //============================================================================
  // Allocate memory for vectors.
  //============================================================================

  p_x = gsl_vector_calloc( i_rows );
  p_y = gsl_vector_calloc( i_rows );
  p_u = gsl_vector_calloc( i_rows * ( (size_t) i_workspace + 1 ) );
  p_output_vector = gsl_vector_alloc( p_input_vector->size );
  p_fpar = gsl_vector_calloc( n_fpar );

  if( !p_x || !p_y || !p_u || !p_output_vector || !p_fpar )
      WNSPARSESOLVE__ERROR( "Couldn't allocate memory in WnSparseSolve" );

  gsl_vector_memcpy( p_output_vector, p_input_vector );

  i_indic = 0;
  i_iter = 0;
  i_ierr = 0;

  while( i_indic != 1 ) {

    if( i_iter == self->iIterMax ) {
      if( self->iDebug == 1 )
        fprintf( stdout, "Exceeded maximum number of iterations.\n" );
      gsl_vector_free( p_x );
      gsl_vector_free( p_y );
      gsl_vector_free( p_u );
      gsl_vector_free( p_output_vector );
      gsl_vector_free( p_fpar );
      return NULL;
    }

    phipro_(
      &i_rows,
      &i_workspace,
      &d_tol,
      &d_tn,
      p_output_vector->data,
      p_const_vector->data,
      p_u->data,
      p_x->data,
      p_y->data,
      p_fpar->data,
      i_ipar,
      &i_indic,
      &i_ierr
    );

    gsl_vector_free( p_y );

#ifndef NO_OPENMP
    #pragma omp parallel for schedule( dynamic, 1 )
#endif
      for( i = 0; i < matrices.size(); i++ )
      {
        xx[i] =
          WnMatrix__computeMatrixTimesVector(
            matrices[i],
            p_x
          );
      }

    size_t i_add = 2;

    size_t i_number = matrices.size();

    while( i_number / 2 > 0 )
    {

      size_t i_c = i_add / 2;

      if( GSL_IS_ODD( i_number ) )
      {
        gsl_vector_add(
          xx[(i_number/2)*i_add - i_c],
          xx[(i_number/2)*i_add]
        );
        gsl_vector_free( xx[(i_number/2)*i_add] );
      }

      i_number /= 2;

#ifndef NO_OPENMP
      #pragma omp parallel for schedule( dynamic, 1 )
#endif
        for( i = 0; i < i_number; i++ )
        {
          gsl_vector_add( xx[i*i_add], xx[i*i_add+i_c] );
          gsl_vector_free( xx[i*i_add+i_c] );
        }

      i_add *= 2;

    }

    p_y = gsl_vector_calloc( i_rows );

    gsl_vector_memcpy( p_y, xx[0] );

    gsl_vector_free( xx[0] );

    gsl_vector_scale( p_y, -1. );

    ++i_iter;

  }

  if( self->iDebug == 1 )
    printf( "%d iterations.\n", i_iter );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  gsl_vector_free( p_x );
  gsl_vector_free( p_y );
  gsl_vector_free( p_u );
  gsl_vector_free( p_fpar );

  return p_output_vector;

}

#endif

} //namespace user
