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
//! \brief A header file for the zone matrix solver routines.
////////////////////////////////////////////////////////////////////////////////

#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP

#include "user/ilu_solvers.h"

namespace user
{

//##############################################################################
// solve_sparse_matrix_with_ilu_preconditioner().
//##############################################################################

template<class Function>
gsl_vector *
solve_sparse_matrix_with_ilu_preconditioner(
  WnMatrix * p_matrix,
  gsl_vector * p_rhs,
  Function param_fn
)
{

  gsl_vector * p_guess_vector, * p_sol;
  WnSparseSolve__Mat * p_solver;
  ilu_structure * p_ilu_structure;
  size_t i_delta;
  double d_droptol;

  //============================================================================
  // Get the guess vector and zero it out.
  //============================================================================

  p_guess_vector =
    gsl_vector_calloc(
      WnMatrix__get_gsl_vector_size( p_rhs )
    );

  //============================================================================
  // Get solver.  Set solver properties.
  //============================================================================

  p_solver = WnSparseSolve__Mat__new();

  boost::tie<size_t, double>( i_delta, d_droptol ) = param_fn( p_solver );

  WnSparseSolve__Mat__updatePreconditionerSolver(
    p_solver,
    (WnSparseSolve__Mat__PreconditionerSolver) ilu_solver
  );

  WnSparseSolve__Mat__updatePreconditionerTransposeSolver(
    p_solver,
    (WnSparseSolve__Mat__PreconditionerTransposeSolver)
      ilu_transpose_solver
  );

  p_ilu_structure =
    ilu_preconditioner_new(
      p_matrix,
      i_delta,
      d_droptol
    );

  WnSparseSolve__Mat__updatePreconditionerUserData(
    p_solver, p_ilu_structure
  );

  //============================================================================
  // Solve the matrix equation.
  //============================================================================

  p_sol =
    WnSparseSolve__Mat__solve(
      p_solver,
      p_matrix,
      p_rhs,
      p_guess_vector
    );

  //============================================================================
  // Clean up and return if valid solution.
  //============================================================================

  WnSparseSolve__Mat__free( p_solver );
  ilu_preconditioner_free( p_ilu_structure );
  gsl_vector_free( p_guess_vector );

  return p_sol;

}

} // namespace user

#endif // MATRIX_SOLVER_HPP
