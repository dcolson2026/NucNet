////////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code to remove duplicate reactions from a network xml file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>

#include "nnt/iter.h"

#include "user/remove_duplicate.h"

int main( int argc, char ** argv )
{

  Libnucnet * p_nucnet;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc!= 3 )
  {
 
    std::cerr << "\nUsage: " << argv[0] << " in_file type out_file\n\n";

    std::cerr << "  in_file = input data xml file\n\n";
    std::cerr << "  out_file = output xml file name\n\n";

    return EXIT_FAILURE;

  }

  //============================================================================
  // Get input.
  //============================================================================

  p_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_nucnet ) );

  //============================================================================
  // Output.  If there are no zones, output as a network file.  If there are,
  // output as a full libnucnet file.
  //============================================================================

  if( Libnucnet__getNumberOfZones( p_nucnet ) > 0 )
  {
    Libnucnet__writeToXmlFile( p_nucnet, argv[2] );
  }
  else
  {
    Libnucnet__Net__writeToXmlFile( Libnucnet__getNet( p_nucnet ), argv[2] );
  }

  //============================================================================
  // Clean up.
  //============================================================================

  Libnucnet__free( p_nucnet );

  return EXIT_SUCCESS;

}
