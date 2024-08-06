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
//! \brief Example code to remove invalid reactions from a network xml file.
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
  // Remove invalid reactions.
  //============================================================================

  Libnucnet__Net * p_net = Libnucnet__getNet( p_nucnet );

  Libnucnet__Reac * p_reac = Libnucnet__Net__getReac( p_net );

  BOOST_FOREACH( nnt::Reaction reaction, nnt::make_reaction_list( p_reac ) )
  {
    if(
      !Libnucnet__Net__isValidReaction( p_net, reaction.getNucnetReaction() )
    )
    {
      Libnucnet__Reac__removeReaction( p_reac, reaction.getNucnetReaction() );
    }
  }
  
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
