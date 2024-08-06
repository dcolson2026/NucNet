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
//! \brief Code for outputting the time, t9, and rho trajectory for a
//!        simple snII zone.
////////////////////////////////////////////////////////////////////////////////

#include "nnt/iter.h"
#include "simple_snII_utilities.h"

int main( int argc, char * argv[] )
{
 
  Libnucnet * p_my_nucnet, * p_special;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 3 )
  {
    fprintf(
      stderr,
      "\nUsage: %s input_xml zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  input_file = input data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_xpath = XPath expression to retrieve zone\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[1],
      NULL
    );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[1], NULL );

  p_special = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_special ) ),
      argv[1],
      NULL
    );

  Libnucnet__assignZoneDataFromXml( p_special, argv[1], argv[2] );

  if( Libnucnet__getNumberOfZones( p_special ) != 1 )
  {
    std::cerr << "Zone xpath must select a single zone." << std::endl;
    return EXIT_FAILURE;
  }
  
  nnt::zone_list_t zone_list = nnt::make_zone_list( p_special );

  boost::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
  > t =
    get_zone_vectors( *zone_list.begin(), p_my_nucnet );

  for( size_t i = 0; i < t.get<0>().size(); i++ )
  {

     std::cout <<
        t.get<0>().at( i ) <<
        "  " <<
        t.get<1>().at( i ) <<
        "  " <<
        pow( 10., t.get<2>().at( i ) ) << std::endl;
  
  }       
    
  Libnucnet__free( p_special );
  Libnucnet__free( p_my_nucnet );
 
  return EXIT_SUCCESS;

}  
