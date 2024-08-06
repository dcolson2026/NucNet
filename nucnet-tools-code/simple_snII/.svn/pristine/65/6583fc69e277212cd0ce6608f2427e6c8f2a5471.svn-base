////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Code for outputting the t9 and rho for all zones for selected times
//!        in a simple snII model.
////////////////////////////////////////////////////////////////////////////////

#include <boost/format.hpp>

#include "simple_snII_utilities.h"

int main( int argc, char * argv[] )
{
 
  Libnucnet * p_my_nucnet, * p_specific;
  nnt::Zone my_zone;

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
      stderr, "  zone_xpath = XPath expression to select zones\n\n"
    );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Get nucnets--the full one and one giving the requested zones.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[1], NULL );

  p_specific = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_specific ) ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( p_specific, argv[1], argv[2] );

  //============================================================================
  // Get zone lists.  Output total number of zones in model.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    ( Libnucnet__Zone__compare_function ) nnt::zone_compare_by_first_label
  );

  Libnucnet__setZoneCompareFunction(
    p_specific,
    ( Libnucnet__Zone__compare_function ) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list1 = nnt::make_zone_list( p_my_nucnet );

  nnt::zone_list_t zone_list2 = nnt::make_zone_list( p_specific );

  std::cout << zone_list1.size() << std::endl;

  //============================================================================
  // If requesting all zones, add an initial output with pre-shock conditions.
  //============================================================================

  if( zone_list1.size() == zone_list2.size() )
  {

    std::cout << boost::format( "t (s) = %.5e\n") % 0.;

    BOOST_FOREACH( nnt::Zone zone, zone_list1 )
    {
      std::cout <<
        zone.getProperty<std::string>( nnt::s_T9_0 ) <<
        "  " <<
        zone.getProperty<std::string>( nnt::s_RHO_0 ) << std::endl;
    }

    std::cout << std::endl;

  }

  //============================================================================
  // Loop over selected zones.
  //============================================================================

  BOOST_FOREACH( nnt::Zone shock_zone, zone_list2 )
  {

    double d_shock_time = shock_zone.getProperty<double>( S_SHOCK_TIME );

    std::cout <<
      boost::format( "t (s) = %.5e\n") % d_shock_time;

    BOOST_FOREACH( nnt::Zone zone, zone_list1 )
    {

      double d_t = zone.getProperty<double>( S_SHOCK_TIME );

      if( d_t > d_shock_time )
      {
        std::cout <<
          zone.getProperty<std::string>( nnt::s_T9_0 ) <<
          "  " <<
          zone.getProperty<std::string>( nnt::s_RHO_0 ) << std::endl;
      }
      else
      {
        double d_rho =
          pow(
            zone.getProperty<double>( S_SHOCK_T9 ) /
            shock_zone.getProperty<double>( S_SHOCK_T9 ),
            3.
          ) *
          shock_zone.getProperty<double>( S_SHOCK_RHO );

        std::cout <<
          shock_zone.getProperty<std::string>( S_SHOCK_T9 ) <<
          "  " <<
          d_rho << std::endl;
      }
          
    }       

    std::cout << std::endl;

  }
    
  //============================================================================
  // Done.  Do cleanup.
  //============================================================================

  Libnucnet__free( p_specific );
  Libnucnet__free( p_my_nucnet );
 
  return EXIT_SUCCESS;

}  
