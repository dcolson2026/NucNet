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
//! \brief Example code to print all properties for a given zone in a given
//!        group.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>

#include <boost/format.hpp>

#include "user/hdf5_routines.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " input.h5 \"Step 00025\" 330 0 0" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 6 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " print properties for zone in a given group." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s hdf5_file group label1 label2 label3\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  hdf5_file = input hdf5 file\n\n"
    );
    fprintf(
      stderr,
      "  group = group label\n\n"
    );
    fprintf(
      stderr,
      "  label1 = first zone label\n\n"
    );
    fprintf(
      stderr,
      "  label2 = second zone label\n\n"
    );
    fprintf(
      stderr,
      "  label3 = third zone label\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  user::hdf5::zone_labels_bimap my_zone_labels;
  user::hdf5::zone_properties_hash my_properties_hash;
  
  check_input( argc, argv );

  my_properties_hash =
    user::hdf5::get_zone_properties(
      argv[1],
      argv[2],
      argv[3],
      argv[4],
      argv[5]
    );

  std::cout <<
    boost::format( "\nProperties for zone (%s,%s,%s) in group \"%s\"\n\n" ) %
      argv[3] %
      argv[4] %
      argv[5] %
      argv[2]; 

  BOOST_FOREACH(
    user::hdf5::zone_properties_hash_entry prop, my_properties_hash
  )
  {

    std::cout << "Property (" <<
                 prop.sName << "," <<
                 prop.sTag1 << "," <<
                 prop.sTag2 << "):  Value: " <<
                 prop.sValue << std::endl;

  }

  std::cout << "\n";

  return EXIT_SUCCESS;

}
