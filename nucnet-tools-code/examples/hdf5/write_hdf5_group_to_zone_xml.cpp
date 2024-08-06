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
//! \brief Example code to write a given group to xml.
////////////////////////////////////////////////////////////////////////////////

/*##############################################################################
// Includes.
//############################################################################*/

#include <iostream>
#include <string>

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
    std::cout << argv[0] << " input.h5 \"Step 00001\" out.xml" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 4 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " write selected group to xml. ";
    fprintf(
      stderr,
      "\nUsage: %s hdf5_file group out_xml\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  hdf5_file = input hdf5 file\n\n"
    );
    fprintf(
      stderr,
      "  group = group to write to xml\n\n"
    );
    fprintf(
      stderr,
      "  out_xml = name of output xml\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// create_network().
//##############################################################################

Libnucnet *
create_network( char ** argv )
{

  Libnucnet * p_nucnet;

  user::hdf5::nuclide_map my_nuclides;

  p_nucnet = Libnucnet__new();

  my_nuclides = user::hdf5::get_network( argv[1] );

  std::vector<user::hdf5::nuclide> my_nuclide_vector( my_nuclides.size() );

  BOOST_FOREACH( user::hdf5::nuclide_map_entry entry, my_nuclides )
  {
    my_nuclide_vector[(size_t) entry.second.iIndex] = entry.second;
  }

  for( size_t i = 0; i < my_nuclide_vector.size(); i++ )
  {
    int i_state = 0;
    if( strcmp( my_nuclide_vector[i].sState, "" ) )
    {
      i_state = 1;
    }
    Libnucnet__Nuc__addSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
      Libnucnet__Species__new(
        my_nuclide_vector[i].iZ,
        my_nuclide_vector[i].iA,
        my_nuclide_vector[i].sSource,
        i_state,
        my_nuclide_vector[i].sState,
        my_nuclide_vector[i].dMassExcess,
        my_nuclide_vector[i].dSpin,
        NULL,
        NULL
      )
    );
  }

  return p_nucnet;

}
   
//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet * p_nucnet;
  user::hdf5::zone_labels_bimap my_zone_labels;

  check_input( argc, argv );

  p_nucnet = create_network( argv ); 

  my_zone_labels = user::hdf5::get_zone_labels( argv[1], argv[2] );

  user::hdf5::mass_fraction_array_t x =
    user::hdf5::get_group_mass_fractions( argv[1], argv[2] );

  for( size_t i = 0; i < my_zone_labels.size(); i++ )
  {

    user::hdf5::zone_labels_bimap::index<user::to>::type::iterator it =
      my_zone_labels.get<user::to>().find( i );

    Libnucnet__Zone * p_zone =
      Libnucnet__Zone__new(
        Libnucnet__getNet( p_nucnet ),
        it->first.tuple.get<0>().c_str(),
        it->first.tuple.get<1>().c_str(),
        it->first.tuple.get<2>().c_str()
      );

    Libnucnet__addZone( p_nucnet, p_zone );

    user::hdf5::zone_properties_hash my_zone_properties =
      user::hdf5::get_zone_properties(
        argv[1],
        argv[2],
        it->first.tuple.get<0>().c_str(),
        it->first.tuple.get<1>().c_str(),
        it->first.tuple.get<2>().c_str()
      );

    BOOST_FOREACH(
      user::hdf5::zone_properties_hash_entry prop, my_zone_properties
    )
    {
      if( prop.sTag1 == "0" && prop.sTag2 == "0" )
      {
        Libnucnet__Zone__updateProperty(
          p_zone,
          prop.sName.c_str(),
          NULL,
          NULL,
          prop.sValue.c_str()
        );
      }
      else if( prop.sTag1 != "0" && prop.sTag2 == "0" )
      {
        Libnucnet__Zone__updateProperty(
          p_zone,
          prop.sName.c_str(),
          prop.sTag1.c_str(),
          NULL,
          prop.sValue.c_str()
        );
      }
      else
      {
        Libnucnet__Zone__updateProperty(
          p_zone,
          prop.sName.c_str(),
          prop.sTag1.c_str(),
          prop.sTag2.c_str(),
          prop.sValue.c_str()
        );
      }
    }

    BOOST_FOREACH(
      nnt::Species species,
      nnt::make_species_list(
        Libnucnet__Net__getNuc(
          Libnucnet__getNet( p_nucnet )
        )
      )
    )
    {
      Libnucnet__Zone__updateSpeciesAbundance(
        p_zone,
        species.getNucnetSpecies(),
        x[i][Libnucnet__Species__getIndex( species.getNucnetSpecies() )] /
           Libnucnet__Species__getA( species.getNucnetSpecies() )
      ); 
    }

  }

  Libnucnet__updateZoneXmlMassFractionFormat( p_nucnet, "%.15e" );

  Libnucnet__writeZoneDataToXmlFile( p_nucnet, argv[3] );

  Libnucnet__free( p_nucnet );

  return EXIT_SUCCESS;

}
