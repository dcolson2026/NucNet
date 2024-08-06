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
//! \brief Code to convert winve file to webnucleo nuclear xml.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <fstream>
#include <Libnucnet.h>

#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

//##############################################################################
// Defines.
//##############################################################################

#define I_VEC  24

//##############################################################################
// get_items_from_line().
//##############################################################################

std::vector<std::string>
get_items_from_line( std::string& line )
{
  
  typedef boost::tokenizer<boost::char_separator<char> > my_tok;
  boost::char_separator<char> sep( " " );

  my_tok tok( line, sep );

  std::vector<std::string> items;

  for( my_tok::iterator it = tok.begin(); it != tok.end();++it )
  {
     items.push_back( *it );
  }

  return items;

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char * argv[] )
{

  Libnucnet__Nuc * p_my_nuclei;
  Libnucnet__Species * p_species;
  std::ifstream input_file;
  std::string line;
  gsl_vector * p_t9, * p_partf;
  std::vector<std::string> items;
  std::vector<double> t9 =
    boost::assign::list_of(0.1)(0.15)(0.2)(0.3)(0.4)(0.5)
                          (0.6)(0.7)(0.8)(0.9)(1.0)(1.5)
                          (2.0)(2.5)(3.0)(3.5)(4.0)(4.5)
                          (5.0)(6.0)(7.0)(8.0)(9.0)(10.0);

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " ../nucnet-tools-code/data_pub/jina_to_webnucleo/nuclide_data/winvne_v2.0.dat nuclides.xml" <<
      std::endl << std::endl;
    return EXIT_FAILURE;
  }

  if( argc != 3 )
  {

    std::cerr << std::endl;
    std::cerr << "Usage: " << argv[0] << " input_text_file output_xml_file " <<
      std::endl << std::endl;
    std::cerr << "  input_text_file = input winve nuclear data file" <<
      std::endl << std::endl;
    std::cerr << "  output_xml_file = output xml nuclear data file" <<
      std::endl << std::endl;
    std::cerr << "For an example usage, type " << std::endl << std::endl;
    std::cerr << argv[0] << " --example" << std::endl << std::endl;

    return EXIT_FAILURE;

  }

  //============================================================================
  // Create nuclear collection.
  //============================================================================

  p_my_nuclei = Libnucnet__Nuc__new();

  //============================================================================
  // Allocate and initialize vector.
  //============================================================================

  p_t9 = gsl_vector_alloc( I_VEC );

  for( size_t i = 0; i < t9.size(); i++ )
    gsl_vector_set( p_t9, i, t9[i] );

  p_partf = gsl_vector_alloc( I_VEC );

  //============================================================================
  // Open file.
  //============================================================================

  input_file.open( argv[1] );

  if( !input_file.good() )
  {
    std::cerr << "Couldn't open input file " << argv[1] << "." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read the first two lines.
  //============================================================================

  std::getline( input_file, line );

  std::getline( input_file, line );

  //============================================================================
  // Loop on the data.
  //============================================================================

  while( std::getline( input_file, line ) )
  {
    if( line.length() > 5 )
    {

      items = get_items_from_line( line );

      p_species =
        Libnucnet__Species__new(
          boost::lexical_cast<unsigned int>( items[2] ),
          boost::lexical_cast<unsigned int>( items[2] ) +
            boost::lexical_cast<unsigned int>( items[3] ),
          items[6].c_str(),
          0,
          "",
          boost::lexical_cast<double>( items[5] ),
          boost::lexical_cast<double>( items[4] ),
          NULL,
          NULL
        );

      for( size_t i = 0; i < 3; i++ )
      {

        std::getline( input_file, line );
        items = get_items_from_line( line );
        for( size_t j = 0; j < items.size(); j++ )
        {
           gsl_vector_set(
             p_partf,
             i * 8 + j,
             log10( boost::lexical_cast<double>( items[j] ) )
           );
        }

      }

      Libnucnet__Species__updatePartitionFunctionData(
        p_species,
        p_t9,
        p_partf
      );
           
      Libnucnet__Nuc__addSpecies( p_my_nuclei, p_species );

    }

  }

  input_file.close();

  //============================================================================
  // Output.
  //============================================================================

  Libnucnet__Nuc__writeToXmlFile( p_my_nuclei, argv[2] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  gsl_vector_free( p_t9 );
  gsl_vector_free( p_partf );
  Libnucnet__Nuc__free( p_my_nuclei );

  return EXIT_SUCCESS;

}

