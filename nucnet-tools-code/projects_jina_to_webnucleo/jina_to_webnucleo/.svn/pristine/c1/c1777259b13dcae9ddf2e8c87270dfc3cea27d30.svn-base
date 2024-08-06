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
//! \brief Code to convert reaclib file to webnucleo reaction xml.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <Libnucnet.h>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <nnt/iter.h>

//##############################################################################
// Defines.
//##############################################################################

#define D_SPINF    0.
#define D_SPINT    0.
#define D_TLOWHF   0.01
#define D_TLOWFIT  0.01
#define D_THIGHFIT 10.
#define D_ACC      0.

//##############################################################################
// Types.
//##############################################################################

struct reaction_data
{
  std::string source;
  std::string flag;
  double * a;

  reaction_data(
    const std::string &s,
    const std::string &f,
    double * p_a
  ) : source( s ), flag( f ), a( p_a ) {}

};

//##############################################################################
// create_reaction().
//##############################################################################

Libnucnet__Reaction *
create_reaction(
  Libnucnet__Nuc * p_nuc,
  std::map<std::string, std::string>& replacements,
  std::set< std::string >& exclusions,
  std::string chapter,
  std::string flag,
  std::vector<std::string>& elements
)
{

  int i_chapter = boost::lexical_cast<int>( chapter );
  size_t i_reac, i_total;

  std::vector<std::string> reactants;
  std::vector<std::string> products;

  for( size_t i = 0; i < 6; i++ )
  {
    if( exclusions.find( elements[i] ) != exclusions.end() )
      return NULL;
  }

  for( size_t i = 0; i < 6; i++ )
  {
    std::map<std::string,std::string>::iterator it =
      replacements.find( elements[i] );
    if( it != replacements.end() )
    {
        elements[i] = it->second;
    }
  }

  switch( i_chapter )
  {

    case 1:
      i_reac = 1;
      i_total = 2;
      break;
    case 2:
      i_reac = 1;
      i_total = 3;
      break;
    case 3:
      i_reac = 1;
      i_total = 4;
      break;
    case 4:
      i_reac = 2;
      i_total = 3;
      break;
    case 5:
      i_reac = 2;
      i_total = 4;
      break;
    case 6:
      i_reac = 2;
      i_total = 5;
      break;
    case 7:
      i_reac = 2;
      i_total = 6;
      break;
    case 8:
      i_reac = 3;
      i_total = 4;
      break;
    case 9:
      i_reac = 3;
      i_total = 5;
      break;
    case 10:
      i_reac = 4;
      i_total = 6;
      break;
    case 11:
      i_reac = 1;
      i_total = 5;
      break;
    default:
      std::cerr << "No such chapter." << std::endl;
      exit( EXIT_FAILURE );
  } 

  for( size_t i = 0; i < i_total; i++ )
  {
    if( i < i_reac )
      reactants.push_back( elements[i] );
    else
      products.push_back( elements[i] );
  }

  if( flag != "w" && flag != "s" )
  {
    switch( i_chapter )
    {

      case 1:
        reactants.push_back( "gamma" );
	break;
      case 2:
        reactants.push_back( "gamma" );
	break;
      case 3:
        reactants.push_back( "gamma" );
	break;
      case 4:
        products.push_back( "gamma" );
	break;
      case 8:
        products.push_back( "gamma" );
	break;
      case 11:
        reactants.push_back( "gamma" );
	break;
      default:
        break;

    }

  }
  else if( flag == "w" )
  {
    unsigned int z_reac = 0;
    unsigned int z_prod = 0;
    for( size_t i = 0; i < reactants.size(); i++ )
    {
      z_reac +=
        Libnucnet__Species__getZ(
          Libnucnet__Nuc__getSpeciesByName( p_nuc, reactants[i].c_str() )
        );
    }
    for( size_t i = 0; i < products.size(); i++ )
    {
      z_prod +=
        Libnucnet__Species__getZ(
          Libnucnet__Nuc__getSpeciesByName( p_nuc, products[i].c_str() )
        );
    }
    
    if( z_reac > z_prod )
    {
      products.push_back( "positron" );
      products.push_back( "neutrino_e" );
    }
    else
    {
      products.push_back( "electron" );
      products.push_back( "anti-neutrino_e" );
    }
  } 

  Libnucnet__Reaction * p_reaction = Libnucnet__Reaction__new();

  for( size_t i = 0; i < reactants.size(); i++ )
  {
    Libnucnet__Reaction__addReactant( p_reaction, reactants[i].c_str() );
  }

  for( size_t i = 0; i < products.size(); i++ )
  {
    Libnucnet__Reaction__addProduct( p_reaction, products[i].c_str() );
  }

  return p_reaction;

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char * argv[] )
{

  Libnucnet__Nuc * p_nuc;
  Libnucnet__Reac * p_my_reactions;
  double aa[8] = {0.};
  std::ifstream input_file;
  std::vector<std::string> elements(7);
  std::string chapter, line, source, flag, sv, sfit;
  std::map<std::string, std::string> replacements;
  std::set<std::string> exclusions;
  std::multimap<std::string, reaction_data> my_reaction_data;
  std::multimap<std::string,reaction_data>::iterator it;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {    std::cout << std::endl;
    std::cout << argv[0] << " ../nucnet-tools-code/data_pub/jina_to_webnucleo/reaction_data/20140918default2 reaction_data.xml" <<
      std::endl << std::endl;
    return EXIT_FAILURE;
  }

  if( argc != 3 )
  {

    std::cerr << std::endl;
    std::cerr << "Usage: " << argv[0] << " input_text_file output_xml_file " <<
      std::endl << std::endl;
    std::cerr << "  input_text_file = input reaclib reaction data file" <<
      std::endl << std::endl;
    std::cerr << "  output_xml_file = output xml reaction data file" <<
      std::endl << std::endl;
    std::cerr << "For an example usage, type " << std::endl << std::endl;
    std::cerr << argv[0] << " --example" << std::endl << std::endl;

    return EXIT_FAILURE;

  }

  //============================================================================
  // Create nuclide collection.
  //============================================================================

  p_nuc = Libnucnet__Nuc__new();

  for( unsigned int iz = 1; iz <= 200; iz++ )
  {
    for( unsigned int ia = 1; ia <= 500; ia++ )
    {
      Libnucnet__Nuc__addSpecies(
        p_nuc,
        Libnucnet__Species__new( iz, ia, "", 0, "", 0., 0., NULL, NULL )
      );
    }
  } 

  //============================================================================
  // Special species.
  //============================================================================

  Libnucnet__Nuc__addSpecies(
    p_nuc,
    Libnucnet__Species__new( 0, 1, "", 0, "", 0., 0., NULL, NULL )
  );

  //============================================================================
  // Replacements.
  //============================================================================

  replacements.insert( std::make_pair( "p", "h1" ) );
  replacements.insert( std::make_pair( "d", "h2" ) );
  replacements.insert( std::make_pair( "t", "h3" ) );

  //============================================================================
  // Exclusions.  Reactions involving these species will not be included.
  //============================================================================

  exclusions.insert( "al-6" );
  exclusions.insert( "al*6" );

  //============================================================================
  // Create reaction collection.
  //============================================================================

  p_my_reactions = Libnucnet__Reac__new();

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
  // Loop on the data.
  //============================================================================

  while( std::getline( input_file, line ) )
  {

    chapter = line;
    
    boost::trim( chapter );

    std::getline( input_file, line );

    for( size_t i = 0; i < 6; i++ )
    {
      elements[i] = line.substr( 5+i*5, 5 );
      boost::trim( elements[i] );
    }

    source = line.substr( 43, 4 );

    flag = line.substr( 47, 1 );

    sv = line.substr( 48, 1 );

    std::getline( input_file, line );

    std::istringstream is1( line );

    is1 >> aa[0] >> aa[1] >> aa[2] >> aa[3];

    std::getline( input_file, line );

    std::istringstream is2( line );

    is2 >> aa[4] >> aa[5] >> aa[6];

    Libnucnet__Reaction * p_reaction =
      create_reaction(
        p_nuc,
        replacements,
        exclusions,
        chapter,
        flag,
        elements
      );

    if( p_reaction )
    {

      double * a = (double *) calloc( 8, sizeof( double ) );

      memcpy( a, aa, sizeof( double ) * 8 );

      Libnucnet__Reaction * p_old_reaction =
	Libnucnet__Reac__getReactionByString(
	  p_my_reactions,
	  Libnucnet__Reaction__getString( p_reaction )
	);

      if( !p_old_reaction )
      {
	Libnucnet__Reac__addReaction( p_my_reactions, p_reaction );
      }
      else
      {
	Libnucnet__Reaction__free( p_reaction );
	p_reaction = p_old_reaction;
      }
	
      my_reaction_data.insert(
	std::pair<std::string,reaction_data>(
	  Libnucnet__Reaction__getString( p_reaction ),
	  reaction_data( source, flag, a )
	)
      );

    }
      
  }

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( p_my_reactions );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    size_t i_count =
      my_reaction_data.count(
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
      );

    if( i_count == 1 )
    {

      it =
        my_reaction_data.find(
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
        );
 
      if(
        it->second.a[1] == 0 &&
        it->second.a[2] == 0 &&
        it->second.a[3] == 0 &&
        it->second.a[4] == 0 &&
        it->second.a[5] == 0 &&
        it->second.a[6] == 0
      )
      {
        Libnucnet__Reaction__updateSingleRate(
          reaction.getNucnetReaction(),
          exp( it->second.a[0] )
        );
      }
      else
      {

        if( it->second.flag == " " )
          sfit = "0";
        else
          sfit = it->second.flag;

        Libnucnet__Reaction__addNonSmokerFit(
          reaction.getNucnetReaction(),
          sfit.c_str(),
          it->second.a,
          D_SPINT,
          D_SPINF,
          D_TLOWHF,
          D_TLOWFIT,
          D_THIGHFIT,
          D_ACC
        );
      }

      Libnucnet__Reaction__updateSource(
        reaction.getNucnetReaction(),
        it->second.source.c_str()
      );

    }
    else
    {

      std::pair<
        std::multimap<std::string,reaction_data>::iterator,
        std::multimap<std::string,reaction_data>::iterator
      > ii =
          my_reaction_data.equal_range(
            Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
          );

      int i_fit = 0;
      int i_wfit = 0;
      int i_nfit = 0;
      int i_rfit = 0;
      int i_sfit = 0;

      for( it = ii.first; it != ii.second; it++ )
      {

        if( it->second.flag == "w" )
          sfit = "w" + boost::lexical_cast<std::string>( i_wfit++ );
        else if( it->second.flag == "n" )
          sfit = "n" + boost::lexical_cast<std::string>( i_nfit++ );
        else if( it->second.flag == "r" )
          sfit = "r" + boost::lexical_cast<std::string>( i_rfit++ );
        else if( it->second.flag == "s" )
          sfit = "s" + boost::lexical_cast<std::string>( i_sfit++ );
        else
          sfit = boost::lexical_cast<std::string>( i_fit++ );

        Libnucnet__Reaction__addNonSmokerFit(
          reaction.getNucnetReaction(),
          sfit.c_str(),
          it->second.a,
          D_SPINT,
          D_SPINF,
          D_TLOWHF,
          D_TLOWFIT,
          D_THIGHFIT,
          D_ACC
        );

        Libnucnet__Reaction__updateSource(
          reaction.getNucnetReaction(),
          it->second.source.c_str()
        );

      }

    }

  }

  input_file.close();

  //============================================================================
  // Output.
  //============================================================================

  Libnucnet__Reac__writeToXmlFile( p_my_reactions, argv[2] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Nuc__free( p_nuc );
  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

