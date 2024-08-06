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
//! \brief Example code to print out properties and mass fractions into a
//!        values delimited file.
////////////////////////////////////////////////////////////////////////////////

/*##############################################################################
// Includes.
//############################################################################*/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <boost/any.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#define S_CUTOFF             "cutoff"
#define S_DELIM              "delimiter"
#define S_EXAMPLE            "example"
#define S_FORMAT             "x_format"
#define S_NUC_XPATH          "nuc_xpath"
#define S_ZONE_LABEL         "zone_label"
#define S_PROPERTY           "property"
#define S_ZONE_XPATH         "zone_xpath"

typedef std::map<std::string, boost::any> my_map_t;

namespace po = boost::program_options;

//##############################################################################
// get_input().
//##############################################################################

my_map_t
get_input( int argc, char **argv )
{

  my_map_t param_map;
  try
  {
    po::variables_map vmap;
    std::string s_nuc_xpath = "none", s_zone_xpath = "";

    po::options_description help( "\nHelp Options" );
    help.add_options()
    // Option to print help message for usage statement
    ( "help", "print out usage statement and exit" )

    // Option to print example usage
    ( S_EXAMPLE, "print out example usage and exit" )

    // Option for value of XPath to select nuclei
    (
      S_NUC_XPATH,
      po::value<std::vector<std::string> >()->multitoken()->composing(),
      "XPath to select nuclei (default: no nuclides)"
    )

    // Option for value of XPath to select nuclei
    (
      S_ZONE_XPATH,
      po::value<std::vector<std::string> >()->multitoken()->composing(),
      "XPath to select zones (default: all zones)"
    )

    // Option for printing zone labels
    (
      S_ZONE_LABEL,
      po::value<std::vector<unsigned int> >()->composing(),
      "Zone label to print out"
    )

    // Option to select output delimiter
    (
      S_DELIM,
      po::value<std::string>()->default_value( " " ),
      "Output delimiter"
    )

    // Option to select output delimiter
    (
      S_CUTOFF,
      po::value<double>(),
      "Mass fraction printout cutoff (default: not set)"
    )

    // Option to select mass fraction format
    (
      S_FORMAT,
      po::value<std::string>()->default_value( "%g" ),
      "Mass fraction output format"
    )

    // Option for property
    (
      S_PROPERTY,
      po::value<std::vector<std::string> >()->composing(),
      "Property (enter for each property)"
    )
    ;

    store(
      po::command_line_parser( argc, argv ).
      options( help ).
      run(),
      vmap
    );

    if( vmap.count( S_EXAMPLE ) )
    {
      std::cout << "\n" <<
        argv[0] << " out.xml --" << S_PROPERTY << " rho " <<
        "--" << S_NUC_XPATH << " \"[z = 2]\" " <<
        "--" << S_PROPERTY << " t9  --" << S_ZONE_LABEL << " 1 " <<
        "--" << S_CUTOFF << " 1.e-25 " <<
        "--" << S_ZONE_XPATH << " \"[position() >= last() - 10]\"\n" <<
        std::endl;
      exit( EXIT_FAILURE );
    }

    if( argc == 1 || vmap.count( "help" ) == 1 )
    {
      std::cout << std::endl << "Usage: " << argv[0]
                << " xml [options]"
                << std::endl;
      std::cout << "\n  xml = input xml file" << std::endl;
      std::cout << help << std::endl;
      exit( EXIT_FAILURE );
    }

    //==========================================================================
    // XPath strings.
    //==========================================================================

    if( vmap.count( S_NUC_XPATH ) )
    {
      BOOST_FOREACH(
        std::string s,
        vmap[S_NUC_XPATH].as<std::vector<std::string> >()
      )
      {
        s_nuc_xpath = s;
      }
    }

    param_map[S_NUC_XPATH] = s_nuc_xpath;

    if( vmap.count( S_ZONE_XPATH ) )
    {
      BOOST_FOREACH(
        std::string s,
        vmap[S_ZONE_XPATH].as<std::vector<std::string> >()
      )
      {
        s_zone_xpath += s + " ";
      }
    }

    param_map[S_ZONE_XPATH] = s_zone_xpath;

    //==========================================================================
    // Properties.
    //==========================================================================

    if( vmap.count( S_PROPERTY ) )
    {
      param_map[S_PROPERTY] = vmap[S_PROPERTY].as<std::vector<std::string> >();
    }

    //==========================================================================
    // Other input.
    //==========================================================================

    param_map[S_DELIM] = vmap[S_DELIM].as<std::string>();

    param_map[S_FORMAT] = vmap[S_FORMAT].as<std::string>();

    if( vmap.count( S_CUTOFF ) )
    {
      param_map[S_CUTOFF] = vmap[S_CUTOFF].as<double>();
    }

    if(
      vmap.count(S_ZONE_LABEL)
    )
    {
      param_map[S_ZONE_LABEL] =
        vmap[S_ZONE_LABEL].as<std::vector<unsigned int> >();
    }

    return param_map;

  }
  catch( std::exception& e )
  {
    std::cerr << "error: " << e.what() << "\n";
    exit( EXIT_FAILURE );
  }
  catch(...)
  {
    std::cerr << "Exception of unknown type!\n";
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// filter_species().
//##############################################################################

void
filter_species(
  Libnucnet * p_nucnet,
  Libnucnet__NucView * p_view,
  boost::unordered_set<std::string>& species_set,
  double d_cutoff
)
{

  nnt::species_list_t species_list =
    nnt::make_species_list( Libnucnet__NucView__getNuc( p_view ) );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    species_set.insert(
      Libnucnet__Species__getName( species.getNucnetSpecies() )
    );
  }

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_nucnet ) )
  {
    Libnucnet__Zone * p_zone = zone.getNucnetZone();

    BOOST_FOREACH( nnt::Species species, species_list )
    {
      Libnucnet__Species * p_species = species.getNucnetSpecies();

      double d_x =
        Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species) *
          Libnucnet__Species__getA( p_species );

      if( d_x > d_cutoff )
      {
        if(
          species_set.find( Libnucnet__Species__getName( p_species ) )
          !=
          species_set.end()
        )
        {
          species_set.erase( Libnucnet__Species__getName( p_species ) );
        }
      }

    }
  }
}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  my_map_t param_map;

  Libnucnet *p_my_nucnet;
  Libnucnet__NucView *p_view;
  std::vector<std::vector<std::string> > v_property;
  std::vector<Libnucnet__Species *> v_species;
  std::vector<unsigned int> v_labels;
  boost::unordered_set<std::string> species_set;

  param_map = get_input( argc, argv );

  if( param_map.find( S_PROPERTY ) != param_map.end() )
  {
    BOOST_FOREACH(
      std::string s,
      boost::any_cast<std::vector<std::string> >( param_map[S_PROPERTY] )
    )
    {

      std::vector<std::string> prop;
      boost::char_separator<char> sep(",");
      tokenizer tok( s, sep );
      int j = 0;
      for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
      {
        std::string s_tmp = *it;
        boost::algorithm::trim( s_tmp );
        prop.push_back( s_tmp );
        j++;
      }

      if( j > 3 )
      {
        std::cerr << "Invalid property and tags." << std::endl;
        return EXIT_FAILURE;
      }

      v_property.push_back( prop );
    }
  }

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    argv[1],
    boost::any_cast<std::string>( param_map[S_ZONE_XPATH] ).c_str()
  );

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  Libnucnet__Nuc__setSpeciesCompareFunction(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    (Libnucnet__Species__compare_function)
      nnt::species_sort_by_z_then_a
  );

  Libnucnet__Nuc__sortSpecies(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
  );

  p_view =
    Libnucnet__NucView__new( 
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      boost::any_cast<std::string>( param_map[S_NUC_XPATH] ).c_str()
    );

  if( param_map.find( S_CUTOFF ) != param_map.end() )
  {
    double d_cutoff = boost::any_cast<double>( param_map[S_CUTOFF] );
    filter_species( p_my_nucnet, p_view, species_set, d_cutoff );
  }

  BOOST_FOREACH(
    nnt::Species species,
    nnt::make_species_list(
      Libnucnet__NucView__getNuc( p_view )
    )
  )
  {
    if(
      species_set.find(
        Libnucnet__Species__getName( species.getNucnetSpecies() )
      ) == species_set.end()
    )
    {
      v_species.push_back( species.getNucnetSpecies() );
    }
  }

  if( param_map.find( S_ZONE_LABEL ) != param_map.end() )
  {
    v_labels =
      boost::any_cast<std::vector<unsigned int> >( param_map[S_ZONE_LABEL] );

    std::sort( v_labels.begin(), v_labels.end() );

    BOOST_FOREACH( unsigned int u, v_labels )
    {
      std::cout << "Label " << u;
      if( u != v_labels.back() )
      {
        std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
      }
    }
    if( v_property.size() > 0 || v_species.size() > 0 )
    {
      std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
    }
  }
      
  BOOST_FOREACH( std::vector<std::string> v, v_property )
  {
    if( v.size() == 1 )
    { 
      std::cout << v[0];
    }
    else
    {
      std::cout << "\"";
      for( size_t i = 0; i < v.size(); i++ )
      {
        std::cout << v[i];
        if( v[i] != v.back() )
        {
          std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
        }
      }
      std::cout << "\"";
    }        
    if( v != v_property.back() )
    {
      std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
    }
    else
    {
      if( v_species.size() > 0 )
      {
        std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
      }
    }
  }

  BOOST_FOREACH( Libnucnet__Species * p, v_species )
  {
    std::cout << Libnucnet__Species__getName( p );
    if( p != v_species.back() )
    { 
      std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
    }
  }

  std::cout << std::endl;

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  boost::format
    fmt( boost::any_cast<std::string>( param_map[S_FORMAT] ).c_str() );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    if( v_labels.size() > 0 )
    {
      BOOST_FOREACH( unsigned int u, v_labels )
      {
        std::cout <<
          Libnucnet__Zone__getLabel( zone.getNucnetZone(), u );
        if( u != v_labels.back() )
        {
          std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
        }
      }
      if( v_property.size() > 0 || v_species.size() > 0 )
      {
        std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
      }
    }
      
    BOOST_FOREACH( std::vector<std::string> v, v_property )
    {
      std::string s;
      switch( v.size() )
      {
        case 1:
          s = zone.getProperty<std::string>( v[0] );
          break;
        case 2:
          s = zone.getProperty<std::string>( v[0], v[1] );
          break;
        case 3:
          s = zone.getProperty<std::string>( v[0], v[1], v[2] );
          break;
        default:
          break;
      }
      std::cout << s;
      if( v != v_property.back() )
      {
        std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
      }
      else
      {
        if( v_species.size() > 0 )
        {
          std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
        }
      }
    }

    BOOST_FOREACH( Libnucnet__Species * p, v_species )
    {
      fmt %
        (
          Libnucnet__Zone__getSpeciesAbundance(
            zone.getNucnetZone(),
            p
          ) * Libnucnet__Species__getA( p )
        );
      std::cout << fmt.str();
      if( p != v_species.back() )
      {
        std::cout << boost::any_cast<std::string>( param_map[S_DELIM] );
      }
    }
          
    std::cout << std::endl;

  }

  Libnucnet__NucView__free( p_view );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
