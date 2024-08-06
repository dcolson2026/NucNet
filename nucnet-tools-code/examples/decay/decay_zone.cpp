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
//! \brief Code to decay input abundances over input times.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>
#include <algorithm>

#include <boost/algorithm/string/find.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>

#include "nnt/auxiliary.h"
#include "user/flow_utilities.h"
#include "user/nuclear_decay_utilities.h"

//##############################################################################
// Types
//##############################################################################

typedef std::map<std::string, boost::any> decay_data_t;

namespace po = boost::program_options;

//##############################################################################
// Defines
//##############################################################################

#define S_NUC_XPATH            "nuc_xpath"
#define S_REAC_XPATH           "reac_xpath"
#define S_ZONE_XPATH           "zone_xpath"
#define S_ZONE_FILE            "zone_file"
#define S_TIMES                "times"
#define S_PRINTOUT             "printout"
#define S_DEBUG                "debug"

//##############################################################################
// at_option_parser().
//##############################################################################

std::pair<std::string,std::string> at_option_parser( std::string const& s )
{
  if ( '@' == s[0] )
    return std::make_pair( std::string( "response-file" ), s.substr( 1 ) );
  else
    return std::pair<std::string,std::string>();
}

//##############################################################################
// response_file().
//##############################################################################

po::variables_map
response_file(
  po::variables_map& vmap,
  po::options_description& all
)
{
  // Load the file and tokenize it
  std::ifstream ifs( vmap["response-file"].as<std::string>().c_str() );
  if( !ifs )
  {
    std::cout << "Could not open the response file\n";
    exit( EXIT_FAILURE );
  }

  // Read the whole file into a string
  std::stringstream ss;
  ss << ifs.rdbuf();

  // Split the file content

  std::string sep1("");
  std::string sep2(" \n\r");
  std::string sep3("\"");
  std::string sstr = ss.str();

  boost::escaped_list_separator<char> els( sep1, sep2, sep3 );
  boost::tokenizer< boost::escaped_list_separator<char> > tok( sstr, els );

  std::vector<std::string> args;
  std::copy( tok.begin(), tok.end(), std::back_inserter( args ) );

  // Parse the file and store the options
  store( po::command_line_parser( args ).options( all ).run(), vmap );

  return vmap;
}

//##############################################################################
// get_input().
//##############################################################################

decay_data_t
get_input( int argc, char **argv )
{

  decay_data_t decay_data;
  try
  {
    std::vector<double> v_times;
    po::variables_map vmap;
    std::string
      s_nuc_xpath = "",
      s_reac_xpath = "[count(reactant) = 1]",
      s_zone_xpath = "[last()]";

    po::options_description all( "\nOptions" );
    all.add_options()

    // Help option
    (
      "help",
      "Print out help and exit"
    )

    // Option for value of XPath to select nuclei
    (
      S_NUC_XPATH,
      po::value<std::vector<std::string> >()->multitoken()->composing(),
      "XPath to select nuclei (default: all nuclides)"
    )

    // Option for value of XPath to select reactions
    (
      S_REAC_XPATH,
      po::value<std::vector<std::string> >()->multitoken()->composing(),
      "XPath to select reactions (default: [count(reactant) = 1])"
    )

    // Option for value of XPath to select zone
    (
      S_ZONE_XPATH,
      po::value<std::vector<std::string> >()->multitoken()->composing(),
      "XPath to select zone (default: [last()])"
    )

    // Option for value of XPath to select zone
    (
      S_ZONE_FILE,
      po::value<std::string>(),
      "XPath to select zone file (if zone data not present in input_xml)"
    )

    // Option to specify response file
    ( "response-file", po::value<std::string>(),
      "can be specified with '@name', too" )

    // Option for times
    (  S_TIMES, po::value<std::vector<double> >()->multitoken()->composing(),
      "Times for decay (in s) as space-delimited list" )

    // Option for times
    (  nnt::s_T9, po::value<double>()->default_value( 1.e-5, "1.e-5" ),
      "Value of T9 to supply in case not present in input" )

    // Option for times
    (  nnt::s_RHO, po::value<double>()->default_value( 1.e-5, "1.e-5" ),
      "Value of density to supply in case not present in input" )

    // Option for debug
    (  S_PRINTOUT, po::value<std::string>()->default_value( "no" ),
      "Print out" )

    // Option for debug
    (  S_DEBUG, po::value<std::string>()->default_value( "no" ),
      "Debug" );

    store(
      po::command_line_parser( argc, argv ).
      options( all ).
      extra_parser( at_option_parser ).
      run(),
      vmap
    );

    if( vmap.count( "response-file" ) )
      vmap = response_file( vmap, all );

    if( argc < 3 || vmap.count("help") == 1 )
    {
      std::cerr <<
        "\nUsage: " << argv[0] << " net_xml zone_xml output_xml [options]" <<
        std::endl;
      std::cout << all << "\n";
      exit( EXIT_FAILURE );
    }

    if(
      boost::algorithm::ifind_first( argv[1], "--" ) ||
      boost::algorithm::ifind_first( argv[2], "--" )
    )
    {
      std::cerr << "Required arguments must not be options." << std::endl;
      exit( EXIT_FAILURE );
    }

    //==========================================================================
    // Get data.
    //==========================================================================

    if( vmap.count( S_ZONE_FILE ) )
    {
      decay_data[S_ZONE_FILE] = vmap[S_ZONE_FILE].as<std::string>();
    }
    else
    {
      decay_data[S_ZONE_FILE] = std::string( argv[1] );
    }

    if( vmap.count( S_NUC_XPATH ) )
    {
      s_nuc_xpath = "";
      BOOST_FOREACH(
        std::string s,
        vmap[S_NUC_XPATH].as<std::vector<std::string> >()
      )
      {
        s_nuc_xpath += s + " ";
      }
    }

    decay_data[S_NUC_XPATH] = s_nuc_xpath;

    if( vmap.count( S_REAC_XPATH ) )
    {
      s_reac_xpath = "";
      BOOST_FOREACH(
        std::string s,
        vmap[S_REAC_XPATH].as<std::vector<std::string> >()
      )
      {
        s_reac_xpath += s + " ";
      }
    }

    decay_data[S_REAC_XPATH] = s_reac_xpath;

    if( vmap.count( S_ZONE_XPATH ) )
    {
      s_zone_xpath = "";
      BOOST_FOREACH(
        std::string s,
        vmap[S_ZONE_XPATH].as<std::vector<std::string> >()
      )
      {
        s_zone_xpath += s + " ";
      }
    }

    decay_data[S_ZONE_XPATH] = s_zone_xpath;

    if( !vmap.count( S_TIMES ) )
    {
      std::cerr <<
        "Must have at least one decay time (set with --times option)." <<
        std::endl;
      exit( EXIT_FAILURE );
    }

    BOOST_FOREACH( double d, vmap[S_TIMES].as<std::vector<double> >() )
    {
      if( d <= 0 )
      {
        std::cerr << d << " is an invalid time." << std::endl;
        exit( EXIT_FAILURE );
      }
      v_times.push_back( d );
    }

    std::sort( v_times.begin(), v_times.end() );

    decay_data[S_TIMES] = v_times;

    decay_data[nnt::s_T9] = vmap[nnt::s_T9].as<double>();

    decay_data[nnt::s_RHO] = vmap[nnt::s_RHO].as<double>();

    decay_data[S_PRINTOUT] = vmap[S_PRINTOUT].as<std::string>();

    decay_data[S_DEBUG] = vmap[S_DEBUG].as<std::string>();

    return decay_data;

  }
  catch( std::exception& e )
  {
    std::cerr << "Error: " << e.what() << "\n";
    exit( EXIT_FAILURE );
  }
  catch(...)
  {
    std::cerr << "Exception of unknown type!\n";
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  decay_data_t decay_data;

  Libnucnet * p_my_nucnet, * p_my_output;
  nnt::Zone zone;

  //============================================================================
  // Get data.
  //============================================================================

  decay_data = get_input( argc, argv );

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Net__updateFromXml(
    Libnucnet__getNet( p_my_nucnet ),
    argv[1],
    boost::any_cast<std::string>( decay_data[S_NUC_XPATH] ).c_str(),
    boost::any_cast<std::string>( decay_data[S_REAC_XPATH] ).c_str()
  );

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    boost::any_cast<std::string>( decay_data[S_ZONE_FILE] ).c_str(),
    boost::any_cast<std::string>( decay_data[S_ZONE_XPATH] ).c_str()
  );

  //============================================================================
  // Get zone.
  //============================================================================

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  if( zone_list.size() != 1 )
  {
    std::cerr << "Invalid number of zones." << std::endl;
    exit( EXIT_FAILURE );
  }

  zone = *(zone_list.begin());

  //============================================================================
  // Create output.
  //============================================================================

  p_my_output = nnt::create_network_copy( p_my_nucnet );

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  //============================================================================
  // Decay abundances.
  //============================================================================

  size_t k = 0;

  double d_t_prev = 0;

  BOOST_FOREACH(
    double d_t, 
    boost::any_cast<std::vector<double> >( decay_data[S_TIMES] )
  )
  {

    double d_dt = d_t - d_t_prev;

    zone.updateProperty( nnt::s_TIME, d_t );
    zone.updateProperty( nnt::s_DTIME, d_dt );

    if( !zone.hasProperty( nnt::s_T9 ) )
    {
      zone.updateProperty(
        nnt::s_T9,
        boost::any_cast<double>( decay_data[nnt::s_T9] )
      );
    }

    if( !zone.hasProperty( nnt::s_RHO ) )
    {
      zone.updateProperty(
        nnt::s_RHO,
        boost::any_cast<double>( decay_data[nnt::s_RHO] )
      );
    }

    user::decay_abundances(
      p_my_nucnet,
      d_dt,
      "",
      boost::any_cast<std::string>( decay_data[S_DEBUG] )
    );

    Libnucnet__relabelZone(
      p_my_nucnet,
      zone.getNucnetZone(),
      ( boost::lexical_cast<std::string>( ++k ) ).c_str(),
      NULL,
      NULL
    );

    nnt::write_xml( p_my_output, zone.getNucnetZone() );
    if( boost::any_cast<std::string> (decay_data[S_PRINTOUT] ) == "yes" )
    {
      nnt::print_zone_abundances( zone );
    }

    d_t_prev = d_t;

  }

  //============================================================================
  // Output.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_output,
    "%.15e"
  );

  Libnucnet__writeToXmlFile( p_my_output, argv[2] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_output );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
