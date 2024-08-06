//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer and Michael J.
// Bojazi.
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
//!
//! \file
//! \brief Code for program options.
//!
////////////////////////////////////////////////////////////////////////////////

#include "program_options.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */
namespace user
{

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
// compose_option_string_from_vector().
//##############################################################################

std::string
compose_option_string_from_vector(
  const po::variables_map& vmap,
  const char * s_option
)
{
  std::string result( "" );
  if( vmap.count( s_option ) )
  {
    BOOST_FOREACH(
      const std::string s,
      vmap[s_option].as<std::vector<std::string> >()
    )
    {
      result += s + " ";
    }
    boost::trim( result );
  }
  return result;
}

//##############################################################################
// compose_option_vector_of_vectors().
//##############################################################################

std::vector<std::vector<std::string> >
compose_option_vector_of_vectors(
  const po::variables_map& vmap,
  const char * s_option,
  std::string s_sep1,
  std::string s_sep2,
  size_t i_expected
)
{

  std::vector<std::vector<std::string> > result;

  std::stringstream ss;
  std::vector<std::string> sv = vmap[s_option].as<std::vector<std::string> >();
  std::copy(
    sv.begin(),
    sv.end(),
    std::ostream_iterator<std::string>( ss )
  );

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep1( s_sep1.c_str() );
  std::string s_s = ss.str();
  tokenizer tokens1( s_s, sep1 );

  std::vector<std::string> v_s;
  std::copy( tokens1.begin(), tokens1.end(), std::back_inserter( v_s ) );

  BOOST_FOREACH( std::string s, v_s )
  {

    boost::char_separator<char> sep2( s_sep2.c_str() );
    tokenizer tokens2( s, sep2 );

    std::vector<std::string> v_x;
    std::copy(
      tokens2.begin(), tokens2.end(), std::back_inserter( v_x )
    );

    if( v_x.size() != i_expected )
    {
      std::cerr << "Invalid input for " << s << std::endl;
      exit( EXIT_FAILURE );
    }

    for( size_t i = 0; i < v_x.size(); i++ )
    {
      boost::trim( v_x[i] );
    }

    result.push_back( v_x );

  }

  return result;

}

std::vector<std::vector<std::string> >
compose_option_vector_of_vectors(
  const po::variables_map& vmap,
  const char * s_option,
  size_t i_expected
)
{
  return
    compose_option_vector_of_vectors( vmap, s_option, "{}", ",", i_expected );
}

std::vector<std::string>
compose_option_vector_of_strings(
  const po::variables_map& vmap,
  const char * s_option,
  std::string s_sep
)
{

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep( s_sep.c_str() );
  std::string s_s = vmap[s_option].as<std::string>();
  tokenizer tokens( s_s, sep );

  std::vector<std::string> v_s;
  for(
    tokenizer::iterator it = tokens.begin();
    it != tokens.end();
    it++
  )
  {
    std::string s = *it;
    boost::algorithm::trim( s );
    v_s.push_back( s );
  }
  return v_s;

}

} // namespace user
