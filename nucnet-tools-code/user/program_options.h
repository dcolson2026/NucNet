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
//! \brief A header file defining program option utilities.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef PROGRAM_OPTIONS_H
#define PROGRAM_OPTIONS_H

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

namespace user
{

std::pair<std::string,std::string>
at_option_parser( std::string const& );

po::variables_map
response_file(
  po::variables_map&,
  po::options_description&
);

std::string
compose_option_string_from_vector(
  const po::variables_map&,
  const char *
);

std::vector<std::vector<std::string> >
compose_option_vector_of_vectors(
  const po::variables_map&,
  const char *,
  std::string,
  std::string,
  size_t
);

std::vector<std::vector<std::string> >
compose_option_vector_of_vectors(
  const po::variables_map&,
  const char *,
  size_t
);

std::vector<std::string>
compose_option_vector_of_strings(
  const po::variables_map&,
  const char *,
  std::string
);

} // user

#endif // PROGRAM_OPTIONS_H
