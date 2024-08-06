//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer and Michael Bojazi.
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
//! \brief Code for useful utilities for the simple snII project.
////////////////////////////////////////////////////////////////////////////////
  
#include "simple_snII_utilities.h"

//##############################################################################
// zone_compare_by_mass().
//##############################################################################

int
zone_compare_by_mass(
  const Libnucnet__Zone *p_zone1,
  const Libnucnet__Zone *p_zone2
)
{

  const char * s_1 =
    Libnucnet__Zone__getProperty( p_zone1, "mass below", NULL, NULL );

  const char * s_2 =
    Libnucnet__Zone__getProperty( p_zone2, "mass below", NULL, NULL );

  if( !s_1 || !s_2 )
  {
    std::cerr << "Zone lacks mass below property." << std::endl;
    exit( EXIT_FAILURE );
  }

  if( atof( s_1 ) < atof( s_2 ) )
     return -1;
  else
     return 1;

}

//##############################################################################
// get_zone_vectors().
//##############################################################################

boost::tuple<std::vector<double>, std::vector<double>, std::vector<double> >
get_zone_vectors(
  nnt::Zone& my_zone,
  Libnucnet * p_my_nucnet
)
{

  std::vector<double> zone_time, zone_t9, zone_log10_rho;

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    ( Libnucnet__Zone__compare_function ) zone_compare_by_mass
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );
  
  zone_time = 
    nnt::get_property_vector_from_zone_list<double>( zone_list, S_SHOCK_TIME );

  zone_time.insert( zone_time.begin(), 0. );
  
  zone_t9.push_back(
    my_zone.getProperty<double>( nnt::s_T9_0 )
  );
  
  zone_log10_rho.push_back(
    log10(
      my_zone.getProperty<double>( nnt::s_RHO_0 )
    )
  );
  
  double d_my_zone_time = 
      my_zone.getProperty<double>( S_SHOCK_TIME ); 
  
  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
    
    double d_t = 
      zone.getProperty<double>( S_SHOCK_TIME ); 

    if( d_t < d_my_zone_time )
    {

      zone_t9.push_back(
        my_zone.getProperty<double>( nnt::s_T9_0 )
      );

      zone_log10_rho.push_back(
        log10(
          my_zone.getProperty<double>( nnt::s_RHO_0 )
        )
      );

    }
    else
    {
      double d_rho = 
        pow(  
          zone.getProperty<double>( S_SHOCK_T9 ) / 
          my_zone.getProperty<double>( S_SHOCK_T9 ),
          3.
        ) * my_zone.getProperty<double>( S_SHOCK_RHO ); 

      zone_t9.push_back(
        zone.getProperty<double>( S_SHOCK_T9 )
      );
   
      zone_log10_rho.push_back(
        log10( d_rho )
      );
    }
  
  }       

  return boost::make_tuple( zone_time, zone_t9, zone_log10_rho );

}

//##############################################################################
// T9_post().
//##############################################################################

double T9_post( double d_t9, nnt::Zone& zone )
{
  
  zone.updateProperty(
    nnt::s_T9,
    d_t9
  );  
  
  double u = zone.getProperty<double>( "u" );

  double E =
    user::compute_thermo_quantity(
      zone,
      "internal energy density",
      zone.getProperty<std::string>( S_PARTICLE_TYPE )
    );

  return u - E;

} 

//##############################################################################
// rho_post().
//##############################################################################

double rho_post( double rho, nnt::Zone& zone )
{
  
  zone.updateProperty(
    nnt::s_RHO,
    rho
  );
 
  double rho0 = zone.getProperty<double>( nnt::s_RHO_0 );
    
  double e0 =  zone.getProperty<double>( S_E0 );

  double p0 = zone.getProperty<double>( S_P0 ); 

  double h0 = e0 + p0 / rho0;
  
  //==========================================================================
  // Find T9 root.          
  //==========================================================================

  double T9_peak =
    compute_1d_root(
      boost::bind( T9_post, _1, boost::ref( zone ) ),
      5. * zone.getProperty<double>( nnt::s_T9_0 ),
      2.
    );

  zone.updateProperty(
    nnt::s_T9,
    T9_peak
  );

  double e = 
    user::compute_thermo_quantity(
      zone,
      "internal energy density", 
      zone.getProperty<std::string>( S_PARTICLE_TYPE )
    ) / rho;
  
  zone.updateProperty(
    S_E,
    e
  );
 
  double p =
    user::compute_thermo_quantity(
      zone,
      nnt::s_PRESSURE,
      zone.getProperty<std::string>( S_PARTICLE_TYPE )
    );
  
  zone.updateProperty(
    S_P,
    p
  );

  double h = e + p / rho;
  
  return
    ( h - h0 )
    -
    ( 1. / 2. ) * ( 1. / rho0 + 1. / rho ) * ( p - p0 );
      
} 

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties(
  nnt::Zone& zone,
  Libnucnet * p_nucnet
)
{

  boost::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
  > t =
    get_zone_vectors( zone, p_nucnet );

  user::update_t9_rho_in_zone_by_interpolation(
    zone,
    "linear",
    t.get<0>(),
    t.get<1>(),
    t.get<2>()
  );

  double d_velocity =
    zone.getProperty<double>( S_SHOCK_SPEED ) -
    zone.getProperty<double>( S_V );

  if( zone.hasProperty( S_VELOCITY_FACTOR ) )
  {
    d_velocity *= zone.getProperty<double>( S_VELOCITY_FACTOR );
  }

  double d_delta_t =
    zone.getProperty<double>( nnt::s_TIME )
    -
    zone.getProperty<double>( S_SHOCK_TIME );

  if( d_delta_t <= 0. )
  {
    zone.updateProperty(
      nnt::s_RADIUS,
      zone.getProperty<std::string>( "cell outer radius" )
    );
  }
  else
  {
    zone.updateProperty(
      nnt::s_RADIUS,
      zone.getProperty<double>( "cell outer radius" ) + d_velocity * d_delta_t
    );
  }

} 

//##############################################################################
// get_evolution_network_and_zones().
//##############################################################################

std::pair<Libnucnet *, std::vector<nnt::Zone> >
get_evolution_network_and_zones( Libnucnet * p_nucnet, char ** argv )
{

  Libnucnet * p_new_nucnet, * p_tmp_nucnet;

  if( !argv )
  {
    std::cerr << "Invalid input." << std::endl;
    exit( EXIT_FAILURE );
  }

  p_new_nucnet = nnt::create_network_copy( p_nucnet );

  p_tmp_nucnet = nnt::create_network_copy( p_nucnet );

  Libnucnet__assignZoneDataFromXml( p_tmp_nucnet, argv[1], argv[3] );

  Libnucnet__setZoneCompareFunction(
    p_tmp_nucnet,
    ( Libnucnet__Zone__compare_function ) zone_compare_by_mass
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_tmp_nucnet );

  std::vector<nnt::Zone> zones( zone_list.size() );

  int i = 0;

  BOOST_FOREACH( nnt::Zone old_zone, zone_list )
  {

    zones[i].setNucnetZone(
      Libnucnet__Zone__new(
        Libnucnet__getNet( p_new_nucnet ),
        "work",
        boost::lexical_cast<std::string>( i ).c_str(),
        "0"
      )
    );

    Libnucnet__Zone__iterateOptionalProperties(
      old_zone.getNucnetZone(),
      NULL,
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function)
         nnt::copy_properties,
      zones[i].getNucnetZone()
    );

    gsl_vector * p_abundances =
      Libnucnet__Zone__getAbundances( old_zone.getNucnetZone() );

    Libnucnet__Zone__updateAbundances( zones[i].getNucnetZone(), p_abundances );

    gsl_vector_free( p_abundances );

    Libnucnet__addZone( p_new_nucnet, zones[i].getNucnetZone() );

    zones[i].updateProperty(
      "original label",
      Libnucnet__Zone__getLabel( old_zone.getNucnetZone(), 1 )
    );

    zones[i].updateProperty(
      S_XML_OUTPUT_FILE,
      argv[4]
    );

    i++;

  }

  Libnucnet__free( p_tmp_nucnet );

  return
    std::make_pair( p_new_nucnet, zones );

}

//##############################################################################
// copy_properties_to_zone().
//##############################################################################

int
copy_properties_to_zone(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  Libnucnet__Zone * p_zone
)
{

  Libnucnet__Zone__updateProperty(
    p_zone,
    s_name,
    s_tag1,
    s_tag2,
    s_value
  );

  return 1;

}

//##############################################################################
// set_low_T_network().
//##############################################################################

void
set_low_T_network( nnt::Zone& zone )
{

  if(
    zone.hasProperty( S_LOW_T_NETWORK_T9 ) &&
    zone.getProperty<double>( nnt::s_T9 ) <
      zone.getProperty<double>( S_LOW_T_NETWORK_T9 )
  )
  {
    zone.updateProperty(
      nnt::s_BASE_EVOLUTION_REAC_XPATH,
      zone.getProperty<std::string>( S_LOW_T_NETWORK_REAC_XPATH )
    );
  }

}
