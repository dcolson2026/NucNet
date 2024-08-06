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
//! \brief Code for running multiple zone simple snII with OpenMP.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <omp.h>

#include "user/remove_duplicate.h"
#include "user/user_rate_functions.h"
#include "user/network_limiter.h"
#include "user/rate_modifiers.h"
#include "user/neutrino_rate_functions.h"
#include "user/evolve.h"

#include "my_user/my_rates.h"
#include "simple_snII_utilities.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_REG_T        0.15    /* Time step regulator */
#define D_REG_Y        0.15    /* Abundance change regulator */
#define D_Y_MIN_DT     1.e-10  /* Abundance minimum for timestep */


//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

#define VALIDATE       "yes"

//##############################################################################
// Strings.
//##############################################################################

#define S_SOLVER_TYPE   nnt::s_ARROW  /* nnt::s_ARROW or nnt::s_GSL */
#define S_DETAILED_WEAK_RATES  "detailed weak rates"
#define S_NEUTRINO_RATES_FILE  "neutrino rates file"
#define S_USE_APPROXIMATE_WEAK_RATES  "use approximate weak rates"
#define S_PARAMS "params"
#define S_UPDATED_RATES_XML "updated rates xml file"

//##############################################################################
// Prototypes.
//##############################################################################

void
set_param_zone_views(
  nnt::Zone&
);

void
set_zone(
  nnt::Zone&,
  nnt::Zone&
);

void
evolve_zone(
  Libnucnet *,
  nnt::Zone&,
  double
);

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i;
  std::vector<nnt::Zone> zone_vector;
  nnt::Zone param_zone;
  Libnucnet *p_my_nucnet, *p_shock_nucnet, * p_tmp;
  Libnucnet__Reac *p_duplicates;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 5 )
  {
    fprintf(
      stderr,
      "\nUsage: %s input_xml param_xml out_xml zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  input_xml = input explosion xml filename\n\n"
    );
    fprintf(
      stderr, "  param_xml = zone parameter xml filename\n\n"
    );
    fprintf(
      stderr, "  out_xml = output xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_xpath = xpath to select zones\n\n"
    );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Validate input file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      return EXIT_FAILURE;
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  p_shock_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_shock_nucnet ) ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( 
    p_shock_nucnet,
    argv[1],
    NULL
  );

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      argv[4]
    );  

  //============================================================================
  // Register user-supplied rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

/*
  my_user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );
*/

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  p_duplicates =
    Libnucnet__Reac__getDuplicateReactions(
      Libnucnet__Net__getReac(
        Libnucnet__getNet( p_my_nucnet )
      )
    );

  Libnucnet__Reac__iterateReactions(
    p_duplicates,
    (Libnucnet__Reaction__iterateFunction) user::remove_duplicate,
    Libnucnet__getNet( p_my_nucnet )
  );

  Libnucnet__Reac__free( p_duplicates );

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if(
    strcmp( S_SOLVER_TYPE, nnt::s_ARROW ) == 0
  )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

  }

  //============================================================================
  // Create parameter zone.
  //============================================================================

  p_tmp = Libnucnet__new();

  Libnucnet__assignZoneDataFromXml( p_tmp, argv[2], NULL );

  param_zone.setNucnetZone(
    Libnucnet__Zone__new(
      Libnucnet__getNet( p_my_nucnet ),
      S_PARAMS,
      "0",
      "0"
    )
  );

  Libnucnet__Zone__iterateOptionalProperties(
    Libnucnet__getZoneByLabels( p_tmp, "0", "0", "0" ),
    NULL,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       copy_properties_to_zone,
    param_zone.getNucnetZone()
  );

  Libnucnet__free( p_tmp );

  //============================================================================
  // Use approximate weak rates or not.
  //============================================================================

  if( param_zone.hasProperty( S_USE_APPROXIMATE_WEAK_RATES ) )
  {
    if(
      param_zone.getProperty<std::string>( S_USE_APPROXIMATE_WEAK_RATES ) ==
      "yes"
    )
      user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );
  }

  //============================================================================
  // Update with detailed weak rates.
  //============================================================================

  if( param_zone.hasProperty( S_DETAILED_WEAK_RATES ) )
  {

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( param_zone.getNucnetZone() )
      ),
      param_zone.getProperty<std::string>( S_DETAILED_WEAK_RATES ).c_str(),
      NULL
    );

    user::set_two_d_weak_rates_hashes(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
    );

  }

  //============================================================================
  // Update with neutrino rates.
  //============================================================================

  if( param_zone.hasProperty( S_NEUTRINO_RATES_FILE ) )
  {

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( param_zone.getNucnetZone() )
      ),
      param_zone.getProperty<std::string>( S_NEUTRINO_RATES_FILE ).c_str(),
      NULL
    );

    user::set_nu_nucl_hash(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( param_zone.getNucnetZone() )
      )
    );

    user::register_neutrino_rate_functions(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( param_zone.getNucnetZone() )
      )
    );

  }

  //============================================================================
  // Update with added rates.
  //============================================================================

  if( param_zone.hasProperty( S_UPDATED_RATES_XML ) )
  {
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      param_zone.getProperty<std::string>( S_UPDATED_RATES_XML ).c_str(),
      NULL
    );
  }

  //============================================================================
  // Create zone vector.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) zone_compare_by_mass
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    zone_vector.push_back( zone );

  }
  
  //============================================================================
  // Set the properties and view in the parameters zone.
  //============================================================================

  set_param_zone_views( param_zone );

  //============================================================================
  // Evolve zones.
  //============================================================================

  #pragma omp parallel for schedule( dynamic, 1 )
    for(
      i = 0;
      i < (int) Libnucnet__getNumberOfZones( p_my_nucnet );
      i++
    )
    {

      set_zone( param_zone, zone_vector[i] );

      evolve_zone(
        p_shock_nucnet,
        zone_vector[i],
        zone_vector[i].getProperty<double>( nnt::s_TEND )
      );

      Libnucnet__Zone__clearNetViews( zone_vector[i].getNucnetZone() );

      Libnucnet__Zone__clearRates( zone_vector[i].getNucnetZone() );

    }
 
  //============================================================================
  // Write out.
  //============================================================================
 
  Libnucnet__writeToXmlFile( p_my_nucnet, argv[3] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Zone__free( param_zone.getNucnetZone() );

  Libnucnet__free( p_my_nucnet );
  Libnucnet__free( p_shock_nucnet );

  return EXIT_SUCCESS;

}

//##############################################################################
// set_zone().
//##############################################################################

void
set_zone(
  nnt::Zone& param_zone,
  nnt::Zone& zone
)
{

  nnt::normalize_zone_abundances( zone );

  Libnucnet__Zone__iterateOptionalProperties(
    param_zone.getNucnetZone(),
    NULL,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       copy_properties_to_zone,
    zone.getNucnetZone()
  );

  Libnucnet__Zone__copy_net_views(
    zone.getNucnetZone(),
    param_zone.getNucnetZone()
  );

  zone.updateProperty(
    nnt::s_RADIUS,
    zone.getProperty<std::string>( "cell outer radius" )
  );

  user::set_rate_data_update_function( zone );

  user::limit_evolution_network( zone );

}

//##############################################################################
// evolve_zone().
//##############################################################################

void
evolve_zone(
  Libnucnet * p_nucnet,
  nnt::Zone& zone,
  double d_tend
)
{  

  int i_steps = 0;
  double d_dt, d_t;

  d_t = zone.getProperty<double>( nnt::s_TIME );

  d_dt = zone.getProperty<double>( nnt::s_DTIME );

     std::cout <<
       "Starting zone " <<
       Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) <<
       "  1 - xsum = " <<
       1. - Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) <<
       std::endl;

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  while ( d_t < zone.getProperty<double>( nnt::s_TEND ) )
  {

    d_t += d_dt;

  //============================================================================
  // Set dt and t.
  //============================================================================

    zone.updateProperty(
      nnt::s_DTIME,
      d_dt
    );

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

  //============================================================================
  // Update temperature and density.  Update dt and time again in case
  // changed in update_zone_properties.
  //============================================================================

    update_zone_properties( zone, p_nucnet );

    d_t = zone.getProperty<double>( nnt::s_TIME );

    d_dt = zone.getProperty<double>( nnt::s_DTIME );

  //============================================================================
  // Evolve abundances.
  //============================================================================

    user::evolve( zone );

  //============================================================================
  // Update exposures.
  //============================================================================

    user::update_exposures( zone );

  //============================================================================
  // Update timestep.
  //============================================================================

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if( zone.getProperty<double>( nnt::s_T9 ) > 10. )
      nnt::normalize_zone_abundances( zone );

    if(
      d_t + d_dt > zone.getProperty<double>( nnt::s_TEND )
    )
    {
      d_dt = zone.getProperty<double>( nnt::s_TEND ) - d_t;
    }

    set_low_T_network( zone );

    user::limit_evolution_network( zone );

    i_steps++;

  }  

  zone.updateProperty(
    "calculation steps",
    i_steps
  );

    std::cout << "   Finishing zone " <<
       Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) <<
       "  1 - xsum = " <<
       1. - Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) <<
       std::endl;

}

//##############################################################################
// set_param_zone_views().
//##############################################################################

void
set_param_zone_views(
  nnt::Zone& param_zone
)
{

  user::view_multi views;

  //============================================================================
  // Get rate modification views.
  //============================================================================

  Libnucnet__Zone__iterateOptionalProperties(
    param_zone.getNucnetZone(),
    NULL,
    nnt::s_RATE_MODIFICATION_VIEW,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       user::modify_rates_views,
    &views
  );

  //============================================================================
  // Assign views.
  //============================================================================

  param_zone.getNetView( "", "" );

  if(
    param_zone.hasProperty( S_LOW_T_NETWORK_T9 ) &&
    param_zone.hasProperty( S_LOW_T_NETWORK_REAC_XPATH )
  )
  {
    param_zone.getNetView(
      "",
      param_zone.getProperty<std::string>(
        S_LOW_T_NETWORK_REAC_XPATH
      ).c_str()
    );
  }

  BOOST_FOREACH( user::view view, views )
  {

    param_zone.getNetView(
      view.nuc_xpath.c_str(),
      view.reac_xpath.c_str()
    );

  }

}
