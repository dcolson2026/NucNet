////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Michael J. Bojazi and Bradley S. Meyer.
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
//! \brief Example code for running a single zone network calculation.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>

#include "nnt/two_d_weak_rates.h"
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

#define D_REG_T        0.15         // Time step change regulator for dt update
#define D_REG_Y        0.15         // Abundance change regulator for dt update 
#define D_Y_MIN_DT     1.e-10       // Smallest y for dt update

#define S_SOLVER       nnt::s_ARROW // Solver type: ARROW or GSL
#define S_DETAILED_WEAK_RATES  "detailed weak rates"
#define S_NEUTRINO_RATES_FILE  "neutrino rates file"

#define S_USE_APPROXIMATE_WEAK_RATES  "use approximate weak rates"
#define S_DETAILED_WEAK_RATES  "detailed weak rates"
#define S_FLOW_CURRENT_XML_FILE  "flow current xml file"
#define S_INTEGRATED_CURRENTS  "integrated currents"
#define S_UPDATED_RATES_XML  "updated rates xml file"

//##############################################################################
// get_nucnet().
//##############################################################################

Libnucnet *
get_nucnet( int argc, char **argv )
{

  Libnucnet * p_nucnet, * p_tmp;
  Libnucnet__Zone * p_new_zone;
  gsl_vector * p_abunds;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "\"[@label1 = '614]\" my_output.xml \"[z <= 20]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if ( argc < 5 || argc > 7 || strcmp( argv[1], "--usage" ) == 0 )
  {

    fprintf(
      stderr,
      "\nUsage: %s expl_xml zone_params zone_xpath out_file xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  expl_xml = input explosion xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_params = xml file with extra data\n\n"
    );
    fprintf(
      stderr, "  zone_xpath = xpath to select zone\n\n"
    );
    fprintf(
      stderr, "  out_file = output xml file\n\n"
    );
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified)\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional)\n\n"
    );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  if( argc == 5 ) {
    p_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );
  } else if( argc == 6 ) {
    p_nucnet = Libnucnet__new_from_xml( argv[1], argv[5], NULL, NULL );
  } else {
    p_nucnet = Libnucnet__new_from_xml( argv[1], argv[5], argv[6], NULL );
  }

  //============================================================================
  // Get parameter zone.
  //============================================================================
  
  p_tmp = Libnucnet__new();

  if( argc == 5 )
  {
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_tmp ) ),
      argv[1],
      NULL
    );
  }
  else
  {
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_tmp ) ),
      argv[1],
      argv[5]
    );
  }

  Libnucnet__assignZoneDataFromXml( p_tmp, argv[2], NULL ); 

  nnt::zone_list_t param_zone_list = nnt::make_zone_list( p_tmp );

  if( param_zone_list.size() != 1 )
  {
    std::cerr << "Incorrect number of parameter zones." << std::endl;
    exit( EXIT_FAILURE );
  }

  p_new_zone =
    Libnucnet__Zone__new(
      Libnucnet__getNet( p_nucnet ),
      "params",
      "0",
      "0"
    );

  Libnucnet__Zone__iterateOptionalProperties(
    param_zone_list.begin()->getNucnetZone(),
    NULL,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
      copy_properties_to_zone,
    p_new_zone
  );

  p_abunds =
    Libnucnet__Zone__getAbundances(
      param_zone_list.begin()->getNucnetZone()
    );

  Libnucnet__Zone__updateAbundances( p_new_zone, p_abunds );

  gsl_vector_free( p_abunds );
    
  Libnucnet__addZone( p_nucnet, p_new_zone );

  Libnucnet__free( p_tmp );

  //============================================================================
  // Done.
  //============================================================================

  return p_nucnet;

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step, k = 0;
  double d_t, d_dt;
  Libnucnet *p_my_nucnet, *p_nucnet_shock, *p_flow_current_nucnet = NULL;
  Libnucnet *p_my_output;
  Libnucnet__Zone * p_zone, * p_zone_params;
  std::vector<nnt::Zone> zones;
  nnt::Zone zone, flow_current_zone;
  Libnucnet__Reac *p_duplicates;
  gsl_vector * p_abunds;

  //============================================================================
  // Get the nucnet and zone.
  //============================================================================

  p_nucnet_shock = get_nucnet( argc, argv );

  boost::tie( p_my_nucnet, zones ) =
     get_evolution_network_and_zones( p_nucnet_shock, argv );

  if( zones.size() != 1 )
  {
    std::cerr << "This routine requires a single zone." << std::endl;
    return EXIT_FAILURE;
  }

  zone = *(zones.begin());

  p_zone_params =
    Libnucnet__getZoneByLabels( p_nucnet_shock, "params", "0", "0" );

  Libnucnet__Zone__iterateOptionalProperties(  
    p_zone_params,
    NULL,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       copy_properties_to_zone,
    zone.getNucnetZone()
  );

  p_abunds = Libnucnet__Zone__getAbundances( p_zone_params );

  if( !gsl_vector_isnull( p_abunds ) )
  {
    Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abunds );
  }

  gsl_vector_free( p_abunds );

  Libnucnet__removeZone( p_nucnet_shock, p_zone_params );

  //============================================================================
  // Update reactions, if desired.
  //============================================================================
  
  if( zone.hasProperty( S_UPDATED_RATES_XML ) )
  {

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      zone.getProperty<std::string>( S_UPDATED_RATES_XML ).c_str(),
      NULL
    );

  }
  
  //============================================================================
  // Register rate functions.
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
  // Initialize time.
  //============================================================================

  if( zone.hasProperty( nnt::s_DTIME ) )
    d_dt = zone.getProperty<double>( nnt::s_DTIME );
  else
  {
    d_dt = 1.e-05;
    zone.updateProperty(
      nnt::s_DTIME,
      d_dt
    );
  }

  if( zone.hasProperty( nnt::s_TIME ) )
    d_t = zone.getProperty<double>( nnt::s_TIME );
  else
  {
    d_t = 0;

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

  }

  //============================================================================
  // Use approximate weak rates or not.
  //============================================================================

  if( zone.hasProperty( S_USE_APPROXIMATE_WEAK_RATES ) )
  {
    if( zone.getProperty<std::string>( S_USE_APPROXIMATE_WEAK_RATES ) == "yes" )
      user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );
  }

  //============================================================================
  // Update with detailed weak rates.
  //============================================================================

  if( zone.hasProperty( S_DETAILED_WEAK_RATES ) )
  {

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      zone.getProperty<std::string>( S_DETAILED_WEAK_RATES ).c_str(),
      NULL
    );

    user::set_two_d_weak_rates_hashes(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
    );

  }

  //============================================================================
  // Update with neutrino rates.
  //============================================================================

  if( zone.hasProperty( S_NEUTRINO_RATES_FILE ) )
  {

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      zone.getProperty<std::string>( S_NEUTRINO_RATES_FILE ).c_str(),
      NULL
    );

    user::set_nu_nucl_hash(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      )
    );

    user::register_neutrino_rate_functions(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      )
    );

  }

  //============================================================================
  // If not set, set the neutrino chemical potential to -inf.
  //============================================================================
  
  if( !zone.hasProperty( nnt::s_MU_NUE_KT ) )
    zone.updateProperty(
      nnt::s_MU_NUE_KT,
      "-inf"
    ); 

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

  if( S_SOLVER == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    zone.updateProperty( nnt::s_SOLVER, nnt::s_ARROW );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "3" );

  }

  //============================================================================
  // Normalize abundances.
  //============================================================================

  nnt::normalize_zone_abundances( zone );

  //============================================================================
  // Set iteration order.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) zone_compare_by_mass
  );

  //============================================================================
  // Print out the rate modification views.
  //============================================================================

  user::print_modified_reactions( zone );

  //============================================================================
  // Limit evolution network.
  //============================================================================

  user::limit_evolution_network( zone );

  //============================================================================
  // Set the rate date update function.
  //============================================================================

  user::set_rate_data_update_function( zone );

  //============================================================================
  // Create output.
  //============================================================================

  p_my_output = nnt::create_network_copy( p_my_nucnet );

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) zone_compare_by_mass
  );

  //============================================================================
  // Create current nucnet.
  //============================================================================

  if( zone.hasProperty( S_FLOW_CURRENT_XML_FILE ) )
  {

    p_flow_current_nucnet = nnt::create_network_copy( p_my_nucnet );

    flow_current_zone.setNucnetZone(
      Libnucnet__Zone__new(
        Libnucnet__getNet( p_flow_current_nucnet ),
        S_INTEGRATED_CURRENTS,
        NULL,
        NULL
      )
    );

    Libnucnet__addZone(
      p_flow_current_nucnet,
      flow_current_zone.getNucnetZone()
    );

    nnt::species_list_t species_list =
      nnt::make_species_list(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
      );

    BOOST_FOREACH( nnt::Species species, species_list )
    {
      flow_current_zone.updateProperty(
        "initial abundance",
        Libnucnet__Species__getName( species.getNucnetSpecies() ),
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies()
        )
      );
    }

  }

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;

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
  // Update temperature, density, and radius.  Update dt and time again in case
  // changed in update_zone_properties.
  //============================================================================

    update_zone_properties( zone, p_nucnet_shock );

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
  // Print out abundances.
  //============================================================================

    if(
       (
         i_step % zone.getProperty<int>( nnt::s_STEPS ) == 0 ||
         d_t >= zone.getProperty<double>( nnt::s_TEND )
       )
    )
    {
      p_zone = Libnucnet__Zone__copy( zone.getNucnetZone() );
      Libnucnet__addZone( p_my_output, p_zone );
      Libnucnet__relabelZone(
        p_my_output,
        p_zone,
        ( boost::lexical_cast<std::string>( ++k ) ).c_str(),
        NULL,
        NULL
      );
      nnt::print_zone_abundances( zone );
    }

    if( zone.hasProperty( S_FLOW_CURRENT_XML_FILE ) )
    {
      user::update_flow_currents( zone, flow_current_zone );
    }

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

    i_step++;

  }  

  //============================================================================
  // Write out currents.
  //============================================================================

  if( zone.hasProperty( S_FLOW_CURRENT_XML_FILE ) )
  {

    nnt::species_list_t species_list =
      nnt::make_species_list(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
      );

    BOOST_FOREACH( nnt::Species species, species_list )
    {
      flow_current_zone.updateProperty(
        "final abundance",
        Libnucnet__Species__getName( species.getNucnetSpecies() ),
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies()
        )
      );
    }

    Libnucnet__writeToXmlFile(
      p_flow_current_nucnet,
      zone.getProperty<std::string>( S_FLOW_CURRENT_XML_FILE ).c_str()
    );

    Libnucnet__free( p_flow_current_nucnet );

  }

  //============================================================================
  // Remove work zone.
  //============================================================================

  Libnucnet__removeZone(
    p_my_nucnet,
    Libnucnet__getZoneByLabels( p_my_nucnet, "work", "0", "0" )
  );

  //============================================================================
  // Write out full output.
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_output,
    "%.15e"
  );

  Libnucnet__writeToXmlFile(
    p_my_output,
    Libnucnet__Zone__getProperty(
      Libnucnet__getZoneByLabels( p_my_output, "1", "0", "0" ),
      S_XML_OUTPUT_FILE,
      NULL,
      NULL
    )
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_nucnet_shock );
  Libnucnet__free( p_my_output );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
