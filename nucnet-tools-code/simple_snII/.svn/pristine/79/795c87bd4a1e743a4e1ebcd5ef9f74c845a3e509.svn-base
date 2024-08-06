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
//! \brief The simple_snII explosion code.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>

#include "nnt/string_defs.h"
#include "nnt/auxiliary.h"

#include "simple_snII_utilities.h"

#define S_MY_PARTICLE_TYPE  "photon"  // "photon" or "total"
#define S_SOUND_SPEED       "sound speed"

#define B_COMPUTE_SOUND_SPEED  false  // Change to true to compute
                                      // pre-supernova sound speed and Mach
                                      // number.

#define B_INCLUDE_PRE_ENERGY   false  // Change to true to include
                                      // pre-supernova energy density.

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] )
{  

  Libnucnet * p_my_nucnet;
  
  double d_time = 0., d_r = 0., EE = 0;
  double f1, f2, f3;
  
  //============================================================================
  // Check input.
  //============================================================================

  if( argc != 4 )
  {

    std::cerr <<
      std::endl <<
      "Usage: " << argv[0] << " input_xml energy output_xml" <<
      std::endl << std::endl <<
      "  input_xml = input presupernova xml file" <<
      std::endl << std::endl <<
      "  energy = energy of shock (ergs)" <<
      std::endl << std::endl <<
      "  output_xml = output xml file" <<
      std::endl << std::endl;

    return EXIT_FAILURE;

  }
 
  //============================================================================
  // Get input.
  //============================================================================

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      NULL
    );

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    ( Libnucnet__Zone__compare_function ) zone_compare_by_mass
  );
     
  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
  
    double E_0 = atof( argv[2] );

    zone.updateProperty( S_SHOCK_ENERGY, argv[2] );

    zone.updateProperty( S_PARTICLE_TYPE, S_MY_PARTICLE_TYPE );

    //==========================================================================
    // Compute and store pre-shock quantities.
    //==========================================================================

    zone.updateProperty(
      nnt::s_T9,
      zone.getProperty<double>( "cell temperature" )
      / GSL_CONST_NUM_GIGA
    );

    zone.updateProperty(
      nnt::s_T9_0,
      zone.getProperty<std::string>( nnt::s_T9 )
    );

    zone.updateProperty(
      nnt::s_RHO,
      zone.getProperty<std::string>( nnt::s_RHO1 )
    );

    zone.updateProperty(
      nnt::s_RHO_0,
      zone.getProperty<std::string>( nnt::s_RHO )
    );
    
    double p0 = user::compute_thermo_quantity( zone, "pressure", "total" );

    zone.updateProperty(
      S_P0,
      p0
    );
    
    double e0 =
      user::compute_thermo_quantity(
        zone, "internal energy density", "total"
      );

    if( B_INCLUDE_PRE_ENERGY )
    {
      EE +=
        e0 * (4./3. ) * ( M_PI ) *
        (
          gsl_pow_3( zone.getProperty<double>( "cell outer radius" ) ) -
          gsl_pow_3( d_r )
        );
    }

    e0 /= zone.getProperty<double>( nnt::s_RHO );

    zone.updateProperty(
      "EE",
      EE
    );
  
    zone.updateProperty(
      S_E0,
      e0
    );
  
    double u = ( E_0 + EE ) / ( ( 4./3. ) * ( M_PI ) * 
      pow( 
        zone.getProperty<double>( "cell outer radius" ),
        3. 
      ) 
    ); 
  
    zone.updateProperty(
      S_U,
      u
    );

    if( B_COMPUTE_SOUND_SPEED )
    {
      zone.updateProperty(
        S_SOUND_SPEED,
        user::compute_sound_speed( zone )
      ); 
    }

    //==========================================================================
    // Find rho root.
    //==========================================================================
    
    double rho = 
      compute_1d_root(
        boost::bind( rho_post, _1, boost::ref( zone ) ),
        zone.getProperty<double>( nnt::s_RHO_0 ),
        2.
      );
        
    zone.updateProperty(
      nnt::s_RHO,
      rho
    );
    
    //==========================================================================
    // Compute and store post shock conditions.
    //==========================================================================

      std::cout <<
        zone.getProperty<std::string>( nnt::s_RHO_0 ) <<
        "  " <<
        zone.getProperty<std::string>( nnt::s_RHO ) <<
        "  " <<
        zone.getProperty<double>( nnt::s_RHO ) /
        zone.getProperty<double>( nnt::s_RHO_0 ) << 
        std::endl;
    
    double e =
      user::compute_thermo_quantity(
        zone,
        "internal energy density",
        S_MY_PARTICLE_TYPE
      ) /
      zone.getProperty<double>( nnt::s_RHO ); 

    zone.updateProperty(
      S_E,
      e
    );
    
    double p =
      user::compute_thermo_quantity(
        zone,
        "pressure",
        S_MY_PARTICLE_TYPE
      );
    
    zone.updateProperty(
      S_P,
      p
    );
    
    double rho0 = zone.getProperty<double>( nnt::s_RHO_0 );
   
    double v_s = 
      pow( ( p - p0 ) / 
        ( rho0 - pow( rho0, 2. ) / 
          zone.getProperty<double>( nnt::s_RHO )
        ),
        1./2.
      ); 

    zone.updateProperty(
      S_SHOCK_T9,
      zone.getProperty<std::string>( nnt::s_T9 )
    );

    zone.updateProperty(
      S_SHOCK_RHO,
      zone.getProperty<std::string>( nnt::s_RHO )
    );

    zone.updateProperty(
      S_SHOCK_SPEED,
      v_s
    );
    
    if( B_COMPUTE_SOUND_SPEED )
    {
      zone.updateProperty(
        nnt::s_MACH,
        v_s / zone.getProperty<double>( S_SOUND_SPEED )
      );
    }
    
    double v =
      v_s *
      ( rho0 / zone.getProperty<double>( nnt::s_RHO ) );

    zone.updateProperty(
      S_V,
      v
    );   
   
    d_time +=
      (
        zone.getProperty<double>( "cell outer radius" )
        -
        d_r
      )
      /
      v_s;

    zone.updateProperty(
      S_SHOCK_TIME,
      d_time
    );

    d_r = zone.getProperty<double>( "cell outer radius" );

    //==========================================================================
    // Compute and store Rankine-Hugioniot equation solutions.
    //==========================================================================
    
    f1 = 
      ( ( rho0 * v_s ) - 
        ( zone.getProperty<double>( nnt::s_RHO ) * v ) 
      ) / ( rho0 * v_s ); 
  
    zone.updateProperty(
      S_F1,
      f1
    );   
   
    f2 = ( ( p0 + ( rho0 * pow( v_s, 2. ) ) ) - ( p + 
           ( zone.getProperty<double>( nnt::s_RHO ) * 
               pow( v, 2. ) ) ) ) /
                 ( p + ( rho0 * pow( v_s, 2. ) ) );    
                                
    zone.updateProperty(
      S_F2,
      f2
    );   

    f3 =
      (
        e0 + ( p0 / rho0 ) + ( 1./2. ) * pow( v_s, 2. ) - 
        (
          e
          +
          p / zone.getProperty<double>( nnt::s_RHO )
          +
          ( 1./2. ) * pow( v, 2. )
        )
      ) / 
      ( e0 + ( p0 / rho0 ) + ( 1./2. ) * pow( v_s, 2. ) );
 
    zone.updateProperty(
      S_F3,
      f3
    );   
 
  } //BOOST_FOREACH

  //============================================================================
  // Write to xml and clean up.
  //============================================================================
    
  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_nucnet,
    "%.15e"
  );

  Libnucnet__writeToXmlFile(
    p_my_nucnet,
    argv[3]
  );

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
