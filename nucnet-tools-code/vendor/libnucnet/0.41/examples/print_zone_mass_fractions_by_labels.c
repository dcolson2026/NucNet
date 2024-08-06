/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
//
//     This is free software; you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     Please see the src/README.txt file in this distribution for more
//     information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use libnucnet routines to create
//       a new Libnucnet structure from an input xml file, print out 
//       abundances for the given zone, and clear the
//       structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet.h>

void print_out_zone( Libnucnet__Zone * );
int print_species_abundance( Libnucnet__Species *, void * );

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;

  if ( argc != 6 ) {
      fprintf(
        stderr,
        "\nUsage: %s nuc_file mass_file label1 label2 label3\n\n", argv[0]
      );
      fprintf(
        stderr, "  nuc_file = input nuclear network data xml filename\n\n"
      );
      fprintf(
        stderr, "  mass_file = input mass fraction xml filename\n\n"
      );
      fprintf(
        stderr, "  label1 = first zone label\n\n"
      );
      fprintf(
        stderr, "  label2 = second zone label\n\n"
      );
      fprintf(
        stderr, "  label3 = third zone label\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_my_nucnet )
    ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  /*============================================================================
  // Retrieve zone and print abundances.
  //==========================================================================*/

  printf( "\n" );

  p_zone = Libnucnet__getZoneByLabels( 
    p_my_nucnet, argv[3], argv[4], argv[5] 
  );

  if( !p_zone ) {

    fprintf( stderr, "Zone not found!\n" );
    return EXIT_FAILURE;

  }

  printf( "For zone with labels %s %s %s:\n\n",
          Libnucnet__Zone__getLabel( p_zone, 1 ),
          Libnucnet__Zone__getLabel( p_zone, 2 ),
          Libnucnet__Zone__getLabel( p_zone, 3 )
  );

  print_out_zone( p_zone );

  printf( "\n" );

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// print_out_zone()
//############################################################################*/

void print_out_zone( Libnucnet__Zone *p_zone ) {

  struct
  {
    Libnucnet__Zone * pZone;
    gsl_vector * pX;
  } user_data;

  user_data.pZone = p_zone;
  user_data.pX = Libnucnet__Zone__getMassFractions( p_zone );

  printf( "Species\t Abundance\tMass Fraction\n" );
  printf( "=======\t============\t=============\n" );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Species__iterateFunction) print_species_abundance,
    &user_data
  );

  printf(
    "\nSum of mass fractions = %e\n",
    Libnucnet__Zone__computeAMoment( p_zone, 1 )
  );

  gsl_vector_free( user_data.pX );

}

int
print_species_abundance(
  Libnucnet__Species *p_species,
  void * p_data
)
{

  typedef struct {
    Libnucnet__Zone * pZone;
    gsl_vector * pX;
  } user_data;

  user_data my_data = *( user_data * ) p_data;

  size_t i_index = Libnucnet__Species__getIndex( p_species );

  double d_x = gsl_vector_get( my_data.pX, i_index );

  if( d_x > 0 ) {
     printf(
       "%s\t%e\t%e\n",
       Libnucnet__Species__getName( p_species ),
       d_x / Libnucnet__Species__getA( p_species ),
       d_x
     ); 
  }

  return 1;

}

