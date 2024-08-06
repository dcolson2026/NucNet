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
//       Example to demonstrate how to use libnucnet routines to parse
//       in a Libnucnet__Nuc structure of nuclear species from an input
//       xml file, extract a subset of those data by an xpath expression,
//       save the subset to a new file,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Nuc.h>

int
reverse_sort_function(
  const Libnucnet__Species *,
  const Libnucnet__Species *
);

/*##############################################################################
// main()
//############################################################################*/

int
main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei, *p_subset;

  if ( argc != 5 ) {
      fprintf(
        stderr, "\nUsage: %s filename xpath indent new_filename\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  xpath = xpath expression for subset.\n\n"
      );
      fprintf(
        stderr, "  indent = indent output (0 = no, 1 = yes).\n\n"
      );
      fprintf(
        stderr, "  new_filename = name of file to dump subset.\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Create nuclear species collection.
  //==========================================================================*/

  p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );
    
  /*============================================================================
  // Set the species comparison function.
  //==========================================================================*/

  Libnucnet__Nuc__setSpeciesCompareFunction(
    p_my_nuclei,
    (Libnucnet__Species__compare_function) reverse_sort_function
  );

  /*============================================================================
  // Create the subset collection.
  //==========================================================================*/

  p_subset = Libnucnet__Nuc__extractSubset( p_my_nuclei, argv[2] );

  /*============================================================================
  // Set the species comparison function in the subset to that in the
  // original collection and sort.
  //==========================================================================*/

  Libnucnet__Nuc__setSpeciesCompareFunction(
    p_subset,
    Libnucnet__Nuc__getSpeciesCompareFunction( p_my_nuclei )
  );

  Libnucnet__Nuc__sortSpecies( p_subset );

  /*============================================================================
  // Free the original nuclear species collection.
  //==========================================================================*/

  Libnucnet__Nuc__free( p_my_nuclei );
  
  /*============================================================================
  // Write the subset to a new xml file.  Use text writer.
  //==========================================================================*/

  Libnucnet__Nuc__writeToXmlFileWithTextWriter(
    p_subset, argv[4], atoi( argv[3] )
  );

  /*============================================================================
  // Free the subset.
  //==========================================================================*/

  Libnucnet__Nuc__free( p_subset );
  
  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// reverse_sort_function()
//############################################################################*/

int
reverse_sort_function(
  const Libnucnet__Species *p_species1,
  const Libnucnet__Species *p_species2
)
{

  int i;

  if( 
      Libnucnet__Species__getZ( p_species1 ) >
      Libnucnet__Species__getZ( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getZ( p_species1 ) <
      Libnucnet__Species__getZ( p_species2 )
  )
    return 1;

  if( 
      Libnucnet__Species__getA( p_species1 ) >
      Libnucnet__Species__getA( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) <
      Libnucnet__Species__getA( p_species2 )
  )
    return 1;

  i =
    strcmp(
      Libnucnet__Species__getName( p_species1 ),
      Libnucnet__Species__getName( p_species2 )
    );

  if( i == 0 ) {
    return 0;
  } else {
    return -GSL_SIGN( i );
  }

}
   
