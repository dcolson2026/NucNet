#///////////////////////////////////////////////////////////////////////////////
#  This file was originally written by Bradley S. Meyer.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file create_reaction_xml.sh
#//! \brief A short shell script to generate reaction xml files from individual
#//!        JINA Reaclib2 format data files.
#//!
#///////////////////////////////////////////////////////////////////////////////

#/bin/bash

#///////////////////////////////////////////////////////////////////////////////
# Check input.
#///////////////////////////////////////////////////////////////////////////////

if [ ! $# -eq 3 -a "$1" == "--example" ]
  then
    echo -e "\n$0 ../nucnet-tools-code/data_pub/jina_to_webnucleo/individual_reaction_data xml reactions.xml\n"
    exit
fi

if [ ! $# -eq 3 ]
  then
    echo -e "\nUsage: $0 input_dir output_dir list_xml\n"
    echo -e "  input_dir = directory containing the reaclib data files\n"
    echo -e "  output_dir = directory to contain output xml file\n"
    echo -e "  list_xml = output xml file that lists the individual files\n"
    echo -e "For an example, type: $0 --example\n"
    exit
fi

#///////////////////////////////////////////////////////////////////////////////
# Output xml files.  First destroy any previous version of the output directory
# to prevent extraneous files from being included in reaction list xml.
#///////////////////////////////////////////////////////////////////////////////

rm -fr $2

mkdir $2 

for file in $1/*
do
  ./reaction_converter $file $2/${file##*/}.xml
done

#///////////////////////////////////////////////////////////////////////////////
# Create xml file with list of individual xml files and move it to the output
# directory.  Note that the xpointer ensures schema validation (use --xinclude
# with xmllint).
#///////////////////////////////////////////////////////////////////////////////

echo "<reaction_data" > $3
echo "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" >> $3
echo "  xmlns:xi=\"http://www.w3.org/2001/XInclude\"" >> $3
echo "  xsi:schemaLocation=\"http://libnucnet.sourceforge.net/xsd_pub/2014-12-13/libnucnet__reac/ http://libnucnet.sourceforge.net/xsd_pub/2014-12-13/libnucnet__reac.xsd\">" >> $3

for file in $2/*
do
  echo "  <xi:include href=\"${file##*/}\" xpointer=\"xpointer(//reaction)\"/>" >> $3
done
echo "</reaction_data>" >> $3

mv -f $3 $2

