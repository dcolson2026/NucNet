<!--/////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2012 Clemson University.
// 
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
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
// 
/////////////////////////////////////////////////////////////////////////////-->

<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:param name="objdir"/>

<xsl:param name="codes"/>

<xsl:output method="text"/>

<xsl:template match="/">
  <xsl:if test="$codes = 'examples'">
    <xsl:apply-templates select="//example"/>
  </xsl:if>
  <xsl:if test="$codes = 'aux'">
    <xsl:apply-templates select="//other_files/file[contains(.,'.c')]"/>
  </xsl:if>
</xsl:template>

<xsl:template match="example">
  <xsl:value-of select="@name"/>
  <xsl:text> </xsl:text>
</xsl:template>

<xsl:template match="other_files/file[contains(.,'.c')]">
  <xsl:value-of select="$objdir"/>
  <xsl:text>/</xsl:text>
  <xsl:value-of select="substring-before(.,'.c')"/><xsl:text>.o</xsl:text>
  <xsl:text> </xsl:text>
</xsl:template>

</xsl:stylesheet>
