<?xml version="1.0"?>

<!--
This file is currently duplicated to each directory containing reference data
XML files. This is to make it compatible with more browsers.
To keep these files in sync, please only modify the version in
  src/testutils/
and use the copy_xsl.sh script to copy it to relevant locations.
-->
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="common-referencedata.xsl"/>

<xsl:template match="AnalysisData">
    <xsl:variable name="has-columnspec"
                  select="DataFrame/Sequence[@Name='Y']/Int[@Name='FirstColumn']"/>
    <table border="1">
        <tr>
            <th>Frame</th>
            <th>X</th>
            <xsl:if test="$has-columnspec">
                <th>Columns</th>
            </xsl:if>
            <th>Values</th>
        </tr>
        <xsl:for-each select="DataFrame/Sequence[@Name='Y']">
        <tr>
            <td><xsl:value-of select="../@Name"/></td>
            <td><xsl:value-of select="../Real[@Name='X']"/></td>
            <xsl:if test="$has-columnspec">
                <td>
                    <xsl:choose>
                        <xsl:when test="Int[@Name='FirstColumn']">
                            <xsl:value-of select="Int[@Name='FirstColumn']"/>
                            <xsl:text>-</xsl:text>
                            <xsl:value-of select="Int[@Name='LastColumn']"/>
                        </xsl:when>
                        <xsl:otherwise>all</xsl:otherwise>
                    </xsl:choose>
                </td>
            </xsl:if>
            <td><xsl:call-template name="SequenceAsCSV"/></td>
        </tr>
        </xsl:for-each>
    </table>
</xsl:template>

<xsl:template match="DataValue">
    <xsl:value-of select="Real[@Name='Value']"/>
</xsl:template>

</xsl:stylesheet>
