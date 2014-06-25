<?xml version="1.0"?>

<!--
This file is currently duplicated to each directory containing reference data
XML files. This is to make it compatible with more browsers.
To keep these files in sync, please only modify the version in
  src/gromacs/analysisdata/tests/refdata/
and use the src/testutils/copy_xsl.sh script to copy it to relevant locations.
-->
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="common-referencedata.xsl"/>

<xsl:template match="AnalysisData">
    <xsl:variable name="has-datasetspec"
                  select="DataFrame/DataValues/Int[@Name='DataSet']"/>
    <xsl:variable name="has-columnspec"
                  select="DataFrame/DataValues/Int[@Name='FirstColumn']"/>
    <table border="1">
        <tr>
            <th>Frame</th>
            <th>X</th>
            <xsl:if test="$has-datasetspec">
                <th>Set</th>
            </xsl:if>
            <xsl:if test="$has-columnspec">
                <th>Columns</th>
            </xsl:if>
            <th>Values</th>
        </tr>
        <xsl:for-each select="DataFrame/DataValues">
        <tr>
            <td><xsl:value-of select="../@Name"/></td>
            <td><xsl:value-of select="../Real[@Name='X']"/></td>
            <xsl:if test="$has-datasetspec">
                <td><xsl:value-of select="Int[@Name='DataSet']"/></td>
            </xsl:if>
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

<xsl:template match="DataValue[Bool[@Name='Present']='false']">
    (
    <xsl:value-of select="Real[@Name='Value']"/>
    <xsl:if test="Real[@Name='Error']">
        &#177; <xsl:value-of select="Real[@Name='Error']"/>
    </xsl:if>
    )
</xsl:template>
<xsl:template match="DataValue">
    <xsl:value-of select="Real[@Name='Value']"/>
    <xsl:if test="Real[@Name='Error']">
        &#177; <xsl:value-of select="Real[@Name='Error']"/>
    </xsl:if>
</xsl:template>

</xsl:stylesheet>
