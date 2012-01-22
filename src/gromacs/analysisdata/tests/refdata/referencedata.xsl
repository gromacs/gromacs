<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="../../../../testutils/common-referencedata.xsl"/>

<xsl:template match="AnalysisData">
    <h2><xsl:value-of select="@Name"/></h2>
    <table border="1">
        <tr>
            <th>Frame</th>
            <th>X</th>
            <th>Values</th>
        </tr>
        <xsl:for-each select="DataFrame/Sequence[@Name='Y']">
        <tr>
            <td><xsl:value-of select="../@Name"/></td>
            <td><xsl:value-of select="../Real[@Name='X']"/></td>
            <td><xsl:call-template name="SequenceAsCSV"/></td>
        </tr>
        </xsl:for-each>
    </table>
</xsl:template>

<xsl:template match="DataValue">
    <xsl:value-of select="Real[@Name='Value']"/>
</xsl:template>

</xsl:stylesheet>
