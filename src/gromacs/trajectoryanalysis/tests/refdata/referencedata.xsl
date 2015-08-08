<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="common-referencedata.xsl"/>
<xsl:import href="analysisdata-referencedata.xsl"/>

<xsl:template match="String[@Name='CommandLine']">
    <xsl:value-of select="."/>
</xsl:template>

<xsl:template match="OutputData">
    <h2>Raw Output Data</h2>
    <xsl:apply-templates />
</xsl:template>

<xsl:template match="AnalysisData">
    <h3><xsl:value-of select="@Name"/></h3>
    <xsl:apply-imports />
</xsl:template>

</xsl:stylesheet>
