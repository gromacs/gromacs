<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="common-referencedata.xsl"/>
<xsl:import href="analysisdata-referencedata.xsl"/>

<xsl:template match="AnalysisData">
    <h2><xsl:value-of select="@Name"/></h2>
    <xsl:apply-imports />
</xsl:template>

</xsl:stylesheet>
