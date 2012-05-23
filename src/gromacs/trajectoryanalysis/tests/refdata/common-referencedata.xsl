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

<xsl:template match="/">
    <html><body>
        <xsl:apply-templates/>
    </body></html>
</xsl:template>

<xsl:template match="/ReferenceData">
    <h1>Test Reference Data</h1>
    <xsl:apply-templates/>
</xsl:template>

<xsl:template match="Vector">
    (<xsl:value-of select="*[@Name='X']"/>;
     <xsl:value-of select="*[@Name='Y']"/>;
     <xsl:value-of select="*[@Name='Z']"/>)
</xsl:template>

<xsl:template name="SequenceAsHorizontalTable">
    <xsl:param name="root" select="."/>
    <table border="1">
        <tr><th>Count</th><th>Items</th></tr>
        <tr>
            <td><xsl:value-of select="$root/Int[@Name='Length']"/></td>
            <td>
                <xsl:call-template name="SequenceAsCSV">
                    <xsl:with-param name="root" select="$root"/>
                </xsl:call-template>
            </td>
        </tr>
    </table>
</xsl:template>

<xsl:template name="SequenceAsCSV">
    <xsl:param name="root" select="."/>
    <xsl:for-each select="$root/*">
        <xsl:if test="not(.[@Name])">
            <xsl:apply-templates select="."/>
            <xsl:if test="position() &lt; last()">, </xsl:if>
        </xsl:if>
    </xsl:for-each>
</xsl:template>

<xsl:template name="Bool">
    <xsl:value-of select="."/>
</xsl:template>

<xsl:template name="String">
    <xsl:value-of select="."/>
</xsl:template>

<xsl:template name="Int">
    <xsl:value-of select="."/>
</xsl:template>

<xsl:template name="Real">
    <xsl:value-of select="."/>
</xsl:template>

</xsl:stylesheet>
