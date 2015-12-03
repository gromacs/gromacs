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

<xsl:template match="OutputFiles">
    <xsl:if test="*/*">
        <h2>Output Files</h2>
        <xsl:apply-templates />
    </xsl:if>
</xsl:template>

<xsl:template match="OutputFiles/File">
    <xsl:if test="*">
        <h3><xsl:value-of select="@Name"/></h3>
        <xsl:apply-templates />
    </xsl:if>
</xsl:template>

<xsl:template match="OutputFiles/File/String[@Name='Contents']">
    <pre>
        <xsl:value-of select="substring(.,2)"/>
    </pre>
</xsl:template>

<xsl:template match="OutputFiles/File/XvgLegend/String[@Name='XvgLegend']">
    <pre>
        <xsl:value-of select="substring(.,2)"/>
    </pre>
</xsl:template>

<xsl:template match="OutputFiles/File/XvgData">
    <xsl:choose>
        <xsl:when test="*">
            <table>
                <xsl:apply-templates />
            </table>
        </xsl:when>
        <xsl:otherwise>Data omitted</xsl:otherwise>
    </xsl:choose>
</xsl:template>

<xsl:template match="OutputFiles/File/XvgData/Sequence">
    <tr>
        <xsl:apply-templates select="Real"/>
    </tr>
</xsl:template>

<xsl:template match="OutputFiles/File/XvgData/Sequence/Real">
    <td><xsl:value-of select="."/></td>
</xsl:template>

<xsl:template match="InteractiveSession">
    <pre>
        <xsl:for-each select="*">
            <xsl:choose>
                <xsl:when test="starts-with(@Name, 'Output')">
                    <xsl:value-of select="substring(.,2)"/>
                </xsl:when>
                <xsl:when test="string-length(.)=1">
                    <xsl:text>&#x25ba;</xsl:text>
                    <xsl:text>&#xb6;</xsl:text>
                </xsl:when>
                <xsl:when test="contains(substring(.,2), '&#10;')">
                    <xsl:text>&#x25ba;</xsl:text>
                    <xsl:value-of select="translate(substring(.,2), '&#10;', '&#x23ce;')"/>
                    <xsl:text>&#10;</xsl:text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:text>&#x25ba;</xsl:text>
                    <xsl:value-of select="substring(.,2)"/>
                    <xsl:text>&#xb6;</xsl:text>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:for-each>
        <xsl:text>[EOF]</xsl:text>
    </pre>
</xsl:template>

</xsl:stylesheet>
