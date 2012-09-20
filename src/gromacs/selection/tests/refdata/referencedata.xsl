<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="common-referencedata.xsl"/>

<xsl:key name="SelectionName" match="ParsedSelections/ParsedSelection" use="@Name"/>

<xsl:template match="ParsedSelections">
    <h2>Parsed Selections</h2>
    <table border="1">
        <tr>
            <th/>
            <th>Input</th>
            <th>Name</th>
            <th>Text</th>
            <th>Dynamic</th>
        </tr>
        <xsl:for-each select="*">
        <tr>
            <td><xsl:value-of select="@Name"/></td>
            <td><xsl:value-of select="String[@Name='Input']"/></td>
            <td><xsl:value-of select="String[@Name='Name']"/></td>
            <td><xsl:value-of select="String[@Name='Text']"/></td>
            <td><xsl:value-of select="Bool[@Name='Dynamic']"/></td>
        </tr>
        </xsl:for-each>
    </table>
</xsl:template>

<xsl:template match="CompiledSelections">
    <h2>Compiled Selections</h2>
    <xsl:apply-templates />
</xsl:template>

<xsl:template match="EvaluatedSelections">
    <h2>Evaluated for <xsl:value-of select="@Name"/></h2>
    <xsl:apply-templates />
</xsl:template>

<xsl:template match="Selection">
    <h3><xsl:value-of select="@Name"/></h3>
    <p>
        Selection text:<br/>
        <xsl:value-of select="key('SelectionName', @Name)/String[@Name='Text']"/>
    </p>
    <xsl:apply-templates />
</xsl:template>

<xsl:template match="Selection/Sequence[@Name='Atoms']">
    <p>
        Atoms:
        <xsl:call-template name="SequenceAsHorizontalTable"/>
    </p>
</xsl:template>

<xsl:template match="Selection/Sequence[@Name='Positions']">
    <p>
        Positions (count: <xsl:value-of select="Int[@Name='Length']"/>):
        <table border="1">
            <tr>
                <xsl:if test="Position/Sequence[@Name='Atoms']">
                    <th>Atom count</th>
                    <th>Atoms</th>
                </xsl:if>
                <xsl:if test="Position/Int[@Name='RefId']">
                    <th>RefId</th>
                    <th>MappedId</th>
                </xsl:if>
                <xsl:if test="Position/Vector[@Name='Coordinates']">
                    <th>Coordinates</th>
                </xsl:if>
            </tr>
            <xsl:for-each select="Position">
            <tr>
                <xsl:if test="Sequence[@Name='Atoms']">
                    <td><xsl:value-of select="Sequence[@Name='Atoms']/Int[@Name='Length']"/></td>
                    <td>
                        <xsl:call-template name="SequenceAsCSV">
                            <xsl:with-param name="root" select="Sequence[@Name='Atoms']"/>
                        </xsl:call-template>
                    </td>
                </xsl:if>
                <xsl:if test="Int[@Name='RefId']">
                    <td><xsl:value-of select="Int[@Name='RefId']"/></td>
                    <td><xsl:value-of select="Int[@Name='MappedId']"/></td>
                </xsl:if>
                <xsl:if test="Vector[@Name='Coordinates']">
                    <td>
                        <xsl:apply-templates select="Vector[@Name='Coordinates']"/>
                    </td>
                </xsl:if>
            </tr>
            </xsl:for-each>
        </table>
    </p>
</xsl:template>

</xsl:stylesheet>
