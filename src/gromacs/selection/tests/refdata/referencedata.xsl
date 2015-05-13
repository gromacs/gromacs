<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="common-referencedata.xsl"/>

<!-- Index handling reference data -->

<xsl:template match="BlockAtoms">
    <xsl:if test="Sequence[@Name='Input']">
        <h2>Input Atoms</h2>
        <xsl:call-template name="SequenceAsHorizontalTable">
            <xsl:with-param name="root" select="Sequence[@Name='Input']"/>
        </xsl:call-template>
    </xsl:if>
    <h2>Blocks</h2>
    <table border="1">
        <tr>
            <th>Atom count</th>
            <th>Atoms</th>
        </tr>
        <xsl:for-each select="Block">
            <tr>
                <td><xsl:value-of select="Sequence[@Name='Atoms']/Int[@Name='Length']"/></td>
                <td>
                    <xsl:call-template name="SequenceAsCSV">
                        <xsl:with-param name="root" select="Sequence[@Name='Atoms']"/>
                    </xsl:call-template>
                </td>
            </tr>
        </xsl:for-each>
    </table>
</xsl:template>

<xsl:template match="IndexMapping">
    <h2><xsl:value-of select="@Name"/></h2>
    <h3>Input Atoms</h3>
    <xsl:call-template name="SequenceAsHorizontalTable">
        <xsl:with-param name="root" select="Sequence[@Name='Input']"/>
    </xsl:call-template>
    <h3>Mapping</h3>
    <table border="1">
        <tr>
            <th>RefId</th>
            <xsl:if test="Block/Int[@Name='MapId']">
                <th>MapId</th>
            </xsl:if>
            <th>Atom count</th>
            <th>Atoms</th>
        </tr>
        <xsl:for-each select="Block">
            <tr>
                <td><xsl:value-of select="Int[@Name='RefId']"/></td>
                <xsl:if test="Int[@Name='MapId']">
                    <td><xsl:value-of select="Int[@Name='MapId']"/></td>
                </xsl:if>
                <td><xsl:value-of select="Sequence[@Name='Atoms']/Int[@Name='Length']"/></td>
                <td>
                    <xsl:call-template name="SequenceAsCSV">
                        <xsl:with-param name="root" select="Sequence[@Name='Atoms']"/>
                    </xsl:call-template>
                </td>
            </tr>
        </xsl:for-each>
    </table>
</xsl:template>

<xsl:template match="OrgIdGroups">
    <h2>Groups: <xsl:value-of select="@Name"/></h2>
    <table>
        <tr>
            <td>Group count:</td>
            <td><xsl:value-of select="Int[@Name='GroupCount']"/></td>
        </tr>
        <tr>
            <td>OrgId</td>
            <td>
                <xsl:call-template name="SequenceAsCSV">
                    <xsl:with-param name="root" select="Sequence[@Name='OrgId']"/>
                </xsl:call-template>
            </td>
        </tr>
    </table>
</xsl:template>

<!-- Position calculation reference data -->

<xsl:template match="InitializedPositions">
    <h2>Initialized Positions</h2>
    <xsl:apply-templates />
</xsl:template>

<xsl:template match="EvaluatedPositions">
    <h2>Evaluated for <xsl:value-of select="@Name"/></h2>
    <xsl:apply-templates />
</xsl:template>

<xsl:template match="Positions">
    <xsl:if test="@Name">
        <h3><xsl:value-of select="@Name"/></h3>
    </xsl:if>
    <table>
        <tr>
            <td>Count:</td>
            <td>
                <xsl:value-of select="Int[@Name='Count']"/>
                (type: <xsl:value-of select="String[@Name='Type']"/>)
            </td>
        </tr>
        <tr>
            <td>Blocks:</td>
            <td>
                <xsl:call-template name="SequenceAsCSV">
                    <xsl:with-param name="root" select="Sequence[@Name='Block']"/>
                </xsl:call-template>
            </td>
        </tr>
    </table>
    <table border="1">
        <tr>
            <th>RefId</th>
            <th>Atom count</th>
            <th>Atoms</th>
            <xsl:if test="Position/Vector[@Name='Coordinates']">
                <th>Coordinates</th>
            </xsl:if>
            <xsl:if test="Position/Vector[@Name='Velocity']">
                <th>Velocity</th>
            </xsl:if>
            <xsl:if test="Position/Vector[@Name='Force']">
                <th>Force</th>
            </xsl:if>
        </tr>
        <xsl:for-each select="Position">
            <tr>
                <td><xsl:value-of select="Int[@Name='RefId']"/></td>
                <td><xsl:value-of select="Sequence[@Name='Atoms']/Int[@Name='Length']"/></td>
                <td>
                    <xsl:call-template name="SequenceAsCSV">
                        <xsl:with-param name="root" select="Sequence[@Name='Atoms']"/>
                    </xsl:call-template>
                </td>
                <xsl:if test="Vector[@Name='Coordinates']">
                    <td>
                        <xsl:apply-templates select="Vector[@Name='Coordinates']"/>
                    </td>
                </xsl:if>
                <xsl:if test="Vector[@Name='Velocity']">
                    <td>
                        <xsl:apply-templates select="Vector[@Name='Velocity']"/>
                    </td>
                </xsl:if>
                <xsl:if test="Vector[@Name='Force']">
                    <td>
                        <xsl:apply-templates select="Vector[@Name='Force']"/>
                    </td>
                </xsl:if>
            </tr>
        </xsl:for-each>
    </table>
</xsl:template>

<!-- Selection reference data -->

<xsl:key name="SelectionName" match="ParsedSelections/ParsedSelection" use="@Name"/>

<xsl:template match="ParsedSelections">
    <h2>Parsed Selections</h2>
    <table border="1">
        <tr>
            <th/>
            <th>Input</th>
            <xsl:if test="*/String[@Name='Name']">
                <th>Name</th>
            </xsl:if>
            <th>Text</th>
            <th>Dynamic</th>
        </tr>
        <xsl:for-each select="*">
        <tr>
            <td><xsl:value-of select="@Name"/></td>
            <td><xsl:value-of select="String[@Name='Input']"/></td>
            <xsl:if test="String[@Name='Name']">
                <td><xsl:value-of select="String[@Name='Name']"/></td>
            </xsl:if>
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
    <table>
        <xsl:if test="String[@Name='Name']">
            <tr>
                <td>Name:</td>
                <td><xsl:value-of select="String[@Name='Name']"/></td>
            </tr>
        </xsl:if>
        <tr>
            <td>Selection text:</td>
            <td>
                <xsl:value-of select="key('SelectionName', @Name)/String[@Name='Text']"/>
            </td>
        </tr>
        <xsl:if test="Sequence[@Name='Atoms']">
            <tr>
                <td>Atoms (<xsl:value-of select="Sequence[@Name='Atoms']/Int[@Name='Length']"/>):</td>
                <td>
                    <xsl:call-template name="SequenceAsCSV">
                        <xsl:with-param name="root" select="Sequence[@Name='Atoms']"/>
                    </xsl:call-template>
                </td>
            </tr>
        </xsl:if>
    </table>
    <xsl:apply-templates select="Sequence[@Name='Positions']"/>
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
                <xsl:if test="Position/Real[@Name='Mass']">
                    <th>Mass</th>
                </xsl:if>
                <xsl:if test="Position/Real[@Name='Charge']">
                    <th>Charge</th>
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
                <xsl:if test="Real[@Name='Mass']">
                    <td><xsl:value-of select="Real[@Name='Mass']"/></td>
                </xsl:if>
                <xsl:if test="Real[@Name='Charge']">
                    <td><xsl:value-of select="Real[@Name='Charge']"/></td>
                </xsl:if>
            </tr>
            </xsl:for-each>
        </table>
    </p>
</xsl:template>

</xsl:stylesheet>
