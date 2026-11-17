<?xml version="1.0"?>
<!--
Combines the CTest results of several compute-sanitizer passes into one JUnit report.

The input is a list of the CTest result files of the passes, one per tool, e.g.
  <runs>
    <run tool="memcheck" href="Testing/ComputeSanitizer/memcheck/20261117-1036/DynamicAnalysis-Test.xml"/>
    ...
  </runs>
where each href is relative to the list file. The output has one <testsuite> per tool,
whose test cases are converted as in CTest2JUnit.xsl.
-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:import href="CTest2JUnit.xsl"/>
    <xsl:output method="xml" indent="yes"/>

    <xsl:template match="/runs">
        <testsuites>
            <xsl:for-each select="run">
                <xsl:variable name="site" select="document(@href, .)/Site"/>
                <testsuite name="{@tool}"
                           hostname="{$site/@Name}"
                           errors="0"
                           failures="{count($site/Testing/Test[@Status='failed'])}"
                           skipped="{count($site/Testing/Test[@Status='notrun'])}"
                           tests="{count($site/Testing/Test)}">
                    <xsl:apply-templates select="$site/Testing/Test"/>
                </testsuite>
            </xsl:for-each>
        </testsuites>
    </xsl:template>
</xsl:stylesheet>
