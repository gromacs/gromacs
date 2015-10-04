/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Definitions for xvgreader.h
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "xvgreader.h"

#include <map>
#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxregex.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

namespace
{

/*! \brief Return maximal substring of \c line enclosed by double quote marks.
 *
 * Nested quote marks are not supported
 *
 * \throws std::bad_alloc if out of memory
 */
std::string findQuotedString(const std::string &line)
{
    std::string toReturn("");
    size_t      firstQuote;

    firstQuote = line.find('"');
    if (firstQuote != std::string::npos)
    {
        size_t secondQuote;

        secondQuote = line.find('"', firstQuote + 1);
        if (secondQuote != std::string::npos)
        {
            toReturn = line.substr(firstQuote + 1, secondQuote - firstQuote - 1);
        }
    }

    return toReturn;
}

}   // namespace

namespace gmx
{

/* TODO this would probably be more readable if implemented using
 * C++11 regex sub-expression matching
 *
 * This implementation tracks the logic of read_xvg_legend in xvgr.cpp
 * fairly closely, but corrects it in a few places and handles more
 * errors.
 *
 * The implementation repeatedly strips strings of leading and
 * trailing whitespace, because we do word-by-word parsing and need to
 * re-start the "parsing" knowing the first word is found at the start
 * of the string.
 *
 * Some of the comments should be removed once read_xvg() is also
 * removed. */
XvgTable readXvgTable(TextReader *reader)
{
    XvgTable                   table;
    std::string                rawLine;
    std::map<int, std::string> dataSetLegendStrings;

    Regex xAxisLabel("\\s*xaxis\\s+label.*\".*\"");
    Regex yAxisLabel("\\s*yaxis\\s+label.*\".*\"");
    Regex dataSetLegendLabelOne("legend string\\s+[0-9]+\\s+.*\".*\"");
    Regex dataSetLegendLabelTwo("s[0-9]+\\s+legend.*\".*\"");

    while (reader->readLine(&rawLine))
    {
        std::string line = stripString(rawLine);
        if (line[0] == '#')
        {
            // Skip comments
        }
        else if (line[0] == '@')
        {
            // We found an .xvg command
            // TODO rebuild this code using proper regular expression submatching, and kill findQuotedString
            std::string xvgCommand = stripString(rawLine.substr(1));
            if (xvgCommand.compare(0, 8, "subtitle") == 0)
            {
                table.subtitle_ = findQuotedString(xvgCommand);
            }
            else if (regexMatch(xvgCommand, xAxisLabel))
            {
                table.xAxisLabel_ = findQuotedString(xvgCommand);
            }
            else if (regexMatch(xvgCommand, yAxisLabel))
            {
                table.yAxisLabel_ = findQuotedString(xvgCommand);
            }
            // Maybe this is a legend string for a data set
            // TODO remove this duplication with a regex matcher
            else if (regexMatch(xvgCommand, dataSetLegendLabelOne))
            {
                int dataSet;
                int returnValue = sscanf(xvgCommand.c_str() + 13, "%d", &dataSet);
                if ((returnValue != EOF) && (returnValue != 0))
                {
                    if (dataSetLegendStrings.count(dataSet) > 0)
                    {
                        GMX_THROW(reader->makeError("Found second legend string for data set %d"));
                    }
                    dataSetLegendStrings[dataSet] = findQuotedString(xvgCommand);
                }
            }
            else if (regexMatch(xvgCommand, dataSetLegendLabelTwo))
            {
                // This contains a legend for data set e.g. 's0 legend "legend string"'
                int         dataSet;
                const char *ptr         = xvgCommand.c_str() + 1;
                int         returnValue = sscanf(ptr, "%d", &dataSet);
                if ((returnValue != EOF) && (returnValue != 0))
                {
                    if (dataSetLegendStrings.count(dataSet) > 0)
                    {
                        GMX_THROW(reader->makeError("Found second legend string for data set %d"));
                    }
                    dataSetLegendStrings[dataSet] = findQuotedString(xvgCommand);
                }
            }
        }
        else if (line[0] == '&')
        {
            // End when we reach an .xvg data-set separator
            break;
        }
        else
        {
            // If we get here, this line must contain columns of data
            if (table.columns_.empty())
            {
                int numColumns = countWords(line);
                if (numColumns == 0)
                {
                    GMX_THROW(reader->makeError("Data contained zero columns"));
                }
                table.columns_.resize(numColumns);
            }

            size_t      numColumnsRead = 0;
            std::string format, baseFormat;
            /* Read each column for this row. Successive format
               strings are constructed so that they ignore the columns
               already read, and return the next data point for the
               next column in valueRead */
            for (auto &column : table.columns_)
            {
                format      = baseFormat + "%lf";
                baseFormat += "%*s";

                double valueRead;
                int    returnValue = sscanf(line.c_str(), format.c_str(), &valueRead);
                if ((returnValue == EOF) || (returnValue == 0))
                {
                    /* We've run out of columns for this
                       row. Previously we filled in blank columns with
                       zero, and put a message to stderr. This means a
                       user might ignore a real problem, because
                       normally they spammed by messages that are not
                       a problem. Instead, make the user give us a
                       well-formed .xvg file. */
                    GMX_THROW(reader->makeError(formatString("Found fewer than the expected %lu data sets", table.columns_.size())));
                }
                column.rows_.push_back(valueRead);
                numColumnsRead++;
            }
            if (!dataSetLegendStrings.empty() && numColumnsRead < dataSetLegendStrings.size())
            {
                GMX_THROW(reader->makeError(formatString("Found fewer data columns than expected for the %lu legend strings found", dataSetLegendStrings.size())));
            }
        }
    }

    if (!dataSetLegendStrings.empty())
    {
        size_t legendOffset = 0;
        if (dataSetLegendStrings.size() + 1 == table.columns_.size())
        {
            /* .xvg data sets are often graph ordinates (unless we're
               abusing .xvg for something else), and since there's a
               single legend string missing, there seems to be an
               abscissa column also. Probably a sensible label is the
               x-axis label. */
            table.columns_[legendOffset++].legend_ = table.xAxisLabel_;
        }
        for (auto &legend : dataSetLegendStrings)
        {
            int column = legend.first + legendOffset;
            /* Even if somehow a legend string was missing, std::map
               will make an empty string spring into existence, and
               that's the best thing we can do for a legend string. */
            table.columns_[column].legend_ = legend.second;
        }
    }

    return table;
}

XvgTable readXvgTable(const char *filename)
{
    try
    {
        TextReader fileText(filename);
        return readXvgTable(&fileText);
    }
    catch (InvalidInputError &ex)
    {
        ex.prependContext(formatString("Error reading .xvg file '%s'\n", filename));
        throw;
    }
}

XvgTable readXvgTableFromLegacyCode(const char *filename)
{
    try
    {
        return readXvgTable(filename);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

} // namespace
