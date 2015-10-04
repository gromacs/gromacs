/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * errors. It specifically targets .xvg files with a single ordinate
 * column (whose label, if present is that of the X axis), and at
 * least one abscissa column (whose labels are in the .xvg legend fields).
 *
 * The implementation repeatedly strips strings of leading and
 * trailing whitespace, because we do word-by-word parsing and need to
 * re-start the "parsing" knowing the first word is found at the start
 * of the string. */
XvgTable readXvgTable(TextReader *reader)
{
    XvgTable                   table;
    std::map<int, std::string> dataSetLegendStrings;

    Regex xAxisLabel("\\s*xaxis\\s+label.*\".*\"");
    Regex yAxisLabel("\\s*yaxis\\s+label.*\".*\"");
    Regex dataSetLegendLabelOne("legend string\\s+[0-9]+\\s+.*\".*\"");
    Regex dataSetLegendLabelTwo("s[0-9]+\\s+legend.*\".*\"");

    while (reader->readLine())
    {
        std::string line = stripString(reader->currentLine());
        if (line[0] == '#')
        {
            // Skip comments
        }
        else if (line[0] == '@')
        {
            // We found an .xvg command
            // TODO rebuild this code using proper regular expression submatching
            // when that is possible (and kill findQuotedString)
            std::string xvgCommand = stripString(line.substr(1));
            if (xvgCommand.compare(0, 8, "subtitle") == 0)
            {
                table.subtitle_ = findQuotedString(xvgCommand);
            }
            else if (regexMatch(xvgCommand, xAxisLabel))
            {
                table.ordinate_.legend_ = findQuotedString(xvgCommand);
            }
            else if (regexMatch(xvgCommand, yAxisLabel))
            {
                table.yAxisLabel_ = findQuotedString(xvgCommand);
            }
            // Maybe this is a legend string for a data set?
            // TODO remove this duplication with a regex matcher when that is possible
            else if (regexMatch(xvgCommand, dataSetLegendLabelOne))
            {
                int dataSet;
                int returnValue = sscanf(xvgCommand.c_str() + 13, "%d", &dataSet);
                if ((returnValue != EOF) && (returnValue != 0))
                {
                    if (dataSetLegendStrings.count(dataSet) > 0)
                    {
                        GMX_THROW(reader->makeError(formatString("Found second legend string for data set %lu", dataSet)));
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
                        GMX_THROW(reader->makeError(formatString("Found second legend string for data set %lu", dataSet)));
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
            std::vector<std::string> xvgColumns = splitString(line);

            if (xvgColumns.empty())
            {
                GMX_THROW(reader->makeError("Data contained zero columns"));
            }
            if (table.columns_.empty())
            {
                size_t numAbscissaColumns = xvgColumns.size() - 1;
                if (!dataSetLegendStrings.empty() && numAbscissaColumns != dataSetLegendStrings.size())
                {
                    GMX_THROW(reader->makeError(formatString("Found %lu data columns, but %lu were expected "
                                                             "from the legend strings plus the ordinate column",
                                                             xvgColumns.size(), dataSetLegendStrings.size())));
                }
                table.columns_.resize(numAbscissaColumns);
            }
            else if (table.columns_.size() + 1 != xvgColumns.size())
            {
                GMX_THROW(reader->makeError(formatString("From the first row of data, %lu columns of data "
                                                         "were expected, but this row had %lu columns",
                                                         table.columns_.size() + 1, xvgColumns.size())));
            }

            /* Read ordinate for this row and check for being a valid floating-point number */
            {
                double valueRead;
                int    returnValue = sscanf(xvgColumns[0].c_str(), "%lf", &valueRead);
                if ((returnValue == EOF) || (returnValue == 0))
                {
                    /* We've found some kind of invalid string. */
                    GMX_THROW(reader->makeError(formatString("Found invalid data %s", xvgColumns[0].c_str())));
                }
                table.ordinate_.rows_.push_back(valueRead);
            }
            /* Read each abscissa column for this row and check for being a valid floating-point number */
            for (size_t i = 1; i != xvgColumns.size(); ++i)
            {
                double valueRead;
                int    returnValue = sscanf(xvgColumns[i].c_str(), "%lf", &valueRead);
                if ((returnValue == EOF) || (returnValue == 0))
                {
                    /* We've found some kind of invalid string. */
                    GMX_THROW(reader->makeError(formatString("Found invalid numerical data '%s'", xvgColumns[i].c_str())));
                }
                table.columns_[i-1].rows_.push_back(valueRead);
            }
        }
    }

    for (size_t i = 0; i != table.columns_.size(); ++i)
    {
        if (dataSetLegendStrings.find(i) == dataSetLegendStrings.end())
        {
            // Note that by design of TextReader, there's no longer a
            // current line, so no context line is associated with
            // this message
            GMX_THROW(reader->makeError(formatString("A legend string is required for data set %lu (in column %lu)", i, i+1)));
        }
        table.columns_[i].legend_ = dataSetLegendStrings[i];
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
        ex.prependContext(formatString("Error reading .xvg file '%s'", filename));
        throw;
    }
}

XvgTable readXvgTableFromLegacyCode(const char *filename)
{
    try
    {
        return readXvgTable(filename);
    }
    // TODO Ideally, we could use prependContext similarly to the
    // above, but since this path is anyway deprecated, it's not a
    // priority.
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

} // namespace
