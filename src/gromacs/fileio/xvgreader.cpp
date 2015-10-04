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
/*! \libinternal \file
 * \brief
 * Definitions for xvgreader.h
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#include "gmxpre.h"

#include "xvgreader.h"

#include <map>
#include <string>
#include <vector>

#include "gromacs/fileio/errormessagemaker.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"

namespace
{

/*! \brief Return substring of \c line beginning at least at \c
 * initialOffset enclosed by double quote marks.
 *
 * Nested quote marks are not supported
 *
 * \throws std::bad_alloc if out of memory
 *
 * \todo Use C++11 regex */
std::string findQuotedString(const std::string &line,
                             size_t             initialOffset = 0)
{
    std::string toReturn("");
    size_t      firstQuote;

    firstQuote = line.find('"', initialOffset);
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
XvgTable readXvgTable(TextInputStream *stream)
{
    XvgTable                   table;
    std::string                rawLine;
    std::map<int, std::string> dataSetLegendStrings;
    ErrorMessageMaker          errorMessage;

    while (stream->readLine(&rawLine))
    {
        errorMessage.setCurrentLine(&rawLine);
        std::string line = stripString(rawLine);
        if (line[0] == '&')
        {
            // Blindly skip .xvg data set separators
            continue;
        }
        if (line[0] == '@')
        {
            // We found an .xvg command
            std::string xvgCommand = stripString(rawLine.substr(1));
            int         dataSet    = -1;
            size_t      offset     = 0;
            if (xvgCommand.compare(offset, 8, "subtitle") == 0)
            {
                offset         += 8;
                table.subtitle_ = findQuotedString(xvgCommand, offset);
            }
            else if (xvgCommand.compare(offset, 12, "xaxis  label") == 0)
            {
                offset           += 12;
                table.xAxisLabel_ = findQuotedString(xvgCommand, offset);
            }
            else if (xvgCommand.compare(offset, 12, "yaxis  label") == 0)
            {
                offset           += 12;
                table.yAxisLabel_ = findQuotedString(xvgCommand, offset);
            }
            else if (xvgCommand.compare(offset, 13, "legend string") == 0)
            {
                offset += 13;
                int nchar;
                int returnValue = sscanf(xvgCommand.c_str() + offset, "%d%n", &dataSet, &nchar);
                if ((returnValue == EOF) || (returnValue == 0))
                {
                    GMX_THROW(errorMessage.make("Cannot parse .xvg"));
                }
                offset += nchar;
            }
            else if (xvgCommand.length() > 1 && xvgCommand[0] == 's')
            {
                // This might contain a legend for data set e.g. 's0'
                offset++;
                int         nchar;
                const char *ptr         = xvgCommand.c_str() + offset;
                int         returnValue = sscanf(ptr, "%d%n", &dataSet, &nchar);
                if ((returnValue == EOF) || (returnValue == 0))
                {
                    GMX_THROW(errorMessage.make("Cannot parse .xvg file"));
                }
                offset += nchar;
                if (xvgCommand.length() > offset)
                {
                    std::string detail = stripString(xvgCommand.substr(offset, std::string::npos));
                    if (detail.compare(0, 6, "legend") == 0)
                    {
                        offset += 6;
                    }
                    else
                    {
                        dataSet = -1;
                    }
                }
            }
            if (dataSet >= 0)
            {
                // We found a legend string for a data set, in one of various ways
                if (dataSetLegendStrings.count(dataSet) > 0)
                {
                    GMX_THROW(errorMessage.make("Found second legend string for data set %d"));
                }
                dataSetLegendStrings[dataSet] = findQuotedString(xvgCommand, offset);
            }
        }
        else if (line[0] != '#')
        {
            // This is not a comment line, so must contain columns of data
            if (table.columns_.empty())
            {
                int numColumns = countWords(line);
                if (numColumns == 0)
                {
                    GMX_THROW(errorMessage.make("Data contained zero columns"));
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
                    GMX_THROW(errorMessage.make(formatString("Found fewer than the expected %lu data sets", table.columns_.size())));
                }
                column.rows_.push_back(valueRead);
                numColumnsRead++;
            }
            if (!dataSetLegendStrings.empty() && numColumnsRead < dataSetLegendStrings.size())
            {
                GMX_THROW(errorMessage.make(formatString("Found fewer data columns than expected for the %lu legend strings found", dataSetLegendStrings.size())));
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
    TextInputFile fileStream(filename);
    try
    {
        return readXvgTable(&fileStream);
    }
    catch (InvalidInputError &ex)
    {
        ex.prependContext(formatString("Error reading .xvg file '%s'\n", filename));
        throw ex;
    }
    catch (const std::exception &ex)
    {
        throw ex;
    }
}

XvgTable readXvgTableFromLegacyCode(const char *filename) noexcept
{
    try
    {
        return readXvgTable(filename);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

} // namespace
