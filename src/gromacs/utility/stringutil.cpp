/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements functions and classes in stringutil.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include "stringutil.h"

#include <cctype>
#include <cstdio>
#include <cstdarg>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

bool endsWith(const std::string &str, const char *suffix)
{
    if (suffix == NULL || suffix[0] == '\0')
    {
        return true;
    }
    size_t length = std::strlen(suffix);
    return (str.length() >= length &&
            str.compare(str.length() - length, length, suffix) == 0);
}

std::string stripSuffixIfPresent(const std::string &str, const char *suffix)
{
    if (suffix != NULL)
    {
        size_t suffixLength = std::strlen(suffix);
        if (suffixLength > 0 && endsWith(str, suffix))
        {
            return str.substr(0, str.length() - suffixLength);
        }
    }
    return str;
}

std::string formatString(const char *fmt, ...)
{
    va_list           ap;
    char              staticBuf[1024];
    int               length = 1024;
    std::vector<char> dynamicBuf;
    char             *buf    = staticBuf;

    // TODO: There may be a better way of doing this on Windows, Microsoft
    // provides their own way of doing things...
    while (1)
    {
        va_start(ap, fmt);
        int n = vsnprintf(buf, length, fmt, ap);
        va_end(ap);
        if (n > -1 && n < length)
        {
            std::string result(buf);
            return result;
        }
        if (n > -1)
        {
            length = n + 1;
        }
        else
        {
            length *= 2;
        }
        dynamicBuf.resize(length);
        buf = &dynamicBuf[0];
    }
}

std::string concatenateStrings(const char *const *sarray, size_t count)
{
    std::string result;

    for (size_t i = 0; i < count && sarray[i] != NULL; ++i)
    {
        if (sarray[i][0] != '\0')
        {
            result.append(sarray[i]);
            char lastchar = sarray[i][std::strlen(sarray[i]) - 1];
            if (!std::isspace(lastchar))
            {
                result.append(" ");
            }
        }
    }
    result.resize(result.find_last_not_of(" \n\r\t") + 1);
    return result;
}

namespace
{

/*! \brief
 * Common implementation for string replacement functions.
 *
 * \param[in] input  Input string.
 * \param[in] from   String to find.
 * \param[in] to     String to use to replace \p from.
 * \param[in] bWholeWords  Whether to only consider matches to whole words.
 * \returns   \p input with all occurrences of \p from replaced with \p to.
 * \throws    std::bad_alloc if out of memory.
 */
std::string
replaceInternal(const std::string &input, const char *from, const char *to,
                bool bWholeWords)
{
    GMX_RELEASE_ASSERT(from != NULL && to != NULL,
                       "Replacement strings must not be NULL");
    size_t      matchLength = std::strlen(from);
    std::string result;
    size_t      inputPos    = 0;
    size_t      matchPos    = input.find(from);
    while (matchPos < input.length())
    {
        size_t matchEnd = matchPos + matchLength;
        if (bWholeWords)
        {
            if (!((matchPos == 0 || !std::isalnum(input[matchPos - 1])) &&
                  (matchEnd == input.length() || !std::isalnum(input[matchEnd]))))
            {
                matchPos = input.find(from, matchPos + 1);
                continue;
            }

        }
        result.append(input, inputPos, matchPos - inputPos);
        result.append(to);
        inputPos = matchEnd;
        matchPos = input.find(from, inputPos);
    }
    result.append(input, inputPos, matchPos - inputPos);
    return result;
}

}   // namespace

std::string
replaceAll(const std::string &input, const char *from, const char *to)
{
    return replaceInternal(input, from, to, false);
}

std::string
replaceAllWords(const std::string &input, const char *from, const char *to)
{
    return replaceInternal(input, from, to, true);
}


/********************************************************************
 * TextLineWrapperSettings
 */

TextLineWrapperSettings::TextLineWrapperSettings()
    : maxLength_(0), indent_(0), firstLineIndent_(-1),
    bStripLeadingWhitespace_(true), continuationChar_('\0')
{
}


/********************************************************************
 * TextLineWrapper
 */

size_t
TextLineWrapper::findNextLine(const char *input, size_t lineStart) const
{
    size_t inputLength = std::strlen(input);
    bool   bFirstLine  = (lineStart == 0 || input[lineStart - 1] == '\n');
    // Ignore leading whitespace if necessary.
    if (!bFirstLine || settings_.bStripLeadingWhitespace_)
    {
        lineStart += std::strspn(input + lineStart, " ");
        if (lineStart >= inputLength)
        {
            return inputLength;
        }
    }

    int    indent = (bFirstLine ? settings_.firstLineIndent() : settings_.indent());
    size_t lastAllowedBreakPoint
        = (settings_.lineLength() > 0
           ? std::min(lineStart + settings_.lineLength() - indent, inputLength)
           : inputLength);
    // Ignore trailing whitespace.
    lastAllowedBreakPoint += std::strspn(input + lastAllowedBreakPoint, " ");
    size_t lineEnd = lineStart;
    do
    {
        const char *nextBreakPtr = std::strpbrk(input + lineEnd, " \n");
        size_t      nextBreak
            = (nextBreakPtr != NULL ? nextBreakPtr - input : inputLength);
        if (nextBreak > lastAllowedBreakPoint && lineEnd > lineStart)
        {
            break;
        }
        lineEnd = nextBreak + 1;
    }
    while (lineEnd < lastAllowedBreakPoint && input[lineEnd - 1] != '\n');
    return (lineEnd < inputLength ? lineEnd : inputLength);
}

size_t
TextLineWrapper::findNextLine(const std::string &input, size_t lineStart) const
{
    return findNextLine(input.c_str(), lineStart);
}

std::string
TextLineWrapper::formatLine(const std::string &input,
                            size_t lineStart, size_t lineEnd) const
{
    size_t inputLength = input.length();
    bool   bFirstLine  = (lineStart == 0 || input[lineStart - 1] == '\n');
    // Strip leading whitespace if necessary.
    if (!bFirstLine || settings_.bStripLeadingWhitespace_)
    {
        lineStart = input.find_first_not_of(' ', lineStart);
        if (lineStart >= inputLength)
        {
            return std::string();
        }
    }
    int  indent        = (bFirstLine ? settings_.firstLineIndent() : settings_.indent());
    bool bContinuation = (lineEnd < inputLength && input[lineEnd - 1] != '\n');
    // Strip trailing whitespace.
    while (lineEnd > lineStart && std::isspace(input[lineEnd - 1]))
    {
        --lineEnd;
    }

    size_t      lineLength = lineEnd - lineStart;
    std::string result(indent, ' ');
    result.append(input, lineStart, lineLength);
    if (bContinuation && settings_.continuationChar_ != '\0')
    {
        result.append(1, ' ');
        result.append(1, settings_.continuationChar_);
    }
    return result;
}

std::string
TextLineWrapper::wrapToString(const std::string &input) const
{
    std::string result;
    size_t      lineStart = 0;
    size_t      length    = input.length();
    while (lineStart < length)
    {
        size_t nextLineStart = findNextLine(input, lineStart);
        result.append(formatLine(input, lineStart, nextLineStart));
        if (nextLineStart < length ||
            (nextLineStart == length && input[length - 1] == '\n'))
        {
            result.append("\n");
        }
        lineStart = nextLineStart;
    }
    return result;
}

std::vector<std::string>
TextLineWrapper::wrapToVector(const std::string &input) const
{
    std::vector<std::string> result;
    size_t lineStart = 0;
    size_t length    = input.length();
    while (lineStart < length)
    {
        size_t nextLineStart = findNextLine(input, lineStart);
        result.push_back(formatLine(input, lineStart, nextLineStart));
        lineStart = nextLineStart;
    }
    return result;
}

} // namespace gmx
