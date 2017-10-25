/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements functions and classes in stringutil.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "stringutil.h"

#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

std::size_t
countWords(const char *s)
{
    std::size_t nWords = 0;
    // Use length variable to avoid N^2 complexity when executing strlen(s) every iteration
    std::size_t length = std::strlen(s);

    for (std::size_t i = 0; i < length; i++)
    {
        // If we found a new word, increase counter and step through the word
        if (std::isalnum(s[i]))
        {
            ++nWords;
            // If we hit string end, '\0' is not alphanumerical
            while (std::isalnum(s[i]))
            {
                // This might increment i to the string end, and then the outer
                // loop will increment i one unit beyond that, but since
                // we compare to the string length in the outer loop this is fine.
                i++;
            }
        }
    }
    return nWords;
}


std::size_t
countWords(const std::string &str)
{
    // Under out beautiful C++ interface hides an ugly c-string implementation :-)
    return countWords(str.c_str());
}

bool endsWith(const char *str, const char *suffix)
{
    if (isNullOrEmpty(suffix))
    {
        return true;
    }
    const size_t strLength    = std::strlen(str);
    const size_t suffixLength = std::strlen(suffix);
    return (strLength >= suffixLength
            && std::strcmp(&str[strLength - suffixLength], suffix) == 0);
}

std::string stripSuffixIfPresent(const std::string &str, const char *suffix)
{
    if (suffix != nullptr)
    {
        size_t suffixLength = std::strlen(suffix);
        if (suffixLength > 0 && endsWith(str, suffix))
        {
            return str.substr(0, str.length() - suffixLength);
        }
    }
    return str;
}

std::string stripString(const std::string &str)
{
    std::string::const_iterator start = str.begin();
    std::string::const_iterator end   = str.end();
    while (start != end && std::isspace(*start))
    {
        ++start;
    }
    while (start != end && std::isspace(*(end - 1)))
    {
        --end;
    }
    return std::string(start, end);
}

std::string formatString(const char *fmt, ...)
{
    va_list     ap;
    va_start(ap, fmt);
    std::string result = formatStringV(fmt, ap);
    va_end(ap);
    return result;
}

std::string formatStringV(const char *fmt, va_list ap)
{
    va_list           ap_copy;
    char              staticBuf[1024];
    int               length = 1024;
    std::vector<char> dynamicBuf;
    char             *buf = staticBuf;

    // TODO: There may be a better way of doing this on Windows, Microsoft
    // provides their own way of doing things...
    while (1)
    {
        va_copy(ap_copy, ap);
        int n = vsnprintf(buf, length, fmt, ap_copy);
        va_end(ap_copy);
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
        buf = dynamicBuf.data();
    }
}

std::vector<std::string> splitString(const std::string &str)
{
    std::vector<std::string>          result;
    std::string::const_iterator       currPos = str.begin();
    const std::string::const_iterator end     = str.end();
    while (currPos != end)
    {
        while (currPos != end && std::isspace(*currPos))
        {
            ++currPos;
        }
        const std::string::const_iterator startPos = currPos;
        while (currPos != end && !std::isspace(*currPos))
        {
            ++currPos;
        }
        if (startPos != end)
        {
            result.emplace_back(startPos, currPos);
        }
    }
    return result;
}

std::vector<std::string> splitDelimitedString(const std::string &str, char delim)
{
    std::vector<std::string> result;
    size_t                   currPos = 0;
    const size_t             len     = str.length();
    if (len > 0)
    {
        size_t nextDelim;
        do
        {
            nextDelim = str.find(delim, currPos);
            result.push_back(str.substr(currPos, nextDelim - currPos));
            currPos = nextDelim < len ? nextDelim + 1 : len;
        }
        while (currPos < len || nextDelim < len);
    }
    return result;
}

std::vector<std::string> splitAndTrimDelimitedString(const std::string &str, char delim)
{
    std::vector<std::string> result;

    result = splitDelimitedString(str, delim);
    std::transform(result.begin(), result.end(), result.begin(), stripString);
    return result;
}

namespace
{

/*! \brief
 * Helper function to identify word boundaries for replaceAllWords().
 *
 * \returns  `true` if the character is considered part of a word.
 *
 * \ingroup module_utility
 */
bool isWordChar(char c)
{
    return std::isalnum(c) || c == '-' || c == '_';
}

/*! \brief
 * Common implementation for string replacement functions.
 *
 * \param[in] input  Input string.
 * \param[in] from   String to find.
 * \param[in] to     String to use to replace \p from.
 * \param[in] bWholeWords  Whether to only consider matches to whole words.
 * \returns   \p input with all occurrences of \p from replaced with \p to.
 * \throws    std::bad_alloc if out of memory.
 *
 * \ingroup module_utility
 */
std::string
replaceInternal(const std::string &input, const char *from, const char *to,
                bool bWholeWords)
{
    GMX_RELEASE_ASSERT(from != nullptr && to != nullptr,
                       "Replacement strings must not be NULL");
    size_t      matchLength = std::strlen(from);
    std::string result;
    size_t      inputPos = 0;
    size_t      matchPos = input.find(from);
    while (matchPos < input.length())
    {
        size_t matchEnd = matchPos + matchLength;
        if (bWholeWords)
        {
            if (!((matchPos == 0 || !isWordChar(input[matchPos-1]))
                  && (matchEnd == input.length() || !isWordChar(input[matchEnd]))))
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
replaceAll(const std::string &input, const std::string &from,
           const std::string &to)
{
    return replaceInternal(input, from.c_str(), to.c_str(), false);
}

std::string
replaceAllWords(const std::string &input, const char *from, const char *to)
{
    return replaceInternal(input, from, to, true);
}

std::string
replaceAllWords(const std::string &input, const std::string &from,
                const std::string &to)
{
    return replaceInternal(input, from.c_str(), to.c_str(), true);
}


/********************************************************************
 * TextLineWrapperSettings
 */

TextLineWrapperSettings::TextLineWrapperSettings()
    : maxLength_(0), indent_(0), firstLineIndent_(-1),
      bKeepFinalSpaces_(false), continuationChar_('\0')
{
}


/********************************************************************
 * TextLineWrapper
 */

bool TextLineWrapper::isTrivial() const
{
    return settings_.lineLength() == 0 && settings_.indent() == 0
           && settings_.firstLineIndent_ <= 0;
}

size_t
TextLineWrapper::findNextLine(const char *input, size_t lineStart) const
{
    size_t inputLength = std::strlen(input);
    bool   bFirstLine  = (lineStart == 0 || input[lineStart - 1] == '\n');
    // Ignore leading whitespace if necessary.
    if (!bFirstLine)
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
            = (nextBreakPtr != nullptr ? nextBreakPtr - input : inputLength);
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
    if (!bFirstLine)
    {
        lineStart = input.find_first_not_of(' ', lineStart);
        if (lineStart >= inputLength)
        {
            return std::string();
        }
    }
    int  indent        = (bFirstLine ? settings_.firstLineIndent() : settings_.indent());
    bool bContinuation = (lineEnd < inputLength && input[lineEnd - 1] != '\n');
    // Remove explicit line breaks in input
    // (the returned line should not contain line breaks).
    while (lineEnd > lineStart && input[lineEnd - 1] == '\n')
    {
        --lineEnd;
    }
    // Strip trailing whitespace, unless they are explicit in the input and it
    // has been requested to keep them.
    if (bContinuation || !settings_.bKeepFinalSpaces_)
    {
        while (lineEnd > lineStart && std::isspace(input[lineEnd - 1]))
        {
            --lineEnd;
        }
    }

    const size_t lineLength = lineEnd - lineStart;
    if (lineLength == 0)
    {
        return std::string();
    }
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
        if (nextLineStart < length
            || (nextLineStart == length && input[length - 1] == '\n'))
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
    size_t                   lineStart = 0;
    size_t                   length    = input.length();
    while (lineStart < length)
    {
        size_t nextLineStart = findNextLine(input, lineStart);
        result.push_back(formatLine(input, lineStart, nextLineStart));
        lineStart = nextLineStart;
    }
    return result;
}

} // namespace gmx
