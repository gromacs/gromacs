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
    return (str.length() >= length
            && str.compare(str.length() - length, length, suffix) == 0);
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
    va_list ap;
    char staticBuf[1024];
    int length = 1024;
    std::vector<char> dynamicBuf;
    char *buf = staticBuf;

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
        if (sarray[i][0] == '\0')
        {
            result.append("\n\n");
        }
        else
        {
            result.append(sarray[i]);
            char lastchar = sarray[i][std::strlen(sarray[i])-1];
            if (!std::isspace(lastchar))
            {
                result.append(" ");
            }
        }
    }
    // Normalize to not have any trailing whitespace.
    result.resize(result.find_last_not_of(" \n\r\t") + 1);
    return result;
}

namespace
{

std::string
replaceInternal(const std::string &input, const char *from, const char *to,
                bool bWholeWords)
{
    GMX_RELEASE_ASSERT(from != NULL && to != NULL,
                       "Replacement strings must not be NULL");
    size_t matchLength = std::strlen(from);
    std::string result;
    size_t inputPos = 0;
    size_t matchPos = input.find(from);
    while (matchPos < input.length())
    {
        size_t matchEnd = matchPos + matchLength;
        if (bWholeWords)
        {
            if (!((matchPos == 0 || !std::isalnum(input[matchPos-1]))
                  && (matchEnd == input.length() || !std::isalnum(input[matchEnd]))))
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

} // namespace

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
      bStripLeadingWhitespace_(false), continuationChar_('\0')
{
}


/********************************************************************
 * TextLineWrapper::Impl
 */

/*! \internal \brief
 * Private implementation class for TextLineWrapper.
 *
 * \todo
 * It is a bit confusing to have this implementation class with only a static
 * method.  This problem may solve itself, though, if this incremental wrapping
 * method also needs to be exposed through TextLineWrapper.
 *
 * \ingroup module_utility
 */
class TextLineWrapper::Impl
{
    public:
        /*! \brief
         * Helper method to find the next wrapped line.
         *
         * \param[in]     settings      Wrapping settings.
         * \param[in]     input         Full input string.
         * \param[in,out] lineStartPtr
         *      Index of first character for the line to wrap.
         *      On exit, index of the first character of the next line.
         * \returns   Next output line.
         */
        static std::string
        wrapNextLine(const TextLineWrapperSettings &settings,
                     const std::string &input, size_t *lineStartPtr);
};

std::string
TextLineWrapper::Impl::wrapNextLine(const TextLineWrapperSettings &settings,
                                    const std::string &input,
                                    size_t *lineStartPtr)
{
    size_t lineStart = *lineStartPtr;
    bool bFirstLine = (lineStart == 0 || input[lineStart - 1] == '\n');
    // Strip leading whitespace.
    if (!bFirstLine || settings.bStripLeadingWhitespace_)
    {
        lineStart = input.find_first_not_of(' ', lineStart);
        if (lineStart == std::string::npos)
        {
            *lineStartPtr = lineStart;
            return std::string();
        }
    }

    int indent = settings.indent_;
    if (bFirstLine && settings.firstLineIndent_ >= 0)
    {
        indent = settings.firstLineIndent_;
    }
    size_t lineEnd = std::string::npos;
    size_t nextNewline
        = std::min(input.find('\n', lineStart), input.length());
    size_t lastAllowedBreakPoint = lineStart + settings.maxLength_ - indent;
    lastAllowedBreakPoint = input.find_first_not_of(' ', lastAllowedBreakPoint);
    if (settings.maxLength_ <= 0 || nextNewline <= lastAllowedBreakPoint)
    {
        lineEnd = nextNewline;
    }
    else
    {
        size_t bestSpace = input.rfind(' ', lastAllowedBreakPoint);
        if (bestSpace < lineStart || bestSpace == std::string::npos)
        {
            bestSpace = input.find(' ', lineStart);
        }
        lineEnd = std::min(bestSpace, nextNewline);
    }
    bool bContinuation = (lineEnd != std::string::npos && lineEnd != nextNewline);

    if (lineEnd == std::string::npos)
    {
        lineEnd = input.length();
    }
    *lineStartPtr = lineEnd + 1;
    // Strip trailing whitespace.
    while (lineEnd > lineStart && std::isspace(input[lineEnd - 1]))
    {
        --lineEnd;
    }

    size_t lineLength = lineEnd - lineStart;
    std::string result(indent, ' ');
    result.append(input, lineStart, lineLength);
    if (bContinuation && settings.continuationChar_ != '\0')
    {
        result.append(1, ' ');
        result.append(1, settings.continuationChar_);
    }
    return result;
}

/********************************************************************
 * TextLineWrapper
 */

std::string
TextLineWrapper::wrapToString(const std::string &input) const
{
    std::string result;
    size_t lineStart = 0;
    size_t length = input.length();
    while (lineStart < length)
    {
        result.append(Impl::wrapNextLine(settings_, input, &lineStart));
        if (lineStart < length
            || (lineStart == length && input[length - 1] == '\n'))
        {
            result.append("\n");
        }
    }
    return result;
}

std::vector<std::string>
TextLineWrapper::wrapToVector(const std::string &input) const
{
    std::vector<std::string> result;
    size_t lineStart = 0;
    while (lineStart < input.length())
    {
        result.push_back(Impl::wrapNextLine(settings_, input, &lineStart));
    }
    return result;
}

} // namespace gmx
