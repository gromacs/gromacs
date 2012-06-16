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
        if (sarray[i][0] != '\0')
        {
            result.append(sarray[i]);
            char lastchar = sarray[i][std::strlen(sarray[i])-1];
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
 * TextLineWrapper::Impl
 */

/*! \internal \brief
 * Private implementation class for TextLineWrapper.
 *
 * \ingroup module_utility
 */
class TextLineWrapper::Impl
{
    public:
        //! Initialize default values for the wrapper.
        Impl() : maxLength_(0) {}

        /*! \brief
         * Helper method to find the next wrapped line.
         *
         * \param[in]     input         Full input string.
         * \param[in,out] lineStartPtr
         *      Index of first character for the line to wrap.
         *      On exit, index of the first character of the next line.
         * \returns   Next output line.
         */
        std::string wrapNextLine(const std::string &input,
                                 size_t *lineStartPtr) const;

        //! Maximum length of output lines, or <= 0 if no limit.
        int                     maxLength_;
};

std::string
TextLineWrapper::Impl::wrapNextLine(const std::string &input,
                                    size_t *lineStartPtr) const
{
    // Strip leading whitespace.
    size_t lineStart = input.find_first_not_of(' ', *lineStartPtr);
    if (lineStart == std::string::npos)
    {
        *lineStartPtr = lineStart;
        return std::string();
    }

    size_t lineEnd = std::string::npos;
    size_t nextNewline
        = std::min(input.find('\n', lineStart), input.length());
    if (maxLength_ <= 0 || nextNewline <= lineStart + maxLength_)
    {
        lineEnd = nextNewline;
    }
    else
    {
        size_t bestSpace = input.rfind(' ', lineStart + maxLength_);
        if (bestSpace < lineStart || bestSpace == std::string::npos)
        {
            bestSpace = input.find(' ', lineStart);
        }
        lineEnd = std::min(bestSpace, nextNewline);
    }

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
    return input.substr(lineStart, lineLength);
}

/********************************************************************
 * TextLineWrapper
 */

TextLineWrapper::TextLineWrapper()
    : impl_(new Impl)
{
}

TextLineWrapper::~TextLineWrapper()
{
}

TextLineWrapper &TextLineWrapper::setLineLength(int length)
{
    impl_->maxLength_ = length;
    return *this;
}

std::string
TextLineWrapper::wrapToString(const std::string &input) const
{
    std::string result;
    size_t lineStart = 0;
    size_t length = input.length();
    while (lineStart < length)
    {
        result.append(impl_->wrapNextLine(input, &lineStart));
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
        result.push_back(impl_->wrapNextLine(input, &lineStart));
    }
    return result;
}

} // namespace gmx
