/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "cstringutil.h"

#include <stdio.h>
#include <stdlib.h>

#include <cassert>
#include <cctype>
#include <cstring>

#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

//! Comment sign to use.
#define COMMENTSIGN ';'

int continuing(char *s)
{
    int sl;
    assert(s);

    rtrim(s);
    sl = strlen(s);
    if ((sl > 0) && (s[sl-1] == CONTINUE))
    {
        s[sl-1] = 0;
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}



char *fgets2(char *line, int n, FILE *stream)
{
    char *c;
    if (fgets(line, n, stream) == nullptr)
    {
        return nullptr;
    }
    if ((c = strchr(line, '\n')) != nullptr)
    {
        *c = '\0';
    }
    else
    {
        /* A line not ending in a newline can only occur at the end of a file,
         * or because of n being too small.
         * Since both cases occur very infrequently, we can check for EOF.
         */
        if (!feof(stream))
        {
            gmx_fatal(FARGS, "An input file contains a line longer than %d characters, while the buffer passed to fgets2 has size %d. The line starts with: '%20.20s'", n, n, line);
        }
    }
    if ((c = strchr(line, '\r')) != nullptr)
    {
        *c = '\0';
    }

    return line;
}

void strip_comment (char *line)
{
    char *c;

    if (!line)
    {
        return;
    }

    /* search for a comment mark and replace it by a zero */
    if ((c = strchr(line, COMMENTSIGN)) != nullptr)
    {
        (*c) = 0;
    }
}

void upstring (char *str)
{
    int i;

    for (i = 0; (i < (int)strlen(str)); i++)
    {
        str[i] = toupper(str[i]);
    }
}

void lowstring (char *str)
{
    int i;
    for (i = 0; (i < (int)strlen(str)); i++)
        str[i] = tolower(str[i]);
}

void ltrim (char *str)
{
    int   i, c;

    if (nullptr == str)
    {
        return;
    }

    c = 0;
    while (('\0' != str[c]) && isspace(str[c]))
    {
        c++;
    }
    if (c > 0)
    {
        for (i = c; ('\0' != str[i]); i++)
        {
            str[i-c] = str[i];
        }
        str[i-c] = '\0';
    }
}

void rtrim (char *str)
{
    int nul;

    if (nullptr == str)
    {
        return;
    }

    nul = strlen(str)-1;
    while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t')) )
    {
        str[nul] = '\0';
        nul--;
    }
}

void trim (char *str)
{
    ltrim (str);
    rtrim (str);
}

int gmx_strcasecmp_min(const char *str1, const char *str2)
{
    char ch1, ch2;

    do
    {
        do
        {
            ch1 = toupper(*(str1++));
        }
        while ((ch1 == '-') || (ch1 == '_'));
        do
        {
            ch2 = toupper(*(str2++));
        }
        while ((ch2 == '-') || (ch2 == '_'));

        if (ch1 != ch2)
        {
            return (ch1-ch2);
        }
    }
    while (ch1);
    return 0;
}

int gmx_strncasecmp_min(const char *str1, const char *str2, int n)
{
    char  ch1, ch2;
    char *stri1, *stri2;

    stri1 = (char *)str1;
    stri2 = (char *)str2;
    do
    {
        do
        {
            ch1 = toupper(*(str1++));
        }
        while ((ch1 == '-') || (ch1 == '_'));
        do
        {
            ch2 = toupper(*(str2++));
        }
        while ((ch2 == '-') || (ch2 == '_'));

        if (ch1 != ch2)
        {
            return (ch1-ch2);
        }
    }
    while (ch1 && (str1-stri1 < n) && (str2-stri2 < n));
    return 0;
}

int gmx_strcasecmp(const char *str1, const char *str2)
{
    char ch1, ch2;

    do
    {
        ch1 = toupper(*(str1++));
        ch2 = toupper(*(str2++));
        if (ch1 != ch2)
        {
            return (ch1-ch2);
        }
    }
    while (ch1);
    return 0;
}

int gmx_strncasecmp(const char *str1, const char *str2, int n)
{
    char ch1, ch2;

    if (n == 0)
    {
        return 0;
    }

    do
    {
        ch1 = toupper(*(str1++));
        ch2 = toupper(*(str2++));
        if (ch1 != ch2)
        {
            return (ch1-ch2);
        }
        n--;
    }
    while (ch1 && n);
    return 0;
}

char *gmx_strdup(const char *src)
{
    char *dest;

    snew(dest, strlen(src)+1);
    strcpy(dest, src);

    return dest;
}

char *
gmx_strndup(const char *src, int n)
{
    int   len;
    char *dest;

    len = strlen(src);
    if (len > n)
    {
        len = n;
    }
    snew(dest, len+1);
    strncpy(dest, src, len);
    dest[len] = 0;
    return dest;
}

/* count the number of text elemets separated by whitespace in a string.
    str = the input string
    maxptr = the maximum number of allowed elements
    ptr = the output array of pointers to the first character of each element
    returns: the number of elements. */
int str_nelem(const char *str, int maxptr, char *ptr[])
{
    int   np = 0;
    char *copy0, *copy;

    copy0 = gmx_strdup(str);
    copy  = copy0;
    ltrim(copy);
    while (*copy != '\0')
    {
        if (np >= maxptr)
        {
            gmx_fatal(FARGS, "Too many groups on line: '%s' (max is %d)",
                      str, maxptr);
        }
        if (ptr)
        {
            ptr[np] = copy;
        }
        np++;
        while ((*copy != '\0') && !isspace(*copy))
        {
            copy++;
        }
        if (*copy != '\0')
        {
            *copy = '\0';
            copy++;
        }
        ltrim(copy);
    }
    if (ptr == nullptr)
    {
        sfree(copy0);
    }

    return np;
}

/* Magic hash init number for Dan J. Bernsteins algorithm.
 * Do NOT use any other value unless you really know what you are doing.
 */
const unsigned int
    gmx_string_hash_init = 5381;


unsigned int
gmx_string_fullhash_func(const char *s, unsigned int hash_init)
{
    int c;

    while ((c = (*s++)) != '\0')
    {
        hash_init = ((hash_init << 5) + hash_init) ^ c; /* (hash * 33) xor c */
    }
    return hash_init;
}

unsigned int
gmx_string_hash_func(const char *s, unsigned int hash_init)
{
    int c;

    while ((c = toupper(*s++)) != '\0')
    {
        if (isalnum(c))
        {
            hash_init = ((hash_init << 5) + hash_init) ^ c;            /* (hash * 33) xor c */
        }
    }
    return hash_init;
}

int
gmx_wcmatch(const char *pattern, const char *str)
{
    while (*pattern)
    {
        if (*pattern == '*')
        {
            /* Skip multiple wildcards in a sequence */
            while (*pattern == '*' || *pattern == '?')
            {
                ++pattern;
                /* For ?, we need to check that there are characters left
                 * in str. */
                if (*pattern == '?')
                {
                    if (*str == 0)
                    {
                        return GMX_NO_WCMATCH;
                    }
                    else
                    {
                        ++str;
                    }
                }
            }
            /* If the pattern ends after the star, we have a match */
            if (*pattern == 0)
            {
                return 0;
            }
            /* Match the rest against each possible suffix of str */
            while (*str)
            {
                /* Only do the recursive call if the first character
                 * matches. We don't have to worry about wildcards here,
                 * since we have processed them above. */
                if (*pattern == *str)
                {
                    int rc;
                    /* Match the suffix, and return if a match or an error */
                    rc = gmx_wcmatch(pattern, str);
                    if (rc != GMX_NO_WCMATCH)
                    {
                        return rc;
                    }
                }
                ++str;
            }
            /* If no suffix of str matches, we don't have a match */
            return GMX_NO_WCMATCH;
        }
        else if ((*pattern == '?' && *str != 0) || *pattern == *str)
        {
            ++str;
        }
        else
        {
            return GMX_NO_WCMATCH;
        }
        ++pattern;
    }
    /* When the pattern runs out, we have a match if the string has ended. */
    return (*str == 0) ? 0 : GMX_NO_WCMATCH;
}

char *wrap_lines(const char *buf, int line_width, int indent, gmx_bool bIndentFirst)
{
    char    *b2;
    int      i, i0, i2, j, b2len, lspace = 0, l2space = 0;
    gmx_bool bFirst, bFitsOnLine;

    /* characters are copied from buf to b2 with possible spaces changed
     * into newlines and extra space added for indentation.
     * i indexes buf (source buffer) and i2 indexes b2 (destination buffer)
     * i0 points to the beginning of the current line (in buf, source)
     * lspace and l2space point to the last space on the current line
     * bFirst is set to prevent indentation of first line
     * bFitsOnLine says if the first space occurred before line_width, if
     * that is not the case, we have a word longer than line_width which
     * will also not fit on the next line, so we might as well keep it on
     * the current line (where it also won't fit, but looks better)
     */

    b2    = nullptr;
    b2len = strlen(buf)+1+indent;
    snew(b2, b2len);
    i0 = i2 = 0;
    if (bIndentFirst)
    {
        for (i2 = 0; (i2 < indent); i2++)
        {
            b2[i2] = ' ';
        }
    }
    bFirst = TRUE;
    do
    {
        l2space = -1;
        /* find the last space before end of line */
        for (i = i0; ((i-i0 < line_width) || (l2space == -1)) && (buf[i]); i++)
        {
            b2[i2++] = buf[i];
            /* remember the position of a space */
            if (buf[i] == ' ')
            {
                lspace  = i;
                l2space = i2-1;
            }
            /* if we have a newline before the line is full, reset counters */
            if (buf[i] == '\n' && buf[i+1])
            {
                i0     = i+1;
                b2len += indent;
                srenew(b2, b2len);
                /* add indentation after the newline */
                for (j = 0; (j < indent); j++)
                {
                    b2[i2++] = ' ';
                }
            }
        }
        /* If we are at the last newline, copy it */
        if (buf[i] == '\n' && !buf[i+1])
        {
            b2[i2++] = buf[i++];
        }
        /* if we're not at the end of the string */
        if (buf[i])
        {
            /* check if one word does not fit on the line */
            bFitsOnLine = (i-i0 <= line_width);
            /* reset line counters to just after the space */
            i0 = lspace+1;
            i2 = l2space+1;
            /* if the words fit on the line, and we're beyond the indentation part */
            if ( (bFitsOnLine) && (l2space >= indent) )
            {
                /* start a new line */
                b2[l2space] = '\n';
                /* and add indentation */
                if (indent)
                {
                    if (bFirst)
                    {
                        line_width -= indent;
                        bFirst      = FALSE;
                    }
                    b2len += indent;
                    srenew(b2, b2len);
                    for (j = 0; (j < indent); j++)
                    {
                        b2[i2++] = ' ';
                    }
                    /* no extra spaces after indent; */
                    while (buf[i0] == ' ')
                    {
                        i0++;
                    }
                }
            }
        }
    }
    while (buf[i]);
    b2[i2] = '\0';

    return b2;
}

gmx_int64_t
str_to_int64_t(const char *str, char **endptr)
{
#ifndef _MSC_VER
    return strtoll(str, endptr, 10);
#else
    return _strtoi64(str, endptr, 10);
#endif
}

char *gmx_step_str(gmx_int64_t i, char *buf)
{
    sprintf(buf, "%" GMX_PRId64, i);
    return buf;
}
