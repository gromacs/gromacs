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
#include "gmxpre.h"

#include "gmxcpp.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    char *name;
    char *def;
} t_define;

static int        ndef   = 0;
static t_define  *defs   = NULL;
static int        nincl  = 0;
static char     **incl   = 0;

/* enum used for handling ifdefs */
enum {
    eifTRUE, eifFALSE, eifIGNORE, eifNR
};

typedef struct gmx_cpp {
    FILE             *fp;
    char             *path, *cwd;
    char             *fn;
    int               line_len;
    char             *line;
    int               line_nr;
    int               nifdef;
    int              *ifdefs;
    struct   gmx_cpp *child, *parent;
} gmx_cpp;

static gmx_bool is_word_end(char c)
{
    return !(isalnum(c) || c == '_');
}

static const char *strstrw(const char *buf, const char *word)
{
    const char *ptr;

    while ((ptr = strstr(buf, word)) != NULL)
    {
        /* Check if we did not find part of a longer word */
        if (ptr &&
            is_word_end(ptr[strlen(word)]) &&
            (((ptr > buf) && is_word_end(ptr[-1])) || (ptr == buf)))
        {
            return ptr;
        }

        buf = ptr + strlen(word);
    }
    return NULL;
}

static gmx_bool find_directive(char *buf, char **name, char **val)
{
    /* Skip initial whitespace */
    while (isspace(*buf))
    {
        ++buf;
    }
    /* Check if this is a directive */
    if (*buf != '#')
    {
        return FALSE;
    }
    /* Skip the hash and any space after it */
    ++buf;
    while (isspace(*buf))
    {
        ++buf;
    }
    /* Set the name pointer and find the next space */
    *name = buf;
    while (*buf != 0 && !isspace(*buf))
    {
        ++buf;
    }
    /* Set the end of the name here, and skip any space */
    if (*buf != 0)
    {
        *buf = 0;
        ++buf;
        while (isspace(*buf))
        {
            ++buf;
        }
    }
    /* Check if anything is remaining */
    *val = (*buf != 0) ? buf : NULL;
    return TRUE;
}

static gmx_bool is_ifdeffed_out(gmx_cpp_t handle)
{
    return ((handle->nifdef > 0) && (handle->ifdefs[handle->nifdef-1] != eifTRUE));
}

static void add_include(const char *include)
{
    int i;

    if (include == NULL)
    {
        return;
    }

    for (i = 0; (i < nincl); i++)
    {
        if (strcmp(incl[i], include) == 0)
        {
            break;
        }
    }
    if (i == nincl)
    {
        nincl++;
        srenew(incl, nincl);
        incl[nincl-1] = gmx_strdup(include);
    }
}

static void done_includes()
{
    int i;
    for (i = 0; (i < nincl); i++)
    {
        sfree(incl[i]);
    }
    sfree(incl);
    incl  = NULL;
    nincl = 0;
}

static void add_define(const char *name, const char *value)
{
    int  i;

    for (i = 0; (i < ndef); i++)
    {
        if (strcmp(defs[i].name, name) == 0)
        {
            break;
        }
    }
    if (i == ndef)
    {
        ndef++;
        srenew(defs, ndef);
        i            = ndef - 1;
        defs[i].name = gmx_strdup(name);
    }
    else if (defs[i].def)
    {
        if (debug)
        {
            fprintf(debug, "Overriding define %s\n", name);
        }
        sfree(defs[i].def);
    }
    if (value && strlen(value) > 0)
    {
        defs[i].def  = gmx_strdup(value);
    }
    else
    {
        defs[i].def  = NULL;
    }
}

static void done_defines()
{
    int i;
    for (i = 0; (i < ndef); i++)
    {
        sfree(defs[i].name);
        sfree(defs[i].def);
    }
    sfree(defs);
    defs = NULL;
    ndef = 0;
}

/* Open the file to be processed. The handle variable holds internal
   info for the cpp emulator. Return integer status */
int cpp_open_file(const char *filenm, gmx_cpp_t *handle, char **cppopts)
{
    gmx_cpp_t    cpp;
    char        *buf, *pdum;
    char        *ptr, *ptr2;
    int          i;
    unsigned int i1;

    /* First process options, they might be necessary for opening files
       (especially include statements). */
    i  = 0;
    if (cppopts)
    {
        while (cppopts[i])
        {
            if (strstr(cppopts[i], "-I") == cppopts[i])
            {
                add_include(cppopts[i]+2);
            }
            if (strstr(cppopts[i], "-D") == cppopts[i])
            {
                /* If the option contains a =, split it into name and value. */
                ptr = strchr(cppopts[i], '=');
                if (ptr)
                {
                    buf = gmx_strndup(cppopts[i] + 2, ptr - cppopts[i] - 2);
                    add_define(buf, ptr + 1);
                    sfree(buf);
                }
                else
                {
                    add_define(cppopts[i] + 2, NULL);
                }
            }
            i++;
        }
    }
    if (debug)
    {
        fprintf(debug, "GMXCPP: added %d command line arguments\n", i);
    }

    snew(cpp, 1);
    *handle      = cpp;
    cpp->fn      = NULL;
    /* Find the file. First check whether it is in the current directory. */
    if (gmx_fexist(filenm))
    {
        cpp->fn = gmx_strdup(filenm);
    }
    else
    {
        /* If not, check all the paths given with -I. */
        for (i = 0; i < nincl; ++i)
        {
            snew(buf, strlen(incl[i]) + strlen(filenm) + 2);
            sprintf(buf, "%s/%s", incl[i], filenm);
            if (gmx_fexist(buf))
            {
                cpp->fn = buf;
                break;
            }
            sfree(buf);
        }
        /* If still not found, check the Gromacs library search path. */
        if (!cpp->fn)
        {
            cpp->fn = low_gmxlibfn(filenm, FALSE, FALSE);
        }
    }
    if (!cpp->fn)
    {
        gmx_fatal(FARGS, "Topology include file \"%s\" not found", filenm);
    }
    if (NULL != debug)
    {
        fprintf(debug, "GMXCPP: cpp file open %s\n", cpp->fn);
    }
    /* If the file name has a path component, we need to change to that
     * directory. Note that we - just as C - always use UNIX path separators
     * internally in include file names.
     */
    ptr  = strrchr(cpp->fn, '/');
    ptr2 = strrchr(cpp->fn, DIR_SEPARATOR);

    if (ptr == NULL || (ptr2 != NULL && ptr2 > ptr))
    {
        ptr = ptr2;
    }
    if (ptr == NULL)
    {
        cpp->path = NULL;
        cpp->cwd  = NULL;
    }
    else
    {
        cpp->path = cpp->fn;
        *ptr      = '\0';
        cpp->fn   = gmx_strdup(ptr+1);
        snew(cpp->cwd, STRLEN);

        gmx_getcwd(cpp->cwd, STRLEN);
        if (NULL != debug)
        {
            fprintf(debug, "GMXCPP: cwd %s\n", cpp->cwd);
        }
        gmx_chdir(cpp->path);

        if (NULL != debug)
        {
            fprintf(debug, "GMXCPP: chdir to %s\n", cpp->path);
        }
    }
    cpp->line_len = 0;
    cpp->line     = NULL;
    cpp->line_nr  = 0;
    cpp->nifdef   = 0;
    cpp->ifdefs   = NULL;
    cpp->child    = NULL;
    cpp->parent   = NULL;
    if (cpp->fp == NULL)
    {
        if (NULL != debug)
        {
            fprintf(debug, "GMXCPP: opening file %s\n", cpp->fn);
        }
        cpp->fp = fopen(cpp->fn, "r");
    }
    if (cpp->fp == NULL)
    {
        switch (errno)
        {
            case EINVAL:
            default:
                return eCPP_UNKNOWN;
        }
    }
    return eCPP_OK;
}

static int
process_directive(gmx_cpp_t *handlep, const char *dname, const char *dval)
{
    gmx_cpp_t    handle = (gmx_cpp_t)*handlep;
    int          i, i0, len, status;
    unsigned int i1;
    char        *inc_fn, *name;
    const char  *ptr;
    int          bIfdef, bIfndef;

    /* #ifdef or ifndef statement */
    bIfdef  = (strcmp(dname, "ifdef") == 0);
    bIfndef = (strcmp(dname, "ifndef") == 0);
    if (bIfdef || bIfndef)
    {
        if ((handle->nifdef > 0) && (handle->ifdefs[handle->nifdef-1] != eifTRUE))
        {
            handle->nifdef++;
            srenew(handle->ifdefs, handle->nifdef);
            handle->ifdefs[handle->nifdef-1] = eifIGNORE;
        }
        else
        {
            snew(name, strlen(dval)+1);
            sscanf(dval, "%s", name);
            for (i = 0; (i < ndef); i++)
            {
                if (strcmp(defs[i].name, name) == 0)
                {
                    break;
                }
            }
            handle->nifdef++;
            srenew(handle->ifdefs, handle->nifdef);
            if ((bIfdef && (i < ndef)) || (bIfndef && (i == ndef)))
            {
                handle->ifdefs[handle->nifdef-1] = eifTRUE;
            }
            else
            {
                handle->ifdefs[handle->nifdef-1] = eifFALSE;
            }
            sfree(name);
        }
        return eCPP_OK;
    }

    /* #else statement */
    if (strcmp(dname, "else") == 0)
    {
        if (handle->nifdef <= 0)
        {
            return eCPP_SYNTAX;
        }
        if (handle->ifdefs[handle->nifdef-1] == eifTRUE)
        {
            handle->ifdefs[handle->nifdef-1] = eifFALSE;
        }
        else if (handle->ifdefs[handle->nifdef-1] == eifFALSE)
        {
            handle->ifdefs[handle->nifdef-1] = eifTRUE;
        }
        return eCPP_OK;
    }

    /* #endif statement */
    if (strcmp(dname, "endif") == 0)
    {
        if (handle->nifdef <= 0)
        {
            return eCPP_SYNTAX;
        }
        handle->nifdef--;
        return eCPP_OK;
    }

    /* Check whether we're not ifdeffed out. The order of this statement
       is important. It has to come after #ifdef, #else and #endif, but
       anything else should be ignored. */
    if (is_ifdeffed_out(handle))
    {
        return eCPP_OK;
    }

    /* Check for include statements */
    if (strcmp(dname, "include") == 0)
    {
        len = -1;
        i0  = 0;
        for (i1 = 0; (i1 < strlen(dval)); i1++)
        {
            if ((dval[i1] == '"') || (dval[i1] == '<') || (dval[i1] == '>'))
            {
                if (len == -1)
                {
                    i0  = i1+1;
                    len = 0;
                }
                else
                {
                    break;
                }
            }
            else if (len >= 0)
            {
                len++;
            }
        }
        if (len == -1)
        {
            return eCPP_SYNTAX;
        }
        snew(inc_fn, len+1);
        strncpy(inc_fn, dval+i0, len);
        inc_fn[len] = '\0';

        if (debug)
        {
            fprintf(debug, "Going to open include file '%s' i0 = %d, strlen = %d\n",
                    inc_fn, i0, len);
        }
        /* Open include file and store it as a child in the handle structure */
        status = cpp_open_file(inc_fn, &(handle->child), NULL);
        sfree(inc_fn);
        if (status != eCPP_OK)
        {
            handle->child = NULL;
            return status;
        }
        /* Make a linked list of open files and move on to the include file */
        handle->child->parent = handle;
        *handlep              = handle->child;
        handle                = *handlep;
        return eCPP_OK;
    }

    /* #define statement */
    if (strcmp(dname, "define") == 0)
    {
        /* Split it into name and value. */
        ptr = dval;
        while ((*ptr != '\0') && !isspace(*ptr))
        {
            ptr++;
        }
        name = gmx_strndup(dval, ptr - dval);

        while ((*ptr != '\0') && isspace(*ptr))
        {
            ptr++;
        }

        add_define(name, ptr);
        sfree(name);
        return eCPP_OK;
    }

    /* #undef statement */
    if (strcmp(dname, "undef") == 0)
    {
        snew(name, strlen(dval)+1);
        sscanf(dval, "%s", name);
        for (i = 0; (i < ndef); i++)
        {
            if (strcmp(defs[i].name, name) == 0)
            {
                sfree(defs[i].name);
                sfree(defs[i].def);
                break;
            }
        }
        sfree(name);
        for (; (i < ndef-1); i++)
        {
            defs[i].name = defs[i+1].name;
            defs[i].def  = defs[i+1].def;
        }
        ndef--;

        return eCPP_OK;
    }

    /* If we haven't matched anything, this is an unknown directive */
    return eCPP_SYNTAX;
}

/* Return one whole line from the file into buf which holds at most n
   characters, for subsequent processing. Returns integer status. This
   routine also does all the "intelligent" work like processing cpp
   directives and so on. Note that often the routine is called
   recursively and no cpp directives are printed. */
int cpp_read_line(gmx_cpp_t *handlep, int n, char buf[])
{
    gmx_cpp_t   handle = (gmx_cpp_t)*handlep;
    int         i, nn, len, status;
    const char *ptr, *ptr2;
    char       *name;
    char       *dname, *dval;
    gmx_bool    bEOF;

    if (!handle)
    {
        return eCPP_INVALID_HANDLE;
    }
    if (!handle->fp)
    {
        return eCPP_FILE_NOT_OPEN;
    }

    bEOF = feof(handle->fp);
    if (!bEOF)
    {
        /* Read the actual line now. */
        if (fgets2(buf, n-1, handle->fp) == NULL)
        {
            /* Recheck EOF, since we could have been at the end before
             * the fgets2 call, but we need to read past the end to know.
             */
            bEOF = feof(handle->fp);
            if (!bEOF)
            {
                /* Something strange happened, fgets returned NULL,
                 * but we are not at EOF.
                 */
                return eCPP_UNKNOWN;
            }
        }
    }

    if (bEOF)
    {
        if (handle->parent == NULL)
        {
            return eCPP_EOF;
        }
        cpp_close_file(handlep);
        *handlep      = handle->parent;
        handle->child = NULL;
        return cpp_read_line(handlep, n, buf);
    }
    else
    {
        if (n > handle->line_len)
        {
            handle->line_len = n;
            srenew(handle->line, n);
        }
        strcpy(handle->line, buf);
        handle->line_nr++;
    }
    /* Now we've read a line! */
    if (debug)
    {
        fprintf(debug, "%s : %4d : %s\n", handle->fn, handle->line_nr, buf);
    }

    /* Process directives if this line contains one */
    if (find_directive(buf, &dname, &dval))
    {
        status = process_directive(handlep, dname, dval);
        if (status != eCPP_OK)
        {
            return status;
        }
        /* Don't print lines with directives, go on to the next */
        return cpp_read_line(handlep, n, buf);
    }

    /* Check whether we're not ifdeffed out. The order of this statement
       is important. It has to come after #ifdef, #else and #endif, but
       anything else should be ignored. */
    if (is_ifdeffed_out(handle))
    {
        return cpp_read_line(handlep, n, buf);
    }

    /* Check whether we have any defines that need to be replaced. Note
       that we have to use a best fit algorithm, rather than first come
       first go. We do this by sorting the defines on length first, and
       then on alphabetical order. */
    for (i = 0; (i < ndef); i++)
    {
        if (defs[i].def)
        {
            nn  = 0;
            ptr = buf;
            while ((ptr = strstrw(ptr, defs[i].name)) != NULL)
            {
                nn++;
                ptr += strlen(defs[i].name);
            }
            if (nn > 0)
            {
                len = strlen(buf) + nn*max(4, 4+strlen(defs[i].def)-strlen(defs[i].name));
                snew(name, len);
                ptr = buf;
                while ((ptr2 = strstrw(ptr, defs[i].name)) != NULL)
                {
                    strncat(name, ptr, (int)(ptr2-ptr));
                    strcat(name, defs[i].def);
                    ptr = ptr2 + strlen(defs[i].name);
                }
                strcat(name, ptr);
                strcpy(buf, name);
                sfree(name);
            }
        }
    }

    return eCPP_OK;
}

char *cpp_cur_file(const gmx_cpp_t *handlep)
{
    return (*handlep)->fn;
}

int cpp_cur_linenr(const gmx_cpp_t *handlep)
{
    return (*handlep)->line_nr;
}

/* Close the file! Return integer status. */
int cpp_close_file(gmx_cpp_t *handlep)
{
    int       i;
    gmx_cpp_t handle = (gmx_cpp_t)*handlep;

    if (!handle)
    {
        return eCPP_INVALID_HANDLE;
    }
    if (!handle->fp)
    {
        return eCPP_FILE_NOT_OPEN;
    }
    if (debug)
    {
        fprintf(debug, "GMXCPP: closing file %s\n", handle->fn);
    }
    fclose(handle->fp);
    if (NULL != handle->cwd)
    {
        if (NULL != debug)
        {
            fprintf(debug, "GMXCPP: chdir to %s\n", handle->cwd);
        }
        gmx_chdir(handle->cwd);
    }

    if (0)
    {
        switch (errno)
        {
            case 0:
                break;
            case ENOENT:
                return eCPP_FILE_NOT_FOUND;
            case EBADF:
                return eCPP_FILE_NOT_OPEN;
            case EINTR:
                return eCPP_INTERRUPT;
            default:
                if (debug)
                {
                    fprintf(debug, "Strange stuff closing file, errno = %d", errno);
                }
                return eCPP_UNKNOWN;
        }
    }
    handle->fp      = NULL;
    handle->line_nr = 0;
    if (NULL != handle->fn)
    {
        sfree(handle->fn);
        handle->fn = NULL;
    }
    if (NULL != handle->line)
    {
        sfree(handle->line);
        handle->line = NULL;
    }
    if (NULL != handle->ifdefs)
    {
        sfree(handle->ifdefs);
    }
    handle->nifdef = 0;
    if (NULL != handle->path)
    {
        sfree(handle->path);
    }
    if (NULL != handle->cwd)
    {
        sfree(handle->cwd);
    }

    return eCPP_OK;
}

void cpp_done()
{
    done_includes();
    done_defines();
}

/* Return a string containing the error message coresponding to status
   variable */
char *cpp_error(gmx_cpp_t *handlep, int status)
{
    char        buf[256];
    const char *ecpp[] = {
        "OK", "File not found", "End of file", "Syntax error", "Interrupted",
        "Invalid file handle",
        "File not open", "Unknown error", "Error status out of range"
    };
    gmx_cpp_t   handle = (gmx_cpp_t)*handlep;

    if (!handle)
    {
        return (char *)ecpp[eCPP_INVALID_HANDLE];
    }

    if ((status < 0) || (status >= eCPP_NR))
    {
        status = eCPP_NR;
    }

    sprintf(buf, "%s - File %s, line %d\nLast line read:\n'%s'",
            ecpp[status],
            (handle && handle->fn) ? handle->fn : "unknown",
            (handle) ? handle->line_nr : -1,
            handle->line ? handle->line : "");

    return gmx_strdup(buf);
}
