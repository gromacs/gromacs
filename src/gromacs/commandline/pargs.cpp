/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 *
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is completely threadsafe - keep it that way! */
#include "gromacs/commandline/pargs.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/string2.h"

#include "gromacs/utility/gmxassert.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

static void usage(const char *type, const char *arg)
{
    GMX_ASSERT(arg != NULL, "NULL command-line argument should not occur");
    gmx_fatal(FARGS, "Expected %s argument for option %s\n", type, arg);
}

/* Scan an int for argument argv[*i] from argument at argv[*i + 1].
 * eg: -p 32.  argv[*i] is only used for error reporting.
 * If there is no value, or the conversion is not successful, the
 * routine exits with an error, otherwise it returns the value found.
 * *i is incremented once.
 */
static int iscan(int argc, char *argv[], int *i)
{
    const char *const arg = argv[*i];
    if (argc <= (*i)+1)
    {
        usage("an integer", arg);
    }
    const char *const value = argv[++(*i)];
    char             *endptr;
    int               var = std::strtol(value, &endptr, 10);
    if (*value == '\0' || *endptr != '\0')
    {
        usage("an integer", arg);
    }
    return var;
}

/* Same as above, but for large integer values */
static gmx_large_int_t istepscan(int argc, char *argv[], int *i)
{
    const char *const arg = argv[*i];
    if (argc <= (*i)+1)
    {
        usage("an integer", arg);
    }
    const char *const value = argv[++(*i)];
    char             *endptr;
    gmx_large_int_t   var = str_to_large_int_t(value, &endptr);
    if (*value == '\0' || *endptr != '\0')
    {
        usage("an integer", arg);
    }
    return var;
}

/* Routine similar to the above, but working on doubles. */
static double dscan(int argc, char *argv[], int *i)
{
    const char *const arg = argv[*i];
    if (argc <= (*i)+1)
    {
        usage("a real", arg);
    }
    const char *const value = argv[++(*i)];
    char             *endptr;
    double            var = std::strtod(value, &endptr);
    if (*value == '\0' || *endptr != '\0')
    {
        usage("a real", arg);
    }
    return var;
}

/* Routine similar to the above, but working on strings. The pointer
 * returned is a pointer to the argv field.
 */
static char *sscan(int argc, char *argv[], int *i)
{
    if (argc > (*i)+1)
    {
        if ( (argv[(*i)+1][0] == '-') && (argc > (*i)+2) &&
             (argv[(*i)+2][0] != '-') )
        {
            fprintf(stderr, "Possible missing string argument for option %s\n\n",
                    argv[*i]);
        }
    }
    else
    {
        usage("a string", argv[*i]);
    }

    return argv[++(*i)];
}

gmx_bool is_hidden(t_pargs *pa)
{
    return ((strstr(pa->desc, "HIDDEN") != NULL) ||
            (strstr(pa->desc, "[hidden]") != NULL));
}

int nenum(const char *const enumc[])
{
    int i;

    i = 1;
    /* we *can* compare pointers directly here! */
    while (enumc[i] && enumc[0] != enumc[i])
    {
        i++;
    }

    return i;
}

void get_pargs(int *argc, char *argv[], int nparg, t_pargs pa[], gmx_bool bKeepArgs)
{
    int       i, j, k, match;
    gmx_bool *bKeep;
    char      buf[32];
    char     *ptr;

    snew(bKeep, *argc+1);
    bKeep[0]     = TRUE;
    bKeep[*argc] = TRUE;

    for (i = 1; (i < *argc); i++)
    {
        bKeep[i] = TRUE;
        for (j = 0; (j < nparg); j++)
        {
            if (pa[j].type == etBOOL)
            {
                sprintf(buf, "-no%s", pa[j].option+1);
                if (strcmp(pa[j].option, argv[i]) == 0)
                {
                    *pa[j].u.b = TRUE;
                    pa[j].bSet = TRUE;
                    bKeep[i]   = FALSE;
                }
                else if (strcmp(buf, argv[i]) == 0)
                {
                    *pa[j].u.b = FALSE;
                    pa[j].bSet = TRUE;
                    bKeep[i]   = FALSE;
                }
            }
            else if (strcmp(pa[j].option, argv[i]) == 0)
            {
                if (pa[j].bSet)
                {
                    fprintf(stderr, "Setting option %s more than once!\n",
                            pa[j].option);
                }
                pa[j].bSet = TRUE;
                bKeep[i]   = FALSE;
                switch (pa[j].type)
                {
                    case etINT:
                        *pa[j].u.i = iscan(*argc, argv, &i);
                        break;
                    case etGMX_LARGE_INT:
                        *pa[j].u.is = istepscan(*argc, argv, &i);
                        break;
                    case etTIME:
                    case etREAL:
                        *pa[j].u.r = dscan(*argc, argv, &i);
                        break;
                    case etSTR:
                        *(pa[j].u.c) = sscan(*argc, argv, &i);
                        break;
                    case etENUM:
                        match = -1;
                        ptr   = sscan(*argc, argv, &i);
                        for (k = 1; (pa[j].u.c[k] != NULL); k++)
                        {
                            /* only check ptr against beginning of
                               pa[j].u.c[k] */
                            if (gmx_strncasecmp(ptr, pa[j].u.c[k], strlen(ptr)) == 0)
                            {
                                if ( ( match == -1 ) ||
                                     ( strlen(pa[j].u.c[k]) <
                                       strlen(pa[j].u.c[match]) ) )
                                {
                                    match = k;
                                }
                            }
                        }
                        if (match != -1)
                        {
                            pa[j].u.c[0] = pa[j].u.c[match];
                        }
                        else
                        {
                            gmx_fatal(FARGS, "Invalid argument %s for option %s",
                                      ptr, pa[j].option);
                        }
                        break;
                    case etRVEC:
                        (*pa[j].u.rv)[0] = dscan(*argc, argv, &i);
                        if ( (i+1 == *argc) ||
                             ( (argv[i+1][0] == '-') &&
                               !isdigit(argv[i+1][1]) ) )
                        {
                            (*pa[j].u.rv)[1]     =
                                (*pa[j].u.rv)[2] =
                                    (*pa[j].u.rv)[0];
                        }
                        else
                        {
                            bKeep[i]         = FALSE;
                            (*pa[j].u.rv)[1] = dscan(*argc, argv, &i);
                            if ( (i+1 == *argc) ||
                                 ( (argv[i+1][0] == '-') &&
                                   !isdigit(argv[i+1][1]) ) )
                            {
                                gmx_fatal(FARGS,
                                          "%s: vector must have 1 or 3 real parameters",
                                          pa[j].option);
                            }
                            bKeep[i]         = FALSE;
                            (*pa[j].u.rv)[2] = dscan(*argc, argv, &i);
                        }
                        break;
                    default:
                        gmx_fatal(FARGS, "Invalid type %d in pargs", pa[j].type);
                }
                /* i may be incremented, so set it to not keep */
                bKeep[i] = FALSE;
            }
        }
    }
    if (!bKeepArgs)
    {
        /* Remove used entries */
        for (i = j = 0; (i <= *argc); i++)
        {
            if (bKeep[i])
            {
                argv[j++] = argv[i];
            }
        }
        (*argc) = j-1;
    }
    sfree(bKeep);
}

int opt2parg_int(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *pa[i].u.i;
        }
    }

    gmx_fatal(FARGS, "No integer option %s in pargs", option);

    return 0;
}

gmx_bool opt2parg_gmx_bool(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *pa[i].u.b;
        }
    }

    gmx_fatal(FARGS, "No boolean option %s in pargs", option);

    return FALSE;
}

real opt2parg_real(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *pa[i].u.r;
        }
    }

    gmx_fatal(FARGS, "No real option %s in pargs", option);

    return 0.0;
}

const char *opt2parg_str(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return *(pa[i].u.c);
        }
    }

    gmx_fatal(FARGS, "No string option %s in pargs", option);

    return NULL;
}

gmx_bool opt2parg_bSet(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return pa[i].bSet;
        }
    }

    gmx_fatal(FARGS, "No such option %s in pargs", option);

    return FALSE; /* Too make some compilers happy */
}

const char *opt2parg_enum(const char *option, int nparg, t_pargs pa[])
{
    int i;

    for (i = 0; (i < nparg); i++)
    {
        if (strcmp(pa[i].option, option) == 0)
        {
            return pa[i].u.c[0];
        }
    }

    gmx_fatal(FARGS, "No such option %s in pargs", option);

    return NULL;
}
