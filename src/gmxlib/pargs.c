/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "typedefs.h"
#include "gmx_fatal.h"
#include "statutil.h"
#include "readinp.h"
#include "smalloc.h"
#include "names.h"
#include "string2.h"
#include "vec.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

gmx_bool is_hidden(t_pargs *pa)
{
    return ((strstr(pa->desc, "HIDDEN") != NULL) ||
            (strstr(pa->desc, "[hidden]") != NULL));
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
                        match = NOTSET;
                        ptr   = sscan(*argc, argv, &i);
                        for (k = 1; (pa[j].u.c[k] != NULL); k++)
                        {
                            /* only check ptr against beginning of
                               pa[j].u.c[k] */
                            if (gmx_strncasecmp(ptr, pa[j].u.c[k], strlen(ptr)) == 0)
                            {
                                if ( ( match == NOTSET ) ||
                                     ( strlen(pa[j].u.c[k]) <
                                       strlen(pa[j].u.c[match]) ) )
                                {
                                    match = k;
                                }
                            }
                        }
                        if (match != NOTSET)
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

    gmx_fatal(FARGS, "No gmx_boolean option %s in pargs", option);

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

char *pa_val(t_pargs *pa, char buf[], int sz)
{
    real r;
    char buf_str[1256]; buf_str[0] = '\0';

    buf[0] = '\0';

    if (sz < 255)
    {
        gmx_fatal(FARGS, "Buffer must be at least 255 chars\n");
    }

    switch (pa->type)
    {
        case etINT:
            sprintf(buf, "%-d", *(pa->u.i));
            break;
        case etGMX_LARGE_INT:
            sprintf(buf, gmx_large_int_pfmt, *(pa->u.is));
            break;
        case etTIME:
        case etREAL:
            r = *(pa->u.r);
            sprintf(buf_str, "%-6g", r);
            strcpy(buf, buf_str);
            break;
        case etBOOL:
            sprintf(buf, "%-6s", *(pa->u.b) ? "yes" : "no");
            break;
        case etSTR:
            if (*(pa->u.c))
            {
                if (strlen(*(pa->u.c)) >= (size_t)sz)
                {
                    gmx_fatal(FARGS, "Argument too long: \"%d\"\n", *(pa->u.c));
                }
                else
                {
                    strcpy(buf, *(pa->u.c));
                }
            }
            break;
        case etENUM:
            strcpy(buf, *(pa->u.c));
            break;
        case etRVEC:
            sprintf(buf, "%g %g %g", (*pa->u.rv)[0],
                    (*pa->u.rv)[1],
                    (*pa->u.rv)[2]);
            break;
    }
    return buf;
}

#define OPTLEN 12
#define TYPELEN 6
#define LONGSTR 1024
char *pargs_print_line(t_pargs *pa, gmx_bool bLeadingSpace)
{
    char buf[LONGSTR], *buf2, *tmp, *desc;

    snew(buf2, LONGSTR+strlen(pa->desc));
    snew(tmp, LONGSTR+strlen(pa->desc));

    if (pa->type == etBOOL)
    {
        sprintf(buf, "-[no]%s", pa->option+1);
    }
    else
    {
        strcpy(buf, pa->option);
    }
    desc = check_tty(pa->desc);
    if (strlen(buf) > ((OPTLEN+TYPELEN)-max(strlen(argtp[pa->type]), 4)))
    {
        sprintf(buf2, "%s%s %-6s %-6s  %-s\n",
                bLeadingSpace ? " " : "", buf,
                argtp[pa->type], pa_val(pa, tmp, LONGSTR-1), desc);
    }
    else if (strlen(buf) > OPTLEN)
    {
        /* so type can be 3 or 4 char's, this fits in the %4s */
        sprintf(buf2, "%s%-14s %-4s %-6s  %-s\n",
                bLeadingSpace ? " " : "", buf, argtp[pa->type],
                pa_val(pa, tmp, LONGSTR-1), desc);
    }
    else
    {
        sprintf(buf2, "%s%-12s %-6s %-6s  %-s\n",
                bLeadingSpace ? " " : "", buf, argtp[pa->type],
                pa_val(pa, tmp, LONGSTR-1), desc);
    }
    sfree(desc);
    sfree(tmp);

    tmp = wrap_lines(buf2, 78, 28, FALSE);

    sfree(buf2);

    return tmp;
}

void print_pargs(FILE *fp, int npargs, t_pargs pa[], gmx_bool bLeadingSpace)
{
    gmx_bool bShowHidden;
    char     buf[32], buf2[256], tmp[256];
    char    *wdesc;
    int      i;

    /* Cannot call opt2parg_bSet here, because it crashes when the option
     * is not in the list (mdrun)
     */
    bShowHidden = FALSE;
    for (i = 0; (i < npargs); i++)
    {
        if ((strcmp(pa[i].option, "-hidden") == 0) && (pa[i].bSet))
        {
            bShowHidden = TRUE;
        }
    }

    if (npargs > 0)
    {
        fprintf(fp, "%s%-12s %-6s %-6s  %-s\n",
                bLeadingSpace ? " " : "", "Option", "Type", "Value", "Description");
        fprintf(fp, "%s------------------------------------------------------\n",
                bLeadingSpace ? " " : "");
        for (i = 0; (i < npargs); i++)
        {
            if (bShowHidden || !is_hidden(&pa[i]))
            {
                wdesc = pargs_print_line(&pa[i], bLeadingSpace);
                fprintf(fp, "%s", wdesc);
                sfree(wdesc);
            }
        }
        fprintf(fp, "\n");
    }
}

void pr_enums(FILE *fp, int npargs, t_pargs pa[], int shell)
{
    int i, j;

    switch (shell)
    {
        case eshellCSH:
            for (i = 0; i < npargs; i++)
            {
                if (pa[i].type == etENUM)
                {
                    fprintf(fp, " \"n/%s/(", pa[i].option);
                    for (j = 1; pa[i].u.c[j]; j++)
                    {
                        fprintf(fp, " %s", pa[i].u.c[j]);
                    }
                    fprintf(fp, ")/\"");
                }
            }
            break;
        case eshellBASH:
            for (i = 0; i < npargs; i++)
            {
                if (pa[i].type == etENUM)
                {
                    fprintf(fp, "%s) COMPREPLY=( $(compgen -W '", pa[i].option);
                    for (j = 1; pa[i].u.c[j]; j++)
                    {
                        fprintf(fp, " %s", pa[i].u.c[j]);
                    }
                    fprintf(fp, " ' -- $c ));;\n");
                }
            }
            break;
        case eshellZSH:
            for (i = 0; i < npargs; i++)
            {
                if (pa[i].type == etENUM)
                {
                    fprintf(fp, "- 'c[-1,%s]' -s \"", pa[i].option);
                    for (j = 1; pa[i].u.c[j]; j++)
                    {
                        fprintf(fp, " %s", pa[i].u.c[j]);
                    }
                    fprintf(fp, "\" ");
                }
            }
            break;
    }
}
