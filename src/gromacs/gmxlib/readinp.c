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

#include "gromacs/legacyheaders/readinp.h"

#include <stdio.h>
#include <stdlib.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/warninp.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/qsort_threadsafe.h"
#include "gromacs/utility/smalloc.h"

t_inpfile *read_inpfile(const char *fn, int *ninp,
                        warninp_t wi)
{
    FILE      *in;
    char       buf[STRLEN], lbuf[STRLEN], rbuf[STRLEN], warn_buf[STRLEN];
    char      *ptr, *cptr;
    t_inpfile *inp = NULL;
    int        nin, lc, i, j, k;
    /* setting cppopts from command-line options would be cooler */
    gmx_bool   allow_override = FALSE;


    if (debug)
    {
        fprintf(debug, "Reading MDP file %s\n", fn);
    }

    in = gmx_ffopen(fn, "r");

    nin = lc  = 0;
    do
    {
        ptr = fgets2(buf, STRLEN-1, in);
        lc++;
        set_warning_line(wi, fn, lc);
        if (ptr)
        {
            /* Strip comment */
            if ((cptr = strchr(buf, COMMENTSIGN)) != NULL)
            {
                *cptr = '\0';
            }
            /* Strip spaces */
            trim(buf);

            for (j = 0; (buf[j] != '=') && (buf[j] != '\0'); j++)
            {
                ;
            }
            if (buf[j] == '\0')
            {
                if (j > 0)
                {
                    if (debug)
                    {
                        fprintf(debug, "No = on line %d in file %s, ignored\n", lc, fn);
                    }
                }
            }
            else
            {
                for (i = 0; (i < j); i++)
                {
                    lbuf[i] = buf[i];
                }
                lbuf[i] = '\0';
                trim(lbuf);
                if (lbuf[0] == '\0')
                {
                    if (debug)
                    {
                        fprintf(debug, "Empty left hand side on line %d in file %s, ignored\n", lc, fn);
                    }
                }
                else
                {
                    for (i = j+1, k = 0; (buf[i] != '\0'); i++, k++)
                    {
                        rbuf[k] = buf[i];
                    }
                    rbuf[k] = '\0';
                    trim(rbuf);
                    if (rbuf[0] == '\0')
                    {
                        if (debug)
                        {
                            fprintf(debug, "Empty right hand side on line %d in file %s, ignored\n", lc, fn);
                        }
                    }
                    else
                    {
                        /* Now finally something sensible */
                        int found_index;

                        /* first check whether we hit the 'multiple_entries' option */
                        if (gmx_strcasecmp_min(eMultentOpt_names[eMultentOptName], lbuf) == 0)
                        {
                            /* we now check whether to allow overrides from here or not */
                            if (gmx_strcasecmp_min(eMultentOpt_names[eMultentOptNo], rbuf) == 0)
                            {
                                allow_override = FALSE;
                            }
                            else if (gmx_strcasecmp_min(eMultentOpt_names[eMultentOptLast], rbuf) == 0)
                            {
                                allow_override = TRUE;
                            }
                            else
                            {
                                sprintf(warn_buf,
                                        "Parameter \"%s\" should either be %s or %s\n",
                                        lbuf,
                                        eMultentOpt_names[eMultentOptNo],
                                        eMultentOpt_names[eMultentOptLast]);
                                warning_error(wi, warn_buf);
                            }
                        }
                        else
                        {
                            /* it is a regular option; check for duplicates */
                            found_index = search_einp(nin, inp, lbuf);

                            if (found_index == -1)
                            {
                                /* add a new item */
                                srenew(inp, ++nin);
                                inp[nin-1].inp_count  = 1;
                                inp[nin-1].count      = 0;
                                inp[nin-1].bObsolete  = FALSE;
                                inp[nin-1].bSet       = FALSE;
                                inp[nin-1].name       = gmx_strdup(lbuf);
                                inp[nin-1].value      = gmx_strdup(rbuf);
                            }
                            else
                            {
                                if (!allow_override)
                                {
                                    sprintf(warn_buf,
                                            "Parameter \"%s\" doubly defined (and multiple assignments not allowed)\n",
                                            lbuf);
                                    warning_error(wi, warn_buf);
                                }
                                else
                                {
                                    /* override */
                                    sfree(inp[found_index].value);
                                    inp[found_index].value = gmx_strdup(rbuf);
                                    sprintf(warn_buf,
                                            "Overriding existing parameter \"%s\" with value \"%s\"\n",
                                            lbuf, rbuf);
                                    warning_note(wi, warn_buf);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    while (ptr);

    gmx_ffclose(in);

    if (debug)
    {
        fprintf(debug, "Done reading MDP file, there were %d entries in there\n",
                nin);
    }

    *ninp = nin;

    return inp;
}




static int inp_comp(const void *a, const void *b)
{
    return ((t_inpfile *)a)->count - ((t_inpfile *)b)->count;
}

static void sort_inp(int ninp, t_inpfile inp[])
{
    int i, mm;

    mm = -1;
    for (i = 0; (i < ninp); i++)
    {
        mm = max(mm, inp[i].count);
    }
    for (i = 0; (i < ninp); i++)
    {
        if (inp[i].count == 0)
        {
            inp[i].count = mm++;
        }
    }
    gmx_qsort(inp, ninp, (size_t)sizeof(inp[0]), inp_comp);
}

void write_inpfile(const char *fn, int ninp, t_inpfile inp[], gmx_bool bHaltOnUnknown,
                   warninp_t wi)
{
    FILE *out;
    int   i;
    char  warn_buf[STRLEN];

    sort_inp(ninp, inp);
    out = gmx_fio_fopen(fn, "w");
    nice_header(out, fn);
    for (i = 0; (i < ninp); i++)
    {
        if (inp[i].bSet)
        {
            if (inp[i].name[0] == ';' || (strlen(inp[i].name) > 2 && inp[i].name[1] == ';'))
            {
                fprintf(out, "%-24s\n", inp[i].name);
            }
            else
            {
                fprintf(out, "%-24s = %s\n", inp[i].name, inp[i].value ? inp[i].value : "");
            }
        }
        else if (!inp[i].bObsolete)
        {
            sprintf(warn_buf, "Unknown left-hand '%s' in parameter file\n",
                    inp[i].name);
            if (bHaltOnUnknown)
            {
                warning_error(wi, warn_buf);
            }
            else
            {
                warning(wi, warn_buf);
            }
        }
    }
    gmx_fio_fclose(out);

    check_warning_error(wi, FARGS);
}

void replace_inp_entry(int ninp, t_inpfile *inp, const char *old_entry, const char *new_entry)
{
    int  i;

    for (i = 0; (i < ninp); i++)
    {
        if (gmx_strcasecmp_min(old_entry, inp[i].name) == 0)
        {
            if (new_entry)
            {
                fprintf(stderr, "Replacing old mdp entry '%s' by '%s'\n",
                        inp[i].name, new_entry);
                sfree(inp[i].name);
                inp[i].name = gmx_strdup(new_entry);
            }
            else
            {
                fprintf(stderr, "Ignoring obsolete mdp entry '%s'\n",
                        inp[i].name);
                inp[i].bObsolete = TRUE;
            }
        }
    }
}

int search_einp(int ninp, const t_inpfile *inp, const char *name)
{
    int i;

    if (inp == NULL)
    {
        return -1;
    }
    for (i = 0; i < ninp; i++)
    {
        if (gmx_strcasecmp_min(name, inp[i].name) == 0)
        {
            return i;
        }
    }
    return -1;
}

static int get_einp(int *ninp, t_inpfile **inp, const char *name)
{
    int    i;
    int    notfound = FALSE;
    char   warn_buf[STRLEN];

/*  if (inp==NULL)
    return -1;
   for(i=0; (i<(*ninp)); i++)
    if (gmx_strcasecmp_min(name,(*inp)[i].name) == 0)
      break;
   if (i == (*ninp)) {*/
    i = search_einp(*ninp, *inp, name);
    if (i == -1)
    {
        notfound = TRUE;
        i        = (*ninp)++;
        srenew(*inp, (*ninp));
        (*inp)[i].name = gmx_strdup(name);
        (*inp)[i].bSet = TRUE;
    }
    (*inp)[i].count = (*inp)[0].inp_count++;
    (*inp)[i].bSet  = TRUE;
    if (debug)
    {
        fprintf(debug, "Inp %d = %s\n", (*inp)[i].count, (*inp)[i].name);
    }

    /*if (i == (*ninp)-1)*/
    if (notfound)
    {
        return -1;
    }
    else
    {
        return i;
    }
}

int get_eint(int *ninp, t_inpfile **inp, const char *name, int def,
             warninp_t wi)
{
    char buf[32], *ptr, warn_buf[STRLEN];
    int  ii;
    int  ret;

    ii = get_einp(ninp, inp, name);

    if (ii == -1)
    {
        sprintf(buf, "%d", def);
        (*inp)[(*ninp)-1].value = gmx_strdup(buf);

        return def;
    }
    else
    {
        ret = strtol((*inp)[ii].value, &ptr, 10);
        if (ptr == (*inp)[ii].value)
        {
            sprintf(warn_buf, "Right hand side '%s' for parameter '%s' in parameter file is not an integer value\n", (*inp)[ii].value, (*inp)[ii].name);
            warning_error(wi, warn_buf);
        }

        return ret;
    }
}

gmx_int64_t get_eint64(int *ninp, t_inpfile **inp,
                       const char *name, gmx_int64_t def,
                       warninp_t wi)
{
    char            buf[32], *ptr, warn_buf[STRLEN];
    int             ii;
    gmx_int64_t     ret;

    ii = get_einp(ninp, inp, name);

    if (ii == -1)
    {
        sprintf(buf, "%"GMX_PRId64, def);
        (*inp)[(*ninp)-1].value = gmx_strdup(buf);

        return def;
    }
    else
    {
        ret = str_to_int64_t((*inp)[ii].value, &ptr);
        if (ptr == (*inp)[ii].value)
        {
            sprintf(warn_buf, "Right hand side '%s' for parameter '%s' in parameter file is not an integer value\n", (*inp)[ii].value, (*inp)[ii].name);
            warning_error(wi, warn_buf);
        }

        return ret;
    }
}

double get_ereal(int *ninp, t_inpfile **inp, const char *name, double def,
                 warninp_t wi)
{
    char   buf[32], *ptr, warn_buf[STRLEN];
    int    ii;
    double ret;

    ii = get_einp(ninp, inp, name);

    if (ii == -1)
    {
        sprintf(buf, "%g", def);
        (*inp)[(*ninp)-1].value = gmx_strdup(buf);

        return def;
    }
    else
    {
        ret = strtod((*inp)[ii].value, &ptr);
        if (ptr == (*inp)[ii].value)
        {
            sprintf(warn_buf, "Right hand side '%s' for parameter '%s' in parameter file is not a real value\n", (*inp)[ii].value, (*inp)[ii].name);
            warning_error(wi, warn_buf);
        }

        return ret;
    }
}

const char *get_estr(int *ninp, t_inpfile **inp, const char *name, const char *def)
{
    char buf[32];
    int  ii;

    ii = get_einp(ninp, inp, name);

    if (ii == -1)
    {
        if (def)
        {
            sprintf(buf, "%s", def);
            (*inp)[(*ninp)-1].value = gmx_strdup(buf);
        }
        else
        {
            (*inp)[(*ninp)-1].value = NULL;
        }

        return def;
    }
    else
    {
        return (*inp)[ii].value;
    }
}

int get_eeenum(int *ninp, t_inpfile **inp, const char *name, const char **defs,
               warninp_t wi)
{
    int  ii, i, j;
    int  n = 0;
    char buf[STRLEN];

    ii = get_einp(ninp, inp, name);

    if (ii == -1)
    {
        (*inp)[(*ninp)-1].value = gmx_strdup(defs[0]);

        return 0;
    }

    for (i = 0; (defs[i] != NULL); i++)
    {
        if (gmx_strcasecmp_min(defs[i], (*inp)[ii].value) == 0)
        {
            break;
        }
    }

    if (defs[i] == NULL)
    {
        n += sprintf(buf, "Invalid enum '%s' for variable %s, using '%s'\n",
                     (*inp)[ii].value, name, defs[0]);
        n += sprintf(buf+n, "Next time use one of:");
        j  = 0;
        while (defs[j])
        {
            n += sprintf(buf+n, " '%s'", defs[j]);
            j++;
        }
        if (wi != NULL)
        {
            warning_error(wi, buf);
        }
        else
        {
            fprintf(stderr, "%s\n", buf);
        }

        (*inp)[ii].value = gmx_strdup(defs[0]);

        return 0;
    }

    return i;
}

int get_eenum(int *ninp, t_inpfile **inp, const char *name, const char **defs)
{
    int dum = 0;

    return get_eeenum(ninp, inp, name, defs, NULL);
}
