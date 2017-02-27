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
#include "gmxpre.h"

#include "filenm.h"

#include <cstdio>
#include <cstring>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

/* Use bitflag ... */
#define IS_SET(fn) ((fn.flag &ffSET) != 0)
#define IS_OPT(fn) ((fn.flag &ffOPT) != 0)

const t_filenm *getFilenm(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            return &fnm[i];
        }
    }

    return nullptr;
}

const char *opt2fn(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (std::strcmp(opt, fnm[i].opt) == 0)
        {
            return fnm[i].fns[0];
        }
    }

    fprintf(stderr, "No option %s\n", opt);

    return nullptr;
}

int opt2fns(char **fns[], const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            *fns = fnm[i].fns;
            return fnm[i].nfiles;
        }
    }

    fprintf(stderr, "No option %s\n", opt);
    return 0;
}

const char *ftp2fn(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            return fnm[i].fns[0];
        }
    }

    fprintf(stderr, "ftp2fn: No filetype %s\n", ftp2ext_with_dot(ftp));
    return nullptr;
}

int ftp2fns(char **fns[], int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            *fns = fnm[i].fns;
            return fnm[i].nfiles;
        }
    }

    fprintf(stderr, "ftp2fn: No filetype %s\n", ftp2ext_with_dot(ftp));
    return 0;
}

gmx_bool ftp2bSet(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            return (gmx_bool) IS_SET(fnm[i]);
        }
    }

    fprintf(stderr, "ftp2fn: No filetype %s\n", ftp2ext_with_dot(ftp));

    return FALSE;
}

gmx_bool opt2bSet(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (std::strcmp(opt, fnm[i].opt) == 0)
        {
            return (gmx_bool) IS_SET(fnm[i]);
        }
    }

    fprintf(stderr, "No option %s\n", opt);

    return FALSE;
}

const char *opt2fn_null(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (std::strcmp(opt, fnm[i].opt) == 0)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
            {
                return nullptr;
            }
            else
            {
                return fnm[i].fns[0];
            }
        }
    }
    fprintf(stderr, "No option %s\n", opt);
    return nullptr;
}

const char *ftp2fn_null(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
            {
                return nullptr;
            }
            else
            {
                return fnm[i].fns[0];
            }
        }
    }
    fprintf(stderr, "ftp2fn: No filetype %s\n", ftp2ext_with_dot(ftp));
    return nullptr;
}

gmx_bool is_optional(const t_filenm *fnm)
{
    return ((fnm->flag & ffOPT) == ffOPT);
}

gmx_bool is_output(const t_filenm *fnm)
{
    return ((fnm->flag & ffWRITE) == ffWRITE);
}

gmx_bool is_set(const t_filenm *fnm)
{
    return ((fnm->flag & ffSET) == ffSET);
}

int add_suffix_to_output_names(t_filenm *fnm, int nfile, const char *suffix)
{
    int   i, j;
    char  buf[STRLEN], newname[STRLEN];
    char *extpos;

    for (i = 0; i < nfile; i++)
    {
        if (is_output(&fnm[i]) && fnm[i].ftp != efCPT)
        {
            /* We never use multiple _outputs_, but we might as well check
               for it, just in case... */
            for (j = 0; j < fnm[i].nfiles; j++)
            {
                std::strncpy(buf, fnm[i].fns[j], STRLEN - 1);
                extpos  = strrchr(buf, '.');
                *extpos = '\0';
                sprintf(newname, "%s%s.%s", buf, suffix, extpos + 1);
                sfree(fnm[i].fns[j]);
                fnm[i].fns[j] = gmx_strdup(newname);
            }
        }
    }
    return 0;
}

t_filenm *dup_tfn(int nf, const t_filenm tfn[])
{
    int       i, j;
    t_filenm *ret;

    snew(ret, nf);
    for (i = 0; i < nf; i++)
    {
        ret[i] = tfn[i]; /* just directly copy all non-string fields */
        if (tfn[i].opt)
        {
            ret[i].opt = gmx_strdup(tfn[i].opt);
        }
        else
        {
            ret[i].opt = nullptr;
        }

        if (tfn[i].fn)
        {
            ret[i].fn = gmx_strdup(tfn[i].fn);
        }
        else
        {
            ret[i].fn = nullptr;
        }

        if (tfn[i].nfiles > 0)
        {
            snew(ret[i].fns, tfn[i].nfiles);
            for (j = 0; j < tfn[i].nfiles; j++)
            {
                ret[i].fns[j] = gmx_strdup(tfn[i].fns[j]);
            }
        }
    }
    return ret;
}

void done_filenms(int nf, t_filenm fnm[])
{
    int i, j;

    for (i = 0; i < nf; ++i)
    {
        for (j = 0; j < fnm[i].nfiles; ++j)
        {
            sfree(fnm[i].fns[j]);
        }
        sfree(fnm[i].fns);
        fnm[i].fns = nullptr;
    }
}
