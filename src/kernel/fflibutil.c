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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "sysstuff.h"
#include "string2.h"
#include "futil.h"
#include "network.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "statutil.h"

#ifdef GMX_NATIVE_WINDOWS
#include <direct.h>
#include <io.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef GMX_THREAD_MPI
#include "thread_mpi.h"
#endif

#include "fflibutil.h"

const char *fflib_forcefield_dir_ext()
{
    return ".ff";
}

const char *fflib_forcefield_itp()
{
    return "forcefield.itp";
}

const char *fflib_forcefield_doc()
{
    return "forcefield.doc";
}

void fflib_filename_base(const char *filename, char *filebase, int maxlen)
{
    const char *cptr;
    char       *ptr;

    cptr = strrchr(filename, DIR_SEPARATOR);
    if (cptr != NULL)
    {
        /* Skip the separator */
        cptr += 1;
    }
    else
    {
        cptr = filename;
    }
    if (strlen(filename) >= (size_t)maxlen)
    {
        gmx_fatal(FARGS, "filename is longer (%d) than maxlen (%d)",
                  strlen(filename), maxlen);
    }
    strcpy(filebase, cptr);
    /* Remove the extension */
    ptr = strrchr(filebase, '.');
    if (ptr != NULL)
    {
        ptr[0] = '\0';
    }
}

static void sort_filenames(int n, char **name, char **name2)
{
    /* Slow sort, but we usually have tens of names */
    int   i, j, f;
    char *tmp;

    for (i = 0; i < n-1; i++)
    {
        f = i;
        for (j = i+1; j < n; j++)
        {
            if (strcmp(name[j], name[f]) < 0)
            {
                f = j;
            }
        }
        if (f > i)
        {
            tmp     = name[i];
            name[i] = name[f];
            name[f] = tmp;
            if (name2 != NULL)
            {
                tmp      = name2[i];
                name2[i] = name2[f];
                name2[f] = tmp;
            }
        }
    }
}

static int low_fflib_search_file_end(const char *ffdir,
                                     gmx_bool    bAddCWD,
                                     const char *file_end,
                                     gmx_bool    bFatalError,
                                     char     ***filenames,
                                     char     ***filenames_short)
{
    char           *ret = NULL;
    char           *lib, *dir;
    char            buf[1024];
    char           *libpath;
    gmx_bool        env_is_set;
    int             len_fe, len_name;
    char          **fns, **fns_short;
    char            dir_print[GMX_PATH_MAX];
    char           *pdum;
    char           *s, fn_dir[GMX_PATH_MAX];
    gmx_directory_t dirhandle;
    char            nextname[STRLEN];
    int             n, n_thisdir, rc;

    len_fe = strlen(file_end);

    env_is_set = FALSE;
    if (ffdir != NULL)
    {
        /* Search in current dir and ffdir */
        libpath = gmxlibfn(ffdir);
    }
    else
    {
        /* GMXLIB can be a path now */
        lib = getenv("GMXLIB");
        snew(libpath, GMX_PATH_MAX);
        if (bAddCWD)
        {
            sprintf(libpath, "%s%s", ".", PATH_SEPARATOR);
        }
        if (lib != NULL)
        {
            env_is_set = TRUE;
            strncat(libpath, lib, GMX_PATH_MAX);
        }
        else if (!get_libdir(libpath+strlen(libpath)))
        {
            strncat(libpath, GMXLIBDIR, GMX_PATH_MAX);
        }
    }
    s         = libpath;
    n         = 0;
    fns       = NULL;
    fns_short = NULL;
    /* Loop over all the entries in libpath */
    while ((dir = gmx_strsep(&s, PATH_SEPARATOR)) != NULL)
    {
        rc = gmx_directory_open(&dirhandle, dir);
        if (rc == 0)
        {
            strcpy(dir_print, dir);

            n_thisdir = 0;
            while (gmx_directory_nextfile(dirhandle, nextname, STRLEN-1) == 0)
            {
                nextname[STRLEN-1] = 0;
                if (debug)
                {
                    fprintf(debug, "dir '%s' %d file '%s'\n",
                            dir, n_thisdir, nextname);
                }
                len_name = strlen(nextname);
                /* What about case sensitivity? */
                if (len_name >= len_fe &&
                    strcmp(nextname+len_name-len_fe, file_end) == 0)
                {
                    /* We have a match */
                    srenew(fns, n+1);
                    sprintf(fn_dir, "%s%c%s",
                            dir_print, DIR_SEPARATOR, nextname);

                    /* Copy the file name, possibly including the path. */
                    fns[n] = strdup(fn_dir);

                    if (ffdir == NULL)
                    {
                        /* We are searching in a path.
                         * Use the relative path when we use share/top
                         * from the installation.
                         * Add the full path when we use the current
                         * working directory of GMXLIB.
                         */
                        srenew(fns_short, n+1);
                        if (strcmp(dir, ".") == 0 || env_is_set)
                        {
                            fns_short[n] = strdup(fn_dir);
                        }
                        else
                        {
                            fns_short[n] = strdup(nextname);
                        }
                    }
                    n++;
                    n_thisdir++;
                }
            }
            gmx_directory_close(dirhandle);

            sort_filenames(n_thisdir,
                           fns+n-n_thisdir,
                           fns_short == NULL ? NULL : fns_short+n-n_thisdir);
        }
    }

    sfree(libpath);

    if (n == 0 && bFatalError)
    {
        if (ffdir != NULL)
        {
            gmx_fatal(FARGS, "Could not find any files ending on '%s' in the force field directory '%s'", file_end, ffdir);
        }
        else
        {
            gmx_fatal(FARGS, "Could not find any files ending on '%s' in the current directory or the GROMACS library search path", file_end);
        }
    }

    *filenames = fns;
    if (ffdir == NULL)
    {
        *filenames_short = fns_short;
    }

    return n;
}

int fflib_search_file_end(const char *ffdir,
                          const char *file_end,
                          gmx_bool    bFatalError,
                          char     ***filenames)
{
    return low_fflib_search_file_end(ffdir, FALSE, file_end, bFatalError,
                                     filenames, NULL);
}

int fflib_search_file_in_dirend(const char *filename, const char *dirend,
                                char ***dirnames)
{
    int             nf, i;
    char          **f, **f_short;
    int             n;
    char          **dns;
    gmx_directory_t dirhandle;
    char            nextname[STRLEN];
    int             rc;

    /* Find all files (not only dir's) ending on dirend */
    nf = low_fflib_search_file_end(NULL, TRUE, dirend, FALSE, &f, &f_short);

    n   = 0;
    dns = NULL;
    for (i = 0; i < nf; i++)
    {
        rc = gmx_directory_open(&dirhandle, f[i]);

        if (rc == 0)
        {
            while (gmx_directory_nextfile(dirhandle, nextname, STRLEN-1) == 0)
            {
                nextname[STRLEN-1] = 0;
                if (strcmp(nextname, filename) == 0)
                {
                    /* We have a match */
                    srenew(dns, n+1);
                    dns[n] = strdup(f_short[i]);
                    n++;
                }
            }
            gmx_directory_close(dirhandle);
        }
        sfree(f[i]);
        sfree(f_short[i]);
    }
    sfree(f);
    sfree(f_short);

    *dirnames = dns;

    return n;
}

gmx_bool fflib_fexist(const char *file)
{
    char *file_fullpath;

    file_fullpath = low_gmxlibfn(file, TRUE, FALSE);

    if (file_fullpath == NULL)
    {
        return FALSE;
    }
    else
    {
        sfree(file_fullpath);

        return TRUE;
    }
}


FILE *fflib_open(const char *file)
{
    char *file_fullpath;
    FILE *fp;

    file_fullpath = gmxlibfn(file);
    fprintf(stderr, "Opening force field file %s\n", file_fullpath);
    fp = ffopen(file_fullpath, "r");
    sfree(file_fullpath);

    return fp;
}
