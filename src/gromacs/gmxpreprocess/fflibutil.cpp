/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012,2013,2014, by the GROMACS development team, led by
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
#include "network.h"
#include "gmx_fatal.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "fflibutil.h"

#include "gromacs/fileio/futil.h"
#include "gromacs/fileio/path.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programcontext.h"

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
    char **fns = NULL, **fns_short = NULL;
    int    n   = 0;
    try
    {
        std::vector<std::string> libPaths;
        bool                     bEnvIsSet = false;

        if (ffdir != NULL)
        {
            /* Search ffdir in current dir and library dirs */
            libPaths.push_back(gmxlibfn(ffdir));
        }
        else
        {
            /* GMXLIB can be a path now */
            if (bAddCWD)
            {
                libPaths.push_back(".");
            }
            const char *lib = getenv("GMXLIB");
            if (lib != NULL)
            {
                bEnvIsSet = true;
                gmx::Path::splitPathEnvironment(lib, &libPaths);
            }
            else
            {
                libPaths.push_back(gmx::getProgramContext().defaultLibraryDataPath());
            }
        }

        const int len_fe = strlen(file_end);

        std::vector<std::string>::const_iterator i;
        for (i = libPaths.begin(); i != libPaths.end(); ++i)
        {
            const char      *dir = i->c_str();
            gmx_directory_t  dirhandle;
            const int        rc  = gmx_directory_open(&dirhandle, dir);
            if (rc == 0)
            {
                char nextname[STRLEN];
                int  n_thisdir = 0;
                while (gmx_directory_nextfile(dirhandle, nextname, STRLEN-1) == 0)
                {
                    nextname[STRLEN-1] = 0;
                    if (debug)
                    {
                        fprintf(debug, "dir '%s' %d file '%s'\n",
                                dir, n_thisdir, nextname);
                    }
                    const int len_name = strlen(nextname);
                    /* What about case sensitivity? */
                    if (len_name >= len_fe &&
                        strcmp(nextname+len_name-len_fe, file_end) == 0)
                    {
                        char fn_dir[GMX_PATH_MAX];
                        /* We have a match */
                        srenew(fns, n+1);
                        sprintf(fn_dir, "%s%c%s", dir, DIR_SEPARATOR, nextname);

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
                            if (strcmp(dir, ".") == 0 || bEnvIsSet)
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
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

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
    fp = gmx_ffopen(file_fullpath, "r");
    sfree(file_fullpath);

    return fp;
}
