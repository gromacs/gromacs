/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define HAVE_DIRENT
#ifdef HAVE_DIRENT
#include <dirent.h>
#endif

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <direct.h>
#include <io.h>
#endif

#include "sysstuff.h"
#include "string2.h"
#include "futil.h"
#include "network.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "statutil.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif


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

void fflib_filename_base(const char *filename,char *filebase,int maxlen)
{
    const char *cptr;
    char *ptr;

    cptr = strrchr(filename,DIR_SEPARATOR);
    if (cptr != NULL)
    {
        /* Skip the separator */
        cptr += 1;
    }
    else
    {
        cptr = filename;
    }
    if (strlen(filename) >= maxlen)
    {
        gmx_fatal(FARGS,"filename is longer (%d) than maxlen (%d)",
                  strlen(filename),maxlen);
    }
    strcpy(filebase,cptr);
    /* Remove the extension */
    ptr = strrchr(filebase,'.');
    if (ptr != NULL)
    {
        ptr[0] = '\0';
    }
}

static void sort_filenames(int n,char **name)
{
    /* Slow sort, but we usually have tens of names */
    int  i,j,f;
    char *tmp;

    for(i=0; i<n-1; i++)
    {
        f = i;
        for(j=i+1; j<n; j++)
        {
            if (strcmp(name[j],name[f]) < 0)
            {
                f = j;
            }
        }
        if (f > i)
        {
            tmp     = name[i];
            name[i] = name[f];
            name[f] = tmp;
        }
    }
}

int fflib_search_file_end(const char *ffdir,const char *file_end,
                          bool bFatalError,
                          char ***filenames)
{
#ifndef HAVE_DIRENT
    gmx_fatal(FARGS,"lib_search_file_end called while the 'dirent' functionality is not available on this system");
    return 0;
#else
    char *ret=NULL;
    char *lib,*dir;
    char buf[1024];
    char libpath[GMX_PATH_MAX];
    bool env_is_set=FALSE;
    int  len_fe,len_name;
    char **fns;
    char *s,tmppath[GMX_PATH_MAX],tmpfn[GMX_PATH_MAX];
    DIR  *dirptr;
    struct dirent *dirent;
    int  n,n_thisdir;

    /* GMXLIB can be a path now */
    lib=getenv("GMXLIB");
    if (lib != NULL)
    {
        env_is_set = TRUE;
        strncpy(libpath,lib,GMX_PATH_MAX);
    } 
    else if (!get_libdir(libpath))
    {
        strncpy(libpath,GMXLIBDIR,GMX_PATH_MAX);
    }

    len_fe = strlen(file_end);

    if (ffdir != NULL)
    {
        /* Search in current dir and ffdir */
        strncpy(tmppath,ffdir,GMX_PATH_MAX);
    }
    else
    {
        /* Search in current dir and libpath */
        strncpy(tmppath,libpath,GMX_PATH_MAX);
    }
    s = tmppath;
    n = 0;
    fns = NULL;
    /* Start with the current directory, continue with libpath */
    dir = ".";
    do
    {
        dirptr = opendir(dir);
        if (dirptr != NULL)
        {
            n_thisdir = 0;
            while ((dirent = readdir(dirptr)) != NULL)
            {
                if (debug)
                {
                    fprintf(debug,"dir '%s' file '%s'\n",dir,dirent->d_name);
                }
                len_name = strlen(dirent->d_name);
                /* What about case sensitivity? */
                if (len_name >= len_fe &&
                    strcmp(dirent->d_name+len_name-len_fe,file_end) == 0)
                {
                    /* We have a match */
                    srenew(fns,n+1);
                    if (strcmp(dir,".") == 0)
                    {
                        fns[n] = strdup(dirent->d_name);
                    }
                    else
                    {
                        sprintf(tmpfn,"%s%c%s",
                                dir,DIR_SEPARATOR,dirent->d_name);
                        fns[n] = strdup(tmpfn);
                    }
                    n++;
                    n_thisdir++;
                }
            }

            sort_filenames(n_thisdir,fns+n-n_thisdir);
        }
    }
    while((dir=gmx_strsep(&s, PATH_SEPARATOR)) != NULL);

    if (n == 0 && bFatalError)
    {
        if (ffdir != NULL)
        {
            gmx_fatal(FARGS,"Could not find any files ending on '%s' in the force field directory '%s'",file_end,ffdir);
        }
        else
        {
            gmx_fatal(FARGS,"Could not find any files ending on '%s' in the current directory or the GROMACS library search path",file_end);
        }
    }

    *filenames = fns;

    return n;
#endif
}

int fflib_search_file_in_dirend(const char *filename,const char *dirend,
                                char ***dirnames)
{
#ifndef HAVE_DIRENT
    gmx_fatal(FARGS,"lib_search_file_in_dirend called while the 'dirent' functionality is not available on this system");
    return 0;
#else
    int  nf,i;
    char **f;
    int  n;
    char **dns;
    DIR  *dirptr;
    struct dirent *dirent;

    /* Find all files (not only dir's) ending on dirend */
    nf = fflib_search_file_end(NULL,dirend,FALSE,&f);

    n = 0;
    dns = NULL;
    for(i=0; i<nf; i++)
    {
        dirptr = opendir(f[i]);
        if (dirptr != NULL)
        {
            while ((dirent = readdir(dirptr)) != NULL)
            {
                if (strcmp(dirent->d_name,filename) == 0)
                {
                    /* We have a match */
                    srenew(dns,n+1);
                    dns[n] = strdup(f[i]);
                    n++;
                }
            }
        }
        sfree(f[i]);
    }
    sfree(f);

    *dirnames = dns;

    return n;
#endif 
}

FILE *fflib_open(const char *file)
{
    FILE *fp;

    fprintf(stderr,"Opening force field file %s\n",file);

    return ffopen(file,"r");
}
