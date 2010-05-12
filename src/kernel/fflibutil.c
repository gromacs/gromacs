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

#ifdef HAVE_DIRENT_H
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

static void sort_filenames(int n,char **name,char **name2)
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
                                     const char *file_end,
                                     bool bFatalError,
                                     char ***filenames,
                                     char ***filenames_short)
{
#ifndef HAVE_DIRENT_H
    gmx_fatal(FARGS,"lib_search_file_end called while the 'dirent' functionality is not available on this system");
    return 0;
#else
    char *ret=NULL;
    char *lib,*dir;
    char buf[1024];
    char *libpath;
    bool env_is_set;
    int  len_fe,len_name;
    char **fns,**fns_short;
    char dir_print[GMX_PATH_MAX];
    char *pdum;
    char *s,fn_dir[GMX_PATH_MAX];
    DIR  *dirptr;
    struct dirent *dirent;
    int  n,n_thisdir;

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
        snew(libpath,GMX_PATH_MAX);
        if (lib != NULL)
        {
            env_is_set = TRUE;
            strncpy(libpath,lib,GMX_PATH_MAX);
        } 
        else if (!get_libdir(libpath))
        {
            strncpy(libpath,GMXLIBDIR,GMX_PATH_MAX);
        }
    }
    s = libpath;
    n = 0;
    fns       = NULL;
    fns_short = NULL;
    /* Start with the current directory, continue with libpath */
    dir = ".";
    do
    {
        dirptr = opendir(dir);
        if (dirptr != NULL)
        {
            if (strcmp(dir,".") == 0)
            {
                /* Print the absolute path to the current working dir. */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
                pdum = _getcwd(dir_print,sizeof(dir_print)-1);
#else
                pdum =  getcwd(dir_print,sizeof(dir_print)-1);
#endif
            }
            else
            {
                /* Print the directory dir, can be relative or absolute. */
                strcpy(dir_print,dir);
            }

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
                    sprintf(fn_dir,"%s%c%s",
                            dir_print,DIR_SEPARATOR,dirent->d_name);

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
                        srenew(fns_short,n+1);
                        if (strcmp(dir,".") == 0 || env_is_set)
                        {
                            fns_short[n] = strdup(fn_dir);
                        }
                        else
                        {
                            fns_short[n] = strdup(dirent->d_name);
                        }
                    }
                    n++;
                    n_thisdir++;
                }
            }

            sort_filenames(n_thisdir,
                           fns+n-n_thisdir,
                           fns_short==NULL ? NULL : fns_short+n-n_thisdir);
        }
    }
    while((dir=gmx_strsep(&s, PATH_SEPARATOR)) != NULL);

    sfree(libpath);

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
    if (ffdir == NULL)
    {
        *filenames_short = fns_short;
    }

    return n;
#endif
}

int fflib_search_file_end(const char *ffdir,const char *file_end,
                          bool bFatalError,
                          char ***filenames)
{
    return low_fflib_search_file_end(ffdir,file_end,bFatalError,filenames,NULL);
}

int fflib_search_file_in_dirend(const char *filename,const char *dirend,
                                char ***dirnames)
{
#ifndef HAVE_DIRENT_H
    gmx_fatal(FARGS,"lib_search_file_in_dirend called while the 'dirent' functionality is not available on this system");
    return 0;
#else
    int  nf,i;
    char **f,**f_short;
    int  n;
    char **dns;
    DIR  *dirptr;
    struct dirent *dirent;

    /* Find all files (not only dir's) ending on dirend */
    nf = low_fflib_search_file_end(NULL,dirend,FALSE,&f,&f_short);

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
                    dns[n] = strdup(f_short[i]);
                    n++;
                }
            }
        }
        sfree(f[i]);
        sfree(f_short[i]);
    }
    sfree(f);
    sfree(f_short);

    *dirnames = dns;

    return n;
#endif 
}

bool fflib_fexist(const char *file)
{
    char *file_fullpath;

    file_fullpath = low_gmxlibfn(file,FALSE);
    
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
    fprintf(stderr,"Opening force field file %s\n",file_fullpath);
    fp = ffopen(file_fullpath,"r");
    sfree(file_fullpath);

    return fp;
}
