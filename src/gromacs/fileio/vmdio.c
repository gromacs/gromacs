/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "vmdio.h"

#include "config.h"

/* Derived from PluginMgr.C and catdcd.c */

/* PluginMgr.C: Copyright: */
/***************************************************************************
 * cr
 * cr            (C) Copyright 1995-2009 The Board of Trustees of the
 * cr                        University of Illinois
 * cr                         All Rights Reserved
 * cr
   Developed by:           Theoretical and Computational Biophysics Group
                        University of Illinois at Urbana-Champaign
                        http://www.ks.uiuc.edu/

   Permission is hereby granted, free of charge, to any person obtaining a copy of
   this software and associated documentation files (the Software), to deal with
   the Software without restriction, including without limitation the rights to
   use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
   of the Software, and to permit persons to whom the Software is furnished to
   do so, subject to the following conditions:

   Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimers.

   Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimers in the documentation
   and/or other materials provided with the distribution.

   Neither the names of Theoretical and Computational Biophysics Group,
   University of Illinois at Urbana-Champaign, nor the names of its contributors
   may be used to endorse or promote products derived from this Software without
   specific prior written permission.

   THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
   THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
   OTHER DEALINGS WITH THE SOFTWARE.
 ***************************************************************************/

/* catdcd.c: Copyright: */
/*****************************************************************************/
/*                                                                           */
/* (C) Copyright 2001-2005 Justin Gullingsrud and the University of Illinois.*/
/*                                                                           */
/*****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Plugin header files; get plugin source from www.ks.uiuc.edu/Research/vmd"
 */
#include "external/vmd_molfile/molfile_plugin.h"
#include "external/vmd_molfile/vmddlopen.h"
#ifndef GMX_NATIVE_WINDOWS
#include <glob.h>
#else
#ifndef _WIN32_IE
#define _WIN32_IE 0x0500 /* SHGetFolderPath is available since WinXP/IE5 */
#endif
#include <shlobj.h>
#include <windows.h>
#endif

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


typedef int (*initfunc)(void);
typedef int (*regfunc)(void *, vmdplugin_register_cb);
typedef int (*finifunc)(void);



static int register_cb(void *v, vmdplugin_t *p)
{
    const char     *key       = p->name;
    t_gmxvmdplugin *vmdplugin = (t_gmxvmdplugin*)v;

    if (strcmp(key, vmdplugin->filetype) == 0)
    {
        vmdplugin->api = (molfile_plugin_t *)p;
    }
    return VMDPLUGIN_SUCCESS;
}

static int load_sharedlibrary_plugins(const char *fullpath, t_gmxvmdplugin* vmdplugin)
{
    /* Open the dll; try to execute the init function. */
    void *handle, *ifunc, *registerfunc;
    handle = vmddlopen(fullpath);
    if (!handle)
    {
        if (debug)
        {
            fprintf(debug, "\nUnable to open dynamic library %s.\n%s\n",  fullpath, vmddlerror());         /*only to debug because of stdc++ erros */
        }
        return 0;
    }

    ifunc = vmddlsym(handle, "vmdplugin_init");
    if (!ifunc || ((initfunc)(ifunc))())
    {
        printf("\nvmdplugin_init() for %s returned an error; plugin(s) not loaded.\n", fullpath);
        vmddlclose(handle);
        return 0;
    }

    registerfunc = vmddlsym(handle, "vmdplugin_register");
    if (!registerfunc)
    {
        printf("\nDidn't find the register function in %s; plugin(s) not loaded.\n", fullpath);
        vmddlclose(handle);
        return 0;
    }
    else
    {
        /* Load plugins from the library.*/
        ((regfunc)registerfunc)(vmdplugin, register_cb);
    }

    /* in case this library does not support the filetype, close it */
    if (vmdplugin->api == NULL)
    {
        vmddlclose(handle);
    }

    return 1;
}

/*return: 1: success, 0: last frame, -1: error*/
gmx_bool read_next_vmd_frame(t_trxframe *fr)
{
    int                rc, i;
    rvec               vec, angle;
    molfile_timestep_t ts;


    fr->bV = fr->vmdplugin->bV;

#ifdef GMX_DOUBLE
    snew(ts.coords, fr->natoms*3);
    if (fr->bV)
    {
        snew(ts.velocities, fr->natoms*3);
    }
#else
    ts.coords = (float*)fr->x;
    if (fr->bV)
    {
        ts.velocities = (float*)fr->v;
    }
#endif

    rc = fr->vmdplugin->api->read_next_timestep(fr->vmdplugin->handle, fr->natoms, &ts);

    if (rc < -1)
    {
        fprintf(stderr, "\nError reading input file (error code %d)\n", rc);
    }
    if (rc < 0)
    {
        fr->vmdplugin->api->close_file_read(fr->vmdplugin->handle);
        return 0;
    }

#ifdef GMX_DOUBLE
    for (i = 0; i < fr->natoms; i++)
    {
        fr->x[i][0] = .1*ts.coords[i*3];
        fr->x[i][1] = .1*ts.coords[i*3+1];
        fr->x[i][2] = .1*ts.coords[i*3+2];
        if (fr->bV)
        {
            fr->v[i][0] = .1*ts.velocities[i*3];
            fr->v[i][1] = .1*ts.velocities[i*3+1];
            fr->v[i][2] = .1*ts.velocities[i*3+2];
        }
    }
    sfree(ts.coords);
    if (fr->bV)
    {
        sfree(ts.velocities);
    }
#else
    for (i = 0; i < fr->natoms; i++)
    {
        svmul(.1, fr->x[i], fr->x[i]);
        if (fr->bV)
        {
            svmul(.1, fr->v[i], fr->v[i]);
        }
    }
#endif

    fr->bX   = 1;
    fr->bBox = 1;
    vec[0]   = .1*ts.A; vec[1] = .1*ts.B; vec[2] = .1*ts.C;
    angle[0] = ts.alpha; angle[1] = ts.beta; angle[2] = ts.gamma;
    matrix_convert(fr->box, vec, angle);
    if (fr->vmdplugin->api->abiversion > 10)
    {
        fr->bTime = TRUE;
        fr->time  = ts.physical_time;
    }
    else
    {
        fr->bTime = FALSE;
    }


    return 1;
}

static int load_vmd_library(const char *fn, t_gmxvmdplugin *vmdplugin)
{
    char            pathname[GMX_PATH_MAX], filename[GMX_PATH_MAX];
    const char     *pathenv;
    const char     *err;
    int             i;
    int             ret = 0;
    char            pathenv_buffer[GMX_PATH_MAX];
#ifndef GMX_NATIVE_WINDOWS
    glob_t          globbuf;
    const char     *defpath_suffix = "/plugins/*/molfile";
    const char     *defpathenv     = GMX_VMD_PLUGIN_PATH;
#else
    WIN32_FIND_DATA ffd;
    HANDLE          hFind = INVALID_HANDLE_VALUE;
    char            progfolder[GMX_PATH_MAX];
    char            defpathenv[GMX_PATH_MAX];
    const char     *defpath_suffix = "\\plugins\\WIN32\\molfile";
    SHGetFolderPath(NULL, CSIDL_PROGRAM_FILES, NULL, SHGFP_TYPE_CURRENT, progfolder);
    sprintf(defpathenv, "%s\\University of Illinois\\VMD\\plugins\\WIN32\\molfile", progfolder);
#endif

    vmdplugin->api      = NULL;
    vmdplugin->filetype = strrchr(fn, '.');
    if (!vmdplugin->filetype)
    {
        return 0;
    }
    vmdplugin->filetype++;

    /* First look for an explicit path given at run time for the
     * plugins, then an implicit run-time path, and finally for one
     * given at configure time. This last might be hard-coded to the
     * default for VMD installs. */
    pathenv = getenv("VMD_PLUGIN_PATH");
    if (pathenv == NULL)
    {
        pathenv = getenv("VMDDIR");
        if (NULL == pathenv)
        {
            printf("\nNeither VMD_PLUGIN_PATH or VMDDIR set. ");
            printf("Using default location:\n%s\n", defpathenv);
            pathenv = defpathenv;
        }
        else
        {
            printf("\nVMD_PLUGIN_PATH no set, but VMDDIR is set. ");
#ifdef _MSC_VER
            _snprintf_s(pathenv_buffer, sizeof(pathenv_buffer), _TRUNCATE, "%s%s", pathenv, defpath_suffix);
#else
            snprintf(pathenv_buffer, sizeof(pathenv_buffer), "%s%s", pathenv, defpath_suffix);
#endif
            printf("Using semi-default location:\n%s\n", pathenv_buffer);
            pathenv = pathenv_buffer;
        }
    }
    strncpy(pathname, pathenv, sizeof(pathname));
#ifndef GMX_NATIVE_WINDOWS
    strcat(pathname, "/*.so");
    glob(pathname, 0, NULL, &globbuf);
    if (globbuf.gl_pathc == 0)
    {
        printf("\nNo VMD Plugins found\n"
               "Set the environment variable VMD_PLUGIN_PATH to the molfile folder within the\n"
               "VMD installation.\n"
               "The architecture (e.g. 32bit versus 64bit) of GROMACS and VMD has to match.\n");
        return 0;
    }
    for (i = 0; i < globbuf.gl_pathc && vmdplugin->api == NULL; i++)
    {
        /* FIXME: Undefined which plugin is chosen if more than one plugin
           can read a certain file ending. Requires some additional command
           line option or enviroment variable to specify which plugin should
           be picked.
         */
        ret |= load_sharedlibrary_plugins(globbuf.gl_pathv[i], vmdplugin);
    }
    globfree(&globbuf);
#else
    strcat(pathname, "\\*.so");
    hFind = FindFirstFile(pathname, &ffd);
    if (INVALID_HANDLE_VALUE == hFind)
    {
        printf("\nNo VMD Plugins found\n");
        return 0;
    }
    do
    {
        sprintf(filename, "%s\\%s", pathenv, ffd.cFileName);
        ret |= load_sharedlibrary_plugins(filename, vmdplugin);
    }
    while (FindNextFile(hFind, &ffd )  != 0 && vmdplugin->api == NULL);
    FindClose(hFind);
#endif

    if (!ret)
    {
        printf("\nCould not open any VMD library.\n");
        err = vmddlerror();
        if (!err)
        {
            printf("Compiled with dlopen?\n");
        }
        else
        {
            printf("Last error:\n%s\n", err);
        }
        return 0;
    }

    if (vmdplugin->api == NULL)
    {
        printf("\nNo plugin for %s found\n", vmdplugin->filetype);
        return 0;
    }

    if (vmdplugin->api->abiversion < 10)
    {
        printf("\nPlugin and/or VMD is too old. At least VMD 1.8.6 is required.\n");
        return 0;
    }

    printf("\nUsing VMD plugin: %s (%s)\n", vmdplugin->api->name, vmdplugin->api->prettyname);

    return 1;

}

int read_first_vmd_frame(const char *fn, t_trxframe *fr)
{
    molfile_timestep_metadata_t *metadata = NULL;

    snew(fr->vmdplugin, 1);
    if (!load_vmd_library(fn, fr->vmdplugin))
    {
        return 0;
    }

    fr->vmdplugin->handle = fr->vmdplugin->api->open_file_read(fn, fr->vmdplugin->filetype, &fr->natoms);

    if (!fr->vmdplugin->handle)
    {
        fprintf(stderr, "\nError: could not open file '%s' for reading.\n",
                fn);
        return 0;
    }

    if (fr->natoms == MOLFILE_NUMATOMS_UNKNOWN)
    {
        fprintf(stderr, "\nFormat of file %s does not record number of atoms.\n", fn);
        return 0;
    }
    else if (fr->natoms == MOLFILE_NUMATOMS_NONE)
    {
        fprintf(stderr, "\nNo atoms found by VMD plugin in file %s.\n", fn );
        return 0;
    }
    else if (fr->natoms < 1)     /*should not be reached*/
    {
        fprintf(stderr, "\nUnknown number of atoms %d for VMD plugin opening file %s.\n",
                fr->natoms, fn );
        return 0;
    }

    snew(fr->x, fr->natoms);

    fr->vmdplugin->bV = 0;
    if (fr->vmdplugin->api->abiversion > 10 && fr->vmdplugin->api->read_timestep_metadata)
    {
        fr->vmdplugin->api->read_timestep_metadata(fr->vmdplugin->handle, metadata);
        assert(metadata);
        fr->vmdplugin->bV = metadata->has_velocities;
        if (fr->vmdplugin->bV)
        {
            snew(fr->v, fr->natoms);
        }
    }
    else
    {
        fprintf(stderr,
                "\nThis trajectory is being read with a VMD plug-in from before VMD"
                "\nversion 1.8, or from a trajectory that lacks time step metadata."
                "\nEither way, GROMACS cannot tell whether the trajectory has velocities.\n");
    }
    return 1;

}
