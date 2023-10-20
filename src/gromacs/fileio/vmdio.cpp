/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "vmdio.h"

#include "config.h"

#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/*
 * Plugin header files; get plugin source from www.ks.uiuc.edu/Research/vmd"
 */
#include "external/vmd_molfile/molfile_plugin.h"
#include "external/vmd_molfile/vmddlopen.h"
#if !GMX_NATIVE_WINDOWS
#    include <glob.h>
#else
#    ifndef _WIN32_IE
#        define _WIN32_IE 0x0500 /* SHGetFolderPath is available since WinXP/IE5 */
#    endif
#    include <shlobj.h>
#    include <windows.h>
#endif

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


typedef int (*initfunc)();
typedef int (*regfunc)(void*, vmdplugin_register_cb);
typedef int (*finifunc)();


static int register_cb(void* v, vmdplugin_t* p)
{
    const char*      key       = p->name;
    gmx_vmdplugin_t* vmdplugin = static_cast<gmx_vmdplugin_t*>(v);

    if (strcmp(key, vmdplugin->filetype.string().c_str()) == 0)
    {
        vmdplugin->api = reinterpret_cast<molfile_plugin_t*>(p);
    }
    return VMDPLUGIN_SUCCESS;
}

static int load_sharedlibrary_plugins(const char* fullpath, gmx_vmdplugin_t* vmdplugin)
{
    /* Open the dll; try to execute the init function. */
    void *handle, *ifunc, *registerfunc;
    handle = vmddlopen(fullpath);
    if (!handle)
    {
        if (debug)
        {
            fprintf(debug, "\nUnable to open dynamic library %s.\n%s\n", fullpath, vmddlerror()); /*only to debug because of stdc++ erros */
        }
        return 0;
    }

    ifunc = vmddlsym(handle, "vmdplugin_init");
    if (!ifunc || (reinterpret_cast<initfunc>(ifunc))())
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
        (reinterpret_cast<regfunc>(registerfunc))(vmdplugin, register_cb);
    }

    /* in case this library does not support the filetype, close it */
    if (vmdplugin->api == nullptr)
    {
        vmddlclose(handle);
    }

    return 1;
}

/*return: 1: success, 0: last frame, -1: error*/
gmx_bool read_next_vmd_frame(gmx_vmdplugin_t* vmdplugin, t_trxframe* fr)
{
    int                rc, i;
    rvec               vec, angle;
    molfile_timestep_t ts;


    fr->bV = vmdplugin->bV;

#if GMX_DOUBLE
    snew(ts.coords, fr->natoms * 3);
    if (fr->bV)
    {
        snew(ts.velocities, fr->natoms * 3);
    }
#else
    ts.coords = reinterpret_cast<float*>(fr->x);
    if (fr->bV)
    {
        ts.velocities = reinterpret_cast<float*>(fr->v);
    }
#endif

    rc = vmdplugin->api->read_next_timestep(vmdplugin->handle, fr->natoms, &ts);

    if (rc < -1)
    {
        fprintf(stderr, "\nError reading input file (error code %d)\n", rc);
    }
    if (rc < 0)
    {
        vmdplugin->api->close_file_read(vmdplugin->handle);
        return false;
    }

#if GMX_DOUBLE
    for (i = 0; i < fr->natoms; i++)
    {
        fr->x[i][0] = .1 * ts.coords[i * 3];
        fr->x[i][1] = .1 * ts.coords[i * 3 + 1];
        fr->x[i][2] = .1 * ts.coords[i * 3 + 2];
        if (fr->bV)
        {
            fr->v[i][0] = .1 * ts.velocities[i * 3];
            fr->v[i][1] = .1 * ts.velocities[i * 3 + 1];
            fr->v[i][2] = .1 * ts.velocities[i * 3 + 2];
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

    fr->bX   = true;
    fr->bBox = true;
    vec[0]   = .1 * ts.A;
    vec[1]   = .1 * ts.B;
    vec[2]   = .1 * ts.C;
    angle[0] = ts.alpha;
    angle[1] = ts.beta;
    angle[2] = ts.gamma;
    matrix_convert(fr->box, vec, angle);
    if (vmdplugin->api->abiversion > 10)
    {
        fr->bTime = TRUE;
        fr->time  = ts.physical_time;
    }
    else
    {
        fr->bTime = FALSE;
    }


    return true;
}

static int load_vmd_library(const std::filesystem::path& fn, gmx_vmdplugin_t* vmdplugin)
{
    const char* err;
    int         ret = 0;
#if !GMX_NATIVE_WINDOWS
    glob_t                      globbuf;
    const std::string           defpath_suffix = "/plugins/*/molfile";
    const std::filesystem::path defpathenv     = GMX_VMD_PLUGIN_PATH;
#else
    WIN32_FIND_DATA       ffd;
    HANDLE                hFind = INVALID_HANDLE_VALUE;
    char                  progfolder[GMX_PATH_MAX];
    std::filesystem::path defpath_suffix = "\\plugins\\WIN32\\molfile";
    SHGetFolderPath(NULL, CSIDL_PROGRAM_FILES, NULL, SHGFP_TYPE_CURRENT, progfolder);
    std::filesystem::path defpathenv =
            progfolder / "University of Illinois" / "VMD" / "plugins" / "WIN32" / "molfile";
#endif

    vmdplugin->api      = nullptr;
    vmdplugin->filetype = fn.extension();
    if (!fn.has_extension())
    {
        return 0;
    }
    vmdplugin->filetype = vmdplugin->filetype.string().substr(1);

    /* First look for an explicit path given at run time for the
     * plugins, then an implicit run-time path, and finally for one
     * given at configure time. This last might be hard-coded to the
     * default for VMD installs. */
    const char*           pathEnvChar = getenv("VMD_PLUGIN_PATH");
    std::filesystem::path pathenv     = pathEnvChar != nullptr ? pathEnvChar : "";
    std::filesystem::path fallBackPathEnv;
    if (pathenv.empty())
    {
        pathenv = getenv("VMDDIR");
        if (pathenv.empty())
        {
            printf("\nNeither VMD_PLUGIN_PATH or VMDDIR set. ");
            printf("Using default location:\n%s\n", defpathenv.c_str());
            pathenv = defpathenv;
        }
        else
        {
            printf("\nVMD_PLUGIN_PATH no set, but VMDDIR is set. ");
            fallBackPathEnv = pathenv / defpath_suffix;
            pathenv         = fallBackPathEnv;
            printf("Using semi-default location:\n%s\n", pathenv.c_str());
        }
    }
#if !GMX_NATIVE_WINDOWS
    auto pathname = pathenv / "/*.so";
    glob(pathname.c_str(), 0, nullptr, &globbuf);
    if (globbuf.gl_pathc == 0)
    {
        printf("\nNo VMD Plugins found\n"
               "Set the environment variable VMD_PLUGIN_PATH to the molfile folder within the\n"
               "VMD installation.\n"
               "The architecture (e.g. 32bit versus 64bit) of GROMACS and VMD has to match.\n");
        return 0;
    }
    for (size_t i = 0; i < globbuf.gl_pathc && vmdplugin->api == nullptr; i++)
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
    auto pathname = pathenv / "*.so";
    hFind         = FindFirstFile(pathname.c_str(), &ffd);
    if (INVALID_HANDLE_VALUE == hFind)
    {
        printf("\nNo VMD Plugins found\n");
        return 0;
    }
    do
    {
        auto filename = pathenv / ffd.cFileName;
        ret |= load_sharedlibrary_plugins(filename.c_str(), vmdplugin);
    } while (FindNextFile(hFind, &ffd) != 0 && vmdplugin->api == NULL);
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

    if (vmdplugin->api == nullptr)
    {
        printf("\nNo plugin for %s found\n", vmdplugin->filetype.c_str());
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

int read_first_vmd_frame(const std::filesystem::path& fn, gmx_vmdplugin_t** vmdpluginp, t_trxframe* fr)
{
    molfile_timestep_metadata_t* metadata = nullptr;
    gmx_vmdplugin_t*             vmdplugin;

    vmdplugin   = new gmx_vmdplugin_t;
    *vmdpluginp = vmdplugin;
    if (!load_vmd_library(fn, vmdplugin))
    {
        return 0;
    }

    vmdplugin->handle =
            vmdplugin->api->open_file_read(fn.c_str(), vmdplugin->filetype.c_str(), &fr->natoms);

    if (!vmdplugin->handle)
    {
        fprintf(stderr, "\nError: could not open file '%s' for reading.\n", fn.c_str());
        return 0;
    }

    if (fr->natoms == MOLFILE_NUMATOMS_UNKNOWN)
    {
        fprintf(stderr, "\nFormat of file %s does not record number of atoms.\n", fn.c_str());
        return 0;
    }
    else if (fr->natoms == MOLFILE_NUMATOMS_NONE)
    {
        fprintf(stderr, "\nNo atoms found by VMD plugin in file %s.\n", fn.c_str());
        return 0;
    }
    else if (fr->natoms < 1) /*should not be reached*/
    {
        fprintf(stderr,
                "\nUnknown number of atoms %d for VMD plugin opening file %s.\n",
                fr->natoms,
                fn.c_str());
        return 0;
    }

    snew(fr->x, fr->natoms);

    vmdplugin->bV = false;
    if (vmdplugin->api->abiversion > 10 && vmdplugin->api->read_timestep_metadata)
    {
        vmdplugin->api->read_timestep_metadata(vmdplugin->handle, metadata);
        assert(metadata);
        vmdplugin->bV = (metadata->has_velocities != 0);
        if (vmdplugin->bV)
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
