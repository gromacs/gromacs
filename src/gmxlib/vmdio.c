/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

/* Derived from PluginMgr.C and catdcd.c */

/* PluginMgr.C: Copyright: */
/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * LICENSE:
 *   UIUC Open Source License
 *   http://www.ks.uiuc.edu/Research/vmd/plugins/pluginlicense.html
 *
 ***************************************************************************/

/* catdcd.c: Copyright: */
/*****************************************************************************/
/*                                                                           */
/* (C) Copyright 2001-2005 Justin Gullingsrud and the University of Illinois.*/
/*                                                                           */
/*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* 
 * Plugin header files; get plugin source from www.ks.uiuc.edu/Research/vmd"
 */
#include "molfile_plugin.h"
#include "vmddlopen.h"
#if !((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include "glob.h"
#else
#include <windows.h>
#include <shlobj.h>
#endif
#include "smalloc.h"
#include "futil.h"


#include "types/simple.h"
#include "vec.h"
#include "gmxfio.h"


typedef int (*initfunc)(void);
typedef int (*regfunc)(void *, vmdplugin_register_cb);
typedef int (*finifunc)(void);



static int register_cb(void *v, vmdplugin_t *p) {
    const char *key = p->name;
    t_gmxvmdplugin *vmdplugin = (t_gmxvmdplugin*)v;

    if (strcmp(key,vmdplugin->filetype)==0)
    {
        vmdplugin->api = (molfile_plugin_t *)p;
    }
    return VMDPLUGIN_SUCCESS;
}

static int load_sharedlibrary_plugins(const char *fullpath,t_gmxvmdplugin* vmdplugin) {
    // Open the dll; try to execute the init function.
    void *handle, *ifunc, *registerfunc; 
    handle = vmddlopen(fullpath);
    if (!handle) {
        if (debug) fprintf(debug, "\nUnable to open dynamic library %s.\n%s\n",  fullpath, vmddlerror());  /*only to debug because of stdc++ erros */
        return 0;
    }

    ifunc = vmddlsym(handle, "vmdplugin_init");
    if (ifunc && ((initfunc)(ifunc))()) {
        printf("\nvmdplugin_init() for %s returned an error; plugin(s) not loaded.\n", fullpath);
        vmddlclose(handle);
        return 0;
    }

    registerfunc = vmddlsym(handle, "vmdplugin_register");
    if (!registerfunc) {
        printf("\nDidn't find the register function in %s; plugin(s) not loaded.\n", fullpath);
        vmddlclose(handle);
        return 0;
    } else {
        // Load plugins from the library.
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
bool read_next_vmd_frame(int status,t_trxframe *fr)
{
    int rc,i;
    rvec vec, angle;
    molfile_timestep_t ts;


    fr->bV = fr->vmdplugin.bV; 
        
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

    rc = fr->vmdplugin.api->read_next_timestep(fr->vmdplugin.handle, fr->natoms, &ts);

    if (rc < -1) {
        fprintf(stderr, "\nError reading input file (error code %d)\n", rc);
    }
    if (rc < 0)
    {
        fr->vmdplugin.api->close_file_read(fr->vmdplugin.handle);
        return 0;
    }

#ifdef GMX_DOUBLE
    for (i=0;i<fr->natoms;i++)
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
    for (i=0;i<fr->natoms;i++)
    {
        svmul(.1,fr->x[i],fr->x[i]);
        if (fr->bV)
        {
            svmul(.1,fr->v[i],fr->v[i]);
        }
    }
#endif

    fr->bX = 1;
    vec[0] = .1*ts.A; vec[1] = .1*ts.B; vec[2] = .1*ts.B;
    angle[0] = ts.alpha; angle[1] = ts.beta; angle[2] = ts.gamma; 
    matrix_convert(fr->box,vec,angle);
    fr->bTime = 1;
    fr->time = ts.physical_time;



    return 1;
}


int read_first_vmd_frame(int *status,const char *fn,t_trxframe *fr,int flags)
{
    char pathname[GMX_PATH_MAX],filename[GMX_PATH_MAX];
    char *pathenv;
    const char *err;
    int i,ret=0;
    molfile_timestep_metadata_t *metadata=NULL;
#if !((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    glob_t globbuf;
    char *defpathenv = "/usr/local/lib/vmd/plugins/*/molfile";
#else
    WIN32_FIND_DATA ffd;
    HANDLE hFind = INVALID_HANDLE_VALUE;
    char progfolder[GMX_PATH_MAX];
    char defpathenv[GMX_PATH_MAX];
    SHGetFolderPath(NULL,CSIDL_PROGRAM_FILES,NULL,SHGFP_TYPE_CURRENT,progfolder);
    sprintf(defpathenv,"%s\\University of Illinois\\VMD\\plugins\\WIN32\\molfile",progfolder);
#endif
  

    fr->vmdplugin.api = NULL;
    fr->vmdplugin.filetype = strrchr(fn,'.')+1;

    pathenv = getenv("VMD_PLUGIN_PATH");
    if (pathenv==NULL) 
    {
        printf("\nVMD_PLUGIN_PATH not set. ");
        printf("Using default location:\n%s\n",defpathenv);
        pathenv=defpathenv;
    }
    strncpy(pathname,pathenv,sizeof(pathname));
#if !((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    strcat(pathname,"/*.so");
    glob(pathname, 0, NULL, &globbuf);
    if (globbuf.gl_pathc == 0)
    {
        printf("\nNo VMD Plugins found\n"
            "Set the environment variable VMD_PLUGIN_PATH to the molfile folder within the\n"
            "VMD installation.\n"
            "The architecture (e.g. 32bit versus 64bit) of Gromacs and VMD has to match.\n");
        return 0;
    }
    for (i=0; i<globbuf.gl_pathc && fr->vmdplugin.api == NULL; i++)
    {
        ret|=load_sharedlibrary_plugins(globbuf.gl_pathv[i],&fr->vmdplugin);
    }
    globfree(&globbuf);
#else
    strcat(pathname,"\\*.so");
    hFind = FindFirstFile(pathname, &ffd);
    if (INVALID_HANDLE_VALUE == hFind) 
    {
        printf("\nNo VMD Plugins found\n");
        return 0;
    } 
    do
    {
        sprintf(filename,"%s\\%s",pathenv,ffd.cFileName);
        ret|=load_sharedlibrary_plugins(filename,&fr->vmdplugin);
    }
    while (FindNextFile(hFind, &ffd )  != 0 && fr->vmdplugin.api == NULL );
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
            printf("Last error:\n%s\n",err);
        }
        return 0;
    }
    if (fr->vmdplugin.api == NULL)
    {
        printf("\nNo plugin for %s found\n",fr->vmdplugin.filetype);
        return 0;
    }

    fr->vmdplugin.handle = fr->vmdplugin.api->open_file_read(fn, fr->vmdplugin.filetype, &fr->natoms);

    if (!fr->vmdplugin.handle) {
        fprintf(stderr, "\nError: could not open file '%s' for reading.\n",
                fn);
        return 0;
    }

    if (fr->natoms < 1) {
        fprintf(stderr, "\nNo atoms found by VMD plugin in %s.\n"
            "Or format does not record number of atoms.\n", fn );
        return 0;
    }


    snew(fr->x,fr->natoms);

    fr->vmdplugin.bV = 0;
    if (fr->vmdplugin.api->read_timestep_metadata) 
    {
        fr->vmdplugin.api->read_timestep_metadata(fr->vmdplugin.handle, metadata);
        fr->vmdplugin.bV = metadata->has_velocities; 
        if (fr->vmdplugin.bV)
        {
            snew(fr->v,fr->natoms);
        }
    }
    return 1;

}



