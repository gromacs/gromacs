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

/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
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

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vmddlopen.c,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.18 $      $Date: 2009/07/07 02:40:05 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Routines for loading dynamic link libraries and shared object files
 *   on various platforms, abstracting from machine dependent APIs.
 *
 ***************************************************************************/

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include "vmddlopen.h"

#if defined(__hpux)

#include <dl.h>
#include <errno.h>
#include <string.h>

void *vmddlopen( const char *path) {
    void *ret;
    ret = shl_load( path, BIND_IMMEDIATE | BIND_FIRST | BIND_VERBOSE, 0);
    return ret;
}

int vmddlclose( void *handle ) {
    return shl_unload( (shl_t) handle );
}

void *vmddlsym( void *handle, const char *sym ) {
    void *value=0;

    if ( shl_findsym( (shl_t*)&handle, sym, TYPE_UNDEFINED, &value ) != 0 ) 
	return 0;
    return value;
}

const char *vmddlerror( void  ) {
    return strerror( errno );
}

#elif 0 && defined(__APPLE__)
/*
 * This is only needed for MacOS X version 10.3 or older
 */
#include <mach-o/dyld.h>

void *vmddlopen( const char *path) {
  NSObjectFileImage image;
  NSObjectFileImageReturnCode retval;
  NSModule module;

  retval = NSCreateObjectFileImageFromFile(path, &image);
  if (retval != NSObjectFileImageSuccess)
    return NULL;

  module = NSLinkModule(image, path,
            NSLINKMODULE_OPTION_BINDNOW | NSLINKMODULE_OPTION_PRIVATE
            | NSLINKMODULE_OPTION_RETURN_ON_ERROR);
  return module;  /* module will be NULL on error */
}

int vmddlclose( void *handle ) {
  NSModule module = (NSModule *)handle;
  NSUnLinkModule(module, NSUNLINKMODULE_OPTION_NONE);
  return 0;
}

void *vmddlsym( void *handle, const char *symname ) {
  char *realsymname;
  NSModule module;
  NSSymbol sym;
  /* Hack around the leading underscore in the symbol name */
  realsymname = (char *)malloc(strlen(symname)+2);
  strcpy(realsymname, "_");
  strcat(realsymname, symname);
  module = (NSModule)handle;
  sym = NSLookupSymbolInModule(module, realsymname);
  free(realsymname);
  if (sym) 
    return (void *)(NSAddressOfSymbol(sym));
  return NULL;
}

const char *vmddlerror( void  ) {
  NSLinkEditErrors c;
  int errorNumber;
  const char *fileName;
  const char *errorString = NULL;
  NSLinkEditError(&c, &errorNumber, &fileName, &errorString);
  return errorString;
}

#elif defined( _WIN32 ) || defined( _WIN64 )

#include <windows.h>

void *vmddlopen(const char *fname) {
  return (void *)LoadLibrary(fname);
}

const char *vmddlerror(void) {
  static CHAR szBuf[80]; 
  DWORD dw = GetLastError(); 
 
  sprintf(szBuf, "vmddlopen failed: GetLastError returned %lu\n", dw);
  return szBuf;
}

void *vmddlsym(void *h, const char *sym) {
  return (void *)GetProcAddress((HINSTANCE)h, sym);
}

int vmddlclose(void *h) {
  /* FreeLibrary returns nonzero on success */
  return !FreeLibrary((HINSTANCE)h);
}

#else

/* All remaining platforms (not Windows, HP-UX, or MacOS X <= 10.3) */
#include <dlfcn.h>

void *vmddlopen(const char *fname) {
  return dlopen(fname, RTLD_NOW);
}
const char *vmddlerror(void) {
  return dlerror();
}
void *vmddlsym(void *h, const char *sym) {
  return dlsym(h, sym);
}
int vmddlclose(void *h) {
  return dlclose(h);
}
#endif 
