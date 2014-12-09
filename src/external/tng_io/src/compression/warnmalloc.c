/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#include <stdio.h>
#include <stdlib.h>
#include "../../include/compression/tng_compress.h"
#include "../../include/compression/warnmalloc.h"

void DECLSPECDLLEXPORT *Ptngc_warnmalloc_x(const size_t size, char *file, const int line)
{
  void *mem=malloc(size);
  if (!mem)
    {
      fprintf(stderr,"TRAJNG ERROR: Could not allocate memory of size %lu at %s:%d\n",(unsigned long) size,file,line);
      exit(EXIT_FAILURE);
    }
  return mem;
}

void DECLSPECDLLEXPORT *Ptngc_warnrealloc_x(void *old, const size_t size, char *file, const int line)
{
  void *mem=realloc(old,size);
  if (!mem)
    {
      fprintf(stderr,"TRAJNG ERROR: Could not allocate memory of size %lu at %s:%d\n",(unsigned long) size,file,line);
      exit(EXIT_FAILURE);
    }
  return mem;
}
