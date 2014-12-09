/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#ifndef WARNMALLOC_H
#define WARNMALLOC_H

#include "../compression/tng_compress.h"

void DECLSPECDLLEXPORT *Ptngc_warnmalloc_x(const size_t size, char *file, const int line);

#define warnmalloc(size) Ptngc_warnmalloc_x(size,__FILE__,__LINE__)

void DECLSPECDLLEXPORT *Ptngc_warnrealloc_x(void *old, const size_t size, char *file, const int line);

#define warnrealloc(old,size) Ptngc_warnrealloc_x(old,size,__FILE__,__LINE__)


#endif
