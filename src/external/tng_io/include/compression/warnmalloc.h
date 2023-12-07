/*
 * This code is part of the tng binary trajectory format.
 *
 * Copyright (c) 2010,2013, The GROMACS development team.
 * Copyright (c) 2020, by the GROMACS development team.
 * TNG was orginally written by Magnus Lundborg, Daniel Spångberg and
 * Rossen Apostolov. The API is implemented mainly by Magnus Lundborg,
 * Daniel Spångberg and Anders Gärdenäs.
 *
 * Please see the AUTHORS file for more information.
 *
 * The TNG library is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 *
 * To help us fund future development, we humbly ask that you cite
 * the research papers on the package.
 *
 * Check out http://www.gromacs.org for more information.
 */

/* This code is part of the tng compression routines
 * Written by Daniel Spangberg
 */

#ifndef WARNMALLOC_H
#define WARNMALLOC_H

#include "../compression/tng_compress.h"

void DECLSPECDLLEXPORT* Ptngc_warnmalloc_x(size_t size, char* file, int line);

#define warnmalloc(size) Ptngc_warnmalloc_x(size, __FILE__, __LINE__)

void DECLSPECDLLEXPORT* Ptngc_warnrealloc_x(void* old, size_t size, char* file, int line);

#define warnrealloc(old, size) Ptngc_warnrealloc_x(old, size, __FILE__, __LINE__)


#endif
