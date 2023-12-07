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

#include <stdio.h>
#include <stdlib.h>
#include "../../include/compression/tng_compress.h"
#include "../../include/compression/warnmalloc.h"

void DECLSPECDLLEXPORT* Ptngc_warnmalloc_x(const size_t size, char* file, const int line)
{
    void* mem = malloc(size);
    if (!mem)
    {
        fprintf(stderr, "TRAJNG ERROR: Could not allocate memory of size %lu at %s:%d\n",
                (unsigned long)size, file, line);
        exit(EXIT_FAILURE);
    }
    return mem;
}

void DECLSPECDLLEXPORT* Ptngc_warnrealloc_x(void* old, const size_t size, char* file, const int line)
{
    void* mem = realloc(old, size);
    if (!mem)
    {
        fprintf(stderr, "TRAJNG ERROR: Could not allocate memory of size %lu at %s:%d\n",
                (unsigned long)size, file, line);
        exit(EXIT_FAILURE);
    }
    return mem;
}
