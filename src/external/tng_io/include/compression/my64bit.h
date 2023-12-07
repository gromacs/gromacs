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

#ifndef MY64BIT_H
#define MY64BIT_H

#ifdef USE_STD_INTTYPES_H
#    include <inttypes.h>
typedef int64_t  my_int64_t;
typedef uint64_t my_uint64_t;
#    define HAVE64BIT
#else /* USE_STD_INTTYPES */
/* The USE_WINDOWS symbol should be automatically defined in tng_compress.h */
#    include "../compression/tng_compress.h"
#    ifdef USE_WINDOWS
typedef __int64          my_int64_t;
typedef unsigned __int64 my_uint64_t;
#        define HAVE64BIT
#    else /* USE_WINDOWS */
/* Fall back to assume that we have unsigned long long */
typedef unsigned long long my_uint64_t;
#        define HAVE64BIT
#    endif /* USE_WINDOWS */
#endif     /* USE_STD_INTTYPES */

#endif
