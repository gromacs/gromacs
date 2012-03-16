/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
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



/*! \file  nb_kernel_ia32_sse2_test_asm.h
 *  \brief Assembly routine to test ia32 SSE2 instructions
 *  \internal
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/*! \brief Try to execute a couple of ia32 SSE2 instructions
 *
 *  \internal
 *
 *  This routine does not produce any real result, but if
 *  ia32 SSE2 support is not present it will trigger an
 *  "illegal instruction" exception, which you should capture
 *  with a signal handling routine (unless you like crashes).
 *
 *  It is meant to be called from the higher level test
 *  routine nb_kernel_ia32_sse2_test().
 */
void
nb_kernel_ia32_sse2_test_asm(void);

#ifdef __cplusplus
}
#endif
