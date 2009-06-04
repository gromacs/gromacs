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
#ifndef _MKNB_DECLARATIONS_H_
#define _MKNB_DECLARATIONS_H_


/*! \file   mknb_declarations.h
 *  \brief Kernel generator (only for compile): Function declarations
 * 
 *  \internal
 *
 *  \note This file is only used to generate the inner loop kernels
 *        at compile time, which in turn are included in Gromacs. 
 *        This code itself is NOT linked into the Gromacs library, so
 *        it does not need to be threadsafe.
 *
 *  This file is only used to generate the inner loop kernels
 *  at compile time, which in turn are included in Gromacs.
 *  This code itself is NOT linked into the Gromacs library, so
 *  it does not need to be threadsafe.
 *
 *  mknb_declarations.h provides routines to generate function
 *  headers for nonbonded interaction kernels in either C or
 *  Fortran, routines to declare local variables and local data
 *  initialization before the actual loop starts.
 *
 *  Definitions and structures used in the nonbonded kernel 
 *  generator itself (i.e. not in the generated code) are in
 *  the file mknb_common.h.
 */



/*! \brief Kernel generator (only for compile): Function definition and header
 *
 *  \internal
 *
 *  \note   Only defined/used in the nonbonded kernel generator 
 *          program mknb. This program is run once at compile
 *          time to create the inner loops, and then discarded.
 *          This source is NOT linked into any Gromacs library.
 *
 *  All level1/level2 Gromacs nonbonded kernels use the same
 *  parameters, so the only difference here is the function 
 *  name and the language used.
 *
 *  \param funcname   Name of the nonbonded kernel to write
 */
void
mknb_function_header       (char *     funcname);



/*! \brief Kernel generator (only for compile): Local function variables 
 *
 *  \internal
 *
 *  \note   Only defined/used in the nonbonded kernel generator 
 *          program mknb. This program is run once at compile
 *          time to create the inner loops, and then discarded.
 *          This source is NOT linked into any Gromacs library.
 *
 *  This routine tries to determine which variables will
 *  be needed based on the current values in the global
 *  variable 'func', which containts the current settings
 *  for Coulomb/L-J interactions, water optimization, etc.
 */
void
mknb_declare_variables     (void);



/*! \brief Kernel generator (only for compile): Init code before loop
 *
 *  \internal
 *
 *  \note   Only defined/used in the nonbonded kernel generator 
 *          program mknb. This program is run once at compile
 *          time to create the inner loops, and then discarded.
 *          This source is NOT linked into any Gromacs library.
 *
 *  Some nonbonded kernels (mainly water optimized ones)
 *  extracts parameters or similar things outside the
 *  outer loop for better performance. This code should
 *  be run after mknb_declare_variables(), but before
 *  the outer loop is started.
 */
void
mknb_initialize_data       (void);



/*! \brief Kernel generator (only for compile): Close a function body.
 *
 *  \internal
 *
 *  \note   Only defined/used in the nonbonded kernel generator 
 *          program mknb. This program is run once at compile
 *          time to create the inner loops, and then discarded.
 *          This source is NOT linked into any Gromacs library.
 *
 *  This finishes a function previously started with
 *  mknb_function_header().
 *
 */
void
mknb_finish_function       (void);




#endif /* _MKNB_DECLARATIONS_H_ */
