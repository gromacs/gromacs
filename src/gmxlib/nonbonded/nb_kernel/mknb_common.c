/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#include <mknb_common.h>

/* This file just contains the global variables for make_nb_kernels,
 * as defined in mknb_common.h. Global variables are B*A*D, but all the 
 * function calls in the generator would be much more complicated 
 * without them, and that feels a bit unnecessary...
 *
 * It's only executed at compile-time anyway :-)
 */




/* Coulomb interaction alternatives */
const char *
mknb_coul_names[MKNB_COUL_NR] = {
	"Not calculated",
	"Normal Coulomb",
	"Reaction field",
	"Tabulated",
	"Generalized-Born"
};


/* VdW interaction alternatives */
const char *
mknb_vdw_names[MKNB_VDW_NR] = {
	"Not calculated",
	"Lennard-Jones",
	"Buckingham",
	"Tabulated"
};


/* Water optimization alternatives */
const char *
mknb_water_names[MKNB_WATER_NR] = {
	"No",
	"SPC/TIP3P - other atoms",
	"pairs of SPC/TIP3P interactions",
	"TIP4P - other atoms",
	"pairs of TIP4P interactions"
};



/* General program options, see mknb_common.h for definition of type. 
 */
struct mknb_options 
mknb_options;



/* Options for the kernel currently being generated. Definition in mknb_common.h.
 */
struct mknb_func
mknb_func;


