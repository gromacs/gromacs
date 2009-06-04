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


