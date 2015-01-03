/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * This file written by Justin A. Lemkul and possibly others.
 * Copyright (c) 2013,2014 Justin A. Lemkul 
 * All rights reserved.

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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

#include "readir.h"
#include "names.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"

void read_drude_opts(int *ninp_p, t_inpfile **inp_p, t_drude *drude, warninp_t wi)
{
    const char *tmp;
    int        ninp;
    t_inpfile *inp;

    ninp   = *ninp_p;
    inp    = *inp_p;

    EETYPE("drude-mode", drude->drudemode, edrude_modes);
    RTYPE ("drude-t", drude->drude_t, 1.0);
    EETYPE("drude-hardwall", drude->bHardWall, yesno_names);
    RTYPE ("drude-r", drude->drude_r, 0.02);
    EETYPE("drude-hyper", drude->bHyper, yesno_names);
    RTYPE ("drude-khyp", drude->drude_khyp, 16736000.00);   /* default CHARMM value: 40,000 kcal/mol/A^2 */
    ITYPE ("drude-pow", drude->drude_hyp_power, 4);
    RTYPE ("nbtholecut", drude->nbtholecut, 0.0);
    ITYPE ("drude-tsteps", drude->tsteps, 20);

    *ninp_p   = ninp;
    *inp_p    = inp;
}
