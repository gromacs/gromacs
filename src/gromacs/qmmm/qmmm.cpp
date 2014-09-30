/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "gromacs/fileio/confio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/ns.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/qmmmxx.h"
#include "gromacs/topology/invblock.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

// QMSystem will eventually take the place of QMrec.
gmx::QMSystem::QMSystem(int grpnr, int nr, int *atomarray, gmx_mtop_t *mtop, t_inputrec *ir) :
    nrQMatoms_(nr),
    xQM_(nr, std::vector<real>(3)),
    indexQM_(nr),
    atomicnumberQM_(nr),
    shiftQM_(nr),
    frontatoms_(nr)
{
    int                   i;
    gmx_mtop_atomlookup_t alook;
    t_atom               *atom;

    alook = gmx_mtop_atomlookup_init(mtop);

    for (i = 0; i < nrQMatoms_; i++)
    {
        gmx_mtop_atomnr_to_atom(alook, indexQM_[i], &atom);
        nelectrons_       += mtop->atomtypes.atomnumber[atom->type];
        atomicnumberQM_[i] = mtop->atomtypes.atomnumber[atom->type];
    }

    gmx_mtop_atomlookup_destroy(alook);

    QMcharge_       = ir->opts.QMcharge[grpnr];
    multiplicity_   = ir->opts.QMmult[grpnr];
    nelectrons_    -= ir->opts.QMcharge[grpnr];

    QMmethod_       = ir->opts.QMmethod[grpnr];
    QMbasis_        = ir->opts.QMbasis[grpnr];

    /* print the current layer to allow users to check their input */
    fprintf(stderr, "Layer %d\nnr of QM atoms %d\n", grpnr, nr);
    fprintf(stderr, "QMlevel: %s/%s\n\n",
            eQMmethod_names[QMmethod_], eQMbasis_names[QMbasis_]);

    bTS_      = ir->opts.bTS[grpnr];
    bOPT_     = ir->opts.bOPT[grpnr];
}

/*
   // HybridQuantumClassical will eventually take the place of QMMMrec.
   gmx::HybridQuantumClassical::HybridQuantumClassical(const t_commrec *cr, const gmx_mtop_t *mtop, const t_inputrec *ir, const t_forcerec *fr)
   {
   }

 */
