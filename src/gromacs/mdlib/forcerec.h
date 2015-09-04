/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_MDLIB_FORCEREC_H
#define GMX_MDLIB_FORCEREC_H

#include "gromacs/legacyheaders/genborn.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/tgroup.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/types/force_flags.h"
#include "gromacs/timing/wallcycle.h"

/*! Create a new forcerec structure */
t_forcerec *mk_forcerec(void);

/*! Print the contents of the forcerec to a file */
void pr_forcerec(FILE *fplog, t_forcerec *fr);

/*! Set the number of charge groups and atoms for the force calculation */
void
forcerec_set_ranges(t_forcerec *fr,
                    int ncg_home, int ncg_force,
                    int natoms_force,
                    int natoms_force_constr, int natoms_f_novirsum);

/*! Initializes the tables in the interaction constant data structure. */
void init_interaction_const_tables(FILE                   *fp,
                                   interaction_const_t    *ic,
                                   real                    rtab);

/*! \brief Initialize force rec structure.
 *
 * The Force rec struct must be created with mk_forcerec
 * The gmx_booleans have the following meaning:
 * bSetQ:    Copy the charges [ only necessary when they change ]
 * bMolEpot: Use the free energy stuff per molecule
 * print_force >= 0: print forces for atoms with force >= print_force
 */
void init_forcerec(FILE                   *fplog,
                   const output_env_t      oenv,
                   t_forcerec             *fr,
                   t_fcdata               *fcd,
                   const t_inputrec       *ir,
                   const gmx_mtop_t       *mtop,
                   const t_commrec        *cr,
                   matrix                  box,
                   const char             *tabfn,
                   const char             *tabafn,
                   const char             *tabpfn,
                   const char             *tabbfn,
                   const char             *nbpu_opt,
                   gmx_bool                bNoSolvOpt,
                   real                    print_force);

/*! Set the exclusion load for the local exclusions and possibly threads */
void forcerec_set_excl_load(t_forcerec           *fr,
                            const gmx_localtop_t *top);

/*! Updates parameters in the forcerec that are time dependent */
void update_forcerec(t_forcerec *fr, matrix box);

#endif
