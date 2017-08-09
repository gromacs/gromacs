/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/genborn.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/timing/wallcycle.h"

struct gmx_device_info_t;
struct t_commrec;
struct t_fcdata;
struct t_filenm;

namespace gmx
{
class MDLogger;
}

/*! \brief Create a new forcerec structure */
t_forcerec *mk_forcerec(void);

/*! \brief Print the contents of the forcerec to a file
 *
 * \param[in] fplog The log file to print to
 * \param[in] fr    The forcerec structure
 */
void pr_forcerec(FILE *fplog, t_forcerec *fr);

/*! \brief Set the number of charge groups and atoms.
 *
 * The force calculation needs information on which atoms it
 * should do work.
 * \param[inout] fr                  The forcerec
 * \param[in]    ncg_home            Number of charge groups on this processor
 * \param[in]    ncg_force           Number of charge groups to compute force on
 * \param[in]    natoms_force        Number of atoms to compute force on
 * \param[in]    natoms_force_constr Number of atoms involved in constraints
 * \param[in]    natoms_f_novirsum   Number of atoms for which
 *                                   force is to be compute but no virial
 */
void
forcerec_set_ranges(t_forcerec *fr,
                    int ncg_home, int ncg_force,
                    int natoms_force,
                    int natoms_force_constr, int natoms_f_novirsum);

/*! \brief Initiate table constants
 *
 * Initializes the tables in the interaction constant data structure.
 * \param[in] fp   File for debugging output
 * \param[in] ic   Structure holding the table constant
 * \param[in] rtab The additional distance to add to tables
 */
void init_interaction_const_tables(FILE                   *fp,
                                   interaction_const_t    *ic,
                                   real                    rtab);

/*! \brief Initialize forcerec structure.
 *
 * The Force rec struct must be created with mk_forcerec.
 * \param[in]  fplog       File for printing
 * \param[in]  mdlog       File for printing
 * \param[out] fr          The forcerec
 * \param[in]  fcd         Force constant data
 * \param[in]  ir          Inputrec structure
 * \param[in]  mtop        Molecular topology
 * \param[in]  cr          Communication structures
 * \param[in]  box         Simulation box
 * \param[in]  tabfn       Table potential file for non-bonded interactions
 * \param[in]  tabpfn      Table potential file for pair interactions
 * \param[in]  tabbfnm     Table potential files for bonded interactions
 * \param[in]  deviceInfo  Info about GPU device to use for short-ranged work
 * \param[in]  bNoSolvOpt  Do not use solvent optimization
 * \param[in]  print_force Print forces for atoms with force >= print_force
 */
void init_forcerec(FILE                   *fplog,
                   const gmx::MDLogger    &mdlog,
                   t_forcerec             *fr,
                   t_fcdata               *fcd,
                   const t_inputrec       *ir,
                   const gmx_mtop_t       *mtop,
                   const t_commrec        *cr,
                   matrix                  box,
                   const char             *tabfn,
                   const char             *tabpfn,
                   const t_filenm         *tabbfnm,
                   gmx_device_info_t      *deviceInfo,
                   gmx_bool                bNoSolvOpt,
                   real                    print_force);

/*! \brief Divide exclusions over threads
 *
 * Set the exclusion load for the local exclusions and possibly threads
 * \param[out] fr  The force record
 * \param[in]  top The topology
 */
void forcerec_set_excl_load(t_forcerec           *fr,
                            const gmx_localtop_t *top);

/*! \brief Update parameters dependent on box
 *
 * Updates parameters in the forcerec that are time dependent
 * \param[out] fr  The force record
 * \param[in]  box The simulation box
 */
void update_forcerec(t_forcerec *fr, matrix box);

#endif
