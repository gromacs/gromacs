/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_MDRUN_H
#define GMX_MDLIB_MDRUN_H

#include <bitset>

#include <stdio.h>
#include <time.h>

#include "gromacs/timing/wallcycle.h"

struct df_history_t;
struct gmx_constr;
struct gmx_edsam;
struct gmx_enerdata_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_expanded;
struct t_extmass;
struct t_inputrec;
struct t_lambda;
struct t_mdatoms;
struct t_simtemp;
class t_state;

enum MdFlags : size_t
{
    polarise             = 2,
    rerun                = 4,
    rerunVSite           = 5,
    ddBondCheck          = 10,
    ddBondComm           = 11,
    confOut              = 12,
    reproducible         = 13,
    appendFiles          = 15,
    appendFilesSet       = 21,
    keepAndNumCpt        = 16,
    readEKin             = 17,
    startFromCpt         = 18,
    resetCountersHalfWay = 19,
    tunePME              = 20,
    ntompSet             = 21,
    imdWait              = 23,
    imdTerm              = 24,
    imdPull              = 25
};
constexpr unsigned long long MD_POLARISE             = 1<< int(MdFlags::polarise);
constexpr unsigned long long MD_RERUN                = 1<< int(MdFlags::rerun);
constexpr unsigned long long MD_RERUN_VSITE          = 1<< int(MdFlags::rerunVSite);
constexpr unsigned long long MD_DDBONDCHECK          = 1<< int(MdFlags::ddBondCheck);
constexpr unsigned long long MD_DDBONDCOMM           = 1<< int(MdFlags::ddBondComm);
constexpr unsigned long long MD_CONFOUT              = 1<< int(MdFlags::confOut);
constexpr unsigned long long MD_REPRODUCIBLE         = 1<< int(MdFlags::reproducible);
constexpr unsigned long long MD_APPENDFILES          = 1<< int(MdFlags::appendFiles);
constexpr unsigned long long MD_APPENDFILESSET       = 1<< int(MdFlags::appendFilesSet);
constexpr unsigned long long MD_KEEPANDNUMCPT        = 1<< int(MdFlags::keepAndNumCpt);
constexpr unsigned long long MD_READ_EKIN            = 1<< int(MdFlags::readEKin);
constexpr unsigned long long MD_STARTFROMCPT         = 1<< int(MdFlags::startFromCpt);
constexpr unsigned long long MD_RESETCOUNTERSHALFWAY = 1<< int(MdFlags::resetCountersHalfWay);
constexpr unsigned long long MD_TUNEPME              = 1<< int(MdFlags::tunePME);
constexpr unsigned long long MD_NTOMPSET             = 1<< int(MdFlags::ntompSet);
constexpr unsigned long long MD_IMDWAIT              = 1<< int(MdFlags::imdWait);
constexpr unsigned long long MD_IMDTERM              = 1<< int(MdFlags::imdTerm);
constexpr unsigned long long MD_IMDPULL              = 1<< int(MdFlags::imdPull);

/* The options for the domain decomposition MPI task ordering */
enum {
    ddrankorderSEL, ddrankorderINTERLEAVE, ddrankorderPP_PME, ddrankorderCARTESIAN, ddrankorderNR
};

void init_npt_masses(t_inputrec *ir, t_state *state, t_extmass *MassQ, gmx_bool bInit);

void init_expanded_ensemble(gmx_bool bStateFromCP, t_inputrec *ir, df_history_t *dfhist);

int ExpandedEnsembleDynamics(FILE *log, t_inputrec *ir, gmx_enerdata_t *enerd,
                             t_state *state, t_extmass *MassQ, int fep_state, df_history_t *dfhist,
                             gmx_int64_t step,
                             rvec *v, t_mdatoms *mdatoms);

void PrintFreeEnergyInfoToFile(FILE *outfile, t_lambda *fep, t_expanded *expand, t_simtemp *simtemp, df_history_t *dfhist,
                               int fep_state, int frequency, gmx_int64_t step);

/* Allocate and initialize node-local state entries. */
void set_state_entries(t_state *state, const t_inputrec *ir);

/* Broadcast the data for a simulation, and allocate node-specific settings */
void init_parallel(t_commrec *cr, t_inputrec *inputrec,
                   gmx_mtop_t *mtop);

void bcast_state(const t_commrec *cr, t_state *state);
/* Broadcasts state from the master to all nodes in cr->mpi_comm_mygroup.
 */

#endif
