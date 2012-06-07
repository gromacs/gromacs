/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#ifndef GMX_IMD_IMD_H
#define GMX_IMD_IMD_H

/* How long shall we wait in seconds until we check for a connection again? */
#define IMDLOOPWAIT 1

/* How long shall we check for the IMD_GO? */
#define IMDCONNECTWAIT 2

/* IMD default port */
#define IMDDEFPORT 8888

#ifdef GMX_NATIVE_WINDOWS
#include <Windows.h>
#define NOFLAGS 0
#endif

/* Put this keyword before all GROMACS/IMD-related output: */
static const char IMDstr[] = "IMD:";

#include <limits.h>

/* GROMACS includes */
#include "typedefs.h"
#include "imdsocket.h"
#include "readinp.h"
#include "types/commrec.h"
#include "oenv.h"
#include "../fileio/filenm.h"

/* We define int32, which is the 32bit integer used in the IMD functions. */
#if ( INT_MAX == 2147483647 )
typedef int int32;
#else
typedef short int32;
#endif

/* IMD Protocol Version: */
#define HEADERSIZE 8
#define IMDVERSION 2

/* We use the same records as the NAMD/VMD IMD implementation. */
typedef enum IMDType_t
{
    IMD_DISCONNECT, /* client disconnect                                      */
    IMD_ENERGIES,   /* energy data                                            */
    IMD_FCOORDS,    /* atomic coordinates                                     */
    IMD_GO,         /* start command for the simulation                       */
    IMD_HANDSHAKE,  /* handshake to determine little/big endian               */
    IMD_KILL,       /* terminates the simulation                              */
    IMD_MDCOMM,     /* force data                                             */
    IMD_PAUSE,      /* pause the simulation                                   */
    IMD_TRATE,      /* sets the IMD transmission, and processing rate         */
    IMD_IOERROR,    /* I/O error                                              */
    IMD_NR
/* <-GROMACS specific extension to access message names */
} IMDMessageType;

/* Macros to access names for the IMDMessageType: */
#define UNDEFINED       "UNDEFINED"

/* Energy record as in the original IMD implementation, energies in kcal/mol. */
/* NOTE: We return the energies in GROMACS / SI units, so they also show up as SI in VMD. */
typedef struct
{
    int32 tstep;    /* time step                                              */
    float T_abs;    /* absolute temperature                                   */
    float E_tot;    /* total energy                                           */
    float E_pot;    /* potential energy                                       */
    float E_vdw;    /* van der Waals energy                                   */
    float E_coul;   /* Coulomb interaction energy                             */
    float E_bond;   /* bonds energy                                           */
    float E_angle;  /* angles energy                                          */
    float E_dihe;   /* dihedrals energy                                       */
    float E_impr;   /* improper dihedrals energy                              */
} IMDEnergyBlock;

/* Definitions for the IMD header & protocol version: */
typedef struct
{
    int32 type;
    int32 length;
} IMDHeader;


/* Public functions, see imd.c */

extern void write_imdatoms(t_inputrec *ir, t_state *state, gmx_mtop_t *sys, const char *fn);

extern void dd_make_local_IMD_atoms(gmx_domdec_t *dd, t_IMD *imd);

extern void init_imd(
        t_inputrec    *ir,
        t_commrec     *cr,
        gmx_mtop_t    *top_global,
        FILE          *fplog,
        int            defnstimd,
        int            nat_total,
        rvec           x[],
        t_mdatoms     *md,
        int            nfile,
        const t_filenm fnm[],
        output_env_t   oenv,
        int            imdport,
        int            imdfreq,
        unsigned long  Flags);

extern void do_imd_prepare_energies(t_gmx_IMD IMDSetup, gmx_enerdata_t *enerd, int step, gmx_bool bHaveNewEnergies);

/* Returns if we should do IMD communication in this step. Also checks for new IMD connection and syncs the nodes. */
extern gmx_bool do_IMD(
        int          step,
        t_commrec   *cr,
        gmx_bool     bNS,       /* Is this a ns step?          */
        matrix       box,
        rvec         x[],       /* The atomic positions local to this node    */
        t_inputrec  *ir,
        double       t);

extern int imd_get_step(t_gmx_IMD IMDSetup);

extern void do_imd_send_positions(
        t_IMD *imd);


extern void imd_apply_forces(t_inputrec *ir, t_commrec *cr, rvec *f);

extern void imd_finalize(t_IMD *imd);

#endif
