/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
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
#include "gmxpre.h"

#include "shellfc.h"

#include <stdlib.h>
#include <string.h>

#include <cstdint>

#include <algorithm>
#include <array>

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/chargegroup.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    int     nnucl;
    int     shell;               /* The shell id				*/
    int     nucl1, nucl2, nucl3; /* The nuclei connected to the shell	*/
    /* gmx_bool    bInterCG; */       /* Coupled to nuclei outside cg?        */
    real    k;                   /* force constant		        */
    real    k_1;                 /* 1 over force constant		*/
    rvec    xold;
    rvec    fold;
    rvec    step;
} t_shell;

struct gmx_shellfc_t {
    /* Shell counts, indices, parameters and working data */
    int          nshell_gl;              /* The number of shells in the system        */
    t_shell     *shell_gl;               /* All the shells (for DD only)              */
    int         *shell_index_gl;         /* Global shell index (for DD only)          */
    gmx_bool     bInterCG;               /* Are there inter charge-group shells?      */
    int          nshell;                 /* The number of local shells                */
    t_shell     *shell;                  /* The local shells                          */
    int          shell_nalloc;           /* The allocation size of shell              */
    gmx_bool     bPredict;               /* Predict shell positions                   */
    gmx_bool     bRequireInit;           /* Require initialization of shell positions */
    int          nflexcon;               /* The number of flexible constraints        */

    /* Temporary arrays, should be fixed size 2 when fully converted to C++ */
    PaddedRVecVector *x;                 /* Array for iterative minimization          */
    PaddedRVecVector *f;                 /* Array for iterative minimization          */

    /* Flexible constraint working data */
    rvec        *acc_dir;                /* Acceleration direction for flexcon        */
    rvec        *x_old;                  /* Old coordinates for flexcon               */
    int          flex_nalloc;            /* The allocation size of acc_dir and x_old  */
    rvec        *adir_xnold;             /* Work space for init_adir                  */
    rvec        *adir_xnew;              /* Work space for init_adir                  */
    int          adir_nalloc;            /* Work space for init_adir                  */
    std::int64_t numForceEvaluations;    /* Total number of force evaluations         */
    int          numConvergedIterations; /* Total number of iterations that converged */
};


static void pr_shell(FILE *fplog, int ns, t_shell s[])
{
    int i;

    fprintf(fplog, "SHELL DATA\n");
    fprintf(fplog, "%5s  %8s  %5s  %5s  %5s\n",
            "Shell", "Force k", "Nucl1", "Nucl2", "Nucl3");
    for (i = 0; (i < ns); i++)
    {
        fprintf(fplog, "%5d  %8.3f  %5d", s[i].shell, 1.0/s[i].k_1, s[i].nucl1);
        if (s[i].nnucl == 2)
        {
            fprintf(fplog, "  %5d\n", s[i].nucl2);
        }
        else if (s[i].nnucl == 3)
        {
            fprintf(fplog, "  %5d  %5d\n", s[i].nucl2, s[i].nucl3);
        }
        else
        {
            fprintf(fplog, "\n");
        }
    }
}

/* TODO The remain call of this function passes non-NULL mass and NULL
 * mtop, so this routine can be simplified.
 *
 * The other code path supported doing prediction before the MD loop
 * started, but even when called, the prediction was always
 * over-written by a subsequent call in the MD loop, so has been
 * removed. */
static void predict_shells(FILE *fplog, rvec x[], rvec v[], real dt,
                           int ns, t_shell s[],
                           real mass[], gmx_mtop_t *mtop, gmx_bool bInit)
{
    int                   i, m, s1, n1, n2, n3;
    real                  dt_1, fudge, tm, m1, m2, m3;
    rvec                 *ptr;

    /* We introduce a fudge factor for performance reasons: with this choice
     * the initial force on the shells is about a factor of two lower than
     * without
     */
    fudge = 1.0;

    if (bInit)
    {
        if (fplog)
        {
            fprintf(fplog, "RELAX: Using prediction for initial shell placement\n");
        }
        ptr  = x;
        dt_1 = 1;
    }
    else
    {
        ptr  = v;
        dt_1 = fudge*dt;
    }

    int molb = 0;
    for (i = 0; (i < ns); i++)
    {
        s1 = s[i].shell;
        if (bInit)
        {
            clear_rvec(x[s1]);
        }
        switch (s[i].nnucl)
        {
            case 1:
                n1 = s[i].nucl1;
                for (m = 0; (m < DIM); m++)
                {
                    x[s1][m] += ptr[n1][m]*dt_1;
                }
                break;
            case 2:
                n1 = s[i].nucl1;
                n2 = s[i].nucl2;
                if (mass)
                {
                    m1 = mass[n1];
                    m2 = mass[n2];
                }
                else
                {
                    /* Not the correct masses with FE, but it is just a prediction... */
                    m1 = mtopGetAtomMass(mtop, n1, &molb);
                    m2 = mtopGetAtomMass(mtop, n2, &molb);
                }
                tm = dt_1/(m1+m2);
                for (m = 0; (m < DIM); m++)
                {
                    x[s1][m] += (m1*ptr[n1][m]+m2*ptr[n2][m])*tm;
                }
                break;
            case 3:
                n1 = s[i].nucl1;
                n2 = s[i].nucl2;
                n3 = s[i].nucl3;
                if (mass)
                {
                    m1 = mass[n1];
                    m2 = mass[n2];
                    m3 = mass[n3];
                }
                else
                {
                    /* Not the correct masses with FE, but it is just a prediction... */
                    m1 = mtopGetAtomMass(mtop, n1, &molb);
                    m2 = mtopGetAtomMass(mtop, n2, &molb);
                    m3 = mtopGetAtomMass(mtop, n3, &molb);
                }
                tm = dt_1/(m1+m2+m3);
                for (m = 0; (m < DIM); m++)
                {
                    x[s1][m] += (m1*ptr[n1][m]+m2*ptr[n2][m]+m3*ptr[n3][m])*tm;
                }
                break;
            default:
                gmx_fatal(FARGS, "Shell %d has %d nuclei!", i, s[i].nnucl);
        }
    }
}

/*! \brief Count the different particle types in a system
 *
 * Routine prints a warning to stderr in case an unknown particle type
 * is encountered.
 * \param[in]  fplog Print what we have found if not NULL
 * \param[in]  mtop  Molecular topology.
 * \returns Array holding the number of particles of a type
 */
static std::array<int, eptNR> countPtypes(FILE       *fplog,
                                          gmx_mtop_t *mtop)
{
    std::array<int, eptNR> nptype = { { 0 } };
    /* Count number of shells, and find their indices */
    for (int i = 0; (i < eptNR); i++)
    {
        nptype[i] = 0;
    }

    gmx_mtop_atomloop_block_t  aloopb = gmx_mtop_atomloop_block_init(mtop);
    int                        nmol;
    const t_atom              *atom;
    while (gmx_mtop_atomloop_block_next(aloopb, &atom, &nmol))
    {
        switch (atom->ptype)
        {
            case eptAtom:
            case eptVSite:
            case eptShell:
                nptype[atom->ptype] += nmol;
                break;
            default:
                fprintf(stderr, "Warning unsupported particle type %d in countPtypes",
                        static_cast<int>(atom->ptype));
        }
    }
    if (fplog)
    {
        /* Print the number of each particle type */
        int n = 0;
        for (const auto &i : nptype)
        {
            if (i != 0)
            {
                fprintf(fplog, "There are: %d %ss\n", i, ptype_str[n]);
            }
            n++;
        }
    }
    return nptype;
}

gmx_shellfc_t *init_shell_flexcon(FILE *fplog,
                                  gmx_mtop_t *mtop, int nflexcon,
                                  int nstcalcenergy,
                                  bool usingDomainDecomposition)
{
    gmx_shellfc_t            *shfc;
    t_shell                  *shell;
    int                      *shell_index = nullptr, *at2cg;
    const t_atom             *atom;

    int                       ns, nshell, nsi;
    int                       i, j, type, mb, a_offset, cg, mol, ftype, nra;
    real                      qS, alpha;
    int                       aS, aN = 0; /* Shell and nucleus */
    int                       bondtypes[] = { F_BONDS, F_HARMONIC, F_CUBICBONDS, F_POLARIZATION, F_ANHARM_POL, F_WATER_POL };
#define NBT asize(bondtypes)
    t_iatom                  *ia;
    gmx_mtop_atomloop_all_t   aloop;
    gmx_ffparams_t           *ffparams;
    gmx_molblock_t           *molb;
    gmx_moltype_t            *molt;
    t_block                  *cgs;

    std::array<int, eptNR>    n = countPtypes(fplog, mtop);
    nshell = n[eptShell];

    if (nshell == 0 && nflexcon == 0)
    {
        /* We're not doing shells or flexible constraints */
        return nullptr;
    }

    snew(shfc, 1);
    shfc->x        = new PaddedRVecVector[2] {};
    shfc->f        = new PaddedRVecVector[2] {};
    shfc->nflexcon = nflexcon;

    if (nshell == 0)
    {
        /* Only flexible constraints, no shells.
         * Note that make_local_shells() does not need to be called.
         */
        shfc->nshell   = 0;
        shfc->bPredict = FALSE;

        return shfc;
    }

    if (nstcalcenergy != 1)
    {
        gmx_fatal(FARGS, "You have nstcalcenergy set to a value (%d) that is different from 1.\nThis is not supported in combination with shell particles.\nPlease make a new tpr file.", nstcalcenergy);
    }
    if (usingDomainDecomposition)
    {
        gmx_fatal(FARGS, "Shell particles are not implemented with domain decomposition, use a single rank");
    }

    /* We have shells: fill the shell data structure */

    /* Global system sized array, this should be avoided */
    snew(shell_index, mtop->natoms);

    aloop  = gmx_mtop_atomloop_all_init(mtop);
    nshell = 0;
    while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
    {
        if (atom->ptype == eptShell)
        {
            shell_index[i] = nshell++;
        }
    }

    snew(shell, nshell);

    /* Initiate the shell structures */
    for (i = 0; (i < nshell); i++)
    {
        shell[i].shell = -1;
        shell[i].nnucl = 0;
        shell[i].nucl1 = -1;
        shell[i].nucl2 = -1;
        shell[i].nucl3 = -1;
        /* shell[i].bInterCG=FALSE; */
        shell[i].k_1   = 0;
        shell[i].k     = 0;
    }

    ffparams = &mtop->ffparams;

    /* Now fill the structures */
    shfc->bInterCG = FALSE;
    ns             = 0;
    a_offset       = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];

        cgs = &molt->cgs;
        snew(at2cg, molt->atoms.nr);
        for (cg = 0; cg < cgs->nr; cg++)
        {
            for (i = cgs->index[cg]; i < cgs->index[cg+1]; i++)
            {
                at2cg[i] = cg;
            }
        }

        atom = molt->atoms.atom;
        for (mol = 0; mol < molb->nmol; mol++)
        {
            for (j = 0; (j < NBT); j++)
            {
                ia = molt->ilist[bondtypes[j]].iatoms;
                for (i = 0; (i < molt->ilist[bondtypes[j]].nr); )
                {
                    type  = ia[0];
                    ftype = ffparams->functype[type];
                    nra   = interaction_function[ftype].nratoms;

                    /* Check whether we have a bond with a shell */
                    aS = -1;

                    switch (bondtypes[j])
                    {
                        case F_BONDS:
                        case F_HARMONIC:
                        case F_CUBICBONDS:
                        case F_POLARIZATION:
                        case F_ANHARM_POL:
                            if (atom[ia[1]].ptype == eptShell)
                            {
                                aS = ia[1];
                                aN = ia[2];
                            }
                            else if (atom[ia[2]].ptype == eptShell)
                            {
                                aS = ia[2];
                                aN = ia[1];
                            }
                            break;
                        case F_WATER_POL:
                            aN    = ia[4]; /* Dummy */
                            aS    = ia[5]; /* Shell */
                            break;
                        default:
                            gmx_fatal(FARGS, "Death Horror: %s, %d", __FILE__, __LINE__);
                    }

                    if (aS != -1)
                    {
                        qS = atom[aS].q;

                        /* Check whether one of the particles is a shell... */
                        nsi = shell_index[a_offset+aS];
                        if ((nsi < 0) || (nsi >= nshell))
                        {
                            gmx_fatal(FARGS, "nsi is %d should be within 0 - %d. aS = %d",
                                      nsi, nshell, aS);
                        }
                        if (shell[nsi].shell == -1)
                        {
                            shell[nsi].shell = a_offset + aS;
                            ns++;
                        }
                        else if (shell[nsi].shell != a_offset+aS)
                        {
                            gmx_fatal(FARGS, "Weird stuff in %s, %d", __FILE__, __LINE__);
                        }

                        if      (shell[nsi].nucl1 == -1)
                        {
                            shell[nsi].nucl1 = a_offset + aN;
                        }
                        else if (shell[nsi].nucl2 == -1)
                        {
                            shell[nsi].nucl2 = a_offset + aN;
                        }
                        else if (shell[nsi].nucl3 == -1)
                        {
                            shell[nsi].nucl3 = a_offset + aN;
                        }
                        else
                        {
                            if (fplog)
                            {
                                pr_shell(fplog, ns, shell);
                            }
                            gmx_fatal(FARGS, "Can not handle more than three bonds per shell\n");
                        }
                        if (at2cg[aS] != at2cg[aN])
                        {
                            /* shell[nsi].bInterCG = TRUE; */
                            shfc->bInterCG = TRUE;
                        }

                        switch (bondtypes[j])
                        {
                            case F_BONDS:
                            case F_HARMONIC:
                                shell[nsi].k    += ffparams->iparams[type].harmonic.krA;
                                break;
                            case F_CUBICBONDS:
                                shell[nsi].k    += ffparams->iparams[type].cubic.kb;
                                break;
                            case F_POLARIZATION:
                            case F_ANHARM_POL:
                                if (!gmx_within_tol(qS, atom[aS].qB, GMX_REAL_EPS*10))
                                {
                                    gmx_fatal(FARGS, "polarize can not be used with qA(%e) != qB(%e) for atom %d of molecule block %d", qS, atom[aS].qB, aS+1, mb+1);
                                }
                                shell[nsi].k    += gmx::square(qS)*ONE_4PI_EPS0/
                                    ffparams->iparams[type].polarize.alpha;
                                break;
                            case F_WATER_POL:
                                if (!gmx_within_tol(qS, atom[aS].qB, GMX_REAL_EPS*10))
                                {
                                    gmx_fatal(FARGS, "water_pol can not be used with qA(%e) != qB(%e) for atom %d of molecule block %d", qS, atom[aS].qB, aS+1, mb+1);
                                }
                                alpha          = (ffparams->iparams[type].wpol.al_x+
                                                  ffparams->iparams[type].wpol.al_y+
                                                  ffparams->iparams[type].wpol.al_z)/3.0;
                                shell[nsi].k  += gmx::square(qS)*ONE_4PI_EPS0/alpha;
                                break;
                            default:
                                gmx_fatal(FARGS, "Death Horror: %s, %d", __FILE__, __LINE__);
                        }
                        shell[nsi].nnucl++;
                    }
                    ia += nra+1;
                    i  += nra+1;
                }
            }
            a_offset += molt->atoms.nr;
        }
        /* Done with this molecule type */
        sfree(at2cg);
    }

    /* Verify whether it's all correct */
    if (ns != nshell)
    {
        gmx_fatal(FARGS, "Something weird with shells. They may not be bonded to something");
    }

    for (i = 0; (i < ns); i++)
    {
        shell[i].k_1 = 1.0/shell[i].k;
    }

    if (debug)
    {
        pr_shell(debug, ns, shell);
    }


    shfc->nshell_gl      = ns;
    shfc->shell_gl       = shell;
    shfc->shell_index_gl = shell_index;

    shfc->bPredict     = (getenv("GMX_NOPREDICT") == nullptr);
    shfc->bRequireInit = FALSE;
    if (!shfc->bPredict)
    {
        if (fplog)
        {
            fprintf(fplog, "\nWill never predict shell positions\n");
        }
    }
    else
    {
        shfc->bRequireInit = (getenv("GMX_REQUIRE_SHELL_INIT") != nullptr);
        if (shfc->bRequireInit && fplog)
        {
            fprintf(fplog, "\nWill always initiate shell positions\n");
        }
    }

    if (shfc->bPredict)
    {
        if (shfc->bInterCG)
        {
            if (fplog)
            {
                fprintf(fplog, "\nNOTE: there all shells that are connected to particles outside thier own charge group, will not predict shells positions during the run\n\n");
            }
            /* Prediction improves performance, so we should implement either:
             * 1. communication for the atoms needed for prediction
             * 2. prediction using the velocities of shells; currently the
             *    shell velocities are zeroed, it's a bit tricky to keep
             *    track of the shell displacements and thus the velocity.
             */
            shfc->bPredict = FALSE;
        }
    }

    return shfc;
}

void make_local_shells(t_commrec *cr, t_mdatoms *md,
                       gmx_shellfc_t *shfc)
{
    t_shell      *shell;
    int           a0, a1, *ind, nshell, i;
    gmx_domdec_t *dd = nullptr;

    if (DOMAINDECOMP(cr))
    {
        dd = cr->dd;
        a0 = 0;
        a1 = dd->nat_home;
    }
    else
    {
        /* Single node: we need all shells, just copy the pointer */
        shfc->nshell = shfc->nshell_gl;
        shfc->shell  = shfc->shell_gl;

        return;
    }

    ind = shfc->shell_index_gl;

    nshell = 0;
    shell  = shfc->shell;
    for (i = a0; i < a1; i++)
    {
        if (md->ptype[i] == eptShell)
        {
            if (nshell+1 > shfc->shell_nalloc)
            {
                shfc->shell_nalloc = over_alloc_dd(nshell+1);
                srenew(shell, shfc->shell_nalloc);
            }
            if (dd)
            {
                shell[nshell] = shfc->shell_gl[ind[dd->gatindex[i]]];
            }
            else
            {
                shell[nshell] = shfc->shell_gl[ind[i]];
            }

            /* With inter-cg shells we can no do shell prediction,
             * so we do not need the nuclei numbers.
             */
            if (!shfc->bInterCG)
            {
                shell[nshell].nucl1   = i + shell[nshell].nucl1 - shell[nshell].shell;
                if (shell[nshell].nnucl > 1)
                {
                    shell[nshell].nucl2 = i + shell[nshell].nucl2 - shell[nshell].shell;
                }
                if (shell[nshell].nnucl > 2)
                {
                    shell[nshell].nucl3 = i + shell[nshell].nucl3 - shell[nshell].shell;
                }
            }
            shell[nshell].shell = i;
            nshell++;
        }
    }

    shfc->nshell = nshell;
    shfc->shell  = shell;
}

static void do_1pos(rvec xnew, const rvec xold, const rvec f, real step)
{
    real xo, yo, zo;
    real dx, dy, dz;

    xo = xold[XX];
    yo = xold[YY];
    zo = xold[ZZ];

    dx = f[XX]*step;
    dy = f[YY]*step;
    dz = f[ZZ]*step;

    xnew[XX] = xo+dx;
    xnew[YY] = yo+dy;
    xnew[ZZ] = zo+dz;
}

static void do_1pos3(rvec xnew, const rvec xold, const rvec f, const rvec step)
{
    real xo, yo, zo;
    real dx, dy, dz;

    xo = xold[XX];
    yo = xold[YY];
    zo = xold[ZZ];

    dx = f[XX]*step[XX];
    dy = f[YY]*step[YY];
    dz = f[ZZ]*step[ZZ];

    xnew[XX] = xo+dx;
    xnew[YY] = yo+dy;
    xnew[ZZ] = zo+dz;
}

static void directional_sd(gmx::ArrayRef<const gmx::RVec> xold,
                           gmx::ArrayRef<gmx::RVec> xnew,
                           const rvec acc_dir[], int homenr, real step)
{
    const rvec *xo = as_rvec_array(xold.data());
    rvec       *xn = as_rvec_array(xnew.data());

    for (int i = 0; i < homenr; i++)
    {
        do_1pos(xn[i], xo[i], acc_dir[i], step);
    }
}

static void shell_pos_sd(gmx::ArrayRef<const gmx::RVec> xcur,
                         gmx::ArrayRef<gmx::RVec> xnew,
                         gmx::ArrayRef<const gmx::RVec> f,
                         int ns, t_shell s[], int count)
{
    const real step_scale_min       = 0.8,
               step_scale_increment = 0.2,
               step_scale_max       = 1.2,
               step_scale_multiple  = (step_scale_max - step_scale_min) / step_scale_increment;
    int        i, shell, d;
    real       dx, df, k_est;
    const real zero = 0;
#ifdef PRINT_STEP
    real       step_min, step_max;

    step_min = 1e30;
    step_max = 0;
#endif
    for (i = 0; (i < ns); i++)
    {
        shell = s[i].shell;
        if (count == 1)
        {
            for (d = 0; d < DIM; d++)
            {
                s[i].step[d] = s[i].k_1;
#ifdef PRINT_STEP
                step_min = std::min(step_min, s[i].step[d]);
                step_max = std::max(step_max, s[i].step[d]);
#endif
            }
        }
        else
        {
            for (d = 0; d < DIM; d++)
            {
                dx = xcur[shell][d] - s[i].xold[d];
                df =    f[shell][d] - s[i].fold[d];
                /* -dx/df gets used to generate an interpolated value, but would
                 * cause a NaN if df were binary-equal to zero. Values close to
                 * zero won't cause problems (because of the min() and max()), so
                 * just testing for binary inequality is OK. */
                if (zero != df)
                {
                    k_est = -dx/df;
                    /* Scale the step size by a factor interpolated from
                     * step_scale_min to step_scale_max, as k_est goes from 0 to
                     * step_scale_multiple * s[i].step[d] */
                    s[i].step[d] =
                        step_scale_min * s[i].step[d] +
                        step_scale_increment * std::min(step_scale_multiple * s[i].step[d], std::max(k_est, zero));
                }
                else
                {
                    /* Here 0 == df */
                    if (gmx_numzero(dx)) /* 0 == dx */
                    {
                        /* Likely this will never happen, but if it does just
                         * don't scale the step. */
                    }
                    else /* 0 != dx */
                    {
                        s[i].step[d] *= step_scale_max;
                    }
                }
#ifdef PRINT_STEP
                step_min = std::min(step_min, s[i].step[d]);
                step_max = std::max(step_max, s[i].step[d]);
#endif
            }
        }
        copy_rvec(xcur [shell], s[i].xold);
        copy_rvec(f[shell],   s[i].fold);

        do_1pos3(xnew[shell], xcur[shell], f[shell], s[i].step);

        if (gmx_debug_at)
        {
            fprintf(debug, "shell[%d] = %d\n", i, shell);
            pr_rvec(debug, 0, "fshell", f[shell], DIM, TRUE);
            pr_rvec(debug, 0, "xold", xcur[shell], DIM, TRUE);
            pr_rvec(debug, 0, "step", s[i].step, DIM, TRUE);
            pr_rvec(debug, 0, "xnew", xnew[shell], DIM, TRUE);
        }
    }
#ifdef PRINT_STEP
    printf("step %.3e %.3e\n", step_min, step_max);
#endif
}

static void decrease_step_size(int nshell, t_shell s[])
{
    int i;

    for (i = 0; i < nshell; i++)
    {
        svmul(0.8, s[i].step, s[i].step);
    }
}

static void print_epot(FILE *fp, gmx_int64_t mdstep, int count, real epot, real df,
                       int ndir, real sf_dir)
{
    char buf[22];

    fprintf(fp, "MDStep=%5s/%2d EPot: %12.8e, rmsF: %6.2e",
            gmx_step_str(mdstep, buf), count, epot, df);
    if (ndir)
    {
        fprintf(fp, ", dir. rmsF: %6.2e\n", std::sqrt(sf_dir/ndir));
    }
    else
    {
        fprintf(fp, "\n");
    }
}


static real rms_force(t_commrec *cr, gmx::ArrayRef<const gmx::RVec> force, int ns, t_shell s[],
                      int ndir, real *sf_dir, real *Epot)
{
    double      buf[4];
    const rvec *f = as_rvec_array(force.data());

    buf[0] = *sf_dir;
    for (int i = 0; i < ns; i++)
    {
        int shell  = s[i].shell;
        buf[0]    += norm2(f[shell]);
    }
    int ntot = ns;

    if (PAR(cr))
    {
        buf[1] = ntot;
        buf[2] = *sf_dir;
        buf[3] = *Epot;
        gmx_sumd(4, buf, cr);
        ntot    = (int)(buf[1] + 0.5);
        *sf_dir = buf[2];
        *Epot   = buf[3];
    }
    ntot += ndir;

    return (ntot ? std::sqrt(buf[0]/ntot) : 0);
}

static void check_pbc(FILE *fp, gmx::ArrayRef<gmx::RVec> x, int shell)
{
    int m, now;

    now = shell-4;
    for (m = 0; (m < DIM); m++)
    {
        if (fabs(x[shell][m]-x[now][m]) > 0.3)
        {
            pr_rvecs(fp, 0, "SHELL-X", as_rvec_array(x.data())+now, 5);
            break;
        }
    }
}

static void dump_shells(FILE *fp, gmx::ArrayRef<gmx::RVec> x, gmx::ArrayRef<gmx::RVec> f, real ftol, int ns, t_shell s[])
{
    int  i, shell;
    real ft2, ff2;

    ft2 = gmx::square(ftol);

    for (i = 0; (i < ns); i++)
    {
        shell = s[i].shell;
        ff2   = iprod(f[shell], f[shell]);
        if (ff2 > ft2)
        {
            fprintf(fp, "SHELL %5d, force %10.5f  %10.5f  %10.5f, |f| %10.5f\n",
                    shell, f[shell][XX], f[shell][YY], f[shell][ZZ], std::sqrt(ff2));
        }
        check_pbc(fp, x, shell);
    }
}

static void init_adir(FILE *log, gmx_shellfc_t *shfc,
                      gmx_constr_t constr, t_idef *idef, t_inputrec *ir,
                      t_commrec *cr, int dd_ac1,
                      gmx_int64_t step, t_mdatoms *md, int end,
                      rvec *x_old, rvec *x_init, rvec *x,
                      rvec *f, rvec *acc_dir,
                      gmx_bool bMolPBC, matrix box,
                      gmx::ArrayRef<const real> lambda, real *dvdlambda,
                      t_nrnb *nrnb)
{
    rvec           *xnold, *xnew;
    double          dt, w_dt;
    int             n, d;
    unsigned short *ptype;

    if (DOMAINDECOMP(cr))
    {
        n = dd_ac1;
    }
    else
    {
        n = end;
    }
    if (n > shfc->adir_nalloc)
    {
        shfc->adir_nalloc = over_alloc_dd(n);
        srenew(shfc->adir_xnold, shfc->adir_nalloc);
        srenew(shfc->adir_xnew, shfc->adir_nalloc);
    }
    xnold = shfc->adir_xnold;
    xnew  = shfc->adir_xnew;

    ptype = md->ptype;

    dt = ir->delta_t;

    /* Does NOT work with freeze or acceleration groups (yet) */
    for (n = 0; n < end; n++)
    {
        w_dt = md->invmass[n]*dt;

        for (d = 0; d < DIM; d++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell))
            {
                xnold[n][d] = x[n][d] - (x_init[n][d] - x_old[n][d]);
                xnew[n][d]  = 2*x[n][d] - x_old[n][d] + f[n][d]*w_dt*dt;
            }
            else
            {
                xnold[n][d] = x[n][d];
                xnew[n][d]  = x[n][d];
            }
        }
    }
    constrain(log, FALSE, FALSE, constr, idef, ir, cr, step, 0, 1.0, md,
              x, xnold, nullptr, bMolPBC, box,
              lambda[efptBONDED], &(dvdlambda[efptBONDED]),
              nullptr, nullptr, nrnb, econqCoord);
    constrain(log, FALSE, FALSE, constr, idef, ir, cr, step, 0, 1.0, md,
              x, xnew, nullptr, bMolPBC, box,
              lambda[efptBONDED], &(dvdlambda[efptBONDED]),
              nullptr, nullptr, nrnb, econqCoord);

    for (n = 0; n < end; n++)
    {
        for (d = 0; d < DIM; d++)
        {
            xnew[n][d] =
                -(2*x[n][d]-xnold[n][d]-xnew[n][d])/gmx::square(dt)
                - f[n][d]*md->invmass[n];
        }
        clear_rvec(acc_dir[n]);
    }

    /* Project the acceleration on the old bond directions */
    constrain(log, FALSE, FALSE, constr, idef, ir, cr, step, 0, 1.0, md,
              x_old, xnew, acc_dir, bMolPBC, box,
              lambda[efptBONDED], &(dvdlambda[efptBONDED]),
              nullptr, nullptr, nrnb, econqDeriv_FlexCon);
}

void relax_shell_flexcon(FILE *fplog, t_commrec *cr, gmx_bool bVerbose,
                         gmx_int64_t mdstep, t_inputrec *inputrec,
                         gmx_bool bDoNS, int force_flags,
                         gmx_localtop_t *top,
                         gmx_constr_t constr,
                         gmx_enerdata_t *enerd, t_fcdata *fcd,
                         t_state *state, PaddedRVecVector *f,
                         tensor force_vir,
                         t_mdatoms *md,
                         t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                         t_graph *graph,
                         gmx_groups_t *groups,
                         gmx_shellfc_t *shfc,
                         t_forcerec *fr,
                         gmx_bool bBornRadii,
                         double t, rvec mu_tot,
                         gmx_vsite_t *vsite,
                         DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
                         DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion)
{
    int        nshell;
    t_shell   *shell;
    t_idef    *idef;
    rvec      *acc_dir = nullptr, *x_old = nullptr;
    real       Epot[2], df[2];
    real       sf_dir, invdt;
    real       ftol, dum = 0;
    char       sbuf[22];
    gmx_bool   bCont, bInit, bConverged;
    int        nat, dd_ac0, dd_ac1 = 0, i;
    int        homenr = md->homenr, end = homenr, cg0, cg1;
    int        nflexcon, number_steps, d, Min = 0, count = 0;
#define  Try (1-Min)             /* At start Try = 1 */

    bCont        = (mdstep == inputrec->init_step) && inputrec->bContinuation;
    bInit        = (mdstep == inputrec->init_step) || shfc->bRequireInit;
    ftol         = inputrec->em_tol;
    number_steps = inputrec->niter;
    nshell       = shfc->nshell;
    shell        = shfc->shell;
    nflexcon     = shfc->nflexcon;

    idef = &top->idef;

    if (DOMAINDECOMP(cr))
    {
        nat = dd_natoms_vsite(cr->dd);
        if (nflexcon > 0)
        {
            dd_get_constraint_range(cr->dd, &dd_ac0, &dd_ac1);
            nat = std::max(nat, dd_ac1);
        }
    }
    else
    {
        nat = state->natoms;
    }

    for (i = 0; (i < 2); i++)
    {
        shfc->x[i].resize(gmx::paddedRVecVectorSize(nat));
        shfc->f[i].resize(gmx::paddedRVecVectorSize(nat));
    }

    /* Create views that we can swap */
    gmx::PaddedArrayRef<gmx::RVec> pos[2];
    gmx::PaddedArrayRef<gmx::RVec> force[2];
    for (i = 0; (i < 2); i++)
    {
        pos[i]   = shfc->x[i];
        force[i] = shfc->f[i];
    }

    if (bDoNS && inputrec->ePBC != epbcNONE && !DOMAINDECOMP(cr))
    {
        /* This is the only time where the coordinates are used
         * before do_force is called, which normally puts all
         * charge groups in the box.
         */
        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            auto xRef = makeArrayRef(state->x);
            put_atoms_in_box_omp(fr->ePBC, state->box, xRef.subArray(0, md->homenr));
        }
        else
        {
            cg0 = 0;
            cg1 = top->cgs.nr;
            put_charge_groups_in_box(fplog, cg0, cg1, fr->ePBC, state->box,
                                     &(top->cgs), as_rvec_array(state->x.data()), fr->cg_cm);
        }

        if (graph)
        {
            mk_mshift(fplog, graph, fr->ePBC, state->box, as_rvec_array(state->x.data()));
        }
    }

    /* After this all coordinate arrays will contain whole charge groups */
    if (graph)
    {
        shift_self(graph, state->box, as_rvec_array(state->x.data()));
    }

    if (nflexcon)
    {
        if (nat > shfc->flex_nalloc)
        {
            shfc->flex_nalloc = over_alloc_dd(nat);
            srenew(shfc->acc_dir, shfc->flex_nalloc);
            srenew(shfc->x_old, shfc->flex_nalloc);
        }
        acc_dir = shfc->acc_dir;
        x_old   = shfc->x_old;
        for (i = 0; i < homenr; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                shfc->x_old[i][d] =
                    state->x[i][d] - state->v[i][d]*inputrec->delta_t;
            }
        }
    }

    /* Do a prediction of the shell positions, when appropriate.
     * Without velocities (EM, NM, BD) we only do initial prediction.
     */
    if (shfc->bPredict && !bCont && (EI_STATE_VELOCITY(inputrec->eI) || bInit))
    {
        predict_shells(fplog, as_rvec_array(state->x.data()), as_rvec_array(state->v.data()), inputrec->delta_t, nshell, shell,
                       md->massT, nullptr, bInit);
    }

    /* do_force expected the charge groups to be in the box */
    if (graph)
    {
        unshift_self(graph, state->box, as_rvec_array(state->x.data()));
    }

    /* Calculate the forces first time around */
    if (gmx_debug_at)
    {
        pr_rvecs(debug, 0, "x b4 do_force", as_rvec_array(state->x.data()), homenr);
    }
    do_force(fplog, cr, inputrec, mdstep, nrnb, wcycle, top, groups,
             state->box, state->x, &state->hist,
             force[Min], force_vir, md, enerd, fcd,
             state->lambda, graph,
             fr, vsite, mu_tot, t, nullptr, bBornRadii,
             (bDoNS ? GMX_FORCE_NS : 0) | force_flags,
             ddOpenBalanceRegion, ddCloseBalanceRegion);

    sf_dir = 0;
    if (nflexcon)
    {
        init_adir(fplog, shfc,
                  constr, idef, inputrec, cr, dd_ac1, mdstep, md, end,
                  shfc->x_old, as_rvec_array(state->x.data()), as_rvec_array(state->x.data()), as_rvec_array(force[Min].data()),
                  shfc->acc_dir,
                  fr->bMolPBC, state->box, state->lambda, &dum, nrnb);

        for (i = 0; i < end; i++)
        {
            sf_dir += md->massT[i]*norm2(shfc->acc_dir[i]);
        }
    }

    Epot[Min] = enerd->term[F_EPOT];

    df[Min] = rms_force(cr, shfc->f[Min], nshell, shell, nflexcon, &sf_dir, &Epot[Min]);
    df[Try] = 0;
    if (debug)
    {
        fprintf(debug, "df = %g  %g\n", df[Min], df[Try]);
    }

    if (gmx_debug_at)
    {
        pr_rvecs(debug, 0, "force0", as_rvec_array(force[Min].data()), md->nr);
    }

    if (nshell+nflexcon > 0)
    {
        /* Copy x to pos[Min] & pos[Try]: during minimization only the
         * shell positions are updated, therefore the other particles must
         * be set here.
         */
        pos[Min] = state->x;
        pos[Try] = state->x;
    }

    if (bVerbose && MASTER(cr))
    {
        print_epot(stdout, mdstep, 0, Epot[Min], df[Min], nflexcon, sf_dir);
    }

    if (debug)
    {
        fprintf(debug, "%17s: %14.10e\n",
                interaction_function[F_EKIN].longname, enerd->term[F_EKIN]);
        fprintf(debug, "%17s: %14.10e\n",
                interaction_function[F_EPOT].longname, enerd->term[F_EPOT]);
        fprintf(debug, "%17s: %14.10e\n",
                interaction_function[F_ETOT].longname, enerd->term[F_ETOT]);
        fprintf(debug, "SHELLSTEP %s\n", gmx_step_str(mdstep, sbuf));
    }

    /* First check whether we should do shells, or whether the force is
     * low enough even without minimization.
     */
    bConverged = (df[Min] < ftol);

    for (count = 1; (!(bConverged) && (count < number_steps)); count++)
    {
        if (vsite)
        {
            construct_vsites(vsite, as_rvec_array(pos[Min].data()),
                             inputrec->delta_t, as_rvec_array(state->v.data()),
                             idef->iparams, idef->il,
                             fr->ePBC, fr->bMolPBC, cr, state->box);
        }

        if (nflexcon)
        {
            init_adir(fplog, shfc,
                      constr, idef, inputrec, cr, dd_ac1, mdstep, md, end,
                      x_old, as_rvec_array(state->x.data()), as_rvec_array(pos[Min].data()), as_rvec_array(force[Min].data()), acc_dir,
                      fr->bMolPBC, state->box, state->lambda, &dum, nrnb);

            directional_sd(pos[Min], pos[Try], acc_dir, end, fr->fc_stepsize);
        }

        /* New positions, Steepest descent */
        shell_pos_sd(pos[Min], pos[Try], force[Min], nshell, shell, count);

        /* do_force expected the charge groups to be in the box */
        if (graph)
        {
            unshift_self(graph, state->box, as_rvec_array(pos[Try].data()));
        }

        if (gmx_debug_at)
        {
            pr_rvecs(debug, 0, "RELAX: pos[Min]  ", as_rvec_array(pos[Min].data()), homenr);
            pr_rvecs(debug, 0, "RELAX: pos[Try]  ", as_rvec_array(pos[Try].data()), homenr);
        }
        /* Try the new positions */
        do_force(fplog, cr, inputrec, 1, nrnb, wcycle,
                 top, groups, state->box, pos[Try], &state->hist,
                 force[Try], force_vir,
                 md, enerd, fcd, state->lambda, graph,
                 fr, vsite, mu_tot, t, nullptr, bBornRadii,
                 force_flags,
                 ddOpenBalanceRegion, ddCloseBalanceRegion);

        if (gmx_debug_at)
        {
            pr_rvecs(debug, 0, "RELAX: force[Min]", as_rvec_array(force[Min].data()), homenr);
            pr_rvecs(debug, 0, "RELAX: force[Try]", as_rvec_array(force[Try].data()), homenr);
        }
        sf_dir = 0;
        if (nflexcon)
        {
            init_adir(fplog, shfc,
                      constr, idef, inputrec, cr, dd_ac1, mdstep, md, end,
                      x_old, as_rvec_array(state->x.data()), as_rvec_array(pos[Try].data()), as_rvec_array(force[Try].data()), acc_dir,
                      fr->bMolPBC, state->box, state->lambda, &dum, nrnb);

            for (i = 0; i < end; i++)
            {
                sf_dir += md->massT[i]*norm2(acc_dir[i]);
            }
        }

        Epot[Try] = enerd->term[F_EPOT];

        df[Try] = rms_force(cr, force[Try], nshell, shell, nflexcon, &sf_dir, &Epot[Try]);

        if (debug)
        {
            fprintf(debug, "df = %g  %g\n", df[Min], df[Try]);
        }

        if (debug)
        {
            if (gmx_debug_at)
            {
                pr_rvecs(debug, 0, "F na do_force", as_rvec_array(force[Try].data()), homenr);
            }
            if (gmx_debug_at)
            {
                fprintf(debug, "SHELL ITER %d\n", count);
                dump_shells(debug, pos[Try], force[Try], ftol, nshell, shell);
            }
        }

        if (bVerbose && MASTER(cr))
        {
            print_epot(stdout, mdstep, count, Epot[Try], df[Try], nflexcon, sf_dir);
        }

        bConverged = (df[Try] < ftol);

        if ((df[Try] < df[Min]))
        {
            if (debug)
            {
                fprintf(debug, "Swapping Min and Try\n");
            }
            if (nflexcon)
            {
                /* Correct the velocities for the flexible constraints */
                invdt = 1/inputrec->delta_t;
                for (i = 0; i < end; i++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        state->v[i][d] += (pos[Try][i][d] - pos[Min][i][d])*invdt;
                    }
                }
            }
            Min  = Try;
        }
        else
        {
            decrease_step_size(nshell, shell);
        }
    }
    shfc->numForceEvaluations += count;
    if (bConverged)
    {
        shfc->numConvergedIterations++;
    }
    if (MASTER(cr) && !(bConverged))
    {
        /* Note that the energies and virial are incorrect when not converged */
        if (fplog)
        {
            fprintf(fplog,
                    "step %s: EM did not converge in %d iterations, RMS force %.3f\n",
                    gmx_step_str(mdstep, sbuf), number_steps, df[Min]);
        }
        fprintf(stderr,
                "step %s: EM did not converge in %d iterations, RMS force %.3f\n",
                gmx_step_str(mdstep, sbuf), number_steps, df[Min]);
    }

    /* Copy back the coordinates and the forces */
    std::copy(pos[Min].begin(), pos[Min].end(), state->x.begin());
    std::copy(force[Min].begin(), force[Min].end(), f->begin());
}

void done_shellfc(FILE *fplog, gmx_shellfc_t *shfc, gmx_int64_t numSteps)
{
    if (shfc && fplog && numSteps > 0)
    {
        double numStepsAsDouble = static_cast<double>(numSteps);
        fprintf(fplog, "Fraction of iterations that converged:           %.2f %%\n",
                (shfc->numConvergedIterations*100.0)/numStepsAsDouble);
        fprintf(fplog, "Average number of force evaluations per MD step: %.2f\n\n",
                shfc->numForceEvaluations/numStepsAsDouble);
    }

    // TODO Deallocate memory in shfc
}
