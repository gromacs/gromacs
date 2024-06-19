/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "shellfc.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <utility>
#include <vector>

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/taskassignment/include/gromacs/taskassignment/decidesimulationworkload.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

using gmx::ArrayRef;
using gmx::ArrayRefWithPadding;
using gmx::RVec;

struct t_shell
{
    int  nnucl      = 0;  /* The number of nuclei */
    int  shellIndex = -1; /* The shell index */
    int  nucl1      = -1; /* The first nuclei connected to the shell	*/
    int  nucl2      = -1; /* The second nuclei connected to the shell	*/
    int  nucl3      = -1; /* The third nuclei connected to the shell	*/
    real k          = 0;  /* force constant		        */
    real k_1        = 0;  /* 1 over force constant		*/
    rvec xold;            /* The old shell coordinates */
    rvec fold;            /* The old force on the shell */
    rvec step;            /* Step size for steepest descents */
};

struct gmx_shellfc_t
{
    /* Shell counts, indices, parameters and working data */
    std::vector<t_shell> shell_gl;              /* All the shells (for DD only)              */
    std::vector<int>     shell_index_gl;        /* Global shell index (for DD only)          */
    gmx_bool             bInterCG;              /* Are there inter charge-group shells?      */
    std::vector<t_shell> shells;                /* The local shells                          */
    bool                 predictShells = false; /* Predict shell positions                   */
    bool                 requireInit   = false; /* Require initialization of shell positions */
    int                  nflexcon      = 0;     /* The number of flexible constraints        */

    std::array<PaddedHostVector<RVec>, 2> x; /* Coordinate buffers for iterative minimization */
    std::array<PaddedHostVector<RVec>, 2> f; /* Force buffers for iterative minimization */

    /* Flexible constraint working data */
    std::vector<RVec>       acc_dir;                /* Acceleration direction for flexcon        */
    gmx::PaddedVector<RVec> x_old;                  /* Old coordinates for flexcon               */
    gmx::PaddedVector<RVec> adir_xnold;             /* Work space for init_adir                  */
    gmx::PaddedVector<RVec> adir_xnew;              /* Work space for init_adir                  */
    std::int64_t            numForceEvaluations;    /* Total number of force evaluations         */
    int                     numConvergedIterations; /* Total number of iterations that converged */
};


static void pr_shell(FILE* fplog, ArrayRef<const t_shell> shells)
{
    fprintf(fplog, "SHELL DATA\n");
    fprintf(fplog, "%5s  %8s  %5s  %5s  %5s\n", "Shell", "Force k", "Nucl1", "Nucl2", "Nucl3");
    for (const t_shell& shell : shells)
    {
        fprintf(fplog, "%5d  %8.3f  %5d", shell.shellIndex, 1.0 / shell.k_1, shell.nucl1);
        if (shell.nnucl == 2)
        {
            fprintf(fplog, "  %5d\n", shell.nucl2);
        }
        else if (shell.nnucl == 3)
        {
            fprintf(fplog, "  %5d  %5d\n", shell.nucl2, shell.nucl3);
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
static void predict_shells(FILE*                     fplog,
                           ArrayRef<RVec>            x,
                           ArrayRef<RVec>            v,
                           real                      dt,
                           ArrayRef<const t_shell>   shells,
                           gmx::ArrayRef<const real> mass,
                           gmx_bool                  bInit)
{
    int  m, n1, n2, n3;
    real dt_1, fudge, tm, m1, m2, m3;

    /* We introduce a fudge factor for performance reasons: with this choice
     * the initial force on the shells is about a factor of two lower than
     * without
     */
    fudge = 1.0;

    ArrayRef<RVec> xOrV;
    if (bInit)
    {
        if (fplog)
        {
            fprintf(fplog, "RELAX: Using prediction for initial shell placement\n");
        }
        xOrV = x;
        dt_1 = 1;
    }
    else
    {
        xOrV = v;
        dt_1 = fudge * dt;
    }

    for (const t_shell& shell : shells)
    {
        const int s1 = shell.shellIndex;
        if (bInit)
        {
            clear_rvec(x[s1]);
        }
        switch (shell.nnucl)
        {
            case 1:
                n1 = shell.nucl1;
                for (m = 0; (m < DIM); m++)
                {
                    x[s1][m] += xOrV[n1][m] * dt_1;
                }
                break;
            case 2:
                n1 = shell.nucl1;
                n2 = shell.nucl2;
                m1 = mass[n1];
                m2 = mass[n2];
                tm = dt_1 / (m1 + m2);
                for (m = 0; (m < DIM); m++)
                {
                    x[s1][m] += (m1 * xOrV[n1][m] + m2 * xOrV[n2][m]) * tm;
                }
                break;
            case 3:
                n1 = shell.nucl1;
                n2 = shell.nucl2;
                n3 = shell.nucl3;
                m1 = mass[n1];
                m2 = mass[n2];
                m3 = mass[n3];
                tm = dt_1 / (m1 + m2 + m3);
                for (m = 0; (m < DIM); m++)
                {
                    x[s1][m] += (m1 * xOrV[n1][m] + m2 * xOrV[n2][m] + m3 * xOrV[n3][m]) * tm;
                }
                break;
            default: gmx_fatal(FARGS, "Shell %d has %d nuclei!", s1, shell.nnucl);
        }
    }
}

gmx_shellfc_t* init_shell_flexcon(FILE*             fplog,
                                  const gmx_mtop_t& mtop,
                                  int               nflexcon,
                                  int               nstcalcenergy,
                                  bool              usingDomainDecomposition,
                                  bool              usingPmeOnGpu)
{
    gmx_shellfc_t* shfc;

    int  ns, nshell, nsi;
    int  i, j, type, a_offset, mol, ftype, nra;
    real qS, alpha;
    int  aS, aN = 0; /* Shell and nucleus */
    int bondtypes[] = { F_BONDS, F_HARMONIC, F_CUBICBONDS, F_POLARIZATION, F_ANHARM_POL, F_WATER_POL };
#define NBT asize(bondtypes)
    const gmx_ffparams_t* ffparams;

    const gmx::EnumerationArray<ParticleType, int> numParticles = gmx_mtop_particletype_count(mtop);
    if (fplog)
    {
        /* Print the number of each particle type */
        for (const auto entry : gmx::keysOf(numParticles))
        {
            const int number = numParticles[entry];
            if (number != 0)
            {
                fprintf(fplog, "There are: %d %ss\n", number, enumValueToString(entry));
            }
        }
    }

    nshell = numParticles[ParticleType::Shell];

    if (nshell == 0 && nflexcon == 0)
    {
        /* We're not doing shells or flexible constraints */
        return nullptr;
    }

    shfc           = new gmx_shellfc_t;
    shfc->nflexcon = nflexcon;

    if (nshell == 0)
    {
        /* Only flexible constraints, no shells.
         * Note that make_local_shells() does not need to be called.
         */
        return shfc;
    }

    if (nstcalcenergy != 1)
    {
        gmx_fatal(FARGS,
                  "You have nstcalcenergy set to a value (%d) that is different from 1.\nThis is "
                  "not supported in combination with shell particles.\nPlease make a new tpr file.",
                  nstcalcenergy);
    }
    if (nshell > 0 && usingDomainDecomposition)
    {
        gmx_fatal(
                FARGS,
                "Shell particles are not implemented with domain decomposition, use a single rank");
    }

    /* We have shells: fill the shell data structure */

    /* Global system sized array, this should be avoided */
    std::vector<int> shell_index(mtop.natoms);

    nshell = 0;
    for (const AtomProxy atomP : AtomRange(mtop))
    {
        const t_atom& local = atomP.atom();
        int           i     = atomP.globalAtomNumber();
        if (local.ptype == ParticleType::Shell)
        {
            shell_index[i] = nshell++;
        }
    }

    std::vector<t_shell> shell(nshell);

    ffparams = &mtop.ffparams;

    /* Now fill the structures */
    /* TODO: See if we can use update groups that cover shell constructions */
    shfc->bInterCG = FALSE;
    ns             = 0;
    a_offset       = 0;
    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb = mtop.molblock[mb];
        const gmx_moltype_t&  molt = mtop.moltype[molb.type];

        const t_atom* atom = molt.atoms.atom;
        for (mol = 0; mol < molb.nmol; mol++)
        {
            for (j = 0; (j < NBT); j++)
            {
                const int* ia = molt.ilist[bondtypes[j]].iatoms.data();
                for (i = 0; (i < molt.ilist[bondtypes[j]].size());)
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
                            if (atom[ia[1]].ptype == ParticleType::Shell)
                            {
                                aS = ia[1];
                                aN = ia[2];
                            }
                            else if (atom[ia[2]].ptype == ParticleType::Shell)
                            {
                                aS = ia[2];
                                aN = ia[1];
                            }
                            break;
                        case F_WATER_POL:
                            aN = ia[4]; /* Dummy */
                            aS = ia[5]; /* Shell */
                            break;
                        default: gmx_fatal(FARGS, "Death Horror: %s, %d", __FILE__, __LINE__);
                    }

                    if (aS != -1)
                    {
                        qS = atom[aS].q;

                        /* Check whether one of the particles is a shell... */
                        nsi = shell_index[a_offset + aS];
                        if ((nsi < 0) || (nsi >= nshell))
                        {
                            gmx_fatal(FARGS, "nsi is %d should be within 0 - %d. aS = %d", nsi, nshell, aS);
                        }
                        if (shell[nsi].shellIndex == -1)
                        {
                            shell[nsi].shellIndex = a_offset + aS;
                            ns++;
                        }
                        else if (shell[nsi].shellIndex != a_offset + aS)
                        {
                            gmx_fatal(FARGS, "Weird stuff in %s, %d", __FILE__, __LINE__);
                        }

                        if (shell[nsi].nucl1 == -1)
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
                                pr_shell(fplog, shell);
                            }
                            gmx_fatal(FARGS, "Can not handle more than three bonds per shell\n");
                        }
                        if (aS != aN)
                        {
                            /* shell[nsi].bInterCG = TRUE; */
                            shfc->bInterCG = TRUE;
                        }

                        switch (bondtypes[j])
                        {
                            case F_BONDS:
                            case F_HARMONIC:
                                shell[nsi].k += ffparams->iparams[type].harmonic.krA;
                                break;
                            case F_CUBICBONDS:
                                shell[nsi].k += ffparams->iparams[type].cubic.kb;
                                break;
                            case F_POLARIZATION:
                            case F_ANHARM_POL:
                                if (!gmx_within_tol(qS, atom[aS].qB, GMX_REAL_EPS * 10))
                                {
                                    gmx_fatal(FARGS,
                                              "polarize can not be used with qA(%e) != qB(%e) for "
                                              "atom %d of molecule block %zu",
                                              qS,
                                              atom[aS].qB,
                                              aS + 1,
                                              mb + 1);
                                }
                                shell[nsi].k += gmx::square(qS) * gmx::c_one4PiEps0
                                                / ffparams->iparams[type].polarize.alpha;
                                break;
                            case F_WATER_POL:
                                if (!gmx_within_tol(qS, atom[aS].qB, GMX_REAL_EPS * 10))
                                {
                                    gmx_fatal(FARGS,
                                              "water_pol can not be used with qA(%e) != qB(%e) for "
                                              "atom %d of molecule block %zu",
                                              qS,
                                              atom[aS].qB,
                                              aS + 1,
                                              mb + 1);
                                }
                                alpha = (ffparams->iparams[type].wpol.al_x
                                         + ffparams->iparams[type].wpol.al_y
                                         + ffparams->iparams[type].wpol.al_z)
                                        / 3.0;
                                shell[nsi].k += gmx::square(qS) * gmx::c_one4PiEps0 / alpha;
                                break;
                            default: gmx_fatal(FARGS, "Death Horror: %s, %d", __FILE__, __LINE__);
                        }
                        shell[nsi].nnucl++;
                    }
                    ia += nra + 1;
                    i += nra + 1;
                }
            }
            a_offset += molt.atoms.nr;
        }
        /* Done with this molecule type */
    }

    /* Verify whether it's all correct */
    if (ns != nshell)
    {
        gmx_fatal(FARGS, "Something weird with shells. They may not be bonded to something");
    }

    for (i = 0; (i < ns); i++)
    {
        shell[i].k_1 = 1.0 / shell[i].k;
    }

    if (debug)
    {
        pr_shell(debug, shell);
    }


    shfc->shell_gl       = shell;
    shfc->shell_index_gl = shell_index;

    shfc->predictShells = (getenv("GMX_NOPREDICT") == nullptr);
    shfc->requireInit   = false;
    if (!shfc->predictShells)
    {
        if (fplog)
        {
            fprintf(fplog, "\nWill never predict shell positions\n");
        }
    }
    else
    {
        shfc->requireInit = (getenv("GMX_REQUIRE_SHELL_INIT") != nullptr);
        if (shfc->requireInit && fplog)
        {
            fprintf(fplog, "\nWill always initiate shell positions\n");
        }
    }

    if (shfc->predictShells)
    {
        if (shfc->bInterCG)
        {
            if (fplog)
            {
                fprintf(fplog,
                        "\nNOTE: in the current version shell prediction during the crun is "
                        "disabled\n\n");
            }
            /* Prediction improves performance, so we should implement either:
             * 1. communication for the atoms needed for prediction
             * 2. prediction using the velocities of shells; currently the
             *    shell velocities are zeroed, it's a bit tricky to keep
             *    track of the shell displacements and thus the velocity.
             */
            shfc->predictShells = false;
        }
    }

    /* shfc->x is used as a coordinate buffer for the sim_util's `do_force` function, and
     * when using PME it must be pinned. */
    if (usingPmeOnGpu)
    {
        for (i = 0; i < 2; i++)
        {
            changePinningPolicy(&shfc->x[i], gmx::PinningPolicy::PinnedIfSupported);
        }
    }

    return shfc;
}

void gmx::make_local_shells(const t_commrec* cr, const t_mdatoms& md, gmx_shellfc_t* shfc)
{
    int           a0, a1;
    gmx_domdec_t* dd = nullptr;

    if (haveDDAtomOrdering(*cr))
    {
        dd = cr->dd;
        a0 = 0;
        a1 = dd_numHomeAtoms(*dd);
    }
    else
    {
        /* Single node: we need all shells, copy them */
        shfc->shells = shfc->shell_gl;

        return;
    }

    ArrayRef<const int> ind = shfc->shell_index_gl;

    std::vector<t_shell>& shells = shfc->shells;
    shells.clear();
    for (int i = a0; i < a1; i++)
    {
        if (md.ptype[i] == ParticleType::Shell)
        {
            if (dd)
            {
                shells.push_back(shfc->shell_gl[ind[dd->globalAtomIndices[i]]]);
            }
            else
            {
                shells.push_back(shfc->shell_gl[ind[i]]);
            }
            t_shell& shell = shells.back();

            /* With inter-cg shells we can no do shell prediction,
             * so we do not need the nuclei numbers.
             */
            if (!shfc->bInterCG)
            {
                shell.nucl1 = i + shell.nucl1 - shell.shellIndex;
                if (shell.nnucl > 1)
                {
                    shell.nucl2 = i + shell.nucl2 - shell.shellIndex;
                }
                if (shell.nnucl > 2)
                {
                    shell.nucl3 = i + shell.nucl3 - shell.shellIndex;
                }
            }
            shell.shellIndex = i;
        }
    }
}

static void do_1pos(rvec xnew, const rvec xold, const rvec f, real step)
{
    real xo, yo, zo;
    real dx, dy, dz;

    xo = xold[XX];
    yo = xold[YY];
    zo = xold[ZZ];

    dx = f[XX] * step;
    dy = f[YY] * step;
    dz = f[ZZ] * step;

    xnew[XX] = xo + dx;
    xnew[YY] = yo + dy;
    xnew[ZZ] = zo + dz;
}

static void do_1pos3(rvec xnew, const rvec xold, const rvec f, const rvec step)
{
    real xo, yo, zo;
    real dx, dy, dz;

    xo = xold[XX];
    yo = xold[YY];
    zo = xold[ZZ];

    dx = f[XX] * step[XX];
    dy = f[YY] * step[YY];
    dz = f[ZZ] * step[ZZ];

    xnew[XX] = xo + dx;
    xnew[YY] = yo + dy;
    xnew[ZZ] = zo + dz;
}

static void directional_sd(ArrayRef<const RVec> xold,
                           ArrayRef<RVec>       xnew,
                           ArrayRef<const RVec> acc_dir,
                           int                  homenr,
                           real                 step)
{
    const rvec* xo = as_rvec_array(xold.data());
    rvec*       xn = as_rvec_array(xnew.data());

    for (int i = 0; i < homenr; i++)
    {
        do_1pos(xn[i], xo[i], acc_dir[i], step);
    }
}

static void shell_pos_sd(ArrayRef<const RVec> xcur,
                         ArrayRef<RVec>       xnew,
                         ArrayRef<const RVec> f,
                         ArrayRef<t_shell>    shells,
                         int                  count)
{
    const real step_scale_min = 0.8, step_scale_increment = 0.2, step_scale_max = 1.2,
               step_scale_multiple = (step_scale_max - step_scale_min) / step_scale_increment;
    int        d;
    real       dx, df, k_est;
    const real zero = 0;
#ifdef PRINT_STEP
    real step_min, step_max;

    step_min = 1e30;
    step_max = 0;
#endif
    for (t_shell& shell : shells)
    {
        const int ind = shell.shellIndex;
        if (count == 1)
        {
            for (d = 0; d < DIM; d++)
            {
                shell.step[d] = shell.k_1;
#ifdef PRINT_STEP
                step_min = std::min(step_min, shell.step[d]);
                step_max = std::max(step_max, shell.step[d]);
#endif
            }
        }
        else
        {
            for (d = 0; d < DIM; d++)
            {
                dx = xcur[ind][d] - shell.xold[d];
                df = f[ind][d] - shell.fold[d];
                /* -dx/df gets used to generate an interpolated value, but would
                 * cause a NaN if df were binary-equal to zero. Values close to
                 * zero won't cause problems (because of the min() and max()), so
                 * just testing for binary inequality is OK. */
                if (zero != df)
                {
                    k_est = -dx / df;
                    /* Scale the step size by a factor interpolated from
                     * step_scale_min to step_scale_max, as k_est goes from 0 to
                     * step_scale_multiple * shell.step[d] */
                    shell.step[d] = step_scale_min * shell.step[d]
                                    + step_scale_increment
                                              * std::min(step_scale_multiple * shell.step[d],
                                                         std::max(k_est, zero));
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
                        shell.step[d] *= step_scale_max;
                    }
                }
#ifdef PRINT_STEP
                step_min = std::min(step_min, shell.step[d]);
                step_max = std::max(step_max, shell.step[d]);
#endif
            }
        }
        copy_rvec(xcur[ind], shell.xold);
        copy_rvec(f[ind], shell.fold);

        do_1pos3(xnew[ind], xcur[ind], f[ind], shell.step);

        if (gmx_debug_at)
        {
            fprintf(debug, "shell = %d\n", ind);
            pr_rvec(debug, 0, "fshell", f[ind], DIM, TRUE);
            pr_rvec(debug, 0, "xold", xcur[ind], DIM, TRUE);
            pr_rvec(debug, 0, "step", shell.step, DIM, TRUE);
            pr_rvec(debug, 0, "xnew", xnew[ind], DIM, TRUE);
        }
    }
#ifdef PRINT_STEP
    printf("step %.3e %.3e\n", step_min, step_max);
#endif
}

static void decrease_step_size(ArrayRef<t_shell> shells)
{
    for (t_shell& shell : shells)
    {
        svmul(0.8, shell.step, shell.step);
    }
}

static void print_epot(FILE* fp, int64_t mdstep, int count, real epot, real df, int ndir, real sf_dir)
{
    char buf[22];

    fprintf(fp, "MDStep=%5s/%2d EPot: %12.8e, rmsF: %6.2e", gmx_step_str(mdstep, buf), count, epot, df);
    if (ndir)
    {
        fprintf(fp, ", dir. rmsF: %6.2e\n", std::sqrt(sf_dir / ndir));
    }
    else
    {
        fprintf(fp, "\n");
    }
}


static real rms_force(const t_commrec*        cr,
                      ArrayRef<const RVec>    force,
                      ArrayRef<const t_shell> shells,
                      int                     ndir,
                      real*                   sf_dir,
                      real*                   Epot)
{
    double      buf[4];
    const rvec* f = as_rvec_array(force.data());

    buf[0] = *sf_dir;
    for (const t_shell& shell : shells)
    {
        buf[0] += norm2(f[shell.shellIndex]);
    }
    int ntot = shells.ssize();

    if (PAR(cr))
    {
        buf[1] = ntot;
        buf[2] = *sf_dir;
        buf[3] = *Epot;
        gmx_sumd(4, buf, cr);
        ntot    = gmx::roundToInt(buf[1]);
        *sf_dir = buf[2];
        *Epot   = buf[3];
    }
    ntot += ndir;

    return (ntot ? std::sqrt(buf[0] / ntot) : 0);
}

static void dump_shells(FILE* fp, ArrayRef<RVec> f, real ftol, ArrayRef<const t_shell> shells)
{
    real ft2, ff2;

    ft2 = gmx::square(ftol);

    for (const t_shell& shell : shells)
    {
        const int ind = shell.shellIndex;
        ff2           = iprod(f[ind], f[ind]);
        if (ff2 > ft2)
        {
            fprintf(fp,
                    "SHELL %5d, force %10.5f  %10.5f  %10.5f, |f| %10.5f\n",
                    ind,
                    f[ind][XX],
                    f[ind][YY],
                    f[ind][ZZ],
                    std::sqrt(ff2));
        }
    }
}

static void init_adir(gmx_shellfc_t*            shfc,
                      gmx::Constraints*         constr,
                      const t_inputrec*         ir,
                      const t_commrec*          cr,
                      int                       dd_ac1,
                      int64_t                   step,
                      const t_mdatoms&          md,
                      int                       end,
                      ArrayRefWithPadding<RVec> xOld,
                      ArrayRef<RVec>            x_init,
                      ArrayRefWithPadding<RVec> xCurrent,
                      ArrayRef<RVec>            f,
                      ArrayRef<RVec>            acc_dir,
                      const matrix              box,
                      ArrayRef<const real>      lambda,
                      real*                     dvdlambda)
{
    double dt, w_dt;
    int    n, d;

    if (haveDDAtomOrdering(*cr))
    {
        n = dd_ac1;
    }
    else
    {
        n = end;
    }
    shfc->adir_xnold.resizeWithPadding(n);
    shfc->adir_xnew.resizeWithPadding(n);
    rvec* xnold = as_rvec_array(shfc->adir_xnold.data());
    rvec* xnew  = as_rvec_array(shfc->adir_xnew.data());
    rvec* x_old = as_rvec_array(xOld.paddedArrayRef().data());
    rvec* x     = as_rvec_array(xCurrent.paddedArrayRef().data());

    dt = ir->delta_t;

    /* Does NOT work with freeze or acceleration groups (yet) */
    for (n = 0; n < end; n++)
    {
        w_dt = md.invmass[n] * dt;

        for (d = 0; d < DIM; d++)
        {
            if ((md.ptype[n] != ParticleType::VSite) && (md.ptype[n] != ParticleType::Shell))
            {
                xnold[n][d] = x[n][d] - (x_init[n][d] - x_old[n][d]);
                xnew[n][d]  = 2 * x[n][d] - x_old[n][d] + f[n][d] * w_dt * dt;
            }
            else
            {
                xnold[n][d] = x[n][d];
                xnew[n][d]  = x[n][d];
            }
        }
    }
    bool computeRmsd   = false;
    bool computeVirial = false;
    constr->apply(computeRmsd,
                  step,
                  0,
                  1.0,
                  xCurrent,
                  shfc->adir_xnold.arrayRefWithPadding(),
                  {},
                  box,
                  lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)],
                  &(dvdlambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)]),
                  {},
                  computeVirial,
                  nullptr,
                  gmx::ConstraintVariable::Positions);
    constr->apply(computeRmsd,
                  step,
                  0,
                  1.0,
                  xCurrent,
                  shfc->adir_xnew.arrayRefWithPadding(),
                  {},
                  box,
                  lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)],
                  &(dvdlambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)]),
                  {},
                  computeVirial,
                  nullptr,
                  gmx::ConstraintVariable::Positions);

    for (n = 0; n < end; n++)
    {
        for (d = 0; d < DIM; d++)
        {
            xnew[n][d] = -(2 * x[n][d] - xnold[n][d] - xnew[n][d]) / gmx::square(dt)
                         - f[n][d] * md.invmass[n];
        }
        clear_rvec(acc_dir[n]);
    }

    /* Project the acceleration on the old bond directions */
    constr->apply(computeRmsd,
                  step,
                  0,
                  1.0,
                  xOld,
                  shfc->adir_xnew.arrayRefWithPadding(),
                  acc_dir,
                  box,
                  lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)],
                  &(dvdlambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)]),
                  {},
                  computeVirial,
                  nullptr,
                  gmx::ConstraintVariable::Deriv_FlexCon);
}

void relax_shell_flexcon(FILE*                             fplog,
                         const t_commrec*                  cr,
                         const gmx_multisim_t*             ms,
                         gmx_bool                          bVerbose,
                         gmx_enfrot*                       enforcedRotation,
                         int64_t                           mdstep,
                         const t_inputrec*                 inputrec,
                         const gmx::MDModulesNotifiers&    mdModulesNotifiers,
                         gmx::ImdSession*                  imdSession,
                         pull_t*                           pull_work,
                         gmx_bool                          bDoNS,
                         const gmx_localtop_t*             top,
                         gmx::Constraints*                 constr,
                         gmx_enerdata_t*                   enerd,
                         int                               natoms,
                         ArrayRefWithPadding<RVec>         xPadded,
                         ArrayRefWithPadding<RVec>         vPadded,
                         const matrix                      box,
                         ArrayRef<real>                    lambda,
                         const history_t*                  hist,
                         gmx::ForceBuffersView*            f,
                         tensor                            force_vir,
                         const t_mdatoms&                  md,
                         CpuPpLongRangeNonbondeds*         longRangeNonbondeds,
                         t_nrnb*                           nrnb,
                         gmx_wallcycle*                    wcycle,
                         gmx_shellfc_t*                    shfc,
                         t_forcerec*                       fr,
                         const gmx::MdrunScheduleWorkload& runScheduleWork,
                         double                            t,
                         rvec                              mu_tot,
                         gmx::VirtualSitesHandler*         vsite,
                         const DDBalanceRegionHandler&     ddBalanceRegionHandler)
{
    real Epot[2], df[2];
    real sf_dir, invdt;
    real dum = 0;
    char sbuf[22];
    int  nat, dd_ac0, dd_ac1 = 0, i;
    int  homenr = md.homenr, end = homenr;
    int  d, Min = 0, count = 0;
#define Try (1 - Min) /* At start Try = 1 */

    const bool        bCont        = (mdstep == inputrec->init_step) && inputrec->bContinuation;
    const bool        bInit        = (mdstep == inputrec->init_step) || shfc->requireInit;
    const real        ftol         = inputrec->em_tol;
    const int         number_steps = inputrec->niter;
    ArrayRef<t_shell> shells       = shfc->shells;
    const int         nflexcon     = shfc->nflexcon;

    if (haveDDAtomOrdering(*cr))
    {
        nat = dd_natoms_vsite(*cr->dd);
        if (nflexcon > 0)
        {
            dd_get_constraint_range(*cr->dd, &dd_ac0, &dd_ac1);
            nat = std::max(nat, dd_ac1);
        }
    }
    else
    {
        nat = natoms;
    }

    for (i = 0; (i < 2); i++)
    {
        shfc->x[i].resizeWithPadding(nat);
        shfc->f[i].resizeWithPadding(nat);
    }

    /* Create views that we can swap for trail and minimum for positions and forces */
    ArrayRefWithPadding<RVec> posWithPadding[2];
    ArrayRefWithPadding<RVec> forceWithPadding[2];
    ArrayRef<RVec>            pos[2];
    ArrayRef<RVec>            force[2];
    for (i = 0; (i < 2); i++)
    {
        posWithPadding[i]   = shfc->x[i].arrayRefWithPadding();
        pos[i]              = posWithPadding[i].paddedArrayRef();
        forceWithPadding[i] = shfc->f[i].arrayRefWithPadding();
        force[i]            = forceWithPadding[i].paddedArrayRef();
    }

    ArrayRef<RVec> x = xPadded.unpaddedArrayRef();
    ArrayRef<RVec> v = vPadded.unpaddedArrayRef();

    if (bDoNS && inputrec->pbcType != PbcType::No && !haveDDAtomOrdering(*cr))
    {
        /* This is the only time where the coordinates are used
         * before do_force is called, which normally puts all
         * charge groups in the box.
         */
        put_atoms_in_box_omp(fr->pbcType,
                             box,
                             fr->haveBoxDeformation,
                             inputrec->deform,
                             x.subArray(0, md.homenr),
                             v.empty() ? ArrayRef<RVec>() : v.subArray(0, md.homenr),
                             gmx_omp_nthreads_get(ModuleMultiThread::Default));
    }

    if (nflexcon)
    {
        shfc->acc_dir.resize(nat);
        shfc->x_old.resizeWithPadding(nat);
        ArrayRef<RVec> x_old = shfc->x_old.arrayRefWithPadding().unpaddedArrayRef();
        for (i = 0; i < homenr; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                x_old[i][d] = x[i][d] - v[i][d] * inputrec->delta_t;
            }
        }
    }

    auto massT = md.massT;
    /* Do a prediction of the shell positions, when appropriate.
     * Without velocities (EM, NM, BD) we only do initial prediction.
     */
    if (shfc->predictShells && !bCont && (EI_STATE_VELOCITY(inputrec->eI) || bInit))
    {
        predict_shells(fplog, x, v, inputrec->delta_t, shells, massT, bInit);
    }

    /* Calculate the forces first time around */
    if (gmx_debug_at)
    {
        pr_rvecs(debug, 0, "x b4 do_force", as_rvec_array(x.data()), homenr);
    }
    gmx::ForceBuffersView forceViewInit = gmx::ForceBuffersView(forceWithPadding[Min], {}, false);

    do_force(fplog,
             cr,
             ms,
             *inputrec,
             mdModulesNotifiers,
             nullptr,
             enforcedRotation,
             imdSession,
             pull_work,
             mdstep,
             nrnb,
             wcycle,
             top,
             box,
             xPadded,
             vPadded.unpaddedArrayRef(),
             hist,
             &forceViewInit,
             force_vir,
             &md,
             enerd,
             lambda,
             fr,
             runScheduleWork,
             vsite,
             mu_tot,
             t,
             nullptr,
             longRangeNonbondeds,
             ddBalanceRegionHandler);

    sf_dir = 0;
    if (nflexcon)
    {
        init_adir(shfc,
                  constr,
                  inputrec,
                  cr,
                  dd_ac1,
                  mdstep,
                  md,
                  end,
                  shfc->x_old.arrayRefWithPadding(),
                  x,
                  xPadded,
                  force[Min],
                  shfc->acc_dir,
                  box,
                  lambda,
                  &dum);

        for (i = 0; i < end; i++)
        {
            sf_dir += massT[i] * norm2(shfc->acc_dir[i]);
        }
    }
    accumulatePotentialEnergies(enerd, lambda, inputrec->fepvals.get());
    Epot[Min] = enerd->term[F_EPOT];

    df[Min] = rms_force(cr, forceWithPadding[Min].paddedArrayRef(), shells, nflexcon, &sf_dir, &Epot[Min]);
    df[Try] = 0;
    if (debug)
    {
        fprintf(debug, "df = %g  %g\n", df[Min], df[Try]);
    }

    if (gmx_debug_at)
    {
        pr_rvecs(debug, 0, "force0", as_rvec_array(force[Min].data()), md.nr);
    }

    if (!shells.empty() || nflexcon > 0)
    {
        /* Copy x to pos[Min] & pos[Try]: during minimization only the
         * shell positions are updated, therefore the other particles must
         * be set here, in advance.
         */
        std::copy(xPadded.paddedArrayRef().begin(),
                  xPadded.paddedArrayRef().end(),
                  posWithPadding[Min].paddedArrayRef().begin());
        std::copy(xPadded.paddedArrayRef().begin(),
                  xPadded.paddedArrayRef().end(),
                  posWithPadding[Try].paddedArrayRef().begin());
    }

    if (bVerbose && MAIN(cr))
    {
        print_epot(stdout, mdstep, 0, Epot[Min], df[Min], nflexcon, sf_dir);
    }

    if (debug)
    {
        fprintf(debug, "%17s: %14.10e\n", interaction_function[F_EKIN].longname, enerd->term[F_EKIN]);
        fprintf(debug, "%17s: %14.10e\n", interaction_function[F_EPOT].longname, enerd->term[F_EPOT]);
        fprintf(debug, "%17s: %14.10e\n", interaction_function[F_ETOT].longname, enerd->term[F_ETOT]);
        fprintf(debug, "SHELLSTEP %s\n", gmx_step_str(mdstep, sbuf));
    }

    // For subsequent calls of do_force() on search steps we need to turn off the
    // the corresponding stepWork flag to avoid executing any (remaining) search-related
    // operations in do_force().
    // This copy can be removed when the doPairSearch() call is moved out of do_force().
    gmx::MdrunScheduleWorkload runScheduleWorkWithoutNS = runScheduleWork;
    runScheduleWorkWithoutNS.stepWork.doNeighborSearch  = false;

    /* First check whether we should do shells, or whether the force is
     * low enough even without minimization.
     */
    bool bConverged = (df[Min] < ftol);

    for (count = 1; (!(bConverged) && (count < number_steps)); count++)
    {
        if (vsite)
        {
            vsite->construct(pos[Min], v, box, gmx::VSiteOperation::PositionsAndVelocities);
        }

        if (nflexcon)
        {
            init_adir(shfc,
                      constr,
                      inputrec,
                      cr,
                      dd_ac1,
                      mdstep,
                      md,
                      end,
                      shfc->x_old.arrayRefWithPadding(),
                      x,
                      posWithPadding[Min],
                      force[Min],
                      shfc->acc_dir,
                      box,
                      lambda,
                      &dum);

            directional_sd(pos[Min], pos[Try], shfc->acc_dir, end, fr->fc_stepsize);
        }

        /* New positions, Steepest descent */
        shell_pos_sd(pos[Min], pos[Try], force[Min], shells, count);

        if (gmx_debug_at)
        {
            pr_rvecs(debug, 0, "RELAX: pos[Min]  ", as_rvec_array(pos[Min].data()), homenr);
            pr_rvecs(debug, 0, "RELAX: pos[Try]  ", as_rvec_array(pos[Try].data()), homenr);
        }
        /* Try the new positions */
        gmx::ForceBuffersView forceViewTry = gmx::ForceBuffersView(forceWithPadding[Try], {}, false);

        do_force(fplog,
                 cr,
                 ms,
                 *inputrec,
                 mdModulesNotifiers,
                 nullptr,
                 enforcedRotation,
                 imdSession,
                 pull_work,
                 1,
                 nrnb,
                 wcycle,
                 top,
                 box,
                 posWithPadding[Try],
                 {},
                 hist,
                 &forceViewTry,
                 force_vir,
                 &md,
                 enerd,
                 lambda,
                 fr,
                 runScheduleWorkWithoutNS,
                 vsite,
                 mu_tot,
                 t,
                 nullptr,
                 longRangeNonbondeds,
                 ddBalanceRegionHandler);
        accumulatePotentialEnergies(enerd, lambda, inputrec->fepvals.get());
        if (gmx_debug_at)
        {
            pr_rvecs(debug, 0, "RELAX: force[Min]", as_rvec_array(force[Min].data()), homenr);
            pr_rvecs(debug, 0, "RELAX: force[Try]", as_rvec_array(force[Try].data()), homenr);
        }
        sf_dir = 0;
        if (nflexcon)
        {
            init_adir(shfc,
                      constr,
                      inputrec,
                      cr,
                      dd_ac1,
                      mdstep,
                      md,
                      end,
                      shfc->x_old.arrayRefWithPadding(),
                      x,
                      posWithPadding[Try],
                      force[Try],
                      shfc->acc_dir,
                      box,
                      lambda,
                      &dum);

            ArrayRef<const RVec> acc_dir = shfc->acc_dir;
            for (i = 0; i < end; i++)
            {
                sf_dir += massT[i] * norm2(acc_dir[i]);
            }
        }

        Epot[Try] = enerd->term[F_EPOT];

        df[Try] = rms_force(cr, force[Try], shells, nflexcon, &sf_dir, &Epot[Try]);

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
                dump_shells(debug, force[Try], ftol, shells);
            }
        }

        if (bVerbose && MAIN(cr))
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
                invdt = 1 / inputrec->delta_t;
                for (i = 0; i < end; i++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        v[i][d] += (pos[Try][i][d] - pos[Min][i][d]) * invdt;
                    }
                }
            }
            Min = Try;
        }
        else
        {
            decrease_step_size(shells);
        }
    }
    shfc->numForceEvaluations += count;
    if (bConverged)
    {
        shfc->numConvergedIterations++;
    }
    if (MAIN(cr) && !(bConverged))
    {
        /* Note that the energies and virial are incorrect when not converged */
        if (fplog)
        {
            fprintf(fplog,
                    "step %s: EM did not converge in %d iterations, RMS force %6.2e\n",
                    gmx_step_str(mdstep, sbuf),
                    number_steps,
                    df[Min]);
        }
        fprintf(stderr,
                "step %s: EM did not converge in %d iterations, RMS force %6.2e\n",
                gmx_step_str(mdstep, sbuf),
                number_steps,
                df[Min]);
    }

    /* Copy back the coordinates and the forces */
    std::copy(pos[Min].begin(), pos[Min].end(), x.data());
    std::copy(force[Min].begin(), force[Min].end(), f->force().begin());
}

void done_shellfc(FILE* fplog, gmx_shellfc_t* shfc, int64_t numSteps)
{
    if (shfc && fplog && numSteps > 0)
    {
        double numStepsAsDouble = static_cast<double>(numSteps);
        fprintf(fplog,
                "Fraction of iterations that converged:           %.2f %%\n",
                (shfc->numConvergedIterations * 100.0) / numStepsAsDouble);
        fprintf(fplog,
                "Average number of force evaluations per MD step: %.2f\n\n",
                shfc->numForceEvaluations / numStepsAsDouble);
    }

    delete shfc;
}
