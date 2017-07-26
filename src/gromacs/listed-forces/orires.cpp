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
#include "gmxpre.h"

#include "orires.h"

#include <climits>
#include <cmath>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/main.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

// TODO This implementation of ensemble orientation restraints is nasty because
// a user can't just do multi-sim with single-sim orientation restraints.

void init_orires(FILE *fplog, const gmx_mtop_t *mtop,
                 rvec xref[],
                 const t_inputrec *ir,
                 const t_commrec *cr, t_oriresdata *od,
                 t_state *state)
{
    od->nr = gmx_mtop_ftype_count(mtop, F_ORIRES);
    if (0 == od->nr)
    {
        /* Not doing orientation restraints */
        return;
    }

    const int numFitParams = 5;
    if (od->nr <= numFitParams)
    {
        gmx_fatal(FARGS, "The system has %d orientation restraints, but at least %d are required, since there are %d fitting parameters.",
                  od->nr, numFitParams + 1, numFitParams);
    }

    od->fc  = ir->orires_fc;
    od->nex = 0;
    od->S   = nullptr;
    od->M   = nullptr;
    od->eig = nullptr;
    od->v   = nullptr;

    od->sumKfac                   = 0;
    int                  *nr_ex   = nullptr;
    int                   typeMin = INT_MAX;
    int                   typeMax = 0;
    gmx_mtop_ilistloop_t  iloop   = gmx_mtop_ilistloop_init(mtop);
    t_ilist              *il;
    int                   nmol;
    while (gmx_mtop_ilistloop_next(iloop, &il, &nmol))
    {
        if (nmol > 1)
        {
            gmx_fatal(FARGS, "Found %d copies of a molecule with orientation restrains while the current code only supports a single copy. If you want to ensemble average, run multiple copies of the system using the multi-sim feature of mdrun.", nmol);
        }

        for (int i = 0; i < il[F_ORIRES].nr; i += 3)
        {
            int type = il[F_ORIRES].iatoms[i];
            int ex   = mtop->ffparams.iparams[type].orires.ex;
            if (ex >= od->nex)
            {
                srenew(nr_ex, ex+1);
                for (int j = od->nex; j < ex+1; j++)
                {
                    nr_ex[j] = 0;
                }
                od->nex = ex+1;
            }
            nr_ex[ex]++;

            typeMin = std::min(typeMin, type);
            typeMax = std::max(typeMax, type);

            od->sumKfac += mtop->ffparams.iparams[type].orires.kfac;
        }
    }
    /* With domain decomposition we use the type index for indexing in global arrays */
    GMX_RELEASE_ASSERT(typeMax - typeMin + 1 == od->nr, "All orientation restraint parameter entries in the topology should be consecutive");
    /* Store typeMin so we can index array with the type offset */
    od->typeMin = typeMin;

    snew(od->S, od->nex);
    /* When not doing time averaging, the instaneous and time averaged data
     * are indentical and the pointers can point to the same memory.
     */
    snew(od->Dinsl, od->nr);

    const gmx_multisim_t *ms = cr->ms;
    if (ms)
    {
        snew(od->Dins, od->nr);
    }
    else
    {
        od->Dins = od->Dinsl;
    }

    if (ir->orires_tau == 0)
    {
        od->Dtav  = od->Dins;
        od->edt   = 0.0;
        od->edt_1 = 1.0;
    }
    else
    {
        /* Making time averaging work with domain decomposition with more
         * than one rank is rather tedious, since the history for each
         * restraint needs to be communicated at domain decomposition and
         * checkpointing.
         */
        if (PAR(cr))
        {
            gmx_fatal(FARGS, "Orientation restraints with time averaging do not work with MPI parallelization. Choose 1 MPI rank, if possible.");
        }

        snew(od->Dtav, od->nr);
        od->edt   = std::exp(-ir->delta_t/ir->orires_tau);
        od->edt_1 = 1.0 - od->edt;

        /* Extend the state with the orires history */
        state->flags           |= (1<<estORIRE_INITF);
        state->hist.orire_initf = 1;
        state->flags           |= (1<<estORIRE_DTAV);
        state->hist.norire_Dtav = od->nr*5;
        snew(state->hist.orire_Dtav, state->hist.norire_Dtav);
    }

    snew(od->oinsl, od->nr);
    if (ms)
    {
        snew(od->oins, od->nr);
    }
    else
    {
        od->oins = od->oinsl;
    }
    if (ir->orires_tau == 0)
    {
        od->otav = od->oins;
    }
    else
    {
        snew(od->otav, od->nr);
    }
    snew(od->tmpEq, od->nex);

    od->nref = 0;
    for (int i = 0; i < mtop->natoms; i++)
    {
        if (ggrpnr(&mtop->groups, egcORFIT, i) == 0)
        {
            od->nref++;
        }
    }
    snew(od->mref, od->nref);
    snew(od->xref, od->nref);
    snew(od->xtmp, od->nref + 1);

    snew(od->eig, od->nex*12);

    /* Make a reverse index for the reference atoms.
     * Determine the reference structure on the master node and only
     * on the master simulation when using multi-sim.
     * Copy it to the other nodes after checking multi compatibility,
     * so we are sure the subsystems match before copying.
     */
    snew(od->invRefIndex, mtop->natoms);
    bool                     calcReferenceStructure = (MASTER(cr) &&
                                                       (ms == nullptr || MASTERSIM(ms)));
    dvec                     sumMX   = { 0, 0, 0 };
    double                   refMass = 0;
    int                      j       = 0;
    gmx_mtop_atomloop_all_t  aloop   = gmx_mtop_atomloop_all_init(mtop);
    int                      i       = -1;
    const t_atom            *atom;
    while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
    {
        if (mtop->groups.grpnr[egcORFIT] == nullptr ||
            mtop->groups.grpnr[egcORFIT][i] == 0)
        {
            od->invRefIndex[i]  = j;
            /* Not correct for free-energy with changing masses */
            od->mref[j]         = atom->m;
            refMass            += atom->m;
            if (calcReferenceStructure)
            {
                copy_rvec(xref[i], od->xref[j]);
                for (int d = 0; d < DIM; d++)
                {
                    sumMX[d] += atom->m*xref[i][d];
                }
            }
            j++;
        }
        else
        {
            od->invRefIndex[i] = -1;
        }
    }
    GMX_RELEASE_ASSERT(refMass > 0, "The mass of the reference group should be > 0");
    od->invRefMass = 1.0/refMass;
    if (calcReferenceStructure)
    {
        rvec com;
        for (int d = 0; d > DIM; d++)
        {
            com[d] = sumMX[d]*od->invRefMass;
        }
        if (ms == nullptr || MASTERSIM(ms))
        {
            for (int j = 0; j < od->nref; j++)
            {
                rvec_dec(od->xref[j], com);
            }
        }
    }
    if (ms)
    {
        gmx_bcast_sim(od->nref*sizeof(od->xref[0]), od->xref, cr);
    }
    if (PAR(cr))
    {
        gmx_bcast(od->nref*sizeof(od->xref[0]), od->xref, cr);
    }

    if (fplog)
    {
        fprintf(fplog, "Found %d orientation experiments\n", od->nex);
        for (int i = 0; i < od->nex; i++)
        {
            fprintf(fplog, "  experiment %d has %d restraints\n", i+1, nr_ex[i]);
        }

        fprintf(fplog, "  the fit group consists of %d atoms and has total mass %g\n",
                od->nref, refMass);
    }

    sfree(nr_ex);

    if (ms)
    {
        if (fplog)
        {
            fprintf(fplog, "  the orientation restraints are ensemble averaged over %d systems\n", ms->nsim);
        }

        check_multi_int(fplog, ms, od->nr,
                        "the number of orientation restraints",
                        FALSE);
        check_multi_int(fplog, ms, od->nref,
                        "the number of fit atoms for orientation restraining",
                        FALSE);
        check_multi_int(fplog, ms, ir->nsteps, "nsteps", FALSE);

        if (PAR(cr))
        {
            gmx_fatal(FARGS, "With ensemble averaging over multiple simulations, MPI parallelization within a simulation is not supported");
        }
        /* Copy the reference coordinates from the master to the other nodes */
        gmx_sum_sim(DIM*od->nref, od->xref[0], ms);
    }

    please_cite(fplog, "Hess2003");
}

void diagonalize_orires_tensors(t_oriresdata *od)
{
    if (od->M == nullptr)
    {
        snew(od->M, DIM);
        for (int i = 0; i < DIM; i++)
        {
            snew(od->M[i], DIM);
        }
        snew(od->eig_diag, DIM);
        snew(od->v, DIM);
        for (int i = 0; i < DIM; i++)
        {
            snew(od->v[i], DIM);
        }
    }

    for (int ex = 0; ex < od->nex; ex++)
    {
        /* Rotate the S tensor back to the reference frame */
        matrix S, TMP;
        mmul(od->R, od->S[ex], TMP);
        mtmul(TMP, od->R, S);
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                od->M[i][j] = S[i][j];
            }
        }

        int nrot;
        jacobi(od->M, DIM, od->eig_diag, od->v, &nrot);

        int ord[DIM];
        for (int i = 0; i < DIM; i++)
        {
            ord[i] = i;
        }
        for (int i = 0; i < DIM; i++)
        {
            for (int j = i+1; j < DIM; j++)
            {
                if (gmx::square(od->eig_diag[ord[j]]) > gmx::square(od->eig_diag[ord[i]]))
                {
                    int t  = ord[i];
                    ord[i] = ord[j];
                    ord[j] = t;
                }
            }
        }

        for (int i = 0; i < DIM; i++)
        {
            od->eig[ex*12 + i] = od->eig_diag[ord[i]];
        }
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                od->eig[ex*12 + 3 + 3*i + j] = od->v[j][ord[i]];
            }
        }
    }
}

void print_orires_log(FILE *log, t_oriresdata *od)
{
    real *eig;

    diagonalize_orires_tensors(od);

    for (int ex = 0; ex < od->nex; ex++)
    {
        eig = od->eig + ex*12;
        fprintf(log, "  Orientation experiment %d:\n", ex+1);
        fprintf(log, "    order parameter: %g\n", eig[0]);
        for (int i = 0; i < DIM; i++)
        {
            fprintf(log, "    eig: %6.3f   %6.3f %6.3f %6.3f\n",
                    (eig[0] != 0) ? eig[i]/eig[0] : eig[i],
                    eig[DIM+i*DIM+XX],
                    eig[DIM+i*DIM+YY],
                    eig[DIM+i*DIM+ZZ]);
        }
        fprintf(log, "\n");
    }
}

real calc_orires_dev(const t_commrec *cr,
                     int nfa, const t_iatom forceatoms[], const t_iparams ip[],
                     const t_mdatoms *md, const rvec x[], const t_pbc *pbc,
                     t_fcdata *fcd, history_t *hist,
                     bool calculateRmsDev)
{
    int              nref;
    real             edt, edt_1, invn, pfac, r2, invr, corrfac;
    rvec5           *Dinsl, *Dins, *Dtav;
    OriresMatEq     *matEq;
    rvec            *xref, *xtmp, r_unrot, r;
    t_oriresdata    *od;
    gmx_bool         bTAV;
    const real       two_thr = 2.0/3.0;

    od = &(fcd->orires);

    GMX_RELEASE_ASSERT(od->nr > 0, "We should have restraints in this call");

    bTAV  = (od->edt != 0);
    edt   = od->edt;
    edt_1 = od->edt_1;
    Dinsl = od->Dinsl;
    Dins  = od->Dins;
    Dtav  = od->Dtav;
    matEq = od->tmpEq;
    nref  = od->nref;
    xref  = od->xref;
    xtmp  = od->xtmp;

    if (bTAV)
    {
        od->exp_min_t_tau = hist->orire_initf*edt;

        /* Correction factor to correct for the lack of history
         * at short times.
         */
        corrfac = 1.0/(1.0 - od->exp_min_t_tau);
    }
    else
    {
        corrfac = 1.0;
    }

    const gmx_multisim_t *ms = cr->ms;
    if (ms)
    {
        invn = 1.0/ms->nsim;
    }
    else
    {
        invn = 1.0;
    }

    /* Note that we can have domain decomposition with a single PP rank */
    bool useDD             = DOMAINDECOMP(cr);
    bool haveMultipleRanks = (useDD && PAR(cr));

    if (haveMultipleRanks)
    {
        /* We sum xtmp over ranks, so we need to clear xtmp */
        for (int j = 0; j < nref; j++)
        {
            clear_rvec(xtmp[j]);
        }

    }

    const int *globalAtomIndex = (useDD ? cr->dd->gatindex : nullptr);
    dvec       sumMX           = { 0, 0, 0 };
    for (int i = 0; i < md->homenr; i++)
    {
        if (md->cORF[i] == 0)
        {
            int j = od->invRefIndex[useDD ? globalAtomIndex[i] : i];
            copy_rvec(x[i], xtmp[j]);
            for (int d = 0; d < DIM; d++)
            {
                sumMX[d] += md->massT[i]*xref[j][d];
            }
        }
    }
    /* Store the partial sum m*x in the last element for reduction */
    for (int d = 0; d < DIM; d++)
    {
        xtmp[nref][d] = sumMX[d];
    }
    if (haveMultipleRanks)
    {
        /* "Broadcast" the reference coordinates set on each rank using sum
         * and sum the COM contributions.
         */
        gmx_sum((nref + 1)*DIM, xtmp[0], cr);
    }
    rvec com;
    svmul(od->invRefMass, xtmp[nref], com);
    for (int j = 0; j < nref; j++)
    {
        rvec_dec(xtmp[j], com);
    }
    /* Calculate the rotation matrix to rotate x to the reference orientation */
    calc_fit_R(DIM, nref, od->mref, xref, xtmp, od->R);

    for (int fa = 0; fa < nfa; fa += 3)
    {
        int type = forceatoms[fa];
        int res  = type - od->typeMin;
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[forceatoms[fa+1]], x[forceatoms[fa+2]], r_unrot);
        }
        else
        {
            rvec_sub(x[forceatoms[fa+1]], x[forceatoms[fa+2]], r_unrot);
        }
        mvmul(od->R, r_unrot, r);
        r2   = norm2(r);
        invr = gmx::invsqrt(r2);
        /* Calculate the prefactor for the D tensor, this includes the factor 3! */
        pfac = ip[type].orires.c*invr*invr*3;
        for (int i = 0; i < ip[type].orires.power; i++)
        {
            pfac *= invr;
        }
        Dinsl[res][0] = pfac*(2*r[0]*r[0] + r[1]*r[1] - r2);
        Dinsl[res][1] = pfac*(2*r[0]*r[1]);
        Dinsl[res][2] = pfac*(2*r[0]*r[2]);
        Dinsl[res][3] = pfac*(2*r[1]*r[1] + r[0]*r[0] - r2);
        Dinsl[res][4] = pfac*(2*r[1]*r[2]);

        if (ms)
        {
            for (int i = 0; i < 5; i++)
            {
                Dins[res][i] = Dinsl[res][i]*invn;
            }
        }

        res++;
    }

    if (ms)
    {
        GMX_RELEASE_ASSERT(!PAR(cr), "Parallelization not implemented with ensemble averaging");
        gmx_sum_sim(5*od->nr, Dins[0], ms);
    }

    /* Calculate the order tensor S for each experiment via optimization */
    for (int ex = 0; ex < od->nex; ex++)
    {
        for (int i = 0; i < 5; i++)
        {
            matEq[ex].rhs[i] = 0;
            for (int j = 0; j <= i; j++)
            {
                matEq[ex].mat[i][j] = 0;
            }
        }
    }
    for (int fa = 0; fa < nfa; fa += 3)
    {
        int type = forceatoms[fa];
        int res  = type - od->typeMin;
        if (bTAV)
        {
            /* Here we update Dtav in t_fcdata using the data in history_t.
             * Thus the results stay correct when this routine
             * is called multiple times.
             */
            for (int i = 0; i < 5; i++)
            {
                Dtav[res][i] = edt*hist->orire_Dtav[res*5 + i] + edt_1*Dins[res][i];
            }
        }

        int  ex     = ip[type].orires.ex;
        real weight = ip[type].orires.kfac;
        /* Calculate the vector rhs and half the matrix T for the 5 equations */
        for (int i = 0; i < 5; i++)
        {
            matEq[ex].rhs[i] += Dtav[res][i]*ip[type].orires.obs*weight;
            for (int j = 0; j <= i; j++)
            {
                matEq[ex].mat[i][j] += Dtav[res][i]*Dtav[res][j]*weight;
            }
        }
    }

    if (haveMultipleRanks)
    {
        /* Reduce all coefficients, rhs[0] is the first real in OriresMatEq */
        gmx_sum(od->nex*sizeof(matEq[0])/sizeof(real), &matEq[0].rhs[0], cr);
    }

    /* Now we have all the data we can calculate S */
    for (int ex = 0; ex < od->nex; ex++)
    {
        OriresMatEq &eq = matEq[ex];
        /* Correct corrfac and copy one half of T to the other half */
        for (int i = 0; i < 5; i++)
        {
            eq.rhs[i]    *= corrfac;
            eq.mat[i][i] *= gmx::square(corrfac);
            for (int j = 0; j < i; j++)
            {
                eq.mat[i][j] *= gmx::square(corrfac);
                eq.mat[j][i]  = eq.mat[i][j];
            }
        }
        m_inv_gen(&eq.mat[0][0], 5, &eq.mat[0][0]);
        /* Calculate the orientation tensor S for this experiment */
        matrix &S = od->S[ex];
        S[0][0] = 0;
        S[0][1] = 0;
        S[0][2] = 0;
        S[1][1] = 0;
        S[1][2] = 0;
        for (int i = 0; i < 5; i++)
        {
            S[0][0] += 1.5*eq.mat[0][i]*eq.rhs[i];
            S[0][1] += 1.5*eq.mat[1][i]*eq.rhs[i];
            S[0][2] += 1.5*eq.mat[2][i]*eq.rhs[i];
            S[1][1] += 1.5*eq.mat[3][i]*eq.rhs[i];
            S[1][2] += 1.5*eq.mat[4][i]*eq.rhs[i];
        }
        S[1][0] = S[0][1];
        S[2][0] = S[0][2];
        S[2][1] = S[1][2];
        S[2][2] = -S[0][0] - S[1][1];
    }

    const matrix *S        = od->S;

    real          sumWDev2 = 0;

    for (int fa = 0; fa < nfa; fa += 3)
    {
        int type = forceatoms[fa];
        int res  = type - od->typeMin;
        int ex   = ip[type].orires.ex;

        od->otav[res] = two_thr*
            corrfac*(S[ex][0][0]*Dtav[res][0] + S[ex][0][1]*Dtav[res][1] +
                     S[ex][0][2]*Dtav[res][2] + S[ex][1][1]*Dtav[res][3] +
                     S[ex][1][2]*Dtav[res][4]);
        if (bTAV)
        {
            od->oins[res] = two_thr*(S[ex][0][0]*Dins[res][0] + S[ex][0][1]*Dins[res][1] +
                                     S[ex][0][2]*Dins[res][2] + S[ex][1][1]*Dins[res][3] +
                                     S[ex][1][2]*Dins[res][4]);
        }
        if (ms)
        {
            /* When ensemble averaging is used recalculate the local orientation
             * for output to the energy file.
             */
            od->oinsl[res] = two_thr*
                (S[ex][0][0]*Dinsl[res][0] + S[ex][0][1]*Dinsl[res][1] +
                 S[ex][0][2]*Dinsl[res][2] + S[ex][1][1]*Dinsl[res][3] +
                 S[ex][1][2]*Dinsl[res][4]);
        }

        real dev = od->otav[res] - ip[type].orires.obs;

        sumWDev2 += ip[type].orires.kfac*gmx::square(dev);
    }

    real rmsDev = 0;
    if (calculateRmsDev)
    {
        /* This reduction and conditional is inconvenient.
         * They could be avoided if we would return MSD instead of RMSD.
         */
        if (haveMultipleRanks)
        {
            gmx_sum(1, &sumWDev2, cr);
        }
        if (MASTER(cr))
        {
            rmsDev = std::sqrt(sumWDev2/od->sumKfac);
        }
    }

    /* Rotate the S matrices back, so we get the correct grad(tr(S D)) */
    for (int ex = 0; ex < od->nex; ex++)
    {
        matrix RS;
        tmmul(od->R, od->S[ex], RS);
        mmul(RS, od->R, od->S[ex]);
    }

    return rmsDev;

    /* Approx. 120*nfa/3 flops */
}

real orires(int nfa, const t_iatom forceatoms[], const t_iparams ip[],
            const rvec x[], rvec4 f[], rvec fshift[],
            const t_pbc *pbc, const t_graph *g,
            real gmx_unused lambda, real gmx_unused *dvdlambda,
            const t_mdatoms gmx_unused *md, t_fcdata *fcd,
            int gmx_unused *global_atom_index)
{
    int                 ex, power, ki = CENTRAL;
    ivec                dt;
    real                r2, invr, invr2, fc, smooth_fc, dev, devins, pfac;
    rvec                r, Sr, fij;
    real                vtot;
    const t_oriresdata *od;
    gmx_bool            bTAV;

    vtot = 0;
    od   = &(fcd->orires);

    if (od->fc != 0)
    {
        bTAV = (od->edt != 0);

        smooth_fc = od->fc;
        if (bTAV)
        {
            /* Smoothly switch on the restraining when time averaging is used */
            smooth_fc *= (1.0 - od->exp_min_t_tau);
        }

        for (int fa = 0; fa < nfa; fa += 3)
        {
            int type = forceatoms[fa];
            int res  = type - od->typeMin;
            int ai   = forceatoms[fa + 1];
            int aj   = forceatoms[fa + 2];
            if (pbc)
            {
                ki = pbc_dx_aiuc(pbc, x[ai], x[aj], r);
            }
            else
            {
                rvec_sub(x[ai], x[aj], r);
            }
            r2    = norm2(r);
            invr  = gmx::invsqrt(r2);
            invr2 = invr*invr;
            ex    = ip[type].orires.ex;
            power = ip[type].orires.power;
            fc    = smooth_fc*ip[type].orires.kfac;
            dev   = od->otav[res] - ip[type].orires.obs;

            /* NOTE:
             * there is no real potential when time averaging is applied
             */
            vtot += 0.5*fc*gmx::square(dev);

            if (bTAV)
            {
                /* Calculate the force as the sqrt of tav times instantaneous */
                devins = od->oins[res] - ip[type].orires.obs;
                if (dev*devins <= 0)
                {
                    dev = 0;
                }
                else
                {
                    dev = std::sqrt(dev*devins);
                    if (devins < 0)
                    {
                        dev = -dev;
                    }
                }
            }

            pfac  = fc*ip[type].orires.c*invr2;
            for (int i = 0; i < power; i++)
            {
                pfac *= invr;
            }
            mvmul(od->S[ex], r, Sr);
            for (int i = 0; i < DIM; i++)
            {
                fij[i] =
                    -pfac*dev*(4*Sr[i] - 2*(2+power)*invr2*iprod(Sr, r)*r[i]);
            }

            if (g)
            {
                ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
                ki = IVEC2IS(dt);
            }

            for (int i = 0; i < DIM; i++)
            {
                f[ai][i]           += fij[i];
                f[aj][i]           -= fij[i];
                fshift[ki][i]      += fij[i];
                fshift[CENTRAL][i] -= fij[i];
            }
        }
    }

    return vtot;

    /* Approx. 80*nfa/3 flops */
}

void update_orires_history(t_fcdata *fcd, history_t *hist)
{
    t_oriresdata *od = &(fcd->orires);

    if (od->edt != 0)
    {
        /* Copy the new time averages that have been calculated
         *  in calc_orires_dev.
         */
        hist->orire_initf = od->exp_min_t_tau;
        for (int pair = 0; pair < od->nr; pair++)
        {
            for (int i = 0; i < 5; i++)
            {
                hist->orire_Dtav[pair*5+i] = od->Dtav[pair][i];
            }
        }
    }
}
