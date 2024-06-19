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

#include "orires.h"

#include <climits>
#include <cmath>

#include <algorithm>
#include <array>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/nrjac.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

using gmx::ArrayRef;
using gmx::RVec;

// TODO This implementation of ensemble orientation restraints is nasty because
// a user can't just do multi-sim with single-sim orientation restraints.

void extendStateWithOriresHistory(const gmx_mtop_t& mtop, const t_inputrec& ir, t_state* globalState)
{
    GMX_RELEASE_ASSERT(globalState != nullptr,
                       "We need a valid global state in extendStateWithOriresHistory()");

    const int numRestraints = gmx_mtop_ftype_count(mtop, F_ORIRES);
    if (numRestraints > 0 && ir.orires_tau > 0)
    {
        /* Extend the state with the orires history */
        globalState->addEntry(StateEntry::OrireInitF);
        globalState->hist.orire_initf = 1;
        globalState->addEntry(StateEntry::OrireDtav);
        globalState->hist.orire_Dtav.resize(numRestraints * 5);
    }
}

namespace
{

//! Creates and returns a list of global atom indices of the orientation restraint fit group
std::vector<gmx::Index> fitGlobalAtomIndices(const gmx_mtop_t& mtop)
{
    std::vector<gmx::Index> indices;

    for (int i = 0; i < mtop.natoms; i++)
    {
        if (getGroupType(mtop.groups, SimulationAtomGroupType::OrientationRestraintsFit, i) == 0)
        {
            indices.push_back(i);
        }
    }

    return indices;
}

} // namespace

t_oriresdata::t_oriresdata(FILE*                     fplog,
                           const gmx_mtop_t&         mtop,
                           const t_inputrec&         ir,
                           const gmx_multisim_t*     ms,
                           t_state*                  globalState,
                           gmx::LocalAtomSetManager* localAtomSetManager) :
    numRestraints(gmx_mtop_ftype_count(mtop, F_ORIRES)),
    fitLocalAtomSet_(localAtomSetManager->add(fitGlobalAtomIndices(mtop)))
{
    GMX_RELEASE_ASSERT(numRestraints > 0,
                       "orires() should only be called with orientation restraints present");

    const int numFitParams = 5;
    if (numRestraints <= numFitParams)
    {
        const std::string mesg = gmx::formatString(
                "The system has %d orientation restraints, but at least %d are required, since "
                "there are %d fitting parameters.",
                numRestraints,
                numFitParams + 1,
                numFitParams);
        GMX_THROW(gmx::InvalidInputError(mesg));
    }

    if (ir.bPeriodicMols)
    {
        /* Since we apply fitting, we need to make molecules whole and this
         * can not be done when periodic molecules are present.
         */
        GMX_THROW(gmx::InvalidInputError(
                "Orientation restraints can not be applied when periodic molecules are present "
                "in the system"));
    }

    GMX_RELEASE_ASSERT(globalState != nullptr, "We need a valid global state in t_oriresdata()");

    fc             = ir.orires_fc;
    numExperiments = 0;

    std::vector<int> nr_ex;
    typeMin     = INT_MAX;
    int typeMax = 0;
    for (const auto il : IListRange(mtop))
    {
        const int numOrires = il.list()[F_ORIRES].size();
        if (il.nmol() > 1 && numOrires > 0)
        {
            const std::string mesg = gmx::formatString(
                    "Found %d copies of a molecule with orientation restrains while the current "
                    "code only supports a single copy. If you want to ensemble average, run "
                    "multiple copies of the system using the multi-sim feature of mdrun.",
                    il.nmol());
            GMX_THROW(gmx::InvalidInputError(mesg));
        }

        for (int i = 0; i < numOrires; i += 3)
        {
            int type = il.list()[F_ORIRES].iatoms[i];
            int ex   = mtop.ffparams.iparams[type].orires.ex;
            if (ex >= numExperiments)
            {
                nr_ex.resize(ex + 1, 0);
                numExperiments = ex + 1;
            }
            nr_ex[ex]++;

            typeMin = std::min(typeMin, type);
            typeMax = std::max(typeMax, type);
        }
    }
    /* With domain decomposition we use the type index for indexing in global arrays */
    GMX_RELEASE_ASSERT(
            typeMax - typeMin + 1 == numRestraints,
            "All orientation restraint parameter entries in the topology should be consecutive");

    snew(orderTensors, numExperiments);
    /* When not doing time averaging, the instaneous and time averaged data
     * are identical and the pointers can point to the same memory.
     */
    snew(DTensors, numRestraints);

    if (ms)
    {
        snew(DTensorsEnsembleAv, numRestraints);
    }
    else
    {
        DTensorsEnsembleAv = DTensors;
    }

    if (ir.orires_tau == 0)
    {
        DTensorsTimeAndEnsembleAv = DTensorsEnsembleAv;
        edt                       = 0.0;
        edt_1                     = 1.0;
    }
    else
    {
        snew(DTensorsTimeAndEnsembleAv, numRestraints);
        edt   = std::exp(-ir.delta_t / ir.orires_tau);
        edt_1 = 1.0 - edt;

        timeAveragingInitFactor_     = std::reference_wrapper<real>(globalState->hist.orire_initf);
        DTensorsTimeAveragedHistory_ = globalState->hist.orire_Dtav;
    }

    orientations.resize(numRestraints);
    if (ms)
    {
        orientationsEnsembleAvBuffer.resize(numRestraints);
        orientationsEnsembleAv = orientationsEnsembleAvBuffer;
    }
    else
    {
        orientationsEnsembleAv = orientations;
    }
    if (ir.orires_tau == 0)
    {
        orientationsTimeAndEnsembleAv = orientationsEnsembleAv;
    }
    else
    {
        orientationsTimeAndEnsembleAvBuffer.resize(numRestraints);
        orientationsTimeAndEnsembleAv = orientationsTimeAndEnsembleAvBuffer;
    }
    tmpEq.resize(numExperiments);

    eigenOutput.resize(numExperiments * c_numEigenRealsPerExperiment);

    /* Determine the reference structure on the main node.
     * Copy it to the other nodes after checking multi compatibility,
     * so we are sure the subsystems match before copying.
     */
    auto   x    = makeConstArrayRef(globalState->x);
    rvec   com  = { 0, 0, 0 };
    double mtot = 0.0;
    for (const AtomProxy atomP : AtomRange(mtop))
    {
        const t_atom& local = atomP.atom();
        int           i     = atomP.globalAtomNumber();
        if (getGroupType(mtop.groups, SimulationAtomGroupType::OrientationRestraintsFit, i) == 0)
        {
            // Not correct for free-energy with changing masses
            const real mass = local.m;
            // Note that only one rank per sim is supported.
            if (isMainSim(ms))
            {
                referenceCoordinates_.push_back(x[i]);
                for (int d = 0; d < DIM; d++)
                {
                    com[d] += mass * x[i][d];
                }
            }
            fitMasses_.push_back(mass);
            mtot += mass;
        }
    }

    svmul(1.0 / mtot, com, com);
    if (isMainSim(ms))
    {
        for (auto& refCoord : referenceCoordinates_)
        {
            refCoord -= com;
        }
    }

    const size_t numFitAtoms = referenceCoordinates_.size();
    xTmp_.resize(numFitAtoms);

    if (fplog)
    {
        fprintf(fplog, "Found %d orientation experiments\n", numExperiments);
        for (int i = 0; i < numExperiments; i++)
        {
            fprintf(fplog, "  experiment %d has %d restraints\n", i + 1, nr_ex[i]);
        }

        fprintf(fplog, "  the fit group consists of %zu atoms and has total mass %g\n", numFitAtoms, mtot);
    }

    if (ms)
    {
        if (fplog)
        {
            fprintf(fplog,
                    "  the orientation restraints are ensemble averaged over %d systems\n",
                    ms->numSimulations_);
        }

        check_multi_int(fplog, ms, numRestraints, "the number of orientation restraints", FALSE);
        check_multi_int(
                fplog, ms, numFitAtoms, "the number of fit atoms for orientation restraining", FALSE);
        check_multi_int(fplog, ms, ir.nsteps, "nsteps", FALSE);
        /* Copy the reference coordinates from the main to the other nodes */
        gmx_sum_sim(DIM * referenceCoordinates_.size(), as_rvec_array(referenceCoordinates_.data())[0], ms);
    }

    please_cite(fplog, "Hess2003");
}

t_oriresdata::~t_oriresdata()
{
    sfree(orderTensors);
    if (DTensorsTimeAndEnsembleAv != DTensorsEnsembleAv)
    {
        sfree(DTensorsTimeAndEnsembleAv);
    }
    if (DTensorsEnsembleAv != DTensors)
    {
        sfree(DTensorsEnsembleAv);
    }
    sfree(DTensors);
}

void diagonalize_orires_tensors(t_oriresdata* od)
{
    for (int ex = 0; ex < od->numExperiments; ex++)
    {
        /* Rotate the S tensor back to the reference frame */
        matrix S, TMP;
        mmul(od->rotationMatrix, od->orderTensors[ex], TMP);
        mtmul(TMP, od->rotationMatrix, S);
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                od->M[i][j] = S[i][j];
            }
        }

        jacobi(od->M, od->eig_diag, od->v);

        int ord[DIM];
        for (int i = 0; i < DIM; i++)
        {
            ord[i] = i;
        }
        for (int i = 0; i < DIM; i++)
        {
            for (int j = i + 1; j < DIM; j++)
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
            od->eigenOutput[ex * t_oriresdata::c_numEigenRealsPerExperiment + i] = od->eig_diag[ord[i]];
        }
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                od->eigenOutput[ex * t_oriresdata::c_numEigenRealsPerExperiment + 3 + 3 * i + j] =
                        od->v[j][ord[i]];
            }
        }
    }
}

void print_orires_log(FILE* log, t_oriresdata* od)
{
    diagonalize_orires_tensors(od);

    for (int ex = 0; ex < od->numExperiments; ex++)
    {
        const real* eig = od->eigenOutput.data() + ex * t_oriresdata::c_numEigenRealsPerExperiment;
        fprintf(log, "  Orientation experiment %d:\n", ex + 1);
        fprintf(log, "    order parameter: %g\n", eig[0]);
        for (int i = 0; i < DIM; i++)
        {
            fprintf(log,
                    "    eig: %6.3f   %6.3f %6.3f %6.3f\n",
                    (eig[0] != 0) ? eig[i] / eig[0] : eig[i],
                    eig[DIM + i * DIM + XX],
                    eig[DIM + i * DIM + YY],
                    eig[DIM + i * DIM + ZZ]);
        }
        fprintf(log, "\n");
    }
}

real calc_orires_dev(const gmx_multisim_t* ms,
                     int                   nfa,
                     const t_iatom         forceatoms[],
                     const t_iparams       ip[],
                     ArrayRef<const RVec>  xWholeMolecules,
                     const rvec            x[],
                     const t_pbc*          pbc,
                     t_oriresdata*         od)
{
    real       invn, pfac, r2, invr, corrfac, wsv2, sw, dev;
    rvec       com, r_unrot, r;
    const real two_thr = 2.0 / 3.0;

    const bool                 bTAV  = (od->edt != 0);
    const real                 edt   = od->edt;
    const real                 edt_1 = od->edt_1;
    gmx::ArrayRef<OriresMatEq> matEq = od->tmpEq;
    gmx::ArrayRef<gmx::RVec>   xFit  = od->xTmp();

    if (bTAV)
    {
        od->exp_min_t_tau = od->timeAveragingInitFactor() * edt;

        /* Correction factor to correct for the lack of history
         * at short times.
         */
        corrfac = 1.0 / (1.0 - od->exp_min_t_tau);
    }
    else
    {
        corrfac = 1.0;
    }

    if (ms)
    {
        invn = 1.0 / ms->numSimulations_;
    }
    else
    {
        invn = 1.0;
    }

    // Extract the local atom indices involved in the fit group
    const auto fitLocalAtomIndices = od->fitLocalAtomSet().localIndex();
    // We need all atoms in the fit group to be local available. This means that
    // orientation restraining is limited to one PP-rank. This should be ensured
    // by the mdrun setup code. We assert here to catch incorrect setup code.
    GMX_RELEASE_ASSERT(fitLocalAtomIndices.size() == od->referenceCoordinates().size(),
                       "All fit atoms should be locally available");

    clear_rvec(com);
    double     mtot               = 0.0;
    gmx::Index referenceListIndex = 0;
    for (const int fitLocalAtomIndex : fitLocalAtomIndices)
    {
        const gmx::RVec& xLocal  = xWholeMolecules[fitLocalAtomIndex];
        const real       mass    = od->fitMasses()[referenceListIndex];
        xFit[referenceListIndex] = xLocal;
        for (int d = 0; d < DIM; d++)
        {
            com[d] += mass * xLocal[d];
        }
        mtot += mass;
        referenceListIndex++;
    }
    svmul(1.0 / mtot, com, com);
    for (auto& xFitCoord : xFit)
    {
        xFitCoord -= com;
    }
    /* Calculate the rotation matrix to rotate x to the reference orientation */
    calc_fit_R(DIM,
               xFit.size(),
               od->fitMasses().data(),
               as_rvec_array(od->referenceCoordinates().data()),
               as_rvec_array(xFit.data()),
               od->rotationMatrix);

    for (int fa = 0; fa < nfa; fa += 3)
    {
        const int type           = forceatoms[fa];
        const int restraintIndex = type - od->typeMin;
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[forceatoms[fa + 1]], x[forceatoms[fa + 2]], r_unrot);
        }
        else
        {
            rvec_sub(x[forceatoms[fa + 1]], x[forceatoms[fa + 2]], r_unrot);
        }
        mvmul(od->rotationMatrix, r_unrot, r);
        r2   = norm2(r);
        invr = gmx::invsqrt(r2);
        /* Calculate the prefactor for the D tensor, this includes the factor 3! */
        pfac = ip[type].orires.c * invr * invr * 3;
        for (int i = 0; i < ip[type].orires.power; i++)
        {
            pfac *= invr;
        }
        rvec5& Dinsl = od->DTensors[restraintIndex];
        Dinsl[0]     = pfac * (2 * r[0] * r[0] + r[1] * r[1] - r2);
        Dinsl[1]     = pfac * (2 * r[0] * r[1]);
        Dinsl[2]     = pfac * (2 * r[0] * r[2]);
        Dinsl[3]     = pfac * (2 * r[1] * r[1] + r[0] * r[0] - r2);
        Dinsl[4]     = pfac * (2 * r[1] * r[2]);

        if (ms)
        {
            for (int i = 0; i < 5; i++)
            {
                od->DTensorsEnsembleAv[restraintIndex][i] = Dinsl[i] * invn;
            }
        }
    }

    if (ms)
    {
        gmx_sum_sim(5 * od->numRestraints, od->DTensorsEnsembleAv[0], ms);
    }

    /* Calculate the order tensor S for each experiment via optimization */
    for (int ex = 0; ex < od->numExperiments; ex++)
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
        const int type           = forceatoms[fa];
        const int restraintIndex = type - od->typeMin;
        rvec5&    Dtav           = od->DTensorsTimeAndEnsembleAv[restraintIndex];
        if (bTAV)
        {
            /* Here we update DTensorsTimeAndEnsembleAv in t_fcdata using the data in history_t.
             * Thus the results stay correct when this routine
             * is called multiple times.
             */
            for (int i = 0; i < 5; i++)
            {
                Dtav[i] = edt * od->DTensorsTimeAveragedHistory()[restraintIndex * 5 + i]
                          + edt_1 * od->DTensorsEnsembleAv[restraintIndex][i];
            }
        }

        int  ex     = ip[type].orires.ex;
        real weight = ip[type].orires.kfac;
        /* Calculate the vector rhs and half the matrix T for the 5 equations */
        for (int i = 0; i < 5; i++)
        {
            matEq[ex].rhs[i] += Dtav[i] * ip[type].orires.obs * weight;
            for (int j = 0; j <= i; j++)
            {
                matEq[ex].mat[i][j] += Dtav[i] * Dtav[j] * weight;
            }
        }
    }
    /* Now we have all the data we can calculate S */
    for (int ex = 0; ex < od->numExperiments; ex++)
    {
        OriresMatEq& eq = matEq[ex];
        /* Correct corrfac and copy one half of T to the other half */
        for (int i = 0; i < 5; i++)
        {
            eq.rhs[i] *= corrfac;
            eq.mat[i][i] *= gmx::square(corrfac);
            for (int j = 0; j < i; j++)
            {
                eq.mat[i][j] *= gmx::square(corrfac);
                eq.mat[j][i] = eq.mat[i][j];
            }
        }
        m_inv_gen(&eq.mat[0][0], 5, &eq.mat[0][0]);
        /* Calculate the orientation tensor S for this experiment */
        matrix& S = od->orderTensors[ex];
        S[0][0]   = 0;
        S[0][1]   = 0;
        S[0][2]   = 0;
        S[1][1]   = 0;
        S[1][2]   = 0;
        for (int i = 0; i < 5; i++)
        {
            S[0][0] += 1.5 * eq.mat[0][i] * eq.rhs[i];
            S[0][1] += 1.5 * eq.mat[1][i] * eq.rhs[i];
            S[0][2] += 1.5 * eq.mat[2][i] * eq.rhs[i];
            S[1][1] += 1.5 * eq.mat[3][i] * eq.rhs[i];
            S[1][2] += 1.5 * eq.mat[4][i] * eq.rhs[i];
        }
        S[1][0] = S[0][1];
        S[2][0] = S[0][2];
        S[2][1] = S[1][2];
        S[2][2] = -S[0][0] - S[1][1];
    }

    const matrix* S = od->orderTensors;

    wsv2 = 0;
    sw   = 0;

    for (int fa = 0; fa < nfa; fa += 3)
    {
        const int type           = forceatoms[fa];
        const int restraintIndex = type - od->typeMin;
        const int ex             = ip[type].orires.ex;

        const rvec5& Dtav = od->DTensorsTimeAndEnsembleAv[restraintIndex];
        od->orientationsTimeAndEnsembleAv[restraintIndex] =
                two_thr * corrfac
                * (S[ex][0][0] * Dtav[0] + S[ex][0][1] * Dtav[1] + S[ex][0][2] * Dtav[2]
                   + S[ex][1][1] * Dtav[3] + S[ex][1][2] * Dtav[4]);
        if (bTAV)
        {
            const rvec5& Dins = od->DTensorsEnsembleAv[restraintIndex];
            od->orientationsEnsembleAv[restraintIndex] =
                    two_thr
                    * (S[ex][0][0] * Dins[0] + S[ex][0][1] * Dins[1] + S[ex][0][2] * Dins[2]
                       + S[ex][1][1] * Dins[3] + S[ex][1][2] * Dins[4]);
        }
        if (ms)
        {
            /* When ensemble averaging is used recalculate the local orientation
             * for output to the energy file.
             */
            const rvec5& Dinsl = od->DTensors[restraintIndex];
            od->orientations[restraintIndex] =
                    two_thr
                    * (S[ex][0][0] * Dinsl[0] + S[ex][0][1] * Dinsl[1] + S[ex][0][2] * Dinsl[2]
                       + S[ex][1][1] * Dinsl[3] + S[ex][1][2] * Dinsl[4]);
        }

        dev = od->orientationsTimeAndEnsembleAv[restraintIndex] - ip[type].orires.obs;

        wsv2 += ip[type].orires.kfac * gmx::square(dev);
        sw += ip[type].orires.kfac;
    }
    od->rmsdev = std::sqrt(wsv2 / sw);

    /* Rotate the S matrices back, so we get the correct grad(tr(S D)) */
    for (int ex = 0; ex < od->numExperiments; ex++)
    {
        matrix RS;
        tmmul(od->rotationMatrix, od->orderTensors[ex], RS);
        mmul(RS, od->rotationMatrix, od->orderTensors[ex]);
    }

    return od->rmsdev;

    /* Approx. 120*nfa/3 flops */
}

real orires(int             nfa,
            const t_iatom   forceatoms[],
            const t_iparams ip[],
            const rvec      x[],
            rvec4           f[],
            rvec            fshift[],
            const t_pbc*    pbc,
            real gmx_unused lambda,
            real gmx_unused* dvdlambda,
            gmx::ArrayRef<const real> /*charge*/,
            t_fcdata gmx_unused* fcd,
            t_disresdata gmx_unused* disresdata,
            t_oriresdata*            oriresdata,
            int gmx_unused* global_atom_index)
{
    int      ex, power, ki = gmx::c_centralShiftIndex;
    real     r2, invr, invr2, fc, smooth_fc, dev, devins, pfac;
    rvec     r, Sr, fij;
    real     vtot;
    gmx_bool bTAV;

    vtot = 0;

    if (oriresdata->fc != 0)
    {
        bTAV = (oriresdata->edt != 0);

        smooth_fc = oriresdata->fc;
        if (bTAV)
        {
            /* Smoothly switch on the restraining when time averaging is used */
            smooth_fc *= (1.0 - oriresdata->exp_min_t_tau);
        }

        for (int fa = 0; fa < nfa; fa += 3)
        {
            const int type           = forceatoms[fa];
            const int ai             = forceatoms[fa + 1];
            const int aj             = forceatoms[fa + 2];
            const int restraintIndex = type - oriresdata->typeMin;
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
            invr2 = invr * invr;
            ex    = ip[type].orires.ex;
            power = ip[type].orires.power;
            fc    = smooth_fc * ip[type].orires.kfac;
            dev   = oriresdata->orientationsTimeAndEnsembleAv[restraintIndex] - ip[type].orires.obs;

            /* NOTE:
             * there is no real potential when time averaging is applied
             */
            vtot += 0.5 * fc * gmx::square(dev);

            if (bTAV)
            {
                /* Calculate the force as the sqrt of tav times instantaneous */
                devins = oriresdata->orientationsEnsembleAv[restraintIndex] - ip[type].orires.obs;
                if (dev * devins <= 0)
                {
                    dev = 0;
                }
                else
                {
                    dev = std::sqrt(dev * devins);
                    if (devins < 0)
                    {
                        dev = -dev;
                    }
                }
            }

            pfac = fc * ip[type].orires.c * invr2;
            for (int i = 0; i < power; i++)
            {
                pfac *= invr;
            }
            mvmul(oriresdata->orderTensors[ex], r, Sr);
            for (int i = 0; i < DIM; i++)
            {
                fij[i] = -pfac * dev * (4 * Sr[i] - 2 * (2 + power) * invr2 * iprod(Sr, r) * r[i]);
            }

            for (int i = 0; i < DIM; i++)
            {
                f[ai][i] += fij[i];
                f[aj][i] -= fij[i];
                if (fshift)
                {
                    fshift[ki][i] += fij[i];
                    fshift[gmx::c_centralShiftIndex][i] -= fij[i];
                }
            }
        }
    }

    return vtot;

    /* Approx. 80*nfa/3 flops */
}

void t_oriresdata::updateHistory()
{
    if (edt != 0)
    {
        /* Copy the new time averages that have been calculated
         *  in calc_orires_dev.
         */
        *timeAveragingInitFactor_ = exp_min_t_tau;
        for (int pair = 0; pair < numRestraints; pair++)
        {
            for (int i = 0; i < 5; i++)
            {
                DTensorsTimeAveragedHistory_[pair * 5 + i] = DTensorsTimeAndEnsembleAv[pair][i];
            }
        }
    }
}
