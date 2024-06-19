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

#include "gromacs/mdtypes/inputrec.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <bitset>
#include <filesystem>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "gromacs/applied_forces/awh/read_params.h"
#include "gromacs/math/veccompare.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/txtdump.h"

//! Macro to select a bool name
#define EBOOL(e) gmx::boolToString(e)

/* The minimum number of integration steps per decay time tau required
 * to stay reasonable close to the target value.
 */
const int c_minimumStepsPerTau = 5;

/* The minimum number of integration steps per period
 * for reasonably accurate integration of second order coupling algorithms.
 */
const int c_minimumStepsPerPeriod = 20;

/* Default values for T- and P- coupling intervals, used when the are no other
 * restrictions.
 */
const int c_defaultNstTCouple = 100;
const int c_defaultNstPCouple = 100;

t_inputrec::t_inputrec() :
    fepvals(std::make_unique<t_lambda>()),
    simtempvals(std::make_unique<t_simtemp>()),
    expandedvals(std::make_unique<t_expanded>())
{
}

t_inputrec::~t_inputrec()
{
    done_inputrec(this);
}

int tcouple_min_integration_steps(TemperatureCoupling etc)
{
    int n;

    switch (etc)
    {
        case TemperatureCoupling::No: n = 0; break;
        case TemperatureCoupling::Berendsen:
        case TemperatureCoupling::Yes:
        case TemperatureCoupling::VRescale:
            /* V-rescale supports instantaneous rescaling */
            n = c_minimumStepsPerTau;
            break;
        case TemperatureCoupling::NoseHoover: n = c_minimumStepsPerPeriod; break;
        case TemperatureCoupling::Andersen:
        case TemperatureCoupling::AndersenMassive: n = 1; break;
        default: gmx_incons("Unknown etc value");
    }

    return n;
}

int ir_optimal_nsttcouple(const t_inputrec* ir)
{
    int  nmin, nwanted, n;
    real tau_min;
    int  g;

    nmin = tcouple_min_integration_steps(ir->etc);

    nwanted = c_defaultNstTCouple;

    tau_min = 1e20;
    if (ir->etc != TemperatureCoupling::No)
    {
        for (g = 0; g < ir->opts.ngtc; g++)
        {
            if (ir->opts.tau_t[g] > 0)
            {
                tau_min = std::min(tau_min, ir->opts.tau_t[g]);
            }
        }
    }

    if (nmin == 0 || ir->delta_t * nwanted <= tau_min / nmin)
    {
        n = nwanted;
    }
    else
    {
        n = static_cast<int>(tau_min / (ir->delta_t * nmin) + 0.001);
        if (n < 1)
        {
            n = 1;
        }
        while (nwanted % n != 0)
        {
            n--;
        }
    }

    return n;
}

int pcouple_min_integration_steps(PressureCoupling epc)
{
    int n;

    switch (epc)
    {
        case PressureCoupling::No: n = 0; break;
        case PressureCoupling::Berendsen:
        case PressureCoupling::CRescale:
        case PressureCoupling::Isotropic: n = c_minimumStepsPerTau; break;
        case PressureCoupling::ParrinelloRahman:
        case PressureCoupling::Mttk: n = c_minimumStepsPerPeriod; break;
        default: gmx_incons("Unknown epc value");
    }

    return n;
}

int ir_optimal_nstpcouple(const t_inputrec* ir)
{
    const int minIntegrationSteps = pcouple_min_integration_steps(ir->pressureCouplingOptions.epc);

    const int nwanted = c_defaultNstPCouple;

    // With multiple time-stepping we can only compute the pressure at slowest steps
    const int minNstPCouple = (ir->useMts ? ir->mtsLevels.back().stepFactor : 1);

    int n;
    if (minIntegrationSteps == 0
        || nwanted * ir->delta_t <= ir->pressureCouplingOptions.tau_p / minIntegrationSteps)
    {
        n = nwanted;
    }
    else
    {
        n = static_cast<int>(ir->pressureCouplingOptions.tau_p / (ir->delta_t * minIntegrationSteps) + 0.001);
        if (n < minNstPCouple)
        {
            n = minNstPCouple;
        }
        // Without MTS we try to make nstpcouple a "nice" number
        if (!ir->useMts)
        {
            while (nwanted % n != 0)
            {
                n--;
            }
        }
    }

    // With MTS, nstpcouple should be a multiple of the slowest MTS interval
    if (ir->useMts)
    {
        n = n - (n % minNstPCouple);
    }

    return n;
}

gmx_bool ir_coulomb_switched(const t_inputrec* ir)
{
    return (ir->coulombtype == CoulombInteractionType::Switch
            || ir->coulombtype == CoulombInteractionType::Shift
            || ir->coulombtype == CoulombInteractionType::PmeSwitch
            || ir->coulombtype == CoulombInteractionType::PmeUserSwitch
            || ir->coulomb_modifier == InteractionModifiers::PotSwitch
            || ir->coulomb_modifier == InteractionModifiers::ForceSwitch);
}

gmx_bool ir_coulomb_is_zero_at_cutoff(const t_inputrec* ir)
{
    return (ir->cutoff_scheme == CutoffScheme::Verlet || ir_coulomb_switched(ir)
            || ir->coulomb_modifier != InteractionModifiers::None
            || ir->coulombtype == CoulombInteractionType::RFZero);
}

gmx_bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec* ir)
{
    return (ir_coulomb_is_zero_at_cutoff(ir) || ir->coulombtype == CoulombInteractionType::User
            || ir->coulombtype == CoulombInteractionType::PmeUser);
}

gmx_bool ir_vdw_switched(const t_inputrec* ir)
{
    return (ir->vdwtype == VanDerWaalsType::Switch || ir->vdwtype == VanDerWaalsType::Shift
            || ir->vdw_modifier == InteractionModifiers::PotSwitch
            || ir->vdw_modifier == InteractionModifiers::ForceSwitch);
}

gmx_bool ir_vdw_is_zero_at_cutoff(const t_inputrec* ir)
{
    return (ir->cutoff_scheme == CutoffScheme::Verlet || ir_vdw_switched(ir)
            || ir->vdw_modifier != InteractionModifiers::None);
}

gmx_bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec* ir)
{
    return (ir_vdw_is_zero_at_cutoff(ir) || ir->vdwtype == VanDerWaalsType::User);
}

bool ir_haveBoxDeformation(const t_inputrec& ir)
{
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            if (ir.deform[i][j] != 0.0_real)
            {
                return true;
            }
        }
    }

    return false;
}

static void done_t_rot(t_rot* rot)
{
    if (rot == nullptr)
    {
        return;
    }
    for (auto& grp : rot->grp)
    {
        sfree(grp.ind);
    }
}

static void done_t_swapCoords(t_swapcoords* swapCoords)
{
    if (swapCoords == nullptr)
    {
        return;
    }
    for (int i = 0; i < swapCoords->ngrp; i++)
    {
        sfree(swapCoords->grp[i].ind);
        sfree(swapCoords->grp[i].molname);
    }
    sfree(swapCoords->grp);
    sfree(swapCoords);
}

void done_inputrec(t_inputrec* ir)
{
    sfree(ir->opts.nrdf);
    sfree(ir->opts.ref_t);
    for (int i = 0; i < ir->opts.ngtc; i++)
    {
        sfree(ir->opts.anneal_time[i]);
        sfree(ir->opts.anneal_temp[i]);
    }
    sfree(ir->opts.annealing);
    sfree(ir->opts.anneal_npoints);
    sfree(ir->opts.anneal_time);
    sfree(ir->opts.anneal_temp);
    sfree(ir->opts.tau_t);
    sfree(ir->opts.acceleration);
    sfree(ir->opts.nFreeze);
    sfree(ir->opts.egp_flags);

    done_t_swapCoords(ir->swap);
    done_t_rot(ir->rot.get());
    delete ir->params;
}

static void pr_grp_opts(FILE* out, int indent, const char* title, const t_grpopts* opts, gmx_bool bMDPformat)
{
    int i, m, j;

    if (!bMDPformat)
    {
        fprintf(out, "%s:\n", title);
    }

    pr_indent(out, indent);
    fprintf(out, "nrdf%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->nrdf[i]);
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "ref-t%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->ref_t[i]);
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "tau-t%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->tau_t[i]);
    }
    fprintf(out, "\n");

    /* Pretty-print the simulated annealing info */
    fprintf(out, "annealing%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10s", enumValueToString(opts->annealing[i]));
    }
    fprintf(out, "\n");

    fprintf(out, "annealing-npoints%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10d", opts->anneal_npoints[i]);
    }
    fprintf(out, "\n");

    for (i = 0; (i < opts->ngtc); i++)
    {
        if (opts->anneal_npoints[i] > 0)
        {
            fprintf(out, "annealing-time [%d]:\t", i);
            for (j = 0; (j < opts->anneal_npoints[i]); j++)
            {
                fprintf(out, "  %10.1f", opts->anneal_time[i][j]);
            }
            fprintf(out, "\n");
            fprintf(out, "annealing-temp [%d]:\t", i);
            for (j = 0; (j < opts->anneal_npoints[i]); j++)
            {
                fprintf(out, "  %10.1f", opts->anneal_temp[i][j]);
            }
            fprintf(out, "\n");
        }
    }

    pr_indent(out, indent);
    fprintf(out, "acc:\t");
    for (i = 0; (i < opts->ngacc); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            fprintf(out, "  %10g", opts->acceleration[i][m]);
        }
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "nfreeze:");
    for (i = 0; (i < opts->ngfrz); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            fprintf(out, "  %10s", opts->nFreeze[i][m] ? "Y" : "N");
        }
    }
    fprintf(out, "\n");


    for (i = 0; (i < opts->ngener); i++)
    {
        pr_indent(out, indent);
        fprintf(out, "energygrp-flags[%3d]:", i);
        for (m = 0; (m < opts->ngener); m++)
        {
            fprintf(out, " %d", opts->egp_flags[opts->ngener * i + m]);
        }
        fprintf(out, "\n");
    }

    fflush(out);
}

static void pr_matrix(FILE* fp, int indent, const char* title, const rvec* m, gmx_bool bMDPformat)
{
    if (bMDPformat)
    {
        fprintf(fp,
                "%-10s    = %g %g %g %g %g %g\n",
                title,
                m[XX][XX],
                m[YY][YY],
                m[ZZ][ZZ],
                m[XX][YY],
                m[XX][ZZ],
                m[YY][ZZ]);
    }
    else
    {
        pr_rvecs(fp, indent, title, m, DIM);
    }
}

#define PS(t, s) pr_str(fp, indent, t, s)
#define PI(t, s) pr_int(fp, indent, t, s)
#define PSTEP(t, s) pr_int64(fp, indent, t, s)
#define PR(t, s) pr_real(fp, indent, t, s)
#define PD(t, s) pr_double(fp, indent, t, s)

static void pr_pull_group(FILE* fp, int indent, int g, const t_pull_group* pgrp)
{
    pr_indent(fp, indent);
    fprintf(fp, "pull-group %d:\n", g);
    indent += 2;
    pr_ivec_block(fp, indent, "atom", pgrp->ind.data(), pgrp->ind.size(), TRUE);
    pr_rvec(fp, indent, "weight", pgrp->weight.data(), pgrp->weight.size(), TRUE);
    PI("pbcatom", pgrp->pbcatom);
}

static void pr_pull_coord(FILE* fp, int indent, int c, const t_pull_coord* pcrd)
{
    pr_indent(fp, indent);
    fprintf(fp, "pull-coord %d:\n", c);
    PS("type", enumValueToString(pcrd->eType));
    if (pcrd->eType == PullingAlgorithm::External)
    {
        PS("potential-provider", pcrd->externalPotentialProvider.c_str());
    }
    PS("geometry", enumValueToString(pcrd->eGeom));
    for (int g = 0; g < pcrd->ngroup; g++)
    {
        std::string buffer = gmx::formatString("group[%d]", g);
        PI(buffer.c_str(), pcrd->group[g]);
    }
    pr_ivec(fp, indent, "dim", pcrd->dim, DIM, TRUE);
    pr_rvec(fp, indent, "origin", pcrd->origin, DIM, TRUE);
    pr_rvec(fp, indent, "vec", pcrd->vec, DIM, TRUE);
    PS("start", EBOOL(pcrd->bStart));
    PR("init", pcrd->init);
    PR("rate", pcrd->rate);
    PR("k", pcrd->k);
    PR("kB", pcrd->kB);
}

static void pr_simtempvals(FILE* fp, int indent, const t_simtemp* simtemp, int n_lambda)
{
    PS("simulated-tempering-scaling", enumValueToString(simtemp->eSimTempScale));
    PR("sim-temp-low", simtemp->simtemp_low);
    PR("sim-temp-high", simtemp->simtemp_high);
    pr_rvec(fp, indent, "simulated tempering temperatures", simtemp->temperatures.data(), n_lambda, TRUE);
}

static void pr_expandedvals(FILE* fp, int indent, const t_expanded* expand, int n_lambda)
{

    PI("nstexpanded", expand->nstexpanded);
    PS("lmc-stats", enumValueToString(expand->elamstats));
    PS("lmc-move", enumValueToString(expand->elmcmove));
    PS("lmc-weights-equil", enumValueToString(expand->elmceq));
    if (expand->elmceq == LambdaWeightWillReachEquilibrium::NumAtLambda)
    {
        PI("weight-equil-number-all-lambda", expand->equil_n_at_lam);
    }
    if (expand->elmceq == LambdaWeightWillReachEquilibrium::Samples)
    {
        PI("weight-equil-number-samples", expand->equil_samples);
    }
    if (expand->elmceq == LambdaWeightWillReachEquilibrium::Steps)
    {
        PI("weight-equil-number-steps", expand->equil_steps);
    }
    if (expand->elmceq == LambdaWeightWillReachEquilibrium::WLDelta)
    {
        PR("weight-equil-wl-delta", expand->equil_wl_delta);
    }
    if (expand->elmceq == LambdaWeightWillReachEquilibrium::Ratio)
    {
        PR("weight-equil-count-ratio", expand->equil_ratio);
    }
    PI("lmc-seed", expand->lmc_seed);
    PR("mc-temperature", expand->mc_temp);
    PI("lmc-repeats", expand->lmc_repeats);
    PI("lmc-gibbsdelta", expand->gibbsdeltalam);
    PI("lmc-forced-nstart", expand->lmc_forced_nstart);
    PS("symmetrized-transition-matrix", EBOOL(expand->bSymmetrizedTMatrix));
    PI("nst-transition-matrix", expand->nstTij);
    PI("mininum-var-min", expand->minvarmin); /*default is reasonable */
    PI("weight-c-range", expand->c_range);    /* default is just C=0 */
    PR("wl-scale", expand->wl_scale);
    PR("wl-ratio", expand->wl_ratio);
    PR("init-wl-delta", expand->init_wl_delta);
    PS("wl-oneovert", EBOOL(expand->bWLoneovert));

    pr_indent(fp, indent);
    pr_rvec(fp, indent, "init-lambda-weights", expand->init_lambda_weights.data(), n_lambda, TRUE);
    PS("init-weights", EBOOL(expand->bInit_weights));
}

static void pr_fepvals(FILE* fp, int indent, const t_lambda* fep, gmx_bool bMDPformat)
{
    int j;

    PR("init-lambda", fep->init_lambda_without_states);
    PI("init-lambda-state", fep->init_fep_state);
    PR("delta-lambda", fep->delta_lambda);
    PI("nstdhdl", fep->nstdhdl);

    if (!bMDPformat)
    {
        PI("n-lambdas", fep->n_lambda);
    }
    if (fep->n_lambda > 0)
    {
        pr_indent(fp, indent);
        fprintf(fp, "separate-dvdl%s\n", bMDPformat ? " = " : ":");
        for (auto i : gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, bool>::keys())
        {
            fprintf(fp, "%18s = ", enumValueToString(i));
            if (fep->separate_dvdl[i])
            {
                fprintf(fp, "  TRUE");
            }
            else
            {
                fprintf(fp, "  FALSE");
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "all-lambdas%s\n", bMDPformat ? " = " : ":");
        for (auto key : gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, bool>::keys())
        {
            fprintf(fp, "%18s = ", enumValueToString(key));
            int i = static_cast<int>(key);
            for (j = 0; j < fep->n_lambda; j++)
            {
                fprintf(fp, "  %10g", fep->all_lambda[i][j]);
            }
            fprintf(fp, "\n");
        }
    }
    PI("calc-lambda-neighbors", fep->lambda_neighbors);
    PS("dhdl-print-energy", enumValueToString(fep->edHdLPrintEnergy));
    PR("sc-alpha", fep->sc_alpha);
    PI("sc-power", fep->sc_power);
    PR("sc-r-power", fep->sc_r_power);
    PR("sc-sigma", fep->sc_sigma);
    PR("sc-sigma-min", fep->sc_sigma_min);
    PS("sc-coul", EBOOL(fep->bScCoul));
    PI("dh-hist-size", fep->dh_hist_size);
    PD("dh-hist-spacing", fep->dh_hist_spacing);
    PS("separate-dhdl-file", enumValueToString(fep->separate_dhdl_file));
    PS("dhdl-derivatives", enumValueToString(fep->dhdl_derivatives));
    PS("sc-function", enumValueToString(fep->softcoreFunction));
    PR("sc-gapsys-scale-linpoint-lj", fep->scGapsysScaleLinpointLJ);
    PR("sc-gapsys-scale-linpoint-q", fep->scGapsysScaleLinpointQ);
    PR("sc-gapsys-sigma-lj", fep->scGapsysSigmaLJ);
};

static void pr_pull(FILE* fp, int indent, const pull_params_t& pull)
{
    int g;

    PR("pull-cylinder-r", pull.cylinder_r);
    PR("pull-constr-tol", pull.constr_tol);
    PS("pull-print-COM", EBOOL(pull.bPrintCOM));
    PS("pull-print-ref-value", EBOOL(pull.bPrintRefValue));
    PS("pull-print-components", EBOOL(pull.bPrintComp));
    PI("pull-nstxout", pull.nstxout);
    PI("pull-nstfout", pull.nstfout);
    PS("pull-pbc-ref-prev-step-com", EBOOL(pull.bSetPbcRefToPrevStepCOM));
    PS("pull-xout-average", EBOOL(pull.bXOutAverage));
    PS("pull-fout-average", EBOOL(pull.bFOutAverage));
    PI("pull-ngroups", pull.ngroup);
    for (g = 0; g < pull.ngroup; g++)
    {
        pr_pull_group(fp, indent, g, &pull.group[g]);
    }
    PI("pull-ncoords", pull.ncoord);
    for (g = 0; g < pull.ncoord; g++)
    {
        pr_pull_coord(fp, indent, g, &pull.coord[g]);
    }
}

static void pr_awh_bias_dim(FILE* fp, int indent, const gmx::AwhDimParams& awhDimParams, const char* prefix)
{
    pr_indent(fp, indent);
    indent++;
    fprintf(fp, "%s:\n", prefix);
    PS("coord-provider", enumValueToString(awhDimParams.coordinateProvider()));
    PI("coord-index", awhDimParams.coordinateIndex() + 1);
    PR("start", awhDimParams.origin());
    PR("end", awhDimParams.end());
    PR("period", awhDimParams.period());
    PR("force-constant", awhDimParams.forceConstant());
    PR("diffusion", awhDimParams.diffusion());
    PR("cover-diameter", awhDimParams.coverDiameter());
}

static void pr_awh_bias(FILE* fp, int indent, const gmx::AwhBiasParams& awhBiasParams, const char* prefix)
{
    char opt[STRLEN];

    sprintf(opt, "%s-error-init", prefix);
    PR(opt, awhBiasParams.initialErrorEstimate());
    sprintf(opt, "%s-growth", prefix);
    PS(opt, enumValueToString(awhBiasParams.growthType()));
    sprintf(opt, "%s-growth-factor", prefix);
    PR(opt, awhBiasParams.growthFactor());
    sprintf(opt, "%s-target", prefix);
    PS(opt, enumValueToString(awhBiasParams.targetDistribution()));
    sprintf(opt, "%s-target-beta-scaling", prefix);
    PR(opt, awhBiasParams.targetBetaScaling());
    sprintf(opt, "%s-target-cutoff", prefix);
    PR(opt, awhBiasParams.targetCutoff());
    sprintf(opt, "%s-target-metric-scaling", prefix);
    PS(opt, EBOOL(awhBiasParams.scaleTargetByMetric()));
    sprintf(opt, "%s-target-metric-scaling-limit", prefix);
    PR(opt, awhBiasParams.targetMetricScalingLimit());
    sprintf(opt, "%s-user-data", prefix);
    PS(opt, EBOOL(awhBiasParams.userPMFEstimate()));
    sprintf(opt, "%s-share-group", prefix);
    PI(opt, awhBiasParams.shareGroup());
    sprintf(opt, "%s-equilibrate-histogram", prefix);
    PS(opt, EBOOL(awhBiasParams.equilibrateHistogram()));
    sprintf(opt, "%s-ndim", prefix);
    PI(opt, awhBiasParams.ndim());

    int d = 0;
    for (const auto& dimParam : awhBiasParams.dimParams())
    {
        char prefixdim[STRLEN];
        sprintf(prefixdim, "%s-dim%d", prefix, d + 1);
        pr_awh_bias_dim(fp, indent, dimParam, prefixdim);
        d++;
    }
}

static void pr_awh(FILE* fp, int indent, gmx::AwhParams* awhParams)
{
    PS("awh-potential", enumValueToString(awhParams->potential()));
    PI("awh-seed", awhParams->seed());
    PI("awh-nstout", awhParams->nstout());
    PI("awh-nstsample", awhParams->nstSampleCoord());
    PI("awh-nsamples-update", awhParams->numSamplesUpdateFreeEnergy());
    PS("awh-share-bias-multisim", EBOOL(awhParams->shareBiasMultisim()));
    PI("awh-nbias", awhParams->numBias());

    int k = 0;
    for (const auto& awhBiasParam : awhParams->awhBiasParams())
    {
        auto prefix = gmx::formatString("awh%d", k + 1);
        pr_awh_bias(fp, indent, awhBiasParam, prefix.c_str());
        k++;
    }
}

static void pr_rotgrp(FILE* fp, int indent, int g, const t_rotgrp* rotg)
{
    pr_indent(fp, indent);
    fprintf(fp, "rot-group %d:\n", g);
    indent += 2;
    PS("rot-type", enumValueToString(rotg->eType));
    PS("rot-massw", EBOOL(rotg->bMassW));
    pr_ivec_block(fp, indent, "atom", rotg->ind, rotg->nat, TRUE);
    pr_rvecs(fp, indent, "x-ref", as_vec_array(rotg->x_ref_original.data()), rotg->x_ref_original.size());
    pr_rvec(fp, indent, "rot-vec", rotg->inputVec, DIM, TRUE);
    pr_rvec(fp, indent, "rot-pivot", rotg->pivot, DIM, TRUE);
    PR("rot-rate", rotg->rate);
    PR("rot-k", rotg->k);
    PR("rot-slab-dist", rotg->slab_dist);
    PR("rot-min-gauss", rotg->min_gaussian);
    PR("rot-eps", rotg->eps);
    PS("rot-fit-method", enumValueToString(rotg->eFittype));
    PI("rot-potfit-nstep", rotg->PotAngle_nstep);
    PR("rot-potfit-step", rotg->PotAngle_step);
}

static void pr_rot(FILE* fp, int indent, const t_rot* rot)
{
    int g;

    PI("rot-nstrout", rot->nstrout);
    PI("rot-nstsout", rot->nstsout);
    PI("rot-ngroups", rot->grp.size());
    for (g = 0; g < gmx::ssize(rot->grp); g++)
    {
        pr_rotgrp(fp, indent, g, &rot->grp[g]);
    }
}


static void pr_swap(FILE* fp, int indent, const t_swapcoords* swap)
{
    char str[STRLEN];

    /* Enums for better readability of the code */
    enum
    {
        eCompA = 0,
        eCompB
    };


    PI("swap-frequency", swap->nstswap);

    /* The split groups that define the compartments */
    for (int j = 0; j < 2; j++)
    {
        snprintf(str, STRLEN, "massw_split%d", j);
        PS(str, EBOOL(swap->massw_split[j]));
        snprintf(str, STRLEN, "split atoms group %d", j);
        pr_ivec_block(fp, indent, str, swap->grp[j].ind, swap->grp[j].nat, TRUE);
    }

    /* The solvent group */
    snprintf(str,
             STRLEN,
             "solvent group %s",
             swap->grp[static_cast<int>(SwapGroupSplittingType::Solvent)].molname);
    pr_ivec_block(fp,
                  indent,
                  str,
                  swap->grp[static_cast<int>(SwapGroupSplittingType::Solvent)].ind,
                  swap->grp[static_cast<int>(SwapGroupSplittingType::Solvent)].nat,
                  TRUE);

    /* Now print the indices for all the ion groups: */
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < swap->ngrp; ig++)
    {
        snprintf(str, STRLEN, "ion group %s", swap->grp[ig].molname);
        pr_ivec_block(fp, indent, str, swap->grp[ig].ind, swap->grp[ig].nat, TRUE);
    }

    PR("cyl0-r", swap->cyl0r);
    PR("cyl0-up", swap->cyl0u);
    PR("cyl0-down", swap->cyl0l);
    PR("cyl1-r", swap->cyl1r);
    PR("cyl1-up", swap->cyl1u);
    PR("cyl1-down", swap->cyl1l);
    PI("coupl-steps", swap->nAverage);

    /* Print the requested ion counts for both compartments */
    for (int ic = eCompA; ic <= eCompB; ic++)
    {
        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < swap->ngrp; ig++)
        {
            snprintf(str, STRLEN, "%s-in-%c", swap->grp[ig].molname, 'A' + ic);
            PI(str, swap->grp[ig].nmolReq[ic]);
        }
    }

    PR("threshold", swap->threshold);
    PR("bulk-offsetA", swap->bulkOffset[eCompA]);
    PR("bulk-offsetB", swap->bulkOffset[eCompB]);
}


static void pr_imd(FILE* fp, int indent, const t_IMD* imd)
{
    PI("IMD-atoms", imd->nat);
    pr_ivec_block(fp, indent, "atom", imd->ind, imd->nat, TRUE);
}


void pr_inputrec(FILE* fp, int indent, const char* title, const t_inputrec* ir, gmx_bool bMDPformat)
{
    const char* infbuf = "inf";

    if (available(fp, ir, indent, title))
    {
        if (!bMDPformat)
        {
            indent = pr_title(fp, indent, title);
        }
        /* Try to make this list appear in the same order as the
         * options are written in the default mdout.mdp, and with
         * the same user-exposed names to facilitate debugging.
         */
        PS("integrator", enumValueToString(ir->eI));
        PR("tinit", ir->init_t);
        PR("dt", ir->delta_t);
        PSTEP("nsteps", ir->nsteps);
        PSTEP("init-step", ir->init_step);
        PI("simulation-part", ir->simulation_part);
        PS("mts", EBOOL(ir->useMts));
        if (ir->useMts)
        {
            for (int mtsIndex = 1; mtsIndex < static_cast<int>(ir->mtsLevels.size()); mtsIndex++)
            {
                const auto&       mtsLevel = ir->mtsLevels[mtsIndex];
                const std::string forceKey = gmx::formatString("mts-level%d-forces", mtsIndex + 1);
                std::string       forceGroups;
                for (int i = 0; i < static_cast<int>(gmx::MtsForceGroups::Count); i++)
                {
                    if (mtsLevel.forceGroups[i])
                    {
                        if (!forceGroups.empty())
                        {
                            forceGroups += " ";
                        }
                        forceGroups += gmx::mtsForceGroupNames[i];
                    }
                }
                PS(forceKey.c_str(), forceGroups.c_str());
                const std::string factorKey = gmx::formatString("mts-level%d-factor", mtsIndex + 1);
                PI(factorKey.c_str(), mtsLevel.stepFactor);
            }
        }
        PR("mass-repartition-factor", ir->massRepartitionFactor);
        PS("comm-mode", enumValueToString(ir->comm_mode));
        PI("nstcomm", ir->nstcomm);

        /* Langevin dynamics */
        PR("bd-fric", ir->bd_fric);
        PSTEP("ld-seed", ir->ld_seed);

        /* Energy minimization */
        PR("emtol", ir->em_tol);
        PR("emstep", ir->em_stepsize);
        PI("niter", ir->niter);
        PR("fcstep", ir->fc_stepsize);
        PI("nstcgsteep", ir->nstcgsteep);
        PI("nbfgscorr", ir->nbfgscorr);

        /* Test particle insertion */
        PR("rtpi", ir->rtpi);

        /* Output control */
        PI("nstxout", ir->nstxout);
        PI("nstvout", ir->nstvout);
        PI("nstfout", ir->nstfout);
        PI("nstlog", ir->nstlog);
        PI("nstcalcenergy", ir->nstcalcenergy);
        PI("nstenergy", ir->nstenergy);
        PI("nstxout-compressed", ir->nstxout_compressed);
        PR("compressed-x-precision", ir->x_compression_precision);

        /* Neighborsearching parameters */
        PS("cutoff-scheme", enumValueToString(ir->cutoff_scheme));
        PI("nstlist", ir->nstlist);
        PS("pbc", c_pbcTypeNames[ir->pbcType].c_str());
        PS("periodic-molecules", EBOOL(ir->bPeriodicMols));
        PR("verlet-buffer-tolerance", ir->verletbuf_tol);
        PR("verlet-buffer-pressure-tolerance", ir->verletBufferPressureTolerance);
        PR("rlist", ir->rlist);

        /* Options for electrostatics and VdW */
        PS("coulombtype", enumValueToString(ir->coulombtype));
        PS("coulomb-modifier", enumValueToString(ir->coulomb_modifier));
        PR("rcoulomb-switch", ir->rcoulomb_switch);
        PR("rcoulomb", ir->rcoulomb);
        if (ir->epsilon_r != 0)
        {
            PR("epsilon-r", ir->epsilon_r);
        }
        else
        {
            PS("epsilon-r", infbuf);
        }
        if (ir->epsilon_rf != 0)
        {
            PR("epsilon-rf", ir->epsilon_rf);
        }
        else
        {
            PS("epsilon-rf", infbuf);
        }
        PS("vdw-type", enumValueToString(ir->vdwtype));
        PS("vdw-modifier", enumValueToString(ir->vdw_modifier));
        PR("rvdw-switch", ir->rvdw_switch);
        PR("rvdw", ir->rvdw);
        PS("DispCorr", enumValueToString(ir->eDispCorr));
        PR("table-extension", ir->tabext);

        PR("fourierspacing", ir->fourier_spacing);
        PI("fourier-nx", ir->nkx);
        PI("fourier-ny", ir->nky);
        PI("fourier-nz", ir->nkz);
        PI("pme-order", ir->pme_order);
        PR("ewald-rtol", ir->ewald_rtol);
        PR("ewald-rtol-lj", ir->ewald_rtol_lj);
        PS("lj-pme-comb-rule", enumValueToString(ir->ljpme_combination_rule));
        PS("ewald-geometry", enumValueToString(ir->ewald_geometry));
        PR("epsilon-surface", ir->epsilon_surface);

        /* Options for weak coupling algorithms */
        PS("ensemble-temperature-setting", enumValueToString(ir->ensembleTemperatureSetting));
        if (ir->ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant)
        {
            PR("ensemble-temperature", ir->ensembleTemperature);
        }
        PS("tcoupl", enumValueToString(ir->etc));
        PI("nsttcouple", ir->nsttcouple);
        PI("nh-chain-length", ir->opts.nhchainlength);
        PS("print-nose-hoover-chain-variables", EBOOL(ir->bPrintNHChains));

        PS("pcoupl", enumValueToString(ir->pressureCouplingOptions.epc));
        if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
        {
            PS("pcoupltype", enumValueToString(ir->pressureCouplingOptions.epct));
            PI("nstpcouple", ir->pressureCouplingOptions.nstpcouple);
            PR("tau-p", ir->pressureCouplingOptions.tau_p);
            pr_matrix(fp, indent, "compressibility", ir->pressureCouplingOptions.compress, bMDPformat);
            pr_matrix(fp, indent, "ref-p", ir->pressureCouplingOptions.ref_p, bMDPformat);
        }
        // Refcoord-scaling is also needed for other algorithms that affect the box
        PS("refcoord-scaling", enumValueToString(ir->pressureCouplingOptions.refcoord_scaling));

        if (bMDPformat)
        {
            fprintf(fp, "posres-com  = %g %g %g\n", ir->posres_com[XX], ir->posres_com[YY], ir->posres_com[ZZ]);
            fprintf(fp, "posres-comB = %g %g %g\n", ir->posres_comB[XX], ir->posres_comB[YY], ir->posres_comB[ZZ]);
        }
        else
        {
            pr_rvec(fp, indent, "posres-com", ir->posres_com, DIM, TRUE);
            pr_rvec(fp, indent, "posres-comB", ir->posres_comB, DIM, TRUE);
        }

        /* QMMM */
        PS("QMMM", EBOOL(ir->bQMMM));
        fprintf(fp, "%s:\n", "qm-opts");
        pr_int(fp, indent, "ngQM", ir->opts.ngQM);

        /* CONSTRAINT OPTIONS */
        PS("constraint-algorithm", enumValueToString(ir->eConstrAlg));
        PS("continuation", EBOOL(ir->bContinuation));

        PS("Shake-SOR", EBOOL(ir->bShakeSOR));
        PR("shake-tol", ir->shake_tol);
        PI("lincs-order", ir->nProjOrder);
        PI("lincs-iter", ir->nLincsIter);
        PR("lincs-warnangle", ir->LincsWarnAngle);

        /* Walls */
        PI("nwall", ir->nwall);
        PS("wall-type", enumValueToString(ir->wall_type));
        PR("wall-r-linpot", ir->wall_r_linpot);
        /* wall-atomtype */
        PI("wall-atomtype[0]", ir->wall_atomtype[0]);
        PI("wall-atomtype[1]", ir->wall_atomtype[1]);
        /* wall-density */
        PR("wall-density[0]", ir->wall_density[0]);
        PR("wall-density[1]", ir->wall_density[1]);
        PR("wall-ewald-zfac", ir->wall_ewald_zfac);

        /* COM PULLING */
        PS("pull", EBOOL(ir->bPull));
        if (ir->bPull)
        {
            pr_pull(fp, indent, *ir->pull);
        }

        /* AWH BIASING */
        PS("awh", EBOOL(ir->bDoAwh));
        if (ir->bDoAwh)
        {
            pr_awh(fp, indent, ir->awhParams.get());
        }

        /* ENFORCED ROTATION */
        PS("rotation", EBOOL(ir->bRot));
        if (ir->bRot)
        {
            pr_rot(fp, indent, ir->rot.get());
        }

        /* INTERACTIVE MD */
        PS("interactiveMD", EBOOL(ir->bIMD));
        if (ir->bIMD)
        {
            pr_imd(fp, indent, ir->imd);
        }

        /* NMR refinement stuff */
        PS("disre", enumValueToString(ir->eDisre));
        PS("disre-weighting", enumValueToString(ir->eDisreWeighting));
        PS("disre-mixed", EBOOL(ir->bDisreMixed));
        PR("dr-fc", ir->dr_fc);
        PR("dr-tau", ir->dr_tau);
        PR("nstdisreout", ir->nstdisreout);

        PR("orire-fc", ir->orires_fc);
        PR("orire-tau", ir->orires_tau);
        PR("nstorireout", ir->nstorireout);

        /* FREE ENERGY VARIABLES */
        PS("free-energy", enumValueToString(ir->efep));
        if (ir->efep != FreeEnergyPerturbationType::No || ir->bSimTemp)
        {
            pr_fepvals(fp, indent, ir->fepvals.get(), bMDPformat);
        }
        if (ir->bExpanded)
        {
            pr_expandedvals(fp, indent, ir->expandedvals.get(), ir->fepvals->n_lambda);
        }

        /* NON-equilibrium MD stuff */
        PR("cos-acceleration", ir->cos_accel);
        pr_matrix(fp, indent, "deform", ir->deform, bMDPformat);

        /* SIMULATED TEMPERING */
        PS("simulated-tempering", EBOOL(ir->bSimTemp));
        if (ir->bSimTemp)
        {
            pr_simtempvals(fp, indent, ir->simtempvals.get(), ir->fepvals->n_lambda);
        }

        /* ION/WATER SWAPPING FOR COMPUTATIONAL ELECTROPHYSIOLOGY */
        PS("swapcoords", enumValueToString(ir->eSwapCoords));
        if (ir->eSwapCoords != SwapType::No)
        {
            pr_swap(fp, indent, ir->swap);
        }

        /* USER-DEFINED THINGIES */
        PI("userint1", ir->userint1);
        PI("userint2", ir->userint2);
        PI("userint3", ir->userint3);
        PI("userint4", ir->userint4);
        PR("userreal1", ir->userreal1);
        PR("userreal2", ir->userreal2);
        PR("userreal3", ir->userreal3);
        PR("userreal4", ir->userreal4);

        if (!bMDPformat)
        {
            gmx::TextWriter writer(fp);
            writer.wrapperSettings().setIndent(indent);
            gmx::dumpKeyValueTree(&writer, *ir->params);
        }

        pr_grp_opts(fp, indent, "grpopts", &(ir->opts), bMDPformat);
    }
}
#undef PS
#undef PR
#undef PI

static void cmpPressureCouplingOptions(FILE*                          fp,
                                       const PressureCouplingOptions& pco1,
                                       const PressureCouplingOptions& pco2,
                                       real                           ftol,
                                       real                           abstol)
{
    GMX_RELEASE_ASSERT(pco1.epc != PressureCoupling::No && pco2.epc != PressureCoupling::No,
                       "Need p-coupling");

    cmpEnum(fp, "inputrec->pressureCouplingOptions.epct", pco1.epct, pco2.epct);
    cmp_real(fp, "inputrec->pressureCouplingOptions.tau_p", -1, pco1.tau_p, pco2.tau_p, ftol, abstol);
    cmp_rvec(fp, "inputrec->pressureCouplingOptions.ref_p(x)", -1, pco1.ref_p[XX], pco2.ref_p[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->pressureCouplingOptions.ref_p(y)", -1, pco1.ref_p[YY], pco2.ref_p[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->pressureCouplingOptions.ref_p(z)", -1, pco1.ref_p[ZZ], pco2.ref_p[ZZ], ftol, abstol);
    cmp_rvec(fp, "inputrec->pressureCouplingOptions.compress(x)", -1, pco1.compress[XX], pco2.compress[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->pressureCouplingOptions.compress(y)", -1, pco1.compress[YY], pco2.compress[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->pressureCouplingOptions.compress(z)", -1, pco1.compress[ZZ], pco2.compress[ZZ], ftol, abstol);
}

static void cmp_grpopts(FILE* fp, const t_grpopts* opt1, const t_grpopts* opt2, real ftol, real abstol)
{
    int  i, j;
    char buf1[256], buf2[256];

    cmp_int(fp, "inputrec->grpopts.ngtc", -1, opt1->ngtc, opt2->ngtc);
    cmp_int(fp, "inputrec->grpopts.ngacc", -1, opt1->ngacc, opt2->ngacc);
    cmp_int(fp, "inputrec->grpopts.ngfrz", -1, opt1->ngfrz, opt2->ngfrz);
    cmp_int(fp, "inputrec->grpopts.ngener", -1, opt1->ngener, opt2->ngener);
    for (i = 0; (i < std::min(opt1->ngtc, opt2->ngtc)); i++)
    {
        cmp_real(fp, "inputrec->grpopts.nrdf", i, opt1->nrdf[i], opt2->nrdf[i], ftol, abstol);
        cmp_real(fp, "inputrec->grpopts.ref_t", i, opt1->ref_t[i], opt2->ref_t[i], ftol, abstol);
        cmp_real(fp, "inputrec->grpopts.tau_t", i, opt1->tau_t[i], opt2->tau_t[i], ftol, abstol);
        cmpEnum(fp, "inputrec->grpopts.annealing", opt1->annealing[i], opt2->annealing[i]);
        cmp_int(fp, "inputrec->grpopts.anneal_npoints", i, opt1->anneal_npoints[i], opt2->anneal_npoints[i]);
        if (opt1->anneal_npoints[i] == opt2->anneal_npoints[i])
        {
            sprintf(buf1, "inputrec->grpopts.anneal_time[%d]", i);
            sprintf(buf2, "inputrec->grpopts.anneal_temp[%d]", i);
            for (j = 0; j < opt1->anneal_npoints[i]; j++)
            {
                cmp_real(fp, buf1, j, opt1->anneal_time[i][j], opt2->anneal_time[i][j], ftol, abstol);
                cmp_real(fp, buf2, j, opt1->anneal_temp[i][j], opt2->anneal_temp[i][j], ftol, abstol);
            }
        }
    }
    if (opt1->ngener == opt2->ngener)
    {
        for (i = 0; i < opt1->ngener; i++)
        {
            for (j = i; j < opt1->ngener; j++)
            {
                sprintf(buf1, "inputrec->grpopts.egp_flags[%d]", i);
                cmp_int(fp, buf1, j, opt1->egp_flags[opt1->ngener * i + j], opt2->egp_flags[opt1->ngener * i + j]);
            }
        }
    }
    for (i = 0; (i < std::min(opt1->ngacc, opt2->ngacc)); i++)
    {
        cmp_rvec(fp, "inputrec->grpopts.acceleration", i, opt1->acceleration[i], opt2->acceleration[i], ftol, abstol);
    }
    for (i = 0; (i < std::min(opt1->ngfrz, opt2->ngfrz)); i++)
    {
        cmp_ivec(fp, "inputrec->grpopts.nFreeze", i, opt1->nFreeze[i], opt2->nFreeze[i]);
    }
}

static void cmp_pull(FILE* fp)
{
    fprintf(fp,
            "WARNING: Both files use COM pulling, but comparing of the pull struct is not "
            "implemented (yet). The pull parameters could be the same or different.\n");
}

static void cmp_awhDimParams(FILE*                    fp,
                             const gmx::AwhDimParams& dimp1,
                             const gmx::AwhDimParams& dimp2,
                             int                      dimIndex,
                             real                     ftol,
                             real                     abstol)
{
    /* Note that we have double index here, but the compare functions only
     * support one index, so here we only print the dim index and not the bias.
     */
    cmp_int(fp,
            "inputrec.awhParams->bias?->dim->coord_index",
            dimIndex,
            dimp1.coordinateIndex(),
            dimp2.coordinateIndex());
    cmp_double(fp, "inputrec->awhParams->bias?->dim->period", dimIndex, dimp1.period(), dimp2.period(), ftol, abstol);
    cmp_double(fp,
               "inputrec->awhParams->bias?->dim->diffusion",
               dimIndex,
               dimp1.diffusion(),
               dimp2.diffusion(),
               ftol,
               abstol);
    cmp_double(fp, "inputrec->awhParams->bias?->dim->origin", dimIndex, dimp1.origin(), dimp2.origin(), ftol, abstol);
    cmp_double(fp, "inputrec->awhParams->bias?->dim->end", dimIndex, dimp1.end(), dimp2.end(), ftol, abstol);
    cmp_double(fp,
               "inputrec->awhParams->bias?->dim->coord_value_init",
               dimIndex,
               dimp1.initialCoordinate(),
               dimp2.initialCoordinate(),
               ftol,
               abstol);
    cmp_double(fp,
               "inputrec->awhParams->bias?->dim->coverDiameter",
               dimIndex,
               dimp1.coverDiameter(),
               dimp2.coverDiameter(),
               ftol,
               abstol);
}

static void cmp_awhBiasParams(FILE*                     fp,
                              const gmx::AwhBiasParams& bias1,
                              const gmx::AwhBiasParams& bias2,
                              int                       biasIndex,
                              real                      ftol,
                              real                      abstol)
{
    cmp_int(fp, "inputrec->awhParams->ndim", biasIndex, bias1.ndim(), bias2.ndim());
    cmpEnum<gmx::AwhTargetType>(
            fp, "inputrec->awhParams->biaseTarget", bias1.targetDistribution(), bias2.targetDistribution());
    cmp_double(fp,
               "inputrec->awhParams->biastargetBetaScaling",
               biasIndex,
               bias1.targetBetaScaling(),
               bias2.targetBetaScaling(),
               ftol,
               abstol);
    cmp_double(fp,
               "inputrec->awhParams->biastargetCutoff",
               biasIndex,
               bias1.targetCutoff(),
               bias2.targetCutoff(),
               ftol,
               abstol);
    cmpEnum<gmx::AwhHistogramGrowthType>(
            fp, "inputrec->awhParams->biaseGrowth", bias1.growthType(), bias2.growthType());
    cmp_double(fp,
               "inputrec->awhParams->biasGrowthFactor",
               biasIndex,
               bias1.growthFactor(),
               bias2.growthFactor(),
               ftol,
               abstol);
    cmp_bool(fp, "inputrec->awhParams->biasbUserData", biasIndex, bias1.userPMFEstimate(), bias2.userPMFEstimate());
    cmp_double(fp,
               "inputrec->awhParams->biaserror_initial",
               biasIndex,
               bias1.initialErrorEstimate(),
               bias2.initialErrorEstimate(),
               ftol,
               abstol);
    cmp_int(fp, "inputrec->awhParams->biasShareGroup", biasIndex, bias1.shareGroup(), bias2.shareGroup());

    const auto dimParams1 = bias1.dimParams();
    const auto dimParams2 = bias2.dimParams();
    for (int dim = 0; dim < std::min(bias1.ndim(), bias2.ndim()); dim++)
    {
        cmp_awhDimParams(fp, dimParams1[dim], dimParams2[dim], dim, ftol, abstol);
    }
}

static void cmp_awhParams(FILE* fp, const gmx::AwhParams& awh1, const gmx::AwhParams& awh2, real ftol, real abstol)
{
    cmp_int(fp, "inputrec->awhParams->nbias", -1, awh1.numBias(), awh2.numBias());
    cmp_int64(fp, "inputrec->awhParams->seed", awh1.seed(), awh2.seed());
    cmp_int(fp, "inputrec->awhParams->nstout", -1, awh1.nstout(), awh2.nstout());
    cmp_int(fp, "inputrec->awhParams->nstsample_coord", -1, awh1.nstSampleCoord(), awh2.nstSampleCoord());
    cmp_int(fp,
            "inputrec->awhParams->nsamples_update_free_energy",
            -1,
            awh1.numSamplesUpdateFreeEnergy(),
            awh2.numSamplesUpdateFreeEnergy());
    cmpEnum<gmx::AwhPotentialType>(
            fp, "inputrec->awhParams->ePotential", awh1.potential(), awh2.potential());
    cmp_bool(fp, "inputrec->awhParams->shareBiasMultisim", -1, awh1.shareBiasMultisim(), awh2.shareBiasMultisim());

    if (awh1.numBias() == awh2.numBias())
    {
        const auto awhBiasParams1 = awh1.awhBiasParams();
        const auto awhBiasParams2 = awh2.awhBiasParams();
        for (int bias = 0; bias < awh1.numBias(); bias++)
        {
            cmp_awhBiasParams(fp, awhBiasParams1[bias], awhBiasParams2[bias], bias, ftol, abstol);
        }
    }
}

static void cmp_simtempvals(FILE*            fp,
                            const t_simtemp* simtemp1,
                            const t_simtemp* simtemp2,
                            int              n_lambda,
                            real             ftol,
                            real             abstol)
{
    int i;
    cmpEnum(fp, "inputrec->simtempvals->eSimTempScale", simtemp1->eSimTempScale, simtemp2->eSimTempScale);
    cmp_real(fp, "inputrec->simtempvals->simtemp_high", -1, simtemp1->simtemp_high, simtemp2->simtemp_high, ftol, abstol);
    cmp_real(fp, "inputrec->simtempvals->simtemp_low", -1, simtemp1->simtemp_low, simtemp2->simtemp_low, ftol, abstol);
    for (i = 0; i < n_lambda; i++)
    {
        cmp_real(fp,
                 "inputrec->simtempvals->temperatures",
                 -1,
                 simtemp1->temperatures[i],
                 simtemp2->temperatures[i],
                 ftol,
                 abstol);
    }
}

static void cmp_expandedvals(FILE*             fp,
                             const t_expanded* expand1,
                             const t_expanded* expand2,
                             int               n_lambda,
                             real              ftol,
                             real              abstol)
{
    int i;

    cmp_bool(fp, "inputrec->fepvals->bInit_weights", -1, expand1->bInit_weights, expand2->bInit_weights);
    cmp_bool(fp, "inputrec->fepvals->bWLoneovert", -1, expand1->bWLoneovert, expand2->bWLoneovert);

    for (i = 0; i < n_lambda; i++)
    {
        cmp_real(fp,
                 "inputrec->expandedvals->init_lambda_weights",
                 -1,
                 expand1->init_lambda_weights[i],
                 expand2->init_lambda_weights[i],
                 ftol,
                 abstol);
    }

    cmpEnum(fp, "inputrec->expandedvals->lambda-stats", expand1->elamstats, expand2->elamstats);
    cmpEnum(fp, "inputrec->expandedvals->lambda-mc-move", expand1->elmcmove, expand2->elmcmove);
    cmp_int(fp, "inputrec->expandedvals->lmc-repeats", -1, expand1->lmc_repeats, expand2->lmc_repeats);
    cmp_int(fp, "inputrec->expandedvals->lmc-gibbsdelta", -1, expand1->gibbsdeltalam, expand2->gibbsdeltalam);
    cmp_int(fp, "inputrec->expandedvals->lmc-forced-nstart", -1, expand1->lmc_forced_nstart, expand2->lmc_forced_nstart);
    cmpEnum(fp, "inputrec->expandedvals->lambda-weights-equil", expand1->elmceq, expand2->elmceq);
    cmp_int(fp,
            "inputrec->expandedvals->,weight-equil-number-all-lambda",
            -1,
            expand1->equil_n_at_lam,
            expand2->equil_n_at_lam);
    cmp_int(fp, "inputrec->expandedvals->weight-equil-number-samples", -1, expand1->equil_samples, expand2->equil_samples);
    cmp_int(fp, "inputrec->expandedvals->weight-equil-number-steps", -1, expand1->equil_steps, expand2->equil_steps);
    cmp_real(fp,
             "inputrec->expandedvals->weight-equil-wl-delta",
             -1,
             expand1->equil_wl_delta,
             expand2->equil_wl_delta,
             ftol,
             abstol);
    cmp_real(fp,
             "inputrec->expandedvals->weight-equil-count-ratio",
             -1,
             expand1->equil_ratio,
             expand2->equil_ratio,
             ftol,
             abstol);
    cmp_bool(fp,
             "inputrec->expandedvals->symmetrized-transition-matrix",
             -1,
             expand1->bSymmetrizedTMatrix,
             expand2->bSymmetrizedTMatrix);
    cmp_int(fp, "inputrec->expandedvals->nstTij", -1, expand1->nstTij, expand2->nstTij);
    cmp_int(fp, "inputrec->expandedvals->mininum-var-min", -1, expand1->minvarmin, expand2->minvarmin); /*default is reasonable */
    cmp_int(fp, "inputrec->expandedvals->weight-c-range", -1, expand1->c_range, expand2->c_range); /* default is just C=0 */
    cmp_real(fp, "inputrec->expandedvals->wl-scale", -1, expand1->wl_scale, expand2->wl_scale, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->init-wl-delta", -1, expand1->init_wl_delta, expand2->init_wl_delta, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->wl-ratio", -1, expand1->wl_ratio, expand2->wl_ratio, ftol, abstol);
    cmp_int(fp, "inputrec->expandedvals->nstexpanded", -1, expand1->nstexpanded, expand2->nstexpanded);
    cmp_int(fp, "inputrec->expandedvals->lmc-seed", -1, expand1->lmc_seed, expand2->lmc_seed);
    cmp_real(fp, "inputrec->expandedvals->mc-temperature", -1, expand1->mc_temp, expand2->mc_temp, ftol, abstol);
}

static void cmp_fepvals(FILE* fp, const t_lambda* fep1, const t_lambda* fep2, real ftol, real abstol)
{
    int i, j;
    cmp_int(fp, "inputrec->nstdhdl", -1, fep1->nstdhdl, fep2->nstdhdl);
    cmp_double(fp, "inputrec->fepvals->init_fep_state", -1, fep1->init_fep_state, fep2->init_fep_state, ftol, abstol);
    cmp_double(fp, "inputrec->fepvals->delta_lambda", -1, fep1->delta_lambda, fep2->delta_lambda, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->n_lambda", -1, fep1->n_lambda, fep2->n_lambda);
    for (i = 0; i < static_cast<int>(FreeEnergyPerturbationCouplingType::Count); i++)
    {
        for (j = 0; j < std::min(fep1->n_lambda, fep2->n_lambda); j++)
        {
            cmp_double(fp,
                       "inputrec->fepvals->all_lambda",
                       -1,
                       fep1->all_lambda[i][j],
                       fep2->all_lambda[i][j],
                       ftol,
                       abstol);
        }
    }
    cmp_int(fp, "inputrec->fepvals->lambda_neighbors", 1, fep1->lambda_neighbors, fep2->lambda_neighbors);
    cmp_real(fp, "inputrec->fepvals->sc_alpha", -1, fep1->sc_alpha, fep2->sc_alpha, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->sc_power", -1, fep1->sc_power, fep2->sc_power);
    cmp_real(fp, "inputrec->fepvals->sc_r_power", -1, fep1->sc_r_power, fep2->sc_r_power, ftol, abstol);
    cmp_real(fp, "inputrec->fepvals->sc_sigma", -1, fep1->sc_sigma, fep2->sc_sigma, ftol, abstol);
    cmpEnum(fp, "inputrec->fepvals->edHdLPrintEnergy", fep1->edHdLPrintEnergy, fep1->edHdLPrintEnergy);
    cmp_bool(fp, "inputrec->fepvals->bScCoul", -1, fep1->bScCoul, fep1->bScCoul);
    cmpEnum(fp, "inputrec->separate_dhdl_file", fep1->separate_dhdl_file, fep2->separate_dhdl_file);
    cmpEnum(fp, "inputrec->dhdl_derivatives", fep1->dhdl_derivatives, fep2->dhdl_derivatives);
    cmp_int(fp, "inputrec->dh_hist_size", -1, fep1->dh_hist_size, fep2->dh_hist_size);
    cmp_double(fp, "inputrec->dh_hist_spacing", -1, fep1->dh_hist_spacing, fep2->dh_hist_spacing, ftol, abstol);
    cmpEnum(fp, "inputrec->fepvals->softcoreFunction", fep1->softcoreFunction, fep2->softcoreFunction);
    cmp_real(fp,
             "inputrec->fepvals->scGapsysScaleLinpointLJ",
             -1,
             fep1->scGapsysScaleLinpointLJ,
             fep2->scGapsysScaleLinpointLJ,
             ftol,
             abstol);
    cmp_real(fp,
             "inputrec->fepvals->scGapsysScaleLinpointQ",
             -1,
             fep1->scGapsysScaleLinpointQ,
             fep2->scGapsysScaleLinpointQ,
             ftol,
             abstol);
    cmp_real(fp, "inputrec->fepvals->scGapsysSigmaLJ", -1, fep1->scGapsysSigmaLJ, fep2->scGapsysSigmaLJ, ftol, abstol);
}

void cmp_inputrec(FILE* fp, const t_inputrec* ir1, const t_inputrec* ir2, real ftol, real abstol)
{
    fprintf(fp, "comparing inputrec\n");

    /* gcc 2.96 doesnt like these defines at all, but issues a huge list
     * of warnings. Maybe it will change in future versions, but for the
     * moment I've spelled them out instead. /EL 000820
     * #define CIB(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
     * #define CII(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
     * #define CIR(s) cmp_real(fp,"inputrec->"#s,0,ir1->##s,ir2->##s,ftol)
     */
    cmpEnum(fp, "inputrec->eI", ir1->eI, ir2->eI);
    cmp_int64(fp, "inputrec->nsteps", ir1->nsteps, ir2->nsteps);
    cmp_int64(fp, "inputrec->init_step", ir1->init_step, ir2->init_step);
    cmp_int(fp, "inputrec->simulation_part", -1, ir1->simulation_part, ir2->simulation_part);
    cmp_int(fp, "inputrec->mts", -1, static_cast<int>(ir1->useMts), static_cast<int>(ir2->useMts));
    if (ir1->useMts && ir2->useMts)
    {
        cmp_int(fp, "inputrec->mts-levels", -1, ir1->mtsLevels.size(), ir2->mtsLevels.size());
        cmp_int(fp,
                "inputrec->mts-level2-forces",
                -1,
                ir1->mtsLevels[1].forceGroups.to_ulong(),
                ir2->mtsLevels[1].forceGroups.to_ulong());
        cmp_int(fp,
                "inputrec->mts-level2-factor",
                -1,
                ir1->mtsLevels[1].stepFactor,
                ir2->mtsLevels[1].stepFactor);
    }
    cmp_real(fp, "inputrec->massRepartitionFactor", -1, ir1->massRepartitionFactor, ir2->massRepartitionFactor, ftol, abstol);
    cmp_int(fp, "inputrec->pbcType", -1, static_cast<int>(ir1->pbcType), static_cast<int>(ir2->pbcType));
    cmp_bool(fp, "inputrec->bPeriodicMols", -1, ir1->bPeriodicMols, ir2->bPeriodicMols);
    cmpEnum(fp, "inputrec->cutoff_scheme", ir1->cutoff_scheme, ir2->cutoff_scheme);
    cmp_int(fp, "inputrec->nstlist", -1, ir1->nstlist, ir2->nstlist);
    cmp_int(fp, "inputrec->nstcomm", -1, ir1->nstcomm, ir2->nstcomm);
    cmpEnum(fp, "inputrec->comm_mode", ir1->comm_mode, ir2->comm_mode);
    cmp_int(fp, "inputrec->nstlog", -1, ir1->nstlog, ir2->nstlog);
    cmp_int(fp, "inputrec->nstxout", -1, ir1->nstxout, ir2->nstxout);
    cmp_int(fp, "inputrec->nstvout", -1, ir1->nstvout, ir2->nstvout);
    cmp_int(fp, "inputrec->nstfout", -1, ir1->nstfout, ir2->nstfout);
    cmp_int(fp, "inputrec->nstcalcenergy", -1, ir1->nstcalcenergy, ir2->nstcalcenergy);
    cmp_int(fp, "inputrec->nstenergy", -1, ir1->nstenergy, ir2->nstenergy);
    cmp_int(fp, "inputrec->nstxout_compressed", -1, ir1->nstxout_compressed, ir2->nstxout_compressed);
    cmp_double(fp, "inputrec->init_t", -1, ir1->init_t, ir2->init_t, ftol, abstol);
    cmp_double(fp, "inputrec->delta_t", -1, ir1->delta_t, ir2->delta_t, ftol, abstol);
    cmp_real(fp,
             "inputrec->x_compression_precision",
             -1,
             ir1->x_compression_precision,
             ir2->x_compression_precision,
             ftol,
             abstol);
    cmp_real(fp, "inputrec->fourierspacing", -1, ir1->fourier_spacing, ir2->fourier_spacing, ftol, abstol);
    cmp_int(fp, "inputrec->nkx", -1, ir1->nkx, ir2->nkx);
    cmp_int(fp, "inputrec->nky", -1, ir1->nky, ir2->nky);
    cmp_int(fp, "inputrec->nkz", -1, ir1->nkz, ir2->nkz);
    cmp_int(fp, "inputrec->pme_order", -1, ir1->pme_order, ir2->pme_order);
    cmp_real(fp, "inputrec->ewald_rtol", -1, ir1->ewald_rtol, ir2->ewald_rtol, ftol, abstol);
    cmpEnum(fp, "inputrec->ewald_geometry", ir1->ewald_geometry, ir2->ewald_geometry);
    cmp_real(fp, "inputrec->epsilon_surface", -1, ir1->epsilon_surface, ir2->epsilon_surface, ftol, abstol);
    cmp_int(fp,
            "inputrec->bContinuation",
            -1,
            static_cast<int>(ir1->bContinuation),
            static_cast<int>(ir2->bContinuation));
    cmp_int(fp, "inputrec->bShakeSOR", -1, static_cast<int>(ir1->bShakeSOR), static_cast<int>(ir2->bShakeSOR));
    cmpEnum(fp, "inputrec->ensembleTemperatureSetting", ir1->ensembleTemperatureSetting, ir2->ensembleTemperatureSetting);
    if (ir1->ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant
        && ir2->ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant)
    {
        cmp_real(fp, "inputrec->ensembleTemperature", -1, ir1->ensembleTemperature, ir2->ensembleTemperature, ftol, abstol);
    }
    cmpEnum(fp, "inputrec->etc", ir1->etc, ir2->etc);
    cmp_int(fp,
            "inputrec->bPrintNHChains",
            -1,
            static_cast<int>(ir1->bPrintNHChains),
            static_cast<int>(ir2->bPrintNHChains));
    cmpEnum(fp,
            "inputrec->pressureCouplingOptions.epc",
            ir1->pressureCouplingOptions.epc,
            ir2->pressureCouplingOptions.epc);
    if (ir1->pressureCouplingOptions.epc != PressureCoupling::No
        && ir2->pressureCouplingOptions.epc != PressureCoupling::No)
    {
        cmpPressureCouplingOptions(
                fp, ir1->pressureCouplingOptions, ir2->pressureCouplingOptions, ftol, abstol);
    }
    // Refcoord-scaling is also needed for other algorithms that affect the box
    cmpEnum(fp,
            "refcoord_scaling",
            ir1->pressureCouplingOptions.refcoord_scaling,
            ir2->pressureCouplingOptions.refcoord_scaling);
    cmp_rvec(fp, "inputrec->posres_com", -1, ir1->posres_com, ir2->posres_com, ftol, abstol);
    cmp_rvec(fp, "inputrec->posres_comB", -1, ir1->posres_comB, ir2->posres_comB, ftol, abstol);
    cmp_real(fp, "inputrec->verletbuf_tol", -1, ir1->verletbuf_tol, ir2->verletbuf_tol, ftol, abstol);
    cmp_real(fp,
             "inputrec->verlet-buffer-pressure-tolerance",
             -1,
             ir1->verletBufferPressureTolerance,
             ir2->verletBufferPressureTolerance,
             ftol,
             abstol);
    cmp_real(fp, "inputrec->rlist", -1, ir1->rlist, ir2->rlist, ftol, abstol);
    cmp_real(fp, "inputrec->rtpi", -1, ir1->rtpi, ir2->rtpi, ftol, abstol);
    cmpEnum(fp, "inputrec->coulombtype", ir1->coulombtype, ir2->coulombtype);
    cmpEnum(fp, "inputrec->coulomb_modifier", ir1->coulomb_modifier, ir2->coulomb_modifier);
    cmp_real(fp, "inputrec->rcoulomb_switch", -1, ir1->rcoulomb_switch, ir2->rcoulomb_switch, ftol, abstol);
    cmp_real(fp, "inputrec->rcoulomb", -1, ir1->rcoulomb, ir2->rcoulomb, ftol, abstol);
    cmpEnum(fp, "inputrec->vdwtype", ir1->vdwtype, ir2->vdwtype);
    cmpEnum(fp, "inputrec->vdw_modifier", ir1->vdw_modifier, ir2->vdw_modifier);
    cmp_real(fp, "inputrec->rvdw_switch", -1, ir1->rvdw_switch, ir2->rvdw_switch, ftol, abstol);
    cmp_real(fp, "inputrec->rvdw", -1, ir1->rvdw, ir2->rvdw, ftol, abstol);
    cmp_real(fp, "inputrec->epsilon_r", -1, ir1->epsilon_r, ir2->epsilon_r, ftol, abstol);
    cmp_real(fp, "inputrec->epsilon_rf", -1, ir1->epsilon_rf, ir2->epsilon_rf, ftol, abstol);
    cmp_real(fp, "inputrec->tabext", -1, ir1->tabext, ir2->tabext, ftol, abstol);

    cmpEnum(fp, "inputrec->eDispCorr", ir1->eDispCorr, ir2->eDispCorr);
    cmp_real(fp, "inputrec->shake_tol", -1, ir1->shake_tol, ir2->shake_tol, ftol, abstol);
    cmpEnum(fp, "inputrec->efep", ir1->efep, ir2->efep);
    cmp_fepvals(fp, ir1->fepvals.get(), ir2->fepvals.get(), ftol, abstol);
    cmp_int(fp, "inputrec->bSimTemp", -1, static_cast<int>(ir1->bSimTemp), static_cast<int>(ir2->bSimTemp));
    if ((ir1->bSimTemp == ir2->bSimTemp) && (ir1->bSimTemp))
    {
        cmp_simtempvals(fp,
                        ir1->simtempvals.get(),
                        ir2->simtempvals.get(),
                        std::min(ir1->fepvals->n_lambda, ir2->fepvals->n_lambda),
                        ftol,
                        abstol);
    }
    cmp_int(fp, "inputrec->bExpanded", -1, static_cast<int>(ir1->bExpanded), static_cast<int>(ir2->bExpanded));
    if ((ir1->bExpanded == ir2->bExpanded) && (ir1->bExpanded))
    {
        cmp_expandedvals(fp,
                         ir1->expandedvals.get(),
                         ir2->expandedvals.get(),
                         std::min(ir1->fepvals->n_lambda, ir2->fepvals->n_lambda),
                         ftol,
                         abstol);
    }
    cmp_int(fp, "inputrec->nwall", -1, ir1->nwall, ir2->nwall);
    cmpEnum(fp, "inputrec->wall_type", ir1->wall_type, ir2->wall_type);
    cmp_int(fp, "inputrec->wall_atomtype[0]", -1, ir1->wall_atomtype[0], ir2->wall_atomtype[0]);
    cmp_int(fp, "inputrec->wall_atomtype[1]", -1, ir1->wall_atomtype[1], ir2->wall_atomtype[1]);
    cmp_real(fp, "inputrec->wall_density[0]", -1, ir1->wall_density[0], ir2->wall_density[0], ftol, abstol);
    cmp_real(fp, "inputrec->wall_density[1]", -1, ir1->wall_density[1], ir2->wall_density[1], ftol, abstol);
    cmp_real(fp, "inputrec->wall_ewald_zfac", -1, ir1->wall_ewald_zfac, ir2->wall_ewald_zfac, ftol, abstol);

    cmp_bool(fp, "inputrec->bPull", -1, ir1->bPull, ir2->bPull);
    if (ir1->bPull && ir2->bPull)
    {
        cmp_pull(fp);
    }

    cmp_bool(fp, "inputrec->bDoAwh", -1, ir1->bDoAwh, ir2->bDoAwh);
    if (ir1->bDoAwh && ir2->bDoAwh)
    {
        cmp_awhParams(fp, *ir1->awhParams, *ir2->awhParams, ftol, abstol);
    }

    cmpEnum(fp, "inputrec->eDisre", ir1->eDisre, ir2->eDisre);
    cmp_real(fp, "inputrec->dr_fc", -1, ir1->dr_fc, ir2->dr_fc, ftol, abstol);
    cmpEnum(fp, "inputrec->eDisreWeighting", ir1->eDisreWeighting, ir2->eDisreWeighting);
    cmp_int(fp, "inputrec->bDisreMixed", -1, static_cast<int>(ir1->bDisreMixed), static_cast<int>(ir2->bDisreMixed));
    cmp_int(fp, "inputrec->nstdisreout", -1, ir1->nstdisreout, ir2->nstdisreout);
    cmp_real(fp, "inputrec->dr_tau", -1, ir1->dr_tau, ir2->dr_tau, ftol, abstol);
    cmp_real(fp, "inputrec->orires_fc", -1, ir1->orires_fc, ir2->orires_fc, ftol, abstol);
    cmp_real(fp, "inputrec->orires_tau", -1, ir1->orires_tau, ir2->orires_tau, ftol, abstol);
    cmp_int(fp, "inputrec->nstorireout", -1, ir1->nstorireout, ir2->nstorireout);
    cmp_real(fp, "inputrec->em_stepsize", -1, ir1->em_stepsize, ir2->em_stepsize, ftol, abstol);
    cmp_real(fp, "inputrec->em_tol", -1, ir1->em_tol, ir2->em_tol, ftol, abstol);
    cmp_int(fp, "inputrec->niter", -1, ir1->niter, ir2->niter);
    cmp_real(fp, "inputrec->fc_stepsize", -1, ir1->fc_stepsize, ir2->fc_stepsize, ftol, abstol);
    cmp_int(fp, "inputrec->nstcgsteep", -1, ir1->nstcgsteep, ir2->nstcgsteep);
    cmp_int(fp, "inputrec->nbfgscorr", 0, ir1->nbfgscorr, ir2->nbfgscorr);
    cmpEnum(fp, "inputrec->eConstrAlg", ir1->eConstrAlg, ir2->eConstrAlg);
    cmp_int(fp, "inputrec->nProjOrder", -1, ir1->nProjOrder, ir2->nProjOrder);
    cmp_real(fp, "inputrec->LincsWarnAngle", -1, ir1->LincsWarnAngle, ir2->LincsWarnAngle, ftol, abstol);
    cmp_int(fp, "inputrec->nLincsIter", -1, ir1->nLincsIter, ir2->nLincsIter);
    cmp_real(fp, "inputrec->bd_fric", -1, ir1->bd_fric, ir2->bd_fric, ftol, abstol);
    cmp_int64(fp, "inputrec->ld_seed", ir1->ld_seed, ir2->ld_seed);
    cmp_real(fp, "inputrec->cos_accel", -1, ir1->cos_accel, ir2->cos_accel, ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(a)", -1, ir1->deform[XX], ir2->deform[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(b)", -1, ir1->deform[YY], ir2->deform[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(c)", -1, ir1->deform[ZZ], ir2->deform[ZZ], ftol, abstol);


    cmp_int(fp, "inputrec->userint1", -1, ir1->userint1, ir2->userint1);
    cmp_int(fp, "inputrec->userint2", -1, ir1->userint2, ir2->userint2);
    cmp_int(fp, "inputrec->userint3", -1, ir1->userint3, ir2->userint3);
    cmp_int(fp, "inputrec->userint4", -1, ir1->userint4, ir2->userint4);
    cmp_real(fp, "inputrec->userreal1", -1, ir1->userreal1, ir2->userreal1, ftol, abstol);
    cmp_real(fp, "inputrec->userreal2", -1, ir1->userreal2, ir2->userreal2, ftol, abstol);
    cmp_real(fp, "inputrec->userreal3", -1, ir1->userreal3, ir2->userreal3, ftol, abstol);
    cmp_real(fp, "inputrec->userreal4", -1, ir1->userreal4, ir2->userreal4, ftol, abstol);
    cmp_grpopts(fp, &(ir1->opts), &(ir2->opts), ftol, abstol);
    gmx::TextWriter writer(fp);
    gmx::compareKeyValueTrees(&writer, *ir1->params, *ir2->params, ftol, abstol);
}

void comp_pull_AB(FILE* fp, const pull_params_t& pull, real ftol, real abstol)
{
    for (int i = 0; i < pull.ncoord; i++)
    {
        fprintf(fp, "comparing pull coord %d\n", i);
        cmp_real(fp, "pull-coord->k", -1, pull.coord[i].k, pull.coord[i].kB, ftol, abstol);
    }
}

gmx_bool inputrecDeform(const t_inputrec* ir)
{
    return (ir->deform[XX][XX] != 0 || ir->deform[YY][YY] != 0 || ir->deform[ZZ][ZZ] != 0
            || ir->deform[YY][XX] != 0 || ir->deform[ZZ][XX] != 0 || ir->deform[ZZ][YY] != 0);
}

gmx_bool inputrecDynamicBox(const t_inputrec* ir)
{
    return (ir->pressureCouplingOptions.epc != PressureCoupling::No
            || ir->eI == IntegrationAlgorithm::TPI || inputrecDeform(ir));
}

bool shouldPreserveBoxShape(const PressureCouplingOptions& pressureCouplingOptions, const tensor deform)
{
    return (pressureCouplingOptions.epc != PressureCoupling::No && deform[XX][XX] == 0
            && (pressureCouplingOptions.epct == PressureCouplingType::Isotropic
                || pressureCouplingOptions.epct == PressureCouplingType::SemiIsotropic));
}

gmx_bool inputrecNeedMutot(const t_inputrec* ir)
{
    return ((ir->coulombtype == CoulombInteractionType::Ewald || usingPme(ir->coulombtype))
            && (ir->ewald_geometry == EwaldGeometry::ThreeDC || ir->epsilon_surface != 0));
}

gmx_bool inputrecExclForces(const t_inputrec* ir)
{
    return (usingFullElectrostatics(ir->coulombtype) || (usingRF(ir->coulombtype)));
}

gmx_bool inputrecNptTrotter(const t_inputrec* ir)
{
    return (((ir->eI == IntegrationAlgorithm::VV) || (ir->eI == IntegrationAlgorithm::VVAK))
            && (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
            && (ir->etc == TemperatureCoupling::NoseHoover));
}

gmx_bool inputrecNvtTrotter(const t_inputrec* ir)
{
    return (((ir->eI == IntegrationAlgorithm::VV) || (ir->eI == IntegrationAlgorithm::VVAK))
            && (ir->pressureCouplingOptions.epc != PressureCoupling::Mttk)
            && (ir->etc == TemperatureCoupling::NoseHoover));
}

gmx_bool inputrecNphTrotter(const t_inputrec* ir)
{
    return (((ir->eI == IntegrationAlgorithm::VV) || (ir->eI == IntegrationAlgorithm::VVAK))
            && (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
            && (ir->etc != TemperatureCoupling::NoseHoover));
}

bool inputrecPbcXY2Walls(const t_inputrec* ir)
{
    return (ir->pbcType == PbcType::XY && ir->nwall == 2);
}

bool inputrecFrozenAtoms(const t_inputrec* ir)
{
    return ((ir->opts.nFreeze != nullptr)
            && (ir->opts.ngfrz > 1 || ir->opts.nFreeze[0][XX] != 0 || ir->opts.nFreeze[0][YY] != 0
                || ir->opts.nFreeze[0][ZZ] != 0));
}

//! Asserts that we have a valid, processed ensemble temperature setting
static void assertValidEnsembleTemperatureSetting(const EnsembleTemperatureSetting rts)
{
    GMX_RELEASE_ASSERT(rts == EnsembleTemperatureSetting::NotAvailable
                               || rts == EnsembleTemperatureSetting::Constant
                               || rts == EnsembleTemperatureSetting::Variable,
                       "At this point in the code the ensemble temperature setting cannot be auto");
}

bool haveConstantEnsembleTemperature(const t_inputrec& ir)
{
    assertValidEnsembleTemperatureSetting(ir.ensembleTemperatureSetting);

    return ir.ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant;
}

real constantEnsembleTemperature(const t_inputrec& ir)
{
    GMX_RELEASE_ASSERT(ir.ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant,
                       "Can only request for a constant ensemble temperature when available");

    return ir.ensembleTemperature;
}

bool haveEnsembleTemperature(const t_inputrec& ir)
{
    assertValidEnsembleTemperatureSetting(ir.ensembleTemperatureSetting);

    return (ir.ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant
            || ir.ensembleTemperatureSetting == EnsembleTemperatureSetting::Variable);
}

bool integratorHasConservedEnergyQuantity(const t_inputrec* ir)
{
    if (!EI_MD(ir->eI))
    { // NOLINT bugprone-branch-clone
        // Energy minimization or stochastic integrator: no conservation
        return false;
    }
    else if (ir->etc == TemperatureCoupling::No && ir->pressureCouplingOptions.epc == PressureCoupling::No)
    {
        // The total energy is conserved, no additional conserved quanitity
        return false;
    }
    else
    {
        // Shear stress with Parrinello-Rahman is not supported (tedious)
        bool shearWithPR = ((ir->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
                             || ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
                            && (ir->pressureCouplingOptions.ref_p[YY][XX] != 0
                                || ir->pressureCouplingOptions.ref_p[ZZ][XX] != 0
                                || ir->pressureCouplingOptions.ref_p[ZZ][YY] != 0));

        return !ETC_ANDERSEN(ir->etc) && !shearWithPR;
    }
}

bool integratorHasReferenceTemperature(const t_inputrec& ir)
{
    return (ir.etc != TemperatureCoupling::No || EI_SD(ir.eI) || ir.eI == IntegrationAlgorithm::BD
            || EI_TPI(ir.eI));
}

bool doSimulatedAnnealing(const t_inputrec& ir)
{
    for (int i = 0; i < ir.opts.ngtc; i++)
    {
        if (ir.opts.annealing[i] != SimulatedAnnealing::No)
        {
            return true;
        }
    }
    return false;
}

int inputrec2nboundeddim(const t_inputrec* ir)
{
    if (inputrecPbcXY2Walls(ir))
    {
        return 3;
    }
    else
    {
        return numPbcDimensions(ir->pbcType);
    }
}

int ndof_com(const t_inputrec* ir)
{
    int n = 0;

    switch (ir->pbcType)
    {
        case PbcType::Xyz:
        case PbcType::No: n = 3; break;
        case PbcType::XY: n = (ir->nwall == 0 ? 3 : 2); break;
        case PbcType::Screw: n = 1; break;
        default: gmx_incons("Unknown pbc in calc_nrdf");
    }

    return n;
}

real maxReferenceTemperature(const t_inputrec& ir)
{
    if (EI_ENERGY_MINIMIZATION(ir.eI) || ir.eI == IntegrationAlgorithm::NM)
    {
        return 0;
    }

    if (EI_MD(ir.eI) && ir.etc == TemperatureCoupling::No)
    {
        return -1;
    }

    /* SD and BD also use ref_t and tau_t for setting the reference temperature.
     * TPI can be treated as MD, since it needs an ensemble temperature.
     */

    real maxTemperature = 0;
    for (int i = 0; i < ir.opts.ngtc; i++)
    {
        if (ir.opts.tau_t[i] >= 0)
        {
            maxTemperature = std::max(maxTemperature, ir.opts.ref_t[i]);
        }
    }

    return maxTemperature;
}

bool haveEwaldSurfaceContribution(const t_inputrec& ir)
{
    return usingPmeOrEwald(ir.coulombtype)
           && (ir.ewald_geometry == EwaldGeometry::ThreeDC || ir.epsilon_surface != 0);
}

bool haveFreeEnergyType(const t_inputrec& ir, const int fepType)
{
    for (int i = 0; i < ir.fepvals->n_lambda; i++)
    {
        if (ir.fepvals->all_lambda[fepType][i] > 0)
        {
            return true;
        }
    }
    return false;
}

bool fepLambdasChangeAtSameRate(
        const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>>& allLambdas)
{
    const int numLambdas = gmx::ssize(allLambdas[FreeEnergyPerturbationCouplingType::Fep]);

    // Check all lambda intervals
    for (int i = 0; i < numLambdas - 1; i++)
    {
        double previousDiff = 0;

        // Check if all lambda components have the same or zero difference
        for (const auto& lambdas : allLambdas)
        {
            // The dhdl code takes the sign into account, so we allow
            // different components to differ by a factor -1
            const double diff = std::abs(lambdas[i + 1] - lambdas[i]);

            // Lambda values are usually of order 1, so here we can,
            // by exception, compare numbers directly to float epsilon
            if (diff > GMX_FLOAT_EPS)
            {
                if (previousDiff > GMX_FLOAT_EPS && std::abs(diff - previousDiff) > GMX_FLOAT_EPS)
                {
                    return false;
                }

                previousDiff = diff;
            }
        }
    }

    return true;
}
