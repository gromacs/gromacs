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

#include "readir.h"

#include <cctype>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>

#include "gromacs/applied_forces/awh/read_params.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/seed.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/ikeyvaluetreeerror.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#define NOGID 255

using gmx::BasicVector;

/* Resource parameters
 * Do not change any of these until you read the instruction
 * in readinp.h. Some cpp's do not take spaces after the backslash
 * (like the c-shell), which will give you a very weird compiler
 * message.
 */

struct gmx_inputrec_strings
{
    char tcgrps[STRLEN], tau_t[STRLEN], ref_t[STRLEN], accelerationGroups[STRLEN],
            acceleration[STRLEN], freeze[STRLEN], frdim[STRLEN], energy[STRLEN], user1[STRLEN],
            user2[STRLEN], vcm[STRLEN], x_compressed_groups[STRLEN], couple_moltype[STRLEN],
            orirefitgrp[STRLEN], egptable[STRLEN], egpexcl[STRLEN], wall_atomtype[STRLEN],
            wall_density[STRLEN], deform[STRLEN], QMMM[STRLEN], imd_grp[STRLEN];
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::string> fep_lambda;
    char                                                                   lambda_weights[STRLEN];
    std::vector<std::string>                                               pullGroupNames;
    std::vector<std::string>                                               rotateGroupNames;
    char anneal[STRLEN], anneal_npoints[STRLEN], anneal_time[STRLEN], anneal_temp[STRLEN];
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static gmx_inputrec_strings* inputrecStrings = nullptr;

void init_inputrec_strings()
{
    if (inputrecStrings)
    {
        gmx_incons(
                "Attempted to call init_inputrec_strings before calling done_inputrec_strings. "
                "Only one inputrec (i.e. .mdp file) can be parsed at a time.");
    }
    inputrecStrings = new gmx_inputrec_strings();
}

void done_inputrec_strings()
{
    delete inputrecStrings;
    inputrecStrings = nullptr;
}


//! How to treat coverage of the whole system for a set of atom groupsx
enum class GroupCoverage
{
    All,             //!< All particles have to be a member of a group
    AllGenerateRest, //<! A rest group with name is generated for particles not part of any group
    Partial,         //<! As \p AllGenerateRest, but no name for the rest group is generated
    OneGroup //<! Merge all selected groups into one group, make a rest group for the remaining particles
};

static const char* const constraints[eshNR + 1] = { "none",     "h-bonds",    "all-bonds",
                                                    "h-angles", "all-angles", nullptr };

static const char* const couple_lam[ecouplamNR + 1] = { "vdw-q", "vdw", "q", "none", nullptr };

static void getSimTemps(int ntemps, t_simtemp* simtemp, gmx::ArrayRef<double> temperature_lambdas)
{

    int i;

    for (i = 0; i < ntemps; i++)
    {
        /* simple linear scaling -- allows more control */
        if (simtemp->eSimTempScale == SimulatedTempering::Linear)
        {
            simtemp->temperatures[i] =
                    simtemp->simtemp_low
                    + (simtemp->simtemp_high - simtemp->simtemp_low) * temperature_lambdas[i];
        }
        else if (simtemp->eSimTempScale
                 == SimulatedTempering::Geometric) /* should give roughly equal acceptance for constant heat capacity . . . */
        {
            simtemp->temperatures[i] = simtemp->simtemp_low
                                       * std::pow(simtemp->simtemp_high / simtemp->simtemp_low,
                                                  static_cast<real>((1.0 * i) / (ntemps - 1)));
        }
        else if (simtemp->eSimTempScale == SimulatedTempering::Exponential)
        {
            simtemp->temperatures[i] = simtemp->simtemp_low
                                       + (simtemp->simtemp_high - simtemp->simtemp_low)
                                                 * (std::expm1(temperature_lambdas[i]) / std::expm1(1.0));
        }
        else
        {
            char errorstr[128];
            sprintf(errorstr, "eSimTempScale=%s not defined", enumValueToString(simtemp->eSimTempScale));
            gmx_fatal(FARGS, "%s", errorstr);
        }
    }
}


static void _low_check(bool b, const char* s, WarningHandler* wi)
{
    if (b)
    {
        wi->addError(s);
    }
}

static void check_nst(const char* desc_nst, int nst, const char* desc_p, int* p, WarningHandler* wi)
{
    char buf[STRLEN];

    if (*p > 0 && *p % nst != 0)
    {
        /* Round up to the next multiple of nst */
        *p = ((*p) / nst + 1) * nst;
        sprintf(buf, "%s should be a multiple of %s, changing %s to %d\n", desc_p, desc_nst, desc_p, *p);
        wi->addWarning(buf);
    }
}

//! Convert legacy mdp entries to modern ones.
static void process_interaction_modifier(InteractionModifiers* eintmod)
{
    if (*eintmod == InteractionModifiers::PotShiftVerletUnsupported)
    {
        *eintmod = InteractionModifiers::PotShift;
    }
}

void check_ir(const char*                    mdparin,
              const gmx::MDModulesNotifiers& mdModulesNotifiers,
              t_inputrec*                    ir,
              t_gromppopts*                  opts,
              WarningHandler*                wi)
/* Check internal consistency.
 * NOTE: index groups are not set here yet, don't check things
 * like temperature coupling group options here, but in triple_check
 */
{
    /* Strange macro: first one fills the err_buf, and then one can check
     * the condition, which will print the message and increase the error
     * counter.
     */
#define CHECK(b) _low_check(b, err_buf, wi)
    char        err_buf[256], warn_buf[STRLEN];
    int         i, j;
    real        dt_pcoupl;
    t_lambda*   fep    = ir->fepvals.get();
    t_expanded* expand = ir->expandedvals.get();

    wi->setFileAndLineNumber(mdparin, -1);

    if (EI_DYNAMICS(ir->eI) && ir->delta_t <= 0)
    {
        wi->addError("dt should be larger than 0");
    }

    /* We cannot check MTS requirements with an invalid MTS setup
     * and we will already have generated errors with an invalid MTS setup.
     */
    if (gmx::haveValidMtsSetup(*ir))
    {
        std::vector<std::string> errorMessages = gmx::checkMtsRequirements(*ir);

        for (const auto& errorMessage : errorMessages)
        {
            wi->addError(errorMessage);
        }
    }

    if (ir->coulombtype == CoulombInteractionType::RFNecUnsupported)
    {
        std::string message =
                gmx::formatString("%s electrostatics is no longer supported",
                                  enumValueToString(CoulombInteractionType::RFNecUnsupported));
        wi->addError(message);
    }

    /* BASIC CUT-OFF STUFF */
    if (ir->rcoulomb < 0)
    {
        wi->addError("rcoulomb should be >= 0");
    }
    if (ir->rvdw < 0)
    {
        wi->addError("rvdw should be >= 0");
    }
    if (ir->rlist < 0 && !(ir->cutoff_scheme == CutoffScheme::Verlet && ir->verletbuf_tol > 0))
    {
        wi->addError("rlist should be >= 0");
    }
    sprintf(err_buf,
            "nstlist can not be smaller than 0. (If you were trying to use the heuristic "
            "neighbour-list update scheme for efficient buffering for improved energy "
            "conservation, please use the Verlet cut-off scheme instead.)");
    CHECK(ir->nstlist < 0);

    process_interaction_modifier(&ir->coulomb_modifier);
    process_interaction_modifier(&ir->vdw_modifier);

    if (ir->cutoff_scheme == CutoffScheme::Group)
    {
        gmx_fatal(FARGS,
                  "The group cutoff scheme has been removed since GROMACS 2020. "
                  "Please use the Verlet cutoff scheme.");
    }
    if (ir->cutoff_scheme == CutoffScheme::Verlet)
    {
        real rc_max;

        /* Normal Verlet type neighbor-list, currently only limited feature support */
        if (inputrec2nboundeddim(ir) < 3)
        {
            wi->addError("With Verlet lists only full pbc or pbc=xy with walls is supported");
        }

        // We don't (yet) have general Verlet kernels for rcoulomb!=rvdw
        if (ir->rcoulomb != ir->rvdw)
        {
            // Since we have PME coulomb + LJ cut-off kernels with rcoulomb>rvdw
            // for PME load balancing, we can support this exception.
            bool bUsesPmeTwinRangeKernel =
                    (usingPmeOrEwald(ir->coulombtype) && ir->vdwtype == VanDerWaalsType::Cut
                     && ir->rcoulomb > ir->rvdw);
            if (!bUsesPmeTwinRangeKernel)
            {
                wi->addError(
                        "With Verlet lists rcoulomb!=rvdw is not supported (except for "
                        "rcoulomb>rvdw with PME electrostatics)");
            }
        }

        if (ir->vdwtype == VanDerWaalsType::Shift || ir->vdwtype == VanDerWaalsType::Switch)
        {
            if (ir->vdw_modifier == InteractionModifiers::None
                || ir->vdw_modifier == InteractionModifiers::PotShift)
            {
                ir->vdw_modifier =
                        (ir->vdwtype == VanDerWaalsType::Shift ? InteractionModifiers::ForceSwitch
                                                               : InteractionModifiers::PotSwitch);

                sprintf(warn_buf,
                        "Replacing vdwtype=%s by the equivalent combination of vdwtype=%s and "
                        "vdw_modifier=%s",
                        enumValueToString(ir->vdwtype),
                        enumValueToString(VanDerWaalsType::Cut),
                        enumValueToString(ir->vdw_modifier));
                wi->addNote(warn_buf);

                ir->vdwtype = VanDerWaalsType::Cut;
            }
            else
            {
                sprintf(warn_buf,
                        "Unsupported combination of vdwtype=%s and vdw_modifier=%s",
                        enumValueToString(ir->vdwtype),
                        enumValueToString(ir->vdw_modifier));
                wi->addError(warn_buf);
            }
        }

        if (!(ir->vdwtype == VanDerWaalsType::Cut || ir->vdwtype == VanDerWaalsType::Pme))
        {
            wi->addError("With Verlet lists only cut-off and PME LJ interactions are supported");
        }
        if (!(ir->coulombtype == CoulombInteractionType::Cut || usingRF(ir->coulombtype)
              || usingPme(ir->coulombtype) || ir->coulombtype == CoulombInteractionType::Ewald))
        {
            wi->addError(
                    "With Verlet lists only cut-off, reaction-field, PME and Ewald "
                    "electrostatics are supported");
        }
        if (!(ir->coulomb_modifier == InteractionModifiers::None
              || ir->coulomb_modifier == InteractionModifiers::PotShift))
        {
            sprintf(warn_buf, "coulomb_modifier=%s is not supported", enumValueToString(ir->coulomb_modifier));
            wi->addError(warn_buf);
        }

        if (usingUserTableElectrostatics(ir->coulombtype))
        {
            sprintf(warn_buf,
                    "Coulomb type %s is not supported with the verlet scheme",
                    enumValueToString(ir->coulombtype));
            wi->addError(warn_buf);
        }

        if (ir->nstlist <= 0)
        {
            wi->addError("With Verlet lists nstlist should be larger than 0");
        }

        if (ir->nstlist < 10)
        {
            wi->addNote(
                    "With Verlet lists the optimal nstlist is >= 10, with GPUs >= 20. Note "
                    "that with the Verlet scheme, nstlist has no effect on the accuracy of "
                    "your simulation.");
        }

        rc_max = std::max(ir->rvdw, ir->rcoulomb);

        if (EI_TPI(ir->eI))
        {
            /* With TPI we set the pairlist cut-off later using the radius of the insterted molecule */
            ir->verletbuf_tol                 = 0;
            ir->verletBufferPressureTolerance = 0;
            ir->rlist                         = rc_max;
        }
        else if (ir->verletbuf_tol <= 0)
        {
            if (ir->verletbuf_tol == 0)
            {
                wi->addError("Can not have Verlet buffer tolerance of exactly 0");
            }

            if (ir->rlist < rc_max)
            {
                wi->addError("With verlet lists rlist can not be smaller than rvdw or rcoulomb");
            }

            if (ir->rlist == rc_max && ir->nstlist > 1)
            {
                wi->addNote(

                        "rlist is equal to rvdw and/or rcoulomb: there is no explicit Verlet "
                        "buffer. The cluster pair list does have a buffering effect, but choosing "
                        "a larger rlist might be necessary for good energy conservation.");
            }

            if (ir->verletBufferPressureTolerance > 0)
            {
                if (ir->nstlist > 1)
                {
                    wi->addNote(
                            "verlet-buffer-pressure-tolerance is ignored when "
                            "verlet-buffer-tolerance < 0");
                }
                ir->verletBufferPressureTolerance = -1;
            }
        }
        else
        {
            if (ir->verletBufferPressureTolerance == 0)
            {
                wi->addError("verlet-buffer-pressure-tolerance cannot be exactly 0");
            }

            if (ir->rlist > rc_max)
            {
                wi->addNote(
                        "You have set rlist larger than the interaction cut-off, but you also "
                        "have verlet-buffer-tolerance > 0. Will set rlist using "
                        "verlet-buffer-tolerance.");
            }

            if (ir->nstlist == 1)
            {
                /* No buffer required */
                ir->rlist = rc_max;
            }
            else
            {
                if (EI_DYNAMICS(ir->eI))
                {
                    if (inputrec2nboundeddim(ir) < 3)
                    {
                        wi->addError(
                                "The box volume is required for calculating rlist from the "
                                "energy drift with verlet-buffer-tolerance > 0. You are "
                                "using at least one unbounded dimension, so no volume can be "
                                "computed. Either use a finite box, or set rlist yourself "
                                "together with verlet-buffer-tolerance = -1.");
                    }
                    /* Set rlist temporarily so we can continue processing */
                    ir->rlist = rc_max;
                }
                else
                {
                    /* Set the buffer to 5% of the cut-off */
                    ir->rlist = (1.0 + verlet_buffer_ratio_nodynamics) * rc_max;
                }
            }
        }
    }

    /* GENERAL INTEGRATOR STUFF */
    if (!EI_MD(ir->eI))
    {
        if (ir->etc != TemperatureCoupling::No)
        {
            if (EI_RANDOM(ir->eI))
            {
                sprintf(warn_buf,
                        "Setting tcoupl from '%s' to 'no'. %s handles temperature coupling "
                        "implicitly. See the documentation for more information on which "
                        "parameters affect temperature for %s.",
                        enumValueToString(ir->etc),
                        enumValueToString(ir->eI),
                        enumValueToString(ir->eI));
            }
            else
            {
                sprintf(warn_buf,
                        "Setting tcoupl from '%s' to 'no'. Temperature coupling does not apply to "
                        "%s.",
                        enumValueToString(ir->etc),
                        enumValueToString(ir->eI));
            }
            wi->addNote(warn_buf);
        }
        ir->etc = TemperatureCoupling::No;
    }
    if (ir->eI == IntegrationAlgorithm::VVAK)
    {
        sprintf(warn_buf,
                "Integrator method %s is implemented primarily for validation purposes; for "
                "molecular dynamics, you should probably be using %s or %s",
                enumValueToString(IntegrationAlgorithm::VVAK),
                enumValueToString(IntegrationAlgorithm::MD),
                enumValueToString(IntegrationAlgorithm::VV));
        wi->addNote(warn_buf);
    }
    if (!EI_DYNAMICS(ir->eI))
    {
        if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
        {
            sprintf(warn_buf,
                    "Setting pcoupl from '%s' to 'no'. Pressure coupling does not apply to %s.",
                    enumValueToString(ir->pressureCouplingOptions.epc),
                    enumValueToString(ir->eI));
            wi->addNote(warn_buf);
        }
        ir->pressureCouplingOptions.epc = PressureCoupling::No;
    }
    if (EI_DYNAMICS(ir->eI))
    {
        // Replace old -1 "automation" values by the default value of 100
        if (ir->nstcalcenergy < 0)
        {
            ir->nstcalcenergy = 100;
        }

        if ((ir->nstenergy > 0 && ir->nstcalcenergy > ir->nstenergy)
            || (ir->efep != FreeEnergyPerturbationType::No && ir->fepvals->nstdhdl > 0
                && (ir->nstcalcenergy > ir->fepvals->nstdhdl)))

        {
            const char* nsten    = "nstenergy";
            const char* nstdh    = "nstdhdl";
            const char* min_name = nsten;
            int         min_nst  = ir->nstenergy;

            /* find the smallest of ( nstenergy, nstdhdl ) */
            if (ir->efep != FreeEnergyPerturbationType::No && ir->fepvals->nstdhdl > 0
                && (ir->nstenergy == 0 || ir->fepvals->nstdhdl < ir->nstenergy))
            {
                min_nst  = ir->fepvals->nstdhdl;
                min_name = nstdh;
            }
            /* If the user sets nstenergy small, we should respect that */
            sprintf(warn_buf, "Setting nstcalcenergy (%d) equal to %s (%d)", ir->nstcalcenergy, min_name, min_nst);
            wi->addNote(warn_buf);
            ir->nstcalcenergy = min_nst;
        }

        if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
        {
            if (ir->pressureCouplingOptions.nstpcouple < 0)
            {
                ir->pressureCouplingOptions.nstpcouple = ir_optimal_nstpcouple(ir);
            }
            if (ir->useMts && ir->pressureCouplingOptions.nstpcouple % ir->mtsLevels.back().stepFactor != 0)
            {
                wi->addError(
                        "With multiple time stepping, nstpcouple should be a multiple of "
                        "mts-factor");
            }
        }

        if (ir->nstcalcenergy > 0)
        {
            if (ir->efep != FreeEnergyPerturbationType::No)
            {
                /* nstdhdl should be a multiple of nstcalcenergy */
                check_nst("nstcalcenergy", ir->nstcalcenergy, "nstdhdl", &ir->fepvals->nstdhdl, wi);
            }
            if (ir->bExpanded)
            {
                /* nstexpanded should be a multiple of nstcalcenergy */
                check_nst("nstcalcenergy", ir->nstcalcenergy, "nstexpanded", &ir->expandedvals->nstexpanded, wi);
            }
            /* for storing exact averages nstenergy should be
             * a multiple of nstcalcenergy
             */
            check_nst("nstcalcenergy", ir->nstcalcenergy, "nstenergy", &ir->nstenergy, wi);
        }

        // Inquire all MDModules, if their parameters match with the energy
        // calculation frequency
        gmx::EnergyCalculationFrequencyErrors energyCalculationFrequencyErrors(ir->nstcalcenergy);
        mdModulesNotifiers.preProcessingNotifier_.notify(&energyCalculationFrequencyErrors);

        // Emit all errors from the energy calculation frequency checks
        for (const std::string& energyFrequencyErrorMessage :
             energyCalculationFrequencyErrors.errorMessages())
        {
            wi->addError(energyFrequencyErrorMessage);
        }
    }

    if (ir->nsteps == 0 && !ir->bContinuation)
    {
        wi->addNote(
                "For a correct single-point energy evaluation with nsteps = 0, use "
                "continuation = yes to avoid constraining the input coordinates.");
    }

    /* LD STUFF */
    if ((EI_SD(ir->eI) || ir->eI == IntegrationAlgorithm::BD) && ir->bContinuation && ir->ld_seed != -1)
    {
        wi->addNote(
                "You are doing a continuation with SD or BD, make sure that ld_seed is "
                "different from the previous run (using ld_seed=-1 will ensure this)");
    }

    /* TPI STUFF */
    if (EI_TPI(ir->eI))
    {
        sprintf(err_buf, "TPI only works with pbc = %s", c_pbcTypeNames[PbcType::Xyz].c_str());
        CHECK(ir->pbcType != PbcType::Xyz);
        sprintf(err_buf, "with TPI nstlist should be larger than zero");
        CHECK(ir->nstlist <= 0);
        sprintf(err_buf, "TPI does not work with full electrostatics other than PME");
        CHECK(usingFullElectrostatics(ir->coulombtype) && !usingPme(ir->coulombtype));
    }

    /* SHAKE / LINCS */
    if ((opts->nshake > 0) && (opts->bMorse))
    {
        sprintf(warn_buf, "Using morse bond-potentials while constraining bonds is useless");
        wi->addWarning(warn_buf);
    }

    if ((EI_SD(ir->eI) || ir->eI == IntegrationAlgorithm::BD) && ir->bContinuation && ir->ld_seed != -1)
    {
        wi->addNote(
                "You are doing a continuation with SD or BD, make sure that ld_seed is "
                "different from the previous run (using ld_seed=-1 will ensure this)");
    }
    /* verify simulated tempering options */

    if (ir->bSimTemp)
    {
        /* We print a warning if users input elements of temperature-lambdas that are > 1 */
        bool bElementsGreaterThanOne = false;
        bool bAllTempZero            = true;
        for (i = 0; i < fep->n_lambda; i++)
        {
            /* We only forbid elements of temperature-lambdas that are < 0 */
            if (fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Temperature)][i] < 0)
            {
                auto message = gmx::formatString(
                        "Entry %d for %s must be greater than 0, instead is %g",
                        i,
                        enumValueToString(FreeEnergyPerturbationCouplingType::Temperature),
                        fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Temperature)][i]);
                wi->addError(message);
            }
            if (fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Temperature)][i] > 1)
            {
                bElementsGreaterThanOne = true;
            }
            if (fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Temperature)][i] > 0)
            {
                bAllTempZero = false;
            }
        }
        if (bElementsGreaterThanOne)
        {
            auto warningText = gmx::formatString(
                    "One or more entries for %s are greater than 1. Please only use this if you "
                    "are aware of "
                    "what you are doing. Check for inconsistencies with "
                    "simulated-tempering-scaling.",
                    enumValueToString(FreeEnergyPerturbationCouplingType::Temperature));
            wi->addWarning(warningText);
        }

        sprintf(err_buf, "if simulated tempering is on, temperature-lambdas may not be all zero");
        CHECK(bAllTempZero == TRUE);

        sprintf(err_buf, "Simulated tempering is currently only compatible with md-vv");
        CHECK(ir->eI != IntegrationAlgorithm::VV);

        /* check compatability of the temperature coupling with simulated tempering */

        if (ir->etc == TemperatureCoupling::NoseHoover)
        {
            sprintf(warn_buf,
                    "Nose-Hoover based temperature control such as [%s] my not be "
                    "entirelyconsistent with simulated tempering",
                    enumValueToString(ir->etc));
            wi->addNote(warn_buf);
        }

        /* check that the temperatures make sense */

        sprintf(err_buf,
                "Higher simulated tempering temperature (%g) must be >= than the simulated "
                "tempering lower temperature (%g)",
                ir->simtempvals->simtemp_high,
                ir->simtempvals->simtemp_low);
        CHECK(ir->simtempvals->simtemp_high <= ir->simtempvals->simtemp_low);

        sprintf(err_buf,
                "Higher simulated tempering temperature (%g) must be >= zero",
                ir->simtempvals->simtemp_high);
        CHECK(ir->simtempvals->simtemp_high <= 0);

        sprintf(err_buf,
                "Lower simulated tempering temperature (%g) must be >= zero",
                ir->simtempvals->simtemp_low);
        CHECK(ir->simtempvals->simtemp_low <= 0);
    }

    /* verify free energy options */

    if (ir->efep != FreeEnergyPerturbationType::No)
    {
        fep = ir->fepvals.get();
        sprintf(err_buf, "The soft-core power is %d and can only be 1 or 2", fep->sc_power);
        CHECK(fep->sc_alpha != 0 && fep->sc_power != 1 && fep->sc_power != 2);

        sprintf(err_buf,
                "The soft-core sc-r-power is %d and can only be 6. (sc-r-power 48 is no longer "
                "supported.)",
                static_cast<int>(fep->sc_r_power));
        CHECK(fep->sc_alpha != 0 && fep->sc_r_power != 6.0);

        /* We need special warnings if init-lambda > 1 && delta-lambda < 0 */
        if (fep->delta_lambda < 0 && fep->init_lambda_without_states > 1)
        {
            if (fep->n_lambda > 1)
            {
                /* warn about capping if lambda vector is provided as user input */
                double stepNumberWhenLambdaIsOne =
                        (1.0 - fep->init_lambda_without_states) / fep->delta_lambda;
                int64_t intStepNumberWhenLambdaIsOne =
                        static_cast<int64_t>(std::round(stepNumberWhenLambdaIsOne));

                auto warningText = gmx::formatString(
                        "With init-lambda = %g and delta_lambda = %g, the lambda "
                        "components won't change before step %" PRId64 " of %" PRId64
                        " simulation steps in total. "
                        "Consider setting init-lambda to a value less or equal to 1.\n",
                        fep->init_lambda_without_states,
                        fep->delta_lambda,
                        intStepNumberWhenLambdaIsOne,
                        ir->nsteps);
                wi->addWarning(warningText);
            }
            else if (fep->n_lambda <= 0)
            {
                /* issue an error if soft-core potentials are used */
                if (fep->sc_alpha > 0 || fep->softcoreFunction == SoftcoreType::Gapsys)
                {
                    auto message = gmx::formatString(
                            "You set init-lambda greater than 1 and provided no lambda vector as "
                            "input. "
                            "Therefore, coul-lambdas and vdw-lambdas will (initially) be greater "
                            "than 1. "
                            "This is not compatible with using soft-core potentials.\n");
                    wi->addError(message);
                }
            }
        }

        if (fep->delta_lambda != 0)
        {
            /* warn about capping */
            int64_t intStepNumberWhenLambdaIsCapped = ir->nsteps;

            if (fep->init_fep_state >= 0 && fep->init_fep_state < fep->n_lambda)
            {
                double deltaLambdaWithMultiplier = ((fep->n_lambda - 1) * fep->delta_lambda);
                if (fep->delta_lambda > 0)
                {
                    double stepNumberWhenLambdaIsCapped =
                            (fep->n_lambda - 1 - fep->init_fep_state) / deltaLambdaWithMultiplier;
                    intStepNumberWhenLambdaIsCapped =
                            static_cast<int64_t>(std::round(stepNumberWhenLambdaIsCapped));
                }
                else if (fep->delta_lambda < 0)
                {
                    double stepNumberWhenLambdaIsCapped =
                            (0 - fep->init_fep_state) / deltaLambdaWithMultiplier;
                    intStepNumberWhenLambdaIsCapped =
                            static_cast<int64_t>(std::round(stepNumberWhenLambdaIsCapped));
                }

                if (intStepNumberWhenLambdaIsCapped < ir->nsteps || ir->nsteps < 0)
                {
                    auto warningText = gmx::formatString(
                            "With init-lambda-state = %d and delta_lambda = %g, the lambda "
                            "components "
                            "won't change anymore after step %" PRId64
                            " until the end of the simulation after %" PRId64 " steps.\n",
                            fep->init_fep_state,
                            fep->delta_lambda,
                            intStepNumberWhenLambdaIsCapped,
                            ir->nsteps);
                    wi->addWarning(warningText);
                }
            }

            else if (fep->init_lambda_without_states >= 0)
            {
                if (fep->delta_lambda > 0)
                {
                    double stepNumberWhenLambdaIsCapped =
                            (1.0 - fep->init_lambda_without_states) / fep->delta_lambda;
                    stepNumberWhenLambdaIsCapped = std::max(stepNumberWhenLambdaIsCapped, 0.0);
                    intStepNumberWhenLambdaIsCapped =
                            static_cast<int64_t>(std::round(stepNumberWhenLambdaIsCapped));

                    /* There's no upper limit (capping) if no lambda value array is specified by the
                     * user. However, soft-core potentials may not be used with coul-lambdas or
                     * vdw-lambdas greater than 1. Make sure to error out.
                     */
                    if ((intStepNumberWhenLambdaIsCapped < ir->nsteps || ir->nsteps < 0)
                        && fep->n_lambda <= 0)
                    {
                        if (fep->sc_alpha > 0 || fep->softcoreFunction == SoftcoreType::Gapsys)
                        {
                            auto message = gmx::formatString(
                                    "With init-lambda = %g and delta_lambda = %g and no explicit "
                                    "input, "
                                    "coul-lambdas and vdw-lambdas will be greater than 1 after "
                                    "step %" PRId64 " of in total %" PRId64
                                    " steps. "
                                    "This is not compatible with using soft-core potentials.\n",
                                    fep->init_lambda_without_states,
                                    fep->delta_lambda,
                                    intStepNumberWhenLambdaIsCapped,
                                    ir->nsteps);
                            wi->addError(message);
                        }
                        /* No capping warning needed. */
                        intStepNumberWhenLambdaIsCapped = ir->nsteps;
                    }
                }
                else if (fep->delta_lambda < 0)
                {
                    double stepNumberWhenLambdaIsCapped =
                            (0.0 - fep->init_lambda_without_states) / fep->delta_lambda;
                    intStepNumberWhenLambdaIsCapped =
                            static_cast<int64_t>(std::round(stepNumberWhenLambdaIsCapped));
                }
                if (intStepNumberWhenLambdaIsCapped < ir->nsteps
                    || (ir->nsteps < 0 && !(fep->delta_lambda > 0 && fep->n_lambda <= 0)))
                {
                    auto warningText = gmx::formatString(
                            "With init-lambda = %g and delta_lambda = %g, the lambda components "
                            "won't change anymore after step %" PRId64
                            " until the end of the simulation after %" PRId64 " steps.\n",
                            fep->init_lambda_without_states,
                            fep->delta_lambda,
                            intStepNumberWhenLambdaIsCapped,
                            ir->nsteps);
                    wi->addWarning(warningText);
                }
            }

            else if (fep->n_lambda == 1)
            {
                auto warningText = gmx::formatString(
                        "You have set delta-lambda non-zero "
                        "while using a lambda vector that has one column. "
                        "The lambda components will therefore stay the same, "
                        "and delta-lambda has no effect.");
                wi->addWarning(warningText);
            }
        }

        sprintf(err_buf,
                "Can't use positive delta-lambda (%g) with expanded ensemble simulations",
                fep->delta_lambda);
        CHECK(fep->delta_lambda > 0 && (ir->efep == FreeEnergyPerturbationType::Expanded));

        sprintf(err_buf, "Can only use expanded ensemble with md-vv (for now)");
        CHECK(!(EI_VV(ir->eI)) && (ir->efep == FreeEnergyPerturbationType::Expanded));

        sprintf(err_buf, "Free-energy not implemented for Ewald");
        CHECK(ir->coulombtype == CoulombInteractionType::Ewald);

        /* check validty of lambda inputs */
        if (fep->n_lambda == 0)
        {
            /* Clear output in case of no states:*/
            sprintf(err_buf, "init-lambda-state set to %d: no lambda states are defined.", fep->init_fep_state);
            CHECK((fep->init_fep_state >= 0) && (fep->n_lambda == 0));
        }
        else
        {
            sprintf(err_buf,
                    "initial thermodynamic state %d does not exist, only goes to %d",
                    fep->init_fep_state,
                    fep->n_lambda - 1);
            CHECK((fep->init_fep_state >= fep->n_lambda));
        }

        sprintf(err_buf,
                "Lambda state must be set, either with init-lambda-state or with init-lambda");
        CHECK((fep->init_fep_state < 0) && (fep->init_lambda_without_states < 0));

        sprintf(err_buf,
                "init-lambda=%g while init-lambda-state=%d. Lambda state must be set either with "
                "init-lambda-state or with init-lambda, but not both",
                fep->init_lambda_without_states,
                fep->init_fep_state);
        CHECK((fep->init_fep_state >= 0) && (fep->init_lambda_without_states >= 0));


        if ((fep->init_lambda_without_states >= 0) && (fep->delta_lambda == 0))
        {
            int n_lambda_terms;
            n_lambda_terms = 0;
            for (i = 0; i < static_cast<int>(FreeEnergyPerturbationCouplingType::Count); i++)
            {
                if (fep->separate_dvdl[i])
                {
                    n_lambda_terms++;
                }
            }
            if (n_lambda_terms > 1)
            {
                sprintf(warn_buf,
                        "If lambda vector states (fep-lambdas, coul-lambdas etc.) are set, don't "
                        "use init-lambda to set lambda state (except for slow growth). Use "
                        "init-lambda-state instead.");
                wi->addWarning(warn_buf);
            }

            if (n_lambda_terms < 2 && fep->n_lambda > 0)
            {
                wi->addNote(
                        "init-lambda is deprecated for setting lambda state (except for slow "
                        "growth). Use init-lambda-state instead.");
            }
        }

        /*  Free Energy Checks -- In an ideal world, slow growth and FEP would
            be treated differently, but that's the next step */

        for (j = 0; j < static_cast<int>(FreeEnergyPerturbationCouplingType::Count); j++)
        {
            auto enumValue = static_cast<FreeEnergyPerturbationCouplingType>(j);
            /* We must restrict elements of the lambda value array for coulomb and vdw to [0,1]
             * if soft-core potentials are used.
             */
            if (((enumValue == FreeEnergyPerturbationCouplingType::Coul)
                 || (enumValue == FreeEnergyPerturbationCouplingType::Vdw))
                && ((fep->sc_alpha > 0) || (fep->softcoreFunction == SoftcoreType::Gapsys)))
            {
                for (i = 0; i < fep->n_lambda; i++)
                {
                    if ((fep->all_lambda[j][i] < 0) || (fep->all_lambda[j][i] > 1))
                    {
                        auto message = gmx::formatString(
                                "As you are using soft-core potentials, entry %d for %s must "
                                "be between 0 and 1, instead is %g",
                                i,
                                enumValueToString(enumValue),
                                fep->all_lambda[j][i]);
                        wi->addError(message);
                    }
                }
            }
            else
            {
                for (i = 0; i < fep->n_lambda; i++)
                {
                    /* We only forbid elements of the lambda value array that are < 0 */
                    if (fep->all_lambda[j][i] < 0)
                    {
                        auto message = gmx::formatString(
                                "Entry %d for %s must be greater than or equal to 0, instead is %g",
                                i,
                                enumValueToString(enumValue),
                                fep->all_lambda[j][i]);
                        wi->addError(message);
                    }
                }
            }
        }

        // the following warning needs to be triggered for the cases where vdw and coul of atoms are changing, but
        //    a) Beutler softcore is used with the softened LJ, yet not softened Coulomb
        bool softcoreConditionCheckBeutler = (fep->softcoreFunction == SoftcoreType::Beutler)
                                             && (fep->sc_alpha > 0) && (!fep->bScCoul);
        //    b) Gapsys softcore is used with the softened LJ, yet not softened Coulomb
        bool softcoreConditionCheckGapsys = (fep->softcoreFunction == SoftcoreType::Gapsys)
                                            && (fep->scGapsysScaleLinpointLJ > 0)
                                            && (fep->scGapsysScaleLinpointQ == 0);
        if (softcoreConditionCheckBeutler || softcoreConditionCheckGapsys)
        {

            for (i = 0; i < fep->n_lambda; i++)
            {
                sprintf(err_buf,
                        "For state %d, vdw-lambdas (%f) is changing with vdw softcore, while "
                        "coul-lambdas (%f) is nonzero without coulomb softcore: this will lead to "
                        "crashes, and is not supported.",
                        i,
                        fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)][i],
                        fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)][i]);
                CHECK((fep->sc_alpha > 0)
                      && (((fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)][i] > 0.0)
                           && (fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)][i] < 1.0))
                          && ((fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)][i] > 0.0)
                              && (fep->all_lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)][i]
                                  < 1.0))));
            }
        }

        if ((fep->softcoreFunction == SoftcoreType::Beutler) && (fep->bScCoul)
            && (usingPme(ir->coulombtype)))
        {
            // PME is formulated for computing 1/r potential, whereas the Coulombic potential softened by means of Beutler softcore no longer has this functional form.
            // This Warning is issued when PME is used for electrostatics and Beutler softcore is applied to the Coulombic interactions.
            // For the Gapsys softcore there is no such issue, because for the interactions more distant than the Coulombic linearization point, the potential takes the standard form of 1/r.
            // It is ensured that this linearization point does not exceed the short range electrostatic cutoff.

            real sigma, lambda, r_sc;

            sigma = 0.34;
            /* Maximum estimate for A and B charges equal with lambda power 1 */
            lambda = 0.5;
            r_sc = std::pow(lambda * fep->sc_alpha * std::pow(sigma / ir->rcoulomb, fep->sc_r_power) + 1.0,
                            1.0 / fep->sc_r_power);
            sprintf(warn_buf,
                    "With PME there is a minor soft core effect present at the cut-off, "
                    "proportional to (LJsigma/rcoulomb)^%g. This could have a minor effect on "
                    "energy conservation, but usually other effects dominate. With a common sigma "
                    "value of %g nm the fraction of the particle-particle potential at the cut-off "
                    "at lambda=%g is around %.1e, while ewald-rtol is %.1e.",
                    fep->sc_r_power,
                    sigma,
                    lambda,
                    r_sc - 1.0,
                    ir->ewald_rtol);
            wi->addNote(warn_buf);
        }

        if (fep->softcoreFunction == SoftcoreType::Gapsys)
        {
            if (fep->scGapsysScaleLinpointQ < 0.0)
            {
                sprintf(warn_buf,
                        "sc_scale_linpoint_Q_gapsys is equal %g but must be >= 0",
                        fep->scGapsysScaleLinpointQ);
                wi->addNote(warn_buf);
            }

            if ((fep->scGapsysScaleLinpointLJ < 0.0) || (fep->scGapsysScaleLinpointLJ >= 1.0))
            {
                sprintf(warn_buf,
                        "sc_scale_linpoint_LJ_gapsys is equal %g but must be in [0,1) when used "
                        "with "
                        "sc_function=gapsys.",
                        fep->scGapsysScaleLinpointLJ);
                wi->addNote(warn_buf);
            }
        }
    }

    if ((ir->bSimTemp) || (ir->efep == FreeEnergyPerturbationType::Expanded))
    {
        fep = ir->fepvals.get();

        /* checking equilibration of weights inputs for validity */

        sprintf(err_buf,
                "weight-equil-number-all-lambda (%d) is ignored if lmc-weights-equil is not equal "
                "to %s",
                expand->equil_n_at_lam,
                enumValueToString(LambdaWeightWillReachEquilibrium::NumAtLambda));
        CHECK((expand->equil_n_at_lam > 0)
              && (expand->elmceq != LambdaWeightWillReachEquilibrium::NumAtLambda));

        sprintf(err_buf,
                "weight-equil-number-samples (%d) is ignored if lmc-weights-equil is not equal to "
                "%s",
                expand->equil_samples,
                enumValueToString(LambdaWeightWillReachEquilibrium::Samples));
        CHECK((expand->equil_samples > 0) && (expand->elmceq != LambdaWeightWillReachEquilibrium::Samples));

        sprintf(err_buf,
                "weight-equil-number-steps (%d) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_steps,
                enumValueToString(LambdaWeightWillReachEquilibrium::Steps));
        CHECK((expand->equil_steps > 0) && (expand->elmceq != LambdaWeightWillReachEquilibrium::Steps));

        sprintf(err_buf,
                "weight-equil-wl-delta (%d) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_samples,
                enumValueToString(LambdaWeightWillReachEquilibrium::WLDelta));
        CHECK((expand->equil_wl_delta > 0) && (expand->elmceq != LambdaWeightWillReachEquilibrium::WLDelta));

        sprintf(err_buf,
                "weight-equil-count-ratio (%f) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_ratio,
                enumValueToString(LambdaWeightWillReachEquilibrium::Ratio));
        CHECK((expand->equil_ratio > 0) && (expand->elmceq != LambdaWeightWillReachEquilibrium::Ratio));

        sprintf(err_buf,
                "weight-equil-number-all-lambda (%d) must be a positive integer if "
                "lmc-weights-equil=%s",
                expand->equil_n_at_lam,
                enumValueToString(LambdaWeightWillReachEquilibrium::NumAtLambda));
        CHECK((expand->equil_n_at_lam <= 0)
              && (expand->elmceq == LambdaWeightWillReachEquilibrium::NumAtLambda));

        sprintf(err_buf,
                "weight-equil-number-samples (%d) must be a positive integer if "
                "lmc-weights-equil=%s",
                expand->equil_samples,
                enumValueToString(LambdaWeightWillReachEquilibrium::Samples));
        CHECK((expand->equil_samples <= 0) && (expand->elmceq == LambdaWeightWillReachEquilibrium::Samples));

        sprintf(err_buf,
                "weight-equil-number-steps (%d) must be a positive integer if lmc-weights-equil=%s",
                expand->equil_steps,
                enumValueToString(LambdaWeightWillReachEquilibrium::Steps));
        CHECK((expand->equil_steps <= 0) && (expand->elmceq == LambdaWeightWillReachEquilibrium::Steps));

        sprintf(err_buf,
                "weight-equil-wl-delta (%f) must be > 0 if lmc-weights-equil=%s",
                expand->equil_wl_delta,
                enumValueToString(LambdaWeightWillReachEquilibrium::WLDelta));
        CHECK((expand->equil_wl_delta <= 0)
              && (expand->elmceq == LambdaWeightWillReachEquilibrium::WLDelta));

        sprintf(err_buf,
                "weight-equil-count-ratio (%f) must be > 0 if lmc-weights-equil=%s",
                expand->equil_ratio,
                enumValueToString(LambdaWeightWillReachEquilibrium::Ratio));
        CHECK((expand->equil_ratio <= 0) && (expand->elmceq == LambdaWeightWillReachEquilibrium::Ratio));

        sprintf(err_buf,
                "lmc-weights-equil=%s only possible when lmc-stats = %s or lmc-stats %s",
                enumValueToString(LambdaWeightWillReachEquilibrium::WLDelta),
                enumValueToString(LambdaWeightCalculation::WL),
                enumValueToString(LambdaWeightCalculation::WWL));
        CHECK((expand->elmceq == LambdaWeightWillReachEquilibrium::WLDelta) && (!EWL(expand->elamstats)));

        sprintf(err_buf, "lmc-repeats (%d) must be greater than 0", expand->lmc_repeats);
        CHECK((expand->lmc_repeats <= 0));
        sprintf(err_buf, "minimum-var-min (%d) must be greater than 0", expand->minvarmin);
        CHECK((expand->minvarmin <= 0));
        sprintf(err_buf, "weight-c-range (%d) must be greater or equal to 0", expand->c_range);
        CHECK((expand->c_range < 0));
        sprintf(err_buf,
                "init-lambda-state (%d) must be zero if lmc-forced-nstart (%d)> 0 and lmc-move != "
                "'no'",
                fep->init_fep_state,
                expand->lmc_forced_nstart);
        CHECK((fep->init_fep_state != 0) && (expand->lmc_forced_nstart > 0)
              && (expand->elmcmove != LambdaMoveCalculation::No));
        sprintf(err_buf, "lmc-forced-nstart (%d) must not be negative", expand->lmc_forced_nstart);
        CHECK((expand->lmc_forced_nstart < 0));
        sprintf(err_buf,
                "init-lambda-state (%d) must be in the interval [0,number of lambdas)",
                fep->init_fep_state);
        CHECK((fep->init_fep_state < 0) || (fep->init_fep_state >= fep->n_lambda));

        sprintf(err_buf, "init-wl-delta (%f) must be greater than or equal to 0", expand->init_wl_delta);
        CHECK((expand->init_wl_delta < 0));
        sprintf(err_buf, "wl-ratio (%f) must be between 0 and 1", expand->wl_ratio);
        CHECK((expand->wl_ratio <= 0) || (expand->wl_ratio >= 1));
        sprintf(err_buf, "wl-scale (%f) must be between 0 and 1", expand->wl_scale);
        CHECK((expand->wl_scale <= 0) || (expand->wl_scale >= 1));

        /* if there is no temperature control, we need to specify an MC temperature */
        if (!integratorHasReferenceTemperature(*ir)
            && (expand->elmcmove != LambdaMoveCalculation::No) && (expand->mc_temp <= 0.0))
        {
            sprintf(err_buf,
                    "If the system has no reference temperature, and lmc-mcmove!='no', mc_temp "
                    "must be set to a positive number");
            wi->addError(err_buf);
        }
        if (expand->nstTij > 0)
        {
            sprintf(err_buf, "nstlog must be non-zero");
            CHECK(ir->nstlog == 0);
            // Avoid modulus by zero in the case that already triggered an error exit.
            if (ir->nstlog != 0)
            {
                sprintf(err_buf,
                        "nst-transition-matrix (%d) must be an integer multiple of nstlog (%d)",
                        expand->nstTij,
                        ir->nstlog);
                CHECK((expand->nstTij % ir->nstlog) != 0);
            }
        }
    }

    /* PBC/WALLS */
    sprintf(err_buf, "walls only work with pbc=%s", c_pbcTypeNames[PbcType::XY].c_str());
    CHECK(ir->nwall && ir->pbcType != PbcType::XY);

    /* VACUUM STUFF */
    if (ir->pbcType != PbcType::Xyz && ir->nwall != 2)
    {
        if (ir->pbcType == PbcType::No)
        {
            if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
            {
                wi->addWarning("Turning off pressure coupling for vacuum system");
                ir->pressureCouplingOptions.epc = PressureCoupling::No;
            }
        }
        else
        {
            sprintf(err_buf,
                    "Can not have pressure coupling with pbc=%s",
                    c_pbcTypeNames[ir->pbcType].c_str());
            CHECK(ir->pressureCouplingOptions.epc != PressureCoupling::No);
        }
        sprintf(err_buf, "Can not have Ewald with pbc=%s", c_pbcTypeNames[ir->pbcType].c_str());
        CHECK(usingFullElectrostatics(ir->coulombtype));

        sprintf(err_buf,
                "Can not have dispersion correction with pbc=%s",
                c_pbcTypeNames[ir->pbcType].c_str());
        CHECK(ir->eDispCorr != DispersionCorrectionType::No);
    }

    if (ir->rlist == 0.0)
    {
        sprintf(err_buf,
                "can only have neighborlist cut-off zero (=infinite)\n"
                "with coulombtype = %s or coulombtype = %s\n"
                "without periodic boundary conditions (pbc = %s) and\n"
                "rcoulomb and rvdw set to zero",
                enumValueToString(CoulombInteractionType::Cut),
                enumValueToString(CoulombInteractionType::User),
                c_pbcTypeNames[PbcType::No].c_str());
        CHECK(((ir->coulombtype != CoulombInteractionType::Cut)
               && (ir->coulombtype != CoulombInteractionType::User))
              || (ir->pbcType != PbcType::No) || (ir->rcoulomb != 0.0) || (ir->rvdw != 0.0));

        if (ir->nstlist > 0)
        {
            wi->addNote(
                    "Simulating without cut-offs can be (slightly) faster with nstlist=0, "
                    "nstype=simple and only one MPI rank");
        }
    }

    /* COMM STUFF */
    if (ir->nstcomm == 0)
    {
        // TODO Change this behaviour. There should be exactly one way
        // to turn off an algorithm.
        ir->comm_mode = ComRemovalAlgorithm::No;
    }
    if (ir->comm_mode != ComRemovalAlgorithm::No)
    {
        if (ir->nstcomm < 0)
        {
            // TODO Such input was once valid. Now that we've been
            // helpful for a few years, we should reject such input,
            // lest we have to support every historical decision
            // forever.
            wi->addWarning(
                    "If you want to remove the rotation around the center of mass, you should set "
                    "comm_mode = Angular instead of setting nstcomm < 0. nstcomm is modified to "
                    "its absolute value");
            ir->nstcomm = std::abs(ir->nstcomm);
        }

        if (ir->nstcalcenergy > 0 && ir->nstcomm < ir->nstcalcenergy
            && ir->comm_mode != ComRemovalAlgorithm::LinearAccelerationCorrection)
        {
            wi->addNote(
                    "nstcomm < nstcalcenergy defeats the purpose of nstcalcenergy, consider "
                    "setting nstcomm equal to nstcalcenergy for less overhead");
        }

        if (ir->comm_mode == ComRemovalAlgorithm::Angular)
        {
            sprintf(err_buf,
                    "Can not remove the rotation around the center of mass with periodic "
                    "molecules");
            CHECK(ir->bPeriodicMols);
            if (ir->pbcType != PbcType::No)
            {
                wi->addWarning(
                        "Removing the rotation around the center of mass in a periodic system, "
                        "this can lead to artifacts. Only use this on a single (cluster of) "
                        "molecules. This cluster should not cross periodic boundaries.");
            }
        }
    }

    if (EI_STATE_VELOCITY(ir->eI) && !EI_SD(ir->eI) && ir->pbcType == PbcType::No
        && ir->comm_mode != ComRemovalAlgorithm::Angular)
    {
        sprintf(warn_buf,
                "Tumbling and flying ice-cubes: We are not removing rotation around center of mass "
                "in a non-periodic system. You should probably set comm_mode = ANGULAR or use "
                "integrator = %s.",
                enumValueToString(IntegrationAlgorithm::SD1));
        wi->addNote(warn_buf);
    }

    /* TEMPERATURE COUPLING */
    if (ir->etc == TemperatureCoupling::Yes)
    {
        ir->etc = TemperatureCoupling::Berendsen;
        wi->addNote(
                "Old option for temperature coupling given: "
                "changing \"yes\" to \"Berendsen\"\n");
    }

    if ((ir->etc == TemperatureCoupling::NoseHoover)
        || (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk))
    {
        if (ir->opts.nhchainlength < 1)
        {
            sprintf(warn_buf,
                    "number of Nose-Hoover chains (currently %d) cannot be less than 1,reset to "
                    "1\n",
                    ir->opts.nhchainlength);
            ir->opts.nhchainlength = 1;
            wi->addWarning(warn_buf);
        }

        if (ir->etc == TemperatureCoupling::NoseHoover && !EI_VV(ir->eI) && ir->opts.nhchainlength > 1)
        {
            wi->addNote(

                    "leapfrog does not yet support Nose-Hoover chains, nhchainlength reset to 1");
            ir->opts.nhchainlength = 1;
        }
    }
    else
    {
        ir->opts.nhchainlength = 0;
    }

    if (ir->eI == IntegrationAlgorithm::VVAK)
    {
        sprintf(err_buf,
                "%s implemented primarily for validation, and requires nsttcouple = 1 and "
                "nstpcouple = 1.",
                enumValueToString(IntegrationAlgorithm::VVAK));
        CHECK((ir->nsttcouple != 1) || (ir->pressureCouplingOptions.nstpcouple != 1));
    }

    if (ETC_ANDERSEN(ir->etc))
    {
        sprintf(err_buf,
                "%s temperature control not supported for integrator %s.",
                enumValueToString(ir->etc),
                enumValueToString(ir->eI));
        CHECK(!(EI_VV(ir->eI)));

        if (ir->nstcomm > 0 && (ir->etc == TemperatureCoupling::Andersen))
        {
            sprintf(warn_buf,
                    "Center of mass removal not necessary for %s.  All velocities of coupled "
                    "groups are rerandomized periodically, so flying ice cube errors will not "
                    "occur.",
                    enumValueToString(ir->etc));
            wi->addNote(warn_buf);
        }

        sprintf(err_buf,
                "nstcomm must be 1, not %d for %s, as velocities of atoms in coupled groups are "
                "randomized every time step",
                ir->nstcomm,
                enumValueToString(ir->etc));
        CHECK(ir->nstcomm > 1 && (ir->etc == TemperatureCoupling::Andersen));
        if (opts->nshake != 0 && ir->etc == TemperatureCoupling::Andersen)
        {
            auto message = gmx::formatString(
                    "%s temperature control does not work with constraints. "
                    "Use %s instead",
                    enumValueToString(TemperatureCoupling::Andersen),
                    enumValueToString(TemperatureCoupling::AndersenMassive));
            wi->addError(message);
        }
    }

    if (ir->etc == TemperatureCoupling::Berendsen)
    {
        sprintf(warn_buf,
                "The %s thermostat does not generate the correct kinetic energy distribution, "
                "and should not be used for new production simulations (in our opinion). "
                "We would recommend the %s thermostat.",
                enumValueToString(ir->etc),
                enumValueToString(TemperatureCoupling::VRescale));
        wi->addWarning(warn_buf);
    }

    /* PRESSURE COUPLING */
    if (ir->pressureCouplingOptions.epc == PressureCoupling::Berendsen)
    {
        sprintf(warn_buf,
                "The %s barostat does not generate any strictly correct ensemble, "
                "and should not be used for new production simulations (in our opinion). "
                "We recommend using the %s barostat instead.",
                enumValueToString(ir->pressureCouplingOptions.epc),
                enumValueToString(PressureCoupling::CRescale));
        wi->addWarning(warn_buf);
    }

    if (ir->pressureCouplingOptions.epc == PressureCoupling::Isotropic)
    {
        ir->pressureCouplingOptions.epc = PressureCoupling::Berendsen;
        wi->addNote(
                "Old option for pressure coupling given: "
                "changing \"Isotropic\" to \"Berendsen\"\n");
    }

    if (ir->pressureCouplingOptions.epc == PressureCoupling::CRescale)
    {
        switch (ir->pressureCouplingOptions.epct)
        {
            case PressureCouplingType::Isotropic:
            case PressureCouplingType::SemiIsotropic:
            case PressureCouplingType::SurfaceTension: break; // supported
            default:
                sprintf(err_buf,
                        "C-rescale does not support pressure coupling type %s yet\n",
                        enumValueToString(ir->pressureCouplingOptions.epct));
                wi->addError(err_buf);
        }
    }

    if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
    {
        dt_pcoupl = ir->pressureCouplingOptions.nstpcouple * ir->delta_t;

        sprintf(err_buf, "tau-p must be > 0 instead of %g\n", ir->pressureCouplingOptions.tau_p);
        CHECK(ir->pressureCouplingOptions.tau_p <= 0);

        if (ir->pressureCouplingOptions.tau_p / dt_pcoupl
            < pcouple_min_integration_steps(ir->pressureCouplingOptions.epc) - 10 * GMX_REAL_EPS)
        {
            sprintf(warn_buf,
                    "For proper integration of the %s barostat, tau-p (%g) should be at least %d "
                    "times larger than nstpcouple*dt (%g)",
                    enumValueToString(ir->pressureCouplingOptions.epc),
                    ir->pressureCouplingOptions.tau_p,
                    pcouple_min_integration_steps(ir->pressureCouplingOptions.epc),
                    dt_pcoupl);
            wi->addWarning(warn_buf);
        }

        sprintf(err_buf,
                "compressibility must be > 0 when using pressure"
                " coupling %s\n",
                enumValueToString(ir->pressureCouplingOptions.epc));
        CHECK(ir->pressureCouplingOptions.compress[XX][XX] < 0
              || ir->pressureCouplingOptions.compress[YY][YY] < 0
              || ir->pressureCouplingOptions.compress[ZZ][ZZ] < 0
              || (trace(ir->pressureCouplingOptions.compress) == 0
                  && ir->pressureCouplingOptions.compress[YY][XX] <= 0
                  && ir->pressureCouplingOptions.compress[ZZ][XX] <= 0
                  && ir->pressureCouplingOptions.compress[ZZ][YY] <= 0));

        if (PressureCoupling::ParrinelloRahman == ir->pressureCouplingOptions.epc && opts->bGenVel)
        {
            sprintf(warn_buf,
                    "You are generating velocities so I am assuming you "
                    "are equilibrating a system. You are using "
                    "%s pressure coupling, but this can be "
                    "unstable for equilibration. If your system crashes, try "
                    "equilibrating first with Berendsen pressure coupling. If "
                    "you are not equilibrating the system, you can probably "
                    "ignore this warning.",
                    enumValueToString(ir->pressureCouplingOptions.epc));
            wi->addWarning(warn_buf);
        }
    }

    if (!EI_VV(ir->eI))
    {
        if (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            wi->addError("MTTK pressure coupling requires a Velocity-verlet integrator");
        }
    }

    /* ELECTROSTATICS */
    /* More checks are in triple check (grompp.c) */

    if (ir->coulombtype == CoulombInteractionType::Switch)
    {
        sprintf(warn_buf,
                "coulombtype = %s is only for testing purposes and can lead to serious "
                "artifacts, advice: use coulombtype = %s",
                enumValueToString(ir->coulombtype),
                enumValueToString(CoulombInteractionType::RFZero));
        wi->addWarning(warn_buf);
    }

    if (usingRF(ir->coulombtype) && ir->epsilon_rf == 1 && ir->epsilon_r != 1)
    {
        sprintf(warn_buf,
                "epsilon-r = %g and epsilon-rf = 1 with reaction field, proceeding assuming old "
                "format and exchanging epsilon-r and epsilon-rf",
                ir->epsilon_r);
        wi->addWarning(warn_buf);
        ir->epsilon_rf = ir->epsilon_r;
        ir->epsilon_r  = 1.0;
    }

    if (ir->epsilon_r == 0)
    {
        sprintf(err_buf,
                "It is pointless to use long-range electrostatics with infinite relative "
                "permittivity."
                "Since you are effectively turning of electrostatics, a plain cutoff will be much "
                "faster.");
        CHECK(usingFullElectrostatics(ir->coulombtype));
    }

    if (getenv("GMX_DO_GALACTIC_DYNAMICS") == nullptr)
    {
        sprintf(err_buf, "epsilon-r must be >= 0 instead of %g\n", ir->epsilon_r);
        CHECK(ir->epsilon_r < 0);
    }

    if (usingRF(ir->coulombtype))
    {
        /* reaction field (at the cut-off) */

        if (ir->coulombtype == CoulombInteractionType::RFZero && ir->epsilon_rf != 0)
        {
            sprintf(warn_buf,
                    "With coulombtype = %s, epsilon-rf must be 0, assuming you meant epsilon_rf=0",
                    enumValueToString(ir->coulombtype));
            wi->addWarning(warn_buf);
            ir->epsilon_rf = 0.0;
        }

        sprintf(err_buf, "epsilon-rf must be >= epsilon-r");
        CHECK((ir->epsilon_rf < ir->epsilon_r && ir->epsilon_rf != 0) || (ir->epsilon_r == 0));
        if (ir->epsilon_rf == ir->epsilon_r)
        {
            sprintf(warn_buf,
                    "Using epsilon-rf = epsilon-r with %s does not make sense",
                    enumValueToString(ir->coulombtype));
            wi->addWarning(warn_buf);
        }
    }
    /* Allow rlist>rcoulomb for tabulated long range stuff. This just
     * means the interaction is zero outside rcoulomb, but it helps to
     * provide accurate energy conservation.
     */
    if (ir_coulomb_might_be_zero_at_cutoff(ir))
    {
        if (ir_coulomb_switched(ir))
        {
            sprintf(err_buf,
                    "With coulombtype = %s rcoulomb_switch must be < rcoulomb. Or, better: Use the "
                    "potential modifier options!",
                    enumValueToString(ir->coulombtype));
            CHECK(ir->rcoulomb_switch >= ir->rcoulomb);
        }
    }

    if (ir->coulombtype == CoulombInteractionType::Switch || ir->coulombtype == CoulombInteractionType::Shift)
    {
        sprintf(err_buf,
                "Explicit switch/shift coulomb interactions cannot be used in combination with a "
                "secondary coulomb-modifier.");
        CHECK(ir->coulomb_modifier != InteractionModifiers::None);
    }
    if (ir->vdwtype == VanDerWaalsType::Switch || ir->vdwtype == VanDerWaalsType::Shift)
    {
        sprintf(err_buf,
                "Explicit switch/shift vdw interactions cannot be used in combination with a "
                "secondary vdw-modifier.");
        CHECK(ir->vdw_modifier != InteractionModifiers::None);
    }

    if (ir->coulombtype == CoulombInteractionType::Switch || ir->coulombtype == CoulombInteractionType::Shift
        || ir->vdwtype == VanDerWaalsType::Switch || ir->vdwtype == VanDerWaalsType::Shift)
    {
        sprintf(warn_buf,
                "The switch/shift interaction settings are just for compatibility; you will get "
                "better "
                "performance from applying potential modifiers to your interactions!\n");
        wi->addNote(warn_buf);
    }

    if (ir->coulombtype == CoulombInteractionType::PmeSwitch
        || ir->coulomb_modifier == InteractionModifiers::PotSwitch)
    {
        if (ir->rcoulomb_switch / ir->rcoulomb < 0.9499)
        {
            real percentage = 100 * (ir->rcoulomb - ir->rcoulomb_switch) / ir->rcoulomb;
            sprintf(warn_buf,
                    "The switching range should be 5%% or less (currently %.2f%% using a switching "
                    "range of %4f-%4f) for accurate electrostatic energies, energy conservation "
                    "will be good regardless, since ewald_rtol = %g.",
                    percentage,
                    ir->rcoulomb_switch,
                    ir->rcoulomb,
                    ir->ewald_rtol);
            wi->addWarning(warn_buf);
        }
    }

    if (ir->vdwtype == VanDerWaalsType::Switch || ir->vdw_modifier == InteractionModifiers::PotSwitch)
    {
        if (ir->rvdw_switch == 0)
        {
            sprintf(warn_buf,
                    "rvdw-switch is equal 0 even though you are using a switched Lennard-Jones "
                    "potential.  This suggests it was not set in the mdp, which can lead to large "
                    "energy errors.  In GROMACS, 0.05 to 0.1 nm is often a reasonable vdw "
                    "switching range.");
            wi->addWarning(warn_buf);
        }
    }

    if (usingFullElectrostatics(ir->coulombtype))
    {
        if (ir->coulombtype == CoulombInteractionType::PmeSwitch
            || ir->coulombtype == CoulombInteractionType::PmeUser
            || ir->coulombtype == CoulombInteractionType::PmeUserSwitch)
        {
            sprintf(err_buf,
                    "With coulombtype = %s, rcoulomb must be <= rlist",
                    enumValueToString(ir->coulombtype));
            CHECK(ir->rcoulomb > ir->rlist);
        }
    }

    if (usingPme(ir->coulombtype) || usingLJPme(ir->vdwtype))
    {
        // TODO: Move these checks into the ewald module with the options class
        int orderMin = 3;
        int orderMax = (ir->coulombtype == CoulombInteractionType::P3mAD ? 8 : 12);

        if (ir->pme_order < orderMin || ir->pme_order > orderMax)
        {
            sprintf(warn_buf,
                    "With coulombtype = %s, you should have %d <= pme-order <= %d",
                    enumValueToString(ir->coulombtype),
                    orderMin,
                    orderMax);
            wi->addError(warn_buf);
        }
    }

    if (ir->nwall == 2 && usingFullElectrostatics(ir->coulombtype))
    {
        if (ir->ewald_geometry == EwaldGeometry::ThreeD)
        {
            sprintf(warn_buf,
                    "With pbc=%s you should use ewald-geometry=%s",
                    c_pbcTypeNames[ir->pbcType].c_str(),
                    enumValueToString(EwaldGeometry::ThreeDC));
            wi->addWarning(warn_buf);
        }
        /* This check avoids extra pbc coding for exclusion corrections */
        sprintf(err_buf, "wall-ewald-zfac should be >= 2");
        CHECK(ir->wall_ewald_zfac < 2);
    }
    if ((ir->ewald_geometry == EwaldGeometry::ThreeDC) && (ir->pbcType != PbcType::XY)
        && usingFullElectrostatics(ir->coulombtype))
    {
        sprintf(warn_buf,
                "With %s and ewald_geometry = %s you should use pbc = %s",
                enumValueToString(ir->coulombtype),
                enumValueToString(EwaldGeometry::ThreeDC),
                c_pbcTypeNames[PbcType::XY].c_str());
        wi->addWarning(warn_buf);
    }
    if ((ir->epsilon_surface != 0) && usingFullElectrostatics(ir->coulombtype))
    {
        sprintf(err_buf, "Cannot have periodic molecules with epsilon_surface > 0");
        CHECK(ir->bPeriodicMols);
        sprintf(warn_buf, "With epsilon_surface > 0 all molecules should be neutral.");
        wi->addNote(warn_buf);
        sprintf(warn_buf,
                "With epsilon_surface > 0 you can only use domain decomposition "
                "when there are only small molecules with all bonds constrained (mdrun will check "
                "for this).");
        wi->addNote(warn_buf);
    }

    if (ir_vdw_switched(ir))
    {
        sprintf(err_buf, "With switched vdw forces or potentials, rvdw-switch must be < rvdw");
        CHECK(ir->rvdw_switch >= ir->rvdw);

        if (ir->rvdw_switch < 0.5 * ir->rvdw)
        {
            sprintf(warn_buf,
                    "You are applying a switch function to vdw forces or potentials from %g to %g "
                    "nm, which is more than half the interaction range, whereas switch functions "
                    "are intended to act only close to the cut-off.",
                    ir->rvdw_switch,
                    ir->rvdw);
            wi->addNote(warn_buf);
        }
    }

    if (ir->vdwtype == VanDerWaalsType::Pme)
    {
        if (!(ir->vdw_modifier == InteractionModifiers::None
              || ir->vdw_modifier == InteractionModifiers::PotShift))
        {
            sprintf(err_buf,
                    "With vdwtype = %s, the only supported modifiers are %s and %s",
                    enumValueToString(ir->vdwtype),
                    enumValueToString(InteractionModifiers::PotShift),
                    enumValueToString(InteractionModifiers::None));
            wi->addError(err_buf);
        }
    }

    if (ir->vdwtype == VanDerWaalsType::User && ir->eDispCorr != DispersionCorrectionType::No)
    {
        wi->addNote(
                "You have selected user tables with dispersion correction, the dispersion "
                "will be corrected to -C6/r^6 beyond rvdw_switch (the tabulated interaction "
                "between rvdw_switch and rvdw will not be double counted). Make sure that you "
                "really want dispersion correction to -C6/r^6.");
    }

    if (ir->eI == IntegrationAlgorithm::LBFGS
        && (ir->coulombtype == CoulombInteractionType::Cut || ir->vdwtype == VanDerWaalsType::Cut)
        && ir->rvdw != 0)
    {
        wi->addWarning("For efficient BFGS minimization, use switch/shift/pme instead of cut-off.");
    }

    if (ir->eI == IntegrationAlgorithm::LBFGS && ir->nbfgscorr <= 0)
    {
        wi->addWarning("Using L-BFGS with nbfgscorr<=0 just gets you steepest descent.");
    }

    /* IMPLICIT SOLVENT */
    if (ir->coulombtype == CoulombInteractionType::GBNotused)
    {
        sprintf(warn_buf, "Invalid option %s for coulombtype", enumValueToString(ir->coulombtype));
        wi->addError(warn_buf);
    }

    if (ir->bQMMM)
    {
        wi->addError("The QMMM integration you are trying to use is no longer supported");
    }

    if (ir->bAdress)
    {
        gmx_fatal(FARGS, "AdResS simulations are no longer supported");
    }

    // cosine acceleration is only supported in leap-frog
    if (ir->cos_accel != 0.0 && ir->eI != IntegrationAlgorithm::MD)
    {
        wi->addError("cos-acceleration is only supported by integrator = md");
    }

    // do checks on the initialization of the flow profile with deform
    if (ir_haveBoxDeformation(*ir) && !opts->deformInitFlow)
    {
        if (opts->bGenVel)
        {
            wi->addError(
                    "When the box is deformed and velocities are generated, the flow profile "
                    "should be initialized by setting deform-init-flow=yes");
        }
        else if (!ir->bContinuation)
        {
            wi->addNote(
                    "Unless the velocities in the initial configuration already obey the flow "
                    "profile, the flow profile should be initialized by setting "
                    "deform-init-flow=yes when using the deform option");
        }
    }
}

/* interpret a number of doubles from a string and put them in an array,
   after allocating space for them.
   str = the input string
   n = the (pre-allocated) number of doubles read
   r = the output array of doubles. */
static std::vector<real> parse_n_real(const std::string& str, int* n, WarningHandler* wi)
{
    auto values = gmx::splitString(str);
    *n          = values.size();

    std::vector<real> r;
    for (int i = 0; i < *n; i++)
    {
        try
        {
            r.emplace_back(gmx::fromString<real>(values[i]));
        }
        catch (gmx::GromacsException&)
        {
            wi->addError("Invalid value " + values[i]
                         + " in string in mdp file. Expected a real number.");
        }
    }
    return r;
}


static void do_fep_params(t_inputrec*                ir,
                          gmx::ArrayRef<std::string> fep_lambda,
                          char                       weights[STRLEN],
                          WarningHandler*            wi)
{

    int         i, j, max_n_lambda, nweights;
    t_lambda*   fep    = ir->fepvals.get();
    t_expanded* expand = ir->expandedvals.get();
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<real>> count_fep_lambdas;
    bool                                                                         bOneLambda = TRUE;
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, int>               nfep;

    /* FEP input processing */
    /* first, identify the number of lambda values for each type.
       All that are nonzero must have the same number */

    for (auto i : keysOf(nfep))
    {
        count_fep_lambdas[i] = parse_n_real(fep_lambda[static_cast<int>(i)], &(nfep[i]), wi);
    }

    /* now, determine the number of components.  All must be either zero, or equal. */

    max_n_lambda = 0;
    for (auto i : keysOf(nfep))
    {
        if (nfep[i] > max_n_lambda)
        {
            max_n_lambda = nfep[i]; /* here's a nonzero one.  All of them
                                       must have the same number if its not zero.*/
            break;
        }
    }

    for (auto i : keysOf(nfep))
    {
        if (nfep[i] == 0)
        {
            ir->fepvals->separate_dvdl[i] = FALSE;
        }
        else if (nfep[i] == max_n_lambda)
        {
            if (i != FreeEnergyPerturbationCouplingType::Temperature) /* we treat this differently -- not really a reason to compute
                                         the derivative with respect to the temperature currently */
            {
                ir->fepvals->separate_dvdl[i] = TRUE;
            }
        }
        else
        {
            gmx_fatal(FARGS,
                      "Number of lambdas (%d) for FEP type %s not equal to number of other types "
                      "(%d)",
                      nfep[i],
                      enumValueToString(i),
                      max_n_lambda);
        }
    }
    /* we don't print out dhdl if the temperature is changing, since we can't correctly define dhdl in this case */
    ir->fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Temperature] = FALSE;

    /* the number of lambdas is the number we've read in, which is either zero
       or the same for all */
    fep->n_lambda = max_n_lambda;

    /* if init_lambda is defined, we need to set lambda */
    if ((fep->init_lambda_without_states > 0) && (fep->n_lambda == 0))
    {
        ir->fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Fep] = TRUE;
    }
    /* otherwise allocate the space for all of the lambdas, and transfer the data */
    for (auto i : keysOf(nfep))
    {
        fep->all_lambda[i].resize(fep->n_lambda);
        if (nfep[i] > 0) /* if it's zero, then the count_fep_lambda arrays
                            are zero */
        {
            for (j = 0; j < fep->n_lambda; j++)
            {
                fep->all_lambda[i][j] = static_cast<double>(count_fep_lambdas[i][j]);
            }
        }
    }

    /* "fep-vals" is either zero or the full number. If zero, we'll need to define fep-lambdas for
       internal bookkeeping -- for now, init_lambda */

    if (nfep[FreeEnergyPerturbationCouplingType::Fep] == 0 && fep->init_lambda_without_states >= 0)
    {
        for (i = 0; i < fep->n_lambda; i++)
        {
            fep->all_lambda[FreeEnergyPerturbationCouplingType::Fep][i] = fep->init_lambda_without_states;
        }
    }

    /* check to see if only a single component lambda is defined, and soft core is defined.
       In this case, turn on coulomb soft core */

    if (max_n_lambda == 0)
    {
        bOneLambda = TRUE;
    }
    else
    {
        for (auto i : keysOf(nfep))
        {
            if ((nfep[i] != 0) && (i != FreeEnergyPerturbationCouplingType::Fep))
            {
                bOneLambda = FALSE;
            }
        }
    }
    if ((bOneLambda) && (fep->sc_alpha > 0))
    {
        fep->bScCoul = TRUE;
    }

    /* Fill in the others with the efptFEP if they are not explicitly
       specified (i.e. nfep[i] == 0).  This means if fep is not defined,
       they are all zero. */

    for (auto i : keysOf(nfep))
    {
        if ((nfep[i] == 0) && (i != FreeEnergyPerturbationCouplingType::Fep))
        {
            for (j = 0; j < fep->n_lambda; j++)
            {
                fep->all_lambda[i][j] = fep->all_lambda[FreeEnergyPerturbationCouplingType::Fep][j];
            }
        }
    }


    /* now read in the weights */
    expand->init_lambda_weights = parse_n_real(weights, &nweights, wi);
    if (nweights == 0)
    {
        expand->init_lambda_weights.resize(fep->n_lambda); /* initialize to zero */
    }
    else if (nweights != fep->n_lambda)
    {
        gmx_fatal(FARGS,
                  "Number of weights (%d) is not equal to number of lambda values (%d)",
                  nweights,
                  fep->n_lambda);
    }
    if ((expand->nstexpanded < 0) && (ir->efep != FreeEnergyPerturbationType::No))
    {
        expand->nstexpanded = fep->nstdhdl;
        /* if you don't specify nstexpanded when doing expanded ensemble free energy calcs, it is set to nstdhdl */
    }
}


static void do_simtemp_params(t_inputrec* ir)
{
    ir->simtempvals->temperatures.resize(ir->fepvals->n_lambda);
    getSimTemps(ir->fepvals->n_lambda,
                ir->simtempvals.get(),
                ir->fepvals->all_lambda[FreeEnergyPerturbationCouplingType::Temperature]);
}

template<typename T>
void convertInts(WarningHandler* wi, gmx::ArrayRef<const std::string> inputs, const char* name, T* outputs)
{
    int i = 0;
    for (const auto& input : inputs)
    {
        try
        {
            outputs[i] = gmx::fromStdString<T>(input);
        }
        catch (gmx::GromacsException&)
        {
            auto message = gmx::formatString(
                    "Invalid value for mdp option %s. %s should only consist of integers separated "
                    "by spaces.",
                    name,
                    name);
            wi->addError(message);
        }
        ++i;
    }
}

static void convertReals(WarningHandler* wi, gmx::ArrayRef<const std::string> inputs, const char* name, real* outputs)
{
    int i = 0;
    for (const auto& input : inputs)
    {
        try
        {
            outputs[i] = gmx::fromString<real>(input);
        }
        catch (gmx::GromacsException&)
        {
            auto message = gmx::formatString(
                    "Invalid value for mdp option %s. %s should only consist of real numbers "
                    "separated by spaces.",
                    name,
                    name);
            wi->addError(message);
        }
        ++i;
    }
}

static void convertRvecs(WarningHandler* wi, gmx::ArrayRef<const std::string> inputs, const char* name, rvec* outputs)
{
    int i = 0, d = 0;
    for (const auto& input : inputs)
    {
        try
        {
            outputs[i][d] = gmx::fromString<real>(input);
        }
        catch (gmx::GromacsException&)
        {
            auto message = gmx::formatString(
                    "Invalid value for mdp option %s. %s should only consist of real numbers "
                    "separated by spaces.",
                    name,
                    name);
            wi->addError(message);
        }
        ++d;
        if (d == DIM)
        {
            d = 0;
            ++i;
        }
    }
}

static void do_wall_params(t_inputrec* ir, char* wall_atomtype, char* wall_density, t_gromppopts* opts, WarningHandler* wi)
{
    opts->wall_atomtype[0] = nullptr;
    opts->wall_atomtype[1] = nullptr;

    ir->wall_atomtype[0] = -1;
    ir->wall_atomtype[1] = -1;
    ir->wall_density[0]  = 0;
    ir->wall_density[1]  = 0;

    if (ir->nwall > 0)
    {
        auto wallAtomTypes = gmx::splitString(wall_atomtype);
        if (wallAtomTypes.size() != size_t(ir->nwall))
        {
            gmx_fatal(FARGS,
                      "Expected %d elements for wall_atomtype, found %zu",
                      ir->nwall,
                      wallAtomTypes.size());
        }
        GMX_RELEASE_ASSERT(ir->nwall < 3, "Invalid number of walls");
        for (int i = 0; i < ir->nwall; i++)
        {
            opts->wall_atomtype[i] = gmx_strdup(wallAtomTypes[i].c_str());
        }

        if (ir->wall_type == WallType::NineThree || ir->wall_type == WallType::TenFour)
        {
            auto wallDensity = gmx::splitString(wall_density);
            if (wallDensity.size() != size_t(ir->nwall))
            {
                gmx_fatal(FARGS,
                          "Expected %d elements for wall-density, found %zu",
                          ir->nwall,
                          wallDensity.size());
            }
            convertReals(wi, wallDensity, "wall-density", ir->wall_density);
            for (int i = 0; i < ir->nwall; i++)
            {
                if (ir->wall_density[i] <= 0)
                {
                    gmx_fatal(FARGS, "wall-density[%d] = %f\n", i, ir->wall_density[i]);
                }
            }
        }
    }
}

static void add_wall_energrps(SimulationGroups* groups, int nwall, t_symtab* symtab)
{
    if (nwall > 0)
    {
        AtomGroupIndices* grps = &(groups->groups[SimulationAtomGroupType::EnergyOutput]);
        for (int i = 0; i < nwall; i++)
        {
            groups->groupNames.emplace_back(put_symtab(symtab, gmx::formatString("wall%d", i).c_str()));
            grps->emplace_back(groups->groupNames.size() - 1);
        }
    }
}

static void read_expandedparams(std::vector<t_inpfile>* inp, t_expanded* expand, WarningHandler* wi)
{
    /* read expanded ensemble parameters */
    printStringNewline(inp, "expanded ensemble variables");
    expand->nstexpanded = get_eint(inp, "nstexpanded", -1, wi);
    expand->elamstats   = getEnum<LambdaWeightCalculation>(inp, "lmc-stats", wi);
    expand->elmcmove    = getEnum<LambdaMoveCalculation>(inp, "lmc-move", wi);
    expand->elmceq      = getEnum<LambdaWeightWillReachEquilibrium>(inp, "lmc-weights-equil", wi);
    expand->equil_n_at_lam = get_eint(inp, "weight-equil-number-all-lambda", -1, wi);
    expand->equil_samples  = get_eint(inp, "weight-equil-number-samples", -1, wi);
    expand->equil_steps    = get_eint(inp, "weight-equil-number-steps", -1, wi);
    expand->equil_wl_delta = get_ereal(inp, "weight-equil-wl-delta", -1, wi);
    expand->equil_ratio    = get_ereal(inp, "weight-equil-count-ratio", -1, wi);
    printStringNewline(inp, "Seed for Monte Carlo in lambda space");
    expand->lmc_seed          = get_eint(inp, "lmc-seed", -1, wi);
    expand->mc_temp           = get_ereal(inp, "mc-temperature", -1, wi);
    expand->lmc_repeats       = get_eint(inp, "lmc-repeats", 1, wi);
    expand->gibbsdeltalam     = get_eint(inp, "lmc-gibbsdelta", -1, wi);
    expand->lmc_forced_nstart = get_eint(inp, "lmc-forced-nstart", 0, wi);
    expand->bSymmetrizedTMatrix =
            (getEnum<Boolean>(inp, "symmetrized-transition-matrix", wi) != Boolean::No);
    expand->nstTij        = get_eint(inp, "nst-transition-matrix", -1, wi);
    expand->minvarmin     = get_eint(inp, "mininum-var-min", 100, wi); /*default is reasonable */
    expand->c_range       = get_eint(inp, "weight-c-range", 0, wi);    /* default is just C=0 */
    expand->wl_scale      = get_ereal(inp, "wl-scale", 0.8, wi);
    expand->wl_ratio      = get_ereal(inp, "wl-ratio", 0.8, wi);
    expand->init_wl_delta = get_ereal(inp, "init-wl-delta", 1.0, wi);
    expand->bWLoneovert   = (getEnum<Boolean>(inp, "wl-oneovert", wi) != Boolean::No);
}

/*! \brief Return whether an end state with the given coupling-lambda
 * value describes fully-interacting VDW.
 *
 * \param[in]  couple_lambda_value  Enumeration ecouplam value describing the end state
 * \return                          Whether VDW is on (i.e. the user chose vdw or vdw-q in the .mdp file)
 */
static bool couple_lambda_has_vdw_on(int couple_lambda_value)
{
    return (couple_lambda_value == ecouplamVDW || couple_lambda_value == ecouplamVDWQ);
}

namespace
{

class MdpErrorHandler : public gmx::IKeyValueTreeErrorHandler
{
public:
    explicit MdpErrorHandler(WarningHandler* wi) : wi_(wi), mapping_(nullptr) {}

    void setBackMapping(const gmx::IKeyValueTreeBackMapping& mapping) { mapping_ = &mapping; }

    bool onError(gmx::UserInputError* ex, const gmx::KeyValueTreePath& context) override
    {
        ex->prependContext(
                gmx::formatString("Error in mdp option \"%s\":", getOptionName(context).c_str()));
        std::string message = gmx::formatExceptionMessageToString(*ex);
        wi_->addError(message);
        return true;
    }

private:
    std::string getOptionName(const gmx::KeyValueTreePath& context)
    {
        if (mapping_ != nullptr)
        {
            gmx::KeyValueTreePath path = mapping_->originalPath(context);
            GMX_ASSERT(path.size() == 1, "Inconsistent mapping back to mdp options");
            return path[0];
        }
        GMX_ASSERT(context.size() == 1, "Inconsistent context for mdp option parsing");
        return context[0];
    }

    WarningHandler*                      wi_;
    const gmx::IKeyValueTreeBackMapping* mapping_;
};

} // namespace

void get_ir(const char*     mdparin,
            const char*     mdparout,
            gmx::MDModules* mdModules,
            t_inputrec*     ir,
            t_gromppopts*   opts,
            WriteMdpHeader  writeMdpHeader,
            WarningHandler* wi)
{
    char*       dumstr[2];
    double      dumdub[2][6];
    int         i, j, m;
    char        warn_buf[STRLEN];
    t_lambda*   fep    = ir->fepvals.get();
    t_expanded* expand = ir->expandedvals.get();

    const char* no_names[] = { "no", nullptr };

    init_inputrec_strings();
    gmx::TextInputFile     stream(mdparin);
    std::vector<t_inpfile> inp = read_inpfile(&stream, mdparin, wi);

    snew(dumstr[0], STRLEN);
    snew(dumstr[1], STRLEN);

    /* ignore the following deprecated commands */
    replace_inp_entry(inp, "title", nullptr);
    replace_inp_entry(inp, "cpp", nullptr);
    replace_inp_entry(inp, "domain-decomposition", nullptr);
    replace_inp_entry(inp, "andersen-seed", nullptr);
    replace_inp_entry(inp, "dihre", nullptr);
    replace_inp_entry(inp, "dihre-fc", nullptr);
    replace_inp_entry(inp, "dihre-tau", nullptr);
    replace_inp_entry(inp, "nstdihreout", nullptr);
    replace_inp_entry(inp, "nstcheckpoint", nullptr);
    replace_inp_entry(inp, "optimize-fft", nullptr);
    replace_inp_entry(inp, "adress_type", nullptr);
    replace_inp_entry(inp, "adress_const_wf", nullptr);
    replace_inp_entry(inp, "adress_ex_width", nullptr);
    replace_inp_entry(inp, "adress_hy_width", nullptr);
    replace_inp_entry(inp, "adress_ex_forcecap", nullptr);
    replace_inp_entry(inp, "adress_interface_correction", nullptr);
    replace_inp_entry(inp, "adress_site", nullptr);
    replace_inp_entry(inp, "adress_reference_coords", nullptr);
    replace_inp_entry(inp, "adress_tf_grp_names", nullptr);
    replace_inp_entry(inp, "adress_cg_grp_names", nullptr);
    replace_inp_entry(inp, "adress_do_hybridpairs", nullptr);
    replace_inp_entry(inp, "rlistlong", nullptr);
    replace_inp_entry(inp, "nstcalclr", nullptr);
    replace_inp_entry(inp, "pull-print-com2", nullptr);
    replace_inp_entry(inp, "gb-algorithm", nullptr);
    replace_inp_entry(inp, "nstgbradii", nullptr);
    replace_inp_entry(inp, "rgbradii", nullptr);
    replace_inp_entry(inp, "gb-epsilon-solvent", nullptr);
    replace_inp_entry(inp, "gb-saltconc", nullptr);
    replace_inp_entry(inp, "gb-obc-alpha", nullptr);
    replace_inp_entry(inp, "gb-obc-beta", nullptr);
    replace_inp_entry(inp, "gb-obc-gamma", nullptr);
    replace_inp_entry(inp, "gb-dielectric-offset", nullptr);
    replace_inp_entry(inp, "sa-algorithm", nullptr);
    replace_inp_entry(inp, "sa-surface-tension", nullptr);
    replace_inp_entry(inp, "ns-type", nullptr);

    /* replace the following commands with the clearer new versions*/
    replace_inp_entry(inp, "unconstrained-start", "continuation");
    replace_inp_entry(inp, "foreign-lambda", "fep-lambdas");
    replace_inp_entry(inp, "verlet-buffer-drift", "verlet-buffer-tolerance");
    replace_inp_entry(inp, "nstxtcout", "nstxout-compressed");
    replace_inp_entry(inp, "xtc-grps", "compressed-x-grps");
    replace_inp_entry(inp, "xtc-precision", "compressed-x-precision");
    replace_inp_entry(inp, "pull-print-com1", "pull-print-com");

    printStringNewline(&inp, "VARIOUS PREPROCESSING OPTIONS");
    printStringNoNewline(&inp, "Preprocessor information: use cpp syntax.");
    printStringNoNewline(&inp, "e.g.: -I/home/joe/doe -I/home/mary/roe");
    setStringEntry(&inp, "include", opts->include, nullptr);
    printStringNoNewline(
            &inp, "e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)");
    setStringEntry(&inp, "define", opts->define, nullptr);

    printStringNewline(&inp, "RUN CONTROL PARAMETERS");
    ir->eI = getEnum<IntegrationAlgorithm>(&inp, "integrator", wi);
    printStringNoNewline(&inp, "Start time and timestep in ps");
    ir->init_t  = get_ereal(&inp, "tinit", 0.0, wi);
    ir->delta_t = get_ereal(&inp, "dt", 0.001, wi);
    ir->nsteps  = get_eint64(&inp, "nsteps", 0, wi);
    printStringNoNewline(&inp, "For exact run continuation or redoing part of a run");
    ir->init_step = get_eint64(&inp, "init-step", 0, wi);
    printStringNoNewline(
            &inp, "Part index is updated automatically on checkpointing (keeps files separate)");
    ir->simulation_part = get_eint(&inp, "simulation-part", 1, wi);
    printStringNoNewline(&inp, "Multiple time-stepping");
    ir->useMts = (getEnum<Boolean>(&inp, "mts", wi) != Boolean::No);
    if (ir->useMts)
    {
        gmx::GromppMtsOpts& mtsOpts = opts->mtsOpts;
        mtsOpts.numLevels           = get_eint(&inp, "mts-levels", 2, wi);
        mtsOpts.level2Forces = setStringEntry(&inp, "mts-level2-forces", "longrange-nonbonded");
        mtsOpts.level2Factor = get_eint(&inp, "mts-level2-factor", 2, wi);

        // We clear after reading without dynamics to not force the user to remove MTS mdp options
        if (!EI_DYNAMICS(ir->eI))
        {
            ir->useMts = false;
        }
    }
    printStringNoNewline(&inp, "factor by which to increase the mass of the lightest atoms");
    ir->massRepartitionFactor = get_ereal(&inp, "mass-repartition-factor", 1.0, wi);
    printStringNoNewline(&inp, "mode for center of mass motion removal");
    ir->comm_mode = getEnum<ComRemovalAlgorithm>(&inp, "comm-mode", wi);
    printStringNoNewline(&inp, "number of steps for center of mass motion removal");
    ir->nstcomm = get_eint(&inp, "nstcomm", 100, wi);
    printStringNoNewline(&inp, "group(s) for center of mass motion removal");
    setStringEntry(&inp, "comm-grps", inputrecStrings->vcm, nullptr);

    printStringNewline(&inp, "LANGEVIN DYNAMICS OPTIONS");
    printStringNoNewline(&inp, "Friction coefficient (amu/ps) and random seed");
    ir->bd_fric = get_ereal(&inp, "bd-fric", 0.0, wi);
    ir->ld_seed = get_eint64(&inp, "ld-seed", -1, wi);

    /* Em stuff */
    printStringNewline(&inp, "ENERGY MINIMIZATION OPTIONS");
    printStringNoNewline(&inp, "Force tolerance and initial step-size");
    ir->em_tol      = get_ereal(&inp, "emtol", 10.0, wi);
    ir->em_stepsize = get_ereal(&inp, "emstep", 0.01, wi);
    printStringNoNewline(&inp, "Max number of iterations in relax-shells");
    ir->niter = get_eint(&inp, "niter", 20, wi);
    printStringNoNewline(&inp, "Step size (ps^2) for minimization of flexible constraints");
    ir->fc_stepsize = get_ereal(&inp, "fcstep", 0, wi);
    printStringNoNewline(&inp, "Frequency of steepest descents steps when doing CG");
    ir->nstcgsteep = get_eint(&inp, "nstcgsteep", 1000, wi);
    ir->nbfgscorr  = get_eint(&inp, "nbfgscorr", 10, wi);

    printStringNewline(&inp, "TEST PARTICLE INSERTION OPTIONS");
    ir->rtpi = get_ereal(&inp, "rtpi", 0.05, wi);

    /* Output options */
    printStringNewline(&inp, "OUTPUT CONTROL OPTIONS");
    printStringNoNewline(&inp, "Output frequency for coords (x), velocities (v) and forces (f)");
    ir->nstxout = get_eint(&inp, "nstxout", 0, wi);
    ir->nstvout = get_eint(&inp, "nstvout", 0, wi);
    ir->nstfout = get_eint(&inp, "nstfout", 0, wi);
    printStringNoNewline(&inp, "Output frequency for energies to log file and energy file");
    ir->nstlog        = get_eint(&inp, "nstlog", 1000, wi);
    ir->nstcalcenergy = get_eint(&inp, "nstcalcenergy", 100, wi);
    ir->nstenergy     = get_eint(&inp, "nstenergy", 1000, wi);
    printStringNoNewline(&inp, "Output frequency and precision for .xtc file");
    ir->nstxout_compressed      = get_eint(&inp, "nstxout-compressed", 0, wi);
    ir->x_compression_precision = get_ereal(&inp, "compressed-x-precision", 1000.0, wi);
    printStringNoNewline(&inp, "This selects the subset of atoms for the compressed");
    printStringNoNewline(&inp, "trajectory file. You can select multiple groups. By");
    printStringNoNewline(&inp, "default, all atoms will be written.");
    setStringEntry(&inp, "compressed-x-grps", inputrecStrings->x_compressed_groups, nullptr);
    printStringNoNewline(&inp, "Selection of energy groups");
    setStringEntry(&inp, "energygrps", inputrecStrings->energy, nullptr);

    /* Neighbor searching */
    printStringNewline(&inp, "NEIGHBORSEARCHING PARAMETERS");
    printStringNoNewline(&inp, "cut-off scheme (Verlet: particle based cut-offs)");
    ir->cutoff_scheme = getEnum<CutoffScheme>(&inp, "cutoff-scheme", wi);
    printStringNoNewline(&inp, "nblist update frequency");
    ir->nstlist = get_eint(&inp, "nstlist", 10, wi);
    printStringNoNewline(&inp, "Periodic boundary conditions: xyz, no, xy");
    // TODO This conversion should be removed when proper std:string handling will be added to get_eeenum(...), etc.
    std::vector<const char*> pbcTypesNamesChar;
    for (const auto& pbcTypeName : c_pbcTypeNames)
    {
        pbcTypesNamesChar.push_back(pbcTypeName.c_str());
    }
    ir->pbcType       = static_cast<PbcType>(get_eeenum(&inp, "pbc", pbcTypesNamesChar.data(), wi));
    ir->bPeriodicMols = getEnum<Boolean>(&inp, "periodic-molecules", wi) != Boolean::No;
    printStringNoNewline(&inp,
                         "Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,");
    printStringNoNewline(&inp, "a value of -1 means: use rlist");
    ir->verletbuf_tol = get_ereal(&inp, "verlet-buffer-tolerance", 0.005, wi);
    printStringNoNewline(
            &inp,
            "Allowed error in the average pressure due to the Verlet buffer for LJ interactions");
    ir->verletBufferPressureTolerance = get_ereal(&inp, "verlet-buffer-pressure-tolerance", 0.5, wi);
    printStringNoNewline(&inp, "nblist cut-off");
    ir->rlist = get_ereal(&inp, "rlist", 1.0, wi);
    printStringNoNewline(&inp, "long-range cut-off for switched potentials");

    /* Electrostatics */
    printStringNewline(&inp, "OPTIONS FOR ELECTROSTATICS AND VDW");
    printStringNoNewline(&inp, "Method for doing electrostatics");
    ir->coulombtype      = getEnum<CoulombInteractionType>(&inp, "coulombtype", wi);
    ir->coulomb_modifier = getEnum<InteractionModifiers>(&inp, "coulomb-modifier", wi);
    printStringNoNewline(&inp, "cut-off lengths");
    ir->rcoulomb_switch = get_ereal(&inp, "rcoulomb-switch", 0.0, wi);
    ir->rcoulomb        = get_ereal(&inp, "rcoulomb", 1.0, wi);
    printStringNoNewline(&inp, "Relative dielectric constant for the medium and the reaction field");
    ir->epsilon_r  = get_ereal(&inp, "epsilon-r", 1.0, wi);
    ir->epsilon_rf = get_ereal(&inp, "epsilon-rf", 0.0, wi);
    printStringNoNewline(&inp, "Method for doing Van der Waals");
    ir->vdwtype      = getEnum<VanDerWaalsType>(&inp, "vdw-type", wi);
    ir->vdw_modifier = getEnum<InteractionModifiers>(&inp, "vdw-modifier", wi);
    printStringNoNewline(&inp, "cut-off lengths");
    ir->rvdw_switch = get_ereal(&inp, "rvdw-switch", 0.0, wi);
    ir->rvdw        = get_ereal(&inp, "rvdw", 1.0, wi);
    printStringNoNewline(&inp, "Apply long range dispersion corrections for Energy and Pressure");
    ir->eDispCorr = getEnum<DispersionCorrectionType>(&inp, "DispCorr", wi);
    printStringNoNewline(&inp, "Extension of the potential lookup tables beyond the cut-off");
    ir->tabext = get_ereal(&inp, "table-extension", 1.0, wi);
    printStringNoNewline(&inp, "Separate tables between energy group pairs");
    setStringEntry(&inp, "energygrp-table", inputrecStrings->egptable, nullptr);
    printStringNoNewline(&inp, "Spacing for the PME/PPPM FFT grid");
    ir->fourier_spacing = get_ereal(&inp, "fourierspacing", 0.12, wi);
    printStringNoNewline(&inp, "FFT grid size, when a value is 0 fourierspacing will be used");
    ir->nkx = get_eint(&inp, "fourier-nx", 0, wi);
    ir->nky = get_eint(&inp, "fourier-ny", 0, wi);
    ir->nkz = get_eint(&inp, "fourier-nz", 0, wi);
    printStringNoNewline(&inp, "EWALD/PME/PPPM parameters");
    ir->pme_order              = get_eint(&inp, "pme-order", 4, wi);
    ir->ewald_rtol             = get_ereal(&inp, "ewald-rtol", 0.00001, wi);
    ir->ewald_rtol_lj          = get_ereal(&inp, "ewald-rtol-lj", 0.001, wi);
    ir->ljpme_combination_rule = getEnum<LongRangeVdW>(&inp, "lj-pme-comb-rule", wi);
    ir->ewald_geometry         = getEnum<EwaldGeometry>(&inp, "ewald-geometry", wi);
    ir->epsilon_surface        = get_ereal(&inp, "epsilon-surface", 0.0, wi);

    /* Implicit solvation is no longer supported, but we need grompp
       to be able to refuse old .mdp files that would have built a tpr
       to run it. Thus, only "no" is accepted. */
    ir->implicit_solvent = (get_eeenum(&inp, "implicit-solvent", no_names, wi) != 0);

    /* Coupling stuff */
    printStringNewline(&inp, "OPTIONS FOR WEAK COUPLING ALGORITHMS");
    ir->ensembleTemperatureSetting =
            getEnum<EnsembleTemperatureSetting>(&inp, "ensemble-temperature-setting", wi);
    ir->ensembleTemperature = get_ereal(&inp, "ensemble-temperature", -1, wi);
    printStringNoNewline(&inp, "Temperature coupling");
    ir->etc                = getEnum<TemperatureCoupling>(&inp, "tcoupl", wi);
    ir->nsttcouple         = get_eint(&inp, "nsttcouple", -1, wi);
    ir->opts.nhchainlength = get_eint(&inp, "nh-chain-length", 10, wi);
    ir->bPrintNHChains = (getEnum<Boolean>(&inp, "print-nose-hoover-chain-variables", wi) != Boolean::No);
    printStringNoNewline(&inp, "Groups to couple separately");
    setStringEntry(&inp, "tc-grps", inputrecStrings->tcgrps, nullptr);
    printStringNoNewline(&inp, "Time constant (ps) and reference temperature (K)");
    setStringEntry(&inp, "tau-t", inputrecStrings->tau_t, nullptr);
    setStringEntry(&inp, "ref-t", inputrecStrings->ref_t, nullptr);
    printStringNoNewline(&inp, "pressure coupling");
    ir->pressureCouplingOptions.epc        = getEnum<PressureCoupling>(&inp, "pcoupl", wi);
    ir->pressureCouplingOptions.epct       = getEnum<PressureCouplingType>(&inp, "pcoupltype", wi);
    ir->pressureCouplingOptions.nstpcouple = get_eint(&inp, "nstpcouple", -1, wi);
    printStringNoNewline(&inp, "Time constant (ps), compressibility (1/bar) and reference P (bar)");
    ir->pressureCouplingOptions.tau_p = get_ereal(&inp, "tau-p", 5.0, wi);
    setStringEntry(&inp, "compressibility", dumstr[0], nullptr);
    setStringEntry(&inp, "ref-p", dumstr[1], nullptr);
    printStringNoNewline(&inp, "Scaling of reference coordinates, No, All or COM");
    ir->pressureCouplingOptions.refcoord_scaling = getEnum<RefCoordScaling>(&inp, "refcoord-scaling", wi);

    /* QMMM */
    printStringNewline(&inp, "OPTIONS FOR QMMM calculations");
    ir->bQMMM = (getEnum<Boolean>(&inp, "QMMM", wi) != Boolean::No);
    printStringNoNewline(&inp, "Groups treated with MiMiC");
    setStringEntry(&inp, "QMMM-grps", inputrecStrings->QMMM, nullptr);

    /* Simulated annealing */
    printStringNewline(&inp, "SIMULATED ANNEALING");
    printStringNoNewline(&inp, "Type of annealing for each temperature group (no/single/periodic)");
    setStringEntry(&inp, "annealing", inputrecStrings->anneal, nullptr);
    printStringNoNewline(&inp,
                         "Number of time points to use for specifying annealing in each group");
    setStringEntry(&inp, "annealing-npoints", inputrecStrings->anneal_npoints, nullptr);
    printStringNoNewline(&inp, "List of times at the annealing points for each group");
    setStringEntry(&inp, "annealing-time", inputrecStrings->anneal_time, nullptr);
    printStringNoNewline(&inp, "Temp. at each annealing point, for each group.");
    setStringEntry(&inp, "annealing-temp", inputrecStrings->anneal_temp, nullptr);

    /* Startup run */
    printStringNewline(&inp, "GENERATE VELOCITIES FOR STARTUP RUN");
    opts->bGenVel = (getEnum<Boolean>(&inp, "gen-vel", wi) != Boolean::No);
    opts->tempi   = get_ereal(&inp, "gen-temp", 300.0, wi);
    opts->seed    = get_eint(&inp, "gen-seed", -1, wi);

    if (opts->seed == -1)
    {
        opts->seed      = static_cast<int>(gmx::makeRandomSeed());
        opts->bMadeSeed = true;
    }

    /* Shake stuff */
    printStringNewline(&inp, "OPTIONS FOR BONDS");
    opts->nshake = get_eeenum(&inp, "constraints", constraints, wi);
    printStringNoNewline(&inp, "Type of constraint algorithm");
    ir->eConstrAlg = getEnum<ConstraintAlgorithm>(&inp, "constraint-algorithm", wi);
    printStringNoNewline(&inp, "Do not constrain the start configuration");
    ir->bContinuation = (getEnum<Boolean>(&inp, "continuation", wi) != Boolean::No);
    printStringNoNewline(&inp,
                         "Use successive overrelaxation to reduce the number of shake iterations");
    ir->bShakeSOR = (getEnum<Boolean>(&inp, "Shake-SOR", wi) != Boolean::No);
    printStringNoNewline(&inp, "Relative tolerance of shake");
    ir->shake_tol = get_ereal(&inp, "shake-tol", 0.0001, wi);
    printStringNoNewline(&inp, "Highest order in the expansion of the constraint coupling matrix");
    ir->nProjOrder = get_eint(&inp, "lincs-order", 4, wi);
    printStringNoNewline(&inp, "Number of iterations in the final step of LINCS. 1 is fine for");
    printStringNoNewline(&inp, "normal simulations, but use 2 to conserve energy in NVE runs.");
    printStringNoNewline(&inp, "For energy minimization with constraints it should be 4 to 8.");
    ir->nLincsIter = get_eint(&inp, "lincs-iter", 1, wi);
    printStringNoNewline(&inp, "Lincs will write a warning to the stderr if in one step a bond");
    printStringNoNewline(&inp, "rotates over more degrees than");
    ir->LincsWarnAngle = get_ereal(&inp, "lincs-warnangle", 30.0, wi);
    printStringNoNewline(&inp, "Convert harmonic bonds to morse potentials");
    opts->bMorse = (getEnum<Boolean>(&inp, "morse", wi) != Boolean::No);

    /* Energy group exclusions */
    printStringNewline(&inp, "ENERGY GROUP EXCLUSIONS");
    printStringNoNewline(
            &inp, "Pairs of energy groups for which all non-bonded interactions are excluded");
    setStringEntry(&inp, "energygrp-excl", inputrecStrings->egpexcl, nullptr);

    /* Walls */
    printStringNewline(&inp, "WALLS");
    printStringNoNewline(
            &inp, "Number of walls, type, atom types, densities and box-z scale factor for Ewald");
    ir->nwall         = get_eint(&inp, "nwall", 0, wi);
    ir->wall_type     = getEnum<WallType>(&inp, "wall-type", wi);
    ir->wall_r_linpot = get_ereal(&inp, "wall-r-linpot", -1, wi);
    setStringEntry(&inp, "wall-atomtype", inputrecStrings->wall_atomtype, nullptr);
    setStringEntry(&inp, "wall-density", inputrecStrings->wall_density, nullptr);
    ir->wall_ewald_zfac = get_ereal(&inp, "wall-ewald-zfac", 3, wi);

    /* COM pulling */
    printStringNewline(&inp, "COM PULLING");
    ir->bPull = (getEnum<Boolean>(&inp, "pull", wi) != Boolean::No);
    if (ir->bPull)
    {
        ir->pull                        = std::make_unique<pull_params_t>();
        inputrecStrings->pullGroupNames = read_pullparams(&inp, ir->pull.get(), wi);

        if (ir->useMts)
        {
            for (int c = 0; c < ir->pull->ncoord; c++)
            {
                if (ir->pull->coord[c].eType == PullingAlgorithm::Constraint)
                {
                    wi->addError(
                            "Constraint COM pulling is not supported in combination with "
                            "multiple time stepping");
                    break;
                }
            }
        }
    }

    /* AWH biasing
       NOTE: needs COM pulling or free energy input */
    printStringNewline(&inp, "AWH biasing");
    ir->bDoAwh = (getEnum<Boolean>(&inp, "awh", wi) != Boolean::No);
    if (ir->bDoAwh)
    {
        ir->awhParams = std::make_unique<gmx::AwhParams>(&inp, wi);
    }

    /* Enforced rotation */
    printStringNewline(&inp, "ENFORCED ROTATION");
    printStringNoNewline(&inp, "Enforced rotation: No or Yes");
    ir->bRot = (getEnum<Boolean>(&inp, "rotation", wi) != Boolean::No);
    if (ir->bRot)
    {
        ir->rot                           = std::make_unique<t_rot>();
        inputrecStrings->rotateGroupNames = read_rotparams(&inp, ir->rot.get(), wi);
    }

    /* Interactive MD */
    ir->bIMD = FALSE;
    printStringNewline(&inp, "Group to display and/or manipulate in interactive MD session");
    setStringEntry(&inp, "IMD-group", inputrecStrings->imd_grp, nullptr);
    if (inputrecStrings->imd_grp[0] != '\0')
    {
        snew(ir->imd, 1);
        ir->bIMD = TRUE;
    }

    /* Refinement */
    printStringNewline(&inp, "NMR refinement stuff");
    printStringNoNewline(&inp, "Distance restraints type: No, Simple or Ensemble");
    ir->eDisre = getEnum<DistanceRestraintRefinement>(&inp, "disre", wi);
    printStringNoNewline(
            &inp, "Force weighting of pairs in one distance restraint: Conservative or Equal");
    ir->eDisreWeighting = getEnum<DistanceRestraintWeighting>(&inp, "disre-weighting", wi);
    printStringNoNewline(&inp, "Use sqrt of the time averaged times the instantaneous violation");
    ir->bDisreMixed = (getEnum<Boolean>(&inp, "disre-mixed", wi) != Boolean::No);
    ir->dr_fc       = get_ereal(&inp, "disre-fc", 1000.0, wi);
    ir->dr_tau      = get_ereal(&inp, "disre-tau", 0.0, wi);
    printStringNoNewline(&inp, "Output frequency for pair distances to energy file");
    ir->nstdisreout = get_eint(&inp, "nstdisreout", 100, wi);
    printStringNoNewline(&inp, "Orientation restraints: No or Yes");
    opts->bOrire = (getEnum<Boolean>(&inp, "orire", wi) != Boolean::No);
    printStringNoNewline(&inp, "Orientation restraints force constant and tau for time averaging");
    ir->orires_fc  = get_ereal(&inp, "orire-fc", 0.0, wi);
    ir->orires_tau = get_ereal(&inp, "orire-tau", 0.0, wi);
    setStringEntry(&inp, "orire-fitgrp", inputrecStrings->orirefitgrp, nullptr);
    printStringNoNewline(&inp, "Output frequency for trace(SD) and S to energy file");
    ir->nstorireout = get_eint(&inp, "nstorireout", 100, wi);

    /* free energy variables */
    printStringNewline(&inp, "Free energy variables");
    ir->efep = getEnum<FreeEnergyPerturbationType>(&inp, "free-energy", wi);
    setStringEntry(&inp, "couple-moltype", inputrecStrings->couple_moltype, nullptr);
    opts->couple_lam0  = get_eeenum(&inp, "couple-lambda0", couple_lam, wi);
    opts->couple_lam1  = get_eeenum(&inp, "couple-lambda1", couple_lam, wi);
    opts->bCoupleIntra = (getEnum<Boolean>(&inp, "couple-intramol", wi) != Boolean::No);

    fep->init_lambda_without_states = get_ereal(&inp, "init-lambda", -1, wi); /* start with -1 so
                                                                                 we can recognize if
                                                                                 it was not entered */
    fep->init_fep_state = get_eint(&inp, "init-lambda-state", -1, wi);
    fep->delta_lambda   = get_ereal(&inp, "delta-lambda", 0.0, wi);
    fep->nstdhdl        = get_eint(&inp, "nstdhdl", 50, wi);
    inputrecStrings->fep_lambda[FreeEnergyPerturbationCouplingType::Fep] =
            setStringEntry(&inp, "fep-lambdas", "");
    inputrecStrings->fep_lambda[FreeEnergyPerturbationCouplingType::Mass] =
            setStringEntry(&inp, "mass-lambdas", "");
    inputrecStrings->fep_lambda[FreeEnergyPerturbationCouplingType::Coul] =
            setStringEntry(&inp, "coul-lambdas", "");
    inputrecStrings->fep_lambda[FreeEnergyPerturbationCouplingType::Vdw] =
            setStringEntry(&inp, "vdw-lambdas", "");
    inputrecStrings->fep_lambda[FreeEnergyPerturbationCouplingType::Bonded] =
            setStringEntry(&inp, "bonded-lambdas", "");
    inputrecStrings->fep_lambda[FreeEnergyPerturbationCouplingType::Restraint] =
            setStringEntry(&inp, "restraint-lambdas", "");
    inputrecStrings->fep_lambda[FreeEnergyPerturbationCouplingType::Temperature] =
            setStringEntry(&inp, "temperature-lambdas", "");
    fep->lambda_neighbors = get_eint(&inp, "calc-lambda-neighbors", 1, wi);
    setStringEntry(&inp, "init-lambda-weights", inputrecStrings->lambda_weights, nullptr);
    fep->edHdLPrintEnergy        = getEnum<FreeEnergyPrintEnergy>(&inp, "dhdl-print-energy", wi);
    fep->softcoreFunction        = getEnum<SoftcoreType>(&inp, "sc-function", wi);
    fep->sc_alpha                = get_ereal(&inp, "sc-alpha", 0.0, wi);
    fep->sc_power                = get_eint(&inp, "sc-power", 1, wi);
    fep->sc_r_power              = get_ereal(&inp, "sc-r-power", 6.0, wi);
    fep->sc_sigma                = get_ereal(&inp, "sc-sigma", 0.3, wi);
    fep->bScCoul                 = (getEnum<Boolean>(&inp, "sc-coul", wi) != Boolean::No);
    fep->scGapsysScaleLinpointLJ = get_ereal(&inp, "sc-gapsys-scale-linpoint-lj", 0.85, wi);
    fep->scGapsysScaleLinpointQ  = get_ereal(&inp, "sc-gapsys-scale-linpoint-q", 0.3, wi);
    fep->scGapsysSigmaLJ         = get_ereal(&inp, "sc-gapsys-sigma-lj", 0.3, wi);
    fep->dh_hist_size            = get_eint(&inp, "dh_hist_size", 0, wi);
    fep->dh_hist_spacing         = get_ereal(&inp, "dh_hist_spacing", 0.1, wi);
    fep->separate_dhdl_file      = getEnum<SeparateDhdlFile>(&inp, "separate-dhdl-file", wi);
    fep->dhdl_derivatives        = getEnum<DhDlDerivativeCalculation>(&inp, "dhdl-derivatives", wi);
    fep->dh_hist_size            = get_eint(&inp, "dh_hist_size", 0, wi);
    fep->dh_hist_spacing         = get_ereal(&inp, "dh_hist_spacing", 0.1, wi);

    /* Non-equilibrium MD stuff */
    printStringNewline(&inp, "Non-equilibrium MD stuff");
    setStringEntry(&inp, "acc-grps", inputrecStrings->accelerationGroups, nullptr);
    setStringEntry(&inp, "accelerate", inputrecStrings->acceleration, nullptr);
    setStringEntry(&inp, "freezegrps", inputrecStrings->freeze, nullptr);
    setStringEntry(&inp, "freezedim", inputrecStrings->frdim, nullptr);
    ir->cos_accel = get_ereal(&inp, "cos-acceleration", 0, wi);
    setStringEntry(&inp, "deform", inputrecStrings->deform, nullptr);
    opts->deformInitFlow = (getEnum<Boolean>(&inp, "deform-init-flow", wi) != Boolean::No);

    /* simulated tempering variables */
    printStringNewline(&inp, "simulated tempering variables");
    ir->bSimTemp = (getEnum<Boolean>(&inp, "simulated-tempering", wi) != Boolean::No);
    ir->simtempvals->eSimTempScale = getEnum<SimulatedTempering>(&inp, "simulated-tempering-scaling", wi);
    ir->simtempvals->simtemp_low  = get_ereal(&inp, "sim-temp-low", 300.0, wi);
    ir->simtempvals->simtemp_high = get_ereal(&inp, "sim-temp-high", 300.0, wi);

    /* expanded ensemble variables */
    if (ir->efep == FreeEnergyPerturbationType::Expanded || ir->bSimTemp)
    {
        read_expandedparams(&inp, expand, wi);
    }

    /* Electric fields */
    {
        gmx::KeyValueTreeObject      convertedValues = flatKeyValueTreeFromInpFile(inp);
        gmx::KeyValueTreeTransformer transform;
        transform.rules()->addRule().keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
        mdModules->initMdpTransform(transform.rules());
        for (const auto& path : transform.mappedPaths())
        {
            GMX_ASSERT(path.size() == 1, "Inconsistent mapping back to mdp options");
            mark_einp_set(inp, path[0].c_str());
        }
        MdpErrorHandler errorHandler(wi);
        auto            result = transform.transform(convertedValues, &errorHandler);
        ir->params             = new gmx::KeyValueTreeObject(result.object());
        mdModules->adjustInputrecBasedOnModules(ir);
        errorHandler.setBackMapping(result.backMapping());
        mdModules->assignOptionsToModules(*ir->params, &errorHandler);
    }

    /* Ion/water position swapping ("computational electrophysiology") */
    printStringNewline(&inp,
                       "Ion/water position swapping for computational electrophysiology setups");
    printStringNoNewline(&inp, "Swap positions along direction: no, X, Y, Z");
    ir->eSwapCoords = getEnum<SwapType>(&inp, "swapcoords", wi);
    if (ir->eSwapCoords != SwapType::No)
    {
        char buf[STRLEN];
        int  nIonTypes;


        snew(ir->swap, 1);
        printStringNoNewline(&inp, "Swap attempt frequency");
        ir->swap->nstswap = get_eint(&inp, "swap-frequency", 1, wi);
        printStringNoNewline(&inp, "Number of ion types to be controlled");
        nIonTypes = get_eint(&inp, "iontypes", 1, wi);
        if (nIonTypes < 1)
        {
            wi->addError("You need to provide at least one ion type for position exchanges.");
        }
        ir->swap->ngrp = nIonTypes + static_cast<int>(SwapGroupSplittingType::Count);
        snew(ir->swap->grp, ir->swap->ngrp);
        for (i = 0; i < ir->swap->ngrp; i++)
        {
            snew(ir->swap->grp[i].molname, STRLEN);
        }
        printStringNoNewline(&inp,
                             "Two index groups that contain the compartment-partitioning atoms");
        setStringEntry(&inp,
                       "split-group0",
                       ir->swap->grp[static_cast<int>(SwapGroupSplittingType::Split0)].molname,
                       nullptr);
        setStringEntry(&inp,
                       "split-group1",
                       ir->swap->grp[static_cast<int>(SwapGroupSplittingType::Split1)].molname,
                       nullptr);
        printStringNoNewline(&inp,
                             "Use center of mass of split groups (yes/no), otherwise center of "
                             "geometry is used");
        ir->swap->massw_split[0] = (getEnum<Boolean>(&inp, "massw-split0", wi) != Boolean::No);
        ir->swap->massw_split[1] = (getEnum<Boolean>(&inp, "massw-split1", wi) != Boolean::No);

        printStringNoNewline(&inp, "Name of solvent molecules");
        setStringEntry(&inp,
                       "solvent-group",
                       ir->swap->grp[static_cast<int>(SwapGroupSplittingType::Solvent)].molname,
                       nullptr);

        printStringNoNewline(&inp,
                             "Split cylinder: radius, upper and lower extension (nm) (this will "
                             "define the channels)");
        printStringNoNewline(&inp,
                             "Note that the split cylinder settings do not have an influence on "
                             "the swapping protocol,");
        printStringNoNewline(
                &inp,
                "however, if correctly defined, the permeation events are recorded per channel");
        ir->swap->cyl0r = get_ereal(&inp, "cyl0-r", 2.0, wi);
        ir->swap->cyl0u = get_ereal(&inp, "cyl0-up", 1.0, wi);
        ir->swap->cyl0l = get_ereal(&inp, "cyl0-down", 1.0, wi);
        ir->swap->cyl1r = get_ereal(&inp, "cyl1-r", 2.0, wi);
        ir->swap->cyl1u = get_ereal(&inp, "cyl1-up", 1.0, wi);
        ir->swap->cyl1l = get_ereal(&inp, "cyl1-down", 1.0, wi);

        printStringNoNewline(
                &inp,
                "Average the number of ions per compartment over these many swap attempt steps");
        ir->swap->nAverage = get_eint(&inp, "coupl-steps", 10, wi);

        printStringNoNewline(
                &inp, "Names of the ion types that can be exchanged with solvent molecules,");
        printStringNoNewline(
                &inp, "and the requested number of ions of this type in compartments A and B");
        printStringNoNewline(&inp, "-1 means fix the numbers as found in step 0");
        for (i = 0; i < nIonTypes; i++)
        {
            int ig = static_cast<int>(SwapGroupSplittingType::Count) + i;

            sprintf(buf, "iontype%d-name", i);
            setStringEntry(&inp, buf, ir->swap->grp[ig].molname, nullptr);
            sprintf(buf, "iontype%d-in-A", i);
            ir->swap->grp[ig].nmolReq[0] = get_eint(&inp, buf, -1, wi);
            sprintf(buf, "iontype%d-in-B", i);
            ir->swap->grp[ig].nmolReq[1] = get_eint(&inp, buf, -1, wi);
        }

        printStringNoNewline(
                &inp,
                "By default (i.e. bulk offset = 0.0), ion/water exchanges happen between layers");
        printStringNoNewline(
                &inp,
                "at maximum distance (= bulk concentration) to the split group layers. However,");
        printStringNoNewline(&inp,
                             "an offset b (-1.0 < b < +1.0) can be specified to offset the bulk "
                             "layer from the middle at 0.0");
        printStringNoNewline(&inp,
                             "towards one of the compartment-partitioning layers (at +/- 1.0).");
        ir->swap->bulkOffset[0] = get_ereal(&inp, "bulk-offsetA", 0.0, wi);
        ir->swap->bulkOffset[1] = get_ereal(&inp, "bulk-offsetB", 0.0, wi);
        if (!(ir->swap->bulkOffset[0] > -1.0 && ir->swap->bulkOffset[0] < 1.0)
            || !(ir->swap->bulkOffset[1] > -1.0 && ir->swap->bulkOffset[1] < 1.0))
        {
            wi->addError("Bulk layer offsets must be > -1.0 and < 1.0 !");
        }

        printStringNoNewline(
                &inp, "Start to swap ions if threshold difference to requested count is reached");
        ir->swap->threshold = get_ereal(&inp, "threshold", 1.0, wi);
    }

    /* AdResS is no longer supported, but we need grompp to be able to
       refuse to process old .mdp files that used it. */
    ir->bAdress = (get_eeenum(&inp, "adress", no_names, wi) != 0);

    /* User defined thingies */
    printStringNewline(&inp, "User defined thingies");
    setStringEntry(&inp, "user1-grps", inputrecStrings->user1, nullptr);
    setStringEntry(&inp, "user2-grps", inputrecStrings->user2, nullptr);
    ir->userint1  = get_eint(&inp, "userint1", 0, wi);
    ir->userint2  = get_eint(&inp, "userint2", 0, wi);
    ir->userint3  = get_eint(&inp, "userint3", 0, wi);
    ir->userint4  = get_eint(&inp, "userint4", 0, wi);
    ir->userreal1 = get_ereal(&inp, "userreal1", 0, wi);
    ir->userreal2 = get_ereal(&inp, "userreal2", 0, wi);
    ir->userreal3 = get_ereal(&inp, "userreal3", 0, wi);
    ir->userreal4 = get_ereal(&inp, "userreal4", 0, wi);
#undef CTYPE

    if (mdparout)
    {
        gmx::TextOutputFile stream(mdparout);

        // Set gen-seed line to actual value instead of -1
        if (opts->bMadeSeed)
        {
            int ii         = search_einp(inp, "gen-seed");
            inp[ii].value_ = std::to_string(opts->seed);
        }

        write_inpfile(&stream, mdparout, &inp, FALSE, writeMdpHeader, wi);

        // Transform module data into a flat key-value tree for output.
        gmx::KeyValueTreeBuilder       builder;
        gmx::KeyValueTreeObjectBuilder builderObject = builder.rootObject();
        mdModules->buildMdpOutput(&builderObject);
        {
            gmx::TextWriter writer(&stream);
            writeKeyValueTreeAsMdp(&writer, builder.build());
        }
        stream.close();
    }

    /* Process options if necessary */
    for (m = 0; m < 2; m++)
    {
        for (i = 0; i < 2 * DIM; i++)
        {
            dumdub[m][i] = 0.0;
        }
        if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
        {
            switch (ir->pressureCouplingOptions.epct)
            {
                case PressureCouplingType::Isotropic:
                    if (sscanf(dumstr[m], "%lf", &(dumdub[m][XX])) != 1)
                    {
                        wi->addError(

                                "Pressure coupling incorrect number of values (I need exactly 1)");
                    }
                    dumdub[m][YY] = dumdub[m][ZZ] = dumdub[m][XX];
                    break;
                case PressureCouplingType::SemiIsotropic:
                case PressureCouplingType::SurfaceTension:
                    if (sscanf(dumstr[m], "%lf%lf", &(dumdub[m][XX]), &(dumdub[m][ZZ])) != 2)
                    {
                        wi->addError(

                                "Pressure coupling incorrect number of values (I need exactly 2)");
                    }
                    dumdub[m][YY] = dumdub[m][XX];
                    break;
                case PressureCouplingType::Anisotropic:
                    if (sscanf(dumstr[m],
                               "%lf%lf%lf%lf%lf%lf",
                               &(dumdub[m][XX]),
                               &(dumdub[m][YY]),
                               &(dumdub[m][ZZ]),
                               &(dumdub[m][3]),
                               &(dumdub[m][4]),
                               &(dumdub[m][5]))
                        != 6)
                    {
                        wi->addError(

                                "Pressure coupling incorrect number of values (I need exactly 6)");
                    }
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Pressure coupling type %s not implemented yet",
                              enumValueToString(ir->pressureCouplingOptions.epct));
            }
        }
    }
    clear_mat(ir->pressureCouplingOptions.ref_p);
    clear_mat(ir->pressureCouplingOptions.compress);
    for (i = 0; i < DIM; i++)
    {
        ir->pressureCouplingOptions.ref_p[i][i]    = dumdub[1][i];
        ir->pressureCouplingOptions.compress[i][i] = dumdub[0][i];
    }
    if (ir->pressureCouplingOptions.epct == PressureCouplingType::Anisotropic)
    {
        ir->pressureCouplingOptions.ref_p[XX][YY] = dumdub[1][3];
        ir->pressureCouplingOptions.ref_p[XX][ZZ] = dumdub[1][4];
        ir->pressureCouplingOptions.ref_p[YY][ZZ] = dumdub[1][5];
        if (ir->pressureCouplingOptions.ref_p[XX][YY] != 0 && ir->pressureCouplingOptions.ref_p[XX][ZZ] != 0
            && ir->pressureCouplingOptions.ref_p[YY][ZZ] != 0)
        {
            wi->addWarning(
                    "All off-diagonal reference pressures are non-zero. Are you sure you want to "
                    "apply a threefold shear stress?\n");
        }
        ir->pressureCouplingOptions.compress[XX][YY] = dumdub[0][3];
        ir->pressureCouplingOptions.compress[XX][ZZ] = dumdub[0][4];
        ir->pressureCouplingOptions.compress[YY][ZZ] = dumdub[0][5];
        for (i = 0; i < DIM; i++)
        {
            for (m = 0; m < i; m++)
            {
                ir->pressureCouplingOptions.ref_p[i][m] = ir->pressureCouplingOptions.ref_p[m][i];
                ir->pressureCouplingOptions.compress[i][m] = ir->pressureCouplingOptions.compress[m][i];
            }
        }
    }

    if (ir->comm_mode == ComRemovalAlgorithm::No)
    {
        ir->nstcomm = 0;
    }

    opts->couple_moltype = nullptr;
    if (strlen(inputrecStrings->couple_moltype) > 0)
    {
        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            opts->couple_moltype = gmx_strdup(inputrecStrings->couple_moltype);
            if (opts->couple_lam0 == opts->couple_lam1)
            {
                wi->addWarning("The lambda=0 and lambda=1 states for coupling are identical");
            }
            if (!EI_SD(ir->eI) && ir->eI != IntegrationAlgorithm::BD
                && (opts->couple_lam0 == ecouplamNONE || opts->couple_lam1 == ecouplamNONE))
            {
                wi->addNote(

                        "For proper sampling of the (nearly) decoupled state, stochastic dynamics "
                        "should be used");
            }
        }
        else
        {
            wi->addNote(
                    "Free energy is turned off, so we will not decouple the molecule listed "
                    "in your input.");
        }
    }
    /* FREE ENERGY AND EXPANDED ENSEMBLE OPTIONS */
    if (ir->efep != FreeEnergyPerturbationType::No)
    {
        if (fep->delta_lambda != 0)
        {
            ir->efep = FreeEnergyPerturbationType::SlowGrowth;
        }
    }

    if (fep->edHdLPrintEnergy == FreeEnergyPrintEnergy::Yes)
    {
        fep->edHdLPrintEnergy = FreeEnergyPrintEnergy::Total;
        wi->addNote(
                "Old option for dhdl-print-energy given: "
                "changing \"yes\" to \"total\"\n");
    }

    if (ir->bSimTemp && (fep->edHdLPrintEnergy == FreeEnergyPrintEnergy::No))
    {
        /* always print out the energy to dhdl if we are doing
           expanded ensemble, since we need the total energy for
           analysis if the temperature is changing. In some
           conditions one may only want the potential energy, so
           we will allow that if the appropriate mdp setting has
           been enabled. Otherwise, total it is:
         */
        fep->edHdLPrintEnergy = FreeEnergyPrintEnergy::Total;
    }

    if ((ir->efep != FreeEnergyPerturbationType::No) || ir->bSimTemp)
    {
        ir->bExpanded = FALSE;
        if ((ir->efep == FreeEnergyPerturbationType::Expanded) || ir->bSimTemp)
        {
            ir->bExpanded = TRUE;
        }
        do_fep_params(ir, inputrecStrings->fep_lambda, inputrecStrings->lambda_weights, wi);
        if (ir->bSimTemp) /* done after fep params */
        {
            do_simtemp_params(ir);
        }

        /* Because sc-coul (=FALSE by default) only acts on the lambda state
         * setup and not on the old way of specifying the free-energy setup,
         * we should check for using soft-core when not needed, since that
         * can complicate the sampling significantly.
         * Note that we only check for the automated coupling setup.
         * If the (advanced) user does FEP through manual topology changes,
         * this check will not be triggered.
         */
        if (ir->efep != FreeEnergyPerturbationType::No && ir->fepvals->n_lambda == 0
            && ir->fepvals->sc_alpha != 0
            && (couple_lambda_has_vdw_on(opts->couple_lam0) && couple_lambda_has_vdw_on(opts->couple_lam1)))
        {
            wi->addWarning(
                    "You are using soft-core interactions while the Van der Waals interactions are "
                    "not decoupled (note that the sc-coul option is only active when using lambda "
                    "states). Although this will not lead to errors, you will need much more "
                    "sampling than without soft-core interactions. Consider using sc-alpha=0.");
        }
    }
    else
    {
        ir->fepvals->n_lambda = 0;
    }

    /* WALL PARAMETERS */

    do_wall_params(ir, inputrecStrings->wall_atomtype, inputrecStrings->wall_density, opts, wi);

    /* ORIENTATION RESTRAINT PARAMETERS */

    if (opts->bOrire && gmx::splitString(inputrecStrings->orirefitgrp).size() != 1)
    {
        wi->addError("ERROR: Need one orientation restraint fit group\n");
    }

    /* DEFORMATION PARAMETERS */

    clear_mat(ir->deform);
    for (i = 0; i < 6; i++)
    {
        dumdub[0][i] = 0;
    }

    double gmx_unused canary;
    int               ndeform = sscanf(inputrecStrings->deform,
                         "%lf %lf %lf %lf %lf %lf %lf",
                         &(dumdub[0][0]),
                         &(dumdub[0][1]),
                         &(dumdub[0][2]),
                         &(dumdub[0][3]),
                         &(dumdub[0][4]),
                         &(dumdub[0][5]),
                         &canary);

    if (strlen(inputrecStrings->deform) > 0 && ndeform != 6)
    {
        wi->addError(gmx::formatString(
                "Cannot parse exactly 6 box deformation velocities from string '%s'",
                inputrecStrings->deform));
    }
    for (i = 0; i < 3; i++)
    {
        ir->deform[i][i] = dumdub[0][i];
    }
    ir->deform[YY][XX] = dumdub[0][3];
    ir->deform[ZZ][XX] = dumdub[0][4];
    ir->deform[ZZ][YY] = dumdub[0][5];
    if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
    {
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j <= i; j++)
            {
                if (ir->deform[i][j] != 0 && ir->pressureCouplingOptions.compress[i][j] != 0)
                {
                    wi->addError("A box element has deform set and compressibility > 0");
                }
            }
        }
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < i; j++)
            {
                if (ir->deform[i][j] != 0)
                {
                    for (m = j; m < DIM; m++)
                    {
                        if (ir->pressureCouplingOptions.compress[m][j] != 0)
                        {
                            sprintf(warn_buf,
                                    "An off-diagonal box element has deform set while "
                                    "compressibility > 0 for the same component of another box "
                                    "vector, this might lead to spurious periodicity effects.");
                            wi->addWarning(warn_buf);
                        }
                    }
                }
            }
        }
    }

    /* Ion/water position swapping checks */
    if (ir->eSwapCoords != SwapType::No)
    {
        if (ir->swap->nstswap < 1)
        {
            wi->addError("swap_frequency must be 1 or larger when ion swapping is requested");
        }
        if (ir->swap->nAverage < 1)
        {
            wi->addError("coupl_steps must be 1 or larger.\n");
        }
        if (ir->swap->threshold < 1.0)
        {
            wi->addError("Ion count threshold must be at least 1.\n");
        }
    }

    /* Set up MTS levels, this needs to happen before checking AWH parameters */
    if (ir->useMts)
    {
        std::vector<std::string> errorMessages;
        ir->mtsLevels = gmx::setupMtsLevels(opts->mtsOpts, &errorMessages);

        for (const auto& errorMessage : errorMessages)
        {
            wi->addError(errorMessage);
        }
    }

    if (ir->bDoAwh)
    {
        gmx::checkAwhParams(*ir->awhParams, *ir, wi);
    }

    sfree(dumstr[0]);
    sfree(dumstr[1]);
}

int getGroupIndex(const std::string& s, gmx::ArrayRef<const IndexGroup> indexGroups)
{
    for (int i = 0; i < gmx::ssize(indexGroups); i++)
    {
        if (gmx_strcasecmp(s.c_str(), indexGroups[i].name.c_str()) == 0)
        {
            return i;
        }
    }

    gmx_fatal(FARGS,
              "Group %s referenced in the .mdp file was not found in the list of index groups.\n"
              "Group names must match either [moleculetype] names or custom index group\n"
              "names, in which case you must supply an index file to the '-n' option\n"
              "of grompp.",
              s.c_str());
}

static void atomGroupRangeValidation(const int natoms, gmx::ArrayRef<const int> particleIndices)
{
    /* Now go over the atoms in the group */
    for (const int aj : particleIndices)
    {
        /* Range checking */
        if ((aj < 0) || (aj >= natoms))
        {
            gmx_fatal(FARGS, "Invalid atom number %d in indexfile", aj + 1);
        }
    }
}

/*! \brief Creates the groups of atom indices for group type \p gtype
 *
 * \param[in] natoms  The total number of atoms in the system
 * \param[in,out] groups  Index \p gtype in this list of list of groups will be set
 * \param[in] groupsFromMdpFile  The list of group names set for \p gtype in the mdp file
 * \param[in] indexGroups The list of all available index groups
 * \param[in] gtype       The group type to creates groups for
 * \param[in] restnm      The index of rest group name in \p gnames
 * \param[in] coverage    How to treat coverage of all atoms in the system
 * \param[in] bVerbose    Whether to print when we make a rest group
 * \param[in,out] wi      List of warnings
 *
 * \returns whether all atoms have been assigned to a group
 */
static bool do_numbering(const int                        natoms,
                         SimulationGroups*                groups,
                         gmx::ArrayRef<const std::string> groupsFromMdpFile,
                         gmx::ArrayRef<const IndexGroup>  indexGroups,
                         const SimulationAtomGroupType    gtype,
                         const int                        restnm,
                         const GroupCoverage              coverage,
                         const bool                       bVerbose,
                         WarningHandler*                  wi)
{
    unsigned short*   cbuf;
    AtomGroupIndices* grps = &(groups->groups[gtype]);
    int               ntot = 0;
    const char*       title;
    char              warn_buf[STRLEN];

    title = shortName(gtype);

    snew(cbuf, natoms);
    /* Mark all id's as not set */
    for (int i = 0; (i < natoms); i++)
    {
        cbuf[i] = NOGID;
    }

    for (int i = 0; i != groupsFromMdpFile.ssize(); ++i)
    {
        /* Lookup the group name in the block structure */
        const int gid = getGroupIndex(groupsFromMdpFile[i], indexGroups);
        if ((coverage != GroupCoverage::OneGroup) || (i == 0))
        {
            grps->emplace_back(gid);
        }
        gmx::ArrayRef<const int> indexGroup = indexGroups[gid].particleIndices;
        atomGroupRangeValidation(natoms, indexGroup);
        /* Now go over the atoms in the group */
        for (const int aj : indexGroup)
        {
            /* Lookup up the old group number */
            const int ognr = cbuf[aj];
            if (ognr != NOGID)
            {
                gmx_fatal(FARGS, "Atom %d in multiple %s groups (%d and %d)", aj + 1, title, ognr + 1, i + 1);
            }
            else
            {
                /* Store the group number in buffer */
                if (coverage == GroupCoverage::OneGroup)
                {
                    cbuf[aj] = 0;
                }
                else
                {
                    cbuf[aj] = i;
                }
                ntot++;
            }
        }
    }

    /* Now check whether we have done all atoms */
    if (ntot != natoms)
    {
        if (coverage == GroupCoverage::All)
        {
            gmx_fatal(FARGS, "%d atoms are not part of any of the %s groups", natoms - ntot, title);
        }
        else if (coverage == GroupCoverage::Partial)
        {
            sprintf(warn_buf, "%d atoms are not part of any of the %s groups", natoms - ntot, title);
            wi->addNote(warn_buf);
        }
        /* Assign all atoms currently unassigned to a rest group */
        for (int j = 0; (j < natoms); j++)
        {
            if (cbuf[j] == NOGID)
            {
                cbuf[j] = grps->size();
            }
        }
        if (coverage != GroupCoverage::Partial)
        {
            if (bVerbose)
            {
                fprintf(stderr, "Making dummy/rest group for %s containing %d elements\n", title, natoms - ntot);
            }
            /* Add group name "rest" */
            grps->emplace_back(restnm);

            /* Assign the rest name to all atoms not currently assigned to a group */
            for (int j = 0; (j < natoms); j++)
            {
                if (cbuf[j] == NOGID)
                {
                    // group size was not updated before this here, so need to use -1.
                    cbuf[j] = grps->size() - 1;
                }
            }
        }
    }

    if (grps->size() == 1 && (ntot == 0 || ntot == natoms))
    {
        /* All atoms are part of one (or no) group, no index required */
        groups->groupNumbers[gtype].clear();
    }
    else
    {
        for (int j = 0; (j < natoms); j++)
        {
            groups->groupNumbers[gtype].emplace_back(cbuf[j]);
        }
    }

    sfree(cbuf);

    return ntot == natoms;
}

static void calc_nrdf(const gmx_mtop_t* mtop, t_inputrec* ir, gmx::ArrayRef<const std::string> gnames)
{
    t_grpopts*     opts;
    pull_params_t* pull;
    int            natoms, imin, jmin;
    int *          nrdf2, *na_vcm, na_tot;
    double *       nrdf_tc, *nrdf_vcm, nrdf_uc, *nrdf_vcm_sub;
    ivec*          dof_vcm;
    int            as;

    /* Calculate nrdf.
     * First calc 3xnr-atoms for each group
     * then subtract half a degree of freedom for each constraint
     *
     * Only atoms and nuclei contribute to the degrees of freedom...
     */

    opts = &ir->opts;

    const SimulationGroups& groups = mtop->groups;
    natoms                         = mtop->natoms;

    /* Allocate one more for a possible rest group */
    /* We need to sum degrees of freedom into doubles,
     * since floats give too low nrdf's above 3 million atoms.
     */
    snew(nrdf_tc, groups.groups[SimulationAtomGroupType::TemperatureCoupling].size() + 1);
    snew(nrdf_vcm, groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval].size() + 1);
    snew(dof_vcm, groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval].size() + 1);
    snew(na_vcm, groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval].size() + 1);
    snew(nrdf_vcm_sub, groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval].size() + 1);

    for (gmx::Index i = 0; i < gmx::ssize(groups.groups[SimulationAtomGroupType::TemperatureCoupling]); i++)
    {
        nrdf_tc[i] = 0;
    }
    for (gmx::Index i = 0;
         i < gmx::ssize(groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval]) + 1;
         i++)
    {
        nrdf_vcm[i] = 0;
        clear_ivec(dof_vcm[i]);
        na_vcm[i]       = 0;
        nrdf_vcm_sub[i] = 0;
    }
    snew(nrdf2, natoms);
    for (const AtomProxy atomP : AtomRange(*mtop))
    {
        const t_atom& local = atomP.atom();
        int           i     = atomP.globalAtomNumber();
        nrdf2[i]            = 0;
        if (local.ptype == ParticleType::Atom || local.ptype == ParticleType::Nucleus)
        {
            int g = getGroupType(groups, SimulationAtomGroupType::Freeze, i);
            for (int d = 0; d < DIM; d++)
            {
                if (opts->nFreeze[g][d] == 0)
                {
                    /* Add one DOF for particle i (counted as 2*1) */
                    nrdf2[i] += 2;
                    /* VCM group i has dim d as a DOF */
                    dof_vcm[getGroupType(groups, SimulationAtomGroupType::MassCenterVelocityRemoval, i)][d] =
                            1;
                }
            }
            nrdf_tc[getGroupType(groups, SimulationAtomGroupType::TemperatureCoupling, i)] +=
                    0.5 * nrdf2[i];
            nrdf_vcm[getGroupType(groups, SimulationAtomGroupType::MassCenterVelocityRemoval, i)] +=
                    0.5 * nrdf2[i];
        }
    }

    as = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const gmx_moltype_t& molt = mtop->moltype[molb.type];
        const t_atom*        atom = molt.atoms.atom;
        for (int mol = 0; mol < molb.nmol; mol++)
        {
            for (int ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
            {
                gmx::ArrayRef<const int> ia = molt.ilist[ftype].iatoms;
                for (int i = 0; i < molt.ilist[ftype].size();)
                {
                    /* Subtract degrees of freedom for the constraints,
                     * if the particles still have degrees of freedom left.
                     * If one of the particles is a vsite or a shell, then all
                     * constraint motion will go there, but since they do not
                     * contribute to the constraints the degrees of freedom do not
                     * change.
                     */
                    int ai = as + ia[i + 1];
                    int aj = as + ia[i + 2];
                    if (((atom[ia[i + 1]].ptype == ParticleType::Nucleus)
                         || (atom[ia[i + 1]].ptype == ParticleType::Atom))
                        && ((atom[ia[i + 2]].ptype == ParticleType::Nucleus)
                            || (atom[ia[i + 2]].ptype == ParticleType::Atom)))
                    {
                        if (nrdf2[ai] > 0)
                        {
                            jmin = 1;
                        }
                        else
                        {
                            jmin = 2;
                        }
                        if (nrdf2[aj] > 0)
                        {
                            imin = 1;
                        }
                        else
                        {
                            imin = 2;
                        }
                        imin = std::min(imin, nrdf2[ai]);
                        jmin = std::min(jmin, nrdf2[aj]);
                        nrdf2[ai] -= imin;
                        nrdf2[aj] -= jmin;
                        nrdf_tc[getGroupType(groups, SimulationAtomGroupType::TemperatureCoupling, ai)] -=
                                0.5 * imin;
                        nrdf_tc[getGroupType(groups, SimulationAtomGroupType::TemperatureCoupling, aj)] -=
                                0.5 * jmin;
                        nrdf_vcm[getGroupType(groups, SimulationAtomGroupType::MassCenterVelocityRemoval, ai)] -=
                                0.5 * imin;
                        nrdf_vcm[getGroupType(groups, SimulationAtomGroupType::MassCenterVelocityRemoval, aj)] -=
                                0.5 * jmin;
                    }
                    i += interaction_function[ftype].nratoms + 1;
                }
            }
            gmx::ArrayRef<const int> ia = molt.ilist[F_SETTLE].iatoms;
            for (int i = 0; i < molt.ilist[F_SETTLE].size();)
            {
                /* Subtract 1 dof from every atom in the SETTLE */
                for (int j = 0; j < 3; j++)
                {
                    int ai = as + ia[i + 1 + j];
                    imin   = std::min(2, nrdf2[ai]);
                    nrdf2[ai] -= imin;
                    nrdf_tc[getGroupType(groups, SimulationAtomGroupType::TemperatureCoupling, ai)] -=
                            0.5 * imin;
                    nrdf_vcm[getGroupType(groups, SimulationAtomGroupType::MassCenterVelocityRemoval, ai)] -=
                            0.5 * imin;
                }
                i += 4;
            }
            as += molt.atoms.nr;
        }
    }

    if (ir->bPull)
    {
        /* Correct nrdf for the COM constraints.
         * We correct using the TC and VCM group of the first atom
         * in the reference and pull group. If atoms in one pull group
         * belong to different TC or VCM groups it is anyhow difficult
         * to determine the optimal nrdf assignment.
         */
        pull = ir->pull.get();

        for (int i = 0; i < pull->ncoord; i++)
        {
            if (pull->coord[i].eType != PullingAlgorithm::Constraint)
            {
                continue;
            }

            imin = 1;

            for (int j = 0; j < 2; j++)
            {
                const t_pull_group* pgrp;

                pgrp = &pull->group[pull->coord[i].group[j]];

                if (!pgrp->ind.empty())
                {
                    /* Subtract 1/2 dof from each group */
                    int ai = pgrp->ind[0];
                    nrdf_tc[getGroupType(groups, SimulationAtomGroupType::TemperatureCoupling, ai)] -=
                            0.5 * imin;
                    nrdf_vcm[getGroupType(groups, SimulationAtomGroupType::MassCenterVelocityRemoval, ai)] -=
                            0.5 * imin;
                    if (nrdf_tc[getGroupType(groups, SimulationAtomGroupType::TemperatureCoupling, ai)] < 0)
                    {
                        gmx_fatal(FARGS,
                                  "Center of mass pulling constraints caused the number of degrees "
                                  "of freedom for temperature coupling group %s to be negative",
                                  gnames[groups.groups[SimulationAtomGroupType::TemperatureCoupling][getGroupType(
                                                 groups, SimulationAtomGroupType::TemperatureCoupling, ai)]]
                                          .c_str());
                    }
                }
                else
                {
                    /* We need to subtract the whole DOF from group j=1 */
                    imin += 1;
                }
            }
        }
    }

    if (ir->nstcomm != 0)
    {
        GMX_RELEASE_ASSERT(!groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval].empty(),
                           "Expect at least one group when removing COM motion");

        /* We remove COM motion up to dim ndof_com() */
        const int ndim_rm_vcm = ndof_com(ir);

        /* Subtract ndim_rm_vcm (or less with frozen dimensions) from
         * the number of degrees of freedom in each vcm group when COM
         * translation is removed and 6 when rotation is removed as well.
         * Note that we do not and should not include the rest group here.
         */
        for (gmx::Index j = 0;
             j < gmx::ssize(groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval]);
             j++)
        {
            switch (ir->comm_mode)
            {
                case ComRemovalAlgorithm::Linear:
                case ComRemovalAlgorithm::LinearAccelerationCorrection:
                    nrdf_vcm_sub[j] = 0;
                    for (int d = 0; d < ndim_rm_vcm; d++)
                    {
                        if (dof_vcm[j][d])
                        {
                            nrdf_vcm_sub[j]++;
                        }
                    }
                    break;
                case ComRemovalAlgorithm::Angular: nrdf_vcm_sub[j] = 6; break;
                default: gmx_incons("Checking comm_mode");
            }
        }

        for (gmx::Index i = 0;
             i < gmx::ssize(groups.groups[SimulationAtomGroupType::TemperatureCoupling]);
             i++)
        {
            /* Count the number of atoms of TC group i for every VCM group */
            for (gmx::Index j = 0;
                 j < gmx::ssize(groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval]) + 1;
                 j++)
            {
                na_vcm[j] = 0;
            }
            na_tot = 0;
            for (int ai = 0; ai < natoms; ai++)
            {
                if (getGroupType(groups, SimulationAtomGroupType::TemperatureCoupling, ai) == i)
                {
                    na_vcm[getGroupType(groups, SimulationAtomGroupType::MassCenterVelocityRemoval, ai)]++;
                    na_tot++;
                }
            }
            /* Correct for VCM removal according to the fraction of each VCM
             * group present in this TC group.
             */
            nrdf_uc    = nrdf_tc[i];
            nrdf_tc[i] = 0;
            for (gmx::Index j = 0;
                 j < gmx::ssize(groups.groups[SimulationAtomGroupType::MassCenterVelocityRemoval]) + 1;
                 j++)
            {
                if (nrdf_vcm[j] > nrdf_vcm_sub[j])
                {
                    nrdf_tc[i] += nrdf_uc * (static_cast<double>(na_vcm[j]) / static_cast<double>(na_tot))
                                  * (nrdf_vcm[j] - nrdf_vcm_sub[j]) / nrdf_vcm[j];
                }
            }
        }
    }
    for (int i = 0; (i < gmx::ssize(groups.groups[SimulationAtomGroupType::TemperatureCoupling])); i++)
    {
        opts->nrdf[i] = nrdf_tc[i];
        if (opts->nrdf[i] < 0)
        {
            opts->nrdf[i] = 0;
        }
        fprintf(stderr,
                "Number of degrees of freedom in T-Coupling group %s is %.2f\n",
                gnames[groups.groups[SimulationAtomGroupType::TemperatureCoupling][i]].c_str(),
                opts->nrdf[i]);
    }

    sfree(nrdf2);
    sfree(nrdf_tc);
    sfree(nrdf_vcm);
    sfree(dof_vcm);
    sfree(na_vcm);
    sfree(nrdf_vcm_sub);
}

static bool do_egp_flag(t_inputrec* ir, SimulationGroups* groups, const char* option, const char* val, int flag)
{
    /* The maximum number of energy group pairs would be MAXPTR*(MAXPTR+1)/2.
     * But since this is much larger than STRLEN, such a line can not be parsed.
     * The real maximum is the number of names that fit in a string: STRLEN/2.
     */
    int  j, k, nr;
    bool bSet;

    auto names = gmx::splitString(val);
    if (names.size() % 2 != 0)
    {
        gmx_fatal(FARGS, "The number of groups for %s is odd", option);
    }
    nr   = groups->groups[SimulationAtomGroupType::EnergyOutput].size();
    bSet = FALSE;
    for (size_t i = 0; i < names.size() / 2; i++)
    {
        // TODO this needs to be replaced by a solution using std::find_if
        j = 0;
        while ((j < nr)
               && gmx_strcasecmp(
                       names[2 * i].c_str(),
                       *(groups->groupNames[groups->groups[SimulationAtomGroupType::EnergyOutput][j]])))
        {
            j++;
        }
        if (j == nr)
        {
            gmx_fatal(FARGS, "%s in %s is not an energy group\n", names[2 * i].c_str(), option);
        }
        k = 0;
        while ((k < nr)
               && gmx_strcasecmp(
                       names[2 * i + 1].c_str(),
                       *(groups->groupNames[groups->groups[SimulationAtomGroupType::EnergyOutput][k]])))
        {
            k++;
        }
        if (k == nr)
        {
            gmx_fatal(FARGS, "%s in %s is not an energy group\n", names[2 * i + 1].c_str(), option);
        }
        if ((j < nr) && (k < nr))
        {
            ir->opts.egp_flags[nr * j + k] |= flag;
            ir->opts.egp_flags[nr * k + j] |= flag;
            bSet = TRUE;
        }
    }

    return bSet;
}


static void make_swap_groups(t_swapcoords* swap, gmx::ArrayRef<const IndexGroup> indexGroups)
{
    int          ig = -1, i = 0, gind;
    t_swapGroup* swapg;


    /* Just a quick check here, more thorough checks are in mdrun */
    if (strcmp(swap->grp[static_cast<int>(SwapGroupSplittingType::Split0)].molname,
               swap->grp[static_cast<int>(SwapGroupSplittingType::Split1)].molname)
        == 0)
    {
        gmx_fatal(FARGS,
                  "The split groups can not both be '%s'.",
                  swap->grp[static_cast<int>(SwapGroupSplittingType::Split0)].molname);
    }

    /* Get the index atoms of the split0, split1, solvent, and swap groups */
    for (ig = 0; ig < swap->ngrp; ig++)
    {
        swapg      = &swap->grp[ig];
        gind       = getGroupIndex(swap->grp[ig].molname, indexGroups);
        swapg->nat = gmx::ssize(indexGroups[gind].particleIndices);

        if (swapg->nat > 0)
        {
            fprintf(stderr,
                    "%s group '%s' contains %d atoms.\n",
                    ig < 3 ? enumValueToString(static_cast<SwapGroupSplittingType>(ig)) : "Swap",
                    swap->grp[ig].molname,
                    swapg->nat);
            snew(swapg->ind, swapg->nat);
            for (i = 0; i < swapg->nat; i++)
            {
                swapg->ind[i] = indexGroups[gind].particleIndices[i];
            }
        }
        else
        {
            gmx_fatal(FARGS, "Swap group %s does not contain any atoms.", swap->grp[ig].molname);
        }
    }
}


static void make_IMD_group(t_IMD* IMDgroup, const char* IMDgname, gmx::ArrayRef<const IndexGroup> indexGroups)
{
    int ig, i;


    ig            = getGroupIndex(IMDgname, indexGroups);
    IMDgroup->nat = gmx::ssize(indexGroups[ig].particleIndices);

    if (IMDgroup->nat > 0)
    {
        fprintf(stderr,
                "Group '%s' with %d atoms can be activated for interactive molecular dynamics "
                "(IMD).\n",
                IMDgname,
                IMDgroup->nat);
        snew(IMDgroup->ind, IMDgroup->nat);
        for (i = 0; i < IMDgroup->nat; i++)
        {
            IMDgroup->ind[i] = indexGroups[ig].particleIndices[i];
        }
    }
}

/* Checks whether atoms are both part of a COM removal group and frozen.
 * If a fully frozen atom is part of a COM removal group, it is removed
 * from the COM removal group. A note is issued if such atoms are present.
 * A warning is issued for atom with one or two dimensions frozen that
 * are part of a COM removal group (mdrun would need to compute COM mass
 * per dimension to handle this correctly).
 * Also issues a warning when non-frozen atoms are not part of a COM
 * removal group while COM removal is active.
 */
static void checkAndUpdateVcmFreezeGroupConsistency(SimulationGroups* groups,
                                                    const int         numAtoms,
                                                    const t_grpopts&  opts,
                                                    WarningHandler*   wi)
{
    const int vcmRestGroup =
            std::max(int(groups->groups[SimulationAtomGroupType::MassCenterVelocityRemoval].size()), 1);

    int numFullyFrozenVcmAtoms     = 0;
    int numPartiallyFrozenVcmAtoms = 0;
    int numNonVcmAtoms             = 0;
    for (int a = 0; a < numAtoms; a++)
    {
        const int freezeGroup   = getGroupType(*groups, SimulationAtomGroupType::Freeze, a);
        int       numFrozenDims = 0;
        for (int d = 0; d < DIM; d++)
        {
            numFrozenDims += opts.nFreeze[freezeGroup][d];
        }

        const int vcmGroup = getGroupType(*groups, SimulationAtomGroupType::MassCenterVelocityRemoval, a);
        if (vcmGroup < vcmRestGroup)
        {
            if (numFrozenDims == DIM)
            {
                /* Do not remove COM motion for this fully frozen atom */
                if (groups->groupNumbers[SimulationAtomGroupType::MassCenterVelocityRemoval].empty())
                {
                    groups->groupNumbers[SimulationAtomGroupType::MassCenterVelocityRemoval].resize(
                            numAtoms, 0);
                }
                groups->groupNumbers[SimulationAtomGroupType::MassCenterVelocityRemoval][a] = vcmRestGroup;
                numFullyFrozenVcmAtoms++;
            }
            else if (numFrozenDims > 0)
            {
                numPartiallyFrozenVcmAtoms++;
            }
        }
        else if (numFrozenDims < DIM)
        {
            numNonVcmAtoms++;
        }
    }

    if (numFullyFrozenVcmAtoms > 0)
    {
        std::string warningText = gmx::formatString(
                "There are %d atoms that are fully frozen and part of COMM removal group(s), "
                "removing these atoms from the COMM removal group(s)",
                numFullyFrozenVcmAtoms);
        wi->addNote(warningText);
    }
    if (numPartiallyFrozenVcmAtoms > 0 && numPartiallyFrozenVcmAtoms < numAtoms)
    {
        std::string warningText = gmx::formatString(
                "There are %d atoms that are frozen along less then %d dimensions and part of COMM "
                "removal group(s), due to limitations in the code these still contribute to the "
                "mass of the COM along frozen dimensions and therefore the COMM correction will be "
                "too small.",
                numPartiallyFrozenVcmAtoms,
                DIM);
        wi->addWarning(warningText);
    }
    if (numNonVcmAtoms > 0)
    {
        std::string warningText = gmx::formatString(
                "%d atoms are not part of any center of mass motion removal group.\n"
                "This may lead to artifacts.\n"
                "In most cases one should use one group for the whole system.",
                numNonVcmAtoms);
        wi->addWarning(warningText);
    }
}

static void processEnsembleTemperature(t_inputrec* ir, const bool allAtomsCoupled, WarningHandler* wi)
{
    if (ir->ensembleTemperatureSetting == EnsembleTemperatureSetting::NotAvailable)
    {
        ir->ensembleTemperature = -1;
    }
    else if (ir->ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant)
    {
        if (ir->ensembleTemperature < 0)
        {
            wi->addError("ensemble-temperature can not be negative");
        }
        else if (integratorHasReferenceTemperature(*ir))
        {
            bool refTEqual = true;
            for (int i = 1; i < ir->opts.ngtc; i++)
            {
                if (ir->opts.ref_t[i] != ir->opts.ref_t[0])
                {
                    refTEqual = false;
                }
            }
            if (refTEqual && ir->ensembleTemperature != ir->opts.ref_t[0])
            {
                wi->addWarning(
                        "The ensemble temperature of the system does not match the reference "
                        "temperature(s) of the T-coupling group(s)");
            }
        }
    }
    else if (ir->ensembleTemperatureSetting == EnsembleTemperatureSetting::Auto)
    {
        if (integratorHasReferenceTemperature(*ir))
        {
            if (!allAtomsCoupled)
            {
                fprintf(stderr,
                        "Not all atoms are temperature coupled: there is no ensemble temperature "
                        "available\n");
                ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::NotAvailable;
                ir->ensembleTemperature        = -1;
            }
            else if (doSimulatedAnnealing(*ir) && ir->opts.ngtc > 1)
            {
                // We could support ensemble temperature if all annealing groups have the same
                // temperature, but that is bug-prone, so we don't implement that.
                fprintf(stderr,
                        "Simulated tempering is used with multiple T-coupling groups: setting the "
                        "ensemble temperature to not available\n");
                ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::NotAvailable;
                ir->ensembleTemperature        = -1;
            }
            else if (doSimulatedAnnealing(*ir) || ir->bSimTemp)
            {
                ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::Variable;
            }
            else
            {
                bool refTEqual = true;
                for (int i = 1; i < ir->opts.ngtc; i++)
                {
                    if (ir->opts.ref_t[i] != ir->opts.ref_t[0])
                    {
                        refTEqual = false;
                    }
                }
                if (refTEqual)
                {
                    ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::Constant;
                    ir->ensembleTemperature        = ir->opts.ref_t[0];
                }
                else
                {
                    ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::NotAvailable;
                    ir->ensembleTemperature        = -1;
                }
            }
        }
        else
        {
            // We do not have an ensemble temperature available
            fprintf(stderr,
                    "The integrator does not provide a ensemble temperature, there is no system "
                    "ensemble temperature\n");
            ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::NotAvailable;
            ir->ensembleTemperature        = -1;
        }
    }
}

void do_index(const char*                                 mdparin,
              const std::optional<std::filesystem::path>& ndx,
              gmx_mtop_t*                                 mtop,
              bool                                        bVerbose,
              const gmx::MDModulesNotifiers&              mdModulesNotifiers,
              t_inputrec*                                 ir,
              WarningHandler*                             wi)
{
    int       natoms;
    t_symtab* symtab;
    t_atoms   atoms_all;
    int       nr;
    real      tau_min;
    int       i, j, k, restnm;
    bool      bExcl, bTable, bAnneal;
    char      warn_buf[STRLEN];

    if (bVerbose)
    {
        fprintf(stderr, "processing index file...\n");
    }
    std::vector<IndexGroup> defaultIndexGroups;
    if (ndx)
    {
        defaultIndexGroups = init_index(ndx.value());
    }
    else
    {
        atoms_all          = gmx_mtop_global_atoms(*mtop);
        defaultIndexGroups = analyse(&atoms_all, false, true);
        done_atom(&atoms_all);
    }

    SimulationGroups* groups = &mtop->groups;
    natoms                   = mtop->natoms;
    symtab                   = &mtop->symtab;

    // We need a temporary list of the group names from the index file plus the rest group
    std::vector<std::string> gnames;
    for (const auto& indexGroup : defaultIndexGroups)
    {
        groups->groupNames.emplace_back(put_symtab(symtab, indexGroup.name.c_str()));
        gnames.emplace_back(*groups->groupNames.back());
    }
    groups->groupNames.emplace_back(put_symtab(symtab, "rest"));
    restnm = groups->groupNames.size() - 1;
    GMX_RELEASE_ASSERT(restnm == gmx::ssize(defaultIndexGroups), "Size of allocations must match");
    gnames.emplace_back(*groups->groupNames.back());

    wi->setFileAndLineNumber(mdparin, -1);

    auto temperatureCouplingTauValues       = gmx::splitString(inputrecStrings->tau_t);
    auto temperatureCouplingReferenceValues = gmx::splitString(inputrecStrings->ref_t);
    auto temperatureCouplingGroupNames      = gmx::splitString(inputrecStrings->tcgrps);
    if (temperatureCouplingTauValues.size() != temperatureCouplingGroupNames.size()
        || temperatureCouplingReferenceValues.size() != temperatureCouplingGroupNames.size())
    {
        gmx_fatal(FARGS,
                  "Invalid T coupling input: %zu groups, %zu ref-t values and "
                  "%zu tau-t values",
                  temperatureCouplingGroupNames.size(),
                  temperatureCouplingReferenceValues.size(),
                  temperatureCouplingTauValues.size());
    }

    const bool useReferenceTemperature = integratorHasReferenceTemperature(*ir);
    const bool allAtomsAreTCoupled =
            do_numbering(natoms,
                         groups,
                         temperatureCouplingGroupNames,
                         defaultIndexGroups,
                         SimulationAtomGroupType::TemperatureCoupling,
                         restnm,
                         useReferenceTemperature ? GroupCoverage::All : GroupCoverage::AllGenerateRest,
                         bVerbose,
                         wi);
    nr            = groups->groups[SimulationAtomGroupType::TemperatureCoupling].size();
    ir->opts.ngtc = nr;
    snew(ir->opts.nrdf, nr);
    snew(ir->opts.tau_t, nr);
    snew(ir->opts.ref_t, nr);
    if (ir->eI == IntegrationAlgorithm::BD && ir->bd_fric == 0)
    {
        fprintf(stderr, "bd-fric=0, so tau-t will be used as the inverse friction constant(s)\n");
    }

    if (useReferenceTemperature)
    {
        if (size_t(nr) != temperatureCouplingReferenceValues.size())
        {
            gmx_fatal(FARGS, "Not enough ref-t and tau-t values!");
        }

        tau_min = 1e20;
        convertReals(wi, temperatureCouplingTauValues, "tau-t", ir->opts.tau_t);
        for (i = 0; (i < nr); i++)
        {
            if ((ir->eI == IntegrationAlgorithm::BD) && ir->opts.tau_t[i] <= 0)
            {
                sprintf(warn_buf,
                        "With integrator %s tau-t should be larger than 0",
                        enumValueToString(ir->eI));
                wi->addError(warn_buf);
            }

            if (ir->etc != TemperatureCoupling::VRescale && ir->opts.tau_t[i] == 0)
            {
                wi->addNote(

                        "tau-t = -1 is the value to signal that a group should not have "
                        "temperature coupling. Treating your use of tau-t = 0 as if you used -1.");
            }

            if (ir->opts.tau_t[i] >= 0)
            {
                tau_min = std::min(tau_min, ir->opts.tau_t[i]);
            }
        }
        if (ir->etc != TemperatureCoupling::No && ir->nsttcouple == -1)
        {
            ir->nsttcouple = ir_optimal_nsttcouple(ir);
        }

        if (EI_VV(ir->eI))
        {
            if ((ir->etc == TemperatureCoupling::NoseHoover)
                && (ir->pressureCouplingOptions.epc == PressureCoupling::Berendsen))
            {
                gmx_fatal(FARGS,
                          "Cannot do Nose-Hoover temperature with Berendsen pressure control with "
                          "md-vv; use either vrescale temperature with berendsen pressure or "
                          "Nose-Hoover temperature with MTTK pressure");
            }
            if (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
            {
                wi->addNote("MTTK coupling is deprecated and will soon be removed");

                if (ir->etc != TemperatureCoupling::NoseHoover)
                {
                    gmx_fatal(FARGS,
                              "Cannot do MTTK pressure coupling without Nose-Hoover temperature "
                              "control");
                }
                else
                {
                    if (ir->pressureCouplingOptions.nstpcouple != ir->nsttcouple)
                    {
                        int mincouple = std::min(ir->pressureCouplingOptions.nstpcouple, ir->nsttcouple);
                        ir->pressureCouplingOptions.nstpcouple = ir->nsttcouple = mincouple;
                        sprintf(warn_buf,
                                "for current Trotter decomposition methods with vv, nsttcouple and "
                                "nstpcouple must be equal.  Both have been reset to "
                                "min(nsttcouple,nstpcouple) = %d",
                                mincouple);
                        wi->addNote(warn_buf);
                    }
                }
            }
        }
        /* velocity verlet with averaged kinetic energy KE = 0.5*(v(t+1/2) - v(t-1/2)) is implemented
           primarily for testing purposes, and does not work with temperature coupling other than 1 */

        if (ETC_ANDERSEN(ir->etc))
        {
            if (ir->nsttcouple != 1)
            {
                ir->nsttcouple = 1;
                sprintf(warn_buf,
                        "Andersen temperature control methods assume nsttcouple = 1; there is no "
                        "need for larger nsttcouple > 1, since no global parameters are computed. "
                        "nsttcouple has been reset to 1");
                wi->addNote(warn_buf);
            }
        }
        const int nstcmin = tcouple_min_integration_steps(ir->etc);
        // V-rescale can act correctly with any coupling interval
        if (nstcmin > 1 && ir->etc != TemperatureCoupling::VRescale)
        {
            if (tau_min / (ir->delta_t * ir->nsttcouple) < nstcmin - 10 * GMX_REAL_EPS)
            {
                sprintf(warn_buf,
                        "For proper integration of the %s thermostat, tau-t (%g) should be at "
                        "least %d times larger than nsttcouple*dt (%g)",
                        enumValueToString(ir->etc),
                        tau_min,
                        nstcmin,
                        ir->nsttcouple * ir->delta_t);
                wi->addWarning(warn_buf);
            }
        }
        convertReals(wi, temperatureCouplingReferenceValues, "ref-t", ir->opts.ref_t);
        for (i = 0; (i < nr); i++)
        {
            if (ir->opts.ref_t[i] < 0)
            {
                gmx_fatal(FARGS, "ref-t for group %d negative", i);
            }
        }
        /* set the lambda mc temperature to the md integrator temperature (which should be defined
           if we are in this conditional) if mc_temp is negative */
        if (ir->expandedvals->mc_temp < 0)
        {
            ir->expandedvals->mc_temp = ir->opts.ref_t[0]; /*for now, set to the first reft */
        }
    }

    /* Simulated annealing for each group. There are nr groups */
    auto simulatedAnnealingGroupNames = gmx::splitString(inputrecStrings->anneal);
    if (simulatedAnnealingGroupNames.size() == 1
        && gmx::equalCaseInsensitive(simulatedAnnealingGroupNames[0], "N", 1))
    {
        simulatedAnnealingGroupNames.resize(0);
    }
    if (!simulatedAnnealingGroupNames.empty() && gmx::ssize(simulatedAnnealingGroupNames) != nr)
    {
        gmx_fatal(FARGS,
                  "Wrong number of annealing values: %zu (for %d groups)\n",
                  simulatedAnnealingGroupNames.size(),
                  nr);
    }
    else
    {
        snew(ir->opts.annealing, nr);
        snew(ir->opts.anneal_npoints, nr);
        snew(ir->opts.anneal_time, nr);
        snew(ir->opts.anneal_temp, nr);
        for (i = 0; i < nr; i++)
        {
            ir->opts.annealing[i]      = SimulatedAnnealing::No;
            ir->opts.anneal_npoints[i] = 0;
            ir->opts.anneal_time[i]    = nullptr;
            ir->opts.anneal_temp[i]    = nullptr;
        }
        if (!simulatedAnnealingGroupNames.empty())
        {
            bAnneal = FALSE;
            for (i = 0; i < nr; i++)
            {
                if (gmx::equalCaseInsensitive(simulatedAnnealingGroupNames[i], "N", 1))
                {
                    ir->opts.annealing[i] = SimulatedAnnealing::No;
                }
                else if (gmx::equalCaseInsensitive(simulatedAnnealingGroupNames[i], "S", 1))
                {
                    ir->opts.annealing[i] = SimulatedAnnealing::Single;
                    bAnneal               = TRUE;
                }
                else if (gmx::equalCaseInsensitive(simulatedAnnealingGroupNames[i], "P", 1))
                {
                    ir->opts.annealing[i] = SimulatedAnnealing::Periodic;
                    bAnneal               = TRUE;
                }
            }
            if (bAnneal)
            {
                /* Read the other fields too */
                auto simulatedAnnealingPoints = gmx::splitString(inputrecStrings->anneal_npoints);
                if (simulatedAnnealingPoints.size() != simulatedAnnealingGroupNames.size())
                {
                    gmx_fatal(FARGS,
                              "Found %zu annealing-npoints values for %zu groups\n",
                              simulatedAnnealingPoints.size(),
                              simulatedAnnealingGroupNames.size());
                }
                convertInts(wi, simulatedAnnealingPoints, "annealing points", ir->opts.anneal_npoints);
                size_t numSimulatedAnnealingFields = 0;
                for (i = 0; i < nr; i++)
                {
                    if (ir->opts.anneal_npoints[i] == 1)
                    {
                        gmx_fatal(
                                FARGS,
                                "Please specify at least a start and an end point for annealing\n");
                    }
                    snew(ir->opts.anneal_time[i], ir->opts.anneal_npoints[i]);
                    snew(ir->opts.anneal_temp[i], ir->opts.anneal_npoints[i]);
                    numSimulatedAnnealingFields += ir->opts.anneal_npoints[i];
                }

                auto simulatedAnnealingTimes = gmx::splitString(inputrecStrings->anneal_time);

                if (simulatedAnnealingTimes.size() != numSimulatedAnnealingFields)
                {
                    gmx_fatal(FARGS,
                              "Found %zu annealing-time values, wanted %zu\n",
                              simulatedAnnealingTimes.size(),
                              numSimulatedAnnealingFields);
                }
                auto simulatedAnnealingTemperatures = gmx::splitString(inputrecStrings->anneal_temp);
                if (simulatedAnnealingTemperatures.size() != numSimulatedAnnealingFields)
                {
                    gmx_fatal(FARGS,
                              "Found %zu annealing-temp values, wanted %zu\n",
                              simulatedAnnealingTemperatures.size(),
                              numSimulatedAnnealingFields);
                }

                std::vector<real> allSimulatedAnnealingTimes(numSimulatedAnnealingFields);
                std::vector<real> allSimulatedAnnealingTemperatures(numSimulatedAnnealingFields);
                convertReals(wi, simulatedAnnealingTimes, "anneal-time", allSimulatedAnnealingTimes.data());
                convertReals(wi,
                             simulatedAnnealingTemperatures,
                             "anneal-temp",
                             allSimulatedAnnealingTemperatures.data());
                for (i = 0, k = 0; i < nr; i++)
                {
                    for (j = 0; j < ir->opts.anneal_npoints[i]; j++)
                    {
                        ir->opts.anneal_time[i][j] = allSimulatedAnnealingTimes[k];
                        ir->opts.anneal_temp[i][j] = allSimulatedAnnealingTemperatures[k];
                        if (j == 0)
                        {
                            if (ir->opts.anneal_time[i][0] > (ir->init_t + GMX_REAL_EPS))
                            {
                                gmx_fatal(FARGS, "First time point for annealing > init_t.\n");
                            }
                        }
                        else
                        {
                            /* j>0 */
                            if (ir->opts.anneal_time[i][j] < ir->opts.anneal_time[i][j - 1])
                            {
                                gmx_fatal(FARGS,
                                          "Annealing timepoints out of order: t=%f comes after "
                                          "t=%f\n",
                                          ir->opts.anneal_time[i][j],
                                          ir->opts.anneal_time[i][j - 1]);
                            }
                        }
                        if (ir->opts.anneal_temp[i][j] < 0)
                        {
                            gmx_fatal(FARGS,
                                      "Found negative temperature in annealing: %f\n",
                                      ir->opts.anneal_temp[i][j]);
                        }
                        k++;
                    }
                }
                /* Print out some summary information, to make sure we got it right */
                for (i = 0; i < nr; i++)
                {
                    if (ir->opts.annealing[i] != SimulatedAnnealing::No)
                    {
                        j = groups->groups[SimulationAtomGroupType::TemperatureCoupling][i];
                        fprintf(stderr,
                                "Simulated annealing for group %s: %s, %d timepoints\n",
                                *(groups->groupNames[j]),
                                enumValueToString(ir->opts.annealing[i]),
                                ir->opts.anneal_npoints[i]);
                        fprintf(stderr, "Time (ps)   Temperature (K)\n");
                        /* All terms except the last one */
                        for (j = 0; j < (ir->opts.anneal_npoints[i] - 1); j++)
                        {
                            fprintf(stderr,
                                    "%9.1f      %5.1f\n",
                                    ir->opts.anneal_time[i][j],
                                    ir->opts.anneal_temp[i][j]);
                        }

                        /* Finally the last one */
                        j = ir->opts.anneal_npoints[i] - 1;
                        if (ir->opts.annealing[i] == SimulatedAnnealing::Single)
                        {
                            fprintf(stderr,
                                    "%9.1f-     %5.1f\n",
                                    ir->opts.anneal_time[i][j],
                                    ir->opts.anneal_temp[i][j]);
                        }
                        else
                        {
                            fprintf(stderr,
                                    "%9.1f      %5.1f\n",
                                    ir->opts.anneal_time[i][j],
                                    ir->opts.anneal_temp[i][j]);
                            if (std::fabs(ir->opts.anneal_temp[i][j] - ir->opts.anneal_temp[i][0]) > GMX_REAL_EPS)
                            {
                                wi->addNote(
                                        "There is a temperature jump when your annealing "
                                        "loops back.\n");
                            }
                        }
                    }
                }
            }
        }
    }

    if (ir->bPull)
    {
        for (int i = 1; i < ir->pull->ngroup; i++)
        {
            const int gid = getGroupIndex(inputrecStrings->pullGroupNames[i], defaultIndexGroups);
            GMX_ASSERT(!defaultIndexGroups.empty(), "Must have initialized default index groups");
            atomGroupRangeValidation(natoms, defaultIndexGroups[gid].particleIndices);
        }

        process_pull_groups(ir->pull->group, inputrecStrings->pullGroupNames, defaultIndexGroups);

        checkPullCoords(ir->pull->group, ir->pull->coord);
    }

    if (ir->bRot)
    {
        make_rotation_groups(ir->rot.get(), inputrecStrings->rotateGroupNames, defaultIndexGroups);
    }

    if (ir->eSwapCoords != SwapType::No)
    {
        make_swap_groups(ir->swap, defaultIndexGroups);
    }

    /* Make indices for IMD session */
    if (ir->bIMD)
    {
        make_IMD_group(ir->imd, inputrecStrings->imd_grp, defaultIndexGroups);
    }

    gmx::IndexGroupsAndNames defaultIndexGroupsAndNames(defaultIndexGroups);
    mdModulesNotifiers.preProcessingNotifier_.notify(defaultIndexGroupsAndNames);

    auto accelerations          = gmx::splitString(inputrecStrings->acceleration);
    auto accelerationGroupNames = gmx::splitString(inputrecStrings->accelerationGroups);
    if (accelerationGroupNames.size() * DIM != accelerations.size())
    {
        gmx_fatal(FARGS,
                  "Invalid Acceleration input: %zu groups and %zu acc. values",
                  accelerationGroupNames.size(),
                  accelerations.size());
    }
    do_numbering(natoms,
                 groups,
                 accelerationGroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::Acceleration,
                 restnm,
                 GroupCoverage::AllGenerateRest,
                 bVerbose,
                 wi);
    nr = groups->groups[SimulationAtomGroupType::Acceleration].size();
    snew(ir->opts.acceleration, nr);
    ir->opts.ngacc = nr;

    convertRvecs(wi, accelerations, "accelerations", ir->opts.acceleration);

    auto freezeDims       = gmx::splitString(inputrecStrings->frdim);
    auto freezeGroupNames = gmx::splitString(inputrecStrings->freeze);
    if (freezeDims.size() != DIM * freezeGroupNames.size())
    {
        gmx_fatal(FARGS,
                  "Invalid Freezing input: %zu groups and %zu freeze values",
                  freezeGroupNames.size(),
                  freezeDims.size());
    }
    do_numbering(natoms,
                 groups,
                 freezeGroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::Freeze,
                 restnm,
                 GroupCoverage::AllGenerateRest,
                 bVerbose,
                 wi);
    nr             = groups->groups[SimulationAtomGroupType::Freeze].size();
    ir->opts.ngfrz = nr;
    snew(ir->opts.nFreeze, nr);
    for (i = k = 0; (size_t(i) < freezeGroupNames.size()); i++)
    {
        for (j = 0; (j < DIM); j++, k++)
        {
            ir->opts.nFreeze[i][j] = static_cast<int>(gmx::equalCaseInsensitive(freezeDims[k], "Y", 1));
            if (!ir->opts.nFreeze[i][j])
            {
                if (!gmx::equalCaseInsensitive(freezeDims[k], "N", 1))
                {
                    sprintf(warn_buf,
                            "Please use Y(ES) or N(O) for freezedim only "
                            "(not %s)",
                            freezeDims[k].c_str());
                    wi->addWarning(warn_buf);
                }
            }
        }
    }
    for (; (i < nr); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            ir->opts.nFreeze[i][j] = 0;
        }
    }

    auto energyGroupNames = gmx::splitString(inputrecStrings->energy);
    do_numbering(natoms,
                 groups,
                 energyGroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::EnergyOutput,
                 restnm,
                 GroupCoverage::AllGenerateRest,
                 bVerbose,
                 wi);
    add_wall_energrps(groups, ir->nwall, symtab);
    ir->opts.ngener    = groups->groups[SimulationAtomGroupType::EnergyOutput].size();
    auto vcmGroupNames = gmx::splitString(inputrecStrings->vcm);
    do_numbering(natoms,
                 groups,
                 vcmGroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::MassCenterVelocityRemoval,
                 restnm,
                 vcmGroupNames.empty() ? GroupCoverage::AllGenerateRest : GroupCoverage::Partial,
                 bVerbose,
                 wi);

    if (ir->comm_mode != ComRemovalAlgorithm::No)
    {
        checkAndUpdateVcmFreezeGroupConsistency(groups, natoms, ir->opts, wi);
    }

    /* Now we have filled the freeze struct, so we can calculate NRDF */
    calc_nrdf(mtop, ir, gnames);

    auto user1GroupNames = gmx::splitString(inputrecStrings->user1);
    do_numbering(natoms,
                 groups,
                 user1GroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::User1,
                 restnm,
                 GroupCoverage::AllGenerateRest,
                 bVerbose,
                 wi);
    auto user2GroupNames = gmx::splitString(inputrecStrings->user2);
    do_numbering(natoms,
                 groups,
                 user2GroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::User2,
                 restnm,
                 GroupCoverage::AllGenerateRest,
                 bVerbose,
                 wi);
    auto compressedXGroupNames = gmx::splitString(inputrecStrings->x_compressed_groups);
    do_numbering(natoms,
                 groups,
                 compressedXGroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::CompressedPositionOutput,
                 restnm,
                 GroupCoverage::OneGroup,
                 bVerbose,
                 wi);
    auto orirefFitGroupNames = gmx::splitString(inputrecStrings->orirefitgrp);
    do_numbering(natoms,
                 groups,
                 orirefFitGroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::OrientationRestraintsFit,
                 restnm,
                 GroupCoverage::AllGenerateRest,
                 bVerbose,
                 wi);

    /* MiMiC QMMM input processing */
    auto qmGroupNames = gmx::splitString(inputrecStrings->QMMM);
    if (qmGroupNames.size() > 1)
    {
        gmx_fatal(FARGS, "Currently, having more than one QM group in MiMiC is not supported");
    }
    /* group rest, if any, is always MM! */
    do_numbering(natoms,
                 groups,
                 qmGroupNames,
                 defaultIndexGroups,
                 SimulationAtomGroupType::QuantumMechanics,
                 restnm,
                 GroupCoverage::AllGenerateRest,
                 bVerbose,
                 wi);
    ir->opts.ngQM = qmGroupNames.size();

    /* end of MiMiC QMMM input */

    if (bVerbose)
    {
        for (auto group : gmx::keysOf(groups->groups))
        {
            fprintf(stderr, "%-16s has %zu element(s):", shortName(group), groups->groups[group].size());
            for (const auto& entry : groups->groups[group])
            {
                fprintf(stderr, " %s", *(groups->groupNames[entry]));
            }
            fprintf(stderr, "\n");
        }
    }

    nr = groups->groups[SimulationAtomGroupType::EnergyOutput].size();
    snew(ir->opts.egp_flags, nr * nr);

    bExcl = do_egp_flag(ir, groups, "energygrp-excl", inputrecStrings->egpexcl, EGP_EXCL);
    if (bExcl && ir->cutoff_scheme == CutoffScheme::Verlet)
    {
        wi->addError("Energy group exclusions are currently not supported");
    }
    if (bExcl && usingFullElectrostatics(ir->coulombtype))
    {
        wi->addWarning("Can not exclude the lattice Coulomb energy between energy groups");
    }

    bTable = do_egp_flag(ir, groups, "energygrp-table", inputrecStrings->egptable, EGP_TABLE);
    if (bTable && !(ir->vdwtype == VanDerWaalsType::User)
        && !(ir->coulombtype == CoulombInteractionType::User)
        && !(ir->coulombtype == CoulombInteractionType::PmeUser)
        && !(ir->coulombtype == CoulombInteractionType::PmeUserSwitch))
    {
        gmx_fatal(FARGS,
                  "Can only have energy group pair tables in combination with user tables for VdW "
                  "and/or Coulomb");
    }

    /* final check before going out of scope if simulated tempering variables
     * need to be set to default values.
     */
    if ((ir->expandedvals->nstexpanded < 0) && ir->bSimTemp)
    {
        ir->expandedvals->nstexpanded = 2 * static_cast<int>(ir->opts.tau_t[0] / ir->delta_t);
        wi->addWarning(gmx::formatString(
                "the value for nstexpanded was not specified for "
                " expanded ensemble simulated tempering. It is set to 2*tau_t (%d) "
                "by default, but it is recommended to set it to an explicit value!",
                ir->expandedvals->nstexpanded));
    }

    // Now that we have the temperature coupling options, we can process the ensemble temperature
    processEnsembleTemperature(ir, allAtomsAreTCoupled, wi);
}


static void check_disre(const gmx_mtop_t& mtop)
{
    if (gmx_mtop_ftype_count(mtop, F_DISRES) > 0)
    {
        const gmx_ffparams_t& ffparams  = mtop.ffparams;
        int                   ndouble   = 0;
        int                   old_label = -1;
        for (int i = 0; i < ffparams.numTypes(); i++)
        {
            int ftype = ffparams.functype[i];
            if (ftype == F_DISRES)
            {
                int label = ffparams.iparams[i].disres.label;
                if (label == old_label)
                {
                    fprintf(stderr, "Distance restraint index %d occurs twice\n", label);
                    ndouble++;
                }
                old_label = label;
            }
        }
        if (ndouble > 0)
        {
            gmx_fatal(FARGS,
                      "Found %d double distance restraint indices,\n"
                      "probably the parameters for multiple pairs in one restraint "
                      "are not identical\n",
                      ndouble);
        }
    }
}

//! Returns whether dimensions have an absolute reference due to walls, pbc or freezing
static BasicVector<bool> haveAbsoluteReference(const t_inputrec& ir)
{
    BasicVector<bool> absRef = { false, false, false };

    /* Check the degrees of freedom of the COM (not taking COMM removal into account) */
    for (int d = 0; d < DIM; d++)
    {
        absRef[d] = (d >= ndof_com(&ir));
    }
    /* Check for freeze groups */
    for (int g = 0; g < ir.opts.ngfrz; g++)
    {
        for (int d = 0; d < DIM; d++)
        {
            if (ir.opts.nFreeze[g][d] != 0)
            {
                absRef[d] = true;
            }
        }
    }

    return absRef;
}

//! Returns whether position restraints are used for dimensions
static BasicVector<bool> havePositionRestraints(const gmx_mtop_t& sys)
{
    BasicVector<bool> havePosres = { false, false, false };

    for (const auto ilists : IListRange(sys))
    {
        const auto& posResList   = ilists.list()[F_POSRES];
        const auto& fbPosResList = ilists.list()[F_FBPOSRES];
        if (ilists.nmol() > 0 && (!havePosres[XX] || !havePosres[YY] || !havePosres[ZZ]))
        {
            for (int i = 0; i < posResList.size(); i += 2)
            {
                const t_iparams& pr = sys.ffparams.iparams[posResList.iatoms[i]];
                for (int d = 0; d < DIM; d++)
                {
                    if (pr.posres.fcA[d] != 0)
                    {
                        havePosres[d] = true;
                    }
                }
            }
            for (int i = 0; i < fbPosResList.size(); i += 2)
            {
                /* Check for flat-bottom posres */
                const t_iparams& pr = sys.ffparams.iparams[fbPosResList.iatoms[i]];
                if (pr.fbposres.k != 0)
                {
                    switch (pr.fbposres.geom)
                    {
                        case efbposresSPHERE: havePosres = { true, true, true }; break;
                        case efbposresCYLINDERX: havePosres[YY] = havePosres[ZZ] = true; break;
                        case efbposresCYLINDERY: havePosres[XX] = havePosres[ZZ] = true; break;
                        case efbposresCYLINDER:
                        /* efbposres is a synonym for efbposresCYLINDERZ for backwards compatibility */
                        case efbposresCYLINDERZ: havePosres[XX] = havePosres[YY] = true; break;
                        case efbposresX: /* d=XX */
                        case efbposresY: /* d=YY */
                        case efbposresZ: /* d=ZZ */
                            havePosres[pr.fbposres.geom - efbposresX] = true;
                            break;
                        default:
                            gmx_fatal(FARGS,
                                      "Invalid geometry for flat-bottom position restraint.\n"
                                      "Expected nr between 1 and %d. Found %d\n",
                                      efbposresNR - 1,
                                      pr.fbposres.geom);
                    }
                }
            }
        }
    }

    return havePosres;
}

static void check_combination_rule_differences(const gmx_mtop_t& mtop,
                                               int               state,
                                               bool* bC6ParametersWorkWithGeometricRules,
                                               bool* bC6ParametersWorkWithLBRules,
                                               bool* bLBRulesPossible)
{
    int         ntypes, tpi, tpj;
    int*        typecount;
    real        tol;
    double      c6i, c6j, c12i, c12j;
    double      c6, c6_geometric, c6_LB;
    double      sigmai, sigmaj, epsi, epsj;
    bool        bCanDoLBRules, bCanDoGeometricRules;
    const char* ptr;

    /* A tolerance of 1e-5 seems reasonable for (possibly hand-typed)
     * force-field floating point parameters.
     */
    tol = 1e-5;
    ptr = getenv("GMX_LJCOMB_TOL");
    if (ptr != nullptr)
    {
        double            dbl;
        double gmx_unused canary;

        if (sscanf(ptr, "%lf%lf", &dbl, &canary) != 1)
        {
            gmx_fatal(
                    FARGS, "Could not parse a single floating-point number from GMX_LJCOMB_TOL (%s)", ptr);
        }
        tol = dbl;
    }

    *bC6ParametersWorkWithLBRules        = TRUE;
    *bC6ParametersWorkWithGeometricRules = TRUE;
    bCanDoLBRules                        = TRUE;
    ntypes                               = mtop.ffparams.atnr;
    snew(typecount, ntypes);
    gmx_mtop_count_atomtypes(mtop, state, typecount);
    *bLBRulesPossible = TRUE;
    for (tpi = 0; tpi < ntypes; ++tpi)
    {
        c6i  = mtop.ffparams.iparams[(ntypes + 1) * tpi].lj.c6;
        c12i = mtop.ffparams.iparams[(ntypes + 1) * tpi].lj.c12;
        for (tpj = tpi; tpj < ntypes; ++tpj)
        {
            c6j          = mtop.ffparams.iparams[(ntypes + 1) * tpj].lj.c6;
            c12j         = mtop.ffparams.iparams[(ntypes + 1) * tpj].lj.c12;
            c6           = mtop.ffparams.iparams[ntypes * tpi + tpj].lj.c6;
            c6_geometric = std::sqrt(c6i * c6j);
            if (!gmx_numzero(c6_geometric))
            {
                if (!gmx_numzero(c12i) && !gmx_numzero(c12j))
                {
                    sigmai = gmx::sixthroot(c12i / c6i);
                    sigmaj = gmx::sixthroot(c12j / c6j);
                    epsi   = c6i * c6i / (4.0 * c12i);
                    epsj   = c6j * c6j / (4.0 * c12j);
                    c6_LB  = 4.0 * std::sqrt(epsi * epsj) * gmx::power6(0.5 * (sigmai + sigmaj));
                }
                else
                {
                    *bLBRulesPossible = FALSE;
                    c6_LB             = c6_geometric;
                }
                bCanDoLBRules = gmx_within_tol(c6_LB, c6, tol);
            }

            if (!bCanDoLBRules)
            {
                *bC6ParametersWorkWithLBRules = FALSE;
            }

            bCanDoGeometricRules = gmx_within_tol(c6_geometric, c6, tol);

            if (!bCanDoGeometricRules)
            {
                *bC6ParametersWorkWithGeometricRules = FALSE;
            }
        }
    }
    sfree(typecount);
}

static void check_combination_rules(const t_inputrec* ir, const gmx_mtop_t& mtop, WarningHandler* wi)
{
    bool bLBRulesPossible, bC6ParametersWorkWithGeometricRules, bC6ParametersWorkWithLBRules;

    check_combination_rule_differences(
            mtop, 0, &bC6ParametersWorkWithGeometricRules, &bC6ParametersWorkWithLBRules, &bLBRulesPossible);
    if (ir->ljpme_combination_rule == LongRangeVdW::LB)
    {
        if (!bC6ParametersWorkWithLBRules || !bLBRulesPossible)
        {
            wi->addWarning(
                    "You are using arithmetic-geometric combination rules "
                    "in LJ-PME, but your non-bonded C6 parameters do not "
                    "follow these rules.");
        }
    }
    else
    {
        if (!bC6ParametersWorkWithGeometricRules)
        {
            if (ir->eDispCorr != DispersionCorrectionType::No)
            {
                wi->addNote(
                        "You are using geometric combination rules in "
                        "LJ-PME, but your non-bonded C6 parameters do "
                        "not follow these rules. "
                        "This will introduce very small errors in the forces and energies in "
                        "your simulations. Dispersion correction will correct total energy "
                        "and/or pressure for isotropic systems, but not forces or surface "
                        "tensions.");
            }
            else
            {
                wi->addNote(
                        "You are using geometric combination rules in "
                        "LJ-PME, but your non-bonded C6 parameters do "
                        "not follow these rules. "
                        "This will introduce very small errors in the forces and energies in "
                        "your simulations. If your system is homogeneous, consider using "
                        "dispersion correction "
                        "for the total energy and pressure.");
            }
        }
    }
}

static bool allTrue(const BasicVector<bool>& boolVector)
{
    return boolVector[0] && boolVector[1] && boolVector[2];
}

/* Generates an error or warning when lambda will become > 1, when appropriate */
static void checksForFepLambaLargerOne(const t_inputrec& ir, const gmx_mtop_t& mtop, WarningHandler* wi)
{
    /* The warnings below are only relevant if we have perturbed atoms */
    if (!haveFepPerturbedNBInteractions(mtop) && !haveFepPerturbedMasses(mtop))
    {
        return;
    }

    /* lambda vector components > 1 at the beginning of the simulation */
    if (ir.fepvals->delta_lambda < 0 && ir.fepvals->init_lambda_without_states > 1
        && ir.fepvals->n_lambda <= 0 && ir.fepvals->sc_alpha <= 0
        && ir.fepvals->softcoreFunction != SoftcoreType::Gapsys)
    {
        auto warningText = gmx::formatString(
                "You set init-lambda greater than 1 such that "
                "all lambdas will (initially) be greater than 1. "
                "Please only use this if you are aware of what you are "
                "doing.\n");
        wi->addWarning(warningText);
    }

    /* lambda vector components become > 1 in the course of the simulation */
    if (ir.fepvals->delta_lambda > 0 && ir.fepvals->init_lambda_without_states >= 0)
    {
        double stepNumberWhenLambdaIsOne =
                (1.0 - ir.fepvals->init_lambda_without_states) / ir.fepvals->delta_lambda;
        stepNumberWhenLambdaIsOne = std::max(stepNumberWhenLambdaIsOne, 0.0);
        int64_t intStepNumberWhenLambdaIsOne =
                static_cast<int64_t>(std::round(stepNumberWhenLambdaIsOne));

        if ((ir.nsteps < 0 || intStepNumberWhenLambdaIsOne < ir.nsteps) && ir.fepvals->n_lambda <= 0
            && ir.fepvals->sc_alpha <= 0 && ir.fepvals->softcoreFunction != SoftcoreType::Gapsys)
        {
            auto warningText = gmx::formatString(
                    "With init-lambda = %g and delta_lambda = %g and no lambda "
                    "vector given, "
                    "all lambdas will be greater than 1 after step %" PRId64 " of in total %" PRId64
                    " steps. "
                    "Please only use this if you are aware of what you are "
                    "doing.\n",
                    ir.fepvals->init_lambda_without_states,
                    ir.fepvals->delta_lambda,
                    intStepNumberWhenLambdaIsOne,
                    ir.nsteps);
            wi->addWarning(warningText);
        }
    }

    /* lambda vector components are set > 1 by the user */
    for (int j = 0; j < static_cast<int>(FreeEnergyPerturbationCouplingType::Count); j++)
    {
        auto enumValue = static_cast<FreeEnergyPerturbationCouplingType>(j);

        /* We have already issued an error for coul-lambdas and vdw-lambdas > 1
         * in combination with soft-core potentials.
         * As we only warn about lambdas > 1 if the non-bonded interactions or
         * masses are perturbed, let's only warn for coul, vdw, mass and bonded.
         */
        if ((((enumValue == FreeEnergyPerturbationCouplingType::Coul)
              || (enumValue == FreeEnergyPerturbationCouplingType::Vdw))
             && ((ir.fepvals->sc_alpha > 0) || (ir.fepvals->softcoreFunction == SoftcoreType::Gapsys)))
            || enumValue == FreeEnergyPerturbationCouplingType::Fep
            || enumValue == FreeEnergyPerturbationCouplingType::Restraint
            || enumValue == FreeEnergyPerturbationCouplingType::Temperature)
        {
            continue;
        }

        bool bElementsGreaterThanOne = false;

        for (int i = 0; i < ir.fepvals->n_lambda; i++)
        {
            if (ir.fepvals->all_lambda[j][i] > 1)
            {
                bElementsGreaterThanOne = true;
                break;
            }
        }

        if (bElementsGreaterThanOne)
        {
            auto warningText = gmx::formatString(
                    "One or more entries for %s are greater than 1. Please only use this "
                    "if you are aware of "
                    "what you are doing.",
                    enumValueToString(enumValue));
            wi->addWarning(warningText);
        }
    }
}

void triple_check(const char* mdparin, t_inputrec* ir, gmx_mtop_t* sys, WarningHandler* wi)
{
    // Not meeting MTS requirements should have resulted in a fatal error, so we can assert here
    GMX_ASSERT(gmx::checkMtsRequirements(*ir).empty(), "All MTS requirements should be met here");

    char                      err_buf[STRLEN];
    int                       i, m, c, nmol;
    bool                      bCharge;
    real *                    mgrp, mt;
    gmx_mtop_atomloop_block_t aloopb;
    char                      warn_buf[STRLEN];

    wi->setFileAndLineNumber(mdparin, -1);

    if (ir->comm_mode != ComRemovalAlgorithm::No && allTrue(havePositionRestraints(*sys)))
    {
        wi->addNote(
                "Removing center of mass motion in the presence of position restraints might "
                "cause artifacts. When you are using position restraints to equilibrate a "
                "macro-molecule, the artifacts are usually negligible.");
    }

    if (ir->cutoff_scheme == CutoffScheme::Verlet && ir->verletbuf_tol > 0 && ir->nstlist > 1
        && ((EI_MD(ir->eI) || EI_SD(ir->eI))
            && (ir->etc == TemperatureCoupling::VRescale || ir->etc == TemperatureCoupling::Berendsen)))
    {
        /* Check if a too small Verlet buffer might potentially
         * cause more drift than the thermostat can couple off.
         */
        /* Temperature error fraction for warning and suggestion */
        const real T_error_warn    = 0.002;
        const real T_error_suggest = 0.001;
        /* For safety: 2 DOF per atom (typical with constraints) */
        const real nrdf_at = 2;
        real       T, tau, max_T_error;
        int        i;

        T   = 0;
        tau = 0;
        for (i = 0; i < ir->opts.ngtc; i++)
        {
            T   = std::max(T, ir->opts.ref_t[i]);
            tau = std::max(tau, ir->opts.tau_t[i]);
        }
        if (T > 0)
        {
            /* This is a worst case estimate of the temperature error,
             * assuming perfect buffer estimation and no cancelation
             * of errors. The factor 0.5 is because energy distributes
             * equally over Ekin and Epot.
             */
            max_T_error = 0.5 * tau * ir->verletbuf_tol / (nrdf_at * gmx::c_boltz * T);
            if (max_T_error > T_error_warn)
            {
                sprintf(warn_buf,
                        "With a verlet-buffer-tolerance of %g kJ/mol/ps, a reference temperature "
                        "of %g and a tau_t of %g, your temperature might be off by up to %.1f%%. "
                        "To ensure the error is below %.1f%%, decrease verlet-buffer-tolerance to "
                        "%.0e or decrease tau_t.",
                        ir->verletbuf_tol,
                        T,
                        tau,
                        100 * max_T_error,
                        100 * T_error_suggest,
                        ir->verletbuf_tol * T_error_suggest / max_T_error);
                wi->addWarning(warn_buf);
            }
        }
    }

    if (ETC_ANDERSEN(ir->etc))
    {
        int i;

        for (i = 0; i < ir->opts.ngtc; i++)
        {
            sprintf(err_buf,
                    "all tau_t must currently be equal using Andersen temperature control, "
                    "violated for group %d",
                    i);
            CHECK(ir->opts.tau_t[0] != ir->opts.tau_t[i]);
            sprintf(err_buf,
                    "all tau_t must be positive using Andersen temperature control, "
                    "tau_t[%d]=%10.6f",
                    i,
                    ir->opts.tau_t[i]);
            CHECK(ir->opts.tau_t[i] < 0);
        }

        if (ir->etc == TemperatureCoupling::AndersenMassive && ir->comm_mode != ComRemovalAlgorithm::No)
        {
            for (i = 0; i < ir->opts.ngtc; i++)
            {
                int nsteps = gmx::roundToInt(ir->opts.tau_t[i] / ir->delta_t);
                sprintf(err_buf,
                        "tau_t/delta_t for group %d for temperature control method %s must be a "
                        "multiple of nstcomm (%d), as velocities of atoms in coupled groups are "
                        "randomized every time step. The input tau_t (%8.3f) leads to %d steps per "
                        "randomization",
                        i,
                        enumValueToString(ir->etc),
                        ir->nstcomm,
                        ir->opts.tau_t[i],
                        nsteps);
                CHECK(nsteps % ir->nstcomm != 0);
            }
        }
    }

    if (EI_DYNAMICS(ir->eI) && !EI_SD(ir->eI) && ir->eI != IntegrationAlgorithm::BD
        && ir->comm_mode == ComRemovalAlgorithm::No
        && !(allTrue(haveAbsoluteReference(*ir)) || allTrue(havePositionRestraints(*sys)) || ir->nsteps <= 10)
        && !ETC_ANDERSEN(ir->etc))
    {
        wi->addWarning(
                "You are not using center of mass motion removal (mdp option comm-mode), numerical "
                "rounding errors can lead to build up of kinetic energy of the center of mass");
    }

    if (ir->pressureCouplingOptions.epc == PressureCoupling::CRescale && !haveEnsembleTemperature(*ir))
    {
        sprintf(warn_buf,
                "Can not use the %s barostat without an ensemble temperature for the system",
                enumValueToString(ir->pressureCouplingOptions.epc));
        wi->addError(warn_buf);
    }

    if (ir->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
        && ir->etc == TemperatureCoupling::NoseHoover)
    {
        real tau_t_max = 0;
        for (int g = 0; g < ir->opts.ngtc; g++)
        {
            tau_t_max = std::max(tau_t_max, ir->opts.tau_t[g]);
        }
        if (ir->pressureCouplingOptions.tau_p < 1.9 * tau_t_max)
        {
            std::string message = gmx::formatString(
                    "With %s T-coupling and %s p-coupling, "
                    "%s (%g) should be at least twice as large as %s (%g) to avoid resonances",
                    enumValueToString(ir->etc),
                    enumValueToString(ir->pressureCouplingOptions.epc),
                    "tau-p",
                    ir->pressureCouplingOptions.tau_p,
                    "tau-t",
                    tau_t_max);
            wi->addWarning(message);
        }
    }

    /* Check for pressure coupling with absolute position restraints */
    if (ir->pressureCouplingOptions.epc != PressureCoupling::No
        && ir->pressureCouplingOptions.refcoord_scaling == RefCoordScaling::No)
    {
        const BasicVector<bool> havePosres = havePositionRestraints(*sys);
        {
            for (m = 0; m < DIM; m++)
            {
                if (havePosres[m] && norm2(ir->pressureCouplingOptions.compress[m]) > 0)
                {
                    wi->addWarning(
                            "You are using pressure coupling with absolute position restraints, "
                            "this will give artifacts. Use the refcoord_scaling option.");
                    break;
                }
            }
        }
    }

    if (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk && !haveConstantEnsembleTemperature(*ir))
    {
        sprintf(warn_buf,
                "The %s barostat requires a constant ensemble temperature for the system",
                enumValueToString(ir->pressureCouplingOptions.epc));
        wi->addError(warn_buf);
    }

    bCharge = FALSE;
    aloopb  = gmx_mtop_atomloop_block_init(*sys);
    const t_atom* atom;
    while (gmx_mtop_atomloop_block_next(aloopb, &atom, &nmol))
    {
        if (atom->q != 0 || atom->qB != 0)
        {
            bCharge = TRUE;
        }
    }

    if (!bCharge)
    {
        if (usingFullElectrostatics(ir->coulombtype))
        {
            sprintf(err_buf,
                    "You are using full electrostatics treatment %s for a system without charges.\n"
                    "This costs a lot of performance for just processing zeros, consider using %s "
                    "instead.\n",
                    enumValueToString(ir->coulombtype),
                    enumValueToString(CoulombInteractionType::Cut));
            wi->addWarning(err_buf);
        }
    }
    else
    {
        if (ir->coulombtype == CoulombInteractionType::Cut && ir->rcoulomb > 0)
        {
            sprintf(err_buf,
                    "You are using a plain Coulomb cut-off, which might produce artifacts.\n"
                    "You might want to consider using %s electrostatics.\n",
                    enumValueToString(CoulombInteractionType::Pme));
            wi->addNote(err_buf);
        }
    }

    /* Check if combination rules used in LJ-PME are the same as in the force field */
    if (usingLJPme(ir->vdwtype))
    {
        check_combination_rules(ir, *sys, wi);
    }

    /* Generalized reaction field */
    if (ir->coulombtype == CoulombInteractionType::GRFNotused)
    {
        wi->addError(
                "Generalized reaction-field electrostatics is no longer supported. "
                "You can use normal reaction-field instead and compute the reaction-field "
                "constant by hand.");
    }

    ir->useConstantAcceleration = false;
    for (int i = 0; (i < gmx::ssize(sys->groups.groups[SimulationAtomGroupType::Acceleration])); i++)
    {
        if (norm2(ir->opts.acceleration[i]) != 0)
        {
            ir->useConstantAcceleration = true;
        }
    }
    if (ir->useConstantAcceleration)
    {
        gmx::RVec acceleration = { 0.0_real, 0.0_real, 0.0_real };
        snew(mgrp, sys->groups.groups[SimulationAtomGroupType::Acceleration].size());
        for (const AtomProxy atomP : AtomRange(*sys))
        {
            const t_atom& local = atomP.atom();
            int           i     = atomP.globalAtomNumber();
            mgrp[getGroupType(sys->groups, SimulationAtomGroupType::Acceleration, i)] += local.m;
        }
        mt = 0.0;
        for (i = 0; (i < gmx::ssize(sys->groups.groups[SimulationAtomGroupType::Acceleration])); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                acceleration[m] += ir->opts.acceleration[i][m] * mgrp[i];
            }
            mt += mgrp[i];
        }
        for (m = 0; (m < DIM); m++)
        {
            if (std::fabs(acceleration[m]) > 1e-6)
            {
                const char* dim[DIM] = { "X", "Y", "Z" };
                fprintf(stderr,
                        "Net Acceleration in %s direction, will %s be corrected\n",
                        dim[m],
                        ir->nstcomm != 0 ? "" : "not");
                if (ir->nstcomm != 0 && m < ndof_com(ir))
                {
                    acceleration[m] /= mt;
                    for (i = 0;
                         (i < gmx::ssize(sys->groups.groups[SimulationAtomGroupType::Acceleration]));
                         i++)
                    {
                        ir->opts.acceleration[i][m] -= acceleration[m];
                    }
                }
            }
        }
        sfree(mgrp);
    }

    /* Checks related to FEP and slow growth */
    if (ir->efep != FreeEnergyPerturbationType::No)
    {
        if (ir->fepvals->sc_alpha != 0 && !gmx_within_tol(sys->ffparams.reppow, 12.0, 10 * GMX_DOUBLE_EPS))
        {
            gmx_fatal(FARGS,
                      "Soft-core interactions are only supported with VdW repulsion power 12");
        }

        checksForFepLambaLargerOne(*ir, *sys, wi);
    }

    if (ir->bPull)
    {
        bool bWarned;

        bWarned = FALSE;
        for (i = 0; i < ir->pull->ncoord && !bWarned; i++)
        {
            if (ir->pull->coord[i].eGeom != PullGroupGeometry::Transformation
                && (ir->pull->coord[i].group[0] == 0 || ir->pull->coord[i].group[1] == 0))
            {
                const auto absRef     = haveAbsoluteReference(*ir);
                const auto havePosres = havePositionRestraints(*sys);
                for (m = 0; m < DIM; m++)
                {
                    if (ir->pull->coord[i].dim[m] && !(absRef[m] || havePosres[m]))
                    {
                        wi->addWarning(
                                "You are using an absolute reference for pulling, but the rest of "
                                "the system does not have an absolute reference. This will lead to "
                                "artifacts.");
                        bWarned = TRUE;
                        break;
                    }
                }
            }
        }

        for (i = 0; i < 3; i++)
        {
            for (m = 0; m <= i; m++)
            {
                if ((ir->pressureCouplingOptions.epc != PressureCoupling::No
                     && ir->pressureCouplingOptions.compress[i][m] != 0)
                    || ir->deform[i][m] != 0)
                {
                    for (c = 0; c < ir->pull->ncoord; c++)
                    {
                        if (ir->pull->coord[c].eGeom == PullGroupGeometry::DirectionPBC
                            && ir->pull->coord[c].vec[m] != 0)
                        {
                            gmx_fatal(FARGS,
                                      "Can not have dynamic box while using pull geometry '%s' "
                                      "(dim %c)",
                                      enumValueToString(ir->pull->coord[c].eGeom),
                                      'x' + m);
                        }
                    }
                }
            }
        }
    }

    if (ir->bDoAwh && !haveConstantEnsembleTemperature(*ir))
    {
        wi->addError("With AWH a constant ensemble temperature is required");
    }

    if (ir_haveBoxDeformation(*ir))
    {
        if (EI_DYNAMICS(ir->eI) && ir->eI != IntegrationAlgorithm::MD
            && (EI_SD(ir->eI) || ir->etc != TemperatureCoupling::No))
        {
            sprintf(warn_buf,
                    "With all integrators except for %s, the whole velocity including the flow "
                    "driven by the deform option is scaled by the thermostat (note that the "
                    "reported kinetic energies and temperature are always computed excluding the "
                    "flow profile)",
                    enumValueToString(IntegrationAlgorithm::MD));
            wi->addNote(warn_buf);
        }

        if (ir->opts.ngtc != 1)
        {
            wi->addError("With box deformation, a single temperature coupling group is required");
        }
    }

    int numAccelerationAlgorithms = 0;
    if (ir->useConstantAcceleration)
    {
        numAccelerationAlgorithms++;
    }
    if (ir->cos_accel != 0)
    {
        numAccelerationAlgorithms++;
    }
    if (ir_haveBoxDeformation(*ir))
    {
        numAccelerationAlgorithms++;
    }
    if (numAccelerationAlgorithms > 1)
    {
        wi->addError(
                "Only one of the following three non-equilibrium methods is supported at a time: "
                "constant acceleration groups, cosine acceleration, box deformation");
    }

    check_disre(*sys);
}

void double_check(t_inputrec* ir, matrix box, bool bHasNormalConstraints, bool bHasAnyConstraints, WarningHandler* wi)
{
    char        warn_buf[STRLEN];
    const char* ptr;

    ptr = check_box(ir->pbcType, box);
    if (ptr)
    {
        wi->addError(ptr);
    }

    if (bHasNormalConstraints && ir->eConstrAlg == ConstraintAlgorithm::Shake)
    {
        if (ir->shake_tol <= 0.0)
        {
            sprintf(warn_buf, "ERROR: shake-tol must be > 0 instead of %g\n", ir->shake_tol);
            wi->addError(warn_buf);
        }
    }

    if ((ir->eConstrAlg == ConstraintAlgorithm::Lincs) && bHasNormalConstraints)
    {
        /* If we have Lincs constraints: */
        if (ir->eI == IntegrationAlgorithm::MD && ir->etc == TemperatureCoupling::No
            && ir->eConstrAlg == ConstraintAlgorithm::Lincs && ir->nLincsIter == 1)
        {
            sprintf(warn_buf,
                    "For energy conservation with LINCS, lincs_iter should be 2 or larger.\n");
            wi->addNote(warn_buf);
        }

        if ((ir->eI == IntegrationAlgorithm::CG || ir->eI == IntegrationAlgorithm::LBFGS)
            && (ir->nProjOrder < 8))
        {
            sprintf(warn_buf,
                    "For accurate %s with LINCS constraints, lincs-order should be 8 or more.",
                    enumValueToString(ir->eI));
            wi->addNote(warn_buf);
        }
        if (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            wi->addError("MTTK not compatible with lincs -- use shake instead.");
        }
    }

    if (bHasAnyConstraints && ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
    {
        wi->addError("Constraints are not implemented with MTTK pressure control.");
    }

    if (ir->LincsWarnAngle > 90.0)
    {
        sprintf(warn_buf, "lincs-warnangle can not be larger than 90 degrees, setting it to 90.\n");
        wi->addWarning(warn_buf);
        ir->LincsWarnAngle = 90.0;
    }

    if (ir->pbcType != PbcType::No)
    {
        if (ir->nstlist == 0)
        {
            wi->addWarning(
                    "With nstlist=0 atoms are only put into the box at step 0, therefore drifting "
                    "atoms might cause the simulation to crash.");
        }
        if (gmx::square(ir->rlist) >= max_cutoff2(ir->pbcType, box))
        {
            sprintf(warn_buf,
                    "ERROR: The cut-off length is longer than half the shortest box vector or "
                    "longer than the smallest box diagonal element. Increase the box size or "
                    "decrease rlist.\n");
            wi->addError(warn_buf);
        }
    }
}
