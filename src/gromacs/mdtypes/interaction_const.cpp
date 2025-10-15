/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#include "interaction_const.h"

#include "config.h"

#include <cmath>
#include <cstdio>

#include <filesystem>
#include <string>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"

interaction_const_t::SoftCoreParameters::SoftCoreParameters(const t_lambda& fepvals) :
    alphaVdw(fepvals.sc_alpha),
    alphaCoulomb(fepvals.bScCoul ? fepvals.sc_alpha : 0),
    lambdaPower(fepvals.sc_power),
    sigma6WithInvalidSigma(gmx::power6(fepvals.sc_sigma)),
    sigma6Minimum(fepvals.bScCoul ? gmx::power6(fepvals.sc_sigma_min) : 0),
    softcoreType(fepvals.softcoreFunction),
    gapsysScaleLinpointVdW(fepvals.scGapsysScaleLinpointLJ),
    gapsysScaleLinpointCoul(fepvals.scGapsysScaleLinpointQ),
    gapsysSigma6VdW(gmx::power6(fepvals.scGapsysSigmaLJ))
{
    // This is checked during tpr reading, so we can assert here
    GMX_RELEASE_ASSERT(fepvals.sc_r_power == 6.0, "We only support soft-core r-power 6");
}

/*! \brief Print Coulomb Ewald citations and set ewald coefficients */
static void initCoulombEwaldParameters(FILE*                                 fp,
                                       const t_inputrec&                     ir,
                                       bool                                  systemHasNetCharge,
                                       interaction_const_t::CoulombSettings* coulombSettings)
{
    if (!usingPmeOrEwald(ir.coulombtype))
    {
        return;
    }

    if (fp)
    {
        fprintf(fp, "Will do PME sum in reciprocal space for electrostatic interactions.\n");

        if (ir.coulombtype == CoulombInteractionType::P3mAD)
        {
            please_cite(fp, "Hockney1988");
            please_cite(fp, "Ballenegger2012");
        }
        else
        {
            please_cite(fp, "Essmann95a");
        }

        if (ir.ewald_geometry == EwaldGeometry::ThreeDC)
        {
            if (fp)
            {
                fprintf(fp,
                        "Using the Ewald3DC correction for systems with a slab geometry%s.\n",
                        systemHasNetCharge ? " and net charge" : "");
            }
            please_cite(fp, "In-Chul99a");
            if (systemHasNetCharge)
            {
                please_cite(fp, "Ballenegger2009");
            }
        }
    }

    coulombSettings->ewaldCoeff = calc_ewaldcoeff_q(ir.rcoulomb, ir.ewald_rtol);
    if (fp)
    {
        fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for Ewald\n", 1 / coulombSettings->ewaldCoeff);
    }

    if (coulombSettings->modifier == InteractionModifiers::PotShift)
    {
        GMX_RELEASE_ASSERT(coulombSettings->cutoff != 0, "Cutoff radius cannot be zero");
        coulombSettings->ewaldShift = std::erfc(coulombSettings->ewaldCoeff * coulombSettings->cutoff)
                                      / coulombSettings->cutoff;
    }
    else
    {
        coulombSettings->ewaldShift = 0;
    }
}

/*! \brief Print Van der Waals Ewald citations and set ewald coefficients */
static void initVdwEwaldParameters(FILE*                                     fp,
                                   const t_inputrec&                         ir,
                                   interaction_const_t::VanDerWaalsSettings* vdwSettings)
{
    if (!usingLJPme(ir.vdwtype))
    {
        return;
    }

    if (fp)
    {
        fprintf(fp, "Will do PME sum in reciprocal space for LJ dispersion interactions.\n");
        please_cite(fp, "Essmann95a");
    }
    vdwSettings->ewaldCoeff = calc_ewaldcoeff_lj(ir.rvdw, ir.ewald_rtol_lj);
    if (fp)
    {
        fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for LJ Ewald\n", 1 / vdwSettings->ewaldCoeff);
    }

    if (vdwSettings->modifier == InteractionModifiers::PotShift)
    {
        real crc2               = gmx::square(vdwSettings->ewaldCoeff * vdwSettings->cutoff);
        vdwSettings->ewaldShift = (std::exp(-crc2) * (1 + crc2 + 0.5 * crc2 * crc2) - 1)
                                  / gmx::power6(vdwSettings->cutoff);
    }
    else
    {
        vdwSettings->ewaldShift = 0;
    }
}

static void clear_force_switch_constants(shift_consts_t* sc)
{
    sc->c2   = 0;
    sc->c3   = 0;
    sc->cpot = 0;
}

static void force_switch_constants(real p, real rsw, real rc, shift_consts_t* sc)
{
    /* Here we determine the coefficient for shifting the force to zero
     * between distance rsw and the cut-off rc.
     * For a potential of r^-p, we have force p*r^-(p+1).
     * But to save flops we absorb p in the coefficient.
     * Thus we get:
     * force/p   = r^-(p+1) + c2*r^2 + c3*r^3
     * potential = r^-p + c2/3*r^3 + c3/4*r^4 + cpot
     */
    sc->c2   = ((p + 1) * rsw - (p + 4) * rc) / (std::pow(rc, p + 2) * gmx::square(rc - rsw));
    sc->c3   = -((p + 1) * rsw - (p + 3) * rc) / (std::pow(rc, p + 2) * gmx::power3(rc - rsw));
    sc->cpot = -std::pow(rc, -p) + p * sc->c2 / 3 * gmx::power3(rc - rsw)
               + p * sc->c3 / 4 * gmx::power4(rc - rsw);
}

static void potential_switch_constants(real rsw, real rc, switch_consts_t* sc)
{
    /* The switch function is 1 at rsw and 0 at rc.
     * The derivative and second derivate are zero at both ends.
     * rsw        = max(r - r_switch, 0)
     * sw         = 1 + c3*rsw^3 + c4*rsw^4 + c5*rsw^5
     * dsw        = 3*c3*rsw^2 + 4*c4*rsw^3 + 5*c5*rsw^4
     * force      = force*dsw - potential*sw
     * potential *= sw
     */
    sc->c3 = -10 / gmx::power3(rc - rsw);
    sc->c4 = 15 / gmx::power4(rc - rsw);
    sc->c5 = -6 / gmx::power5(rc - rsw);
}


interaction_const_t init_interaction_const(FILE*               fp,
                                           const t_inputrec&   ir,
                                           const gmx_mtop_t&   mtop,
                                           bool                systemHasNetCharge,
                                           std::optional<bool> anMDModuleProvidesDirectCoulomb)
{
    interaction_const_t interactionConst;

    interactionConst.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
    interactionConst.vdwEwaldTables     = std::make_unique<EwaldCorrectionTables>();

    /* Lennard-Jones */
    interaction_const_t::VanDerWaalsSettings& vdw = interactionConst.vdw;

    vdw.type               = ir.vdwtype;
    vdw.modifier           = ir.vdw_modifier;
    vdw.repulsionPower     = mtop.ffparams.reppow;
    vdw.cutoff             = cutoff_inf(ir.rvdw);
    vdw.switchDistance     = ir.rvdw_switch;
    vdw.pmeCombinationRule = ir.ljpme_combination_rule;

    initVdwEwaldParameters(fp, ir, &interactionConst.vdw);

    clear_force_switch_constants(&vdw.dispersionShift);
    clear_force_switch_constants(&vdw.repulsionShift);

    switch (vdw.modifier)
    {
        case InteractionModifiers::PotShift:
            /* Only shift the potential, don't touch the force */
            vdw.dispersionShift.cpot = -1.0 / gmx::power6(vdw.cutoff);
            vdw.repulsionShift.cpot  = -1.0 / gmx::power12(vdw.cutoff);
            break;
        case InteractionModifiers::ForceSwitch:
            /* Switch the force, switch and shift the potential */
            force_switch_constants(6.0, vdw.switchDistance, vdw.cutoff, &vdw.dispersionShift);
            force_switch_constants(12.0, vdw.switchDistance, vdw.cutoff, &vdw.repulsionShift);
            break;
        case InteractionModifiers::PotSwitch:
            /* Switch the potential and force */
            potential_switch_constants(vdw.switchDistance, vdw.cutoff, &vdw.switchConstants);
            break;
        case InteractionModifiers::None:
        case InteractionModifiers::ExactCutoff:
            /* Nothing to do here */
            break;
        default: gmx_incons("unimplemented potential modifier");
    }

    /* Electrostatics */
    interaction_const_t::CoulombSettings& coulomb = interactionConst.coulomb;

    coulomb.type           = ir.coulombtype;
    coulomb.modifier       = ir.coulomb_modifier;
    coulomb.cutoff         = cutoff_inf(ir.rcoulomb);
    coulomb.switchDistance = ir.rcoulomb_switch;
    coulomb.epsilon_r      = ir.epsilon_r;

    /* Set the Coulomb energy conversion factor */
    if (coulomb.epsilon_r != 0)
    {
        coulomb.epsfac = gmx::c_one4PiEps0 / coulomb.epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        coulomb.epsfac = 0;
    }

    /* Reaction-field */
    if (usingRF(coulomb.type))
    {
        GMX_RELEASE_ASSERT(coulomb.type != CoulombInteractionType::GRFNotused,
                           "GRF is no longer supported");
        calc_rffac(fp,
                   coulomb.epsilon_r,
                   ir.epsilon_rf,
                   coulomb.cutoff,
                   &coulomb.reactionFieldCoefficient,
                   &coulomb.reactionFieldShift);
    }
    else
    {
        /* For plain cut-off we might use the reaction-field kernels */
        coulomb.reactionFieldCoefficient = 0;
        if (ir.coulomb_modifier == InteractionModifiers::PotShift)
        {
            coulomb.reactionFieldShift = 1 / coulomb.cutoff;
        }
        else
        {
            coulomb.reactionFieldShift = 0;
        }
    }

    initCoulombEwaldParameters(fp, ir, systemHasNetCharge, &interactionConst.coulomb);

    if (fp != nullptr)
    {
        real dispersionShift;

        dispersionShift = vdw.dispersionShift.cpot;
        if (usingLJPme(vdw.type))
        {
            dispersionShift -= vdw.ewaldShift;
        }
        fprintf(fp, "Potential shift: LJ r^-12: %.3e r^-6: %.3e", vdw.repulsionShift.cpot, dispersionShift);

        if (coulomb.type == CoulombInteractionType::Cut)
        {
            fprintf(fp, ", Coulomb %.e", -coulomb.reactionFieldShift);
        }
        else if (usingPme(coulomb.type))
        {
            fprintf(fp, ", Ewald %.3e", -coulomb.ewaldShift);
        }
        fprintf(fp, "\n");
    }

    if (ir.efep != FreeEnergyPerturbationType::No)
    {
        GMX_RELEASE_ASSERT(ir.fepvals, "ir.fepvals should be set with free-energy");
        interactionConst.softCoreParameters =
                std::make_unique<interaction_const_t::SoftCoreParameters>(*ir.fepvals);
    }

    if (GMX_USE_EXT_FMM)
    {
        // Override the default only if the MDModule has specified it
        if (anMDModuleProvidesDirectCoulomb.has_value())
        {
            interactionConst.nbnxmIsDirectCoulombProvider = !anMDModuleProvidesDirectCoulomb.value();
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(!anMDModuleProvidesDirectCoulomb.has_value(),
                           "No MDModule may specify a direct Coulomb provider when an external FMM "
                           "module is not enabled");
    }

    return interactionConst;
}
