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

#include "md_enums.h"

#include "gromacs/utility/enumerationhelpers.h"

const char* enum_name(int index, int max_index, const char* const names[])
{
    if (index < 0 || index >= max_index)
    {
        static const char* undef = "no name defined";
        return undef;
    }
    else
    {
        return names[index];
    }
}

const char* enumValueToString(IntegrationAlgorithm enumValue)
{
    static constexpr gmx::EnumerationArray<IntegrationAlgorithm, const char*> interationAlgorithmNames = {
        "md",  "steep", "cg", "bd",    "sd2 - removed", "nm",   "l-bfgs",
        "tpi", "tpic",  "sd", "md-vv", "md-vv-avek",    "mimic"
    };
    return interationAlgorithmNames[enumValue];
}

const char* enumValueToString(CoulombInteractionType enumValue)
{
    static constexpr gmx::EnumerationArray<CoulombInteractionType, const char*> coloumbTreatmentNames = {
        "Cut-off",
        "Reaction-Field",
        "Generalized-Reaction-Field (unused)",
        "PME",
        "Ewald",
        "P3M-AD",
        "Poisson",
        "Switch",
        "Shift",
        "User",
        "Generalized-Born (unused)",
        "Reaction-Field-nec (unsupported)",
        "Encad-shift (unused)",
        "PME-User",
        "PME-Switch",
        "PME-User-Switch",
        "Reaction-Field-zero"
    };
    return coloumbTreatmentNames[enumValue];
}

const char* enumValueToString(EwaldGeometry enumValue)
{
    static constexpr gmx::EnumerationArray<EwaldGeometry, const char*> ewaldGeometryNames = {
        "3d", "3dc"
    };
    return ewaldGeometryNames[enumValue];
}

const char* enumValueToString(LongRangeVdW enumValue)
{
    static constexpr gmx::EnumerationArray<LongRangeVdW, const char*> longRangeVdWNames = {
        "Geometric", "Lorentz-Berthelot"
    };
    return longRangeVdWNames[enumValue];
}

const char* enumValueToString(VanDerWaalsType enumValue)
{
    static constexpr gmx::EnumerationArray<VanDerWaalsType, const char*> vanDerWaalsTypeNames = {
        "Cut-off", "Switch", "Shift", "User", "Encad-shift (unused)", "PME"
    };
    return vanDerWaalsTypeNames[enumValue];
}

const char* enumValueToString(ConstraintAlgorithm enumValue)
{
    static constexpr gmx::EnumerationArray<ConstraintAlgorithm, const char*> constraintAlgorithmNames = {
        "Lincs", "Shake"
    };
    return constraintAlgorithmNames[enumValue];
}

const char* enumValueToString(InteractionModifiers enumValue)
{
    static constexpr gmx::EnumerationArray<InteractionModifiers, const char*> interactionModifierNames = {
        "Potential-shift-Verlet", "Potential-shift", "None",
        "Potential-switch",       "Exact-cutoff",    "Force-switch"
    };
    return interactionModifierNames[enumValue];
}

const char* enumValueToString(TemperatureCoupling enumValue)
{
    static constexpr gmx::EnumerationArray<TemperatureCoupling, const char*> temperatureCouplingNames = {
        "No", "Berendsen", "Nose-Hoover", "yes", "Andersen", "Andersen-massive", "V-rescale"
    }; /* yes is alias for berendsen */
    return temperatureCouplingNames[enumValue];
}

const char* enumValueToString(PressureCoupling enumValue)
{
    static constexpr gmx::EnumerationArray<PressureCoupling, const char*> pressureCouplingNames = {
        "No", "Berendsen", "Parrinello-Rahman", "Isotropic", "MTTK", "C-rescale"
    }; /* isotropic is alias for berendsen */
    return pressureCouplingNames[enumValue];
}

const char* enumValueToString(Boolean enumValue)
{
    static constexpr gmx::EnumerationArray<Boolean, const char*> booleanNames = { "no", "yes" };
    return booleanNames[enumValue];
}

const char* booleanValueToString(bool value)
{
    Boolean enumValue = value ? Boolean::Yes : Boolean::No;
    return enumValueToString(enumValue);
}

const char* enumValueToString(RefCoordScaling enumValue)
{
    static constexpr gmx::EnumerationArray<RefCoordScaling, const char*> refCoordScalingNames = {
        "No", "All", "COM"
    };
    return refCoordScalingNames[enumValue];
}

const char* enumValueToString(CutoffScheme enumValue)
{
    static constexpr gmx::EnumerationArray<CutoffScheme, const char*> cutoffSchemeNames = {
        "Verlet", "Group"
    };
    return cutoffSchemeNames[enumValue];
}

const char* enumValueToString(PressureCouplingType enumValue)
{
    static constexpr gmx::EnumerationArray<PressureCouplingType, const char*> pressureCouplingTypeNames = {
        "Isotropic", "Semiisotropic", "Anisotropic", "Surface-Tension"
    };
    return pressureCouplingTypeNames[enumValue];
}

const char* enumValueToString(DistanceRestraintRefinement enumValue)
{
    static constexpr gmx::EnumerationArray<DistanceRestraintRefinement, const char*> distanceRestraintRefinementNames = {
        "No", "Simple", "Ensemble"
    };
    return distanceRestraintRefinementNames[enumValue];
}

const char* enumValueToString(DistanceRestraintWeighting enumValue)
{
    static constexpr gmx::EnumerationArray<DistanceRestraintWeighting, const char*> distanceRestraintWeightingNames = {
        "Conservative", "Equal"
    };
    return distanceRestraintWeightingNames[enumValue];
}

const char* enumValueToString(VanDerWaalsPotential enumValue)
{
    static constexpr gmx::EnumerationArray<VanDerWaalsPotential, const char*> vanDerWaalsPotentialNames = {
        "None", "LJ", "Buckingham"
    };
    return vanDerWaalsPotentialNames[enumValue];
}

const char* enumValueToString(CombinationRule enumValue)
{
    static constexpr gmx::EnumerationArray<CombinationRule, const char*> combinationRuleNames = {
        "None", "Geometric", "Arithmetic", "GeomSigEps"
    };
    return combinationRuleNames[enumValue];
}

const char* enumValueToString(SimulatedTempering enumValue)
{
    static constexpr gmx::EnumerationArray<SimulatedTempering, const char*> simulatedTemperingNames = {
        "geometric", "exponential", "linear"
    };
    return simulatedTemperingNames[enumValue];
}

const char* enumValueToString(FreeEnergyPerturbationType enumValue)
{
    static constexpr gmx::EnumerationArray<FreeEnergyPerturbationType, const char*> freeEnergyPerturbationTypeNames = {
        "no", "yes", "static", "slow-growth", "expanded"
    };
    return freeEnergyPerturbationTypeNames[enumValue];
}

const char* enumValueToString(FreeEnergyPerturbationCouplingType enumValue)
{
    static constexpr gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, const char*> freeEnergyPerturbationCouplingTypeNames = {
        "fep-lambdas",    "mass-lambdas",      "coul-lambdas",       "vdw-lambdas",
        "bonded-lambdas", "restraint-lambdas", "temperature-lambdas"
    };
    return freeEnergyPerturbationCouplingTypeNames[enumValue];
}

const char* enumValueToStringSingular(FreeEnergyPerturbationCouplingType enumValue)
{
    static constexpr gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, const char*> freeEnergyPerturbationCouplingTypeNames = {
        "fep-lambda",    "mass-lambda",      "coul-lambda",       "vdw-lambda",
        "bonded-lambda", "restraint-lambda", "temperature-lambda"
    };
    return freeEnergyPerturbationCouplingTypeNames[enumValue];
}

const char* enumValueToString(FreeEnergyPrintEnergy enumValue)
{
    static constexpr gmx::EnumerationArray<FreeEnergyPrintEnergy, const char*> freeEnergyPrintNames = {
        "no", "total", "potential", "yes"
    };
    return freeEnergyPrintNames[enumValue];
}

const char* enumValueToString(LambdaWeightCalculation enumValue)
{
    static constexpr gmx::EnumerationArray<LambdaWeightCalculation, const char*> lambdaWeightCalculationNames = {
        "no",     "metropolis-transition", "barker-transition",
        "minvar", "wang-landau",           "weighted-wang-landau"
    };
    return lambdaWeightCalculationNames[enumValue];
}

const char* enumValueToString(LambdaMoveCalculation enumValue)
{
    static constexpr gmx::EnumerationArray<LambdaMoveCalculation, const char*> lambdaMoveCalculationNames = {
        "no", "metropolis", "barker", "gibbs", "metropolized-gibbs"
    };
    return lambdaMoveCalculationNames[enumValue];
}

const char* enumValueToString(LambdaWeightWillReachEquilibrium enumValue)
{
    static constexpr gmx::EnumerationArray<LambdaWeightWillReachEquilibrium, const char*> lambdaWeightEquilibriumNames = {
        "no",         "yes", "wl-delta", "number-all-lambda", "number-steps", "number-samples",
        "count-ratio"
    };
    return lambdaWeightEquilibriumNames[enumValue];
}

const char* enumValueToString(SeparateDhdlFile enumValue)
{
    static constexpr gmx::EnumerationArray<SeparateDhdlFile, const char*> separateDhdlFileNames = {
        "yes", "no"
    };
    return separateDhdlFileNames[enumValue];
}

const char* enumValueToString(SoftcoreType enumValue)
{
    static constexpr gmx::EnumerationArray<SoftcoreType, const char*> softcoreTypeNames = {
        "beutler", "gapsys"
    };
    return softcoreTypeNames[enumValue];
}

const char* enumValueToString(KernelSoftcoreType enumValue)
{
    static constexpr gmx::EnumerationArray<KernelSoftcoreType, const char*> softcoreTypeNames = {
        "beutler", "gapsys", "none"
    };
    return softcoreTypeNames[enumValue];
}

const char* enumValueToString(DhDlDerivativeCalculation enumValue)
{
    static constexpr gmx::EnumerationArray<DhDlDerivativeCalculation, const char*> dhdlDerivativeCalculationNames = {
        "yes", "no"
    };
    return dhdlDerivativeCalculationNames[enumValue];
}

const char* enumValueToString(SolventModel enumValue)
{
    static constexpr gmx::EnumerationArray<SolventModel, const char*> solventModelNames = {
        "No", "SPC", "TIP4p"
    };
    return solventModelNames[enumValue];
}

const char* enumValueToString(DispersionCorrectionType enumValue)
{
    static constexpr gmx::EnumerationArray<DispersionCorrectionType, const char*> dispersionCorrectionTypeNames = {
        "No", "EnerPres", "Ener", "AllEnerPres", "AllEner"
    };
    return dispersionCorrectionTypeNames[enumValue];
}

const char* enumValueToString(SimulatedAnnealing enumValue)
{
    static constexpr gmx::EnumerationArray<SimulatedAnnealing, const char*> simulatedAnnealingNames = {
        "No", "Single", "Periodic"
    };
    return simulatedAnnealingNames[enumValue];
}

const char* enumValueToString(WallType enumValue)
{
    static constexpr gmx::EnumerationArray<WallType, const char*> wallTypeNames = {
        "9-3", "10-4", "table", "12-6"
    };
    return wallTypeNames[enumValue];
}

const char* enumValueToString(PullingAlgorithm enumValue)
{
    static constexpr gmx::EnumerationArray<PullingAlgorithm, const char*> pullAlgorithmNames = {
        "umbrella",    "constraint",       "constant-force",
        "flat-bottom", "flat-bottom-high", "external-potential"
    };
    return pullAlgorithmNames[enumValue];
}

const char* enumValueToString(PullGroupGeometry enumValue)
{
    static constexpr gmx::EnumerationArray<PullGroupGeometry, const char*> pullGroupControlNames = {
        "distance", "direction", "cylinder",   "direction-periodic", "direction-relative",
        "angle",    "dihedral",  "angle-axis", "transformation"
    };
    return pullGroupControlNames[enumValue];
}

const char* enumValueToString(EnforcedRotationGroupType enumValue)
{
    static constexpr gmx::EnumerationArray<EnforcedRotationGroupType, const char*> enforcedRotationGroupNames = {
        "iso", "iso-pf", "pm",   "pm-pf",  "rm",    "rm-pf",
        "rm2", "rm2-pf", "flex", "flex-t", "flex2", "flex2-t"
    };
    return enforcedRotationGroupNames[enumValue];
}

const char* enumValueToString(RotationGroupFitting enumValue)
{
    static constexpr gmx::EnumerationArray<RotationGroupFitting, const char*> rotationGroupFittingNames = {
        "rmsd", "norm", "potential"
    };
    return rotationGroupFittingNames[enumValue];
}

const char* enumValueToString(SwapType enumValue)
{
    static constexpr gmx::EnumerationArray<SwapType, const char*> swapTypeNames = {
        "no", "X", "Y", "Z"
    };
    return swapTypeNames[enumValue];
}

const char* enumValueToString(SwapGroupSplittingType enumValue)
{
    static constexpr gmx::EnumerationArray<SwapGroupSplittingType, const char*> swapGroupSplittingTypeNames = {
        "Split0", "Split1", "Solvent"
    };
    return swapGroupSplittingTypeNames[enumValue];
}

const char* enumValueToString(NbkernelElecType enumValue)
{
    static constexpr gmx::EnumerationArray<NbkernelElecType, const char*> nbkernelElecTypeNames = {
        "None", "Coulomb", "Reaction-Field", "Cubic-Spline-Table", "Ewald"
    };
    return nbkernelElecTypeNames[enumValue];
}

const char* enumValueToString(NbkernelVdwType enumValue)
{
    static constexpr gmx::EnumerationArray<NbkernelVdwType, const char*> nbkernelVdwTypeNames = {
        "None", "Lennard-Jones", "Buckingham", "Cubic-Spline-Table", "LJEwald"
    };
    return nbkernelVdwTypeNames[enumValue];
}
