#!/usr/bin/env python3
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2013- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# This script is used by the GROMACS developers to generate the .cpp
# files that instiantiate all the template specializations of the kernel.
# This script also generates .h files with a list of external template
# declarations and a table of function points. This script is not called
# at CMake time, and users should never need to use it. The generated
# files are versions of the *.pre files in this directory, customized
# for the kernel structure type and/or the detailed kernel type. These
# are:
#
#   A single header file that declares all the kernel function template
#   specializations for this NBNxM kernel layout and a function pointer table.
#
#   Many C++ kernel files, each defining a single kernel function template
#   specialization. These functions can take a noticeable time to compile,
#   and should be in separate files to take advantage of make-time parallelism.
#
# These files are written to two separate directories for the two kernel
# layouts 2xMM and 4xM.
#
# Note that while functions for both NBNxM kernel layouts are compiled
# and built into an mdrun executable, because that executable
# is not portable, only the functions for the useful NBNxM kernel
# structure for the hardware selected at CMake time contain real
# kernel logic. A run-time error occurs if an inappropriate kernel
# dispatcher function is called (but that is normally impossible).

import re
import sys
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
import collections  # Requires Python 2.7

sys.path.append("../../../../admin")
from copyright import create_copyright_header

FileHeader = create_copyright_header("2012")
FileHeader += """/*
 * Note: this file was generated by the Verlet kernel generator for
 * kernel type {0}.
 */

"""


def read_kernel_template(filename):
    with open(filename, "r") as TemplateFile:
        TemplateText = TemplateFile.read()
    copyright_re = r"/\*\n \* This file is part of the GROMACS molecular simulation package\.\n( \*.*\n)* \*/\n"
    match = re.match(copyright_re, TemplateText)
    if match:
        TemplateText = TemplateText[match.end() :]
    return TemplateText


# The dict order must match the order of an enumeration in
# kernel_simd_template.c.pre
ElectrostaticsDict = collections.OrderedDict()
ElectrostaticsDict["ElecRF"] = {"param": "KernelCoulombType::RF, VdwCutoffCheck::No"}
ElectrostaticsDict["ElecQSTab"] = {
    "param": "KernelCoulombType::EwaldTabulated, VdwCutoffCheck::No"
}
ElectrostaticsDict["ElecQSTabTwinCut"] = {
    "param": "KernelCoulombType::EwaldTabulated, VdwCutoffCheck::Yes"
}
ElectrostaticsDict["ElecEw"] = {
    "param": "KernelCoulombType::EwaldAnalytical, VdwCutoffCheck::No"
}
ElectrostaticsDict["ElecEwTwinCut"] = {
    "param": "KernelCoulombType::EwaldAnalytical, VdwCutoffCheck::Yes"
}
ElectrostaticsDict["ElecNone"] = {
    "param": "KernelCoulombType::None, VdwCutoffCheck::Yes"
}


# The dict order must match the order of a C enumeration.
VdwTreatmentDict = collections.OrderedDict()
VdwTreatmentDict["VdwLJCombGeom"] = {
    "param": "LJCombinationRule::Geometric, InteractionModifiers::PotShift, LJEwald::None"
}
VdwTreatmentDict["VdwLJCombLB"] = {
    "param": "LJCombinationRule::LorentzBerthelot, InteractionModifiers::PotShift, LJEwald::None"
}
VdwTreatmentDict["VdwLJ"] = {
    "param": "LJCombinationRule::None, InteractionModifiers::PotShift, LJEwald::None"
}
VdwTreatmentDict["VdwLJFSw"] = {
    "param": "LJCombinationRule::None, InteractionModifiers::ForceSwitch, LJEwald::None"
}
VdwTreatmentDict["VdwLJPSw"] = {
    "param": "LJCombinationRule::None, InteractionModifiers::PotSwitch, LJEwald::None"
}
VdwTreatmentDict["VdwLJEwCombGeom"] = {
    "param": "LJCombinationRule::None, InteractionModifiers::PotShift, LJEwald::CombGeometric"
}

# This is OK as an unordered dict
EnergiesComputationDict = {
    "F": {
        "function type": "void",
        "param": "EnergyOutput::None",
    },
    "VF": {
        "function type": "void",
        "param": "EnergyOutput::System",
    },
    "VgrpF": {
        "function type": "void",
        "param": "EnergyOutput::GroupPairs",
    },
}

# This is OK as an unordered dict
VerletKernelTypeDict = {
    "2xmm": {
        "param": "KernelLayout::r2xMM",
        "define": "GMX_HAVE_NBNXM_SIMD_2XMM",
    },
    "4xm": {
        "param": "KernelLayout::r4xM",
        "define": "GMX_HAVE_NBNXM_SIMD_4XM",
    },
}

KernelsHeaderTemplate = read_kernel_template("kernel_simd_template.h.pre")

# For each Verlet kernel type, write two kinds of files:
#   a header file defining the functions for all the kernels and
#     the kernel function lookup table
#   for each kernel, a file defining the single C function for that kernel
for type in VerletKernelTypeDict:
    DirName = "../kernels_simd_{0}".format(type)
    KernelNamePrefix = "nbnxmKernelSimd"
    KernelFileNamePrefix = "nbnxm_kernel"
    KernelFileNamePrefix = "kernel"
    KernelsFileName = "{0}_simd_{1}".format(KernelFileNamePrefix, type)
    KernelsHeaderFileName = "kernels.h"
    KernelsHeaderPathName = "{1}".format(type, KernelsHeaderFileName)
    KernelFunctionLookupTable = {}
    KernelDeclarations = ""
    KernelTemplate = read_kernel_template("kernel_simd_template.cpp.pre")

    # Loop over all kernels
    for ener in EnergiesComputationDict:
        KernelFunctionLookupTable[ener] = ""
        for elec in ElectrostaticsDict:
            if elec == "ElecNone":
                KernelFunctionLookupTable[ener] += "#    if GMX_USE_EXT_FMM\n"
            KernelFunctionLookupTable[ener] += "    {\n"
            for ljtreat in VdwTreatmentDict:
                KernelName = "{0}<{1}, {2}, {3}, {4}>".format(
                    KernelNamePrefix,
                    VerletKernelTypeDict[type]["param"],
                    ElectrostaticsDict[elec]["param"],
                    VdwTreatmentDict[ljtreat]["param"],
                    EnergiesComputationDict[ener]["param"],
                )
                KernelFileName = "{0}_{1}_{2}_{3}".format(
                    KernelFileNamePrefix, elec, ljtreat, ener, type
                )

                # Declare the kernel function
                KernelDeclarations += "extern template {1}\n        {0}(\n".format(
                    KernelName, EnergiesComputationDict[ener]["function type"]
                )
                KernelDeclarations += "        const NbnxnPairlistCpu&    pairlist,\n        const nbnxn_atomdata_t&    nbat,\n        const interaction_const_t& ic,\n        const rvec*                shift_vec,\n        nbnxn_atomdata_output_t*   out);\n\n"

                # Write the file with the kernel definition
                with open(
                    "{0}/{1}.cpp".format(DirName, KernelFileName), "w"
                ) as kernelfp:
                    kernelfp.write(FileHeader.format(type))
                    kernelfp.write(
                        KernelTemplate.format(
                            VerletKernelTypeDict[type]["param"],
                            ElectrostaticsDict[elec]["param"],
                            VdwTreatmentDict[ljtreat]["param"],
                            EnergiesComputationDict[ener]["param"],
                            VerletKernelTypeDict[type]["define"],
                        )
                    )

                # Enter the kernel function in the lookup table
                KernelFunctionLookupTable[ener] += "            {0},\n".format(
                    KernelName
                )

            KernelFunctionLookupTable[ener] += "    },\n"
            if elec == "ElecNone":
                KernelFunctionLookupTable[ener] += "#    endif\n"

        KernelDeclarations += "\n"

    # Write the header file that declares all the kernel
    # functions for this type
    with open("{0}/{1}".format(DirName, KernelsHeaderFileName), "w") as fp:
        fp.write(FileHeader.format(type))
        fp.write(
            KernelsHeaderTemplate.format(
                KernelDeclarations,
                type,
                KernelFunctionLookupTable["F"],
                KernelFunctionLookupTable["VF"],
                KernelFunctionLookupTable["VgrpF"],
            )
        )

sys.exit()
