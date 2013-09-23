#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org

# This script is used by the GROMACS developers to build most of the
# files from which the nbnxn kernels are compiled. It is not called at
# CMake time, and users should never need to use it. It currently
# works for nbnxn kernel structure types 2xnn and 4xn. The generated
# files are versions of the *.pre files in this directory, customized
# for the kernel structure type and/or the detailed kernel type. These
# are:
#
#   A single header file that declares all the kernel functions for
#   this nbnxn kernel structure type, including the function that does
#   the dispatch via the function pointer table.
#
#   A single C kernel dispatcher file that defines the function that
#   decides at run time which kernel to call.
#
#   Many C kernel files, each defining a single kernel function. These
#   functions can take a noticeable time to compile, and should tend
#   to be in seperate files to take advantage of make-time
#   parallelism.
#
# This script should be run from the directory in which it is
# located. The generated files are located in ../simd_<type>. There
# are three other files in those locations that are not generated. These
# contain:
#
#   setup logic peculiar to the kernel structure type but common to
#   all the kernels within that type, and
#
#   the logic for the outer and inner loops of the kernels, as
#   customized by numerous preprocessor defines to suit the hardware
#   and kernel type.
#
# Note that while functions for both nbnxn kernel structures are
# compiled and built into an mdrun executable, because that executable
# is not portable, only the functions for the useful nbnxn kernel
# structure for the hardware selected at CMake time contain real
# kernel logic. A run-time error occurs if an inappropriate kernel
# dispatcher function is called (but that is normally impossible).

import sys
import os
import collections # Requires Python 2.7

FileHeader = \
'/*\n' \
' * This file is part of the GROMACS molecular simulation package.\n' \
' *\n' \
' * Copyright (c) 2012,2013, by the GROMACS development team, led by\n' \
' * David van der Spoel, Berk Hess, Erik Lindahl, and including many\n' \
' * others, as listed in the AUTHORS file in the top-level source\n' \
' * directory and at http://www.gromacs.org.\n' \
' *\n' \
' * GROMACS is free software; you can redistribute it and/or\n' \
' * modify it under the terms of the GNU Lesser General Public License\n' \
' * as published by the Free Software Foundation; either version 2.1\n' \
' * of the License, or (at your option) any later version.\n' \
' *\n' \
' * GROMACS is distributed in the hope that it will be useful,\n' \
' * but WITHOUT ANY WARRANTY; without even the implied warranty of\n' \
' * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n' \
' * Lesser General Public License for more details.\n' \
' *\n' \
' * You should have received a copy of the GNU Lesser General Public\n' \
' * License along with GROMACS; if not, see\n' \
' * http://www.gnu.org/licenses, or write to the Free Software Foundation,\n' \
' * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.\n' \
' *\n' \
' * If you want to redistribute modifications to GROMACS, please\n' \
' * consider that scientific software is very special. Version\n' \
' * control is crucial - bugs must be traceable. We will be happy to\n' \
' * consider code for inclusion in the official distribution, but\n' \
' * derived work must not be called official GROMACS. Details are found\n' \
' * in the README & COPYING files - if they are missing, get the\n' \
' * official version at http://www.gromacs.org.\n' \
' *\n' \
' * To help us fund GROMACS development, we humbly ask that you cite\n' \
' * the research papers on the package. Check out http://www.gromacs.org.\n' \
' */\n' \
'/*\n' \
' * Note: this file was generated by the Verlet kernel generator for\n' \
' * kernel type {0}.\n' \
' */\n\n'

# The dict order must match the order of an enumeration in
# nbnxn_kernel_simd_template.c.pre
ElectrostaticsDict = collections.OrderedDict()
ElectrostaticsDict['rf'] = { 'define' : '#define CALC_COUL_RF' }
ElectrostaticsDict['tab'] = { 'define' : '#define CALC_COUL_TAB' }
ElectrostaticsDict['tab_twin'] = { 'define' : '#define CALC_COUL_TAB\n#define VDW_CUTOFF_CHECK /* Use twin-range cut-off */' }
ElectrostaticsDict['ewald'] = { 'define' : '#define CALC_COUL_EWALD' }
ElectrostaticsDict['ewald_twin'] = { 'define' : '#define CALC_COUL_EWALD\n#define VDW_CUTOFF_CHECK /* Use twin-range cut-off */' }

# The dict order must match the order of a C enumeration.
LJCombinationRuleDict = collections.OrderedDict()
LJCombinationRuleDict['geom'] = { 'define' : '#define LJ_COMB_GEOM' }
LJCombinationRuleDict['lb'] = { 'define' : '#define LJ_COMB_LB' }
LJCombinationRuleDict['none'] = { 'define' : '/* Use no LJ combination rule */' }

# This is OK as an unordered dict
EnergiesComputationDict = {
    'ener'    : {
        'function type' : 'nbk_func_ener',
        'define' : '#define CALC_ENERGIES',
    },
    'energrp' : {
        'function type' : 'nbk_func_ener',
        'define' : '#define CALC_ENERGIES\n#define ENERGY_GROUPS',
    },
    'noener'  : {
        'function type' : 'nbk_func_noener',
        'define' : '/* Will not calculate energies */',
    },
}

# This is OK as an unordered dict
VerletKernelTypeDict = {
    '2xnn' : {
        'Define' : 'GMX_NBNXN_SIMD_2XNN',
        'WidthSetup' : '/* Include the full-width SIMD macros */\n',
        'WidthCheck' : ('#if !(GMX_SIMD_WIDTH_HERE == 8 || GMX_SIMD_WIDTH_HERE == 16)\n' \
                        '#error "unsupported SIMD width"\n' \
                        '#endif\n'),
        'UnrollSize' : 2,
    },
    '4xn' : {
        'Define' : 'GMX_NBNXN_SIMD_4XN',
        'WidthSetup' : ('#ifdef GMX_NBNXN_HALF_WIDTH_SIMD\n' \
                        '#define GMX_USE_HALF_WIDTH_SIMD_HERE\n' \
                        '#endif\n'),
        'WidthCheck' : ('#if !(GMX_SIMD_WIDTH_HERE == 2 || GMX_SIMD_WIDTH_HERE == 4 || GMX_SIMD_WIDTH_HERE == 8)\n' \
                        '#error "unsupported SIMD width"\n' \
                        '#endif\n'),
        'UnrollSize' : 1,
    },
}

with open ("nbnxn_kernel_simd_template.c.pre", "r") as KernelDispatcherTemplateFile:
    KernelDispatcherTemplate = KernelDispatcherTemplateFile.read()

with open ("nbnxn_kernel_simd_template.h.pre", "r") as KernelsHeaderTemplateFile:
    KernelsHeaderTemplate =  KernelsHeaderTemplateFile.read()

# For each Verlet kernel type, write three kinds of files:
#   a header file defining the functions for all the kernels,
#   a code file containing the kernel function lookup table and
#     the kernel dispatcher function
#   for each kernel, a file defining the single C function for that kernel
for type in VerletKernelTypeDict:
    DirName = "../simd_{0}".format(type)
    KernelNamePrefix = 'nbnxn_kernel_simd_{0}'.format(type)
    KernelsHeaderFileName = "{0}.h".format(KernelNamePrefix)
    KernelFunctionLookupTable = {}
    KernelDeclarations = ''
    with open ("{1}_kernel.c.pre".format(DirName,KernelNamePrefix), "r") as KernelTemplateFile:
        KernelTemplate =  KernelTemplateFile.read()

    # Loop over all kernels
    for ener in EnergiesComputationDict:
        KernelFunctionLookupTable[ener] = '{\n'
        for elec in ElectrostaticsDict:
            KernelFunctionLookupTable[ener] += '    {\n'
            for ljcomb in LJCombinationRuleDict:
                KernelName = ('{0}_{1}_comb_{2}_{3}'
                              .format(KernelNamePrefix,elec,ljcomb,ener))

                # Declare the kernel function
                KernelDeclarations += ('{1:21} {0};\n'
                                       .format(KernelName,
                                               EnergiesComputationDict[ener]['function type']))

                # Write the file with the kernel definition
                with open('{0}/{1}.c'.format(DirName,KernelName), 'w') as kernelfp:
                    kernelfp.write(FileHeader.format(type))
                    kernelfp.write(KernelTemplate
                                   .format(VerletKernelTypeDict[type]['Define'],
                                           ElectrostaticsDict[elec]['define'],
                                           LJCombinationRuleDict[ljcomb]['define'],
                                           EnergiesComputationDict[ener]['define'],
                                           KernelsHeaderFileName,
                                           KernelName,
                                           " " * (len(KernelName) + 1),
                                           VerletKernelTypeDict[type]['UnrollSize'],
                                       )
                               )

                # Enter the kernel function in the lookup table
                KernelFunctionLookupTable[ener] += '        {0},\n'.format(KernelName)

            KernelFunctionLookupTable[ener] += '    },\n'
        KernelFunctionLookupTable[ener] += '};\n'
        KernelDeclarations += '\n'

    # Write the header file that declares all the kernel
    # functions for this type
    with open('{0}/{1}'.format(DirName,KernelsHeaderFileName),'w') as fp:
        fp.write(FileHeader.format(type))
        fp.write(KernelsHeaderTemplate
                 .format(KernelNamePrefix,
                         " " * (len(KernelNamePrefix) + 1),
                         KernelDeclarations))

    # Write the file defining the kernel dispatcher
    # function for this type
    with open('{0}/{1}'.format(DirName,"{0}.c".format(KernelNamePrefix)),'w') as fp:
        fp.write(FileHeader.format(type))
        fp.write(KernelDispatcherTemplate
                 .format(VerletKernelTypeDict[type]['Define'],
                         VerletKernelTypeDict[type]['WidthSetup'],
                         VerletKernelTypeDict[type]['WidthCheck'],
                         VerletKernelTypeDict[type]['UnrollSize'],
                         KernelsHeaderFileName,
                         KernelNamePrefix,
                         ' ' * (len(KernelNamePrefix)+1),
                         KernelFunctionLookupTable['ener'],
                         KernelFunctionLookupTable['energrp'],
                         KernelFunctionLookupTable['noener'],
                     )
             )

sys.exit()
