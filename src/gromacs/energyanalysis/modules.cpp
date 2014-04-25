/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gromacs/energyanalysis/modules.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"

#include "gromacs/energyanalysis/viscosity.h"
#include "gromacs/energyanalysis/fluctprops.h"
#include "gromacs/energyanalysis/dhdl.h"
#include "gromacs/energyanalysis/freeenergydifference.h"
#include "gromacs/energyanalysis/simple.h"
#include "gromacs/energyanalysis/handler.h"

static int gmx_etool(int                   *argc,
                     char                  *argv[],
                     gmx::EnergyAnalysisPtr ptr)
{
    // Energy handler utility
    gmx::EnergyHandler eh;

    eh.addAnalysisTool(ptr);

    eh.prepare(argc, argv);

    return eh.readFiles();
}

static int gmx_viscosity(int argc, char *argv[])
{
    return gmx_etool(&argc, argv, new(gmx::Viscosity));
}

static int gmx_fluctprops(int argc, char *argv[])
{
    return gmx_etool(&argc, argv, new(gmx::FluctProps));
}

static int gmx_dhdl(int argc, char *argv[])
{
    return gmx_etool(&argc, argv, new(gmx::DhdlEnergy));
}

static int gmx_fediff(int argc, char *argv[])
{
    return gmx_etool(&argc, argv, new(gmx::FreeEnergyDifference));
}

static int gmx_energy(int argc, char *argv[])
{
    return gmx_etool(&argc, argv, new(gmx::SimpleEnergy));
}

void registerEnergyAnalysisModules(gmx::CommandLineModuleManager *manager)
{
    manager->addModuleCMain("dhdl",
                            "Extract dH/dlambda terms from energy file",
                            &gmx_dhdl);
    manager->addModuleCMain("energy",
                            "Writes energies to xvg files and display averages",
                            &gmx_energy);
    manager->addModuleCMain("fediff",
                            "Compute free energy estimate from two energy files",
                            &gmx_fediff);
    manager->addModuleCMain("fluctprops",
                            "Compute heat capacities, and other fluctuation properties",
                            &gmx_fluctprops);
    manager->addModuleCMain("viscosity",
                            "Determines the shear- and bulk viscosity",
                            &gmx_viscosity);
}
