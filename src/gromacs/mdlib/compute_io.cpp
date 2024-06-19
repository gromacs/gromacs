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

#include "compute_io.h"

#include <csignal>
#include <cstdlib>

#include <memory>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

static int div_nsteps(int nsteps, int nst)
{
    if (nst > 0)
    {
        return gmx::divideRoundUp(1 + nsteps, nst);
    }
    else
    {
        return 0;
    }
}

double compute_io(const t_inputrec* ir, int natoms, const SimulationGroups& groups, int nrener, int nrepl)
{

    int    nsteps    = ir->nsteps;
    int    nxtcatoms = 0;
    int    nstx, nstv, nstf, nste, nstlog, nstxtc;
    double cio;

    nstx   = div_nsteps(nsteps, ir->nstxout);
    nstv   = div_nsteps(nsteps, ir->nstvout);
    nstf   = div_nsteps(nsteps, ir->nstfout);
    nstxtc = div_nsteps(nsteps, ir->nstxout_compressed);
    if (ir->nstxout_compressed > 0)
    {
        for (int i = 0; i < natoms; i++)
        {
            if (groups.groupNumbers[SimulationAtomGroupType::CompressedPositionOutput].empty()
                || groups.groupNumbers[SimulationAtomGroupType::CompressedPositionOutput][i] == 0)
            {
                nxtcatoms++;
            }
        }
    }
    nstlog = div_nsteps(nsteps, ir->nstlog);
    /* We add 2 for the header */
    nste = div_nsteps(2 + nsteps, ir->nstenergy);

    cio = 80 * natoms;
    cio += (nstx + nstf + nstv) * sizeof(real) * (3.0 * natoms);
    cio += nstxtc * (14 * 4 + nxtcatoms * 5.0); /* roughly 5 bytes per atom */
    cio += nstlog * (nrener * 16 * 2.0);        /* 16 bytes per energy term plus header */
    /* t_energy contains doubles, but real is written to edr */
    cio += (1.0 * nste) * nrener * 3 * sizeof(real);

    if ((ir->efep != FreeEnergyPerturbationType::No || ir->bSimTemp) && (ir->fepvals->nstdhdl > 0))
    {
        int ndh    = ir->fepvals->n_lambda;
        int ndhdl  = 0;
        int nchars = 0;

        for (auto i : keysOf(ir->fepvals->separate_dvdl))
        {
            if (ir->fepvals->separate_dvdl[i])
            {
                ndhdl += 1;
            }
        }

        if (ir->fepvals->separate_dhdl_file == SeparateDhdlFile::Yes)
        {
            nchars = 8 + ndhdl * 8 + ndh * 10; /* time data ~8 chars/entry, dH data ~10 chars/entry */
            if (ir->expandedvals->elmcmove > LambdaMoveCalculation::No)
            {
                nchars += 5; /* alchemical state */
            }

            if (ir->fepvals->edHdLPrintEnergy != FreeEnergyPrintEnergy::No)
            {
                nchars += 12; /* energy for dhdl */
            }
            cio += div_nsteps(nsteps, ir->fepvals->nstdhdl) * nchars;
        }
        else
        {
            /* dH output to ener.edr: */
            if (ir->fepvals->dh_hist_size <= 0)
            {
                int ndh_tot = ndh + ndhdl;
                if (ir->expandedvals->elmcmove > LambdaMoveCalculation::No)
                {
                    ndh_tot += 1;
                }
                if (ir->fepvals->edHdLPrintEnergy != FreeEnergyPrintEnergy::No)
                {
                    ndh_tot += 1;
                }
                /* as data blocks: 1 real per dH point */
                cio += div_nsteps(nsteps, ir->fepvals->nstdhdl) * ndh_tot * sizeof(real);
            }
            else
            {
                /* as histograms: dh_hist_size ints per histogram */
                cio += div_nsteps(nsteps, ir->nstenergy) * sizeof(int) * ir->fepvals->dh_hist_size * ndh;
            }
        }
    }
    if (ir->pull != nullptr)
    {
        cio += div_nsteps(nsteps, ir->pull->nstxout) * 20; /* roughly 20 chars per line */
        cio += div_nsteps(nsteps, ir->pull->nstfout) * 20; /* roughly 20 chars per line */
    }

    return cio * nrepl / (1024 * 1024);
}
