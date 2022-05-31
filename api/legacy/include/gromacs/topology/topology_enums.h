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
#ifndef GMX_TOPOLOGY_TOPOLOGY_ENUMS_H
#define GMX_TOPOLOGY_TOPOLOGY_ENUMS_H

enum class SimulationAtomGroupType : int
{
    TemperatureCoupling,
    EnergyOutput,
    Acceleration,
    Freeze,
    User1,
    User2,
    MassCenterVelocityRemoval,
    CompressedPositionOutput,
    OrientationRestraintsFit,
    QuantumMechanics,
    Count
};

//! Short strings used for describing atom groups in log and energy files
const char* shortName(SimulationAtomGroupType type);

/* The particle type */
enum class ParticleType : int
{
    Atom,
    Nucleus,
    Shell,
    Bond,
    VSite,
    Count
};

/* The particle type names */
const char* enumValueToString(ParticleType enumValue);

/* Enumerated type for pdb records. The other entries are ignored
 * when reading a pdb file
 */
enum class PdbRecordType : int
{
    Atom,
    Hetatm,
    Anisou,
    Cryst1,
    Compound,
    Model,
    EndModel,
    Ter,
    Header,
    Title,
    Remark,
    Conect,
    Count
};

const char* enumValueToString(PdbRecordType enumValue);

#endif
