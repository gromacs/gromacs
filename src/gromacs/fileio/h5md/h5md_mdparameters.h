/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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

/*! \internal \file
 * \brief Declares the description of the MD parameters written to H5MD files.
 *
 * The MD parameters are described with engine-agnostic names and units, following the
 * MD parameters of the H5MD biophysics module, as one object per category (run control,
 * output control, temperature and pressure coupling, Coulomb and van der Waals interactions,
 * boundary conditions and constraints). Each parameter is an object holding a "value" and,
 * for a quantity with a unit, a "unit". A tensor is held row by row, with a "shape" of [3, 3].
 *
 * Only the settings that act in the simulation are described: a category or subcategory
 * that does not apply is omitted, and so is a parameter without a value.
 *
 * The topology provides what the input record does not hold: the names of the simulation
 * groups and whether any bond is constrained.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_H5MD_MDPARAMETERS_H
#define GMX_FILEIO_H5MD_MDPARAMETERS_H

struct gmx_mtop_t;
struct t_inputrec;

namespace gmx
{

class KeyValueTreeObjectBuilder;

/*! \brief Describe the MD parameters of a simulation in \p builder.
 *
 * The parameters are described as set out in the file documentation.
 *
 * \param[in,out] builder     The object to add the categories to.
 * \param[in]     inputRecord The input record of the simulation.
 * \param[in]     mtop        The topology of the simulation.
 * \param[in]     isRerun     Whether the energies of an existing trajectory are recomputed,
 *                            rather than a simulation being run.
 *
 * \throws std::bad_alloc If out of memory.
 */
void buildMdParameterTree(KeyValueTreeObjectBuilder* builder,
                          const t_inputrec&          inputRecord,
                          const gmx_mtop_t&          mtop,
                          bool                       isRerun);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_MDPARAMETERS_H
