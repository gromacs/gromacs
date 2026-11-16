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

/*! \libinternal \file
 * \brief Declares the writing of a KeyValueTree to H5MD.
 *
 * A key-value tree has no notion of HDF5 datasets and attributes, so the writer maps it
 * onto HDF5 by the structure of each object and a few reserved keys:
 *  - An object with a "value" becomes a dataset of that name, holding the value.
 *  - Any other object becomes a group of that name.
 *  - In either, a reserved attribute key becomes an attribute of the dataset or group.
 *  - In a group, a nested object is written in the same way, to any depth. A dataset cannot
 *    hold children, so an object with a "value" cannot hold nested objects.
 *  - Any other key in an object is an error. At the top level of the tree, a scalar or an
 *    array is ignored.
 *
 * Reserved keys that describe the dataset:
 *  - "value": the data, a scalar or an array;
 *  - "shape": the dimensions of an array "value" in row-major order, e.g. [3, 3] for a
 *    tensor, as positive integers whose product equals the number of values. Without it,
 *    an array has one dimension.
 *
 * Reserved keys written as attributes:
 *  - "unit":        physical unit of the value;
 *  - "description": description of the element;
 *  - "doi":         DOI of the publication that describes the method;
 *  - "url":         URL of further documentation;
 *  - "version":     version of the method or of its implementation.
 *
 * \verbatim
 * Example key-value tree:
 * {
 *     temperature_coupling: {
 *         reference_temperature: {
 *             value: [300, 300],
 *             unit: "K",
 *             description: "Simulation temperature"
 *         },
 *         stochastic_velocity_rescaling: {
 *             doi: "xxx.xxx",
 *             coupling_time: {
 *                 value: [0.1, 0.1],
 *                 unit: "ps",
 *             }
 *         }
 *     }
 * }
 * is written as the HDF5 structure:
 * container/
 *     \-- temperature_coupling (group)
 *         \-- reference_temperature (dataset, value=[300, 300])
 *         |       +-- unit="K" (attribute)
 *         |       +-- description="Simulation temperature" (attribute)
 *         \-- stochastic_velocity_rescaling (group)
 *                 +-- doi="xxx.xxx" (attribute)
 *                 \-- coupling_time (dataset, value=[0.1, 0.1])
 *                         +-- unit="ps" (attribute)
 * \endverbatim
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \inlibraryapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_H5MD_KVTREE_WRITER_H
#define GMX_FILEIO_H5MD_KVTREE_WRITER_H

#include <hdf5.h>

namespace gmx
{

class KeyValueTreeObject;

/*! \brief Write a KeyValueTree to H5MD.
 *
 * The properties of \p tree are written as set out in the file documentation.
 *
 * \param[in] container HDF5 group or file handle to write into
 * \param[in] tree      KeyValueTree to write
 *
 * \throws H5mdError If an object holds a key that is not allowed, if an object with a "value"
 *                   holds a nested object, if a "value" is of a type that is not supported,
 *                   if an array mixes element types, if a "shape" is invalid or does not
 *                   match its "value", or if writing to HDF5 fails.
 */
void writeKeyValueTreeToH5md(hid_t container, const KeyValueTreeObject& tree);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_KVTREE_WRITER_H
