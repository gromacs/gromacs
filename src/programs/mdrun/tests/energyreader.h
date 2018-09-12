/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief Interfaces of related classes for tests that want to
 * inspect energies produced by mdrun.
 *
 * The responsibilities for reading and sharing information from .edr
 * files in an exception-safe manner are shared between two
 * classes. Intended usage is to call openEnergyFileToReadFields() to
 * return an EnergyFrameReader that knows it expects to read the named
 * .edr fields from each frame. Successive calls to its
 * EnergyFrameReader::readNextFrame() and EnergyFrameReader::frame()
 * methods will return EnergyFrame objects that contain the values for
 * the .edr fields registered earlier.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_PROGRAMS_MDRUN_TESTS_ENERGYREADER_H
#define GMX_PROGRAMS_MDRUN_TESTS_ENERGYREADER_H

#include <cstdint>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "gromacs/fileio/enxio.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/testasserts.h"

namespace gmx
{

class EnergyFrame;

namespace test
{

/*! \brief Convenience function to get std::string keys from a map.
 *
 * This function can be used to provide an input for
 * openEnergyFileToReadFields().
 *
 * \todo This returns a copy of the keys, which is convenient, but
 * inefficient. Alternatively, this could return a view of the keys
 * from a range rather than a container, but there's no implementation
 * of that in C++11 at the moment. */
template <typename Map>
std::vector<std::string> getKeys(const Map &m)
{
    std::vector<std::string> keys;
    keys.reserve(m.size());
    for (const auto &it : m)
    {
        keys.push_back(it.first);
    }
    return keys;
}

class EnergyFrameReader;
//! Convenience smart pointer typedef
typedef std::unique_ptr<EnergyFrameReader> EnergyFrameReaderPtr;

/*! \brief Open the file and return an object that can read the required fields from frames of an .edr file.
 *
 * \param[in] filename                 Name of the energy file to use
 * \param[in] requiredEnergyFieldNames Names of the energy fields that the caller requires to
 *                                     be present for an .edr file frame to be considered valid
 * \throws    FileIOError              If the .edr file cannot be opened
 * \throws    APIError                 If any required energy field is not present in the file
 * \throws    std::bad_alloc           When out of memory
 *
 * This function is intended to have the main responsibility for
 * making EnergyFrameReader objects. */
EnergyFrameReaderPtr openEnergyFileToReadFields(const std::string              &filename,
                                                const std::vector<std::string> &requiredEnergyFieldNames);

//! Convenience smart pointer typedef
typedef unique_cptr<ener_file, done_ener_file> ener_file_ptr;
//! Helper function to free resources (NB free_enxframe only frees the contents, not the pointer itself)
void done_enxframe(t_enxframe *fr);
//! Convenience smart pointer typedef
typedef unique_cptr<t_enxframe, done_enxframe> enxframe_ptr;

/*! \internal
 * \brief Manages returning an EnergyFrame containing required energy
 * field values read from successive frames of an .edr file. */
class EnergyFrameReader
{
    public:
        /*! \brief Attempt to read the next frame from the energy file.
         *
         * \return Whether a frame was available to read.
         *
         * If true is returned, then frame() should be called
         * to get access to the data. If false is returned, then no
         * further data exists and no further call to
         * readNextFrame() or frame() should occur.
         *
         * \throws APIError  if an earlier probe has not been properly handled
         *                   (by calling frame(), or stopping trying to read
         *                   from the file). */
        bool readNextFrame();
        /*! \brief Make an EnergyFrame from the contents of the next frame in the energy file.
         *
         * If the next frame has not been probed for, then probe for
         * it. If no next frame exists, then throw APIError, because
         * user code should have called readNextFrame() itself if this
         * is possible. (This permits user code to avoid making calls
         * to readNextFrame() in a case where it already knows that
         * the frame exists.)
         *
         * \throws APIError        if no next frame exists.
         * \throws std::bad_alloc  when out of memory. */
        EnergyFrame frame();
        /*! \brief Constructor
         *
         * \param[in] indicesOfEnergyFields  Looks up energy fields by name to get the index into a t_enxframe structure read by the legacy API.
         * \param[in] energyFile             Open energy file object to manage, and from which to read frames */
        explicit EnergyFrameReader(const std::map<std::string, int> &indicesOfEnergyFields,
                                   ener_file *energyFile);
    private:
        //! Convert energy field name to its index within a t_enxframe from this file.
        std::map<std::string, int> indicesOfEnergyFields_;
        //! Owning handle of an open energy file ready to read frames.
        const ener_file_ptr        energyFileGuard_;
        //! Owning handle of contents of .edr file frame after reading.
        const enxframe_ptr         enxframeGuard_;
        //! Whether the API has been used properly (ie. probe before reading).
        bool                       haveProbedForNextFrame_;
        //! Whether there has been a probe that found a next frame.
        bool                       nextFrameExists_;

        // Multiple owners of these resources isn't very sensible, so prevent it
        GMX_DISALLOW_COPY_AND_ASSIGN(EnergyFrameReader);
};

}  // namespace test
}  // namespace gmx

#endif
