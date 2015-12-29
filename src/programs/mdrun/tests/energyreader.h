/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * The responsibilities for reading information from .edr files are
 * shared between several classes. Typical usage is to create an
 * EnergyFileReader, and call its EnergyFileReader::openToReadFields()
 * method to return an EnergyFrameReader that knows it expects to read
 * the named .edr frame fields. Successive calls to its
 * EnergyFrameReader::probeForNextFrame() and
 * EnergyFrameReader::makeFrameInfo() methods will return
 * EnergyFrameInfo objects that contain the values for the .edr fields
 * registered earlier.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_PROGRAMS_MDRUN_TESTS_ENERGYREADER_H
#define GMX_PROGRAMS_MDRUN_TESTS_ENERGYREADER_H

#include <map>
#include <string>
#include <vector>

#include "gromacs/fileio/enxio.h"
#include "gromacs/utility/scoped_cptr.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

class EnergyFrameReader;

/*! \internal
 * \brief Manages opening an .edr file and reading required energy fields from it safely.
 *
 * This class manages resources for reading .edr file frames, to
 * ensure that objects of EnergyFrameReader and EnergyFrameInfo can be
 * used with as much exception safety as possible (given the
 * limitations of the legacy API it wraps).
 *
 * As such, it has sole responsibility for making EnergyFrameReader objects. */
class EnergyFileReader
{
    public:
        /*! \brief Constructor, taking the name of the .edr file to read. */
        explicit EnergyFileReader(const std::string &filename);
        /*! \brief Open the file and return an object that can read the required fields from an .edr file.
         *
         * \param[in] requiredEnergyFieldNames Names of the energy fields that the caller requires to
         *                                     be present for an .edr file frame to be considered valid
         * \throws    FileIOError              If the .edr file cannot be opened
         * \throws    APIError                 Upon any subsequent call to this function
         * \throws    APIError                 If any required energy field is not present in the file
         * \throws    std::bad_alloc           When out of memory */
        EnergyFrameReader openToReadFields(const std::vector<std::string> &requiredEnergyFieldNames);

    private:
        //! Name of .edr file.
        std::string                       filename_;
        //! Legacy file object managed for exception safety.
        scoped_cptr<ener_file, close_enx> energyFile_;
        //! Object to contain data read from .edr file frames, used by associated EnergyFrameReader object.
        t_enxframe frame_;
        //! RAII-aware handle to contents of .edr file frames, to permit exception-safe use of EnergyFrameReader.
        scoped_cptr<t_enxframe, free_enxframe> frameContents_;
        // Multiple owners of these resources isn't very sensible, so prevent it
        GMX_DISALLOW_COPY_AND_ASSIGN(EnergyFileReader);
};

class EnergyFrameInfo;

/*! \internal
 * \brief Manages returning an EnergyFrameInfo object containing
 * required energy field values read from successive frames of an .edr
 * file.
 *
 * The scope of such an object should not exceed that of the
 * EnergyFileReader object that created it. */
class EnergyFrameReader
{
    public:
        /*! \brief Return whether a next frame exists in the energy file.
         *
         * \return Whether a next frame exists.
         *
         * If true is returned, then makeFrameInfo() should be called
         * to get access to the data. If false is returned, then no
         * further data exists and no further call to
         * probeForNextFrame() or makeFrameInfo() should occur.
         *
         * \throws APIError  if an earlier probe has not been properly handled
         *                   (by calling makeFrameInfo, or stopping trying to read
         *                   from the file). */
        bool probeForNextFrame();
        /*! \brief Make an EnergyFrameInfo from the contents of the next frame in the energy file.
         *
         * If the next frame has not been probed for, then probe for
         * it. If no next frame exists, then throw APIError, because
         * user code should have called probeForNextFrame() itself if
         * this is possible. (This permits user code to avoid making
         * calls to probeForNextFrame() in cases where it already
         * knows that the frame exists.)
         *
         * \throws APIError  If no next frame exists.
         * \throws std::bad_alloc  when out of memory */
        EnergyFrameInfo makeFrameInfo();
    private:
        /*! \brief Constructor (private, so only EnergyFileReader can call it)
         *
         * \param[in] indicesOfEnergyFields  Looks up energy fields by name to get the index into a t_enxframe structure read by the legacy API.
         * \param[in] energyFile             Energy file object from which to read frames
         * \param[in] frame                  Handle to object into the legacy energy-frame reading API stores temporary data */
        explicit EnergyFrameReader(std::map<std::string, int> indicesOfEnergyFields,
                                   ener_file *energyFile,
                                   t_enxframe *frame);
    private:
        //! Convert energy field name to its index within a t_enxframe.
        std::map<std::string, int> indicesOfEnergyFields_;
        //! Handle to open energy file
        ener_file                 *energyFile_;
        //! Handle to object with data read from the current .edr file frame.
        t_enxframe                *frame_;
        //! Whether the API has been used properly (ie. probe before reading).
        bool                       haveProbedForNextFrame_;
        //! Whether there has been a probe that found a next frame.
        bool                       nextFrameExists_;

        friend class EnergyFileReader;
};

/*! \internal
 * \brief Contains the content of an .edr frame read by an EnergyFrameReader
 *
 * Objects of this type are intended to be constructed by
 * EnergyFrameReader objects, and as such will always contain valid
 * data from an .edr file frame. */
class EnergyFrameInfo
{
    public:
        /*! \brief Return string that helps users identify this frame, containing time and step number.
         *
         * \throws std::bad_alloc  when out of memory */
        std::string getFrameName() const;
        /*! \brief Return the value read for energy \c name.
         *
         * \throws APIError  if \c name was not registered with EnergyFileReader. */
        real getValue(const std::string &name) const;
    private:
        //! Constructor (private, so that only EnergyFrameReader can call it).
        EnergyFrameInfo();
    private:
        //! Container for energy values, indexed by name
        std::map<std::string, real> values_;
        //! Step number read from the .edr file frame
        gmx_int64_t                 step_;
        //! Time read from the .edr file frame
        double time_;

        friend class EnergyFrameReader;
};

} // namespace
} // namespace

#endif
