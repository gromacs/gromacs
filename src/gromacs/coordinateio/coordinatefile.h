/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \file
 * \brief
 * CoordinateFile takes care of opening files and writing output to them.
 *
 * \author
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_COORDINATEFILE_H
#define GMX_COORDINATEIO_COORDINATEFILE_H

#include <string>
#include <utility>

#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/outputadaptercontainer.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"

struct gmx_mtop_t;
struct t_trxstatus;

namespace gmx
{

class Selection;
class TrajectoryFrameWriter;
struct OutputRequirements;

/*! \brief
 * Factory function for TrajectoryFrameWriter.
 *
 * Used to initialize a new instance of TrajectoryFrameWriter with the user supplied information
 * for writing trajectory data to disk. Information needed is the file type, file name
 * corresponding to the type, if available topology information and selection information.
 *
 * If supplied, the modules contained within \p adapters are registered on
 * the TrajectoryFrameWriter if possible.
 *
 * The factory function is responsible for the initial santity checks concerning file types and
 * availability of topology information, with the registration of modules being the second part.
 *
 * \param[in] top                    Pointer to full topology or null.
 * \param[in] sel                    Reference to global selection used to construct the object.
 * \param[in] filename               Name of new output file, used to deduce file type.
 * \param[in] atoms                  Smart Pointer to atoms data or null.
 * \param[in] requirements           Container for settings obtained to specify which
 *                                   OutputAdapters should be registered.
 * \throws    InconsistentInputError When user input and requirements don't match.
 */
std::unique_ptr<TrajectoryFrameWriter> createTrajectoryFrameWriter(const gmx_mtop_t*  top,
                                                                   const Selection&   sel,
                                                                   const std::string& filename,
                                                                   AtomsDataPtr       atoms,
                                                                   OutputRequirements requirements);

/*!\brief
 * Low level method to take care of only file opening and closing.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class TrajectoryFileOpener
{
public:
    /*! \brief
     * Constructor, taking all arguments needed to open valid coordinate files of any type.
     *
     * \param[in] name Name of the file to create.
     * \param[in] filetype Internal filetype to know which kind we are going to have.
     * \param[in] sel Reference to selection of atoms to write. Needs to be valid for
     *                longer than the lifetime of the object created here.
     * \param[in] mtop Topology used to create TNG output. Needs to be valid for longer
     *                 than the object created here.
     */
    TrajectoryFileOpener(const std::string& name, int filetype, const Selection& sel, const gmx_mtop_t* mtop) :
        outputFileName_(name), outputFile_(nullptr), filetype_(filetype), sel_(sel), mtop_(mtop)
    {
    }

    /*! \brief
     * Closes new trajectory file after finishing the writing to it.
     */
    ~TrajectoryFileOpener();

    /*! \brief
     * Get access to initialized output file object.
     *
     * Performs lazy initialization if needed.
     */
    t_trxstatus* outputFile();

private:
    //! Name for the new coordinate file.
    std::string outputFileName_;

    //! File pointer to the coordinate file being written.
    t_trxstatus* outputFile_;

    //! Storage of file type for determing what kind of file will be written to disk.
    int filetype_;

    /*! \brief
     * Selection of atoms to write out.
     *
     * Currently, CoordinateFile expects that the lifetime of the selection is longer
     * than that of itself, and that the selection remains unchanged during this lifetime.
     * A better approach will be to pass the selection to it and expect it to
     * manage the lifetime instead.
     */
    const Selection& sel_;

    //! Pointer to topology information if available.
    const gmx_mtop_t* mtop_;
};

/*!\brief
 * Writes coordinate frames to a sink, e.g. a trajectory file.
 *
 * Writes all frames passed to the trajectory file handle it
 * manages. It can use IOutputAdapter modules to transform the
 * frame data before writing. If any transform modules are used,
 * makes a deep copy of the frame contents.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class TrajectoryFrameWriter
{
public:
    friend std::unique_ptr<TrajectoryFrameWriter> createTrajectoryFrameWriter(const gmx_mtop_t* top,
                                                                              const Selection&  sel,
                                                                              const std::string& filename,
                                                                              AtomsDataPtr atoms,
                                                                              OutputRequirements requirements);

    /*! \brief
     * Writes the input frame, after applying any IOutputAdapters.
     *
     * \param[in] framenumber Number of frame being currently processed.
     * \param[in] input View of the constant t_trxframe object provided by the
     *                  method that calls the output manager.
     */
    void prepareAndWriteFrame(int framenumber, const t_trxframe& input);

private:
    /*! \brief
     * Creates fully initialized object.
     *
     * \param[in] name Name of the file to create.
     * \param[in] filetype Internal filetype to know which kind we are going to have.
     * \param[in] sel Reference to selection of atoms to write. Needs to be valid for
     *                longer than the lifetime of the object created here.
     * \param[in] mtop Topology used to create TNG output. Needs to be valid for longer
     *                 than the object created here.
     * \param[in] adapters Container of methods that can modify the information written
     *                     to the new file.
     */
    TrajectoryFrameWriter(const std::string&     name,
                          int                    filetype,
                          const Selection&       sel,
                          const gmx_mtop_t*      mtop,
                          OutputAdapterContainer adapters) :
        file_(name, filetype, sel, mtop), outputAdapters_(std::move(adapters))
    {
    }

    //! Underlying object for open/writing to file.
    TrajectoryFileOpener file_;

    //! Storage for list of output adapters.
    OutputAdapterContainer outputAdapters_;

    //! Local storage for modified positions.
    std::vector<RVec> localX_;
    //! Local storage for modified velocities.
    std::vector<RVec> localV_;
    //! Local storage for modified forces.
    std::vector<RVec> localF_;
    //! Local storage for modified indices.
    std::vector<int> localIndex_;
};

//! Smart pointer to manage the TrajectoryFrameWriter object.
using TrajectoryFrameWriterPointer = std::unique_ptr<TrajectoryFrameWriter>;

} // namespace gmx

#endif
