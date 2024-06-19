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
/*! \defgroup module_coordinateio Handling of writing new coordinate files
 * \ingroup group_analysismodules
 * \brief
 * Provides basic functions to handle writing of new coordinate files.
 *
 * The methods included in the coordinateio module implement the basics
 * for manipulating and writing coordinate trajectory files
 * and changing metadata in the underlying datastructures. Included are a container for storing
 * modules that change trajectory data, as well as a manager to write output files that uses
 * those methods. It can be used from within \ref module_trajectoryanalysis, and uses
 * methods from:
 * - \ref module_options
 * - \ref module_selection
 *
 * <H3>Overview</H3>
 * The methods in coordinateio provide the infrastructure to perform operations on coordinate data files
 * and structures during data analysis. It implements ways to change the information
 * in coordinate data structures as well as checking that both input data and output
 * method are matching for a given coordinate file writing operation. For this
 * components verify first that all the requirements can be satisfied. Then
 * components are build that will change the coordinate information accordingly.
 *
 * The main parts are the outputadapters implemented using the
 * IOutputAdapter interface to change information in a local (deep) copy of t_trxframes
 * stored in the coordinatefile.
 *
 * <H3>Outputadapter</H3>
 * Each OutputAdapter module implements the same IOutputAdapter interface and
 * has to set its requirements for final
 * processing as a flag from the enum in requirementflags. During processing, they implement a custom
 * version of the processFrame directive that modifies the stored trajectory data before writing
 * a new file to disk.
 *
 *
 * The interaction between the CoordinateFile and the OutputAdapter modules derived from
 * IOutputAdapter is shown in the diagram below.
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   analysistool,
   builder [ label="CoordinateFileBuilder" ],
   coordinatefile [ label="CoordinateFile" ],
   container [ label="OutputAdapterStorage" ],
   outputadapters [ label="OutputAdapters" ];

   analysistool => builder [ label="Requests new coordinate output" ];
   analysistool => builder [ label="Specifies required OutputAdapters" ];
   builder => outputadapters [ label="Tries to construct new outputadapters" ];
   outputadapters => builder [ label="Return or give error for wrong preconditions" ];
   outputadapters => container [ label="Outputadapters are stored" ];
   container => builder [ label="Gives error if storage conditions are violated" ];
   builder => coordinatefile [ label="Constructs new manager according to specifications" ];
   builder => container [ label="Requests storage object with registered outputadapters" ];
   container => builder [ label="Gives ownership of stored outputadapters" ];
   builder box builder [ label="Tests preconditions of storage object and new coordinatefile" ];
   builder => analysistool [ label="Raise error if preconditions don't match" ];
   builder => coordinatefile [ label="Add storage object to new coordinatefile" ];
   coordinatefile => analysistool [ label="Returns finished coordinatefile" ];
   builder box builder [ label="coordinatefile created, can start to work on input data" ];

 * \endmsc
 *
 * Once the CoordinateFile object and its registered modules are created, they can be used to
 * iterate over input data to write new coordinate frames.
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   analysistool,
   analysisloop,
   coordinatefile [ label="CoordinateFile" ],
   outputadapters [ label="OutputAdapters" ] ,
   filewriting;

   --- [ label="Setup of coordinatefile complete, analysis phase begins" ];
    analysistool   => analysisloop [ label="Starts iteration over frames" ];
    analysisloop   => coordinatefile [ label="Provides coordinates" ];
    coordinatefile  => outputadapters [ label="Provide coordinate frames for changing" ];
    outputadapters => coordinatefile [ label="Return after changing data" ];
    coordinatefile  => filewriting [ label="Send new coordinates for writing" ];
    filewriting    => coordinatefile [ label="Continue after writing to disk" ];
    coordinatefile  => analysisloop [ label="Returns after writing" ];
    analysisloop box analysisloop [ label="Iterates over frames" ];
    --- [ label="Analysis complete, object is destructed and files are closed" ];

 *  \endmsc
 *
 *
 *
 * \if libapi
 * <H3>Preparing new OutputAdapters</H3>
 *
 * If additional methods are needed to perform changes to the t_trxframe metadata,
 * new OutputAdapters can be written that again implement the IOutputAdapter interface.
 * The new method should follow the approach of the other modules that are present
 * in changing a minimal set of t_trxframe data.
 * \endif
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
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

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/outputadaptercontainer.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"

struct gmx_mtop_t;
struct t_trxstatus;
struct t_trxframe;

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
