/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 *
 * \brief This file declares implementation types for handleRestart()
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrunutility
 */

#ifndef GMX_MDRUN_HANDLERESTART_IMPL_H
#define GMX_MDRUN_HANDLERESTART_IMPL_H

#include <string.h>

#include <string>
#include <vector>

#include "gromacs/fileio/fileiohandlerinterface.h"
//#include "gromacs/fileio/gmxfio.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/mpihandlerinterface.h"

namespace gmx
{

namespace HandleRestart
{

//! Convenience typedef
typedef ConstArrayRef<t_filenm> ConstCommandLineFilenames;

/*! \internal
 * \brief POD class for data from checkpoint file for restarting
 *
 * Used internally by handleRestart() and methods of Impl. In
 * particular, its constructor implements default initialization of
 * values that might be read from the checkpoint file, which is useful
 * while testing Impl::getDataFromCheckpoint().
 */
class DataFromCheckpoint
{
    public:
        //! Constructor
        DataFromCheckpoint() : partNumber_(0),
                               expectedOutputFilenames_() {}

        //! The number of the simulation part in the checkpoint file
        int                      partNumber_;
        //! Container of filenames from the checkpoint
        std::vector<std::string> expectedOutputFilenames_;
};


/*! \internal
 * \brief Implementation class for handleRestart()
 *
 * A class is used so that functions that need to be mocked while
 * testing functions that call them can be virtual, and thus easy to
 * mock.
 */
class Impl
{
    public:
        /*! \brief Constructor
         *
         * \param[in]    fileIOHandler   Handle to object that injects behaviour dependent on file systems or filenames. See gmx::FileIOHandler.
         * \param[in]    mpiHandler      Handle to object that injects behaviour dependent on MPI. See gmx::MpiHandlerHandler.
         * \param[inout] mdrunFilenames  Container of filenames parsed from mdrun command line
         */
        Impl(FileIOHandlerInterface   *fileIOHandler,
             MpiHandlerInterface      *mpiHandler,
             ConstCommandLineFilenames mdrunFilenames);

        //! Virtual destructor required
        virtual ~Impl();

        /*! \brief Get data from mdrun -cpi if appropriate
         *
         * If a checkpoint file was supplied to mdrun and this is a
         * master rank, read that file (if it exists and can be read),
         * to get the old simulation part number and names of files
         * used in that run. This anticipates communication of
         * partNumber_ to achieve consistency within and between
         * running simulations.
         *
         * \return An object whose fields contain valid data read
         * from a checkpoint file, or default data if reading did not
         * take place.
         *
         * \param[in]    bIsSimMaster    Flag to permit only one rank to write to read the checkpoint file; should normally be set to SIMMASTER(cr)
         * \param[in]    bIsMultiMaster  Flag to permit only one rank to write to stdout; should normally be set to MULTIMASTER(cr)
         *
         * \throws std::bad_alloc          if out of memory
         * \throws InconsistentInputError  if checkpoint file is missing and a log file exists
         * \throws InternalError           if the checkpoint file contains no output filenames (ie. is malformed)
         */
        DataFromCheckpoint getDataFromCheckpoint(const bool bIsSimMaster,
                                                 const bool bIsMultiMaster) const;
        /*! \brief Search for \p filename in \p mdrunFilenames_ and return true iff it is found, is an outputfile, and exists
         *
         * Does not throw */
        virtual bool
        outputFileExists(const std::string &filename) const;
        /*! \brief Construct an error message when the files named in the checkpoint is only a partial match with what is on disk
         *
         * When only some of the set of files in the checkpoint are
         * present on disk, we cannot proceed with appending, and wish
         * to issue a detailed error message about what is wrong.
         *
         * \param[in] expectedOutputFilenames  Container of filenames from the checkpoint
         * \param[in] numFilesThatExist        Number of files found on disk
         * \return                             The error message
         *
         * \throws std::bad_alloc              if out of memory
         */
        virtual std::string
        makeMessageWhenSetsOfFilesDontMatch(const std::vector<std::string> &expectedOutputFilenames,
                                            size_t                          numFilesThatExist) const;
        /*! \brief Throw unless all the mdrun files for appending exist
         *
         * Check if the output filenames stored in the checkpoint file
         * match the output filenames of mdrun, and actually exist.
         *
         * \throws std::bad_alloc          if out of memory
         * \throws InconsistentInputError  if not all of the output filenames in the checkpoint files are found on disk
         */
        virtual void
        checkAllFilesForAppendingExist(const std::vector<std::string> &expectedOutputFilenames) const;
        /*! \brief Return whether the conditions on this rank permit appending
         *
         * If using mdrun -append and the old simulation part number
         * is > 0, appending might occur. Check that the files on disk
         * form a complete/empty set, and prepare to return true/false
         * respectively. If only some output files are present, give a
         * fatal error. However, if the .log filename in the
         * checkpoint file contains ".part", then the files on
         * disk have part numbers, so we cannot append any files
         * because that is not supported behaviour.
         *
         * \return On master rank, when the conditions are correct for
         * the run to actually do appending, return true; otherwise
         * return false. It is expected that non-master ranks will get
         * the information via broadcast.
         *
         * \param[in]  dataFromCheckpoint  Data read from the checkpoint
         * \param[in]  bTryToAppendFiles   Whether mdrun -append was used
         *
         * \throws std::bad_alloc  if out of memory
         * \throws InternalError   if the contents of the list of output filenames in the checkpoint file do not conform to expectations from how they should be written
         * Also throws anything thrown by checkAllFilesForAppendingExist() */
        bool canAppend(const DataFromCheckpoint &dataFromCheckpoint,
                       const bool                bTryToAppendFiles) const;

        /*! \brief Injected dependency for external functions that
         * refer to file names or the file system. Mock object used in
         * testing. */
        FileIOHandlerInterface *fileIOHandler_;
        /*! \brief Injected dependency for external functions that
         * refer to MPI functionality controlled by the contents of a
         * t_commrec. Mock object used in testing. */
        MpiHandlerInterface      *mpiHandler_;
        /*! \brief Container of filenames on the mdrun comand line */
        ConstCommandLineFilenames mdrunFilenames_;
};


} // namespace

} // namespace

#endif
