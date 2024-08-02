/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 *
 * \brief This file declares functions for mdrun to call to manage the
 * details of doing a restart (ie. reading checkpoints, appending
 * output files).
 *
 * \todo Clean up the error-prone logic here. Add doxygen.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrunutility
 */

#include "gmxpre.h"

#include "handlerestart.h"

#include "config.h"

#include <fcntl.h>

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <array>
#include <filesystem>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"
#if GMX_NATIVE_WINDOWS
#    include <io.h>

#    include <sys/locking.h>
#endif

#include <algorithm>
#include <exception>
#include <functional>
#include <optional>
#include <tuple>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{
namespace
{

/*! \brief Search for \p fnm_cp in fnm and return true iff found
 *
 * \todo This could be implemented sanely with a for loop. */
gmx_bool exist_output_file(const std::filesystem::path& fnm_cp, int nfile, const t_filenm fnm[])
{
    int i;

    /* Check if the output file name stored in the checkpoint file
     * is one of the output file names of mdrun.
     */
    i = 0;
    while (i < nfile && !(is_output(&fnm[i]) && fnm_cp == fnm[i].filenames[0]))
    {
        i++;
    }

    return (i < nfile && gmx_fexist(fnm_cp));
}

/*! \brief Throw when mdrun -cpi fails because previous output files are missing.
 *
 * If we get here, the user requested restarting from a checkpoint file, that checkpoint
 * file was found (so it is not the first part of a new run), but we are still missing
 * some or all checkpoint files. In this case we issue a fatal error since there are
 * so many special cases we cannot keep track of, and better safe than sorry. */
[[noreturn]] void throwBecauseOfMissingOutputFiles(const std::filesystem::path& checkpointFilename,
                                                   ArrayRef<const gmx_file_position_t> outputfiles,
                                                   int                                 nfile,
                                                   const t_filenm                      fnm[],
                                                   size_t numFilesMissing)
{
    StringOutputStream stream;
    TextWriter         writer(&stream);
    writer.writeLineFormatted(
            "Some output files listed in the checkpoint file %s are not present or not named "
            "as the output files by the current program:)",
            checkpointFilename.c_str());
    auto& settings  = writer.wrapperSettings();
    auto  oldIndent = settings.indent(), newIndent = 2;

    writer.writeLine("Expected output files that are present:");
    settings.setIndent(newIndent);
    settings.setLineLength(78);
    for (const auto& outputfile : outputfiles)
    {
        if (exist_output_file(outputfile.filename, nfile, fnm))
        {
            writer.writeLine(outputfile.filename);
        }
    }
    settings.setIndent(oldIndent);
    writer.ensureEmptyLine();

    // The implementation of -deffnm does not handle properly the
    // naming of output files that share a common suffix, such as
    // pullx.xvg and pullf.xvg from the pull module. Such output files
    // will be sought by the wrong name by the code that handles the
    // restart, even though the pull module would later work out what
    // they should have been called. Since there is a straightforward
    // way to work around that, we help the user with that. This can
    // be removed when gitlab issue #3875 is resolved.
    bool missingFilesIncludedPullOutputFiles = false;
    writer.writeLine("Expected output files that are not present or named differently:");
    settings.setIndent(newIndent);
    for (const auto& outputfile : outputfiles)
    {
        if (!exist_output_file(outputfile.filename, nfile, fnm))
        {
            writer.writeLine(outputfile.filename);
            // If this was a pull file, then we have a known issue and
            // work-around (See gitlab issue #3442).
            if (!missingFilesIncludedPullOutputFiles
                && (contains(outputfile.filename, "pullx")
                    || contains(outputfile.filename, "pullf")))
            {
                missingFilesIncludedPullOutputFiles = true;
            }
        }
    }
    if (missingFilesIncludedPullOutputFiles)
    {
        writer.ensureEmptyLine();
        writer.writeLineFormatted(
                "It appears that pull output files were not found. It is known that "
                "using gmx mdrun -deffnm test with pulling and later "
                "gmx mdrun -deffnm test -cpi will fail to consider the changed default "
                "filename when checking the pull output files for restarting with "
                "appending. You may be able to work around this by using a command like "
                "gmx mdrun -deffnm test -px test_pullx -pf test_pullf -cpi.");
    }
    settings.setIndent(oldIndent);

    writer.ensureEmptyLine();
    writer.writeLineFormatted(
            "To keep your simulation files safe, this simulation will not restart. "
            "Either name your output files exactly the same as the previous simulation "
            "part (e.g. with -deffnm or explicit naming), or make sure all the output "
            "files are present (e.g. run from the same directory as the previous simulation "
            "part), or instruct mdrun to write new output files with mdrun -noappend. In "
            "the last case, you will not be able to use appending in future for this "
            "simulation.",
            numFilesMissing,
            outputfiles.size());
    GMX_THROW(InconsistentInputError(stream.toString()));
}

//! Return a string describing the precision of a build of GROMACS.
const char* precisionToString(bool isDoublePrecision)
{
    return isDoublePrecision ? "double" : "mixed";
}

/*! \brief Describes how mdrun will (re)start and provides supporting
 * functionality based on that data. */
class StartingBehaviorHandler
{
public:
    /*! \brief Throw unless all simulations in a multi-sim restart the same way.
     *
     * The restart could differ if checkpoint or other output files are
     * not found in a consistent way across the set of multi-simulations,
     * or are from different simulation parts.
     *
     * \param[in]  ms      Multi-sim handler.
     *
     * May only be called from the main rank of each simulation.
     *
     * \throws InconsistentInputError if either simulations restart
     * differently, or from checkpoints from different simulation parts.
     */
    void ensureMultiSimBehaviorsMatch(const gmx_multisim_t* ms);
    /*! \brief Return an optional value that describes the index
     * of the next simulation part when not appending.
     *
     * \param[in] appendingBehavior  Whether this simulation is appending
     *                                (relevant only when restarting)
     *
     * Does not throw */
    std::optional<int> makeIndexOfNextPart(AppendingBehavior appendingBehavior) const;

    //! Describes how mdrun will (re)start
    StartingBehavior startingBehavior = StartingBehavior::NewSimulation;
    //! When restarting from a checkpoint, contains the contents of its header
    std::optional<CheckpointHeaderContents> headerContents;
    //! When restarting from a checkpoint, contains the names of expected output files
    std::optional<std::vector<gmx_file_position_t>> outputFiles;
};

/*! \brief Choose the starting behaviour for this simulation
 *
 * This routine cannot print tons of data, since it is called before
 * the log file is opened.
 *
 * Note that different simulations in a multi-simulation can return
 * values that depend on whether the respective checkpoint files are
 * found (and other files found, when appending), and so can differ
 * between multi-simulations. It is the caller's responsibility to
 * detect this and react accordingly. */
StartingBehaviorHandler chooseStartingBehavior(const AppendingBehavior appendingBehavior,
                                               const int               nfile,
                                               t_filenm                fnm[])
{
    StartingBehaviorHandler handler;
    if (!opt2bSet("-cpi", nfile, fnm))
    {
        // No need to tell the user anything
        handler.startingBehavior = StartingBehavior::NewSimulation;
        return handler;
    }

    // A -cpi option was provided, do a restart if there is an input checkpoint file available
    const char* checkpointFilename = opt2fn("-cpi", nfile, fnm);
    if (!gmx_fexist(checkpointFilename))
    {
        // This is interpreted as the user intending a new
        // simulation, so that scripts can call "gmx mdrun -cpi"
        // for all simulation parts. Thus, appending cannot occur.
        if (appendingBehavior == AppendingBehavior::Appending)
        {
            GMX_THROW(InconsistentInputError(
                    "Could not do a restart with appending because the checkpoint file "
                    "was not found. Either supply the name of the right checkpoint file "
                    "or do not use -append"));
        }
        // No need to tell the user that mdrun -cpi without a file means a new simulation
        handler.startingBehavior = StartingBehavior::NewSimulation;
        return handler;
    }

    t_fileio* fp = gmx_fio_open(checkpointFilename, "r");
    if (fp == nullptr)
    {
        GMX_THROW(FileIOError(
                formatString("Checkpoint file '%s' was found but could not be opened for "
                             "reading. Check the file permissions.",
                             checkpointFilename)));
    }

    std::vector<gmx_file_position_t> outputFiles;
    CheckpointHeaderContents         headerContents =
            read_checkpoint_simulation_part_and_filenames(fp, &outputFiles);

    GMX_RELEASE_ASSERT(!outputFiles.empty(),
                       "The checkpoint file or its reading is broken, as no output "
                       "file information is stored in it");
    const auto& logFilename = outputFiles[0].filename;
    GMX_RELEASE_ASSERT(fn2ftp(logFilename) == efLOG,
                       formatString("The checkpoint file or its reading is broken, the first "
                                    "output file '%s' must be a log file with extension '%s'",
                                    logFilename,
                                    ftp2ext(efLOG))
                               .c_str());

    if (appendingBehavior != AppendingBehavior::NoAppending)
    {
        // See whether appending can be done.

        size_t numFilesMissing = std::count_if(
                std::begin(outputFiles), std::end(outputFiles), [nfile, fnm](const auto& outputFile) {
                    return !exist_output_file(outputFile.filename, nfile, fnm);
                });
        if (numFilesMissing != 0)
        {
            // Appending is not possible, because not all previous
            // output files are present. We don't automatically switch
            // to numbered output files either, because that prevents
            // the user from using appending in future. If they want
            // to restart with missing files, they need to use
            // -noappend.
            throwBecauseOfMissingOutputFiles(checkpointFilename, outputFiles, nfile, fnm, numFilesMissing);
        }

        for (const auto& outputFile : outputFiles)
        {
            if (outputFile.offset < 0)
            {
                // Appending of large files is not possible unless
                // mdrun and the filesystem can do a correct job. We
                // don't automatically switch to numbered output files
                // either, because the user can benefit from
                // understanding that their infrastructure is not very
                // suitable for running a simulation producing lots of
                // output.
                auto message = formatString(
                        "The original mdrun wrote a file called '%s' which "
                        "is larger than 2 GB, but that mdrun or the filesystem "
                        "it ran on (e.g FAT32) did not support such large files. "
                        "This simulation cannot be restarted with appending. It will "
                        "be easier for you to use mdrun on a 64-bit filesystem, but "
                        "if you choose not to, then you must run mdrun with "
                        "-noappend once your output gets large enough.",
                        outputFile.filename);
                GMX_THROW(InconsistentInputError(message));
            }
        }

        const char* logFilename = outputFiles[0].filename;
        // If the precision does not match, we cannot continue with
        // appending, and will switch to not appending unless
        // instructed otherwise.
        if (headerContents.file_version >= CheckPointVersion::DoublePrecisionBuild
            && headerContents.double_prec != GMX_DOUBLE)
        {
            if (appendingBehavior == AppendingBehavior::Appending)
            {
                GMX_THROW(InconsistentInputError(formatString(
                        "Cannot restart with appending because the previous simulation part used "
                        "%s precision which does not match the %s precision used by this build "
                        "of GROMACS. Either use matching precision or use mdrun -noappend.",
                        precisionToString(headerContents.double_prec),
                        precisionToString(GMX_DOUBLE))));
            }
        }
        // If the previous log filename had a part number, then we
        // cannot continue with appending, and will continue without
        // appending.
        else if (hasSuffixFromNoAppend(logFilename))
        {
            if (appendingBehavior == AppendingBehavior::Appending)
            {
                GMX_THROW(InconsistentInputError(
                        "Cannot restart with appending because the previous simulation "
                        "part did not use appending. Either do not use mdrun -append, or "
                        "provide the correct checkpoint file."));
            }
        }
        else
        {
            // Everything is perfect - we can do an appending restart.
            handler = { StartingBehavior::RestartWithAppending, headerContents, outputFiles };
            return handler;
        }

        // No need to tell the user anything because the previous
        // simulation part also didn't append and that can only happen
        // when they ask for it.
    }

    GMX_RELEASE_ASSERT(appendingBehavior != AppendingBehavior::Appending,
                       "Logic error in appending");
    handler = { StartingBehavior::RestartWithoutAppending, headerContents, outputFiles };
    return handler;
}

//! Check whether chksum_file output file has a checksum that matches \c outputfile from the checkpoint.
void checkOutputFile(t_fileio* fileToCheck, const gmx_file_position_t& outputfile)
{
    /* compute md5 chksum */
    std::array<unsigned char, 16> digest;
    if (outputfile.checksumSize != -1)
    {
        if (gmx_fio_get_file_md5(fileToCheck, outputfile.offset, &digest)
            != outputfile.checksumSize) /*at the end of the call the file position is at the end of the file*/
        {
            auto message = formatString(
                    "Can't read %d bytes of '%s' to compute checksum. The file "
                    "has been replaced or its contents have been modified. Cannot "
                    "do appending because of this condition.",
                    outputfile.checksumSize,
                    outputfile.filename);
            GMX_THROW(InconsistentInputError(message));
        }
    }

    /* compare md5 chksum */
    if (outputfile.checksumSize != -1 && digest != outputfile.checksum)
    {
        if (debug)
        {
            fprintf(debug, "chksum for %s: ", outputfile.filename);
            for (int j = 0; j < 16; j++)
            {
                fprintf(debug, "%02x", digest[j]);
            }
            fprintf(debug, "\n");
        }
        auto message = formatString(
                "Checksum wrong for '%s'. The file has been replaced "
                "or its contents have been modified. Cannot do appending "
                "because of this condition.",
                outputfile.filename);
        GMX_THROW(InconsistentInputError(message));
    }
}

/*! \brief If supported, obtain a write lock on the log file.
 *
 * This wil prevent e.g. other mdrun instances from changing it while
 * we attempt to restart with appending. */
void lockLogFile(t_fileio* logfio, const std::filesystem::path& logFilename)
{
    /* Note that there are systems where the lock operation
     * will succeed, but a second process can also lock the file.
     * We should probably try to detect this.
     */
#if defined __native_client__
    errno = ENOSYS;
    if (true)
#elif GMX_NATIVE_WINDOWS
    if (_locking(fileno(gmx_fio_getfp(logfio)), _LK_NBLCK, LONG_MAX) == -1)
#else
    // don't initialize here: the struct order is OS dependent!
    struct flock fl;
    fl.l_type   = F_WRLCK;
    fl.l_whence = SEEK_SET;
    fl.l_start  = 0;
    fl.l_len    = 0;
    fl.l_pid    = 0;

    if (fcntl(fileno(gmx_fio_getfp(logfio)), F_SETLK, &fl) == -1)
#endif
    {
        if (errno == ENOSYS)
        {
            std::string message =
                    "File locking is not supported on this system. "
                    "Use mdrun -noappend to restart.";
            GMX_THROW(FileIOError(message));
        }
        else if (errno == EACCES || errno == EAGAIN)
        {
            auto message = formatString(
                    "Failed to lock: %s. Already running "
                    "simulation?",
                    logFilename.string().c_str());
            GMX_THROW(FileIOError(message));
        }
        else
        {
            auto message = formatString(
                    "Failed to lock: %s. %s.", logFilename.string().c_str(), std::strerror(errno));
            GMX_THROW(FileIOError(message));
        }
    }
}

/*! \brief Prepare to append to output files.
 *
 * We use the file pointer positions of the output files stored in the
 * checkpoint file and truncate the files such that any frames written
 * after the checkpoint time are removed.  All files are md5sum
 * checked such that we can be sure that we do not truncate other
 * (maybe important) files. The log file is locked so that we can
 * avoid cases where another mdrun instance might still be writing to
 * the file. */
void prepareForAppending(const ArrayRef<const gmx_file_position_t> outputFiles, t_fileio* logfio)
{
    if (GMX_FAHCORE)
    {
        // Can't check or truncate output files in general
        // TODO do we do this elsewhere for GMX_FAHCORE?
        return;
    }

    // Handle the log file separately - it comes first in the list
    // because we have already opened the log file. This ensures that
    // we retain a lock on the open file that is never lifted after
    // the checksum is calculated.
    const gmx_file_position_t& logOutputFile = outputFiles[0];
    lockLogFile(logfio, logOutputFile.filename);
    checkOutputFile(logfio, logOutputFile);

    if (gmx_fio_seek(logfio, logOutputFile.offset) != 0)
    {
        auto message = formatString("Seek error! Failed to truncate log file: %s.", std::strerror(errno));
        GMX_THROW(FileIOError(message));
    }

    // Now handle the remaining outputFiles
    for (const auto& outputFile : outputFiles.subArray(1, outputFiles.size() - 1))
    {
        t_fileio* fileToCheck = gmx_fio_open(outputFile.filename, "r+");
        checkOutputFile(fileToCheck, outputFile);
        gmx_fio_close(fileToCheck);

        if (gmx_truncate(outputFile.filename, outputFile.offset) != 0)
        {
            auto message = formatString(
                    "Truncation of file %s failed. Cannot do appending "
                    "because of this failure.",
                    outputFile.filename);
            GMX_THROW(FileIOError(message));
        }
    }
}

void StartingBehaviorHandler::ensureMultiSimBehaviorsMatch(const gmx_multisim_t* ms)
{
    if (!isMultiSim(ms))
    {
        // Trivially, the multi-sim behaviors match
        return;
    }

    auto startingBehaviors = gatherIntFromMultiSimulation(ms, static_cast<int>(startingBehavior));
    bool identicalStartingBehaviors =
            (std::adjacent_find(
                     std::begin(startingBehaviors), std::end(startingBehaviors), std::not_equal_to<>())
             == std::end(startingBehaviors));

    const EnumerationArray<StartingBehavior, std::string> behaviorStrings = {
        { "restart with appending", "restart without appending", "new simulation" }
    };
    if (!identicalStartingBehaviors)
    {
        std::string message = formatString(R"(
Multi-simulations must all start in the same way, either a new
simulation, a restart with appending, or a restart without appending.
However, the contents of the multi-simulation directories you specified
were inconsistent with each other. Either remove the checkpoint file
from each directory, or ensure each directory has a checkpoint file from
the same simulation part (and, if you want to append to output files,
ensure the old output files are present and named as they were when the
checkpoint file was written).

To help you identify which directories need attention, the %d
simulations wanted the following respective behaviors:
)",
                                           ms->numSimulations_);
        for (Index simIndex = 0; simIndex != gmx::ssize(startingBehaviors); ++simIndex)
        {
            auto behavior = static_cast<StartingBehavior>(startingBehaviors[simIndex]);
            message += formatString(
                    "  Simulation %6zd: %s\n", simIndex, behaviorStrings[behavior].c_str());
        }
        GMX_THROW(InconsistentInputError(message));
    }

    if (startingBehavior == StartingBehavior::NewSimulation)
    {
        // When not restarting, the behaviors are now known to match
        return;
    }

    // Multi-simulation restarts require that each checkpoint
    // describes the same simulation part. If those don't match, then
    // the simulation cannot proceed.
    auto simulationParts = gatherIntFromMultiSimulation(ms, headerContents->simulation_part);
    bool identicalSimulationParts =
            (std::adjacent_find(
                     std::begin(simulationParts), std::end(simulationParts), std::not_equal_to<>())
             == std::end(simulationParts));

    if (!identicalSimulationParts)
    {
        std::string message = formatString(R"(
Multi-simulations must all start in the same way, either a new
simulation, a restart with appending, or a restart without appending.
However, the checkpoint files you specified were from different
simulation parts. Either remove the checkpoint file from each directory,
or ensure each directory has a checkpoint file from the same simulation
part (and, if you want to append to output files, ensure the old output
files are present and named as they were when the checkpoint file was
written).

To help you identify which directories need attention, the %d
simulation checkpoint files were from the following respective
simulation parts:
)",
                                           ms->numSimulations_);
        for (Index partIndex = 0; partIndex != gmx::ssize(simulationParts); ++partIndex)
        {
            message += formatString("  Simulation %6zd: %d\n", partIndex, simulationParts[partIndex]);
        }
        GMX_THROW(InconsistentInputError(message));
    }
}

std::optional<int> StartingBehaviorHandler::makeIndexOfNextPart(const AppendingBehavior appendingBehavior) const
{
    std::optional<int> indexOfNextPart;

    if (startingBehavior == StartingBehavior::RestartWithoutAppending)
    {
        indexOfNextPart = headerContents->simulation_part + 1;
    }
    else if (startingBehavior == StartingBehavior::NewSimulation
             && appendingBehavior == AppendingBehavior::NoAppending)
    {
        indexOfNextPart = 1;
    }

    return indexOfNextPart;
}

} // namespace

std::tuple<StartingBehavior, LogFilePtr> handleRestart(const bool              isSimulationMain,
                                                       MPI_Comm                communicator,
                                                       const gmx_multisim_t*   ms,
                                                       const AppendingBehavior appendingBehavior,
                                                       const int               nfile,
                                                       t_filenm                fnm[])
{
    StartingBehaviorHandler handler;
    LogFilePtr              logFileGuard = nullptr;

    // Make sure all ranks agree on whether the (multi-)simulation can
    // proceed.
    int                numErrorsFound = 0;
    std::exception_ptr exceptionPtr;

    // Only the main rank of each simulation can do anything with
    // output files, so it is the only one that needs to consider
    // whether a restart might take place, and how to implement it.
    if (isSimulationMain)
    {
        try
        {
            handler = chooseStartingBehavior(appendingBehavior, nfile, fnm);

            handler.ensureMultiSimBehaviorsMatch(ms);

            // When not appending, prepare a suffix for the part number
            std::optional<int> indexOfNextPart = handler.makeIndexOfNextPart(appendingBehavior);

            // If a part suffix is used, change the file names accordingly.
            if (indexOfNextPart)
            {
                std::string suffix = formatString(".part%04d", *indexOfNextPart);
                add_suffix_to_output_names(fnm, nfile, suffix.c_str());
            }

            // Open the log file, now that it has the right name
            logFileGuard = openLogFile(ftp2fn(efLOG, nfile, fnm),
                                       handler.startingBehavior == StartingBehavior::RestartWithAppending);

            // When appending, the other output files need special handling before opening
            if (handler.startingBehavior == StartingBehavior::RestartWithAppending)
            {
                prepareForAppending(*handler.outputFiles, logFileGuard.get());
            }
        }
        catch (const std::exception& /*ex*/)
        {
            exceptionPtr   = std::current_exception();
            numErrorsFound = 1;
        }
    }
    // Since the main rank (perhaps of only one simulation) may have
    // found an error condition, we now coordinate the behavior across
    // all ranks. However, only the applicable ranks will throw a
    // non-default exception.
    //
    // TODO Evolve some re-usable infrastructure for this, because it
    // will be needed in many places while setting up simulations.
#if GMX_LIB_MPI
    int reducedNumErrorsFound;
    MPI_Allreduce(&numErrorsFound, &reducedNumErrorsFound, 1, MPI_INT, MPI_SUM, communicator);
    numErrorsFound = reducedNumErrorsFound;
#else
    // There is nothing to do with no MPI or thread-MPI, as there is
    // only one rank at this point.
    GMX_RELEASE_ASSERT(communicator == MPI_COMM_NULL, "Must have null communicator at this point");
#endif

    // Throw in a globally coordinated way, if needed
    if (numErrorsFound > 0)
    {
        if (exceptionPtr)
        {
            std::rethrow_exception(exceptionPtr);
        }
        else
        {
            GMX_THROW(ParallelConsistencyError("Another MPI rank encountered an exception"));
        }
    }

    // Ensure all ranks agree on the starting behavior, which is easy
    // because all simulations in a multi-simulation already agreed on
    // the starting behavior. There is nothing to do with no
    // MPI or thread-MPI.
#if GMX_LIB_MPI
    MPI_Bcast(&handler.startingBehavior, 1, MPI_INT, 0, communicator);
#endif

    return std::make_tuple(handler.startingBehavior, std::move(logFileGuard));
}

} // namespace gmx
