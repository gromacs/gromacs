/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Defines the CheckpointHelper class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "checkpointhelper.h"

#include "buildinfo.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/edsamhistory.h"        // TODO: Only needed for readLegacy / writeLegacy
#include "gromacs/mdtypes/energyhistory.h"       // TODO: Only needed for readLegacy / writeLegacy
#include "gromacs/mdtypes/inputrec.h"            // TODO: Only needed for readLegacy / writeLegacy
#include "gromacs/mdtypes/observableshistory.h"  // TODO: Only needed for readLegacy / writeLegacy
#include "gromacs/mdtypes/state.h"               // TODO: Only needed for readLegacy / writeLegacy
#include "gromacs/mdtypes/swaphistory.h"         // TODO: Only needed for readLegacy / writeLegacy
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/programcontext.h"

namespace gmx
{
const int CheckpointHelper::CPT_MAGIC1 = 171817;
const int CheckpointHelper::CPT_MAGIC2 = 171819;

CheckpointHelper::CheckpointHelper(
        const std::string &checkpointFilename) :
    cptFileVersion_(-1),
    doublePrecision_(-1),
    nnodes_(-1),
    npme_(-1),
    simulationPart_(-1),
    step_(-1),
    time_(-1),
    checkpointFilename_(checkpointFilename)
{}

CheckpointHelper::CheckpointHelper(
        const std::string               &gmxVersion,
        bool                             isDouble,
        const std::string               &programPath,
        const std::string               &generationTime,
        const std::string               &buildTime,
        const std::string               &buildUser,
        const std::string               &buildHost,
        int                              nnodes,
        int                              npmeNodes,
        const ivec                       domdecCells,
        int                              simulationPart,
        int64_t                          step,
        double                           time,
        std::vector<gmx_file_position_t> outputfiles,
        const std::string               &checkpointFilename) :
    cptFileVersion_(-1),
    doublePrecision_(static_cast<int>(isDouble)),
    nnodes_(nnodes),
    npme_(npmeNodes),
    simulationPart_(simulationPart),
    step_(step),
    time_(time),
    outputfiles_(std::move(outputfiles)),
    checkpointFilename_(checkpointFilename)
{
    std::strcpy(gmxVersion_, gmxVersion.c_str());
    std::strcpy(currentBinary_, programPath.c_str());
    std::strcpy(currentTime_, generationTime.c_str());
    std::strcpy(buildTime_, buildTime.c_str());
    std::strcpy(buildUser_, buildUser.c_str());
    std::strcpy(buildHost_, buildHost.c_str());
    copy_ivec(domdecCells, domdecNumCells_);
}

CheckpointVersion CheckpointHelper::getCheckpointVersion(const std::string &checkpointFilename)
{
    // Create helper & read header
    CheckpointHelper checkpointHelper = CheckpointHelper(checkpointFilename);
    auto             fp               = checkpointHelper.openReadCheckpointFile();
    checkpointHelper.doHeader(compat::make_unique<XDRWrapper>(fp), false);
    checkpointHelper.closeReadCheckpointFile(fp);

    if (checkpointHelper.cptFileVersion_ <= static_cast<int>(CheckpointVersion::legacy))
    {
        return CheckpointVersion::legacy;
    }
    if (checkpointHelper.cptFileVersion_ > static_cast<int>(CheckpointVersion::current))
    {
        gmx_fatal(FARGS, "Attempting to read a checkpoint file of version %d with code of version %d\n",
                  checkpointHelper.cptFileVersion_, static_cast<int>(CheckpointVersion::current));
    }
    return static_cast<CheckpointVersion>(checkpointHelper.cptFileVersion_);
}

void CheckpointHelper::doHeader(std::unique_ptr<XDRWrapper> xd, bool bRead)
{
    bool_t res = 0;
    int    magic;

    if (bRead)
    {
        magic = -1;
    }
    else
    {
        magic = CPT_MAGIC1;
    }
    res = xdr_int(xd->getXDR(), &magic);
    if (res == 0)
    {
        gmx_fatal(FARGS, "The checkpoint file is empty/corrupted, or maybe you are out of disk space?");
    }
    if (magic != CPT_MAGIC1)
    {
        gmx_fatal(FARGS, "Start of file magic number mismatch, checkpoint file has %d, should be %d\n"
                  "The checkpoint file is corrupted or not a checkpoint file",
                  magic, CPT_MAGIC1);
    }

    xd->doCptStringErr(gmxVersion_);
    xd->doCptStringErr(currentBinary_);
    xd->doCptStringErr(currentTime_);
    cptFileVersion_ = static_cast<int>(CheckpointVersion::current);
    xd->doCptIntErr(&cptFileVersion_);
    if (cptFileVersion_ > static_cast<int>(CheckpointVersion::current))
    {
        gmx_fatal(FARGS, "Attempting to read a checkpoint file of version %d with code of version %d\n",
                  cptFileVersion_, static_cast<int>(CheckpointVersion::current));
    }
    xd->doCptStringErr(buildTime_);
    xd->doCptStringErr(buildUser_);
    xd->doCptStringErr(buildHost_);
    xd->doCptIntErr(&doublePrecision_);

    xd->doCptIntErr(&nnodes_);
    xd->doCptIntErr(&npme_);
    xd->doCptIntErr(&domdecNumCells_[XX]);
    xd->doCptIntErr(&domdecNumCells_[YY]);
    xd->doCptIntErr(&domdecNumCells_[ZZ]);

    xd->doCptIntErr(&simulationPart_);
    xd->doCptInt64Err(&step_);
    xd->doCptDoubleErr(&time_);
}

void CheckpointHelper::doFooter(std::unique_ptr<XDRWrapper> xd, gmx_bool bRead)
{
    bool_t res = 0;
    int    magic;

    if (bRead)
    {
        magic = -1;
    }
    else
    {
        magic = CPT_MAGIC2;
    }
    res = xdr_int(xd->getXDR(), &magic);
    if (res == 0)
    {
        gmx_fatal(FARGS, "The checkpoint file is empty/corrupted, or maybe you are out of disk space?");
    }
    if (magic != CPT_MAGIC2)
    {
        gmx_fatal(FARGS, "End of file magic number mismatch, checkpoint file has %d, should be %d\n"
                  "The checkpoint file is corrupted or not a checkpoint file",
                  magic, CPT_MAGIC2);
    }
}

t_fileio* CheckpointHelper::openWriteCheckpointFile(int64_t step)
{
    std::string temporaryFilename = checkpointFilename_;
    /* if we can't rename, we just overwrite the cpt file.
     * dangerous if interrupted.
     */
    if (!GMX_NO_RENAME)
    {
        /* make the new temporary filename */
        std::string suffix      = "_" + std::to_string(step);
        size_t      dotPosition = checkpointFilename_.size() -
            std::strlen(ftp2ext(fn2ftp(checkpointFilename_.c_str()))) - 1;
        temporaryFilename.insert(dotPosition, suffix);
    }
    return gmx_fio_open(temporaryFilename.c_str(), "w");
}

void CheckpointHelper::closeWriteCheckpointFile(t_fileio* fp, bool bNumberAndKeep)
{
    /* we really, REALLY, want to make sure to physically write the checkpoint,
       and all the files it depends on, out to disk. Because we've
       opened the checkpoint with gmx_fio_open(), it's in our list
       of open files.  */
    t_fileio* ret = gmx_fio_all_output_fsync();

    if (ret)
    {
        char buf[STRLEN];
        sprintf(buf,
                "Cannot fsync '%s'; maybe you are out of disk space?",
                gmx_fio_getname(ret));

        if (getenv(GMX_IGNORE_FSYNC_FAILURE_ENV) == nullptr)
        {
            gmx_file(buf);
        }
        else
        {
            gmx_warning("%s", buf);
        }
    }

    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    /* we don't move the checkpoint if the user specified they didn't want it,
       or if the fsyncs failed */
    if (!GMX_NO_RENAME)
    {
        auto temporaryFilename = std::string(gmx_fio_getname(fp));
        if (!bNumberAndKeep && !ret)
        {
            if (gmx_fexist(checkpointFilename_.c_str()))
            {
                /* Rename the previous checkpoint file */
                std::string suffix      = "_prev";
                size_t      dotPosition = checkpointFilename_.size() -
                    std::strlen(ftp2ext(fn2ftp(checkpointFilename_.c_str()))) - 1;
                auto        oldCheckpointFilename = checkpointFilename_;
                oldCheckpointFilename.insert(dotPosition, suffix);
#ifndef GMX_FAHCORE
                {
                    /* we copy here so that if something goes wrong between now and
                     * the rename below, there's always a state.cpt.
                     * If renames are atomic (such as in POSIX systems),
                     * this copying should be unneccesary.
                     */
                    gmx_file_copy(checkpointFilename_.c_str(), oldCheckpointFilename.c_str(), FALSE);
                    /* We don't really care if this fails:
                     * there's already a new checkpoint.
                     */
                }
#else
                {
                    gmx_file_rename(checkpointFilename.c_str(), oldCheckpointFilename.c_str());
                }
#endif
            }
            if (gmx_file_rename(temporaryFilename.c_str(), checkpointFilename_.c_str()) != 0)
            {
                gmx_file("Cannot rename checkpoint file; maybe you are out of disk space?");
            }
        }
    }
}

t_fileio* CheckpointHelper::openReadCheckpointFile()
{
    return gmx_fio_open(checkpointFilename_.c_str(), "r");
}

void CheckpointHelper::closeReadCheckpointFile(t_fileio* fp)
{
    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }
}

int CheckpointHelper::doFiles(
        std::unique_ptr<XDRWrapper> xd, gmx_bool bRead)
{
    gmx_off_t                   offset;
    gmx_off_t                   mask = 0xFFFFFFFFL;
    int                         offset_high, offset_low;
    std::array<char, CPTSTRLEN> buf;

    // Ensure that reading pre-allocates outputfiles, while writing
    // writes what is already there.
    int nfiles = outputfiles_.size();
    if (xd->doCptInt(&nfiles) != 0)
    {
        return -1;
    }
    if (bRead)
    {
        outputfiles_.resize(nfiles);
    }

    for (auto &outputfile : outputfiles_)
    {
        /* 64-bit XDR numbers are not portable, so it is stored as separate high/low fractions */
        if (bRead)
        {
            xd->doCptStringErr(buf);
            std::strncpy(outputfile.filename, buf.data(), buf.size()-1);

            if (xd->doCptInt(&offset_high) != 0)
            {
                return -1;
            }
            if (xd->doCptInt(&offset_low) != 0)
            {
                return -1;
            }
            outputfile.offset = (static_cast<gmx_off_t>(offset_high) << 32 ) | ( static_cast<gmx_off_t>(offset_low) & mask );
        }
        else
        {
            xd->doCptStringErr(outputfile.filename);
            /* writing */
            offset      = outputfile.offset;
            if (offset == -1)
            {
                offset_low  = -1;
                offset_high = -1;
            }
            else
            {
                offset_low  = static_cast<int>(offset & mask);
                offset_high = static_cast<int>((offset >> 32) & mask);
            }
            if (xd->doCptInt(&offset_high) != 0)
            {
                return -1;
            }
            if (xd->doCptInt(&offset_low) != 0)
            {
                return -1;
            }
        }
        if (cptFileVersion_ >= 8)
        {
            if (xd->doCptInt(&outputfile.chksum_size) != 0)
            {
                return -1;
            }
            if (xd->doCptUChars(16, outputfile.chksum) != 0)
            {
                return -1;
            }
        }
        else
        {
            outputfile.chksum_size = -1;
        }
    }
    return 0;
}

void CheckpointHelper::printHeader(FILE *fplog)
{
    fprintf(fplog, "\n");
    fprintf(fplog, "Reading checkpoint file %s\n", checkpointFilename_.c_str());
    fprintf(fplog, "  file generated by:     %s\n", currentBinary_);
    fprintf(fplog, "  file generated at:     %s\n", currentTime_);
    fprintf(fplog, "  GROMACS build time:    %s\n", buildTime_);
    fprintf(fplog, "  GROMACS build user:    %s\n", buildUser_);
    fprintf(fplog, "  GROMACS build host:    %s\n", buildHost_);
    fprintf(fplog, "  GROMACS double prec.:  %d\n", doublePrecision_);
    fprintf(fplog, "  simulation part #:     %d\n", simulationPart_);
    fprintf(fplog, "  step:                  %s\n", std::to_string(step_).c_str());
    fprintf(fplog, "  time:                  %f\n", time_);
    fprintf(fplog, "\n");
}

void CheckpointHelper::checkInt(FILE *fplog, const char *type, int p, int f, bool *mm)
{
    bool foundMismatch = (p != f);
    if (!foundMismatch)
    {
        return;
    }
    *mm = TRUE;
    if (fplog)
    {
        fprintf(fplog, "  %s mismatch,\n", type);
        fprintf(fplog, "    current program: %d\n", p);
        fprintf(fplog, "    checkpoint file: %d\n", f);
        fprintf(fplog, "\n");
    }
}

void CheckpointHelper::checkString(
        FILE *fplog, const char *type, const char *p,
        const char *f, bool *mm)
{
    bool foundMismatch = (std::strcmp(p, f) != 0);
    if (!foundMismatch)
    {
        return;
    }
    *mm = TRUE;
    if (fplog)
    {
        fprintf(fplog, "  %s mismatch,\n", type);
        fprintf(fplog, "    current program: %s\n", p);
        fprintf(fplog, "    checkpoint file: %s\n", f);
        fprintf(fplog, "\n");
    }
}

void CheckpointHelper::checkMatch(
        FILE* fplog, bool reproducibilityRequested,
        const t_commrec* cr, const DomdecOptions &domdecOptions)
{
    /* Note that this checkString on the version will also print a message
     * when only the minor version differs. But we only print a warning
     * message further down with reproducibilityRequested=TRUE.
     */
    bool versionDiffers = false;
    checkString(fplog, "Version", gmx_version(), gmxVersion_, &versionDiffers);

    bool precisionDiffers = FALSE;
    checkInt(fplog, "Double prec.", GMX_DOUBLE, doublePrecision_, &precisionDiffers);
    if (precisionDiffers)
    {
        const char msg_precision_difference[] =
            "You are continuing a simulation with a different precision. Not matching\n"
            "single/double precision will lead to precision or performance loss.\n";
        if (fplog)
        {
            fprintf(fplog, "%s\n", msg_precision_difference);
        }
    }

    bool mm = (versionDiffers || precisionDiffers);

    if (reproducibilityRequested)
    {
        checkString(fplog, "Build time", BUILD_TIME, buildTime_, &mm);
        checkString(fplog, "Build user", BUILD_USER, buildUser_, &mm);
        checkString(fplog, "Build host", BUILD_HOST, buildHost_, &mm);
        checkString(fplog, "Program name", getProgramContext().fullBinaryPath(), currentBinary_, &mm);

        checkInt(fplog, "#ranks", cr->nnodes, nnodes_, &mm);
    }

    if (cr->nnodes > 1 && reproducibilityRequested)
    {
        checkInt(fplog, "#PME-ranks", cr->npmenodes, npme_, &mm);

        int npp = cr->nnodes;
        if (cr->npmenodes >= 0)
        {
            npp -= cr->npmenodes;
        }
        int npp_f = nnodes_;
        if (npme_ >= 0)
        {
            npp_f -= npme_;
        }
        if (npp == npp_f)
        {
            checkInt(fplog, "#DD-cells[x]", domdecOptions.numCells[XX], domdecNumCells_[XX], &mm);
            checkInt(fplog, "#DD-cells[y]", domdecOptions.numCells[YY], domdecNumCells_[YY], &mm);
            checkInt(fplog, "#DD-cells[z]", domdecOptions.numCells[ZZ], domdecNumCells_[ZZ], &mm);
        }
    }

    if (mm)
    {
        /* Gromacs should be able to continue from checkpoints between
         * different patch level versions, but we do not guarantee
         * compatibility between different major/minor versions - check this.
         */
        int        gmx_major;
        int        cpt_major;
        sscanf(gmx_version(), "%5d", &gmx_major);
        int        ret                 = sscanf(gmxVersion_, "%5d", &cpt_major);
        gmx_bool   majorVersionDiffers = (ret < 1 || gmx_major != cpt_major);

        const char msg_major_version_difference[] =
            "The current GROMACS major version is not identical to the one that\n"
            "generated the checkpoint file. In principle GROMACS does not support\n"
            "continuation from checkpoints between different versions, so we advise\n"
            "against this. If you still want to try your luck we recommend that you use\n"
            "the -noappend flag to keep your output files from the two versions separate.\n"
            "This might also work around errors where the output fields in the energy\n"
            "file have changed between the different versions.\n";

        const char msg_mismatch_notice[] =
            "GROMACS patchlevel, binary or parallel settings differ from previous run.\n"
            "Continuation is exact, but not guaranteed to be binary identical.\n";

        if (majorVersionDiffers)
        {
            if (fplog)
            {
                fprintf(fplog, "%s\n", msg_major_version_difference);
            }
        }
        else if (reproducibilityRequested)
        {
            /* Major & minor versions match at least, but something is different. */
            if (fplog)
            {
                fprintf(fplog, "%s\n", msg_mismatch_notice);
            }
        }
    }
}

void CheckpointHelper::appendOutputFiles(FILE *fplog gmx_unused, FILE **pfplog, bool bForceAppend gmx_unused)
{
    t_fileio     *chksum_file;
    unsigned char digest[16];
    if (outputfiles_.empty())
    {
        gmx_fatal(FARGS, "No names of output files were recorded in the checkpoint");
    }
    if (fn2ftp(outputfiles_[0].filename) != efLOG)
    {
        /* make sure first file is log file so that it is OK to use it for
         * locking
         */
        gmx_fatal(FARGS, "The first output file should always be the log "
                  "file but instead is: %s. Cannot do appending because of this condition.", outputfiles_[0].filename);
    }
    bool firstFile = true;
    for (const auto &outputfile : outputfiles_)
    {
        if (outputfile.offset < 0)
        {
            gmx_fatal(FARGS, "The original run wrote a file called '%s' which "
                      "is larger than 2 GB, but mdrun did not support large file"
                      " offsets. Can not append. Run mdrun with -noappend",
                      outputfile.filename);
        }
#ifndef GMX_FAHCORE
        {
            chksum_file = gmx_fio_open(outputfile.filename, "a");
        }
#else
        {
            chksum_file = gmx_fio_open(outputfile.filename, "r+");

            /* lock log file */
            if (firstFile)
            {
                /* Note that there are systems where the lock operation
                 * will succeed, but a second process can also lock the file.
                 * We should probably try to detect this.
                 */

#if !defined __native_client__ && !GMX_NATIVE_WINDOWS
                struct flock         fl; /* don't initialize here: the struct order is OS
                                            dependent! */
                fl.l_type   = F_WRLCK;
                fl.l_whence = SEEK_SET;
                fl.l_start  = 0;
                fl.l_len    = 0;
                fl.l_pid    = 0;
#endif
#if defined __native_client__
                errno = ENOSYS;
                if (1)

#elif GMX_NATIVE_WINDOWS
                if (_locking(fileno(gmx_fio_getfp(chksum_file)), _LK_NBLCK, LONG_MAX) == -1)
#else
                if (fcntl(fileno(gmx_fio_getfp(chksum_file)), F_SETLK, &fl) == -1)
#endif
                {
                    if (errno == ENOSYS)
                    {
                        if (!bForceAppend)
                        {
                            gmx_fatal(FARGS, "File locking is not supported on this system. Use -noappend or specify -append explicitly to append anyhow.");
                        }
                        else
                        {
                            if (fplog)
                            {
                                fprintf(fplog, "\nNOTE: File locking not supported on this system, will not lock %s\n\n", outputfile.filename);
                            }
                        }
                    }
                    else if (errno == EACCES || errno == EAGAIN)
                    {
                        gmx_fatal(FARGS, "Failed to lock: %s. Already running "
                                  "simulation?", outputfile.filename);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Failed to lock: %s. %s.",
                                  outputfile.filename, std::strerror(errno));
                    }
                }
            }

            /* compute md5 chksum */
            if (outputfile.chksum_size != -1)
            {
                if (gmx_fio_get_file_md5(chksum_file, outputfile.offset,
                                         digest) != outputfile.chksum_size) /*at the end of the call the file position is at the end of the file*/
                {
                    gmx_fatal(FARGS, "Can't read %d bytes of '%s' to compute checksum. The file has been replaced or its contents have been modified. Cannot do appending because of this condition.",
                              outputfile.chksum_size,
                              outputfile.filename);
                }
            }
            /* log file needs to be seeked in case we need to truncate
             * (other files are truncated below)*/
            if (firstFile)
            {
                if (gmx_fio_seek(chksum_file, outputfile.offset))
                {
                    gmx_fatal(FARGS, "Seek error! Failed to truncate log-file: %s.", std::strerror(errno));
                }
            }
        }
#endif

        /* open log file here - so that lock is never lifted
         * after chksum is calculated */
        if (firstFile)
        {
            *pfplog = gmx_fio_getfp(chksum_file);
        }
        else
        {
            gmx_fio_close(chksum_file);
        }
#ifndef GMX_FAHCORE
        /* compare md5 chksum */
        if (outputfile.chksum_size != -1 &&
            memcmp(digest, outputfile.chksum, 16) != 0)
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
            gmx_fatal(FARGS, "Checksum wrong for '%s'. The file has been replaced or its contents have been modified. Cannot do appending because of this condition.",
                      outputfile.filename);
        }
#endif

        // We could preprocess less, but then checkers complain
        // about possible bugs when using logic on constant
        // expressions.
#if !GMX_NATIVE_WINDOWS || !GMX_FAHCORE
        // The log file is already seeked to correct position, but
        // others must be truncated.
        if (!firstFile)
        {
            /* For FAHCORE, we do this elsewhere*/
            int rc = gmx_truncate(outputfile.filename, outputfile.offset);
            if (rc != 0)
            {
                gmx_fatal(FARGS, "Truncation of file %s failed. Cannot do appending because of this failure.", outputfile.filename);
            }
        }
#endif
        firstFile = false;
    }
}

void CheckpointHelper::readLegacy(
        std::unique_ptr<XDRWrapper> xd, const std::string &checkpointFilename, t_inputrec *ir,
        t_state *state, gmx_bool *bReadEkin, ObservablesHistory *observablesHistory)
{
    using legacy::cp_error;

    // Data that used to be in the header but doesn't really belong there in the new scheme

    int eIntegrator;
    int nlambda;

    int flags_eks;
    int flags_enh;
    int flags_dfh;
    int flags_awhh;

    int nED;
    int eSwapCoords;

    xd->doCptIntErr(&eIntegrator);

    xd->doCptIntErr(&nlambda);

    xd->doCptIntErr(&flags_eks);
    xd->doCptIntErr(&flags_enh);
    xd->doCptIntErr(&flags_dfh);
    xd->doCptIntErr(&flags_awhh);

    xd->doCptIntErr(&nED);
    xd->doCptIntErr(&eSwapCoords);

    int nlambdaHistory = (state->dfhist ? state->dfhist->nlambda : 0);
    if (nlambda != nlambdaHistory)
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system with %d lambda states, while the current system consists of %d lambda states", nlambda, nlambdaHistory);
    }

    if (eIntegrator != ir->eI)
    {
        gmx_fatal(FARGS, "Cannot change integrator during a checkpoint restart. Perhaps you should make a new .tpr with grompp -f new.mdp -t %s", checkpointFilename.c_str());
    }

    int ret;

    ret = legacy::do_cpt_ekinstate(xd->getXDR(), flags_eks, &state->ekinstate, nullptr);
    if (ret)
    {
        cp_error();
    }
    *bReadEkin = (((flags_eks & (1<<legacy::eeksEKINH)) != 0) ||
                  ((flags_eks & (1<<legacy::eeksEKINF)) != 0) ||
                  ((flags_eks & (1<<legacy::eeksEKINO)) != 0) ||
                  (((flags_eks & (1<<legacy::eeksEKINSCALEF)) |
                    (flags_eks & (1<<legacy::eeksEKINSCALEH)) |
                    (flags_eks & (1<<legacy::eeksVSCALE))) != 0));

    if (flags_enh && observablesHistory->energyHistory == nullptr)
    {
        observablesHistory->energyHistory = gmx::compat::make_unique<energyhistory_t>();
    }
    ret = legacy::do_cpt_enerhist(xd->getXDR(), TRUE,
                                  flags_enh, observablesHistory->energyHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = legacy::do_cpt_df_hist(xd->getXDR(), flags_dfh, nlambda, &state->dfhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    if (nED > 0 && observablesHistory->edsamHistory == nullptr)
    {
        observablesHistory->edsamHistory = gmx::compat::make_unique<edsamhistory_t>(edsamhistory_t {});
    }
    ret = legacy::do_cpt_EDstate(xd->getXDR(), TRUE, nED, observablesHistory->edsamHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    if (flags_awhh != 0 && state->awhHistory == nullptr)
    {
        state->awhHistory = std::make_shared<gmx::AwhHistory>();
    }
    ret = legacy::do_cpt_awh(xd->getXDR(), TRUE,
                             flags_awhh, state->awhHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    if (eSwapCoords != eswapNO && observablesHistory->swapHistory == nullptr)
    {
        observablesHistory->swapHistory = gmx::compat::make_unique<swaphistory_t>(swaphistory_t {});
    }
    ret = legacy::do_cpt_swapstate(xd->getXDR(), TRUE, eSwapCoords, observablesHistory->swapHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }
}

void CheckpointHelper::writeLegacy(
        std::unique_ptr<XDRWrapper> xd,
        bool bExpanded, int elamstats, int eIntegrator,
        t_state *state, ObservablesHistory *observablesHistory)
{
    int  flags_eks;
    if (state->ekinstate.bUpToDate)
    {
        flags_eks =
            ((1<<legacy::eeksEKIN_N) | (1<<legacy::eeksEKINH) | (1<<legacy::eeksEKINF) |
             (1<<legacy::eeksEKINO) | (1<<legacy::eeksEKINSCALEF) | (1<<legacy::eeksEKINSCALEH) |
             (1<<legacy::eeksVSCALE) | (1<<legacy::eeksDEKINDL) | (1<<legacy::eeksMVCOS));
    }
    else
    {
        flags_eks = 0;
    }

    energyhistory_t *enerhist  = observablesHistory->energyHistory.get();
    int              flags_enh = 0;
    if (enerhist != nullptr && (enerhist->nsum > 0 || enerhist->nsum_sim > 0))
    {
        flags_enh |= (1<<legacy::eenhENERGY_N) | (1<<legacy::eenhENERGY_NSTEPS) | (1<<legacy::eenhENERGY_NSTEPS_SIM);
        if (enerhist->nsum > 0)
        {
            flags_enh |= ((1<<legacy::eenhENERGY_AVER) | (1<<legacy::eenhENERGY_SUM) |
                          (1<<legacy::eenhENERGY_NSUM));
        }
        if (enerhist->nsum_sim > 0)
        {
            flags_enh |= ((1<<legacy::eenhENERGY_SUM_SIM) | (1<<legacy::eenhENERGY_NSUM_SIM));
        }
        if (enerhist->deltaHForeignLambdas != nullptr)
        {
            flags_enh |= ( (1<< legacy::eenhENERGY_DELTA_H_NN) |
                           (1<< legacy::eenhENERGY_DELTA_H_LIST) |
                           (1<< legacy::eenhENERGY_DELTA_H_STARTTIME) |
                           (1<< legacy::eenhENERGY_DELTA_H_STARTLAMBDA) );
        }
    }

    int flags_dfh;
    if (bExpanded)
    {
        flags_dfh = ((1<<legacy::edfhBEQUIL) | (1<<legacy::edfhNATLAMBDA) | (1<<legacy::edfhSUMWEIGHTS) |  (1<<legacy::edfhSUMDG)  |
                     (1<<legacy::edfhTIJ) | (1<<legacy::edfhTIJEMP));
        if (EWL(elamstats))
        {
            flags_dfh |= ((1<<legacy::edfhWLDELTA) | (1<<legacy::edfhWLHISTO));
        }
        if ((elamstats == elamstatsMINVAR) || (elamstats == elamstatsBARKER) || (elamstats == elamstatsMETROPOLIS))
        {
            flags_dfh |= ((1<<legacy::edfhACCUMP) | (1<<legacy::edfhACCUMM) | (1<<legacy::edfhACCUMP2) | (1<<legacy::edfhACCUMM2)
                          | (1<<legacy::edfhSUMMINVAR) | (1<<legacy::edfhSUMVAR));
        }
    }
    else
    {
        flags_dfh = 0;
    }

    int flags_awhh = 0;
    if (state->awhHistory != nullptr && !state->awhHistory->bias.empty())
    {
        flags_awhh |= ( (1 << legacy::eawhhIN_INITIAL) |
                        (1 << legacy::eawhhEQUILIBRATEHISTOGRAM) |
                        (1 << legacy::eawhhHISTSIZE) |
                        (1 << legacy::eawhhNPOINTS) |
                        (1 << legacy::eawhhCOORDPOINT) |
                        (1 << legacy::eawhhUMBRELLAGRIDPOINT) |
                        (1 << legacy::eawhhUPDATELIST) |
                        (1 << legacy::eawhhLOGSCALEDSAMPLEWEIGHT) |
                        (1 << legacy::eawhhNUMUPDATES) |
                        (1 << legacy::eawhhFORCECORRELATIONGRID));
    }

    /* We can check many more things now (CPU, acceleration, etc), but
     * it is highly unlikely to have two separate builds with exactly
     * the same version, user, time, and build host!
     */

    int                      nlambda     = (state->dfhist ? state->dfhist->nlambda : 0);

    edsamhistory_t          *edsamhist   = observablesHistory->edsamHistory.get();
    int                      nED         = (edsamhist ? edsamhist->nED : 0);

    swaphistory_t           *swaphist    = observablesHistory->swapHistory.get();
    int                      eSwapCoords = (swaphist ? swaphist->eSwapCoords : eswapNO);

    xd->doCptIntErr(&eIntegrator);
    xd->doCptIntErr(&nlambda);

    xd->doCptIntErr(&flags_eks);
    xd->doCptIntErr(&flags_enh);
    xd->doCptIntErr(&flags_dfh);
    xd->doCptIntErr(&flags_awhh);

    xd->doCptIntErr(&nED);
    xd->doCptIntErr(&eSwapCoords);

    if ((legacy::do_cpt_ekinstate(xd->getXDR(), flags_eks, &state->ekinstate, nullptr) < 0) ||
        (legacy::do_cpt_enerhist(xd->getXDR(), FALSE, flags_enh, enerhist, nullptr) < 0)  ||
        (legacy::do_cpt_df_hist(xd->getXDR(), flags_dfh, nlambda, &state->dfhist, nullptr) < 0)  ||
        (legacy::do_cpt_EDstate(xd->getXDR(), FALSE, nED, edsamhist, nullptr) < 0)      ||
        (legacy::do_cpt_awh(xd->getXDR(), FALSE, flags_awhh, state->awhHistory.get(), nullptr) < 0) ||
        (legacy::do_cpt_swapstate(xd->getXDR(), FALSE, eSwapCoords, swaphist, nullptr) < 0))
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }
}

}  // namespace gmx
