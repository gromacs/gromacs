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
 * Declares the CheckpointHelper class.
 *
 * This contains the header data of the checkpoint and helper functions to
 * read and write the checkpoints.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */
#ifndef GROMACS_CHECKPOINTHELPER_H
#define GROMACS_CHECKPOINTHELPER_H

#include "gromacs/fileio/xdrwrapper.h"

struct DomdecOptions;
struct ObservablesHistory;
struct t_commrec;
struct t_inputrec;
class t_state;

namespace gmx
{
/*! \brief Enum of values that describe the contents of a cpt file
 * whose format matches a version number
 *
 * The enum helps the code be more self-documenting and ensure merges
 * do not silently resolve when two patches make the same bump. When
 * adding new functionality, add a new element just above current
 * in this enumeration, and write code below that does the right thing
 * according to the value of file_version.
 */
enum class CheckpointVersion : int
{
    legacy = 19,       /**< The legacy version */
    current            /**< the total number of versions */
};

// TODO: Put that in the right spot - for now here to make our life easier
#define CPTSTRLEN 1024

class CheckpointHandler;

class CheckpointHelper
{
    public:
        void doHeader(std::unique_ptr<XDRWrapper> xd, bool bRead);
        void doFooter(std::unique_ptr<XDRWrapper> xd, gmx_bool bRead);
        int doFiles(std::unique_ptr<XDRWrapper> xd, gmx_bool bRead);

        void printHeader(FILE *fplog);
        void checkMatch(FILE* fplog, bool reproducibilityRequested,
                        const t_commrec* cr, const DomdecOptions &domdecOptions);
        void appendOutputFiles(FILE *fplog, FILE **pfplog, bool bForceAppend);

        t_fileio* openWriteCheckpointFile(int64_t step);
        void closeWriteCheckpointFile(t_fileio* fp, bool bNumberAndKeep);
        t_fileio* openReadCheckpointFile();
        void closeReadCheckpointFile(t_fileio* fp);

        static CheckpointVersion getCheckpointVersion(const std::string &checkpointFilename);

        friend class CheckpointHandler;

    private:
        CheckpointHelper(
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
            const std::string               &checkpointFilename);

        explicit CheckpointHelper(const std::string &checkpointFilename);

        static void checkInt(FILE *fplog, const char *type, int p, int f, bool *mm);
        static void checkString(FILE *fplog, const char *type, const char *p,
                                const char *f, bool *mm);

        static void readLegacy(
            std::unique_ptr<XDRWrapper> xd, const std::string &checkpointFilename, t_inputrec *ir,
            t_state *state, gmx_bool *bReadEkin, ObservablesHistory *observablesHistory);
        static void writeLegacy(
            std::unique_ptr<XDRWrapper> xd,
            bool bExpanded, int elamstats, int eIntegrator,
            t_state *state, ObservablesHistory *observablesHistory);

        const static int CPT_MAGIC1;
        const static int CPT_MAGIC2;

        /* Header information */
        int     cptFileVersion_;           // Checkpoint version
        char    gmxVersion_[CPTSTRLEN];    // GROMACS version
        char    buildTime_[CPTSTRLEN];     // Build time
        char    buildUser_[CPTSTRLEN];     // Build user
        char    buildHost_[CPTSTRLEN];     // Build host
        int     doublePrecision_;          // Double precision binary?
        char    currentBinary_[CPTSTRLEN]; // Generating program
        char    currentTime_[CPTSTRLEN];   // Generation time
        int     nnodes_;                   // Number of total nodes (cr)
        int     npme_;                     // Number of PME-only ranks (cr)
        ivec    domdecNumCells_;           // Domain decomposition vector (cr)
        int     simulationPart_;           // Simulation part (ir, gets incremented at checkpoint continuation)
        int64_t step_;
        double  time_;

        std::vector<gmx_file_position_t> outputfiles_;        // file positions (not part of header)
        const std::string                checkpointFilename_; // The filename (not part of header)
};

}      // namespace gmx

#endif //GROMACS_CHECKPOINTHELPER_H
