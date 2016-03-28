/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 *  \brief
 *  Tool for extracting AWH (Accelerated Weight Histogram) data from energy files.
 *
 *  \author Viveca Lindahl
 *  \author Berk Hess
 */

#include "gmxpre.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>
#include <string>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

using gmx::AwhParams;
using gmx::AwhBiasParams;

/*! \brief All meta-data that is shared for one output file type for one bias */
struct OutFile
{
    std::string              baseFileName;       /**< Base of the output file name. */
    std::string              title;              /**< Title for the graph. */
    int                      numDim;             /**< Number of dimensions. */
    int                      firstGraphSubblock; /**< Index in the energy sub-blocks for the first graph */
    int                      numGraph;           /**< Number of actual graphs. */
    bool                     kTUnit;             /**< Whether to use kT instead of kJ/mol */
    real                     kT;                 /**< k_B*T */

    std::vector<std::string> legend;             /**< Legends for the output */
    std::string              xLabel;             /**< Label for the x-axis. */
    std::string              yLabel;             /**< Label for the y-axis. */
};

/*! \brief All meta-data that is shared for all output files for one bias */
class BiasReader
{
    public:
        //! Constructor.
        BiasReader(int                      subblockStart,
                   int                      numSubblocks,
                   std::unique_ptr<OutFile> awh) :
            subblockStart_(subblockStart),
            numSubblocks_(numSubblocks),
            awh_(std::move(awh))
        {
        }

        //! Return the AWH output file data.
        const OutFile &awh() const
        {
            return *awh_.get();
        }

        //! Return the starting subblock.
        int subblockStart() const
        {
            return subblockStart_;
        }

    private:
        const int                subblockStart_;  /**< The start index of the subblocks to read. */
        const int                numSubblocks_;   /**< Number of subblocks to read. */
        std::unique_ptr<OutFile> awh_;            /**< The standard AWH output file data. */
        /* NOTE: A second OutFile will be added soon, this will also make numSubblocks_ useful. */
};

/*! \brief All options and meta-data needed for the AWH output */
class AwhReader
{
    public:
        //! Constructor
        AwhReader(const AwhParams  *awhParams,
                  int               nfile,
                  const t_filenm    fnm[],
                  bool              moreGraphs,
                  bool              kTUnit,
                  real              kT,
                  const t_enxblock *block);

        //! Extract and print AWH data for an AWH block of one energy frame
        void processAwhFrame(const t_enxblock       *block,
                             double                  t,
                             const gmx_output_env_t *oenv);

    private:
        std::vector<BiasReader> biasReader_; /**< The readers, one for each AWH bias. */
};

/* NOTE: A second value will be added soon. */
enum {
    awhGraphTypeAwh
};

/*! \brief The maximum number of observables per subblock, not including the full friction tensor.
 *
 * This number, as well as the legends in make_legend() should match
 * the output that mdrun writes. It would be better to define these
 * values in a single location.
 */
static const int maxAwhGraphs = 5;

/*! \brief Constructs a legend for a standard awh output file */
static std::vector<std::string>makeLegend(const AwhBiasParams *awhBiasParams,
                                          int                  awhGraphType,
                                          size_t               numLegend)
{
    const char              *legendBase[maxAwhGraphs] = {
        "PMF",
        "Coord bias",
        "Coord distr",
        "Ref value distr",
        "Target ref value distr"
    };

    std::vector<std::string> legend;
    for (int d = 1; d < awhBiasParams->ndim; d++)
    {
        /* Give legends to higher dimensions */
        legend.push_back(gmx::formatString("dim%d", d));
    }

    switch (awhGraphType)
    {
        case awhGraphTypeAwh:
        {
            /* Add as many legends as possible from the "base" legend list */
            int legendBaseIndex = 0;
            while (legend.size() < numLegend && legendBaseIndex < asize(legendBase))
            {
                legend.push_back(legendBase[legendBaseIndex]);
                legendBaseIndex++;
            }
        }
        break;
    }

    GMX_RELEASE_ASSERT(legend.size() == numLegend, "The number of legends requested for printing and the number generated should be the same");

    return legend;
}

/*! \brief Set the output base file name and title for one output file type */
static void setOutfileNameAndTitle(OutFile           *outFile,
                                   const std::string  name,
                                   const std::string  title,
                                   int                nbias,
                                   int                biasIndex)
{
    char buf[STRLEN];

    /* Store the input file name without extension in buf */
    outFile->baseFileName = name.substr(0, name.find('.'));

    if (nbias == 1)
    {
        outFile->baseFileName = buf;
        outFile->title        = title;
    }
    else
    {
        outFile->baseFileName = gmx::formatString("%s%d", buf, biasIndex + 1);
        outFile->title        = gmx::formatString("%s %d", title, biasIndex + 1);
    }
}

/*! \brief Initializes the output file setup for the AWH output (note that the filename is not set here). */
static void initOutfileSetupAwh(OutFile             *outf,
                                int                  subblockStart,
                                int                  numSubblocks,
                                const AwhBiasParams *awhBiasParams,
                                bool                 moreGraphs,
                                bool                 kTUnit,
                                real                 kT)
{
    /* The first subblock with actual graph y-values is index 1 + ndim */
    int ndim                 = awhBiasParams->ndim;
    outf->numDim             = ndim;
    outf->firstGraphSubblock = subblockStart + 1 + ndim;
    if (moreGraphs)
    {
        /* There are one metadata and ndim coordinate blocks */
        outf->numGraph       = std::min(numSubblocks - 1 - ndim,
                                        maxAwhGraphs);
    }
    else
    {
        outf->numGraph       = 1;
    }
    outf->kTUnit             = kTUnit;
    outf->kT                 = kT;
    int numLegend            = ndim - 1 + outf->numGraph;
    outf->legend = makeLegend(awhBiasParams, awhGraphTypeAwh, numLegend);
    outf->xLabel = "(nm or deg)";
    outf->yLabel = outf->kTUnit ? "(k\\sB\\NT)" : "(kJ/mol)";
    if (moreGraphs)
    {
        outf->yLabel += ", (nm\\S-%d\\N or rad\\S-%d\\N), (-)";
    }
}

AwhReader::AwhReader(const AwhParams  *awhParams,
                     int               nfile,
                     const t_filenm    fnm[],
                     bool              moreGraphs,
                     bool              kTUnit,
                     real              kT,
                     const t_enxblock *block)
{
    GMX_RELEASE_ASSERT(block != nullptr, "NULL pointer passed to initAwhReader");

    /* The first subblock of each AWH block has metadata about
     * the number of subblocks belonging to that AWH block.
     * (block->nsub gives only the total number of subblocks and not how
     * they are distributed between the AWH biases).
     */

    /* Keep track of the first subblock of this AWH */
    int subblockStart = 0;
    for (int k = 0; k < awhParams->numBias; k++)
    {
        AwhBiasParams *awhBiasParams = &awhParams->awhBiasParams[k];

        int            numSubblocks  = (int)block->sub[subblockStart].fval[0];

        std::unique_ptr<OutFile> awh(new OutFile());
        setOutfileNameAndTitle(awh.get(),
                               opt2fn("-o", nfile, fnm), "AWH",
                               awhParams->numBias, k);

        initOutfileSetupAwh(awh.get(), subblockStart, numSubblocks,
                            awhBiasParams, moreGraphs, kTUnit, kT);

        biasReader_.emplace_back(BiasReader(subblockStart, numSubblocks, std::move(awh)));

        subblockStart += numSubblocks;
    }
}

/*! \brief Opens a single output file for a bias, prints title and legends */
static FILE *openBiasOutputFile(const OutFile          &outFile,
                                double                  t,
                                const gmx_output_env_t *oenv)
{
    std::string filename = outFile.baseFileName + gmx::formatString("_t%g.xvg", t);

    FILE       *fp = xvgropen(filename.c_str(), outFile.title.c_str(),
                              outFile.xLabel, outFile.yLabel, oenv);
    xvgrLegend(fp, outFile.legend, oenv);

    return fp;
}

/*! \brief Prints data selected by \p outf from \p block to \p fp */
static void printOutfileData(const t_enxblock       *block,
                             int                     subStart,
                             const OutFile          &outf,
                             FILE                   *fp)
{
    int numPoints = block->sub[subStart + 1].nr;
    for (int j = 0; j < numPoints; j++)
    {
        /* Print the coordinates for numDim dimensions */
        for (int d = 0; d < outf.numDim; d++)
        {
            fprintf(fp, "  %8.4f", block->sub[subStart + 1 + d].fval[j]);
        }

        /* Print numGraph observables */
        for (int i = outf.firstGraphSubblock; i < outf.firstGraphSubblock + outf.numGraph; i++)
        {
            real scaleFactor = (outf.kTUnit || i > outf.firstGraphSubblock + 1) ? 1 : outf.kT;
            fprintf(fp, "  %g", block->sub[i].fval[j]*scaleFactor);
        }

        fprintf(fp, "\n");
    }
}

/*! \brief Extract and print AWH data for an AWH block of one energy frame */
void AwhReader::processAwhFrame(const t_enxblock       *block,
                                double                  t,
                                const gmx_output_env_t *oenv)
{
    /* We look for AWH data every energy frame and count the no of AWH frames found. We only extract every 'skip' AWH frame. */

    for (auto &biasReader : biasReader_)
    {
        /* Each frame and AWH instance extracted generates one xvg file. */
        {
            const OutFile &outf     = biasReader.awh();

            const int      subStart = biasReader.subblockStart();

            FILE          *fp       = openBiasOutputFile(outf, t, oenv);

            /* Now do the actual printing. Metadata in first subblock is treated separately. */
            fprintf(fp, "# AWH metadata: target error = %4.2f %s\n",
                    block->sub[subStart].fval[1]*(outf.kTUnit ? 1 : outf.kT),
                    outf.kTUnit ? "kT" : "kJ/mol");

            fprintf(fp, "# AWH metadata: log sample weight = %4.2f\n",
                    block->sub[subStart].fval[2]);

            printOutfileData(block, subStart, outf, fp);

            gmx_ffclose(fp);
        }
    }
}

/*! \brief The main function for the AWH tool */
int gmx_awh(int argc, char *argv[])
{
    const char         *desc[] = {
        "[THISMODULE] extracts AWH data from an energy file.",
        "One file is written per AWH bias per time frame.",
        "The bias index, if more than one, is appended to the file, as well as",
        "the time of the frame. By default only the PMF is printed.",
        "With [TT]-more[tt] the bias, target and coordinate distributions",
        "are also printed.",
    };
    static gmx_bool     moreGraphs = FALSE;
    static int          skip       = 0;
    static gmx_bool     kTUnit     = FALSE;
    t_pargs             pa[]       = {
        { "-skip", FALSE, etINT,  {&skip},
          "Skip number of frames between data points" },
        { "-more", FALSE, etBOOL, {&moreGraphs},
          "Print more output" },
        { "-kt", FALSE, etBOOL, {&kTUnit},
          "Print free energy output in units of kT instead of kJ/mol" }
    };

    ener_file_t         fp;
    t_inputrec          ir;
    gmx_enxnm_t        *enm = nullptr;
    t_enxframe         *frame;
    int                 nre;
    gmx_output_env_t   *oenv;

    t_filenm            fnm[] = {
        { efEDR, "-f", nullptr,           ffREAD  },
        { efTPR, "-s", nullptr,           ffREAD  },
        { efXVG, "-o",    "awh",          ffWRITE }
    };
    const int           nfile  = asize(fnm);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END,
                           nfile, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    snew(frame, 1);
    fp = open_enx(ftp2fn(efEDR, nfile, fnm), "r");
    do_enxnms(fp, &nre, &enm);

    /* We just need the AWH parameters from inputrec. These are used to initialize
       the AWH reader when we have a frame to read later on. */
    gmx_mtop_t mtop;
    int        natoms;
    matrix     box;
    read_tpx(ftp2fn(efTPR, nfile, fnm), &ir, box, &natoms, nullptr, nullptr, &mtop);

    if (!ir.bDoAwh)
    {
        gmx_fatal(FARGS, "No AWH data in %s\n", opt2fn("-f", nfile, fnm));
    }

    std::unique_ptr<AwhReader> awhReader;

    /* Initiate counters */
    gmx_bool   haveFrame;
    int        awhFrameCounter = 0;
    int        timeCheck       = 0;
    do
    {
        haveFrame = do_enx(fp, frame);

        bool        useFrame = false;

        t_enxblock *block    = nullptr;

        if (haveFrame)
        {
            timeCheck = check_times(frame->t);

            if (timeCheck == 0)
            {
                block = find_block_id_enxframe(frame, enxAWH, nullptr);

                if (block != nullptr)
                {
                    if ((skip == 0) || (awhFrameCounter % skip == 0))
                    {
                        useFrame = true;
                    }
                    awhFrameCounter++;
                }
            }
        }

        if (useFrame)
        {
            /* We read a valid frame with an AWH block, so we can use it */

            /* Currently we have the number of subblocks per AWH stored
             * in the first subblock (we cannot get this directly from the tpr),
             * so we have to read an AWH block before making the legends.
             */
            if (awhReader == nullptr)
            {
                awhReader =
                    std::unique_ptr<AwhReader>(new AwhReader(ir.awhParams,
                                                             nfile, fnm, moreGraphs,
                                                             kTUnit, BOLTZ*ir.opts.ref_t[0],
                                                             block));
            }

            awhReader->processAwhFrame(block, frame->t, oenv);
        }
    }
    while (haveFrame && (timeCheck <= 0));

    fprintf(stderr, "\n");

    close_enx(fp);

    return 0;
}
