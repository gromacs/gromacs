/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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

#include <string>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/*! \brief All meta-data that is shared for one output file type for one bias */
struct OutFile
{
    char          baseFileName[STRLEN]; /**< Base output file name. */
    char          title[STRLEN];        /**< Title for the graph. */
    int           numDim;               /**< Number of dimensions. */
    int           firstGraphSubblock;   /**< Index in the energy sub-blocks for the first graph */
    int           numGraph;             /**< Number of graphs. */
    char        **legend;               /**< Legends for the output. */
    const char   *xLabel;               /**< Label for the x-axis. */
    const char   *yLabel;               /**< Label for the y-axis. */
};

/*! \brief All meta-data that is shared for all output files for one bias */
struct BiasReader
{
    int           subblockStart;        /**< The start index of the subblocks to read. */
    int           numSubblocks;         /**< Number of subblocks to read. */
    OutFile       awh;                  /**< The standard AWH output file data */
};

/*! \brief All options and meta-data needed for the AWH output */
struct AwhReader
{
    int           nbias;                /**< The number of AWH biases. */
    BiasReader   *biasReader;           /**< The readers, one for each AWH bias. */
};

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
static void make_legend(const awh_bias_params_t *awh_bias_params,
                        int awhGraphType,
                        int numGraph, char **legend)
{
    char               buf[256];

    const char        *legendBase[maxAwhGraphs] = {
        "PMF",
        "Coord bias",
        "Coord distr",
        "Ref value distr",
        "Target ref value distr"
    };

    int                numLegend   = awh_bias_params->ndim - 1 + numGraph;

    int                legendIndex = 0;
    for (int d = 1; d < awh_bias_params->ndim; d++)
    {
        /* Give legends to higher dimensions */
        sprintf(buf, "dim%d", d);
        legend[legendIndex] = strdup(buf);
        legendIndex++;
    }

    switch (awhGraphType)
    {
        case awhGraphTypeAwh:
        {
            /* Add as many legends as possible from the "base" legend list */
            int legendBaseIndex = 0;
            while (legendIndex < numLegend && legendBaseIndex < asize(legendBase))
            {
                legend[legendIndex] = strdup(legendBase[legendBaseIndex]);
                legendIndex++;
                legendBaseIndex++;
            }
        }
        break;
    }

    if (legendIndex != numLegend)
    {
        gmx_incons("Mismatch between the number of legends requested for printing and the number present!");
    }
}

/*! \brief Set the output base file name and title for one output file type */
static void setOutfileNameAndTitle(OutFile    *outFile,
                                   const char *name,
                                   const char *title,
                                   int         nbias,
                                   int         biasIndex)
{
    char buf[STRLEN];

    /* Store the input file name without extension in buf */
    sprintf(buf, "%s", name);
    char *ptr = strrchr(buf, '.');
    if (ptr != NULL)
    {
        *ptr = '\0';
    }

    if (nbias == 1)
    {
        sprintf(outFile->baseFileName, "%s", buf);
        sprintf(outFile->title, "%s", title);
    }
    else
    {
        sprintf(outFile->baseFileName, "%s%d", buf, biasIndex + 1);
        sprintf(outFile->title, "%s %d", title, biasIndex + 1);
    }
}

/*! \brief Initializes the output file setup for the AWH output (note that the filename is not set here). */
static void initOutfileSetupAwh(OutFile                 *outf,
                                const BiasReader        *biasReader,
                                const awh_bias_params_t *awh_bias_params,
                                gmx_bool                 moreGraphs)
{
    /* The first subblock with actual graph y-values is index 1 + ndim */
    int ndim                 = awh_bias_params->ndim;
    outf->numDim             = ndim;
    outf->firstGraphSubblock = biasReader->subblockStart + 1 + ndim;
    if (moreGraphs)
    {
        outf->numGraph       = std::min(biasReader->numSubblocks - 1 - ndim,
                                        maxAwhGraphs);
    }
    else
    {
        outf->numGraph       = 1;
    }

    snew(outf->legend, awh_bias_params->ndim - 1 + outf->numGraph);
    make_legend(awh_bias_params, awhGraphTypeAwh,
                outf->numGraph, outf->legend);
    outf->xLabel = "(nm or deg)";

    char buf[STRLEN];

    if (moreGraphs)
    {
        sprintf(buf, "(k\\sB\\NT), (nm\\S-%d\\N or rad\\S-%d\\N), (-)",
                awh_bias_params->ndim, awh_bias_params->ndim);
    }
    else
    {
        sprintf(buf, "(k\\sB\\NT)");
    }
    outf->yLabel = strdup(buf);
}

/*! \brief Initializes the reader, so we can start processing frames. */
static AwhReader *initAwhReader(const awh_params_t *awh_params,
                                int                 nfile,
                                const t_filenm      fnm[],
                                gmx_bool            moreGraphs,
                                const t_enxblock   *block)
{
    GMX_RELEASE_ASSERT(block != NULL, "NULL pointer passed to initAwhReader");

    AwhReader       *awhReader;

    snew(awhReader, 1);
    awhReader->nbias = awh_params->nbias;
    snew(awhReader->biasReader, awhReader->nbias);

    /* The first subblock of each AWH block has metadata has metadata about
     * the number of subblocks belonging to that AWH block.
     * (block->nsub gives only the total number of subblocks and not how
     * they are distributed between the AWH biases).
     */

    /* Keep track of the first subblock of this AWH */
    int subStart = 0;
    for (int k = 0; k < awhReader->nbias; k++)
    {
        awh_bias_params_t *awh_bias_params = &awh_params->awh_bias_params[k];
        BiasReader        *biasReader      = &awhReader->biasReader[k];

        biasReader->subblockStart = subStart;
        biasReader->numSubblocks  = (int)block->sub[subStart].fval[0];

        OutFile *outf = &biasReader->awh;
        setOutfileNameAndTitle(outf,
                               opt2fn("-o", nfile, fnm), "AWH",
                               awhReader->nbias, k);

        initOutfileSetupAwh(outf, biasReader, awh_bias_params, moreGraphs);

        subStart += biasReader->numSubblocks;
    }

    return awhReader;
}

/*! \brief Opens a single output file for a bias, prints title and legends */
static FILE *openBiasOutputFile(const OutFile          *outFile,
                                double                  t,
                                const gmx_output_env_t *oenv)
{
    char filename[STRLEN];

    sprintf(filename, "%s_t%g.xvg", outFile->baseFileName, t);

    FILE *fp = xvgropen(filename, outFile->title,
                        outFile->xLabel, outFile->yLabel, oenv);
    xvgr_legend(fp, outFile->numGraph, (const char**)outFile->legend, oenv);

    return fp;
}

/*! \brief Prints data selected by \p outf from \p block to \p fp */
static void printOutfileData(const t_enxblock       *block,
                             int                     subStart,
                             const OutFile          *outf,
                             FILE                   *fp)
{
    int numPoints = block->sub[subStart + 1].nr;
    for (int j = 0; j < numPoints; j++)
    {
        /* Print the coordinates for numDim dimensions */
        for (int d = 0; d < outf->numDim; d++)
        {
            fprintf(fp, "  %8.4f", block->sub[subStart + 1 + d].fval[j]);
        }

        /* Print numGraph observables */
        for (int i = outf->firstGraphSubblock; i < outf->firstGraphSubblock + outf->numGraph; i++)
        {
            fprintf(fp, "  %g", block->sub[i].fval[j]);
        }

        fprintf(fp, "\n");
    }
}

/*! \brief Extract and print AWH data for an AWH block of one energy frame */
static void processAwhFrame(const t_enxblock       *block,
                            double                  t,
                            const AwhReader        *awhReader,
                            const gmx_output_env_t *oenv)
{
    /* We look for AWH data every energy frame and count the no of AWH frames found. We only extract every 'skip' AWH frame. */

    /* Keep track of the starting index for each AWH */
    int subStart = 0;
    for (int k = 0; k < awhReader->nbias; k++)
    {
        /* Each frame and AWH instance extracted generates one xvg file. */
        {
            OutFile *outf = &awhReader->biasReader[k].awh;

            FILE    *fp   = openBiasOutputFile(outf, t, oenv);

            /* Now do the actual printing. Metadata in first subblock is treated separately. */
            fprintf(fp, "# AWH metadata: target error = %4.2f kT\n",
                    block->sub[subStart].fval[1]);

            fprintf(fp, "# AWH metadata: log sample weight = %4.2f\n",
                    block->sub[subStart].fval[2]);

            printOutfileData(block, subStart, outf, fp);

            gmx_ffclose(fp);
        }

        subStart += awhReader->biasReader[k].numSubblocks;
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
        "With [TT]-more[tt] the bias, and target and coordinate distributions",
        "are also printed.",
    };
    static gmx_bool     moreGraphs = FALSE;
    static int          skip       = 0;
    t_pargs             pa[]       = {
        { "-skip", FALSE, etINT,  {&skip},
          "Skip number of frames between data points" },
        { "-more", FALSE, etBOOL, {&moreGraphs},
          "Print more output" }
    };

    ener_file_t         fp;
    t_inputrec          ir;
    gmx_enxnm_t        *enm = NULL;
    t_enxframe         *frame;
    int                 nre;
    gmx_output_env_t   *oenv;

    t_filenm            fnm[] = {
        { efEDR, "-f",    NULL,           ffREAD  },
        { efTPR, "-s",    NULL,           ffREAD  },
        { efXVG, "-o",    "awh",          ffWRITE }
    };
    const int           nfile  = asize(fnm);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END,
                           nfile, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
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
    read_tpx(ftp2fn(efTPR, nfile, fnm), &ir, box, &natoms, NULL, NULL, &mtop);

    if (!ir.bDoAwh)
    {
        gmx_fatal(FARGS, "No AWH data in %s\n", opt2fn("-f", nfile, fnm));
    }

    AwhReader *awhReader = NULL;

    /* Initiate counters */
    gmx_bool   haveFrame;
    int        awhFrameCounter = 0;
    int        timeCheck       = 0;
    do
    {
        haveFrame = do_enx(fp, frame);

        bool        useFrame = false;

        t_enxblock *block = NULL;

        if (haveFrame)
        {
            timeCheck = check_times(frame->t);

            if (timeCheck == 0)
            {
                block = find_block_id_enxframe(frame, enxAWH, NULL);

                if (block != NULL)
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
            if (awhReader == NULL)
            {
                awhReader = initAwhReader(ir.awh_params,
                                          nfile, fnm, moreGraphs,
                                          block);
            }

            processAwhFrame(block, frame->t, awhReader, oenv);
        }
    }
    while (haveFrame && (timeCheck <= 0));

    fprintf(stderr, "\n");

    close_enx(fp);

    return 0;
}
