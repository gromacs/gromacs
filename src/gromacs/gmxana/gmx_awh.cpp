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
#include "gromacs/utility/smalloc.h"

struct outfile_t {
    char          baseFileName[STRLEN];      /* Base output file name. */
    char          title[STRLEN];             /* Title for the graph. */
    int           startGraph;                /* Start index for the first graph */
    int           numGraph;                  /* Number of graphs. */
    char        **legend;                    /* Legends for the output. */
    const char   *ylabel;                    /* Label for the y-axis. */
};

struct bias_reader_t {
    int           nsubblocks;                /* Number of subblocks to read. */
    outfile_t     awh;                       /* The standard AWH output file data */
    outfile_t     tensor;                    /* The metric tensor output file data */
};

struct awh_reader_t {
    int               nbias;                  /* The number of AWH biases. */
    bias_reader_t    *bias_reader;            /* The readers, one for each AWH bias. */
    bool              outputMetricTensor;     /* Write metric tensor files */
};


static void make_awh_legend(const awh_bias_params_t *awh_bias_params,
                            int nleg, char **leg)
{
    int                i, d, ileg;
    char               buf[256];

    const char        *leg_base[] = {
        "PMF",
        "Coord bias",
        "Coord distr",
        "Ref value distr",
        "Target ref value distr"
    };

    ileg = 0;
    for (d = 1; d < awh_bias_params->ndim; d++)
    {
        /* Give legends to higher dimensions */
        sprintf(buf, "dim%d", d);
        leg[ileg] = strdup(buf);
        ileg++;
    }

    /* Add as many legends as possible from the "base" legend list */
    i = 0;
    while (ileg < nleg && i < asize(leg_base))
    {
        leg[ileg] = strdup(leg_base[i]);
        ileg++;
        i++;
    }

    /* If there are still more legends to add, do it. */
    if (ileg < nleg)
    {
        leg[ileg] =  strdup("Force corr distr.");
        ileg++;
    }

    if (ileg != nleg)
    {
        gmx_incons("Mismatch between the number of legends requested for printing and the number present!");
    }
}

static void make_tensor_legend(const awh_bias_params_t *awh_bias_params,
                               int nleg, char **leg)
{
    /* Want the 2 indices of the symmetric correlation matrix.
     * Only print half the elements */
    int ileg = 0;
    for (int i0 = 0; i0 < awh_bias_params->ndim; i0++)
    {
        for (int i1 = 0; i1 <= i0; i1++)
        {
            char buf[256];

            sprintf(buf, "%d,%d", i0, i1);

            leg[ileg] = strdup(buf);
            ileg++;
        }
    }

    if (ileg != nleg)
    {
        gmx_incons("Mismatch between the number of legends requested for printing and the number present!");
    }
}

static void set_outfile_name_and_title(outfile_t *outfile,
                                       const char *name,
                                       const char *title,
                                       int nbias, int biasIndex)
{
    char buf[STRLEN];

    /* Store the input file name without extension in buf buf */
    sprintf(buf, name);
    char *ptr = strrchr(buf, '.');
    if (ptr != NULL)
    {
        *ptr = '\0';
    }

    if (nbias == 1)
    {
        sprintf(outfile->baseFileName, "%s", buf);
        sprintf(outfile->title, "%s", title);
    }
    else
    {
        sprintf(outfile->baseFileName, "%s%d", buf, biasIndex + 1);
        sprintf(outfile->title, "%s %d", title, biasIndex + 1);
    }
}

static awh_reader_t *init_awh_reader(awh_params_t *awh_params,
                                     int nfile, const t_filenm fnm[],
                                     gmx_bool bMore,
                                     t_enxframe *fr)
{
    awh_reader_t       *awh_reader       = NULL;
    t_enxblock         *block            = NULL;

    /* The first subblock of each AWH has metadata has metadata about the number of subblocks belonging to that AWH
     * (block->nsub gives only the total number of subblocks and not how they are distributed between the AWHs). */
    block = find_block_id_enxframe(fr, enxAWH, NULL);
    if (block != NULL)
    {
        int imin_sub;

        snew(awh_reader, 1);
        awh_reader->nbias = awh_params->nbias;
        snew(awh_reader->bias_reader, awh_reader->nbias);

        awh_reader->outputMetricTensor = opt2bSet("-mt", nfile, fnm);

        /* Keep track of the first subblock of this AWH */
        imin_sub = 0;
        for (int k = 0; k < awh_reader->nbias; k++)
        {
            awh_bias_params_t *awh_bias_params = &awh_params->awh_bias_params[k];
            int                ndim            = awh_bias_params->ndim;
            int                numTensorElem   = (ndim*(ndim + 1))/2;

            bias_reader_t     *bias_reader     = &awh_reader->bias_reader[k];

            bias_reader->nsubblocks = (int)block->sub[imin_sub].dval[0];

            outfile_t *outf = &bias_reader->awh;
            set_outfile_name_and_title(outf,
                                       opt2fn("-o", nfile, fnm), "AWH",
                                       awh_reader->nbias, k);

            /* The metadata and the grid point values (for the first dimension) do not get legends (-2) */
            outf->startGraph = 2;
            outf->numGraph   = bMore ? bias_reader->nsubblocks - 2 - numTensorElem : 1;

            snew(outf->legend, outf->numGraph);
            make_awh_legend(awh_bias_params, outf->numGraph, outf->legend);
            outf->ylabel = bMore ? "(k\\sB\\NT), (nm\\S-1\\N or rad\\S-1\\N)" : "(k\\sB\\NT)";

            if (awh_reader->outputMetricTensor)
            {
                outfile_t *outf = &bias_reader->tensor;
                set_outfile_name_and_title(outf,
                                           opt2fn("-mt", nfile, fnm), "Metric tensor",
                                           awh_reader->nbias, k);

                outf->startGraph = bias_reader->nsubblocks - numTensorElem;
                outf->numGraph   = numTensorElem;
                snew(outf->legend, outf->numGraph);
                make_tensor_legend(awh_bias_params, outf->numGraph, outf->legend);
                outf->ylabel     = "(nm\\S2\\Nps or rad\\S2\\Nps)";
            }

            imin_sub               += bias_reader->nsubblocks;
        }
    }

    return awh_reader;
}

static FILE *open_awh_outfile(const outfile_t        *outfile,
                              double                  t,
                              const gmx_output_env_t *oenv)
{
    char filename[STRLEN];

    sprintf(filename, "%s_t%g.xvg", outfile->baseFileName, t);

    FILE *fp = xvgropen(filename, outfile->title,
                        "(nm or deg)", outfile->ylabel, oenv);
    xvgr_legend(fp, outfile->numGraph, (const char**)outfile->legend, oenv);

    return fp;
}

static void do_awh_frame(t_enxframe *fr, awh_reader_t *awh_reader, int *frame_counter, int skip, const gmx_output_env_t *oenv)
{
    t_enxblock *block = NULL;

    /* We look for AWH data every energy frame and count the no of AWH frames found. We only extract every 'skip' AWH frame. */
    block = find_block_id_enxframe(fr, enxAWH, NULL);

    if (block != NULL)
    {
        if (!skip || (*frame_counter) % skip == 0)
        {
            int    imin_sub;

            /* Keep track of the starting index for each AWH */
            imin_sub = 0;
            for (int k = 0; k < awh_reader->nbias; k++)
            {
                outfile_t *outf = &awh_reader->bias_reader[k].awh;

                /* Each frame and AWH instance extracted generates one xvg file. */
                FILE *out = open_awh_outfile(outf, fr->t, oenv);

                /* Now do the actual printing. Metadata in first subblock is treated separately. */
                fprintf(out, "# AWH metadata: target error = %4.2f kT\n",
                        block->sub[imin_sub].dval[1]);

                fprintf(out, "# AWH metadata: log sample weight = %4.2f\n",
                        block->sub[imin_sub].dval[2]);

                int npoints = block->sub[imin_sub + 1].nr;
                for (int j = 0; j < npoints; j++)
                {
                    fprintf(out, "  %8.4f", block->sub[imin_sub + 1].dval[j]);

                    for (int i = outf->startGraph; i < outf->startGraph + outf->numGraph; i++)
                    {
                        fprintf(out, "  %8.4f", block->sub[i].dval[j]);
                    }
                    fprintf(out, "\n");
                }

                gmx_ffclose(out);

                if (awh_reader->outputMetricTensor)
                {
                    outfile_t *outf = &awh_reader->bias_reader[k].tensor;

                    FILE      *out = open_awh_outfile(outf, fr->t, oenv);
                    for (int j = 0; j < npoints; j++)
                    {
                        fprintf(out, "  %8.4f", block->sub[imin_sub + 1].dval[j]);

                        for (int i = outf->startGraph; i < outf->startGraph + outf->numGraph; i++)
                        {
                            fprintf(out, "  %g", block->sub[i].dval[j]);
                        }
                        fprintf(out, "\n");
                    }

                    gmx_ffclose(out);
                }

                imin_sub += awh_reader->bias_reader[k].nsubblocks;
            }
        }
        (*frame_counter)++;
    }
}

int gmx_awh(int argc, char *argv[])
{
    const char         *desc[] = {
        "[THISMODULE] extracts AWH data from an energy file.",
        "One file is output per AWH bias per time frame.",
        "By default only the PMF is printed. With [TT]-more[tt] the bias,",
        "and target, coordinate distributions and friction are also printed.",
        "Option [TT]-mt[tt] prints all components of the friction metric",
        "tensor to another set of files."
    };
    static gmx_bool     bMore_output = FALSE;
    static int          skip         = 0;
    t_pargs             pa[]         = {
        { "-skip", FALSE, etINT,  {&skip},
          "Skip number of frames between data points" },
        { "-more", FALSE, etBOOL, {&bMore_output},
          "Print more output" }
    };

    ener_file_t         fp;
    int                 timecheck = 0;
    t_inputrec          ir;
    gmx_enxnm_t        *enm = NULL;
    t_enxframe         *frame, *fr = NULL;
    int                 cur = 0;
#define NEXT (1-cur)
    int                 nre, awh_frame_counter;
    gmx_bool            bCont;
    gmx_output_env_t   *oenv;
    awh_reader_t       *awh_reader = NULL;

    /* For reading the tpr containing AWH parameters*/
    gmx_mtop_t         mtop;
    int                natoms;
    matrix             box;

    t_filenm           fnm[] = {
        { efEDR, "-f",    NULL,           ffREAD  },
        { efTPR, "-s",    NULL,           ffOPTRD },
        { efXVG, "-o",    "awh",          ffWRITE },
        { efXVG, "-mt",   "metrictensor", ffOPTWR }
    };
#define NFILE asize(fnm)
    int                npargs;

    npargs = asize(pa);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END,
                           NFILE, fnm, npargs, pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    snew(frame, 2);
    fp = open_enx(ftp2fn(efEDR, NFILE, fnm), "r");
    do_enxnms(fp, &nre, &enm);

    /* We just need the AWH parameters from inputrec. These are used to initialize
       the AWH reader when we have a frame to read later on. */
    read_tpx(ftp2fn(efTPR, NFILE, fnm), &ir, box, &natoms, NULL, NULL, &mtop);

    if (!ir.bDoAwh)
    {
        gmx_fatal(FARGS, "No AWH data in %s\n", opt2fn("-f", NFILE, fnm));
    }

    /* Initiate counters */
    awh_frame_counter   = 0;
    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(fp, &(frame[NEXT]));
            if (bCont)
            {
                timecheck = check_times(frame[NEXT].t);
            }
        }
        while (bCont && (timecheck < 0));

        if ((timecheck == 0) && bCont)
        {
            /* We read a valid frame, so we can use it */
            fr = &(frame[NEXT]);

            if (fr->nre > 0)
            {
                /* The frame contains energies, so update cur */
                cur  = NEXT;
            }

            /* Currently we have the number of subblocks per AWH stored in the first subblock (we cannot get this directly from the tpr)
               so we have to read an AWH block before making the legends. Possibly its better to call init_awh and get the number
               of subblocks and therefore the number of legends that way instead (TODO). */
            if (awh_reader == NULL)
            {
                awh_reader = init_awh_reader(ir.awh_params,
                                             NFILE, fnm, bMore_output,
                                             fr);
            }

            do_awh_frame(fr, awh_reader, &awh_frame_counter, skip, oenv);
        }
    }
    while (bCont && (timecheck == 0));

    fprintf(stderr, "\n");
    close_enx(fp);

    return 0;
}
