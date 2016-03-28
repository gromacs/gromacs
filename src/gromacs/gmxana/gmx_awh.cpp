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

typedef struct bias_reader_t {
    int           nsubblocks;                /* Number of subblocks to read. */
    int           nsubblocks_out;            /* Number of subblocks to print. */
    int           nlegend;                   /* Number of legend strings. */
    char        **legend;                    /* Legends for the output. */
    char         *ylabel;                    /* Label for the y-axis. */
} bias_reader_t;

typedef struct awh_reader_t {
    int               nbias;                  /* The number of AWH biases. */
    bias_reader_t    *bias_reader;            /* The readers, one for each AWH bias. */
} awh_reader_t;


static void make_legend(awh_bias_params_t *awh_bias_params, int nleg, char **leg)
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

    if (ileg != nleg)
    {
        gmx_incons("Mismatch between the number of legends requested for printing and the number present!");
    }
}

static awh_reader_t *init_awh_reader(awh_params_t *awh_params, t_enxframe *fr, gmx_bool bMore)
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

        /* Keep track of the first subblock of this AWH */
        imin_sub = 0;
        for (int k = 0; k < awh_reader->nbias; k++)
        {
            awh_bias_params_t *awh_bias_params = &awh_params->awh_bias_params[k];
            bias_reader_t     *bias_reader     = &awh_reader->bias_reader[k];

            bias_reader->nsubblocks = (int)block->sub[imin_sub].dval[0];

            /* We at least print the grid point values (ndim) and the PMF (1) and the metadata (1). */
            bias_reader->nsubblocks_out = bMore ? bias_reader->nsubblocks : awh_params->awh_bias_params[k].ndim + 2;

            /* The metadata and the grid point values (for the first dimension) do not get legends (-2) */
            bias_reader->nlegend = bias_reader->nsubblocks_out - 2;

            snew(bias_reader->legend, bias_reader->nlegend);
            make_legend(awh_bias_params, bias_reader->nlegend, bias_reader->legend);

            bias_reader->ylabel = strdup("(k\\sB\\NT)");
            imin_sub           += bias_reader->nsubblocks;
        }
    }

    return awh_reader;
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
                int            i, j, npoints;
                FILE          *out = NULL;
                char           filename[STRLEN], title[STRLEN];
                bias_reader_t *bias_reader = &awh_reader->bias_reader[k];

                /* Each frame and AWH instance extracted generates one xvg file. */
                if (awh_reader->nbias > 1)
                {
                    sprintf(filename, "awh%d_t%g.xvg", k + 1, fr->t);
                    sprintf(title, "AWH %d",  k + 1);
                }
                else
                {
                    sprintf(filename, "awh_t%g.xvg", fr->t);
                    sprintf(title, "AWH");
                }

                out = xvgropen(filename, title,
                               "(nm or deg)", bias_reader->ylabel, oenv);
                xvgr_legend(out, bias_reader->nlegend, (const char**)bias_reader->legend, oenv);

                /* Now do the actual printing. Metadata in first subblock is treated separately. */
                fprintf(out, "# AWH metadata: target error = %3.2f kT\n",
                        block->sub[imin_sub].dval[1]);

                fprintf(out, "# AWH metadata: log sample weight = %3.2f \n",
                        block->sub[imin_sub].dval[2]);

                npoints = block->sub[imin_sub + 1].nr;
                for (j = 0; j < npoints; j++)
                {
                    for (i = imin_sub + 1; i < imin_sub + bias_reader->nsubblocks_out; i++)
                    {
                        fprintf(out, "  %8.4f", block->sub[i].dval[j]);
                    }

                    fprintf(out, "\n");
                }

                gmx_ffclose(out);
                imin_sub += bias_reader->nsubblocks;
            }
        }
        (*frame_counter)++;
    }
}

int gmx_awh(int argc, char *argv[])
{
    const char         *desc[] = {
        "[THISMODULE] extracts AWH data from an energy file.[PAR]",
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
        { efEDR, "-f",    NULL,      ffREAD  },
        { efTPR, "-s",    NULL,      ffOPTRD },
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
                awh_reader = init_awh_reader(ir.awh_params, fr, bMore_output);
            }

            do_awh_frame(fr, awh_reader, &awh_frame_counter, skip, oenv);
        }
    }
    while (bCont && (timecheck == 0));

    fprintf(stderr, "\n");
    close_enx(fp);

    return 0;
}
