/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include <stdio.h>
#include <string.h>
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/options.h"
#include "dhdl.h"

namespace gmx
{

void DhdlEnergy::setParameters(const char *topnm)
{
    gmx_mtop_t  mtop;
    int         natoms;
    matrix      box;

    /* all we need is the ir to be able to write the label */
    read_tpx(topnm, &ir, box, &natoms, NULL, NULL, NULL, &mtop);

    // Implement done mtop? There is no such routine!
}

DhdlEnergy::DhdlEnergy()
{
    init_inputrec(&ir);
    fp_dhdl  = NULL;
    filename = NULL;
    blocks   = NULL;
    hists    = NULL;
    samples  = NULL;
    nlambdas = NULL;
}

void DhdlEnergy::initOptions(Options *options)
{
    static const char *const desc[] = {
        "[THISMODULE] extracts and plots the free energy data",
        "(Hamiltoian differences and/or the Hamiltonian derivative dhdl)",
        "from the [TT]ener.edr[tt] file.[PAR]"
    };
    options->setDescription(desc);
}

bool DhdlEnergy::initAnalysis(std::vector<std::string> eName,
                              std::vector<std::string> eUnit)
{
    return SimpleEnergy::initAnalysis(eName, eUnit);
}

bool DhdlEnergy::addAnalysisFrame(t_enxframe *fr)

{
    const char  *dhdl = "dH/d\\lambda", *deltag = "\\DeltaH", *lambda = "\\lambda";
    char         title[STRLEN], label_x[STRLEN], label_y[STRLEN], legend[STRLEN];
    char         buf[STRLEN];
    int          nblock_hist = 0, nblock_dh = 0, nblock_dhcoll = 0;
    int          i, j, k;
    /* coll data */
    double       temp              = 0, start_time = 0, delta_time = 0, start_lambda = 0; //, delta_lambda = 0;
    static int   setnr             = 0;
    double      *native_lambda_vec = NULL;
    const char **lambda_components = NULL;
    int          n_lambda_vec      = 0;
    //int          lambda_fep_state;

    /* now count the blocks & handle the global dh data */
    for (i = 0; i < fr->nblock; i++)
    {
        if (fr->block[i].id == enxDHHIST)
        {
            nblock_hist++;
        }
        else if (fr->block[i].id == enxDH)
        {
            nblock_dh++;
        }
        else if (fr->block[i].id == enxDHCOLL)
        {
            nblock_dhcoll++;
            if ( (fr->block[i].nsub < 1) ||
                 (fr->block[i].sub[0].type != xdr_datatype_double) ||
                 (fr->block[i].sub[0].nr < 5))
            {
                gmx_fatal(FARGS, "Unexpected block data");
            }

            /* read the data from the DHCOLL block */
            temp            = fr->block[i].sub[0].dval[0];
            start_time      = fr->block[i].sub[0].dval[1];
            delta_time      = fr->block[i].sub[0].dval[2];
            start_lambda    = fr->block[i].sub[0].dval[3];
            //delta_lambda    = fr->block[i].sub[0].dval[4];
            if (fr->block[i].nsub > 1)
            {
                //lambda_fep_state = fr->block[i].sub[1].ival[0];
                if (n_lambda_vec == 0)
                {
                    n_lambda_vec = fr->block[i].sub[1].ival[1];
                }
                else
                {
                    if (n_lambda_vec != fr->block[i].sub[1].ival[1])
                    {
                        gmx_fatal(FARGS,
                                  "Unexpected change of basis set in lambda");
                    }
                }
                if (lambda_components == NULL)
                {
                    snew(lambda_components, n_lambda_vec);
                }
                if (native_lambda_vec == NULL)
                {
                    snew(native_lambda_vec, n_lambda_vec);
                }
                for (j = 0; j < n_lambda_vec; j++)
                {
                    native_lambda_vec[j] = fr->block[i].sub[0].dval[5+j];
                    lambda_components[j] =
                        efpt_singular_names[fr->block[i].sub[1].ival[2+j]];
                }
            }
        }
    }

    if (nblock_hist == 0 && nblock_dh == 0)
    {
        /* don't do anything */
        return false;
    }
    if (nblock_hist > 0 && nblock_dh > 0)
    {
        gmx_fatal(FARGS, "This energy file contains both histogram dhdl data and non-histogram dhdl data. Don't know what to do.");
    }
    if (!*fp_dhdl)
    {
        if (nblock_dh > 0)
        {
            /* we have standard, non-histogram data --
               call open_dhdl to open the file */
            /* TODO this is an ugly hack that needs to be fixed: this will only
               work if the order of data is always the same and if we're
               only using the g_energy compiled with the mdrun that produced
               the ener.edr. */
            *fp_dhdl = open_dhdl(filename, &ir, helper()->getOutputEnvironment());
        }
        else
        {
            sprintf(title, "N(%s)", deltag);
            sprintf(label_x, "%s (%s)", deltag, unit_energy);
            sprintf(label_y, "Samples");
            *fp_dhdl = xvgropen_type(filename, title, label_x, label_y, exvggtXNY,
                                     helper()->getOutputEnvironment());
            sprintf(buf, "T = %g (K), %s = %g", temp, lambda, start_lambda);
            xvgr_subtitle(*fp_dhdl, buf,
                          helper()->getOutputEnvironment());
        }
    }

    (*hists)   += nblock_hist;
    (*blocks)  += nblock_dh;
    (*nlambdas) = nblock_hist+nblock_dh;

    /* write the data */
    if (nblock_hist > 0)
    {
        gmx_int64_t sum = 0;
        /* histograms */
        for (i = 0; i < fr->nblock; i++)
        {
            t_enxblock *blk = &(fr->block[i]);
            if (blk->id == enxDHHIST)
            {
                double          foreign_lambda, dx;
                gmx_int64_t     x0;
                int             nhist, derivative;

                /* check the block types etc. */
                if ( (blk->nsub < 2) ||
                     (blk->sub[0].type != xdr_datatype_double) ||
                     (blk->sub[1].type != xdr_datatype_int64) ||
                     (blk->sub[0].nr < 2)  ||
                     (blk->sub[1].nr < 2) )
                {
                    gmx_fatal(FARGS, "Unexpected block data in file");
                }
                foreign_lambda = blk->sub[0].dval[0];
                dx             = blk->sub[0].dval[1];
                nhist          = blk->sub[1].lval[0];
                derivative     = blk->sub[1].lval[1];
                for (j = 0; j < nhist; j++)
                {
                    const char *lg[1];
                    x0 = blk->sub[1].lval[2+j];

                    if (!derivative)
                    {
                        sprintf(legend, "N(%s(%s=%g) | %s=%g)",
                                deltag, lambda, foreign_lambda,
                                lambda, start_lambda);
                    }
                    else
                    {
                        sprintf(legend, "N(%s | %s=%g)",
                                dhdl, lambda, start_lambda);
                    }

                    lg[0] = legend;
                    xvgr_new_dataset(*fp_dhdl, setnr, 1, lg, helper()->getOutputEnvironment());
                    setnr++;
                    for (k = 0; k < blk->sub[j+2].nr; k++)
                    {
                        int    hist;
                        double xmin, xmax;

                        hist = blk->sub[j+2].ival[k];
                        xmin = (x0+k)*dx;
                        xmax = (x0+k+1)*dx;
                        fprintf(*fp_dhdl, "%g %d\n%g %d\n", xmin, hist,
                                xmax, hist);
                        sum += hist;
                    }
                    /* multiple histogram data blocks in one histogram
                       mean that the second one is the reverse of the first one:
                       for dhdl derivatives, it's important to know both the
                       maximum and minimum values */
                    dx = -dx;
                }
            }
        }
        (*samples) += (int)(sum/nblock_hist);
    }
    else
    {
        /* raw dh */
        int    len      = 0;

        for (i = 0; i < fr->nblock; i++)
        {
            t_enxblock *blk = &(fr->block[i]);
            if (blk->id == enxDH)
            {
                if (len == 0)
                {
                    len = blk->sub[2].nr;
                }
                else
                {
                    if (len != blk->sub[2].nr)
                    {
                        gmx_fatal(FARGS, "Length inconsistency in dhdl data");
                    }
                }
            }
        }
        (*samples) += len;

        for (i = 0; i < len; i++)
        {
            double time = start_time + delta_time*i;

            fprintf(*fp_dhdl, "%.4f ", time);

            for (j = 0; j < fr->nblock; j++)
            {
                t_enxblock *blk = &(fr->block[j]);
                if (blk->id == enxDH)
                {
                    double value;
                    if (blk->sub[2].type == xdr_datatype_float)
                    {
                        value = blk->sub[2].fval[i];
                    }
                    else
                    {
                        value = blk->sub[2].dval[i];
                    }
                    /* we need to decide which data type it is based on the count*/

                    if (j == 1 && ir.bExpanded)
                    {
                        fprintf(*fp_dhdl, "%4d", (int)value);   /* if expanded ensembles and zero, this is a state value, it's an integer. We need a cleaner conditional than if j==1! */
                    }
                    else
                    {
                        if (helper()->getDoublePrecision())
                        {
                            fprintf(*fp_dhdl, " %#.12g", value);   /* print normal precision */
                        }
                        else
                        {
                            fprintf(*fp_dhdl, " %#.8g", value);   /* print normal precision */
                        }
                    }
                }
            }
            fprintf(*fp_dhdl, "\n");
        }
    }
    return true;
}

bool DhdlEnergy::finalizeAnalysis()
{
    return true;
}


}
