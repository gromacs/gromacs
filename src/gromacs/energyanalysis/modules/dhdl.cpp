/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 * \brief
 * Implements classes in dhdl.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "dhdl.h"

#include <cstdio>
#include <cstring>

#include "gromacs/options.h"
#include "gromacs/energyanalysis/analysismodule.h"
#include "gromacs/energyanalysis/energytermcontainer.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace energyanalysis
{

class DhdlEnergyModule : public EnergyAnalysisModule
{
    public:
        //! Constructor
        DhdlEnergyModule();

        virtual void setOutputEnvironment(gmx_output_env_t *oenv)
        {
            ehelper_.setOutputEnvironment(oenv);
        }

        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings);

        virtual void initAnalysis(const std::vector<EnergyNameUnit> &eNU);

        virtual void analyzeFrame(t_enxframe *fr);

        virtual void finalizeAnalysis();

    private:
        //! MD inputrec to extract free energy options from
        t_inputrec                           *inputrec_;
        //! File pointers for output dhdl
        std::unique_ptr<FILE, void(*)(FILE*)> fp_dhdl_;
        //! Energy output xvg file
        std::string                           fnEnergy_;
        //! Topology input file
        std::string                           fnTopology_;
        //! Bookkeeping stuff for histogramming
        int                                  *blocks, *hists, *samples, *nlambdas;
        //! Energy helper class for low level stuff
        EnergyTermContainer                   ehelper_;

};

DhdlEnergyModule::DhdlEnergyModule() : fp_dhdl_(nullptr, xvgrclose)
{
    gmx::MDModules mdModules;
    inputrec_  = mdModules.inputrec();
    blocks     = nullptr;
    hists      = nullptr;
    samples    = nullptr;
    nlambdas   = nullptr;
}

void DhdlEnergyModule::initOptions(IOptionsContainer                 *options,
                                   ICommandLineOptionsModuleSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] extracts and plots the free energy data",
        "(Hamiltonian differences and/or the Hamiltonian derivative dhdl)",
        "from the [TT]ener.edr[tt] file.[PAR]"
    };
    settings->setHelpText(desc);
    // Add input file for topology
    options->addOption(FileNameOption("s")
                           .filetype(eftTopology)
                           .inputFile()
                           .store(&fnTopology_)
                           .defaultBasename("topol")
                           .required());
    // Add option for output files
    options->addOption(FileNameOption("o")
                           .filetype(eftPlot)
                           .outputFile()
                           .store(&fnEnergy_)
                           .defaultBasename("dhdl"));
}

void DhdlEnergyModule::initAnalysis(const std::vector<EnergyNameUnit> &)
{
    // Read the topology
    {
        gmx_mtop_t  mtop;
        int         natoms;
        matrix      box;

        /* all we need is the ir to be able to write the label */
        read_tpx(fnTopology_.c_str(), inputrec_, box, &natoms, NULL, NULL, &mtop);
    }
}

void DhdlEnergyModule::analyzeFrame(t_enxframe *fr)

{
    /* Data to be collected from the enxframe. */
    double temp         = 0, start_time = 0, delta_time = 0, start_lambda = 0;
    int    n_lambda_vec = 0;

    /* now count the blocks & handle the global dh data */
    int nblock_hist   = 0;
    int nblock_dh     = 0;
    int nblock_dhcoll = 0;
    for (int i = 0; i < fr->nblock; i++)
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
                 (fr->block[i].sub[0].nr < 5) )
            {
                GMX_THROW(InvalidInputError("Unexpected block data in energy file"));
            }

            /* read the data from the DHCOLL block */
            temp            = fr->block[i].sub[0].dval[0];
            start_time      = fr->block[i].sub[0].dval[1];
            delta_time      = fr->block[i].sub[0].dval[2];
            start_lambda    = fr->block[i].sub[0].dval[3];

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
                        GMX_THROW(InvalidInputError("Unexpected change of basis set in lambda"));
                    }
                }
            }
        }
    }

    if (nblock_hist == 0 && nblock_dh == 0)
    {
        return;
    }
    static int   setnr             = 0;
    std::string  dhdl("dH/d\\lambda"), deltag("\\DeltaH"), lambda("\\lambda");

    if (nblock_hist > 0 && nblock_dh > 0)
    {
        GMX_THROW(InvalidInputError("This energy file contains both histogram dhdl data and non-histogram dhdl data. Don't know what to do."));
    }
    if (!fp_dhdl_)
    {
        if (nblock_dh > 0)
        {
            /* we have standard, non-histogram data --
               call open_dhdl to open the file */
            /* TODO this is an ugly hack that needs to be fixed: this will only
               work if the order of data is always the same and if we're
               only using the gmx dhdl compiled with the mdrun that produced
               the ener.edr. */
            fp_dhdl_.reset(open_dhdl(fnEnergy_.c_str(),
                                     inputrec_,
                                     ehelper_.outputEnvironment()));
        }
        else
        {
            std::string           label_x = deltag + " (" unit_energy + ")";
            std::string           title   = "N(" + deltag + ")";
            std::string           label_y("Samples");

            fp_dhdl_.reset(xvgropen_type(fnEnergy_.c_str(), title.c_str(),
                                         label_x.c_str(), label_y.c_str(), exvggtXNY,
                                         ehelper_.outputEnvironment()));
            std::string buf = gmx::formatString("T = %g (K), %s = %g",
                                                temp, lambda.c_str(), start_lambda);
            xvgr_subtitle(fp_dhdl_.get(), buf.c_str(),
                          ehelper_.outputEnvironment());
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
        for (int i = 0; i < fr->nblock; i++)
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
                    GMX_THROW(InvalidInputError("Unexpected block data in file"));
                }
                foreign_lambda = blk->sub[0].dval[0];
                dx             = blk->sub[0].dval[1];
                nhist          = blk->sub[1].lval[0];
                derivative     = blk->sub[1].lval[1];
                for (int j = 0; j < nhist; j++)
                {
                    const char *lg[1];
                    std::string legend;

                    if (!derivative)
                    {
                        legend = formatString("N(%s(%s=%g) | %s=%g)",
                                              deltag.c_str(), lambda.c_str(), foreign_lambda,
                                              lambda.c_str(), start_lambda);
                    }
                    else
                    {
                        legend = formatString("N(%s | %s=%g)",
                                              dhdl.c_str(), lambda.c_str(), start_lambda);
                    }

                    lg[0] = legend.c_str();
                    xvgr_new_dataset(fp_dhdl_.get(), setnr, 1, lg, ehelper_.outputEnvironment());
                    setnr++;
                    x0 = blk->sub[1].lval[2+j];
                    for (int k = 0; k < blk->sub[j+2].nr; k++)
                    {
                        int    hist;
                        double xmin, xmax;

                        hist = blk->sub[j+2].ival[k];
                        xmin = (x0+k)*dx;
                        xmax = (x0+k+1)*dx;
                        fprintf(fp_dhdl_.get(), "%g %d\n%g %d\n", xmin, hist,
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

        for (int i = 0; i < fr->nblock; i++)
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
                        GMX_THROW(InvalidInputError("Length inconsistency in dhdl data"));
                    }
                }
            }
        }
        (*samples) += len;

        for (int i = 0; i < len; i++)
        {
            double time = start_time + delta_time*i;

            fprintf(fp_dhdl_.get(), "%.4f ", time);

            for (int j = 0; j < fr->nblock; j++)
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

                    if (j == 1 && inputrec_->bExpanded)
                    {
                        fprintf(fp_dhdl_.get(), "%4d", (int)value);   /* if expanded ensembles and zero, this is a state value, it's an integer. We need a cleaner conditional than if j==1! */
                    }
                    else
                    {
                        if (ehelper_.doublePrecision())
                        {
                            fprintf(fp_dhdl_.get(), " %#.12g", value);   /* print normal precision */
                        }
                        else
                        {
                            fprintf(fp_dhdl_.get(), " %#.8g", value);   /* print normal precision */
                        }
                    }
                }
            }
            fprintf(fp_dhdl_.get(), "\n");
        }
    }
}

void DhdlEnergyModule::finalizeAnalysis()
{
}

const char DhdlInfo::name[] = "dhdl";

const char DhdlInfo::shortDescription[] = "Extract dH/dlambda terms from energy file";

EnergyAnalysisModulePointer DhdlInfo::create()
{
    return EnergyAnalysisModulePointer(new DhdlEnergyModule());
}

} // namespace energyanalysis

} // namespace gmx
