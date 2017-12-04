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

#include "gmxpre.h"
#include "config.h"

#include "string.h"
#include "densfit.h"
#include "localdensfitdata.h"

#include "gromacs/math/vec.h"

#include "gromacs/fileio/xvgr.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/mdlib/sim_util.h"

#include "gromacs/timing/cyclecounter.h"

#include "gromacs/domdec/domdec.h"

#include "gromacs/fileio/griddataio.h"

#include "gromacs/mdlib/broadcaststructs.h"

#define MAXPTR 254

const std::string Densfit::name_ = "Density fitting:";

Densfit::Densfit() {}

Densfit::Densfit(const gmx::DensfitData &parameters)
{
    parameters_ = std::unique_ptr<gmx::DensfitData>(new gmx::DensfitData(parameters));
}

Densfit::~Densfit(){};

const gmx::GridDataReal3D &Densfit::simulatedMap() { return df_->simulatedMap(); }

void Densfit::init(int ePBC, gmx_mtop_t *mtop, rvec *x, matrix box, gmx::LocalAtomSetManager * atomsets, const t_commrec *cr)
{

    pbc_ = ePBC;
    if ((PAR(cr)) && !DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS, "%sOnly implemented for domain decomposition.\n", name_.c_str());
    }

    auto densfitAtoms = atomsets->add(parameters_->fittingGroup().ind_);

    /* To be able to make the density fitting molecule whole: */
    rvec *x_densfit_pbc = nullptr;
    snew(x_densfit_pbc, densfitAtoms.numAtomsGlobal()); /* There ... */
    if (MASTER(cr))
    {
        rvec *x_pbc = nullptr;
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        snew(x_pbc, mtop->natoms); /* There ... */
        //        m_rveccopy(mtop->natoms,x,x_pbc);
        copy_rvecn(x, x_pbc, 0, mtop->natoms);
        do_pbc_first_mtop(nullptr, pbc_, box, mtop, x_pbc);
        for (size_t i = 0; i < densfitAtoms.numAtomsGlobal(); i++)
        {
            copy_rvec(x_pbc[densfitAtoms.globalIndex()[i]], x_densfit_pbc[i]);
        }
        sfree(x_pbc); /* ... and back again */
    }
    if (PAR(cr))
    {
        gmx_bcast(parameters_->fittingGroup().ind_.size() * sizeof(rvec), x_densfit_pbc, cr);
    }

    df_ = std::unique_ptr<gmx::LocalDensfitData>(new gmx::LocalDensfitData(densfitAtoms, *parameters_, x_densfit_pbc));

    sfree(x_densfit_pbc); /* ... and back again */

}

void Densfit::initOutputFromCommandLineParameters(FILE *fplog, int nfile, const t_filenm *fnm, const gmx_output_env_t *oenv, gmx_bool bAppend, gmx_bool bVerbose)
{
    output_ = std::unique_ptr<gmx::DensfitOutput> (new gmx::DensfitOutput(fplog, opt2fn("-d", nfile, fnm), opt2fn("-mo", nfile, fnm), oenv, bAppend));
    if (bVerbose && output_ && parameters_ != nullptr)
    {
        parameters_->printToOut(output_->logFile(), oenv);
    }
}

void Densfit::triggerShiftUpdate()
{
    df_->triggerShiftUpdate();
}

void Densfit::add_forces(rvec *f, gmx_int64_t step, real time)
{
    if (!parameters_->fitThisStep(step))
    {
        return;
    }

    df_->add_forces(f);

    if (output_ && output_->logFile() && parameters_->outputMapThisStep(step))
    {
        fprintf(output_->logFile(), "%12.5e%s\n", time, df_->infoString().c_str());
    }

}

void Densfit::do_densfit(real t, gmx_int64_t step,
                         const t_commrec *cr, gmx::PaddedArrayRef<gmx::RVec> x, matrix box,
                         gmx_wallcycle *wcycle)
{

    /* If the user requested to have the density fitting forces calculated only
     * every N steps, we can skip the expensive do_densfit routine. However, we
     * must always calculate the fitting forces
     * - at the first step (obviously)
     * - after domain decompositions since the local atoms have changed
     */
    if ((!parameters_->fitThisStep(step) && (!df_->afterDomainDecomposition())))
    {
        return;
    }

    /* Update the time-dependent parameters like k, sigma, etc. */
    if (parameters_->timeDependent().vary())
    {
        df_->updateTimeDependentData(t, *parameters_);
    }


    /* ASSEMBLE THE POSITIONS OF THE ATOMS USED FOR DENSITY FITTING           */
    /* Broadcast the DF positions such that every node has all of them
     * Every node contributes its local positions x and stores it in
     * the collective df->x_assembled array. Do all the communication first!  */
    wallcycle_start(wcycle, ewcDENSFIT_COMM);
    df_->communicate(x, cr, box, pbc_);
    wallcycle_stop(wcycle, ewcDENSFIT_COMM);
    /* Now all nodes have all of the DF positions in df->x_assembled */

    /* SPREAD THE ATOMS ON THE DENSITY GRID                                   */
    /* Produces the density map. Each node spreads an equal part of all atoms,
     * this way we do not need load balancing here                            */
    wallcycle_start(wcycle, ewcDENSFIT_SPREAD);
    df_->spreadAtoms(cr);
    wallcycle_stop(wcycle, ewcDENSFIT_SPREAD);

    if (PAR(cr))
    {
        wallcycle_start(wcycle, ewcDENSFIT_SUM_GRID);
        df_->sumSimulatedMap(cr);
        wallcycle_stop(wcycle, ewcDENSFIT_SUM_GRID);
    }


    if (output_ && parameters_->outputMapThisStep(step))
    {
        if (parameters_->keepAndNumberMaps())
        {
            gmx::MrcFile().write(output_->outputMapFileNameWithStep(step), df_->simulatedMap());
        }
        else
        {
            gmx::MrcFile().write(output_->outputMapFileName(), df_->simulatedMap());
        }
    }

    if (MASTER(cr))
    {
        df_->calculateCorrelationCoefficient(parameters_->referenceMap());
    }
    /**************************************************************************/
    /* CALCULATE THE FORCES FROM THE DENSITY FITTING POTENTIAL                */
    /* Each rank calculates the forces on its home atoms                      */

    /* Start to count cycles for the dynamic load balancing now */
    gmx_cycles_t cycles_comp = gmx_cycles_read(); /* Cycles for the density fitting computations
                                                     only, does not count communication. This
                                                     counter is used for load-balancing */

    //    long int c1 = gmx_cycles_read();
    wallcycle_start(wcycle, ewcDENSFIT_FORCES);
    df_->do_forces(box, parameters_->referenceMap());
    wallcycle_stop(wcycle, ewcDENSFIT_FORCES);
    //    fprintf(stderr, "--- Forces took %g M cycles\n",
    //    0.000001*(gmx_cycles_read()-c1));

    /* Stop the density fitting cycle counter and add the computation-only
     * cycles to the force cycles for load balancing */
    cycles_comp = gmx_cycles_read() - cycles_comp;

    if (DOMAINDECOMP(cr) && wcycle)
    {
        dd_cycles_add(cr->dd, cycles_comp, ddCyclF);
    }
}

gmx::GridDataReal3D Densfit::referenceMapCopy() const
{
    return parameters_->referenceMap();
};

const gmx::DensfitData &Densfit::parameters() const { return *parameters_; }

namespace gmx
{

void make_filename(const char *outf_name, int ndigit, int file_nr,
                   char newname[STRLEN])
{
    char  fmt[128];
    char  nbuf[128];

    /* Strip away the extension */
    char *outf_base = gmx_strdup(outf_name);
    char *extpos    = strrchr(outf_base, '.');
    if (extpos)
    {
        *extpos = '\0';
    }

    int fnr = file_nr;
    int nd  = 0;
    do
    {
        fnr /= 10;
        nd++;
    }
    while (fnr > 0);

    if (ndigit > 0)
    {
        snprintf(fmt, 128, "%%0%dd", ndigit);
        snprintf(nbuf, 128, fmt, file_nr);
    }
    else
    {
        snprintf(nbuf, 128, "%d", file_nr);
    }

    if (extpos)
    {
        snprintf(newname, STRLEN, "%s%s.%s", outf_base, nbuf, extpos + 1);
    }
    else
    {
        snprintf(newname, STRLEN, "%s%s", outf_base, nbuf);
    }

    free(outf_base);
}


DensfitOutput::DensfitOutput(FILE * fplog, const char *fnLogFile, const char * fnDensity, const gmx_output_env_t *oenv, bool bAppend)
{
    if (bAppend)
    {
        logFile_ = gmx_fio_fopen(fnLogFile, "a");
    }
    else
    {
        please_cite(fplog, "Tama2008");
        logFile_ = xvgropen(fnLogFile, "Density fitting correlation coefficient", "Time (ps)",
                            "correlation coefficient", oenv);
        fflush(logFile_);
    }
    fnDensity_ = fnDensity;
}

FILE* DensfitOutput::logFile()
{
    return logFile_;
}

DensfitOutput::~DensfitOutput()
{
    if (logFile_) // Close the density fitting output file
    {
        if (bAppend_)
        {
            gmx_fio_fclose(logFile_);
        }
        else
        {
            xvgrclose(logFile_);
        }
    }
}

std::string DensfitOutput::outputMapFileName() const
{
    return fnDensity_;
}

std::string DensfitOutput::outputMapFileNameWithStep(gmx_int64_t step) const
{
    char         fnWithStep[STRLEN];
    make_filename(fnDensity_.c_str(), 0, step, fnWithStep);
    return fnWithStep;
}
}
