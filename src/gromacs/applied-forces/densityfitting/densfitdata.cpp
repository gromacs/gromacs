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

#include "densfitdata.h"

#include <cassert>
#include <cstdio>
#include <vector>

#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdlib/sim_util.h"

#include "gromacs/mdlib/groupcoord.h"

#include "gromacs/domdec/domdec_struct.h"

#include "gromacs/timing/cyclecounter.h"

#include "gromacs/domdec/domdec.h"

#include "gromacs/fileio/griddataio.h"

#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/pbcutil/pbc.h"

#include "gromacs/fileio/mrcmetadata.h"

#include "gromacs/utility/txtdump.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/utility/strconvert.h"

#include "gromacs/fileio/gmxfio-xdr.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/fileio/readinp.h"

static const char *FitStr = "";                       //!< Used in gmx map output
//! \brief Environment variable for setting nstfit
static const char *NSTFIT_ENVVAR = "GMX_NSTFIT";

#define MAXPTR 254

namespace gmx
{

namespace
{
//! \brief Adds 'buf' to 'str'
static void add_to_string(char **str, char *buf)
{
    int len;

    len = strlen(*str) + strlen(buf) + 1;
    srenew(*str, len);
    strcat(*str, buf);
}

//! \brief TODO
static void add_to_string_aligned(char **str, const char *buf)
{
    char buf_aligned[STRLEN];

    sprintf(buf_aligned, "%12s", buf);

    add_to_string(str, buf_aligned);
}
}

DensfitData::DensfitData() {}

DensfitData::DensfitData(int npoint, real sigma, real k,
                         real sigma_dist, int nAtoms, int *index,
                         bool bVerbose)
{
    npoints_ = npoint;
    snew(time_values_, npoints_);
    snew(temp_values_, npoints_);
    snew(k_values_, npoints_);
    snew(sigma_values_, npoints_);
    sigma_values_[0] = sigma;
    k_values_[0]     = k;
    dist_            = sigma_dist;
    nat_             = nAtoms;

    if (NULL == ind_)
    {
        if (bVerbose)
        {
            fprintf(stdout, "%sYou did not provide an index group - will use all "
                    "atoms for spreading.\n",
                    FitStr);
        }

        snew(ind_, nat_);
        for (int i = 0; i < nat_; i++)
        {
            ind_[i] = i;
        }
    }
    else
    {
        ind_ = index;
    }
}

void DensfitData::makeGroups(char *groupname, t_blocka *grps,
                             char **gnames)
{
    /* Group index from group name */
    int ig = search_string(groupname, grps->nr, gnames);
    nat_ = grps->index[ig + 1] - grps->index[ig];

    fprintf(stderr, "Density fitting group '%s' has %d atoms\n", groupname, nat_);
    snew(ind_, nat_);
    for (int i = 0; i < nat_; i++)
    {
        ind_[i] = grps->a[grps->index[ig] + i];
    }
}

real DensfitData::currentK(real time) const
{
    return interpolateLinearly(time, npoints_, time_values_, k_values_);
}
real DensfitData::currentT(real time) const
{
    return interpolateLinearly(time, npoints_, time_values_, temp_values_);
}
real DensfitData::currentSigma(real time) const
{
    return interpolateLinearly(time, npoints_, time_values_, sigma_values_);
}

//! \brief TODO
real DensfitData::interpolateLinearly(real time, int npoints, real *time_values,
                                      real *y_values) const
{
    real x, y;
    int  j;

    /* If we are left of any time point, return the zero'th y value */
    if (time <= time_values[0])
    {
        return y_values[0];
    }

    /* If we are to the right of any time point, return the last y value */
    if (time >= time_values[npoints - 1])
    {
        return y_values[npoints - 1];
    }

    /* We are somewhere in between, that means we need to interpolate y linearly
     */
    j = 0;
    while (time > time_values[j + 1])
    {
        j++;
        assert(j < npoints);
    }

    /* Found our position between points j and j+1.
     * Interpolate: x is the amount from j+1, (1-x) from point j.
     */
    x = ((time - time_values[j]) / (time_values[j + 1] - time_values[j]));
    y = x * y_values[j + 1] + (1 - x) * y_values[j];

    return y;
}

int DensfitData::nAtoms() const
{
    return nat_;
}

void DensfitData::printToOut(FILE *fp, const gmx_output_env_t * oenv) const
{
    fprintf(fp, "# Density fitting forces for %d atoms are updated every %d "
            "time step%s.\n",
            nat_, nstfit_,
            nstfit_ != 1 ? "s" : "");
    fprintf(fp, "#   Note that the force update frequency can be overwritten "
            "by the environment variable %s.\n",
            NSTFIT_ENVVAR);
    fprintf(fp, "# Output is written in intervals of %d time step%s.\n",
            nstout_, nstout_ > 1 ? "s" : "");
    fprintf(fp, "# Number of refinement steps is %d.\n", npoints_);
    fprintf(fp, "# Parameter sigma_dist = %g\n", dist_);
    if (0 != nstmapout_)
    {
        fprintf(fp, "# Calculated map will be saved every %d steps\n",
                nstmapout_);
    }
    int          nsets;
    const char **setname;
    char         buf[50];
    char        *LegendStr = NULL;
    /* Print a nice legend */
    snew(LegendStr, 1);
    LegendStr[0] = '\0';
    sprintf(buf, "#     %6s", "time");
    add_to_string_aligned(&LegendStr, buf);

    nsets = 0;
    snew(setname, 5);
    setname[nsets++] = strdup("temperature (K)");
    add_to_string_aligned(&LegendStr, "T");
    setname[nsets++] = strdup("k (kJ/mol)");
    add_to_string_aligned(&LegendStr, "k");
    setname[nsets++] = strdup("sigma (nm)");
    add_to_string_aligned(&LegendStr, "sigma");
    setname[nsets++] = strdup("correlation coefficient");
    add_to_string_aligned(&LegendStr, "c.c.");
    setname[nsets++] = strdup("V_dens (kJ/mol)");
    add_to_string_aligned(&LegendStr, "V_dens");
    xvgr_legend(fp, nsets, setname, oenv);
    sfree(setname);

    fprintf(fp, "#\n# Legend for the following data columns:\n");
    fprintf(fp, "%s\n", LegendStr);
    sfree(LegendStr);

}

void DensfitData::setNstFitFromEnv(const t_commrec *cr, int nstlist, FILE *fp)
{
    char *env;
    char  buf[STRLEN];
    int   env_nstfit;

    if (MASTER(cr))
    {
        env = getenv(NSTFIT_ENVVAR);

        if (env != NULL)
        {
            sprintf(buf, "Getting nstfit from environment variable %s=%s",
                    NSTFIT_ENVVAR, env);
            fprintf(stderr, "%s%s\n", FitStr, buf);

            sscanf(env, "%d", &env_nstfit);

            if (env_nstfit >= 0)
            {
                nstfit_ = env_nstfit;
                if (fp)
                {
                    fprintf(fp, "# %s\n", buf);
                }
            }
            else
            {
                fprintf(stderr, "WARNING: Could not get a meaningful value for "
                        "environment variable %s!\n",
                        NSTFIT_ENVVAR);
            }
        }
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(nstfit_), &(nstfit_), cr);
    }

    /* We check nstfit versus nstlist here so we can output a warning on stderr
     * AND a note in the density fitting output file */
    if (MASTER(cr))
    {
        if (nstfit_ > 1 && (nstfit_ % nstlist != 0))
        {
            fprintf(stderr, "\n%s\n"
                    "WARNING: The frequency at which the fitting forces are "
                    "calculated (nstfit=%d)\n"
                    "         is not a multiple of the neighbor searching "
                    "frequency (nstlist=%d).\n"
                    "         Fitting forces will be calculated when the "
                    "time step is a multiple\n"
                    "         of nstfit or nstlist, that could also be "
                    "irregular intervals.\n\n",
                    FitStr, nstfit_, nstlist);

            fprintf(fp, "# Note: will update fitting forces when the time step is a "
                    "multiple of nstlist=%d or nstfit=%d.\n",
                    nstlist, nstfit_);
        }

        if (nstfit_ < 1)
        {
            fprintf(stderr, "\n%s\n"
                    "WARNING: The frequency at which the density fitting "
                    "potential is calculated (nstfit=%d) is < 1!\n"
                    "         No density fitting will be performed and no "
                    "fitting output will be written.\n\n",
                    FitStr, nstfit_);
            fprintf(fp, "# No output follows since nstfit < 1\n");
        }
    }

}

void DensfitData::broadcast(const t_commrec *cr)
{
    for (int i = 0; i < npoints_; i++)
    {
        snew_bc(cr, time_values_, npoints_);
        snew_bc(cr, temp_values_, npoints_);
        snew_bc(cr, k_values_, npoints_);
        snew_bc(cr, sigma_values_, npoints_);
    }
    nblock_bc(cr, npoints_, time_values_);
    nblock_bc(cr, npoints_, temp_values_);
    nblock_bc(cr, npoints_, k_values_);
    nblock_bc(cr, npoints_, sigma_values_);

    snew_bc(cr, ind_, nat_);
    nblock_bc(cr, nat_, ind_);
    block_bc(cr, map_ref_);

    int grid1d = map_ref_.grid[XX] * map_ref_.grid[YY] * map_ref_.grid[ZZ];
    if (!MASTER(cr))
    {
        /* For now it is sufficient to have the map title on the master node only */
        map_ref_.title = nullptr;
        map_ref_.vox.resize(grid1d);
    }
    nblock_bc(cr, grid1d, map_ref_.vox.data());
}
void DensfitData::setReferenceMap(const gmx::t_mapdata &reference)
{
    map_ref_ = reference;
}

const t_mapdata &DensfitData::referenceMap()
{
    return map_ref_;
}

bool DensfitData::keepAndNumberMaps() const
{
    return bKeepAndNumberMaps_;
};

int DensfitData::nStepsMapOutput() const
{return nstmapout_; }


int DensfitData::timePoints() const
{
    return npoints_;
}

int DensfitData::nStepsFit() const
{
    return nstfit_;
}

int * DensfitData::indices() const
{
    return ind_;
}

void DensfitData::print(FILE * fp, int indent) const
{
    fprintf(fp, "density fitting options:\n");
    pr_int(fp, indent, "densfit-npoints", npoints_ );
    pr_reals(fp, indent, "densfit-time ", time_values_, npoints_);
    pr_reals(fp, indent, "densfit-sigma", sigma_values_, npoints_);
    pr_reals(fp, indent, "densfit-k    ", k_values_, npoints_);
    pr_reals(fp, indent, "densfit-temp ", temp_values_, npoints_);
    pr_real(fp, indent, "densfit-dist", dist_);
    pr_int(fp, indent, "densfit-nstfit", nstfit_);
    pr_int(fp, indent, "densfit-nstout", nstout_);
    pr_int(fp, indent, "densfit-nstmapout", nstmapout_);
    pr_str(fp, indent, "densfit-bKeepMaps",  boolToString(bKeepAndNumberMaps_));
    pr_ivec_block(fp, indent, "atom", ind_, nat_, true);
}

void DensfitData::do_fio(t_fileio *fio, bool bRead)
{

    gmx_fio_do_int(fio, npoints_);
    gmx_fio_do_int     (fio, nat_);

    if (bRead)
    {
        snew(time_values_, npoints_);
        snew(k_values_, npoints_);
        snew(sigma_values_, npoints_);
        snew(temp_values_, npoints_);

        snew(ind_, nat_);
    }

    gmx_fio_ndo_real   (fio, time_values_, npoints_);
    gmx_fio_ndo_real   (fio, k_values_, npoints_);
    gmx_fio_ndo_real   (fio, sigma_values_, npoints_);
    gmx_fio_ndo_real   (fio, temp_values_, npoints_);

    gmx_fio_do_real    (fio, dist_              );
    gmx_fio_do_int     (fio, nstfit_            );
    gmx_fio_do_int     (fio, nstout_            );
    gmx_fio_do_int     (fio, nstmapout_         );
    gmx_fio_do_gmx_bool(fio, bKeepAndNumberMaps_);

    gmx_fio_ndo_int    (fio, ind_, nat_);
}

void DensfitData::read_densparams(int *ninp_p, t_inpfile **inp_p, warninp * wi, gmx_inputrec_strings *is)
{

    int         i, nRead;
    int         ninp;
    t_inpfile  *inp;
    const char *tmp;
    char       *ptr1[MAXPTR];

    ninp   = *ninp_p;
    inp    = *inp_p;

    snew(is->dens_grp, STRLEN);
    CTYPE ("Index group of atoms used for spreading");
    STYPE ("spread-group", is->dens_grp, "");
    /* We need at least one value for k and sigma. But if the user wishes
     * to do multiple refinement steps in an automated way, additional
     * values can be added. Interpolation between these values will be
     * done like in the simulated annealing case. */
    CTYPE("Number of time points to use for specifying the refinement with different values for k, sigma, and temperature");
    ITYPE("densfit-npoints", npoints_, 1);
    if (npoints_ < 1 || npoints_ >= MAXPTR)
    {
        gmx_fatal(FARGS, "The number of refinement points should be at least 1 and at most %d\n", MAXPTR);
    }
    snew(time_values_, npoints_);
    snew(sigma_values_, npoints_);
    snew(k_values_, npoints_);
    snew(temp_values_, npoints_);

    CTYPE ("List of times (ps) at which density fitting / refinement parameters are given and exactly matched,");
    CTYPE ("interpolate linearly in between");
    STYPE ("densfit-time", is->densfit_time, NULL);
    CTYPE ("Gaussian width(s) used to transform the discrete atomic positions into a density (one sigma per refinement point)");
    STYPE ("densfit-sigma", is->densfit_sigma, NULL);
    CTYPE ("Strength constant(s) k (kJ/mol) of map potential V_fit = k*(1 - corr.coeff.) (one k per refinement point)");
    STYPE ("densfit-k", is->densfit_k, NULL);
    CTYPE ("Temperature(s) (K) (one T per refinement point) - leave empty if you do not want to change the temperature!");
    STYPE ("densfit-temp", is->densfit_temp,  NULL);

    CTYPE ("Cutoff distance (in sigmas) for density spreading, spread only within +/- sigma*dist");
    RTYPE ("dist", dist_, 4.0);
    CTYPE ("Calculate V_fit every nstfit time steps");
    ITYPE ("nstfit", nstfit_, 1);
    CTYPE ("Write diagnostic output every nstout time steps");
    ITYPE ("nstout", nstout_, 1);
    CTYPE ("Write simulated density map to file every nstmapout time steps");
    ITYPE ("nstmapout", nstmapout_, 1000);
    CTYPE ("Keep and number simulated density maps (results in a bunch of map files)");
    EETYPE("KeepMaps", bKeepAndNumberMaps_, yesno_names);

    /* Get values for the refinement points from input strings */
    /* Read time values from input string*/
    nRead = str_nelem(is->densfit_time, MAXPTR, ptr1);
    if (nRead != npoints_)
    {
        gmx_fatal(FARGS, "Found %d time values for density fitting refinement,\n"
                  "but expected %d (as given in densfit-npoints)\n", nRead, npoints_);
    }
    for (i = 0; i < npoints_; i++)
    {
        time_values_[i] = strtod(ptr1[i], NULL);
    }
    for (i = 1; i < npoints_; i++)
    {
        if (time_values_[i] <= time_values_[i-1])
        {
            gmx_fatal(FARGS, "Density fitting time points must be given in increasing order!\n");
        }
    }

    /* Read sigma values from input string */
    nRead = str_nelem(is->densfit_sigma, MAXPTR, ptr1);
    if (nRead != npoints_)
    {
        gmx_fatal(FARGS, "Found %d sigma values for density fitting refinement,\n"
                  "but expected %d (as given in densfit-npoints)\n", nRead, npoints_);
    }
    for (i = 0; i < npoints_; i++)
    {
        sigma_values_[i] = strtod(ptr1[i], NULL);
    }

    /* Read k values from input string */
    nRead = str_nelem(is->densfit_k, MAXPTR, ptr1);
    if (nRead != npoints_)
    {
        gmx_fatal(FARGS, "Found %d k values for density fitting refinement,\n"
                  "but expected %d (as given in densfit-npoints)\n", nRead, npoints_);
    }
    for (i = 0; i < npoints_; i++)
    {
        k_values_[i] = strtod(ptr1[i], NULL);
        if (k_values_[i] < 0)
        {
            warning(wi, "Got a negative value for a density fitting k value. This will expel atoms from the map!");
        }
    }

    /* Read temperature values from string */
    nRead = str_nelem(is->densfit_temp, MAXPTR, ptr1);
    if (0 == nRead)     /* A user can omit the temperature specification altogether if T should not be touched */
    {
        for (i = 0; i < npoints_; i++)
        {
            temp_values_[i] = -1;      /* negative value signals: not set */
        }
    }
    else
    {
        if (nRead != npoints_)
        {
            gmx_fatal(FARGS, "Found %d temperature values for density fitting refinement,\n"
                      "but expected %d (as given in densfit-npoints)\n", nRead, npoints_);
        }
        for (i = 0; i < npoints_; i++)
        {
            temp_values_[i] = strtod(ptr1[i], NULL);
            if (temp_values_[i] < 0.0)
            {
                gmx_fatal(FARGS, "Sorry, I don't know what a negative temperature means.\n"
                          "Density fitting temperature point %d is %f K\n", i, temp_values_[i]);
            }
        }
    }

    *ninp_p   = ninp;
    *inp_p    = inp;
}


real DensfitData::dist() const
{
    return dist_;
};

} // namespace gmx
