/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 *  \brief Implementation of the Weighted Histogram Analysis Method (WHAM)
 *
 *  \author Jochen Hub <jhub@gwdg.de>
 */
#include "gmxpre.h"

#include "config.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <sstream>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

//! longest file names allowed in input files
#define WHAM_MAXFILELEN 2048

/*! \brief
 * x-axis legend for output files
 */
static const char *xlabel = "\\xx\\f{} (nm)";

/*! \brief
 * enum for energy units
 */
enum {
    enSel, en_kJ, en_kCal, en_kT, enNr
};
/*! \brief
 * enum for type of input files (pdos, tpr, or pullf)
 */
enum {
    whamin_unknown, whamin_tpr, whamin_pullxf, whamin_pdo
};

/*! \brief enum for bootstrapping method
 *
 * These bootstrap methods are supported:
 *  - bootstrap complete histograms with continuous weights (Bayesian bootstrap)
 *    (bsMethod_BayesianHist)
 *  - bootstrap complete histograms (bsMethod_hist)
 *  - bootstrap trajectories from given umbrella histograms. This generates new
 *    "synthetic" histograms (bsMethod_traj)
 *  - bootstrap trajectories from Gaussian with mu/sigma computed from
 *    the respective histogram (bsMethod_trajGauss). This gives very similar
 *    results compared to bsMethod_traj.
 *
 *  ********************************************************************
 *  FOR MORE DETAILS ON THE BOOTSTRAP METHODS (INCLUDING EXAMPLES), SEE
 *  JS Hub, BL de Groot, D van der Spoel
 *  g_wham - A free weighted histogram analysis implementation including
 *  robust error and autocorrelation estimates,
 *  J Chem Theory Comput, 6(12), 3713-3720 (2010)
 *  ********************************************************************
 */
enum {
    bsMethod_unknown, bsMethod_BayesianHist, bsMethod_hist,
    bsMethod_traj, bsMethod_trajGauss
};


//! Parameters of the umbrella potentials
typedef struct
{
    /*!
     * \name Using umbrella pull code since gromacs 4.x
     */
    /*!\{*/
    int      npullcrds_tot; //!< nr of pull coordinates in tpr file
    int      npullcrds;     //!< nr of umbrella pull coordinates for reading
    int      pull_geometry; //!< such as distance, direction
    ivec     pull_dim;      //!< pull dimension with geometry distance
    int      pull_ndim;     //!< nr of pull_dim != 0
    gmx_bool bPrintRef;     //!< Coordinates of reference groups written to pullx.xvg ?
    gmx_bool bPrintComp;    //!< Components of pull distance written to pullx.xvg ?
    real    *k;             //!< force constants in tpr file
    real    *init_dist;     //!< reference displacements
    real    *umbInitDist;   //!< reference displacement in umbrella direction
    /*!\}*/
    /*!
     * \name Using PDO files common until gromacs 3.x
     */
    /*!\{*/
    int    nSkip;
    char   Reference[256];
    int    nPull;
    int    nDim;
    ivec   Dims;
    char   PullName[4][256];
    double UmbPos[4][3];
    double UmbCons[4][3];
    /*!\}*/
} t_UmbrellaHeader;

//! Data in the umbrella histograms
typedef struct
{
    int      nPull;       //!< nr of pull groups in this pdo or pullf/x file
    double **Histo;       //!< nPull histograms
    double **cum;         //!< nPull cumulative distribution functions
    int      nBin;        //!< nr of bins. identical to opt->bins
    double  *k;           //!< force constants for the nPull groups
    double  *pos;         //!< umbrella positions for the nPull groups
    double  *z;           //!< z=(-Fi/kT) for the nPull groups. These values are iteratively computed during wham
    int     *N;           //!< nr of data points in nPull histograms.
    int     *Ntot;        //!< also nr of data points. N and Ntot only differ if bHistEq==TRUE

    /*! \brief  g = 1 + 2*tau[int]/dt where tau is the integrated autocorrelation time.
     *
     * Compare, e.g. Ferrenberg/Swendsen, PRL 63:1195 (1989),
     * Kumar et al, J Comp Chem 13, 1011-1021 (1992), eq. 28
     */
    double *g;
    double *tau;         //!< intetrated autocorrelation time (IACT)
    double *tausmooth;   //!< smoothed IACT

    double  dt;          //!< timestep in the input data. Can be adapted with g_wham option -dt

    /*! \brief TRUE, if any data point of the histogram is within min and max, otherwise FALSE */
    gmx_bool **bContrib;
    real     **ztime;     //!< input data z(t) as a function of time. Required to compute ACTs

    /*! \brief average force estimated from average displacement, fAv=dzAv*k
     *
     *  Used for integration to guess the potential.
     */
    real   *forceAv;
    real   *aver;         //!< average of histograms
    real   *sigma;        //!< stddev of histograms
    double *bsWeight;     //!< for bootstrapping complete histograms with continuous weights
} t_UmbrellaWindow;

//! Selection of pull groups to be used in WHAM (one structure for each tpr file)
typedef struct
{
    int       n;         //!< total nr of pull groups in this tpr file
    int       nUse;      //!< nr of pull groups used
    gmx_bool *bUse;      //!< boolean array of size n. =1 if used, =0 if not
} t_groupselection;

//! Parameters of WHAM
typedef struct
{
    /*!
     * \name Input stuff
     */
    /*!\{*/
    const char       *fnTpr, *fnPullf, *fnGroupsel;
    const char       *fnPdo, *fnPullx;            //!< file names of input
    gmx_bool          bTpr, bPullf, bPdo, bPullx; //!< input file types given?
    real              tmin, tmax, dt;             //!< only read input within tmin and tmax with dt

    gmx_bool          bInitPotByIntegration;      //!< before WHAM, guess potential by force integration. Yields 1.5 to 2 times faster convergence
    int               stepUpdateContrib;          //!< update contribution table every ... iterations. Accelerates WHAM.
    int               nGroupsel;                  //!< if >0: use only certain group in WHAM, if ==0: use all groups
    t_groupselection *groupsel;                   //!< for each tpr file: which pull groups to use in WHAM?
    /*!\}*/
    /*!
     * \name Basic WHAM options
     */
    /*!\{*/
    int      bins;                   //!< nr of bins, min, max, and dz of profile
    real     min, max, dz;
    real     Temperature, Tolerance; //!< temperature, converged when probability changes less than Tolerance
    gmx_bool bCycl;                  //!< generate cyclic (periodic) PMF
    /*!\}*/
    /*!
     * \name Output control
     */
    /*!\{*/
    gmx_bool bLog;                   //!< energy output (instead of probability) for profile
    int      unit;                   //!< unit for PMF output kJ/mol or kT or kCal/mol
    gmx_bool bSym;                   //!< symmetrize PMF around z=0 after WHAM, useful for membranes etc.
    /*! \brief after wham, set prof to zero at this z-position.
     * When bootstrapping, set zProf0 to a "stable" reference position.
     */
    real         zProf0;
    gmx_bool     bProf0Set;              //!< setting profile to 0 at zProf0?

    gmx_bool     bBoundsOnly, bHistOnly; //!< determine min and max, or write histograms and exit
    gmx_bool     bAuto;                  //!< determine min and max automatically but do not exit

    gmx_bool     verbose;                //!< more noisy wham mode
    int          stepchange;             //!< print maximum change in prof after how many interations
    output_env_t oenv;                   //!< xvgr options
    /*!\}*/
    /*!
     * \name Autocorrelation stuff
     */
    /*!\{*/
    gmx_bool bTauIntGiven, bCalcTauInt; //!< IACT given or should be calculated?
    real     sigSmoothIact;             //!< sigma of Gaussian to smooth ACTs
    gmx_bool bAllowReduceIact;          //!< Allow to reduce ACTs during smoothing. Otherwise ACT are only increased during smoothing
    real     acTrestart;                //!< when computing ACT, time between restarting points

    /* \brief Enforce the same weight for each umbella window, that is
     *  calculate with the same number of data points for
     *  each window. That can be reasonable, if the histograms
     *  have different length, but due to autocorrelation,
     *  a longer simulation should not have larger weightin wham.
     */
    gmx_bool bHistEq;
    /*!\}*/

    /*!
     * \name Bootstrapping stuff
     */
    /*!\{*/
    int nBootStrap;              //!< nr of bootstraps (50 is usually enough)

    /* \brief bootstrap method
     *
     * if == bsMethod_hist, consider complete histograms as independent
     * data points and, hence, only mix complete histograms.
     * if == bsMethod_BayesianHist, consider complete histograms
     * as independent data points, but assign random weights
     * to the histograms during the bootstrapping ("Bayesian bootstrap")
     * In case of long correlations (e.g., inside a channel), these
     * will yield a more realistic error.
     * if == bsMethod_traj(Gauss), generate synthetic histograms
     * for each given
     * histogram by generating an autocorrelated random sequence
     * that is distributed according to the respective given
     * histogram. With bsMethod_trajGauss, bootstrap from a Gaussian
     * (instead of from the umbrella histogram) to generate a new
     * histogram.
     */
    int bsMethod;

    /* \brief  autocorrelation time (ACT) used to generate synthetic histograms. If ==0, use calculated ACF */
    real tauBootStrap;

    /* \brief when mixing histograms, mix only histograms withing blocks
              long the reaction coordinate xi. Avoids gaps along xi. */
    int histBootStrapBlockLength;

    int bsSeed;                    //!< random seed for bootstrapping

    /* \brief Write cumulative distribution functions (CDFs) of histograms
              and write the generated histograms for each bootstrap */
    gmx_bool bs_verbose;
    /*!\}*/
    /*!
     * \name tabulated umbrella potential stuff
     */
    /*!\{*/
    gmx_bool  bTab;
    double   *tabX, *tabY, tabMin, tabMax, tabDz;
    int       tabNbins;
    /*!\}*/
    gmx_rng_t rng;                  //!< gromacs random number generator
} t_UmbrellaOptions;

//! Make an umbrella window (may contain several histograms)
t_UmbrellaWindow * initUmbrellaWindows(int nwin)
{
    t_UmbrellaWindow *win;
    int               i;
    snew(win, nwin);
    for (i = 0; i < nwin; i++)
    {
        win[i].Histo    = win[i].cum  = 0;
        win[i].k        = win[i].pos  = win[i].z = 0;
        win[i].N        = win[i].Ntot = 0;
        win[i].g        = win[i].tau  = win[i].tausmooth = 0;
        win[i].bContrib = 0;
        win[i].ztime    = 0;
        win[i].forceAv  = 0;
        win[i].aver     = win[i].sigma = 0;
        win[i].bsWeight = 0;
    }
    return win;
}

//! Delete an umbrella window (may contain several histograms)
void freeUmbrellaWindows(t_UmbrellaWindow *win, int nwin)
{
    int i, j;
    for (i = 0; i < nwin; i++)
    {
        if (win[i].Histo)
        {
            for (j = 0; j < win[i].nPull; j++)
            {
                sfree(win[i].Histo[j]);
            }
        }
        if (win[i].cum)
        {
            for (j = 0; j < win[i].nPull; j++)
            {
                sfree(win[i].cum[j]);
            }
        }
        if (win[i].bContrib)
        {
            for (j = 0; j < win[i].nPull; j++)
            {
                sfree(win[i].bContrib[j]);
            }
        }
        sfree(win[i].Histo);
        sfree(win[i].cum);
        sfree(win[i].k);
        sfree(win[i].pos);
        sfree(win[i].z);
        sfree(win[i].N);
        sfree(win[i].Ntot);
        sfree(win[i].g);
        sfree(win[i].tau);
        sfree(win[i].tausmooth);
        sfree(win[i].bContrib);
        sfree(win[i].ztime);
        sfree(win[i].forceAv);
        sfree(win[i].aver);
        sfree(win[i].sigma);
        sfree(win[i].bsWeight);
    }
    sfree(win);
}

/*! \brief
 * Read and setup tabulated umbrella potential
 */
void setup_tab(const char *fn, t_UmbrellaOptions *opt)
{
    int      i, ny, nl;
    double **y;

    printf("Setting up tabulated potential from file %s\n", fn);
    nl            = read_xvg(fn, &y, &ny);
    opt->tabNbins = nl;
    if (ny != 2)
    {
        gmx_fatal(FARGS, "Found %d columns in %s. Expected 2.\n", ny, fn);
    }
    opt->tabMin = y[0][0];
    opt->tabMax = y[0][nl-1];
    opt->tabDz  = (opt->tabMax-opt->tabMin)/(nl-1);
    if (opt->tabDz <= 0)
    {
        gmx_fatal(FARGS, "The tabulated potential in %s must be provided in \n"
                  "ascending z-direction", fn);
    }
    for (i = 0; i < nl-1; i++)
    {
        if  (fabs(y[0][i+1]-y[0][i]-opt->tabDz) > opt->tabDz/1e6)
        {
            gmx_fatal(FARGS, "z-values in %s are not equally spaced.\n", ny, fn);
        }
    }
    snew(opt->tabY, nl);
    snew(opt->tabX, nl);
    for (i = 0; i < nl; i++)
    {
        opt->tabX[i] = y[0][i];
        opt->tabY[i] = y[1][i];
    }
    printf("Found equally spaced tabulated potential from %g to %g, spacing %g\n",
           opt->tabMin, opt->tabMax, opt->tabDz);
}

//! Read the header of an PDO file (position, force const, nr of groups)
void read_pdo_header(FILE * file, t_UmbrellaHeader * header, t_UmbrellaOptions *opt)
{
    char               line[2048];
    char               Buffer0[256], Buffer1[256], Buffer2[256], Buffer3[256], Buffer4[256];
    int                i;
    std::istringstream ist;

    /*  line 1 */
    if (fgets(line, 2048, file) == NULL)
    {
        gmx_fatal(FARGS, "Error reading header from pdo file\n");
    }
    ist.str(line);
    ist >> Buffer0 >> Buffer1 >> Buffer2;
    if (strcmp(Buffer1, "UMBRELLA"))
    {
        gmx_fatal(FARGS, "This does not appear to be a valid pdo file. Found %s, expected %s\n"
                  "(Found in first line: `%s')\n",
                  Buffer1, "UMBRELLA", line);
    }
    if (strcmp(Buffer2, "3.0"))
    {
        gmx_fatal(FARGS, "This does not appear to be a version 3.0 pdo file");
    }

    /*  line 2 */
    if (fgets(line, 2048, file) == NULL)
    {
        gmx_fatal(FARGS, "Error reading header from pdo file\n");
    }
    ist.str(line);
    ist >> Buffer0 >> Buffer1 >> Buffer2 >> header->Dims[0] >> header->Dims[1] >> header->Dims[2];
    /* printf("%d %d %d\n", header->Dims[0],header->Dims[1],header->Dims[2]); */

    header->nDim = header->Dims[0] + header->Dims[1] + header->Dims[2];
    if (header->nDim != 1)
    {
        gmx_fatal(FARGS, "Currently only supports one dimension");
    }

    /* line3 */
    if (fgets(line, 2048, file) == NULL)
    {
        gmx_fatal(FARGS, "Error reading header from pdo file\n");
    }
    ist.str(line);
    ist >> Buffer0 >> Buffer1 >> header->nSkip;

    /* line 4 */
    if (fgets(line, 2048, file) == NULL)
    {
        gmx_fatal(FARGS, "Error reading header from pdo file\n");
    }
    ist.str(line);
    ist >> Buffer0 >> Buffer1 >> Buffer2 >> header->Reference;

    /* line 5 */
    if (fgets(line, 2048, file) == NULL)
    {
        gmx_fatal(FARGS, "Error reading header from pdo file\n");
    }
    ist.str(line);
    ist >> Buffer0 >> Buffer1 >> Buffer2 >> Buffer3 >> Buffer4 >> header->nPull;

    if (opt->verbose)
    {
        printf("\tFound nPull=%d , nSkip=%d, ref=%s\n", header->nPull, header->nSkip,
               header->Reference);
    }

    for (i = 0; i < header->nPull; ++i)
    {
        if (fgets(line, 2048, file) == NULL)
        {
            gmx_fatal(FARGS, "Error reading header from pdo file\n");
        }
        ist.str(line);
        ist >> Buffer0 >> Buffer1 >> Buffer2 >> header->PullName[i];
        ist >> Buffer0 >> Buffer1 >> header->UmbPos[i][0];
        ist >> Buffer0 >> Buffer1 >> header->UmbCons[i][0];

        if (opt->verbose)
        {
            printf("\tpullgroup %d, pullname = %s, UmbPos = %g, UmbConst = %g\n",
                   i, header->PullName[i], header->UmbPos[i][0], header->UmbCons[i][0]);
        }
    }

    if (fgets(line, 2048, file) == NULL)
    {
        gmx_fatal(FARGS, "Cannot read from file\n");
    }
    ist.str(line);
    ist >> Buffer3;
    if (strcmp(Buffer3, "#####") != 0)
    {
        gmx_fatal(FARGS, "Expected '#####', found %s. Hick.\n", Buffer3);
    }
}

//! smarter fgets
static char *fgets3(FILE *fp, char ptr[], int *len)
{
    char *p;
    int   slen;

    if (fgets(ptr, *len-1, fp) == NULL)
    {
        return NULL;
    }
    p = ptr;
    while ((strchr(ptr, '\n') == NULL) && (!feof(fp)))
    {
        /* This line is longer than len characters, let's increase len! */
        *len += STRLEN;
        p    += STRLEN;
        srenew(ptr, *len);
        if (fgets(p-1, STRLEN, fp) == NULL)
        {
            break;
        }
    }
    slen = strlen(ptr);
    if (ptr[slen-1] == '\n')
    {
        ptr[slen-1] = '\0';
    }

    return ptr;
}

/*! \brief Read the data columns of and PDO file.
 *
 *  TO DO: Get rid of the scanf function to avoid the clang warning.
 *         At the moment, this warning is avoided by hiding the format string
 *         the variable fmtlf.
 */
void read_pdo_data(FILE * file, t_UmbrellaHeader * header,
                   int fileno, t_UmbrellaWindow * win,
                   t_UmbrellaOptions *opt,
                   gmx_bool bGetMinMax, real *mintmp, real *maxtmp)
{
    int                i, inttemp, bins, count, ntot;
    real               min, max, minfound = 1e20, maxfound = -1e20;
    double             temp, time, time0 = 0, dt;
    char              *ptr    = 0;
    t_UmbrellaWindow * window = 0;
    gmx_bool           timeok, dt_ok = 1;
    char              *tmpbuf   = 0, fmt[256], fmtign[256], fmtlf[5] = "%lf";
    int                len      = STRLEN, dstep = 1;
    const int          blocklen = 4096;
    int               *lennow   = 0;

    if (!bGetMinMax)
    {
        bins = opt->bins;
        min  = opt->min;
        max  = opt->max;

        window = win+fileno;
        /* Need to alocate memory and set up structure */
        window->nPull = header->nPull;
        window->nBin  = bins;

        snew(window->Histo, window->nPull);
        snew(window->z, window->nPull);
        snew(window->k, window->nPull);
        snew(window->pos, window->nPull);
        snew(window->N, window->nPull);
        snew(window->Ntot, window->nPull);
        snew(window->g, window->nPull);
        snew(window->bsWeight, window->nPull);

        window->bContrib = 0;

        if (opt->bCalcTauInt)
        {
            snew(window->ztime, window->nPull);
        }
        else
        {
            window->ztime = 0;
        }
        snew(lennow, window->nPull);

        for (i = 0; i < window->nPull; ++i)
        {
            window->z[i]        = 1;
            window->bsWeight[i] = 1.;
            snew(window->Histo[i], bins);
            window->k[i]    = header->UmbCons[i][0];
            window->pos[i]  = header->UmbPos[i][0];
            window->N[i]    = 0;
            window->Ntot[i] = 0;
            window->g[i]    = 1.;
            if (opt->bCalcTauInt)
            {
                window->ztime[i] = 0;
            }
        }

        /* Done with setup */
    }
    else
    {
        minfound = 1e20;
        maxfound = -1e20;
        min      = max = bins = 0; /* Get rid of warnings */
    }

    count = 0;
    snew(tmpbuf, len);
    while ( (ptr = fgets3(file, tmpbuf, &len)) != NULL)
    {
        trim(ptr);

        if (ptr[0] == '#' || strlen(ptr) < 2)
        {
            continue;
        }

        /* Initiate format string */
        fmtign[0] = '\0';
        strcat(fmtign, "%*s");

        sscanf(ptr, fmtlf, &time); /* printf("Time %f\n",time); */
        /* Round time to fs */
        time = 1.0/1000*( static_cast<int> (time*1000+0.5) );

        /* get time step of pdo file */
        if (count == 0)
        {
            time0 = time;
        }
        else if (count == 1)
        {
            dt = time-time0;
            if (opt->dt > 0.0)
            {
                dstep = static_cast<int>(opt->dt/dt+0.5);
                if (dstep == 0)
                {
                    dstep = 1;
                }
            }
            if (!bGetMinMax)
            {
                window->dt = dt*dstep;
            }
        }
        count++;

        dt_ok  = ((count-1)%dstep == 0);
        timeok = (dt_ok && time >= opt->tmin && time <= opt->tmax);
        /* if (opt->verbose)
           printf(" time = %f, (tmin,tmax)=(%e,%e), dt_ok=%d timeok=%d\n",
           time,opt->tmin, opt->tmax, dt_ok,timeok); */

        if (timeok)
        {
            for (i = 0; i < header->nPull; ++i)
            {
                strcpy(fmt, fmtign);
                strcat(fmt, "%lf");      /* Creating a format stings such as "%*s...%*s%lf" */
                strcat(fmtign, "%*s");   /* ignoring one more entry in the next loop */
                if (sscanf(ptr, fmt, &temp))
                {
                    temp += header->UmbPos[i][0];
                    if (bGetMinMax)
                    {
                        if (temp < minfound)
                        {
                            minfound = temp;
                        }
                        if (temp > maxfound)
                        {
                            maxfound = temp;
                        }
                    }
                    else
                    {
                        if (opt->bCalcTauInt)
                        {
                            /* save time series for autocorrelation analysis */
                            ntot = window->Ntot[i];
                            if (ntot >= lennow[i])
                            {
                                lennow[i] += blocklen;
                                srenew(window->ztime[i], lennow[i]);
                            }
                            window->ztime[i][ntot] = temp;
                        }

                        temp -= min;
                        temp /= (max-min);
                        temp *= bins;
                        temp  = floor(temp);

                        inttemp = static_cast<int> (temp);
                        if (opt->bCycl)
                        {
                            if (inttemp < 0)
                            {
                                inttemp += bins;
                            }
                            else if (inttemp >= bins)
                            {
                                inttemp -= bins;
                            }
                        }

                        if (inttemp >= 0 && inttemp < bins)
                        {
                            window->Histo[i][inttemp] += 1.;
                            window->N[i]++;
                        }
                        window->Ntot[i]++;
                    }
                }
            }
        }
        if (time > opt->tmax)
        {
            if (opt->verbose)
            {
                printf("time %f larger than tmax %f, stop reading pdo file\n", time, opt->tmax);
            }
            break;
        }
    }

    if (bGetMinMax)
    {
        *mintmp = minfound;
        *maxtmp = maxfound;
    }

    sfree(lennow);
    sfree(tmpbuf);
}

/*! \brief Set identical weights for all histograms
 *
 * Normally, the weight is given by the number data points in each
 * histogram, together with the autocorrelation time. This can be overwritten
 * by this routine (not recommended). Since we now support autocorrelations, it is better to set
 * an appropriate autocorrelation times instead of using this function.
 */
void enforceEqualWeights(t_UmbrellaWindow * window, int nWindows)
{
    int    i, k, j, NEnforced;
    double ratio;

    NEnforced = window[0].Ntot[0];
    printf("\nFound -hist-eq. Enforcing equal weights for all histograms, \ni.e. doing a "
           "non-weighted histogram analysis method. Ndata = %d\n", NEnforced);
    /* enforce all histograms to have the same weight as the very first histogram */

    for (j = 0; j < nWindows; ++j)
    {
        for (k = 0; k < window[j].nPull; ++k)
        {
            ratio = 1.0*NEnforced/window[j].Ntot[k];
            for (i = 0; i < window[0].nBin; ++i)
            {
                window[j].Histo[k][i] *= ratio;
            }
            window[j].N[k] = static_cast<int>(ratio*window[j].N[k] + 0.5);
        }
    }
}

/*! \brief Simple linear interpolation between two given tabulated points
 */
double tabulated_pot(double dist, t_UmbrellaOptions *opt)
{
    int    jl, ju;
    double pl, pu, dz, dp;

    jl = static_cast<int> (floor((dist-opt->tabMin)/opt->tabDz));
    ju = jl+1;
    if (jl < 0 || ju >= opt->tabNbins)
    {
        gmx_fatal(FARGS, "Distance %f out of bounds of tabulated potential (jl=%d, ju=%d).\n"
                  "Provide an extended table.", dist, jl, ju);
    }
    pl = opt->tabY[jl];
    pu = opt->tabY[ju];
    dz = dist-opt->tabX[jl];
    dp = (pu-pl)*dz/opt->tabDz;
    return pl+dp;
}


/*! \brief
 * Check which bins substiantially contribute (accelerates WHAM)
 *
 * Don't worry, that routine does not mean we compute the PMF in limited precision.
 * After rapid convergence (using only substiantal contributions), we always switch to
 * full precision.
 */
void setup_acc_wham(double *profile, t_UmbrellaWindow * window, int nWindows,
                    t_UmbrellaOptions *opt)
{
    int           i, j, k, nGrptot = 0, nContrib = 0, nTot = 0;
    double        U, min = opt->min, dz = opt->dz, temp, ztot_half, distance, ztot, contrib1, contrib2;
    gmx_bool      bAnyContrib;
    static int    bFirst = 1;
    static double wham_contrib_lim;

    if (bFirst)
    {
        for (i = 0; i < nWindows; ++i)
        {
            nGrptot += window[i].nPull;
        }
        wham_contrib_lim = opt->Tolerance/nGrptot;
    }

    ztot      = opt->max-opt->min;
    ztot_half = ztot/2;

    for (i = 0; i < nWindows; ++i)
    {
        if (!window[i].bContrib)
        {
            snew(window[i].bContrib, window[i].nPull);
        }
        for (j = 0; j < window[i].nPull; ++j)
        {
            if (!window[i].bContrib[j])
            {
                snew(window[i].bContrib[j], opt->bins);
            }
            bAnyContrib = FALSE;
            for (k = 0; k < opt->bins; ++k)
            {
                temp     = (1.0*k+0.5)*dz+min;
                distance = temp - window[i].pos[j];   /* distance to umbrella center */
                if (opt->bCycl)
                {                                     /* in cyclic wham:             */
                    if (distance > ztot_half)         /*    |distance| < ztot_half   */
                    {
                        distance -= ztot;
                    }
                    else if (distance < -ztot_half)
                    {
                        distance += ztot;
                    }
                }
                /* Note: there are two contributions to bin k in the wham equations:
                   i)  N[j]*exp(- U/(8.314e-3*opt->Temperature) + window[i].z[j])
                   ii) exp(- U/(8.314e-3*opt->Temperature))
                   where U is the umbrella potential
                   If any of these number is larger wham_contrib_lim, I set contrib=TRUE
                 */

                if (!opt->bTab)
                {
                    U = 0.5*window[i].k[j]*sqr(distance);       /* harmonic potential assumed. */
                }
                else
                {
                    U = tabulated_pot(distance, opt);            /* Use tabulated potential     */
                }
                contrib1                 = profile[k]*exp(-U/(8.314e-3*opt->Temperature));
                contrib2                 = window[i].N[j]*exp(-U/(8.314e-3*opt->Temperature) + window[i].z[j]);
                window[i].bContrib[j][k] = (contrib1 > wham_contrib_lim || contrib2 > wham_contrib_lim);
                bAnyContrib              = (bAnyContrib | window[i].bContrib[j][k]);
                if (window[i].bContrib[j][k])
                {
                    nContrib++;
                }
                nTot++;
            }
            /* If this histo is far outside min and max all bContrib may be FALSE,
               causing a floating point exception later on. To avoid that, switch
               them all to true.*/
            if (!bAnyContrib)
            {
                for (k = 0; k < opt->bins; ++k)
                {
                    window[i].bContrib[j][k] = TRUE;
                }
            }
        }
    }
    if (bFirst)
    {
        printf("Initialized rapid wham stuff (contrib tolerance %g)\n"
               "Evaluating only %d of %d expressions.\n\n", wham_contrib_lim, nContrib, nTot);
    }

    if (opt->verbose)
    {
        printf("Updated rapid wham stuff. (evaluating only %d of %d contributions)\n",
               nContrib, nTot);
    }
    bFirst = 0;
}

//! Compute the PMF (one of the two main WHAM routines)
void calc_profile(double *profile, t_UmbrellaWindow * window, int nWindows,
                  t_UmbrellaOptions *opt, gmx_bool bExact)
{
    double ztot_half, ztot, min = opt->min, dz = opt->dz;

    ztot      = opt->max-opt->min;
    ztot_half = ztot/2;

#pragma omp parallel
    {
        int nthreads  = gmx_omp_get_max_threads();
        int thread_id = gmx_omp_get_thread_num();
        int i;
        int i0        = thread_id*opt->bins/nthreads;
        int i1        = std::min(opt->bins, ((thread_id+1)*opt->bins)/nthreads);

        for (i = i0; i < i1; ++i)
        {
            int    j, k;
            double num, denom, invg, temp = 0, distance, U = 0;
            num = denom = 0.;
            for (j = 0; j < nWindows; ++j)
            {
                for (k = 0; k < window[j].nPull; ++k)
                {
                    invg = 1.0/window[j].g[k] * window[j].bsWeight[k];
                    temp = (1.0*i+0.5)*dz+min;
                    num += invg*window[j].Histo[k][i];

                    if (!(bExact || window[j].bContrib[k][i]))
                    {
                        continue;
                    }
                    distance = temp - window[j].pos[k];   /* distance to umbrella center */
                    if (opt->bCycl)
                    {                                     /* in cyclic wham:             */
                        if (distance > ztot_half)         /*    |distance| < ztot_half   */
                        {
                            distance -= ztot;
                        }
                        else if (distance < -ztot_half)
                        {
                            distance += ztot;
                        }
                    }

                    if (!opt->bTab)
                    {
                        U = 0.5*window[j].k[k]*sqr(distance);       /* harmonic potential assumed. */
                    }
                    else
                    {
                        U = tabulated_pot(distance, opt);            /* Use tabulated potential     */
                    }
                    denom += invg*window[j].N[k]*exp(-U/(8.314e-3*opt->Temperature) + window[j].z[k]);
                }
            }
            profile[i] = num/denom;
        }
    }
}

//! Compute the free energy offsets z (one of the two main WHAM routines)
double calc_z(double * profile, t_UmbrellaWindow * window, int nWindows,
              t_UmbrellaOptions *opt, gmx_bool bExact)
{
    double min     = opt->min, dz = opt->dz, ztot_half, ztot;
    double maxglob = -1e20;

    ztot      = opt->max-opt->min;
    ztot_half = ztot/2;

#pragma omp parallel
    {
        int    nthreads  = gmx_omp_get_max_threads();
        int    thread_id = gmx_omp_get_thread_num();
        int    i;
        int    i0        = thread_id*nWindows/nthreads;
        int    i1        = std::min(nWindows, ((thread_id+1)*nWindows)/nthreads);
        double maxloc    = -1e20;

        for (i = i0; i < i1; ++i)
        {
            double total     = 0, temp, distance, U = 0;
            int    j, k;

            for (j = 0; j < window[i].nPull; ++j)
            {
                total = 0;
                for (k = 0; k < window[i].nBin; ++k)
                {
                    if (!(bExact || window[i].bContrib[j][k]))
                    {
                        continue;
                    }
                    temp     = (1.0*k+0.5)*dz+min;
                    distance = temp - window[i].pos[j];   /* distance to umbrella center */
                    if (opt->bCycl)
                    {                                     /* in cyclic wham:             */
                        if (distance > ztot_half)         /*    |distance| < ztot_half   */
                        {
                            distance -= ztot;
                        }
                        else if (distance < -ztot_half)
                        {
                            distance += ztot;
                        }
                    }

                    if (!opt->bTab)
                    {
                        U = 0.5*window[i].k[j]*sqr(distance);       /* harmonic potential assumed. */
                    }
                    else
                    {
                        U = tabulated_pot(distance, opt);            /* Use tabulated potential     */
                    }
                    total += profile[k]*exp(-U/(8.314e-3*opt->Temperature));
                }
                /* Avoid floating point exception if window is far outside min and max */
                if (total != 0.0)
                {
                    total = -log(total);
                }
                else
                {
                    total = 1000.0;
                }
                temp = fabs(total - window[i].z[j]);
                if (temp > maxloc)
                {
                    maxloc = temp;
                }
                window[i].z[j] = total;
            }
        }
        /* Now get maximum maxloc from the threads and put in maxglob */
        if (maxloc > maxglob)
        {
#pragma omp critical
            {
                if (maxloc > maxglob)
                {
                    maxglob = maxloc;
                }
            }
        }
    }

    return maxglob;
}

//! Make PMF symmetric around 0 (useful e.g. for membranes)
void symmetrizeProfile(double* profile, t_UmbrellaOptions *opt)
{
    int     i, j, bins = opt->bins;
    double *prof2, min = opt->min, max = opt->max, dz = opt->dz, zsym, deltaz, profsym;
    double  z, z1;

    if (min > 0. || max < 0.)
    {
        gmx_fatal(FARGS, "Cannot symmetrize profile around z=0 with min=%f and max=%f\n",
                  opt->min, opt->max);
    }

    snew(prof2, bins);

    for (i = 0; i < bins; i++)
    {
        z    = min+(i+0.5)*dz;
        zsym = -z;
        /* bin left of zsym */
        j = static_cast<int> (floor((zsym-min)/dz-0.5));
        if (j >= 0 && (j+1) < bins)
        {
            /* interpolate profile linearly between bins j and j+1 */
            z1      = min+(j+0.5)*dz;
            deltaz  = zsym-z1;
            profsym = profile[j] + (profile[j+1]-profile[j])/dz*deltaz;
            /* average between left and right */
            prof2[i] = 0.5*(profsym+profile[i]);
        }
        else
        {
            prof2[i] = profile[i];
        }
    }

    memcpy(profile, prof2, bins*sizeof(double));
    sfree(prof2);
}

//! Set energy unit (kJ/mol,kT,kCal/mol) and set it to zero at opt->zProf0
void prof_normalization_and_unit(double * profile, t_UmbrellaOptions *opt)
{
    int    i, bins, imin;
    double unit_factor = 1., R_MolarGasConst, diff;

    R_MolarGasConst = 8.314472e-3; /* in kJ/(mol*K) */
    bins            = opt->bins;

    /* Not log? Nothing to do! */
    if (!opt->bLog)
    {
        return;
    }

    /* Get profile in units of kT, kJ/mol, or kCal/mol */
    if (opt->unit == en_kT)
    {
        unit_factor = 1.0;
    }
    else if (opt->unit == en_kJ)
    {
        unit_factor = R_MolarGasConst*opt->Temperature;
    }
    else if (opt->unit == en_kCal)
    {
        unit_factor = R_MolarGasConst*opt->Temperature/4.1868;
    }
    else
    {
        gmx_fatal(FARGS, "Sorry, I don't know this energy unit.");
    }

    for (i = 0; i < bins; i++)
    {
        if (profile[i] > 0.0)
        {
            profile[i] = -log(profile[i])*unit_factor;
        }
    }

    /* shift to zero at z=opt->zProf0 */
    if (!opt->bProf0Set)
    {
        diff = profile[0];
    }
    else
    {
        /* Get bin with shortest distance to opt->zProf0
           (-0.5 from bin position and +0.5 from rounding cancel) */
        imin = static_cast<int>((opt->zProf0-opt->min)/opt->dz);
        if (imin < 0)
        {
            imin = 0;
        }
        else if (imin >= bins)
        {
            imin = bins-1;
        }
        diff = profile[imin];
    }

    /* Shift to zero */
    for (i = 0; i < bins; i++)
    {
        profile[i] -= diff;
    }
}

//! Make an array of random integers (used for bootstrapping)
void getRandomIntArray(int nPull, int blockLength, int* randomArray, gmx_rng_t rng)
{
    int ipull, blockBase, nr, ipullRandom;

    if (blockLength == 0)
    {
        blockLength = nPull;
    }

    for (ipull = 0; ipull < nPull; ipull++)
    {
        blockBase = (ipull/blockLength)*blockLength;
        do
        {      /* make sure nothing bad happens in the last block */
            nr          = static_cast<int>(gmx_rng_uniform_real(rng)*blockLength);
            ipullRandom = blockBase + nr;
        }
        while (ipullRandom >= nPull);
        if (ipullRandom < 0 || ipullRandom >= nPull)
        {
            gmx_fatal(FARGS, "Ups, random iWin = %d, nPull = %d, nr = %d, "
                      "blockLength = %d, blockBase = %d\n",
                      ipullRandom, nPull, nr, blockLength, blockBase);
        }
        randomArray[ipull] = ipullRandom;
    }
    /*for (ipull=0; ipull<nPull; ipull++)
       printf("%d ",randomArray[ipull]); printf("\n"); */
}

/*! \brief Set pull group information of a synthetic histogram
 *
 * This is used when bootstapping new trajectories and thereby create new histogtrams,
 * but it is not required if we bootstrap complete histograms.
 */
void copy_pullgrp_to_synthwindow(t_UmbrellaWindow *synthWindow,
                                 t_UmbrellaWindow *thisWindow, int pullid)
{
    synthWindow->N       [0] = thisWindow->N        [pullid];
    synthWindow->Histo   [0] = thisWindow->Histo    [pullid];
    synthWindow->pos     [0] = thisWindow->pos      [pullid];
    synthWindow->z       [0] = thisWindow->z        [pullid];
    synthWindow->k       [0] = thisWindow->k        [pullid];
    synthWindow->bContrib[0] = thisWindow->bContrib [pullid];
    synthWindow->g       [0] = thisWindow->g        [pullid];
    synthWindow->bsWeight[0] = thisWindow->bsWeight [pullid];
}

/*! \brief Calculate cumulative distribution function of of all histograms.
 *
 * This allow to create random number sequences
 * which are distributed according to the histograms. Required to generate
 * the "synthetic" histograms for the Bootstrap method
 */
void calc_cumulatives(t_UmbrellaWindow *window, int nWindows,
                      t_UmbrellaOptions *opt, const char *fnhist)
{
    int    i, j, k, nbin;
    double last;
    char  *fn = 0, *buf = 0;
    FILE  *fp = 0;

    if (opt->bs_verbose)
    {
        snew(fn, strlen(fnhist)+10);
        snew(buf, strlen(fnhist)+10);
        sprintf(fn, "%s_cumul.xvg", strncpy(buf, fnhist, strlen(fnhist)-4));
        fp = xvgropen(fn, "CDFs of umbrella windows", xlabel, "CDF", opt->oenv);
    }

    nbin = opt->bins;
    for (i = 0; i < nWindows; i++)
    {
        snew(window[i].cum, window[i].nPull);
        for (j = 0; j < window[i].nPull; j++)
        {
            snew(window[i].cum[j], nbin+1);
            window[i].cum[j][0] = 0.;
            for (k = 1; k <= nbin; k++)
            {
                window[i].cum[j][k] = window[i].cum[j][k-1]+window[i].Histo[j][k-1];
            }

            /* normalize CDFs. Ensure cum[nbin]==1 */
            last = window[i].cum[j][nbin];
            for (k = 0; k <= nbin; k++)
            {
                window[i].cum[j][k] /= last;
            }
        }
    }

    printf("Cumulative distriubtion functions of all histograms created.\n");
    if (opt->bs_verbose)
    {
        for (k = 0; k <= nbin; k++)
        {
            fprintf(fp, "%g\t", opt->min+k*opt->dz);
            for (i = 0; i < nWindows; i++)
            {
                for (j = 0; j < window[i].nPull; j++)
                {
                    fprintf(fp, "%g\t", window[i].cum[j][k]);
                }
            }
            fprintf(fp, "\n");
        }
        printf("Wrote cumulative distribution functions to %s\n", fn);
        xvgrclose(fp);
        sfree(fn);
        sfree(buf);
    }
}


/*! \brief Return j such that xx[j] <= x < xx[j+1]
 *
 *  This is used to generate a random sequence distributed according to a histogram
 */
void searchCumulative(double xx[], int n, double x, int *j)
{
    int ju, jm, jl;

    jl = -1;
    ju = n;
    while (ju-jl > 1)
    {
        jm = (ju+jl) >> 1;
        if (x >= xx[jm])
        {
            jl = jm;
        }
        else
        {
            ju = jm;
        }
    }
    if (x == xx[0])
    {
        *j = 0;
    }
    else if (x == xx[n-1])
    {
        *j = n-2;
    }
    else
    {
        *j = jl;
    }
}

//! Bootstrap new trajectories and thereby generate new (bootstrapped) histograms
void create_synthetic_histo(t_UmbrellaWindow *synthWindow, t_UmbrellaWindow *thisWindow,
                            int pullid, t_UmbrellaOptions *opt)
{
    int    N, i, nbins, r_index, ibin;
    double r, tausteps = 0.0, a, ap, dt, x, invsqrt2, g, y, sig = 0., z, mu = 0.;
    char   errstr[1024];

    N     = thisWindow->N[pullid];
    dt    = thisWindow->dt;
    nbins = thisWindow->nBin;

    /* tau = autocorrelation time */
    if (opt->tauBootStrap > 0.0)
    {
        tausteps = opt->tauBootStrap/dt;
    }
    else if (opt->bTauIntGiven || opt->bCalcTauInt)
    {
        /* calc tausteps from g=1+2tausteps */
        g        = thisWindow->g[pullid];
        tausteps = (g-1)/2;
    }
    else
    {
        sprintf(errstr,
                "When generating hypothetical trajctories from given umbrella histograms,\n"
                "autocorrelation times (ACTs) are required. Otherwise the statistical error\n"
                "cannot be predicted. You have 3 options:\n"
                "1) Make g_wham estimate the ACTs (options -ac and -acsig).\n"
                "2) Calculate the ACTs by yourself (e.g. with g_analyze) and provide them\n");
        strcat(errstr,
               "   with option -iiact for all umbrella windows.\n"
               "3) If all ACTs are identical and know, you can define them with -bs-tau.\n"
               "   Use option (3) only if you are sure what you're doing, you may severely\n"
               "   underestimate the error if a too small ACT is given.\n");
        gmx_fatal(FARGS, errstr);
    }

    synthWindow->N       [0] = N;
    synthWindow->pos     [0] = thisWindow->pos[pullid];
    synthWindow->z       [0] = thisWindow->z[pullid];
    synthWindow->k       [0] = thisWindow->k[pullid];
    synthWindow->bContrib[0] = thisWindow->bContrib[pullid];
    synthWindow->g       [0] = thisWindow->g       [pullid];
    synthWindow->bsWeight[0] = thisWindow->bsWeight[pullid];

    for (i = 0; i < nbins; i++)
    {
        synthWindow->Histo[0][i] = 0.;
    }

    if (opt->bsMethod == bsMethod_trajGauss)
    {
        sig = thisWindow->sigma [pullid];
        mu  = thisWindow->aver  [pullid];
    }

    /* Genrate autocorrelated Gaussian random variable with autocorrelation time tau
       Use the following:
       If x and y are random numbers from N(0,1) (Gaussian with average 0 and sigma=1),
       then
       z = a*x + sqrt(1-a^2)*y
       is also from N(0,1), and cov(z,x) = a. Thus, by gerenating a sequence
       x' = a*x + sqrt(1-a^2)*y, the sequnce x(t) is from N(0,1) and has an autocorrelation
       function
       C(t) = exp(-t/tau) with tau=-1/ln(a)

       Then, use error function to turn the Gaussian random variable into a uniformly
       distributed one in [0,1]. Eventually, use cumulative distribution function of
       histogram to get random variables distributed according to histogram.
       Note: The ACT of the flat distribution and of the generated histogram is not
       100% exactly tau, but near tau (my test was 3.8 instead of 4).
     */
    a        = exp(-1.0/tausteps);
    ap       = sqrt(1-a*a);
    invsqrt2 = 1./sqrt(2.0);

    /* init random sequence */
    x = gmx_rng_gaussian_table(opt->rng);

    if (opt->bsMethod == bsMethod_traj)
    {
        /* bootstrap points from the umbrella histograms */
        for (i = 0; i < N; i++)
        {
            y = gmx_rng_gaussian_table(opt->rng);
            x = a*x+ap*y;
            /* get flat distribution in [0,1] using cumulative distribution function of Gauusian
               Note: CDF(Gaussian) = 0.5*{1+erf[x/sqrt(2)]}
             */
            r = 0.5*(1+gmx_erf(x*invsqrt2));
            searchCumulative(thisWindow->cum[pullid], nbins+1, r, &r_index);
            synthWindow->Histo[0][r_index] += 1.;
        }
    }
    else if (opt->bsMethod == bsMethod_trajGauss)
    {
        /* bootstrap points from a Gaussian with the same average and sigma
           as the respective umbrella histogram. The idea was, that -given
           limited sampling- the bootstrapped histograms are otherwise biased
           from the limited sampling of the US histos. However, bootstrapping from
           the Gaussian seems to yield a similar estimate. */
        i = 0;
        while (i < N)
        {
            y    = gmx_rng_gaussian_table(opt->rng);
            x    = a*x+ap*y;
            z    = x*sig+mu;
            ibin = static_cast<int> (floor((z-opt->min)/opt->dz));
            if (opt->bCycl)
            {
                if (ibin < 0)
                {
                    while ( (ibin += nbins) < 0)
                    {
                        ;
                    }
                }
                else if (ibin >= nbins)
                {
                    while ( (ibin -= nbins) >= nbins)
                    {
                        ;
                    }
                }
            }

            if (ibin >= 0 && ibin < nbins)
            {
                synthWindow->Histo[0][ibin] += 1.;
                i++;
            }
        }
    }
    else
    {
        gmx_fatal(FARGS, "Unknown bsMethod (id %d). That should not happen.\n", opt->bsMethod);
    }
}

/*! \brief Write all histograms to a file
 *
 * If bs_index>=0, a number is added to the output file name to allow the ouput of all
 * sets of bootstrapped histograms.
 */
void print_histograms(const char *fnhist, t_UmbrellaWindow * window, int nWindows,
                      int bs_index, t_UmbrellaOptions *opt)
{
    char *fn = 0, *buf = 0, title[256];
    FILE *fp;
    int   bins, l, i, j;

    if (bs_index >= 0)
    {
        snew(fn, strlen(fnhist)+10);
        snew(buf, strlen(fnhist)+1);
        sprintf(fn, "%s_bs%d.xvg", strncpy(buf, fnhist, strlen(fnhist)-4), bs_index);
        sprintf(title, "Umbrella histograms. Bootstrap #%d", bs_index);
    }
    else
    {
        fn = gmx_strdup(fnhist);
        strcpy(title, "Umbrella histograms");
    }

    fp   = xvgropen(fn, title, xlabel, "count", opt->oenv);
    bins = opt->bins;

    /* Write histograms */
    for (l = 0; l < bins; ++l)
    {
        fprintf(fp, "%e\t", (double)(l+0.5)*opt->dz+opt->min);
        for (i = 0; i < nWindows; ++i)
        {
            for (j = 0; j < window[i].nPull; ++j)
            {
                fprintf(fp, "%e\t", window[i].Histo[j][l]);
            }
        }
        fprintf(fp, "\n");
    }

    xvgrclose(fp);
    printf("Wrote %s\n", fn);
    if (bs_index >= 0)
    {
        sfree(buf);
    }
    sfree(fn);
}

//! Used for qsort to sort random numbers
int func_wham_is_larger(const void *a, const void *b)
{
    double *aa, *bb;
    aa = (double*)a;
    bb = (double*)b;
    if (*aa < *bb)
    {
        return -1;
    }
    else if (*aa > *bb)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

//! Make random weights for histograms for the Bayesian bootstrap of complete histograms)
void setRandomBsWeights(t_UmbrellaWindow *synthwin, int nAllPull, t_UmbrellaOptions *opt)
{
    int     i;
    double *r;

    snew(r, nAllPull);

    /* generate ordered random numbers between 0 and nAllPull  */
    for (i = 0; i < nAllPull-1; i++)
    {
        r[i] = gmx_rng_uniform_real(opt->rng) * nAllPull;
    }
    qsort((void *)r, nAllPull-1, sizeof(double), &func_wham_is_larger);
    r[nAllPull-1] = 1.0*nAllPull;

    synthwin[0].bsWeight[0] = r[0];
    for (i = 1; i < nAllPull; i++)
    {
        synthwin[i].bsWeight[0] = r[i]-r[i-1];
    }

    /* avoid to have zero weight by adding a tiny value */
    for (i = 0; i < nAllPull; i++)
    {
        if (synthwin[i].bsWeight[0] < 1e-5)
        {
            synthwin[i].bsWeight[0] = 1e-5;
        }
    }

    sfree(r);
}

//! The main bootstrapping routine
void do_bootstrapping(const char *fnres, const char* fnprof, const char *fnhist,
                      char* ylabel, double *profile,
                      t_UmbrellaWindow * window, int nWindows, t_UmbrellaOptions *opt)
{
    t_UmbrellaWindow * synthWindow;
    double            *bsProfile, *bsProfiles_av, *bsProfiles_av2, maxchange = 1e20, tmp, stddev;
    int                i, j, *randomArray = 0, winid, pullid, ib;
    int                iAllPull, nAllPull, *allPull_winId, *allPull_pullId;
    FILE              *fp;
    gmx_bool           bExact = FALSE;

    /* init random generator */
    if (opt->bsSeed == -1)
    {
        opt->rng = gmx_rng_init(gmx_rng_make_seed());
    }
    else
    {
        opt->rng = gmx_rng_init(opt->bsSeed);
    }

    snew(bsProfile,     opt->bins);
    snew(bsProfiles_av, opt->bins);
    snew(bsProfiles_av2, opt->bins);

    /* Create array of all pull groups. Note that different windows
       may have different nr of pull groups
       First: Get total nr of pull groups */
    nAllPull = 0;
    for (i = 0; i < nWindows; i++)
    {
        nAllPull += window[i].nPull;
    }
    snew(allPull_winId, nAllPull);
    snew(allPull_pullId, nAllPull);
    iAllPull = 0;
    /* Setup one array of all pull groups */
    for (i = 0; i < nWindows; i++)
    {
        for (j = 0; j < window[i].nPull; j++)
        {
            allPull_winId[iAllPull]  = i;
            allPull_pullId[iAllPull] = j;
            iAllPull++;
        }
    }

    /* setup stuff for synthetic windows */
    snew(synthWindow, nAllPull);
    for (i = 0; i < nAllPull; i++)
    {
        synthWindow[i].nPull = 1;
        synthWindow[i].nBin  = opt->bins;
        snew(synthWindow[i].Histo, 1);
        if (opt->bsMethod == bsMethod_traj || opt->bsMethod == bsMethod_trajGauss)
        {
            snew(synthWindow[i].Histo[0], opt->bins);
        }
        snew(synthWindow[i].N, 1);
        snew(synthWindow[i].pos, 1);
        snew(synthWindow[i].z, 1);
        snew(synthWindow[i].k, 1);
        snew(synthWindow[i].bContrib, 1);
        snew(synthWindow[i].g, 1);
        snew(synthWindow[i].bsWeight, 1);
    }

    switch (opt->bsMethod)
    {
        case bsMethod_hist:
            snew(randomArray, nAllPull);
            printf("\n\nWhen computing statistical errors by bootstrapping entire histograms:\n");
            please_cite(stdout, "Hub2006");
            break;
        case bsMethod_BayesianHist:
            /* just copy all histogams into synthWindow array */
            for (i = 0; i < nAllPull; i++)
            {
                winid  = allPull_winId [i];
                pullid = allPull_pullId[i];
                copy_pullgrp_to_synthwindow(synthWindow+i, window+winid, pullid);
            }
            break;
        case bsMethod_traj:
        case bsMethod_trajGauss:
            calc_cumulatives(window, nWindows, opt, fnhist);
            break;
        default:
            gmx_fatal(FARGS, "Unknown bootstrap method. That should not have happened.\n");
    }

    /* do bootstrapping */
    fp = xvgropen(fnprof, "Boot strap profiles", xlabel, ylabel, opt->oenv);
    for (ib = 0; ib < opt->nBootStrap; ib++)
    {
        printf("  *******************************************\n"
               "  ******** Start bootstrap nr %d ************\n"
               "  *******************************************\n", ib+1);

        switch (opt->bsMethod)
        {
            case bsMethod_hist:
                /* bootstrap complete histograms from given histograms */
                getRandomIntArray(nAllPull, opt->histBootStrapBlockLength, randomArray, opt->rng);
                for (i = 0; i < nAllPull; i++)
                {
                    winid  = allPull_winId [randomArray[i]];
                    pullid = allPull_pullId[randomArray[i]];
                    copy_pullgrp_to_synthwindow(synthWindow+i, window+winid, pullid);
                }
                break;
            case bsMethod_BayesianHist:
                /* keep histos, but assign random weights ("Bayesian bootstrap") */
                setRandomBsWeights(synthWindow, nAllPull, opt);
                break;
            case bsMethod_traj:
            case bsMethod_trajGauss:
                /* create new histos from given histos, that is generate new hypothetical
                   trajectories */
                for (i = 0; i < nAllPull; i++)
                {
                    winid  = allPull_winId[i];
                    pullid = allPull_pullId[i];
                    create_synthetic_histo(synthWindow+i, window+winid, pullid, opt);
                }
                break;
        }

        /* write histos in case of verbose output */
        if (opt->bs_verbose)
        {
            print_histograms(fnhist, synthWindow, nAllPull, ib, opt);
        }

        /* do wham */
        i         = 0;
        bExact    = FALSE;
        maxchange = 1e20;
        memcpy(bsProfile, profile, opt->bins*sizeof(double)); /* use profile as guess */
        do
        {
            if ( (i%opt->stepUpdateContrib) == 0)
            {
                setup_acc_wham(bsProfile, synthWindow, nAllPull, opt);
            }
            if (maxchange < opt->Tolerance)
            {
                bExact = TRUE;
            }
            if (((i%opt->stepchange) == 0 || i == 1) && i != 0)
            {
                printf("\t%4d) Maximum change %e\n", i, maxchange);
            }
            calc_profile(bsProfile, synthWindow, nAllPull, opt, bExact);
            i++;
        }
        while ( (maxchange = calc_z(bsProfile, synthWindow, nAllPull, opt, bExact)) > opt->Tolerance || !bExact);
        printf("\tConverged in %d iterations. Final maximum change %g\n", i, maxchange);

        if (opt->bLog)
        {
            prof_normalization_and_unit(bsProfile, opt);
        }

        /* symmetrize profile around z=0 */
        if (opt->bSym)
        {
            symmetrizeProfile(bsProfile, opt);
        }

        /* save stuff to get average and stddev */
        for (i = 0; i < opt->bins; i++)
        {
            tmp                = bsProfile[i];
            bsProfiles_av[i]  += tmp;
            bsProfiles_av2[i] += tmp*tmp;
            fprintf(fp, "%e\t%e\n", (i+0.5)*opt->dz+opt->min, tmp);
        }
        fprintf(fp, "%s\n", output_env_get_print_xvgr_codes(opt->oenv) ? "&" : "");
    }
    xvgrclose(fp);

    /* write average and stddev */
    fp = xvgropen(fnres, "Average and stddev from bootstrapping", xlabel, ylabel, opt->oenv);
    if (output_env_get_print_xvgr_codes(opt->oenv))
    {
        fprintf(fp, "@TYPE xydy\n");
    }
    for (i = 0; i < opt->bins; i++)
    {
        bsProfiles_av [i] /= opt->nBootStrap;
        bsProfiles_av2[i] /= opt->nBootStrap;
        tmp                = bsProfiles_av2[i]-sqr(bsProfiles_av[i]);
        stddev             = (tmp >= 0.) ? sqrt(tmp) : 0.; /* Catch rouding errors */
        fprintf(fp, "%e\t%e\t%e\n", (i+0.5)*opt->dz+opt->min, bsProfiles_av [i], stddev);
    }
    xvgrclose(fp);
    printf("Wrote boot strap result to %s\n", fnres);
}

//! Return type of input file based on file extension (xvg, pdo, or tpr)
int whaminFileType(char *fn)
{
    int len;
    len = strlen(fn);
    if (strcmp(fn+len-3, "tpr") == 0)
    {
        return whamin_tpr;
    }
    else if (strcmp(fn+len-3, "xvg") == 0 || strcmp(fn+len-6, "xvg.gz") == 0)
    {
        return whamin_pullxf;
    }
    else if (strcmp(fn+len-3, "pdo") == 0 || strcmp(fn+len-6, "pdo.gz") == 0)
    {
        return whamin_pdo;
    }
    else
    {
        gmx_fatal(FARGS, "Unknown file type of %s. Should be tpr, xvg, or pdo.\n", fn);
    }
    return whamin_unknown;
}

//! Read the files names in pdo-files.dat, pullf/x-files.dat, tpr-files.dat
void read_wham_in(const char *fn, char ***filenamesRet, int *nfilesRet,
                  t_UmbrellaOptions *opt)
{
    char **filename = 0, tmp[WHAM_MAXFILELEN+2];
    int    nread, sizenow, i, block = 1;
    FILE  *fp;

    fp      = gmx_ffopen(fn, "r");
    nread   = 0;
    sizenow = 0;
    while (fgets(tmp, sizeof(tmp), fp) != NULL)
    {
        if (strlen(tmp) >= WHAM_MAXFILELEN)
        {
            gmx_fatal(FARGS, "Filename too long in %s. Only %d characters allowed.\n", fn, WHAM_MAXFILELEN);
        }
        if (nread >= sizenow)
        {
            sizenow += block;
            srenew(filename, sizenow);
            for (i = sizenow-block; i < sizenow; i++)
            {
                snew(filename[i], WHAM_MAXFILELEN);
            }
        }
        /* remove newline if there is one */
        if (tmp[strlen(tmp)-1] == '\n')
        {
            tmp[strlen(tmp)-1] = '\0';
        }
        strcpy(filename[nread], tmp);
        if (opt->verbose)
        {
            printf("Found file %s in %s\n", filename[nread], fn);
        }
        nread++;
    }
    *filenamesRet = filename;
    *nfilesRet    = nread;
}

//! Open a file or a pipe to a gzipped file
FILE *open_pdo_pipe(const char *fn, t_UmbrellaOptions *opt, gmx_bool *bPipeOpen)
{
    char            Buffer[1024], gunzip[1024], *Path = 0;
    FILE           *pipe   = 0;
    static gmx_bool bFirst = 1;

    /* gzipped pdo file? */
    if ((strcmp(fn+strlen(fn)-3, ".gz") == 0))
    {
        /* search gunzip executable */
        if (!(Path = getenv("GMX_PATH_GZIP")))
        {
            if (gmx_fexist("/bin/gunzip"))
            {
                sprintf(gunzip, "%s", "/bin/gunzip");
            }
            else if (gmx_fexist("/usr/bin/gunzip"))
            {
                sprintf(gunzip, "%s", "/usr/bin/gunzip");
            }
            else
            {
                gmx_fatal(FARGS, "Cannot find executable gunzip in /bin or /usr/bin.\n"
                          "You may want to define the path to gunzip "
                          "with the environment variable GMX_PATH_GZIP.", gunzip);
            }
        }
        else
        {
            sprintf(gunzip, "%s/gunzip", Path);
            if (!gmx_fexist(gunzip))
            {
                gmx_fatal(FARGS, "Cannot find executable %s. Please define the path to gunzip"
                          " in the environmental varialbe GMX_PATH_GZIP.", gunzip);
            }
        }
        if (bFirst)
        {
            printf("Using gunzig executable %s\n", gunzip);
            bFirst = 0;
        }
        if (!gmx_fexist(fn))
        {
            gmx_fatal(FARGS, "File %s does not exist.\n", fn);
        }
        sprintf(Buffer, "%s -c < %s", gunzip, fn);
        if (opt->verbose)
        {
            printf("Executing command '%s'\n", Buffer);
        }
#ifdef HAVE_PIPES
        if ((pipe = popen(Buffer, "r")) == NULL)
        {
            gmx_fatal(FARGS, "Unable to open pipe to `%s'\n", Buffer);
        }
#else
        gmx_fatal(FARGS, "Cannot open a compressed file on platform without pipe support");
#endif
        *bPipeOpen = TRUE;
    }
    else
    {
        pipe       = gmx_ffopen(fn, "r");
        *bPipeOpen = FALSE;
    }

    return pipe;
}

//! Close file or pipe
void pdo_close_file(FILE *fp)
{
#ifdef HAVE_PIPES
    pclose(fp);
#else
    gmx_ffclose(fp);
#endif
}

//! Reading all pdo files
void read_pdo_files(char **fn, int nfiles, t_UmbrellaHeader* header,
                    t_UmbrellaWindow *window, t_UmbrellaOptions *opt)
{
    FILE    *file;
    real     mintmp, maxtmp, done = 0.;
    int      i;
    gmx_bool bPipeOpen;
    /* char Buffer0[1000]; */

    if (nfiles < 1)
    {
        gmx_fatal(FARGS, "No files found. Hick.");
    }

    /* if min and max are not given, get min and max from the input files */
    if (opt->bAuto)
    {
        printf("Automatic determination of boundaries from %d pdo files...\n", nfiles);
        opt->min = 1e20;
        opt->max = -1e20;
        for (i = 0; i < nfiles; ++i)
        {
            file = open_pdo_pipe(fn[i], opt, &bPipeOpen);
            /*fgets(Buffer0,999,file);
               fprintf(stderr,"First line '%s'\n",Buffer0); */
            done = 100.0*(i+1)/nfiles;
            printf("\rOpening %s ... [%2.0f%%]", fn[i], done); fflush(stdout);
            if (opt->verbose)
            {
                printf("\n");
            }
            read_pdo_header(file, header, opt);
            /* here only determine min and max of this window */
            read_pdo_data(file, header, i, NULL, opt, TRUE, &mintmp, &maxtmp);
            if (maxtmp > opt->max)
            {
                opt->max = maxtmp;
            }
            if (mintmp < opt->min)
            {
                opt->min = mintmp;
            }
            if (bPipeOpen)
            {
                pdo_close_file(file);
            }
            else
            {
                gmx_ffclose(file);
            }
        }
        printf("\n");
        printf("\nDetermined boundaries to %f and %f\n\n", opt->min, opt->max);
        if (opt->bBoundsOnly)
        {
            printf("Found option -boundsonly, now exiting.\n");
            exit (0);
        }
    }
    /* store stepsize in profile */
    opt->dz = (opt->max-opt->min)/opt->bins;

    /* Having min and max, we read in all files */
    /* Loop over all files */
    for (i = 0; i < nfiles; ++i)
    {
        done = 100.0*(i+1)/nfiles;
        printf("\rOpening %s ... [%2.0f%%]", fn[i], done); fflush(stdout);
        if (opt->verbose)
        {
            printf("\n");
        }
        file = open_pdo_pipe(fn[i], opt, &bPipeOpen);
        read_pdo_header(file, header, opt);
        /* load data into window */
        read_pdo_data(file, header, i, window, opt, FALSE, NULL, NULL);
        if ((window+i)->Ntot[0] == 0)
        {
            fprintf(stderr, "\nWARNING, no data points read from file %s (check -b option)\n", fn[i]);
        }
        if (bPipeOpen)
        {
            pdo_close_file(file);
        }
        else
        {
            gmx_ffclose(file);
        }
    }
    printf("\n");
    for (i = 0; i < nfiles; ++i)
    {
        sfree(fn[i]);
    }
    sfree(fn);
}

//! translate 0/1 to N/Y to write pull dimensions
#define int2YN(a) (((a) == 0) ? ("N") : ("Y"))

//! Read pull groups from a tpr file (including position, force const, geometry, number of groups)
void read_tpr_header(const char *fn, t_UmbrellaHeader* header, t_UmbrellaOptions *opt)
{
    t_inputrec  ir;
    int         i, ncrd;
    t_state     state;
    static int  first = 1;

    /* printf("Reading %s \n",fn); */
    read_tpx_state(fn, &ir, &state, NULL, NULL);

    if (!ir.bPull)
    {
        gmx_fatal(FARGS, "This is not a tpr with COM pulling");
    }

    /* nr of pull groups */
    ncrd = 0;
    for (i = 0; i < ir.pull->ncoord; i++)
    {
        if (ir.pull->coord[i].eType == epullUMBRELLA)
        {
            int m;

            if (ncrd == 0)
            {
                header->pull_geometry = ir.pull->coord[i].eGeom;
                copy_ivec(ir.pull->coord[i].dim, header->pull_dim);
            }

            if (ncrd != i)
            {
                /* TODO: remove this restriction */
                gmx_fatal(FARGS, "Currently tpr files where a non-umbrella pull coordinate is followed by an umbrella pull coordinate are not supported");
            }

            for (m = 0; m < DIM; m++)
            {
                if (ir.pull->coord[i].eGeom != header->pull_geometry)
                {
                    /* TODO: remove the restriction */
                    gmx_fatal(FARGS, "Currently all umbrella pull coordinates should use the same geometry");
                }

                if (ir.pull->coord[i].dim[m] != header->pull_dim[m])
                {
                    /* TODO: remove the restriction */
                    gmx_fatal(FARGS, "Currently all umbrella pull coordinates should operate on the same dimensions");
                }
            }

            ncrd++;
        }
    }

    if (ncrd < 1)
    {
        gmx_fatal(FARGS, "This is not a tpr of umbrella simulation. Found only %d umbrella pull coordinates\n", ncrd);
    }

    header->npullcrds_tot = ir.pull->ncoord;
    header->npullcrds     = ncrd;
    header->bPrintRef     = ir.pull->bPrintCOM1;
    if (ir.pull->bPrintCOM2)
    {
        gmx_fatal(FARGS, "Reading pull output with printing of the COM of group 2 is currently not supported");
    }
    if (ir.pull->bPrintRefValue)
    {
        gmx_fatal(FARGS, "Reading pull output with printing of the reference value of the coordinates is currently not supported");
    }
    header->bPrintComp    = ir.pull->bPrintComp;
    header->pull_ndim     = header->pull_dim[0]+header->pull_dim[1]+header->pull_dim[2];
    /* We should add a struct for each pull coord with all pull coord data
       instead of allocating many arrays for each property */
    snew(header->k, ncrd);
    snew(header->init_dist, ncrd);
    snew(header->umbInitDist, ncrd);

    /* only z-direction with epullgCYL? */
    if (header->pull_geometry == epullgCYL)
    {
        if (header->pull_dim[XX] || header->pull_dim[YY] || (!header->pull_dim[ZZ]))
        {
            gmx_fatal(FARGS, "With pull geometry 'cylinder', expected pulling in Z direction only.\n"
                      "However, found dimensions [%s %s %s]\n",
                      int2YN(header->pull_dim[XX]), int2YN(header->pull_dim[YY]),
                      int2YN(header->pull_dim[ZZ]));
        }
    }

    for (i = 0; i < ncrd; i++)
    {
        header->k[i] = ir.pull->coord[i].k;
        if (header->k[i] == 0.0)
        {
            gmx_fatal(FARGS, "Pull coordinate %d has force constant of of 0.0 in %s.\n"
                      "That doesn't seem to be an Umbrella tpr.\n",
                      i, fn);
        }
        header->init_dist[i] =  ir.pull->coord[i].init;

        /* initial distance to reference */
        switch (header->pull_geometry)
        {
            case epullgCYL:
            /* umbrella distance stored in init_dist[i] for geometry cylinder (not in ...[i][ZZ]) */
            case epullgDIST:
            case epullgDIR:
            case epullgDIRPBC:
                header->umbInitDist[i] = header->init_dist[i];
                break;
            default:
                gmx_fatal(FARGS, "Pull geometry %s not supported\n", epullg_names[header->pull_geometry]);
        }
    }

    if (opt->verbose || first)
    {
        printf("File %s, %d coordinates, geometry \"%s\", dimensions [%s %s %s], (%d dimensions)\n"
               "\tPull group coordinates%s expected in pullx files.\n",
               fn, header->npullcrds, epullg_names[header->pull_geometry],
               int2YN(header->pull_dim[0]), int2YN(header->pull_dim[1]), int2YN(header->pull_dim[2]),
               header->pull_ndim, (header->bPrintRef ? "" : " not"));
        for (i = 0; i < ncrd; i++)
        {
            printf("\tcrd %d) k = %-5g  position = %g\n", i, header->k[i], header->umbInitDist[i]);
        }
    }
    if (!opt->verbose && first)
    {
        printf("\tUse option -v to see this output for all input tpr files\n\n");
    }

    first = 0;
}

//! 2-norm in a ndim-dimensional space
double dist_ndim(double **dx, int ndim, int line)
{
    int    i;
    double r2 = 0.;
    for (i = 0; i < ndim; i++)
    {
        r2 += sqr(dx[i][line]);
    }
    return sqrt(r2);
}

//! Read pullx.xvg or pullf.xvg
void read_pull_xf(const char *fn, const char *fntpr, t_UmbrellaHeader * header,
                  t_UmbrellaWindow * window,
                  t_UmbrellaOptions *opt,
                  gmx_bool bGetMinMax, real *mintmp, real *maxtmp,
                  t_groupselection *groupsel)
{
    double        **y = 0, pos = 0., t, force, time0 = 0., dt;
    int             ny, nt, bins, ibin, i, g, gUsed, dstep = 1, nColPerCrd, nColRefEachCrd, nColExpect, ntot, column;
    real            min, max, minfound = 1e20, maxfound = -1e20;
    gmx_bool        dt_ok, timeok, bHaveForce;
    const char     *quantity;
    const int       blocklen = 4096;
    int            *lennow   = 0;
    static gmx_bool bFirst   = TRUE;

    /*
     * Data columns in pull output:
     *  - in force output pullf.xvg:
     *    No reference columns, one column per pull coordinate
     *
     *  - in position output pullx.xvg
     *    bPrintRef == TRUE:  for each pull coordinate: ndim reference columns, and ndim dx columns
     *    bPrintRef == FALSE: for each pull coordinate: no   reference columns, but ndim dx columns
     */

    nColPerCrd = 1;
    if (opt->bPullx && header->bPrintComp)
    {
        nColPerCrd += header->pull_ndim;
    }
    quantity   = opt->bPullx ? "position" : "force";

    if (opt->bPullx && header->bPrintRef)
    {
        nColRefEachCrd = header->pull_ndim;
    }
    else
    {
        nColRefEachCrd = 0;
    }

    nColExpect = 1 + header->npullcrds*(nColRefEachCrd+nColPerCrd);
    bHaveForce = opt->bPullf;

    /* With geometry "distance" or "distance_periodic", only force reading is supported so far.
       That avoids the somewhat tedious extraction of the reaction coordinate from the pullx files.
       Sorry for the laziness, this is a To-do. */
    if  ( (header->pull_geometry == epullgDIR || header->pull_geometry == epullgDIRPBC)
          && opt->bPullx)
    {
        gmx_fatal(FARGS, "With pull geometries \"direction\" and \"direction_periodic\", only pull force "
                  "reading \n(option -if) is supported at present, "
                  "not pull position reading (options -ix).\nMake sure mdrun writes the pull "
                  "forces (pullf.xvg files)\nand provide them to g_wham with option -if.",
                  epullg_names[header->pull_geometry]);
    }

    nt = read_xvg(fn, &y, &ny);

    /* Check consistency */
    if (nt < 1)
    {
        gmx_fatal(FARGS, "Empty pull %s file %s\n", quantity, fn);
    }
    if (bFirst)
    {
        printf("Reading pull %s file with pull geometry %s and %d pull dimensions\n",
               bHaveForce ? "force" : "position", epullg_names[header->pull_geometry],
               header->pull_ndim);
        /* Since this tool code has not updated yet to allow difference pull
         * dimensions per pull coordinate, we can't easily check the exact
         * number of expected columns, so we only check what we expect for
         * the pull coordinates selected for the WHAM analysis.
         */
        printf("Expecting these columns in pull file:\n"
               "\t%d reference columns for each individual pull coordinate\n"
               "\t%d data columns for each pull coordinate\n", nColRefEachCrd, nColPerCrd);
        printf("With %d pull groups, expect %s%d columns (including the time column)\n",
               header->npullcrds,
               header->npullcrds < header->npullcrds_tot ? "at least " : "",
               nColExpect);
        bFirst = FALSE;
    }
    if (ny < nColExpect ||
        (header->npullcrds == header->npullcrds_tot && ny > nColExpect))
    {
        gmx_fatal(FARGS, "Using %d pull coodinates from %s,\n but found %d data columns in %s (expected %s%d)\n"
                  "\nMaybe you confused options -ix and -if ?\n",
                  header->npullcrds, fntpr, ny-1, fn,
                  header->npullcrds < header->npullcrds_tot ? "at least " : "",
                  nColExpect-1);
    }

    if (opt->verbose)
    {
        printf("Found %d times and %d %s sets %s\n", nt, (ny-1)/nColPerCrd, quantity, fn);
    }

    if (!bGetMinMax)
    {
        bins = opt->bins;
        min  = opt->min;
        max  = opt->max;
        if (nt > 1)
        {
            window->dt = y[0][1]-y[0][0];
        }
        else if (opt->nBootStrap && opt->tauBootStrap != 0.0)
        {
            fprintf(stderr, "\n *** WARNING, Could not determine time step in %s\n", fn);
        }

        /* Need to alocate memory and set up structure */

        if (groupsel)
        {
            /* Use only groups selected with option -is file */
            if (header->npullcrds != groupsel->n)
            {
                gmx_fatal(FARGS, "tpr file contains %d pull groups, but expected %d from group selection file\n",
                          header->npullcrds, groupsel->n);
            }
            window->nPull = groupsel->nUse;
        }
        else
        {
            window->nPull = header->npullcrds;
        }

        window->nBin = bins;
        snew(window->Histo, window->nPull);
        snew(window->z, window->nPull);
        snew(window->k, window->nPull);
        snew(window->pos, window->nPull);
        snew(window->N, window->nPull);
        snew(window->Ntot, window->nPull);
        snew(window->g, window->nPull);
        snew(window->bsWeight, window->nPull);
        window->bContrib = 0;

        if (opt->bCalcTauInt)
        {
            snew(window->ztime, window->nPull);
        }
        else
        {
            window->ztime = NULL;
        }
        snew(lennow, window->nPull);

        for (g = 0; g < window->nPull; ++g)
        {
            window->z[g]        = 1;
            window->bsWeight[g] = 1.;
            snew(window->Histo[g], bins);
            window->N[g]    = 0;
            window->Ntot[g] = 0;
            window->g[g]    = 1.;
            if (opt->bCalcTauInt)
            {
                window->ztime[g] = NULL;
            }
        }

        /* Copying umbrella center and force const is more involved since not
           all pull groups from header (tpr file) may be used in window variable */
        for (g = 0, gUsed = 0; g < header->npullcrds; ++g)
        {
            if (groupsel && (groupsel->bUse[g] == FALSE))
            {
                continue;
            }
            window->k[gUsed]   = header->k[g];
            window->pos[gUsed] = header->umbInitDist[g];
            gUsed++;
        }
    }
    else
    {   /* only determine min and max */
        minfound = 1e20;
        maxfound = -1e20;
        min      = max = bins = 0; /* Get rid of warnings */
    }


    for (i = 0; i < nt; i++)
    {
        /* Do you want that time frame? */
        t = 1.0/1000*( static_cast<int> ((y[0][i]*1000) + 0.5)); /* round time to fs */

        /* get time step of pdo file and get dstep from opt->dt */
        if (i == 0)
        {
            time0 = t;
        }
        else if (i == 1)
        {
            dt = t-time0;
            if (opt->dt > 0.0)
            {
                dstep = static_cast<int>(opt->dt/dt+0.5);
                if (dstep == 0)
                {
                    dstep = 1;
                }
            }
            if (!bGetMinMax)
            {
                window->dt = dt*dstep;
            }
        }

        dt_ok  = (i%dstep == 0);
        timeok = (dt_ok && t >= opt->tmin && t <= opt->tmax);
        /*if (opt->verbose)
           printf(" time = %f, (tmin,tmax)=(%e,%e), dt_ok=%d timeok=%d\n",
           t,opt->tmin, opt->tmax, dt_ok,timeok); */
        if (timeok)
        {
            /* Note: if groupsel == NULL:
             *          all groups in pullf/x file are stored in this window, and gUsed == g
             *       if groupsel != NULL:
             *          only groups with groupsel.bUse[g]==TRUE are stored. gUsed is not always equal g
             */
            gUsed = -1;
            for (g = 0; g < header->npullcrds; ++g)
            {
                /* was this group selected for application in WHAM? */
                if (groupsel && (groupsel->bUse[g] == FALSE))
                {
                    continue;
                }

                gUsed++;

                if (bHaveForce)
                {
                    /* y has 1 time column y[0] and one column per force y[1],...,y[nCrds] */
                    force = y[g+1][i];
                    pos   = -force/header->k[g] + header->umbInitDist[g];
                }
                else
                {
                    /* pick the right column from:
                     * time ref1[ndim] coord1[1(ndim)] ref2[ndim] coord2[1(ndim)] ...
                     */
                    column = 1 +  nColRefEachCrd + g*(nColRefEachCrd+nColPerCrd);
                    pos    = y[column][i];
                }

                /* printf("crd %d dpos %f poseq %f pos %f \n",g,dpos,poseq,pos); */
                if (bGetMinMax)
                {
                    if (pos < minfound)
                    {
                        minfound = pos;
                    }
                    if (pos > maxfound)
                    {
                        maxfound = pos;
                    }
                }
                else
                {
                    if (gUsed >= window->nPull)
                    {
                        gmx_fatal(FARGS, "gUsed too large (%d, nPull=%d). This error should have been catched before.\n",
                                  gUsed, window->nPull);
                    }

                    if (opt->bCalcTauInt && !bGetMinMax)
                    {
                        /* save time series for autocorrelation analysis */
                        ntot = window->Ntot[gUsed];
                        /* printf("i %d, ntot %d, lennow[g] = %d\n",i,ntot,lennow[g]); */
                        if (ntot >= lennow[gUsed])
                        {
                            lennow[gUsed] += blocklen;
                            srenew(window->ztime[gUsed], lennow[gUsed]);
                        }
                        window->ztime[gUsed][ntot] = pos;
                    }

                    ibin = static_cast<int> (floor((pos-min)/(max-min)*bins));
                    if (opt->bCycl)
                    {
                        if (ibin < 0)
                        {
                            while ( (ibin += bins) < 0)
                            {
                                ;
                            }
                        }
                        else if (ibin >= bins)
                        {
                            while ( (ibin -= bins) >= bins)
                            {
                                ;
                            }
                        }
                    }
                    if (ibin >= 0 && ibin < bins)
                    {
                        window->Histo[gUsed][ibin] += 1.;
                        window->N[gUsed]++;
                    }
                    window->Ntot[gUsed]++;
                }
            }
        }
        else if (t > opt->tmax)
        {
            if (opt->verbose)
            {
                printf("time %f larger than tmax %f, stop reading this pullx/pullf file\n", t, opt->tmax);
            }
            break;
        }
    }

    if (bGetMinMax)
    {
        *mintmp = minfound;
        *maxtmp = maxfound;
    }
    sfree(lennow);
    for (i = 0; i < ny; i++)
    {
        sfree(y[i]);
    }
}

//! read pullf-files.dat or pullx-files.dat and tpr-files.dat
void read_tpr_pullxf_files(char **fnTprs, char **fnPull, int nfiles,
                           t_UmbrellaHeader* header,
                           t_UmbrellaWindow *window, t_UmbrellaOptions *opt)
{
    int  i;
    real mintmp, maxtmp;

    printf("Reading %d tpr and pullf files\n", nfiles/2);

    /* min and max not given? */
    if (opt->bAuto)
    {
        printf("Automatic determination of boundaries...\n");
        opt->min = 1e20;
        opt->max = -1e20;
        for (i = 0; i < nfiles; i++)
        {
            if (whaminFileType(fnTprs[i]) != whamin_tpr)
            {
                gmx_fatal(FARGS, "Expected the %d'th file in input file to be a tpr file\n", i);
            }
            read_tpr_header(fnTprs[i], header, opt);
            if (whaminFileType(fnPull[i]) != whamin_pullxf)
            {
                gmx_fatal(FARGS, "Expected the %d'th file in input file to be a xvg (pullx/pullf) file\n", i);
            }
            read_pull_xf(fnPull[i], fnTprs[i], header, NULL, opt, TRUE, &mintmp, &maxtmp,
                         (opt->nGroupsel > 0) ? &opt->groupsel[i] : NULL);
            if (maxtmp > opt->max)
            {
                opt->max = maxtmp;
            }
            if (mintmp < opt->min)
            {
                opt->min = mintmp;
            }
        }
        printf("\nDetermined boundaries to %f and %f\n\n", opt->min, opt->max);
        if (opt->bBoundsOnly)
        {
            printf("Found option -boundsonly, now exiting.\n");
            exit (0);
        }
    }
    /* store stepsize in profile */
    opt->dz = (opt->max-opt->min)/opt->bins;

    for (i = 0; i < nfiles; i++)
    {
        if (whaminFileType(fnTprs[i]) != whamin_tpr)
        {
            gmx_fatal(FARGS, "Expected the %d'th file in input file to be a tpr file\n", i);
        }
        read_tpr_header(fnTprs[i], header, opt);
        if (whaminFileType(fnPull[i]) != whamin_pullxf)
        {
            gmx_fatal(FARGS, "Expected the %d'th file in input file to be a xvg (pullx/pullf) file\n", i);
        }
        read_pull_xf(fnPull[i], fnTprs[i], header, window+i, opt, FALSE, NULL, NULL,
                     (opt->nGroupsel > 0) ? &opt->groupsel[i] : NULL);
        if (window[i].Ntot[0] == 0)
        {
            fprintf(stderr, "\nWARNING, no data points read from file %s (check -b option)\n", fnPull[i]);
        }
    }

    for (i = 0; i < nfiles; i++)
    {
        sfree(fnTprs[i]);
        sfree(fnPull[i]);
    }
    sfree(fnTprs);
    sfree(fnPull);
}

/*! \brief Read integrated autocorrelation times from input file (option -iiact)
 *
 * Note: Here we consider tau[int] := int_0^inf ACF(t) as the integrated autocorrelation times.
 * The factor `g := 1 + 2*tau[int]` subsequently enters the uncertainty.
 */
void readIntegratedAutocorrelationTimes(t_UmbrellaWindow *window, int nwins, const char* fn)
{
    int      nlines, ny, i, ig;
    double **iact;

    printf("Readging Integrated autocorrelation times from %s ...\n", fn);
    nlines = read_xvg(fn, &iact, &ny);
    if (nlines != nwins)
    {
        gmx_fatal(FARGS, "Found %d lines with integrated autocorrelation times in %s.\nExpected %d",
                  nlines, fn, nwins);
    }
    for (i = 0; i < nlines; i++)
    {
        if (window[i].nPull != ny)
        {
            gmx_fatal(FARGS, "You are providing autocorrelation times with option -iiact and the\n"
                      "number of pull groups is different in different simulations. That is not\n"
                      "supported yet. Sorry.\n");
        }
        for (ig = 0; ig < window[i].nPull; ig++)
        {
            /* compare Kumar et al, J Comp Chem 13, 1011-1021 (1992) */
            window[i].g[ig] = 1+2*iact[ig][i]/window[i].dt;

            if (iact[ig][i] <= 0.0)
            {
                fprintf(stderr, "\nWARNING, IACT = %f (window %d, group %d)\n", iact[ig][i], i, ig);
            }
        }
    }
}


/*! \brief Smooth autocorreltion times along the reaction coordinate.
 *
 * This is useful
 * if the ACT is subject to high uncertainty in case if limited sampling. Note
 * that -in case of limited sampling- the ACT may be severely underestimated.
 * Note: the g=1+2tau are overwritten.
 * If opt->bAllowReduceIact==FALSE, the ACTs are never reduced, only increased
 * by the smoothing
 */
void smoothIact(t_UmbrellaWindow *window, int nwins, t_UmbrellaOptions *opt)
{
    int    i, ig, j, jg;
    double pos, dpos2, siglim, siglim2, gaufact, invtwosig2, w, weight, tausm;

    /* only evaluate within +- 3sigma of the Gausian */
    siglim  = 3.0*opt->sigSmoothIact;
    siglim2 = dsqr(siglim);
    /* pre-factor of Gaussian */
    gaufact    = 1.0/(sqrt(2*M_PI)*opt->sigSmoothIact);
    invtwosig2 = 0.5/dsqr(opt->sigSmoothIact);

    for (i = 0; i < nwins; i++)
    {
        snew(window[i].tausmooth, window[i].nPull);
        for (ig = 0; ig < window[i].nPull; ig++)
        {
            tausm  = 0.;
            weight = 0;
            pos    = window[i].pos[ig];
            for (j = 0; j < nwins; j++)
            {
                for (jg = 0; jg < window[j].nPull; jg++)
                {
                    dpos2 = dsqr(window[j].pos[jg]-pos);
                    if (dpos2 < siglim2)
                    {
                        w       = gaufact*exp(-dpos2*invtwosig2);
                        weight += w;
                        tausm  += w*window[j].tau[jg];
                        /*printf("Weight %g dpos2=%g pos=%g gaufact=%g invtwosig2=%g\n",
                           w,dpos2,pos,gaufact,invtwosig2); */
                    }
                }
            }
            tausm /= weight;
            if (opt->bAllowReduceIact || tausm > window[i].tau[ig])
            {
                window[i].tausmooth[ig] = tausm;
            }
            else
            {
                window[i].tausmooth[ig] = window[i].tau[ig];
            }
            window[i].g[ig] = 1+2*tausm/window[i].dt;
        }
    }
}

//! Stop integrating autoccorelation function when ACF drops under this value
#define WHAM_AC_ZERO_LIMIT 0.05

/*! \brief Try to compute the autocorrelation time for each umbrealla window
 */
void calcIntegratedAutocorrelationTimes(t_UmbrellaWindow *window, int nwins,
                                        t_UmbrellaOptions *opt, const char *fn)
{
    int   i, ig, ncorr, ntot, j, k, *count, restart;
    real *corr, c0, dt, tmp;
    real *ztime, av, tausteps;
    FILE *fp, *fpcorr = 0;

    if (opt->verbose)
    {
        fpcorr = xvgropen("hist_autocorr.xvg", "Autocorrelation functions of umbrella windows",
                          "time [ps]", "autocorrelation function", opt->oenv);
    }

    printf("\n");
    for (i = 0; i < nwins; i++)
    {
        printf("\rEstimating integrated autocorreltion times ... [%2.0f%%] ...", 100.*(i+1)/nwins);
        fflush(stdout);
        ntot = window[i].Ntot[0];

        /* using half the maximum time as length of autocorrelation function */
        ncorr = ntot/2;
        if (ntot < 10)
        {
            gmx_fatal(FARGS, "Tryig to estimtate autocorrelation time from only %d"
                      " points. Provide more pull data!", ntot);
        }
        snew(corr, ncorr);
        /* snew(corrSq,ncorr); */
        snew(count, ncorr);
        dt = window[i].dt;
        snew(window[i].tau, window[i].nPull);
        restart = static_cast<int>(opt->acTrestart/dt+0.5);
        if (restart == 0)
        {
            restart = 1;
        }

        for (ig = 0; ig < window[i].nPull; ig++)
        {
            if (ntot != window[i].Ntot[ig])
            {
                gmx_fatal(FARGS, "Encountered different nr of frames in different pull groups.\n"
                          "That should not happen. (%d and %d)\n", ntot, window[i].Ntot[ig]);
            }
            ztime = window[i].ztime[ig];

            /* calc autocorrelation function C(t) = < [z(tau)-<z>]*[z(tau+t)-<z>]> */
            for (j = 0, av = 0; (j < ntot); j++)
            {
                av += ztime[j];
            }
            av /= ntot;
            for (k = 0; (k < ncorr); k++)
            {
                corr[k]  = 0.;
                count[k] = 0;
            }
            for (j = 0; (j < ntot); j += restart)
            {
                for (k = 0; (k < ncorr) && (j+k < ntot); k++)
                {
                    tmp        = (ztime[j]-av)*(ztime[j+k]-av);
                    corr  [k] += tmp;
                    /* corrSq[k] += tmp*tmp; */
                    count[k]++;
                }
            }
            /* divide by nr of frames for each time displacement */
            for (k = 0; (k < ncorr); k++)
            {
                /* count probably = (ncorr-k+(restart-1))/restart; */
                corr[k] = corr[k]/count[k];
                /* variance of autocorrelation function */
                /* corrSq[k]=corrSq[k]/count[k]; */
            }
            /* normalize such that corr[0] == 0 */
            c0 = 1./corr[0];
            for (k = 0; (k < ncorr); k++)
            {
                corr[k] *= c0;
                /* corrSq[k]*=c0*c0; */
            }

            /* write ACFs in verbose mode */
            if (fpcorr)
            {
                for (k = 0; (k < ncorr); k++)
                {
                    fprintf(fpcorr, "%g  %g\n", k*dt, corr[k]);
                }
                fprintf(fpcorr, "%s\n", output_env_get_print_xvgr_codes(opt->oenv) ? "&" : "");
            }

            /* esimate integrated correlation time, fitting is too unstable */
            tausteps = 0.5*corr[0];
            /* consider corr below WHAM_AC_ZERO_LIMIT as noise */
            for (j = 1; (j < ncorr) && (corr[j] > WHAM_AC_ZERO_LIMIT); j++)
            {
                tausteps += corr[j];
            }

            /* g = 1+2*tau, see. Ferrenberg/Swendsen, PRL 63:1195 (1989) or
               Kumar et al, eq. 28 ff. */
            window[i].tau[ig] = tausteps*dt;
            window[i].g[ig]   = 1+2*tausteps;
            /* printf("win %d, group %d, estimated correlation time = %g ps\n",i,ig,window[i].tau[ig]); */
        } /* ig loop */
        sfree(corr);
        sfree(count);
    }
    printf(" done\n");
    if (fpcorr)
    {
        xvgrclose(fpcorr);
    }

    /* plot IACT along reaction coordinate */
    fp = xvgropen(fn, "Integrated autocorrelation times", xlabel, "IACT [ps]", opt->oenv);
    if (output_env_get_print_xvgr_codes(opt->oenv))
    {
        fprintf(fp, "@    s0 symbol 1\n@    s0 symbol size 0.5\n@    s0 line linestyle 0\n");
        fprintf(fp, "#  WIN   tau(gr1)  tau(gr2) ...\n");
        for (i = 0; i < nwins; i++)
        {
            fprintf(fp, "# %3d   ", i);
            for (ig = 0; ig < window[i].nPull; ig++)
            {
                fprintf(fp, " %11g", window[i].tau[ig]);
            }
            fprintf(fp, "\n");
        }
    }
    for (i = 0; i < nwins; i++)
    {
        for (ig = 0; ig < window[i].nPull; ig++)
        {
            fprintf(fp, "%8g %8g\n", window[i].pos[ig], window[i].tau[ig]);
        }
    }
    if (opt->sigSmoothIact > 0.0)
    {
        printf("Smoothing autocorrelation times along reaction coordinate with Gaussian of sig = %g\n",
               opt->sigSmoothIact);
        /* smooth IACT along reaction coordinate and overwrite g=1+2tau */
        smoothIact(window, nwins, opt);
        fprintf(fp, "%s\n", output_env_get_print_xvgr_codes(opt->oenv) ? "&" : "");
        if (output_env_get_print_xvgr_codes(opt->oenv))
        {
            fprintf(fp, "@    s1 symbol 1\n@    s1 symbol size 0.5\n@    s1 line linestyle 0\n");
            fprintf(fp, "@    s1 symbol color 2\n");
        }
        for (i = 0; i < nwins; i++)
        {
            for (ig = 0; ig < window[i].nPull; ig++)
            {
                fprintf(fp, "%8g %8g\n", window[i].pos[ig], window[i].tausmooth[ig]);
            }
        }
    }
    xvgrclose(fp);
    printf("Wrote %s\n", fn);
}

/*! \brief
 * compute average and sigma of each umbrella histogram
 */
void averageSigma(t_UmbrellaWindow *window, int nwins)
{
    int  i, ig, ntot, k;
    real av, sum2, sig, diff, *ztime, nSamplesIndep;

    for (i = 0; i < nwins; i++)
    {
        snew(window[i].aver, window[i].nPull);
        snew(window[i].sigma, window[i].nPull);

        ntot = window[i].Ntot[0];
        for (ig = 0; ig < window[i].nPull; ig++)
        {
            ztime = window[i].ztime[ig];
            for (k = 0, av = 0.; k < ntot; k++)
            {
                av += ztime[k];
            }
            av /= ntot;
            for (k = 0, sum2 = 0.; k < ntot; k++)
            {
                diff  = ztime[k]-av;
                sum2 += diff*diff;
            }
            sig                = sqrt(sum2/ntot);
            window[i].aver[ig] = av;

            /* Note: This estimate for sigma is biased from the limited sampling.
               Correct sigma by n/(n-1) where n = number of independent
               samples. Only possible if IACT is known.
             */
            if (window[i].tau)
            {
                nSamplesIndep       = window[i].N[ig]/(window[i].tau[ig]/window[i].dt);
                window[i].sigma[ig] = sig * nSamplesIndep/(nSamplesIndep-1);
            }
            else
            {
                window[i].sigma[ig] = sig;
            }
            printf("win %d, aver = %f  sig = %f\n", i, av, window[i].sigma[ig]);
        }
    }
}


/*! \brief
 * Use histograms to compute average force on pull group.
 */
void computeAverageForce(t_UmbrellaWindow *window, int nWindows, t_UmbrellaOptions *opt)
{
    int    i, j, bins = opt->bins, k;
    double dz, min = opt->min, max = opt->max, displAv, temp, distance, ztot, ztot_half, w, weight;
    double posmirrored;

    dz        = (max-min)/bins;
    ztot      = opt->max-min;
    ztot_half = ztot/2;

    /* Compute average displacement from histograms */
    for (j = 0; j < nWindows; ++j)
    {
        snew(window[j].forceAv, window[j].nPull);
        for (k = 0; k < window[j].nPull; ++k)
        {
            displAv = 0.0;
            weight  = 0.0;
            for (i = 0; i < opt->bins; ++i)
            {
                temp     = (1.0*i+0.5)*dz+min;
                distance = temp - window[j].pos[k];
                if (opt->bCycl)
                {                                       /* in cyclic wham:             */
                    if (distance > ztot_half)           /*    |distance| < ztot_half   */
                    {
                        distance -= ztot;
                    }
                    else if (distance < -ztot_half)
                    {
                        distance += ztot;
                    }
                }
                w         = window[j].Histo[k][i]/window[j].g[k];
                displAv  += w*distance;
                weight   += w;
                /* Are we near min or max? We are getting wrong forces from the histgrams since
                   the histograms are zero outside [min,max). Therefore, assume that the position
                   on the other side of the histomgram center is equally likely. */
                if (!opt->bCycl)
                {
                    posmirrored = window[j].pos[k]-distance;
                    if (posmirrored >= max || posmirrored < min)
                    {
                        displAv  += -w*distance;
                        weight   += w;
                    }
                }
            }
            displAv  /= weight;

            /* average force from average displacement */
            window[j].forceAv[k] = displAv*window[j].k[k];
            /* sigma from average square displacement */
            /* window[j].sigma  [k] = sqrt(displAv2); */
            /* printf("Win %d, sigma = %f\n",j,sqrt(displAv2)); */
        }
    }
}

/*! \brief
 * Check if the complete reaction coordinate is covered by the histograms
 */
void  checkReactionCoordinateCovered(t_UmbrellaWindow *window, int nwins,
                                     t_UmbrellaOptions *opt)
{
    int  i, ig, j, bins = opt->bins, bBoundary;
    real avcount = 0, z, relcount, *count;
    snew(count, opt->bins);

    for (j = 0; j < opt->bins; ++j)
    {
        for (i = 0; i < nwins; i++)
        {
            for (ig = 0; ig < window[i].nPull; ig++)
            {
                count[j] += window[i].Histo[ig][j];
            }
        }
        avcount += 1.0*count[j];
    }
    avcount /= bins;
    for (j = 0; j < bins; ++j)
    {
        relcount  = count[j]/avcount;
        z         = (j+0.5)*opt->dz+opt->min;
        bBoundary = ( j<bins/20 || (bins-j)>bins/20 );
        /* check for bins with no data */
        if (count[j] == 0)
        {
            fprintf(stderr, "\nWARNING, no data point in bin %d (z=%g) !\n"
                    "You may not get a reasonable profile. Check your histograms!\n", j, z);
        }
        /* and check for poor sampling */
        else if (relcount < 0.005 && !bBoundary)
        {
            fprintf(stderr, "Warning, poor sampling bin %d (z=%g). Check your histograms!\n", j, z);
        }
    }
    sfree(count);
}

/*! \brief Compute initial potential by integrating the average force
 *
 * This speeds up the convergence by roughly a factor of 2
 */
void guessPotByIntegration(t_UmbrellaWindow *window, int nWindows, t_UmbrellaOptions *opt)
{
    int    i, j, ig, bins = opt->bins, nHist, winmin, groupmin;
    double dz, min = opt->min, *pot, pos, hispos, dist, diff, fAv, distmin, *f;
    FILE  *fp;

    dz = (opt->max-min)/bins;

    printf("Getting initial potential by integration.\n");

    /* Compute average displacement from histograms */
    computeAverageForce(window, nWindows, opt);

    /* Get force for each bin from all histograms in this bin, or, alternatively,
       if no histograms are inside this bin, from the closest histogram */
    snew(pot, bins);
    snew(f, bins);
    for (j = 0; j < opt->bins; ++j)
    {
        pos      = (1.0*j+0.5)*dz+min;
        nHist    = 0;
        fAv      = 0.;
        distmin  = 1e20;
        groupmin = winmin = 0;
        for (i = 0; i < nWindows; i++)
        {
            for (ig = 0; ig < window[i].nPull; ig++)
            {
                hispos = window[i].pos[ig];
                dist   = fabs(hispos-pos);
                /* average force within bin */
                if (dist < dz/2)
                {
                    nHist++;
                    fAv += window[i].forceAv[ig];
                }
                /* at the same time, remember closest histogram */
                if (dist < distmin)
                {
                    winmin   = i;
                    groupmin = ig;
                    distmin  = dist;
                }
            }
        }
        /* if no histogram found in this bin, use closest histogram */
        if (nHist > 0)
        {
            fAv = fAv/nHist;
        }
        else
        {
            fAv = window[winmin].forceAv[groupmin];
        }
        f[j] = fAv;
    }
    for (j = 1; j < opt->bins; ++j)
    {
        pot[j] = pot[j-1] - 0.5*dz*(f[j-1]+f[j]);
    }

    /* cyclic wham: linearly correct possible offset */
    if (opt->bCycl)
    {
        diff = (pot[bins-1]-pot[0])/(bins-1);
        for (j = 1; j < opt->bins; ++j)
        {
            pot[j] -= j*diff;
        }
    }
    if (opt->verbose)
    {
        fp = xvgropen("pmfintegrated.xvg", "PMF from force integration", xlabel, "PMF (kJ/mol)", opt->oenv);
        for (j = 0; j < opt->bins; ++j)
        {
            fprintf(fp, "%g  %g\n", (j+0.5)*dz+opt->min, pot[j]);
        }
        xvgrclose(fp);
        printf("verbose mode: wrote %s with PMF from interated forces\n", "pmfintegrated.xvg");
    }

    /* get initial z=exp(-F[i]/kT) from integrated potential, where F[i] denote the free
       energy offsets which are usually determined by wham
       First: turn pot into probabilities:
     */
    for (j = 0; j < opt->bins; ++j)
    {
        pot[j] = exp(-pot[j]/(8.314e-3*opt->Temperature));
    }
    calc_z(pot, window, nWindows, opt, TRUE);

    sfree(pot);
    sfree(f);
}

//! Count number of words in an line
static int wordcount(char *ptr)
{
    int i, n, is[2];
    int cur = 0;

    if (strlen(ptr) == 0)
    {
        return 0;
    }
    /* fprintf(stderr,"ptr='%s'\n",ptr); */
    n = 1;
    for (i = 0; (ptr[i] != '\0'); i++)
    {
        is[cur] = isspace(ptr[i]);
        if ((i > 0)  && (is[cur] && !is[1-cur]))
        {
            n++;
        }
        cur = 1-cur;
    }
    return n;
}

/*! \brief Read input file for pull group selection (option -is)
 *
 * TO DO: ptr=fgets(...) is never freed (small memory leak)
 */
void readPullGroupSelection(t_UmbrellaOptions *opt, char **fnTpr, int nTpr)
{
    FILE *fp;
    int   i, iline, n, len = STRLEN, temp;
    char *ptr = 0, *tmpbuf = 0;
    char  fmt[1024], fmtign[1024];
    int   block = 1, sizenow;

    fp            = gmx_ffopen(opt->fnGroupsel, "r");
    opt->groupsel = NULL;

    snew(tmpbuf, len);
    sizenow = 0;
    iline   = 0;
    while ( (ptr = fgets3(fp, tmpbuf, &len)) != NULL)
    {
        trim(ptr);
        n = wordcount(ptr);

        if (iline >= sizenow)
        {
            sizenow += block;
            srenew(opt->groupsel, sizenow);
        }
        opt->groupsel[iline].n    = n;
        opt->groupsel[iline].nUse = 0;
        snew(opt->groupsel[iline].bUse, n);

        fmtign[0] = '\0';
        for (i = 0; i < n; i++)
        {
            strcpy(fmt, fmtign);
            strcat(fmt, "%d");
            if (sscanf(ptr, fmt, &temp))
            {
                opt->groupsel[iline].bUse[i] = (temp > 0);
                if (opt->groupsel[iline].bUse[i])
                {
                    opt->groupsel[iline].nUse++;
                }
            }
            strcat(fmtign, "%*s");
        }
        iline++;
    }
    opt->nGroupsel = iline;
    if (nTpr != opt->nGroupsel)
    {
        gmx_fatal(FARGS, "Found %d tpr files but %d lines in %s\n", nTpr, opt->nGroupsel,
                  opt->fnGroupsel);
    }

    printf("\nUse only these pull groups:\n");
    for (iline = 0; iline < nTpr; iline++)
    {
        printf("%s (%d of %d groups):", fnTpr[iline], opt->groupsel[iline].nUse, opt->groupsel[iline].n);
        for (i = 0; i < opt->groupsel[iline].n; i++)
        {
            if (opt->groupsel[iline].bUse[i])
            {
                printf(" %d", i+1);
            }
        }
        printf("\n");
    }
    printf("\n");

    sfree(tmpbuf);
}

//! Boolean XOR
#define WHAMBOOLXOR(a, b) ( ((!(a)) && (b)) || ((a) && (!(b))))

//! Number of elements in fnm (used for command line parsing)
#define NFILE asize(fnm)

//! The main g_wham routine
int gmx_wham(int argc, char *argv[])
{
    const char              *desc[] = {
        "[THISMODULE] is an analysis program that implements the Weighted",
        "Histogram Analysis Method (WHAM). It is intended to analyze",
        "output files generated by umbrella sampling simulations to ",
        "compute a potential of mean force (PMF).[PAR]",
        "",
        "[THISMODULE] is currently not fully up to date. It only supports pull setups",
        "where the first pull coordinate(s) is/are umbrella pull coordinates",
        "and, if multiple coordinates need to be analyzed, all used the same",
        "geometry and dimensions. In most cases this is not an issue.[PAR]",
        "",
        "At present, three input modes are supported.",
        "",
        "* With option [TT]-it[tt], the user provides a file which contains the",
        "  file names of the umbrella simulation run-input files ([REF].tpr[ref] files),",
        "  AND, with option [TT]-ix[tt], a file which contains file names of",
        "  the pullx [TT]mdrun[tt] output files. The [REF].tpr[ref] and pullx files must",
        "  be in corresponding order, i.e. the first [REF].tpr[ref] created the",
        "  first pullx, etc.",
        "* Same as the previous input mode, except that the the user",
        "  provides the pull force output file names ([TT]pullf.xvg[tt]) with option [TT]-if[tt].",
        "  From the pull force the position in the umbrella potential is",
        "  computed. This does not work with tabulated umbrella potentials.",
        "* With option [TT]-ip[tt], the user provides file names of (gzipped) [REF].pdo[ref] files, i.e.",
        "  the GROMACS 3.3 umbrella output files. If you have some unusual"
        "  reaction coordinate you may also generate your own [REF].pdo[ref] files and",
        "  feed them with the [TT]-ip[tt] option into to [THISMODULE]. The [REF].pdo[ref] file header",
        "  must be similar to the following::",
        "",
        "    # UMBRELLA      3.0",
        "    # Component selection: 0 0 1",
        "    # nSkip 1",
        "    # Ref. Group 'TestAtom'",
        "    # Nr. of pull groups 2",
        "    # Group 1 'GR1'  Umb. Pos. 5.0 Umb. Cons. 1000.0",
        "    # Group 2 'GR2'  Umb. Pos. 2.0 Umb. Cons. 500.0",
        "    #####",
        "",
        "  The number of pull groups, umbrella positions, force constants, and names ",
        "  may (of course) differ. Following the header, a time column and ",
        "  a data column for each pull group follows (i.e. the displacement",
        "  with respect to the umbrella center). Up to four pull groups are possible ",
        "  per [REF].pdo[ref] file at present.",
        "",
        "By default, all pull groups found in all pullx/pullf files are used in WHAM. If only ",
        "some of the pull groups should be used, a pull group selection file (option [TT]-is[tt]) can ",
        "be provided. The selection file must contain one line for each tpr file in tpr-files.dat.",
        "Each of these lines must contain one digit (0 or 1) for each pull group in the tpr file. ",
        "Here, 1 indicates that the pull group is used in WHAM, and 0 means it is omitted. Example:",
        "If you have three tpr files, each containing 4 pull groups, but only pull group 1 and 2 should be ",
        "used, groupsel.dat looks like this::",
        "",
        "    1 1 0 0",
        "    1 1 0 0",
        "    1 1 0 0",
        "",
        "By default, the output files are",
        "",
        "* [TT]-o[tt]      PMF output file",
        "* [TT]-hist[tt]   Histograms output file",
        "",
        "Always check whether the histograms sufficiently overlap.[PAR]",
        "The umbrella potential is assumed to be harmonic and the force constants are ",
        "read from the [REF].tpr[ref] or [REF].pdo[ref] files. If a non-harmonic umbrella force was applied ",
        "a tabulated potential can be provided with [TT]-tab[tt].",
        "",
        "WHAM options",
        "^^^^^^^^^^^^",
        "",
        "* [TT]-bins[tt]   Number of bins used in analysis",
        "* [TT]-temp[tt]   Temperature in the simulations",
        "* [TT]-tol[tt]    Stop iteration if profile (probability) changed less than tolerance",
        "* [TT]-auto[tt]   Automatic determination of boundaries",
        "* [TT]-min,-max[tt]   Boundaries of the profile",
        "",
        "The data points that are used to compute the profile",
        "can be restricted with options [TT]-b[tt], [TT]-e[tt], and [TT]-dt[tt]. ",
        "Adjust [TT]-b[tt] to ensure sufficient equilibration in each ",
        "umbrella window.[PAR]",
        "With [TT]-log[tt] (default) the profile is written in energy units, otherwise ",
        "(with [TT]-nolog[tt]) as probability. The unit can be specified with [TT]-unit[tt]. ",
        "With energy output, the energy in the first bin is defined to be zero. ",
        "If you want the free energy at a different ",
        "position to be zero, set [TT]-zprof0[tt] (useful with bootstrapping, see below).[PAR]",
        "For cyclic or periodic reaction coordinates (dihedral angle, channel PMF",
        "without osmotic gradient), the option [TT]-cycl[tt] is useful.",
        "[THISMODULE] will make use of the",
        "periodicity of the system and generate a periodic PMF. The first and the last bin of the",
        "reaction coordinate will assumed be be neighbors.[PAR]",
        "Option [TT]-sym[tt] symmetrizes the profile around z=0 before output, ",
        "which may be useful for, e.g. membranes.",
        "",
        "Parallelization",
        "^^^^^^^^^^^^^^^",
        "",
        "If available, the number of OpenMP threads used by g_wham is controlled with [TT]-nt[tt].",
        "",
        "Autocorrelations",
        "^^^^^^^^^^^^^^^^",
        "",
        "With [TT]-ac[tt], [THISMODULE] estimates the integrated autocorrelation ",
        "time (IACT) [GRK]tau[grk] for each umbrella window and weights the respective ",
        "window with 1/[1+2*[GRK]tau[grk]/dt]. The IACTs are written ",
        "to the file defined with [TT]-oiact[tt]. In verbose mode, all ",
        "autocorrelation functions (ACFs) are written to [TT]hist_autocorr.xvg[tt]. ",
        "Because the IACTs can be severely underestimated in case of limited ",
        "sampling, option [TT]-acsig[tt] allows one to smooth the IACTs along the ",
        "reaction coordinate with a Gaussian ([GRK]sigma[grk] provided with [TT]-acsig[tt], ",
        "see output in [TT]iact.xvg[tt]). Note that the IACTs are estimated by simple ",
        "integration of the ACFs while the ACFs are larger 0.05.",
        "If you prefer to compute the IACTs by a more sophisticated (but possibly ",
        "less robust) method such as fitting to a double exponential, you can ",
        "compute the IACTs with [gmx-analyze] and provide them to [THISMODULE] with the file ",
        "[TT]iact-in.dat[tt] (option [TT]-iiact[tt]), which should contain one line per ",
        "input file ([REF].pdo[ref] or pullx/f file) and one column per pull group in the respective file.",
        "",
        "Error analysis",
        "^^^^^^^^^^^^^^",
        "",
        "Statistical errors may be estimated with bootstrap analysis. Use it with care, ",
        "otherwise the statistical error may be substantially underestimated. ",
        "More background and examples for the bootstrap technique can be found in ",
        "Hub, de Groot and Van der Spoel, JCTC (2010) 6: 3713-3720.",
        "[TT]-nBootstrap[tt] defines the number of bootstraps (use, e.g., 100). ",
        "Four bootstrapping methods are supported and ",
        "selected with [TT]-bs-method[tt].",
        "",
        "* [TT]b-hist[tt]   Default: complete histograms are considered as independent ",
        "  data points, and the bootstrap is carried out by assigning random weights to the ",
        "  histograms (\"Bayesian bootstrap\"). Note that each point along the reaction coordinate",
        "  must be covered by multiple independent histograms (e.g. 10 histograms), otherwise the ",
        "  statistical error is underestimated.",
        "* [TT]hist[tt]    Complete histograms are considered as independent data points. ",
        "  For each bootstrap, N histograms are randomly chosen from the N given histograms ",
        "  (allowing duplication, i.e. sampling with replacement).",
        "  To avoid gaps without data along the reaction coordinate blocks of histograms ",
        "  ([TT]-histbs-block[tt]) may be defined. In that case, the given histograms are ",
        "  divided into blocks and only histograms within each block are mixed. Note that ",
        "  the histograms within each block must be representative for all possible histograms, ",
        "  otherwise the statistical error is underestimated.",
        "* [TT]traj[tt]  The given histograms are used to generate new random trajectories,",
        "  such that the generated data points are distributed according the given histograms ",
        "  and properly autocorrelated. The autocorrelation time (ACT) for each window must be ",
        "  known, so use [TT]-ac[tt] or provide the ACT with [TT]-iiact[tt]. If the ACT of all ",
        "  windows are identical (and known), you can also provide them with [TT]-bs-tau[tt]. ",
        "  Note that this method may severely underestimate the error in case of limited sampling, ",
        "  that is if individual histograms do not represent the complete phase space at ",
        "  the respective positions.",
        "* [TT]traj-gauss[tt]  The same as method [TT]traj[tt], but the trajectories are ",
        "  not bootstrapped from the umbrella histograms but from Gaussians with the average ",
        "  and width of the umbrella histograms. That method yields similar error estimates ",
        "  like method [TT]traj[tt].",
        "",
        "Bootstrapping output:",
        "",
        "* [TT]-bsres[tt]   Average profile and standard deviations",
        "* [TT]-bsprof[tt]  All bootstrapping profiles",
        "",
        "With [TT]-vbs[tt] (verbose bootstrapping), the histograms of each bootstrap are written, ",
        "and, with bootstrap method [TT]traj[tt], the cumulative distribution functions of ",
        "the histograms."
    };

    const char              *en_unit[]       = {NULL, "kJ", "kCal", "kT", NULL};
    const char              *en_unit_label[] = {"", "E (kJ mol\\S-1\\N)", "E (kcal mol\\S-1\\N)", "E (kT)", NULL};
    const char              *en_bsMethod[]   = { NULL, "b-hist", "hist", "traj", "traj-gauss", NULL };
    static t_UmbrellaOptions opt;
    static int               nthreads = -1;

    t_pargs                  pa[] = {
#ifdef GMX_OPENMP
        { "-nt", FALSE, etINT, {&nthreads},
          "Number of threads used by g_wham (if -1, all threads will be used or what is specified by the environment variable OMP_NUM_THREADS)"},
#endif
        { "-min", FALSE, etREAL, {&opt.min},
          "Minimum coordinate in profile"},
        { "-max", FALSE, etREAL, {&opt.max},
          "Maximum coordinate in profile"},
        { "-auto", FALSE, etBOOL, {&opt.bAuto},
          "Determine min and max automatically"},
        { "-bins", FALSE, etINT, {&opt.bins},
          "Number of bins in profile"},
        { "-temp", FALSE, etREAL, {&opt.Temperature},
          "Temperature"},
        { "-tol", FALSE, etREAL, {&opt.Tolerance},
          "Tolerance"},
        { "-v", FALSE, etBOOL, {&opt.verbose},
          "Verbose mode"},
        { "-b", FALSE, etREAL, {&opt.tmin},
          "First time to analyse (ps)"},
        { "-e", FALSE, etREAL, {&opt.tmax},
          "Last time to analyse (ps)"},
        { "-dt", FALSE, etREAL, {&opt.dt},
          "Analyse only every dt ps"},
        { "-histonly", FALSE, etBOOL, {&opt.bHistOnly},
          "Write histograms and exit"},
        { "-boundsonly", FALSE, etBOOL, {&opt.bBoundsOnly},
          "Determine min and max and exit (with [TT]-auto[tt])"},
        { "-log", FALSE, etBOOL, {&opt.bLog},
          "Calculate the log of the profile before printing"},
        { "-unit", FALSE,  etENUM, {en_unit},
          "Energy unit in case of log output" },
        { "-zprof0", FALSE, etREAL, {&opt.zProf0},
          "Define profile to 0.0 at this position (with [TT]-log[tt])"},
        { "-cycl", FALSE, etBOOL, {&opt.bCycl},
          "Create cyclic/periodic profile. Assumes min and max are the same point."},
        { "-sym", FALSE, etBOOL, {&opt.bSym},
          "Symmetrize profile around z=0"},
        { "-hist-eq", FALSE, etBOOL, {&opt.bHistEq},
          "HIDDENEnforce equal weight for all histograms. (Non-Weighed-HAM)"},
        { "-ac", FALSE, etBOOL, {&opt.bCalcTauInt},
          "Calculate integrated autocorrelation times and use in wham"},
        { "-acsig", FALSE, etREAL, {&opt.sigSmoothIact},
          "Smooth autocorrelation times along reaction coordinate with Gaussian of this [GRK]sigma[grk]"},
        { "-ac-trestart", FALSE, etREAL, {&opt.acTrestart},
          "When computing autocorrelation functions, restart computing every .. (ps)"},
        { "-acred", FALSE, etBOOL, {&opt.bAllowReduceIact},
          "HIDDENWhen smoothing the ACTs, allow to reduce ACTs. Otherwise, only increase ACTs "
          "during smoothing"},
        { "-nBootstrap", FALSE,  etINT, {&opt.nBootStrap},
          "nr of bootstraps to estimate statistical uncertainty (e.g., 200)" },
        { "-bs-method", FALSE,  etENUM, {en_bsMethod},
          "Bootstrap method" },
        { "-bs-tau", FALSE, etREAL, {&opt.tauBootStrap},
          "Autocorrelation time (ACT) assumed for all histograms. Use option [TT]-ac[tt] if ACT is unknown."},
        { "-bs-seed", FALSE, etINT, {&opt.bsSeed},
          "Seed for bootstrapping. (-1 = use time)"},
        { "-histbs-block", FALSE, etINT, {&opt.histBootStrapBlockLength},
          "When mixing histograms only mix within blocks of [TT]-histbs-block[tt]."},
        { "-vbs", FALSE, etBOOL, {&opt.bs_verbose},
          "Verbose bootstrapping. Print the CDFs and a histogram file for each bootstrap."},
        { "-stepout", FALSE, etINT, {&opt.stepchange},
          "HIDDENWrite maximum change every ... (set to 1 with [TT]-v[tt])"},
        { "-updateContr", FALSE, etINT, {&opt.stepUpdateContrib},
          "HIDDENUpdate table with significan contributions to WHAM every ... iterations"},
    };

    t_filenm                 fnm[] = {
        { efDAT, "-ix", "pullx-files", ffOPTRD}, /* wham input: pullf.xvg's and tprs           */
        { efDAT, "-if", "pullf-files", ffOPTRD}, /* wham input: pullf.xvg's and tprs           */
        { efDAT, "-it", "tpr-files", ffOPTRD},   /* wham input: tprs                           */
        { efDAT, "-ip", "pdo-files", ffOPTRD},   /* wham input: pdo files (gmx3 style)         */
        { efDAT, "-is", "groupsel", ffOPTRD},    /* input: select pull groups to use           */
        { efXVG, "-o", "profile", ffWRITE },     /* output file for profile                     */
        { efXVG, "-hist", "histo", ffWRITE},     /* output file for histograms                  */
        { efXVG, "-oiact", "iact", ffOPTWR},     /* writing integrated autocorrelation times    */
        { efDAT, "-iiact", "iact-in", ffOPTRD},  /* reading integrated autocorrelation times   */
        { efXVG, "-bsres", "bsResult", ffOPTWR}, /* average and errors of bootstrap analysis    */
        { efXVG, "-bsprof", "bsProfs", ffOPTWR}, /* output file for bootstrap profiles          */
        { efDAT, "-tab", "umb-pot", ffOPTRD},    /* Tabulated umbrella potential (if not harmonic) */
    };

    int                      i, j, l, nfiles, nwins, nfiles2;
    t_UmbrellaHeader         header;
    t_UmbrellaWindow       * window = NULL;
    double                  *profile, maxchange = 1e20;
    gmx_bool                 bMinSet, bMaxSet, bAutoSet, bExact = FALSE;
    char                   **fninTpr, **fninPull, **fninPdo;
    const char              *fnPull;
    FILE                    *histout, *profout;
    char                     ylabel[256], title[256];

    opt.bins      = 200;
    opt.verbose   = FALSE;
    opt.bHistOnly = FALSE;
    opt.bCycl     = FALSE;
    opt.tmin      = 50;
    opt.tmax      = 1e20;
    opt.dt        = 0.0;
    opt.min       = 0;
    opt.max       = 0;
    opt.bAuto     = TRUE;
    opt.nGroupsel = 0;

    /* bootstrapping stuff */
    opt.nBootStrap               = 0;
    opt.bsMethod                 = bsMethod_hist;
    opt.tauBootStrap             = 0.0;
    opt.bsSeed                   = -1;
    opt.histBootStrapBlockLength = 8;
    opt.bs_verbose               = FALSE;

    opt.bLog                  = TRUE;
    opt.unit                  = en_kJ;
    opt.zProf0                = 0.;
    opt.Temperature           = 298;
    opt.Tolerance             = 1e-6;
    opt.bBoundsOnly           = FALSE;
    opt.bSym                  = FALSE;
    opt.bCalcTauInt           = FALSE;
    opt.sigSmoothIact         = 0.0;
    opt.bAllowReduceIact      = TRUE;
    opt.bInitPotByIntegration = TRUE;
    opt.acTrestart            = 1.0;
    opt.stepchange            = 100;
    opt.stepUpdateContrib     = 100;

    if (!parse_common_args(&argc, argv, 0,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &opt.oenv))
    {
        return 0;
    }

    opt.unit     = nenum(en_unit);
    opt.bsMethod = nenum(en_bsMethod);

    opt.bProf0Set = opt2parg_bSet("-zprof0",  asize(pa), pa);

    opt.bTab         = opt2bSet("-tab", NFILE, fnm);
    opt.bPdo         = opt2bSet("-ip", NFILE, fnm);
    opt.bTpr         = opt2bSet("-it", NFILE, fnm);
    opt.bPullx       = opt2bSet("-ix", NFILE, fnm);
    opt.bPullf       = opt2bSet("-if", NFILE, fnm);
    opt.bTauIntGiven = opt2bSet("-iiact", NFILE, fnm);
    if  (opt.bTab && opt.bPullf)
    {
        gmx_fatal(FARGS, "Force input does not work with tabulated potentials. "
                  "Provide pullx.xvg or pdo files!");
    }

    if (!opt.bPdo && !WHAMBOOLXOR(opt.bPullx, opt.bPullf))
    {
        gmx_fatal(FARGS, "Give either pullx (-ix) OR pullf (-if) data. Not both.");
    }
    if (!opt.bPdo && !(opt.bTpr || opt.bPullf || opt.bPullx))
    {
        gmx_fatal(FARGS, "g_wham supports three input modes, pullx, pullf, or pdo file input."
                  "\n\n Check g_wham -h !");
    }

    opt.fnPdo      = opt2fn("-ip", NFILE, fnm);
    opt.fnTpr      = opt2fn("-it", NFILE, fnm);
    opt.fnPullf    = opt2fn("-if", NFILE, fnm);
    opt.fnPullx    = opt2fn("-ix", NFILE, fnm);
    opt.fnGroupsel = opt2fn_null("-is", NFILE, fnm);

    bMinSet  = opt2parg_bSet("-min",  asize(pa), pa);
    bMaxSet  = opt2parg_bSet("-max",  asize(pa), pa);
    bAutoSet = opt2parg_bSet("-auto",  asize(pa), pa);
    if ( (bMinSet || bMaxSet) && bAutoSet)
    {
        gmx_fatal(FARGS, "With -auto, do not give -min or -max\n");
    }

    if ( (bMinSet && !bMaxSet) || (!bMinSet && bMaxSet))
    {
        gmx_fatal(FARGS, "When giving -min, you must give -max (and vice versa), too\n");
    }

    if (bMinSet && opt.bAuto)
    {
        printf("Note: min and max given, switching off -auto.\n");
        opt.bAuto = FALSE;
    }

    if (opt.bTauIntGiven && opt.bCalcTauInt)
    {
        gmx_fatal(FARGS, "Either read (option -iiact) or calculate (option -ac) the\n"
                  "the autocorrelation times. Not both.");
    }

    if (opt.tauBootStrap > 0.0 && opt2parg_bSet("-ac", asize(pa), pa))
    {
        gmx_fatal(FARGS, "Either compute autocorrelation times (ACTs) (option -ac) or "
                  "provide it with -bs-tau for bootstrapping. Not Both.\n");
    }
    if (opt.tauBootStrap > 0.0 && opt2bSet("-iiact", NFILE, fnm))
    {
        gmx_fatal(FARGS, "Either provide autocorrelation times (ACTs) with file iact-in.dat "
                  "(option -iiact) or define all ACTs with -bs-tau for bootstrapping\n. Not Both.");
    }

    /* Set # of OpenMP threads */
    if (nthreads > 0)
    {
        gmx_omp_set_num_threads(nthreads);
    }
    else
    {
        nthreads = gmx_omp_get_max_threads();
    }
    if (nthreads > 1)
    {
        printf("\nNote: Will use %d OpenMP threads.\n\n", nthreads);
    }

    /* Reading gmx4 pull output and tpr files */
    if (opt.bTpr || opt.bPullf || opt.bPullx)
    {
        read_wham_in(opt.fnTpr, &fninTpr, &nfiles, &opt);

        fnPull = opt.bPullf ? opt.fnPullf : opt.fnPullx;
        read_wham_in(fnPull, &fninPull, &nfiles2, &opt);
        printf("Found %d tpr and %d pull %s files in %s and %s, respectively\n",
               nfiles, nfiles2, opt.bPullf ? "force" : "position", opt.fnTpr, fnPull);
        if (nfiles != nfiles2)
        {
            gmx_fatal(FARGS, "Found %d file names in %s, but %d in %s\n", nfiles,
                      opt.fnTpr, nfiles2, fnPull);
        }

        /* Read file that selects the pull group to be used */
        if (opt.fnGroupsel != NULL)
        {
            readPullGroupSelection(&opt, fninTpr, nfiles);
        }

        window = initUmbrellaWindows(nfiles);
        read_tpr_pullxf_files(fninTpr, fninPull, nfiles, &header, window, &opt);
    }
    else
    {   /* reading pdo files */
        if  (opt.fnGroupsel != NULL)
        {
            gmx_fatal(FARGS, "Reading a -is file is not supported with PDO input files.\n"
                      "Use awk or a similar tool to pick the required pull groups from your PDO files\n");
        }
        read_wham_in(opt.fnPdo, &fninPdo, &nfiles, &opt);
        printf("Found %d pdo files in %s\n", nfiles, opt.fnPdo);
        window = initUmbrellaWindows(nfiles);
        read_pdo_files(fninPdo, nfiles, &header, window, &opt);
    }
    nwins = nfiles;

    /* enforce equal weight for all histograms? */
    if (opt.bHistEq)
    {
        enforceEqualWeights(window, nwins);
    }

    /* write histograms */
    histout = xvgropen(opt2fn("-hist", NFILE, fnm), "Umbrella histograms",
                       xlabel, "count", opt.oenv);
    for (l = 0; l < opt.bins; ++l)
    {
        fprintf(histout, "%e\t", (double)(l+0.5)/opt.bins*(opt.max-opt.min)+opt.min);
        for (i = 0; i < nwins; ++i)
        {
            for (j = 0; j < window[i].nPull; ++j)
            {
                fprintf(histout, "%e\t", window[i].Histo[j][l]);
            }
        }
        fprintf(histout, "\n");
    }
    xvgrclose(histout);
    printf("Wrote %s\n", opt2fn("-hist", NFILE, fnm));
    if (opt.bHistOnly)
    {
        printf("Wrote histograms to %s, now exiting.\n", opt2fn("-hist", NFILE, fnm));
        return 0;
    }

    /* Using tabulated umbrella potential */
    if (opt.bTab)
    {
        setup_tab(opt2fn("-tab", NFILE, fnm), &opt);
    }

    /* Integrated autocorrelation times provided ? */
    if (opt.bTauIntGiven)
    {
        readIntegratedAutocorrelationTimes(window, nwins, opt2fn("-iiact", NFILE, fnm));
    }

    /* Compute integrated autocorrelation times */
    if (opt.bCalcTauInt)
    {
        calcIntegratedAutocorrelationTimes(window, nwins, &opt, opt2fn("-oiact", NFILE, fnm));
    }

    /* calc average and sigma for each histogram
       (maybe required for bootstrapping. If not, this is fast anyhow) */
    if (opt.nBootStrap && opt.bsMethod == bsMethod_trajGauss)
    {
        averageSigma(window, nwins);
    }

    /* Get initial potential by simple integration */
    if (opt.bInitPotByIntegration)
    {
        guessPotByIntegration(window, nwins, &opt);
    }

    /* Check if complete reaction coordinate is covered */
    checkReactionCoordinateCovered(window, nwins, &opt);

    /* Calculate profile */
    snew(profile, opt.bins);
    if (opt.verbose)
    {
        opt.stepchange = 1;
    }
    i = 0;
    do
    {
        if ( (i%opt.stepUpdateContrib) == 0)
        {
            setup_acc_wham(profile, window, nwins, &opt);
        }
        if (maxchange < opt.Tolerance)
        {
            bExact = TRUE;
            /* if (opt.verbose) */
            printf("Switched to exact iteration in iteration %d\n", i);
        }
        calc_profile(profile, window, nwins, &opt, bExact);
        if (((i%opt.stepchange) == 0 || i == 1) && i != 0)
        {
            printf("\t%4d) Maximum change %e\n", i, maxchange);
        }
        i++;
    }
    while ( (maxchange = calc_z(profile, window, nwins, &opt, bExact)) > opt.Tolerance || !bExact);
    printf("Converged in %d iterations. Final maximum change %g\n", i, maxchange);

    /* calc error from Kumar's formula */
    /* Unclear how the error propagates along reaction coordinate, therefore
       commented out  */
    /* calc_error_kumar(profile,window, nwins,&opt); */

    /* Write profile in energy units? */
    if (opt.bLog)
    {
        prof_normalization_and_unit(profile, &opt);
        strcpy(ylabel, en_unit_label[opt.unit]);
        strcpy(title, "Umbrella potential");
    }
    else
    {
        strcpy(ylabel, "Density of states");
        strcpy(title, "Density of states");
    }

    /* symmetrize profile around z=0? */
    if (opt.bSym)
    {
        symmetrizeProfile(profile, &opt);
    }

    /* write profile or density of states */
    profout = xvgropen(opt2fn("-o", NFILE, fnm), title, xlabel, ylabel, opt.oenv);
    for (i = 0; i < opt.bins; ++i)
    {
        fprintf(profout, "%e\t%e\n", (double)(i+0.5)/opt.bins*(opt.max-opt.min)+opt.min, profile[i]);
    }
    xvgrclose(profout);
    printf("Wrote %s\n", opt2fn("-o", NFILE, fnm));

    /* Bootstrap Method */
    if (opt.nBootStrap)
    {
        do_bootstrapping(opt2fn("-bsres", NFILE, fnm), opt2fn("-bsprof", NFILE, fnm),
                         opt2fn("-hist", NFILE, fnm),
                         ylabel, profile, window, nwins, &opt);
    }

    sfree(profile);
    freeUmbrellaWindows(window, nfiles);

    printf("\nIn case you use results from g_wham for a publication, please cite:\n");
    please_cite(stdout, "Hub2010");

    return 0;
}
