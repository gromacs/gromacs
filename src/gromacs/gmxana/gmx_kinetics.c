/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "gromacs/commandline/pargs.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "gromacs/fileio/futil.h"
#include "readinp.h"
#include "txtdump.h"
#include "gstat.h"
#include "xvgr.h"
#include "physics.h"
#include "gmx_ana.h"

enum {
    epAuf, epEuf, epAfu, epEfu, epNR
};
enum {
    eqAif, eqEif, eqAfi, eqEfi, eqAui, eqEui, eqAiu, eqEiu, eqNR
};
static char *eep[epNR] = { "Af", "Ef", "Au", "Eu" };
static char *eeq[eqNR] = { "Aif", "Eif", "Afi", "Efi", "Aui", "Eui", "Aiu", "Eiu" };

typedef struct {
    int       nreplica;  /* Number of replicas in the calculation                   */
    int       nframe;    /* Number of time frames                                   */
    int       nstate;    /* Number of states the system can be in, e.g. F,I,U       */
    int       nparams;   /* Is 2, 4 or 8                                            */
    gmx_bool *bMask;     /* Determine whether this replica is part of the d2 comp.  */
    gmx_bool  bSum;
    gmx_bool  bDiscrete; /* Use either discrete folding (0/1) or a continuous       */
    /* criterion */
    int       nmask;     /* Number of replicas taken into account                   */
    real      dt;        /* Timestep between frames                                 */
    int       j0, j1;    /* Range of frames used in calculating delta               */
    real    **temp, **data, **data2;
    int     **state;     /* State index running from 0 (F) to nstate-1 (U)          */
    real    **beta, **fcalt, **icalt;
    real     *time, *sumft, *sumit, *sumfct, *sumict;
    real     *params;
    real     *d2_replica;
} t_remd_data;


int gmx_kinetics(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] reads two [TT].xvg[tt] files, each one containing data for N replicas.",
        "The first file contains the temperature of each replica at each timestep,",
        "and the second contains real values that can be interpreted as",
        "an indicator for folding. If the value in the file is larger than",
        "the cutoff it is taken to be unfolded and the other way around.[PAR]",
        "From these data an estimate of the forward and backward rate constants",
        "for folding is made at a reference temperature. In addition,",
        "a theoretical melting curve and free energy as a function of temperature",
        "are printed in an [TT].xvg[tt] file.[PAR]",
        "The user can give a max value to be regarded as intermediate",
        "([TT]-ucut[tt]), which, when given will trigger the use of an intermediate state",
        "in the algorithm to be defined as those structures that have",
        "cutoff < DATA < ucut. Structures with DATA values larger than ucut will",
        "not be regarded as potential folders. In this case 8 parameters are optimized.[PAR]",
        "The average fraction foled is printed in an [TT].xvg[tt] file together with the fit to it.",
        "If an intermediate is used a further file will show the build of the intermediate and the fit to that process.[PAR]",
        "The program can also be used with continuous variables (by setting",
        "[TT]-nodiscrete[tt]). In this case kinetics of other processes can be",
        "studied. This is very much a work in progress and hence the manual",
        "(this information) is lagging behind somewhat.[PAR]",
        "In order to run [THISMODULE], GROMACS must be compiled with the GNU",
        "scientific library."
    };
    static int      nreplica  = 1;
    static real     tref      = 298.15;
    static real     cutoff    = 0.2;
    static real     ucut      = 0.0;
    static real     Euf       = 10;
    static real     Efu       = 30;
    static real     Ei        = 10;
    static gmx_bool bHaveT    = TRUE;
    static real     t0        = -1;
    static real     t1        = -1;
    static real     tb        = 0;
    static real     te        = 0;
    static real     tol       = 1e-3;
    static int      maxiter   = 100;
    static int      skip      = 0;
    static int      nmult     = 1;
    static gmx_bool bBack     = TRUE;
    static gmx_bool bSplit    = TRUE;
    static gmx_bool bSum      = TRUE;
    static gmx_bool bDiscrete = TRUE;
    t_pargs         pa[]      = {
        { "-time",    FALSE, etBOOL, {&bHaveT},
          "Expect a time in the input" },
        { "-b",       FALSE, etREAL, {&tb},
          "First time to read from set" },
        { "-e",       FALSE, etREAL, {&te},
          "Last time to read from set" },
        { "-bfit",    FALSE, etREAL, {&t0},
          "Time to start the fit from" },
        { "-efit",    FALSE, etREAL, {&t1},
          "Time to end the fit" },
        { "-T",       FALSE, etREAL, {&tref},
          "Reference temperature for computing rate constants" },
        { "-n",       FALSE, etINT, {&nreplica},
          "Read data for this number of replicas. Only necessary when files are written in xmgrace format using @type and & as delimiters." },
        { "-cut",     FALSE, etREAL, {&cutoff},
          "Cut-off (max) value for regarding a structure as folded" },
        { "-ucut",    FALSE, etREAL, {&ucut},
          "Cut-off (max) value for regarding a structure as intermediate (if not folded)" },
        { "-euf",     FALSE, etREAL, {&Euf},
          "Initial guess for energy of activation for folding (kJ/mol)" },
        { "-efu",     FALSE, etREAL, {&Efu},
          "Initial guess for energy of activation for unfolding (kJ/mol)" },
        { "-ei",      FALSE, etREAL, {&Ei},
          "Initial guess for energy of activation for intermediates (kJ/mol)" },
        { "-maxiter", FALSE, etINT, {&maxiter},
          "Max number of iterations" },
        { "-back",    FALSE, etBOOL, {&bBack},
          "Take the back reaction into account" },
        { "-tol",     FALSE, etREAL, {&tol},
          "Absolute tolerance for convergence of the Nelder and Mead simplex algorithm" },
        { "-skip",    FALSE, etINT, {&skip},
          "Skip points in the output [TT].xvg[tt] file" },
        { "-split",   FALSE, etBOOL, {&bSplit},
          "Estimate error by splitting the number of replicas in two and refitting" },
        { "-sum",     FALSE, etBOOL, {&bSum},
          "Average folding before computing [GRK]chi[grk]^2" },
        { "-discrete", FALSE, etBOOL, {&bDiscrete},
          "Use a discrete folding criterion (F <-> U) or a continuous one" },
        { "-mult",    FALSE, etINT, {&nmult},
          "Factor to multiply the data with before discretization" }
    };
#define NPA asize(pa)

    FILE        *fp;
    real         dt_t, dt_d, dt_d2;
    int          nset_t, nset_d, nset_d2, n_t, n_d, n_d2, i;
    const char  *tfile, *dfile, *dfile2;
    t_remd_data  remd;
    output_env_t oenv;

    t_filenm     fnm[] = {
        { efXVG, "-f",    "temp",    ffREAD   },
        { efXVG, "-d",    "data",    ffREAD   },
        { efXVG, "-d2",   "data2",   ffOPTRD  },
        { efXVG, "-o",    "ft_all",  ffWRITE  },
        { efXVG, "-o2",   "it_all",  ffOPTWR  },
        { efXVG, "-o3",   "ft_repl", ffOPTWR  },
        { efXVG, "-ee",   "err_est", ffOPTWR  },
        { efLOG, "-g",    "remd",    ffWRITE  },
        { efXVG, "-m",    "melt",    ffWRITE  }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_BE_NICE | PCA_TIME_UNIT,
                           NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    fprintf(stderr, "You have requested code to run that is deprecated.\n");
    fprintf(stderr, "Revert to an older GROMACS version or help in porting the code.\n");

    return 0;
}
