/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/* The make_edi program was generously contributed by Oliver Lange, based
 * on the code from g_anaeig. You can reach him as olange@gwdg.de. He
 * probably also holds copyright to the following code.
 */
#include "gmxpre.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/eigio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct
{
    real     deltaF0;
    gmx_bool bHarmonic;
    gmx_bool bConstForce; /* Do constant force flooding instead of
                             evaluating a flooding potential             */
    real tau;
    real deltaF;
    real kT;
    real constEfl;
    real alpha2;
} t_edflood;


/* This type is for the average, reference, target, and origin structure   */
struct edix
{
    int   nr;   /* number of atoms this structure contains */
    int*  anrs; /* atom index numbers                      */
    rvec* x;    /* positions                               */
};


typedef struct edipar
{
    int      nini;          /* total Nr of atoms                    */
    gmx_bool fitmas;        /* true if trans fit with cm            */
    gmx_bool pcamas;        /* true if mass-weighted PCA            */
    int      presteps;      /* number of steps to run without any
                             *    perturbations ... just monitoring */
    int         outfrq;     /* freq (in steps) of writing to edo    */
    int         maxedsteps; /* max nr of steps per cycle            */
    struct edix sref;       /* reference positions, to these fitting
                             * will be done                         */
    struct edix sav;        /* average positions                    */
    struct edix star;       /* target positions                     */
    struct edix sori;       /* origin positions                     */
    real        slope;      /* minimal slope in acceptance radexp   */
    int         ned;        /* Nr of atoms in essdyn buffer         */
    t_edflood   flood;      /* parameters especially for flooding   */
} t_edipar;


static void make_t_edx(struct edix* edx, int natoms, rvec* pos, int index[])
{
    edx->nr   = natoms;
    edx->anrs = index;
    edx->x    = pos;
}

static void write_t_edx(FILE* fp, struct edix edx, const char* comment)
{
    /*here we copy only the pointers into the t_edx struct
       no data is copied and edx.box is ignored  */
    int i;
    fprintf(fp, "#%s \n %d \n", comment, edx.nr);
    for (i = 0; i < edx.nr; i++)
    {
        fprintf(fp, "%d  %f  %f  %f\n", (edx.anrs)[i] + 1, (edx.x)[i][XX], (edx.x)[i][YY], (edx.x)[i][ZZ]);
    }
}

static int sscan_list(int* list[], const char* str, const char* listname)
{
    /*this routine scans a string of the form 1,3-6,9 and returns the
       selected numbers (in this case 1 3 4 5 6 9) in NULL-terminated array of integers.
       memory for this list will be allocated  in this routine -- sscan_list expects *list to
       be a NULL-Pointer

       listname is a string used in the errormessage*/


    int   i, istep;
    char  c;
    char *pos, *startpos, *step;
    int   n = std::strlen(str);

    /*enums to define the different lexical stati */
    enum
    {
        sBefore,
        sNumber,
        sMinus,
        sRange,
        sZero,
        sSmaller,
        sError,
        sSteppedRange
    };

    int status     = sBefore; /*status of the deterministic automat to scan str   */
    int number     = 0;
    int end_number = 0;

    char* start = nullptr; /*holds the string of the number behind a ','*/
    char* end   = nullptr; /*holds the string of the number behind a '-' */

    int nvecs = 0; /* counts the number of vectors in the list*/

    step = nullptr;
    snew(pos, n + 4);
    startpos = pos;
    std::strcpy(pos, str);
    pos[n]     = ',';
    pos[n + 1] = '1';
    pos[n + 2] = '\0';

    *list = nullptr;

    while ((c = *pos) != 0)
    {
        switch (status)
        {
            /* expect a number */
            case sBefore:
                if (std::isdigit(c))
                {
                    start  = pos;
                    status = sNumber;
                    break;
                }
                else
                {
                    status = sError;
                }
                break;

            /* have read a number, expect ',' or '-' */
            case sNumber:
                if (c == ',')
                {
                    /*store number*/
                    srenew(*list, nvecs + 1);
                    (*list)[nvecs++] = number = std::strtol(start, nullptr, 10);
                    status                    = sBefore;
                    if (number == 0)
                    {
                        status = sZero;
                    }
                    break;
                }
                else if (c == '-')
                {
                    status = sMinus;
                    break;
                }
                else if (std::isdigit(c))
                {
                    break;
                }
                else
                {
                    status = sError;
                }
                break;

            /* have read a '-' -> expect a number */
            case sMinus:
                if (std::isdigit(c))
                {
                    end    = pos;
                    status = sRange;
                    break;
                }
                else
                {
                    status = sError;
                }
                break;

            case sSteppedRange:
                if (std::isdigit(c))
                {
                    if (step)
                    {
                        status = sError;
                        break;
                    }
                    else
                    {
                        step = pos;
                    }
                    status = sRange;
                    break;
                }
                else
                {
                    status = sError;
                }
                break;

            /* have read the number after a minus, expect ',' or ':' */
            case sRange:
                if (c == ',')
                {
                    /*store numbers*/
                    end_number = std::strtol(end, nullptr, 10);
                    number     = std::strtol(start, nullptr, 10);
                    status     = sBefore;
                    if (number == 0)
                    {
                        status = sZero;
                        break;
                    }
                    if (end_number <= number)
                    {
                        status = sSmaller;
                        break;
                    }
                    srenew(*list, nvecs + end_number - number + 1);
                    if (step)
                    {
                        istep = strtol(step, nullptr, 10);
                        step  = nullptr;
                    }
                    else
                    {
                        istep = 1;
                    }
                    for (i = number; i <= end_number; i += istep)
                    {
                        (*list)[nvecs++] = i;
                    }
                    break;
                }
                else if (c == ':')
                {
                    status = sSteppedRange;
                    break;
                }
                else if (std::isdigit(c))
                {
                    break;
                }
                else
                {
                    status = sError;
                }
                break;

            /* format error occurred */
            case sError:
                gmx_fatal(FARGS,
                          "Error in the list of eigenvectors for %s at pos %td with char %c",
                          listname,
                          pos - startpos,
                          *(pos - 1));
            /* logical error occurred */
            case sZero:
                gmx_fatal(FARGS,
                          "Error in the list of eigenvectors for %s at pos %td: eigenvector 0 is "
                          "not valid",
                          listname,
                          pos - startpos);
            case sSmaller:
                gmx_fatal(FARGS,
                          "Error in the list of eigenvectors for %s at pos %td: second index %d is "
                          "not bigger than %d",
                          listname,
                          pos - startpos,
                          end_number,
                          number);
        }
        ++pos; /* read next character */
    }          /*scanner has finished */

    /* append zero to list of eigenvectors */
    srenew(*list, nvecs + 1);
    (*list)[nvecs] = 0;
    sfree(startpos);
    return nvecs;
} /*sscan_list*/

static void
write_eigvec(FILE* fp, int natoms, int eig_list[], rvec** eigvecs, int nvec, const char* grouptitle, real steps[])
{
    /* eig_list is a zero-terminated list of indices into the eigvecs array.
       eigvecs are coordinates of eigenvectors
       grouptitle to write in the comment line
       steps  -- array with stepsizes for evLINFIX, evLINACC and evRADACC
     */

    int  n = 0, i;
    rvec x;
    while (eig_list[n++])
    {
        /*count selected eigenvecs*/
    }
    fprintf(fp, "# NUMBER OF EIGENVECTORS + %s\n %d\n", grouptitle, n - 1);

    /* write list of eigenvector indicess */
    for (n = 0; eig_list[n]; n++)
    {
        if (steps)
        {
            fprintf(fp, "%8d   %g\n", eig_list[n], steps[n]);
        }
        else
        {
            fprintf(fp, "%8d   %g\n", eig_list[n], 1.0);
        }
    }
    n = 0;

    /* dump coordinates of the selected eigenvectors */
    while (eig_list[n])
    {
        for (i = 0; i < natoms; i++)
        {
            if (eig_list[n] > nvec)
            {
                gmx_fatal(FARGS,
                          "Selected eigenvector %d is higher than maximum number %d of available "
                          "eigenvectors",
                          eig_list[n],
                          nvec);
            }
            copy_rvec(eigvecs[eig_list[n] - 1][i], x);
            fprintf(fp, "%8.5f %8.5f %8.5f\n", x[XX], x[YY], x[ZZ]);
        }
        n++;
    }
}


/*enum referring to the different lists of eigenvectors*/
enum
{
    evLINFIX,
    evLINACC,
    evFLOOD,
    evRADFIX,
    evRADACC,
    evRADCON,
    evMON,
    evNr
};
#define MAGIC 670


static void write_the_whole_thing(const char* filename,
                                  t_edipar*   edpars,
                                  rvec**      eigvecs,
                                  int         nvec,
                                  int*        eig_listen[],
                                  real*       evStepList[])
{
    FILE* fp = gmx_ffopen(filename, "w");

    /* write edi-file */

    /*Header*/
    fprintf(fp,
            "#MAGIC\n %d \n#NINI\n %d\n#FITMAS\n %d\n#ANALYSIS_MAS\n %d\n",
            MAGIC,
            edpars->nini,
            int(edpars->fitmas),
            int(edpars->pcamas));
    fprintf(fp, "#OUTFRQ\n %d\n#MAXLEN\n %d\n#SLOPECRIT\n %f\n", edpars->outfrq, edpars->maxedsteps, edpars->slope);
    fprintf(fp,
            "#PRESTEPS\n %d\n#DELTA_F0\n %f\n#INIT_DELTA_F\n %f\n#TAU\n %f\n#EFL_NULL\n "
            "%f\n#ALPHA2\n %f\n#KT\n %f\n#HARMONIC\n %d\n#CONST_FORCE_FLOODING\n %d\n",
            edpars->presteps,
            edpars->flood.deltaF0,
            edpars->flood.deltaF,
            edpars->flood.tau,
            edpars->flood.constEfl,
            edpars->flood.alpha2,
            edpars->flood.kT,
            int(edpars->flood.bHarmonic),
            int(edpars->flood.bConstForce));

    /* Average and reference positions */
    write_t_edx(fp, edpars->sref, "NREF, XREF");
    write_t_edx(fp, edpars->sav, "NAV, XAV");

    /*Eigenvectors */

    write_eigvec(fp, edpars->ned, eig_listen[evMON], eigvecs, nvec, "COMPONENTS GROUP 1", nullptr);
    write_eigvec(fp, edpars->ned, eig_listen[evLINFIX], eigvecs, nvec, "COMPONENTS GROUP 2", evStepList[evLINFIX]);
    write_eigvec(fp, edpars->ned, eig_listen[evLINACC], eigvecs, nvec, "COMPONENTS GROUP 3", evStepList[evLINACC]);
    write_eigvec(fp, edpars->ned, eig_listen[evRADFIX], eigvecs, nvec, "COMPONENTS GROUP 4", evStepList[evRADFIX]);
    write_eigvec(fp, edpars->ned, eig_listen[evRADACC], eigvecs, nvec, "COMPONENTS GROUP 5", nullptr);
    write_eigvec(fp, edpars->ned, eig_listen[evRADCON], eigvecs, nvec, "COMPONENTS GROUP 6", nullptr);
    write_eigvec(fp, edpars->ned, eig_listen[evFLOOD], eigvecs, nvec, "COMPONENTS GROUP 7", evStepList[evFLOOD]);


    /*Target and Origin positions */
    write_t_edx(fp, edpars->star, "NTARGET, XTARGET");
    write_t_edx(fp, edpars->sori, "NORIGIN, XORIGIN");

    gmx_ffclose(fp);
}

static int read_conffile(const char* confin, rvec** x)
{
    t_topology top;
    matrix     box;
    printf("read coordnumber from file %s\n", confin);
    read_tps_conf(confin, &top, nullptr, x, nullptr, box, FALSE);
    printf("number of coordinates in file %d\n", top.atoms.nr);
    return top.atoms.nr;
}


static void read_eigenvalues(const int vecs[], const char* eigfile, real values[], gmx_bool bHesse, real kT, int natoms_average_struct)
{
    int      neig, nrow, i;
    double** eigval;

    neig = read_xvg(eigfile, &eigval, &nrow);

    fprintf(stderr, "Read %d eigenvalues\n", neig);
    for (i = bHesse ? 6 : 0; i < neig; i++)
    {
        if (eigval[1][i] < -0.001 && bHesse)
        {
            fprintf(stderr,
                    "WARNING: The Hessian Matrix has negative eigenvalue %f, we set it to zero (no "
                    "flooding in this direction)\n\n",
                    eigval[1][i]);
        }

        if (eigval[1][i] < 0)
        {
            eigval[1][i] = 0;
        }
    }
    if (bHesse)
    {
        for (i = 0; vecs[i]; i++)
        {
            if (vecs[i] < 7)
            {
                gmx_fatal(FARGS,
                          "ERROR: You have chosen one of the first 6 eigenvectors of the HESSE "
                          "Matrix. That does not make sense, since they correspond to the 6 "
                          "rotational and translational degrees of freedom.\n\n");
            }
            values[i] = eigval[1][vecs[i] - 1] / kT;
        }
    }
    else
    {
        for (i = 0; vecs[i]; i++)
        {
            /* Make sure this eigenvalue does not correspond to one of the last 6 eigenvectors of the
             * covariance matrix. These correspond to the rotational and translational degrees of
             * freedom and will be zero within numerical accuracy.
             *
             * Note that the total number of eigenvectors produced by gmx covar depends on:
             * 1) the total number of degrees of freedom of the system (3N, with N the number of atoms)
             * 2) the number S of independent configurations fed into gmx covar.
             * For long trajectories with lots of frames, usually S >= 3N + 1, so that one indeed gets
             * 3N eigenvalues (of which the last 6 will have zero eigenvalues).
             * For S < 3N + 1, however, the covariance matrix becomes rank deficient, and the number
             * of possible eigenvalues is just S - 1. Since in make_edi we only know N but not S, we can
             * only warn the user if he picked one of the last 6 of 3N.
             */
            if (vecs[i] > (3 * natoms_average_struct - 6))
            {
                gmx_fatal(FARGS,
                          "ERROR: You have chosen one of the last 6 eigenvectors of the COVARIANCE "
                          "Matrix. That does not make sense, since they correspond to the 6 "
                          "rotational and translational degrees of freedom.\n\n");
            }
            values[i] = 1 / eigval[1][vecs[i] - 1];
        }
    }
    /* free memory */
    for (i = 0; i < nrow; i++)
    {
        sfree(eigval[i]);
    }
    sfree(eigval);
}


static real* scan_vecparams(const char* str, const char* par, int nvecs)
{
    char   f0[256], f1[256]; /*format strings adapted every pass of the loop*/
    double d;
    int    i;
    real*  vec_params;

    snew(vec_params, nvecs);
    if (str)
    {
        f0[0] = '\0';
        for (i = 0; (i < nvecs); i++)
        {
            std::strcpy(f1, f0);    /*f0 is the format string for the "to-be-ignored" numbers*/
            std::strcat(f1, "%lf"); /*and f1 to read the actual number in this pass of the loop*/
            if (sscanf(str, f1, &d) != 1)
            {
                gmx_fatal(FARGS, "Not enough elements for %s parameter (I need %d)", par, nvecs);
            }
            vec_params[i] = d;
            std::strcat(f0, "%*s");
        }
    }
    return vec_params;
}


static void init_edx(struct edix* edx)
{
    edx->nr = 0;
    snew(edx->x, 1);
    snew(edx->anrs, 1);
}

static void
filter2edx(struct edix* edx, int nindex, int index[], int ngro, const int igro[], const rvec* x, const char* structure)
{
    /* filter2edx copies coordinates from x to edx which are given in index
     */

    int pos, i;
    int ix = edx->nr;
    edx->nr += nindex;
    srenew(edx->x, edx->nr);
    srenew(edx->anrs, edx->nr);
    for (i = 0; i < nindex; i++, ix++)
    {
        for (pos = 0; pos < ngro - 1 && igro[pos] != index[i]; ++pos) {} /*search element in igro*/
        if (igro[pos] != index[i])
        {
            gmx_fatal(FARGS, "Couldn't find atom with index %d in structure %s", index[i], structure);
        }
        edx->anrs[ix] = index[i];
        copy_rvec(x[pos], edx->x[ix]);
    }
}

static void get_structure(const t_atoms* atoms,
                          const char*    IndexFile,
                          const char*    StructureFile,
                          struct edix*   edx,
                          int            nfit,
                          int            ifit[],
                          int            nav,
                          int            index[])
{
    int*  igro; /*index corresponding to target or origin structure*/
    int   ngro;
    int   ntar;
    rvec* xtar;
    char* grpname;


    ntar = read_conffile(StructureFile, &xtar);
    printf("Select an index group of %d elements that corresponds to the atoms in the structure "
           "file %s\n",
           ntar,
           StructureFile);
    get_index(atoms, IndexFile, 1, &ngro, &igro, &grpname);
    if (ngro != ntar)
    {
        gmx_fatal(FARGS, "You selected an index group with %d elements instead of %d", ngro, ntar);
    }
    init_edx(edx);
    filter2edx(edx, nfit, ifit, ngro, igro, xtar, StructureFile);

    /* If average and reference/fitting structure differ, append the average structure as well */
    if (ifit != index) /*if fit structure is different append these coordinates, too -- don't mind duplicates*/
    {
        filter2edx(edx, nav, index, ngro, igro, xtar, StructureFile);
    }
}

int gmx_make_edi(int argc, char* argv[])
{

    static const char* desc[] = {
        "[THISMODULE] generates an essential dynamics (ED) sampling input file to be used with ",
        "[TT]mdrun[tt] based on eigenvectors of a covariance matrix ([gmx-covar]) or from a",
        "normal modes analysis ([gmx-nmeig]).",
        "ED sampling can be used to manipulate the position along collective coordinates",
        "(eigenvectors) of (biological) macromolecules during a simulation. Particularly,",
        "it may be used to enhance the sampling efficiency of MD simulations by stimulating",
        "the system to explore new regions along these collective coordinates. A number",
        "of different algorithms are implemented to drive the system along the eigenvectors",
        "([TT]-linfix[tt], [TT]-linacc[tt], [TT]-radfix[tt], [TT]-radacc[tt], [TT]-radcon[tt]),",
        "to keep the position along a certain (set of) coordinate(s) fixed ([TT]-linfix[tt]),",
        "or to only monitor the projections of the positions onto",
        "these coordinates ([TT]-mon[tt]).[PAR]",
        "References:[PAR]",
        "A. Amadei, A.B.M. Linssen, B.L. de Groot, D.M.F. van Aalten and ",
        "H.J.C. Berendsen; An efficient method for sampling the essential subspace ",
        "of proteins., J. Biomol. Struct. Dyn. 13:615-626 (1996)[PAR]",
        "B.L. de Groot, A. Amadei, D.M.F. van Aalten and H.J.C. Berendsen; ",
        "Towards an exhaustive sampling of the configurational spaces of the ",
        "two forms of the peptide hormone guanylin,",
        "J. Biomol. Struct. Dyn. 13 : 741-751 (1996)[PAR]",
        "B.L. de Groot, A.Amadei, R.M. Scheek, N.A.J. van Nuland and H.J.C. Berendsen; ",
        "An extended sampling of the configurational space of HPr from E. coli",
        "Proteins: Struct. Funct. Gen. 26: 314-322 (1996)",
        "[PAR]You will be prompted for one or more index groups that correspond to the ",
        "eigenvectors,",
        "reference structure, target positions, etc.[PAR]",

        "[TT]-mon[tt]: monitor projections of the coordinates onto selected eigenvectors.[PAR]",
        "[TT]-linfix[tt]: perform fixed-step linear expansion along selected eigenvectors.[PAR]",
        "[TT]-linacc[tt]: perform acceptance linear expansion along selected eigenvectors.",
        "(steps in the desired directions will be accepted, others will be rejected).[PAR]",
        "[TT]-radfix[tt]: perform fixed-step radius expansion along selected eigenvectors.[PAR]",
        "[TT]-radacc[tt]: perform acceptance radius expansion along selected eigenvectors.",
        "(steps in the desired direction will be accepted, others will be rejected).",
        "[BB]Note:[bb] by default the starting MD structure will be taken as origin of the first",
        "expansion cycle for radius expansion. If [TT]-ori[tt] is specified, you will be able",
        "to read in a structure file that defines an external origin.[PAR]",
        "[TT]-radcon[tt]: perform acceptance radius contraction along selected eigenvectors",
        "towards a target structure specified with [TT]-tar[tt].[PAR]",
        "NOTE: each eigenvector can be selected only once. [PAR]",
        "[TT]-outfrq[tt]: frequency (in steps) of writing out projections etc. to [REF].xvg[ref] ",
        "file",
        "[PAR]",
        "[TT]-slope[tt]: minimal slope in acceptance radius expansion. A new expansion",
        "cycle will be started if the spontaneous increase of the radius (in nm/step)",
        "is less than the value specified.[PAR]",
        "[TT]-maxedsteps[tt]: maximum number of steps per cycle in radius expansion",
        "before a new cycle is started.[PAR]",
        "Note on the parallel implementation: since ED sampling is a 'global' thing",
        "(collective coordinates etc.), at least on the 'protein' side, ED sampling",
        "is not very parallel-friendly from an implementation point of view. Because",
        "parallel ED requires some extra communication, expect the performance to be",
        "lower as in a free MD simulation, especially on a large number of ranks and/or",
        "when the ED group contains a lot of atoms. [PAR]",
        "Please also note that if your ED group contains more than a single protein,",
        "then the [REF].tpr[ref] file must contain the correct PBC representation of the ED group.",
        "Take a look on the initial RMSD from the reference structure, which is printed",
        "out at the start of the simulation; if this is much higher than expected, one",
        "of the ED molecules might be shifted by a box vector. [PAR]",
        "All ED-related output of [TT]mdrun[tt] (specify with [TT]-eo[tt]) is written to a ",
        "[REF].xvg[ref] file as a function of time in intervals of OUTFRQ steps.[PAR]",
        "[BB]Note[bb] that you can impose multiple ED constraints and flooding potentials in",
        "a single simulation (on different molecules) if several [REF].edi[ref] files were ",
        "concatenated first. The constraints are applied in the order they appear in ",
        "the [REF].edi[ref] file. Depending on what was specified in the [REF].edi[ref] ",
        "input file, the output file contains for each ED dataset",
        "",
        " * the RMSD of the fitted molecule to the reference structure (for atoms involved in ",
        "   fitting prior to calculating the ED constraints)",
        " * projections of the positions onto selected eigenvectors",
        "",
        "FLOODING:[PAR]",
        "with [TT]-flood[tt], you can specify which eigenvectors are used to compute a flooding ",
        "potential,",
        "which will lead to extra forces expelling the structure out of the region described",
        "by the covariance matrix. If you switch -restrain the potential is inverted and the ",
        "structure is kept in that region.",
        "[PAR]",
        "The origin is normally the average structure stored in the [TT]eigvec.trr[tt] file.",
        "It can be changed with [TT]-ori[tt] to an arbitrary position in configuration space.",
        "With [TT]-tau[tt], [TT]-deltaF0[tt], and [TT]-Eflnull[tt] you control the flooding ",
        "behaviour. Efl is the flooding strength, it is updated according to the rule of ",
        "adaptive flooding. Tau is the time constant of adaptive flooding, high ",
        "[GRK]tau[grk] means slow adaption (i.e. growth). ",
        "DeltaF0 is the flooding strength you want to reach after tau ps of simulation.",
        "To use constant Efl set [TT]-tau[tt] to zero.",
        "[PAR]",
        "[TT]-alpha[tt] is a fudge parameter to control the width of the flooding potential. A ",
        "value of 2 has been found",
        "to give good results for most standard cases in flooding of proteins.",
        "[GRK]alpha[grk] basically accounts for incomplete sampling, if you sampled further the ",
        "width of the ensemble would",
        "increase, this is mimicked by [GRK]alpha[grk] > 1.",
        "For restraining, [GRK]alpha[grk] < 1 can give you smaller width in the restraining ",
        "potential.",
        "[PAR]",
        "RESTART and FLOODING:",
        "If you want to restart a crashed flooding simulation please find the values deltaF and ",
        "Efl in",
        "the output file and manually put them into the [REF].edi[ref] file under DELTA_F0 and ",
        "EFL_NULL."
    };

    /* Save all the params in this struct and then save it in an edi file.
     * ignoring fields nmass,massnrs,mass,tmass,nfit,fitnrs,edo
     */
    static t_edipar edi_params;

    enum
    {
        evStepNr = evRADFIX + 1
    };
    static const char* evSelections[evNr] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
    static const char* evOptions[evNr]    = { "-linfix", "-linacc", "-flood", "-radfix",
                                           "-radacc", "-radcon", "-mon" };
    static const char* evParams[evStepNr]      = { nullptr, nullptr };
    static const char* evStepOptions[evStepNr] = { "-linstep", "-accdir", "-not_used", "-radstep" };
    static const char* ConstForceStr;
    static real*       evStepList[evStepNr];
    static real        radstep  = 0.0;
    static real        deltaF0  = 150;
    static real        deltaF   = 0;
    static real        tau      = .1;
    static real        constEfl = 0.0;
    static real        alpha    = 1;
    static int         eqSteps  = 0;
    static int*        listen[evNr];
    static real        T         = 300.0;
    const real         kB        = 2.5 / 300.0; /* k_boltzmann in MD units */
    static gmx_bool    bRestrain = FALSE;
    static gmx_bool    bHesse    = FALSE;
    static gmx_bool    bHarmonic = FALSE;
    t_pargs            pa[]      = {
        { "-mon",
          FALSE,
          etSTR,
          { &evSelections[evMON] },
          "Indices of eigenvectors for projections of x (e.g. 1,2-5,9) or 1-100:10 means 1 11 21 "
          "31 ... 91" },
        { "-linfix",
          FALSE,
          etSTR,
          { &evSelections[0] },
          "Indices of eigenvectors for fixed increment linear sampling" },
        { "-linacc",
          FALSE,
          etSTR,
          { &evSelections[1] },
          "Indices of eigenvectors for acceptance linear sampling" },
        { "-radfix",
          FALSE,
          etSTR,
          { &evSelections[3] },
          "Indices of eigenvectors for fixed increment radius expansion" },
        { "-radacc",
          FALSE,
          etSTR,
          { &evSelections[4] },
          "Indices of eigenvectors for acceptance radius expansion" },
        { "-radcon",
          FALSE,
          etSTR,
          { &evSelections[5] },
          "Indices of eigenvectors for acceptance radius contraction" },
        { "-flood", FALSE, etSTR, { &evSelections[2] }, "Indices of eigenvectors for flooding" },
        { "-outfrq",
          FALSE,
          etINT,
          { &edi_params.outfrq },
          "Frequency (in steps) of writing output in [REF].xvg[ref] file" },
        { "-slope",
          FALSE,
          etREAL,
          { &edi_params.slope },
          "Minimal slope in acceptance radius expansion" },
        { "-linstep",
          FALSE,
          etSTR,
          { &evParams[0] },
          "Stepsizes (nm/step) for fixed increment linear sampling (put in quotes! \"1.0 2.3 5.1 "
          "-3.1\")" },
        { "-accdir",
          FALSE,
          etSTR,
          { &evParams[1] },
          "Directions for acceptance linear sampling - only sign counts! (put in quotes! \"-1 +1 "
          "-1.1\")" },
        { "-radstep",
          FALSE,
          etREAL,
          { &radstep },
          "Stepsize (nm/step) for fixed increment radius expansion" },
        { "-maxedsteps",
          FALSE,
          etINT,
          { &edi_params.maxedsteps },
          "Maximum number of steps per cycle" },
        { "-eqsteps",
          FALSE,
          etINT,
          { &eqSteps },
          "Number of steps to run without any perturbations " },
        { "-deltaF0", FALSE, etREAL, { &deltaF0 }, "Target destabilization energy for flooding" },
        { "-deltaF",
          FALSE,
          etREAL,
          { &deltaF },
          "Start deltaF with this parameter - default 0, nonzero values only needed for restart" },
        { "-tau",
          FALSE,
          etREAL,
          { &tau },
          "Coupling constant for adaption of flooding strength according to deltaF0, 0 = infinity "
          "i.e. constant flooding strength" },
        { "-Eflnull",
          FALSE,
          etREAL,
          { &constEfl },
          "The starting value of the flooding strength. The flooding strength is updated "
          "according to the adaptive flooding scheme. For a constant flooding strength use "
          "[TT]-tau[tt] 0. " },
        { "-T",
          FALSE,
          etREAL,
          { &T },
          "T is temperature, the value is needed if you want to do flooding " },
        { "-alpha",
          FALSE,
          etREAL,
          { &alpha },
          "Scale width of gaussian flooding potential with alpha^2 " },
        { "-restrain",
          FALSE,
          etBOOL,
          { &bRestrain },
          "Use the flooding potential with inverted sign -> effects as quasiharmonic restraining "
          "potential" },
        { "-hessian",
          FALSE,
          etBOOL,
          { &bHesse },
          "The eigenvectors and eigenvalues are from a Hessian matrix" },
        { "-harmonic",
          FALSE,
          etBOOL,
          { &bHarmonic },
          "The eigenvalues are interpreted as spring constant" },
        { "-constF",
          FALSE,
          etSTR,
          { &ConstForceStr },
          "Constant force flooding: manually set the forces for the eigenvectors selected with "
          "-flood "
          "(put in quotes! \"1.0 2.3 5.1 -3.1\"). No other flooding parameters are needed when "
          "specifying the forces directly." }
    };
#define NPA asize(pa)

    rvec*       xref1;
    int         nvec1, *eignr1  = nullptr;
    rvec *      xav1, **eigvec1 = nullptr;
    t_atoms*    atoms = nullptr;
    int         nav; /* Number of atoms in the average structure */
    char*       grpname;
    const char* indexfile;
    int         i;
    int *       index, *ifit;
    int         nfit;     /* Number of atoms in the reference/fit structure */
    int         ev_class; /* parameter _class i.e. evMON, evRADFIX etc. */
    int         nvecs;
    real*       eigval1 = nullptr; /* in V3.3 this is parameter of read_eigenvectors */

    const char* EdiFile;
    const char* TargetFile;
    const char* OriginFile;
    const char* EigvecFile;

    gmx_output_env_t* oenv;

    /*to read topology file*/
    t_topology top;
    PbcType    pbcType;
    matrix     topbox;
    rvec*      xtop;
    gmx_bool   bFit1;

    t_filenm fnm[] = { { efTRN, "-f", "eigenvec", ffREAD },  { efXVG, "-eig", "eigenval", ffOPTRD },
                       { efTPS, nullptr, nullptr, ffREAD },  { efNDX, nullptr, nullptr, ffOPTRD },
                       { efSTX, "-tar", "target", ffOPTRD }, { efSTX, "-ori", "origin", ffOPTRD },
                       { efEDI, "-o", "sam", ffWRITE } };
#define NFILE asize(fnm)
    edi_params.outfrq     = 100;
    edi_params.slope      = 0.0;
    edi_params.maxedsteps = 0;
    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, NPA, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    indexfile  = ftp2fn_null(efNDX, NFILE, fnm);
    EdiFile    = ftp2fn(efEDI, NFILE, fnm);
    TargetFile = opt2fn_null("-tar", NFILE, fnm);
    OriginFile = opt2fn_null("-ori", NFILE, fnm);


    for (ev_class = 0; ev_class < evNr; ++ev_class)
    {
        if (opt2parg_bSet(evOptions[ev_class], NPA, pa))
        {
            /*get list of eigenvectors*/
            nvecs = sscan_list(
                    &(listen[ev_class]), opt2parg_str(evOptions[ev_class], NPA, pa), evOptions[ev_class]);
            if (ev_class < evStepNr - 2)
            {
                /*if apropriate get list of stepsizes for these eigenvectors*/
                if (opt2parg_bSet(evStepOptions[ev_class], NPA, pa))
                {
                    evStepList[ev_class] = scan_vecparams(
                            opt2parg_str(evStepOptions[ev_class], NPA, pa), evStepOptions[ev_class], nvecs);
                }
                else /*if list is not given fill with zeros */
                {
                    snew(evStepList[ev_class], nvecs);
                    for (i = 0; i < nvecs; i++)
                    {
                        evStepList[ev_class][i] = 0.0;
                    }
                }
            }
            else if (ev_class == evRADFIX)
            {
                snew(evStepList[ev_class], nvecs);
                for (i = 0; i < nvecs; i++)
                {
                    evStepList[ev_class][i] = radstep;
                }
            }
            else if (ev_class == evFLOOD)
            {
                snew(evStepList[ev_class], nvecs);

                /* Are we doing constant force flooding? In that case, we read in
                 * the fproj values from the command line */
                if (opt2parg_bSet("-constF", NPA, pa))
                {
                    evStepList[ev_class] =
                            scan_vecparams(opt2parg_str("-constF", NPA, pa), "-constF", nvecs);
                }
            }
            else
            {
            } /*to avoid ambiguity   */
        }
        else /* if there are no eigenvectors for this option set list to zero */
        {
            listen[ev_class] = nullptr;
            snew(listen[ev_class], 1);
            listen[ev_class][0] = 0;
        }
    }

    /* print the interpreted list of eigenvectors - to give some feedback*/
    for (ev_class = 0; ev_class < evNr; ++ev_class)
    {
        printf("Eigenvector list %7s consists of the indices: ", evOptions[ev_class]);
        i = 0;
        while (listen[ev_class][i])
        {
            printf("%d ", listen[ev_class][i++]);
        }
        printf("\n");
    }

    EigvecFile = opt2fn("-f", NFILE, fnm);

    /*read eigenvectors from eigvec.trr*/
    read_eigenvectors(
            EigvecFile, &nav, &bFit1, &xref1, &edi_params.fitmas, &xav1, &edi_params.pcamas, &nvec1, &eignr1, &eigvec1, &eigval1);

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &xtop, nullptr, topbox, false);
    atoms = &top.atoms;


    printf("\nSelect an index group of %d elements that corresponds to the eigenvectors\n", nav);
    get_index(atoms, indexfile, 1, &i, &index, &grpname); /*if indexfile != NULL parameter 'atoms' is ignored */
    if (i != nav)
    {
        gmx_fatal(FARGS, "you selected a group with %d elements instead of %d", i, nav);
    }
    printf("\n");


    if (xref1 == nullptr)
    {
        if (bFit1)
        {
            /* if g_covar used different coordinate groups to fit and to do the PCA */
            printf("\nNote: the structure in %s should be the same\n"
                   "      as the one used for the fit in g_covar\n",
                   ftp2fn(efTPS, NFILE, fnm));
            printf("\nSelect the index group that was used for the least squares fit in g_covar\n");
        }
        else
        {
            printf("\nNote: Apparently no fitting was done in g_covar.\n"
                   "      However, you need to select a reference group for fitting in mdrun\n");
        }
        get_index(atoms, indexfile, 1, &nfit, &ifit, &grpname);
        snew(xref1, nfit);
        for (i = 0; i < nfit; i++)
        {
            copy_rvec(xtop[ifit[i]], xref1[i]);
        }
    }
    else
    {
        nfit = nav;
        ifit = index;
    }

    if (opt2parg_bSet("-constF", NPA, pa))
    {
        /* Constant force flooding is special: Most of the normal flooding
         * options are not needed. */
        edi_params.flood.bConstForce = TRUE;
    }
    else
    {
        /* For normal flooding read eigenvalues and store them in evSteplist[evFLOOD] */

        if (listen[evFLOOD][0] != 0)
        {
            read_eigenvalues(
                    listen[evFLOOD], opt2fn("-eig", NFILE, fnm), evStepList[evFLOOD], bHesse, kB * T, nav);
        }

        edi_params.flood.tau       = tau;
        edi_params.flood.deltaF0   = deltaF0;
        edi_params.flood.deltaF    = deltaF;
        edi_params.presteps        = eqSteps;
        edi_params.flood.kT        = kB * T;
        edi_params.flood.bHarmonic = bHarmonic;
        if (bRestrain)
        {
            /* Trick: invert sign of Efl and alpha2 then this will give the same sign in the exponential and inverted sign outside */
            edi_params.flood.constEfl = -constEfl;
            edi_params.flood.alpha2   = -gmx::square(alpha);
        }
        else
        {
            edi_params.flood.constEfl = constEfl;
            edi_params.flood.alpha2   = gmx::square(alpha);
        }
    }

    edi_params.ned = nav;

    /*number of system atoms  */
    edi_params.nini = atoms->nr;


    /*store reference and average structure in edi_params*/
    make_t_edx(&edi_params.sref, nfit, xref1, ifit);
    make_t_edx(&edi_params.sav, nav, xav1, index);


    /* Store target positions in edi_params */
    if (opt2bSet("-tar", NFILE, fnm))
    {
        if (0 != listen[evFLOOD][0])
        {
            fprintf(stderr,
                    "\nNote: Providing a TARGET structure has no effect when using flooding.\n"
                    "      You may want to use -ori to define the flooding potential center.\n\n");
        }
        get_structure(atoms, indexfile, TargetFile, &edi_params.star, nfit, ifit, nav, index);
    }
    else
    {
        make_t_edx(&edi_params.star, 0, nullptr, index);
    }

    /* Store origin positions */
    if (opt2bSet("-ori", NFILE, fnm))
    {
        get_structure(atoms, indexfile, OriginFile, &edi_params.sori, nfit, ifit, nav, index);
    }
    else
    {
        make_t_edx(&edi_params.sori, 0, nullptr, index);
    }

    /* Write edi-file */
    write_the_whole_thing(EdiFile, &edi_params, eigvec1, nvec1, listen, evStepList);

    return 0;
}
