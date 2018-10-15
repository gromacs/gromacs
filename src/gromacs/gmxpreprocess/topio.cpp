/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "topio.h"

#include <cassert>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include <sys/types.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/gmxcpp.h"
#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/topshake.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/vsite_parm.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

#define OPENDIR     '[' /* starting sign for directive */
#define CLOSEDIR    ']' /* ending sign for directive   */

static void free_nbparam(t_nbparam **param, int nr)
{
    int i;

    assert(param);
    for (i = 0; i < nr; i++)
    {
        assert(param[i]);
        sfree(param[i]);
    }
    sfree(param);
}

static int copy_nbparams(t_nbparam **param, int ftype, t_params *plist, int nr)
{
    int i, j, f;
    int nrfp, ncopy;

    nrfp = NRFP(ftype);

    ncopy = 0;
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j <= i; j++)
        {
            assert(param);
            if (param[i][j].bSet)
            {
                for (f = 0; f < nrfp; f++)
                {
                    plist->param[nr*i+j].c[f] = param[i][j].c[f];
                    plist->param[nr*j+i].c[f] = param[i][j].c[f];
                }
                ncopy++;
            }
        }
    }

    return ncopy;
}

static void gen_pairs(t_params *nbs, t_params *pairs, real fudge, int comb)
{
    int     i, j, ntp, nrfp, nrfpA, nrfpB, nnn;
    real    scaling;
    ntp       = nbs->nr;
    nnn       = static_cast<int>(std::sqrt(static_cast<double>(ntp)));
    GMX_ASSERT(nnn * nnn == ntp, "Number of pairs of generated non-bonded parameters should be a perfect square");
    nrfp      = NRFP(F_LJ);
    nrfpA     = interaction_function[F_LJ14].nrfpA;
    nrfpB     = interaction_function[F_LJ14].nrfpB;
    pairs->nr = ntp;

    if ((nrfp  != nrfpA) || (nrfpA != nrfpB))
    {
        gmx_incons("Number of force parameters in gen_pairs wrong");
    }

    fprintf(stderr, "Generating 1-4 interactions: fudge = %g\n", fudge);
    snew(pairs->param, pairs->nr);
    for (i = 0; (i < ntp); i++)
    {
        /* Copy param.a */
        pairs->param[i].a[0] = i / nnn;
        pairs->param[i].a[1] = i % nnn;
        /* Copy normal and FEP parameters and multiply by fudge factor */



        for (j = 0; (j < nrfp); j++)
        {
            /* If we are using sigma/epsilon values, only the epsilon values
             * should be scaled, but not sigma.
             * The sigma values have even indices 0,2, etc.
             */
            if ((comb == eCOMB_ARITHMETIC || comb == eCOMB_GEOM_SIG_EPS) && (j%2 == 0))
            {
                scaling = 1.0;
            }
            else
            {
                scaling = fudge;
            }

            pairs->param[i].c[j]      = scaling*nbs->param[i].c[j];
            /* NOTE: this should be cleat to the compiler, but some gcc 5.2 versions
             *  issue false positive warnings for the pairs->param.c[] indexing below.
             */
            assert(2*nrfp <= MAXFORCEPARAM);
            pairs->param[i].c[nrfp+j] = scaling*nbs->param[i].c[j];
        }
    }
}

double check_mol(const gmx_mtop_t *mtop, warninp_t wi)
{
    char     buf[256];
    int      i, ri, pt;
    double   q;
    real     m, mB;

    /* Check mass and charge */
    q = 0.0;

    for (const gmx_molblock_t &molb : mtop->molblock)
    {
        const t_atoms *atoms = &mtop->moltype[molb.type].atoms;
        for (i = 0; (i < atoms->nr); i++)
        {
            q += molb.nmol*atoms->atom[i].q;
            m  = atoms->atom[i].m;
            mB = atoms->atom[i].mB;
            pt = atoms->atom[i].ptype;
            /* If the particle is an atom or a nucleus it must have a mass,
             * else, if it is a shell, a vsite or a bondshell it can have mass zero
             */
            if (((m <= 0.0) || (mB <= 0.0)) && ((pt == eptAtom) || (pt == eptNucleus)))
            {
                ri = atoms->atom[i].resind;
                sprintf(buf, "atom %s (Res %s-%d) has mass %g (state A) / %g (state B)\n",
                        *(atoms->atomname[i]),
                        *(atoms->resinfo[ri].name),
                        atoms->resinfo[ri].nr,
                        m, mB);
                warning_error(wi, buf);
            }
            else
            if (((m != 0) || (mB != 0)) && (pt == eptVSite))
            {
                ri = atoms->atom[i].resind;
                sprintf(buf, "virtual site %s (Res %s-%d) has non-zero mass %g (state A) / %g (state B)\n"
                        "     Check your topology.\n",
                        *(atoms->atomname[i]),
                        *(atoms->resinfo[ri].name),
                        atoms->resinfo[ri].nr,
                        m, mB);
                warning_error(wi, buf);
                /* The following statements make LINCS break! */
                /* atoms->atom[i].m=0; */
            }
        }
    }
    return q;
}

/*! \brief Returns the rounded charge of a molecule, when close to integer, otherwise returns the original charge.
 *
 * The results of this routine are only used for checking and for
 * printing warning messages. Thus we can assume that charges of molecules
 * should be integer. If the user wanted non-integer molecular charge,
 * an undesired warning is printed and the user should use grompp -maxwarn 1.
 *
 * \param qMol     The total, unrounded, charge of the molecule
 * \param sumAbsQ  The sum of absolute values of the charges, used for determining the tolerance for the rounding.
 */
static double roundedMoleculeCharge(double qMol, double sumAbsQ)
{
    /* We use a tolerance of 1e-6 for inaccuracies beyond the 6th decimal
     * of the charges for ascii float truncation in the topology files.
     * Although the summation here uses double precision, the charges
     * are read and stored in single precision when real=float. This can
     * lead to rounding errors of half the least significant bit.
     * Note that, unfortunately, we can not assume addition of random
     * rounding errors. It is not entirely unlikely that many charges
     * have a near half-bit rounding error with the same sign.
     */
    double tolAbs = 1e-6;
    double tol    = std::max(tolAbs, 0.5*GMX_REAL_EPS*sumAbsQ);
    double qRound = std::round(qMol);
    if (std::abs(qMol - qRound) <= tol)
    {
        return qRound;
    }
    else
    {
        return qMol;
    }
}

static void sum_q(const t_atoms *atoms, int numMols,
                  double *qTotA, double *qTotB)
{
    /* sum charge */
    double qmolA    = 0;
    double qmolB    = 0;
    double sumAbsQA = 0;
    double sumAbsQB = 0;
    for (int i = 0; i < atoms->nr; i++)
    {
        qmolA    += atoms->atom[i].q;
        qmolB    += atoms->atom[i].qB;
        sumAbsQA += std::abs(atoms->atom[i].q);
        sumAbsQB += std::abs(atoms->atom[i].qB);
    }

    *qTotA += numMols*roundedMoleculeCharge(qmolA, sumAbsQA);
    *qTotB += numMols*roundedMoleculeCharge(qmolB, sumAbsQB);
}

static void get_nbparm(char *nb_str, char *comb_str, int *nb, int *comb,
                       warninp_t wi)
{
    int  i;
    char warn_buf[STRLEN];

    *nb   = -1;
    for (i = 1; (i < eNBF_NR); i++)
    {
        if (gmx_strcasecmp(nb_str, enbf_names[i]) == 0)
        {
            *nb = i;
        }
    }
    if (*nb == -1)
    {
        *nb = strtol(nb_str, nullptr, 10);
    }
    if ((*nb < 1) || (*nb >= eNBF_NR))
    {
        sprintf(warn_buf, "Invalid nonbond function selector '%s' using %s",
                nb_str, enbf_names[1]);
        warning_error(wi, warn_buf);
        *nb = 1;
    }
    *comb = -1;
    for (i = 1; (i < eCOMB_NR); i++)
    {
        if (gmx_strcasecmp(comb_str, ecomb_names[i]) == 0)
        {
            *comb = i;
        }
    }
    if (*comb == -1)
    {
        *comb = strtol(comb_str, nullptr, 10);
    }
    if ((*comb < 1) || (*comb >= eCOMB_NR))
    {
        sprintf(warn_buf, "Invalid combination rule selector '%s' using %s",
                comb_str, ecomb_names[1]);
        warning_error(wi, warn_buf);
        *comb = 1;
    }
}

static char ** cpp_opts(const char *define, const char *include,
                        warninp_t wi)
{
    int         n, len;
    int         ncppopts = 0;
    const char *cppadds[2];
    char      **cppopts   = nullptr;
    const char *option[2] = { "-D", "-I" };
    const char *nopt[2]   = { "define", "include" };
    const char *ptr;
    const char *rptr;
    char       *buf;
    char        warn_buf[STRLEN];

    cppadds[0] = define;
    cppadds[1] = include;
    for (n = 0; (n < 2); n++)
    {
        if (cppadds[n])
        {
            ptr = cppadds[n];
            while (*ptr != '\0')
            {
                while ((*ptr != '\0') && isspace(*ptr))
                {
                    ptr++;
                }
                rptr = ptr;
                while ((*rptr != '\0') && !isspace(*rptr))
                {
                    rptr++;
                }
                len = (rptr - ptr);
                if (len > 2)
                {
                    snew(buf, (len+1));
                    strncpy(buf, ptr, len);
                    if (strstr(ptr, option[n]) != ptr)
                    {
                        set_warning_line(wi, "mdp file", -1);
                        sprintf(warn_buf, "Malformed %s option %s", nopt[n], buf);
                        warning(wi, warn_buf);
                    }
                    else
                    {
                        srenew(cppopts, ++ncppopts);
                        cppopts[ncppopts-1] = gmx_strdup(buf);
                    }
                    sfree(buf);
                    ptr = rptr;
                }
            }
        }
    }
    srenew(cppopts, ++ncppopts);
    cppopts[ncppopts-1] = nullptr;

    return cppopts;
}


static void make_atoms_sys(const std::vector<gmx_molblock_t> &molblock,
                           const t_molinfo                   *molinfo,
                           t_atoms                           *atoms)
{
    atoms->nr   = 0;
    atoms->atom = nullptr;

    for (const gmx_molblock_t &molb : molblock)
    {
        const t_atoms &mol_atoms = molinfo[molb.type].atoms;

        srenew(atoms->atom, atoms->nr + molb.nmol*mol_atoms.nr);

        for (int m = 0; m < molb.nmol; m++)
        {
            for (int a = 0; a < mol_atoms.nr; a++)
            {
                atoms->atom[atoms->nr++] = mol_atoms.atom[a];
            }
        }
    }
}


static char **read_topol(const char *infile, const char *outfile,
                         const char *define, const char *include,
                         t_symtab    *symtab,
                         gpp_atomtype_t atype,
                         int         *nrmols,
                         t_molinfo   **molinfo,
                         t_molinfo   **intermolecular_interactions,
                         t_params    plist[],
                         int         *combination_rule,
                         double      *reppow,
                         t_gromppopts *opts,
                         real        *fudgeQQ,
                         std::vector<gmx_molblock_t> *molblock,
                         bool       *ffParametrizedWithHBondConstraints,
                         bool        bFEP,
                         bool        bZero,
                         bool        usingFullRangeElectrostatics,
                         warninp_t       wi)
{
    FILE                 *out;
    int                   i, sl, nb_funct;
    char                 *pline = nullptr, **title = nullptr;
    char                  line[STRLEN], errbuf[256], comb_str[256], nb_str[256];
    char                  genpairs[32];
    char                 *dirstr, *dummy2;
    int                   nrcopies, nmol, nscan, ncombs, ncopy;
    double                fLJ, fQQ, fPOW;
    t_molinfo            *mi0   = nullptr;
    DirStack             *DS;
    directive             d, newd;
    t_nbparam           **nbparam, **pair;
    gmx::ExclusionBlocks *exclusionBlocks;
    real                  fudgeLJ = -1;    /* Multiplication factor to generate 1-4 from LJ */
    bool                  bReadDefaults, bReadMolType, bGenPairs, bWarn_copy_A_B;
    double                qt = 0, qBt = 0; /* total charge */
    t_bond_atomtype       batype;
    int                   lastcg = -1;
    int                   dcatt  = -1, nmol_couple;
    /* File handling variables */
    int                   status;
    bool                  done;
    gmx_cpp_t             handle;
    char                 *tmp_line = nullptr;
    char                  warn_buf[STRLEN];
    const char           *floating_point_arithmetic_tip =
        "Total charge should normally be an integer. See\n"
        "http://www.gromacs.org/Documentation/Floating_Point_Arithmetic\n"
        "for discussion on how close it should be to an integer.\n";
    /* We need to open the output file before opening the input file,
     * because cpp_open_file can change the current working directory.
     */
    if (outfile)
    {
        out = gmx_fio_fopen(outfile, "w");
    }
    else
    {
        out = nullptr;
    }

    /* open input file */
    auto cpp_opts_return = cpp_opts(define, include, wi);
    status = cpp_open_file(infile, &handle, cpp_opts_return);
    if (status != 0)
    {
        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
    }

    /* some local variables */
    DS_Init(&DS);                /* directive stack	*/
    nmol            = 0;         /* no molecules yet...	*/
    d               = d_invalid; /* first thing should be a directive */
    nbparam         = nullptr;   /* The temporary non-bonded matrix */
    pair            = nullptr;   /* The temporary pair interaction matrix */
    exclusionBlocks = nullptr;   /* the extra exclusions */
    nb_funct        = F_LJ;

    *reppow  = 12.0;      /* Default value for repulsion power     */

    *intermolecular_interactions = nullptr;

    /* Init the number of CMAP torsion angles  and grid spacing */
    plist[F_CMAP].grid_spacing = 0;
    plist[F_CMAP].nc           = 0;

    bWarn_copy_A_B = bFEP;

    batype = init_bond_atomtype();
    /* parse the actual file */
    bReadDefaults = FALSE;
    bGenPairs     = FALSE;
    bReadMolType  = FALSE;
    nmol_couple   = 0;

    do
    {
        status = cpp_read_line(&handle, STRLEN, line);
        done   = (status == eCPP_EOF);
        if (!done)
        {
            if (status != eCPP_OK)
            {
                gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
            }
            else if (out)
            {
                fprintf(out, "%s\n", line);
            }

            set_warning_line(wi, cpp_cur_file(&handle), cpp_cur_linenr(&handle));

            pline = gmx_strdup(line);

            /* Strip trailing '\' from pline, if it exists */
            sl = strlen(pline);
            if ((sl > 0) && (pline[sl-1] == CONTINUE))
            {
                pline[sl-1] = ' ';
            }

            /* build one long line from several fragments - necessary for CMAP */
            while (continuing(line))
            {
                status = cpp_read_line(&handle, STRLEN, line);
                set_warning_line(wi, cpp_cur_file(&handle), cpp_cur_linenr(&handle));

                /* Since we depend on the '\' being present to continue to read, we copy line
                 * to a tmp string, strip the '\' from that string, and cat it to pline
                 */
                tmp_line = gmx_strdup(line);

                sl = strlen(tmp_line);
                if ((sl > 0) && (tmp_line[sl-1] == CONTINUE))
                {
                    tmp_line[sl-1] = ' ';
                }

                done = (status == eCPP_EOF);
                if (!done)
                {
                    if (status != eCPP_OK)
                    {
                        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
                    }
                    else if (out)
                    {
                        fprintf(out, "%s\n", line);
                    }
                }

                srenew(pline, strlen(pline)+strlen(tmp_line)+1);
                strcat(pline, tmp_line);
                sfree(tmp_line);
            }

            /* skip trailing and leading spaces and comment text */
            strip_comment (pline);
            trim (pline);

            /* if there is something left... */
            if (static_cast<int>(strlen(pline)) > 0)
            {
                if (pline[0] == OPENDIR)
                {
                    /* A directive on this line: copy the directive
                     * without the brackets into dirstr, then
                     * skip spaces and tabs on either side of directive
                     */
                    dirstr = gmx_strdup((pline+1));
                    if ((dummy2 = strchr (dirstr, CLOSEDIR)) != nullptr)
                    {
                        (*dummy2) = 0;
                    }
                    trim (dirstr);

                    if ((newd = str2dir(dirstr)) == d_invalid)
                    {
                        sprintf(errbuf, "Invalid directive %s", dirstr);
                        warning_error(wi, errbuf);
                    }
                    else
                    {
                        /* Directive found */
                        if (DS_Check_Order (DS, newd))
                        {
                            DS_Push (&DS, newd);
                            d = newd;
                        }
                        else
                        {
                            /* we should print here which directives should have
                               been present, and which actually are */
                            gmx_fatal(FARGS, "%s\nInvalid order for directive %s",
                                      cpp_error(&handle, eCPP_SYNTAX), dir2str(newd));
                            /* d = d_invalid; */
                        }

                        if (d == d_intermolecular_interactions)
                        {
                            if (*intermolecular_interactions == nullptr)
                            {
                                /* We (mis)use the moleculetype processing
                                 * to process the intermolecular interactions
                                 * by making a "molecule" of the size of the system.
                                 */
                                snew(*intermolecular_interactions, 1);
                                init_molinfo(*intermolecular_interactions);
                                mi0 = *intermolecular_interactions;
                                make_atoms_sys(*molblock, *molinfo,
                                               &mi0->atoms);
                            }
                        }
                    }
                    sfree(dirstr);
                }
                else if (d != d_invalid)
                {
                    /* Not a directive, just a plain string
                     * use a gigantic switch to decode,
                     * if there is a valid directive!
                     */
                    switch (d)
                    {
                        case d_defaults:
                            if (bReadDefaults)
                            {
                                gmx_fatal(FARGS, "%s\nFound a second defaults directive.\n",
                                          cpp_error(&handle, eCPP_SYNTAX));
                            }
                            bReadDefaults = TRUE;
                            nscan         = sscanf(pline, "%s%s%s%lf%lf%lf",
                                                   nb_str, comb_str, genpairs, &fLJ, &fQQ, &fPOW);
                            if (nscan < 2)
                            {
                                too_few(wi);
                            }
                            else
                            {
                                bGenPairs = FALSE;
                                fudgeLJ   = 1.0;
                                *fudgeQQ  = 1.0;

                                get_nbparm(nb_str, comb_str, &nb_funct, combination_rule, wi);
                                if (nscan >= 3)
                                {
                                    bGenPairs = (gmx_strncasecmp(genpairs, "Y", 1) == 0);
                                    if (nb_funct != eNBF_LJ && bGenPairs)
                                    {
                                        gmx_fatal(FARGS, "Generating pair parameters is only supported with LJ non-bonded interactions");
                                    }
                                }
                                if (nscan >= 4)
                                {
                                    fudgeLJ   = fLJ;
                                }
                                if (nscan >= 5)
                                {
                                    *fudgeQQ  = fQQ;
                                }
                                if (nscan >= 6)
                                {
                                    *reppow   = fPOW;
                                }
                            }
                            nb_funct = ifunc_index(d_nonbond_params, nb_funct);

                            break;
                        case d_atomtypes:
                            push_at(symtab, atype, batype, pline, nb_funct,
                                    &nbparam, bGenPairs ? &pair : nullptr, wi);
                            break;

                        case d_bondtypes:
                            push_bt(d, plist, 2, nullptr, batype, pline, wi);
                            break;
                        case d_constrainttypes:
                            push_bt(d, plist, 2, nullptr, batype, pline, wi);
                            break;
                        case d_pairtypes:
                            if (bGenPairs)
                            {
                                push_nbt(d, pair, atype, pline, F_LJ14, wi);
                            }
                            else
                            {
                                push_bt(d, plist, 2, atype, nullptr, pline, wi);
                            }
                            break;
                        case d_angletypes:
                            push_bt(d, plist, 3, nullptr, batype, pline, wi);
                            break;
                        case d_dihedraltypes:
                            /* Special routine that can read both 2 and 4 atom dihedral definitions. */
                            push_dihedraltype(d, plist, batype, pline, wi);
                            break;

                        case d_nonbond_params:
                            push_nbt(d, nbparam, atype, pline, nb_funct, wi);
                            break;
                        /*
                           case d_blocktype:
                           nblock++;
                           srenew(block,nblock);
                           srenew(blockinfo,nblock);
                           blk0=&(block[nblock-1]);
                           bi0=&(blockinfo[nblock-1]);
                           init_top(blk0);
                           init_molinfo(bi0);
                           push_molt(symtab,bi0,pline);
                           break;
                         */

                        case d_implicit_genborn_params:
                            // Skip this line, so old topologies with
                            // GB parameters can be read.
                            break;

                        case d_implicit_surface_params:
                            // Skip this line, so that any topologies
                            // with surface parameters can be read
                            // (even though these were never formally
                            // supported).
                            break;

                        case d_cmaptypes:
                            push_cmaptype(d, plist, 5, atype, batype, pline, wi);
                            break;

                        case d_moleculetype:
                        {
                            if (!bReadMolType)
                            {
                                int ntype;
                                if (opts->couple_moltype != nullptr &&
                                    (opts->couple_lam0 == ecouplamNONE ||
                                     opts->couple_lam0 == ecouplamQ ||
                                     opts->couple_lam1 == ecouplamNONE ||
                                     opts->couple_lam1 == ecouplamQ))
                                {
                                    dcatt = add_atomtype_decoupled(symtab, atype,
                                                                   &nbparam, bGenPairs ? &pair : nullptr);
                                }
                                ntype  = get_atomtype_ntypes(atype);
                                ncombs = (ntype*(ntype+1))/2;
                                generate_nbparams(*combination_rule, nb_funct, &(plist[nb_funct]), atype, wi);
                                ncopy = copy_nbparams(nbparam, nb_funct, &(plist[nb_funct]),
                                                      ntype);
                                fprintf(stderr, "Generated %d of the %d non-bonded parameter combinations\n", ncombs-ncopy, ncombs);
                                free_nbparam(nbparam, ntype);
                                if (bGenPairs)
                                {
                                    gen_pairs(&(plist[nb_funct]), &(plist[F_LJ14]), fudgeLJ, *combination_rule);
                                    ncopy = copy_nbparams(pair, nb_funct, &(plist[F_LJ14]),
                                                          ntype);
                                    fprintf(stderr, "Generated %d of the %d 1-4 parameter combinations\n", ncombs-ncopy, ncombs);
                                    free_nbparam(pair, ntype);
                                }
                                /* Copy GBSA parameters to atomtype array? */

                                bReadMolType = TRUE;
                            }

                            push_molt(symtab, &nmol, molinfo, pline, wi);
                            srenew(exclusionBlocks, nmol);
                            exclusionBlocks[nmol-1].nr      = 0;
                            mi0                             = &((*molinfo)[nmol-1]);
                            mi0->atoms.haveMass             = TRUE;
                            mi0->atoms.haveCharge           = TRUE;
                            mi0->atoms.haveType             = TRUE;
                            mi0->atoms.haveBState           = TRUE;
                            mi0->atoms.havePdbInfo          = FALSE;
                            break;
                        }
                        case d_atoms:
                            push_atom(symtab, &(mi0->cgs), &(mi0->atoms), atype, pline, &lastcg, wi);
                            break;

                        case d_pairs:
                            push_bond(d, plist, mi0->plist, &(mi0->atoms), atype, pline, FALSE,
                                      bGenPairs, *fudgeQQ, bZero, &bWarn_copy_A_B, wi);
                            break;
                        case d_pairs_nb:
                            push_bond(d, plist, mi0->plist, &(mi0->atoms), atype, pline, FALSE,
                                      FALSE, 1.0, bZero, &bWarn_copy_A_B, wi);
                            break;

                        case d_vsites2:
                        case d_vsites3:
                        case d_vsites4:
                        case d_bonds:
                        case d_angles:
                        case d_constraints:
                        case d_settles:
                        case d_position_restraints:
                        case d_angle_restraints:
                        case d_angle_restraints_z:
                        case d_distance_restraints:
                        case d_orientation_restraints:
                        case d_dihedral_restraints:
                        case d_dihedrals:
                        case d_polarization:
                        case d_water_polarization:
                        case d_thole_polarization:
                            push_bond(d, plist, mi0->plist, &(mi0->atoms), atype, pline, TRUE,
                                      bGenPairs, *fudgeQQ, bZero, &bWarn_copy_A_B, wi);
                            break;
                        case d_cmap:
                            push_cmap(d, plist, mi0->plist, &(mi0->atoms), atype, pline, wi);
                            break;

                        case d_vsitesn:
                            push_vsitesn(d, mi0->plist, &(mi0->atoms), pline, wi);
                            break;
                        case d_exclusions:
                            GMX_ASSERT(exclusionBlocks, "exclusionBlocks must always be allocated so exclusions can be processed");
                            if (!exclusionBlocks[nmol-1].nr)
                            {
                                initExclusionBlocks(&(exclusionBlocks[nmol-1]), mi0->atoms.nr);
                            }
                            push_excl(pline, &(exclusionBlocks[nmol-1]), wi);
                            break;
                        case d_system:
                            trim(pline);
                            title = put_symtab(symtab, pline);
                            break;
                        case d_molecules:
                        {
                            int      whichmol;
                            bool     bCouple;

                            push_mol(nmol, *molinfo, pline, &whichmol, &nrcopies, wi);
                            mi0 = &((*molinfo)[whichmol]);
                            molblock->resize(molblock->size() + 1);
                            molblock->back().type = whichmol;
                            molblock->back().nmol = nrcopies;

                            bCouple = (opts->couple_moltype != nullptr &&
                                       (gmx_strcasecmp("system", opts->couple_moltype) == 0 ||
                                        strcmp(*(mi0->name), opts->couple_moltype) == 0));
                            if (bCouple)
                            {
                                nmol_couple += nrcopies;
                            }

                            if (mi0->atoms.nr == 0)
                            {
                                gmx_fatal(FARGS, "Molecule type '%s' contains no atoms",
                                          *mi0->name);
                            }
                            fprintf(stderr,
                                    "Excluding %d bonded neighbours molecule type '%s'\n",
                                    mi0->nrexcl, *mi0->name);
                            sum_q(&mi0->atoms, nrcopies, &qt, &qBt);
                            if (!mi0->bProcessed)
                            {
                                t_nextnb nnb;
                                generate_excl(mi0->nrexcl,
                                              mi0->atoms.nr,
                                              mi0->plist,
                                              &nnb,
                                              &(mi0->excls));
                                gmx::mergeExclusions(&(mi0->excls), &(exclusionBlocks[whichmol]));
                                gmx::doneExclusionBlocks(&(exclusionBlocks[whichmol]));
                                make_shake(mi0->plist, &mi0->atoms, opts->nshake);



                                done_nnb(&nnb);

                                if (bCouple)
                                {
                                    convert_moltype_couple(mi0, dcatt, *fudgeQQ,
                                                           opts->couple_lam0, opts->couple_lam1,
                                                           opts->bCoupleIntra,
                                                           nb_funct, &(plist[nb_funct]), wi);
                                }
                                stupid_fill_block(&mi0->mols, mi0->atoms.nr, TRUE);
                                mi0->bProcessed = TRUE;
                            }
                            break;
                        }
                        default:
                            fprintf (stderr, "case: %d\n", static_cast<int>(d));
                            gmx_incons("unknown directive");
                    }
                }
            }
            sfree(pline);
            pline = nullptr;
        }
    }
    while (!done);
    sfree(cpp_opts_return);

    if (out)
    {
        gmx_fio_fclose(out);
    }

    /* List of GROMACS define names for force fields that have been
     * parametrized using constraints involving hydrogens only.
     *
     * We should avoid hardcoded names, but this is hopefully only
     * needed temparorily for discouraging use of constraints=all-bonds.
     */
    const std::array<std::string, 3> ffDefines = {
        "_FF_AMBER",
        "_FF_CHARMM",
        "_FF_OPLSAA"
    };
    *ffParametrizedWithHBondConstraints = false;
    for (const std::string &ffDefine : ffDefines)
    {
        if (cpp_find_define(&handle, ffDefine))
        {
            *ffParametrizedWithHBondConstraints = true;
        }
    }

    cpp_done(handle);

    if (opts->couple_moltype)
    {
        if (nmol_couple == 0)
        {
            gmx_fatal(FARGS, "Did not find any molecules of type '%s' for coupling",
                      opts->couple_moltype);
        }
        fprintf(stderr, "Coupling %d copies of molecule type '%s'\n",
                nmol_couple, opts->couple_moltype);
    }

    /* this is not very clean, but fixes core dump on empty system name */
    if (!title)
    {
        title = put_symtab(symtab, "");
    }

    if (fabs(qt) > 1e-4)
    {
        sprintf(warn_buf, "System has non-zero total charge: %.6f\n%s\n", qt, floating_point_arithmetic_tip);
        warning_note(wi, warn_buf);
    }
    if (fabs(qBt) > 1e-4 && !gmx_within_tol(qBt, qt, 1e-6))
    {
        sprintf(warn_buf, "State B has non-zero total charge: %.6f\n%s\n", qBt, floating_point_arithmetic_tip);
        warning_note(wi, warn_buf);
    }
    if (usingFullRangeElectrostatics && (fabs(qt) > 1e-4 || fabs(qBt) > 1e-4))
    {
        warning(wi, "You are using Ewald electrostatics in a system with net charge. This can lead to severe artifacts, such as ions moving into regions with low dielectric, due to the uniform background charge. We suggest to neutralize your system with counter ions, possibly in combination with a physiological salt concentration.");
        please_cite(stdout, "Hub2014a");
    }

    DS_Done (&DS);
    for (i = 0; i < nmol; i++)
    {
        gmx::doneExclusionBlocks(&(exclusionBlocks[i]));
    }
    free(exclusionBlocks);

    done_bond_atomtype(&batype);

    if (*intermolecular_interactions != nullptr)
    {
        sfree(mi0->atoms.atom);
    }

    *nrmols = nmol;

    return title;
}

char **do_top(bool                          bVerbose,
              const char                   *topfile,
              const char                   *topppfile,
              t_gromppopts                 *opts,
              bool                          bZero,
              t_symtab                     *symtab,
              t_params                      plist[],
              int                          *combination_rule,
              double                       *repulsion_power,
              real                         *fudgeQQ,
              gpp_atomtype_t                atype,
              int                          *nrmols,
              t_molinfo                   **molinfo,
              t_molinfo                   **intermolecular_interactions,
              const t_inputrec             *ir,
              std::vector<gmx_molblock_t>  *molblock,
              bool                         *ffParametrizedWithHBondConstraints,
              warninp_t                     wi)
{
    /* Tmpfile might contain a long path */
    const char *tmpfile;
    char      **title;

    if (topppfile)
    {
        tmpfile = topppfile;
    }
    else
    {
        tmpfile = nullptr;
    }

    if (bVerbose)
    {
        printf("processing topology...\n");
    }
    title = read_topol(topfile, tmpfile, opts->define, opts->include,
                       symtab, atype,
                       nrmols, molinfo, intermolecular_interactions,
                       plist, combination_rule, repulsion_power,
                       opts, fudgeQQ, molblock,
                       ffParametrizedWithHBondConstraints,
                       ir->efep != efepNO, bZero,
                       EEL_FULL(ir->coulombtype), wi);

    if ((*combination_rule != eCOMB_GEOMETRIC) &&
        (ir->vdwtype == evdwUSER))
    {
        warning(wi, "Using sigma/epsilon based combination rules with"
                " user supplied potential function may produce unwanted"
                " results");
    }

    return title;
}

/*! \brief
 * Generate exclusion lists for QM/MM.
 *
 * This routine updates the exclusion lists for QM atoms in order to include all other QM
 * atoms of this molecule. Moreover, this routine replaces bonds between QM atoms with
 * CONNBOND and, when MiMiC is not used, removes bonded interactions between QM and link atoms.
 * Finally, in case if MiMiC QM/MM is used - charges of QM atoms are set to 0
 *
 * @param molt molecule type with QM atoms
 * @param grpnr group informatio
 * @param ir input record
 * @param qmmmMode QM/MM mode switch: original/MiMiC
 */
static void generate_qmexcl_moltype(gmx_moltype_t *molt, const unsigned char *grpnr,
                                    t_inputrec *ir, GmxQmmmMode qmmmMode)
{
    /* This routine expects molt->ilist to be of size F_NRE and ordered. */

    /* generates the exclusions between the individual QM atoms, as
     * these interactions should be handled by the QM subroutines and
     * not by the gromacs routines
     */
    int       qm_max = 0, qm_nr = 0, link_nr = 0, link_max = 0;
    int      *qm_arr = nullptr, *link_arr = nullptr;
    bool     *bQMMM, *blink;

    /* First we search and select the QM atoms in an qm_arr array that
     * we use to create the exclusions.
     *
     * we take the possibility into account that a user has defined more
     * than one QM group:
     *
     * for that we also need to do this an ugly work-about just in case
     * the QM group contains the entire system...
     */

    /* we first search for all the QM atoms and put them in an array
     */
    for (int j = 0; j < ir->opts.ngQM; j++)
    {
        for (int i = 0; i < molt->atoms.nr; i++)
        {
            if (qm_nr >= qm_max)
            {
                qm_max += 100;
                srenew(qm_arr, qm_max);
            }
            if ((grpnr ? grpnr[i] : 0) == j)
            {
                qm_arr[qm_nr++]        = i;
                molt->atoms.atom[i].q  = 0.0;
                molt->atoms.atom[i].qB = 0.0;
            }
        }
    }
    /* bQMMM[..] is an array containin TRUE/FALSE for atoms that are
     * QM/not QM. We first set all elements to false. Afterwards we use
     * the qm_arr to change the elements corresponding to the QM atoms
     * to TRUE.
     */
    snew(bQMMM, molt->atoms.nr);
    for (int i = 0; i < molt->atoms.nr; i++)
    {
        bQMMM[i] = FALSE;
    }
    for (int i = 0; i < qm_nr; i++)
    {
        bQMMM[qm_arr[i]] = TRUE;
    }

    /* We remove all bonded interactions (i.e. bonds,
     * angles, dihedrals, 1-4's), involving the QM atoms. The way they
     * are removed is as follows: if the interaction invloves 2 atoms,
     * it is removed if both atoms are QMatoms. If it involves 3 atoms,
     * it is removed if at least two of the atoms are QM atoms, if the
     * interaction involves 4 atoms, it is removed if there are at least
     * 2 QM atoms.  Since this routine is called once before any forces
     * are computed, the top->idef.il[N].iatom[] array (see idef.h) can
     * be rewritten at this poitn without any problem. 25-9-2002 */

    /* first check whether we already have CONNBONDS.
     * Note that if we don't, we don't add a param entry and set ftype=0,
     * which is ok, since CONNBONDS does not use parameters.
     */
    int ftype_connbond = 0;
    int ind_connbond   = 0;
    if (molt->ilist[F_CONNBONDS].size() != 0)
    {
        fprintf(stderr, "nr. of CONNBONDS present already: %d\n",
                molt->ilist[F_CONNBONDS].size()/3);
        ftype_connbond = molt->ilist[F_CONNBONDS].iatoms[0];
        ind_connbond   = molt->ilist[F_CONNBONDS].size();
    }
    /* now we delete all bonded interactions, except the ones describing
     * a chemical bond. These are converted to CONNBONDS
     */
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (!(interaction_function[ftype].flags & IF_BOND) ||
            ftype == F_CONNBONDS)
        {
            continue;
        }
        int nratoms = interaction_function[ftype].nratoms;
        int j       = 0;
        while (j < molt->ilist[ftype].size())
        {
            bool bexcl;

            if (nratoms == 2)
            {
                /* Remove an interaction between two atoms when both are
                 * in the QM region. Note that we don't have to worry about
                 * link atoms here, as they won't have 2-atom interactions.
                 */
                int a1 = molt->ilist[ftype].iatoms[1 + j + 0];
                int a2 = molt->ilist[ftype].iatoms[1 + j + 1];
                bexcl  = (bQMMM[a1] && bQMMM[a2]);
                /* A chemical bond between two QM atoms will be copied to
                 * the F_CONNBONDS list, for reasons mentioned above.
                 */
                if (bexcl && IS_CHEMBOND(ftype))
                {
                    InteractionList &ilist = molt->ilist[F_CONNBONDS];
                    ilist.iatoms.resize(ind_connbond + 3);
                    ilist.iatoms[ind_connbond++]  = ftype_connbond;
                    ilist.iatoms[ind_connbond++]  = a1;
                    ilist.iatoms[ind_connbond++]  = a2;
                }
            }
            else
            {
                /* MM interactions have to be excluded if they are included
                 * in the QM already. Because we use a link atom (H atom)
                 * when the QM/MM boundary runs through a chemical bond, this
                 * means that as long as one atom is MM, we still exclude,
                 * as the interaction is included in the QM via:
                 * QMatom1-QMatom2-QMatom-3-Linkatom.
                 */
                int numQmAtoms = 0;
                for (int jj = j + 1; jj < j + 1 + nratoms; jj++)
                {
                    if (bQMMM[molt->ilist[ftype].iatoms[jj]])
                    {
                        numQmAtoms++;
                    }
                }

                /* MiMiC treats link atoms as quantum atoms - therefore
                 * we do not need do additional exclusions here */
                if (qmmmMode == GmxQmmmMode::GMX_QMMM_MIMIC)
                {
                    bexcl = numQmAtoms == nratoms;
                }
                else
                {
                    bexcl = (numQmAtoms >= nratoms - 1);
                }

                if (bexcl && ftype == F_SETTLE)
                {
                    gmx_fatal(FARGS, "Can not apply QM to molecules with SETTLE, replace the moleculetype using QM and SETTLE by one without SETTLE");
                }
            }
            if (bexcl)
            {
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                InteractionList &ilist = molt->ilist[ftype];
                for (int k = j; k < ilist.size() - (nratoms + 1); k++)
                {
                    ilist.iatoms[k] = ilist.iatoms[k + (nratoms + 1)];
                }
                ilist.iatoms.resize(ilist.size() - (nratoms + 1));
            }
            else
            {
                j += nratoms+1; /* the +1 is for the functype */
            }
        }
    }
    /* Now, we search for atoms bonded to a QM atom because we also want
     * to exclude their nonbonded interactions with the QM atoms. The
     * reason for this is that this interaction is accounted for in the
     * linkatoms interaction with the QMatoms and would be counted
     * twice.  */

    if (qmmmMode != GmxQmmmMode::GMX_QMMM_MIMIC)
    {
        for (int i = 0; i < F_NRE; i++)
        {
            if (IS_CHEMBOND(i))
            {
                int j = 0;
                while (j < molt->ilist[i].size())
                {
                    int a1 = molt->ilist[i].iatoms[j + 1];
                    int a2 = molt->ilist[i].iatoms[j + 2];
                    if ((bQMMM[a1] && !bQMMM[a2]) || (!bQMMM[a1] && bQMMM[a2]))
                    {
                        if (link_nr >= link_max)
                        {
                            link_max += 10;
                            srenew(link_arr, link_max);
                        }
                        if (bQMMM[a1])
                        {
                            link_arr[link_nr++] = a2;
                        }
                        else
                        {
                            link_arr[link_nr++] = a1;
                        }
                    }
                    j += 3;
                }
            }
        }
    }
    snew(blink, molt->atoms.nr);
    for (int i = 0; i < molt->atoms.nr; i++)
    {
        blink[i] = FALSE;
    }

    if (qmmmMode != GmxQmmmMode::GMX_QMMM_MIMIC)
    {
        for (int i = 0; i < link_nr; i++)
        {
            blink[link_arr[i]] = TRUE;
        }
    }
    /* creating the exclusion block for the QM atoms. Each QM atom has
     * as excluded elements all the other QMatoms (and itself).
     */
    t_blocka  qmexcl;
    qmexcl.nr  = molt->atoms.nr;
    qmexcl.nra = qm_nr*(qm_nr+link_nr)+link_nr*qm_nr;
    snew(qmexcl.index, qmexcl.nr+1);
    snew(qmexcl.a, qmexcl.nra);
    int j = 0;
    for (int i = 0; i < qmexcl.nr; i++)
    {
        qmexcl.index[i] = j;
        if (bQMMM[i])
        {
            for (int k = 0; k < qm_nr; k++)
            {
                qmexcl.a[k+j] = qm_arr[k];
            }
            for (int k = 0; k < link_nr; k++)
            {
                qmexcl.a[qm_nr+k+j] = link_arr[k];
            }
            j += (qm_nr+link_nr);
        }
        if (blink[i])
        {
            for (int k = 0; k < qm_nr; k++)
            {
                qmexcl.a[k+j] = qm_arr[k];
            }
            j += qm_nr;
        }
    }
    qmexcl.index[qmexcl.nr] = j;

    /* and merging with the exclusions already present in sys.
     */

    gmx::ExclusionBlocks  qmexcl2;
    initExclusionBlocks(&qmexcl2, molt->atoms.nr);
    gmx::blockaToExclusionBlocks(&qmexcl, &qmexcl2);
    gmx::mergeExclusions(&(molt->excls), &qmexcl2);
    gmx::doneExclusionBlocks(&qmexcl2);

    /* Finally, we also need to get rid of the pair interactions of the
     * classical atom bonded to the boundary QM atoms with the QMatoms,
     * as this interaction is already accounted for by the QM, so also
     * here we run the risk of double counting! We proceed in a similar
     * way as we did above for the other bonded interactions: */
    for (int i = F_LJ14; i < F_COUL14; i++)
    {
        int nratoms = interaction_function[i].nratoms;
        int j       = 0;
        while (j < molt->ilist[i].size())
        {
            int  a1    = molt->ilist[i].iatoms[j+1];
            int  a2    = molt->ilist[i].iatoms[j+2];
            bool bexcl = ((bQMMM[a1] && bQMMM[a2]) ||
                          (blink[a1] && bQMMM[a2]) ||
                          (bQMMM[a1] && blink[a2]));
            if (bexcl)
            {
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                InteractionList &ilist = molt->ilist[i];
                for (int k = j; k < ilist.size() - (nratoms + 1); k++)
                {
                    ilist.iatoms[k] = ilist.iatoms[k + (nratoms + 1)];
                }
                ilist.iatoms.resize(ilist.size() - (nratoms + 1));
            }
            else
            {
                j += nratoms+1; /* the +1 is for the functype */
            }
        }
    }

    free(qm_arr);
    free(bQMMM);
    free(link_arr);
    free(blink);
} /* generate_qmexcl */

void generate_qmexcl(gmx_mtop_t *sys, t_inputrec *ir, warninp_t wi, GmxQmmmMode qmmmMode)
{
    /* This routine expects molt->molt[m].ilist to be of size F_NRE and ordered.
     */

    unsigned char   *grpnr;
    int              mol, nat_mol, nr_mol_with_qm_atoms = 0;
    gmx_molblock_t  *molb;
    bool             bQMMM;
    int              index_offset = 0;
    int              qm_nr        = 0;

    grpnr = sys->groups.grpnr[egcQMMM];

    for (size_t mb = 0; mb < sys->molblock.size(); mb++)
    {
        molb    = &sys->molblock[mb];
        nat_mol = sys->moltype[molb->type].atoms.nr;
        for (mol = 0; mol < molb->nmol; mol++)
        {
            bQMMM = FALSE;
            for (int i = 0; i < nat_mol; i++)
            {
                if ((grpnr ? grpnr[i] : 0) < (ir->opts.ngQM))
                {
                    bQMMM                    = TRUE;
                    qm_nr++;
                }
            }

            if (bQMMM)
            {
                nr_mol_with_qm_atoms++;
                if (molb->nmol > 1)
                {
                    /* We need to split this molblock */
                    if (mol > 0)
                    {
                        /* Split the molblock at this molecule */
                        auto pos = sys->molblock.begin() + mb + 1;
                        sys->molblock.insert(pos, sys->molblock[mb]);
                        sys->molblock[mb  ].nmol  = mol;
                        sys->molblock[mb+1].nmol -= mol;
                        mb++;
                        molb = &sys->molblock[mb];
                    }
                    if (molb->nmol > 1)
                    {
                        /* Split the molblock after this molecule */
                        auto pos = sys->molblock.begin() + mb + 1;
                        sys->molblock.insert(pos, sys->molblock[mb]);
                        molb = &sys->molblock[mb];
                        sys->molblock[mb  ].nmol  = 1;
                        sys->molblock[mb+1].nmol -= 1;
                    }

                    /* Create a copy of a moltype for a molecule
                     * containing QM atoms and append it in the end of the list
                     */
                    std::vector<gmx_moltype_t> temp(sys->moltype.size());
                    for (size_t i = 0; i < sys->moltype.size(); ++i)
                    {
                        copy_moltype(&sys->moltype[i], &temp[i]);
                    }
                    sys->moltype.resize(sys->moltype.size() + 1);
                    for (size_t i = 0; i < temp.size(); ++i)
                    {
                        copy_moltype(&temp[i], &sys->moltype[i]);
                    }
                    copy_moltype(&sys->moltype[molb->type], &sys->moltype.back());
                    /* Copy the exclusions to a new array, since this is the only
                     * thing that needs to be modified for QMMM.
                     */
                    copy_blocka(&sys->moltype[molb->type].excls,
                                &sys->moltype.back().excls);
                    /* Set the molecule type for the QMMM molblock */
                    molb->type = sys->moltype.size() - 1;
                }
                generate_qmexcl_moltype(&sys->moltype[molb->type], grpnr, ir, qmmmMode);
            }
            if (grpnr)
            {
                grpnr += nat_mol;
            }
            index_offset += nat_mol;
        }
    }
    if (qmmmMode == GmxQmmmMode::GMX_QMMM_ORIGINAL &&
        nr_mol_with_qm_atoms > 1)
    {
        /* generate a warning is there are QM atoms in different
         * topologies. In this case it is not possible at this stage to
         * mutualy exclude the non-bonded interactions via the
         * exclusions (AFAIK). Instead, the user is advised to use the
         * energy group exclusions in the mdp file
         */
        warning_note(wi,
                     "\nThe QM subsystem is divided over multiple topologies. "
                     "The mutual non-bonded interactions cannot be excluded. "
                     "There are two ways to achieve this:\n\n"
                     "1) merge the topologies, such that the atoms of the QM "
                     "subsystem are all present in one single topology file. "
                     "In this case this warning will dissappear\n\n"
                     "2) exclude the non-bonded interactions explicitly via the "
                     "energygrp-excl option in the mdp file. if this is the case "
                     "this warning may be ignored"
                     "\n\n");
    }
}
