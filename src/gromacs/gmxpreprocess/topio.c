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
#include "gmxpre.h"

#include "topio.h"

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxpreprocess/gmxcpp.h"
#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/topshake.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/vsite_parm.h"
#include "gromacs/legacyheaders/genborn.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/warninp.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#define OPENDIR     '[' /* starting sign for directive */
#define CLOSEDIR    ']' /* ending sign for directive   */

static void free_nbparam(t_nbparam **param, int nr)
{
    int i;

    for (i = 0; i < nr; i++)
    {
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
    nnn       = sqrt(ntp);
    nrfp      = NRFP(F_LJ);
    nrfpA     = interaction_function[F_LJ14].nrfpA;
    nrfpB     = interaction_function[F_LJ14].nrfpB;
    pairs->nr = ntp;

    if ((nrfp  != nrfpA) || (nrfpA != nrfpB))
    {
        gmx_incons("Number of force parameters in gen_pairs wrong");
    }

    fprintf(stderr, "Generating 1-4 interactions: fudge = %g\n", fudge);
    if (debug)
    {
        fprintf(debug, "Fudge factor for 1-4 interactions: %g\n", fudge);
        fprintf(debug, "Holy Cow! there are %d types\n", ntp);
    }
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
            pairs->param[i].c[nrfp+j] = scaling*nbs->param[i].c[j];
        }
    }
}

double check_mol(gmx_mtop_t *mtop, warninp_t wi)
{
    char     buf[256];
    int      i, mb, nmol, ri, pt;
    double   q;
    real     m, mB;
    t_atoms *atoms;

    /* Check mass and charge */
    q = 0.0;

    for (mb = 0; mb < mtop->nmoltype; mb++)
    {
        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
        nmol  = mtop->molblock[mb].nmol;
        for (i = 0; (i < atoms->nr); i++)
        {
            q += nmol*atoms->atom[i].q;
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

static void sum_q(t_atoms *atoms, int n, double *qt, double *qBt)
{
    double  qmolA, qmolB;
    int     i;

    /* sum charge */
    qmolA = 0;
    qmolB = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        qmolA += atoms->atom[i].q;
        qmolB += atoms->atom[i].qB;
    }
    /* Unfortunately an absolute comparison,
     * but this avoids unnecessary warnings and gmx-users mails.
     */
    if (fabs(qmolA) >= 1e-6 || fabs(qmolB) >= 1e-6)
    {
        *qt  += n*qmolA;
        *qBt += n*qmolB;
    }
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
        *nb = strtol(nb_str, NULL, 10);
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
        *comb = strtol(comb_str, NULL, 10);
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
    char      **cppopts   = NULL;
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
    cppopts[ncppopts-1] = NULL;

    return cppopts;
}


int
find_gb_bondlength(t_params *plist, int ai, int aj, real *length)
{
    int i, j, a1, a2;

    int found = 0;
    int status;

    for (i = 0; i < F_NRE && !found; i++)
    {
        if (IS_CHEMBOND(i))
        {
            for (j = 0; j < plist[i].nr; j++)
            {
                a1 = plist[i].param[j].a[0];
                a2 = plist[i].param[j].a[1];

                if ( (a1 == ai && a2 == aj) || (a1 == aj && a2 == ai))
                {
                    /* Equilibrium bond distance */
                    *length = plist[i].param[j].c[0];
                    found   = 1;
                }
            }
        }
    }
    status = !found;

    return status;
}


int
find_gb_anglelength(t_params *plist, int ai, int ak, real *length)
{
    int  i, j, a1, a2, a3;
    real r12, r23, a123;
    int  found = 0;
    int  status, status1, status2;

    r12 = r23 = 0;

    for (i = 0; i < F_NRE && !found; i++)
    {
        if (IS_ANGLE(i))
        {
            for (j = 0; j < plist[i].nr; j++)
            {
                a1 = plist[i].param[j].a[0];
                a2 = plist[i].param[j].a[1];
                a3 = plist[i].param[j].a[2];

                /* We dont care what the middle atom is, but use it below */
                if ( (a1 == ai && a3 == ak) || (a1 == ak && a3 == ai) )
                {
                    /* Equilibrium bond distance */
                    a123 = plist[i].param[j].c[0];
                    /* Use middle atom to find reference distances r12 and r23 */
                    status1 = find_gb_bondlength(plist, a1, a2, &r12);
                    status2 = find_gb_bondlength(plist, a2, a3, &r23);

                    if (status1 == 0 && status2 == 0)
                    {
                        /* cosine theorem to get r13 */
                        *length = sqrt(r12*r12+r23*r23-(2*r12*r23*cos(a123/RAD2DEG)));
                        found   = 1;
                    }
                }
            }
        }
    }
    status = !found;

    return status;
}

int
generate_gb_exclusion_interactions(t_molinfo *mi, gpp_atomtype_t atype, t_nextnb *nnb)
{
    int          i, j, k, n, ai, aj, ti, tj;
    int          n12, n13, n14;
    int          ftype;
    t_param      param;
    t_params *   plist;
    t_atoms *    at;
    real         radiusi, radiusj;
    real         gb_radiusi, gb_radiusj;
    real         param_c2, param_c4;
    real         distance;

    plist = mi->plist;
    at    = &mi->atoms;

    for (n = 1; n <= nnb->nrex; n++)
    {
        switch (n)
        {
            case 1:
                ftype    = F_GB12;
                param_c2 = STILL_P2;
                param_c4 = 0.8875;
                break;
            case 2:
                ftype    = F_GB13;
                param_c2 = STILL_P3;
                param_c4 = 0.3516;
                break;
            default:
                /* Put all higher-order exclusions into 1,4 list so we dont miss them */
                ftype    = F_GB14;
                param_c2 = STILL_P3;
                param_c4 = 0.3516;
                break;
        }

        for (ai = 0; ai < nnb->nr; ai++)
        {
            ti         = at->atom[ai].type;
            radiusi    = get_atomtype_radius(ti, atype);
            gb_radiusi = get_atomtype_gb_radius(ti, atype);

            for (j = 0; j < nnb->nrexcl[ai][n]; j++)
            {
                aj = nnb->a[ai][n][j];

                /* Only add the interactions once */
                if (aj > ai)
                {
                    tj         = at->atom[aj].type;
                    radiusj    = get_atomtype_radius(tj, atype);
                    gb_radiusj = get_atomtype_gb_radius(tj, atype);

                    /* There is an exclusion of type "ftype" between atoms ai and aj */
                    param.a[0] = ai;
                    param.a[1] = aj;

                    /* Reference distance, not used for 1-4 interactions */
                    switch (ftype)
                    {
                        case F_GB12:
                            if (find_gb_bondlength(plist, ai, aj, &distance) != 0)
                            {
                                gmx_fatal(FARGS, "Cannot find bond length for atoms %d-%d", ai, aj);
                            }
                            break;
                        case F_GB13:
                            if (find_gb_anglelength(plist, ai, aj, &distance) != 0)
                            {
                                gmx_fatal(FARGS, "Cannot find length for atoms %d-%d involved in angle", ai, aj);
                            }
                            break;
                        default:
                            distance = -1;
                            break;
                    }
                    /* Assign GB parameters */
                    /* Sum of radii */
                    param.c[0] = radiusi+radiusj;
                    /* Reference distance distance */
                    param.c[1] = distance;
                    /* Still parameter */
                    param.c[2] = param_c2;
                    /* GB radius */
                    param.c[3] = gb_radiusi+gb_radiusj;
                    /* Parameter */
                    param.c[4] = param_c4;

                    /* Add it to the parameter list */
                    add_param_to_list(&plist[ftype], &param);
                }
            }
        }
    }
    return 0;
}


static void make_atoms_sys(int nmolb, const gmx_molblock_t *molb,
                           const t_molinfo *molinfo,
                           t_atoms *atoms)
{
    int            mb, m, a;
    const t_atoms *mol_atoms;

    atoms->nr   = 0;
    atoms->atom = NULL;

    for (mb = 0; mb < nmolb; mb++)
    {
        mol_atoms = &molinfo[molb[mb].type].atoms;

        srenew(atoms->atom, atoms->nr + molb[mb].nmol*mol_atoms->nr);

        for (m = 0; m < molb[mb].nmol; m++)
        {
            for (a = 0; a < mol_atoms->nr; a++)
            {
                atoms->atom[atoms->nr++] = mol_atoms->atom[a];
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
                         int         *nmolblock,
                         gmx_molblock_t **molblock,
                         gmx_bool        bFEP,
                         gmx_bool        bGenborn,
                         gmx_bool        bZero,
                         warninp_t   wi)
{
    FILE           *out;
    int             i, sl, nb_funct, comb;
    char           *pline = NULL, **title = NULL;
    char            line[STRLEN], errbuf[256], comb_str[256], nb_str[256];
    char            genpairs[32];
    char           *dirstr, *dummy2;
    int             nrcopies, nmol, nmolb = 0, nscan, ncombs, ncopy;
    double          fLJ, fQQ, fPOW;
    gmx_molblock_t *molb  = NULL;
    t_topology     *block = NULL;
    t_molinfo      *mi0   = NULL;
    DirStack       *DS;
    directive       d, newd;
    t_nbparam     **nbparam, **pair;
    gmx_bool        bIntermolecularInteractions;
    t_block2       *block2;
    real            fudgeLJ = -1;    /* Multiplication factor to generate 1-4 from LJ */
    gmx_bool        bReadDefaults, bReadMolType, bGenPairs, bWarn_copy_A_B;
    double          qt = 0, qBt = 0; /* total charge */
    t_bond_atomtype batype;
    int             lastcg = -1;
    int             dcatt  = -1, nmol_couple;
    /* File handling variables */
    int             status, done;
    gmx_cpp_t       handle;
    char           *tmp_line = NULL;
    char            warn_buf[STRLEN];
    const char     *floating_point_arithmetic_tip =
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
        out = NULL;
    }

    /* open input file */
    status = cpp_open_file(infile, &handle, cpp_opts(define, include, wi));
    if (status != 0)
    {
        gmx_fatal(FARGS, cpp_error(&handle, status));
    }

    /* some local variables */
    DS_Init(&DS);         /* directive stack			 */
    nmol     = 0;         /* no molecules yet...			 */
    d        = d_invalid; /* first thing should be a directive   */
    nbparam  = NULL;      /* The temporary non-bonded matrix       */
    pair     = NULL;      /* The temporary pair interaction matrix */
    block2   = NULL;      /* the extra exclusions			 */
    nb_funct = F_LJ;

    *reppow  = 12.0;      /* Default value for repulsion power     */

    *intermolecular_interactions = NULL;

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
                gmx_fatal(FARGS, cpp_error(&handle, status));
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
                        gmx_fatal(FARGS, cpp_error(&handle, status));
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
            if ((int)strlen(pline) > 0)
            {
                if (pline[0] == OPENDIR)
                {
                    /* A directive on this line: copy the directive
                     * without the brackets into dirstr, then
                     * skip spaces and tabs on either side of directive
                     */
                    dirstr = gmx_strdup((pline+1));
                    if ((dummy2 = strchr (dirstr, CLOSEDIR)) != NULL)
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
                        if (debug)
                        {
                            fprintf(debug, "found directive '%s'\n", dir2str(newd));
                        }
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
                            if (*intermolecular_interactions == NULL)
                            {
                                /* We (mis)use the moleculetype processing
                                 * to process the intermolecular interactions
                                 * by making a "molecule" of the size of the system.
                                 */
                                snew(*intermolecular_interactions, 1);
                                init_molinfo(*intermolecular_interactions);
                                mi0 = *intermolecular_interactions;
                                make_atoms_sys(nmolb, molb, *molinfo,
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
                                    &nbparam, bGenPairs ? &pair : NULL, wi);
                            break;

                        case d_bondtypes:
                            push_bt(d, plist, 2, NULL, batype, pline, wi);
                            break;
                        case d_constrainttypes:
                            push_bt(d, plist, 2, NULL, batype, pline, wi);
                            break;
                        case d_pairtypes:
                            if (bGenPairs)
                            {
                                push_nbt(d, pair, atype, pline, F_LJ14, wi);
                            }
                            else
                            {
                                push_bt(d, plist, 2, atype, NULL, pline, wi);
                            }
                            break;
                        case d_angletypes:
                            push_bt(d, plist, 3, NULL, batype, pline, wi);
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
                            push_gb_params(atype, pline, wi);
                            break;

                        case d_implicit_surface_params:
                            gmx_fatal(FARGS, "Implicit surface directive not supported yet.");
                            break;

                        case d_cmaptypes:
                            push_cmaptype(d, plist, 5, atype, batype, pline, wi);
                            break;

                        case d_moleculetype:
                        {
                            if (!bReadMolType)
                            {
                                int ntype;
                                if (opts->couple_moltype != NULL &&
                                    (opts->couple_lam0 == ecouplamNONE ||
                                     opts->couple_lam0 == ecouplamQ ||
                                     opts->couple_lam1 == ecouplamNONE ||
                                     opts->couple_lam1 == ecouplamQ))
                                {
                                    dcatt = add_atomtype_decoupled(symtab, atype,
                                                                   &nbparam, bGenPairs ? &pair : NULL);
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
                            srenew(block2, nmol);
                            block2[nmol-1].nr = 0;
                            mi0               = &((*molinfo)[nmol-1]);
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
                            assert(block2);
                            if (!block2[nmol-1].nr)
                            {
                                init_block2(&(block2[nmol-1]), mi0->atoms.nr);
                            }
                            push_excl(pline, &(block2[nmol-1]));
                            break;
                        case d_system:
                            trim(pline);
                            title = put_symtab(symtab, pline);
                            break;
                        case d_molecules:
                        {
                            int      whichmol;
                            gmx_bool bCouple;

                            push_mol(nmol, *molinfo, pline, &whichmol, &nrcopies, wi);
                            mi0 = &((*molinfo)[whichmol]);
                            srenew(molb, nmolb+1);
                            molb[nmolb].type = whichmol;
                            molb[nmolb].nmol = nrcopies;
                            nmolb++;

                            bCouple = (opts->couple_moltype != NULL &&
                                       (gmx_strcasecmp("system", opts->couple_moltype) == 0 ||
                                        gmx_strcasecmp(*(mi0->name), opts->couple_moltype) == 0));
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
                                merge_excl(&(mi0->excls), &(block2[whichmol]));
                                done_block2(&(block2[whichmol]));
                                make_shake(mi0->plist, &mi0->atoms, opts->nshake);



                                /* nnb contains information about first,2nd,3rd bonded neighbors.
                                 * Use this to generate GB 1-2,1-3,1-4 interactions when necessary.
                                 */
                                if (bGenborn == TRUE)
                                {
                                    generate_gb_exclusion_interactions(mi0, atype, &nnb);
                                }

                                done_nnb(&nnb);

                                if (bCouple)
                                {
                                    convert_moltype_couple(mi0, dcatt, *fudgeQQ,
                                                           opts->couple_lam0, opts->couple_lam1,
                                                           opts->bCoupleIntra,
                                                           nb_funct, &(plist[nb_funct]));
                                }
                                stupid_fill_block(&mi0->mols, mi0->atoms.nr, TRUE);
                                mi0->bProcessed = TRUE;
                            }
                            break;
                        }
                        default:
                            fprintf (stderr, "case: %d\n", (int)d);
                            gmx_incons("unknown directive");
                    }
                }
            }
            sfree(pline);
            pline = NULL;
        }
    }
    while (!done);
    status = cpp_close_file(&handle);
    if (status != eCPP_OK)
    {
        gmx_fatal(FARGS, cpp_error(&handle, status));
    }
    cpp_done();
    if (out)
    {
        gmx_fio_fclose(out);
    }

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
    DS_Done (&DS);
    for (i = 0; i < nmol; i++)
    {
        done_block2(&(block2[i]));
    }
    free(block2);

    done_bond_atomtype(&batype);

    if (*intermolecular_interactions != NULL)
    {
        sfree(mi0->atoms.atom);
    }

    *nrmols = nmol;

    *nmolblock = nmolb;
    *molblock  = molb;

    return title;
}

char **do_top(gmx_bool          bVerbose,
              const char       *topfile,
              const char       *topppfile,
              t_gromppopts     *opts,
              gmx_bool          bZero,
              t_symtab         *symtab,
              t_params          plist[],
              int              *combination_rule,
              double           *repulsion_power,
              real             *fudgeQQ,
              gpp_atomtype_t    atype,
              int              *nrmols,
              t_molinfo       **molinfo,
              t_molinfo       **intermolecular_interactions,
              t_inputrec       *ir,
              int              *nmolblock,
              gmx_molblock_t  **molblock,
              gmx_bool          bGenborn,
              warninp_t         wi)
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
        tmpfile = NULL;
    }

    if (bVerbose)
    {
        printf("processing topology...\n");
    }
    title = read_topol(topfile, tmpfile, opts->define, opts->include,
                       symtab, atype,
                       nrmols, molinfo, intermolecular_interactions,
                       plist, combination_rule, repulsion_power,
                       opts, fudgeQQ, nmolblock, molblock,
                       ir->efep != efepNO, bGenborn, bZero, wi);
    if ((*combination_rule != eCOMB_GEOMETRIC) &&
        (ir->vdwtype == evdwUSER))
    {
        warning(wi, "Using sigma/epsilon based combination rules with"
                " user supplied potential function may produce unwanted"
                " results");
    }

    return title;
}


static void generate_qmexcl_moltype(gmx_moltype_t *molt, unsigned char *grpnr,
                                    t_inputrec *ir)
{
    /* This routine expects molt->ilist to be of size F_NRE and ordered. */

    /* generates the exclusions between the individual QM atoms, as
     * these interactions should be handled by the QM subroutines and
     * not by the gromacs routines
     */
    int
        i, j, l, k = 0, jmax, qm_max = 0, qm_nr = 0, nratoms = 0, link_nr = 0, link_max = 0;
    atom_id
       *qm_arr = NULL, *link_arr = NULL, a1, a2, a3, a4, ftype = 0;
    t_blocka
        qmexcl;
    t_block2
        qmexcl2;
    gmx_bool
       *bQMMM, *blink, bexcl;

    /* First we search and select the QM atoms in an qm_arr array that
     * we use to create the exclusions.
     *
     * we take the possibility into account that a user has defined more
     * than one QM group:
     *
     * for that we also need to do this an ugly work-about just in case
     * the QM group contains the entire system...
     */
    jmax = ir->opts.ngQM;

    /* we first search for all the QM atoms and put them in an array
     */
    for (j = 0; j < jmax; j++)
    {
        for (i = 0; i < molt->atoms.nr; i++)
        {
            if (qm_nr >= qm_max)
            {
                qm_max += 100;
                srenew(qm_arr, qm_max);
            }
            if ((grpnr ? grpnr[i] : 0) == j)
            {
                qm_arr[qm_nr++] = i;
            }
        }
    }
    /* bQMMM[..] is an array containin TRUE/FALSE for atoms that are
     * QM/not QM. We first set all elements to false. Afterwards we use
     * the qm_arr to change the elements corresponding to the QM atoms
     * to TRUE.
     */
    snew(bQMMM, molt->atoms.nr);
    for (i = 0; i < molt->atoms.nr; i++)
    {
        bQMMM[i] = FALSE;
    }
    for (i = 0; i < qm_nr; i++)
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

    /* first check weter we already have CONNBONDS: */
    if (molt->ilist[F_CONNBONDS].nr != 0)
    {
        fprintf(stderr, "nr. of CONNBONDS present already: %d\n",
                molt->ilist[F_CONNBONDS].nr/3);
        ftype = molt->ilist[F_CONNBONDS].iatoms[0];
        k     = molt->ilist[F_CONNBONDS].nr;
    }
    /* now we delete all bonded interactions, except the ones describing
     * a chemical bond. These are converted to CONNBONDS
     */
    for (i = 0; i < F_LJ; i++)
    {
        if (i == F_CONNBONDS)
        {
            continue;
        }
        nratoms = interaction_function[i].nratoms;
        j       = 0;
        while (j < molt->ilist[i].nr)
        {
            bexcl = FALSE;
            switch (nratoms)
            {
                case 2:
                    a1    = molt->ilist[i].iatoms[j+1];
                    a2    = molt->ilist[i].iatoms[j+2];
                    bexcl = (bQMMM[a1] && bQMMM[a2]);
                    /* a bonded beteen two QM atoms will be copied to the
                     * CONNBONDS list, for reasons mentioned above
                     */
                    if (bexcl && i < F_ANGLES)
                    {
                        srenew(molt->ilist[F_CONNBONDS].iatoms, k+3);
                        molt->ilist[F_CONNBONDS].nr         += 3;
                        molt->ilist[F_CONNBONDS].iatoms[k++] = ftype;
                        molt->ilist[F_CONNBONDS].iatoms[k++] = a1;
                        molt->ilist[F_CONNBONDS].iatoms[k++] = a2;
                    }
                    break;
                case 3:
                    a1    = molt->ilist[i].iatoms[j+1];
                    a2    = molt->ilist[i].iatoms[j+2];
                    a3    = molt->ilist[i].iatoms[j+3];
                    bexcl = ((bQMMM[a1] && bQMMM[a2]) ||
                             (bQMMM[a1] && bQMMM[a3]) ||
                             (bQMMM[a2] && bQMMM[a3]));
                    break;
                case 4:
                    a1    = molt->ilist[i].iatoms[j+1];
                    a2    = molt->ilist[i].iatoms[j+2];
                    a3    = molt->ilist[i].iatoms[j+3];
                    a4    = molt->ilist[i].iatoms[j+4];
                    bexcl = ((bQMMM[a1] && bQMMM[a2] && bQMMM[a3]) ||
                             (bQMMM[a1] && bQMMM[a2] && bQMMM[a4]) ||
                             (bQMMM[a1] && bQMMM[a3] && bQMMM[a4]) ||
                             (bQMMM[a2] && bQMMM[a3] && bQMMM[a4]));
                    break;
                default:
                    gmx_fatal(FARGS, "no such bonded interactions with %d atoms\n", nratoms);
            }
            if (bexcl)
            {
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                molt->ilist[i].nr -= (nratoms+1);
                for (l = j; l < molt->ilist[i].nr; l++)
                {
                    molt->ilist[i].iatoms[l] = molt->ilist[i].iatoms[l+(nratoms+1)];
                }
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

    for (i = 0; i < F_NRE; i++)
    {
        if (IS_CHEMBOND(i))
        {
            j = 0;
            while (j < molt->ilist[i].nr)
            {
                a1 = molt->ilist[i].iatoms[j+1];
                a2 = molt->ilist[i].iatoms[j+2];
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
    snew(blink, molt->atoms.nr);
    for (i = 0; i < molt->atoms.nr; i++)
    {
        blink[i] = FALSE;
    }
    for (i = 0; i < link_nr; i++)
    {
        blink[link_arr[i]] = TRUE;
    }
    /* creating the exclusion block for the QM atoms. Each QM atom has
     * as excluded elements all the other QMatoms (and itself).
     */
    qmexcl.nr  = molt->atoms.nr;
    qmexcl.nra = qm_nr*(qm_nr+link_nr)+link_nr*qm_nr;
    snew(qmexcl.index, qmexcl.nr+1);
    snew(qmexcl.a, qmexcl.nra);
    j = 0;
    for (i = 0; i < qmexcl.nr; i++)
    {
        qmexcl.index[i] = j;
        if (bQMMM[i])
        {
            for (k = 0; k < qm_nr; k++)
            {
                qmexcl.a[k+j] = qm_arr[k];
            }
            for (k = 0; k < link_nr; k++)
            {
                qmexcl.a[qm_nr+k+j] = link_arr[k];
            }
            j += (qm_nr+link_nr);
        }
        if (blink[i])
        {
            for (k = 0; k < qm_nr; k++)
            {
                qmexcl.a[k+j] = qm_arr[k];
            }
            j += qm_nr;
        }
    }
    qmexcl.index[qmexcl.nr] = j;

    /* and merging with the exclusions already present in sys.
     */

    init_block2(&qmexcl2, molt->atoms.nr);
    b_to_b2(&qmexcl, &qmexcl2);
    merge_excl(&(molt->excls), &qmexcl2);
    done_block2(&qmexcl2);

    /* Finally, we also need to get rid of the pair interactions of the
     * classical atom bonded to the boundary QM atoms with the QMatoms,
     * as this interaction is already accounted for by the QM, so also
     * here we run the risk of double counting! We proceed in a similar
     * way as we did above for the other bonded interactions: */
    for (i = F_LJ14; i < F_COUL14; i++)
    {
        nratoms = interaction_function[i].nratoms;
        j       = 0;
        while (j < molt->ilist[i].nr)
        {
            a1    = molt->ilist[i].iatoms[j+1];
            a2    = molt->ilist[i].iatoms[j+2];
            bexcl = ((bQMMM[a1] && bQMMM[a2]) ||
                     (blink[a1] && bQMMM[a2]) ||
                     (bQMMM[a1] && blink[a2]));
            if (bexcl)
            {
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                molt->ilist[i].nr -= (nratoms+1);
                for (k = j; k < molt->ilist[i].nr; k++)
                {
                    molt->ilist[i].iatoms[k] = molt->ilist[i].iatoms[k+(nratoms+1)];
                }
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

void generate_qmexcl(gmx_mtop_t *sys, t_inputrec *ir, warninp_t    wi)
{
    /* This routine expects molt->molt[m].ilist to be of size F_NRE and ordered.
     */

    unsigned char  *grpnr;
    int             mb, mol, nat_mol, i, nr_mol_with_qm_atoms = 0;
    gmx_molblock_t *molb;
    gmx_bool        bQMMM;

    grpnr = sys->groups.grpnr[egcQMMM];

    for (mb = 0; mb < sys->nmolblock; mb++)
    {
        molb    = &sys->molblock[mb];
        nat_mol = sys->moltype[molb->type].atoms.nr;
        for (mol = 0; mol < molb->nmol; mol++)
        {
            bQMMM = FALSE;
            for (i = 0; i < nat_mol; i++)
            {
                if ((grpnr ? grpnr[i] : 0) < ir->opts.ngQM)
                {
                    bQMMM = TRUE;
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
                        sys->nmolblock++;
                        srenew(sys->molblock, sys->nmolblock);
                        for (i = sys->nmolblock-2; i >= mb; i--)
                        {
                            sys->molblock[i+1] = sys->molblock[i];
                        }
                        sys->molblock[mb  ].nmol  = mol;
                        sys->molblock[mb+1].nmol -= mol;
                        mb++;
                        molb = &sys->molblock[mb];
                    }
                    if (molb->nmol > 1)
                    {
                        /* Split the molblock after this molecule */
                        sys->nmolblock++;
                        srenew(sys->molblock, sys->nmolblock);
                        molb = &sys->molblock[mb];
                        for (i = sys->nmolblock-2; i >= mb; i--)
                        {
                            sys->molblock[i+1] = sys->molblock[i];
                        }
                        sys->molblock[mb  ].nmol  = 1;
                        sys->molblock[mb+1].nmol -= 1;
                    }

                    /* Add a moltype for the QMMM molecule */
                    sys->nmoltype++;
                    srenew(sys->moltype, sys->nmoltype);
                    /* Copy the moltype struct */
                    sys->moltype[sys->nmoltype-1] = sys->moltype[molb->type];
                    /* Copy the exclusions to a new array, since this is the only
                     * thing that needs to be modified for QMMM.
                     */
                    copy_blocka(&sys->moltype[molb->type     ].excls,
                                &sys->moltype[sys->nmoltype-1].excls);
                    /* Set the molecule type for the QMMM molblock */
                    molb->type = sys->nmoltype - 1;
                }
                generate_qmexcl_moltype(&sys->moltype[molb->type], grpnr, ir);
            }
            if (grpnr)
            {
                grpnr += nat_mol;
            }
        }
    }
    if (nr_mol_with_qm_atoms > 1)
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
