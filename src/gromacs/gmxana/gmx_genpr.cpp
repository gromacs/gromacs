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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "config.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/waxsdebye/scattering_factors.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/math/vec.h"
#include "gmx_ana.h"

static void gen_macros(FILE *out, gmx::ScatteringFactorTable &sft)
{
    fprintf(out, "\n; Macros to translate the types used in the pair section to something interpretable for poor old grompp.\n");

    for (std::vector<gmx::ScatteringFactorPointer>::iterator sfp = sft.beginSFP(); (sfp < sft.endSFP()); ++sfp)
    {
        fprintf(out, "#define SF_%s_%s %d\n",
                (*sfp)->residue().c_str(),
                (*sfp)->atom().c_str(),
                (*sfp)->type());
    }
    fprintf(out, "\n");
}

static void com()
{
    //
}

static void gen_residue_matrix(const char *fn,
                               gmx::ScatteringFactorTable &sft,
                               t_atoms *myatoms, rvec **x, int resmax)
{
    FILE     *out;
    int       vsmin;
    int       nres  = 0;
    int       i, j, resind, rn2, vsnum, atomcount;
    rvec      cog;
    t_symtab *symtab;

    out = gmx_fio_fopen(fn, "w");

    /* warn about running without the -resmax option */
    if (resmax <= 0)
    {
        printf("\nNOTE: Generating Debye scattering terms as bonded interactions without a specified maximum number of residues (option -resmax). If the input file contains multiple molecules, this is not meaningful.\n");
    }

    nres = myatoms->nres;
    /* special treatment if a maximum number of residues is specified: */
    if ((resmax > 0) && (resmax <= nres))
    {
        /* find vsmin, the atom number of the first virtual site: */
        nres = resmax;
        for (i = 0; (i < myatoms->nr); i++)
        {
            resind = myatoms->atom[i].resind;
            if (resind >= nres)
            {
                break;
            }
        }
        vsmin = i;
    }
    else
    {
        resmax = nres;
        vsmin = myatoms->nr;
    }

    /* extend the topology to make space for the new MW atoms: */
    snew(symtab, 1);
    open_symtab(symtab);
    srenew((*x), myatoms->nr+nres);
    add_t_atoms(myatoms, nres, nres);

    /* shift all atoms which are not included in the matrix upwards: */
    for (i = myatoms->nr-1; (i >= vsmin); i--)
    {
        copy_rvec((*x)[i-nres], (*x)[i]);
        memcpy(&myatoms->atom[i], &myatoms->atom[i-nres], sizeof(myatoms->atom[i]));
        myatoms->atom[i].resind += nres;
        myatoms->atomname[i]     = myatoms->atomname[i-nres];
    }
    /* shift all residues in the same way: */
    for (i = myatoms->nres-1; (i >= 2*resmax); i--)
    {
        memcpy(&myatoms->resinfo[i], &myatoms->resinfo[i-nres], sizeof(myatoms->resinfo[i]));
    }

    /* Calculate the positions of the new MW atoms */
    vsnum = vsmin;
    for (i = 0; (i < myatoms->nr); )
    {
        atomcount = 0;
        clear_rvec(cog);
        resind = myatoms->atom[i].resind;
        if (resind >= nres)
        {
            break;
        }
        for (; (i < myatoms->nr) && (myatoms->atom[i].resind == resind); i++)
        {
            for (j = 0; (j<3); j++)
            {
                cog[j] += (*x)[i][j];
            }
            atomcount++;
        }
        for (j = 0; (j<3); j++)
        {
            cog[j] /= atomcount;
        }
        copy_rvec(cog, (*x)[vsnum++]);
    }

    /* construct and write out the topology section */
    fprintf(out, "; Include this file straight after the [ atoms ] section\n");
    fprintf(out, "; since it adds additional atoms to the topology.\n\n");
    for (i = 0; (i < nres); i++)
    {
        resind = nres+i+1;
        fprintf(out, "    %d    MW     %d   %s    MW     %d  0\n",
                vsmin+i+1, resind, *(myatoms->resinfo[i].name),
                vsmin+i+1);
        myatoms->atom[vsmin+i].m          = 0;
        myatoms->atom[vsmin+i].q          = 0;
        myatoms->atom[vsmin+i].mB         = 0;
        myatoms->atom[vsmin+i].qB         = 0;
        myatoms->atom[vsmin+i].type       = 0;
        myatoms->atom[vsmin+i].typeB      = 0;
        myatoms->atom[vsmin+i].ptype      = eptVSite;
        myatoms->atom[vsmin+i].resind     = resind-1;
        myatoms->atom[vsmin+i].atomnumber = 0;
        strcpy(myatoms->atom[vsmin+i].elem, "0");
        myatoms->atomname[vsmin+i] = put_symtab(symtab, "MW");
        t_atoms_set_resinfo(myatoms, vsmin+i, symtab, "MW", resind, ' ', 0, ' ');
    }
    fprintf(out, "\n");
    fprintf(out, "[ virtual_sitesn ]\n");
    vsnum = vsmin+1;
    printf("NOTE: Using center-of-geometry (COG) vsites\n");
    for (i = 0; (i < myatoms->nr); )
    {
        resind = myatoms->atom[i].resind;
        if (resind >= nres)
        {
            break;
        }
        /* Using center-of-geometry (COG) vsites: */
        fprintf(out, "%5d  1", vsnum++);
        for (; (i < myatoms->nr) && (myatoms->atom[i].resind == resind); i++)
        {
            fprintf(out, "  %d", i+1);
        }
        fprintf(out, "\n");
    }

    gen_macros(out, sft);

    fprintf(out, "\n[ pairs ]\n");
    for (resind = 1; (resind <= nres); resind++)
    {
        for (rn2 = resind; (rn2 <= nres); rn2++)
        {
            fprintf(out, "%5d  %5d  3  SF_%s_MW  SF_%s_MW\n",
                    vsmin+resind, vsmin+rn2,
                    *myatoms->resinfo[resind-1].name,
                    *myatoms->resinfo[rn2-1].name);
        }
    }
    gmx_fio_fclose(out);
}

static void gen_atom_matrix(const char                 *fn,
                            gmx::ScatteringFactorTable &sft,
                            t_atoms                    *myatoms)
{
    FILE                    *out;
    int                      ai, aj;
    std::vector<std::string> res;

    out = gmx_fio_fopen(fn, "w");
    if (NULL == out)
    {
        return;
    }

    fprintf(out, "; Include this anywhere after the atoms section in your topology file\n");

    gen_macros(out, sft);

    gmx_atomprop_t aps = gmx_atomprop_init();
    for (ai = 0; (ai < myatoms->nr); ai++)
    {
        int   ri = myatoms->atom[ai].resind;
        char  ibuf[32];
        char *element;

        ibuf[0] = '\0';

        // First try full residue- and atom names
        if (sft.computeScatteringFactor((*myatoms->resinfo[ri].name),
                                        (*myatoms->atomname[ai]),
                                        1.0) != -1)
        {
            sprintf(ibuf, "SF_%s_%s", (*myatoms->resinfo[ri].name),
                    (*myatoms->atomname[ai]));
        }
        // Then try full residue- and element name if present
        else
        {
            element = NULL;
            if ((0 != myatoms->atom[ai].atomnumber) &&
                (strlen(myatoms->atom[ai].elem) > 0))
            {
                element = myatoms->atom[ai].elem;
            }
            else
            {
                real value;
                if (gmx_atomprop_query(aps, epropElement,
                                       (*myatoms->resinfo[ri].name),
                                       (*myatoms->atomname[ai]),
                                       &value))
                {
                    int atomnumber = (int)floor(value+0.001);
                    element = gmx_atomprop_element(aps, atomnumber);
                }
            }

            if (NULL != element)
            {
                if (sft.computeScatteringFactor((*myatoms->resinfo[ri].name),
                                                element, 1.0) != -1)
                {
                    sprintf(ibuf, "SF_%s_%s", (*myatoms->resinfo[ri].name),
                            element);
                }
                else if (sft.computeScatteringFactor("*", element, 1.0) != -1)
                {
                    sprintf(ibuf, "SF_*_%s", element);
                }
            }
        }

        if (0 == strlen(ibuf))
        {
            char errbuf[STRLEN];
            sprintf(errbuf, "No scattering factor for res %s atom %s",
                    (*myatoms->resinfo[ri].name),
                    (*myatoms->atomname[ai]));
            GMX_THROW(gmx::InvalidInputError(errbuf));
        }

        res.push_back(ibuf);
    }
    gmx_atomprop_destroy(aps);

    fprintf(out, "\n[ pairs ]\n");
    for (ai = 0; (ai < myatoms->nr); ai++)
    {
        for (aj = ai; (aj < myatoms->nr); aj++)
        {
            fprintf(out, "%5d  %5d  3  %s  %s\n",
                    ai+1, aj+1, res[ai].c_str(), res[aj].c_str());
        }
    }
    gmx_fio_fclose(out);
}

static void modify_topology(const char *top, const char *itp)
{
    if (NULL != top)
    {
        FILE  *fp;
        char **strings = NULL;
        int    i, nlines;

        nlines = get_file(top, &strings);
        if (0 == nlines)
        {
            char errbuf[100];
            sprintf(errbuf, "Topology file %s is empty", top);
            GMX_THROW(gmx::InvalidInputError(errbuf));
        }
        fp = gmx_fio_fopen(top, "w");
        for (i = 0; (i < nlines); i++)
        {
            if (NULL != strstr(strings[i], "bonds"))
            {
                fprintf(fp, "\n; Line added by %s\n#include \"%s\"\n\n",
                        ShortProgram(), itp);
            }
            fprintf(fp, "%s\n", strings[i]);
            sfree(strings[i]);
        }
        sfree(strings);
        gmx_fio_fclose(fp);
    }
}

int gmx_genpr(int argc, char *argv[])
{
    const char        *desc[] = {
        "[TT]genrestr[tt] produces an include file for a topology containing",
        "restraints of different kinds.[PAR]",
        "[BB]Position restraints[bb] can be generated as",
        "a list of atom numbers and three force constants for the",
        "[IT]x[it]-, [IT]y[it]-, and [IT]z[it]-direction. A single isotropic force constant may",
        "be given on the command line instead of three components.[PAR]",
        "WARNING: position restraints only work for the one molecule at a time.",
        "Position restraints are interactions within molecules, therefore",
        "they should be included within the correct [TT][ moleculetype ][tt]",
        "block in the topology. Since the atom numbers in every moleculetype",
        "in the topology start at 1 and the numbers in the input file for",
        "[TT]genrestr[tt] number consecutively from 1, [TT]genrestr[tt] will only",
        "produce a useful file for the first molecule.[PAR]",

        "The [TT]-of[tt] option produces an index file that can be used for",
        "[BB]freezing atoms[bb]. In this case, the input file must be a [TT].pdb[tt] file.[PAR]",

        "With the [TT]-disre[tt] option, half a matrix of [BB]distance restraints[bb]",
        "is generated instead of position restraints. With this matrix, that",
        "one typically would apply to C[GRK]alpha[grk] atoms in a protein, one can",
        "maintain the overall conformation of a protein without tieing it to",
        "a specific position (as with position restraints).[PAR]",

        "The [TT]-matrix[tt] option generates half a matrix of pairs of type 3",
        "between either atoms or residues, for use with [BB]SAXS/WAXS refinement[bb]",
        "based on the Debye equation (double sum over atom pairs). Citations",
        "will be added soon. For this purpose an additional file can be read",
        "containing scattering factor tables using the [TT]-d[tt] option, which is used",
        "to generate additional input to the topology. This file is expected",
        "to be compatible with the [TT]mdrun[tt] WAXS/SAXS refinement option.",
        "No checking for consistency is done, but [TT]grompp[tt] will crash if",
        "there is an inconsistency between the scattering factor file and",
        "the generated itp file. The topology file can optionially be modified in",
        "order to include the restraint file."
    };
    static rvec        fc            = {1000.0, 1000.0, 1000.0};
    static real        freeze_level  = 0.0;
    static real        disre_dist    = 0.1;
    static real        disre_frac    = 0.0;
    static real        disre_up2     = 1.0;
    static gmx_bool    bDisre        = FALSE;
    static gmx_bool    bConstr       = FALSE;
    static real        cutoff        = -1.0;
    static int         resmax        = 0;
    const char        *enum_matrix[] = {NULL, "none", "residue", "atom", NULL};

    t_pargs            pa[] = {
        { "-fc", FALSE, etRVEC, {fc},
          "Force constants (kJ/mol nm^2)" },
        { "-freeze", FALSE, etREAL, {&freeze_level},
          "If the [TT]-of[tt] option or this one is given an index file will be written containing atom numbers of all atoms that have a B-factor less than the level given here" },
        { "-disre", FALSE, etBOOL, {&bDisre},
          "Generate a distance restraint matrix for all the atoms in index" },
        { "-disre_dist", FALSE, etREAL, {&disre_dist},
          "Distance range around the actual distance for generating distance restraints" },
        { "-disre_frac", FALSE, etREAL, {&disre_frac},
          "Fraction of distance to be used as interval rather than a fixed distance. If the fraction of the distance that you specify here is less than the distance given in the previous option, that one is used instead." },
        { "-disre_up2", FALSE, etREAL, {&disre_up2},
          "Distance between upper bound for distance restraints, and the distance at which the force becomes constant (see manual)" },
        { "-cutoff", FALSE, etREAL, {&cutoff},
          "Only generate distance restraints for atoms pairs within cutoff (nm)" },
        { "-constr", FALSE, etBOOL, {&bConstr},
          "Generate a constraint matrix rather than distance restraints. Constraints of type 2 will be generated that do generate exclusions." },
        { "-matrix", FALSE, etENUM, {enum_matrix},
          "Generate a matrix of restraints between either atoms (option atom) or one center of mass vsite per residue (option residue) or none (default) for use with SAXS/WAXS refinement. See also the [TT]-d[tt] option," },
        { "-resmax", FALSE, etINT, {&resmax},
          "Highest residue number to take into account when making a matrix with the residue option. If <= 0 all residues will be taken into account." }
    };
#define npargs ((sizeof(pa)/sizeof(pa[0])))

    output_env_t   oenv;
    int            i, j, k;
    FILE          *out;
    int            igrp;
    real           d, dd;
    real           lo, hi;
    atom_id       *ind_grp;
    const char    *xfn, *nfn;
    char          *gn_grp;
    char           title[STRLEN];
    matrix         box;
    t_atoms       *myatoms;
    gmx_bool       bFreeze;
    rvec           dx, *x = NULL, *v = NULL;
    int            nmatrix;
    t_filenm       fnm[] = {
        { efSTX, "-f",  NULL,    ffREAD },
        { efNDX, "-n",  NULL,    ffOPTRD },
        { efITP, "-o",  "posre", ffWRITE },
        { efNDX, "-of", "freeze", ffOPTWR },
        { efXML, "-d", "sfactor", ffOPTRD },
        { efTOP, "-p", "topol",   ffOPTWR },
        { efSTO, "-oc",  NULL,    ffOPTWR },
    };
#define NFILE ((sizeof(fnm)/sizeof(fnm[0])))

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, npargs, pa,
                           (sizeof(desc)/sizeof(desc[0])), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }

    bFreeze = opt2bSet("-of", NFILE, fnm) || opt2parg_bSet("-freeze",
                                                           npargs, pa);
    bDisre  = bDisre || opt2parg_bSet("-disre_dist", npargs, pa);
    xfn     = opt2fn_null("-f", NFILE, fnm);
    nfn     = opt2fn_null("-n", NFILE, fnm);
    nmatrix = nenum(enum_matrix)-1;

    char errbuf[STRLEN];
    if (( nfn == NULL ) && ( xfn == NULL))
    {
        sprintf(errbuf, "no index file and no structure file suplied");
        GMX_THROW(gmx::InvalidInputError(errbuf));
    }

    if ((disre_frac < 0) || (disre_frac >= 1))
    {
        sprintf(errbuf, "disre_frac should be between 0 and 1 instead of %f", disre_frac);
        GMX_THROW(gmx::InvalidInputError(errbuf));
    }
    if (disre_dist < 0)
    {
        sprintf(errbuf, "disre_dist should be >= 0 instead of %f", disre_dist);
        GMX_THROW(gmx::InvalidInputError(errbuf));
    }
    snew(myatoms, 1);
    if (xfn != NULL)
    {
        int natoms = 0;
        get_stx_coordnum(xfn, &natoms);
        if (natoms <= 0)
        {
            sprintf(errbuf, "Something weird with input structure in %20s", xfn);
            GMX_THROW(gmx::InvalidInputError(errbuf));
        }
        init_t_atoms(myatoms, natoms, TRUE);
        snew(x, natoms);
        snew(v, natoms);
        read_stx_conf(xfn, title, myatoms, x, v, NULL, box);
    }

    if (bFreeze)
    {
        if (NULL == myatoms->pdbinfo)
        {
            sprintf(errbuf, "No B-factors in input file %s, use a pdb file next time.", xfn);
            GMX_THROW(gmx::InvalidInputError(errbuf));
        }

        out = opt2FILE("-of", NFILE, fnm, "w");
        fprintf(out, "[ freeze ]\n");
        for (i = 0; (i < myatoms->nr); i++)
        {
            if (myatoms->pdbinfo[i].bfac <= freeze_level)
            {
                fprintf(out, "%d\n", i+1);
            }
        }
        gmx_fio_fclose(out);
    }
    else if ((bDisre || bConstr) && x)
    {
        printf("Select group to generate %s matrix from\n",
               bConstr ? "constraint" : "distance restraint");
        get_index(myatoms, nfn, 1, &igrp, &ind_grp, &gn_grp);

        out = ftp2FILE(efITP, NFILE, fnm, "w");
        if (bConstr)
        {
            fprintf(out, "; constraints for %s of %s\n\n", gn_grp, title);
            fprintf(out, "[ constraints ]\n");
            fprintf(out, ";%4s %5s %1s %10s\n", "i", "j", "tp", "dist");
        }
        else
        {
            fprintf(out, "; distance restraints for %s of %s\n\n", gn_grp, title);
            fprintf(out, "[ distance_restraints ]\n");
            fprintf(out, ";%4s %5s %1s %5s %10s %10s %10s %10s %10s\n", "i", "j", "?",
                    "label", "funct", "lo", "up1", "up2", "weight");
        }
        for (i = k = 0; i < igrp; i++)
        {
            for (j = i+1; j < igrp; j++, k++)
            {
                rvec_sub(x[ind_grp[i]], x[ind_grp[j]], dx);
                d = norm(dx);
                if (bConstr)
                {
                    fprintf(out, "%5d %5d %1d %10g\n", ind_grp[i]+1, ind_grp[j]+1, 2, d);
                }
                else
                {
                    if (cutoff < 0 || d < cutoff)
                    {
                        if (disre_frac > 0)
                        {
                            dd = std::min(disre_dist, disre_frac*d);
                        }
                        else
                        {
                            dd = disre_dist;
                        }
                        lo = std::max(0.0, (double)(d-dd));
                        hi = d+dd;
                        fprintf(out, "%5d %5d %1d %5d %10d %10g %10g %10g %10g\n",
                                ind_grp[i]+1, ind_grp[j]+1, 1, k, 1,
                                lo, hi, hi+1, 1.0);
                    }
                }
            }
        }
        gmx_fio_fclose(out);
    }
    else if ((1 == nmatrix) || (2 == nmatrix))
    {
        gmx::ScatteringFactorTable sft;
        const char                *dat = opt2fn("-d", NFILE, fnm);

        if (!sft.read(dat))
        {
            GMX_THROW(gmx::FileIOError(dat));
        }
        if (1 == nmatrix)
        {
            gen_residue_matrix(ftp2fn(efITP, NFILE, fnm),
                               sft, myatoms, &x, resmax);

            write_sto_conf(opt2fn("-oc", NFILE, fnm), title,
                           myatoms, x, NULL, epbcXYZ, box);
        }
        else if (2 == nmatrix)
        {
            gen_atom_matrix(ftp2fn(efITP, NFILE, fnm),
                            sft, myatoms);
        }
    }
    else
    {
        printf("Select group to position restrain\n");
        get_index(myatoms, nfn, 1, &igrp, &ind_grp, &gn_grp);

        out = ftp2FILE(efITP, NFILE, fnm, "w");
        fprintf(out, "; position restraints for %s of %s\n\n", gn_grp, title);
        fprintf(out, "[ position_restraints ]\n");
        fprintf(out, ";%3s %5s %9s %10s %10s\n", "i", "funct", "fcx", "fcy", "fcz");
        for (i = 0; i < igrp; i++)
        {
            fprintf(out, "%4d %4d %10g %10g %10g\n",
                    ind_grp[i]+1, 1, fc[XX], fc[YY], fc[ZZ]);
        }
        gmx_fio_fclose(out);
    }
    modify_topology(opt2fn_null("-p", NFILE, fnm), opt2fn("-o", NFILE, fnm));

    if (xfn)
    {
        sfree(x);
        sfree(v);
    }
    printf("Successfully generated %s.\n", ftp2fn(efITP, NFILE, fnm));

    return 0;
}
