/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "pdb2gmx.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/genhydro.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/hizzie.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/specbond.h"
#include "gromacs/gmxpreprocess/ter_db.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/xlate.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

#define RTP_MAXCHAR 5
typedef struct {
    char gmx[RTP_MAXCHAR+2];
    char main[RTP_MAXCHAR+2];
    char nter[RTP_MAXCHAR+2];
    char cter[RTP_MAXCHAR+2];
    char bter[RTP_MAXCHAR+2];
} rtprename_t;


static const char *res2bb_notermini(const char *name,
                                    int nrr, const rtprename_t *rr)
{
    /* NOTE: This function returns the main building block name,
     *       it does not take terminal renaming into account.
     */
    int i;

    i = 0;
    while (i < nrr && gmx_strcasecmp(name, rr[i].gmx) != 0)
    {
        i++;
    }

    return (i < nrr ? rr[i].main : name);
}

static const char *select_res(int nr, int resnr,
                              const char *name[], const char *expl[],
                              const char *title,
                              int nrr, const rtprename_t *rr)
{
    int sel = 0;

    printf("Which %s type do you want for residue %d\n", title, resnr+1);
    for (sel = 0; (sel < nr); sel++)
    {
        printf("%d. %s (%s)\n",
               sel, expl[sel], res2bb_notermini(name[sel], nrr, rr));
    }
    printf("\nType a number:"); fflush(stdout);

    if (scanf("%d", &sel) != 1)
    {
        gmx_fatal(FARGS, "Answer me for res %s %d!", title, resnr+1);
    }

    return name[sel];
}

static const char *get_asptp(int resnr, int nrr, const rtprename_t *rr)
{
    enum {
        easp, easpH, easpNR
    };
    const char *lh[easpNR]   = { "ASP", "ASPH" };
    const char *expl[easpNR] = {
        "Not protonated (charge -1)",
        "Protonated (charge 0)"
    };

    return select_res(easpNR, resnr, lh, expl, "ASPARTIC ACID", nrr, rr);
}

static const char *get_glutp(int resnr, int nrr, const rtprename_t *rr)
{
    enum {
        eglu, egluH, egluNR
    };
    const char *lh[egluNR]   = { "GLU", "GLUH" };
    const char *expl[egluNR] = {
        "Not protonated (charge -1)",
        "Protonated (charge 0)"
    };

    return select_res(egluNR, resnr, lh, expl, "GLUTAMIC ACID", nrr, rr);
}

static const char *get_glntp(int resnr, int nrr, const rtprename_t *rr)
{
    enum {
        egln, eglnH, eglnNR
    };
    const char *lh[eglnNR]   = { "GLN", "QLN" };
    const char *expl[eglnNR] = {
        "Not protonated (charge 0)",
        "Protonated (charge +1)"
    };

    return select_res(eglnNR, resnr, lh, expl, "GLUTAMINE", nrr, rr);
}

static const char *get_lystp(int resnr, int nrr, const rtprename_t *rr)
{
    enum {
        elys, elysH, elysNR
    };
    const  char *lh[elysNR]   = { "LYSN", "LYS" };
    const char  *expl[elysNR] = {
        "Not protonated (charge 0)",
        "Protonated (charge +1)"
    };

    return select_res(elysNR, resnr, lh, expl, "LYSINE", nrr, rr);
}

static const char *get_argtp(int resnr, int nrr, const rtprename_t *rr)
{
    enum {
        earg, eargH, eargNR
    };
    const  char *lh[eargNR]   = { "ARGN", "ARG" };
    const char  *expl[eargNR] = {
        "Not protonated (charge 0)",
        "Protonated (charge +1)"
    };

    return select_res(eargNR, resnr, lh, expl, "ARGININE", nrr, rr);
}

static const char *get_histp(int resnr, int nrr, const rtprename_t *rr)
{
    const char *expl[ehisNR] = {
        "H on ND1 only",
        "H on NE2 only",
        "H on ND1 and NE2",
        "Coupled to Heme"
    };

    return select_res(ehisNR, resnr, hh, expl, "HISTIDINE", nrr, rr);
}

static void read_rtprename(const char *fname, FILE *fp,
                           int *nrtprename, rtprename_t **rtprename)
{
    char         line[STRLEN], buf[STRLEN];
    int          n;
    rtprename_t *rr;
    int          ncol, nc;

    n  = *nrtprename;
    rr = *rtprename;

    ncol = 0;
    while (get_a_line(fp, line, STRLEN))
    {
        srenew(rr, n+1);
        /* line is NULL-terminated and length<STRLEN, so final arg cannot overflow.
         * For other args, we read up to 6 chars (so we can detect if the length is > 5).
         * Note that the buffer length has been increased to 7 to allow this,
         * so we just need to make sure the strings have zero-length initially.
         */
        rr[n].gmx[0]  = '\0';
        rr[n].main[0] = '\0';
        rr[n].nter[0] = '\0';
        rr[n].cter[0] = '\0';
        rr[n].bter[0] = '\0';
        nc            = sscanf(line, "%6s %6s %6s %6s %6s %s",
                               rr[n].gmx, rr[n].main, rr[n].nter, rr[n].cter, rr[n].bter, buf);
        if (strlen(rr[n].gmx) > RTP_MAXCHAR || strlen(rr[n].main) > RTP_MAXCHAR ||
            strlen(rr[n].nter) > RTP_MAXCHAR || strlen(rr[n].cter) > RTP_MAXCHAR || strlen(rr[n].bter) > RTP_MAXCHAR)
        {
            gmx_fatal(FARGS, "Residue renaming database '%s' has strings longer than %d chars in first 5 columns:\n%s",
                      fname, RTP_MAXCHAR, line);
        }

        if (ncol == 0)
        {
            if (nc != 2 && nc != 5)
            {
                gmx_fatal(FARGS, "Residue renaming database '%s' has %d columns instead of %d, %d or %d", fname, ncol, 2, 5);
            }
            ncol = nc;
        }
        else if (nc != ncol)
        {
            gmx_fatal(FARGS, "A line in residue renaming database '%s' has %d columns, while previous lines have %d columns", fname, nc, ncol);
        }

        if (nc == 2)
        {
            /* This file does not have special termini names, copy them from main */
            strcpy(rr[n].nter, rr[n].main);
            strcpy(rr[n].cter, rr[n].main);
            strcpy(rr[n].bter, rr[n].main);
        }

        n++;
    }

    *nrtprename = n;
    *rtprename  = rr;
}

static char *search_resrename(int nrr, rtprename_t *rr,
                              const char *name,
                              gmx_bool bStart, gmx_bool bEnd,
                              gmx_bool bCompareFFRTPname)
{
    char *nn;
    int   i;

    nn = nullptr;

    i = 0;
    while (i < nrr && ((!bCompareFFRTPname && strcmp(name, rr[i].gmx)  != 0) ||
                       ( bCompareFFRTPname && strcmp(name, rr[i].main) != 0)))
    {
        i++;
    }

    /* If found in the database, rename this residue's rtp building block,
     * otherwise keep the old name.
     */
    if (i < nrr)
    {
        if (bStart && bEnd)
        {
            nn = rr[i].bter;
        }
        else if (bStart)
        {
            nn = rr[i].nter;
        }
        else if (bEnd)
        {
            nn = rr[i].cter;
        }
        else
        {
            nn = rr[i].main;
        }

        if (nn[0] == '-')
        {
            gmx_fatal(FARGS, "In the chosen force field there is no residue type for '%s'%s", name, bStart ? ( bEnd ? " as a standalone (starting & ending) residue" : " as a starting terminus") : (bEnd ? " as an ending terminus" : ""));
        }
    }

    return nn;
}


static void rename_resrtp(t_atoms *pdba, int nterpairs, int *r_start, int *r_end,
                          int nrr, rtprename_t *rr, t_symtab *symtab,
                          gmx_bool bVerbose)
{
    int      r, j;
    gmx_bool bStart, bEnd;
    char    *nn;
    gmx_bool bFFRTPTERRNM;

    bFFRTPTERRNM = (getenv("GMX_NO_FFRTP_TER_RENAME") == nullptr);

    for (r = 0; r < pdba->nres; r++)
    {
        bStart = FALSE;
        bEnd   = FALSE;
        for (j = 0; j < nterpairs; j++)
        {
            if (r == r_start[j])
            {
                bStart = TRUE;
            }
        }
        for (j = 0; j < nterpairs; j++)
        {
            if (r == r_end[j])
            {
                bEnd = TRUE;
            }
        }

        nn = search_resrename(nrr, rr, *pdba->resinfo[r].rtp, bStart, bEnd, FALSE);

        if (bFFRTPTERRNM && nn == nullptr && (bStart || bEnd))
        {
            /* This is a terminal residue, but the residue name,
             * currently stored in .rtp, is not a standard residue name,
             * but probably a force field specific rtp name.
             * Check if we need to rename it because it is terminal.
             */
            nn = search_resrename(nrr, rr,
                                  *pdba->resinfo[r].rtp, bStart, bEnd, TRUE);
        }

        if (nn != nullptr && strcmp(*pdba->resinfo[r].rtp, nn) != 0)
        {
            if (bVerbose)
            {
                printf("Changing rtp entry of residue %d %s to '%s'\n",
                       pdba->resinfo[r].nr, *pdba->resinfo[r].name, nn);
            }
            pdba->resinfo[r].rtp = put_symtab(symtab, nn);
        }
    }
}

static void pdbres_to_gmxrtp(t_atoms *pdba)
{
    int i;

    for (i = 0; (i < pdba->nres); i++)
    {
        if (pdba->resinfo[i].rtp == nullptr)
        {
            pdba->resinfo[i].rtp = pdba->resinfo[i].name;
        }
    }
}

static void rename_pdbres(t_atoms *pdba, const char *oldnm, const char *newnm,
                          gmx_bool bFullCompare, t_symtab *symtab)
{
    char *resnm;
    int   i;

    for (i = 0; (i < pdba->nres); i++)
    {
        resnm = *pdba->resinfo[i].name;
        if ((bFullCompare && (gmx_strcasecmp(resnm, oldnm) == 0)) ||
            (!bFullCompare && strstr(resnm, oldnm) != nullptr))
        {
            /* Rename the residue name (not the rtp name) */
            pdba->resinfo[i].name = put_symtab(symtab, newnm);
        }
    }
}

static void rename_bb(t_atoms *pdba, const char *oldnm, const char *newnm,
                      gmx_bool bFullCompare, t_symtab *symtab)
{
    char *bbnm;
    int   i;

    for (i = 0; (i < pdba->nres); i++)
    {
        /* We have not set the rtp name yes, use the residue name */
        bbnm = *pdba->resinfo[i].name;
        if ((bFullCompare && (gmx_strcasecmp(bbnm, oldnm) == 0)) ||
            (!bFullCompare && strstr(bbnm, oldnm) != nullptr))
        {
            /* Change the rtp builing block name */
            pdba->resinfo[i].rtp = put_symtab(symtab, newnm);
        }
    }
}

static void rename_bbint(t_atoms *pdba, const char *oldnm,
                         const char *gettp(int, int, const rtprename_t *),
                         gmx_bool bFullCompare,
                         t_symtab *symtab,
                         int nrr, const rtprename_t *rr)
{
    int         i;
    const char *ptr;
    char       *bbnm;

    for (i = 0; i < pdba->nres; i++)
    {
        /* We have not set the rtp name yes, use the residue name */
        bbnm = *pdba->resinfo[i].name;
        if ((bFullCompare && (strcmp(bbnm, oldnm) == 0)) ||
            (!bFullCompare && strstr(bbnm, oldnm) != nullptr))
        {
            ptr                  = gettp(i, nrr, rr);
            pdba->resinfo[i].rtp = put_symtab(symtab, ptr);
        }
    }
}

static void check_occupancy(t_atoms *atoms, const char *filename, gmx_bool bVerbose)
{
    int i, ftp;
    int nzero   = 0;
    int nnotone = 0;

    ftp = fn2ftp(filename);
    if (!atoms->pdbinfo || ((ftp != efPDB) && (ftp != efBRK) && (ftp != efENT)))
    {
        fprintf(stderr, "No occupancies in %s\n", filename);
    }
    else
    {
        for (i = 0; (i < atoms->nr); i++)
        {
            if (atoms->pdbinfo[i].occup != 1)
            {
                if (bVerbose)
                {
                    fprintf(stderr, "Occupancy for atom %s%d-%s is %f rather than 1\n",
                            *atoms->resinfo[atoms->atom[i].resind].name,
                            atoms->resinfo[ atoms->atom[i].resind].nr,
                            *atoms->atomname[i],
                            atoms->pdbinfo[i].occup);
                }
                if (atoms->pdbinfo[i].occup == 0)
                {
                    nzero++;
                }
                else
                {
                    nnotone++;
                }
            }
        }
        if (nzero == atoms->nr)
        {
            fprintf(stderr, "All occupancy fields zero. This is probably not an X-Ray structure\n");
        }
        else if ((nzero > 0) || (nnotone > 0))
        {
            fprintf(stderr,
                    "\n"
                    "WARNING: there were %d atoms with zero occupancy and %d atoms with\n"
                    "         occupancy unequal to one (out of %d atoms). Check your pdb file.\n"
                    "\n",
                    nzero, nnotone, atoms->nr);
        }
        else
        {
            fprintf(stderr, "All occupancies are one\n");
        }
    }
}

static void write_posres(char *fn, t_atoms *pdba, real fc)
{
    FILE *fp;
    int   i;

    fp = gmx_fio_fopen(fn, "w");
    fprintf(fp,
            "; In this topology include file, you will find position restraint\n"
            "; entries for all the heavy atoms in your original pdb file.\n"
            "; This means that all the protons which were added by pdb2gmx are\n"
            "; not restrained.\n"
            "\n"
            "[ position_restraints ]\n"
            "; %4s%6s%8s%8s%8s\n", "atom", "type", "fx", "fy", "fz"
            );
    for (i = 0; (i < pdba->nr); i++)
    {
        if (!is_hydrogen(*pdba->atomname[i]) && !is_dummymass(*pdba->atomname[i]))
        {
            fprintf(fp, "%6d%6d  %g  %g  %g\n", i+1, 1, fc, fc, fc);
        }
    }
    gmx_fio_fclose(fp);
}

static int read_pdball(const char *inf, const char *outf, char *title,
                       t_atoms *atoms, rvec **x,
                       int *ePBC, matrix box, gmx_bool bRemoveH,
                       t_symtab *symtab, gmx_residuetype_t *rt, const char *watres,
                       gmx_atomprop_t aps, gmx_bool bVerbose)
/* Read a pdb file. (containing proteins) */
{
    int  natom, new_natom, i;

    /* READ IT */
    printf("Reading %s...\n", inf);
    t_topology *top;
    snew(top, 1);
    read_tps_conf(inf, top, ePBC, x, nullptr, box, FALSE);
    strncpy(title, *top->name, STRLEN);
    title[STRLEN-1] = '\0';
    *atoms          = top->atoms;
    sfree(top);
    natom = atoms->nr;
    if (atoms->pdbinfo == nullptr)
    {
        snew(atoms->pdbinfo, atoms->nr);
    }
    if (fn2ftp(inf) == efPDB)
    {
        get_pdb_atomnumber(atoms, aps);
    }
    if (bRemoveH)
    {
        new_natom = 0;
        for (i = 0; i < atoms->nr; i++)
        {
            if (!is_hydrogen(*atoms->atomname[i]))
            {
                atoms->atom[new_natom]     = atoms->atom[i];
                atoms->atomname[new_natom] = atoms->atomname[i];
                atoms->pdbinfo[new_natom]  = atoms->pdbinfo[i];
                copy_rvec((*x)[i], (*x)[new_natom]);
                new_natom++;
            }
        }
        atoms->nr = new_natom;
        natom     = new_natom;
    }

    printf("Read");
    if (title[0])
    {
        printf(" '%s',", title);
    }
    printf(" %d atoms\n", natom);

    /* Rename residues */
    rename_pdbres(atoms, "HOH", watres, FALSE, symtab);
    rename_pdbres(atoms, "SOL", watres, FALSE, symtab);
    rename_pdbres(atoms, "WAT", watres, FALSE, symtab);

    rename_atoms("xlateat.dat", nullptr,
                 atoms, symtab, nullptr, TRUE, rt, TRUE, bVerbose);

    if (natom == 0)
    {
        return 0;
    }

    if (outf)
    {
        write_sto_conf(outf, title, atoms, *x, nullptr, *ePBC, box);
    }

    return natom;
}

static void process_chain(t_atoms *pdba, rvec *x,
                          gmx_bool bTrpU, gmx_bool bPheU, gmx_bool bTyrU,
                          gmx_bool bLysMan, gmx_bool bAspMan, gmx_bool bGluMan,
                          gmx_bool bHisMan, gmx_bool bArgMan, gmx_bool bGlnMan,
                          real angle, real distance, t_symtab *symtab,
                          int nrr, const rtprename_t *rr)
{
    /* Rename aromatics, lys, asp and histidine */
    if (bTyrU)
    {
        rename_bb(pdba, "TYR", "TYRU", FALSE, symtab);
    }
    if (bTrpU)
    {
        rename_bb(pdba, "TRP", "TRPU", FALSE, symtab);
    }
    if (bPheU)
    {
        rename_bb(pdba, "PHE", "PHEU", FALSE, symtab);
    }
    if (bLysMan)
    {
        rename_bbint(pdba, "LYS", get_lystp, FALSE, symtab, nrr, rr);
    }
    if (bArgMan)
    {
        rename_bbint(pdba, "ARG", get_argtp, FALSE, symtab, nrr, rr);
    }
    if (bGlnMan)
    {
        rename_bbint(pdba, "GLN", get_glntp, FALSE, symtab, nrr, rr);
    }
    if (bAspMan)
    {
        rename_bbint(pdba, "ASP", get_asptp, FALSE, symtab, nrr, rr);
    }
    else
    {
        rename_bb(pdba, "ASPH", "ASP", FALSE, symtab);
    }
    if (bGluMan)
    {
        rename_bbint(pdba, "GLU", get_glutp, FALSE, symtab, nrr, rr);
    }
    else
    {
        rename_bb(pdba, "GLUH", "GLU", FALSE, symtab);
    }

    if (!bHisMan)
    {
        set_histp(pdba, x, angle, distance);
    }
    else
    {
        rename_bbint(pdba, "HIS", get_histp, TRUE, symtab, nrr, rr);
    }

    /* Initialize the rtp builing block names with the residue names
     * for the residues that have not been processed above.
     */
    pdbres_to_gmxrtp(pdba);

    /* Now we have all rtp names set.
     * The rtp names will conform to Gromacs naming,
     * unless the input pdb file contained one or more force field specific
     * rtp names as residue names.
     */
}

/* struct for sorting the atoms from the pdb file */
typedef struct {
    int  resnr;  /* residue number               */
    int  j;      /* database order index         */
    int  index;  /* original atom number         */
    char anm1;   /* second letter of atom name   */
    char altloc; /* alternate location indicator */
} t_pdbindex;

static int pdbicomp(const void *a, const void *b)
{
    t_pdbindex *pa, *pb;
    int         d;

    pa = (t_pdbindex *)a;
    pb = (t_pdbindex *)b;

    d = (pa->resnr - pb->resnr);
    if (d == 0)
    {
        d = (pa->j - pb->j);
        if (d == 0)
        {
            d = (pa->anm1 - pb->anm1);
            if (d == 0)
            {
                d = (pa->altloc - pb->altloc);
            }
        }
    }

    return d;
}

static void sort_pdbatoms(t_restp restp[],
                          int natoms, t_atoms **pdbaptr, rvec **x,
                          t_blocka *block, char ***gnames)
{
    t_atoms     *pdba, *pdbnew;
    rvec       **xnew;
    int          i, j;
    t_restp     *rptr;
    t_pdbindex  *pdbi;
    int         *a;
    char        *atomnm;

    pdba   = *pdbaptr;
    natoms = pdba->nr;
    pdbnew = nullptr;
    snew(xnew, 1);
    snew(pdbi, natoms);

    for (i = 0; i < natoms; i++)
    {
        atomnm = *pdba->atomname[i];
        rptr   = &restp[pdba->atom[i].resind];
        for (j = 0; (j < rptr->natom); j++)
        {
            if (gmx_strcasecmp(atomnm, *(rptr->atomname[j])) == 0)
            {
                break;
            }
        }
        if (j == rptr->natom)
        {
            char buf[STRLEN];

            sprintf(buf,
                    "Atom %s in residue %s %d was not found in rtp entry %s with %d atoms\n"
                    "while sorting atoms.\n%s", atomnm,
                    *pdba->resinfo[pdba->atom[i].resind].name,
                    pdba->resinfo[pdba->atom[i].resind].nr,
                    rptr->resname,
                    rptr->natom,
                    is_hydrogen(atomnm) ?
                    "\nFor a hydrogen, this can be a different protonation state, or it\n"
                    "might have had a different number in the PDB file and was rebuilt\n"
                    "(it might for instance have been H3, and we only expected H1 & H2).\n"
                    "Note that hydrogens might have been added to the entry for the N-terminus.\n"
                    "Remove this hydrogen or choose a different protonation state to solve it.\n"
                    "Option -ignh will ignore all hydrogens in the input." : ".");
            gmx_fatal(FARGS, buf);
        }
        /* make shadow array to be sorted into indexgroup */
        pdbi[i].resnr  = pdba->atom[i].resind;
        pdbi[i].j      = j;
        pdbi[i].index  = i;
        pdbi[i].anm1   = atomnm[1];
        pdbi[i].altloc = pdba->pdbinfo[i].altloc;
    }
    qsort(pdbi, natoms, (size_t)sizeof(pdbi[0]), pdbicomp);

    /* pdba is sorted in pdbnew using the pdbi index */
    snew(a, natoms);
    snew(pdbnew, 1);
    init_t_atoms(pdbnew, natoms, TRUE);
    snew(*xnew, natoms);
    pdbnew->nr   = pdba->nr;
    pdbnew->nres = pdba->nres;
    sfree(pdbnew->resinfo);
    pdbnew->resinfo = pdba->resinfo;
    for (i = 0; i < natoms; i++)
    {
        pdbnew->atom[i]     = pdba->atom[pdbi[i].index];
        pdbnew->atomname[i] = pdba->atomname[pdbi[i].index];
        pdbnew->pdbinfo[i]  = pdba->pdbinfo[pdbi[i].index];
        copy_rvec((*x)[pdbi[i].index], (*xnew)[i]);
        /* make indexgroup in block */
        a[i] = pdbi[i].index;
    }
    /* clean up */
    sfree(pdba->atomname);
    sfree(pdba->atom);
    sfree(pdba->pdbinfo);
    sfree(pdba);
    sfree(*x);
    /* copy the sorted pdbnew back to pdba */
    *pdbaptr = pdbnew;
    *x       = *xnew;
    add_grp(block, gnames, natoms, a, "prot_sort");
    sfree(xnew);
    sfree(a);
    sfree(pdbi);
}

static int remove_duplicate_atoms(t_atoms *pdba, rvec x[], gmx_bool bVerbose)
{
    int        i, j, oldnatoms, ndel;
    t_resinfo *ri;

    printf("Checking for duplicate atoms....\n");
    oldnatoms    = pdba->nr;
    ndel         = 0;
    /* NOTE: pdba->nr is modified inside the loop */
    for (i = 1; (i < pdba->nr); i++)
    {
        /* compare 'i' and 'i-1', throw away 'i' if they are identical
           this is a 'while' because multiple alternate locations can be present */
        while ( (i < pdba->nr) &&
                (pdba->atom[i-1].resind == pdba->atom[i].resind) &&
                (strcmp(*pdba->atomname[i-1], *pdba->atomname[i]) == 0) )
        {
            ndel++;
            if (bVerbose)
            {
                ri = &pdba->resinfo[pdba->atom[i].resind];
                printf("deleting duplicate atom %4s  %s%4d%c",
                       *pdba->atomname[i], *ri->name, ri->nr, ri->ic);
                if (ri->chainid && (ri->chainid != ' '))
                {
                    printf(" ch %c", ri->chainid);
                }
                if (pdba->pdbinfo)
                {
                    if (pdba->pdbinfo[i].atomnr)
                    {
                        printf("  pdb nr %4d", pdba->pdbinfo[i].atomnr);
                    }
                    if (pdba->pdbinfo[i].altloc && (pdba->pdbinfo[i].altloc != ' '))
                    {
                        printf("  altloc %c", pdba->pdbinfo[i].altloc);
                    }
                }
                printf("\n");
            }
            pdba->nr--;
            /* We can not free, since it might be in the symtab */
            /* sfree(pdba->atomname[i]); */
            for (j = i; j < pdba->nr; j++)
            {
                pdba->atom[j]     = pdba->atom[j+1];
                pdba->atomname[j] = pdba->atomname[j+1];
                if (pdba->pdbinfo)
                {
                    pdba->pdbinfo[j]  = pdba->pdbinfo[j+1];
                }
                copy_rvec(x[j+1], x[j]);
            }
            srenew(pdba->atom,     pdba->nr);
            /* srenew(pdba->atomname, pdba->nr); */
            srenew(pdba->pdbinfo,  pdba->nr);
        }
    }
    if (pdba->nr != oldnatoms)
    {
        printf("Now there are %d atoms. Deleted %d duplicates.\n", pdba->nr, ndel);
    }

    return pdba->nr;
}

static void find_nc_ter(t_atoms *pdba, int r0, int r1, int *r_start, int *r_end,
                        gmx_residuetype_t *rt)
{
    int         i;
    const char *p_startrestype;
    const char *p_restype;
    int         nstartwarn, nendwarn;

    *r_start = -1;
    *r_end   = -1;

    nstartwarn = 0;
    nendwarn   = 0;

    /* Find the starting terminus (typially N or 5') */
    for (i = r0; i < r1 && *r_start == -1; i++)
    {
        gmx_residuetype_get_type(rt, *pdba->resinfo[i].name, &p_startrestype);
        if (!gmx_strcasecmp(p_startrestype, "Protein") || !gmx_strcasecmp(p_startrestype, "DNA") || !gmx_strcasecmp(p_startrestype, "RNA") )
        {
            printf("Identified residue %s%d as a starting terminus.\n", *pdba->resinfo[i].name, pdba->resinfo[i].nr);
            *r_start = i;
        }
        else
        {
            if (nstartwarn < 5)
            {
                printf("Warning: Starting residue %s%d in chain not identified as Protein/RNA/DNA.\n", *pdba->resinfo[i].name, pdba->resinfo[i].nr);
            }
            if (nstartwarn == 5)
            {
                printf("More than 5 unidentified residues at start of chain - disabling further warnings.\n");
            }
            nstartwarn++;
        }
    }

    if (*r_start >= 0)
    {
        /* Go through the rest of the residues, check that they are the same class, and identify the ending terminus. */
        for (i = *r_start; i < r1; i++)
        {
            gmx_residuetype_get_type(rt, *pdba->resinfo[i].name, &p_restype);
            if (!gmx_strcasecmp(p_restype, p_startrestype) && nendwarn == 0)
            {
                *r_end = i;
            }
            else
            {
                if (nendwarn < 5)
                {
                    printf("Warning: Residue %s%d in chain has different type (%s) from starting residue %s%d (%s).\n",
                           *pdba->resinfo[i].name, pdba->resinfo[i].nr, p_restype,
                           *pdba->resinfo[*r_start].name, pdba->resinfo[*r_start].nr, p_startrestype);
                }
                if (nendwarn == 5)
                {
                    printf("More than 5 unidentified residues at end of chain - disabling further warnings.\n");
                }
                nendwarn++;
            }
        }
    }

    if (*r_end >= 0)
    {
        printf("Identified residue %s%d as a ending terminus.\n", *pdba->resinfo[*r_end].name, pdba->resinfo[*r_end].nr);
    }
}


enum SplittingType
{
    SPLIT_ID_OR_TER,
    SPLIT_ID_AND_TER,
    SPLIT_ID_ONLY,
    SPLIT_TER_ONLY,
    SPLIT_INTERACTIVE
};

static SplittingType getSplittingType(const char *chainsep)
{
    SplittingType splitting = SPLIT_TER_ONLY; /* keep compiler happy */

    /* Be a bit flexible to catch typos */
    if (!strncmp(chainsep, "id_o", 4))
    {
        /* For later interactive splitting we tentatively assign new chain numbers at either changing id or ter records */
        splitting = SPLIT_ID_OR_TER;
        printf("Splitting chemical chains based on TER records or chain id changing.\n");
    }
    else if (!strncmp(chainsep, "int", 3))
    {
        /* For later interactive splitting we tentatively assign new chain numbers at either changing id or ter records */
        splitting = SPLIT_INTERACTIVE;
        printf("Splitting chemical chains interactively.\n");
    }
    else if (!strncmp(chainsep, "id_a", 4))
    {
        splitting = SPLIT_ID_AND_TER;
        printf("Splitting chemical chains based on TER records and chain id changing.\n");
    }
    else if (strlen(chainsep) == 2 && !strncmp(chainsep, "id", 4))
    {
        splitting = SPLIT_ID_ONLY;
        printf("Splitting chemical chains based on changing chain id only (ignoring TER records).\n");
    }
    else if (chainsep[0] == 't')
    {
        splitting = SPLIT_TER_ONLY;
        printf("Splitting chemical chains based on TER records only (ignoring chain id).\n");
    }
    else
    {
        gmx_fatal(FARGS, "Unidentified setting for chain separation: %s\n", chainsep);
    }
    return splitting;
}

static void
modify_chain_numbers(t_atoms *       pdba,
                     const char *    chainsep)
{
    int           i;
    char          old_prev_chainid;
    char          old_this_chainid;
    int           old_prev_chainnum;
    int           old_this_chainnum;
    t_resinfo    *ri;
    char          select[STRLEN];
    int           new_chainnum;
    int           this_atomnum;
    int           prev_atomnum;
    const char *  prev_atomname;
    const char *  this_atomname;
    const char *  prev_resname;
    const char *  this_resname;
    int           prev_resnum;
    int           this_resnum;
    char          prev_chainid;
    char          this_chainid;

    SplittingType splitting = getSplittingType(chainsep);

    /* The default chain enumeration is based on TER records only, which is reflected in chainnum below */

    old_prev_chainid  = '?';
    old_prev_chainnum = -1;
    new_chainnum      = -1;

    this_atomname       = nullptr;
    this_atomnum        = -1;
    this_resname        = nullptr;
    this_resnum         = -1;
    this_chainid        = '?';

    for (i = 0; i < pdba->nres; i++)
    {
        ri                 = &pdba->resinfo[i];
        old_this_chainid   = ri->chainid;
        old_this_chainnum  = ri->chainnum;

        prev_atomname      = this_atomname;
        prev_atomnum       = this_atomnum;
        prev_resname       = this_resname;
        prev_resnum        = this_resnum;
        prev_chainid       = this_chainid;

        this_atomname      = *(pdba->atomname[i]);
        this_atomnum       = (pdba->pdbinfo != nullptr) ? pdba->pdbinfo[i].atomnr : i+1;
        this_resname       = *ri->name;
        this_resnum        = ri->nr;
        this_chainid       = ri->chainid;

        switch (splitting)
        {
            case SPLIT_ID_OR_TER:
                if (old_this_chainid != old_prev_chainid || old_this_chainnum != old_prev_chainnum)
                {
                    new_chainnum++;
                }
                break;

            case SPLIT_ID_AND_TER:
                if (old_this_chainid != old_prev_chainid && old_this_chainnum != old_prev_chainnum)
                {
                    new_chainnum++;
                }
                break;

            case SPLIT_ID_ONLY:
                if (old_this_chainid != old_prev_chainid)
                {
                    new_chainnum++;
                }
                break;

            case SPLIT_TER_ONLY:
                if (old_this_chainnum != old_prev_chainnum)
                {
                    new_chainnum++;
                }
                break;
            case SPLIT_INTERACTIVE:
                if (old_this_chainid != old_prev_chainid || old_this_chainnum != old_prev_chainnum)
                {
                    if (i > 0)
                    {
                        printf("Split the chain (and introduce termini) between residue %s%d (chain id '%c', atom %d %s)\n"
                               "and residue %s%d (chain id '%c', atom %d %s) ? [n/y]\n",
                               prev_resname, prev_resnum, prev_chainid, prev_atomnum, prev_atomname,
                               this_resname, this_resnum, this_chainid, this_atomnum, this_atomname);

                        if (nullptr == fgets(select, STRLEN-1, stdin))
                        {
                            gmx_fatal(FARGS, "Error reading from stdin");
                        }
                    }
                    if (i == 0 || select[0] == 'y')
                    {
                        new_chainnum++;
                    }
                }
                break;
            default:
                gmx_fatal(FARGS, "Internal inconsistency - this shouldn't happen...");
                break;
        }
        old_prev_chainid  = old_this_chainid;
        old_prev_chainnum = old_this_chainnum;

        ri->chainnum = new_chainnum;
    }
}


typedef struct {
    char     chainid;
    char     chainnum;
    int      start;
    int      natom;
    gmx_bool bAllWat;
    int      nterpairs;
    int     *chainstart;
} t_pdbchain;

typedef struct {
    char          chainid;
    int           chainnum;
    gmx_bool      bAllWat;
    int           nterpairs;
    int          *chainstart;
    t_hackblock **ntdb;
    t_hackblock **ctdb;
    int          *r_start;
    int          *r_end;
    t_atoms      *pdba;
    rvec         *x;
} t_chain;

int gmx_pdb2gmx(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] reads a [REF].pdb[ref] (or [REF].gro[ref]) file, reads",
        "some database files, adds hydrogens to the molecules and generates",
        "coordinates in GROMACS (GROMOS), or optionally [REF].pdb[ref], format",
        "and a topology in GROMACS format.",
        "These files can subsequently be processed to generate a run input file.",
        "[PAR]",
        "[THISMODULE] will search for force fields by looking for",
        "a [TT]forcefield.itp[tt] file in subdirectories [TT]<forcefield>.ff[tt]",
        "of the current working directory and of the GROMACS library directory",
        "as inferred from the path of the binary or the [TT]GMXLIB[tt] environment",
        "variable.",
        "By default the forcefield selection is interactive,",
        "but you can use the [TT]-ff[tt] option to specify one of the short names",
        "in the list on the command line instead. In that case [THISMODULE] just looks",
        "for the corresponding [TT]<forcefield>.ff[tt] directory.",
        "[PAR]",
        "After choosing a force field, all files will be read only from",
        "the corresponding force field directory.",
        "If you want to modify or add a residue types, you can copy the force",
        "field directory from the GROMACS library directory to your current",
        "working directory. If you want to add new protein residue types,",
        "you will need to modify [TT]residuetypes.dat[tt] in the library directory",
        "or copy the whole library directory to a local directory and set",
        "the environment variable [TT]GMXLIB[tt] to the name of that directory.",
        "Check Chapter 5 of the manual for more information about file formats.",
        "[PAR]",

        "Note that a [REF].pdb[ref] file is nothing more than a file format, and it",
        "need not necessarily contain a protein structure. Every kind of",
        "molecule for which there is support in the database can be converted.",
        "If there is no support in the database, you can add it yourself.[PAR]",

        "The program has limited intelligence, it reads a number of database",
        "files, that allow it to make special bonds (Cys-Cys, Heme-His, etc.),",
        "if necessary this can be done manually. The program can prompt the",
        "user to select which kind of LYS, ASP, GLU, CYS or HIS residue is",
        "desired. For Lys the choice is between neutral (two protons on NZ) or",
        "protonated (three protons, default), for Asp and Glu unprotonated",
        "(default) or protonated, for His the proton can be either on ND1,",
        "on NE2 or on both. By default these selections are done automatically.",
        "For His, this is based on an optimal hydrogen bonding",
        "conformation. Hydrogen bonds are defined based on a simple geometric",
        "criterion, specified by the maximum hydrogen-donor-acceptor angle",
        "and donor-acceptor distance, which are set by [TT]-angle[tt] and",
        "[TT]-dist[tt] respectively.[PAR]",

        "The protonation state of N- and C-termini can be chosen interactively",
        "with the [TT]-ter[tt] flag.  Default termini are ionized (NH3+ and COO-),",
        "respectively.  Some force fields support zwitterionic forms for chains of",
        "one residue, but for polypeptides these options should NOT be selected.",
        "The AMBER force fields have unique forms for the terminal residues,",
        "and these are incompatible with the [TT]-ter[tt] mechanism. You need",
        "to prefix your N- or C-terminal residue names with \"N\" or \"C\"",
        "respectively to use these forms, making sure you preserve the format",
        "of the coordinate file. Alternatively, use named terminating residues",
        "(e.g. ACE, NME).[PAR]",

        "The separation of chains is not entirely trivial since the markup",
        "in user-generated PDB files frequently varies and sometimes it",
        "is desirable to merge entries across a TER record, for instance",
        "if you want a disulfide bridge or distance restraints between",
        "two protein chains or if you have a HEME group bound to a protein.",
        "In such cases multiple chains should be contained in a single",
        "[TT]moleculetype[tt] definition.",
        "To handle this, [THISMODULE] uses two separate options.",
        "First, [TT]-chainsep[tt] allows you to choose when a new chemical chain should",
        "start, and termini added when applicable. This can be done based on the",
        "existence of TER records, when the chain id changes, or combinations of either",
        "or both of these. You can also do the selection fully interactively.",
        "In addition, there is a [TT]-merge[tt] option that controls how multiple chains",
        "are merged into one moleculetype, after adding all the chemical termini (or not).",
        "This can be turned off (no merging), all non-water chains can be merged into a",
        "single molecule, or the selection can be done interactively.[PAR]",

        "[THISMODULE] will also check the occupancy field of the [REF].pdb[ref] file.",
        "If any of the occupancies are not one, indicating that the atom is",
        "not resolved well in the structure, a warning message is issued.",
        "When a [REF].pdb[ref] file does not originate from an X-ray structure determination",
        "all occupancy fields may be zero. Either way, it is up to the user",
        "to verify the correctness of the input data (read the article!).[PAR]",

        "During processing the atoms will be reordered according to GROMACS",
        "conventions. With [TT]-n[tt] an index file can be generated that",
        "contains one group reordered in the same way. This allows you to",
        "convert a GROMOS trajectory and coordinate file to GROMOS. There is",
        "one limitation: reordering is done after the hydrogens are stripped",
        "from the input and before new hydrogens are added. This means that",
        "you should not use [TT]-ignh[tt].[PAR]",

        "The [REF].gro[ref] and [TT].g96[tt] file formats do not support chain",
        "identifiers. Therefore it is useful to enter a [REF].pdb[ref] file name at",
        "the [TT]-o[tt] option when you want to convert a multi-chain [REF].pdb[ref] file.",
        "[PAR]",

        "The option [TT]-vsite[tt] removes hydrogen and fast improper dihedral",
        "motions. Angular and out-of-plane motions can be removed by changing",
        "hydrogens into virtual sites and fixing angles, which fixes their",
        "position relative to neighboring atoms. Additionally, all atoms in the",
        "aromatic rings of the standard amino acids (i.e. PHE, TRP, TYR and HIS)",
        "can be converted into virtual sites, eliminating the fast improper dihedral",
        "fluctuations in these rings. [BB]Note[bb] that in this case all other hydrogen",
        "atoms are also converted to virtual sites. The mass of all atoms that are",
        "converted into virtual sites, is added to the heavy atoms.[PAR]",
        "Also slowing down of dihedral motion can be done with [TT]-heavyh[tt]",
        "done by increasing the hydrogen-mass by a factor of 4. This is also",
        "done for water hydrogens to slow down the rotational motion of water.",
        "The increase in mass of the hydrogens is subtracted from the bonded",
        "(heavy) atom so that the total mass of the system remains the same."
    };


    FILE             *fp, *top_file, *top_file2, *itp_file = nullptr;
    int               natom, nres;
    t_atoms           pdba_all, *pdba;
    t_atoms          *atoms;
    t_resinfo        *ri;
    t_blocka         *block;
    int               chain, nch, maxch, nwaterchain;
    t_pdbchain       *pdb_ch;
    t_chain          *chains, *cc;
    char              select[STRLEN];
    int               nincl, nmol;
    char            **incls;
    t_mols           *mols;
    char            **gnames;
    int               ePBC;
    matrix            box;
    rvec              box_space;
    int               i, j, k, l, nrtp;
    int              *swap_index, si;
    t_restp          *restp;
    t_hackblock      *ah;
    t_symtab          symtab;
    gpp_atomtype_t    atype;
    gmx_residuetype_t*rt;
    const char       *top_fn;
    char              fn[256], itp_fn[STRLEN], posre_fn[STRLEN], buf_fn[STRLEN];
    char              molname[STRLEN], title[STRLEN];
    char             *c, forcefield[STRLEN], ffdir[STRLEN];
    char              ffname[STRLEN], suffix[STRLEN], buf[STRLEN];
    char             *watermodel;
    const char       *watres;
    int               nrtpf;
    char            **rtpf;
    int               nrrn;
    char            **rrn;
    int               nrtprename;
    rtprename_t      *rtprename = nullptr;
    int               nah, nNtdb, nCtdb, ntdblist;
    t_hackblock      *ntdb, *ctdb, **tdblist;
    int               nssbonds;
    t_ssbond         *ssbonds;
    rvec             *pdbx, *x;
    gmx_bool          bVsites = FALSE, bWat, bPrevWat = FALSE, bITP, bVsiteAromatics = FALSE;
    real              mHmult  = 0;
    t_hackblock      *hb_chain;
    t_restp          *restp_chain;
    gmx_output_env_t *oenv;
    const char       *p_restype;
    int               rc;
    int               this_atomnum;
    int               prev_atomnum;
    const char     *  prev_atomname;
    const char     *  this_atomname;
    const char     *  prev_resname;
    const char     *  this_resname;
    int               prev_resnum;
    int               this_resnum;
    char              prev_chainid;
    char              this_chainid;
    int               prev_chainnumber;
    int               this_chainnumber;
    int               nid_used;
    int               this_chainstart;
    int               prev_chainstart;
    gmx_bool          bMerged;
    int               nchainmerges;

    gmx_atomprop_t    aps;

    t_filenm          fnm[] = {
        { efSTX, "-f", "eiwit.pdb", ffREAD  },
        { efSTO, "-o", "conf",      ffWRITE },
        { efTOP, nullptr, nullptr,        ffWRITE },
        { efITP, "-i", "posre",     ffWRITE },
        { efNDX, "-n", "clean",     ffOPTWR },
        { efSTO, "-q", "clean.pdb", ffOPTWR }
    };
#define NFILE asize(fnm)

    gmx_bool           bNewRTP        = FALSE;
    gmx_bool           bInter         = FALSE, bCysMan = FALSE;
    gmx_bool           bLysMan        = FALSE, bAspMan = FALSE, bGluMan = FALSE, bHisMan = FALSE;
    gmx_bool           bGlnMan        = FALSE, bArgMan = FALSE;
    gmx_bool           bTerMan        = FALSE, bUnA = FALSE, bHeavyH = FALSE;
    gmx_bool           bSort          = TRUE, bAllowMissing = FALSE, bRemoveH = FALSE;
    gmx_bool           bDeuterate     = FALSE, bVerbose = FALSE, bChargeGroups = TRUE, bCmap = TRUE;
    gmx_bool           bRenumRes      = FALSE, bRTPresname = FALSE;
    real               angle          = 135.0, distance = 0.3, posre_fc = 1000;
    real               long_bond_dist = 0.25, short_bond_dist = 0.05;
    const char        *vsitestr[]     = { nullptr, "none", "hydrogens", "aromatics", nullptr };
    const char        *watstr[]       = { nullptr, "select", "none", "spc", "spce", "tip3p", "tip4p", "tip5p", "tips3p", nullptr };
    const char        *chainsep[]     = { nullptr, "id_or_ter", "id_and_ter", "ter", "id", "interactive", nullptr };
    const char        *merge[]        = {nullptr, "no", "all", "interactive", nullptr };
    const char        *ff             = "select";

    t_pargs            pa[] = {
        { "-newrtp", FALSE, etBOOL, {&bNewRTP},
          "HIDDENWrite the residue database in new format to [TT]new.rtp[tt]"},
        { "-lb",     FALSE, etREAL, {&long_bond_dist},
          "HIDDENLong bond warning distance" },
        { "-sb",     FALSE, etREAL, {&short_bond_dist},
          "HIDDENShort bond warning distance" },
        { "-chainsep", FALSE, etENUM, {chainsep},
          "Condition in PDB files when a new chain should be started (adding termini)" },
        { "-merge",  FALSE, etENUM, {&merge},
          "Merge multiple chains into a single [moleculetype]" },
        { "-ff",     FALSE, etSTR,  {&ff},
          "Force field, interactive by default. Use [TT]-h[tt] for information." },
        { "-water",  FALSE, etENUM, {watstr},
          "Water model to use" },
        { "-inter",  FALSE, etBOOL, {&bInter},
          "Set the next 8 options to interactive"},
        { "-ss",     FALSE, etBOOL, {&bCysMan},
          "Interactive SS bridge selection" },
        { "-ter",    FALSE, etBOOL, {&bTerMan},
          "Interactive termini selection, instead of charged (default)" },
        { "-lys",    FALSE, etBOOL, {&bLysMan},
          "Interactive lysine selection, instead of charged" },
        { "-arg",    FALSE, etBOOL, {&bArgMan},
          "Interactive arginine selection, instead of charged" },
        { "-asp",    FALSE, etBOOL, {&bAspMan},
          "Interactive aspartic acid selection, instead of charged" },
        { "-glu",    FALSE, etBOOL, {&bGluMan},
          "Interactive glutamic acid selection, instead of charged" },
        { "-gln",    FALSE, etBOOL, {&bGlnMan},
          "Interactive glutamine selection, instead of neutral" },
        { "-his",    FALSE, etBOOL, {&bHisMan},
          "Interactive histidine selection, instead of checking H-bonds" },
        { "-angle",  FALSE, etREAL, {&angle},
          "Minimum hydrogen-donor-acceptor angle for a H-bond (degrees)" },
        { "-dist",   FALSE, etREAL, {&distance},
          "Maximum donor-acceptor distance for a H-bond (nm)" },
        { "-una",    FALSE, etBOOL, {&bUnA},
          "Select aromatic rings with united CH atoms on phenylalanine, "
          "tryptophane and tyrosine" },
        { "-sort",   FALSE, etBOOL, {&bSort},
          "HIDDENSort the residues according to database, turning this off is dangerous as charge groups might be broken in parts" },
        { "-ignh",   FALSE, etBOOL, {&bRemoveH},
          "Ignore hydrogen atoms that are in the coordinate file" },
        { "-missing", FALSE, etBOOL, {&bAllowMissing},
          "Continue when atoms are missing, dangerous" },
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Be slightly more verbose in messages" },
        { "-posrefc", FALSE, etREAL, {&posre_fc},
          "Force constant for position restraints" },
        { "-vsite",  FALSE, etENUM, {vsitestr},
          "Convert atoms to virtual sites" },
        { "-heavyh", FALSE, etBOOL, {&bHeavyH},
          "Make hydrogen atoms heavy" },
        { "-deuterate", FALSE, etBOOL, {&bDeuterate},
          "Change the mass of hydrogens to 2 amu" },
        { "-chargegrp", TRUE, etBOOL, {&bChargeGroups},
          "Use charge groups in the [REF].rtp[ref] file"  },
        { "-cmap", TRUE, etBOOL, {&bCmap},
          "Use cmap torsions (if enabled in the [REF].rtp[ref] file)"  },
        { "-renum", TRUE, etBOOL, {&bRenumRes},
          "Renumber the residues consecutively in the output"  },
        { "-rtpres", TRUE, etBOOL, {&bRTPresname},
          "Use [REF].rtp[ref] entry names as residue names"  }
    };
#define NPARGS asize(pa)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }

    /* Force field selection, interactive or direct */
    choose_ff(strcmp(ff, "select") == 0 ? nullptr : ff,
              forcefield, sizeof(forcefield),
              ffdir, sizeof(ffdir));

    if (strlen(forcefield) > 0)
    {
        strcpy(ffname, forcefield);
        ffname[0] = toupper(ffname[0]);
    }
    else
    {
        gmx_fatal(FARGS, "Empty forcefield string");
    }

    printf("\nUsing the %s force field in directory %s\n\n",
           ffname, ffdir);

    choose_watermodel(watstr[0], ffdir, &watermodel);

    if (bInter)
    {
        /* if anything changes here, also change description of -inter */
        bCysMan = TRUE;
        bTerMan = TRUE;
        bLysMan = TRUE;
        bArgMan = TRUE;
        bAspMan = TRUE;
        bGluMan = TRUE;
        bGlnMan = TRUE;
        bHisMan = TRUE;
    }

    if (bHeavyH)
    {
        mHmult = 4.0;
    }
    else if (bDeuterate)
    {
        mHmult = 2.0;
    }
    else
    {
        mHmult = 1.0;
    }

    /* parse_common_args ensures vsitestr has been selected, but
       clang-static-analyzer needs clues to know that */
    GMX_ASSERT(vsitestr[0], "-vsite default wasn't processed correctly");
    switch (vsitestr[0][0])
    {
        case 'n': /* none */
            bVsites         = FALSE;
            bVsiteAromatics = FALSE;
            break;
        case 'h': /* hydrogens */
            bVsites         = TRUE;
            bVsiteAromatics = FALSE;
            break;
        case 'a': /* aromatics */
            bVsites         = TRUE;
            bVsiteAromatics = TRUE;
            break;
        default:
            gmx_fatal(FARGS, "DEATH HORROR in $s (%d): vsitestr[0]='%s'",
                      __FILE__, __LINE__, vsitestr[0]);
    } /* end switch */

    /* Open the symbol table */
    open_symtab(&symtab);

    /* Residue type database */
    gmx_residuetype_init(&rt);

    /* Read residue renaming database(s), if present */
    nrrn = fflib_search_file_end(ffdir, ".r2b", FALSE, &rrn);

    nrtprename = 0;
    rtprename  = nullptr;
    for (i = 0; i < nrrn; i++)
    {
        fp = fflib_open(rrn[i]);
        read_rtprename(rrn[i], fp, &nrtprename, &rtprename);
        gmx_ffclose(fp);
        sfree(rrn[i]);
    }
    sfree(rrn);

    /* Add all alternative names from the residue renaming database to the list of recognized amino/nucleic acids. */
    for (i = 0; i < nrtprename; i++)
    {
        rc = gmx_residuetype_get_type(rt, rtprename[i].gmx, &p_restype);

        /* Only add names if the 'standard' gromacs/iupac base name was found */
        if (rc == 0)
        {
            gmx_residuetype_add(rt, rtprename[i].main, p_restype);
            gmx_residuetype_add(rt, rtprename[i].nter, p_restype);
            gmx_residuetype_add(rt, rtprename[i].cter, p_restype);
            gmx_residuetype_add(rt, rtprename[i].bter, p_restype);
        }
    }

    clear_mat(box);
    if (watermodel != nullptr && (strstr(watermodel, "4p") ||
                                  strstr(watermodel, "4P")))
    {
        watres = "HO4";
    }
    else if (watermodel != nullptr && (strstr(watermodel, "5p") ||
                                       strstr(watermodel, "5P")))
    {
        watres = "HO5";
    }
    else
    {
        watres = "HOH";
    }

    aps   = gmx_atomprop_init();
    natom = read_pdball(opt2fn("-f", NFILE, fnm), opt2fn_null("-q", NFILE, fnm), title,
                        &pdba_all, &pdbx, &ePBC, box, bRemoveH, &symtab, rt, watres,
                        aps, bVerbose);

    if (natom == 0)
    {
        gmx_fatal(FARGS, "No atoms found in pdb file %s\n", opt2fn("-f", NFILE, fnm));
    }

    printf("Analyzing pdb file\n");
    nwaterchain = 0;

    modify_chain_numbers(&pdba_all, chainsep[0]);

    nchainmerges        = 0;

    this_atomname       = nullptr;
    this_atomnum        = -1;
    this_resname        = nullptr;
    this_resnum         = -1;
    this_chainid        = '?';
    this_chainnumber    = -1;
    this_chainstart     = 0;
    /* Keep the compiler happy */
    prev_chainstart     = 0;

    nch   = 0;
    maxch = 16;
    snew(pdb_ch, maxch);

    bMerged = FALSE;
    for (i = 0; (i < natom); i++)
    {
        ri = &pdba_all.resinfo[pdba_all.atom[i].resind];

        /* TODO this should live in a helper object, and consolidate
           that with code in modify_chain_numbers */
        prev_atomname      = this_atomname;
        prev_atomnum       = this_atomnum;
        prev_resname       = this_resname;
        prev_resnum        = this_resnum;
        prev_chainid       = this_chainid;
        prev_chainnumber   = this_chainnumber;
        if (!bMerged)
        {
            prev_chainstart    = this_chainstart;
        }

        this_atomname      = *pdba_all.atomname[i];
        this_atomnum       = (pdba_all.pdbinfo != nullptr) ? pdba_all.pdbinfo[i].atomnr : i+1;
        this_resname       = *ri->name;
        this_resnum        = ri->nr;
        this_chainid       = ri->chainid;
        this_chainnumber   = ri->chainnum;

        bWat = gmx_strcasecmp(*ri->name, watres) == 0;
        if ((i == 0) || (this_chainnumber != prev_chainnumber) || (bWat != bPrevWat))
        {
            this_chainstart = pdba_all.atom[i].resind;

            bMerged = FALSE;
            if (i > 0 && !bWat)
            {
                if (!strncmp(merge[0], "int", 3))
                {
                    printf("Merge chain ending with residue %s%d (chain id '%c', atom %d %s) and chain starting with\n"
                           "residue %s%d (chain id '%c', atom %d %s) into a single moleculetype (keeping termini)? [n/y]\n",
                           prev_resname, prev_resnum, prev_chainid, prev_atomnum, prev_atomname,
                           this_resname, this_resnum, this_chainid, this_atomnum, this_atomname);

                    if (nullptr == fgets(select, STRLEN-1, stdin))
                    {
                        gmx_fatal(FARGS, "Error reading from stdin");
                    }
                    bMerged = (select[0] == 'y');
                }
                else if (!strncmp(merge[0], "all", 3))
                {
                    bMerged = TRUE;
                }
            }

            if (bMerged)
            {
                pdb_ch[nch-1].chainstart[pdb_ch[nch-1].nterpairs] =
                    pdba_all.atom[i].resind - prev_chainstart;
                pdb_ch[nch-1].nterpairs++;
                srenew(pdb_ch[nch-1].chainstart, pdb_ch[nch-1].nterpairs+1);
                nchainmerges++;
            }
            else
            {
                /* set natom for previous chain */
                if (nch > 0)
                {
                    pdb_ch[nch-1].natom = i-pdb_ch[nch-1].start;
                }
                if (bWat)
                {
                    nwaterchain++;
                    ri->chainid = ' ';
                }
                /* check if chain identifier was used before */
                for (j = 0; (j < nch); j++)
                {
                    if (pdb_ch[j].chainid != ' ' && pdb_ch[j].chainid == ri->chainid)
                    {
                        printf("WARNING: Chain identifier '%c' is used in two non-sequential blocks.\n"
                               "They will be treated as separate chains unless you reorder your file.\n",
                               ri->chainid);
                    }
                }
                // TODO This is too convoluted. Use a std::vector
                if (nch == maxch)
                {
                    maxch += 16;
                    srenew(pdb_ch, maxch);
                }
                pdb_ch[nch].chainid  = ri->chainid;
                pdb_ch[nch].chainnum = ri->chainnum;
                pdb_ch[nch].start    = i;
                pdb_ch[nch].bAllWat  = bWat;
                if (bWat)
                {
                    pdb_ch[nch].nterpairs = 0;
                }
                else
                {
                    pdb_ch[nch].nterpairs = 1;
                }
                snew(pdb_ch[nch].chainstart, pdb_ch[nch].nterpairs+1);
                /* modified [nch] to [0] below */
                pdb_ch[nch].chainstart[0] = 0;
                nch++;
            }
        }
        bPrevWat = bWat;
    }
    pdb_ch[nch-1].natom = natom-pdb_ch[nch-1].start;

    /* set all the water blocks at the end of the chain */
    snew(swap_index, nch);
    j = 0;
    for (i = 0; i < nch; i++)
    {
        if (!pdb_ch[i].bAllWat)
        {
            swap_index[j] = i;
            j++;
        }
    }
    for (i = 0; i < nch; i++)
    {
        if (pdb_ch[i].bAllWat)
        {
            swap_index[j] = i;
            j++;
        }
    }
    if (nwaterchain > 1)
    {
        printf("Moved all the water blocks to the end\n");
    }

    snew(chains, nch);
    /* copy pdb data and x for all chains */
    for (i = 0; (i < nch); i++)
    {
        si                   = swap_index[i];
        chains[i].chainid    = pdb_ch[si].chainid;
        chains[i].chainnum   = pdb_ch[si].chainnum;
        chains[i].bAllWat    = pdb_ch[si].bAllWat;
        chains[i].nterpairs  = pdb_ch[si].nterpairs;
        chains[i].chainstart = pdb_ch[si].chainstart;
        snew(chains[i].ntdb, pdb_ch[si].nterpairs);
        snew(chains[i].ctdb, pdb_ch[si].nterpairs);
        snew(chains[i].r_start, pdb_ch[si].nterpairs);
        snew(chains[i].r_end, pdb_ch[si].nterpairs);

        snew(chains[i].pdba, 1);
        init_t_atoms(chains[i].pdba, pdb_ch[si].natom, TRUE);
        snew(chains[i].x, chains[i].pdba->nr);
        for (j = 0; j < chains[i].pdba->nr; j++)
        {
            chains[i].pdba->atom[j] = pdba_all.atom[pdb_ch[si].start+j];
            snew(chains[i].pdba->atomname[j], 1);
            *chains[i].pdba->atomname[j] =
                gmx_strdup(*pdba_all.atomname[pdb_ch[si].start+j]);
            chains[i].pdba->pdbinfo[j] = pdba_all.pdbinfo[pdb_ch[si].start+j];
            copy_rvec(pdbx[pdb_ch[si].start+j], chains[i].x[j]);
        }
        /* Re-index the residues assuming that the indices are continuous */
        k                    = chains[i].pdba->atom[0].resind;
        nres                 = chains[i].pdba->atom[chains[i].pdba->nr-1].resind - k + 1;
        chains[i].pdba->nres = nres;
        for (j = 0; j < chains[i].pdba->nr; j++)
        {
            chains[i].pdba->atom[j].resind -= k;
        }
        srenew(chains[i].pdba->resinfo, nres);
        for (j = 0; j < nres; j++)
        {
            chains[i].pdba->resinfo[j] = pdba_all.resinfo[k+j];
            snew(chains[i].pdba->resinfo[j].name, 1);
            *chains[i].pdba->resinfo[j].name = gmx_strdup(*pdba_all.resinfo[k+j].name);
            /* make all chain identifiers equal to that of the chain */
            chains[i].pdba->resinfo[j].chainid = pdb_ch[si].chainid;
        }
    }

    if (nchainmerges > 0)
    {
        printf("\nMerged chains into joint molecule definitions at %d places.\n\n",
               nchainmerges);
    }

    printf("There are %d chains and %d blocks of water and "
           "%d residues with %d atoms\n",
           nch-nwaterchain, nwaterchain,
           pdba_all.nres, natom);

    printf("\n  %5s  %4s %6s\n", "chain", "#res", "#atoms");
    for (i = 0; (i < nch); i++)
    {
        printf("  %d '%c' %5d %6d  %s\n",
               i+1, chains[i].chainid ? chains[i].chainid : '-',
               chains[i].pdba->nres, chains[i].pdba->nr,
               chains[i].bAllWat ? "(only water)" : "");
    }
    printf("\n");

    check_occupancy(&pdba_all, opt2fn("-f", NFILE, fnm), bVerbose);

    /* Read atomtypes... */
    atype = read_atype(ffdir, &symtab);

    /* read residue database */
    printf("Reading residue database... (%s)\n", forcefield);
    nrtpf = fflib_search_file_end(ffdir, ".rtp", TRUE, &rtpf);
    nrtp  = 0;
    restp = nullptr;
    for (i = 0; i < nrtpf; i++)
    {
        read_resall(rtpf[i], &nrtp, &restp, atype, &symtab, FALSE);
        sfree(rtpf[i]);
    }
    sfree(rtpf);
    if (bNewRTP)
    {
        /* Not correct with multiple rtp input files with different bonded types */
        fp = gmx_fio_fopen("new.rtp", "w");
        print_resall(fp, nrtp, restp, atype);
        gmx_fio_fclose(fp);
    }

    /* read hydrogen database */
    nah = read_h_db(ffdir, &ah);

    /* Read Termini database... */
    nNtdb = read_ter_db(ffdir, 'n', &ntdb, atype);
    nCtdb = read_ter_db(ffdir, 'c', &ctdb, atype);

    top_fn   = ftp2fn(efTOP, NFILE, fnm);
    top_file = gmx_fio_fopen(top_fn, "w");

    print_top_header(top_file, top_fn, FALSE, ffdir, mHmult);

    nincl = 0;
    nmol  = 0;
    incls = nullptr;
    mols  = nullptr;
    for (chain = 0; (chain < nch); chain++)
    {
        cc = &(chains[chain]);

        /* set pdba, natom and nres to the current chain */
        pdba  = cc->pdba;
        x     = cc->x;
        natom = cc->pdba->nr;
        nres  = cc->pdba->nres;

        if (cc->chainid && ( cc->chainid != ' ' ) )
        {
            printf("Processing chain %d '%c' (%d atoms, %d residues)\n",
                   chain+1, cc->chainid, natom, nres);
        }
        else
        {
            printf("Processing chain %d (%d atoms, %d residues)\n",
                   chain+1, natom, nres);
        }

        process_chain(pdba, x, bUnA, bUnA, bUnA, bLysMan, bAspMan, bGluMan,
                      bHisMan, bArgMan, bGlnMan, angle, distance, &symtab,
                      nrtprename, rtprename);

        cc->chainstart[cc->nterpairs] = pdba->nres;
        j = 0;
        for (i = 0; i < cc->nterpairs; i++)
        {
            find_nc_ter(pdba, cc->chainstart[i], cc->chainstart[i+1],
                        &(cc->r_start[j]), &(cc->r_end[j]), rt);

            if (cc->r_start[j] >= 0 && cc->r_end[j] >= 0)
            {
                j++;
            }
        }
        cc->nterpairs = j;
        if (cc->nterpairs == 0)
        {
            printf("Problem with chain definition, or missing terminal residues.\n"
                   "This chain does not appear to contain a recognized chain molecule.\n"
                   "If this is incorrect, you can edit residuetypes.dat to modify the behavior.\n");
        }

        /* Check for disulfides and other special bonds */
        nssbonds = mk_specbonds(pdba, x, bCysMan, &ssbonds, bVerbose);

        if (nrtprename > 0)
        {
            rename_resrtp(pdba, cc->nterpairs, cc->r_start, cc->r_end, nrtprename, rtprename,
                          &symtab, bVerbose);
        }

        if (debug)
        {
            if (nch == 1)
            {
                sprintf(fn, "chain.pdb");
            }
            else
            {
                sprintf(fn, "chain_%c%d.pdb", cc->chainid, cc->chainnum);
            }
            write_sto_conf(fn, title, pdba, x, nullptr, ePBC, box);
        }


        for (i = 0; i < cc->nterpairs; i++)
        {

            /* Set termini.
             * We first apply a filter so we only have the
             * termini that can be applied to the residue in question
             * (or a generic terminus if no-residue specific is available).
             */
            /* First the N terminus */
            if (nNtdb > 0)
            {
                tdblist = filter_ter(nrtp, restp, nNtdb, ntdb,
                                     *pdba->resinfo[cc->r_start[i]].name,
                                     *pdba->resinfo[cc->r_start[i]].rtp,
                                     &ntdblist);
                if (ntdblist == 0)
                {
                    printf("No suitable end (N or 5') terminus found in database - assuming this residue\n"
                           "is already in a terminus-specific form and skipping terminus selection.\n");
                    cc->ntdb[i] = nullptr;
                }
                else
                {
                    if (bTerMan && ntdblist > 1)
                    {
                        sprintf(select, "Select start terminus type for %s-%d",
                                *pdba->resinfo[cc->r_start[i]].name,
                                pdba->resinfo[cc->r_start[i]].nr);
                        cc->ntdb[i] = choose_ter(ntdblist, tdblist, select);
                    }
                    else
                    {
                        cc->ntdb[i] = tdblist[0];
                    }

                    printf("Start terminus %s-%d: %s\n",
                           *pdba->resinfo[cc->r_start[i]].name,
                           pdba->resinfo[cc->r_start[i]].nr,
                           (cc->ntdb[i])->name);
                    sfree(tdblist);
                }
            }
            else
            {
                cc->ntdb[i] = nullptr;
            }

            /* And the C terminus */
            if (nCtdb > 0)
            {
                tdblist = filter_ter(nrtp, restp, nCtdb, ctdb,
                                     *pdba->resinfo[cc->r_end[i]].name,
                                     *pdba->resinfo[cc->r_end[i]].rtp,
                                     &ntdblist);
                if (ntdblist == 0)
                {
                    printf("No suitable end (C or 3') terminus found in database - assuming this residue\n"
                           "is already in a terminus-specific form and skipping terminus selection.\n");
                    cc->ctdb[i] = nullptr;
                }
                else
                {
                    if (bTerMan && ntdblist > 1)
                    {
                        sprintf(select, "Select end terminus type for %s-%d",
                                *pdba->resinfo[cc->r_end[i]].name,
                                pdba->resinfo[cc->r_end[i]].nr);
                        cc->ctdb[i] = choose_ter(ntdblist, tdblist, select);
                    }
                    else
                    {
                        cc->ctdb[i] = tdblist[0];
                    }
                    printf("End terminus %s-%d: %s\n",
                           *pdba->resinfo[cc->r_end[i]].name,
                           pdba->resinfo[cc->r_end[i]].nr,
                           (cc->ctdb[i])->name);
                    sfree(tdblist);
                }
            }
            else
            {
                cc->ctdb[i] = nullptr;
            }
        }
        /* lookup hackblocks and rtp for all residues */
        get_hackblocks_rtp(&hb_chain, &restp_chain,
                           nrtp, restp, pdba->nres, pdba->resinfo,
                           cc->nterpairs, cc->ntdb, cc->ctdb, cc->r_start, cc->r_end);
        /* ideally, now we would not need the rtp itself anymore, but do
           everything using the hb and restp arrays. Unfortunately, that
           requires some re-thinking of code in gen_vsite.c, which I won't
           do now :( AF 26-7-99 */

        rename_atoms(nullptr, ffdir,
                     pdba, &symtab, restp_chain, FALSE, rt, FALSE, bVerbose);

        match_atomnames_with_rtp(restp_chain, hb_chain, pdba, x, bVerbose);

        if (bSort)
        {
            block = new_blocka();
            snew(gnames, 1);
            sort_pdbatoms(restp_chain, natom, &pdba, &x, block, &gnames);
            remove_duplicate_atoms(pdba, x, bVerbose);
            if (ftp2bSet(efNDX, NFILE, fnm))
            {
                if (bRemoveH)
                {
                    fprintf(stderr, "WARNING: with the -remh option the generated "
                            "index file (%s) might be useless\n"
                            "(the index file is generated before hydrogens are added)",
                            ftp2fn(efNDX, NFILE, fnm));
                }
                write_index(ftp2fn(efNDX, NFILE, fnm), block, gnames, FALSE, 0);
            }
            for (i = 0; i < block->nr; i++)
            {
                sfree(gnames[i]);
            }
            sfree(gnames);
            done_blocka(block);
        }
        else
        {
            fprintf(stderr, "WARNING: "
                    "without sorting no check for duplicate atoms can be done\n");
        }

        /* Generate Hydrogen atoms (and termini) in the sequence */
        printf("Generating any missing hydrogen atoms and/or adding termini.\n");
        natom = add_h(&pdba, &x, nah, ah,
                      cc->nterpairs, cc->ntdb, cc->ctdb, cc->r_start, cc->r_end, bAllowMissing,
                      nullptr, nullptr, TRUE, FALSE);
        printf("Now there are %d residues with %d atoms\n",
               pdba->nres, pdba->nr);
        if (debug)
        {
            write_pdbfile(debug, title, pdba, x, ePBC, box, ' ', 0, nullptr, TRUE);
        }

        if (debug)
        {
            for (i = 0; (i < natom); i++)
            {
                fprintf(debug, "Res %s%d atom %d %s\n",
                        *(pdba->resinfo[pdba->atom[i].resind].name),
                        pdba->resinfo[pdba->atom[i].resind].nr, i+1, *pdba->atomname[i]);
            }
        }

        strcpy(posre_fn, ftp2fn(efITP, NFILE, fnm));

        /* make up molecule name(s) */

        k = (cc->nterpairs > 0 && cc->r_start[0] >= 0) ? cc->r_start[0] : 0;

        gmx_residuetype_get_type(rt, *pdba->resinfo[k].name, &p_restype);

        suffix[0] = '\0';

        if (cc->bAllWat)
        {
            sprintf(molname, "Water");
        }
        else
        {
            this_chainid = cc->chainid;

            /* Add the chain id if we have one */
            if (this_chainid != ' ')
            {
                sprintf(buf, "_chain_%c", this_chainid);
                strcat(suffix, buf);
            }

            /* Check if there have been previous chains with the same id */
            nid_used = 0;
            for (k = 0; k < chain; k++)
            {
                if (cc->chainid == chains[k].chainid)
                {
                    nid_used++;
                }
            }
            /* Add the number for this chain identifier if there are multiple copies */
            if (nid_used > 0)
            {

                sprintf(buf, "%d", nid_used+1);
                strcat(suffix, buf);
            }

            if (strlen(suffix) > 0)
            {
                sprintf(molname, "%s%s", p_restype, suffix);
            }
            else
            {
                strcpy(molname, p_restype);
            }
        }

        if ((nch-nwaterchain > 1) && !cc->bAllWat)
        {
            bITP = TRUE;
            strcpy(itp_fn, top_fn);
            printf("Chain time...\n");
            c = strrchr(itp_fn, '.');
            sprintf(c, "_%s.itp", molname);
            c = strrchr(posre_fn, '.');
            sprintf(c, "_%s.itp", molname);
            if (strcmp(itp_fn, posre_fn) == 0)
            {
                strcpy(buf_fn, posre_fn);
                c  = strrchr(buf_fn, '.');
                *c = '\0';
                sprintf(posre_fn, "%s_pr.itp", buf_fn);
            }

            nincl++;
            srenew(incls, nincl);
            incls[nincl-1] = gmx_strdup(itp_fn);
            itp_file       = gmx_fio_fopen(itp_fn, "w");
        }
        else
        {
            bITP = FALSE;
        }

        srenew(mols, nmol+1);
        if (cc->bAllWat)
        {
            mols[nmol].name = gmx_strdup("SOL");
            mols[nmol].nr   = pdba->nres;
        }
        else
        {
            mols[nmol].name = gmx_strdup(molname);
            mols[nmol].nr   = 1;
        }
        nmol++;

        if (bITP)
        {
            print_top_comment(itp_file, itp_fn, ffdir, TRUE);
        }

        if (cc->bAllWat)
        {
            top_file2 = nullptr;
        }
        else
        if (bITP)
        {
            top_file2 = itp_file;
        }
        else
        {
            top_file2 = top_file;
        }

        pdb2top(top_file2, posre_fn, molname, pdba, &x, atype, &symtab,
                nrtp, restp,
                restp_chain, hb_chain,
                bAllowMissing,
                bVsites, bVsiteAromatics, ffdir,
                mHmult, nssbonds, ssbonds,
                long_bond_dist, short_bond_dist, bDeuterate, bChargeGroups, bCmap,
                bRenumRes, bRTPresname);

        if (!cc->bAllWat)
        {
            write_posres(posre_fn, pdba, posre_fc);
        }

        if (bITP)
        {
            gmx_fio_fclose(itp_file);
        }

        /* pdba and natom have been reassigned somewhere so: */
        cc->pdba = pdba;
        cc->x    = x;

        if (debug)
        {
            if (cc->chainid == ' ')
            {
                sprintf(fn, "chain.pdb");
            }
            else
            {
                sprintf(fn, "chain_%c.pdb", cc->chainid);
            }
            write_sto_conf(fn, "", pdba, x, nullptr, ePBC, box);
        }
    }

    if (watermodel == nullptr)
    {
        for (chain = 0; chain < nch; chain++)
        {
            if (chains[chain].bAllWat)
            {
                gmx_fatal(FARGS, "You have chosen not to include a water model, but there is water in the input file. Select a water model or remove the water from your input file.");
            }
        }
    }
    else
    {
        sprintf(buf_fn, "%s%c%s.itp", ffdir, DIR_SEPARATOR, watermodel);
        if (!fflib_fexist(buf_fn))
        {
            gmx_fatal(FARGS, "The topology file '%s' for the selected water model '%s' can not be found in the force field directory. Select a different water model.",
                      buf_fn, watermodel);
        }
    }

    print_top_mols(top_file, title, ffdir, watermodel, nincl, incls, nmol, mols);
    gmx_fio_fclose(top_file);

    gmx_residuetype_destroy(rt);

    /* now merge all chains back together */
    natom = 0;
    nres  = 0;
    for (i = 0; (i < nch); i++)
    {
        natom += chains[i].pdba->nr;
        nres  += chains[i].pdba->nres;
    }
    snew(atoms, 1);
    init_t_atoms(atoms, natom, FALSE);
    for (i = 0; i < atoms->nres; i++)
    {
        sfree(atoms->resinfo[i].name);
    }
    sfree(atoms->resinfo);
    atoms->nres = nres;
    snew(atoms->resinfo, nres);
    snew(x, natom);
    k = 0;
    l = 0;
    for (i = 0; (i < nch); i++)
    {
        if (nch > 1)
        {
            printf("Including chain %d in system: %d atoms %d residues\n",
                   i+1, chains[i].pdba->nr, chains[i].pdba->nres);
        }
        for (j = 0; (j < chains[i].pdba->nr); j++)
        {
            atoms->atom[k]         = chains[i].pdba->atom[j];
            atoms->atom[k].resind += l; /* l is processed nr of residues */
            atoms->atomname[k]     = chains[i].pdba->atomname[j];
            atoms->resinfo[atoms->atom[k].resind].chainid = chains[i].chainid;
            copy_rvec(chains[i].x[j], x[k]);
            k++;
        }
        for (j = 0; (j < chains[i].pdba->nres); j++)
        {
            atoms->resinfo[l] = chains[i].pdba->resinfo[j];
            if (bRTPresname)
            {
                atoms->resinfo[l].name = atoms->resinfo[l].rtp;
            }
            l++;
        }
    }

    if (nch > 1)
    {
        fprintf(stderr, "Now there are %d atoms and %d residues\n", k, l);
        print_sums(atoms, TRUE);
    }

    fprintf(stderr, "\nWriting coordinate file...\n");
    clear_rvec(box_space);
    if (box[0][0] == 0)
    {
        make_new_box(atoms->nr, x, box, box_space, FALSE);
    }
    write_sto_conf(ftp2fn(efSTO, NFILE, fnm), title, atoms, x, nullptr, ePBC, box);

    printf("\t\t--------- PLEASE NOTE ------------\n");
    printf("You have successfully generated a topology from: %s.\n",
           opt2fn("-f", NFILE, fnm));
    if (watermodel != nullptr)
    {
        printf("The %s force field and the %s water model are used.\n",
               ffname, watermodel);
    }
    else
    {
        printf("The %s force field is used.\n",
               ffname);
    }
    printf("\t\t--------- ETON ESAELP ------------\n");

    return 0;
}
