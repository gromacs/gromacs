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
#include "solvate.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include "sysstuff.h"
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/fileio/confio.h"
#include "macros.h"
#include "gromacs/fileio/futil.h"
#include "atomprop.h"
#include "names.h"
#include "vec.h"
#include "gmx_fatal.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "addconf.h"
#include "read-conformation.h"
#include "gromacs/fileio/pdbio.h"
#include "pbc.h"

#ifdef DEBUG
static void print_stat(rvec *x, int natoms, matrix box)
{
    int  i, m;
    rvec xmin, xmax;
    for (m = 0; (m < DIM); m++)
    {
        xmin[m] = x[0][m];
        xmax[m] = x[0][m];
    }
    for (i = 0; (i < natoms); i++)
    {
        for (m = 0; m < DIM; m++)
        {
            xmin[m] = min(xmin[m], x[i][m]);
            xmax[m] = max(xmax[m], x[i][m]);
        }
    }
    for (m = 0; (m < DIM); m++)
    {
        fprintf(stderr, "DIM %d XMIN %8.3f XMAX %8.3f BOX %8.3f\n",
                m, xmin[m], xmax[m], box[m][m]);
    }
}
#endif

typedef struct {
    char *name;
    int   natoms;
    int   nmol;
    int   i, i0;
    int   res0;
} t_moltypes;

static void sort_molecule(t_atoms **atoms_solvt, rvec *x, rvec *v, real *r)
{
    int         atnr, i, j, moltp = 0, nrmoltypes, resi_o, resi_n, resnr;
    t_moltypes *moltypes;
    t_atoms    *atoms, *newatoms;
    rvec       *newx, *newv = NULL;
    real       *newr;

    fprintf(stderr, "Sorting configuration\n");

    atoms = *atoms_solvt;

    /* copy each residue from *atoms to a molecule in *molecule */
    moltypes   = NULL;
    nrmoltypes = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        if ( (i == 0) || (atoms->atom[i].resind != atoms->atom[i-1].resind) )
        {
            /* see if this was a molecule type we haven't had yet: */
            moltp = NOTSET;
            for (j = 0; (j < nrmoltypes) && (moltp == NOTSET); j++)
            {
                if (strcmp(*(atoms->resinfo[atoms->atom[i].resind].name),
                           moltypes[j].name) == 0)
                {
                    moltp = j;
                }
            }
            if (moltp == NOTSET)
            {
                moltp = nrmoltypes;
                nrmoltypes++;
                srenew(moltypes, nrmoltypes);
                moltypes[moltp].name = *(atoms->resinfo[atoms->atom[i].resind].name);
                atnr                 = 0;
                while ((i+atnr < atoms->nr) &&
                       (atoms->atom[i].resind == atoms->atom[i+atnr].resind))
                {
                    atnr++;
                }
                moltypes[moltp].natoms = atnr;
                moltypes[moltp].nmol   = 0;
            }
            moltypes[moltp].nmol++;
        }
    }

    fprintf(stderr, "Found %d%s molecule type%s:\n",
            nrmoltypes, nrmoltypes == 1 ? "" : " different", nrmoltypes == 1 ? "" : "s");
    for (j = 0; j < nrmoltypes; j++)
    {
        if (j == 0)
        {
            moltypes[j].res0 = 0;
        }
        else
        {
            moltypes[j].res0 = moltypes[j-1].res0+moltypes[j-1].nmol;
        }
        fprintf(stderr, "%7s (%4d atoms): %5d residues\n",
                moltypes[j].name, moltypes[j].natoms, moltypes[j].nmol);
    }

    /* if we have only 1 moleculetype, we don't have to sort */
    if (nrmoltypes > 1)
    {
        /* find out which molecules should go where: */
        moltypes[0].i = moltypes[0].i0 = 0;
        for (j = 1; j < nrmoltypes; j++)
        {
            moltypes[j].i      =
                moltypes[j].i0 =
                    moltypes[j-1].i0+moltypes[j-1].natoms*moltypes[j-1].nmol;
        }

        /* now put them there: */
        snew(newatoms, 1);
        init_t_atoms(newatoms, atoms->nr, FALSE);
        newatoms->nres = atoms->nres;
        snew(newatoms->resinfo, atoms->nres);
        snew(newx, atoms->nr);
        if (v)
        {
            snew(newv, atoms->nr);
        }
        snew(newr, atoms->nr);

        resi_n = 0;
        resnr  = 1;
        j      = 0;
        for (moltp = 0; moltp < nrmoltypes; moltp++)
        {
            i = 0;
            while (i < atoms->nr)
            {
                resi_o = atoms->atom[i].resind;
                if (strcmp(*atoms->resinfo[resi_o].name, moltypes[moltp].name) == 0)
                {
                    /* Copy the residue info */
                    newatoms->resinfo[resi_n]    = atoms->resinfo[resi_o];
                    newatoms->resinfo[resi_n].nr = resnr;
                    /* Copy the atom info */
                    do
                    {
                        newatoms->atom[j]        = atoms->atom[i];
                        newatoms->atomname[j]    = atoms->atomname[i];
                        newatoms->atom[j].resind = resi_n;
                        copy_rvec(x[i], newx[j]);
                        if (v != NULL)
                        {
                            copy_rvec(v[i], newv[j]);
                        }
                        newr[j] = r[i];
                        i++;
                        j++;
                    }
                    while (i < atoms->nr && atoms->atom[i].resind == resi_o);
                    /* Increase the new residue counters */
                    resi_n++;
                    resnr++;
                }
                else
                {
                    /* Skip this residue */
                    do
                    {
                        i++;
                    }
                    while (i < atoms->nr && atoms->atom[i].resind == resi_o);
                }
            }
        }

        /* put them back into the original arrays and throw away temporary arrays */
        sfree(atoms->atomname);
        sfree(atoms->resinfo);
        sfree(atoms->atom);
        sfree(atoms);
        *atoms_solvt = newatoms;
        for (i = 0; i < (*atoms_solvt)->nr; i++)
        {
            copy_rvec(newx[i], x[i]);
            if (v)
            {
                copy_rvec(newv[i], v[i]);
            }
            r[i] = newr[i];
        }
        sfree(newx);
        if (v)
        {
            sfree(newv);
        }
        sfree(newr);
    }
    sfree(moltypes);
}

static void rm_res_pbc(t_atoms *atoms, rvec *x, matrix box)
{
    int  i, start, n, d, nat;
    rvec xcg;

    start = 0;
    nat   = 0;
    clear_rvec(xcg);
    for (n = 0; n < atoms->nr; n++)
    {
        if (!is_hydrogen(*(atoms->atomname[n])))
        {
            nat++;
            rvec_inc(xcg, x[n]);
        }
        if ( (n+1 == atoms->nr) ||
             (atoms->atom[n+1].resind != atoms->atom[n].resind) )
        {
            /* if nat==0 we have only hydrogens in the solvent,
               we take last coordinate as cg */
            if (nat == 0)
            {
                nat = 1;
                copy_rvec(x[n], xcg);
            }
            svmul(1.0/nat, xcg, xcg);
            for (d = 0; d < DIM; d++)
            {
                while (xcg[d] < 0)
                {
                    for (i = start; i <= n; i++)
                    {
                        x[i][d] += box[d][d];
                    }
                    xcg[d] += box[d][d];
                }
                while (xcg[d] >= box[d][d])
                {
                    for (i = start; i <= n; i++)
                    {
                        x[i][d] -= box[d][d];
                    }
                    xcg[d] -= box[d][d];
                }
            }
            start = n+1;
            nat   = 0;
            clear_rvec(xcg);
        }
    }
}

/* Make a new configuration by adding boxes*/
static void make_new_conformation(t_atoms *atoms, rvec *x, rvec *v, real *r, matrix box, ivec n_box)
{
    int     i, ix, iy, iz, m, j, imol, offset;
    rvec    delta;
    int     nmol;

    nmol = n_box[XX]*n_box[YY]*n_box[ZZ];

    /*print message*/
    fprintf(stderr, "Generating configuration\n");
    imol = 0;
    for (ix = 0; (ix < n_box[XX]); ix++)
    {
        delta[XX] = ix*box[XX][XX];
        for (iy = 0; (iy < n_box[YY]); iy++)
        {
            delta[YY] = iy*box[YY][YY];
            for (iz = 0; (iz < n_box[ZZ]); iz++)
            {
                delta[ZZ] = iz*box[ZZ][ZZ];
                offset    = imol*atoms->nr;
                for (i = 0; (i < atoms->nr); i++)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        x[offset+i][m] = delta[m]+x[i][m];
                    }
                    if (v)
                    {
                        for (m = 0; (m < DIM); m++)
                        {
                            v[offset+i][m] = v[i][m];
                        }
                    }
                    r[offset+i] = r[i];
                }
                imol++;
            }
        }
    }
    for (i = 1; (i < nmol); i++)
    {
        int offs    = i*atoms->nr;
        int offsres = i*atoms->nres;
        for (j = 0; (j < atoms->nr); j++)
        {
            atoms->atomname[offs+j]                    = atoms->atomname[j];
            atoms->atom[offs+j].resind                 = atoms->atom[j].resind + offsres;
            atoms->resinfo[atoms->atom[offs+j].resind] =
                atoms->resinfo[atoms->atom[j].resind];
            atoms->resinfo[atoms->atom[offs+j].resind].nr += offsres;
        }
    }
    atoms->nr   *= nmol;
    atoms->nres *= nmol;
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            box[j][i] *= n_box[j];
        }
    }
}

static void add_solv(const char *fn, t_atoms *atoms, rvec **x, rvec **v, real **exclusionDistances,
                     int ePBC, matrix box,
                     gmx_atomprop_t aps,
                     real defaultDistance, real scaleFactor, int *atoms_added,
                     int *residues_added, real rshell, int max_sol,
                     const output_env_t oenv)
{
    int      i, nmol;
    ivec     n_box;
    char     filename[STRLEN];
    char     title_solvt[STRLEN];
    t_atoms *atoms_solvt;
    rvec    *x_solvt, *v_solvt = NULL;
    real    *exclusionDistances_solvt;
    int      ePBC_solvt;
    matrix   box_solvt;
    int      onr, onres;
    char    *lfn;

    lfn = gmxlibfn(fn);
    strncpy(filename, lfn, STRLEN);
    sfree(lfn);
    {
        int natoms;
        get_stx_coordnum(filename, &natoms);
        if (0 == natoms)
        {
            gmx_fatal(FARGS, "No solvent in %s, please check your input\n", filename);
        }
        snew(atoms_solvt, 1);
        init_t_atoms(atoms_solvt, natoms, FALSE);
    }
    snew(x_solvt, atoms_solvt->nr);
    if (v)
    {
        snew(v_solvt, atoms_solvt->nr);
    }
    snew(exclusionDistances_solvt, atoms_solvt->nr);
    snew(atoms_solvt->resinfo, atoms_solvt->nr);
    snew(atoms_solvt->atomname, atoms_solvt->nr);
    snew(atoms_solvt->atom, atoms_solvt->nr);
    atoms_solvt->pdbinfo = NULL;
    fprintf(stderr, "Reading solvent configuration%s\n",
            v_solvt ? " and velocities" : "");
    read_stx_conf(filename, title_solvt, atoms_solvt, x_solvt, v_solvt,
                  &ePBC_solvt, box_solvt);
    fprintf(stderr, "\"%s\"\n", title_solvt);
    fprintf(stderr, "solvent configuration contains %d atoms in %d residues\n",
            atoms_solvt->nr, atoms_solvt->nres);
    fprintf(stderr, "\n");

    /* apply pbc for solvent configuration for whole molecules */
    rm_res_pbc(atoms_solvt, x_solvt, box_solvt);

    /* initialise distance arrays for solvent configuration */
    exclusionDistances_solvt = makeExclusionDistances(atoms_solvt, aps, defaultDistance, scaleFactor);

    /* calculate the box multiplication factors n_box[0...DIM] */
    nmol = 1;
    for (i = 0; (i < DIM); i++)
    {
        n_box[i] = 1;
        while (n_box[i]*box_solvt[i][i] < box[i][i])
        {
            n_box[i]++;
        }
        nmol *= n_box[i];
    }
    fprintf(stderr, "Will generate new solvent configuration of %dx%dx%d boxes\n",
            n_box[XX], n_box[YY], n_box[ZZ]);

    /* realloc atoms_solvt for the new solvent configuration */
    srenew(atoms_solvt->resinfo, atoms_solvt->nres*nmol);
    srenew(atoms_solvt->atomname, atoms_solvt->nr*nmol);
    srenew(atoms_solvt->atom, atoms_solvt->nr*nmol);
    srenew(x_solvt, atoms_solvt->nr*nmol);
    if (v_solvt)
    {
        srenew(v_solvt, atoms_solvt->nr*nmol);
    }
    srenew(exclusionDistances_solvt, atoms_solvt->nr*nmol);

    /* generate a new solvent configuration */
    make_new_conformation(atoms_solvt, x_solvt, v_solvt, exclusionDistances_solvt, box_solvt, n_box);

#ifdef DEBUG
    print_stat(x_solvt, atoms_solvt->nr, box_solvt);
#endif

#ifdef DEBUG
    print_stat(x_solvt, atoms_solvt->nr, box_solvt);
#endif
    /* Sort the solvent mixture, not the protein... */
    sort_molecule(&atoms_solvt, x_solvt, v_solvt, exclusionDistances_solvt);

    /* add the two configurations */
    onr   = atoms->nr;
    onres = atoms->nres;
    add_conf(atoms, x, v, exclusionDistances, TRUE, ePBC, box, FALSE,
             atoms_solvt, x_solvt, v_solvt, exclusionDistances_solvt, TRUE, rshell, max_sol, oenv);
    *atoms_added    = atoms->nr-onr;
    *residues_added = atoms->nres-onres;

    sfree(x_solvt);
    sfree(exclusionDistances_solvt);
    done_atom(atoms_solvt);
    sfree(atoms_solvt);

    fprintf(stderr, "Generated solvent containing %d atoms in %d residues\n",
            *atoms_added, *residues_added);
}

static void update_top(t_atoms *atoms, matrix box, int NFILE, t_filenm fnm[],
                       gmx_atomprop_t aps)
{
#define TEMP_FILENM "temp.top"
    FILE       *fpin, *fpout;
    char        buf[STRLEN], buf2[STRLEN], *temp;
    const char *topinout;
    int         line;
    gmx_bool    bSystem, bMolecules, bSkip;
    int         i, nsol = 0;
    double      mtot;
    real        vol, mm;

    for (i = 0; (i < atoms->nres); i++)
    {
        /* calculate number of SOLvent molecules */
        if ( (strcmp(*atoms->resinfo[i].name, "SOL") == 0) ||
             (strcmp(*atoms->resinfo[i].name, "WAT") == 0) ||
             (strcmp(*atoms->resinfo[i].name, "HOH") == 0) )
        {
            nsol++;
        }
    }
    mtot = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        gmx_atomprop_query(aps, epropMass,
                           *atoms->resinfo[atoms->atom[i].resind].name,
                           *atoms->atomname[i], &mm);
        mtot += mm;
    }

    vol = det(box);

    fprintf(stderr, "Volume                 :  %10g (nm^3)\n", vol);
    fprintf(stderr, "Density                :  %10g (g/l)\n",
            (mtot*1e24)/(AVOGADRO*vol));
    fprintf(stderr, "Number of SOL molecules:  %5d   \n\n", nsol);

    /* open topology file and append sol molecules */
    topinout  = ftp2fn(efTOP, NFILE, fnm);
    if (ftp2bSet(efTOP, NFILE, fnm) )
    {
        fprintf(stderr, "Processing topology\n");
        fpin    = gmx_ffopen(topinout, "r");
        fpout   = gmx_ffopen(TEMP_FILENM, "w");
        line    = 0;
        bSystem = bMolecules = FALSE;
        while (fgets(buf, STRLEN, fpin))
        {
            bSkip = FALSE;
            line++;
            strcpy(buf2, buf);
            if ((temp = strchr(buf2, '\n')) != NULL)
            {
                temp[0] = '\0';
            }
            ltrim(buf2);
            if (buf2[0] == '[')
            {
                buf2[0] = ' ';
                if ((temp = strchr(buf2, '\n')) != NULL)
                {
                    temp[0] = '\0';
                }
                rtrim(buf2);
                if (buf2[strlen(buf2)-1] == ']')
                {
                    buf2[strlen(buf2)-1] = '\0';
                    ltrim(buf2);
                    rtrim(buf2);
                    bSystem    = (gmx_strcasecmp(buf2, "system") == 0);
                    bMolecules = (gmx_strcasecmp(buf2, "molecules") == 0);
                }
            }
            else if (bSystem && nsol && (buf[0] != ';') )
            {
                /* if sol present, append "in water" to system name */
                rtrim(buf2);
                if (buf2[0] && (!strstr(buf2, " water")) )
                {
                    sprintf(buf, "%s in water\n", buf2);
                    bSystem = FALSE;
                }
            }
            else if (bMolecules)
            {
                /* check if this is a line with solvent molecules */
                sscanf(buf, "%4095s", buf2);
                if (strcmp(buf2, "SOL") == 0)
                {
                    sscanf(buf, "%*4095s %20d", &i);
                    nsol -= i;
                    if (nsol < 0)
                    {
                        bSkip = TRUE;
                        nsol += i;
                    }
                }
            }
            if (bSkip)
            {
                if ((temp = strchr(buf, '\n')) != NULL)
                {
                    temp[0] = '\0';
                }
                fprintf(stdout, "Removing line #%d '%s' from topology file (%s)\n",
                        line, buf, topinout);
            }
            else
            {
                fprintf(fpout, "%s", buf);
            }
        }
        gmx_ffclose(fpin);
        if (nsol)
        {
            fprintf(stdout, "Adding line for %d solvent molecules to "
                    "topology file (%s)\n", nsol, topinout);
            fprintf(fpout, "%-15s %5d\n", "SOL", nsol);
        }
        gmx_ffclose(fpout);
        /* use gmx_ffopen to generate backup of topinout */
        fpout = gmx_ffopen(topinout, "w");
        gmx_ffclose(fpout);
        rename(TEMP_FILENM, topinout);
    }
#undef TEMP_FILENM
}

int gmx_solvate(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] can do one of 2 things:[PAR]",

        "1) Generate a box of solvent. Specify [TT]-cs[tt] and [TT]-box[tt].",
        "Or specify [TT]-cs[tt] and [TT]-cp[tt] with a structure file with",
        "a box, but without atoms.[PAR]",

        "2) Solvate a solute configuration, e.g. a protein, in a bath of solvent ",
        "molecules. Specify [TT]-cp[tt] (solute) and [TT]-cs[tt] (solvent). ",
        "The box specified in the solute coordinate file ([TT]-cp[tt]) is used,",
        "unless [TT]-box[tt] is set.",
        "If you want the solute to be centered in the box,",
        "the program [gmx-editconf] has sophisticated options",
        "to change the box dimensions and center the solute.",
        "Solvent molecules are removed from the box where the ",
        "distance between any atom of the solute molecule(s) and any atom of ",
        "the solvent molecule is less than the sum of the scaled van der Waals",
        "radii of both atoms. A database ([TT]vdwradii.dat[tt]) of van der",
        "Waals radii is read by the program, and the resulting radii scaled",
        "by [TT]-scale[tt]. If radii are not found in the database, those"
        "atoms are assigned the (pre-scaled) distance [TT]-radius[tt].[PAR]",

        "The default solvent is Simple Point Charge water (SPC), with coordinates ",
        "from [TT]$GMXLIB/spc216.gro[tt]. These coordinates can also be used",
        "for other 3-site water models, since a short equibilibration will remove",
        "the small differences between the models.",
        "Other solvents are also supported, as well as mixed solvents. The",
        "only restriction to solvent types is that a solvent molecule consists",
        "of exactly one residue. The residue information in the coordinate",
        "files is used, and should therefore be more or less consistent.",
        "In practice this means that two subsequent solvent molecules in the ",
        "solvent coordinate file should have different residue number.",
        "The box of solute is built by stacking the coordinates read from",
        "the coordinate file. This means that these coordinates should be ",
        "equlibrated in periodic boundary conditions to ensure a good",
        "alignment of molecules on the stacking interfaces.",
        "The [TT]-maxsol[tt] option simply adds only the first [TT]-maxsol[tt]",
        "solvent molecules and leaves out the rest that would have fitted",
        "into the box. This can create a void that can cause problems later.",
        "Choose your volume wisely.[PAR]",

        "The program can optionally rotate the solute molecule to align the",
        "longest molecule axis along a box edge. This way the amount of solvent",
        "molecules necessary is reduced.",
        "It should be kept in mind that this only works for",
        "short simulations, as e.g. an alpha-helical peptide in solution can ",
        "rotate over 90 degrees, within 500 ps. In general it is therefore ",
        "better to make a more or less cubic box.[PAR]",

        "Setting [TT]-shell[tt] larger than zero will place a layer of water of",
        "the specified thickness (nm) around the solute. Hint: it is a good",
        "idea to put the protein in the center of a box first (using [gmx-editconf]).",
        "[PAR]",

        "Finally, [THISMODULE] will optionally remove lines from your topology file in ",
        "which a number of solvent molecules is already added, and adds a ",
        "line with the total number of solvent molecules in your coordinate file."
    };

    const char *bugs[] = {
        "Molecules must be whole in the initial configurations.",
    };

    /* parameter data */
    gmx_bool       bProt, bBox;
    const char    *conf_prot, *confout;
    real          *exclusionDistances = NULL;
    gmx_atomprop_t aps;

    /* protein configuration data */
    char    *title = NULL;
    t_atoms *atoms;
    rvec    *x    = NULL, *v = NULL;
    int      ePBC = -1;
    matrix   box;

    /* other data types */
    int      atoms_added, residues_added;

    t_filenm fnm[] = {
        { efSTX, "-cp", "protein", ffOPTRD },
        { efSTX, "-cs", "spc216",  ffLIBRD},
        { efSTO, NULL,  NULL,      ffWRITE},
        { efTOP, NULL,  NULL,      ffOPTRW},
    };
#define NFILE asize(fnm)

    static real     defaultDistance = 0.105, r_shell = 0, scaleFactor = 0.57;
    static rvec     new_box         = {0.0, 0.0, 0.0};
    static gmx_bool bReadV          = FALSE;
    static int      max_sol         = 0;
    output_env_t    oenv;
    t_pargs         pa[]              = {
        { "-box",    FALSE, etRVEC, {new_box},
          "Box size (in nm)" },
        { "-radius",   FALSE, etREAL, {&defaultDistance},
          "Default van der Waals distance"},
        { "-scale", FALSE, etREAL, {&scaleFactor},
          "Scale factor to multiply Van der Waals radii from the database in share/gromacs/top/vdwradii.dat. The default value of 0.57 yields density close to 1000 g/l for proteins in water." },
        { "-shell",  FALSE, etREAL, {&r_shell},
          "Thickness of optional water layer around solute" },
        { "-maxsol", FALSE, etINT,  {&max_sol},
          "Maximum number of solvent molecules to add if they fit in the box. If zero (default) this is ignored" },
        { "-vel",    FALSE, etBOOL, {&bReadV},
          "Keep velocities from input solute and solvent" },
    };

    if (!parse_common_args(&argc, argv, PCA_BE_NICE, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    const char *solventFileName = opt2fn("-cs", NFILE, fnm);
    bProt     = opt2bSet("-cp", NFILE, fnm);
    bBox      = opt2parg_bSet("-box", asize(pa), pa);

    /* check input */
    if (!bProt && !bBox)
    {
        gmx_fatal(FARGS, "When no solute (-cp) is specified, "
                  "a box size (-box) must be specified");
    }

    aps = gmx_atomprop_init();

    snew(atoms, 1);
    init_t_atoms(atoms, 0, FALSE);
    if (bProt)
    {
        /* Generate a solute configuration */
        conf_prot = opt2fn("-cp", NFILE, fnm);
        title     = readConformation(conf_prot, atoms, &x,
                                     bReadV ? &v : NULL, &ePBC, box);
        exclusionDistances = makeExclusionDistances(atoms, aps, defaultDistance, scaleFactor);

        if (bReadV && !v)
        {
            fprintf(stderr, "Note: no velocities found\n");
        }
        if (atoms->nr == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", conf_prot);
            bProt = FALSE;
        }
    }
    if (bBox)
    {
        ePBC = epbcXYZ;
        clear_mat(box);
        box[XX][XX] = new_box[XX];
        box[YY][YY] = new_box[YY];
        box[ZZ][ZZ] = new_box[ZZ];
    }
    if (det(box) == 0)
    {
        gmx_fatal(FARGS, "Undefined solute box.\nCreate one with gmx editconf "
                  "or give explicit -box command line option");
    }

    add_solv(solventFileName, atoms, &x, v ? &v : NULL, &exclusionDistances, ePBC, box,
             aps, defaultDistance, scaleFactor, &atoms_added, &residues_added, r_shell, max_sol,
             oenv);

    /* write new configuration 1 to file confout */
    confout = ftp2fn(efSTO, NFILE, fnm);
    fprintf(stderr, "Writing generated configuration to %s\n", confout);
    if (bProt)
    {
        write_sto_conf(confout, title, atoms, x, v, ePBC, box);
        /* print box sizes and box type to stderr */
        fprintf(stderr, "%s\n", title);
    }
    else
    {
        write_sto_conf(confout, "Generated by gmx solvate", atoms, x, v, ePBC, box);
    }

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n",
            atoms->nr, atoms->nres);
    update_top(atoms, box, NFILE, fnm, aps);

    gmx_atomprop_destroy(aps);
    sfree(exclusionDistances);
    sfree(x);
    sfree(v);
    done_atom(atoms);
    sfree(atoms);

    return 0;
}
