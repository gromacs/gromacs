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
#include "filenm.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "futil.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "types/simple.h"
#include "types/commrec.h"

/* XDR should be available on all platforms now,
 * but we keep the possibility of turning it off...
 */
#define USE_XDR

/* Use bitflag ... */
#define IS_SET(fn) ((fn.flag & ffSET) != 0)
#define IS_OPT(fn) ((fn.flag & ffOPT) != 0)
#define IS_MULT(fn) ((fn.flag & ffMULT) != 0)
#define UN_SET(fn) (fn.flag = (fn.flag & ~ffSET))
#define DO_SET(fn) (fn.flag = (fn.flag |  ffSET))

enum
{
    eftASC, eftBIN, eftXDR, eftTNG, eftGEN, eftNR
};

/* To support multiple file types with one general (eg TRX) we have
 * these arrays.
 */
static const int trxs[] =
{
#ifdef USE_XDR
    efXTC, efTRR, efCPT,
#endif
    efTRJ, efGRO, efG96, efPDB, efTNG
};
#define NTRXS asize(trxs)

static const int trcompressed[] =
{
#ifdef USE_XDR
    efXTC,
#endif
    efTNG
};
#define NTRCOMPRESSED asize(trcompressed)

static const int tros[] =
{
#ifdef USE_XDR
    efXTC, efTRR,
#endif
    efTRJ, efGRO, efG96, efPDB, efTNG
};
#define NTROS asize(tros)

static const int trns[] =
{
#ifdef USE_XDR
    efTRR, efCPT,
#endif
    efTRJ, efTNG
};
#define NTRNS asize(trns)

static const int stos[] =
{ efGRO, efG96, efPDB, efBRK, efENT, efESP };
#define NSTOS asize(stos)

static const int stxs[] =
{
    efGRO, efG96, efPDB, efBRK, efENT, efESP,
#ifdef USE_XDR
    efTPR,
#endif
    efTPB, efTPA
};
#define NSTXS asize(stxs)

static const int tpxs[] =
{
#ifdef USE_XDR
    efTPR,
#endif
    efTPB, efTPA
};
#define NTPXS asize(tpxs)

static const int tpss[] =
{
#ifdef USE_XDR
    efTPR,
#endif
    efTPB, efTPA, efGRO, efG96, efPDB, efBRK, efENT
};
#define NTPSS asize(tpss)

typedef struct
{
    int         ftype;
    const char *ext;
    const char *defnm;
    const char *defopt;
    const char *descr;
    int         ntps;
    const int  *tps;
} t_deffile;

/* this array should correspond to the enum in filenm.h */
static const t_deffile
    deffile[efNR] =
{
    { eftASC, ".mdp", "grompp", "-f", "grompp input file with MD parameters" },
    { eftGEN, ".???", "traj", "-f", "Trajectory", NTRXS, trxs },
    { eftGEN, ".???", "trajout", "-f", "Trajectory", NTROS, tros },
    { eftGEN, ".???", "traj", NULL,
      "Full precision trajectory", NTRNS, trns },
    { eftXDR, ".trr", "traj", NULL, "Trajectory in portable xdr format" },
    { eftBIN, ".trj", "traj", NULL, "Trajectory file (architecture specific)" },
    { eftGEN, ".???", "traj_comp", NULL,
      "Compressed trajectory (tng format or portable xdr format)", NTRCOMPRESSED, trcompressed},
    { eftXDR, ".xtc", "traj", NULL,
      "Compressed trajectory (portable xdr format): xtc" },
    { eftTNG, ".tng", "traj", NULL,
      "Trajectory file (tng format)" },
    { eftXDR, ".edr", "ener",   NULL, "Energy file"},
    { eftGEN, ".???", "conf", "-c", "Structure file", NSTXS, stxs },
    { eftGEN, ".???", "out", "-o", "Structure file", NSTOS, stos },
    { eftASC, ".gro", "conf", "-c", "Coordinate file in Gromos-87 format" },
    { eftASC, ".g96", "conf", "-c", "Coordinate file in Gromos-96 format" },
    { eftASC, ".pdb", "eiwit",  "-f", "Protein data bank file"},
    { eftASC, ".brk", "eiwit",  "-f", "Brookhaven data bank file"},
    { eftASC, ".ent", "eiwit", "-f", "Entry in the protein date bank" },
    { eftASC, ".esp", "conf", "-f", "Coordinate file in Espresso format" },
    { eftASC, ".pqr", "state",  "-o", "Coordinate file for MEAD"},
    { eftXDR, ".cpt", "state",  "-cp", "Checkpoint file"},
    { eftASC, ".log", "run",    "-l", "Log file"},
    { eftASC, ".xvg", "graph",  "-o", "xvgr/xmgr file"},
    { eftASC, ".out", "hello",  "-o", "Generic output file"},
    { eftASC, ".ndx", "index",  "-n", "Index file", },
    { eftASC, ".top", "topol",  "-p", "Topology file"},
    { eftASC, ".itp", "topinc", NULL, "Include file for topology"},
    { eftGEN, ".???", "topol", "-s", "Run input file", NTPXS, tpxs },
    { eftGEN, ".???", "topol", "-s", "Structure+mass(db)", NTPSS, tpss },
    { eftXDR, ".tpr", "topol",  "-s", "Portable xdr run input file"},
    { eftASC, ".tpa", "topol",  "-s", "Ascii run input file"},
    { eftBIN, ".tpb", "topol",  "-s", "Binary run input file"},
    { eftASC, ".tex", "doc",    "-o", "LaTeX file"},
    { eftASC, ".rtp", "residue", NULL, "Residue Type file used by pdb2gmx" },
    { eftASC, ".atp", "atomtp", NULL, "Atomtype file used by pdb2gmx" },
    { eftASC, ".hdb", "polar",  NULL, "Hydrogen data base"},
    { eftASC, ".dat", "nnnice", NULL, "Generic data file"},
    { eftASC, ".dlg", "user",   NULL, "Dialog Box data for ngmx"},
    { eftASC, ".map", "ss", NULL, "File that maps matrix data to colors" },
    { eftASC, ".eps", "plot", NULL, "Encapsulated PostScript (tm) file" },
    { eftASC, ".mat", "ss",     NULL, "Matrix Data file"},
    { eftASC, ".m2p", "ps",     NULL, "Input file for mat2ps"},
    { eftXDR, ".mtx", "hessian", "-m", "Hessian matrix"},
    { eftASC, ".edi", "sam",    NULL, "ED sampling input"},
    { eftASC, ".cub", "pot",  NULL, "Gaussian cube file" },
    { eftASC, ".xpm", "root", NULL, "X PixMap compatible matrix file" },
    { eftASC, "", "rundir", NULL, "Run directory" }
};

#define NZEXT 2
static const char *z_ext[NZEXT] =
{ ".gz", ".Z" };

const char *ftp2ext(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].ext[0] != '\0' ? deffile[ftp].ext + 1 : "";
    }
    else
    {
        return "unknown";
    }
}

const char *ftp2ext_generic(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        switch (ftp)
        {
            case efTRX:
                return "trx";
            case efTRN:
                return "trn";
            case efSTO:
                return "sto";
            case efSTX:
                return "stx";
            case efTPX:
                return "tpx";
            case efTPS:
                return "tps";
            default:
                return ftp2ext(ftp);
        }
    }
    else
    {
        return "unknown";
    }
}

const char *ftp2ext_with_dot(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].ext;
    }
    else
    {
        return "unknown";
    }
}

int ftp2generic_count(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].ntps;
    }
    else
    {
        return 0;
    }
}

const int *ftp2generic_list(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].tps;
    }
    else
    {
        return 0;
    }
}

const char *ftp2desc(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].descr;
    }
    else
    {
        return "unknown filetype";
    }
}

const char *ftp2ftype(int ftp)
{
    if ((ftp >= 0) && (ftp < efNR))
    {
        switch (deffile[ftp].ftype)
        {
            case eftASC:
                return "ASCII";
            case eftBIN:
                return "Binary";
            case eftXDR:
                return "XDR portable";
            case eftTNG:
                return "TNG";
            case eftGEN:
                return "";
            default:
                gmx_fatal(FARGS, "Unknown filetype %d in ftp2ftype", deffile[ftp].ftype);
                break;
        }
    }
    return "unknown";
}

const char *ftp2defnm(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].defnm;
    }
    else
    {
        return NULL;
    }
}

static void check_opts(int nf, t_filenm fnm[])
{
    int              i;
    const t_deffile *df;

    for (i = 0; (i < nf); i++)
    {
        df = &(deffile[fnm[i].ftp]);
        if (fnm[i].opt == NULL)
        {
            if (df->defopt == NULL)
            {
                gmx_fatal(FARGS, "No default cmd-line option for %s (type %d)\n",
                          deffile[fnm[i].ftp].ext, fnm[i].ftp);
            }
            else
            {
                fnm[i].opt = df->defopt;
            }
        }
    }
}

int fn2ftp(const char *fn)
{
    int         i, len;
    const char *feptr;
    const char *eptr;

    if (!fn)
    {
        return efNR;
    }

    len = strlen(fn);
    if ((len >= 4) && (fn[len - 4] == '.'))
    {
        feptr = &(fn[len - 4]);
    }
    else
    {
        return efNR;
    }

    for (i = 0; (i < efNR); i++)
    {
        if ((eptr = deffile[i].ext) != NULL)
        {
            if (gmx_strcasecmp(feptr, eptr) == 0)
            {
                break;
            }
        }
    }

    return i;
}

static void set_extension(char *buf, int ftp)
{
    int              len, extlen;
    const t_deffile *df;

    /* check if extension is already at end of filename */
    df     = &(deffile[ftp]);
    len    = strlen(buf);
    extlen = strlen(df->ext);
    if ((len <= extlen) || (gmx_strcasecmp(&(buf[len - extlen]), df->ext) != 0))
    {
        strcat(buf, df->ext);
    }
}

static void add_filenm(t_filenm *fnm, const char *filenm)
{
    srenew(fnm->fns, fnm->nfiles+1);
    fnm->fns[fnm->nfiles] = strdup(filenm);
    fnm->nfiles++;
}

static void set_grpfnm(t_filenm *fnm, const char *name, const char *deffnm)
{
    char       buf[256], buf2[256];
    int        i, type;
    gmx_bool   bValidExt;
    int        nopts;
    const int *ftps;

    nopts = deffile[fnm->ftp].ntps;
    ftps  = deffile[fnm->ftp].tps;
    if ((nopts == 0) || (ftps == NULL))
    {
        gmx_fatal(FARGS, "nopts == 0 || ftps == NULL");
    }

    bValidExt = FALSE;
    if (name && deffnm == NULL)
    {
        strcpy(buf, name);
        /* First check whether we have a valid filename already */
        type = fn2ftp(name);
        if ((fnm->flag & ffREAD) && (fnm->ftp == efTRX))
        {
            /*if file exist don't add an extension for trajectory reading*/
            bValidExt = gmx_fexist(name);
        }
        for (i = 0; (i < nopts) && !bValidExt; i++)
        {
            if (type == ftps[i])
            {
                bValidExt = TRUE;
            }
        }
    }
    else if (deffnm != NULL)
    {
        strcpy(buf, deffnm);
    }
    else
    {
        /* No name given, set the default name */
        strcpy(buf, ftp2defnm(fnm->ftp));
    }

    if (!bValidExt && (fnm->flag & ffREAD))
    {
        /* for input-files only: search for filenames in the directory */
        for (i = 0; (i < nopts) && !bValidExt; i++)
        {
            type = ftps[i];
            strcpy(buf2, buf);
            set_extension(buf2, type);
            if (gmx_fexist(buf2))
            {
                bValidExt = TRUE;
                strcpy(buf, buf2);
            }
        }
    }

    if (!bValidExt)
    {
        /* Use the first extension type */
        set_extension(buf, ftps[0]);
    }

    add_filenm(fnm, buf);
}

static void set_filenm(t_filenm *fnm, const char *name, const char *deffnm,
                       gmx_bool bReadNode)
{
    /* Set the default filename, extension and option for those fields that
     * are not already set. An extension is added if not present, if fn = NULL
     * or empty, the default filename is given.
     */
    char buf[256];
    int  i, len, extlen;

    if ((fnm->flag & ffREAD) && !bReadNode)
    {
        return;
    }

    if ((fnm->ftp < 0) || (fnm->ftp >= efNR))
    {
        gmx_fatal(FARGS, "file type out of range (%d)", fnm->ftp);
    }

    if (name)
    {
        strcpy(buf, name);
    }
    if ((fnm->flag & ffREAD) && name && gmx_fexist(name))
    {
        /* check if filename ends in .gz or .Z, if so remove that: */
        len = strlen(name);
        for (i = 0; i < NZEXT; i++)
        {
            extlen = strlen(z_ext[i]);
            if (len > extlen)
            {
                if (gmx_strcasecmp(name+len-extlen, z_ext[i]) == 0)
                {
                    buf[len-extlen] = '\0';
                    break;
                }
            }
        }
    }

    if (deffile[fnm->ftp].ntps)
    {
        set_grpfnm(fnm, name ? buf : NULL, deffnm);
    }
    else
    {
        if (name == NULL || deffnm != NULL)
        {
            if (deffnm != NULL)
            {
                strcpy(buf, deffnm);
            }
            else
            {
                strcpy(buf, ftp2defnm(fnm->ftp));
            }
        }
        set_extension(buf, fnm->ftp);

        add_filenm(fnm, buf);
    }
}

static void set_filenms(int nf, t_filenm fnm[], const char *deffnm, gmx_bool bReadNode)
{
    int i;

    for (i = 0; (i < nf); i++)
    {
        if (!IS_SET(fnm[i]))
        {
            set_filenm(&(fnm[i]), fnm[i].fn, deffnm, bReadNode);
        }
    }
}

void parse_file_args(int *argc, char *argv[], int nf, t_filenm fnm[],
                     const char *deffnm, gmx_bool bReadNode)
{
    int       i, j;
    gmx_bool *bRemove;

    check_opts(nf, fnm);

    for (i = 0; (i < nf); i++)
    {
        UN_SET(fnm[i]);
    }

    if (*argc > 1)
    {
        snew(bRemove, (*argc)+1);
        i = 1;
        do
        {
            for (j = 0; (j < nf); j++)
            {
                if (strcmp(argv[i], fnm[j].opt) == 0)
                {
                    DO_SET(fnm[j]);
                    bRemove[i] = TRUE;
                    i++;
                    /* check if we are out of arguments for this option */
                    if ((i >= *argc) || (argv[i][0] == '-'))
                    {
                        set_filenm(&fnm[j], fnm[j].fn, deffnm, bReadNode);
                    }
                    /* sweep up all file arguments for this option */
                    while ((i < *argc) && (argv[i][0] != '-'))
                    {
                        set_filenm(&fnm[j], argv[i], NULL, bReadNode);
                        bRemove[i] = TRUE;
                        i++;
                        /* only repeat for 'multiple' file options: */
                        if (!IS_MULT(fnm[j]))
                        {
                            break;
                        }
                    }

                    break; /* jump out of 'j' loop */
                }
            }
            /* No file found corresponding to option argv[i] */
            if (j == nf)
            {
                i++;
            }
        }
        while (i < *argc);

        /* Remove used entries */
        for (i = j = 0; (i <= *argc); i++)
        {
            if (!bRemove[i])
            {
                argv[j++] = argv[i];
            }
        }
        (*argc) = j - 1;
        sfree(bRemove);
    }

    set_filenms(nf, fnm, deffnm, bReadNode);

}

const char *opt2fn(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            return fnm[i].fns[0];
        }
    }

    fprintf(stderr, "No option %s\n", opt);

    return NULL;
}

const char *opt2fn_master(const char *opt, int nfile, const t_filenm fnm[],
                          t_commrec *cr)
{
    return SIMMASTER(cr) ? opt2fn(opt, nfile, fnm) : NULL;
}

int opt2fns(char **fns[], const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            *fns = fnm[i].fns;
            return fnm[i].nfiles;
        }
    }

    fprintf(stderr, "No option %s\n", opt);
    return 0;
}

const char *ftp2fn(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            return fnm[i].fns[0];
        }
    }

    fprintf(stderr, "ftp2fn: No filetype %s\n", deffile[ftp].ext);
    return NULL;
}

int ftp2fns(char **fns[], int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            *fns = fnm[i].fns;
            return fnm[i].nfiles;
        }
    }

    fprintf(stderr, "ftp2fn: No filetype %s\n", deffile[ftp].ext);
    return 0;
}

gmx_bool ftp2bSet(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            return (gmx_bool) IS_SET(fnm[i]);
        }
    }

    fprintf(stderr, "ftp2bSet: No filetype %s\n", deffile[ftp].ext);

    return FALSE;
}

gmx_bool opt2bSet(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            return (gmx_bool) IS_SET(fnm[i]);
        }
    }

    fprintf(stderr, "No option %s\n", opt);

    return FALSE;
}

const char *opt2fn_null(const char *opt, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (strcmp(opt, fnm[i].opt) == 0)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
            {
                return NULL;
            }
            else
            {
                return fnm[i].fns[0];
            }
        }
    }
    fprintf(stderr, "No option %s\n", opt);
    return NULL;
}

const char *ftp2fn_null(int ftp, int nfile, const t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (ftp == fnm[i].ftp)
        {
            if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
            {
                return NULL;
            }
            else
            {
                return fnm[i].fns[0];
            }
        }
    }
    fprintf(stderr, "ftp2fn: No filetype %s\n", deffile[ftp].ext);
    return NULL;
}

gmx_bool is_optional(const t_filenm *fnm)
{
    return ((fnm->flag & ffOPT) == ffOPT);
}

gmx_bool is_output(const t_filenm *fnm)
{
    return ((fnm->flag & ffWRITE) == ffWRITE);
}

gmx_bool is_set(const t_filenm *fnm)
{
    return ((fnm->flag & ffSET) == ffSET);
}

int add_suffix_to_output_names(t_filenm *fnm, int nfile, const char *suffix)
{
    int   i, j, pos;
    char  buf[STRLEN], newname[STRLEN];
    char *extpos;

    for (i = 0; i < nfile; i++)
    {
        if (is_output(&fnm[i]) && fnm[i].ftp != efCPT)
        {
            /* We never use multiple _outputs_, but we might as well check
               for it, just in case... */
            for (j = 0; j < fnm[i].nfiles; j++)
            {
                strncpy(buf, fnm[i].fns[j], STRLEN - 1);
                extpos  = strrchr(buf, '.');
                *extpos = '\0';
                sprintf(newname, "%s%s.%s", buf, suffix, extpos + 1);
                free(fnm[i].fns[j]);
                fnm[i].fns[j] = strdup(newname);
            }
        }
    }
    return 0;
}

t_filenm *dup_tfn(int nf, const t_filenm tfn[])
{
    int       i, j;
    t_filenm *ret;

    snew(ret, nf);
    for (i = 0; i < nf; i++)
    {
        ret[i] = tfn[i]; /* just directly copy all non-string fields */
        if (tfn[i].opt)
        {
            ret[i].opt = strdup(tfn[i].opt);
        }
        else
        {
            ret[i].opt = NULL;
        }

        if (tfn[i].fn)
        {
            ret[i].fn = strdup(tfn[i].fn);
        }
        else
        {
            ret[i].fn = NULL;
        }

        if (tfn[i].nfiles > 0)
        {
            snew(ret[i].fns, tfn[i].nfiles);
            for (j = 0; j < tfn[i].nfiles; j++)
            {
                ret[i].fns[j] = strdup(tfn[i].fns[j]);
            }
        }
    }
    return ret;
}

void done_filenms(int nf, t_filenm fnm[])
{
    int i, j;

    for (i = 0; i < nf; ++i)
    {
        for (j = 0; j < fnm[i].nfiles; ++j)
        {
            sfree(fnm[i].fns[j]);
        }
        sfree(fnm[i].fns);
        fnm[i].fns = NULL;
    }
}
