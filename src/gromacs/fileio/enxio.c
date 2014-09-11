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
#include "gmxpre.h"

#include "enxio.h"

#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/* The source code in this file should be thread-safe.
         Please keep it that way. */

/* This number should be increased whenever the file format changes! */
static const int enx_version = 5;

const char      *enx_block_id_name[] = {
    "Averaged orientation restraints",
    "Instantaneous orientation restraints",
    "Orientation restraint order tensor(s)",
    "Distance restraints",
    "Free energy data",
    "BAR histogram",
    "Delta H raw data"
};


/* Stuff for reading pre 4.1 energy files */
typedef struct {
    gmx_bool     bOldFileOpen;   /* Is this an open old file? */
    gmx_bool     bReadFirstStep; /* Did we read the first step? */
    int          first_step;     /* First step in the energy file */
    int          step_prev;      /* Previous step */
    int          nsum_prev;      /* Previous step sum length */
    t_energy    *ener_prev;      /* Previous energy sums */
} ener_old_t;

struct ener_file
{
    ener_old_t eo;
    t_fileio  *fio;
    int        framenr;
    real       frametime;
};

static void enxsubblock_init(t_enxsubblock *sb)
{
    sb->nr = 0;
#ifdef GMX_DOUBLE
    sb->type = xdr_datatype_double;
#else
    sb->type = xdr_datatype_float;
#endif
    sb->fval       = NULL;
    sb->dval       = NULL;
    sb->ival       = NULL;
    sb->lval       = NULL;
    sb->cval       = NULL;
    sb->sval       = NULL;
    sb->fval_alloc = 0;
    sb->dval_alloc = 0;
    sb->ival_alloc = 0;
    sb->lval_alloc = 0;
    sb->cval_alloc = 0;
    sb->sval_alloc = 0;
}

static void enxsubblock_free(t_enxsubblock *sb)
{
    if (sb->fval_alloc)
    {
        sfree(sb->fval);
        sb->fval_alloc = 0;
        sb->fval       = NULL;
    }
    if (sb->dval_alloc)
    {
        sfree(sb->dval);
        sb->dval_alloc = 0;
        sb->dval       = NULL;
    }
    if (sb->ival_alloc)
    {
        sfree(sb->ival);
        sb->ival_alloc = 0;
        sb->ival       = NULL;
    }
    if (sb->lval_alloc)
    {
        sfree(sb->lval);
        sb->lval_alloc = 0;
        sb->lval       = NULL;
    }
    if (sb->cval_alloc)
    {
        sfree(sb->cval);
        sb->cval_alloc = 0;
        sb->cval       = NULL;
    }
    if (sb->sval_alloc)
    {
        int i;

        for (i = 0; i < sb->sval_alloc; i++)
        {
            if (sb->sval[i])
            {
                sfree(sb->sval[i]);
            }
        }
        sfree(sb->sval);
        sb->sval_alloc = 0;
        sb->sval       = NULL;
    }
}

/* allocate the appropriate amount of memory for the given type and nr */
static void enxsubblock_alloc(t_enxsubblock *sb)
{
    /* allocate the appropriate amount of memory */
    switch (sb->type)
    {
        case xdr_datatype_float:
            if (sb->nr > sb->fval_alloc)
            {
                srenew(sb->fval, sb->nr);
                sb->fval_alloc = sb->nr;
            }
            break;
        case xdr_datatype_double:
            if (sb->nr > sb->dval_alloc)
            {
                srenew(sb->dval, sb->nr);
                sb->dval_alloc = sb->nr;
            }
            break;
        case xdr_datatype_int:
            if (sb->nr > sb->ival_alloc)
            {
                srenew(sb->ival, sb->nr);
                sb->ival_alloc = sb->nr;
            }
            break;
        case xdr_datatype_int64:
            if (sb->nr > sb->lval_alloc)
            {
                srenew(sb->lval, sb->nr);
                sb->lval_alloc = sb->nr;
            }
            break;
        case xdr_datatype_char:
            if (sb->nr > sb->cval_alloc)
            {
                srenew(sb->cval, sb->nr);
                sb->cval_alloc = sb->nr;
            }
            break;
        case xdr_datatype_string:
            if (sb->nr > sb->sval_alloc)
            {
                int i;

                srenew(sb->sval, sb->nr);
                for (i = sb->sval_alloc; i < sb->nr; i++)
                {
                    sb->sval[i] = NULL;
                }
                sb->sval_alloc = sb->nr;
            }
            break;
        default:
            gmx_incons("Unknown block type: this file is corrupted or from the future");
    }
}

static void enxblock_init(t_enxblock *eb)
{
    eb->id         = enxOR;
    eb->nsub       = 0;
    eb->sub        = NULL;
    eb->nsub_alloc = 0;
}

static void enxblock_free(t_enxblock *eb)
{
    if (eb->nsub_alloc > 0)
    {
        int i;
        for (i = 0; i < eb->nsub_alloc; i++)
        {
            enxsubblock_free(&(eb->sub[i]));
        }
        sfree(eb->sub);
        eb->nsub_alloc = 0;
        eb->sub        = NULL;
    }
}

void init_enxframe(t_enxframe *fr)
{
    fr->e_alloc = 0;
    fr->ener    = NULL;

    /*fr->d_alloc=0;*/
    fr->ener = NULL;

    /*fr->ndisre=0;*/

    fr->nblock       = 0;
    fr->nblock_alloc = 0;
    fr->block        = NULL;
}


void free_enxframe(t_enxframe *fr)
{
    int b;

    if (fr->e_alloc)
    {
        sfree(fr->ener);
    }
    for (b = 0; b < fr->nblock_alloc; b++)
    {
        enxblock_free(&(fr->block[b]));
    }
    sfree(fr->block);
}

void add_blocks_enxframe(t_enxframe *fr, int n)
{
    fr->nblock = n;
    if (n > fr->nblock_alloc)
    {
        int b;

        srenew(fr->block, n);
        for (b = fr->nblock_alloc; b < fr->nblock; b++)
        {
            enxblock_init(&(fr->block[b]));
        }
        fr->nblock_alloc = n;
    }
}

t_enxblock *find_block_id_enxframe(t_enxframe *ef, int id, t_enxblock *prev)
{
    gmx_off_t starti = 0;
    gmx_off_t i;

    if (prev)
    {
        starti = (prev - ef->block) + 1;
    }
    for (i = starti; i < ef->nblock; i++)
    {
        if (ef->block[i].id == id)
        {
            return &(ef->block[i]);
        }
    }
    return NULL;
}

void add_subblocks_enxblock(t_enxblock *eb, int n)
{
    eb->nsub = n;
    if (eb->nsub > eb->nsub_alloc)
    {
        int b;

        srenew(eb->sub, n);
        for (b = eb->nsub_alloc; b < n; b++)
        {
            enxsubblock_init(&(eb->sub[b]));
        }
        eb->nsub_alloc = n;
    }
}

static void enx_warning(const char *msg)
{
    if (getenv("GMX_ENX_NO_FATAL") != NULL)
    {
        gmx_warning(msg);
    }
    else
    {
        gmx_fatal(FARGS, "%s\n%s",
                  msg,
                  "If you want to use the correct frames before the corrupted frame and avoid this fatal error set the env.var. GMX_ENX_NO_FATAL");
    }
}

static void edr_strings(XDR *xdr, gmx_bool bRead, int file_version,
                        int n, gmx_enxnm_t **nms)
{
    int          i;
    gmx_enxnm_t *nm;

    if (*nms == NULL)
    {
        snew(*nms, n);
    }
    for (i = 0; i < n; i++)
    {
        nm = &(*nms)[i];
        if (bRead)
        {
            if (nm->name)
            {
                sfree(nm->name);
                nm->name = NULL;
            }
            if (nm->unit)
            {
                sfree(nm->unit);
                nm->unit = NULL;
            }
        }
        if (!xdr_string(xdr, &(nm->name), STRLEN))
        {
            gmx_file("Cannot write energy names to file; maybe you are out of disk space?");
        }
        if (file_version >= 2)
        {
            if (!xdr_string(xdr, &(nm->unit), STRLEN))
            {
                gmx_file("Cannot write energy names to file; maybe you are out of disk space?");
            }
        }
        else
        {
            nm->unit = gmx_strdup("kJ/mol");
        }
    }
}

void do_enxnms(ener_file_t ef, int *nre, gmx_enxnm_t **nms)
{
    int      magic = -55555;
    XDR     *xdr;
    gmx_bool bRead = gmx_fio_getread(ef->fio);
    int      file_version;
    int      i;

    gmx_fio_checktype(ef->fio);

    xdr = gmx_fio_getxdr(ef->fio);

    if (!xdr_int(xdr, &magic))
    {
        if (!bRead)
        {
            gmx_file("Cannot write energy names to file; maybe you are out of disk space?");
        }
        *nre = 0;
        return;
    }
    if (magic > 0)
    {
        /* Assume this is an old edr format */
        file_version          = 1;
        *nre                  = magic;
        ef->eo.bOldFileOpen   = TRUE;
        ef->eo.bReadFirstStep = FALSE;
        srenew(ef->eo.ener_prev, *nre);
    }
    else
    {
        ef->eo.bOldFileOpen = FALSE;

        if (magic != -55555)
        {
            gmx_fatal(FARGS, "Energy names magic number mismatch, this is not a GROMACS edr file");
        }
        file_version = enx_version;
        xdr_int(xdr, &file_version);
        if (file_version > enx_version)
        {
            gmx_fatal(FARGS, "reading tpx file (%s) version %d with version %d program", gmx_fio_getname(ef->fio), file_version, enx_version);
        }
        xdr_int(xdr, nre);
    }
    if (file_version != enx_version)
    {
        fprintf(stderr, "Note: enx file_version %d, software version %d\n",
                file_version, enx_version);
    }

    edr_strings(xdr, bRead, file_version, *nre, nms);
}

static gmx_bool do_eheader(ener_file_t ef, int *file_version, t_enxframe *fr,
                           int nre_test, gmx_bool *bWrongPrecision, gmx_bool *bOK)
{
    int          magic = -7777777;
    real         first_real_to_check;
    int          b, i, zero = 0, dum = 0;
    gmx_bool     bRead      = gmx_fio_getread(ef->fio);
    int          tempfix_nr = 0;
    int          ndisre     = 0;
    int          startb     = 0;
#ifndef GMX_DOUBLE
    xdr_datatype dtreal = xdr_datatype_float;
#else
    xdr_datatype dtreal = xdr_datatype_double;
#endif

    if (bWrongPrecision)
    {
        *bWrongPrecision = FALSE;
    }

    *bOK = TRUE;
    /* The original energy frame started with a real,
     * so we have to use a real for compatibility.
     * This is VERY DIRTY code, since do_eheader can be called
     * with the wrong precision set and then we could read r > -1e10,
     * while actually the intention was r < -1e10.
     * When nre_test >= 0, do_eheader should therefore terminate
     * before the number of i/o calls starts depending on what has been read
     * (which is the case for for instance the block sizes for variable
     * number of blocks, where this number is read before).
     */
    first_real_to_check = -2e10;
    if (!gmx_fio_do_real(ef->fio, first_real_to_check))
    {
        return FALSE;
    }
    if (first_real_to_check > -1e10)
    {
        /* Assume we are reading an old format */
        *file_version = 1;
        fr->t         = first_real_to_check;
        if (!gmx_fio_do_int(ef->fio, dum))
        {
            *bOK = FALSE;
        }
        fr->step = dum;
    }
    else
    {
        if (!gmx_fio_do_int(ef->fio, magic))
        {
            *bOK = FALSE;
        }
        if (magic != -7777777)
        {
            enx_warning("Energy header magic number mismatch, this is not a GROMACS edr file");
            *bOK = FALSE;
            return FALSE;
        }
        *file_version = enx_version;
        if (!gmx_fio_do_int(ef->fio, *file_version))
        {
            *bOK = FALSE;
        }
        if (*bOK && *file_version > enx_version)
        {
            gmx_fatal(FARGS, "reading tpx file (%s) version %d with version %d program", gmx_fio_getname(ef->fio), file_version, enx_version);
        }
        if (!gmx_fio_do_double(ef->fio, fr->t))
        {
            *bOK = FALSE;
        }
        if (!gmx_fio_do_int64(ef->fio, fr->step))
        {
            *bOK = FALSE;
        }
        if (!bRead && fr->nsum == 1)
        {
            /* Do not store sums of length 1,
             * since this does not add information.
             */
            if (!gmx_fio_do_int(ef->fio, zero))
            {
                *bOK = FALSE;
            }
        }
        else
        {
            if (!gmx_fio_do_int(ef->fio, fr->nsum))
            {
                *bOK = FALSE;
            }
        }
        if (*file_version >= 3)
        {
            if (!gmx_fio_do_int64(ef->fio, fr->nsteps))
            {
                *bOK = FALSE;
            }
        }
        else
        {
            fr->nsteps = max(1, fr->nsum);
        }
        if (*file_version >= 5)
        {
            if (!gmx_fio_do_double(ef->fio, fr->dt))
            {
                *bOK = FALSE;
            }
        }
        else
        {
            fr->dt = 0;
        }
    }
    if (!gmx_fio_do_int(ef->fio, fr->nre))
    {
        *bOK = FALSE;
    }
    if (*file_version < 4)
    {
        if (!gmx_fio_do_int(ef->fio, ndisre))
        {
            *bOK = FALSE;
        }
    }
    else
    {
        /* now reserved for possible future use */
        if (!gmx_fio_do_int(ef->fio, dum))
        {
            *bOK = FALSE;
        }
    }

    if (!gmx_fio_do_int(ef->fio, fr->nblock))
    {
        *bOK = FALSE;
    }
    if (fr->nblock < 0)
    {
        *bOK = FALSE;
    }

    if (ndisre != 0)
    {
        if (*file_version >= 4)
        {
            enx_warning("Distance restraint blocks in old style in new style file");
            *bOK = FALSE;
            return FALSE;
        }
        fr->nblock += 1;
    }


    /* Frames could have nre=0, so we can not rely only on the fr->nre check */
    if (bRead && nre_test >= 0 &&
        ((fr->nre > 0 && fr->nre != nre_test) ||
         fr->nre < 0 || ndisre < 0 || fr->nblock < 0))
    {
        *bWrongPrecision = TRUE;
        return *bOK;
    }

    /* we now know what these should be, or we've already bailed out because
       of wrong precision */
    if (*file_version == 1 && (fr->t < 0 || fr->t > 1e20 || fr->step < 0 ) )
    {
        enx_warning("edr file with negative step number or unreasonable time (and without version number).");
        *bOK = FALSE;
        return FALSE;
    }


    if (*bOK && bRead)
    {
        add_blocks_enxframe(fr, fr->nblock);
    }

    startb = 0;
    if (ndisre > 0)
    {
        /* sub[0] is the instantaneous data, sub[1] is time averaged */
        add_subblocks_enxblock(&(fr->block[0]), 2);
        fr->block[0].id          = enxDISRE;
        fr->block[0].sub[0].nr   = ndisre;
        fr->block[0].sub[1].nr   = ndisre;
        fr->block[0].sub[0].type = dtreal;
        fr->block[0].sub[1].type = dtreal;
        startb++;
    }

    /* read block header info */
    for (b = startb; b < fr->nblock; b++)
    {
        if (*file_version < 4)
        {
            /* blocks in old version files always have 1 subblock that
               consists of reals. */
            int nrint;

            if (bRead)
            {
                add_subblocks_enxblock(&(fr->block[b]), 1);
            }
            else
            {
                if (fr->block[b].nsub != 1)
                {
                    gmx_incons("Writing an old version .edr file with too many subblocks");
                }
                if (fr->block[b].sub[0].type != dtreal)
                {
                    gmx_incons("Writing an old version .edr file the wrong subblock type");
                }
            }
            nrint = fr->block[b].sub[0].nr;

            if (!gmx_fio_do_int(ef->fio, nrint))
            {
                *bOK = FALSE;
            }
            fr->block[b].id          = b - startb;
            fr->block[b].sub[0].nr   = nrint;
            fr->block[b].sub[0].type = dtreal;
        }
        else
        {
            int i;
            /* in the new version files, the block header only contains
               the ID and the number of subblocks */
            int nsub = fr->block[b].nsub;
            *bOK = *bOK && gmx_fio_do_int(ef->fio, fr->block[b].id);
            *bOK = *bOK && gmx_fio_do_int(ef->fio, nsub);

            fr->block[b].nsub = nsub;
            if (bRead)
            {
                add_subblocks_enxblock(&(fr->block[b]), nsub);
            }

            /* read/write type & size for each subblock */
            for (i = 0; i < nsub; i++)
            {
                t_enxsubblock *sub    = &(fr->block[b].sub[i]); /* shortcut */
                int            typenr = sub->type;

                *bOK = *bOK && gmx_fio_do_int(ef->fio, typenr);
                *bOK = *bOK && gmx_fio_do_int(ef->fio, sub->nr);

                sub->type = (xdr_datatype)typenr;
            }
        }
    }
    if (!gmx_fio_do_int(ef->fio, fr->e_size))
    {
        *bOK = FALSE;
    }

    /* now reserved for possible future use */
    if (!gmx_fio_do_int(ef->fio, dum))
    {
        *bOK = FALSE;
    }

    /* Do a dummy int to keep the format compatible with the old code */
    if (!gmx_fio_do_int(ef->fio, dum))
    {
        *bOK = FALSE;
    }

    if (*bOK && *file_version == 1 && nre_test < 0)
    {
        if (!ef->eo.bReadFirstStep)
        {
            ef->eo.bReadFirstStep = TRUE;
            ef->eo.first_step     = fr->step;
            ef->eo.step_prev      = fr->step;
            ef->eo.nsum_prev      = 0;
        }

        fr->nsum   = fr->step - ef->eo.first_step + 1;
        fr->nsteps = fr->step - ef->eo.step_prev;
        fr->dt     = 0;
    }

    return *bOK;
}

void free_enxnms(int n, gmx_enxnm_t *nms)
{
    int i;

    for (i = 0; i < n; i++)
    {
        sfree(nms[i].name);
        sfree(nms[i].unit);
    }

    sfree(nms);
}

void close_enx(ener_file_t ef)
{
    if (gmx_fio_close(ef->fio) != 0)
    {
        gmx_file("Cannot close energy file; it might be corrupt, or maybe you are out of disk space?");
    }
}

static gmx_bool empty_file(const char *fn)
{
    FILE    *fp;
    char     dum;
    int      ret;
    gmx_bool bEmpty;

    fp     = gmx_fio_fopen(fn, "r");
    ret    = fread(&dum, sizeof(dum), 1, fp);
    bEmpty = feof(fp);
    gmx_fio_fclose(fp);

    return bEmpty;
}


ener_file_t open_enx(const char *fn, const char *mode)
{
    int               nre, i;
    gmx_enxnm_t      *nms          = NULL;
    int               file_version = -1;
    t_enxframe       *fr;
    gmx_bool          bWrongPrecision, bOK = TRUE;
    struct ener_file *ef;

    snew(ef, 1);

    if (mode[0] == 'r')
    {
        ef->fio = gmx_fio_open(fn, mode);
        gmx_fio_checktype(ef->fio);
        gmx_fio_setprecision(ef->fio, FALSE);
        do_enxnms(ef, &nre, &nms);
        snew(fr, 1);
        do_eheader(ef, &file_version, fr, nre, &bWrongPrecision, &bOK);
        if (!bOK)
        {
            gmx_file("Cannot read energy file header. Corrupt file?");
        }

        /* Now check whether this file is in single precision */
        if (!bWrongPrecision &&
            ((fr->e_size && (fr->nre == nre) &&
              (nre*4*(long int)sizeof(float) == fr->e_size)) ) )
        {
            fprintf(stderr, "Opened %s as single precision energy file\n", fn);
            free_enxnms(nre, nms);
        }
        else
        {
            gmx_fio_rewind(ef->fio);
            gmx_fio_checktype(ef->fio);
            gmx_fio_setprecision(ef->fio, TRUE);
            do_enxnms(ef, &nre, &nms);
            do_eheader(ef, &file_version, fr, nre, &bWrongPrecision, &bOK);
            if (!bOK)
            {
                gmx_file("Cannot write energy file header; maybe you are out of disk space?");
            }

            if (((fr->e_size && (fr->nre == nre) &&
                  (nre*4*(long int)sizeof(double) == fr->e_size)) ))
            {
                fprintf(stderr, "Opened %s as double precision energy file\n",
                        fn);
            }
            else
            {
                if (empty_file(fn))
                {
                    gmx_fatal(FARGS, "File %s is empty", fn);
                }
                else
                {
                    gmx_fatal(FARGS, "Energy file %s not recognized, maybe different CPU?",
                              fn);
                }
            }
            free_enxnms(nre, nms);
        }
        free_enxframe(fr);
        sfree(fr);
        gmx_fio_rewind(ef->fio);
    }
    else
    {
        ef->fio = gmx_fio_open(fn, mode);
    }

    ef->framenr   = 0;
    ef->frametime = 0;
    return ef;
}

t_fileio *enx_file_pointer(const ener_file_t ef)
{
    return ef->fio;
}

static void convert_full_sums(ener_old_t *ener_old, t_enxframe *fr)
{
    int    nstep_all;
    int    ne, ns, i;
    double esum_all, eav_all;

    if (fr->nsum > 0)
    {
        ne = 0;
        ns = 0;
        for (i = 0; i < fr->nre; i++)
        {
            if (fr->ener[i].e    != 0)
            {
                ne++;
            }
            if (fr->ener[i].esum != 0)
            {
                ns++;
            }
        }
        if (ne > 0 && ns == 0)
        {
            /* We do not have all energy sums */
            fr->nsum = 0;
        }
    }

    /* Convert old full simulation sums to sums between energy frames */
    nstep_all = fr->step - ener_old->first_step + 1;
    if (fr->nsum > 1 && fr->nsum == nstep_all && ener_old->nsum_prev > 0)
    {
        /* Set the new sum length: the frame step difference */
        fr->nsum = fr->step - ener_old->step_prev;
        for (i = 0; i < fr->nre; i++)
        {
            esum_all         = fr->ener[i].esum;
            eav_all          = fr->ener[i].eav;
            fr->ener[i].esum = esum_all - ener_old->ener_prev[i].esum;
            fr->ener[i].eav  = eav_all  - ener_old->ener_prev[i].eav
                - dsqr(ener_old->ener_prev[i].esum/(nstep_all - fr->nsum)
                       - esum_all/nstep_all)*
                (nstep_all - fr->nsum)*nstep_all/(double)fr->nsum;
            ener_old->ener_prev[i].esum = esum_all;
            ener_old->ener_prev[i].eav  = eav_all;
        }
        ener_old->nsum_prev = nstep_all;
    }
    else if (fr->nsum > 0)
    {
        if (fr->nsum != nstep_all)
        {
            fprintf(stderr, "\nWARNING: something is wrong with the energy sums, will not use exact averages\n");
            ener_old->nsum_prev = 0;
        }
        else
        {
            ener_old->nsum_prev = nstep_all;
        }
        /* Copy all sums to ener_prev */
        for (i = 0; i < fr->nre; i++)
        {
            ener_old->ener_prev[i].esum = fr->ener[i].esum;
            ener_old->ener_prev[i].eav  = fr->ener[i].eav;
        }
    }

    ener_old->step_prev = fr->step;
}

gmx_bool do_enx(ener_file_t ef, t_enxframe *fr)
{
    int           file_version = -1;
    int           i, b;
    gmx_bool      bRead, bOK, bOK1, bSane;
    real          tmp1, tmp2, rdum;
    /*int       d_size;*/

    bOK   = TRUE;
    bRead = gmx_fio_getread(ef->fio);
    if (!bRead)
    {
        fr->e_size = fr->nre*sizeof(fr->ener[0].e)*4;
        /*d_size = fr->ndisre*(sizeof(real)*2);*/
    }
    gmx_fio_checktype(ef->fio);

    if (!do_eheader(ef, &file_version, fr, -1, NULL, &bOK))
    {
        if (bRead)
        {
            fprintf(stderr, "\rLast energy frame read %d time %8.3f         ",
                    ef->framenr-1, ef->frametime);
            if (!bOK)
            {
                fprintf(stderr,
                        "\nWARNING: Incomplete energy frame: nr %d time %8.3f\n",
                        ef->framenr, fr->t);
            }
        }
        else
        {
            gmx_file("Cannot write energy file header; maybe you are out of disk space?");
        }
        return FALSE;
    }
    if (bRead)
    {
        if ((ef->framenr <   20 || ef->framenr %   10 == 0) &&
            (ef->framenr <  200 || ef->framenr %  100 == 0) &&
            (ef->framenr < 2000 || ef->framenr % 1000 == 0))
        {
            fprintf(stderr, "\rReading energy frame %6d time %8.3f         ",
                    ef->framenr, fr->t);
        }
        ef->framenr++;
        ef->frametime = fr->t;
    }
    /* Check sanity of this header */
    bSane = fr->nre > 0;
    for (b = 0; b < fr->nblock; b++)
    {
        bSane = bSane || (fr->block[b].nsub > 0);
    }
    if (!((fr->step >= 0) && bSane))
    {
        fprintf(stderr, "\nWARNING: there may be something wrong with energy file %s\n",
                gmx_fio_getname(ef->fio));
        fprintf(stderr, "Found: step=%"GMX_PRId64 ", nre=%d, nblock=%d, time=%g.\n"
                "Trying to skip frame expect a crash though\n",
                fr->step, fr->nre, fr->nblock, fr->t);
    }
    if (bRead && fr->nre > fr->e_alloc)
    {
        srenew(fr->ener, fr->nre);
        for (i = fr->e_alloc; (i < fr->nre); i++)
        {
            fr->ener[i].e    = 0;
            fr->ener[i].eav  = 0;
            fr->ener[i].esum = 0;
        }
        fr->e_alloc = fr->nre;
    }

    for (i = 0; i < fr->nre; i++)
    {
        bOK = bOK && gmx_fio_do_real(ef->fio, fr->ener[i].e);

        /* Do not store sums of length 1,
         * since this does not add information.
         */
        if (file_version == 1 ||
            (bRead && fr->nsum > 0) || fr->nsum > 1)
        {
            tmp1 = fr->ener[i].eav;
            bOK  = bOK && gmx_fio_do_real(ef->fio, tmp1);
            if (bRead)
            {
                fr->ener[i].eav = tmp1;
            }

            /* This is to save only in single precision (unless compiled in DP) */
            tmp2 = fr->ener[i].esum;
            bOK  = bOK && gmx_fio_do_real(ef->fio, tmp2);
            if (bRead)
            {
                fr->ener[i].esum = tmp2;
            }

            if (file_version == 1)
            {
                /* Old, unused real */
                rdum = 0;
                bOK  = bOK && gmx_fio_do_real(ef->fio, rdum);
            }
        }
    }

    /* Here we can not check for file_version==1, since one could have
     * continued an old format simulation with a new one with mdrun -append.
     */
    if (bRead && ef->eo.bOldFileOpen)
    {
        /* Convert old full simulation sums to sums between energy frames */
        convert_full_sums(&(ef->eo), fr);
    }
    /* read the blocks */
    for (b = 0; b < fr->nblock; b++)
    {
        /* now read the subblocks. */
        int nsub = fr->block[b].nsub; /* shortcut */
        int i;

        for (i = 0; i < nsub; i++)
        {
            t_enxsubblock *sub = &(fr->block[b].sub[i]); /* shortcut */

            if (bRead)
            {
                enxsubblock_alloc(sub);
            }

            /* read/write data */
            bOK1 = TRUE;
            switch (sub->type)
            {
                case xdr_datatype_float:
                    bOK1 = gmx_fio_ndo_float(ef->fio, sub->fval, sub->nr);
                    break;
                case xdr_datatype_double:
                    bOK1 = gmx_fio_ndo_double(ef->fio, sub->dval, sub->nr);
                    break;
                case xdr_datatype_int:
                    bOK1 = gmx_fio_ndo_int(ef->fio, sub->ival, sub->nr);
                    break;
                case xdr_datatype_int64:
                    bOK1 = gmx_fio_ndo_int64(ef->fio, sub->lval, sub->nr);
                    break;
                case xdr_datatype_char:
                    bOK1 = gmx_fio_ndo_uchar(ef->fio, sub->cval, sub->nr);
                    break;
                case xdr_datatype_string:
                    bOK1 = gmx_fio_ndo_string(ef->fio, sub->sval, sub->nr);
                    break;
                default:
                    gmx_incons("Reading unknown block data type: this file is corrupted or from the future");
            }
            bOK = bOK && bOK1;
        }
    }

    if (!bRead)
    {
        if (gmx_fio_flush(ef->fio) != 0)
        {
            gmx_file("Cannot write energy file; maybe you are out of disk space?");
        }
    }

    if (!bOK)
    {
        if (bRead)
        {
            fprintf(stderr, "\nLast energy frame read %d",
                    ef->framenr-1);
            fprintf(stderr, "\nWARNING: Incomplete energy frame: nr %d time %8.3f\n",
                    ef->framenr, fr->t);
        }
        else
        {
            gmx_fatal(FARGS, "could not write energies");
        }
        return FALSE;
    }

    return TRUE;
}

static real find_energy(const char *name, int nre, gmx_enxnm_t *enm,
                        t_enxframe *fr)
{
    int i;

    for (i = 0; i < nre; i++)
    {
        if (strcmp(enm[i].name, name) == 0)
        {
            return fr->ener[i].e;
        }
    }

    gmx_fatal(FARGS, "Could not find energy term named '%s'", name);

    return 0;
}


void get_enx_state(const char *fn, real t, gmx_groups_t *groups, t_inputrec *ir,
                   t_state *state)
{
    /* Should match the names in mdebin.c */
    static const char *boxvel_nm[] = {
        "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
        "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY"
    };

    static const char *pcouplmu_nm[] = {
        "Pcoupl-Mu-XX", "Pcoupl-Mu-YY", "Pcoupl-Mu-ZZ",
        "Pcoupl-Mu-YX", "Pcoupl-Mu-ZX", "Pcoupl-Mu-ZY"
    };
    static const char *baro_nm[] = {
        "Barostat"
    };


    int          ind0[] = { XX, YY, ZZ, YY, ZZ, ZZ };
    int          ind1[] = { XX, YY, ZZ, XX, XX, YY };
    int          nre, nfr, i, j, ni, npcoupl;
    char         buf[STRLEN];
    const char  *bufi;
    gmx_enxnm_t *enm = NULL;
    t_enxframe  *fr;
    ener_file_t  in;

    in = open_enx(fn, "r");
    do_enxnms(in, &nre, &enm);
    snew(fr, 1);
    nfr = 0;
    while ((nfr == 0 || fr->t != t) && do_enx(in, fr))
    {
        nfr++;
    }
    close_enx(in);
    fprintf(stderr, "\n");

    if (nfr == 0 || fr->t != t)
    {
        gmx_fatal(FARGS, "Could not find frame with time %f in '%s'", t, fn);
    }

    npcoupl = TRICLINIC(ir->compress) ? 6 : 3;
    if (ir->epc == epcPARRINELLORAHMAN)
    {
        clear_mat(state->boxv);
        for (i = 0; i < npcoupl; i++)
        {
            state->boxv[ind0[i]][ind1[i]] =
                find_energy(boxvel_nm[i], nre, enm, fr);
        }
        fprintf(stderr, "\nREAD %d BOX VELOCITIES FROM %s\n\n", npcoupl, fn);
    }

    if (ir->etc == etcNOSEHOOVER)
    {
        char cns[20];

        cns[0] = '\0';

        for (i = 0; i < state->ngtc; i++)
        {
            ni   = groups->grps[egcTC].nm_ind[i];
            bufi = *(groups->grpname[ni]);
            for (j = 0; (j < state->nhchainlength); j++)
            {
                if (IR_NVT_TROTTER(ir))
                {
                    sprintf(cns, "-%d", j);
                }
                sprintf(buf, "Xi%s-%s", cns, bufi);
                state->nosehoover_xi[i] = find_energy(buf, nre, enm, fr);
                sprintf(buf, "vXi%s-%s", cns, bufi);
                state->nosehoover_vxi[i] = find_energy(buf, nre, enm, fr);
            }

        }
        fprintf(stderr, "\nREAD %d NOSE-HOOVER Xi chains FROM %s\n\n", state->ngtc, fn);

        if (IR_NPT_TROTTER(ir) || IR_NPH_TROTTER(ir))
        {
            for (i = 0; i < state->nnhpres; i++)
            {
                bufi = baro_nm[0]; /* All barostat DOF's together for now */
                for (j = 0; (j < state->nhchainlength); j++)
                {
                    sprintf(buf, "Xi-%d-%s", j, bufi);
                    state->nhpres_xi[i] = find_energy(buf, nre, enm, fr);
                    sprintf(buf, "vXi-%d-%s", j, bufi);
                    state->nhpres_vxi[i] = find_energy(buf, nre, enm, fr);
                }
            }
            fprintf(stderr, "\nREAD %d NOSE-HOOVER BAROSTAT Xi chains FROM %s\n\n", state->nnhpres, fn);
        }
    }

    free_enxnms(nre, enm);
    free_enxframe(fr);
    sfree(fr);
}
