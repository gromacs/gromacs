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

#include "trnio.h"

#include <string.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#define BUFSIZE     128
#define GROMACS_MAGIC   1993

static int nFloatSize(t_trnheader *sh)
{
    int nflsize = 0;

    if (sh->box_size)
    {
        nflsize = sh->box_size/(DIM*DIM);
    }
    else if (sh->x_size)
    {
        nflsize = sh->x_size/(sh->natoms*DIM);
    }
    else if (sh->v_size)
    {
        nflsize = sh->v_size/(sh->natoms*DIM);
    }
    else if (sh->f_size)
    {
        nflsize = sh->f_size/(sh->natoms*DIM);
    }
    else
    {
        gmx_file("Can not determine precision of trn file");
    }

    if (((nflsize != sizeof(float)) && (nflsize != sizeof(double))))
    {
        gmx_fatal(FARGS, "Float size %d. Maybe different CPU?", nflsize);
    }

    return nflsize;
}

static gmx_bool do_trnheader(t_fileio *fio, gmx_bool bRead, t_trnheader *sh, gmx_bool *bOK)
{
    int             magic  = GROMACS_MAGIC;
    static gmx_bool bFirst = TRUE;
    char            buf[256];

    *bOK = TRUE;

    gmx_fio_checktype(fio);

    if (!gmx_fio_do_int(fio, magic) || magic != GROMACS_MAGIC)
    {
        return FALSE;
    }

    if (bRead)
    {
        *bOK = *bOK && gmx_fio_do_string(fio, buf);
        if (bFirst)
        {
            fprintf(stderr, "trn version: %s ", buf);
        }
    }
    else
    {
        sprintf(buf, "GMX_trn_file");
        *bOK = *bOK && gmx_fio_do_string(fio, buf);
    }
    *bOK = *bOK && gmx_fio_do_int(fio, sh->ir_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->e_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->box_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->vir_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->pres_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->top_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->sym_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->x_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->v_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->f_size);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->natoms);

    if (!*bOK)
    {
        return *bOK;
    }
    sh->bDouble = (nFloatSize(sh) == sizeof(double));
    gmx_fio_setprecision(fio, sh->bDouble);

    if (bRead && bFirst)
    {
        fprintf(stderr, "(%s precision)\n", sh->bDouble ? "double" : "single");
        bFirst = FALSE;
    }

    *bOK = *bOK && gmx_fio_do_int(fio, sh->step);
    *bOK = *bOK && gmx_fio_do_int(fio, sh->nre);
    *bOK = *bOK && gmx_fio_do_real(fio, sh->t);
    *bOK = *bOK && gmx_fio_do_real(fio, sh->lambda);

    return *bOK;
}

void pr_trnheader(FILE *fp, int indent, char *title, t_trnheader *sh)
{
    if (sh)
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "box_size    = %d\n", sh->box_size);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "x_size      = %d\n", sh->x_size);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "v_size      = %d\n", sh->v_size);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "f_size      = %d\n", sh->f_size);

        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "natoms      = %d\n", sh->natoms);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "step        = %d\n", sh->step);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "t           = %e\n", sh->t);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "lambda      = %e\n", sh->lambda);
    }
}

static gmx_bool do_htrn(t_fileio *fio, t_trnheader *sh,
                        rvec *box, rvec *x, rvec *v, rvec *f)
{
    matrix   pv;
    gmx_bool bOK;

    bOK = TRUE;
    if (sh->box_size != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, box, DIM);
    }
    if (sh->vir_size != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, pv, DIM);
    }
    if (sh->pres_size != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, pv, DIM);
    }
    if (sh->x_size   != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, x, sh->natoms);
    }
    if (sh->v_size   != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, v, sh->natoms);
    }
    if (sh->f_size   != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, f, sh->natoms);
    }

    return bOK;
}

static gmx_bool do_trn(t_fileio *fio, gmx_bool bRead, int *step, real *t, real *lambda,
                       rvec *box, int *natoms, rvec *x, rvec *v, rvec *f)
{
    t_trnheader *sh;
    gmx_bool     bOK;

    snew(sh, 1);
    if (!bRead)
    {
        sh->box_size = (box) ? sizeof(matrix) : 0;
        sh->x_size   = ((x) ? (*natoms*sizeof(x[0])) : 0);
        sh->v_size   = ((v) ? (*natoms*sizeof(v[0])) : 0);
        sh->f_size   = ((f) ? (*natoms*sizeof(f[0])) : 0);
        sh->natoms   = *natoms;
        sh->step     = *step;
        sh->nre      = 0;
        sh->t        = *t;
        sh->lambda   = *lambda;
    }
    if (!do_trnheader(fio, bRead, sh, &bOK))
    {
        return FALSE;
    }
    if (bRead)
    {
        *natoms = sh->natoms;
        *step   = sh->step;
        *t      = sh->t;
        *lambda = sh->lambda;
        if (sh->ir_size)
        {
            gmx_file("inputrec in trn file");
        }
        if (sh->e_size)
        {
            gmx_file("energies in trn file");
        }
        if (sh->top_size)
        {
            gmx_file("topology in trn file");
        }
        if (sh->sym_size)
        {
            gmx_file("symbol table in trn file");
        }
    }
    bOK = do_htrn(fio, sh, box, x, v, f);

    sfree(sh);

    return bOK;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/

void read_trnheader(const char *fn, t_trnheader *trn)
{
    t_fileio *fio;
    gmx_bool  bOK;

    fio = open_trn(fn, "r");
    if (!do_trnheader(fio, TRUE, trn, &bOK))
    {
        gmx_fatal(FARGS, "Empty file %s", fn);
    }
    close_trn(fio);
}

gmx_bool fread_trnheader(t_fileio *fio, t_trnheader *trn, gmx_bool *bOK)
{
    return do_trnheader(fio, TRUE, trn, bOK);
}

void write_trn(const char *fn, int step, real t, real lambda,
               rvec *box, int natoms, rvec *x, rvec *v, rvec *f)
{
    t_fileio *fio;

    fio = open_trn(fn, "w");
    do_trn(fio, FALSE, &step, &t, &lambda, box, &natoms, x, v, f);
    close_trn(fio);
}

void read_trn(const char *fn, int *step, real *t, real *lambda,
              rvec *box, int *natoms, rvec *x, rvec *v, rvec *f)
{
    t_fileio *fio;

    fio = open_trn(fn, "r");
    (void) do_trn(fio, TRUE, step, t, lambda, box, natoms, x, v, f);
    close_trn(fio);
}

void fwrite_trn(t_fileio *fio, int step, real t, real lambda,
                rvec *box, int natoms, rvec *x, rvec *v, rvec *f)
{
    if (do_trn(fio, FALSE, &step, &t, &lambda, box, &natoms, x, v, f) == FALSE)
    {
        gmx_file("Cannot write trajectory frame; maybe you are out of disk space?");
    }
}


gmx_bool fread_trn(t_fileio *fio, int *step, real *t, real *lambda,
                   rvec *box, int *natoms, rvec *x, rvec *v, rvec *f)
{
    return do_trn(fio, TRUE, step, t, lambda, box, natoms, x, v, f);
}

gmx_bool fread_htrn(t_fileio *fio, t_trnheader *trn, rvec *box, rvec *x, rvec *v,
                    rvec *f)
{
    return do_htrn(fio, trn, box, x, v, f);
}

t_fileio *open_trn(const char *fn, const char *mode)
{
    return gmx_fio_open(fn, mode);
}

void close_trn(t_fileio *fio)
{
    gmx_fio_close(fio);
}
