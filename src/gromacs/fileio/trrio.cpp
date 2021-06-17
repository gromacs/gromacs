/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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

#include "trrio.h"

#include <cstring>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static int nFloatSize(gmx_trr_header_t* sh)
{
    int nflsize = 0;

    if (sh->box_size)
    {
        nflsize = sh->box_size / (DIM * DIM);
    }
    else if (sh->x_size)
    {
        nflsize = sh->x_size / (sh->natoms * DIM);
    }
    else if (sh->v_size)
    {
        nflsize = sh->v_size / (sh->natoms * DIM);
    }
    else if (sh->f_size)
    {
        nflsize = sh->f_size / (sh->natoms * DIM);
    }
    else
    {
        gmx_file("Can not determine precision of trr file");
    }

    if (((nflsize != sizeof(float)) && (nflsize != sizeof(double))))
    {
        gmx_fatal(FARGS, "Float size %d. Maybe different CPU?", nflsize);
    }

    return nflsize;
}

/* Returns whether a valid frame header was read. Upon exit, *bOK is
   TRUE if a normal outcome resulted. Usually that is the same thing,
   but encountering the end of the file before reading the magic
   integer is a normal outcome for TRR reading, and it does not
   produce a valid frame header, so the values differ in that case.
   That does not exclude the possibility of a reading error between
   frames, but the trajectory-handling infrastructure needs an
   overhaul before we can handle that. */
static gmx_bool do_trr_frame_header(t_fileio* fio, bool bRead, gmx_trr_header_t* sh, gmx_bool* bOK)
{
    const int       magicValue = 1993;
    int             magic      = magicValue;
    static gmx_bool bFirst     = TRUE;
    char            buf[256];

    *bOK = TRUE;

    if (!gmx_fio_do_int(fio, magic))
    {
        /* Failed to read an integer, which should be the magic
           number, which usually means we've reached the end
           of the file (but could be an I/O error that we now
           might mishandle). */
        return FALSE;
    }
    if (magic != magicValue)
    {
        *bOK = FALSE;
        gmx_fatal(FARGS,
                  "Failed to find GROMACS magic number in trr frame header, so this is not a trr "
                  "file!\n");
    }

    if (bRead)
    {
        *bOK = *bOK && gmx_fio_do_string(fio, buf);
        if (bFirst)
        {
            fprintf(stderr, "trr version: %s ", buf);
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

    /* Note that TRR wasn't defined to be extensible, so we can't fix
     * the fact that we used a default int for the step number, which
     * is typically defined to be signed and 32 bit. */
    int intStep = sh->step;
    *bOK        = *bOK && gmx_fio_do_int(fio, intStep);
    sh->step    = intStep;
    *bOK        = *bOK && gmx_fio_do_int(fio, sh->nre);
    *bOK        = *bOK && gmx_fio_do_real(fio, sh->t);
    *bOK        = *bOK && gmx_fio_do_real(fio, sh->lambda);

    return *bOK;
}

static gmx_bool do_trr_frame_data(t_fileio* fio, gmx_trr_header_t* sh, rvec* box, rvec* x, rvec* v, rvec* f)
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
    if (sh->x_size != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, x, sh->natoms);
    }
    if (sh->v_size != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, v, sh->natoms);
    }
    if (sh->f_size != 0)
    {
        bOK = bOK && gmx_fio_ndo_rvec(fio, f, sh->natoms);
    }

    return bOK;
}

static gmx_bool do_trr_frame(t_fileio* fio,
                             bool      bRead,
                             int64_t*  step,
                             real*     t,
                             real*     lambda,
                             rvec*     box,
                             int*      natoms,
                             rvec*     x,
                             rvec*     v,
                             rvec*     f)
{
    gmx_trr_header_t* sh;
    gmx_bool          bOK;

    snew(sh, 1);
    if (!bRead)
    {
        sh->box_size = (box) ? sizeof(matrix) : 0;
        sh->x_size   = ((x) ? (*natoms * sizeof(x[0])) : 0);
        sh->v_size   = ((v) ? (*natoms * sizeof(v[0])) : 0);
        sh->f_size   = ((f) ? (*natoms * sizeof(f[0])) : 0);
        sh->natoms   = *natoms;
        sh->step     = *step;
        sh->nre      = 0;
        sh->t        = *t;
        sh->lambda   = *lambda;
    }
    if (!do_trr_frame_header(fio, bRead, sh, &bOK))
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
            gmx_file("inputrec in trr file");
        }
        if (sh->e_size)
        {
            gmx_file("energies in trr file");
        }
        if (sh->top_size)
        {
            gmx_file("topology in trr file");
        }
        if (sh->sym_size)
        {
            gmx_file("symbol table in trr file");
        }
    }
    bOK = do_trr_frame_data(fio, sh, box, x, v, f);

    sfree(sh);

    return bOK;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/

void gmx_trr_read_single_header(const char* fn, gmx_trr_header_t* header)
{
    t_fileio* fio = gmx_trr_open(fn, "r");
    gmx_bool  bOK;
    if (!do_trr_frame_header(fio, true, header, &bOK))
    {
        gmx_fatal(FARGS, "Empty file %s", fn);
    }
    gmx_trr_close(fio);
}

gmx_bool gmx_trr_read_frame_header(t_fileio* fio, gmx_trr_header_t* header, gmx_bool* bOK)
{
    return do_trr_frame_header(fio, true, header, bOK);
}

void gmx_trr_write_single_frame(const char* fn,
                                int64_t     step,
                                real        t,
                                real        lambda,
                                const rvec* box,
                                int         natoms,
                                const rvec* x,
                                const rvec* v,
                                const rvec* f)
{
    t_fileio* fio = gmx_trr_open(fn, "w");
    do_trr_frame(fio,
                 false,
                 &step,
                 &t,
                 &lambda,
                 const_cast<rvec*>(box),
                 &natoms,
                 const_cast<rvec*>(x),
                 const_cast<rvec*>(v),
                 const_cast<rvec*>(f));
    gmx_trr_close(fio);
}

void gmx_trr_read_single_frame(const char* fn,
                               int64_t*    step,
                               real*       t,
                               real*       lambda,
                               rvec*       box,
                               int*        natoms,
                               rvec*       x,
                               rvec*       v,
                               rvec*       f)
{
    t_fileio* fio = gmx_trr_open(fn, "r");
    do_trr_frame(fio, true, step, t, lambda, box, natoms, x, v, f);
    gmx_trr_close(fio);
}

void gmx_trr_write_frame(t_fileio*   fio,
                         int64_t     step,
                         real        t,
                         real        lambda,
                         const rvec* box,
                         int         natoms,
                         const rvec* x,
                         const rvec* v,
                         const rvec* f)
{
    if (!do_trr_frame(fio,
                      false,
                      &step,
                      &t,
                      &lambda,
                      const_cast<rvec*>(box),
                      &natoms,
                      const_cast<rvec*>(x),
                      const_cast<rvec*>(v),
                      const_cast<rvec*>(f)))
    {
        gmx_file("Cannot write trajectory frame; maybe you are out of disk space?");
    }
}


gmx_bool gmx_trr_read_frame(t_fileio* fio,
                            int64_t*  step,
                            real*     t,
                            real*     lambda,
                            rvec*     box,
                            int*      natoms,
                            rvec*     x,
                            rvec*     v,
                            rvec*     f)
{
    return do_trr_frame(fio, true, step, t, lambda, box, natoms, x, v, f);
}

gmx_bool gmx_trr_read_frame_data(t_fileio* fio, gmx_trr_header_t* header, rvec* box, rvec* x, rvec* v, rvec* f)
{
    return do_trr_frame_data(fio, header, box, x, v, f);
}

t_fileio* gmx_trr_open(const char* fn, const char* mode)
{
    return gmx_fio_open(fn, mode);
}

void gmx_trr_close(t_fileio* fio)
{
    gmx_fio_close(fio);
}
