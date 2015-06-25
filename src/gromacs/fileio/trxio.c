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

#include "trxio.h"

#include "config.h"

#include <assert.h>
#include <math.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/timecontrol.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/tngio_for_tools.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#ifdef GMX_USE_PLUGINS
#include "gromacs/fileio/vmdio.h"
#endif

/* defines for frame counter output */
#define SKIP1   10
#define SKIP2  100
#define SKIP3 1000

struct t_trxstatus
{
    int                     __frame;
    t_trxframe             *xframe;
    int                     nxframe;
    t_fileio               *fio;
    tng_trajectory_t        tng;
    int                     natoms;
    double                  DT, BOX[3];
    gmx_bool                bReadBox;
    char                   *persistent_line; /* Persistent line for reading g96 trajectories */
};

/* utility functions */

gmx_bool bRmod_fd(double a, double b, double c, gmx_bool bDouble)
{
    int    iq;
    double tol;

    tol = 2*(bDouble ? GMX_DOUBLE_EPS : GMX_FLOAT_EPS);

    iq = (int)((a - b + tol*a)/c);

    if (fabs(a - b - c*iq) <= tol*fabs(a))
    {
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}



int check_times2(real t, real t0, gmx_bool bDouble)
{
    int  r;

#ifndef GMX_DOUBLE
    /* since t is float, we can not use double precision for bRmod */
    bDouble = FALSE;
#endif

    r = -1;
    if ((!bTimeSet(TBEGIN) || (t >= rTimeValue(TBEGIN)))  &&
        (!bTimeSet(TEND)   || (t <= rTimeValue(TEND))))
    {
        if (bTimeSet(TDELTA) && !bRmod_fd(t, t0, rTimeValue(TDELTA), bDouble))
        {
            r = -1;
        }
        else
        {
            r = 0;
        }
    }
    else if (bTimeSet(TEND) && (t >= rTimeValue(TEND)))
    {
        r = 1;
    }
    if (debug)
    {
        fprintf(debug, "t=%g, t0=%g, b=%g, e=%g, dt=%g: r=%d\n",
                t, t0, rTimeValue(TBEGIN), rTimeValue(TEND), rTimeValue(TDELTA), r);
    }
    return r;
}

int check_times(real t)
{
    return check_times2(t, t, FALSE);
}

static void initcount(t_trxstatus *status)
{
    status->__frame = -1;
}

static void status_init(t_trxstatus *status)
{
    status->nxframe         = 0;
    status->xframe          = NULL;
    status->fio             = NULL;
    status->__frame         = -1;
    status->persistent_line = NULL;
    status->tng             = NULL;
}


int nframes_read(t_trxstatus *status)
{
    return status->__frame;
}

static void printcount_(t_trxstatus *status, const output_env_t oenv,
                        const char *l, real t)
{
    if ((status->__frame < 2*SKIP1 || status->__frame % SKIP1 == 0) &&
        (status->__frame < 2*SKIP2 || status->__frame % SKIP2 == 0) &&
        (status->__frame < 2*SKIP3 || status->__frame % SKIP3 == 0))
    {
        fprintf(stderr, "\r%-14s %6d time %8.3f   ", l, status->__frame,
                output_env_conv_time(oenv, t));
    }
}

static void printcount(t_trxstatus *status, const output_env_t oenv, real t,
                       gmx_bool bSkip)
{
    status->__frame++;
    printcount_(status, oenv, bSkip ? "Skipping frame" : "Reading frame", t);
}

static void printlast(t_trxstatus *status, const output_env_t oenv, real t)
{
    printcount_(status, oenv, "Last frame", t);
    fprintf(stderr, "\n");
}

static void printincomp(t_trxstatus *status, t_trxframe *fr)
{
    if (fr->not_ok & HEADER_NOT_OK)
    {
        fprintf(stderr, "WARNING: Incomplete header: nr %d time %g\n",
                status->__frame+1, fr->time);
    }
    else if (fr->not_ok)
    {
        fprintf(stderr, "WARNING: Incomplete frame: nr %d time %g\n",
                status->__frame+1, fr->time);
    }
}

int prec2ndec(real prec)
{
    if (prec <= 0)
    {
        gmx_fatal(FARGS, "DEATH HORROR prec (%g) <= 0 in prec2ndec", prec);
    }

    return (int)(log(prec)/log(10.0)+0.5);
}

real ndec2prec(int ndec)
{
    return pow(10.0, ndec);
}

t_fileio *trx_get_fileio(t_trxstatus *status)
{
    return status->fio;
}

float trx_get_time_of_final_frame(t_trxstatus *status)
{
    t_fileio *stfio    = trx_get_fileio(status);
    int       filetype = gmx_fio_getftp(stfio);
    int       bOK;
    float     lasttime = -1;

    if (filetype == efXTC)
    {
        lasttime =
            xdr_xtc_get_last_frame_time(gmx_fio_getfp(stfio),
                                        gmx_fio_getxdr(stfio),
                                        status->natoms, &bOK);
        if (!bOK)
        {
            gmx_fatal(FARGS, "Error reading last frame. Maybe seek not supported." );
        }
    }
    else if (filetype == efTNG)
    {
        tng_trajectory_t tng = status->tng;
        if (!tng)
        {
            gmx_fatal(FARGS, "Error opening TNG file.");
        }
        lasttime = gmx_tng_get_time_of_final_frame(tng);
    }
    else
    {
        gmx_incons("Only supported for TNG and XTC");
    }
    return lasttime;
}

void clear_trxframe(t_trxframe *fr, gmx_bool bFirst)
{
    fr->not_ok    = 0;
    fr->bTitle    = FALSE;
    fr->bStep     = FALSE;
    fr->bTime     = FALSE;
    fr->bLambda   = FALSE;
    fr->bFepState = FALSE;
    fr->bAtoms    = FALSE;
    fr->bPrec     = FALSE;
    fr->bX        = FALSE;
    fr->bV        = FALSE;
    fr->bF        = FALSE;
    fr->bBox      = FALSE;
    if (bFirst)
    {
        fr->flags     = 0;
        fr->bDouble   = FALSE;
        fr->natoms    = -1;
        fr->t0        = 0;
        fr->tf        = 0;
        fr->tpf       = 0;
        fr->tppf      = 0;
        fr->title     = NULL;
        fr->step      = 0;
        fr->time      = 0;
        fr->lambda    = 0;
        fr->fep_state = 0;
        fr->atoms     = NULL;
        fr->prec      = 0;
        fr->x         = NULL;
        fr->v         = NULL;
        fr->f         = NULL;
        clear_mat(fr->box);
        fr->bPBC   = FALSE;
        fr->ePBC   = -1;
    }
}

void set_trxframe_ePBC(t_trxframe *fr, int ePBC)
{
    fr->bPBC = (ePBC == -1);
    fr->ePBC = ePBC;
}

int write_trxframe_indexed(t_trxstatus *status, t_trxframe *fr, int nind,
                           const atom_id *ind, gmx_conect gc)
{
    char  title[STRLEN];
    rvec *xout = NULL, *vout = NULL, *fout = NULL;
    int   i, ftp = -1;
    real  prec;

    if (fr->bPrec)
    {
        prec = fr->prec;
    }
    else
    {
        prec = 1000.0;
    }

    if (status->tng)
    {
        ftp = efTNG;
    }
    else if (status->fio)
    {
        ftp = gmx_fio_getftp(status->fio);
    }
    else
    {
        gmx_incons("No input file available");
    }

    switch (ftp)
    {
        case efTRR:
        case efTNG:
            break;
        default:
            if (!fr->bX)
            {
                gmx_fatal(FARGS, "Need coordinates to write a %s trajectory",
                          ftp2ext(ftp));
            }
            break;
    }

    switch (ftp)
    {
        case efTRR:
        case efTNG:
            if (fr->bV)
            {
                snew(vout, nind);
                for (i = 0; i < nind; i++)
                {
                    copy_rvec(fr->v[ind[i]], vout[i]);
                }
            }
            if (fr->bF)
            {
                snew(fout, nind);
                for (i = 0; i < nind; i++)
                {
                    copy_rvec(fr->f[ind[i]], fout[i]);
                }
            }
        /* no break */
        case efXTC:
            if (fr->bX)
            {
                snew(xout, nind);
                for (i = 0; i < nind; i++)
                {
                    copy_rvec(fr->x[ind[i]], xout[i]);
                }
            }
            break;
        default:
            break;
    }

    switch (ftp)
    {
        case efTNG:
            gmx_write_tng_from_trxframe(status->tng, fr, nind);
            break;
        case efXTC:
            write_xtc(status->fio, nind, fr->step, fr->time, fr->box, xout, prec);
            break;
        case efTRR:
            fwrite_trn(status->fio, nframes_read(status),
                       fr->time, fr->step, fr->box, nind, xout, vout, fout);
            break;
        case efGRO:
        case efPDB:
        case efBRK:
        case efENT:
            if (!fr->bAtoms)
            {
                gmx_fatal(FARGS, "Can not write a %s file without atom names",
                          ftp2ext(ftp));
            }
            sprintf(title, "frame t= %.3f", fr->time);
            if (ftp == efGRO)
            {
                write_hconf_indexed_p(gmx_fio_getfp(status->fio), title, fr->atoms, nind, ind,
                                      prec2ndec(prec),
                                      fr->x, fr->bV ? fr->v : NULL, fr->box);
            }
            else
            {
                write_pdbfile_indexed(gmx_fio_getfp(status->fio), title, fr->atoms,
                                      fr->x, -1, fr->box, ' ', fr->step, nind, ind, gc, TRUE);
            }
            break;
        case efG96:
            write_g96_conf(gmx_fio_getfp(status->fio), fr, nind, ind);
            break;
        default:
            gmx_fatal(FARGS, "Sorry, write_trxframe_indexed can not write %s",
                      ftp2ext(ftp));
            break;
    }

    switch (ftp)
    {
        case efTRR:
        case efTNG:
            if (vout)
            {
                sfree(vout);
            }
            if (fout)
            {
                sfree(fout);
            }
        /* no break */
        case efXTC:
            sfree(xout);
            break;
        default:
            break;
    }

    return 0;
}

void trjtools_gmx_prepare_tng_writing(const char       *filename,
                                      char              filemode,
                                      t_trxstatus      *in,
                                      t_trxstatus     **out,
                                      const char       *infile,
                                      const int         natoms,
                                      const gmx_mtop_t *mtop,
                                      const atom_id    *index,
                                      const char       *index_group_name)
{
    if (filemode != 'w' && filemode != 'a')
    {
        gmx_incons("Sorry, can only prepare for TNG output.");
    }

    if (*out == NULL)
    {
        snew((*out), 1);
    }
    status_init(*out);

    if (in != NULL)
    {
        gmx_prepare_tng_writing(filename,
                                filemode,
                                &in->tng,
                                &(*out)->tng,
                                natoms,
                                mtop,
                                index,
                                index_group_name);
    }
    else if (efTNG == fn2ftp(infile))
    {
        tng_trajectory_t tng_in;
        gmx_tng_open(infile, 'r', &tng_in);

        gmx_prepare_tng_writing(filename,
                                filemode,
                                &tng_in,
                                &(*out)->tng,
                                natoms,
                                mtop,
                                index,
                                index_group_name);
    }
}

void write_tng_frame(t_trxstatus *status,
                     t_trxframe  *frame)
{
    gmx_write_tng_from_trxframe(status->tng, frame, -1);
}

int write_trxframe(t_trxstatus *status, t_trxframe *fr, gmx_conect gc)
{
    char title[STRLEN];
    real prec;

    if (fr->bPrec)
    {
        prec = fr->prec;
    }
    else
    {
        prec = 1000.0;
    }

    if (status->tng)
    {
        gmx_tng_set_compression_precision(status->tng, prec);
        write_tng_frame(status, fr);

        return 0;
    }

    switch (gmx_fio_getftp(status->fio))
    {
        case efTRR:
            break;
        default:
            if (!fr->bX)
            {
                gmx_fatal(FARGS, "Need coordinates to write a %s trajectory",
                          ftp2ext(gmx_fio_getftp(status->fio)));
            }
            break;
    }

    switch (gmx_fio_getftp(status->fio))
    {
        case efXTC:
            write_xtc(status->fio, fr->natoms, fr->step, fr->time, fr->box, fr->x, prec);
            break;
        case efTRR:
            fwrite_trn(status->fio, fr->step, fr->time, fr->lambda, fr->box, fr->natoms,
                       fr->bX ? fr->x : NULL, fr->bV ? fr->v : NULL, fr->bF ? fr->f : NULL);
            break;
        case efGRO:
        case efPDB:
        case efBRK:
        case efENT:
            if (!fr->bAtoms)
            {
                gmx_fatal(FARGS, "Can not write a %s file without atom names",
                          ftp2ext(gmx_fio_getftp(status->fio)));
            }
            sprintf(title, "frame t= %.3f", fr->time);
            if (gmx_fio_getftp(status->fio) == efGRO)
            {
                write_hconf_p(gmx_fio_getfp(status->fio), title, fr->atoms,
                              prec2ndec(prec), fr->x, fr->bV ? fr->v : NULL, fr->box);
            }
            else
            {
                write_pdbfile(gmx_fio_getfp(status->fio), title,
                              fr->atoms, fr->x, fr->bPBC ? fr->ePBC : -1, fr->box,
                              ' ', fr->step, gc, TRUE);
            }
            break;
        case efG96:
            write_g96_conf(gmx_fio_getfp(status->fio), fr, -1, NULL);
            break;
        default:
            gmx_fatal(FARGS, "Sorry, write_trxframe can not write %s",
                      ftp2ext(gmx_fio_getftp(status->fio)));
            break;
    }

    return 0;
}

int write_trx(t_trxstatus *status, int nind, const atom_id *ind, t_atoms *atoms,
              int step, real time, matrix box, rvec x[], rvec *v,
              gmx_conect gc)
{
    t_trxframe fr;

    clear_trxframe(&fr, TRUE);
    fr.bStep  = TRUE;
    fr.step   = step;
    fr.bTime  = TRUE;
    fr.time   = time;
    fr.bAtoms = atoms != NULL;
    fr.atoms  = atoms;
    fr.bX     = TRUE;
    fr.x      = x;
    fr.bV     = v != NULL;
    fr.v      = v;
    fr.bBox   = TRUE;
    copy_mat(box, fr.box);

    return write_trxframe_indexed(status, &fr, nind, ind, gc);
}

void close_trx(t_trxstatus *status)
{
    gmx_tng_close(&status->tng);
    if (status->fio)
    {
        gmx_fio_close(status->fio);
    }
    sfree(status);
}

t_trxstatus *open_trx(const char *outfile, const char *filemode)
{
    t_trxstatus *stat;
    if (filemode[0] != 'w' && filemode[0] != 'a' && filemode[1] != '+')
    {
        gmx_fatal(FARGS, "Sorry, write_trx can only write");
    }

    snew(stat, 1);
    status_init(stat);

    stat->fio = gmx_fio_open(outfile, filemode);
    return stat;
}

static gmx_bool gmx_next_frame(t_trxstatus *status, t_trxframe *fr)
{
    t_trnheader sh;
    gmx_bool    bOK, bRet;

    bRet = FALSE;

    if (fread_trnheader(status->fio, &sh, &bOK))
    {
        fr->bDouble   = sh.bDouble;
        fr->natoms    = sh.natoms;
        fr->bStep     = TRUE;
        fr->step      = sh.step;
        fr->bTime     = TRUE;
        fr->time      = sh.t;
        fr->bLambda   = TRUE;
        fr->bFepState = TRUE;
        fr->lambda    = sh.lambda;
        fr->bBox      = sh.box_size > 0;
        if (fr->flags & (TRX_READ_X | TRX_NEED_X))
        {
            if (fr->x == NULL)
            {
                snew(fr->x, sh.natoms);
            }
            fr->bX = sh.x_size > 0;
        }
        if (fr->flags & (TRX_READ_V | TRX_NEED_V))
        {
            if (fr->v == NULL)
            {
                snew(fr->v, sh.natoms);
            }
            fr->bV = sh.v_size > 0;
        }
        if (fr->flags & (TRX_READ_F | TRX_NEED_F))
        {
            if (fr->f == NULL)
            {
                snew(fr->f, sh.natoms);
            }
            fr->bF = sh.f_size > 0;
        }
        if (fread_htrn(status->fio, &sh, fr->box, fr->x, fr->v, fr->f))
        {
            bRet = TRUE;
        }
        else
        {
            fr->not_ok = DATA_NOT_OK;
        }
    }
    else
    if (!bOK)
    {
        fr->not_ok = HEADER_NOT_OK;
    }

    return bRet;
}

static gmx_bool pdb_next_x(t_trxstatus *status, FILE *fp, t_trxframe *fr)
{
    t_atoms   atoms;
    matrix    boxpdb;
    int       ePBC, model_nr, na;
    char      title[STRLEN], *time;
    double    dbl;

    atoms.nr      = fr->natoms;
    atoms.atom    = NULL;
    atoms.pdbinfo = NULL;
    /* the other pointers in atoms should not be accessed if these are NULL */
    model_nr = NOTSET;
    na       = read_pdbfile(fp, title, &model_nr, &atoms, fr->x, &ePBC, boxpdb, TRUE, NULL);
    set_trxframe_ePBC(fr, ePBC);
    if (nframes_read(status) == 0)
    {
        fprintf(stderr, " '%s', %d atoms\n", title, fr->natoms);
    }
    fr->bPrec = TRUE;
    fr->prec  = 10000;
    fr->bX    = TRUE;
    fr->bBox  = (boxpdb[XX][XX] != 0.0);
    if (fr->bBox)
    {
        copy_mat(boxpdb, fr->box);
    }

    if (model_nr != NOTSET)
    {
        fr->bStep = TRUE;
        fr->step  = model_nr;
    }
    time = strstr(title, " t= ");
    if (time)
    {
        fr->bTime = TRUE;
        sscanf(time+4, "%lf", &dbl);
        fr->time = (real)dbl;
    }
    else
    {
        fr->bTime = FALSE;
        /* this is a bit dirty, but it will work: if no time is read from
           comment line in pdb file, set time to current frame number */
        if (fr->bStep)
        {
            fr->time = (real)fr->step;
        }
        else
        {
            fr->time = (real)nframes_read(status);
        }
    }
    if (na == 0)
    {
        return FALSE;
    }
    else
    {
        if (na != fr->natoms)
        {
            gmx_fatal(FARGS, "Number of atoms in pdb frame %d is %d instead of %d",
                      nframes_read(status), na, fr->natoms);
        }
        return TRUE;
    }
}

static int pdb_first_x(t_trxstatus *status, FILE *fp, t_trxframe *fr)
{
    initcount(status);

    fprintf(stderr, "Reading frames from pdb file");
    frewind(fp);
    get_pdb_coordnum(fp, &fr->natoms);
    if (fr->natoms == 0)
    {
        gmx_fatal(FARGS, "\nNo coordinates in pdb file\n");
    }
    frewind(fp);
    snew(fr->x, fr->natoms);
    pdb_next_x(status, fp, fr);

    return fr->natoms;
}

gmx_bool read_next_frame(const output_env_t oenv, t_trxstatus *status, t_trxframe *fr)
{
    real     pt;
    int      ct;
    gmx_bool bOK, bRet, bMissingData = FALSE, bSkip = FALSE;
    int      dummy = 0;
    int      ftp;

    bRet = FALSE;
    pt   = fr->tf;

    do
    {
        clear_trxframe(fr, FALSE);
        fr->tppf = fr->tpf;
        fr->tpf  = fr->tf;

        if (status->tng)
        {
            /* Special treatment for TNG files */
            ftp = efTNG;
        }
        else
        {
            ftp = gmx_fio_getftp(status->fio);
        }
        switch (ftp)
        {
            case efTRR:
                bRet = gmx_next_frame(status, fr);
                break;
            case efCPT:
                /* Checkpoint files can not contain mulitple frames */
                break;
            case efG96:
                read_g96_conf(gmx_fio_getfp(status->fio), NULL, fr,
                              status->persistent_line);
                bRet = (fr->natoms > 0);
                break;
            case efXTC:
                /* B. Hess 2005-4-20
                 * Sometimes is off by one frame
                 * and sometimes reports frame not present/file not seekable
                 */
                /* DvdS 2005-05-31: this has been fixed along with the increased
                 * accuracy of the control over -b and -e options.
                 */
                if (bTimeSet(TBEGIN) && (fr->tf < rTimeValue(TBEGIN)))
                {
                    if (xtc_seek_time(status->fio, rTimeValue(TBEGIN), fr->natoms, TRUE))
                    {
                        gmx_fatal(FARGS, "Specified frame (time %f) doesn't exist or file corrupt/inconsistent.",
                                  rTimeValue(TBEGIN));
                    }
                    initcount(status);
                }
                bRet = read_next_xtc(status->fio, fr->natoms, &fr->step, &fr->time, fr->box,
                                     fr->x, &fr->prec, &bOK);
                fr->bPrec = (bRet && fr->prec > 0);
                fr->bStep = bRet;
                fr->bTime = bRet;
                fr->bX    = bRet;
                fr->bBox  = bRet;
                if (!bOK)
                {
                    /* Actually the header could also be not ok,
                       but from bOK from read_next_xtc this can't be distinguished */
                    fr->not_ok = DATA_NOT_OK;
                }
                break;
            case efTNG:
                bRet = gmx_read_next_tng_frame(status->tng, fr, NULL, 0);
                break;
            case efPDB:
                bRet = pdb_next_x(status, gmx_fio_getfp(status->fio), fr);
                break;
            case efGRO:
                bRet = gro_next_x_or_v(gmx_fio_getfp(status->fio), fr);
                break;
            default:
#ifdef GMX_USE_PLUGINS
                bRet = read_next_vmd_frame(fr);
#else
                gmx_fatal(FARGS, "DEATH HORROR in read_next_frame ftp=%s,status=%s",
                          ftp2ext(gmx_fio_getftp(status->fio)),
                          gmx_fio_getname(status->fio));
#endif
        }
        fr->tf = fr->time;

        if (bRet)
        {
            bMissingData = (((fr->flags & TRX_NEED_X) && !fr->bX) ||
                            ((fr->flags & TRX_NEED_V) && !fr->bV) ||
                            ((fr->flags & TRX_NEED_F) && !fr->bF));
            bSkip = FALSE;
            if (!bMissingData)
            {
                ct = check_times2(fr->time, fr->t0, fr->bDouble);
                if (ct == 0 || ((fr->flags & TRX_DONT_SKIP) && ct < 0))
                {
                    printcount(status, oenv, fr->time, FALSE);
                }
                else if (ct > 0)
                {
                    bRet = FALSE;
                }
                else
                {
                    printcount(status, oenv, fr->time, TRUE);
                    bSkip = TRUE;
                }
            }
        }

    }
    while (bRet && (bMissingData || bSkip));

    if (!bRet)
    {
        printlast(status, oenv, pt);
        if (fr->not_ok)
        {
            printincomp(status, fr);
        }
    }

    return bRet;
}

int read_first_frame(const output_env_t oenv, t_trxstatus **status,
                     const char *fn, t_trxframe *fr, int flags)
{
    t_fileio      *fio;
    gmx_bool       bFirst, bOK;
    int            dummy = 0;
    int            ftp   = fn2ftp(fn);
    gmx_int64_t   *tng_ids;

    clear_trxframe(fr, TRUE);
    fr->flags = flags;

    bFirst = TRUE;

    snew((*status), 1);

    status_init( *status );
    (*status)->nxframe = 1;
    initcount(*status);

    if (efTNG == ftp)
    {
        /* Special treatment for TNG files */
        gmx_tng_open(fn, 'r', &(*status)->tng);
    }
    else
    {
        fio = (*status)->fio = gmx_fio_open(fn, "r");
    }
    switch (ftp)
    {
        case efTRR:
            break;
        case efCPT:
            read_checkpoint_trxframe(fio, fr);
            bFirst = FALSE;
            break;
        case efG96:
            /* Can not rewind a compressed file, so open it twice */
            if (!(*status)->persistent_line)
            {
                /* allocate the persistent line */
                snew((*status)->persistent_line, STRLEN+1);
            }
            read_g96_conf(gmx_fio_getfp(fio), fn, fr, (*status)->persistent_line);
            gmx_fio_close(fio);
            clear_trxframe(fr, FALSE);
            if (flags & (TRX_READ_X | TRX_NEED_X))
            {
                snew(fr->x, fr->natoms);
            }
            if (flags & (TRX_READ_V | TRX_NEED_V))
            {
                snew(fr->v, fr->natoms);
            }
            fio = (*status)->fio = gmx_fio_open(fn, "r");
            break;
        case efXTC:
            if (read_first_xtc(fio, &fr->natoms, &fr->step, &fr->time, fr->box, &fr->x,
                               &fr->prec, &bOK) == 0)
            {
                assert(!bOK);
                fr->not_ok = DATA_NOT_OK;
            }
            if (fr->not_ok)
            {
                fr->natoms = 0;
                printincomp(*status, fr);
            }
            else
            {
                fr->bPrec = (fr->prec > 0);
                fr->bStep = TRUE;
                fr->bTime = TRUE;
                fr->bX    = TRUE;
                fr->bBox  = TRUE;
                printcount(*status, oenv, fr->time, FALSE);
            }
            bFirst = FALSE;
            break;
        case efTNG:
            fr->step = -1;
            if (!gmx_read_next_tng_frame((*status)->tng, fr, NULL, 0))
            {
                fr->not_ok = DATA_NOT_OK;
                fr->natoms = 0;
                printincomp(*status, fr);
            }
            else
            {
                printcount(*status, oenv, fr->time, FALSE);
            }
            bFirst = FALSE;
            break;
        case efPDB:
            pdb_first_x(*status, gmx_fio_getfp(fio), fr);
            if (fr->natoms)
            {
                printcount(*status, oenv, fr->time, FALSE);
            }
            bFirst = FALSE;
            break;
        case efGRO:
            if (gro_first_x_or_v(gmx_fio_getfp(fio), fr))
            {
                printcount(*status, oenv, fr->time, FALSE);
            }
            bFirst = FALSE;
            break;
        default:
#ifdef GMX_USE_PLUGINS
            fprintf(stderr, "The file format of %s is not a known trajectory format to GROMACS.\n"
                    "Please make sure that the file is a trajectory!\n"
                    "GROMACS will now assume it to be a trajectory and will try to open it using the VMD plug-ins.\n"
                    "This will only work in case the VMD plugins are found and it is a trajectory format supported by VMD.\n", fn);
            gmx_fio_fp_close(fio); /*only close the file without removing FIO entry*/
            if (!read_first_vmd_frame(fn, fr))
            {
                gmx_fatal(FARGS, "Not supported in read_first_frame: %s", fn);
            }
#else
            gmx_fatal(FARGS, "Not supported in read_first_frame: %s. Please make sure that the file is a trajectory.\n"
                      "GROMACS is not compiled with plug-in support. Thus it cannot read non-GROMACS trajectory formats using the VMD plug-ins.\n"
                      "Please compile with plug-in support if you want to read non-GROMACS trajectory formats.\n", fn);
#endif
            break;
    }
    fr->tf = fr->time;

    /* Return FALSE if we read a frame that's past the set ending time. */
    if (!bFirst && (!(fr->flags & TRX_DONT_SKIP) && check_times(fr->time) > 0))
    {
        fr->t0 = fr->time;
        return FALSE;
    }

    if (bFirst ||
        (!(fr->flags & TRX_DONT_SKIP) && check_times(fr->time) < 0))
    {
        /* Read a frame when no frame was read or the first was skipped */
        if (!read_next_frame(oenv, *status, fr))
        {
            return FALSE;
        }
    }
    fr->t0 = fr->time;

    /* We need the number of atoms for random-access XTC searching, even when
     * we don't have access to the actual frame data.
     */
    (*status)->natoms = fr->natoms;

    return (fr->natoms > 0);
}

/***** C O O R D I N A T E   S T U F F *****/

int read_first_x(const output_env_t oenv, t_trxstatus **status, const char *fn,
                 real *t, rvec **x, matrix box)
{
    t_trxframe fr;

    read_first_frame(oenv, status, fn, &fr, TRX_NEED_X);

    snew((*status)->xframe, 1);
    (*status)->nxframe   = 1;
    (*(*status)->xframe) = fr;
    *t                   = (*status)->xframe->time;
    *x                   = (*status)->xframe->x;
    copy_mat((*status)->xframe->box, box);

    return (*status)->xframe->natoms;
}

gmx_bool read_next_x(const output_env_t oenv, t_trxstatus *status, real *t,
                     rvec x[], matrix box)
{
    gmx_bool bRet;

    status->xframe->x = x;
    /*xframe[status].x = x;*/
    bRet = read_next_frame(oenv, status, status->xframe);
    *t   = status->xframe->time;
    copy_mat(status->xframe->box, box);

    return bRet;
}

void close_trj(t_trxstatus *status)
{
    gmx_tng_close(&status->tng);
    if (status->fio)
    {
        gmx_fio_close(status->fio);
    }

    /* The memory in status->xframe is lost here,
     * but the read_first_x/read_next_x functions are deprecated anyhow.
     * read_first_frame/read_next_frame and close_trx should be used.
     */
    sfree(status);
}

void rewind_trj(t_trxstatus *status)
{
    initcount(status);

    gmx_fio_rewind(status->fio);
}

/***** T O P O L O G Y   S T U F F ******/

t_topology *read_top(const char *fn, int *ePBC)
{
    int         epbc, natoms;
    t_topology *top;

    snew(top, 1);
    epbc = read_tpx_top(fn, NULL, NULL, &natoms, NULL, NULL, NULL, top);
    if (ePBC)
    {
        *ePBC = epbc;
    }

    return top;
}
