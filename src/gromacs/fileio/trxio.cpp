/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "gromacs/fileio/trxio.h"

#include "config.h"

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <filesystem>
#include <optional>
#include <string>

#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/g96io.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/timecontrol.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;

#if GMX_USE_PLUGINS
#    include "gromacs/fileio/vmdio.h"
#endif

/* defines for frame counter output */
#define SKIP1 10
#define SKIP2 100
#define SKIP3 1000

struct t_trxstatus
{
    int  flags; /* flags for read_first/next_frame  */
    int  currentFrame;
    real t0;                 /* time of the first frame, needed  *
                              * for skipping frames with -dt     */
    real                 tf; /* internal frame time              */
    t_trxframe*          xframe;
    t_fileio*            fio;
    gmx_tng_trajectory_t tng;
    int                  natoms;
    char*                persistent_line; /* Persistent line for reading g96 trajectories */
#if GMX_USE_PLUGINS
    gmx_vmdplugin_t* vmdplugin;
#endif
};

/* utility functions */

gmx_bool bRmod_fd(double a, double b, double c, gmx_bool bDouble)
{
    int    iq;
    double tol;

    tol = 2 * (bDouble ? GMX_DOUBLE_EPS : GMX_FLOAT_EPS);

    iq = static_cast<int>((a - b + tol * a) / c);

    return std::fabs(a - b - c * iq) <= tol * std::fabs(a);
}


int check_times2(real t, real t0, gmx_bool bDouble)
{
    int r;

#if !GMX_DOUBLE
    /* since t is float, we can not use double precision for bRmod */
    bDouble = FALSE;
#endif

    r              = -1;
    auto startTime = timeValue(TimeControl::Begin);
    auto endTime   = timeValue(TimeControl::End);
    auto deltaTime = timeValue(TimeControl::Delta);
    if ((!startTime.has_value() || (t >= startTime.value()))
        && (!endTime.has_value() || (t <= endTime.value())))
    {
        if (deltaTime.has_value() && !bRmod_fd(t, t0, deltaTime.value(), bDouble))
        {
            r = -1;
        }
        else
        {
            r = 0;
        }
    }
    else if (endTime.has_value() && (t >= endTime.value()))
    {
        r = 1;
    }
    if (debug)
    {
        fprintf(debug,
                "t=%g, t0=%g, b=%g, e=%g, dt=%g: r=%d\n",
                t,
                t0,
                startTime.value_or(0),
                endTime.value_or(0),
                deltaTime.value_or(0),
                r);
    }
    return r;
}

int check_times(real t)
{
    return check_times2(t, t, FALSE);
}

static void initcount(t_trxstatus* status)
{
    status->currentFrame = -1;
}

static void status_init(t_trxstatus* status)
{
    status->flags           = 0;
    status->xframe          = nullptr;
    status->fio             = nullptr;
    status->currentFrame    = -1;
    status->t0              = 0;
    status->tf              = 0;
    status->persistent_line = nullptr;
    status->tng             = nullptr;
}


int nframes_read(t_trxstatus* status)
{
    return status->currentFrame;
}

bool trxio_should_print_count(const gmx_output_env_t* oenv, t_trxstatus* status)
{
    return ((status->currentFrame < 2 * SKIP1 || status->currentFrame % SKIP1 == 0)
            && (status->currentFrame < 2 * SKIP2 || status->currentFrame % SKIP2 == 0)
            && (status->currentFrame < 2 * SKIP3 || status->currentFrame % SKIP3 == 0)
            && output_env_get_trajectory_io_verbosity(oenv) != 0);
}

static void printcount_(t_trxstatus* status, const gmx_output_env_t* oenv, const char* l, real t)
{
    if (trxio_should_print_count(oenv, status))
    {
        fprintf(stderr, "\r%-14s %6d time %8.3f   ", l, status->currentFrame, output_env_conv_time(oenv, t));
        fflush(stderr);
    }
}

static void printcount(t_trxstatus* status, const gmx_output_env_t* oenv, real t, gmx_bool bSkip)
{
    status->currentFrame++;
    printcount_(status, oenv, bSkip ? "Skipping frame" : "Reading frame", t);
}

static void printlast(t_trxstatus* status, const gmx_output_env_t* oenv, real t)
{
    printcount_(status, oenv, "Last frame", t);
    fprintf(stderr, "\n");
    fflush(stderr);
}

static void printincomp(t_trxstatus* status, t_trxframe* fr)
{
    if (fr->not_ok & HEADER_NOT_OK)
    {
        fprintf(stderr, "WARNING: Incomplete header: nr %d time %g\n", status->currentFrame + 1, fr->time);
    }
    else if (fr->not_ok)
    {
        fprintf(stderr, "WARNING: Incomplete frame: nr %d time %g\n", status->currentFrame + 1, fr->time);
    }
    fflush(stderr);
}

int prec2ndec(real prec)
{
    if (prec <= 0)
    {
        gmx_fatal(FARGS, "DEATH HORROR prec (%g) <= 0 in prec2ndec", prec);
    }

    return gmx::roundToInt(std::log(prec) / std::log(10.0));
}

real ndec2prec(int ndec)
{
    return std::pow(10.0, ndec);
}

t_fileio* trx_get_fileio(t_trxstatus* status)
{
    return status->fio;
}

float trx_get_time_of_final_frame(t_trxstatus* status)
{
    t_fileio* stfio    = trx_get_fileio(status);
    int       filetype = gmx_fio_getftp(stfio);
    gmx_bool  bOK;
    float     lasttime = -1;

    if (filetype == efXTC)
    {
        lasttime = xdr_xtc_get_last_frame_time(
                gmx_fio_getfp(stfio), gmx_fio_getxdr(stfio), status->natoms, &bOK);
        if (!bOK)
        {
            gmx_fatal(FARGS, "Error reading last frame. Maybe seek not supported.");
        }
    }
    else if (filetype == efTNG)
    {
        gmx_tng_trajectory_t tng = status->tng;
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

void clear_trxframe(t_trxframe* fr, gmx_bool bFirst)
{
    fr->not_ok    = 0;
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
        fr->bDouble   = FALSE;
        fr->natoms    = -1;
        fr->step      = 0;
        fr->time      = 0;
        fr->lambda    = 0;
        fr->fep_state = 0;
        fr->atoms     = nullptr;
        fr->prec      = 0;
        fr->x         = nullptr;
        fr->v         = nullptr;
        fr->f         = nullptr;
        clear_mat(fr->box);
        fr->bPBC    = FALSE;
        fr->pbcType = PbcType::Unset;
        fr->bIndex  = false;
        fr->index   = nullptr;
    }
}

void setTrxFramePbcType(t_trxframe* fr, PbcType pbcType)
{
    fr->bPBC    = (pbcType == PbcType::Unset);
    fr->pbcType = pbcType;
}

int write_trxframe_indexed(t_trxstatus* status, const t_trxframe* fr, int nind, const int* ind, gmx_conect gc)
{
    char  title[STRLEN];
    rvec *xout = nullptr, *vout = nullptr, *fout = nullptr;
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
        case efTNG: break;
        default:
            if (!fr->bX)
            {
                gmx_fatal(FARGS, "Need coordinates to write a %s trajectory", ftp2ext(ftp));
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
            if (fr->bX)
            {
                snew(xout, nind);
                for (i = 0; i < nind; i++)
                {
                    copy_rvec(fr->x[ind[i]], xout[i]);
                }
            }
            break;
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
        default: break;
    }

    switch (ftp)
    {
        case efTNG: gmx_write_tng_from_trxframe(status->tng, fr, nind); break;
        case efXTC: write_xtc(status->fio, nind, fr->step, fr->time, fr->box, xout, prec); break;
        case efTRR:
            gmx_trr_write_frame(
                    status->fio, nframes_read(status), fr->time, fr->step, fr->box, nind, xout, vout, fout);
            break;
        case efGRO:
        case efPDB:
        case efBRK:
        case efENT:
            if (!fr->bAtoms)
            {
                gmx_fatal(FARGS, "Can not write a %s file without atom names", ftp2ext(ftp));
            }
            sprintf(title, "frame t= %.3f", fr->time);
            if (ftp == efGRO)
            {
                write_hconf_indexed_p(gmx_fio_getfp(status->fio),
                                      title,
                                      fr->atoms,
                                      nind,
                                      ind,
                                      fr->x,
                                      fr->bV ? fr->v : nullptr,
                                      fr->box);
            }
            else
            {
                write_pdbfile_indexed(gmx_fio_getfp(status->fio),
                                      title,
                                      fr->atoms,
                                      fr->x,
                                      PbcType::Unset,
                                      fr->box,
                                      ' ',
                                      fr->step,
                                      nind,
                                      ind,
                                      gc,
                                      FALSE);
            }
            break;
        case efG96:
            sprintf(title, "frame t= %.3f", fr->time);
            write_g96_conf(gmx_fio_getfp(status->fio), title, fr, nind, ind);
            break;
        default: gmx_fatal(FARGS, "Sorry, write_trxframe_indexed can not write %s", ftp2ext(ftp));
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
            sfree(xout);
            break;
        case efXTC: sfree(xout); break;
        default: break;
    }

    return 0;
}

t_trxstatus* trjtools_gmx_prepare_tng_writing(const std::filesystem::path& filename,
                                              char                         filemode,
                                              t_trxstatus*                 in,
                                              const std::filesystem::path& infile,
                                              const int                    natoms,
                                              const gmx_mtop_t*            mtop,
                                              gmx::ArrayRef<const int>     index,
                                              const char*                  index_group_name)
{
    if (filemode != 'w' && filemode != 'a')
    {
        gmx_incons("Sorry, can only prepare for TNG output.");
    }
    t_trxstatus* out;
    snew(out, 1);
    status_init(out);

    if (in != nullptr)
    {
        gmx_prepare_tng_writing(
                filename, filemode, &in->tng, &out->tng, natoms, mtop, index, index_group_name);
    }
    else if (fn2ftp(infile) == efTNG)
    {
        gmx_tng_trajectory_t tng_in;
        gmx_tng_open(infile, 'r', &tng_in);

        gmx_prepare_tng_writing(
                filename, filemode, &tng_in, &out->tng, natoms, mtop, index, index_group_name);
    }
    else
    {
        // we start from a file that is not a tng file or have been unable to load the
        // input file, so we need to populate the fields independently of it
        gmx_prepare_tng_writing(
                filename, filemode, nullptr, &out->tng, natoms, mtop, index, index_group_name);
    }
    return out;
}

void write_tng_frame(t_trxstatus* status, t_trxframe* frame)
{
    gmx_write_tng_from_trxframe(status->tng, frame, -1);
}

int write_trxframe(t_trxstatus* status, t_trxframe* fr, gmx_conect gc)
{
    char title[STRLEN];
    title[0] = '\0';
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
        case efTRR: break;
        default:
            if (!fr->bX)
            {
                gmx_fatal(FARGS,
                          "Need coordinates to write a %s trajectory",
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
            gmx_trr_write_frame(status->fio,
                                fr->step,
                                fr->time,
                                fr->lambda,
                                fr->box,
                                fr->natoms,
                                fr->bX ? fr->x : nullptr,
                                fr->bV ? fr->v : nullptr,
                                fr->bF ? fr->f : nullptr);
            break;
        case efGRO:
        case efPDB:
        case efBRK:
        case efENT:
            if (!fr->bAtoms)
            {
                gmx_fatal(FARGS,
                          "Can not write a %s file without atom names",
                          ftp2ext(gmx_fio_getftp(status->fio)));
            }
            sprintf(title, "frame t= %.3f", fr->time);
            if (gmx_fio_getftp(status->fio) == efGRO)
            {
                write_hconf_p(
                        gmx_fio_getfp(status->fio), title, fr->atoms, fr->x, fr->bV ? fr->v : nullptr, fr->box);
            }
            else
            {
                write_pdbfile(gmx_fio_getfp(status->fio),
                              title,
                              fr->atoms,
                              fr->x,
                              fr->bPBC ? fr->pbcType : PbcType::Unset,
                              fr->box,
                              ' ',
                              fr->step,
                              gc);
            }
            break;
        case efG96: write_g96_conf(gmx_fio_getfp(status->fio), title, fr, -1, nullptr); break;
        default:
            gmx_fatal(FARGS, "Sorry, write_trxframe can not write %s", ftp2ext(gmx_fio_getftp(status->fio)));
    }

    return 0;
}

int write_trx(t_trxstatus*   status,
              int            nind,
              const int*     ind,
              const t_atoms* atoms,
              int            step,
              real           time,
              matrix         box,
              rvec           x[],
              rvec*          v,
              gmx_conect     gc)
{
    t_trxframe fr;

    clear_trxframe(&fr, TRUE);
    fr.bStep  = TRUE;
    fr.step   = step;
    fr.bTime  = TRUE;
    fr.time   = time;
    fr.bAtoms = atoms != nullptr;
    fr.atoms  = const_cast<t_atoms*>(atoms);
    fr.bX     = TRUE;
    fr.x      = x;
    fr.bV     = v != nullptr;
    fr.v      = v;
    fr.bBox   = TRUE;
    copy_mat(box, fr.box);

    return write_trxframe_indexed(status, &fr, nind, ind, gc);
}

void close_trx(t_trxstatus* status)
{
    if (status == nullptr)
    {
        return;
    }
    gmx_tng_close(&status->tng);
    if (status->fio)
    {
        gmx_fio_close(status->fio);
    }
    sfree(status->persistent_line);
#if GMX_USE_PLUGINS
    delete status->vmdplugin;
#endif
    /* The memory in status->xframe is lost here,
     * but the read_first_x/read_next_x functions are deprecated anyhow.
     * read_first_frame/read_next_frame and close_trx should be used.
     */
    sfree(status);
}

void done_trx_xframe(t_trxstatus* status)
{
    done_frame(status->xframe);
    sfree(status->xframe);
}

t_trxstatus* open_trx(const std::filesystem::path& outfile, const char* filemode)
{
    t_trxstatus* stat;
    if (filemode[0] != 'w' && filemode[0] != 'a' && filemode[1] != '+')
    {
        gmx_fatal(FARGS, "Sorry, write_trx can only write");
    }

    snew(stat, 1);
    status_init(stat);

    stat->fio = gmx_fio_open(outfile, filemode);
    return stat;
}

static gmx_bool gmx_next_frame(t_trxstatus* status, t_trxframe* fr)
{
    gmx_trr_header_t sh;
    gmx_bool         bOK, bRet;

    bRet = FALSE;

    if (gmx_trr_read_frame_header(status->fio, &sh, &bOK))
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
        if (status->flags & (TRX_READ_X | TRX_NEED_X))
        {
            if (fr->x == nullptr)
            {
                snew(fr->x, sh.natoms);
            }
            fr->bX = sh.x_size > 0;
        }
        if (status->flags & (TRX_READ_V | TRX_NEED_V))
        {
            if (fr->v == nullptr)
            {
                snew(fr->v, sh.natoms);
            }
            fr->bV = sh.v_size > 0;
        }
        if (status->flags & (TRX_READ_F | TRX_NEED_F))
        {
            if (fr->f == nullptr)
            {
                snew(fr->f, sh.natoms);
            }
            fr->bF = sh.f_size > 0;
        }
        if (gmx_trr_read_frame_data(status->fio, &sh, fr->box, fr->x, fr->v, fr->f))
        {
            bRet = TRUE;
        }
        else
        {
            fr->not_ok = DATA_NOT_OK;
        }
    }
    else if (!bOK)
    {
        fr->not_ok = HEADER_NOT_OK;
    }

    return bRet;
}

static gmx_bool pdb_next_x(t_trxstatus* status, FILE* fp, t_trxframe* fr)
{
    t_atoms   atoms;
    t_symtab* symtab;
    matrix    boxpdb;
    // Initiate model_nr to -1 rather than NOTSET.
    // It is not worthwhile introducing extra variables in the
    // read_pdbfile call to verify that a model_nr was read.
    PbcType pbcType;
    int     model_nr = -1, na;
    char    title[STRLEN], *time, *step;
    double  dbl;

    atoms.nr      = fr->natoms;
    atoms.atom    = nullptr;
    atoms.pdbinfo = nullptr;
    /* the other pointers in atoms should not be accessed if these are NULL */
    snew(symtab, 1);
    open_symtab(symtab);
    na = read_pdbfile(fp, title, &model_nr, &atoms, symtab, fr->x, &pbcType, boxpdb, nullptr);
    free_symtab(symtab);
    sfree(symtab);
    setTrxFramePbcType(fr, pbcType);
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

    fr->step  = 0;
    step      = std::strstr(title, " step= ");
    fr->bStep = ((step != nullptr) && sscanf(step + 7, "%" SCNd64, &fr->step) == 1);

    dbl       = 0.0;
    time      = std::strstr(title, " t= ");
    fr->bTime = ((time != nullptr) && sscanf(time + 4, "%lf", &dbl) == 1);
    fr->time  = dbl;

    if (na == 0)
    {
        return FALSE;
    }
    else
    {
        if (na != fr->natoms)
        {
            gmx_fatal(FARGS,
                      "Number of atoms in pdb frame %d is %d instead of %d",
                      nframes_read(status),
                      na,
                      fr->natoms);
        }
        return TRUE;
    }
}

static int pdb_first_x(t_trxstatus* status, FILE* fp, t_trxframe* fr)
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

bool read_next_frame(const gmx_output_env_t* oenv, t_trxstatus* status, t_trxframe* fr)
{
    real     pt;
    int      ct;
    gmx_bool bOK, bMissingData = FALSE, bSkip = FALSE;
    bool     bRet = false;
    int      ftp;

    pt = status->tf;

    do
    {
        clear_trxframe(fr, FALSE);

        if (status->tng)
        {
            /* Special treatment for TNG files */
            ftp = efTNG;
        }
        else
        {
            ftp = gmx_fio_getftp(status->fio);
        }
        auto startTime = timeValue(TimeControl::Begin);
        switch (ftp)
        {
            case efTRR: bRet = gmx_next_frame(status, fr); break;
            case efCPT:
                /* Checkpoint files can not contain mulitple frames */
                break;
            case efG96:
            {
                t_symtab* symtab = nullptr;
                read_g96_conf(gmx_fio_getfp(status->fio), {}, nullptr, fr, symtab, status->persistent_line);
                bRet = (fr->natoms > 0);
                break;
            }
            case efXTC:
                if (startTime.has_value() && (status->tf < startTime.value()))
                {
                    if (xtc_seek_time(status->fio, startTime.value(), fr->natoms, TRUE))
                    {
                        gmx_fatal(FARGS,
                                  "Specified frame (time %f) doesn't exist or file "
                                  "corrupt/inconsistent.",
                                  startTime.value());
                    }
                    initcount(status);
                }
                bRet      = (read_next_xtc(
                                status->fio, fr->natoms, &fr->step, &fr->time, fr->box, fr->x, &fr->prec, &bOK)
                        != 0);
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
            case efTNG: bRet = gmx_read_next_tng_frame(status->tng, fr, nullptr, 0); break;
            case efPDB: bRet = pdb_next_x(status, gmx_fio_getfp(status->fio), fr); break;
            case efGRO: bRet = gro_next_x_or_v(gmx_fio_getfp(status->fio), fr); break;
            default:
#if GMX_USE_PLUGINS
                bRet = read_next_vmd_frame(status->vmdplugin, fr);
#else
                gmx_fatal(FARGS,
                          "DEATH HORROR in read_next_frame ftp=%s,status=%s",
                          ftp2ext(gmx_fio_getftp(status->fio)),
                          gmx_fio_getname(status->fio).string().c_str());
#endif
        }
        status->tf = fr->time;

        if (bRet)
        {
            bMissingData = ((((status->flags & TRX_NEED_X) != 0) && !fr->bX)
                            || (((status->flags & TRX_NEED_V) != 0) && !fr->bV)
                            || (((status->flags & TRX_NEED_F) != 0) && !fr->bF));
            bSkip        = FALSE;
            if (!bMissingData)
            {
                ct = check_times2(fr->time, status->t0, fr->bDouble);
                if (ct == 0 || ((status->flags & TRX_DONT_SKIP) && ct < 0))
                {
                    printcount(status, oenv, fr->time, FALSE);
                }
                else if (ct > 0)
                {
                    bRet = false;
                }
                else
                {
                    printcount(status, oenv, fr->time, TRUE);
                    bSkip = TRUE;
                }
            }
        }

    } while (bRet && (bMissingData || bSkip));

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

bool read_first_frame(const gmx_output_env_t*      oenv,
                      t_trxstatus**                status,
                      const std::filesystem::path& fn,
                      t_trxframe*                  fr,
                      int                          flags)
{
    t_fileio* fio = nullptr;
    gmx_bool  bFirst, bOK;
    int       ftp = fn2ftp(fn);

    clear_trxframe(fr, TRUE);

    bFirst = TRUE;

    snew((*status), 1);

    status_init(*status);
    initcount(*status);
    (*status)->flags = flags;

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
        case efTRR: break;
        case efCPT:
            read_checkpoint_trxframe(fio, fr);
            bFirst = FALSE;
            break;
        case efG96:
        {
            /* Can not rewind a compressed file, so open it twice */
            if (!(*status)->persistent_line)
            {
                /* allocate the persistent line */
                snew((*status)->persistent_line, STRLEN + 1);
            }
            t_symtab* symtab = nullptr;
            read_g96_conf(gmx_fio_getfp(fio), fn, nullptr, fr, symtab, (*status)->persistent_line);
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
            (*status)->fio = gmx_fio_open(fn, "r");
            break;
        }
        case efXTC:
            if (read_first_xtc(fio, &fr->natoms, &fr->step, &fr->time, fr->box, &fr->x, &fr->prec, &bOK) == 0)
            {
                GMX_RELEASE_ASSERT(!bOK,
                                   "Inconsistent results - OK status from read_first_xtc, but 0 "
                                   "atom coords read");
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
            if (!gmx_read_next_tng_frame((*status)->tng, fr, nullptr, 0))
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
#if GMX_USE_PLUGINS
            fprintf(stderr,
                    "The file format of %s is not a known trajectory format to GROMACS.\n"
                    "Please make sure that the file is a trajectory!\n"
                    "GROMACS will now assume it to be a trajectory and will try to open it using "
                    "the VMD plug-ins.\n"
                    "This will only work in case the VMD plugins are found and it is a trajectory "
                    "format supported by VMD.\n",
                    fn.string().c_str());
            gmx_fio_fp_close(fio); /*only close the file without removing FIO entry*/
            if (!read_first_vmd_frame(fn, &(*status)->vmdplugin, fr))
            {
                gmx_fatal(FARGS, "Not supported in read_first_frame: %s", fn.string().c_str());
            }
#else
            gmx_fatal(FARGS,
                      "Not supported in read_first_frame: %s. Please make sure that the file is a "
                      "trajectory.\n"
                      "GROMACS is not compiled with plug-in support. Thus it cannot read "
                      "non-GROMACS trajectory formats using the VMD plug-ins.\n"
                      "Please compile with plug-in support if you want to read non-GROMACS "
                      "trajectory formats.\n",
                      fn.string().c_str());
#endif
    }
    (*status)->tf = fr->time;

    /* Return FALSE if we read a frame that's past the set ending time. */
    if (!bFirst && (!(flags & TRX_DONT_SKIP) && check_times(fr->time) > 0))
    {
        (*status)->t0 = fr->time;
        return FALSE;
    }

    if (bFirst || (!(flags & TRX_DONT_SKIP) && check_times(fr->time) < 0))
    {
        /* Read a frame when no frame was read or the first was skipped */
        if (!read_next_frame(oenv, *status, fr))
        {
            return FALSE;
        }
    }
    (*status)->t0 = fr->time;

    /* We need the number of atoms for random-access XTC searching, even when
     * we don't have access to the actual frame data.
     */
    (*status)->natoms = fr->natoms;

    return (fr->natoms > 0);
}

/***** C O O R D I N A T E   S T U F F *****/

int read_first_x(const gmx_output_env_t*      oenv,
                 t_trxstatus**                status,
                 const std::filesystem::path& fn,
                 real*                        t,
                 rvec**                       x,
                 matrix                       box)
{
    t_trxframe fr;

    read_first_frame(oenv, status, fn, &fr, TRX_NEED_X);

    snew((*status)->xframe, 1);
    (*(*status)->xframe) = fr;
    *t                   = (*status)->xframe->time;
    *x                   = (*status)->xframe->x;
    copy_mat((*status)->xframe->box, box);

    return (*status)->xframe->natoms;
}

gmx_bool read_next_x(const gmx_output_env_t* oenv, t_trxstatus* status, real* t, rvec x[], matrix box)
{
    gmx_bool bRet;

    status->xframe->x = x;
    /*xframe[status].x = x;*/
    bRet = read_next_frame(oenv, status, status->xframe);
    *t   = status->xframe->time;
    copy_mat(status->xframe->box, box);

    return bRet;
}

void rewind_trj(t_trxstatus* status)
{
    initcount(status);

    gmx_fio_rewind(status->fio);
}

/***** T O P O L O G Y   S T U F F ******/

t_topology* read_top(const std::filesystem::path& fn, PbcType* pbcType)
{
    int         natoms;
    PbcType     pbcTypeFile;
    t_topology* top;

    snew(top, 1);
    pbcTypeFile = read_tpx_top(fn, nullptr, nullptr, &natoms, nullptr, nullptr, top);
    if (pbcType)
    {
        *pbcType = pbcTypeFile;
    }

    return top;
}
