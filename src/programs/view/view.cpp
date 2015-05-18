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

#include "config.h"

#include <stdio.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#ifdef GMX_X11

#include "gromacs/fileio/writeps.h"

#include "Xstuff.h"
#include "dialogs.h"
#include "gromacs.bm"
#include "molps.h"
#include "nmol.h"
#include "xutil.h"

static void dump_it(t_manager *man)
{
    t_psdata ps;

    ps = ps_open("view.ps", 0, 0, man->molw->wd.width, man->molw->wd.height);
    ps_draw_mol(ps, man);
    ps_close(ps);
}

static void done_gmx(t_x11 *x11, t_gmx *gmx)
{
    done_logo(x11, gmx->logo);
    done_pd(x11, gmx->pd);
    done_man(x11, gmx->man);
    done_dlgs(gmx);
    x11->UnRegisterCallback(x11, gmx->wd->self);
}

static void move_gmx(t_x11 *x11, t_gmx *gmx, int width, int height,
                     bool bSizePD)
{
    int y0, wl, hl;
#ifdef DEBUG
    fprintf(stderr, "Move gmx %dx%d\n", width, height);
#endif
    y0 = XTextHeight(x11->font);
    /* Resize PD-Menu */
    if (bSizePD)
    {
        XResizeWindow(x11->disp, gmx->pd->wd.self, width, y0);
    }

    XMoveWindow(x11->disp, gmx->man->wd.self, 0, y0+1);
    XResizeWindow(x11->disp, gmx->man->wd.self, width, height-y0-1);

    wl = gmx->logo->wd.width;
    hl = gmx->logo->wd.height;
    XMoveWindow(x11->disp, gmx->logo->wd.self, (width-wl)/2, (height-y0-hl)/2);
}

static bool HandleClient(t_x11 *x11, int ID, t_gmx *gmx)
{
    t_pulldown *pd;

    pd = gmx->pd;

    switch (ID)
    {
        /* File Menu */
        case IDDUMPWIN:
            write_gmx(x11, gmx, IDDODUMP);
            break;
        case IDDODUMP:
            if (gmx->man->bAnimate)
            {
                hide_but(x11, gmx->man->vbox);
            }
            dump_it(gmx->man);
            if (gmx->man->bAnimate)
            {
                show_but(x11, gmx->man->vbox);
            }
            break;
        case IDCLOSE:
        case IDIMPORT:
        case IDEXPORT:
            ShowDlg(gmx->dlgs[edExport]);
            break;
        case IDDOEXPORT:
            write_sto_conf(gmx->confout, *gmx->man->top.name,
                           &(gmx->man->top.atoms),
                           gmx->man->x, NULL, gmx->man->molw->ePBC, gmx->man->box);
            break;
        case IDQUIT:
            show_mb(gmx, emQuit);
            break;
        case IDTERM:
            done_gmx(x11, gmx);
            return true;

        /* Edit Menu */
        case IDEDITTOP:
            edit_file("topol.gmx");
            break;
        case IDEDITCOORDS:
            edit_file("confin.gmx");
            break;
        case IDEDITPARAMS:
            edit_file("mdparin.gmx");
            break;

        /* Display Menu */
        case IDFILTER:
            if (gmx->filter)
            {
                ShowDlg(gmx->dlgs[edFilter]);
            }
            break;
        case IDDOFILTER:
            do_filter(x11, gmx->man, gmx->filter);
            break;
        case IDANIMATE:
            check_pd_item(pd, IDANIMATE, toggle_animate(x11, gmx->man));
            break;
        case IDLABELSOFF:
            no_labels(x11, gmx->man);
            break;
        case IDRESETVIEW:
            reset_view(gmx->man->view);
            ExposeWin(x11->disp, gmx->man->molw->wd.self);
            break;
        case IDPHOTO:
            show_mb(gmx, emNotImplemented);
            break;
        case IDBONDOPTS:
            ShowDlg(gmx->dlgs[edBonds]);
            break;
        case IDTHIN:
            set_bond_type(x11, gmx->man->molw, eBThin);
            break;
        case IDFAT:
            set_bond_type(x11, gmx->man->molw, eBFat);
            break;
        case IDVERYFAT:
            set_bond_type(x11, gmx->man->molw, eBVeryFat);
            break;
        case IDBALLS:
            set_bond_type(x11, gmx->man->molw, eBSpheres);
            break;
        case IDNOBOX:
            set_box_type(x11, gmx->man->molw, esbNone);
            break;
        case IDRECTBOX:
            set_box_type(x11, gmx->man->molw, esbRect);
            break;
        case IDTRIBOX:
            set_box_type(x11, gmx->man->molw, esbTri);
            break;
        case IDTOBOX:
            set_box_type(x11, gmx->man->molw, esbTrunc);
            break;

        /* Analysis Menu */
        case IDBOND:
        case IDANGLE:
        case IDDIH:
        case IDRMS:
        case IDRDF:
        case IDENERGIES:
        case IDCORR:
            show_mb(gmx, emNotImplemented);
            break;

        /* Help Menu */
        case IDHELP:
            show_mb(gmx, emHelp);
            break;
        case IDABOUT:
            show_logo(x11, gmx->logo);
            break;

        default:
            break;
    }
    return false;
}

static bool MainCallBack(t_x11 *x11, XEvent *event, Window /*w*/, void *data)
{
    t_gmx    *gmx;
    int       nsel, width, height;
    bool      result;

    result = false;
    gmx    = (t_gmx *)data;
    switch (event->type)
    {
        case ButtonRelease:
            hide_pd(x11, gmx->pd);
            break;
        case ConfigureNotify:
            width  = event->xconfigure.width;
            height = event->xconfigure.height;
            if ((width != gmx->wd->width) || (height != gmx->wd->height))
            {
                move_gmx(x11, gmx, width, height, true);
            }
            break;
        case ClientMessage:
            hide_pd(x11, gmx->pd);
            nsel   = event->xclient.data.l[0];
            result = HandleClient(x11, nsel, gmx);
            break;
        default:
            break;
    }
    return result;
}

static t_mentry  FileMenu[] = {
    { 0,  IDEXPORT,   false,  "Export..." },
    { 0,  IDDUMPWIN,  false,  "Print"     },
    { 0,  IDQUIT,     false,  "Quit"      }
};

static t_mentry  DispMenu[] = {
    { 0,  IDFILTER,   false,  "Filter..." },
    { 0,  IDANIMATE,  false,  "Animate"   },
    { 0,  IDLABELSOFF,    false,  "Labels Off"},
    { 0,  IDRESETVIEW,    false,  "Reset View"},
    { 0,  IDBONDOPTS,     false,  "Options..."}
};

static t_mentry  HelpMenu[] = {
    { 0,  IDHELP,     false,  "Help"             },
    { 0,  IDABOUT,    false,  "About GROMACS..." }
};

static t_mentry *gmx_pd[] = { FileMenu, DispMenu, HelpMenu };

#define MSIZE asize(gmx_pd)

static int         gmx_pd_size[MSIZE] = {
    asize(FileMenu), asize(DispMenu), asize(HelpMenu)
};

static const char *MenuTitle[MSIZE] = {
    "File", "Display", "Help"
};

static void init_gmx(t_x11 *x11, char *program, int nfile, t_filenm fnm[],
                     const output_env_t oenv)
{
    Pixmap                pm;
    t_gmx                *gmx;
    XSizeHints            hints;
    int                   w0, h0;
    int                   natom, natom_trx;
    t_topology            top;
    int                   ePBC;
    matrix                box;
    t_trxframe            fr;
    t_trxstatus          *status;
    char                  quote[256];

    snew(gmx, 1);
    snew(gmx->wd, 1);

    ePBC = read_tpx_top(ftp2fn(efTPR, nfile, fnm),
                        NULL, box, &natom, NULL, NULL, NULL, &top);

    read_first_frame(oenv, &status, ftp2fn(efTRX, nfile, fnm), &fr, TRX_DONT_SKIP);
    close_trx(status);
    natom_trx = fr.natoms;

    /* Creates a simple window */
    w0 = DisplayWidth(x11->disp, x11->screen)-132;
    h0 = DisplayHeight(x11->disp, x11->screen)-140;
    bromacs(quote, 255);
    InitWin(gmx->wd, 0, 0, w0, h0, 3, quote);
    gmx->wd->self = XCreateSimpleWindow(x11->disp, x11->root,
                                        gmx->wd->x, gmx->wd->y,
                                        gmx->wd->width, gmx->wd->height,
                                        gmx->wd->bwidth, WHITE, BLACK);
    pm = XCreatePixmapFromBitmapData(x11->disp, x11->root,
                                     (char *)gromacs_bits, gromacs_width,
                                     gromacs_height,
                                     WHITE, BLACK, 1);
    hints.flags      = PMinSize;
    hints.min_width  = 2*EWIDTH+40;
    hints.min_height = EHEIGHT+LDHEIGHT+LEGHEIGHT+40;
    XSetStandardProperties(x11->disp, gmx->wd->self, gmx->wd->text, program,
                           pm, NULL, 0, &hints);

    x11->RegisterCallback(x11, gmx->wd->self, x11->root, MainCallBack, gmx);
    x11->SetInputMask(x11, gmx->wd->self,
                      ButtonPressMask     | ButtonReleaseMask |
                      OwnerGrabButtonMask | ExposureMask      |
                      StructureNotifyMask);

    /* The order of creating windows is important here! */
    /* Manager */
    gmx->man  = init_man(x11, gmx->wd->self, 0, 0, 1, 1, WHITE, BLACK, ePBC, box, oenv);
    gmx->logo = init_logo(x11, gmx->wd->self, false);

    /* Now put all windows in the proper place */
    move_gmx(x11, gmx, w0, h0, false);

    XMapWindow(x11->disp, gmx->wd->self);
    map_man(x11, gmx->man);

    /* Pull Down menu */
    gmx->pd = init_pd(x11, gmx->wd->self, gmx->wd->width,
                      x11->fg, x11->bg,
                      MSIZE, gmx_pd_size, gmx_pd, MenuTitle);

    /* Dialogs & Filters */

    gmx->filter = init_filter(&(top.atoms), ftp2fn_null(efNDX, nfile, fnm),
                              natom_trx);

    init_dlgs(x11, gmx);

    /* Now do file operations */
    set_file(x11, gmx->man, ftp2fn(efTRX, nfile, fnm), ftp2fn(efTPR, nfile, fnm));

    ShowDlg(gmx->dlgs[edFilter]);
}
#endif

int gmx_view(int argc, char *argv[])
{
    const char  *desc[] = {
        "[THISMODULE] is the GROMACS trajectory viewer. This program reads a",
        "trajectory file, a run input file and an index file and plots a",
        "3D structure of your molecule on your standard X Window",
        "screen. No need for a high end graphics workstation, it even",
        "works on Monochrome screens.[PAR]",
        "The following features have been implemented:",
        "3D view, rotation, translation and scaling of your molecule(s),",
        "labels on atoms, animation of trajectories,",
        "hardcopy in PostScript format, user defined atom-filters",
        "runs on MIT-X (real X), open windows and motif,",
        "user friendly menus, option to remove periodicity, option to",
        "show computational box.[PAR]",
        "Some of the more common X command line options can be used: ",
        "[TT]-bg[tt], [TT]-fg[tt] change colors, [TT]-font fontname[tt] changes the font."
    };
    const char  *bugs[] = {
        "Balls option does not work",
        "Some times dumps core without a good reason"
    };

    output_env_t oenv;
    t_filenm     fnm[] = {
        { efTRX, "-f", NULL, ffREAD },
        { efTPR, NULL, NULL, ffREAD },
        { efNDX, NULL, NULL, ffOPTRD }
    };
#define NFILE asize(fnm)

    if (parse_common_args(&argc, argv, PCA_CAN_TIME, NFILE, fnm,
                          0, NULL, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
#ifndef GMX_X11
        fprintf(stderr, "Compiled without X-Windows - can not run viewer.\n");
#else
        t_x11 *x11;

        if ((x11 = GetX11(&argc, argv)) == NULL)
        {
            fprintf(stderr, "Can't connect to X Server.\n"
                    "Check your DISPLAY environment variable\n");
        }
        else
        {
            init_gmx(x11, argv[0], NFILE, fnm, oenv);
            x11->MainLoop(x11);
            x11->CleanUp(x11);
        }
#endif
    }
    return 0;
}
