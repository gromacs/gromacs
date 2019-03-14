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

#include "manager.h"

#include "config.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#ifdef HAVE_UNISTD_H
#include <unistd.h> // for usleep()
#endif

#include "gromacs/fileio/tpxio.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "3dview.h"
#include "nmol.h"

static void add_object(t_manager *man, eObject eO, int ai, int aj)
{
    srenew(man->obj, ++man->nobj);
    man->obj[man->nobj-1].eO    = eO;
    man->obj[man->nobj-1].eV    = eVNormal;
    man->obj[man->nobj-1].color = WHITE;
    man->obj[man->nobj-1].ai    = ai;
    man->obj[man->nobj-1].aj    = aj;
    man->obj[man->nobj-1].z     = 0.0;
}

static void add_bonds(t_manager *man, t_functype func[],
                      t_ilist *b, bool bB[])
{
    bool        *bH = man->bHydro;
    t_iatom     *ia;
    t_iatom      type, ai, aj, ak;
    int          i, delta, ftype;

#ifdef DEBUG
    std::fprintf(stderr, "Going to make bonds from an ilist with %d entries\n", b->nr);
#endif
    ia = b->iatoms;
    for (i = 0; (i < b->nr); )
    {
        type  = ia[0];
        ai    = ia[1];
        ftype = func[type];
        delta = interaction_function[ftype].nratoms;

        if (ftype == F_SETTLE)
        {
            aj     = ia[2];
            ak     = ia[3];
            bB[ai] = bB[aj] = bB[ak] = true;
            add_object(man, eOHBond, ai, aj);
            add_object(man, eOHBond, ai, ak);
        }
        else if (IS_CHEMBOND(ftype))
        {
            aj = ia[2];
#ifdef DEBUG
            std::fprintf(stderr, "Adding bond from %d to %d\n", ai, aj);
#endif
            bB[ai] = bB[aj] = true;
            if (!(bH[ai] == bH[aj]))
            {
                add_object(man, eOHBond, ai, aj);
            }
            else if (!bH[ai] && !bH[aj])
            {
                add_object(man, eOBond, ai, aj);
            }
        }
#ifdef DEBUG
        std::fprintf(stderr, "Type: %5d, delta: %5d\n", type, delta);
#endif
        ia += delta+1;
        i  += delta+1;
    }
}

static void add_bpl(t_manager *man, t_idef *idef, bool bB[])
{
    int ftype;

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (IS_CHEMBOND(ftype) || ftype == F_SETTLE)
        {
            add_bonds(man, idef->functype, &idef->il[ftype], bB);
        }
    }
}

static int which_atom(t_manager *man, int x, int y)
{
#define DELTA 5
    int  i;
    iv2 *ix = man->ix;

    for (i = 0; (i < man->natom); i++)
    {
        if ((std::abs(ix[i][XX]-x) < DELTA) && (std::abs(ix[i][YY]-y) < DELTA))
        {
            if (man->bVis[i])
            {
                return (int) i;
            }
        }
    }
    return -1;
}

static void do_label(t_x11 *x11, t_manager *man, int x, int y, bool bSet)
{
    int             ai;
    unsigned long   col;

    if ((ai = which_atom(man, x, y)) != -1)
    {
        x = man->ix[ai][XX];
        y = man->ix[ai][YY];
        if (bSet && !man->bLabel[ai])
        {
            col             = WHITE;
            man->bLabel[ai] = true;
        }
        else if (!bSet && man->bLabel[ai])
        {
            col             = BLUE;
            man->bLabel[ai] = false;
        }
        else
        {
            return;
        }
        XSetForeground(x11->disp, x11->gc, col);
        XDrawString(x11->disp, man->molw->wd.self, x11->gc, x+2, y-2, man->szLab[ai],
                    std::strlen(man->szLab[ai]));
        XSetForeground(x11->disp, x11->gc, x11->fg);
    }
}

static void show_label(t_x11 *x11, t_manager *man, int x, int y)
{
    do_label(x11, man, x, y, true);
}

static void hide_label(t_x11 *x11, t_manager *man, int x, int y)
{
    do_label(x11, man, x, y, false);
}

void set_file(t_x11 *x11, t_manager *man, const char *trajectory,
              const char *status)
{
    gmx_atomprop_t    aps;
    t_tpxheader       sh;
    t_atoms          *at;
    bool             *bB;
    int               i;

    read_tpxheader(status, &sh, true);
    snew(man->ix, sh.natoms);
    snew(man->zz, sh.natoms);
    snew(man->col, sh.natoms);
    snew(man->size, sh.natoms);
    snew(man->vdw, sh.natoms);
    snew(man->bLabel, sh.natoms);
    snew(man->bVis, sh.natoms);
    for (i = 0; (i < sh.natoms); i++)
    {
        man->bVis[i] = false;
    }

    man->bPbc = false;

    snew(man->szLab, sh.natoms);
    snew(man->bHydro, sh.natoms);
    snew(bB, sh.natoms);
    read_tpx_top(status, nullptr, man->box, &man->natom, nullptr, nullptr, &man->top);
    man->gpbc = gmx_rmpbc_init(&man->top.idef, -1, man->natom);

    man->natom =
        read_first_x(man->oenv, &man->status, trajectory, &(man->time), &(man->x),
                     man->box);
    man->trajfile = gmx_strdup(trajectory);
    if (man->natom > man->top.atoms.nr)
    {
        gmx_fatal(FARGS, "Topology %s (%d atoms) and trajectory %s (%d atoms) "
                  "do not match", status, man->top.atoms.nr,
                  trajectory, man->natom);
    }

    man->title.text = gmx_strdup(gmx::formatString("%s: %s", *man->top.name, gmx::getCoolQuote().c_str()).c_str());
    man->view       = init_view(man->box);
    at              = &(man->top.atoms);
    aps             = gmx_atomprop_init();
    for (i = 0; (i < man->natom); i++)
    {
        char      *aname = *(at->atomname[i]);
        t_resinfo *ri    = &at->resinfo[at->atom[i].resind];

        man->col[i] = Type2Color(aname);
        snew(man->szLab[i], 20);
        if (ri->ic != ' ')
        {
            std::sprintf(man->szLab[i], "%s%d%c, %s", *ri->name, ri->nr, ri->ic, aname);
        }
        else
        {
            std::sprintf(man->szLab[i], "%s%d, %s", *ri->name, ri->nr, aname);
        }
        man->bHydro[i] = (toupper(aname[0]) == 'H');
        if (man->bHydro[i])
        {
            man->vdw[i] = 0;
        }
        else if (!gmx_atomprop_query(aps, epropVDW, *ri->name, aname, &(man->vdw[i])))
        {
            man->vdw[i] = 0;
        }
    }
    gmx_atomprop_destroy(aps);
    add_bpl(man, &(man->top.idef), bB);
    for (i = 0; (i < man->natom); i++)
    {
        if (!bB[i])
        {
            add_object(man, eOSingle, (int) i, 0);
        }
    }
    sfree(bB);

    ExposeWin(x11->disp, man->molw->wd.self);
}

void step_message(t_x11 *x11, t_manager *man)
{
    XEvent letter;

    letter.type                 = ClientMessage;
    letter.xclient.display      = x11->disp;
    letter.xclient.window       = man->wd.self;
    letter.xclient.message_type = 0;
    letter.xclient.format       = 32;
    letter.xclient.data.l[0]    = IDSTEP;
    letter.xclient.data.l[1]    = Button1;
    XSendEvent(x11->disp, letter.xclient.window, True, 0, &letter);
}

static void reset_mols(t_block *mols, matrix box, rvec x[])
{
    int  i, m0, m1, j, m;
    rvec xcm, icm;
    real ix, iy, iz;

    for (i = 0; (i < mols->nr); i++)
    {
        m0 = mols->index[i];
        m1 = mols->index[i+1];

        clear_rvec(xcm);
        clear_rvec(icm);

        for (j = m0; (j < m1); j++)
        {
            rvec_inc(xcm, x[j]);
        }
        for (m = 0; (m < DIM); m++)
        {
            xcm[m] /= (m1-m0);
        }
        for (m = 0; (m < DIM); m++)
        {
            if (xcm[m] < 0)
            {
                icm[m] = box[m][m];
            }
            else if (xcm[m] >= box[m][m])
            {
                icm[m] = -box[m][m];
            }
        }
        ix = icm[XX], iy = icm[YY], iz = icm[ZZ];

        if ((ix != 0) || (iy != 0) || (iz != 0))
        {
            for (j = m0; (j < m1); j++)
            {
                x[j][XX] += ix;
                x[j][YY] += iy;
                x[j][ZZ] += iz;
            }
        }
    }
}

static bool step_man(t_manager *man, int *nat)
{
    static int      ncount = 0;
    bool            bEof;

    if (!man->natom)
    {
        std::fprintf(stderr, "Not initiated yet!");
        std::exit(1);
    }
    bEof = read_next_x(man->oenv, man->status, &man->time, man->x, man->box);
    *nat = man->natom;
    if (ncount == man->nSkip)
    {
        auto atomsArrayRef = gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec *>(man->x), man->natom);
        switch (man->molw->boxtype)
        {
            case esbTri:
                put_atoms_in_triclinic_unitcell(ecenterDEF, man->box, atomsArrayRef);
                break;
            case esbTrunc:
                put_atoms_in_compact_unitcell(man->molw->ePBC, ecenterDEF, man->box,
                                              atomsArrayRef);
                break;
            case esbRect:
            case esbNone:
            default:
                break;
        }
        if (man->bPbc)
        {
            gmx_rmpbc(man->gpbc, man->natom, man->box, man->x);
            reset_mols(&(man->top.mols), man->box, man->x);
        }
        ncount = 0;
    }
    else
    {
        if (man->nSkip > 0)
        {
            ncount++;
            return step_man(man, nat);
        }
    }

    return bEof;
}

static void HandleClient(t_x11 *x11, t_manager *man, long data[])
{
    int  ID, button, x, y;
    bool bPos;
    real fac;

    ID     = data[0];
    button = data[1];
    x      = data[2];
    y      = data[3];
    bPos   = (button == Button1);
    switch (ID)
    {
        case IDROTX:
        case IDROTY:
        case IDROTZ:
            rotate_3d(man->view, ID-IDROTX, bPos);
            draw_mol(x11, man);
            break;
        case IDZOOM:
            if (bPos)
            {
                fac = 0.8; /* Reduce distance between eye and origin */
            }
            else
            {
                fac = 1.25;
            }

            /*  zoom changed to scale by Berk Hess 3-7-96
               if (zoom_3d(man->view,fac))
               draw_mol(x11,man); */
            man->view->sc_x /= fac;
            man->view->sc_y /= fac;
            draw_mol(x11, man);
            break;
        case IDTRANSX:
        case IDTRANSY:
        case IDTRANSZ:
            translate_view(man->view, ID-IDTRANSX, bPos);
            draw_mol(x11, man);
            break;
        case IDREWIND:
            if (man->status)
            {
                rewind_trj(man->status);
                read_next_x(man->oenv, man->status, &(man->time), man->x,
                            man->box);
                man->bEof = false;
                draw_mol(x11, man);
            }
            break;
        case IDSTEP:
        {
            int      nat;

            nat = 0;
            if (!step_man(man, &nat))
            {
                man->bEof  = true;
                man->bStop = true;
            }
            else
            {
                if (nat > 0)
                {
                    draw_mol(x11, man);
                    usleep(man->nWait*1000);
                }
            }
            break;
        }
        case IDFF:
            man->bStop = false;
            break;
        case IDSTOP_ANI:
            man->bStop = true;
            break;
        case IDDRAWMOL:
            draw_mol(x11, man);
            break;
        case IDLABEL:
            switch (button)
            {
                case Button1:
                case Button2:
                    show_label(x11, man, x, y);
                    break;
                case Button3:
                    hide_label(x11, man, x, y);
                    break;
            }
            break;
        default:
            break;
    }
    if (man->bAnimate && !man->bEof && !man->bStop)
    {
        step_message(x11, man);
    }
}

static bool TitleCallBack(t_x11 *x11, XEvent *event, Window /*w*/, void *data)
{
    t_windata *wd;

    wd = static_cast<t_windata *>(data);
    switch (event->type)
    {
        case Expose:
            if (wd->text && (wd->width > 10))
            {
                XSetForeground(x11->disp, x11->gc, WHITE);
                TextInWin(x11, wd, wd->text, eXCenter, eYCenter);
                XDrawLine(x11->disp, wd->self, x11->gc, 0, wd->height,
                          wd->width, wd->height);
            }
            break;
        case ConfigureNotify:
            wd->width  = event->xconfigure.width;
            wd->height = event->xconfigure.height;
            break;
    }
    return false;
}

static bool ManCallBack(t_x11 *x11, XEvent *event, Window /*w*/, void *data)
{
    t_manager *man;
    int        width, height;

    man = static_cast<t_manager *>(data);
    switch (event->type)
    {
        case ConfigureNotify:
            width  = event->xconfigure.width;
            height = event->xconfigure.height;
            if ((width != man->wd.width) || (height != man->wd.height))
            {
                move_man(x11, man, width, height);
            }
            break;
        case ClientMessage:
            HandleClient(x11, man, event->xclient.data.l);
            break;
        default:
            break;
    }
    return false;
}

void no_labels(t_x11 *x11, t_manager *man)
{
    int i;

    for (i = 0; (i < man->natom); i++)
    {
        man->bLabel[i] = false;
    }
    draw_mol(x11, man);
}

void move_man(t_x11 *x11, t_manager *man, int width, int height)
{
    int x0, y0, mw, mh, hb;
    int th;

#ifdef DEBUG
    std::fprintf(stderr, "Move manager %dx%d\n", width, height);
#endif
    man->wd.width  = width;
    man->wd.height = height;

    /* Move all subwindows, resize only Mol window */
    x0 = width-EWIDTH-AIR-4*BORDER;           /* Starting of ewin etc. */
    y0 = AIR;

    /* Mol Window */
    mw = x0-2*AIR-4*BORDER;
    mh = height-y0-AIR-2*BORDER;
    XMoveResizeWindow(x11->disp, man->molw->wd.self, AIR, y0, mw, mh);

    /* Title Window */
    th = XTextHeight(x11->font);
    XMoveResizeWindow(x11->disp, man->title.self, 0, 0, mw, th+AIR);

    /* Legend Window */
    XMoveResizeWindow(x11->disp, man->legw->wd.self, x0, y0, EWIDTH, LEGHEIGHT);
    y0 += LEGHEIGHT+AIR+2*BORDER;

    if (y0 > height)
    {
        std::printf("Error: Windows falling out of main window!\n");
    }

    /* Button Box */
    hb = height-y0-AIR-2*BORDER;
    XMoveResizeWindow(x11->disp, man->bbox->wd.self, x0, y0, EWIDTH, hb);

    /* Video Box */
    x0 = (mw-man->vbox->wd.width)/2;
    y0 = (mh-2-AIR-man->vbox->wd.height);
    XMoveWindow(x11->disp, man->vbox->wd.self, x0, y0);
}

void map_man(t_x11 *x11, t_manager *man)
{
    XMapWindow(x11->disp, man->wd.self);
    map_mw(x11, man->molw);
    XMapWindow(x11->disp, man->title.self);
    map_legw(x11, man->legw);
    show_but(x11, man->bbox);
}

bool toggle_animate (t_x11 *x11, t_manager *man)
{
    if (man->status)
    {
        man->bAnimate = !man->bAnimate;
        man->bStop    = true;
        man->bEof     = false;
        if (man->bAnimate)
        {
            show_but(x11, man->vbox);
        }
        else
        {
            hide_but(x11, man->vbox);
        }
    }
    return man->bAnimate;
}

bool toggle_pbc (t_manager *man)
{
    man->bPbc = !man->bPbc;

    return man->bPbc;
}


t_manager *init_man(t_x11 *x11, Window Parent,
                    int x, int y, int width, int height,
                    unsigned long fg, unsigned long bg,
                    int ePBC, matrix box,
                    gmx_output_env_t *oenv)
{
    t_manager *man;

    snew(man, 1);
    man->status = nullptr;
    man->bPlus  = true;
    man->bSort  = true;
    man->oenv   = oenv;
    InitWin(&(man->wd), x, y, width, height, 0, "Manager");
    man->wd.self = XCreateSimpleWindow(x11->disp, Parent, man->wd.x, man->wd.y,
                                       man->wd.width, man->wd.height,
                                       man->wd.bwidth, fg, bg);
    x11->RegisterCallback(x11, man->wd.self, Parent, ManCallBack, man);
    x11->SetInputMask(x11, man->wd.self, StructureNotifyMask |
                      ExposureMask | ButtonPressMask);

    /* The order of creating windows is important for the stacking order */
    /* Mol Window */
    man->molw = init_mw(x11, man->wd.self, 0, 0, 1, 1, WHITE, BLUE, ePBC, box);

    /* Title Window */
    InitWin(&(man->title), 0, 0, 1, 1, 0, nullptr);
    man->title.self = XCreateSimpleWindow(x11->disp, man->molw->wd.self,
                                          man->title.x, man->title.y,
                                          man->title.width, man->title.height,
                                          man->title.bwidth, WHITE, BLUE);
    x11->RegisterCallback(x11, man->title.self, man->molw->wd.self,
                          TitleCallBack, &(man->title));
    x11->SetInputMask(x11, man->title.self, ExposureMask | StructureNotifyMask);

    /* Button box */
    man->bbox = init_bbox(x11, man->wd.self, man->wd.self, 1, WHITE, BLUE);

    /* Legend Window */
    man->legw = init_legw(x11, man->wd.self, 0, 0, EWIDTH, LEGHEIGHT, WHITE, BLUE);

    /* Video Box */
    man->vbox = init_vbox(x11, man->molw->wd.self, man->wd.self, WHITE, BLUE);

    return man;
}

void done_man(t_x11 *x11, t_manager *man)
{
    done_bbox(x11, man->vbox);
    done_bbox(x11, man->bbox);
    done_mw(x11, man->molw);
    done_legw(x11, man->legw);
    x11->UnRegisterCallback(x11, man->title.self);
    x11->UnRegisterCallback(x11, man->wd.self);
    sfree(man->x);
    sfree(man->obj);
    sfree(man->bHydro);
    sfree(man->bLabel);
    sfree(man->szLab);
    sfree(man->col);
    sfree(man);
}

void do_filter(t_x11 *x11, t_manager *man, t_filter *filter)
{
    int      i;
    int      j;

    for (i = 0; (i < man->natom); i++)
    {
        man->bVis[i] = false;
    }
    for (i = 0; (i < filter->grps->nr); i++)
    {
        if (filter->bShow[i])
        {
            for (j = filter->grps->index[i]; (j < filter->grps->index[i+1]); j++)
            {
                man->bVis[filter->grps->a[j]] = true;
            }
        }
    }

    ExposeWin(x11->disp, man->wd.self);
}
