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

#include "dialogs.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef HAVE_UNISTD_H
#    include <unistd.h> // for fork()
#endif

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "manager.h"
#include "nmol.h"
#include "x11.h"
#include "xdlghi.h"
#include "xmb.h"

#define MBFLAGS /* MB_APPLMODAL | */ MB_DONTSHOW

void write_gmx(t_x11* x11, t_gmx* gmx, int mess)
{
    XEvent letter;

    letter.type                 = ClientMessage;
    letter.xclient.display      = x11->disp;
    letter.xclient.window       = gmx->wd->self;
    letter.xclient.message_type = 0;
    letter.xclient.format       = 32;
    letter.xclient.data.l[0]    = mess;
    letter.xclient.data.l[1]    = Button1;
    XSendEvent(x11->disp, letter.xclient.window, True, 0, &letter);
}

static void shell_comm(const char* title, const char* script, int nsleep)
{
    FILE* tfil;
    char  command[STRLEN];
    char  tmp[32];

    std::strcpy(tmp, "dialogXXXXXX");
    tfil = gmx_fopen_temporary(tmp);

    fprintf(tfil, "%s\n", script);
    fprintf(tfil, "sleep %d\n", nsleep);
    gmx_ffclose(tfil);

    std::sprintf(command, "xterm -title %s -e sh %s", title, tmp);
#ifdef DEBUG
    std::fprintf(stderr, "command: %s\n", command);
#endif

    if (0 != std::system(command))
    {
        gmx_fatal(FARGS, "Failed to execute command: %s", command);
    }

#ifdef DEBUG
    unlink(tmp)
#endif
}

void show_mb(t_gmx* gmx, int mb)
{
    if (mb >= 0 && mb < emNR)
    {
        gmx->which_mb = mb;
        ShowDlg(gmx->mboxes[mb]);
    }
}

static void hide_mb(t_gmx* gmx)
{
    if (gmx->which_mb >= 0 && gmx->which_mb < emNR)
    {
        HideDlg(gmx->mboxes[gmx->which_mb]);
        gmx->which_mb = -1;
    }
}

static void MBCallback(t_x11* /*x11*/, int dlg_mess, int /*item_id*/, char* /*set*/, void* data)
{
    t_gmx* gmx;

    gmx = static_cast<t_gmx*>(data);
    if (dlg_mess == DLG_EXIT)
    {
        hide_mb(gmx);
    }
}

static t_dlg* about_mb(t_x11* x11, t_gmx* gmx)
{
    const char* lines[] = { "         G R O M A C S",
                            " Machine for Simulating Chemistry",
                            "       Copyright (c) 1992-2013",
                            "  Berk Hess, David van der Spoel, Erik Lindahl",
                            "        and many collaborators!" };

    return MessageBox(
            x11, gmx->wd->self, gmx->wd->text, asize(lines), lines, MB_OK | MB_ICONGMX | MBFLAGS, MBCallback, gmx);
}

static void QuitCB(t_x11* x11, int dlg_mess, int /*item_id*/, char* set, void* data)
{
    t_gmx* gmx;
    gmx = static_cast<t_gmx*>(data);

    hide_mb(gmx);
    if (dlg_mess == DLG_EXIT)
    {
        if (gmx_strcasecmp("yes", set) == 0)
        {
            write_gmx(x11, gmx, IDTERM);
        }
    }
}

static t_dlg* quit_mb(t_x11* x11, t_gmx* gmx)
{
    const char* lines[] = { " Do you really want to Quit ?" };

    return MessageBox(
            x11, gmx->wd->self, gmx->wd->text, asize(lines), lines, MB_YESNO | MB_ICONSTOP | MBFLAGS, QuitCB, gmx);
}

static t_dlg* help_mb(t_x11* x11, t_gmx* gmx)
{
    const char* lines[] = { " Help will soon be added" };

    return MessageBox(
            x11, gmx->wd->self, gmx->wd->text, asize(lines), lines, MB_OK | MB_ICONINFORMATION | MBFLAGS, MBCallback, gmx);
}

static t_dlg* ni_mb(t_x11* x11, t_gmx* gmx)
{
    const char* lines[] = { " This feature has not been", " implemented yet." };

    return MessageBox(
            x11, gmx->wd->self, gmx->wd->text, asize(lines), lines, MB_OK | MB_ICONEXCLAMATION | MBFLAGS, MBCallback, gmx);
}

enum
{
    eExE,
    eExGrom,
    eExPdb,
    eExConf,
    eExNR
};

static void ExportCB(t_x11* x11, int dlg_mess, int item_id, char* set, void* data)
{
    bool   bOk;
    t_gmx* gmx;
    t_dlg* dlg;

    gmx = static_cast<t_gmx*>(data);
    dlg = gmx->dlgs[edExport];
    switch (dlg_mess)
    {
        case DLG_SET:
            switch (item_id)
            {
                case eExGrom: gmx->ExpMode = eExpGromos; break;
                case eExPdb: gmx->ExpMode = eExpPDB; break;
                default: break;
            }
#ifdef DEBUG
            std::fprintf(stderr, "exportcb: item_id=%d\n", item_id);
#endif
            break;
        case DLG_EXIT:
            if ((bOk = gmx_strcasecmp("ok", set)) == 0)
            {
                std::strcpy(gmx->confout, EditText(dlg, eExConf));
            }
            HideDlg(dlg);
            if (bOk)
            {
                write_gmx(x11, gmx, IDDOEXPORT);
            }
            break;
    }
}

enum
{
    eg0,
    egTOPOL,
    egCONFIN,
    egPARAM,
    eg1,
    eg1PROC,
    eg32PROC
};

enum bond_set
{
    ebShowH = 11,
    ebDPlus,
    ebRMPBC,
    ebCue,
    ebSkip,
    ebWait
};

static void BondsCB(t_x11* x11, int dlg_mess, int item_id, char* set, void* data)
{
    static int ebond = -1;
    static int ebox  = -1;
    bool       bOk, bBond = false;
    int        nskip, nwait;
    t_gmx*     gmx;
    char*      endptr;

    gmx = static_cast<t_gmx*>(data);
    if (ebond == -1)
    {
        ebond = gmx->man->molw->bond_type;
        ebox  = gmx->man->molw->boxtype;
    }
    switch (dlg_mess)
    {
        case DLG_SET:
            if (item_id <= eBNR)
            {
                ebond = item_id - 1;
                bBond = false;
            }
            else if (item_id <= eBNR + esbNR + 1)
            {
                ebox  = item_id - eBNR - 2;
                bBond = true;
            }
            else
            {

#define DO_NOT(b) (b) = (!(b))

                switch (item_id)
                {
                    case ebShowH: toggle_hydrogen(x11, gmx->man->molw); break;
                    case ebDPlus: DO_NOT(gmx->man->bPlus);
#ifdef DEBUG
                        std::fprintf(stderr, "gmx->man->bPlus=%s\n", gmx->man->bPlus ? "true" : "false");
#endif
                        break;
                    case ebRMPBC: toggle_pbc(gmx->man); break;
                    case ebCue: DO_NOT(gmx->man->bSort);
#ifdef DEBUG
                        std::fprintf(stderr, "gmx->man->bSort=%s\n", gmx->man->bSort ? "true" : "false");
#endif
                        break;
                    case ebSkip:
                        nskip = std::strtol(set, &endptr, 10);
                        if (endptr != set)
                        {
#ifdef DEBUG
                            std::fprintf(stderr, "nskip: %d frames\n", nskip);
#endif
                            if (nskip >= 0)
                            {
                                gmx->man->nSkip = nskip;
                            }
                        }
                        break;
                    case ebWait:
                        nwait = std::strtol(set, &endptr, 10);
                        if (endptr != set)
                        {
#ifdef DEBUG
                            std::fprintf(stderr, "wait: %d ms\n", nwait);
#endif
                            if (nwait >= 0)
                            {
                                gmx->man->nWait = nwait;
                            }
                        }
                    default:
#ifdef DEBUG
                        std::fprintf(stderr, "item_id: %d, set: %s\n", item_id, set);
#endif
                        break;
                }
            }
            break;
        case DLG_EXIT:
            bOk = (gmx_strcasecmp("ok", set) == 0);
            HideDlg(gmx->dlgs[edBonds]);
            if (bOk)
            {
                if (bBond)
                {
                    switch (ebond)
                    {
                        case eBThin: write_gmx(x11, gmx, IDTHIN); break;
                        case eBFat: write_gmx(x11, gmx, IDFAT); break;
                        case eBVeryFat: write_gmx(x11, gmx, IDVERYFAT); break;
                        case eBSpheres: write_gmx(x11, gmx, IDBALLS); break;
                        default:
                            gmx_fatal(FARGS, "Invalid bond type %d at %s, %d", ebond, __FILE__, __LINE__);
                    }
                }
                else
                {
                    switch (ebox)
                    {
                        case esbNone: write_gmx(x11, gmx, IDNOBOX); break;
                        case esbRect: write_gmx(x11, gmx, IDRECTBOX); break;
                        case esbTri: write_gmx(x11, gmx, IDTRIBOX); break;
                        case esbTrunc: write_gmx(x11, gmx, IDTOBOX); break;
                        default:
                            gmx_fatal(FARGS, "Invalid box type %d at %s, %d", ebox, __FILE__, __LINE__);
                    }
                }
            }
            break;
    }
}

enum
{
    esFUNCT = 1,
    esBSHOW,
    esINFIL,
    esINDEXFIL,
    esLSQ,
    esSHOW,
    esPLOTFIL
};

typedef t_dlg* t_mmb(t_x11* x11, t_gmx* gmx);

typedef struct
{
    const char*  dlgfile;
    DlgCallback* cb;
} t_dlginit;

void init_dlgs(t_x11* x11, t_gmx* gmx)
{
    static t_dlginit di[]     = { { "export.dlg", ExportCB }, { "bonds.dlg", BondsCB } };
    static t_mmb*    mi[emNR] = { quit_mb, help_mb, about_mb, ni_mb };

    snew(gmx->dlgs, edNR);
    for (int i = 0; (i < asize(di)); i++)
    {
        gmx->dlgs[i] = ReadDlg(
                x11, gmx->wd->self, di[i].dlgfile, di[i].dlgfile, 0, 0, true, false, di[i].cb, gmx);
    }

    gmx->dlgs[edFilter] = select_filter(x11, gmx);

    snew(gmx->mboxes, emNR);
    for (int i = 0; (i < emNR); i++)
    {
        gmx->mboxes[i] = mi[i](x11, gmx);
    }
    gmx->which_mb = -1;
}

void done_dlgs(t_gmx* gmx)
{
    int i;

    for (i = 0; (i < edNR); i++)
    {
        FreeDlg(gmx->dlgs[i]);
    }
    for (i = 0; (i < emNR); i++)
    {
        FreeDlg(gmx->mboxes[i]);
    }
}

void edit_file(const char* fn)
{
    if (fork() == 0)
    {
        char script[256];

        std::sprintf(script, "vi  %s", fn);
        shell_comm(fn, script, 0);
        std::exit(0);
    }
}
