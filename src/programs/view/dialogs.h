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

#ifndef _dialogs_h
#define _dialogs_h

#include "logo.h"
#include "manager.h"
#include "pulldown.h"
#include "xdlg.h"

typedef struct {
    bool      bMapped;
    t_dlg    *dlg;
} t_dialogs;

enum eDialogs {
    edExport, edBonds, edFilter, edNR
};

enum eMBoxes {
    emQuit, emHelp, emAbout, emNotImplemented, emNR
};

enum eExport {
    eExpGromos, eExpPDB, eExpNR
};

typedef struct {
    char           confout[256]; /* Export file			*/
    int            ExpMode;      /* Export mode			*/
    t_dlg        **dlgs;         /* Temporary storage for dlgs	*/
    int            which_mb;     /* Which mb is visible          */
    t_dlg        **mboxes;       /* id for message boxes         */
    t_filter      *filter;       /* Filter for visibility etc.	*/
    t_windata     *wd;           /* The main window		*/
    t_pulldown    *pd;           /* The pull-down menu		*/
    t_manager     *man;          /* The manager			*/
    /*t_statrec    *sr;*/		/* The statistics dlg		*/
    t_logo        *logo;         /* The gromacs logo             */
} t_gmx;

enum {
    IDNEW, IDOPEN, IDOPENED, IDCLOSE, IDIMPORT, IDEXPORT, IDDOEXPORT, IDQUIT, IDTERM,
    IDEDITTOP, IDEDITCOORDS, IDEDITPARAMS,
    IDGROMPP, IDRUNMD, IDDOGROMPP, IDGSTAT, IDDOGSTAT, IDDORUNMD,
    IDFILTER, IDDOFILTER,
    IDANIMATE, IDSHOWBOX, IDRMPBC, IDHYDROGEN, IDLABELSOFF, IDRESETVIEW, IDPHOTO,
    IDDUMPWIN, IDDODUMP,
    IDBONDOPTS, IDTHIN, IDFAT, IDVERYFAT, IDBALLS,
    IDNOBOX, IDRECTBOX, IDTRIBOX, IDTOBOX,
    IDBOND, IDANGLE, IDDIH, IDRMS, IDRDF, IDENERGIES, IDCORR,
    IDHELP, IDABOUT,

    /* Last line specifies how many IDs there are */
    IDMENUNR
};

extern void run_grompp(t_gmx *gmx);

extern void run_mdrun(t_gmx *gmx);

extern void write_gmx(t_x11 *x11, t_gmx *gmx, int mess);

/*extern void run_sr(t_statrec *sr);

   extern t_statrec *init_sr();*/

extern void init_dlgs(t_x11 *x11, t_gmx *gmx);

extern void show_mb(t_gmx *gmx, int mb);

extern void done_dlgs(t_gmx *gmx);

extern void edit_file(const char *fn);

extern t_filter *init_filter(t_atoms *atoms, const char *fn, int natom_trx);

extern t_dlg *select_filter(t_x11 *x11, t_gmx *gmx);

#endif
