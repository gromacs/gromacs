/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 *
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef VMDIO_H_
#define VMDIO_H_

#include "external/vmd_molfile/molfile_plugin.h"
#include "types/trx.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmxvmdplugin
{
    molfile_plugin_t *api;
    const char      * filetype;
    void            * handle;
    gmx_bool          bV;
};

int read_first_vmd_frame(int *status, const char *fn, struct trxframe *fr, int flags);
gmx_bool read_next_vmd_frame(int status, struct trxframe *fr);

#ifdef __cplusplus
}
#endif

#endif /* VMDIO_H_ */
