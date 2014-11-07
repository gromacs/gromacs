/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GENTOP_QGEN_H
#define GENTOP_QGEN_H

#include <stdio.h>
#include "poldata.h"
#include "gmx_resp.h"

enum {
    eQGEN_OK, eQGEN_NOTCONVERGED, eQGEN_NOSUPPORT, eQGEN_ERROR, eQGEN_NR
};

typedef struct gentop_qgen *gentop_qgen_t;

extern gentop_qgen_t
gentop_qgen_init(gmx_poldata_t pd, t_atoms *atoms,
                 gmx_atomprop_t aps,
                 rvec *x, ChargeGenerationModel eqg_model,
                 real hfac, int qtotal,
                 real epsr);

extern void
gentop_qgen_done(gentop_qgen_t qgen);

extern int
generate_charges_sm(FILE *fp, gentop_qgen_t qgen,
                    gmx_poldata_t pd, t_atoms *atoms,
                    real tol, int maxiter, gmx_atomprop_t aps,
                    real *chieq);

extern int
generate_charges(FILE *fp,
                 gentop_qgen_t qgen,
                 gmx_resp_t gr, const char *molname,
                 gmx_poldata_t pd,
                 t_atoms *atoms,
                 real tol, int maxiter, int maxcycle,
                 gmx_atomprop_t aps);

extern void
qgen_message(gentop_qgen_t qgen, int len, char buf[], gmx_resp_t gr);

extern gmx_bool
bSplitQ(ChargeGenerationModel iModel);

/* The routines below return NOTSET if something is out of the ordinary */
extern int gentop_qgen_get_nzeta(gentop_qgen_t qgen, int atom);

extern int gentop_qgen_get_row(gentop_qgen_t qgen, int atom, int z);

extern double gentop_qgen_get_q(gentop_qgen_t qgen, int atom, int z);

extern double gentop_qgen_get_zeta(gentop_qgen_t qgen, int atom, int z);

#endif
