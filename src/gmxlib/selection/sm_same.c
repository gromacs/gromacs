/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief Implementation of the \p same selection method.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <macros.h>
#include <smalloc.h>

#include <indexutil.h>
#include <selmethod.h>

/*! \internal \brief
 * Data structure for the \p same selection method.
 *
 * All angle values are in the units of radians.
 */
typedef struct
{
    /*! \brief Input group. */
    gmx_ana_index_t       g;
} t_methoddata_same;

/*! \brief Allocates data for the \p same selection method. */
static void *
init_data_same(int npar, gmx_ana_selparam_t *param);
/*! \brief Initializes the \p same selection method. */
static int
init_same(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/*! \brief Frees the data allocated for the \p same selection method. */
static void
free_data_same(void *data);
/*! \brief Initializes the evaluation of the \p same selection method for a frame. */
static int
init_frame_same(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data);
/*! \brief Evaluates the \p same selection method. */
static int
evaluate_same(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);

/*! \brief Enum array for the \p same selection method. */
static char *same_enum[] = { NULL, "residue", NULL };

/*! \brief Parameters for the \p same selection method. */
static gmx_ana_selparam_t smparams_same[] = {
    {NULL, {STR_VALUE,   1, {same_enum}}, NULL, SPAR_ENUMVAL},
    {"as", {GROUP_VALUE, 1, {NULL}},      NULL, SPAR_DYNAMIC},
};

/*! \internal \brief Selection method data for the \p same method. */
gmx_ana_selmethod_t sm_same = {
    "same", GROUP_VALUE, SMETH_REQTOP,
    asize(smparams_same), smparams_same,
    &init_data_same,
    NULL,
    NULL,
    NULL,
    &free_data_same,
    NULL,
    &evaluate_same,
    NULL,
};

/*!
 * \param[in]     npar  Not used (should be 2).
 * \param[in,out] param Method parameters (should point to 
 *   \ref smparams_same).
 * \returns Pointer to the allocated data (\ref t_methoddata_same).
 */
static void *
init_data_same(int npar, gmx_ana_selparam_t *param)
{
    t_methoddata_same *data;

    snew(data, 1);
    gmx_ana_index_clear(&data->g);
    param[1].val.u.g = &data->g;
    return data;
}

/*!
 * \param data Data to free (should point to a \ref t_methoddata_same).
 */
static void
free_data_same(void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;

    gmx_ana_index_deinit(&d->g);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_same.
 *
 * Calculates which atoms in \p g are in the same residues as the atoms in
 * \c t_methoddata_same::g.
 */
static int
evaluate_same(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;
    int                    i, j, resind;

    out->u.g->isize = 0;
    i = j = 0;
    while (i < d->g.isize)
    {
        resind = top->atoms.atom[d->g.index[i]].resind;
        /* Find the first atom in the current residue */
        while (top->atoms.atom[g->index[j]].resind < resind) ++j;
        /* Copy all the atoms in the residue to the output */
        while (j < g->isize && top->atoms.atom[g->index[j]].resind == resind)
        {
            out->u.g->index[out->u.g->isize++] = g->index[j];
            ++j;
        }
        /* Skip the rest of the atoms in the residue */
        ++i;
        while (i < d->g.isize && top->atoms.atom[d->g.index[i]].resind == resind) ++i;
    }
    return 0;
}

