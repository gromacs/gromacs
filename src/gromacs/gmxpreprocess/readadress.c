/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

#include <stdlib.h>
#include <string.h>

#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#define MAXPTR 254

static char adress_refs[STRLEN], adress_tf_grp_names[STRLEN], adress_cg_grp_names[STRLEN];

void read_adressparams(int *ninp_p, t_inpfile **inp_p, t_adress *adress, warninp_t wi)
{
    int         nadress_refs, i;
    const char *tmp;
    char       *ptr1[MAXPTR];


    int        ninp;
    t_inpfile *inp;

    ninp   = *ninp_p;
    inp    = *inp_p;

    EETYPE("adress_type",                adress->type,         eAdresstype_names);
    RTYPE ("adress_const_wf",            adress->const_wf,     1);
    RTYPE ("adress_ex_width",            adress->ex_width,     0);
    RTYPE ("adress_hy_width",            adress->hy_width,     0);
    RTYPE ("adress_ex_forcecap",         adress->ex_forcecap,     0);
    EETYPE("adress_interface_correction", adress->icor,         eAdressICtype_names);
    EETYPE("adress_site",                adress->site,         eAdressSITEtype_names);
    STYPE ("adress_reference_coords",    adress_refs,             NULL);
    STYPE ("adress_tf_grp_names",        adress_tf_grp_names,     NULL);
    STYPE ("adress_cg_grp_names",        adress_cg_grp_names,     NULL);
    EETYPE("adress_do_hybridpairs",      adress->do_hybridpairs, yesno_names);

    nadress_refs = str_nelem(adress_refs, MAXPTR, ptr1);

    for (i = 0; (i < nadress_refs); i++) /*read vector components*/
    {
        adress->refs[i] = strtod(ptr1[i], NULL);
    }
    for (; (i < DIM); i++) /*remaining undefined components of the vector set to zero*/
    {
        adress->refs[i] = 0;
    }

    *ninp_p   = ninp;
    *inp_p    = inp;
}

void do_adress_index(t_adress *adress, gmx_groups_t *groups, char **gnames, t_grpopts *opts, warninp_t wi)
{
    int      nr, i, j, k;
    char    *ptr1[MAXPTR];
    int      nadress_cg_grp_names, nadress_tf_grp_names;

    /* AdResS coarse grained groups input */

    nr = groups->grps[egcENER].nr;
    snew(adress->group_explicit, nr);
    for (i = 0; i < nr; i++)
    {
        adress->group_explicit[i] = TRUE;
    }
    adress->n_energy_grps = nr;

    nadress_cg_grp_names = str_nelem(adress_cg_grp_names, MAXPTR, ptr1);

    if (nadress_cg_grp_names > 0)
    {
        for (i = 0; i < nadress_cg_grp_names; i++)
        {
            /* search for the group name mathching the tf group name */
            k = 0;
            while ((k < nr) &&
                   gmx_strcasecmp(ptr1[i], (char*)(gnames[groups->grps[egcENER].nm_ind[k]])))
            {
                k++;
            }
            if (k == nr)
            {
                gmx_fatal(FARGS, "Adress cg energy group %s not found\n", ptr1[i]);
            }
            adress->group_explicit[k] = FALSE;
            printf ("AdResS: Energy group %s is treated as coarse-grained \n",
                    (char*)(gnames[groups->grps[egcENER].nm_ind[k]]));
        }
        /* set energy group exclusions between all coarse-grained and explicit groups */
        for (j = 0; j < nr; j++)
        {
            for (k = 0; k < nr; k++)
            {
                if (!(adress->group_explicit[k] == adress->group_explicit[j]))
                {
                    opts->egp_flags[nr * j + k] |= EGP_EXCL;
                    if (debug)
                    {
                        fprintf(debug, "AdResS excl %s %s \n",
                                (char*)(gnames[groups->grps[egcENER].nm_ind[j]]),
                                (char*)(gnames[groups->grps[egcENER].nm_ind[k]]));
                    }
                }
            }
        }
    }
    else
    {
        warning(wi, "For an AdResS simulation at least one coarse-grained energy group has to be specified in adress_cg_grp_names");
    }


    /* AdResS multiple tf tables input */
    nadress_tf_grp_names = str_nelem(adress_tf_grp_names, MAXPTR, ptr1);
    adress->n_tf_grps    = nadress_tf_grp_names;
    snew(adress->tf_table_index, nadress_tf_grp_names);

    nr = groups->grps[egcENER].nr;

    if (nadress_tf_grp_names > 0)
    {
        for (i = 0; i < nadress_tf_grp_names; i++)
        {
            /* search for the group name mathching the tf group name */
            k = 0;
            while ((k < nr) &&
                   gmx_strcasecmp(ptr1[i], (char*)(gnames[groups->grps[egcENER].nm_ind[k]])))
            {
                k++;
            }
            if (k == nr)
            {
                gmx_fatal(FARGS, "Adress tf energy group %s not found\n", ptr1[i]);
            }

            adress->tf_table_index[i] = k;
            if (debug)
            {
                fprintf(debug, "found tf group %s id %d \n", ptr1[i], k);
            }
            if (adress->group_explicit[k])
            {
                gmx_fatal(FARGS, "Thermodynamic force group %s is not a coarse-grained group in adress_cg_grp_names. The thermodynamic force has to act on the coarse-grained vsite of a molecule.\n", ptr1[i]);
            }

        }
    }
    /* end AdResS multiple tf tables input */


}
