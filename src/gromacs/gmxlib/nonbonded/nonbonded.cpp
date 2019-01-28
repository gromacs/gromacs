/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "nonbonded.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "thread_mpi/threads.h"

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_free_energy.h"
#include "gromacs/gmxlib/nonbonded/nb_generic.h"
#include "gromacs/gmxlib/nonbonded/nb_generic_cg.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/listed_forces/bonded.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void
gmx_nonbonded_set_kernel_pointers(t_nblist *nl)
{
    const char *     elec;
    const char *     elec_mod;
    const char *     vdw;
    const char *     vdw_mod;
    const char *     geom;

    nl->kernelptr_vf = nullptr;
    nl->kernelptr_v  = nullptr;
    nl->kernelptr_f  = nullptr;

    elec     = gmx_nbkernel_elec_names[nl->ielec];
    elec_mod = eintmod_names[nl->ielecmod];
    vdw      = gmx_nbkernel_vdw_names[nl->ivdw];
    vdw_mod  = eintmod_names[nl->ivdwmod];
    geom     = gmx_nblist_geometry_names[nl->igeometry];

    if (nl->type == GMX_NBLIST_INTERACTION_FREE_ENERGY)
    {
        nl->kernelptr_vf       = reinterpret_cast<void *>(gmx_nb_free_energy_kernel);
        nl->kernelptr_f        = reinterpret_cast<void *>(gmx_nb_free_energy_kernel);
        nl->simd_padding_width = 1;
    }
    else if (!gmx_strcasecmp_min(geom, "CG-CG"))
    {
        nl->kernelptr_vf       = reinterpret_cast<void *>(gmx_nb_generic_cg_kernel);
        nl->kernelptr_f        = reinterpret_cast<void *>(gmx_nb_generic_cg_kernel);
        nl->simd_padding_width = 1;
    }
    else
    {
        /* "Pick" the only remaining kernel, the generic one.
         * We only do this for particle-particle kernels; by leaving the water-optimized kernel
         * pointers to NULL, the water optimization will automatically be disabled for this interaction.
         */
        if (!gmx_strcasecmp_min(geom, "Particle-Particle"))
        {
            nl->kernelptr_vf       = reinterpret_cast<void *>(gmx_nb_generic_kernel);
            nl->kernelptr_f        = reinterpret_cast<void *>(gmx_nb_generic_kernel);
            nl->simd_padding_width = 1;
            if (debug)
            {
                fprintf(debug,
                        "WARNING - Slow generic NB kernel used for neighborlist with\n"
                        "    Elec: '%s', Modifier: '%s'\n"
                        "    Vdw:  '%s', Modifier: '%s'\n",
                        elec, elec_mod, vdw, vdw_mod);
            }
        }
    }
}

void do_nonbonded(const t_forcerec  *fr,
                  rvec               x[],
                  rvec               f_shortrange[],
                  const t_mdatoms   *mdatoms,
                  const t_blocka    *excl,
                  gmx_grppairener_t *grppener,
                  t_nrnb            *nrnb,
                  real              *lambda,
                  real              *dvdl,
                  int                nls,
                  int                eNL,
                  int                flags)
{
    t_nblist *        nlist;
    int               n, n0, n1, i, i0, i1;
    t_nblists *       nblists;
    nb_kernel_data_t  kernel_data;
    nb_kernel_t *     kernelptr = nullptr;
    rvec *            f;

    kernel_data.flags                   = flags;
    kernel_data.exclusions              = excl;
    kernel_data.lambda                  = lambda;
    kernel_data.dvdl                    = dvdl;

    if (eNL >= 0)
    {
        i0 = eNL;
        i1 = i0+1;
    }
    else
    {
        i0 = 0;
        i1 = eNL_NR;
    }

    if (nls >= 0)
    {
        n0 = nls;
        n1 = nls+1;
    }
    else
    {
        n0 = 0;
        n1 = fr->nnblists;
    }

    for (n = n0; (n < n1); n++)
    {
        nblists = &fr->nblists[n];

        /* Tabulated kernels hard-code a lot of assumptions about the
         * structure of these tables, but that's not worth fixing with
         * the group scheme due for removal soon. As a token
         * improvement, this assertion will stop code segfaulting if
         * someone assumes that extending the group-scheme table-type
         * enumeration is something that GROMACS supports. */
        static_assert(etiNR == 3, "");

        kernel_data.table_elec              = nblists->table_elec;
        kernel_data.table_vdw               = nblists->table_vdw;
        kernel_data.table_elec_vdw          = nblists->table_elec_vdw;

        {
            {
                /* Short-range */
                if (!(flags & GMX_NONBONDED_DO_SR))
                {
                    continue;
                }
                kernel_data.energygrp_elec          = grppener->ener[egCOULSR];
                kernel_data.energygrp_vdw           = grppener->ener[fr->bBHAM ? egBHAMSR : egLJSR];
                nlist = nblists->nlist_sr;
                f                                   = f_shortrange;
            }

            for (i = i0; (i < i1); i++)
            {
                if (nlist[i].nri > 0)
                {
                    if (flags & GMX_NONBONDED_DO_POTENTIAL)
                    {
                        /* Potential and force */
                        kernelptr = reinterpret_cast<nb_kernel_t *>(nlist[i].kernelptr_vf);
                    }
                    else
                    {
                        /* Force only, no potential */
                        kernelptr = reinterpret_cast<nb_kernel_t *>(nlist[i].kernelptr_f);
                    }

                    if (nlist[i].type != GMX_NBLIST_INTERACTION_FREE_ENERGY && (flags & GMX_NONBONDED_DO_FOREIGNLAMBDA))
                    {
                        /* We don't need the non-perturbed interactions */
                        continue;
                    }
                    /* Neighborlists whose kernelptr==NULL will always be empty */
                    if (kernelptr != nullptr)
                    {
                        (*kernelptr)(&(nlist[i]), x, f, const_cast<t_forcerec*>(fr),
                                     const_cast<t_mdatoms*>(mdatoms), &kernel_data, nrnb);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Non-empty neighborlist does not have any kernel pointer assigned.");
                    }
                }
            }
        }
    }
}
