/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "do_fit.h"

#include <stdio.h>

#include <cmath>

#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{
namespace
{
/*!\brief Mass weighted root mean square deviation kernel.
 * This struct aids templating the loop for over the atoms for structure simlarity measurements.
 * The total rmsd = sqrt(sum(m*(x-y)*(x-y)/sum(m)).
 */
struct RMSDKernel
{
    /*!\brief The nominator for the rmsd calculation.
     * \param[in] m mass of the atom
     * \param[in] x coordinate
     * \param[in] y coordinate
     * \returns contribution of this atom pair to nominator of rmsd calculation
     */
    static inline real nominator(real m, const gmx::RVec &x, const gmx::RVec &y)
    {
        return m * distance2(x, y);
    }
    /*!\brief The denominator for the rmsd calculation.
     * \param[in] m mass of the atom
     * \returns contribution of this atom pair to denominator of rmsd calculation
     */
    static inline real denominator(real m, const gmx::RVec & /*x*/, const gmx::RVec & /*y*/)
    {
        return m;
    }
};
/*!\brief Mass weighted rho-measure kernel.
 * This struct aids templating the loop for over the atoms for structure simlarity measurements.
 * The total rho = sqrt(sum(m*(x-y)*(x-y)/sum(m*(x+y)*(x+y)/4)).
 */
struct RhoMeasureKernel
{
    /*!\brief The nominator for the rho-measure calculation.
     * \param[in] m mass of the atom
     * \param[in] x coordinate
     * \param[in] y coordinate
     * \returns contribution of this atom pair to nominator of rho-measure calculation
     */
    static inline real nominator(real m, const gmx::RVec &x, const gmx::RVec &y)
    {
        return RMSDKernel::nominator(m, x, y);
    }
    /*!\brief The denominator for the rho-measure calculation.
     * \param[in] m mass of the atom
     * \param[in] x coordinate
     * \param[in] y coordinate
     * \returns contribution of this atom pair to denominator of rho-measure calculation
     */
    static inline real denominator(real m, const gmx::RVec &x, const gmx::RVec &y)
    {
        rvec sum;
        rvec_add(x, y, sum);
        return m * norm2(sum) / 4.;
    }
};

/*!\brief Evaluate structural simlarity in the form sqrt(sum(nominator(mass,x,y))/sum(denominator(mass,x,y))).
   *\tparam F provides the nominator and denominator function for the similarity measure
   *\param[in] nAtoms number of atoms of structures to compare
   *\param[in] index Index for selecting a sub-set of atoms, that is applied to mass, x, and xp
   *\param[in] mass Masses of the atoms for similarity comparison
   *\param[in] x the coordinates of the reference structure
   *\param[in] xp the coordinates of the structure to compare
   *\returns Mass-weighted measure of similarity between two structures
 */
template <typename F>
real structureSimilarity(int nAtoms, const int *index, const real *mass, const rvec *x, const rvec *xp)
{
    real denominator = 0;
    real nominator   = 0;
    if (index != nullptr)
    {
        gmx::ArrayRef<const int> indexRef = {index, index+nAtoms};
        for (int i : indexRef)
        {
            nominator   += F::nominator(mass[i], x[i], xp[i]);
            denominator += F::denominator(mass[i], x[i], xp[i]);
        }
    }
    else
    {
        for (int i = 0; i < nAtoms; i++)
        {
            nominator   += F::nominator(mass[i], x[i], xp[i]);
            denominator += F::denominator(mass[i], x[i], xp[i]);
        }
    }
    return std::sqrt(nominator/denominator);
}
}   // namespace

real StructureSimilarityMeasure::rmsd(int nAtoms, const real *mass, const rvec *x, const rvec *xp, const int * ind )
{
    return structureSimilarity<RMSDKernel>(nAtoms, ind, mass, x, xp);
}

real StructureSimilarityMeasure::sizeIndependentRho(int nAtoms, const real *mass, const rvec *x, const rvec *xp, const int *index)
{
    return structureSimilarity<RhoMeasureKernel>(nAtoms, index, mass, x, xp);
}
} // namespace gmx

void calc_fit_R(int ndim, int natoms, const real *w_rls, const rvec *xp, rvec *x, matrix R)
{
    int      c, r, n, j, i, irot, s;
    double **omega, **om;
    double   d[2*DIM], xnr, xpc;
    matrix   vh, vk, u;
    real     mn;
    int      index;
    real     max_d;

    if (ndim != 3 && ndim != 2)
    {
        gmx_fatal(FARGS, "calc_fit_R called with ndim=%d instead of 3 or 2", ndim);
    }

    snew(omega, 2*ndim);
    snew(om, 2*ndim);
    for (i = 0; i < 2*ndim; i++)
    {
        snew(omega[i], 2*ndim);
        snew(om[i], 2*ndim);
    }

    for (i = 0; i < 2*ndim; i++)
    {
        d[i] = 0;
        for (j = 0; j < 2*ndim; j++)
        {
            omega[i][j] = 0;
            om[i][j]    = 0;
        }
    }

    /*calculate the matrix U*/
    clear_mat(u);
    for (n = 0; (n < natoms); n++)
    {
        if ((mn = w_rls[n]) != 0.0)
        {
            for (c = 0; (c < ndim); c++)
            {
                xpc = xp[n][c];
                for (r = 0; (r < ndim); r++)
                {
                    xnr      = x[n][r];
                    u[c][r] += mn*xnr*xpc;
                }
            }
        }
    }

    /*construct omega*/
    /*omega is symmetric -> omega==omega' */
    for (r = 0; r < 2*ndim; r++)
    {
        for (c = 0; c <= r; c++)
        {
            if (r >= ndim && c < ndim)
            {
                omega[r][c] = u[r-ndim][c];
                omega[c][r] = u[r-ndim][c];
            }
            else
            {
                omega[r][c] = 0;
                omega[c][r] = 0;
            }
        }
    }

    /*determine h and k*/
    jacobi(omega, 2*ndim, d, om, &irot);
    /*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
     * int     natoms = number of rows and columns
     * real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
     * real       **v = v[0..n-1][0..n-1] contains the vectors in columns
     * int      *irot = number of jacobi rotations
     */

    if (debug && irot == 0)
    {
        fprintf(debug, "IROT=0\n");
    }

    index = 0; /* For the compiler only */

    /* Copy only the first ndim-1 eigenvectors */
    for (j = 0; j < ndim-1; j++)
    {
        max_d = -1000;
        for (i = 0; i < 2*ndim; i++)
        {
            if (d[i] > max_d)
            {
                max_d = d[i];
                index = i;
            }
        }
        d[index] = -10000;
        for (i = 0; i < ndim; i++)
        {
            vh[j][i] = M_SQRT2*om[i][index];
            vk[j][i] = M_SQRT2*om[i+ndim][index];
        }
    }
    if (ndim == 3)
    {
        /* Calculate the last eigenvector as the outer-product of the first two.
         * This insures that the conformation is not mirrored and
         * prevents problems with completely flat reference structures.
         */
        cprod(vh[0], vh[1], vh[2]);
        cprod(vk[0], vk[1], vk[2]);
    }
    else if (ndim == 2)
    {
        /* Calculate the last eigenvector from the first one */
        vh[1][XX] = -vh[0][YY];
        vh[1][YY] =  vh[0][XX];
        vk[1][XX] = -vk[0][YY];
        vk[1][YY] =  vk[0][XX];
    }

    /* determine R */
    clear_mat(R);
    for (r = 0; r < ndim; r++)
    {
        for (c = 0; c < ndim; c++)
        {
            for (s = 0; s < ndim; s++)
            {
                R[r][c] += vk[s][r]*vh[s][c];
            }
        }
    }
    for (r = ndim; r < DIM; r++)
    {
        R[r][r] = 1;
    }

    for (i = 0; i < 2*ndim; i++)
    {
        sfree(omega[i]);
        sfree(om[i]);
    }
    sfree(omega);
    sfree(om);
}

void do_fit_ndim(int ndim, int natoms, real *w_rls, const rvec *xp, rvec *x)
{
    int    j, m, r, c;
    matrix R;
    rvec   x_old;

    /* Calculate the rotation matrix R */
    calc_fit_R(ndim, natoms, w_rls, xp, x, R);

    /*rotate X*/
    for (j = 0; j < natoms; j++)
    {
        for (m = 0; m < DIM; m++)
        {
            x_old[m] = x[j][m];
        }
        for (r = 0; r < DIM; r++)
        {
            x[j][r] = 0;
            for (c = 0; c < DIM; c++)
            {
                x[j][r] += R[r][c]*x_old[c];
            }
        }
    }
}

void do_fit(int natoms, real *w_rls, const rvec *xp, rvec *x)
{
    do_fit_ndim(3, natoms, w_rls, xp, x);
}

void reset_x_ndim(int ndim, int ncm, const int *ind_cm,
                  int nreset, const int *ind_reset,
                  rvec x[], const real mass[])
{
    int  i, m, ai;
    rvec xcm;
    real tm, mm;

    if (ndim > DIM)
    {
        gmx_incons("More than 3 dimensions not supported.");
    }
    tm = 0.0;
    clear_rvec(xcm);
    if (ind_cm != nullptr)
    {
        for (i = 0; i < ncm; i++)
        {
            ai = ind_cm[i];
            mm = mass[ai];
            for (m = 0; m < ndim; m++)
            {
                xcm[m] += mm*x[ai][m];
            }
            tm += mm;
        }
    }
    else
    {
        for (i = 0; i < ncm; i++)
        {
            mm = mass[i];
            for (m = 0; m < ndim; m++)
            {
                xcm[m] += mm*x[i][m];
            }
            tm += mm;
        }
    }
    for (m = 0; m < ndim; m++)
    {
        xcm[m] /= tm;
    }

    if (ind_reset != nullptr)
    {
        for (i = 0; i < nreset; i++)
        {
            rvec_dec(x[ind_reset[i]], xcm);
        }
    }
    else
    {
        for (i = 0; i < nreset; i++)
        {
            rvec_dec(x[i], xcm);
        }
    }
}

void reset_x(int ncm, const int *ind_cm,
             int nreset, const int *ind_reset,
             rvec x[], const real mass[])
{
    reset_x_ndim(3, ncm, ind_cm, nreset, ind_reset, x, mass);
}
