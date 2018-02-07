/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "eigenvec.h"

#include "config.h"

#include <memory>
#include <vector>

#include "gromacs/fileio/confio.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

void FcaEigenVec::internal_fitit(const int nr, rvec* x, const real* w_rls, rvec* transvec,
                                 matrix &rmat, const int nfit, const int* ifit, const rvec* xref)
{
    std::unique_ptr<rvec[]> xdum(new rvec[nr]());
    copy_rvecn(x, xdum.get(), 0, nr);
    if (nfit > 3)
    {
        reset_x(nfit, ifit, nr, nullptr, x, const_cast<real*>(w_rls)); //reset nr first atoms in x to mass-center of nfit/ifit atoms
    }
    for (int i = 0; (i < nr); i++) \
    {
        rvec_sub(x[i], xdum[i], transvec[i]); // zeigt von xfit nach xread
    }
    if (nfit > 3)
    {
        calc_fit_R(DIM, nr, const_cast<real*>(w_rls), const_cast<rvec*>(xref), x, rmat); //auswahl der atome in x und x_ref geschieht durch w_rls, 0 -> atom is not choosen
        fca_utils::rotate_vec(nr, x, rmat);
    }
}

FcaEigenVec::FcaEigenVec(const char* EigFile, const char* TpsFile, const char* IndexFile,
                         const char* AverageFile)
{
    read_eigenvectors(EigFile, &this->natoms, &this->bFit,
                      fca_utils::make_auto_move_ptr(this->xref),
                      &this->bDMR, fca_utils::make_auto_move_ptr(this->xav),
                      &this->bDMA, &this->nvec,
                      fca_utils::make_auto_move_ptr(this->eignr), &this->vec,
                      fca_utils::make_auto_move_ptr(this->val));

    gmx_bool bM = TRUE;
    if ((!bFit || this->xref) && !bDMR && !bDMA)
    {
        bM = FALSE;
    }

    gmx_bool bTPS = TRUE; // from g_anaeig
    if (this->xref == nullptr)
    {
        bTPS = TRUE;
    }

    this->nfit  = 0;
    gmx_bool bTop = TRUE;
    t_atoms  atoms;
    atoms.nr = 0;

    if (!bTPS)
    {
        bTop = FALSE;
    }
    else
    {
        get_stx_coordnum(TpsFile, &natoms);

        this->xtop.reset(new rvec[natoms]());
        init_t_atoms(&atoms, natoms, TRUE);

        char   title[STRLEN] = "";
        matrix topbox;
        int    ePBC;
        read_stx_conf(TpsFile, nullptr, (char**)&title, &atoms, this->xtop.get(), nullptr, &ePBC,
                      topbox);

        rm_gropbc(&atoms, this->xtop.get(), topbox);

        /* Fitting is only needed when we need to read a trajectory */
        char* grpname      = nullptr;

        if (this->xref == nullptr)
        {
            if (bFit)
            {
                printf("\nNote: the structure in %s should be the same\n"
                       "      as the one used for the fit in g_covar\n",
                       TpsFile);
                printf( "\nSelect the index group that was used for the least squares fit in g_covar\n");
                get_index((&atoms), IndexFile, 1, &this->nfit, fca_utils::make_auto_move_ptr(this->ifit),
                          &grpname);
                this->w_rls.reset(new real[(&atoms)->nr]());
                if (bM && bDMR)
                {
                    for (int i = 0; (i < this->nfit); i++)
                    {
                        this->w_rls[this->ifit[i]] = (&atoms)->atom[this->ifit[i]].m;
                    }
                }
                else
                {
                    for (int i = 0; (i < this->nfit); i++)
                    {
                        this->w_rls[this->ifit[i]] = 1.0;
                    }
                }
            }
            else /*no fit*/
            {
                fprintf(stdout,
                        "\
                            *********************************************************************************\n\
                            no fit was used in g_covar, this seems odd... and we normally do not support such behaviour\n\
                            ------------- however, we assume for now that your trajectory has already been fitted and proceed ------------\n\
                            *********************************************************************************************\n");
                this->nfit  = 0;
                this->ifit  = nullptr;
                this->w_rls = nullptr;
            }
        } // now we want to first load the this->index ---> and then go thru the *have xref* code

        int i;
        get_index((&atoms), IndexFile, 1, &i, fca_utils::make_auto_move_ptr(this->index), &grpname);

        if (this->xref)
        {   /* we have xref */
            /* make the fit index in xref instead of xtop */
            this->nfit = this->natoms;
            int* new_ifit;
            snew(new_ifit, this->nfit);
            this->ifit.reset(new_ifit);
            this->w_rls.reset(new real[this->nfit]);
            for (int i = 0; (i < this->nfit); i++)
            {
                this->ifit[i] = i;
                if (bDMR && bM)
                {
                    this->w_rls[i] = (&atoms)->atom[this->index[i]].m;
                }
                else
                {
                    this->w_rls[i] = 1.0;
                }
            }
        }
    }

    if (AverageFile && bTop)
    {
        fca_utils::InputFunctions::get_structure((&atoms), IndexFile, AverageFile, this->xav.get(),
                                                 this->natoms, this->index.get());
    }

    this->sqrtm.reset(new real[this->natoms]);
    if (bM && bDMA)
    {
        for (int i = 0; (i < this->natoms); i++)
        {
            this->sqrtm[i] = sqrt((&atoms)->atom[this->index[i]].m);
        }
    }
    else
    {
        for (int i = 0; (i < this->natoms); i++)
        {
            this->sqrtm[i] = 1.0;
        }
    }

    if (this->xref)
    {
        this->ned = this->index[this->natoms - 1] + 1;
    }
    else
    {
        int nfit = 0;
        if (this->ifit)
        {
            nfit = this->ifit[this->nfit - 1] + 1;
        }
        this->ned = fca_utils::maxv(this->index[this->natoms - 1] + 1, nfit);
    }
} //init_eigenvectors

FcaEigenVec::~FcaEigenVec()
{
    for (int idxvec = 0; idxvec < nvec; ++idxvec)
    {
        sfree(vec[idxvec]);
    }
    sfree(vec);
}

void FcaEigenVec::setVal(const real newVal[], const int nbVal)
{
    val.reset(new real[nbVal]());
    for (int idx = 0; idx < nbVal; ++idx)
    {
        val[idx] = newVal[idx];
    }
}

void FcaEigenVec::fitit(const rvec* xread, rvec* x, rvec* transvec, matrix &rotmat) const
{
    if (xref)
    {
        for (int i = 0; i < natoms; i++)
        {
            copy_rvec(xread[index[i]], x[i]);
        }
    }
    if (bFit)
    {
        if (xref)
        {
            internal_fitit(natoms, x, w_rls.get(), transvec, rotmat,
                           nfit, ifit.get(), xref.get());
        }
        else
        {
            internal_fitit(ned, const_cast<rvec*>(xread), w_rls.get(), transvec, rotmat, // TODO
                           nfit, ifit.get(), xtop.get());
        }
    }
    if (!xref)
    {
        for (int i = 0; i < natoms; i++)
        {
            copy_rvec(xread[index[i]], x[i]);
        }
    }
}

void FcaEigenVec::fitit(rvec* xread, rvec* transvec, matrix &rotmat) const
{
    std::unique_ptr<rvec[]> x(new rvec[natoms]());
    fitit(xread, x.get(), transvec, rotmat);
}

std::unique_ptr<real[]> FcaEigenVec::project_frame(rvec* x, rvec* v, real* vproj,
                                                   rvec* f, real* fproj, int natoms,
                                                   matrix box, const real t,
                                                   const int nout, const int outvec[]) const
{
    fca_utils::rm_pbc_plain(natoms, x, box);

    std::unique_ptr<rvec[]> old_transvec(new rvec[ned]());
    std::unique_ptr<rvec[]> older_transvec(new rvec[ned]());
    std::unique_ptr<rvec[]> old_vcorr(new rvec[ned]());
    matrix                  old_rotmat;

    return project_vector(x, v, vproj, f, fproj, t, nout, outvec,
                          old_transvec.get(), older_transvec.get(),
                          old_vcorr.get(), &old_rotmat, true);
}

std::vector< std::unique_ptr< real[] > > FcaEigenVec::project_vectors(const std::vector< std::unique_ptr< rvec[] > > &xread, rvec* velread,
                                                                      real* vproj, rvec* fread, real* fproj, const real input_dt,
                                                                      const int nout, const int outvec[]) const
{
    std::vector< std::unique_ptr< real[] > > result;

    std::unique_ptr<rvec[]>                  old_transvec(new rvec[ned]());
    std::unique_ptr<rvec[]>                  older_transvec(new rvec[ned]());
    std::unique_ptr<rvec[]>                  old_vcorr(new rvec[ned]());
    matrix old_rotmat;

    bool   resetBuffer = true;

    for (const auto &vec : xread)
    {
        result.emplace_back(project_vector(vec.get(), velread,
                                           vproj, fread, fproj, input_dt,
                                           nout, outvec, old_transvec.get(),
                                           older_transvec.get(), old_vcorr.get(),
                                           &old_rotmat, resetBuffer));
        // reset buffer only the first loop
        resetBuffer = false;
    }
    return result;
}


std::unique_ptr<real[]> FcaEigenVec::project_vector(const rvec* xread, rvec* velread,
                                                    real* vproj, rvec* fread, real* fproj, const real input_dt,
                                                    const int nout, const int outvec[],
                                                    rvec old_transvec[], rvec older_transvec[],
                                                    rvec old_vcorr[], matrix* old_rotmat,
                                                    const bool resetBuffer) const
{
    const real              dt          = fca_utils::maxv(input_dt, real(1));
    const gmx_bool          bCorrection = (vproj || fproj);
    std::unique_ptr<rvec[]> vel;
    std::unique_ptr<rvec[]> force;
    std::unique_ptr<rvec[]> xdum;
    /*------------------------- start computation ----------------------*/
    /* get everything in sequence of just projection atoms */
    if (bCorrection)
    {
        vel.reset(new rvec[ned]());
        xdum.reset(new rvec[ned]());

        assert(velread && xread);
        for (int i = 0; i < natoms; i++)
        {
            copy_rvec(velread[index[i]], vel[i]);
            copy_rvec(xread[index[i]], xdum[i]); /* get unfitted copy to rotate with alternative rotation matrix */
        }

        if (fread)
        {
            force.reset(new rvec[ned]());
            for (int i = 0; i < natoms; i++)
            {
                copy_rvec(fread[index[i]], force[i]);
            }
        }
    }

    matrix                  rotmat;
    std::unique_ptr<rvec[]> transvec(new rvec[ned]());
    std::unique_ptr<rvec[]> x(new rvec[natoms]());
    fitit(xread, x.get(), transvec.get(), rotmat);
    /* now x is only analysis atoms and fitted, xread might be fitted or not*/

    /* rotate velocities */
    if (bCorrection && bFit)
    {
        assert(vel != nullptr);
        if (resetBuffer)
        {
            copy_rvecn(transvec.get(), old_transvec, 0, ned);
            copy_mat(rotmat, (*old_rotmat));
            if (force)
            {
                copy_rvecn(old_transvec, older_transvec, 0, ned);
            }
        }
        std::unique_ptr<rvec[]> vcorr(new rvec[ned]());
        fca_utils::correct_vel(natoms, vel.get(), transvec.get(), old_transvec,
                               -1.0 * dt, vcorr.get());
        if (force)
        {   /* correct forces to incorporate effect of corrected velocities */
            fca_utils::correct_force(natoms, sqrtm.get(), force.get(), transvec.get(),
                                     old_transvec, -1.0 * dt * dt);
            fca_utils::correct_force(natoms, sqrtm.get(), force.get(), old_transvec,
                                     older_transvec, 1.0 * dt * dt); /* opposite sign's here are correct */
            copy_rvecn(old_transvec, older_transvec, 0, natoms);
        }
        copy_rvecn(transvec.get(), old_transvec, 0, natoms);

        /*correct for rotational motion */
        std::unique_ptr<rvec[]> vdum(new rvec[ned]());
        copy_rvecn(vel.get(), vdum.get(), 0, natoms);
        fca_utils::rotate_vec(natoms, vdum.get(), (*old_rotmat)); /* two copies of velocities rotated with different matrices */
        fca_utils::rotate_vec(natoms, vel.get(), rotmat);         /* to get rotational motion in velocities -- > correct forces */

        if (force)
        {
            fca_utils::rotate_vec(natoms, force.get(), rotmat);                                           /* rotate forces before you correct them */
            fca_utils::correct_force(natoms, sqrtm.get(), force.get(), vel.get(), vdum.get(), -1.0 * dt); /* correct with rotation in velocities */
        }

        fca_utils::rvecadd(natoms, xdum.get(), transvec.get(), xdum.get()); /*move copy of actual structure into origin */
        fca_utils::rotate_vec(natoms, xdum.get(), (*old_rotmat));           /* rotate positions with old matrix -> rotation in positions can be extracted */
        fca_utils::correct_vel(natoms, vel.get(), x.get(), xdum.get(), -1.0 * dt, vcorr.get());

        if (!resetBuffer && force)
        {
            fca_utils::correct_force(natoms, sqrtm.get(), force.get(), vcorr.get(), old_vcorr, -1.0 * dt); /* add some forces to account for corrected velocities */
        }
        if (force)
        {
            copy_rvecn(vcorr.get(), old_vcorr, 0, natoms);
        }
        copy_mat(rotmat, (*old_rotmat));
    }

    std::unique_ptr<real[]> xproj(new real[nout]());
    /*  calculate (mass-weighted) projection */
    for (int v = 0; v < nout; v++)
    {
        const int vecidx = outvec[v];
        real      xinp   = 0;
        for (int i = 0; i < natoms; i++)
        {
            xinp += (vec[vecidx][i][0] * (x[i][0] - xav[i][0]) + vec[vecidx][i][1] * (x[i][1] - xav[i][1]) + vec[vecidx][i][2] * (x[i][2] - xav[i][2])) * sqrtm[i];
        }
        xproj[v] = xinp;
    }

    if (bCorrection && vel)
    {
        assert(vproj);
        for (int v = 0; v < nout; v++)
        {
            const int vecidx = outvec[v];
            real      vinp   = 0;
            for (int i = 0; i < natoms; i++)
            {
                if (bCorrection)
                {
                    vinp += (vec[vecidx][i][0] * (vel[i][0]) + vec[vecidx][i][1] * (vel[i][1]) + vec[vecidx][i][2] * (vel[i][2])) * sqrtm[i];
                }
            }
            vproj[v] = vinp;
        }
    }

    if (bCorrection && force)
    {
        assert(fproj);
        for (int v = 0; v < nout; v++)
        {
            const int vecidx = outvec[v];
            real      finp   = 0;
            for (int i = 0; i < natoms; i++)
            {
                finp += (vec[vecidx][i][0] * (force[i][0]) + vec[vecidx][i][1] * (force[i][1]) + vec[vecidx][i][2] * (force[i][2])) / sqrtm[i];
            }
            fproj[v] = finp;
        }
    }

    return xproj;
}

std::unique_ptr< real[] > FcaEigenVec::produce_new_eigvecs(const std::unique_ptr< real[] >* ic_vec, const int dim) const
{
    const int                 ndim = natoms * DIM;
    std::unique_ptr< real[] > mat(new real[nvec * ndim]());
#if GMX_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int k = 0; k < nvec; k++)
    {
        for (int j = 0; j < natoms; j++)
        {
            for (int d = 0; d < DIM; d++)
            {
                if (k < dim)
                {
                    mat[k * ndim + j * DIM + d] = 0; //< run multiplication of v_ijd w_ik
                    for (int i = 0; i < nvec; i++)
                    {
                        if (i < dim)
                        {
                            mat[k * ndim + j * DIM + d] += vec[i][j][d] * ic_vec[k][i];
                        }
                    }
                }
                else
                {
                    mat[k * ndim + j * DIM + d] = vec[k][j][d]; //< mode number higher then fca.dim, just copy PCA modes
                }
            }
        }
    }
    return mat;
}

} //gmx namespace
