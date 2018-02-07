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
/*! \internal \file
 * \brief
 * Declares FcaEigenVec.
 */
#ifndef FCA_EIGENVEC_H
#define FCA_EIGENVEC_H

#include <memory>
#include <vector>

#include "gromacs/gmxana/eigio.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "utils.h"

namespace gmx
{

class FcaEigenVec
{
    int    nvec;                                                     /* number eigenvectors saved in vec */
    rvec** vec;                                                      /* eigenvectors */
    std::unique_ptr < int[], fca_utils::sfree_deleter < int>> eignr; /* indices of eigenvectors */
    int    ned;                                                      /* highest atomindex used */
    int    natoms;                                                   /* number of atoms used in PCA */
    std::unique_ptr < int[], fca_utils::sfree_deleter < int>> index; /* indices of atoms used in PCA */
    std::unique_ptr < rvec[], fca_utils::sfree_deleter < rvec>> xav;
    std::unique_ptr < rvec[], fca_utils::sfree_deleter < rvec>> xref;
    std::unique_ptr<rvec[]> xtop;
    gmx_bool                bDMA;
    gmx_bool                bDMR;
    gmx_bool                bFit;
    std::unique_ptr<real[]> w_rls;
    std::unique_ptr<real[]> sqrtm;
    int                     nfit;                                   /* number of atoms used in fit */
    std::unique_ptr < int[], fca_utils::sfree_deleter < int>> ifit; /* indices of atoms used for fit */
    std::unique_ptr<rvec[]> x;                                      /* after project this is the fitted structure */
    std::unique_ptr < real[], fca_utils::sfree_deleter < real>> val;


    static void internal_fitit(const int nr, rvec* x, const real* w_rls, rvec* transvec,
                               matrix &rmat, const int nfit, const int* ifit, const rvec* xref);

    public:
        FcaEigenVec(const char* EigFile, const char* TpsFile, const char* IndexFile,
                    const char* AverageFile);

        ~FcaEigenVec();

        int getNbAtoms() const { return natoms; }

        int getNbVec() const{ return nvec; }

        const int* geteignr() const{ return eignr.get(); }

        void setVal(const real newVal[], const int nbVal);

        void fitit(const rvec* xread, rvec* x, rvec* transvec, matrix &rotmat) const;

        void fitit(rvec* xread, rvec* transvec, matrix &rotmat) const;

        std::unique_ptr<real[]> project_frame(rvec* x, rvec* v, real* vproj,
                                              rvec* f, real* fproj, int natoms,
                                              matrix box, const real t,
                                              const int nout, const int outvec[]) const;

        std::vector< std::unique_ptr< real[] > > project_vectors(const std::vector< std::unique_ptr< rvec[] > > &xread, rvec* velread,
                                                                 real* vproj, rvec* fread, real* fproj, const real input_dt,
                                                                 const int nout, const int outvec[]) const;


        std::unique_ptr<real[]> project_vector(const rvec* xread, rvec* velread,
                                               real* vproj, rvec* fread, real* fproj,
                                               const real input_dt, const int nout, const int outvec[],
                                               rvec old_transvec[], rvec older_transvec[],
                                               rvec old_vcorr[], matrix* old_rotmat,
                                               const bool resetBuffer) const;

        /*! \brief
         * FCA-mode k is generetad from multiplication of x_jd with PCA eigenvector matrix v_ijd
         * and subsequently rotating the first fca.dim PCA-modes using w_ik.
         * this matrix multiplication is performed as y_k= x_jd v_ijd w_ik
         * output is sum_i v_ijd w_ik. All k>fca.dim are the old PCA modes i=k.
         */
        std::unique_ptr< real[] > produce_new_eigvecs(const std::unique_ptr< real[] >* ic_vec, const int dim) const;


        void writeEigvecsWithMat(const real mat[], const char* EigFile)
        {
            write_eigenvectors(EigFile, getNbAtoms(), mat, FALSE, 1, nvec,
                               (xref != nullptr ? eWXR_YES : eWXR_NOFIT),
                               xref.get(), bDMR, xav.get(), bDMA, val.get());
        }
};

} //gmx namespace

#endif
