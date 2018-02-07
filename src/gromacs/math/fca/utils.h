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
 * Implements utility functions and classes for Fca.
 */
#ifndef FCA_UTILS_H
#define FCA_UTILS_H

#define FCA_MASTER_RANK 0

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctime>

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if GMX_OPENMP
#include <omp.h>
#endif

#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

namespace fca_utils
{

template <class NumType>
int fscanf_floating_point(FILE* fl, NumType* ptr);

template <> inline
int fscanf_floating_point<float>(FILE* fl, float* ptr) { return fscanf(fl, "%f", ptr); }

template <> inline
int fscanf_floating_point<double>(FILE* fl, double* ptr) { return fscanf(fl, "%lf", ptr); }


inline static double PI_constant() { return 3.14159265358979323846; }

/* MPI Wrapper section */

#if GMX_LIB_MPI
inline void AssertMpi(const int returnedValue, const std::string &filename, const int lineNumber)
{
    if (MPI_SUCCESS != returnedValue)
    {
        std::string messageErrorMpi = "MPI Error at line " + std::to_string(lineNumber) + " in file "
            + filename + ", value "  + std::to_string(returnedValue);
        throw std::runtime_error(messageErrorMpi);
    }
}
#endif

class mpi
{
    int my_rank;
    int num_nodes;

    public:
        mpi()
        {
#if GMX_LIB_MPI
            my_rank   = gmx_node_rank();
            num_nodes = gmx_node_num();
#else
            my_rank   = 0;
            num_nodes = 1;
#endif
        }

        virtual ~mpi() { }

        mpi(const mpi &)            = delete;
        mpi &operator=(const mpi &) = delete;

        bool isMaster() const { return ((my_rank == FCA_MASTER_RANK)); }

        bool isPar() const { return ((num_nodes > 1)); }

        int getMyRank() const { return my_rank; }

        int getNumNodes() const { return num_nodes; }
};

/* Vec/mat functions section */

/*! \brief
 * Check periodic boundary
 */
inline void rm_pbc_plain(const int natoms, rvec* x, const matrix box)
{
    for (int n = 1; n < natoms; n++)
    {
        for (int m = DIM - 1; m >= 0; m--)
        {
            real dist = x[n][m] - x[n - 1][m];
            if (fabs(dist) > 0.9 * box[m][m])
            {
                if (dist > 0)
                {
                    for (int d = 0; d <= m; d++)
                    {
                        x[n][d] -= box[m][d];
                    }
                }
                else
                {
                    for (int d = 0; d <= m; d++)
                    {
                        x[n][d] += box[m][d];
                    }
                }
            }
        }
    }
}

inline void apply_rm_pbc_plain_to_all(const int natoms, std::vector < std::unique_ptr < rvec[]>>* vecs, const matrix &box)
{
#if GMX_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int idxVec = 0; idxVec < int(vecs->size()); ++idxVec)
    {
        rm_pbc_plain(natoms, (*vecs)[idxVec].get(), box);
    }
}


/*! \brief
 * Rotate X
 */
inline void rotate_vec(const int natoms, rvec* x, const matrix &R)
{
    for (int j = 0; j < natoms; j++)
    {
        rvec x_old;
        for (int m = 0; m < DIM; m++)
        {
            x_old[m] = x[j][m];
        }
        for (int r = 0; r < DIM; r++)
        {
            x[j][r] = 0;
            for (int c = 0; c < DIM; c++)
            {
                x[j][r] += R[r][c] * x_old[c];
            }
        }
    }
}

inline void correct_vel(const int nr, rvec* vel, const rvec* v1, const rvec* v2, const real dt,
                        rvec* vcorr)
{
    const real fact = 1.0 / dt;
    for (int i = 0; i < nr; i++)
    {
        for (int m = 0; m < DIM; m++)     /* x3-2*x2+x1/(dt^2*m) is force */
        {
            vcorr[i][m] = fact * v2[i][m] - fact * v1[i][m];
            vel[i][m]  += vcorr[i][m];
        }
    }
}

inline void correct_force(const int nr, const real* eig_sqrtm, rvec* force, const rvec* v1, const rvec* v2,
                          const real dt)
{
    const real fact = 1.0 / dt;
    for (int i = 0; i < nr; i++)
    {
        for (int m = 0; m < DIM; m++)     /* x3-2*x2+x1/(dt^2*m) is force */
        {
            force[i][m] += eig_sqrtm[i] * eig_sqrtm[i] * (fact * v2[i][m] - fact * v1[i][m]);
        }
    }
}

template <class ObjType>
inline ObjType squareof(const ObjType &obj) { return obj*obj; }

/* Ascii eig vev output section */

inline void write_eigvecs_ascii_fca(real mat[], const int natoms, const int nvec, char* fname)
{
    FILE    * fp = gmx_ffopen(fname, "w");
    assert(fp);
    const int ndim = natoms * DIM;
    for (int i = 0; i < ndim; i++)
    {
        for (int j = 0; j < nvec; j++)
        {
            fprintf(fp, "%g ", mat[j * ndim + i]);
        }
        fprintf(fp, "\n");
    }
}

/*! \brief
 * rvec extra method: c=a+b
 */
inline void rvecadd(const int dim, const rvec a[], const rvec b[], rvec c[])
{
    for (int i = 0; i < dim; i++)
    {
        rvec_add(a[i], b[i], c[i]);
    }
}

/* Ascii matrix output section */

inline void dump_sqrmatrix(FILE* out, const double* data, const int dim)
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            fprintf(out, "%e ", data[i + dim * j]);
        }
        fprintf(out, "\n");
    }
}

inline void dump_sqrmatrixr(FILE* out, const real* data, const int dim)
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            fprintf(out, "%e ", data[i + dim * j]);
        }
        fprintf(out, "\n");
    }
}

/* Ascii eig vev output section */

template <class T>
struct sfree_deleter
{
    void operator()(T* b) { sfree(b); }
};

template <class T>
class auto_move_ptr
{
    std::unique_ptr < T[], sfree_deleter < T>>&uptr;
    T* temp_ptr;
    public:
        explicit auto_move_ptr(std::unique_ptr < T[], sfree_deleter < T>> &&in_uptr)
            : uptr(in_uptr), temp_ptr(nullptr)
        {
        }

        auto_move_ptr(const auto_move_ptr &)            = delete;
        auto_move_ptr &operator=(const auto_move_ptr &) = delete;

        auto_move_ptr(auto_move_ptr &&)            = default;
        auto_move_ptr &operator=(auto_move_ptr &&) = default;

        ~auto_move_ptr()
        {
            if (temp_ptr)
            {
                uptr.reset(temp_ptr);
            }
        }

        operator T**(){
            return &temp_ptr;
        }
};

template <class T>
auto_move_ptr<T> make_auto_move_ptr(std::unique_ptr < T[], sfree_deleter < T>> &ptr)
{
    return auto_move_ptr<T>(std::move(ptr));
}

/* MinMax section */

template <class Type>
inline constexpr Type limits_min()
{
    return std::numeric_limits<Type>::min();
}

template <class Type>
inline constexpr Type limits_max()
{
    return std::numeric_limits<Type>::max();
}

template <class Type>
inline constexpr const Type &minv(const Type &v1, const Type &v2)
{
    return v1 < v2 ? v1 : v2;
}

template <class Type>
inline constexpr const Type &maxv(const Type &v1, const Type &v2)
{
    return v1 > v2 ? v1 : v2;
}

/*! \internal \brief
 * Implements Log structure.
 */
struct Log
{
    Log() : MI(nullptr),
            ic_matrix(nullptr),
            move_log(nullptr),
            move_matrix(nullptr),
            icx(nullptr),
            MI_matrix(nullptr)
    {
    }
    ~Log()
    {
        if (MI)
        {
            fclose(MI);
        }
        if (ic_matrix)
        {
            fclose(ic_matrix);
        }
        if (move_log)
        {
            fclose(move_log);
        }
        if (move_matrix)
        {
            fclose(move_matrix);
        }
        if (icx)
        {
            fclose(icx);
        }
        if (MI_matrix)
        {
            fclose(MI_matrix);
        }
    }

    Log(const Log &)            = delete;
    Log &operator=(const Log &) = delete;

    FILE* MI;
    FILE* ic_matrix;
    FILE* move_log;
    FILE* move_matrix;
    FILE* icx;
    FILE* MI_matrix;
};

/*! \internal \brief
 * Implements TrxInfo class.
 */
class TrxInfo
{
    const char* trj_file;
    int         dimpca;
    const char* pca_file;

    public:
        TrxInfo() : trj_file(nullptr), dimpca(0), pca_file(nullptr) { }

        void setTrj(const char* inTrj) { trj_file = inTrj; }

        void setPca(const char* inPca) { pca_file = inPca; }
        void setDimPca(const int inDimPca) { dimpca = inDimPca; }

        template <class FcaMasterClass>
        void write_fcalog(const char* fname, const FcaMasterClass &fca,
                          const char* posfile, const char* eigvecfile)
        {
            FILE * out;
            time_t now;
            out = gmx_ffopen(fname, "w");
            now = time(nullptr);
            fprintf(out, "FCA - Analysis log, written %s\n", ctime(&now));

            if (trj_file)
            {
                fprintf(out,
                        "Read %d principal components for intial projection from %s\n",
                        dimpca, pca_file);
            }
            else
            {
                fprintf(out, "Read data from %s\n", posfile);
            }

            fprintf(out, "\n");
            fprintf(out, "did FCA-Minimization in %d dimensions \n", fca.getDim());
            fprintf(out, "finally reached MI-sum: %e \n", fca.getMIsum());
            if (eigvecfile)
            {
                fprintf(out, "Wrote fca-vectors to %s\n", eigvecfile);
            }
            fclose(out);
        }
};

/*! \internal \brief
 * Implements parallel (multi thread) sort class.
 */
class PSort
{
#if GMX_OPENMP && _OPENMP >= 200805
    template <class ArrayType, class CompareType>
    static ArrayType selectPivot(ArrayType begin, ArrayType end, CompareType &&comparator)
    {
        ArrayType middle = std::next(begin, std::distance(begin, end)/2);

        std::advance(end, -1);

        if ((comparator(*begin, *middle) && comparator(*middle, *end))
            || (comparator(*end, *middle) && comparator(*middle, *begin)))
        {
            return middle;
        }
        else if ((comparator(*middle, *begin) && comparator(*begin, *end))
                 || (comparator(*end, *begin) && comparator(*begin, *middle)))
        {
            return begin;
        }
        else
        {
            return end;
        }
    }

    template <class ArrayType, class CompareType>
    static void sortCore(ArrayType begin, ArrayType end, CompareType &&comparator, const int deepTask)
    {
        if (std::distance(begin, end) == 0)
        {
            return;
        }

        if (deepTask == 0)
        {
            std::sort(begin, end, std::forward<CompareType>(comparator));
        }
        else
        {
            ArrayType pivot = selectPivot(begin, end, std::forward<CompareType>(comparator));
            using ValueType = decltype(*pivot);
            auto      it = std::partition(begin, end, [pivot, &comparator](const ValueType &item){return comparator(item, *pivot); });

#pragma omp task default(shared)
            sortCore(begin, it, std::forward<CompareType>(comparator), deepTask-1);

            sortCore(it, end, std::forward<CompareType>(comparator), deepTask-1);

#pragma omp taskwait
        }
    }

    public:
        template <class ArrayType, class CompareType>
        static void sort(ArrayType begin, ArrayType end, CompareType &&comparator)
        {
            const int nbTasks = std::min(static_cast<decltype(std::distance(begin, end)/512)>(omp_get_max_threads()),
                                         std::distance(begin, end)/512);
            int       deepTask = 1;
            while ((1 << deepTask) < nbTasks) { deepTask += 1; }

#pragma omp parallel default(shared)
#pragma omp master
            sortCore(begin, end, std::forward<CompareType>(comparator), deepTask);
        }
#else
    public:
        template <class ArrayType, class CompareType>
        static void sort(ArrayType begin, ArrayType end, CompareType &&comparator)
        {
            std::sort(begin, end, std::forward<CompareType>(comparator));
        }
#endif

        template <class ArrayType>
        static void sort(ArrayType begin, ArrayType end)
        {
            using ValueType = decltype(*begin);
            std::sort(begin, end, [](const ValueType &v1, const ValueType &v2){return v1 < v2; });
        }
};

/*! \internal \brief
 * Declares input functions.
 */
class InputFunctions
{
    /*! \brief
     * File loading function, read coordinates out of STX file.
     */
    static int read_conffile(const char* confin, char* title, rvec* x[]);

    /*! \brief
     * Copies coordinates from x to edx which are given in index
     */
    static void filter2x(rvec* x, int nindex, int index[], int ngro,
                         int igro[], rvec* xori, const char* structure);

    /*! \brief
     * Reading frames.
     */
    static int count_number_columns(FILE* fp);

    public:

        static bool get_structure(t_atoms* atoms, const char* IndexFile,
                                  const char* StructureFile, rvec* x,
                                  int natoms, int index[]);

        static std::vector< std::unique_ptr< real[] > > read_fca_proj(const char pos_file[],
                                                                      int      * dim);
};

} //fca_utils namespace

} //gmx namespace

#endif
