#ifndef UTILS_HPP
#define UTILS_HPP

#ifdef GMX_MPI
#define MPICC
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#else
#define MPI_Bcast(A, B, C, D, E)
#define MPI_COMM_WORLD 0
#endif
#define MASTER_RANK 0

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#include <limits>
#include <algorithm>
#include <memory>

#include <vector>
#include <cstdio>
#include <cassert>
#include <unistd.h>

#include "gmxextern.hpp"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/cstringutil.h"

namespace FCA{

namespace utils{

/////////////////////////////////////////////////////////////////////////////////////
/// MPI Wrapper
/////////////////////////////////////////////////////////////////////////////////////

#ifndef AssertMpi
#define AssertMpi(X) if(MPI_SUCCESS != (X)) { printf("MPI Error at line %d\n",__LINE__); fflush(stdout) ; throw std::runtime_error("Stop from from mpi erro"); }
#endif

class mpi {
    int my_rank;
    int num_nodes;

public:
#ifndef MPICC
    mpi(int* /*argc*/, char*** /*argv*/) {
        my_rank   = 0;
        num_nodes = 1;
    }

    virtual ~mpi() {
    }
#else
    mpi(int* argc, char*** argv) {
        AssertMpi(MPI_Init(argc, argv)); /*  already includes MPE_Init_log(); */
        AssertMpi(MPI_Comm_size(MPI_COMM_WORLD, &num_nodes));
        AssertMpi(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));
        /*fprintf(stderr,"Initialized MPI on node %d for %d nodes\n",  mpi->my_rank,mpi->num_nodes);*/
    }

    virtual ~mpi() {
        AssertMpi(MPI_Finalize());
    }
#endif

    mpi(const mpi&) = delete;
    mpi& operator=(const mpi&) = delete;

    bool isMaster() const {
        return ((my_rank == MASTER_RANK));
    }

    bool isPar() const {
        return ((num_nodes > 1));
    }

    int getMyRank() const {
        return my_rank;
    }

    int getNumNodes() const {
        return num_nodes;
    }
};

/////////////////////////////////////////////////////////////////////////////////////
/// Vec/mat functions
/////////////////////////////////////////////////////////////////////////////////////

inline void rm_pbc_plain(const int natoms, rvec* x, const matrix box) {
    /* check periodic boundary */
    for(int n = 1; n < natoms; n++) {
        for(int m = DIM - 1; m >= 0; m--) {
            real dist = x[n][m] - x[n - 1][m];
            if(fabs(dist) > 0.9 * box[m][m]) {
                if(dist > 0) {
                    for(int d = 0; d <= m; d++) {
                        x[n][d] -= box[m][d];
                    }
                } else {
                    for(int d = 0; d <= m; d++) {
                        x[n][d] += box[m][d];
                    }
                }
            }
        }
    }
}

inline void apply_rm_pbc_plain_to_all(const int natoms, std::vector<std::unique_ptr<rvec[]>>* vecs, const matrix& box){
    for(size_t idxVec = 0 ; idxVec < vecs->size() ; ++idxVec){
        rm_pbc_plain(natoms, (*vecs)[idxVec].get(), box);
    }
}


inline void rotate_vec(const int natoms, rvec* x, const matrix& R) {
    /*rotate X*/
    for(int j = 0; j < natoms; j++) {
        rvec x_old;
        for(int m = 0; m < DIM; m++){
            x_old[m] = x[j][m];
        }
        for(int r = 0; r < DIM; r++) {
            x[j][r] = 0;
            for(int c = 0; c < DIM; c++){
                x[j][r] += R[r][c] * x_old[c];
            }
        }
    }
}

inline void correct_vel(const int nr, rvec* vel, const rvec* v1, const rvec* v2, const real dt,
                        rvec* vcorr) {
    const real fact = 1.0 / dt;
    for(int i = 0; i < nr; i++) {
        for(int m = 0; m < DIM; m++) /* x3-2*x2+x1/(dt^2*m) is force */ {
            vcorr[i][m] = fact * v2[i][m] - fact * v1[i][m];
            vel[i][m] += vcorr[i][m];
        }
    }
}

inline void correct_force(const int nr, const real* eig_sqrtm, rvec* force, const rvec* v1, const rvec* v2,
                          const real dt) {
    const real fact = 1.0 / dt;
    for(int i = 0; i < nr; i++){
        for(int m = 0; m < DIM; m++){ /* x3-2*x2+x1/(dt^2*m) is force */
            force[i][m] += eig_sqrtm[i] * eig_sqrtm[i] * (fact * v2[i][m] - fact * v1[i][m]);
        }
    }
}

template <class ObjType>
inline ObjType squareof(const ObjType& obj){
    return obj*obj;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Ascii eig vev output
/////////////////////////////////////////////////////////////////////////////////////

inline void write_eigvecs_ascii_fca(real mat[], const int natoms, const int nvec, char* fname) {
    FILE* fp = gmx_ffopen(fname, "w");
    assert(fp);
    const int ndim = natoms * DIM;
    for(int i = 0; i < ndim; i++) {
        for(int j = 0; j < nvec; j++){
            fprintf(fp, "%g ", mat[j * ndim + i]);
        }
        fprintf(fp, "\n");
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/// rvec extra method
/////////////////////////////////////////////////////////////////////////////////////

inline void rvecadd(const int dim, const rvec a[], const rvec b[], rvec c[]) {
    //c=a+b;
    for(int i = 0; i < dim; i++) {
        rvec_add(a[i], b[i], c[i]);
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/// Ascii matrix output
/////////////////////////////////////////////////////////////////////////////////////

inline void dump_sqrmatrix(FILE* out, const double* data, const int dim) {
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            fprintf(out, "%10.5f", data[i + dim * j]);
        }
        fprintf(out, "\n");
    }
}

inline void dump_sqrmatrixr(FILE* out, const real* data, const int dim) {
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            fprintf(out, "%10.5f", data[i + dim * j]);
        }
        fprintf(out, "\n");
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/// Ascii eig vev output
/////////////////////////////////////////////////////////////////////////////////////

template <class T>
struct sfree_deleter {
    void operator()(T* b) { sfree(b); }
};

template <class T>
class auto_move_ptr{
    std::unique_ptr<T[], sfree_deleter<T>>& uptr;
    T* temp_ptr;
public:
    explicit auto_move_ptr(std::unique_ptr<T[], sfree_deleter<T>>& in_uptr)
        : uptr(in_uptr), temp_ptr(nullptr){
    }

    ~auto_move_ptr(){
        if(temp_ptr){
            uptr.reset(temp_ptr);
        }
    }

    operator T**(){
        return &temp_ptr;
    }
};

template <class T>
auto_move_ptr<T> make_auto_move_ptr(std::unique_ptr<T[], sfree_deleter<T>>& ptr){
    return auto_move_ptr<T>(ptr);
}


/////////////////////////////////////////////////////////////////////////////////////
/// MinMax
/////////////////////////////////////////////////////////////////////////////////////

template <class Type>
inline constexpr Type limits_min(){
    return std::numeric_limits<Type>::min();
}

template <class Type>
inline constexpr Type limits_max(){
    return std::numeric_limits<Type>::max();
}

template <class Type>
inline constexpr const Type& minv(const Type& v1, const Type& v2){
    return v1 < v2 ? v1 : v2;
}

template <class Type>
inline constexpr const Type& maxv(const Type& v1, const Type& v2){
    return v1 > v2 ? v1 : v2;
}


/////////////////////////////////////////////////////////////////////////////////////
/// Log
/////////////////////////////////////////////////////////////////////////////////////

struct Log {
    Log(): MI(nullptr),
        ic_matrix(nullptr),
        move_log(nullptr),
        move_matrix(nullptr),
        icx(nullptr),
        MI_matrix(nullptr){
    }
    ~Log(){
        if(MI){
            fclose(MI);
        }
        if(ic_matrix){
            fclose(ic_matrix);
        }
        if(move_log){
            fclose(move_log);
        }
        if(move_matrix){
            fclose(move_matrix);
        }
        if(icx){
            fclose(icx);
        }
        if(MI_matrix){
            fclose(MI_matrix);
        }
    }

    Log(const Log&) = delete;
    Log& operator=(const Log&) = delete;

    FILE* MI;
    FILE* ic_matrix;
    FILE* move_log;
    FILE* move_matrix;
    FILE* icx;
    FILE* MI_matrix;
};

/////////////////////////////////////////////////////////////////////////////////////
/// TrxInfo
/////////////////////////////////////////////////////////////////////////////////////


class TrxInfo {
    const char* trj_file;
    int dimpca;
    const char* pca_file;

public:
    TrxInfo() : trj_file(nullptr), dimpca(0), pca_file(nullptr){
    }

    void setTrj(const char* inTrj){
        trj_file = inTrj;
    }

    void setPca(const char* inPca){
        pca_file = inPca;
    }
    void setDimPca(const int inDimPca){
        dimpca = inDimPca;
    }

    template <class FcaMasterClass>
    void write_fcalog(const char* fname, const FcaMasterClass& fca,
                      const char* posfile, const char* eigvecfile) {
        char str[STRLEN];
        FILE* out;
        time_t now;
        out = gmx_ffopen(fname, "w");
        now = time(nullptr);
        fprintf(out, "FCA - Analysis log, written %s\n", ctime(&now));
        char* retptr = getcwd(str, STRLEN);
        assert(retptr != nullptr);
        fprintf(out, "Working directory: %s\n\n", str);

        if(trj_file) {
            fprintf(out,
                    "Read %d principal components for intial projection from %s\n",
                    dimpca, pca_file);
        } else {
            fprintf(out, "Read data from %s\n", posfile);
        }

        fprintf(out, "\n");
        fprintf(out, "did FCA-Minimization in %d dimensions \n", fca.getDim());
        fprintf(out, "finally reached MI-sum: %10.5f \n", fca.getMIsum());
        if(eigvecfile){
            fprintf(out, "Wrote fca-vectors to %s\n", eigvecfile);
        }
        fclose(out);
    }
};

}
}

#endif
