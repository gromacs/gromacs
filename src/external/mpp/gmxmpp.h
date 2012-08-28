#include "gromacs/utility/gmx_header_config.h"
// For GMX_LIB_MPI
#ifdef GMX_CXX11
#define MPP_CXX11_RVALREF
#endif
#ifdef GMX_LIB_MPI
#include "mpp.h"
#else
#ifdef GMX_THREAD_MPI
#include <tmpi.h>
#else
typedef void* MPI_Datatype;
extern const MPI_Datatype MPI_DOUBLE;
extern const MPI_Datatype MPI_INT;
extern const MPI_Datatype MPI_CHAR;
extern const MPI_Datatype MPI_FLOAT;
extern const MPI_Datatype MPI_LONG;
extern const MPI_Datatype MPI_BYTE;
int MPI_Type_commit(MPI_Datatype *datatype);
#endif //GMX_THREAD_MPI
#include <cstddef>
typedef std::size_t MPI_Aint;
extern const MPI_Datatype MPI_UNSIGNED_LONG;
int MPI_Get_address(void*,MPI_Aint*);
int MPI_Type_indexed(int, int*, int*, MPI_Datatype, MPI_Datatype*);
#define MPP_NO_MPI_INCL
#include "type_traits.h"
#endif //GMX_LIB_MPI
