#ifndef GMX_RAMD_READPARAMS_H
#define GMX_RAMD_READPARAMS_H

#include "gromacs/fileio/readinp.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

#include "ramdparameters.h"

struct t_grpopts;
struct t_inputrec;
struct gmx_mtop_t;
struct pull_params_t;
struct pull_t;
enum class PbcType : int;
class WarningHandler;

namespace gmx
{

class ISerializer;

void read_params(std::vector<t_inpfile>* inp, RAMDParameters* ramdparams, WarningHandler* wi);

} // namespace gmx

#endif // GMX_RAMD_READPARAMS_H
