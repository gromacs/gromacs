/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmx_resp_atom.h"

#include "gromacs/topology/atoms.h"

#include "gentop_qgen.h"
#include "gmx_resp.h"
#include "poldata.h"
#include "stringutil.h"

namespace alexandria
{

RespAtomType::RespAtomType(int atype,
                           const char *atomtype, 
                           const Poldata &pd,
                           ChargeDistributionModel iDistributionModel,
                           const std::vector<std::string> &dzatoms)
{
    int         k, shell;
    std::string atomtype_new;

    bool bRestr = false;
    if (!dzatoms.empty())
    {
        k = 0;
        while (dzatoms[k].size() > 0 && !bRestr)
        {
            bRestr = (strcasecmp(atomtype, dzatoms[k].c_str()) == 0);
            k++;
        }
    }
    atype_       = atype;
    atomtype_    = atomtype;
    atomtype_new = atomtype_;
    shell        = 0;

    bRestrained_ = bRestr;

    int nZeta    = std::max(1, pd.getNzeta(iDistributionModel, atomtype_));
    for(int i = 0; i < nZeta-1; i++)
    {
        rz_.push_back(RespZeta(0, 0, 0));
    }

    size_t shell_name = atomtype_.find("_s");
    if (shell_name != std::string::npos)
    {
        shell        = 1;
        atomtype_new = atomtype_.substr(0, shell_name);
    }
    rz_.push_back(RespZeta(pd.getRow(iDistributionModel, atomtype_new, shell),
                           pd.getQ(iDistributionModel, atomtype_new, shell),
                           pd.getZeta(iDistributionModel, atomtype_new, shell)));
}

} // namespace
