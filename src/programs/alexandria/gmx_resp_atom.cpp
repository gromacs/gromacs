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

RespAtom::RespAtom(int atomnumber, 
                   int atype,
                   const char *atomtype, 
                   const Poldata &pd,
                   ChargeDistributionModel iDistributionModel,
                   std::vector<std::string> dzatoms)
{
    int         k, zz = 0, shell;
    bool        bRestr;
    std::string atomtype_new;
    size_t      shell_name;

    bRestr = false;
    if (!dzatoms.empty())
    {
        k = 0;
        while (("" != dzatoms[k]) && !bRestr)
        {
            bRestr = (strcasecmp(atomtype, dzatoms[k].c_str()) == 0);
            k++;
        }
    }
    // _nZeta       = pd->getNzeta(iDistributionModel, _atomtype.c_str());
    _atomnumber  = atomnumber;
    _atype       = atype;
    _atomtype    = (atomtype);
    atomtype_new = _atomtype;
    _nZeta       = 1;
    shell        = 0;

    /*if (_nZeta <= 0)
       {
        _bSetUpcorrectly = false;
       return;
       }*/

    _bRestrained = bRestr;

    _zeta.resize(_nZeta);
    _zetaRef.resize(_nZeta);
    _q.resize(_nZeta);
    _iz.resize(_nZeta);
    _iq.resize(_nZeta);
    _row.resize(_nZeta);

    //for (zz = 0; (zz < _nZeta); zz++)
    //{	}

    _iq[zz]        = -1;
    _iz[zz]        = -1;

    shell_name = _atomtype.find("_s");
    if (shell_name != std::string::npos)
    {
        shell        = 1;
        atomtype_new = _atomtype.substr(0, shell_name);
    }

    _q[zz]            = pd.getQ( iDistributionModel, atomtype_new.c_str(), shell);
    _zetaRef[zz]      =
        _zeta[zz]     = pd.getZeta( iDistributionModel, atomtype_new.c_str(), shell);
    _row[zz]          = pd.getRow( iDistributionModel, atomtype_new.c_str(), shell);

    _bSetUpcorrectly = true;
}


bool RespAtom::setUpcorrectly()
{
    return _bSetUpcorrectly;
}

RespAtom::~RespAtom()
{
    /* sfree(_iz);
       sfree(_iq);
       sfree(_zeta);
       sfree(_q);
       sfree(_row);*/
    // sfree(_atomtype);
}

real RespAtom::getQ()
{
    int  i;
    real q = 0;

    for (i = 0; (i < _nZeta); i++)
    {
        q += _q[i];
    }

    return q;
}
}
