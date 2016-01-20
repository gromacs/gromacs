/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMX_RESP_ATOM_H
#define GMX_RESP_ATOM_H

#include <vector>

#include "gromacs/utility/real.h"

#include "poldata.h"

namespace alexandria
{
class RespAtom
{
    public:
        RespAtom(int atomnumber, 
                 int atype,
                 const char *atomtype, 
                 const Poldata &pd,
                 ChargeDistributionModel iDistributionModel, 
                 std::vector<std::string>  dzatoms);
        ~RespAtom();

        real getQ();

        bool setUpcorrectly();

        real getZeta(int i)
        {
            return _zeta[i];
        }

        void setZeta(int i, real value)
        {
            _zeta[i] = value;
        }

        int getNZeta()
        {
            return _nZeta;
        }

        int getIq(int i)
        {
            return _iq[i];
        }

        void setIq(int i, int value)
        {
            _iq[i] = value;
        }

        int getIz(int i)
        {
            return _iz[i];
        }

        void setIz(int i, int value)
        {
            _iz[i] = value;
        }


        real getQ(int i)
        {
            return _q[i];
        }

        void setQ(int i, real value)
        {
            _q[i] = value;
        }

        int getRow(int i)
        {
            return _row[i];
        }


        int getAtomnumber()
        {
            return _atomnumber;
        }

        void setAtomnumber(int value)
        {
            _atomnumber = value;
        }

        real getZetaRef(int i)
        {
            return _zetaRef[i];
        }

        std::string getAtomtype()
        {
            return _atomtype;
        }

        bool getBRestrained()
        {
            return _bRestrained;
        }

        int getAtype()
        {
            return _atype;
        }

        void setAtype(int i )
        {
            _atype = i;
        }

    private:

        std::vector<int>  _row;
        int               _atomnumber, _atype;
        bool              _bRestrained;
        std::string       _atomtype;
        int               _nZeta;
        std::vector<real> _q, _zeta, _zetaRef;
        std::vector<int>  _iq, _iz;
        bool              _bSetUpcorrectly;

};

} // namespace

#endif
