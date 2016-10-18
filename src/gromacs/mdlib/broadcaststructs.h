#ifndef GMX_MDLIB_BROADCASTSTRUCTS_H
#define GMX_MDLIB_BROADCASTSTRUCTS_H

#include <cmath>
#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/smalloc.h"

template <typename T>
void block_bc(const t_commrec *cr, T &data)
{
    gmx_bcast(sizeof(T), static_cast<void *>(&data), cr);
}
template <typename T>
void nblock_bc(const t_commrec *cr, int numElements, T *data)
{
    gmx_bcast(numElements * sizeof(T), static_cast<void *>(data), cr);
}
template <typename T>
void snew_bc(const t_commrec *cr, T * &data, int numElements)
{
    if (!MASTER(cr))
    {
        snew(data, numElements);
    }
}
template <typename T>
static void nblock_abc(const t_commrec *cr, int numElements, T **v)
{
    snew_bc(cr, v, numElements);
    nblock_bc(cr, numElements, *v);
}

template <typename T>
static void nblock_abc(const t_commrec *cr, int numElements, std::vector<T> *v)
{
    if (!MASTER(cr))
    {
        v->resize(numElements);
    }
    gmx_bcast(numElements*sizeof(T), v->data(), cr);
}

#endif
