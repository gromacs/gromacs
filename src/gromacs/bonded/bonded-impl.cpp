#include "bonded-impl.h"

#include <cstdlib>

int glatnr(int *global_atom_index, int i)
{
    int atnr;

    if (global_atom_index == NULL)
    {
        atnr = i + 1;
    }
    else
    {
        atnr = global_atom_index[i] + 1;
    }

    return atnr;
}

