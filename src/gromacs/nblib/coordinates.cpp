#include "coordinates.h"

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"


void generateCoordinates(int                     multiplicationFactor,
                        std::vector<gmx::RVec> *coordinates,
                        matrix                  box)
{
    if (multiplicationFactor < 1 ||
        (multiplicationFactor & (multiplicationFactor - 1)) != 0)
    {
        gmx_fatal(FARGS, "The size factor has to be a power of 2");
    }

    if (multiplicationFactor == 1)
    {
        *coordinates = coordinates12;
        copy_mat(box12, box);

        return;
    }

    ivec factors = { 1, 1, 1 };

    int  dim = 0;
    while (multiplicationFactor > 1)
    {
        factors[dim]         *= 2;
        multiplicationFactor /= 2;
        dim++;
        if (dim == DIM)
        {
            dim = 0;
        }
    }
    printf("Stacking a box of %zu atoms %d x %d x %d times\n",
           coordinates12.size(), factors[XX], factors[YY], factors[ZZ]);

    coordinates->resize(factors[XX]*factors[YY]*factors[ZZ]*coordinates12.size());

    int       i = 0;
    gmx::RVec shift;
    for (int x = 0; x < factors[XX]; x++)
    {
        shift[XX] = x*box12[XX][XX];
        for (int y = 0; y < factors[YY]; y++)
        {
            shift[YY] = y*box12[YY][YY];
            for (int z = 0; z < factors[ZZ]; z++)
            {
                shift[ZZ] = z*box12[ZZ][ZZ];

                for (const gmx::RVec &coordOrig : coordinates12)
                {
                    (*coordinates)[i] = coordOrig + shift;
                    i++;
                }
            }
        }
    }

    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            box[d1][d2] = factors[d1]*box12[d1][d2];
        }
    }
}
