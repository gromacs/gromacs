#include <stdio.h>
#include "physics.h"

int main(int argc, char *argv[])
{
    int    i;
    double x, y, z;

    x = 3.25;
    for (i = 0; (i < eg2cNR); i++)
    {
        y = gmx2convert(x, i);
        z = convert2gmx(y, i);
        printf("Converted %g [type %d] to %g and back to %g. Diff %g\n",
               x, i, y, z, x-z);
    }
}
