/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "binsearch.h"

#include "gromacs/utility/real.h"

/*Make range-array (Permutation identity) for sorting */
void rangeArray(int* ar, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        ar[i] = i;
    }
}

static void pswap(int* v1, int* v2)
{
    int temp;
    temp = *v1;
    *v1  = *v2;
    *v2  = temp;
}


static void Swap(real* v1, real* v2)
{
    real temp;
    temp = *v1;
    *v1  = *v2;
    *v2  = temp;
}


void insertionSort(real* arr, int* perm, int startndx, int endndx, int direction)
{
    int i, j;

    if (direction >= 0)
    {
        for (i = startndx; i <= endndx; i++)
        {
            j = i;

            while (j > startndx && arr[j - 1] > arr[j])
            {
                Swap(&arr[j], &arr[j - 1]);
                pswap(&perm[j], &perm[j - 1]);
                j--;
            }
        }
    }

    if (direction < 0)
    {
        for (i = startndx; i <= endndx; i++)
        {
            j = i;

            while (j > startndx && arr[j - 1] < arr[j])
            {
                Swap(&arr[j], &arr[j - 1]);
                pswap(&perm[j], &perm[j - 1]);
                j--;
            }
        }
    }
}


int BinarySearch(const real* array, int low, int high, real key, int direction)
{
    int iMid, iMax, iMin;
    iMax = high + 2;
    iMin = low + 1;

    /*Iterative implementation*/

    if (direction >= 0)
    {
        while (iMax - iMin > 1)
        {
            iMid = (iMin + iMax) >> 1;
            if (key < array[iMid - 1])
            {
                iMax = iMid;
            }
            else
            {
                iMin = iMid;
            }
        }
        return iMin;
    }
    else
    {
        while (iMax - iMin > 1)
        {
            iMid = (iMin + iMax) >> 1;
            if (key > array[iMid - 1])
            {
                iMax = iMid;
            }
            else
            {
                iMin = iMid;
            }
        }
        return iMin - 1;
    }
}


int start_binsearch(real* array, int* perm, int low, int high, real key, int direction)
{
    insertionSort(array, perm, low, high, direction);
    return BinarySearch(array, low, high, key, direction);
}
