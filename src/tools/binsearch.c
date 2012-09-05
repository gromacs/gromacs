/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: densorder.c,v 0.9
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "types/simple.h"
#include "gmx_fatal.h"

/*Make range-array (Permutation identity) for sorting */
void rangeArray(int *ar,int size)
{
    int i;	
	for (i=0;i<size;i++){
		ar[i]=i;
	}
}

void pswap(int *v1, int *v2)
{
	int temp;
	temp=*v1;
	*v1=*v2;
	*v2=temp;
}


void Swap (real *v1, real *v2)
{
    real temp;
    temp = *v1;
    *v1 = *v2;
    *v2 = temp;
}



void insertionSort(real *arr, int *perm, int startndx, int endndx, int direction) {
    int i, j;

	if(direction>=0){
        for (i = startndx; i <=endndx; i++) {
            j = i;

            while (j > startndx && arr[j - 1] > arr[j]) {
                Swap(&arr[j],&arr[j-1]);
				pswap(&perm[j],&perm[j-1]);
                j--;
            }

        }

	}

	if(direction<0){
        for (i = startndx; i <=endndx; i++) {
            j = i;

            while (j > startndx && arr[j - 1] < arr[j]) {
                Swap(&arr[j],&arr[j-1]);
				pswap(&perm[j],&perm[j-1]);
                j--;
            }

        }
	}
}


int BinarySearch (real *array, int low, int high, real key,int direction)
{
    int mid, max, min;
    max=high+2;
    min=low+1;

/*Iterative implementation*/

    if (direction>=0){
        while (max-min>1){
            mid=(min+max)>>1;
            if(key<array[mid-1]) max=mid;
            else min=mid;
        }
        return min;
    }

    else if (direction<0){
        while(max-min>1){
            mid=(min+max)>>1;
            if(key>array[mid-1]) max=mid;
            else min=mid;
        }
        return min-1;

    }/*end -ifelse direction*/
   return -1;
}


int start_binsearch(real *array, int *perm, int low, int high, 
                    real key, int direction)
{
	insertionSort(array,perm,low,high,direction);
	return BinarySearch(array,low,high,key,direction);
}

int LinearSearch (double *array,int startindx, int stopindx, 
                  double key,int *count, int direction)
{
    /*Iterative implementation - assume elements sorted*/
    int i;
    int keyindex;

    if(direction>=0){
        keyindex=startindx;
        for (i=startindx;i<=stopindx;i++){
            (*count)++;
            if(array[i]>key) {	
                keyindex=i-1;
                return keyindex;
            }
        }
    }
    else if(direction<0 ){
        keyindex=stopindx;
        for(i=stopindx;i>=startindx;i--){
            (*count)++;
            if (array[i]>key){
                keyindex=i+1;
                return keyindex;
            }
        }
    }

    else 
        gmx_fatal(FARGS,"Startindex=stopindex=%d\n",startindx);
        
    return -1;
}
		
