/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*
 * TOMS Transpose.  Revised version of algorithm 380.
 * 
 * These routines do in-place transposes of arrays.
 * 
 * [ Cate, E.G. and Twigg, D.W., ACM Transactions on Mathematical Software, 
 *   vol. 3, no. 1, 104-110 (1977) ]
 * 
 * C version by Steven G. Johnson. February 1997.
 */

#include <stdlib.h>
#include <string.h>

#include "TOMS_transpose.h"

static int TOMS_gcd(int a, int b);

/*
 * "a" is a 1D array of length ny*nx which constains the nx x ny matrix to be
 * transposed.  "a" is stored in C order (last index varies fastest).  move
 * is a 1D array of length move_size used to store information to speed up
 * the process.  The value move_size=(ny+nx)/2 is recommended.
 * 
 * The return value indicates the success or failure of the routine. Returns 0
 * if okay, -1 if ny or nx < 0, and -2 if move_size < 1. The return value
 * should never be positive, but it it is, it is set to the final position in
 * a when the search is completed but some elements have not been moved.
 * 
 * Note: move[i] will stay zero for fixed points.
 */

short TOMS_transpose_2d(TOMS_el_type * a,
			int nx, int ny,
			char *move,
			int move_size)
{
	int             i, j, im, mn;
	TOMS_el_type    b, c, d;
	int             ncount;
	int             k;

	/* check arguments and initialize: */
	if (ny < 0 || nx < 0)
		return -1;
	if (ny < 2 || nx < 2)
		return 0;
	if (move_size < 1)
		return -2;

	if (ny == nx) {
		/*
		 * if matrix is square, exchange elements a(i,j) and a(j,i):
		 */
		for (i = 0; i < nx; ++i)
			for (j = i + 1; j < nx; ++j) {
				b = a[i + j * nx];
				a[i + j * nx] = a[j + i * nx];
				a[j + i * nx] = b;
			}
		return 0;
	}
	ncount = 2;		/* always at least 2 fixed points */
	k = (mn = ny * nx) - 1;

	for (i = 0; i < move_size; ++i)
		move[i] = 0;

	if (ny >= 3 && nx >= 3)
		ncount += TOMS_gcd(ny - 1, nx - 1) - 1;	/* # fixed points */

	i = 1;
	im = ny;

	while (1) {
		int             i1, i2, i1c, i2c;
		int             kmi;

		/** Rearrange the elements of a loop
	              and its companion loop: **/

		i1 = i;
		kmi = k - i;
		b = a[i1];
		i1c = kmi;
		c = a[i1c];

		while (1) {
			i2 = ny * i1 - k * (i1 / nx);
			i2c = k - i2;
			if (i1 < move_size)
				move[i1] = 1;
			if (i1c < move_size)
				move[i1c] = 1;
			ncount += 2;
			if (i2 == i)
				break;
			if (i2 == kmi) {
				d = b;
				b = c;
				c = d;
				break;
			}
			a[i1] = a[i2];
			a[i1c] = a[i2c];
			i1 = i2;
			i1c = i2c;
		}
		a[i1] = b;
		a[i1c] = c;

		if (ncount >= mn)
			break;	/* we've moved all elements */

		/** Search for loops to rearrange: **/

		while (1) {
			int             max;

			max = k - i;
			++i;
			if (i > max)
				return i;
			im += ny;
			if (im > k)
				im -= k;
			i2 = im;
			if (i == i2)
				continue;
			if (i >= move_size) {
				while (i2 > i && i2 < max) {
					i1 = i2;
					i2 = ny * i1 - k * (i1 / nx);
				}
				if (i2 == i)
					break;
			} else if (!move[i])
				break;
		}
	}

	return 0;
}

/*
 * "a" is a 1D array of length ny*nx which constains the nx x ny matrix to be
 * transposed.  "a" is stored in C order (last index varies fastest).  move
 * is a 1D array of length move_size used to store information to speed up
 * the process.  The value move_size=(ny+nx)/2 is recommended.
 * 
 * Here, instead of each element of "a" being a single value of type
 * TOMS_el_type, each element is el_size values of type TOMS_el_type.
 * 
 * The return value indicates the success or failure of the routine. Returns 0
 * if okay, -1 if ny or nx < 0, and -2 if move_size < 1. Also, returns -3 if
 * it ran out of memory.  The return value should never be positive, but it
 * it is, it is set to the final position in a when the search is completed
 * but some elements have not been moved.
 * 
 * Note: move[i] will stay zero for fixed points.
 */
short TOMS_transpose_2d_arbitrary(TOMS_el_type * a,
				  int nx, int ny,
				  int el_size,
				  char *move,
				  int move_size)
{
	int             i, j, im, mn;
	TOMS_el_type   *b, *c, *d;
	int             ncount;
	int             k;

	/* check arguments and initialize: */
	if (ny < 0 || nx < 0)
		return -1;
	if (ny < 2 || nx < 2 || el_size < 1)
		return 0;
	if (move_size < 1)
		return -2;

	b = (TOMS_el_type *) malloc(sizeof(TOMS_el_type) * el_size);
	if (!b)
		return -3;

	if (ny == nx) {
		/*
		 * if matrix is square, exchange elements a(i,j) and a(j,i):
		 */
		for (i = 0; i < nx; ++i)
			for (j = i + 1; j < nx; ++j) {
				memcpy(b, &a[el_size * (i + j * nx)], el_size * sizeof(TOMS_el_type));
				memcpy(&a[el_size * (i + j * nx)], &a[el_size * (j + i * nx)], el_size * sizeof(TOMS_el_type));
				memcpy(&a[el_size * (j + i * nx)], b, el_size * sizeof(TOMS_el_type));
			}
		free(b);
		return 0;
	}
	c = (TOMS_el_type *) malloc(sizeof(TOMS_el_type) * el_size);
	if (!c) {
		free(b);
		return -3;
	}
	ncount = 2;		/* always at least 2 fixed points */
	k = (mn = ny * nx) - 1;

	for (i = 0; i < move_size; ++i)
		move[i] = 0;

	if (ny >= 3 && nx >= 3)
		ncount += TOMS_gcd(ny - 1, nx - 1) - 1;	/* # fixed points */

	i = 1;
	im = ny;

	while (1) {
		int             i1, i2, i1c, i2c;
		int             kmi;

		/** Rearrange the elements of a loop
	              and its companion loop: **/

		i1 = i;
		kmi = k - i;
		memcpy(b, &a[el_size * i1], el_size * sizeof(TOMS_el_type));
		i1c = kmi;
		memcpy(c, &a[el_size * i1c], el_size * sizeof(TOMS_el_type));

		while (1) {
			i2 = ny * i1 - k * (i1 / nx);
			i2c = k - i2;
			if (i1 < move_size)
				move[i1] = 1;
			if (i1c < move_size)
				move[i1c] = 1;
			ncount += 2;
			if (i2 == i)
				break;
			if (i2 == kmi) {
				d = b;
				b = c;
				c = d;
				break;
			}
			memcpy(&a[el_size * i1], &a[el_size * i2], 
			       el_size * sizeof(TOMS_el_type));
			memcpy(&a[el_size * i1c], &a[el_size * i2c], 
			       el_size * sizeof(TOMS_el_type));
			i1 = i2;
			i1c = i2c;
		}
		memcpy(&a[el_size * i1], b, el_size * sizeof(TOMS_el_type));
		memcpy(&a[el_size * i1c], c, el_size * sizeof(TOMS_el_type));

		if (ncount >= mn)
			break;	/* we've moved all elements */

		/** Search for loops to rearrange: **/

		while (1) {
			int             max;

			max = k - i;
			++i;
			if (i > max) {
				free(b);
				free(c);
				return i;
			}
			im += ny;
			if (im > k)
				im -= k;
			i2 = im;
			if (i == i2)
				continue;
			if (i >= move_size) {
				while (i2 > i && i2 < max) {
					i1 = i2;
					i2 = ny * i1 - k * (i1 / nx);
				}
				if (i2 == i)
					break;
			} else if (!move[i])
				break;
		}
	}

	free(b);
	free(c);
	return 0;
}

static int TOMS_gcd(int a, int b)
{
	int r;
	do {
		r = a % b;
		a = b;
		b = r;
	} while (r != 0);

	return a;
}
