/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#ifndef MERGE_SORT_H
#define MERGE_SORT_H

void Ptngc_merge_sort(void *base, const size_t nmemb, const size_t size,
		int (*compar)(const void *v1,const void *v2,const void *private),
		void *private);


#endif
