/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */

#ifndef TNG_COMPRESS_H
#define TNG_COMPRESS_H

#ifndef USE_WINDOWS
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define USE_WINDOWS
#endif /* win32... */
#endif /* not defined USE_WINDOWS */

#ifndef DECLSPECDLLEXPORT
#ifdef USE_WINDOWS
#define DECLSPECDLLEXPORT __declspec(dllexport)
#else /* USE_WINDOWS */
#define DECLSPECDLLEXPORT
#endif /* USE_WINDOWS */
#endif /* DECLSPECDLLEXPORT */

#ifdef __cplusplus
 extern "C" {
#endif

/* tng_compress_pos expects positions to have the order:
   first xyz, then sorted in atom order
   then all the frames repeated, i.e.:
   nframes * [
    natoms* [
      x, y, z
    ]
   ]
   desired_precision what to round the numbers to, i.e. integers will be created as:
   round(pos[]/desired_precision).

   algo should first be determined by calling
   tng_compress_pos_find_algo

   The compressed data is returned in a malloced pointer (so free can
   be called to free the memory), the number of chars in the compressed
   data is put into *nitems.

   If too large values are input (compared to the precision), NULL is returned.
*/

char DECLSPECDLLEXPORT *tng_compress_pos(double *pos, const int natoms, const int nframes,
					 const double desired_precision,
					 const int speed, int *algo,
					 int *nitems);

char DECLSPECDLLEXPORT *tng_compress_pos_float(float *pos, const int natoms, const int nframes,
					       const float desired_precision,
					       const int speed, int *algo,
					       int *nitems);

char DECLSPECDLLEXPORT *tng_compress_pos_int(int *pos, const int natoms, const int nframes,
					     const unsigned long prec_hi, const unsigned long prec_lo,
					     int speed,int *algo,
					     int *nitems);

/* The tng_compress_pos_find_algo works the same as tng_compress_pos, but
   it performs benchmarking to find the algorithms with the best
   compression ratio.
   The search is controlled by giving speed:
   speed=1:  Fast algorithms only. This excludes all BWLZH algorithms and
             the XTC3 algorithm.
   speed=2:  Same as 1 and also includes the XTC3 algorithm using base compression
             only.
   speed=3:  Same as 2 and also includes the XTC3 algorithm which will use BWLZH
             compression when it seems likely to give better
             compression. Also includes the interframe BWLZH algorithm for
             coordinates and velocities.
   speed=4:  Enable the inter frame BWLZH algorithm for the coordinates.
             The one-to-one BWLZH algorithm is enabled for velocities.
   speed=5:  Enable the LZ77 part of the BWLZH algorithm.
   speed=6:  Enable the intra frame BWLZH algorithm for the coordinates. Always try
             the BWLZH compression in the XTC3 algorithm.

   Set speed=0 to allow tng_compression to set the default speed (which is currently 2).
   For very good compression it makes sense to choose speed=4 or speed=5

   The number of items required in the algorithm array can be found
   by calling tng_compress_nalgo
*/

char DECLSPECDLLEXPORT *tng_compress_pos_find_algo(double *pos, const int natoms, const int nframes,
						   const double desired_precision,
						   const int speed,
						   int *algo,
						   int *nitems);

char DECLSPECDLLEXPORT *tng_compress_pos_float_find_algo(float *pos, const int natoms, const int nframes,
							 const float desired_precision,
							 const int speed,
							 int *algo,
							 int *nitems);

char DECLSPECDLLEXPORT *tng_compress_pos_int_find_algo(int *pos, const int natoms, const int nframes,
						       const unsigned long prec_hi, const unsigned long prec_lo,
						       const int speed, int *algo,
						       int *nitems);

/* This returns the number of integers required for the storage of the algorithm
   with the best compression ratio. */
int DECLSPECDLLEXPORT tng_compress_nalgo(void);

/* The following two routines does the same as the compression of the
   positions, but compresses velocities instead. The algorithm
   selection for velocities is different, so the position and
   velocities routines should not be mixed. */

char DECLSPECDLLEXPORT *tng_compress_vel(double *vel, const int natoms, const int nframes,
					 const double desired_precision,
					 const int speed, int *algo,
					 int *nitems);

char DECLSPECDLLEXPORT *tng_compress_vel_float(float *vel, const int natoms, const int nframes,
					       const float desired_precision,
					       const int speed, int *algo,
					       int *nitems);

char DECLSPECDLLEXPORT *tng_compress_vel_int(int *vel, const int natoms, const int nframes,
					     const unsigned long prec_hi, const unsigned long prec_lo,
					     int speed, int *algo,
					     int *nitems);

char DECLSPECDLLEXPORT *tng_compress_vel_find_algo(double *vel, const int natoms, const int nframes,
						   const double desired_precision,
						   const int speed,
						   int *algo,
						   int *nitems);

char DECLSPECDLLEXPORT *tng_compress_vel_float_find_algo(float *vel, const int natoms, const int nframes,
							 const float desired_precision,
							 const int speed,
							 int *algo,
							 int *nitems);

char DECLSPECDLLEXPORT *tng_compress_vel_int_find_algo(int *vel, const int natoms, const int nframes,
						       const unsigned long prec_hi, const unsigned long prec_lo,
						       const int speed,
						       int *algo,
						       int *nitems);

/* From a compressed block, obtain information about
   whether it is a position or velocity block:
   *vel=1 means velocity block, *vel=0 means position block.
   It also gives info about the number of atoms,
   frames, and the precision used to compress the block, and the algorithms used to
   compress the block. The return value=0 if the block looks like a tng compressed block,
   and 1 otherwise. If the return value is 1 no information is returned. */
int DECLSPECDLLEXPORT tng_compress_inquire(char *data,int *vel, int *natoms,
					   int *nframes, double *precision,
					   int *algo);

/* Uncompresses any tng compress block, positions or velocities. It determines whether it is positions or velocities from the data buffer. The return value is 0 if ok, and 1 if not.
*/
int DECLSPECDLLEXPORT tng_compress_uncompress(char *data,double *posvel);

int DECLSPECDLLEXPORT tng_compress_uncompress_float(char *data,float *posvel);

int DECLSPECDLLEXPORT tng_compress_uncompress_int(char *data,int *posvel, unsigned long *prec_hi, unsigned long *prec_lo);

/* This converts a block of integers, as obtained from tng_compress_uncompress_int, to floating point values
   either double precision or single precision. */
void DECLSPECDLLEXPORT tng_compress_int_to_double(int *posvel_int, const unsigned long prec_hi, const unsigned long prec_lo,
						  const int natoms, const int nframes,
						  double *posvel_double);

void DECLSPECDLLEXPORT tng_compress_int_to_float(int *posvel_int, const unsigned long prec_hi, const unsigned long prec_lo,
						 const int natoms, const int nframes,
						 float *posvel_float);


/* Compression algorithms (matching the original trajng
   assignments) The compression backends require that some of the
   algorithms must have the same value. */

#define TNG_COMPRESS_ALGO_STOPBIT 1
#define TNG_COMPRESS_ALGO_TRIPLET 2
#define TNG_COMPRESS_ALGO_BWLZH1  8
#define TNG_COMPRESS_ALGO_BWLZH2  9

#define TNG_COMPRESS_ALGO_POS_STOPBIT_INTER     TNG_COMPRESS_ALGO_STOPBIT
#define TNG_COMPRESS_ALGO_POS_TRIPLET_INTER     TNG_COMPRESS_ALGO_TRIPLET
#define TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA     3
#define TNG_COMPRESS_ALGO_POS_XTC2              5
#define TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE  7
#define TNG_COMPRESS_ALGO_POS_BWLZH_INTER       TNG_COMPRESS_ALGO_BWLZH1
#define TNG_COMPRESS_ALGO_POS_BWLZH_INTRA       TNG_COMPRESS_ALGO_BWLZH2
#define TNG_COMPRESS_ALGO_POS_XTC3              10
#define TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE  TNG_COMPRESS_ALGO_STOPBIT
#define TNG_COMPRESS_ALGO_VEL_TRIPLET_INTER     TNG_COMPRESS_ALGO_TRIPLET
#define TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE  3
#define TNG_COMPRESS_ALGO_VEL_STOPBIT_INTER     6
#define TNG_COMPRESS_ALGO_VEL_BWLZH_INTER       TNG_COMPRESS_ALGO_BWLZH1
#define TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE    TNG_COMPRESS_ALGO_BWLZH2
#define TNG_COMPRESS_ALGO_MAX 11



/* Obtain strings describing the actual algorithms. These point to static memory, so should
   not be freed. */
char DECLSPECDLLEXPORT *tng_compress_initial_pos_algo(int *algo);
char DECLSPECDLLEXPORT *tng_compress_pos_algo(int *algo);
char DECLSPECDLLEXPORT *tng_compress_initial_vel_algo(int *algo);
char DECLSPECDLLEXPORT *tng_compress_vel_algo(int *algo);



#ifdef __cplusplus
 }
#endif


#endif
