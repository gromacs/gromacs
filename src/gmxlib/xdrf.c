/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */



#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <stdio.h>
#include <stdlib.h>
#include "xdrf.h"

int ftocstr(char *, int, char *, int);
int ctofstr(char *, int, char *);

#define MAXID 20
static FILE *xdrfiles[MAXID];
static XDR *xdridptr[MAXID];
static char xdrmodes[MAXID];
static unsigned int cnt;

typedef void (* xdrfproc_) (int *, void *, int *);

void
xdrfbool_ (xdrid, pb, ret)
int *xdrid, *ret;
int *pb;
{
	xdr_bool(xdridptr[*xdrid], pb);
	cnt += sizeof(int);
}

void
xdrfchar_ (xdrid, cp, ret)
int *xdrid, *ret;
char *cp;
{
	*ret = xdr_char(xdridptr[*xdrid], cp);
	cnt += sizeof(char);
}

void
xdrfdouble_ (xdrid, dp, ret)
int *xdrid, *ret;
double *dp;
{
	*ret = xdr_double(xdridptr[*xdrid], dp);
	cnt += sizeof(double);
}

void
xdrffloat_ (xdrid, fp, ret)
int *xdrid, *ret;
float *fp;
{
	*ret = xdr_float(xdridptr[*xdrid], fp);
	cnt += sizeof(float);
}

void
xdrfint_ (xdrid, ip, ret)
int *xdrid, *ret;
int *ip;
{
	*ret = xdr_int(xdridptr[*xdrid], ip);
	cnt += sizeof(int);
}

void
xdrflong_ (xdrid, lp, ret)
int *xdrid, *ret;
long *lp;
{
	*ret = xdr_long(xdridptr[*xdrid], lp);
	cnt += sizeof(long);
}

void
xdrfshort_ (xdrid, sp, ret)
int *xdrid, *ret;
short *sp;
{
	*ret = xdr_short(xdridptr[*xdrid], sp);
	cnt += sizeof(sp);
}

void
xdrfuchar_ (xdrid, ucp, ret)
int *xdrid, *ret;
char *ucp;
{
	*ret = xdr_u_char(xdridptr[*xdrid], (u_char *) ucp);
	cnt += sizeof(char);
}

void
xdrfulong_ (xdrid, ulp, ret)
int *xdrid, *ret;
unsigned long *ulp;
{
	*ret = xdr_u_long(xdridptr[*xdrid], ulp);
	cnt += sizeof(unsigned long);
}

void
xdrfushort_ (xdrid, usp, ret)
int *xdrid, *ret;
unsigned short *usp;
{
	*ret = xdr_u_short(xdridptr[*xdrid], usp);
	cnt += sizeof(unsigned short);
}

void 
xdrf3dfcoord_ (xdrid, fp, size, precision, ret)
int *xdrid, *ret;
float *fp;
int *size;
float *precision;
{
	*ret = xdr3dfcoord(xdridptr[*xdrid], fp, size, precision);
}

void
xdrfstring_ (xdrid, sp_ptr, maxsize, ret, sp_len)
int *xdrid, *ret;
char * sp_ptr; int sp_len;
int *maxsize;
{
	char *tsp;

	tsp = (char*) malloc(((sp_len) + 1) * sizeof(char));
	if (tsp == NULL) {
	    *ret = -1;
	    return;
	}
	if (ftocstr(tsp, *maxsize+1, sp_ptr, sp_len)) {
	    *ret = -1;
	    free(tsp);
	    return;
	}
	*ret = xdr_string(xdridptr[*xdrid], &tsp, *maxsize);
	ctofstr( sp_ptr, sp_len, tsp);
	cnt += *maxsize;
	free(tsp);
}

void
xdrfwrapstring_ (xdrid,  sp_ptr, ret, sp_len)
int *xdrid, *ret;
char * sp_ptr; int sp_len;
{
	char *tsp;
	int maxsize;
	maxsize = (sp_len) + 1;
	tsp = (char*) malloc(maxsize * sizeof(char));
	if (tsp == NULL) {
	    *ret = -1;
	    return;
	}
	if (ftocstr(tsp, maxsize, sp_ptr, sp_len)) {
	    *ret = -1;
	    free(tsp);
	    return;
	}
	*ret = xdr_string(xdridptr[*xdrid], &tsp, maxsize);
	ctofstr( sp_ptr, sp_len, tsp);
	cnt += maxsize;
	free(tsp);
}

void
xdrfopaque_ (xdrid, cp, ccnt, ret)
int *xdrid, *ret;
caddr_t *cp;
int *ccnt;
{
	*ret = xdr_opaque(xdridptr[*xdrid], *cp, *ccnt);
	cnt += *ccnt;
}

void
xdrfsetpos_ (xdrid, pos, ret)
int *xdrid, *ret;
int *pos;
{
	*ret = xdr_setpos(xdridptr[*xdrid], *pos);
}

void
xdrf_ (xdrid, pos)
int *xdrid, *pos;
{
	*pos = xdr_getpos(xdridptr[*xdrid]);
}

void
xdrfvector_ (xdrid, cp, size, elproc, ret)
int *xdrid, *ret;
char *cp;
int *size;
xdrfproc_ elproc;
{
	int lcnt;
	cnt = 0;
	for (lcnt = 0; lcnt < *size; lcnt++) {
		elproc(xdrid, (cp+cnt) , ret);
	}
}


void
xdrfclose_ (xdrid, ret)
int *xdrid;
int *ret;
{
	*ret = xdrclose(xdridptr[*xdrid]);
	cnt = 0;
}

void
xdrfopen_ (xdrid,  fp_ptr, mode_ptr, ret, fp_len, mode_len)
int *xdrid;
char * fp_ptr; int fp_len;
char * mode_ptr; int mode_len;
int *ret;
{
	char fname[512];
	char fmode[3];

	if (ftocstr(fname, sizeof(fname), fp_ptr, fp_len)) {
		*ret = 0;
	}
	if (ftocstr(fmode, sizeof(fmode), mode_ptr,
			mode_len)) {
		*ret = 0;
	}

	*xdrid = xdropen(NULL, fname, fmode);
	if (*xdrid == 0)
		*ret = 0;
	else 
		*ret = 1;	
}

/*___________________________________________________________________________
 |
 | what follows are the C routines for opening, closing xdr streams
 | and the routine to read/write compressed coordinates together
 | with some routines to assist in this task (those are marked
 | static and cannot be called from user programs)
*/
#define MAXABS INT_MAX-2

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x):(y))
#endif
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
static int magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
    8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
    82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))


/*__________________________________________________________________________
 |
 | xdropen - open xdr file
 |
 | This versions differs from xdrstdio_create, because I need to know
 | the state of the file (read or write) so I can use xdr3dfcoord
 | in eigther read or write mode, and the file descriptor
 | so I can close the file (something xdr_destroy doesn't do).
 |
*/

int xdropen(XDR *xdrs, const char *filename, const char *type) {
    static int init_done = 0;
    enum xdr_op lmode;
    int xdrid;
    
    if (init_done == 0) {
	for (xdrid = 1; xdrid < MAXID; xdrid++) {
	    xdridptr[xdrid] = NULL;
	}
	init_done = 1;
    }
    xdrid = 1;
    while (xdrid < MAXID && xdridptr[xdrid] != NULL) {
	xdrid++;
    }
    if (xdrid == MAXID) {
	fprintf(stderr,"xdropen: couldn't open file '%s'\n", filename);
	return 0;
    }
    if (*type == 'w' || *type == 'W') {
      type = "w+";
      lmode = XDR_ENCODE;
    } else if (*type == 'a' || *type == "A") {
      type = "a+";
      lmode = XDR_ENCODE;
    } else {
      type = "r";
      lmode = XDR_DECODE;
    }
    if ((xdrfiles[xdrid] = fopen(filename, type))==NULL) {
      fprintf(stderr,"Sorry, I couldn't open %s\n",filename);
      exit(0);
    }
    if (xdrfiles[xdrid] == NULL) {
	xdrs = NULL;
	return 0;
    }
    xdrmodes[xdrid] = *type;
    /* next test isn't usefull in the case of C language
     * but is used for the Fortran interface
     * (C users are expected to pass the address of an already allocated
     * XDR staructure)
     */
    if (xdrs == NULL) {
	xdridptr[xdrid] = (XDR *) malloc(sizeof(XDR));
	xdrstdio_create(xdridptr[xdrid], xdrfiles[xdrid], lmode);
    } else {
	xdridptr[xdrid] = xdrs;
	xdrstdio_create(xdrs, xdrfiles[xdrid], lmode);
    }
    return xdrid;
}

/*_________________________________________________________________________
 |
 | xdrclose - close a xdr file
 |
 | This will flush the xdr buffers, and destroy the xdr stream.
 | It also closes the associated file descriptor (this is *not*
 | done by xdr_destroy).
 |
*/
 
int xdrclose(XDR *xdrs) {
    int xdrid;
    
    if (xdrs == NULL) {
	fprintf(stderr, "xdrclose: passed a NULL pointer\n");
	return 0;
    }
    for (xdrid = 1; xdrid < MAXID; xdrid++) {
	if (xdridptr[xdrid] == xdrs) {
	    
	    xdr_destroy(xdrs);
	    fclose(xdrfiles[xdrid]);
	    xdridptr[xdrid] = NULL;
	    return 1;
	}
    } 
    fprintf(stderr, "xdrclose: no such open xdr file\n");
    return 0;
    
}

/*____________________________________________________________________________
 |
 | sendbits - encode num into buf using the specified number of bits
 |
 | This routines appends the value of num to the bits already present in
 | the array buf. You need to give it the number of bits to use and you
 | better make sure that thi snumber if bits is enough to hold the value
 | Also num must be positive.
 |
*/

static void sendbits(int buf[], int num_of_bits, int num) {
    
    unsigned int cnt, lastbyte;
    int lastbits;
    unsigned char * cbuf;
    
    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte =(unsigned int) buf[2];
    while (num_of_bits >= 8) {
	lastbyte = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
	cbuf[cnt++] = lastbyte >> lastbits;
	num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
	lastbyte = (lastbyte << num_of_bits) | num;
	lastbits += num_of_bits;
	if (lastbits >= 8) {
	    lastbits -= 8;
	    cbuf[cnt++] = lastbyte >> lastbits;
	}
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits>0) {
	cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}

/*___________________________________________________________________________
 |
 | sizeofints - calculate 'bitsize' of compressed ints
 |
 | given the number of small unsigned integers and the maximum value
 | return the number of bits needed to read or write them with the
 | routines receiveints and sendints. You need this parameter when
 | calling these routines. Note that for many calls I can use
 | the variable 'smallidx' which is exactly the number of bits, and
 | So I don't need to call 'sizeofints for those calls.
*/

static int sizeofints( const int num_of_ints, unsigned int sizes[]) {
    int bytes[32];
    int num_of_bits, i, num;
    int num_of_bytes, bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++) {	
	tmp = 0;
	for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
	    tmp = bytes[bytecnt] * sizes[i] + tmp;
	    bytes[bytecnt] = tmp & 0xff;
	    tmp >>= 8;
	}
	while (tmp != 0) {
	    bytes[bytecnt++] = tmp & 0xff;
	    tmp >>= 8;
	}
	num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
	num_of_bits++;
	num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;
}
    
/*____________________________________________________________________________
 |
 | sendints - send a small set of small integers in compressed format
 |
 | this routine is used internally by xdr3dfcoord, to send a set of
 | small integers to the buffer. 
 | Multiplication with fixed (specified maximum ) sizes is used to get
 | to one big, multibyte interger. Allthough the routine could be
 | modified to handle sizes bigger than 16777216, ot more than just
 | a few integers, this is not done, because the gain in compression
 | isn't worth the effort. Note that overflowing the multiplication
 | or the byte buffer (32 bytes) is unchecked and causes bad results.
 |
 */
 
static void sendints(int buf[], const int num_of_ints, const int num_of_bits,
	unsigned int sizes[], int nums[]) {
    unsigned int bytes[32], bytecnt;
    int i, num_of_bytes;
    unsigned int tmp;
    bytes[0] = 0;
    num_of_bytes = 0;
    for (i = 0; i < num_of_ints; i++) {
	if (nums[i] < 0 || nums[i]>= sizes[i]) {
	    fprintf(stderr,"major breakdown in sendints num %d doesn't "
		    "match size %d\n", nums[i], sizes[i]);
	    exit(1);
	}
	bytecnt = 0;
	tmp = nums[i];
	for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
	    tmp = bytes[bytecnt] * sizes[i] + tmp;
	    bytes[bytecnt] = tmp & 0xff;
	    tmp >>= 8;
	}
	while (tmp != 0) {
	    bytes[bytecnt++] = tmp & 0xff;
	    tmp >>= 8;
	}
	num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8) {
	for (i = 0; i < num_of_bytes; i++) {
	    sendbits(buf, 8, bytes[i]);
	}
	sendbits(buf, num_of_bits - num_of_bytes * 8, 0);
    } else {
	for (i = 0; i < num_of_bytes-1; i++) {
	    sendbits(buf, 8, bytes[i]);
	}
	sendbits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
}


/*___________________________________________________________________________
 |
 | receivebits - decode number from buf using specified number of bits
 | 
 | extract the number of bits from the array buf and construct an integer
 | from it. Return that value.
 |
*/

static int receivebits(int buf[], int num_of_bits) {

    int cnt, num; 
    unsigned int lastbits, lastbyte;
    unsigned char * cbuf;
    int mask = (1 << num_of_bits) -1;
    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];
    
    num = 0;
    while (num_of_bits >= 8) {
	lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
	num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
	num_of_bits -=8;
    }
    if (num_of_bits > 0) {
	if (lastbits < num_of_bits) {
	    lastbits += 8;
	    lastbyte = (lastbyte << 8) | cbuf[cnt++];
	}
	lastbits -= num_of_bits;
	num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num; 
}

/*____________________________________________________________________________
 |
 | receiveints - decode 'small' integers from the buf array
 |
 | this routine is the inverse from sendints() and decodes the small integers
 | written to buf by calculating the remainder and doing divisions with
 | the given sizes[]. You need to specify the total number of bits to be
 | used from buf in num_of_bits.
 |
*/

static void receiveints(int buf[], const int num_of_ints, int num_of_bits,
	unsigned int sizes[], int nums[]) {
    int bytes[32];
    int i, j, num_of_bytes, p, num;
    
    bytes[0] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8) {
	bytes[num_of_bytes++] = receivebits(buf, 8);
	num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
	bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
    }
    for (i = num_of_ints-1; i >=0; i--) {
	num = 0;
	for (j = num_of_bytes-1; j >=0; j--) {
	    num = (num << 8) | bytes[j];
	    p = num / sizes[i];
	    bytes[j] = p;
	    num = num - p * sizes[i];
	}
	nums[i] = num;
    }
}
    
/*____________________________________________________________________________
 |
 | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
 |
 | this routine reads or writes (depending on how you opened the file with
 | xdropen() ) a large number of 3d coordinates (stored in *fp).
 | The number of coordinates triplets to write is given by *size. On
 | read this number may be zero, in which case it reads as many as were written
 | or it may specify the number if triplets to read (which should match the
 | number written).
 | Compression is achieved by first converting all floating numbers to integer
 | using multiplication by *precision and rounding to the nearest integer.
 | Then the minimum and maximum value are calculated to determine the range.
 | The limited range of integers so found, is used to compress the coordinates.
 | In addition the differences between succesive coordinates is calculated.
 | If the difference happens to be 'small' then only the difference is saved,
 | compressing the data even more. The notion of 'small' is changed dynamically
 | and is enlarged or reduced whenever needed or possible.
 | Extra compression is achieved in the case of GROMOS and coordinates of
 | water molecules. GROMOS first writes out the Oxygen possition, followed by
 | the two hydrogens. In order to make the differences smaller (and there by
 | compression the data better) the order is changed into first one hydrogen
 | then the oxygen, followed by the other hydrogen. The is rather special, but
 | it shouldn't harm in the general case.
 |
 */
 
int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision) {
    

    static int *ip = NULL;
    static int oldsize;
    static int *buf;

    int minint[3], maxint[3], mindiff, *lip, diff;
    int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int minidx, maxidx;
    unsigned sizeint[3], sizesmall[3], size3, *luip;
    int flag, k;
    int small, smaller, larger, i, is_small, is_smaller, run, prevrun;
    float * lfp;
    int tmp, *thiscoord, tmpcoord[30], prevcoord[3];
    int bufsize, xdrid, lsize;
    unsigned int bitsize;
    float inv_precision;
    
    /* find out if xdrs is opened for reading or for writing */
    xdrid = 0;
    while (xdridptr[xdrid] != xdrs) {
	xdrid++;
	if (xdrid >= MAXID) {
	    fprintf(stderr, "xdr error. no open xdr stream\n");
	    exit (1);
	}
    }
    if (xdrmodes[xdrid] == 'w') {

	/* xdrs is open for writing */

	size3 = *size * 3;
	xdr_int(xdrs, size);
	/* when the number of coordinates is small, don't try to compress; just
	 * write them as floats using xdr_vector
	 */
	if (*size <= 9 ) {
	    return (xdr_vector(xdrs, (char *) fp, size3, sizeof(*fp),
		(xdrproc_t)xdr_float));
	}
	
	xdr_float(xdrs, precision);
	if (ip == NULL) {
	    ip = (int *)malloc(size3 * sizeof(*ip));
	    if (ip == NULL) {
		return 0;
	    }
	    bufsize = size3 * 1.2;
	    buf = (int *)malloc(bufsize * sizeof(*buf));
	    if (buf == NULL) {
		return 0;
	    }
	    oldsize = *size;
	} else if (*size > oldsize) {
	    ip = (int *)realloc(ip, size3 * sizeof(*ip));
	    bufsize = size3 * 1.2;
	    buf = (int *)realloc(buf, bufsize * sizeof(*buf));
	    oldsize = *size;
	}
	/* buf[0-2] are special and do not contain actual data */
	buf[0] = buf[1] = buf[2] = 0;
	minint[0] = minint[1] = minint[2] = INT_MAX;
	maxint[0] = maxint[1] = maxint[2] = INT_MIN;
	prevrun = -1;
	lfp = fp;
	lip = ip;
	mindiff = INT_MAX;
	if (fabs(*lfp) > MAXABS) {
	    /* scaling would cause overflow */
	    return 0;
	}
	oldlint1 = oldlint2 = oldlint3 = 0;
	while(lfp < fp + size3 ) {
	    if (fabs(*lfp) > MAXABS) {
		/* scaling would cause overflow */
		return 0;
	    }
	    /* find nearest integer */
	    if (*lfp >= 0.0)
		lint1 = *lfp * *precision + 0.5;
	    else
		lint1 = *lfp * *precision - 0.5;
	    if (lint1 < minint[0]) minint[0] = lint1;
	    if (lint1 > maxint[0]) maxint[0] = lint1;
	    *lip++ = lint1;
	    lfp++;
	    if (fabs(*lfp) > MAXABS) {
		/* scaling would cause overflow */
		return 0;
	    }
	    /* find nearest integer */
	    if (*lfp >= 0.0)
		lint2 = *lfp * *precision + 0.5;
	    else
		lint2 = *lfp * *precision - 0.5;
	    if (lint2 < minint[1]) minint[1] = lint2;
	    if (lint2 > maxint[1]) maxint[1] = lint2;
	    *lip++ = lint2;
	    lfp++;
	    if (fabs(*lfp) > MAXABS) {
		/* scaling would cause overflow */
		return 0;
	    }
	    /* find nearest integer */
	    if (*lfp >= 0.0)
		lint3 = *lfp * *precision + 0.5;
	    else
		lint3 = *lfp * *precision - 0.5;
	    if (lint3 < minint[2]) minint[2] = lint3;
	    if (lint3 > maxint[2]) maxint[2] = lint3;
	    *lip++ = lint3;
	    lfp++;
	    diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
	    if (diff < mindiff && lfp > fp + 3)
		mindiff = diff;
	    oldlint1 = lint1;
	    oldlint2 = lint2;
	    oldlint3 = lint3;
	}
	xdr_int(xdrs, &(minint[0]));
	xdr_int(xdrs, &(minint[1]));
	xdr_int(xdrs, &(minint[2]));
	
	xdr_int(xdrs, &(maxint[0]));
	xdr_int(xdrs, &(maxint[1]));
	xdr_int(xdrs, &(maxint[2]));
	
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	bitsize = sizeofints(3, sizeint);
	lip = ip;
	luip = (unsigned int *) ip;
	smallidx = FIRSTIDX;
	while (smallidx < LASTIDX && magicints[smallidx] < mindiff) {
	    smallidx++;
	}
	xdr_int(xdrs, &smallidx);
	maxidx = MIN(LASTIDX, smallidx + 8) ;
	minidx = maxidx - 8; /* often this equal smallidx */
	smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
	small = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
	larger = magicints[maxidx] / 2;
	i = 0;
	while (i < *size) {
	    is_small = 0;
	    thiscoord = (int *)(luip) + i * 3;
	    if (smallidx < maxidx && i >= 1 &&
		    abs(thiscoord[0] - prevcoord[0]) < larger &&
		    abs(thiscoord[1] - prevcoord[1]) < larger &&
		    abs(thiscoord[2] - prevcoord[2]) < larger) {
		is_smaller = 1;
	    } else if (smallidx > minidx) {
		is_smaller = -1;
	    } else {
		is_smaller = 0;
	    }
	    if (i + 1 < *size) {
		if (abs(thiscoord[0] - thiscoord[3]) < small &&
			abs(thiscoord[1] - thiscoord[4]) < small &&
			abs(thiscoord[2] - thiscoord[5]) < small) {
		    /* interchange first with second atom for better
		     * compression of water molecules
		     */
		    tmp = thiscoord[0]; thiscoord[0] = thiscoord[3];
			thiscoord[3] = tmp;
		    tmp = thiscoord[1]; thiscoord[1] = thiscoord[4];
			thiscoord[4] = tmp;
		    tmp = thiscoord[2]; thiscoord[2] = thiscoord[5];
			thiscoord[5] = tmp;
		    is_small = 1;
		}
    
	    }
	    tmpcoord[0] = thiscoord[0] - minint[0];
	    tmpcoord[1] = thiscoord[1] - minint[1];
	    tmpcoord[2] = thiscoord[2] - minint[2];
	    sendints(buf, 3, bitsize, sizeint, tmpcoord);
	    prevcoord[0] = thiscoord[0];
	    prevcoord[1] = thiscoord[1];
	    prevcoord[2] = thiscoord[2];
	    thiscoord = thiscoord + 3;
	    i++;
	    
	    run = 0;
	    if (is_small == 0 && is_smaller == -1)
		is_smaller = 0;
	    while (is_small && run < 8*3) {
		if (is_smaller == -1 && (
			SQR(thiscoord[0] - prevcoord[0]) +
			SQR(thiscoord[1] - prevcoord[1]) +
			SQR(thiscoord[2] - prevcoord[2]) >= smaller * smaller)) {
		    is_smaller = 0;
		}

		tmpcoord[run++] = thiscoord[0] - prevcoord[0] + small;
		tmpcoord[run++] = thiscoord[1] - prevcoord[1] + small;
		tmpcoord[run++] = thiscoord[2] - prevcoord[2] + small;
		
		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];

		i++;
		thiscoord = thiscoord + 3;
		is_small = 0;
		if (i < *size &&
			abs(thiscoord[0] - prevcoord[0]) < small &&
			abs(thiscoord[1] - prevcoord[1]) < small &&
			abs(thiscoord[2] - prevcoord[2]) < small) {
		    is_small = 1;
		}
	    }
	    if (run != prevrun || is_smaller != 0) {
		prevrun = run;
		sendbits(buf, 1, 1); /* flag the change in run-length */
		sendbits(buf, 5, run+is_smaller+1);
	    } else {
		sendbits(buf, 1, 0); /* flag the fact that runlength did not change */
	    }
	    for (k=0; k < run; k+=3) {
		sendints(buf, 3, smallidx, sizesmall, &tmpcoord[k]);	
	    }
	    if (is_smaller != 0) {
		smallidx += is_smaller;
		if (is_smaller < 0) {
		    small = smaller;
		    smaller = magicints[smallidx-1] / 2;
		} else {
		    smaller = small;
		    small = magicints[smallidx] / 2;
		}
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
	    }
	}
	if (buf[1] != 0) buf[0]++;;
	xdr_int(xdrs, &(buf[0])); /* buf[0] holds the length in bytes */
	return (xdr_opaque(xdrs, &(buf[3]), buf[0]));
    } else {
	
	/* xdrs is open for reading */
	
	xdr_int(xdrs, &lsize);
	if (*size != 0 && lsize != *size) {
	    fprintf(stderr, "wrong number of coordinates in xdr3dfcoor; "
		    "%d arg vs %d in file", *size, lsize);
	}
	*size = lsize;
	if (*size <= 9) {
	    return (xdr_vector(xdrs, (char *) fp, size3, sizeof(*fp),
		(xdrproc_t)xdr_float));
	}
	size3 = *size * 3;
	xdr_float(xdrs, precision);
	if (ip == NULL) {
	    ip = (int *)malloc(size3 * sizeof(*ip));
	    if (ip == NULL) {
		return 0;
	    }
	    bufsize = size3 * 1.2;
	    buf = (int *)malloc(bufsize * sizeof(*buf));
	    if (buf == NULL) {
		return 0;
	    }
	    oldsize = *size;
	} else if (*size > oldsize) {
	    ip = (int *)realloc(ip, size3 * sizeof(*ip));
	    bufsize = size3 * 1.2;
	    buf = (int *)realloc(buf, bufsize * sizeof(*buf));
	    oldsize = *size;
	}
	buf[0] = buf[1] = buf[2] = 0;
	
	xdr_int(xdrs, &(minint[0]));
	xdr_int(xdrs, &(minint[1]));
	xdr_int(xdrs, &(minint[2]));

	xdr_int(xdrs, &(maxint[0]));
	xdr_int(xdrs, &(maxint[1]));
	xdr_int(xdrs, &(maxint[2]));
		
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	bitsize = sizeofints(3, sizeint);
	
	xdr_int(xdrs, &smallidx);
	maxidx = MIN(LASTIDX, smallidx + 8) ;
	minidx = maxidx - 8; /* often this equal smallidx */
	smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
	small = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	larger = magicints[maxidx];
    	
	xdr_int(xdrs, &(buf[0])); /* buf[0] holds the length in bytes */
	xdr_opaque(xdrs, &(buf[3]), buf[0]);
	buf[0] = buf[1] = buf[2] = 0;
	
	lfp = fp;
	inv_precision = 1.0 / * precision;
	i = 0;
	lip = ip;
	while ( i < lsize ) {
	    thiscoord = (int *)(lip) + i * 3;

	    receiveints(buf, 3, bitsize, sizeint, thiscoord);
	    i++;
	    thiscoord[0] += minint[0];
	    thiscoord[1] += minint[1];
	    thiscoord[2] += minint[2];
	    
	    prevcoord[0] = thiscoord[0];
	    prevcoord[1] = thiscoord[1];
	    prevcoord[2] = thiscoord[2];
	    
	   
	    flag = receivebits(buf, 1);
	    is_smaller = 0;
	    if (flag == 1) {
		run = receivebits(buf, 5);
		is_smaller = run % 3;
		run -= is_smaller;
		is_smaller--;
	    }
	    if (run > 0) {
		thiscoord += 3;
		for (k = 0; k < run; k+=3) {
		    receiveints(buf, 3, smallidx, sizesmall, thiscoord);
		    i++;
		    thiscoord[0] += prevcoord[0] - small;
		    thiscoord[1] += prevcoord[1] - small;
		    thiscoord[2] += prevcoord[2] - small;
		    if (k == 0) {
			/* interchange first with second atom for better
			 * compression of water molecules
			 */
			tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
				prevcoord[0] = tmp;
			tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
				prevcoord[1] = tmp;
			tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
				prevcoord[2] = tmp;
			*lfp++ = prevcoord[0] * inv_precision;
			*lfp++ = prevcoord[1] * inv_precision;
			*lfp++ = prevcoord[2] * inv_precision;
		    } else {
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
		    }
		    *lfp++ = thiscoord[0] * inv_precision;
		    *lfp++ = thiscoord[1] * inv_precision;
		    *lfp++ = thiscoord[2] * inv_precision;
		}
	    } else {
		*lfp++ = thiscoord[0] * inv_precision;
		*lfp++ = thiscoord[1] * inv_precision;
		*lfp++ = thiscoord[2] * inv_precision;		
	    }
	    smallidx += is_smaller;
	    if (is_smaller < 0) {
		small = smaller;
		if (smallidx > FIRSTIDX) {
		    smaller = magicints[smallidx - 1] /2;
		} else {
		    smaller = 0;
		}
	    } else if (is_smaller > 0) {
		smaller = small;
		small = magicints[smallidx] / 2;
	    }
	    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	}
	return 1;
    }
}


   
