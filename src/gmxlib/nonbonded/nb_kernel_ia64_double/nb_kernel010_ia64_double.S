/*
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/*
 * The ia64-assembly Gromacs inner loops would not have been
 * possible without a lot of support and optimization suggestions 
 * from John Worley at Hewlett-Packard.
 */

/* Each thread locks a counter and grabs a couple of neighborlists.
 * Available sizes for this chunk: 1,2,4,8, or 16 
 */
#define THREAD_CHUNK_SIZE       8
#define JJNR_PREFETCH_DISTANCE  128

#define	zero	r0	/* permanent zero					*/
#define	gp		r1	/* global data pointer				*/

#define	at0		r2	/* temp, target of addi				*/
#define	at1		r3	/* temp, target of addi				*/

#define	S0		r4	/* callee saves register			*/
#define	S1		r5	/* callee saves register			*/
#define	S2		r6	/* callee saves register			*/
#define	S3		r7	/* callee saves register			*/

#define	v0		r8	/* 1st fixed point return value/ptr	*/
#define	v1		r9	/* 2nd fixed return value/ptr		*/
#define	v2		r10	/* 3rd fixed return value/ptr		*/
#define	v3		r11	/* 4th fixed return value/ptr		*/

#define	sp		r12	/* memory stack pointer				*/
#define	tp		r13	/* thread pointer					*/

#define	t0		r14	/* caller saves register			*/
#define	t1		r15	/* caller saves register			*/
#define	t2		r16	/* caller saves register			*/
#define	t3		r17	/* caller saves register			*/
#define	t4		r18	/* caller saves register			*/
#define	t5		r19	/* caller saves register			*/

#define	t6		r20	/* caller saves register			*/
#define	t7		r21	/* caller saves register			*/
#define	t8		r22	/* caller saves register			*/
#define	t9		r23	/* caller saves register			*/
#define	t10		r24	/* caller saves register			*/
#define	t11		r25	/* caller saves register			*/
#define	t12		r26	/* caller saves register			*/
#define	t13		r27	/* caller saves register			*/
#define	t14		r28	/* caller saves register			*/
#define	t15		r29	/* caller saves register			*/
#define	t16		r30	/* caller saves register			*/
#define	t17		r31	/* caller saves register			*/


//	ia64 Floating-point register definitions
#define	fZero	f0	/* permanent floating point 0.0		*/
#define	fOne	f1	/* permanent floating point 1.0		*/

#define	fs0		f2	/* callee saves register			*/
#define	fs1		f3	/* callee saves register			*/
#define	fs2		f4	/* callee saves register			*/
#define	fs3		f5	/* callee saves register			*/
	
#define	ft0		f6	/* caller saves register			*/
#define	ft1		f7	/* caller saves register			*/

#define	fa0		f8	/* argument register 0				*/
#define	fa1		f9	/* argument register 1				*/
#define	fa2		f10	/* argument register 2				*/
#define	fa3		f11	/* argument register 3				*/
#define	fa4		f12	/* argument register 4				*/
#define	fa5		f13	/* argument register 5				*/
#define	fa6		f14	/* argument register 6				*/
#define	fa7		f15	/* argument register 7				*/

#define	fv0		f8	/* return value register 0			*/
#define	fv1		f9	/* return value register 1			*/
#define	fv2		f10	/* return value register 2			*/
#define	fv3		f11	/* return value register 3			*/
#define	fv4		f12	/* return value register 4			*/
#define	fv5		f13	/* return value register 5			*/
#define	fv6		f14	/* return value register 6			*/
#define	fv7		f15	/* return value register 7			*/

#define	fs4		f16	/* callee saves register			*/
#define	fs5		f17	/* callee saves register			*/
#define	fs6		f18	/* callee saves register			*/
#define	fs7		f19	/* callee saves register			*/
#define	fs8		f20	/* callee saves register			*/
#define	fs9		f21	/* callee saves register			*/
#define	fs10	f22	/* callee saves register			*/
#define	fs11	f23	/* callee saves register			*/

#define	fs12	f24	/* callee saves register			*/
#define	fs13	f25	/* callee saves register			*/
#define	fs14	f26	/* callee saves register			*/
#define	fs15	f27	/* callee saves register			*/
#define	fs16	f28	/* callee saves register			*/
#define	fs17	f29	/* callee saves register			*/
#define	fs18	f30	/* callee saves register			*/
#define	fs19	f31	/* callee saves register			*/

// ia64 predicate register definitions
#define	pone	p0	/* permanent one predicate			*/
#define	pTrue	p0	/* permanent one predicate			*/

#define	ps0		p1	/* callee saves predicate			*/
#define	ps1		p2	/* callee saves predicate			*/
#define	ps2		p3	/* callee saves predicate			*/
#define	ps3		p4	/* callee saves predicate			*/
#define	ps4		p5	/* callee saves predicate			*/

#define	pt0		p6	/* caller saves predicate			*/
#define	pt1		p7	/* caller saves predicate			*/
#define	pt2		p8	/* caller saves predicate			*/
#define	pt3		p9	/* caller saves predicate			*/
#define	pt4		p10	/* caller saves predicate			*/
#define	pt5		p11	/* caller saves predicate			*/
#define	pt6		p12	/* caller saves predicate			*/
#define	pt7		p13	/* caller saves predicate			*/
#define	pt8		p14	/* caller saves predicate			*/
#define	pt9		p15	/* caller saves predicate			*/

// ia64 branch register definitions
#define	rb		b0	/* return link						*/

#define	bs0		b1	/* callee saves branch register		*/
#define	bs1		b2	/* callee saves branch register		*/
#define	bs2		b3	/* callee saves branch register		*/
#define	bs3		b4	/* callee saves branch register		*/
#define	bs4		b5	/* callee saves branch register		*/
	
#define	bt0		b6	/* caller saves branch register		*/
#define	bt1		b7	/* caller saves branch register		*/
	
		
.text

#define	CHARGE		t10
#define	FACTION		t9
#define	FActII		loc8
#define	FActIX		fs1
#define	FActIY		fs2
#define	FActIZ		fs3
#define	FIX			fa6
#define	FIY			fa7
#define	FIZ			fs4
#define	FSHIFT		t6
#define	FShiftIS	loc9
#define	FShiftX		fa3
#define	FShiftY		fa4
#define	FShiftZ		fa5
#define	Facel		fa3
#define	ICharge		fa7
#define	InnerCnt	t17
#define	II			t13
#define	II3			in7
#define	IQ			fa7
#define	IS			t12
#define	IS3			in6
#define	IX			fa0
#define	IY			fa1
#define	IZ			fa2
#define	In_FSHIFT	in6
#define	In_GID		in7
#define	In_IINR		in1
#define	In_JINDEX	in2
#define	In_JJNR		in3
#define	In_NRI		in0
#define	In_SHIFT	in4
#define	In_SHIFTVEC	in5
#define NRI			loc12
#define IINR		loc13
#define JINDEX		loc14
#define JJNR		loc15
#define SHIFT		loc16
#define GID			loc17
#define COUNT		loc18
#define	JX			DX[0]
#define	JY			DY[0]
#define	JZ			DZ[0]
#define	LCSave		at0
#define	NJ0			t14
#define	NJ1			t15
#define	POSITION	t8
#define	PRSave		at1
#define	PosX		f88
#define	PosY		f89
#define	PosZ		f90
#define	QCharge		Charge[0]
#define	QQ			Charge[3]
#define	SHIFTVEC	t5
#define	VC			t11
#define	VCPtr		ggid
#define	VNBTotal	fs0
#define	VCoul		Charge[3]
#define	VCoul2		Charge[4]
#define	argPtr		loc23
#define	chargePtr	v2
#define	Tmp1		t0
#define	Tmp2		t17
#define	Tmp3		loc11
#define	Tmp4		t2
#define	fHALF		ft0
#define fSIX        ft0
#define fTWELVE     ft1
#define	fillP0		v0
#define	fillP1		v1
#define NN0			t0
#define NN1			loc11
#define	ggid		loc10
#define	gidPtr		t7
#define	iinrPtr		t1
#define	jindexPtr	t2
#define	jjnrPtr		t3
#define	jnr			t16
#define	jnr3		v0
#define	nriCount	t0
#define	pCont		pt0
#define	pDone		pt1
#define	pJJNR		pt2
#define	pMore		pt3
#define	pLast		pt4
#define	posPtr		v3
#define	shX			fa0
#define	shY			fa1
#define	shZ			fa2
#define	shiftPtr	t4
#define	shiftVPtr	v1
#define	spillPtr	v0
#define	spillPtr2   t0
#define	xPFS		at0
#define TYPE        loc19
#define NTYPE       loc20
#define typePtr     loc21
#define NBFP		loc22
#define NTI     	loc24
#define VNBPtr      loc25
#define Nouter      loc26
#define Ninner      loc27
#define OuterIter   loc28
#define InnerIter   loc29
#define VNB         loc30
		
#define	_NINPUTS	8
#define	_NLOCALS	32
#define	_NOUTPUT	0
#define	_NROTATE	16

	.regstk	8, 32, 0, 16
	.rotr	FActPtr[7], TypeJ[3]
	.rotf	DX[7], DY[7], DZ[7], FActX[7], FActY[7], FActZ[7], RSqr[3], RInv2Err[2], FScalar[3], RInv2[4], RInv6[2], C6[4], C12[4]
	.rotp	pPipe[7]

#define	PIPE_DEPTH	7

#define	EXP(n)					(0xffff + (n))

#define	POS_STK_OFFSET			0x10
#define	FACTION_STK_OFFSET		0x18
#define	CHARGE_STK_OFFSET		0x20
#define	FACEL_STK_OFFSET		0x28
#define	KRF_STK_OFFSET			0x30
#define	CRF_STK_OFFSET			0x38
#define	VC_STK_OFFSET			0x40
#define	TYPE_STK_OFFSET			0x48
#define	NTYPE_STK_OFFSET		0x50
#define	NBFP_STK_OFFSET			0x58
#define	VNB_STK_OFFSET			0x60
#define	TABSCALE_STK_OFFSET		0x68
#define	VFTAB_STK_OFFSET		0x70
#define	INVSQRTA_STK_OFFSET		0x78
#define	DVDA_STK_OFFSET			0x80
#define	GBTABSCALE_STK_OFFSET		0x88
#define	GBTAB_STK_OFFSET		0x90
#define	NTHREADS_STK_OFFSET		0x98
#define	COUNT_STK_OFFSET		0xA0
#define MTX_STK_OFFSET			0xA8
#define OUTERITER_STK_OFFSET		0xB0
#define INNERITER_STK_OFFSET		0xB8
#define WORK_STK_OFFSET     		0xC0


	.global nb_kernel010_ia64_double
	.proc	nb_kernel010_ia64_double
	.align	32

nb_kernel010_ia64_double:
//	INIT 1
	{	.mmi
		alloc			xPFS = ar.pfs, _NINPUTS, _NLOCALS, _NOUTPUT, _NROTATE
		mov				spillPtr = sp
		mov				Tmp1 = EXP(-1)
	}
	{	.mfi
		nop				0x0
		fadd			fSIX = fOne, fOne
		add				argPtr = COUNT_STK_OFFSET, sp
	} ;;
//	INIT 2
	{	.mlx									
		nop				0x0
		movl			Tmp2 = 0x40c00000   // 6.0
	}
	{	.mii
		stf.spill		[spillPtr] = fs0, -16
		mov				LCSave = ar.lc
		add				sp = -5 * 16, sp
	} ;;
//	INIT 3	
	{	.mmi					
		stf.spill		[spillPtr] = fs1, -16
		ld8				COUNT = [argPtr], POS_STK_OFFSET - COUNT_STK_OFFSET
		mov				PRSave	= pr
	} ;;
//	INIT 4	
	{	.mfi
		ld8				POSITION = [argPtr], FACTION_STK_OFFSET - POS_STK_OFFSET
		fma				fSIX = fSIX, fOne, fOne  // 3.0
		nop				0x0
	}
	{	.mmi									
		nop				0x0
		stf.spill		[spillPtr] = fs2, -16
		mov				pr.rot	= 0x0
	} ;;
//  INIT 5	
	{	.mmi								
		ld8			FACTION = [argPtr], NTYPE_STK_OFFSET - FACTION_STK_OFFSET
		nop			0x0
		mov			SHIFTVEC   = In_SHIFTVEC
	}
   	{	.mfi
		stf.spill		[spillPtr] = fs3, -16
		nop			0x0
		nop			0x0
	} ;;
//  INIT 6	
	{	.mfi
		ld8			NTYPE = [argPtr], TYPE_STK_OFFSET - NTYPE_STK_OFFSET
		fadd			fSIX = fSIX, fSIX
		mov			FSHIFT   = In_FSHIFT
	} 
	{	.mfi
		stf.spill		[spillPtr] = fs4, -16
		nop			0x0
		mov			GID = In_GID
	} ;;
//  INIT 7	
	{	.mmf
		ld8			TYPE = [argPtr], NBFP_STK_OFFSET - TYPE_STK_OFFSET
		mov			SHIFT = In_SHIFT
		nop			0x0
	}
	{	.mmi
		stf.spill		[spillPtr] = fs5, -16
		ld4			NRI = [In_NRI]
		nop			0x0
	} ;;
//  INIT 8	
	{	.mmf
		ld8			NBFP = [argPtr], VNB_STK_OFFSET - NBFP_STK_OFFSET
		mov			JJNR = In_JJNR
		nop			0x0
	}	
	{	.mii
		nop			0x0
		mov			IINR = In_IINR
		nop			0x0
	} ;;
//  INIT 9		
	{	.mmf
		ld8			VNB = [argPtr],	OUTERITER_STK_OFFSET - VNB_STK_OFFSET
		mov			JINDEX = In_JINDEX
		fadd			fTWELVE = fSIX, fSIX
	} ;; 
//  INIT 10
	{	.mfi
		ld8			OuterIter = [argPtr], INNERITER_STK_OFFSET - OUTERITER_STK_OFFSET
		nop			0x0
		mov			Nouter = 0
	} ;;
//  INIT 11
	{	.mmi
		ld8			InnerIter = [argPtr]
		ld4			NTYPE = [NTYPE]
		mov			Ninner = 0
	} ;;
//  18 bundles used for init - still aligned.




	
threadLoop:
//  THREAD PROLOGUE 1	
	{	.mfi		
		fetchadd4.rel	NN0 = [COUNT], THREAD_CHUNK_SIZE
		nop				0x0
		nop				0x0
	}
	{	.mfi
		setf.sig		f33 = NTYPE
		nop				0x0
		nop				0x0
	} ;;    
//  THREAD PROLOGUE 2 - at least 12 cycle latency hole before this bundle (fetchadd4)
	{	.mmi		
		cmp.lt			pCont, pDone = NN0, NRI
		shladd			gidPtr = NN0, 2, GID
		adds			NN1 = THREAD_CHUNK_SIZE, NN0
	}
	{	.mmi
		shladd			jindexPtr = NN0, 2, JINDEX
		shladd   		shiftPtr  = NN0, 2, SHIFT
		shladd			iinrPtr   = NN0, 2, IINR
	} ;; 
//  THREAD PROLOGUE 3 	
	{ .mmi				
	(pCont) ld4			II = [iinrPtr], 4
	(pCont) ld4			IS = [shiftPtr], 4
		cmp.ge			pLast, pMore = NN1, NRI
	}
	{ .mib
	(pCont) ld4			NJ0 = [jindexPtr], 4
	(pCont) adds		Tmp2 = 1, NN0
	(pDone) br.cond.spnt.few finish
	} ;; 		
//  THREAD PROLOGUE 4	
	{ .mmi				
		ld4				ggid = [gidPtr], 4
		shladd			II3 = II, 1, II
		shladd			IS3 = IS, 1, IS
	}
	{ .mfi
		ld4				NJ1 = [jindexPtr], 4
		nop				0x0
		shladd			jjnrPtr = NJ0, 2, JJNR
	} ;;
//  THREAD PROLOGUE 5	
	{ .mmi
		cmp.lt			pCont, pDone = Tmp2, NRI						
		shladd			FShiftIS  = IS3, 3, FSHIFT
		shladd			typePtr = II, 2, TYPE
	}	
	{ .mmi
		shladd			posPtr    = II3, 3, POSITION
		shladd			FActII    = II3, 3, FACTION
		shladd			shiftVPtr = IS3, 3, SHIFTVEC	
	} ;;
//  THREAD PROLOGUE 6	
	{ .mmi				
		ld4				jnr = [jjnrPtr], 4
	(pCont)	ld4			IS = [shiftPtr], 4
		nop				0x0
	}	
	{ .mmi
	(pCont) ld4			II = [iinrPtr], 4
			ld4			NTI = [typePtr]
	(pLast)	mov			NN1 = NRI
	} ;;
//  12 bundles in thread prologue - still aligned







outerLoop:
	//	At this point in the outer loop, the following values are ready
	//
	//		FActII		Pointer to FACTION XYZ for II
	//		FShiftIS	Pointer to FSHIFT XYZ for IS
	//		shiftVPtr	Pointer to current shift XYZ values
	//		posPtr		Pointer to current XYZ position
	//		chargePtr	Pointer to current atom charge
	//		ggid		Index for Vc array
	//		jjnr		Pointer to next neighbor index
	//		jnr			Current jnr value
	//		NJ0, NJ1	Bounds of current neighbor list
	//
	//	Load up all the floating-point values (yes, McKinley can do 4 FP loads
	//	per cycle) and initialize the loop counters and predicates. Compute
	//	the initial position <x, y, z> and charge. If this isn't the last time
	//	through the loop, start loading the next value for NJ1 - we already
	//	moved the previous NJ1 -> NJ0.
//	OUTER PROLOGUE 1
	{	.mfi						
		nop 		0x0
		mov			FIX = f0
		add		Nouter = 1, Nouter
	}
	{	.mmf
		ldfd		shX = [shiftVPtr], 8
		ldfd		PosX = [posPtr], 8
		mov			FIY = f0
	} ;;
//	OUTER PROLOGUE 2
	{	.mmf
		setf.sig	f32 = NTI
		ldfd		shY = [shiftVPtr], 8
		nop			0x0
	}
	{	.mfi
		ldfd		PosY = [posPtr], 8
		nop			0x0
		sub			InnerCnt = NJ1, NJ0, 1
	} ;;

	{	.mmf						
		ldfd		shZ = [shiftVPtr]
		ldfd		PosZ = [posPtr]
		mov			FIZ = f0
	}
	{	.mmi
		ldfd		FShiftX = [FShiftIS], 8
		ldfd		FActIX = [FActII], 8
		shladd		VNBPtr = ggid, 3, VNB
	} ;;
//	OUTER PROLOGUE 4
	{	.mmf	
		ldfd		FShiftY = [FShiftIS], 8
		ldfd		FActIY = [FActII], 8
		xma.l		f32 = f32, f33, fZero
	};;
//	OUTER PROLOGUE 5
	{	.mmi
		ldfd		FActIZ = [FActII], -16
		ldfd		FShiftZ = [FShiftIS], -16
		mov			NJ0 = NJ1
	} ;;
//	OUTER PROLOGUE 6
	{	.mfi
		ldfd		VNBTotal = [VNBPtr]
		fadd		IX = shX, PosX
		add			NN0 = 1, NN0
	}
	{	.mfi
	(pCont)	ld4		NJ1 = [jindexPtr], 4
		//	This may seem strange, but we set the first stage of the
		//	pipe to execute this way because setting pr.rot doesn't take
		//	into account how much the predicates have rotated. If this is
		//	the first time through, we cleared all the pipeline predicates
		//	in the initialization. If not, flushing the pipeline set all
		//	the pipeline predicates to 0

		cmp.eq		pPipe[0], p0 = zero, zero
	} ;;
//	OUTER PROLOGUE 7
	{	.mfi		
		cmp.lt		pCont, pDone = NN0, NN1
		fadd		IY = shY, PosY
		mov		    ar.lc = InnerCnt
	} ;;
//	OUTER PROLOGUE 8
	{	.mfi		
		getf.sig	NTI = f32
		fadd		IZ = shZ, PosZ
		mov			ar.ec = PIPE_DEPTH
	} ;;
// 12 bundles in outer loop - still aligned.

	//	The inner loop is a 6-stage pipeline. The serial sequence of float ops
	//	is folded into a 17-cycle loop (17 * 2 = 34 float ops, one empty), 
    //  then divided
	//	into 5 stages.


	

innerLoop:
//	INNER LOOP 1
	{	.mfi	
	(pPipe[0])	shladd	chargePtr = jnr, 3, CHARGE
	(pPipe[1])	fsub	DX[1] = IX, DX[1]
	(pPipe[0])	shladd	jnr3 = jnr, 1, jnr
	}
	//	We march through jjnr[] sequentially, so it's usually a good idea
	//	to preload the next value. However, we don't want to do this if
	//	(1) we're in the epilogue or (2) this is the last time through and
	//	there are no more atoms to inspect. Thus, we keep track of the loop
	//	trip and use the logic below to see if we should load ahead

	.pred.rel "mutex", pCont, pDone
	{	.mfi
	(pCont)		cmp.ge	pJJNR, p0 = InnerCnt, zero
	(pPipe[3])	fma		RInv2[1] = RInv2[1], RInv2Err[1], RInv2[1] 
	(pDone)		cmp.gt	pJJNR, p0 = InnerCnt, zero
	} ;;
//	INNER LOOP 2
	{	.mfi	
				nop		0x0
	(pPipe[1])	fsub	DY[1] = IY, DY[1]
	(pPipe[0])	add		InnerCnt = -1, InnerCnt
	}
	{	.mfi
	(pPipe[0])	shladd	posPtr = jnr3, 3, POSITION
	(pPipe[5])	fma		FScalar[1] = fTWELVE, C12[3], FScalar[1]
	(pPipe[0])	shladd	FActPtr[0] = jnr3, 3, FACTION
	} ;;
//	INNER LOOP 3
	{	.mfi									
	(pPipe[0])	ldfd	JX = [posPtr], 8
	(pPipe[1])	fsub	DZ[1] = IZ, DZ[1]
	(pPipe[0])	shladd  TypeJ[0] = jnr, 2, TYPE
	}
	{  	.mfi
	(pJJNR)		ld4	jnr = [jjnrPtr], 4
	(pPipe[2])	frcpa 	RInv2[0], p0 = fOne, RSqr[1]	
	(pPipe[2])	add		TypeJ[2] = NTI, TypeJ[2]	
	} ;;
//	INNER LOOP 4
	{	.mfi	
	(pPipe[0])	ldfd	JY = [posPtr], 8
	(pPipe[4])	fmpy	RInv6[1] = RInv6[1], RInv2[2]
	(pPipe[0])	add	Ninner = 1, Ninner
	}
	{	.mfi
				nop		0x0
	(pPipe[5])	fsub	VNBTotal = VNBTotal, C6[3]
				nop		0x0
	} ;;
//	INNER LOOP 5
	{	.mfi									
	(pPipe[0])	ldfd	JZ = [posPtr], 8
	(pPipe[1])	fmpy	RSqr[0] = DX[1], DX[1]
	(pJJNR)     add     jjnrPtr = JJNR_PREFETCH_DISTANCE, jjnrPtr
	}
	{	.mfi
				nop		0x0
	(pPipe[3])	fnma	RInv2Err[1] = RSqr[2], RInv2[1], fOne
				nop		0x0
	} ;;
//	INNER LOOP 6
	{	.mfi	
	(pPipe[0])	ldfd	FActX[0] = [FActPtr[0]], 8
	(pPipe[5])	fmpy	FScalar[1] = RInv2[3], FScalar[1]	
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[6])	fma 	FIX = DX[6], FScalar[2], FIX
				nop		0x0
	} ;;
//	INNER LOOP 7
	{	.mfi									
	(pPipe[0])	ldfd	FActY[0] = [FActPtr[0]], 8
	(pPipe[2])	fnma	RInv2Err[0] = RInv2[0], RSqr[1], fOne
				nop		0x0
	}
	{	.mfi
	(pJJNR)     lfetch.nta  [jjnrPtr]
	(pPipe[6])	fma 	FIY = DY[6], FScalar[2], FIY
				nop		0x0
	} ;;
//	INNER LOOP 8
	{	.mfi	
	(pPipe[0])	ldfd	FActZ[0] = [FActPtr[0]], -16
	(pPipe[4])	fmpy	C6[2] = C6[2], RInv6[1]
	(pPipe[2])	shladd	TypeJ[2] = TypeJ[2], 4, NBFP
	}
	{	.mfi
				nop		0x0
	(pPipe[4])	fmpy	RInv6[1] = RInv6[1], RInv6[1]
				nop		0x0
	} ;;
//	INNER LOOP 9
	{	.mfi									
	(pPipe[2])	ldfd	C6[0] = [TypeJ[2]], 8
	(pPipe[1])	fma	RSqr[0] = DY[1], DY[1], RSqr[0]
	(pJJNR)     add     jjnrPtr = -JJNR_PREFETCH_DISTANCE, jjnrPtr
	}
	{	.mfi
	(pPipe[0])	ld4 	TypeJ[0] = [TypeJ[0]]
	(pPipe[3])	fma		RInv2[1] = RInv2[1], RInv2Err[1], RInv2[1] 
				nop		0x0
	} ;;
//	INNER LOOP 10
	{	.mfi		
	(pPipe[2])	ldfd	C12[0] = [TypeJ[2]]
	(pPipe[5])	fadd	VNBTotal = VNBTotal, C12[3]
				nop		0x0
	}
	{	.mfi
				nop		0x0				
	(pPipe[6])	fma 	FIZ = DZ[6], FScalar[2], FIZ
				nop		0x0
	} ;;
//	INNER LOOP 11
	{	.mfi		
				nop		0x0
	(pPipe[2])	fma		RInv2Err[0] = RInv2Err[0], RInv2Err[0], RInv2Err[0]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[5])	fnma	FActX[5] = FScalar[1], DX[5], FActX[5]
				nop		0x0
	} ;;
//	INNER LOOP 12
	{	.mfi
	(pPipe[6])	stfd	[FActPtr[6]] = FActX[6], 8	
	(pPipe[4])	fmpy	C12[2] = C12[2], RInv6[1]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[4])	fnma	FScalar[0] = fSIX, C6[2], fZero
				nop		0x0
	} ;;
//	INNER LOOP 13
	{	.mfi
	(pPipe[6])	stfd	[FActPtr[6]] = FActY[6], 8	
	(pPipe[1])	fma	RSqr[0] = DZ[1], DZ[1], RSqr[0]
				nop		0x0
	}
	{	.mfi
				nop		0x0
	(pPipe[3])	fmpy	RInv6[0] = RInv2[1], RInv2[1]
				nop		0x0
	} ;;
//	INNER LOOP 14
	{	.mfi
	(pPipe[6])	stfd	[FActPtr[6]] = FActZ[6], 8
	(pPipe[5])	fnma	FActY[5] = FScalar[1], DY[5], FActY[5]
				nop		0x0
	}
	{	.mfb
				nop		0x0
	(pPipe[5])	fnma	FActZ[5] = FScalar[1], DZ[5], FActZ[5]
				br.ctop.sptk.many 	innerLoop
	} ;;



// 	End of modulo-scheduled inner loop

	//	Having finshed the loop, we now compute various quantities to
	//	store. In paralllel, start computing computing some of the values
	//	for the next loop trip, if we're going there.

//	OUTER EPILOGUE 1
    {   .mfi
	(pCont)	shladd	typePtr = II, 2, TYPE
			nop		0x0
	(pCont)	shladd	II3 = II, 1, II
    }
	{	.mfi								
			nop		0x0
			nop		0x0
	(pCont)	shladd	IS3 = IS, 1, IS
    } ;;
//	OUTER EPILOGUE 2
    {   .mfi
	(pCont)	ld4	IS = [shiftPtr], 4
		fadd	FActIX = FActIX, FIX
		nop		0x0
	}
    {   .mmf
	(pCont)	setf.sig	f33 = NTYPE
	(pCont)	ld4	II = [iinrPtr] ,4
		fadd	FShiftX = FShiftX, FIX
	} ;;
// 	OUTER EPILOGUE 3
    {   .mfi
	(pCont)	ld4				NTI = [typePtr]	  	
		fadd	FActIY = FActIY, FIY
	(pCont)	shladd	shiftVPtr = IS3, 3, SHIFTVEC						
	} 
    {   .mfi
		nop 0x0
		fadd	FShiftY = FShiftY, FIY
	(pCont)	shladd	posPtr = II3, 3, POSITION
	} ;;
//	OUTER EPILOGUE 4
    {   .mfi
		nop 	0x0
		fadd	FActIZ = FActIZ, FIZ
		nop 	0x0
	} 
    {   .mfi
		nop 	0x0
		fadd	FShiftZ = FShiftZ, FIZ		
		nop 	0x0
	} ;;
//	OUTER EPILOGUE 5
	{	.mmi
		stfd	[FActII] = FActIX, 8
		stfd	[FShiftIS] = FShiftX, 8
		nop 	0x0
	}
    {   .mmi
		stfd    [VNBPtr] = VNBTotal
	(pCont)		ld4     ggid = [gidPtr], 4
		nop 	0x0
	} ;;
//	OUTER EPILOGUE 6
	{	.mmi
		stfd	[FActII] = FActIY, 8
		stfd	[FShiftIS] = FShiftY, 8
		nop 0x0
	} ;;
//	OUTER EPILOGUE 7
	{	.mmi
		stfd	[FActII] = FActIZ
		nop		0x0
	(pCont)	shladd	FActII = II3, 3, FACTION
	}
	{	.mib
		stfd	[FShiftIS] = FShiftZ
	(pCont)	shladd	FShiftIS = IS3, 3, FSHIFT
	(pCont)	br.cond.sptk.many	outerLoop
	} ;;




	// Finish if this was the last chunk, or do another thread-loop iteration
//  THREAD EPILOGUE 1
	{ .mib				
		nop				0x0
		nop				0x0
	(pMore) br.cond.sptk.many threadLoop
	} ;;
	






	//	Ready to exit - restore the floating-point registers we saved, the
	//	loop counter, and the predicates, then we're done. Note that the
	//	stack pointer has the address of the last saved FP register.

finish:
//  EXIT 1
	{	.mmi
		mov			fillP0 = sp
		add			fillP1 = 16, sp
		mov			ar.lc = LCSave
	}  
	{	.mmi
		st4			[OuterIter] = Nouter
		st4			[InnerIter] = Ninner
		nop			0x0
	} ;;
//  EXIT 2
	{	.mmi
		ldf.fill		fs5 = [fillP0], 32
		ldf.fill		fs4 = [fillP1], 32
		mov				pr = PRSave, 0x1ffff
	} ;;
//  EXIT 3
	{	.mmi
		ldf.fill		fs3 = [fillP0], 32
		ldf.fill		fs2 = [fillP1], 32
		add				sp = 5 * 16, sp
	} ;;
//  EXIT 4
	{	.mmb
		ldf.fill		fs1 = [fillP0]
		ldf.fill		fs0 = [fillP1]
		br.ret.sptk.few	rp
	} ;;

	.endp	 nb_kernel010_ia64_double



