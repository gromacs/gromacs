;#
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as

.section .text

	
.globl nb_kernel114_ia32_sse2
.globl _nb_kernel114_ia32_sse2
nb_kernel114_ia32_sse2:	
_nb_kernel114_ia32_sse2:	
.equiv          nb114_p_nri,            8
.equiv          nb114_iinr,             12
.equiv          nb114_jindex,           16
.equiv          nb114_jjnr,             20
.equiv          nb114_shift,            24
.equiv          nb114_shiftvec,         28
.equiv          nb114_fshift,           32
.equiv          nb114_gid,              36
.equiv          nb114_pos,              40
.equiv          nb114_faction,          44
.equiv          nb114_charge,           48
.equiv          nb114_p_facel,          52
.equiv          nb114_argkrf,           56
.equiv          nb114_argcrf,           60
.equiv          nb114_Vc,               64
.equiv          nb114_type,             68
.equiv          nb114_p_ntype,          72
.equiv          nb114_vdwparam,         76
.equiv          nb114_Vvdw,             80
.equiv          nb114_p_tabscale,       84
.equiv          nb114_VFtab,            88
.equiv          nb114_invsqrta,         92
.equiv          nb114_dvda,             96
.equiv          nb114_p_gbtabscale,     100
.equiv          nb114_GBtab,            104
.equiv          nb114_p_nthreads,       108
.equiv          nb114_count,            112
.equiv          nb114_mtx,              116
.equiv          nb114_outeriter,        120
.equiv          nb114_inneriter,        124
.equiv          nb114_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb114_ixO,              0
.equiv          nb114_iyO,              16
.equiv          nb114_izO,              32
.equiv          nb114_ixH1,             48
.equiv          nb114_iyH1,             64
.equiv          nb114_izH1,             80
.equiv          nb114_ixH2,             96
.equiv          nb114_iyH2,             112
.equiv          nb114_izH2,             128
.equiv          nb114_ixM,              144
.equiv          nb114_iyM,              160
.equiv          nb114_izM,              176
.equiv          nb114_jxO,              192
.equiv          nb114_jyO,              208
.equiv          nb114_jzO,              224
.equiv          nb114_jxH1,             240
.equiv          nb114_jyH1,             256
.equiv          nb114_jzH1,             272
.equiv          nb114_jxH2,             288
.equiv          nb114_jyH2,             304
.equiv          nb114_jzH2,             320
.equiv          nb114_jxM,              336
.equiv          nb114_jyM,              352
.equiv          nb114_jzM,              368
.equiv          nb114_dxOO,             384
.equiv          nb114_dyOO,             400
.equiv          nb114_dzOO,             416
.equiv          nb114_dxH1H1,           432
.equiv          nb114_dyH1H1,           448
.equiv          nb114_dzH1H1,           464
.equiv          nb114_dxH1H2,           480
.equiv          nb114_dyH1H2,           496
.equiv          nb114_dzH1H2,           512
.equiv          nb114_dxH1M,            528
.equiv          nb114_dyH1M,            544
.equiv          nb114_dzH1M,            560
.equiv          nb114_dxH2H1,           576
.equiv          nb114_dyH2H1,           592
.equiv          nb114_dzH2H1,           608
.equiv          nb114_dxH2H2,           624
.equiv          nb114_dyH2H2,           640
.equiv          nb114_dzH2H2,           656
.equiv          nb114_dxH2M,            672
.equiv          nb114_dyH2M,            688
.equiv          nb114_dzH2M,            704
.equiv          nb114_dxMH1,            720
.equiv          nb114_dyMH1,            736
.equiv          nb114_dzMH1,            752
.equiv          nb114_dxMH2,            768
.equiv          nb114_dyMH2,            784
.equiv          nb114_dzMH2,            800
.equiv          nb114_dxMM,             816
.equiv          nb114_dyMM,             832
.equiv          nb114_dzMM,             848
.equiv          nb114_qqMM,             864
.equiv          nb114_qqMH,             880
.equiv          nb114_qqHH,             896
.equiv          nb114_two,              912
.equiv          nb114_c6,               944
.equiv          nb114_c12,              960
.equiv          nb114_vctot,            976
.equiv          nb114_Vvdwtot,          992
.equiv          nb114_fixO,             1008
.equiv          nb114_fiyO,             1024
.equiv          nb114_fizO,             1040
.equiv          nb114_fixH1,            1056
.equiv          nb114_fiyH1,            1072
.equiv          nb114_fizH1,            1088
.equiv          nb114_fixH2,            1104
.equiv          nb114_fiyH2,            1120
.equiv          nb114_fizH2,            1136
.equiv          nb114_fixM,             1152
.equiv          nb114_fiyM,             1168
.equiv          nb114_fizM,             1184
.equiv          nb114_fjxO,             1200
.equiv          nb114_fjyO,             1216
.equiv          nb114_fjzO,             1232
.equiv          nb114_fjxH1,            1248
.equiv          nb114_fjyH1,            1264
.equiv          nb114_fjzH1,            1280
.equiv          nb114_fjxH2,            1296
.equiv          nb114_fjyH2,            1312
.equiv          nb114_fjzH2,            1328
.equiv          nb114_fjxM,             1344
.equiv          nb114_fjyM,             1360
.equiv          nb114_fjzM,             1376
.equiv          nb114_half,             1392
.equiv          nb114_three,            1408
.equiv          nb114_six,              1424
.equiv          nb114_twelve,           1440
.equiv          nb114_rsqOO,            1456
.equiv          nb114_rsqH1H1,          1472
.equiv          nb114_rsqH1H2,          1488
.equiv          nb114_rsqH1M,           1504
.equiv          nb114_rsqH2H1,          1520
.equiv          nb114_rsqH2H2,          1536
.equiv          nb114_rsqH2M,           1552
.equiv          nb114_rsqMH1,           1568
.equiv          nb114_rsqMH2,           1584
.equiv          nb114_rsqMM,            1600
.equiv          nb114_rinvsqOO,         1616
.equiv          nb114_rinvH1H1,         1632
.equiv          nb114_rinvH1H2,         1648
.equiv          nb114_rinvH1M,          1664
.equiv          nb114_rinvH2H1,         1680
.equiv          nb114_rinvH2H2,         1696
.equiv          nb114_rinvH2M,          1712
.equiv          nb114_rinvMH1,          1728
.equiv          nb114_rinvMH2,          1744
.equiv          nb114_rinvMM,           1760
.equiv          nb114_is3,              1776
.equiv          nb114_ii3,              1780
.equiv          nb114_innerjjnr,        1784
.equiv          nb114_innerk,           1788
.equiv          nb114_n,                1792
.equiv          nb114_nn1,              1796
.equiv          nb114_nri,              1800
.equiv          nb114_nouter,           1804
.equiv          nb114_ninner,           1808
.equiv          nb114_salign,           1812
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1816		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb114_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb114_p_nri]
	mov ecx, [ecx]
	mov [esp + nb114_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb114_nouter], eax
	mov [esp + nb114_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb114_half], eax
	mov [esp + nb114_half+4], ebx
	movsd xmm1, [esp + nb114_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# 6.0
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# 12.0
	movapd [esp + nb114_half], xmm1
	movapd [esp + nb114_two], xmm2
	movapd [esp + nb114_three], xmm3
	movapd [esp + nb114_six], xmm4
	movapd [esp + nb114_twelve], xmm5
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb114_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb114_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb114_p_facel]
	movsd xmm6, [esi]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb114_qqMM], xmm3
	movapd [esp + nb114_qqMH], xmm4
	movapd [esp + nb114_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb114_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb114_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb114_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movhpd xmm0, [eax + edx*8 + 8]
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb114_c6], xmm0
	movapd [esp + nb114_c12], xmm1

.nb114_threadloop:
        mov   esi, [ebp + nb114_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb114_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb114_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb114_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb114_n], eax
        mov [esp + nb114_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb114_outerstart
        jmp .nb114_end

.nb114_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb114_nouter]
	mov [esp + nb114_nouter], ebx

.nb114_outer:
	mov   eax, [ebp + nb114_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb114_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb114_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb114_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb114_pos]    ;# eax = base of pos[]  
	mov   [esp + nb114_ii3], ebx		

	addsd xmm3, [eax + ebx*8] 	;# ox
	addsd xmm4, [eax + ebx*8 + 8] 	;# oy
	addsd xmm5, [eax + ebx*8 + 16]	;# oz	
	addsd xmm6, [eax + ebx*8 + 24] 	;# h1x
	addsd xmm7, [eax + ebx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [esp + nb114_ixO], xmm3
	movapd [esp + nb114_iyO], xmm4
	movapd [esp + nb114_izO], xmm5
	movapd [esp + nb114_ixH1], xmm6
	movapd [esp + nb114_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [eax + ebx*8 + 40] ;# h1z
	addsd xmm0, [eax + ebx*8 + 48] ;# h2x
	addsd xmm1, [eax + ebx*8 + 56] ;# h2y
	addsd xmm2, [eax + ebx*8 + 64] ;# h2z
	addsd xmm3, [eax + ebx*8 + 72] ;# mx
	addsd xmm4, [eax + ebx*8 + 80] ;# my
	addsd xmm5, [eax + ebx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb114_izH1], xmm6
	movapd [esp + nb114_ixH2], xmm0
	movapd [esp + nb114_iyH2], xmm1
	movapd [esp + nb114_izH2], xmm2
	movapd [esp + nb114_ixM], xmm3
	movapd [esp + nb114_iyM], xmm4
	movapd [esp + nb114_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb114_vctot], xmm4
	movapd [esp + nb114_Vvdwtot], xmm4
	movapd [esp + nb114_fixO], xmm4
	movapd [esp + nb114_fiyO], xmm4
	movapd [esp + nb114_fizO], xmm4
	movapd [esp + nb114_fixH1], xmm4
	movapd [esp + nb114_fiyH1], xmm4
	movapd [esp + nb114_fizH1], xmm4
	movapd [esp + nb114_fixH2], xmm4
	movapd [esp + nb114_fiyH2], xmm4
	movapd [esp + nb114_fizH2], xmm4
	movapd [esp + nb114_fixM], xmm4
	movapd [esp + nb114_fiyM], xmm4
	movapd [esp + nb114_fizM], xmm4
	
	mov   eax, [ebp + nb114_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb114_pos] 
	mov   edi, [ebp + nb114_faction]	
	mov   eax, [ebp + nb114_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb114_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb114_ninner]
	mov   [esp + nb114_ninner], ecx
	add   edx, 0
	mov   [esp + nb114_innerk], edx    ;# number of innerloop atoms 
	jge   .nb114_unroll_loop
	jmp   .nb114_checksingle
.nb114_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb114_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb114_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb114_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8]
	movlpd xmm2, [esi + ebx*8]
	movhpd xmm0, [esi + eax*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 8]
	movlpd xmm3, [esi + eax*8 + 16]
	movlpd xmm5, [esi + ebx*8 + 16]
	movhpd xmm3, [esi + eax*8 + 24]
	movhpd xmm5, [esi + ebx*8 + 24]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# ox 
	unpckhpd xmm1, xmm2 ;# oy
	unpcklpd xmm3, xmm5 ;# ox 
	unpckhpd xmm4, xmm5 ;# oy
	movapd 	[esp + nb114_jxO], xmm0
	movapd 	[esp + nb114_jyO], xmm1
	movapd 	[esp + nb114_jzO], xmm3
	movapd 	[esp + nb114_jxH1], xmm4
	
	;# load h1y, h1z, h2x, h2y 
	movlpd xmm0, [esi + eax*8 + 32]
	movlpd xmm2, [esi + ebx*8 + 32]
	movhpd xmm0, [esi + eax*8 + 40]
	movhpd xmm2, [esi + ebx*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movlpd xmm5, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movhpd xmm5, [esi + ebx*8 + 56]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h1y
	unpckhpd xmm1, xmm2 ;# h1z
	unpcklpd xmm3, xmm5 ;# h2x
	unpckhpd xmm4, xmm5 ;# h2y
	movapd 	[esp + nb114_jyH1], xmm0
	movapd 	[esp + nb114_jzH1], xmm1
	movapd 	[esp + nb114_jxH2], xmm3
	movapd 	[esp + nb114_jyH2], xmm4
	
	;# load h2z, mx, my, mz
	movlpd xmm0, [esi + eax*8 + 64]
	movlpd xmm2, [esi + ebx*8 + 64]
	movhpd xmm0, [esi + eax*8 + 72]
	movhpd xmm2, [esi + ebx*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm5, [esi + ebx*8 + 80]
	movhpd xmm3, [esi + eax*8 + 88]
	movhpd xmm5, [esi + ebx*8 + 88]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h2z
	unpckhpd xmm1, xmm2 ;# mx
	unpcklpd xmm3, xmm5 ;# my
	unpckhpd xmm4, xmm5 ;# mz
	movapd 	[esp + nb114_jzH2], xmm0
	movapd 	[esp + nb114_jxM], xmm1
	movapd 	[esp + nb114_jyM], xmm3
	movapd 	[esp + nb114_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [esp + nb114_ixO]
	movapd xmm1, [esp + nb114_iyO]
	movapd xmm2, [esp + nb114_izO]
	movapd xmm3, [esp + nb114_ixH1]
	movapd xmm4, [esp + nb114_iyH1]
	movapd xmm5, [esp + nb114_izH1]
	subpd  xmm0, [esp + nb114_jxO]
	subpd  xmm1, [esp + nb114_jyO]
	subpd  xmm2, [esp + nb114_jzO]
	subpd  xmm3, [esp + nb114_jxH1]
	subpd  xmm4, [esp + nb114_jyH1]
	subpd  xmm5, [esp + nb114_jzH1]
	movapd [esp + nb114_dxOO], xmm0
	movapd [esp + nb114_dyOO], xmm1
	movapd [esp + nb114_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb114_dxH1H1], xmm3
	movapd [esp + nb114_dyH1H1], xmm4
	movapd [esp + nb114_dzH1H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb114_rsqOO], xmm0
	movapd [esp + nb114_rsqH1H1], xmm3

	movapd xmm0, [esp + nb114_ixH1]
	movapd xmm1, [esp + nb114_iyH1]
	movapd xmm2, [esp + nb114_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb114_jxH2]
	subpd  xmm1, [esp + nb114_jyH2]
	subpd  xmm2, [esp + nb114_jzH2]
	subpd  xmm3, [esp + nb114_jxM]
	subpd  xmm4, [esp + nb114_jyM]
	subpd  xmm5, [esp + nb114_jzM]
	movapd [esp + nb114_dxH1H2], xmm0
	movapd [esp + nb114_dyH1H2], xmm1
	movapd [esp + nb114_dzH1H2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb114_dxH1M], xmm3
	movapd [esp + nb114_dyH1M], xmm4
	movapd [esp + nb114_dzH1M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb114_rsqH1H2], xmm0
	movapd [esp + nb114_rsqH1M], xmm3

	movapd xmm0, [esp + nb114_ixH2]
	movapd xmm1, [esp + nb114_iyH2]
	movapd xmm2, [esp + nb114_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb114_jxH1]
	subpd  xmm1, [esp + nb114_jyH1]
	subpd  xmm2, [esp + nb114_jzH1]
	subpd  xmm3, [esp + nb114_jxH2]
	subpd  xmm4, [esp + nb114_jyH2]
	subpd  xmm5, [esp + nb114_jzH2]
	movapd [esp + nb114_dxH2H1], xmm0
	movapd [esp + nb114_dyH2H1], xmm1
	movapd [esp + nb114_dzH2H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb114_dxH2H2], xmm3
	movapd [esp + nb114_dyH2H2], xmm4
	movapd [esp + nb114_dzH2H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb114_rsqH2H1], xmm0
	movapd [esp + nb114_rsqH2H2], xmm3

	movapd xmm0, [esp + nb114_ixH2]
	movapd xmm1, [esp + nb114_iyH2]
	movapd xmm2, [esp + nb114_izH2]
	movapd xmm3, [esp + nb114_ixM]
	movapd xmm4, [esp + nb114_iyM]
	movapd xmm5, [esp + nb114_izM]
	subpd  xmm0, [esp + nb114_jxM]
	subpd  xmm1, [esp + nb114_jyM]
	subpd  xmm2, [esp + nb114_jzM]
	subpd  xmm3, [esp + nb114_jxH1]
	subpd  xmm4, [esp + nb114_jyH1]
	subpd  xmm5, [esp + nb114_jzH1]
	movapd [esp + nb114_dxH2M], xmm0
	movapd [esp + nb114_dyH2M], xmm1
	movapd [esp + nb114_dzH2M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb114_dxMH1], xmm3
	movapd [esp + nb114_dyMH1], xmm4
	movapd [esp + nb114_dzMH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb114_rsqH2M], xmm0
	movapd [esp + nb114_rsqMH1], xmm4

	movapd xmm0, [esp + nb114_ixM]
	movapd xmm1, [esp + nb114_iyM]
	movapd xmm2, [esp + nb114_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb114_jxH2]
	subpd  xmm1, [esp + nb114_jyH2]
	subpd  xmm2, [esp + nb114_jzH2]
	subpd  xmm3, [esp + nb114_jxM]
	subpd  xmm4, [esp + nb114_jyM]
	subpd  xmm5, [esp + nb114_jzM]
	movapd [esp + nb114_dxMH2], xmm0
	movapd [esp + nb114_dyMH2], xmm1
	movapd [esp + nb114_dzMH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb114_dxMM], xmm3
	movapd [esp + nb114_dyMM], xmm4
	movapd [esp + nb114_dzMM], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb114_rsqMH2], xmm0
	movapd [esp + nb114_rsqMM], xmm4
	
	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtpd2ps xmm1, xmm0
	cvtpd2ps xmm5, xmm4
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1   ;# luA
	cvtps2pd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114_half] ;# iter1 
	mulpd   xmm7, [esp + nb114_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114_half] ;# rinv 
	mulpd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvMH2], xmm1
	movapd [esp + nb114_rinvMM], xmm5

	movapd xmm0, [esp + nb114_rsqOO]
	movapd xmm4, [esp + nb114_rsqH1H1]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb114_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114_half] ;# rinv 
	mulpd   xmm5, [esp + nb114_half] ;# rinv
	mulpd   xmm1, xmm1
	movapd [esp + nb114_rinvsqOO], xmm1
	movapd [esp + nb114_rinvH1H1], xmm5

	movapd xmm0, [esp + nb114_rsqH1H2]
	movapd xmm4, [esp + nb114_rsqH1M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114_half] ;# iter1 
	mulpd   xmm7, [esp + nb114_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114_half] ;# rinv 
	mulpd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvH1H2], xmm1
	movapd [esp + nb114_rinvH1M], xmm5

	movapd xmm0, [esp + nb114_rsqH2H1]
	movapd xmm4, [esp + nb114_rsqH2H2]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114_half] ;# iter1a 
	mulpd   xmm7, [esp + nb114_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114_half] ;# rinv 
	mulpd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvH2H1], xmm1
	movapd [esp + nb114_rinvH2H2], xmm5

	movapd xmm0, [esp + nb114_rsqMH1]
	movapd xmm4, [esp + nb114_rsqH2M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114_half] ;# iter1a 
	mulpd   xmm7, [esp + nb114_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114_half] ;# rinv 
	mulpd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvMH1], xmm1
	movapd [esp + nb114_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb114_rinvsqOO] ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulpd   xmm1, xmm1 ;# rinv4
	mulpd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulpd	xmm2, xmm2 ;# rinvtwelve
	mulpd  xmm1, [esp + nb114_c6]
	mulpd  xmm2, [esp + nb114_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb114_Vvdwtot]
	mulpd  xmm1, [esp + nb114_six]
	mulpd  xmm2, [esp + nb114_twelve]
	subpd  xmm2, xmm1
	mulpd  xmm2, xmm0
	movapd [esp + nb114_Vvdwtot], xmm3

	movapd xmm0, xmm2
	movapd xmm1, xmm2 ;# fscal

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb114_dxOO]
	mulpd xmm1, [esp + nb114_dyOO]
	mulpd xmm2, [esp + nb114_dzOO]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixO]
	addpd xmm1, [esp + nb114_fiyO]
	addpd xmm2, [esp + nb114_fizO]
	movapd [esp + nb114_fjxO], xmm3
	movapd [esp + nb114_fjyO], xmm4
	movapd [esp + nb114_fjzO], xmm5
	movapd [esp + nb114_fixO], xmm0
	movapd [esp + nb114_fiyO], xmm1
	movapd [esp + nb114_fizO], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb114_rinvH1H1]
	movapd xmm7, xmm0		;# xmm7=rinv 
	mulpd  xmm0, xmm0		;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqHH] ;# xmm7=vcoul
	mulpd  xmm0, xmm7		;# fscal
	addpd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movapd [esp + nb114_vctot], xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb114_dxH1H1]
	mulpd xmm1, [esp + nb114_dyH1H1]
	mulpd xmm2, [esp + nb114_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixH1]
	addpd xmm1, [esp + nb114_fiyH1]
	addpd xmm2, [esp + nb114_fizH1]
	movapd [esp + nb114_fjxH1], xmm3
	movapd [esp + nb114_fjyH1], xmm4
	movapd [esp + nb114_fjzH1], xmm5
	movapd [esp + nb114_fixH1], xmm0
	movapd [esp + nb114_fiyH1], xmm1
	movapd [esp + nb114_fizH1], xmm2

	;# H1-H2 interaction  
	movapd xmm0, [esp + nb114_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqHH] ;# vcoul
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movapd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb114_dxH1H2]
	mulpd xmm1, [esp + nb114_dyH1H2]
	mulpd xmm2, [esp + nb114_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixH1]
	addpd xmm1, [esp + nb114_fiyH1]
	addpd xmm2, [esp + nb114_fizH1]
	movapd [esp + nb114_fjxH2], xmm3
	movapd [esp + nb114_fjyH2], xmm4
	movapd [esp + nb114_fjzH2], xmm5
	movapd [esp + nb114_fixH1], xmm0
	movapd [esp + nb114_fiyH1], xmm1
	movapd [esp + nb114_fizH1], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb114_rinvH1M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqMH] ;# vcoul
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movapd [esp + nb114_vctot], xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb114_dxH1M]
	mulpd xmm1, [esp + nb114_dyH1M]
	mulpd xmm2, [esp + nb114_dzH1M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixH1]
	addpd xmm1, [esp + nb114_fiyH1]
	addpd xmm2, [esp + nb114_fizH1]
	movapd [esp + nb114_fjxM], xmm3
	movapd [esp + nb114_fjyM], xmm4
	movapd [esp + nb114_fjzM], xmm5
	movapd [esp + nb114_fixH1], xmm0
	movapd [esp + nb114_fiyH1], xmm1
	movapd [esp + nb114_fizH1], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb114_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqHH] 
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movapd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb114_fjxH1]
	movapd xmm4, [esp + nb114_fjyH1]
	movapd xmm5, [esp + nb114_fjzH1]
	mulpd xmm0, [esp + nb114_dxH2H1]
	mulpd xmm1, [esp + nb114_dyH2H1]
	mulpd xmm2, [esp + nb114_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixH2]
	addpd xmm1, [esp + nb114_fiyH2]
	addpd xmm2, [esp + nb114_fizH2]
	movapd [esp + nb114_fjxH1], xmm3
	movapd [esp + nb114_fjyH1], xmm4
	movapd [esp + nb114_fjzH1], xmm5
	movapd [esp + nb114_fixH2], xmm0
	movapd [esp + nb114_fiyH2], xmm1
	movapd [esp + nb114_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb114_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqHH] ;# vcoul
	
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movapd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb114_fjxH2]
	movapd xmm4, [esp + nb114_fjyH2]
	movapd xmm5, [esp + nb114_fjzH2]
	mulpd xmm0, [esp + nb114_dxH2H2]
	mulpd xmm1, [esp + nb114_dyH2H2]
	mulpd xmm2, [esp + nb114_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixH2]
	addpd xmm1, [esp + nb114_fiyH2]
	addpd xmm2, [esp + nb114_fizH2]
	movapd [esp + nb114_fjxH2], xmm3
	movapd [esp + nb114_fjyH2], xmm4
	movapd [esp + nb114_fjzH2], xmm5
	movapd [esp + nb114_fixH2], xmm0
	movapd [esp + nb114_fiyH2], xmm1
	movapd [esp + nb114_fizH2], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb114_rinvH2M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqMH] ;# vcoul
	
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot] 
	movapd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb114_fjxM]
	movapd xmm4, [esp + nb114_fjyM]
	movapd xmm5, [esp + nb114_fjzM]
	mulpd xmm0, [esp + nb114_dxH2M]
	mulpd xmm1, [esp + nb114_dyH2M]
	mulpd xmm2, [esp + nb114_dzH2M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixH2]
	addpd xmm1, [esp + nb114_fiyH2]
	addpd xmm2, [esp + nb114_fizH2]
	movapd [esp + nb114_fjxM], xmm3
	movapd [esp + nb114_fjyM], xmm4
	movapd [esp + nb114_fjzM], xmm5
	movapd [esp + nb114_fixH2], xmm0
	movapd [esp + nb114_fiyH2], xmm1
	movapd [esp + nb114_fizH2], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb114_rinvMH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqMH]
	
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movapd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb114_fjxH1]
	movapd xmm4, [esp + nb114_fjyH1]
	movapd xmm5, [esp + nb114_fjzH1]
	mulpd xmm0, [esp + nb114_dxMH1]
	mulpd xmm1, [esp + nb114_dyMH1]
	mulpd xmm2, [esp + nb114_dzMH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixM]
	addpd xmm1, [esp + nb114_fiyM]
	addpd xmm2, [esp + nb114_fizM]
	movapd [esp + nb114_fjxH1], xmm3
	movapd [esp + nb114_fjyH1], xmm4
	movapd [esp + nb114_fjzH1], xmm5
	movapd [esp + nb114_fixM], xmm0
	movapd [esp + nb114_fiyM], xmm1
	movapd [esp + nb114_fizM], xmm2

	;# M-H2 interaction 
	movapd xmm0, [esp + nb114_rinvMH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqMH] 
	
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot]
	movapd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb114_fjxH2]
	movapd xmm4, [esp + nb114_fjyH2]
	movapd xmm5, [esp + nb114_fjzH2]
	mulpd xmm0, [esp + nb114_dxMH2]
	mulpd xmm1, [esp + nb114_dyMH2]
	mulpd xmm2, [esp + nb114_dzMH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixM]
	addpd xmm1, [esp + nb114_fiyM]
	addpd xmm2, [esp + nb114_fizM]
	movapd [esp + nb114_fjxH2], xmm3
	movapd [esp + nb114_fjyH2], xmm4
	movapd [esp + nb114_fjzH2], xmm5
	movapd [esp + nb114_fixM], xmm0
	movapd [esp + nb114_fiyM], xmm1
	movapd [esp + nb114_fizM], xmm2

	;# M-M interaction 
	movapd xmm0, [esp + nb114_rinvMM]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	mulpd  xmm7, [esp + nb114_qqMM] ;# vcoul
	
	mulpd  xmm0, xmm7
	addpd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movapd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb114_fjxM]
	movapd xmm4, [esp + nb114_fjyM] 
	movapd xmm5, [esp + nb114_fjzM]
	mulpd xmm0, [esp + nb114_dxMM]
	mulpd xmm1, [esp + nb114_dyMM]
	mulpd xmm2, [esp + nb114_dzMM]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb114_fixM]
	addpd xmm1, [esp + nb114_fiyM]
	addpd xmm2, [esp + nb114_fizM]
	movapd [esp + nb114_fjxM], xmm3
	movapd [esp + nb114_fjyM], xmm4
	movapd [esp + nb114_fjzM], xmm5
	movapd [esp + nb114_fixM], xmm0
	movapd [esp + nb114_fiyM], xmm1
	movapd [esp + nb114_fizM], xmm2

	mov edi, [ebp + nb114_faction]
	
	;# Did all interactions - now update j forces 
	;# Step1 - transpose fjxO, fjyO and fjzO, fjxH1
	movapd xmm0, [esp + nb114_fjxO]
	movapd xmm1, [esp + nb114_fjyO]
	movapd xmm2, [esp + nb114_fjzO]
	movapd xmm3, [esp + nb114_fjxH1]
	movapd xmm4, xmm0
	movapd xmm5, xmm2
	unpcklpd xmm0, xmm1   ;# fjOxA fjOyA
	unpckhpd xmm4, xmm1   ;# fjOxB fjOyB
	unpcklpd xmm2, xmm3   ;# fjOzA fjH1xA
	unpckhpd xmm5, xmm3   ;# fjOzB fjH1xB
	movlpd xmm1, [edi + eax*8]
	movlpd xmm3, [edi + ebx*8]
	movhpd xmm1, [edi + eax*8 + 8]
	movhpd xmm3, [edi + ebx*8 + 8]
	movlpd xmm6, [edi + eax*8 + 16]
	movlpd xmm7, [edi + ebx*8 + 16]
	movhpd xmm6, [edi + eax*8 + 24]
	movhpd xmm7, [edi + ebx*8 + 24]
	addpd  xmm1, xmm0
	addpd  xmm3, xmm4
	addpd  xmm6, xmm2
	addpd  xmm7, xmm5
	movlpd [edi + eax*8],      xmm1
	movlpd [edi + ebx*8],      xmm3
	movhpd [edi + eax*8 + 8],  xmm1
	movhpd [edi + ebx*8 + 8],  xmm3
	movlpd [edi + eax*8 + 16], xmm6
	movlpd [edi + ebx*8 + 16], xmm7
	movhpd [edi + eax*8 + 24], xmm6
	movhpd [edi + ebx*8 + 24], xmm7

	;# Step2 - transpose fjyH1, fjzH1 and fjxH2, fjyH2
	movapd xmm0, [esp + nb114_fjyH1]
	movapd xmm1, [esp + nb114_fjzH1]
	movapd xmm2, [esp + nb114_fjxH2]
	movapd xmm3, [esp + nb114_fjyH2]
	movapd xmm4, xmm0
	movapd xmm5, xmm2
	unpcklpd xmm0, xmm1   ;# fjOxA fjOyA
	unpckhpd xmm4, xmm1   ;# fjOxB fjOyB
	unpcklpd xmm2, xmm3   ;# fjOzA fjH1xA
	unpckhpd xmm5, xmm3   ;# fjOzB fjH1xB
	movlpd xmm1, [edi + eax*8 + 32]
	movlpd xmm3, [edi + ebx*8 + 32]
	movhpd xmm1, [edi + eax*8 + 40]
	movhpd xmm3, [edi + ebx*8 + 40]
	movlpd xmm6, [edi + eax*8 + 48]
	movlpd xmm7, [edi + ebx*8 + 48]
	movhpd xmm6, [edi + eax*8 + 56]
	movhpd xmm7, [edi + ebx*8 + 56]
	addpd  xmm1, xmm0
	addpd  xmm3, xmm4
	addpd  xmm6, xmm2
	addpd  xmm7, xmm5
	movlpd [edi + eax*8 + 32], xmm1
	movlpd [edi + ebx*8 + 32], xmm3
	movhpd [edi + eax*8 + 40], xmm1
	movhpd [edi + ebx*8 + 40], xmm3
	movlpd [edi + eax*8 + 48], xmm6
	movlpd [edi + ebx*8 + 48], xmm7
	movhpd [edi + eax*8 + 56], xmm6
	movhpd [edi + ebx*8 + 56], xmm7

	;# Step3 - transpose fjzH2, fjxM and fjyM, fjzM
	movapd xmm0, [esp + nb114_fjzH2]
	movapd xmm1, [esp + nb114_fjxM]
	movapd xmm2, [esp + nb114_fjyM]
	movapd xmm3, [esp + nb114_fjzM]
	movapd xmm4, xmm0
	movapd xmm5, xmm2
	unpcklpd xmm0, xmm1   ;# fjOxA fjOyA
	unpckhpd xmm4, xmm1   ;# fjOxB fjOyB
	unpcklpd xmm2, xmm3   ;# fjOzA fjH1xA
	unpckhpd xmm5, xmm3   ;# fjOzB fjH1xB
	movlpd xmm1, [edi + eax*8 + 64]
	movlpd xmm3, [edi + ebx*8 + 64]
	movhpd xmm1, [edi + eax*8 + 72]
	movhpd xmm3, [edi + ebx*8 + 72]
	movlpd xmm6, [edi + eax*8 + 80]
	movlpd xmm7, [edi + ebx*8 + 80]
	movhpd xmm6, [edi + eax*8 + 88]
	movhpd xmm7, [edi + ebx*8 + 88]
	addpd  xmm1, xmm0
	addpd  xmm3, xmm4
	addpd  xmm6, xmm2
	addpd  xmm7, xmm5
	movlpd [edi + eax*8 + 64], xmm1
	movlpd [edi + ebx*8 + 64], xmm3
	movhpd [edi + eax*8 + 72], xmm1
	movhpd [edi + ebx*8 + 72], xmm3
	movlpd [edi + eax*8 + 80], xmm6
	movlpd [edi + ebx*8 + 80], xmm7
	movhpd [edi + eax*8 + 88], xmm6
	movhpd [edi + ebx*8 + 88], xmm7

	;# should we do one more iteration? 
	sub dword ptr [esp + nb114_innerk],  2
	jl    .nb114_checksingle
	jmp   .nb114_unroll_loop
.nb114_checksingle:
	mov   edx, [esp + nb114_innerk]
	and   edx, 1
	jnz   .nb114_dosingle
	jmp   .nb114_updateouterdata
.nb114_dosingle:
	mov   edx, [esp + nb114_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb114_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8]
	movhpd xmm0, [esi + eax*8 + 8]
	movlpd xmm1, [esi + eax*8 + 16]
	movhpd xmm1, [esi + eax*8 + 24]
	movlpd xmm2, [esi + eax*8 + 32]
	movhpd xmm2, [esi + eax*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm4, [esi + eax*8 + 72]
	movlpd xmm5, [esi + eax*8 + 80]
	movhpd xmm5, [esi + eax*8 + 88]
	movsd  [esp + nb114_jxO], xmm0
	movsd  [esp + nb114_jzO], xmm1
	movsd  [esp + nb114_jyH1], xmm2
	movsd  [esp + nb114_jxH2], xmm3
	movsd  [esp + nb114_jzH2], xmm4
	movsd  [esp + nb114_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [esp + nb114_jyO], xmm0
	movsd  [esp + nb114_jxH1], xmm1
	movsd  [esp + nb114_jzH1], xmm2
	movsd  [esp + nb114_jyH2], xmm3
	movsd  [esp + nb114_jxM], xmm4
	movsd  [esp + nb114_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [esp + nb114_ixO]
	movapd xmm1, [esp + nb114_iyO]
	movapd xmm2, [esp + nb114_izO]
	movapd xmm3, [esp + nb114_ixH1]
	movapd xmm4, [esp + nb114_iyH1]
	movapd xmm5, [esp + nb114_izH1]
	subsd  xmm0, [esp + nb114_jxO]
	subsd  xmm1, [esp + nb114_jyO]
	subsd  xmm2, [esp + nb114_jzO]
	subsd  xmm3, [esp + nb114_jxH1]
	subsd  xmm4, [esp + nb114_jyH1]
	subsd  xmm5, [esp + nb114_jzH1]
	movapd [esp + nb114_dxOO], xmm0
	movapd [esp + nb114_dyOO], xmm1
	movapd [esp + nb114_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb114_dxH1H1], xmm3
	movapd [esp + nb114_dyH1H1], xmm4
	movapd [esp + nb114_dzH1H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb114_rsqOO], xmm0
	movapd [esp + nb114_rsqH1H1], xmm3

	movapd xmm0, [esp + nb114_ixH1]
	movapd xmm1, [esp + nb114_iyH1]
	movapd xmm2, [esp + nb114_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb114_jxH2]
	subsd  xmm1, [esp + nb114_jyH2]
	subsd  xmm2, [esp + nb114_jzH2]
	subsd  xmm3, [esp + nb114_jxM]
	subsd  xmm4, [esp + nb114_jyM]
	subsd  xmm5, [esp + nb114_jzM]
	movapd [esp + nb114_dxH1H2], xmm0
	movapd [esp + nb114_dyH1H2], xmm1
	movapd [esp + nb114_dzH1H2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb114_dxH1M], xmm3
	movapd [esp + nb114_dyH1M], xmm4
	movapd [esp + nb114_dzH1M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb114_rsqH1H2], xmm0
	movapd [esp + nb114_rsqH1M], xmm3

	movapd xmm0, [esp + nb114_ixH2]
	movapd xmm1, [esp + nb114_iyH2]
	movapd xmm2, [esp + nb114_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb114_jxH1]
	subsd  xmm1, [esp + nb114_jyH1]
	subsd  xmm2, [esp + nb114_jzH1]
	subsd  xmm3, [esp + nb114_jxH2]
	subsd  xmm4, [esp + nb114_jyH2]
	subsd  xmm5, [esp + nb114_jzH2]
	movapd [esp + nb114_dxH2H1], xmm0
	movapd [esp + nb114_dyH2H1], xmm1
	movapd [esp + nb114_dzH2H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb114_dxH2H2], xmm3
	movapd [esp + nb114_dyH2H2], xmm4
	movapd [esp + nb114_dzH2H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb114_rsqH2H1], xmm0
	movapd [esp + nb114_rsqH2H2], xmm3

	movapd xmm0, [esp + nb114_ixH2]
	movapd xmm1, [esp + nb114_iyH2]
	movapd xmm2, [esp + nb114_izH2]
	movapd xmm3, [esp + nb114_ixM]
	movapd xmm4, [esp + nb114_iyM]
	movapd xmm5, [esp + nb114_izM]
	subsd  xmm0, [esp + nb114_jxM]
	subsd  xmm1, [esp + nb114_jyM]
	subsd  xmm2, [esp + nb114_jzM]
	subsd  xmm3, [esp + nb114_jxH1]
	subsd  xmm4, [esp + nb114_jyH1]
	subsd  xmm5, [esp + nb114_jzH1]
	movapd [esp + nb114_dxH2M], xmm0
	movapd [esp + nb114_dyH2M], xmm1
	movapd [esp + nb114_dzH2M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb114_dxMH1], xmm3
	movapd [esp + nb114_dyMH1], xmm4
	movapd [esp + nb114_dzMH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb114_rsqH2M], xmm0
	movapd [esp + nb114_rsqMH1], xmm4

	movapd xmm0, [esp + nb114_ixM]
	movapd xmm1, [esp + nb114_iyM]
	movapd xmm2, [esp + nb114_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb114_jxH2]
	subsd  xmm1, [esp + nb114_jyH2]
	subsd  xmm2, [esp + nb114_jzH2]
	subsd  xmm3, [esp + nb114_jxM]
	subsd  xmm4, [esp + nb114_jyM]
	subsd  xmm5, [esp + nb114_jzM]
	movapd [esp + nb114_dxMH2], xmm0
	movapd [esp + nb114_dyMH2], xmm1
	movapd [esp + nb114_dzMH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb114_dxMM], xmm3
	movapd [esp + nb114_dyMM], xmm4
	movapd [esp + nb114_dzMM], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb114_rsqMH2], xmm0
	movapd [esp + nb114_rsqMM], xmm4

	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtsd2ss xmm1, xmm0
	cvtsd2ss xmm5, xmm4
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1   ;# luA
	cvtss2sd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114_half] ;# iter1 
	mulsd   xmm7, [esp + nb114_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114_half] ;# rinv 
	mulsd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvMH2], xmm1
	movapd [esp + nb114_rinvMM], xmm5

	movapd xmm0, [esp + nb114_rsqOO]
	movapd xmm4, [esp + nb114_rsqH1H1]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb114_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114_half] ;# rinv 
	mulsd   xmm5, [esp + nb114_half] ;# rinv
	mulpd   xmm1, xmm1
	movapd [esp + nb114_rinvsqOO], xmm1
	movapd [esp + nb114_rinvH1H1], xmm5

	movapd xmm0, [esp + nb114_rsqH1H2]
	movapd xmm4, [esp + nb114_rsqH1M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114_half] ;# iter1 
	mulsd   xmm7, [esp + nb114_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114_half] ;# rinv 
	mulsd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvH1H2], xmm1
	movapd [esp + nb114_rinvH1M], xmm5

	movapd xmm0, [esp + nb114_rsqH2H1]
	movapd xmm4, [esp + nb114_rsqH2H2]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114_half] ;# iter1a 
	mulsd   xmm7, [esp + nb114_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114_half] ;# rinv 
	mulsd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvH2H1], xmm1
	movapd [esp + nb114_rinvH2H2], xmm5

	movapd xmm0, [esp + nb114_rsqMH1]
	movapd xmm4, [esp + nb114_rsqH2M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114_half] ;# iter1a 
	mulsd   xmm7, [esp + nb114_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114_half] ;# rinv 
	mulsd   xmm5, [esp + nb114_half] ;# rinv 
	movapd [esp + nb114_rinvMH1], xmm1
	movapd [esp + nb114_rinvH2M], xmm5

	;# start with OO interaction 
	movsd xmm0, [esp + nb114_rinvsqOO] ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movsd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [esp + nb114_c6]
	mulsd  xmm2, [esp + nb114_c12]
	movsd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb114_Vvdwtot]
	mulsd  xmm1, [esp + nb114_six]
	mulsd  xmm2, [esp + nb114_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm2, xmm0
	movsd [esp + nb114_Vvdwtot], xmm3

	movapd xmm0, xmm2
	movapd xmm1, xmm2

	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb114_dxOO]
	mulsd xmm1, [esp + nb114_dyOO]
	mulsd xmm2, [esp + nb114_dzOO]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixO]
	addsd xmm1, [esp + nb114_fiyO]
	addsd xmm2, [esp + nb114_fizO]
	movsd [esp + nb114_fjxO], xmm3
	movsd [esp + nb114_fjyO], xmm4
	movsd [esp + nb114_fjzO], xmm5
	movsd [esp + nb114_fixO], xmm0
	movsd [esp + nb114_fiyO], xmm1
	movsd [esp + nb114_fizO], xmm2

	;# H1-H1 interaction 
	movsd xmm0, [esp + nb114_rinvH1H1]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	mulsd  xmm7, [esp + nb114_qqHH] 
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot]
	movsd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb114_dxH1H1]
	mulsd xmm1, [esp + nb114_dyH1H1]
	mulsd xmm2, [esp + nb114_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixH1]
	addsd xmm1, [esp + nb114_fiyH1]
	addsd xmm2, [esp + nb114_fizH1]
	movsd [esp + nb114_fjxH1], xmm3
	movsd [esp + nb114_fjyH1], xmm4
	movsd [esp + nb114_fjzH1], xmm5
	movsd [esp + nb114_fixH1], xmm0
	movsd [esp + nb114_fiyH1], xmm1
	movsd [esp + nb114_fizH1], xmm2

	;# H1-H2 interaction  
	movsd xmm0, [esp + nb114_rinvH1H2]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	mulsd  xmm7, [esp + nb114_qqHH] ;# vcoul
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot]
	movsd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb114_dxH1H2]
	mulsd xmm1, [esp + nb114_dyH1H2]
	mulsd xmm2, [esp + nb114_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixH1]
	addsd xmm1, [esp + nb114_fiyH1]
	addsd xmm2, [esp + nb114_fizH1]
	movsd [esp + nb114_fjxH2], xmm3
	movsd [esp + nb114_fjyH2], xmm4
	movsd [esp + nb114_fjzH2], xmm5
	movsd [esp + nb114_fixH1], xmm0
	movsd [esp + nb114_fiyH1], xmm1
	movsd [esp + nb114_fizH1], xmm2

	;# H1-M interaction 
	movsd xmm0, [esp + nb114_rinvH1M]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	mulsd  xmm7, [esp + nb114_qqMH] ;# vcoul
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot] 
	movsd [esp + nb114_vctot], xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb114_dxH1M]
	mulsd xmm1, [esp + nb114_dyH1M]
	mulsd xmm2, [esp + nb114_dzH1M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixH1]
	addsd xmm1, [esp + nb114_fiyH1]
	addsd xmm2, [esp + nb114_fizH1]
	movsd [esp + nb114_fjxM], xmm3
	movsd [esp + nb114_fjyM], xmm4
	movsd [esp + nb114_fjzM], xmm5
	movsd [esp + nb114_fixH1], xmm0
	movsd [esp + nb114_fiyH1], xmm1
	movsd [esp + nb114_fizH1], xmm2

	;# H2-H1 interaction 
	movsd xmm0, [esp + nb114_rinvH2H1]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	mulsd  xmm7, [esp + nb114_qqHH] ;# vcoul
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot] 
	movsd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb114_fjxH1]
	movapd xmm4, [esp + nb114_fjyH1]
	movapd xmm5, [esp + nb114_fjzH1]
	mulsd xmm0, [esp + nb114_dxH2H1]
	mulsd xmm1, [esp + nb114_dyH2H1]
	mulsd xmm2, [esp + nb114_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixH2]
	addsd xmm1, [esp + nb114_fiyH2]
	addsd xmm2, [esp + nb114_fizH2]
	movsd [esp + nb114_fjxH1], xmm3
	movsd [esp + nb114_fjyH1], xmm4
	movsd [esp + nb114_fjzH1], xmm5
	movsd [esp + nb114_fixH2], xmm0
	movsd [esp + nb114_fiyH2], xmm1
	movsd [esp + nb114_fizH2], xmm2

	;# H2-H2 interaction 
	movsd xmm0, [esp + nb114_rinvH2H2]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	mulsd  xmm7, [esp + nb114_qqHH] ;# vcoul
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot] 
	movsd [esp + nb114_vctot], xmm7
	
	movsd xmm1, xmm0
	movsd xmm2, xmm0

	movsd xmm3, [esp + nb114_fjxH2]
	movsd xmm4, [esp + nb114_fjyH2]
	movsd xmm5, [esp + nb114_fjzH2]
	mulsd xmm0, [esp + nb114_dxH2H2]
	mulsd xmm1, [esp + nb114_dyH2H2]
	mulsd xmm2, [esp + nb114_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixH2]
	addsd xmm1, [esp + nb114_fiyH2]
	addsd xmm2, [esp + nb114_fizH2]
	movsd [esp + nb114_fjxH2], xmm3
	movsd [esp + nb114_fjyH2], xmm4
	movsd [esp + nb114_fjzH2], xmm5
	movsd [esp + nb114_fixH2], xmm0
	movsd [esp + nb114_fiyH2], xmm1
	movsd [esp + nb114_fizH2], xmm2

	;# H2-M interaction 
	movsd xmm0, [esp + nb114_rinvH2M]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	mulsd  xmm7, [esp + nb114_qqMH] 	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot] 
	movsd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb114_fjxM]
	movapd xmm4, [esp + nb114_fjyM]
	movapd xmm5, [esp + nb114_fjzM]
	mulsd xmm0, [esp + nb114_dxH2M]
	mulsd xmm1, [esp + nb114_dyH2M]
	mulsd xmm2, [esp + nb114_dzH2M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixH2]
	addsd xmm1, [esp + nb114_fiyH2]
	addsd xmm2, [esp + nb114_fizH2]
	movsd [esp + nb114_fjxM], xmm3
	movsd [esp + nb114_fjyM], xmm4
	movsd [esp + nb114_fjzM], xmm5
	movsd [esp + nb114_fixH2], xmm0
	movsd [esp + nb114_fiyH2], xmm1
	movsd [esp + nb114_fizH2], xmm2

	;# M-H1 interaction 
	movsd xmm0, [esp + nb114_rinvMH1]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 

	mulsd  xmm7, [esp + nb114_qqMH] 
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot] ;# local vctot summation variable 
	movsd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb114_fjxH1]
	movapd xmm4, [esp + nb114_fjyH1]
	movapd xmm5, [esp + nb114_fjzH1]
	mulsd xmm0, [esp + nb114_dxMH1]
	mulsd xmm1, [esp + nb114_dyMH1]
	mulsd xmm2, [esp + nb114_dzMH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixM]
	addsd xmm1, [esp + nb114_fiyM]
	addsd xmm2, [esp + nb114_fizM]
	movsd [esp + nb114_fjxH1], xmm3
	movsd [esp + nb114_fjyH1], xmm4
	movsd [esp + nb114_fjzH1], xmm5
	movsd [esp + nb114_fixM], xmm0
	movsd [esp + nb114_fiyM], xmm1
	movsd [esp + nb114_fizM], xmm2

	;# M-H2 interaction 
	movsd xmm0, [esp + nb114_rinvMH2]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 

	mulsd  xmm7, [esp + nb114_qqMH] 
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot]
	movsd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb114_fjxH2]
	movapd xmm4, [esp + nb114_fjyH2]
	movapd xmm5, [esp + nb114_fjzH2]
	mulsd xmm0, [esp + nb114_dxMH2]
	mulsd xmm1, [esp + nb114_dyMH2]
	mulsd xmm2, [esp + nb114_dzMH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixM]
	addsd xmm1, [esp + nb114_fiyM]
	addsd xmm2, [esp + nb114_fizM]
	movsd [esp + nb114_fjxH2], xmm3
	movsd [esp + nb114_fjyH2], xmm4
	movsd [esp + nb114_fjzH2], xmm5
	movsd [esp + nb114_fixM], xmm0
	movsd [esp + nb114_fiyM], xmm1
	movsd [esp + nb114_fizM], xmm2

	;# M-M interaction 
	movsd xmm0, [esp + nb114_rinvMM]
	movsd xmm7, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 

	mulsd  xmm7, [esp + nb114_qqMM] 
	
	mulsd  xmm0, xmm7
	addsd  xmm7, [esp + nb114_vctot] 
	movsd [esp + nb114_vctot], xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb114_fjxM]
	movapd xmm4, [esp + nb114_fjyM]
	movapd xmm5, [esp + nb114_fjzM]
	mulsd xmm0, [esp + nb114_dxMM]
	mulsd xmm1, [esp + nb114_dyMM]
	mulsd xmm2, [esp + nb114_dzMM]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb114_fixM]
	addsd xmm1, [esp + nb114_fiyM]
	addsd xmm2, [esp + nb114_fizM]
	movsd [esp + nb114_fjxM], xmm3
	movsd [esp + nb114_fjyM], xmm4
	movsd [esp + nb114_fjzM], xmm5
	movsd [esp + nb114_fixM], xmm0
	movsd [esp + nb114_fiyM], xmm1
	movsd [esp + nb114_fizM], xmm2

	mov edi, [ebp + nb114_faction]

	;# Did all interactions - now update j forces 
	;# Step1 - merge forces
	movlpd xmm0, [esp + nb114_fjxO]
	movlpd xmm1, [esp + nb114_fjzO]
	movlpd xmm2, [esp + nb114_fjyH1]
	movlpd xmm3, [esp + nb114_fjxH2]
	movlpd xmm4, [esp + nb114_fjzH2]
	movlpd xmm5, [esp + nb114_fjyM]

	movhpd xmm0, [esp + nb114_fjyO]
	movhpd xmm1, [esp + nb114_fjxH1]
	movhpd xmm2, [esp + nb114_fjzH1]
	movhpd xmm3, [esp + nb114_fjyH2]
	movhpd xmm4, [esp + nb114_fjxM]
	movhpd xmm5, [esp + nb114_fjzM]

	movlps xmm6, [edi + eax*8]
	movhps xmm6, [edi + eax*8 + 8]
	movlps xmm7, [edi + eax*8 + 16]
	movhps xmm7, [edi + eax*8 + 24]
	addpd  xmm0, xmm6
	addpd  xmm1, xmm7
	movlps xmm6, [edi + eax*8 + 32]
	movhps xmm6, [edi + eax*8 + 40]
	movlps xmm7, [edi + eax*8 + 48]
	movhps xmm7, [edi + eax*8 + 56]
	addpd  xmm2, xmm6
	addpd  xmm3, xmm7
	movlps xmm6, [edi + eax*8 + 64]
	movhps xmm6, [edi + eax*8 + 72]
	movlps xmm7, [edi + eax*8 + 80]
	movhps xmm7, [edi + eax*8 + 88]
	addpd  xmm4, xmm6
	addpd  xmm5, xmm7
	
	movlps [edi + eax*8],      xmm0
	movhps [edi + eax*8 + 8],  xmm0
	movlps [edi + eax*8 + 16], xmm1
	movhps [edi + eax*8 + 24], xmm1
	movlps [edi + eax*8 + 32], xmm2
	movhps [edi + eax*8 + 40], xmm2
	movlps [edi + eax*8 + 48], xmm3
	movhps [edi + eax*8 + 56], xmm3
	movlps [edi + eax*8 + 64], xmm4
	movhps [edi + eax*8 + 72], xmm4
	movlps [edi + eax*8 + 80], xmm5
	movhps [edi + eax*8 + 88], xmm5
	
.nb114_updateouterdata:
	mov   ecx, [esp + nb114_ii3]
	mov   edi, [ebp + nb114_faction]
	mov   esi, [ebp + nb114_fshift]
	mov   edx, [esp + nb114_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb114_fixO]
	movapd xmm1, [esp + nb114_fiyO]
	movapd xmm2, [esp + nb114_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [edi + ecx*8]
	movsd  xmm4, [edi + ecx*8 + 8]
	movsd  xmm5, [edi + ecx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8],     xmm3
	movsd  [edi + ecx*8 + 8], xmm4
	movsd  [edi + ecx*8 + 16], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb114_fixH1]
	movapd xmm1, [esp + nb114_fiyH1]
	movapd xmm2, [esp + nb114_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 24]
	movsd  xmm4, [edi + ecx*8 + 32]
	movsd  xmm5, [edi + ecx*8 + 40]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 24], xmm3
	movsd  [edi + ecx*8 + 32], xmm4
	movsd  [edi + ecx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb114_fixH2]
	movapd xmm1, [esp + nb114_fiyH2]
	movapd xmm2, [esp + nb114_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 48]
	movsd  xmm4, [edi + ecx*8 + 56]
	movsd  xmm5, [edi + ecx*8 + 64]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 48], xmm3
	movsd  [edi + ecx*8 + 56], xmm4
	movsd  [edi + ecx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb114_fixM]
	movapd xmm1, [esp + nb114_fiyM]
	movapd xmm2, [esp + nb114_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 72]
	movsd  xmm4, [edi + ecx*8 + 80]
	movsd  xmm5, [edi + ecx*8 + 88]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 72], xmm3
	movsd  [edi + ecx*8 + 80], xmm4
	movsd  [edi + ecx*8 + 88], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [esi + edx*8]
	movhpd xmm3, [esi + edx*8 + 8]
	movsd  xmm4, [esi + edx*8 + 16]
	addpd  xmm3, xmm6
	addsd  xmm4, xmm7
	movlpd [esi + edx*8],      xmm3
	movhpd [esi + edx*8 + 8],  xmm3
	movsd  [esi + edx*8 + 16], xmm4

	;# get n from stack
	mov esi, [esp + nb114_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb114_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb114_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb114_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb114_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb114_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb114_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb114_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb114_n], esi
        jmp .nb114_outer
.nb114_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb114_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb114_end
        ;# non-zero, do one more workunit
        jmp   .nb114_threadloop
.nb114_end:
	emms

	mov eax, [esp + nb114_nouter]
	mov ebx, [esp + nb114_ninner]
	mov ecx, [ebp + nb114_outeriter]
	mov edx, [ebp + nb114_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb114_salign]
	add esp, eax
	add esp, 1816
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




	
.globl nb_kernel114nf_ia32_sse2
.globl _nb_kernel114nf_ia32_sse2
nb_kernel114nf_ia32_sse2:	
_nb_kernel114nf_ia32_sse2:	
.equiv          nb114nf_p_nri,          8
.equiv          nb114nf_iinr,           12
.equiv          nb114nf_jindex,         16
.equiv          nb114nf_jjnr,           20
.equiv          nb114nf_shift,          24
.equiv          nb114nf_shiftvec,       28
.equiv          nb114nf_fshift,         32
.equiv          nb114nf_gid,            36
.equiv          nb114nf_pos,            40
.equiv          nb114nf_faction,        44
.equiv          nb114nf_charge,         48
.equiv          nb114nf_p_facel,        52
.equiv          nb114nf_argkrf,         56
.equiv          nb114nf_argcrf,         60
.equiv          nb114nf_Vc,             64
.equiv          nb114nf_type,           68
.equiv          nb114nf_p_ntype,        72
.equiv          nb114nf_vdwparam,       76
.equiv          nb114nf_Vvdw,           80
.equiv          nb114nf_p_tabscale,     84
.equiv          nb114nf_VFtab,          88
.equiv          nb114nf_invsqrta,       92
.equiv          nb114nf_dvda,           96
.equiv          nb114nf_p_gbtabscale,   100
.equiv          nb114nf_GBtab,          104
.equiv          nb114nf_p_nthreads,     108
.equiv          nb114nf_count,          112
.equiv          nb114nf_mtx,            116
.equiv          nb114nf_outeriter,      120
.equiv          nb114nf_inneriter,      124
.equiv          nb114nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb114nf_ixO,            0
.equiv          nb114nf_iyO,            16
.equiv          nb114nf_izO,            32
.equiv          nb114nf_ixH1,           48
.equiv          nb114nf_iyH1,           64
.equiv          nb114nf_izH1,           80
.equiv          nb114nf_ixH2,           96
.equiv          nb114nf_iyH2,           112
.equiv          nb114nf_izH2,           128
.equiv          nb114nf_ixM,            144
.equiv          nb114nf_iyM,            160
.equiv          nb114nf_izM,            176
.equiv          nb114nf_jxO,            192
.equiv          nb114nf_jyO,            208
.equiv          nb114nf_jzO,            224
.equiv          nb114nf_jxH1,           240
.equiv          nb114nf_jyH1,           256
.equiv          nb114nf_jzH1,           272
.equiv          nb114nf_jxH2,           288
.equiv          nb114nf_jyH2,           304
.equiv          nb114nf_jzH2,           320
.equiv          nb114nf_jxM,            336
.equiv          nb114nf_jyM,            352
.equiv          nb114nf_jzM,            368
.equiv          nb114nf_qqMM,           384
.equiv          nb114nf_qqMH,           400
.equiv          nb114nf_qqHH,           416
.equiv          nb114nf_two,            432
.equiv          nb114nf_c6,             448
.equiv          nb114nf_c12,            464
.equiv          nb114nf_vctot,          480
.equiv          nb114nf_Vvdwtot,        496
.equiv          nb114nf_half,           512
.equiv          nb114nf_three,          528
.equiv          nb114nf_rsqOO,          544
.equiv          nb114nf_rsqH1H1,        560
.equiv          nb114nf_rsqH1H2,        576
.equiv          nb114nf_rsqH1M,         592
.equiv          nb114nf_rsqH2H1,        608
.equiv          nb114nf_rsqH2H2,        624
.equiv          nb114nf_rsqH2M,         640
.equiv          nb114nf_rsqMH1,         656
.equiv          nb114nf_rsqMH2,         672
.equiv          nb114nf_rsqMM,          688
.equiv          nb114nf_rinvsqOO,       704
.equiv          nb114nf_rinvH1H1,       720
.equiv          nb114nf_rinvH1H2,       736
.equiv          nb114nf_rinvH1M,        752
.equiv          nb114nf_rinvH2H1,       768
.equiv          nb114nf_rinvH2H2,       784
.equiv          nb114nf_rinvH2M,        800
.equiv          nb114nf_rinvMH1,        816
.equiv          nb114nf_rinvMH2,        832
.equiv          nb114nf_rinvMM,         848
.equiv          nb114nf_is3,            864
.equiv          nb114nf_ii3,            868
.equiv          nb114nf_innerjjnr,      872
.equiv          nb114nf_innerk,         876
.equiv          nb114nf_n,              880
.equiv          nb114nf_nn1,            884
.equiv          nb114nf_nri,            888
.equiv          nb114nf_nouter,         892
.equiv          nb114nf_ninner,         896
.equiv          nb114nf_salign,         900
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 904		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb114nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb114nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb114nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb114nf_nouter], eax
	mov [esp + nb114nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb114nf_half], eax
	mov [esp + nb114nf_half+4], ebx
	movsd xmm1, [esp + nb114nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb114nf_half], xmm1
	movapd [esp + nb114nf_two], xmm2
	movapd [esp + nb114nf_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb114nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb114nf_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb114nf_p_facel]
	movsd xmm6, [esi]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb114nf_qqMM], xmm3
	movapd [esp + nb114nf_qqMH], xmm4
	movapd [esp + nb114nf_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb114nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb114nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb114nf_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movhpd xmm0, [eax + edx*8 + 8]
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb114nf_c6], xmm0
	movapd [esp + nb114nf_c12], xmm1

.nb114nf_threadloop:
        mov   esi, [ebp + nb114nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb114nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb114nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb114nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb114nf_n], eax
        mov [esp + nb114nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb114nf_outerstart
        jmp .nb114nf_end

.nb114nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb114nf_nouter]
	mov [esp + nb114nf_nouter], ebx

.nb114nf_outer:
	mov   eax, [ebp + nb114nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb114nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb114nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb114nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb114nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb114nf_ii3], ebx		

	addsd xmm3, [eax + ebx*8] 	;# ox
	addsd xmm4, [eax + ebx*8 + 8] 	;# oy
	addsd xmm5, [eax + ebx*8 + 16]	;# oz	
	addsd xmm6, [eax + ebx*8 + 24] 	;# h1x
	addsd xmm7, [eax + ebx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [esp + nb114nf_ixO], xmm3
	movapd [esp + nb114nf_iyO], xmm4
	movapd [esp + nb114nf_izO], xmm5
	movapd [esp + nb114nf_ixH1], xmm6
	movapd [esp + nb114nf_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [eax + ebx*8 + 40] ;# h1z
	addsd xmm0, [eax + ebx*8 + 48] ;# h2x
	addsd xmm1, [eax + ebx*8 + 56] ;# h2y
	addsd xmm2, [eax + ebx*8 + 64] ;# h2z
	addsd xmm3, [eax + ebx*8 + 72] ;# mx
	addsd xmm4, [eax + ebx*8 + 80] ;# my
	addsd xmm5, [eax + ebx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb114nf_izH1], xmm6
	movapd [esp + nb114nf_ixH2], xmm0
	movapd [esp + nb114nf_iyH2], xmm1
	movapd [esp + nb114nf_izH2], xmm2
	movapd [esp + nb114nf_ixM], xmm3
	movapd [esp + nb114nf_iyM], xmm4
	movapd [esp + nb114nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb114nf_vctot], xmm4
	movapd [esp + nb114nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb114nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb114nf_pos] 
	mov   edi, [ebp + nb114nf_faction]	
	mov   eax, [ebp + nb114nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb114nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb114nf_ninner]
	mov   [esp + nb114nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb114nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb114nf_unroll_loop
	jmp   .nb114nf_checksingle
.nb114nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb114nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb114nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb114nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8]
	movlpd xmm2, [esi + ebx*8]
	movhpd xmm0, [esi + eax*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 8]
	movlpd xmm3, [esi + eax*8 + 16]
	movlpd xmm5, [esi + ebx*8 + 16]
	movhpd xmm3, [esi + eax*8 + 24]
	movhpd xmm5, [esi + ebx*8 + 24]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# ox 
	unpckhpd xmm1, xmm2 ;# oy
	unpcklpd xmm3, xmm5 ;# ox 
	unpckhpd xmm4, xmm5 ;# oy
	movapd 	[esp + nb114nf_jxO], xmm0
	movapd 	[esp + nb114nf_jyO], xmm1
	movapd 	[esp + nb114nf_jzO], xmm3
	movapd 	[esp + nb114nf_jxH1], xmm4
	
	;# load h1y, h1z, h2x, h2y 
	movlpd xmm0, [esi + eax*8 + 32]
	movlpd xmm2, [esi + ebx*8 + 32]
	movhpd xmm0, [esi + eax*8 + 40]
	movhpd xmm2, [esi + ebx*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movlpd xmm5, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movhpd xmm5, [esi + ebx*8 + 56]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h1y
	unpckhpd xmm1, xmm2 ;# h1z
	unpcklpd xmm3, xmm5 ;# h2x
	unpckhpd xmm4, xmm5 ;# h2y
	movapd 	[esp + nb114nf_jyH1], xmm0
	movapd 	[esp + nb114nf_jzH1], xmm1
	movapd 	[esp + nb114nf_jxH2], xmm3
	movapd 	[esp + nb114nf_jyH2], xmm4
	
	;# load h2z, mx, my, mz
	movlpd xmm0, [esi + eax*8 + 64]
	movlpd xmm2, [esi + ebx*8 + 64]
	movhpd xmm0, [esi + eax*8 + 72]
	movhpd xmm2, [esi + ebx*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm5, [esi + ebx*8 + 80]
	movhpd xmm3, [esi + eax*8 + 88]
	movhpd xmm5, [esi + ebx*8 + 88]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h2z
	unpckhpd xmm1, xmm2 ;# mx
	unpcklpd xmm3, xmm5 ;# my
	unpckhpd xmm4, xmm5 ;# mz
	movapd 	[esp + nb114nf_jzH2], xmm0
	movapd 	[esp + nb114nf_jxM], xmm1
	movapd 	[esp + nb114nf_jyM], xmm3
	movapd 	[esp + nb114nf_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [esp + nb114nf_ixO]
	movapd xmm1, [esp + nb114nf_iyO]
	movapd xmm2, [esp + nb114nf_izO]
	movapd xmm3, [esp + nb114nf_ixH1]
	movapd xmm4, [esp + nb114nf_iyH1]
	movapd xmm5, [esp + nb114nf_izH1]
	subpd  xmm0, [esp + nb114nf_jxO]
	subpd  xmm1, [esp + nb114nf_jyO]
	subpd  xmm2, [esp + nb114nf_jzO]
	subpd  xmm3, [esp + nb114nf_jxH1]
	subpd  xmm4, [esp + nb114nf_jyH1]
	subpd  xmm5, [esp + nb114nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb114nf_rsqOO], xmm0
	movapd [esp + nb114nf_rsqH1H1], xmm3

	movapd xmm0, [esp + nb114nf_ixH1]
	movapd xmm1, [esp + nb114nf_iyH1]
	movapd xmm2, [esp + nb114nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb114nf_jxH2]
	subpd  xmm1, [esp + nb114nf_jyH2]
	subpd  xmm2, [esp + nb114nf_jzH2]
	subpd  xmm3, [esp + nb114nf_jxM]
	subpd  xmm4, [esp + nb114nf_jyM]
	subpd  xmm5, [esp + nb114nf_jzM]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb114nf_rsqH1H2], xmm0
	movapd [esp + nb114nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb114nf_ixH2]
	movapd xmm1, [esp + nb114nf_iyH2]
	movapd xmm2, [esp + nb114nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb114nf_jxH1]
	subpd  xmm1, [esp + nb114nf_jyH1]
	subpd  xmm2, [esp + nb114nf_jzH1]
	subpd  xmm3, [esp + nb114nf_jxH2]
	subpd  xmm4, [esp + nb114nf_jyH2]
	subpd  xmm5, [esp + nb114nf_jzH2]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb114nf_rsqH2H1], xmm0
	movapd [esp + nb114nf_rsqH2H2], xmm3

	movapd xmm0, [esp + nb114nf_ixH2]
	movapd xmm1, [esp + nb114nf_iyH2]
	movapd xmm2, [esp + nb114nf_izH2]
	movapd xmm3, [esp + nb114nf_ixM]
	movapd xmm4, [esp + nb114nf_iyM]
	movapd xmm5, [esp + nb114nf_izM]
	subpd  xmm0, [esp + nb114nf_jxM]
	subpd  xmm1, [esp + nb114nf_jyM]
	subpd  xmm2, [esp + nb114nf_jzM]
	subpd  xmm3, [esp + nb114nf_jxH1]
	subpd  xmm4, [esp + nb114nf_jyH1]
	subpd  xmm5, [esp + nb114nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb114nf_rsqH2M], xmm0
	movapd [esp + nb114nf_rsqMH1], xmm4

	movapd xmm0, [esp + nb114nf_ixM]
	movapd xmm1, [esp + nb114nf_iyM]
	movapd xmm2, [esp + nb114nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb114nf_jxH2]
	subpd  xmm1, [esp + nb114nf_jyH2]
	subpd  xmm2, [esp + nb114nf_jzH2]
	subpd  xmm3, [esp + nb114nf_jxM]
	subpd  xmm4, [esp + nb114nf_jyM]
	subpd  xmm5, [esp + nb114nf_jzM]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb114nf_rsqMH2], xmm0
	movapd [esp + nb114nf_rsqMM], xmm4
	
	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtpd2ps xmm1, xmm0
	cvtpd2ps xmm5, xmm4
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1   ;# luA
	cvtps2pd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb114nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvMH2], xmm1
	movapd [esp + nb114nf_rinvMM], xmm5

	movapd xmm0, [esp + nb114nf_rsqOO]
	movapd xmm4, [esp + nb114nf_rsqH1H1]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb114nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb114nf_half] ;# rinv
	mulpd   xmm1, xmm1
	movapd [esp + nb114nf_rinvsqOO], xmm1
	movapd [esp + nb114nf_rinvH1H1], xmm5

	movapd xmm0, [esp + nb114nf_rsqH1H2]
	movapd xmm4, [esp + nb114nf_rsqH1M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb114nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvH1H2], xmm1
	movapd [esp + nb114nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb114nf_rsqH2H1]
	movapd xmm4, [esp + nb114nf_rsqH2H2]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb114nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvH2H1], xmm1
	movapd [esp + nb114nf_rinvH2H2], xmm5

	movapd xmm0, [esp + nb114nf_rsqMH1]
	movapd xmm4, [esp + nb114nf_rsqH2M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb114nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb114nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvMH1], xmm1
	movapd [esp + nb114nf_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb114nf_rinvsqOO] ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulpd   xmm1, xmm1 ;# rinv4
	mulpd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulpd	xmm2, xmm2 ;# rinvtwelve
	mulpd  xmm1, [esp + nb114nf_c6]
	mulpd  xmm2, [esp + nb114nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb114nf_Vvdwtot]
	movapd [esp + nb114nf_Vvdwtot], xmm3

	;# Coulomb interactions.
	;# Sum rinv for H-H interactions in xmm0, H-M in xmm1.
	;# Use xmm2 for the M-M pair
	movapd xmm0, [esp + nb114nf_rinvH1H1]
	movapd xmm1, [esp + nb114nf_rinvH1M]
	movapd xmm2, [esp + nb114nf_rinvMM]

	addpd  xmm0, [esp + nb114nf_rinvH1H2]
	addpd  xmm1, [esp + nb114nf_rinvH2M]
	addpd  xmm0, [esp + nb114nf_rinvH2H1]
	addpd  xmm1, [esp + nb114nf_rinvMH1]
	addpd  xmm0, [esp + nb114nf_rinvH2H2]
	addpd  xmm1, [esp + nb114nf_rinvMH2]
	mulpd  xmm0, [esp + nb114nf_qqHH]
	mulpd  xmm1, [esp + nb114nf_qqMH]
	mulpd  xmm2, [esp + nb114nf_qqMM]
	addpd  xmm0, xmm1
	addpd  xmm2, [esp + nb114nf_vctot]
	addpd  xmm2, xmm0
	movapd [esp + nb114nf_vctot], xmm2

	;# should we do one more iteration? 
	sub dword ptr [esp + nb114nf_innerk],  2
	jl    .nb114nf_checksingle
	jmp   .nb114nf_unroll_loop
.nb114nf_checksingle:
	mov   edx, [esp + nb114nf_innerk]
	and   edx, 1
	jnz   .nb114nf_dosingle
	jmp   .nb114nf_updateouterdata
.nb114nf_dosingle:
	mov   edx, [esp + nb114nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb114nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [esi + eax*8] 
	movhpd xmm0, [esi + eax*8 + 8]
	movlpd xmm1, [esi + eax*8 + 16]
	movhpd xmm1, [esi + eax*8 + 24]
	movlpd xmm2, [esi + eax*8 + 32]
	movhpd xmm2, [esi + eax*8 + 40]
	movlpd xmm3, [esi + eax*8 + 48]
	movhpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm4, [esi + eax*8 + 72]
	movlpd xmm5, [esi + eax*8 + 80]
	movhpd xmm5, [esi + eax*8 + 88]
	movsd  [esp + nb114nf_jxO], xmm0
	movsd  [esp + nb114nf_jzO], xmm1
	movsd  [esp + nb114nf_jyH1], xmm2
	movsd  [esp + nb114nf_jxH2], xmm3
	movsd  [esp + nb114nf_jzH2], xmm4
	movsd  [esp + nb114nf_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [esp + nb114nf_jyO], xmm0
	movsd  [esp + nb114nf_jxH1], xmm1
	movsd  [esp + nb114nf_jzH1], xmm2
	movsd  [esp + nb114nf_jyH2], xmm3
	movsd  [esp + nb114nf_jxM], xmm4
	movsd  [esp + nb114nf_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [esp + nb114nf_ixO]
	movapd xmm1, [esp + nb114nf_iyO]
	movapd xmm2, [esp + nb114nf_izO]
	movapd xmm3, [esp + nb114nf_ixH1]
	movapd xmm4, [esp + nb114nf_iyH1]
	movapd xmm5, [esp + nb114nf_izH1]
	subsd  xmm0, [esp + nb114nf_jxO]
	subsd  xmm1, [esp + nb114nf_jyO]
	subsd  xmm2, [esp + nb114nf_jzO]
	subsd  xmm3, [esp + nb114nf_jxH1]
	subsd  xmm4, [esp + nb114nf_jyH1]
	subsd  xmm5, [esp + nb114nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb114nf_rsqOO], xmm0
	movapd [esp + nb114nf_rsqH1H1], xmm3

	movapd xmm0, [esp + nb114nf_ixH1]
	movapd xmm1, [esp + nb114nf_iyH1]
	movapd xmm2, [esp + nb114nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb114nf_jxH2]
	subsd  xmm1, [esp + nb114nf_jyH2]
	subsd  xmm2, [esp + nb114nf_jzH2]
	subsd  xmm3, [esp + nb114nf_jxM]
	subsd  xmm4, [esp + nb114nf_jyM]
	subsd  xmm5, [esp + nb114nf_jzM]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb114nf_rsqH1H2], xmm0
	movapd [esp + nb114nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb114nf_ixH2]
	movapd xmm1, [esp + nb114nf_iyH2]
	movapd xmm2, [esp + nb114nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb114nf_jxH1]
	subsd  xmm1, [esp + nb114nf_jyH1]
	subsd  xmm2, [esp + nb114nf_jzH1]
	subsd  xmm3, [esp + nb114nf_jxH2]
	subsd  xmm4, [esp + nb114nf_jyH2]
	subsd  xmm5, [esp + nb114nf_jzH2]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb114nf_rsqH2H1], xmm0
	movapd [esp + nb114nf_rsqH2H2], xmm3

	movapd xmm0, [esp + nb114nf_ixH2]
	movapd xmm1, [esp + nb114nf_iyH2]
	movapd xmm2, [esp + nb114nf_izH2]
	movapd xmm3, [esp + nb114nf_ixM]
	movapd xmm4, [esp + nb114nf_iyM]
	movapd xmm5, [esp + nb114nf_izM]
	subsd  xmm0, [esp + nb114nf_jxM]
	subsd  xmm1, [esp + nb114nf_jyM]
	subsd  xmm2, [esp + nb114nf_jzM]
	subsd  xmm3, [esp + nb114nf_jxH1]
	subsd  xmm4, [esp + nb114nf_jyH1]
	subsd  xmm5, [esp + nb114nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb114nf_rsqH2M], xmm0
	movapd [esp + nb114nf_rsqMH1], xmm4

	movapd xmm0, [esp + nb114nf_ixM]
	movapd xmm1, [esp + nb114nf_iyM]
	movapd xmm2, [esp + nb114nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb114nf_jxH2]
	subsd  xmm1, [esp + nb114nf_jyH2]
	subsd  xmm2, [esp + nb114nf_jzH2]
	subsd  xmm3, [esp + nb114nf_jxM]
	subsd  xmm4, [esp + nb114nf_jyM]
	subsd  xmm5, [esp + nb114nf_jzM]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb114nf_rsqMH2], xmm0
	movapd [esp + nb114nf_rsqMM], xmm4

	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtsd2ss xmm1, xmm0
	cvtsd2ss xmm5, xmm4
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1   ;# luA
	cvtss2sd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb114nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvMH2], xmm1
	movapd [esp + nb114nf_rinvMM], xmm5

	movapd xmm0, [esp + nb114nf_rsqOO]
	movapd xmm4, [esp + nb114nf_rsqH1H1]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb114nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb114nf_half] ;# rinv
	mulpd   xmm1, xmm1
	movapd [esp + nb114nf_rinvsqOO], xmm1
	movapd [esp + nb114nf_rinvH1H1], xmm5

	movapd xmm0, [esp + nb114nf_rsqH1H2]
	movapd xmm4, [esp + nb114nf_rsqH1M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb114nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvH1H2], xmm1
	movapd [esp + nb114nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb114nf_rsqH2H1]
	movapd xmm4, [esp + nb114nf_rsqH2H2]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb114nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvH2H1], xmm1
	movapd [esp + nb114nf_rinvH2H2], xmm5

	movapd xmm0, [esp + nb114nf_rsqMH1]
	movapd xmm4, [esp + nb114nf_rsqH2M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [esp + nb114nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb114nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb114nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb114nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb114nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb114nf_half] ;# rinv 
	movapd [esp + nb114nf_rinvMH1], xmm1
	movapd [esp + nb114nf_rinvH2M], xmm5

	;# start with OO interaction 
	movsd xmm0, [esp + nb114nf_rinvsqOO] ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movsd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [esp + nb114nf_c6]
	mulsd  xmm2, [esp + nb114nf_c12]
	movsd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb114nf_Vvdwtot]
	movsd [esp + nb114nf_Vvdwtot], xmm3

	;# Coulomb interactions.
	;# Sum rinv for H-H interactions in xmm0, H-M in xmm1.
	;# Use xmm2 for the M-M pair
	movsd xmm0, [esp + nb114nf_rinvH1H1]
	movsd xmm1, [esp + nb114nf_rinvH1M]
	movsd xmm2, [esp + nb114nf_rinvMM]

	addsd  xmm0, [esp + nb114nf_rinvH1H2]
	addsd  xmm1, [esp + nb114nf_rinvH2M]
	addsd  xmm0, [esp + nb114nf_rinvH2H1]
	addsd  xmm1, [esp + nb114nf_rinvMH1]
	addsd  xmm0, [esp + nb114nf_rinvH2H2]
	addsd  xmm1, [esp + nb114nf_rinvMH2]
	mulsd  xmm0, [esp + nb114nf_qqHH]
	mulsd  xmm1, [esp + nb114nf_qqMH]
	mulsd  xmm2, [esp + nb114nf_qqMM]
	addsd  xmm0, xmm1
	addsd  xmm2, [esp + nb114nf_vctot]
	addsd  xmm2, xmm0
	movsd [esp + nb114nf_vctot], xmm2

	
.nb114nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb114nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb114nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb114nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb114nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb114nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb114nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb114nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb114nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb114nf_n], esi
        jmp .nb114nf_outer
.nb114nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb114nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb114nf_end
        ;# non-zero, do one more workunit
        jmp   .nb114nf_threadloop
.nb114nf_end:
	emms

	mov eax, [esp + nb114nf_nouter]
	mov ebx, [esp + nb114nf_ninner]
	mov ecx, [ebp + nb114nf_outeriter]
	mov edx, [ebp + nb114nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb114nf_salign]
	add esp, eax
	add esp, 904
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


