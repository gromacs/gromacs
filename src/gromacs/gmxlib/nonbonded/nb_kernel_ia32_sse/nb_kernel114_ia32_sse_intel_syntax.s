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


	
.globl nb_kernel114_ia32_sse 
.globl _nb_kernel114_ia32_sse
nb_kernel114_ia32_sse:	
_nb_kernel114_ia32_sse:	
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
.equiv          nb114_p_krf,            56
.equiv          nb114_p_crf,            60
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
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb114_c6,               928
.equiv          nb114_c12,              944
.equiv          nb114_six,              960
.equiv          nb114_twelve,           976
.equiv          nb114_vctot,            992
.equiv          nb114_Vvdwtot,          1008
.equiv          nb114_fixO,             1024
.equiv          nb114_fiyO,             1040
.equiv          nb114_fizO,             1056
.equiv          nb114_fixH1,            1072
.equiv          nb114_fiyH1,            1088
.equiv          nb114_fizH1,            1104
.equiv          nb114_fixH2,            1120
.equiv          nb114_fiyH2,            1136
.equiv          nb114_fizH2,            1152
.equiv          nb114_fixM,             1168
.equiv          nb114_fiyM,             1184
.equiv          nb114_fizM,             1200
.equiv          nb114_fjxO,             1216
.equiv          nb114_fjyO,             1232
.equiv          nb114_fjzO,             1248
.equiv          nb114_fjxH1,            1264
.equiv          nb114_fjyH1,            1280
.equiv          nb114_fjzH1,            1296
.equiv          nb114_fjxH2,            1312
.equiv          nb114_fjyH2,            1328
.equiv          nb114_fjzH2,            1344
.equiv          nb114_fjxM,             1360
.equiv          nb114_fjyM,             1376
.equiv          nb114_fjzM,             1392
.equiv          nb114_half,             1408
.equiv          nb114_three,            1424
.equiv          nb114_rsqOO,            1440
.equiv          nb114_rsqH1H1,          1456
.equiv          nb114_rsqH1H2,          1472
.equiv          nb114_rsqH1M,           1488
.equiv          nb114_rsqH2H1,          1504
.equiv          nb114_rsqH2H2,          1520
.equiv          nb114_rsqH2M,           1536
.equiv          nb114_rsqMH1,           1552
.equiv          nb114_rsqMH2,           1568
.equiv          nb114_rsqMM,            1584
.equiv          nb114_rinvsqOO,         1600
.equiv          nb114_rinvH1H1,         1616
.equiv          nb114_rinvH1H2,         1632
.equiv          nb114_rinvH1M,          1648
.equiv          nb114_rinvH2H1,         1664
.equiv          nb114_rinvH2H2,         1680
.equiv          nb114_rinvH2M,          1696
.equiv          nb114_rinvMH1,          1712
.equiv          nb114_rinvMH2,          1728
.equiv          nb114_rinvMM,           1744
.equiv          nb114_fstmp,            1760
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
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb114_half], eax
	movss xmm1, [esp + nb114_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# 6.0
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# constant 12.0
	movaps [esp + nb114_half],  xmm1
	movaps [esp + nb114_two],  xmm2
	movaps [esp + nb114_three],  xmm3
	movaps [esp + nb114_six],  xmm4
	movaps [esp + nb114_twelve],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb114_iinr]   ;# ecx = pointer into iinr[]
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb114_charge]
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	movss xmm4, xmm3	
	mov esi, [ebp + nb114_p_facel]
	movss xmm6, [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb114_qqMM], xmm3
	movaps [esp + nb114_qqMH], xmm4
	movaps [esp + nb114_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb114_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb114_p_ntype]
	imul  ecx, [edi]  ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb114_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [esp + nb114_c6], xmm0
	movaps [esp + nb114_c12], xmm1

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
	mov   eax, [ebp + nb114_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb114_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb114_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb114_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb114_pos]	;# eax = base of pos[]  
	mov   [esp + nb114_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	movaps xmm6, xmm0
	movaps xmm7, xmm1

	addss xmm3, [eax + ebx*4]  	;# ox
	addss xmm4, [eax + ebx*4 + 4]  ;# oy
	addss xmm5, [eax + ebx*4 + 8]  ;# oz
	addss xmm6, [eax + ebx*4 + 12] ;# h1x
	addss xmm7, [eax + ebx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [esp + nb114_ixO], xmm3
	movaps [esp + nb114_iyO], xmm4
	movaps [esp + nb114_izO], xmm5
	movaps [esp + nb114_ixH1], xmm6
	movaps [esp + nb114_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [eax + ebx*4 + 20] ;# h1z
	addss xmm0, [eax + ebx*4 + 24] ;# h2x
	addss xmm1, [eax + ebx*4 + 28] ;# h2y
	addss xmm2, [eax + ebx*4 + 32] ;# h2z
	addss xmm3, [eax + ebx*4 + 36] ;# mx
	addss xmm4, [eax + ebx*4 + 40] ;# my
	addss xmm5, [eax + ebx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb114_izH1], xmm6
	movaps [esp + nb114_ixH2], xmm0
	movaps [esp + nb114_iyH2], xmm1
	movaps [esp + nb114_izH2], xmm2
	movaps [esp + nb114_ixM], xmm3
	movaps [esp + nb114_iyM], xmm4
	movaps [esp + nb114_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb114_vctot], xmm4
	movaps [esp + nb114_Vvdwtot], xmm4
	movaps [esp + nb114_fixO], xmm4
	movaps [esp + nb114_fiyO], xmm4
	movaps [esp + nb114_fizO], xmm4
	movaps [esp + nb114_fixH1], xmm4
	movaps [esp + nb114_fiyH1], xmm4
	movaps [esp + nb114_fizH1], xmm4
	movaps [esp + nb114_fixH2], xmm4
	movaps [esp + nb114_fiyH2], xmm4
	movaps [esp + nb114_fizH2], xmm4
	movaps [esp + nb114_fixM], xmm4
	movaps [esp + nb114_fiyM], xmm4
	movaps [esp + nb114_fizM], xmm4
	
	mov   eax, [ebp + nb114_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb114_pos]
	mov   edi, [ebp + nb114_faction]	
	mov   eax, [ebp + nb114_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb114_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb114_ninner]
	mov   [esp + nb114_ninner], ecx
	add   edx, 0
	mov   [esp + nb114_innerk], edx	;# number of innerloop atoms 
	jge   .nb114_unroll_loop
	jmp   .nb114_single_check
.nb114_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb114_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb114_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov esi, [ebp + nb114_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables.
	;# Load Ox Oy Oz H1x
	movlps xmm1, [esi + eax*4]	;#  Oxa   Oya    -    -
	movlps xmm4, [esi + ecx*4]	;#  Oxc   Oyc    -    -
	movhps xmm1, [esi + ebx*4]	;#  Oxa   Oya   Oxb   Oyb 
	movhps xmm4, [esi + edx*4]	;#  Oxc   Oyc   Oxd   Oyd 
	movaps xmm0, xmm1		;#  Oxa   Oya   Oxb   Oyb 
	shufps xmm0, xmm4, 0x88		;#  Oxa   Oxb   Oxc   Oxd
	shufps xmm1, xmm4, 0xDD		;#  Oya   Oyb   Oyc   Oyd
	movlps xmm3, [esi + eax*4 + 8]	;#  Oza  H1xa    -    -
	movlps xmm5, [esi + ecx*4 + 8]	;#  Ozc  H1xc    -    -
	movhps xmm3, [esi + ebx*4 + 8]	;#  Oza  H1xa   Ozb  H1xb 
	movhps xmm5, [esi + edx*4 + 8]	;#  Ozc  H1xc   Ozd  H1xd 
	movaps xmm2, xmm3		;#  Oza  H1xa   Ozb  H1xb 
	shufps xmm2, xmm5, 0x88		;#  Oza   Ozb   Ozc   Ozd
	shufps xmm3, xmm5, 0xDD		;# H1xa  H1xb  H1xc  H1xd
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb114_jxO], xmm0
	movaps [esp + nb114_jyO], xmm1
	movaps [esp + nb114_jzO], xmm2
	movaps [esp + nb114_jxH1], xmm3

	;# Load H1y H1z H2x H2y 
	movlps xmm1, [esi + eax*4 + 16]	
	movlps xmm4, [esi + ecx*4 + 16]	
	movhps xmm1, [esi + ebx*4 + 16]	
	movhps xmm4, [esi + edx*4 + 16]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 24]	
	movlps xmm5, [esi + ecx*4 + 24]	
	movhps xmm3, [esi + ebx*4 + 24]	
	movhps xmm5, [esi + edx*4 + 24]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb114_jyH1], xmm0
	movaps [esp + nb114_jzH1], xmm1
	movaps [esp + nb114_jxH2], xmm2
	movaps [esp + nb114_jyH2], xmm3

	;# Load H2z Mx My Mz 
	movlps xmm1, [esi + eax*4 + 32]	
	movlps xmm4, [esi + ecx*4 + 32]	
	movhps xmm1, [esi + ebx*4 + 32]	
	movhps xmm4, [esi + edx*4 + 32]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 40]	
	movlps xmm5, [esi + ecx*4 + 40]	
	movhps xmm3, [esi + ebx*4 + 40]	
	movhps xmm5, [esi + edx*4 + 40]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb114_jzH2], xmm0
	movaps [esp + nb114_jxM], xmm1
	movaps [esp + nb114_jyM], xmm2
	movaps [esp + nb114_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [esp + nb114_ixO]
	movaps xmm1, [esp + nb114_iyO]
	movaps xmm2, [esp + nb114_izO]
	movaps xmm3, [esp + nb114_ixH1]
	movaps xmm4, [esp + nb114_iyH1]
	movaps xmm5, [esp + nb114_izH1]
	subps  xmm0, [esp + nb114_jxO]
	subps  xmm1, [esp + nb114_jyO]
	subps  xmm2, [esp + nb114_jzO]
	subps  xmm3, [esp + nb114_jxH1]
	subps  xmm4, [esp + nb114_jyH1]
	subps  xmm5, [esp + nb114_jzH1]
	movaps [esp + nb114_dxOO], xmm0
	movaps [esp + nb114_dyOO], xmm1
	movaps [esp + nb114_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb114_dxH1H1], xmm3
	movaps [esp + nb114_dyH1H1], xmm4
	movaps [esp + nb114_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb114_rsqOO], xmm0
	movaps [esp + nb114_rsqH1H1], xmm3

	movaps xmm0, [esp + nb114_ixH1]
	movaps xmm1, [esp + nb114_iyH1]
	movaps xmm2, [esp + nb114_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb114_jxH2]
	subps  xmm1, [esp + nb114_jyH2]
	subps  xmm2, [esp + nb114_jzH2]
	subps  xmm3, [esp + nb114_jxM]
	subps  xmm4, [esp + nb114_jyM]
	subps  xmm5, [esp + nb114_jzM]
	movaps [esp + nb114_dxH1H2], xmm0
	movaps [esp + nb114_dyH1H2], xmm1
	movaps [esp + nb114_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb114_dxH1M], xmm3
	movaps [esp + nb114_dyH1M], xmm4
	movaps [esp + nb114_dzH1M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb114_rsqH1H2], xmm0
	movaps [esp + nb114_rsqH1M], xmm3

	movaps xmm0, [esp + nb114_ixH2]
	movaps xmm1, [esp + nb114_iyH2]
	movaps xmm2, [esp + nb114_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb114_jxH1]
	subps  xmm1, [esp + nb114_jyH1]
	subps  xmm2, [esp + nb114_jzH1]
	subps  xmm3, [esp + nb114_jxH2]
	subps  xmm4, [esp + nb114_jyH2]
	subps  xmm5, [esp + nb114_jzH2]
	movaps [esp + nb114_dxH2H1], xmm0
	movaps [esp + nb114_dyH2H1], xmm1
	movaps [esp + nb114_dzH2H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb114_dxH2H2], xmm3
	movaps [esp + nb114_dyH2H2], xmm4
	movaps [esp + nb114_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb114_rsqH2H1], xmm0
	movaps [esp + nb114_rsqH2H2], xmm3

	movaps xmm0, [esp + nb114_ixH2]
	movaps xmm1, [esp + nb114_iyH2]
	movaps xmm2, [esp + nb114_izH2]
	movaps xmm3, [esp + nb114_ixM]
	movaps xmm4, [esp + nb114_iyM]
	movaps xmm5, [esp + nb114_izM]
	subps  xmm0, [esp + nb114_jxM]
	subps  xmm1, [esp + nb114_jyM]
	subps  xmm2, [esp + nb114_jzM]
	subps  xmm3, [esp + nb114_jxH1]
	subps  xmm4, [esp + nb114_jyH1]
	subps  xmm5, [esp + nb114_jzH1]
	movaps [esp + nb114_dxH2M], xmm0
	movaps [esp + nb114_dyH2M], xmm1
	movaps [esp + nb114_dzH2M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb114_dxMH1], xmm3
	movaps [esp + nb114_dyMH1], xmm4
	movaps [esp + nb114_dzMH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb114_rsqH2M], xmm0
	movaps [esp + nb114_rsqMH1], xmm4

	movaps xmm0, [esp + nb114_ixM]
	movaps xmm1, [esp + nb114_iyM]
	movaps xmm2, [esp + nb114_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb114_jxH2]
	subps  xmm1, [esp + nb114_jyH2]
	subps  xmm2, [esp + nb114_jzH2]
	subps  xmm3, [esp + nb114_jxM]
	subps  xmm4, [esp + nb114_jyM]
	subps  xmm5, [esp + nb114_jzM]
	movaps [esp + nb114_dxMH2], xmm0
	movaps [esp + nb114_dyMH2], xmm1
	movaps [esp + nb114_dzMH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb114_dxMM], xmm3
	movaps [esp + nb114_dyMM], xmm4
	movaps [esp + nb114_dzMM], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb114_rsqMH2], xmm0
	movaps [esp + nb114_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [esp + nb114_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [esp + nb114_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [esp + nb114_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [esp + nb114_rsqH1H1]
	rsqrtps xmm5, [esp + nb114_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114_rsqH1H1]
	mulps   xmm5, [esp + nb114_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114_half] ;# rinvH1H1 
	mulps   xmm7, [esp + nb114_half] ;# rinvH1H2 
	movaps  [esp + nb114_rinvH1H1], xmm3
	movaps  [esp + nb114_rinvH1H2], xmm7
			
	rsqrtps xmm1, [esp + nb114_rsqH1M]
	rsqrtps xmm5, [esp + nb114_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114_rsqH1M]
	mulps   xmm5, [esp + nb114_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114_half] 
	mulps   xmm7, [esp + nb114_half]
	movaps  [esp + nb114_rinvH1M], xmm3
	movaps  [esp + nb114_rinvH2H1], xmm7
			
	rsqrtps xmm1, [esp + nb114_rsqH2H2]
	rsqrtps xmm5, [esp + nb114_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114_rsqH2H2]
	mulps   xmm5, [esp + nb114_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114_half] 
	mulps   xmm7, [esp + nb114_half]
	movaps  [esp + nb114_rinvH2H2], xmm3
	movaps  [esp + nb114_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb114_rsqMH1]
	rsqrtps xmm5, [esp + nb114_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114_rsqMH1]
	mulps   xmm5, [esp + nb114_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114_half] 
	mulps   xmm7, [esp + nb114_half]
	movaps  [esp + nb114_rinvMH1], xmm3
	movaps  [esp + nb114_rinvMH2], xmm7
        		
	rsqrtps xmm1, [esp + nb114_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb114_three]
	mulps   xmm1, [esp + nb114_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb114_half] 
	movaps  [esp + nb114_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [esp + nb114_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb114_c6]
	mulps  xmm2, [esp + nb114_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [esp + nb114_Vvdwtot]
	mulps  xmm1, [esp + nb114_six]
	mulps  xmm2, [esp + nb114_twelve]
	movaps [esp + nb114_Vvdwtot], xmm4
	subps  xmm2, xmm1
	mulps  xmm0, xmm2 	;# fscal 
	movaps xmm1, xmm0
	movaps xmm2, xmm0		

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb114_dxOO]
	mulps xmm1, [esp + nb114_dyOO]
	mulps xmm2, [esp + nb114_dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixO]
	addps xmm1, [esp + nb114_fiyO]
	addps xmm2, [esp + nb114_fizO]
	movaps [esp + nb114_fjxO], xmm3
	movaps [esp + nb114_fjyO], xmm4
	movaps [esp + nb114_fjzO], xmm5
	movaps [esp + nb114_fixO], xmm0
	movaps [esp + nb114_fiyO], xmm1
	movaps [esp + nb114_fizO], xmm2

	;# Coulomb interactions 
	;# start with H1-H1 interaction 
	movaps xmm0, [esp + nb114_rinvH1H1]
	movaps xmm7, xmm0
	mulps  xmm0, xmm0
	mulps  xmm7, [esp + nb114_qqHH]
	mulps  xmm0, xmm7	
	addps  xmm7, [esp + nb114_vctot] 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb114_dxH1H1]
	mulps xmm1, [esp + nb114_dyH1H1]
	mulps xmm2, [esp + nb114_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixH1]
	addps xmm1, [esp + nb114_fiyH1]
	addps xmm2, [esp + nb114_fizH1]
	movaps [esp + nb114_fjxH1], xmm3
	movaps [esp + nb114_fjyH1], xmm4
	movaps [esp + nb114_fjzH1], xmm5
	movaps [esp + nb114_fixH1], xmm0
	movaps [esp + nb114_fiyH1], xmm1
	movaps [esp + nb114_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [esp + nb114_rinvH1H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqHH]
	mulps xmm0, xmm1	;# fs H1-H2  
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb114_dxH1H2]
	mulps xmm1, [esp + nb114_dyH1H2]
	mulps xmm2, [esp + nb114_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixH1]
	addps xmm1, [esp + nb114_fiyH1]
	addps xmm2, [esp + nb114_fizH1]
	movaps [esp + nb114_fjxH2], xmm3
	movaps [esp + nb114_fjyH2], xmm4
	movaps [esp + nb114_fjzH2], xmm5
	movaps [esp + nb114_fixH1], xmm0
	movaps [esp + nb114_fiyH1], xmm1
	movaps [esp + nb114_fizH1], xmm2

	;# H1-M interaction  
	movaps xmm0, [esp + nb114_rinvH1M]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqMH]
	mulps xmm0, xmm1	;# fs H1-M  
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb114_dxH1M]
	mulps xmm1, [esp + nb114_dyH1M]
	mulps xmm2, [esp + nb114_dzH1M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixH1]
	addps xmm1, [esp + nb114_fiyH1]
	addps xmm2, [esp + nb114_fizH1]
	movaps [esp + nb114_fjxM], xmm3
	movaps [esp + nb114_fjyM], xmm4
	movaps [esp + nb114_fjzM], xmm5
	movaps [esp + nb114_fixH1], xmm0
	movaps [esp + nb114_fiyH1], xmm1
	movaps [esp + nb114_fizH1], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [esp + nb114_rinvH2H1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqHH]
	mulps xmm0, xmm1	;# fs H2-H1 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + nb114_fjxH1]
	movaps xmm4, [esp + nb114_fjyH1]
	movaps xmm5, [esp + nb114_fjzH1]
	mulps xmm0, [esp + nb114_dxH2H1]
	mulps xmm1, [esp + nb114_dyH2H1]
	mulps xmm2, [esp + nb114_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixH2]
	addps xmm1, [esp + nb114_fiyH2]
	addps xmm2, [esp + nb114_fizH2]
	movaps [esp + nb114_fjxH1], xmm3
	movaps [esp + nb114_fjyH1], xmm4
	movaps [esp + nb114_fjzH1], xmm5
	movaps [esp + nb114_fixH2], xmm0
	movaps [esp + nb114_fiyH2], xmm1
	movaps [esp + nb114_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [esp + nb114_rinvH2H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqHH]
	mulps xmm0, xmm1	;# fsH2H2
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + nb114_fjxH2]
	movaps xmm4, [esp + nb114_fjyH2]
	movaps xmm5, [esp + nb114_fjzH2]
	mulps xmm0, [esp + nb114_dxH2H2]
	mulps xmm1, [esp + nb114_dyH2H2]
	mulps xmm2, [esp + nb114_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixH2]
	addps xmm1, [esp + nb114_fiyH2]
	addps xmm2, [esp + nb114_fizH2]
	movaps [esp + nb114_fjxH2], xmm3
	movaps [esp + nb114_fjyH2], xmm4
	movaps [esp + nb114_fjzH2], xmm5
	movaps [esp + nb114_fixH2], xmm0
	movaps [esp + nb114_fiyH2], xmm1
	movaps [esp + nb114_fizH2], xmm2

	;# H2-M interaction 
	movaps xmm0, [esp + nb114_rinvH2M]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqMH]
	mulps xmm0, xmm1	;# fs H2-M  
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + nb114_fjxM]
	movaps xmm4, [esp + nb114_fjyM]
	movaps xmm5, [esp + nb114_fjzM]
	mulps xmm0, [esp + nb114_dxH2M]
	mulps xmm1, [esp + nb114_dyH2M]
	mulps xmm2, [esp + nb114_dzH2M]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixH2]
	addps xmm1, [esp + nb114_fiyH2]
	addps xmm2, [esp + nb114_fizH2]
	movaps [esp + nb114_fjxM], xmm3
	movaps [esp + nb114_fjyM], xmm4
	movaps [esp + nb114_fjzM], xmm5
	movaps [esp + nb114_fixH2], xmm0
	movaps [esp + nb114_fiyH2], xmm1
	movaps [esp + nb114_fizH2], xmm2

	;# M-H1 interaction 
	movaps xmm0, [esp + nb114_rinvMH1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqMH]
	mulps xmm0, xmm1	;# fs M-H1 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + nb114_fjxH1]
	movaps xmm4, [esp + nb114_fjyH1]
	movaps xmm5, [esp + nb114_fjzH1]
	mulps xmm0, [esp + nb114_dxMH1]
	mulps xmm1, [esp + nb114_dyMH1]
	mulps xmm2, [esp + nb114_dzMH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixM]
	addps xmm1, [esp + nb114_fiyM]
	addps xmm2, [esp + nb114_fizM]
	movaps [esp + nb114_fjxH1], xmm3
	movaps [esp + nb114_fjyH1], xmm4
	movaps [esp + nb114_fjzH1], xmm5
	movaps [esp + nb114_fixM], xmm0
	movaps [esp + nb114_fiyM], xmm1
	movaps [esp + nb114_fizM], xmm2
	
	;# M-H2 interaction 
	movaps xmm0, [esp + nb114_rinvMH2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqMH]
	mulps xmm0, xmm1	;# fs M-H2 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + nb114_fjxH2]
	movaps xmm4, [esp + nb114_fjyH2]
	movaps xmm5, [esp + nb114_fjzH2]
	mulps xmm0, [esp + nb114_dxMH2]
	mulps xmm1, [esp + nb114_dyMH2]
	mulps xmm2, [esp + nb114_dzMH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixM]
	addps xmm1, [esp + nb114_fiyM]
	addps xmm2, [esp + nb114_fizM]
	movaps [esp + nb114_fjxH2], xmm3
	movaps [esp + nb114_fjyH2], xmm4
	movaps [esp + nb114_fjzH2], xmm5
	movaps [esp + nb114_fixM], xmm0
	movaps [esp + nb114_fiyM], xmm1
	movaps [esp + nb114_fizM], xmm2

	;# M-M interaction 
	movaps xmm0, [esp + nb114_rinvMM]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + nb114_qqMM]
	mulps xmm0, xmm1	;# fs M-M 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps [esp + nb114_vctot], xmm7
	movaps xmm2, xmm0
	movaps xmm3, [esp + nb114_fjxM]
	movaps xmm4, [esp + nb114_fjyM]
	movaps xmm5, [esp + nb114_fjzM]
	mulps xmm0, [esp + nb114_dxMM]
	mulps xmm1, [esp + nb114_dyMM]
	mulps xmm2, [esp + nb114_dzMM]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb114_fixM]
	addps xmm1, [esp + nb114_fiyM]
	addps xmm2, [esp + nb114_fizM]
	movaps [esp + nb114_fjxM], xmm3
	movaps [esp + nb114_fjyM], xmm4
	movaps [esp + nb114_fjzM], xmm5
	movaps [esp + nb114_fixM], xmm0
	movaps [esp + nb114_fiyM], xmm1
	movaps [esp + nb114_fizM], xmm2

	;# Did all interactions - now update j forces 
	mov edi, [ebp + nb114_faction]
	;# 4 j waters with four atoms each.
	;# step 1 : transpose fjxO, fjyO, fjzO, fjxH1
	movaps xmm0, [esp + nb114_fjxO]
	movaps xmm1, [esp + nb114_fjyO]
	movaps xmm2, [esp + nb114_fjzO]
	movaps xmm3, [esp + nb114_fjxH1]
	movaps xmm4, xmm0
	movaps xmm5, xmm1
	unpcklps xmm4, xmm2 	
	unpcklps xmm5, xmm3	
	unpckhps xmm0, xmm2 
	unpckhps xmm1, xmm3 
	movaps xmm2, xmm4
	movaps xmm3, xmm0

	unpcklps xmm4, xmm5 
	unpckhps xmm2, xmm5 
	unpcklps xmm0, xmm1 
	unpckhps xmm3, xmm1
	
	;# results are now in xmm4, xmm2, xmm0, xmm3
	;# load the corresponding j forces from memory
	movlps   xmm1, [edi + eax*4]
	movlps   xmm5, [edi + ebx*4]
	movlps   xmm6, [edi + ecx*4]
	movlps   xmm7, [edi + edx*4]
	movhps   xmm1, [edi + eax*4 + 8]
	movhps   xmm5, [edi + ebx*4 + 8]
	movhps   xmm6, [edi + ecx*4 + 8]
	movhps   xmm7, [edi + edx*4 + 8]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [edi + eax*4], xmm1
	movlps   [edi + ebx*4], xmm5
	movlps   [edi + ecx*4], xmm6
	movlps   [edi + edx*4], xmm7
	movhps   [edi + eax*4 + 8], xmm1
	movhps   [edi + ebx*4 + 8], xmm5
	movhps   [edi + ecx*4 + 8], xmm6
	movhps   [edi + edx*4 + 8], xmm7

	;# step 2 : transpose fjyH1, fjzH1, fjxH2, fjyH2
	movaps xmm0, [esp + nb114_fjyH1]
	movaps xmm1, [esp + nb114_fjzH1]
	movaps xmm2, [esp + nb114_fjxH2]
	movaps xmm3, [esp + nb114_fjyH2]
	movaps xmm4, xmm0
	movaps xmm5, xmm1
	unpcklps xmm4, xmm2 	
	unpcklps xmm5, xmm3	
	unpckhps xmm0, xmm2 
	unpckhps xmm1, xmm3 
	movaps xmm2, xmm4
	movaps xmm3, xmm0
	unpcklps xmm4, xmm5 
	unpckhps xmm2, xmm5 

	unpcklps xmm0, xmm1 
	unpckhps xmm3, xmm1 
	;# results are now in xmm4, xmm2, xmm0, xmm3
	;# load the corresponding j forces from memory
	movlps   xmm1, [edi + eax*4 + 16]
	movlps   xmm5, [edi + ebx*4 + 16]
	movlps   xmm6, [edi + ecx*4 + 16]
	movlps   xmm7, [edi + edx*4 + 16]
	movhps   xmm1, [edi + eax*4 + 24]
	movhps   xmm5, [edi + ebx*4 + 24]
	movhps   xmm6, [edi + ecx*4 + 24]
	movhps   xmm7, [edi + edx*4 + 24]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [edi + eax*4 + 16], xmm1
	movlps   [edi + ebx*4 + 16], xmm5
	movlps   [edi + ecx*4 + 16], xmm6
	movlps   [edi + edx*4 + 16], xmm7
	movhps   [edi + eax*4 + 24], xmm1
	movhps   [edi + ebx*4 + 24], xmm5
	movhps   [edi + ecx*4 + 24], xmm6
	movhps   [edi + edx*4 + 24], xmm7

	;# step 3 : transpose fjzH2, fjxM, fjyM, fjzM. xmm4 is scratch
	movaps xmm0, [esp + nb114_fjzH2]
	movaps xmm1, [esp + nb114_fjxM]
	movaps xmm2, [esp + nb114_fjyM]
	movaps xmm3, [esp + nb114_fjzM]

        movaps xmm4, xmm0
        movaps xmm5, xmm1
        unpcklps xmm4, xmm2
        unpcklps xmm5, xmm3
        unpckhps xmm0, xmm2
        unpckhps xmm1, xmm3
        movaps xmm2, xmm4
        movaps xmm3, xmm0
        unpcklps xmm4, xmm5
        unpckhps xmm2, xmm5
	
        unpcklps xmm0, xmm1
        unpckhps xmm3, xmm1

	;# results are now in xmm0, xmm1, xmm2, xmm3
	;# load the corresponding j forces from memory
	movlps   xmm1, [edi + eax*4 + 32]
	movlps   xmm5, [edi + ebx*4 + 32]
	movlps   xmm6, [edi + ecx*4 + 32]
	movlps   xmm7, [edi + edx*4 + 32]
	movhps   xmm1, [edi + eax*4 + 40]
	movhps   xmm5, [edi + ebx*4 + 40]
	movhps   xmm6, [edi + ecx*4 + 40]
	movhps   xmm7, [edi + edx*4 + 40]
	;# add
	addps    xmm1, xmm4
	addps    xmm5, xmm2
	addps    xmm6, xmm0
	addps    xmm7, xmm3
	;# store back
	movlps   [edi + eax*4 + 32], xmm1
	movlps   [edi + ebx*4 + 32], xmm5
	movlps   [edi + ecx*4 + 32], xmm6
	movlps   [edi + edx*4 + 32], xmm7
	movhps   [edi + eax*4 + 40], xmm1
	movhps   [edi + ebx*4 + 40], xmm5
	movhps   [edi + ecx*4 + 40], xmm6
	movhps   [edi + edx*4 + 40], xmm7

	;# should we do one more iteration? 
	sub dword ptr [esp + nb114_innerk],  4
	jl    .nb114_single_check
	jmp   .nb114_unroll_loop
.nb114_single_check:
	add dword ptr [esp + nb114_innerk],  4
	jnz   .nb114_single_loop
	jmp   .nb114_updateouterdata
.nb114_single_loop:
	mov   edx, [esp + nb114_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb114_innerjjnr],  4	

	mov esi, [ebp + nb114_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates	
	movlps xmm3,  [esi + eax*4]		;#  Ox  Oy  
	movlps xmm4,  [esi + eax*4 + 16]	;# H1y H1z 
	movlps xmm5,  [esi + eax*4 + 32]	;# H2z  Mx 
	movhps xmm3,  [esi + eax*4 + 8]   	;#  Ox  Oy  Oz H1x
	movhps xmm4,  [esi + eax*4 + 24]	;# H1y H1z H2x H2y
	movhps xmm5,  [esi + eax*4 + 40]	;# H2z  Mx  My  Mz
	;# transpose
	movaps xmm0, xmm4
	movaps xmm1, xmm3
	movaps xmm2, xmm4
	movaps xmm6, xmm3
	shufps xmm4, xmm5, 18  ;# (00010010)  h2x - Mx  - 
	shufps xmm3, xmm0, 193 ;# (11000001)  Oy  - H1y - 
	shufps xmm2, xmm5, 35  ;# (00100011) H2y - My  - 
 	shufps xmm1, xmm0, 18  ;# (00010010)  Oz  - H1z - 
	;#  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
	shufps xmm6, xmm4, 140 ;# (10001100) Ox H1x H2x Mx 
	shufps xmm3, xmm2, 136 ;# (10001000) Oy H1y H2y My 
	shufps xmm1, xmm5, 200 ;# (11001000) Oz H1z H2z Mz

	;# store all j coordinates in jO  
	movaps [esp + nb114_jxO], xmm6
	movaps [esp + nb114_jyO], xmm3
	movaps [esp + nb114_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [esp + nb114_ixO]
	movaps xmm1, [esp + nb114_iyO]
	movaps xmm2, [esp + nb114_izO]
	movaps xmm3, [esp + nb114_ixM]
	movaps xmm4, [esp + nb114_iyM]
	movaps xmm5, [esp + nb114_izM]
	subps  xmm0, [esp + nb114_jxO]
	subps  xmm1, [esp + nb114_jyO]
	subps  xmm2, [esp + nb114_jzO]
	subps  xmm3, [esp + nb114_jxO]
	subps  xmm4, [esp + nb114_jyO]
	subps  xmm5, [esp + nb114_jzO]
	
	movaps [esp + nb114_dxOO], xmm0
	movaps [esp + nb114_dyOO], xmm1
	movaps [esp + nb114_dzOO], xmm2
	movaps [esp + nb114_dxMM], xmm3
	movaps [esp + nb114_dyMM], xmm4
	movaps [esp + nb114_dzMM], xmm5
	
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm4, xmm3
	addps xmm4, xmm5	;# have rsq in xmm4
	;# Save M data 
	movaps [esp + nb114_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [esp + nb114_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [esp + nb114_three]
	mulss  xmm2, xmm1 	;# constant 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# constant 1/r6
	mulps   xmm7, [esp + nb114_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [esp + nb114_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# constant 1/r12
	mulss  xmm1, [esp + nb114_c6]
	mulss  xmm2, [esp + nb114_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [esp + nb114_Vvdwtot]
	movss  [esp + nb114_Vvdwtot], xmm3
	mulss  xmm1, [esp + nb114_six]
	mulss  xmm2, [esp + nb114_twelve]
	subss  xmm2, xmm1
	mulss  xmm0, xmm2 	;# fscal
	movss  xmm1, xmm0
	movss  xmm2, xmm0
	mulss  xmm0, [esp + nb114_dxOO]
	mulss  xmm1, [esp + nb114_dyOO]
	mulss  xmm2, [esp + nb114_dzOO]
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movaps  [esp + nb114_fjxO], xmm3
	movaps  [esp + nb114_fjyO], xmm4
	movaps  [esp + nb114_fjzO], xmm5
	addss   xmm0, [esp + nb114_fixO]
	addss   xmm1, [esp + nb114_fiyO]
	addss   xmm2, [esp + nb114_fizO]
	movss  [esp + nb114_fixO], xmm0
	movss  [esp + nb114_fiyO], xmm1
	movss  [esp + nb114_fizO], xmm2

	;# do  M coulomb interaction
	movaps xmm0, [esp + nb114_rinvMM]
	movaps xmm4, xmm0	;# xmm4=rinv
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb114_qqMH]
	movhps  xmm3, [esp + nb114_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm4, xmm3 	;# xmm4=voul	
	mulps  xmm0, xmm4
	addps  xmm4, [esp + nb114_vctot] 
	movaps [esp + nb114_vctot], xmm4

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + nb114_dxMM]
	mulps   xmm1, [esp + nb114_dyMM]
	mulps   xmm2, [esp + nb114_dzMM]
	;# update forces M - j water 
	movaps  xmm3, [esp + nb114_fjxO]
	movaps  xmm4, [esp + nb114_fjyO]
	movaps  xmm5, [esp + nb114_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb114_fjxO], xmm3
	movaps  [esp + nb114_fjyO], xmm4
	movaps  [esp + nb114_fjzO], xmm5
	addps   xmm0, [esp + nb114_fixM]
	addps   xmm1, [esp + nb114_fiyM]
	addps   xmm2, [esp + nb114_fizM]
	movaps  [esp + nb114_fixM], xmm0
	movaps  [esp + nb114_fiyM], xmm1
	movaps  [esp + nb114_fizM], xmm2	
	
	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb114_ixH1]
	movaps  xmm1, [esp + nb114_iyH1]
	movaps  xmm2, [esp + nb114_izH1]	
	movaps  xmm3, [esp + nb114_ixH2] 
	movaps  xmm4, [esp + nb114_iyH2] 
	movaps  xmm5, [esp + nb114_izH2] 
	subps   xmm0, [esp + nb114_jxO]
	subps   xmm1, [esp + nb114_jyO]
	subps   xmm2, [esp + nb114_jzO]
	subps   xmm3, [esp + nb114_jxO]
	subps   xmm4, [esp + nb114_jyO]
	subps   xmm5, [esp + nb114_jzO]
	movaps [esp + nb114_dxH1H1], xmm0
	movaps [esp + nb114_dyH1H1], xmm1
	movaps [esp + nb114_dzH1H1], xmm2
	movaps [esp + nb114_dxH2H2], xmm3
	movaps [esp + nb114_dyH2H2], xmm4
	movaps [esp + nb114_dzH2H2], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	movaps  [esp + nb114_rsqH1H1], xmm0
	movaps  [esp + nb114_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114_half] ;# rinvH1H1
	mulps   xmm7, [esp + nb114_half] ;# rinvH2H2
	movaps  [esp + nb114_rinvH1H1], xmm3
	movaps  [esp + nb114_rinvH2H2], xmm7
	
	;# Do H1 coulomb interaction
	movaps xmm0, [esp + nb114_rinvH1H1]
	movaps xmm4, xmm0	;# xmm4=rinv 
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb114_qqHH]
	movhps  xmm3, [esp + nb114_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm4, xmm3 	;# xmm4=voul
	mulps  xmm0, xmm4
	addps  xmm4, [esp + nb114_vctot] 
	movaps [esp + nb114_vctot], xmm4


	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + nb114_dxH1H1]
	mulps   xmm1, [esp + nb114_dyH1H1]
	mulps   xmm2, [esp + nb114_dzH1H1]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + nb114_fjxO]
	movaps  xmm4, [esp + nb114_fjyO]
	movaps  xmm5, [esp + nb114_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb114_fjxO], xmm3
	movaps  [esp + nb114_fjyO], xmm4
	movaps  [esp + nb114_fjzO], xmm5
	addps   xmm0, [esp + nb114_fixH1]
	addps   xmm1, [esp + nb114_fiyH1]
	addps   xmm2, [esp + nb114_fizH1]
	movaps  [esp + nb114_fixH1], xmm0
	movaps  [esp + nb114_fiyH1], xmm1
	movaps  [esp + nb114_fizH1], xmm2	

	;# H2 Coulomb
	movaps xmm0, [esp + nb114_rinvH2H2]
	movaps xmm4, xmm0	;# xmm4=rinv 
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb114_qqHH]
	movhps  xmm3, [esp + nb114_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm4, xmm3       ;# xmm4=voul
	mulps  xmm0, xmm4
	addps  xmm4, [esp + nb114_vctot] ;# local vctot summation variable
	movaps [esp + nb114_vctot], xmm4

	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + nb114_dxH2H2]
	mulps   xmm1, [esp + nb114_dyH2H2]
	mulps   xmm2, [esp + nb114_dzH2H2]
	;# update forces H2 - j water 
	movaps  xmm3, [esp + nb114_fjxO]
	movaps  xmm4, [esp + nb114_fjyO]
	movaps  xmm5, [esp + nb114_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb114_fjxO], xmm3
	movaps  [esp + nb114_fjyO], xmm4
	movaps  [esp + nb114_fjzO], xmm5
	addps   xmm0, [esp + nb114_fixH2]
	addps   xmm1, [esp + nb114_fiyH2]
	addps   xmm2, [esp + nb114_fizH2]
	movaps  [esp + nb114_fixH2], xmm0
	movaps  [esp + nb114_fiyH2], xmm1
	movaps  [esp + nb114_fizH2], xmm2			
	
	mov     esi, [ebp + nb114_faction]
	;# update j water forces from local variables.
	;# transpose back first
	movaps  xmm0, [esp + nb114_fjxO] ;# Ox H1x H2x Mx 
	movaps  xmm1, [esp + nb114_fjyO] ;# Oy H1y H2y My
	movaps  xmm2, [esp + nb114_fjzO] ;# Oz H1z H2z Mz
	 
	movaps  xmm3, xmm0
	movaps  xmm4, xmm0
	unpcklps xmm3, xmm1       	;# Ox Oy - -
	shufps  xmm4, xmm2, 0x1	  	;# h1x - Oz -
	movaps  xmm5, xmm1
	movaps  xmm6, xmm0
	unpcklps xmm5, xmm2 	  	;# - - H1y H1z
	unpckhps xmm6, xmm1	  	;# h2x h2y - - 
	unpckhps xmm1, xmm2	  	;# - - My Mz
	
	shufps   xmm2, xmm0, 0x32  ;# (00110010) h2z - Mx -
	shufps   xmm3, xmm4, 36  ;# constant 00100100 ;# Ox Oy Oz H1x 
	shufps   xmm5, xmm6, 78  ;# constant 01001110 ;# h1y h1z h2x h2y
	shufps   xmm2, xmm1, 232  ;# constant 11101000 ;# h2z mx my mz

	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 16]
	movlps  xmm4, [esi + eax*4 + 32]
	movhps  xmm0, [esi + eax*4 + 8]
	movhps  xmm1, [esi + eax*4 + 24]
	movhps  xmm4, [esi + eax*4 + 40]
	addps   xmm0, xmm3
	addps   xmm1, xmm5
	addps   xmm4, xmm2
	movlps   [esi + eax*4], xmm0
	movlps   [esi + eax*4 + 16], xmm1
	movlps   [esi + eax*4 + 32], xmm4
	movhps   [esi + eax*4 + 8], xmm0
	movhps   [esi + eax*4 + 24], xmm1
	movhps   [esi + eax*4 + 40], xmm4

	dec dword ptr [esp + nb114_innerk]
	jz    .nb114_updateouterdata
	jmp   .nb114_single_loop
.nb114_updateouterdata:
	mov   ecx, [esp + nb114_ii3]
	mov   edi, [ebp + nb114_faction]
	mov   esi, [ebp + nb114_fshift]
	mov   edx, [esp + nb114_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb114_fixO]
	movaps xmm1, [esp + nb114_fiyO] 
	movaps xmm2, [esp + nb114_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb114_fixH1]
	movaps xmm1, [esp + nb114_fiyH1]
	movaps xmm2, [esp + nb114_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb114_fixH2]
	movaps xmm1, [esp + nb114_fiyH2]
	movaps xmm2, [esp + nb114_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb114_fixM]
	movaps xmm1, [esp + nb114_fiyM]
	movaps xmm2, [esp + nb114_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 36]
	movss  xmm4, [edi + ecx*4 + 40]
	movss  xmm5, [edi + ecx*4 + 44]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 36], xmm3
	movss  [edi + ecx*4 + 40], xmm4
	movss  [edi + ecx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb114_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb114_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb114_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb114_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb114_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb114_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
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




	
.globl nb_kernel114nf_ia32_sse 
.globl _nb_kernel114nf_ia32_sse
nb_kernel114nf_ia32_sse:	
_nb_kernel114nf_ia32_sse:	
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
.equiv          nb114nf_p_krf,          56
.equiv          nb114nf_p_crf,          60
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
	;# bottom of stack is cache-aligned for sse use 
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
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb114nf_half], eax
	movss xmm1, [esp + nb114nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb114nf_half],  xmm1
	movaps [esp + nb114nf_two],  xmm2
	movaps [esp + nb114nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb114nf_iinr]   ;# ecx = pointer into iinr[]
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb114nf_charge]
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	movss xmm4, xmm3	
	mov esi, [ebp + nb114nf_p_facel]
	movss xmm6, [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb114nf_qqMM], xmm3
	movaps [esp + nb114nf_qqMH], xmm4
	movaps [esp + nb114nf_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + nb114nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb114nf_p_ntype]
	imul  ecx, [edi]  ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb114nf_vdwparam]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0x55
	movaps [esp + nb114nf_c6], xmm0
	movaps [esp + nb114nf_c12], xmm1

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
	mov   eax, [ebp + nb114nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb114nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb114nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb114nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb114nf_pos]	;# eax = base of pos[]  
	mov   [esp + nb114nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	movaps xmm6, xmm0
	movaps xmm7, xmm1

	addss xmm3, [eax + ebx*4]  	;# ox
	addss xmm4, [eax + ebx*4 + 4]  ;# oy
	addss xmm5, [eax + ebx*4 + 8]  ;# oz
	addss xmm6, [eax + ebx*4 + 12] ;# h1x
	addss xmm7, [eax + ebx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [esp + nb114nf_ixO], xmm3
	movaps [esp + nb114nf_iyO], xmm4
	movaps [esp + nb114nf_izO], xmm5
	movaps [esp + nb114nf_ixH1], xmm6
	movaps [esp + nb114nf_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [eax + ebx*4 + 20] ;# h1z
	addss xmm0, [eax + ebx*4 + 24] ;# h2x
	addss xmm1, [eax + ebx*4 + 28] ;# h2y
	addss xmm2, [eax + ebx*4 + 32] ;# h2z
	addss xmm3, [eax + ebx*4 + 36] ;# mx
	addss xmm4, [eax + ebx*4 + 40] ;# my
	addss xmm5, [eax + ebx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb114nf_izH1], xmm6
	movaps [esp + nb114nf_ixH2], xmm0
	movaps [esp + nb114nf_iyH2], xmm1
	movaps [esp + nb114nf_izH2], xmm2
	movaps [esp + nb114nf_ixM], xmm3
	movaps [esp + nb114nf_iyM], xmm4
	movaps [esp + nb114nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [esp + nb114nf_vctot], xmm4
	movaps [esp + nb114nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb114nf_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb114nf_pos]
	mov   edi, [ebp + nb114nf_faction]	
	mov   eax, [ebp + nb114nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb114nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb114nf_ninner]
	mov   [esp + nb114nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb114nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb114nf_unroll_loop
	jmp   .nb114nf_single_check
.nb114nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb114nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb114nf_innerjjnr],  16 ;# advance pointer (unroll 4) 

	mov esi, [ebp + nb114nf_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables
	;# Load Ox, Oy, Oz, H1x 
	movlps xmm1, [esi + eax*4]	;#  Oxa   Oya    -    -
	movlps xmm4, [esi + ecx*4]	;#  Oxc   Oyc    -    -
	movhps xmm1, [esi + ebx*4]	;#  Oxa   Oya   Oxb   Oyb 
	movhps xmm4, [esi + edx*4]	;#  Oxc   Oyc   Oxd   Oyd 
	movaps xmm0, xmm1		;#  Oxa   Oya   Oxb   Oyb 
	shufps xmm0, xmm4, 0x88		;#  Oxa   Oxb   Oxc   Oxd
	shufps xmm1, xmm4, 0xDD		;#  Oya   Oyb   Oyc   Oyd
	movlps xmm3, [esi + eax*4 + 8]	;#  Oza  H1xa    -    -
	movlps xmm5, [esi + ecx*4 + 8]	;#  Ozc  H1xc    -    -
	movhps xmm3, [esi + ebx*4 + 8]	;#  Oza  H1xa   Ozb  H1xb 
	movhps xmm5, [esi + edx*4 + 8]	;#  Ozc  H1xc   Ozd  H1xd 
	movaps xmm2, xmm3		;#  Oza  H1xa   Ozb  H1xb 
	shufps xmm2, xmm5, 0x88		;#  Oza   Ozb   Ozc   Ozd
	shufps xmm3, xmm5, 0xDD		;# H1xa  H1xb  H1xc  H1xd
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb114nf_jxO], xmm0
	movaps [esp + nb114nf_jyO], xmm1
	movaps [esp + nb114nf_jzO], xmm2
	movaps [esp + nb114nf_jxH1], xmm3

	;# Load H1y H1z H2x H2y 
	movlps xmm1, [esi + eax*4 + 16]	
	movlps xmm4, [esi + ecx*4 + 16]	
	movhps xmm1, [esi + ebx*4 + 16]	
	movhps xmm4, [esi + edx*4 + 16]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 24]	
	movlps xmm5, [esi + ecx*4 + 24]	
	movhps xmm3, [esi + ebx*4 + 24]	
	movhps xmm5, [esi + edx*4 + 24]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb114nf_jyH1], xmm0
	movaps [esp + nb114nf_jzH1], xmm1
	movaps [esp + nb114nf_jxH2], xmm2
	movaps [esp + nb114nf_jyH2], xmm3

	;# Load H2z Mx My Mz 
	movlps xmm1, [esi + eax*4 + 32]	
	movlps xmm4, [esi + ecx*4 + 32]	
	movhps xmm1, [esi + ebx*4 + 32]	
	movhps xmm4, [esi + edx*4 + 32]	
	movaps xmm0, xmm1		
	shufps xmm0, xmm4, 0x88		
	shufps xmm1, xmm4, 0xDD		
	movlps xmm3, [esi + eax*4 + 40]	
	movlps xmm5, [esi + ecx*4 + 40]	
	movhps xmm3, [esi + ebx*4 + 40]	
	movhps xmm5, [esi + edx*4 + 40]	
	movaps xmm2, xmm3		
	shufps xmm2, xmm5, 0x88		
	shufps xmm3, xmm5, 0xDD		
	;# coordinates in xmm0-xmm3	
	;# store
	movaps [esp + nb114nf_jzH2], xmm0
	movaps [esp + nb114nf_jxM], xmm1
	movaps [esp + nb114nf_jyM], xmm2
	movaps [esp + nb114nf_jzM], xmm3
	
	;# start calculating pairwise distances
	movaps xmm0, [esp + nb114nf_ixO]
	movaps xmm1, [esp + nb114nf_iyO]
	movaps xmm2, [esp + nb114nf_izO]
	movaps xmm3, [esp + nb114nf_ixH1]
	movaps xmm4, [esp + nb114nf_iyH1]
	movaps xmm5, [esp + nb114nf_izH1]
	subps  xmm0, [esp + nb114nf_jxO]
	subps  xmm1, [esp + nb114nf_jyO]
	subps  xmm2, [esp + nb114nf_jzO]
	subps  xmm3, [esp + nb114nf_jxH1]
	subps  xmm4, [esp + nb114nf_jyH1]
	subps  xmm5, [esp + nb114nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb114nf_rsqOO], xmm0
	movaps [esp + nb114nf_rsqH1H1], xmm3

	movaps xmm0, [esp + nb114nf_ixH1]
	movaps xmm1, [esp + nb114nf_iyH1]
	movaps xmm2, [esp + nb114nf_izH1]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb114nf_jxH2]
	subps  xmm1, [esp + nb114nf_jyH2]
	subps  xmm2, [esp + nb114nf_jzH2]
	subps  xmm3, [esp + nb114nf_jxM]
	subps  xmm4, [esp + nb114nf_jyM]
	subps  xmm5, [esp + nb114nf_jzM]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb114nf_rsqH1H2], xmm0
	movaps [esp + nb114nf_rsqH1M], xmm3

	movaps xmm0, [esp + nb114nf_ixH2]
	movaps xmm1, [esp + nb114nf_iyH2]
	movaps xmm2, [esp + nb114nf_izH2]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb114nf_jxH1]
	subps  xmm1, [esp + nb114nf_jyH1]
	subps  xmm2, [esp + nb114nf_jzH1]
	subps  xmm3, [esp + nb114nf_jxH2]
	subps  xmm4, [esp + nb114nf_jyH2]
	subps  xmm5, [esp + nb114nf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb114nf_rsqH2H1], xmm0
	movaps [esp + nb114nf_rsqH2H2], xmm3

	movaps xmm0, [esp + nb114nf_ixH2]
	movaps xmm1, [esp + nb114nf_iyH2]
	movaps xmm2, [esp + nb114nf_izH2]
	movaps xmm3, [esp + nb114nf_ixM]
	movaps xmm4, [esp + nb114nf_iyM]
	movaps xmm5, [esp + nb114nf_izM]
	subps  xmm0, [esp + nb114nf_jxM]
	subps  xmm1, [esp + nb114nf_jyM]
	subps  xmm2, [esp + nb114nf_jzM]
	subps  xmm3, [esp + nb114nf_jxH1]
	subps  xmm4, [esp + nb114nf_jyH1]
	subps  xmm5, [esp + nb114nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb114nf_rsqH2M], xmm0
	movaps [esp + nb114nf_rsqMH1], xmm4

	movaps xmm0, [esp + nb114nf_ixM]
	movaps xmm1, [esp + nb114nf_iyM]
	movaps xmm2, [esp + nb114nf_izM]
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	subps  xmm0, [esp + nb114nf_jxH2]
	subps  xmm1, [esp + nb114nf_jyH2]
	subps  xmm2, [esp + nb114nf_jzH2]
	subps  xmm3, [esp + nb114nf_jxM]
	subps  xmm4, [esp + nb114nf_jyM]
	subps  xmm5, [esp + nb114nf_jzM]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb114nf_rsqMH2], xmm0
	movaps [esp + nb114nf_rsqMM], xmm4

	;# start by doing reciprocal for OO
	movaps  xmm7, [esp + nb114nf_rsqOO]
	rcpps   xmm2, xmm7
	movaps  xmm1, [esp + nb114nf_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps [esp + nb114nf_rinvsqOO], xmm2
	
	;# next step is invsqrt - do two at a time.
	rsqrtps xmm1, [esp + nb114nf_rsqH1H1]
	rsqrtps xmm5, [esp + nb114nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114nf_rsqH1H1]
	mulps   xmm5, [esp + nb114nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114nf_half] ;# rinvH1H1 
	mulps   xmm7, [esp + nb114nf_half] ;# rinvH1H2 
	movaps  [esp + nb114nf_rinvH1H1], xmm3
	movaps  [esp + nb114nf_rinvH1H2], xmm7
			
	rsqrtps xmm1, [esp + nb114nf_rsqH1M]
	rsqrtps xmm5, [esp + nb114nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114nf_rsqH1M]
	mulps   xmm5, [esp + nb114nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114nf_half] 
	mulps   xmm7, [esp + nb114nf_half]
	movaps  [esp + nb114nf_rinvH1M], xmm3
	movaps  [esp + nb114nf_rinvH2H1], xmm7
			
	rsqrtps xmm1, [esp + nb114nf_rsqH2H2]
	rsqrtps xmm5, [esp + nb114nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114nf_rsqH2H2]
	mulps   xmm5, [esp + nb114nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114nf_half] 
	mulps   xmm7, [esp + nb114nf_half]
	movaps  [esp + nb114nf_rinvH2H2], xmm3
	movaps  [esp + nb114nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [esp + nb114nf_rsqMH1]
	rsqrtps xmm5, [esp + nb114nf_rsqMH2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb114nf_rsqMH1]
	mulps   xmm5, [esp + nb114nf_rsqMH2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114nf_half] 
	mulps   xmm7, [esp + nb114nf_half]
	movaps  [esp + nb114nf_rinvMH1], xmm3
	movaps  [esp + nb114nf_rinvMH2], xmm7
        		
	rsqrtps xmm1, [esp + nb114nf_rsqMM]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb114nf_three]
	mulps   xmm1, [esp + nb114nf_rsqMM]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb114nf_half] 
	movaps  [esp + nb114nf_rinvMM], xmm3
	
	;# start with OO LJ interaction
	movaps xmm0, [esp + nb114nf_rinvsqOO]
	movaps xmm1, xmm0
	mulps  xmm1, xmm1	;# rinv4
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb114nf_c6]
	mulps  xmm2, [esp + nb114nf_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [esp + nb114nf_Vvdwtot]
	movaps [esp + nb114nf_Vvdwtot], xmm4

	;# Coulomb interactions 
	;# all H-H interactions
	movaps xmm0, [esp + nb114nf_rinvH1H1]
	addps  xmm0, [esp + nb114nf_rinvH1H2]
	addps  xmm0, [esp + nb114nf_rinvH2H1]
	addps  xmm0, [esp + nb114nf_rinvH2H2]
	mulps  xmm0, [esp + nb114nf_qqHH]
	;# all M-H interactions
	movaps xmm1, [esp + nb114nf_rinvH1M]
	addps  xmm1, [esp + nb114nf_rinvH2M]
	addps  xmm1, [esp + nb114nf_rinvMH1]
	addps  xmm1, [esp + nb114nf_rinvMH2]
	mulps  xmm1, [esp + nb114nf_qqMH]
	;# The M-M interaction
	movaps xmm2, [esp + nb114nf_rinvMM]
	mulps  xmm2, [esp + nb114nf_qqMM]
	addps  xmm0, xmm1
	addps  xmm2, [esp + nb114nf_vctot] 
	addps  xmm0, xmm2
	movaps [esp + nb114nf_vctot], xmm0 
	;# should we do one more iteration? 
	sub dword ptr [esp + nb114nf_innerk],  4
	jl    .nb114nf_single_check
	jmp   .nb114nf_unroll_loop
.nb114nf_single_check:
	add dword ptr [esp + nb114nf_innerk],  4
	jnz   .nb114nf_single_loop
	jmp   .nb114nf_updateouterdata
.nb114nf_single_loop:
	mov   edx, [esp + nb114nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb114nf_innerjjnr],  4	

	mov esi, [ebp + nb114nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates
	movlps xmm3,  [esi + eax*4]		;#  Ox  Oy  
	movlps xmm4,  [esi + eax*4 + 16]	;# H1y H1z 
	movlps xmm5,  [esi + eax*4 + 32]	;# H2z  Mx 
	movhps xmm3,  [esi + eax*4 + 8]   	;#  Ox  Oy  Oz H1x
	movhps xmm4,  [esi + eax*4 + 24]	;# H1y H1z H2x H2y
	movhps xmm5,  [esi + eax*4 + 40]	;# H2z  Mx  My  Mz
	;# transpose
	movaps xmm0, xmm4
	movaps xmm1, xmm3
	movaps xmm2, xmm4
	movaps xmm6, xmm3
	shufps xmm4, xmm5, 18  ;# (00010010)  h2x - Mx  - 
	shufps xmm3, xmm0, 193 ;# (11000001)  Oy  - H1y - 
	shufps xmm2, xmm5, 35  ;# (00100011) H2y - My  - 
 	shufps xmm1, xmm0, 18  ;# (00010010)  Oz  - H1z - 
	;#  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
	shufps xmm6, xmm4, 140 ;# (10001100) Ox H1x H2x Mx 
	shufps xmm3, xmm2, 136 ;# (10001000) Oy H1y H2y My 
	shufps xmm1, xmm5, 200 ;# (11001000) Oz H1z H2z Mz

	;# store all j coordinates in jO  
	movaps [esp + nb114nf_jxO], xmm6
	movaps [esp + nb114nf_jyO], xmm3
	movaps [esp + nb114nf_jzO], xmm1

	;# do O and M in parallel
	movaps xmm0, [esp + nb114nf_ixO]
	movaps xmm1, [esp + nb114nf_iyO]
	movaps xmm2, [esp + nb114nf_izO]
	movaps xmm3, [esp + nb114nf_ixM]
	movaps xmm4, [esp + nb114nf_iyM]
	movaps xmm5, [esp + nb114nf_izM]
	subps  xmm0, [esp + nb114nf_jxO]
	subps  xmm1, [esp + nb114nf_jyO]
	subps  xmm2, [esp + nb114nf_jzO]
	subps  xmm3, [esp + nb114nf_jxO]
	subps  xmm4, [esp + nb114nf_jyO]
	subps  xmm5, [esp + nb114nf_jzO]
		
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm4, xmm3
	addps xmm4, xmm5	;# have rsq in xmm4
	;# Save M data 
	movaps [esp + nb114nf_rsqMM], xmm4
	
	;# do 1/x for O and 1/sqrt(x) for M
	rcpss  xmm1, xmm0
	rsqrtps xmm5, xmm4
	movss  xmm2, [esp + nb114nf_two]
	movaps  xmm6, xmm5	
	mulss  xmm0, xmm1
	mulps   xmm5, xmm5
	subss  xmm2, xmm0
	movaps  xmm7, [esp + nb114nf_three]
	mulss  xmm2, xmm1 	;# constant 1/r2
	
	mulps   xmm5, xmm4
	movss  xmm0, xmm2
	subps   xmm7, xmm5
	mulss  xmm2, xmm2
	mulps   xmm7, xmm6
	mulss  xmm2, xmm0 	;# constant 1/r6
	mulps   xmm7, [esp + nb114nf_half] ;# rinv iH1 - j water 
	movss  xmm1, xmm2
	movaps [esp + nb114nf_rinvMM], xmm7

	mulss  xmm2, xmm2 	;# constant 1/r12
	mulss  xmm1, [esp + nb114nf_c6]
	mulss  xmm2, [esp + nb114nf_c12]
	movss  xmm3, xmm2
	subss  xmm3, xmm1
	addss  xmm3, [esp + nb114nf_Vvdwtot]
	movss  [esp + nb114nf_Vvdwtot], xmm3

	;# do  M coulomb interaction
	movaps xmm0, [esp + nb114nf_rinvMM]
	movaps xmm4, xmm0	;# xmm4=rinv

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb114nf_qqMH]
	movhps  xmm3, [esp + nb114nf_qqMM]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm4, xmm3 	;# xmm4=voul	
	addps  xmm4, [esp + nb114nf_vctot] 
	movaps [esp + nb114nf_vctot], xmm4

	;# i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb114nf_ixH1]
	movaps  xmm1, [esp + nb114nf_iyH1]
	movaps  xmm2, [esp + nb114nf_izH1]	
	movaps  xmm3, [esp + nb114nf_ixH2] 
	movaps  xmm4, [esp + nb114nf_iyH2] 
	movaps  xmm5, [esp + nb114nf_izH2] 
	subps   xmm0, [esp + nb114nf_jxO]
	subps   xmm1, [esp + nb114nf_jyO]
	subps   xmm2, [esp + nb114nf_jzO]
	subps   xmm3, [esp + nb114nf_jxO]
	subps   xmm4, [esp + nb114nf_jyO]
	subps   xmm5, [esp + nb114nf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	movaps  [esp + nb114nf_rsqH1H1], xmm0
	movaps  [esp + nb114nf_rsqH2H2], xmm4
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb114nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb114nf_half] ;# rinvH1H1
	mulps   xmm7, [esp + nb114nf_half] ;# rinvH2H2
	movaps  [esp + nb114nf_rinvH1H1], xmm3
	movaps  [esp + nb114nf_rinvH2H2], xmm7
	
	;# Do H1 & H2 coulomb interaction
	movaps xmm0, [esp + nb114nf_rinvH1H1]
        addps  xmm0, [esp + nb114nf_rinvH2H2]

	;# fetch charges to xmm3 (temporary) 
	xorps  xmm3, xmm3
	movss   xmm3, [esp + nb114nf_qqHH]
	movhps  xmm3, [esp + nb114nf_qqMH]
	shufps  xmm3, xmm3, 193 ;# constant 11000001 

	mulps  xmm0, xmm3 	;# xmm4=voul
	addps  xmm0, [esp + nb114nf_vctot] 
	movaps [esp + nb114nf_vctot], xmm0

	dec dword ptr [esp + nb114nf_innerk]
	jz    .nb114nf_updateouterdata
	jmp   .nb114nf_single_loop
.nb114nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb114nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb114nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb114nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb114nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb114nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb114nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
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

