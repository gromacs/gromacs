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

	
.globl nb_kernel134_ia32_sse2
.globl _nb_kernel134_ia32_sse2
nb_kernel134_ia32_sse2:	
_nb_kernel134_ia32_sse2:	
.equiv          nb134_p_nri,            8
.equiv          nb134_iinr,             12
.equiv          nb134_jindex,           16
.equiv          nb134_jjnr,             20
.equiv          nb134_shift,            24
.equiv          nb134_shiftvec,         28
.equiv          nb134_fshift,           32
.equiv          nb134_gid,              36
.equiv          nb134_pos,              40
.equiv          nb134_faction,          44
.equiv          nb134_charge,           48
.equiv          nb134_p_facel,          52
.equiv          nb134_argkrf,           56
.equiv          nb134_argcrf,           60
.equiv          nb134_Vc,               64
.equiv          nb134_type,             68
.equiv          nb134_p_ntype,          72
.equiv          nb134_vdwparam,         76
.equiv          nb134_Vvdw,             80
.equiv          nb134_p_tabscale,       84
.equiv          nb134_VFtab,            88
.equiv          nb134_invsqrta,         92
.equiv          nb134_dvda,             96
.equiv          nb134_p_gbtabscale,     100
.equiv          nb134_GBtab,            104
.equiv          nb134_p_nthreads,       108
.equiv          nb134_count,            112
.equiv          nb134_mtx,              116
.equiv          nb134_outeriter,        120
.equiv          nb134_inneriter,        124
.equiv          nb134_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb134_ixO,              0
.equiv          nb134_iyO,              16
.equiv          nb134_izO,              32
.equiv          nb134_ixH1,             48
.equiv          nb134_iyH1,             64
.equiv          nb134_izH1,             80
.equiv          nb134_ixH2,             96
.equiv          nb134_iyH2,             112
.equiv          nb134_izH2,             128
.equiv          nb134_ixM,              144
.equiv          nb134_iyM,              160
.equiv          nb134_izM,              176
.equiv          nb134_jxO,              192
.equiv          nb134_jyO,              208
.equiv          nb134_jzO,              224
.equiv          nb134_jxH1,             240
.equiv          nb134_jyH1,             256
.equiv          nb134_jzH1,             272
.equiv          nb134_jxH2,             288
.equiv          nb134_jyH2,             304
.equiv          nb134_jzH2,             320
.equiv          nb134_jxM,              336
.equiv          nb134_jyM,              352
.equiv          nb134_jzM,              368
.equiv          nb134_dxOO,             384
.equiv          nb134_dyOO,             400
.equiv          nb134_dzOO,             416
.equiv          nb134_dxH1H1,           432
.equiv          nb134_dyH1H1,           448
.equiv          nb134_dzH1H1,           464
.equiv          nb134_dxH1H2,           480
.equiv          nb134_dyH1H2,           496
.equiv          nb134_dzH1H2,           512
.equiv          nb134_dxH1M,            528
.equiv          nb134_dyH1M,            544
.equiv          nb134_dzH1M,            560
.equiv          nb134_dxH2H1,           576
.equiv          nb134_dyH2H1,           592
.equiv          nb134_dzH2H1,           608
.equiv          nb134_dxH2H2,           624
.equiv          nb134_dyH2H2,           640
.equiv          nb134_dzH2H2,           656
.equiv          nb134_dxH2M,            672
.equiv          nb134_dyH2M,            688
.equiv          nb134_dzH2M,            704
.equiv          nb134_dxMH1,            720
.equiv          nb134_dyMH1,            736
.equiv          nb134_dzMH1,            752
.equiv          nb134_dxMH2,            768
.equiv          nb134_dyMH2,            784
.equiv          nb134_dzMH2,            800
.equiv          nb134_dxMM,             816
.equiv          nb134_dyMM,             832
.equiv          nb134_dzMM,             848
.equiv          nb134_qqMM,             864
.equiv          nb134_qqMH,             880
.equiv          nb134_qqHH,             896
.equiv          nb134_two,              912
.equiv          nb134_c6,               944
.equiv          nb134_c12,              960
.equiv          nb134_vctot,            976
.equiv          nb134_Vvdwtot,          992
.equiv          nb134_fixO,             1008
.equiv          nb134_fiyO,             1024
.equiv          nb134_fizO,             1040
.equiv          nb134_fixH1,            1056
.equiv          nb134_fiyH1,            1072
.equiv          nb134_fizH1,            1088
.equiv          nb134_fixH2,            1104
.equiv          nb134_fiyH2,            1120
.equiv          nb134_fizH2,            1136
.equiv          nb134_fixM,             1152
.equiv          nb134_fiyM,             1168
.equiv          nb134_fizM,             1184
.equiv          nb134_fjxO,             1200
.equiv          nb134_fjyO,             1216
.equiv          nb134_fjzO,             1232
.equiv          nb134_fjxH1,            1248
.equiv          nb134_fjyH1,            1264
.equiv          nb134_fjzH1,            1280
.equiv          nb134_fjxH2,            1296
.equiv          nb134_fjyH2,            1312
.equiv          nb134_fjzH2,            1328
.equiv          nb134_fjxM,             1344
.equiv          nb134_fjyM,             1360
.equiv          nb134_fjzM,             1376
.equiv          nb134_half,             1392
.equiv          nb134_three,            1408
.equiv          nb134_tsc,              1424
.equiv          nb134_fstmp,            1440
.equiv          nb134_rsqOO,            1456
.equiv          nb134_rsqH1H1,          1472
.equiv          nb134_rsqH1H2,          1488
.equiv          nb134_rsqH1M,           1504
.equiv          nb134_rsqH2H1,          1520
.equiv          nb134_rsqH2H2,          1536
.equiv          nb134_rsqH2M,           1552
.equiv          nb134_rsqMH1,           1568
.equiv          nb134_rsqMH2,           1584
.equiv          nb134_rsqMM,            1600
.equiv          nb134_rinvOO,           1616
.equiv          nb134_rinvH1H1,         1632
.equiv          nb134_rinvH1H2,         1648
.equiv          nb134_rinvH1M,          1664
.equiv          nb134_rinvH2H1,         1680
.equiv          nb134_rinvH2H2,         1696
.equiv          nb134_rinvH2M,          1712
.equiv          nb134_rinvMH1,          1728
.equiv          nb134_rinvMH2,          1744
.equiv          nb134_rinvMM,           1760
.equiv          nb134_is3,              1808
.equiv          nb134_ii3,              1812
.equiv          nb134_innerjjnr,        1816
.equiv          nb134_innerk,           1820
.equiv          nb134_n,                1824
.equiv          nb134_nn1,              1828
.equiv          nb134_nri,              1832
.equiv          nb134_nouter,           1836
.equiv          nb134_ninner,           1840
.equiv          nb134_salign,           1844
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1848		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb134_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb134_p_nri]
	mov ecx, [ecx]
	mov [esp + nb134_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb134_nouter], eax
	mov [esp + nb134_ninner], eax

	mov eax, [ebp + nb134_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb134_tsc], xmm3


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb134_half], eax
	mov [esp + nb134_half+4], ebx
	movsd xmm1, [esp + nb134_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb134_half], xmm1
	movapd [esp + nb134_two], xmm2
	movapd [esp + nb134_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb134_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb134_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb134_p_facel]
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
	movapd [esp + nb134_qqMM], xmm3
	movapd [esp + nb134_qqMH], xmm4
	movapd [esp + nb134_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb134_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb134_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb134_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movhpd xmm0, [eax + edx*8 + 8]
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb134_c6], xmm0
	movapd [esp + nb134_c12], xmm1

.nb134_threadloop:
        mov   esi, [ebp + nb134_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb134_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb134_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb134_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb134_n], eax
        mov [esp + nb134_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb134_outerstart
        jmp .nb134_end

.nb134_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb134_nouter]
	mov [esp + nb134_nouter], ebx

.nb134_outer:
	mov   eax, [ebp + nb134_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb134_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb134_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb134_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb134_pos]    ;# eax = base of pos[]  
	mov   [esp + nb134_ii3], ebx		

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
	movapd [esp + nb134_ixO], xmm3
	movapd [esp + nb134_iyO], xmm4
	movapd [esp + nb134_izO], xmm5
	movapd [esp + nb134_ixH1], xmm6
	movapd [esp + nb134_iyH1], xmm7

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
	movapd [esp + nb134_izH1], xmm6
	movapd [esp + nb134_ixH2], xmm0
	movapd [esp + nb134_iyH2], xmm1
	movapd [esp + nb134_izH2], xmm2
	movapd [esp + nb134_ixM], xmm3
	movapd [esp + nb134_iyM], xmm4
	movapd [esp + nb134_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb134_vctot], xmm4
	movapd [esp + nb134_Vvdwtot], xmm4
	movapd [esp + nb134_fixO], xmm4
	movapd [esp + nb134_fiyO], xmm4
	movapd [esp + nb134_fizO], xmm4
	movapd [esp + nb134_fixH1], xmm4
	movapd [esp + nb134_fiyH1], xmm4
	movapd [esp + nb134_fizH1], xmm4
	movapd [esp + nb134_fixH2], xmm4
	movapd [esp + nb134_fiyH2], xmm4
	movapd [esp + nb134_fizH2], xmm4
	movapd [esp + nb134_fixM], xmm4
	movapd [esp + nb134_fiyM], xmm4
	movapd [esp + nb134_fizM], xmm4
	
	mov   eax, [ebp + nb134_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb134_pos] 
	mov   edi, [ebp + nb134_faction]	
	mov   eax, [ebp + nb134_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb134_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb134_ninner]
	mov   [esp + nb134_ninner], ecx
	add   edx, 0
	mov   [esp + nb134_innerk], edx    ;# number of innerloop atoms 
	jge   .nb134_unroll_loop
	jmp   .nb134_checksingle
.nb134_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb134_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb134_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb134_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb134_jxO], xmm0
	movapd 	[esp + nb134_jyO], xmm1
	movapd 	[esp + nb134_jzO], xmm3
	movapd 	[esp + nb134_jxH1], xmm4
	
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
	movapd 	[esp + nb134_jyH1], xmm0
	movapd 	[esp + nb134_jzH1], xmm1
	movapd 	[esp + nb134_jxH2], xmm3
	movapd 	[esp + nb134_jyH2], xmm4
	
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
	movapd 	[esp + nb134_jzH2], xmm0
	movapd 	[esp + nb134_jxM], xmm1
	movapd 	[esp + nb134_jyM], xmm3
	movapd 	[esp + nb134_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [esp + nb134_ixO]
	movapd xmm1, [esp + nb134_iyO]
	movapd xmm2, [esp + nb134_izO]
	movapd xmm3, [esp + nb134_ixH1]
	movapd xmm4, [esp + nb134_iyH1]
	movapd xmm5, [esp + nb134_izH1]
	subpd  xmm0, [esp + nb134_jxO]
	subpd  xmm1, [esp + nb134_jyO]
	subpd  xmm2, [esp + nb134_jzO]
	subpd  xmm3, [esp + nb134_jxH1]
	subpd  xmm4, [esp + nb134_jyH1]
	subpd  xmm5, [esp + nb134_jzH1]
	movapd [esp + nb134_dxOO], xmm0
	movapd [esp + nb134_dyOO], xmm1
	movapd [esp + nb134_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb134_dxH1H1], xmm3
	movapd [esp + nb134_dyH1H1], xmm4
	movapd [esp + nb134_dzH1H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb134_rsqOO], xmm0
	movapd [esp + nb134_rsqH1H1], xmm3

	movapd xmm0, [esp + nb134_ixH1]
	movapd xmm1, [esp + nb134_iyH1]
	movapd xmm2, [esp + nb134_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb134_jxH2]
	subpd  xmm1, [esp + nb134_jyH2]
	subpd  xmm2, [esp + nb134_jzH2]
	subpd  xmm3, [esp + nb134_jxM]
	subpd  xmm4, [esp + nb134_jyM]
	subpd  xmm5, [esp + nb134_jzM]
	movapd [esp + nb134_dxH1H2], xmm0
	movapd [esp + nb134_dyH1H2], xmm1
	movapd [esp + nb134_dzH1H2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb134_dxH1M], xmm3
	movapd [esp + nb134_dyH1M], xmm4
	movapd [esp + nb134_dzH1M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb134_rsqH1H2], xmm0
	movapd [esp + nb134_rsqH1M], xmm3

	movapd xmm0, [esp + nb134_ixH2]
	movapd xmm1, [esp + nb134_iyH2]
	movapd xmm2, [esp + nb134_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb134_jxH1]
	subpd  xmm1, [esp + nb134_jyH1]
	subpd  xmm2, [esp + nb134_jzH1]
	subpd  xmm3, [esp + nb134_jxH2]
	subpd  xmm4, [esp + nb134_jyH2]
	subpd  xmm5, [esp + nb134_jzH2]
	movapd [esp + nb134_dxH2H1], xmm0
	movapd [esp + nb134_dyH2H1], xmm1
	movapd [esp + nb134_dzH2H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb134_dxH2H2], xmm3
	movapd [esp + nb134_dyH2H2], xmm4
	movapd [esp + nb134_dzH2H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb134_rsqH2H1], xmm0
	movapd [esp + nb134_rsqH2H2], xmm3

	movapd xmm0, [esp + nb134_ixH2]
	movapd xmm1, [esp + nb134_iyH2]
	movapd xmm2, [esp + nb134_izH2]
	movapd xmm3, [esp + nb134_ixM]
	movapd xmm4, [esp + nb134_iyM]
	movapd xmm5, [esp + nb134_izM]
	subpd  xmm0, [esp + nb134_jxM]
	subpd  xmm1, [esp + nb134_jyM]
	subpd  xmm2, [esp + nb134_jzM]
	subpd  xmm3, [esp + nb134_jxH1]
	subpd  xmm4, [esp + nb134_jyH1]
	subpd  xmm5, [esp + nb134_jzH1]
	movapd [esp + nb134_dxH2M], xmm0
	movapd [esp + nb134_dyH2M], xmm1
	movapd [esp + nb134_dzH2M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb134_dxMH1], xmm3
	movapd [esp + nb134_dyMH1], xmm4
	movapd [esp + nb134_dzMH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb134_rsqH2M], xmm0
	movapd [esp + nb134_rsqMH1], xmm4

	movapd xmm0, [esp + nb134_ixM]
	movapd xmm1, [esp + nb134_iyM]
	movapd xmm2, [esp + nb134_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb134_jxH2]
	subpd  xmm1, [esp + nb134_jyH2]
	subpd  xmm2, [esp + nb134_jzH2]
	subpd  xmm3, [esp + nb134_jxM]
	subpd  xmm4, [esp + nb134_jyM]
	subpd  xmm5, [esp + nb134_jzM]
	movapd [esp + nb134_dxMH2], xmm0
	movapd [esp + nb134_dyMH2], xmm1
	movapd [esp + nb134_dzMH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb134_dxMM], xmm3
	movapd [esp + nb134_dyMM], xmm4
	movapd [esp + nb134_dzMM], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb134_rsqMH2], xmm0
	movapd [esp + nb134_rsqMM], xmm4
	
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
	movapd  xmm3, [esp + nb134_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134_half] ;# iter1 
	mulpd   xmm7, [esp + nb134_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134_half] ;# rinv 
	mulpd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvMH2], xmm1
	movapd [esp + nb134_rinvMM], xmm5

	movapd xmm0, [esp + nb134_rsqOO]
	movapd xmm4, [esp + nb134_rsqH1H1]	
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
	movapd  xmm3, [esp + nb134_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb134_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134_half] ;# rinv 
	mulpd   xmm5, [esp + nb134_half] ;# rinv
	movapd [esp + nb134_rinvOO], xmm1
	movapd [esp + nb134_rinvH1H1], xmm5

	movapd xmm0, [esp + nb134_rsqH1H2]
	movapd xmm4, [esp + nb134_rsqH1M]	
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
	movapd  xmm3, [esp + nb134_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134_half] ;# iter1 
	mulpd   xmm7, [esp + nb134_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134_half] ;# rinv 
	mulpd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvH1H2], xmm1
	movapd [esp + nb134_rinvH1M], xmm5

	movapd xmm0, [esp + nb134_rsqH2H1]
	movapd xmm4, [esp + nb134_rsqH2H2]	
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
	movapd  xmm3, [esp + nb134_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134_half] ;# iter1a 
	mulpd   xmm7, [esp + nb134_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134_half] ;# rinv 
	mulpd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvH2H1], xmm1
	movapd [esp + nb134_rinvH2H2], xmm5

	movapd xmm0, [esp + nb134_rsqMH1]
	movapd xmm4, [esp + nb134_rsqH2M]	
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
	movapd  xmm3, [esp + nb134_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134_half] ;# iter1a 
	mulpd   xmm7, [esp + nb134_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134_half] ;# rinv 
	mulpd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvMH1], xmm1
	movapd [esp + nb134_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb134_rinvOO] 
	movapd xmm4, [esp + nb134_rsqOO]
	
		mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb134_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb134_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb134_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb134_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb134_Vvdwtot]
	xorpd  xmm3, xmm3
	mulpd  xmm7, [esp + nb134_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb134_fstmp], xmm3
	movapd [esp + nb134_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb134_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb134_c12]
	mulpd  xmm7, xmm4 
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb134_Vvdwtot]
	movapd xmm3, [esp + nb134_fstmp]
	mulpd  xmm7, [esp + nb134_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb134_Vvdwtot], xmm5

	mulpd  xmm0, xmm3
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movd eax, mm0	
	movd ebx, mm1

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb134_dxOO]
	mulpd xmm1, [esp + nb134_dyOO]
	mulpd xmm2, [esp + nb134_dzOO]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixO]
	addpd xmm1, [esp + nb134_fiyO]
	addpd xmm2, [esp + nb134_fizO]
	movapd [esp + nb134_fjxO], xmm3
	movapd [esp + nb134_fjyO], xmm4
	movapd [esp + nb134_fjzO], xmm5
	movapd [esp + nb134_fixO], xmm0
	movapd [esp + nb134_fiyO], xmm1
	movapd [esp + nb134_fizO], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb134_rinvH1H1]
	movapd xmm6, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm6, [esp + nb134_qqHH] ;# vcoul
	mulpd xmm0, xmm6
	
	addpd  xmm6, [esp + nb134_vctot] ;# local vctot summation variable 

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb134_dxH1H1]
	mulpd xmm1, [esp + nb134_dyH1H1]
	mulpd xmm2, [esp + nb134_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixH1]
	addpd xmm1, [esp + nb134_fiyH1]
	addpd xmm2, [esp + nb134_fizH1]
	movapd [esp + nb134_fjxH1], xmm3
	movapd [esp + nb134_fjyH1], xmm4
	movapd [esp + nb134_fjzH1], xmm5
	movapd [esp + nb134_fixH1], xmm0
	movapd [esp + nb134_fiyH1], xmm1
	movapd [esp + nb134_fizH1], xmm2

	;# H1-H2 interaction  
	movapd xmm0, [esp + nb134_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqHH] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb134_dxH1H2]
	mulpd xmm1, [esp + nb134_dyH1H2]
	mulpd xmm2, [esp + nb134_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixH1]
	addpd xmm1, [esp + nb134_fiyH1]
	addpd xmm2, [esp + nb134_fizH1]
	movapd [esp + nb134_fjxH2], xmm3
	movapd [esp + nb134_fjyH2], xmm4
	movapd [esp + nb134_fjzH2], xmm5
	movapd [esp + nb134_fixH1], xmm0
	movapd [esp + nb134_fiyH1], xmm1
	movapd [esp + nb134_fizH1], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb134_rinvH1M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqMH] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb134_dxH1M]
	mulpd xmm1, [esp + nb134_dyH1M]
	mulpd xmm2, [esp + nb134_dzH1M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixH1]
	addpd xmm1, [esp + nb134_fiyH1]
	addpd xmm2, [esp + nb134_fizH1]
	movapd [esp + nb134_fjxM], xmm3
	movapd [esp + nb134_fjyM], xmm4
	movapd [esp + nb134_fjzM], xmm5
	movapd [esp + nb134_fixH1], xmm0
	movapd [esp + nb134_fiyH1], xmm1
	movapd [esp + nb134_fizH1], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb134_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqHH] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxH1]
	movapd xmm4, [esp + nb134_fjyH1]
	movapd xmm5, [esp + nb134_fjzH1]
	mulpd xmm0, [esp + nb134_dxH2H1]
	mulpd xmm1, [esp + nb134_dyH2H1]
	mulpd xmm2, [esp + nb134_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixH2]
	addpd xmm1, [esp + nb134_fiyH2]
	addpd xmm2, [esp + nb134_fizH2]
	movapd [esp + nb134_fjxH1], xmm3
	movapd [esp + nb134_fjyH1], xmm4
	movapd [esp + nb134_fjzH1], xmm5
	movapd [esp + nb134_fixH2], xmm0
	movapd [esp + nb134_fiyH2], xmm1
	movapd [esp + nb134_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb134_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqHH] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxH2]
	movapd xmm4, [esp + nb134_fjyH2]
	movapd xmm5, [esp + nb134_fjzH2]
	mulpd xmm0, [esp + nb134_dxH2H2]
	mulpd xmm1, [esp + nb134_dyH2H2]
	mulpd xmm2, [esp + nb134_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixH2]
	addpd xmm1, [esp + nb134_fiyH2]
	addpd xmm2, [esp + nb134_fizH2]
	movapd [esp + nb134_fjxH2], xmm3
	movapd [esp + nb134_fjyH2], xmm4
	movapd [esp + nb134_fjzH2], xmm5
	movapd [esp + nb134_fixH2], xmm0
	movapd [esp + nb134_fiyH2], xmm1
	movapd [esp + nb134_fizH2], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb134_rinvH2M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqMH] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxM]
	movapd xmm4, [esp + nb134_fjyM]
	movapd xmm5, [esp + nb134_fjzM]
	mulpd xmm0, [esp + nb134_dxH2M]
	mulpd xmm1, [esp + nb134_dyH2M]
	mulpd xmm2, [esp + nb134_dzH2M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixH2]
	addpd xmm1, [esp + nb134_fiyH2]
	addpd xmm2, [esp + nb134_fizH2]
	movapd [esp + nb134_fjxM], xmm3
	movapd [esp + nb134_fjyM], xmm4
	movapd [esp + nb134_fjzM], xmm5
	movapd [esp + nb134_fixH2], xmm0
	movapd [esp + nb134_fiyH2], xmm1
	movapd [esp + nb134_fizH2], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb134_rinvMH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqMH] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxH1]
	movapd xmm4, [esp + nb134_fjyH1]
	movapd xmm5, [esp + nb134_fjzH1]
	mulpd xmm0, [esp + nb134_dxMH1]
	mulpd xmm1, [esp + nb134_dyMH1]
	mulpd xmm2, [esp + nb134_dzMH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixM]
	addpd xmm1, [esp + nb134_fiyM]
	addpd xmm2, [esp + nb134_fizM]
	movapd [esp + nb134_fjxH1], xmm3
	movapd [esp + nb134_fjyH1], xmm4
	movapd [esp + nb134_fjzH1], xmm5
	movapd [esp + nb134_fixM], xmm0
	movapd [esp + nb134_fiyM], xmm1
	movapd [esp + nb134_fizM], xmm2

	;# M-H2 interaction 
	movapd xmm0, [esp + nb134_rinvMH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqMH] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7

	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxH2]
	movapd xmm4, [esp + nb134_fjyH2]
	movapd xmm5, [esp + nb134_fjzH2]
	mulpd xmm0, [esp + nb134_dxMH2]
	mulpd xmm1, [esp + nb134_dyMH2]
	mulpd xmm2, [esp + nb134_dzMH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixM]
	addpd xmm1, [esp + nb134_fiyM]
	addpd xmm2, [esp + nb134_fizM]
	movapd [esp + nb134_fjxH2], xmm3
	movapd [esp + nb134_fjyH2], xmm4
	movapd [esp + nb134_fjzH2], xmm5
	movapd [esp + nb134_fixM], xmm0
	movapd [esp + nb134_fiyM], xmm1
	movapd [esp + nb134_fizM], xmm2

	;# M-M interaction 
	movapd xmm0, [esp + nb134_rinvMM]
	movapd xmm7, xmm0	;# xmm7=rinv 
	mulpd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulpd  xmm7, [esp + nb134_qqMM] ;# vcoul
	mulpd xmm0, xmm7
	addpd  xmm6, xmm7
	movapd [esp + nb134_vctot], xmm6
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxM]
	movapd xmm4, [esp + nb134_fjyM] 
	movapd xmm5, [esp + nb134_fjzM]
	mulpd xmm0, [esp + nb134_dxMM]
	mulpd xmm1, [esp + nb134_dyMM]
	mulpd xmm2, [esp + nb134_dzMM]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb134_fixM]
	addpd xmm1, [esp + nb134_fiyM]
	addpd xmm2, [esp + nb134_fizM]
	movapd [esp + nb134_fjxM], xmm3
	movapd [esp + nb134_fjyM], xmm4
	movapd [esp + nb134_fjzM], xmm5
	movapd [esp + nb134_fixM], xmm0
	movapd [esp + nb134_fiyM], xmm1
	movapd [esp + nb134_fizM], xmm2

	mov edi, [ebp + nb134_faction]
	
	;# Did all interactions - now update j forces 
	;# Step1 - transpose fjxO, fjyO and fjzO, fjxH1
	movapd xmm0, [esp + nb134_fjxO]
	movapd xmm1, [esp + nb134_fjyO]
	movapd xmm2, [esp + nb134_fjzO]
	movapd xmm3, [esp + nb134_fjxH1]
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
	movapd xmm0, [esp + nb134_fjyH1]
	movapd xmm1, [esp + nb134_fjzH1]
	movapd xmm2, [esp + nb134_fjxH2]
	movapd xmm3, [esp + nb134_fjyH2]
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
	movapd xmm0, [esp + nb134_fjzH2]
	movapd xmm1, [esp + nb134_fjxM]
	movapd xmm2, [esp + nb134_fjyM]
	movapd xmm3, [esp + nb134_fjzM]
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
	sub dword ptr [esp + nb134_innerk],  2
	jl    .nb134_checksingle
	jmp   .nb134_unroll_loop
.nb134_checksingle:
	mov   edx, [esp + nb134_innerk]
	and   edx, 1
	jnz   .nb134_dosingle
	jmp   .nb134_updateouterdata
.nb134_dosingle:
	mov   edx, [esp + nb134_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb134_pos]       ;# base of pos[] 

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
	movsd  [esp + nb134_jxO], xmm0
	movsd  [esp + nb134_jzO], xmm1
	movsd  [esp + nb134_jyH1], xmm2
	movsd  [esp + nb134_jxH2], xmm3
	movsd  [esp + nb134_jzH2], xmm4
	movsd  [esp + nb134_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [esp + nb134_jyO], xmm0
	movsd  [esp + nb134_jxH1], xmm1
	movsd  [esp + nb134_jzH1], xmm2
	movsd  [esp + nb134_jyH2], xmm3
	movsd  [esp + nb134_jxM], xmm4
	movsd  [esp + nb134_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [esp + nb134_ixO]
	movapd xmm1, [esp + nb134_iyO]
	movapd xmm2, [esp + nb134_izO]
	movapd xmm3, [esp + nb134_ixH1]
	movapd xmm4, [esp + nb134_iyH1]
	movapd xmm5, [esp + nb134_izH1]
	subsd  xmm0, [esp + nb134_jxO]
	subsd  xmm1, [esp + nb134_jyO]
	subsd  xmm2, [esp + nb134_jzO]
	subsd  xmm3, [esp + nb134_jxH1]
	subsd  xmm4, [esp + nb134_jyH1]
	subsd  xmm5, [esp + nb134_jzH1]
	movapd [esp + nb134_dxOO], xmm0
	movapd [esp + nb134_dyOO], xmm1
	movapd [esp + nb134_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb134_dxH1H1], xmm3
	movapd [esp + nb134_dyH1H1], xmm4
	movapd [esp + nb134_dzH1H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb134_rsqOO], xmm0
	movapd [esp + nb134_rsqH1H1], xmm3

	movapd xmm0, [esp + nb134_ixH1]
	movapd xmm1, [esp + nb134_iyH1]
	movapd xmm2, [esp + nb134_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb134_jxH2]
	subsd  xmm1, [esp + nb134_jyH2]
	subsd  xmm2, [esp + nb134_jzH2]
	subsd  xmm3, [esp + nb134_jxM]
	subsd  xmm4, [esp + nb134_jyM]
	subsd  xmm5, [esp + nb134_jzM]
	movapd [esp + nb134_dxH1H2], xmm0
	movapd [esp + nb134_dyH1H2], xmm1
	movapd [esp + nb134_dzH1H2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb134_dxH1M], xmm3
	movapd [esp + nb134_dyH1M], xmm4
	movapd [esp + nb134_dzH1M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb134_rsqH1H2], xmm0
	movapd [esp + nb134_rsqH1M], xmm3

	movapd xmm0, [esp + nb134_ixH2]
	movapd xmm1, [esp + nb134_iyH2]
	movapd xmm2, [esp + nb134_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb134_jxH1]
	subsd  xmm1, [esp + nb134_jyH1]
	subsd  xmm2, [esp + nb134_jzH1]
	subsd  xmm3, [esp + nb134_jxH2]
	subsd  xmm4, [esp + nb134_jyH2]
	subsd  xmm5, [esp + nb134_jzH2]
	movapd [esp + nb134_dxH2H1], xmm0
	movapd [esp + nb134_dyH2H1], xmm1
	movapd [esp + nb134_dzH2H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb134_dxH2H2], xmm3
	movapd [esp + nb134_dyH2H2], xmm4
	movapd [esp + nb134_dzH2H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb134_rsqH2H1], xmm0
	movapd [esp + nb134_rsqH2H2], xmm3

	movapd xmm0, [esp + nb134_ixH2]
	movapd xmm1, [esp + nb134_iyH2]
	movapd xmm2, [esp + nb134_izH2]
	movapd xmm3, [esp + nb134_ixM]
	movapd xmm4, [esp + nb134_iyM]
	movapd xmm5, [esp + nb134_izM]
	subsd  xmm0, [esp + nb134_jxM]
	subsd  xmm1, [esp + nb134_jyM]
	subsd  xmm2, [esp + nb134_jzM]
	subsd  xmm3, [esp + nb134_jxH1]
	subsd  xmm4, [esp + nb134_jyH1]
	subsd  xmm5, [esp + nb134_jzH1]
	movapd [esp + nb134_dxH2M], xmm0
	movapd [esp + nb134_dyH2M], xmm1
	movapd [esp + nb134_dzH2M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb134_dxMH1], xmm3
	movapd [esp + nb134_dyMH1], xmm4
	movapd [esp + nb134_dzMH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb134_rsqH2M], xmm0
	movapd [esp + nb134_rsqMH1], xmm4

	movapd xmm0, [esp + nb134_ixM]
	movapd xmm1, [esp + nb134_iyM]
	movapd xmm2, [esp + nb134_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb134_jxH2]
	subsd  xmm1, [esp + nb134_jyH2]
	subsd  xmm2, [esp + nb134_jzH2]
	subsd  xmm3, [esp + nb134_jxM]
	subsd  xmm4, [esp + nb134_jyM]
	subsd  xmm5, [esp + nb134_jzM]
	movapd [esp + nb134_dxMH2], xmm0
	movapd [esp + nb134_dyMH2], xmm1
	movapd [esp + nb134_dzMH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb134_dxMM], xmm3
	movapd [esp + nb134_dyMM], xmm4
	movapd [esp + nb134_dzMM], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb134_rsqMH2], xmm0
	movapd [esp + nb134_rsqMM], xmm4

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
	movapd  xmm3, [esp + nb134_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134_half] ;# iter1 
	mulsd   xmm7, [esp + nb134_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134_half] ;# rinv 
	mulsd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvMH2], xmm1
	movapd [esp + nb134_rinvMM], xmm5

	movapd xmm0, [esp + nb134_rsqOO]
	movapd xmm4, [esp + nb134_rsqH1H1]	
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
	movapd  xmm3, [esp + nb134_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb134_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134_half] ;# rinv 
	mulsd   xmm5, [esp + nb134_half] ;# rinv
	movapd [esp + nb134_rinvOO], xmm1
	movapd [esp + nb134_rinvH1H1], xmm5

	movapd xmm0, [esp + nb134_rsqH1H2]
	movapd xmm4, [esp + nb134_rsqH1M]	
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
	movapd  xmm3, [esp + nb134_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134_half] ;# iter1 
	mulsd   xmm7, [esp + nb134_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134_half] ;# rinv 
	mulsd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvH1H2], xmm1
	movapd [esp + nb134_rinvH1M], xmm5

	movapd xmm0, [esp + nb134_rsqH2H1]
	movapd xmm4, [esp + nb134_rsqH2H2]	
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
	movapd  xmm3, [esp + nb134_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134_half] ;# iter1a 
	mulsd   xmm7, [esp + nb134_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134_half] ;# rinv 
	mulsd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvH2H1], xmm1
	movapd [esp + nb134_rinvH2H2], xmm5

	movapd xmm0, [esp + nb134_rsqMH1]
	movapd xmm4, [esp + nb134_rsqH2M]	
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
	movapd  xmm3, [esp + nb134_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134_half] ;# iter1a 
	mulsd   xmm7, [esp + nb134_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134_half] ;# rinv 
	mulsd   xmm5, [esp + nb134_half] ;# rinv 
	movapd [esp + nb134_rinvMH1], xmm1
	movapd [esp + nb134_rinvH2M], xmm5

	;# start with OO interaction 
	movsd xmm0, [esp + nb134_rinvOO] 
	movsd xmm4, [esp + nb134_rsqOO]
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb134_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb134_VFtab]

	;# dispersion 
	movlpd xmm4, [esi + ebx*8]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 16]	;# G1
	movhpd xmm6, [esi + ebx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb134_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb134_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb134_Vvdwtot]
	xorpd  xmm3, xmm3
	mulsd  xmm7, [esp + nb134_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb134_fstmp], xmm3
	movsd [esp + nb134_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + ebx*8 + 32]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 48]	;# G1
	movhpd xmm6, [esi + ebx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb134_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb134_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb134_Vvdwtot]
	movsd xmm3, [esp + nb134_fstmp]
	mulsd  xmm7, [esp + nb134_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb134_Vvdwtot], xmm5

	mulsd  xmm0, xmm3
	movsd xmm1, xmm0
	movsd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb134_dxOO]
	mulsd xmm1, [esp + nb134_dyOO]
	mulsd xmm2, [esp + nb134_dzOO]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixO]
	addsd xmm1, [esp + nb134_fiyO]
	addsd xmm2, [esp + nb134_fizO]
	movsd [esp + nb134_fjxO], xmm3
	movsd [esp + nb134_fjyO], xmm4
	movsd [esp + nb134_fjzO], xmm5
	movsd [esp + nb134_fixO], xmm0
	movsd [esp + nb134_fiyO], xmm1
	movsd [esp + nb134_fizO], xmm2

	;# H1-H1 interaction 
	movsd  xmm0, [esp + nb134_rinvH1H1]
	movsd  xmm6, xmm0	;# xmm6=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm6, [esp + nb134_qqHH] ;# vcoul
	mulsd xmm0, xmm6
	addsd  xmm6, [esp + nb134_vctot]

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb134_dxH1H1]
	mulsd xmm1, [esp + nb134_dyH1H1]
	mulsd xmm2, [esp + nb134_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixH1]
	addsd xmm1, [esp + nb134_fiyH1]
	addsd xmm2, [esp + nb134_fizH1]
	movsd [esp + nb134_fjxH1], xmm3
	movsd [esp + nb134_fjyH1], xmm4
	movsd [esp + nb134_fjzH1], xmm5
	movsd [esp + nb134_fixH1], xmm0
	movsd [esp + nb134_fiyH1], xmm1
	movsd [esp + nb134_fizH1], xmm2

	;# H1-H2 interaction  
	movsd  xmm0, [esp + nb134_rinvH1H2]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqHH] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb134_dxH1H2]
	mulsd xmm1, [esp + nb134_dyH1H2]
	mulsd xmm2, [esp + nb134_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixH1]
	addsd xmm1, [esp + nb134_fiyH1]
	addsd xmm2, [esp + nb134_fizH1]
	movsd [esp + nb134_fjxH2], xmm3
	movsd [esp + nb134_fjyH2], xmm4
	movsd [esp + nb134_fjzH2], xmm5
	movsd [esp + nb134_fixH1], xmm0
	movsd [esp + nb134_fiyH1], xmm1
	movsd [esp + nb134_fizH1], xmm2

	;# H1-M interaction 
	movsd  xmm0, [esp + nb134_rinvH1M]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqMH] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb134_dxH1M]
	mulsd xmm1, [esp + nb134_dyH1M]
	mulsd xmm2, [esp + nb134_dzH1M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixH1]
	addsd xmm1, [esp + nb134_fiyH1]
	addsd xmm2, [esp + nb134_fizH1]
	movsd [esp + nb134_fjxM], xmm3
	movsd [esp + nb134_fjyM], xmm4
	movsd [esp + nb134_fjzM], xmm5
	movsd [esp + nb134_fixH1], xmm0
	movsd [esp + nb134_fiyH1], xmm1
	movsd [esp + nb134_fizH1], xmm2

	;# H2-H1 interaction 
	movsd  xmm0, [esp + nb134_rinvH2H1]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqHH] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxH1]
	movapd xmm4, [esp + nb134_fjyH1]
	movapd xmm5, [esp + nb134_fjzH1]
	mulsd xmm0, [esp + nb134_dxH2H1]
	mulsd xmm1, [esp + nb134_dyH2H1]
	mulsd xmm2, [esp + nb134_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixH2]
	addsd xmm1, [esp + nb134_fiyH2]
	addsd xmm2, [esp + nb134_fizH2]
	movsd [esp + nb134_fjxH1], xmm3
	movsd [esp + nb134_fjyH1], xmm4
	movsd [esp + nb134_fjzH1], xmm5
	movsd [esp + nb134_fixH2], xmm0
	movsd [esp + nb134_fiyH2], xmm1
	movsd [esp + nb134_fizH2], xmm2

	;# H2-H2 interaction 
	movsd  xmm0, [esp + nb134_rinvH2H2]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqHH] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	movsd xmm3, [esp + nb134_fjxH2]
	movsd xmm4, [esp + nb134_fjyH2]
	movsd xmm5, [esp + nb134_fjzH2]
	mulsd xmm0, [esp + nb134_dxH2H2]
	mulsd xmm1, [esp + nb134_dyH2H2]
	mulsd xmm2, [esp + nb134_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixH2]
	addsd xmm1, [esp + nb134_fiyH2]
	addsd xmm2, [esp + nb134_fizH2]
	movsd [esp + nb134_fjxH2], xmm3
	movsd [esp + nb134_fjyH2], xmm4
	movsd [esp + nb134_fjzH2], xmm5
	movsd [esp + nb134_fixH2], xmm0
	movsd [esp + nb134_fiyH2], xmm1
	movsd [esp + nb134_fizH2], xmm2

	;# H2-M interaction 
	movsd  xmm0, [esp + nb134_rinvH2M]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqMH] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxM]
	movapd xmm4, [esp + nb134_fjyM]
	movapd xmm5, [esp + nb134_fjzM]
	mulsd xmm0, [esp + nb134_dxH2M]
	mulsd xmm1, [esp + nb134_dyH2M]
	mulsd xmm2, [esp + nb134_dzH2M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixH2]
	addsd xmm1, [esp + nb134_fiyH2]
	addsd xmm2, [esp + nb134_fizH2]
	movsd [esp + nb134_fjxM], xmm3
	movsd [esp + nb134_fjyM], xmm4
	movsd [esp + nb134_fjzM], xmm5
	movsd [esp + nb134_fixH2], xmm0
	movsd [esp + nb134_fiyH2], xmm1
	movsd [esp + nb134_fizH2], xmm2

	;# M-H1 interaction 
	movsd  xmm0, [esp + nb134_rinvMH1]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqMH] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxH1]
	movapd xmm4, [esp + nb134_fjyH1]
	movapd xmm5, [esp + nb134_fjzH1]
	mulsd xmm0, [esp + nb134_dxMH1]
	mulsd xmm1, [esp + nb134_dyMH1]
	mulsd xmm2, [esp + nb134_dzMH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixM]
	addsd xmm1, [esp + nb134_fiyM]
	addsd xmm2, [esp + nb134_fizM]
	movsd [esp + nb134_fjxH1], xmm3
	movsd [esp + nb134_fjyH1], xmm4
	movsd [esp + nb134_fjzH1], xmm5
	movsd [esp + nb134_fixM], xmm0
	movsd [esp + nb134_fiyM], xmm1
	movsd [esp + nb134_fizM], xmm2

	;# M-H2 interaction 
	movsd  xmm0, [esp + nb134_rinvMH2]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqMH] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4

	movsd xmm1, xmm0
	movsd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxH2]
	movapd xmm4, [esp + nb134_fjyH2]
	movapd xmm5, [esp + nb134_fjzH2]
	mulsd xmm0, [esp + nb134_dxMH2]
	mulsd xmm1, [esp + nb134_dyMH2]
	mulsd xmm2, [esp + nb134_dzMH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixM]
	addsd xmm1, [esp + nb134_fiyM]
	addsd xmm2, [esp + nb134_fizM]
	movsd [esp + nb134_fjxH2], xmm3
	movsd [esp + nb134_fjyH2], xmm4
	movsd [esp + nb134_fjzH2], xmm5
	movsd [esp + nb134_fixM], xmm0
	movsd [esp + nb134_fiyM], xmm1
	movsd [esp + nb134_fizM], xmm2

	;# M-M interaction 
	movsd  xmm0, [esp + nb134_rinvMM]
	movsd  xmm4, xmm0	;# xmm7=rinv 
	mulsd  xmm0, xmm0	;# xmm0=rinvsq 
	
	mulsd  xmm4, [esp + nb134_qqMM] ;# vcoul
	mulsd xmm0, xmm4
	addsd  xmm6, xmm4
	movsd [esp + nb134_vctot], xmm6
	movsd xmm1, xmm0
	movsd xmm2, xmm0

	movapd xmm3, [esp + nb134_fjxM]
	movapd xmm4, [esp + nb134_fjyM]
	movapd xmm5, [esp + nb134_fjzM]
	mulsd xmm0, [esp + nb134_dxMM]
	mulsd xmm1, [esp + nb134_dyMM]
	mulsd xmm2, [esp + nb134_dzMM]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb134_fixM]
	addsd xmm1, [esp + nb134_fiyM]
	addsd xmm2, [esp + nb134_fizM]
	movsd [esp + nb134_fjxM], xmm3
	movsd [esp + nb134_fjyM], xmm4
	movsd [esp + nb134_fjzM], xmm5
	movsd [esp + nb134_fixM], xmm0
	movsd [esp + nb134_fiyM], xmm1
	movsd [esp + nb134_fizM], xmm2

	mov edi, [ebp + nb134_faction]

	;# Did all interactions - now update j forces 
	;# Step1 - merge forces
	movlpd xmm0, [esp + nb134_fjxO]
	movlpd xmm1, [esp + nb134_fjzO]
	movlpd xmm2, [esp + nb134_fjyH1]
	movlpd xmm3, [esp + nb134_fjxH2]
	movlpd xmm4, [esp + nb134_fjzH2]
	movlpd xmm5, [esp + nb134_fjyM]

	movhpd xmm0, [esp + nb134_fjyO]
	movhpd xmm1, [esp + nb134_fjxH1]
	movhpd xmm2, [esp + nb134_fjzH1]
	movhpd xmm3, [esp + nb134_fjyH2]
	movhpd xmm4, [esp + nb134_fjxM]
	movhpd xmm5, [esp + nb134_fjzM]

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
	
	movlpd [edi + eax*8],      xmm0
	movhpd [edi + eax*8 + 8],  xmm0
	movlpd [edi + eax*8 + 16], xmm1
	movhpd [edi + eax*8 + 24], xmm1
	movlpd [edi + eax*8 + 32], xmm2
	movhpd [edi + eax*8 + 40], xmm2
	movlpd [edi + eax*8 + 48], xmm3
	movhpd [edi + eax*8 + 56], xmm3
	movlpd [edi + eax*8 + 64], xmm4
	movhpd [edi + eax*8 + 72], xmm4
	movlpd [edi + eax*8 + 80], xmm5
	movhpd [edi + eax*8 + 88], xmm5
	
.nb134_updateouterdata:
	mov   ecx, [esp + nb134_ii3]
	mov   edi, [ebp + nb134_faction]
	mov   esi, [ebp + nb134_fshift]
	mov   edx, [esp + nb134_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb134_fixO]
	movapd xmm1, [esp + nb134_fiyO]
	movapd xmm2, [esp + nb134_fizO]

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
	movapd xmm0, [esp + nb134_fixH1]
	movapd xmm1, [esp + nb134_fiyH1]
	movapd xmm2, [esp + nb134_fizH1]

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
	movapd xmm0, [esp + nb134_fixH2]
	movapd xmm1, [esp + nb134_fiyH2]
	movapd xmm2, [esp + nb134_fizH2]

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
	movapd xmm0, [esp + nb134_fixM]
	movapd xmm1, [esp + nb134_fiyM]
	movapd xmm2, [esp + nb134_fizM]

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
	mov esi, [esp + nb134_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb134_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb134_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb134_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb134_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb134_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb134_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb134_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb134_n], esi
        jmp .nb134_outer
.nb134_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb134_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb134_end
        ;# non-zero, do one more workunit
        jmp   .nb134_threadloop
.nb134_end:
	emms

	mov eax, [esp + nb134_nouter]
	mov ebx, [esp + nb134_ninner]
	mov ecx, [ebp + nb134_outeriter]
	mov edx, [ebp + nb134_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb134_salign]
	add esp, eax
	add esp, 1848
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel134nf_ia32_sse2
.globl _nb_kernel134nf_ia32_sse2
nb_kernel134nf_ia32_sse2:	
_nb_kernel134nf_ia32_sse2:	
.equiv          nb134nf_p_nri,            8
.equiv          nb134nf_iinr,             12
.equiv          nb134nf_jindex,           16
.equiv          nb134nf_jjnr,             20
.equiv          nb134nf_shift,            24
.equiv          nb134nf_shiftvec,         28
.equiv          nb134nf_fshift,           32
.equiv          nb134nf_gid,              36
.equiv          nb134nf_pos,              40
.equiv          nb134nf_faction,          44
.equiv          nb134nf_charge,           48
.equiv          nb134nf_p_facel,          52
.equiv          nb134nf_argkrf,           56
.equiv          nb134nf_argcrf,           60
.equiv          nb134nf_Vc,               64
.equiv          nb134nf_type,             68
.equiv          nb134nf_p_ntype,          72
.equiv          nb134nf_vdwparam,         76
.equiv          nb134nf_Vvdw,             80
.equiv          nb134nf_p_tabscale,       84
.equiv          nb134nf_VFtab,            88
.equiv          nb134nf_invsqrta,         92
.equiv          nb134nf_dvda,             96
.equiv          nb134nf_p_gbtabscale,     100
.equiv          nb134nf_GBtab,            104
.equiv          nb134nf_p_nthreads,       108
.equiv          nb134nf_count,            112
.equiv          nb134nf_mtx,              116
.equiv          nb134nf_outeriter,        120
.equiv          nb134nf_inneriter,        124
.equiv          nb134nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb134nf_ixO,              0
.equiv          nb134nf_iyO,              16
.equiv          nb134nf_izO,              32
.equiv          nb134nf_ixH1,             48
.equiv          nb134nf_iyH1,             64
.equiv          nb134nf_izH1,             80
.equiv          nb134nf_ixH2,             96
.equiv          nb134nf_iyH2,             112
.equiv          nb134nf_izH2,             128
.equiv          nb134nf_ixM,              144
.equiv          nb134nf_iyM,              160
.equiv          nb134nf_izM,              176
.equiv          nb134nf_jxO,              192
.equiv          nb134nf_jyO,              208
.equiv          nb134nf_jzO,              224
.equiv          nb134nf_jxH1,             240
.equiv          nb134nf_jyH1,             256
.equiv          nb134nf_jzH1,             272
.equiv          nb134nf_jxH2,             288
.equiv          nb134nf_jyH2,             304
.equiv          nb134nf_jzH2,             320
.equiv          nb134nf_jxM,              336
.equiv          nb134nf_jyM,              352
.equiv          nb134nf_jzM,              368
.equiv          nb134nf_qqMM,             384
.equiv          nb134nf_qqMH,             400
.equiv          nb134nf_qqHH,             416
.equiv          nb134nf_two,              432
.equiv          nb134nf_c6,               448
.equiv          nb134nf_c12,              464
.equiv          nb134nf_vctot,            480
.equiv          nb134nf_Vvdwtot,          496
.equiv          nb134nf_half,             512
.equiv          nb134nf_three,            528
.equiv          nb134nf_tsc,              544
.equiv          nb134nf_rsqOO,            560
.equiv          nb134nf_rsqH1H1,          576
.equiv          nb134nf_rsqH1H2,          592
.equiv          nb134nf_rsqH1M,           608
.equiv          nb134nf_rsqH2H1,          624
.equiv          nb134nf_rsqH2H2,          640
.equiv          nb134nf_rsqH2M,           656
.equiv          nb134nf_rsqMH1,           672
.equiv          nb134nf_rsqMH2,           688
.equiv          nb134nf_rsqMM,            704
.equiv          nb134nf_rinvOO,           720
.equiv          nb134nf_rinvH1H1,         736
.equiv          nb134nf_rinvH1H2,         752
.equiv          nb134nf_rinvH1M,          768
.equiv          nb134nf_rinvH2H1,         784
.equiv          nb134nf_rinvH2H2,         800
.equiv          nb134nf_rinvH2M,          816
.equiv          nb134nf_rinvMH1,          832
.equiv          nb134nf_rinvMH2,          848
.equiv          nb134nf_rinvMM,           864
.equiv          nb134nf_is3,              912
.equiv          nb134nf_ii3,              916
.equiv          nb134nf_innerjjnr,        920
.equiv          nb134nf_innerk,           924
.equiv          nb134nf_n,                928
.equiv          nb134nf_nn1,              932
.equiv          nb134nf_nri,              936
.equiv          nb134nf_nouter,           940
.equiv          nb134nf_ninner,           944
.equiv          nb134nf_salign,           948
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 952		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb134nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb134nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb134nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb134nf_nouter], eax
	mov [esp + nb134nf_ninner], eax

	mov eax, [ebp + nb134nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb134nf_tsc], xmm3


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb134nf_half], eax
	mov [esp + nb134nf_half+4], ebx
	movsd xmm1, [esp + nb134nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb134nf_half], xmm1
	movapd [esp + nb134nf_two], xmm2
	movapd [esp + nb134nf_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb134nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb134nf_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb134nf_p_facel]
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
	movapd [esp + nb134nf_qqMM], xmm3
	movapd [esp + nb134nf_qqMH], xmm4
	movapd [esp + nb134nf_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb134nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb134nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb134nf_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movhpd xmm0, [eax + edx*8 + 8]
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb134nf_c6], xmm0
	movapd [esp + nb134nf_c12], xmm1

.nb134nf_threadloop:
        mov   esi, [ebp + nb134nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb134nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb134nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb134nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb134nf_n], eax
        mov [esp + nb134nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb134nf_outerstart
        jmp .nb134nf_end

.nb134nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb134nf_nouter]
	mov [esp + nb134nf_nouter], ebx

.nb134nf_outer:
	mov   eax, [ebp + nb134nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb134nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb134nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb134nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb134nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb134nf_ii3], ebx		

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
	movapd [esp + nb134nf_ixO], xmm3
	movapd [esp + nb134nf_iyO], xmm4
	movapd [esp + nb134nf_izO], xmm5
	movapd [esp + nb134nf_ixH1], xmm6
	movapd [esp + nb134nf_iyH1], xmm7

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
	movapd [esp + nb134nf_izH1], xmm6
	movapd [esp + nb134nf_ixH2], xmm0
	movapd [esp + nb134nf_iyH2], xmm1
	movapd [esp + nb134nf_izH2], xmm2
	movapd [esp + nb134nf_ixM], xmm3
	movapd [esp + nb134nf_iyM], xmm4
	movapd [esp + nb134nf_izM], xmm5

	;# clear vctot
	xorpd xmm4, xmm4
	movapd [esp + nb134nf_vctot], xmm4
	movapd [esp + nb134nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb134nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb134nf_pos] 
	mov   edi, [ebp + nb134nf_faction]	
	mov   eax, [ebp + nb134nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb134nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb134nf_ninner]
	mov   [esp + nb134nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb134nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb134nf_unroll_loop
	jmp   .nb134nf_checksingle
.nb134nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb134nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb134nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb134nf_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb134nf_jxO], xmm0
	movapd 	[esp + nb134nf_jyO], xmm1
	movapd 	[esp + nb134nf_jzO], xmm3
	movapd 	[esp + nb134nf_jxH1], xmm4
	
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
	movapd 	[esp + nb134nf_jyH1], xmm0
	movapd 	[esp + nb134nf_jzH1], xmm1
	movapd 	[esp + nb134nf_jxH2], xmm3
	movapd 	[esp + nb134nf_jyH2], xmm4
	
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
	movapd 	[esp + nb134nf_jzH2], xmm0
	movapd 	[esp + nb134nf_jxM], xmm1
	movapd 	[esp + nb134nf_jyM], xmm3
	movapd 	[esp + nb134nf_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [esp + nb134nf_ixO]
	movapd xmm1, [esp + nb134nf_iyO]
	movapd xmm2, [esp + nb134nf_izO]
	movapd xmm3, [esp + nb134nf_ixH1]
	movapd xmm4, [esp + nb134nf_iyH1]
	movapd xmm5, [esp + nb134nf_izH1]
	subpd  xmm0, [esp + nb134nf_jxO]
	subpd  xmm1, [esp + nb134nf_jyO]
	subpd  xmm2, [esp + nb134nf_jzO]
	subpd  xmm3, [esp + nb134nf_jxH1]
	subpd  xmm4, [esp + nb134nf_jyH1]
	subpd  xmm5, [esp + nb134nf_jzH1]
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
	movapd [esp + nb134nf_rsqOO], xmm0
	movapd [esp + nb134nf_rsqH1H1], xmm3

	movapd xmm0, [esp + nb134nf_ixH1]
	movapd xmm1, [esp + nb134nf_iyH1]
	movapd xmm2, [esp + nb134nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb134nf_jxH2]
	subpd  xmm1, [esp + nb134nf_jyH2]
	subpd  xmm2, [esp + nb134nf_jzH2]
	subpd  xmm3, [esp + nb134nf_jxM]
	subpd  xmm4, [esp + nb134nf_jyM]
	subpd  xmm5, [esp + nb134nf_jzM]
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
	movapd [esp + nb134nf_rsqH1H2], xmm0
	movapd [esp + nb134nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb134nf_ixH2]
	movapd xmm1, [esp + nb134nf_iyH2]
	movapd xmm2, [esp + nb134nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb134nf_jxH1]
	subpd  xmm1, [esp + nb134nf_jyH1]
	subpd  xmm2, [esp + nb134nf_jzH1]
	subpd  xmm3, [esp + nb134nf_jxH2]
	subpd  xmm4, [esp + nb134nf_jyH2]
	subpd  xmm5, [esp + nb134nf_jzH2]
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
	movapd [esp + nb134nf_rsqH2H1], xmm0
	movapd [esp + nb134nf_rsqH2H2], xmm3

	movapd xmm0, [esp + nb134nf_ixH2]
	movapd xmm1, [esp + nb134nf_iyH2]
	movapd xmm2, [esp + nb134nf_izH2]
	movapd xmm3, [esp + nb134nf_ixM]
	movapd xmm4, [esp + nb134nf_iyM]
	movapd xmm5, [esp + nb134nf_izM]
	subpd  xmm0, [esp + nb134nf_jxM]
	subpd  xmm1, [esp + nb134nf_jyM]
	subpd  xmm2, [esp + nb134nf_jzM]
	subpd  xmm3, [esp + nb134nf_jxH1]
	subpd  xmm4, [esp + nb134nf_jyH1]
	subpd  xmm5, [esp + nb134nf_jzH1]
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
	movapd [esp + nb134nf_rsqH2M], xmm0
	movapd [esp + nb134nf_rsqMH1], xmm4

	movapd xmm0, [esp + nb134nf_ixM]
	movapd xmm1, [esp + nb134nf_iyM]
	movapd xmm2, [esp + nb134nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [esp + nb134nf_jxH2]
	subpd  xmm1, [esp + nb134nf_jyH2]
	subpd  xmm2, [esp + nb134nf_jzH2]
	subpd  xmm3, [esp + nb134nf_jxM]
	subpd  xmm4, [esp + nb134nf_jyM]
	subpd  xmm5, [esp + nb134nf_jzM]
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
	movapd [esp + nb134nf_rsqMH2], xmm0
	movapd [esp + nb134nf_rsqMM], xmm4
	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvMH2], xmm1
	movapd [esp + nb134nf_rinvMM], xmm5

	movapd xmm0, [esp + nb134nf_rsqOO]
	movapd xmm4, [esp + nb134nf_rsqH1H1]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb134nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb134nf_half] ;# rinv
	movapd [esp + nb134nf_rinvOO], xmm1
	movapd [esp + nb134nf_rinvH1H1], xmm5

	movapd xmm0, [esp + nb134nf_rsqH1H2]
	movapd xmm4, [esp + nb134nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvH1H2], xmm1
	movapd [esp + nb134nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb134nf_rsqH2H1]
	movapd xmm4, [esp + nb134nf_rsqH2H2]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvH2H1], xmm1
	movapd [esp + nb134nf_rinvH2H2], xmm5

	movapd xmm0, [esp + nb134nf_rsqMH1]
	movapd xmm4, [esp + nb134nf_rsqH2M]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb134nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvMH1], xmm1
	movapd [esp + nb134nf_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [esp + nb134nf_rinvOO] 
	movapd xmm4, [esp + nb134nf_rsqOO]
	
		mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb134nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb134nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb134nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [esp + nb134nf_Vvdwtot]
	movapd [esp + nb134nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb134nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb134nf_Vvdwtot]
	movapd [esp + nb134nf_Vvdwtot], xmm5

	;# All Coulomb interactions
	movapd xmm0, [esp + nb134nf_rinvH1H1]
	movapd xmm1, [esp + nb134nf_rinvH1M]
	addpd  xmm0, [esp + nb134nf_rinvH1H2]
	addpd  xmm1, [esp + nb134nf_rinvH2M]
	addpd  xmm0, [esp + nb134nf_rinvH2H1]
	addpd  xmm1, [esp + nb134nf_rinvMH1]
	addpd  xmm0, [esp + nb134nf_rinvH2H2]
	addpd  xmm1, [esp + nb134nf_rinvMH2]
	movapd xmm2, [esp + nb134nf_rinvMM]
	
	mulpd  xmm0, [esp + nb134nf_qqHH]
	mulpd  xmm1, [esp + nb134nf_qqMH]
	mulpd  xmm2, [esp + nb134nf_qqMM]
	addpd  xmm0, xmm1
	addpd  xmm2, [esp + nb134nf_vctot]
	addpd  xmm0, xmm2
	movapd [esp + nb134nf_vctot], xmm0

	;# should we do one more iteration? 
	sub dword ptr [esp + nb134nf_innerk],  2
	jl    .nb134nf_checksingle
	jmp   .nb134nf_unroll_loop
.nb134nf_checksingle:
	mov   edx, [esp + nb134nf_innerk]
	and   edx, 1
	jnz   .nb134nf_dosingle
	jmp   .nb134nf_updateouterdata
.nb134nf_dosingle:
	mov   edx, [esp + nb134nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb134nf_pos]       ;# base of pos[] 

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
	movsd  [esp + nb134nf_jxO], xmm0
	movsd  [esp + nb134nf_jzO], xmm1
	movsd  [esp + nb134nf_jyH1], xmm2
	movsd  [esp + nb134nf_jxH2], xmm3
	movsd  [esp + nb134nf_jzH2], xmm4
	movsd  [esp + nb134nf_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [esp + nb134nf_jyO], xmm0
	movsd  [esp + nb134nf_jxH1], xmm1
	movsd  [esp + nb134nf_jzH1], xmm2
	movsd  [esp + nb134nf_jyH2], xmm3
	movsd  [esp + nb134nf_jxM], xmm4
	movsd  [esp + nb134nf_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [esp + nb134nf_ixO]
	movapd xmm1, [esp + nb134nf_iyO]
	movapd xmm2, [esp + nb134nf_izO]
	movapd xmm3, [esp + nb134nf_ixH1]
	movapd xmm4, [esp + nb134nf_iyH1]
	movapd xmm5, [esp + nb134nf_izH1]
	subsd  xmm0, [esp + nb134nf_jxO]
	subsd  xmm1, [esp + nb134nf_jyO]
	subsd  xmm2, [esp + nb134nf_jzO]
	subsd  xmm3, [esp + nb134nf_jxH1]
	subsd  xmm4, [esp + nb134nf_jyH1]
	subsd  xmm5, [esp + nb134nf_jzH1]
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
	movapd [esp + nb134nf_rsqOO], xmm0
	movapd [esp + nb134nf_rsqH1H1], xmm3

	movapd xmm0, [esp + nb134nf_ixH1]
	movapd xmm1, [esp + nb134nf_iyH1]
	movapd xmm2, [esp + nb134nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb134nf_jxH2]
	subsd  xmm1, [esp + nb134nf_jyH2]
	subsd  xmm2, [esp + nb134nf_jzH2]
	subsd  xmm3, [esp + nb134nf_jxM]
	subsd  xmm4, [esp + nb134nf_jyM]
	subsd  xmm5, [esp + nb134nf_jzM]
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
	movapd [esp + nb134nf_rsqH1H2], xmm0
	movapd [esp + nb134nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb134nf_ixH2]
	movapd xmm1, [esp + nb134nf_iyH2]
	movapd xmm2, [esp + nb134nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb134nf_jxH1]
	subsd  xmm1, [esp + nb134nf_jyH1]
	subsd  xmm2, [esp + nb134nf_jzH1]
	subsd  xmm3, [esp + nb134nf_jxH2]
	subsd  xmm4, [esp + nb134nf_jyH2]
	subsd  xmm5, [esp + nb134nf_jzH2]
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
	movapd [esp + nb134nf_rsqH2H1], xmm0
	movapd [esp + nb134nf_rsqH2H2], xmm3

	movapd xmm0, [esp + nb134nf_ixH2]
	movapd xmm1, [esp + nb134nf_iyH2]
	movapd xmm2, [esp + nb134nf_izH2]
	movapd xmm3, [esp + nb134nf_ixM]
	movapd xmm4, [esp + nb134nf_iyM]
	movapd xmm5, [esp + nb134nf_izM]
	subsd  xmm0, [esp + nb134nf_jxM]
	subsd  xmm1, [esp + nb134nf_jyM]
	subsd  xmm2, [esp + nb134nf_jzM]
	subsd  xmm3, [esp + nb134nf_jxH1]
	subsd  xmm4, [esp + nb134nf_jyH1]
	subsd  xmm5, [esp + nb134nf_jzH1]
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
	movapd [esp + nb134nf_rsqH2M], xmm0
	movapd [esp + nb134nf_rsqMH1], xmm4

	movapd xmm0, [esp + nb134nf_ixM]
	movapd xmm1, [esp + nb134nf_iyM]
	movapd xmm2, [esp + nb134nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [esp + nb134nf_jxH2]
	subsd  xmm1, [esp + nb134nf_jyH2]
	subsd  xmm2, [esp + nb134nf_jzH2]
	subsd  xmm3, [esp + nb134nf_jxM]
	subsd  xmm4, [esp + nb134nf_jyM]
	subsd  xmm5, [esp + nb134nf_jzM]
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
	movapd [esp + nb134nf_rsqMH2], xmm0
	movapd [esp + nb134nf_rsqMM], xmm4

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
	movapd  xmm3, [esp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvMH2], xmm1
	movapd [esp + nb134nf_rinvMM], xmm5

	movapd xmm0, [esp + nb134nf_rsqOO]
	movapd xmm4, [esp + nb134nf_rsqH1H1]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb134nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb134nf_half] ;# rinv
	movapd [esp + nb134nf_rinvOO], xmm1
	movapd [esp + nb134nf_rinvH1H1], xmm5

	movapd xmm0, [esp + nb134nf_rsqH1H2]
	movapd xmm4, [esp + nb134nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvH1H2], xmm1
	movapd [esp + nb134nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb134nf_rsqH2H1]
	movapd xmm4, [esp + nb134nf_rsqH2H2]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvH2H1], xmm1
	movapd [esp + nb134nf_rinvH2H2], xmm5

	movapd xmm0, [esp + nb134nf_rsqMH1]
	movapd xmm4, [esp + nb134nf_rsqH2M]	
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
	movapd  xmm3, [esp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb134nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb134nf_half] ;# rinv 
	movapd [esp + nb134nf_rinvMH1], xmm1
	movapd [esp + nb134nf_rinvH2M], xmm5

	;# start with OO interaction 
	movsd xmm0, [esp + nb134nf_rinvOO] 
	movsd xmm4, [esp + nb134nf_rsqOO]
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb134nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb134nf_VFtab]

	;# dispersion 
	movlpd xmm4, [esi + ebx*8]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 16]	;# G1
	movhpd xmm6, [esi + ebx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb134nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [esp + nb134nf_Vvdwtot]
	movsd [esp + nb134nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + ebx*8 + 32]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 48]	;# G1
	movhpd xmm6, [esi + ebx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb134nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb134nf_Vvdwtot]
	movsd [esp + nb134nf_Vvdwtot], xmm5

	;# All Coulomb interactions
	movsd xmm0, [esp + nb134nf_rinvH1H1]
	movsd xmm1, [esp + nb134nf_rinvH1M]
	addsd  xmm0, [esp + nb134nf_rinvH1H2]
	addsd  xmm1, [esp + nb134nf_rinvH2M]
	addsd  xmm0, [esp + nb134nf_rinvH2H1]
	addsd  xmm1, [esp + nb134nf_rinvMH1]
	addsd  xmm0, [esp + nb134nf_rinvH2H2]
	addsd  xmm1, [esp + nb134nf_rinvMH2]
	movsd xmm2, [esp + nb134nf_rinvMM]
	
	mulsd  xmm0, [esp + nb134nf_qqHH]
	mulsd  xmm1, [esp + nb134nf_qqMH]
	mulsd  xmm2, [esp + nb134nf_qqMM]
	addsd  xmm0, xmm1
	addsd  xmm2, [esp + nb134nf_vctot]
	addsd  xmm0, xmm2
	movsd [esp + nb134nf_vctot], xmm0
	
.nb134nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb134nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb134nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb134nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb134nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb134nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb134nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb134nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb134nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb134nf_n], esi
        jmp .nb134nf_outer
.nb134nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb134nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb134nf_end
        ;# non-zero, do one more workunit
        jmp   .nb134nf_threadloop
.nb134nf_end:
	emms

	mov eax, [esp + nb134nf_nouter]
	mov ebx, [esp + nb134nf_ninner]
	mov ecx, [ebp + nb134nf_outeriter]
	mov edx, [ebp + nb134nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb134nf_salign]
	add esp, eax
	add esp, 952
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


