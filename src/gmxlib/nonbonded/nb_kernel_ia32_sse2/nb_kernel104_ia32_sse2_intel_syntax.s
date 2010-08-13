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



.globl nb_kernel104_ia32_sse2
.globl _nb_kernel104_ia32_sse2
nb_kernel104_ia32_sse2:	
_nb_kernel104_ia32_sse2:	
.equiv          nb104_p_nri,            8
.equiv          nb104_iinr,             12
.equiv          nb104_jindex,           16
.equiv          nb104_jjnr,             20
.equiv          nb104_shift,            24
.equiv          nb104_shiftvec,         28
.equiv          nb104_fshift,           32
.equiv          nb104_gid,              36
.equiv          nb104_pos,              40
.equiv          nb104_faction,          44
.equiv          nb104_charge,           48
.equiv          nb104_p_facel,          52
.equiv          nb104_argkrf,           56
.equiv          nb104_argcrf,           60
.equiv          nb104_Vc,               64
.equiv          nb104_type,             68
.equiv          nb104_p_ntype,          72
.equiv          nb104_vdwparam,         76
.equiv          nb104_Vvdw,             80
.equiv          nb104_p_tabscale,       84
.equiv          nb104_VFtab,            88
.equiv          nb104_invsqrta,         92
.equiv          nb104_dvda,             96
.equiv          nb104_p_gbtabscale,     100
.equiv          nb104_GBtab,            104
.equiv          nb104_p_nthreads,       108
.equiv          nb104_count,            112
.equiv          nb104_mtx,              116
.equiv          nb104_outeriter,        120
.equiv          nb104_inneriter,        124
.equiv          nb104_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 	
.equiv          nb104_ixM,              0
.equiv          nb104_iyM,              16
.equiv          nb104_izM,              32
.equiv          nb104_ixH1,             48
.equiv          nb104_iyH1,             64
.equiv          nb104_izH1,             80
.equiv          nb104_ixH2,             96
.equiv          nb104_iyH2,             112
.equiv          nb104_izH2,             128
.equiv          nb104_jxM,              144
.equiv          nb104_jyM,              160
.equiv          nb104_jzM,              176
.equiv          nb104_jxH1,             192
.equiv          nb104_jyH1,             208
.equiv          nb104_jzH1,             224
.equiv          nb104_jxH2,             240
.equiv          nb104_jyH2,             256
.equiv          nb104_jzH2,             272
.equiv          nb104_dxMM,             288
.equiv          nb104_dyMM,             304
.equiv          nb104_dzMM,             320
.equiv          nb104_dxMH1,            336
.equiv          nb104_dyMH1,            352
.equiv          nb104_dzMH1,            368
.equiv          nb104_dxMH2,            384
.equiv          nb104_dyMH2,            400
.equiv          nb104_dzMH2,            416
.equiv          nb104_dxH1M,            432
.equiv          nb104_dyH1M,            448
.equiv          nb104_dzH1M,            464
.equiv          nb104_dxH1H1,           480
.equiv          nb104_dyH1H1,           496
.equiv          nb104_dzH1H1,           512
.equiv          nb104_dxH1H2,           528
.equiv          nb104_dyH1H2,           544
.equiv          nb104_dzH1H2,           560
.equiv          nb104_dxH2M,            576
.equiv          nb104_dyH2M,            592
.equiv          nb104_dzH2M,            608
.equiv          nb104_dxH2H1,           624
.equiv          nb104_dyH2H1,           640
.equiv          nb104_dzH2H1,           656
.equiv          nb104_dxH2H2,           672
.equiv          nb104_dyH2H2,           688
.equiv          nb104_dzH2H2,           704
.equiv          nb104_qqMM,             720
.equiv          nb104_qqMH,             736
.equiv          nb104_qqHH,             752
.equiv          nb104_vctot,            768
.equiv          nb104_fixM,             784
.equiv          nb104_fiyM,             800
.equiv          nb104_fizM,             816
.equiv          nb104_fixH1,            832
.equiv          nb104_fiyH1,            848
.equiv          nb104_fizH1,            864
.equiv          nb104_fixH2,            880
.equiv          nb104_fiyH2,            896
.equiv          nb104_fizH2,            912
.equiv          nb104_fjxM,             928
.equiv          nb104_fjyM,             944
.equiv          nb104_fjzM,             960
.equiv          nb104_fjxH1,            976
.equiv          nb104_fjyH1,            992
.equiv          nb104_fjzH1,            1008
.equiv          nb104_fjxH2,            1024
.equiv          nb104_fjyH2,            1040
.equiv          nb104_fjzH2,            1056
.equiv          nb104_half,             1072
.equiv          nb104_three,            1088
.equiv          nb104_rsqMM,            1104
.equiv          nb104_rsqMH1,           1120
.equiv          nb104_rsqMH2,           1136
.equiv          nb104_rsqH1M,           1152
.equiv          nb104_rsqH1H1,          1168
.equiv          nb104_rsqH1H2,          1184
.equiv          nb104_rsqH2M,           1200
.equiv          nb104_rsqH2H1,          1216
.equiv          nb104_rsqH2H2,          1232
.equiv          nb104_rinvMM,           1248
.equiv          nb104_rinvMH1,          1264
.equiv          nb104_rinvMH2,          1280
.equiv          nb104_rinvH1M,          1296
.equiv          nb104_rinvH1H1,         1312
.equiv          nb104_rinvH1H2,         1328
.equiv          nb104_rinvH2M,          1344
.equiv          nb104_rinvH2H1,         1360
.equiv          nb104_rinvH2H2,         1376
.equiv          nb104_is3,              1392
.equiv          nb104_ii3,              1396
.equiv          nb104_innerjjnr,        1400
.equiv          nb104_innerk,           1404
.equiv          nb104_n,                1408
.equiv          nb104_nn1,              1412
.equiv          nb104_nri,              1416
.equiv          nb104_nouter,           1420
.equiv          nb104_ninner,           1424
.equiv          nb104_salign,           1428
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 1432		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb104_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb104_p_nri]
	mov ecx, [ecx]
	mov [esp + nb104_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb104_nouter], eax
	mov [esp + nb104_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb104_half], eax
	mov [esp + nb104_half+4], ebx
	movsd xmm1, [esp + nb104_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb104_half], xmm1
	movapd [esp + nb104_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb104_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb104_charge]
	movsd xmm3, [edx + ebx*8 + 24]	;# qM 
	movsd xmm4, xmm3		;# qM 
	movsd xmm5, [edx + ebx*8 + 8]	;# qH 
	mov esi, [ebp + nb104_p_facel]
	movsd xmm6, [esi]	;# facel 
	mulsd  xmm3, xmm3		;# qM*qM 
	mulsd  xmm4, xmm5		;# qM*qH 
	mulsd  xmm5, xmm5		;# qH*qH 
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb104_qqMM], xmm3
	movapd [esp + nb104_qqMH], xmm4
	movapd [esp + nb104_qqHH], xmm5

.nb104_threadloop:
        mov   esi, [ebp + nb104_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb104_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb104_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb104_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb104_n], eax
        mov [esp + nb104_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb104_outerstart
        jmp .nb104_end

.nb104_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb104_nouter]
	mov [esp + nb104_nouter], ebx

.nb104_outer:
	mov   eax, [ebp + nb104_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb104_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb104_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb104_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb104_pos]    ;# eax = base of pos[]  
	mov   [esp + nb104_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb104_ixH1], xmm3
	movapd [esp + nb104_iyH1], xmm4
	movapd [esp + nb104_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 48]
	addsd xmm1, [eax + ebx*8 + 56]
	addsd xmm2, [eax + ebx*8 + 64]		
	addsd xmm3, [eax + ebx*8 + 72]
	addsd xmm4, [eax + ebx*8 + 80]
	addsd xmm5, [eax + ebx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb104_ixH2], xmm0
	movapd [esp + nb104_iyH2], xmm1
	movapd [esp + nb104_izH2], xmm2
	movapd [esp + nb104_ixM], xmm3
	movapd [esp + nb104_iyM], xmm4
	movapd [esp + nb104_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb104_vctot], xmm4
	movapd [esp + nb104_fixM], xmm4
	movapd [esp + nb104_fiyM], xmm4
	movapd [esp + nb104_fizM], xmm4
	movapd [esp + nb104_fixH1], xmm4
	movapd [esp + nb104_fiyH1], xmm4
	movapd [esp + nb104_fizH1], xmm4
	movapd [esp + nb104_fixH2], xmm4
	movapd [esp + nb104_fiyH2], xmm4
	movapd [esp + nb104_fizH2], xmm4
	
	mov   eax, [ebp + nb104_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax +esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb104_pos]
	mov   edi, [ebp + nb104_faction]	
	mov   eax, [ebp + nb104_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb104_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb104_ninner]
	mov   [esp + nb104_ninner], ecx
	add   edx, 0
	mov   [esp + nb104_innerk], edx    ;# number of innerloop atoms 
	jge   .nb104_unroll_loop
	jmp   .nb104_checksingle
.nb104_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb104_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb104_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb104_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 24]
	movhpd xmm3, [esi + ebx*8 + 32]
	movhpd xmm4, [esi + ebx*8 + 40]
	movhpd xmm5, [esi + ebx*8 + 48]
	movhpd xmm6, [esi + ebx*8 + 56]
	movhpd xmm7, [esi + ebx*8 + 64]
	movapd 	[esp + nb104_jxH1], xmm2
	movapd 	[esp + nb104_jyH1], xmm3
	movapd 	[esp + nb104_jzH1], xmm4
	movapd 	[esp + nb104_jxH2], xmm5
	movapd 	[esp + nb104_jyH2], xmm6
	movapd 	[esp + nb104_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movhpd xmm2, [esi + ebx*8 + 72]
	movhpd xmm3, [esi + ebx*8 + 80]
	movhpd xmm4, [esi + ebx*8 + 88]
	movapd 	[esp + nb104_jxM], xmm2
	movapd 	[esp + nb104_jyM], xmm3
	movapd 	[esp + nb104_jzM], xmm4
	
	movapd xmm0, [esp + nb104_ixM]
	movapd xmm1, [esp + nb104_iyM]
	movapd xmm2, [esp + nb104_izM]
	movapd xmm3, [esp + nb104_ixM]
	movapd xmm4, [esp + nb104_iyM]
	movapd xmm5, [esp + nb104_izM]
	subpd  xmm0, [esp + nb104_jxM]
	subpd  xmm1, [esp + nb104_jyM]
	subpd  xmm2, [esp + nb104_jzM]
	subpd  xmm3, [esp + nb104_jxH1]
	subpd  xmm4, [esp + nb104_jyH1]
	subpd  xmm5, [esp + nb104_jzH1]
	movapd [esp + nb104_dxMM], xmm0
	movapd [esp + nb104_dyMM], xmm1
	movapd [esp + nb104_dzMM], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb104_dxMH1], xmm3
	movapd [esp + nb104_dyMH1], xmm4
	movapd [esp + nb104_dzMH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb104_rsqMM], xmm0
	movapd [esp + nb104_rsqMH1], xmm3

	movapd xmm0, [esp + nb104_ixM]
	movapd xmm1, [esp + nb104_iyM]
	movapd xmm2, [esp + nb104_izM]
	movapd xmm3, [esp + nb104_ixH1]
	movapd xmm4, [esp + nb104_iyH1]
	movapd xmm5, [esp + nb104_izH1]
	subpd  xmm0, [esp + nb104_jxH2]
	subpd  xmm1, [esp + nb104_jyH2]
	subpd  xmm2, [esp + nb104_jzH2]
	subpd  xmm3, [esp + nb104_jxM]
	subpd  xmm4, [esp + nb104_jyM]
	subpd  xmm5, [esp + nb104_jzM]
	movapd [esp + nb104_dxMH2], xmm0
	movapd [esp + nb104_dyMH2], xmm1
	movapd [esp + nb104_dzMH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb104_dxH1M], xmm3
	movapd [esp + nb104_dyH1M], xmm4
	movapd [esp + nb104_dzH1M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb104_rsqMH2], xmm0
	movapd [esp + nb104_rsqH1M], xmm3

	movapd xmm0, [esp + nb104_ixH1]
	movapd xmm1, [esp + nb104_iyH1]
	movapd xmm2, [esp + nb104_izH1]
	movapd xmm3, [esp + nb104_ixH1]
	movapd xmm4, [esp + nb104_iyH1]
	movapd xmm5, [esp + nb104_izH1]
	subpd  xmm0, [esp + nb104_jxH1]
	subpd  xmm1, [esp + nb104_jyH1]
	subpd  xmm2, [esp + nb104_jzH1]
	subpd  xmm3, [esp + nb104_jxH2]
	subpd  xmm4, [esp + nb104_jyH2]
	subpd  xmm5, [esp + nb104_jzH2]
	movapd [esp + nb104_dxH1H1], xmm0
	movapd [esp + nb104_dyH1H1], xmm1
	movapd [esp + nb104_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb104_dxH1H2], xmm3
	movapd [esp + nb104_dyH1H2], xmm4
	movapd [esp + nb104_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb104_rsqH1H1], xmm0
	movapd [esp + nb104_rsqH1H2], xmm3

	movapd xmm0, [esp + nb104_ixH2]
	movapd xmm1, [esp + nb104_iyH2]
	movapd xmm2, [esp + nb104_izH2]
	movapd xmm3, [esp + nb104_ixH2]
	movapd xmm4, [esp + nb104_iyH2]
	movapd xmm5, [esp + nb104_izH2]
	subpd  xmm0, [esp + nb104_jxM]
	subpd  xmm1, [esp + nb104_jyM]
	subpd  xmm2, [esp + nb104_jzM]
	subpd  xmm3, [esp + nb104_jxH1]
	subpd  xmm4, [esp + nb104_jyH1]
	subpd  xmm5, [esp + nb104_jzH1]
	movapd [esp + nb104_dxH2M], xmm0
	movapd [esp + nb104_dyH2M], xmm1
	movapd [esp + nb104_dzH2M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb104_dxH2H1], xmm3
	movapd [esp + nb104_dyH2H1], xmm4
	movapd [esp + nb104_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb104_rsqH2M], xmm0
	movapd [esp + nb104_rsqH2H1], xmm4

	movapd xmm0, [esp + nb104_ixH2]
	movapd xmm1, [esp + nb104_iyH2]
	movapd xmm2, [esp + nb104_izH2]
	subpd  xmm0, [esp + nb104_jxH2]
	subpd  xmm1, [esp + nb104_jyH2]
	subpd  xmm2, [esp + nb104_jzH2]
	movapd [esp + nb104_dxH2H2], xmm0
	movapd [esp + nb104_dyH2H2], xmm1
	movapd [esp + nb104_dzH2H2], xmm2
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb104_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0 (h2h2) , xmm4 (h2h1) 
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
	movapd  xmm3, [esp + nb104_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104_half] ;# iter1 
	mulpd   xmm7, [esp + nb104_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104_half] ;# rinv 
	mulpd   xmm5, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvH2H2], xmm1
	movapd [esp + nb104_rinvH2H1], xmm5

	movapd xmm0, [esp + nb104_rsqMM]
	movapd xmm4, [esp + nb104_rsqMH1]	
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
	movapd  xmm3, [esp + nb104_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb104_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104_half] ;# rinv 
	mulpd   xmm5, [esp + nb104_half] ;# rinv
	movapd [esp + nb104_rinvMM], xmm1
	movapd [esp + nb104_rinvMH1], xmm5

	movapd xmm0, [esp + nb104_rsqMH2]
	movapd xmm4, [esp + nb104_rsqH1M]	
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
	movapd  xmm3, [esp + nb104_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104_half] ;# iter1 
	mulpd   xmm7, [esp + nb104_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104_half] ;# rinv 
	mulpd   xmm5, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvMH2], xmm1
	movapd [esp + nb104_rinvH1M], xmm5

	movapd xmm0, [esp + nb104_rsqH1H1]
	movapd xmm4, [esp + nb104_rsqH1H2]	
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
	movapd  xmm3, [esp + nb104_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104_half] ;# iter1a 
	mulpd   xmm7, [esp + nb104_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104_half] ;# rinv 
	mulpd   xmm5, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvH1H1], xmm1
	movapd [esp + nb104_rinvH1H2], xmm5

	movapd xmm0, [esp + nb104_rsqH2M]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb104_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb104_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb104_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvH2M], xmm1

	;# start with MM interaction 
	movapd xmm0, [esp + nb104_rinvMM]
	movapd xmm7, xmm0
	mulpd  xmm0, xmm0		;# rinvsq 
	mulpd  xmm7, [esp + nb104_qqMM]	
	mulpd  xmm0, xmm7	
	addpd  xmm7, [esp + nb104_vctot] 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb104_dxMM]
	mulpd xmm1, [esp + nb104_dyMM]
	mulpd xmm2, [esp + nb104_dzMM]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixM]
	addpd xmm1, [esp + nb104_fiyM]
	addpd xmm2, [esp + nb104_fizM]
	movapd [esp + nb104_fjxM], xmm3
	movapd [esp + nb104_fjyM], xmm4
	movapd [esp + nb104_fjzM], xmm5
	movapd [esp + nb104_fixM], xmm0
	movapd [esp + nb104_fiyM], xmm1
	movapd [esp + nb104_fizM], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb104_rinvMH1]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqMH]
	mulpd xmm0, xmm1	;# fsMH1  
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb104_dxMH1]
	mulpd xmm1, [esp + nb104_dyMH1]
	mulpd xmm2, [esp + nb104_dzMH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixM]
	addpd xmm1, [esp + nb104_fiyM]
	addpd xmm2, [esp + nb104_fizM]
	movapd [esp + nb104_fjxH1], xmm3
	movapd [esp + nb104_fjyH1], xmm4
	movapd [esp + nb104_fjzH1], xmm5
	movapd [esp + nb104_fixM], xmm0
	movapd [esp + nb104_fiyM], xmm1
	movapd [esp + nb104_fizM], xmm2

	;# M-H2 interaction  
	movapd xmm0, [esp + nb104_rinvMH2]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqMH]
	mulpd xmm0, xmm1	;# fsMH2  
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb104_dxMH2]
	mulpd xmm1, [esp + nb104_dyMH2]
	mulpd xmm2, [esp + nb104_dzMH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixM]
	addpd xmm1, [esp + nb104_fiyM]
	addpd xmm2, [esp + nb104_fizM]
	movapd [esp + nb104_fjxH2], xmm3
	movapd [esp + nb104_fjyH2], xmm4
	movapd [esp + nb104_fjzH2], xmm5
	movapd [esp + nb104_fixM], xmm0
	movapd [esp + nb104_fiyM], xmm1
	movapd [esp + nb104_fizM], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb104_rinvH1M]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqMH]
	mulpd xmm0, xmm1	;# fsH1M 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxM]
	movapd xmm4, [esp + nb104_fjyM]
	movapd xmm5, [esp + nb104_fjzM]
	mulpd xmm0, [esp + nb104_dxH1M]
	mulpd xmm1, [esp + nb104_dyH1M]
	mulpd xmm2, [esp + nb104_dzH1M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixH1]
	addpd xmm1, [esp + nb104_fiyH1]
	addpd xmm2, [esp + nb104_fizH1]
	movapd [esp + nb104_fjxM], xmm3
	movapd [esp + nb104_fjyM], xmm4
	movapd [esp + nb104_fjzM], xmm5
	movapd [esp + nb104_fixH1], xmm0
	movapd [esp + nb104_fiyH1], xmm1
	movapd [esp + nb104_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb104_rinvH1H1]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqHH]
	mulpd xmm0, xmm1	;# fsH1H1 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH1]
	movapd xmm4, [esp + nb104_fjyH1]
	movapd xmm5, [esp + nb104_fjzH1]
	mulpd xmm0, [esp + nb104_dxH1H1]
	mulpd xmm1, [esp + nb104_dyH1H1]
	mulpd xmm2, [esp + nb104_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixH1]
	addpd xmm1, [esp + nb104_fiyH1]
	addpd xmm2, [esp + nb104_fizH1]
	movapd [esp + nb104_fjxH1], xmm3
	movapd [esp + nb104_fjyH1], xmm4
	movapd [esp + nb104_fjzH1], xmm5
	movapd [esp + nb104_fixH1], xmm0
	movapd [esp + nb104_fiyH1], xmm1
	movapd [esp + nb104_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb104_rinvH1H2]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqHH]
	mulpd xmm0, xmm1	;# fsMH2  
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH2]
	movapd xmm4, [esp + nb104_fjyH2]
	movapd xmm5, [esp + nb104_fjzH2]
	mulpd xmm0, [esp + nb104_dxH1H2]
	mulpd xmm1, [esp + nb104_dyH1H2]
	mulpd xmm2, [esp + nb104_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixH1]
	addpd xmm1, [esp + nb104_fiyH1]
	addpd xmm2, [esp + nb104_fizH1]
	movapd [esp + nb104_fjxH2], xmm3
	movapd [esp + nb104_fjyH2], xmm4
	movapd [esp + nb104_fjzH2], xmm5
	movapd [esp + nb104_fixH1], xmm0
	movapd [esp + nb104_fiyH1], xmm1
	movapd [esp + nb104_fizH1], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb104_rinvH2M]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqMH]
	mulpd xmm0, xmm1	;# fsH2M 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxM]
	movapd xmm4, [esp + nb104_fjyM]
	movapd xmm5, [esp + nb104_fjzM]
	mulpd xmm0, [esp + nb104_dxH2M]
	mulpd xmm1, [esp + nb104_dyH2M]
	mulpd xmm2, [esp + nb104_dzH2M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixH2]
	addpd xmm1, [esp + nb104_fiyH2]
	addpd xmm2, [esp + nb104_fizH2]
	movapd [esp + nb104_fjxM], xmm3
	movapd [esp + nb104_fjyM], xmm4
	movapd [esp + nb104_fjzM], xmm5
	movapd [esp + nb104_fixH2], xmm0
	movapd [esp + nb104_fiyH2], xmm1
	movapd [esp + nb104_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb104_rinvH2H1]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqHH]
	mulpd xmm0, xmm1	;# fsH2H1 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH1]
	movapd xmm4, [esp + nb104_fjyH1]
	movapd xmm5, [esp + nb104_fjzH1]
	mulpd xmm0, [esp + nb104_dxH2H1]
	mulpd xmm1, [esp + nb104_dyH2H1]
	mulpd xmm2, [esp + nb104_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixH2]
	addpd xmm1, [esp + nb104_fiyH2]
	addpd xmm2, [esp + nb104_fizH2]
	movapd [esp + nb104_fjxH1], xmm3
	movapd [esp + nb104_fjyH1], xmm4
	movapd [esp + nb104_fjzH1], xmm5
	movapd [esp + nb104_fixH2], xmm0
	movapd [esp + nb104_fiyH2], xmm1
	movapd [esp + nb104_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb104_rinvH2H2]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb104_qqHH]
	mulpd xmm0, xmm1	;# fsH2H2 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd [esp + nb104_vctot], xmm7
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH2]
	movapd xmm4, [esp + nb104_fjyH2]
	movapd xmm5, [esp + nb104_fjzH2]
	mulpd xmm0, [esp + nb104_dxH2H2]
	mulpd xmm1, [esp + nb104_dyH2H2]
	mulpd xmm2, [esp + nb104_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb104_fixH2]
	addpd xmm1, [esp + nb104_fiyH2]
	addpd xmm2, [esp + nb104_fizH2]
	movapd [esp + nb104_fjxH2], xmm3
	movapd [esp + nb104_fjyH2], xmm4
	movapd [esp + nb104_fjzH2], xmm5
	movapd [esp + nb104_fixH2], xmm0
	movapd [esp + nb104_fiyH2], xmm1
	movapd [esp + nb104_fizH2], xmm2

	mov edi, [ebp + nb104_faction]
		
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8 + 24]
	movlpd xmm1, [edi + eax*8 + 32]
	movlpd xmm2, [edi + eax*8 + 40]
	movlpd xmm3, [edi + eax*8 + 48]
	movlpd xmm4, [edi + eax*8 + 56]
	movlpd xmm5, [edi + eax*8 + 64]
	movlpd xmm6, [edi + eax*8 + 72]
	movlpd xmm7, [edi + eax*8 + 80]
	movhpd xmm0, [edi + ebx*8 + 24]
	movhpd xmm1, [edi + ebx*8 + 32]
	movhpd xmm2, [edi + ebx*8 + 40]
	movhpd xmm3, [edi + ebx*8 + 48]
	movhpd xmm4, [edi + ebx*8 + 56]
	movhpd xmm5, [edi + ebx*8 + 64]
	movhpd xmm6, [edi + ebx*8 + 72]
	movhpd xmm7, [edi + ebx*8 + 80]
	addpd xmm0, [esp + nb104_fjxH1]
	addpd xmm1, [esp + nb104_fjyH1]
	addpd xmm2, [esp + nb104_fjzH1]
	addpd xmm3, [esp + nb104_fjxH2]
	addpd xmm4, [esp + nb104_fjyH2]
	addpd xmm5, [esp + nb104_fjzH2]
	addpd xmm6, [esp + nb104_fjxM]
	addpd xmm7, [esp + nb104_fjyM]
	movlpd [edi + eax*8 + 24], xmm0
	movlpd [edi + eax*8 + 32], xmm1
	movlpd [edi + eax*8 + 40], xmm2
	movlpd [edi + eax*8 + 48], xmm3
	movlpd [edi + eax*8 + 56], xmm4
	movlpd [edi + eax*8 + 64], xmm5
	movlpd [edi + eax*8 + 72], xmm6
	movlpd [edi + eax*8 + 80], xmm7
	movhpd [edi + ebx*8 + 24], xmm0
	movhpd [edi + ebx*8 + 32], xmm1
	movhpd [edi + ebx*8 + 40], xmm2
	movhpd [edi + ebx*8 + 48], xmm3
	movhpd [edi + ebx*8 + 56], xmm4
	movhpd [edi + ebx*8 + 64], xmm5
	movhpd [edi + ebx*8 + 72], xmm6
	movhpd [edi + ebx*8 + 80], xmm7

	movlpd xmm0, [edi + eax*8 + 88]
	movhpd xmm0, [edi + ebx*8 + 88]
	addpd xmm0, [esp + nb104_fjzM]
	movlpd [edi + eax*8 + 88], xmm0
	movhpd [edi + ebx*8 + 88], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb104_innerk],  2
	jl    .nb104_checksingle
	jmp   .nb104_unroll_loop
.nb104_checksingle:
	mov   edx, [esp + nb104_innerk]
	and   edx, 1
	jnz   .nb104_dosingle
	jmp   .nb104_updateouterdata
.nb104_dosingle:
	mov   edx, [esp + nb104_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	mov esi, [ebp + nb104_pos]
	lea   eax, [eax + eax*2]  

	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movapd 	[esp + nb104_jxH1], xmm2
	movapd 	[esp + nb104_jyH1], xmm3
	movapd 	[esp + nb104_jzH1], xmm4
	movapd 	[esp + nb104_jxH2], xmm5
	movapd 	[esp + nb104_jyH2], xmm6
	movapd 	[esp + nb104_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movapd 	[esp + nb104_jxM], xmm2
	movapd 	[esp + nb104_jyM], xmm3
	movapd 	[esp + nb104_jzM], xmm4
	
	movapd xmm0, [esp + nb104_ixM]
	movapd xmm1, [esp + nb104_iyM]
	movapd xmm2, [esp + nb104_izM]
	movapd xmm3, [esp + nb104_ixM]
	movapd xmm4, [esp + nb104_iyM]
	movapd xmm5, [esp + nb104_izM]
	subsd  xmm0, [esp + nb104_jxM]
	subsd  xmm1, [esp + nb104_jyM]
	subsd  xmm2, [esp + nb104_jzM]
	subsd  xmm3, [esp + nb104_jxH1]
	subsd  xmm4, [esp + nb104_jyH1]
	subsd  xmm5, [esp + nb104_jzH1]
	movapd [esp + nb104_dxMM], xmm0
	movapd [esp + nb104_dyMM], xmm1
	movapd [esp + nb104_dzMM], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb104_dxMH1], xmm3
	movapd [esp + nb104_dyMH1], xmm4
	movapd [esp + nb104_dzMH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb104_rsqMM], xmm0
	movapd [esp + nb104_rsqMH1], xmm3

	movapd xmm0, [esp + nb104_ixM]
	movapd xmm1, [esp + nb104_iyM]
	movapd xmm2, [esp + nb104_izM]
	movapd xmm3, [esp + nb104_ixH1]
	movapd xmm4, [esp + nb104_iyH1]
	movapd xmm5, [esp + nb104_izH1]
	subsd  xmm0, [esp + nb104_jxH2]
	subsd  xmm1, [esp + nb104_jyH2]
	subsd  xmm2, [esp + nb104_jzH2]
	subsd  xmm3, [esp + nb104_jxM]
	subsd  xmm4, [esp + nb104_jyM]
	subsd  xmm5, [esp + nb104_jzM]
	movapd [esp + nb104_dxMH2], xmm0
	movapd [esp + nb104_dyMH2], xmm1
	movapd [esp + nb104_dzMH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb104_dxH1M], xmm3
	movapd [esp + nb104_dyH1M], xmm4
	movapd [esp + nb104_dzH1M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb104_rsqMH2], xmm0
	movapd [esp + nb104_rsqH1M], xmm3

	movapd xmm0, [esp + nb104_ixH1]
	movapd xmm1, [esp + nb104_iyH1]
	movapd xmm2, [esp + nb104_izH1]
	movapd xmm3, [esp + nb104_ixH1]
	movapd xmm4, [esp + nb104_iyH1]
	movapd xmm5, [esp + nb104_izH1]
	subsd  xmm0, [esp + nb104_jxH1]
	subsd  xmm1, [esp + nb104_jyH1]
	subsd  xmm2, [esp + nb104_jzH1]
	subsd  xmm3, [esp + nb104_jxH2]
	subsd  xmm4, [esp + nb104_jyH2]
	subsd  xmm5, [esp + nb104_jzH2]
	movapd [esp + nb104_dxH1H1], xmm0
	movapd [esp + nb104_dyH1H1], xmm1
	movapd [esp + nb104_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb104_dxH1H2], xmm3
	movapd [esp + nb104_dyH1H2], xmm4
	movapd [esp + nb104_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb104_rsqH1H1], xmm0
	movapd [esp + nb104_rsqH1H2], xmm3

	movapd xmm0, [esp + nb104_ixH2]
	movapd xmm1, [esp + nb104_iyH2]
	movapd xmm2, [esp + nb104_izH2]
	movapd xmm3, [esp + nb104_ixH2]
	movapd xmm4, [esp + nb104_iyH2]
	movapd xmm5, [esp + nb104_izH2]
	subsd  xmm0, [esp + nb104_jxM]
	subsd  xmm1, [esp + nb104_jyM]
	subsd  xmm2, [esp + nb104_jzM]
	subsd  xmm3, [esp + nb104_jxH1]
	subsd  xmm4, [esp + nb104_jyH1]
	subsd  xmm5, [esp + nb104_jzH1]
	movapd [esp + nb104_dxH2M], xmm0
	movapd [esp + nb104_dyH2M], xmm1
	movapd [esp + nb104_dzH2M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb104_dxH2H1], xmm3
	movapd [esp + nb104_dyH2H1], xmm4
	movapd [esp + nb104_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb104_rsqH2M], xmm0
	movapd [esp + nb104_rsqH2H1], xmm4

	movapd xmm0, [esp + nb104_ixH2]
	movapd xmm1, [esp + nb104_iyH2]
	movapd xmm2, [esp + nb104_izH2]
	subsd  xmm0, [esp + nb104_jxH2]
	subsd  xmm1, [esp + nb104_jyH2]
	subsd  xmm2, [esp + nb104_jzH2]
	movapd [esp + nb104_dxH2H2], xmm0
	movapd [esp + nb104_dyH2H2], xmm1
	movapd [esp + nb104_dzH2H2], xmm2
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb104_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
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
	movapd  xmm3, [esp + nb104_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104_half] ;# iter1 
	mulsd   xmm7, [esp + nb104_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104_half] ;# rinv 
	mulsd   xmm5, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvH2H2], xmm1
	movapd [esp + nb104_rinvH2H1], xmm5

	movapd xmm0, [esp + nb104_rsqMM]
	movapd xmm4, [esp + nb104_rsqMH1]	
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
	movapd  xmm3, [esp + nb104_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb104_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104_half] ;# rinv 
	mulsd   xmm5, [esp + nb104_half] ;# rinv
	movapd [esp + nb104_rinvMM], xmm1
	movapd [esp + nb104_rinvMH1], xmm5

	movapd xmm0, [esp + nb104_rsqMH2]
	movapd xmm4, [esp + nb104_rsqH1M]	
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
	movapd  xmm3, [esp + nb104_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104_half] ;# iter1 
	mulsd   xmm7, [esp + nb104_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104_half] ;# rinv 
	mulsd   xmm5, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvMH2], xmm1
	movapd [esp + nb104_rinvH1M], xmm5

	movapd xmm0, [esp + nb104_rsqH1H1]
	movapd xmm4, [esp + nb104_rsqH1H2]	
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
	movapd  xmm3, [esp + nb104_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104_half] ;# iter1a 
	mulsd   xmm7, [esp + nb104_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104_half] ;# rinv 
	mulsd   xmm5, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvH1H1], xmm1
	movapd [esp + nb104_rinvH1H2], xmm5

	movapd xmm0, [esp + nb104_rsqH2M]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb104_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb104_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb104_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb104_half] ;# rinv 
	movapd [esp + nb104_rinvH2M], xmm1

	;# start with MM interaction 
	movapd xmm0, [esp + nb104_rinvMM]
	movapd xmm7, xmm0
	mulsd  xmm0, xmm0
	mulsd  xmm7, [esp + nb104_qqMM]
	mulsd  xmm0, xmm7	
	addsd  xmm7, [esp + nb104_vctot] 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb104_dxMM]
	mulsd xmm1, [esp + nb104_dyMM]
	mulsd xmm2, [esp + nb104_dzMM]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixM]
	addsd xmm1, [esp + nb104_fiyM]
	addsd xmm2, [esp + nb104_fizM]
	movlpd [esp + nb104_fjxM], xmm3
	movlpd [esp + nb104_fjyM], xmm4
	movlpd [esp + nb104_fjzM], xmm5
	movlpd [esp + nb104_fixM], xmm0
	movlpd [esp + nb104_fiyM], xmm1
	movlpd [esp + nb104_fizM], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb104_rinvMH1]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqMH]
	mulsd xmm0, xmm1	;# fsMH1  
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb104_dxMH1]
	mulsd xmm1, [esp + nb104_dyMH1]
	mulsd xmm2, [esp + nb104_dzMH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixM]
	addsd xmm1, [esp + nb104_fiyM]
	addsd xmm2, [esp + nb104_fizM]
	movlpd [esp + nb104_fjxH1], xmm3
	movlpd [esp + nb104_fjyH1], xmm4
	movlpd [esp + nb104_fjzH1], xmm5
	movlpd [esp + nb104_fixM], xmm0
	movlpd [esp + nb104_fiyM], xmm1
	movlpd [esp + nb104_fizM], xmm2

	;# M-H2 interaction  
	movapd xmm0, [esp + nb104_rinvMH2]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqMH]
	mulsd xmm0, xmm1	;# fsMH2  
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb104_dxMH2]
	mulsd xmm1, [esp + nb104_dyMH2]
	mulsd xmm2, [esp + nb104_dzMH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixM]
	addsd xmm1, [esp + nb104_fiyM]
	addsd xmm2, [esp + nb104_fizM]
	movlpd [esp + nb104_fjxH2], xmm3
	movlpd [esp + nb104_fjyH2], xmm4
	movlpd [esp + nb104_fjzH2], xmm5
	movlpd [esp + nb104_fixM], xmm0
	movlpd [esp + nb104_fiyM], xmm1
	movlpd [esp + nb104_fizM], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb104_rinvH1M]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqMH]
	mulsd xmm0, xmm1	;# fsH1M 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxM]
	movapd xmm4, [esp + nb104_fjyM]
	movapd xmm5, [esp + nb104_fjzM]
	mulsd xmm0, [esp + nb104_dxH1M]
	mulsd xmm1, [esp + nb104_dyH1M]
	mulsd xmm2, [esp + nb104_dzH1M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixH1]
	addsd xmm1, [esp + nb104_fiyH1]
	addsd xmm2, [esp + nb104_fizH1]
	movlpd [esp + nb104_fjxM], xmm3
	movlpd [esp + nb104_fjyM], xmm4
	movlpd [esp + nb104_fjzM], xmm5
	movlpd [esp + nb104_fixH1], xmm0
	movlpd [esp + nb104_fiyH1], xmm1
	movlpd [esp + nb104_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb104_rinvH1H1]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqHH]
	mulsd xmm0, xmm1	;# fsH1H1 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH1]
	movapd xmm4, [esp + nb104_fjyH1]
	movapd xmm5, [esp + nb104_fjzH1]
	mulsd xmm0, [esp + nb104_dxH1H1]
	mulsd xmm1, [esp + nb104_dyH1H1]
	mulsd xmm2, [esp + nb104_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixH1]
	addsd xmm1, [esp + nb104_fiyH1]
	addsd xmm2, [esp + nb104_fizH1]
	movlpd [esp + nb104_fjxH1], xmm3
	movlpd [esp + nb104_fjyH1], xmm4
	movlpd [esp + nb104_fjzH1], xmm5
	movlpd [esp + nb104_fixH1], xmm0
	movlpd [esp + nb104_fiyH1], xmm1
	movlpd [esp + nb104_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb104_rinvH1H2]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqHH]
	mulsd xmm0, xmm1	;# fsMH2  
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH2]
	movapd xmm4, [esp + nb104_fjyH2]
	movapd xmm5, [esp + nb104_fjzH2]
	mulsd xmm0, [esp + nb104_dxH1H2]
	mulsd xmm1, [esp + nb104_dyH1H2]
	mulsd xmm2, [esp + nb104_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixH1]
	addsd xmm1, [esp + nb104_fiyH1]
	addsd xmm2, [esp + nb104_fizH1]
	movlpd [esp + nb104_fjxH2], xmm3
	movlpd [esp + nb104_fjyH2], xmm4
	movlpd [esp + nb104_fjzH2], xmm5
	movlpd [esp + nb104_fixH1], xmm0
	movlpd [esp + nb104_fiyH1], xmm1
	movlpd [esp + nb104_fizH1], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb104_rinvH2M]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqMH]
	mulsd xmm0, xmm1	;# fsH2M 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxM]
	movapd xmm4, [esp + nb104_fjyM]
	movapd xmm5, [esp + nb104_fjzM]
	mulsd xmm0, [esp + nb104_dxH2M]
	mulsd xmm1, [esp + nb104_dyH2M]
	mulsd xmm2, [esp + nb104_dzH2M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixH2]
	addsd xmm1, [esp + nb104_fiyH2]
	addsd xmm2, [esp + nb104_fizH2]
	movlpd [esp + nb104_fjxM], xmm3
	movlpd [esp + nb104_fjyM], xmm4
	movlpd [esp + nb104_fjzM], xmm5
	movlpd [esp + nb104_fixH2], xmm0
	movlpd [esp + nb104_fiyH2], xmm1
	movlpd [esp + nb104_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb104_rinvH2H1]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqHH]
	mulsd xmm0, xmm1	;# fsH2H1 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH1]
	movapd xmm4, [esp + nb104_fjyH1]
	movapd xmm5, [esp + nb104_fjzH1]
	mulsd xmm0, [esp + nb104_dxH2H1]
	mulsd xmm1, [esp + nb104_dyH2H1]
	mulsd xmm2, [esp + nb104_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixH2]
	addsd xmm1, [esp + nb104_fiyH2]
	addsd xmm2, [esp + nb104_fizH2]
	movlpd [esp + nb104_fjxH1], xmm3
	movlpd [esp + nb104_fjyH1], xmm4
	movlpd [esp + nb104_fjzH1], xmm5
	movlpd [esp + nb104_fixH2], xmm0
	movlpd [esp + nb104_fiyH2], xmm1
	movlpd [esp + nb104_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb104_rinvH2H2]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb104_qqHH]
	mulsd xmm0, xmm1	;# fsH2H2 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movsd [esp + nb104_vctot], xmm7
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb104_fjxH2]
	movapd xmm4, [esp + nb104_fjyH2]
	movapd xmm5, [esp + nb104_fjzH2]
	mulsd xmm0, [esp + nb104_dxH2H2]
	mulsd xmm1, [esp + nb104_dyH2H2]
	mulsd xmm2, [esp + nb104_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb104_fixH2]
	addsd xmm1, [esp + nb104_fiyH2]
	addsd xmm2, [esp + nb104_fizH2]
	movlpd [esp + nb104_fjxH2], xmm3
	movlpd [esp + nb104_fjyH2], xmm4
	movlpd [esp + nb104_fjzH2], xmm5
	movlpd [esp + nb104_fixH2], xmm0
	movlpd [esp + nb104_fiyH2], xmm1
	movlpd [esp + nb104_fizH2], xmm2

	mov edi, [ebp + nb104_faction]
		
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8 + 24]
	movlpd xmm1, [edi + eax*8 + 32]
	movlpd xmm2, [edi + eax*8 + 40]
	movlpd xmm3, [edi + eax*8 + 48]
	movlpd xmm4, [edi + eax*8 + 56]
	movlpd xmm5, [edi + eax*8 + 64]
	movlpd xmm6, [edi + eax*8 + 72]
	movlpd xmm7, [edi + eax*8 + 80]
	addsd xmm0, [esp + nb104_fjxH1]
	addsd xmm1, [esp + nb104_fjyH1]
	addsd xmm2, [esp + nb104_fjzH1]
	addsd xmm3, [esp + nb104_fjxH2]
	addsd xmm4, [esp + nb104_fjyH2]
	addsd xmm5, [esp + nb104_fjzH2]
	addsd xmm6, [esp + nb104_fjxM]
	addsd xmm7, [esp + nb104_fjyM]
	movlpd [edi + eax*8 + 24], xmm0
	movlpd [edi + eax*8 + 32], xmm1
	movlpd [edi + eax*8 + 40], xmm2
	movlpd [edi + eax*8 + 48], xmm3
	movlpd [edi + eax*8 + 56], xmm4
	movlpd [edi + eax*8 + 64], xmm5
	movlpd [edi + eax*8 + 72], xmm6
	movlpd [edi + eax*8 + 80], xmm7

	movlpd xmm0, [edi + eax*8 + 88]
	addsd xmm0, [esp + nb104_fjzM]
	movlpd [edi + eax*8 + 88], xmm0
	
.nb104_updateouterdata:
	mov   ecx, [esp + nb104_ii3]
	mov   edi, [ebp + nb104_faction]
	mov   esi, [ebp + nb104_fshift]
	mov   edx, [esp + nb104_is3]

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb104_fixH1]
	movapd xmm1, [esp + nb104_fiyH1]
	movapd xmm2, [esp + nb104_fizH1]

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
	movsd  [edi + ecx*8 + 24],     xmm3
	movsd  [edi + ecx*8 + 32], xmm4
	movsd  [edi + ecx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb104_fixH2]
	movapd xmm1, [esp + nb104_fiyH2]
	movapd xmm2, [esp + nb104_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb104_fixM]
	movapd xmm1, [esp + nb104_fiyM]
	movapd xmm2, [esp + nb104_fizM]

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
	mov esi, [esp + nb104_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb104_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb104_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   eax, [ebp + nb104_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb104_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb104_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb104_n], esi
        jmp .nb104_outer
.nb104_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb104_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb104_end
        ;# non-zero, do one more workunit
        jmp   .nb104_threadloop
.nb104_end:
	emms

	mov eax, [esp + nb104_nouter]
	mov ebx, [esp + nb104_ninner]
	mov ecx, [ebp + nb104_outeriter]
	mov edx, [ebp + nb104_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb104_salign]
	add esp, eax
	add esp, 1432
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


	
.globl nb_kernel104nf_ia32_sse2
.globl _nb_kernel104nf_ia32_sse2
nb_kernel104nf_ia32_sse2:	
_nb_kernel104nf_ia32_sse2:	
.equiv          nb104nf_p_nri,          8
.equiv          nb104nf_iinr,           12
.equiv          nb104nf_jindex,         16
.equiv          nb104nf_jjnr,           20
.equiv          nb104nf_shift,          24
.equiv          nb104nf_shiftvec,       28
.equiv          nb104nf_fshift,         32
.equiv          nb104nf_gid,            36
.equiv          nb104nf_pos,            40
.equiv          nb104nf_faction,        44
.equiv          nb104nf_charge,         48
.equiv          nb104nf_p_facel,        52
.equiv          nb104nf_argkrf,         56
.equiv          nb104nf_argcrf,         60
.equiv          nb104nf_Vc,             64
.equiv          nb104nf_type,           68
.equiv          nb104nf_p_ntype,        72
.equiv          nb104nf_vdwparam,       76
.equiv          nb104nf_Vvdw,           80
.equiv          nb104nf_p_tabscale,     84
.equiv          nb104nf_VFtab,          88
.equiv          nb104nf_invsqrta,       92
.equiv          nb104nf_dvda,           96
.equiv          nb104nf_p_gbtabscale,   100
.equiv          nb104nf_GBtab,          104
.equiv          nb104nf_p_nthreads,     108
.equiv          nb104nf_count,          112
.equiv          nb104nf_mtx,            116
.equiv          nb104nf_outeriter,      120
.equiv          nb104nf_inneriter,      124
.equiv          nb104nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb104nf_ixM,            0
.equiv          nb104nf_iyM,            16
.equiv          nb104nf_izM,            32
.equiv          nb104nf_ixH1,           48
.equiv          nb104nf_iyH1,           64
.equiv          nb104nf_izH1,           80
.equiv          nb104nf_ixH2,           96
.equiv          nb104nf_iyH2,           112
.equiv          nb104nf_izH2,           128
.equiv          nb104nf_jxM,            144
.equiv          nb104nf_jyM,            160
.equiv          nb104nf_jzM,            176
.equiv          nb104nf_jxH1,           192
.equiv          nb104nf_jyH1,           208
.equiv          nb104nf_jzH1,           224
.equiv          nb104nf_jxH2,           240
.equiv          nb104nf_jyH2,           256
.equiv          nb104nf_jzH2,           272
.equiv          nb104nf_qqMM,           288
.equiv          nb104nf_qqMH,           304
.equiv          nb104nf_qqHH,           320
.equiv          nb104nf_vctot,          336
.equiv          nb104nf_half,           352
.equiv          nb104nf_three,          368
.equiv          nb104nf_rsqMM,          384
.equiv          nb104nf_rsqMH1,         400
.equiv          nb104nf_rsqMH2,         416
.equiv          nb104nf_rsqH1M,         432
.equiv          nb104nf_rsqH1H1,        448
.equiv          nb104nf_rsqH1H2,        464
.equiv          nb104nf_rsqH2M,         480
.equiv          nb104nf_rsqH2H1,        496
.equiv          nb104nf_rsqH2H2,        512
.equiv          nb104nf_rinvMM,         528
.equiv          nb104nf_rinvMH1,        544
.equiv          nb104nf_rinvMH2,        560
.equiv          nb104nf_rinvH1M,        576
.equiv          nb104nf_rinvH1H1,       592
.equiv          nb104nf_rinvH1H2,       608
.equiv          nb104nf_rinvH2M,        624
.equiv          nb104nf_rinvH2H1,       640
.equiv          nb104nf_rinvH2H2,       656
.equiv          nb104nf_is3,            672
.equiv          nb104nf_ii3,            676
.equiv          nb104nf_innerjjnr,      680
.equiv          nb104nf_innerk,         684
.equiv          nb104nf_n,              688
.equiv          nb104nf_nn1,            692
.equiv          nb104nf_nri,            696
.equiv          nb104nf_nouter,         700
.equiv          nb104nf_ninner,         704
.equiv          nb104nf_salign,         708
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 712		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb104nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb104nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb104nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb104nf_nouter], eax
	mov [esp + nb104nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb104nf_half], eax
	mov [esp + nb104nf_half+4], ebx
	movsd xmm1, [esp + nb104nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb104nf_half], xmm1
	movapd [esp + nb104nf_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb104nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb104nf_charge]
	movsd xmm3, [edx + ebx*8 + 24]	;# qM 
	movsd xmm4, xmm3		;# qM 
	movsd xmm5, [edx + ebx*8 + 8]	;# qH 
	mov esi, [ebp + nb104nf_p_facel]
	movsd xmm6, [esi]	;# facel 
	mulsd  xmm3, xmm3		;# qM*qM 
	mulsd  xmm4, xmm5		;# qM*qH 
	mulsd  xmm5, xmm5		;# qH*qH 
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb104nf_qqMM], xmm3
	movapd [esp + nb104nf_qqMH], xmm4
	movapd [esp + nb104nf_qqHH], xmm5

.nb104nf_threadloop:
        mov   esi, [ebp + nb104nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb104nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb104nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb104nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb104nf_n], eax
        mov [esp + nb104nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb104nf_outerstart
        jmp .nb104nf_end

.nb104nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb104nf_nouter]
	mov [esp + nb104nf_nouter], ebx

.nb104nf_outer:
	mov   eax, [ebp + nb104nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb104nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb104nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb104nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb104nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb104nf_ixH1], xmm3
	movapd [esp + nb104nf_iyH1], xmm4
	movapd [esp + nb104nf_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 48]
	addsd xmm1, [eax + ebx*8 + 56]
	addsd xmm2, [eax + ebx*8 + 64]		
	addsd xmm3, [eax + ebx*8 + 72]
	addsd xmm4, [eax + ebx*8 + 80]
	addsd xmm5, [eax + ebx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb104nf_ixH2], xmm0
	movapd [esp + nb104nf_iyH2], xmm1
	movapd [esp + nb104nf_izH2], xmm2
	movapd [esp + nb104nf_ixM], xmm3
	movapd [esp + nb104nf_iyM], xmm4
	movapd [esp + nb104nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb104nf_vctot], xmm4
	
	mov   eax, [ebp + nb104nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb104nf_pos]
	mov   eax, [ebp + nb104nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb104nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb104nf_ninner]
	mov   [esp + nb104nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb104nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb104nf_unroll_loop
	jmp   .nb104nf_checksingle
.nb104nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb104nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb104nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb104nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 24]
	movhpd xmm3, [esi + ebx*8 + 32]
	movhpd xmm4, [esi + ebx*8 + 40]
	movhpd xmm5, [esi + ebx*8 + 48]
	movhpd xmm6, [esi + ebx*8 + 56]
	movhpd xmm7, [esi + ebx*8 + 64]
	movapd 	[esp + nb104nf_jxH1], xmm2
	movapd 	[esp + nb104nf_jyH1], xmm3
	movapd 	[esp + nb104nf_jzH1], xmm4
	movapd 	[esp + nb104nf_jxH2], xmm5
	movapd 	[esp + nb104nf_jyH2], xmm6
	movapd 	[esp + nb104nf_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movhpd xmm2, [esi + ebx*8 + 72]
	movhpd xmm3, [esi + ebx*8 + 80]
	movhpd xmm4, [esi + ebx*8 + 88]
	movapd 	[esp + nb104nf_jxM], xmm2
	movapd 	[esp + nb104nf_jyM], xmm3
	movapd 	[esp + nb104nf_jzM], xmm4
	
	movapd xmm0, [esp + nb104nf_ixM]
	movapd xmm1, [esp + nb104nf_iyM]
	movapd xmm2, [esp + nb104nf_izM]
	movapd xmm3, [esp + nb104nf_ixM]
	movapd xmm4, [esp + nb104nf_iyM]
	movapd xmm5, [esp + nb104nf_izM]
	subpd  xmm0, [esp + nb104nf_jxM]
	subpd  xmm1, [esp + nb104nf_jyM]
	subpd  xmm2, [esp + nb104nf_jzM]
	subpd  xmm3, [esp + nb104nf_jxH1]
	subpd  xmm4, [esp + nb104nf_jyH1]
	subpd  xmm5, [esp + nb104nf_jzH1]
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
	movapd [esp + nb104nf_rsqMM], xmm0
	movapd [esp + nb104nf_rsqMH1], xmm3

	movapd xmm0, [esp + nb104nf_ixM]
	movapd xmm1, [esp + nb104nf_iyM]
	movapd xmm2, [esp + nb104nf_izM]
	movapd xmm3, [esp + nb104nf_ixH1]
	movapd xmm4, [esp + nb104nf_iyH1]
	movapd xmm5, [esp + nb104nf_izH1]
	subpd  xmm0, [esp + nb104nf_jxH2]
	subpd  xmm1, [esp + nb104nf_jyH2]
	subpd  xmm2, [esp + nb104nf_jzH2]
	subpd  xmm3, [esp + nb104nf_jxM]
	subpd  xmm4, [esp + nb104nf_jyM]
	subpd  xmm5, [esp + nb104nf_jzM]
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
	movapd [esp + nb104nf_rsqMH2], xmm0
	movapd [esp + nb104nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb104nf_ixH1]
	movapd xmm1, [esp + nb104nf_iyH1]
	movapd xmm2, [esp + nb104nf_izH1]
	movapd xmm3, [esp + nb104nf_ixH1]
	movapd xmm4, [esp + nb104nf_iyH1]
	movapd xmm5, [esp + nb104nf_izH1]
	subpd  xmm0, [esp + nb104nf_jxH1]
	subpd  xmm1, [esp + nb104nf_jyH1]
	subpd  xmm2, [esp + nb104nf_jzH1]
	subpd  xmm3, [esp + nb104nf_jxH2]
	subpd  xmm4, [esp + nb104nf_jyH2]
	subpd  xmm5, [esp + nb104nf_jzH2]
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
	movapd [esp + nb104nf_rsqH1H1], xmm0
	movapd [esp + nb104nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb104nf_ixH2]
	movapd xmm1, [esp + nb104nf_iyH2]
	movapd xmm2, [esp + nb104nf_izH2]
	movapd xmm3, [esp + nb104nf_ixH2]
	movapd xmm4, [esp + nb104nf_iyH2]
	movapd xmm5, [esp + nb104nf_izH2]
	subpd  xmm0, [esp + nb104nf_jxM]
	subpd  xmm1, [esp + nb104nf_jyM]
	subpd  xmm2, [esp + nb104nf_jzM]
	subpd  xmm3, [esp + nb104nf_jxH1]
	subpd  xmm4, [esp + nb104nf_jyH1]
	subpd  xmm5, [esp + nb104nf_jzH1]
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
	movapd [esp + nb104nf_rsqH2M], xmm0
	movapd [esp + nb104nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb104nf_ixH2]
	movapd xmm1, [esp + nb104nf_iyH2]
	movapd xmm2, [esp + nb104nf_izH2]
	subpd  xmm0, [esp + nb104nf_jxH2]
	subpd  xmm1, [esp + nb104nf_jyH2]
	subpd  xmm2, [esp + nb104nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb104nf_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0 (h2h2) , xmm4 (h2h1) 
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
	movapd  xmm3, [esp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvH2H2], xmm1
	movapd [esp + nb104nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb104nf_rsqMM]
	movapd xmm4, [esp + nb104nf_rsqMH1]	
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
	movapd  xmm3, [esp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb104nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb104nf_half] ;# rinv
	movapd [esp + nb104nf_rinvMM], xmm1
	movapd [esp + nb104nf_rinvMH1], xmm5

	movapd xmm0, [esp + nb104nf_rsqMH2]
	movapd xmm4, [esp + nb104nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvMH2], xmm1
	movapd [esp + nb104nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb104nf_rsqH1H1]
	movapd xmm4, [esp + nb104nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb104nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb104nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvH1H1], xmm1
	movapd [esp + nb104nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb104nf_rsqH2M]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvH2M], xmm1

	;# start with MM interaction 
	movapd xmm0, [esp + nb104nf_rinvMM]
	mulpd  xmm0, [esp + nb104nf_qqMM]	
	addpd  xmm0, [esp + nb104nf_vctot]
	
	;# other interactions 
	movapd xmm1, [esp + nb104nf_rinvMH1]
	movapd xmm2, [esp + nb104nf_rinvH1H1]
	
	addpd xmm1, [esp + nb104nf_rinvMH2]
	addpd xmm2, [esp + nb104nf_rinvH1H2]
	
	addpd xmm1, [esp + nb104nf_rinvH1M]
	addpd xmm2, [esp + nb104nf_rinvH2H1]

	addpd xmm1, [esp + nb104nf_rinvH2M]
	addpd xmm2, [esp + nb104nf_rinvH2H2]

	mulpd xmm1, [esp + nb104nf_qqMH]
	mulpd xmm2, [esp + nb104nf_qqHH]
	
	addpd xmm0, xmm1	
	addpd xmm0, xmm2

	movapd [esp + nb104nf_vctot], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb104nf_innerk],  2
	jl    .nb104nf_checksingle
	jmp   .nb104nf_unroll_loop
.nb104nf_checksingle:
	mov   edx, [esp + nb104nf_innerk]
	and   edx, 1
	jnz   .nb104nf_dosingle
	jmp   .nb104nf_updateouterdata
.nb104nf_dosingle:
	mov   edx, [esp + nb104nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	mov esi, [ebp + nb104nf_pos]
	lea   eax, [eax + eax*2]  

	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movapd 	[esp + nb104nf_jxH1], xmm2
	movapd 	[esp + nb104nf_jyH1], xmm3
	movapd 	[esp + nb104nf_jzH1], xmm4
	movapd 	[esp + nb104nf_jxH2], xmm5
	movapd 	[esp + nb104nf_jyH2], xmm6
	movapd 	[esp + nb104nf_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movapd 	[esp + nb104nf_jxM], xmm2
	movapd 	[esp + nb104nf_jyM], xmm3
	movapd 	[esp + nb104nf_jzM], xmm4
	
	movapd xmm0, [esp + nb104nf_ixM]
	movapd xmm1, [esp + nb104nf_iyM]
	movapd xmm2, [esp + nb104nf_izM]
	movapd xmm3, [esp + nb104nf_ixM]
	movapd xmm4, [esp + nb104nf_iyM]
	movapd xmm5, [esp + nb104nf_izM]
	subsd  xmm0, [esp + nb104nf_jxM]
	subsd  xmm1, [esp + nb104nf_jyM]
	subsd  xmm2, [esp + nb104nf_jzM]
	subsd  xmm3, [esp + nb104nf_jxH1]
	subsd  xmm4, [esp + nb104nf_jyH1]
	subsd  xmm5, [esp + nb104nf_jzH1]
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
	movapd [esp + nb104nf_rsqMM], xmm0
	movapd [esp + nb104nf_rsqMH1], xmm3

	movapd xmm0, [esp + nb104nf_ixM]
	movapd xmm1, [esp + nb104nf_iyM]
	movapd xmm2, [esp + nb104nf_izM]
	movapd xmm3, [esp + nb104nf_ixH1]
	movapd xmm4, [esp + nb104nf_iyH1]
	movapd xmm5, [esp + nb104nf_izH1]
	subsd  xmm0, [esp + nb104nf_jxH2]
	subsd  xmm1, [esp + nb104nf_jyH2]
	subsd  xmm2, [esp + nb104nf_jzH2]
	subsd  xmm3, [esp + nb104nf_jxM]
	subsd  xmm4, [esp + nb104nf_jyM]
	subsd  xmm5, [esp + nb104nf_jzM]
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
	movapd [esp + nb104nf_rsqMH2], xmm0
	movapd [esp + nb104nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb104nf_ixH1]
	movapd xmm1, [esp + nb104nf_iyH1]
	movapd xmm2, [esp + nb104nf_izH1]
	movapd xmm3, [esp + nb104nf_ixH1]
	movapd xmm4, [esp + nb104nf_iyH1]
	movapd xmm5, [esp + nb104nf_izH1]
	subsd  xmm0, [esp + nb104nf_jxH1]
	subsd  xmm1, [esp + nb104nf_jyH1]
	subsd  xmm2, [esp + nb104nf_jzH1]
	subsd  xmm3, [esp + nb104nf_jxH2]
	subsd  xmm4, [esp + nb104nf_jyH2]
	subsd  xmm5, [esp + nb104nf_jzH2]
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
	movapd [esp + nb104nf_rsqH1H1], xmm0
	movapd [esp + nb104nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb104nf_ixH2]
	movapd xmm1, [esp + nb104nf_iyH2]
	movapd xmm2, [esp + nb104nf_izH2]
	movapd xmm3, [esp + nb104nf_ixH2]
	movapd xmm4, [esp + nb104nf_iyH2]
	movapd xmm5, [esp + nb104nf_izH2]
	subsd  xmm0, [esp + nb104nf_jxM]
	subsd  xmm1, [esp + nb104nf_jyM]
	subsd  xmm2, [esp + nb104nf_jzM]
	subsd  xmm3, [esp + nb104nf_jxH1]
	subsd  xmm4, [esp + nb104nf_jyH1]
	subsd  xmm5, [esp + nb104nf_jzH1]
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
	movapd [esp + nb104nf_rsqH2M], xmm0
	movapd [esp + nb104nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb104nf_ixH2]
	movapd xmm1, [esp + nb104nf_iyH2]
	movapd xmm2, [esp + nb104nf_izH2]
	subsd  xmm0, [esp + nb104nf_jxH2]
	subsd  xmm1, [esp + nb104nf_jyH2]
	subsd  xmm2, [esp + nb104nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb104nf_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
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
	movapd  xmm3, [esp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvH2H2], xmm1
	movapd [esp + nb104nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb104nf_rsqMM]
	movapd xmm4, [esp + nb104nf_rsqMH1]	
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
	movapd  xmm3, [esp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb104nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb104nf_half] ;# rinv
	movapd [esp + nb104nf_rinvMM], xmm1
	movapd [esp + nb104nf_rinvMH1], xmm5

	movapd xmm0, [esp + nb104nf_rsqMH2]
	movapd xmm4, [esp + nb104nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvMH2], xmm1
	movapd [esp + nb104nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb104nf_rsqH1H1]
	movapd xmm4, [esp + nb104nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb104nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb104nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvH1H1], xmm1
	movapd [esp + nb104nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb104nf_rsqH2M]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb104nf_half] ;# rinv 
	movapd [esp + nb104nf_rinvH2M], xmm1

	;# start with MM interaction 
	movapd xmm0, [esp + nb104nf_rinvMM]
	mulpd  xmm0, [esp + nb104nf_qqMM]	
	addpd  xmm0, [esp + nb104nf_vctot]
	
	;# other interactions 
	movapd xmm1, [esp + nb104nf_rinvMH1]
	movapd xmm2, [esp + nb104nf_rinvH1H1]
	
	addsd xmm1, [esp + nb104nf_rinvMH2]
	addsd xmm2, [esp + nb104nf_rinvH1H2]
	
	addsd xmm1, [esp + nb104nf_rinvH1M]
	addsd xmm2, [esp + nb104nf_rinvH2H1]

	addsd xmm1, [esp + nb104nf_rinvH2M]
	addsd xmm2, [esp + nb104nf_rinvH2H2]

	mulsd xmm1, [esp + nb104nf_qqMH]
	mulsd xmm2, [esp + nb104nf_qqHH]
	
	addsd xmm0, xmm1	
	addsd xmm0, xmm2

	movlpd [esp + nb104nf_vctot], xmm0
	
.nb104nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb104nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb104nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movapd xmm7, [esp + nb104nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   eax, [ebp + nb104nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb104nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb104nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb104nf_n], esi
        jmp .nb104nf_outer
.nb104nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb104nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb104nf_end
        ;# non-zero, do one more workunit
        jmp   .nb104nf_threadloop
.nb104nf_end:
	emms

	mov eax, [esp + nb104nf_nouter]
	mov ebx, [esp + nb104nf_ninner]
	mov ecx, [ebp + nb104nf_outeriter]
	mov edx, [ebp + nb104nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb104nf_salign]
	add esp, eax
	add esp, 712
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


