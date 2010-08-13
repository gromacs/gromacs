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


	
.globl nb_kernel204_ia32_sse2
.globl _nb_kernel204_ia32_sse2
nb_kernel204_ia32_sse2:	
_nb_kernel204_ia32_sse2:	
.equiv          nb204_p_nri,            8
.equiv          nb204_iinr,             12
.equiv          nb204_jindex,           16
.equiv          nb204_jjnr,             20
.equiv          nb204_shift,            24
.equiv          nb204_shiftvec,         28
.equiv          nb204_fshift,           32
.equiv          nb204_gid,              36
.equiv          nb204_pos,              40
.equiv          nb204_faction,          44
.equiv          nb204_charge,           48
.equiv          nb204_p_facel,          52
.equiv          nb204_argkrf,           56
.equiv          nb204_argcrf,           60
.equiv          nb204_Vc,               64
.equiv          nb204_type,             68
.equiv          nb204_p_ntype,          72
.equiv          nb204_vdwparam,         76
.equiv          nb204_Vvdw,             80
.equiv          nb204_p_tabscale,       84
.equiv          nb204_VFtab,            88
.equiv          nb204_invsqrta,         92
.equiv          nb204_dvda,             96
.equiv          nb204_p_gbtabscale,     100
.equiv          nb204_GBtab,            104
.equiv          nb204_p_nthreads,       108
.equiv          nb204_count,            112
.equiv          nb204_mtx,              116
.equiv          nb204_outeriter,        120
.equiv          nb204_inneriter,        124
.equiv          nb204_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb204_ixM,              0
.equiv          nb204_iyM,              16
.equiv          nb204_izM,              32
.equiv          nb204_ixH1,             48
.equiv          nb204_iyH1,             64
.equiv          nb204_izH1,             80
.equiv          nb204_ixH2,             96
.equiv          nb204_iyH2,             112
.equiv          nb204_izH2,             128
.equiv          nb204_jxM,              144
.equiv          nb204_jyM,              160
.equiv          nb204_jzM,              176
.equiv          nb204_jxH1,             192
.equiv          nb204_jyH1,             208
.equiv          nb204_jzH1,             224
.equiv          nb204_jxH2,             240
.equiv          nb204_jyH2,             256
.equiv          nb204_jzH2,             272
.equiv          nb204_dxMM,             288
.equiv          nb204_dyMM,             304
.equiv          nb204_dzMM,             320
.equiv          nb204_dxMH1,            336
.equiv          nb204_dyMH1,            352
.equiv          nb204_dzMH1,            368
.equiv          nb204_dxMH2,            384
.equiv          nb204_dyMH2,            400
.equiv          nb204_dzMH2,            416
.equiv          nb204_dxH1M,            432
.equiv          nb204_dyH1M,            448
.equiv          nb204_dzH1M,            464
.equiv          nb204_dxH1H1,           480
.equiv          nb204_dyH1H1,           496
.equiv          nb204_dzH1H1,           512
.equiv          nb204_dxH1H2,           528
.equiv          nb204_dyH1H2,           544
.equiv          nb204_dzH1H2,           560
.equiv          nb204_dxH2M,            576
.equiv          nb204_dyH2M,            592
.equiv          nb204_dzH2M,            608
.equiv          nb204_dxH2H1,           624
.equiv          nb204_dyH2H1,           640
.equiv          nb204_dzH2H1,           656
.equiv          nb204_dxH2H2,           672
.equiv          nb204_dyH2H2,           688
.equiv          nb204_dzH2H2,           704
.equiv          nb204_qqMM,             720
.equiv          nb204_qqMH,             736
.equiv          nb204_qqHH,             752
.equiv          nb204_vctot,            768
.equiv          nb204_fixM,             784
.equiv          nb204_fiyM,             800
.equiv          nb204_fizM,             816
.equiv          nb204_fixH1,            832
.equiv          nb204_fiyH1,            848
.equiv          nb204_fizH1,            864
.equiv          nb204_fixH2,            880
.equiv          nb204_fiyH2,            896
.equiv          nb204_fizH2,            912
.equiv          nb204_fjxM,             928
.equiv          nb204_fjyM,             944
.equiv          nb204_fjzM,             960
.equiv          nb204_fjxH1,            976
.equiv          nb204_fjyH1,            992
.equiv          nb204_fjzH1,            1008
.equiv          nb204_fjxH2,            1024
.equiv          nb204_fjyH2,            1040
.equiv          nb204_fjzH2,            1056
.equiv          nb204_half,             1072
.equiv          nb204_three,            1088
.equiv          nb204_rsqMM,            1104
.equiv          nb204_rsqMH1,           1120
.equiv          nb204_rsqMH2,           1136
.equiv          nb204_rsqH1M,           1152
.equiv          nb204_rsqH1H1,          1168
.equiv          nb204_rsqH1H2,          1184
.equiv          nb204_rsqH2M,           1200
.equiv          nb204_rsqH2H1,          1216
.equiv          nb204_rsqH2H2,          1232
.equiv          nb204_rinvMM,           1248
.equiv          nb204_rinvMH1,          1264
.equiv          nb204_rinvMH2,          1280
.equiv          nb204_rinvH1M,          1296
.equiv          nb204_rinvH1H1,         1312
.equiv          nb204_rinvH1H2,         1328
.equiv          nb204_rinvH2M,          1344
.equiv          nb204_rinvH2H1,         1360
.equiv          nb204_rinvH2H2,         1376
.equiv          nb204_two,              1392
.equiv          nb204_krf,              1408
.equiv          nb204_crf,              1424
.equiv          nb204_is3,              1440
.equiv          nb204_ii3,              1444
.equiv          nb204_innerjjnr,        1448
.equiv          nb204_innerk,           1452
.equiv          nb204_n,                1456
.equiv          nb204_nn1,              1460
.equiv          nb204_nri,              1464
.equiv          nb204_nouter,           1468
.equiv          nb204_ninner,           1472
.equiv          nb204_salign,           1476
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 1480		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb204_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb204_p_nri]
	mov ecx, [ecx]
	mov [esp + nb204_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb204_nouter], eax
	mov [esp + nb204_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb204_half], eax
	mov [esp + nb204_half+4], ebx
	movsd xmm1, [esp + nb204_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb204_half], xmm1
	movapd [esp + nb204_two], xmm2
	movapd [esp + nb204_three], xmm3

	mov esi, [ebp + nb204_argkrf]
	mov edi, [ebp + nb204_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb204_krf], xmm5
	movapd [esp + nb204_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb204_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb204_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb204_p_facel]
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
	movapd [esp + nb204_qqMM], xmm3
	movapd [esp + nb204_qqMH], xmm4
	movapd [esp + nb204_qqHH], xmm5
	
.nb204_threadloop:
        mov   esi, [ebp + nb204_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb204_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb204_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb204_n], eax
        mov [esp + nb204_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204_outerstart
        jmp .nb204_end

.nb204_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb204_nouter]
	mov [esp + nb204_nouter], ebx

.nb204_outer:
	mov   eax, [ebp + nb204_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb204_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb204_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb204_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb204_pos]    ;# eax = base of pos[]  
	mov   [esp + nb204_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb204_ixH1], xmm3
	movapd [esp + nb204_iyH1], xmm4
	movapd [esp + nb204_izH1], xmm5

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
	movapd [esp + nb204_ixH2], xmm0
	movapd [esp + nb204_iyH2], xmm1
	movapd [esp + nb204_izH2], xmm2
	movapd [esp + nb204_ixM], xmm3
	movapd [esp + nb204_iyM], xmm4
	movapd [esp + nb204_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb204_vctot], xmm4
	movapd [esp + nb204_fixM], xmm4
	movapd [esp + nb204_fiyM], xmm4
	movapd [esp + nb204_fizM], xmm4
	movapd [esp + nb204_fixH1], xmm4
	movapd [esp + nb204_fiyH1], xmm4
	movapd [esp + nb204_fizH1], xmm4
	movapd [esp + nb204_fixH2], xmm4
	movapd [esp + nb204_fiyH2], xmm4
	movapd [esp + nb204_fizH2], xmm4
	
	mov   eax, [ebp + nb204_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb204_pos]
	mov   edi, [ebp + nb204_faction]	
	mov   eax, [ebp + nb204_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb204_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb204_ninner]
	mov   [esp + nb204_ninner], ecx
	add   edx, 0
	mov   [esp + nb204_innerk], edx    ;# number of innerloop atoms 
	jge   .nb204_unroll_loop
	jmp   .nb204_checksingle
.nb204_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb204_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb204_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb204_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb204_jxH1], xmm2
	movapd 	[esp + nb204_jyH1], xmm3
	movapd 	[esp + nb204_jzH1], xmm4
	movapd 	[esp + nb204_jxH2], xmm5
	movapd 	[esp + nb204_jyH2], xmm6
	movapd 	[esp + nb204_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movhpd xmm2, [esi + ebx*8 + 72]
	movhpd xmm3, [esi + ebx*8 + 80]
	movhpd xmm4, [esi + ebx*8 + 88]
	movapd 	[esp + nb204_jxM], xmm2
	movapd 	[esp + nb204_jyM], xmm3
	movapd 	[esp + nb204_jzM], xmm4
	
	movapd xmm0, [esp + nb204_ixM]
	movapd xmm1, [esp + nb204_iyM]
	movapd xmm2, [esp + nb204_izM]
	movapd xmm3, [esp + nb204_ixM]
	movapd xmm4, [esp + nb204_iyM]
	movapd xmm5, [esp + nb204_izM]
	subpd  xmm0, [esp + nb204_jxM]
	subpd  xmm1, [esp + nb204_jyM]
	subpd  xmm2, [esp + nb204_jzM]
	subpd  xmm3, [esp + nb204_jxH1]
	subpd  xmm4, [esp + nb204_jyH1]
	subpd  xmm5, [esp + nb204_jzH1]
	movapd [esp + nb204_dxMM], xmm0
	movapd [esp + nb204_dyMM], xmm1
	movapd [esp + nb204_dzMM], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb204_dxMH1], xmm3
	movapd [esp + nb204_dyMH1], xmm4
	movapd [esp + nb204_dzMH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb204_rsqMM], xmm0
	movapd [esp + nb204_rsqMH1], xmm3

	movapd xmm0, [esp + nb204_ixM]
	movapd xmm1, [esp + nb204_iyM]
	movapd xmm2, [esp + nb204_izM]
	movapd xmm3, [esp + nb204_ixH1]
	movapd xmm4, [esp + nb204_iyH1]
	movapd xmm5, [esp + nb204_izH1]
	subpd  xmm0, [esp + nb204_jxH2]
	subpd  xmm1, [esp + nb204_jyH2]
	subpd  xmm2, [esp + nb204_jzH2]
	subpd  xmm3, [esp + nb204_jxM]
	subpd  xmm4, [esp + nb204_jyM]
	subpd  xmm5, [esp + nb204_jzM]
	movapd [esp + nb204_dxMH2], xmm0
	movapd [esp + nb204_dyMH2], xmm1
	movapd [esp + nb204_dzMH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb204_dxH1M], xmm3
	movapd [esp + nb204_dyH1M], xmm4
	movapd [esp + nb204_dzH1M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb204_rsqMH2], xmm0
	movapd [esp + nb204_rsqH1M], xmm3

	movapd xmm0, [esp + nb204_ixH1]
	movapd xmm1, [esp + nb204_iyH1]
	movapd xmm2, [esp + nb204_izH1]
	movapd xmm3, [esp + nb204_ixH1]
	movapd xmm4, [esp + nb204_iyH1]
	movapd xmm5, [esp + nb204_izH1]
	subpd  xmm0, [esp + nb204_jxH1]
	subpd  xmm1, [esp + nb204_jyH1]
	subpd  xmm2, [esp + nb204_jzH1]
	subpd  xmm3, [esp + nb204_jxH2]
	subpd  xmm4, [esp + nb204_jyH2]
	subpd  xmm5, [esp + nb204_jzH2]
	movapd [esp + nb204_dxH1H1], xmm0
	movapd [esp + nb204_dyH1H1], xmm1
	movapd [esp + nb204_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb204_dxH1H2], xmm3
	movapd [esp + nb204_dyH1H2], xmm4
	movapd [esp + nb204_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb204_rsqH1H1], xmm0
	movapd [esp + nb204_rsqH1H2], xmm3

	movapd xmm0, [esp + nb204_ixH2]
	movapd xmm1, [esp + nb204_iyH2]
	movapd xmm2, [esp + nb204_izH2]
	movapd xmm3, [esp + nb204_ixH2]
	movapd xmm4, [esp + nb204_iyH2]
	movapd xmm5, [esp + nb204_izH2]
	subpd  xmm0, [esp + nb204_jxM]
	subpd  xmm1, [esp + nb204_jyM]
	subpd  xmm2, [esp + nb204_jzM]
	subpd  xmm3, [esp + nb204_jxH1]
	subpd  xmm4, [esp + nb204_jyH1]
	subpd  xmm5, [esp + nb204_jzH1]
	movapd [esp + nb204_dxH2M], xmm0
	movapd [esp + nb204_dyH2M], xmm1
	movapd [esp + nb204_dzH2M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb204_dxH2H1], xmm3
	movapd [esp + nb204_dyH2H1], xmm4
	movapd [esp + nb204_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb204_rsqH2M], xmm0
	movapd [esp + nb204_rsqH2H1], xmm4

	movapd xmm0, [esp + nb204_ixH2]
	movapd xmm1, [esp + nb204_iyH2]
	movapd xmm2, [esp + nb204_izH2]
	subpd  xmm0, [esp + nb204_jxH2]
	subpd  xmm1, [esp + nb204_jyH2]
	subpd  xmm2, [esp + nb204_jzH2]
	movapd [esp + nb204_dxH2H2], xmm0
	movapd [esp + nb204_dyH2H2], xmm1
	movapd [esp + nb204_dzH2H2], xmm2
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb204_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
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
	movapd  xmm3, [esp + nb204_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204_half] ;# iter1 
	mulpd   xmm7, [esp + nb204_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204_half] ;# rinv 
	mulpd   xmm5, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvH2H2], xmm1
	movapd [esp + nb204_rinvH2H1], xmm5

	movapd xmm0, [esp + nb204_rsqMM]
	movapd xmm4, [esp + nb204_rsqMH1]	
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
	movapd  xmm3, [esp + nb204_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb204_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204_half] ;# rinv 
	mulpd   xmm5, [esp + nb204_half] ;# rinv
	movapd [esp + nb204_rinvMM], xmm1
	movapd [esp + nb204_rinvMH1], xmm5

	movapd xmm0, [esp + nb204_rsqMH2]
	movapd xmm4, [esp + nb204_rsqH1M]	
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
	movapd  xmm3, [esp + nb204_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204_half] ;# iter1 
	mulpd   xmm7, [esp + nb204_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204_half] ;# rinv 
	mulpd   xmm5, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvMH2], xmm1
	movapd [esp + nb204_rinvH1M], xmm5

	movapd xmm0, [esp + nb204_rsqH1H1]
	movapd xmm4, [esp + nb204_rsqH1H2]	
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
	movapd  xmm3, [esp + nb204_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204_half] ;# iter1a 
	mulpd   xmm7, [esp + nb204_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204_half] ;# rinv 
	mulpd   xmm5, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvH1H1], xmm1
	movapd [esp + nb204_rinvH1H2], xmm5

	movapd xmm0, [esp + nb204_rsqH2M]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb204_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb204_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb204_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvH2M], xmm1
	
	;# start with MM interaction 
	movapd xmm0, [esp + nb204_rinvMM]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]
	mulpd  xmm0, xmm0	;# rinvsq 
	mulpd  xmm5, [esp + nb204_rsqMM] ;# xmm5=krsq 
	movapd xmm6, xmm5
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subpd  xmm6, [esp + nb204_crf]
	
	mulpd  xmm6, [esp + nb204_qqMM] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqMM] ;# xmm7 = coul part of fscal 
	
	addpd  xmm6, [esp + nb204_vctot] ;# local vctot summation variable 
	mulpd  xmm0, xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb204_dxMM]
	mulpd xmm1, [esp + nb204_dyMM]
	mulpd xmm2, [esp + nb204_dzMM]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixM]
	addpd xmm1, [esp + nb204_fiyM]
	addpd xmm2, [esp + nb204_fizM]
	movapd [esp + nb204_fjxM], xmm3
	movapd [esp + nb204_fjyM], xmm4
	movapd [esp + nb204_fjzM], xmm5
	movapd [esp + nb204_fixM], xmm0
	movapd [esp + nb204_fiyM], xmm1
	movapd [esp + nb204_fizM], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb204_rinvMH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqMH1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulpd  xmm0, xmm0
	subpd  xmm4, [esp + nb204_crf]
	mulpd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH1  
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb204_dxMH1]
	mulpd xmm1, [esp + nb204_dyMH1]
	mulpd xmm2, [esp + nb204_dzMH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixM]
	addpd xmm1, [esp + nb204_fiyM]
	addpd xmm2, [esp + nb204_fizM]
	movapd [esp + nb204_fjxH1], xmm3
	movapd [esp + nb204_fjyH1], xmm4
	movapd [esp + nb204_fjzH1], xmm5
	movapd [esp + nb204_fixM], xmm0
	movapd [esp + nb204_fiyM], xmm1
	movapd [esp + nb204_fizM], xmm2

	;# M-H2 interaction  
	movapd xmm0, [esp + nb204_rinvMH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqMH2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb204_crf]
	mulpd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb204_dxMH2]
	mulpd xmm1, [esp + nb204_dyMH2]
	mulpd xmm2, [esp + nb204_dzMH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixM]
	addpd xmm1, [esp + nb204_fiyM]
	addpd xmm2, [esp + nb204_fizM]
	movapd [esp + nb204_fjxH2], xmm3
	movapd [esp + nb204_fjyH2], xmm4
	movapd [esp + nb204_fjzH2], xmm5
	movapd [esp + nb204_fixM], xmm0
	movapd [esp + nb204_fiyM], xmm1
	movapd [esp + nb204_fizM], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb204_rinvH1M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqH1M] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb204_crf]
	mulpd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxM]
	movapd xmm4, [esp + nb204_fjyM]
	movapd xmm5, [esp + nb204_fjzM]
	mulpd xmm0, [esp + nb204_dxH1M]
	mulpd xmm1, [esp + nb204_dyH1M]
	mulpd xmm2, [esp + nb204_dzH1M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixH1]
	addpd xmm1, [esp + nb204_fiyH1]
	addpd xmm2, [esp + nb204_fizH1]
	movapd [esp + nb204_fjxM], xmm3
	movapd [esp + nb204_fjyM], xmm4
	movapd [esp + nb204_fjzM], xmm5
	movapd [esp + nb204_fixH1], xmm0
	movapd [esp + nb204_fiyH1], xmm1
	movapd [esp + nb204_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb204_rinvH1H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqH1H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb204_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxH1]
	movapd xmm4, [esp + nb204_fjyH1]
	movapd xmm5, [esp + nb204_fjzH1]
	mulpd xmm0, [esp + nb204_dxH1H1]
	mulpd xmm1, [esp + nb204_dyH1H1]
	mulpd xmm2, [esp + nb204_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixH1]
	addpd xmm1, [esp + nb204_fiyH1]
	addpd xmm2, [esp + nb204_fizH1]
	movapd [esp + nb204_fjxH1], xmm3
	movapd [esp + nb204_fjyH1], xmm4
	movapd [esp + nb204_fjzH1], xmm5
	movapd [esp + nb204_fixH1], xmm0
	movapd [esp + nb204_fiyH1], xmm1
	movapd [esp + nb204_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb204_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqH1H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb204_crf]
	mulpd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb204_fjxH2]
	movapd xmm4, [esp + nb204_fjyH2]
	movapd xmm5, [esp + nb204_fjzH2]
	mulpd xmm0, [esp + nb204_dxH1H2]
	mulpd xmm1, [esp + nb204_dyH1H2]
	mulpd xmm2, [esp + nb204_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixH1]
	addpd xmm1, [esp + nb204_fiyH1]
	addpd xmm2, [esp + nb204_fizH1]
	movapd [esp + nb204_fjxH2], xmm3
	movapd [esp + nb204_fjyH2], xmm4
	movapd [esp + nb204_fjzH2], xmm5
	movapd [esp + nb204_fixH1], xmm0
	movapd [esp + nb204_fiyH1], xmm1
	movapd [esp + nb204_fizH1], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb204_rinvH2M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqH2M] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb204_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxM]
	movapd xmm4, [esp + nb204_fjyM]
	movapd xmm5, [esp + nb204_fjzM]
	mulpd xmm0, [esp + nb204_dxH2M]
	mulpd xmm1, [esp + nb204_dyH2M]
	mulpd xmm2, [esp + nb204_dzH2M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixH2]
	addpd xmm1, [esp + nb204_fiyH2]
	addpd xmm2, [esp + nb204_fizH2]
	movapd [esp + nb204_fjxM], xmm3
	movapd [esp + nb204_fjyM], xmm4
	movapd [esp + nb204_fjzM], xmm5
	movapd [esp + nb204_fixH2], xmm0
	movapd [esp + nb204_fiyH2], xmm1
	movapd [esp + nb204_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb204_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqH2H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb204_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxH1]
	movapd xmm4, [esp + nb204_fjyH1]
	movapd xmm5, [esp + nb204_fjzH1]
	mulpd xmm0, [esp + nb204_dxH2H1]
	mulpd xmm1, [esp + nb204_dyH2H1]
	mulpd xmm2, [esp + nb204_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixH2]
	addpd xmm1, [esp + nb204_fiyH2]
	addpd xmm2, [esp + nb204_fizH2]
	movapd [esp + nb204_fjxH1], xmm3
	movapd [esp + nb204_fjyH1], xmm4
	movapd [esp + nb204_fjzH1], xmm5
	movapd [esp + nb204_fixH2], xmm0
	movapd [esp + nb204_fiyH2], xmm1
	movapd [esp + nb204_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb204_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb204_rsqH2H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb204_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb204_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm1, xmm0
	movapd [esp + nb204_vctot], xmm6
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb204_fjxH2]
	movapd xmm4, [esp + nb204_fjyH2]
	movapd xmm5, [esp + nb204_fjzH2]
	mulpd xmm0, [esp + nb204_dxH2H2]
	mulpd xmm1, [esp + nb204_dyH2H2]
	mulpd xmm2, [esp + nb204_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb204_fixH2]
	addpd xmm1, [esp + nb204_fiyH2]
	addpd xmm2, [esp + nb204_fizH2]
	movapd [esp + nb204_fjxH2], xmm3
	movapd [esp + nb204_fjyH2], xmm4
	movapd [esp + nb204_fjzH2], xmm5
	movapd [esp + nb204_fixH2], xmm0
	movapd [esp + nb204_fiyH2], xmm1
	movapd [esp + nb204_fizH2], xmm2

	mov edi, [ebp + nb204_faction]
		
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
	addpd xmm0, [esp + nb204_fjxH1]
	addpd xmm1, [esp + nb204_fjyH1]
	addpd xmm2, [esp + nb204_fjzH1]
	addpd xmm3, [esp + nb204_fjxH2]
	addpd xmm4, [esp + nb204_fjyH2]
	addpd xmm5, [esp + nb204_fjzH2]
	addpd xmm6, [esp + nb204_fjxM]
	addpd xmm7, [esp + nb204_fjyM]
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
	addpd xmm0, [esp + nb204_fjzM]
	movlpd [edi + eax*8 + 88], xmm0
	movhpd [edi + ebx*8 + 88], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb204_innerk],  2
	jl    .nb204_checksingle
	jmp   .nb204_unroll_loop
.nb204_checksingle:
	mov   edx, [esp + nb204_innerk]
	and   edx, 1
	jnz   .nb204_dosingle
	jmp   .nb204_updateouterdata
.nb204_dosingle:
	mov   edx, [esp + nb204_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]
	
	mov esi, [ebp + nb204_pos]
	lea   eax, [eax + eax*2]  

	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movapd 	[esp + nb204_jxH1], xmm2
	movapd 	[esp + nb204_jyH1], xmm3
	movapd 	[esp + nb204_jzH1], xmm4
	movapd 	[esp + nb204_jxH2], xmm5
	movapd 	[esp + nb204_jyH2], xmm6
	movapd 	[esp + nb204_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movapd 	[esp + nb204_jxM], xmm2
	movapd 	[esp + nb204_jyM], xmm3
	movapd 	[esp + nb204_jzM], xmm4
	
	movapd xmm0, [esp + nb204_ixM]
	movapd xmm1, [esp + nb204_iyM]
	movapd xmm2, [esp + nb204_izM]
	movapd xmm3, [esp + nb204_ixM]
	movapd xmm4, [esp + nb204_iyM]
	movapd xmm5, [esp + nb204_izM]
	subsd  xmm0, [esp + nb204_jxM]
	subsd  xmm1, [esp + nb204_jyM]
	subsd  xmm2, [esp + nb204_jzM]
	subsd  xmm3, [esp + nb204_jxH1]
	subsd  xmm4, [esp + nb204_jyH1]
	subsd  xmm5, [esp + nb204_jzH1]
	movapd [esp + nb204_dxMM], xmm0
	movapd [esp + nb204_dyMM], xmm1
	movapd [esp + nb204_dzMM], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb204_dxMH1], xmm3
	movapd [esp + nb204_dyMH1], xmm4
	movapd [esp + nb204_dzMH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb204_rsqMM], xmm0
	movapd [esp + nb204_rsqMH1], xmm3

	movapd xmm0, [esp + nb204_ixM]
	movapd xmm1, [esp + nb204_iyM]
	movapd xmm2, [esp + nb204_izM]
	movapd xmm3, [esp + nb204_ixH1]
	movapd xmm4, [esp + nb204_iyH1]
	movapd xmm5, [esp + nb204_izH1]
	subsd  xmm0, [esp + nb204_jxH2]
	subsd  xmm1, [esp + nb204_jyH2]
	subsd  xmm2, [esp + nb204_jzH2]
	subsd  xmm3, [esp + nb204_jxM]
	subsd  xmm4, [esp + nb204_jyM]
	subsd  xmm5, [esp + nb204_jzM]
	movapd [esp + nb204_dxMH2], xmm0
	movapd [esp + nb204_dyMH2], xmm1
	movapd [esp + nb204_dzMH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb204_dxH1M], xmm3
	movapd [esp + nb204_dyH1M], xmm4
	movapd [esp + nb204_dzH1M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb204_rsqMH2], xmm0
	movapd [esp + nb204_rsqH1M], xmm3

	movapd xmm0, [esp + nb204_ixH1]
	movapd xmm1, [esp + nb204_iyH1]
	movapd xmm2, [esp + nb204_izH1]
	movapd xmm3, [esp + nb204_ixH1]
	movapd xmm4, [esp + nb204_iyH1]
	movapd xmm5, [esp + nb204_izH1]
	subsd  xmm0, [esp + nb204_jxH1]
	subsd  xmm1, [esp + nb204_jyH1]
	subsd  xmm2, [esp + nb204_jzH1]
	subsd  xmm3, [esp + nb204_jxH2]
	subsd  xmm4, [esp + nb204_jyH2]
	subsd  xmm5, [esp + nb204_jzH2]
	movapd [esp + nb204_dxH1H1], xmm0
	movapd [esp + nb204_dyH1H1], xmm1
	movapd [esp + nb204_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb204_dxH1H2], xmm3
	movapd [esp + nb204_dyH1H2], xmm4
	movapd [esp + nb204_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb204_rsqH1H1], xmm0
	movapd [esp + nb204_rsqH1H2], xmm3

	movapd xmm0, [esp + nb204_ixH2]
	movapd xmm1, [esp + nb204_iyH2]
	movapd xmm2, [esp + nb204_izH2]
	movapd xmm3, [esp + nb204_ixH2]
	movapd xmm4, [esp + nb204_iyH2]
	movapd xmm5, [esp + nb204_izH2]
	subsd  xmm0, [esp + nb204_jxM]
	subsd  xmm1, [esp + nb204_jyM]
	subsd  xmm2, [esp + nb204_jzM]
	subsd  xmm3, [esp + nb204_jxH1]
	subsd  xmm4, [esp + nb204_jyH1]
	subsd  xmm5, [esp + nb204_jzH1]
	movapd [esp + nb204_dxH2M], xmm0
	movapd [esp + nb204_dyH2M], xmm1
	movapd [esp + nb204_dzH2M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb204_dxH2H1], xmm3
	movapd [esp + nb204_dyH2H1], xmm4
	movapd [esp + nb204_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb204_rsqH2M], xmm0
	movapd [esp + nb204_rsqH2H1], xmm4

	movapd xmm0, [esp + nb204_ixH2]
	movapd xmm1, [esp + nb204_iyH2]
	movapd xmm2, [esp + nb204_izH2]
	subsd  xmm0, [esp + nb204_jxH2]
	subsd  xmm1, [esp + nb204_jyH2]
	subsd  xmm2, [esp + nb204_jzH2]
	movapd [esp + nb204_dxH2H2], xmm0
	movapd [esp + nb204_dyH2H2], xmm1
	movapd [esp + nb204_dzH2H2], xmm2
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb204_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb204_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204_half] ;# iter1 
	mulsd   xmm7, [esp + nb204_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204_half] ;# rinv 
	mulsd   xmm5, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvH2H2], xmm1
	movapd [esp + nb204_rinvH2H1], xmm5

	movapd xmm0, [esp + nb204_rsqMM]
	movapd xmm4, [esp + nb204_rsqMH1]	
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
	movapd  xmm3, [esp + nb204_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb204_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204_half] ;# rinv 
	mulsd   xmm5, [esp + nb204_half] ;# rinv
	movapd [esp + nb204_rinvMM], xmm1
	movapd [esp + nb204_rinvMH1], xmm5

	movapd xmm0, [esp + nb204_rsqMH2]
	movapd xmm4, [esp + nb204_rsqH1M]	
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
	movapd  xmm3, [esp + nb204_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204_half] ;# iter1 
	mulsd   xmm7, [esp + nb204_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204_half] ;# rinv 
	mulsd   xmm5, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvMH2], xmm1
	movapd [esp + nb204_rinvH1M], xmm5

	movapd xmm0, [esp + nb204_rsqH1H1]
	movapd xmm4, [esp + nb204_rsqH1H2]	
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
	movapd  xmm3, [esp + nb204_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204_half] ;# iter1a 
	mulsd   xmm7, [esp + nb204_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204_half] ;# rinv 
	mulsd   xmm5, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvH1H1], xmm1
	movapd [esp + nb204_rinvH1H2], xmm5

	movapd xmm0, [esp + nb204_rsqH2M]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb204_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb204_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb204_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb204_half] ;# rinv 
	movapd [esp + nb204_rinvH2M], xmm1
	
	;# start with MM interaction 
	movapd xmm0, [esp + nb204_rinvMM]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]
	mulsd  xmm0, xmm0
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	mulsd  xmm5, [esp + nb204_rsqMM] ;# xmm5=krsq 
	movapd xmm6, xmm5
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subsd  xmm6, [esp + nb204_crf]
	
	mulsd  xmm6, [esp + nb204_qqMM] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqMM] ;# xmm7 = coul part of fscal 
	
	addsd  xmm6, [esp + nb204_vctot] ;# local vctot summation variable 
	mulsd  xmm0, xmm7
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb204_dxMM]
	mulsd xmm1, [esp + nb204_dyMM]
	mulsd xmm2, [esp + nb204_dzMM]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixM]
	addsd xmm1, [esp + nb204_fiyM]
	addsd xmm2, [esp + nb204_fizM]
	movlpd [esp + nb204_fjxM], xmm3
	movlpd [esp + nb204_fjyM], xmm4
	movlpd [esp + nb204_fjzM], xmm5
	movlpd [esp + nb204_fixM], xmm0
	movlpd [esp + nb204_fiyM], xmm1
	movlpd [esp + nb204_fizM], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb204_rinvMH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqMH1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulsd  xmm0, xmm0
	subsd  xmm4, [esp + nb204_crf]
	mulsd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH1  
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb204_dxMH1]
	mulsd xmm1, [esp + nb204_dyMH1]
	mulsd xmm2, [esp + nb204_dzMH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixM]
	addsd xmm1, [esp + nb204_fiyM]
	addsd xmm2, [esp + nb204_fizM]
	movlpd [esp + nb204_fjxH1], xmm3
	movlpd [esp + nb204_fjyH1], xmm4
	movlpd [esp + nb204_fjzH1], xmm5
	movlpd [esp + nb204_fixM], xmm0
	movlpd [esp + nb204_fiyM], xmm1
	movlpd [esp + nb204_fizM], xmm2

	;# M-H2 interaction  
	movapd xmm0, [esp + nb204_rinvMH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqMH2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulsd  xmm0, xmm0
	subsd  xmm4, [esp + nb204_crf]
	mulsd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb204_dxMH2]
	mulsd xmm1, [esp + nb204_dyMH2]
	mulsd xmm2, [esp + nb204_dzMH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixM]
	addsd xmm1, [esp + nb204_fiyM]
	addsd xmm2, [esp + nb204_fizM]
	movlpd [esp + nb204_fjxH2], xmm3
	movlpd [esp + nb204_fjyH2], xmm4
	movlpd [esp + nb204_fjzH2], xmm5
	movlpd [esp + nb204_fixM], xmm0
	movlpd [esp + nb204_fiyM], xmm1
	movlpd [esp + nb204_fizM], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb204_rinvH1M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqH1M] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulsd xmm0, xmm0
	subsd  xmm4, [esp + nb204_crf]
	mulsd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxM]
	movapd xmm4, [esp + nb204_fjyM]
	movapd xmm5, [esp + nb204_fjzM]
	mulsd xmm0, [esp + nb204_dxH1M]
	mulsd xmm1, [esp + nb204_dyH1M]
	mulsd xmm2, [esp + nb204_dzH1M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixH1]
	addsd xmm1, [esp + nb204_fiyH1]
	addsd xmm2, [esp + nb204_fizH1]
	movlpd [esp + nb204_fjxM], xmm3
	movlpd [esp + nb204_fjyM], xmm4
	movlpd [esp + nb204_fjzM], xmm5
	movlpd [esp + nb204_fixH1], xmm0
	movlpd [esp + nb204_fiyH1], xmm1
	movlpd [esp + nb204_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb204_rinvH1H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqH1H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb204_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxH1]
	movapd xmm4, [esp + nb204_fjyH1]
	movapd xmm5, [esp + nb204_fjzH1]
	mulsd xmm0, [esp + nb204_dxH1H1]
	mulsd xmm1, [esp + nb204_dyH1H1]
	mulsd xmm2, [esp + nb204_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixH1]
	addsd xmm1, [esp + nb204_fiyH1]
	addsd xmm2, [esp + nb204_fizH1]
	movlpd [esp + nb204_fjxH1], xmm3
	movlpd [esp + nb204_fjyH1], xmm4
	movlpd [esp + nb204_fjzH1], xmm5
	movlpd [esp + nb204_fixH1], xmm0
	movlpd [esp + nb204_fiyH1], xmm1
	movlpd [esp + nb204_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb204_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqH1H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulsd xmm0, xmm0
	subsd  xmm4, [esp + nb204_crf]
	mulsd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb204_fjxH2]
	movapd xmm4, [esp + nb204_fjyH2]
	movapd xmm5, [esp + nb204_fjzH2]
	mulsd xmm0, [esp + nb204_dxH1H2]
	mulsd xmm1, [esp + nb204_dyH1H2]
	mulsd xmm2, [esp + nb204_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixH1]
	addsd xmm1, [esp + nb204_fiyH1]
	addsd xmm2, [esp + nb204_fizH1]
	movlpd [esp + nb204_fjxH2], xmm3
	movlpd [esp + nb204_fjyH2], xmm4
	movlpd [esp + nb204_fjzH2], xmm5
	movlpd [esp + nb204_fixH1], xmm0
	movlpd [esp + nb204_fiyH1], xmm1
	movlpd [esp + nb204_fizH1], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb204_rinvH2M]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqH2M] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb204_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb204_qqMH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqMH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxM]
	movapd xmm4, [esp + nb204_fjyM]
	movapd xmm5, [esp + nb204_fjzM]
	mulsd xmm0, [esp + nb204_dxH2M]
	mulsd xmm1, [esp + nb204_dyH2M]
	mulsd xmm2, [esp + nb204_dzH2M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixH2]
	addsd xmm1, [esp + nb204_fiyH2]
	addsd xmm2, [esp + nb204_fizH2]
	movlpd [esp + nb204_fjxM], xmm3
	movlpd [esp + nb204_fjyM], xmm4
	movlpd [esp + nb204_fjzM], xmm5
	movlpd [esp + nb204_fixH2], xmm0
	movlpd [esp + nb204_fiyH2], xmm1
	movlpd [esp + nb204_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb204_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqH2H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb204_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb204_fjxH1]
	movapd xmm4, [esp + nb204_fjyH1]
	movapd xmm5, [esp + nb204_fjzH1]
	mulsd xmm0, [esp + nb204_dxH2H1]
	mulsd xmm1, [esp + nb204_dyH2H1]
	mulsd xmm2, [esp + nb204_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixH2]
	addsd xmm1, [esp + nb204_fiyH2]
	addsd xmm2, [esp + nb204_fizH2]
	movlpd [esp + nb204_fjxH1], xmm3
	movlpd [esp + nb204_fjyH1], xmm4
	movlpd [esp + nb204_fjzH1], xmm5
	movlpd [esp + nb204_fixH2], xmm0
	movlpd [esp + nb204_fiyH2], xmm1
	movlpd [esp + nb204_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb204_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb204_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb204_rsqH2H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb204_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb204_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb204_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb204_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsMH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm1, xmm0
	movlpd [esp + nb204_vctot], xmm6
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb204_fjxH2]
	movapd xmm4, [esp + nb204_fjyH2]
	movapd xmm5, [esp + nb204_fjzH2]
	mulsd xmm0, [esp + nb204_dxH2H2]
	mulsd xmm1, [esp + nb204_dyH2H2]
	mulsd xmm2, [esp + nb204_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb204_fixH2]
	addsd xmm1, [esp + nb204_fiyH2]
	addsd xmm2, [esp + nb204_fizH2]
	movlpd [esp + nb204_fjxH2], xmm3
	movlpd [esp + nb204_fjyH2], xmm4
	movlpd [esp + nb204_fjzH2], xmm5
	movlpd [esp + nb204_fixH2], xmm0
	movlpd [esp + nb204_fiyH2], xmm1
	movlpd [esp + nb204_fizH2], xmm2

	mov edi, [ebp + nb204_faction]
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8 + 24]
	movlpd xmm1, [edi + eax*8 + 32]
	movlpd xmm2, [edi + eax*8 + 40]
	movlpd xmm3, [edi + eax*8 + 48]
	movlpd xmm4, [edi + eax*8 + 56]
	movlpd xmm5, [edi + eax*8 + 64]
	movlpd xmm6, [edi + eax*8 + 72]
	movlpd xmm7, [edi + eax*8 + 80]
	addsd xmm0, [esp + nb204_fjxH1]
	addsd xmm1, [esp + nb204_fjyH1]
	addsd xmm2, [esp + nb204_fjzH1]
	addsd xmm3, [esp + nb204_fjxH2]
	addsd xmm4, [esp + nb204_fjyH2]
	addsd xmm5, [esp + nb204_fjzH2]
	addsd xmm6, [esp + nb204_fjxM]
	addsd xmm7, [esp + nb204_fjyM]
	movlpd [edi + eax*8 + 24], xmm0
	movlpd [edi + eax*8 + 32], xmm1
	movlpd [edi + eax*8 + 40], xmm2
	movlpd [edi + eax*8 + 48], xmm3
	movlpd [edi + eax*8 + 56], xmm4
	movlpd [edi + eax*8 + 64], xmm5
	movlpd [edi + eax*8 + 72], xmm6
	movlpd [edi + eax*8 + 80], xmm7

	movlpd xmm0, [edi + eax*8 + 88]
	addsd xmm0, [esp + nb204_fjzM]
	movlpd [edi + eax*8 + 88], xmm0
	
.nb204_updateouterdata:
	mov   ecx, [esp + nb204_ii3]
	mov   edi, [ebp + nb204_faction]
	mov   esi, [ebp + nb204_fshift]
	mov   edx, [esp + nb204_is3]

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb204_fixH1]
	movapd xmm1, [esp + nb204_fiyH1]
	movapd xmm2, [esp + nb204_fizH1]

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
	movapd xmm0, [esp + nb204_fixH2]
	movapd xmm1, [esp + nb204_fiyH2]
	movapd xmm2, [esp + nb204_fizH2]

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
	movapd xmm0, [esp + nb204_fixM]
	movapd xmm1, [esp + nb204_fiyM]
	movapd xmm2, [esp + nb204_fizM]

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
	mov esi, [esp + nb204_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb204_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb204_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb204_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb204_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb204_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb204_n], esi
        jmp .nb204_outer
.nb204_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb204_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb204_end
        ;# non-zero, do one more workunit
        jmp   .nb204_threadloop
.nb204_end:
	emms

	mov eax, [esp + nb204_nouter]
	mov ebx, [esp + nb204_ninner]
	mov ecx, [ebp + nb204_outeriter]
	mov edx, [ebp + nb204_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb204_salign]
	add esp, eax
	add esp, 1480
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	

	
.globl nb_kernel204nf_ia32_sse2
.globl _nb_kernel204nf_ia32_sse2
nb_kernel204nf_ia32_sse2:	
_nb_kernel204nf_ia32_sse2:	
.equiv          nb204nf_p_nri,          8
.equiv          nb204nf_iinr,           12
.equiv          nb204nf_jindex,         16
.equiv          nb204nf_jjnr,           20
.equiv          nb204nf_shift,          24
.equiv          nb204nf_shiftvec,       28
.equiv          nb204nf_fshift,         32
.equiv          nb204nf_gid,            36
.equiv          nb204nf_pos,            40
.equiv          nb204nf_faction,        44
.equiv          nb204nf_charge,         48
.equiv          nb204nf_p_facel,        52
.equiv          nb204nf_argkrf,         56
.equiv          nb204nf_argcrf,         60
.equiv          nb204nf_Vc,             64
.equiv          nb204nf_type,           68
.equiv          nb204nf_p_ntype,        72
.equiv          nb204nf_vdwparam,       76
.equiv          nb204nf_Vvdw,           80
.equiv          nb204nf_p_tabscale,     84
.equiv          nb204nf_VFtab,          88
.equiv          nb204nf_invsqrta,       92
.equiv          nb204nf_dvda,           96
.equiv          nb204nf_p_gbtabscale,   100
.equiv          nb204nf_GBtab,          104
.equiv          nb204nf_p_nthreads,     108
.equiv          nb204nf_count,          112
.equiv          nb204nf_mtx,            116
.equiv          nb204nf_outeriter,      120
.equiv          nb204nf_inneriter,      124
.equiv          nb204nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb204nf_ixM,            0
.equiv          nb204nf_iyM,            16
.equiv          nb204nf_izM,            32
.equiv          nb204nf_ixH1,           48
.equiv          nb204nf_iyH1,           64
.equiv          nb204nf_izH1,           80
.equiv          nb204nf_ixH2,           96
.equiv          nb204nf_iyH2,           112
.equiv          nb204nf_izH2,           128
.equiv          nb204nf_jxM,            144
.equiv          nb204nf_jyM,            160
.equiv          nb204nf_jzM,            176
.equiv          nb204nf_jxH1,           192
.equiv          nb204nf_jyH1,           208
.equiv          nb204nf_jzH1,           224
.equiv          nb204nf_jxH2,           240
.equiv          nb204nf_jyH2,           256
.equiv          nb204nf_jzH2,           272
.equiv          nb204nf_qqMM,           288
.equiv          nb204nf_qqMH,           304
.equiv          nb204nf_qqHH,           320
.equiv          nb204nf_vctot,          336
.equiv          nb204nf_half,           352
.equiv          nb204nf_three,          368
.equiv          nb204nf_rsqMM,          384
.equiv          nb204nf_rsqMH1,         400
.equiv          nb204nf_rsqMH2,         416
.equiv          nb204nf_rsqH1M,         432
.equiv          nb204nf_rsqH1H1,        448
.equiv          nb204nf_rsqH1H2,        464
.equiv          nb204nf_rsqH2M,         480
.equiv          nb204nf_rsqH2H1,        496
.equiv          nb204nf_rsqH2H2,        512
.equiv          nb204nf_rinvMM,         528
.equiv          nb204nf_rinvMH1,        544
.equiv          nb204nf_rinvMH2,        560
.equiv          nb204nf_rinvH1M,        576
.equiv          nb204nf_rinvH1H1,       592
.equiv          nb204nf_rinvH1H2,       608
.equiv          nb204nf_rinvH2M,        624
.equiv          nb204nf_rinvH2H1,       640
.equiv          nb204nf_rinvH2H2,       656
.equiv          nb204nf_krf,            672
.equiv          nb204nf_crf,            688
.equiv          nb204nf_is3,            704
.equiv          nb204nf_ii3,            708
.equiv          nb204nf_innerjjnr,      712
.equiv          nb204nf_innerk,         716
.equiv          nb204nf_n,              720
.equiv          nb204nf_nn1,            724
.equiv          nb204nf_nri,            728
.equiv          nb204nf_nouter,         732
.equiv          nb204nf_ninner,         736
.equiv          nb204nf_salign,         740
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 744		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb204nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb204nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb204nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb204nf_nouter], eax
	mov [esp + nb204nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb204nf_half], eax
	mov [esp + nb204nf_half+4], ebx
	movsd xmm1, [esp + nb204nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb204nf_half], xmm1
	movapd [esp + nb204nf_three], xmm3

	mov esi, [ebp + nb204nf_argkrf]
	mov edi, [ebp + nb204nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb204nf_krf], xmm5
	movapd [esp + nb204nf_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb204nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb204nf_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb204nf_p_facel]
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
	movapd [esp + nb204nf_qqMM], xmm3
	movapd [esp + nb204nf_qqMH], xmm4
	movapd [esp + nb204nf_qqHH], xmm5
	
.nb204nf_threadloop:
        mov   esi, [ebp + nb204nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb204nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb204nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb204nf_n], eax
        mov [esp + nb204nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204nf_outerstart
        jmp .nb204nf_end

.nb204nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb204nf_nouter]
	mov [esp + nb204nf_nouter], ebx

.nb204nf_outer:
	mov   eax, [ebp + nb204nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb204nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb204nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb204nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb204nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb204nf_ixH1], xmm3
	movapd [esp + nb204nf_iyH1], xmm4
	movapd [esp + nb204nf_izH1], xmm5

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
	movapd [esp + nb204nf_ixH2], xmm0
	movapd [esp + nb204nf_iyH2], xmm1
	movapd [esp + nb204nf_izH2], xmm2
	movapd [esp + nb204nf_ixM], xmm3
	movapd [esp + nb204nf_iyM], xmm4
	movapd [esp + nb204nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb204nf_vctot], xmm4
	
	mov   eax, [ebp + nb204nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb204nf_pos]
	mov   eax, [ebp + nb204nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb204nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb204nf_ninner]
	mov   [esp + nb204nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb204nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb204nf_unroll_loop
	jmp   .nb204nf_checksingle
.nb204nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb204nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb204nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb204nf_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb204nf_jxH1], xmm2
	movapd 	[esp + nb204nf_jyH1], xmm3
	movapd 	[esp + nb204nf_jzH1], xmm4
	movapd 	[esp + nb204nf_jxH2], xmm5
	movapd 	[esp + nb204nf_jyH2], xmm6
	movapd 	[esp + nb204nf_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movhpd xmm2, [esi + ebx*8 + 72]
	movhpd xmm3, [esi + ebx*8 + 80]
	movhpd xmm4, [esi + ebx*8 + 88]
	movapd 	[esp + nb204nf_jxM], xmm2
	movapd 	[esp + nb204nf_jyM], xmm3
	movapd 	[esp + nb204nf_jzM], xmm4
	
	movapd xmm0, [esp + nb204nf_ixM]
	movapd xmm1, [esp + nb204nf_iyM]
	movapd xmm2, [esp + nb204nf_izM]
	movapd xmm3, [esp + nb204nf_ixM]
	movapd xmm4, [esp + nb204nf_iyM]
	movapd xmm5, [esp + nb204nf_izM]
	subpd  xmm0, [esp + nb204nf_jxM]
	subpd  xmm1, [esp + nb204nf_jyM]
	subpd  xmm2, [esp + nb204nf_jzM]
	subpd  xmm3, [esp + nb204nf_jxH1]
	subpd  xmm4, [esp + nb204nf_jyH1]
	subpd  xmm5, [esp + nb204nf_jzH1]
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
	movapd [esp + nb204nf_rsqMM], xmm0
	movapd [esp + nb204nf_rsqMH1], xmm3

	movapd xmm0, [esp + nb204nf_ixM]
	movapd xmm1, [esp + nb204nf_iyM]
	movapd xmm2, [esp + nb204nf_izM]
	movapd xmm3, [esp + nb204nf_ixH1]
	movapd xmm4, [esp + nb204nf_iyH1]
	movapd xmm5, [esp + nb204nf_izH1]
	subpd  xmm0, [esp + nb204nf_jxH2]
	subpd  xmm1, [esp + nb204nf_jyH2]
	subpd  xmm2, [esp + nb204nf_jzH2]
	subpd  xmm3, [esp + nb204nf_jxM]
	subpd  xmm4, [esp + nb204nf_jyM]
	subpd  xmm5, [esp + nb204nf_jzM]
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
	movapd [esp + nb204nf_rsqMH2], xmm0
	movapd [esp + nb204nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb204nf_ixH1]
	movapd xmm1, [esp + nb204nf_iyH1]
	movapd xmm2, [esp + nb204nf_izH1]
	movapd xmm3, [esp + nb204nf_ixH1]
	movapd xmm4, [esp + nb204nf_iyH1]
	movapd xmm5, [esp + nb204nf_izH1]
	subpd  xmm0, [esp + nb204nf_jxH1]
	subpd  xmm1, [esp + nb204nf_jyH1]
	subpd  xmm2, [esp + nb204nf_jzH1]
	subpd  xmm3, [esp + nb204nf_jxH2]
	subpd  xmm4, [esp + nb204nf_jyH2]
	subpd  xmm5, [esp + nb204nf_jzH2]
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
	movapd [esp + nb204nf_rsqH1H1], xmm0
	movapd [esp + nb204nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb204nf_ixH2]
	movapd xmm1, [esp + nb204nf_iyH2]
	movapd xmm2, [esp + nb204nf_izH2]
	movapd xmm3, [esp + nb204nf_ixH2]
	movapd xmm4, [esp + nb204nf_iyH2]
	movapd xmm5, [esp + nb204nf_izH2]
	subpd  xmm0, [esp + nb204nf_jxM]
	subpd  xmm1, [esp + nb204nf_jyM]
	subpd  xmm2, [esp + nb204nf_jzM]
	subpd  xmm3, [esp + nb204nf_jxH1]
	subpd  xmm4, [esp + nb204nf_jyH1]
	subpd  xmm5, [esp + nb204nf_jzH1]
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
	movapd [esp + nb204nf_rsqH2M], xmm0
	movapd [esp + nb204nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb204nf_ixH2]
	movapd xmm1, [esp + nb204nf_iyH2]
	movapd xmm2, [esp + nb204nf_izH2]
	subpd  xmm0, [esp + nb204nf_jxH2]
	subpd  xmm1, [esp + nb204nf_jyH2]
	subpd  xmm2, [esp + nb204nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb204nf_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
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
	movapd  xmm3, [esp + nb204nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb204nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvH2H2], xmm1
	movapd [esp + nb204nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb204nf_rsqMM]
	movapd xmm4, [esp + nb204nf_rsqMH1]	
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
	movapd  xmm3, [esp + nb204nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb204nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb204nf_half] ;# rinv
	movapd [esp + nb204nf_rinvMM], xmm1
	movapd [esp + nb204nf_rinvMH1], xmm5

	movapd xmm0, [esp + nb204nf_rsqMH2]
	movapd xmm4, [esp + nb204nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb204nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb204nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvMH2], xmm1
	movapd [esp + nb204nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb204nf_rsqH1H1]
	movapd xmm4, [esp + nb204nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb204nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb204nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb204nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvH1H1], xmm1
	movapd [esp + nb204nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb204nf_rsqH2M]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb204nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb204nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb204nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvH2M], xmm1
	
	;# start with MM interaction 
	movapd xmm6, [esp + nb204nf_krf]
	mulpd  xmm6, [esp + nb204nf_rsqMM]	;# xmm5=krsq 
	addpd  xmm6, [esp + nb204nf_rinvMM]	;# xmm6=rinv+ krsq 
	subpd  xmm6, [esp + nb204nf_crf]
	
	mulpd  xmm6, [esp + nb204nf_qqMM] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, [esp + nb204nf_vctot] ;# local vctot summation variable 

	;# M-H1 interaction 
	movapd xmm5, [esp + nb204nf_krf]
	mulpd  xmm5, [esp + nb204nf_rsqMH1]	;# xmm5=krsq 
	addpd  xmm5, [esp + nb204nf_rinvMH1]	;# xmm6=rinv+ krsq 
	subpd  xmm5, [esp + nb204nf_crf]
	
	mulpd  xmm5, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm5 ;# local vctot summation variable 

	;# M-H2 interaction 
	movapd xmm7, [esp + nb204nf_krf]
	mulpd  xmm7, [esp + nb204nf_rsqMH2]	;# xmm5=krsq 
	addpd  xmm7, [esp + nb204nf_rinvMH2]	;# xmm6=rinv+ krsq 
	subpd  xmm7, [esp + nb204nf_crf]
	
	mulpd  xmm7, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm7 ;# local vctot summation variable 

	;# H1-M interaction 
	movapd xmm4, [esp + nb204nf_krf]
	mulpd  xmm4, [esp + nb204nf_rsqH1M]	;# xmm5=krsq 
	addpd  xmm4, [esp + nb204nf_rinvH1M]	;# xmm6=rinv+ krsq 
	subpd  xmm4, [esp + nb204nf_crf]
	
	mulpd  xmm4, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H1-H1 interaction 
	movapd xmm5, [esp + nb204nf_krf]
	mulpd  xmm5, [esp + nb204nf_rsqH1H1]	;# xmm5=krsq 
	addpd  xmm5, [esp + nb204nf_rinvH1H1]	;# xmm6=rinv+ krsq 
	subpd  xmm5, [esp + nb204nf_crf]
	
	mulpd  xmm5, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm5 ;# local vctot summation variable 

	;# H1-H2 interaction 
	movapd xmm7, [esp + nb204nf_krf]
	mulpd  xmm7, [esp + nb204nf_rsqH1H2]	;# xmm5=krsq 
	addpd  xmm7, [esp + nb204nf_rinvH1H2]	;# xmm6=rinv+ krsq 
	subpd  xmm7, [esp + nb204nf_crf]
	
	mulpd  xmm7, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm7 ;# local vctot summation variable 

	;# H2-M interaction 
	movapd xmm4, [esp + nb204nf_krf]
	mulpd  xmm4, [esp + nb204nf_rsqH2M]	;# xmm5=krsq 
	addpd  xmm4, [esp + nb204nf_rinvH2M]	;# xmm6=rinv+ krsq 
	subpd  xmm4, [esp + nb204nf_crf]
	
	mulpd  xmm4, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H2-H1 interaction 
	movapd xmm5, [esp + nb204nf_krf]
	mulpd  xmm5, [esp + nb204nf_rsqH2H1]	;# xmm5=krsq 
	addpd  xmm5, [esp + nb204nf_rinvH2H1]	;# xmm6=rinv+ krsq 
	subpd  xmm5, [esp + nb204nf_crf]
	
	mulpd  xmm5, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm5 ;# local vctot summation variable 

	;# H2-H2 interaction 
	movapd xmm7, [esp + nb204nf_krf]
	mulpd  xmm7, [esp + nb204nf_rsqH2H2]	;# xmm5=krsq 
	addpd  xmm7, [esp + nb204nf_rinvH2H2]	;# xmm6=rinv+ krsq 
	subpd  xmm7, [esp + nb204nf_crf]
	
	mulpd  xmm7, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm7 ;# local vctot summation variable 
	movapd [esp + nb204nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [esp + nb204nf_innerk],  2
	jl    .nb204nf_checksingle
	jmp   .nb204nf_unroll_loop
.nb204nf_checksingle:
	mov   edx, [esp + nb204nf_innerk]
	and   edx, 1
	jnz   .nb204nf_dosingle
	jmp   .nb204nf_updateouterdata
.nb204nf_dosingle:
	mov   edx, [esp + nb204nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]
	
	mov esi, [ebp + nb204nf_pos]
	lea   eax, [eax + eax*2]  

	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movapd 	[esp + nb204nf_jxH1], xmm2
	movapd 	[esp + nb204nf_jyH1], xmm3
	movapd 	[esp + nb204nf_jzH1], xmm4
	movapd 	[esp + nb204nf_jxH2], xmm5
	movapd 	[esp + nb204nf_jyH2], xmm6
	movapd 	[esp + nb204nf_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movapd 	[esp + nb204nf_jxM], xmm2
	movapd 	[esp + nb204nf_jyM], xmm3
	movapd 	[esp + nb204nf_jzM], xmm4
	
	movapd xmm0, [esp + nb204nf_ixM]
	movapd xmm1, [esp + nb204nf_iyM]
	movapd xmm2, [esp + nb204nf_izM]
	movapd xmm3, [esp + nb204nf_ixM]
	movapd xmm4, [esp + nb204nf_iyM]
	movapd xmm5, [esp + nb204nf_izM]
	subsd  xmm0, [esp + nb204nf_jxM]
	subsd  xmm1, [esp + nb204nf_jyM]
	subsd  xmm2, [esp + nb204nf_jzM]
	subsd  xmm3, [esp + nb204nf_jxH1]
	subsd  xmm4, [esp + nb204nf_jyH1]
	subsd  xmm5, [esp + nb204nf_jzH1]
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
	movapd [esp + nb204nf_rsqMM], xmm0
	movapd [esp + nb204nf_rsqMH1], xmm3

	movapd xmm0, [esp + nb204nf_ixM]
	movapd xmm1, [esp + nb204nf_iyM]
	movapd xmm2, [esp + nb204nf_izM]
	movapd xmm3, [esp + nb204nf_ixH1]
	movapd xmm4, [esp + nb204nf_iyH1]
	movapd xmm5, [esp + nb204nf_izH1]
	subsd  xmm0, [esp + nb204nf_jxH2]
	subsd  xmm1, [esp + nb204nf_jyH2]
	subsd  xmm2, [esp + nb204nf_jzH2]
	subsd  xmm3, [esp + nb204nf_jxM]
	subsd  xmm4, [esp + nb204nf_jyM]
	subsd  xmm5, [esp + nb204nf_jzM]
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
	movapd [esp + nb204nf_rsqMH2], xmm0
	movapd [esp + nb204nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb204nf_ixH1]
	movapd xmm1, [esp + nb204nf_iyH1]
	movapd xmm2, [esp + nb204nf_izH1]
	movapd xmm3, [esp + nb204nf_ixH1]
	movapd xmm4, [esp + nb204nf_iyH1]
	movapd xmm5, [esp + nb204nf_izH1]
	subsd  xmm0, [esp + nb204nf_jxH1]
	subsd  xmm1, [esp + nb204nf_jyH1]
	subsd  xmm2, [esp + nb204nf_jzH1]
	subsd  xmm3, [esp + nb204nf_jxH2]
	subsd  xmm4, [esp + nb204nf_jyH2]
	subsd  xmm5, [esp + nb204nf_jzH2]
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
	movapd [esp + nb204nf_rsqH1H1], xmm0
	movapd [esp + nb204nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb204nf_ixH2]
	movapd xmm1, [esp + nb204nf_iyH2]
	movapd xmm2, [esp + nb204nf_izH2]
	movapd xmm3, [esp + nb204nf_ixH2]
	movapd xmm4, [esp + nb204nf_iyH2]
	movapd xmm5, [esp + nb204nf_izH2]
	subsd  xmm0, [esp + nb204nf_jxM]
	subsd  xmm1, [esp + nb204nf_jyM]
	subsd  xmm2, [esp + nb204nf_jzM]
	subsd  xmm3, [esp + nb204nf_jxH1]
	subsd  xmm4, [esp + nb204nf_jyH1]
	subsd  xmm5, [esp + nb204nf_jzH1]
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
	movapd [esp + nb204nf_rsqH2M], xmm0
	movapd [esp + nb204nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb204nf_ixH2]
	movapd xmm1, [esp + nb204nf_iyH2]
	movapd xmm2, [esp + nb204nf_izH2]
	subsd  xmm0, [esp + nb204nf_jxH2]
	subsd  xmm1, [esp + nb204nf_jyH2]
	subsd  xmm2, [esp + nb204nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb204nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb204nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb204nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvH2H2], xmm1
	movapd [esp + nb204nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb204nf_rsqMM]
	movapd xmm4, [esp + nb204nf_rsqMH1]	
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
	movapd  xmm3, [esp + nb204nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb204nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb204nf_half] ;# rinv
	movapd [esp + nb204nf_rinvMM], xmm1
	movapd [esp + nb204nf_rinvMH1], xmm5

	movapd xmm0, [esp + nb204nf_rsqMH2]
	movapd xmm4, [esp + nb204nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb204nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb204nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvMH2], xmm1
	movapd [esp + nb204nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb204nf_rsqH1H1]
	movapd xmm4, [esp + nb204nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb204nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb204nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb204nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb204nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb204nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvH1H1], xmm1
	movapd [esp + nb204nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb204nf_rsqH2M]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb204nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb204nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb204nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb204nf_half] ;# rinv 
	movapd [esp + nb204nf_rinvH2M], xmm1
	
	;# start with MM interaction 
	movapd xmm6, [esp + nb204nf_krf]
	mulsd  xmm6, [esp + nb204nf_rsqMM]	;# xmm5=krsq 
	addsd  xmm6, [esp + nb204nf_rinvMM]	;# xmm6=rinv+ krsq 
	subsd  xmm6, [esp + nb204nf_crf]
	
	mulsd  xmm6, [esp + nb204nf_qqMM] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, [esp + nb204nf_vctot] ;# local vctot summation variable 

	;# M-H1 interaction 
	movapd xmm5, [esp + nb204nf_krf]
	mulsd  xmm5, [esp + nb204nf_rsqMH1]	;# xmm5=krsq 
	addsd  xmm5, [esp + nb204nf_rinvMH1]	;# xmm6=rinv+ krsq 
	subsd  xmm5, [esp + nb204nf_crf]
	
	mulsd  xmm5, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm5 ;# local vctot summation variable 

	;# M-H2 interaction 
	movapd xmm7, [esp + nb204nf_krf]
	mulsd  xmm7, [esp + nb204nf_rsqMH2]	;# xmm5=krsq 
	addsd  xmm7, [esp + nb204nf_rinvMH2]	;# xmm6=rinv+ krsq 
	subsd  xmm7, [esp + nb204nf_crf]
	
	mulsd  xmm7, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm7 ;# local vctot summation variable 

	;# H1-M interaction 
	movapd xmm4, [esp + nb204nf_krf]
	mulsd  xmm4, [esp + nb204nf_rsqH1M]	;# xmm5=krsq 
	addsd  xmm4, [esp + nb204nf_rinvH1M]	;# xmm6=rinv+ krsq 
	subsd  xmm4, [esp + nb204nf_crf]
	
	mulsd  xmm4, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H1-H1 interaction 
	movapd xmm5, [esp + nb204nf_krf]
	mulsd  xmm5, [esp + nb204nf_rsqH1H1]	;# xmm5=krsq 
	addsd  xmm5, [esp + nb204nf_rinvH1H1]	;# xmm6=rinv+ krsq 
	subsd  xmm5, [esp + nb204nf_crf]
	
	mulsd  xmm5, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm5 ;# local vctot summation variable 

	;# H1-H2 interaction 
	movapd xmm7, [esp + nb204nf_krf]
	mulsd  xmm7, [esp + nb204nf_rsqH1H2]	;# xmm5=krsq 
	addsd  xmm7, [esp + nb204nf_rinvH1H2]	;# xmm6=rinv+ krsq 
	subsd  xmm7, [esp + nb204nf_crf]
	
	mulsd  xmm7, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm7 ;# local vctot summation variable 

	;# H2-M interaction 
	movapd xmm4, [esp + nb204nf_krf]
	mulsd  xmm4, [esp + nb204nf_rsqH2M]	;# xmm5=krsq 
	addsd  xmm4, [esp + nb204nf_rinvH2M]	;# xmm6=rinv+ krsq 
	subsd  xmm4, [esp + nb204nf_crf]
	
	mulsd  xmm4, [esp + nb204nf_qqMH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H2-H1 interaction 
	movapd xmm5, [esp + nb204nf_krf]
	mulsd  xmm5, [esp + nb204nf_rsqH2H1]	;# xmm5=krsq 
	addsd  xmm5, [esp + nb204nf_rinvH2H1]	;# xmm6=rinv+ krsq 
	subsd  xmm5, [esp + nb204nf_crf]
	
	mulsd  xmm5, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm5 ;# local vctot summation variable 

	;# H2-H2 interaction 
	movapd xmm7, [esp + nb204nf_krf]
	mulsd  xmm7, [esp + nb204nf_rsqH2H2]	;# xmm5=krsq 
	addsd  xmm7, [esp + nb204nf_rinvH2H2]	;# xmm6=rinv+ krsq 
	subsd  xmm7, [esp + nb204nf_crf]
	
	mulsd  xmm7, [esp + nb204nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm7 ;# local vctot summation variable 
	movlpd [esp + nb204nf_vctot], xmm6
	
.nb204nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb204nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb204nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb204nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb204nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb204nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb204nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb204nf_n], esi
        jmp .nb204nf_outer
.nb204nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb204nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb204nf_end
        ;# non-zero, do one more workunit
        jmp   .nb204nf_threadloop
.nb204nf_end:
	emms

	mov eax, [esp + nb204nf_nouter]
	mov ebx, [esp + nb204nf_ninner]
	mov ecx, [ebp + nb204nf_outeriter]
	mov edx, [ebp + nb204nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb204nf_salign]
	add esp, eax
	add esp, 744
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

	
