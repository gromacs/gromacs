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


		
.globl nb_kernel304_ia32_sse2
.globl _nb_kernel304_ia32_sse2
nb_kernel304_ia32_sse2:	
_nb_kernel304_ia32_sse2:	
.equiv          nb304_p_nri,            8
.equiv          nb304_iinr,             12
.equiv          nb304_jindex,           16
.equiv          nb304_jjnr,             20
.equiv          nb304_shift,            24
.equiv          nb304_shiftvec,         28
.equiv          nb304_fshift,           32
.equiv          nb304_gid,              36
.equiv          nb304_pos,              40
.equiv          nb304_faction,          44
.equiv          nb304_charge,           48
.equiv          nb304_p_facel,          52
.equiv          nb304_argkrf,           56
.equiv          nb304_argcrf,           60
.equiv          nb304_Vc,               64
.equiv          nb304_type,             68
.equiv          nb304_p_ntype,          72
.equiv          nb304_vdwparam,         76
.equiv          nb304_Vvdw,             80
.equiv          nb304_p_tabscale,       84
.equiv          nb304_VFtab,            88
.equiv          nb304_invsqrta,         92
.equiv          nb304_dvda,             96
.equiv          nb304_p_gbtabscale,     100
.equiv          nb304_GBtab,            104
.equiv          nb304_p_nthreads,       108
.equiv          nb304_count,            112
.equiv          nb304_mtx,              116
.equiv          nb304_outeriter,        120
.equiv          nb304_inneriter,        124
.equiv          nb304_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb304_ixM,              0
.equiv          nb304_iyM,              16
.equiv          nb304_izM,              32
.equiv          nb304_ixH1,             48
.equiv          nb304_iyH1,             64
.equiv          nb304_izH1,             80
.equiv          nb304_ixH2,             96
.equiv          nb304_iyH2,             112
.equiv          nb304_izH2,             128
.equiv          nb304_jxM,              144
.equiv          nb304_jyM,              160
.equiv          nb304_jzM,              176
.equiv          nb304_jxH1,             192
.equiv          nb304_jyH1,             208
.equiv          nb304_jzH1,             224
.equiv          nb304_jxH2,             240
.equiv          nb304_jyH2,             256
.equiv          nb304_jzH2,             272
.equiv          nb304_dxMM,             288
.equiv          nb304_dyMM,             304
.equiv          nb304_dzMM,             320
.equiv          nb304_dxMH1,            336
.equiv          nb304_dyMH1,            352
.equiv          nb304_dzMH1,            368
.equiv          nb304_dxMH2,            384
.equiv          nb304_dyMH2,            400
.equiv          nb304_dzMH2,            416
.equiv          nb304_dxH1M,            432
.equiv          nb304_dyH1M,            448
.equiv          nb304_dzH1M,            464
.equiv          nb304_dxH1H1,           480
.equiv          nb304_dyH1H1,           496
.equiv          nb304_dzH1H1,           512
.equiv          nb304_dxH1H2,           528
.equiv          nb304_dyH1H2,           544
.equiv          nb304_dzH1H2,           560
.equiv          nb304_dxH2M,            576
.equiv          nb304_dyH2M,            592
.equiv          nb304_dzH2M,            608
.equiv          nb304_dxH2H1,           624
.equiv          nb304_dyH2H1,           640
.equiv          nb304_dzH2H1,           656
.equiv          nb304_dxH2H2,           672
.equiv          nb304_dyH2H2,           688
.equiv          nb304_dzH2H2,           704
.equiv          nb304_qqMM,             720
.equiv          nb304_qqMH,             736
.equiv          nb304_qqHH,             752
.equiv          nb304_two,              768
.equiv          nb304_tsc,              784
.equiv          nb304_vctot,            800
.equiv          nb304_fixM,             816
.equiv          nb304_fiyM,             832
.equiv          nb304_fizM,             848
.equiv          nb304_fixH1,            864
.equiv          nb304_fiyH1,            880
.equiv          nb304_fizH1,            896
.equiv          nb304_fixH2,            912
.equiv          nb304_fiyH2,            928
.equiv          nb304_fizH2,            944
.equiv          nb304_fjxM,             960
.equiv          nb304_fjyM,             976
.equiv          nb304_fjzM,             992
.equiv          nb304_fjxH1,            1008
.equiv          nb304_fjyH1,            1024
.equiv          nb304_fjzH1,            1040
.equiv          nb304_fjxH2,            1056
.equiv          nb304_fjyH2,            1072
.equiv          nb304_fjzH2,            1088
.equiv          nb304_half,             1104
.equiv          nb304_three,            1120
.equiv          nb304_rsqMM,            1136
.equiv          nb304_rsqMH1,           1152
.equiv          nb304_rsqMH2,           1168
.equiv          nb304_rsqH1M,           1184
.equiv          nb304_rsqH1H1,          1200
.equiv          nb304_rsqH1H2,          1216
.equiv          nb304_rsqH2M,           1232
.equiv          nb304_rsqH2H1,          1248
.equiv          nb304_rsqH2H2,          1264
.equiv          nb304_rinvMM,           1280
.equiv          nb304_rinvMH1,          1296
.equiv          nb304_rinvMH2,          1312
.equiv          nb304_rinvH1M,          1328
.equiv          nb304_rinvH1H1,         1344
.equiv          nb304_rinvH1H2,         1360
.equiv          nb304_rinvH2M,          1376
.equiv          nb304_rinvH2H1,         1392
.equiv          nb304_rinvH2H2,         1408
.equiv          nb304_is3,              1424
.equiv          nb304_ii3,              1428
.equiv          nb304_innerjjnr,        1432
.equiv          nb304_innerk,           1436
.equiv          nb304_n,                1440
.equiv          nb304_nn1,              1444
.equiv          nb304_nri,              1448
.equiv          nb304_nouter,           1452
.equiv          nb304_ninner,           1456
.equiv          nb304_salign,           1460
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 1464		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb304_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb304_p_nri]
	mov ecx, [ecx]
	mov [esp + nb304_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb304_nouter], eax
	mov [esp + nb304_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb304_half], eax
	mov [esp + nb304_half+4], ebx
	movsd xmm1, [esp + nb304_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb304_half], xmm1
	movapd [esp + nb304_two], xmm2
	movapd [esp + nb304_three], xmm3

	mov eax, [ebp + nb304_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb304_tsc],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb304_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb304_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb304_p_facel]
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
	movapd [esp + nb304_qqMM], xmm3
	movapd [esp + nb304_qqMH], xmm4
	movapd [esp + nb304_qqHH], xmm5		

.nb304_threadloop:
        mov   esi, [ebp + nb304_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb304_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb304_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb304_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb304_n], eax
        mov [esp + nb304_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb304_outerstart
        jmp .nb304_end

.nb304_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb304_nouter]
	mov [esp + nb304_nouter], ebx

.nb304_outer:
	mov   eax, [ebp + nb304_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax +esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb304_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb304_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb304_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb304_pos]    ;# eax = base of pos[]  
	mov   [esp + nb304_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb304_ixH1], xmm3
	movapd [esp + nb304_iyH1], xmm4
	movapd [esp + nb304_izH1], xmm5

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
	movapd [esp + nb304_ixH2], xmm0
	movapd [esp + nb304_iyH2], xmm1
	movapd [esp + nb304_izH2], xmm2
	movapd [esp + nb304_ixM], xmm3
	movapd [esp + nb304_iyM], xmm4
	movapd [esp + nb304_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb304_vctot], xmm4
	movapd [esp + nb304_fixM], xmm4
	movapd [esp + nb304_fiyM], xmm4
	movapd [esp + nb304_fizM], xmm4
	movapd [esp + nb304_fixH1], xmm4
	movapd [esp + nb304_fiyH1], xmm4
	movapd [esp + nb304_fizH1], xmm4
	movapd [esp + nb304_fixH2], xmm4
	movapd [esp + nb304_fiyH2], xmm4
	movapd [esp + nb304_fizH2], xmm4
	
	mov   eax, [ebp + nb304_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb304_pos]
	mov   edi, [ebp + nb304_faction]	
	mov   eax, [ebp + nb304_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb304_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb304_ninner]
	mov   [esp + nb304_ninner], ecx
	add   edx, 0
	mov   [esp + nb304_innerk], edx    ;# number of innerloop atoms 
	jge   .nb304_unroll_loop
	jmp   .nb304_checksingle
.nb304_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb304_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb304_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb304_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb304_jxH1], xmm2
	movapd 	[esp + nb304_jyH1], xmm3
	movapd 	[esp + nb304_jzH1], xmm4
	movapd 	[esp + nb304_jxH2], xmm5
	movapd 	[esp + nb304_jyH2], xmm6
	movapd 	[esp + nb304_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movhpd xmm2, [esi + ebx*8 + 72]
	movhpd xmm3, [esi + ebx*8 + 80]
	movhpd xmm4, [esi + ebx*8 + 88]
	movapd 	[esp + nb304_jxM], xmm2
	movapd 	[esp + nb304_jyM], xmm3
	movapd 	[esp + nb304_jzM], xmm4
	
	movapd xmm0, [esp + nb304_ixM]
	movapd xmm1, [esp + nb304_iyM]
	movapd xmm2, [esp + nb304_izM]
	movapd xmm3, [esp + nb304_ixM]
	movapd xmm4, [esp + nb304_iyM]
	movapd xmm5, [esp + nb304_izM]
	subpd  xmm0, [esp + nb304_jxM]
	subpd  xmm1, [esp + nb304_jyM]
	subpd  xmm2, [esp + nb304_jzM]
	subpd  xmm3, [esp + nb304_jxH1]
	subpd  xmm4, [esp + nb304_jyH1]
	subpd  xmm5, [esp + nb304_jzH1]
	movapd [esp + nb304_dxMM], xmm0
	movapd [esp + nb304_dyMM], xmm1
	movapd [esp + nb304_dzMM], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb304_dxMH1], xmm3
	movapd [esp + nb304_dyMH1], xmm4
	movapd [esp + nb304_dzMH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb304_rsqMM], xmm0
	movapd [esp + nb304_rsqMH1], xmm3

	movapd xmm0, [esp + nb304_ixM]
	movapd xmm1, [esp + nb304_iyM]
	movapd xmm2, [esp + nb304_izM]
	movapd xmm3, [esp + nb304_ixH1]
	movapd xmm4, [esp + nb304_iyH1]
	movapd xmm5, [esp + nb304_izH1]
	subpd  xmm0, [esp + nb304_jxH2]
	subpd  xmm1, [esp + nb304_jyH2]
	subpd  xmm2, [esp + nb304_jzH2]
	subpd  xmm3, [esp + nb304_jxM]
	subpd  xmm4, [esp + nb304_jyM]
	subpd  xmm5, [esp + nb304_jzM]
	movapd [esp + nb304_dxMH2], xmm0
	movapd [esp + nb304_dyMH2], xmm1
	movapd [esp + nb304_dzMH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb304_dxH1M], xmm3
	movapd [esp + nb304_dyH1M], xmm4
	movapd [esp + nb304_dzH1M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb304_rsqMH2], xmm0
	movapd [esp + nb304_rsqH1M], xmm3

	movapd xmm0, [esp + nb304_ixH1]
	movapd xmm1, [esp + nb304_iyH1]
	movapd xmm2, [esp + nb304_izH1]
	movapd xmm3, [esp + nb304_ixH1]
	movapd xmm4, [esp + nb304_iyH1]
	movapd xmm5, [esp + nb304_izH1]
	subpd  xmm0, [esp + nb304_jxH1]
	subpd  xmm1, [esp + nb304_jyH1]
	subpd  xmm2, [esp + nb304_jzH1]
	subpd  xmm3, [esp + nb304_jxH2]
	subpd  xmm4, [esp + nb304_jyH2]
	subpd  xmm5, [esp + nb304_jzH2]
	movapd [esp + nb304_dxH1H1], xmm0
	movapd [esp + nb304_dyH1H1], xmm1
	movapd [esp + nb304_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb304_dxH1H2], xmm3
	movapd [esp + nb304_dyH1H2], xmm4
	movapd [esp + nb304_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb304_rsqH1H1], xmm0
	movapd [esp + nb304_rsqH1H2], xmm3

	movapd xmm0, [esp + nb304_ixH2]
	movapd xmm1, [esp + nb304_iyH2]
	movapd xmm2, [esp + nb304_izH2]
	movapd xmm3, [esp + nb304_ixH2]
	movapd xmm4, [esp + nb304_iyH2]
	movapd xmm5, [esp + nb304_izH2]
	subpd  xmm0, [esp + nb304_jxM]
	subpd  xmm1, [esp + nb304_jyM]
	subpd  xmm2, [esp + nb304_jzM]
	subpd  xmm3, [esp + nb304_jxH1]
	subpd  xmm4, [esp + nb304_jyH1]
	subpd  xmm5, [esp + nb304_jzH1]
	movapd [esp + nb304_dxH2M], xmm0
	movapd [esp + nb304_dyH2M], xmm1
	movapd [esp + nb304_dzH2M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb304_dxH2H1], xmm3
	movapd [esp + nb304_dyH2H1], xmm4
	movapd [esp + nb304_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb304_rsqH2M], xmm0
	movapd [esp + nb304_rsqH2H1], xmm4

	movapd xmm0, [esp + nb304_ixH2]
	movapd xmm1, [esp + nb304_iyH2]
	movapd xmm2, [esp + nb304_izH2]
	subpd  xmm0, [esp + nb304_jxH2]
	subpd  xmm1, [esp + nb304_jyH2]
	subpd  xmm2, [esp + nb304_jzH2]
	movapd [esp + nb304_dxH2H2], xmm0
	movapd [esp + nb304_dyH2H2], xmm1
	movapd [esp + nb304_dzH2H2], xmm2
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb304_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb304_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304_half] ;# iter1 
	mulpd   xmm7, [esp + nb304_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304_half] ;# rinv 
	mulpd   xmm5, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvH2H2], xmm1
	movapd [esp + nb304_rinvH2H1], xmm5

	movapd xmm0, [esp + nb304_rsqMM]
	movapd xmm4, [esp + nb304_rsqMH1]	
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
	movapd  xmm3, [esp + nb304_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb304_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304_half] ;# rinv 
	mulpd   xmm5, [esp + nb304_half] ;# rinv
	movapd [esp + nb304_rinvMM], xmm1
	movapd [esp + nb304_rinvMH1], xmm5

	movapd xmm0, [esp + nb304_rsqMH2]
	movapd xmm4, [esp + nb304_rsqH1M]	
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
	movapd  xmm3, [esp + nb304_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304_half] ;# iter1 
	mulpd   xmm7, [esp + nb304_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304_half] ;# rinv 
	mulpd   xmm5, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvMH2], xmm1
	movapd [esp + nb304_rinvH1M], xmm5

	movapd xmm0, [esp + nb304_rsqH1H1]
	movapd xmm4, [esp + nb304_rsqH1H2]	
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
	movapd  xmm3, [esp + nb304_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304_half] ;# iter1a 
	mulpd   xmm7, [esp + nb304_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304_half] ;# rinv 
	mulpd   xmm5, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvH1H1], xmm1
	movapd [esp + nb304_rinvH1H2], xmm5

	movapd xmm0, [esp + nb304_rsqH2M]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb304_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb304_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb304_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvH2M], xmm1
	
	;# start with MM interaction 
	movapd xmm0, [esp + nb304_rinvMM]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqMM] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	movd mm0, eax	
	movd mm1, ebx
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMM]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addpd  xmm5, [esp + nb304_vctot]
	xorpd  xmm2, xmm2
    movapd [esp + nb304_vctot], xmm5
	mulpd  xmm3, [esp + nb304_tsc]
	
	subpd  xmm2, xmm3
	mulpd  xmm0, xmm2	;# mult by rinv 
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0		

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb304_dxMM]
	mulpd xmm1, [esp + nb304_dyMM]
	mulpd xmm2, [esp + nb304_dzMM]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixM]
	addpd xmm1, [esp + nb304_fiyM]
	addpd xmm2, [esp + nb304_fizM]
	movapd [esp + nb304_fjxM], xmm3
	movapd [esp + nb304_fjyM], xmm4
	movapd [esp + nb304_fjzM], xmm5
	movapd [esp + nb304_fixM], xmm0
	movapd [esp + nb304_fiyM], xmm1
	movapd [esp + nb304_fizM], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb304_rinvMH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqMH1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb304_dxMH1]
	mulpd xmm1, [esp + nb304_dyMH1]
	mulpd xmm2, [esp + nb304_dzMH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixM]
	addpd xmm1, [esp + nb304_fiyM]
	addpd xmm2, [esp + nb304_fizM]
	movapd [esp + nb304_fjxH1], xmm3
	movapd [esp + nb304_fjyH1], xmm4
	movapd [esp + nb304_fjzH1], xmm5
	movapd [esp + nb304_fixM], xmm0
	movapd [esp + nb304_fiyM], xmm1
	movapd [esp + nb304_fizM], xmm2

	;# M-H2 interaction  
	movapd xmm0, [esp + nb304_rinvMH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqMH2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb304_dxMH2]
	mulpd xmm1, [esp + nb304_dyMH2]
	mulpd xmm2, [esp + nb304_dzMH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixM]
	addpd xmm1, [esp + nb304_fiyM]
	addpd xmm2, [esp + nb304_fizM]
	movapd [esp + nb304_fjxH2], xmm3
	movapd [esp + nb304_fjyH2], xmm4
	movapd [esp + nb304_fjzH2], xmm5
	movapd [esp + nb304_fixM], xmm0
	movapd [esp + nb304_fiyM], xmm1
	movapd [esp + nb304_fizM], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb304_rinvH1M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqH1M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxM]
	movapd xmm4, [esp + nb304_fjyM]
	movapd xmm5, [esp + nb304_fjzM]
	mulpd xmm0, [esp + nb304_dxH1M]
	mulpd xmm1, [esp + nb304_dyH1M]
	mulpd xmm2, [esp + nb304_dzH1M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixH1]
	addpd xmm1, [esp + nb304_fiyH1]
	addpd xmm2, [esp + nb304_fizH1]
	movapd [esp + nb304_fjxM], xmm3
	movapd [esp + nb304_fjyM], xmm4
	movapd [esp + nb304_fjzM], xmm5
	movapd [esp + nb304_fixH1], xmm0
	movapd [esp + nb304_fiyH1], xmm1
	movapd [esp + nb304_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb304_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH1]
	movapd xmm4, [esp + nb304_fjyH1]
	movapd xmm5, [esp + nb304_fjzH1]
	mulpd xmm0, [esp + nb304_dxH1H1]
	mulpd xmm1, [esp + nb304_dyH1H1]
	mulpd xmm2, [esp + nb304_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixH1]
	addpd xmm1, [esp + nb304_fiyH1]
	addpd xmm2, [esp + nb304_fizH1]
	movapd [esp + nb304_fjxH1], xmm3
	movapd [esp + nb304_fjyH1], xmm4
	movapd [esp + nb304_fjzH1], xmm5
	movapd [esp + nb304_fixH1], xmm0
	movapd [esp + nb304_fiyH1], xmm1
	movapd [esp + nb304_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb304_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH2]
	movapd xmm4, [esp + nb304_fjyH2]
	movapd xmm5, [esp + nb304_fjzH2]
	mulpd xmm0, [esp + nb304_dxH1H2]
	mulpd xmm1, [esp + nb304_dyH1H2]
	mulpd xmm2, [esp + nb304_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixH1]
	addpd xmm1, [esp + nb304_fiyH1]
	addpd xmm2, [esp + nb304_fizH1]
	movapd [esp + nb304_fjxH2], xmm3
	movapd [esp + nb304_fjyH2], xmm4
	movapd [esp + nb304_fjzH2], xmm5
	movapd [esp + nb304_fixH1], xmm0
	movapd [esp + nb304_fiyH1], xmm1
	movapd [esp + nb304_fizH1], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb304_rinvH2M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqH2M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1

	movapd xmm3, [esp + nb304_fjxM]
	movapd xmm4, [esp + nb304_fjyM]
	movapd xmm5, [esp + nb304_fjzM]
	mulpd xmm0, [esp + nb304_dxH2M]
	mulpd xmm1, [esp + nb304_dyH2M]
	mulpd xmm2, [esp + nb304_dzH2M]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixH2]
	addpd xmm1, [esp + nb304_fiyH2]
	addpd xmm2, [esp + nb304_fizH2]
	movapd [esp + nb304_fjxM], xmm3
	movapd [esp + nb304_fjyM], xmm4
	movapd [esp + nb304_fjzM], xmm5
	movapd [esp + nb304_fixH2], xmm0
	movapd [esp + nb304_fiyH2], xmm1
	movapd [esp + nb304_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb304_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH1]
	movapd xmm4, [esp + nb304_fjyH1]
	movapd xmm5, [esp + nb304_fjzH1]
	mulpd xmm0, [esp + nb304_dxH2H1]
	mulpd xmm1, [esp + nb304_dyH2H1]
	mulpd xmm2, [esp + nb304_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixH2]
	addpd xmm1, [esp + nb304_fiyH2]
	addpd xmm2, [esp + nb304_fizH2]
	movapd [esp + nb304_fjxH1], xmm3
	movapd [esp + nb304_fjyH1], xmm4
	movapd [esp + nb304_fjzH1], xmm5
	movapd [esp + nb304_fixH2], xmm0
	movapd [esp + nb304_fiyH2], xmm1
	movapd [esp + nb304_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb304_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
	
    addpd  xmm5, [esp + nb304_vctot]
    movapd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb304_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH2]
	movapd xmm4, [esp + nb304_fjyH2]
	movapd xmm5, [esp + nb304_fjzH2]
	mulpd xmm0, [esp + nb304_dxH2H2]
	mulpd xmm1, [esp + nb304_dyH2H2]
	mulpd xmm2, [esp + nb304_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb304_fixH2]
	addpd xmm1, [esp + nb304_fiyH2]
	addpd xmm2, [esp + nb304_fizH2]
	movapd [esp + nb304_fjxH2], xmm3
	movapd [esp + nb304_fjyH2], xmm4
	movapd [esp + nb304_fjzH2], xmm5
	movapd [esp + nb304_fixH2], xmm0
	movapd [esp + nb304_fiyH2], xmm1
	movapd [esp + nb304_fizH2], xmm2

	mov edi, [ebp + nb304_faction]

	movd eax, mm0
	movd ebx, mm1
	
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
	addpd xmm0, [esp + nb304_fjxH1]
	addpd xmm1, [esp + nb304_fjyH1]
	addpd xmm2, [esp + nb304_fjzH1]
	addpd xmm3, [esp + nb304_fjxH2]
	addpd xmm4, [esp + nb304_fjyH2]
	addpd xmm5, [esp + nb304_fjzH2]
	addpd xmm6, [esp + nb304_fjxM]
	addpd xmm7, [esp + nb304_fjyM]
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
	addpd xmm0, [esp + nb304_fjzM]
	movlpd [edi + eax*8 + 88], xmm0
	movhpd [edi + ebx*8 + 88], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb304_innerk],  2
	jl    .nb304_checksingle
	jmp   .nb304_unroll_loop
.nb304_checksingle:
	mov   edx, [esp + nb304_innerk]
	and   edx, 1
	jnz   .nb304_dosingle
	jmp   .nb304_updateouterdata
.nb304_dosingle:
	mov   edx, [esp + nb304_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb304_pos]
	lea   eax, [eax + eax*2]  

	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movapd 	[esp + nb304_jxH1], xmm2
	movapd 	[esp + nb304_jyH1], xmm3
	movapd 	[esp + nb304_jzH1], xmm4
	movapd 	[esp + nb304_jxH2], xmm5
	movapd 	[esp + nb304_jyH2], xmm6
	movapd 	[esp + nb304_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movapd 	[esp + nb304_jxM], xmm2
	movapd 	[esp + nb304_jyM], xmm3
	movapd 	[esp + nb304_jzM], xmm4
	
	movapd xmm0, [esp + nb304_ixM]
	movapd xmm1, [esp + nb304_iyM]
	movapd xmm2, [esp + nb304_izM]
	movapd xmm3, [esp + nb304_ixM]
	movapd xmm4, [esp + nb304_iyM]
	movapd xmm5, [esp + nb304_izM]
	subsd  xmm0, [esp + nb304_jxM]
	subsd  xmm1, [esp + nb304_jyM]
	subsd  xmm2, [esp + nb304_jzM]
	subsd  xmm3, [esp + nb304_jxH1]
	subsd  xmm4, [esp + nb304_jyH1]
	subsd  xmm5, [esp + nb304_jzH1]
	movapd [esp + nb304_dxMM], xmm0
	movapd [esp + nb304_dyMM], xmm1
	movapd [esp + nb304_dzMM], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb304_dxMH1], xmm3
	movapd [esp + nb304_dyMH1], xmm4
	movapd [esp + nb304_dzMH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb304_rsqMM], xmm0
	movapd [esp + nb304_rsqMH1], xmm3

	movapd xmm0, [esp + nb304_ixM]
	movapd xmm1, [esp + nb304_iyM]
	movapd xmm2, [esp + nb304_izM]
	movapd xmm3, [esp + nb304_ixH1]
	movapd xmm4, [esp + nb304_iyH1]
	movapd xmm5, [esp + nb304_izH1]
	subsd  xmm0, [esp + nb304_jxH2]
	subsd  xmm1, [esp + nb304_jyH2]
	subsd  xmm2, [esp + nb304_jzH2]
	subsd  xmm3, [esp + nb304_jxM]
	subsd  xmm4, [esp + nb304_jyM]
	subsd  xmm5, [esp + nb304_jzM]
	movapd [esp + nb304_dxMH2], xmm0
	movapd [esp + nb304_dyMH2], xmm1
	movapd [esp + nb304_dzMH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb304_dxH1M], xmm3
	movapd [esp + nb304_dyH1M], xmm4
	movapd [esp + nb304_dzH1M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb304_rsqMH2], xmm0
	movapd [esp + nb304_rsqH1M], xmm3

	movapd xmm0, [esp + nb304_ixH1]
	movapd xmm1, [esp + nb304_iyH1]
	movapd xmm2, [esp + nb304_izH1]
	movapd xmm3, [esp + nb304_ixH1]
	movapd xmm4, [esp + nb304_iyH1]
	movapd xmm5, [esp + nb304_izH1]
	subsd  xmm0, [esp + nb304_jxH1]
	subsd  xmm1, [esp + nb304_jyH1]
	subsd  xmm2, [esp + nb304_jzH1]
	subsd  xmm3, [esp + nb304_jxH2]
	subsd  xmm4, [esp + nb304_jyH2]
	subsd  xmm5, [esp + nb304_jzH2]
	movapd [esp + nb304_dxH1H1], xmm0
	movapd [esp + nb304_dyH1H1], xmm1
	movapd [esp + nb304_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb304_dxH1H2], xmm3
	movapd [esp + nb304_dyH1H2], xmm4
	movapd [esp + nb304_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb304_rsqH1H1], xmm0
	movapd [esp + nb304_rsqH1H2], xmm3

	movapd xmm0, [esp + nb304_ixH2]
	movapd xmm1, [esp + nb304_iyH2]
	movapd xmm2, [esp + nb304_izH2]
	movapd xmm3, [esp + nb304_ixH2]
	movapd xmm4, [esp + nb304_iyH2]
	movapd xmm5, [esp + nb304_izH2]
	subsd  xmm0, [esp + nb304_jxM]
	subsd  xmm1, [esp + nb304_jyM]
	subsd  xmm2, [esp + nb304_jzM]
	subsd  xmm3, [esp + nb304_jxH1]
	subsd  xmm4, [esp + nb304_jyH1]
	subsd  xmm5, [esp + nb304_jzH1]
	movapd [esp + nb304_dxH2M], xmm0
	movapd [esp + nb304_dyH2M], xmm1
	movapd [esp + nb304_dzH2M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb304_dxH2H1], xmm3
	movapd [esp + nb304_dyH2H1], xmm4
	movapd [esp + nb304_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb304_rsqH2M], xmm0
	movapd [esp + nb304_rsqH2H1], xmm4

	movapd xmm0, [esp + nb304_ixH2]
	movapd xmm1, [esp + nb304_iyH2]
	movapd xmm2, [esp + nb304_izH2]
	subsd  xmm0, [esp + nb304_jxH2]
	subsd  xmm1, [esp + nb304_jyH2]
	subsd  xmm2, [esp + nb304_jzH2]
	movapd [esp + nb304_dxH2H2], xmm0
	movapd [esp + nb304_dyH2H2], xmm1
	movapd [esp + nb304_dzH2H2], xmm2
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb304_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb304_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304_half] ;# iter1 
	mulsd   xmm7, [esp + nb304_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304_half] ;# rinv 
	mulsd   xmm5, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvH2H2], xmm1
	movapd [esp + nb304_rinvH2H1], xmm5

	movapd xmm0, [esp + nb304_rsqMM]
	movapd xmm4, [esp + nb304_rsqMH1]	
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
	movapd  xmm3, [esp + nb304_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb304_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304_half] ;# rinv 
	mulsd   xmm5, [esp + nb304_half] ;# rinv
	movapd [esp + nb304_rinvMM], xmm1
	movapd [esp + nb304_rinvMH1], xmm5

	movapd xmm0, [esp + nb304_rsqMH2]
	movapd xmm4, [esp + nb304_rsqH1M]	
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
	movapd  xmm3, [esp + nb304_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304_half] ;# iter1 
	mulsd   xmm7, [esp + nb304_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304_half] ;# rinv 
	mulsd   xmm5, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvMH2], xmm1
	movapd [esp + nb304_rinvH1M], xmm5

	movapd xmm0, [esp + nb304_rsqH1H1]
	movapd xmm4, [esp + nb304_rsqH1H2]	
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
	movapd  xmm3, [esp + nb304_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304_half] ;# iter1a 
	mulsd   xmm7, [esp + nb304_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304_half] ;# rinv 
	mulsd   xmm5, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvH1H1], xmm1
	movapd [esp + nb304_rinvH1H2], xmm5

	movapd xmm0, [esp + nb304_rsqH2M]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb304_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb304_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb304_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb304_half] ;# rinv 
	movapd [esp + nb304_rinvH2M], xmm1
	
	movd mm0, eax	
	;# start with MM interaction 
	movapd xmm0, [esp + nb304_rinvMM]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqMM] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMM]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addsd  xmm5, [esp + nb304_vctot]
	xorpd  xmm2, xmm2
    movlpd [esp + nb304_vctot], xmm5
	mulsd  xmm3, [esp + nb304_tsc]
	
	subsd  xmm2, xmm3
	mulsd  xmm0, xmm2
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0		

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb304_dxMM]
	mulsd xmm1, [esp + nb304_dyMM]
	mulsd xmm2, [esp + nb304_dzMM]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixM]
	addsd xmm1, [esp + nb304_fiyM]
	addsd xmm2, [esp + nb304_fizM]
	movlpd [esp + nb304_fjxM], xmm3
	movlpd [esp + nb304_fjyM], xmm4
	movlpd [esp + nb304_fjzM], xmm5
	movlpd [esp + nb304_fixM], xmm0
	movlpd [esp + nb304_fiyM], xmm1
	movlpd [esp + nb304_fizM], xmm2

	;# M-H1 interaction 
	movapd xmm0, [esp + nb304_rinvMH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqMH1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb304_dxMH1]
	mulsd xmm1, [esp + nb304_dyMH1]
	mulsd xmm2, [esp + nb304_dzMH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixM]
	addsd xmm1, [esp + nb304_fiyM]
	addsd xmm2, [esp + nb304_fizM]
	movlpd [esp + nb304_fjxH1], xmm3
	movlpd [esp + nb304_fjyH1], xmm4
	movlpd [esp + nb304_fjzH1], xmm5
	movlpd [esp + nb304_fixM], xmm0
	movlpd [esp + nb304_fiyM], xmm1
	movlpd [esp + nb304_fizM], xmm2

	;# M-H2 interaction  
	movapd xmm0, [esp + nb304_rinvMH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqMH2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb304_dxMH2]
	mulsd xmm1, [esp + nb304_dyMH2]
	mulsd xmm2, [esp + nb304_dzMH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixM]
	addsd xmm1, [esp + nb304_fiyM]
	addsd xmm2, [esp + nb304_fizM]
	movlpd [esp + nb304_fjxH2], xmm3
	movlpd [esp + nb304_fjyH2], xmm4
	movlpd [esp + nb304_fjzH2], xmm5
	movlpd [esp + nb304_fixM], xmm0
	movlpd [esp + nb304_fiyM], xmm1
	movlpd [esp + nb304_fizM], xmm2

	;# H1-M interaction 
	movapd xmm0, [esp + nb304_rinvH1M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqH1M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxM]
	movapd xmm4, [esp + nb304_fjyM]
	movapd xmm5, [esp + nb304_fjzM]
	mulsd xmm0, [esp + nb304_dxH1M]
	mulsd xmm1, [esp + nb304_dyH1M]
	mulsd xmm2, [esp + nb304_dzH1M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixH1]
	addsd xmm1, [esp + nb304_fiyH1]
	addsd xmm2, [esp + nb304_fizH1]
	movlpd [esp + nb304_fjxM], xmm3
	movlpd [esp + nb304_fjyM], xmm4
	movlpd [esp + nb304_fjzM], xmm5
	movlpd [esp + nb304_fixH1], xmm0
	movlpd [esp + nb304_fiyH1], xmm1
	movlpd [esp + nb304_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb304_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH1]
	movapd xmm4, [esp + nb304_fjyH1]
	movapd xmm5, [esp + nb304_fjzH1]
	mulsd xmm0, [esp + nb304_dxH1H1]
	mulsd xmm1, [esp + nb304_dyH1H1]
	mulsd xmm2, [esp + nb304_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixH1]
	addsd xmm1, [esp + nb304_fiyH1]
	addsd xmm2, [esp + nb304_fizH1]
	movlpd [esp + nb304_fjxH1], xmm3
	movlpd [esp + nb304_fjyH1], xmm4
	movlpd [esp + nb304_fjzH1], xmm5
	movlpd [esp + nb304_fixH1], xmm0
	movlpd [esp + nb304_fiyH1], xmm1
	movlpd [esp + nb304_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb304_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH2]
	movapd xmm4, [esp + nb304_fjyH2]
	movapd xmm5, [esp + nb304_fjzH2]
	mulsd xmm0, [esp + nb304_dxH1H2]
	mulsd xmm1, [esp + nb304_dyH1H2]
	mulsd xmm2, [esp + nb304_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixH1]
	addsd xmm1, [esp + nb304_fiyH1]
	addsd xmm2, [esp + nb304_fizH1]
	movlpd [esp + nb304_fjxH2], xmm3
	movlpd [esp + nb304_fjyH2], xmm4
	movlpd [esp + nb304_fjzH2], xmm5
	movlpd [esp + nb304_fixH1], xmm0
	movlpd [esp + nb304_fiyH1], xmm1
	movlpd [esp + nb304_fizH1], xmm2

	;# H2-M interaction 
	movapd xmm0, [esp + nb304_rinvH2M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqH2M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqMH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1

	movapd xmm3, [esp + nb304_fjxM]
	movapd xmm4, [esp + nb304_fjyM]
	movapd xmm5, [esp + nb304_fjzM]
	mulsd xmm0, [esp + nb304_dxH2M]
	mulsd xmm1, [esp + nb304_dyH2M]
	mulsd xmm2, [esp + nb304_dzH2M]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixH2]
	addsd xmm1, [esp + nb304_fiyH2]
	addsd xmm2, [esp + nb304_fizH2]
	movlpd [esp + nb304_fjxM], xmm3
	movlpd [esp + nb304_fjyM], xmm4
	movlpd [esp + nb304_fjzM], xmm5
	movlpd [esp + nb304_fixH2], xmm0
	movlpd [esp + nb304_fiyH2], xmm1
	movlpd [esp + nb304_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb304_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH1]
	movapd xmm4, [esp + nb304_fjyH1]
	movapd xmm5, [esp + nb304_fjzH1]
	mulsd xmm0, [esp + nb304_dxH2H1]
	mulsd xmm1, [esp + nb304_dyH2H1]
	mulsd xmm2, [esp + nb304_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixH2]
	addsd xmm1, [esp + nb304_fiyH2]
	addsd xmm2, [esp + nb304_fizH2]
	movlpd [esp + nb304_fjxH1], xmm3
	movlpd [esp + nb304_fjyH1], xmm4
	movlpd [esp + nb304_fjzH1], xmm5
	movlpd [esp + nb304_fixH2], xmm0
	movlpd [esp + nb304_fiyH2], xmm1
	movlpd [esp + nb304_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb304_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb304_two]	;# two*Heps2 
	movapd xmm3, [esp + nb304_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb304_vctot]
    movlpd [esp + nb304_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb304_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb304_fjxH2]
	movapd xmm4, [esp + nb304_fjyH2]
	movapd xmm5, [esp + nb304_fjzH2]
	mulsd xmm0, [esp + nb304_dxH2H2]
	mulsd xmm1, [esp + nb304_dyH2H2]
	mulsd xmm2, [esp + nb304_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb304_fixH2]
	addsd xmm1, [esp + nb304_fiyH2]
	addsd xmm2, [esp + nb304_fizH2]
	movlpd [esp + nb304_fjxH2], xmm3
	movlpd [esp + nb304_fjyH2], xmm4
	movlpd [esp + nb304_fjzH2], xmm5
	movlpd [esp + nb304_fixH2], xmm0
	movlpd [esp + nb304_fiyH2], xmm1
	movlpd [esp + nb304_fizH2], xmm2

	mov edi, [ebp + nb304_faction]

	movd eax, mm0
	
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8 + 24]
	movlpd xmm1, [edi + eax*8 + 32]
	movlpd xmm2, [edi + eax*8 + 40]
	movlpd xmm3, [edi + eax*8 + 48]
	movlpd xmm4, [edi + eax*8 + 56]
	movlpd xmm5, [edi + eax*8 + 64]
	movlpd xmm6, [edi + eax*8 + 72]
	movlpd xmm7, [edi + eax*8 + 80]
	addsd xmm0, [esp + nb304_fjxH1]
	addsd xmm1, [esp + nb304_fjyH1]
	addsd xmm2, [esp + nb304_fjzH1]
	addsd xmm3, [esp + nb304_fjxH2]
	addsd xmm4, [esp + nb304_fjyH2]
	addsd xmm5, [esp + nb304_fjzH2]
	addsd xmm6, [esp + nb304_fjxM]
	addsd xmm7, [esp + nb304_fjyM]
	movlpd [edi + eax*8 + 24], xmm0
	movlpd [edi + eax*8 + 32], xmm1
	movlpd [edi + eax*8 + 40], xmm2
	movlpd [edi + eax*8 + 48], xmm3
	movlpd [edi + eax*8 + 56], xmm4
	movlpd [edi + eax*8 + 64], xmm5
	movlpd [edi + eax*8 + 72], xmm6
	movlpd [edi + eax*8 + 80], xmm7

	movlpd xmm0, [edi + eax*8 + 88]
	addsd xmm0, [esp + nb304_fjzM]
	movlpd [edi + eax*8 + 88], xmm0
	
.nb304_updateouterdata:
	mov   ecx, [esp + nb304_ii3]
	mov   edi, [ebp + nb304_faction]
	mov   esi, [ebp + nb304_fshift]
	mov   edx, [esp + nb304_is3]

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb304_fixH1]
	movapd xmm1, [esp + nb304_fiyH1]
	movapd xmm2, [esp + nb304_fizH1]

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
	movapd xmm0, [esp + nb304_fixH2]
	movapd xmm1, [esp + nb304_fiyH2]
	movapd xmm2, [esp + nb304_fizH2]

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
	movapd xmm0, [esp + nb304_fixM]
	movapd xmm1, [esp + nb304_fiyM]
	movapd xmm2, [esp + nb304_fizM]

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
	mov esi, [esp + nb304_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb304_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb304_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb304_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb304_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb304_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb304_n], esi
        jmp .nb304_outer
.nb304_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb304_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb304_end
        ;# non-zero, do one more workunit
        jmp   .nb304_threadloop
.nb304_end:
	emms

	mov eax, [esp + nb304_nouter]
	mov ebx, [esp + nb304_ninner]
	mov ecx, [ebp + nb304_outeriter]
	mov edx, [ebp + nb304_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb304_salign]
	add esp, eax
	add esp, 1464
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl nb_kernel304nf_ia32_sse2
.globl _nb_kernel304nf_ia32_sse2
nb_kernel304nf_ia32_sse2:	
_nb_kernel304nf_ia32_sse2:	
.equiv          nb304nf_p_nri,          8
.equiv          nb304nf_iinr,           12
.equiv          nb304nf_jindex,         16
.equiv          nb304nf_jjnr,           20
.equiv          nb304nf_shift,          24
.equiv          nb304nf_shiftvec,       28
.equiv          nb304nf_fshift,         32
.equiv          nb304nf_gid,            36
.equiv          nb304nf_pos,            40
.equiv          nb304nf_faction,        44
.equiv          nb304nf_charge,         48
.equiv          nb304nf_p_facel,        52
.equiv          nb304nf_argkrf,         56
.equiv          nb304nf_argcrf,         60
.equiv          nb304nf_Vc,             64
.equiv          nb304nf_type,           68
.equiv          nb304nf_p_ntype,        72
.equiv          nb304nf_vdwparam,       76
.equiv          nb304nf_Vvdw,           80
.equiv          nb304nf_p_tabscale,     84
.equiv          nb304nf_VFtab,          88
.equiv          nb304nf_invsqrta,       92
.equiv          nb304nf_dvda,           96
.equiv          nb304nf_p_gbtabscale,   100
.equiv          nb304nf_GBtab,          104
.equiv          nb304nf_p_nthreads,     108
.equiv          nb304nf_count,          112
.equiv          nb304nf_mtx,            116
.equiv          nb304nf_outeriter,      120
.equiv          nb304nf_inneriter,      124
.equiv          nb304nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb304nf_ixM,            0
.equiv          nb304nf_iyM,            16
.equiv          nb304nf_izM,            32
.equiv          nb304nf_ixH1,           48
.equiv          nb304nf_iyH1,           64
.equiv          nb304nf_izH1,           80
.equiv          nb304nf_ixH2,           96
.equiv          nb304nf_iyH2,           112
.equiv          nb304nf_izH2,           128
.equiv          nb304nf_jxM,            144
.equiv          nb304nf_jyM,            160
.equiv          nb304nf_jzM,            176
.equiv          nb304nf_jxH1,           192
.equiv          nb304nf_jyH1,           208
.equiv          nb304nf_jzH1,           224
.equiv          nb304nf_jxH2,           240
.equiv          nb304nf_jyH2,           256
.equiv          nb304nf_jzH2,           272
.equiv          nb304nf_qqMM,           288
.equiv          nb304nf_qqMH,           304
.equiv          nb304nf_qqHH,           320
.equiv          nb304nf_tsc,            336
.equiv          nb304nf_vctot,          352
.equiv          nb304nf_half,           368
.equiv          nb304nf_three,          384
.equiv          nb304nf_rsqMM,          400
.equiv          nb304nf_rsqMH1,         416
.equiv          nb304nf_rsqMH2,         432
.equiv          nb304nf_rsqH1M,         448
.equiv          nb304nf_rsqH1H1,        464
.equiv          nb304nf_rsqH1H2,        480
.equiv          nb304nf_rsqH2M,         496
.equiv          nb304nf_rsqH2H1,        512
.equiv          nb304nf_rsqH2H2,        528
.equiv          nb304nf_rinvMM,         544
.equiv          nb304nf_rinvMH1,        560
.equiv          nb304nf_rinvMH2,        576
.equiv          nb304nf_rinvH1M,        592
.equiv          nb304nf_rinvH1H1,       608
.equiv          nb304nf_rinvH1H2,       624
.equiv          nb304nf_rinvH2M,        640
.equiv          nb304nf_rinvH2H1,       656
.equiv          nb304nf_rinvH2H2,       672
.equiv          nb304nf_is3,            688
.equiv          nb304nf_ii3,            692
.equiv          nb304nf_innerjjnr,      696
.equiv          nb304nf_innerk,         700
.equiv          nb304nf_n,              704
.equiv          nb304nf_nn1,            708
.equiv          nb304nf_nri,            712
.equiv          nb304nf_nouter,         716
.equiv          nb304nf_ninner,         720
.equiv          nb304nf_salign,         724
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 728		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb304nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb304nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb304nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb304nf_nouter], eax
	mov [esp + nb304nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb304nf_half], eax
	mov [esp + nb304nf_half+4], ebx
	movsd xmm1, [esp + nb304nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb304nf_half], xmm1
	movapd [esp + nb304nf_three], xmm3

	mov eax, [ebp + nb304nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb304nf_tsc],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb304nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb304nf_charge]
	movsd xmm3, [edx + ebx*8 + 24]	
	movsd xmm4, xmm3
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb304nf_p_facel]
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
	movapd [esp + nb304nf_qqMM], xmm3
	movapd [esp + nb304nf_qqMH], xmm4
	movapd [esp + nb304nf_qqHH], xmm5		

.nb304nf_threadloop:
        mov   esi, [ebp + nb304nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb304nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb304nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb304nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb304nf_n], eax
        mov [esp + nb304nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb304nf_outerstart
        jmp .nb304nf_end

.nb304nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb304nf_nouter]
	mov [esp + nb304nf_nouter], ebx

.nb304nf_outer:
	mov   eax, [ebp + nb304nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb304nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb304nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb304nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb304nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb304nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb304nf_ixH1], xmm3
	movapd [esp + nb304nf_iyH1], xmm4
	movapd [esp + nb304nf_izH1], xmm5

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
	movapd [esp + nb304nf_ixH2], xmm0
	movapd [esp + nb304nf_iyH2], xmm1
	movapd [esp + nb304nf_izH2], xmm2
	movapd [esp + nb304nf_ixM], xmm3
	movapd [esp + nb304nf_iyM], xmm4
	movapd [esp + nb304nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb304nf_vctot], xmm4
	
	mov   eax, [ebp + nb304nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb304nf_pos]
	mov   eax, [ebp + nb304nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb304nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb304nf_ninner]
	mov   [esp + nb304nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb304nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb304nf_unroll_loop
	jmp   .nb304nf_checksingle
.nb304nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb304nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb304nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb304nf_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb304nf_jxH1], xmm2
	movapd 	[esp + nb304nf_jyH1], xmm3
	movapd 	[esp + nb304nf_jzH1], xmm4
	movapd 	[esp + nb304nf_jxH2], xmm5
	movapd 	[esp + nb304nf_jyH2], xmm6
	movapd 	[esp + nb304nf_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movhpd xmm2, [esi + ebx*8 + 72]
	movhpd xmm3, [esi + ebx*8 + 80]
	movhpd xmm4, [esi + ebx*8 + 88]
	movapd 	[esp + nb304nf_jxM], xmm2
	movapd 	[esp + nb304nf_jyM], xmm3
	movapd 	[esp + nb304nf_jzM], xmm4
	
	movapd xmm0, [esp + nb304nf_ixM]
	movapd xmm1, [esp + nb304nf_iyM]
	movapd xmm2, [esp + nb304nf_izM]
	movapd xmm3, [esp + nb304nf_ixM]
	movapd xmm4, [esp + nb304nf_iyM]
	movapd xmm5, [esp + nb304nf_izM]
	subpd  xmm0, [esp + nb304nf_jxM]
	subpd  xmm1, [esp + nb304nf_jyM]
	subpd  xmm2, [esp + nb304nf_jzM]
	subpd  xmm3, [esp + nb304nf_jxH1]
	subpd  xmm4, [esp + nb304nf_jyH1]
	subpd  xmm5, [esp + nb304nf_jzH1]
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
	movapd [esp + nb304nf_rsqMM], xmm0
	movapd [esp + nb304nf_rsqMH1], xmm3

	movapd xmm0, [esp + nb304nf_ixM]
	movapd xmm1, [esp + nb304nf_iyM]
	movapd xmm2, [esp + nb304nf_izM]
	movapd xmm3, [esp + nb304nf_ixH1]
	movapd xmm4, [esp + nb304nf_iyH1]
	movapd xmm5, [esp + nb304nf_izH1]
	subpd  xmm0, [esp + nb304nf_jxH2]
	subpd  xmm1, [esp + nb304nf_jyH2]
	subpd  xmm2, [esp + nb304nf_jzH2]
	subpd  xmm3, [esp + nb304nf_jxM]
	subpd  xmm4, [esp + nb304nf_jyM]
	subpd  xmm5, [esp + nb304nf_jzM]
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
	movapd [esp + nb304nf_rsqMH2], xmm0
	movapd [esp + nb304nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb304nf_ixH1]
	movapd xmm1, [esp + nb304nf_iyH1]
	movapd xmm2, [esp + nb304nf_izH1]
	movapd xmm3, [esp + nb304nf_ixH1]
	movapd xmm4, [esp + nb304nf_iyH1]
	movapd xmm5, [esp + nb304nf_izH1]
	subpd  xmm0, [esp + nb304nf_jxH1]
	subpd  xmm1, [esp + nb304nf_jyH1]
	subpd  xmm2, [esp + nb304nf_jzH1]
	subpd  xmm3, [esp + nb304nf_jxH2]
	subpd  xmm4, [esp + nb304nf_jyH2]
	subpd  xmm5, [esp + nb304nf_jzH2]
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
	movapd [esp + nb304nf_rsqH1H1], xmm0
	movapd [esp + nb304nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb304nf_ixH2]
	movapd xmm1, [esp + nb304nf_iyH2]
	movapd xmm2, [esp + nb304nf_izH2]
	movapd xmm3, [esp + nb304nf_ixH2]
	movapd xmm4, [esp + nb304nf_iyH2]
	movapd xmm5, [esp + nb304nf_izH2]
	subpd  xmm0, [esp + nb304nf_jxM]
	subpd  xmm1, [esp + nb304nf_jyM]
	subpd  xmm2, [esp + nb304nf_jzM]
	subpd  xmm3, [esp + nb304nf_jxH1]
	subpd  xmm4, [esp + nb304nf_jyH1]
	subpd  xmm5, [esp + nb304nf_jzH1]
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
	movapd [esp + nb304nf_rsqH2M], xmm0
	movapd [esp + nb304nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb304nf_ixH2]
	movapd xmm1, [esp + nb304nf_iyH2]
	movapd xmm2, [esp + nb304nf_izH2]
	subpd  xmm0, [esp + nb304nf_jxH2]
	subpd  xmm1, [esp + nb304nf_jyH2]
	subpd  xmm2, [esp + nb304nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb304nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb304nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb304nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvH2H2], xmm1
	movapd [esp + nb304nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb304nf_rsqMM]
	movapd xmm4, [esp + nb304nf_rsqMH1]	
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
	movapd  xmm3, [esp + nb304nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb304nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb304nf_half] ;# rinv
	movapd [esp + nb304nf_rinvMM], xmm1
	movapd [esp + nb304nf_rinvMH1], xmm5

	movapd xmm0, [esp + nb304nf_rsqMH2]
	movapd xmm4, [esp + nb304nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb304nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb304nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvMH2], xmm1
	movapd [esp + nb304nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb304nf_rsqH1H1]
	movapd xmm4, [esp + nb304nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb304nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb304nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb304nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvH1H1], xmm1
	movapd [esp + nb304nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb304nf_rsqH2M]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb304nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb304nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb304nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvH2M], xmm1
	
	;# start with MM interaction 
	movapd xmm0, [esp + nb304nf_rinvMM]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqMM] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMM]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5

	;# M-H1 interaction 
	movapd xmm0, [esp + nb304nf_rinvMH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqMH1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5

	;# M-H2 interaction  
	movapd xmm0, [esp + nb304nf_rinvMH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqMH2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5

	;# H1-M interaction 
	movapd xmm0, [esp + nb304nf_rinvH1M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqH1M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb304nf_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5
	
	;# H1-H2 interaction 
	movapd xmm0, [esp + nb304nf_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5

	;# H2-M interaction 
	movapd xmm0, [esp + nb304nf_rinvH2M]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqH2M] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb304nf_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5	

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb304nf_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb304nf_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb304nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
	
    addpd  xmm5, [esp + nb304nf_vctot]
    movapd [esp + nb304nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb304nf_innerk],  2
	jl    .nb304nf_checksingle
	jmp   .nb304nf_unroll_loop
.nb304nf_checksingle:
	mov   edx, [esp + nb304nf_innerk]
	and   edx, 1
	jnz   .nb304nf_dosingle
	jmp   .nb304nf_updateouterdata
.nb304nf_dosingle:
	mov   edx, [esp + nb304nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb304nf_pos]
	lea   eax, [eax + eax*2]  

	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8 + 24]
	movlpd xmm3, [esi + eax*8 + 32]
	movlpd xmm4, [esi + eax*8 + 40]
	movlpd xmm5, [esi + eax*8 + 48]
	movlpd xmm6, [esi + eax*8 + 56]
	movlpd xmm7, [esi + eax*8 + 64]
	movapd 	[esp + nb304nf_jxH1], xmm2
	movapd 	[esp + nb304nf_jyH1], xmm3
	movapd 	[esp + nb304nf_jzH1], xmm4
	movapd 	[esp + nb304nf_jxH2], xmm5
	movapd 	[esp + nb304nf_jyH2], xmm6
	movapd 	[esp + nb304nf_jzH2], xmm7
	movlpd xmm2, [esi + eax*8 + 72]
	movlpd xmm3, [esi + eax*8 + 80]
	movlpd xmm4, [esi + eax*8 + 88]
	movapd 	[esp + nb304nf_jxM], xmm2
	movapd 	[esp + nb304nf_jyM], xmm3
	movapd 	[esp + nb304nf_jzM], xmm4
	
	movapd xmm0, [esp + nb304nf_ixM]
	movapd xmm1, [esp + nb304nf_iyM]
	movapd xmm2, [esp + nb304nf_izM]
	movapd xmm3, [esp + nb304nf_ixM]
	movapd xmm4, [esp + nb304nf_iyM]
	movapd xmm5, [esp + nb304nf_izM]
	subsd  xmm0, [esp + nb304nf_jxM]
	subsd  xmm1, [esp + nb304nf_jyM]
	subsd  xmm2, [esp + nb304nf_jzM]
	subsd  xmm3, [esp + nb304nf_jxH1]
	subsd  xmm4, [esp + nb304nf_jyH1]
	subsd  xmm5, [esp + nb304nf_jzH1]
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
	movapd [esp + nb304nf_rsqMM], xmm0
	movapd [esp + nb304nf_rsqMH1], xmm3

	movapd xmm0, [esp + nb304nf_ixM]
	movapd xmm1, [esp + nb304nf_iyM]
	movapd xmm2, [esp + nb304nf_izM]
	movapd xmm3, [esp + nb304nf_ixH1]
	movapd xmm4, [esp + nb304nf_iyH1]
	movapd xmm5, [esp + nb304nf_izH1]
	subsd  xmm0, [esp + nb304nf_jxH2]
	subsd  xmm1, [esp + nb304nf_jyH2]
	subsd  xmm2, [esp + nb304nf_jzH2]
	subsd  xmm3, [esp + nb304nf_jxM]
	subsd  xmm4, [esp + nb304nf_jyM]
	subsd  xmm5, [esp + nb304nf_jzM]
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
	movapd [esp + nb304nf_rsqMH2], xmm0
	movapd [esp + nb304nf_rsqH1M], xmm3

	movapd xmm0, [esp + nb304nf_ixH1]
	movapd xmm1, [esp + nb304nf_iyH1]
	movapd xmm2, [esp + nb304nf_izH1]
	movapd xmm3, [esp + nb304nf_ixH1]
	movapd xmm4, [esp + nb304nf_iyH1]
	movapd xmm5, [esp + nb304nf_izH1]
	subsd  xmm0, [esp + nb304nf_jxH1]
	subsd  xmm1, [esp + nb304nf_jyH1]
	subsd  xmm2, [esp + nb304nf_jzH1]
	subsd  xmm3, [esp + nb304nf_jxH2]
	subsd  xmm4, [esp + nb304nf_jyH2]
	subsd  xmm5, [esp + nb304nf_jzH2]
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
	movapd [esp + nb304nf_rsqH1H1], xmm0
	movapd [esp + nb304nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb304nf_ixH2]
	movapd xmm1, [esp + nb304nf_iyH2]
	movapd xmm2, [esp + nb304nf_izH2]
	movapd xmm3, [esp + nb304nf_ixH2]
	movapd xmm4, [esp + nb304nf_iyH2]
	movapd xmm5, [esp + nb304nf_izH2]
	subsd  xmm0, [esp + nb304nf_jxM]
	subsd  xmm1, [esp + nb304nf_jyM]
	subsd  xmm2, [esp + nb304nf_jzM]
	subsd  xmm3, [esp + nb304nf_jxH1]
	subsd  xmm4, [esp + nb304nf_jyH1]
	subsd  xmm5, [esp + nb304nf_jzH1]
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
	movapd [esp + nb304nf_rsqH2M], xmm0
	movapd [esp + nb304nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb304nf_ixH2]
	movapd xmm1, [esp + nb304nf_iyH2]
	movapd xmm2, [esp + nb304nf_izH2]
	subsd  xmm0, [esp + nb304nf_jxH2]
	subsd  xmm1, [esp + nb304nf_jyH2]
	subsd  xmm2, [esp + nb304nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb304nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb304nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb304nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvH2H2], xmm1
	movapd [esp + nb304nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb304nf_rsqMM]
	movapd xmm4, [esp + nb304nf_rsqMH1]	
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
	movapd  xmm3, [esp + nb304nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb304nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb304nf_half] ;# rinv
	movapd [esp + nb304nf_rinvMM], xmm1
	movapd [esp + nb304nf_rinvMH1], xmm5

	movapd xmm0, [esp + nb304nf_rsqMH2]
	movapd xmm4, [esp + nb304nf_rsqH1M]	
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
	movapd  xmm3, [esp + nb304nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb304nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvMH2], xmm1
	movapd [esp + nb304nf_rinvH1M], xmm5

	movapd xmm0, [esp + nb304nf_rsqH1H1]
	movapd xmm4, [esp + nb304nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb304nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb304nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb304nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb304nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb304nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvH1H1], xmm1
	movapd [esp + nb304nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb304nf_rsqH2M]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb304nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb304nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb304nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb304nf_half] ;# rinv 
	movapd [esp + nb304nf_rinvH2M], xmm1
	
	;# start with MM interaction 
	movapd xmm0, [esp + nb304nf_rinvMM]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqMM] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMM]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# M-H1 interaction 
	movapd xmm0, [esp + nb304nf_rinvMH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqMH1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# M-H2 interaction  
	movapd xmm0, [esp + nb304nf_rinvMH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqMH2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# H1-M interaction 
	movapd xmm0, [esp + nb304nf_rinvH1M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqH1M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb304nf_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb304nf_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# H2-M interaction 
	movapd xmm0, [esp + nb304nf_rinvH2M]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqH2M] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqMH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb304nf_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb304nf_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb304nf_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb304nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb304nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb304nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb304nf_vctot]
    movlpd [esp + nb304nf_vctot], xmm5
	
.nb304nf_updateouterdata:
	;# get group index for i particle 
	;# get n from stack
	mov esi, [esp + nb304nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb304nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb304nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb304nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb304nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb304nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb304nf_n], esi
        jmp .nb304nf_outer
.nb304nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb304nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb304nf_end
        ;# non-zero, do one more workunit
        jmp   .nb304nf_threadloop
.nb304nf_end:
	emms

	mov eax, [esp + nb304nf_nouter]
	mov ebx, [esp + nb304nf_ninner]
	mov ecx, [ebp + nb304nf_outeriter]
	mov edx, [ebp + nb304nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb304nf_salign]
	add esp, eax
	add esp, 728
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
	

