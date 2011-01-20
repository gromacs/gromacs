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

	
.globl nb_kernel312_ia32_sse2
.globl _nb_kernel312_ia32_sse2
nb_kernel312_ia32_sse2:	
_nb_kernel312_ia32_sse2:	
.equiv          nb312_p_nri,            8
.equiv          nb312_iinr,             12
.equiv          nb312_jindex,           16
.equiv          nb312_jjnr,             20
.equiv          nb312_shift,            24
.equiv          nb312_shiftvec,         28
.equiv          nb312_fshift,           32
.equiv          nb312_gid,              36
.equiv          nb312_pos,              40
.equiv          nb312_faction,          44
.equiv          nb312_charge,           48
.equiv          nb312_p_facel,          52
.equiv          nb312_argkrf,           56
.equiv          nb312_argcrf,           60
.equiv          nb312_Vc,               64
.equiv          nb312_type,             68
.equiv          nb312_p_ntype,          72
.equiv          nb312_vdwparam,         76
.equiv          nb312_Vvdw,             80
.equiv          nb312_p_tabscale,       84
.equiv          nb312_VFtab,            88
.equiv          nb312_invsqrta,         92
.equiv          nb312_dvda,             96
.equiv          nb312_p_gbtabscale,     100
.equiv          nb312_GBtab,            104
.equiv          nb312_p_nthreads,       108
.equiv          nb312_count,            112
.equiv          nb312_mtx,              116
.equiv          nb312_outeriter,        120
.equiv          nb312_inneriter,        124
.equiv          nb312_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb312_ixO,              0
.equiv          nb312_iyO,              16
.equiv          nb312_izO,              32
.equiv          nb312_ixH1,             48
.equiv          nb312_iyH1,             64
.equiv          nb312_izH1,             80
.equiv          nb312_ixH2,             96
.equiv          nb312_iyH2,             112
.equiv          nb312_izH2,             128
.equiv          nb312_jxO,              144
.equiv          nb312_jyO,              160
.equiv          nb312_jzO,              176
.equiv          nb312_jxH1,             192
.equiv          nb312_jyH1,             208
.equiv          nb312_jzH1,             224
.equiv          nb312_jxH2,             240
.equiv          nb312_jyH2,             256
.equiv          nb312_jzH2,             272
.equiv          nb312_dxOO,             288
.equiv          nb312_dyOO,             304
.equiv          nb312_dzOO,             320
.equiv          nb312_dxOH1,            336
.equiv          nb312_dyOH1,            352
.equiv          nb312_dzOH1,            368
.equiv          nb312_dxOH2,            384
.equiv          nb312_dyOH2,            400
.equiv          nb312_dzOH2,            416
.equiv          nb312_dxH1O,            432
.equiv          nb312_dyH1O,            448
.equiv          nb312_dzH1O,            464
.equiv          nb312_dxH1H1,           480
.equiv          nb312_dyH1H1,           496
.equiv          nb312_dzH1H1,           512
.equiv          nb312_dxH1H2,           528
.equiv          nb312_dyH1H2,           544
.equiv          nb312_dzH1H2,           560
.equiv          nb312_dxH2O,            576
.equiv          nb312_dyH2O,            592
.equiv          nb312_dzH2O,            608
.equiv          nb312_dxH2H1,           624
.equiv          nb312_dyH2H1,           640
.equiv          nb312_dzH2H1,           656
.equiv          nb312_dxH2H2,           672
.equiv          nb312_dyH2H2,           688
.equiv          nb312_dzH2H2,           704
.equiv          nb312_qqOO,             720
.equiv          nb312_qqOH,             736
.equiv          nb312_qqHH,             752
.equiv          nb312_two,              768
.equiv          nb312_tsc,              784
.equiv          nb312_c6,               800
.equiv          nb312_c12,              816
.equiv          nb312_six,              832
.equiv          nb312_twelve,           848
.equiv          nb312_vctot,            864
.equiv          nb312_Vvdwtot,          880
.equiv          nb312_fixO,             896
.equiv          nb312_fiyO,             912
.equiv          nb312_fizO,             928
.equiv          nb312_fixH1,            944
.equiv          nb312_fiyH1,            960
.equiv          nb312_fizH1,            976
.equiv          nb312_fixH2,            992
.equiv          nb312_fiyH2,            1008
.equiv          nb312_fizH2,            1024
.equiv          nb312_fjxO,             1040
.equiv          nb312_fjyO,             1056
.equiv          nb312_fjzO,             1072
.equiv          nb312_fjxH1,            1088
.equiv          nb312_fjyH1,            1104
.equiv          nb312_fjzH1,            1120
.equiv          nb312_fjxH2,            1136
.equiv          nb312_fjyH2,            1152
.equiv          nb312_fjzH2,            1168
.equiv          nb312_half,             1184
.equiv          nb312_three,            1200
.equiv          nb312_rsqOO,            1216
.equiv          nb312_rsqOH1,           1232
.equiv          nb312_rsqOH2,           1248
.equiv          nb312_rsqH1O,           1264
.equiv          nb312_rsqH1H1,          1280
.equiv          nb312_rsqH1H2,          1296
.equiv          nb312_rsqH2O,           1312
.equiv          nb312_rsqH2H1,          1328
.equiv          nb312_rsqH2H2,          1344
.equiv          nb312_rinvOO,           1360
.equiv          nb312_rinvOH1,          1376
.equiv          nb312_rinvOH2,          1392
.equiv          nb312_rinvH1O,          1408
.equiv          nb312_rinvH1H1,         1424
.equiv          nb312_rinvH1H2,         1440
.equiv          nb312_rinvH2O,          1456
.equiv          nb312_rinvH2H1,         1472
.equiv          nb312_rinvH2H2,         1488
.equiv          nb312_fstmp,            1504
.equiv          nb312_is3,              1520
.equiv          nb312_ii3,              1524
.equiv          nb312_innerjjnr,        1528
.equiv          nb312_innerk,           1532
.equiv          nb312_n,                1536
.equiv          nb312_nn1,              1540
.equiv          nb312_nri,              1544
.equiv          nb312_nouter,           1548
.equiv          nb312_ninner,           1552
.equiv          nb312_salign,           1556
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 1560		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb312_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb312_p_nri]
	mov ecx, [ecx]
	mov [esp + nb312_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb312_nouter], eax
	mov [esp + nb312_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb312_half], eax
	mov [esp + nb312_half+4], ebx
	movsd xmm1, [esp + nb312_half]
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
	movapd [esp + nb312_half], xmm1
	movapd [esp + nb312_two], xmm2
	movapd [esp + nb312_three], xmm3
	movapd [esp + nb312_six], xmm4
	movapd [esp + nb312_twelve], xmm5

	mov eax, [ebp + nb312_p_tabscale]
	movsd xmm5, [eax]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb312_tsc],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb312_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb312_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb312_p_facel]
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
	movapd [esp + nb312_qqOO], xmm3
	movapd [esp + nb312_qqOH], xmm4
	movapd [esp + nb312_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb312_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb312_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb312_vdwparam]
	movlpd xmm0, [eax + edx*8] 
	movlpd xmm1, [eax + edx*8 + 8] 
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb312_c6], xmm0
	movapd [esp + nb312_c12], xmm1

.nb312_threadloop:
        mov   esi, [ebp + nb312_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb312_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb312_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb312_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb312_n], eax
        mov [esp + nb312_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb312_outerstart
        jmp .nb312_end

.nb312_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb312_nouter]
	mov [esp + nb312_nouter], ebx

.nb312_outer:
	mov   eax, [ebp + nb312_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb312_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb312_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb312_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb312_pos]    ;# eax = base of pos[]  
	mov   [esp + nb312_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb312_ixO], xmm3
	movapd [esp + nb312_iyO], xmm4
	movapd [esp + nb312_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 24]
	addsd xmm1, [eax + ebx*8 + 32]
	addsd xmm2, [eax + ebx*8 + 40]		
	addsd xmm3, [eax + ebx*8 + 48]
	addsd xmm4, [eax + ebx*8 + 56]
	addsd xmm5, [eax + ebx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb312_ixH1], xmm0
	movapd [esp + nb312_iyH1], xmm1
	movapd [esp + nb312_izH1], xmm2
	movapd [esp + nb312_ixH2], xmm3
	movapd [esp + nb312_iyH2], xmm4
	movapd [esp + nb312_izH2], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb312_vctot], xmm4
	movapd [esp + nb312_Vvdwtot], xmm4
	movapd [esp + nb312_fixO], xmm4
	movapd [esp + nb312_fiyO], xmm4
	movapd [esp + nb312_fizO], xmm4
	movapd [esp + nb312_fixH1], xmm4
	movapd [esp + nb312_fiyH1], xmm4
	movapd [esp + nb312_fizH1], xmm4
	movapd [esp + nb312_fixH2], xmm4
	movapd [esp + nb312_fiyH2], xmm4
	movapd [esp + nb312_fizH2], xmm4
	
	mov   eax, [ebp + nb312_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb312_pos]
	mov   edi, [ebp + nb312_faction]	
	mov   eax, [ebp + nb312_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb312_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb312_ninner]
	mov   [esp + nb312_ninner], ecx
	add   edx, 0
	mov   [esp + nb312_innerk], edx    ;# number of innerloop atoms 
	jge   .nb312_unroll_loop
	jmp   .nb312_checksingle
.nb312_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb312_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb312_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb312_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movhpd xmm2, [esi + ebx*8]
	movhpd xmm3, [esi + ebx*8 + 8]
	movhpd xmm4, [esi + ebx*8 + 16]
	movhpd xmm5, [esi + ebx*8 + 24]
	movhpd xmm6, [esi + ebx*8 + 32]
	movhpd xmm7, [esi + ebx*8 + 40]
	movapd 	[esp + nb312_jxO], xmm2
	movapd 	[esp + nb312_jyO], xmm3
	movapd 	[esp + nb312_jzO], xmm4
	movapd 	[esp + nb312_jxH1], xmm5
	movapd 	[esp + nb312_jyH1], xmm6
	movapd 	[esp + nb312_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + ebx*8 + 56]
	movhpd xmm4, [esi + ebx*8 + 64]
	movapd 	[esp + nb312_jxH2], xmm2
	movapd 	[esp + nb312_jyH2], xmm3
	movapd 	[esp + nb312_jzH2], xmm4
	
	movapd xmm0, [esp + nb312_ixO]
	movapd xmm1, [esp + nb312_iyO]
	movapd xmm2, [esp + nb312_izO]
	movapd xmm3, [esp + nb312_ixO]
	movapd xmm4, [esp + nb312_iyO]
	movapd xmm5, [esp + nb312_izO]
	subpd  xmm0, [esp + nb312_jxO]
	subpd  xmm1, [esp + nb312_jyO]
	subpd  xmm2, [esp + nb312_jzO]
	subpd  xmm3, [esp + nb312_jxH1]
	subpd  xmm4, [esp + nb312_jyH1]
	subpd  xmm5, [esp + nb312_jzH1]
	movapd [esp + nb312_dxOO], xmm0
	movapd [esp + nb312_dyOO], xmm1
	movapd [esp + nb312_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb312_dxOH1], xmm3
	movapd [esp + nb312_dyOH1], xmm4
	movapd [esp + nb312_dzOH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb312_rsqOO], xmm0
	movapd [esp + nb312_rsqOH1], xmm3

	movapd xmm0, [esp + nb312_ixO]
	movapd xmm1, [esp + nb312_iyO]
	movapd xmm2, [esp + nb312_izO]
	movapd xmm3, [esp + nb312_ixH1]
	movapd xmm4, [esp + nb312_iyH1]
	movapd xmm5, [esp + nb312_izH1]
	subpd  xmm0, [esp + nb312_jxH2]
	subpd  xmm1, [esp + nb312_jyH2]
	subpd  xmm2, [esp + nb312_jzH2]
	subpd  xmm3, [esp + nb312_jxO]
	subpd  xmm4, [esp + nb312_jyO]
	subpd  xmm5, [esp + nb312_jzO]
	movapd [esp + nb312_dxOH2], xmm0
	movapd [esp + nb312_dyOH2], xmm1
	movapd [esp + nb312_dzOH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb312_dxH1O], xmm3
	movapd [esp + nb312_dyH1O], xmm4
	movapd [esp + nb312_dzH1O], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb312_rsqOH2], xmm0
	movapd [esp + nb312_rsqH1O], xmm3

	movapd xmm0, [esp + nb312_ixH1]
	movapd xmm1, [esp + nb312_iyH1]
	movapd xmm2, [esp + nb312_izH1]
	movapd xmm3, [esp + nb312_ixH1]
	movapd xmm4, [esp + nb312_iyH1]
	movapd xmm5, [esp + nb312_izH1]
	subpd  xmm0, [esp + nb312_jxH1]
	subpd  xmm1, [esp + nb312_jyH1]
	subpd  xmm2, [esp + nb312_jzH1]
	subpd  xmm3, [esp + nb312_jxH2]
	subpd  xmm4, [esp + nb312_jyH2]
	subpd  xmm5, [esp + nb312_jzH2]
	movapd [esp + nb312_dxH1H1], xmm0
	movapd [esp + nb312_dyH1H1], xmm1
	movapd [esp + nb312_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb312_dxH1H2], xmm3
	movapd [esp + nb312_dyH1H2], xmm4
	movapd [esp + nb312_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb312_rsqH1H1], xmm0
	movapd [esp + nb312_rsqH1H2], xmm3

	movapd xmm0, [esp + nb312_ixH2]
	movapd xmm1, [esp + nb312_iyH2]
	movapd xmm2, [esp + nb312_izH2]
	movapd xmm3, [esp + nb312_ixH2]
	movapd xmm4, [esp + nb312_iyH2]
	movapd xmm5, [esp + nb312_izH2]
	subpd  xmm0, [esp + nb312_jxO]
	subpd  xmm1, [esp + nb312_jyO]
	subpd  xmm2, [esp + nb312_jzO]
	subpd  xmm3, [esp + nb312_jxH1]
	subpd  xmm4, [esp + nb312_jyH1]
	subpd  xmm5, [esp + nb312_jzH1]
	movapd [esp + nb312_dxH2O], xmm0
	movapd [esp + nb312_dyH2O], xmm1
	movapd [esp + nb312_dzH2O], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb312_dxH2H1], xmm3
	movapd [esp + nb312_dyH2H1], xmm4
	movapd [esp + nb312_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb312_rsqH2O], xmm0
	movapd [esp + nb312_rsqH2H1], xmm4

	movapd xmm0, [esp + nb312_ixH2]
	movapd xmm1, [esp + nb312_iyH2]
	movapd xmm2, [esp + nb312_izH2]
	subpd  xmm0, [esp + nb312_jxH2]
	subpd  xmm1, [esp + nb312_jyH2]
	subpd  xmm2, [esp + nb312_jzH2]
	movapd [esp + nb312_dxH2H2], xmm0
	movapd [esp + nb312_dyH2H2], xmm1
	movapd [esp + nb312_dzH2H2], xmm2
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb312_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb312_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312_half] ;# iter1 
	mulpd   xmm7, [esp + nb312_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312_half] ;# rinv 
	mulpd   xmm5, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvH2H2], xmm1
	movapd [esp + nb312_rinvH2H1], xmm5

	movapd xmm0, [esp + nb312_rsqOO]
	movapd xmm4, [esp + nb312_rsqOH1]	
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
	movapd  xmm3, [esp + nb312_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb312_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312_half] ;# rinv 
	mulpd   xmm5, [esp + nb312_half] ;# rinv
	movapd [esp + nb312_rinvOO], xmm1
	movapd [esp + nb312_rinvOH1], xmm5

	movapd xmm0, [esp + nb312_rsqOH2]
	movapd xmm4, [esp + nb312_rsqH1O]	
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
	movapd  xmm3, [esp + nb312_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312_half] ;# iter1 
	mulpd   xmm7, [esp + nb312_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312_half] ;# rinv 
	mulpd   xmm5, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvOH2], xmm1
	movapd [esp + nb312_rinvH1O], xmm5

	movapd xmm0, [esp + nb312_rsqH1H1]
	movapd xmm4, [esp + nb312_rsqH1H2]	
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
	movapd  xmm3, [esp + nb312_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312_half] ;# iter1a 
	mulpd   xmm7, [esp + nb312_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312_half] ;# rinv 
	mulpd   xmm5, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvH1H1], xmm1
	movapd [esp + nb312_rinvH1H2], xmm5

	movapd xmm0, [esp + nb312_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb312_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb312_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb312_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm0, [esp + nb312_rinvOO]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqOO] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	movd mm0, eax	
	movd mm1, ebx
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOO]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addpd  xmm5, [esp + nb312_vctot]
	xorpd  xmm2, xmm2
    movapd [esp + nb312_vctot], xmm5
	mulpd  xmm3, [esp + nb312_tsc]

	;# start doing lj 
	movapd xmm2, xmm0
	mulpd  xmm2, xmm2
	movapd xmm1, xmm2
	mulpd  xmm1, xmm2
	mulpd  xmm1, xmm2	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm1, [esp + nb312_c6]
	mulpd  xmm2, [esp + nb312_c12]
	movapd xmm4, xmm2
	subpd  xmm4, xmm1
	addpd  xmm4, [esp + nb312_Vvdwtot]
	mulpd  xmm1, [esp + nb312_six]
	mulpd  xmm2, [esp + nb312_twelve]
	movapd [esp + nb312_Vvdwtot], xmm4
	subpd  xmm2, xmm1
	mulpd  xmm2, xmm0

	subpd  xmm2, xmm3
	mulpd  xmm0, xmm2
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0		

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb312_dxOO]
	mulpd xmm1, [esp + nb312_dyOO]
	mulpd xmm2, [esp + nb312_dzOO]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixO]
	addpd xmm1, [esp + nb312_fiyO]
	addpd xmm2, [esp + nb312_fizO]
	movapd [esp + nb312_fjxO], xmm3
	movapd [esp + nb312_fjyO], xmm4
	movapd [esp + nb312_fjzO], xmm5
	movapd [esp + nb312_fixO], xmm0
	movapd [esp + nb312_fiyO], xmm1
	movapd [esp + nb312_fizO], xmm2

	;# O-H1 interaction 
	movapd xmm0, [esp + nb312_rinvOH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqOH1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb312_dxOH1]
	mulpd xmm1, [esp + nb312_dyOH1]
	mulpd xmm2, [esp + nb312_dzOH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixO]
	addpd xmm1, [esp + nb312_fiyO]
	addpd xmm2, [esp + nb312_fizO]
	movapd [esp + nb312_fjxH1], xmm3
	movapd [esp + nb312_fjyH1], xmm4
	movapd [esp + nb312_fjzH1], xmm5
	movapd [esp + nb312_fixO], xmm0
	movapd [esp + nb312_fiyO], xmm1
	movapd [esp + nb312_fizO], xmm2

	;# O-H2 interaction  
	movapd xmm0, [esp + nb312_rinvOH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqOH2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb312_dxOH2]
	mulpd xmm1, [esp + nb312_dyOH2]
	mulpd xmm2, [esp + nb312_dzOH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixO]
	addpd xmm1, [esp + nb312_fiyO]
	addpd xmm2, [esp + nb312_fizO]
	movapd [esp + nb312_fjxH2], xmm3
	movapd [esp + nb312_fjyH2], xmm4
	movapd [esp + nb312_fjzH2], xmm5
	movapd [esp + nb312_fixO], xmm0
	movapd [esp + nb312_fiyO], xmm1
	movapd [esp + nb312_fizO], xmm2

	;# H1-O interaction 
	movapd xmm0, [esp + nb312_rinvH1O]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqH1O] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxO]
	movapd xmm4, [esp + nb312_fjyO]
	movapd xmm5, [esp + nb312_fjzO]
	mulpd xmm0, [esp + nb312_dxH1O]
	mulpd xmm1, [esp + nb312_dyH1O]
	mulpd xmm2, [esp + nb312_dzH1O]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixH1]
	addpd xmm1, [esp + nb312_fiyH1]
	addpd xmm2, [esp + nb312_fizH1]
	movapd [esp + nb312_fjxO], xmm3
	movapd [esp + nb312_fjyO], xmm4
	movapd [esp + nb312_fjzO], xmm5
	movapd [esp + nb312_fixH1], xmm0
	movapd [esp + nb312_fiyH1], xmm1
	movapd [esp + nb312_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb312_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH1]
	movapd xmm4, [esp + nb312_fjyH1]
	movapd xmm5, [esp + nb312_fjzH1]
	mulpd xmm0, [esp + nb312_dxH1H1]
	mulpd xmm1, [esp + nb312_dyH1H1]
	mulpd xmm2, [esp + nb312_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixH1]
	addpd xmm1, [esp + nb312_fiyH1]
	addpd xmm2, [esp + nb312_fizH1]
	movapd [esp + nb312_fjxH1], xmm3
	movapd [esp + nb312_fjyH1], xmm4
	movapd [esp + nb312_fjzH1], xmm5
	movapd [esp + nb312_fixH1], xmm0
	movapd [esp + nb312_fiyH1], xmm1
	movapd [esp + nb312_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb312_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH2]
	movapd xmm4, [esp + nb312_fjyH2]
	movapd xmm5, [esp + nb312_fjzH2]
	mulpd xmm0, [esp + nb312_dxH1H2]
	mulpd xmm1, [esp + nb312_dyH1H2]
	mulpd xmm2, [esp + nb312_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixH1]
	addpd xmm1, [esp + nb312_fiyH1]
	addpd xmm2, [esp + nb312_fizH1]
	movapd [esp + nb312_fjxH2], xmm3
	movapd [esp + nb312_fjyH2], xmm4
	movapd [esp + nb312_fjzH2], xmm5
	movapd [esp + nb312_fixH1], xmm0
	movapd [esp + nb312_fiyH1], xmm1
	movapd [esp + nb312_fizH1], xmm2

	;# H2-O interaction 
	movapd xmm0, [esp + nb312_rinvH2O]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqH2O] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1

	movapd xmm3, [esp + nb312_fjxO]
	movapd xmm4, [esp + nb312_fjyO]
	movapd xmm5, [esp + nb312_fjzO]
	mulpd xmm0, [esp + nb312_dxH2O]
	mulpd xmm1, [esp + nb312_dyH2O]
	mulpd xmm2, [esp + nb312_dzH2O]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixH2]
	addpd xmm1, [esp + nb312_fiyH2]
	addpd xmm2, [esp + nb312_fizH2]
	movapd [esp + nb312_fjxO], xmm3
	movapd [esp + nb312_fjyO], xmm4
	movapd [esp + nb312_fjzO], xmm5
	movapd [esp + nb312_fixH2], xmm0
	movapd [esp + nb312_fiyH2], xmm1
	movapd [esp + nb312_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb312_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH1]
	movapd xmm4, [esp + nb312_fjyH1]
	movapd xmm5, [esp + nb312_fjzH1]
	mulpd xmm0, [esp + nb312_dxH2H1]
	mulpd xmm1, [esp + nb312_dyH2H1]
	mulpd xmm2, [esp + nb312_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixH2]
	addpd xmm1, [esp + nb312_fiyH2]
	addpd xmm2, [esp + nb312_fizH2]
	movapd [esp + nb312_fjxH1], xmm3
	movapd [esp + nb312_fjyH1], xmm4
	movapd [esp + nb312_fjzH1], xmm5
	movapd [esp + nb312_fixH2], xmm0
	movapd [esp + nb312_fiyH2], xmm1
	movapd [esp + nb312_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb312_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]
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
	mulpd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
	
    addpd  xmm5, [esp + nb312_vctot]
    movapd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulpd  xmm3,  [esp + nb312_tsc]
	mulpd  xmm3, xmm0
	subpd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH2]
	movapd xmm4, [esp + nb312_fjyH2]
	movapd xmm5, [esp + nb312_fjzH2]
	mulpd xmm0, [esp + nb312_dxH2H2]
	mulpd xmm1, [esp + nb312_dyH2H2]
	mulpd xmm2, [esp + nb312_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb312_fixH2]
	addpd xmm1, [esp + nb312_fiyH2]
	addpd xmm2, [esp + nb312_fizH2]
	movapd [esp + nb312_fjxH2], xmm3
	movapd [esp + nb312_fjyH2], xmm4
	movapd [esp + nb312_fjzH2], xmm5
	movapd [esp + nb312_fixH2], xmm0
	movapd [esp + nb312_fiyH2], xmm1
	movapd [esp + nb312_fizH2], xmm2

	mov edi, [ebp + nb312_faction]

	movd eax, mm0
	movd ebx, mm1
	
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]
	movlpd xmm3, [edi + eax*8 + 24]
	movlpd xmm4, [edi + eax*8 + 32]
	movlpd xmm5, [edi + eax*8 + 40]
	movlpd xmm6, [edi + eax*8 + 48]
	movlpd xmm7, [edi + eax*8 + 56]
	movhpd xmm0, [edi + ebx*8]
	movhpd xmm1, [edi + ebx*8 + 8]
	movhpd xmm2, [edi + ebx*8 + 16]
	movhpd xmm3, [edi + ebx*8 + 24]
	movhpd xmm4, [edi + ebx*8 + 32]
	movhpd xmm5, [edi + ebx*8 + 40]
	movhpd xmm6, [edi + ebx*8 + 48]
	movhpd xmm7, [edi + ebx*8 + 56]
	addpd xmm0, [esp + nb312_fjxO]
	addpd xmm1, [esp + nb312_fjyO]
	addpd xmm2, [esp + nb312_fjzO]
	addpd xmm3, [esp + nb312_fjxH1]
	addpd xmm4, [esp + nb312_fjyH1]
	addpd xmm5, [esp + nb312_fjzH1]
	addpd xmm6, [esp + nb312_fjxH2]
	addpd xmm7, [esp + nb312_fjyH2]
	movlpd [edi + eax*8], xmm0
	movlpd [edi + eax*8 + 8], xmm1
	movlpd [edi + eax*8 + 16], xmm2
	movlpd [edi + eax*8 + 24], xmm3
	movlpd [edi + eax*8 + 32], xmm4
	movlpd [edi + eax*8 + 40], xmm5
	movlpd [edi + eax*8 + 48], xmm6
	movlpd [edi + eax*8 + 56], xmm7
	movhpd [edi + ebx*8], xmm0
	movhpd [edi + ebx*8 + 8], xmm1
	movhpd [edi + ebx*8 + 16], xmm2
	movhpd [edi + ebx*8 + 24], xmm3
	movhpd [edi + ebx*8 + 32], xmm4
	movhpd [edi + ebx*8 + 40], xmm5
	movhpd [edi + ebx*8 + 48], xmm6
	movhpd [edi + ebx*8 + 56], xmm7
	
	movlpd xmm0, [edi + eax*8 + 64]
	movhpd xmm0, [edi + ebx*8 + 64]
	addpd xmm0, [esp + nb312_fjzH2]
	movlpd [edi + eax*8 + 64], xmm0
	movhpd [edi + ebx*8 + 64], xmm0

	;# should we do one more iteration? 
	sub dword ptr [esp + nb312_innerk],  2
	jl    .nb312_checksingle
	jmp   .nb312_unroll_loop
.nb312_checksingle:
	mov   edx, [esp + nb312_innerk]
	and   edx, 1
	jnz   .nb312_dosingle
	jmp   .nb312_updateouterdata
.nb312_dosingle:
	mov   edx, [esp + nb312_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb312_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movapd 	[esp + nb312_jxO], xmm2
	movapd 	[esp + nb312_jyO], xmm3
	movapd 	[esp + nb312_jzO], xmm4
	movapd 	[esp + nb312_jxH1], xmm5
	movapd 	[esp + nb312_jyH1], xmm6
	movapd 	[esp + nb312_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movapd 	[esp + nb312_jxH2], xmm2
	movapd 	[esp + nb312_jyH2], xmm3
	movapd 	[esp + nb312_jzH2], xmm4
	
	movapd xmm0, [esp + nb312_ixO]
	movapd xmm1, [esp + nb312_iyO]
	movapd xmm2, [esp + nb312_izO]
	movapd xmm3, [esp + nb312_ixO]
	movapd xmm4, [esp + nb312_iyO]
	movapd xmm5, [esp + nb312_izO]
	subsd  xmm0, [esp + nb312_jxO]
	subsd  xmm1, [esp + nb312_jyO]
	subsd  xmm2, [esp + nb312_jzO]
	subsd  xmm3, [esp + nb312_jxH1]
	subsd  xmm4, [esp + nb312_jyH1]
	subsd  xmm5, [esp + nb312_jzH1]
	movapd [esp + nb312_dxOO], xmm0
	movapd [esp + nb312_dyOO], xmm1
	movapd [esp + nb312_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb312_dxOH1], xmm3
	movapd [esp + nb312_dyOH1], xmm4
	movapd [esp + nb312_dzOH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb312_rsqOO], xmm0
	movapd [esp + nb312_rsqOH1], xmm3

	movapd xmm0, [esp + nb312_ixO]
	movapd xmm1, [esp + nb312_iyO]
	movapd xmm2, [esp + nb312_izO]
	movapd xmm3, [esp + nb312_ixH1]
	movapd xmm4, [esp + nb312_iyH1]
	movapd xmm5, [esp + nb312_izH1]
	subsd  xmm0, [esp + nb312_jxH2]
	subsd  xmm1, [esp + nb312_jyH2]
	subsd  xmm2, [esp + nb312_jzH2]
	subsd  xmm3, [esp + nb312_jxO]
	subsd  xmm4, [esp + nb312_jyO]
	subsd  xmm5, [esp + nb312_jzO]
	movapd [esp + nb312_dxOH2], xmm0
	movapd [esp + nb312_dyOH2], xmm1
	movapd [esp + nb312_dzOH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb312_dxH1O], xmm3
	movapd [esp + nb312_dyH1O], xmm4
	movapd [esp + nb312_dzH1O], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb312_rsqOH2], xmm0
	movapd [esp + nb312_rsqH1O], xmm3

	movapd xmm0, [esp + nb312_ixH1]
	movapd xmm1, [esp + nb312_iyH1]
	movapd xmm2, [esp + nb312_izH1]
	movapd xmm3, [esp + nb312_ixH1]
	movapd xmm4, [esp + nb312_iyH1]
	movapd xmm5, [esp + nb312_izH1]
	subsd  xmm0, [esp + nb312_jxH1]
	subsd  xmm1, [esp + nb312_jyH1]
	subsd  xmm2, [esp + nb312_jzH1]
	subsd  xmm3, [esp + nb312_jxH2]
	subsd  xmm4, [esp + nb312_jyH2]
	subsd  xmm5, [esp + nb312_jzH2]
	movapd [esp + nb312_dxH1H1], xmm0
	movapd [esp + nb312_dyH1H1], xmm1
	movapd [esp + nb312_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb312_dxH1H2], xmm3
	movapd [esp + nb312_dyH1H2], xmm4
	movapd [esp + nb312_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb312_rsqH1H1], xmm0
	movapd [esp + nb312_rsqH1H2], xmm3

	movapd xmm0, [esp + nb312_ixH2]
	movapd xmm1, [esp + nb312_iyH2]
	movapd xmm2, [esp + nb312_izH2]
	movapd xmm3, [esp + nb312_ixH2]
	movapd xmm4, [esp + nb312_iyH2]
	movapd xmm5, [esp + nb312_izH2]
	subsd  xmm0, [esp + nb312_jxO]
	subsd  xmm1, [esp + nb312_jyO]
	subsd  xmm2, [esp + nb312_jzO]
	subsd  xmm3, [esp + nb312_jxH1]
	subsd  xmm4, [esp + nb312_jyH1]
	subsd  xmm5, [esp + nb312_jzH1]
	movapd [esp + nb312_dxH2O], xmm0
	movapd [esp + nb312_dyH2O], xmm1
	movapd [esp + nb312_dzH2O], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb312_dxH2H1], xmm3
	movapd [esp + nb312_dyH2H1], xmm4
	movapd [esp + nb312_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb312_rsqH2O], xmm0
	movapd [esp + nb312_rsqH2H1], xmm4

	movapd xmm0, [esp + nb312_ixH2]
	movapd xmm1, [esp + nb312_iyH2]
	movapd xmm2, [esp + nb312_izH2]
	subsd  xmm0, [esp + nb312_jxH2]
	subsd  xmm1, [esp + nb312_jyH2]
	subsd  xmm2, [esp + nb312_jzH2]
	movapd [esp + nb312_dxH2H2], xmm0
	movapd [esp + nb312_dyH2H2], xmm1
	movapd [esp + nb312_dzH2H2], xmm2
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb312_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb312_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312_half] ;# iter1 
	mulsd   xmm7, [esp + nb312_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312_half] ;# rinv 
	mulsd   xmm5, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvH2H2], xmm1
	movapd [esp + nb312_rinvH2H1], xmm5

	movapd xmm0, [esp + nb312_rsqOO]
	movapd xmm4, [esp + nb312_rsqOH1]	
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
	movapd  xmm3, [esp + nb312_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb312_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312_half] ;# rinv 
	mulsd   xmm5, [esp + nb312_half] ;# rinv
	movapd [esp + nb312_rinvOO], xmm1
	movapd [esp + nb312_rinvOH1], xmm5

	movapd xmm0, [esp + nb312_rsqOH2]
	movapd xmm4, [esp + nb312_rsqH1O]	
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
	movapd  xmm3, [esp + nb312_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312_half] ;# iter1 
	mulsd   xmm7, [esp + nb312_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312_half] ;# rinv 
	mulsd   xmm5, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvOH2], xmm1
	movapd [esp + nb312_rinvH1O], xmm5

	movapd xmm0, [esp + nb312_rsqH1H1]
	movapd xmm4, [esp + nb312_rsqH1H2]	
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
	movapd  xmm3, [esp + nb312_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312_half] ;# iter1a 
	mulsd   xmm7, [esp + nb312_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312_half] ;# rinv 
	mulsd   xmm5, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvH1H1], xmm1
	movapd [esp + nb312_rinvH1H2], xmm5

	movapd xmm0, [esp + nb312_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb312_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb312_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb312_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb312_half] ;# rinv 
	movapd [esp + nb312_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm0, [esp + nb312_rinvOO]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqOO] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOO]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addsd  xmm5, [esp + nb312_vctot]
	xorpd  xmm2, xmm2
    movlpd [esp + nb312_vctot], xmm5
	mulsd  xmm3, [esp + nb312_tsc]

	;# start doing lj 
	movapd xmm2, xmm0
	mulsd  xmm2, xmm2
	movapd xmm1, xmm2
	mulsd  xmm1, xmm2
	mulsd  xmm1, xmm2	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm1, [esp + nb312_c6]
	mulsd  xmm2, [esp + nb312_c12]
	movapd xmm4, xmm2
	subsd  xmm4, xmm1
	addsd  xmm4, [esp + nb312_Vvdwtot]
	mulsd  xmm1, [esp + nb312_six]
	mulsd  xmm2, [esp + nb312_twelve]
	movlpd [esp + nb312_Vvdwtot], xmm4
	subsd  xmm2, xmm1
	mulsd  xmm2, xmm0

	subsd  xmm2, xmm3
	mulsd  xmm0, xmm2
	
	movapd xmm1, xmm0
	movapd xmm2, xmm0		

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb312_dxOO]
	mulpd xmm1, [esp + nb312_dyOO]
	mulpd xmm2, [esp + nb312_dzOO]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixO]
	addsd xmm1, [esp + nb312_fiyO]
	addsd xmm2, [esp + nb312_fizO]
	movlpd [esp + nb312_fjxO], xmm3
	movlpd [esp + nb312_fjyO], xmm4
	movlpd [esp + nb312_fjzO], xmm5
	movlpd [esp + nb312_fixO], xmm0
	movlpd [esp + nb312_fiyO], xmm1
	movlpd [esp + nb312_fizO], xmm2

	;# O-H1 interaction 
	movapd xmm0, [esp + nb312_rinvOH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqOH1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb312_dxOH1]
	mulsd xmm1, [esp + nb312_dyOH1]
	mulsd xmm2, [esp + nb312_dzOH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixO]
	addsd xmm1, [esp + nb312_fiyO]
	addsd xmm2, [esp + nb312_fizO]
	movlpd [esp + nb312_fjxH1], xmm3
	movlpd [esp + nb312_fjyH1], xmm4
	movlpd [esp + nb312_fjzH1], xmm5
	movlpd [esp + nb312_fixO], xmm0
	movlpd [esp + nb312_fiyO], xmm1
	movlpd [esp + nb312_fizO], xmm2

	;# O-H2 interaction  
	movapd xmm0, [esp + nb312_rinvOH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqOH2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb312_dxOH2]
	mulsd xmm1, [esp + nb312_dyOH2]
	mulsd xmm2, [esp + nb312_dzOH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixO]
	addsd xmm1, [esp + nb312_fiyO]
	addsd xmm2, [esp + nb312_fizO]
	movlpd [esp + nb312_fjxH2], xmm3
	movlpd [esp + nb312_fjyH2], xmm4
	movlpd [esp + nb312_fjzH2], xmm5
	movlpd [esp + nb312_fixO], xmm0
	movlpd [esp + nb312_fiyO], xmm1
	movlpd [esp + nb312_fizO], xmm2

	;# H1-O interaction 
	movapd xmm0, [esp + nb312_rinvH1O]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqH1O] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxO]
	movapd xmm4, [esp + nb312_fjyO]
	movapd xmm5, [esp + nb312_fjzO]
	mulsd xmm0, [esp + nb312_dxH1O]
	mulsd xmm1, [esp + nb312_dyH1O]
	mulsd xmm2, [esp + nb312_dzH1O]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixH1]
	addsd xmm1, [esp + nb312_fiyH1]
	addsd xmm2, [esp + nb312_fizH1]
	movlpd [esp + nb312_fjxO], xmm3
	movlpd [esp + nb312_fjyO], xmm4
	movlpd [esp + nb312_fjzO], xmm5
	movlpd [esp + nb312_fixH1], xmm0
	movlpd [esp + nb312_fiyH1], xmm1
	movlpd [esp + nb312_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb312_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH1]
	movapd xmm4, [esp + nb312_fjyH1]
	movapd xmm5, [esp + nb312_fjzH1]
	mulsd xmm0, [esp + nb312_dxH1H1]
	mulsd xmm1, [esp + nb312_dyH1H1]
	mulsd xmm2, [esp + nb312_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixH1]
	addsd xmm1, [esp + nb312_fiyH1]
	addsd xmm2, [esp + nb312_fizH1]
	movlpd [esp + nb312_fjxH1], xmm3
	movlpd [esp + nb312_fjyH1], xmm4
	movlpd [esp + nb312_fjzH1], xmm5
	movlpd [esp + nb312_fixH1], xmm0
	movlpd [esp + nb312_fiyH1], xmm1
	movlpd [esp + nb312_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb312_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH2]
	movapd xmm4, [esp + nb312_fjyH2]
	movapd xmm5, [esp + nb312_fjzH2]
	mulsd xmm0, [esp + nb312_dxH1H2]
	mulsd xmm1, [esp + nb312_dyH1H2]
	mulsd xmm2, [esp + nb312_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixH1]
	addsd xmm1, [esp + nb312_fiyH1]
	addsd xmm2, [esp + nb312_fizH1]
	movlpd [esp + nb312_fjxH2], xmm3
	movlpd [esp + nb312_fjyH2], xmm4
	movlpd [esp + nb312_fjzH2], xmm5
	movlpd [esp + nb312_fixH1], xmm0
	movlpd [esp + nb312_fiyH1], xmm1
	movlpd [esp + nb312_fizH1], xmm2

	;# H2-O interaction 
	movapd xmm0, [esp + nb312_rinvH2O]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqH2O] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqOH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1

	movapd xmm3, [esp + nb312_fjxO]
	movapd xmm4, [esp + nb312_fjyO]
	movapd xmm5, [esp + nb312_fjzO]
	mulsd xmm0, [esp + nb312_dxH2O]
	mulsd xmm1, [esp + nb312_dyH2O]
	mulsd xmm2, [esp + nb312_dzH2O]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixH2]
	addsd xmm1, [esp + nb312_fiyH2]
	addsd xmm2, [esp + nb312_fizH2]
	movlpd [esp + nb312_fjxO], xmm3
	movlpd [esp + nb312_fjyO], xmm4
	movlpd [esp + nb312_fjzO], xmm5
	movlpd [esp + nb312_fixH2], xmm0
	movlpd [esp + nb312_fiyH2], xmm1
	movlpd [esp + nb312_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb312_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH1]
	movapd xmm4, [esp + nb312_fjyH1]
	movapd xmm5, [esp + nb312_fjzH1]
	mulsd xmm0, [esp + nb312_dxH2H1]
	mulsd xmm1, [esp + nb312_dyH2H1]
	mulsd xmm2, [esp + nb312_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixH2]
	addsd xmm1, [esp + nb312_fiyH2]
	addsd xmm2, [esp + nb312_fizH2]
	movlpd [esp + nb312_fjxH1], xmm3
	movlpd [esp + nb312_fjyH1], xmm4
	movlpd [esp + nb312_fjzH1], xmm5
	movlpd [esp + nb312_fixH2], xmm0
	movlpd [esp + nb312_fiyH2], xmm1
	movlpd [esp + nb312_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb312_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb312_two]	;# two*Heps2 
	movapd xmm3, [esp + nb312_qqHH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 

    addsd  xmm5, [esp + nb312_vctot]
    movlpd [esp + nb312_vctot], xmm5
	xorpd  xmm1, xmm1
	mulsd  xmm3,  [esp + nb312_tsc]
	mulsd  xmm3, xmm0
	subsd  xmm1, xmm3

	movapd xmm0, xmm1
	movapd xmm2, xmm1
	
	movapd xmm3, [esp + nb312_fjxH2]
	movapd xmm4, [esp + nb312_fjyH2]
	movapd xmm5, [esp + nb312_fjzH2]
	mulsd xmm0, [esp + nb312_dxH2H2]
	mulsd xmm1, [esp + nb312_dyH2H2]
	mulsd xmm2, [esp + nb312_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb312_fixH2]
	addsd xmm1, [esp + nb312_fiyH2]
	addsd xmm2, [esp + nb312_fizH2]
	movlpd [esp + nb312_fjxH2], xmm3
	movlpd [esp + nb312_fjyH2], xmm4
	movlpd [esp + nb312_fjzH2], xmm5
	movlpd [esp + nb312_fixH2], xmm0
	movlpd [esp + nb312_fiyH2], xmm1
	movlpd [esp + nb312_fizH2], xmm2

	mov edi, [ebp + nb312_faction]

	movd eax, mm0
	
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]
	movlpd xmm3, [edi + eax*8 + 24]
	movlpd xmm4, [edi + eax*8 + 32]
	movlpd xmm5, [edi + eax*8 + 40]
	movlpd xmm6, [edi + eax*8 + 48]
	movlpd xmm7, [edi + eax*8 + 56]
	addsd xmm0, [esp + nb312_fjxO]
	addsd xmm1, [esp + nb312_fjyO]
	addsd xmm2, [esp + nb312_fjzO]
	addsd xmm3, [esp + nb312_fjxH1]
	addsd xmm4, [esp + nb312_fjyH1]
	addsd xmm5, [esp + nb312_fjzH1]
	addsd xmm6, [esp + nb312_fjxH2]
	addsd xmm7, [esp + nb312_fjyH2]
	movlpd [edi + eax*8], xmm0
	movlpd [edi + eax*8 + 8], xmm1
	movlpd [edi + eax*8 + 16], xmm2
	movlpd [edi + eax*8 + 24], xmm3
	movlpd [edi + eax*8 + 32], xmm4
	movlpd [edi + eax*8 + 40], xmm5
	movlpd [edi + eax*8 + 48], xmm6
	movlpd [edi + eax*8 + 56], xmm7

	movlpd xmm0, [edi + eax*8 + 64]
	addsd xmm0, [esp + nb312_fjzH2]
	movlpd [edi + eax*8 + 64], xmm0
	
.nb312_updateouterdata:
	mov   ecx, [esp + nb312_ii3]
	mov   edi, [ebp + nb312_faction]
	mov   esi, [ebp + nb312_fshift]
	mov   edx, [esp + nb312_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb312_fixO]
	movapd xmm1, [esp + nb312_fiyO]
	movapd xmm2, [esp + nb312_fizO]

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
	movapd xmm0, [esp + nb312_fixH1]
	movapd xmm1, [esp + nb312_fiyH1]
	movapd xmm2, [esp + nb312_fizH1]

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
	movapd xmm0, [esp + nb312_fixH2]
	movapd xmm1, [esp + nb312_fiyH2]
	movapd xmm2, [esp + nb312_fizH2]

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
	mov esi, [esp + nb312_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb312_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb312_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb312_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb312_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb312_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb312_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb312_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb312_n], esi
        jmp .nb312_outer
.nb312_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb312_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb312_end
        ;# non-zero, do one more workunit
        jmp   .nb312_threadloop
.nb312_end:
	emms

	mov eax, [esp + nb312_nouter]
	mov ebx, [esp + nb312_ninner]
	mov ecx, [ebp + nb312_outeriter]
	mov edx, [ebp + nb312_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb312_salign]
	add esp, eax
	add esp, 1560
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	



	
.globl nb_kernel312nf_ia32_sse2
.globl _nb_kernel312nf_ia32_sse2
nb_kernel312nf_ia32_sse2:	
_nb_kernel312nf_ia32_sse2:	
.equiv          nb312nf_p_nri,          8
.equiv          nb312nf_iinr,           12
.equiv          nb312nf_jindex,         16
.equiv          nb312nf_jjnr,           20
.equiv          nb312nf_shift,          24
.equiv          nb312nf_shiftvec,       28
.equiv          nb312nf_fshift,         32
.equiv          nb312nf_gid,            36
.equiv          nb312nf_pos,            40
.equiv          nb312nf_faction,        44
.equiv          nb312nf_charge,         48
.equiv          nb312nf_p_facel,        52
.equiv          nb312nf_argkrf,         56
.equiv          nb312nf_argcrf,         60
.equiv          nb312nf_Vc,             64
.equiv          nb312nf_type,           68
.equiv          nb312nf_p_ntype,        72
.equiv          nb312nf_vdwparam,       76
.equiv          nb312nf_Vvdw,           80
.equiv          nb312nf_p_tabscale,     84
.equiv          nb312nf_VFtab,          88
.equiv          nb312nf_invsqrta,       92
.equiv          nb312nf_dvda,           96
.equiv          nb312nf_p_gbtabscale,   100
.equiv          nb312nf_GBtab,          104
.equiv          nb312nf_p_nthreads,     108
.equiv          nb312nf_count,          112
.equiv          nb312nf_mtx,            116
.equiv          nb312nf_outeriter,      120
.equiv          nb312nf_inneriter,      124
.equiv          nb312nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb312nf_ixO,            0
.equiv          nb312nf_iyO,            16
.equiv          nb312nf_izO,            32
.equiv          nb312nf_ixH1,           48
.equiv          nb312nf_iyH1,           64
.equiv          nb312nf_izH1,           80
.equiv          nb312nf_ixH2,           96
.equiv          nb312nf_iyH2,           112
.equiv          nb312nf_izH2,           128
.equiv          nb312nf_jxO,            144
.equiv          nb312nf_jyO,            160
.equiv          nb312nf_jzO,            176
.equiv          nb312nf_jxH1,           192
.equiv          nb312nf_jyH1,           208
.equiv          nb312nf_jzH1,           224
.equiv          nb312nf_jxH2,           240
.equiv          nb312nf_jyH2,           256
.equiv          nb312nf_jzH2,           272
.equiv          nb312nf_qqOO,           288
.equiv          nb312nf_qqOH,           304
.equiv          nb312nf_qqHH,           320
.equiv          nb312nf_tsc,            336
.equiv          nb312nf_c6,             352
.equiv          nb312nf_c12,            368
.equiv          nb312nf_vctot,          384
.equiv          nb312nf_Vvdwtot,        400
.equiv          nb312nf_half,           416
.equiv          nb312nf_three,          432
.equiv          nb312nf_rsqOO,          448
.equiv          nb312nf_rsqOH1,         464
.equiv          nb312nf_rsqOH2,         480
.equiv          nb312nf_rsqH1O,         496
.equiv          nb312nf_rsqH1H1,        512
.equiv          nb312nf_rsqH1H2,        528
.equiv          nb312nf_rsqH2O,         544
.equiv          nb312nf_rsqH2H1,        560
.equiv          nb312nf_rsqH2H2,        576
.equiv          nb312nf_rinvOO,         592
.equiv          nb312nf_rinvOH1,        608
.equiv          nb312nf_rinvOH2,        624
.equiv          nb312nf_rinvH1O,        640
.equiv          nb312nf_rinvH1H1,       656
.equiv          nb312nf_rinvH1H2,       672
.equiv          nb312nf_rinvH2O,        688
.equiv          nb312nf_rinvH2H1,       704
.equiv          nb312nf_rinvH2H2,       720
.equiv          nb312nf_is3,            736
.equiv          nb312nf_ii3,            740
.equiv          nb312nf_innerjjnr,      744
.equiv          nb312nf_innerk,         748
.equiv          nb312nf_n,              752
.equiv          nb312nf_nn1,            756
.equiv          nb312nf_nri,            760
.equiv          nb312nf_nouter,         764
.equiv          nb312nf_ninner,         768
.equiv          nb312nf_salign,         772
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 776		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb312nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb312nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb312nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb312nf_nouter], eax
	mov [esp + nb312nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb312nf_half], eax
	mov [esp + nb312nf_half+4], ebx
	movsd xmm1, [esp + nb312nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb312nf_half], xmm1
	movapd [esp + nb312nf_three], xmm3

	mov eax, [ebp + nb312nf_p_tabscale]
	movsd xmm5, [eax]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb312nf_tsc],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb312nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb312nf_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb312nf_p_facel]
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
	movapd [esp + nb312nf_qqOO], xmm3
	movapd [esp + nb312nf_qqOH], xmm4
	movapd [esp + nb312nf_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb312nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb312nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb312nf_vdwparam]
	movlpd xmm0, [eax + edx*8] 
	movlpd xmm1, [eax + edx*8 + 8] 
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb312nf_c6], xmm0
	movapd [esp + nb312nf_c12], xmm1

.nb312nf_threadloop:
        mov   esi, [ebp + nb312nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb312nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb312nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb312nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb312nf_n], eax
        mov [esp + nb312nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb312nf_outerstart
        jmp .nb312nf_end

.nb312nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb312nf_nouter]
	mov [esp + nb312nf_nouter], ebx

.nb312nf_outer:
	mov   eax, [ebp + nb312nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb312nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb312nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb312nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb312nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb312nf_ixO], xmm3
	movapd [esp + nb312nf_iyO], xmm4
	movapd [esp + nb312nf_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 24]
	addsd xmm1, [eax + ebx*8 + 32]
	addsd xmm2, [eax + ebx*8 + 40]		
	addsd xmm3, [eax + ebx*8 + 48]
	addsd xmm4, [eax + ebx*8 + 56]
	addsd xmm5, [eax + ebx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb312nf_ixH1], xmm0
	movapd [esp + nb312nf_iyH1], xmm1
	movapd [esp + nb312nf_izH1], xmm2
	movapd [esp + nb312nf_ixH2], xmm3
	movapd [esp + nb312nf_iyH2], xmm4
	movapd [esp + nb312nf_izH2], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb312nf_vctot], xmm4
	movapd [esp + nb312nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb312nf_jindex]
	mov   ecx, [eax +esi*4]	     ;# jindex[n] 
	mov   edx, [eax +esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb312nf_pos]
	mov   eax, [ebp + nb312nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb312nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb312nf_ninner]
	mov   [esp + nb312nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb312nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb312nf_unroll_loop
	jmp   .nb312nf_checksingle
.nb312nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb312nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb312nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb312nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	
	;# move j coordinates to local temp variables 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movhpd xmm2, [esi + ebx*8]
	movhpd xmm3, [esi + ebx*8 + 8]
	movhpd xmm4, [esi + ebx*8 + 16]
	movhpd xmm5, [esi + ebx*8 + 24]
	movhpd xmm6, [esi + ebx*8 + 32]
	movhpd xmm7, [esi + ebx*8 + 40]
	movapd 	[esp + nb312nf_jxO], xmm2
	movapd 	[esp + nb312nf_jyO], xmm3
	movapd 	[esp + nb312nf_jzO], xmm4
	movapd 	[esp + nb312nf_jxH1], xmm5
	movapd 	[esp + nb312nf_jyH1], xmm6
	movapd 	[esp + nb312nf_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + ebx*8 + 56]
	movhpd xmm4, [esi + ebx*8 + 64]
	movapd 	[esp + nb312nf_jxH2], xmm2
	movapd 	[esp + nb312nf_jyH2], xmm3
	movapd 	[esp + nb312nf_jzH2], xmm4
	
	movapd xmm0, [esp + nb312nf_ixO]
	movapd xmm1, [esp + nb312nf_iyO]
	movapd xmm2, [esp + nb312nf_izO]
	movapd xmm3, [esp + nb312nf_ixO]
	movapd xmm4, [esp + nb312nf_iyO]
	movapd xmm5, [esp + nb312nf_izO]
	subpd  xmm0, [esp + nb312nf_jxO]
	subpd  xmm1, [esp + nb312nf_jyO]
	subpd  xmm2, [esp + nb312nf_jzO]
	subpd  xmm3, [esp + nb312nf_jxH1]
	subpd  xmm4, [esp + nb312nf_jyH1]
	subpd  xmm5, [esp + nb312nf_jzH1]
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
	movapd [esp + nb312nf_rsqOO], xmm0
	movapd [esp + nb312nf_rsqOH1], xmm3

	movapd xmm0, [esp + nb312nf_ixO]
	movapd xmm1, [esp + nb312nf_iyO]
	movapd xmm2, [esp + nb312nf_izO]
	movapd xmm3, [esp + nb312nf_ixH1]
	movapd xmm4, [esp + nb312nf_iyH1]
	movapd xmm5, [esp + nb312nf_izH1]
	subpd  xmm0, [esp + nb312nf_jxH2]
	subpd  xmm1, [esp + nb312nf_jyH2]
	subpd  xmm2, [esp + nb312nf_jzH2]
	subpd  xmm3, [esp + nb312nf_jxO]
	subpd  xmm4, [esp + nb312nf_jyO]
	subpd  xmm5, [esp + nb312nf_jzO]
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
	movapd [esp + nb312nf_rsqOH2], xmm0
	movapd [esp + nb312nf_rsqH1O], xmm3

	movapd xmm0, [esp + nb312nf_ixH1]
	movapd xmm1, [esp + nb312nf_iyH1]
	movapd xmm2, [esp + nb312nf_izH1]
	movapd xmm3, [esp + nb312nf_ixH1]
	movapd xmm4, [esp + nb312nf_iyH1]
	movapd xmm5, [esp + nb312nf_izH1]
	subpd  xmm0, [esp + nb312nf_jxH1]
	subpd  xmm1, [esp + nb312nf_jyH1]
	subpd  xmm2, [esp + nb312nf_jzH1]
	subpd  xmm3, [esp + nb312nf_jxH2]
	subpd  xmm4, [esp + nb312nf_jyH2]
	subpd  xmm5, [esp + nb312nf_jzH2]
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
	movapd [esp + nb312nf_rsqH1H1], xmm0
	movapd [esp + nb312nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb312nf_ixH2]
	movapd xmm1, [esp + nb312nf_iyH2]
	movapd xmm2, [esp + nb312nf_izH2]
	movapd xmm3, [esp + nb312nf_ixH2]
	movapd xmm4, [esp + nb312nf_iyH2]
	movapd xmm5, [esp + nb312nf_izH2]
	subpd  xmm0, [esp + nb312nf_jxO]
	subpd  xmm1, [esp + nb312nf_jyO]
	subpd  xmm2, [esp + nb312nf_jzO]
	subpd  xmm3, [esp + nb312nf_jxH1]
	subpd  xmm4, [esp + nb312nf_jyH1]
	subpd  xmm5, [esp + nb312nf_jzH1]
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
	movapd [esp + nb312nf_rsqH2O], xmm0
	movapd [esp + nb312nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb312nf_ixH2]
	movapd xmm1, [esp + nb312nf_iyH2]
	movapd xmm2, [esp + nb312nf_izH2]
	subpd  xmm0, [esp + nb312nf_jxH2]
	subpd  xmm1, [esp + nb312nf_jyH2]
	subpd  xmm2, [esp + nb312nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb312nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb312nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb312nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvH2H2], xmm1
	movapd [esp + nb312nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb312nf_rsqOO]
	movapd xmm4, [esp + nb312nf_rsqOH1]	
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
	movapd  xmm3, [esp + nb312nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb312nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb312nf_half] ;# rinv
	movapd [esp + nb312nf_rinvOO], xmm1
	movapd [esp + nb312nf_rinvOH1], xmm5

	movapd xmm0, [esp + nb312nf_rsqOH2]
	movapd xmm4, [esp + nb312nf_rsqH1O]	
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
	movapd  xmm3, [esp + nb312nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb312nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvOH2], xmm1
	movapd [esp + nb312nf_rinvH1O], xmm5

	movapd xmm0, [esp + nb312nf_rsqH1H1]
	movapd xmm4, [esp + nb312nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb312nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb312nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb312nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvH1H1], xmm1
	movapd [esp + nb312nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb312nf_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb312nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb312nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb312nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm0, [esp + nb312nf_rinvOO]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqOO] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqOO]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# start doing lj 
	movapd xmm2, xmm0
	mulpd  xmm2, xmm2
	movapd xmm1, xmm2
	mulpd  xmm1, xmm2
	mulpd  xmm1, xmm2	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm1, [esp + nb312nf_c6]
	mulpd  xmm2, [esp + nb312nf_c12]
	movapd xmm4, xmm2
	subpd  xmm4, xmm1
	addpd  xmm4, [esp + nb312nf_Vvdwtot]
	movapd [esp + nb312nf_Vvdwtot], xmm4

	;# O-H1 interaction 
	movapd xmm0, [esp + nb312nf_rinvOH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqOH1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# O-H2 interaction  
	movapd xmm0, [esp + nb312nf_rinvOH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqOH2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# H1-O interaction 
	movapd xmm0, [esp + nb312nf_rinvH1O]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqH1O] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb312nf_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb312nf_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# H2-O interaction 
	movapd xmm0, [esp + nb312nf_rinvH2O]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqH2O] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb312nf_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb312nf_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [esp + nb312nf_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [esp + nb312nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]
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
	movapd xmm3, [esp + nb312nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
	
    addpd  xmm5, [esp + nb312nf_vctot]
    movapd [esp + nb312nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb312nf_innerk],  2
	jl    .nb312nf_checksingle
	jmp   .nb312nf_unroll_loop
.nb312nf_checksingle:
	mov   edx, [esp + nb312nf_innerk]
	and   edx, 1
	jnz   .nb312nf_dosingle
	jmp   .nb312nf_updateouterdata
.nb312nf_dosingle:
	mov   edx, [esp + nb312nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]

	mov esi, [ebp + nb312nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movapd 	[esp + nb312nf_jxO], xmm2
	movapd 	[esp + nb312nf_jyO], xmm3
	movapd 	[esp + nb312nf_jzO], xmm4
	movapd 	[esp + nb312nf_jxH1], xmm5
	movapd 	[esp + nb312nf_jyH1], xmm6
	movapd 	[esp + nb312nf_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movapd 	[esp + nb312nf_jxH2], xmm2
	movapd 	[esp + nb312nf_jyH2], xmm3
	movapd 	[esp + nb312nf_jzH2], xmm4
	
	movapd xmm0, [esp + nb312nf_ixO]
	movapd xmm1, [esp + nb312nf_iyO]
	movapd xmm2, [esp + nb312nf_izO]
	movapd xmm3, [esp + nb312nf_ixO]
	movapd xmm4, [esp + nb312nf_iyO]
	movapd xmm5, [esp + nb312nf_izO]
	subsd  xmm0, [esp + nb312nf_jxO]
	subsd  xmm1, [esp + nb312nf_jyO]
	subsd  xmm2, [esp + nb312nf_jzO]
	subsd  xmm3, [esp + nb312nf_jxH1]
	subsd  xmm4, [esp + nb312nf_jyH1]
	subsd  xmm5, [esp + nb312nf_jzH1]
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
	movapd [esp + nb312nf_rsqOO], xmm0
	movapd [esp + nb312nf_rsqOH1], xmm3

	movapd xmm0, [esp + nb312nf_ixO]
	movapd xmm1, [esp + nb312nf_iyO]
	movapd xmm2, [esp + nb312nf_izO]
	movapd xmm3, [esp + nb312nf_ixH1]
	movapd xmm4, [esp + nb312nf_iyH1]
	movapd xmm5, [esp + nb312nf_izH1]
	subsd  xmm0, [esp + nb312nf_jxH2]
	subsd  xmm1, [esp + nb312nf_jyH2]
	subsd  xmm2, [esp + nb312nf_jzH2]
	subsd  xmm3, [esp + nb312nf_jxO]
	subsd  xmm4, [esp + nb312nf_jyO]
	subsd  xmm5, [esp + nb312nf_jzO]
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
	movapd [esp + nb312nf_rsqOH2], xmm0
	movapd [esp + nb312nf_rsqH1O], xmm3

	movapd xmm0, [esp + nb312nf_ixH1]
	movapd xmm1, [esp + nb312nf_iyH1]
	movapd xmm2, [esp + nb312nf_izH1]
	movapd xmm3, [esp + nb312nf_ixH1]
	movapd xmm4, [esp + nb312nf_iyH1]
	movapd xmm5, [esp + nb312nf_izH1]
	subsd  xmm0, [esp + nb312nf_jxH1]
	subsd  xmm1, [esp + nb312nf_jyH1]
	subsd  xmm2, [esp + nb312nf_jzH1]
	subsd  xmm3, [esp + nb312nf_jxH2]
	subsd  xmm4, [esp + nb312nf_jyH2]
	subsd  xmm5, [esp + nb312nf_jzH2]
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
	movapd [esp + nb312nf_rsqH1H1], xmm0
	movapd [esp + nb312nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb312nf_ixH2]
	movapd xmm1, [esp + nb312nf_iyH2]
	movapd xmm2, [esp + nb312nf_izH2]
	movapd xmm3, [esp + nb312nf_ixH2]
	movapd xmm4, [esp + nb312nf_iyH2]
	movapd xmm5, [esp + nb312nf_izH2]
	subsd  xmm0, [esp + nb312nf_jxO]
	subsd  xmm1, [esp + nb312nf_jyO]
	subsd  xmm2, [esp + nb312nf_jzO]
	subsd  xmm3, [esp + nb312nf_jxH1]
	subsd  xmm4, [esp + nb312nf_jyH1]
	subsd  xmm5, [esp + nb312nf_jzH1]
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
	movapd [esp + nb312nf_rsqH2O], xmm0
	movapd [esp + nb312nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb312nf_ixH2]
	movapd xmm1, [esp + nb312nf_iyH2]
	movapd xmm2, [esp + nb312nf_izH2]
	subsd  xmm0, [esp + nb312nf_jxH2]
	subsd  xmm1, [esp + nb312nf_jyH2]
	subsd  xmm2, [esp + nb312nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb312nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb312nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb312nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvH2H2], xmm1
	movapd [esp + nb312nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb312nf_rsqOO]
	movapd xmm4, [esp + nb312nf_rsqOH1]	
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
	movapd  xmm3, [esp + nb312nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb312nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb312nf_half] ;# rinv
	movapd [esp + nb312nf_rinvOO], xmm1
	movapd [esp + nb312nf_rinvOH1], xmm5

	movapd xmm0, [esp + nb312nf_rsqOH2]
	movapd xmm4, [esp + nb312nf_rsqH1O]	
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
	movapd  xmm3, [esp + nb312nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb312nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvOH2], xmm1
	movapd [esp + nb312nf_rinvH1O], xmm5

	movapd xmm0, [esp + nb312nf_rsqH1H1]
	movapd xmm4, [esp + nb312nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb312nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb312nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb312nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb312nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb312nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvH1H1], xmm1
	movapd [esp + nb312nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb312nf_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb312nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb312nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb312nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb312nf_half] ;# rinv 
	movapd [esp + nb312nf_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm0, [esp + nb312nf_rinvOO]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqOO] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqOO]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5

	;# start doing lj 
	movapd xmm2, xmm0
	mulsd  xmm2, xmm2
	movapd xmm1, xmm2
	mulsd  xmm1, xmm2
	mulsd  xmm1, xmm2	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm1, [esp + nb312nf_c6]
	mulsd  xmm2, [esp + nb312nf_c12]
	movapd xmm4, xmm2
	subsd  xmm4, xmm1
	addsd  xmm4, [esp + nb312nf_Vvdwtot]
	movlpd [esp + nb312nf_Vvdwtot], xmm4

	;# O-H1 interaction 
	movapd xmm0, [esp + nb312nf_rinvOH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqOH1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5

	;# O-H2 interaction  
	movapd xmm0, [esp + nb312nf_rinvOH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqOH2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5

	;# H1-O interaction 
	movapd xmm0, [esp + nb312nf_rinvH1O]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqH1O] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb312nf_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb312nf_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5

	;# H2-O interaction 
	movapd xmm0, [esp + nb312nf_rinvH2O]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqH2O] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	

	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb312nf_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5
	
	;# H2-H2 interaction 
	movapd xmm0, [esp + nb312nf_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [esp + nb312nf_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [esp + nb312nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb312nf_VFtab]

	movsd xmm4, [esi + eax*8]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 8]	;# F1 	
	movsd xmm6, [esi + eax*8 + 16]	;# G1 	
	movsd xmm7, [esi + eax*8 + 24]	;# H1 	
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb312nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [esp + nb312nf_vctot]
    movlpd [esp + nb312nf_vctot], xmm5
	
.nb312nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb312nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb312nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb312nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb312nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb312nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb312nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb312nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb312nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb312nf_n], esi
        jmp .nb312nf_outer
.nb312nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb312nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb312nf_end
        ;# non-zero, do one more workunit
        jmp   .nb312nf_threadloop
.nb312nf_end:
	emms

	mov eax, [esp + nb312nf_nouter]
	mov ebx, [esp + nb312nf_ninner]
	mov ecx, [ebp + nb312nf_outeriter]
	mov edx, [ebp + nb312nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb312nf_salign]
	add esp, eax
	add esp, 776
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
