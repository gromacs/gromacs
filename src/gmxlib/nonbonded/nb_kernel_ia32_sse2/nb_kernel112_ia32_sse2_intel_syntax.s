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



	
.globl nb_kernel112_ia32_sse2
.globl _nb_kernel112_ia32_sse2
nb_kernel112_ia32_sse2:	
_nb_kernel112_ia32_sse2:	
.equiv          nb112_p_nri,            8
.equiv          nb112_iinr,             12
.equiv          nb112_jindex,           16
.equiv          nb112_jjnr,             20
.equiv          nb112_shift,            24
.equiv          nb112_shiftvec,         28
.equiv          nb112_fshift,           32
.equiv          nb112_gid,              36
.equiv          nb112_pos,              40
.equiv          nb112_faction,          44
.equiv          nb112_charge,           48
.equiv          nb112_p_facel,          52
.equiv          nb112_argkrf,           56
.equiv          nb112_argcrf,           60
.equiv          nb112_Vc,               64
.equiv          nb112_type,             68
.equiv          nb112_p_ntype,          72
.equiv          nb112_vdwparam,         76
.equiv          nb112_Vvdw,             80
.equiv          nb112_p_tabscale,       84
.equiv          nb112_VFtab,            88
.equiv          nb112_invsqrta,         92
.equiv          nb112_dvda,             96
.equiv          nb112_p_gbtabscale,     100
.equiv          nb112_GBtab,            104
.equiv          nb112_p_nthreads,       108
.equiv          nb112_count,            112
.equiv          nb112_mtx,              116
.equiv          nb112_outeriter,        120
.equiv          nb112_inneriter,        124
.equiv          nb112_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb112_ixO,              0
.equiv          nb112_iyO,              16
.equiv          nb112_izO,              32
.equiv          nb112_ixH1,             48
.equiv          nb112_iyH1,             64
.equiv          nb112_izH1,             80
.equiv          nb112_ixH2,             96
.equiv          nb112_iyH2,             112
.equiv          nb112_izH2,             128
.equiv          nb112_jxO,              144
.equiv          nb112_jyO,              160
.equiv          nb112_jzO,              176
.equiv          nb112_jxH1,             192
.equiv          nb112_jyH1,             208
.equiv          nb112_jzH1,             224
.equiv          nb112_jxH2,             240
.equiv          nb112_jyH2,             256
.equiv          nb112_jzH2,             272
.equiv          nb112_dxOO,             288
.equiv          nb112_dyOO,             304
.equiv          nb112_dzOO,             320
.equiv          nb112_dxOH1,            336
.equiv          nb112_dyOH1,            352
.equiv          nb112_dzOH1,            368
.equiv          nb112_dxOH2,            384
.equiv          nb112_dyOH2,            400
.equiv          nb112_dzOH2,            416
.equiv          nb112_dxH1O,            432
.equiv          nb112_dyH1O,            448
.equiv          nb112_dzH1O,            464
.equiv          nb112_dxH1H1,           480
.equiv          nb112_dyH1H1,           496
.equiv          nb112_dzH1H1,           512
.equiv          nb112_dxH1H2,           528
.equiv          nb112_dyH1H2,           544
.equiv          nb112_dzH1H2,           560
.equiv          nb112_dxH2O,            576
.equiv          nb112_dyH2O,            592
.equiv          nb112_dzH2O,            608
.equiv          nb112_dxH2H1,           624
.equiv          nb112_dyH2H1,           640
.equiv          nb112_dzH2H1,           656
.equiv          nb112_dxH2H2,           672
.equiv          nb112_dyH2H2,           688
.equiv          nb112_dzH2H2,           704
.equiv          nb112_qqOO,             720
.equiv          nb112_qqOH,             736
.equiv          nb112_qqHH,             752
.equiv          nb112_c6,               768
.equiv          nb112_c12,              784
.equiv          nb112_six,              800
.equiv          nb112_twelve,           816
.equiv          nb112_vctot,            832
.equiv          nb112_Vvdwtot,          848
.equiv          nb112_fixO,             864
.equiv          nb112_fiyO,             880
.equiv          nb112_fizO,             896
.equiv          nb112_fixH1,            912
.equiv          nb112_fiyH1,            928
.equiv          nb112_fizH1,            944
.equiv          nb112_fixH2,            960
.equiv          nb112_fiyH2,            976
.equiv          nb112_fizH2,            992
.equiv          nb112_fjxO,             1008
.equiv          nb112_fjyO,             1024
.equiv          nb112_fjzO,             1040
.equiv          nb112_fjxH1,            1056
.equiv          nb112_fjyH1,            1072
.equiv          nb112_fjzH1,            1088
.equiv          nb112_fjxH2,            1104
.equiv          nb112_fjyH2,            1120
.equiv          nb112_fjzH2,            1136
.equiv          nb112_half,             1152
.equiv          nb112_three,            1168
.equiv          nb112_rsqOO,            1184
.equiv          nb112_rsqOH1,           1200
.equiv          nb112_rsqOH2,           1216
.equiv          nb112_rsqH1O,           1232
.equiv          nb112_rsqH1H1,          1248
.equiv          nb112_rsqH1H2,          1264
.equiv          nb112_rsqH2O,           1280
.equiv          nb112_rsqH2H1,          1296
.equiv          nb112_rsqH2H2,          1312
.equiv          nb112_rinvOO,           1328
.equiv          nb112_rinvOH1,          1344
.equiv          nb112_rinvOH2,          1360
.equiv          nb112_rinvH1O,          1376
.equiv          nb112_rinvH1H1,         1392
.equiv          nb112_rinvH1H2,         1408
.equiv          nb112_rinvH2O,          1424
.equiv          nb112_rinvH2H1,         1440
.equiv          nb112_rinvH2H2,         1456
.equiv          nb112_is3,              1472
.equiv          nb112_ii3,              1476
.equiv          nb112_innerjjnr,        1480
.equiv          nb112_innerk,           1484
.equiv          nb112_n,                1488
.equiv          nb112_nn1,              1492
.equiv          nb112_nri,              1496
.equiv          nb112_nouter,           1500
.equiv          nb112_ninner,           1504
.equiv          nb112_salign,           1508
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1512		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb112_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb112_p_nri]
	mov ecx, [ecx]
	mov [esp + nb112_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb112_nouter], eax
	mov [esp + nb112_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb112_half], eax
	mov [esp + nb112_half+4], ebx
	movsd xmm1, [esp + nb112_half]
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
	movapd [esp + nb112_half], xmm1
	movapd [esp + nb112_three], xmm3
	movapd [esp + nb112_six], xmm4
	movapd [esp + nb112_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb112_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb112_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb112_p_facel]
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
	movapd [esp + nb112_qqOO], xmm3
	movapd [esp + nb112_qqOH], xmm4
	movapd [esp + nb112_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb112_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb112_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb112_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movhpd xmm0, [eax + edx*8 + 8]
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0	
	movapd [esp + nb112_c6], xmm0
	movapd [esp + nb112_c12], xmm1

.nb112_threadloop:
        mov   esi, [ebp + nb112_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb112_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb112_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb112_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb112_n], eax
        mov [esp + nb112_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb112_outerstart
        jmp .nb112_end

.nb112_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb112_nouter]
	mov [esp + nb112_nouter], ebx

.nb112_outer:
	mov   eax, [ebp + nb112_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb112_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb112_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb112_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb112_pos]    ;# eax = base of pos[]  
	mov   [esp + nb112_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb112_ixO], xmm3
	movapd [esp + nb112_iyO], xmm4
	movapd [esp + nb112_izO], xmm5

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
	movapd [esp + nb112_ixH1], xmm0
	movapd [esp + nb112_iyH1], xmm1
	movapd [esp + nb112_izH1], xmm2
	movapd [esp + nb112_ixH2], xmm3
	movapd [esp + nb112_iyH2], xmm4
	movapd [esp + nb112_izH2], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb112_vctot], xmm4
	movapd [esp + nb112_Vvdwtot], xmm4
	movapd [esp + nb112_fixO], xmm4
	movapd [esp + nb112_fiyO], xmm4
	movapd [esp + nb112_fizO], xmm4
	movapd [esp + nb112_fixH1], xmm4
	movapd [esp + nb112_fiyH1], xmm4
	movapd [esp + nb112_fizH1], xmm4
	movapd [esp + nb112_fixH2], xmm4
	movapd [esp + nb112_fiyH2], xmm4
	movapd [esp + nb112_fizH2], xmm4
	
	mov   eax, [ebp + nb112_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb112_pos]
	mov   edi, [ebp + nb112_faction]	
	mov   eax, [ebp + nb112_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb112_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb112_ninner]
	mov   [esp + nb112_ninner], ecx
	add   edx, 0
	mov   [esp + nb112_innerk], edx    ;# number of innerloop atoms 
	jge  .nb112_unroll_loop
	jmp  .nb112_checksingle
.nb112_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb112_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb112_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb112_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb112_jxO], xmm2
	movapd 	[esp + nb112_jyO], xmm3
	movapd 	[esp + nb112_jzO], xmm4
	movapd 	[esp + nb112_jxH1], xmm5
	movapd 	[esp + nb112_jyH1], xmm6
	movapd 	[esp + nb112_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + ebx*8 + 56]
	movhpd xmm4, [esi + ebx*8 + 64]
	movapd 	[esp + nb112_jxH2], xmm2
	movapd 	[esp + nb112_jyH2], xmm3
	movapd 	[esp + nb112_jzH2], xmm4
	
	movapd xmm0, [esp + nb112_ixO]
	movapd xmm1, [esp + nb112_iyO]
	movapd xmm2, [esp + nb112_izO]
	movapd xmm3, [esp + nb112_ixO]
	movapd xmm4, [esp + nb112_iyO]
	movapd xmm5, [esp + nb112_izO]
	subpd  xmm0, [esp + nb112_jxO]
	subpd  xmm1, [esp + nb112_jyO]
	subpd  xmm2, [esp + nb112_jzO]
	subpd  xmm3, [esp + nb112_jxH1]
	subpd  xmm4, [esp + nb112_jyH1]
	subpd  xmm5, [esp + nb112_jzH1]
	movapd [esp + nb112_dxOO], xmm0
	movapd [esp + nb112_dyOO], xmm1
	movapd [esp + nb112_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb112_dxOH1], xmm3
	movapd [esp + nb112_dyOH1], xmm4
	movapd [esp + nb112_dzOH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb112_rsqOO], xmm0
	movapd [esp + nb112_rsqOH1], xmm3

	movapd xmm0, [esp + nb112_ixO]
	movapd xmm1, [esp + nb112_iyO]
	movapd xmm2, [esp + nb112_izO]
	movapd xmm3, [esp + nb112_ixH1]
	movapd xmm4, [esp + nb112_iyH1]
	movapd xmm5, [esp + nb112_izH1]
	subpd  xmm0, [esp + nb112_jxH2]
	subpd  xmm1, [esp + nb112_jyH2]
	subpd  xmm2, [esp + nb112_jzH2]
	subpd  xmm3, [esp + nb112_jxO]
	subpd  xmm4, [esp + nb112_jyO]
	subpd  xmm5, [esp + nb112_jzO]
	movapd [esp + nb112_dxOH2], xmm0
	movapd [esp + nb112_dyOH2], xmm1
	movapd [esp + nb112_dzOH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb112_dxH1O], xmm3
	movapd [esp + nb112_dyH1O], xmm4
	movapd [esp + nb112_dzH1O], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb112_rsqOH2], xmm0
	movapd [esp + nb112_rsqH1O], xmm3

	movapd xmm0, [esp + nb112_ixH1]
	movapd xmm1, [esp + nb112_iyH1]
	movapd xmm2, [esp + nb112_izH1]
	movapd xmm3, [esp + nb112_ixH1]
	movapd xmm4, [esp + nb112_iyH1]
	movapd xmm5, [esp + nb112_izH1]
	subpd  xmm0, [esp + nb112_jxH1]
	subpd  xmm1, [esp + nb112_jyH1]
	subpd  xmm2, [esp + nb112_jzH1]
	subpd  xmm3, [esp + nb112_jxH2]
	subpd  xmm4, [esp + nb112_jyH2]
	subpd  xmm5, [esp + nb112_jzH2]
	movapd [esp + nb112_dxH1H1], xmm0
	movapd [esp + nb112_dyH1H1], xmm1
	movapd [esp + nb112_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb112_dxH1H2], xmm3
	movapd [esp + nb112_dyH1H2], xmm4
	movapd [esp + nb112_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb112_rsqH1H1], xmm0
	movapd [esp + nb112_rsqH1H2], xmm3

	movapd xmm0, [esp + nb112_ixH2]
	movapd xmm1, [esp + nb112_iyH2]
	movapd xmm2, [esp + nb112_izH2]
	movapd xmm3, [esp + nb112_ixH2]
	movapd xmm4, [esp + nb112_iyH2]
	movapd xmm5, [esp + nb112_izH2]
	subpd  xmm0, [esp + nb112_jxO]
	subpd  xmm1, [esp + nb112_jyO]
	subpd  xmm2, [esp + nb112_jzO]
	subpd  xmm3, [esp + nb112_jxH1]
	subpd  xmm4, [esp + nb112_jyH1]
	subpd  xmm5, [esp + nb112_jzH1]
	movapd [esp + nb112_dxH2O], xmm0
	movapd [esp + nb112_dyH2O], xmm1
	movapd [esp + nb112_dzH2O], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb112_dxH2H1], xmm3
	movapd [esp + nb112_dyH2H1], xmm4
	movapd [esp + nb112_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb112_rsqH2O], xmm0
	movapd [esp + nb112_rsqH2H1], xmm4

	movapd xmm0, [esp + nb112_ixH2]
	movapd xmm1, [esp + nb112_iyH2]
	movapd xmm2, [esp + nb112_izH2]
	subpd  xmm0, [esp + nb112_jxH2]
	subpd  xmm1, [esp + nb112_jyH2]
	subpd  xmm2, [esp + nb112_jzH2]
	movapd [esp + nb112_dxH2H2], xmm0
	movapd [esp + nb112_dyH2H2], xmm1
	movapd [esp + nb112_dzH2H2], xmm2
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb112_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb112_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112_half] ;# iter1 
	mulpd   xmm7, [esp + nb112_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112_half] ;# rinv 
	mulpd   xmm5, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvH2H2], xmm1
	movapd [esp + nb112_rinvH2H1], xmm5

	movapd xmm0, [esp + nb112_rsqOO]
	movapd xmm4, [esp + nb112_rsqOH1]	
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
	movapd  xmm3, [esp + nb112_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb112_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112_half] ;# rinv 
	mulpd   xmm5, [esp + nb112_half] ;# rinv
	movapd [esp + nb112_rinvOO], xmm1
	movapd [esp + nb112_rinvOH1], xmm5

	movapd xmm0, [esp + nb112_rsqOH2]
	movapd xmm4, [esp + nb112_rsqH1O]	
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
	movapd  xmm3, [esp + nb112_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112_half] ;# iter1 
	mulpd   xmm7, [esp + nb112_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112_half] ;# rinv 
	mulpd   xmm5, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvOH2], xmm1
	movapd [esp + nb112_rinvH1O], xmm5

	movapd xmm0, [esp + nb112_rsqH1H1]
	movapd xmm4, [esp + nb112_rsqH1H2]	
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
	movapd  xmm3, [esp + nb112_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112_half] ;# iter1a 
	mulpd   xmm7, [esp + nb112_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112_half] ;# rinv 
	mulpd   xmm5, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvH1H1], xmm1
	movapd [esp + nb112_rinvH1H2], xmm5

	movapd xmm0, [esp + nb112_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb112_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb112_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb112_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvH2O], xmm1

	;# start with OO interaction 
	movapd xmm0, [esp + nb112_rinvOO]
	movapd xmm7, xmm0
	mulpd  xmm0, xmm0
	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	mulpd  xmm7, [esp + nb112_qqOO]
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm1, [esp + nb112_c6]	
	mulpd  xmm2, [esp + nb112_c12]	
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addpd  xmm3, [esp + nb112_Vvdwtot]
	mulpd  xmm1, [esp + nb112_six]
	mulpd  xmm2, [esp + nb112_twelve]
	movapd [esp + nb112_Vvdwtot], xmm3
	subpd  xmm2, xmm1
	addpd  xmm2, xmm7
	addpd  xmm7, [esp + nb112_vctot]
	mulpd  xmm0, xmm2	
 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb112_dxOO]
	mulpd xmm1, [esp + nb112_dyOO]
	mulpd xmm2, [esp + nb112_dzOO]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixO]
	addpd xmm1, [esp + nb112_fiyO]
	addpd xmm2, [esp + nb112_fizO]
	movapd [esp + nb112_fjxO], xmm3
	movapd [esp + nb112_fjyO], xmm4
	movapd [esp + nb112_fjzO], xmm5
	movapd [esp + nb112_fixO], xmm0
	movapd [esp + nb112_fiyO], xmm1
	movapd [esp + nb112_fizO], xmm2

	;# O-H1 interaction 
	movapd xmm0, [esp + nb112_rinvOH1]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqOH]
	mulpd xmm0, xmm1	;# fsOH1  
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb112_dxOH1]
	mulpd xmm1, [esp + nb112_dyOH1]
	mulpd xmm2, [esp + nb112_dzOH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixO]
	addpd xmm1, [esp + nb112_fiyO]
	addpd xmm2, [esp + nb112_fizO]
	movapd [esp + nb112_fjxH1], xmm3
	movapd [esp + nb112_fjyH1], xmm4
	movapd [esp + nb112_fjzH1], xmm5
	movapd [esp + nb112_fixO], xmm0
	movapd [esp + nb112_fiyO], xmm1
	movapd [esp + nb112_fizO], xmm2

	;# O-H2 interaction  
	movapd xmm0, [esp + nb112_rinvOH2]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqOH]
	mulpd xmm0, xmm1	;# fsOH2  
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb112_dxOH2]
	mulpd xmm1, [esp + nb112_dyOH2]
	mulpd xmm2, [esp + nb112_dzOH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixO]
	addpd xmm1, [esp + nb112_fiyO]
	addpd xmm2, [esp + nb112_fizO]
	movapd [esp + nb112_fjxH2], xmm3
	movapd [esp + nb112_fjyH2], xmm4
	movapd [esp + nb112_fjzH2], xmm5
	movapd [esp + nb112_fixO], xmm0
	movapd [esp + nb112_fiyO], xmm1
	movapd [esp + nb112_fizO], xmm2

	;# H1-O interaction 
	movapd xmm0, [esp + nb112_rinvH1O]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqOH]
	mulpd xmm0, xmm1	;# fsH1O 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxO]
	movapd xmm4, [esp + nb112_fjyO]
	movapd xmm5, [esp + nb112_fjzO]
	mulpd xmm0, [esp + nb112_dxH1O]
	mulpd xmm1, [esp + nb112_dyH1O]
	mulpd xmm2, [esp + nb112_dzH1O]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixH1]
	addpd xmm1, [esp + nb112_fiyH1]
	addpd xmm2, [esp + nb112_fizH1]
	movapd [esp + nb112_fjxO], xmm3
	movapd [esp + nb112_fjyO], xmm4
	movapd [esp + nb112_fjzO], xmm5
	movapd [esp + nb112_fixH1], xmm0
	movapd [esp + nb112_fiyH1], xmm1
	movapd [esp + nb112_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb112_rinvH1H1]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqHH]
	mulpd xmm0, xmm1	;# fsH1H1 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxH1]
	movapd xmm4, [esp + nb112_fjyH1]
	movapd xmm5, [esp + nb112_fjzH1]
	mulpd xmm0, [esp + nb112_dxH1H1]
	mulpd xmm1, [esp + nb112_dyH1H1]
	mulpd xmm2, [esp + nb112_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixH1]
	addpd xmm1, [esp + nb112_fiyH1]
	addpd xmm2, [esp + nb112_fizH1]
	movapd [esp + nb112_fjxH1], xmm3
	movapd [esp + nb112_fjyH1], xmm4
	movapd [esp + nb112_fjzH1], xmm5
	movapd [esp + nb112_fixH1], xmm0
	movapd [esp + nb112_fiyH1], xmm1
	movapd [esp + nb112_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb112_rinvH1H2]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqHH]
	mulpd xmm0, xmm1	;# fsOH2  
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxH2]
	movapd xmm4, [esp + nb112_fjyH2]
	movapd xmm5, [esp + nb112_fjzH2]
	mulpd xmm0, [esp + nb112_dxH1H2]
	mulpd xmm1, [esp + nb112_dyH1H2]
	mulpd xmm2, [esp + nb112_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixH1]
	addpd xmm1, [esp + nb112_fiyH1]
	addpd xmm2, [esp + nb112_fizH1]
	movapd [esp + nb112_fjxH2], xmm3
	movapd [esp + nb112_fjyH2], xmm4
	movapd [esp + nb112_fjzH2], xmm5
	movapd [esp + nb112_fixH1], xmm0
	movapd [esp + nb112_fiyH1], xmm1
	movapd [esp + nb112_fizH1], xmm2

	;# H2-O interaction 
	movapd xmm0, [esp + nb112_rinvH2O]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqOH]
	mulpd xmm0, xmm1	;# fsH2O 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxO]
	movapd xmm4, [esp + nb112_fjyO]
	movapd xmm5, [esp + nb112_fjzO]
	mulpd xmm0, [esp + nb112_dxH2O]
	mulpd xmm1, [esp + nb112_dyH2O]
	mulpd xmm2, [esp + nb112_dzH2O]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixH2]
	addpd xmm1, [esp + nb112_fiyH2]
	addpd xmm2, [esp + nb112_fizH2]
	movapd [esp + nb112_fjxO], xmm3
	movapd [esp + nb112_fjyO], xmm4
	movapd [esp + nb112_fjzO], xmm5
	movapd [esp + nb112_fixH2], xmm0
	movapd [esp + nb112_fiyH2], xmm1
	movapd [esp + nb112_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb112_rinvH2H1]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqHH]
	mulpd xmm0, xmm1	;# fsH2H1 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxH1]
	movapd xmm4, [esp + nb112_fjyH1]
	movapd xmm5, [esp + nb112_fjzH1]
	mulpd xmm0, [esp + nb112_dxH2H1]
	mulpd xmm1, [esp + nb112_dyH2H1]
	mulpd xmm2, [esp + nb112_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixH2]
	addpd xmm1, [esp + nb112_fiyH2]
	addpd xmm2, [esp + nb112_fizH2]
	movapd [esp + nb112_fjxH1], xmm3
	movapd [esp + nb112_fjyH1], xmm4
	movapd [esp + nb112_fjzH1], xmm5
	movapd [esp + nb112_fixH2], xmm0
	movapd [esp + nb112_fiyH2], xmm1
	movapd [esp + nb112_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb112_rinvH2H2]
	movapd xmm1, xmm0
	mulpd xmm0, xmm0
	mulpd xmm1, [esp + nb112_qqHH]
	mulpd xmm0, xmm1	;# fsH2H2 
	addpd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd [esp + nb112_vctot], xmm7
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxH2]
	movapd xmm4, [esp + nb112_fjyH2]
	movapd xmm5, [esp + nb112_fjzH2]
	mulpd xmm0, [esp + nb112_dxH2H2]
	mulpd xmm1, [esp + nb112_dyH2H2]
	mulpd xmm2, [esp + nb112_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb112_fixH2]
	addpd xmm1, [esp + nb112_fiyH2]
	addpd xmm2, [esp + nb112_fizH2]
	movapd [esp + nb112_fjxH2], xmm3
	movapd [esp + nb112_fjyH2], xmm4
	movapd [esp + nb112_fjzH2], xmm5
	movapd [esp + nb112_fixH2], xmm0
	movapd [esp + nb112_fiyH2], xmm1
	movapd [esp + nb112_fizH2], xmm2

	mov edi, [ebp + nb112_faction]
		
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
	addpd xmm0, [esp + nb112_fjxO]
	addpd xmm1, [esp + nb112_fjyO]
	addpd xmm2, [esp + nb112_fjzO]
	addpd xmm3, [esp + nb112_fjxH1]
	addpd xmm4, [esp + nb112_fjyH1]
	addpd xmm5, [esp + nb112_fjzH1]
	addpd xmm6, [esp + nb112_fjxH2]
	addpd xmm7, [esp + nb112_fjyH2]
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
	addpd xmm0, [esp + nb112_fjzH2]
	movlpd [edi + eax*8 + 64], xmm0
	movhpd [edi + ebx*8 + 64], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb112_innerk],  2
	jl   .nb112_checksingle
	jmp  .nb112_unroll_loop
.nb112_checksingle:
	mov  edx, [esp + nb112_innerk]
	and  edx, 1
	jnz  .nb112_dosingle
	jmp  .nb112_updateouterdata
.nb112_dosingle:
	mov   edx, [esp + nb112_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb112_innerjjnr],  4	

	mov esi, [ebp + nb112_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movapd 	[esp + nb112_jxO], xmm2
	movapd 	[esp + nb112_jyO], xmm3
	movapd 	[esp + nb112_jzO], xmm4
	movapd 	[esp + nb112_jxH1], xmm5
	movapd 	[esp + nb112_jyH1], xmm6
	movapd 	[esp + nb112_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movapd 	[esp + nb112_jxH2], xmm2
	movapd 	[esp + nb112_jyH2], xmm3
	movapd 	[esp + nb112_jzH2], xmm4
	
	movapd xmm0, [esp + nb112_ixO]
	movapd xmm1, [esp + nb112_iyO]
	movapd xmm2, [esp + nb112_izO]
	movapd xmm3, [esp + nb112_ixO]
	movapd xmm4, [esp + nb112_iyO]
	movapd xmm5, [esp + nb112_izO]
	subsd  xmm0, [esp + nb112_jxO]
	subsd  xmm1, [esp + nb112_jyO]
	subsd  xmm2, [esp + nb112_jzO]
	subsd  xmm3, [esp + nb112_jxH1]
	subsd  xmm4, [esp + nb112_jyH1]
	subsd  xmm5, [esp + nb112_jzH1]
	movapd [esp + nb112_dxOO], xmm0
	movapd [esp + nb112_dyOO], xmm1
	movapd [esp + nb112_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb112_dxOH1], xmm3
	movapd [esp + nb112_dyOH1], xmm4
	movapd [esp + nb112_dzOH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb112_rsqOO], xmm0
	movapd [esp + nb112_rsqOH1], xmm3

	movapd xmm0, [esp + nb112_ixO]
	movapd xmm1, [esp + nb112_iyO]
	movapd xmm2, [esp + nb112_izO]
	movapd xmm3, [esp + nb112_ixH1]
	movapd xmm4, [esp + nb112_iyH1]
	movapd xmm5, [esp + nb112_izH1]
	subsd  xmm0, [esp + nb112_jxH2]
	subsd  xmm1, [esp + nb112_jyH2]
	subsd  xmm2, [esp + nb112_jzH2]
	subsd  xmm3, [esp + nb112_jxO]
	subsd  xmm4, [esp + nb112_jyO]
	subsd  xmm5, [esp + nb112_jzO]
	movapd [esp + nb112_dxOH2], xmm0
	movapd [esp + nb112_dyOH2], xmm1
	movapd [esp + nb112_dzOH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb112_dxH1O], xmm3
	movapd [esp + nb112_dyH1O], xmm4
	movapd [esp + nb112_dzH1O], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb112_rsqOH2], xmm0
	movapd [esp + nb112_rsqH1O], xmm3

	movapd xmm0, [esp + nb112_ixH1]
	movapd xmm1, [esp + nb112_iyH1]
	movapd xmm2, [esp + nb112_izH1]
	movapd xmm3, [esp + nb112_ixH1]
	movapd xmm4, [esp + nb112_iyH1]
	movapd xmm5, [esp + nb112_izH1]
	subsd  xmm0, [esp + nb112_jxH1]
	subsd  xmm1, [esp + nb112_jyH1]
	subsd  xmm2, [esp + nb112_jzH1]
	subsd  xmm3, [esp + nb112_jxH2]
	subsd  xmm4, [esp + nb112_jyH2]
	subsd  xmm5, [esp + nb112_jzH2]
	movapd [esp + nb112_dxH1H1], xmm0
	movapd [esp + nb112_dyH1H1], xmm1
	movapd [esp + nb112_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb112_dxH1H2], xmm3
	movapd [esp + nb112_dyH1H2], xmm4
	movapd [esp + nb112_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb112_rsqH1H1], xmm0
	movapd [esp + nb112_rsqH1H2], xmm3

	movapd xmm0, [esp + nb112_ixH2]
	movapd xmm1, [esp + nb112_iyH2]
	movapd xmm2, [esp + nb112_izH2]
	movapd xmm3, [esp + nb112_ixH2]
	movapd xmm4, [esp + nb112_iyH2]
	movapd xmm5, [esp + nb112_izH2]
	subsd  xmm0, [esp + nb112_jxO]
	subsd  xmm1, [esp + nb112_jyO]
	subsd  xmm2, [esp + nb112_jzO]
	subsd  xmm3, [esp + nb112_jxH1]
	subsd  xmm4, [esp + nb112_jyH1]
	subsd  xmm5, [esp + nb112_jzH1]
	movapd [esp + nb112_dxH2O], xmm0
	movapd [esp + nb112_dyH2O], xmm1
	movapd [esp + nb112_dzH2O], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb112_dxH2H1], xmm3
	movapd [esp + nb112_dyH2H1], xmm4
	movapd [esp + nb112_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb112_rsqH2O], xmm0
	movapd [esp + nb112_rsqH2H1], xmm4

	movapd xmm0, [esp + nb112_ixH2]
	movapd xmm1, [esp + nb112_iyH2]
	movapd xmm2, [esp + nb112_izH2]
	subsd  xmm0, [esp + nb112_jxH2]
	subsd  xmm1, [esp + nb112_jyH2]
	subsd  xmm2, [esp + nb112_jzH2]
	movapd [esp + nb112_dxH2H2], xmm0
	movapd [esp + nb112_dyH2H2], xmm1
	movapd [esp + nb112_dzH2H2], xmm2
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb112_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb112_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112_half] ;# iter1 
	mulsd   xmm7, [esp + nb112_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112_half] ;# rinv 
	mulsd   xmm5, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvH2H2], xmm1
	movapd [esp + nb112_rinvH2H1], xmm5

	movapd xmm0, [esp + nb112_rsqOO]
	movapd xmm4, [esp + nb112_rsqOH1]	
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
	movapd  xmm3, [esp + nb112_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb112_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112_half] ;# rinv 
	mulsd   xmm5, [esp + nb112_half] ;# rinv
	movapd [esp + nb112_rinvOO], xmm1
	movapd [esp + nb112_rinvOH1], xmm5

	movapd xmm0, [esp + nb112_rsqOH2]
	movapd xmm4, [esp + nb112_rsqH1O]	
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
	movapd  xmm3, [esp + nb112_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112_half] ;# iter1 
	mulsd   xmm7, [esp + nb112_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112_half] ;# rinv 
	mulsd   xmm5, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvOH2], xmm1
	movapd [esp + nb112_rinvH1O], xmm5

	movapd xmm0, [esp + nb112_rsqH1H1]
	movapd xmm4, [esp + nb112_rsqH1H2]	
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
	movapd  xmm3, [esp + nb112_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112_half] ;# iter1a 
	mulsd   xmm7, [esp + nb112_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112_half] ;# rinv 
	mulsd   xmm5, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvH1H1], xmm1
	movapd [esp + nb112_rinvH1H2], xmm5

	movapd xmm0, [esp + nb112_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb112_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb112_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb112_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb112_half] ;# rinv 
	movapd [esp + nb112_rinvH2O], xmm1

	;# start with OO interaction 
	movapd xmm0, [esp + nb112_rinvOO]
	movapd xmm7, xmm0
	mulsd  xmm0, xmm0
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	mulsd  xmm7, [esp + nb112_qqOO]
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm1, [esp + nb112_c6]	
	mulsd  xmm2, [esp + nb112_c12]	
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addsd  xmm3, [esp + nb112_Vvdwtot]
	mulsd  xmm1, [esp + nb112_six]
	mulsd  xmm2, [esp + nb112_twelve]
	movlpd [esp + nb112_Vvdwtot], xmm3
	subsd  xmm2, xmm1
	addsd  xmm2, xmm7
	addsd  xmm7, [esp + nb112_vctot]
	mulsd  xmm0, xmm2	
 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb112_dxOO]
	mulsd xmm1, [esp + nb112_dyOO]
	mulsd xmm2, [esp + nb112_dzOO]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixO]
	addsd xmm1, [esp + nb112_fiyO]
	addsd xmm2, [esp + nb112_fizO]
	movlpd [esp + nb112_fjxO], xmm3
	movlpd [esp + nb112_fjyO], xmm4
	movlpd [esp + nb112_fjzO], xmm5
	movlpd [esp + nb112_fixO], xmm0
	movlpd [esp + nb112_fiyO], xmm1
	movlpd [esp + nb112_fizO], xmm2

	;# O-H1 interaction 
	movapd xmm0, [esp + nb112_rinvOH1]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqOH]
	mulsd xmm0, xmm1	;# fsOH1  
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb112_dxOH1]
	mulsd xmm1, [esp + nb112_dyOH1]
	mulsd xmm2, [esp + nb112_dzOH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixO]
	addsd xmm1, [esp + nb112_fiyO]
	addsd xmm2, [esp + nb112_fizO]
	movlpd [esp + nb112_fjxH1], xmm3
	movlpd [esp + nb112_fjyH1], xmm4
	movlpd [esp + nb112_fjzH1], xmm5
	movlpd [esp + nb112_fixO], xmm0
	movlpd [esp + nb112_fiyO], xmm1
	movlpd [esp + nb112_fizO], xmm2

	;# O-H2 interaction  
	movapd xmm0, [esp + nb112_rinvOH2]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqOH]
	mulsd xmm0, xmm1	;# fsOH2  
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb112_dxOH2]
	mulsd xmm1, [esp + nb112_dyOH2]
	mulsd xmm2, [esp + nb112_dzOH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixO]
	addsd xmm1, [esp + nb112_fiyO]
	addsd xmm2, [esp + nb112_fizO]
	movlpd [esp + nb112_fjxH2], xmm3
	movlpd [esp + nb112_fjyH2], xmm4
	movlpd [esp + nb112_fjzH2], xmm5
	movlpd [esp + nb112_fixO], xmm0
	movlpd [esp + nb112_fiyO], xmm1
	movlpd [esp + nb112_fizO], xmm2

	;# H1-O interaction 
	movapd xmm0, [esp + nb112_rinvH1O]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqOH]
	mulsd xmm0, xmm1	;# fsH1O 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxO]
	movapd xmm4, [esp + nb112_fjyO]
	movapd xmm5, [esp + nb112_fjzO]
	mulsd xmm0, [esp + nb112_dxH1O]
	mulsd xmm1, [esp + nb112_dyH1O]
	mulsd xmm2, [esp + nb112_dzH1O]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixH1]
	addsd xmm1, [esp + nb112_fiyH1]
	addsd xmm2, [esp + nb112_fizH1]
	movlpd [esp + nb112_fjxO], xmm3
	movlpd [esp + nb112_fjyO], xmm4
	movlpd [esp + nb112_fjzO], xmm5
	movlpd [esp + nb112_fixH1], xmm0
	movlpd [esp + nb112_fiyH1], xmm1
	movlpd [esp + nb112_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb112_rinvH1H1]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqHH]
	mulsd xmm0, xmm1	;# fsH1H1 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxH1]
	movapd xmm4, [esp + nb112_fjyH1]
	movapd xmm5, [esp + nb112_fjzH1]
	mulsd xmm0, [esp + nb112_dxH1H1]
	mulsd xmm1, [esp + nb112_dyH1H1]
	mulsd xmm2, [esp + nb112_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixH1]
	addsd xmm1, [esp + nb112_fiyH1]
	addsd xmm2, [esp + nb112_fizH1]
	movlpd [esp + nb112_fjxH1], xmm3
	movlpd [esp + nb112_fjyH1], xmm4
	movlpd [esp + nb112_fjzH1], xmm5
	movlpd [esp + nb112_fixH1], xmm0
	movlpd [esp + nb112_fiyH1], xmm1
	movlpd [esp + nb112_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb112_rinvH1H2]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqHH]
	mulsd xmm0, xmm1	;# fsOH2  
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxH2]
	movapd xmm4, [esp + nb112_fjyH2]
	movapd xmm5, [esp + nb112_fjzH2]
	mulsd xmm0, [esp + nb112_dxH1H2]
	mulsd xmm1, [esp + nb112_dyH1H2]
	mulsd xmm2, [esp + nb112_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixH1]
	addsd xmm1, [esp + nb112_fiyH1]
	addsd xmm2, [esp + nb112_fizH1]
	movlpd [esp + nb112_fjxH2], xmm3
	movlpd [esp + nb112_fjyH2], xmm4
	movlpd [esp + nb112_fjzH2], xmm5
	movlpd [esp + nb112_fixH1], xmm0
	movlpd [esp + nb112_fiyH1], xmm1
	movlpd [esp + nb112_fizH1], xmm2

	;# H2-O interaction 
	movapd xmm0, [esp + nb112_rinvH2O]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqOH]
	mulsd xmm0, xmm1	;# fsH2O 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxO]
	movapd xmm4, [esp + nb112_fjyO]
	movapd xmm5, [esp + nb112_fjzO]
	mulsd xmm0, [esp + nb112_dxH2O]
	mulsd xmm1, [esp + nb112_dyH2O]
	mulsd xmm2, [esp + nb112_dzH2O]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixH2]
	addsd xmm1, [esp + nb112_fiyH2]
	addsd xmm2, [esp + nb112_fizH2]
	movlpd [esp + nb112_fjxO], xmm3
	movlpd [esp + nb112_fjyO], xmm4
	movlpd [esp + nb112_fjzO], xmm5
	movlpd [esp + nb112_fixH2], xmm0
	movlpd [esp + nb112_fiyH2], xmm1
	movlpd [esp + nb112_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb112_rinvH2H1]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqHH]
	mulsd xmm0, xmm1	;# fsH2H1 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	movapd xmm3, [esp + nb112_fjxH1]
	movapd xmm4, [esp + nb112_fjyH1]
	movapd xmm5, [esp + nb112_fjzH1]
	mulsd xmm0, [esp + nb112_dxH2H1]
	mulsd xmm1, [esp + nb112_dyH2H1]
	mulsd xmm2, [esp + nb112_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixH2]
	addsd xmm1, [esp + nb112_fiyH2]
	addsd xmm2, [esp + nb112_fizH2]
	movlpd [esp + nb112_fjxH1], xmm3
	movlpd [esp + nb112_fjyH1], xmm4
	movlpd [esp + nb112_fjzH1], xmm5
	movlpd [esp + nb112_fixH2], xmm0
	movlpd [esp + nb112_fiyH2], xmm1
	movlpd [esp + nb112_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb112_rinvH2H2]
	movapd xmm1, xmm0
	mulsd xmm0, xmm0
	mulsd xmm1, [esp + nb112_qqHH]
	mulsd xmm0, xmm1	;# fsH2H2 
	addsd xmm7, xmm1	;# add to local vctot 
	movapd xmm1, xmm0
	movlpd [esp + nb112_vctot], xmm7
	movapd xmm2, xmm0
	movlpd xmm3, [esp + nb112_fjxH2]
	movlpd xmm4, [esp + nb112_fjyH2]
	movlpd xmm5, [esp + nb112_fjzH2]
	mulsd xmm0, [esp + nb112_dxH2H2]
	mulsd xmm1, [esp + nb112_dyH2H2]
	mulsd xmm2, [esp + nb112_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb112_fixH2]
	addsd xmm1, [esp + nb112_fiyH2]
	addsd xmm2, [esp + nb112_fizH2]
	movlpd [esp + nb112_fjxH2], xmm3
	movlpd [esp + nb112_fjyH2], xmm4
	movlpd [esp + nb112_fjzH2], xmm5
	movlpd [esp + nb112_fixH2], xmm0
	movlpd [esp + nb112_fiyH2], xmm1
	movlpd [esp + nb112_fizH2], xmm2

	mov edi, [ebp + nb112_faction]
		
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]
	movlpd xmm3, [edi + eax*8 + 24]
	movlpd xmm4, [edi + eax*8 + 32]
	movlpd xmm5, [edi + eax*8 + 40]
	movlpd xmm6, [edi + eax*8 + 48]
	movlpd xmm7, [edi + eax*8 + 56]
	addsd xmm0, [esp + nb112_fjxO]
	addsd xmm1, [esp + nb112_fjyO]
	addsd xmm2, [esp + nb112_fjzO]
	addsd xmm3, [esp + nb112_fjxH1]
	addsd xmm4, [esp + nb112_fjyH1]
	addsd xmm5, [esp + nb112_fjzH1]
	addsd xmm6, [esp + nb112_fjxH2]
	addsd xmm7, [esp + nb112_fjyH2]
	movlpd [edi + eax*8], xmm0
	movlpd [edi + eax*8 + 8], xmm1
	movlpd [edi + eax*8 + 16], xmm2
	movlpd [edi + eax*8 + 24], xmm3
	movlpd [edi + eax*8 + 32], xmm4
	movlpd [edi + eax*8 + 40], xmm5
	movlpd [edi + eax*8 + 48], xmm6
	movlpd [edi + eax*8 + 56], xmm7

	movlpd xmm0, [edi + eax*8 + 64]
	addsd xmm0, [esp + nb112_fjzH2]
	movlpd [edi + eax*8 + 64], xmm0
.nb112_updateouterdata:
	mov   ecx, [esp + nb112_ii3]
	mov   edi, [ebp + nb112_faction]
	mov   esi, [ebp + nb112_fshift]
	mov   edx, [esp + nb112_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb112_fixO]
	movapd xmm1, [esp + nb112_fiyO]
	movapd xmm2, [esp + nb112_fizO]

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
	movapd xmm0, [esp + nb112_fixH1]
	movapd xmm1, [esp + nb112_fiyH1]
	movapd xmm2, [esp + nb112_fizH1]

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
	movapd xmm0, [esp + nb112_fixH2]
	movapd xmm1, [esp + nb112_fiyH2]
	movapd xmm2, [esp + nb112_fizH2]

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
	mov esi, [esp + nb112_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb112_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb112_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb112_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb112_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb112_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb112_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb112_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb112_n], esi
        jmp .nb112_outer
.nb112_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb112_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb112_end
        ;# non-zero, do one more workunit
        jmp   .nb112_threadloop
.nb112_end:
	emms

	mov eax, [esp + nb112_nouter]
	mov ebx, [esp + nb112_ninner]
	mov ecx, [ebp + nb112_outeriter]
	mov edx, [ebp + nb112_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb112_salign]
	add esp, eax
	add esp, 1512
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




	
.globl nb_kernel112nf_ia32_sse2
.globl _nb_kernel112nf_ia32_sse2
nb_kernel112nf_ia32_sse2:	
_nb_kernel112nf_ia32_sse2:	
.equiv          nb112nf_p_nri,          8
.equiv          nb112nf_iinr,           12
.equiv          nb112nf_jindex,         16
.equiv          nb112nf_jjnr,           20
.equiv          nb112nf_shift,          24
.equiv          nb112nf_shiftvec,       28
.equiv          nb112nf_fshift,         32
.equiv          nb112nf_gid,            36
.equiv          nb112nf_pos,            40
.equiv          nb112nf_faction,        44
.equiv          nb112nf_charge,         48
.equiv          nb112nf_p_facel,        52
.equiv          nb112nf_argkrf,         56
.equiv          nb112nf_argcrf,         60
.equiv          nb112nf_Vc,             64
.equiv          nb112nf_type,           68
.equiv          nb112nf_p_ntype,        72
.equiv          nb112nf_vdwparam,       76
.equiv          nb112nf_Vvdw,           80
.equiv          nb112nf_p_tabscale,     84
.equiv          nb112nf_VFtab,          88
.equiv          nb112nf_invsqrta,       92
.equiv          nb112nf_dvda,           96
.equiv          nb112nf_p_gbtabscale,   100
.equiv          nb112nf_GBtab,          104
.equiv          nb112nf_p_nthreads,     108
.equiv          nb112nf_count,          112
.equiv          nb112nf_mtx,            116
.equiv          nb112nf_outeriter,      120
.equiv          nb112nf_inneriter,      124
.equiv          nb112nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb112nf_ixO,            0
.equiv          nb112nf_iyO,            16
.equiv          nb112nf_izO,            32
.equiv          nb112nf_ixH1,           48
.equiv          nb112nf_iyH1,           64
.equiv          nb112nf_izH1,           80
.equiv          nb112nf_ixH2,           96
.equiv          nb112nf_iyH2,           112
.equiv          nb112nf_izH2,           128
.equiv          nb112nf_jxO,            144
.equiv          nb112nf_jyO,            160
.equiv          nb112nf_jzO,            176
.equiv          nb112nf_jxH1,           192
.equiv          nb112nf_jyH1,           208
.equiv          nb112nf_jzH1,           224
.equiv          nb112nf_jxH2,           240
.equiv          nb112nf_jyH2,           256
.equiv          nb112nf_jzH2,           272
.equiv          nb112nf_qqOO,           288
.equiv          nb112nf_qqOH,           304
.equiv          nb112nf_qqHH,           320
.equiv          nb112nf_c6,             336
.equiv          nb112nf_c12,            352
.equiv          nb112nf_vctot,          368
.equiv          nb112nf_Vvdwtot,        384
.equiv          nb112nf_half,           400
.equiv          nb112nf_three,          416
.equiv          nb112nf_rsqOO,          432
.equiv          nb112nf_rsqOH1,         448
.equiv          nb112nf_rsqOH2,         464
.equiv          nb112nf_rsqH1O,         480
.equiv          nb112nf_rsqH1H1,        496
.equiv          nb112nf_rsqH1H2,        512
.equiv          nb112nf_rsqH2O,         528
.equiv          nb112nf_rsqH2H1,        544
.equiv          nb112nf_rsqH2H2,        560
.equiv          nb112nf_rinvOO,         576
.equiv          nb112nf_rinvOH1,        592
.equiv          nb112nf_rinvOH2,        608
.equiv          nb112nf_rinvH1O,        624
.equiv          nb112nf_rinvH1H1,       640
.equiv          nb112nf_rinvH1H2,       656
.equiv          nb112nf_rinvH2O,        672
.equiv          nb112nf_rinvH2H1,       688
.equiv          nb112nf_rinvH2H2,       704
.equiv          nb112nf_is3,            720
.equiv          nb112nf_ii3,            724
.equiv          nb112nf_innerjjnr,      728
.equiv          nb112nf_innerk,         732
.equiv          nb112nf_n,              736
.equiv          nb112nf_nn1,            740
.equiv          nb112nf_nri,            744
.equiv          nb112nf_nouter,         748
.equiv          nb112nf_ninner,         752
.equiv          nb112nf_salign,         756
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 760		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb112nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb112nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb112nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb112nf_nouter], eax
	mov [esp + nb112nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb112nf_half], eax
	mov [esp + nb112nf_half+4], ebx
	movsd xmm1, [esp + nb112nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb112nf_half], xmm1
	movapd [esp + nb112nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb112nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb112nf_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb112nf_p_facel]
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
	movapd [esp + nb112nf_qqOO], xmm3
	movapd [esp + nb112nf_qqOH], xmm4
	movapd [esp + nb112nf_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb112nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb112nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb112nf_vdwparam]
	movlpd xmm0, [eax + edx*8]
	movhpd xmm0, [eax + edx*8 + 8]
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0	
	movapd [esp + nb112nf_c6], xmm0
	movapd [esp + nb112nf_c12], xmm1

.nb112nf_threadloop:
        mov   esi, [ebp + nb112nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb112nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb112nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb112nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb112nf_n], eax
        mov [esp + nb112nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb112nf_outerstart
        jmp .nb112nf_end

.nb112nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb112nf_nouter]
	mov [esp + nb112nf_nouter], ebx

.nb112nf_outer:
	mov   eax, [ebp + nb112nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb112nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb112nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb112nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb112nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb112nf_ixO], xmm3
	movapd [esp + nb112nf_iyO], xmm4
	movapd [esp + nb112nf_izO], xmm5

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
	movapd [esp + nb112nf_ixH1], xmm0
	movapd [esp + nb112nf_iyH1], xmm1
	movapd [esp + nb112nf_izH1], xmm2
	movapd [esp + nb112nf_ixH2], xmm3
	movapd [esp + nb112nf_iyH2], xmm4
	movapd [esp + nb112nf_izH2], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb112nf_vctot], xmm4
	movapd [esp + nb112nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb112nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb112nf_pos]
	mov   eax, [ebp + nb112nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb112nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb112nf_ninner]
	mov   [esp + nb112nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb112nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb112nf_unroll_loop
	jmp   .nb112nf_checksingle
.nb112nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb112nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb112nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb112nf_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb112nf_jxO], xmm2
	movapd 	[esp + nb112nf_jyO], xmm3
	movapd 	[esp + nb112nf_jzO], xmm4
	movapd 	[esp + nb112nf_jxH1], xmm5
	movapd 	[esp + nb112nf_jyH1], xmm6
	movapd 	[esp + nb112nf_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + ebx*8 + 56]
	movhpd xmm4, [esi + ebx*8 + 64]
	movapd 	[esp + nb112nf_jxH2], xmm2
	movapd 	[esp + nb112nf_jyH2], xmm3
	movapd 	[esp + nb112nf_jzH2], xmm4
	
	movapd xmm0, [esp + nb112nf_ixO]
	movapd xmm1, [esp + nb112nf_iyO]
	movapd xmm2, [esp + nb112nf_izO]
	movapd xmm3, [esp + nb112nf_ixO]
	movapd xmm4, [esp + nb112nf_iyO]
	movapd xmm5, [esp + nb112nf_izO]
	subpd  xmm0, [esp + nb112nf_jxO]
	subpd  xmm1, [esp + nb112nf_jyO]
	subpd  xmm2, [esp + nb112nf_jzO]
	subpd  xmm3, [esp + nb112nf_jxH1]
	subpd  xmm4, [esp + nb112nf_jyH1]
	subpd  xmm5, [esp + nb112nf_jzH1]
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
	movapd [esp + nb112nf_rsqOO], xmm0
	movapd [esp + nb112nf_rsqOH1], xmm3

	movapd xmm0, [esp + nb112nf_ixO]
	movapd xmm1, [esp + nb112nf_iyO]
	movapd xmm2, [esp + nb112nf_izO]
	movapd xmm3, [esp + nb112nf_ixH1]
	movapd xmm4, [esp + nb112nf_iyH1]
	movapd xmm5, [esp + nb112nf_izH1]
	subpd  xmm0, [esp + nb112nf_jxH2]
	subpd  xmm1, [esp + nb112nf_jyH2]
	subpd  xmm2, [esp + nb112nf_jzH2]
	subpd  xmm3, [esp + nb112nf_jxO]
	subpd  xmm4, [esp + nb112nf_jyO]
	subpd  xmm5, [esp + nb112nf_jzO]
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
	movapd [esp + nb112nf_rsqOH2], xmm0
	movapd [esp + nb112nf_rsqH1O], xmm3

	movapd xmm0, [esp + nb112nf_ixH1]
	movapd xmm1, [esp + nb112nf_iyH1]
	movapd xmm2, [esp + nb112nf_izH1]
	movapd xmm3, [esp + nb112nf_ixH1]
	movapd xmm4, [esp + nb112nf_iyH1]
	movapd xmm5, [esp + nb112nf_izH1]
	subpd  xmm0, [esp + nb112nf_jxH1]
	subpd  xmm1, [esp + nb112nf_jyH1]
	subpd  xmm2, [esp + nb112nf_jzH1]
	subpd  xmm3, [esp + nb112nf_jxH2]
	subpd  xmm4, [esp + nb112nf_jyH2]
	subpd  xmm5, [esp + nb112nf_jzH2]
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
	movapd [esp + nb112nf_rsqH1H1], xmm0
	movapd [esp + nb112nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb112nf_ixH2]
	movapd xmm1, [esp + nb112nf_iyH2]
	movapd xmm2, [esp + nb112nf_izH2]
	movapd xmm3, [esp + nb112nf_ixH2]
	movapd xmm4, [esp + nb112nf_iyH2]
	movapd xmm5, [esp + nb112nf_izH2]
	subpd  xmm0, [esp + nb112nf_jxO]
	subpd  xmm1, [esp + nb112nf_jyO]
	subpd  xmm2, [esp + nb112nf_jzO]
	subpd  xmm3, [esp + nb112nf_jxH1]
	subpd  xmm4, [esp + nb112nf_jyH1]
	subpd  xmm5, [esp + nb112nf_jzH1]
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
	movapd [esp + nb112nf_rsqH2O], xmm0
	movapd [esp + nb112nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb112nf_ixH2]
	movapd xmm1, [esp + nb112nf_iyH2]
	movapd xmm2, [esp + nb112nf_izH2]
	subpd  xmm0, [esp + nb112nf_jxH2]
	subpd  xmm1, [esp + nb112nf_jyH2]
	subpd  xmm2, [esp + nb112nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb112nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvH2H2], xmm1
	movapd [esp + nb112nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb112nf_rsqOO]
	movapd xmm4, [esp + nb112nf_rsqOH1]	
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
	movapd  xmm3, [esp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb112nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb112nf_half] ;# rinv
	movapd [esp + nb112nf_rinvOO], xmm1
	movapd [esp + nb112nf_rinvOH1], xmm5

	movapd xmm0, [esp + nb112nf_rsqOH2]
	movapd xmm4, [esp + nb112nf_rsqH1O]	
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
	movapd  xmm3, [esp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvOH2], xmm1
	movapd [esp + nb112nf_rinvH1O], xmm5

	movapd xmm0, [esp + nb112nf_rsqH1H1]
	movapd xmm4, [esp + nb112nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb112nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb112nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvH1H1], xmm1
	movapd [esp + nb112nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb112nf_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvH2O], xmm1

	;# start with OO interaction 
	movapd xmm0, [esp + nb112nf_rinvOO]
	movapd xmm7, xmm0
	mulpd  xmm0, xmm0
	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	mulpd  xmm7, [esp + nb112nf_qqOO]
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm1, [esp + nb112nf_c6]	
	mulpd  xmm2, [esp + nb112nf_c12]	
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addpd  xmm3, [esp + nb112nf_Vvdwtot]
	movapd [esp + nb112nf_Vvdwtot], xmm3
	addpd  xmm7, [esp + nb112nf_vctot]

	;# other interactions 
	movapd xmm1, [esp + nb112nf_rinvOH1]
	movapd xmm2, [esp + nb112nf_rinvH1H1]
	
	addpd xmm1, [esp + nb112nf_rinvOH2]
	addpd xmm2, [esp + nb112nf_rinvH1H2]
	
	addpd xmm1, [esp + nb112nf_rinvH1O]
	addpd xmm2, [esp + nb112nf_rinvH2H1]

	addpd xmm1, [esp + nb112nf_rinvH2O]
	addpd xmm2, [esp + nb112nf_rinvH2H2]

	mulpd xmm1, [esp + nb112nf_qqOH]
	mulpd xmm2, [esp + nb112nf_qqHH]
	
	addpd xmm7, xmm1	
	addpd xmm7, xmm2

	movapd [esp + nb112nf_vctot], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb112nf_innerk],  2
	jl    .nb112nf_checksingle
	jmp   .nb112nf_unroll_loop
.nb112nf_checksingle:
	mov   edx, [esp + nb112nf_innerk]
	and   edx, 1
	jnz   .nb112nf_dosingle
	jmp   .nb112nf_updateouterdata
.nb112nf_dosingle:
	mov   edx, [esp + nb112nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb112nf_innerjjnr],  4	

	mov esi, [ebp + nb112nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movapd 	[esp + nb112nf_jxO], xmm2
	movapd 	[esp + nb112nf_jyO], xmm3
	movapd 	[esp + nb112nf_jzO], xmm4
	movapd 	[esp + nb112nf_jxH1], xmm5
	movapd 	[esp + nb112nf_jyH1], xmm6
	movapd 	[esp + nb112nf_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movapd 	[esp + nb112nf_jxH2], xmm2
	movapd 	[esp + nb112nf_jyH2], xmm3
	movapd 	[esp + nb112nf_jzH2], xmm4
	
	movapd xmm0, [esp + nb112nf_ixO]
	movapd xmm1, [esp + nb112nf_iyO]
	movapd xmm2, [esp + nb112nf_izO]
	movapd xmm3, [esp + nb112nf_ixO]
	movapd xmm4, [esp + nb112nf_iyO]
	movapd xmm5, [esp + nb112nf_izO]
	subsd  xmm0, [esp + nb112nf_jxO]
	subsd  xmm1, [esp + nb112nf_jyO]
	subsd  xmm2, [esp + nb112nf_jzO]
	subsd  xmm3, [esp + nb112nf_jxH1]
	subsd  xmm4, [esp + nb112nf_jyH1]
	subsd  xmm5, [esp + nb112nf_jzH1]
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
	movapd [esp + nb112nf_rsqOO], xmm0
	movapd [esp + nb112nf_rsqOH1], xmm3

	movapd xmm0, [esp + nb112nf_ixO]
	movapd xmm1, [esp + nb112nf_iyO]
	movapd xmm2, [esp + nb112nf_izO]
	movapd xmm3, [esp + nb112nf_ixH1]
	movapd xmm4, [esp + nb112nf_iyH1]
	movapd xmm5, [esp + nb112nf_izH1]
	subsd  xmm0, [esp + nb112nf_jxH2]
	subsd  xmm1, [esp + nb112nf_jyH2]
	subsd  xmm2, [esp + nb112nf_jzH2]
	subsd  xmm3, [esp + nb112nf_jxO]
	subsd  xmm4, [esp + nb112nf_jyO]
	subsd  xmm5, [esp + nb112nf_jzO]
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
	movapd [esp + nb112nf_rsqOH2], xmm0
	movapd [esp + nb112nf_rsqH1O], xmm3

	movapd xmm0, [esp + nb112nf_ixH1]
	movapd xmm1, [esp + nb112nf_iyH1]
	movapd xmm2, [esp + nb112nf_izH1]
	movapd xmm3, [esp + nb112nf_ixH1]
	movapd xmm4, [esp + nb112nf_iyH1]
	movapd xmm5, [esp + nb112nf_izH1]
	subsd  xmm0, [esp + nb112nf_jxH1]
	subsd  xmm1, [esp + nb112nf_jyH1]
	subsd  xmm2, [esp + nb112nf_jzH1]
	subsd  xmm3, [esp + nb112nf_jxH2]
	subsd  xmm4, [esp + nb112nf_jyH2]
	subsd  xmm5, [esp + nb112nf_jzH2]
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
	movapd [esp + nb112nf_rsqH1H1], xmm0
	movapd [esp + nb112nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb112nf_ixH2]
	movapd xmm1, [esp + nb112nf_iyH2]
	movapd xmm2, [esp + nb112nf_izH2]
	movapd xmm3, [esp + nb112nf_ixH2]
	movapd xmm4, [esp + nb112nf_iyH2]
	movapd xmm5, [esp + nb112nf_izH2]
	subsd  xmm0, [esp + nb112nf_jxO]
	subsd  xmm1, [esp + nb112nf_jyO]
	subsd  xmm2, [esp + nb112nf_jzO]
	subsd  xmm3, [esp + nb112nf_jxH1]
	subsd  xmm4, [esp + nb112nf_jyH1]
	subsd  xmm5, [esp + nb112nf_jzH1]
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
	movapd [esp + nb112nf_rsqH2O], xmm0
	movapd [esp + nb112nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb112nf_ixH2]
	movapd xmm1, [esp + nb112nf_iyH2]
	movapd xmm2, [esp + nb112nf_izH2]
	subsd  xmm0, [esp + nb112nf_jxH2]
	subsd  xmm1, [esp + nb112nf_jyH2]
	subsd  xmm2, [esp + nb112nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb112nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvH2H2], xmm1
	movapd [esp + nb112nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb112nf_rsqOO]
	movapd xmm4, [esp + nb112nf_rsqOH1]	
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
	movapd  xmm3, [esp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb112nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb112nf_half] ;# rinv
	movapd [esp + nb112nf_rinvOO], xmm1
	movapd [esp + nb112nf_rinvOH1], xmm5

	movapd xmm0, [esp + nb112nf_rsqOH2]
	movapd xmm4, [esp + nb112nf_rsqH1O]	
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
	movapd  xmm3, [esp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvOH2], xmm1
	movapd [esp + nb112nf_rinvH1O], xmm5

	movapd xmm0, [esp + nb112nf_rsqH1H1]
	movapd xmm4, [esp + nb112nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb112nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb112nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvH1H1], xmm1
	movapd [esp + nb112nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb112nf_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb112nf_half] ;# rinv 
	movapd [esp + nb112nf_rinvH2O], xmm1

	;# start with OO interaction 
	movapd xmm0, [esp + nb112nf_rinvOO]
	movapd xmm7, xmm0
	mulsd  xmm0, xmm0
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	mulsd  xmm7, [esp + nb112nf_qqOO]
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm1, [esp + nb112nf_c6]	
	mulsd  xmm2, [esp + nb112nf_c12]	
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addsd  xmm3, [esp + nb112nf_Vvdwtot]
	movlpd [esp + nb112nf_Vvdwtot], xmm3
	addsd  xmm7, [esp + nb112nf_vctot]

	;# other interactions 
	movapd xmm1, [esp + nb112nf_rinvOH1]
	movapd xmm2, [esp + nb112nf_rinvH1H1]
	
	addsd xmm1, [esp + nb112nf_rinvOH2]
	addsd xmm2, [esp + nb112nf_rinvH1H2]
	
	addsd xmm1, [esp + nb112nf_rinvH1O]
	addsd xmm2, [esp + nb112nf_rinvH2H1]

	addsd xmm1, [esp + nb112nf_rinvH2O]
	addsd xmm2, [esp + nb112nf_rinvH2H2]

	mulsd xmm1, [esp + nb112nf_qqOH]
	mulsd xmm2, [esp + nb112nf_qqHH]
	
	addsd xmm7, xmm1	
	addsd xmm7, xmm2

	movlpd [esp + nb112nf_vctot], xmm7

.nb112nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb112nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb112nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb112nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb112nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb112nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb112nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb112nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb112nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb112nf_n], esi
        jmp .nb112nf_outer
.nb112nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb112nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb112nf_end
        ;# non-zero, do one more workunit
        jmp   .nb112nf_threadloop
.nb112nf_end:
	emms

	mov eax, [esp + nb112nf_nouter]
	mov ebx, [esp + nb112nf_ninner]
	mov ecx, [ebp + nb112nf_outeriter]
	mov edx, [ebp + nb112nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb112nf_salign]
	add esp, eax
	add esp, 760
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
