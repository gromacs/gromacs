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


.globl nb_kernel232_ia32_sse2
.globl _nb_kernel232_ia32_sse2
nb_kernel232_ia32_sse2:	
_nb_kernel232_ia32_sse2:	
.equiv          nb232_p_nri,            8
.equiv          nb232_iinr,             12
.equiv          nb232_jindex,           16
.equiv          nb232_jjnr,             20
.equiv          nb232_shift,            24
.equiv          nb232_shiftvec,         28
.equiv          nb232_fshift,           32
.equiv          nb232_gid,              36
.equiv          nb232_pos,              40
.equiv          nb232_faction,          44
.equiv          nb232_charge,           48
.equiv          nb232_p_facel,          52
.equiv          nb232_argkrf,           56
.equiv          nb232_argcrf,           60
.equiv          nb232_Vc,               64
.equiv          nb232_type,             68
.equiv          nb232_p_ntype,          72
.equiv          nb232_vdwparam,         76
.equiv          nb232_Vvdw,             80
.equiv          nb232_p_tabscale,       84
.equiv          nb232_VFtab,            88
.equiv          nb232_invsqrta,         92
.equiv          nb232_dvda,             96
.equiv          nb232_p_gbtabscale,     100
.equiv          nb232_GBtab,            104
.equiv          nb232_p_nthreads,       108
.equiv          nb232_count,            112
.equiv          nb232_mtx,              116
.equiv          nb232_outeriter,        120
.equiv          nb232_inneriter,        124
.equiv          nb232_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb232_ixO,              0
.equiv          nb232_iyO,              16
.equiv          nb232_izO,              32
.equiv          nb232_ixH1,             48
.equiv          nb232_iyH1,             64
.equiv          nb232_izH1,             80
.equiv          nb232_ixH2,             96
.equiv          nb232_iyH2,             112
.equiv          nb232_izH2,             128
.equiv          nb232_jxO,              144
.equiv          nb232_jyO,              160
.equiv          nb232_jzO,              176
.equiv          nb232_jxH1,             192
.equiv          nb232_jyH1,             208
.equiv          nb232_jzH1,             224
.equiv          nb232_jxH2,             240
.equiv          nb232_jyH2,             256
.equiv          nb232_jzH2,             272
.equiv          nb232_dxOO,             288
.equiv          nb232_dyOO,             304
.equiv          nb232_dzOO,             320
.equiv          nb232_dxOH1,            336
.equiv          nb232_dyOH1,            352
.equiv          nb232_dzOH1,            368
.equiv          nb232_dxOH2,            384
.equiv          nb232_dyOH2,            400
.equiv          nb232_dzOH2,            416
.equiv          nb232_dxH1O,            432
.equiv          nb232_dyH1O,            448
.equiv          nb232_dzH1O,            464
.equiv          nb232_dxH1H1,           480
.equiv          nb232_dyH1H1,           496
.equiv          nb232_dzH1H1,           512
.equiv          nb232_dxH1H2,           528
.equiv          nb232_dyH1H2,           544
.equiv          nb232_dzH1H2,           560
.equiv          nb232_dxH2O,            576
.equiv          nb232_dyH2O,            592
.equiv          nb232_dzH2O,            608
.equiv          nb232_dxH2H1,           624
.equiv          nb232_dyH2H1,           640
.equiv          nb232_dzH2H1,           656
.equiv          nb232_dxH2H2,           672
.equiv          nb232_dyH2H2,           688
.equiv          nb232_dzH2H2,           704
.equiv          nb232_qqOO,             720
.equiv          nb232_qqOH,             736
.equiv          nb232_qqHH,             752
.equiv          nb232_c6,               768
.equiv          nb232_c12,              784
.equiv          nb232_tsc,              800
.equiv          nb232_fstmp,            816
.equiv          nb232_vctot,            832
.equiv          nb232_Vvdwtot,          848
.equiv          nb232_fixO,             864
.equiv          nb232_fiyO,             880
.equiv          nb232_fizO,             896
.equiv          nb232_fixH1,            912
.equiv          nb232_fiyH1,            928
.equiv          nb232_fizH1,            944
.equiv          nb232_fixH2,            960
.equiv          nb232_fiyH2,            976
.equiv          nb232_fizH2,            992
.equiv          nb232_fjxO,             1008
.equiv          nb232_fjyO,             1024
.equiv          nb232_fjzO,             1040
.equiv          nb232_fjxH1,            1056
.equiv          nb232_fjyH1,            1072
.equiv          nb232_fjzH1,            1088
.equiv          nb232_fjxH2,            1104
.equiv          nb232_fjyH2,            1120
.equiv          nb232_fjzH2,            1136
.equiv          nb232_half,             1152
.equiv          nb232_three,            1168
.equiv          nb232_rsqOO,            1184
.equiv          nb232_rsqOH1,           1200
.equiv          nb232_rsqOH2,           1216
.equiv          nb232_rsqH1O,           1232
.equiv          nb232_rsqH1H1,          1248
.equiv          nb232_rsqH1H2,          1264
.equiv          nb232_rsqH2O,           1280
.equiv          nb232_rsqH2H1,          1296
.equiv          nb232_rsqH2H2,          1312
.equiv          nb232_rinvOO,           1328
.equiv          nb232_rinvOH1,          1344
.equiv          nb232_rinvOH2,          1360
.equiv          nb232_rinvH1O,          1376
.equiv          nb232_rinvH1H1,         1392
.equiv          nb232_rinvH1H2,         1408
.equiv          nb232_rinvH2O,          1424
.equiv          nb232_rinvH2H1,         1440
.equiv          nb232_rinvH2H2,         1456
.equiv          nb232_two,              1472
.equiv          nb232_krf,              1488
.equiv          nb232_crf,              1504
.equiv          nb232_is3,              1520
.equiv          nb232_ii3,              1524
.equiv          nb232_innerjjnr,        1528
.equiv          nb232_innerk,           1532
.equiv          nb232_n,                1536
.equiv          nb232_nn1,              1540
.equiv          nb232_nri,              1544
.equiv          nb232_nouter,           1548
.equiv          nb232_ninner,           1552
.equiv          nb232_salign,           1556
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
	mov [esp + nb232_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb232_p_nri]
	mov ecx, [ecx]
	mov [esp + nb232_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb232_nouter], eax
	mov [esp + nb232_ninner], eax
	
	mov eax, [ebp + nb232_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb232_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb232_half], eax
	mov [esp + nb232_half+4], ebx
	movsd xmm1, [esp + nb232_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb232_half], xmm1
	movapd [esp + nb232_two], xmm2
	movapd [esp + nb232_three], xmm3

	mov esi, [ebp + nb232_argkrf]
	mov edi, [ebp + nb232_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb232_krf], xmm5
	movapd [esp + nb232_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb232_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb232_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb232_p_facel]
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
	movapd [esp + nb232_qqOO], xmm3
	movapd [esp + nb232_qqOH], xmm4
	movapd [esp + nb232_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb232_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb232_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb232_vdwparam]
	movlpd xmm0, [eax + edx*8] 
	movhpd xmm0, [eax + edx*8 + 8] 
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb232_c6], xmm0
	movapd [esp + nb232_c12], xmm1

.nb232_threadloop:
        mov   esi, [ebp + nb232_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb232_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb232_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb232_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb232_n], eax
        mov [esp + nb232_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb232_outerstart
        jmp .nb232_end

.nb232_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb232_nouter]
	mov [esp + nb232_nouter], ebx

.nb232_outer:
	mov   eax, [ebp + nb232_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb232_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb232_shiftvec]   ;# eax = base of shiftvec[] 

	movlpd xmm0, [eax + ebx*8]
	movlpd xmm1, [eax + ebx*8 + 8]
	movlpd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb232_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb232_pos]    ;# eax = base of pos[]  
	mov   [esp + nb232_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb232_ixO], xmm3
	movapd [esp + nb232_iyO], xmm4
	movapd [esp + nb232_izO], xmm5

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
	movapd [esp + nb232_ixH1], xmm0
	movapd [esp + nb232_iyH1], xmm1
	movapd [esp + nb232_izH1], xmm2
	movapd [esp + nb232_ixH2], xmm3
	movapd [esp + nb232_iyH2], xmm4
	movapd [esp + nb232_izH2], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb232_vctot], xmm4
	movapd [esp + nb232_Vvdwtot], xmm4
	movapd [esp + nb232_fixO], xmm4
	movapd [esp + nb232_fiyO], xmm4
	movapd [esp + nb232_fizO], xmm4
	movapd [esp + nb232_fixH1], xmm4
	movapd [esp + nb232_fiyH1], xmm4
	movapd [esp + nb232_fizH1], xmm4
	movapd [esp + nb232_fixH2], xmm4
	movapd [esp + nb232_fiyH2], xmm4
	movapd [esp + nb232_fizH2], xmm4
	
	mov   eax, [ebp + nb232_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb232_pos]
	mov   edi, [ebp + nb232_faction]	
	mov   eax, [ebp + nb232_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb232_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb232_ninner]
	mov   [esp + nb232_ninner], ecx
	add   edx, 0
	mov   [esp + nb232_innerk], edx    ;# number of innerloop atoms 
	jge   .nb232_unroll_loop
	jmp   .nb232_checksingle
.nb232_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb232_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb232_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb232_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb232_jxO], xmm2
	movapd 	[esp + nb232_jyO], xmm3
	movapd 	[esp + nb232_jzO], xmm4
	movapd 	[esp + nb232_jxH1], xmm5
	movapd 	[esp + nb232_jyH1], xmm6
	movapd 	[esp + nb232_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + ebx*8 + 56]
	movhpd xmm4, [esi + ebx*8 + 64]
	movapd 	[esp + nb232_jxH2], xmm2
	movapd 	[esp + nb232_jyH2], xmm3
	movapd 	[esp + nb232_jzH2], xmm4
	
	movapd xmm0, [esp + nb232_ixO]
	movapd xmm1, [esp + nb232_iyO]
	movapd xmm2, [esp + nb232_izO]
	movapd xmm3, [esp + nb232_ixO]
	movapd xmm4, [esp + nb232_iyO]
	movapd xmm5, [esp + nb232_izO]
	subpd  xmm0, [esp + nb232_jxO]
	subpd  xmm1, [esp + nb232_jyO]
	subpd  xmm2, [esp + nb232_jzO]
	subpd  xmm3, [esp + nb232_jxH1]
	subpd  xmm4, [esp + nb232_jyH1]
	subpd  xmm5, [esp + nb232_jzH1]
	movapd [esp + nb232_dxOO], xmm0
	movapd [esp + nb232_dyOO], xmm1
	movapd [esp + nb232_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb232_dxOH1], xmm3
	movapd [esp + nb232_dyOH1], xmm4
	movapd [esp + nb232_dzOH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb232_rsqOO], xmm0
	movapd [esp + nb232_rsqOH1], xmm3

	movapd xmm0, [esp + nb232_ixO]
	movapd xmm1, [esp + nb232_iyO]
	movapd xmm2, [esp + nb232_izO]
	movapd xmm3, [esp + nb232_ixH1]
	movapd xmm4, [esp + nb232_iyH1]
	movapd xmm5, [esp + nb232_izH1]
	subpd  xmm0, [esp + nb232_jxH2]
	subpd  xmm1, [esp + nb232_jyH2]
	subpd  xmm2, [esp + nb232_jzH2]
	subpd  xmm3, [esp + nb232_jxO]
	subpd  xmm4, [esp + nb232_jyO]
	subpd  xmm5, [esp + nb232_jzO]
	movapd [esp + nb232_dxOH2], xmm0
	movapd [esp + nb232_dyOH2], xmm1
	movapd [esp + nb232_dzOH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb232_dxH1O], xmm3
	movapd [esp + nb232_dyH1O], xmm4
	movapd [esp + nb232_dzH1O], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb232_rsqOH2], xmm0
	movapd [esp + nb232_rsqH1O], xmm3

	movapd xmm0, [esp + nb232_ixH1]
	movapd xmm1, [esp + nb232_iyH1]
	movapd xmm2, [esp + nb232_izH1]
	movapd xmm3, [esp + nb232_ixH1]
	movapd xmm4, [esp + nb232_iyH1]
	movapd xmm5, [esp + nb232_izH1]
	subpd  xmm0, [esp + nb232_jxH1]
	subpd  xmm1, [esp + nb232_jyH1]
	subpd  xmm2, [esp + nb232_jzH1]
	subpd  xmm3, [esp + nb232_jxH2]
	subpd  xmm4, [esp + nb232_jyH2]
	subpd  xmm5, [esp + nb232_jzH2]
	movapd [esp + nb232_dxH1H1], xmm0
	movapd [esp + nb232_dyH1H1], xmm1
	movapd [esp + nb232_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb232_dxH1H2], xmm3
	movapd [esp + nb232_dyH1H2], xmm4
	movapd [esp + nb232_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [esp + nb232_rsqH1H1], xmm0
	movapd [esp + nb232_rsqH1H2], xmm3

	movapd xmm0, [esp + nb232_ixH2]
	movapd xmm1, [esp + nb232_iyH2]
	movapd xmm2, [esp + nb232_izH2]
	movapd xmm3, [esp + nb232_ixH2]
	movapd xmm4, [esp + nb232_iyH2]
	movapd xmm5, [esp + nb232_izH2]
	subpd  xmm0, [esp + nb232_jxO]
	subpd  xmm1, [esp + nb232_jyO]
	subpd  xmm2, [esp + nb232_jzO]
	subpd  xmm3, [esp + nb232_jxH1]
	subpd  xmm4, [esp + nb232_jyH1]
	subpd  xmm5, [esp + nb232_jzH1]
	movapd [esp + nb232_dxH2O], xmm0
	movapd [esp + nb232_dyH2O], xmm1
	movapd [esp + nb232_dzH2O], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [esp + nb232_dxH2H1], xmm3
	movapd [esp + nb232_dyH2H1], xmm4
	movapd [esp + nb232_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [esp + nb232_rsqH2O], xmm0
	movapd [esp + nb232_rsqH2H1], xmm4

	movapd xmm0, [esp + nb232_ixH2]
	movapd xmm1, [esp + nb232_iyH2]
	movapd xmm2, [esp + nb232_izH2]
	subpd  xmm0, [esp + nb232_jxH2]
	subpd  xmm1, [esp + nb232_jyH2]
	subpd  xmm2, [esp + nb232_jzH2]
	movapd [esp + nb232_dxH2H2], xmm0
	movapd [esp + nb232_dyH2H2], xmm1
	movapd [esp + nb232_dzH2H2], xmm2
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb232_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb232_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232_half] ;# iter1 
	mulpd   xmm7, [esp + nb232_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232_half] ;# rinv 
	mulpd   xmm5, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvH2H2], xmm1
	movapd [esp + nb232_rinvH2H1], xmm5

	movapd xmm0, [esp + nb232_rsqOO]
	movapd xmm4, [esp + nb232_rsqOH1]	
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
	movapd  xmm3, [esp + nb232_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb232_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232_half] ;# rinv 
	mulpd   xmm5, [esp + nb232_half] ;# rinv
	movapd [esp + nb232_rinvOO], xmm1
	movapd [esp + nb232_rinvOH1], xmm5

	movapd xmm0, [esp + nb232_rsqOH2]
	movapd xmm4, [esp + nb232_rsqH1O]	
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
	movapd  xmm3, [esp + nb232_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232_half] ;# iter1 
	mulpd   xmm7, [esp + nb232_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232_half] ;# rinv 
	mulpd   xmm5, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvOH2], xmm1
	movapd [esp + nb232_rinvH1O], xmm5

	movapd xmm0, [esp + nb232_rsqH1H1]
	movapd xmm4, [esp + nb232_rsqH1H2]	
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
	movapd  xmm3, [esp + nb232_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232_half] ;# iter1a 
	mulpd   xmm7, [esp + nb232_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232_half] ;# rinv 
	mulpd   xmm5, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvH1H1], xmm1
	movapd [esp + nb232_rinvH1H2], xmm5

	movapd xmm0, [esp + nb232_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb232_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb232_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb232_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvH2O], xmm1
	
	;# start with OO interaction 	
	movapd xmm0, [esp + nb232_rinvOO]
	movapd xmm7, xmm0		;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]
	mulpd  xmm5, [esp + nb232_rsqOO]	;# xmm5=krsq 
	movapd xmm6, xmm5
	addpd  xmm6, xmm7		;# xmm6=rinv+ krsq 
	subpd  xmm6, [esp + nb232_crf]	;# rinv+krsq-crf 	
	mulpd  xmm6, [esp + nb232_qqOO]	;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulpd  xmm5, [esp + nb232_two]	;# 2*krsq 
	subpd  xmm7, xmm5		;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqOO]	;# xmm7 = qq*(rinv-2*krsq) 
	mulpd  xmm7, xmm0
	
	addpd  xmm6, [esp + nb232_vctot]
	movapd [esp + nb232_vctot], xmm6
	movapd [esp + nb232_fstmp], xmm7

	;# LJ table interaction
	movapd xmm4, [esp + nb232_rsqOO]
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb232_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb232_VFtab]
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
	mulpd  xmm7, [esp + nb232_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb232_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb232_Vvdwtot]
	movapd xmm3, [esp + nb232_fstmp]
	mulpd  xmm7, [esp + nb232_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb232_fstmp], xmm3
	movapd [esp + nb232_Vvdwtot], xmm5

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
	mulpd  xmm7, [esp + nb232_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb232_c12]
	mulpd  xmm7, xmm4 
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb232_Vvdwtot]
	movapd xmm3, [esp + nb232_fstmp]
	mulpd  xmm7, [esp + nb232_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb232_Vvdwtot], xmm5

	mulpd  xmm0, xmm3
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb232_faction]

	;# update forces
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb232_dxOO]
	mulpd xmm1, [esp + nb232_dyOO]
	mulpd xmm2, [esp + nb232_dzOO]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixO]
	addpd xmm1, [esp + nb232_fiyO]
	addpd xmm2, [esp + nb232_fizO]
	movapd [esp + nb232_fjxO], xmm3
	movapd [esp + nb232_fjyO], xmm4
	movapd [esp + nb232_fjzO], xmm5
	movapd [esp + nb232_fixO], xmm0
	movapd [esp + nb232_fiyO], xmm1
	movapd [esp + nb232_fizO], xmm2

	;# O-H1 interaction 
	movapd xmm0, [esp + nb232_rinvOH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqOH1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulpd  xmm0, xmm0
	subpd  xmm4, [esp + nb232_crf]
	mulpd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	movapd xmm6, [esp + nb232_vctot]
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH1  
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb232_dxOH1]
	mulpd xmm1, [esp + nb232_dyOH1]
	mulpd xmm2, [esp + nb232_dzOH1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixO]
	addpd xmm1, [esp + nb232_fiyO]
	addpd xmm2, [esp + nb232_fizO]
	movapd [esp + nb232_fjxH1], xmm3
	movapd [esp + nb232_fjyH1], xmm4
	movapd [esp + nb232_fjzH1], xmm5
	movapd [esp + nb232_fixO], xmm0
	movapd [esp + nb232_fiyO], xmm1
	movapd [esp + nb232_fizO], xmm2

	;# O-H2 interaction  
	movapd xmm0, [esp + nb232_rinvOH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqOH2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb232_crf]
	mulpd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulpd xmm0, [esp + nb232_dxOH2]
	mulpd xmm1, [esp + nb232_dyOH2]
	mulpd xmm2, [esp + nb232_dzOH2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixO]
	addpd xmm1, [esp + nb232_fiyO]
	addpd xmm2, [esp + nb232_fizO]
	movapd [esp + nb232_fjxH2], xmm3
	movapd [esp + nb232_fjyH2], xmm4
	movapd [esp + nb232_fjzH2], xmm5
	movapd [esp + nb232_fixO], xmm0
	movapd [esp + nb232_fiyO], xmm1
	movapd [esp + nb232_fizO], xmm2

	;# H1-O interaction 
	movapd xmm0, [esp + nb232_rinvH1O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqH1O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb232_crf]
	mulpd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxO]
	movapd xmm4, [esp + nb232_fjyO]
	movapd xmm5, [esp + nb232_fjzO]
	mulpd xmm0, [esp + nb232_dxH1O]
	mulpd xmm1, [esp + nb232_dyH1O]
	mulpd xmm2, [esp + nb232_dzH1O]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixH1]
	addpd xmm1, [esp + nb232_fiyH1]
	addpd xmm2, [esp + nb232_fizH1]
	movapd [esp + nb232_fjxO], xmm3
	movapd [esp + nb232_fjyO], xmm4
	movapd [esp + nb232_fjzO], xmm5
	movapd [esp + nb232_fixH1], xmm0
	movapd [esp + nb232_fiyH1], xmm1
	movapd [esp + nb232_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb232_rinvH1H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqH1H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxH1]
	movapd xmm4, [esp + nb232_fjyH1]
	movapd xmm5, [esp + nb232_fjzH1]
	mulpd xmm0, [esp + nb232_dxH1H1]
	mulpd xmm1, [esp + nb232_dyH1H1]
	mulpd xmm2, [esp + nb232_dzH1H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixH1]
	addpd xmm1, [esp + nb232_fiyH1]
	addpd xmm2, [esp + nb232_fizH1]
	movapd [esp + nb232_fjxH1], xmm3
	movapd [esp + nb232_fjyH1], xmm4
	movapd [esp + nb232_fjzH1], xmm5
	movapd [esp + nb232_fixH1], xmm0
	movapd [esp + nb232_fiyH1], xmm1
	movapd [esp + nb232_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb232_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqH1H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb232_crf]
	mulpd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb232_fjxH2]
	movapd xmm4, [esp + nb232_fjyH2]
	movapd xmm5, [esp + nb232_fjzH2]
	mulpd xmm0, [esp + nb232_dxH1H2]
	mulpd xmm1, [esp + nb232_dyH1H2]
	mulpd xmm2, [esp + nb232_dzH1H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixH1]
	addpd xmm1, [esp + nb232_fiyH1]
	addpd xmm2, [esp + nb232_fizH1]
	movapd [esp + nb232_fjxH2], xmm3
	movapd [esp + nb232_fjyH2], xmm4
	movapd [esp + nb232_fjzH2], xmm5
	movapd [esp + nb232_fixH1], xmm0
	movapd [esp + nb232_fiyH1], xmm1
	movapd [esp + nb232_fizH1], xmm2

	;# H2-O interaction 
	movapd xmm0, [esp + nb232_rinvH2O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqH2O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxO]
	movapd xmm4, [esp + nb232_fjyO]
	movapd xmm5, [esp + nb232_fjzO]
	mulpd xmm0, [esp + nb232_dxH2O]
	mulpd xmm1, [esp + nb232_dyH2O]
	mulpd xmm2, [esp + nb232_dzH2O]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixH2]
	addpd xmm1, [esp + nb232_fiyH2]
	addpd xmm2, [esp + nb232_fizH2]
	movapd [esp + nb232_fjxO], xmm3
	movapd [esp + nb232_fjyO], xmm4
	movapd [esp + nb232_fjzO], xmm5
	movapd [esp + nb232_fixH2], xmm0
	movapd [esp + nb232_fiyH2], xmm1
	movapd [esp + nb232_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb232_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqH2H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxH1]
	movapd xmm4, [esp + nb232_fjyH1]
	movapd xmm5, [esp + nb232_fjzH1]
	mulpd xmm0, [esp + nb232_dxH2H1]
	mulpd xmm1, [esp + nb232_dyH2H1]
	mulpd xmm2, [esp + nb232_dzH2H1]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixH2]
	addpd xmm1, [esp + nb232_fiyH2]
	addpd xmm2, [esp + nb232_fizH2]
	movapd [esp + nb232_fjxH1], xmm3
	movapd [esp + nb232_fjyH1], xmm4
	movapd [esp + nb232_fjzH1], xmm5
	movapd [esp + nb232_fixH2], xmm0
	movapd [esp + nb232_fiyH2], xmm1
	movapd [esp + nb232_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb232_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232_rsqH2H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulpd  xmm5, [esp + nb232_two]
	subpd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulpd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addpd  xmm6, xmm4	;# add to local vctot 
	mulpd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm1, xmm0
	movapd [esp + nb232_vctot], xmm6
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb232_fjxH2]
	movapd xmm4, [esp + nb232_fjyH2]
	movapd xmm5, [esp + nb232_fjzH2]
	mulpd xmm0, [esp + nb232_dxH2H2]
	mulpd xmm1, [esp + nb232_dyH2H2]
	mulpd xmm2, [esp + nb232_dzH2H2]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	addpd xmm0, [esp + nb232_fixH2]
	addpd xmm1, [esp + nb232_fiyH2]
	addpd xmm2, [esp + nb232_fizH2]
	movapd [esp + nb232_fjxH2], xmm3
	movapd [esp + nb232_fjyH2], xmm4
	movapd [esp + nb232_fjzH2], xmm5
	movapd [esp + nb232_fixH2], xmm0
	movapd [esp + nb232_fiyH2], xmm1
	movapd [esp + nb232_fizH2], xmm2

	mov edi, [ebp + nb232_faction]
		
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
	addpd xmm0, [esp + nb232_fjxO]
	addpd xmm1, [esp + nb232_fjyO]
	addpd xmm2, [esp + nb232_fjzO]
	addpd xmm3, [esp + nb232_fjxH1]
	addpd xmm4, [esp + nb232_fjyH1]
	addpd xmm5, [esp + nb232_fjzH1]
	addpd xmm6, [esp + nb232_fjxH2]
	addpd xmm7, [esp + nb232_fjyH2]
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
	addpd xmm0, [esp + nb232_fjzH2]
	movlpd [edi + eax*8 + 64], xmm0
	movhpd [edi + ebx*8 + 64], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb232_innerk],  2
	jl    .nb232_checksingle
	jmp   .nb232_unroll_loop
.nb232_checksingle:
	mov   edx, [esp + nb232_innerk]
	and   edx, 1
	jnz   .nb232_dosingle
	jmp   .nb232_updateouterdata
.nb232_dosingle:
	mov   edx, [esp + nb232_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]
	
	mov esi, [ebp + nb232_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movapd 	[esp + nb232_jxO], xmm2
	movapd 	[esp + nb232_jyO], xmm3
	movapd 	[esp + nb232_jzO], xmm4
	movapd 	[esp + nb232_jxH1], xmm5
	movapd 	[esp + nb232_jyH1], xmm6
	movapd 	[esp + nb232_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movapd 	[esp + nb232_jxH2], xmm2
	movapd 	[esp + nb232_jyH2], xmm3
	movapd 	[esp + nb232_jzH2], xmm4
	
	movapd xmm0, [esp + nb232_ixO]
	movapd xmm1, [esp + nb232_iyO]
	movapd xmm2, [esp + nb232_izO]
	movapd xmm3, [esp + nb232_ixO]
	movapd xmm4, [esp + nb232_iyO]
	movapd xmm5, [esp + nb232_izO]
	subsd  xmm0, [esp + nb232_jxO]
	subsd  xmm1, [esp + nb232_jyO]
	subsd  xmm2, [esp + nb232_jzO]
	subsd  xmm3, [esp + nb232_jxH1]
	subsd  xmm4, [esp + nb232_jyH1]
	subsd  xmm5, [esp + nb232_jzH1]
	movapd [esp + nb232_dxOO], xmm0
	movapd [esp + nb232_dyOO], xmm1
	movapd [esp + nb232_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb232_dxOH1], xmm3
	movapd [esp + nb232_dyOH1], xmm4
	movapd [esp + nb232_dzOH1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb232_rsqOO], xmm0
	movapd [esp + nb232_rsqOH1], xmm3

	movapd xmm0, [esp + nb232_ixO]
	movapd xmm1, [esp + nb232_iyO]
	movapd xmm2, [esp + nb232_izO]
	movapd xmm3, [esp + nb232_ixH1]
	movapd xmm4, [esp + nb232_iyH1]
	movapd xmm5, [esp + nb232_izH1]
	subsd  xmm0, [esp + nb232_jxH2]
	subsd  xmm1, [esp + nb232_jyH2]
	subsd  xmm2, [esp + nb232_jzH2]
	subsd  xmm3, [esp + nb232_jxO]
	subsd  xmm4, [esp + nb232_jyO]
	subsd  xmm5, [esp + nb232_jzO]
	movapd [esp + nb232_dxOH2], xmm0
	movapd [esp + nb232_dyOH2], xmm1
	movapd [esp + nb232_dzOH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb232_dxH1O], xmm3
	movapd [esp + nb232_dyH1O], xmm4
	movapd [esp + nb232_dzH1O], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb232_rsqOH2], xmm0
	movapd [esp + nb232_rsqH1O], xmm3

	movapd xmm0, [esp + nb232_ixH1]
	movapd xmm1, [esp + nb232_iyH1]
	movapd xmm2, [esp + nb232_izH1]
	movapd xmm3, [esp + nb232_ixH1]
	movapd xmm4, [esp + nb232_iyH1]
	movapd xmm5, [esp + nb232_izH1]
	subsd  xmm0, [esp + nb232_jxH1]
	subsd  xmm1, [esp + nb232_jyH1]
	subsd  xmm2, [esp + nb232_jzH1]
	subsd  xmm3, [esp + nb232_jxH2]
	subsd  xmm4, [esp + nb232_jyH2]
	subsd  xmm5, [esp + nb232_jzH2]
	movapd [esp + nb232_dxH1H1], xmm0
	movapd [esp + nb232_dyH1H1], xmm1
	movapd [esp + nb232_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb232_dxH1H2], xmm3
	movapd [esp + nb232_dyH1H2], xmm4
	movapd [esp + nb232_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [esp + nb232_rsqH1H1], xmm0
	movapd [esp + nb232_rsqH1H2], xmm3

	movapd xmm0, [esp + nb232_ixH2]
	movapd xmm1, [esp + nb232_iyH2]
	movapd xmm2, [esp + nb232_izH2]
	movapd xmm3, [esp + nb232_ixH2]
	movapd xmm4, [esp + nb232_iyH2]
	movapd xmm5, [esp + nb232_izH2]
	subsd  xmm0, [esp + nb232_jxO]
	subsd  xmm1, [esp + nb232_jyO]
	subsd  xmm2, [esp + nb232_jzO]
	subsd  xmm3, [esp + nb232_jxH1]
	subsd  xmm4, [esp + nb232_jyH1]
	subsd  xmm5, [esp + nb232_jzH1]
	movapd [esp + nb232_dxH2O], xmm0
	movapd [esp + nb232_dyH2O], xmm1
	movapd [esp + nb232_dzH2O], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [esp + nb232_dxH2H1], xmm3
	movapd [esp + nb232_dyH2H1], xmm4
	movapd [esp + nb232_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [esp + nb232_rsqH2O], xmm0
	movapd [esp + nb232_rsqH2H1], xmm4

	movapd xmm0, [esp + nb232_ixH2]
	movapd xmm1, [esp + nb232_iyH2]
	movapd xmm2, [esp + nb232_izH2]
	subsd  xmm0, [esp + nb232_jxH2]
	subsd  xmm1, [esp + nb232_jyH2]
	subsd  xmm2, [esp + nb232_jzH2]
	movapd [esp + nb232_dxH2H2], xmm0
	movapd [esp + nb232_dyH2H2], xmm1
	movapd [esp + nb232_dzH2H2], xmm2
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb232_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb232_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232_half] ;# iter1 
	mulsd   xmm7, [esp + nb232_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232_half] ;# rinv 
	mulsd   xmm5, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvH2H2], xmm1
	movapd [esp + nb232_rinvH2H1], xmm5

	movapd xmm0, [esp + nb232_rsqOO]
	movapd xmm4, [esp + nb232_rsqOH1]	
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
	movapd  xmm3, [esp + nb232_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb232_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232_half] ;# rinv 
	mulsd   xmm5, [esp + nb232_half] ;# rinv
	movapd [esp + nb232_rinvOO], xmm1
	movapd [esp + nb232_rinvOH1], xmm5

	movapd xmm0, [esp + nb232_rsqOH2]
	movapd xmm4, [esp + nb232_rsqH1O]	
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
	movapd  xmm3, [esp + nb232_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232_half] ;# iter1 
	mulsd   xmm7, [esp + nb232_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232_half] ;# rinv 
	mulsd   xmm5, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvOH2], xmm1
	movapd [esp + nb232_rinvH1O], xmm5

	movapd xmm0, [esp + nb232_rsqH1H1]
	movapd xmm4, [esp + nb232_rsqH1H2]	
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
	movapd  xmm3, [esp + nb232_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232_half] ;# iter1a 
	mulsd   xmm7, [esp + nb232_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232_half] ;# rinv 
	mulsd   xmm5, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvH1H1], xmm1
	movapd [esp + nb232_rinvH1H2], xmm5

	movapd xmm0, [esp + nb232_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb232_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb232_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb232_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb232_half] ;# rinv 
	movapd [esp + nb232_rinvH2O], xmm1
	
	;# start with OO interaction 	
	movsd xmm0, [esp + nb232_rinvOO]
	movsd xmm7, xmm0		;# xmm7=rinv 
	movsd xmm5, [esp + nb232_krf]
	mulsd  xmm5, [esp + nb232_rsqOO]	;# xmm5=krsq 
	movsd xmm6, xmm5		;# krsq 
	addsd  xmm6, xmm7		;# xmm6=rinv+ krsq 
	subsd  xmm6, [esp + nb232_crf]	;# rinv+krsq-crf 	
	mulsd  xmm6, [esp + nb232_qqOO]	;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulsd  xmm5, [esp + nb232_two]	;# 2*krsq 
	subsd  xmm7, xmm5		;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqOO]	;# xmm7 = qq*(rinv-2*krsq) 
	mulsd  xmm7, xmm0
	
	addsd  xmm6, [esp + nb232_vctot]
	movsd [esp + nb232_vctot], xmm6
	movsd [esp + nb232_fstmp], xmm7

	mulsd xmm7, xmm0
	
	;# LJ table interaction
	movsd xmm4, [esp + nb232_rsqOO]
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb232_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb232_VFtab]

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
	mulsd  xmm7, [esp + nb232_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb232_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb232_Vvdwtot]
	movsd xmm3, [esp + nb232_fstmp]
	mulsd  xmm7, [esp + nb232_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb232_fstmp], xmm3
	movsd [esp + nb232_Vvdwtot], xmm5

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
	mulsd  xmm7, [esp + nb232_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb232_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb232_Vvdwtot]
	movsd xmm3, [esp + nb232_fstmp]
	mulsd  xmm7, [esp + nb232_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb232_Vvdwtot], xmm5

	mulsd  xmm0, xmm3
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	mov    edi, [ebp + nb232_faction]

	;# update forces
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb232_dxOO]
	mulsd xmm1, [esp + nb232_dyOO]
	mulsd xmm2, [esp + nb232_dzOO]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixO]
	addsd xmm1, [esp + nb232_fiyO]
	addsd xmm2, [esp + nb232_fizO]
	movlpd [esp + nb232_fjxO], xmm3
	movlpd [esp + nb232_fjyO], xmm4
	movlpd [esp + nb232_fjzO], xmm5
	movlpd [esp + nb232_fixO], xmm0
	movlpd [esp + nb232_fiyO], xmm1
	movlpd [esp + nb232_fizO], xmm2

	;# O-H1 interaction 
	movapd xmm0, [esp + nb232_rinvOH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqOH1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulsd  xmm0, xmm0
	subsd  xmm4, [esp + nb232_crf]
	mulsd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	movsd  xmm6, [esp + nb232_vctot]
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH1  
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb232_dxOH1]
	mulsd xmm1, [esp + nb232_dyOH1]
	mulsd xmm2, [esp + nb232_dzOH1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixO]
	addsd xmm1, [esp + nb232_fiyO]
	addsd xmm2, [esp + nb232_fizO]
	movlpd [esp + nb232_fjxH1], xmm3
	movlpd [esp + nb232_fjyH1], xmm4
	movlpd [esp + nb232_fjzH1], xmm5
	movlpd [esp + nb232_fixO], xmm0
	movlpd [esp + nb232_fiyO], xmm1
	movlpd [esp + nb232_fizO], xmm2

	;# O-H2 interaction  
	movapd xmm0, [esp + nb232_rinvOH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqOH2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulsd  xmm0, xmm0
	subsd  xmm4, [esp + nb232_crf]
	mulsd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	xorpd xmm3, xmm3
	movapd xmm4, xmm3
	movapd xmm5, xmm3
	mulsd xmm0, [esp + nb232_dxOH2]
	mulsd xmm1, [esp + nb232_dyOH2]
	mulsd xmm2, [esp + nb232_dzOH2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixO]
	addsd xmm1, [esp + nb232_fiyO]
	addsd xmm2, [esp + nb232_fizO]
	movlpd [esp + nb232_fjxH2], xmm3
	movlpd [esp + nb232_fjyH2], xmm4
	movlpd [esp + nb232_fjzH2], xmm5
	movlpd [esp + nb232_fixO], xmm0
	movlpd [esp + nb232_fiyO], xmm1
	movlpd [esp + nb232_fizO], xmm2

	;# H1-O interaction 
	movapd xmm0, [esp + nb232_rinvH1O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqH1O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulsd xmm0, xmm0
	subsd  xmm4, [esp + nb232_crf]
	mulsd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxO]
	movapd xmm4, [esp + nb232_fjyO]
	movapd xmm5, [esp + nb232_fjzO]
	mulsd xmm0, [esp + nb232_dxH1O]
	mulsd xmm1, [esp + nb232_dyH1O]
	mulsd xmm2, [esp + nb232_dzH1O]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixH1]
	addsd xmm1, [esp + nb232_fiyH1]
	addsd xmm2, [esp + nb232_fizH1]
	movlpd [esp + nb232_fjxO], xmm3
	movlpd [esp + nb232_fjyO], xmm4
	movlpd [esp + nb232_fjzO], xmm5
	movlpd [esp + nb232_fixH1], xmm0
	movlpd [esp + nb232_fiyH1], xmm1
	movlpd [esp + nb232_fizH1], xmm2

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb232_rinvH1H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqH1H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxH1]
	movapd xmm4, [esp + nb232_fjyH1]
	movapd xmm5, [esp + nb232_fjzH1]
	mulsd xmm0, [esp + nb232_dxH1H1]
	mulsd xmm1, [esp + nb232_dyH1H1]
	mulsd xmm2, [esp + nb232_dzH1H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixH1]
	addsd xmm1, [esp + nb232_fiyH1]
	addsd xmm2, [esp + nb232_fizH1]
	movlpd [esp + nb232_fjxH1], xmm3
	movlpd [esp + nb232_fjyH1], xmm4
	movlpd [esp + nb232_fjzH1], xmm5
	movlpd [esp + nb232_fixH1], xmm0
	movlpd [esp + nb232_fiyH1], xmm1
	movlpd [esp + nb232_fizH1], xmm2

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb232_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqH1H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulsd xmm0, xmm0
	subsd  xmm4, [esp + nb232_crf]
	mulsd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb232_fjxH2]
	movapd xmm4, [esp + nb232_fjyH2]
	movapd xmm5, [esp + nb232_fjzH2]
	mulsd xmm0, [esp + nb232_dxH1H2]
	mulsd xmm1, [esp + nb232_dyH1H2]
	mulsd xmm2, [esp + nb232_dzH1H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixH1]
	addsd xmm1, [esp + nb232_fiyH1]
	addsd xmm2, [esp + nb232_fizH1]
	movlpd [esp + nb232_fjxH2], xmm3
	movlpd [esp + nb232_fjyH2], xmm4
	movlpd [esp + nb232_fjzH2], xmm5
	movlpd [esp + nb232_fixH1], xmm0
	movlpd [esp + nb232_fiyH1], xmm1
	movlpd [esp + nb232_fizH1], xmm2

	;# H2-O interaction 
	movapd xmm0, [esp + nb232_rinvH2O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqH2O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqOH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxO]
	movapd xmm4, [esp + nb232_fjyO]
	movapd xmm5, [esp + nb232_fjzO]
	mulsd xmm0, [esp + nb232_dxH2O]
	mulsd xmm1, [esp + nb232_dyH2O]
	mulsd xmm2, [esp + nb232_dzH2O]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixH2]
	addsd xmm1, [esp + nb232_fiyH2]
	addsd xmm2, [esp + nb232_fizH2]
	movlpd [esp + nb232_fjxO], xmm3
	movlpd [esp + nb232_fjyO], xmm4
	movlpd [esp + nb232_fjzO], xmm5
	movlpd [esp + nb232_fixH2], xmm0
	movlpd [esp + nb232_fiyH2], xmm1
	movlpd [esp + nb232_fizH2], xmm2

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb232_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqH2H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm3, [esp + nb232_fjxH1]
	movapd xmm4, [esp + nb232_fjyH1]
	movapd xmm5, [esp + nb232_fjzH1]
	mulsd xmm0, [esp + nb232_dxH2H1]
	mulsd xmm1, [esp + nb232_dyH2H1]
	mulsd xmm2, [esp + nb232_dzH2H1]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixH2]
	addsd xmm1, [esp + nb232_fiyH2]
	addsd xmm2, [esp + nb232_fizH2]
	movlpd [esp + nb232_fjxH1], xmm3
	movlpd [esp + nb232_fjyH1], xmm4
	movlpd [esp + nb232_fjzH1], xmm5
	movlpd [esp + nb232_fixH2], xmm0
	movlpd [esp + nb232_fiyH2], xmm1
	movlpd [esp + nb232_fizH2], xmm2

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb232_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232_rsqH2H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	mulsd  xmm5, [esp + nb232_two]
	subsd  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulsd  xmm7, [esp + nb232_qqHH] ;# xmm7 = coul part of fscal 
	addsd  xmm6, xmm4	;# add to local vctot 
	mulsd xmm0, xmm7	;# fsOH2 
	movapd xmm1, xmm0
	movapd xmm2, xmm0

	movapd xmm1, xmm0
	movsd [esp + nb232_vctot], xmm6
	movapd xmm2, xmm0
	
	movapd xmm3, [esp + nb232_fjxH2]
	movapd xmm4, [esp + nb232_fjyH2]
	movapd xmm5, [esp + nb232_fjzH2]
	mulsd xmm0, [esp + nb232_dxH2H2]
	mulsd xmm1, [esp + nb232_dyH2H2]
	mulsd xmm2, [esp + nb232_dzH2H2]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	addsd xmm0, [esp + nb232_fixH2]
	addsd xmm1, [esp + nb232_fiyH2]
	addsd xmm2, [esp + nb232_fizH2]
	movlpd [esp + nb232_fjxH2], xmm3
	movlpd [esp + nb232_fjyH2], xmm4
	movlpd [esp + nb232_fjzH2], xmm5
	movlpd [esp + nb232_fixH2], xmm0
	movlpd [esp + nb232_fiyH2], xmm1
	movlpd [esp + nb232_fizH2], xmm2

	mov edi, [ebp + nb232_faction]
	;# Did all interactions - now update j forces 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]
	movlpd xmm3, [edi + eax*8 + 24]
	movlpd xmm4, [edi + eax*8 + 32]
	movlpd xmm5, [edi + eax*8 + 40]
	movlpd xmm6, [edi + eax*8 + 48]
	movlpd xmm7, [edi + eax*8 + 56]
	addsd xmm0, [esp + nb232_fjxO]
	addsd xmm1, [esp + nb232_fjyO]
	addsd xmm2, [esp + nb232_fjzO]
	addsd xmm3, [esp + nb232_fjxH1]
	addsd xmm4, [esp + nb232_fjyH1]
	addsd xmm5, [esp + nb232_fjzH1]
	addsd xmm6, [esp + nb232_fjxH2]
	addsd xmm7, [esp + nb232_fjyH2]
	movlpd [edi + eax*8], xmm0
	movlpd [edi + eax*8 + 8], xmm1
	movlpd [edi + eax*8 + 16], xmm2
	movlpd [edi + eax*8 + 24], xmm3
	movlpd [edi + eax*8 + 32], xmm4
	movlpd [edi + eax*8 + 40], xmm5
	movlpd [edi + eax*8 + 48], xmm6
	movlpd [edi + eax*8 + 56], xmm7

	movlpd xmm0, [edi + eax*8 + 64]
	addsd xmm0, [esp + nb232_fjzH2]
	movlpd [edi + eax*8 + 64], xmm0
	
.nb232_updateouterdata:
	mov   ecx, [esp + nb232_ii3]
	mov   edi, [ebp + nb232_faction]
	mov   esi, [ebp + nb232_fshift]
	mov   edx, [esp + nb232_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb232_fixO]
	movapd xmm1, [esp + nb232_fiyO]
	movapd xmm2, [esp + nb232_fizO]

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
	movapd xmm0, [esp + nb232_fixH1]
	movapd xmm1, [esp + nb232_fiyH1]
	movapd xmm2, [esp + nb232_fizH1]

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
	movapd xmm0, [esp + nb232_fixH2]
	movapd xmm1, [esp + nb232_fiyH2]
	movapd xmm2, [esp + nb232_fizH2]

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
	mov esi, [esp + nb232_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb232_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb232_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb232_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb232_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb232_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb232_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb232_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb232_n], esi
        jmp .nb232_outer
.nb232_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb232_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb232_end
        ;# non-zero, do one more workunit
        jmp   .nb232_threadloop
.nb232_end:
	emms

	mov eax, [esp + nb232_nouter]
	mov ebx, [esp + nb232_ninner]
	mov ecx, [ebp + nb232_outeriter]
	mov edx, [ebp + nb232_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb232_salign]
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




.globl nb_kernel232nf_ia32_sse2
.globl _nb_kernel232nf_ia32_sse2
nb_kernel232nf_ia32_sse2:	
_nb_kernel232nf_ia32_sse2:	
.equiv          nb232nf_p_nri,            8
.equiv          nb232nf_iinr,             12
.equiv          nb232nf_jindex,           16
.equiv          nb232nf_jjnr,             20
.equiv          nb232nf_shift,            24
.equiv          nb232nf_shiftvec,         28
.equiv          nb232nf_fshift,           32
.equiv          nb232nf_gid,              36
.equiv          nb232nf_pos,              40
.equiv          nb232nf_faction,          44
.equiv          nb232nf_charge,           48
.equiv          nb232nf_p_facel,          52
.equiv          nb232nf_argkrf,           56
.equiv          nb232nf_argcrf,           60
.equiv          nb232nf_Vc,               64
.equiv          nb232nf_type,             68
.equiv          nb232nf_p_ntype,          72
.equiv          nb232nf_vdwparam,         76
.equiv          nb232nf_Vvdw,             80
.equiv          nb232nf_p_tabscale,       84
.equiv          nb232nf_VFtab,            88
.equiv          nb232nf_invsqrta,         92
.equiv          nb232nf_dvda,             96
.equiv          nb232nf_p_gbtabscale,     100
.equiv          nb232nf_GBtab,            104
.equiv          nb232nf_p_nthreads,       108
.equiv          nb232nf_count,            112
.equiv          nb232nf_mtx,              116
.equiv          nb232nf_outeriter,        120
.equiv          nb232nf_inneriter,        124
.equiv          nb232nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb232nf_ixO,              0
.equiv          nb232nf_iyO,              16
.equiv          nb232nf_izO,              32
.equiv          nb232nf_ixH1,             48
.equiv          nb232nf_iyH1,             64
.equiv          nb232nf_izH1,             80
.equiv          nb232nf_ixH2,             96
.equiv          nb232nf_iyH2,             112
.equiv          nb232nf_izH2,             128
.equiv          nb232nf_jxO,              144
.equiv          nb232nf_jyO,              160
.equiv          nb232nf_jzO,              176
.equiv          nb232nf_jxH1,             192
.equiv          nb232nf_jyH1,             208
.equiv          nb232nf_jzH1,             224
.equiv          nb232nf_jxH2,             240
.equiv          nb232nf_jyH2,             256
.equiv          nb232nf_jzH2,             272
.equiv          nb232nf_qqOO,             288
.equiv          nb232nf_qqOH,             304
.equiv          nb232nf_qqHH,             320
.equiv          nb232nf_c6,               336
.equiv          nb232nf_c12,              352
.equiv          nb232nf_tsc,              368
.equiv          nb232nf_vctot,            384
.equiv          nb232nf_Vvdwtot,          400
.equiv          nb232nf_half,             416
.equiv          nb232nf_three,            432
.equiv          nb232nf_rsqOO,            448
.equiv          nb232nf_rsqOH1,           464
.equiv          nb232nf_rsqOH2,           480
.equiv          nb232nf_rsqH1O,           496
.equiv          nb232nf_rsqH1H1,          512
.equiv          nb232nf_rsqH1H2,          528
.equiv          nb232nf_rsqH2O,           544
.equiv          nb232nf_rsqH2H1,          560
.equiv          nb232nf_rsqH2H2,          576
.equiv          nb232nf_rinvOO,           592
.equiv          nb232nf_rinvOH1,          608
.equiv          nb232nf_rinvOH2,          624
.equiv          nb232nf_rinvH1O,          640
.equiv          nb232nf_rinvH1H1,         656
.equiv          nb232nf_rinvH1H2,         672
.equiv          nb232nf_rinvH2O,          688
.equiv          nb232nf_rinvH2H1,         704
.equiv          nb232nf_rinvH2H2,         720
.equiv          nb232nf_two,              736
.equiv          nb232nf_krf,              752
.equiv          nb232nf_crf,              768
.equiv          nb232nf_is3,              784
.equiv          nb232nf_ii3,              788
.equiv          nb232nf_innerjjnr,        792
.equiv          nb232nf_innerk,           796
.equiv          nb232nf_n,                800
.equiv          nb232nf_nn1,              804
.equiv          nb232nf_nri,              808
.equiv          nb232nf_nouter,           812
.equiv          nb232nf_ninner,           816
.equiv          nb232nf_salign,           820
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 824		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb232nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb232nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb232nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb232nf_nouter], eax
	mov [esp + nb232nf_ninner], eax
	
	mov eax, [ebp + nb232nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb232nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb232nf_half], eax
	mov [esp + nb232nf_half+4], ebx
	movsd xmm1, [esp + nb232nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb232nf_half], xmm1
	movapd [esp + nb232nf_two], xmm2
	movapd [esp + nb232nf_three], xmm3

	mov esi, [ebp + nb232nf_argkrf]
	mov edi, [ebp + nb232nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb232nf_krf], xmm5
	movapd [esp + nb232nf_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb232nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb232nf_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb232nf_p_facel]
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
	movapd [esp + nb232nf_qqOO], xmm3
	movapd [esp + nb232nf_qqOH], xmm4
	movapd [esp + nb232nf_qqHH], xmm5
		
	xorpd xmm0, xmm0
	mov   edx, [ebp + nb232nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + nb232nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + nb232nf_vdwparam]
	movlpd xmm0, [eax + edx*8] 
	movhpd xmm0, [eax + edx*8 + 8] 
	movhlps xmm1, xmm0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [esp + nb232nf_c6], xmm0
	movapd [esp + nb232nf_c12], xmm1

.nb232nf_threadloop:
        mov   esi, [ebp + nb232nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb232nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb232nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb232nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb232nf_n], eax
        mov [esp + nb232nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb232nf_outerstart
        jmp .nb232nf_end

.nb232nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb232nf_nouter]
	mov [esp + nb232nf_nouter], ebx

.nb232nf_outer:
	mov   eax, [ebp + nb232nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb232nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb232nf_shiftvec]   ;# eax = base of shiftvec[] 

	movlpd xmm0, [eax + ebx*8]
	movlpd xmm1, [eax + ebx*8 + 8]
	movlpd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb232nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb232nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb232nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb232nf_ixO], xmm3
	movapd [esp + nb232nf_iyO], xmm4
	movapd [esp + nb232nf_izO], xmm5

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
	movapd [esp + nb232nf_ixH1], xmm0
	movapd [esp + nb232nf_iyH1], xmm1
	movapd [esp + nb232nf_izH1], xmm2
	movapd [esp + nb232nf_ixH2], xmm3
	movapd [esp + nb232nf_iyH2], xmm4
	movapd [esp + nb232nf_izH2], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb232nf_vctot], xmm4
	movapd [esp + nb232nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb232nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb232nf_pos]
	mov   edi, [ebp + nb232nf_faction]	
	mov   eax, [ebp + nb232nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb232nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb232nf_ninner]
	mov   [esp + nb232nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb232nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb232nf_unroll_loop
	jmp   .nb232nf_checksingle
.nb232nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb232nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	
	add dword ptr [esp + nb232nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb232nf_pos]       ;# base of pos[] 

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
	movapd 	[esp + nb232nf_jxO], xmm2
	movapd 	[esp + nb232nf_jyO], xmm3
	movapd 	[esp + nb232nf_jzO], xmm4
	movapd 	[esp + nb232nf_jxH1], xmm5
	movapd 	[esp + nb232nf_jyH1], xmm6
	movapd 	[esp + nb232nf_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movhpd xmm2, [esi + ebx*8 + 48]
	movhpd xmm3, [esi + ebx*8 + 56]
	movhpd xmm4, [esi + ebx*8 + 64]
	movapd 	[esp + nb232nf_jxH2], xmm2
	movapd 	[esp + nb232nf_jyH2], xmm3
	movapd 	[esp + nb232nf_jzH2], xmm4
	
	movapd xmm0, [esp + nb232nf_ixO]
	movapd xmm1, [esp + nb232nf_iyO]
	movapd xmm2, [esp + nb232nf_izO]
	movapd xmm3, [esp + nb232nf_ixO]
	movapd xmm4, [esp + nb232nf_iyO]
	movapd xmm5, [esp + nb232nf_izO]
	subpd  xmm0, [esp + nb232nf_jxO]
	subpd  xmm1, [esp + nb232nf_jyO]
	subpd  xmm2, [esp + nb232nf_jzO]
	subpd  xmm3, [esp + nb232nf_jxH1]
	subpd  xmm4, [esp + nb232nf_jyH1]
	subpd  xmm5, [esp + nb232nf_jzH1]
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
	movapd [esp + nb232nf_rsqOO], xmm0
	movapd [esp + nb232nf_rsqOH1], xmm3

	movapd xmm0, [esp + nb232nf_ixO]
	movapd xmm1, [esp + nb232nf_iyO]
	movapd xmm2, [esp + nb232nf_izO]
	movapd xmm3, [esp + nb232nf_ixH1]
	movapd xmm4, [esp + nb232nf_iyH1]
	movapd xmm5, [esp + nb232nf_izH1]
	subpd  xmm0, [esp + nb232nf_jxH2]
	subpd  xmm1, [esp + nb232nf_jyH2]
	subpd  xmm2, [esp + nb232nf_jzH2]
	subpd  xmm3, [esp + nb232nf_jxO]
	subpd  xmm4, [esp + nb232nf_jyO]
	subpd  xmm5, [esp + nb232nf_jzO]
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
	movapd [esp + nb232nf_rsqOH2], xmm0
	movapd [esp + nb232nf_rsqH1O], xmm3

	movapd xmm0, [esp + nb232nf_ixH1]
	movapd xmm1, [esp + nb232nf_iyH1]
	movapd xmm2, [esp + nb232nf_izH1]
	movapd xmm3, [esp + nb232nf_ixH1]
	movapd xmm4, [esp + nb232nf_iyH1]
	movapd xmm5, [esp + nb232nf_izH1]
	subpd  xmm0, [esp + nb232nf_jxH1]
	subpd  xmm1, [esp + nb232nf_jyH1]
	subpd  xmm2, [esp + nb232nf_jzH1]
	subpd  xmm3, [esp + nb232nf_jxH2]
	subpd  xmm4, [esp + nb232nf_jyH2]
	subpd  xmm5, [esp + nb232nf_jzH2]
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
	movapd [esp + nb232nf_rsqH1H1], xmm0
	movapd [esp + nb232nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb232nf_ixH2]
	movapd xmm1, [esp + nb232nf_iyH2]
	movapd xmm2, [esp + nb232nf_izH2]
	movapd xmm3, [esp + nb232nf_ixH2]
	movapd xmm4, [esp + nb232nf_iyH2]
	movapd xmm5, [esp + nb232nf_izH2]
	subpd  xmm0, [esp + nb232nf_jxO]
	subpd  xmm1, [esp + nb232nf_jyO]
	subpd  xmm2, [esp + nb232nf_jzO]
	subpd  xmm3, [esp + nb232nf_jxH1]
	subpd  xmm4, [esp + nb232nf_jyH1]
	subpd  xmm5, [esp + nb232nf_jzH1]
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
	movapd [esp + nb232nf_rsqH2O], xmm0
	movapd [esp + nb232nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb232nf_ixH2]
	movapd xmm1, [esp + nb232nf_iyH2]
	movapd xmm2, [esp + nb232nf_izH2]
	subpd  xmm0, [esp + nb232nf_jxH2]
	subpd  xmm1, [esp + nb232nf_jyH2]
	subpd  xmm2, [esp + nb232nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [esp + nb232nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb232nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb232nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvH2H2], xmm1
	movapd [esp + nb232nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb232nf_rsqOO]
	movapd xmm4, [esp + nb232nf_rsqOH1]	
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
	movapd  xmm3, [esp + nb232nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232nf_half] ;# iter1 of  
	mulpd   xmm7, [esp + nb232nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb232nf_half] ;# rinv
	movapd [esp + nb232nf_rinvOO], xmm1
	movapd [esp + nb232nf_rinvOH1], xmm5

	movapd xmm0, [esp + nb232nf_rsqOH2]
	movapd xmm4, [esp + nb232nf_rsqH1O]	
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
	movapd  xmm3, [esp + nb232nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232nf_half] ;# iter1 
	mulpd   xmm7, [esp + nb232nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvOH2], xmm1
	movapd [esp + nb232nf_rinvH1O], xmm5

	movapd xmm0, [esp + nb232nf_rsqH1H1]
	movapd xmm4, [esp + nb232nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb232nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [esp + nb232nf_half] ;# iter1a 
	mulpd   xmm7, [esp + nb232nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulpd   xmm5, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvH1H1], xmm1
	movapd [esp + nb232nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb232nf_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb232nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [esp + nb232nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb232nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvH2O], xmm1
	
	;# start with OO interaction 	
	movapd xmm0, [esp + nb232nf_rinvOO]
	movapd xmm7, xmm0		;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]
	mulpd  xmm5, [esp + nb232nf_rsqOO]	;# xmm5=krsq 
	movapd xmm6, xmm5
	addpd  xmm6, xmm7		;# xmm6=rinv+ krsq 
	subpd  xmm6, [esp + nb232nf_crf]	;# rinv+krsq-crf 	
	mulpd  xmm6, [esp + nb232nf_qqOO]	;# xmm6=voul=qq*(rinv+ krsq-crf) 
	
	addpd  xmm6, [esp + nb232nf_vctot]
	movapd [esp + nb232nf_vctot], xmm6

	;# LJ table interaction
	movapd xmm4, [esp + nb232nf_rsqOO]
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb232nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  esi, [ebp + nb232nf_VFtab]
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

	movapd xmm4, [esp + nb232nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [esp + nb232nf_Vvdwtot]
	movapd [esp + nb232nf_Vvdwtot], xmm5

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
	
	movapd xmm4, [esp + nb232nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb232nf_Vvdwtot]
	movapd [esp + nb232nf_Vvdwtot], xmm5

	;# O-H1 interaction 
	movapd xmm0, [esp + nb232nf_rinvOH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqOH1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulpd  xmm0, xmm0
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	
	movapd xmm6, [esp + nb232nf_vctot]
	addpd  xmm6, xmm4
	
	;# O-H2 interaction  
	movapd xmm0, [esp + nb232nf_rinvOH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqOH2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addpd  xmm6, xmm4	;# add to local vctot 

	;# H1-O interaction 
	movapd xmm0, [esp + nb232nf_rinvH1O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqH1O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addpd  xmm6, xmm4	;# add to local vctot 

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb232nf_rinvH1H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqH1H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addpd  xmm6, xmm4	;# add to local vctot 

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb232nf_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqH1H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulpd xmm0, xmm0
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addpd  xmm6, xmm4	;# add to local vctot 

	;# H2-O interaction 
	movapd xmm0, [esp + nb232nf_rinvH2O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqH2O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addpd  xmm6, xmm4	;# add to local vctot 

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb232nf_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqH2H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addpd  xmm6, xmm4	;# add to local vctot 

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb232nf_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulpd  xmm5, [esp + nb232nf_rsqH2H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addpd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subpd  xmm4, [esp + nb232nf_crf]
	mulpd xmm0, xmm0
	mulpd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addpd  xmm6, xmm4	;# add to local vctot 
	movapd [esp + nb232nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb232nf_innerk],  2
	jl    .nb232nf_checksingle
	jmp   .nb232nf_unroll_loop
.nb232nf_checksingle:
	mov   edx, [esp + nb232nf_innerk]
	and   edx, 1
	jnz   .nb232nf_dosingle
	jmp   .nb232nf_updateouterdata
.nb232nf_dosingle:
	mov   edx, [esp + nb232nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]
	
	mov esi, [ebp + nb232nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [esi + eax*8]
	movlpd xmm3, [esi + eax*8 + 8]
	movlpd xmm4, [esi + eax*8 + 16]
	movlpd xmm5, [esi + eax*8 + 24]
	movlpd xmm6, [esi + eax*8 + 32]
	movlpd xmm7, [esi + eax*8 + 40]
	movapd 	[esp + nb232nf_jxO], xmm2
	movapd 	[esp + nb232nf_jyO], xmm3
	movapd 	[esp + nb232nf_jzO], xmm4
	movapd 	[esp + nb232nf_jxH1], xmm5
	movapd 	[esp + nb232nf_jyH1], xmm6
	movapd 	[esp + nb232nf_jzH1], xmm7
	movlpd xmm2, [esi + eax*8 + 48]
	movlpd xmm3, [esi + eax*8 + 56]
	movlpd xmm4, [esi + eax*8 + 64]
	movapd 	[esp + nb232nf_jxH2], xmm2
	movapd 	[esp + nb232nf_jyH2], xmm3
	movapd 	[esp + nb232nf_jzH2], xmm4
	
	movapd xmm0, [esp + nb232nf_ixO]
	movapd xmm1, [esp + nb232nf_iyO]
	movapd xmm2, [esp + nb232nf_izO]
	movapd xmm3, [esp + nb232nf_ixO]
	movapd xmm4, [esp + nb232nf_iyO]
	movapd xmm5, [esp + nb232nf_izO]
	subsd  xmm0, [esp + nb232nf_jxO]
	subsd  xmm1, [esp + nb232nf_jyO]
	subsd  xmm2, [esp + nb232nf_jzO]
	subsd  xmm3, [esp + nb232nf_jxH1]
	subsd  xmm4, [esp + nb232nf_jyH1]
	subsd  xmm5, [esp + nb232nf_jzH1]
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
	movapd [esp + nb232nf_rsqOO], xmm0
	movapd [esp + nb232nf_rsqOH1], xmm3

	movapd xmm0, [esp + nb232nf_ixO]
	movapd xmm1, [esp + nb232nf_iyO]
	movapd xmm2, [esp + nb232nf_izO]
	movapd xmm3, [esp + nb232nf_ixH1]
	movapd xmm4, [esp + nb232nf_iyH1]
	movapd xmm5, [esp + nb232nf_izH1]
	subsd  xmm0, [esp + nb232nf_jxH2]
	subsd  xmm1, [esp + nb232nf_jyH2]
	subsd  xmm2, [esp + nb232nf_jzH2]
	subsd  xmm3, [esp + nb232nf_jxO]
	subsd  xmm4, [esp + nb232nf_jyO]
	subsd  xmm5, [esp + nb232nf_jzO]
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
	movapd [esp + nb232nf_rsqOH2], xmm0
	movapd [esp + nb232nf_rsqH1O], xmm3

	movapd xmm0, [esp + nb232nf_ixH1]
	movapd xmm1, [esp + nb232nf_iyH1]
	movapd xmm2, [esp + nb232nf_izH1]
	movapd xmm3, [esp + nb232nf_ixH1]
	movapd xmm4, [esp + nb232nf_iyH1]
	movapd xmm5, [esp + nb232nf_izH1]
	subsd  xmm0, [esp + nb232nf_jxH1]
	subsd  xmm1, [esp + nb232nf_jyH1]
	subsd  xmm2, [esp + nb232nf_jzH1]
	subsd  xmm3, [esp + nb232nf_jxH2]
	subsd  xmm4, [esp + nb232nf_jyH2]
	subsd  xmm5, [esp + nb232nf_jzH2]
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
	movapd [esp + nb232nf_rsqH1H1], xmm0
	movapd [esp + nb232nf_rsqH1H2], xmm3

	movapd xmm0, [esp + nb232nf_ixH2]
	movapd xmm1, [esp + nb232nf_iyH2]
	movapd xmm2, [esp + nb232nf_izH2]
	movapd xmm3, [esp + nb232nf_ixH2]
	movapd xmm4, [esp + nb232nf_iyH2]
	movapd xmm5, [esp + nb232nf_izH2]
	subsd  xmm0, [esp + nb232nf_jxO]
	subsd  xmm1, [esp + nb232nf_jyO]
	subsd  xmm2, [esp + nb232nf_jzO]
	subsd  xmm3, [esp + nb232nf_jxH1]
	subsd  xmm4, [esp + nb232nf_jyH1]
	subsd  xmm5, [esp + nb232nf_jzH1]
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
	movapd [esp + nb232nf_rsqH2O], xmm0
	movapd [esp + nb232nf_rsqH2H1], xmm4

	movapd xmm0, [esp + nb232nf_ixH2]
	movapd xmm1, [esp + nb232nf_iyH2]
	movapd xmm2, [esp + nb232nf_izH2]
	subsd  xmm0, [esp + nb232nf_jxH2]
	subsd  xmm1, [esp + nb232nf_jyH2]
	subsd  xmm2, [esp + nb232nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [esp + nb232nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [esp + nb232nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb232nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvH2H2], xmm1
	movapd [esp + nb232nf_rinvH2H1], xmm5

	movapd xmm0, [esp + nb232nf_rsqOO]
	movapd xmm4, [esp + nb232nf_rsqOH1]	
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
	movapd  xmm3, [esp + nb232nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232nf_half] ;# iter1 of  
	mulsd   xmm7, [esp + nb232nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb232nf_half] ;# rinv
	movapd [esp + nb232nf_rinvOO], xmm1
	movapd [esp + nb232nf_rinvOH1], xmm5

	movapd xmm0, [esp + nb232nf_rsqOH2]
	movapd xmm4, [esp + nb232nf_rsqH1O]	
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
	movapd  xmm3, [esp + nb232nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232nf_half] ;# iter1 
	mulsd   xmm7, [esp + nb232nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvOH2], xmm1
	movapd [esp + nb232nf_rinvH1O], xmm5

	movapd xmm0, [esp + nb232nf_rsqH1H1]
	movapd xmm4, [esp + nb232nf_rsqH1H2]	
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
	movapd  xmm3, [esp + nb232nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [esp + nb232nf_half] ;# iter1a 
	mulsd   xmm7, [esp + nb232nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [esp + nb232nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [esp + nb232nf_half] ;# rinv 
	mulsd   xmm5, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvH1H1], xmm1
	movapd [esp + nb232nf_rinvH1H2], xmm5

	movapd xmm0, [esp + nb232nf_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [esp + nb232nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [esp + nb232nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [esp + nb232nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [esp + nb232nf_half] ;# rinv 
	movapd [esp + nb232nf_rinvH2O], xmm1
	
	;# start with OO interaction 	
	movsd xmm0, [esp + nb232nf_rinvOO]
	movsd xmm7, xmm0		;# xmm7=rinv 
	movsd xmm5, [esp + nb232nf_krf]
	mulsd  xmm5, [esp + nb232nf_rsqOO]	;# xmm5=krsq 
	movsd xmm6, xmm5		;# krsq 
	addsd  xmm6, xmm7		;# xmm6=rinv+ krsq 
	subsd  xmm6, [esp + nb232nf_crf]	;# rinv+krsq-crf 	
	mulsd  xmm6, [esp + nb232nf_qqOO]	;# xmm6=voul=qq*(rinv+ krsq-crf) 
	
	addsd  xmm6, [esp + nb232nf_vctot]
	movsd [esp + nb232nf_vctot], xmm6
	
	;# LJ table interaction
	movsd xmm4, [esp + nb232nf_rsqOO]
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb232nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb232nf_VFtab]

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

	movsd xmm4, [esp + nb232nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [esp + nb232nf_Vvdwtot]
	movsd [esp + nb232nf_Vvdwtot], xmm5

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
	
	movsd xmm4, [esp + nb232nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb232nf_Vvdwtot]
	movsd [esp + nb232nf_Vvdwtot], xmm5

	;# O-H1 interaction 
	movapd xmm0, [esp + nb232nf_rinvOH1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqOH1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulsd  xmm0, xmm0
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 

	movsd  xmm6, [esp + nb232nf_vctot]
	addsd  xmm6, xmm4

	;# O-H2 interaction  
	movapd xmm0, [esp + nb232nf_rinvOH2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqOH2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulsd  xmm0, xmm0
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addsd  xmm6, xmm4	;# add to local vctot 

	;# H1-O interaction 
	movapd xmm0, [esp + nb232nf_rinvH1O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqH1O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=rinv+ krsq 
	mulsd xmm0, xmm0
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addsd  xmm6, xmm4	;# add to local vctot 

	;# H1-H1 interaction 
	movapd xmm0, [esp + nb232nf_rinvH1H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqH1H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addsd  xmm6, xmm4	;# add to local vctot 

	;# H1-H2 interaction 
	movapd xmm0, [esp + nb232nf_rinvH1H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqH1H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	mulsd xmm0, xmm0
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addsd  xmm6, xmm4	;# add to local vctot 

	;# H2-O interaction 
	movapd xmm0, [esp + nb232nf_rinvH2O]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqH2O] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232nf_qqOH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addsd  xmm6, xmm4	;# add to local vctot 

	;# H2-H1 interaction 
	movapd xmm0, [esp + nb232nf_rinvH2H1]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqH2H1] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addsd  xmm6, xmm4	;# add to local vctot 

	;# H2-H2 interaction 
	movapd xmm0, [esp + nb232nf_rinvH2H2]
	movapd xmm7, xmm0	;# xmm7=rinv 
	movapd xmm5, [esp + nb232nf_krf]	
	movapd xmm1, xmm0
	mulsd  xmm5, [esp + nb232nf_rsqH2H2] ;# xmm5=krsq 
	movapd xmm4, xmm5
	addsd  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subsd  xmm4, [esp + nb232nf_crf]
	mulsd xmm0, xmm0
	mulsd  xmm4, [esp + nb232nf_qqHH] ;# xmm4=voul=qq*(rinv+ krsq) 
	addsd  xmm6, xmm4	;# add to local vctot 

	movsd [esp + nb232nf_vctot], xmm6
	
.nb232nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb232nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb232nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb232nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb232nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb232nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb232nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb232nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb232nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb232nf_n], esi
        jmp .nb232nf_outer
.nb232nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb232nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb232nf_end
        ;# non-zero, do one more workunit
        jmp   .nb232nf_threadloop
.nb232nf_end:
	emms

	mov eax, [esp + nb232nf_nouter]
	mov ebx, [esp + nb232nf_ninner]
	mov ecx, [ebp + nb232nf_outeriter]
	mov edx, [ebp + nb232nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb232nf_salign]
	add esp, eax
	add esp, 824
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


