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


.globl nb_kernel113_ia32_sse2
.globl _nb_kernel113_ia32_sse2
nb_kernel113_ia32_sse2:	
_nb_kernel113_ia32_sse2:	
.equiv          nb113_p_nri,            8
.equiv          nb113_iinr,             12
.equiv          nb113_jindex,           16
.equiv          nb113_jjnr,             20
.equiv          nb113_shift,            24
.equiv          nb113_shiftvec,         28
.equiv          nb113_fshift,           32
.equiv          nb113_gid,              36
.equiv          nb113_pos,              40
.equiv          nb113_faction,          44
.equiv          nb113_charge,           48
.equiv          nb113_p_facel,          52
.equiv          nb113_argkrf,           56
.equiv          nb113_argcrf,           60
.equiv          nb113_Vc,               64
.equiv          nb113_type,             68
.equiv          nb113_p_ntype,          72
.equiv          nb113_vdwparam,         76
.equiv          nb113_Vvdw,             80
.equiv          nb113_p_tabscale,       84
.equiv          nb113_VFtab,            88
.equiv          nb113_invsqrta,         92
.equiv          nb113_dvda,             96
.equiv          nb113_p_gbtabscale,     100
.equiv          nb113_GBtab,            104
.equiv          nb113_p_nthreads,       108
.equiv          nb113_count,            112
.equiv          nb113_mtx,              116
.equiv          nb113_outeriter,        120
.equiv          nb113_inneriter,        124
.equiv          nb113_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb113_ixO,              0
.equiv          nb113_iyO,              16
.equiv          nb113_izO,              32
.equiv          nb113_ixH1,             48
.equiv          nb113_iyH1,             64
.equiv          nb113_izH1,             80
.equiv          nb113_ixH2,             96
.equiv          nb113_iyH2,             112
.equiv          nb113_izH2,             128
.equiv          nb113_ixM,              144
.equiv          nb113_iyM,              160
.equiv          nb113_izM,              176
.equiv          nb113_iqH,              192
.equiv          nb113_iqM,              208
.equiv          nb113_dxO,              224
.equiv          nb113_dyO,              240
.equiv          nb113_dzO,              256
.equiv          nb113_dxH1,             272
.equiv          nb113_dyH1,             288
.equiv          nb113_dzH1,             304
.equiv          nb113_dxH2,             320
.equiv          nb113_dyH2,             336
.equiv          nb113_dzH2,             352
.equiv          nb113_dxM,              368
.equiv          nb113_dyM,              384
.equiv          nb113_dzM,              400
.equiv          nb113_qqH,              416
.equiv          nb113_qqM,              432
.equiv          nb113_c6,               448
.equiv          nb113_c12,              464
.equiv          nb113_six,              480
.equiv          nb113_twelve,           496
.equiv          nb113_vctot,            512
.equiv          nb113_Vvdwtot,          528
.equiv          nb113_fixO,             544
.equiv          nb113_fiyO,             560
.equiv          nb113_fizO,             576
.equiv          nb113_fixH1,            592
.equiv          nb113_fiyH1,            608
.equiv          nb113_fizH1,            624
.equiv          nb113_fixH2,            640
.equiv          nb113_fiyH2,            656
.equiv          nb113_fizH2,            672
.equiv          nb113_fixM,             688
.equiv          nb113_fiyM,             704
.equiv          nb113_fizM,             720
.equiv          nb113_fjx,              736
.equiv          nb113_fjy,              752
.equiv          nb113_fjz,              768
.equiv          nb113_half,             784
.equiv          nb113_three,            800
.equiv          nb113_two,              816
.equiv          nb113_rinvH1,           832
.equiv          nb113_rinvH2,           848
.equiv          nb113_rinvM,            864
.equiv          nb113_is3,              880
.equiv          nb113_ii3,              884
.equiv          nb113_ntia,             888
.equiv          nb113_innerjjnr,        892
.equiv          nb113_innerk,           896
.equiv          nb113_n,                900
.equiv          nb113_nn1,              904
.equiv          nb113_nri,              908
.equiv          nb113_nouter,           912
.equiv          nb113_ninner,           916
.equiv          nb113_salign,           920
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 924		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb113_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb113_p_nri]
	mov ecx, [ecx]
	mov [esp + nb113_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb113_nouter], eax
	mov [esp + nb113_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000     ;# upper half of 0.5
	mov [esp + nb113_half], eax
	mov [esp + nb113_half+4], ebx
	movsd xmm1, [esp + nb113_half]
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
	movapd [esp + nb113_half], xmm1
	movapd [esp + nb113_two], xmm2
	movapd [esp + nb113_three], xmm3
	movapd [esp + nb113_six], xmm4
	movapd [esp + nb113_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb113_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb113_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb113_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb113_iqH], xmm3
	movapd [esp + nb113_iqM], xmm4
	
	mov   edx, [ebp + nb113_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb113_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb113_ntia], ecx		
.nb113_threadloop:
        mov   esi, [ebp + nb113_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb113_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb113_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb113_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb113_n], eax
        mov [esp + nb113_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb113_outerstart
        jmp .nb113_end

.nb113_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb113_nouter]
	mov [esp + nb113_nouter], ebx

.nb113_outer:
	mov   eax, [ebp + nb113_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb113_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb113_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb113_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb113_pos]    ;# eax = base of pos[]  
	mov   [esp + nb113_ii3], ebx

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
	movapd [esp + nb113_ixO], xmm3
	movapd [esp + nb113_iyO], xmm4
	movapd [esp + nb113_izO], xmm5
	movapd [esp + nb113_ixH1], xmm6
	movapd [esp + nb113_iyH1], xmm7

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
	movapd [esp + nb113_izH1], xmm6
	movapd [esp + nb113_ixH2], xmm0
	movapd [esp + nb113_iyH2], xmm1
	movapd [esp + nb113_izH2], xmm2
	movapd [esp + nb113_ixM], xmm3
	movapd [esp + nb113_iyM], xmm4
	movapd [esp + nb113_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb113_vctot], xmm4
	movapd [esp + nb113_Vvdwtot], xmm4
	movapd [esp + nb113_fixO], xmm4
	movapd [esp + nb113_fiyO], xmm4
	movapd [esp + nb113_fizO], xmm4
	movapd [esp + nb113_fixH1], xmm4
	movapd [esp + nb113_fiyH1], xmm4
	movapd [esp + nb113_fizH1], xmm4
	movapd [esp + nb113_fixH2], xmm4
	movapd [esp + nb113_fiyH2], xmm4
	movapd [esp + nb113_fizH2], xmm4
	movapd [esp + nb113_fixM], xmm4
	movapd [esp + nb113_fiyM], xmm4
	movapd [esp + nb113_fizM], xmm4
	
	mov   eax, [ebp + nb113_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb113_pos]
	mov   edi, [ebp + nb113_faction]	
	mov   eax, [ebp + nb113_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb113_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb113_ninner]
	mov   [esp + nb113_ninner], ecx
	add   edx, 0
	mov   [esp + nb113_innerk], edx    ;# number of innerloop atoms 
	jge   .nb113_unroll_loop
	jmp   .nb113_checksingle
.nb113_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb113_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb113_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb113_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb113_iqM]
	mulpd  xmm4, [esp + nb113_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb113_qqM], xmm3
	movapd  [esp + nb113_qqH], xmm4
	
	mov esi, [ebp + nb113_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb113_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb113_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movlpd xmm7, [esi + ebx*8]	;# c6b
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + ebx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb113_c6], xmm4
	movapd [esp + nb113_c12], xmm6
	
	mov esi, [ebp + nb113_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb113_ixO]
	movapd xmm5, [esp + nb113_iyO]
	movapd xmm6, [esp + nb113_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb113_dxO], xmm4
	movapd [esp + nb113_dyO], xmm5
	movapd [esp + nb113_dzO], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb113_ixH1]
	movapd xmm5, [esp + nb113_iyH1]
	movapd xmm6, [esp + nb113_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb113_dxH1], xmm4
	movapd [esp + nb113_dyH1], xmm5
	movapd [esp + nb113_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb113_ixH2]
	movapd xmm4, [esp + nb113_iyH2]
	movapd xmm5, [esp + nb113_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb113_dxH2], xmm3
	movapd [esp + nb113_dyH2], xmm4
	movapd [esp + nb113_dzH2], xmm5
	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3

	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb113_iyM]
	movapd xmm4, [esp + nb113_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb113_ixM]
	subpd  xmm2, xmm0	

	;# store dr 
	movapd [esp + nb113_dxM], xmm2
	movapd [esp + nb113_dyM], xmm3
	movapd [esp + nb113_dzM], xmm4
	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb113_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb113_half] ;# rinv 
	movapd  [esp + nb113_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb113_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb113_half] ;# rinv 
	movapd  [esp + nb113_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb113_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb113_half] ;# rinv 
	movapd  [esp + nb113_rinvM], xmm1	

	;# do O interactions directly - rsqO is in xmm7
	cvtpd2ps xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtps2pd xmm2, xmm2
	movapd   xmm1, [esp + nb113_two]
	movapd   xmm0, xmm1
	mulpd   xmm7, xmm2
	subpd   xmm1, xmm7
	mulpd   xmm2, xmm1 ;# iter1 
	mulpd   xmm6, xmm2
	subpd   xmm0, xmm6
	mulpd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulpd   xmm1, xmm1 ;# rinv4
	mulpd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulpd	xmm2, xmm2 ;# rinvtwelve
	mulpd  xmm1, [esp + nb113_c6]
	mulpd  xmm2, [esp + nb113_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb113_Vvdwtot]
	mulpd  xmm1, [esp + nb113_six]
	mulpd  xmm2, [esp + nb113_twelve]
	subpd  xmm2, xmm1
	mulpd  xmm2, xmm0
	movapd xmm4, xmm2 ;# total fsO 
	movapd [esp + nb113_Vvdwtot], xmm3

	movapd xmm0, [esp + nb113_dxO]
	movapd xmm1, [esp + nb113_dyO]
	movapd xmm2, [esp + nb113_dzO]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [esp + nb113_fixO]
	movapd xmm4, [esp + nb113_fiyO]
	movapd xmm7, [esp + nb113_fizO]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb113_fixO], xmm3
	movapd [esp + nb113_fiyO], xmm4
	movapd [esp + nb113_fizO], xmm7
	;# update j forces with water O 
	movapd [esp + nb113_fjx], xmm0
	movapd [esp + nb113_fjy], xmm1
	movapd [esp + nb113_fjz], xmm2

	;# H1 interactions
	movapd  xmm6, [esp + nb113_rinvH1] 
	movapd  xmm4, xmm6	
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulpd  xmm6, [esp + nb113_qqH]	;# xmm6=vcoul 
	mulpd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addpd  xmm6, [esp + nb113_vctot]

	movapd xmm0, [esp + nb113_dxH1]
	movapd xmm1, [esp + nb113_dyH1]
	movapd xmm2, [esp + nb113_dzH1]
	movapd [esp + nb113_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb113_fixH1]
	movapd xmm4, [esp + nb113_fiyH1]
	movapd xmm7, [esp + nb113_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb113_fixH1], xmm3
	movapd [esp + nb113_fiyH1], xmm4
	movapd [esp + nb113_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb113_fjx]
	addpd  xmm1, [esp + nb113_fjy]
	addpd  xmm2, [esp + nb113_fjz]
	movapd [esp + nb113_fjx], xmm0
	movapd [esp + nb113_fjy], xmm1
	movapd [esp + nb113_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [esp + nb113_rinvH2] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulpd  xmm5, [esp + nb113_qqH]	;# xmm5=vcoul 
	mulpd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addpd  xmm5, [esp + nb113_vctot]

	movapd xmm0, [esp + nb113_dxH2]
	movapd xmm1, [esp + nb113_dyH2]
	movapd xmm2, [esp + nb113_dzH2]
	movapd [esp + nb113_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb113_fixH2]
	movapd xmm4, [esp + nb113_fiyH2]
	movapd xmm7, [esp + nb113_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb113_fixH2], xmm3
	movapd [esp + nb113_fiyH2], xmm4
	movapd [esp + nb113_fizH2], xmm7
	;# update j forces with water H2
	addpd  xmm0, [esp + nb113_fjx]
	addpd  xmm1, [esp + nb113_fjy]
	addpd  xmm2, [esp + nb113_fjz]
	movapd [esp + nb113_fjx], xmm0
	movapd [esp + nb113_fjy], xmm1
	movapd [esp + nb113_fjz], xmm2

	;# M interactions 
	movapd  xmm5, [esp + nb113_rinvM] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulpd  xmm5, [esp + nb113_qqM]	;# xmm5=vcoul 
	mulpd  xmm4, xmm5		;# total fsM in xmm4 
	
	addpd  xmm5, [esp + nb113_vctot]

	movapd xmm0, [esp + nb113_dxM]
	movapd xmm1, [esp + nb113_dyM]
	movapd xmm2, [esp + nb113_dzM]
	movapd [esp + nb113_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [esp + nb113_fixM]
	movapd xmm4, [esp + nb113_fiyM]
	movapd xmm7, [esp + nb113_fizM]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb113_fixM], xmm3
	movapd [esp + nb113_fiyM], xmm4
	movapd [esp + nb113_fizM], xmm7

	mov edi, [ebp + nb113_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb113_fjx]
	addpd  xmm1, [esp + nb113_fjy]
	addpd  xmm2, [esp + nb113_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	movhpd xmm3, [edi + ebx*8]
	movhpd xmm4, [edi + ebx*8 + 8]
	movhpd xmm5, [edi + ebx*8 + 16]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5
	movhpd [edi + ebx*8], xmm3
	movhpd [edi + ebx*8 + 8], xmm4
	movhpd [edi + ebx*8 + 16], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb113_innerk],  2
	jl   .nb113_checksingle
	jmp  .nb113_unroll_loop
.nb113_checksingle:	
	mov   edx, [esp + nb113_innerk]
	and   edx, 1
	jnz  .nb113_dosingle
	jmp  .nb113_updateouterdata
.nb113_dosingle:
	mov   edx, [esp + nb113_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb113_innerjjnr],  4	

	mov esi, [ebp + nb113_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb113_iqM]
	mulpd  xmm4, [esp + nb113_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb113_qqM], xmm3
	movapd  [esp + nb113_qqH], xmm4
	
	mov esi, [ebp + nb113_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb113_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb113_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb113_c6], xmm4
	movapd [esp + nb113_c12], xmm6
	
	mov esi, [ebp + nb113_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb113_ixO]
	movapd xmm5, [esp + nb113_iyO]
	movapd xmm6, [esp + nb113_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb113_dxO], xmm4
	movapd [esp + nb113_dyO], xmm5
	movapd [esp + nb113_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb113_ixH1]
	movapd xmm5, [esp + nb113_iyH1]
	movapd xmm6, [esp + nb113_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb113_dxH1], xmm4
	movapd [esp + nb113_dyH1], xmm5
	movapd [esp + nb113_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb113_ixH2]
	movapd xmm4, [esp + nb113_iyH2]
	movapd xmm5, [esp + nb113_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb113_dxH2], xmm3
	movapd [esp + nb113_dyH2], xmm4
	movapd [esp + nb113_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb113_iyM]
	movapd xmm4, [esp + nb113_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb113_ixM]
	subpd  xmm2, xmm0	

	;# store dr 
	movapd [esp + nb113_dxM], xmm2
	movapd [esp + nb113_dyM], xmm3
	movapd [esp + nb113_dzM], xmm4
	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb113_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb113_half] ;# rinv 
	movapd [esp + nb113_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb113_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb113_half] ;# rinv 
	movapd [esp + nb113_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb113_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb113_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb113_half] ;# rinv 
	movapd [esp + nb113_rinvM], xmm1

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [esp + nb113_two]
	movapd   xmm0, xmm1
	mulsd   xmm7, xmm2
	subsd   xmm1, xmm7
	mulsd   xmm2, xmm1 ;# iter1 
	mulsd   xmm6, xmm2
	subsd   xmm0, xmm6
	mulsd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [esp + nb113_c6]
	mulsd  xmm2, [esp + nb113_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb113_Vvdwtot]
	mulsd  xmm1, [esp + nb113_six]
	mulsd  xmm2, [esp + nb113_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm2, xmm0
	movapd xmm4, xmm2 ;# total fsO 
	movsd [esp + nb113_Vvdwtot], xmm3

	movapd xmm0, [esp + nb113_dxO]
	movapd xmm1, [esp + nb113_dyO]
	movapd xmm2, [esp + nb113_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [esp + nb113_fixO]
	movapd xmm4, [esp + nb113_fiyO]
	movapd xmm7, [esp + nb113_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb113_fixO], xmm3
	movsd [esp + nb113_fiyO], xmm4
	movsd [esp + nb113_fizO], xmm7
	;# update j forces with water O 
	movsd [esp + nb113_fjx], xmm0
	movsd [esp + nb113_fjy], xmm1
	movsd [esp + nb113_fjz], xmm2

	;# H1 interactions
	movapd  xmm6, [esp + nb113_rinvH1] 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd  xmm6, [esp + nb113_qqH]	;# xmm6=vcoul 
	mulsd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addsd  xmm6, [esp + nb113_vctot]

	movapd xmm0, [esp + nb113_dxH1]
	movapd xmm1, [esp + nb113_dyH1]
	movapd xmm2, [esp + nb113_dzH1]
	movsd [esp + nb113_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb113_fixH1]
	movapd xmm4, [esp + nb113_fiyH1]
	movapd xmm7, [esp + nb113_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb113_fixH1], xmm3
	movsd [esp + nb113_fiyH1], xmm4
	movsd [esp + nb113_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb113_fjx]
	addsd  xmm1, [esp + nb113_fjy]
	addsd  xmm2, [esp + nb113_fjz]
	movsd [esp + nb113_fjx], xmm0
	movsd [esp + nb113_fjy], xmm1
	movsd [esp + nb113_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [esp + nb113_rinvH2] 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [esp + nb113_qqH]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [esp + nb113_vctot]

	movapd xmm0, [esp + nb113_dxH2]
	movapd xmm1, [esp + nb113_dyH2]
	movapd xmm2, [esp + nb113_dzH2]
	movsd [esp + nb113_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb113_fixH2]
	movapd xmm4, [esp + nb113_fiyH2]
	movapd xmm7, [esp + nb113_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb113_fixH2], xmm3
	movsd [esp + nb113_fiyH2], xmm4
	movsd [esp + nb113_fizH2], xmm7
	;# update j forces with water H2 
	addsd  xmm0, [esp + nb113_fjx]
	addsd  xmm1, [esp + nb113_fjy]
	addsd  xmm2, [esp + nb113_fjz]
	movsd [esp + nb113_fjx], xmm0
	movsd [esp + nb113_fjy], xmm1
	movsd [esp + nb113_fjz], xmm2

	;# M interactions 
	movapd  xmm5, [esp + nb113_rinvM] 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [esp + nb113_qqM]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [esp + nb113_vctot]

	movapd xmm0, [esp + nb113_dxM]
	movapd xmm1, [esp + nb113_dyM]
	movapd xmm2, [esp + nb113_dzM]
	movsd [esp + nb113_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [esp + nb113_fixM]
	movapd xmm4, [esp + nb113_fiyM]
	movapd xmm7, [esp + nb113_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb113_fixM], xmm3
	movsd [esp + nb113_fiyM], xmm4
	movsd [esp + nb113_fizM], xmm7

	mov edi, [ebp + nb113_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb113_fjx]
	addsd  xmm1, [esp + nb113_fjy]
	addsd  xmm2, [esp + nb113_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb113_updateouterdata:
	mov   ecx, [esp + nb113_ii3]
	mov   edi, [ebp + nb113_faction]
	mov   esi, [ebp + nb113_fshift]
	mov   edx, [esp + nb113_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb113_fixO]
	movapd xmm1, [esp + nb113_fiyO]
	movapd xmm2, [esp + nb113_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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
	unpcklpd xmm6,xmm1 

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb113_fixH1]
	movapd xmm1, [esp + nb113_fiyH1]
	movapd xmm2, [esp + nb113_fizH1]

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
	movapd xmm0, [esp + nb113_fixH2]
	movapd xmm1, [esp + nb113_fiyH2]
	movapd xmm2, [esp + nb113_fizH2]

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

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb113_fixM]
	movapd xmm1, [esp + nb113_fiyM]
	movapd xmm2, [esp + nb113_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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
	mov esi, [esp + nb113_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb113_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb113_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb113_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb113_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb113_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb113_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb113_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb113_n], esi
        jmp .nb113_outer
.nb113_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb113_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb113_end
        ;# non-zero, do one more workunit
        jmp   .nb113_threadloop
.nb113_end:
	emms

	mov eax, [esp + nb113_nouter]
	mov ebx, [esp + nb113_ninner]
	mov ecx, [ebp + nb113_outeriter]
	mov edx, [ebp + nb113_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb113_salign]
	add esp, eax
	add esp, 924
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret






.globl nb_kernel113nf_ia32_sse2
.globl _nb_kernel113nf_ia32_sse2
nb_kernel113nf_ia32_sse2:	
_nb_kernel113nf_ia32_sse2:	
.equiv          nb113nf_p_nri,          8
.equiv          nb113nf_iinr,           12
.equiv          nb113nf_jindex,         16
.equiv          nb113nf_jjnr,           20
.equiv          nb113nf_shift,          24
.equiv          nb113nf_shiftvec,       28
.equiv          nb113nf_fshift,         32
.equiv          nb113nf_gid,            36
.equiv          nb113nf_pos,            40
.equiv          nb113nf_faction,        44
.equiv          nb113nf_charge,         48
.equiv          nb113nf_p_facel,        52
.equiv          nb113nf_argkrf,         56
.equiv          nb113nf_argcrf,         60
.equiv          nb113nf_Vc,             64
.equiv          nb113nf_type,           68
.equiv          nb113nf_p_ntype,        72
.equiv          nb113nf_vdwparam,       76
.equiv          nb113nf_Vvdw,           80
.equiv          nb113nf_p_tabscale,     84
.equiv          nb113nf_VFtab,          88
.equiv          nb113nf_invsqrta,       92
.equiv          nb113nf_dvda,           96
.equiv          nb113nf_p_gbtabscale,   100
.equiv          nb113nf_GBtab,          104
.equiv          nb113nf_p_nthreads,     108
.equiv          nb113nf_count,          112
.equiv          nb113nf_mtx,            116
.equiv          nb113nf_outeriter,      120
.equiv          nb113nf_inneriter,      124
.equiv          nb113nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb113nf_ixO,            0
.equiv          nb113nf_iyO,            16
.equiv          nb113nf_izO,            32
.equiv          nb113nf_ixH1,           48
.equiv          nb113nf_iyH1,           64
.equiv          nb113nf_izH1,           80
.equiv          nb113nf_ixH2,           96
.equiv          nb113nf_iyH2,           112
.equiv          nb113nf_izH2,           128
.equiv          nb113nf_ixM,            144
.equiv          nb113nf_iyM,            160
.equiv          nb113nf_izM,            176
.equiv          nb113nf_iqH,            192
.equiv          nb113nf_iqM,            208
.equiv          nb113nf_qqH,            224
.equiv          nb113nf_qqM,            240
.equiv          nb113nf_c6,             256
.equiv          nb113nf_c12,            272
.equiv          nb113nf_vctot,          288
.equiv          nb113nf_Vvdwtot,        304
.equiv          nb113nf_half,           320
.equiv          nb113nf_three,          336
.equiv          nb113nf_two,            352
.equiv          nb113nf_is3,            368
.equiv          nb113nf_ii3,            372
.equiv          nb113nf_ntia,           376
.equiv          nb113nf_innerjjnr,      380
.equiv          nb113nf_innerk,         384
.equiv          nb113nf_n,              388
.equiv          nb113nf_nn1,            392
.equiv          nb113nf_nri,            396
.equiv          nb113nf_nouter,         400
.equiv          nb113nf_ninner,         404
.equiv          nb113nf_salign,         408
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 412		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb113nf_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb113nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb113nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb113nf_nouter], eax
	mov [esp + nb113nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb113nf_half], eax
	mov [esp + nb113nf_half+4], ebx
	movsd xmm1, [esp + nb113nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb113nf_half], xmm1
	movapd [esp + nb113nf_two], xmm2
	movapd [esp + nb113nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb113nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb113nf_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb113nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb113nf_iqH], xmm3
	movapd [esp + nb113nf_iqM], xmm4
	
	mov   edx, [ebp + nb113nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb113nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb113nf_ntia], ecx		

.nb113nf_threadloop:
        mov   esi, [ebp + nb113nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb113nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb113nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb113nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb113nf_n], eax
        mov [esp + nb113nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb113nf_outerstart
        jmp .nb113nf_end

.nb113nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb113nf_nouter]
	mov [esp + nb113nf_nouter], ebx

.nb113nf_outer:
	mov   eax, [ebp + nb113nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb113nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb113nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb113nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb113nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb113nf_ii3], ebx

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
	movapd [esp + nb113nf_ixO], xmm3
	movapd [esp + nb113nf_iyO], xmm4
	movapd [esp + nb113nf_izO], xmm5
	movapd [esp + nb113nf_ixH1], xmm6
	movapd [esp + nb113nf_iyH1], xmm7

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
	movapd [esp + nb113nf_izH1], xmm6
	movapd [esp + nb113nf_ixH2], xmm0
	movapd [esp + nb113nf_iyH2], xmm1
	movapd [esp + nb113nf_izH2], xmm2
	movapd [esp + nb113nf_ixM], xmm3
	movapd [esp + nb113nf_iyM], xmm4
	movapd [esp + nb113nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb113nf_vctot], xmm4
	movapd [esp + nb113nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb113nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb113nf_pos]
	mov   edi, [ebp + nb113nf_faction]	
	mov   eax, [ebp + nb113nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb113nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb113nf_ninner]
	mov   [esp + nb113nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb113nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb113nf_unroll_loop
	jmp   .nb113nf_checksingle
.nb113nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb113nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb113nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb113nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb113nf_iqM]
	mulpd  xmm4, [esp + nb113nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb113nf_qqM], xmm3
	movapd  [esp + nb113nf_qqH], xmm4
	
	mov esi, [ebp + nb113nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb113nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb113nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movlpd xmm7, [esi + ebx*8]	;# c6b
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + ebx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb113nf_c6], xmm4
	movapd [esp + nb113nf_c12], xmm6
	
	mov esi, [ebp + nb113nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb113nf_ixO]
	movapd xmm5, [esp + nb113nf_iyO]
	movapd xmm6, [esp + nb113nf_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb113nf_ixH1]
	movapd xmm5, [esp + nb113nf_iyH1]
	movapd xmm6, [esp + nb113nf_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb113nf_ixH2]
	movapd xmm4, [esp + nb113nf_iyH2]
	movapd xmm5, [esp + nb113nf_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3

	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb113nf_iyM]
	movapd xmm4, [esp + nb113nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb113nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb113nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113nf_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb113nf_half] ;# rinv 
	movapd  xmm6, xmm1	;# rinvH1

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb113nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113nf_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb113nf_half] ;# rinv 
	movapd  xmm5, xmm1	;# rinvH2
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb113nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113nf_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb113nf_half] ;# rinv 
	movapd  xmm4, xmm1	;# rinvM

	;# calculate coulomb potentials from rinv.
	addpd   xmm6, xmm5	;# rinvH1+rinvH2
	mulpd	xmm4, [esp + nb113nf_qqM]
	mulpd	xmm6, [esp + nb113nf_qqH]
	addpd   xmm4, xmm6
	addpd   xmm4, [esp + nb113nf_vctot]
	movapd  [esp + nb113nf_vctot], xmm4

	;# do O interactions - rsqO is in xmm7
	cvtpd2ps xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtps2pd xmm2, xmm2
	movapd   xmm1, [esp + nb113nf_two]
	movapd   xmm0, xmm1
	mulpd   xmm7, xmm2
	subpd   xmm1, xmm7
	mulpd   xmm2, xmm1 ;# iter1 
	mulpd   xmm6, xmm2
	subpd   xmm0, xmm6
	mulpd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulpd   xmm1, xmm1 ;# rinv4
	mulpd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulpd	xmm2, xmm2 ;# rinvtwelve
	mulpd  xmm1, [esp + nb113nf_c6]
	mulpd  xmm2, [esp + nb113nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb113nf_Vvdwtot]
	movapd [esp + nb113nf_Vvdwtot], xmm3
	;# should we do one more iteration? 
	sub dword ptr [esp + nb113nf_innerk],  2
	jl   .nb113nf_checksingle
	jmp  .nb113nf_unroll_loop
.nb113nf_checksingle:	
	add dword ptr [esp + nb113nf_innerk],  2
	jnz  .nb113nf_dosingle
	jmp  .nb113nf_updateouterdata
.nb113nf_dosingle:
	mov   edx, [esp + nb113nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb113nf_innerjjnr],  4	

	mov esi, [ebp + nb113nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb113nf_iqM]
	mulpd  xmm4, [esp + nb113nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb113nf_qqM], xmm3
	movapd  [esp + nb113nf_qqH], xmm4
	
	mov esi, [ebp + nb113nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb113nf_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb113nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb113nf_c6], xmm4
	movapd [esp + nb113nf_c12], xmm6
	
	mov esi, [ebp + nb113nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb113nf_ixO]
	movapd xmm5, [esp + nb113nf_iyO]
	movapd xmm6, [esp + nb113nf_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb113nf_ixH1]
	movapd xmm5, [esp + nb113nf_iyH1]
	movapd xmm6, [esp + nb113nf_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb113nf_ixH2]
	movapd xmm4, [esp + nb113nf_iyH2]
	movapd xmm5, [esp + nb113nf_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb113nf_iyM]
	movapd xmm4, [esp + nb113nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb113nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb113nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113nf_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb113nf_half] ;# rinv 
	movapd xmm6, xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb113nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113nf_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb113nf_half] ;# rinv 
	movapd xmm5, xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb113nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb113nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb113nf_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb113nf_half] ;# rinv 
	movapd xmm4, xmm1

	;# Calculate coulomb potential
	addsd  xmm6, xmm5 	;# rinvH1+rinvH2
	mulsd  xmm4, [esp + nb113nf_qqM]
	mulsd  xmm6, [esp + nb113nf_qqH]
	addsd  xmm4, xmm6
	addsd xmm4, [esp + nb113nf_vctot]
	movsd [esp + nb113nf_vctot], xmm4

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [esp + nb113nf_two]
	movapd   xmm0, xmm1
	mulsd   xmm7, xmm2
	subsd   xmm1, xmm7
	mulsd   xmm2, xmm1 ;# iter1 
	mulsd   xmm6, xmm2
	subsd   xmm0, xmm6
	mulsd   xmm0, xmm2 ;# xmm0=rinvsq
	movapd  xmm1, xmm0	
	mulsd   xmm1, xmm1 ;# rinv4
	mulsd   xmm1, xmm0 ;#rinvsix
	movapd  xmm2, xmm1
	mulsd	xmm2, xmm2 ;# rinvtwelve
	mulsd  xmm1, [esp + nb113nf_c6]
	mulsd  xmm2, [esp + nb113nf_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb113nf_Vvdwtot]
	movsd [esp + nb113nf_Vvdwtot], xmm3

.nb113nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb113nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb113nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb113nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb113nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb113nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb113nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb113nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb113nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb113nf_n], esi
        jmp .nb113nf_outer
.nb113nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb113nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb113nf_end
        ;# non-zero, do one more workunit
        jmp   .nb113nf_threadloop
.nb113nf_end:
	emms

	mov eax, [esp + nb113nf_nouter]
	mov ebx, [esp + nb113nf_ninner]
	mov ecx, [ebp + nb113nf_outeriter]
	mov edx, [ebp + nb113nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb113nf_salign]
	add esp, eax
	add esp, 412
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



