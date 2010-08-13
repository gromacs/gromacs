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


.globl nb_kernel233_ia32_sse2
.globl _nb_kernel233_ia32_sse2
nb_kernel233_ia32_sse2:	
_nb_kernel233_ia32_sse2:	
.equiv          nb233_p_nri,            8
.equiv          nb233_iinr,             12
.equiv          nb233_jindex,           16
.equiv          nb233_jjnr,             20
.equiv          nb233_shift,            24
.equiv          nb233_shiftvec,         28
.equiv          nb233_fshift,           32
.equiv          nb233_gid,              36
.equiv          nb233_pos,              40
.equiv          nb233_faction,          44
.equiv          nb233_charge,           48
.equiv          nb233_p_facel,          52
.equiv          nb233_argkrf,           56
.equiv          nb233_argcrf,           60
.equiv          nb233_Vc,               64
.equiv          nb233_type,             68
.equiv          nb233_p_ntype,          72
.equiv          nb233_vdwparam,         76
.equiv          nb233_Vvdw,             80
.equiv          nb233_p_tabscale,       84
.equiv          nb233_VFtab,            88
.equiv          nb233_invsqrta,         92
.equiv          nb233_dvda,             96
.equiv          nb233_p_gbtabscale,     100
.equiv          nb233_GBtab,            104
.equiv          nb233_p_nthreads,       108
.equiv          nb233_count,            112
.equiv          nb233_mtx,              116
.equiv          nb233_outeriter,        120
.equiv          nb233_inneriter,        124
.equiv          nb233_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb233_ixO,              0
.equiv          nb233_iyO,              16
.equiv          nb233_izO,              32
.equiv          nb233_ixH1,             48
.equiv          nb233_iyH1,             64
.equiv          nb233_izH1,             80
.equiv          nb233_ixH2,             96
.equiv          nb233_iyH2,             112
.equiv          nb233_izH2,             128
.equiv          nb233_ixM,              144
.equiv          nb233_iyM,              160
.equiv          nb233_izM,              176
.equiv          nb233_iqH,              192
.equiv          nb233_iqM,              208
.equiv          nb233_dxO,              224
.equiv          nb233_dyO,              240
.equiv          nb233_dzO,              256
.equiv          nb233_dxH1,             272
.equiv          nb233_dyH1,             288
.equiv          nb233_dzH1,             304
.equiv          nb233_dxH2,             320
.equiv          nb233_dyH2,             336
.equiv          nb233_dzH2,             352
.equiv          nb233_dxM,              368
.equiv          nb233_dyM,              384
.equiv          nb233_dzM,              400
.equiv          nb233_qqH,              416
.equiv          nb233_qqM,              432
.equiv          nb233_c6,               448
.equiv          nb233_c12,              464
.equiv          nb233_tsc,              480
.equiv          nb233_fstmp,            496
.equiv          nb233_vctot,            512
.equiv          nb233_Vvdwtot,          528
.equiv          nb233_fixO,             544
.equiv          nb233_fiyO,             560
.equiv          nb233_fizO,             576
.equiv          nb233_fixH1,            592
.equiv          nb233_fiyH1,            608
.equiv          nb233_fizH1,            624
.equiv          nb233_fixH2,            640
.equiv          nb233_fiyH2,            656
.equiv          nb233_fizH2,            672
.equiv          nb233_fixM,             688
.equiv          nb233_fiyM,             704
.equiv          nb233_fizM,             720
.equiv          nb233_fjx,              736
.equiv          nb233_fjy,              752
.equiv          nb233_fjz,              768
.equiv          nb233_half,             784
.equiv          nb233_three,            800
.equiv          nb233_two,              816
.equiv          nb233_rinvH1,           832
.equiv          nb233_rinvH2,           848
.equiv          nb233_rinvM,            864
.equiv          nb233_krsqH1,           880
.equiv          nb233_krsqH2,           896
.equiv          nb233_krsqM,            912
.equiv          nb233_krf,              928
.equiv          nb233_crf,              944
.equiv          nb233_rsqO,             960
.equiv          nb233_is3,              976
.equiv          nb233_ii3,              980
.equiv          nb233_ntia,             984
.equiv          nb233_innerjjnr,        988
.equiv          nb233_innerk,           992
.equiv          nb233_n,                996
.equiv          nb233_nn1,              1000
.equiv          nb233_nri,              1004
.equiv          nb233_nouter,           1008
.equiv          nb233_ninner,           1012
.equiv          nb233_salign,           1016
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1020		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb233_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb233_p_nri]
	mov ecx, [ecx]
	mov [esp + nb233_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb233_nouter], eax
	mov [esp + nb233_ninner], eax

	mov eax, [ebp + nb233_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb233_tsc], xmm3

	mov esi, [ebp + nb233_argkrf]
	mov edi, [ebp + nb233_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb233_krf], xmm5
	movapd [esp + nb233_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb233_half], eax
	mov [esp + nb233_half+4], ebx
	movsd xmm1, [esp + nb233_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb233_half], xmm1
	movapd [esp + nb233_two], xmm2
	movapd [esp + nb233_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb233_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb233_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb233_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb233_iqH], xmm3
	movapd [esp + nb233_iqM], xmm4
	
	mov   edx, [ebp + nb233_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb233_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb233_ntia], ecx		
.nb233_threadloop:
        mov   esi, [ebp + nb233_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb233_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb233_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb233_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb233_n], eax
        mov [esp + nb233_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb233_outerstart
        jmp .nb233_end

.nb233_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb233_nouter]
	mov [esp + nb233_nouter], ebx

.nb233_outer:
	mov   eax, [ebp + nb233_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb233_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb233_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb233_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb233_pos]    ;# eax = base of pos[]  
	mov   [esp + nb233_ii3], ebx

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
	movapd [esp + nb233_ixO], xmm3
	movapd [esp + nb233_iyO], xmm4
	movapd [esp + nb233_izO], xmm5
	movapd [esp + nb233_ixH1], xmm6
	movapd [esp + nb233_iyH1], xmm7

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
	movapd [esp + nb233_izH1], xmm6
	movapd [esp + nb233_ixH2], xmm0
	movapd [esp + nb233_iyH2], xmm1
	movapd [esp + nb233_izH2], xmm2
	movapd [esp + nb233_ixM], xmm3
	movapd [esp + nb233_iyM], xmm4
	movapd [esp + nb233_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb233_vctot], xmm4
	movapd [esp + nb233_Vvdwtot], xmm4
	movapd [esp + nb233_fixO], xmm4
	movapd [esp + nb233_fiyO], xmm4
	movapd [esp + nb233_fizO], xmm4
	movapd [esp + nb233_fixH1], xmm4
	movapd [esp + nb233_fiyH1], xmm4
	movapd [esp + nb233_fizH1], xmm4
	movapd [esp + nb233_fixH2], xmm4
	movapd [esp + nb233_fiyH2], xmm4
	movapd [esp + nb233_fizH2], xmm4
	movapd [esp + nb233_fixM], xmm4
	movapd [esp + nb233_fiyM], xmm4
	movapd [esp + nb233_fizM], xmm4
	
	mov   eax, [ebp + nb233_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb233_pos]
	mov   edi, [ebp + nb233_faction]	
	mov   eax, [ebp + nb233_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb233_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb233_ninner]
	mov   [esp + nb233_ninner], ecx
	add   edx, 0
	mov   [esp + nb233_innerk], edx    ;# number of innerloop atoms 
	jge   .nb233_unroll_loop
	jmp   .nb233_checksingle
.nb233_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb233_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb233_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb233_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb233_iqM]
	mulpd  xmm4, [esp + nb233_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb233_qqM], xmm3
	movapd  [esp + nb233_qqH], xmm4
	
	mov esi, [ebp + nb233_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb233_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb233_ntia]
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
	movapd [esp + nb233_c6], xmm4
	movapd [esp + nb233_c12], xmm6
	
	mov esi, [ebp + nb233_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb233_ixO]
	movapd xmm5, [esp + nb233_iyO]
	movapd xmm6, [esp + nb233_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb233_dxO], xmm4
	movapd [esp + nb233_dyO], xmm5
	movapd [esp + nb233_dzO], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb233_ixH1]
	movapd xmm5, [esp + nb233_iyH1]
	movapd xmm6, [esp + nb233_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb233_dxH1], xmm4
	movapd [esp + nb233_dyH1], xmm5
	movapd [esp + nb233_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb233_ixH2]
	movapd xmm4, [esp + nb233_iyH2]
	movapd xmm5, [esp + nb233_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb233_dxH2], xmm3
	movapd [esp + nb233_dyH2], xmm4
	movapd [esp + nb233_dzH2], xmm5
	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3

	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb233_iyM]
	movapd xmm4, [esp + nb233_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb233_ixM]
	subpd  xmm2, xmm0	

	;# store dr 
	movapd [esp + nb233_dxM], xmm2
	movapd [esp + nb233_dyM], xmm3
	movapd [esp + nb233_dzM], xmm4
	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	movapd [esp + nb233_rsqO], xmm7
	
	;# calculate krsq
	movapd xmm0, [esp + nb233_krf]
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	mulpd xmm0, xmm4  
	mulpd xmm1, xmm5
	mulpd xmm2, xmm6
	movapd [esp + nb233_krsqM], xmm0
	movapd [esp + nb233_krsqH2], xmm1
	movapd [esp + nb233_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb233_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb233_half] ;# rinv 
	movapd  [esp + nb233_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb233_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb233_half] ;# rinv 
	movapd  [esp + nb233_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb233_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb233_half] ;# rinv 
	movapd  [esp + nb233_rinvM], xmm1	

		
	;# rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb233_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb233_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb233_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	
	
	movapd xmm4, [esp + nb233_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulpd xmm4, xmm7	;# xmm4=r 
	mulpd xmm4, [esp + nb233_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb233_VFtab]
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
	mulpd  xmm7, [esp + nb233_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb233_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb233_Vvdwtot]
	xorpd  xmm3, xmm3
	mulpd  xmm7, [esp + nb233_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb233_fstmp], xmm3
	movapd [esp + nb233_Vvdwtot], xmm5

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
	mulpd  xmm7, [esp + nb233_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb233_c12]
	mulpd  xmm7, xmm4 
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb233_Vvdwtot]
	movapd xmm3, [esp + nb233_fstmp]
	mulpd  xmm7, [esp + nb233_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb233_Vvdwtot], xmm5

	mulpd  xmm3, xmm0
		
		
	movapd xmm0, [esp + nb233_dxO]
	movapd xmm1, [esp + nb233_dyO]
	movapd xmm2, [esp + nb233_dzO]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb233_faction]
	mulpd  xmm0, xmm3
	mulpd  xmm1, xmm3
	mulpd  xmm2, xmm3
	
	;# update O forces 
	movapd xmm3, [esp + nb233_fixO]
	movapd xmm4, [esp + nb233_fiyO]
	movapd xmm7, [esp + nb233_fizO]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb233_fixO], xmm3
	movapd [esp + nb233_fiyO], xmm4
	movapd [esp + nb233_fizO], xmm7
	;# update j forces with water O 
	movapd [esp + nb233_fjx], xmm0
	movapd [esp + nb233_fjy], xmm1
	movapd [esp + nb233_fjz], xmm2

	;# H1 interactions 
	movapd  xmm6, [esp + nb233_rinvH1] 
	movapd  xmm4, xmm6
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [esp + nb233_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulpd   xmm0, [esp + nb233_two]
	subpd   xmm6, [esp + nb233_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm6, [esp + nb233_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb233_qqH]
	mulpd  xmm4, xmm7		;# total fsH1 in xmm4 
	addpd  xmm6, [esp + nb233_vctot]

	movapd xmm0, [esp + nb233_dxH1]
	movapd xmm1, [esp + nb233_dyH1]
	movapd xmm2, [esp + nb233_dzH1]
	movapd [esp + nb233_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb233_fixH1]
	movapd xmm4, [esp + nb233_fiyH1]
	movapd xmm7, [esp + nb233_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb233_fixH1], xmm3
	movapd [esp + nb233_fiyH1], xmm4
	movapd [esp + nb233_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb233_fjx]
	addpd  xmm1, [esp + nb233_fjy]
	addpd  xmm2, [esp + nb233_fjz]
	movapd [esp + nb233_fjx], xmm0
	movapd [esp + nb233_fjy], xmm1
	movapd [esp + nb233_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [esp + nb233_rinvH2] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb233_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulpd   xmm0, [esp + nb233_two]
	subpd   xmm5, [esp + nb233_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm5, [esp + nb233_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb233_qqH]
	mulpd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addpd  xmm5, [esp + nb233_vctot]

	movapd xmm0, [esp + nb233_dxH2]
	movapd xmm1, [esp + nb233_dyH2]
	movapd xmm2, [esp + nb233_dzH2]
	movapd [esp + nb233_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb233_fixH2]
	movapd xmm4, [esp + nb233_fiyH2]
	movapd xmm7, [esp + nb233_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb233_fixH2], xmm3
	movapd [esp + nb233_fiyH2], xmm4
	movapd [esp + nb233_fizH2], xmm7
	;# update j forces with water H2
	addpd  xmm0, [esp + nb233_fjx]
	addpd  xmm1, [esp + nb233_fjy]
	addpd  xmm2, [esp + nb233_fjz]
	movapd [esp + nb233_fjx], xmm0
	movapd [esp + nb233_fjy], xmm1
	movapd [esp + nb233_fjz], xmm2

	;# M interactions 
	movapd  xmm5, [esp + nb233_rinvM] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb233_krsqM]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulpd   xmm0, [esp + nb233_two]
	subpd   xmm5, [esp + nb233_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm5, [esp + nb233_qqM] ;# vcoul 
	mulpd   xmm7, [esp + nb233_qqM]
	mulpd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addpd  xmm5, [esp + nb233_vctot]

	movapd xmm0, [esp + nb233_dxM]
	movapd xmm1, [esp + nb233_dyM]
	movapd xmm2, [esp + nb233_dzM]
	movapd [esp + nb233_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb233_fixM]
	movapd xmm4, [esp + nb233_fiyM]
	movapd xmm7, [esp + nb233_fizM]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb233_fixM], xmm3
	movapd [esp + nb233_fiyM], xmm4
	movapd [esp + nb233_fizM], xmm7

	mov edi, [ebp + nb233_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb233_fjx]
	addpd  xmm1, [esp + nb233_fjy]
	addpd  xmm2, [esp + nb233_fjz]
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
	sub dword ptr [esp + nb233_innerk],  2
	jl   .nb233_checksingle
	jmp  .nb233_unroll_loop
.nb233_checksingle:	
	mov   edx, [esp + nb233_innerk]
	and   edx, 1
	jnz  .nb233_dosingle
	jmp  .nb233_updateouterdata
.nb233_dosingle:
	mov   edx, [esp + nb233_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb233_innerjjnr],  4	

	mov esi, [ebp + nb233_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [esp + nb233_iqM]
	mulsd  xmm4, [esp + nb233_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb233_qqM], xmm3
	movapd  [esp + nb233_qqH], xmm4
	
	mov esi, [ebp + nb233_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb233_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb233_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb233_c6], xmm4
	movapd [esp + nb233_c12], xmm6
	
	mov esi, [ebp + nb233_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb233_ixO]
	movapd xmm5, [esp + nb233_iyO]
	movapd xmm6, [esp + nb233_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb233_dxO], xmm4
	movapd [esp + nb233_dyO], xmm5
	movapd [esp + nb233_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [esp + nb233_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb233_ixH1]
	movapd xmm5, [esp + nb233_iyH1]
	movapd xmm6, [esp + nb233_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb233_dxH1], xmm4
	movapd [esp + nb233_dyH1], xmm5
	movapd [esp + nb233_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb233_ixH2]
	movapd xmm4, [esp + nb233_iyH2]
	movapd xmm5, [esp + nb233_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb233_dxH2], xmm3
	movapd [esp + nb233_dyH2], xmm4
	movapd [esp + nb233_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb233_iyM]
	movapd xmm4, [esp + nb233_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb233_ixM]
	subpd  xmm2, xmm0	

	;# store dr 
	movapd [esp + nb233_dxM], xmm2
	movapd [esp + nb233_dyM], xmm3
	movapd [esp + nb233_dzM], xmm4
	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [esp + nb233_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [esp + nb233_krsqM], xmm0
	movsd [esp + nb233_krsqH2], xmm1
	movsd [esp + nb233_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb233_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb233_half] ;# rinv 
	movapd [esp + nb233_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb233_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb233_half] ;# rinv 
	movapd [esp + nb233_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb233_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb233_half] ;# rinv 
	movapd [esp + nb233_rinvM], xmm1

	;# rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movsd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movsd  xmm4, [esp + nb233_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb233_half] ;# iter1 ( new lu) 

	movsd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movsd xmm4, [esp + nb233_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb233_half] ;# rinv 
	movsd  xmm7, xmm4	;# rinvO in xmm7 
	
	movsd xmm4, [esp + nb233_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulsd xmm4, xmm7	;# xmm4=r 
	mulsd xmm4, [esp + nb233_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3

	mov  esi, [ebp + nb233_VFtab]

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
	mulsd  xmm7, [esp + nb233_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb233_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb233_Vvdwtot]
	xorpd  xmm3, xmm3
	mulsd  xmm7, [esp + nb233_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb233_fstmp], xmm3
	movsd [esp + nb233_Vvdwtot], xmm5

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
	mulsd  xmm7, [esp + nb233_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb233_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb233_Vvdwtot]
	movsd xmm3, [esp + nb233_fstmp]
	mulsd  xmm7, [esp + nb233_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb233_Vvdwtot], xmm5

	mulsd  xmm3, xmm0
		
		
	movsd xmm0, [esp + nb233_dxO]
	movsd xmm1, [esp + nb233_dyO]
	movsd xmm2, [esp + nb233_dzO]

	mov    edi, [ebp + nb233_faction]
	mulsd  xmm0, xmm3
	mulsd  xmm1, xmm3
	mulsd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [esp + nb233_fixO]
	movapd xmm4, [esp + nb233_fiyO]
	movapd xmm7, [esp + nb233_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb233_fixO], xmm3
	movsd [esp + nb233_fiyO], xmm4
	movsd [esp + nb233_fizO], xmm7
	;# update j forces with water O 
	movsd [esp + nb233_fjx], xmm0
	movsd [esp + nb233_fjy], xmm1
	movsd [esp + nb233_fjz], xmm2

	;# H1 interactions
	movsd  xmm6, [esp + nb233_rinvH1] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm6
	movsd  xmm0, [esp + nb233_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [esp + nb233_two]
	subsd   xmm6, [esp + nb233_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [esp + nb233_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb233_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [esp + nb233_vctot]

	movapd xmm0, [esp + nb233_dxH1]
	movapd xmm1, [esp + nb233_dyH1]
	movapd xmm2, [esp + nb233_dzH1]
	movsd [esp + nb233_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb233_fixH1]
	movapd xmm4, [esp + nb233_fiyH1]
	movapd xmm7, [esp + nb233_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb233_fixH1], xmm3
	movsd [esp + nb233_fiyH1], xmm4
	movsd [esp + nb233_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb233_fjx]
	addsd  xmm1, [esp + nb233_fjy]
	addsd  xmm2, [esp + nb233_fjz]
	movsd [esp + nb233_fjx], xmm0
	movsd [esp + nb233_fjy], xmm1
	movsd [esp + nb233_fjz], xmm2

	;# H2 interactions 
	movsd  xmm5, [esp + nb233_rinvH2] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [esp + nb233_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [esp + nb233_two]
	subsd   xmm5, [esp + nb233_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [esp + nb233_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb233_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [esp + nb233_vctot]

	movapd xmm0, [esp + nb233_dxH2]
	movapd xmm1, [esp + nb233_dyH2]
	movapd xmm2, [esp + nb233_dzH2]
	movsd [esp + nb233_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb233_fixH2]
	movapd xmm4, [esp + nb233_fiyH2]
	movapd xmm7, [esp + nb233_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb233_fixH2], xmm3
	movsd [esp + nb233_fiyH2], xmm4
	movsd [esp + nb233_fizH2], xmm7
	;# update j forces with water H2 
	addsd  xmm0, [esp + nb233_fjx]
	addsd  xmm1, [esp + nb233_fjy]
	addsd  xmm2, [esp + nb233_fjz]
	movsd [esp + nb233_fjx], xmm0
	movsd [esp + nb233_fjy], xmm1
	movsd [esp + nb233_fjz], xmm2

	;# M interactions 
	movsd  xmm5, [esp + nb233_rinvM] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [esp + nb233_krsqM]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [esp + nb233_two]
	subsd   xmm5, [esp + nb233_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [esp + nb233_qqM] ;# vcoul 
	mulsd   xmm7, [esp + nb233_qqM]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [esp + nb233_vctot]

	movapd xmm0, [esp + nb233_dxM]
	movapd xmm1, [esp + nb233_dyM]
	movapd xmm2, [esp + nb233_dzM]
	movsd [esp + nb233_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [esp + nb233_fixM]
	movapd xmm4, [esp + nb233_fiyM]
	movapd xmm7, [esp + nb233_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb233_fixM], xmm3
	movsd [esp + nb233_fiyM], xmm4
	movsd [esp + nb233_fizM], xmm7

	mov edi, [ebp + nb233_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb233_fjx]
	addsd  xmm1, [esp + nb233_fjy]
	addsd  xmm2, [esp + nb233_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb233_updateouterdata:
	mov   ecx, [esp + nb233_ii3]
	mov   edi, [ebp + nb233_faction]
	mov   esi, [ebp + nb233_fshift]
	mov   edx, [esp + nb233_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb233_fixO]
	movapd xmm1, [esp + nb233_fiyO]
	movapd xmm2, [esp + nb233_fizO]

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
	movapd xmm0, [esp + nb233_fixH1]
	movapd xmm1, [esp + nb233_fiyH1]
	movapd xmm2, [esp + nb233_fizH1]

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
	movapd xmm0, [esp + nb233_fixH2]
	movapd xmm1, [esp + nb233_fiyH2]
	movapd xmm2, [esp + nb233_fizH2]

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
	movapd xmm0, [esp + nb233_fixM]
	movapd xmm1, [esp + nb233_fiyM]
	movapd xmm2, [esp + nb233_fizM]

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
	mov esi, [esp + nb233_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb233_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb233_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb233_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb233_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb233_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb233_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb233_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb233_n], esi
        jmp .nb233_outer
.nb233_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb233_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb233_end
        ;# non-zero, do one more workunit
        jmp   .nb233_threadloop
.nb233_end:
	emms

	mov eax, [esp + nb233_nouter]
	mov ebx, [esp + nb233_ninner]
	mov ecx, [ebp + nb233_outeriter]
	mov edx, [ebp + nb233_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb233_salign]
	add esp, eax
	add esp, 1020
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel233nf_ia32_sse2
.globl _nb_kernel233nf_ia32_sse2
nb_kernel233nf_ia32_sse2:	
_nb_kernel233nf_ia32_sse2:	
.equiv          nb233nf_p_nri,            8
.equiv          nb233nf_iinr,             12
.equiv          nb233nf_jindex,           16
.equiv          nb233nf_jjnr,             20
.equiv          nb233nf_shift,            24
.equiv          nb233nf_shiftvec,         28
.equiv          nb233nf_fshift,           32
.equiv          nb233nf_gid,              36
.equiv          nb233nf_pos,              40
.equiv          nb233nf_faction,          44
.equiv          nb233nf_charge,           48
.equiv          nb233nf_p_facel,          52
.equiv          nb233nf_argkrf,           56
.equiv          nb233nf_argcrf,           60
.equiv          nb233nf_Vc,               64
.equiv          nb233nf_type,             68
.equiv          nb233nf_p_ntype,          72
.equiv          nb233nf_vdwparam,         76
.equiv          nb233nf_Vvdw,             80
.equiv          nb233nf_p_tabscale,       84
.equiv          nb233nf_VFtab,            88
.equiv          nb233nf_invsqrta,         92
.equiv          nb233nf_dvda,             96
.equiv          nb233nf_p_gbtabscale,     100
.equiv          nb233nf_GBtab,            104
.equiv          nb233nf_p_nthreads,       108
.equiv          nb233nf_count,            112
.equiv          nb233nf_mtx,              116
.equiv          nb233nf_outeriter,        120
.equiv          nb233nf_inneriter,        124
.equiv          nb233nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb233nf_ixO,              0
.equiv          nb233nf_iyO,              16
.equiv          nb233nf_izO,              32
.equiv          nb233nf_ixH1,             48
.equiv          nb233nf_iyH1,             64
.equiv          nb233nf_izH1,             80
.equiv          nb233nf_ixH2,             96
.equiv          nb233nf_iyH2,             112
.equiv          nb233nf_izH2,             128
.equiv          nb233nf_ixM,              144
.equiv          nb233nf_iyM,              160
.equiv          nb233nf_izM,              176
.equiv          nb233nf_iqH,              192
.equiv          nb233nf_iqM,              208
.equiv          nb233nf_qqH,              416
.equiv          nb233nf_qqM,              432
.equiv          nb233nf_c6,               448
.equiv          nb233nf_c12,              464
.equiv          nb233nf_tsc,              480
.equiv          nb233nf_vctot,            512
.equiv          nb233nf_Vvdwtot,          528
.equiv          nb233nf_half,             784
.equiv          nb233nf_three,            800
.equiv          nb233nf_two,              816
.equiv          nb233nf_rinvH1,           832
.equiv          nb233nf_rinvH2,           848
.equiv          nb233nf_rinvM,            864
.equiv          nb233nf_krsqH1,           880
.equiv          nb233nf_krsqH2,           896
.equiv          nb233nf_krsqM,            912
.equiv          nb233nf_krf,              928
.equiv          nb233nf_crf,              944
.equiv          nb233nf_rsqO,             960
.equiv          nb233nf_is3,              976
.equiv          nb233nf_ii3,              980
.equiv          nb233nf_ntia,             984
.equiv          nb233nf_innerjjnr,        988
.equiv          nb233nf_innerk,           992
.equiv          nb233nf_n,                996
.equiv          nb233nf_nn1,              1000
.equiv          nb233nf_nri,              1004
.equiv          nb233nf_nouter,           1008
.equiv          nb233nf_ninner,           1012
.equiv          nb233nf_salign,           1016
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1020		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb233nf_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb233nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb233nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb233nf_nouter], eax
	mov [esp + nb233nf_ninner], eax

	mov eax, [ebp + nb233nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb233nf_tsc], xmm3

	mov esi, [ebp + nb233nf_argkrf]
	mov edi, [ebp + nb233nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb233nf_krf], xmm5
	movapd [esp + nb233nf_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb233nf_half], eax
	mov [esp + nb233nf_half+4], ebx
	movsd xmm1, [esp + nb233nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb233nf_half], xmm1
	movapd [esp + nb233nf_two], xmm2
	movapd [esp + nb233nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb233nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb233nf_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb233nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb233nf_iqH], xmm3
	movapd [esp + nb233nf_iqM], xmm4
	
	mov   edx, [ebp + nb233nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb233nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb233nf_ntia], ecx		
.nb233nf_threadloop:
        mov   esi, [ebp + nb233nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb233nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb233nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb233nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb233nf_n], eax
        mov [esp + nb233nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb233nf_outerstart
        jmp .nb233nf_end

.nb233nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb233nf_nouter]
	mov [esp + nb233nf_nouter], ebx

.nb233nf_outer:
	mov   eax, [ebp + nb233nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb233nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb233nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb233nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb233nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb233nf_ii3], ebx

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
	movapd [esp + nb233nf_ixO], xmm3
	movapd [esp + nb233nf_iyO], xmm4
	movapd [esp + nb233nf_izO], xmm5
	movapd [esp + nb233nf_ixH1], xmm6
	movapd [esp + nb233nf_iyH1], xmm7

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
	movapd [esp + nb233nf_izH1], xmm6
	movapd [esp + nb233nf_ixH2], xmm0
	movapd [esp + nb233nf_iyH2], xmm1
	movapd [esp + nb233nf_izH2], xmm2
	movapd [esp + nb233nf_ixM], xmm3
	movapd [esp + nb233nf_iyM], xmm4
	movapd [esp + nb233nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb233nf_vctot], xmm4
	movapd [esp + nb233nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb233nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb233nf_pos]
	mov   edi, [ebp + nb233nf_faction]	
	mov   eax, [ebp + nb233nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb233nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb233nf_ninner]
	mov   [esp + nb233nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb233nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb233nf_unroll_loop
	jmp   .nb233nf_checksingle
.nb233nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb233nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb233nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb233nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb233nf_iqM]
	mulpd  xmm4, [esp + nb233nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb233nf_qqM], xmm3
	movapd  [esp + nb233nf_qqH], xmm4
	
	mov esi, [ebp + nb233nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb233nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb233nf_ntia]
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
	movapd [esp + nb233nf_c6], xmm4
	movapd [esp + nb233nf_c12], xmm6
	
	mov esi, [ebp + nb233nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb233nf_ixO]
	movapd xmm5, [esp + nb233nf_iyO]
	movapd xmm6, [esp + nb233nf_izO]

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
	movapd xmm4, [esp + nb233nf_ixH1]
	movapd xmm5, [esp + nb233nf_iyH1]
	movapd xmm6, [esp + nb233nf_izH1]

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
	movapd xmm3, [esp + nb233nf_ixH2]
	movapd xmm4, [esp + nb233nf_iyH2]
	movapd xmm5, [esp + nb233nf_izH2]

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
	movapd xmm3, [esp + nb233nf_iyM]
	movapd xmm4, [esp + nb233nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb233nf_ixM]
	subpd  xmm2, xmm0	


	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	movapd [esp + nb233nf_rsqO], xmm7
	
	;# calculate krsq
	movapd xmm0, [esp + nb233nf_krf]
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	mulpd xmm0, xmm4  
	mulpd xmm1, xmm5
	mulpd xmm2, xmm6
	movapd [esp + nb233nf_krsqM], xmm0
	movapd [esp + nb233nf_krsqH2], xmm1
	movapd [esp + nb233nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb233nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233nf_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb233nf_half] ;# rinv 
	movapd  [esp + nb233nf_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb233nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233nf_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb233nf_half] ;# rinv 
	movapd  [esp + nb233nf_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb233nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233nf_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb233nf_half] ;# rinv 
	movapd  [esp + nb233nf_rinvM], xmm1	

		
	;# rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb233nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb233nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb233nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	
	
	movapd xmm4, [esp + nb233nf_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulpd xmm4, xmm7	;# xmm4=r 
	mulpd xmm4, [esp + nb233nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb233nf_VFtab]
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

	movapd xmm4, [esp + nb233nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [esp + nb233nf_Vvdwtot]
	movapd [esp + nb233nf_Vvdwtot], xmm5

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
	
	movapd xmm4, [esp + nb233nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb233nf_Vvdwtot]
	movapd [esp + nb233nf_Vvdwtot], xmm5

	;# H1 interactions 
	movapd  xmm6, [esp + nb233nf_rinvH1] 
	movapd  xmm4, xmm6
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [esp + nb233nf_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subpd   xmm6, [esp + nb233nf_crf]
	mulpd   xmm6, [esp + nb233nf_qqH] ;# vcoul 
	addpd  xmm6, [esp + nb233nf_vctot]
	movapd [esp + nb233nf_vctot], xmm6

	;# H2 interactions 
	movapd  xmm5, [esp + nb233nf_rinvH2] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb233nf_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [esp + nb233nf_crf]
	mulpd   xmm5, [esp + nb233nf_qqH] ;# vcoul 
	addpd  xmm5, [esp + nb233nf_vctot]
	movapd [esp + nb233nf_vctot], xmm5

	;# M interactions 
	movapd  xmm5, [esp + nb233nf_rinvM] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb233nf_krsqM]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [esp + nb233nf_crf]
	mulpd   xmm5, [esp + nb233nf_qqM] ;# vcoul 
	addpd  xmm5, [esp + nb233nf_vctot]
	movapd [esp + nb233nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb233nf_innerk],  2
	jl   .nb233nf_checksingle
	jmp  .nb233nf_unroll_loop
.nb233nf_checksingle:	
	mov   edx, [esp + nb233nf_innerk]
	and   edx, 1
	jnz  .nb233nf_dosingle
	jmp  .nb233nf_updateouterdata
.nb233nf_dosingle:
	mov   edx, [esp + nb233nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb233nf_innerjjnr],  4	

	mov esi, [ebp + nb233nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [esp + nb233nf_iqM]
	mulsd  xmm4, [esp + nb233nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb233nf_qqM], xmm3
	movapd  [esp + nb233nf_qqH], xmm4
	
	mov esi, [ebp + nb233nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb233nf_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb233nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb233nf_c6], xmm4
	movapd [esp + nb233nf_c12], xmm6
	
	mov esi, [ebp + nb233nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb233nf_ixO]
	movapd xmm5, [esp + nb233nf_iyO]
	movapd xmm6, [esp + nb233nf_izO]

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
	movapd [esp + nb233nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb233nf_ixH1]
	movapd xmm5, [esp + nb233nf_iyH1]
	movapd xmm6, [esp + nb233nf_izH1]

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
	movapd xmm3, [esp + nb233nf_ixH2]
	movapd xmm4, [esp + nb233nf_iyH2]
	movapd xmm5, [esp + nb233nf_izH2]

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
	movapd xmm3, [esp + nb233nf_iyM]
	movapd xmm4, [esp + nb233nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb233nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [esp + nb233nf_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [esp + nb233nf_krsqM], xmm0
	movsd [esp + nb233nf_krsqH2], xmm1
	movsd [esp + nb233nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb233nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233nf_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb233nf_half] ;# rinv 
	movapd [esp + nb233nf_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb233nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233nf_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb233nf_half] ;# rinv 
	movapd [esp + nb233nf_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb233nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb233nf_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb233nf_half] ;# rinv 
	movapd [esp + nb233nf_rinvM], xmm1

	;# rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movsd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movsd  xmm4, [esp + nb233nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb233nf_half] ;# iter1 ( new lu) 

	movsd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movsd xmm4, [esp + nb233nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb233nf_half] ;# rinv 
	movsd  xmm7, xmm4	;# rinvO in xmm7 
	
	movsd xmm4, [esp + nb233nf_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulsd xmm4, xmm7	;# xmm4=r 
	mulsd xmm4, [esp + nb233nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3

	mov  esi, [ebp + nb233nf_VFtab]

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

	movsd xmm4, [esp + nb233nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [esp + nb233nf_Vvdwtot]
	movsd [esp + nb233nf_Vvdwtot], xmm5

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
	
	movsd xmm4, [esp + nb233nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb233nf_Vvdwtot]
	movsd [esp + nb233nf_Vvdwtot], xmm5

	;# H1 interactions
	movsd  xmm6, [esp + nb233nf_rinvH1] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm6
	movsd  xmm0, [esp + nb233nf_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subsd   xmm6, [esp + nb233nf_crf]
	mulsd   xmm6, [esp + nb233nf_qqH] ;# vcoul 
	addsd  xmm6, [esp + nb233nf_vctot]
	movsd [esp + nb233nf_vctot], xmm6

	;# H2 interactions 
	movsd  xmm5, [esp + nb233nf_rinvH2] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [esp + nb233nf_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [esp + nb233nf_crf]
	mulsd   xmm5, [esp + nb233nf_qqH] ;# vcoul 
	addsd  xmm5, [esp + nb233nf_vctot]
	movsd [esp + nb233nf_vctot], xmm5

	;# M interactions 
	movsd  xmm5, [esp + nb233nf_rinvM] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [esp + nb233nf_krsqM]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [esp + nb233nf_crf]
	mulsd   xmm5, [esp + nb233nf_qqM] ;# vcoul 
	addsd  xmm5, [esp + nb233nf_vctot]
	movsd [esp + nb233nf_vctot], xmm5

.nb233nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb233nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb233nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb233nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb233nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb233nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb233nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb233nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb233nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb233nf_n], esi
        jmp .nb233nf_outer
.nb233nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb233nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb233nf_end
        ;# non-zero, do one more workunit
        jmp   .nb233nf_threadloop
.nb233nf_end:
	emms

	mov eax, [esp + nb233nf_nouter]
	mov ebx, [esp + nb233nf_ninner]
	mov ecx, [ebp + nb233nf_outeriter]
	mov edx, [ebp + nb233nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb233nf_salign]
	add esp, eax
	add esp, 1020
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

