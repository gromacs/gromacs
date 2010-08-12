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


.globl nb_kernel213_ia32_sse2
.globl _nb_kernel213_ia32_sse2
nb_kernel213_ia32_sse2:	
_nb_kernel213_ia32_sse2:	
.equiv          nb213_p_nri,            8
.equiv          nb213_iinr,             12
.equiv          nb213_jindex,           16
.equiv          nb213_jjnr,             20
.equiv          nb213_shift,            24
.equiv          nb213_shiftvec,         28
.equiv          nb213_fshift,           32
.equiv          nb213_gid,              36
.equiv          nb213_pos,              40
.equiv          nb213_faction,          44
.equiv          nb213_charge,           48
.equiv          nb213_p_facel,          52
.equiv          nb213_argkrf,           56
.equiv          nb213_argcrf,           60
.equiv          nb213_Vc,               64
.equiv          nb213_type,             68
.equiv          nb213_p_ntype,          72
.equiv          nb213_vdwparam,         76
.equiv          nb213_Vvdw,             80
.equiv          nb213_p_tabscale,       84
.equiv          nb213_VFtab,            88
.equiv          nb213_invsqrta,         92
.equiv          nb213_dvda,             96
.equiv          nb213_p_gbtabscale,     100
.equiv          nb213_GBtab,            104
.equiv          nb213_p_nthreads,       108
.equiv          nb213_count,            112
.equiv          nb213_mtx,              116
.equiv          nb213_outeriter,        120
.equiv          nb213_inneriter,        124
.equiv          nb213_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb213_ixO,              0
.equiv          nb213_iyO,              16
.equiv          nb213_izO,              32
.equiv          nb213_ixH1,             48
.equiv          nb213_iyH1,             64
.equiv          nb213_izH1,             80
.equiv          nb213_ixH2,             96
.equiv          nb213_iyH2,             112
.equiv          nb213_izH2,             128
.equiv          nb213_ixM,              144
.equiv          nb213_iyM,              160
.equiv          nb213_izM,              176
.equiv          nb213_iqH,              192
.equiv          nb213_iqM,              208
.equiv          nb213_dxO,              224
.equiv          nb213_dyO,              240
.equiv          nb213_dzO,              256
.equiv          nb213_dxH1,             272
.equiv          nb213_dyH1,             288
.equiv          nb213_dzH1,             304
.equiv          nb213_dxH2,             320
.equiv          nb213_dyH2,             336
.equiv          nb213_dzH2,             352
.equiv          nb213_dxM,              368
.equiv          nb213_dyM,              384
.equiv          nb213_dzM,              400
.equiv          nb213_qqH,              416
.equiv          nb213_qqM,              432
.equiv          nb213_c6,               448
.equiv          nb213_c12,              464
.equiv          nb213_six,              480
.equiv          nb213_twelve,           496
.equiv          nb213_vctot,            512
.equiv          nb213_Vvdwtot,          528
.equiv          nb213_fixO,             544
.equiv          nb213_fiyO,             560
.equiv          nb213_fizO,             576
.equiv          nb213_fixH1,            592
.equiv          nb213_fiyH1,            608
.equiv          nb213_fizH1,            624
.equiv          nb213_fixH2,            640
.equiv          nb213_fiyH2,            656
.equiv          nb213_fizH2,            672
.equiv          nb213_fixM,             688
.equiv          nb213_fiyM,             704
.equiv          nb213_fizM,             720
.equiv          nb213_fjx,              736
.equiv          nb213_fjy,              752
.equiv          nb213_fjz,              768
.equiv          nb213_half,             784
.equiv          nb213_three,            800
.equiv          nb213_two,              816
.equiv          nb213_rinvH1,           832
.equiv          nb213_rinvH2,           848
.equiv          nb213_rinvM,            864
.equiv          nb213_krsqH1,           880
.equiv          nb213_krsqH2,           896
.equiv          nb213_krsqM,            912
.equiv          nb213_krf,              928
.equiv          nb213_crf,              944
.equiv          nb213_is3,              960
.equiv          nb213_ii3,              964
.equiv          nb213_ntia,             968
.equiv          nb213_innerjjnr,        972
.equiv          nb213_innerk,           976
.equiv          nb213_n,                980
.equiv          nb213_nn1,              984
.equiv          nb213_nri,              988
.equiv          nb213_nouter,           992
.equiv          nb213_ninner,           996
.equiv          nb213_salign,           1000
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1004		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb213_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb213_p_nri]
	mov ecx, [ecx]
	mov [esp + nb213_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb213_nouter], eax
	mov [esp + nb213_ninner], eax


	mov esi, [ebp + nb213_argkrf]
	mov edi, [ebp + nb213_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb213_krf], xmm5
	movapd [esp + nb213_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb213_half], eax
	mov [esp + nb213_half+4], ebx
	movsd xmm1, [esp + nb213_half]
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
	movapd [esp + nb213_half], xmm1
	movapd [esp + nb213_two], xmm2
	movapd [esp + nb213_three], xmm3
	movapd [esp + nb213_six], xmm4
	movapd [esp + nb213_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb213_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb213_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb213_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb213_iqH], xmm3
	movapd [esp + nb213_iqM], xmm4
	
	mov   edx, [ebp + nb213_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb213_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb213_ntia], ecx		
.nb213_threadloop:
        mov   esi, [ebp + nb213_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb213_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb213_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb213_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb213_n], eax
        mov [esp + nb213_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb213_outerstart
        jmp .nb213_end

.nb213_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb213_nouter]
	mov [esp + nb213_nouter], ebx

.nb213_outer:
	mov   eax, [ebp + nb213_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb213_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb213_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb213_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb213_pos]    ;# eax = base of pos[]  
	mov   [esp + nb213_ii3], ebx

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
	movapd [esp + nb213_ixO], xmm3
	movapd [esp + nb213_iyO], xmm4
	movapd [esp + nb213_izO], xmm5
	movapd [esp + nb213_ixH1], xmm6
	movapd [esp + nb213_iyH1], xmm7

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
	movapd [esp + nb213_izH1], xmm6
	movapd [esp + nb213_ixH2], xmm0
	movapd [esp + nb213_iyH2], xmm1
	movapd [esp + nb213_izH2], xmm2
	movapd [esp + nb213_ixM], xmm3
	movapd [esp + nb213_iyM], xmm4
	movapd [esp + nb213_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb213_vctot], xmm4
	movapd [esp + nb213_Vvdwtot], xmm4
	movapd [esp + nb213_fixO], xmm4
	movapd [esp + nb213_fiyO], xmm4
	movapd [esp + nb213_fizO], xmm4
	movapd [esp + nb213_fixH1], xmm4
	movapd [esp + nb213_fiyH1], xmm4
	movapd [esp + nb213_fizH1], xmm4
	movapd [esp + nb213_fixH2], xmm4
	movapd [esp + nb213_fiyH2], xmm4
	movapd [esp + nb213_fizH2], xmm4
	movapd [esp + nb213_fixM], xmm4
	movapd [esp + nb213_fiyM], xmm4
	movapd [esp + nb213_fizM], xmm4
	
	mov   eax, [ebp + nb213_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb213_pos]
	mov   edi, [ebp + nb213_faction]	
	mov   eax, [ebp + nb213_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb213_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb213_ninner]
	mov   [esp + nb213_ninner], ecx
	add   edx, 0
	mov   [esp + nb213_innerk], edx    ;# number of innerloop atoms 
	jge   .nb213_unroll_loop
	jmp   .nb213_checksingle
.nb213_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb213_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb213_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb213_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb213_iqM]
	mulpd  xmm4, [esp + nb213_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb213_qqM], xmm3
	movapd  [esp + nb213_qqH], xmm4
	
	mov esi, [ebp + nb213_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb213_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb213_ntia]
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
	movapd [esp + nb213_c6], xmm4
	movapd [esp + nb213_c12], xmm6
	
	mov esi, [ebp + nb213_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb213_ixO]
	movapd xmm5, [esp + nb213_iyO]
	movapd xmm6, [esp + nb213_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb213_dxO], xmm4
	movapd [esp + nb213_dyO], xmm5
	movapd [esp + nb213_dzO], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb213_ixH1]
	movapd xmm5, [esp + nb213_iyH1]
	movapd xmm6, [esp + nb213_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb213_dxH1], xmm4
	movapd [esp + nb213_dyH1], xmm5
	movapd [esp + nb213_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb213_ixH2]
	movapd xmm4, [esp + nb213_iyH2]
	movapd xmm5, [esp + nb213_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb213_dxH2], xmm3
	movapd [esp + nb213_dyH2], xmm4
	movapd [esp + nb213_dzH2], xmm5
	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3

	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb213_iyM]
	movapd xmm4, [esp + nb213_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb213_ixM]
	subpd  xmm2, xmm0	

	;# store dr 
	movapd [esp + nb213_dxM], xmm2
	movapd [esp + nb213_dyM], xmm3
	movapd [esp + nb213_dzM], xmm4
	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movapd xmm0, [esp + nb213_krf]
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	mulpd xmm0, xmm4  
	mulpd xmm1, xmm5
	mulpd xmm2, xmm6
	movapd [esp + nb213_krsqM], xmm0
	movapd [esp + nb213_krsqH2], xmm1
	movapd [esp + nb213_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb213_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb213_half] ;# rinv 
	movapd  [esp + nb213_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb213_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb213_half] ;# rinv 
	movapd  [esp + nb213_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb213_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb213_half] ;# rinv 
	movapd  [esp + nb213_rinvM], xmm1	

	;# do O interactions directly - rsqO is in xmm7
	cvtpd2ps xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtps2pd xmm2, xmm2
	movapd   xmm1, [esp + nb213_two]
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
	mulpd  xmm1, [esp + nb213_c6]
	mulpd  xmm2, [esp + nb213_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb213_Vvdwtot]
	mulpd  xmm1, [esp + nb213_six]
	mulpd  xmm2, [esp + nb213_twelve]
	subpd  xmm2, xmm1
	mulpd  xmm2, xmm0
	movapd xmm4, xmm2 ;# total fsO 
	movapd [esp + nb213_Vvdwtot], xmm3

	movapd xmm0, [esp + nb213_dxO]
	movapd xmm1, [esp + nb213_dyO]
	movapd xmm2, [esp + nb213_dzO]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [esp + nb213_fixO]
	movapd xmm4, [esp + nb213_fiyO]
	movapd xmm7, [esp + nb213_fizO]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb213_fixO], xmm3
	movapd [esp + nb213_fiyO], xmm4
	movapd [esp + nb213_fizO], xmm7
	;# update j forces with water O 
	movapd [esp + nb213_fjx], xmm0
	movapd [esp + nb213_fjy], xmm1
	movapd [esp + nb213_fjz], xmm2

	;# H1 interactions 
	movapd  xmm6, [esp + nb213_rinvH1] 
	movapd  xmm4, xmm6
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [esp + nb213_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulpd   xmm0, [esp + nb213_two]
	subpd   xmm6, [esp + nb213_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm6, [esp + nb213_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb213_qqH]
	mulpd  xmm4, xmm7		;# total fsH1 in xmm4 
	addpd  xmm6, [esp + nb213_vctot]

	movapd xmm0, [esp + nb213_dxH1]
	movapd xmm1, [esp + nb213_dyH1]
	movapd xmm2, [esp + nb213_dzH1]
	movapd [esp + nb213_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb213_fixH1]
	movapd xmm4, [esp + nb213_fiyH1]
	movapd xmm7, [esp + nb213_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb213_fixH1], xmm3
	movapd [esp + nb213_fiyH1], xmm4
	movapd [esp + nb213_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb213_fjx]
	addpd  xmm1, [esp + nb213_fjy]
	addpd  xmm2, [esp + nb213_fjz]
	movapd [esp + nb213_fjx], xmm0
	movapd [esp + nb213_fjy], xmm1
	movapd [esp + nb213_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [esp + nb213_rinvH2] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb213_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulpd   xmm0, [esp + nb213_two]
	subpd   xmm5, [esp + nb213_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm5, [esp + nb213_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb213_qqH]
	mulpd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addpd  xmm5, [esp + nb213_vctot]

	movapd xmm0, [esp + nb213_dxH2]
	movapd xmm1, [esp + nb213_dyH2]
	movapd xmm2, [esp + nb213_dzH2]
	movapd [esp + nb213_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb213_fixH2]
	movapd xmm4, [esp + nb213_fiyH2]
	movapd xmm7, [esp + nb213_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb213_fixH2], xmm3
	movapd [esp + nb213_fiyH2], xmm4
	movapd [esp + nb213_fizH2], xmm7
	;# update j forces with water H2
	addpd  xmm0, [esp + nb213_fjx]
	addpd  xmm1, [esp + nb213_fjy]
	addpd  xmm2, [esp + nb213_fjz]
	movapd [esp + nb213_fjx], xmm0
	movapd [esp + nb213_fjy], xmm1
	movapd [esp + nb213_fjz], xmm2

	;# M interactions 
	movapd  xmm5, [esp + nb213_rinvM] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb213_krsqM]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulpd   xmm0, [esp + nb213_two]
	subpd   xmm5, [esp + nb213_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm5, [esp + nb213_qqM] ;# vcoul 
	mulpd   xmm7, [esp + nb213_qqM]
	mulpd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addpd  xmm5, [esp + nb213_vctot]

	movapd xmm0, [esp + nb213_dxM]
	movapd xmm1, [esp + nb213_dyM]
	movapd xmm2, [esp + nb213_dzM]
	movapd [esp + nb213_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb213_fixM]
	movapd xmm4, [esp + nb213_fiyM]
	movapd xmm7, [esp + nb213_fizM]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb213_fixM], xmm3
	movapd [esp + nb213_fiyM], xmm4
	movapd [esp + nb213_fizM], xmm7

	mov edi, [ebp + nb213_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb213_fjx]
	addpd  xmm1, [esp + nb213_fjy]
	addpd  xmm2, [esp + nb213_fjz]
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
	sub dword ptr [esp + nb213_innerk],  2
	jl   .nb213_checksingle
	jmp  .nb213_unroll_loop
.nb213_checksingle:	
	mov   edx, [esp + nb213_innerk]
	and   edx, 1
	jnz  .nb213_dosingle
	jmp  .nb213_updateouterdata
.nb213_dosingle:
	mov   edx, [esp + nb213_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb213_innerjjnr],  4	

	mov esi, [ebp + nb213_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [esp + nb213_iqM]
	mulsd  xmm4, [esp + nb213_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb213_qqM], xmm3
	movapd  [esp + nb213_qqH], xmm4
	
	mov esi, [ebp + nb213_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb213_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb213_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb213_c6], xmm4
	movapd [esp + nb213_c12], xmm6
	
	mov esi, [ebp + nb213_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb213_ixO]
	movapd xmm5, [esp + nb213_iyO]
	movapd xmm6, [esp + nb213_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb213_dxO], xmm4
	movapd [esp + nb213_dyO], xmm5
	movapd [esp + nb213_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb213_ixH1]
	movapd xmm5, [esp + nb213_iyH1]
	movapd xmm6, [esp + nb213_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb213_dxH1], xmm4
	movapd [esp + nb213_dyH1], xmm5
	movapd [esp + nb213_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb213_ixH2]
	movapd xmm4, [esp + nb213_iyH2]
	movapd xmm5, [esp + nb213_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb213_dxH2], xmm3
	movapd [esp + nb213_dyH2], xmm4
	movapd [esp + nb213_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# move ixM-izM to xmm2-xmm4  
	movapd xmm3, [esp + nb213_iyM]
	movapd xmm4, [esp + nb213_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb213_ixM]
	subpd  xmm2, xmm0	

	;# store dr 
	movapd [esp + nb213_dxM], xmm2
	movapd [esp + nb213_dyM], xmm3
	movapd [esp + nb213_dzM], xmm4
	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [esp + nb213_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [esp + nb213_krsqM], xmm0
	movsd [esp + nb213_krsqH2], xmm1
	movsd [esp + nb213_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb213_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb213_half] ;# rinv 
	movapd [esp + nb213_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb213_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb213_half] ;# rinv 
	movapd [esp + nb213_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb213_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb213_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb213_half] ;# rinv 
	movapd [esp + nb213_rinvM], xmm1

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [esp + nb213_two]
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
	mulsd  xmm1, [esp + nb213_c6]
	mulsd  xmm2, [esp + nb213_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb213_Vvdwtot]
	mulsd  xmm1, [esp + nb213_six]
	mulsd  xmm2, [esp + nb213_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm2, xmm0
	movapd xmm4, xmm2 ;# total fsO 
	movsd [esp + nb213_Vvdwtot], xmm3

	movapd xmm0, [esp + nb213_dxO]
	movapd xmm1, [esp + nb213_dyO]
	movapd xmm2, [esp + nb213_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [esp + nb213_fixO]
	movapd xmm4, [esp + nb213_fiyO]
	movapd xmm7, [esp + nb213_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb213_fixO], xmm3
	movsd [esp + nb213_fiyO], xmm4
	movsd [esp + nb213_fizO], xmm7
	;# update j forces with water O 
	movsd [esp + nb213_fjx], xmm0
	movsd [esp + nb213_fjy], xmm1
	movsd [esp + nb213_fjz], xmm2

	;# H1 interactions
	movsd  xmm6, [esp + nb213_rinvH1] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm6
	movsd  xmm0, [esp + nb213_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [esp + nb213_two]
	subsd   xmm6, [esp + nb213_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [esp + nb213_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb213_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [esp + nb213_vctot]

	movapd xmm0, [esp + nb213_dxH1]
	movapd xmm1, [esp + nb213_dyH1]
	movapd xmm2, [esp + nb213_dzH1]
	movsd [esp + nb213_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb213_fixH1]
	movapd xmm4, [esp + nb213_fiyH1]
	movapd xmm7, [esp + nb213_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb213_fixH1], xmm3
	movsd [esp + nb213_fiyH1], xmm4
	movsd [esp + nb213_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb213_fjx]
	addsd  xmm1, [esp + nb213_fjy]
	addsd  xmm2, [esp + nb213_fjz]
	movsd [esp + nb213_fjx], xmm0
	movsd [esp + nb213_fjy], xmm1
	movsd [esp + nb213_fjz], xmm2

	;# H2 interactions 
	movsd  xmm5, [esp + nb213_rinvH2] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [esp + nb213_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [esp + nb213_two]
	subsd   xmm5, [esp + nb213_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [esp + nb213_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb213_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [esp + nb213_vctot]

	movapd xmm0, [esp + nb213_dxH2]
	movapd xmm1, [esp + nb213_dyH2]
	movapd xmm2, [esp + nb213_dzH2]
	movsd [esp + nb213_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb213_fixH2]
	movapd xmm4, [esp + nb213_fiyH2]
	movapd xmm7, [esp + nb213_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb213_fixH2], xmm3
	movsd [esp + nb213_fiyH2], xmm4
	movsd [esp + nb213_fizH2], xmm7
	;# update j forces with water H2 
	addsd  xmm0, [esp + nb213_fjx]
	addsd  xmm1, [esp + nb213_fjy]
	addsd  xmm2, [esp + nb213_fjz]
	movsd [esp + nb213_fjx], xmm0
	movsd [esp + nb213_fjy], xmm1
	movsd [esp + nb213_fjz], xmm2

	;# M interactions 
	movsd  xmm5, [esp + nb213_rinvM] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [esp + nb213_krsqM]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [esp + nb213_two]
	subsd   xmm5, [esp + nb213_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [esp + nb213_qqM] ;# vcoul 
	mulsd   xmm7, [esp + nb213_qqM]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [esp + nb213_vctot]

	movapd xmm0, [esp + nb213_dxM]
	movapd xmm1, [esp + nb213_dyM]
	movapd xmm2, [esp + nb213_dzM]
	movsd [esp + nb213_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [esp + nb213_fixM]
	movapd xmm4, [esp + nb213_fiyM]
	movapd xmm7, [esp + nb213_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb213_fixM], xmm3
	movsd [esp + nb213_fiyM], xmm4
	movsd [esp + nb213_fizM], xmm7

	mov edi, [ebp + nb213_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb213_fjx]
	addsd  xmm1, [esp + nb213_fjy]
	addsd  xmm2, [esp + nb213_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb213_updateouterdata:
	mov   ecx, [esp + nb213_ii3]
	mov   edi, [ebp + nb213_faction]
	mov   esi, [ebp + nb213_fshift]
	mov   edx, [esp + nb213_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb213_fixO]
	movapd xmm1, [esp + nb213_fiyO]
	movapd xmm2, [esp + nb213_fizO]

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
	movapd xmm0, [esp + nb213_fixH1]
	movapd xmm1, [esp + nb213_fiyH1]
	movapd xmm2, [esp + nb213_fizH1]

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
	movapd xmm0, [esp + nb213_fixH2]
	movapd xmm1, [esp + nb213_fiyH2]
	movapd xmm2, [esp + nb213_fizH2]

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
	movapd xmm0, [esp + nb213_fixM]
	movapd xmm1, [esp + nb213_fiyM]
	movapd xmm2, [esp + nb213_fizM]

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
	mov esi, [esp + nb213_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb213_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb213_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb213_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb213_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb213_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb213_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb213_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb213_n], esi
        jmp .nb213_outer
.nb213_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb213_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb213_end
        ;# non-zero, do one more workunit
        jmp   .nb213_threadloop
.nb213_end:
	emms

	mov eax, [esp + nb213_nouter]
	mov ebx, [esp + nb213_ninner]
	mov ecx, [ebp + nb213_outeriter]
	mov edx, [ebp + nb213_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb213_salign]
	add esp, eax
	add esp, 1004
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel213nf_ia32_sse2
.globl _nb_kernel213nf_ia32_sse2
nb_kernel213nf_ia32_sse2:	
_nb_kernel213nf_ia32_sse2:	
.equiv          nb213nf_p_nri,          8
.equiv          nb213nf_iinr,           12
.equiv          nb213nf_jindex,         16
.equiv          nb213nf_jjnr,           20
.equiv          nb213nf_shift,          24
.equiv          nb213nf_shiftvec,       28
.equiv          nb213nf_fshift,         32
.equiv          nb213nf_gid,            36
.equiv          nb213nf_pos,            40
.equiv          nb213nf_faction,        44
.equiv          nb213nf_charge,         48
.equiv          nb213nf_p_facel,        52
.equiv          nb213nf_argkrf,         56
.equiv          nb213nf_argcrf,         60
.equiv          nb213nf_Vc,             64
.equiv          nb213nf_type,           68
.equiv          nb213nf_p_ntype,        72
.equiv          nb213nf_vdwparam,       76
.equiv          nb213nf_Vvdw,           80
.equiv          nb213nf_p_tabscale,     84
.equiv          nb213nf_VFtab,          88
.equiv          nb213nf_invsqrta,       92
.equiv          nb213nf_dvda,           96
.equiv          nb213nf_p_gbtabscale,   100
.equiv          nb213nf_GBtab,          104
.equiv          nb213nf_p_nthreads,     108
.equiv          nb213nf_count,          112
.equiv          nb213nf_mtx,            116
.equiv          nb213nf_outeriter,      120
.equiv          nb213nf_inneriter,      124
.equiv          nb213nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb213nf_ixO,            0
.equiv          nb213nf_iyO,            16
.equiv          nb213nf_izO,            32
.equiv          nb213nf_ixH1,           48
.equiv          nb213nf_iyH1,           64
.equiv          nb213nf_izH1,           80
.equiv          nb213nf_ixH2,           96
.equiv          nb213nf_iyH2,           112
.equiv          nb213nf_izH2,           128
.equiv          nb213nf_ixM,            144
.equiv          nb213nf_iyM,            160
.equiv          nb213nf_izM,            176
.equiv          nb213nf_iqH,            192
.equiv          nb213nf_iqM,            208
.equiv          nb213nf_qqH,            224
.equiv          nb213nf_qqM,            240
.equiv          nb213nf_c6,             256
.equiv          nb213nf_c12,            272
.equiv          nb213nf_vctot,          288
.equiv          nb213nf_Vvdwtot,        304
.equiv          nb213nf_half,           320
.equiv          nb213nf_three,          336
.equiv          nb213nf_two,            352
.equiv          nb213nf_rinvH1,         368
.equiv          nb213nf_rinvH2,         384
.equiv          nb213nf_rinvM,          400
.equiv          nb213nf_krsqH1,         416
.equiv          nb213nf_krsqH2,         432
.equiv          nb213nf_krsqM,          448
.equiv          nb213nf_krf,            464
.equiv          nb213nf_crf,            480
.equiv          nb213nf_is3,            496
.equiv          nb213nf_ii3,            500
.equiv          nb213nf_ntia,           504
.equiv          nb213nf_innerjjnr,      508
.equiv          nb213nf_innerk,         512
.equiv          nb213nf_n,              516
.equiv          nb213nf_nn1,            520
.equiv          nb213nf_nri,            524
.equiv          nb213nf_nouter,         528
.equiv          nb213nf_ninner,         532
.equiv          nb213nf_salign,         536
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 540		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb213nf_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb213nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb213nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb213nf_nouter], eax
	mov [esp + nb213nf_ninner], eax


	mov esi, [ebp + nb213nf_argkrf]
	mov edi, [ebp + nb213nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb213nf_krf], xmm5
	movapd [esp + nb213nf_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb213nf_half], eax
	mov [esp + nb213nf_half+4], ebx
	movsd xmm1, [esp + nb213nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb213nf_half], xmm1
	movapd [esp + nb213nf_two], xmm2
	movapd [esp + nb213nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb213nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb213nf_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb213nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb213nf_iqH], xmm3
	movapd [esp + nb213nf_iqM], xmm4
	
	mov   edx, [ebp + nb213nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb213nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb213nf_ntia], ecx		
.nb213nf_threadloop:
        mov   esi, [ebp + nb213nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb213nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb213nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb213nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb213nf_n], eax
        mov [esp + nb213nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb213nf_outerstart
        jmp .nb213nf_end

.nb213nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb213nf_nouter]
	mov [esp + nb213nf_nouter], ebx

.nb213nf_outer:
	mov   eax, [ebp + nb213nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb213nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb213nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb213nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb213nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb213nf_ii3], ebx

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
	movapd [esp + nb213nf_ixO], xmm3
	movapd [esp + nb213nf_iyO], xmm4
	movapd [esp + nb213nf_izO], xmm5
	movapd [esp + nb213nf_ixH1], xmm6
	movapd [esp + nb213nf_iyH1], xmm7

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
	movapd [esp + nb213nf_izH1], xmm6
	movapd [esp + nb213nf_ixH2], xmm0
	movapd [esp + nb213nf_iyH2], xmm1
	movapd [esp + nb213nf_izH2], xmm2
	movapd [esp + nb213nf_ixM], xmm3
	movapd [esp + nb213nf_iyM], xmm4
	movapd [esp + nb213nf_izM], xmm5

	;# clear vctot
	xorpd xmm4, xmm4
	movapd [esp + nb213nf_vctot], xmm4
	movapd [esp + nb213nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb213nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb213nf_pos]
	mov   edi, [ebp + nb213nf_faction]	
	mov   eax, [ebp + nb213nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb213nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb213nf_ninner]
	mov   [esp + nb213nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb213nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb213nf_unroll_loop
	jmp   .nb213nf_checksingle
.nb213nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb213nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb213nf_innerjjnr],  8	
	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb213nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb213nf_iqM]
	mulpd  xmm4, [esp + nb213nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb213nf_qqM], xmm3
	movapd  [esp + nb213nf_qqH], xmm4
	
	mov esi, [ebp + nb213nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb213nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb213nf_ntia]
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
	movapd [esp + nb213nf_c6], xmm4
	movapd [esp + nb213nf_c12], xmm6
	
	mov esi, [ebp + nb213nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb213nf_ixO]
	movapd xmm5, [esp + nb213nf_iyO]
	movapd xmm6, [esp + nb213nf_izO]

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
	movapd xmm4, [esp + nb213nf_ixH1]
	movapd xmm5, [esp + nb213nf_iyH1]
	movapd xmm6, [esp + nb213nf_izH1]

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
	movapd xmm3, [esp + nb213nf_ixH2]
	movapd xmm4, [esp + nb213nf_iyH2]
	movapd xmm5, [esp + nb213nf_izH2]

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
	movapd xmm3, [esp + nb213nf_iyM]
	movapd xmm4, [esp + nb213nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb213nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movapd xmm0, [esp + nb213nf_krf]
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	mulpd xmm0, xmm4  
	mulpd xmm1, xmm5
	mulpd xmm2, xmm6
	movapd [esp + nb213nf_krsqM], xmm0
	movapd [esp + nb213nf_krsqH2], xmm1
	movapd [esp + nb213nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb213nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213nf_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb213nf_half] ;# rinv 
	movapd  [esp + nb213nf_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb213nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213nf_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb213nf_half] ;# rinv 
	movapd  [esp + nb213nf_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [esp + nb213nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [esp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213nf_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [esp + nb213nf_half] ;# rinv 
	movapd  [esp + nb213nf_rinvM], xmm1	

	;# do O interactions directly - rsqO is in xmm7
	cvtpd2ps xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtps2pd xmm2, xmm2
	movapd   xmm1, [esp + nb213nf_two]
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
	mulpd  xmm1, [esp + nb213nf_c6]
	mulpd  xmm2, [esp + nb213nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb213nf_Vvdwtot]
	movapd [esp + nb213nf_Vvdwtot], xmm3

	;# H1 interactions 
	movapd  xmm6, [esp + nb213nf_rinvH1]
	addpd   xmm6, [esp + nb213nf_krsqH1]
 	subpd   xmm6, [esp + nb213nf_crf]
	mulpd   xmm6, [esp + nb213nf_qqH] ;# vcoul 
	addpd   xmm6, [esp + nb213nf_vctot]
	movapd [esp + nb213nf_vctot], xmm6
	
	;# H2 interactions 
	movapd  xmm6, [esp + nb213nf_rinvH2]
	addpd   xmm6, [esp + nb213nf_krsqH2]
 	subpd   xmm6, [esp + nb213nf_crf]
	mulpd   xmm6, [esp + nb213nf_qqH] ;# vcoul 
	addpd   xmm6, [esp + nb213nf_vctot]
	movapd [esp + nb213nf_vctot], xmm6

	;# M interactions 
	movapd  xmm6, [esp + nb213nf_rinvM]
	addpd   xmm6, [esp + nb213nf_krsqM]
 	subpd   xmm6, [esp + nb213nf_crf]
	mulpd   xmm6, [esp + nb213nf_qqM] ;# vcoul 
	addpd   xmm6, [esp + nb213nf_vctot]
	movapd [esp + nb213nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb213nf_innerk],  2
	jl   .nb213nf_checksingle
	jmp  .nb213nf_unroll_loop
.nb213nf_checksingle:	
	mov   edx, [esp + nb213nf_innerk]
	and   edx, 1
	jnz  .nb213nf_dosingle
	jmp  .nb213nf_updateouterdata
.nb213nf_dosingle:
	mov   edx, [esp + nb213nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb213nf_innerjjnr],  4	

	mov esi, [ebp + nb213nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [esp + nb213nf_iqM]
	mulsd  xmm4, [esp + nb213nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb213nf_qqM], xmm3
	movapd  [esp + nb213nf_qqH], xmm4
	
	mov esi, [ebp + nb213nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb213nf_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb213nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb213nf_c6], xmm4
	movapd [esp + nb213nf_c12], xmm6
	
	mov esi, [ebp + nb213nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb213nf_ixO]
	movapd xmm5, [esp + nb213nf_iyO]
	movapd xmm6, [esp + nb213nf_izO]

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
	movapd xmm4, [esp + nb213nf_ixH1]
	movapd xmm5, [esp + nb213nf_iyH1]
	movapd xmm6, [esp + nb213nf_izH1]

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
	movapd xmm3, [esp + nb213nf_ixH2]
	movapd xmm4, [esp + nb213nf_iyH2]
	movapd xmm5, [esp + nb213nf_izH2]

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
	movapd xmm3, [esp + nb213nf_iyM]
	movapd xmm4, [esp + nb213nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [esp + nb213nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [esp + nb213nf_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [esp + nb213nf_krsqM], xmm0
	movsd [esp + nb213nf_krsqH2], xmm1
	movsd [esp + nb213nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb213nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213nf_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb213nf_half] ;# rinv 
	movapd [esp + nb213nf_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb213nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213nf_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb213nf_half] ;# rinv 
	movapd [esp + nb213nf_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [esp + nb213nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [esp + nb213nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [esp + nb213nf_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [esp + nb213nf_half] ;# rinv 
	movapd [esp + nb213nf_rinvM], xmm1

	;# do O interactions directly. xmm7=rsq
	cvtsd2ss xmm2, xmm7
	movapd   xmm6, xmm7
	rcpps    xmm2, xmm2
	cvtss2sd xmm2, xmm2
	movapd   xmm1, [esp + nb213nf_two]
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
	mulsd  xmm1, [esp + nb213nf_c6]
	mulsd  xmm2, [esp + nb213nf_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb213nf_Vvdwtot]
	movsd [esp + nb213nf_Vvdwtot], xmm3

	;# H1 interactions 
	movsd  xmm6, [esp + nb213nf_rinvH1]
	addsd   xmm6, [esp + nb213nf_krsqH1]
 	subsd   xmm6, [esp + nb213nf_crf]
	mulsd   xmm6, [esp + nb213nf_qqH] ;# vcoul 
	addsd   xmm6, [esp + nb213nf_vctot]
	movsd [esp + nb213nf_vctot], xmm6
	
	;# H2 interactions 
	movsd  xmm6, [esp + nb213nf_rinvH2]
	addsd   xmm6, [esp + nb213nf_krsqH2]
 	subsd   xmm6, [esp + nb213nf_crf]
	mulsd   xmm6, [esp + nb213nf_qqH] ;# vcoul 
	addsd   xmm6, [esp + nb213nf_vctot]
	movsd [esp + nb213nf_vctot], xmm6

	;# M interactions 
	movsd  xmm6, [esp + nb213nf_rinvM]
	addsd   xmm6, [esp + nb213nf_krsqM]
 	subsd   xmm6, [esp + nb213nf_crf]
	mulsd   xmm6, [esp + nb213nf_qqM] ;# vcoul 
	addsd   xmm6, [esp + nb213nf_vctot]
	movsd [esp + nb213nf_vctot], xmm6
	
.nb213nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb213nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb213nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb213nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb213nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb213nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb213nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb213nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb213nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb213nf_n], esi
        jmp .nb213nf_outer
.nb213nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb213nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb213nf_end
        ;# non-zero, do one more workunit
        jmp   .nb213nf_threadloop
.nb213nf_end:
	emms

	mov eax, [esp + nb213nf_nouter]
	mov ebx, [esp + nb213nf_ninner]
	mov ecx, [ebp + nb213nf_outeriter]
	mov edx, [ebp + nb213nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb213nf_salign]
	add esp, eax
	add esp, 540
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

