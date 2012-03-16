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


	

.globl nb_kernel203_ia32_sse2
.globl _nb_kernel203_ia32_sse2
nb_kernel203_ia32_sse2:	
_nb_kernel203_ia32_sse2:	
.equiv          nb203_p_nri,            8
.equiv          nb203_iinr,             12
.equiv          nb203_jindex,           16
.equiv          nb203_jjnr,             20
.equiv          nb203_shift,            24
.equiv          nb203_shiftvec,         28
.equiv          nb203_fshift,           32
.equiv          nb203_gid,              36
.equiv          nb203_pos,              40
.equiv          nb203_faction,          44
.equiv          nb203_charge,           48
.equiv          nb203_p_facel,          52
.equiv          nb203_argkrf,           56
.equiv          nb203_argcrf,           60
.equiv          nb203_Vc,               64
.equiv          nb203_type,             68
.equiv          nb203_p_ntype,          72
.equiv          nb203_vdwparam,         76
.equiv          nb203_Vvdw,             80
.equiv          nb203_p_tabscale,       84
.equiv          nb203_VFtab,            88
.equiv          nb203_invsqrta,         92
.equiv          nb203_dvda,             96
.equiv          nb203_p_gbtabscale,     100
.equiv          nb203_GBtab,            104
.equiv          nb203_p_nthreads,       108
.equiv          nb203_count,            112
.equiv          nb203_mtx,              116
.equiv          nb203_outeriter,        120
.equiv          nb203_inneriter,        124
.equiv          nb203_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb203_ixM,              0
.equiv          nb203_iyM,              16
.equiv          nb203_izM,              32
.equiv          nb203_ixH1,             48
.equiv          nb203_iyH1,             64
.equiv          nb203_izH1,             80
.equiv          nb203_ixH2,             96
.equiv          nb203_iyH2,             112
.equiv          nb203_izH2,             128
.equiv          nb203_iqM,              144
.equiv          nb203_iqH,              160
.equiv          nb203_dxM,              176
.equiv          nb203_dyM,              192
.equiv          nb203_dzM,              208
.equiv          nb203_dxH1,             224
.equiv          nb203_dyH1,             240
.equiv          nb203_dzH1,             256
.equiv          nb203_dxH2,             272
.equiv          nb203_dyH2,             288
.equiv          nb203_dzH2,             304
.equiv          nb203_qqM,              320
.equiv          nb203_qqH,              336
.equiv          nb203_vctot,            352
.equiv          nb203_fixM,             384
.equiv          nb203_fiyM,             400
.equiv          nb203_fizM,             416
.equiv          nb203_fixH1,            432
.equiv          nb203_fiyH1,            448
.equiv          nb203_fizH1,            464
.equiv          nb203_fixH2,            480
.equiv          nb203_fiyH2,            496
.equiv          nb203_fizH2,            512
.equiv          nb203_fjx,              528
.equiv          nb203_fjy,              544
.equiv          nb203_fjz,              560
.equiv          nb203_half,             576
.equiv          nb203_three,            592
.equiv          nb203_two,              608
.equiv          nb203_krf,              624
.equiv          nb203_crf,              640
.equiv          nb203_krsqM,            656
.equiv          nb203_krsqH1,           672
.equiv          nb203_krsqH2,           688
.equiv          nb203_is3,              704
.equiv          nb203_ii3,              708
.equiv          nb203_innerjjnr,        712
.equiv          nb203_innerk,           716
.equiv          nb203_n,                720
.equiv          nb203_nn1,              724
.equiv          nb203_nri,              728
.equiv          nb203_nouter,           732
.equiv          nb203_ninner,           736
.equiv          nb203_salign,           740
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
	mov [esp + nb203_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb203_p_nri]
	mov ecx, [ecx]
	mov [esp + nb203_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb203_nouter], eax
	mov [esp + nb203_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb203_half], eax
	mov [esp + nb203_half+4], ebx
	movsd xmm1, [esp + nb203_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb203_half], xmm1
	movapd [esp + nb203_two], xmm2
	movapd [esp + nb203_three], xmm3

	mov esi, [ebp + nb203_argkrf]
	mov edi, [ebp + nb203_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb203_krf], xmm5
	movapd [esp + nb203_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb203_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb203_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb203_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb203_iqH], xmm3
	movapd [esp + nb203_iqM], xmm4
			
.nb203_threadloop:
        mov   esi, [ebp + nb203_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb203_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb203_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb203_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb203_n], eax
        mov [esp + nb203_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb203_outerstart
        jmp .nb203_end

.nb203_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb203_nouter]
	mov [esp + nb203_nouter], ebx

.nb203_outer:
	mov   eax, [ebp + nb203_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb203_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb203_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb203_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb203_pos]    ;# eax = base of pos[]  
	mov   [esp + nb203_ii3], ebx

	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb203_ixH1], xmm3
	movapd [esp + nb203_iyH1], xmm4
	movapd [esp + nb203_izH1], xmm5

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
	movapd [esp + nb203_ixH2], xmm0
	movapd [esp + nb203_iyH2], xmm1
	movapd [esp + nb203_izH2], xmm2
	movapd [esp + nb203_ixM], xmm3
	movapd [esp + nb203_iyM], xmm4
	movapd [esp + nb203_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb203_vctot], xmm4
	movapd [esp + nb203_fixM], xmm4
	movapd [esp + nb203_fiyM], xmm4
	movapd [esp + nb203_fizM], xmm4
	movapd [esp + nb203_fixH1], xmm4
	movapd [esp + nb203_fiyH1], xmm4
	movapd [esp + nb203_fizH1], xmm4
	movapd [esp + nb203_fixH2], xmm4
	movapd [esp + nb203_fiyH2], xmm4
	movapd [esp + nb203_fizH2], xmm4
	
	mov   eax, [ebp + nb203_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb203_pos]
	mov   edi, [ebp + nb203_faction]	
	mov   eax, [ebp + nb203_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb203_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb203_ninner]
	mov   [esp + nb203_ninner], ecx
	add   edx, 0
	mov   [esp + nb203_innerk], edx    ;# number of innerloop atoms 
	jge   .nb203_unroll_loop
	jmp   .nb203_checksingle
.nb203_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb203_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb203_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb203_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb203_iqM]
	mulpd  xmm4, [esp + nb203_iqH]
	movapd  [esp + nb203_qqM], xmm3
	movapd  [esp + nb203_qqH], xmm4	

	mov esi, [ebp + nb203_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb203_ixM]
	movapd xmm5, [esp + nb203_iyM]
	movapd xmm6, [esp + nb203_izM]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb203_dxM], xmm4
	movapd [esp + nb203_dyM], xmm5
	movapd [esp + nb203_dzM], xmm6

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb203_ixH1]
	movapd xmm5, [esp + nb203_iyH1]
	movapd xmm6, [esp + nb203_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb203_dxH1], xmm4
	movapd [esp + nb203_dyH1], xmm5
	movapd [esp + nb203_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb203_ixH2]
	movapd xmm4, [esp + nb203_iyH2]
	movapd xmm5, [esp + nb203_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb203_dxH2], xmm3
	movapd [esp + nb203_dyH2], xmm4
	movapd [esp + nb203_dzH2], xmm5
	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulpd  xmm0, [esp + nb203_krf]	
	mulpd  xmm1, [esp + nb203_krf]	
	mulpd  xmm2, [esp + nb203_krf]	

	movapd [esp + nb203_krsqH2], xmm0
	movapd [esp + nb203_krsqH1], xmm1
	movapd [esp + nb203_krsqM], xmm2
	
	;# start with rsqM - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb203_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb203_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb203_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb203_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb203_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb203_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb203_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb203_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb203_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	movapd  xmm4, xmm7	
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm3, xmm7
	movapd  xmm0, [esp + nb203_krsqM]
	addpd   xmm7, xmm0	;# xmm6=rinv+ krsq 
	mulpd   xmm0, [esp + nb203_two]
	subpd   xmm7, [esp + nb203_crf]
	subpd   xmm3, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm7, [esp + nb203_qqM] ;# vcoul 
	mulpd   xmm3, [esp + nb203_qqM]
	mulpd  xmm4, xmm3	;# total fsH1 in xmm4 
	
	addpd  xmm7, [esp + nb203_vctot]

	movapd xmm0, [esp + nb203_dxM]
	movapd xmm1, [esp + nb203_dyM]
	movapd xmm2, [esp + nb203_dzM]
	movapd [esp + nb203_vctot], xmm7
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	
	;# update M forces 
	movapd xmm3, [esp + nb203_fixM]
	movapd xmm4, [esp + nb203_fiyM]
	movapd xmm7, [esp + nb203_fizM]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb203_fixM], xmm3
	movapd [esp + nb203_fiyM], xmm4
	movapd [esp + nb203_fizM], xmm7
	;# update j forces with water M 
	movapd [esp + nb203_fjx], xmm0
	movapd [esp + nb203_fjy], xmm1
	movapd [esp + nb203_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [esp + nb203_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulpd   xmm0, [esp + nb203_two]
	subpd   xmm6, [esp + nb203_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm6, [esp + nb203_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb203_qqH]
	mulpd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addpd  xmm6, [esp + nb203_vctot]

	movapd xmm0, [esp + nb203_dxH1]
	movapd xmm1, [esp + nb203_dyH1]
	movapd xmm2, [esp + nb203_dzH1]
	movapd [esp + nb203_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb203_fixH1]
	movapd xmm4, [esp + nb203_fiyH1]
	movapd xmm7, [esp + nb203_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb203_fixH1], xmm3
	movapd [esp + nb203_fiyH1], xmm4
	movapd [esp + nb203_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb203_fjx]
	addpd  xmm1, [esp + nb203_fjy]
	addpd  xmm2, [esp + nb203_fjz]
	movapd [esp + nb203_fjx], xmm0
	movapd [esp + nb203_fjy], xmm1
	movapd [esp + nb203_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb203_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulpd   xmm0, [esp + nb203_two]
	subpd   xmm5, [esp + nb203_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm5, [esp + nb203_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb203_qqH]
	mulpd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addpd  xmm5, [esp + nb203_vctot]

	movapd xmm0, [esp + nb203_dxH2]
	movapd xmm1, [esp + nb203_dyH2]
	movapd xmm2, [esp + nb203_dzH2]
	movapd [esp + nb203_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb203_fixH2]
	movapd xmm4, [esp + nb203_fiyH2]
	movapd xmm7, [esp + nb203_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb203_fixH2], xmm3
	movapd [esp + nb203_fiyH2], xmm4
	movapd [esp + nb203_fizH2], xmm7

	mov edi, [ebp + nb203_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb203_fjx]
	addpd  xmm1, [esp + nb203_fjy]
	addpd  xmm2, [esp + nb203_fjz]
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
	sub dword ptr [esp + nb203_innerk],  2
	jl    .nb203_checksingle
	jmp   .nb203_unroll_loop
.nb203_checksingle:	
	mov   edx, [esp + nb203_innerk]
	and   edx, 1
	jnz   .nb203_dosingle
	jmp   .nb203_updateouterdata
.nb203_dosingle:
	mov   edx, [esp + nb203_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb203_innerjjnr],  4	

	mov esi, [ebp + nb203_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb203_iqM]
	mulpd  xmm4, [esp + nb203_iqH]
	movapd  [esp + nb203_qqM], xmm3
	movapd  [esp + nb203_qqH], xmm4
	
	mov esi, [ebp + nb203_pos]       ;# base of pos[] 
	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb203_ixM]
	movapd xmm5, [esp + nb203_iyM]
	movapd xmm6, [esp + nb203_izM]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb203_dxM], xmm4
	movapd [esp + nb203_dyM], xmm5
	movapd [esp + nb203_dzM], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb203_ixH1]
	movapd xmm5, [esp + nb203_iyH1]
	movapd xmm6, [esp + nb203_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb203_dxH1], xmm4
	movapd [esp + nb203_dyH1], xmm5
	movapd [esp + nb203_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb203_ixH2]
	movapd xmm4, [esp + nb203_iyH2]
	movapd xmm5, [esp + nb203_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb203_dxH2], xmm3
	movapd [esp + nb203_dyH2], xmm4
	movapd [esp + nb203_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 
	
	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulsd  xmm0, [esp + nb203_krf]	
	mulsd  xmm1, [esp + nb203_krf]	
	mulsd  xmm2, [esp + nb203_krf]	

	movapd [esp + nb203_krsqH2], xmm0
	movapd [esp + nb203_krsqH1], xmm1
	movapd [esp + nb203_krsqM], xmm2
	
	;# start with rsqM - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb203_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb203_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb203_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb203_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb203_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb203_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb203_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb203_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb203_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm3, xmm7
	movapd  xmm0, [esp + nb203_krsqM]
	addsd   xmm7, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [esp + nb203_two]
	subsd   xmm7, [esp + nb203_crf]
	subsd   xmm3, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm7, [esp + nb203_qqM] ;# vcoul 
	mulsd   xmm3, [esp + nb203_qqM]
	mulsd  xmm4, xmm3	;# total fsH1 in xmm4 
	
	addsd  xmm7, [esp + nb203_vctot]

	movapd xmm0, [esp + nb203_dxM]
	movapd xmm1, [esp + nb203_dyM]
	movapd xmm2, [esp + nb203_dzM]
	movlpd [esp + nb203_vctot], xmm7
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [esp + nb203_fixM]
	movapd xmm4, [esp + nb203_fiyM]
	movapd xmm7, [esp + nb203_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb203_fixM], xmm3
	movlpd [esp + nb203_fiyM], xmm4
	movlpd [esp + nb203_fizM], xmm7
	;# update j forces with water M 
	movlpd [esp + nb203_fjx], xmm0
	movlpd [esp + nb203_fjy], xmm1
	movlpd [esp + nb203_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [esp + nb203_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [esp + nb203_two]
	subsd   xmm6, [esp + nb203_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [esp + nb203_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb203_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [esp + nb203_vctot]

	movapd xmm0, [esp + nb203_dxH1]
	movapd xmm1, [esp + nb203_dyH1]
	movapd xmm2, [esp + nb203_dzH1]
	movlpd [esp + nb203_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb203_fixH1]
	movapd xmm4, [esp + nb203_fiyH1]
	movapd xmm7, [esp + nb203_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb203_fixH1], xmm3
	movlpd [esp + nb203_fiyH1], xmm4
	movlpd [esp + nb203_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb203_fjx]
	addsd  xmm1, [esp + nb203_fjy]
	addsd  xmm2, [esp + nb203_fjz]
	movlpd [esp + nb203_fjx], xmm0
	movlpd [esp + nb203_fjy], xmm1
	movlpd [esp + nb203_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb203_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [esp + nb203_two]
	subsd   xmm5, [esp + nb203_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [esp + nb203_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb203_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [esp + nb203_vctot]

	movapd xmm0, [esp + nb203_dxH2]
	movapd xmm1, [esp + nb203_dyH2]
	movapd xmm2, [esp + nb203_dzH2]
	movlpd [esp + nb203_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb203_fixH2]
	movapd xmm4, [esp + nb203_fiyH2]
	movapd xmm7, [esp + nb203_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb203_fixH2], xmm3
	movlpd [esp + nb203_fiyH2], xmm4
	movlpd [esp + nb203_fizH2], xmm7

	mov edi, [ebp + nb203_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb203_fjx]
	addsd  xmm1, [esp + nb203_fjy]
	addsd  xmm2, [esp + nb203_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5	

.nb203_updateouterdata:
	mov   ecx, [esp + nb203_ii3]
	mov   edi, [ebp + nb203_faction]
	mov   esi, [ebp + nb203_fshift]
	mov   edx, [esp + nb203_is3]

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb203_fixH1]
	movapd xmm1, [esp + nb203_fiyH1]
	movapd xmm2, [esp + nb203_fizH1]

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
	movapd xmm0, [esp + nb203_fixH2]
	movapd xmm1, [esp + nb203_fiyH2]
	movapd xmm2, [esp + nb203_fizH2]

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
	movapd xmm0, [esp + nb203_fixM]
	movapd xmm1, [esp + nb203_fiyM]
	movapd xmm2, [esp + nb203_fizM]

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
	mov esi, [esp + nb203_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb203_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb203_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb203_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb203_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb203_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb203_n], esi
        jmp .nb203_outer
.nb203_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb203_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb203_end
        ;# non-zero, do one more workunit
        jmp   .nb203_threadloop
.nb203_end:
	emms

	mov eax, [esp + nb203_nouter]
	mov ebx, [esp + nb203_ninner]
	mov ecx, [ebp + nb203_outeriter]
	mov edx, [ebp + nb203_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb203_salign]
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




.globl nb_kernel203nf_ia32_sse2
.globl _nb_kernel203nf_ia32_sse2
nb_kernel203nf_ia32_sse2:	
_nb_kernel203nf_ia32_sse2:	
.equiv          nb203nf_p_nri,          8
.equiv          nb203nf_iinr,           12
.equiv          nb203nf_jindex,         16
.equiv          nb203nf_jjnr,           20
.equiv          nb203nf_shift,          24
.equiv          nb203nf_shiftvec,       28
.equiv          nb203nf_fshift,         32
.equiv          nb203nf_gid,            36
.equiv          nb203nf_pos,            40
.equiv          nb203nf_faction,        44
.equiv          nb203nf_charge,         48
.equiv          nb203nf_p_facel,        52
.equiv          nb203nf_argkrf,         56
.equiv          nb203nf_argcrf,         60
.equiv          nb203nf_Vc,             64
.equiv          nb203nf_type,           68
.equiv          nb203nf_p_ntype,        72
.equiv          nb203nf_vdwparam,       76
.equiv          nb203nf_Vvdw,           80
.equiv          nb203nf_p_tabscale,     84
.equiv          nb203nf_VFtab,          88
.equiv          nb203nf_invsqrta,       92
.equiv          nb203nf_dvda,           96
.equiv          nb203nf_p_gbtabscale,   100
.equiv          nb203nf_GBtab,          104
.equiv          nb203nf_p_nthreads,     108
.equiv          nb203nf_count,          112
.equiv          nb203nf_mtx,            116
.equiv          nb203nf_outeriter,      120
.equiv          nb203nf_inneriter,      124
.equiv          nb203nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb203nf_ixM,            0
.equiv          nb203nf_iyM,            16
.equiv          nb203nf_izM,            32
.equiv          nb203nf_ixH1,           48
.equiv          nb203nf_iyH1,           64
.equiv          nb203nf_izH1,           80
.equiv          nb203nf_ixH2,           96
.equiv          nb203nf_iyH2,           112
.equiv          nb203nf_izH2,           128
.equiv          nb203nf_iqM,            144
.equiv          nb203nf_iqH,            160
.equiv          nb203nf_qqM,            176
.equiv          nb203nf_qqH,            192
.equiv          nb203nf_vctot,          208
.equiv          nb203nf_half,           224
.equiv          nb203nf_three,          240
.equiv          nb203nf_krf,            256
.equiv          nb203nf_crf,            272
.equiv          nb203nf_krsqM,          288
.equiv          nb203nf_krsqH1,         304
.equiv          nb203nf_krsqH2,         320
.equiv          nb203nf_is3,            336
.equiv          nb203nf_ii3,            340
.equiv          nb203nf_innerjjnr,      344
.equiv          nb203nf_innerk,         348
.equiv          nb203nf_n,              352
.equiv          nb203nf_nn1,            356
.equiv          nb203nf_nri,            360
.equiv          nb203nf_nouter,         364
.equiv          nb203nf_ninner,         368
.equiv          nb203nf_salign,         372
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 376		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb203nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb203nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb203nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb203nf_nouter], eax
	mov [esp + nb203nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb203nf_half], eax
	mov [esp + nb203nf_half+4], ebx
	movsd xmm1, [esp + nb203nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb203nf_half], xmm1
	movapd [esp + nb203nf_three], xmm3

	mov esi, [ebp + nb203nf_argkrf]
	mov edi, [ebp + nb203nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb203nf_krf], xmm5
	movapd [esp + nb203nf_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb203nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb203nf_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb203nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb203nf_iqH], xmm3
	movapd [esp + nb203nf_iqM], xmm4
			
.nb203nf_threadloop:
        mov   esi, [ebp + nb203nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb203nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb203nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb203nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb203nf_n], eax
        mov [esp + nb203nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb203nf_outerstart
        jmp .nb203nf_end

.nb203nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb203nf_nouter]
	mov [esp + nb203nf_nouter], ebx

.nb203nf_outer:
	mov   eax, [ebp + nb203nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb203nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb203nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb203nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb203nf_ii3], ebx

	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb203nf_ixH1], xmm3
	movapd [esp + nb203nf_iyH1], xmm4
	movapd [esp + nb203nf_izH1], xmm5

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
	movapd [esp + nb203nf_ixH2], xmm0
	movapd [esp + nb203nf_iyH2], xmm1
	movapd [esp + nb203nf_izH2], xmm2
	movapd [esp + nb203nf_ixM], xmm3
	movapd [esp + nb203nf_iyM], xmm4
	movapd [esp + nb203nf_izM], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb203nf_vctot], xmm4
	
	mov   eax, [ebp + nb203nf_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb203nf_pos]
	mov   eax, [ebp + nb203nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb203nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb203nf_ninner]
	mov   [esp + nb203nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb203nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb203nf_unroll_loop
	jmp   .nb203nf_checksingle
.nb203nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb203nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb203nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb203nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb203nf_iqM]
	mulpd  xmm4, [esp + nb203nf_iqH]
	movapd  [esp + nb203nf_qqM], xmm3
	movapd  [esp + nb203nf_qqH], xmm4	

	mov esi, [ebp + nb203nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb203nf_ixM]
	movapd xmm5, [esp + nb203nf_iyM]
	movapd xmm6, [esp + nb203nf_izM]

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
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb203nf_ixH1]
	movapd xmm5, [esp + nb203nf_iyH1]
	movapd xmm6, [esp + nb203nf_izH1]

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
	movapd xmm3, [esp + nb203nf_ixH2]
	movapd xmm4, [esp + nb203nf_iyH2]
	movapd xmm5, [esp + nb203nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulpd  xmm0, [esp + nb203nf_krf]	
	mulpd  xmm1, [esp + nb203nf_krf]	
	mulpd  xmm2, [esp + nb203nf_krf]	

	movapd [esp + nb203nf_krsqH2], xmm0
	movapd [esp + nb203nf_krsqH1], xmm1
	movapd [esp + nb203nf_krsqM], xmm2
	
	;# start with rsqM - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb203nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb203nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb203nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb203nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb203nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb203nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb203nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb203nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb203nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	movapd  xmm0, [esp + nb203nf_krsqM]
	addpd   xmm7, xmm0	;# xmm7=rinv+ krsq 
	subpd   xmm7, [esp + nb203nf_crf]
	mulpd   xmm7, [esp + nb203nf_qqM] ;# vcoul 	
	addpd  xmm7, [esp + nb203nf_vctot]

	;# H1 interactions 
	movapd  xmm0, [esp + nb203nf_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subpd   xmm6, [esp + nb203nf_crf]
	mulpd   xmm6, [esp + nb203nf_qqH] ;# vcoul 
	addpd  xmm6, xmm7

	;# H2 interactions 
	movapd  xmm0, [esp + nb203nf_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [esp + nb203nf_crf]
	mulpd   xmm5, [esp + nb203nf_qqH] ;# vcoul 
	addpd  xmm5, xmm6
	movapd [esp + nb203nf_vctot], xmm5
		
	;# should we do one more iteration? 
	sub dword ptr [esp + nb203nf_innerk],  2
	jl    .nb203nf_checksingle
	jmp   .nb203nf_unroll_loop
.nb203nf_checksingle:	
	mov   edx, [esp + nb203nf_innerk]
	and   edx, 1
	jnz   .nb203nf_dosingle
	jmp   .nb203nf_updateouterdata
.nb203nf_dosingle:
	mov   edx, [esp + nb203nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb203nf_innerjjnr],  4	

	mov esi, [ebp + nb203nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb203nf_iqM]
	mulpd  xmm4, [esp + nb203nf_iqH]
	movapd  [esp + nb203nf_qqM], xmm3
	movapd  [esp + nb203nf_qqH], xmm4
	
	mov esi, [ebp + nb203nf_pos]       ;# base of pos[] 
	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb203nf_ixM]
	movapd xmm5, [esp + nb203nf_iyM]
	movapd xmm6, [esp + nb203nf_izM]

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
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb203nf_ixH1]
	movapd xmm5, [esp + nb203nf_iyH1]
	movapd xmm6, [esp + nb203nf_izH1]

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
	movapd xmm3, [esp + nb203nf_ixH2]
	movapd xmm4, [esp + nb203nf_iyH2]
	movapd xmm5, [esp + nb203nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 
	
	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulsd  xmm0, [esp + nb203nf_krf]	
	mulsd  xmm1, [esp + nb203nf_krf]	
	mulsd  xmm2, [esp + nb203nf_krf]	

	movapd [esp + nb203nf_krsqH2], xmm0
	movapd [esp + nb203nf_krsqH1], xmm1
	movapd [esp + nb203nf_krsqM], xmm2
	
	;# start with rsqM - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb203nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb203nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb203nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb203nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb203nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb203nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb203nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb203nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb203nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb203nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	movapd  xmm0, [esp + nb203nf_krsqM]
	addsd   xmm7, xmm0	;# xmm7=rinv+ krsq 
	subsd   xmm7, [esp + nb203nf_crf]
	mulsd   xmm7, [esp + nb203nf_qqM] ;# vcoul 	
	addsd  xmm7, [esp + nb203nf_vctot]

	;# H1 interactions 
	movapd  xmm0, [esp + nb203nf_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subsd   xmm6, [esp + nb203nf_crf]
	mulsd   xmm6, [esp + nb203nf_qqH] ;# vcoul 
	addsd  xmm6, xmm7

	;# H2 interactions 
	movapd  xmm0, [esp + nb203nf_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [esp + nb203nf_crf]
	mulsd   xmm5, [esp + nb203nf_qqH] ;# vcoul 
	addsd  xmm5, xmm6
	movlpd [esp + nb203nf_vctot], xmm5
		
.nb203nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb203nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb203nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb203nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb203nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb203nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb203nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb203nf_n], esi
        jmp .nb203nf_outer
.nb203nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb203nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb203nf_end
        ;# non-zero, do one more workunit
        jmp   .nb203nf_threadloop
.nb203nf_end:
	emms

	mov eax, [esp + nb203nf_nouter]
	mov ebx, [esp + nb203nf_ninner]
	mov ecx, [ebp + nb203nf_outeriter]
	mov edx, [ebp + nb203nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb203nf_salign]
	add esp, eax
	add esp, 376
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
