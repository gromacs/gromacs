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


.globl nb_kernel230_ia32_sse2
.globl _nb_kernel230_ia32_sse2
nb_kernel230_ia32_sse2:	
_nb_kernel230_ia32_sse2:	
.equiv          nb230_p_nri,            8
.equiv          nb230_iinr,             12
.equiv          nb230_jindex,           16
.equiv          nb230_jjnr,             20
.equiv          nb230_shift,            24
.equiv          nb230_shiftvec,         28
.equiv          nb230_fshift,           32
.equiv          nb230_gid,              36
.equiv          nb230_pos,              40
.equiv          nb230_faction,          44
.equiv          nb230_charge,           48
.equiv          nb230_p_facel,          52
.equiv          nb230_argkrf,           56
.equiv          nb230_argcrf,           60
.equiv          nb230_Vc,               64
.equiv          nb230_type,             68
.equiv          nb230_p_ntype,          72
.equiv          nb230_vdwparam,         76
.equiv          nb230_Vvdw,             80
.equiv          nb230_p_tabscale,       84
.equiv          nb230_VFtab,            88
.equiv          nb230_invsqrta,         92
.equiv          nb230_dvda,             96
.equiv          nb230_p_gbtabscale,     100
.equiv          nb230_GBtab,            104
.equiv          nb230_p_nthreads,       108
.equiv          nb230_count,            112
.equiv          nb230_mtx,              116
.equiv          nb230_outeriter,        120
.equiv          nb230_inneriter,        124
.equiv          nb230_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb230_ix,               0
.equiv          nb230_iy,               16
.equiv          nb230_iz,               32
.equiv          nb230_iq,               48
.equiv          nb230_dx,               64
.equiv          nb230_dy,               80
.equiv          nb230_dz,               96
.equiv          nb230_c6,               112
.equiv          nb230_c12,              128
.equiv          nb230_tsc,              144
.equiv          nb230_fstmp,            160
.equiv          nb230_vctot,            176
.equiv          nb230_Vvdwtot,          192
.equiv          nb230_fix,              208
.equiv          nb230_fiy,              224
.equiv          nb230_fiz,              240
.equiv          nb230_half,             256
.equiv          nb230_three,            272
.equiv          nb230_two,              288
.equiv          nb230_krf,              304
.equiv          nb230_crf,              320
.equiv          nb230_is3,              336
.equiv          nb230_ii3,              340
.equiv          nb230_ntia,             344
.equiv          nb230_innerjjnr,        348
.equiv          nb230_innerk,           352
.equiv          nb230_n,                356
.equiv          nb230_nn1,              360
.equiv          nb230_nri,              364
.equiv          nb230_facel,            368   ;# uses 8 bytes
.equiv          nb230_ntype,            376
.equiv          nb230_nouter,           380
.equiv          nb230_ninner,           384
.equiv          nb230_salign,           388
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  392		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb230_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb230_p_nri]
	mov esi, [ebp + nb230_p_facel]
	mov edi, [ebp + nb230_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb230_nri], ecx
	movsd [esp + nb230_facel], xmm7
	mov [esp + nb230_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb230_nouter], eax
	mov [esp + nb230_ninner], eax

	mov eax, [ebp + nb230_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb230_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb230_half], eax
	mov [esp + nb230_half+4], ebx
	movsd xmm1, [esp + nb230_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb230_half], xmm1
	movapd [esp + nb230_two], xmm2
	movapd [esp + nb230_three], xmm3

	mov esi, [ebp + nb230_argkrf]
	mov edi, [ebp + nb230_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb230_krf], xmm5
	movapd [esp + nb230_crf], xmm6

.nb230_threadloop:
        mov   esi, [ebp + nb230_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb230_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb230_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb230_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb230_n], eax
        mov [esp + nb230_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb230_outerstart
        jmp .nb230_end

.nb230_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb230_nouter]
	mov [esp + nb230_nouter], ebx

.nb230_outer:
	mov   eax, [ebp + nb230_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb230_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb230_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb230_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb230_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb230_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb230_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb230_ntype]
    	shl   edx, 1
    	mov   [esp + nb230_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb230_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb230_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb230_ix], xmm0
	movapd [esp + nb230_iy], xmm1
	movapd [esp + nb230_iz], xmm2

	mov   [esp + nb230_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb230_vctot], xmm4
	movapd [esp + nb230_Vvdwtot], xmm4
	movapd [esp + nb230_fix], xmm4
	movapd [esp + nb230_fiy], xmm4
	movapd [esp + nb230_fiz], xmm4
	
	mov   eax, [ebp + nb230_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb230_pos]
	mov   edi, [ebp + nb230_faction]	
	mov   eax, [ebp + nb230_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb230_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb230_ninner]
	mov   [esp + nb230_ninner], ecx
	add   edx, 0
	mov   [esp + nb230_innerk], edx    ;# number of innerloop atoms 
	jge   .nb230_unroll_loop
	jmp   .nb230_checksingle
.nb230_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb230_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb230_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb230_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb230_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb230_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb230_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb230_ntia]
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
	movapd [esp + nb230_c6], xmm4
	movapd [esp + nb230_c12], xmm6
	
	mov esi, [ebp + nb230_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb230_ix]
	movapd xmm5, [esp + nb230_iy]
	movapd xmm6, [esp + nb230_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb230_dx], xmm4
	movapd [esp + nb230_dy], xmm5
	movapd [esp + nb230_dz], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	movapd xmm7, [esp + nb230_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb230_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb230_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [esp + nb230_crf]
	mulpd  xmm6, xmm3
	mulpd  xmm7, [esp + nb230_two]
	movapd xmm1, xmm0
	subpd  xmm1, xmm7   ;# rinv-2*krsq
	mulpd  xmm1, xmm0   ;# (rinv-2*krsq)*rinv
	mulpd  xmm1, xmm3   ;# qq*(rinv-2*krsq)*rinv
	
	movapd [esp + nb230_fstmp], xmm1
	
	addpd  xmm6, [esp + nb230_vctot]
	movapd [esp + nb230_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb230_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb230_VFtab]
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
	mulpd  xmm7, [esp + nb230_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb230_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb230_Vvdwtot]
	movapd xmm3, [esp + nb230_fstmp]
	mulpd  xmm7, [esp + nb230_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb230_fstmp], xmm3
	movapd [esp + nb230_Vvdwtot], xmm5

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
	mulpd  xmm7, [esp + nb230_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb230_c12]
	mulpd  xmm7, xmm4 
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb230_Vvdwtot]
	movapd xmm3, [esp + nb230_fstmp]
	mulpd  xmm7, [esp + nb230_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb230_Vvdwtot], xmm5

	mulpd  xmm3, xmm0
		
	movapd xmm0, [esp + nb230_dx]
	movapd xmm1, [esp + nb230_dy]
	movapd xmm2, [esp + nb230_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb230_faction]
	mulpd  xmm0, xmm3
	mulpd  xmm1, xmm3
	mulpd  xmm2, xmm3
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb230_fix]
	movapd xmm4, [esp + nb230_fiy]
	movapd xmm5, [esp + nb230_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb230_fix], xmm3
	movapd [esp + nb230_fiy], xmm4
	movapd [esp + nb230_fiz], xmm5
	;# the fj's - start by accumulating forces from memory 
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
	sub dword ptr [esp + nb230_innerk],  2
	jl    .nb230_checksingle
	jmp   .nb230_unroll_loop

.nb230_checksingle:				
	mov   edx, [esp + nb230_innerk]
	and   edx, 1
	jnz    .nb230_dosingle
	jmp    .nb230_updateouterdata
.nb230_dosingle:			
	mov esi, [ebp + nb230_charge]
	mov edi, [ebp + nb230_pos]
	mov   ecx, [esp + nb230_innerjjnr]
	xorpd xmm3, xmm3
	mov   eax, [ecx]

	movlpd xmm3, [esi + eax*8]
	movapd xmm5, [esp + nb230_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb230_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb230_vdwparam]
	shl eax, 1
	mov edi, [esp + nb230_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		
	movapd [esp + nb230_c6], xmm4
	movapd [esp + nb230_c12], xmm6
	
	mov esi, [ebp + nb230_pos]       ;# base of pos[] 

	lea eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb230_ix]
	movapd xmm5, [esp + nb230_iy]
	movapd xmm6, [esp + nb230_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb230_dx], xmm4
	movapd [esp + nb230_dy], xmm5
	movapd [esp + nb230_dz], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	movapd xmm7, [esp + nb230_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb230_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb230_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [esp + nb230_crf]
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulsd  xmm7, [esp + nb230_two]
	movapd xmm1, xmm0
	subsd  xmm1, xmm7   ;# rinv-2*krsq
	mulsd  xmm1, xmm0   ;# (rinv-2*krsq)*rinv
	mulsd  xmm1, xmm3
	movsd [esp + nb230_fstmp], xmm1
	
	addsd  xmm6, [esp + nb230_vctot]
	movsd [esp + nb230_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, cmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb230_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl   ebx, 3

	mov  esi, [ebp + nb230_VFtab]

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
	mulsd  xmm7, [esp + nb230_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb230_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb230_Vvdwtot]
	movsd xmm3, [esp + nb230_fstmp]
	mulsd  xmm7, [esp + nb230_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb230_fstmp], xmm3
	movsd [esp + nb230_Vvdwtot], xmm5

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
	mulsd  xmm7, [esp + nb230_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb230_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb230_Vvdwtot]
	movsd xmm3, [esp + nb230_fstmp]
	mulsd  xmm7, [esp + nb230_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb230_Vvdwtot], xmm5

	mulsd  xmm3, xmm0
		
	movsd xmm0, [esp + nb230_dx]
	movsd xmm1, [esp + nb230_dy]
	movsd xmm2, [esp + nb230_dz]

	mov    edi, [ebp + nb230_faction]
	mulsd  xmm0, xmm3
	mulsd  xmm1, xmm3
	mulsd  xmm2, xmm3
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movlpd xmm3, [esp + nb230_fix]
	movlpd xmm4, [esp + nb230_fiy]
	movlpd xmm5, [esp + nb230_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb230_fix], xmm3
	movlpd [esp + nb230_fiy], xmm4
	movlpd [esp + nb230_fiz], xmm5
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5
	
.nb230_updateouterdata:
	mov   ecx, [esp + nb230_ii3]
	mov   edi, [ebp + nb230_faction]
	mov   esi, [ebp + nb230_fshift]
	mov   edx, [esp + nb230_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb230_fix]
	movapd xmm1, [esp + nb230_fiy]
	movapd xmm2, [esp + nb230_fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addpd  xmm0, xmm3
	addpd  xmm1, xmm4
	addpd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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

	;# increment fshift force  
	movsd  xmm3, [esi + edx*8]
	movsd  xmm4, [esi + edx*8 + 8]
	movsd  xmm5, [esi + edx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [esi + edx*8],     xmm3
	movsd  [esi + edx*8 + 8], xmm4
	movsd  [esi + edx*8 + 16], xmm5

	;# get n from stack
	mov esi, [esp + nb230_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb230_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb230_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb230_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb230_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb230_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb230_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb230_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb230_n], esi
        jmp .nb230_outer
.nb230_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb230_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb230_end
        ;# non-zero, do one more workunit
        jmp   .nb230_threadloop
.nb230_end:
	emms

	mov eax, [esp + nb230_nouter]
	mov ebx, [esp + nb230_ninner]
	mov ecx, [ebp + nb230_outeriter]
	mov edx, [ebp + nb230_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb230_salign]
	add esp, eax
	add esp,  392
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


.globl nb_kernel230nf_ia32_sse2
.globl _nb_kernel230nf_ia32_sse2
nb_kernel230nf_ia32_sse2:	
_nb_kernel230nf_ia32_sse2:	
.equiv          nb230nf_p_nri,            8
.equiv          nb230nf_iinr,             12
.equiv          nb230nf_jindex,           16
.equiv          nb230nf_jjnr,             20
.equiv          nb230nf_shift,            24
.equiv          nb230nf_shiftvec,         28
.equiv          nb230nf_fshift,           32
.equiv          nb230nf_gid,              36
.equiv          nb230nf_pos,              40
.equiv          nb230nf_faction,          44
.equiv          nb230nf_charge,           48
.equiv          nb230nf_p_facel,          52
.equiv          nb230nf_argkrf,           56
.equiv          nb230nf_argcrf,           60
.equiv          nb230nf_Vc,               64
.equiv          nb230nf_type,             68
.equiv          nb230nf_p_ntype,          72
.equiv          nb230nf_vdwparam,         76
.equiv          nb230nf_Vvdw,             80
.equiv          nb230nf_p_tabscale,       84
.equiv          nb230nf_VFtab,            88
.equiv          nb230nf_invsqrta,         92
.equiv          nb230nf_dvda,             96
.equiv          nb230nf_p_gbtabscale,     100
.equiv          nb230nf_GBtab,            104
.equiv          nb230nf_p_nthreads,       108
.equiv          nb230nf_count,            112
.equiv          nb230nf_mtx,              116
.equiv          nb230nf_outeriter,        120
.equiv          nb230nf_inneriter,        124
.equiv          nb230nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb230nf_ix,               0
.equiv          nb230nf_iy,               16
.equiv          nb230nf_iz,               32
.equiv          nb230nf_iq,               48
.equiv          nb230nf_dx,               64
.equiv          nb230nf_dy,               80
.equiv          nb230nf_dz,               96
.equiv          nb230nf_c6,               112
.equiv          nb230nf_c12,              128
.equiv          nb230nf_tsc,              144
.equiv          nb230nf_vctot,            176
.equiv          nb230nf_Vvdwtot,          192
.equiv          nb230nf_half,             256
.equiv          nb230nf_three,            272
.equiv          nb230nf_krf,              304
.equiv          nb230nf_crf,              320
.equiv          nb230nf_is3,              336
.equiv          nb230nf_ii3,              340
.equiv          nb230nf_ntia,             344
.equiv          nb230nf_innerjjnr,        348
.equiv          nb230nf_innerk,           352
.equiv          nb230nf_n,                356
.equiv          nb230nf_nn1,              360
.equiv          nb230nf_nri,              364
.equiv          nb230nf_facel,            368   ;# uses 8 bytes
.equiv          nb230nf_ntype,            376
.equiv          nb230nf_nouter,           380
.equiv          nb230nf_ninner,           384
.equiv          nb230nf_salign,           388
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  392		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb230nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb230nf_p_nri]
	mov esi, [ebp + nb230nf_p_facel]
	mov edi, [ebp + nb230nf_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb230nf_nri], ecx
	movsd [esp + nb230nf_facel], xmm7
	mov [esp + nb230nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb230nf_nouter], eax
	mov [esp + nb230nf_ninner], eax

	mov eax, [ebp + nb230nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb230nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb230nf_half], eax
	mov [esp + nb230nf_half+4], ebx
	movsd xmm1, [esp + nb230nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb230nf_half], xmm1
	movapd [esp + nb230nf_three], xmm3

	mov esi, [ebp + nb230nf_argkrf]
	mov edi, [ebp + nb230nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb230nf_krf], xmm5
	movapd [esp + nb230nf_crf], xmm6

.nb230nf_threadloop:
        mov   esi, [ebp + nb230nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb230nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb230nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb230nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb230nf_n], eax
        mov [esp + nb230nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb230nf_outerstart
        jmp .nb230nf_end

.nb230nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb230nf_nouter]
	mov [esp + nb230nf_nouter], ebx

.nb230nf_outer:
	mov   eax, [ebp + nb230nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb230nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb230nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb230nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb230nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb230nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb230nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb230nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb230nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb230nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb230nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb230nf_ix], xmm0
	movapd [esp + nb230nf_iy], xmm1
	movapd [esp + nb230nf_iz], xmm2

	mov   [esp + nb230nf_ii3], ebx
	
	;# clear vctot
	xorpd xmm4, xmm4
	movapd [esp + nb230nf_vctot], xmm4
	movapd [esp + nb230nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb230nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb230nf_pos]
	mov   edi, [ebp + nb230nf_faction]	
	mov   eax, [ebp + nb230nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb230nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb230nf_ninner]
	mov   [esp + nb230nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb230nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb230nf_unroll_loop
	jmp   .nb230nf_checksingle
.nb230nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb230nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb230nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb230nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb230nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb230nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb230nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb230nf_ntia]
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
	movapd [esp + nb230nf_c6], xmm4
	movapd [esp + nb230nf_c12], xmm6
	
	mov esi, [ebp + nb230nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb230nf_ix]
	movapd xmm5, [esp + nb230nf_iy]
	movapd xmm6, [esp + nb230nf_iz]

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
	;# rsq in xmm4 

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	movapd xmm7, [esp + nb230nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb230nf_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb230nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [esp + nb230nf_crf]
	mulpd  xmm6, xmm3
	
	addpd  xmm6, [esp + nb230nf_vctot]
	movapd [esp + nb230nf_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb230nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  esi, [ebp + nb230nf_VFtab]
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

	movapd xmm4, [esp + nb230nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;#Update Vvdwtot directly 
	addpd  xmm5, [esp + nb230nf_Vvdwtot]
	movapd [esp + nb230nf_Vvdwtot], xmm5

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
	
	movapd xmm4, [esp + nb230nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb230nf_Vvdwtot]
	movapd [esp + nb230nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb230nf_innerk],  2
	jl    .nb230nf_checksingle
	jmp   .nb230nf_unroll_loop

.nb230nf_checksingle:				
	mov   edx, [esp + nb230nf_innerk]
	and   edx, 1
	jnz    .nb230nf_dosingle
	jmp    .nb230nf_updateouterdata
.nb230nf_dosingle:			
	mov esi, [ebp + nb230nf_charge]
	mov edi, [ebp + nb230nf_pos]
	mov   ecx, [esp + nb230nf_innerjjnr]
	xorpd xmm3, xmm3
	mov   eax, [ecx]

	movlpd xmm3, [esi + eax*8]
	movapd xmm5, [esp + nb230nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb230nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb230nf_vdwparam]
	shl eax, 1
	mov edi, [esp + nb230nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		
	movapd [esp + nb230nf_c6], xmm4
	movapd [esp + nb230nf_c12], xmm6
	
	mov esi, [ebp + nb230nf_pos]       ;# base of pos[] 

	lea eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb230nf_ix]
	movapd xmm5, [esp + nb230nf_iy]
	movapd xmm6, [esp + nb230nf_iz]

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
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	movapd xmm7, [esp + nb230nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb230nf_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb230nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb230nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [esp + nb230nf_crf]
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	
	addsd  xmm6, [esp + nb230nf_vctot]
	movsd [esp + nb230nf_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, cmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb230nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl   ebx, 3

	mov  esi, [ebp + nb230nf_VFtab]

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

	movsd xmm4, [esp + nb230nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb230nf_Vvdwtot]
	movsd [esp + nb230nf_Vvdwtot], xmm5

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
	
	movsd xmm4, [esp + nb230nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb230nf_Vvdwtot]
	movsd [esp + nb230nf_Vvdwtot], xmm5
	
.nb230nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb230nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb230nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb230nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb230nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb230nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb230nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb230nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb230nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb230nf_n], esi
        jmp .nb230nf_outer
.nb230nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb230nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb230nf_end
        ;# non-zero, do one more workunit
        jmp   .nb230nf_threadloop
.nb230nf_end:
	emms

	mov eax, [esp + nb230nf_nouter]
	mov ebx, [esp + nb230nf_ninner]
	mov ecx, [ebp + nb230nf_outeriter]
	mov edx, [ebp + nb230nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb230nf_salign]
	add esp, eax
	add esp,  392
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

