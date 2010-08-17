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




.globl nb_kernel310_ia32_sse2
.globl _nb_kernel310_ia32_sse2
nb_kernel310_ia32_sse2:	
_nb_kernel310_ia32_sse2:	
.equiv          nb310_p_nri,            8
.equiv          nb310_iinr,             12
.equiv          nb310_jindex,           16
.equiv          nb310_jjnr,             20
.equiv          nb310_shift,            24
.equiv          nb310_shiftvec,         28
.equiv          nb310_fshift,           32
.equiv          nb310_gid,              36
.equiv          nb310_pos,              40
.equiv          nb310_faction,          44
.equiv          nb310_charge,           48
.equiv          nb310_p_facel,          52
.equiv          nb310_argkrf,           56
.equiv          nb310_argcrf,           60
.equiv          nb310_Vc,               64
.equiv          nb310_type,             68
.equiv          nb310_p_ntype,          72
.equiv          nb310_vdwparam,         76
.equiv          nb310_Vvdw,             80
.equiv          nb310_p_tabscale,       84
.equiv          nb310_VFtab,            88
.equiv          nb310_invsqrta,         92
.equiv          nb310_dvda,             96
.equiv          nb310_p_gbtabscale,     100
.equiv          nb310_GBtab,            104
.equiv          nb310_p_nthreads,       108
.equiv          nb310_count,            112
.equiv          nb310_mtx,              116
.equiv          nb310_outeriter,        120
.equiv          nb310_inneriter,        124
.equiv          nb310_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb310_ix,               0
.equiv          nb310_iy,               16
.equiv          nb310_iz,               32
.equiv          nb310_iq,               48
.equiv          nb310_dx,               64
.equiv          nb310_dy,               80
.equiv          nb310_dz,               96
.equiv          nb310_two,              112
.equiv          nb310_six,              128
.equiv          nb310_twelve,           144
.equiv          nb310_tsc,              160
.equiv          nb310_qq,               176
.equiv          nb310_c6,               192
.equiv          nb310_c12,              208
.equiv          nb310_fscal,            224
.equiv          nb310_vctot,            240
.equiv          nb310_Vvdwtot,          256
.equiv          nb310_fix,              272
.equiv          nb310_fiy,              288
.equiv          nb310_fiz,              304
.equiv          nb310_half,             320
.equiv          nb310_three,            336
.equiv          nb310_is3,              352
.equiv          nb310_ii3,              356
.equiv          nb310_ntia,             360
.equiv          nb310_innerjjnr,        364
.equiv          nb310_innerk,           368
.equiv          nb310_n,                372
.equiv          nb310_nn1,              376
.equiv          nb310_nri,              380
.equiv          nb310_facel,            384   ;# uses 8 bytes
.equiv          nb310_ntype,            392
.equiv          nb310_nouter,           396
.equiv          nb310_ninner,           400
.equiv          nb310_salign,           404
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 408		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb310_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb310_p_nri]
	mov esi, [ebp + nb310_p_facel]
	mov edi, [ebp + nb310_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb310_nri], ecx
	movsd [esp + nb310_facel], xmm7
	mov [esp + nb310_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb310_nouter], eax
	mov [esp + nb310_ninner], eax


	mov eax, [ebp + nb310_p_tabscale]
	movsd xmm5, [eax]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb310_tsc], xmm5
	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb310_half], eax
	mov [esp + nb310_half+4], ebx
	movsd xmm1, [esp + nb310_half]
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
	movapd [esp + nb310_half], xmm1
	movapd [esp + nb310_two], xmm2
	movapd [esp + nb310_three], xmm3
	movapd [esp + nb310_six], xmm4
	movapd [esp + nb310_twelve], xmm5

.nb310_threadloop:
        mov   esi, [ebp + nb310_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb310_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb310_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb310_n], eax
        mov [esp + nb310_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310_outerstart
        jmp .nb310_end

.nb310_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb310_nouter]
	mov [esp + nb310_nouter], ebx

.nb310_outer:
	mov   eax, [ebp + nb310_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb310_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb310_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb310_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb310_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb310_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb310_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb310_ntype]
    	shl   edx, 1
    	mov   [esp + nb310_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb310_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb310_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb310_ix], xmm0
	movapd [esp + nb310_iy], xmm1
	movapd [esp + nb310_iz], xmm2

	mov   [esp + nb310_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb310_vctot], xmm4
	movapd [esp + nb310_Vvdwtot], xmm4
	movapd [esp + nb310_fix], xmm4
	movapd [esp + nb310_fiy], xmm4
	movapd [esp + nb310_fiz], xmm4
	
	mov   eax, [ebp + nb310_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb310_pos]
	mov   edi, [ebp + nb310_faction]	
	mov   eax, [ebp + nb310_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb310_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb310_ninner]
	mov   [esp + nb310_ninner], ecx
	add   edx, 0
	mov   [esp + nb310_innerk], edx    ;# number of innerloop atoms 
	jge   .nb310_unroll_loop
	jmp   .nb310_checksingle
.nb310_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb310_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb310_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm2, [esp + nb310_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb310_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb310_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb310_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb310_ntia]
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
	movapd [esp + nb310_c6], xmm4
	movapd [esp + nb310_c12], xmm6
	
	mov esi, [ebp + nb310_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb310_ix]
	movapd xmm5, [esp + nb310_iy]
	movapd xmm6, [esp + nb310_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb310_dx], xmm4
	movapd [esp + nb310_dy], xmm5
	movapd [esp + nb310_dz], xmm6
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

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb310_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb310_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb310_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb310_VFtab]
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
	mulpd  xmm7, [esp + nb310_two]	;# two*Heps2 
	movapd xmm3, [esp + nb310_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	
	;# L-J 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [esp + nb310_vctot]

	movapd xmm6, xmm4
	mulpd  xmm6, xmm4

	movapd [esp + nb310_vctot], xmm5 

	mulpd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulpd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulpd  xmm6, [esp + nb310_c6]
	mulpd  xmm4, [esp + nb310_c12]
	movapd xmm7, [esp + nb310_Vvdwtot]
	addpd  xmm7, xmm4
	mulpd  xmm4, [esp + nb310_twelve]
	subpd  xmm7, xmm6
	mulpd  xmm3, [esp + nb310_tsc]
	mulpd  xmm6, [esp + nb310_six]
	movapd [esp + nb310_Vvdwtot], xmm7
	subpd  xmm4, xmm6
	mulpd  xmm4, xmm0
	subpd  xmm4, xmm3
	mulpd  xmm4, xmm0

	movapd xmm0, [esp + nb310_dx]
	movapd xmm1, [esp + nb310_dy]
	movapd xmm2, [esp + nb310_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb310_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb310_fix]
	movapd xmm4, [esp + nb310_fiy]
	movapd xmm5, [esp + nb310_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb310_fix], xmm3
	movapd [esp + nb310_fiy], xmm4
	movapd [esp + nb310_fiz], xmm5
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
	sub dword ptr [esp + nb310_innerk],  2
	jl    .nb310_checksingle
	jmp   .nb310_unroll_loop
.nb310_checksingle:
	mov   edx, [esp + nb310_innerk]
	and   edx, 1
	jnz    .nb310_dosingle
	jmp    .nb310_updateouterdata
.nb310_dosingle:
	mov esi, [ebp + nb310_charge]
	mov edi, [ebp + nb310_pos]
	mov   ecx, [esp + nb310_innerjjnr]
	mov   eax, [ecx]
	
	xorpd  xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm2, [esp + nb310_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb310_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb310_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb310_vdwparam]
	shl eax, 1
	mov edi, [esp + nb310_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb310_c6], xmm4
	movapd [esp + nb310_c12], xmm6
	
	mov esi, [ebp + nb310_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb310_ix]
	movapd xmm5, [esp + nb310_iy]
	movapd xmm6, [esp + nb310_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb310_dx], xmm4
	movapd [esp + nb310_dy], xmm5
	movapd [esp + nb310_dz], xmm6
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

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb310_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb310_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb310_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb310_VFtab]

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
	mulsd  xmm7, [esp + nb310_two]	;# two*Heps2 
	movapd xmm3, [esp + nb310_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	
	;# L-J 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [esp + nb310_vctot]

	movapd xmm6, xmm4
	mulsd  xmm6, xmm4

	movlpd [esp + nb310_vctot], xmm5 

	mulsd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulsd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulsd  xmm6, [esp + nb310_c6]
	mulsd  xmm4, [esp + nb310_c12]
	movapd xmm7, [esp + nb310_Vvdwtot]
	addsd  xmm7, xmm4
	mulsd  xmm4, [esp + nb310_twelve]
	subsd  xmm7, xmm6
	mulsd  xmm3, [esp + nb310_tsc]
	mulsd  xmm6, [esp + nb310_six]
	movlpd [esp + nb310_Vvdwtot], xmm7
	subsd  xmm4, xmm6
	mulsd  xmm4, xmm0
	subsd  xmm4, xmm3
	mulsd  xmm4, xmm0

	movapd xmm0, [esp + nb310_dx]
	movapd xmm1, [esp + nb310_dy]
	movapd xmm2, [esp + nb310_dz]

	movd eax, mm0	

	mov    edi, [ebp + nb310_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb310_fix]
	movapd xmm4, [esp + nb310_fiy]
	movapd xmm5, [esp + nb310_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb310_fix], xmm3
	movlpd [esp + nb310_fiy], xmm4
	movlpd [esp + nb310_fiz], xmm5
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
		
.nb310_updateouterdata:
	mov   ecx, [esp + nb310_ii3]
	mov   edi, [ebp + nb310_faction]
	mov   esi, [ebp + nb310_fshift]
	mov   edx, [esp + nb310_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb310_fix]
	movapd xmm1, [esp + nb310_fiy]
	movapd xmm2, [esp + nb310_fiz]

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
	mov esi, [esp + nb310_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb310_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb310_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb310_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb310_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb310_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb310_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb310_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb310_n], esi
        jmp .nb310_outer
.nb310_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb310_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb310_end
        ;# non-zero, do one more workunit
        jmp   .nb310_threadloop
.nb310_end:
	emms

	mov eax, [esp + nb310_nouter]
	mov ebx, [esp + nb310_ninner]
	mov ecx, [ebp + nb310_outeriter]
	mov edx, [ebp + nb310_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb310_salign]
	add esp, eax
	add esp, 408
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel310nf_ia32_sse2
.globl _nb_kernel310nf_ia32_sse2
nb_kernel310nf_ia32_sse2:	
_nb_kernel310nf_ia32_sse2:	
.equiv          nb310nf_p_nri,          8
.equiv          nb310nf_iinr,           12
.equiv          nb310nf_jindex,         16
.equiv          nb310nf_jjnr,           20
.equiv          nb310nf_shift,          24
.equiv          nb310nf_shiftvec,       28
.equiv          nb310nf_fshift,         32
.equiv          nb310nf_gid,            36
.equiv          nb310nf_pos,            40
.equiv          nb310nf_faction,        44
.equiv          nb310nf_charge,         48
.equiv          nb310nf_p_facel,        52
.equiv          nb310nf_argkrf,         56
.equiv          nb310nf_argcrf,         60
.equiv          nb310nf_Vc,             64
.equiv          nb310nf_type,           68
.equiv          nb310nf_p_ntype,        72
.equiv          nb310nf_vdwparam,       76
.equiv          nb310nf_Vvdw,           80
.equiv          nb310nf_p_tabscale,     84
.equiv          nb310nf_VFtab,          88
.equiv          nb310nf_invsqrta,       92
.equiv          nb310nf_dvda,           96
.equiv          nb310nf_p_gbtabscale,   100
.equiv          nb310nf_GBtab,          104
.equiv          nb310nf_p_nthreads,     108
.equiv          nb310nf_count,          112
.equiv          nb310nf_mtx,            116
.equiv          nb310nf_outeriter,      120
.equiv          nb310nf_inneriter,      124
.equiv          nb310nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb310nf_ix,             0
.equiv          nb310nf_iy,             16
.equiv          nb310nf_iz,             32
.equiv          nb310nf_iq,             48
.equiv          nb310nf_tsc,            64
.equiv          nb310nf_qq,             80
.equiv          nb310nf_c6,             96
.equiv          nb310nf_c12,            112
.equiv          nb310nf_vctot,          128
.equiv          nb310nf_Vvdwtot,        144
.equiv          nb310nf_half,           160
.equiv          nb310nf_three,          176
.equiv          nb310nf_is3,            192
.equiv          nb310nf_ii3,            196
.equiv          nb310nf_ntia,           200
.equiv          nb310nf_innerjjnr,      204
.equiv          nb310nf_innerk,         208
.equiv          nb310nf_n,              212
.equiv          nb310nf_nn1,            216
.equiv          nb310nf_nri,            220
.equiv          nb310nf_facel,          224   ;# uses 8 bytes
.equiv          nb310nf_ntype,          232
.equiv          nb310nf_nouter,         236
.equiv          nb310nf_ninner,         240
.equiv          nb310nf_salign,         244
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 248		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb310nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb310nf_p_nri]
	mov esi, [ebp + nb310nf_p_facel]
	mov edi, [ebp + nb310nf_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb310nf_nri], ecx
	movsd [esp + nb310nf_facel], xmm7
	mov [esp + nb310nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb310nf_nouter], eax
	mov [esp + nb310nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb310nf_half], eax
	mov [esp + nb310nf_half+4], ebx
	movsd xmm1, [esp + nb310nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb310nf_half], xmm1
	movapd [esp + nb310nf_three], xmm3
	mov eax, [ebp + nb310nf_p_tabscale]
	movsd xmm5, [eax]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb310nf_tsc], xmm5

.nb310nf_threadloop:
        mov   esi, [ebp + nb310nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb310nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb310nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb310nf_n], eax
        mov [esp + nb310nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310nf_outerstart
        jmp .nb310nf_end

.nb310nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb310nf_nouter]
	mov [esp + nb310nf_nouter], ebx

.nb310nf_outer:
	mov   eax, [ebp + nb310nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb310nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb310nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb310nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb310nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb310nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb310nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb310nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb310nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb310nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb310nf_ix], xmm0
	movapd [esp + nb310nf_iy], xmm1
	movapd [esp + nb310nf_iz], xmm2

	mov   [esp + nb310nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb310nf_vctot], xmm4
	movapd [esp + nb310nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb310nf_jindex]
	mov   ecx, [eax+esi*4]	     		;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     	;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb310nf_pos]
	mov   eax, [ebp + nb310nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb310nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb310nf_ninner]
	mov   [esp + nb310nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb310nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb310nf_unroll_loop
	jmp   .nb310nf_checksingle
.nb310nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb310nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb310nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb310nf_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm2, [esp + nb310nf_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb310nf_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb310nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb310nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb310nf_ntia]
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
	movapd [esp + nb310nf_c6], xmm4
	movapd [esp + nb310nf_c12], xmm6
	
	mov esi, [ebp + nb310nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb310nf_ix]
	movapd xmm5, [esp + nb310nf_iy]
	movapd xmm6, [esp + nb310nf_iz]

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

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb310nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb310nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb310nf_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb310nf_VFtab]
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
	movapd xmm3, [esp + nb310nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	
	;# L-J 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [esp + nb310nf_vctot]

	movapd xmm6, xmm4
	mulpd  xmm6, xmm4

	movapd [esp + nb310nf_vctot], xmm5 

	mulpd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulpd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulpd  xmm6, [esp + nb310nf_c6]
	mulpd  xmm4, [esp + nb310nf_c12]
	movapd xmm7, [esp + nb310nf_Vvdwtot]
	addpd  xmm7, xmm4
	subpd  xmm7, xmm6
	movapd [esp + nb310nf_Vvdwtot], xmm7
		
	;# should we do one more iteration? 
	sub dword ptr [esp + nb310nf_innerk],  2
	jl    .nb310nf_checksingle
	jmp   .nb310nf_unroll_loop
.nb310nf_checksingle:
	mov   edx, [esp + nb310nf_innerk]
	and   edx, 1
	jnz    .nb310nf_dosingle
	jmp    .nb310nf_updateouterdata
.nb310nf_dosingle:
	mov esi, [ebp + nb310nf_charge]
	mov edi, [ebp + nb310nf_pos]
	mov   ecx, [esp + nb310nf_innerjjnr]
	mov   eax, [ecx]
	
	xorpd  xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm2, [esp + nb310nf_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb310nf_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb310nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb310nf_vdwparam]
	shl eax, 1
	mov edi, [esp + nb310nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb310nf_c6], xmm4
	movapd [esp + nb310nf_c12], xmm6
	
	mov esi, [ebp + nb310nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb310nf_ix]
	movapd xmm5, [esp + nb310nf_iy]
	movapd xmm6, [esp + nb310nf_iz]

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

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb310nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb310nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb310nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb310nf_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb310nf_VFtab]

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
	movapd xmm3, [esp + nb310nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	
	;# L-J 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [esp + nb310nf_vctot]

	movapd xmm6, xmm4
	mulsd  xmm6, xmm4

	movlpd [esp + nb310nf_vctot], xmm5 

	mulsd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulsd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulsd  xmm6, [esp + nb310nf_c6]
	mulsd  xmm4, [esp + nb310nf_c12]
	movapd xmm7, [esp + nb310nf_Vvdwtot]
	addsd  xmm7, xmm4
	subsd  xmm7, xmm6
	movlpd [esp + nb310nf_Vvdwtot], xmm7
		
.nb310nf_updateouterdata:
	;# get group index for i particle 
	;# get n from stack
	mov esi, [esp + nb310nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb310nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb310nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb310nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb310nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb310nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb310nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb310nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb310nf_n], esi
        jmp .nb310nf_outer
.nb310nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb310nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb310nf_end
        ;# non-zero, do one more workunit
        jmp   .nb310nf_threadloop
.nb310nf_end:
	emms

	mov eax, [esp + nb310nf_nouter]
	mov ebx, [esp + nb310nf_ninner]
	mov ecx, [ebp + nb310nf_outeriter]
	mov edx, [ebp + nb310nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb310nf_salign]
	add esp, eax
	add esp, 248
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


