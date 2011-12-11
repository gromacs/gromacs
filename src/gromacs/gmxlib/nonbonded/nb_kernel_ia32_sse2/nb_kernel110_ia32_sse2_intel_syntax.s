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





.globl nb_kernel110_ia32_sse2
.globl _nb_kernel110_ia32_sse2
nb_kernel110_ia32_sse2:	
_nb_kernel110_ia32_sse2:	
.equiv          nb110_p_nri,            8
.equiv          nb110_iinr,             12
.equiv          nb110_jindex,           16
.equiv          nb110_jjnr,             20
.equiv          nb110_shift,            24
.equiv          nb110_shiftvec,         28
.equiv          nb110_fshift,           32
.equiv          nb110_gid,              36
.equiv          nb110_pos,              40
.equiv          nb110_faction,          44
.equiv          nb110_charge,           48
.equiv          nb110_p_facel,          52
.equiv          nb110_argkrf,           56
.equiv          nb110_argcrf,           60
.equiv          nb110_Vc,               64
.equiv          nb110_type,             68
.equiv          nb110_p_ntype,          72
.equiv          nb110_vdwparam,         76
.equiv          nb110_Vvdw,             80
.equiv          nb110_p_tabscale,       84
.equiv          nb110_VFtab,            88
.equiv          nb110_invsqrta,         92
.equiv          nb110_dvda,             96
.equiv          nb110_p_gbtabscale,     100
.equiv          nb110_GBtab,            104
.equiv          nb110_p_nthreads,       108
.equiv          nb110_count,            112
.equiv          nb110_mtx,              116
.equiv          nb110_outeriter,        120
.equiv          nb110_inneriter,        124
.equiv          nb110_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb110_ix,               0
.equiv          nb110_iy,               16
.equiv          nb110_iz,               32
.equiv          nb110_iq,               48
.equiv          nb110_dx,               64
.equiv          nb110_dy,               80
.equiv          nb110_dz,               96
.equiv          nb110_c6,               112
.equiv          nb110_c12,              128
.equiv          nb110_six,              144
.equiv          nb110_twelve,           160
.equiv          nb110_vctot,            176
.equiv          nb110_Vvdwtot,          192
.equiv          nb110_fix,              208
.equiv          nb110_fiy,              224
.equiv          nb110_fiz,              240
.equiv          nb110_half,             256
.equiv          nb110_three,            272
.equiv          nb110_facel,            288
.equiv          nb110_is3,              304
.equiv          nb110_ii3,              308
.equiv          nb110_ntia,             312
.equiv          nb110_innerjjnr,        316
.equiv          nb110_innerk,           320
.equiv          nb110_n,                324
.equiv          nb110_nn1,              328
.equiv          nb110_nri,              332
.equiv          nb110_ntype,            336
.equiv          nb110_nouter,           340
.equiv          nb110_ninner,           344
.equiv          nb110_salign,           348
	
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  352		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb110_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb110_p_nri]
	mov esi, [ebp + nb110_p_facel]
	mov edi, [ebp + nb110_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb110_nri], ecx
	movsd [esp + nb110_facel], xmm7
	mov [esp + nb110_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb110_nouter], eax
	mov [esp + nb110_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb110_half], eax
	mov [esp + nb110_half+4], ebx
	movsd xmm1, [esp + nb110_half]
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
	movapd [esp + nb110_half], xmm1
	movapd [esp + nb110_three], xmm3
	movapd [esp + nb110_six], xmm4
	movapd [esp + nb110_twelve], xmm5

.nb110_threadloop:
        mov   esi, [ebp + nb110_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb110_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb110_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb110_n], eax
        mov [esp + nb110_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110_outerstart
        jmp .nb110_end

.nb110_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb110_nouter]
	mov [esp + nb110_nouter], ebx

.nb110_outer:
	mov   eax, [ebp + nb110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb110_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb110_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb110_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb110_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb110_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb110_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb110_ntype]
    	shl   edx, 1
    	mov   [esp + nb110_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb110_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb110_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb110_ix], xmm0
	movapd [esp + nb110_iy], xmm1
	movapd [esp + nb110_iz], xmm2

	mov   [esp + nb110_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb110_vctot], xmm4
	movapd [esp + nb110_Vvdwtot], xmm4
	movapd [esp + nb110_fix], xmm4
	movapd [esp + nb110_fiy], xmm4
	movapd [esp + nb110_fiz], xmm4
	
	mov   eax, [ebp + nb110_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb110_pos]
	mov   edi, [ebp + nb110_faction]	
	mov   eax, [ebp + nb110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb110_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb110_ninner]
	mov   [esp + nb110_ninner], ecx
	add   edx, 0
	mov   [esp + nb110_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110_unroll_loop
	jmp   .nb110_checksingle
.nb110_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb110_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb110_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb110_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb110_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb110_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb110_ntia]
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
	movapd [esp + nb110_c6], xmm4
	movapd [esp + nb110_c12], xmm6
	
	mov esi, [ebp + nb110_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb110_ix]
	movapd xmm5, [esp + nb110_iy]
	movapd xmm6, [esp + nb110_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb110_dx], xmm4
	movapd [esp + nb110_dy], xmm5
	movapd [esp + nb110_dz], xmm6
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
	movapd xmm1, [esp + nb110_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb110_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	mulpd  xmm1, [esp + nb110_c6]
	mulpd  xmm2, [esp + nb110_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [esp + nb110_Vvdwtot]
	mulpd  xmm1, [esp + nb110_six]
	mulpd  xmm2, [esp + nb110_twelve]
	subpd  xmm2, xmm1
	addpd  xmm2, xmm3
	mulpd  xmm4, xmm2	;# xmm4=total fscal 
	addpd  xmm3, [esp + nb110_vctot]

	movapd xmm0, [esp + nb110_dx]
	movapd xmm1, [esp + nb110_dy]
	movapd xmm2, [esp + nb110_dz]

	movapd [esp + nb110_vctot], xmm3
	movapd [esp + nb110_Vvdwtot], xmm5

	mov    edi, [ebp + nb110_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb110_fix]
	movapd xmm4, [esp + nb110_fiy]
	movapd xmm5, [esp + nb110_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb110_fix], xmm3
	movapd [esp + nb110_fiy], xmm4
	movapd [esp + nb110_fiz], xmm5
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
	sub dword ptr [esp + nb110_innerk],  2
	jl    .nb110_checksingle
	jmp   .nb110_unroll_loop	
.nb110_checksingle:
	mov   edx, [esp + nb110_innerk]
	and   edx, 1
	jnz    .nb110_dosingle
	jmp    .nb110_updateouterdata
.nb110_dosingle:
	mov esi, [ebp + nb110_charge]
	mov edi, [ebp + nb110_pos]
	mov ecx, [esp + nb110_innerjjnr]
	mov   eax, [ecx]
	
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]

	movapd xmm5, [esp + nb110_iq]
	mulsd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov esi, [ebp + nb110_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb110_vdwparam]
	shl eax, 1
	mov edi, [esp + nb110_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	
	movapd [esp + nb110_c6], xmm4
	movapd [esp + nb110_c12], xmm6
	
	mov esi, [ebp + nb110_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb110_ix]
	movapd xmm5, [esp + nb110_iy]
	movapd xmm6, [esp + nb110_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb110_dx], xmm4
	movapd [esp + nb110_dy], xmm5
	movapd [esp + nb110_dz], xmm6
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
	movapd xmm1, [esp + nb110_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb110_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	mulsd  xmm1, [esp + nb110_c6]
	mulsd  xmm2, [esp + nb110_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [esp + nb110_Vvdwtot]
	mulsd  xmm1, [esp + nb110_six]
	mulsd  xmm2, [esp + nb110_twelve]
	subsd  xmm2, xmm1
	addsd  xmm2, xmm3
	mulsd  xmm4, xmm2	;# xmm4=total fscal 
	addsd  xmm3, [esp + nb110_vctot]

	movapd xmm0, [esp + nb110_dx]
	movapd xmm1, [esp + nb110_dy]
	movapd xmm2, [esp + nb110_dz]

	movlpd [esp + nb110_vctot], xmm3
	movlpd [esp + nb110_Vvdwtot], xmm5

	mov    edi, [ebp + nb110_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movlpd xmm3, [esp + nb110_fix]
	movlpd xmm4, [esp + nb110_fiy]
	movlpd xmm5, [esp + nb110_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb110_fix], xmm3
	movlpd [esp + nb110_fiy], xmm4
	movlpd [esp + nb110_fiz], xmm5
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
	
.nb110_updateouterdata:
	mov   ecx, [esp + nb110_ii3]
	mov   edi, [ebp + nb110_faction]
	mov   esi, [ebp + nb110_fshift]
	mov   edx, [esp + nb110_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb110_fix]
	movapd xmm1, [esp + nb110_fiy]
	movapd xmm2, [esp + nb110_fiz]

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
	mov esi, [esp + nb110_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb110_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb110_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb110_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb110_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb110_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb110_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb110_n], esi
        jmp .nb110_outer
.nb110_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb110_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110_end
        ;# non-zero, do one more workunit
        jmp   .nb110_threadloop
.nb110_end:
	emms

	mov eax, [esp + nb110_nouter]
	mov ebx, [esp + nb110_ninner]
	mov ecx, [ebp + nb110_outeriter]
	mov edx, [ebp + nb110_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb110_salign]
	add esp, eax
	add esp,  352
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret






.globl nb_kernel110nf_ia32_sse2
.globl _nb_kernel110nf_ia32_sse2
nb_kernel110nf_ia32_sse2:	
_nb_kernel110nf_ia32_sse2:	
.equiv          nb110nf_p_nri,          8
.equiv          nb110nf_iinr,           12
.equiv          nb110nf_jindex,         16
.equiv          nb110nf_jjnr,           20
.equiv          nb110nf_shift,          24
.equiv          nb110nf_shiftvec,       28
.equiv          nb110nf_fshift,         32
.equiv          nb110nf_gid,            36
.equiv          nb110nf_pos,            40
.equiv          nb110nf_faction,        44
.equiv          nb110nf_charge,         48
.equiv          nb110nf_p_facel,        52
.equiv          nb110nf_argkrf,         56
.equiv          nb110nf_argcrf,         60
.equiv          nb110nf_Vc,             64
.equiv          nb110nf_type,           68
.equiv          nb110nf_p_ntype,        72
.equiv          nb110nf_vdwparam,       76
.equiv          nb110nf_Vvdw,           80
.equiv          nb110nf_p_tabscale,     84
.equiv          nb110nf_VFtab,          88
.equiv          nb110nf_invsqrta,       92
.equiv          nb110nf_dvda,           96
.equiv          nb110nf_p_gbtabscale,   100
.equiv          nb110nf_GBtab,          104
.equiv          nb110nf_p_nthreads,     108
.equiv          nb110nf_count,          112
.equiv          nb110nf_mtx,            116
.equiv          nb110nf_outeriter,      120
.equiv          nb110nf_inneriter,      124
.equiv          nb110nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb110nf_ix,             0
.equiv          nb110nf_iy,             16
.equiv          nb110nf_iz,             32
.equiv          nb110nf_iq,             48
.equiv          nb110nf_c6,             64
.equiv          nb110nf_c12,            80
.equiv          nb110nf_vctot,          96
.equiv          nb110nf_Vvdwtot,        112
.equiv          nb110nf_half,           128
.equiv          nb110nf_three,          144
.equiv          nb110nf_is3,            160
.equiv          nb110nf_ii3,            164
.equiv          nb110nf_ntia,           168
.equiv          nb110nf_innerjjnr,      172
.equiv          nb110nf_innerk,         176
.equiv          nb110nf_n,              180
.equiv          nb110nf_nn1,            184
.equiv          nb110nf_nri,            188
.equiv          nb110nf_facel,          192   ;# uses 8 bytes
.equiv          nb110nf_ntype,          200
.equiv          nb110nf_nouter,         204
.equiv          nb110nf_ninner,         208
.equiv          nb110nf_salign,         212
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  192		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb110nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb110nf_p_nri]
	mov esi, [ebp + nb110nf_p_facel]
	mov edi, [ebp + nb110nf_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb110nf_nri], ecx
	movsd [esp + nb110nf_facel], xmm7
	mov [esp + nb110nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb110nf_nouter], eax
	mov [esp + nb110nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb110nf_half], eax
	mov [esp + nb110nf_half+4], ebx
	movsd xmm1, [esp + nb110nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb110nf_half], xmm1
	movapd [esp + nb110nf_three], xmm3

.nb110nf_threadloop:
        mov   esi, [ebp + nb110nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb110nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb110nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb110nf_n], eax
        mov [esp + nb110nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110nf_outerstart
        jmp .nb110nf_end

.nb110nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb110nf_nouter]
	mov [esp + nb110nf_nouter], ebx

.nb110nf_outer:
	mov   eax, [ebp + nb110nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb110nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb110nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb110nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb110nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb110nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb110nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb110nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb110nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb110nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb110nf_ix], xmm0
	movapd [esp + nb110nf_iy], xmm1
	movapd [esp + nb110nf_iz], xmm2

	mov   [esp + nb110nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb110nf_vctot], xmm4
	movapd [esp + nb110nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb110nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb110nf_pos]
	mov   eax, [ebp + nb110nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb110nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb110nf_ninner]
	mov   [esp + nb110nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb110nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110nf_unroll_loop
	jmp   .nb110nf_checksingle
.nb110nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb110nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb110nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb110nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb110nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb110nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb110nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb110nf_ntia]
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
	movapd [esp + nb110nf_c6], xmm4
	movapd [esp + nb110nf_c12], xmm6
	
	mov esi, [ebp + nb110nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb110nf_ix]
	movapd xmm5, [esp + nb110nf_iy]
	movapd xmm6, [esp + nb110nf_iz]

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
	movapd xmm1, [esp + nb110nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb110nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	mulpd  xmm1, [esp + nb110nf_c6]
	mulpd  xmm2, [esp + nb110nf_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [esp + nb110nf_Vvdwtot]
	addpd  xmm3, [esp + nb110nf_vctot]
	movapd [esp + nb110nf_vctot], xmm3
	movapd [esp + nb110nf_Vvdwtot], xmm5
		
	;# should we do one more iteration? 
	sub dword ptr [esp + nb110nf_innerk],  2
	jl    .nb110nf_checksingle
	jmp   .nb110nf_unroll_loop	
.nb110nf_checksingle:
	mov   edx, [esp + nb110nf_innerk]
	and   edx, 1
	jnz   .nb110nf_dosingle
	jmp   .nb110nf_updateouterdata
.nb110nf_dosingle:
	mov esi, [ebp + nb110nf_charge]
	mov edi, [ebp + nb110nf_pos]
	mov ecx, [esp + nb110nf_innerjjnr]
	mov   eax, [ecx]
	
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]

	movapd xmm5, [esp + nb110nf_iq]
	mulsd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov esi, [ebp + nb110nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb110nf_vdwparam]
	shl eax, 1
	mov edi, [esp + nb110nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	
	movapd [esp + nb110nf_c6], xmm4
	movapd [esp + nb110nf_c12], xmm6
	
	mov esi, [ebp + nb110nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb110nf_ix]
	movapd xmm5, [esp + nb110nf_iy]
	movapd xmm6, [esp + nb110nf_iz]

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
	movapd xmm1, [esp + nb110nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb110nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb110nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	mulsd  xmm1, [esp + nb110nf_c6]
	mulsd  xmm2, [esp + nb110nf_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [esp + nb110nf_Vvdwtot]
	addsd  xmm3, [esp + nb110nf_vctot]
	movlpd [esp + nb110nf_vctot], xmm3
	movlpd [esp + nb110nf_Vvdwtot], xmm5
	
.nb110nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb110nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb110nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb110nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb110nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb110nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb110nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb110nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb110nf_n], esi
        jmp .nb110nf_outer
.nb110nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb110nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110nf_end
        ;# non-zero, do one more workunit
        jmp   .nb110nf_threadloop
.nb110nf_end:
	emms

	mov eax, [esp + nb110nf_nouter]
	mov ebx, [esp + nb110nf_ninner]
	mov ecx, [ebp + nb110nf_outeriter]
	mov edx, [ebp + nb110nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb110nf_salign]
	add esp, eax
	add esp,  192
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


