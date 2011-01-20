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




.globl nb_kernel100_ia32_sse2
.globl _nb_kernel100_ia32_sse2
nb_kernel100_ia32_sse2:	
_nb_kernel100_ia32_sse2:	
.equiv          nb100_p_nri,            8
.equiv          nb100_iinr,             12
.equiv          nb100_jindex,           16
.equiv          nb100_jjnr,             20
.equiv          nb100_shift,            24
.equiv          nb100_shiftvec,         28
.equiv          nb100_fshift,           32
.equiv          nb100_gid,              36
.equiv          nb100_pos,              40
.equiv          nb100_faction,          44
.equiv          nb100_charge,           48
.equiv          nb100_p_facel,          52
.equiv          nb100_argkrf,           56
.equiv          nb100_argcrf,           60
.equiv          nb100_Vc,               64
.equiv          nb100_type,             68
.equiv          nb100_p_ntype,          72
.equiv          nb100_vdwparam,         76
.equiv          nb100_Vvdw,             80
.equiv          nb100_p_tabscale,       84
.equiv          nb100_VFtab,            88
.equiv          nb100_invsqrta,         92
.equiv          nb100_dvda,             96
.equiv          nb100_p_gbtabscale,     100
.equiv          nb100_GBtab,            104
.equiv          nb100_p_nthreads,       108
.equiv          nb100_count,            112
.equiv          nb100_mtx,              116
.equiv          nb100_outeriter,        120
.equiv          nb100_inneriter,        124
.equiv          nb100_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb100_ix,               0
.equiv          nb100_iy,               16
.equiv          nb100_iz,               32
.equiv          nb100_iq,               48
.equiv          nb100_dx,               64
.equiv          nb100_dy,               80
.equiv          nb100_dz,               96
.equiv          nb100_vctot,            112
.equiv          nb100_fix,              128
.equiv          nb100_fiy,              144
.equiv          nb100_fiz,              160
.equiv          nb100_half,             176
.equiv          nb100_three,            192
.equiv          nb100_is3,              208
.equiv          nb100_ii3,              212
.equiv          nb100_innerjjnr,        216
.equiv          nb100_innerk,           220
.equiv          nb100_n,                224
.equiv          nb100_nn1,              228
.equiv          nb100_nri,              232
.equiv          nb100_facel,            240   ;# uses 8 bytes
.equiv          nb100_nouter,           248
.equiv          nb100_ninner,           252
.equiv          nb100_salign,           256
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 260		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb100_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb100_p_nri]
	mov esi, [ebp + nb100_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb100_nri], ecx
	movsd [esp + nb100_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb100_nouter], eax
	mov [esp + nb100_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb100_half], eax
	mov [esp + nb100_half+4], ebx
	movsd xmm1, [esp + nb100_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb100_half], xmm1
	movapd [esp + nb100_three], xmm3

.nb100_threadloop:
        mov   esi, [ebp + nb100_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb100_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb100_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb100_n], eax
        mov [esp + nb100_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb100_outerstart
        jmp .nb100_end

.nb100_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb100_nouter]
	mov [esp + nb100_nouter], ebx

.nb100_outer:
	mov   eax, [ebp + nb100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb100_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb100_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb100_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb100_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb100_facel]
	shufpd xmm3, xmm3, 0	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb100_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb100_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb100_ix], xmm0
	movapd [esp + nb100_iy], xmm1
	movapd [esp + nb100_iz], xmm2

	mov   [esp + nb100_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb100_vctot], xmm4
	movapd [esp + nb100_fix], xmm4
	movapd [esp + nb100_fiy], xmm4
	movapd [esp + nb100_fiz], xmm4
	
	mov   eax, [ebp + nb100_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb100_pos]
	mov   edi, [ebp + nb100_faction]	
	mov   eax, [ebp + nb100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb100_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb100_ninner]
	mov   [esp + nb100_ninner], ecx
	add   edx, 0
	mov   [esp + nb100_innerk], edx    ;# number of innerloop atoms 
	jge   .nb100_unroll_loop
	jmp   .nb100_checksingle
.nb100_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb100_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb100_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]	;# jq A 
	movhpd xmm3, [esi + ebx*8]	;# jq B 

	movapd xmm5, [esp + nb100_iq]
	
	mulpd xmm3, xmm5		;# qq 
	
	mov esi, [ebp + nb100_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	mov    edi, [ebp + nb100_faction]
	
	;# move nb100_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb100_ix]
	movapd xmm5, [esp + nb100_iy]
	movapd xmm6, [esp + nb100_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb100_dx], xmm4
	movapd [esp + nb100_dy], xmm5
	movapd [esp + nb100_dz], xmm6
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
	movapd xmm1, [esp + nb100_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb100_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 

	
	movapd xmm5, [esp + nb100_vctot]
	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	mulpd  xmm4, xmm3	;# xmm4=fscal 
	addpd  xmm5, xmm3

	movapd xmm0, [esp + nb100_dx]
	movapd xmm1, [esp + nb100_dy]
	movapd xmm2, [esp + nb100_dz]

	movapd [esp + nb100_vctot], xmm5

	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb100_fix]
	movapd xmm4, [esp + nb100_fiy]
	movapd xmm5, [esp + nb100_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb100_fix], xmm3
	movapd [esp + nb100_fiy], xmm4
	movapd [esp + nb100_fiz], xmm5
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
	sub dword ptr [esp + nb100_innerk],  2
	jl    .nb100_checksingle
	jmp   .nb100_unroll_loop

.nb100_checksingle:				
	mov   edx, [esp + nb100_innerk]
	and   edx, 1
	jnz    .nb100_dosingle
	jmp    .nb100_updateouterdata
.nb100_dosingle:			
	mov esi, [ebp + nb100_charge]
	mov edi, [ebp + nb100_pos]

	mov edx, [esp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov eax, [edx]	

	xorpd xmm3, xmm3
	movsd xmm3, [esi + eax*8]	;# jq A 
	movapd xmm5, [esp + nb100_iq]
	unpcklpd xmm3, xmm6
	mulpd xmm3, xmm5		;# qq 
	
	mov esi, [ebp + nb100_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	mov    edi, [ebp + nb100_faction]

	;# move nb100_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb100_ix]
	movapd xmm5, [esp + nb100_iy]
	movapd xmm6, [esp + nb100_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movlpd [esp + nb100_dx], xmm4
	movlpd [esp + nb100_dy], xmm5
	movlpd [esp + nb100_dz], xmm6
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
	movapd xmm1, [esp + nb100_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb100_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 

	movlpd xmm5, [esp + nb100_vctot]
	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	mulsd  xmm4, xmm3	;# xmm4=fscal 
	addsd  xmm5, xmm3

	movapd xmm0, [esp + nb100_dx]
	movapd xmm1, [esp + nb100_dy]
	movapd xmm2, [esp + nb100_dz]

	movlpd [esp + nb100_vctot], xmm5

	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movlpd xmm3, [esp + nb100_fix]
	movlpd xmm4, [esp + nb100_fiy]
	movlpd xmm5, [esp + nb100_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb100_fix], xmm3
	movlpd [esp + nb100_fiy], xmm4
	movlpd [esp + nb100_fiz], xmm5
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

.nb100_updateouterdata:
	mov   ecx, [esp + nb100_ii3]
	mov   edi, [ebp + nb100_faction]
	mov   esi, [ebp + nb100_fshift]
	mov   edx, [esp + nb100_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb100_fix]
	movapd xmm1, [esp + nb100_fiy]
	movapd xmm2, [esp + nb100_fiz]

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
	mov esi, [esp + nb100_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb100_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb100_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb100_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb100_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb100_n], esi
        jmp .nb100_outer
.nb100_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb100_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb100_end
        ;# non-zero, do one more workunit
        jmp   .nb100_threadloop
.nb100_end:
	emms

	mov eax, [esp + nb100_nouter]
	mov ebx, [esp + nb100_ninner]
	mov ecx, [ebp + nb100_outeriter]
	mov edx, [ebp + nb100_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb100_salign]
	add esp, eax
	add esp, 260
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel100nf_ia32_sse2
.globl _nb_kernel100nf_ia32_sse2
nb_kernel100nf_ia32_sse2:	
_nb_kernel100nf_ia32_sse2:	
.equiv          nb100nf_p_nri,          8
.equiv          nb100nf_iinr,           12
.equiv          nb100nf_jindex,         16
.equiv          nb100nf_jjnr,           20
.equiv          nb100nf_shift,          24
.equiv          nb100nf_shiftvec,       28
.equiv          nb100nf_fshift,         32
.equiv          nb100nf_gid,            36
.equiv          nb100nf_pos,            40
.equiv          nb100nf_faction,        44
.equiv          nb100nf_charge,         48
.equiv          nb100nf_p_facel,        52
.equiv          nb100nf_argkrf,         56
.equiv          nb100nf_argcrf,         60
.equiv          nb100nf_Vc,             64
.equiv          nb100nf_type,           68
.equiv          nb100nf_p_ntype,        72
.equiv          nb100nf_vdwparam,       76
.equiv          nb100nf_Vvdw,           80
.equiv          nb100nf_p_tabscale,     84
.equiv          nb100nf_VFtab,          88
.equiv          nb100nf_invsqrta,       92
.equiv          nb100nf_dvda,           96
.equiv          nb100nf_p_gbtabscale,   100
.equiv          nb100nf_GBtab,          104
.equiv          nb100nf_p_nthreads,     108
.equiv          nb100nf_count,          112
.equiv          nb100nf_mtx,            116
.equiv          nb100nf_outeriter,      120
.equiv          nb100nf_inneriter,      124
.equiv          nb100nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb100nf_ix,             0
.equiv          nb100nf_iy,             16
.equiv          nb100nf_iz,             32
.equiv          nb100nf_iq,             48
.equiv          nb100nf_vctot,          64
.equiv          nb100nf_half,           80
.equiv          nb100nf_three,          96
.equiv          nb100nf_is3,            112
.equiv          nb100nf_ii3,            116
.equiv          nb100nf_innerjjnr,      120
.equiv          nb100nf_innerk,         124
.equiv          nb100nf_n,              128
.equiv          nb100nf_nn1,            132
.equiv          nb100nf_nri,            136
.equiv          nb100nf_facel,          144   ;# uses 8 bytes
.equiv          nb100nf_nouter,         152
.equiv          nb100nf_ninner,         156
.equiv          nb100nf_salign,         160
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 164		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb100nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb100nf_p_nri]
	mov esi, [ebp + nb100nf_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb100nf_nri], ecx
	movsd [esp + nb100nf_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb100nf_nouter], eax
	mov [esp + nb100nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb100nf_half], eax
	mov [esp + nb100nf_half+4], ebx
	movsd xmm1, [esp + nb100nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb100nf_half], xmm1
	movapd [esp + nb100nf_three], xmm3

.nb100nf_threadloop:
        mov   esi, [ebp + nb100nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb100nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb100nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb100nf_n], eax
        mov [esp + nb100nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb100nf_outerstart
        jmp .nb100nf_end

.nb100nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb100nf_nouter]
	mov [esp + nb100nf_nouter], ebx

.nb100nf_outer:
	mov   eax, [ebp + nb100nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb100nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb100nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb100nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb100nf_facel]
	shufpd xmm3, xmm3, 0	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb100nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb100nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb100nf_ix], xmm0
	movapd [esp + nb100nf_iy], xmm1
	movapd [esp + nb100nf_iz], xmm2

	mov   [esp + nb100nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb100nf_vctot], xmm4
	
	mov   eax, [ebp + nb100nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb100nf_pos]
	mov   eax, [ebp + nb100nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb100nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb100nf_ninner]
	mov   [esp + nb100nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb100nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb100nf_unroll_loop
	jmp   .nb100nf_checksingle
.nb100nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb100nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb100nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb100nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]	;# jq A 
	movhpd xmm3, [esi + ebx*8]	;# jq B 

	movapd xmm5, [esp + nb100nf_iq]
	
	mulpd xmm3, xmm5		;# qq 
	
	mov esi, [ebp + nb100nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move nb100nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb100nf_ix]
	movapd xmm5, [esp + nb100nf_iy]
	movapd xmm6, [esp + nb100nf_iz]

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
	movapd xmm1, [esp + nb100nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb100nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	movapd xmm5, [esp + nb100nf_vctot]
	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	addpd  xmm5, xmm3
	movapd [esp + nb100nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb100nf_innerk],  2
	jl    .nb100nf_checksingle
	jmp   .nb100nf_unroll_loop

.nb100nf_checksingle:				
	mov   edx, [esp + nb100nf_innerk]
	and   edx, 1
	jnz    .nb100nf_dosingle
	jmp    .nb100nf_updateouterdata
.nb100nf_dosingle:			
	mov esi, [ebp + nb100nf_charge]
	mov edi, [ebp + nb100nf_pos]

	mov edx, [esp + nb100nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov eax, [edx]	

	xorpd xmm3, xmm3
	movsd xmm3, [esi + eax*8]	;# jq A 
	movapd xmm5, [esp + nb100nf_iq]
	unpcklpd xmm3, xmm6
	mulpd xmm3, xmm5		;# qq 
	
	mov esi, [ebp + nb100nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move nb100nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb100nf_ix]
	movapd xmm5, [esp + nb100nf_iy]
	movapd xmm6, [esp + nb100nf_iz]

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
	movapd xmm1, [esp + nb100nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb100nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb100nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 

	movlpd xmm5, [esp + nb100nf_vctot]
	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	addsd  xmm5, xmm3
	movlpd [esp + nb100nf_vctot], xmm5

.nb100nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb100nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb100nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb100nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb100nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb100nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb100nf_n], esi
        jmp .nb100nf_outer
.nb100nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb100nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb100nf_end
        ;# non-zero, do one more workunit
        jmp   .nb100nf_threadloop
.nb100nf_end:
	emms

	mov eax, [esp + nb100nf_nouter]
	mov ebx, [esp + nb100nf_ninner]
	mov ecx, [ebp + nb100nf_outeriter]
	mov edx, [ebp + nb100nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb100nf_salign]
	add esp, eax
	add esp, 164
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




