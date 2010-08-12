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


.globl nb_kernel200_ia32_sse2
.globl _nb_kernel200_ia32_sse2
nb_kernel200_ia32_sse2:	
_nb_kernel200_ia32_sse2:	
.equiv          nb200_p_nri,            8
.equiv          nb200_iinr,             12
.equiv          nb200_jindex,           16
.equiv          nb200_jjnr,             20
.equiv          nb200_shift,            24
.equiv          nb200_shiftvec,         28
.equiv          nb200_fshift,           32
.equiv          nb200_gid,              36
.equiv          nb200_pos,              40
.equiv          nb200_faction,          44
.equiv          nb200_charge,           48
.equiv          nb200_p_facel,          52
.equiv          nb200_argkrf,           56
.equiv          nb200_argcrf,           60
.equiv          nb200_Vc,               64
.equiv          nb200_type,             68
.equiv          nb200_p_ntype,          72
.equiv          nb200_vdwparam,         76
.equiv          nb200_Vvdw,             80
.equiv          nb200_p_tabscale,       84
.equiv          nb200_VFtab,            88
.equiv          nb200_invsqrta,         92
.equiv          nb200_dvda,             96
.equiv          nb200_p_gbtabscale,     100
.equiv          nb200_GBtab,            104
.equiv          nb200_p_nthreads,       108
.equiv          nb200_count,            112
.equiv          nb200_mtx,              116
.equiv          nb200_outeriter,        120
.equiv          nb200_inneriter,        124
.equiv          nb200_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb200_ix,               0
.equiv          nb200_iy,               16
.equiv          nb200_iz,               32
.equiv          nb200_iq,               48
.equiv          nb200_dx,               64
.equiv          nb200_dy,               80
.equiv          nb200_dz,               96
.equiv          nb200_vctot,            112
.equiv          nb200_fix,              128
.equiv          nb200_fiy,              144
.equiv          nb200_fiz,              160
.equiv          nb200_half,             176
.equiv          nb200_three,            192
.equiv          nb200_two,              208
.equiv          nb200_krf,              224
.equiv          nb200_crf,              240
.equiv          nb200_facel,            256 ;# uses 8 bytes
.equiv          nb200_is3,              264
.equiv          nb200_ii3,              268
.equiv          nb200_innerjjnr,        272
.equiv          nb200_innerk,           276
.equiv          nb200_n,                280
.equiv          nb200_nn1,              284
.equiv          nb200_nri,              288
.equiv          nb200_nouter,           292
.equiv          nb200_ninner,           296
.equiv          nb200_salign,           300
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  304		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb200_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb200_p_nri]
	mov esi, [ebp + nb200_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb200_nri], ecx
	movsd [esp + nb200_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb200_nouter], eax
	mov [esp + nb200_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb200_half], eax
	mov [esp + nb200_half+4], ebx
	movsd xmm1, [esp + nb200_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb200_half], xmm1
	movapd [esp + nb200_two], xmm2
	movapd [esp + nb200_three], xmm3

	mov esi, [ebp + nb200_argkrf]
	mov edi, [ebp + nb200_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb200_krf], xmm5
	shufpd xmm6, xmm6, 0
	movapd [esp + nb200_crf], xmm6

.nb200_threadloop:
        mov   esi, [ebp + nb200_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb200_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb200_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb200_n], eax
        mov [esp + nb200_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200_outerstart
        jmp .nb200_end

.nb200_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb200_nouter]
	mov [esp + nb200_nouter], ebx

.nb200_outer:
	mov   eax, [ebp + nb200_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb200_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb200_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb200_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb200_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb200_facel]
	shufpd xmm3, xmm3, 0
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb200_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb200_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb200_ix], xmm0
	movapd [esp + nb200_iy], xmm1
	movapd [esp + nb200_iz], xmm2

	mov   [esp + nb200_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb200_vctot], xmm4
	movapd [esp + nb200_fix], xmm4
	movapd [esp + nb200_fiy], xmm4
	movapd [esp + nb200_fiz], xmm4
	
	mov   eax, [ebp + nb200_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb200_pos]
	mov   edi, [ebp + nb200_faction]	
	mov   eax, [ebp + nb200_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb200_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb200_ninner]
	mov   [esp + nb200_ninner], ecx
	add   edx, 0
	mov   [esp + nb200_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200_unroll_loop
	jmp   .nb200_checksingle
.nb200_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb200_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb200_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb200_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb200_iq]
	mulpd xmm3, xmm5		;# qq 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mov esi, [ebp + nb200_pos]       ;# base of pos[] 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb200_ix]
	movapd xmm5, [esp + nb200_iy]
	movapd xmm6, [esp + nb200_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb200_dx], xmm4
	movapd [esp + nb200_dy], xmm5
	movapd [esp + nb200_dz], xmm6
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

	movapd xmm7, [esp + nb200_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb200_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb200_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [esp + nb200_crf]
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	
	mulpd  xmm7, [esp + nb200_two]

	subpd  xmm0, xmm7
	mulpd  xmm3, xmm0	
	mulpd  xmm4, xmm3	;# xmm4=total fscal 
	addpd  xmm6, [esp + nb200_vctot]

	movapd xmm0, [esp + nb200_dx]
	movapd xmm1, [esp + nb200_dy]
	movapd xmm2, [esp + nb200_dz]

	movapd [esp + nb200_vctot], xmm6

	mov    edi, [ebp + nb200_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb200_fix]
	movapd xmm4, [esp + nb200_fiy]
	movapd xmm5, [esp + nb200_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb200_fix], xmm3
	movapd [esp + nb200_fiy], xmm4
	movapd [esp + nb200_fiz], xmm5
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
	sub dword ptr [esp + nb200_innerk],  2
	jl    .nb200_checksingle
	jmp   .nb200_unroll_loop

.nb200_checksingle:				
	mov   edx, [esp + nb200_innerk]
	and   edx, 1
	jnz    .nb200_dosingle
	jmp    .nb200_updateouterdata
.nb200_dosingle:			
	mov esi, [ebp + nb200_charge]
	mov edi, [ebp + nb200_pos]
	mov   ecx, [esp + nb200_innerjjnr]
	
	xorpd xmm3, xmm3
	mov   eax, [ecx]

	movlpd xmm3, [esi + eax*8]
	movapd xmm5, [esp + nb200_iq]
	mulpd xmm3, xmm5		;# qq 
	
	mov esi, [ebp + nb200_pos]       ;# base of pos[] 

	lea eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb200_ix]
	movapd xmm5, [esp + nb200_iy]
	movapd xmm6, [esp + nb200_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb200_dx], xmm4
	movapd [esp + nb200_dy], xmm5
	movapd [esp + nb200_dz], xmm6
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

	movapd xmm7, [esp + nb200_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb200_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb200_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [esp + nb200_crf]
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulsd  xmm7, [esp + nb200_two]

	subsd  xmm0, xmm7
	mulsd  xmm3, xmm0
	mulsd  xmm4, xmm3	;# xmm4=total fscal 
	addsd  xmm6, [esp + nb200_vctot]

	movlpd xmm0, [esp + nb200_dx]
	movlpd xmm1, [esp + nb200_dy]
	movlpd xmm2, [esp + nb200_dz]

	movlpd [esp + nb200_vctot], xmm6

	mov    edi, [ebp + nb200_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movlpd xmm3, [esp + nb200_fix]
	movlpd xmm4, [esp + nb200_fiy]
	movlpd xmm5, [esp + nb200_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb200_fix], xmm3
	movlpd [esp + nb200_fiy], xmm4
	movlpd [esp + nb200_fiz], xmm5
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

.nb200_updateouterdata:
	mov   ecx, [esp + nb200_ii3]
	mov   edi, [ebp + nb200_faction]
	mov   esi, [ebp + nb200_fshift]
	mov   edx, [esp + nb200_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb200_fix]
	movapd xmm1, [esp + nb200_fiy]
	movapd xmm2, [esp + nb200_fiz]

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
	mov esi, [esp + nb200_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb200_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb200_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb200_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb200_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb200_n], esi
        jmp .nb200_outer
.nb200_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb200_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200_end
        ;# non-zero, do one more workunit
        jmp   .nb200_threadloop
.nb200_end:
	emms

	mov eax, [esp + nb200_nouter]
	mov ebx, [esp + nb200_ninner]
	mov ecx, [ebp + nb200_outeriter]
	mov edx, [ebp + nb200_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb200_salign]
	add esp, eax
	add esp,  304
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



.globl nb_kernel200nf_ia32_sse2
.globl _nb_kernel200nf_ia32_sse2
nb_kernel200nf_ia32_sse2:	
_nb_kernel200nf_ia32_sse2:	
.equiv          nb200nf_p_nri,          8
.equiv          nb200nf_iinr,           12
.equiv          nb200nf_jindex,         16
.equiv          nb200nf_jjnr,           20
.equiv          nb200nf_shift,          24
.equiv          nb200nf_shiftvec,       28
.equiv          nb200nf_fshift,         32
.equiv          nb200nf_gid,            36
.equiv          nb200nf_pos,            40
.equiv          nb200nf_faction,        44
.equiv          nb200nf_charge,         48
.equiv          nb200nf_p_facel,        52
.equiv          nb200nf_argkrf,         56
.equiv          nb200nf_argcrf,         60
.equiv          nb200nf_Vc,             64
.equiv          nb200nf_type,           68
.equiv          nb200nf_p_ntype,        72
.equiv          nb200nf_vdwparam,       76
.equiv          nb200nf_Vvdw,           80
.equiv          nb200nf_p_tabscale,     84
.equiv          nb200nf_VFtab,          88
.equiv          nb200nf_invsqrta,       92
.equiv          nb200nf_dvda,           96
.equiv          nb200nf_p_gbtabscale,   100
.equiv          nb200nf_GBtab,          104
.equiv          nb200nf_p_nthreads,     108
.equiv          nb200nf_count,          112
.equiv          nb200nf_mtx,            116
.equiv          nb200nf_outeriter,      120
.equiv          nb200nf_inneriter,      124
.equiv          nb200nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb200nf_ix,             0
.equiv          nb200nf_iy,             16
.equiv          nb200nf_iz,             32
.equiv          nb200nf_iq,             48
.equiv          nb200nf_vctot,          64
.equiv          nb200nf_half,           80
.equiv          nb200nf_three,          96
.equiv          nb200nf_krf,            112
.equiv          nb200nf_crf,            128
.equiv          nb200nf_is3,            144
.equiv          nb200nf_ii3,            148
.equiv          nb200nf_innerjjnr,      152
.equiv          nb200nf_innerk,         156
.equiv          nb200nf_n,              160
.equiv          nb200nf_nn1,            164
.equiv          nb200nf_nri,            168
.equiv          nb200nf_facel,          176   ;# uses 8 bytes
.equiv          nb200nf_nouter,         184
.equiv          nb200nf_ninner,         188
.equiv          nb200nf_salign,         192
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  172		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb200nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb200nf_p_nri]
	mov esi, [ebp + nb200nf_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb200nf_nri], ecx
	movsd [esp + nb200nf_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb200nf_nouter], eax
	mov [esp + nb200nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb200nf_half], eax
	mov [esp + nb200nf_half+4], ebx
	movsd xmm1, [esp + nb200nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb200nf_half], xmm1
	movapd [esp + nb200nf_three], xmm3

	mov esi, [ebp + nb200nf_argkrf]
	mov edi, [ebp + nb200nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb200nf_krf], xmm5
	shufpd xmm6, xmm6, 0
	movapd [esp + nb200nf_crf], xmm6

.nb200nf_threadloop:
        mov   esi, [ebp + nb200nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb200nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb200nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb200nf_n], eax
        mov [esp + nb200nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200nf_outerstart
        jmp .nb200nf_end

.nb200nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb200nf_nouter]
	mov [esp + nb200nf_nouter], ebx

.nb200nf_outer:
	mov   eax, [ebp + nb200nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb200nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb200nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb200nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb200nf_facel]
	shufpd xmm3, xmm3, 0
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb200nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb200nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb200nf_ix], xmm0
	movapd [esp + nb200nf_iy], xmm1
	movapd [esp + nb200nf_iz], xmm2

	mov   [esp + nb200nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb200nf_vctot], xmm4
	
	mov   eax, [ebp + nb200nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb200nf_pos]
	mov   eax, [ebp + nb200nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb200nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb200nf_ninner]
	mov   [esp + nb200nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb200nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200nf_unroll_loop
	jmp   .nb200nf_checksingle
.nb200nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb200nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb200nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb200nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb200nf_iq]
	mulpd xmm3, xmm5		;# qq 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mov esi, [ebp + nb200nf_pos]       ;# base of pos[] 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb200nf_ix]
	movapd xmm5, [esp + nb200nf_iy]
	movapd xmm6, [esp + nb200nf_iz]

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

	movapd xmm7, [esp + nb200nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb200nf_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb200nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [esp + nb200nf_crf]
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addpd  xmm6, [esp + nb200nf_vctot]
	movapd [esp + nb200nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [esp + nb200nf_innerk],  2
	jl    .nb200nf_checksingle
	jmp   .nb200nf_unroll_loop

.nb200nf_checksingle:				
	mov   edx, [esp + nb200nf_innerk]
	and   edx, 1
	jnz    .nb200nf_dosingle
	jmp    .nb200nf_updateouterdata
.nb200nf_dosingle:			
	mov esi, [ebp + nb200nf_charge]
	mov edi, [ebp + nb200nf_pos]
	mov   ecx, [esp + nb200nf_innerjjnr]
	
	xorpd xmm3, xmm3
	mov   eax, [ecx]

	movlpd xmm3, [esi + eax*8]
	movapd xmm5, [esp + nb200nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	mov esi, [ebp + nb200nf_pos]       ;# base of pos[] 

	lea eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb200nf_ix]
	movapd xmm5, [esp + nb200nf_iy]
	movapd xmm6, [esp + nb200nf_iz]

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

	movapd xmm7, [esp + nb200nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb200nf_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb200nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb200nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [esp + nb200nf_crf]
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addsd  xmm6, [esp + nb200nf_vctot]
	movlpd [esp + nb200nf_vctot], xmm6

.nb200nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb200nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb200nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb200nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb200nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb200nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb200nf_n], esi
        jmp .nb200nf_outer
.nb200nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb200nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200nf_end
        ;# non-zero, do one more workunit
        jmp   .nb200nf_threadloop
.nb200nf_end:
	emms

	mov eax, [esp + nb200nf_nouter]
	mov ebx, [esp + nb200nf_ninner]
	mov ecx, [ebp + nb200nf_outeriter]
	mov edx, [ebp + nb200nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb200nf_salign]
	add esp, eax
	add esp,  172
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

