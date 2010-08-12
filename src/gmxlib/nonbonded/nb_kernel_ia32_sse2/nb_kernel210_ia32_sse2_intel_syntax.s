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


.globl nb_kernel210_ia32_sse2
.globl _nb_kernel210_ia32_sse2
nb_kernel210_ia32_sse2:	
_nb_kernel210_ia32_sse2:	
.equiv          nb210_p_nri,            8
.equiv          nb210_iinr,             12
.equiv          nb210_jindex,           16
.equiv          nb210_jjnr,             20
.equiv          nb210_shift,            24
.equiv          nb210_shiftvec,         28
.equiv          nb210_fshift,           32
.equiv          nb210_gid,              36
.equiv          nb210_pos,              40
.equiv          nb210_faction,          44
.equiv          nb210_charge,           48
.equiv          nb210_p_facel,          52
.equiv          nb210_argkrf,           56
.equiv          nb210_argcrf,           60
.equiv          nb210_Vc,               64
.equiv          nb210_type,             68
.equiv          nb210_p_ntype,          72
.equiv          nb210_vdwparam,         76
.equiv          nb210_Vvdw,             80
.equiv          nb210_p_tabscale,       84
.equiv          nb210_VFtab,            88
.equiv          nb210_invsqrta,         92
.equiv          nb210_dvda,             96
.equiv          nb210_p_gbtabscale,     100
.equiv          nb210_GBtab,            104
.equiv          nb210_p_nthreads,       108
.equiv          nb210_count,            112
.equiv          nb210_mtx,              116
.equiv          nb210_outeriter,        120
.equiv          nb210_inneriter,        124
.equiv          nb210_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb210_ix,               0
.equiv          nb210_iy,               16
.equiv          nb210_iz,               32
.equiv          nb210_iq,               48
.equiv          nb210_dx,               64
.equiv          nb210_dy,               80
.equiv          nb210_dz,               96
.equiv          nb210_c6,               112
.equiv          nb210_c12,              128
.equiv          nb210_six,              144
.equiv          nb210_twelve,           160
.equiv          nb210_vctot,            176
.equiv          nb210_Vvdwtot,          192
.equiv          nb210_fix,              208
.equiv          nb210_fiy,              224
.equiv          nb210_fiz,              240
.equiv          nb210_half,             256
.equiv          nb210_three,            272
.equiv          nb210_two,              288
.equiv          nb210_krf,              304
.equiv          nb210_crf,              320
.equiv          nb210_facel,            336 ;# uses 8 bytes
.equiv          nb210_is3,              344
.equiv          nb210_ii3,              348
.equiv          nb210_ntia,             352
.equiv          nb210_innerjjnr,        356
.equiv          nb210_innerk,           360
.equiv          nb210_n,                364
.equiv          nb210_nn1,              368
.equiv          nb210_nri,              372
.equiv          nb210_ntype,            376
.equiv          nb210_nouter,           380
.equiv          nb210_ninner,           384
.equiv          nb210_salign,           388
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
	mov [esp + nb210_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb210_p_nri]
	mov esi, [ebp + nb210_p_facel]
	mov edi, [ebp + nb210_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb210_nri], ecx
	movsd [esp + nb210_facel], xmm7
	mov [esp + nb210_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb210_nouter], eax
	mov [esp + nb210_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb210_half], eax
	mov [esp + nb210_half+4], ebx
	movsd xmm1, [esp + nb210_half]
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
	movapd [esp + nb210_half], xmm1
	movapd [esp + nb210_two], xmm2
	movapd [esp + nb210_three], xmm3
	movapd [esp + nb210_six], xmm4
	movapd [esp + nb210_twelve], xmm5

	mov esi, [ebp + nb210_argkrf]
	mov edi, [ebp + nb210_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb210_krf], xmm5
	movapd [esp + nb210_crf], xmm6

.nb210_threadloop:
        mov   esi, [ebp + nb210_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb210_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb210_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb210_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb210_n], eax
        mov [esp + nb210_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb210_outerstart
        jmp .nb210_end

.nb210_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb210_nouter]
	mov [esp + nb210_nouter], ebx

.nb210_outer:
	mov   eax, [ebp + nb210_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb210_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb210_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb210_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb210_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb210_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb210_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb210_ntype]
    	shl   edx, 1
    	mov   [esp + nb210_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb210_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb210_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb210_ix], xmm0
	movapd [esp + nb210_iy], xmm1
	movapd [esp + nb210_iz], xmm2

	mov   [esp + nb210_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb210_vctot], xmm4
	movapd [esp + nb210_Vvdwtot], xmm4
	movapd [esp + nb210_fix], xmm4
	movapd [esp + nb210_fiy], xmm4
	movapd [esp + nb210_fiz], xmm4
	
	mov   eax, [ebp + nb210_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb210_pos]
	mov   edi, [ebp + nb210_faction]	
	mov   eax, [ebp + nb210_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb210_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb210_ninner]
	mov   [esp + nb210_ninner], ecx
	add   edx, 0
	mov   [esp + nb210_innerk], edx    ;# number of innerloop atoms 
	jge   .nb210_unroll_loop
	jmp   .nb210_checksingle
.nb210_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb210_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb210_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb210_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb210_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb210_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb210_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb210_ntia]
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
	movapd [esp + nb210_c6], xmm4
	movapd [esp + nb210_c12], xmm6
	
	mov esi, [ebp + nb210_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb210_ix]
	movapd xmm5, [esp + nb210_iy]
	movapd xmm6, [esp + nb210_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb210_dx], xmm4
	movapd [esp + nb210_dy], xmm5
	movapd [esp + nb210_dz], xmm6
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

	movapd xmm7, [esp + nb210_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb210_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb210_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [esp + nb210_crf]
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulpd  xmm7, [esp + nb210_two]
	mulpd  xmm1, [esp + nb210_c6]
	mulpd  xmm2, [esp + nb210_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [esp + nb210_Vvdwtot]
	mulpd  xmm1, [esp + nb210_six]
	mulpd  xmm2, [esp + nb210_twelve]
	subpd  xmm2, xmm1
	subpd  xmm0, xmm7
	mulpd  xmm3, xmm0
	addpd  xmm2, xmm3
	mulpd  xmm4, xmm2	;# xmm4=total fscal 
	addpd  xmm6, [esp + nb210_vctot]

	movapd xmm0, [esp + nb210_dx]
	movapd xmm1, [esp + nb210_dy]
	movapd xmm2, [esp + nb210_dz]

	movapd [esp + nb210_vctot], xmm6
	movapd [esp + nb210_Vvdwtot], xmm5

	mov    edi, [ebp + nb210_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb210_fix]
	movapd xmm4, [esp + nb210_fiy]
	movapd xmm5, [esp + nb210_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb210_fix], xmm3
	movapd [esp + nb210_fiy], xmm4
	movapd [esp + nb210_fiz], xmm5
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
	sub dword ptr [esp + nb210_innerk],  2
	jl    .nb210_checksingle
	jmp   .nb210_unroll_loop

.nb210_checksingle:				
	mov   edx, [esp + nb210_innerk]
	and   edx, 1
	jnz    .nb210_dosingle
	jmp    .nb210_updateouterdata
.nb210_dosingle:			
	mov esi, [ebp + nb210_charge]
	mov edi, [ebp + nb210_pos]
	mov   ecx, [esp + nb210_innerjjnr]
	xorpd xmm3, xmm3
	mov   eax, [ecx]

	movlpd xmm3, [esi + eax*8]
	movapd xmm5, [esp + nb210_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb210_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb210_vdwparam]
	shl eax, 1
	mov edi, [esp + nb210_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		
	movapd [esp + nb210_c6], xmm4
	movapd [esp + nb210_c12], xmm6
	
	mov esi, [ebp + nb210_pos]       ;# base of pos[] 

	lea eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb210_ix]
	movapd xmm5, [esp + nb210_iy]
	movapd xmm6, [esp + nb210_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb210_dx], xmm4
	movapd [esp + nb210_dy], xmm5
	movapd [esp + nb210_dz], xmm6
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

	movapd xmm7, [esp + nb210_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb210_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb210_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [esp + nb210_crf]
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulsd  xmm7, [esp + nb210_two]
	mulsd  xmm1, [esp + nb210_c6]
	mulsd  xmm2, [esp + nb210_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [esp + nb210_Vvdwtot]
	mulsd  xmm1, [esp + nb210_six]
	mulsd  xmm2, [esp + nb210_twelve]
	subsd  xmm2, xmm1
	subsd  xmm0, xmm7
	mulsd  xmm3, xmm0
	addsd  xmm2, xmm3
	mulsd  xmm4, xmm2	;# xmm4=total fscal 
	addsd  xmm6, [esp + nb210_vctot]

	movlpd xmm0, [esp + nb210_dx]
	movlpd xmm1, [esp + nb210_dy]
	movlpd xmm2, [esp + nb210_dz]

	movlpd [esp + nb210_vctot], xmm6
	movlpd [esp + nb210_Vvdwtot], xmm5

	mov    edi, [ebp + nb210_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movlpd xmm3, [esp + nb210_fix]
	movlpd xmm4, [esp + nb210_fiy]
	movlpd xmm5, [esp + nb210_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb210_fix], xmm3
	movlpd [esp + nb210_fiy], xmm4
	movlpd [esp + nb210_fiz], xmm5
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
	
.nb210_updateouterdata:
	mov   ecx, [esp + nb210_ii3]
	mov   edi, [ebp + nb210_faction]
	mov   esi, [ebp + nb210_fshift]
	mov   edx, [esp + nb210_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb210_fix]
	movapd xmm1, [esp + nb210_fiy]
	movapd xmm2, [esp + nb210_fiz]

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
	mov esi, [esp + nb210_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb210_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb210_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb210_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb210_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb210_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb210_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb210_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb210_n], esi
        jmp .nb210_outer
.nb210_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb210_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb210_end
        ;# non-zero, do one more workunit
        jmp   .nb210_threadloop
.nb210_end:
	emms

	mov eax, [esp + nb210_nouter]
	mov ebx, [esp + nb210_ninner]
	mov ecx, [ebp + nb210_outeriter]
	mov edx, [ebp + nb210_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb210_salign]
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


.globl nb_kernel210nf_ia32_sse2
.globl _nb_kernel210nf_ia32_sse2
nb_kernel210nf_ia32_sse2:	
_nb_kernel210nf_ia32_sse2:	
.equiv          nb210nf_p_nri,          8
.equiv          nb210nf_iinr,           12
.equiv          nb210nf_jindex,         16
.equiv          nb210nf_jjnr,           20
.equiv          nb210nf_shift,          24
.equiv          nb210nf_shiftvec,       28
.equiv          nb210nf_fshift,         32
.equiv          nb210nf_gid,            36
.equiv          nb210nf_pos,            40
.equiv          nb210nf_faction,        44
.equiv          nb210nf_charge,         48
.equiv          nb210nf_p_facel,        52
.equiv          nb210nf_argkrf,         56
.equiv          nb210nf_argcrf,         60
.equiv          nb210nf_Vc,             64
.equiv          nb210nf_type,           68
.equiv          nb210nf_p_ntype,        72
.equiv          nb210nf_vdwparam,       76
.equiv          nb210nf_Vvdw,           80
.equiv          nb210nf_p_tabscale,     84
.equiv          nb210nf_VFtab,          88
.equiv          nb210nf_invsqrta,       92
.equiv          nb210nf_dvda,           96
.equiv          nb210nf_p_gbtabscale,   100
.equiv          nb210nf_GBtab,          104
.equiv          nb210nf_p_nthreads,     108
.equiv          nb210nf_count,          112
.equiv          nb210nf_mtx,            116
.equiv          nb210nf_outeriter,      120
.equiv          nb210nf_inneriter,      124
.equiv          nb210nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb210nf_ix,             0
.equiv          nb210nf_iy,             16
.equiv          nb210nf_iz,             32
.equiv          nb210nf_iq,             48
.equiv          nb210nf_c6,             64
.equiv          nb210nf_c12,            80
.equiv          nb210nf_vctot,          96
.equiv          nb210nf_Vvdwtot,        112
.equiv          nb210nf_half,           128
.equiv          nb210nf_three,          144
.equiv          nb210nf_krf,            160
.equiv          nb210nf_crf,            176
.equiv          nb210nf_is3,            192
.equiv          nb210nf_ii3,            196
.equiv          nb210nf_ntia,           200
.equiv          nb210nf_innerjjnr,      204
.equiv          nb210nf_innerk,         208
.equiv          nb210nf_n,              212
.equiv          nb210nf_nn1,            216
.equiv          nb210nf_nri,            220
.equiv          nb210nf_facel,          224   ;# uses 8 bytes
.equiv          nb210nf_ntype,          232
.equiv          nb210nf_nouter,         236
.equiv          nb210nf_ninner,         240
.equiv          nb210nf_salign,         244
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp,  224		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb210nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb210nf_p_nri]
	mov esi, [ebp + nb210nf_p_facel]
	mov edi, [ebp + nb210nf_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb210nf_nri], ecx
	movsd [esp + nb210nf_facel], xmm7
	mov [esp + nb210nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb210nf_nouter], eax
	mov [esp + nb210nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb210nf_half], eax
	mov [esp + nb210nf_half+4], ebx
	movsd xmm1, [esp + nb210nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb210nf_half], xmm1
	movapd [esp + nb210nf_three], xmm3

	mov esi, [ebp + nb210nf_argkrf]
	mov edi, [ebp + nb210nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb210nf_krf], xmm5
	movapd [esp + nb210nf_crf], xmm6

.nb210nf_threadloop:
        mov   esi, [ebp + nb210nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb210nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb210nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb210nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb210nf_n], eax
        mov [esp + nb210nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb210nf_outerstart
        jmp .nb210nf_end

.nb210nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb210nf_nouter]
	mov [esp + nb210nf_nouter], ebx

.nb210nf_outer:
	mov   eax, [ebp + nb210nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb210nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb210nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb210nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb210nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb210nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb210nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb210nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb210nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb210nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb210nf_ix], xmm0
	movapd [esp + nb210nf_iy], xmm1
	movapd [esp + nb210nf_iz], xmm2

	mov   [esp + nb210nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb210nf_vctot], xmm4
	movapd [esp + nb210nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb210nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb210nf_pos]
	mov   eax, [ebp + nb210nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb210nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb210nf_ninner]
	mov   [esp + nb210nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb210nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb210nf_unroll_loop
	jmp   .nb210nf_checksingle
.nb210nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb210nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb210nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb210nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm5, [esp + nb210nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb210nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb210nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb210nf_ntia]
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
	movapd [esp + nb210nf_c6], xmm4
	movapd [esp + nb210nf_c12], xmm6
	
	mov esi, [ebp + nb210nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb210nf_ix]
	movapd xmm5, [esp + nb210nf_iy]
	movapd xmm6, [esp + nb210nf_iz]

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

	movapd xmm7, [esp + nb210nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb210nf_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb210nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [esp + nb210nf_crf]
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulpd  xmm1, [esp + nb210nf_c6]
	mulpd  xmm2, [esp + nb210nf_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [esp + nb210nf_Vvdwtot]
	addpd  xmm6, [esp + nb210nf_vctot]
	movapd [esp + nb210nf_vctot], xmm6
	movapd [esp + nb210nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb210nf_innerk],  2
	jl    .nb210nf_checksingle
	jmp   .nb210nf_unroll_loop

.nb210nf_checksingle:				
	mov   edx, [esp + nb210nf_innerk]
	and   edx, 1
	jnz    .nb210nf_dosingle
	jmp    .nb210nf_updateouterdata
.nb210nf_dosingle:			
	mov esi, [ebp + nb210nf_charge]
	mov edi, [ebp + nb210nf_pos]
	mov   ecx, [esp + nb210nf_innerjjnr]
	xorpd xmm3, xmm3
	mov   eax, [ecx]

	movlpd xmm3, [esi + eax*8]
	movapd xmm5, [esp + nb210nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb210nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb210nf_vdwparam]
	shl eax, 1
	mov edi, [esp + nb210nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		
	movapd [esp + nb210nf_c6], xmm4
	movapd [esp + nb210nf_c12], xmm6
	
	mov esi, [ebp + nb210nf_pos]       ;# base of pos[] 

	lea eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb210nf_ix]
	movapd xmm5, [esp + nb210nf_iy]
	movapd xmm6, [esp + nb210nf_iz]

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

	movapd xmm7, [esp + nb210nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb210nf_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb210nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb210nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [esp + nb210nf_crf]
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulsd  xmm1, [esp + nb210nf_c6]
	mulsd  xmm2, [esp + nb210nf_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [esp + nb210nf_Vvdwtot]
	addsd  xmm6, [esp + nb210nf_vctot]
	movlpd [esp + nb210nf_vctot], xmm6
	movlpd [esp + nb210nf_Vvdwtot], xmm5
	
.nb210nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb210nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb210nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb210nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb210nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb210nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb210nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb210nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb210nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb210nf_n], esi
        jmp .nb210nf_outer
.nb210nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb210nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb210nf_end
        ;# non-zero, do one more workunit
        jmp   .nb210nf_threadloop
.nb210nf_end:
	emms

	mov eax, [esp + nb210nf_nouter]
	mov ebx, [esp + nb210nf_ninner]
	mov ecx, [ebp + nb210nf_outeriter]
	mov edx, [ebp + nb210nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb210nf_salign]
	add esp, eax
	add esp,  224
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


