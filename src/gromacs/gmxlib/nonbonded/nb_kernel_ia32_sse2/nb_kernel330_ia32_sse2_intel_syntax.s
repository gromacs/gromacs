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

.globl nb_kernel330_ia32_sse2
.globl _nb_kernel330_ia32_sse2
nb_kernel330_ia32_sse2:	
_nb_kernel330_ia32_sse2:	
.equiv          nb330_p_nri,            8
.equiv          nb330_iinr,             12
.equiv          nb330_jindex,           16
.equiv          nb330_jjnr,             20
.equiv          nb330_shift,            24
.equiv          nb330_shiftvec,         28
.equiv          nb330_fshift,           32
.equiv          nb330_gid,              36
.equiv          nb330_pos,              40
.equiv          nb330_faction,          44
.equiv          nb330_charge,           48
.equiv          nb330_p_facel,          52
.equiv          nb330_argkrf,           56
.equiv          nb330_argcrf,           60
.equiv          nb330_Vc,               64
.equiv          nb330_type,             68
.equiv          nb330_p_ntype,          72
.equiv          nb330_vdwparam,         76
.equiv          nb330_Vvdw,             80
.equiv          nb330_p_tabscale,       84
.equiv          nb330_VFtab,            88
.equiv          nb330_invsqrta,         92
.equiv          nb330_dvda,             96
.equiv          nb330_p_gbtabscale,     100
.equiv          nb330_GBtab,            104
.equiv          nb330_p_nthreads,       108
.equiv          nb330_count,            112
.equiv          nb330_mtx,              116
.equiv          nb330_outeriter,        120
.equiv          nb330_inneriter,        124
.equiv          nb330_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb330_ix,               0
.equiv          nb330_iy,               16
.equiv          nb330_iz,               32
.equiv          nb330_iq,               48
.equiv          nb330_dx,               64
.equiv          nb330_dy,               80
.equiv          nb330_dz,               96
.equiv          nb330_two,              112
.equiv          nb330_tsc,              128
.equiv          nb330_qq,               144
.equiv          nb330_c6,               160
.equiv          nb330_c12,              176
.equiv          nb330_fscal,            192
.equiv          nb330_vctot,            208
.equiv          nb330_Vvdwtot,          224
.equiv          nb330_fix,              240
.equiv          nb330_fiy,              256
.equiv          nb330_fiz,              272
.equiv          nb330_half,             288
.equiv          nb330_three,            304
.equiv          nb330_is3,              320
.equiv          nb330_ii3,              324
.equiv          nb330_ntia,             328
.equiv          nb330_innerjjnr,        332
.equiv          nb330_innerk,           336
.equiv          nb330_n,                340
.equiv          nb330_nn1,              344
.equiv          nb330_nri,              348
.equiv          nb330_facel,            352   ;# uses 8 bytes
.equiv          nb330_ntype,            360
.equiv          nb330_nouter,           364
.equiv          nb330_ninner,           368
.equiv          nb330_salign,           372
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
	mov [esp + nb330_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb330_p_nri]
	mov esi, [ebp + nb330_p_facel]
	mov edi, [ebp + nb330_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb330_nri], ecx
	movsd [esp + nb330_facel], xmm7
	mov [esp + nb330_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb330_nouter], eax
	mov [esp + nb330_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb330_half], eax
	mov [esp + nb330_half+4], ebx
	movsd xmm1, [esp + nb330_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb330_half], xmm1
	movapd [esp + nb330_two], xmm2
	movapd [esp + nb330_three], xmm3
	mov eax, [ebp + nb330_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb330_tsc], xmm3

.nb330_threadloop:
        mov   esi, [ebp + nb330_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb330_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb330_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb330_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb330_n], eax
        mov [esp + nb330_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb330_outerstart
        jmp .nb330_end

.nb330_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb330_nouter]
	mov [esp + nb330_nouter], ebx

.nb330_outer:
	mov   eax, [ebp + nb330_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb330_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb330_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb330_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb330_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb330_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb330_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb330_ntype]
    	shl   edx, 1
    	mov   [esp + nb330_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb330_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb330_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb330_ix], xmm0
	movapd [esp + nb330_iy], xmm1
	movapd [esp + nb330_iz], xmm2

	mov   [esp + nb330_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb330_vctot], xmm4
	movapd [esp + nb330_Vvdwtot], xmm4
	movapd [esp + nb330_fix], xmm4
	movapd [esp + nb330_fiy], xmm4
	movapd [esp + nb330_fiz], xmm4
	
	mov   eax, [ebp + nb330_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb330_pos]
	mov   edi, [ebp + nb330_faction]	
	mov   eax, [ebp + nb330_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb330_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb330_ninner]
	mov   [esp + nb330_ninner], ecx
	add   edx, 0
	mov   [esp + nb330_innerk], edx    ;# number of innerloop atoms 
	jge   .nb330_unroll_loop
	jmp   .nb330_checksingle
.nb330_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb330_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb330_innerjjnr], 8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb330_charge]    ;# base of charge[] 

	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm2, [esp + nb330_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb330_qq], xmm3	

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb330_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb330_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb330_ntia]
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
	movapd [esp + nb330_c6], xmm4
	movapd [esp + nb330_c12], xmm6
	
	mov esi, [ebp + nb330_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	mov    edi, [ebp + nb330_faction]
	
	;# move nb330_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb330_ix]
	movapd xmm5, [esp + nb330_iy]
	movapd xmm6, [esp + nb330_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb330_dx], xmm4
	movapd [esp + nb330_dy], xmm5
	movapd [esp + nb330_dz], xmm6
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
	movapd xmm1, [esp + nb330_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb330_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb330_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 

	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb330_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	;# Coulomb 
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
	mulpd  xmm7, [esp + nb330_two]	;# two*Heps2 
	movapd xmm3, [esp + nb330_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [esp + nb330_vctot]
	movapd [esp + nb330_vctot], xmm5 

	;# put scalar force on stack temporarily 
	movapd [esp + nb330_fscal], xmm3

	;# Dispersion 
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
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb330_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb330_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 
	addpd  xmm7, [esp + nb330_fscal] ;# add to fscal 

	;# put scalar force back on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb330_Vvdwtot]
	movapd [esp + nb330_fscal], xmm7
	movapd [esp + nb330_Vvdwtot], xmm5

	;# Repulsion 
	movlpd xmm4, [esi + eax*8 + 64]	;# Y1	
	movlpd xmm3, [esi + ebx*8 + 64]	;# Y2
	movhpd xmm4, [esi + eax*8 + 72]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 72]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 80]	;# G1
	movlpd xmm3, [esi + ebx*8 + 80]	;# G2
	movhpd xmm6, [esi + eax*8 + 88]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 88]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb330_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb330_c12]
	mulpd  xmm7, xmm4 ;# fijR 
	mulpd  xmm5, xmm4 ;# Vvdw12 
	addpd  xmm7, [esp + nb330_fscal] 
	
	addpd  xmm5, [esp + nb330_Vvdwtot]
	movapd [esp + nb330_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulpd xmm7, [esp + nb330_tsc]
	mulpd xmm7, xmm0
	subpd xmm4, xmm7

	movapd xmm0, [esp + nb330_dx]
	movapd xmm1, [esp + nb330_dy]
	movapd xmm2, [esp + nb330_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb330_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb330_fix]
	movapd xmm4, [esp + nb330_fiy]
	movapd xmm5, [esp + nb330_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb330_fix], xmm3
	movapd [esp + nb330_fiy], xmm4
	movapd [esp + nb330_fiz], xmm5
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
	sub dword ptr [esp + nb330_innerk],  2
	jl    .nb330_checksingle
	jmp   .nb330_unroll_loop
.nb330_checksingle:
	mov   edx, [esp + nb330_innerk]
	and   edx, 1
	jnz    .nb330_dosingle
	jmp    .nb330_updateouterdata
.nb330_dosingle:
	mov esi, [ebp + nb330_charge]
	mov edi, [ebp + nb330_pos]
	mov   ecx, [esp + nb330_innerjjnr]
	mov   eax, [ecx]	
	xorpd  xmm3, xmm3
	movlpd xmm3, [esi + eax*8]	;# xmm6(0) has the charge 	
	mulpd  xmm3, [esp + nb330_iq]
	movapd [esp + nb330_qq], xmm3

	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov esi, [ebp + nb330_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb330_vdwparam]
	shl eax, 1
	mov edi, [esp + nb330_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb330_c6], xmm4
	movapd [esp + nb330_c12], xmm6
	
	mov esi, [ebp + nb330_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	mov    edi, [ebp + nb330_faction]

	;# move nb330_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb330_ix]
	movapd xmm5, [esp + nb330_iy]
	movapd xmm6, [esp + nb330_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb330_dx], xmm4
	movapd [esp + nb330_dy], xmm5
	movapd [esp + nb330_dz], xmm6
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
	movapd xmm1, [esp + nb330_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb330_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb330_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb330_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	;# Coulomb 
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
	mulsd  xmm7, [esp + nb330_two]	;# two*Heps2 
	movapd xmm3, [esp + nb330_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [esp + nb330_vctot]
	movlpd [esp + nb330_vctot], xmm5 

	;# put scalar force on stack temporarily 
	movapd [esp + nb330_fscal], xmm3

	;# Dispersion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1	
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb330_two]	;# two*Heps2 
	movapd xmm3, [esp + nb330_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb330_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 
	addsd  xmm7, [esp + nb330_fscal] ;# add to fscal 

	;# put scalar force back on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb330_Vvdwtot]
	movlpd [esp + nb330_fscal], xmm7
	movlpd [esp + nb330_Vvdwtot], xmm5

	;# Repulsion 
	movlpd xmm4, [esi + eax*8 + 64]	;# Y1 	
	movhpd xmm4, [esi + eax*8 + 72]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movlpd xmm6, [esi + eax*8 + 80]	;# G1	
	movhpd xmm6, [esi + eax*8 + 88]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb330_two]	;# two*Heps2 
	movapd xmm3, [esp + nb330_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb330_c12]
	mulsd  xmm7, xmm4 ;# fijR 
	mulsd  xmm5, xmm4 ;# Vvdw12 
	addsd  xmm7, [esp + nb330_fscal] 
	
	addsd  xmm5, [esp + nb330_Vvdwtot]
	movlpd [esp + nb330_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulsd xmm7, [esp + nb330_tsc]
	mulsd xmm7, xmm0
	subsd xmm4, xmm7

	movapd xmm0, [esp + nb330_dx]
	movapd xmm1, [esp + nb330_dy]
	movapd xmm2, [esp + nb330_dz]

	movd eax, mm0	

	mov    edi, [ebp + nb330_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb330_fix]
	movapd xmm4, [esp + nb330_fiy]
	movapd xmm5, [esp + nb330_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb330_fix], xmm3
	movlpd [esp + nb330_fiy], xmm4
	movlpd [esp + nb330_fiz], xmm5
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
.nb330_updateouterdata:
	mov   ecx, [esp + nb330_ii3]
	mov   edi, [ebp + nb330_faction]
	mov   esi, [ebp + nb330_fshift]
	mov   edx, [esp + nb330_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb330_fix]
	movapd xmm1, [esp + nb330_fiy]
	movapd xmm2, [esp + nb330_fiz]

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
	mov esi, [esp + nb330_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb330_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb330_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb330_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb330_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb330_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb330_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb330_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb330_n], esi
        jmp .nb330_outer
.nb330_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb330_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb330_end
        ;# non-zero, do one more workunit
        jmp   .nb330_threadloop
.nb330_end:
	emms

	mov eax, [esp + nb330_nouter]
	mov ebx, [esp + nb330_ninner]
	mov ecx, [ebp + nb330_outeriter]
	mov edx, [ebp + nb330_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb330_salign]
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



.globl nb_kernel330nf_ia32_sse2
.globl _nb_kernel330nf_ia32_sse2
nb_kernel330nf_ia32_sse2:	
_nb_kernel330nf_ia32_sse2:	
.equiv          nb330nf_p_nri,          8
.equiv          nb330nf_iinr,           12
.equiv          nb330nf_jindex,         16
.equiv          nb330nf_jjnr,           20
.equiv          nb330nf_shift,          24
.equiv          nb330nf_shiftvec,       28
.equiv          nb330nf_fshift,         32
.equiv          nb330nf_gid,            36
.equiv          nb330nf_pos,            40
.equiv          nb330nf_faction,        44
.equiv          nb330nf_charge,         48
.equiv          nb330nf_p_facel,        52
.equiv          nb330nf_argkrf,         56
.equiv          nb330nf_argcrf,         60
.equiv          nb330nf_Vc,             64
.equiv          nb330nf_type,           68
.equiv          nb330nf_p_ntype,        72
.equiv          nb330nf_vdwparam,       76
.equiv          nb330nf_Vvdw,           80
.equiv          nb330nf_p_tabscale,     84
.equiv          nb330nf_VFtab,          88
.equiv          nb330nf_invsqrta,       92
.equiv          nb330nf_dvda,           96
.equiv          nb330nf_p_gbtabscale,   100
.equiv          nb330nf_GBtab,          104
.equiv          nb330nf_p_nthreads,     108
.equiv          nb330nf_count,          112
.equiv          nb330nf_mtx,            116
.equiv          nb330nf_outeriter,      120
.equiv          nb330nf_inneriter,      124
.equiv          nb330nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb330nf_ix,             0
.equiv          nb330nf_iy,             16
.equiv          nb330nf_iz,             32
.equiv          nb330nf_iq,             48
.equiv          nb330nf_tsc,            64
.equiv          nb330nf_qq,             80
.equiv          nb330nf_c6,             96
.equiv          nb330nf_c12,            112
.equiv          nb330nf_vctot,          128
.equiv          nb330nf_Vvdwtot,        144
.equiv          nb330nf_half,           160
.equiv          nb330nf_three,          176
.equiv          nb330nf_is3,            192
.equiv          nb330nf_ii3,            196
.equiv          nb330nf_ntia,           200
.equiv          nb330nf_innerjjnr,      204
.equiv          nb330nf_innerk,         208
.equiv          nb330nf_n,              212
.equiv          nb330nf_nn1,            216
.equiv          nb330nf_nri,            220
.equiv          nb330nf_facel,          224   ;# uses 8 bytes
.equiv          nb330nf_ntype,          232
.equiv          nb330nf_nouter,         236
.equiv          nb330nf_ninner,         240
.equiv          nb330nf_salign,         244
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
	mov [esp + nb330nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb330nf_p_nri]
	mov esi, [ebp + nb330nf_p_facel]
	mov edi, [ebp + nb330nf_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb330nf_nri], ecx
	movsd [esp + nb330nf_facel], xmm7
	mov [esp + nb330nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb330nf_nouter], eax
	mov [esp + nb330nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb330nf_half], eax
	mov [esp + nb330nf_half+4], ebx
	movsd xmm1, [esp + nb330nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb330nf_half], xmm1
	movapd [esp + nb330nf_three], xmm3
	mov eax, [ebp + nb330nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb330nf_tsc], xmm3

.nb330nf_threadloop:
        mov   esi, [ebp + nb330nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb330nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb330nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb330nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb330nf_n], eax
        mov [esp + nb330nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb330nf_outerstart
        jmp .nb330nf_end

.nb330nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb330nf_nouter]
	mov [esp + nb330nf_nouter], ebx

.nb330nf_outer:
	mov   eax, [ebp + nb330nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb330nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb330nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb330nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb330nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   edx, [ebp + nb330nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb330nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb330nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb330nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb330nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb330nf_ix], xmm0
	movapd [esp + nb330nf_iy], xmm1
	movapd [esp + nb330nf_iz], xmm2

	mov   [esp + nb330nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb330nf_vctot], xmm4
	movapd [esp + nb330nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb330nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb330nf_pos]
	mov   eax, [ebp + nb330nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb330nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb330nf_ninner]
	mov   [esp + nb330nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb330nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb330nf_unroll_loop
	jmp   .nb330nf_checksingle
.nb330nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb330nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb330nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb330nf_charge]    ;# base of charge[] 

	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm2, [esp + nb330nf_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb330nf_qq], xmm3	

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb330nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb330nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb330nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [esi + eax*8]     ;# c6a 
	movlpd xmm7, [esi + ebx*8]     ; # c6b
	movhpd xmm6, [esi + eax*8 + 8] ; # c6a c12a
	movhpd xmm7, [esi + ebx*8 + 8] ; # c6b c12b
	
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb330nf_c6], xmm4
	movapd [esp + nb330nf_c12], xmm6
	
	mov esi, [ebp + nb330nf_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move nb330nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb330nf_ix]
	movapd xmm5, [esp + nb330nf_iy]
	movapd xmm6, [esp + nb330nf_iz]

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
	movapd xmm1, [esp + nb330nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb330nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb330nf_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 

	mov  esi, [ebp + nb330nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 
	lea   ebx, [ebx + ebx*2]	;# idx*=3 (12 total now) 

	;# Coulomb 
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
	movapd xmm3, [esp + nb330nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV 
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [esp + nb330nf_vctot]
	movapd [esp + nb330nf_vctot], xmm5 

	;# Dispersion 
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
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [esp + nb330nf_c6] ;# Vvdw6 

	addpd  xmm5, [esp + nb330nf_Vvdwtot]
	movapd [esp + nb330nf_Vvdwtot], xmm5

	;# Repulsion 
	movlpd xmm4, [esi + eax*8 + 64]	;# Y1	
	movlpd xmm3, [esi + ebx*8 + 64]	;# Y2
	movhpd xmm4, [esi + eax*8 + 72]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 72]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 80]	;# G1
	movlpd xmm3, [esi + ebx*8 + 80]	;# G2
	movhpd xmm6, [esi + eax*8 + 88]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 88]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [esp + nb330nf_c12] ;# Vvdw12 
	
	addpd  xmm5, [esp + nb330nf_Vvdwtot]
	movapd [esp + nb330nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb330nf_innerk],  2
	jl    .nb330nf_checksingle
	jmp   .nb330nf_unroll_loop
.nb330nf_checksingle:
	mov   edx, [esp + nb330nf_innerk]
	and   edx, 1
	jnz    .nb330nf_dosingle
	jmp    .nb330nf_updateouterdata
.nb330nf_dosingle:
	mov esi, [ebp + nb330nf_charge]
	mov edi, [ebp + nb330nf_pos]
	mov   ecx, [esp + nb330nf_innerjjnr]
	mov   eax, [ecx]	
	xorpd  xmm3, xmm3
	movlpd xmm3, [esi + eax*8]	;# xmm6(0) has the charge 	
	mulpd  xmm3, [esp + nb330nf_iq]
	movapd [esp + nb330nf_qq], xmm3

	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov esi, [ebp + nb330nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb330nf_vdwparam]
	shl eax, 1
	mov edi, [esp + nb330nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb330nf_c6], xmm4
	movapd [esp + nb330nf_c12], xmm6
	
	mov esi, [ebp + nb330nf_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move nb330nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb330nf_ix]
	movapd xmm5, [esp + nb330nf_iy]
	movapd xmm6, [esp + nb330nf_iz]

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
	movapd xmm1, [esp + nb330nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb330nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb330nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb330nf_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb330nf_VFtab]
	lea   eax, [eax + eax*2]	;# idx*=3 (12 total now) 

	;# Coulomb 
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
	movapd xmm3, [esp + nb330nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [esp + nb330nf_vctot]
	movlpd [esp + nb330nf_vctot], xmm5 

	;# Dispersion 
	movsd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 40]	;# Y1 F1 	
	movsd xmm6, [esi + eax*8 + 48]	;# G1	
	movsd xmm7, [esi + eax*8 + 56]	;# G1 H1 	
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb330nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [esp + nb330nf_c6]	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [esp + nb330nf_Vvdwtot]
	movlpd [esp + nb330nf_Vvdwtot], xmm5

	;# Repulsion 
	movsd xmm4, [esi + eax*8 + 64]	;# Y1 	
	movsd xmm5, [esi + eax*8 + 72]	;# Y1 F1 	
	movsd xmm6, [esi + eax*8 + 80]	;# G1	
	movsd xmm7, [esi + eax*8 + 88]	;# G1 H1 	
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb330nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [esp + nb330nf_c12] ;# Vvdw12 
	
	addsd  xmm5, [esp + nb330nf_Vvdwtot]
	movlpd [esp + nb330nf_Vvdwtot], xmm5
	
.nb330nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb330nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb330nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb330nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb330nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb330nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb330nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb330nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb330nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb330nf_n], esi
        jmp .nb330nf_outer
.nb330nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb330nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb330nf_end
        ;# non-zero, do one more workunit
        jmp   .nb330nf_threadloop
.nb330nf_end:
	emms

	mov eax, [esp + nb330nf_nouter]
	mov ebx, [esp + nb330nf_ninner]
	mov ecx, [ebp + nb330nf_outeriter]
	mov edx, [ebp + nb330nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb330nf_salign]
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

