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

.globl nb_kernel233_x86_64_sse2
.globl _nb_kernel233_x86_64_sse2
nb_kernel233_x86_64_sse2:	
_nb_kernel233_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb233_fshift,           16
.equiv          nb233_gid,              24
.equiv          nb233_pos,              32
.equiv          nb233_faction,          40
.equiv          nb233_charge,           48
.equiv          nb233_p_facel,          56
.equiv          nb233_argkrf,           64
.equiv          nb233_argcrf,           72
.equiv          nb233_Vc,               80
.equiv          nb233_type,             88
.equiv          nb233_p_ntype,          96
.equiv          nb233_vdwparam,         104
.equiv          nb233_Vvdw,             112
.equiv          nb233_p_tabscale,       120
.equiv          nb233_VFtab,            128
.equiv          nb233_invsqrta,         136
.equiv          nb233_dvda,             144
.equiv          nb233_p_gbtabscale,     152
.equiv          nb233_GBtab,            160
.equiv          nb233_p_nthreads,       168
.equiv          nb233_count,            176
.equiv          nb233_mtx,              184
.equiv          nb233_outeriter,        192
.equiv          nb233_inneriter,        200
.equiv          nb233_work,             208
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
.equiv          nb233_facel,            976
.equiv          nb233_innerjjnr,        984	
.equiv          nb233_iinr,             992
.equiv          nb233_jindex,           1000
.equiv          nb233_jjnr,             1008
.equiv          nb233_shift,            1016
.equiv          nb233_shiftvec,         1024
.equiv          nb233_is3,              1032
.equiv          nb233_ii3,              1036
.equiv          nb233_nri,              1040
.equiv          nb233_ntia,             1044
.equiv          nb233_innerk,           1048
.equiv          nb233_n,                1052
.equiv          nb233_nn1,              1056
.equiv          nb233_nouter,           1060
.equiv          nb233_ninner,           1064

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
	sub rsp, 1072		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb233_nouter], eax
	mov [rsp + nb233_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb233_nri], edi
	mov [rsp + nb233_iinr], rsi
	mov [rsp + nb233_jindex], rdx
	mov [rsp + nb233_jjnr], rcx
	mov [rsp + nb233_shift], r8
	mov [rsp + nb233_shiftvec], r9
	mov rsi, [rbp + nb233_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb233_facel], xmm0

	mov rax, [rbp + nb233_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb233_tsc], xmm3

	mov rsi, [rbp + nb233_argkrf]
	mov rdi, [rbp + nb233_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb233_krf], xmm1
	movapd [rsp + nb233_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb233_half], eax
	mov [rsp + nb233_half+4], ebx
	movsd xmm1, [rsp + nb233_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb233_half], xmm1
	movapd [rsp + nb233_two], xmm2
	movapd [rsp + nb233_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb233_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb233_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	

	movsd xmm5, [rsp + nb233_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb233_iqH], xmm3
	movapd [rsp + nb233_iqM], xmm4
	
	mov   rdx, [rbp + nb233_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb233_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb233_ntia], ecx		
.nb233_threadloop:
        mov   rsi, [rbp + nb233_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb233_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb233_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb233_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb233_n], eax
        mov [rsp + nb233_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb233_outerstart
        jmp .nb233_end

.nb233_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb233_nouter]
	mov [rsp + nb233_nouter], ebx

.nb233_outer:
	mov   rax, [rsp + nb233_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb233_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb233_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb233_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb233_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb233_ii3], ebx

	addsd xmm3, [rax + rbx*8] 	;# ox
	addsd xmm4, [rax + rbx*8 + 8] 	;# oy
	addsd xmm5, [rax + rbx*8 + 16]	;# oz	
	addsd xmm6, [rax + rbx*8 + 24] 	;# h1x
	addsd xmm7, [rax + rbx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [rsp + nb233_ixO], xmm3
	movapd [rsp + nb233_iyO], xmm4
	movapd [rsp + nb233_izO], xmm5
	movapd [rsp + nb233_ixH1], xmm6
	movapd [rsp + nb233_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [rax + rbx*8 + 40] ;# h1z
	addsd xmm0, [rax + rbx*8 + 48] ;# h2x
	addsd xmm1, [rax + rbx*8 + 56] ;# h2y
	addsd xmm2, [rax + rbx*8 + 64] ;# h2z
	addsd xmm3, [rax + rbx*8 + 72] ;# mx
	addsd xmm4, [rax + rbx*8 + 80] ;# my
	addsd xmm5, [rax + rbx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb233_izH1], xmm6
	movapd [rsp + nb233_ixH2], xmm0
	movapd [rsp + nb233_iyH2], xmm1
	movapd [rsp + nb233_izH2], xmm2
	movapd [rsp + nb233_ixM], xmm3
	movapd [rsp + nb233_iyM], xmm4
	movapd [rsp + nb233_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb233_vctot], xmm4
	movapd [rsp + nb233_Vvdwtot], xmm4
	movapd [rsp + nb233_fixO], xmm4
	movapd [rsp + nb233_fiyO], xmm4
	movapd [rsp + nb233_fizO], xmm4
	movapd [rsp + nb233_fixH1], xmm4
	movapd [rsp + nb233_fiyH1], xmm4
	movapd [rsp + nb233_fizH1], xmm4
	movapd [rsp + nb233_fixH2], xmm4
	movapd [rsp + nb233_fiyH2], xmm4
	movapd [rsp + nb233_fizH2], xmm4
	movapd [rsp + nb233_fixM], xmm4
	movapd [rsp + nb233_fiyM], xmm4
	movapd [rsp + nb233_fizM], xmm4
	
	mov   rax, [rsp + nb233_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb233_pos]
	mov   rdi, [rbp + nb233_faction]	
	mov   rax, [rsp + nb233_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb233_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb233_ninner]
	mov   [rsp + nb233_ninner], ecx
	add   edx, 0
	mov   [rsp + nb233_innerk], edx    ;# number of innerloop atoms 
	jge   .nb233_unroll_loop
	jmp   .nb233_checksingle
.nb233_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb233_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb233_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb233_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb233_iqM]
	mulpd  xmm4, [rsp + nb233_iqH]

	movapd  [rsp + nb233_qqM], xmm3
	movapd  [rsp + nb233_qqH], xmm4
	
	mov rsi, [rbp + nb233_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb233_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb233_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r9*8]	;# c6b
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + r9*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movapd [rsp + nb233_c6], xmm4
	movapd [rsp + nb233_c12], xmm6
	
	mov rsi, [rbp + nb233_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load j coordinates 
    movlpd xmm0, [rsi + rax*8] 
    movlpd xmm1, [rsi + rax*8 + 8] 
    movlpd xmm2, [rsi + rax*8 + 16] 
    movhpd xmm0, [rsi + rbx*8] 
    movhpd xmm1, [rsi + rbx*8 + 8] 
    movhpd xmm2, [rsi + rbx*8 + 16] 

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    ;# O interaction
    ;# copy to xmm3-xmm5
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    
    subpd xmm3, [rsp + nb233_ixO]
    subpd xmm4, [rsp + nb233_iyO]
    subpd xmm5, [rsp + nb233_izO]
    
    movapd [rsp + nb233_dxO], xmm3
    movapd [rsp + nb233_dyO], xmm4
    movapd [rsp + nb233_dzO], xmm5
    
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5

	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    ;# xmm3=rsq

    cvtpd2ps xmm5, xmm3     
    rsqrtps xmm5, xmm5
    cvtps2pd xmm15, xmm5     ;# lu in low xmm2 

    ;# lookup seed in xmm2 
    movapd xmm5, xmm15       ;# copy of lu 
    mulpd xmm15, xmm15        ;# lu*lu 
    movapd xmm7, [rsp + nb233_three]
    mulpd xmm15, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb233_half]
    subpd xmm7, xmm15        ;# 30-rsq*lu*lu 
    mulpd xmm7, xmm5        
    mulpd xmm7, xmm6        ;# xmm0=iter1 of rinv (new lu) 

    movapd xmm5, xmm7       ;# copy of lu 
    mulpd xmm7, xmm7        ;# lu*lu 
    movapd xmm15, [rsp + nb233_three]
    mulpd xmm7, xmm3        ;# rsq*lu*lu                    
    movapd xmm6, [rsp + nb233_half]
    subpd xmm15, xmm7        ;# 30-rsq*lu*lu 
    mulpd xmm15, xmm5        
    mulpd xmm15, xmm6        ;# xmm15=rinv
        
    mulpd xmm3, xmm15        ;# xmm3=r 

    ;# xmm15=rinv
    ;# xmm3=r

    mulpd xmm3, [rsp + nb233_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttpd2pi mm6, xmm3
    
    ;# convert back to float
    cvtpi2pd  xmm4, mm6
    
    ;# multiply by 8
    pslld   mm6, 3

    ;# calculate eps
    subpd     xmm3, xmm4    ;# xmm3=eps
    
    ;# move to integer registers
    movd r10d, mm6
    psrlq mm6, 32
    movd r11d, mm6
    
    ;# xmm3=eps
    ;# xmm15=rinv

	mov rsi, [rbp + nb233_VFtab]
    ;# indices in r10, r11. Load dispersion and repulsion tables in parallel.
    movapd xmm4, [rsi + r10*8]          ;# Y1d F1d  
    movapd xmm12, [rsi + r11*8]         ;# Y2d F2d 
    movapd xmm8, [rsi + r10*8 + 32]     ;# Y1r F1r  
    movapd xmm13, [rsi + r11*8 + 32]    ;# Y2r F2r 
    movapd xmm5, xmm4
    movapd xmm9, xmm8
    unpcklpd xmm4, xmm12    ;# Y1d Y2d 
    unpckhpd xmm5, xmm12    ;# F1d F2d 
    unpcklpd xmm8, xmm13    ;# Y1r Y2r 
    unpckhpd xmm9, xmm13    ;# F1r F2r 

    movapd xmm6, [rsi + r10*8 + 16]     ;# G1d H1d  
    movapd xmm12, [rsi + r11*8 + 16]        ;# G2d H2d 
    movapd xmm10, [rsi + r10*8 + 48]        ;# G1r H1r      
    movapd xmm13, [rsi + r11*8 + 48]        ;# G2r H2r 
    movapd xmm7, xmm6
    movapd xmm11, xmm10
    unpcklpd xmm6, xmm12    ;# G1d G2d 
    unpckhpd xmm7, xmm12    ;# H1d H2d 
    unpcklpd xmm10, xmm13   ;# G1r G2r 
    unpckhpd xmm11, xmm13   ;# H1r H2r 
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    mulpd  xmm7, xmm3    ;# Heps
    mulpd  xmm11, xmm3 
    mulpd  xmm6, xmm3   ;# Geps
    mulpd  xmm10, xmm3 
    mulpd  xmm7, xmm3   ;# Heps2
    mulpd  xmm11, xmm3 
    addpd  xmm5, xmm6  ;# F+Geps
    addpd  xmm9, xmm10 
    addpd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addpd  xmm9, xmm11 
    addpd  xmm7, xmm7    ;# 2*Heps2
    addpd  xmm11, xmm11
    addpd  xmm7, xmm6   ;# 2*Heps2+Geps
    addpd  xmm11, xmm10
    
    addpd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm11, xmm9
    mulpd  xmm5, xmm3  ;# eps*Fp
    mulpd  xmm9, xmm3
    movapd xmm12, [rsp + nb233_c6]
    movapd xmm13, [rsp + nb233_c12]
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb233_Vvdwtot]
    movapd [rsp + nb233_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    
    mulpd  xmm7, [rsp + nb233_tsc]
    mulpd  xmm7, xmm15   ;# -fscal
    xorpd  xmm9, xmm9
    
    subpd  xmm9, xmm7     ;# fscal
    movapd xmm10, xmm9
    movapd xmm11, xmm9

    mulpd  xmm9,  [rsp + nb233_dxO] ;# fx/fy/fz
    mulpd  xmm10, [rsp + nb233_dyO]
    mulpd  xmm11, [rsp + nb233_dzO]

    ;# save j force temporarily
    movapd [rsp + nb233_fjx], xmm9
    movapd [rsp + nb233_fjy], xmm10
    movapd [rsp + nb233_fjz], xmm11
    
    ;# increment i O force
    addpd xmm9, [rsp + nb233_fixO]
    addpd xmm10, [rsp + nb233_fiyO]
    addpd xmm11, [rsp + nb233_fizO]
    movapd [rsp + nb233_fixO], xmm9
    movapd [rsp + nb233_fiyO], xmm10
    movapd [rsp + nb233_fizO], xmm11
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.                
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb233_ixH1]
    subpd xmm1, [rsp + nb233_iyH1]
    subpd xmm2, [rsp + nb233_izH1]
    subpd xmm3, [rsp + nb233_ixH2]
    subpd xmm4, [rsp + nb233_iyH2]
    subpd xmm5, [rsp + nb233_izH2]
    subpd xmm6, [rsp + nb233_ixM]
    subpd xmm7, [rsp + nb233_iyM]
    subpd xmm8, [rsp + nb233_izM]
    
	movapd [rsp + nb233_dxH1], xmm0
	movapd [rsp + nb233_dyH1], xmm1
	movapd [rsp + nb233_dzH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb233_dxH2], xmm3
	movapd [rsp + nb233_dyH2], xmm4
	movapd [rsp + nb233_dzH2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb233_dxM], xmm6
	movapd [rsp + nb233_dyM], xmm7
	movapd [rsp + nb233_dzM], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for j atoms
    cvtpd2ps xmm1, xmm0
    cvtpd2ps xmm4, xmm3
    cvtpd2ps xmm7, xmm6
	rsqrtps xmm1, xmm1
	rsqrtps xmm4, xmm4
    rsqrtps xmm7, xmm7
    cvtps2pd xmm1, xmm1
    cvtps2pd xmm4, xmm4
    cvtps2pd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulpd   xmm1, xmm1 ;# lu*lu
	mulpd   xmm4, xmm4 ;# lu*lu
    mulpd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb233_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulpd   xmm1, xmm0 ;# rsq*lu*lu
	mulpd   xmm4, xmm3 ;# rsq*lu*lu 
    mulpd   xmm7, xmm6 ;# rsq*lu*lu
	
	subpd   xmm9, xmm1
	subpd   xmm10, xmm4
    subpd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm2
	mulpd   xmm10, xmm5
    mulpd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb233_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2
    mulpd   xmm11, xmm15 ;# first iteration for rinvM

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb233_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulpd   xmm2, xmm0 ;# rsq*lu*lu
	mulpd   xmm5, xmm3 ;# rsq*lu*lu 
    mulpd   xmm8, xmm6 ;# rsq*lu*lu
	
	subpd   xmm1, xmm2
	subpd   xmm4, xmm5
    subpd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm1
	mulpd   xmm10, xmm4
    mulpd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb233_half]
	mulpd   xmm9, xmm15  ;#  rinvH1
	mulpd   xmm10, xmm15 ;#   rinvH2
    mulpd   xmm11, xmm15 ;#   rinvM
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, xmm9 ;# copy of rinv
    movapd xmm4, xmm10
    movapd xmm7, xmm11
    movapd xmm2, [rsp + nb233_krf]    
    mulpd  xmm9, xmm9   ;# rinvsq
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, xmm2  ;# k*rsq
    mulpd  xmm3, xmm2
    mulpd  xmm6, xmm2
    movapd xmm2, xmm0 ;# copy of k*rsq
    movapd xmm5, xmm3
    movapd xmm8, xmm6
    addpd  xmm2, xmm1  ;# rinv+krsq
    addpd  xmm5, xmm4
    addpd  xmm8, xmm7
    movapd xmm14, [rsp + nb233_crf]
    subpd  xmm2, xmm14   ;# rinv+krsq-crf
    subpd  xmm5, xmm14
    subpd  xmm8, xmm14
    movapd xmm12, [rsp + nb233_qqH]
    movapd xmm13, [rsp + nb233_qqM]    
    mulpd  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm5, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addpd  xmm0, xmm0 ;# 2*krsq
    addpd  xmm3, xmm3 
    addpd  xmm6, xmm6 
    subpd  xmm1, xmm0 ;# rinv-2*krsq
    subpd  xmm4, xmm3
    subpd  xmm7, xmm6
    mulpd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulpd  xmm4, xmm12
    mulpd  xmm7, xmm13
    addpd  xmm2, [rsp + nb233_vctot]
    addpd  xmm5, xmm8
    addpd  xmm2, xmm5
    movapd [rsp + nb233_vctot], xmm2
    
    mulpd  xmm9, xmm1   ;# fscal
    mulpd  xmm10, xmm4
    mulpd  xmm11, xmm7

    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb233_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

    ;# add forces from O interaction
    addpd xmm0, [rsp + nb233_fjx]
    addpd xmm1, [rsp + nb233_fjy]
    addpd xmm2, [rsp + nb233_fjz]

	mulpd xmm7, [rsp + nb233_dxH1]
	mulpd xmm8, [rsp + nb233_dyH1]
	mulpd xmm9, [rsp + nb233_dzH1]
	mulpd xmm10, [rsp + nb233_dxH2]
	mulpd xmm11, [rsp + nb233_dyH2]
	mulpd xmm12, [rsp + nb233_dzH2]
	mulpd xmm13, [rsp + nb233_dxM]
	mulpd xmm14, [rsp + nb233_dyM]
	mulpd xmm15, [rsp + nb233_dzM]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb233_fixH1]
    addpd xmm8, [rsp + nb233_fiyH1]
    addpd xmm9, [rsp + nb233_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb233_fixH2]
    addpd xmm11, [rsp + nb233_fiyH2]
    addpd xmm12, [rsp + nb233_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb233_fixM]
    addpd xmm14, [rsp + nb233_fiyM]
    addpd xmm15, [rsp + nb233_fizM]

    movapd [rsp + nb233_fixH1], xmm7
    movapd [rsp + nb233_fiyH1], xmm8
    movapd [rsp + nb233_fizH1], xmm9
    movapd [rsp + nb233_fixH2], xmm10
    movapd [rsp + nb233_fiyH2], xmm11
    movapd [rsp + nb233_fizH2], xmm12
    movapd [rsp + nb233_fixM], xmm13
    movapd [rsp + nb233_fiyM], xmm14
    movapd [rsp + nb233_fizM], xmm15
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb233_innerk],  2
	jl   .nb233_checksingle
	jmp  .nb233_unroll_loop
.nb233_checksingle:	
	mov   edx, [rsp + nb233_innerk]
	and   edx, 1
	jnz  .nb233_dosingle
	jmp  .nb233_updateouterdata
.nb233_dosingle:
	mov   rdx, [rsp + nb233_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb233_innerjjnr],  4

	mov rsi, [rbp + nb233_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [rsp + nb233_iqM]
	mulsd  xmm4, [rsp + nb233_iqH]

	movapd  [rsp + nb233_qqM], xmm3
	movapd  [rsp + nb233_qqH], xmm4
	
	mov rsi, [rbp + nb233_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb233_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb233_ntia]
	add r8d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb233_c6], xmm4
	movapd [rsp + nb233_c12], xmm6
	
	mov rsi, [rbp + nb233_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2  and xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb233_ixO]
	subsd xmm5, [rsp + nb233_iyO]
	subsd xmm6, [rsp + nb233_izO]

	;# store dr 
	movapd [rsp + nb233_dxO], xmm4
	movapd [rsp + nb233_dyO], xmm5
	movapd [rsp + nb233_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [rsp + nb233_rsqO], xmm7
	
	;# move j coords to xmm4-xmm6 
	movapd xmm4, xmm0
	movapd xmm5, xmm1
	movapd xmm6, xmm2

	;# calc dr 
	subsd xmm4, [rsp + nb233_ixH1]
	subsd xmm5, [rsp + nb233_iyH1]
	subsd xmm6, [rsp + nb233_izH1]

	;# store dr 
	movapd [rsp + nb233_dxH1], xmm4
	movapd [rsp + nb233_dyH1], xmm5
	movapd [rsp + nb233_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move j coords to xmm3-xmm5 
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	;# calc dr 
	subsd xmm3, [rsp + nb233_ixH2]
	subsd xmm4, [rsp + nb233_iyH2]
	subsd xmm5, [rsp + nb233_izH2]

	;# store dr 
	movapd [rsp + nb233_dxH2], xmm3
	movapd [rsp + nb233_dyH2], xmm4
	movapd [rsp + nb233_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3

	;# move j coords to xmm4-xmm2
	movapd xmm4, xmm0
	movapd xmm3, xmm1
    ;# xmm2 already contains z

	;# calc dr 
	subsd xmm4, [rsp + nb233_ixM]
	subsd xmm3, [rsp + nb233_iyM]
	subsd xmm2, [rsp + nb233_izM]

	;# store dr 
	movapd [rsp + nb233_dxM], xmm4
	movapd [rsp + nb233_dyM], xmm3
	movapd [rsp + nb233_dzM], xmm2

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [rsp + nb233_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [rsp + nb233_krsqM], xmm0
	movsd [rsp + nb233_krsqH2], xmm1
	movsd [rsp + nb233_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb233_half] ;# rinv 
	movapd [rsp + nb233_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb233_half] ;# rinv 
	movapd [rsp + nb233_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb233_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb233_half] ;# rinv 
	movapd [rsp + nb233_rinvM], xmm1

	;# rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movsd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movsd  xmm4, [rsp + nb233_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb233_half] ;# iter1 ( new lu) 

	movsd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movsd xmm4, [rsp + nb233_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb233_half] ;# rinv 
	movsd  xmm7, xmm4	;# rinvO in xmm7 
	
	movsd xmm4, [rsp + nb233_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulsd xmm4, xmm7	;# xmm4=r 
	mulsd xmm4, [rsp + nb233_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3

	mov  rsi, [rbp + nb233_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb233_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [rsp + nb233_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb233_Vvdwtot]
	xorpd  xmm3, xmm3
	mulsd  xmm7, [rsp + nb233_tsc]
	subsd  xmm3, xmm7
	movsd [rsp + nb233_fstmp], xmm3
	movsd [rsp + nb233_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb233_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [rsp + nb233_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb233_Vvdwtot]
	movsd xmm3, [rsp + nb233_fstmp]
	mulsd  xmm7, [rsp + nb233_tsc]
	subsd  xmm3, xmm7
	movsd [rsp + nb233_Vvdwtot], xmm5

	mulsd  xmm3, xmm0
		
		
	movsd xmm0, [rsp + nb233_dxO]
	movsd xmm1, [rsp + nb233_dyO]
	movsd xmm2, [rsp + nb233_dzO]

	mov    rdi, [rbp + nb233_faction]
	mulsd  xmm0, xmm3
	mulsd  xmm1, xmm3
	mulsd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [rsp + nb233_fixO]
	movapd xmm4, [rsp + nb233_fiyO]
	movapd xmm7, [rsp + nb233_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb233_fixO], xmm3
	movsd [rsp + nb233_fiyO], xmm4
	movsd [rsp + nb233_fizO], xmm7
	;# update j forces with water O 
	movsd [rsp + nb233_fjx], xmm0
	movsd [rsp + nb233_fjy], xmm1
	movsd [rsp + nb233_fjz], xmm2

	;# H1 interactions
	movsd  xmm6, [rsp + nb233_rinvH1] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm6
	movsd  xmm0, [rsp + nb233_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [rsp + nb233_two]
	subsd   xmm6, [rsp + nb233_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [rsp + nb233_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb233_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [rsp + nb233_vctot]

	movapd xmm0, [rsp + nb233_dxH1]
	movapd xmm1, [rsp + nb233_dyH1]
	movapd xmm2, [rsp + nb233_dzH1]
	movsd [rsp + nb233_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb233_fixH1]
	movapd xmm4, [rsp + nb233_fiyH1]
	movapd xmm7, [rsp + nb233_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb233_fixH1], xmm3
	movsd [rsp + nb233_fiyH1], xmm4
	movsd [rsp + nb233_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb233_fjx]
	addsd  xmm1, [rsp + nb233_fjy]
	addsd  xmm2, [rsp + nb233_fjz]
	movsd [rsp + nb233_fjx], xmm0
	movsd [rsp + nb233_fjy], xmm1
	movsd [rsp + nb233_fjz], xmm2

	;# H2 interactions 
	movsd  xmm5, [rsp + nb233_rinvH2] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [rsp + nb233_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [rsp + nb233_two]
	subsd   xmm5, [rsp + nb233_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [rsp + nb233_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb233_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [rsp + nb233_vctot]

	movapd xmm0, [rsp + nb233_dxH2]
	movapd xmm1, [rsp + nb233_dyH2]
	movapd xmm2, [rsp + nb233_dzH2]
	movsd [rsp + nb233_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb233_fixH2]
	movapd xmm4, [rsp + nb233_fiyH2]
	movapd xmm7, [rsp + nb233_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb233_fixH2], xmm3
	movsd [rsp + nb233_fiyH2], xmm4
	movsd [rsp + nb233_fizH2], xmm7
	;# update j forces with water H2 
	addsd  xmm0, [rsp + nb233_fjx]
	addsd  xmm1, [rsp + nb233_fjy]
	addsd  xmm2, [rsp + nb233_fjz]
	movsd [rsp + nb233_fjx], xmm0
	movsd [rsp + nb233_fjy], xmm1
	movsd [rsp + nb233_fjz], xmm2

	;# M interactions 
	movsd  xmm5, [rsp + nb233_rinvM] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [rsp + nb233_krsqM]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [rsp + nb233_two]
	subsd   xmm5, [rsp + nb233_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [rsp + nb233_qqM] ;# vcoul 
	mulsd   xmm7, [rsp + nb233_qqM]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [rsp + nb233_vctot]

	movapd xmm0, [rsp + nb233_dxM]
	movapd xmm1, [rsp + nb233_dyM]
	movapd xmm2, [rsp + nb233_dzM]
	movsd [rsp + nb233_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [rsp + nb233_fixM]
	movapd xmm4, [rsp + nb233_fiyM]
	movapd xmm7, [rsp + nb233_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb233_fixM], xmm3
	movsd [rsp + nb233_fiyM], xmm4
	movsd [rsp + nb233_fizM], xmm7

	mov rdi, [rbp + nb233_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb233_fjx]
	addsd  xmm1, [rsp + nb233_fjy]
	addsd  xmm2, [rsp + nb233_fjz]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb233_updateouterdata:
	mov   ecx, [rsp + nb233_ii3]
	mov   rdi, [rbp + nb233_faction]
	mov   rsi, [rbp + nb233_fshift]
	mov   edx, [rsp + nb233_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb233_fixO]
	movapd xmm1, [rsp + nb233_fiyO]
	movapd xmm2, [rsp + nb233_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6,xmm1 

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb233_fixH1]
	movapd xmm1, [rsp + nb233_fiyH1]
	movapd xmm2, [rsp + nb233_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 24]
	movsd  xmm4, [rdi + rcx*8 + 32]
	movsd  xmm5, [rdi + rcx*8 + 40]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 24], xmm3
	movsd  [rdi + rcx*8 + 32], xmm4
	movsd  [rdi + rcx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb233_fixH2]
	movapd xmm1, [rsp + nb233_fiyH2]
	movapd xmm2, [rsp + nb233_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 48]
	movsd  xmm4, [rdi + rcx*8 + 56]
	movsd  xmm5, [rdi + rcx*8 + 64]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 48], xmm3
	movsd  [rdi + rcx*8 + 56], xmm4
	movsd  [rdi + rcx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb233_fixM]
	movapd xmm1, [rsp + nb233_fiyM]
	movapd xmm2, [rsp + nb233_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 72]
	movsd  xmm4, [rdi + rcx*8 + 80]
	movsd  xmm5, [rdi + rcx*8 + 88]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 72], xmm3
	movsd  [rdi + rcx*8 + 80], xmm4
	movsd  [rdi + rcx*8 + 88], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [rsi + rdx*8]
	movhpd xmm3, [rsi + rdx*8 + 8]
	movsd  xmm4, [rsi + rdx*8 + 16]
	subpd  xmm3, xmm6
	subsd  xmm4, xmm7
	movlpd [rsi + rdx*8],      xmm3
	movhpd [rsi + rdx*8 + 8],  xmm3
	movsd  [rsi + rdx*8 + 16], xmm4

	;# get n from stack
	mov esi, [rsp + nb233_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb233_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb233_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb233_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb233_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb233_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb233_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb233_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb233_n], esi
        jmp .nb233_outer
.nb233_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb233_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb233_end
        ;# non-zero, do one more workunit
        jmp   .nb233_threadloop
.nb233_end:
	mov eax, [rsp + nb233_nouter]
	mov ebx, [rsp + nb233_ninner]
	mov rcx, [rbp + nb233_outeriter]
	mov rdx, [rbp + nb233_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1072
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret




.globl nb_kernel233nf_x86_64_sse2
.globl _nb_kernel233nf_x86_64_sse2
nb_kernel233nf_x86_64_sse2:	
_nb_kernel233nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb233nf_fshift,         16
.equiv          nb233nf_gid,            24
.equiv          nb233nf_pos,            32
.equiv          nb233nf_faction,        40
.equiv          nb233nf_charge,         48
.equiv          nb233nf_p_facel,        56
.equiv          nb233nf_argkrf,         64
.equiv          nb233nf_argcrf,         72
.equiv          nb233nf_Vc,             80
.equiv          nb233nf_type,           88
.equiv          nb233nf_p_ntype,        96
.equiv          nb233nf_vdwparam,       104
.equiv          nb233nf_Vvdw,           112
.equiv          nb233nf_p_tabscale,     120
.equiv          nb233nf_VFtab,          128
.equiv          nb233nf_invsqrta,       136
.equiv          nb233nf_dvda,           144
.equiv          nb233nf_p_gbtabscale,   152
.equiv          nb233nf_GBtab,          160
.equiv          nb233nf_p_nthreads,     168
.equiv          nb233nf_count,          176
.equiv          nb233nf_mtx,            184
.equiv          nb233nf_outeriter,      192
.equiv          nb233nf_inneriter,      200
.equiv          nb233nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb233nf_ixO,            0
.equiv          nb233nf_iyO,            16
.equiv          nb233nf_izO,            32
.equiv          nb233nf_ixH1,           48
.equiv          nb233nf_iyH1,           64
.equiv          nb233nf_izH1,           80
.equiv          nb233nf_ixH2,           96
.equiv          nb233nf_iyH2,           112
.equiv          nb233nf_izH2,           128
.equiv          nb233nf_ixM,            144
.equiv          nb233nf_iyM,            160
.equiv          nb233nf_izM,            176
.equiv          nb233nf_iqH,            192
.equiv          nb233nf_iqM,            208
.equiv          nb233nf_qqH,            224
.equiv          nb233nf_qqM,            240
.equiv          nb233nf_c6,             256
.equiv          nb233nf_c12,            272
.equiv          nb233nf_vctot,          288
.equiv          nb233nf_Vvdwtot,        304
.equiv          nb233nf_half,           320
.equiv          nb233nf_three,          336
.equiv          nb233nf_tsc,            352
.equiv          nb233nf_rinvH1,         368
.equiv          nb233nf_rinvH2,         384
.equiv          nb233nf_rinvM,          400
.equiv          nb233nf_krsqH1,         416
.equiv          nb233nf_krsqH2,         432
.equiv          nb233nf_krsqM,          448
.equiv          nb233nf_krf,            464
.equiv          nb233nf_crf,            480
.equiv          nb233nf_rsqO,           496
.equiv          nb233nf_facel,          512
.equiv          nb233nf_iinr,           520
.equiv          nb233nf_jindex,         528
.equiv          nb233nf_jjnr,           536
.equiv          nb233nf_shift,          544
.equiv          nb233nf_shiftvec,       552
.equiv          nb233nf_innerjjnr,      560
.equiv          nb233nf_nri,            568
.equiv          nb233nf_is3,            572
.equiv          nb233nf_ii3,            576
.equiv          nb233nf_ntia,           580
.equiv          nb233nf_innerk,         584
.equiv          nb233nf_n,              588
.equiv          nb233nf_nn1,            592
.equiv          nb233nf_nouter,         596
.equiv          nb233nf_ninner,         600

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
	sub rsp, 608		;# local variable stack space (n*16+8)
	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb233nf_nouter], eax
	mov [rsp + nb233nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb233nf_nri], edi
	mov [rsp + nb233nf_iinr], rsi
	mov [rsp + nb233nf_jindex], rdx
	mov [rsp + nb233nf_jjnr], rcx
	mov [rsp + nb233nf_shift], r8
	mov [rsp + nb233nf_shiftvec], r9
	mov rsi, [rbp + nb233nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb233nf_facel], xmm0

	mov rax, [rbp + nb233nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb233nf_tsc], xmm3

	mov rsi, [rbp + nb233nf_argkrf]
	mov rdi, [rbp + nb233nf_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb233nf_krf], xmm1
	movapd [rsp + nb233nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb233nf_half], eax
	mov [rsp + nb233nf_half+4], ebx
	movsd xmm1, [rsp + nb233nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb233nf_half], xmm1
	movapd [rsp + nb233nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb233nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb233nf_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	

	movsd xmm5, [rsp + nb233nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb233nf_iqH], xmm3
	movapd [rsp + nb233nf_iqM], xmm4
	
	mov   rdx, [rbp + nb233nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb233nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb233nf_ntia], ecx		
.nb233nf_threadloop:
        mov   rsi, [rbp + nb233nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb233nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb233nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb233nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb233nf_n], eax
        mov [rsp + nb233nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb233nf_outerstart
        jmp .nb233nf_end

.nb233nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb233nf_nouter]
	mov [rsp + nb233nf_nouter], ebx

.nb233nf_outer:
	mov   rax, [rsp + nb233nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb233nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb233nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb233nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb233nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb233nf_ii3], ebx

	addsd xmm3, [rax + rbx*8] 	;# ox
	addsd xmm4, [rax + rbx*8 + 8] 	;# oy
	addsd xmm5, [rax + rbx*8 + 16]	;# oz	
	addsd xmm6, [rax + rbx*8 + 24] 	;# h1x
	addsd xmm7, [rax + rbx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [rsp + nb233nf_ixO], xmm3
	movapd [rsp + nb233nf_iyO], xmm4
	movapd [rsp + nb233nf_izO], xmm5
	movapd [rsp + nb233nf_ixH1], xmm6
	movapd [rsp + nb233nf_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [rax + rbx*8 + 40] ;# h1z
	addsd xmm0, [rax + rbx*8 + 48] ;# h2x
	addsd xmm1, [rax + rbx*8 + 56] ;# h2y
	addsd xmm2, [rax + rbx*8 + 64] ;# h2z
	addsd xmm3, [rax + rbx*8 + 72] ;# mx
	addsd xmm4, [rax + rbx*8 + 80] ;# my
	addsd xmm5, [rax + rbx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb233nf_izH1], xmm6
	movapd [rsp + nb233nf_ixH2], xmm0
	movapd [rsp + nb233nf_iyH2], xmm1
	movapd [rsp + nb233nf_izH2], xmm2
	movapd [rsp + nb233nf_ixM], xmm3
	movapd [rsp + nb233nf_iyM], xmm4
	movapd [rsp + nb233nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb233nf_vctot], xmm4
	movapd [rsp + nb233nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb233nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb233nf_pos]
	mov   rdi, [rbp + nb233nf_faction]	
	mov   rax, [rsp + nb233nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb233nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb233nf_ninner]
	mov   [rsp + nb233nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb233nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb233nf_unroll_loop
	jmp   .nb233nf_checksingle
.nb233nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb233nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb233nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb233nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb233nf_iqM]
	mulpd  xmm4, [rsp + nb233nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [rsp + nb233nf_qqM], xmm3
	movapd  [rsp + nb233nf_qqH], xmm4
	
	mov rsi, [rbp + nb233nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb233nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb233nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb233nf_c6], xmm4
	movapd [rsp + nb233nf_c12], xmm6
	
	mov rsi, [rbp + nb233nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb233nf_ixO]
	movapd xmm5, [rsp + nb233nf_iyO]
	movapd xmm6, [rsp + nb233nf_izO]

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
	movapd xmm4, [rsp + nb233nf_ixH1]
	movapd xmm5, [rsp + nb233nf_iyH1]
	movapd xmm6, [rsp + nb233nf_izH1]

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
	movapd xmm3, [rsp + nb233nf_ixH2]
	movapd xmm4, [rsp + nb233nf_iyH2]
	movapd xmm5, [rsp + nb233nf_izH2]

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
	movapd xmm3, [rsp + nb233nf_iyM]
	movapd xmm4, [rsp + nb233nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb233nf_ixM]
	subpd  xmm2, xmm0	


	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	movapd [rsp + nb233nf_rsqO], xmm7
	
	;# calculate krsq
	movapd xmm0, [rsp + nb233nf_krf]
	movapd xmm1, xmm0
	movapd xmm2, xmm0
	mulpd xmm0, xmm4  
	mulpd xmm1, xmm5
	mulpd xmm2, xmm6
	movapd [rsp + nb233nf_krsqM], xmm0
	movapd [rsp + nb233nf_krsqH2], xmm1
	movapd [rsp + nb233nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2
	
	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233nf_three]
	subpd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb233nf_half] ;# rinv 
	movapd  [rsp + nb233nf_rinvH1], xmm1	

	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233nf_three]
	subpd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb233nf_half] ;# rinv 
	movapd  [rsp + nb233nf_rinvH2], xmm1	
	
	;# rsqM - seed in xmm2 
	cvtpd2ps xmm2, xmm4	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233nf_three]
	mulpd   xmm2, xmm4	;# rsq*lu*lu 
	subpd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm1, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulpd xmm1, xmm1	;# lu*lu 
	mulpd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233nf_three]
	subpd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulpd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm1, [rsp + nb233nf_half] ;# rinv 
	movapd  [rsp + nb233nf_rinvM], xmm1	

		
	;# rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb233nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb233nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb233nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	
	
	movapd xmm4, [rsp + nb233nf_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulpd xmm4, xmm7	;# xmm4=r 
	mulpd xmm4, [rsp + nb233nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  rsi, [rbp + nb233nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [rsi + rax*8]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 16]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 16]	;# G2
	movhpd xmm6, [rsi + rax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 24]	;# G2 H2 
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

	movapd xmm4, [rsp + nb233nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [rsp + nb233nf_Vvdwtot]
	movapd [rsp + nb233nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rax*8 + 32]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8 + 32]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 48]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 48]	;# G2
	movhpd xmm6, [rsi + rax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 56]	;# G2 H2 

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
	
	movapd xmm4, [rsp + nb233nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [rsp + nb233nf_Vvdwtot]
	movapd [rsp + nb233nf_Vvdwtot], xmm5

	;# H1 interactions 
	movapd  xmm6, [rsp + nb233nf_rinvH1] 
	movapd  xmm4, xmm6
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [rsp + nb233nf_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subpd   xmm6, [rsp + nb233nf_crf]
	mulpd   xmm6, [rsp + nb233nf_qqH] ;# vcoul 
	addpd  xmm6, [rsp + nb233nf_vctot]
	movapd [rsp + nb233nf_vctot], xmm6

	;# H2 interactions 
	movapd  xmm5, [rsp + nb233nf_rinvH2] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [rsp + nb233nf_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [rsp + nb233nf_crf]
	mulpd   xmm5, [rsp + nb233nf_qqH] ;# vcoul 
	addpd  xmm5, [rsp + nb233nf_vctot]
	movapd [rsp + nb233nf_vctot], xmm5

	;# M interactions 
	movapd  xmm5, [rsp + nb233nf_rinvM] 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [rsp + nb233nf_krsqM]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [rsp + nb233nf_crf]
	mulpd   xmm5, [rsp + nb233nf_qqM] ;# vcoul 
	addpd  xmm5, [rsp + nb233nf_vctot]
	movapd [rsp + nb233nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb233nf_innerk],  2
	jl   .nb233nf_checksingle
	jmp  .nb233nf_unroll_loop
.nb233nf_checksingle:	
	mov   edx, [rsp + nb233nf_innerk]
	and   edx, 1
	jnz  .nb233nf_dosingle
	jmp  .nb233nf_updateouterdata
.nb233nf_dosingle:
	mov   rdx, [rsp + nb233nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb233nf_innerjjnr],  4

	mov rsi, [rbp + nb233nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulsd  xmm3, [rsp + nb233nf_iqM]
	mulsd  xmm4, [rsp + nb233nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [rsp + nb233nf_qqM], xmm3
	movapd  [rsp + nb233nf_qqH], xmm4
	
	mov rsi, [rbp + nb233nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb233nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb233nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb233nf_c6], xmm4
	movapd [rsp + nb233nf_c12], xmm6
	
	mov rsi, [rbp + nb233nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb233nf_ixO]
	movapd xmm5, [rsp + nb233nf_iyO]
	movapd xmm6, [rsp + nb233nf_izO]

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
	movapd [rsp + nb233nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb233nf_ixH1]
	movapd xmm5, [rsp + nb233nf_iyH1]
	movapd xmm6, [rsp + nb233nf_izH1]

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
	movapd xmm3, [rsp + nb233nf_ixH2]
	movapd xmm4, [rsp + nb233nf_iyH2]
	movapd xmm5, [rsp + nb233nf_izH2]

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
	movapd xmm3, [rsp + nb233nf_iyM]
	movapd xmm4, [rsp + nb233nf_izM]
	subpd  xmm3, xmm1
	subpd  xmm4, xmm2
	movapd xmm2, [rsp + nb233nf_ixM]
	subpd  xmm2, xmm0	

	;# square it 
	mulpd xmm2,xmm2
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	addpd xmm4, xmm3
	addpd xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# calculate krsq
	movsd xmm0, [rsp + nb233nf_krf]
	movsd xmm1, xmm0
	movsd xmm2, xmm0
	mulsd xmm0, xmm4  
	mulsd xmm1, xmm5
	mulsd xmm2, xmm6
	movsd [rsp + nb233nf_krsqM], xmm0
	movsd [rsp + nb233nf_krsqH2], xmm1
	movsd [rsp + nb233nf_krsqH1], xmm2

	;# start with rsqH1 - put seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm6, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233nf_three]
	subsd xmm1, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb233nf_half] ;# rinv 
	movapd [rsp + nb233nf_rinvH1], xmm1
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm5, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233nf_three]
	subsd xmm1, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb233nf_half] ;# rinv 
	movapd [rsp + nb233nf_rinvH2], xmm1
	
	;# rsqM - seed in xmm2 
	cvtsd2ss xmm2, xmm4
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm1, [rsp + nb233nf_three]
	mulsd   xmm2, xmm4	;# rsq*lu*lu 
	subsd   xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm1, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm1, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm1
	mulsd xmm1, xmm1	;# lu*lu 
	mulsd xmm4, xmm1	;# rsq*lu*lu 
	movapd xmm1, [rsp + nb233nf_three]
	subsd xmm1, xmm4	;# 3-rsq*lu*lu 
	mulsd xmm1, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm1, [rsp + nb233nf_half] ;# rinv 
	movapd [rsp + nb233nf_rinvM], xmm1

	;# rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movsd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movsd  xmm4, [rsp + nb233nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb233nf_half] ;# iter1 ( new lu) 

	movsd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movsd xmm4, [rsp + nb233nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb233nf_half] ;# rinv 
	movsd  xmm7, xmm4	;# rinvO in xmm7 
	
	movsd xmm4, [rsp + nb233nf_rsqO]
	movapd xmm0, xmm7
	;# LJ table interaction.
	mulsd xmm4, xmm7	;# xmm4=r 
	mulsd xmm4, [rsp + nb233nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3

	mov  rsi, [rbp + nb233nf_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
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

	movsd xmm4, [rsp + nb233nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb233nf_Vvdwtot]
	movsd [rsp + nb233nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

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
	
	movsd xmm4, [rsp + nb233nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb233nf_Vvdwtot]
	movsd [rsp + nb233nf_Vvdwtot], xmm5

	;# H1 interactions
	movsd  xmm6, [rsp + nb233nf_rinvH1] 
	movsd  xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm6
	movsd  xmm0, [rsp + nb233nf_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subsd   xmm6, [rsp + nb233nf_crf]
	mulsd   xmm6, [rsp + nb233nf_qqH] ;# vcoul 
	addsd  xmm6, [rsp + nb233nf_vctot]
	movsd [rsp + nb233nf_vctot], xmm6

	;# H2 interactions 
	movsd  xmm5, [rsp + nb233nf_rinvH2] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [rsp + nb233nf_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [rsp + nb233nf_crf]
	mulsd   xmm5, [rsp + nb233nf_qqH] ;# vcoul 
	addsd  xmm5, [rsp + nb233nf_vctot]
	movsd [rsp + nb233nf_vctot], xmm5

	;# M interactions 
	movsd  xmm5, [rsp + nb233nf_rinvM] 
	movsd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movsd  xmm7, xmm5
	movsd  xmm0, [rsp + nb233nf_krsqM]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [rsp + nb233nf_crf]
	mulsd   xmm5, [rsp + nb233nf_qqM] ;# vcoul 
	addsd  xmm5, [rsp + nb233nf_vctot]
	movsd [rsp + nb233nf_vctot], xmm5

.nb233nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb233nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb233nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb233nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb233nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb233nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb233nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb233nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb233nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb233nf_n], esi
        jmp .nb233nf_outer
.nb233nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb233nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb233nf_end
        ;# non-zero, do one more workunit
        jmp   .nb233nf_threadloop
.nb233nf_end:
	mov eax, [rsp + nb233nf_nouter]
	mov ebx, [rsp + nb233nf_ninner]
	mov rcx, [rbp + nb233nf_outeriter]
	mov rdx, [rbp + nb233nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 608
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret
