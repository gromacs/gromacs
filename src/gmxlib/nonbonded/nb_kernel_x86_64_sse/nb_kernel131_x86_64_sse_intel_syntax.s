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


.globl nb_kernel131_x86_64_sse
.globl _nb_kernel131_x86_64_sse
nb_kernel131_x86_64_sse:	
_nb_kernel131_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb131_fshift,           16
.equiv          nb131_gid,              24
.equiv          nb131_pos,              32
.equiv          nb131_faction,          40
.equiv          nb131_charge,           48
.equiv          nb131_p_facel,          56
.equiv          nb131_argkrf,           64
.equiv          nb131_argcrf,           72
.equiv          nb131_Vc,               80
.equiv          nb131_type,             88
.equiv          nb131_p_ntype,          96
.equiv          nb131_vdwparam,         104
.equiv          nb131_Vvdw,             112
.equiv          nb131_p_tabscale,       120
.equiv          nb131_VFtab,            128
.equiv          nb131_invsqrta,         136
.equiv          nb131_dvda,             144
.equiv          nb131_p_gbtabscale,     152
.equiv          nb131_GBtab,            160
.equiv          nb131_p_nthreads,       168
.equiv          nb131_count,            176
.equiv          nb131_mtx,              184
.equiv          nb131_outeriter,        192
.equiv          nb131_inneriter,        200
.equiv          nb131_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb131_ixO,              0
.equiv          nb131_iyO,              16
.equiv          nb131_izO,              32
.equiv          nb131_ixH1,             48
.equiv          nb131_iyH1,             64
.equiv          nb131_izH1,             80
.equiv          nb131_ixH2,             96
.equiv          nb131_iyH2,             112
.equiv          nb131_izH2,             128
.equiv          nb131_iqO,              144
.equiv          nb131_iqH,              160
.equiv          nb131_dxO,              176
.equiv          nb131_dyO,              192
.equiv          nb131_dzO,              208
.equiv          nb131_dxH1,             224
.equiv          nb131_dyH1,             240
.equiv          nb131_dzH1,             256
.equiv          nb131_dxH2,             272
.equiv          nb131_dyH2,             288
.equiv          nb131_dzH2,             304
.equiv          nb131_qqO,              320
.equiv          nb131_qqH,              336
.equiv          nb131_c6,               352
.equiv          nb131_c12,              368
.equiv          nb131_tsc,              384
.equiv          nb131_fstmp,            400
.equiv          nb131_vctot,            416
.equiv          nb131_Vvdwtot,          432
.equiv          nb131_fixO,             448
.equiv          nb131_fiyO,             464
.equiv          nb131_fizO,             480
.equiv          nb131_fixH1,            496
.equiv          nb131_fiyH1,            512
.equiv          nb131_fizH1,            528
.equiv          nb131_fixH2,            544
.equiv          nb131_fiyH2,            560
.equiv          nb131_fizH2,            576
.equiv          nb131_fjx,              592
.equiv          nb131_fjy,              608
.equiv          nb131_fjz,              624
.equiv          nb131_half,             640
.equiv          nb131_three,            656
.equiv          nb131_two,              672
.equiv          nb131_krf,              688
.equiv          nb131_crf,              704
.equiv          nb131_rsqO,             720
.equiv          nb131_rsqH1,            736
.equiv          nb131_rsqH2,            752
.equiv          nb131_rinvO,            768
.equiv          nb131_rinvH1,           784
.equiv          nb131_rinvH2,           800
.equiv          nb131_facel,            816
.equiv          nb131_iinr,             824
.equiv          nb131_jindex,           832
.equiv          nb131_jjnr,             840
.equiv          nb131_shift,            848
.equiv          nb131_shiftvec,         856
.equiv          nb131_innerjjnr,        864
.equiv          nb131_nri,              872
.equiv          nb131_is3,              876
.equiv          nb131_ii3,              880
.equiv          nb131_ntia,             884
.equiv          nb131_innerk,           888
.equiv          nb131_n,                892
.equiv          nb131_nn1,              896
.equiv          nb131_nouter,           900
.equiv          nb131_ninner,           904

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
    
	emms
	sub rsp, 912		;# local variable stack space (n*16+8)
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

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb131_nouter], eax
	mov [rsp + nb131_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb131_nri], edi
	mov [rsp + nb131_iinr], rsi
	mov [rsp + nb131_jindex], rdx
	mov [rsp + nb131_jjnr], rcx
	mov [rsp + nb131_shift], r8
	mov [rsp + nb131_shiftvec], r9
	mov rsi, [rbp + nb131_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb131_facel], xmm0

	mov rax, [rbp + nb131_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb131_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb131_half], eax
	movss xmm1, [rsp + nb131_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb131_half],  xmm1
	movaps [rsp + nb131_two],  xmm2
	movaps [rsp + nb131_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb131_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb131_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb131_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb131_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb131_iqO], xmm3
	movaps [rsp + nb131_iqH], xmm4
	
	mov   rdx, [rbp + nb131_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb131_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb131_ntia], ecx		

.nb131_threadloop:
        mov   rsi, [rbp + nb131_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb131_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb131_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb131_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb131_n], eax
        mov [rsp + nb131_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb131_outerstart
        jmp .nb131_end

.nb131_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb131_nouter]
	mov [rsp + nb131_nouter], ebx

.nb131_outer:
	mov   rax, [rsp + nb131_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax +rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb131_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb131_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb131_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb131_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb131_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb131_ixO], xmm3
	movaps [rsp + nb131_iyO], xmm4
	movaps [rsp + nb131_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb131_ixH1], xmm0
	movaps [rsp + nb131_iyH1], xmm1
	movaps [rsp + nb131_izH1], xmm2
	movaps [rsp + nb131_ixH2], xmm3
	movaps [rsp + nb131_iyH2], xmm4
	movaps [rsp + nb131_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb131_vctot], xmm4
	movaps [rsp + nb131_Vvdwtot], xmm4
	movaps [rsp + nb131_fixO], xmm4
	movaps [rsp + nb131_fiyO], xmm4
	movaps [rsp + nb131_fizO], xmm4
	movaps [rsp + nb131_fixH1], xmm4
	movaps [rsp + nb131_fiyH1], xmm4
	movaps [rsp + nb131_fizH1], xmm4
	movaps [rsp + nb131_fixH2], xmm4
	movaps [rsp + nb131_fiyH2], xmm4
	movaps [rsp + nb131_fizH2], xmm4
	
	mov   rax, [rsp + nb131_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb131_pos]
	mov   rdi, [rbp + nb131_faction]	
	mov   rax, [rsp + nb131_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb131_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb131_ninner]
	mov   [rsp + nb131_ninner], ecx
	add   edx, 0
	mov   [rsp + nb131_innerk], edx    ;# number of innerloop atoms 
	jge   .nb131_unroll_loop
	jmp   .nb131_odd_inner
.nb131_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb131_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb131_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb131_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb131_iqO]
	mulps  xmm4, [rsp + nb131_iqH]

	movaps  [rsp + nb131_qqO], xmm3
	movaps  [rsp + nb131_qqH], xmm4
	
	mov rsi, [rbp + nb131_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov r10d, [rsi + rcx*4]
	mov r11d, [rsi + rdx*4]
	mov rsi, [rbp + nb131_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	
	mov edi, [rsp + nb131_ntia]
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

	movlps xmm6, [rsi + r8*4]
	movlps xmm7, [rsi + r10*4]
	movhps xmm6, [rsi + r9*4]
	movhps xmm7, [rsi + r11*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101

	movaps [rsp + nb131_c6], xmm4
	movaps [rsp + nb131_c12], xmm6

	mov rsi, [rbp + nb131_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four j coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2

    subps xmm0, [rsp + nb131_ixO]
    subps xmm1, [rsp + nb131_iyO]
    subps xmm2, [rsp + nb131_izO]
    subps xmm3, [rsp + nb131_ixH1]
    subps xmm4, [rsp + nb131_iyH1]
    subps xmm5, [rsp + nb131_izH1]
    subps xmm6, [rsp + nb131_ixH2]
    subps xmm7, [rsp + nb131_iyH2]
    subps xmm8, [rsp + nb131_izH2]
    
	movaps [rsp + nb131_dxO], xmm0
	movaps [rsp + nb131_dyO], xmm1
	movaps [rsp + nb131_dzO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb131_dxH1], xmm3
	movaps [rsp + nb131_dyH1], xmm4
	movaps [rsp + nb131_dzH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb131_dxH2], xmm6
	movaps [rsp + nb131_dyH2], xmm7
	movaps [rsp + nb131_dzH2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb131_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm2, [rsp + nb131_half]
	mulps   xmm9, xmm2  ;# rinvO 
	mulps   xmm10, xmm2 ;# rinvH1
    mulps   xmm11, xmm2 ;# rinvH2
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11
    movaps [rsp + nb131_rsqO], xmm0
    movaps [rsp + nb131_rsqH1], xmm3
    movaps [rsp + nb131_rsqH2], xmm6
    movaps [rsp + nb131_rinvO], xmm9
    movaps [rsp + nb131_rinvH1], xmm10
    movaps [rsp + nb131_rinvH2], xmm11

    ;# table LJ interaction
    mulps  xmm0, xmm9
    mulps  xmm0, [rsp + nb131_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0

    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
         
    ;# multiply by 8
    pslld   xmm1, 3

    ;# move to integer registers
    movhlps xmm13, xmm1
    movd    r8d, xmm1
    movd    r10d, xmm13
    shufps  xmm1, xmm1, 1
    shufps  xmm13, xmm13, 1
    movd    r9d, xmm1
    movd    r11d, xmm13
    
    ;# calculate eps
    subps     xmm0, xmm2
    mov  rsi, [rbp + nb131_VFtab]
            
    movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r8*4 + 16]

	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r10*4 + 16]

	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r9*4 + 16]

	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r11*4 + 16]

    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101

	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r8*4 + 24]
    
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r10*4 + 24]

	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r9*4 + 24]
    
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r11*4 + 24]

    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulps  xmm7, xmm0    ;# Heps
    mulps  xmm11, xmm0 
    mulps  xmm6, xmm0   ;# Geps
    mulps  xmm10, xmm0 
    mulps  xmm7, xmm0   ;# Heps2
    mulps  xmm11, xmm0 
    addps  xmm5, xmm6  ;# F+Geps
    addps  xmm9, xmm10 
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm9, xmm11 
    addps  xmm7, xmm7    ;# 2*Heps2
    addps  xmm11, xmm11
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm11, xmm10
    
    addps  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm11, xmm9
    mulps  xmm5, xmm0  ;# eps*Fp
    mulps  xmm9, xmm0
    movaps xmm12, [rsp + nb131_c6]
    movaps xmm13, [rsp + nb131_c12]
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

    mulps  xmm5, xmm12  ;# VV*c6 = vnb6
    mulps  xmm9, xmm13  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb131_Vvdwtot]
    movaps [rsp + nb131_Vvdwtot], xmm5
        
    mulps  xmm7, xmm12   ;# FF*c6 = fnb6
    mulps  xmm11, xmm13   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    mulps  xmm7, [rsp + nb131_tsc]
     
    movaps xmm9, [rsp + nb131_rinvO]
    movaps xmm10, [rsp + nb131_rinvH1]
    movaps xmm11, [rsp + nb131_rinvH2]
    movaps xmm0, xmm9   ;# rinv
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    
    mulps  xmm10, xmm10  ;# rinvsq
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb131_qqO] 
    mulps  xmm1, [rsp + nb131_qqH] 
    mulps  xmm2, [rsp + nb131_qqH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    subps  xmm9, xmm7
    mulps  xmm9, [rsp + nb131_rinvO]
    
    addps xmm0, [rsp + nb131_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb131_vctot], xmm0

	;# move j  forces to local temp variables 
	mov rdi, [rbp + nb131_faction]
    movlps xmm0, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm1, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm0, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm1, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm2, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm3, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm4, [rdi + rbx*4 + 8] ;# jzb
    movss  xmm5, [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm2, xmm4 ;# jza  -  jzb  -
    movlhps xmm3, xmm5 ;# jzc  -  jzd -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm0: jxa jya jxb jyb 
    ;# xmm1: jxc jyc jxd jyd
    ;# xmm2: jza jzb jzc jzd

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb131_dxO]
	mulps xmm8, [rsp + nb131_dyO]
	mulps xmm9, [rsp + nb131_dzO]
	mulps xmm10, [rsp + nb131_dxH1]
	mulps xmm11, [rsp + nb131_dyH1]
	mulps xmm12, [rsp + nb131_dzH1]
	mulps xmm13, [rsp + nb131_dxH2]
	mulps xmm14, [rsp + nb131_dyH2]
	mulps xmm15, [rsp + nb131_dzH2]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb131_fixO]
    addps xmm8, [rsp + nb131_fiyO]
    addps xmm9, [rsp + nb131_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb131_fixH1]
    addps xmm11, [rsp + nb131_fiyH1]
    addps xmm12, [rsp + nb131_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb131_fixH2]
    addps xmm14, [rsp + nb131_fiyH2]
    addps xmm15, [rsp + nb131_fizH2]

    movaps [rsp + nb131_fixO], xmm7
    movaps [rsp + nb131_fiyO], xmm8
    movaps [rsp + nb131_fizO], xmm9
    movaps [rsp + nb131_fixH1], xmm10
    movaps [rsp + nb131_fiyH1], xmm11
    movaps [rsp + nb131_fizH1], xmm12
    movaps [rsp + nb131_fixH2], xmm13
    movaps [rsp + nb131_fiyH2], xmm14
    movaps [rsp + nb131_fizH2], xmm15
    
    ;# xmm0 = fjx
    ;# xmm1 = fjy
    ;# xmm2 = fjz
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 
    
    movlps [rdi + rax*4], xmm0
    movhps [rdi + rbx*4], xmm0
    movlps [rdi + rcx*4], xmm1
    movhps [rdi + rdx*4], xmm1
    movss  [rdi + rax*4 + 8], xmm2
    movss  [rdi + rcx*4 + 8], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 8], xmm2
    movss  [rdi + rdx*4 + 8], xmm3

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb131_innerk],  4
	jl    .nb131_odd_inner
	jmp   .nb131_unroll_loop
.nb131_odd_inner:	
	add dword ptr [rsp + nb131_innerk],  4
	jnz   .nb131_odd_loop
	jmp   .nb131_updateouterdata
.nb131_odd_loop:
	mov   rdx, [rsp + nb131_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb131_innerjjnr],  4

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb131_iqO]
	mov rsi, [rbp + nb131_charge] 
	movhps xmm4, [rsp + nb131_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb131_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb131_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb131_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb131_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [rsp + nb131_c6], xmm6
	movaps [rsp + nb131_c12], xmm7

	mov rsi, [rbp + nb131_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	movss xmm0, [rsp + nb131_ixO]
	movss xmm1, [rsp + nb131_iyO]
	movss xmm2, [rsp + nb131_izO]
	
	movlps xmm6, [rsp + nb131_ixH1]
	movlps xmm7, [rsp + nb131_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm0, xmm6
	movlps xmm6, [rsp + nb131_iyH1]
	movlps xmm7, [rsp + nb131_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm1, xmm6
	movlps xmm6, [rsp + nb131_izH1]
	movlps xmm7, [rsp + nb131_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm2, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [rsp + nb131_dxO], xmm3
	movaps [rsp + nb131_dyO], xmm4
	movaps [rsp + nb131_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb131_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb131_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000
		
	;# rsq still in xmm4, rinv in xmm0.
	mulps  xmm4, xmm0
	mulps  xmm4, [rsp + nb131_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3	

	mov  rsi, [rbp + nb131_VFtab]
	movd r12d, mm6
	
	;# dispersion 
	movlps xmm5, [rsi + r12*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + r12*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [rsp + nb131_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [rsp + nb131_c6]
	mulss  xmm7, xmm4	 ;# fijD 
	mulss  xmm5, xmm4	 ;# Vvdw6 
	mulss  xmm7, [rsp + nb131_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [rsp + nb131_Vvdwtot]
	movss [rsp + nb131_fstmp], xmm7
	movss [rsp + nb131_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + r12*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + r12*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [rsp + nb131_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [rsp + nb131_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	mulss  xmm7, [rsp + nb131_tsc]
	addss  xmm7, [rsp + nb131_fstmp]
	movss [rsp + nb131_fstmp], xmm7

	addss  xmm5, [rsp + nb131_Vvdwtot]
	movss [rsp + nb131_Vvdwtot], xmm5

	;# coulomb in analytical form. xmm0=rinv
	movaps xmm2, xmm0               ;# rinv
	mulps  xmm2, [rsp + nb131_qqO]  ;# vcoulO
	movaps xmm3, xmm2               ;# vcoulO
	
	addps  xmm2, [rsp + nb131_vctot]
	movaps [rsp + nb131_vctot], xmm2
	
	mulps  xmm3, xmm0
	subss  xmm3, [rsp + nb131_fstmp]
	mulps  xmm3, xmm0
	
	movaps xmm0, [rsp + nb131_dxO]
	movaps xmm1, [rsp + nb131_dyO]
	movaps xmm2, [rsp + nb131_dzO]

	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb131_fixO]	
	movss  xmm4, [rsp + nb131_fiyO]	
	movss  xmm5, [rsp + nb131_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb131_fixO], xmm3	
	movss  [rsp + nb131_fiyO], xmm4	
	movss  [rsp + nb131_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [rsp + nb131_fixH1]
	addss  xmm4, [rsp + nb131_fiyH1]
	addss  xmm5, [rsp + nb131_fizH1]
	movss  [rsp + nb131_fixH1], xmm3	
	movss  [rsp + nb131_fiyH1], xmm4	
	movss  [rsp + nb131_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb131_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [rsp + nb131_fixH2]
	addss  xmm4, [rsp + nb131_fiyH2]
	addss  xmm5, [rsp + nb131_fizH2]
	movss  [rsp + nb131_fixH2], xmm3	
	movss  [rsp + nb131_fiyH2], xmm4	
	movss  [rsp + nb131_fizH2], xmm5	;# updated the H2 force 

	;# the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [rdi + rax*4]
	movss  xmm7, [rdi + rax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	;# x,y sum in xmm0 

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1 
	addss    xmm2, xmm1    ;# z sum in xmm2 
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb131_innerk]
	jz    .nb131_updateouterdata
	jmp   .nb131_odd_loop
.nb131_updateouterdata:
	mov   ecx, [rsp + nb131_ii3]
	mov   rdi, [rbp + nb131_faction]
	mov   rsi, [rbp + nb131_fshift]
	mov   edx, [rsp + nb131_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb131_fixO]
	movaps xmm1, [rsp + nb131_fiyO]
	movaps xmm2, [rsp + nb131_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4]
	movss  xmm4, [rdi + rcx*4 + 4]
	movss  xmm5, [rdi + rcx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb131_fixH1]
	movaps xmm1, [rsp + nb131_fiyH1]
	movaps xmm2, [rsp + nb131_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb131_fixH2]
	movaps xmm1, [rsp + nb131_fiyH2]
	movaps xmm2, [rsp + nb131_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb131_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb131_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb131_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb131_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb131_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb131_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb131_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb131_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb131_n], esi
        jmp .nb131_outer
.nb131_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb131_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb131_end
        ;# non-zero, do one more workunit
        jmp   .nb131_threadloop
.nb131_end:
	mov eax, [rsp + nb131_nouter]
	mov ebx, [rsp + nb131_ninner]
	mov rcx, [rbp + nb131_outeriter]
	mov rdx, [rbp + nb131_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 912
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



.globl nb_kernel131nf_x86_64_sse
.globl _nb_kernel131nf_x86_64_sse
nb_kernel131nf_x86_64_sse:	
_nb_kernel131nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb131nf_fshift,         16
.equiv          nb131nf_gid,            24
.equiv          nb131nf_pos,            32
.equiv          nb131nf_faction,        40
.equiv          nb131nf_charge,         48
.equiv          nb131nf_p_facel,        56
.equiv          nb131nf_argkrf,         64
.equiv          nb131nf_argcrf,         72
.equiv          nb131nf_Vc,             80
.equiv          nb131nf_type,           88
.equiv          nb131nf_p_ntype,        96
.equiv          nb131nf_vdwparam,       104
.equiv          nb131nf_Vvdw,           112
.equiv          nb131nf_p_tabscale,     120
.equiv          nb131nf_VFtab,          128
.equiv          nb131nf_invsqrta,       136
.equiv          nb131nf_dvda,           144
.equiv          nb131nf_p_gbtabscale,   152
.equiv          nb131nf_GBtab,          160
.equiv          nb131nf_p_nthreads,     168
.equiv          nb131nf_count,          176
.equiv          nb131nf_mtx,            184
.equiv          nb131nf_outeriter,      192
.equiv          nb131nf_inneriter,      200
.equiv          nb131nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb131nf_ixO,            0
.equiv          nb131nf_iyO,            16
.equiv          nb131nf_izO,            32
.equiv          nb131nf_ixH1,           48
.equiv          nb131nf_iyH1,           64
.equiv          nb131nf_izH1,           80
.equiv          nb131nf_ixH2,           96
.equiv          nb131nf_iyH2,           112
.equiv          nb131nf_izH2,           128
.equiv          nb131nf_iqO,            144
.equiv          nb131nf_iqH,            160
.equiv          nb131nf_qqO,            176
.equiv          nb131nf_qqH,            192
.equiv          nb131nf_c6,             208
.equiv          nb131nf_c12,            224
.equiv          nb131nf_vctot,          240
.equiv          nb131nf_Vvdwtot,        256
.equiv          nb131nf_half,           272
.equiv          nb131nf_three,          288
.equiv          nb131nf_rinvO,          304
.equiv          nb131nf_rinvH1,         320
.equiv          nb131nf_rinvH2,	        336
.equiv          nb131nf_krsqO,          352
.equiv          nb131nf_krsqH1,         368
.equiv          nb131nf_krsqH2,         384
.equiv          nb131nf_tsc,            400
.equiv          nb131nf_facel,          416
.equiv          nb131nf_iinr,           424
.equiv          nb131nf_jindex,         432
.equiv          nb131nf_jjnr,           440
.equiv          nb131nf_shift,          448
.equiv          nb131nf_shiftvec,       456
.equiv          nb131nf_innerjjnr,      464
.equiv          nb131nf_nri,            472
.equiv          nb131nf_ntia,           476
.equiv          nb131nf_is3,            480
.equiv          nb131nf_ii3,            484
.equiv          nb131nf_innerk,         488
.equiv          nb131nf_n,              492
.equiv          nb131nf_nn1,            496
.equiv          nb131nf_nouter,         500
.equiv          nb131nf_ninner,         504

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
    
	emms
	sub rsp, 512		;# local variable stack space (n*16+8)
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

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb131nf_nouter], eax
	mov [rsp + nb131nf_ninner], eax
	
	mov edi, [rdi]
	mov [rsp + nb131nf_nri], edi
	mov [rsp + nb131nf_iinr], rsi
	mov [rsp + nb131nf_jindex], rdx
	mov [rsp + nb131nf_jjnr], rcx
	mov [rsp + nb131nf_shift], r8
	mov [rsp + nb131nf_shiftvec], r9
	mov rsi, [rbp + nb131nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb131nf_facel], xmm0

	mov rax, [rbp + nb131nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb131nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb131nf_half], eax
	movss xmm1, [rsp + nb131nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb131nf_half],  xmm1
	movaps [rsp + nb131nf_three],  xmm3	

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb131nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb131nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb131nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb131nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb131nf_iqO], xmm3
	movaps [rsp + nb131nf_iqH], xmm4
	
	mov   rdx, [rbp + nb131nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb131nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb131nf_ntia], ecx		

.nb131nf_threadloop:
        mov   rsi, [rbp + nb131nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb131nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb131nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb131nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb131nf_n], eax
        mov [rsp + nb131nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb131nf_outerstart
        jmp .nb131nf_end

.nb131nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb131nf_nouter]
	mov [rsp + nb131nf_nouter], ebx

.nb131nf_outer:
	mov   rax, [rsp + nb131nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax +rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb131nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb131nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb131nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb131nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb131nf_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb131nf_ixO], xmm3
	movaps [rsp + nb131nf_iyO], xmm4
	movaps [rsp + nb131nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb131nf_ixH1], xmm0
	movaps [rsp + nb131nf_iyH1], xmm1
	movaps [rsp + nb131nf_izH1], xmm2
	movaps [rsp + nb131nf_ixH2], xmm3
	movaps [rsp + nb131nf_iyH2], xmm4
	movaps [rsp + nb131nf_izH2], xmm5
	
	;# clear vctot
	xorps xmm4, xmm4
	movaps [rsp + nb131nf_vctot], xmm4
	movaps [rsp + nb131nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb131nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb131nf_pos]
	mov   rax, [rsp + nb131nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb131nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb131nf_ninner]
	mov   [rsp + nb131nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb131nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb131nf_unroll_loop
	jmp   .nb131nf_odd_inner
.nb131nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb131nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb131nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb131nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb131nf_iqO]
	mulps  xmm4, [rsp + nb131nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb131nf_qqO], xmm3
	movaps  [rsp + nb131nf_qqH], xmm4
	
	mov rsi, [rbp + nb131nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb131nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb131nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb131nf_c6], xmm4
	movaps [rsp + nb131nf_c12], xmm6

	mov rsi, [rbp + nb131nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb131nf_ixO]
	movaps xmm5, [rsp + nb131nf_iyO]
	movaps xmm6, [rsp + nb131nf_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb131nf_ixH1]
	movaps xmm5, [rsp + nb131nf_iyH1]
	movaps xmm6, [rsp + nb131nf_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [rsp + nb131nf_ixH2]
	movaps xmm4, [rsp + nb131nf_iyH2]
	movaps xmm5, [rsp + nb131nf_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	;# start with rsqO (still in xmm7) - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb131nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb131nf_half]
	movaps  [rsp + nb131nf_rinvO], xmm4	;# rinvO in xmm0, rsqO in xmm7
	
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb131nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb131nf_half]
	movaps  [rsp + nb131nf_rinvH1], xmm4	;# rinvH1 in xmm6 

	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb131nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb131nf_half]
	movaps  [rsp + nb131nf_rinvH2], xmm4	;# rinvH2 in xmm5 
	
	
	;# do O table interactions - rsqO in xmm7.
	mulps xmm7, [rsp + nb131nf_rinvO]
	mulps xmm7, [rsp + nb131nf_tsc] ;# rtab
	
	movhlps xmm5, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps  xmm7, xmm6	
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	mov  rsi, [rbp + nb131nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb131nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 

	addps  xmm5, [rsp + nb131nf_Vvdwtot]
	movaps [rsp + nb131nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movlps xmm7, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb131nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 

	addps  xmm5, [rsp + nb131nf_Vvdwtot]
	movaps [rsp + nb131nf_Vvdwtot], xmm5
	
	movaps xmm2, [rsp + nb131nf_rinvO]
	mulps  xmm2, [rsp + nb131nf_qqO]

	;# H1 & H2 interactions 
	movaps xmm3, [rsp + nb131nf_rinvH1]
	addps  xmm3, [rsp + nb131nf_rinvH2]
	
	mulps  xmm3, [rsp + nb131nf_qqH]
	
	addps  xmm2, xmm3
	addps  xmm2, [rsp + nb131nf_vctot]
	movaps [rsp + nb131nf_vctot], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb131nf_innerk],  4
	jl    .nb131nf_odd_inner
	jmp   .nb131nf_unroll_loop
.nb131nf_odd_inner:	
	add dword ptr [rsp + nb131nf_innerk],  4
	jnz   .nb131nf_odd_loop
	jmp   .nb131nf_updateouterdata
.nb131nf_odd_loop:
	mov   rdx, [rsp + nb131nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb131nf_innerjjnr],  4

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb131nf_iqO]
	mov rsi, [rbp + nb131nf_charge] 
	movhps xmm4, [rsp + nb131nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb131nf_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb131nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb131nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb131nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [rsp + nb131nf_c6], xmm6
	movaps [rsp + nb131nf_c12], xmm7

	mov rsi, [rbp + nb131nf_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb131nf_ixO]
	movss xmm4, [rsp + nb131nf_iyO]
	movss xmm5, [rsp + nb131nf_izO]
		
	movlps xmm6, [rsp + nb131nf_ixH1]
	movlps xmm7, [rsp + nb131nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb131nf_iyH1]
	movlps xmm7, [rsp + nb131nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb131nf_izH1]
	movlps xmm7, [rsp + nb131nf_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb131nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb131nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000
		
	;# rsq still in xmm4, rinv in xmm0.
	mulps  xmm4, xmm0
	mulps  xmm4, [rsp + nb131nf_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	mov  rsi, [rbp + nb131nf_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [rsp + nb131nf_c6]
	mulss  xmm5, xmm4	 ;# Vvdw6 

	addss  xmm5, [rsp + nb131nf_Vvdwtot]
	movss [rsp + nb131nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [rsp + nb131nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 
	addss  xmm5, [rsp + nb131nf_Vvdwtot]
	movss [rsp + nb131nf_Vvdwtot], xmm5

	;# xmm0=rinv
	mulps  xmm0, [rsp + nb131nf_qqO]	;# xmm0=vcoul 
	
	addps  xmm0, [rsp + nb131nf_vctot]
	movaps [rsp + nb131nf_vctot], xmm0
	
	dec dword ptr [rsp + nb131nf_innerk]
	jz    .nb131nf_updateouterdata
	jmp   .nb131nf_odd_loop
.nb131nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb131nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb131nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb131nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb131nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb131nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb131nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb131nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb131nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb131nf_n], esi
        jmp .nb131nf_outer
.nb131nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb131nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb131nf_end
        ;# non-zero, do one more workunit
        jmp   .nb131nf_threadloop
.nb131nf_end:
	mov eax, [rsp + nb131nf_nouter]
	mov ebx, [rsp + nb131nf_ninner]
	mov rcx, [rbp + nb131nf_outeriter]
	mov rdx, [rbp + nb131nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 512
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
