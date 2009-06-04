##
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##




.globl nb_kernel430_x86_64_sse2
.globl _nb_kernel430_x86_64_sse2
nb_kernel430_x86_64_sse2:       
_nb_kernel430_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb430_fshift, 16
.set nb430_gid, 24
.set nb430_pos, 32
.set nb430_faction, 40
.set nb430_charge, 48
.set nb430_p_facel, 56
.set nb430_argkrf, 64
.set nb430_argcrf, 72
.set nb430_Vc, 80
.set nb430_type, 88
.set nb430_p_ntype, 96
.set nb430_vdwparam, 104
.set nb430_Vvdw, 112
.set nb430_p_tabscale, 120
.set nb430_VFtab, 128
.set nb430_invsqrta, 136
.set nb430_dvda, 144
.set nb430_p_gbtabscale, 152
.set nb430_GBtab, 160
.set nb430_p_nthreads, 168
.set nb430_count, 176
.set nb430_mtx, 184
.set nb430_outeriter, 192
.set nb430_inneriter, 200
.set nb430_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb430_ix, 0
.set nb430_iy, 16
.set nb430_iz, 32
.set nb430_iq, 48
.set nb430_dx, 64
.set nb430_dy, 80
.set nb430_dz, 96
.set nb430_eps, 112
.set nb430_gbtsc, 128
.set nb430_tsc, 144
.set nb430_qq, 160
.set nb430_c6, 176
.set nb430_c12, 192
.set nb430_epsgb, 208
.set nb430_vctot, 224
.set nb430_Vvdwtot, 240
.set nb430_fix, 256
.set nb430_fiy, 272
.set nb430_fiz, 288
.set nb430_half, 304
.set nb430_three, 320
.set nb430_r, 336
.set nb430_isai, 352
.set nb430_isaprod, 368
.set nb430_dvdasum, 384
.set nb430_gbscale, 400
.set nb430_rinv, 416
.set nb430_nri, 432
.set nb430_iinr, 440
.set nb430_jindex, 448
.set nb430_jjnr, 456
.set nb430_shift, 464
.set nb430_shiftvec, 472
.set nb430_facel, 480
.set nb430_innerjjnr, 488
.set nb430_ii, 496
.set nb430_is3, 500
.set nb430_ii3, 504
.set nb430_ntia, 508
.set nb430_innerk, 512
.set nb430_n, 516
.set nb430_nn1, 520
.set nb430_ntype, 524
.set nb430_nouter, 528
.set nb430_ninner, 532

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $536,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb430_nouter(%rsp)
        movl %eax,nb430_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb430_nri(%rsp)
        movq %rsi,nb430_iinr(%rsp)
        movq %rdx,nb430_jindex(%rsp)
        movq %rcx,nb430_jjnr(%rsp)
        movq %r8,nb430_shift(%rsp)
        movq %r9,nb430_shiftvec(%rsp)
        movq nb430_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb430_ntype(%rsp)
        movq nb430_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb430_facel(%rsp)

        movq nb430_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb430_tsc(%rsp)

        movq nb430_p_gbtabscale(%rbp),%rbx
        movsd (%rbx),%xmm4
        shufpd $0,%xmm4,%xmm4
        movapd %xmm4,nb430_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb430_half(%rsp)
        movl %ebx,nb430_half+4(%rsp)
        movsd nb430_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb430_half(%rsp)
        movapd %xmm3,nb430_three(%rsp)

_nb_kernel430_x86_64_sse2.nb430_threadloop: 
        movq  nb430_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel430_x86_64_sse2.nb430_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430_x86_64_sse2.nb430_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb430_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb430_n(%rsp)
        movl %ebx,nb430_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel430_x86_64_sse2.nb430_outerstart
        jmp _nb_kernel430_x86_64_sse2.nb430_end

_nb_kernel430_x86_64_sse2.nb430_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430_nouter(%rsp),%ebx
        movl %ebx,nb430_nouter(%rsp)

_nb_kernel430_x86_64_sse2.nb430_outer: 
        movq  nb430_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb430_is3(%rsp)      ## store is3 

        movq  nb430_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb430_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb430_ii(%rsp)

        movq  nb430_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb430_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb430_invsqrta(%rbp),%rdx         ## load invsqrta[ii]
        movsd (%rdx,%rbx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movq  nb430_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb430_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb430_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb430_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb430_iq(%rsp)
        movapd %xmm4,nb430_isai(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb430_ix(%rsp)
        movapd %xmm1,nb430_iy(%rsp)
        movapd %xmm2,nb430_iz(%rsp)

        movl  %ebx,nb430_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb430_vctot(%rsp)
        movapd %xmm4,nb430_Vvdwtot(%rsp)
        movapd %xmm4,nb430_dvdasum(%rsp)
        movapd %xmm4,nb430_fix(%rsp)
        movapd %xmm4,nb430_fiy(%rsp)
        movapd %xmm4,nb430_fiz(%rsp)

        movq  nb430_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb430_pos(%rbp),%rsi
        movq  nb430_faction(%rbp),%rdi
        movq  nb430_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb430_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb430_ninner(%rsp),%ecx
        movl  %ecx,nb430_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb430_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel430_x86_64_sse2.nb430_unroll_loop
        jmp   _nb_kernel430_x86_64_sse2.nb430_checksingle
_nb_kernel430_x86_64_sse2.nb430_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb430_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb430_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 


        movq nb430_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%rax,%rax,2),%r10     ## j3 
        lea  (%rbx,%rbx,2),%r11

        ## move two coordinates to xmm4-xmm6 
        movlpd (%rsi,%r10,8),%xmm4
        movlpd 8(%rsi,%r10,8),%xmm5
        movlpd 16(%rsi,%r10,8),%xmm6
        movhpd (%rsi,%r11,8),%xmm4
        movhpd 8(%rsi,%r11,8),%xmm5
        movhpd 16(%rsi,%r11,8),%xmm6

        ## calc dr 
        subpd nb430_ix(%rsp),%xmm4
        subpd nb430_iy(%rsp),%xmm5
        subpd nb430_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb430_dx(%rsp)
        movapd %xmm5,nb430_dy(%rsp)
        movapd %xmm6,nb430_dz(%rsp)

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        ## load isaj
        movq nb430_invsqrta(%rbp),%rsi
        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        mulpd  nb430_isai(%rsp),%xmm3
        movapd %xmm3,nb430_isaprod(%rsp)
        movapd %xmm3,%xmm6
        mulpd nb430_gbtsc(%rsp),%xmm3
        movapd %xmm3,nb430_gbscale(%rsp)

    ##invsqrt
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb430_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

    mulpd  nb430_iq(%rsp),%xmm6
        movq nb430_charge(%rbp),%rsi     ## base of charge[] 
        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        mulpd  %xmm6,%xmm3
        movapd %xmm3,nb430_qq(%rsp)

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb430_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb430_r(%rsp)
        movapd %xmm0,nb430_rinv(%rsp)

        movq nb430_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        shll %r8d
        shll %r9d
        movl nb430_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

    movapd %xmm4,%xmm8 ## r
        mulpd nb430_gbscale(%rsp),%xmm4
        mulpd nb430_tsc(%rsp),%xmm8

    ## truncate and convert to integers
    cvttpd2pi %xmm4,%mm0 ## gb
    cvttpd2pi %xmm8,%mm1 ## lj

    ## convert back to float
    cvtpi2pd  %mm0,%xmm6  ## gb
    cvtpi2pd  %mm1,%xmm10 ## lj

    ## multiply by 4 and 8, respectively
    pslld   $2,%mm0  ## gb
    pslld   $3,%mm1  ## lj

    ## move to integer registers
    movd    %mm0,%r12d      ## gb
    movd    %mm1,%r14d     ## lj
        psrlq $32,%mm0
        psrlq $32,%mm1
    movd    %mm0,%r13d     ## gb
    movd    %mm1,%r15d    ## lj
    ## GB indices: r10-11   LJ indices: r12-r13

    ## calculate eps
    subpd     %xmm6,%xmm4  ## gb
    subpd     %xmm10,%xmm8 ## lj
    movapd    %xmm4,nb430_epsgb(%rsp)   ## gb eps
    movapd    %xmm8,nb430_eps(%rsp)   ## lj eps

        movq nb430_GBtab(%rbp),%rsi
        movq nb430_VFtab(%rbp),%rdi

    ## load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
    movapd (%rsi,%r12,8),%xmm0         ## Y1c F1c
    movapd (%rsi,%r13,8),%xmm12        ## Y2c F2c
    movapd (%rdi,%r14,8),%xmm4         ## Y1d F1d
    movapd (%rdi,%r15,8),%xmm13        ## Y2d F2d
    movapd 32(%rdi,%r14,8),%xmm8       ## Y1r F1r
    movapd 32(%rdi,%r15,8),%xmm14      ## Y2r F2r
        movapd %xmm0,%xmm1
        movapd %xmm4,%xmm5
        movapd %xmm8,%xmm9
        unpcklpd %xmm12,%xmm0   ## Y1c Y2c 
        unpckhpd %xmm12,%xmm1   ## F1c F2c 
        unpcklpd %xmm13,%xmm4   ## Y1d Y2d 
        unpckhpd %xmm13,%xmm5   ## F1d F2d 
        unpcklpd %xmm14,%xmm8   ## Y1r Y2r 
        unpckhpd %xmm14,%xmm9   ## F1r F2r 

    movapd 16(%rsi,%r12,8),%xmm2       ## G1c H1c
    movapd 16(%rsi,%r13,8),%xmm12      ## G2c H2c
    movapd 16(%rdi,%r14,8),%xmm6       ## G1d H1d
    movapd 16(%rdi,%r15,8),%xmm13      ## G2d H2d
    movapd 48(%rdi,%r14,8),%xmm10      ## G1r H1r
    movapd 48(%rdi,%r15,8),%xmm14      ## G2r H2r
        movapd %xmm2,%xmm3
        movapd %xmm6,%xmm7
        movapd %xmm10,%xmm11
        unpcklpd %xmm12,%xmm2   ## G1c G2c 
        unpckhpd %xmm12,%xmm3   ## H1c H2c 
        unpcklpd %xmm13,%xmm6   ## G1d G2d 
        unpckhpd %xmm13,%xmm7   ## H1d H2d 
        unpcklpd %xmm14,%xmm10  ## G1r G2r 
        unpckhpd %xmm14,%xmm11  ## H1r H2r 
    ## table data ready. Coul GB in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11
    movq nb430_vdwparam(%rbp),%rdi

    movapd nb430_epsgb(%rsp),%xmm12
    movapd nb430_eps(%rsp),%xmm13

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm13,%xmm11
    mulpd  %xmm12,%xmm2    ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm13,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm13,%xmm11

    movlpd (%rdi,%r8,8),%xmm14
    movlpd 8(%rdi,%r8,8),%xmm15

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    movhpd (%rdi,%r9,8),%xmm14
    movhpd 8(%rdi,%r9,8),%xmm15

    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm13,%xmm9
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  nb430_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulpd  %xmm14,%xmm5  ## vnb6
    mulpd  %xmm15,%xmm9  ## vnb12
    mulpd  nb430_qq(%rsp),%xmm3      ## FF*qq = fij
    mulpd  %xmm14,%xmm7  ## fijD
    mulpd  %xmm15,%xmm11  ##fijR

    addpd  %xmm7,%xmm11 ## fijD+fijR
    mulpd  nb430_tsc(%rsp),%xmm11   ## (fijD+fijR)*tabscale

    ## accumulate Vvdwtot
    addpd  nb430_Vvdwtot(%rsp),%xmm5
    addpd  %xmm9,%xmm5
    movapd %xmm5,nb430_Vvdwtot(%rsp)

        movq nb430_dvda(%rbp),%rsi

        ## Calculate dVda
        mulpd nb430_gbscale(%rsp),%xmm3     ## fijC=qq*FF*gbscale
        movapd %xmm3,%xmm6
        mulpd  nb430_r(%rsp),%xmm6
        addpd  %xmm1,%xmm6  ## vcoul+fijC*r

    addpd  %xmm11,%xmm3 ## fijC+fijD+fijR

    ## increment vctot
        addpd  nb430_vctot(%rsp),%xmm1
    movapd %xmm1,nb430_vctot(%rsp)

        ## xmm6=(vcoul+fijC*r)
        xorpd  %xmm7,%xmm7
        subpd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

    ## the fj's - start by combiningg forces from memory 
    movq nb430_faction(%rbp),%rdi
        movlpd (%rdi,%r10,8),%xmm0
        movlpd 8(%rdi,%r10,8),%xmm1
        movlpd 16(%rdi,%r10,8),%xmm2
        movhpd (%rdi,%r11,8),%xmm0
        movhpd 8(%rdi,%r11,8),%xmm1
        movhpd 16(%rdi,%r11,8),%xmm2

        ## update dvdasum 
        addpd  nb430_dvdasum(%rsp),%xmm7
    movapd %xmm7,nb430_dvdasum(%rsp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        addsd  (%rsi,%rax,8),%xmm6
        addsd  (%rsi,%rbx,8),%xmm7
        movsd  %xmm6,(%rsi,%rax,8)
        movsd  %xmm7,(%rsi,%rbx,8)

        xorpd  %xmm4,%xmm4
        mulpd nb430_rinv(%rsp),%xmm3
        subpd  %xmm3,%xmm4

    movapd  %xmm4,%xmm9
    movapd  %xmm4,%xmm10
    movapd  %xmm4,%xmm11

    mulpd  nb430_dx(%rsp),%xmm9
    mulpd  nb430_dy(%rsp),%xmm10
    mulpd  nb430_dz(%rsp),%xmm11

        addpd %xmm9,%xmm0
        addpd %xmm10,%xmm1
        addpd %xmm11,%xmm2

        ## accumulate i forces
    addpd nb430_fix(%rsp),%xmm9
    addpd nb430_fiy(%rsp),%xmm10
    addpd nb430_fiz(%rsp),%xmm11

        movlpd %xmm0,(%rdi,%r10,8)
        movlpd %xmm1,8(%rdi,%r10,8)
        movlpd %xmm2,16(%rdi,%r10,8)

    movapd %xmm9,nb430_fix(%rsp)
    movapd %xmm10,nb430_fiy(%rsp)
    movapd %xmm11,nb430_fiz(%rsp)

        movhpd %xmm0,(%rdi,%r11,8)
        movhpd %xmm1,8(%rdi,%r11,8)
        movhpd %xmm2,16(%rdi,%r11,8)

    ## should we do one more iteration? 
        subl $2,nb430_innerk(%rsp)
        jl    _nb_kernel430_x86_64_sse2.nb430_checksingle
        jmp   _nb_kernel430_x86_64_sse2.nb430_unroll_loop
_nb_kernel430_x86_64_sse2.nb430_checksingle: 
        movl  nb430_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430_x86_64_sse2.nb430_dosingle
        jmp    _nb_kernel430_x86_64_sse2.nb430_updateouterdata
_nb_kernel430_x86_64_sse2.nb430_dosingle: 
        movq nb430_charge(%rbp),%rsi
        movq nb430_invsqrta(%rbp),%rdx
        movq nb430_pos(%rbp),%rdi
        movq  nb430_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        ## load isaj
        movq nb430_invsqrta(%rbp),%rsi
        movsd (%rsi,%rax,8),%xmm2
        mulsd  nb430_isai(%rsp),%xmm2
        movapd %xmm2,nb430_isaprod(%rsp)
        movapd %xmm2,%xmm1
        mulsd nb430_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb430_gbscale(%rsp)

    mulsd nb430_iq(%rsp),%xmm2
        movq nb430_charge(%rbp),%rsi     ## base of charge[] 
        movsd (%rsi,%rax,8),%xmm3
        mulsd  %xmm2,%xmm3
        movapd %xmm3,nb430_qq(%rsp)

        movq nb430_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb430_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb430_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb430_c6(%rsp)
        movapd %xmm6,nb430_c12(%rsp)

        movq nb430_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%rax,%rax,2),%r10     ## j3 

        ## move coordinate to xmm4-xmm6 
        movsd (%rsi,%r10,8),%xmm4
        movsd 8(%rsi,%r10,8),%xmm5
        movsd 16(%rsi,%r10,8),%xmm6

        movq   nb430_faction(%rbp),%rdi

        ## calc dr 
        subsd nb430_ix(%rsp),%xmm4
        subsd nb430_iy(%rsp),%xmm5
        subsd nb430_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb430_dx(%rsp)
        movapd %xmm5,nb430_dy(%rsp)
        movapd %xmm6,nb430_dz(%rsp)

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb430_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb430_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb430_r(%rsp)
        movapd %xmm0,nb430_rinv(%rsp)

    movapd %xmm4,%xmm8 ## r
        mulsd nb430_gbscale(%rsp),%xmm4
        mulsd nb430_tsc(%rsp),%xmm8

    ## truncate and convert to integers
    cvttsd2si %xmm4,%r12d ## gb
    cvttsd2si %xmm8,%r14d ## lj

    ## convert back to float
    cvtsi2sd  %r12d,%xmm6  ## gb
    cvtsi2sd  %r14d,%xmm10 ## lj

    ## multiply by 4 and 8, respectively
    shll   $2,%r12d  ## gb
    shll   $3,%r14d  ## lj

    ## GB indices: r10   LJ indices: r12

    ## calculate eps
    subsd     %xmm6,%xmm4  ## gb
    subsd     %xmm10,%xmm8 ## lj
    movapd    %xmm4,nb430_epsgb(%rsp)   ## gb eps
    movapd    %xmm8,nb430_eps(%rsp)   ## lj eps

        movq nb430_GBtab(%rbp),%rsi
        movq nb430_VFtab(%rbp),%rdi

    ## load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
    movapd (%rsi,%r12,8),%xmm0         ## Y1c F1c
    movapd (%rdi,%r14,8),%xmm4         ## Y1d F1d
    movapd 32(%rdi,%r14,8),%xmm8       ## Y1r F1r
        movhlps %xmm0,%xmm1
        movhlps %xmm4,%xmm5
        movhlps %xmm8,%xmm9

    movapd 16(%rsi,%r12,8),%xmm2       ## G1c H1c
    movapd 16(%rdi,%r14,8),%xmm6       ## G1d H1d
    movapd 48(%rdi,%r14,8),%xmm10      ## G1r H1r
        movhlps %xmm2,%xmm3
        movhlps %xmm6,%xmm7
        movhlps %xmm10,%xmm11
    ## table data ready. Coul GB in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movapd nb430_epsgb(%rsp),%xmm12
    movapd nb430_eps(%rsp),%xmm13

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm13,%xmm11
    mulsd  %xmm12,%xmm2    ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm13,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm13,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm13,%xmm9
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  nb430_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulsd  nb430_c6(%rsp),%xmm5     ## vnb6
    mulsd  nb430_c12(%rsp),%xmm9     ## vnb12
    mulsd  nb430_qq(%rsp),%xmm3      ## FF*qq = fij
    mulsd  nb430_c6(%rsp),%xmm7     ## fijD
    mulsd  nb430_c12(%rsp),%xmm11     ##fijR

    addsd  %xmm7,%xmm11 ## fijD+fijR
    mulsd  nb430_tsc(%rsp),%xmm11   ## (fijD+fijR)*tabscale

    ## accumulate Vvdwtot
    addsd  nb430_Vvdwtot(%rsp),%xmm5
    addsd  %xmm9,%xmm5
    movsd %xmm5,nb430_Vvdwtot(%rsp)

        movq nb430_dvda(%rbp),%rsi

        ## Calculate dVda
        mulsd nb430_gbscale(%rsp),%xmm3     ## fijC=qq*FF*gbscale
        movapd %xmm3,%xmm6
        mulsd  nb430_r(%rsp),%xmm6
        addsd  %xmm1,%xmm6  ## vcoul+fijC*r

    addsd  %xmm11,%xmm3 ## fijC+fijD+fijR

    ## increment vctot
        addsd  nb430_vctot(%rsp),%xmm1
    movsd %xmm1,nb430_vctot(%rsp)

        ## xmm6=(vcoul+fijC*r)
        xorpd  %xmm7,%xmm7
        subsd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        ## update dvdasum 
        addsd  nb430_dvdasum(%rsp),%xmm7
    movsd %xmm7,nb430_dvdasum(%rsp)

        ## update j atoms dvdaj
        addsd  (%rsi,%rax,8),%xmm6
        movsd  %xmm6,(%rsi,%rax,8)

        xorpd  %xmm4,%xmm4
        mulsd nb430_rinv(%rsp),%xmm3
        subsd  %xmm3,%xmm4

    movapd  %xmm4,%xmm9
    movapd  %xmm4,%xmm10
    movapd  %xmm4,%xmm11

    mulsd  nb430_dx(%rsp),%xmm9
    mulsd  nb430_dy(%rsp),%xmm10
    mulsd  nb430_dz(%rsp),%xmm11

    movapd %xmm9,%xmm3
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm5

        ## accumulate i forces
    addsd nb430_fix(%rsp),%xmm9
    addsd nb430_fiy(%rsp),%xmm10
    addsd nb430_fiz(%rsp),%xmm11
    movsd %xmm9,nb430_fix(%rsp)
    movsd %xmm10,nb430_fiy(%rsp)
    movsd %xmm11,nb430_fiz(%rsp)

    movq nb430_faction(%rbp),%rdi
        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%r10,8),%xmm3
        addsd 8(%rdi,%r10,8),%xmm4
        addsd 16(%rdi,%r10,8),%xmm5
        movsd %xmm3,(%rdi,%r10,8)
        movsd %xmm4,8(%rdi,%r10,8)
        movsd %xmm5,16(%rdi,%r10,8)

_nb_kernel430_x86_64_sse2.nb430_updateouterdata: 
        movl  nb430_ii3(%rsp),%ecx
        movq  nb430_faction(%rbp),%rdi
        movq  nb430_fshift(%rbp),%rsi
        movl  nb430_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb430_fix(%rsp),%xmm0
        movapd nb430_fiy(%rsp),%xmm1
        movapd nb430_fiz(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## increment fshift force  
        movsd  (%rsi,%rdx,8),%xmm3
        movsd  8(%rsi,%rdx,8),%xmm4
        movsd  16(%rsi,%rdx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rsi,%rdx,8)
        movsd  %xmm4,8(%rsi,%rdx,8)
        movsd  %xmm5,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb430_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb430_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb430_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb430_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb430_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb430_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate dVda and update it 
        movapd nb430_dvdasum(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        movl nb430_ii(%rsp),%edx
        movq nb430_dvda(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb430_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel430_x86_64_sse2.nb430_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430_n(%rsp)
        jmp _nb_kernel430_x86_64_sse2.nb430_outer
_nb_kernel430_x86_64_sse2.nb430_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel430_x86_64_sse2.nb430_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430_x86_64_sse2.nb430_threadloop
_nb_kernel430_x86_64_sse2.nb430_end: 
        movl nb430_nouter(%rsp),%eax
        movl nb430_ninner(%rsp),%ebx
        movq nb430_outeriter(%rbp),%rcx
        movq nb430_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $536,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel430nf_x86_64_sse2
.globl _nb_kernel430nf_x86_64_sse2
nb_kernel430nf_x86_64_sse2:     
_nb_kernel430nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb430nf_fshift, 16
.set nb430nf_gid, 24
.set nb430nf_pos, 32
.set nb430nf_faction, 40
.set nb430nf_charge, 48
.set nb430nf_p_facel, 56
.set nb430nf_argkrf, 64
.set nb430nf_argcrf, 72
.set nb430nf_Vc, 80
.set nb430nf_type, 88
.set nb430nf_p_ntype, 96
.set nb430nf_vdwparam, 104
.set nb430nf_Vvdw, 112
.set nb430nf_p_tabscale, 120
.set nb430nf_VFtab, 128
.set nb430nf_invsqrta, 136
.set nb430nf_dvda, 144
.set nb430nf_p_gbtabscale, 152
.set nb430nf_GBtab, 160
.set nb430nf_p_nthreads, 168
.set nb430nf_count, 176
.set nb430nf_mtx, 184
.set nb430nf_outeriter, 192
.set nb430nf_inneriter, 200
.set nb430nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb430nf_ix, 0
.set nb430nf_iy, 16
.set nb430nf_iz, 32
.set nb430nf_iq, 48
.set nb430nf_gbtsc, 64
.set nb430nf_tsc, 80
.set nb430nf_qq, 96
.set nb430nf_c6, 112
.set nb430nf_c12, 128
.set nb430nf_vctot, 144
.set nb430nf_Vvdwtot, 160
.set nb430nf_half, 176
.set nb430nf_three, 192
.set nb430nf_r, 208
.set nb430nf_isai, 224
.set nb430nf_isaprod, 240
.set nb430nf_gbscale, 256
.set nb430nf_nri, 272
.set nb430nf_iinr, 280
.set nb430nf_jindex, 288
.set nb430nf_jjnr, 296
.set nb430nf_shift, 304
.set nb430nf_shiftvec, 312
.set nb430nf_facel, 320
.set nb430nf_innerjjnr, 328
.set nb430nf_is3, 336
.set nb430nf_ii3, 340
.set nb430nf_ntia, 344
.set nb430nf_innerk, 348
.set nb430nf_n, 352
.set nb430nf_nn1, 356
.set nb430nf_ntype, 360
.set nb430nf_nouter, 364
.set nb430nf_ninner, 368
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $392,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb430nf_nouter(%rsp)
        movl %eax,nb430nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb430nf_nri(%rsp)
        movq %rsi,nb430nf_iinr(%rsp)
        movq %rdx,nb430nf_jindex(%rsp)
        movq %rcx,nb430nf_jjnr(%rsp)
        movq %r8,nb430nf_shift(%rsp)
        movq %r9,nb430nf_shiftvec(%rsp)
        movq nb430nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb430nf_ntype(%rsp)
        movq nb430nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb430nf_facel(%rsp)

        movq nb430nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb430nf_tsc(%rsp)

        movq nb430nf_p_gbtabscale(%rbp),%rbx
        movsd (%rbx),%xmm4
        shufpd $0,%xmm4,%xmm4
        movapd %xmm4,nb430nf_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb430nf_half(%rsp)
        movl %ebx,nb430nf_half+4(%rsp)
        movsd nb430nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb430nf_half(%rsp)
        movapd %xmm3,nb430nf_three(%rsp)

_nb_kernel430nf_x86_64_sse2.nb430nf_threadloop: 
        movq  nb430nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel430nf_x86_64_sse2.nb430nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430nf_x86_64_sse2.nb430nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb430nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb430nf_n(%rsp)
        movl %ebx,nb430nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel430nf_x86_64_sse2.nb430nf_outerstart
        jmp _nb_kernel430nf_x86_64_sse2.nb430nf_end

_nb_kernel430nf_x86_64_sse2.nb430nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430nf_nouter(%rsp),%ebx
        movl %ebx,nb430nf_nouter(%rsp)

_nb_kernel430nf_x86_64_sse2.nb430nf_outer: 
        movq  nb430nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb430nf_is3(%rsp)            ## store is3 

        movq  nb430nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb430nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb430nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb430nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb430nf_invsqrta(%rbp),%rdx       ## load invsqrta[ii]
        movsd (%rdx,%rbx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movq  nb430nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb430nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb430nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb430nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb430nf_iq(%rsp)
        movapd %xmm4,nb430nf_isai(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb430nf_ix(%rsp)
        movapd %xmm1,nb430nf_iy(%rsp)
        movapd %xmm2,nb430nf_iz(%rsp)

        movl  %ebx,nb430nf_ii3(%rsp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb430nf_vctot(%rsp)
        movapd %xmm4,nb430nf_Vvdwtot(%rsp)

        movq  nb430nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb430nf_pos(%rbp),%rsi
        movq  nb430nf_faction(%rbp),%rdi
        movq  nb430nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb430nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb430nf_ninner(%rsp),%ecx
        movl  %ecx,nb430nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb430nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel430nf_x86_64_sse2.nb430nf_unroll_loop
        jmp   _nb_kernel430nf_x86_64_sse2.nb430nf_checksingle
_nb_kernel430nf_x86_64_sse2.nb430nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb430nf_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb430nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        ## load isaj
        movq nb430nf_invsqrta(%rbp),%rsi
        movlpd (%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm2
        mulpd  nb430nf_isai(%rsp),%xmm2
        movapd %xmm2,nb430nf_isaprod(%rsp)
        movapd %xmm2,%xmm1
        mulpd nb430nf_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb430nf_gbscale(%rsp)

        movq nb430nf_charge(%rbp),%rsi     ## base of charge[] 
        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        mulpd nb430nf_iq(%rsp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb430nf_qq(%rsp)

        movq nb430nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ecx
        movl (%rsi,%rbx,4),%edx
        movq nb430nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb430nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx

        movlpd (%rsi,%rcx,8),%xmm6      ## c6a
        movlpd (%rsi,%rdx,8),%xmm7      ## c6b
        movhpd 8(%rsi,%rcx,8),%xmm6     ## c6a c12a 
        movhpd 8(%rsi,%rdx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb430nf_c6(%rsp)
        movapd %xmm6,nb430nf_c12(%rsp)

        movq nb430nf_pos(%rbp),%rsi             ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        movq   nb430nf_faction(%rbp),%rdi

        ## move nb430nf_ix-iz to xmm4-xmm6 
        movapd nb430nf_ix(%rsp),%xmm4
        movapd nb430nf_iy(%rsp),%xmm5
        movapd nb430nf_iz(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb430nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb430nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb430nf_r(%rsp)
        mulpd nb430nf_gbscale(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movq nb430nf_GBtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx          ## indices in eax/ebx 

        ## Coulomb 
        movapd (%rsi,%rcx,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rdx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rcx,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rdx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb430nf_qq(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb430nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb430nf_vctot(%rsp)

        movapd nb430nf_r(%rsp),%xmm4
        mulpd  nb430nf_tsc(%rsp),%xmm4
        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8

        movq nb430nf_VFtab(%rbp),%rsi

        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx          ## indices in eax/ebx 

        ## Dispersion 
        movapd (%rsi,%rcx,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rdx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rcx,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rdx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb430nf_c6(%rsp),%xmm5    ## Vvdw6
        addpd  nb430nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb430nf_Vvdwtot(%rsp)

        ## Repulsion 
        movapd 32(%rsi,%rcx,8),%xmm4    ## Y1 F1        
        movapd 32(%rsi,%rdx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 48(%rsi,%rcx,8),%xmm6    ## G1 H1        
        movapd 48(%rsi,%rdx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb430nf_c12(%rsp),%xmm5   ## Vvdw12 
        addpd  nb430nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb430nf_Vvdwtot(%rsp)
        xorpd  %xmm4,%xmm4

        ## should we do one more iteration? 
        subl $2,nb430nf_innerk(%rsp)
        jl    _nb_kernel430nf_x86_64_sse2.nb430nf_checksingle
        jmp   _nb_kernel430nf_x86_64_sse2.nb430nf_unroll_loop
_nb_kernel430nf_x86_64_sse2.nb430nf_checksingle: 
        movl  nb430nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430nf_x86_64_sse2.nb430nf_dosingle
        jmp    _nb_kernel430nf_x86_64_sse2.nb430nf_updateouterdata
_nb_kernel430nf_x86_64_sse2.nb430nf_dosingle: 
        movq nb430nf_charge(%rbp),%rsi
        movq nb430nf_invsqrta(%rbp),%rdx
        movq nb430nf_pos(%rbp),%rdi
        movq  nb430nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%rdx,%rax,8),%xmm7
        movlpd (%rsi,%rax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb430nf_isai(%rsp),%xmm7
        movapd %xmm7,nb430nf_isaprod(%rsp)
        movapd %xmm7,%xmm1
        mulpd nb430nf_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb430nf_gbscale(%rsp)

        mulsd  nb430nf_iq(%rsp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb430nf_qq(%rsp)

        movq nb430nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%edx
        movq nb430nf_vdwparam(%rbp),%rsi
        shll %edx
        movl nb430nf_ntia(%rsp),%edi
        addl %edi,%edx

        movlpd (%rsi,%rdx,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rdx,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb430nf_c6(%rsp)
        movapd %xmm6,nb430nf_c12(%rsp)

        movq nb430nf_pos(%rbp),%rsi             ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        movq   nb430nf_faction(%rbp),%rdi

        ## move nb430nf_ix-iz to xmm4-xmm6 
        movapd nb430nf_ix(%rsp),%xmm4
        movapd nb430nf_iy(%rsp),%xmm5
        movapd nb430nf_iz(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb430nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb430nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        movsd %xmm4,nb430nf_r(%rsp)
        mulsd nb430nf_gbscale(%rsp),%xmm4

        cvttsd2si %xmm4,%edx    ## mm6 = lu idx 
        cvtsi2sd %edx,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%edx            ## idx *= 4 
        movq nb430nf_GBtab(%rbp),%rsi

        ## Coulomb 
        movapd (%rsi,%rdx,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rdx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb430nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addsd  nb430nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb430nf_vctot(%rsp)

        movsd nb430nf_r(%rsp),%xmm4
        mulsd  nb430nf_tsc(%rsp),%xmm4
        cvttsd2si %xmm4,%edx    ## mm6 = lu idx 
        cvtsi2sd %edx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2

        shll $3,%edx

        movq nb430nf_VFtab(%rbp),%rsi

        ## Dispersion 
        movapd (%rsi,%rdx,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rdx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb430nf_c6(%rsp),%xmm5    ## Vvdw6
        addsd  nb430nf_Vvdwtot(%rsp),%xmm5
        movlpd %xmm5,nb430nf_Vvdwtot(%rsp)

        ## Repulsion 
        movapd 32(%rsi,%rdx,8),%xmm4    ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 48(%rsi,%rdx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  nb430nf_c12(%rsp),%xmm5   ## Vvdw12 
        addsd  nb430nf_Vvdwtot(%rsp),%xmm5
        movlpd %xmm5,nb430nf_Vvdwtot(%rsp)
_nb_kernel430nf_x86_64_sse2.nb430nf_updateouterdata: 
        ## get n from stack
        movl nb430nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb430nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb430nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb430nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb430nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb430nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb430nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel430nf_x86_64_sse2.nb430nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430nf_n(%rsp)
        jmp _nb_kernel430nf_x86_64_sse2.nb430nf_outer
_nb_kernel430nf_x86_64_sse2.nb430nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel430nf_x86_64_sse2.nb430nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430nf_x86_64_sse2.nb430nf_threadloop
_nb_kernel430nf_x86_64_sse2.nb430nf_end: 
        movl nb430nf_nouter(%rsp),%eax
        movl nb430nf_ninner(%rsp),%ebx
        movq nb430nf_outeriter(%rbp),%rcx
        movq nb430nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $392,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




