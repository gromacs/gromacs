##
## $Id$
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







.globl nb_kernel430_x86_64_sse
.globl _nb_kernel430_x86_64_sse
nb_kernel430_x86_64_sse:        
_nb_kernel430_x86_64_sse:       
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
        ## bottom of stack is cache-aligned for sse use 
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

        subq $552,%rsp          ## local variable stack space (n*16+8)

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
        movss (%rsi),%xmm0
        movss %xmm0,nb430_facel(%rsp)

        movq nb430_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb430_tsc(%rsp)

        movq nb430_p_gbtabscale(%rbp),%rbx
        movss (%rbx),%xmm4
        shufps $0,%xmm4,%xmm4
        movaps %xmm4,nb430_gbtsc(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb430_half(%rsp)
        movss nb430_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb430_half(%rsp)
        movaps %xmm3,nb430_three(%rsp)

_nb_kernel430_x86_64_sse.nb430_threadloop: 
        movq  nb430_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel430_x86_64_sse.nb430_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430_x86_64_sse.nb430_spinlock

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
        jg  _nb_kernel430_x86_64_sse.nb430_outerstart
        jmp _nb_kernel430_x86_64_sse.nb430_end

_nb_kernel430_x86_64_sse.nb430_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430_nouter(%rsp),%ebx
        movl %ebx,nb430_nouter(%rsp)

_nb_kernel430_x86_64_sse.nb430_outer: 
        movq  nb430_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb430_is3(%rsp)      ## store is3 

        movq  nb430_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb430_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 
        movl  %ebx,nb430_ii(%rsp)

        movq  nb430_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb430_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb430_invsqrta(%rbp),%rdx         ## load invsqrta[ii]
        movss (%rdx,%rbx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        movq  nb430_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb430_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb430_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb430_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb430_iq(%rsp)
        movaps %xmm4,nb430_isai(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb430_ix(%rsp)
        movaps %xmm1,nb430_iy(%rsp)
        movaps %xmm2,nb430_iz(%rsp)

        movl  %ebx,nb430_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb430_vctot(%rsp)
        movaps %xmm4,nb430_Vvdwtot(%rsp)
        movaps %xmm4,nb430_dvdasum(%rsp)
        movaps %xmm4,nb430_fix(%rsp)
        movaps %xmm4,nb430_fiy(%rsp)
        movaps %xmm4,nb430_fiz(%rsp)

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
        subl  $4,%edx
        addl  nb430_ninner(%rsp),%ecx
        movl  %ecx,nb430_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb430_innerk(%rsp)      ## number of innerloop atoms

        jge   _nb_kernel430_x86_64_sse.nb430_unroll_loop
        jmp   _nb_kernel430_x86_64_sse.nb430_finish_inner
_nb_kernel430_x86_64_sse.nb430_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb430_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb430_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        ## load isaj
        movq nb430_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7
        movaps nb430_isai(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all isaj in xmm3 
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb430_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb430_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb430_gbscale(%rsp)

        movq nb430_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        mulps nb430_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430_qq(%rsp)

    ## vdw parameters
        movq nb430_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        movl (%rsi,%rcx,4),%r14d
        movl (%rsi,%rdx,4),%r15d
        shll %r12d
        shll %r13d
        shll %r14d
        shll %r15d
    movl nb430_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d
        addl %edi,%r14d
        addl %edi,%r15d

        movq nb430_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movlps (%rsi,%r14,4),%xmm7
        movhps (%rsi,%r13,4),%xmm3
        movhps (%rsi,%r15,4),%xmm7

        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb430_c6(%rsp)
    movaps %xmm3,nb430_c12(%rsp)

        movq nb430_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## jnr
        lea  (%rbx,%rbx,2),%r9
        lea  (%rcx,%rcx,2),%r10
        lea  (%rdx,%rdx,2),%r11

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%r8,4),%xmm4
        movlps (%rsi,%r10,4),%xmm5
        movss 8(%rsi,%r8,4),%xmm2
        movss 8(%rsi,%r10,4),%xmm6

        movhps (%rsi,%r9,4),%xmm4
        movhps (%rsi,%r11,4),%xmm5

        movss 8(%rsi,%r9,4),%xmm0
        movss 8(%rsi,%r11,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## calc dr 
        subps nb430_ix(%rsp),%xmm0
        subps nb430_iy(%rsp),%xmm1
        subps nb430_iz(%rsp),%xmm2

        ## store dr 
        movaps %xmm0,nb430_dx(%rsp)
        movaps %xmm1,nb430_dy(%rsp)
        movaps %xmm2,nb430_dz(%rsp)

    movd %r8,%mm0 ## store j3
    movd %r9,%mm1
    movd %r10,%mm2
    movd %r11,%mm3

        ## square it 
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
    movaps %xmm0,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb430_r(%rsp)
    movaps %xmm0,nb430_rinv(%rsp)

    movaps %xmm4,%xmm8   ## r
        mulps nb430_gbscale(%rsp),%xmm4   ## rgbtab
    mulps nb430_tsc(%rsp),%xmm8      ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm4,%xmm5 ## gb
    cvttps2dq %xmm8,%xmm9 ## lj

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm6  ## gb
    cvtdq2ps  %xmm9,%xmm10 ## lj

    ## multiply by 4 and 8, respectively
    pslld   $2,%xmm5  ## gb
    pslld   $3,%xmm9  ## lj

    ## move to integer registers
    movhlps %xmm5,%xmm7    ## gb
    movhlps %xmm9,%xmm11   ## lj
    movd    %xmm5,%r8d      ## gb
    movd    %xmm9,%r12d     ## lj
    movd    %xmm7,%r10d     ## gb
    movd    %xmm11,%r14d    ## lj
    pshufd $1,%xmm5,%xmm5 ## gb
    pshufd $1,%xmm9,%xmm9 ## lj
    pshufd $1,%xmm7,%xmm7 ## gb
    pshufd $1,%xmm11,%xmm11 ## lj
    movd    %xmm5,%r9d      ## gb
    movd    %xmm9,%r13d     ## lj
    movd    %xmm7,%r11d     ## gb
    movd    %xmm11,%r15d    ## lj
    ## GB indices: r8-r11   LJ indices: r12-r15

    ## calculate eps
    subps     %xmm6,%xmm4  ## gb
    subps     %xmm10,%xmm8 ## lj
    movaps    %xmm4,nb430_epsgb(%rsp)   ## gb eps
    movaps    %xmm8,nb430_eps(%rsp)   ## lj eps

        movq nb430_GBtab(%rbp),%rsi
        movq nb430_VFtab(%rbp),%rdi

    ## load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
        movlps (%rsi,%r8,4),%xmm1         ## Y1c F1c 
        movlps (%rdi,%r12,4),%xmm5        ## Y1d F1d 
        movlps 16(%rdi,%r12,4),%xmm9      ## Y1r F1r 

        movlps (%rsi,%r10,4),%xmm3        ## Y3c F3c 
        movlps (%rdi,%r14,4),%xmm7        ## Y3d F3d 
        movlps 16(%rdi,%r14,4),%xmm11     ## Y3r F3r 

        movhps (%rsi,%r9,4),%xmm1         ## Y1c F1c Y2c F2c
        movhps (%rdi,%r13,4),%xmm5        ## Y1d F1d Y2d F2d
        movhps 16(%rdi,%r13,4),%xmm9      ## Y1r F1r Y2r F2r

        movhps (%rsi,%r11,4),%xmm3        ## Y3c F3c Y4c F4c
        movhps (%rdi,%r15,4),%xmm7        ## Y3d F3d Y4d F4d
        movhps 16(%rdi,%r15,4),%xmm11     ## Y3r F3r Y4r F4r

    movaps %xmm1,%xmm0
    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm3,%xmm0 ## 10001000   => Y1c Y2c Y3c Y4c
        shufps $136,%xmm7,%xmm4 ## 10001000   => Y1d Y2d Y3d Y4d
        shufps $136,%xmm11,%xmm8 ## 10001000  => Y1r Y2r Y3r Y4r
        shufps $221,%xmm3,%xmm1 ## 11011101   => F1c F2c F3c F4c
        shufps $221,%xmm7,%xmm5 ## 11011101   => F1d F2d F3d F4d
        shufps $221,%xmm11,%xmm9 ## 11011101  => F1r F2r F3r F4r

        movlps 8(%rsi,%r8,4),%xmm3         ## G1c H1c 
        movlps 8(%rdi,%r12,4),%xmm7        ## G1d H1d 
        movlps 24(%rdi,%r12,4),%xmm11      ## G1r H1r 

        movlps 8(%rsi,%r10,4),%xmm12       ## G3c H3c 
        movlps 8(%rdi,%r14,4),%xmm13       ## G3d H3d 
        movlps 24(%rdi,%r14,4),%xmm14      ## G3r H3r 

        movhps 8(%rsi,%r9,4),%xmm3         ## G1c H1c G2c H2c
        movhps 8(%rdi,%r13,4),%xmm7        ## G1d H1d G2d H2d
        movhps 24(%rdi,%r13,4),%xmm11      ## G1r H1r G2r H2r

        movhps 8(%rsi,%r11,4),%xmm12       ## G3c H3c G4c H4c
        movhps 8(%rdi,%r15,4),%xmm13       ## G3d H3d G4d H4d
        movhps 24(%rdi,%r15,4),%xmm14      ## G3r H3r G4r H4r
    movaps %xmm3,%xmm2
    movaps %xmm7,%xmm6
    movaps %xmm11,%xmm10

        shufps $136,%xmm12,%xmm2 ## 10001000  => G1c G2c G3c G4c
        shufps $136,%xmm13,%xmm6 ## 10001000  => G1d G2d G3d G4d
        shufps $136,%xmm14,%xmm10 ## 10001000 => G1r G2r G3r G4r
        shufps $221,%xmm12,%xmm3 ## 11011101  => H1c H2c H3c H4c
        shufps $221,%xmm13,%xmm7 ## 11011101  => H1d H2d H3d H4d
        shufps $221,%xmm14,%xmm11 ## 11011101 => H1r H2r H3r H4r
    ## table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps nb430_epsgb(%rsp),%xmm12
    movaps nb430_eps(%rsp),%xmm13

    mulps  %xmm12,%xmm3  ## Heps
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11
    mulps  %xmm12,%xmm2    ## Geps
    mulps  %xmm13,%xmm6
    mulps  %xmm13,%xmm10
    mulps  %xmm12,%xmm3  ## Heps2
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11

    addps  %xmm2,%xmm1  ## F+Geps
    addps  %xmm6,%xmm5
    addps  %xmm10,%xmm9
    addps  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addps  %xmm7,%xmm5
    addps  %xmm11,%xmm9
    addps  %xmm3,%xmm3   ## 2*Heps2
    addps  %xmm7,%xmm7
    addps  %xmm11,%xmm11
    addps  %xmm2,%xmm3   ## 2*Heps2+Geps
    addps  %xmm6,%xmm7
    addps  %xmm10,%xmm11
    addps  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addps  %xmm5,%xmm7
    addps  %xmm9,%xmm11
    mulps  %xmm12,%xmm1  ## eps*Fp
    mulps  %xmm13,%xmm5
    mulps  %xmm13,%xmm9
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  nb430_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulps  nb430_c6(%rsp),%xmm5     ## vnb6
    mulps  nb430_c12(%rsp),%xmm9     ## vnb12
    mulps  nb430_qq(%rsp),%xmm3      ## FF*qq = fij
    mulps  nb430_c6(%rsp),%xmm7     ## fijD
    mulps  nb430_c12(%rsp),%xmm11     ##fijR

    addps  %xmm7,%xmm11 ## fijD+fijR
    mulps  nb430_tsc(%rsp),%xmm11   ## (fijD+fijR)*tabscale

    ## accumulate Vvdwtot
    addps  nb430_Vvdwtot(%rsp),%xmm5
    addps  %xmm9,%xmm5
    movaps %xmm5,nb430_Vvdwtot(%rsp)

        movq nb430_dvda(%rbp),%rsi

        ## Calculate dVda
        mulps nb430_gbscale(%rsp),%xmm3     ## fijC=qq*FF*gbscale
        movaps %xmm3,%xmm6
        mulps  nb430_r(%rsp),%xmm6
        addps  %xmm1,%xmm6  ## vcoul+fijC*r

    addps  %xmm11,%xmm3 ## fijC+fijD+fijR

    ## increment vctot
        addps  nb430_vctot(%rsp),%xmm1
    movaps %xmm1,nb430_vctot(%rsp)

        ## xmm6=(vcoul+fijC*r)
        xorps  %xmm7,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum 
        addps  nb430_dvdasum(%rsp),%xmm7
    movaps %xmm7,nb430_dvdasum(%rsp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        movaps  %xmm6,%xmm5
        movaps  %xmm7,%xmm4
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm4,%xmm4

        ## xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
        addss  (%rsi,%rax,4),%xmm6
        addss  (%rsi,%rbx,4),%xmm5
        addss  (%rsi,%rcx,4),%xmm7
        addss  (%rsi,%rdx,4),%xmm4
        movss  %xmm6,(%rsi,%rax,4)
        movss  %xmm5,(%rsi,%rbx,4)
        movss  %xmm7,(%rsi,%rcx,4)
        movss  %xmm4,(%rsi,%rdx,4)

        xorps  %xmm4,%xmm4
        mulps nb430_rinv(%rsp),%xmm3
        subps  %xmm3,%xmm4

    movd %mm0,%r8  ## fetch j3
    movd %mm1,%r9
    movd %mm2,%r10
    movd %mm3,%r11

    movaps  %xmm4,%xmm9
    movaps  %xmm4,%xmm10
    movaps  %xmm4,%xmm11

    mulps  nb430_dx(%rsp),%xmm9
    mulps  nb430_dy(%rsp),%xmm10
    mulps  nb430_dz(%rsp),%xmm11

        ## accumulate i forces
    movaps nb430_fix(%rsp),%xmm12
    movaps nb430_fiy(%rsp),%xmm13
    movaps nb430_fiz(%rsp),%xmm14
    addps %xmm9,%xmm12
    addps %xmm10,%xmm13
    addps %xmm11,%xmm14
    movaps %xmm12,nb430_fix(%rsp)
    movaps %xmm13,nb430_fiy(%rsp)
    movaps %xmm14,nb430_fiz(%rsp)

        movq nb430_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%r10,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%r11,4),%xmm1 ## x3 y3 x4 y4

    movaps %xmm9,%xmm8
    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    unpckhps %xmm10,%xmm8 ## x3 y3 x4 y4

    ## update fjx and fjy
        addps  %xmm9,%xmm0
        addps  %xmm8,%xmm1

        movlps %xmm0,(%rsi,%r8,4)
        movlps %xmm1,(%rsi,%r10,4)
        movhps %xmm0,(%rsi,%r9,4)
        movhps %xmm1,(%rsi,%r11,4)

    ## xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd $1,%xmm11,%xmm10 ## fjz2 - - -
    movhlps %xmm11,%xmm9     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm8  ## fjz4 - - -

        addss  8(%rsi,%r8,4),%xmm11
        addss  8(%rsi,%r9,4),%xmm10
        addss  8(%rsi,%r10,4),%xmm9
        addss  8(%rsi,%r11,4),%xmm8
        movss  %xmm11,8(%rsi,%r8,4)
        movss  %xmm10,8(%rsi,%r9,4)
        movss  %xmm9,8(%rsi,%r10,4)
        movss  %xmm8,8(%rsi,%r11,4)

        ## should we do one more iteration? 
        subl $4,nb430_innerk(%rsp)
        jl    _nb_kernel430_x86_64_sse.nb430_finish_inner
        jmp   _nb_kernel430_x86_64_sse.nb430_unroll_loop
_nb_kernel430_x86_64_sse.nb430_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb430_innerk(%rsp)
        movl  nb430_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel430_x86_64_sse.nb430_dopair
        jmp   _nb_kernel430_x86_64_sse.nb430_checksingle
_nb_kernel430_x86_64_sse.nb430_dopair: 
        movq  nb430_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb430_innerjjnr(%rsp)

        ## load isaj
        movq nb430_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        movaps nb430_isai(%rsp),%xmm2
    unpcklps %xmm6,%xmm3
        mulps  %xmm3,%xmm2
    movaps %xmm2,nb430_isaprod(%rsp)

        movaps %xmm2,%xmm1
        mulps nb430_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb430_gbscale(%rsp)

        movq nb430_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
    unpcklps %xmm6,%xmm3
        mulps nb430_iq(%rsp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430_qq(%rsp)

    ## vdw parameters
        movq nb430_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        shll %r12d
        shll %r13d
    movl nb430_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb430_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb430_c6(%rsp)
    movaps %xmm3,nb430_c12(%rsp)

        movq nb430_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## j3
        lea  (%rbx,%rbx,2),%r9

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%r8,4),%xmm0       ## x1 y1 - - 
        movlps (%rsi,%r9,4),%xmm1       ## x2 y2 - - 

        movss 8(%rsi,%r8,4),%xmm2       ## z1 - - - 
        movss 8(%rsi,%r9,4),%xmm7       ## z2 - - - 

    unpcklps %xmm1,%xmm0 ## x1 x2 y1 y2
    movhlps  %xmm0,%xmm1 ## y1 y2 -  -
    unpcklps %xmm7,%xmm2 ## z1 z2 -  -

        ## calc dr 
        subps nb430_ix(%rsp),%xmm0
        subps nb430_iy(%rsp),%xmm1
        subps nb430_iz(%rsp),%xmm2

        ## store dr 
        movaps %xmm0,nb430_dx(%rsp)
        movaps %xmm1,nb430_dy(%rsp)
        movaps %xmm2,nb430_dz(%rsp)

        ## square it 
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
    movaps %xmm0,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb430_r(%rsp)
    movaps %xmm0,nb430_rinv(%rsp)

    movaps %xmm4,%xmm8   ## r
        mulps nb430_gbscale(%rsp),%xmm4   ## rgbtab
    mulps nb430_tsc(%rsp),%xmm8      ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm4,%xmm5 ## gb
    cvttps2dq %xmm8,%xmm9 ## lj

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm6  ## gb
    cvtdq2ps  %xmm9,%xmm10 ## lj

    ## multiply by 4 and 8, respectively
    pslld   $2,%xmm5  ## gb
    pslld   $3,%xmm9  ## lj

    ## move to integer registers
    movd    %xmm5,%r12d      ## gb
    movd    %xmm9,%r14d     ## lj
    pshufd $1,%xmm5,%xmm5  ## gb
    pshufd $1,%xmm9,%xmm9  ## lj
    movd    %xmm5,%r13d      ## gb
    movd    %xmm9,%r15d     ## lj
    ## GB indices: r12-r13   LJ indices: r14-r15

    ## calculate eps
    subps     %xmm6,%xmm4  ## gb
    subps     %xmm10,%xmm8 ## lj
    movaps    %xmm4,nb430_epsgb(%rsp)   ## gb eps
    movaps    %xmm8,nb430_eps(%rsp)   ## lj eps

        movq nb430_GBtab(%rbp),%rsi
        movq nb430_VFtab(%rbp),%rdi

    ## load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
        movlps (%rsi,%r12,4),%xmm0       ## Y1c F1c
        movlps (%rsi,%r13,4),%xmm1       ## Y2c F2c
        movlps (%rdi,%r14,4),%xmm4       ## Y1d F1d  
        movlps (%rdi,%r15,4),%xmm5       ## Y2d F2d
        movlps 16(%rdi,%r14,4),%xmm8     ## Y1r F1r
        movlps 16(%rdi,%r15,4),%xmm9     ## Y2r F2r

    unpcklps %xmm1,%xmm0
    movhlps  %xmm0,%xmm1
    unpcklps %xmm5,%xmm4
    movhlps  %xmm4,%xmm5
    unpcklps %xmm9,%xmm8
    movhlps  %xmm8,%xmm9
        movlps 8(%rsi,%r12,4),%xmm2       ## G1c H1c
        movlps 8(%rsi,%r13,4),%xmm3       ## G2c H2c
        movlps 8(%rdi,%r14,4),%xmm6       ## G1d H1d  
        movlps 8(%rdi,%r15,4),%xmm7       ## G2d H2d
        movlps 24(%rdi,%r14,4),%xmm10     ## G1r H1r
        movlps 24(%rdi,%r15,4),%xmm11     ## G2r H2r
    unpcklps %xmm3,%xmm2
    movhlps  %xmm2,%xmm3
    unpcklps %xmm7,%xmm6
    movhlps  %xmm6,%xmm7
    unpcklps %xmm11,%xmm10
    movhlps  %xmm10,%xmm11
    ## table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps nb430_epsgb(%rsp),%xmm12
    movaps nb430_eps(%rsp),%xmm13

    mulps  %xmm12,%xmm3  ## Heps
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11
    mulps  %xmm12,%xmm2    ## Geps
    mulps  %xmm13,%xmm6
    mulps  %xmm13,%xmm10
    mulps  %xmm12,%xmm3  ## Heps2
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11

    addps  %xmm2,%xmm1  ## F+Geps
    addps  %xmm6,%xmm5
    addps  %xmm10,%xmm9
    addps  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addps  %xmm7,%xmm5
    addps  %xmm11,%xmm9
    addps  %xmm3,%xmm3   ## 2*Heps2
    addps  %xmm7,%xmm7
    addps  %xmm11,%xmm11
    addps  %xmm2,%xmm3   ## 2*Heps2+Geps
    addps  %xmm6,%xmm7
    addps  %xmm10,%xmm11
    addps  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addps  %xmm5,%xmm7
    addps  %xmm9,%xmm11
    mulps  %xmm12,%xmm1  ## eps*Fp
    mulps  %xmm13,%xmm5
    mulps  %xmm13,%xmm9
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  nb430_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulps  nb430_c6(%rsp),%xmm5     ## vnb6
    mulps  nb430_c12(%rsp),%xmm9     ## vnb12
    mulps  nb430_qq(%rsp),%xmm3      ## FF*qq = fij
    mulps  nb430_c6(%rsp),%xmm7     ## fijD
    mulps  nb430_c12(%rsp),%xmm11     ##fijR

    addps  %xmm7,%xmm11 ## fijD+fijR
    mulps  nb430_tsc(%rsp),%xmm11   ## (fijD+fijR)*tabscale

    ## accumulate Vvdwtot
    addps  nb430_Vvdwtot(%rsp),%xmm5
    addps  %xmm9,%xmm5
    movlps %xmm5,nb430_Vvdwtot(%rsp)

        movq nb430_dvda(%rbp),%rsi

        ## Calculate dVda
        mulps nb430_gbscale(%rsp),%xmm3     ## fijC=qq*FF*gbscale
        movaps %xmm3,%xmm6
        mulps  nb430_r(%rsp),%xmm6
        addps  %xmm1,%xmm6  ## vcoul+fijC*r

    addps  %xmm11,%xmm3 ## fijC+fijD+fijR

    ## increment vctot
        addps  nb430_vctot(%rsp),%xmm1
    movlps %xmm1,nb430_vctot(%rsp)

        ## xmm6=(vcoul+fijC*r)
        xorps  %xmm7,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum 
        addps  nb430_dvdasum(%rsp),%xmm7
    movlps %xmm7,nb430_dvdasum(%rsp)

        ## update j atoms dvdaj
        movaps  %xmm6,%xmm5
        shufps $0x1,%xmm5,%xmm5

        ## xmm6=dvdaj1 xmm5=dvdaj2 
        addss  (%rsi,%rax,4),%xmm6
        addss  (%rsi,%rbx,4),%xmm5
        movss  %xmm6,(%rsi,%rax,4)
        movss  %xmm5,(%rsi,%rbx,4)

        xorps  %xmm4,%xmm4
        mulps nb430_rinv(%rsp),%xmm3
        subps  %xmm3,%xmm4

    movaps  %xmm4,%xmm9
    movaps  %xmm4,%xmm10
    movaps  %xmm4,%xmm11

    mulps  nb430_dx(%rsp),%xmm9
    mulps  nb430_dy(%rsp),%xmm10
    mulps  nb430_dz(%rsp),%xmm11


        ## accumulate i forces
    movaps nb430_fix(%rsp),%xmm12
    movaps nb430_fiy(%rsp),%xmm13
    movaps nb430_fiz(%rsp),%xmm14
    addps %xmm9,%xmm12
    addps %xmm10,%xmm13
    addps %xmm11,%xmm14
    movlps %xmm12,nb430_fix(%rsp)
    movlps %xmm13,nb430_fiy(%rsp)
    movlps %xmm14,nb430_fiz(%rsp)

        movq nb430_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2

    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    addps    %xmm9,%xmm0

        movlps %xmm0,(%rsi,%r8,4)
        movhps %xmm0,(%rsi,%r9,4)

    ## z forces
    pshufd $1,%xmm11,%xmm8
    addss  8(%rsi,%r8,4),%xmm11
    addss  8(%rsi,%r9,4),%xmm8
    movss  %xmm11,8(%rsi,%r8,4)
    movss  %xmm8,8(%rsi,%r9,4)

_nb_kernel430_x86_64_sse.nb430_checksingle:     
        movl  nb430_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430_x86_64_sse.nb430_dosingle
        jmp    _nb_kernel430_x86_64_sse.nb430_updateouterdata
_nb_kernel430_x86_64_sse.nb430_dosingle: 
        movq nb430_charge(%rbp),%rsi
        movq nb430_invsqrta(%rbp),%rdx
        movq nb430_pos(%rbp),%rdi
        movq  nb430_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        ## load isaj
        movq nb430_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movaps nb430_isai(%rsp),%xmm2
        mulss  %xmm3,%xmm2
    movaps %xmm2,nb430_isaprod(%rsp)

        movaps %xmm2,%xmm1
        mulss nb430_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb430_gbscale(%rsp)

        movq nb430_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        mulss nb430_iq(%rsp),%xmm2
        mulss  %xmm2,%xmm3
        movaps %xmm3,nb430_qq(%rsp)

    ## vdw parameters
        movq nb430_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb430_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb430_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
        movss 4(%rsi,%r12,4),%xmm3
    movaps %xmm0,nb430_c6(%rsp)
    movaps %xmm3,nb430_c12(%rsp)

        movq nb430_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## j3

        ## move four coordinates to xmm0-xmm2   
    movss  (%rsi,%r8,4),%xmm0
    movss  4(%rsi,%r8,4),%xmm1
    movss  8(%rsi,%r8,4),%xmm2

        ## calc dr 
        subss nb430_ix(%rsp),%xmm0
        subss nb430_iy(%rsp),%xmm1
        subss nb430_iz(%rsp),%xmm2

        ## store dr 
        movaps %xmm0,nb430_dx(%rsp)
        movaps %xmm1,nb430_dy(%rsp)
        movaps %xmm2,nb430_dz(%rsp)

        ## square it 
        mulss %xmm0,%xmm0
        mulss %xmm1,%xmm1
        mulss %xmm2,%xmm2
        addss %xmm1,%xmm0
        addss %xmm2,%xmm0
    movaps %xmm0,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movaps nb430_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430_half(%rsp),%xmm0
        subss %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 
        mulss %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb430_r(%rsp)
    movaps %xmm0,nb430_rinv(%rsp)

    movaps %xmm4,%xmm8   ## r
        mulss nb430_gbscale(%rsp),%xmm4   ## rgbtab
    mulss nb430_tsc(%rsp),%xmm8      ## rtab

    ## truncate and convert to integers
    cvttss2si %xmm4,%r12d ## gb
    cvttss2si %xmm8,%r14d ## lj

    ## convert back to float
    cvtsi2ss  %r12d,%xmm6  ## gb
    cvtsi2ss  %r14d,%xmm10 ## lj

    ## multiply by 4 and 8, respectively
    shll  $2,%r12d  ## gb
    shll  $3,%r14d  ## lj

    ## GB index: r12   LJ indices: r14

    ## calculate eps
    subss     %xmm6,%xmm4  ## gb
    subss     %xmm10,%xmm8 ## lj
    movaps    %xmm4,nb430_epsgb(%rsp)   ## gb eps
    movaps    %xmm8,nb430_eps(%rsp)   ## lj eps

        movq nb430_GBtab(%rbp),%rsi
        movq nb430_VFtab(%rbp),%rdi

    ## load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
    movss  (%rsi,%r12,4),%xmm0
    movss  4(%rsi,%r12,4),%xmm1
    movss  8(%rsi,%r12,4),%xmm2
    movss  12(%rsi,%r12,4),%xmm3
    movss  (%rdi,%r14,4),%xmm4
    movss  4(%rdi,%r14,4),%xmm5
    movss  8(%rdi,%r14,4),%xmm6
    movss  12(%rdi,%r14,4),%xmm7
    movss  16(%rdi,%r14,4),%xmm8
    movss  20(%rdi,%r14,4),%xmm9
    movss  24(%rdi,%r14,4),%xmm10
    movss  28(%rdi,%r14,4),%xmm11
    ## table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps nb430_epsgb(%rsp),%xmm12
    movaps nb430_eps(%rsp),%xmm13

    mulss  %xmm12,%xmm3  ## Heps
    mulss  %xmm13,%xmm7
    mulss  %xmm13,%xmm11
    mulss  %xmm12,%xmm2    ## Geps
    mulss  %xmm13,%xmm6
    mulss  %xmm13,%xmm10
    mulss  %xmm12,%xmm3  ## Heps2
    mulss  %xmm13,%xmm7
    mulss  %xmm13,%xmm11

    addss  %xmm2,%xmm1  ## F+Geps
    addss  %xmm6,%xmm5
    addss  %xmm10,%xmm9
    addss  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addss  %xmm7,%xmm5
    addss  %xmm11,%xmm9
    addss  %xmm3,%xmm3   ## 2*Heps2
    addss  %xmm7,%xmm7
    addss  %xmm11,%xmm11
    addss  %xmm2,%xmm3   ## 2*Heps2+Geps
    addss  %xmm6,%xmm7
    addss  %xmm10,%xmm11
    addss  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addss  %xmm5,%xmm7
    addss  %xmm9,%xmm11
    mulss  %xmm12,%xmm1  ## eps*Fp
    mulss  %xmm13,%xmm5
    mulss  %xmm13,%xmm9
    addss  %xmm0,%xmm1    ## VV
    addss  %xmm4,%xmm5
    addss  %xmm8,%xmm9
    mulss  nb430_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulss  nb430_c6(%rsp),%xmm5     ## vnb6
    mulss  nb430_c12(%rsp),%xmm9     ## vnb12
    mulss  nb430_qq(%rsp),%xmm3      ## FF*qq = fij
    mulss  nb430_c6(%rsp),%xmm7     ## fijD
    mulss  nb430_c12(%rsp),%xmm11     ##fijR

    addss  %xmm7,%xmm11 ## fijD+fijR
    mulss  nb430_tsc(%rsp),%xmm11   ## (fijD+fijR)*tabscale

    ## accumulate Vvdwtot
    addss  nb430_Vvdwtot(%rsp),%xmm5
    addss  %xmm9,%xmm5
    movss %xmm5,nb430_Vvdwtot(%rsp)

        movq nb430_dvda(%rbp),%rsi

        ## Calculate dVda
        mulss nb430_gbscale(%rsp),%xmm3     ## fijC=qq*FF*gbscale
        movaps %xmm3,%xmm6
        mulss  nb430_r(%rsp),%xmm6
        addss  %xmm1,%xmm6  ## vcoul+fijC*r

    addss  %xmm11,%xmm3 ## fijC+fijD+fijR

    ## increment vctot
        addss  nb430_vctot(%rsp),%xmm1
    movss %xmm1,nb430_vctot(%rsp)

        ## xmm6=(vcoul+fijC*r)
        xorps  %xmm7,%xmm7
        subss  %xmm6,%xmm7
        movaps %xmm7,%xmm6

        ## update dvdasum 
        addss  nb430_dvdasum(%rsp),%xmm7
    movss %xmm7,nb430_dvdasum(%rsp)

        ## update j atoms dvdaj

        ## xmm6=dvdaj1
        addss  (%rsi,%rax,4),%xmm6
        movss  %xmm6,(%rsi,%rax,4)

        xorps  %xmm4,%xmm4
        mulss nb430_rinv(%rsp),%xmm3
        subss  %xmm3,%xmm4

    movss  %xmm4,%xmm9
    movss  %xmm4,%xmm10
    movss  %xmm4,%xmm11

    mulss  nb430_dx(%rsp),%xmm9
    mulss  nb430_dy(%rsp),%xmm10
    mulss  nb430_dz(%rsp),%xmm11

        ## accumulate i forces
    movaps nb430_fix(%rsp),%xmm12
    movaps nb430_fiy(%rsp),%xmm13
    movaps nb430_fiz(%rsp),%xmm14
    addss %xmm9,%xmm12
    addss %xmm10,%xmm13
    addss %xmm11,%xmm14
    movss %xmm12,nb430_fix(%rsp)
    movss %xmm13,nb430_fiy(%rsp)
    movss %xmm14,nb430_fiz(%rsp)

        movq nb430_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%r8,4),%xmm9
    addss  4(%rsi,%r8,4),%xmm10
    addss  8(%rsi,%r8,4),%xmm11
    movss  %xmm9,(%rsi,%r8,4)
    movss  %xmm10,4(%rsi,%r8,4)
    movss  %xmm11,8(%rsi,%r8,4)

_nb_kernel430_x86_64_sse.nb430_updateouterdata: 
        movl  nb430_ii3(%rsp),%ecx
        movq  nb430_faction(%rbp),%rdi
        movq  nb430_fshift(%rbp),%rsi
        movl  nb430_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb430_fix(%rsp),%xmm0
        movaps nb430_fiy(%rsp),%xmm1
        movaps nb430_fiz(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%rdi,%rcx,4),%xmm3
        movss  4(%rdi,%rcx,4),%xmm4
        movss  8(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rdi,%rcx,4)
        movss  %xmm4,4(%rdi,%rcx,4)
        movss  %xmm5,8(%rdi,%rcx,4)

        ## increment fshift force  
        movss  (%rsi,%rdx,4),%xmm3
        movss  4(%rsi,%rdx,4),%xmm4
        movss  8(%rsi,%rdx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,4(%rsi,%rdx,4)
        movss  %xmm5,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb430_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb430_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb430_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb430_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb430_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb430_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate dVda and update it 
        movaps nb430_dvdasum(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        movl nb430_ii(%rsp),%edx
        movq nb430_dvda(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb430_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel430_x86_64_sse.nb430_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430_n(%rsp)
        jmp _nb_kernel430_x86_64_sse.nb430_outer
_nb_kernel430_x86_64_sse.nb430_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel430_x86_64_sse.nb430_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430_x86_64_sse.nb430_threadloop
_nb_kernel430_x86_64_sse.nb430_end: 
        movl nb430_nouter(%rsp),%eax
        movl nb430_ninner(%rsp),%ebx
        movq nb430_outeriter(%rbp),%rcx
        movq nb430_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $552,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel430nf_x86_64_sse
.globl _nb_kernel430nf_x86_64_sse
nb_kernel430nf_x86_64_sse:      
_nb_kernel430nf_x86_64_sse:     
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
        ## bottom of stack is cache-aligned for sse use 
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
.set nb430nf_isai, 208
.set nb430nf_isaprod, 224
.set nb430nf_gbscale, 240
.set nb430nf_r, 256
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
        movss (%rsi),%xmm0
        movss %xmm0,nb430nf_facel(%rsp)

        movq nb430nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb430nf_tsc(%rsp)

        movq nb430nf_p_gbtabscale(%rbp),%rbx
        movss (%rbx),%xmm4
        shufps $0,%xmm4,%xmm4
        movaps %xmm4,nb430nf_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb430nf_half(%rsp)
        movss nb430nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb430nf_half(%rsp)
        movaps %xmm3,nb430nf_three(%rsp)

_nb_kernel430nf_x86_64_sse.nb430nf_threadloop: 
        movq  nb430nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel430nf_x86_64_sse.nb430nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430nf_x86_64_sse.nb430nf_spinlock

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
        jg  _nb_kernel430nf_x86_64_sse.nb430nf_outerstart
        jmp _nb_kernel430nf_x86_64_sse.nb430nf_end

_nb_kernel430nf_x86_64_sse.nb430nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430nf_nouter(%rsp),%ebx
        movl %ebx,nb430nf_nouter(%rsp)

_nb_kernel430nf_x86_64_sse.nb430nf_outer: 
        movq  nb430nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb430nf_is3(%rsp)            ## store is3 

        movq  nb430nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb430nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb430nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb430nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb430nf_invsqrta(%rbp),%rdx       ## load invsqrta[ii]
        movss (%rdx,%rbx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        movq  nb430nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb430nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb430nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb430nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb430nf_iq(%rsp)
        movaps %xmm4,nb430nf_isai(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb430nf_ix(%rsp)
        movaps %xmm1,nb430nf_iy(%rsp)
        movaps %xmm2,nb430nf_iz(%rsp)

        movl  %ebx,nb430nf_ii3(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb430nf_vctot(%rsp)
        movaps %xmm4,nb430nf_Vvdwtot(%rsp)

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
        subl  $4,%edx
        addl  nb430nf_ninner(%rsp),%ecx
        movl  %ecx,nb430nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb430nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel430nf_x86_64_sse.nb430nf_unroll_loop
        jmp   _nb_kernel430nf_x86_64_sse.nb430nf_finish_inner
_nb_kernel430nf_x86_64_sse.nb430nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb430nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb430nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        ## load isa2
        movq nb430nf_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7
        movaps nb430nf_isai(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb430nf_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb430nf_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb430nf_gbscale(%rsp)

        movq nb430nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        mulps nb430nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430nf_qq(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movq nb430nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb430nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb430nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%rsi,%rax,4),%xmm6
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb430nf_c6(%rsp)
        movaps %xmm6,nb430nf_c12(%rsp)

        movq nb430nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   

        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ix-iz to xmm4-xmm6 
        movaps nb430nf_ix(%rsp),%xmm4
        movaps nb430nf_iy(%rsp),%xmm5
        movaps nb430nf_iz(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb430nf_r(%rsp)
        mulps nb430nf_gbscale(%rsp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movq nb430nf_GBtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## load coulomb table
        movaps (%rsi,%rax,4),%xmm4
        movaps (%rsi,%rbx,4),%xmm5
        movaps (%rsi,%rcx,4),%xmm6
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7            

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb430nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb430nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb430nf_vctot(%rsp)


        movaps nb430nf_r(%rsp),%xmm4
        mulps nb430nf_tsc(%rsp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6
        pslld $3,%mm7

        movq nb430nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movaps (%rsi,%rax,4),%xmm4
        movaps (%rsi,%rbx,4),%xmm5
        movaps (%rsi,%rcx,4),%xmm6
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  nb430nf_c6(%rsp),%xmm5    ## Vvdw6
        addps  nb430nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%rsp)

        ## repulsion 
        movaps 16(%rsi,%rax,4),%xmm4
        movaps 16(%rsi,%rbx,4),%xmm5
        movaps 16(%rsi,%rcx,4),%xmm6
        movaps 16(%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        mulps  nb430nf_c12(%rsp),%xmm5   ## Vvdw12
        addps  nb430nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb430nf_innerk(%rsp)
        jl    _nb_kernel430nf_x86_64_sse.nb430nf_finish_inner
        jmp   _nb_kernel430nf_x86_64_sse.nb430nf_unroll_loop
_nb_kernel430nf_x86_64_sse.nb430nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb430nf_innerk(%rsp)
        movl  nb430nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel430nf_x86_64_sse.nb430nf_dopair
        jmp   _nb_kernel430nf_x86_64_sse.nb430nf_checksingle
_nb_kernel430nf_x86_64_sse.nb430nf_dopair: 

        movq  nb430nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb430nf_innerjjnr(%rsp)

        xorps %xmm2,%xmm2
        movaps %xmm2,%xmm6

        ## load isa2
        movq nb430nf_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm2
        movss (%rsi,%rbx,4),%xmm3
        unpcklps %xmm3,%xmm2    ## isa2 in xmm3(0,1)
        mulps  nb430nf_isai(%rsp),%xmm2
        movaps %xmm2,nb430nf_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb430nf_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb430nf_gbscale(%rsp)

        movq nb430nf_charge(%rbp),%rsi     ## base of charge[]  
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        unpcklps %xmm6,%xmm3 ## 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb430nf_iq(%rsp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb430nf_qq(%rsp)

        movq nb430nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb430nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb430nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb430nf_pos(%rbp),%rdi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb430nf_c6(%rsp)
        movaps %xmm6,nb430nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        ## move coordinates to xmm0-xmm2 
        movlps (%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2
        movhps (%rdi,%rbx,4),%xmm1
        movss 8(%rdi,%rbx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## 10001000

        shufps $136,%xmm0,%xmm0 ## 10001000
        shufps $221,%xmm1,%xmm1 ## 11011101

        movq   nb430nf_faction(%rbp),%rdi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb430nf_ix(%rsp),%xmm4
        movaps nb430nf_iy(%rsp),%xmm5
        movaps nb430nf_iz(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb430nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb430nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm4,nb430nf_r(%rsp)
        mulps nb430nf_gbscale(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb430nf_GBtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## load coulomb table
        movaps (%rsi,%rcx,4),%xmm4
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb430nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb430nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb430nf_vctot(%rsp)

        movaps nb430nf_r(%rsp),%xmm4
        mulps nb430nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movq nb430nf_VFtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## dispersion 
        movaps (%rsi,%rcx,4),%xmm4
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        mulps  nb430nf_c6(%rsp),%xmm5    ## Vvdw6 
        addps  nb430nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%rsp)

        ## repulsion 
        movaps 16(%rsi,%rcx,4),%xmm4
        movaps 16(%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        mulps  nb430nf_c12(%rsp),%xmm5   ## Vvdw12 

        addps  nb430nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb430nf_Vvdwtot(%rsp)
_nb_kernel430nf_x86_64_sse.nb430nf_checksingle: 
        movl  nb430nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430nf_x86_64_sse.nb430nf_dosingle
        jmp    _nb_kernel430nf_x86_64_sse.nb430nf_updateouterdata
_nb_kernel430nf_x86_64_sse.nb430nf_dosingle: 
        movq nb430nf_charge(%rbp),%rsi
        movq nb430nf_invsqrta(%rbp),%rdx
        movq nb430nf_pos(%rbp),%rdi
        movq  nb430nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorps  %xmm2,%xmm2
        movaps %xmm2,%xmm6
        movss (%rdx,%rax,4),%xmm2       ## isa2
        mulss nb430nf_isai(%rsp),%xmm2
        movss %xmm2,nb430nf_isaprod(%rsp)
        movss %xmm2,%xmm1
        mulss nb430nf_gbtsc(%rsp),%xmm1
        movss %xmm1,nb430nf_gbscale(%rsp)

        mulss  nb430nf_iq(%rsp),%xmm2
        movss (%rsi,%rax,4),%xmm6       ## xmm6(0) has the charge       
        mulss  %xmm2,%xmm6
        movss %xmm6,nb430nf_qq(%rsp)

        movq nb430nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb430nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb430nf_ntia(%rsp),%ecx
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100     
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movss %xmm4,nb430nf_c6(%rsp)
        movss %xmm6,nb430nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        movss nb430nf_ix(%rsp),%xmm4
        movss nb430nf_iy(%rsp),%xmm5
        movss nb430nf_iz(%rsp),%xmm6

        ## calc dr 
        subss %xmm0,%xmm4
        subss %xmm1,%xmm5
        subss %xmm2,%xmm6

        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb430nf_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb430nf_half(%rsp),%xmm0
        subss %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 

        mulss %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm4,nb430nf_r(%rsp)
        mulss nb430nf_gbscale(%rsp),%xmm4

        cvttss2si %xmm4,%ebx    ## mm6 contain lu indices 
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%ebx

        movq nb430nf_GBtab(%rbp),%rsi

        movaps (%rsi,%rbx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        movss nb430nf_qq(%rsp),%xmm3
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        addss  nb430nf_vctot(%rsp),%xmm5
        movss %xmm5,nb430nf_vctot(%rsp)

        movss nb430nf_r(%rsp),%xmm4
        mulps nb430nf_tsc(%rsp),%xmm4

        cvttss2si %xmm4,%ebx
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx
        movq nb430nf_VFtab(%rbp),%rsi

        ## dispersion 
        movaps (%rsi,%rbx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  nb430nf_c6(%rsp),%xmm5    ## Vvdw6
        addss  nb430nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb430nf_Vvdwtot(%rsp)

        ## repulsion 
        movaps 16(%rsi,%rbx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        mulss  nb430nf_c12(%rsp),%xmm5   ## Vvdw12 

        addss  nb430nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb430nf_Vvdwtot(%rsp)

_nb_kernel430nf_x86_64_sse.nb430nf_updateouterdata: 
        ## get n from stack
        movl nb430nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb430nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb430nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb430nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb430nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb430nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb430nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel430nf_x86_64_sse.nb430nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430nf_n(%rsp)
        jmp _nb_kernel430nf_x86_64_sse.nb430nf_outer
_nb_kernel430nf_x86_64_sse.nb430nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel430nf_x86_64_sse.nb430nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430nf_x86_64_sse.nb430nf_threadloop
_nb_kernel430nf_x86_64_sse.nb430nf_end: 

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






