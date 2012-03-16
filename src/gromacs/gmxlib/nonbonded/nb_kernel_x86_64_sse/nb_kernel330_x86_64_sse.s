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







.globl nb_kernel330_x86_64_sse
.globl _nb_kernel330_x86_64_sse
nb_kernel330_x86_64_sse:        
_nb_kernel330_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb330_fshift, 16
.set nb330_gid, 24
.set nb330_pos, 32
.set nb330_faction, 40
.set nb330_charge, 48
.set nb330_p_facel, 56
.set nb330_argkrf, 64
.set nb330_argcrf, 72
.set nb330_Vc, 80
.set nb330_type, 88
.set nb330_p_ntype, 96
.set nb330_vdwparam, 104
.set nb330_Vvdw, 112
.set nb330_p_tabscale, 120
.set nb330_VFtab, 128
.set nb330_invsqrta, 136
.set nb330_dvda, 144
.set nb330_p_gbtabscale, 152
.set nb330_GBtab, 160
.set nb330_p_nthreads, 168
.set nb330_count, 176
.set nb330_mtx, 184
.set nb330_outeriter, 192
.set nb330_inneriter, 200
.set nb330_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb330_ix, 0
.set nb330_iy, 16
.set nb330_iz, 32
.set nb330_iq, 48
.set nb330_dx, 64
.set nb330_dy, 80
.set nb330_dz, 96
.set nb330_rinv, 112
.set nb330_tsc, 128
.set nb330_qq, 144
.set nb330_c6, 160
.set nb330_c12, 176
.set nb330_eps, 192
.set nb330_vctot, 208
.set nb330_Vvdwtot, 224
.set nb330_fix, 240
.set nb330_fiy, 256
.set nb330_fiz, 272
.set nb330_half, 288
.set nb330_three, 304
.set nb330_nri, 320
.set nb330_iinr, 328
.set nb330_jindex, 336
.set nb330_jjnr, 344
.set nb330_shift, 352
.set nb330_shiftvec, 360
.set nb330_facel, 368
.set nb330_innerjjnr, 376
.set nb330_is3, 384
.set nb330_ii3, 388
.set nb330_ntia, 392
.set nb330_innerk, 396
.set nb330_n, 400
.set nb330_nn1, 404
.set nb330_ntype, 408
.set nb330_nouter, 412
.set nb330_ninner, 416

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $440,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb330_nouter(%rsp)
        movl %eax,nb330_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb330_nri(%rsp)
        movq %rsi,nb330_iinr(%rsp)
        movq %rdx,nb330_jindex(%rsp)
        movq %rcx,nb330_jjnr(%rsp)
        movq %r8,nb330_shift(%rsp)
        movq %r9,nb330_shiftvec(%rsp)
        movq nb330_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb330_ntype(%rsp)
        movq nb330_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb330_facel(%rsp)

        movq nb330_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb330_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb330_half(%rsp)
        movss nb330_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb330_half(%rsp)
        movaps %xmm3,nb330_three(%rsp)

_nb_kernel330_x86_64_sse.nb330_threadloop: 
        movq  nb330_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel330_x86_64_sse.nb330_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330_x86_64_sse.nb330_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330_n(%rsp)
        movl %ebx,nb330_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330_x86_64_sse.nb330_outerstart
        jmp _nb_kernel330_x86_64_sse.nb330_end

_nb_kernel330_x86_64_sse.nb330_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330_nouter(%rsp),%ebx
        movl %ebx,nb330_nouter(%rsp)

_nb_kernel330_x86_64_sse.nb330_outer: 
        movq  nb330_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb330_is3(%rsp)      ## store is3 

        movq  nb330_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb330_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb330_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb330_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb330_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb330_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb330_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb330_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb330_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb330_ix(%rsp)
        movaps %xmm1,nb330_iy(%rsp)
        movaps %xmm2,nb330_iz(%rsp)

        movl  %ebx,nb330_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb330_vctot(%rsp)
        movaps %xmm4,nb330_Vvdwtot(%rsp)
        movaps %xmm4,nb330_fix(%rsp)
        movaps %xmm4,nb330_fiy(%rsp)
        movaps %xmm4,nb330_fiz(%rsp)

        movq  nb330_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl nb330_charge(%rsp),%esi
        movl nb330_ntia(%rsp),%edi

        movq  nb330_pos(%rbp),%rsi
        movq  nb330_faction(%rbp),%rdi
        movq  nb330_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb330_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb330_ninner(%rsp),%ecx
        movl  %ecx,nb330_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb330_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel330_x86_64_sse.nb330_unroll_loop
        jmp   _nb_kernel330_x86_64_sse.nb330_finish_inner
_nb_kernel330_x86_64_sse.nb330_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb330_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r12d
        movl  4(%rdx),%r13d
        movl  8(%rdx),%r14d
        movl  12(%rdx),%r15d           ## eax-edx=jnr1-4 
        addq $16,nb330_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%r12,%r12,2),%rax     ## replace jnr with j3 
        lea  (%r13,%r13,2),%rbx
        lea  (%r14,%r14,2),%rcx
        lea  (%r15,%r15,2),%rdx

        movq nb330_pos(%rbp),%rdi
        ## load coordinates
        movlps (%rdi,%rax,4),%xmm1      ## x1 y1 - - 
        movlps (%rdi,%rcx,4),%xmm2      ## x3 y3 - - 
        movhps (%rdi,%rbx,4),%xmm1      ## x2 y2 - -
        movhps (%rdi,%rdx,4),%xmm2      ## x4 y4 - -

        movss 8(%rdi,%rax,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%rcx,4),%xmm6      ## z2 - - - 
        movss 8(%rdi,%rbx,4),%xmm7      ## z3 - - - 
        movss 8(%rdi,%rdx,4),%xmm8      ## z4 - - - 
    movlhps %xmm7,%xmm5 ## jzOa  -  jzOb  -
    movlhps %xmm8,%xmm6 ## jzOc  -  jzOd -

    movaps %xmm1,%xmm4
    unpcklps %xmm2,%xmm1 ## jxa jxc jya jyc        
    unpckhps %xmm2,%xmm4 ## jxb jxd jyb jyd
    movaps %xmm1,%xmm2
    unpcklps %xmm4,%xmm1 ## x
    unpckhps %xmm4,%xmm2 ## y
    shufps  $136,%xmm6,%xmm5  ## 10001000 => jzH2a jzH2b jzH2c jzH2d
        movq nb330_charge(%rbp),%rsi

        ## calc dr  
        subps nb330_ix(%rsp),%xmm1
        subps nb330_iy(%rsp),%xmm2
        subps nb330_iz(%rsp),%xmm5

        ## store dr
    movaps %xmm1,nb330_dx(%rsp)
    movaps %xmm2,nb330_dy(%rsp)
    movaps %xmm5,nb330_dz(%rsp)

        movss (%rsi,%r12,4),%xmm0
        movss (%rsi,%r14,4),%xmm3
        movss (%rsi,%r13,4),%xmm4
        movss (%rsi,%r15,4),%xmm6

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

    unpcklps %xmm3,%xmm0 ## jqa jqc - -
    unpcklps %xmm6,%xmm4 ## jqb jqd - -
    unpcklps %xmm4,%xmm0 ## jqa jqb jqc jqd
        mulps nb330_iq(%rsp),%xmm0
    movaps %xmm0,nb330_qq(%rsp)

        movq nb330_type(%rbp),%rsi
    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm12
        mulps %xmm5,%xmm5
        movaps nb330_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm12,%xmm4
        mulps nb330_half(%rsp),%xmm4
        movaps %xmm4,%xmm15
        mulps  %xmm4,%xmm1
    ## xmm15=rinv
    ## xmm1=r

    ## vdw parameters
        movl (%rsi,%r12,4),%r12d
        movl (%rsi,%r13,4),%r13d
        movl (%rsi,%r14,4),%r14d
        movl (%rsi,%r15,4),%r15d

    mulps nb330_tsc(%rsp),%xmm1   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm1,%xmm5

        shll %r12d
        shll %r13d
        shll %r14d
        shll %r15d
    movl nb330_ntia(%rsp),%edi

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm4

    ## multiply by 4
    pslld   $2,%xmm5

    ## multiply by three (copy, mult. by two, add back)
    movaps  %xmm5,%xmm6
    pslld   $1,%xmm5
    paddd   %xmm6,%xmm5

    ## calculate eps
    subps     %xmm4,%xmm1

        addl %edi,%r12d
        addl %edi,%r13d
        addl %edi,%r14d
        addl %edi,%r15d

    ## move to integer registers
    movhlps %xmm5,%xmm6
    movd    %xmm5,%r8d
    movd    %xmm6,%r10d
    pshufd $1,%xmm5,%xmm5
    pshufd $1,%xmm6,%xmm6
    movd    %xmm5,%r9d
    movd    %xmm6,%r11d

        movq nb330_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movlps (%rsi,%r14,4),%xmm7
        movhps (%rsi,%r13,4),%xmm3
        movhps (%rsi,%r15,4),%xmm7

        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb330_c6(%rsp)
    movaps %xmm3,nb330_c12(%rsp)

    movaps %xmm1,nb330_eps(%rsp)
    ## xmm15=rinv

        movq nb330_VFtab(%rbp),%rsi
    ## load Coulomb and LJ table data in parallel
        movlps (%rsi,%r8,4),%xmm1        ## Y1c F1c 
        movlps 16(%rsi,%r8,4),%xmm5      ## Y1d F1d 
        movlps 32(%rsi,%r8,4),%xmm9      ## Y1r F1r 

        movlps (%rsi,%r10,4),%xmm3       ## Y3c F3c 
        movlps 16(%rsi,%r10,4),%xmm7     ## Y3d F3d 
        movlps 32(%rsi,%r10,4),%xmm11    ## Y3r F3r 

        movhps (%rsi,%r9,4),%xmm1        ## Y1c F1c Y2c F2c
        movhps 16(%rsi,%r9,4),%xmm5      ## Y1d F1d Y2d F2d
        movhps 32(%rsi,%r9,4),%xmm9      ## Y1r F1r Y2r F2r

        movhps (%rsi,%r11,4),%xmm3       ## Y3c F3c Y4c F4c
        movhps 16(%rsi,%r11,4),%xmm7     ## Y3d F3d Y4d F4d
        movhps 32(%rsi,%r11,4),%xmm11    ## Y3r F3r Y4r F4r

    movaps %xmm1,%xmm0
    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm3,%xmm0 ## 10001000   => Y1c Y2c Y3c Y4c
        shufps $136,%xmm7,%xmm4 ## 10001000   => Y1d Y2d Y3d Y4d
        shufps $136,%xmm11,%xmm8 ## 10001000  => Y1r Y2r Y3r Y4r
        shufps $221,%xmm3,%xmm1 ## 11011101   => F1c F2c F3c F4c
        shufps $221,%xmm7,%xmm5 ## 11011101   => F1d F2d F3d F4d
        shufps $221,%xmm11,%xmm9 ## 11011101  => F1r F2r F3r F4r

        movlps 8(%rsi,%r8,4),%xmm3        ## G1c H1c 
        movlps 24(%rsi,%r8,4),%xmm7       ## G1d H1d 
        movlps 40(%rsi,%r8,4),%xmm11      ## G1r H1r 

        movlps 8(%rsi,%r10,4),%xmm12      ## G3c H3c 
        movlps 24(%rsi,%r10,4),%xmm13     ## G3d H3d 
        movlps 40(%rsi,%r10,4),%xmm14     ## G3r H3r 

        movhps 8(%rsi,%r9,4),%xmm3        ## G1c H1c G2c H2c
        movhps 24(%rsi,%r9,4),%xmm7       ## G1d H1d G2d H2d
        movhps 40(%rsi,%r9,4),%xmm11      ## G1r H1r G2r H2r

        movhps 8(%rsi,%r11,4),%xmm12      ## G3c H3c G4c H4c
        movhps 24(%rsi,%r11,4),%xmm13     ## G3d H3d G4d H4d
        movhps 40(%rsi,%r11,4),%xmm14     ## G3r H3r G4r H4r


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

    movaps nb330_eps(%rsp),%xmm12

    mulps  %xmm12,%xmm3  ## Heps
    mulps  %xmm12,%xmm7
    mulps  %xmm12,%xmm11
    mulps  %xmm12,%xmm2    ## Geps
    mulps  %xmm12,%xmm6
    mulps  %xmm12,%xmm10
    mulps  %xmm12,%xmm3  ## Heps2
    mulps  %xmm12,%xmm7
    mulps  %xmm12,%xmm11

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
    mulps  %xmm12,%xmm5
    mulps  %xmm12,%xmm9
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  nb330_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulps  nb330_c6(%rsp),%xmm5     ## vnb6
    mulps  nb330_c12(%rsp),%xmm9     ## vnb12
    mulps  nb330_qq(%rsp),%xmm3      ## FF*qq = fij
    mulps  nb330_c6(%rsp),%xmm7     ## fijD
    mulps  nb330_c12(%rsp),%xmm11     ##fijR

    ## accumulate vctot
    addps  nb330_vctot(%rsp),%xmm1

    ## accumulate Vvdwtot
    addps  nb330_Vvdwtot(%rsp),%xmm5
    addps  %xmm9,%xmm5
    xorps  %xmm9,%xmm9

    movaps %xmm1,nb330_vctot(%rsp)
    movaps %xmm5,nb330_Vvdwtot(%rsp)


        movq nb330_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%rax,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%rcx,4),%xmm1 ## x3 y3 - -
    addps  %xmm7,%xmm3
    addps  %xmm11,%xmm3
    mulps  %xmm15,%xmm3
    mulps  nb330_tsc(%rsp),%xmm3    ## fscal

    subps  %xmm3,%xmm9
    movaps %xmm9,%xmm10
    movaps %xmm9,%xmm11

        movhps (%rsi,%rbx,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%rdx,4),%xmm1 ## x3 y3 x4 y4

    movaps nb330_fix(%rsp),%xmm12
    movaps nb330_fiy(%rsp),%xmm13
    movaps nb330_fiz(%rsp),%xmm14

    mulps  nb330_dx(%rsp),%xmm9
    mulps  nb330_dy(%rsp),%xmm10
    mulps  nb330_dz(%rsp),%xmm11

    ## accumulate i forces
    addps %xmm9,%xmm12
    addps %xmm10,%xmm13
    addps %xmm11,%xmm14
    movaps %xmm12,nb330_fix(%rsp)
    movaps %xmm13,nb330_fiy(%rsp)
    movaps %xmm14,nb330_fiz(%rsp)

    movaps %xmm9,%xmm8
    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    unpckhps %xmm10,%xmm8 ## x3 y3 x4 y4

    ## update fjx and fjy
        addps  %xmm9,%xmm0
        addps  %xmm8,%xmm1

        movlps %xmm0,(%rsi,%rax,4)
        movlps %xmm1,(%rsi,%rcx,4)
        movhps %xmm0,(%rsi,%rbx,4)
        movhps %xmm1,(%rsi,%rdx,4)

    ## xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd $1,%xmm11,%xmm10 ## fjz2 - - -
    movhlps %xmm11,%xmm9     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm8  ## fjz4 - - -

        addss  8(%rsi,%rax,4),%xmm11
        addss  8(%rsi,%rbx,4),%xmm10
        addss  8(%rsi,%rcx,4),%xmm9
        addss  8(%rsi,%rdx,4),%xmm8
        movss  %xmm11,8(%rsi,%rax,4)
        movss  %xmm10,8(%rsi,%rbx,4)
        movss  %xmm9,8(%rsi,%rcx,4)
        movss  %xmm8,8(%rsi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb330_innerk(%rsp)
        jl    _nb_kernel330_x86_64_sse.nb330_finish_inner
        jmp   _nb_kernel330_x86_64_sse.nb330_unroll_loop
_nb_kernel330_x86_64_sse.nb330_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb330_innerk(%rsp)
        movl  nb330_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel330_x86_64_sse.nb330_dopair
        jmp   _nb_kernel330_x86_64_sse.nb330_checksingle
_nb_kernel330_x86_64_sse.nb330_dopair: 
    movq  nb330_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb330_innerjjnr(%rsp)

        movq nb330_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0
        movss (%rsi,%rbx,4),%xmm2

    unpcklps %xmm2,%xmm0 ## jqa jqb 
        mulps nb330_iq(%rsp),%xmm0
    movaps %xmm0,nb330_qq(%rsp)

        movq nb330_type(%rbp),%rsi
    ## vdw parameters
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        shll %r12d
        shll %r13d
    movl nb330_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb330_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps  %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb330_c6(%rsp)
    movaps %xmm3,nb330_c12(%rsp)

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load coordinates
        movq nb330_pos(%rbp),%rdi

        movlps (%rdi,%rax,4),%xmm1      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm2      ## x2 y2 - - 

        movss 8(%rdi,%rax,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm6      ## z2 - - - 

    unpcklps %xmm2,%xmm1 ## x1 x2 y1 y2
    movhlps  %xmm1,%xmm2 ## y1 y2 -  -
    unpcklps %xmm6,%xmm5 ## z1 z2 -  -

        ## calc dr  
        subps nb330_ix(%rsp),%xmm1
        subps nb330_iy(%rsp),%xmm2
        subps nb330_iz(%rsp),%xmm5

        ## store dr
    movaps %xmm1,nb330_dx(%rsp)
    movaps %xmm2,nb330_dy(%rsp)
    movaps %xmm5,nb330_dz(%rsp)

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1

        ## rsq in xmm1

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb330_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm4
        mulps nb330_half(%rsp),%xmm4
        movaps %xmm4,%xmm15
        mulps  %xmm4,%xmm1
    ## xmm15=rinv
    ## xmm1=r

    mulps nb330_tsc(%rsp),%xmm1   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm1,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm4

    ## multiply by 4
    pslld   $2,%xmm5

    ## multiply by three (copy, mult. by two, add back)
    movaps  %xmm5,%xmm6
    pslld   $1,%xmm5
    paddd   %xmm6,%xmm5

    ## calculate eps
    subps     %xmm4,%xmm1

    ## move to integer registers
    movd    %xmm5,%r8d
    pshufd $1,%xmm5,%xmm5
    movd    %xmm5,%r9d

    movaps %xmm1,nb330_eps(%rsp)
    ## xmm15=rinv

        movq nb330_VFtab(%rbp),%rsi
    ## load Coulomb and LJ table data in parallel
        movlps (%rsi,%r8,4),%xmm0        ## Y1c F1c 
        movlps 16(%rsi,%r8,4),%xmm4      ## Y1d F1d 
        movlps 32(%rsi,%r8,4),%xmm8      ## Y1r F1r 

        movlps (%rsi,%r9,4),%xmm1        ## Y2c F2c
        movlps 16(%rsi,%r9,4),%xmm5      ## Y2d F2d
        movlps 32(%rsi,%r9,4),%xmm9     ## Y2r F2r

    unpcklps %xmm1,%xmm0 ## Y1c Y2c F1c F2c
    unpcklps %xmm5,%xmm4 ## Y1d Y2d F1d F2d
    unpcklps %xmm9,%xmm8 ## Y1r Y2r F1r F2r
    movhlps  %xmm0,%xmm1 ## F1c F2c
    movhlps  %xmm4,%xmm5 ## F1d F2d
    movhlps  %xmm8,%xmm9 ## F1r F2r

        movlps 8(%rsi,%r8,4),%xmm2       ## G1c H1c 
        movlps 24(%rsi,%r8,4),%xmm6      ## G1d H1d 
        movlps 40(%rsi,%r8,4),%xmm10     ## G1r H1r 

        movlps 8(%rsi,%r9,4),%xmm3       ## G2c H2c
        movlps 24(%rsi,%r9,4),%xmm7      ## G2d H2d
        movlps 40(%rsi,%r9,4),%xmm11     ## G2r H2r

    unpcklps %xmm3,%xmm2  ## G1c G2c H1c H2c
    unpcklps %xmm7,%xmm6  ## G1d G2d H1d H2d
    unpcklps %xmm11,%xmm10 ## G1r G2r H1r H2r
    movhlps  %xmm2,%xmm3  ## H1c H2c
    movhlps  %xmm6,%xmm7  ## H1d H2d
    movhlps  %xmm10,%xmm11 ## H1r H2r
    ## table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps nb330_eps(%rsp),%xmm12

    mulps  %xmm12,%xmm3  ## Heps
    mulps  %xmm12,%xmm7
    mulps  %xmm12,%xmm11
    mulps  %xmm12,%xmm2    ## Geps
    mulps  %xmm12,%xmm6
    mulps  %xmm12,%xmm10
    mulps  %xmm12,%xmm3  ## Heps2
    mulps  %xmm12,%xmm7
    mulps  %xmm12,%xmm11

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
    mulps  %xmm12,%xmm5
    mulps  %xmm12,%xmm9
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  nb330_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulps  nb330_c6(%rsp),%xmm5     ## vnb6
    mulps  nb330_c12(%rsp),%xmm9     ## vnb12
    mulps  nb330_qq(%rsp),%xmm3      ## FF*qq = fij
    mulps  nb330_c6(%rsp),%xmm7     ## fijD
    mulps  nb330_c12(%rsp),%xmm11     ##fijR

    ## accumulate vctot
    addps  nb330_vctot(%rsp),%xmm1
    movlps %xmm1,nb330_vctot(%rsp)

    ## accumulate Vvdwtot
    addps  nb330_Vvdwtot(%rsp),%xmm5
    addps  %xmm9,%xmm5
    movlps %xmm5,nb330_Vvdwtot(%rsp)

    xorps  %xmm9,%xmm9

    addps  %xmm7,%xmm3
    addps  %xmm11,%xmm3
    mulps  %xmm15,%xmm3
    mulps  nb330_tsc(%rsp),%xmm3    ## fscal

    subps  %xmm3,%xmm9
    movaps %xmm9,%xmm10
    movaps %xmm9,%xmm11

    movaps nb330_fix(%rsp),%xmm12
    movaps nb330_fiy(%rsp),%xmm13
    movaps nb330_fiz(%rsp),%xmm14

    mulps  nb330_dx(%rsp),%xmm9
    mulps  nb330_dy(%rsp),%xmm10
    mulps  nb330_dz(%rsp),%xmm11

    ## accumulate i forces
    addps %xmm9,%xmm12
    addps %xmm10,%xmm13
    addps %xmm11,%xmm14
    movlps %xmm12,nb330_fix(%rsp)
    movlps %xmm13,nb330_fiy(%rsp)
    movlps %xmm14,nb330_fiz(%rsp)

        movq nb330_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%rax,4),%xmm0 ## x1 y1 - -
        movhps (%rsi,%rbx,4),%xmm0 ## x1 y1 x2 y2

    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    addps    %xmm9,%xmm0

        movlps %xmm0,(%rsi,%rax,4)
        movhps %xmm0,(%rsi,%rbx,4)

    ## z forces
    pshufd $1,%xmm11,%xmm8
    addss  8(%rsi,%rax,4),%xmm11
    addss  8(%rsi,%rbx,4),%xmm8
    movss  %xmm11,8(%rsi,%rax,4)
    movss  %xmm8,8(%rsi,%rbx,4)

_nb_kernel330_x86_64_sse.nb330_checksingle:     
        movl  nb330_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330_x86_64_sse.nb330_dosingle
        jmp    _nb_kernel330_x86_64_sse.nb330_updateouterdata
_nb_kernel330_x86_64_sse.nb330_dosingle: 
        movq nb330_pos(%rbp),%rdi
        movq  nb330_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb330_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0

        mulss nb330_iq(%rsp),%xmm0
    movaps %xmm0,nb330_qq(%rsp)

    ## vdw parameters
        movq nb330_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb330_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb330_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
        movss 4(%rsi,%r12,4),%xmm3

    movaps %xmm0,nb330_c6(%rsp)
    movaps %xmm3,nb330_c12(%rsp)

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        movq nb330_pos(%rbp),%rdi
        ## load coordinates
        movss (%rdi,%rax,4),%xmm1
        movss 4(%rdi,%rax,4),%xmm2
        movss 8(%rdi,%rax,4),%xmm5

        ## calc dr  
        subss nb330_ix(%rsp),%xmm1
        subss nb330_iy(%rsp),%xmm2
        subss nb330_iz(%rsp),%xmm5

        ## store dr
    movaps %xmm1,nb330_dx(%rsp)
    movaps %xmm2,nb330_dy(%rsp)
    movaps %xmm5,nb330_dz(%rsp)

        ## square it 
        mulss %xmm1,%xmm1
        mulss %xmm2,%xmm2
        mulss %xmm5,%xmm5
        addss %xmm2,%xmm1
        addss %xmm5,%xmm1

        ## rsq in xmm1

    ## calculate rinv=1/sqrt(rsq)
        rsqrtss %xmm1,%xmm5
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movaps nb330_three(%rsp),%xmm4
        mulss %xmm1,%xmm5       ## rsq*lu*lu    
    subss %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm4
        mulss nb330_half(%rsp),%xmm4
        movaps %xmm4,%xmm15
        mulss  %xmm4,%xmm1
    ## xmm15=rinv
    ## xmm1=r

    mulss nb330_tsc(%rsp),%xmm1   ## rtab

    ## truncate and convert to integers
    cvttss2si %xmm1,%r8d

    ## convert back to float
    cvtsi2ss  %r8d,%xmm4

    ## multiply by 4
    shll      $2,%r8d

    ## calculate eps
    subss     %xmm4,%xmm1

    ## mult. by 3
        lea  (%r8,%r8,2),%r8

    movaps %xmm1,nb330_eps(%rsp)
    ## xmm15=rinv

        movq nb330_VFtab(%rbp),%rsi
    ## load Coulomb and LJ table data in parallel
    movss  (%rsi,%r8,4),%xmm0
    movss  4(%rsi,%r8,4),%xmm1
    movss  8(%rsi,%r8,4),%xmm2
    movss  12(%rsi,%r8,4),%xmm3
    movss  16(%rsi,%r8,4),%xmm4
    movss  20(%rsi,%r8,4),%xmm5
    movss  24(%rsi,%r8,4),%xmm6
    movss  28(%rsi,%r8,4),%xmm7
    movss  32(%rsi,%r8,4),%xmm8
    movss  36(%rsi,%r8,4),%xmm9
    movss  40(%rsi,%r8,4),%xmm10
    movss  44(%rsi,%r8,4),%xmm11
    ## table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps nb330_eps(%rsp),%xmm12

    mulps  %xmm12,%xmm3  ## Heps
    mulps  %xmm12,%xmm7
    mulps  %xmm12,%xmm11
    mulps  %xmm12,%xmm2    ## Geps
    mulps  %xmm12,%xmm6
    mulps  %xmm12,%xmm10
    mulps  %xmm12,%xmm3  ## Heps2
    mulps  %xmm12,%xmm7
    mulps  %xmm12,%xmm11

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
    mulps  nb330_eps(%rsp),%xmm1     ## eps*Fp
    mulps  nb330_eps(%rsp),%xmm5
    mulps  nb330_eps(%rsp),%xmm9
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  nb330_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulps  nb330_c6(%rsp),%xmm5     ## vnb6
    mulps  nb330_c12(%rsp),%xmm9     ## vnb12
    mulps  nb330_qq(%rsp),%xmm3      ## FF*qq = fij
    mulps  nb330_c6(%rsp),%xmm7     ## fijD
    mulps  nb330_c12(%rsp),%xmm11     ##fijR

    ## accumulate vctot
    addps  nb330_vctot(%rsp),%xmm1
    movaps %xmm1,nb330_vctot(%rsp)

    ## accumulate Vvdwtot
    addps  nb330_Vvdwtot(%rsp),%xmm5
    addps  %xmm9,%xmm5
    movaps %xmm5,nb330_Vvdwtot(%rsp)

    xorps  %xmm9,%xmm9

    addps  %xmm7,%xmm3
    addps  %xmm11,%xmm3
    mulss  %xmm15,%xmm3
    mulps  nb330_tsc(%rsp),%xmm3    ## fscal

    subps  %xmm3,%xmm9
    movaps %xmm9,%xmm10
    movaps %xmm9,%xmm11

    movaps nb330_fix(%rsp),%xmm12
    movaps nb330_fiy(%rsp),%xmm13
    movaps nb330_fiz(%rsp),%xmm14

    mulss  nb330_dx(%rsp),%xmm9
    mulss  nb330_dy(%rsp),%xmm10
    mulss  nb330_dz(%rsp),%xmm11

    ## accumulate i forces
    addss %xmm9,%xmm12
    addss %xmm10,%xmm13
    addss %xmm11,%xmm14
    movaps %xmm12,nb330_fix(%rsp)
    movaps %xmm13,nb330_fiy(%rsp)
    movaps %xmm14,nb330_fiz(%rsp)

        movq nb330_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel330_x86_64_sse.nb330_updateouterdata: 
        movl  nb330_ii3(%rsp),%ecx
        movq  nb330_faction(%rbp),%rdi
        movq  nb330_fshift(%rbp),%rsi
        movl  nb330_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb330_fix(%rsp),%xmm0
        movaps nb330_fiy(%rsp),%xmm1
        movaps nb330_fiz(%rsp),%xmm2

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
        movl nb330_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb330_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb330_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb330_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb330_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb330_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb330_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330_x86_64_sse.nb330_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330_n(%rsp)
        jmp _nb_kernel330_x86_64_sse.nb330_outer
_nb_kernel330_x86_64_sse.nb330_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330_x86_64_sse.nb330_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330_x86_64_sse.nb330_threadloop
_nb_kernel330_x86_64_sse.nb330_end: 

        movl nb330_nouter(%rsp),%eax
        movl nb330_ninner(%rsp),%ebx
        movq nb330_outeriter(%rbp),%rcx
        movq nb330_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $440,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel330nf_x86_64_sse
.globl _nb_kernel330nf_x86_64_sse
nb_kernel330nf_x86_64_sse:      
_nb_kernel330nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb330nf_fshift, 16
.set nb330nf_gid, 24
.set nb330nf_pos, 32
.set nb330nf_faction, 40
.set nb330nf_charge, 48
.set nb330nf_p_facel, 56
.set nb330nf_argkrf, 64
.set nb330nf_argcrf, 72
.set nb330nf_Vc, 80
.set nb330nf_type, 88
.set nb330nf_p_ntype, 96
.set nb330nf_vdwparam, 104
.set nb330nf_Vvdw, 112
.set nb330nf_p_tabscale, 120
.set nb330nf_VFtab, 128
.set nb330nf_invsqrta, 136
.set nb330nf_dvda, 144
.set nb330nf_p_gbtabscale, 152
.set nb330nf_GBtab, 160
.set nb330nf_p_nthreads, 168
.set nb330nf_count, 176
.set nb330nf_mtx, 184
.set nb330nf_outeriter, 192
.set nb330nf_inneriter, 200
.set nb330nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb330nf_ix, 0
.set nb330nf_iy, 16
.set nb330nf_iz, 32
.set nb330nf_iq, 48
.set nb330nf_tsc, 64
.set nb330nf_qq, 80
.set nb330nf_c6, 96
.set nb330nf_c12, 112
.set nb330nf_vctot, 128
.set nb330nf_Vvdwtot, 144
.set nb330nf_half, 160
.set nb330nf_three, 176
.set nb330nf_nri, 192
.set nb330nf_iinr, 200
.set nb330nf_jindex, 208
.set nb330nf_jjnr, 216
.set nb330nf_shift, 224
.set nb330nf_shiftvec, 232
.set nb330nf_facel, 240
.set nb330nf_innerjjnr, 248
.set nb330nf_is3, 256
.set nb330nf_ii3, 260
.set nb330nf_ntia, 264
.set nb330nf_innerk, 268
.set nb330nf_n, 272
.set nb330nf_nn1, 276
.set nb330nf_ntype, 280
.set nb330nf_nouter, 284
.set nb330nf_ninner, 288

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $312,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb330nf_nouter(%rsp)
        movl %eax,nb330nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb330nf_nri(%rsp)
        movq %rsi,nb330nf_iinr(%rsp)
        movq %rdx,nb330nf_jindex(%rsp)
        movq %rcx,nb330nf_jjnr(%rsp)
        movq %r8,nb330nf_shift(%rsp)
        movq %r9,nb330nf_shiftvec(%rsp)
        movq nb330nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb330nf_ntype(%rsp)
        movq nb330nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb330nf_facel(%rsp)


        movq nb330nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb330nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb330nf_half(%rsp)
        movss nb330nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb330nf_half(%rsp)
        movaps %xmm3,nb330nf_three(%rsp)

_nb_kernel330nf_x86_64_sse.nb330nf_threadloop: 
        movq  nb330nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel330nf_x86_64_sse.nb330nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330nf_x86_64_sse.nb330nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330nf_n(%rsp)
        movl %ebx,nb330nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330nf_x86_64_sse.nb330nf_outerstart
        jmp _nb_kernel330nf_x86_64_sse.nb330nf_end

_nb_kernel330nf_x86_64_sse.nb330nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330nf_nouter(%rsp),%ebx
        movl %ebx,nb330nf_nouter(%rsp)

_nb_kernel330nf_x86_64_sse.nb330nf_outer: 
        movq  nb330nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb330nf_is3(%rsp)            ## store is3 

        movq  nb330nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb330nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb330nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb330nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb330nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb330nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb330nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb330nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb330nf_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb330nf_ix(%rsp)
        movaps %xmm1,nb330nf_iy(%rsp)
        movaps %xmm2,nb330nf_iz(%rsp)

        movl  %ebx,nb330nf_ii3(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb330nf_vctot(%rsp)
        movaps %xmm4,nb330nf_Vvdwtot(%rsp)

        movq  nb330nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb330nf_pos(%rbp),%rsi
        movq  nb330nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb330nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb330nf_ninner(%rsp),%ecx
        movl  %ecx,nb330nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb330nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel330nf_x86_64_sse.nb330nf_unroll_loop
        jmp   _nb_kernel330nf_x86_64_sse.nb330nf_finish_inner
_nb_kernel330nf_x86_64_sse.nb330nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb330nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb330nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb330nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        movaps nb330nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        mulps  %xmm2,%xmm3
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm3,nb330nf_qq(%rsp)

        movq nb330nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb330nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb330nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb330nf_c6(%rsp)
        movaps %xmm6,nb330nf_c12(%rsp)

        movq nb330nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb330nf_ix(%rsp),%xmm4
        movaps nb330nf_iy(%rsp),%xmm5
        movaps nb330nf_iz(%rsp),%xmm6

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
        movaps nb330nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330nf_tsc(%rsp),%xmm4

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

        movq nb330nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movlps (%rsi,%rax,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rbx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb330nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb330nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb330nf_vctot(%rsp)

        ## dispersion 
        movlps 16(%rsi,%rax,4),%xmm5
        movlps 16(%rsi,%rcx,4),%xmm7
        movhps 16(%rsi,%rbx,4),%xmm5
        movhps 16(%rsi,%rdx,4),%xmm7    ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 24(%rsi,%rax,4),%xmm7
        movlps 24(%rsi,%rcx,4),%xmm3
        movhps 24(%rsi,%rbx,4),%xmm7
        movhps 24(%rsi,%rdx,4),%xmm3    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 
        ## put scalar force on stack 
        addps  nb330nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%rsp)

        ## repulsion 
        movlps 32(%rsi,%rax,4),%xmm5
        movlps 32(%rsi,%rcx,4),%xmm7
        movhps 32(%rsi,%rbx,4),%xmm5
        movhps 32(%rsi,%rdx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 40(%rsi,%rax,4),%xmm7
        movlps 40(%rsi,%rcx,4),%xmm3
        movhps 40(%rsi,%rbx,4),%xmm7
        movhps 40(%rsi,%rdx,4),%xmm3    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb330nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb330nf_innerk(%rsp)
        jl    _nb_kernel330nf_x86_64_sse.nb330nf_finish_inner
        jmp   _nb_kernel330nf_x86_64_sse.nb330nf_unroll_loop
_nb_kernel330nf_x86_64_sse.nb330nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb330nf_innerk(%rsp)
        movl  nb330nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel330nf_x86_64_sse.nb330nf_dopair
        jmp   _nb_kernel330nf_x86_64_sse.nb330nf_checksingle
_nb_kernel330nf_x86_64_sse.nb330nf_dopair: 
        movq nb330nf_charge(%rbp),%rsi

        movq  nb330nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb330nf_innerjjnr(%rsp)
        xorps %xmm7,%xmm7
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb330nf_iq(%rsp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb330nf_qq(%rsp)

        movq nb330nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb330nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb330nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb330nf_pos(%rbp),%rdi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb330nf_c6(%rsp)
        movaps %xmm6,nb330nf_c12(%rsp)

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


        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb330nf_ix(%rsp),%xmm4
        movaps nb330nf_iy(%rsp),%xmm5
        movaps nb330nf_iz(%rsp),%xmm6

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
        movaps nb330nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb330nf_VFtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movlps (%rsi,%rcx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm5 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rcx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## 10001000
        shufps $221,%xmm7,%xmm7 ## 11011101
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb330nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb330nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb330nf_vctot(%rsp)

        ## dispersion 
        movlps 16(%rsi,%rcx,4),%xmm5
        movhps 16(%rsi,%rdx,4),%xmm5   ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## 10001000
        shufps $221,%xmm5,%xmm5 ## 11011101

        movlps 24(%rsi,%rcx,4),%xmm7
        movhps 24(%rsi,%rdx,4),%xmm7    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## 10001000
        shufps $221,%xmm7,%xmm7 ## 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 
        ## put scalar force on stack 
        addps  nb330nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%rsp)

        ## repulsion 
        movlps 32(%rsi,%rcx,4),%xmm5
        movhps 32(%rsi,%rdx,4),%xmm5    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 40(%rsi,%rcx,4),%xmm7
        movhps 40(%rsi,%rdx,4),%xmm7    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb330nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%rsp)

_nb_kernel330nf_x86_64_sse.nb330nf_checksingle: 
        movl  nb330nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330nf_x86_64_sse.nb330nf_dosingle
        jmp    _nb_kernel330nf_x86_64_sse.nb330nf_updateouterdata
_nb_kernel330nf_x86_64_sse.nb330nf_dosingle: 
        movq nb330nf_charge(%rbp),%rsi
        movq nb330nf_pos(%rbp),%rdi
        movq  nb330nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorps  %xmm6,%xmm6
        movss (%rsi,%rax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb330nf_iq(%rsp),%xmm6
        movaps %xmm6,nb330nf_qq(%rsp)

        movq nb330nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb330nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb330nf_ntia(%rsp),%ecx
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100     
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movaps %xmm4,nb330nf_c6(%rsp)
        movaps %xmm6,nb330nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        movaps nb330nf_ix(%rsp),%xmm4
        movaps nb330nf_iy(%rsp),%xmm5
        movaps nb330nf_iz(%rsp),%xmm6

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
        movaps nb330nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb330nf_VFtab(%rbp),%rsi
        movd %mm6,%ebx

        lea (%rbx,%rbx,2),%rbx

        movlps (%rsi,%rbx,4),%xmm4
        movlps 8(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb330nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb330nf_vctot(%rsp),%xmm5
        movss %xmm5,nb330nf_vctot(%rsp)

        ## dispersion 
        movlps 16(%rsi,%rbx,4),%xmm4
        movlps 24(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack  
        addss  nb330nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb330nf_Vvdwtot(%rsp)

        ## repulsion 
        movlps 32(%rsi,%rbx,4),%xmm4
        movlps 40(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addss  nb330nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb330nf_Vvdwtot(%rsp)

_nb_kernel330nf_x86_64_sse.nb330nf_updateouterdata: 
        ## get n from stack
        movl nb330nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb330nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb330nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb330nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb330nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb330nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb330nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330nf_x86_64_sse.nb330nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330nf_n(%rsp)
        jmp _nb_kernel330nf_x86_64_sse.nb330nf_outer
_nb_kernel330nf_x86_64_sse.nb330nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330nf_x86_64_sse.nb330nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330nf_x86_64_sse.nb330nf_threadloop
_nb_kernel330nf_x86_64_sse.nb330nf_end: 

        movl nb330nf_nouter(%rsp),%eax
        movl nb330nf_ninner(%rsp),%ebx
        movq nb330nf_outeriter(%rbp),%rcx
        movq nb330nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $312,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





