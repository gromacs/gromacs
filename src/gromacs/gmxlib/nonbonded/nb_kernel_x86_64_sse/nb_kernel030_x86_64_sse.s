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





.globl nb_kernel030_x86_64_sse
.globl _nb_kernel030_x86_64_sse
nb_kernel030_x86_64_sse:        
_nb_kernel030_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb030_fshift, 16
.set nb030_gid, 24
.set nb030_pos, 32
.set nb030_faction, 40
.set nb030_charge, 48
.set nb030_p_facel, 56
.set nb030_argkrf, 64
.set nb030_argcrf, 72
.set nb030_Vc, 80
.set nb030_type, 88
.set nb030_p_ntype, 96
.set nb030_vdwparam, 104
.set nb030_Vvdw, 112
.set nb030_p_tabscale, 120
.set nb030_VFtab, 128
.set nb030_invsqrta, 136
.set nb030_dvda, 144
.set nb030_p_gbtabscale, 152
.set nb030_GBtab, 160
.set nb030_p_nthreads, 168
.set nb030_count, 176
.set nb030_mtx, 184
.set nb030_outeriter, 192
.set nb030_inneriter, 200
.set nb030_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb030_ix, 0
.set nb030_iy, 16
.set nb030_iz, 32
.set nb030_dx, 48
.set nb030_dy, 64
.set nb030_dz, 80
.set nb030_two, 96
.set nb030_tsc, 112
.set nb030_c6, 128
.set nb030_c12, 144
.set nb030_fscal, 160
.set nb030_Vvdwtot, 176
.set nb030_fix, 192
.set nb030_fiy, 208
.set nb030_fiz, 224
.set nb030_half, 240
.set nb030_three, 256
.set nb030_is3, 272
.set nb030_ii3, 276
.set nb030_nri, 280
.set nb030_iinr, 288
.set nb030_jindex, 296
.set nb030_jjnr, 304
.set nb030_shift, 312
.set nb030_shiftvec, 320
.set nb030_innerjjnr, 328
.set nb030_ntia, 336
.set nb030_innerk, 340
.set nb030_n, 344
.set nb030_nn1, 348
.set nb030_ntype, 352
.set nb030_nouter, 356
.set nb030_ninner, 360

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

    subq $376,%rsp              # # local variable stack space (n*16+8)
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb030_nouter(%rsp)
        movl %eax,nb030_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb030_nri(%rsp)
        movq %rsi,nb030_iinr(%rsp)
        movq %rdx,nb030_jindex(%rsp)
        movq %rcx,nb030_jjnr(%rsp)
        movq %r8,nb030_shift(%rsp)
        movq %r9,nb030_shiftvec(%rsp)
        movq nb030_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb030_ntype(%rsp)

        movq nb030_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb030_tsc(%rsp)

    movq nb030_pos(%rbp),%r8
    movq nb030_faction(%rbp),%r9
    movq nb030_type(%rbp),%r10
    movq nb030_vdwparam(%rbp),%r11


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb030_half(%rsp)
        movss nb030_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb030_half(%rsp)
        movaps %xmm2,nb030_two(%rsp)
        movaps %xmm3,nb030_three(%rsp)

_nb_kernel030_x86_64_sse.nb030_threadloop: 
        movq  nb030_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel030_x86_64_sse.nb030_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030_x86_64_sse.nb030_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030_n(%rsp)
        movl %ebx,nb030_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030_x86_64_sse.nb030_outerstart
        jmp _nb_kernel030_x86_64_sse.nb030_end

_nb_kernel030_x86_64_sse.nb030_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030_nouter(%rsp),%ebx
        movl %ebx,nb030_nouter(%rsp)

_nb_kernel030_x86_64_sse.nb030_outer: 
        movq  nb030_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb030_is3(%rsp)      ## store is3 

        movq  nb030_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb030_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb030_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb030_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb030_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb030_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb030_ix(%rsp)
        movaps %xmm1,nb030_iy(%rsp)
        movaps %xmm2,nb030_iz(%rsp)

        movl  %ebx,nb030_ii3(%rsp)

        ## clear tot potential and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb030_Vvdwtot(%rsp)
        movaps %xmm4,nb030_fix(%rsp)
        movaps %xmm4,nb030_fiy(%rsp)
        movaps %xmm4,nb030_fiz(%rsp)

        movq  nb030_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 


        movq  nb030_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb030_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb030_ninner(%rsp),%ecx
        movl  %ecx,nb030_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb030_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel030_x86_64_sse.nb030_unroll_loop
        jmp   _nb_kernel030_x86_64_sse.nb030_finish_inner
_nb_kernel030_x86_64_sse.nb030_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb030_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb030_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%rax,%rax,2),%r12     ##  j3 
        lea  (%rbx,%rbx,2),%r13
        lea  (%rcx,%rcx,2),%r14
        lea  (%rdx,%rdx,2),%r15

        movq nb030_pos(%rbp),%rdi
        ## load coordinates
        movlps (%rdi,%r12,4),%xmm1      ## x1 y1 - - 
        movlps (%rdi,%r14,4),%xmm2      ## x3 y3 - - 
        movhps (%rdi,%r13,4),%xmm1      ## x2 y2 - -
        movhps (%rdi,%r15,4),%xmm2      ## x4 y4 - -

        movss 8(%rdi,%r12,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%r14,4),%xmm6      ## z2 - - - 
        movss 8(%rdi,%r13,4),%xmm7      ## z3 - - - 
        movss 8(%rdi,%r15,4),%xmm8      ## z4 - - - 
    movlhps %xmm7,%xmm5 ## jzOa  -  jzOb  -
    movlhps %xmm8,%xmm6 ## jzOc  -  jzOd -

    movaps %xmm1,%xmm4
    unpcklps %xmm2,%xmm1 ## jxa jxc jya jyc        
    unpckhps %xmm2,%xmm4 ## jxb jxd jyb jyd
    movaps %xmm1,%xmm2
    unpcklps %xmm4,%xmm1 ## x
    unpckhps %xmm4,%xmm2 ## y
    shufps  $136,%xmm6,%xmm5  ## 10001000 => jzH2a jzH2b jzH2c jzH2d

        movq nb030_type(%rbp),%rsi

        ## calc dr  
        subps nb030_ix(%rsp),%xmm1
        subps nb030_iy(%rsp),%xmm2
        subps nb030_iz(%rsp),%xmm5

        ## store dr
    movaps %xmm1,nb030_dx(%rsp)
    movaps %xmm2,nb030_dy(%rsp)
    movaps %xmm5,nb030_dz(%rsp)


        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

    ## load vdw types
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb030_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm4
        mulps nb030_half(%rsp),%xmm4
        movaps %xmm4,%xmm2
        mulps  %xmm4,%xmm1
    ## xmm2=rinv
    ## xmm1=r

    movl nb030_ntia(%rsp),%edi
    ## type*=2
    shll %eax
        shll %ebx
        shll %ecx
        shll %edx

    mulps nb030_tsc(%rsp),%xmm1   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm1,%xmm5

        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm4

    ## multiply by 8
    pslld   $3,%xmm5

    ## calculate eps
    subps     %xmm4,%xmm1

    ## move to integer registers
    movhlps %xmm5,%xmm6
    movd    %xmm5,%r8d
    movd    %xmm6,%r10d
    pshufd $1,%xmm5,%xmm5
    pshufd $1,%xmm6,%xmm6
    movd    %xmm5,%r9d
    movd    %xmm6,%r11d

    ## xmm1=eps
    ## xmm2=rinv

        movq nb030_VFtab(%rbp),%rsi
    ## calculate LJ table
    movlps (%rsi,%r8,4),%xmm5
        movlps 16(%rsi,%r8,4),%xmm9

        movlps (%rsi,%r10,4),%xmm7
        movlps 16(%rsi,%r10,4),%xmm11

        movhps (%rsi,%r9,4),%xmm5
        movhps 16(%rsi,%r9,4),%xmm9

        movhps (%rsi,%r11,4),%xmm7
        movhps 16(%rsi,%r11,4),%xmm11

    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $136,%xmm11,%xmm8 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101
        shufps $221,%xmm11,%xmm9 ## 11011101

        movlps 8(%rsi,%r8,4),%xmm7
        movlps 24(%rsi,%r8,4),%xmm11

        movlps 8(%rsi,%r10,4),%xmm13
        movlps 24(%rsi,%r10,4),%xmm14

        movhps 8(%rsi,%r9,4),%xmm7
        movhps 24(%rsi,%r9,4),%xmm11

        movhps 8(%rsi,%r11,4),%xmm13
        movhps 24(%rsi,%r11,4),%xmm14

    movaps %xmm7,%xmm6
    movaps %xmm11,%xmm10

        shufps $136,%xmm13,%xmm6 ## 10001000
        shufps $136,%xmm14,%xmm10 ## 10001000
        shufps $221,%xmm13,%xmm7 ## 11011101
        shufps $221,%xmm14,%xmm11 ## 11011101
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

        movq nb030_vdwparam(%rbp),%rsi

    mulps  %xmm1,%xmm7   ## Heps
    mulps  %xmm1,%xmm11
    mulps  %xmm1,%xmm6  ## Geps
    mulps  %xmm1,%xmm10
    mulps  %xmm1,%xmm7  ## Heps2
    mulps  %xmm1,%xmm11

    ## load c6/c12
    movlps (%rsi,%rax,4),%xmm13
        movlps (%rsi,%rcx,4),%xmm14

    addps  %xmm6,%xmm5 ## F+Geps
    addps  %xmm10,%xmm9
    addps  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addps  %xmm11,%xmm9

        movhps (%rsi,%rbx,4),%xmm13
        movhps (%rsi,%rdx,4),%xmm14

    addps  %xmm7,%xmm7   ## 2*Heps2
    addps  %xmm11,%xmm11
    addps  %xmm6,%xmm7  ## 2*Heps2+Geps
    addps  %xmm10,%xmm11

    ## shuffle c6/c12
    movaps %xmm13,%xmm12
        shufps $136,%xmm14,%xmm12 ## 10001000  => c6
        shufps $221,%xmm14,%xmm13 ## 11011101  => c12

    addps  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addps  %xmm9,%xmm11
    mulps  %xmm1,%xmm5 ## eps*Fp
    mulps  %xmm1,%xmm9
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

        movq nb030_faction(%rbp),%rsi

    mulps  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulps  %xmm13,%xmm9 ## VV*c12 = vnb12
    addps  %xmm9,%xmm5
    addps  nb030_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb030_Vvdwtot(%rsp)

    mulps  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulps  %xmm13,%xmm11  ## FF*c12  = fnb12
    addps  %xmm11,%xmm7

    mulps  nb030_tsc(%rsp),%xmm7
    mulps  %xmm2,%xmm7
    xorps  %xmm9,%xmm9

    subps  %xmm7,%xmm9
    movaps %xmm9,%xmm10
    movaps %xmm9,%xmm11

    movaps nb030_fix(%rsp),%xmm12
    movaps nb030_fiy(%rsp),%xmm13
    movaps nb030_fiz(%rsp),%xmm14

    mulps  nb030_dx(%rsp),%xmm9
    mulps  nb030_dy(%rsp),%xmm10
    mulps  nb030_dz(%rsp),%xmm11

    ## accumulate i forces
    addps %xmm9,%xmm12
    addps %xmm10,%xmm13
    addps %xmm11,%xmm14
    movaps %xmm12,nb030_fix(%rsp)
    movaps %xmm13,nb030_fiy(%rsp)
    movaps %xmm14,nb030_fiz(%rsp)

        ## the fj's - start by combining x & y forces from memory 
        movlps (%rsi,%r12,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%r14,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%r13,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%r15,4),%xmm1 ## x3 y3 x4 y4

    movaps %xmm9,%xmm8
    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    unpckhps %xmm10,%xmm8 ## x3 y3 x4 y4

    ## update fjx and fjy
        addps  %xmm9,%xmm0
        addps  %xmm8,%xmm1

        movlps %xmm0,(%rsi,%r12,4)
        movlps %xmm1,(%rsi,%r14,4)
        movhps %xmm0,(%rsi,%r13,4)
        movhps %xmm1,(%rsi,%r15,4)

    ## xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd $1,%xmm11,%xmm5 ## fjz2 - - -
    movhlps %xmm11,%xmm4     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm3  ## fjz4 - - -

        addss  8(%rsi,%r12,4),%xmm11
        addss  8(%rsi,%r13,4),%xmm5
        addss  8(%rsi,%r14,4),%xmm4
        addss  8(%rsi,%r15,4),%xmm3
        movss  %xmm11,8(%rsi,%r12,4)
        movss  %xmm5,8(%rsi,%r13,4)
        movss  %xmm4,8(%rsi,%r14,4)
        movss  %xmm3,8(%rsi,%r15,4)

        ## should we do one more iteration? 
        subl $4,nb030_innerk(%rsp)
        jl    _nb_kernel030_x86_64_sse.nb030_finish_inner
        jmp   _nb_kernel030_x86_64_sse.nb030_unroll_loop
_nb_kernel030_x86_64_sse.nb030_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb030_innerk(%rsp)
        movl  nb030_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel030_x86_64_sse.nb030_dopair
        jmp   _nb_kernel030_x86_64_sse.nb030_checksingle
_nb_kernel030_x86_64_sse.nb030_dopair: 
    movq  nb030_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb030_innerjjnr(%rsp)

        movq nb030_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        shll %r12d
        shll %r13d
    movl nb030_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb030_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps  %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb030_c6(%rsp)
    movaps %xmm3,nb030_c12(%rsp)

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        movq nb030_pos(%rbp),%rdi
        ## load coordinates
        movlps (%rdi,%rax,4),%xmm1      ## x1 y1  -  - 
        movlps (%rdi,%rbx,4),%xmm4      ## x2 y2  -  -

        movss 8(%rdi,%rax,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm7      ## z2 - - - 

    unpcklps %xmm4,%xmm1 ## x1 x2 y1 y2
    movhlps  %xmm1,%xmm2 ## y1 y2 -  -
    unpcklps %xmm7,%xmm5 ## z1 z2 -  - 

        ## calc dr  
        subps nb030_ix(%rsp),%xmm1
        subps nb030_iy(%rsp),%xmm2
        subps nb030_iz(%rsp),%xmm5

        ## store dr
    movaps %xmm1,nb030_dx(%rsp)
    movaps %xmm2,nb030_dy(%rsp)
    movaps %xmm5,nb030_dz(%rsp)

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
        movaps nb030_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm4
        mulps nb030_half(%rsp),%xmm4
        movaps %xmm4,%xmm2
        mulps  %xmm4,%xmm1
    ## xmm2=rinv
    ## xmm1=r

    mulps nb030_tsc(%rsp),%xmm1   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm1,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm4

    ## multiply by 8
    pslld   $3,%xmm5

    ## calculate eps
    subps     %xmm4,%xmm1

    ## move to integer registers
    movd    %xmm5,%r8d
    pshufd $1,%xmm5,%xmm5
    movd    %xmm5,%r9d

    ## xmm1=eps
    ## xmm2=rinv

        movq nb030_VFtab(%rbp),%rsi
    ## calculate LJ table
    movlps (%rsi,%r8,4),%xmm4
        movlps (%rsi,%r9,4),%xmm5

    unpcklps %xmm5,%xmm4
    movhlps  %xmm4,%xmm5

    movlps 8(%rsi,%r8,4),%xmm6
        movlps 8(%rsi,%r9,4),%xmm7

    unpcklps %xmm7,%xmm6
    movhlps  %xmm6,%xmm7

    movlps 16(%rsi,%r8,4),%xmm8
        movlps 16(%rsi,%r9,4),%xmm9

    unpcklps %xmm9,%xmm8
    movhlps  %xmm8,%xmm9

    movlps 24(%rsi,%r8,4),%xmm10
        movlps 24(%rsi,%r9,4),%xmm11

    unpcklps %xmm11,%xmm10
    movhlps  %xmm10,%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulps  %xmm1,%xmm7   ## Heps
    mulps  %xmm1,%xmm11
    mulps  %xmm1,%xmm6  ## Geps
    mulps  %xmm1,%xmm10
    mulps  %xmm1,%xmm7  ## Heps2
    mulps  %xmm1,%xmm11
    addps  %xmm6,%xmm5 ## F+Geps
    addps  %xmm10,%xmm9
    addps  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addps  %xmm11,%xmm9
    addps  %xmm7,%xmm7   ## 2*Heps2
    addps  %xmm11,%xmm11
    addps  %xmm6,%xmm7  ## 2*Heps2+Geps
    addps  %xmm10,%xmm11

    addps  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addps  %xmm9,%xmm11
    mulps  %xmm1,%xmm5 ## eps*Fp
    mulps  %xmm1,%xmm9
    movaps nb030_c6(%rsp),%xmm12
    movaps nb030_c12(%rsp),%xmm13
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

    mulps  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulps  %xmm13,%xmm9 ## VV*c12 = vnb12
    addps  %xmm9,%xmm5

    xorps  %xmm8,%xmm8
    movlhps %xmm8,%xmm5

    addps  nb030_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb030_Vvdwtot(%rsp)

    mulps  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulps  %xmm13,%xmm11  ## FF*c12  = fnb12
    addps  %xmm11,%xmm7

    mulps  nb030_tsc(%rsp),%xmm7
    mulps  %xmm2,%xmm7
    xorps  %xmm9,%xmm9

    subps  %xmm7,%xmm9
    movaps %xmm9,%xmm10
    movaps %xmm9,%xmm11

    movaps nb030_fix(%rsp),%xmm12
    movaps nb030_fiy(%rsp),%xmm13
    movaps nb030_fiz(%rsp),%xmm14

    mulps  nb030_dx(%rsp),%xmm9
    mulps  nb030_dy(%rsp),%xmm10
    mulps  nb030_dz(%rsp),%xmm11

    movlhps %xmm8,%xmm9
    movlhps %xmm8,%xmm10
    movlhps %xmm8,%xmm11

    ## accumulate i forces
    addps %xmm9,%xmm12
    addps %xmm10,%xmm13
    addps %xmm11,%xmm14
    movaps %xmm12,nb030_fix(%rsp)
    movaps %xmm13,nb030_fiy(%rsp)
    movaps %xmm14,nb030_fiz(%rsp)

        ## the fj's - start by accumulating x & y forces from memory 
        movq nb030_faction(%rbp),%rsi
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

_nb_kernel030_x86_64_sse.nb030_checksingle:     
        movl  nb030_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030_x86_64_sse.nb030_dosingle
        jmp    _nb_kernel030_x86_64_sse.nb030_updateouterdata
_nb_kernel030_x86_64_sse.nb030_dosingle: 
        movq nb030_pos(%rbp),%rdi
        movq  nb030_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb030_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb030_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb030_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
        movss 4(%rsi,%r12,4),%xmm3

    movaps %xmm0,nb030_c6(%rsp)
    movaps %xmm3,nb030_c12(%rsp)

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        movq nb030_pos(%rbp),%rdi
        ## load coordinates
        movss (%rdi,%rax,4),%xmm1
        movss 4(%rdi,%rax,4),%xmm2
        movss 8(%rdi,%rax,4),%xmm5

        ## calc dr  
        subss nb030_ix(%rsp),%xmm1
        subss nb030_iy(%rsp),%xmm2
        subss nb030_iz(%rsp),%xmm5

        ## store dr
    movaps %xmm1,nb030_dx(%rsp)
    movaps %xmm2,nb030_dy(%rsp)
    movaps %xmm5,nb030_dz(%rsp)

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
        movaps nb030_three(%rsp),%xmm4
        mulss %xmm1,%xmm5       ## rsq*lu*lu    
    subss %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm4
        mulss nb030_half(%rsp),%xmm4
        movaps %xmm4,%xmm2
        mulss  %xmm4,%xmm1
    ## xmm2=rinv
    ## xmm1=r

    mulss nb030_tsc(%rsp),%xmm1   ## rtab

    ## truncate and convert to integers
    cvttss2si %xmm1,%r8d

    ## convert back to float
    cvtsi2ss  %r8d,%xmm4

    ## multiply by 8
    shll  $3,%r8d

    ## calculate eps
    subss     %xmm4,%xmm1

    ## xmm1=eps
    ## xmm2=rinv

        movq nb030_VFtab(%rbp),%rsi
    ## calculate LJ table
    movss (%rsi,%r8,4),%xmm4
        movss 4(%rsi,%r8,4),%xmm5
    movss 8(%rsi,%r8,4),%xmm6
        movss 12(%rsi,%r8,4),%xmm7
    movss 16(%rsi,%r8,4),%xmm8
        movss 20(%rsi,%r8,4),%xmm9
    movss 24(%rsi,%r8,4),%xmm10
        movss 28(%rsi,%r8,4),%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulss  %xmm1,%xmm7   ## Heps
    mulss  %xmm1,%xmm11
    mulss  %xmm1,%xmm6  ## Geps
    mulss  %xmm1,%xmm10
    mulss  %xmm1,%xmm7  ## Heps2
    mulss  %xmm1,%xmm11
    addss  %xmm6,%xmm5 ## F+Geps
    addss  %xmm10,%xmm9
    addss  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addss  %xmm11,%xmm9
    addss  %xmm7,%xmm7   ## 2*Heps2
    addss  %xmm11,%xmm11
    addss  %xmm6,%xmm7  ## 2*Heps2+Geps
    addss  %xmm10,%xmm11

    addss  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addss  %xmm9,%xmm11
    mulss  %xmm1,%xmm5 ## eps*Fp
    mulss  %xmm1,%xmm9
    movaps nb030_c6(%rsp),%xmm12
    movaps nb030_c12(%rsp),%xmm13
    addss  %xmm4,%xmm5 ## VV
    addss  %xmm8,%xmm9

    mulss  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulss  %xmm13,%xmm9 ## VV*c12 = vnb12
    addss  %xmm9,%xmm5

    addss  nb030_Vvdwtot(%rsp),%xmm5
    movss %xmm5,nb030_Vvdwtot(%rsp)

    mulss  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulss  %xmm13,%xmm11  ## FF*c12  = fnb12
    addss  %xmm11,%xmm7

    mulss  nb030_tsc(%rsp),%xmm7
    mulss  %xmm2,%xmm7
    xorps  %xmm9,%xmm9

    subss  %xmm7,%xmm9
    movaps %xmm9,%xmm10
    movaps %xmm9,%xmm11

    movaps nb030_fix(%rsp),%xmm12
    movaps nb030_fiy(%rsp),%xmm13
    movaps nb030_fiz(%rsp),%xmm14

    mulss  nb030_dx(%rsp),%xmm9
    mulss  nb030_dy(%rsp),%xmm10
    mulss  nb030_dz(%rsp),%xmm11

    ## accumulate i forces
    addss %xmm9,%xmm12
    addss %xmm10,%xmm13
    addss %xmm11,%xmm14
    movss %xmm12,nb030_fix(%rsp)
    movss %xmm13,nb030_fiy(%rsp)
    movss %xmm14,nb030_fiz(%rsp)

        movq nb030_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel030_x86_64_sse.nb030_updateouterdata: 
        movl  nb030_ii3(%rsp),%ecx
        movq  nb030_faction(%rbp),%rdi
        movq  nb030_fshift(%rbp),%rsi
        movl  nb030_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb030_fix(%rsp),%xmm0
        movaps nb030_fiy(%rsp),%xmm1
        movaps nb030_fiz(%rsp),%xmm2

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
        movl nb030_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb030_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movaps nb030_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb030_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb030_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030_x86_64_sse.nb030_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030_n(%rsp)
        jmp _nb_kernel030_x86_64_sse.nb030_outer
_nb_kernel030_x86_64_sse.nb030_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030_x86_64_sse.nb030_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030_x86_64_sse.nb030_threadloop
_nb_kernel030_x86_64_sse.nb030_end: 

        movl nb030_nouter(%rsp),%eax
        movl nb030_ninner(%rsp),%ebx
        movq nb030_outeriter(%rbp),%rcx
        movq nb030_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $376,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel030nf_x86_64_sse
.globl _nb_kernel030nf_x86_64_sse
nb_kernel030nf_x86_64_sse:      
_nb_kernel030nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb030nf_fshift, 16
.set nb030nf_gid, 24
.set nb030nf_pos, 32
.set nb030nf_faction, 40
.set nb030nf_charge, 48
.set nb030nf_p_facel, 56
.set nb030nf_argkrf, 64
.set nb030nf_argcrf, 72
.set nb030nf_Vc, 80
.set nb030nf_type, 88
.set nb030nf_p_ntype, 96
.set nb030nf_vdwparam, 104
.set nb030nf_Vvdw, 112
.set nb030nf_p_tabscale, 120
.set nb030nf_VFtab, 128
.set nb030nf_invsqrta, 136
.set nb030nf_dvda, 144
.set nb030nf_p_gbtabscale, 152
.set nb030nf_GBtab, 160
.set nb030nf_p_nthreads, 168
.set nb030nf_count, 176
.set nb030nf_mtx, 184
.set nb030nf_outeriter, 192
.set nb030nf_inneriter, 200
.set nb030nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb030nf_ix, 0
.set nb030nf_iy, 16
.set nb030nf_iz, 32
.set nb030nf_tsc, 48
.set nb030nf_c6, 64
.set nb030nf_c12, 80
.set nb030nf_Vvdwtot, 96
.set nb030nf_half, 112
.set nb030nf_three, 128
.set nb030nf_is3, 144
.set nb030nf_ii3, 148
.set nb030nf_innerjjnr, 152
.set nb030nf_nri, 160
.set nb030nf_iinr, 168
.set nb030nf_jindex, 176
.set nb030nf_jjnr, 184
.set nb030nf_shift, 192
.set nb030nf_shiftvec, 200
.set nb030nf_ntia, 208
.set nb030nf_innerk, 212
.set nb030nf_n, 216
.set nb030nf_nn1, 220
.set nb030nf_ntype, 224
.set nb030nf_nouter, 228
.set nb030nf_ninner, 232

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $248,%rsp          # # local variable stack space (n*16+8)
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb030nf_nouter(%rsp)
        movl %eax,nb030nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb030nf_nri(%rsp)
        movq %rsi,nb030nf_iinr(%rsp)
        movq %rdx,nb030nf_jindex(%rsp)
        movq %rcx,nb030nf_jjnr(%rsp)
        movq %r8,nb030nf_shift(%rsp)
        movq %r9,nb030nf_shiftvec(%rsp)
        movq nb030nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb030nf_ntype(%rsp)

        movq nb030nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb030nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb030nf_half(%rsp)
        movss nb030nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb030nf_half(%rsp)
        movaps %xmm3,nb030nf_three(%rsp)

_nb_kernel030nf_x86_64_sse.nb030nf_threadloop: 
        movq  nb030nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel030nf_x86_64_sse.nb030nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030nf_x86_64_sse.nb030nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030nf_n(%rsp)
        movl %ebx,nb030nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030nf_x86_64_sse.nb030nf_outerstart
        jmp _nb_kernel030nf_x86_64_sse.nb030nf_end

_nb_kernel030nf_x86_64_sse.nb030nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030nf_nouter(%rsp),%ebx
        movl %ebx,nb030nf_nouter(%rsp)

_nb_kernel030nf_x86_64_sse.nb030nf_outer: 
        movq  nb030nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb030nf_is3(%rsp)            ## store is3 

        movq  nb030nf_shiftvec(%rsp),%rax       ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb030nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movq  nb030nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb030nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb030nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx                ## rbx = 3*ii=ii3 
        movq  nb030nf_pos(%rbp),%rax            ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb030nf_ix(%rsp)
        movaps %xmm1,nb030nf_iy(%rsp)
        movaps %xmm2,nb030nf_iz(%rsp)

        movl  %ebx,nb030nf_ii3(%rsp)

        ## clear tot potential and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb030nf_Vvdwtot(%rsp)

        movq  nb030nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movq  nb030nf_pos(%rbp),%rsi
        movq  nb030nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb030nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb030nf_ninner(%rsp),%ecx
        movl  %ecx,nb030nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb030nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel030nf_x86_64_sse.nb030nf_unroll_loop
        jmp   _nb_kernel030nf_x86_64_sse.nb030nf_finish_inner
_nb_kernel030nf_x86_64_sse.nb030nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb030nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb030nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movq nb030nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb030nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb030nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb030nf_c6(%rsp)
        movaps %xmm6,nb030nf_c12(%rsp)

        movq nb030nf_pos(%rbp),%rsi        ## base of pos[] 

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

        ## move nb030nf_ix-iz to xmm4-xmm6 
        movaps nb030nf_ix(%rsp),%xmm4
        movaps nb030nf_iy(%rsp),%xmm5
        movaps nb030nf_iz(%rsp),%xmm6

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
        movaps nb030nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030nf_tsc(%rsp),%xmm4

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

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movq nb030nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps 0(%rsi,%rax,4),%xmm5
        movlps 0(%rsi,%rcx,4),%xmm7
        movhps 0(%rsi,%rbx,4),%xmm5
        movhps 0(%rsi,%rdx,4),%xmm7    ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rbx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm3    ## other half of dispersion table 
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

        movaps nb030nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot 
        addps  nb030nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%rsp)

        ## repulsion 
        movlps 16(%rsi,%rax,4),%xmm5
        movlps 16(%rsi,%rcx,4),%xmm7
        movhps 16(%rsi,%rbx,4),%xmm5
        movhps 16(%rsi,%rdx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 24(%rsi,%rax,4),%xmm7
        movlps 24(%rsi,%rcx,4),%xmm3
        movhps 24(%rsi,%rbx,4),%xmm7
        movhps 24(%rsi,%rdx,4),%xmm3    ## other half of repulsion table 
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

        movaps nb030nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5
        addps  nb030nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%rsp)
        xorps  %xmm4,%xmm4

        ## should we do one more iteration? 
        subl $4,nb030nf_innerk(%rsp)
        jl    _nb_kernel030nf_x86_64_sse.nb030nf_finish_inner
        jmp   _nb_kernel030nf_x86_64_sse.nb030nf_unroll_loop
_nb_kernel030nf_x86_64_sse.nb030nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb030nf_innerk(%rsp)
        movl  nb030nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel030nf_x86_64_sse.nb030nf_dopair
        jmp   _nb_kernel030nf_x86_64_sse.nb030nf_checksingle
_nb_kernel030nf_x86_64_sse.nb030nf_dopair: 
        movq  nb030nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb030nf_innerjjnr(%rsp)
        xorps %xmm7,%xmm7

        movq nb030nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb030nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb030nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb030nf_pos(%rbp),%rdi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb030nf_c6(%rsp)
        movaps %xmm6,nb030nf_c12(%rsp)

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

        ## move nb030nf_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb030nf_ix(%rsp),%xmm4
        movaps nb030nf_iy(%rsp),%xmm5
        movaps nb030nf_iz(%rsp),%xmm6

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
        movaps nb030nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6

        movq nb030nf_VFtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## dispersion 
        movlps 0(%rsi,%rcx,4),%xmm5
        movhps 0(%rsi,%rdx,4),%xmm5   ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## 10001000
        shufps $221,%xmm5,%xmm5 ## 11011101

        movlps 8(%rsi,%rcx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm7    ## other half of dispersion table 
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

        movaps nb030nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ##  Update Vvdwtot  
        addps  nb030nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%rsp)

        ## repulsion 
        movlps 16(%rsi,%rcx,4),%xmm5
        movhps 16(%rsi,%rdx,4),%xmm5    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 24(%rsi,%rcx,4),%xmm7
        movhps 24(%rsi,%rdx,4),%xmm7    ## other half of repulsion table 
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

        movaps nb030nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb030nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%rsp)

_nb_kernel030nf_x86_64_sse.nb030nf_checksingle: 
        movl  nb030nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030nf_x86_64_sse.nb030nf_dosingle
        jmp    _nb_kernel030nf_x86_64_sse.nb030nf_updateouterdata
_nb_kernel030nf_x86_64_sse.nb030nf_dosingle: 
        movq nb030nf_pos(%rbp),%rdi
        movq  nb030nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorps  %xmm6,%xmm6

        movq nb030nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb030nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb030nf_ntia(%rsp),%ecx
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100     
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movaps %xmm4,nb030nf_c6(%rsp)
        movaps %xmm6,nb030nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        movaps nb030nf_ix(%rsp),%xmm4
        movaps nb030nf_iy(%rsp),%xmm5
        movaps nb030nf_iz(%rsp),%xmm6

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
        movaps nb030nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6

        movq nb030nf_VFtab(%rbp),%rsi
        movd %mm6,%ebx

        ## dispersion 
        movlps 0(%rsi,%rbx,4),%xmm4
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
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb030nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb030nf_Vvdwtot(%rsp)

        ## repulsion 
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

        movaps nb030nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addss  nb030nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb030nf_Vvdwtot(%rsp)

_nb_kernel030nf_x86_64_sse.nb030nf_updateouterdata: 
        ## get n from stack
        movl nb030nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb030nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movaps nb030nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb030nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb030nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030nf_x86_64_sse.nb030nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030nf_n(%rsp)
        jmp _nb_kernel030nf_x86_64_sse.nb030nf_outer
_nb_kernel030nf_x86_64_sse.nb030nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030nf_x86_64_sse.nb030nf_end ## non-zero, do one more workunit
        jmp   _nb_kernel030nf_x86_64_sse.nb030nf_threadloop
_nb_kernel030nf_x86_64_sse.nb030nf_end: 

        movl nb030nf_nouter(%rsp),%eax
        movl nb030nf_ninner(%rsp),%ebx
        movq nb030nf_outeriter(%rbp),%rcx
        movq nb030nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $248,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret

