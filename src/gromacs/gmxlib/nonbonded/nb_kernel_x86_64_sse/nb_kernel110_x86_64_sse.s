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







.globl nb_kernel110_x86_64_sse
.globl _nb_kernel110_x86_64_sse
nb_kernel110_x86_64_sse:        
_nb_kernel110_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb110_fshift, 16
.set nb110_gid, 24
.set nb110_pos, 32
.set nb110_faction, 40
.set nb110_charge, 48
.set nb110_p_facel, 56
.set nb110_argkrf, 64
.set nb110_argcrf, 72
.set nb110_Vc, 80
.set nb110_type, 88
.set nb110_p_ntype, 96
.set nb110_vdwparam, 104
.set nb110_Vvdw, 112
.set nb110_p_tabscale, 120
.set nb110_VFtab, 128
.set nb110_invsqrta, 136
.set nb110_dvda, 144
.set nb110_p_gbtabscale, 152
.set nb110_GBtab, 160
.set nb110_p_nthreads, 168
.set nb110_count, 176
.set nb110_mtx, 184
.set nb110_outeriter, 192
.set nb110_inneriter, 200
.set nb110_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb110_ix, 0
.set nb110_iy, 16
.set nb110_iz, 32
.set nb110_iq, 48
.set nb110_qq, 64
.set nb110_c6, 112
.set nb110_c12, 128
.set nb110_six, 144
.set nb110_twelve, 160
.set nb110_vctot, 176
.set nb110_Vvdwtot, 192
.set nb110_fix, 208
.set nb110_fiy, 224
.set nb110_fiz, 240
.set nb110_half, 256
.set nb110_three, 272
.set nb110_nri, 288
.set nb110_iinr, 296
.set nb110_jindex, 304
.set nb110_jjnr, 312
.set nb110_shift, 320
.set nb110_shiftvec, 328
.set nb110_facel, 336
.set nb110_innerjjnr, 344
.set nb110_is3, 352
.set nb110_ii3, 356
.set nb110_ntia, 360
.set nb110_innerk, 364
.set nb110_n, 368
.set nb110_nn1, 372
.set nb110_ntype, 376
.set nb110_nouter, 380
.set nb110_ninner, 384



        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        emms
        subq $408,%rsp

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb110_nouter(%rsp)
        movl %eax,nb110_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb110_nri(%rsp)
        movq %rsi,nb110_iinr(%rsp)
        movq %rdx,nb110_jindex(%rsp)
        movq %rcx,nb110_jjnr(%rsp)
        movq %r8,nb110_shift(%rsp)
        movq %r9,nb110_shiftvec(%rsp)
        movq nb110_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb110_ntype(%rsp)
        movq nb110_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb110_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb110_half(%rsp)
        movss nb110_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## six
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## twelve
        movaps %xmm1,nb110_half(%rsp)
        movaps %xmm3,nb110_three(%rsp)
        movaps %xmm4,nb110_six(%rsp)
        movaps %xmm5,nb110_twelve(%rsp)

_nb_kernel110_x86_64_sse.nb110_threadloop: 
        movq  nb110_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel110_x86_64_sse.nb110_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110_x86_64_sse.nb110_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110_n(%rsp)
        movl %ebx,nb110_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel110_x86_64_sse.nb110_outerstart
        jmp _nb_kernel110_x86_64_sse.nb110_end

_nb_kernel110_x86_64_sse.nb110_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110_nouter(%rsp),%ebx
        movl %ebx,nb110_nouter(%rsp)

_nb_kernel110_x86_64_sse.nb110_outer: 
        movq  nb110_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb110_is3(%rsp)      ## store is3 

        movq  nb110_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb110_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb110_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb110_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb110_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb110_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb110_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb110_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb110_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb110_ix(%rsp)
        movaps %xmm1,nb110_iy(%rsp)
        movaps %xmm2,nb110_iz(%rsp)

        movl  %ebx,nb110_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm12,%xmm12
        movaps %xmm12,nb110_vctot(%rsp)
        movaps %xmm12,%xmm13
        movaps %xmm12,%xmm14
        movaps %xmm12,%xmm15

        movq  nb110_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb110_jjnr(%rsp),%rax

        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb110_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb110_ninner(%rsp),%ecx
        movl  %ecx,nb110_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb110_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel110_x86_64_sse.nb110_unroll_loop
        jmp   _nb_kernel110_x86_64_sse.nb110_finish_inner
_nb_kernel110_x86_64_sse.nb110_unroll_loop: 
        ## quad-unrolled innerloop here 
        movq  nb110_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d           ## eax-edx=jnr1-4 

        addq $16,nb110_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx
        lea  (%r10,%r10,2),%rcx
        lea  (%r11,%r11,2),%rdx

        movq nb110_pos(%rbp),%rdi
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

        movq nb110_charge(%rbp),%rsi
        movq nb110_type(%rbp),%rdi

    movaps %xmm1,%xmm4
    unpcklps %xmm2,%xmm1 ## jxa jxc jya jyc        
    unpckhps %xmm2,%xmm4 ## jxb jxd jyb jyd

        movss (%rsi,%r8,4),%xmm7
        movss (%rsi,%r9,4),%xmm8
        movss (%rsi,%r10,4),%xmm10
        movss (%rsi,%r11,4),%xmm11

    movaps %xmm1,%xmm2
    unpcklps %xmm4,%xmm1 ## x
    unpckhps %xmm4,%xmm2 ## y
    shufps  $136,%xmm6,%xmm5  ## 10001000 => jzH2a jzH2b jzH2c jzH2d

        movl (%rdi,%r8,4),%r8d
        movl (%rdi,%r9,4),%r9d
        movl (%rdi,%r10,4),%r10d
        movl (%rdi,%r11,4),%r11d

    unpcklps %xmm10,%xmm7
    unpcklps %xmm11,%xmm8

        ## calc dr  
        subps nb110_ix(%rsp),%xmm1
        subps nb110_iy(%rsp),%xmm2
        subps nb110_iz(%rsp),%xmm5


    unpcklps %xmm8,%xmm7 ## jq

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm2,%xmm10
    movaps %xmm5,%xmm11

    mulps nb110_iq(%rsp),%xmm7

        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
    movl nb110_ntia(%rsp),%edi

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1
        movq nb110_vdwparam(%rbp),%rsi

        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm6
        mulps %xmm5,%xmm5
        movaps nb110_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm6,%xmm4
        mulps nb110_half(%rsp),%xmm4
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

        movlps (%rsi,%r8,4),%xmm3
        movlps (%rsi,%r10,4),%xmm2
        movhps (%rsi,%r9,4),%xmm3
        movhps (%rsi,%r11,4),%xmm2

        movaps %xmm4,%xmm5
        mulps  %xmm4,%xmm5  ## rinv4
        mulps  %xmm4,%xmm5      ## rinv6
        movaps %xmm5,%xmm6
        mulps  %xmm5,%xmm5      ## xmm5=rinv12

    ## coulomb stuff
    mulps  %xmm7,%xmm1 ## vcoul=rinv*qq
    movaps %xmm1,%xmm8 ## fijC

        movaps %xmm3,%xmm0
        shufps $136,%xmm2,%xmm0 ## 10001000
        shufps $221,%xmm2,%xmm3 ## 11011101

    ## add to vctot
    addps  nb110_vctot(%rsp),%xmm1
    movaps %xmm1,nb110_vctot(%rsp)

        mulps  %xmm0,%xmm6  ## vvdw6=c6*rinv6
        mulps  %xmm3,%xmm5  ## vvdw12=c12*rinv12     
        movaps %xmm5,%xmm7
        subps  %xmm6,%xmm5      ## Vvdw=Vvdw12-Vvdw6

        mulps  nb110_six(%rsp),%xmm6
        mulps  nb110_twelve(%rsp),%xmm7
        subps  %xmm6,%xmm7
    addps  %xmm7,%xmm8
        mulps  %xmm8,%xmm4      ## xmm4=total fscal 

        movq nb110_faction(%rbp),%rsi
        ## the fj's - start by combining x & y forces from memory 
        movlps (%rsi,%rax,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%rcx,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%rbx,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%rdx,4),%xmm1 ## x3 y3 x4 y4

    ## add potential to Vvdwtot (sum in xmm12)
        addps  %xmm5,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

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
        subl $4,nb110_innerk(%rsp)
        jl    _nb_kernel110_x86_64_sse.nb110_finish_inner
        jmp   _nb_kernel110_x86_64_sse.nb110_unroll_loop
_nb_kernel110_x86_64_sse.nb110_finish_inner: 
    ## check if at least two particles remain 
    addl $4,nb110_innerk(%rsp)
    movl  nb110_innerk(%rsp),%edx
    andl  $2,%edx
    jnz   _nb_kernel110_x86_64_sse.nb110_dopair
    jmp   _nb_kernel110_x86_64_sse.nb110_checksingle
_nb_kernel110_x86_64_sse.nb110_dopair: 
        ## twice-unrolled innerloop here 
        movq  nb110_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d

        addq $8,nb110_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb110_charge(%rbp),%rsi
        movss (%rsi,%r8,4),%xmm0
        movss (%rsi,%r9,4),%xmm2

    unpcklps %xmm2,%xmm0 ## jqa jqb - -
        mulps nb110_iq(%rsp),%xmm0
    movaps %xmm0,nb110_qq(%rsp)

        movq nb110_type(%rbp),%rsi
        movl (%rsi,%r8,4),%r12d
        movl (%rsi,%r9,4),%r13d
        shll %r12d
        shll %r13d
    movl nb110_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb110_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps  %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    ## xmm0=c6
    ## xmm3=c12

        lea  (%r8,%r8,2),%rax     ## j3 
        lea  (%r9,%r9,2),%rbx

        ## load coordinates
        movq nb110_pos(%rbp),%rdi
        movlps (%rdi,%rax,4),%xmm1      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm2      ## x2 y2 - - 

        movss 8(%rdi,%rax,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm6      ## z2 - - - 

    unpcklps %xmm2,%xmm1 ## x1 x2 y1 y2
    movhlps  %xmm1,%xmm2 ## y1 y2 -  -
    unpcklps %xmm6,%xmm5 ## z1 z2 -  -

        ## calc dr  
        subps nb110_ix(%rsp),%xmm1
        subps nb110_iy(%rsp),%xmm2
        subps nb110_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm2,%xmm10
    movaps %xmm5,%xmm11

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm6
        mulps %xmm5,%xmm5
        movaps nb110_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm6,%xmm4
        mulps nb110_half(%rsp),%xmm4
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

        movaps %xmm4,%xmm5
        mulps  %xmm4,%xmm5  ## rinv4
        mulps  %xmm4,%xmm5      ## rinv6
        movaps %xmm5,%xmm6
        mulps  %xmm5,%xmm5      ## xmm5=rinv12

    ## coulomb stuff
    mulps  nb110_qq(%rsp),%xmm1    ## vcoul=rinv*qq
    movaps %xmm1,%xmm8 ## fijC

    ## add to vctot
    addps  nb110_vctot(%rsp),%xmm1
    movlps %xmm1,nb110_vctot(%rsp)

        mulps  %xmm0,%xmm6  ## vvdw6=c6*rinv6
        mulps  %xmm3,%xmm5  ## vvdw12=c12*rinv12     
        movaps %xmm5,%xmm7
        subps  %xmm6,%xmm5      ## Vvdw=Vvdw12-Vvdw6

        mulps  nb110_six(%rsp),%xmm6
        mulps  nb110_twelve(%rsp),%xmm7
        subps  %xmm6,%xmm7
    addps  %xmm7,%xmm8
        mulps  %xmm8,%xmm4      ## xmm4=total fscal 

    xorps  %xmm7,%xmm7
    movlhps %xmm7,%xmm5

    ## add potential to Vvdwtot (sum in xmm12)
        addps  %xmm5,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

    movlhps %xmm7,%xmm9
    movlhps %xmm7,%xmm10
    movlhps %xmm7,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb110_faction(%rbp),%rsi
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

_nb_kernel110_x86_64_sse.nb110_checksingle:     
    movl  nb110_innerk(%rsp),%edx
    andl  $1,%edx
    jnz    _nb_kernel110_x86_64_sse.nb110_dosingle
    jmp    _nb_kernel110_x86_64_sse.nb110_updateouterdata

_nb_kernel110_x86_64_sse.nb110_dosingle: 
    movq nb110_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb110_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0
    mulss nb110_iq(%rsp),%xmm0
    movaps %xmm0,nb110_qq(%rsp)

        movq nb110_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb110_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb110_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
        movss 4(%rsi,%r12,4),%xmm3

    ## xmm0=c6
    ## xmm3=c12

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        movq nb110_pos(%rbp),%rdi
        ## load coordinates
        movss (%rdi,%rax,4),%xmm1           ## x1 - - - 
        movss 4(%rdi,%rax,4),%xmm2       ## y2 - - - 
        movss 8(%rdi,%rax,4),%xmm5       ## 13 - - - 

        ## calc dr  
        subss nb110_ix(%rsp),%xmm1
        subss nb110_iy(%rsp),%xmm2
        subss nb110_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm2,%xmm10
    movaps %xmm5,%xmm11

        ## square it 
        mulss %xmm1,%xmm1
        mulss %xmm2,%xmm2
        mulss %xmm5,%xmm5
        addss %xmm2,%xmm1
        addss %xmm5,%xmm1
        ## rsq in xmm1

    ## calculate rinv=1/sqrt(rsq)
        rsqrtss %xmm1,%xmm5
        movaps %xmm5,%xmm6
        mulss %xmm5,%xmm5
        movaps nb110_three(%rsp),%xmm4
        mulss %xmm1,%xmm5       ## rsq*lu*lu    
    subss %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulss %xmm6,%xmm4
        mulss nb110_half(%rsp),%xmm4
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

        movaps %xmm4,%xmm5
        mulss  %xmm4,%xmm5  ## rinv4
        mulss  %xmm4,%xmm5      ## rinv6
        movaps %xmm5,%xmm6
        mulss  %xmm5,%xmm5      ## xmm5=rinv12

    ## coulomb stuff
    mulss  nb110_qq(%rsp),%xmm1    ## vcoul=rinv*qq
    movaps %xmm1,%xmm8 ## fijC

    ## add to vctot
    addss  nb110_vctot(%rsp),%xmm1
    movss %xmm1,nb110_vctot(%rsp)

        mulss  %xmm0,%xmm6  ## vvdw6=c6*rinv6
        mulss  %xmm3,%xmm5  ## vvdw12=c12*rinv12     
        movaps %xmm5,%xmm7
        subss  %xmm6,%xmm5      ## Vvdw=Vvdw12-Vvdw6

        mulss  nb110_six(%rsp),%xmm6
        mulss  nb110_twelve(%rsp),%xmm7
        subss  %xmm6,%xmm7
    addss  %xmm7,%xmm8
        mulss  %xmm8,%xmm4      ## xmm4=total fscal 

    ## add potential to Vvdwtot (sum in xmm12)
        addss  %xmm5,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulss  %xmm4,%xmm9
        mulss  %xmm4,%xmm10
        mulss  %xmm4,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb110_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel110_x86_64_sse.nb110_updateouterdata: 
        movl  nb110_ii3(%rsp),%ecx
        movq  nb110_faction(%rbp),%rdi
        movq  nb110_fshift(%rbp),%rsi
        movl  nb110_is3(%rsp),%edx

        ## accumulate i forces in xmm13, xmm14, xmm15
        movhlps %xmm13,%xmm0
        movhlps %xmm14,%xmm1
        movhlps %xmm15,%xmm2
        addps  %xmm13,%xmm0
        addps  %xmm14,%xmm1
        addps  %xmm15,%xmm2
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
        movl nb110_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb110_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb110_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb110_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total potential energy and update it 
        ## accumulate 
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb110_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)

        ## finish if last 
        movl nb110_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110_x86_64_sse.nb110_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110_n(%rsp)
        jmp _nb_kernel110_x86_64_sse.nb110_outer
_nb_kernel110_x86_64_sse.nb110_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110_x86_64_sse.nb110_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110_x86_64_sse.nb110_threadloop
_nb_kernel110_x86_64_sse.nb110_end: 


        movl nb110_nouter(%rsp),%eax
        movl nb110_ninner(%rsp),%ebx
        movq nb110_outeriter(%rbp),%rcx
        movq nb110_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)
        addq $408,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


.globl nb_kernel110nf_x86_64_sse
.globl _nb_kernel110nf_x86_64_sse
nb_kernel110nf_x86_64_sse:      
_nb_kernel110nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb110nf_fshift, 16
.set nb110nf_gid, 24
.set nb110nf_pos, 32
.set nb110nf_faction, 40
.set nb110nf_charge, 48
.set nb110nf_p_facel, 56
.set nb110nf_argkrf, 64
.set nb110nf_argcrf, 72
.set nb110nf_Vc, 80
.set nb110nf_type, 88
.set nb110nf_p_ntype, 96
.set nb110nf_vdwparam, 104
.set nb110nf_Vvdw, 112
.set nb110nf_p_tabscale, 120
.set nb110nf_VFtab, 128
.set nb110nf_invsqrta, 136
.set nb110nf_dvda, 144
.set nb110nf_p_gbtabscale, 152
.set nb110nf_GBtab, 160
.set nb110nf_p_nthreads, 168
.set nb110nf_count, 176
.set nb110nf_mtx, 184
.set nb110nf_outeriter, 192
.set nb110nf_inneriter, 200
.set nb110nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb110nf_ix, 0
.set nb110nf_iy, 16
.set nb110nf_iz, 32
.set nb110nf_iq, 48
.set nb110nf_c6, 64
.set nb110nf_c12, 80
.set nb110nf_vctot, 96
.set nb110nf_Vvdwtot, 112
.set nb110nf_half, 128
.set nb110nf_three, 144
.set nb110nf_nri, 160
.set nb110nf_iinr, 168
.set nb110nf_jindex, 176
.set nb110nf_jjnr, 184
.set nb110nf_shift, 192
.set nb110nf_shiftvec, 200
.set nb110nf_facel, 208
.set nb110nf_innerjjnr, 216
.set nb110nf_is3, 224
.set nb110nf_ii3, 228
.set nb110nf_ntia, 232
.set nb110nf_innerk, 236
.set nb110nf_n, 240
.set nb110nf_nn1, 244
.set nb110nf_ntype, 248
.set nb110nf_nouter, 252
.set nb110nf_ninner, 256

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms
        subq $280,%rsp
        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb110nf_nouter(%rsp)
        movl %eax,nb110nf_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb110nf_nri(%rsp)
        movq %rsi,nb110nf_iinr(%rsp)
        movq %rdx,nb110nf_jindex(%rsp)
        movq %rcx,nb110nf_jjnr(%rsp)
        movq %r8,nb110nf_shift(%rsp)
        movq %r9,nb110nf_shiftvec(%rsp)
        movq nb110nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb110nf_ntype(%rsp)
        movq nb110nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb110nf_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb110nf_half(%rsp)
        movss nb110nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb110nf_half(%rsp)
        movaps %xmm3,nb110nf_three(%rsp)

_nb_kernel110nf_x86_64_sse.nb110nf_threadloop: 
        movq  nb110nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel110nf_x86_64_sse.nb110nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110nf_x86_64_sse.nb110nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110nf_n(%rsp)
        movl %ebx,nb110nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel110nf_x86_64_sse.nb110nf_outerstart
        jmp _nb_kernel110nf_x86_64_sse.nb110nf_end

_nb_kernel110nf_x86_64_sse.nb110nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110nf_nouter(%rsp),%ebx
        movl %ebx,nb110nf_nouter(%rsp)

_nb_kernel110nf_x86_64_sse.nb110nf_outer: 
        movq  nb110nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb110nf_is3(%rsp)            ## store is3 

        movq  nb110nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb110nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb110nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb110nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb110nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb110nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb110nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb110nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb110nf_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb110nf_ix(%rsp)
        movaps %xmm1,nb110nf_iy(%rsp)
        movaps %xmm2,nb110nf_iz(%rsp)

        movl  %ebx,nb110nf_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb110nf_vctot(%rsp)
        movaps %xmm4,nb110nf_Vvdwtot(%rsp)

        movq  nb110nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb110nf_pos(%rbp),%rsi
        movq  nb110nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb110nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb110nf_ninner(%rsp),%ecx
        movl  %ecx,nb110nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb110nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel110nf_x86_64_sse.nb110nf_unroll_loop
        jmp   _nb_kernel110nf_x86_64_sse.nb110nf_finish_inner
_nb_kernel110nf_x86_64_sse.nb110nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb110nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb110nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb110nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        movaps nb110nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movq nb110nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb110nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb110nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb110nf_c6(%rsp)
        movaps %xmm6,nb110nf_c12(%rsp)

        movq nb110nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        mulps %xmm2,%xmm3
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
        movaps nb110nf_ix(%rsp),%xmm4
        movaps nb110nf_iy(%rsp),%xmm5
        movaps nb110nf_iz(%rsp),%xmm6

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
        movaps nb110nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110nf_c6(%rsp),%xmm1
        mulps  nb110nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb110nf_Vvdwtot(%rsp),%xmm5
        addps  nb110nf_vctot(%rsp),%xmm3
        movaps %xmm3,nb110nf_vctot(%rsp)
        movaps %xmm5,nb110nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb110nf_innerk(%rsp)
        jl    _nb_kernel110nf_x86_64_sse.nb110nf_finish_inner
        jmp   _nb_kernel110nf_x86_64_sse.nb110nf_unroll_loop
_nb_kernel110nf_x86_64_sse.nb110nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb110nf_innerk(%rsp)
        movl  nb110nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel110nf_x86_64_sse.nb110nf_dopair
        jmp   _nb_kernel110nf_x86_64_sse.nb110nf_checksingle
_nb_kernel110nf_x86_64_sse.nb110nf_dopair: 
        movq nb110nf_charge(%rbp),%rsi

        movq  nb110nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb110nf_innerjjnr(%rsp)

        xorps %xmm3,%xmm3
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## 00001100 
        shufps $88,%xmm3,%xmm3 ## 01011000 ;# xmm3(0,1) has the charges 

        movq nb110nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb110nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb110nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb110nf_pos(%rbp),%rdi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb110nf_c6(%rsp)
        movaps %xmm6,nb110nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        ## move coordinates to xmm0-xmm2 
        movlps (%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2
        movhps (%rdi,%rbx,4),%xmm1
        movss 8(%rdi,%rbx,4),%xmm0

        mulps  nb110nf_iq(%rsp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## 10001000

        shufps $136,%xmm0,%xmm0 ## 10001000
        shufps $221,%xmm1,%xmm1 ## 11011101

        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb110nf_ix(%rsp),%xmm4
        movaps nb110nf_iy(%rsp),%xmm5
        movaps nb110nf_iz(%rsp),%xmm6

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
        movaps nb110nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110nf_c6(%rsp),%xmm1
        mulps  nb110nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb110nf_Vvdwtot(%rsp),%xmm5
        addps  nb110nf_vctot(%rsp),%xmm3
        movaps %xmm3,nb110nf_vctot(%rsp)
        movaps %xmm5,nb110nf_Vvdwtot(%rsp)

_nb_kernel110nf_x86_64_sse.nb110nf_checksingle: 
        movl  nb110nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel110nf_x86_64_sse.nb110nf_dosingle
        jmp    _nb_kernel110nf_x86_64_sse.nb110nf_updateouterdata
_nb_kernel110nf_x86_64_sse.nb110nf_dosingle: 
        movq nb110nf_charge(%rbp),%rsi
        movq nb110nf_pos(%rbp),%rdi
        movq  nb110nf_innerjjnr(%rsp),%rcx
        xorps %xmm3,%xmm3
        movl  (%rcx),%eax
        movss (%rsi,%rax,4),%xmm3       ## xmm3(0) has the charge       

        movq nb110nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb110nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb110nf_ntia(%rsp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100     
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movaps %xmm4,nb110nf_c6(%rsp)
        movaps %xmm6,nb110nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        mulps  nb110nf_iq(%rsp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb110nf_ix(%rsp),%xmm4
        movaps nb110nf_iy(%rsp),%xmm5
        movaps nb110nf_iz(%rsp),%xmm6

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
        movaps nb110nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb110nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb110nf_c6(%rsp),%xmm1
        mulps  nb110nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb110nf_Vvdwtot(%rsp),%xmm5
        addss  nb110nf_vctot(%rsp),%xmm3
        movss %xmm3,nb110nf_vctot(%rsp)
        movss %xmm5,nb110nf_Vvdwtot(%rsp)

_nb_kernel110nf_x86_64_sse.nb110nf_updateouterdata: 
        ## get n from stack
        movl nb110nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb110nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb110nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb110nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb110nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb110nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb110nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110nf_x86_64_sse.nb110nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110nf_n(%rsp)
        jmp _nb_kernel110nf_x86_64_sse.nb110nf_outer
_nb_kernel110nf_x86_64_sse.nb110nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110nf_x86_64_sse.nb110nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110nf_x86_64_sse.nb110nf_threadloop
_nb_kernel110nf_x86_64_sse.nb110nf_end: 


        movl nb110nf_nouter(%rsp),%eax
        movl nb110nf_ninner(%rsp),%ebx
        movq nb110nf_outeriter(%rbp),%rcx
        movq nb110nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $280,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret



