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






.globl nb_kernel210_x86_64_sse
.globl _nb_kernel210_x86_64_sse
nb_kernel210_x86_64_sse:        
_nb_kernel210_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb210_fshift, 16
.set nb210_gid, 24
.set nb210_pos, 32
.set nb210_faction, 40
.set nb210_charge, 48
.set nb210_p_facel, 56
.set nb210_argkrf, 64
.set nb210_argcrf, 72
.set nb210_Vc, 80
.set nb210_type, 88
.set nb210_p_ntype, 96
.set nb210_vdwparam, 104
.set nb210_Vvdw, 112
.set nb210_p_tabscale, 120
.set nb210_VFtab, 128
.set nb210_invsqrta, 136
.set nb210_dvda, 144
.set nb210_p_gbtabscale, 152
.set nb210_GBtab, 160
.set nb210_p_nthreads, 168
.set nb210_count, 176
.set nb210_mtx, 184
.set nb210_outeriter, 192
.set nb210_inneriter, 200
.set nb210_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb210_ix, 0
.set nb210_iy, 16
.set nb210_iz, 32
.set nb210_iq, 48
.set nb210_c6, 64
.set nb210_c12, 128
.set nb210_six, 144
.set nb210_twelve, 160
.set nb210_vctot, 176
.set nb210_Vvdwtot, 192
.set nb210_fix, 208
.set nb210_fiy, 224
.set nb210_fiz, 240
.set nb210_half, 256
.set nb210_three, 272
.set nb210_two, 288
.set nb210_krf, 304
.set nb210_crf, 320
.set nb210_nri, 336
.set nb210_iinr, 344
.set nb210_jindex, 352
.set nb210_jjnr, 360
.set nb210_shift, 368
.set nb210_shiftvec, 376
.set nb210_facel, 384
.set nb210_innerjjnr, 392
.set nb210_is3, 400
.set nb210_ii3, 404
.set nb210_ntia, 408
.set nb210_innerk, 412
.set nb210_n, 416
.set nb210_nn1, 420
.set nb210_ntype, 424
.set nb210_nouter, 428
.set nb210_ninner, 432


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
        movl %eax,nb210_nouter(%rsp)
        movl %eax,nb210_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb210_nri(%rsp)
        movq %rsi,nb210_iinr(%rsp)
        movq %rdx,nb210_jindex(%rsp)
        movq %rcx,nb210_jjnr(%rsp)
        movq %r8,nb210_shift(%rsp)
        movq %r9,nb210_shiftvec(%rsp)
        movq nb210_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb210_ntype(%rsp)
        movq nb210_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb210_facel(%rsp)


        movq nb210_argkrf(%rbp),%rsi
        movq nb210_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb210_krf(%rsp)
        movaps %xmm2,nb210_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb210_half(%rsp)
        movss nb210_half(%rsp),%xmm1
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
        movaps %xmm1,nb210_half(%rsp)
        movaps %xmm2,nb210_two(%rsp)
        movaps %xmm3,nb210_three(%rsp)
        movaps %xmm4,nb210_six(%rsp)
        movaps %xmm5,nb210_twelve(%rsp)

_nb_kernel210_x86_64_sse.nb210_threadloop: 
        movq  nb210_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel210_x86_64_sse.nb210_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210_x86_64_sse.nb210_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210_n(%rsp)
        movl %ebx,nb210_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210_x86_64_sse.nb210_outerstart
        jmp _nb_kernel210_x86_64_sse.nb210_end

_nb_kernel210_x86_64_sse.nb210_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210_nouter(%rsp),%ebx
        movl %ebx,nb210_nouter(%rsp)

_nb_kernel210_x86_64_sse.nb210_outer: 
        movq  nb210_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb210_is3(%rsp)      ## store is3 

        movq  nb210_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb210_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb210_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb210_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb210_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb210_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb210_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb210_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb210_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb210_ix(%rsp)
        movaps %xmm1,nb210_iy(%rsp)
        movaps %xmm2,nb210_iz(%rsp)

        movl  %ebx,nb210_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm12,%xmm12
        movaps %xmm12,%xmm13
        movaps %xmm12,%xmm14
        movaps %xmm12,%xmm15
        movaps %xmm12,nb210_vctot(%rsp)
        movaps %xmm12,nb210_Vvdwtot(%rsp)

        movq  nb210_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb210_pos(%rbp),%rsi
        movq  nb210_faction(%rbp),%rdi
        movq  nb210_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb210_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb210_ninner(%rsp),%ecx
        movl  %ecx,nb210_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb210_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel210_x86_64_sse.nb210_unroll_loop
        jmp   _nb_kernel210_x86_64_sse.nb210_finish_inner
_nb_kernel210_x86_64_sse.nb210_unroll_loop: 
        ## quad-unrolled innerloop here 
        movq  nb210_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d           ## eax-edx=jnr1-4 

        addq $16,nb210_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx
        lea  (%r10,%r10,2),%rcx
        lea  (%r11,%r11,2),%rdx

        movq nb210_pos(%rbp),%rdi
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

        movq nb210_charge(%rbp),%rsi

    movaps %xmm1,%xmm4
    unpcklps %xmm2,%xmm1 ## jxa jxc jya jyc        
    unpckhps %xmm2,%xmm4 ## jxb jxd jyb jyd
    movaps %xmm1,%xmm2
    unpcklps %xmm4,%xmm1 ## x
    unpckhps %xmm4,%xmm2 ## y
    shufps  $136,%xmm6,%xmm5  ## 10001000 => jzH2a jzH2b jzH2c jzH2d

        ## calc dr  
        subps nb210_ix(%rsp),%xmm1
        subps nb210_iy(%rsp),%xmm2
        subps nb210_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm2,%xmm10
    movaps %xmm5,%xmm11

    ## load charges
        movss (%rsi,%r8,4),%xmm12
        movss (%rsi,%r10,4),%xmm3
        movss (%rsi,%r9,4),%xmm4
        movss (%rsi,%r11,4),%xmm6
    unpcklps %xmm3,%xmm12 ## jqa jqc - -
    unpcklps %xmm6,%xmm4 ## jqb jqd - -
    unpcklps %xmm4,%xmm12 ## jqa jqb jqc jqd
        mulps nb210_iq(%rsp),%xmm12    ## qq


        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

        movq nb210_type(%rbp),%rsi
        movl (%rsi,%r8,4),%r12d
        movl (%rsi,%r9,4),%r13d
        movl (%rsi,%r10,4),%r14d
        movl (%rsi,%r11,4),%r15d

    movaps nb210_krf(%rsp),%xmm7
    mulps  %xmm1,%xmm7   ## krsq

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm6
        mulps %xmm5,%xmm5

        shll %r12d
        shll %r13d
        shll %r14d
        shll %r15d

        movaps nb210_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        movq nb210_vdwparam(%rbp),%rsi

    movl nb210_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d
        addl %edi,%r14d
        addl %edi,%r15d

        mulps %xmm6,%xmm4
        mulps nb210_half(%rsp),%xmm4
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

    ## load c6/c12
        movlps (%rsi,%r12,4),%xmm0
        movlps (%rsi,%r14,4),%xmm2
        movhps (%rsi,%r13,4),%xmm0
        movhps (%rsi,%r15,4),%xmm2

    movaps %xmm1,%xmm3

    subps  %xmm7,%xmm1
    subps  %xmm7,%xmm1  ## rinv-2*krsq
    addps  %xmm7,%xmm3  ## rinv+krsq

    subps  nb210_crf(%rsp),%xmm3   ## rinv+krsq-crf
    mulps  %xmm12,%xmm3 ## vcoul=qq*(rinv+krsq-crf)
    mulps  %xmm12,%xmm1 ## fijC

        movaps %xmm4,%xmm5
        mulps  %xmm4,%xmm5  ## rinv4
        mulps  %xmm4,%xmm5      ## rinv6
        movaps %xmm5,%xmm6
        mulps  %xmm5,%xmm5      ## xmm5=rinv12

        movaps %xmm0,%xmm8
        shufps $136,%xmm2,%xmm0 ## 10001000
        shufps $221,%xmm2,%xmm8 ## 11011101     

    ## add to vctot
    addps  nb210_vctot(%rsp),%xmm3
    movaps %xmm3,nb210_vctot(%rsp)

        mulps  %xmm0,%xmm6  ## vvdw6=c6*rinv6
        mulps  %xmm8,%xmm5  ## vvdw12=c12*rinv12     
        movaps %xmm5,%xmm7
        subps  %xmm6,%xmm5      ## Vvdw=Vvdw12-Vvdw6

        mulps  nb210_six(%rsp),%xmm6
        mulps  nb210_twelve(%rsp),%xmm7
        subps  %xmm6,%xmm7
    addps  %xmm7,%xmm1
        mulps  %xmm1,%xmm4      ## xmm4=total fscal 

        movq nb210_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%rax,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%rcx,4),%xmm1 ## x3 y3 - -

    ## add potential to Vvdwtot 
        addps  nb210_Vvdwtot(%rsp),%xmm5
        movhps (%rsi,%rbx,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%rdx,4),%xmm1 ## x3 y3 x4 y4

    movaps %xmm5,nb210_Vvdwtot(%rsp)

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
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
        subl $4,nb210_innerk(%rsp)
        jl    _nb_kernel210_x86_64_sse.nb210_finish_inner
        jmp   _nb_kernel210_x86_64_sse.nb210_unroll_loop
_nb_kernel210_x86_64_sse.nb210_finish_inner: 
    ## check if at least two particles remain 
    addl $4,nb210_innerk(%rsp)
    movl  nb210_innerk(%rsp),%edx
    andl  $2,%edx
    jnz   _nb_kernel210_x86_64_sse.nb210_dopair
    jmp   _nb_kernel210_x86_64_sse.nb210_checksingle
_nb_kernel210_x86_64_sse.nb210_dopair: 
        ## twice-unrolled innerloop here 
        movq  nb210_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb210_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb210_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm12
        movss (%rsi,%rbx,4),%xmm2

    unpcklps %xmm2,%xmm12 ## jqa jqb - -
        mulps nb210_iq(%rsp),%xmm12     ## qq

        movq nb210_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        shll %r12d
        shll %r13d
    movl nb210_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb210_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps  %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb210_c6(%rsp)
    movaps %xmm3,nb210_c12(%rsp)

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load coordinates
        movq nb210_pos(%rbp),%rdi

        movlps (%rdi,%rax,4),%xmm1      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm2      ## x2 y2 - - 

        movss 8(%rdi,%rax,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm6      ## z2 - - - 

    unpcklps %xmm2,%xmm1 ## x1 x2 y1 y2
    movhlps  %xmm1,%xmm2 ## y1 y2 -  -
    unpcklps %xmm6,%xmm5 ## z1 z2 -  -

        ## calc dr  
        subps nb210_ix(%rsp),%xmm1
        subps nb210_iy(%rsp),%xmm2
        subps nb210_iz(%rsp),%xmm5

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

    movaps nb210_krf(%rsp),%xmm7
    mulps  %xmm1,%xmm7   ## krsq

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm6
        mulps %xmm5,%xmm5
        movaps nb210_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm6,%xmm4
        mulps nb210_half(%rsp),%xmm4
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

    movaps %xmm1,%xmm3

    subps  %xmm7,%xmm1
    subps  %xmm7,%xmm1  ## rinv-2*krsq
    addps  %xmm7,%xmm3  ## rinv+krsq

    subps  nb210_crf(%rsp),%xmm3   ## rinv+krsq-crf
    mulps  %xmm12,%xmm3 ## vcoul=qq*(rinv+krsq-crf)
    mulps  %xmm12,%xmm1 ## fijC

        movaps %xmm4,%xmm5
        mulps  %xmm4,%xmm5  ## rinv4
        mulps  %xmm4,%xmm5      ## rinv6
        movaps %xmm5,%xmm6
        mulps  %xmm5,%xmm5      ## xmm5=rinv12


    ## add to vctot
    addps  nb210_vctot(%rsp),%xmm3
    movlps %xmm3,nb210_vctot(%rsp)

        mulps  nb210_c6(%rsp),%xmm6     ## vvdw6=c6*rinv6
        mulps  nb210_c12(%rsp),%xmm5    ## vvdw12=c12*rinv12     
        movaps %xmm5,%xmm7
        subps  %xmm6,%xmm5      ## Vvdw=Vvdw12-Vvdw6

        mulps  nb210_six(%rsp),%xmm6
        mulps  nb210_twelve(%rsp),%xmm7
        subps  %xmm6,%xmm7
    addps  %xmm7,%xmm1
        mulps  %xmm1,%xmm4      ## xmm4=total fscal 

    xorps  %xmm7,%xmm7
    movlhps %xmm7,%xmm5

    ## add potential to Vvdwtot 
        addps  nb210_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb210_Vvdwtot(%rsp)

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

    movlhps %xmm7,%xmm9
    movlhps %xmm7,%xmm10
    movlhps %xmm7,%xmm11

        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb210_faction(%rbp),%rsi
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

_nb_kernel210_x86_64_sse.nb210_checksingle:     
    movl  nb210_innerk(%rsp),%edx
    andl  $1,%edx
    jnz    _nb_kernel210_x86_64_sse.nb210_dosingle
    jmp    _nb_kernel210_x86_64_sse.nb210_updateouterdata

_nb_kernel210_x86_64_sse.nb210_dosingle: 
    movq nb210_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb210_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm12
    mulss nb210_iq(%rsp),%xmm12   ## qq

        movq nb210_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb210_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb210_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
        movss 4(%rsi,%r12,4),%xmm3

    movaps %xmm0,nb210_c6(%rsp)
    movaps %xmm3,nb210_c12(%rsp)

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        movq nb210_pos(%rbp),%rdi
        ## load coordinates
        movss (%rdi,%rax,4),%xmm1           ## x1 - - - 
        movss 4(%rdi,%rax,4),%xmm2       ## y2 - - - 
        movss 8(%rdi,%rax,4),%xmm5       ## 13 - - - 

        ## calc dr  
        subss nb210_ix(%rsp),%xmm1
        subss nb210_iy(%rsp),%xmm2
        subss nb210_iz(%rsp),%xmm5

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

    movaps nb210_krf(%rsp),%xmm7
    mulss  %xmm1,%xmm7   ## krsq

    ## calculate rinv=1/sqrt(rsq)
        rsqrtss %xmm1,%xmm5
        movaps %xmm5,%xmm6
        mulss %xmm5,%xmm5
        movaps nb210_three(%rsp),%xmm4
        mulss %xmm1,%xmm5       ## rsq*lu*lu    
    subss %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulss %xmm6,%xmm4
        mulss nb210_half(%rsp),%xmm4
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

    movaps %xmm1,%xmm3

    subss  %xmm7,%xmm1
    subss  %xmm7,%xmm1  ## rinv-2*krsq
    addss  %xmm7,%xmm3  ## rinv+krsq

    subss  nb210_crf(%rsp),%xmm3   ## rinv+krsq-crf
    mulss  %xmm12,%xmm3 ## vcoul=qq*(rinv+krsq-crf)
    mulss  %xmm12,%xmm1 ## fijC

        movaps %xmm4,%xmm5
        mulss  %xmm4,%xmm5  ## rinv4
        mulss  %xmm4,%xmm5      ## rinv6
        movaps %xmm5,%xmm6
        mulss  %xmm5,%xmm5      ## xmm5=rinv12


    ## add to vctot
    addss  nb210_vctot(%rsp),%xmm3
    movss %xmm3,nb210_vctot(%rsp)

        mulss  nb210_c6(%rsp),%xmm6     ## vvdw6=c6*rinv6
        mulss  nb210_c12(%rsp),%xmm5     ## vvdw12=c12*rinv12     
        movaps %xmm5,%xmm7
        subss  %xmm6,%xmm5      ## Vvdw=Vvdw12-Vvdw6

        mulss  nb210_six(%rsp),%xmm6
        mulss  nb210_twelve(%rsp),%xmm7
        subss  %xmm6,%xmm7
    addss  %xmm7,%xmm1
        mulss  %xmm1,%xmm4      ## xmm4=total fscal 

    ## add potential to Vvdwtot 
        addss  nb210_Vvdwtot(%rsp),%xmm5
    movss  %xmm5,nb210_Vvdwtot(%rsp)

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulss  %xmm4,%xmm9
        mulss  %xmm4,%xmm10
        mulss  %xmm4,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb210_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel210_x86_64_sse.nb210_updateouterdata: 
        movl  nb210_ii3(%rsp),%ecx
        movq  nb210_faction(%rbp),%rdi
        movq  nb210_fshift(%rbp),%rsi
        movl  nb210_is3(%rsp),%edx

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
        movl nb210_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb210_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb210_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb210_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total potential energy and update it 
        movaps nb210_Vvdwtot(%rsp),%xmm12
    ## accumulate
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb210_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)

        ## finish if last 
        movl nb210_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210_x86_64_sse.nb210_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210_n(%rsp)
        jmp _nb_kernel210_x86_64_sse.nb210_outer
_nb_kernel210_x86_64_sse.nb210_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210_x86_64_sse.nb210_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210_x86_64_sse.nb210_threadloop
_nb_kernel210_x86_64_sse.nb210_end: 
        movl nb210_nouter(%rsp),%eax
        movl nb210_ninner(%rsp),%ebx
        movq nb210_outeriter(%rbp),%rcx
        movq nb210_inneriter(%rbp),%rdx
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







.globl nb_kernel210nf_x86_64_sse
.globl _nb_kernel210nf_x86_64_sse
nb_kernel210nf_x86_64_sse:      
_nb_kernel210nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb210nf_fshift, 16
.set nb210nf_gid, 24
.set nb210nf_pos, 32
.set nb210nf_faction, 40
.set nb210nf_charge, 48
.set nb210nf_p_facel, 56
.set nb210nf_argkrf, 64
.set nb210nf_argcrf, 72
.set nb210nf_Vc, 80
.set nb210nf_type, 88
.set nb210nf_p_ntype, 96
.set nb210nf_vdwparam, 104
.set nb210nf_Vvdw, 112
.set nb210nf_p_tabscale, 120
.set nb210nf_VFtab, 128
.set nb210nf_invsqrta, 136
.set nb210nf_dvda, 144
.set nb210nf_p_gbtabscale, 152
.set nb210nf_GBtab, 160
.set nb210nf_p_nthreads, 168
.set nb210nf_count, 176
.set nb210nf_mtx, 184
.set nb210nf_outeriter, 192
.set nb210nf_inneriter, 200
.set nb210nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb210nf_ix, 0
.set nb210nf_iy, 16
.set nb210nf_iz, 32
.set nb210nf_iq, 48
.set nb210nf_c6, 64
.set nb210nf_c12, 80
.set nb210nf_vctot, 96
.set nb210nf_Vvdwtot, 112
.set nb210nf_half, 128
.set nb210nf_three, 144
.set nb210nf_krf, 160
.set nb210nf_crf, 176
.set nb210nf_nri, 192
.set nb210nf_iinr, 200
.set nb210nf_jindex, 208
.set nb210nf_jjnr, 216
.set nb210nf_shift, 224
.set nb210nf_shiftvec, 232
.set nb210nf_facel, 240
.set nb210nf_innerjjnr, 248
.set nb210nf_is3, 256
.set nb210nf_ii3, 260
.set nb210nf_ntia, 264
.set nb210nf_innerk, 268
.set nb210nf_n, 272
.set nb210nf_nn1, 276
.set nb210nf_ntype, 280
.set nb210nf_nouter, 284
.set nb210nf_ninner, 288

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
        movl %eax,nb210nf_nouter(%rsp)
        movl %eax,nb210nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb210nf_nri(%rsp)
        movq %rsi,nb210nf_iinr(%rsp)
        movq %rdx,nb210nf_jindex(%rsp)
        movq %rcx,nb210nf_jjnr(%rsp)
        movq %r8,nb210nf_shift(%rsp)
        movq %r9,nb210nf_shiftvec(%rsp)
        movq nb210nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb210nf_ntype(%rsp)
        movq nb210nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb210nf_facel(%rsp)


        movq nb210nf_argkrf(%rbp),%rsi
        movq nb210nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb210nf_krf(%rsp)
        movaps %xmm2,nb210nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb210nf_half(%rsp)
        movss nb210nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb210nf_half(%rsp)
        movaps %xmm3,nb210nf_three(%rsp)


_nb_kernel210nf_x86_64_sse.nb210nf_threadloop: 
        movq  nb210nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel210nf_x86_64_sse.nb210nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210nf_x86_64_sse.nb210nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210nf_n(%rsp)
        movl %ebx,nb210nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210nf_x86_64_sse.nb210nf_outerstart
        jmp _nb_kernel210nf_x86_64_sse.nb210nf_end

_nb_kernel210nf_x86_64_sse.nb210nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210nf_nouter(%rsp),%ebx
        movl %ebx,nb210nf_nouter(%rsp)

_nb_kernel210nf_x86_64_sse.nb210nf_outer: 
        movq  nb210nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb210nf_is3(%rsp)            ## store is3 

        movq  nb210nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb210nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb210nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb210nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb210nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb210nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb210nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb210nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb210nf_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb210nf_ix(%rsp)
        movaps %xmm1,nb210nf_iy(%rsp)
        movaps %xmm2,nb210nf_iz(%rsp)

        movl  %ebx,nb210nf_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb210nf_vctot(%rsp)
        movaps %xmm4,nb210nf_Vvdwtot(%rsp)

        movq  nb210nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb210nf_pos(%rbp),%rsi
        movq  nb210nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb210nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb210nf_ninner(%rsp),%ecx
        movl  %ecx,nb210nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb210nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel210nf_x86_64_sse.nb210nf_unroll_loop
        jmp   _nb_kernel210nf_x86_64_sse.nb210nf_finish_inner
_nb_kernel210nf_x86_64_sse.nb210nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb210nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb210nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb210nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        movaps nb210nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movq nb210nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb210nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb210nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb210nf_c6(%rsp)
        movaps %xmm6,nb210nf_c12(%rsp)

        movq nb210nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb210nf_ix(%rsp),%xmm4
        movaps nb210nf_iy(%rsp),%xmm5
        movaps nb210nf_iz(%rsp),%xmm6

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

        movaps nb210nf_krf(%rsp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210nf_half(%rsp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210nf_crf(%rsp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulps  nb210nf_c6(%rsp),%xmm1
        mulps  nb210nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb210nf_Vvdwtot(%rsp),%xmm5
        addps  nb210nf_vctot(%rsp),%xmm6
        movaps %xmm6,nb210nf_vctot(%rsp)
        movaps %xmm5,nb210nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb210nf_innerk(%rsp)
        jl    _nb_kernel210nf_x86_64_sse.nb210nf_finish_inner
        jmp   _nb_kernel210nf_x86_64_sse.nb210nf_unroll_loop
_nb_kernel210nf_x86_64_sse.nb210nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb210nf_innerk(%rsp)
        movl  nb210nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel210nf_x86_64_sse.nb210nf_dopair
        jmp   _nb_kernel210nf_x86_64_sse.nb210nf_checksingle
_nb_kernel210nf_x86_64_sse.nb210nf_dopair: 
        movq nb210nf_charge(%rbp),%rsi

        movq  nb210nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb210nf_innerjjnr(%rsp)

        xorps %xmm3,%xmm3
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## 00001100 
        shufps $88,%xmm3,%xmm3 ## 01011000 ;# xmm3(0,1) has the charges 

        movq nb210nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb210nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb210nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb210nf_pos(%rbp),%rdi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb210nf_c6(%rsp)
        movaps %xmm6,nb210nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        ## move coordinates to xmm0-xmm2 
        movlps (%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2
        movhps (%rdi,%rbx,4),%xmm1
        movss 8(%rdi,%rbx,4),%xmm0

        mulps  nb210nf_iq(%rsp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## 10001000

        shufps $136,%xmm0,%xmm0 ## 10001000
        shufps $221,%xmm1,%xmm1 ## 11011101

        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb210nf_ix(%rsp),%xmm4
        movaps nb210nf_iy(%rsp),%xmm5
        movaps nb210nf_iz(%rsp),%xmm6

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

        movaps nb210nf_krf(%rsp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210nf_half(%rsp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210nf_crf(%rsp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq-crf) 
        mulps  nb210nf_c6(%rsp),%xmm1
        mulps  nb210nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb210nf_Vvdwtot(%rsp),%xmm5
        addps  nb210nf_vctot(%rsp),%xmm6
        movaps %xmm6,nb210nf_vctot(%rsp)
        movaps %xmm5,nb210nf_Vvdwtot(%rsp)

_nb_kernel210nf_x86_64_sse.nb210nf_checksingle: 
        movl  nb210nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel210nf_x86_64_sse.nb210nf_dosingle
        jmp    _nb_kernel210nf_x86_64_sse.nb210nf_updateouterdata
_nb_kernel210nf_x86_64_sse.nb210nf_dosingle: 
        movq nb210nf_charge(%rbp),%rsi
        movq nb210nf_pos(%rbp),%rdi
        movq  nb210nf_innerjjnr(%rsp),%rcx
        xorps %xmm3,%xmm3
        movl  (%rcx),%eax
        movss (%rsi,%rax,4),%xmm3       ## xmm3(0) has the charge       

        movq nb210nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb210nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb210nf_ntia(%rsp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100     
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movaps %xmm4,nb210nf_c6(%rsp)
        movaps %xmm6,nb210nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        mulps  nb210nf_iq(%rsp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb210nf_ix(%rsp),%xmm4
        movaps nb210nf_iy(%rsp),%xmm5
        movaps nb210nf_iz(%rsp),%xmm6

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

        movaps nb210nf_krf(%rsp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb210nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb210nf_half(%rsp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movaps %xmm4,%xmm1
        subps  nb210nf_crf(%rsp),%xmm6
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul 
        mulps  nb210nf_c6(%rsp),%xmm1
        mulps  nb210nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb210nf_Vvdwtot(%rsp),%xmm5
        addss  nb210nf_vctot(%rsp),%xmm6
        movss %xmm6,nb210nf_vctot(%rsp)
        movss %xmm5,nb210nf_Vvdwtot(%rsp)

_nb_kernel210nf_x86_64_sse.nb210nf_updateouterdata: 
        ## get n from stack
        movl nb210nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb210nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb210nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb210nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb210nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb210nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb210nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210nf_x86_64_sse.nb210nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210nf_n(%rsp)
        jmp _nb_kernel210nf_x86_64_sse.nb210nf_outer
_nb_kernel210nf_x86_64_sse.nb210nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210nf_x86_64_sse.nb210nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210nf_x86_64_sse.nb210nf_threadloop
_nb_kernel210nf_x86_64_sse.nb210nf_end: 

        movl nb210nf_nouter(%rsp),%eax
        movl nb210nf_ninner(%rsp),%ebx
        movq nb210nf_outeriter(%rbp),%rcx
        movq nb210nf_inneriter(%rbp),%rdx
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

