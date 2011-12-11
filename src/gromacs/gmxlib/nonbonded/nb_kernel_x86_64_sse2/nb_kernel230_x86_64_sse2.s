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






.globl nb_kernel230_x86_64_sse2
.globl _nb_kernel230_x86_64_sse2
nb_kernel230_x86_64_sse2:       
_nb_kernel230_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb230_fshift, 16
.set nb230_gid, 24
.set nb230_pos, 32
.set nb230_faction, 40
.set nb230_charge, 48
.set nb230_p_facel, 56
.set nb230_argkrf, 64
.set nb230_argcrf, 72
.set nb230_Vc, 80
.set nb230_type, 88
.set nb230_p_ntype, 96
.set nb230_vdwparam, 104
.set nb230_Vvdw, 112
.set nb230_p_tabscale, 120
.set nb230_VFtab, 128
.set nb230_invsqrta, 136
.set nb230_dvda, 144
.set nb230_p_gbtabscale, 152
.set nb230_GBtab, 160
.set nb230_p_nthreads, 168
.set nb230_count, 176
.set nb230_mtx, 184
.set nb230_outeriter, 192
.set nb230_inneriter, 200
.set nb230_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb230_ix, 0
.set nb230_iy, 16
.set nb230_iz, 32
.set nb230_iq, 48
.set nb230_dx, 64
.set nb230_dy, 80
.set nb230_dz, 96
.set nb230_c6, 112
.set nb230_c12, 128
.set nb230_tsc, 144
.set nb230_fstmp, 160
.set nb230_vctot, 176
.set nb230_Vvdwtot, 192
.set nb230_fix, 208
.set nb230_fiy, 224
.set nb230_fiz, 240
.set nb230_half, 256
.set nb230_three, 272
.set nb230_two, 288
.set nb230_krf, 304
.set nb230_crf, 320
.set nb230_nri, 336
.set nb230_iinr, 344
.set nb230_jindex, 352
.set nb230_jjnr, 360
.set nb230_shift, 368
.set nb230_shiftvec, 376
.set nb230_facel, 384
.set nb230_innerjjnr, 392
.set nb230_is3, 400
.set nb230_ii3, 404
.set nb230_ntia, 408
.set nb230_innerk, 412
.set nb230_n, 416
.set nb230_nn1, 420
.set nb230_ntype, 424
.set nb230_nouter, 428
.set nb230_ninner, 432
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $456,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb230_nouter(%rsp)
        movl %eax,nb230_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb230_nri(%rsp)
        movq %rsi,nb230_iinr(%rsp)
        movq %rdx,nb230_jindex(%rsp)
        movq %rcx,nb230_jjnr(%rsp)
        movq %r8,nb230_shift(%rsp)
        movq %r9,nb230_shiftvec(%rsp)
        movq nb230_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb230_ntype(%rsp)
        movq nb230_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb230_facel(%rsp)

        movq nb230_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb230_tsc(%rsp)

        movq nb230_argkrf(%rbp),%rsi
        movq nb230_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb230_krf(%rsp)
        movapd %xmm2,nb230_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb230_half(%rsp)
        movl %ebx,nb230_half+4(%rsp)
        movsd nb230_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb230_half(%rsp)
        movapd %xmm2,nb230_two(%rsp)
        movapd %xmm3,nb230_three(%rsp)

_nb_kernel230_x86_64_sse2.nb230_threadloop: 
        movq  nb230_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel230_x86_64_sse2.nb230_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel230_x86_64_sse2.nb230_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb230_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb230_n(%rsp)
        movl %ebx,nb230_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel230_x86_64_sse2.nb230_outerstart
        jmp _nb_kernel230_x86_64_sse2.nb230_end

_nb_kernel230_x86_64_sse2.nb230_outerstart: 
        ## ebx contains number of outer iterations
        addl nb230_nouter(%rsp),%ebx
        movl %ebx,nb230_nouter(%rsp)

_nb_kernel230_x86_64_sse2.nb230_outer: 
        movq  nb230_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb230_is3(%rsp)      ## store is3 

        movq  nb230_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb230_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb230_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb230_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb230_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb230_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb230_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb230_pos(%rbp),%rax      ## eax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb230_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb230_ix(%rsp)
        movapd %xmm1,nb230_iy(%rsp)
        movapd %xmm2,nb230_iz(%rsp)

        movl  %ebx,nb230_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb230_vctot(%rsp)
        movapd %xmm4,nb230_Vvdwtot(%rsp)
        movapd %xmm4,nb230_fix(%rsp)
        movapd %xmm4,nb230_fiy(%rsp)
        movapd %xmm4,nb230_fiz(%rsp)

        movq  nb230_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb230_pos(%rbp),%rsi
        movq  nb230_faction(%rbp),%rdi
        movq  nb230_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb230_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb230_ninner(%rsp),%ecx
        movl  %ecx,nb230_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb230_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel230_x86_64_sse2.nb230_unroll_loop
        jmp   _nb_kernel230_x86_64_sse2.nb230_checksingle
_nb_kernel230_x86_64_sse2.nb230_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb230_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb230_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb230_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        mulpd nb230_iq(%rsp),%xmm3              ## qq 

        movq nb230_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb230_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb230_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb230_c6(%rsp)
        movapd %xmm6,nb230_c12(%rsp)

        movq nb230_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        ## calc dr 
        subpd nb230_ix(%rsp),%xmm4
        subpd nb230_iy(%rsp),%xmm5
        subpd nb230_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb230_dx(%rsp)
        movapd %xmm5,nb230_dy(%rsp)
        movapd %xmm6,nb230_dz(%rsp)
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

        movapd nb230_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb230_three(%rsp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb230_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb230_crf(%rsp),%xmm6
        mulpd  %xmm3,%xmm6
        mulpd  nb230_two(%rsp),%xmm7
        movapd %xmm0,%xmm1
        subpd  %xmm7,%xmm1  ## rinv-2*krsq
        mulpd  %xmm0,%xmm1  ## (rinv-2*krsq)*rinv
        mulpd  %xmm1,%xmm3  ## qq*(rinv-2*krsq)*rinv

    ## xmm3=fstmp

        addpd  nb230_vctot(%rsp),%xmm6
        movapd %xmm6,nb230_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb230_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movq nb230_VFtab(%rbp),%rsi
        movd %mm6,%r8d
        psrlq $32,%mm6
        movd %mm6,%r9d

        ## load both disp. and rep. tables in parallel
        movlpd (%rsi,%r8,8),%xmm4
        movlpd 8(%rsi,%r8,8),%xmm5
        movlpd 16(%rsi,%r8,8),%xmm6
        movlpd 24(%rsi,%r8,8),%xmm7
        movlpd 32(%rsi,%r8,8),%xmm8
        movlpd 40(%rsi,%r8,8),%xmm9
        movlpd 48(%rsi,%r8,8),%xmm10
        movlpd 56(%rsi,%r8,8),%xmm11
        movhpd (%rsi,%r9,8),%xmm4
        movhpd 8(%rsi,%r9,8),%xmm5
        movhpd 16(%rsi,%r9,8),%xmm6
        movhpd 24(%rsi,%r9,8),%xmm7
        movhpd 32(%rsi,%r9,8),%xmm8
        movhpd 40(%rsi,%r9,8),%xmm9
        movhpd 48(%rsi,%r9,8),%xmm10
        movhpd 56(%rsi,%r9,8),%xmm11
        ## dispersion table ready in xmm4-xmm7, repulsion in xmm8-xmm11

    mulpd  %xmm1,%xmm7   ## Heps
    mulpd  %xmm1,%xmm11
    mulpd  %xmm1,%xmm6  ## Geps
    mulpd  %xmm1,%xmm10
    mulpd  %xmm1,%xmm7  ## Heps2
    mulpd  %xmm1,%xmm11
    addpd  %xmm6,%xmm5 ## F+Geps
    addpd  %xmm10,%xmm9
    addpd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addpd  %xmm11,%xmm9
    addpd  %xmm7,%xmm7   ## 2*Heps2
    addpd  %xmm11,%xmm11
    addpd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addpd  %xmm10,%xmm11

    addpd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm9,%xmm11
    mulpd  %xmm1,%xmm5 ## eps*Fp
    mulpd  %xmm1,%xmm9
    movapd nb230_c6(%rsp),%xmm12
    movapd nb230_c12(%rsp),%xmm13
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb230_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb230_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7

    mulpd  nb230_tsc(%rsp),%xmm7
    subpd  %xmm7,%xmm3
    mulpd  %xmm0,%xmm3  ## fscal

    movapd %xmm3,%xmm9
    movapd %xmm3,%xmm10
    movapd %xmm3,%xmm11

    movapd nb230_fix(%rsp),%xmm12
    movapd nb230_fiy(%rsp),%xmm13
    movapd nb230_fiz(%rsp),%xmm14

    mulpd  nb230_dx(%rsp),%xmm9
    mulpd  nb230_dy(%rsp),%xmm10
    mulpd  nb230_dz(%rsp),%xmm11

    ## accumulate i forces
    addpd %xmm9,%xmm12
    addpd %xmm10,%xmm13
    addpd %xmm11,%xmm14
    movapd %xmm12,nb230_fix(%rsp)
    movapd %xmm13,nb230_fiy(%rsp)
    movapd %xmm14,nb230_fiz(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb230_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        movhpd (%rdi,%rbx,8),%xmm3
        movhpd 8(%rdi,%rbx,8),%xmm4
        movhpd 16(%rdi,%rbx,8),%xmm5
        addpd %xmm9,%xmm3
        addpd %xmm10,%xmm4
        addpd %xmm11,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)
        movhpd %xmm3,(%rdi,%rbx,8)
        movhpd %xmm4,8(%rdi,%rbx,8)
        movhpd %xmm5,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb230_innerk(%rsp)
        jl    _nb_kernel230_x86_64_sse2.nb230_checksingle
        jmp   _nb_kernel230_x86_64_sse2.nb230_unroll_loop

_nb_kernel230_x86_64_sse2.nb230_checksingle:    
        movl  nb230_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel230_x86_64_sse2.nb230_dosingle
        jmp    _nb_kernel230_x86_64_sse2.nb230_updateouterdata
_nb_kernel230_x86_64_sse2.nb230_dosingle: 
        movq nb230_charge(%rbp),%rsi
        movq nb230_pos(%rbp),%rdi
        movq  nb230_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb230_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
        mulsd nb230_iq(%rsp),%xmm3              ## qq 

        movq nb230_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb230_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb230_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb230_c6(%rsp)
        movapd %xmm6,nb230_c12(%rsp)

        movq nb230_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb230_ix(%rsp),%xmm4
        subsd nb230_iy(%rsp),%xmm5
        subsd nb230_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb230_dx(%rsp)
        movapd %xmm5,nb230_dy(%rsp)
        movapd %xmm6,nb230_dz(%rsp)
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

        movapd nb230_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb230_three(%rsp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb230_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb230_crf(%rsp),%xmm6
        mulsd  %xmm3,%xmm6
        mulsd  nb230_two(%rsp),%xmm7
        movapd %xmm0,%xmm1
        subsd  %xmm7,%xmm1  ## rinv-2*krsq
        mulsd  %xmm0,%xmm1  ## (rinv-2*krsq)*rinv
        mulsd  %xmm1,%xmm3  ## qq*(rinv-2*krsq)*rinv

    ## xmm3=fstmp

        addsd  nb230_vctot(%rsp),%xmm6
        movsd %xmm6,nb230_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq
        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb230_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%r8d            ## idx *= 8 

        movq nb230_VFtab(%rbp),%rsi


        ## load both disp. and rep. tables in parallel
    movsd (%rsi,%r8,8),%xmm4
    movsd 8(%rsi,%r8,8),%xmm5
    movsd 16(%rsi,%r8,8),%xmm6
    movsd 24(%rsi,%r8,8),%xmm7
    movsd 32(%rsi,%r8,8),%xmm8
    movsd 40(%rsi,%r8,8),%xmm9
    movsd 48(%rsi,%r8,8),%xmm10
    movsd 56(%rsi,%r8,8),%xmm11
        ## dispersion table ready in xmm4-xmm7, repulsion in xmm8-xmm11

    mulsd  %xmm1,%xmm7   ## Heps
    mulsd  %xmm1,%xmm11
    mulsd  %xmm1,%xmm6  ## Geps
    mulsd  %xmm1,%xmm10
    mulsd  %xmm1,%xmm7  ## Heps2
    mulsd  %xmm1,%xmm11
    addsd  %xmm6,%xmm5 ## F+Geps
    addsd  %xmm10,%xmm9
    addsd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addsd  %xmm11,%xmm9
    addsd  %xmm7,%xmm7   ## 2*Heps2
    addsd  %xmm11,%xmm11
    addsd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addsd  %xmm10,%xmm11

    addsd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm9,%xmm11
    mulsd  %xmm1,%xmm5 ## eps*Fp
    mulsd  %xmm1,%xmm9
    movapd nb230_c6(%rsp),%xmm12
    movapd nb230_c12(%rsp),%xmm13
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulsd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb230_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb230_Vvdwtot(%rsp)

    mulsd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulsd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7

    mulsd  nb230_tsc(%rsp),%xmm7
    subsd  %xmm7,%xmm3
    mulsd  %xmm0,%xmm3  ## fscal

    movapd %xmm3,%xmm9
    movapd %xmm3,%xmm10
    movapd %xmm3,%xmm11

    movapd nb230_fix(%rsp),%xmm12
    movapd nb230_fiy(%rsp),%xmm13
    movapd nb230_fiz(%rsp),%xmm14

    mulsd  nb230_dx(%rsp),%xmm9
    mulsd  nb230_dy(%rsp),%xmm10
    mulsd  nb230_dz(%rsp),%xmm11

    ## accumulate i forces
    addsd %xmm9,%xmm12
    addsd %xmm10,%xmm13
    addsd %xmm11,%xmm14
    movsd %xmm12,nb230_fix(%rsp)
    movsd %xmm13,nb230_fiy(%rsp)
    movsd %xmm14,nb230_fiz(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb230_faction(%rbp),%rdi
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel230_x86_64_sse2.nb230_updateouterdata: 
        movl  nb230_ii3(%rsp),%ecx
        movq  nb230_faction(%rbp),%rdi
        movq  nb230_fshift(%rbp),%rsi
        movl  nb230_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb230_fix(%rsp),%xmm0
        movapd nb230_fiy(%rsp),%xmm1
        movapd nb230_fiz(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addpd  %xmm3,%xmm0
        addpd  %xmm4,%xmm1
        addpd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movl nb230_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb230_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb230_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb230_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb230_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb230_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb230_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel230_x86_64_sse2.nb230_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb230_n(%rsp)
        jmp _nb_kernel230_x86_64_sse2.nb230_outer
_nb_kernel230_x86_64_sse2.nb230_outerend: 
        ## check if more outer neighborlists remain
        movl  nb230_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel230_x86_64_sse2.nb230_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel230_x86_64_sse2.nb230_threadloop
_nb_kernel230_x86_64_sse2.nb230_end: 
        movl nb230_nouter(%rsp),%eax
        movl nb230_ninner(%rsp),%ebx
        movq nb230_outeriter(%rbp),%rcx
        movq nb230_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $456,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel230nf_x86_64_sse2
.globl _nb_kernel230nf_x86_64_sse2
nb_kernel230nf_x86_64_sse2:     
_nb_kernel230nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb230nf_fshift, 16
.set nb230nf_gid, 24
.set nb230nf_pos, 32
.set nb230nf_faction, 40
.set nb230nf_charge, 48
.set nb230nf_p_facel, 56
.set nb230nf_argkrf, 64
.set nb230nf_argcrf, 72
.set nb230nf_Vc, 80
.set nb230nf_type, 88
.set nb230nf_p_ntype, 96
.set nb230nf_vdwparam, 104
.set nb230nf_Vvdw, 112
.set nb230nf_p_tabscale, 120
.set nb230nf_VFtab, 128
.set nb230nf_invsqrta, 136
.set nb230nf_dvda, 144
.set nb230nf_p_gbtabscale, 152
.set nb230nf_GBtab, 160
.set nb230nf_p_nthreads, 168
.set nb230nf_count, 176
.set nb230nf_mtx, 184
.set nb230nf_outeriter, 192
.set nb230nf_inneriter, 200
.set nb230nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb230nf_ix, 0
.set nb230nf_iy, 16
.set nb230nf_iz, 32
.set nb230nf_iq, 48
.set nb230nf_c6, 64
.set nb230nf_c12, 80
.set nb230nf_vctot, 96
.set nb230nf_Vvdwtot, 112
.set nb230nf_half, 128
.set nb230nf_three, 144
.set nb230nf_krf, 160
.set nb230nf_crf, 176
.set nb230nf_tsc, 192
.set nb230nf_nri, 208
.set nb230nf_iinr, 216
.set nb230nf_jindex, 224
.set nb230nf_jjnr, 232
.set nb230nf_shift, 240
.set nb230nf_shiftvec, 248
.set nb230nf_facel, 256
.set nb230nf_innerjjnr, 264
.set nb230nf_is3, 272
.set nb230nf_ii3, 280
.set nb230nf_ntia, 284
.set nb230nf_innerk, 288
.set nb230nf_n, 292
.set nb230nf_nn1, 296
.set nb230nf_ntype, 300
.set nb230nf_nouter, 304
.set nb230nf_ninner, 308


        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $328,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb230nf_nouter(%rsp)
        movl %eax,nb230nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb230nf_nri(%rsp)
        movq %rsi,nb230nf_iinr(%rsp)
        movq %rdx,nb230nf_jindex(%rsp)
        movq %rcx,nb230nf_jjnr(%rsp)
        movq %r8,nb230nf_shift(%rsp)
        movq %r9,nb230nf_shiftvec(%rsp)
        movq nb230nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb230nf_ntype(%rsp)
        movq nb230nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb230nf_facel(%rsp)

        movq nb230nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb230nf_tsc(%rsp)

        movq nb230nf_argkrf(%rbp),%rsi
        movq nb230nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb230nf_krf(%rsp)
        movapd %xmm2,nb230nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb230nf_half(%rsp)
        movl %ebx,nb230nf_half+4(%rsp)
        movsd nb230nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb230nf_half(%rsp)
        movapd %xmm3,nb230nf_three(%rsp)

_nb_kernel230nf_x86_64_sse2.nb230nf_threadloop: 
        movq  nb230nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel230nf_x86_64_sse2.nb230nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel230nf_x86_64_sse2.nb230nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb230nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb230nf_n(%rsp)
        movl %ebx,nb230nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel230nf_x86_64_sse2.nb230nf_outerstart
        jmp _nb_kernel230nf_x86_64_sse2.nb230nf_end

_nb_kernel230nf_x86_64_sse2.nb230nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb230nf_nouter(%rsp),%ebx
        movl %ebx,nb230nf_nouter(%rsp)

_nb_kernel230nf_x86_64_sse2.nb230nf_outer: 
        movq  nb230nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb230nf_is3(%rsp)            ## store is3 

        movq  nb230nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb230nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb230nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb230nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb230nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb230nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb230nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb230nf_pos(%rbp),%rax      ## eax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb230nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb230nf_ix(%rsp)
        movapd %xmm1,nb230nf_iy(%rsp)
        movapd %xmm2,nb230nf_iz(%rsp)

        movl  %ebx,nb230nf_ii3(%rsp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb230nf_vctot(%rsp)
        movapd %xmm4,nb230nf_Vvdwtot(%rsp)

        movq  nb230nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb230nf_pos(%rbp),%rsi
        movq  nb230nf_faction(%rbp),%rdi
        movq  nb230nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb230nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb230nf_ninner(%rsp),%ecx
        movl  %ecx,nb230nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb230nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel230nf_x86_64_sse2.nb230nf_unroll_loop
        jmp   _nb_kernel230nf_x86_64_sse2.nb230nf_checksingle
_nb_kernel230nf_x86_64_sse2.nb230nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb230nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb230nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb230nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb230nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb230nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb230nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb230nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movlpd (%rsi,%rbx,8),%xmm7      ## c6b
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        movhpd 8(%rsi,%rbx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb230nf_c6(%rsp)
        movapd %xmm6,nb230nf_c12(%rsp)

        movq nb230nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb230nf_ix(%rsp),%xmm4
        movapd nb230nf_iy(%rsp),%xmm5
        movapd nb230nf_iz(%rsp),%xmm6

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

        movapd nb230nf_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb230nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb230nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb230nf_crf(%rsp),%xmm6
        mulpd  %xmm3,%xmm6

        addpd  nb230nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb230nf_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb230nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movq nb230nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%rsi,%rax,8),%xmm4      ## Y1   
        movlpd (%rsi,%rbx,8),%xmm3      ## Y2 
        movhpd 8(%rsi,%rax,8),%xmm4     ## Y1 F1        
        movhpd 8(%rsi,%rbx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rax,8),%xmm6    ## G1
        movlpd 16(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 24(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 24(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb230nf_c6(%rsp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ##Update Vvdwtot directly 
        addpd  nb230nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb230nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rax,8),%xmm4    ## Y1   
        movlpd 32(%rsi,%rbx,8),%xmm3    ## Y2 
        movhpd 40(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movhpd 40(%rsi,%rbx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rax,8),%xmm6    ## G1
        movlpd 48(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 56(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 56(%rsi,%rbx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb230nf_c12(%rsp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb230nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb230nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb230nf_innerk(%rsp)
        jl    _nb_kernel230nf_x86_64_sse2.nb230nf_checksingle
        jmp   _nb_kernel230nf_x86_64_sse2.nb230nf_unroll_loop

_nb_kernel230nf_x86_64_sse2.nb230nf_checksingle: 
        movl  nb230nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel230nf_x86_64_sse2.nb230nf_dosingle
        jmp    _nb_kernel230nf_x86_64_sse2.nb230nf_updateouterdata
_nb_kernel230nf_x86_64_sse2.nb230nf_dosingle: 
        movq nb230nf_charge(%rbp),%rsi
        movq nb230nf_pos(%rbp),%rdi
        movq  nb230nf_innerjjnr(%rsp),%rcx
        xorpd %xmm3,%xmm3
        movl  (%rcx),%eax

        movlpd (%rsi,%rax,8),%xmm3
        movapd nb230nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movq nb230nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb230nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb230nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb230nf_c6(%rsp)
        movapd %xmm6,nb230nf_c12(%rsp)

        movq nb230nf_pos(%rbp),%rsi        ## base of pos[] 

        lea (%rax,%rax,2),%rax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb230nf_ix(%rsp),%xmm4
        movapd nb230nf_iy(%rsp),%xmm5
        movapd nb230nf_iz(%rsp),%xmm6

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

        movapd nb230nf_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb230nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb230nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb230nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb230nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb230nf_crf(%rsp),%xmm6
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 

        addsd  nb230nf_vctot(%rsp),%xmm6
        movsd %xmm6,nb230nf_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, cmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb230nf_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll  $3,%ebx

        movq nb230nf_VFtab(%rbp),%rsi

        ## dispersion 
        movlpd (%rsi,%rbx,8),%xmm4      ## Y1   
        movhpd 8(%rsi,%rbx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 24(%rsi,%rbx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb230nf_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb230nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb230nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rbx,8),%xmm4    ## Y1   
        movhpd 40(%rsi,%rbx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 56(%rsi,%rbx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb230nf_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb230nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb230nf_Vvdwtot(%rsp)

_nb_kernel230nf_x86_64_sse2.nb230nf_updateouterdata: 
        ## get n from stack
        movl nb230nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb230nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb230nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb230nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb230nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb230nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb230nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel230nf_x86_64_sse2.nb230nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb230nf_n(%rsp)
        jmp _nb_kernel230nf_x86_64_sse2.nb230nf_outer
_nb_kernel230nf_x86_64_sse2.nb230nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb230nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel230nf_x86_64_sse2.nb230nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel230nf_x86_64_sse2.nb230nf_threadloop
_nb_kernel230nf_x86_64_sse2.nb230nf_end: 
        movl nb230nf_nouter(%rsp),%eax
        movl nb230nf_ninner(%rsp),%ebx
        movq nb230nf_outeriter(%rbp),%rcx
        movq nb230nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $328,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

