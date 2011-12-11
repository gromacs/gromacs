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






.globl nb_kernel130_x86_64_sse2
.globl _nb_kernel130_x86_64_sse2
nb_kernel130_x86_64_sse2:       
_nb_kernel130_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb130_fshift, 16
.set nb130_gid, 24
.set nb130_pos, 32
.set nb130_faction, 40
.set nb130_charge, 48
.set nb130_p_facel, 56
.set nb130_argkrf, 64
.set nb130_argcrf, 72
.set nb130_Vc, 80
.set nb130_type, 88
.set nb130_p_ntype, 96
.set nb130_vdwparam, 104
.set nb130_Vvdw, 112
.set nb130_p_tabscale, 120
.set nb130_VFtab, 128
.set nb130_invsqrta, 136
.set nb130_dvda, 144
.set nb130_p_gbtabscale, 152
.set nb130_GBtab, 160
.set nb130_p_nthreads, 168
.set nb130_count, 176
.set nb130_mtx, 184
.set nb130_outeriter, 192
.set nb130_inneriter, 200
.set nb130_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb130_ix, 0
.set nb130_iy, 16
.set nb130_iz, 32
.set nb130_iq, 48
.set nb130_dx, 64
.set nb130_dy, 80
.set nb130_dz, 96
.set nb130_c6, 112
.set nb130_c12, 128
.set nb130_tsc, 144
.set nb130_fstmp, 160
.set nb130_vctot, 176
.set nb130_Vvdwtot, 192
.set nb130_fix, 208
.set nb130_fiy, 224
.set nb130_fiz, 240
.set nb130_half, 256
.set nb130_three, 272
.set nb130_two, 288
.set nb130_krf, 304
.set nb130_crf, 320
.set nb130_nri, 336
.set nb130_iinr, 344
.set nb130_jindex, 352
.set nb130_jjnr, 360
.set nb130_shift, 368
.set nb130_shiftvec, 376
.set nb130_facel, 384
.set nb130_innerjjnr, 392
.set nb130_is3, 400
.set nb130_ii3, 404
.set nb130_ntia, 408
.set nb130_innerk, 412
.set nb130_n, 416
.set nb130_nn1, 420
.set nb130_ntype, 424
.set nb130_nouter, 428
.set nb130_ninner, 432
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
        movl %eax,nb130_nouter(%rsp)
        movl %eax,nb130_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb130_nri(%rsp)
        movq %rsi,nb130_iinr(%rsp)
        movq %rdx,nb130_jindex(%rsp)
        movq %rcx,nb130_jjnr(%rsp)
        movq %r8,nb130_shift(%rsp)
        movq %r9,nb130_shiftvec(%rsp)
        movq nb130_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb130_ntype(%rsp)
        movq nb130_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb130_facel(%rsp)

        movq nb130_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb130_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb130_half(%rsp)
        movl %ebx,nb130_half+4(%rsp)
        movsd nb130_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb130_half(%rsp)
        movapd %xmm2,nb130_two(%rsp)
        movapd %xmm3,nb130_three(%rsp)

_nb_kernel130_x86_64_sse2.nb130_threadloop: 
        movq  nb130_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel130_x86_64_sse2.nb130_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel130_x86_64_sse2.nb130_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb130_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb130_n(%rsp)
        movl %ebx,nb130_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel130_x86_64_sse2.nb130_outerstart
        jmp _nb_kernel130_x86_64_sse2.nb130_end

_nb_kernel130_x86_64_sse2.nb130_outerstart: 
        ## ebx contains number of outer iterations
        addl nb130_nouter(%rsp),%ebx
        movl %ebx,nb130_nouter(%rsp)

_nb_kernel130_x86_64_sse2.nb130_outer: 
        movq  nb130_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb130_is3(%rsp)      ## store is3 

        movq  nb130_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb130_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb130_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb130_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb130_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb130_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb130_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb130_pos(%rbp),%rax      ## eax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb130_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb130_ix(%rsp)
        movapd %xmm1,nb130_iy(%rsp)
        movapd %xmm2,nb130_iz(%rsp)

        movl  %ebx,nb130_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb130_vctot(%rsp)
        movapd %xmm4,nb130_Vvdwtot(%rsp)
        movapd %xmm4,nb130_fix(%rsp)
        movapd %xmm4,nb130_fiy(%rsp)
        movapd %xmm4,nb130_fiz(%rsp)

        movq  nb130_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb130_pos(%rbp),%rsi
        movq  nb130_faction(%rbp),%rdi
        movq  nb130_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb130_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb130_ninner(%rsp),%ecx
        movl  %ecx,nb130_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb130_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel130_x86_64_sse2.nb130_unroll_loop
        jmp   _nb_kernel130_x86_64_sse2.nb130_checksingle
_nb_kernel130_x86_64_sse2.nb130_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb130_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb130_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb130_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        mulpd nb130_iq(%rsp),%xmm3              ## qq 

        movq nb130_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb130_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb130_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6
        movhpd (%rsi,%r9,8),%xmm6
        movlpd 8(%rsi,%r8,8),%xmm7      ## c12 
        movhpd 8(%rsi,%r9,8),%xmm7
        movapd %xmm6,nb130_c6(%rsp)
        movapd %xmm7,nb130_c12(%rsp)

        movq nb130_pos(%rbp),%rsi        ## base of pos[] 

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
        subpd nb130_ix(%rsp),%xmm4
        subpd nb130_iy(%rsp),%xmm5
        subpd nb130_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb130_dx(%rsp)
        movapd %xmm5,nb130_dy(%rsp)
        movapd %xmm6,nb130_dz(%rsp)
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
        movapd nb130_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb130_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        movapd %xmm4,%xmm1
        mulpd  %xmm3,%xmm6  ## vcoul = rinv*qq
        movapd %xmm6,%xmm3
        mulpd  %xmm0,%xmm3

    ## fstmp in xmm3

        addpd  nb130_vctot(%rsp),%xmm6
        movapd %xmm6,nb130_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb130_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 

        pslld $3,%mm6           ## idx *= 8 

        movq nb130_VFtab(%rbp),%rsi
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
    movapd nb130_c6(%rsp),%xmm12
    movapd nb130_c12(%rsp),%xmm13
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb130_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb130_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7

    mulpd  nb130_tsc(%rsp),%xmm7
    subpd  %xmm7,%xmm3
    mulpd  %xmm0,%xmm3  ## fscal

    movapd %xmm3,%xmm9
    movapd %xmm3,%xmm10
    movapd %xmm3,%xmm11

    movapd nb130_fix(%rsp),%xmm12
    movapd nb130_fiy(%rsp),%xmm13
    movapd nb130_fiz(%rsp),%xmm14

    mulpd  nb130_dx(%rsp),%xmm9
    mulpd  nb130_dy(%rsp),%xmm10
    mulpd  nb130_dz(%rsp),%xmm11

    ## accumulate i forces
    addpd %xmm9,%xmm12
    addpd %xmm10,%xmm13
    addpd %xmm11,%xmm14
    movapd %xmm12,nb130_fix(%rsp)
    movapd %xmm13,nb130_fiy(%rsp)
    movapd %xmm14,nb130_fiz(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb130_faction(%rbp),%rdi
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
        subl $2,nb130_innerk(%rsp)
        jl    _nb_kernel130_x86_64_sse2.nb130_checksingle
        jmp   _nb_kernel130_x86_64_sse2.nb130_unroll_loop

_nb_kernel130_x86_64_sse2.nb130_checksingle:    
        movl  nb130_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel130_x86_64_sse2.nb130_dosingle
        jmp    _nb_kernel130_x86_64_sse2.nb130_updateouterdata
_nb_kernel130_x86_64_sse2.nb130_dosingle: 
        movq nb130_charge(%rbp),%rsi
        movq nb130_pos(%rbp),%rdi
        movq  nb130_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb130_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
        mulsd nb130_iq(%rsp),%xmm3              ## qq 

        movq nb130_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb130_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb130_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb130_c6(%rsp)
        movapd %xmm6,nb130_c12(%rsp)

        movq nb130_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinate to xmm4-xmm6         
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb130_ix(%rsp),%xmm4
        subsd nb130_iy(%rsp),%xmm5
        subsd nb130_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb130_dx(%rsp)
        movapd %xmm5,nb130_dy(%rsp)
        movapd %xmm6,nb130_dz(%rsp)
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
        movapd nb130_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb130_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        movapd %xmm4,%xmm1
        mulsd  %xmm3,%xmm6  ## vcoul = rinv*qq
        movapd %xmm6,%xmm3
        mulsd  %xmm0,%xmm3

    ## fstmp in xmm3

        addsd  nb130_vctot(%rsp),%xmm6
        movsd %xmm6,nb130_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb130_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 

        shll    $3,%r8d         ## idx *= 8 

        movq nb130_VFtab(%rbp),%rsi

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
    movapd nb130_c6(%rsp),%xmm12
    movapd nb130_c12(%rsp),%xmm13
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulsd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb130_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb130_Vvdwtot(%rsp)

    mulsd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulsd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7

    mulsd  nb130_tsc(%rsp),%xmm7
    subsd  %xmm7,%xmm3
    mulsd  %xmm0,%xmm3  ## fscal

    movapd %xmm3,%xmm9
    movapd %xmm3,%xmm10
    movapd %xmm3,%xmm11

    movapd nb130_fix(%rsp),%xmm12
    movapd nb130_fiy(%rsp),%xmm13
    movapd nb130_fiz(%rsp),%xmm14

    mulsd  nb130_dx(%rsp),%xmm9
    mulsd  nb130_dy(%rsp),%xmm10
    mulsd  nb130_dz(%rsp),%xmm11

    ## accumulate i forces
    addsd %xmm9,%xmm12
    addsd %xmm10,%xmm13
    addsd %xmm11,%xmm14
    movsd %xmm12,nb130_fix(%rsp)
    movsd %xmm13,nb130_fiy(%rsp)
    movsd %xmm14,nb130_fiz(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb130_faction(%rbp),%rdi
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel130_x86_64_sse2.nb130_updateouterdata: 
        movl  nb130_ii3(%rsp),%ecx
        movq  nb130_faction(%rbp),%rdi
        movq  nb130_fshift(%rbp),%rsi
        movl  nb130_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb130_fix(%rsp),%xmm0
        movapd nb130_fiy(%rsp),%xmm1
        movapd nb130_fiz(%rsp),%xmm2

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
        movl nb130_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb130_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb130_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb130_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb130_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb130_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb130_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel130_x86_64_sse2.nb130_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb130_n(%rsp)
        jmp _nb_kernel130_x86_64_sse2.nb130_outer
_nb_kernel130_x86_64_sse2.nb130_outerend: 
        ## check if more outer neighborlists remain
        movl  nb130_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel130_x86_64_sse2.nb130_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel130_x86_64_sse2.nb130_threadloop
_nb_kernel130_x86_64_sse2.nb130_end: 
        movl nb130_nouter(%rsp),%eax
        movl nb130_ninner(%rsp),%ebx
        movq nb130_outeriter(%rbp),%rcx
        movq nb130_inneriter(%rbp),%rdx
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




.globl nb_kernel130nf_x86_64_sse2
.globl _nb_kernel130nf_x86_64_sse2
nb_kernel130nf_x86_64_sse2:     
_nb_kernel130nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb130nf_fshift, 16
.set nb130nf_gid, 24
.set nb130nf_pos, 32
.set nb130nf_faction, 40
.set nb130nf_charge, 48
.set nb130nf_p_facel, 56
.set nb130nf_argkrf, 64
.set nb130nf_argcrf, 72
.set nb130nf_Vc, 80
.set nb130nf_type, 88
.set nb130nf_p_ntype, 96
.set nb130nf_vdwparam, 104
.set nb130nf_Vvdw, 112
.set nb130nf_p_tabscale, 120
.set nb130nf_VFtab, 128
.set nb130nf_invsqrta, 136
.set nb130nf_dvda, 144
.set nb130nf_p_gbtabscale, 152
.set nb130nf_GBtab, 160
.set nb130nf_p_nthreads, 168
.set nb130nf_count, 176
.set nb130nf_mtx, 184
.set nb130nf_outeriter, 192
.set nb130nf_inneriter, 200
.set nb130nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb130nf_ix, 0
.set nb130nf_iy, 16
.set nb130nf_iz, 32
.set nb130nf_iq, 48
.set nb130nf_c6, 64
.set nb130nf_c12, 80
.set nb130nf_vctot, 96
.set nb130nf_Vvdwtot, 112
.set nb130nf_half, 128
.set nb130nf_three, 144
.set nb130nf_krf, 160
.set nb130nf_crf, 176
.set nb130nf_tsc, 192
.set nb130nf_nri, 208
.set nb130nf_iinr, 216
.set nb130nf_jindex, 224
.set nb130nf_jjnr, 232
.set nb130nf_shift, 240
.set nb130nf_shiftvec, 248
.set nb130nf_facel, 256
.set nb130nf_innerjjnr, 264
.set nb130nf_is3, 272
.set nb130nf_ii3, 280
.set nb130nf_ntia, 284
.set nb130nf_innerk, 288
.set nb130nf_n, 292
.set nb130nf_nn1, 296
.set nb130nf_ntype, 300
.set nb130nf_nouter, 304
.set nb130nf_ninner, 308

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
        movl %eax,nb130nf_nouter(%rsp)
        movl %eax,nb130nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb130nf_nri(%rsp)
        movq %rsi,nb130nf_iinr(%rsp)
        movq %rdx,nb130nf_jindex(%rsp)
        movq %rcx,nb130nf_jjnr(%rsp)
        movq %r8,nb130nf_shift(%rsp)
        movq %r9,nb130nf_shiftvec(%rsp)
        movq nb130nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb130nf_ntype(%rsp)
        movq nb130nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb130nf_facel(%rsp)

        movq nb130nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb130nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb130nf_half(%rsp)
        movl %ebx,nb130nf_half+4(%rsp)
        movsd nb130nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb130nf_half(%rsp)
        movapd %xmm3,nb130nf_three(%rsp)

_nb_kernel130nf_x86_64_sse2.nb130nf_threadloop: 
        movq  nb130nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel130nf_x86_64_sse2.nb130nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel130nf_x86_64_sse2.nb130nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb130nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb130nf_n(%rsp)
        movl %ebx,nb130nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel130nf_x86_64_sse2.nb130nf_outerstart
        jmp _nb_kernel130nf_x86_64_sse2.nb130nf_end

_nb_kernel130nf_x86_64_sse2.nb130nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb130nf_nouter(%rsp),%ebx
        movl %ebx,nb130nf_nouter(%rsp)

_nb_kernel130nf_x86_64_sse2.nb130nf_outer: 
        movq  nb130nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb130nf_is3(%rsp)            ## store is3 

        movq  nb130nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb130nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb130nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb130nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb130nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb130nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb130nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb130nf_pos(%rbp),%rax      ## eax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb130nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb130nf_ix(%rsp)
        movapd %xmm1,nb130nf_iy(%rsp)
        movapd %xmm2,nb130nf_iz(%rsp)

        movl  %ebx,nb130nf_ii3(%rsp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb130nf_vctot(%rsp)
        movapd %xmm4,nb130nf_Vvdwtot(%rsp)

        movq  nb130nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb130nf_pos(%rbp),%rsi
        movq  nb130nf_faction(%rbp),%rdi
        movq  nb130nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb130nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb130nf_ninner(%rsp),%ecx
        movl  %ecx,nb130nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb130nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel130nf_x86_64_sse2.nb130nf_unroll_loop
        jmp   _nb_kernel130nf_x86_64_sse2.nb130nf_checksingle
_nb_kernel130nf_x86_64_sse2.nb130nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb130nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb130nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb130nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb130nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb130nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb130nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb130nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb130nf_c6(%rsp)
        movapd %xmm6,nb130nf_c12(%rsp)

        movq nb130nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb130nf_ix(%rsp),%xmm4
        movapd nb130nf_iy(%rsp),%xmm5
        movapd nb130nf_iz(%rsp),%xmm6

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
        movapd nb130nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb130nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        mulpd  %xmm3,%xmm6  ## vcoul

        addpd  nb130nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb130nf_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb130nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movq nb130nf_VFtab(%rbp),%rsi
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

        movapd nb130nf_c6(%rsp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ##Update Vvdwtot directly 
        addpd  nb130nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb130nf_Vvdwtot(%rsp)

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

        movapd nb130nf_c12(%rsp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb130nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb130nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb130nf_innerk(%rsp)
        jl    _nb_kernel130nf_x86_64_sse2.nb130nf_checksingle
        jmp   _nb_kernel130nf_x86_64_sse2.nb130nf_unroll_loop

_nb_kernel130nf_x86_64_sse2.nb130nf_checksingle: 
        movl  nb130nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel130nf_x86_64_sse2.nb130nf_dosingle
        jmp    _nb_kernel130nf_x86_64_sse2.nb130nf_updateouterdata
_nb_kernel130nf_x86_64_sse2.nb130nf_dosingle: 
        movq nb130nf_charge(%rbp),%rsi
        movq nb130nf_pos(%rbp),%rdi
        movq  nb130nf_innerjjnr(%rsp),%rcx
        xorpd %xmm3,%xmm3
        movl  (%rcx),%eax

        movlpd (%rsi,%rax,8),%xmm3
        movapd nb130nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movq nb130nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb130nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb130nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb130nf_c6(%rsp)
        movapd %xmm6,nb130nf_c12(%rsp)

        movq nb130nf_pos(%rbp),%rsi        ## base of pos[] 

        lea (%rax,%rax,2),%rax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb130nf_ix(%rsp),%xmm4
        movapd nb130nf_iy(%rsp),%xmm5
        movapd nb130nf_iz(%rsp),%xmm6

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
        movapd nb130nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb130nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb130nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb130nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm6
        movapd %xmm4,%xmm1
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*rinv

        addsd  nb130nf_vctot(%rsp),%xmm6
        movsd %xmm6,nb130nf_vctot(%rsp)

        ## LJ table interaction. xmm0=rinv, cmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb130nf_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll  $3,%ebx

        movq nb130nf_VFtab(%rbp),%rsi

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

        movsd nb130nf_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb130nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb130nf_Vvdwtot(%rsp)

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

        movsd nb130nf_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb130nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb130nf_Vvdwtot(%rsp)

_nb_kernel130nf_x86_64_sse2.nb130nf_updateouterdata: 
        ## get n from stack
        movl nb130nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb130nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb130nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb130nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb130nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb130nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb130nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel130nf_x86_64_sse2.nb130nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb130nf_n(%rsp)
        jmp _nb_kernel130nf_x86_64_sse2.nb130nf_outer
_nb_kernel130nf_x86_64_sse2.nb130nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb130nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel130nf_x86_64_sse2.nb130nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel130nf_x86_64_sse2.nb130nf_threadloop
_nb_kernel130nf_x86_64_sse2.nb130nf_end: 
        movl nb130nf_nouter(%rsp),%eax
        movl nb130nf_ninner(%rsp),%ebx
        movq nb130nf_outeriter(%rbp),%rcx
        movq nb130nf_inneriter(%rbp),%rdx
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

