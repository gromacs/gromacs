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








.globl nb_kernel110_x86_64_sse2
.globl _nb_kernel110_x86_64_sse2
nb_kernel110_x86_64_sse2:       
_nb_kernel110_x86_64_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
.set nb110_ix, 0
.set nb110_iy, 16
.set nb110_iz, 32
.set nb110_iq, 48
.set nb110_dx, 64
.set nb110_dy, 80
.set nb110_dz, 96
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
.set nb110_is3, 288
.set nb110_ii3, 292
.set nb110_nri, 296
.set nb110_iinr, 304
.set nb110_jindex, 312
.set nb110_jjnr, 320
.set nb110_shift, 328
.set nb110_shiftvec, 336
.set nb110_facel, 344
.set nb110_innerjjnr, 352
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

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $408,%rsp          ## local variable stack space (n*16+8)

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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb110_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb110_half(%rsp)
        movl %ebx,nb110_half+4(%rsp)
        movsd nb110_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## six
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## twelve
        movapd %xmm1,nb110_half(%rsp)
        movapd %xmm3,nb110_three(%rsp)
        movapd %xmm4,nb110_six(%rsp)
        movapd %xmm5,nb110_twelve(%rsp)

_nb_kernel110_x86_64_sse2.nb110_threadloop: 
        movq  nb110_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel110_x86_64_sse2.nb110_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110_x86_64_sse2.nb110_spinlock

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
        jg  _nb_kernel110_x86_64_sse2.nb110_outerstart
        jmp _nb_kernel110_x86_64_sse2.nb110_end

_nb_kernel110_x86_64_sse2.nb110_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110_nouter(%rsp),%ebx
        movl %ebx,nb110_nouter(%rsp)

_nb_kernel110_x86_64_sse2.nb110_outer: 
        movq  nb110_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb110_is3(%rsp)      ## store is3 

        movq  nb110_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb110_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb110_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb110_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb110_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb110_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb110_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb110_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb110_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb110_ix(%rsp)
        movapd %xmm1,nb110_iy(%rsp)
        movapd %xmm2,nb110_iz(%rsp)

        movl  %ebx,nb110_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm12,%xmm12
        movapd %xmm12,nb110_Vvdwtot(%rsp)
        movapd %xmm12,%xmm13
        movapd %xmm12,%xmm14
        movapd %xmm12,%xmm15

        movq  nb110_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb110_pos(%rbp),%rsi
        movq  nb110_faction(%rbp),%rdi
        movq  nb110_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb110_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb110_ninner(%rsp),%ecx
        movl  %ecx,nb110_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb110_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel110_x86_64_sse2.nb110_unroll_loop
        jmp   _nb_kernel110_x86_64_sse2.nb110_checksingle
_nb_kernel110_x86_64_sse2.nb110_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb110_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        addq $8,nb110_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb110_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx

        ## move two coordinates to xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        ## calc dr 
        subpd nb110_ix(%rsp),%xmm4
        subpd nb110_iy(%rsp),%xmm5
        subpd nb110_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        movq nb110_charge(%rbp),%rsi     ## base of charge[] 

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        movq nb110_type(%rbp),%rdi
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        movlpd (%rsi,%r8,8),%xmm3
        movhpd (%rsi,%r9,8),%xmm3

        cvtpd2ps %xmm4,%xmm5
        movl (%rdi,%r8,4),%r8d
        movl (%rdi,%r9,4),%r9d
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        shll %r8d
        shll %r9d

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movl nb110_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        movapd nb110_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movq nb110_vdwparam(%rbp),%rdi
        movapd nb110_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 
        mulpd nb110_iq(%rsp),%xmm3              ## qq 

    movlpd (%rdi,%r8,8),%xmm6
    movlpd 8(%rdi,%r8,8),%xmm7

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb110_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

    movhpd (%rdi,%r9,8),%xmm6
    movhpd 8(%rdi,%r9,8),%xmm7

        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulpd  %xmm6,%xmm1
        mulpd  %xmm7,%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb110_Vvdwtot(%rsp),%xmm5
        mulpd  nb110_six(%rsp),%xmm1
        mulpd  nb110_twelve(%rsp),%xmm2
        subpd  %xmm1,%xmm2
        addpd  %xmm3,%xmm2
        mulpd  %xmm2,%xmm4      ## xmm4=total fscal 
        addpd  %xmm3,%xmm12 ## add to vctot
        movq   nb110_faction(%rbp),%rdi

        ## the fj's - start by accumulating forces from memory 
        movlpd (%rdi,%rax,8),%xmm6
        movlpd 8(%rdi,%rax,8),%xmm7
        movlpd 16(%rdi,%rax,8),%xmm8
        movhpd (%rdi,%rbx,8),%xmm6
        movhpd 8(%rdi,%rbx,8),%xmm7
        movhpd 16(%rdi,%rbx,8),%xmm8

        movapd %xmm5,nb110_Vvdwtot(%rsp)

        mulpd  %xmm4,%xmm9
        mulpd  %xmm4,%xmm10
        mulpd  %xmm4,%xmm11

        addpd %xmm9,%xmm6
        addpd %xmm10,%xmm7
        addpd %xmm11,%xmm8

        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15

        movlpd %xmm6,(%rdi,%rax,8)
        movlpd %xmm7,8(%rdi,%rax,8)
        movlpd %xmm8,16(%rdi,%rax,8)
        movhpd %xmm6,(%rdi,%rbx,8)
        movhpd %xmm7,8(%rdi,%rbx,8)
        movhpd %xmm8,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb110_innerk(%rsp)
        jl    _nb_kernel110_x86_64_sse2.nb110_checksingle
        jmp   _nb_kernel110_x86_64_sse2.nb110_unroll_loop
_nb_kernel110_x86_64_sse2.nb110_checksingle: 
        movl  nb110_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel110_x86_64_sse2.nb110_dosingle
        jmp    _nb_kernel110_x86_64_sse2.nb110_updateouterdata
_nb_kernel110_x86_64_sse2.nb110_dosingle: 
        movq nb110_charge(%rbp),%rsi
        movq nb110_pos(%rbp),%rdi
        movq nb110_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb110_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
    mulsd nb110_iq(%rsp),%xmm3   ## qq

        movq nb110_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb110_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb110_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4            ## c6
        movsd 8(%rsi,%r8,8),%xmm6       ## c12
        movapd %xmm4,nb110_c6(%rsp)
        movapd %xmm6,nb110_c12(%rsp)

        movq nb110_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb110_ix(%rsp),%xmm4
        subsd nb110_iy(%rsp),%xmm5
        subsd nb110_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

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
        movapd nb110_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb110_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb110_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulsd  nb110_c6(%rsp),%xmm1
        mulsd  nb110_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb110_Vvdwtot(%rsp),%xmm5
        mulsd  nb110_six(%rsp),%xmm1
        mulsd  nb110_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2
        addsd  %xmm3,%xmm2
        mulsd  %xmm2,%xmm4      ## xmm4=total fscal 
        addsd  %xmm3,%xmm12  ## add to vctot

        movsd %xmm5,nb110_Vvdwtot(%rsp)

        movq   nb110_faction(%rbp),%rdi
        mulsd  %xmm4,%xmm9
        mulsd  %xmm4,%xmm10
        mulsd  %xmm4,%xmm11

        ## now update f_i 
        addsd  %xmm9,%xmm13
        addsd  %xmm10,%xmm14
        addsd  %xmm11,%xmm15
        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel110_x86_64_sse2.nb110_updateouterdata: 
        movl  nb110_ii3(%rsp),%ecx
        movq  nb110_faction(%rbp),%rdi
        movq  nb110_fshift(%rbp),%rsi
        movl  nb110_is3(%rsp),%edx

        ## accumulate i forces in xmm13, xmm14, xmm15
        movhlps %xmm13,%xmm3
        movhlps %xmm14,%xmm4
        movhlps %xmm15,%xmm5
        addsd  %xmm3,%xmm13
        addsd  %xmm4,%xmm14
        addsd  %xmm5,%xmm15 ## sum is in low xmm13-xmm15

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## increment fshift force  
        movsd  (%rsi,%rdx,8),%xmm3
        movsd  8(%rsi,%rdx,8),%xmm4
        movsd  16(%rsi,%rdx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rsi,%rdx,8)
        movsd  %xmm4,8(%rsi,%rdx,8)
        movsd  %xmm5,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb110_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb110_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movhlps %xmm12,%xmm6
        addsd  %xmm6,%xmm12     ## low xmm12 has the sum now 

        ## add earlier value from mem 
        movq  nb110_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm12
        ## move back to mem 
        movsd %xmm12,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb110_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb110_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb110_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110_x86_64_sse2.nb110_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110_n(%rsp)
        jmp _nb_kernel110_x86_64_sse2.nb110_outer
_nb_kernel110_x86_64_sse2.nb110_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110_x86_64_sse2.nb110_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110_x86_64_sse2.nb110_threadloop
_nb_kernel110_x86_64_sse2.nb110_end: 
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






.globl nb_kernel110nf_x86_64_sse2
.globl _nb_kernel110nf_x86_64_sse2
nb_kernel110nf_x86_64_sse2:     
_nb_kernel110nf_x86_64_sse2:    
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb110nf_is3, 160
.set nb110nf_ii3, 164
.set nb110nf_nri, 168
.set nb110nf_iinr, 176
.set nb110nf_jindex, 184
.set nb110nf_jjnr, 192
.set nb110nf_shift, 200
.set nb110nf_shiftvec, 208
.set nb110nf_facel, 216
.set nb110nf_innerjjnr, 224
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

        push %r12
        push %r13
        push %r14
        push %r15

        subq $280,%rsp          ## local variable stack space (n*16+8)

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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb110nf_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb110nf_half(%rsp)
        movl %ebx,nb110nf_half+4(%rsp)
        movsd nb110nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb110nf_half(%rsp)
        movapd %xmm3,nb110nf_three(%rsp)

_nb_kernel110nf_x86_64_sse2.nb110nf_threadloop: 
        movq  nb110nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel110nf_x86_64_sse2.nb110nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110nf_x86_64_sse2.nb110nf_spinlock

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
        jg  _nb_kernel110nf_x86_64_sse2.nb110nf_outerstart
        jmp _nb_kernel110nf_x86_64_sse2.nb110nf_end

_nb_kernel110nf_x86_64_sse2.nb110nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110nf_nouter(%rsp),%ebx
        movl %ebx,nb110nf_nouter(%rsp)

_nb_kernel110nf_x86_64_sse2.nb110nf_outer: 
        movq  nb110nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb110nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb110nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb110nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb110nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb110nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb110nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb110nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb110nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb110nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb110nf_ix(%rsp)
        movapd %xmm1,nb110nf_iy(%rsp)
        movapd %xmm2,nb110nf_iz(%rsp)

        movl  %ebx,nb110nf_ii3(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb110nf_vctot(%rsp)
        movapd %xmm4,nb110nf_Vvdwtot(%rsp)

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
        subl  $2,%edx
        addl  nb110nf_ninner(%rsp),%ecx
        movl  %ecx,nb110nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb110nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel110nf_x86_64_sse2.nb110nf_unroll_loop
        jmp   _nb_kernel110nf_x86_64_sse2.nb110nf_checksingle
_nb_kernel110nf_x86_64_sse2.nb110nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb110nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb110nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb110nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb110nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb110nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb110nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb110nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb110nf_c6(%rsp)
        movapd %xmm6,nb110nf_c12(%rsp)

        movq nb110nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb110nf_ix(%rsp),%xmm4
        movapd nb110nf_iy(%rsp),%xmm5
        movapd nb110nf_iz(%rsp),%xmm6

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
        movapd nb110nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb110nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb110nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulpd  nb110nf_c6(%rsp),%xmm1
        mulpd  nb110nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb110nf_Vvdwtot(%rsp),%xmm5
        addpd  nb110nf_vctot(%rsp),%xmm3
        movapd %xmm3,nb110nf_vctot(%rsp)
        movapd %xmm5,nb110nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb110nf_innerk(%rsp)
        jl    _nb_kernel110nf_x86_64_sse2.nb110nf_checksingle
        jmp   _nb_kernel110nf_x86_64_sse2.nb110nf_unroll_loop
_nb_kernel110nf_x86_64_sse2.nb110nf_checksingle: 
        movl  nb110nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel110nf_x86_64_sse2.nb110nf_dosingle
        jmp   _nb_kernel110nf_x86_64_sse2.nb110nf_updateouterdata
_nb_kernel110nf_x86_64_sse2.nb110nf_dosingle: 
        movq nb110nf_charge(%rbp),%rsi
        movq nb110nf_pos(%rbp),%rdi
        movq nb110nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3

        movapd nb110nf_iq(%rsp),%xmm5
        mulsd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movq nb110nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb110nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb110nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb110nf_c6(%rsp)
        movapd %xmm6,nb110nf_c12(%rsp)

        movq nb110nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb110nf_ix(%rsp),%xmm4
        movapd nb110nf_iy(%rsp),%xmm5
        movapd nb110nf_iz(%rsp),%xmm6

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
        movapd nb110nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb110nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb110nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulsd  nb110nf_c6(%rsp),%xmm1
        mulsd  nb110nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb110nf_Vvdwtot(%rsp),%xmm5
        addsd  nb110nf_vctot(%rsp),%xmm3
        movlpd %xmm3,nb110nf_vctot(%rsp)
        movlpd %xmm5,nb110nf_Vvdwtot(%rsp)

_nb_kernel110nf_x86_64_sse2.nb110nf_updateouterdata: 
        ## get n from stack
        movl nb110nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb110nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb110nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb110nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb110nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb110nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb110nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110nf_x86_64_sse2.nb110nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110nf_n(%rsp)
        jmp _nb_kernel110nf_x86_64_sse2.nb110nf_outer
_nb_kernel110nf_x86_64_sse2.nb110nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110nf_x86_64_sse2.nb110nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110nf_x86_64_sse2.nb110nf_threadloop
_nb_kernel110nf_x86_64_sse2.nb110nf_end: 
        movl nb110nf_nouter(%rsp),%eax
        movl nb110nf_ninner(%rsp),%ebx
        movq nb110nf_outeriter(%rbp),%rcx
        movq nb110nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $280,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


