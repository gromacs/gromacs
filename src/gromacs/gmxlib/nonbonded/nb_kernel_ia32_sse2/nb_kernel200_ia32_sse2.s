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



.globl nb_kernel200_ia32_sse2
.globl _nb_kernel200_ia32_sse2
nb_kernel200_ia32_sse2: 
_nb_kernel200_ia32_sse2:        
.set nb200_p_nri, 8
.set nb200_iinr, 12
.set nb200_jindex, 16
.set nb200_jjnr, 20
.set nb200_shift, 24
.set nb200_shiftvec, 28
.set nb200_fshift, 32
.set nb200_gid, 36
.set nb200_pos, 40
.set nb200_faction, 44
.set nb200_charge, 48
.set nb200_p_facel, 52
.set nb200_argkrf, 56
.set nb200_argcrf, 60
.set nb200_Vc, 64
.set nb200_type, 68
.set nb200_p_ntype, 72
.set nb200_vdwparam, 76
.set nb200_Vvdw, 80
.set nb200_p_tabscale, 84
.set nb200_VFtab, 88
.set nb200_invsqrta, 92
.set nb200_dvda, 96
.set nb200_p_gbtabscale, 100
.set nb200_GBtab, 104
.set nb200_p_nthreads, 108
.set nb200_count, 112
.set nb200_mtx, 116
.set nb200_outeriter, 120
.set nb200_inneriter, 124
.set nb200_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb200_ix, 0
.set nb200_iy, 16
.set nb200_iz, 32
.set nb200_iq, 48
.set nb200_dx, 64
.set nb200_dy, 80
.set nb200_dz, 96
.set nb200_vctot, 112
.set nb200_fix, 128
.set nb200_fiy, 144
.set nb200_fiz, 160
.set nb200_half, 176
.set nb200_three, 192
.set nb200_two, 208
.set nb200_krf, 224
.set nb200_crf, 240
.set nb200_facel, 256                       ## uses 8 bytes
.set nb200_is3, 264
.set nb200_ii3, 268
.set nb200_innerjjnr, 272
.set nb200_innerk, 276
.set nb200_n, 280
.set nb200_nn1, 284
.set nb200_nri, 288
.set nb200_nouter, 292
.set nb200_ninner, 296
.set nb200_salign, 300
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $304,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb200_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb200_p_nri(%ebp),%ecx
        movl nb200_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb200_nri(%esp)
        movsd %xmm7,nb200_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb200_nouter(%esp)
        movl %eax,nb200_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb200_half(%esp)
        movl %ebx,nb200_half+4(%esp)
        movsd nb200_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb200_half(%esp)
        movapd %xmm2,nb200_two(%esp)
        movapd %xmm3,nb200_three(%esp)

        movl nb200_argkrf(%ebp),%esi
        movl nb200_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb200_krf(%esp)
        shufpd $0,%xmm6,%xmm6
        movapd %xmm6,nb200_crf(%esp)

_nb_kernel200_ia32_sse2.nb200_threadloop: 
        movl  nb200_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel200_ia32_sse2.nb200_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200_ia32_sse2.nb200_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200_n(%esp)
        movl %ebx,nb200_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200_ia32_sse2.nb200_outerstart
        jmp _nb_kernel200_ia32_sse2.nb200_end

_nb_kernel200_ia32_sse2.nb200_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200_nouter(%esp),%ebx
        movl %ebx,nb200_nouter(%esp)

_nb_kernel200_ia32_sse2.nb200_outer: 
        movl  nb200_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb200_is3(%esp)      ## store is3 

        movl  nb200_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb200_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb200_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb200_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb200_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb200_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb200_ix(%esp)
        movapd %xmm1,nb200_iy(%esp)
        movapd %xmm2,nb200_iz(%esp)

        movl  %ebx,nb200_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb200_vctot(%esp)
        movapd %xmm4,nb200_fix(%esp)
        movapd %xmm4,nb200_fiy(%esp)
        movapd %xmm4,nb200_fiz(%esp)

        movl  nb200_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb200_pos(%ebp),%esi
        movl  nb200_faction(%ebp),%edi
        movl  nb200_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb200_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb200_ninner(%esp),%ecx
        movl  %ecx,nb200_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb200_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel200_ia32_sse2.nb200_unroll_loop
        jmp   _nb_kernel200_ia32_sse2.nb200_checksingle
_nb_kernel200_ia32_sse2.nb200_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb200_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb200_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb200_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb200_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl nb200_pos(%ebp),%esi        ## base of pos[] 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb200_ix(%esp),%xmm4
        movapd nb200_iy(%esp),%xmm5
        movapd nb200_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb200_dx(%esp)
        movapd %xmm5,nb200_dy(%esp)
        movapd %xmm6,nb200_dz(%esp)
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

        movapd nb200_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb200_three(%esp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb200_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb200_crf(%esp),%xmm6
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 

        mulpd  nb200_two(%esp),%xmm7

        subpd  %xmm7,%xmm0
        mulpd  %xmm0,%xmm3
        mulpd  %xmm3,%xmm4      ## xmm4=total fscal 
        addpd  nb200_vctot(%esp),%xmm6

        movapd nb200_dx(%esp),%xmm0
        movapd nb200_dy(%esp),%xmm1
        movapd nb200_dz(%esp),%xmm2

        movapd %xmm6,nb200_vctot(%esp)

        movl   nb200_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb200_fix(%esp),%xmm3
        movapd nb200_fiy(%esp),%xmm4
        movapd nb200_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb200_fix(%esp)
        movapd %xmm4,nb200_fiy(%esp)
        movapd %xmm5,nb200_fiz(%esp)
        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        movhpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%ebx,8),%xmm4
        movhpd 16(%edi,%ebx,8),%xmm5
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
        movhpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm4,8(%edi,%ebx,8)
        movhpd %xmm5,16(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb200_innerk(%esp)
        jl    _nb_kernel200_ia32_sse2.nb200_checksingle
        jmp   _nb_kernel200_ia32_sse2.nb200_unroll_loop

_nb_kernel200_ia32_sse2.nb200_checksingle:      
        movl  nb200_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200_ia32_sse2.nb200_dosingle
        jmp    _nb_kernel200_ia32_sse2.nb200_updateouterdata
_nb_kernel200_ia32_sse2.nb200_dosingle: 
        movl nb200_charge(%ebp),%esi
        movl nb200_pos(%ebp),%edi
        movl  nb200_innerjjnr(%esp),%ecx

        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb200_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movl nb200_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb200_ix(%esp),%xmm4
        movapd nb200_iy(%esp),%xmm5
        movapd nb200_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb200_dx(%esp)
        movapd %xmm5,nb200_dy(%esp)
        movapd %xmm6,nb200_dz(%esp)
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

        movapd nb200_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb200_three(%esp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb200_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb200_crf(%esp),%xmm6
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulsd  nb200_two(%esp),%xmm7

        subsd  %xmm7,%xmm0
        mulsd  %xmm0,%xmm3
        mulsd  %xmm3,%xmm4      ## xmm4=total fscal 
        addsd  nb200_vctot(%esp),%xmm6

        movlpd nb200_dx(%esp),%xmm0
        movlpd nb200_dy(%esp),%xmm1
        movlpd nb200_dz(%esp),%xmm2

        movlpd %xmm6,nb200_vctot(%esp)

        movl   nb200_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movlpd nb200_fix(%esp),%xmm3
        movlpd nb200_fiy(%esp),%xmm4
        movlpd nb200_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb200_fix(%esp)
        movlpd %xmm4,nb200_fiy(%esp)
        movlpd %xmm5,nb200_fiz(%esp)
        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel200_ia32_sse2.nb200_updateouterdata: 
        movl  nb200_ii3(%esp),%ecx
        movl  nb200_faction(%ebp),%edi
        movl  nb200_fshift(%ebp),%esi
        movl  nb200_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb200_fix(%esp),%xmm0
        movapd nb200_fiy(%esp),%xmm1
        movapd nb200_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## increment fshift force  
        movsd  (%esi,%edx,8),%xmm3
        movsd  8(%esi,%edx,8),%xmm4
        movsd  16(%esi,%edx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%esi,%edx,8)
        movsd  %xmm4,8(%esi,%edx,8)
        movsd  %xmm5,16(%esi,%edx,8)

        ## get n from stack
        movl nb200_n(%esp),%esi
        ## get group index for i particle 
        movl  nb200_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb200_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb200_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb200_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200_ia32_sse2.nb200_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200_n(%esp)
        jmp _nb_kernel200_ia32_sse2.nb200_outer
_nb_kernel200_ia32_sse2.nb200_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200_ia32_sse2.nb200_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200_ia32_sse2.nb200_threadloop
_nb_kernel200_ia32_sse2.nb200_end: 
        emms

        movl nb200_nouter(%esp),%eax
        movl nb200_ninner(%esp),%ebx
        movl nb200_outeriter(%ebp),%ecx
        movl nb200_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb200_salign(%esp),%eax
        addl %eax,%esp
        addl $304,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret



.globl nb_kernel200nf_ia32_sse2
.globl _nb_kernel200nf_ia32_sse2
nb_kernel200nf_ia32_sse2:       
_nb_kernel200nf_ia32_sse2:      
.set nb200nf_p_nri, 8
.set nb200nf_iinr, 12
.set nb200nf_jindex, 16
.set nb200nf_jjnr, 20
.set nb200nf_shift, 24
.set nb200nf_shiftvec, 28
.set nb200nf_fshift, 32
.set nb200nf_gid, 36
.set nb200nf_pos, 40
.set nb200nf_faction, 44
.set nb200nf_charge, 48
.set nb200nf_p_facel, 52
.set nb200nf_argkrf, 56
.set nb200nf_argcrf, 60
.set nb200nf_Vc, 64
.set nb200nf_type, 68
.set nb200nf_p_ntype, 72
.set nb200nf_vdwparam, 76
.set nb200nf_Vvdw, 80
.set nb200nf_p_tabscale, 84
.set nb200nf_VFtab, 88
.set nb200nf_invsqrta, 92
.set nb200nf_dvda, 96
.set nb200nf_p_gbtabscale, 100
.set nb200nf_GBtab, 104
.set nb200nf_p_nthreads, 108
.set nb200nf_count, 112
.set nb200nf_mtx, 116
.set nb200nf_outeriter, 120
.set nb200nf_inneriter, 124
.set nb200nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb200nf_ix, 0
.set nb200nf_iy, 16
.set nb200nf_iz, 32
.set nb200nf_iq, 48
.set nb200nf_vctot, 64
.set nb200nf_half, 80
.set nb200nf_three, 96
.set nb200nf_krf, 112
.set nb200nf_crf, 128
.set nb200nf_is3, 144
.set nb200nf_ii3, 148
.set nb200nf_innerjjnr, 152
.set nb200nf_innerk, 156
.set nb200nf_n, 160
.set nb200nf_nn1, 164
.set nb200nf_nri, 168
.set nb200nf_facel, 176                       ## uses 8 bytes
.set nb200nf_nouter, 184
.set nb200nf_ninner, 188
.set nb200nf_salign, 192
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $172,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb200nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb200nf_p_nri(%ebp),%ecx
        movl nb200nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb200nf_nri(%esp)
        movsd %xmm7,nb200nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb200nf_nouter(%esp)
        movl %eax,nb200nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb200nf_half(%esp)
        movl %ebx,nb200nf_half+4(%esp)
        movsd nb200nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb200nf_half(%esp)
        movapd %xmm3,nb200nf_three(%esp)

        movl nb200nf_argkrf(%ebp),%esi
        movl nb200nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb200nf_krf(%esp)
        shufpd $0,%xmm6,%xmm6
        movapd %xmm6,nb200nf_crf(%esp)

_nb_kernel200nf_ia32_sse2.nb200nf_threadloop: 
        movl  nb200nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel200nf_ia32_sse2.nb200nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200nf_ia32_sse2.nb200nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200nf_n(%esp)
        movl %ebx,nb200nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200nf_ia32_sse2.nb200nf_outerstart
        jmp _nb_kernel200nf_ia32_sse2.nb200nf_end

_nb_kernel200nf_ia32_sse2.nb200nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200nf_nouter(%esp),%ebx
        movl %ebx,nb200nf_nouter(%esp)

_nb_kernel200nf_ia32_sse2.nb200nf_outer: 
        movl  nb200nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb200nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb200nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb200nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb200nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb200nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb200nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb200nf_ix(%esp)
        movapd %xmm1,nb200nf_iy(%esp)
        movapd %xmm2,nb200nf_iz(%esp)

        movl  %ebx,nb200nf_ii3(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb200nf_vctot(%esp)

        movl  nb200nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb200nf_pos(%ebp),%esi
        movl  nb200nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb200nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb200nf_ninner(%esp),%ecx
        movl  %ecx,nb200nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb200nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel200nf_ia32_sse2.nb200nf_unroll_loop
        jmp   _nb_kernel200nf_ia32_sse2.nb200nf_checksingle
_nb_kernel200nf_ia32_sse2.nb200nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb200nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb200nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb200nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb200nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        movl nb200nf_pos(%ebp),%esi        ## base of pos[] 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb200nf_ix(%esp),%xmm4
        movapd nb200nf_iy(%esp),%xmm5
        movapd nb200nf_iz(%esp),%xmm6

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

        movapd nb200nf_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb200nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb200nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb200nf_crf(%esp),%xmm6
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addpd  nb200nf_vctot(%esp),%xmm6
        movapd %xmm6,nb200nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb200nf_innerk(%esp)
        jl    _nb_kernel200nf_ia32_sse2.nb200nf_checksingle
        jmp   _nb_kernel200nf_ia32_sse2.nb200nf_unroll_loop

_nb_kernel200nf_ia32_sse2.nb200nf_checksingle:  
        movl  nb200nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200nf_ia32_sse2.nb200nf_dosingle
        jmp    _nb_kernel200nf_ia32_sse2.nb200nf_updateouterdata
_nb_kernel200nf_ia32_sse2.nb200nf_dosingle: 
        movl nb200nf_charge(%ebp),%esi
        movl nb200nf_pos(%ebp),%edi
        movl  nb200nf_innerjjnr(%esp),%ecx

        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb200nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movl nb200nf_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb200nf_ix(%esp),%xmm4
        movapd nb200nf_iy(%esp),%xmm5
        movapd nb200nf_iz(%esp),%xmm6

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

        movapd nb200nf_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb200nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb200nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb200nf_crf(%esp),%xmm6
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addsd  nb200nf_vctot(%esp),%xmm6
        movlpd %xmm6,nb200nf_vctot(%esp)

_nb_kernel200nf_ia32_sse2.nb200nf_updateouterdata: 
        ## get n from stack
        movl nb200nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb200nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb200nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb200nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb200nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200nf_ia32_sse2.nb200nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200nf_n(%esp)
        jmp _nb_kernel200nf_ia32_sse2.nb200nf_outer
_nb_kernel200nf_ia32_sse2.nb200nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200nf_ia32_sse2.nb200nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200nf_ia32_sse2.nb200nf_threadloop
_nb_kernel200nf_ia32_sse2.nb200nf_end: 
        emms

        movl nb200nf_nouter(%esp),%eax
        movl nb200nf_ninner(%esp),%ebx
        movl nb200nf_outeriter(%ebp),%ecx
        movl nb200nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb200nf_salign(%esp),%eax
        addl %eax,%esp
        addl $172,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


