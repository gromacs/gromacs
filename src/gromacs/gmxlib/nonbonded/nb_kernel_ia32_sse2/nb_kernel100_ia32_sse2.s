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



.globl nb_kernel100_ia32_sse2
.globl _nb_kernel100_ia32_sse2
nb_kernel100_ia32_sse2: 
_nb_kernel100_ia32_sse2:        
.set nb100_p_nri, 8
.set nb100_iinr, 12
.set nb100_jindex, 16
.set nb100_jjnr, 20
.set nb100_shift, 24
.set nb100_shiftvec, 28
.set nb100_fshift, 32
.set nb100_gid, 36
.set nb100_pos, 40
.set nb100_faction, 44
.set nb100_charge, 48
.set nb100_p_facel, 52
.set nb100_argkrf, 56
.set nb100_argcrf, 60
.set nb100_Vc, 64
.set nb100_type, 68
.set nb100_p_ntype, 72
.set nb100_vdwparam, 76
.set nb100_Vvdw, 80
.set nb100_p_tabscale, 84
.set nb100_VFtab, 88
.set nb100_invsqrta, 92
.set nb100_dvda, 96
.set nb100_p_gbtabscale, 100
.set nb100_GBtab, 104
.set nb100_p_nthreads, 108
.set nb100_count, 112
.set nb100_mtx, 116
.set nb100_outeriter, 120
.set nb100_inneriter, 124
.set nb100_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb100_ix, 0
.set nb100_iy, 16
.set nb100_iz, 32
.set nb100_iq, 48
.set nb100_dx, 64
.set nb100_dy, 80
.set nb100_dz, 96
.set nb100_vctot, 112
.set nb100_fix, 128
.set nb100_fiy, 144
.set nb100_fiz, 160
.set nb100_half, 176
.set nb100_three, 192
.set nb100_is3, 208
.set nb100_ii3, 212
.set nb100_innerjjnr, 216
.set nb100_innerk, 220
.set nb100_n, 224
.set nb100_nn1, 228
.set nb100_nri, 232
.set nb100_facel, 240                         ## uses 8 bytes
.set nb100_nouter, 248
.set nb100_ninner, 252
.set nb100_salign, 256
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $260,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb100_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb100_p_nri(%ebp),%ecx
        movl nb100_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb100_nri(%esp)
        movsd %xmm7,nb100_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb100_nouter(%esp)
        movl %eax,nb100_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb100_half(%esp)
        movl %ebx,nb100_half+4(%esp)
        movsd nb100_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb100_half(%esp)
        movapd %xmm3,nb100_three(%esp)

_nb_kernel100_ia32_sse2.nb100_threadloop: 
        movl  nb100_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel100_ia32_sse2.nb100_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100_ia32_sse2.nb100_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100_n(%esp)
        movl %ebx,nb100_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel100_ia32_sse2.nb100_outerstart
        jmp _nb_kernel100_ia32_sse2.nb100_end

_nb_kernel100_ia32_sse2.nb100_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100_nouter(%esp),%ebx
        movl %ebx,nb100_nouter(%esp)

_nb_kernel100_ia32_sse2.nb100_outer: 
        movl  nb100_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb100_is3(%esp)      ## store is3 

        movl  nb100_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb100_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb100_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb100_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3
        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb100_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb100_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb100_ix(%esp)
        movapd %xmm1,nb100_iy(%esp)
        movapd %xmm2,nb100_iz(%esp)

        movl  %ebx,nb100_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb100_vctot(%esp)
        movapd %xmm4,nb100_fix(%esp)
        movapd %xmm4,nb100_fiy(%esp)
        movapd %xmm4,nb100_fiz(%esp)

        movl  nb100_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb100_pos(%ebp),%esi
        movl  nb100_faction(%ebp),%edi
        movl  nb100_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb100_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb100_ninner(%esp),%ecx
        movl  %ecx,nb100_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb100_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel100_ia32_sse2.nb100_unroll_loop
        jmp   _nb_kernel100_ia32_sse2.nb100_checksingle
_nb_kernel100_ia32_sse2.nb100_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb100_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb100_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb100_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3      ## jq A 
        movhpd (%esi,%ebx,8),%xmm3      ## jq B 

        movapd nb100_iq(%esp),%xmm5

        mulpd %xmm5,%xmm3               ## qq 

        movl nb100_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        movl   nb100_faction(%ebp),%edi

        ## move nb100_ix-iz to xmm4-xmm6 
        movapd nb100_ix(%esp),%xmm4
        movapd nb100_iy(%esp),%xmm5
        movapd nb100_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb100_dx(%esp)
        movapd %xmm5,nb100_dy(%esp)
        movapd %xmm6,nb100_dz(%esp)
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
        movapd nb100_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb100_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 


        movapd nb100_vctot(%esp),%xmm5
        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulpd  %xmm3,%xmm4      ## xmm4=fscal 
        addpd  %xmm3,%xmm5

        movapd nb100_dx(%esp),%xmm0
        movapd nb100_dy(%esp),%xmm1
        movapd nb100_dz(%esp),%xmm2

        movapd %xmm5,nb100_vctot(%esp)

        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb100_fix(%esp),%xmm3
        movapd nb100_fiy(%esp),%xmm4
        movapd nb100_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb100_fix(%esp)
        movapd %xmm4,nb100_fiy(%esp)
        movapd %xmm5,nb100_fiz(%esp)
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
        subl $2,nb100_innerk(%esp)
        jl    _nb_kernel100_ia32_sse2.nb100_checksingle
        jmp   _nb_kernel100_ia32_sse2.nb100_unroll_loop

_nb_kernel100_ia32_sse2.nb100_checksingle:      
        movl  nb100_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel100_ia32_sse2.nb100_dosingle
        jmp    _nb_kernel100_ia32_sse2.nb100_updateouterdata
_nb_kernel100_ia32_sse2.nb100_dosingle: 
        movl nb100_charge(%ebp),%esi
        movl nb100_pos(%ebp),%edi

        movl nb100_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl (%edx),%eax

        xorpd %xmm3,%xmm3
        movsd (%esi,%eax,8),%xmm3       ## jq A 
        movapd nb100_iq(%esp),%xmm5
        unpcklpd %xmm6,%xmm3
        mulpd %xmm5,%xmm3               ## qq 

        movl nb100_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        movl   nb100_faction(%ebp),%edi

        ## move nb100_ix-iz to xmm4-xmm6 
        movapd nb100_ix(%esp),%xmm4
        movapd nb100_iy(%esp),%xmm5
        movapd nb100_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movlpd %xmm4,nb100_dx(%esp)
        movlpd %xmm5,nb100_dy(%esp)
        movlpd %xmm6,nb100_dz(%esp)
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
        movapd nb100_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb100_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 

        movlpd nb100_vctot(%esp),%xmm5
        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulsd  %xmm3,%xmm4      ## xmm4=fscal 
        addsd  %xmm3,%xmm5

        movapd nb100_dx(%esp),%xmm0
        movapd nb100_dy(%esp),%xmm1
        movapd nb100_dz(%esp),%xmm2

        movlpd %xmm5,nb100_vctot(%esp)

        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movlpd nb100_fix(%esp),%xmm3
        movlpd nb100_fiy(%esp),%xmm4
        movlpd nb100_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb100_fix(%esp)
        movlpd %xmm4,nb100_fiy(%esp)
        movlpd %xmm5,nb100_fiz(%esp)
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

_nb_kernel100_ia32_sse2.nb100_updateouterdata: 
        movl  nb100_ii3(%esp),%ecx
        movl  nb100_faction(%ebp),%edi
        movl  nb100_fshift(%ebp),%esi
        movl  nb100_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb100_fix(%esp),%xmm0
        movapd nb100_fiy(%esp),%xmm1
        movapd nb100_fiz(%esp),%xmm2

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
        movl nb100_n(%esp),%esi
        ## get group index for i particle 
        movl  nb100_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb100_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb100_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb100_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100_ia32_sse2.nb100_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100_n(%esp)
        jmp _nb_kernel100_ia32_sse2.nb100_outer
_nb_kernel100_ia32_sse2.nb100_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel100_ia32_sse2.nb100_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100_ia32_sse2.nb100_threadloop
_nb_kernel100_ia32_sse2.nb100_end: 
        emms

        movl nb100_nouter(%esp),%eax
        movl nb100_ninner(%esp),%ebx
        movl nb100_outeriter(%ebp),%ecx
        movl nb100_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb100_salign(%esp),%eax
        addl %eax,%esp
        addl $260,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel100nf_ia32_sse2
.globl _nb_kernel100nf_ia32_sse2
nb_kernel100nf_ia32_sse2:       
_nb_kernel100nf_ia32_sse2:      
.set nb100nf_p_nri, 8
.set nb100nf_iinr, 12
.set nb100nf_jindex, 16
.set nb100nf_jjnr, 20
.set nb100nf_shift, 24
.set nb100nf_shiftvec, 28
.set nb100nf_fshift, 32
.set nb100nf_gid, 36
.set nb100nf_pos, 40
.set nb100nf_faction, 44
.set nb100nf_charge, 48
.set nb100nf_p_facel, 52
.set nb100nf_argkrf, 56
.set nb100nf_argcrf, 60
.set nb100nf_Vc, 64
.set nb100nf_type, 68
.set nb100nf_p_ntype, 72
.set nb100nf_vdwparam, 76
.set nb100nf_Vvdw, 80
.set nb100nf_p_tabscale, 84
.set nb100nf_VFtab, 88
.set nb100nf_invsqrta, 92
.set nb100nf_dvda, 96
.set nb100nf_p_gbtabscale, 100
.set nb100nf_GBtab, 104
.set nb100nf_p_nthreads, 108
.set nb100nf_count, 112
.set nb100nf_mtx, 116
.set nb100nf_outeriter, 120
.set nb100nf_inneriter, 124
.set nb100nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb100nf_ix, 0
.set nb100nf_iy, 16
.set nb100nf_iz, 32
.set nb100nf_iq, 48
.set nb100nf_vctot, 64
.set nb100nf_half, 80
.set nb100nf_three, 96
.set nb100nf_is3, 112
.set nb100nf_ii3, 116
.set nb100nf_innerjjnr, 120
.set nb100nf_innerk, 124
.set nb100nf_n, 128
.set nb100nf_nn1, 132
.set nb100nf_nri, 136
.set nb100nf_facel, 144                       ## uses 8 bytes
.set nb100nf_nouter, 152
.set nb100nf_ninner, 156
.set nb100nf_salign, 160
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $164,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb100nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb100nf_p_nri(%ebp),%ecx
        movl nb100nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb100nf_nri(%esp)
        movsd %xmm7,nb100nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb100nf_nouter(%esp)
        movl %eax,nb100nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb100nf_half(%esp)
        movl %ebx,nb100nf_half+4(%esp)
        movsd nb100nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb100nf_half(%esp)
        movapd %xmm3,nb100nf_three(%esp)

_nb_kernel100nf_ia32_sse2.nb100nf_threadloop: 
        movl  nb100nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel100nf_ia32_sse2.nb100nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100nf_ia32_sse2.nb100nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100nf_n(%esp)
        movl %ebx,nb100nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel100nf_ia32_sse2.nb100nf_outerstart
        jmp _nb_kernel100nf_ia32_sse2.nb100nf_end

_nb_kernel100nf_ia32_sse2.nb100nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100nf_nouter(%esp),%ebx
        movl %ebx,nb100nf_nouter(%esp)

_nb_kernel100nf_ia32_sse2.nb100nf_outer: 
        movl  nb100nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb100nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb100nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb100nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb100nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3
        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb100nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb100nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb100nf_ix(%esp)
        movapd %xmm1,nb100nf_iy(%esp)
        movapd %xmm2,nb100nf_iz(%esp)

        movl  %ebx,nb100nf_ii3(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb100nf_vctot(%esp)

        movl  nb100nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb100nf_pos(%ebp),%esi
        movl  nb100nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb100nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb100nf_ninner(%esp),%ecx
        movl  %ecx,nb100nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb100nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel100nf_ia32_sse2.nb100nf_unroll_loop
        jmp   _nb_kernel100nf_ia32_sse2.nb100nf_checksingle
_nb_kernel100nf_ia32_sse2.nb100nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb100nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb100nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb100nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3      ## jq A 
        movhpd (%esi,%ebx,8),%xmm3      ## jq B 

        movapd nb100nf_iq(%esp),%xmm5

        mulpd %xmm5,%xmm3               ## qq 

        movl nb100nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move nb100nf_ix-iz to xmm4-xmm6 
        movapd nb100nf_ix(%esp),%xmm4
        movapd nb100nf_iy(%esp),%xmm5
        movapd nb100nf_iz(%esp),%xmm6

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
        movapd nb100nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb100nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        movapd nb100nf_vctot(%esp),%xmm5
        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        addpd  %xmm3,%xmm5
        movapd %xmm5,nb100nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb100nf_innerk(%esp)
        jl    _nb_kernel100nf_ia32_sse2.nb100nf_checksingle
        jmp   _nb_kernel100nf_ia32_sse2.nb100nf_unroll_loop

_nb_kernel100nf_ia32_sse2.nb100nf_checksingle:  
        movl  nb100nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel100nf_ia32_sse2.nb100nf_dosingle
        jmp    _nb_kernel100nf_ia32_sse2.nb100nf_updateouterdata
_nb_kernel100nf_ia32_sse2.nb100nf_dosingle: 
        movl nb100nf_charge(%ebp),%esi
        movl nb100nf_pos(%ebp),%edi

        movl nb100nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl (%edx),%eax

        xorpd %xmm3,%xmm3
        movsd (%esi,%eax,8),%xmm3       ## jq A 
        movapd nb100nf_iq(%esp),%xmm5
        unpcklpd %xmm6,%xmm3
        mulpd %xmm5,%xmm3               ## qq 

        movl nb100nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move nb100nf_ix-iz to xmm4-xmm6 
        movapd nb100nf_ix(%esp),%xmm4
        movapd nb100nf_iy(%esp),%xmm5
        movapd nb100nf_iz(%esp),%xmm6

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
        movapd nb100nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb100nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        movlpd nb100nf_vctot(%esp),%xmm5
        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        addsd  %xmm3,%xmm5
        movlpd %xmm5,nb100nf_vctot(%esp)

_nb_kernel100nf_ia32_sse2.nb100nf_updateouterdata: 
        ## get n from stack
        movl nb100nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb100nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb100nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb100nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb100nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100nf_ia32_sse2.nb100nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100nf_n(%esp)
        jmp _nb_kernel100nf_ia32_sse2.nb100nf_outer
_nb_kernel100nf_ia32_sse2.nb100nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel100nf_ia32_sse2.nb100nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100nf_ia32_sse2.nb100nf_threadloop
_nb_kernel100nf_ia32_sse2.nb100nf_end: 
        emms

        movl nb100nf_nouter(%esp),%eax
        movl nb100nf_ninner(%esp),%ebx
        movl nb100nf_outeriter(%ebp),%ecx
        movl nb100nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb100nf_salign(%esp),%eax
        addl %eax,%esp
        addl $164,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





