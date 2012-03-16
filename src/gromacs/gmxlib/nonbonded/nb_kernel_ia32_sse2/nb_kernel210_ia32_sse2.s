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




.globl nb_kernel210_ia32_sse2
.globl _nb_kernel210_ia32_sse2
nb_kernel210_ia32_sse2: 
_nb_kernel210_ia32_sse2:        
.set nb210_p_nri, 8
.set nb210_iinr, 12
.set nb210_jindex, 16
.set nb210_jjnr, 20
.set nb210_shift, 24
.set nb210_shiftvec, 28
.set nb210_fshift, 32
.set nb210_gid, 36
.set nb210_pos, 40
.set nb210_faction, 44
.set nb210_charge, 48
.set nb210_p_facel, 52
.set nb210_argkrf, 56
.set nb210_argcrf, 60
.set nb210_Vc, 64
.set nb210_type, 68
.set nb210_p_ntype, 72
.set nb210_vdwparam, 76
.set nb210_Vvdw, 80
.set nb210_p_tabscale, 84
.set nb210_VFtab, 88
.set nb210_invsqrta, 92
.set nb210_dvda, 96
.set nb210_p_gbtabscale, 100
.set nb210_GBtab, 104
.set nb210_p_nthreads, 108
.set nb210_count, 112
.set nb210_mtx, 116
.set nb210_outeriter, 120
.set nb210_inneriter, 124
.set nb210_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb210_ix, 0
.set nb210_iy, 16
.set nb210_iz, 32
.set nb210_iq, 48
.set nb210_dx, 64
.set nb210_dy, 80
.set nb210_dz, 96
.set nb210_c6, 112
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
.set nb210_facel, 336                       ## uses 8 bytes
.set nb210_is3, 344
.set nb210_ii3, 348
.set nb210_ntia, 352
.set nb210_innerjjnr, 356
.set nb210_innerk, 360
.set nb210_n, 364
.set nb210_nn1, 368
.set nb210_nri, 372
.set nb210_ntype, 376
.set nb210_nouter, 380
.set nb210_ninner, 384
.set nb210_salign, 388
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $392,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb210_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb210_p_nri(%ebp),%ecx
        movl nb210_p_facel(%ebp),%esi
        movl nb210_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb210_nri(%esp)
        movsd %xmm7,nb210_facel(%esp)
        movl %edi,nb210_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb210_nouter(%esp)
        movl %eax,nb210_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb210_half(%esp)
        movl %ebx,nb210_half+4(%esp)
        movsd nb210_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## 6.0
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## 12.0
        movapd %xmm1,nb210_half(%esp)
        movapd %xmm2,nb210_two(%esp)
        movapd %xmm3,nb210_three(%esp)
        movapd %xmm4,nb210_six(%esp)
        movapd %xmm5,nb210_twelve(%esp)

        movl nb210_argkrf(%ebp),%esi
        movl nb210_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb210_krf(%esp)
        movapd %xmm6,nb210_crf(%esp)

_nb_kernel210_ia32_sse2.nb210_threadloop: 
        movl  nb210_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel210_ia32_sse2.nb210_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210_ia32_sse2.nb210_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210_n(%esp)
        movl %ebx,nb210_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210_ia32_sse2.nb210_outerstart
        jmp _nb_kernel210_ia32_sse2.nb210_end

_nb_kernel210_ia32_sse2.nb210_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210_nouter(%esp),%ebx
        movl %ebx,nb210_nouter(%esp)

_nb_kernel210_ia32_sse2.nb210_outer: 
        movl  nb210_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb210_is3(%esp)      ## store is3 

        movl  nb210_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb210_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb210_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb210_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb210_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb210_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb210_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb210_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb210_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb210_ix(%esp)
        movapd %xmm1,nb210_iy(%esp)
        movapd %xmm2,nb210_iz(%esp)

        movl  %ebx,nb210_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb210_vctot(%esp)
        movapd %xmm4,nb210_Vvdwtot(%esp)
        movapd %xmm4,nb210_fix(%esp)
        movapd %xmm4,nb210_fiy(%esp)
        movapd %xmm4,nb210_fiz(%esp)

        movl  nb210_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb210_pos(%ebp),%esi
        movl  nb210_faction(%ebp),%edi
        movl  nb210_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb210_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb210_ninner(%esp),%ecx
        movl  %ecx,nb210_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb210_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel210_ia32_sse2.nb210_unroll_loop
        jmp   _nb_kernel210_ia32_sse2.nb210_checksingle
_nb_kernel210_ia32_sse2.nb210_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb210_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb210_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb210_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb210_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb210_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb210_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb210_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb210_c6(%esp)
        movapd %xmm6,nb210_c12(%esp)

        movl nb210_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb210_ix(%esp),%xmm4
        movapd nb210_iy(%esp),%xmm5
        movapd nb210_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb210_dx(%esp)
        movapd %xmm5,nb210_dy(%esp)
        movapd %xmm6,nb210_dz(%esp)
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

        movapd nb210_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb210_three(%esp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb210_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb210_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb210_crf(%esp),%xmm6
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulpd  nb210_two(%esp),%xmm7
        mulpd  nb210_c6(%esp),%xmm1
        mulpd  nb210_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb210_Vvdwtot(%esp),%xmm5
        mulpd  nb210_six(%esp),%xmm1
        mulpd  nb210_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        subpd  %xmm7,%xmm0
        mulpd  %xmm0,%xmm3
        addpd  %xmm3,%xmm2
        mulpd  %xmm2,%xmm4      ## xmm4=total fscal 
        addpd  nb210_vctot(%esp),%xmm6

        movapd nb210_dx(%esp),%xmm0
        movapd nb210_dy(%esp),%xmm1
        movapd nb210_dz(%esp),%xmm2

        movapd %xmm6,nb210_vctot(%esp)
        movapd %xmm5,nb210_Vvdwtot(%esp)

        movl   nb210_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb210_fix(%esp),%xmm3
        movapd nb210_fiy(%esp),%xmm4
        movapd nb210_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb210_fix(%esp)
        movapd %xmm4,nb210_fiy(%esp)
        movapd %xmm5,nb210_fiz(%esp)
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
        subl $2,nb210_innerk(%esp)
        jl    _nb_kernel210_ia32_sse2.nb210_checksingle
        jmp   _nb_kernel210_ia32_sse2.nb210_unroll_loop

_nb_kernel210_ia32_sse2.nb210_checksingle:      
        movl  nb210_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel210_ia32_sse2.nb210_dosingle
        jmp    _nb_kernel210_ia32_sse2.nb210_updateouterdata
_nb_kernel210_ia32_sse2.nb210_dosingle: 
        movl nb210_charge(%ebp),%esi
        movl nb210_pos(%ebp),%edi
        movl  nb210_innerjjnr(%esp),%ecx
        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb210_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb210_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb210_vdwparam(%ebp),%esi
        shll %eax
        movl nb210_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb210_c6(%esp)
        movapd %xmm6,nb210_c12(%esp)

        movl nb210_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb210_ix(%esp),%xmm4
        movapd nb210_iy(%esp),%xmm5
        movapd nb210_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb210_dx(%esp)
        movapd %xmm5,nb210_dy(%esp)
        movapd %xmm6,nb210_dz(%esp)
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

        movapd nb210_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb210_three(%esp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb210_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb210_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb210_crf(%esp),%xmm6
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulsd  nb210_two(%esp),%xmm7
        mulsd  nb210_c6(%esp),%xmm1
        mulsd  nb210_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb210_Vvdwtot(%esp),%xmm5
        mulsd  nb210_six(%esp),%xmm1
        mulsd  nb210_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        subsd  %xmm7,%xmm0
        mulsd  %xmm0,%xmm3
        addsd  %xmm3,%xmm2
        mulsd  %xmm2,%xmm4      ## xmm4=total fscal 
        addsd  nb210_vctot(%esp),%xmm6

        movlpd nb210_dx(%esp),%xmm0
        movlpd nb210_dy(%esp),%xmm1
        movlpd nb210_dz(%esp),%xmm2

        movlpd %xmm6,nb210_vctot(%esp)
        movlpd %xmm5,nb210_Vvdwtot(%esp)

        movl   nb210_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movlpd nb210_fix(%esp),%xmm3
        movlpd nb210_fiy(%esp),%xmm4
        movlpd nb210_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb210_fix(%esp)
        movlpd %xmm4,nb210_fiy(%esp)
        movlpd %xmm5,nb210_fiz(%esp)
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

_nb_kernel210_ia32_sse2.nb210_updateouterdata: 
        movl  nb210_ii3(%esp),%ecx
        movl  nb210_faction(%ebp),%edi
        movl  nb210_fshift(%ebp),%esi
        movl  nb210_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb210_fix(%esp),%xmm0
        movapd nb210_fiy(%esp),%xmm1
        movapd nb210_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addpd  %xmm3,%xmm0
        addpd  %xmm4,%xmm1
        addpd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movl nb210_n(%esp),%esi
        ## get group index for i particle 
        movl  nb210_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb210_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb210_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb210_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb210_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb210_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210_ia32_sse2.nb210_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210_n(%esp)
        jmp _nb_kernel210_ia32_sse2.nb210_outer
_nb_kernel210_ia32_sse2.nb210_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210_ia32_sse2.nb210_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210_ia32_sse2.nb210_threadloop
_nb_kernel210_ia32_sse2.nb210_end: 
        emms

        movl nb210_nouter(%esp),%eax
        movl nb210_ninner(%esp),%ebx
        movl nb210_outeriter(%ebp),%ecx
        movl nb210_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb210_salign(%esp),%eax
        addl %eax,%esp
        addl $392,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel210nf_ia32_sse2
.globl _nb_kernel210nf_ia32_sse2
nb_kernel210nf_ia32_sse2:       
_nb_kernel210nf_ia32_sse2:      
.set nb210nf_p_nri, 8
.set nb210nf_iinr, 12
.set nb210nf_jindex, 16
.set nb210nf_jjnr, 20
.set nb210nf_shift, 24
.set nb210nf_shiftvec, 28
.set nb210nf_fshift, 32
.set nb210nf_gid, 36
.set nb210nf_pos, 40
.set nb210nf_faction, 44
.set nb210nf_charge, 48
.set nb210nf_p_facel, 52
.set nb210nf_argkrf, 56
.set nb210nf_argcrf, 60
.set nb210nf_Vc, 64
.set nb210nf_type, 68
.set nb210nf_p_ntype, 72
.set nb210nf_vdwparam, 76
.set nb210nf_Vvdw, 80
.set nb210nf_p_tabscale, 84
.set nb210nf_VFtab, 88
.set nb210nf_invsqrta, 92
.set nb210nf_dvda, 96
.set nb210nf_p_gbtabscale, 100
.set nb210nf_GBtab, 104
.set nb210nf_p_nthreads, 108
.set nb210nf_count, 112
.set nb210nf_mtx, 116
.set nb210nf_outeriter, 120
.set nb210nf_inneriter, 124
.set nb210nf_work, 128
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
.set nb210nf_is3, 192
.set nb210nf_ii3, 196
.set nb210nf_ntia, 200
.set nb210nf_innerjjnr, 204
.set nb210nf_innerk, 208
.set nb210nf_n, 212
.set nb210nf_nn1, 216
.set nb210nf_nri, 220
.set nb210nf_facel, 224                       ## uses 8 bytes
.set nb210nf_ntype, 232
.set nb210nf_nouter, 236
.set nb210nf_ninner, 240
.set nb210nf_salign, 244
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $224,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb210nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb210nf_p_nri(%ebp),%ecx
        movl nb210nf_p_facel(%ebp),%esi
        movl nb210nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb210nf_nri(%esp)
        movsd %xmm7,nb210nf_facel(%esp)
        movl %edi,nb210nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb210nf_nouter(%esp)
        movl %eax,nb210nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb210nf_half(%esp)
        movl %ebx,nb210nf_half+4(%esp)
        movsd nb210nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb210nf_half(%esp)
        movapd %xmm3,nb210nf_three(%esp)

        movl nb210nf_argkrf(%ebp),%esi
        movl nb210nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb210nf_krf(%esp)
        movapd %xmm6,nb210nf_crf(%esp)

_nb_kernel210nf_ia32_sse2.nb210nf_threadloop: 
        movl  nb210nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel210nf_ia32_sse2.nb210nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210nf_ia32_sse2.nb210nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210nf_n(%esp)
        movl %ebx,nb210nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210nf_ia32_sse2.nb210nf_outerstart
        jmp _nb_kernel210nf_ia32_sse2.nb210nf_end

_nb_kernel210nf_ia32_sse2.nb210nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210nf_nouter(%esp),%ebx
        movl %ebx,nb210nf_nouter(%esp)

_nb_kernel210nf_ia32_sse2.nb210nf_outer: 
        movl  nb210nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb210nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb210nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb210nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb210nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb210nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb210nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb210nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb210nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb210nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb210nf_ix(%esp)
        movapd %xmm1,nb210nf_iy(%esp)
        movapd %xmm2,nb210nf_iz(%esp)

        movl  %ebx,nb210nf_ii3(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb210nf_vctot(%esp)
        movapd %xmm4,nb210nf_Vvdwtot(%esp)

        movl  nb210nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb210nf_pos(%ebp),%esi
        movl  nb210nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb210nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb210nf_ninner(%esp),%ecx
        movl  %ecx,nb210nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb210nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel210nf_ia32_sse2.nb210nf_unroll_loop
        jmp   _nb_kernel210nf_ia32_sse2.nb210nf_checksingle
_nb_kernel210nf_ia32_sse2.nb210nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb210nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb210nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb210nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb210nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb210nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb210nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb210nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb210nf_c6(%esp)
        movapd %xmm6,nb210nf_c12(%esp)

        movl nb210nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb210nf_ix(%esp),%xmm4
        movapd nb210nf_iy(%esp),%xmm5
        movapd nb210nf_iz(%esp),%xmm6

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

        movapd nb210nf_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb210nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb210nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb210nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb210nf_crf(%esp),%xmm6
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulpd  nb210nf_c6(%esp),%xmm1
        mulpd  nb210nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb210nf_Vvdwtot(%esp),%xmm5
        addpd  nb210nf_vctot(%esp),%xmm6
        movapd %xmm6,nb210nf_vctot(%esp)
        movapd %xmm5,nb210nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $2,nb210nf_innerk(%esp)
        jl    _nb_kernel210nf_ia32_sse2.nb210nf_checksingle
        jmp   _nb_kernel210nf_ia32_sse2.nb210nf_unroll_loop

_nb_kernel210nf_ia32_sse2.nb210nf_checksingle:  
        movl  nb210nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel210nf_ia32_sse2.nb210nf_dosingle
        jmp    _nb_kernel210nf_ia32_sse2.nb210nf_updateouterdata
_nb_kernel210nf_ia32_sse2.nb210nf_dosingle: 
        movl nb210nf_charge(%ebp),%esi
        movl nb210nf_pos(%ebp),%edi
        movl  nb210nf_innerjjnr(%esp),%ecx
        xorpd %xmm3,%xmm3
        movl  (%ecx),%eax

        movlpd (%esi,%eax,8),%xmm3
        movapd nb210nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb210nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb210nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb210nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb210nf_c6(%esp)
        movapd %xmm6,nb210nf_c12(%esp)

        movl nb210nf_pos(%ebp),%esi        ## base of pos[] 

        leal (%eax,%eax,2),%eax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb210nf_ix(%esp),%xmm4
        movapd nb210nf_iy(%esp),%xmm5
        movapd nb210nf_iz(%esp),%xmm6

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

        movapd nb210nf_krf(%esp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb210nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb210nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb210nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb210nf_crf(%esp),%xmm6
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulsd  nb210nf_c6(%esp),%xmm1
        mulsd  nb210nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb210nf_Vvdwtot(%esp),%xmm5
        addsd  nb210nf_vctot(%esp),%xmm6
        movlpd %xmm6,nb210nf_vctot(%esp)
        movlpd %xmm5,nb210nf_Vvdwtot(%esp)

_nb_kernel210nf_ia32_sse2.nb210nf_updateouterdata: 
        ## get n from stack
        movl nb210nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb210nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb210nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb210nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb210nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb210nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb210nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210nf_ia32_sse2.nb210nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210nf_n(%esp)
        jmp _nb_kernel210nf_ia32_sse2.nb210nf_outer
_nb_kernel210nf_ia32_sse2.nb210nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210nf_ia32_sse2.nb210nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210nf_ia32_sse2.nb210nf_threadloop
_nb_kernel210nf_ia32_sse2.nb210nf_end: 
        emms

        movl nb210nf_nouter(%esp),%eax
        movl nb210nf_ninner(%esp),%ebx
        movl nb210nf_outeriter(%ebp),%ecx
        movl nb210nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb210nf_salign(%esp),%eax
        addl %eax,%esp
        addl $224,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



